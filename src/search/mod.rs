// SPDX-License-Identifier: Apache-2.0
//! Monte-Carlo global search (docking).
//!
//! Each task runs an independent Metropolis Monte-Carlo over conformation space:
//! randomize a start pose, then repeatedly mutate one degree of freedom, locally
//! optimize on the grid cache (BFGS), and accept by the Metropolis criterion,
//! saving improving local minima. Tasks' minima are merged, re-optimized off the
//! grid (`non_cache`), rescored, and ranked.
//!
//! The RNG is not Boost's, so the search trajectory is not bit-identical to
//! Vina's; correctness is validated by convergence to the same minima.

pub mod mutate;
pub mod output;

use crate::model::conf::Conf;
use crate::model::Model;
use crate::random::Rng;
use crate::scoring::eval::{
    assemble_breakdown, default_max_steps, local_optimize, local_optimize_cache, pair_deriv_interp,
    score_components, ScoreBreakdown, SearchBox,
};
use crate::scoring::grid::Cache;
use crate::scoring::ScoringFunction;

use mutate::mutate_conf;
use output::{OutputContainer, OutputType};

/// The `hunt_cap` energy bound used during search BFGS.
const HUNT_CAP_V: f64 = 10.0;
/// `authentic_v` used for the final re-optimization of accepted minima.
const AUTHENTIC_V: f64 = 1000.0;

/// Monte-Carlo parameters.
#[derive(Debug, Clone, Copy)]
pub struct MonteCarloParams {
    pub global_steps: usize,
    pub local_steps: usize,
    pub temperature: f64,
    pub mutation_amplitude: f64,
    pub num_saved_mins: usize,
    pub min_rmsd: f64,
}

impl MonteCarloParams {
    /// The default search parameters computed for a model.
    pub fn for_model(m: &Model, num_saved_mins: usize, min_rmsd: f64) -> Self {
        let heuristic = m.num_movable_atoms() + 10 * m.get_size().num_degrees_of_freedom();
        let global_steps = 70 * 3 * (50 + heuristic) / 2;
        MonteCarloParams {
            global_steps,
            local_steps: default_max_steps(m),
            temperature: 1.2,
            mutation_amplitude: 2.0,
            num_saved_mins,
            min_rmsd,
        }
    }
}

/// The Metropolis acceptance test.
fn metropolis_accept(old_f: f64, new_f: f64, temperature: f64, rng: &mut Rng) -> bool {
    if new_f < old_f {
        return true;
    }
    let p = ((old_f - new_f) / temperature).exp();
    rng.random_fl(0.0, 1.0) < p
}

/// One Monte-Carlo run.
fn monte_carlo_run(
    m: &mut Model,
    cache: &Cache,
    sf: &ScoringFunction,
    params: &MonteCarloParams,
    corner1: &crate::math::Vec3,
    corner2: &crate::math::Vec3,
    rng: &mut Rng,
) -> OutputContainer {
    let s = m.get_size();
    let mut out = OutputContainer::new();
    let mut tmp = OutputType::new(Conf::new(&s), 0.0);
    tmp.conf.randomize(corner1, corner2, rng);
    let mut best_e = f64::MAX;

    for step in 0..params.global_steps {
        let mut candidate = tmp.clone();
        mutate_conf(&mut candidate.conf, m, params.mutation_amplitude, rng);
        let (e, _) = local_optimize_cache(
            m,
            sf,
            cache,
            HUNT_CAP_V,
            &mut candidate.conf,
            params.local_steps,
        );
        candidate.e = e;

        if step == 0 || metropolis_accept(tmp.e, candidate.e, params.temperature, rng) {
            tmp = candidate;
            if tmp.e < best_e || out.len() < params.num_saved_mins {
                // Re-optimize accepted minima with the authentic cap.
                let (e, _) = local_optimize_cache(
                    m,
                    sf,
                    cache,
                    AUTHENTIC_V,
                    &mut tmp.conf,
                    params.local_steps,
                );
                tmp.e = e;
                tmp.coords = m.heavy_atom_movable_coords();
                out.add(tmp.clone(), params.min_rmsd, params.num_saved_mins);
                if tmp.e < best_e {
                    best_e = tmp.e;
                }
            }
        }
    }
    out
}

/// A finalized docked pose: the full score breakdown plus its conformation.
#[derive(Debug, Clone)]
pub struct DockedPose {
    pub breakdown: ScoreBreakdown,
    pub conf: Conf,
}

/// Run the full docking search: independent MC tasks, merged minima, non_cache
/// refinement, rescoring, and ranking. Returns poses sorted best-first.
#[allow(clippy::too_many_arguments)]
pub fn dock(
    base: &Model,
    cache: &Cache,
    sf: &ScoringFunction,
    gbox: &SearchBox,
    exhaustiveness: usize,
    num_modes: usize,
    min_rmsd: f64,
    seed: u32,
    global_steps_override: Option<usize>,
) -> Vec<DockedPose> {
    let corner1 = crate::math::Vec3::new(gbox.begin[0], gbox.begin[1], gbox.begin[2]);
    let corner2 = crate::math::Vec3::new(gbox.end[0], gbox.end[1], gbox.end[2]);
    let mut params = MonteCarloParams::for_model(base, num_modes, min_rmsd);
    if let Some(gs) = global_steps_override {
        params.global_steps = gs;
    }

    // Per-task seeds drawn from the main generator.
    let mut main_rng = Rng::seed(seed);
    let mut merged = OutputContainer::new();
    const MERGE_MIN_RMSD: f64 = 2.0; // merge cutoff for combining task minima
    for _ in 0..exhaustiveness {
        let task_seed = main_rng.random_int(0, 1_000_000) as u32;
        let mut task_rng = Rng::seed(task_seed);
        let mut task_model = base.clone();
        let task_out = monte_carlo_run(
            &mut task_model,
            cache,
            sf,
            &params,
            &corner1,
            &corner2,
            &mut task_rng,
        );
        for pose in task_out.poses {
            merged.add(pose, MERGE_MIN_RMSD, params.num_saved_mins);
        }
    }
    merged.sort();

    // Remove redundant poses by the output min_rmsd.
    let mut deduped = OutputContainer::new();
    for pose in merged.poses {
        deduped.add(pose, min_rmsd, merged_capacity(exhaustiveness, num_modes));
    }

    // Refine each pose off-grid (non_cache) and record its raw energy components.
    // The canonical refinement uses a 5-step out-of-box penalty ramp (slope
    // 1e2..1e10, breaking once the pose is back in the box); for a pose that stays
    // inside the box the penalty and its gradient are identically zero at every
    // slope, so a single pass at the fixed slope is equivalent. Out-of-box
    // stragglers (atypical after the grid-clamped search) would differ — not
    // handled here.
    let mut model = base.clone();
    let refine_steps = default_max_steps(&model);
    let pair = pair_deriv_interp(sf);
    // (inter, intra, conf) per pose.
    let mut refined: Vec<(f64, f64, Conf)> = Vec::new();
    for pose in &deduped.poses {
        let mut conf = pose.conf.clone();
        let _ = local_optimize(
            &mut model,
            sf,
            gbox,
            AUTHENTIC_V,
            &pair,
            &mut conf,
            refine_steps,
        );
        // local_optimize leaves the model `set` to `conf`.
        let (inter, intra) = score_components(&model, sf, gbox);
        refined.push((inter, intra, conf));
    }

    // Rank by the raw search energy `inter + intra`, then take the unbound
    // reference from the best pose's intra and report every pose against that
    // shared reference.
    refined.sort_by(|a, b| {
        (a.0 + a.1)
            .partial_cmp(&(b.0 + b.1))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let num_tors = model.num_tors();
    let unbound = refined.first().map(|p| p.1).unwrap_or(0.0);
    let mut results: Vec<DockedPose> = refined
        .into_iter()
        .map(|(inter, intra, conf)| DockedPose {
            breakdown: assemble_breakdown(inter, intra, unbound, num_tors, sf),
            conf,
        })
        .collect();
    results.truncate(num_modes);
    results
}

fn merged_capacity(exhaustiveness: usize, num_modes: usize) -> usize {
    (exhaustiveness * num_modes).max(num_modes)
}
