// SPDX-License-Identifier: Apache-2.0
//! Single-point scoring: the `--score_only` energy breakdown.
//!
//! Scores the Vina scoring function with a rigid receptor: the intermolecular
//! term via direct pairwise summation, the intramolecular term over the
//! ligand's interacting pairs, and the conformation-independent torsional term.
//!
//! The 3-element coordinate loops index several parallel arrays, so the
//! `needless_range_loop` lint (which assumes a single iterable) is allowed here.
#![allow(clippy::needless_range_loop)]

use crate::atom::constants::{is_glue_type, XS_TYPE_SIZE};
use crate::math::{distance_sqr, Vec3, ZERO};
use crate::model::conf::{Change, Conf};
use crate::model::Model;

use super::{curl, curl_deriv, ScoringFunction};

/// The energy cap used for scoring.
const AUTHENTIC_V: f64 = 1000.0;

/// The axis-aligned search box. For an in-box pose the box only gates the
/// out-of-bounds penalty (which is then zero).
#[derive(Debug, Clone, Copy)]
pub struct SearchBox {
    pub begin: [f64; 3],
    pub end: [f64; 3],
    /// Number of voxels per axis (sample points = `n_voxels + 1`).
    pub n_voxels: [usize; 3],
    /// Out-of-bounds penalty slope (irrelevant for in-box poses).
    pub slope: f64,
}

/// The default grid spacing for the Vina CLI (`--spacing 0.375`).
pub const DEFAULT_GRANULARITY: f64 = 0.375;

impl SearchBox {
    /// Construct the box: `n_voxels = ceil(size/granularity)` (rounded up to
    /// even if `force_even_voxels`), `real_span = granularity * n_voxels`,
    /// `begin = center - real_span/2`, `end = begin + real_span`.
    pub fn new(
        center: [f64; 3],
        size: [f64; 3],
        granularity: f64,
        force_even_voxels: bool,
    ) -> Self {
        let mut begin = [0.0; 3];
        let mut end = [0.0; 3];
        let mut nv = [0usize; 3];
        for j in 0..3 {
            let mut n_voxels = (size[j] / granularity).ceil() as usize;
            if force_even_voxels && n_voxels % 2 == 1 {
                n_voxels += 1;
            }
            let real_span = granularity * n_voxels as f64;
            begin[j] = center[j] - real_span / 2.0;
            end[j] = begin[j] + real_span;
            nv[j] = n_voxels;
        }
        SearchBox {
            begin,
            end,
            n_voxels: nv,
            slope: 1.0e6,
        }
    }

    /// Span (`end - begin`) per axis.
    pub fn span(&self, j: usize) -> f64 {
        self.end[j] - self.begin[j]
    }

    /// Convenience: the box with the Vina CLI defaults (granularity 0.375, no
    /// forced-even voxels).
    pub fn from_center_size(center: [f64; 3], size: [f64; 3]) -> Self {
        Self::new(center, size, DEFAULT_GRANULARITY, false)
    }
}

/// Vina's reported single-pose energy breakdown.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScoreBreakdown {
    /// Estimated free energy of binding (the headline affinity).
    pub estimated_free_energy: f64,
    /// (1) Final intermolecular energy.
    pub intermolecular: f64,
    /// (2) Final total internal energy.
    pub internal: f64,
    /// (3) Torsional free energy.
    pub torsional: f64,
    /// (4) Unbound system's energy.
    pub unbound: f64,
}

/// Intermolecular ligand-receptor energy by direct pairwise summation over
/// receptor atoms. A neighbour list would be a pure speed optimization;
/// iterating all receptor atoms with the cutoff filter yields the same pairs.
fn eval_inter(model: &Model, sf: &ScoringFunction, gbox: &SearchBox, v: f64) -> f64 {
    let cutoff_sqr = sf.cutoff_sqr();
    let mut e = 0.0;
    for i in 0..model.num_movable_atoms() {
        let a = &model.atoms[i];
        let t1 = a.xs();
        if t1 >= XS_TYPE_SIZE || is_glue_type(t1) {
            continue;
        }
        let a_coords = model.coords[i];

        // Clamp to the box, accumulating the out-of-bounds penalty.
        let mut adjusted = a_coords;
        let mut penalty = 0.0;
        for j in 0..3 {
            if a_coords[j] < gbox.begin[j] {
                adjusted[j] = gbox.begin[j];
                penalty += (a_coords[j] - gbox.begin[j]).abs();
            } else if a_coords[j] > gbox.end[j] {
                adjusted[j] = gbox.end[j];
                penalty += (a_coords[j] - gbox.end[j]).abs();
            }
        }
        penalty *= gbox.slope;

        let mut this_e = 0.0;
        for b in &model.grid_atoms {
            let t2 = b.xs();
            if t2 >= XS_TYPE_SIZE {
                continue;
            }
            let r_ba = adjusted - b.coords;
            let r2 = r_ba.norm_sqr();
            if r2 < cutoff_sqr {
                this_e += sf.eval_fast(t1, t2, r2);
            }
        }
        curl(&mut this_e, v);
        e += this_e + penalty;
    }
    e
}

/// Sum over each ligand's interacting pairs (the intramolecular energy, which
/// is also Vina's unbound reference).
fn eval_intra(model: &Model, sf: &ScoringFunction, v: f64) -> f64 {
    let cutoff_sqr = sf.cutoff_sqr();
    let mut e = 0.0;
    for lig in &model.ligands {
        for &(a, b) in &lig.pairs {
            let r2 = distance_sqr(&model.coords[a], &model.coords[b]);
            if r2 < cutoff_sqr {
                let mut tmp = sf.eval_fast(model.atoms[a].xs(), model.atoms[b].xs(), r2);
                curl(&mut tmp, v);
                e += tmp;
            }
        }
    }
    e
}

/// The raw `(intermolecular, intramolecular)` energies of the current pose,
/// before the conf-independent term.
pub fn score_components(model: &Model, sf: &ScoringFunction, gbox: &SearchBox) -> (f64, f64) {
    (
        eval_inter(model, sf, gbox, AUTHENTIC_V),
        eval_intra(model, sf, AUTHENTIC_V),
    )
}

/// Assemble the reported breakdown from raw components and an `unbound`
/// reference: `EFE = conf_independent(num_tors, inter + intra - unbound)`.
/// For `--score_only` the unbound reference is the pose's own intra; for a
/// dock it is the best pose's intra.
pub fn assemble_breakdown(
    inter: f64,
    intra: f64,
    unbound: f64,
    num_tors: f64,
    sf: &ScoringFunction,
) -> ScoreBreakdown {
    let base = inter + intra - unbound;
    let total = sf.conf_independent(num_tors, base);
    let torsional = total - base;
    ScoreBreakdown {
        estimated_free_energy: total,
        intermolecular: inter,
        internal: intra,
        torsional,
        unbound,
    }
}

/// Compute the full `--score_only` breakdown for a rigid-receptor Vina dock
/// (the unbound reference is the input pose's own intramolecular energy).
pub fn score(model: &Model, sf: &ScoringFunction, gbox: &SearchBox) -> ScoreBreakdown {
    let (inter, intra) = score_components(model, sf, gbox);
    assemble_breakdown(inter, intra, intra, model.num_tors(), sf)
}

/// Convenience: a `[f64; 3]` from a [`Vec3`].
pub fn vec3_to_array(v: Vec3) -> [f64; 3] {
    [v[0], v[1], v[2]]
}

// ---------------------------------------------------------------------------
// Derivatives (forces + conformation-space gradient)
// ---------------------------------------------------------------------------

/// A per-pair derivative provider: given XS types and squared distance, returns
/// `(energy, derivative-over-r)`. Two implementations: the engine's interpolated
/// path ([`pair_deriv_interp`]) and the exact analytic path ([`pair_deriv_exact`]).
pub type PairDeriv<'a> = dyn Fn(usize, usize, f64) -> (f64, f64) + 'a;

/// The engine's interpolated `(energy, dor)` — what the grid uses.
pub fn pair_deriv_interp(sf: &ScoringFunction) -> impl Fn(usize, usize, f64) -> (f64, f64) + '_ {
    move |t1, t2, r2| sf.eval_deriv(t1, t2, r2)
}

/// The exact analytic `(energy, dor)` where `dor = (dE/dr) / r` — energy and
/// derivative are consistent by construction (for finite-difference validation).
pub fn pair_deriv_exact(sf: &ScoringFunction) -> impl Fn(usize, usize, f64) -> (f64, f64) + '_ {
    move |t1, t2, r2| {
        let r = r2.sqrt();
        let (e, dedr) = sf.eval_exact_deriv(t1, t2, r);
        let dor = if r < crate::math::EPSILON_FL {
            0.0
        } else {
            dedr / r
        };
        (e, dor)
    }
}

/// Intermolecular energy and per-atom forces, using the supplied per-pair
/// derivative provider. Writes `forces[i]` for every movable atom.
fn eval_inter_deriv_with(
    model: &Model,
    sf: &ScoringFunction,
    gbox: &SearchBox,
    v: f64,
    pair: &PairDeriv,
    forces: &mut [Vec3],
) -> f64 {
    let cutoff_sqr = sf.cutoff_sqr();
    let mut e = 0.0;
    for i in 0..model.num_movable_atoms() {
        let a = &model.atoms[i];
        let t1 = a.xs();
        if t1 >= XS_TYPE_SIZE || is_glue_type(t1) {
            forces[i] = ZERO;
            continue;
        }
        let a_coords = model.coords[i];

        let mut adjusted = a_coords;
        let mut out_of_bounds_deriv = ZERO;
        let mut penalty = 0.0;
        for j in 0..3 {
            if a_coords[j] < gbox.begin[j] {
                adjusted[j] = gbox.begin[j];
                out_of_bounds_deriv[j] = -1.0;
                penalty += (a_coords[j] - gbox.begin[j]).abs();
            } else if a_coords[j] > gbox.end[j] {
                adjusted[j] = gbox.end[j];
                out_of_bounds_deriv[j] = 1.0;
                penalty += (a_coords[j] - gbox.end[j]).abs();
            }
        }
        penalty *= gbox.slope;
        out_of_bounds_deriv = gbox.slope * out_of_bounds_deriv;

        let mut this_e = 0.0;
        let mut deriv = ZERO;
        for b in &model.grid_atoms {
            let t2 = b.xs();
            if t2 >= XS_TYPE_SIZE {
                continue;
            }
            let r_ba = adjusted - b.coords;
            let r2 = r_ba.norm_sqr();
            if r2 < cutoff_sqr {
                let (e_pair, dor) = pair(a.xs(), b.xs(), r2);
                this_e += e_pair;
                deriv += dor * r_ba;
            }
        }
        curl_deriv(&mut this_e, &mut deriv, v);
        forces[i] = deriv + out_of_bounds_deriv;
        e += this_e + penalty;
    }
    e
}

/// Intramolecular energy and forces over the ligand interacting pairs.
fn eval_intra_deriv_with(
    model: &Model,
    sf: &ScoringFunction,
    v: f64,
    pair: &PairDeriv,
    forces: &mut [Vec3],
) -> f64 {
    let cutoff_sqr = sf.cutoff_sqr();
    let mut e = 0.0;
    for lig in &model.ligands {
        for &(a, b) in &lig.pairs {
            let r = model.coords[b] - model.coords[a]; // a -> b
            let r2 = r.norm_sqr();
            if r2 < cutoff_sqr {
                let (mut e_pair, dor) = pair(model.atoms[a].xs(), model.atoms[b].xs(), r2);
                let mut force = dor * r;
                curl_deriv(&mut e_pair, &mut force, v);
                e += e_pair;
                forces[a] -= force;
                forces[b] += force;
            }
        }
    }
    e
}

/// Total energy with a chosen per-pair derivative provider, plus the
/// conformation-space gradient `g`. Returns the energy.
pub fn eval_deriv_with(
    model: &Model,
    sf: &ScoringFunction,
    gbox: &SearchBox,
    v: f64,
    pair: &PairDeriv,
    g: &mut Change,
) -> f64 {
    let mut forces = vec![ZERO; model.atoms.len()];
    let mut e = eval_inter_deriv_with(model, sf, gbox, v, pair, &mut forces);
    e += eval_intra_deriv_with(model, sf, v, pair, &mut forces);
    for (i, lig) in model.ligands.iter().enumerate() {
        lig.body
            .derivative(&model.coords, &forces, &mut g.ligands[i]);
    }
    e
}

/// Total energy (interpolated) plus the conformation-space gradient `g`, using
/// direct pairwise summation for the intermolecular forces. This is the path
/// the optimizer follows.
pub fn eval_deriv(
    model: &Model,
    sf: &ScoringFunction,
    gbox: &SearchBox,
    v: f64,
    g: &mut Change,
) -> f64 {
    eval_deriv_with(model, sf, gbox, v, &pair_deriv_interp(sf), g)
}

/// BFGS local optimization of `conf` (updated in place) against the chosen
/// per-pair derivative provider: `set` the model from the conformation,
/// evaluate energy + gradient, step. Returns `(final_energy, eval_count)`.
pub fn local_optimize(
    model: &mut Model,
    sf: &ScoringFunction,
    gbox: &SearchBox,
    v: f64,
    pair: &PairDeriv,
    conf: &mut Conf,
    max_steps: usize,
) -> (f64, i64) {
    let size = model.get_size();
    let mut g = Change::new(&size);
    let mut evalcount: i64 = 0;
    let res = {
        let mut f = |c: &Conf, grad: &mut Change| {
            model.set(c);
            eval_deriv_with(model, sf, gbox, v, pair, grad)
        };
        crate::optimize::bfgs(&mut f, conf, &mut g, max_steps, &mut evalcount)
    };
    model.set(conf);
    (res, evalcount)
}

/// `max_steps = (25 + num_movable_atoms) / 3` (the Vina default).
pub fn default_max_steps(model: &Model) -> usize {
    (25 + model.num_movable_atoms()) / 3
}

/// Total energy using the precomputed grid [`Cache`](super::grid::Cache) for the intermolecular
/// forces. The intramolecular term still uses the interpolated pairwise
/// potential.
pub fn eval_deriv_cache(
    model: &Model,
    sf: &ScoringFunction,
    cache: &super::grid::Cache,
    v: f64,
    g: &mut Change,
) -> f64 {
    let mut forces = vec![ZERO; model.atoms.len()];
    let mut e = cache.eval_deriv(model, v, &mut forces);
    e += eval_intra_deriv_with(model, sf, v, &pair_deriv_interp(sf), &mut forces);
    for (i, lig) in model.ligands.iter().enumerate() {
        lig.body
            .derivative(&model.coords, &forces, &mut g.ligands[i]);
    }
    e
}

/// BFGS local optimization against the grid cache.
pub fn local_optimize_cache(
    model: &mut Model,
    sf: &ScoringFunction,
    cache: &super::grid::Cache,
    v: f64,
    conf: &mut Conf,
    max_steps: usize,
) -> (f64, i64) {
    let size = model.get_size();
    let mut g = Change::new(&size);
    let mut evalcount: i64 = 0;
    let res = {
        let mut f = |c: &Conf, grad: &mut Change| {
            model.set(c);
            eval_deriv_cache(model, sf, cache, v, grad)
        };
        crate::optimize::bfgs(&mut f, conf, &mut g, max_steps, &mut evalcount)
    };
    model.set(conf);
    (res, evalcount)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn box_quantizes_like_vina() {
        // size 20 with granularity 0.375: n_voxels = ceil(53.33) = 54,
        // real_span = 0.375*54 = 20.25, so the box is slightly larger than 20.
        let b = SearchBox::from_center_size([0.0, 0.0, 0.0], [20.0, 20.0, 20.0]);
        let real_span = 0.375 * 54.0;
        assert!((b.end[0] - b.begin[0] - real_span).abs() < 1e-12);
        assert!((b.begin[0] - (-real_span / 2.0)).abs() < 1e-12);

        // A non-divisible size: 22.3 -> ceil(22.3/0.375)=ceil(59.47)=60 ->
        // real_span = 22.5 (differs from the naive 22.3 box).
        let b = SearchBox::new([1.0, 2.0, 3.0], [22.3, 22.3, 22.3], 0.375, false);
        let span = b.end[1] - b.begin[1];
        assert!((span - 22.5).abs() < 1e-12, "span = {span}");
        assert!((b.begin[1] - (2.0 - 22.5 / 2.0)).abs() < 1e-12);
    }

    #[test]
    fn force_even_voxels_rounds_up_odd() {
        // size that yields an odd voxel count: 19.0/0.5 = 38 (even already);
        // use 0.4 -> 19/0.4 = 47.5 -> ceil 48 (even). Pick 18.6/0.4 = 46.5 ->
        // ceil 47 (odd) -> forced to 48.
        let odd = SearchBox::new([0.0; 3], [18.6; 3], 0.4, false);
        let even = SearchBox::new([0.0; 3], [18.6; 3], 0.4, true);
        assert!((odd.end[0] - odd.begin[0] - 0.4 * 47.0).abs() < 1e-12);
        assert!((even.end[0] - even.begin[0] - 0.4 * 48.0).abs() < 1e-12);
    }
}
