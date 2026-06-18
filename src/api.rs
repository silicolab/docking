// SPDX-License-Identifier: Apache-2.0
//! High-level docking API — the stable entry point for downstream consumers.
//!
//! Given an already-prepared receptor and ligand (PDBQT) plus a search box, it
//! produces ranked docked poses with Vina scores, or a single-point score for an
//! input pose. This wraps parsing, grid construction, the Monte-Carlo search, and
//! scoring behind a small, stable surface.
//!
//! ```no_run
//! use docking::api::{dock, DockConfig};
//!
//! let receptor = std::fs::read_to_string("receptor.pdbqt")?;
//! let ligand = std::fs::read_to_string("ligand.pdbqt")?;
//! let poses = dock(
//!     &receptor,
//!     &ligand,
//!     [15.190, 53.903, 16.917], // box center (Å)
//!     [20.0, 20.0, 20.0],       // box size (Å)
//!     &DockConfig { seed: 42, ..Default::default() },
//! )?;
//! for (rank, pose) in poses.iter().enumerate() {
//!     println!("mode {}: {:.3} kcal/mol", rank + 1, pose.affinity);
//!     // `pose.pdbqt` is the docked conformation, ready to write to a file.
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::atom::AtomTyping;
use crate::pdbqt::{parse_ligand_pdbqt_from_string, parse_receptor_pdbqt_from_string, PdbqtError};
use crate::scoring::eval::{score, SearchBox};
use crate::scoring::grid::Cache;
use crate::scoring::ScoringFunction;
use crate::search::dock as run_dock;

/// Errors from the docking API.
#[derive(Debug)]
#[non_exhaustive]
pub enum Error {
    /// A PDBQT input could not be parsed.
    Pdbqt(PdbqtError),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Pdbqt(e) => write!(f, "PDBQT parse error: {e}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<PdbqtError> for Error {
    fn from(e: PdbqtError) -> Self {
        Error::Pdbqt(e)
    }
}

/// Search configuration. `Default` matches the AutoDock Vina CLI defaults.
#[derive(Debug, Clone)]
pub struct DockConfig {
    /// Number of independent Monte-Carlo runs (search thoroughness).
    pub exhaustiveness: usize,
    /// Maximum number of binding modes to return.
    pub num_modes: usize,
    /// Minimum RMSD between reported modes (Å).
    pub min_rmsd: f64,
    /// Random seed (the search is deterministic for a fixed seed).
    pub seed: u32,
    /// Override the per-run Monte-Carlo step count. `None` uses Vina's heuristic
    /// (`70*3*(50 + num_movable_atoms + 10*num_dof)/2`); set a smaller value for
    /// a faster, lighter search.
    pub max_global_steps: Option<usize>,
}

impl Default for DockConfig {
    fn default() -> Self {
        DockConfig {
            exhaustiveness: 8,
            num_modes: 9,
            min_rmsd: 1.0,
            seed: 0,
            max_global_steps: None,
        }
    }
}

/// A docked pose: its Vina energies and the conformation as PDBQT text.
#[derive(Debug, Clone)]
pub struct Pose {
    /// Estimated free energy of binding (kcal/mol) — the headline affinity.
    pub affinity: f64,
    /// Final intermolecular energy (kcal/mol).
    pub intermolecular: f64,
    /// Final total internal energy (kcal/mol).
    pub internal: f64,
    /// Torsional free energy penalty (kcal/mol).
    pub torsional: f64,
    /// The pose as a PDBQT block (a `REMARK VINA RESULT` line plus the ligand).
    pub pdbqt: String,
}

/// Dock a ligand into a receptor within the given box; returns ranked poses
/// (best affinity first).
pub fn dock(
    receptor_pdbqt: &str,
    ligand_pdbqt: &str,
    center: [f64; 3],
    size: [f64; 3],
    config: &DockConfig,
) -> Result<Vec<Pose>, Error> {
    let mut model = parse_ligand_pdbqt_from_string(ligand_pdbqt, AtomTyping::Xs)?;
    let receptor = parse_receptor_pdbqt_from_string(receptor_pdbqt, AtomTyping::Xs)?;
    model.attach_receptor(receptor);

    let gbox = SearchBox::from_center_size(center, size);
    let sf = ScoringFunction::vina();
    let cache = Cache::populate(&model, &sf, &gbox);

    let docked = run_dock(
        &model,
        &cache,
        &sf,
        &gbox,
        config.exhaustiveness,
        config.num_modes,
        config.min_rmsd,
        config.seed,
        config.max_global_steps,
    );

    // Render each pose to PDBQT, with RMSD reported relative to the best mode.
    let mut out_model = model.clone();
    let els = out_model.heavy_atom_movable_els();
    let mut best_coords: Vec<crate::math::Vec3> = Vec::new();
    let mut poses = Vec::with_capacity(docked.len());
    for (i, d) in docked.iter().enumerate() {
        out_model.set(&d.conf);
        let coords = out_model.heavy_atom_movable_coords();
        if i == 0 {
            best_coords = coords.clone();
        }
        // Vina reports two RMSDs to the best mode: a symmetry-aware lower bound
        // and a same-order upper bound (both zero for the best mode itself).
        let (rmsd_lb, rmsd_ub) = if i == 0 {
            (0.0, 0.0)
        } else {
            (
                crate::search::output::rmsd_lower_bound(&els, &coords, &best_coords),
                crate::search::output::rmsd_upper_bound(&coords, &best_coords),
            )
        };
        let remark = format!(
            "REMARK VINA RESULT:  {:8.3} {:8.3} {:8.3}\n",
            d.breakdown.estimated_free_energy, rmsd_lb, rmsd_ub
        );
        poses.push(Pose {
            affinity: d.breakdown.estimated_free_energy,
            intermolecular: d.breakdown.intermolecular,
            internal: d.breakdown.internal,
            torsional: d.breakdown.torsional,
            pdbqt: out_model.write_ligand_pdbqt(0, &remark),
        });
    }
    Ok(poses)
}

/// Single-point score of the ligand's input pose (no search) — the
/// `--score_only` breakdown.
pub fn score_only(
    receptor_pdbqt: &str,
    ligand_pdbqt: &str,
    center: [f64; 3],
    size: [f64; 3],
) -> Result<crate::scoring::eval::ScoreBreakdown, Error> {
    let mut model = parse_ligand_pdbqt_from_string(ligand_pdbqt, AtomTyping::Xs)?;
    let receptor = parse_receptor_pdbqt_from_string(receptor_pdbqt, AtomTyping::Xs)?;
    model.attach_receptor(receptor);
    let gbox = SearchBox::from_center_size(center, size);
    let sf = ScoringFunction::vina();
    Ok(score(&model, &sf, &gbox))
}
