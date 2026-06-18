// SPDX-License-Identifier: Apache-2.0
//! Ground-truth diff harness shared by integration tests.
//!
//! These helpers parse the artifacts produced by official AutoDock Vina v1.2.7
//! (see `tests/golden/`) and compare them against values produced by this crate
//! within explicit tolerances. They are deliberately
//! independent of the engine's own parser so that a bug in the engine cannot
//! also corrupt the reference it is checked against.
//!
//! This module is compiled into each integration-test binary that declares
//! `mod common;`. Some helpers are not yet used by any test;
//! `#[allow(dead_code)]` keeps the harness complete without warnings meanwhile.
#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// Absolute path to `tests/golden/<system>/<file>`.
pub fn golden_path(system: &str, file: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("golden")
        .join(system)
        .join(file)
}

/// Absolute path to `tests/fixtures/<rel>`.
pub fn fixture_path(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(rel)
}

// ---------------------------------------------------------------------------
// Energy breakdown (from `--score_only` / `--local_only` stdout)
// ---------------------------------------------------------------------------

/// The energy breakdown Vina prints for a single pose. Fields map to the
/// numbered lines Vina writes, e.g.:
///
/// ```text
/// Estimated Free Energy of Binding   : -12.513 (kcal/mol) [=(1)+(2)+(3)-(4)]
/// (1) Final Intermolecular Energy    : -17.634 (kcal/mol)
/// (2) Final Total Internal Energy    : -0.485 (kcal/mol)
/// (3) Torsional Free Energy          : 5.121 (kcal/mol)
/// (4) Unbound System's Energy        : -0.485 (kcal/mol)
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScoreBreakdown {
    pub estimated_free_energy: f64,
    pub intermolecular: f64,
    pub internal: f64,
    pub torsional: f64,
    pub unbound: f64,
}

impl ScoreBreakdown {
    /// Parse the breakdown out of Vina's stdout text.
    pub fn parse(text: &str) -> Self {
        let find = |needle: &str| -> f64 {
            let line = text
                .lines()
                .find(|l| l.contains(needle))
                .unwrap_or_else(|| panic!("no line containing {needle:?} in vina output"));
            parse_first_float(line)
                .unwrap_or_else(|| panic!("no float on line for {needle:?}: {line:?}"))
        };
        ScoreBreakdown {
            estimated_free_energy: find("Estimated Free Energy of Binding"),
            intermolecular: find("(1) Final Intermolecular Energy"),
            internal: find("(2) Final Total Internal Energy"),
            torsional: find("(3) Torsional Free Energy"),
            unbound: find("(4) Unbound System"),
        }
    }

    /// Load and parse a golden `*.out.txt` for a system.
    pub fn load_golden(system: &str, file: &str) -> Self {
        let p = golden_path(system, file);
        let text =
            std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()));
        Self::parse(&text)
    }
}

// ---------------------------------------------------------------------------
// Dock result table (from full-dock stdout)
// ---------------------------------------------------------------------------

/// One row of Vina's ranked result table.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DockMode {
    pub rank: usize,
    pub affinity: f64,
    pub rmsd_lb: f64,
    pub rmsd_ub: f64,
}

/// Parse the ranked mode table from a full-dock stdout, e.g. rows like
/// `   1        -13.2          0          0`.
pub fn parse_dock_table(text: &str) -> Vec<DockMode> {
    let mut modes = Vec::new();
    for line in text.lines() {
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() != 4 {
            continue;
        }
        // A data row is "<int> <float> <float> <float>".
        let (Ok(rank), Ok(affinity), Ok(rmsd_lb), Ok(rmsd_ub)) = (
            cols[0].parse::<usize>(),
            cols[1].parse::<f64>(),
            cols[2].parse::<f64>(),
            cols[3].parse::<f64>(),
        ) else {
            continue;
        };
        modes.push(DockMode {
            rank,
            affinity,
            rmsd_lb,
            rmsd_ub,
        });
    }
    modes
}

// ---------------------------------------------------------------------------
// Poses (from PDBQT output, possibly multi-model)
// ---------------------------------------------------------------------------

/// A single pose parsed from a PDBQT model: the `REMARK VINA RESULT` triple
/// (affinity, rmsd l.b., rmsd u.b.) if present, plus per-atom data.
#[derive(Debug, Clone, Default)]
pub struct Pose {
    pub vina_result: Option<[f64; 3]>,
    pub atom_names: Vec<String>,
    pub elements: Vec<String>,
    pub coords: Vec<[f64; 3]>,
}

impl Pose {
    pub fn heavy_atom_indices(&self) -> impl Iterator<Item = usize> + '_ {
        // PDBQT autodock types: hydrogens are "H" (nonpolar, usually merged) and
        // "HD" (polar). Heavy atoms are everything else.
        self.elements
            .iter()
            .enumerate()
            .filter(|(_, e)| *e != "H" && *e != "HD")
            .map(|(i, _)| i)
    }
}

/// Parse all poses (MODEL...ENDMDL blocks, or a single implicit model) from
/// PDBQT text. Coordinates are read from fixed PDB columns.
pub fn parse_poses(text: &str) -> Vec<Pose> {
    let mut poses = Vec::new();
    let mut cur = Pose::default();
    let mut started = false;

    let flush = |cur: &mut Pose, poses: &mut Vec<Pose>| {
        if !cur.coords.is_empty() || cur.vina_result.is_some() {
            poses.push(std::mem::take(cur));
        }
    };

    for line in text.lines() {
        if line.starts_with("MODEL") {
            flush(&mut cur, &mut poses);
            started = true;
        } else if line.starts_with("ENDMDL") {
            flush(&mut cur, &mut poses);
        } else if line.starts_with("REMARK VINA RESULT") {
            let nums = parse_floats(line);
            if nums.len() >= 3 {
                cur.vina_result = Some([nums[0], nums[1], nums[2]]);
            }
        } else if line.starts_with("ATOM") || line.starts_with("HETATM") {
            started = true;
            let x = col_f64(line, 30, 38);
            let y = col_f64(line, 38, 46);
            let z = col_f64(line, 46, 54);
            cur.coords.push([x, y, z]);
            cur.atom_names.push(col_str(line, 12, 16));
            // PDBQT autodock atom type lives in the last whitespace-separated token.
            let element = line.split_whitespace().last().unwrap_or("").to_string();
            cur.elements.push(element);
        }
    }
    // Single-model files have no MODEL/ENDMDL.
    if !started {
        return poses;
    }
    flush(&mut cur, &mut poses);
    poses
}

/// Load and parse a golden PDBQT file for a system.
pub fn load_golden_poses(system: &str, file: &str) -> Vec<Pose> {
    let p = golden_path(system, file);
    let text = std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()));
    parse_poses(&text)
}

// ---------------------------------------------------------------------------
// RMSD
// ---------------------------------------------------------------------------

/// Coordinate (no-superposition) RMSD between two equal-length atom sets, in the
/// same atom order. This matches Vina's "rmsd l.b." style of comparison on poses
/// already in the receptor frame.
pub fn rmsd(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    assert_eq!(
        a.len(),
        b.len(),
        "rmsd: mismatched atom counts {} vs {}",
        a.len(),
        b.len()
    );
    assert!(!a.is_empty(), "rmsd: empty coordinate set");
    let mut sum = 0.0;
    for (p, q) in a.iter().zip(b.iter()) {
        let dx = p[0] - q[0];
        let dy = p[1] - q[1];
        let dz = p[2] - q[2];
        sum += dx * dx + dy * dy + dz * dz;
    }
    (sum / a.len() as f64).sqrt()
}

// ---------------------------------------------------------------------------
// Tolerant comparison
// ---------------------------------------------------------------------------

/// Initial tolerances. Centralized so changes are deliberate.
pub mod tol {
    /// Each reported energy component of `--score_only` / `--local_only`, kcal/mol.
    pub const SCORE_COMPONENT_KCAL: f64 = 1e-3;
    /// Final energy of a local optimization, kcal/mol.
    pub const LOCAL_ENERGY_KCAL: f64 = 1e-2;
    /// Per-atom coordinate after local optimization, Angstrom.
    pub const LOCAL_COORD_ANGSTROM: f64 = 1e-2;
    /// Best-mode affinity of a full dock, kcal/mol.
    pub const DOCK_BEST_AFFINITY_KCAL: f64 = 0.3;
    /// Relative tolerance for analytic gradient vs. finite differences.
    pub const GRADIENT_REL: f64 = 1e-5;
}

/// Assert `got` is within `tol` of `want` (absolute), with a descriptive message.
#[track_caller]
pub fn assert_close(label: &str, got: f64, want: f64, tol: f64) {
    let d = (got - want).abs();
    assert!(
        d <= tol,
        "{label}: |{got} - {want}| = {d:.3e} exceeds tolerance {tol:.3e}"
    );
}

// ---------------------------------------------------------------------------
// Small parsing helpers
// ---------------------------------------------------------------------------

fn col_str(line: &str, start: usize, end: usize) -> String {
    line.get(start..end).unwrap_or("").trim().to_string()
}

fn col_f64(line: &str, start: usize, end: usize) -> f64 {
    let s = line.get(start..end).unwrap_or("").trim();
    s.parse::<f64>()
        .unwrap_or_else(|_| panic!("bad float in columns {start}..{end}: {s:?} (line {line:?})"))
}

/// First float anywhere on the line (after the first ':' if present, to skip
/// numbers embedded in labels like "(1)").
fn parse_first_float(line: &str) -> Option<f64> {
    let scan = line.split_once(':').map(|(_, r)| r).unwrap_or(line);
    parse_floats(scan).into_iter().next()
}

/// Every float-looking token on the line.
fn parse_floats(line: &str) -> Vec<f64> {
    line.split(|c: char| {
        !(c.is_ascii_digit() || c == '.' || c == '-' || c == '+' || c == 'e' || c == 'E')
    })
    .filter_map(|tok| {
        if tok.is_empty() || tok == "-" || tok == "+" || tok == "." {
            None
        } else {
            tok.parse::<f64>().ok()
        }
    })
    .collect()
}
