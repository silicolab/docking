// SPDX-License-Identifier: Apache-2.0
//! Self-tests for the ground-truth diff harness.
//!
//! These prove the harness correctly reads the archived official-Vina golden
//! artifacts and that its math (RMSD, energy-identity) is sound. Engine output
//! is compared against these same parsed golden values.

mod common;
use common::*;

#[test]
fn score_only_breakdown_parses_with_expected_values() {
    let s = ScoreBreakdown::load_golden("1iep", "score_only.out.txt");
    // Values printed by official Vina v1.2.7 for the 1iep input pose.
    assert_close("estimated_fe", s.estimated_free_energy, -12.513, 1e-3);
    assert_close("intermolecular", s.intermolecular, -17.634, 1e-3);
    assert_close("internal", s.internal, -0.485, 1e-3);
    assert_close("torsional", s.torsional, 5.121, 1e-3);
    assert_close("unbound", s.unbound, -0.485, 1e-3);
}

#[test]
fn score_breakdown_satisfies_vina_energy_identity() {
    // Vina reports EFE = (1) + (2) + (3) - (4).
    let s = ScoreBreakdown::load_golden("1iep", "score_only.out.txt");
    let recomputed = s.intermolecular + s.internal + s.torsional - s.unbound;
    assert_close("EFE identity", recomputed, s.estimated_free_energy, 1e-3);
}

#[test]
fn local_only_breakdown_parses() {
    let s = ScoreBreakdown::load_golden("1iep", "local_only.out.txt");
    assert_close("local estimated_fe", s.estimated_free_energy, -13.170, 1e-3);
    assert_close("local intermolecular", s.intermolecular, -18.559, 1e-3);
}

#[test]
fn dock_table_parses_ranked_modes() {
    let text = std::fs::read_to_string(golden_path("1iep", "dock_seed42.out.txt")).unwrap();
    let modes = parse_dock_table(&text);
    assert_eq!(modes.len(), 9, "expected 9 ranked modes in stdout table");
    assert_eq!(modes[0].rank, 1);
    assert_close("best affinity", modes[0].affinity, -13.2, 0.05);
    // Rank 1 is the reference for rmsd, so its rmsd columns are zero.
    assert_eq!(modes[0].rmsd_lb, 0.0);
    assert_eq!(modes[0].rmsd_ub, 0.0);
    // Affinities are sorted ascending in energy (best/most-negative first).
    for w in modes.windows(2) {
        assert!(
            w[0].affinity <= w[1].affinity + 1e-9,
            "modes not sorted by affinity"
        );
    }
}

#[test]
fn dock_poses_parse_with_results_and_atom_counts() {
    let poses = load_golden_poses("1iep", "dock_seed42.pdbqt");
    // Output file keeps only modes within the default 3 kcal/mol energy range.
    assert_eq!(poses.len(), 4, "expected 4 written poses");

    let best = &poses[0];
    let result = best
        .vina_result
        .expect("first pose has a VINA RESULT remark");
    assert_close("best pose affinity", result[0], -13.201, 1e-3);

    // 1iep ligand: 40 atoms total, 3 polar hydrogens (HD), rest heavy.
    assert_eq!(best.coords.len(), 40);
    assert_eq!(best.elements.iter().filter(|e| *e == "HD").count(), 3);
    assert_eq!(best.heavy_atom_indices().count(), 37);

    // Every written pose has the same atom count.
    for (i, p) in poses.iter().enumerate() {
        assert_eq!(p.coords.len(), 40, "pose {i} atom count");
    }
}

#[test]
fn local_only_pose_parses_and_differs_from_input() {
    // `--local_only` writes only SMILES remarks (no VINA RESULT); its energy is
    // reported on stdout. Here we check the optimized pose itself: it parses with
    // the full atom count and BFGS actually moved atoms away from the input pose.
    let breakdown = ScoreBreakdown::load_golden("1iep", "local_only.out.txt");
    assert!(
        breakdown.estimated_free_energy < 0.0,
        "expected a favorable energy"
    );

    let opt = load_golden_poses("1iep", "local_only.pdbqt");
    assert_eq!(opt.len(), 1, "local_only writes a single pose");
    assert_eq!(opt[0].coords.len(), 40);
    assert!(
        opt[0].vina_result.is_none(),
        "local_only has no VINA RESULT remark"
    );

    let input = parse_poses(&std::fs::read_to_string(golden_path("1iep", "ligand.pdbqt")).unwrap());
    assert_eq!(input.len(), 1);
    let moved = rmsd(&opt[0].coords, &input[0].coords);
    assert!(
        moved > 1e-3,
        "optimized pose should differ from input (rmsd={moved:.4} A)"
    );
}

#[test]
fn rmsd_is_zero_for_identical_and_correct_for_offset() {
    let a = [[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-1.0, 5.0, 2.0]];
    assert_eq!(rmsd(&a, &a), 0.0);

    // Shift every atom by (3,4,0): each squared displacement is 25, RMSD = 5.
    let b: Vec<[f64; 3]> = a.iter().map(|p| [p[0] + 3.0, p[1] + 4.0, p[2]]).collect();
    assert_close("offset rmsd", rmsd(&a, &b), 5.0, 1e-12);
}

#[test]
fn pose_self_rmsd_against_golden_is_zero() {
    // Sanity: a golden pose compared to itself is exactly zero, and the parser
    // produced usable coordinates.
    let poses = load_golden_poses("1iep", "dock_seed42.pdbqt");
    let coords = &poses[0].coords;
    assert!(coords.len() >= 37);
    assert_eq!(rmsd(coords, coords), 0.0);
}
