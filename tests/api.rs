// SPDX-License-Identifier: Apache-2.0
//! Exercise the public `docking::api` surface that downstream consumers call.

mod common;
use common::*;

use docking::api::{dock, score_only, DockConfig};

fn read(p: std::path::PathBuf) -> String {
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

const CENTER: [f64; 3] = [15.190, 53.903, 16.917];
const SIZE: [f64; 3] = [20.0, 20.0, 20.0];

#[test]
fn api_score_only_matches_golden() {
    let receptor = read(golden_path("1iep", "receptor.pdbqt"));
    let ligand = read(golden_path("1iep", "ligand.pdbqt"));
    let got = score_only(&receptor, &ligand, CENTER, SIZE).expect("score_only");
    let want = ScoreBreakdown::load_golden("1iep", "score_only.out.txt");
    assert_close(
        "affinity",
        got.estimated_free_energy,
        want.estimated_free_energy,
        1e-3,
    );
    assert_close(
        "intermolecular",
        got.intermolecular,
        want.intermolecular,
        1e-3,
    );
}

#[test]
fn api_dock_returns_ranked_poses() {
    let receptor = read(golden_path("1iep", "receptor.pdbqt"));
    let ligand = read(golden_path("1iep", "ligand.pdbqt"));

    // Light config — this test validates the API plumbing (ranking, PDBQT
    // emission); the strict convergence-to-Vina check lives in the docking tests.
    let config = DockConfig {
        exhaustiveness: 3,
        num_modes: 9,
        min_rmsd: 1.0,
        seed: 42,
        max_global_steps: Some(120),
    };
    let poses = dock(&receptor, &ligand, CENTER, SIZE, &config).expect("dock");

    assert!(!poses.is_empty(), "no poses returned");
    assert!(poses.len() <= 9);

    // Ranked best-first.
    for w in poses.windows(2) {
        assert!(
            w[0].affinity <= w[1].affinity + 1e-9,
            "poses not ranked by affinity"
        );
    }

    // The search reaches a credible binding pose.
    eprintln!("api dock best affinity = {:.3}", poses[0].affinity);
    assert!(
        poses[0].affinity < -10.0,
        "best affinity {} not a binding pose",
        poses[0].affinity
    );

    // The emitted PDBQT is well-formed: a VINA RESULT remark + the right atom count.
    let best_pdbqt = &poses[0].pdbqt;
    assert!(best_pdbqt.contains("REMARK VINA RESULT:"));
    let parsed = parse_poses(best_pdbqt);
    assert_eq!(parsed.len(), 1);
    assert_eq!(parsed[0].coords.len(), 40);
    let result = parsed[0].vina_result.expect("remark parsed");
    assert_close("pdbqt remark affinity", result[0], poses[0].affinity, 1e-3);
}
