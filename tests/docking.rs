// SPDX-License-Identifier: Apache-2.0
//! Monte-Carlo global search (end-to-end docking) against official
//! Vina's `--seed 42 --exhaustiveness 8` golden result.
//!
//! The search uses a different RNG than Boost, so the trajectory is not
//! bit-identical; we validate **convergence**:
//! the search must find the same global minimum (binding pose) as Vina, with the
//! top affinity within tolerance and a pose close to Vina's best mode. A reduced
//! step count is used to keep the test fast — it already converges to the global
//! minimum (full exhaustiveness only adds redundant sampling).

mod common;
use common::*;

use docking::atom::AtomTyping;
use docking::pdbqt::{parse_ligand_pdbqt_from_string, parse_receptor_pdbqt_from_string};
use docking::scoring::eval::SearchBox;
use docking::scoring::grid::Cache;
use docking::scoring::ScoringFunction;
use docking::search::dock;

fn read(p: std::path::PathBuf) -> String {
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

#[test]
fn docking_converges_to_vina_global_minimum() {
    let lig = read(golden_path("1iep", "ligand.pdbqt"));
    let rec = read(golden_path("1iep", "receptor.pdbqt"));
    let mut m = parse_ligand_pdbqt_from_string(&lig, AtomTyping::Xs).unwrap();
    m.attach_receptor(parse_receptor_pdbqt_from_string(&rec, AtomTyping::Xs).unwrap());
    let gbox = SearchBox::from_center_size([15.190, 53.903, 16.917], [20.0, 20.0, 20.0]);
    let sf = ScoringFunction::vina();
    let cache = Cache::populate(&m, &sf, &gbox);

    // Reduced steps (vs Vina's ~23100) — already converges to the global min
    // and stays light enough to run in a debug build.
    let poses = dock(&m, &cache, &sf, &gbox, 4, 9, 1.0, 42, Some(150));
    assert!(!poses.is_empty(), "search produced no poses");

    // Vina's golden best mode (mode 1) affinity for seed 42.
    let golden_table = parse_dock_table(&read(golden_path("1iep", "dock_seed42.out.txt")));
    let vina_best = golden_table[0].affinity; // -13.2
    let my_best = poses[0].breakdown.estimated_free_energy;
    eprintln!("dock: my best EFE = {my_best:.3}, Vina best = {vina_best:.3}");

    // Top affinity within tolerance of Vina's.
    assert_close(
        "best affinity",
        my_best,
        vina_best,
        tol::DOCK_BEST_AFFINITY_KCAL,
    );

    for w in poses.windows(2) {
        assert!(
            w[0].breakdown.estimated_free_energy <= w[1].breakdown.estimated_free_energy + 1e-9,
            "poses not sorted by affinity"
        );
    }

    // The best pose should match Vina's best mode (the binding pose) closely.
    let mut model = m.clone();
    model.set(&poses[0].conf);
    let my_out = model.write_ligand_pdbqt(0, "");
    let my_pose = parse_poses(&my_out);
    let golden_best = load_golden_poses("1iep", "dock_seed42.pdbqt");
    let pose_rmsd = rmsd(&my_pose[0].coords, &golden_best[0].coords);
    eprintln!("dock: best-pose RMSD vs Vina's best mode = {pose_rmsd:.3} A");
    assert!(
        pose_rmsd < 2.0,
        "best pose RMSD {pose_rmsd:.3} A not close to Vina's best mode"
    );
}
