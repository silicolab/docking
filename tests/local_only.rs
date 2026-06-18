// SPDX-License-Identifier: Apache-2.0
//! Reproduce official Vina `--local_only` (grid-based BFGS local
//! optimization of the input pose), comparing the optimized energy and pose.

mod common;
use common::*;

use docking::atom::AtomTyping;
use docking::pdbqt::{parse_ligand_pdbqt_from_string, parse_receptor_pdbqt_from_string};
use docking::scoring::eval::{default_max_steps, local_optimize_cache, score, SearchBox};
use docking::scoring::grid::Cache;
use docking::scoring::ScoringFunction;

const V: f64 = 1000.0;

fn read(p: std::path::PathBuf) -> String {
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

#[test]
fn local_only_matches_official_vina() {
    let lig = read(golden_path("1iep", "ligand.pdbqt"));
    let rec = read(golden_path("1iep", "receptor.pdbqt"));
    let mut m = parse_ligand_pdbqt_from_string(&lig, AtomTyping::Xs).unwrap();
    m.attach_receptor(parse_receptor_pdbqt_from_string(&rec, AtomTyping::Xs).unwrap());
    let gbox = SearchBox::from_center_size([15.190, 53.903, 16.917], [20.0, 20.0, 20.0]);
    let sf = ScoringFunction::vina();

    // Build the affinity grids (as Vina does before optimizing).
    let cache = Cache::populate(&m, &sf, &gbox);

    // Optimize from the input pose.
    let mut conf = m.get_initial_conf();
    let steps = default_max_steps(&m);
    let _ = local_optimize_cache(&mut m, &sf, &cache, V, &mut conf, steps);

    // Vina reports the optimized energy via the non_cache score path.
    let got = score(&m, &sf, &gbox);
    let want = ScoreBreakdown::load_golden("1iep", "local_only.out.txt");

    eprintln!(
        "local_only: EFE got {:.3} want {:.3} | inter got {:.3} want {:.3}",
        got.estimated_free_energy,
        want.estimated_free_energy,
        got.intermolecular,
        want.intermolecular
    );

    assert_close(
        "EFE",
        got.estimated_free_energy,
        want.estimated_free_energy,
        tol::LOCAL_ENERGY_KCAL,
    );
    assert_close(
        "intermolecular",
        got.intermolecular,
        want.intermolecular,
        tol::LOCAL_ENERGY_KCAL,
    );

    // Compare the optimized pose to Vina's. Re-emit in file order (the engine
    // stores movable atoms in tree order) so atoms correspond to the golden pose.
    let my_out = m.write_ligand_pdbqt(0, "");
    let my_pose = parse_poses(&my_out);
    let golden = load_golden_poses("1iep", "local_only.pdbqt");
    let rmsd = rmsd(&my_pose[0].coords, &golden[0].coords);
    eprintln!("local_only pose RMSD vs Vina = {rmsd:.4} A");
    assert!(
        rmsd < tol::LOCAL_COORD_ANGSTROM * 5.0,
        "optimized pose RMSD {rmsd:.4} A too large"
    );
}
