// SPDX-License-Identifier: Apache-2.0
//! Single-point scoring against official Vina `--score_only`.
//!
//! Validates that the Rust engine reproduces every component of Vina's energy
//! breakdown for the 1iep input pose within tolerance.

mod common;
use common::*;

use docking::pdbqt::{parse_ligand_pdbqt_from_string, parse_receptor_pdbqt_from_string};
use docking::scoring::eval::{score, SearchBox};
use docking::scoring::ScoringFunction;

fn read(p: std::path::PathBuf) -> String {
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()))
}

fn build_1iep() -> (docking::model::Model, SearchBox) {
    use docking::atom::AtomTyping;
    let lig_text = read(golden_path("1iep", "ligand.pdbqt"));
    let rec_text = read(golden_path("1iep", "receptor.pdbqt"));
    let mut lig = parse_ligand_pdbqt_from_string(&lig_text, AtomTyping::Xs).expect("ligand");
    let rec = parse_receptor_pdbqt_from_string(&rec_text, AtomTyping::Xs).expect("receptor");
    lig.attach_receptor(rec);
    // Box from tests/golden/1iep/config.txt.
    let gbox = SearchBox::from_center_size([15.190, 53.903, 16.917], [20.0, 20.0, 20.0]);
    (lig, gbox)
}

#[test]
fn num_tors_is_seven() {
    let (m, _) = build_1iep();
    assert!(
        (m.num_tors() - 7.0).abs() < 1e-9,
        "num_tors = {}",
        m.num_tors()
    );
}

#[test]
fn score_only_breakdown_matches_official_vina() {
    let (m, gbox) = build_1iep();
    let sf = ScoringFunction::vina();
    let got = score(&m, &sf, &gbox);

    let want = ScoreBreakdown::load_golden("1iep", "score_only.out.txt");

    assert_close(
        "intermolecular (1)",
        got.intermolecular,
        want.intermolecular,
        tol::SCORE_COMPONENT_KCAL,
    );
    assert_close(
        "internal (2)",
        got.internal,
        want.internal,
        tol::SCORE_COMPONENT_KCAL,
    );
    assert_close(
        "torsional (3)",
        got.torsional,
        want.torsional,
        tol::SCORE_COMPONENT_KCAL,
    );
    assert_close(
        "unbound (4)",
        got.unbound,
        want.unbound,
        tol::SCORE_COMPONENT_KCAL,
    );
    assert_close(
        "estimated free energy",
        got.estimated_free_energy,
        want.estimated_free_energy,
        tol::SCORE_COMPONENT_KCAL,
    );
}
