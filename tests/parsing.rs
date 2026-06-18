// SPDX-License-Identifier: Apache-2.0
//! PDBQT parsing, the torsion-tree model, and coordinate generation.
//!
//! Fidelity checks here:
//! - structural: atom counts, torsion counts, AD types;
//! - round-trip: parse -> re-emit is byte-identical to the input PDBQT;
//! - coordinate generation: applying the initial conformation reproduces the
//!   parsed coordinates, rigid motions preserve internal distances, and a full
//!   2*pi torsion returns to the start.

mod common;

use docking::atom::constants;
use docking::atom::AtomTyping;
use docking::math::{distance_sqr, Vec3, PI};
use docking::pdbqt::{parse_ligand_pdbqt_from_string, parse_receptor_pdbqt_from_string};

fn read(path: std::path::PathBuf) -> String {
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
}

/// Reconstruct the byte-exact text that a faithful re-emission should produce:
/// lines split on '\n' (CR stripped), trailing empty record dropped, each line
/// '\n'-terminated. This matches the parser's line handling.
fn canonical(text: &str) -> String {
    let mut lines: Vec<&str> = text.split('\n').map(|l| l.trim_end_matches('\r')).collect();
    if lines.last() == Some(&"") {
        lines.pop();
    }
    let mut out = String::new();
    for l in lines {
        out.push_str(l);
        out.push('\n');
    }
    out
}

#[test]
fn parse_1iep_ligand_structure() {
    let text = read(common::golden_path("1iep", "ligand.pdbqt"));
    let m = parse_ligand_pdbqt_from_string(&text, AtomTyping::Xs).expect("parse 1iep ligand");
    assert_eq!(m.num_ligands(), 1);
    assert_eq!(m.num_movable_atoms(), 40, "1iep ligand has 40 atoms");
    assert_eq!(m.num_atoms(), 40);
    assert_eq!(m.coords.len(), 40);
    let lig = &m.ligands[0];
    assert_eq!(lig.degrees_of_freedom, 7, "TORSDOF 7");
    assert_eq!(
        lig.body.count_torsions(),
        7,
        "7 rotatable bonds (BRANCH records)"
    );
    // Ligand spans the whole movable-atom range.
    assert_eq!(lig.begin, 0);
    assert_eq!(lig.end, 40);
}

#[test]
fn ad_type_histogram_matches_input() {
    let text = read(common::golden_path("1iep", "ligand.pdbqt"));
    let m = parse_ligand_pdbqt_from_string(&text, AtomTyping::Xs).unwrap();
    let mut counts = std::collections::BTreeMap::<usize, usize>::new();
    for a in &m.atoms {
        *counts.entry(a.ad()).or_default() += 1;
    }
    // From the 1iep ligand PDBQT: 21 A, 8 C, 3 HD, 3 N, 4 NA, 1 OA.
    assert_eq!(counts[&constants::AD_TYPE_A], 21);
    assert_eq!(counts[&constants::AD_TYPE_C], 8);
    assert_eq!(counts[&constants::AD_TYPE_HD], 3);
    assert_eq!(counts[&constants::AD_TYPE_N], 3);
    assert_eq!(counts[&constants::AD_TYPE_NA], 4);
    assert_eq!(counts[&constants::AD_TYPE_OA], 1);
    assert_eq!(counts.values().sum::<usize>(), 40);
}

#[test]
fn roundtrip_1iep_ligand_is_byte_identical() {
    let text = read(common::golden_path("1iep", "ligand.pdbqt"));
    let m = parse_ligand_pdbqt_from_string(&text, AtomTyping::Xs).unwrap();
    let out = m.write_ligand_pdbqt(0, "");
    assert_eq!(
        out,
        canonical(&text),
        "re-emitted ligand differs from input"
    );
}

#[test]
fn roundtrip_all_fixtures_byte_identical() {
    // (file, atom count, branch/torsion count)
    let cases = [
        ("1iep_ligand.pdbqt", 40, 7),
        ("1s63_ligand_zinc.pdbqt", 30, 6),
        ("1uw6_ligand_hydrated.pdbqt", 15, 2),
        ("5x72_ligand_p59.pdbqt", 25, 2),
        ("BACE_1_ligand_macrocycle.pdbqt", 43, 22),
    ];
    for (file, atoms, torsions) in cases {
        let text = read(common::fixture_path(&format!("pdbqt/{file}")));
        let m = parse_ligand_pdbqt_from_string(&text, AtomTyping::Xs)
            .unwrap_or_else(|e| panic!("parse {file}: {e}"));
        assert_eq!(m.num_movable_atoms(), atoms, "{file} atom count");
        assert_eq!(
            m.ligands[0].body.count_torsions(),
            torsions,
            "{file} torsion count"
        );
        assert_eq!(
            m.write_ligand_pdbqt(0, ""),
            canonical(&text),
            "{file} round-trip"
        );
    }
}

#[test]
fn parse_1iep_receptor() {
    let text = read(common::golden_path("1iep", "receptor.pdbqt"));
    let m = parse_receptor_pdbqt_from_string(&text, AtomTyping::Xs).expect("parse 1iep receptor");
    assert_eq!(m.grid_atoms.len(), 2702, "1iep receptor atom count");
    // Receptor atoms carry absolute coordinates and valid AD types.
    assert!(m
        .grid_atoms
        .iter()
        .all(|a| a.ad() < constants::AD_TYPE_SIZE));
}

#[test]
fn initial_conf_reproduces_parsed_coordinates() {
    let text = read(common::golden_path("1iep", "ligand.pdbqt"));
    let mut m = parse_ligand_pdbqt_from_string(&text, AtomTyping::Xs).unwrap();
    let parsed: Vec<Vec3> = m.coords.clone();

    let initial = m.get_initial_conf();
    m.set(&initial);

    let mut max_dev = 0.0_f64;
    for (a, b) in parsed.iter().zip(m.coords.iter()) {
        max_dev = max_dev.max(distance_sqr(a, b).sqrt());
    }
    // Identity orientation + zero torsions should reproduce the parsed pose to
    // within floating-point round-off (origin + (abs - origin)).
    assert!(
        max_dev < 1e-9,
        "initial conf reproduction max deviation = {max_dev:.3e} A"
    );
}

#[test]
fn rigid_motion_preserves_internal_distances() {
    let text = read(common::golden_path("1iep", "ligand.pdbqt"));
    let mut m = parse_ligand_pdbqt_from_string(&text, AtomTyping::Xs).unwrap();
    let before: Vec<Vec3> = m.coords.clone();

    // A pure rigid-body move: translate the root and rotate ~0.7 rad about a
    // tilted axis, leaving all torsions at zero.
    let mut c = m.get_initial_conf();
    c.ligands[0].rigid.position += Vec3::new(5.0, -3.0, 2.0);
    let axis = Vec3::new(1.0, 2.0, -2.0);
    let axis = (1.0 / axis.norm()) * axis;
    c.ligands[0].rigid.orientation = docking::math::angle_to_quaternion_axis(&axis, 0.7);
    m.set(&c);

    // Internal pairwise distances must be invariant under rigid motion.
    let n = before.len();
    let mut max_diff = 0.0_f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let d0 = distance_sqr(&before[i], &before[j]).sqrt();
            let d1 = distance_sqr(&m.coords[i], &m.coords[j]).sqrt();
            max_diff = max_diff.max((d0 - d1).abs());
        }
    }
    assert!(
        max_diff < 1e-9,
        "rigid motion changed internal distances by {max_diff:.3e} A"
    );
    // And it actually moved the molecule.
    let moved = distance_sqr(&before[0], &m.coords[0]).sqrt();
    assert!(moved > 1.0, "rigid motion should displace atoms");
}

#[test]
fn full_2pi_torsion_returns_to_start() {
    let text = read(common::golden_path("1iep", "ligand.pdbqt"));
    let mut m = parse_ligand_pdbqt_from_string(&text, AtomTyping::Xs).unwrap();
    let before: Vec<Vec3> = m.coords.clone();

    let mut c = m.get_initial_conf();
    // Rotate the first torsion by a full turn: geometry must be unchanged.
    c.ligands[0].torsions[0] = 2.0 * PI;
    m.set(&c);

    let mut max_dev = 0.0_f64;
    for (a, b) in before.iter().zip(m.coords.iter()) {
        max_dev = max_dev.max(distance_sqr(a, b).sqrt());
    }
    assert!(
        max_dev < 1e-9,
        "2*pi torsion changed coordinates by {max_dev:.3e} A"
    );
}

#[test]
fn nonzero_torsion_moves_only_the_distal_subtree() {
    let text = read(common::golden_path("1iep", "ligand.pdbqt"));
    let mut m = parse_ligand_pdbqt_from_string(&text, AtomTyping::Xs).unwrap();
    let before: Vec<Vec3> = m.coords.clone();

    let mut c = m.get_initial_conf();
    c.ligands[0].torsions[0] = 1.0; // rotate one bond
    m.set(&c);

    // Some atoms move (the rotated subtree), some do not (root + other subtrees).
    let moved = (0..before.len())
        .filter(|&i| distance_sqr(&before[i], &m.coords[i]).sqrt() > 1e-6)
        .count();
    assert!(
        moved > 0 && moved < before.len(),
        "expected a partial rotation, moved {moved}/40"
    );
}
