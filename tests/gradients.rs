// SPDX-License-Identifier: Apache-2.0
//! Validate the analytic conformation-space gradient against finite
//! differences of the (interpolated) energy the optimizer follows.
//!
//! This is an independent correctness check of the whole derivative path:
//! per-atom forces (`eval_deriv`) and their propagation through the torsion tree
//! into the rigid-body + torsion gradient.

mod common;
use common::*;

use docking::atom::AtomTyping;
use docking::model::conf::{Change, Conf};
use docking::model::Model;
use docking::pdbqt::{parse_ligand_pdbqt_from_string, parse_receptor_pdbqt_from_string};
use docking::scoring::eval::{
    default_max_steps, eval_deriv_with, local_optimize, pair_deriv_exact, SearchBox,
};
use docking::scoring::ScoringFunction;

const V: f64 = 1000.0;

fn build() -> (Model, SearchBox, ScoringFunction) {
    let lig = read_to_string(golden_path("1iep", "ligand.pdbqt"));
    let rec = read_to_string(golden_path("1iep", "receptor.pdbqt"));
    let mut m = parse_ligand_pdbqt_from_string(&lig, AtomTyping::Xs).unwrap();
    m.attach_receptor(parse_receptor_pdbqt_from_string(&rec, AtomTyping::Xs).unwrap());
    let gbox = SearchBox::from_center_size([15.190, 53.903, 16.917], [20.0, 20.0, 20.0]);
    (m, gbox, ScoringFunction::vina())
}

fn read_to_string(p: std::path::PathBuf) -> String {
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

fn energy_at(m: &mut Model, c: &Conf, sf: &ScoringFunction, gbox: &SearchBox) -> f64 {
    m.set(c);
    let mut throwaway = Change::new(&m.get_size());
    // Exact continuous energy (consistent with the exact gradient).
    eval_deriv_with(m, sf, gbox, V, &pair_deriv_exact(sf), &mut throwaway)
}

#[test]
fn analytic_gradient_matches_finite_differences() {
    let (mut m, gbox, sf) = build();

    // A generic, off-minimum conformation so the gradient is non-trivial and the
    // pose is away from interpolation sample boundaries.
    let mut c = m.get_initial_conf();
    c.ligands[0].rigid.position += docking::math::Vec3::new(0.21, -0.13, 0.07);
    let axis = docking::math::Vec3::new(0.3, -0.5, 0.8);
    let axis = (1.0 / axis.norm()) * axis;
    c.ligands[0].rigid.orientation = docking::math::angle_to_quaternion_axis(&axis, 0.17);
    for (k, t) in c.ligands[0].torsions.iter_mut().enumerate() {
        *t = 0.1 * (k as f64 + 1.0);
    }

    // Exact analytic gradient at c (consistent with the exact energy above).
    m.set(&c);
    let mut g = Change::new(&m.get_size());
    let _e = eval_deriv_with(&m, &sf, &gbox, V, &pair_deriv_exact(&sf), &mut g);

    let n = g.num_floats();
    let h = 1e-5;
    let mut max_abs_err = 0.0_f64;
    let mut worst = 0usize;
    for idx in 0..n {
        // Unit direction for this degree of freedom.
        let mut dir = Change::new(&m.get_size());
        *dir.get_mut(idx) = 1.0;

        let mut cp = c.clone();
        cp.increment(&dir, h);
        let ep = energy_at(&mut m, &cp, &sf, &gbox);

        let mut cm = c.clone();
        cm.increment(&dir, -h);
        let em = energy_at(&mut m, &cm, &sf, &gbox);

        let numeric = (ep - em) / (2.0 * h);
        let err = (numeric - g.get(idx)).abs();
        if err > max_abs_err {
            max_abs_err = err;
            worst = idx;
        }
    }

    eprintln!("max |analytic - FD| = {max_abs_err:.3e} at dof {worst} (n={n})");
    assert!(
        max_abs_err < 1e-5,
        "analytic gradient disagrees with finite differences: max err {max_abs_err:.3e} at dof {worst}"
    );
}

#[test]
fn bfgs_reduces_energy_and_converges() {
    let (mut m, gbox, sf) = build();

    // Start from an off-minimum conformation.
    let mut c = m.get_initial_conf();
    c.ligands[0].rigid.position += docking::math::Vec3::new(0.4, -0.3, 0.2);
    for (k, t) in c.ligands[0].torsions.iter_mut().enumerate() {
        *t = 0.15 * (k as f64 + 1.0);
    }

    // Energy before, on the exact potential.
    m.set(&c);
    let mut g0 = Change::new(&m.get_size());
    let e_before = eval_deriv_with(&m, &sf, &gbox, V, &pair_deriv_exact(&sf), &mut g0);

    // Optimize on the exact (consistent) potential.
    let pair = pair_deriv_exact(&sf);
    let steps = default_max_steps(&m);
    let (e_after, _evals) = local_optimize(&mut m, &sf, &gbox, V, &pair, &mut c, steps);

    // The gradient norm at the optimized point.
    let mut g1 = Change::new(&m.get_size());
    let _ = eval_deriv_with(&m, &sf, &gbox, V, &pair_deriv_exact(&sf), &mut g1);
    let gnorm = {
        let n = g1.num_floats();
        let mut s = 0.0;
        for i in 0..n {
            s += g1.get(i) * g1.get(i);
        }
        s.sqrt()
    };

    eprintln!("BFGS: e_before={e_before:.4} e_after={e_after:.4} steps={steps} |g|={gnorm:.4}");
    assert!(
        e_after <= e_before,
        "BFGS increased energy: {e_before} -> {e_after}"
    );
    assert!(
        e_after < e_before - 0.1,
        "BFGS barely moved: {e_before} -> {e_after}"
    );
    // Converged toward a local minimum (gradient much smaller than at the start).
    assert!(
        gnorm < 50.0,
        "gradient still large after optimization: {gnorm}"
    );
}
