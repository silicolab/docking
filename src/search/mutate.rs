// SPDX-License-Identifier: Apache-2.0
//! Conformation mutation for the Monte-Carlo search.

use crate::math::quaternion_increment;
use crate::math::EPSILON_FL;
use crate::model::conf::Conf;
use crate::model::Model;
use crate::random::Rng;

/// Per ligand: 2 (position + orientation) plus one per torsion; per flex
/// residue: one per torsion.
fn count_mutable_entities(c: &Conf) -> usize {
    c.ligands
        .iter()
        .map(|l| 2 + l.torsions.len())
        .sum::<usize>()
        + c.flex.iter().map(|f| f.torsions.len()).sum::<usize>()
}

/// Randomly perturb exactly one degree of freedom: the rigid position (by
/// `amplitude`), the orientation (by `amplitude / gyration_radius`), or one
/// torsion (re-randomized in `[-pi, pi]`). Does not update the model.
pub fn mutate_conf(c: &mut Conf, m: &Model, amplitude: f64, rng: &mut Rng) {
    let n = count_mutable_entities(c);
    if n == 0 {
        return;
    }
    let mut which = rng.random_sz(0, n - 1);

    for i in 0..c.ligands.len() {
        if which == 0 {
            let delta = amplitude * rng.random_inside_sphere();
            c.ligands[i].rigid.position += delta;
            return;
        }
        which -= 1;
        if which == 0 {
            let gr = m.gyration_radius(i);
            if gr > EPSILON_FL {
                let rotation = (amplitude / gr) * rng.random_inside_sphere();
                quaternion_increment(&mut c.ligands[i].rigid.orientation, &rotation);
            }
            return;
        }
        which -= 1;
        let nt = c.ligands[i].torsions.len();
        if which < nt {
            c.ligands[i].torsions[which] = rng.random_fl(-crate::math::PI, crate::math::PI);
            return;
        }
        which -= nt;
    }
    for f in &mut c.flex {
        let nt = f.torsions.len();
        if which < nt {
            f.torsions[which] = rng.random_fl(-crate::math::PI, crate::math::PI);
            return;
        }
        which -= nt;
    }
}
