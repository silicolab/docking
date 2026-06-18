// SPDX-License-Identifier: Apache-2.0
//! Conformation types.
//!
//! A [`Conf`] is the full set of degrees of freedom: per-ligand rigid-body
//! position + orientation and a torsion angle per rotatable bond, plus per-flex
//! residue torsions.

use crate::math::{normalized_angle, quaternion_increment, Quat, Vec3, IDENTITY, PI, ZERO};
use crate::random::Rng;

/// The number of torsions per ligand and per flex residue.
#[derive(Debug, Clone, Default)]
pub struct ConfSize {
    pub ligands: Vec<usize>,
    pub flex: Vec<usize>,
}

impl ConfSize {
    /// `num_degrees_of_freedom() = sum(ligands) + sum(flex) + 6 * ligands.len()`.
    pub fn num_degrees_of_freedom(&self) -> usize {
        self.ligands.iter().sum::<usize>()
            + self.flex.iter().sum::<usize>()
            + 6 * self.ligands.len()
    }
}

/// Rigid-body position and orientation.
#[derive(Debug, Clone, Copy)]
pub struct RigidConf {
    pub position: Vec3,
    pub orientation: Quat,
}

impl Default for RigidConf {
    fn default() -> Self {
        RigidConf {
            position: ZERO,
            orientation: IDENTITY,
        }
    }
}

impl RigidConf {
    /// Reset to the origin position and identity orientation.
    pub fn set_to_null(&mut self) {
        self.position = ZERO;
        self.orientation = IDENTITY;
    }

    /// Advance by `factor * change`. The orientation is updated via
    /// `quaternion_increment` and is **not** renormalized beyond the approximate
    /// normalization.
    pub fn increment(&mut self, c: &RigidChange, factor: f64) {
        self.position += factor * c.position;
        let rotation = factor * c.orientation;
        quaternion_increment(&mut self.orientation, &rotation);
    }

    /// Random position in the box, random orientation.
    pub fn randomize(&mut self, corner1: &Vec3, corner2: &Vec3, rng: &mut Rng) {
        self.position = rng.random_in_box(corner1, corner2);
        self.orientation = rng.random_orientation();
    }
}

/// Each torsion uniform in `[-pi, pi]`.
fn torsions_randomize(torsions: &mut [f64], rng: &mut Rng) {
    for t in torsions.iter_mut() {
        *t = rng.random_fl(-PI, PI);
    }
}

/// A ligand's rigid body plus its torsion angles.
#[derive(Debug, Clone, Default)]
pub struct LigandConf {
    pub rigid: RigidConf,
    pub torsions: Vec<f64>,
}

impl LigandConf {
    pub fn set_to_null(&mut self) {
        self.rigid.set_to_null();
        for t in &mut self.torsions {
            *t = 0.0;
        }
    }

    /// Advance the rigid body and torsions by `factor * change`.
    pub fn increment(&mut self, c: &LigandChange, factor: f64) {
        self.rigid.increment(&c.rigid, factor);
        torsions_increment(&mut self.torsions, &c.torsions, factor);
    }

    /// Randomize the rigid body and torsions.
    pub fn randomize(&mut self, corner1: &Vec3, corner2: &Vec3, rng: &mut Rng) {
        self.rigid.randomize(corner1, corner2, rng);
        torsions_randomize(&mut self.torsions, rng);
    }
}

/// Advance each torsion by `factor * change`; the result is normalized.
fn torsions_increment(torsions: &mut [f64], c: &[f64], factor: f64) {
    for (t, &ci) in torsions.iter_mut().zip(c.iter()) {
        *t += normalized_angle(factor * ci);
        *t = normalized_angle(*t);
    }
}

/// A flexible residue's torsion angles.
#[derive(Debug, Clone, Default)]
pub struct ResidueConf {
    pub torsions: Vec<f64>,
}

impl ResidueConf {
    pub fn set_to_null(&mut self) {
        for t in &mut self.torsions {
            *t = 0.0;
        }
    }
}

/// The full conformation.
#[derive(Debug, Clone, Default)]
pub struct Conf {
    pub ligands: Vec<LigandConf>,
    pub flex: Vec<ResidueConf>,
}

impl Conf {
    /// Build a zeroed conformation with torsions sized to the per-body counts.
    pub fn new(size: &ConfSize) -> Self {
        let ligands = size
            .ligands
            .iter()
            .map(|&n| LigandConf {
                rigid: RigidConf::default(),
                torsions: vec![0.0; n],
            })
            .collect();
        let flex = size
            .flex
            .iter()
            .map(|&n| ResidueConf {
                torsions: vec![0.0; n],
            })
            .collect();
        Conf { ligands, flex }
    }

    /// Reset every ligand and flex residue to its null pose.
    pub fn set_to_null(&mut self) {
        for l in &mut self.ligands {
            l.set_to_null();
        }
        for f in &mut self.flex {
            f.set_to_null();
        }
    }

    /// Advance by `factor * change` (torsions normalized, orientations not).
    pub fn increment(&mut self, c: &Change, factor: f64) {
        for (l, lc) in self.ligands.iter_mut().zip(&c.ligands) {
            l.increment(lc, factor);
        }
        for (f, fc) in self.flex.iter_mut().zip(&c.flex) {
            torsions_increment(&mut f.torsions, &fc.torsions, factor);
        }
    }

    /// Random rigid body + torsions for every ligand and flex residue.
    pub fn randomize(&mut self, corner1: &Vec3, corner2: &Vec3, rng: &mut Rng) {
        for l in &mut self.ligands {
            l.randomize(corner1, corner2, rng);
        }
        for f in &mut self.flex {
            torsions_randomize(&mut f.torsions, rng);
        }
    }
}

// ---------------------------------------------------------------------------
// Change — the gradient / step in conformation space
// ---------------------------------------------------------------------------

/// Gradient of the rigid-body degrees of freedom.
#[derive(Debug, Clone, Copy, Default)]
pub struct RigidChange {
    pub position: Vec3,
    pub orientation: Vec3,
}

/// Gradient of a ligand's rigid body and torsions.
#[derive(Debug, Clone, Default)]
pub struct LigandChange {
    pub rigid: RigidChange,
    pub torsions: Vec<f64>,
}

/// Gradient of a flexible residue's torsions.
#[derive(Debug, Clone, Default)]
pub struct ResidueChange {
    pub torsions: Vec<f64>,
}

/// The full conformation-space vector (gradient or BFGS direction).
#[derive(Debug, Clone, Default)]
pub struct Change {
    pub ligands: Vec<LigandChange>,
    pub flex: Vec<ResidueChange>,
}

impl Change {
    /// Build a zeroed change with torsions sized to the per-body counts.
    pub fn new(size: &ConfSize) -> Self {
        let ligands = size
            .ligands
            .iter()
            .map(|&n| LigandChange {
                rigid: RigidChange::default(),
                torsions: vec![0.0; n],
            })
            .collect();
        let flex = size
            .flex
            .iter()
            .map(|&n| ResidueChange {
                torsions: vec![0.0; n],
            })
            .collect();
        Change { ligands, flex }
    }

    /// Total degrees of freedom (6 + torsions per ligand, plus flex torsions).
    pub fn num_floats(&self) -> usize {
        self.ligands
            .iter()
            .map(|l| 6 + l.torsions.len())
            .sum::<usize>()
            + self.flex.iter().map(|f| f.torsions.len()).sum::<usize>()
    }

    /// Flat read of the `index`-th degree of freedom (per ligand:
    /// position[0..3], orientation[0..3], torsions; then flex torsions).
    pub fn get(&self, mut index: usize) -> f64 {
        for lig in &self.ligands {
            if index < 3 {
                return lig.rigid.position[index];
            }
            index -= 3;
            if index < 3 {
                return lig.rigid.orientation[index];
            }
            index -= 3;
            if index < lig.torsions.len() {
                return lig.torsions[index];
            }
            index -= lig.torsions.len();
        }
        for res in &self.flex {
            if index < res.torsions.len() {
                return res.torsions[index];
            }
            index -= res.torsions.len();
        }
        panic!("change index out of range");
    }

    /// Mutable flat access.
    pub fn get_mut(&mut self, mut index: usize) -> &mut f64 {
        for lig in &mut self.ligands {
            if index < 3 {
                return &mut lig.rigid.position[index];
            }
            index -= 3;
            if index < 3 {
                return &mut lig.rigid.orientation[index];
            }
            index -= 3;
            if index < lig.torsions.len() {
                return &mut lig.torsions[index];
            }
            index -= lig.torsions.len();
        }
        for res in &mut self.flex {
            if index < res.torsions.len() {
                return &mut res.torsions[index];
            }
            index -= res.torsions.len();
        }
        panic!("change index out of range");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dof_count() {
        let cs = ConfSize {
            ligands: vec![7],
            flex: vec![],
        };
        // 7 torsions + 6 rigid dof = 13
        assert_eq!(cs.num_degrees_of_freedom(), 13);
    }

    #[test]
    fn conf_new_sizes_torsions() {
        let cs = ConfSize {
            ligands: vec![7, 2],
            flex: vec![3],
        };
        let c = Conf::new(&cs);
        assert_eq!(c.ligands.len(), 2);
        assert_eq!(c.ligands[0].torsions.len(), 7);
        assert_eq!(c.ligands[1].torsions.len(), 2);
        assert_eq!(c.flex[0].torsions.len(), 3);
        assert!(c.ligands[0].torsions.iter().all(|&t| t == 0.0));
    }
}
