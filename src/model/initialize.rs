// SPDX-License-Identifier: Apache-2.0
//! `model::initialize` — bond detection, XS atom typing, intramolecular pair
//! setup, and the `num_tors` conf-independent input.
//!
//! Notes on simpler choices:
//! - Bond detection uses a brute-force neighbour scan instead of a `beads`
//!   spatial structure; the resulting bond set is identical (the beads
//!   structure is purely a speed optimization), and bond *order* does not affect
//!   any scoring-relevant query (all are set/any/count reductions).
//! - Macrocycle glue/closure special cases (`is_glue_pair`, `is_closure_clash`,
//!   `is_unmatched_closure_dummy`) are treated as inactive when the model has no
//!   glue/CG atoms (true for all non-macrocycle systems); full macrocycle
//!   support is not yet implemented.

use crate::atom::constants::*;
use crate::atom::{Atom, AtomIndex, Bond};
use crate::math::{distance_sqr, sqr, Vec3};

use super::matrix::{DistanceType, DistanceTypeMatrix};
use super::Model;

const BOND_LENGTH_ALLOWANCE_FACTOR: f64 = 1.1;

impl Model {
    /// Assign bonds and XS types, build interacting pairs, and compute
    /// `num_tors`.
    pub fn initialize(&mut self, mobility: &DistanceTypeMatrix) {
        for lig in &mut self.ligands {
            lig.set_range();
        }
        self.assign_bonds(mobility);
        self.assign_types();
        self.initialize_pairs(mobility);
        self.compute_num_tors();
    }

    // -- atom indexing helpers (grid atoms first, then movable atoms) ----------

    fn total_atoms(&self) -> usize {
        self.grid_atoms.len() + self.atoms.len()
    }

    fn sz_to_atom_index(&self, i: usize) -> AtomIndex {
        if i < self.grid_atoms.len() {
            AtomIndex::new(i, true)
        } else {
            AtomIndex::new(i - self.grid_atoms.len(), false)
        }
    }

    fn get_atom(&self, ai: AtomIndex) -> &Atom {
        if ai.in_grid {
            &self.grid_atoms[ai.i]
        } else {
            &self.atoms[ai.i]
        }
    }

    fn atom_coords(&self, ai: AtomIndex) -> Vec3 {
        if ai.in_grid {
            self.grid_atoms[ai.i].coords
        } else {
            self.coords[ai.i]
        }
    }

    fn distance_sqr_between(&self, a: AtomIndex, b: AtomIndex) -> f64 {
        distance_sqr(&self.atom_coords(a), &self.atom_coords(b))
    }

    /// Relative mobility of two atoms. Grid-grid pairs and grid-to-immovable
    /// pairs are FIXED; grid-to-movable is VARIABLE; movable-movable uses the
    /// mobility matrix.
    fn distance_type_between(
        &self,
        mobility: &DistanceTypeMatrix,
        i: AtomIndex,
        j: AtomIndex,
    ) -> DistanceType {
        if i.in_grid && j.in_grid {
            return DistanceType::Fixed;
        }
        if i.in_grid {
            return if j.i < self.num_movable_atoms {
                DistanceType::Variable
            } else {
                DistanceType::Fixed
            };
        }
        if j.in_grid {
            return if i.i < self.num_movable_atoms {
                DistanceType::Variable
            } else {
                DistanceType::Fixed
            };
        }
        let (a, b) = (i.i, j.i);
        if a == b {
            return DistanceType::Fixed;
        }
        if a < b {
            *mobility.get(a, b)
        } else {
            *mobility.get(b, a)
        }
    }

    // -- bond detection --------------------------------------------------------

    fn covalent_radius_or_max(atom: &Atom, max_cov: f64) -> f64 {
        if atom.ad() < AD_TYPE_SIZE {
            ATOM_KIND_DATA[atom.ad()].covalent_radius
        } else {
            max_cov
        }
    }

    /// True if some third (non-variable, closer) atom lies between `a` and `b`,
    /// which would make their bond spurious.
    fn atom_exists_between(
        &self,
        mobility: &DistanceTypeMatrix,
        a: AtomIndex,
        b: AtomIndex,
        relevant: &[usize],
    ) -> bool {
        let r2 = self.distance_sqr_between(a, b);
        for &i in relevant {
            let c = self.sz_to_atom_index(i);
            if a == c || b == c {
                continue;
            }
            let ac = self.distance_type_between(mobility, a, c);
            let bc = self.distance_type_between(mobility, b, c);
            if ac != DistanceType::Variable
                && bc != DistanceType::Variable
                && self.distance_sqr_between(a, c) < r2
                && self.distance_sqr_between(b, c) < r2
            {
                return true;
            }
        }
        false
    }

    /// Assign covalent bonds from relative mobility, distance, and covalent
    /// radii.
    fn assign_bonds(&mut self, mobility: &DistanceTypeMatrix) {
        let n = self.total_atoms();
        let max_cov = max_covalent_radius();
        // Collected as (owner atom, bond) so we can push to both endpoints without
        // aliasing two mutable atom borrows.
        let mut to_add: Vec<(AtomIndex, Bond)> = Vec::new();

        for i in 0..n {
            let ai = self.sz_to_atom_index(i);
            let i_coords = self.atom_coords(ai);
            let i_cov = Self::covalent_radius_or_max(self.get_atom(ai), max_cov);

            // Relevant atoms: those close enough to possibly bond (distance type
            // not VARIABLE).
            let bond_cutoff = BOND_LENGTH_ALLOWANCE_FACTOR * (i_cov + max_cov);
            let bond_cutoff_sqr = sqr(bond_cutoff);
            let mut relevant: Vec<usize> = Vec::new();
            for j in 0..n {
                if i == j {
                    continue;
                }
                let aj = self.sz_to_atom_index(j);
                if self.distance_type_between(mobility, ai, aj) != DistanceType::Variable {
                    let r2 = distance_sqr(&i_coords, &self.atom_coords(aj));
                    if r2 < bond_cutoff_sqr {
                        relevant.push(j);
                    }
                }
            }

            for &j in &relevant {
                if j <= i {
                    continue;
                }
                let aj = self.sz_to_atom_index(j);
                let bond_length = self
                    .get_atom(ai)
                    .ty
                    .optimal_covalent_bond_length(&self.get_atom(aj).ty);
                let dt = self.distance_type_between(mobility, ai, aj);
                let r2 = self.distance_sqr_between(ai, aj);
                if r2 < sqr(BOND_LENGTH_ALLOWANCE_FACTOR * bond_length)
                    && !self.atom_exists_between(mobility, ai, aj, &relevant)
                {
                    let rotatable = dt == DistanceType::Rotor;
                    let length = r2.sqrt();
                    to_add.push((
                        ai,
                        Bond {
                            connected_atom_index: aj,
                            length,
                            rotatable,
                        },
                    ));
                    to_add.push((
                        aj,
                        Bond {
                            connected_atom_index: ai,
                            length,
                            rotatable,
                        },
                    ));
                }
            }
        }

        for (owner, bond) in to_add {
            if owner.in_grid {
                self.grid_atoms[owner.i].bonds.push(bond);
            } else {
                self.atoms[owner.i].bonds.push(bond);
            }
        }
    }

    fn bonded_to_hd(&self, ai: AtomIndex) -> bool {
        self.get_atom(ai)
            .bonds
            .iter()
            .any(|b| self.get_atom(b.connected_atom_index).ad() == AD_TYPE_HD)
    }

    fn bonded_to_heteroatom(&self, ai: AtomIndex) -> bool {
        self.get_atom(ai)
            .bonds
            .iter()
            .any(|b| self.get_atom(b.connected_atom_index).ty.is_heteroatom())
    }

    // -- XS typing -------------------------------------------------------------

    /// Assign EL and XS types from AD types + bonding.
    fn assign_types(&mut self) {
        let n = self.total_atoms();
        for idx in 0..n {
            let ai = self.sz_to_atom_index(idx);
            // assign_el first; donor_NorO reads the element.
            {
                let a = if ai.in_grid {
                    &mut self.grid_atoms[ai.i]
                } else {
                    &mut self.atoms[ai.i]
                };
                a.ty.assign_el();
            }
            let ad = self.get_atom(ai).ad();
            let el = self.get_atom(ai).el();
            let acceptor = ad == AD_TYPE_OA || ad == AD_TYPE_NA; // SA deliberately ignored
            let donor_n_or_o = el == EL_TYPE_MET || self.bonded_to_hd(ai);

            let xs = match el {
                EL_TYPE_H => None,
                EL_TYPE_C => Some(if ad == AD_TYPE_CG0 {
                    if self.bonded_to_heteroatom(ai) {
                        XS_TYPE_C_P_CG0
                    } else {
                        XS_TYPE_C_H_CG0
                    }
                } else if ad == AD_TYPE_CG1 {
                    if self.bonded_to_heteroatom(ai) {
                        XS_TYPE_C_P_CG1
                    } else {
                        XS_TYPE_C_H_CG1
                    }
                } else if ad == AD_TYPE_CG2 {
                    if self.bonded_to_heteroatom(ai) {
                        XS_TYPE_C_P_CG2
                    } else {
                        XS_TYPE_C_H_CG2
                    }
                } else if ad == AD_TYPE_CG3 {
                    if self.bonded_to_heteroatom(ai) {
                        XS_TYPE_C_P_CG3
                    } else {
                        XS_TYPE_C_H_CG3
                    }
                } else if self.bonded_to_heteroatom(ai) {
                    XS_TYPE_C_P
                } else {
                    XS_TYPE_C_H
                }),
                EL_TYPE_N => Some(match (acceptor, donor_n_or_o) {
                    (true, true) => XS_TYPE_N_DA,
                    (true, false) => XS_TYPE_N_A,
                    (false, true) => XS_TYPE_N_D,
                    (false, false) => XS_TYPE_N_P,
                }),
                EL_TYPE_O => Some(match (acceptor, donor_n_or_o) {
                    (true, true) => XS_TYPE_O_DA,
                    (true, false) => XS_TYPE_O_A,
                    (false, true) => XS_TYPE_O_D,
                    (false, false) => XS_TYPE_O_P,
                }),
                EL_TYPE_S => Some(XS_TYPE_S_P),
                EL_TYPE_P => Some(XS_TYPE_P_P),
                EL_TYPE_F => Some(XS_TYPE_F_H),
                EL_TYPE_CL => Some(XS_TYPE_CL_H),
                EL_TYPE_BR => Some(XS_TYPE_BR_H),
                EL_TYPE_I => Some(XS_TYPE_I_H),
                EL_TYPE_SI => Some(XS_TYPE_SI),
                EL_TYPE_AT => Some(XS_TYPE_AT),
                EL_TYPE_MET => Some(XS_TYPE_MET_D),
                EL_TYPE_DUMMY => Some(match ad {
                    AD_TYPE_G0 => XS_TYPE_G0,
                    AD_TYPE_G1 => XS_TYPE_G1,
                    AD_TYPE_G2 => XS_TYPE_G2,
                    AD_TYPE_G3 => XS_TYPE_G3,
                    AD_TYPE_W => XS_TYPE_SIZE, // no W in XS types
                    _ => unreachable!("dummy element with unexpected AD type {ad}"),
                }),
                _ => None,
            };

            if let Some(xs) = xs {
                let a = if ai.in_grid {
                    &mut self.grid_atoms[ai.i]
                } else {
                    &mut self.atoms[ai.i]
                };
                a.ty.xs = xs;
            }
        }
    }

    // -- interacting pairs -----------------------------------------------------

    /// Movable atoms within `n` bonds of `a` (inclusive), as a set.
    fn bonded_to_set(&self, a: usize, n: usize) -> Vec<usize> {
        let mut out = Vec::new();
        self.bonded_to_recur(a, n, &mut out);
        out
    }

    fn bonded_to_recur(&self, a: usize, n: usize, out: &mut Vec<usize>) {
        if !out.contains(&a) {
            out.push(a);
            if n > 0 {
                for b in &self.atoms[a].bonds {
                    if !b.connected_atom_index.in_grid {
                        self.bonded_to_recur(b.connected_atom_index.i, n - 1, out);
                    }
                }
            }
        }
    }

    fn find_ligand(&self, a: usize) -> Option<usize> {
        self.ligands.iter().position(|l| a >= l.begin && a < l.end)
    }

    fn is_atom_in_ligand(&self, a: usize) -> bool {
        self.ligands.iter().any(|l| a >= l.begin && a < l.end)
    }

    /// Build intramolecular interacting pairs (ligand-internal 1-4+ pairs that
    /// are mobile relative to each other).
    fn initialize_pairs(&mut self, mobility: &DistanceTypeMatrix) {
        let typing = self.atom_typing_used;
        let n_types = num_atom_types(typing);
        let na = self.atoms.len();

        // Accumulate then assign, to avoid borrowing ligands while reading atoms.
        let mut lig_pairs: Vec<Vec<(usize, usize)>> = vec![Vec::new(); self.ligands.len()];
        let mut other_pairs: Vec<(usize, usize)> = Vec::new();

        for i in 0..na {
            let i_lig = self.find_ligand(i);
            let bonded = self.bonded_to_set(i, 3); // 1-2, 1-3, 1-4 exclusion
            for j in (i + 1)..na {
                if *mobility.get(i, j) != DistanceType::Variable || bonded.contains(&j) {
                    continue;
                }
                // Macrocycle closure exclusions are inactive without glue atoms.
                let t1 = self.atoms[i].ty.get(typing);
                let t2 = self.atoms[j].ty.get(typing);
                if t1 < n_types && t2 < n_types {
                    // is_glue_pair is false without glue atoms.
                    match (i_lig, self.find_ligand(j)) {
                        (Some(li), Some(lj)) if li == lj => lig_pairs[li].push((i, j)),
                        _ if !self.is_atom_in_ligand(i) && !self.is_atom_in_ligand(j) => {
                            other_pairs.push((i, j))
                        }
                        _ => {}
                    }
                }
            }
        }

        for (lig, pairs) in self.ligands.iter_mut().zip(lig_pairs) {
            lig.pairs = pairs;
        }
        self.other_pairs = other_pairs;
    }

    // -- num_tors (conf-independent input) -------------------------------------

    fn num_bonded_heavy_atoms(&self, ai: AtomIndex) -> usize {
        self.get_atom(ai)
            .bonds
            .iter()
            .filter(|b| !self.get_atom(b.connected_atom_index).ty.is_hydrogen())
            .count()
    }

    /// Rotatable bonds to heavy ligand atoms that themselves have more than one
    /// heavy neighbour (excludes methyls etc.).
    fn atom_rotors(&self, a: usize) -> usize {
        self.atoms[a]
            .bonds
            .iter()
            .filter(|b| {
                b.rotatable
                    && !self.get_atom(b.connected_atom_index).ty.is_hydrogen()
                    && self.num_bonded_heavy_atoms(b.connected_atom_index) > 1
            })
            .count()
    }

    /// `num_tors = sum over ligand heavy atoms of 0.5 * atom_rotors`.
    fn compute_num_tors(&mut self) {
        let mut num_tors = 0.0;
        for li in 0..self.ligands.len() {
            let (begin, end) = (self.ligands[li].begin, self.ligands[li].end);
            for j in begin..end {
                if self.atoms[j].el() != EL_TYPE_H {
                    num_tors += 0.5 * self.atom_rotors(j) as f64;
                }
            }
        }
        self.set_num_tors(num_tors);
    }
}
