// SPDX-License-Identifier: Apache-2.0
//! The molecular model: atoms, coordinates, the ligand torsion tree, and PDBQT
//! re-emission.
//!
//! Bond detection, XS typing, and interaction-pair setup (`model::initialize`)
//! are scoring concerns; this module covers the structure, conformation
//! application, and round-trip output.

pub mod conf;
mod initialize;
pub mod matrix;
pub mod tree;

use crate::atom::{Atom, AtomTyping};
use crate::math::Vec3;

use conf::{Conf, ConfSize};
use tree::FlexibleBody;

/// One re-emission line: the original text plus, for `ATOM`/`HETATM` lines, the
/// index into [`Model::coords`] whose value should overwrite the coordinate
/// columns on output.
pub type ParsedLine = (String, Option<usize>);
/// The full ordered set of lines for a ligand or flex block.
pub type Context = Vec<ParsedLine>;

/// A ligand: its torsion tree, atom range, torsional DOF, and the PDBQT context
/// for output.
#[derive(Debug, Clone)]
pub struct Ligand {
    pub body: FlexibleBody,
    pub begin: usize,
    pub end: usize,
    pub degrees_of_freedom: u32,
    pub cont: Context,
    /// Intramolecular interacting pairs (ligand_i - ligand_i, 1-4+), as movable
    /// atom-index pairs.
    pub pairs: Vec<(usize, usize)>,
}

impl Ligand {
    /// Derive `begin`/`end` from the tree atom ranges.
    pub fn set_range(&mut self) {
        let (b, e) = self.body.atom_range();
        self.begin = b;
        self.end = e;
    }
}

/// The molecular model.
#[derive(Debug, Clone)]
pub struct Model {
    /// Absolute atom coordinates (lab frame).
    pub coords: Vec<Vec3>,
    /// Movable + inflex atoms. For movable atoms, `atom.coords` holds the
    /// **frame-relative** coordinates used by `set_conf`.
    pub atoms: Vec<Atom>,
    /// Rigid receptor atoms.
    pub grid_atoms: Vec<Atom>,
    pub ligands: Vec<Ligand>,
    pub flex_context: Context,
    /// INTRAmolecular flex-flex pairs (empty without flexible residues).
    pub other_pairs: Vec<(usize, usize)>,
    /// INTRAmolecular macrocycle glue pairs (empty without glue atoms).
    pub glue_pairs: Vec<(usize, usize)>,
    /// `num_tors` conf-independent input, computed by [`Model::initialize`].
    num_tors: f64,
    num_movable_atoms: usize,
    atom_typing_used: AtomTyping,
}

impl Model {
    pub(crate) fn new(atom_typing_used: AtomTyping) -> Self {
        Model {
            coords: Vec::new(),
            atoms: Vec::new(),
            grid_atoms: Vec::new(),
            ligands: Vec::new(),
            flex_context: Vec::new(),
            other_pairs: Vec::new(),
            glue_pairs: Vec::new(),
            num_tors: 0.0,
            num_movable_atoms: 0,
            atom_typing_used,
        }
    }

    /// `num_tors` conf-independent input (set by [`Model::initialize`]).
    pub fn num_tors(&self) -> f64 {
        self.num_tors
    }

    pub(crate) fn set_num_tors(&mut self, n: f64) {
        self.num_tors = n;
    }

    pub(crate) fn set_num_movable_atoms(&mut self, n: usize) {
        self.num_movable_atoms = n;
    }

    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }
    pub fn num_movable_atoms(&self) -> usize {
        self.num_movable_atoms
    }
    pub fn num_ligands(&self) -> usize {
        self.ligands.len()
    }
    pub fn atom_typing_used(&self) -> AtomTyping {
        self.atom_typing_used
    }

    /// Attach an already-initialized rigid receptor's atoms to this ligand model
    /// as `grid_atoms`. Ligand and receptor are typed independently (there are no
    /// covalent ligand-receptor bonds), so combining is a move of the grid atoms.
    pub fn attach_receptor(&mut self, receptor: Model) {
        self.grid_atoms = receptor.grid_atoms;
    }

    /// Torsions per ligand and flex residue. `flex` is empty because flexible
    /// receptor residues are not yet modelled.
    pub fn get_size(&self) -> ConfSize {
        ConfSize {
            ligands: self
                .ligands
                .iter()
                .map(|l| l.body.count_torsions())
                .collect(),
            flex: Vec::new(),
        }
    }

    /// Torsions = 0, orientations = identity, ligand positions = current tree
    /// origins. Applying this reproduces the parsed coordinates.
    pub fn get_initial_conf(&self) -> Conf {
        let mut tmp = Conf::new(&self.get_size());
        tmp.set_to_null();
        for (i, lig) in self.ligands.iter().enumerate() {
            tmp.ligands[i].rigid.position = lig.body.node.frame.origin;
        }
        tmp
    }

    /// RMS distance of the ligand's heavy atoms from its rigid-body origin
    /// (uses current coordinates).
    pub fn gyration_radius(&self, ligand_number: usize) -> f64 {
        let lig = &self.ligands[ligand_number];
        let origin = lig.body.node.frame.origin;
        let mut acc = 0.0;
        let mut counter = 0usize;
        for i in lig.begin..lig.end {
            if self.atoms[i].el() != crate::atom::constants::EL_TYPE_H {
                acc += crate::math::distance_sqr(&self.coords[i], &origin);
                counter += 1;
            }
        }
        if counter > 0 {
            (acc / counter as f64).sqrt()
        } else {
            0.0
        }
    }

    /// Current coordinates of the movable heavy atoms (the pose representation
    /// used for clustering/RMSD).
    pub fn heavy_atom_movable_coords(&self) -> Vec<Vec3> {
        (0..self.num_movable_atoms)
            .filter(|&i| self.atoms[i].el() != crate::atom::constants::EL_TYPE_H)
            .map(|i| self.coords[i])
            .collect()
    }

    /// Element types of the movable heavy atoms, in the same order as
    /// [`Self::heavy_atom_movable_coords`]. Used to match atoms by element for
    /// the symmetry-aware lower-bound RMSD.
    pub fn heavy_atom_movable_els(&self) -> Vec<usize> {
        (0..self.num_movable_atoms)
            .map(|i| self.atoms[i].el())
            .filter(|&el| el != crate::atom::constants::EL_TYPE_H)
            .collect()
    }

    /// Regenerate all coordinates from a conformation.
    pub fn set(&mut self, c: &Conf) {
        for (i, lig) in self.ligands.iter_mut().enumerate() {
            lig.body
                .set_conf(&self.atoms, &mut self.coords, &c.ligands[i]);
        }
    }

    /// Render a ligand's PDBQT block with current coordinates substituted into
    /// the coordinate columns. `remark` is prepended verbatim (already
    /// newline-terminated by the caller).
    pub fn write_ligand_pdbqt(&self, ligand_number: usize, remark: &str) -> String {
        let mut out = String::new();
        out.push_str(remark);
        write_context(&self.ligands[ligand_number].cont, &self.coords, &mut out);
        out
    }

    /// Render a full `MODEL ... ENDMDL` block.
    pub fn write_model(&self, model_number: usize, remark: &str) -> String {
        let mut out = String::new();
        out.push_str(&format!("MODEL {model_number}\n"));
        out.push_str(remark);
        for lig in &self.ligands {
            write_context(&lig.cont, &self.coords, &mut out);
        }
        out.push_str("ENDMDL\n");
        out
    }
}

/// Overwrite the three 8-column coordinate fields (PDB columns 31-38, 39-46,
/// 47-54) with `%8.3f` values.
fn coords_to_pdbqt_string(coords: &Vec3, line: &str) -> String {
    let mut bytes: Vec<u8> = line.as_bytes().to_vec();
    write_coord_field(&mut bytes, 30, coords[0]);
    write_coord_field(&mut bytes, 38, coords[1]);
    write_coord_field(&mut bytes, 46, coords[2]);
    String::from_utf8(bytes).expect("pdbqt line is ascii")
}

/// Write an 8-wide `%8.3f` field at 0-based byte offset `start`.
fn write_coord_field(bytes: &mut [u8], start: usize, value: f64) {
    let s = format!("{value:8.3}");
    debug_assert_eq!(s.len(), 8, "coordinate {value} does not fit in 8 columns");
    let field = s.as_bytes();
    for (k, &b) in field.iter().enumerate() {
        bytes[start + k] = b;
    }
}

fn write_context(c: &Context, coords: &[Vec3], out: &mut String) {
    for (str_line, idx) in c {
        match idx {
            Some(i) => {
                out.push_str(&coords_to_pdbqt_string(&coords[*i], str_line));
                out.push('\n');
            }
            None => {
                out.push_str(str_line);
                out.push('\n');
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Vec3;

    #[test]
    fn coord_substitution_overwrites_only_coordinate_columns() {
        // A representative ATOM line from 1iep (79 chars).
        let line =
            "ATOM      1  N   UNL     1      16.600  51.810  14.798  1.00  0.00    -0.322 N ";
        let new = coords_to_pdbqt_string(&Vec3::new(1.0, -2.5, 100.125), line);
        // Columns 31-54 (0-based 30..54) replaced; the rest untouched.
        assert_eq!(&new[0..30], &line[0..30]);
        assert_eq!(&new[30..38], "   1.000");
        assert_eq!(&new[38..46], "  -2.500");
        assert_eq!(&new[46..54], " 100.125");
        assert_eq!(&new[54..], &line[54..]);
    }
}
