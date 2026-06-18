// SPDX-License-Identifier: Apache-2.0
//! Atom types and the atom data model.

pub mod constants;

use crate::math::Vec3;
use constants::*;

/// The four typing schemes: `EL`, `AD`, `XS`, `SY`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomTyping {
    El,
    Ad,
    Xs,
    Sy,
}

/// An atom's type under each scheme. Unassigned fields hold the corresponding
/// `*_TYPE_SIZE` sentinel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomType {
    pub el: usize,
    pub ad: usize,
    pub xs: usize,
    pub sy: usize,
}

impl Default for AtomType {
    fn default() -> Self {
        AtomType {
            el: EL_TYPE_SIZE,
            ad: AD_TYPE_SIZE,
            xs: XS_TYPE_SIZE,
            sy: SY_TYPE_SIZE,
        }
    }
}

impl AtomType {
    /// Returns the type index for the given scheme.
    #[inline]
    pub fn get(&self, typing: AtomTyping) -> usize {
        match typing {
            AtomTyping::El => self.el,
            AtomTyping::Ad => self.ad,
            AtomTyping::Xs => self.xs,
            AtomTyping::Sy => self.sy,
        }
    }

    /// Whether the atom is hydrogen, by its AD type.
    #[inline]
    pub fn is_hydrogen(&self) -> bool {
        ad_is_hydrogen(self.ad)
    }

    /// Whether the atom is a heteroatom (or a metal donor).
    #[inline]
    pub fn is_heteroatom(&self) -> bool {
        ad_is_heteroatom(self.ad) || self.xs == XS_TYPE_MET_D
    }

    /// Whether the atom has an assigned, usable type.
    #[inline]
    pub fn acceptable_type(&self) -> bool {
        self.ad < AD_TYPE_SIZE || self.xs == XS_TYPE_MET_D
    }

    /// Sets the element type from the AD type.
    pub fn assign_el(&mut self) {
        self.el = ad_type_to_el_type(self.ad);
        if self.ad == AD_TYPE_SIZE && self.xs == XS_TYPE_MET_D {
            self.el = EL_TYPE_MET;
        }
    }

    /// Whether two atoms share an element type; does not distinguish metals or
    /// unassigned types.
    #[inline]
    pub fn same_element(&self, a: &AtomType) -> bool {
        self.el == a.el
    }

    /// The atom's covalent radius.
    ///
    /// # Panics
    /// Panics if the atom has no assigned AD type and is not a metal donor —
    /// i.e. it was never typed. Atoms produced by the parser are always typed.
    pub fn covalent_radius(&self) -> f64 {
        if self.ad < AD_TYPE_SIZE {
            ATOM_KIND_DATA[self.ad].covalent_radius
        } else if self.xs == XS_TYPE_MET_D {
            METAL_COVALENT_RADIUS
        } else {
            unreachable!("covalent_radius: unassigned atom type");
        }
    }

    /// The optimal covalent bond length between this atom and `x`.
    #[inline]
    pub fn optimal_covalent_bond_length(&self, x: &AtomType) -> f64 {
        self.covalent_radius() + x.covalent_radius()
    }
}

/// Reference to an atom in either the receptor grid set or the movable set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomIndex {
    pub i: usize,
    pub in_grid: bool,
}

impl AtomIndex {
    #[inline]
    pub fn new(i: usize, in_grid: bool) -> Self {
        AtomIndex { i, in_grid }
    }
}

/// A covalent bond record.
#[derive(Debug, Clone, Copy)]
pub struct Bond {
    pub connected_atom_index: AtomIndex,
    pub length: f64,
    pub rotatable: bool,
}

/// An atom: its types, partial charge, coordinates, and bonds.
///
/// For movable atoms `coords` holds the **frame-relative** coordinates;
/// absolute coordinates live in `Model::coords`.
#[derive(Debug, Clone, Default)]
pub struct Atom {
    pub ty: AtomType,
    pub charge: f64,
    pub coords: Vec3,
    pub bonds: Vec<Bond>,
}

impl Atom {
    /// AD type index (convenience).
    #[inline]
    pub fn ad(&self) -> usize {
        self.ty.ad
    }
    /// Element type index (convenience).
    #[inline]
    pub fn el(&self) -> usize {
        self.ty.el
    }
    /// XS type index (convenience).
    #[inline]
    pub fn xs(&self) -> usize {
        self.ty.xs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_type_is_unassigned() {
        let t = AtomType::default();
        assert_eq!(t.ad, AD_TYPE_SIZE);
        assert_eq!(t.xs, XS_TYPE_SIZE);
        assert!(!t.acceptable_type());
    }

    #[test]
    fn assign_el_from_ad() {
        let mut t = AtomType {
            ad: AD_TYPE_OA,
            ..Default::default()
        };
        t.assign_el();
        assert_eq!(t.el, EL_TYPE_O);
    }

    #[test]
    fn covalent_bond_length() {
        let c = AtomType {
            ad: AD_TYPE_C,
            ..Default::default()
        };
        let n = AtomType {
            ad: AD_TYPE_N,
            ..Default::default()
        };
        assert_eq!(c.optimal_covalent_bond_length(&n), 0.77 + 0.75);
    }
}
