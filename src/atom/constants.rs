// SPDX-License-Identifier: Apache-2.0
//! AutoDock atom-type constants and typing tables.
//!
//! Type indices matter: they are used as array offsets throughout (e.g. into
//! [`ATOM_KIND_DATA`] and [`XS_VDW_RADII`]).

// ----------------------------------------------------------------------------
// Element (EL) types — "based on SY_TYPE_* but includes H".
// ----------------------------------------------------------------------------
pub const EL_TYPE_H: usize = 0;
pub const EL_TYPE_C: usize = 1;
pub const EL_TYPE_N: usize = 2;
pub const EL_TYPE_O: usize = 3;
pub const EL_TYPE_S: usize = 4;
pub const EL_TYPE_P: usize = 5;
pub const EL_TYPE_F: usize = 6;
pub const EL_TYPE_CL: usize = 7;
pub const EL_TYPE_BR: usize = 8;
pub const EL_TYPE_I: usize = 9;
pub const EL_TYPE_SI: usize = 10;
pub const EL_TYPE_AT: usize = 11;
pub const EL_TYPE_MET: usize = 12;
pub const EL_TYPE_DUMMY: usize = 13;
pub const EL_TYPE_SIZE: usize = 14;

// ----------------------------------------------------------------------------
// AutoDock4 (AD) types.
// ----------------------------------------------------------------------------
pub const AD_TYPE_C: usize = 0;
pub const AD_TYPE_A: usize = 1;
pub const AD_TYPE_N: usize = 2;
pub const AD_TYPE_O: usize = 3;
pub const AD_TYPE_P: usize = 4;
pub const AD_TYPE_S: usize = 5;
pub const AD_TYPE_H: usize = 6; // non-polar hydrogen
pub const AD_TYPE_F: usize = 7;
pub const AD_TYPE_I: usize = 8;
pub const AD_TYPE_NA: usize = 9;
pub const AD_TYPE_OA: usize = 10;
pub const AD_TYPE_SA: usize = 11;
pub const AD_TYPE_HD: usize = 12;
pub const AD_TYPE_MG: usize = 13;
pub const AD_TYPE_MN: usize = 14;
pub const AD_TYPE_ZN: usize = 15;
pub const AD_TYPE_CA: usize = 16;
pub const AD_TYPE_FE: usize = 17;
pub const AD_TYPE_CL: usize = 18;
pub const AD_TYPE_BR: usize = 19;
pub const AD_TYPE_SI: usize = 20;
pub const AD_TYPE_AT: usize = 21;
pub const AD_TYPE_G0: usize = 22; // closure of cyclic molecules
pub const AD_TYPE_G1: usize = 23;
pub const AD_TYPE_G2: usize = 24;
pub const AD_TYPE_G3: usize = 25;
pub const AD_TYPE_CG0: usize = 26;
pub const AD_TYPE_CG1: usize = 27;
pub const AD_TYPE_CG2: usize = 28;
pub const AD_TYPE_CG3: usize = 29;
pub const AD_TYPE_W: usize = 30; // hydrated ligand
pub const AD_TYPE_SIZE: usize = 31;

// ----------------------------------------------------------------------------
// X-Score (XS) types.
// ----------------------------------------------------------------------------
pub const XS_TYPE_C_H: usize = 0;
pub const XS_TYPE_C_P: usize = 1;
pub const XS_TYPE_N_P: usize = 2;
pub const XS_TYPE_N_D: usize = 3;
pub const XS_TYPE_N_A: usize = 4;
pub const XS_TYPE_N_DA: usize = 5;
pub const XS_TYPE_O_P: usize = 6;
pub const XS_TYPE_O_D: usize = 7;
pub const XS_TYPE_O_A: usize = 8;
pub const XS_TYPE_O_DA: usize = 9;
pub const XS_TYPE_S_P: usize = 10;
pub const XS_TYPE_P_P: usize = 11;
pub const XS_TYPE_F_H: usize = 12;
pub const XS_TYPE_CL_H: usize = 13;
pub const XS_TYPE_BR_H: usize = 14;
pub const XS_TYPE_I_H: usize = 15;
pub const XS_TYPE_SI: usize = 16;
pub const XS_TYPE_AT: usize = 17;
pub const XS_TYPE_MET_D: usize = 18;
pub const XS_TYPE_C_H_CG0: usize = 19; // closure of cyclic molecules
pub const XS_TYPE_C_P_CG0: usize = 20;
pub const XS_TYPE_G0: usize = 21;
pub const XS_TYPE_C_H_CG1: usize = 22;
pub const XS_TYPE_C_P_CG1: usize = 23;
pub const XS_TYPE_G1: usize = 24;
pub const XS_TYPE_C_H_CG2: usize = 25;
pub const XS_TYPE_C_P_CG2: usize = 26;
pub const XS_TYPE_G2: usize = 27;
pub const XS_TYPE_C_H_CG3: usize = 28;
pub const XS_TYPE_C_P_CG3: usize = 29;
pub const XS_TYPE_G3: usize = 30;
pub const XS_TYPE_W: usize = 31; // hydrated ligand
pub const XS_TYPE_SIZE: usize = 32;

/// DrugScore-CSD (SY) types — unused by the Vina scoring path; only the size is
/// needed for `num_atom_types`.
pub const SY_TYPE_SIZE: usize = 18;

/// Per-AD-type physical properties.
#[derive(Debug, Clone, Copy)]
pub struct AtomKind {
    pub name: &'static str,
    pub radius: f64,
    pub depth: f64,
    /// Pair `(i, j)` is an H-bond if `hb_depth[i] * hb_depth[j] < 0`.
    pub hb_depth: f64,
    pub hb_radius: f64,
    pub solvation: f64,
    pub volume: f64,
    pub covalent_radius: f64,
}

/// Per-AD-type physical-property table, indexed by AD type.
/// 31 entries (indices 0..=30), matching [`AD_TYPE_SIZE`].
#[rustfmt::skip]
pub const ATOM_KIND_DATA: [AtomKind; 31] = [
    // name, radius, depth, hb_depth, hb_radius, solvation, volume, covalent_radius
    AtomKind { name: "C",   radius: 2.00000, depth: 0.15000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00143, volume: 33.51030, covalent_radius: 0.77 }, //  0
    AtomKind { name: "A",   radius: 2.00000, depth: 0.15000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00052, volume: 33.51030, covalent_radius: 0.77 }, //  1
    AtomKind { name: "N",   radius: 1.75000, depth: 0.16000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00162, volume: 22.44930, covalent_radius: 0.75 }, //  2
    AtomKind { name: "O",   radius: 1.60000, depth: 0.20000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00251, volume: 17.15730, covalent_radius: 0.73 }, //  3
    AtomKind { name: "P",   radius: 2.10000, depth: 0.20000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00110, volume: 38.79240, covalent_radius: 1.06 }, //  4
    AtomKind { name: "S",   radius: 2.00000, depth: 0.20000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00214, volume: 33.51030, covalent_radius: 1.02 }, //  5
    AtomKind { name: "H",   radius: 1.00000, depth: 0.02000, hb_depth:  0.0, hb_radius: 0.0, solvation:  0.00051, volume:  0.00000, covalent_radius: 0.37 }, //  6
    AtomKind { name: "F",   radius: 1.54500, depth: 0.08000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00110, volume: 15.44800, covalent_radius: 0.71 }, //  7
    AtomKind { name: "I",   radius: 2.36000, depth: 0.55000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00110, volume: 55.05850, covalent_radius: 1.33 }, //  8
    AtomKind { name: "NA",  radius: 1.75000, depth: 0.16000, hb_depth: -5.0, hb_radius: 1.9, solvation: -0.00162, volume: 22.44930, covalent_radius: 0.75 }, //  9
    AtomKind { name: "OA",  radius: 1.60000, depth: 0.20000, hb_depth: -5.0, hb_radius: 1.9, solvation: -0.00251, volume: 17.15730, covalent_radius: 0.73 }, // 10
    AtomKind { name: "SA",  radius: 2.00000, depth: 0.20000, hb_depth: -1.0, hb_radius: 2.5, solvation: -0.00214, volume: 33.51030, covalent_radius: 1.02 }, // 11
    AtomKind { name: "HD",  radius: 1.00000, depth: 0.02000, hb_depth:  1.0, hb_radius: 0.0, solvation:  0.00051, volume:  0.00000, covalent_radius: 0.37 }, // 12
    AtomKind { name: "Mg",  radius: 0.65000, depth: 0.87500, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00110, volume:  1.56000, covalent_radius: 1.30 }, // 13
    AtomKind { name: "Mn",  radius: 0.65000, depth: 0.87500, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00110, volume:  2.14000, covalent_radius: 1.39 }, // 14
    AtomKind { name: "Zn",  radius: 0.74000, depth: 0.55000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00110, volume:  1.70000, covalent_radius: 1.31 }, // 15
    AtomKind { name: "Ca",  radius: 0.99000, depth: 0.55000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00110, volume:  2.77000, covalent_radius: 1.74 }, // 16
    AtomKind { name: "Fe",  radius: 0.65000, depth: 0.01000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00110, volume:  1.84000, covalent_radius: 1.25 }, // 17
    AtomKind { name: "Cl",  radius: 2.04500, depth: 0.27600, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00110, volume: 35.82350, covalent_radius: 0.99 }, // 18
    AtomKind { name: "Br",  radius: 2.16500, depth: 0.38900, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00110, volume: 42.56610, covalent_radius: 1.14 }, // 19
    AtomKind { name: "Si",  radius: 2.30000, depth: 0.20000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00143, volume: 50.96500, covalent_radius: 1.11 }, // 20
    AtomKind { name: "At",  radius: 2.40000, depth: 0.55000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00110, volume: 57.90580, covalent_radius: 1.44 }, // 21
    AtomKind { name: "G0",  radius: 0.00000, depth: 0.00000, hb_depth:  0.0, hb_radius: 0.0, solvation:  0.00000, volume:  0.00000, covalent_radius: 0.77 }, // 22
    AtomKind { name: "G1",  radius: 0.00000, depth: 0.00000, hb_depth:  0.0, hb_radius: 0.0, solvation:  0.00000, volume:  0.00000, covalent_radius: 0.77 }, // 23
    AtomKind { name: "G2",  radius: 0.00000, depth: 0.00000, hb_depth:  0.0, hb_radius: 0.0, solvation:  0.00000, volume:  0.00000, covalent_radius: 0.77 }, // 24
    AtomKind { name: "G3",  radius: 0.00000, depth: 0.00000, hb_depth:  0.0, hb_radius: 0.0, solvation:  0.00000, volume:  0.00000, covalent_radius: 0.77 }, // 25
    AtomKind { name: "CG0", radius: 2.00000, depth: 0.15000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00143, volume: 33.51030, covalent_radius: 0.77 }, // 26
    AtomKind { name: "CG1", radius: 2.00000, depth: 0.15000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00143, volume: 33.51030, covalent_radius: 0.77 }, // 27
    AtomKind { name: "CG2", radius: 2.00000, depth: 0.15000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00143, volume: 33.51030, covalent_radius: 0.77 }, // 28
    AtomKind { name: "CG3", radius: 2.00000, depth: 0.15000, hb_depth:  0.0, hb_radius: 0.0, solvation: -0.00143, volume: 33.51030, covalent_radius: 0.77 }, // 29
    AtomKind { name: "W",   radius: 0.00000, depth: 0.00000, hb_depth:  0.0, hb_radius: 0.0, solvation:  0.00000, volume:  0.00000, covalent_radius: 0.00 }, // 30
];

/// Covalent radius used for metals not in the AD table.
pub const METAL_COVALENT_RADIUS: f64 = 1.75;

/// Name aliases resolved during typing.
pub const ATOM_EQUIVALENCE_DATA: [(&str, &str); 1] = [("Se", "S")];

/// XS van der Waals radii (Vina), indexed by XS type.
#[rustfmt::skip]
pub const XS_VDW_RADII: [f64; XS_TYPE_SIZE] = [
    1.9, // C_H
    1.9, // C_P
    1.8, // N_P
    1.8, // N_D
    1.8, // N_A
    1.8, // N_DA
    1.7, // O_P
    1.7, // O_D
    1.7, // O_A
    1.7, // O_DA
    2.0, // S_P
    2.1, // P_P
    1.5, // F_H
    1.8, // Cl_H
    2.0, // Br_H
    2.2, // I_H
    2.2, // Si
    2.3, // At
    1.2, // Met_D
    1.9, // C_H_CG0
    1.9, // C_P_CG0
    1.9, // C_H_CG1
    1.9, // C_P_CG1
    1.9, // C_H_CG2
    1.9, // C_P_CG2
    1.9, // C_H_CG3
    1.9, // C_P_CG3
    0.0, // G0
    0.0, // G1
    0.0, // G2
    0.0, // G3
    0.0, // W
];

/// XS radii for the Vinardo variant.
#[rustfmt::skip]
pub const XS_VINARDO_VDW_RADII: [f64; XS_TYPE_SIZE] = [
    2.0, 2.0, 1.7, 1.7, 1.7, 1.7, 1.6, 1.6, 1.6, 1.6, 2.0, 2.1, 1.5, 1.8, 2.0, 2.2,
    2.2, 2.3, 1.2, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0,
];

/// Metal element names not present in the AD4 type table.
pub const NON_AD_METAL_NAMES: [&str; 9] = ["Cu", "Fe", "Na", "K", "Hg", "Co", "U", "Cd", "Ni"];

/// XS vdW radius for XS type `t`.
#[inline]
pub fn xs_radius(t: usize) -> f64 {
    XS_VDW_RADII[t]
}

/// Whether the AD type is a hydrogen (polar or non-polar).
#[inline]
pub fn ad_is_hydrogen(ad: usize) -> bool {
    ad == AD_TYPE_H || ad == AD_TYPE_HD
}

/// Whether the AD type is a heteroatom; returns false for `ad >= AD_TYPE_SIZE`.
#[inline]
pub fn ad_is_heteroatom(ad: usize) -> bool {
    ad != AD_TYPE_A && ad != AD_TYPE_C && ad != AD_TYPE_H && ad != AD_TYPE_HD && ad < AD_TYPE_SIZE
}

/// Maps an AD type to its element (EL) type.
pub fn ad_type_to_el_type(t: usize) -> usize {
    match t {
        AD_TYPE_C => EL_TYPE_C,
        AD_TYPE_A => EL_TYPE_C,
        AD_TYPE_N => EL_TYPE_N,
        AD_TYPE_O => EL_TYPE_O,
        AD_TYPE_P => EL_TYPE_P,
        AD_TYPE_S => EL_TYPE_S,
        AD_TYPE_H => EL_TYPE_H,
        AD_TYPE_F => EL_TYPE_F,
        AD_TYPE_I => EL_TYPE_I,
        AD_TYPE_NA => EL_TYPE_N,
        AD_TYPE_OA => EL_TYPE_O,
        AD_TYPE_SA => EL_TYPE_S,
        AD_TYPE_HD => EL_TYPE_H,
        AD_TYPE_MG => EL_TYPE_MET,
        AD_TYPE_MN => EL_TYPE_MET,
        AD_TYPE_ZN => EL_TYPE_MET,
        AD_TYPE_CA => EL_TYPE_MET,
        AD_TYPE_FE => EL_TYPE_MET,
        AD_TYPE_CL => EL_TYPE_CL,
        AD_TYPE_BR => EL_TYPE_BR,
        AD_TYPE_SI => EL_TYPE_SI,
        AD_TYPE_AT => EL_TYPE_AT,
        AD_TYPE_CG0 | AD_TYPE_CG1 | AD_TYPE_CG2 | AD_TYPE_CG3 => EL_TYPE_C,
        AD_TYPE_G0 | AD_TYPE_G1 | AD_TYPE_G2 | AD_TYPE_G3 => EL_TYPE_DUMMY,
        AD_TYPE_W => EL_TYPE_DUMMY,
        AD_TYPE_SIZE => EL_TYPE_SIZE,
        _ => unreachable!("ad_type_to_el_type: invalid AD type {t}"),
    }
}

/// Whether `name` is a metal element not present in the AD4 type table.
#[inline]
pub fn is_non_ad_metal_name(name: &str) -> bool {
    NON_AD_METAL_NAMES.contains(&name)
}

/// Whether the XS type is one of the macrocycle-closure glue dummies G0-G3.
#[inline]
pub fn is_glue_type(xs: usize) -> bool {
    xs == XS_TYPE_G0 || xs == XS_TYPE_G1 || xs == XS_TYPE_G2 || xs == XS_TYPE_G3
}

/// Whether the pair is a glue dummy paired with its matching CG carbon,
/// used by the `linearattraction` term.
pub fn is_glued(t1: usize, t2: usize) -> bool {
    const PAIRS: [(usize, usize); 4] = [
        (XS_TYPE_G0, XS_TYPE_C_H_CG0),
        (XS_TYPE_G1, XS_TYPE_C_H_CG1),
        (XS_TYPE_G2, XS_TYPE_C_H_CG2),
        (XS_TYPE_G3, XS_TYPE_C_H_CG3),
    ];
    const PAIRS_P: [(usize, usize); 4] = [
        (XS_TYPE_G0, XS_TYPE_C_P_CG0),
        (XS_TYPE_G1, XS_TYPE_C_P_CG1),
        (XS_TYPE_G2, XS_TYPE_C_P_CG2),
        (XS_TYPE_G3, XS_TYPE_C_P_CG3),
    ];
    for (g, c) in PAIRS.into_iter().chain(PAIRS_P) {
        if (t1 == g && t2 == c) || (t2 == g && t1 == c) {
            return true;
        }
    }
    false
}

/// Whether the XS type is hydrophobic.
#[inline]
pub fn xs_is_hydrophobic(xs: usize) -> bool {
    xs == XS_TYPE_C_H
        || xs == XS_TYPE_F_H
        || xs == XS_TYPE_CL_H
        || xs == XS_TYPE_BR_H
        || xs == XS_TYPE_I_H
}

/// Whether the XS type is a hydrogen-bond acceptor.
#[inline]
pub fn xs_is_acceptor(xs: usize) -> bool {
    xs == XS_TYPE_N_A || xs == XS_TYPE_N_DA || xs == XS_TYPE_O_A || xs == XS_TYPE_O_DA
}

/// Whether the XS type is a hydrogen-bond donor.
#[inline]
pub fn xs_is_donor(xs: usize) -> bool {
    xs == XS_TYPE_N_D
        || xs == XS_TYPE_N_DA
        || xs == XS_TYPE_O_D
        || xs == XS_TYPE_O_DA
        || xs == XS_TYPE_MET_D
}

/// Whether `t1` is a donor and `t2` is an acceptor.
#[inline]
pub fn xs_donor_acceptor(t1: usize, t2: usize) -> bool {
    xs_is_donor(t1) && xs_is_acceptor(t2)
}

/// Whether a hydrogen bond is possible between the two XS types (either order).
#[inline]
pub fn xs_h_bond_possible(t1: usize, t2: usize) -> bool {
    xs_donor_acceptor(t1, t2) || xs_donor_acceptor(t2, t1)
}

/// Resolves an atom-type name to its AD type, returning [`AD_TYPE_SIZE`] if not
/// found (no error: metals unknown to AD4 are not exceptional). Resolves
/// equivalence aliases (e.g. `Se -> S`).
pub fn string_to_ad_type(name: &str) -> usize {
    if let Some(i) = ATOM_KIND_DATA.iter().position(|k| k.name == name) {
        return i;
    }
    for (alias, to) in ATOM_EQUIVALENCE_DATA {
        if alias == name {
            return string_to_ad_type(to);
        }
    }
    AD_TYPE_SIZE
}

/// Maximum covalent radius over the AD table.
pub fn max_covalent_radius() -> f64 {
    ATOM_KIND_DATA.iter().fold(0.0, |acc, k| {
        if k.covalent_radius > acc {
            k.covalent_radius
        } else {
            acc
        }
    })
}

/// Number of atom types for each of the four typing schemes.
pub fn num_atom_types(typing: super::AtomTyping) -> usize {
    match typing {
        super::AtomTyping::El => EL_TYPE_SIZE,
        super::AtomTyping::Ad => AD_TYPE_SIZE,
        super::AtomTyping::Xs => XS_TYPE_SIZE,
        super::AtomTyping::Sy => SY_TYPE_SIZE,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ad_table_has_31_entries_indexed_correctly() {
        assert_eq!(ATOM_KIND_DATA.len(), AD_TYPE_SIZE);
        assert_eq!(ATOM_KIND_DATA[AD_TYPE_C].name, "C");
        assert_eq!(ATOM_KIND_DATA[AD_TYPE_HD].name, "HD");
        assert_eq!(ATOM_KIND_DATA[AD_TYPE_W].name, "W");
    }

    #[test]
    fn string_to_ad_type_lookups() {
        assert_eq!(string_to_ad_type("C"), AD_TYPE_C);
        assert_eq!(string_to_ad_type("A"), AD_TYPE_A);
        assert_eq!(string_to_ad_type("OA"), AD_TYPE_OA);
        assert_eq!(string_to_ad_type("HD"), AD_TYPE_HD);
        assert_eq!(string_to_ad_type("Se"), AD_TYPE_S, "Se aliases to S");
        assert_eq!(string_to_ad_type("Xx"), AD_TYPE_SIZE, "unknown -> size");
        // case-sensitive: lowercase 'c' is not a type
        assert_eq!(string_to_ad_type("c"), AD_TYPE_SIZE);
    }

    #[test]
    fn xs_radii_table_size_matches() {
        assert_eq!(XS_VDW_RADII.len(), XS_TYPE_SIZE);
        assert_eq!(XS_VINARDO_VDW_RADII.len(), XS_TYPE_SIZE);
        assert_eq!(xs_radius(XS_TYPE_C_H), 1.9);
        assert_eq!(xs_radius(XS_TYPE_O_A), 1.7);
    }

    #[test]
    fn predicates() {
        assert!(ad_is_hydrogen(AD_TYPE_H));
        assert!(ad_is_hydrogen(AD_TYPE_HD));
        assert!(!ad_is_hydrogen(AD_TYPE_C));
        assert!(!ad_is_heteroatom(AD_TYPE_C));
        assert!(!ad_is_heteroatom(AD_TYPE_A));
        assert!(ad_is_heteroatom(AD_TYPE_NA));
        assert!(xs_is_hydrophobic(XS_TYPE_C_H));
        assert!(!xs_is_hydrophobic(XS_TYPE_C_P));
        assert!(xs_h_bond_possible(XS_TYPE_N_D, XS_TYPE_O_A));
        assert!(!xs_h_bond_possible(XS_TYPE_C_H, XS_TYPE_C_H));
    }

    #[test]
    fn max_covalent_radius_is_calcium() {
        // Ca has the largest covalent radius (1.74) in the table.
        assert_eq!(max_covalent_radius(), 1.74);
    }

    #[test]
    fn el_mapping() {
        assert_eq!(ad_type_to_el_type(AD_TYPE_A), EL_TYPE_C);
        assert_eq!(ad_type_to_el_type(AD_TYPE_OA), EL_TYPE_O);
        assert_eq!(ad_type_to_el_type(AD_TYPE_ZN), EL_TYPE_MET);
        assert_eq!(ad_type_to_el_type(AD_TYPE_W), EL_TYPE_DUMMY);
    }
}
