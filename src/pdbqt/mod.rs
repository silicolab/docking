// SPDX-License-Identifier: Apache-2.0
//! PDBQT parsing.
//!
//! Parses a prepared ligand (or rigid receptor) PDBQT into a [`Model`]: the
//! torsion tree, frame-relative + absolute coordinates, and the line `context`
//! used for byte-faithful re-emission. Bond detection / XS typing are handled
//! separately (they are scoring concerns).

use crate::atom::constants::{is_non_ad_metal_name, string_to_ad_type, XS_TYPE_MET_D};
use crate::atom::{Atom, AtomType, AtomTyping};
use crate::math::Vec3;
use crate::model::matrix::{DistanceType, DistanceTypeMatrix};
use crate::model::tree::{Branch, FlexibleBody, Frame, RigidBody, Segment};
use crate::model::{Context, Ligand, Model};

/// A PDBQT parse error, with the offending line when available.
#[derive(Debug, Clone)]
pub struct PdbqtError {
    pub message: String,
    pub line: Option<String>,
}

impl PdbqtError {
    fn new(message: impl Into<String>) -> Self {
        PdbqtError {
            message: message.into(),
            line: None,
        }
    }
    fn at(message: impl Into<String>, line: &str) -> Self {
        PdbqtError {
            message: message.into(),
            line: Some(line.to_string()),
        }
    }
}

impl std::fmt::Display for PdbqtError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.line {
            Some(l) => write!(f, "{} (line: {:?})", self.message, l),
            None => write!(f, "{}", self.message),
        }
    }
}

impl std::error::Error for PdbqtError {}

type Result<T> = std::result::Result<T, PdbqtError>;

// ---------------------------------------------------------------------------
// Atom-line parsing (fixed PDB columns, 1-based inclusive ranges)
// ---------------------------------------------------------------------------

/// A parsed atom plus its PDBQT serial number.
#[derive(Debug, Clone)]
struct ParsedAtom {
    atom: Atom,
    number: u32,
}

/// Error if the field's end column `j` (1-based, inclusive) exceeds the line
/// length. A too-short line is rejected with "This line is too short." rather
/// than truncated, so a malformed line errors instead of yielding a
/// silently-wrong value.
fn require_len(line: &str, j: usize) -> Result<()> {
    if j > line.len() {
        Err(PdbqtError::at("This line is too short.", line))
    } else {
        Ok(())
    }
}

/// Trim-extract 1-based inclusive columns `[i, j]` after a length check.
/// PDBQT is ASCII, so byte ranges equal column ranges.
fn checked_col(line: &str, i: usize, j: usize) -> Result<&str> {
    require_len(line, j)?;
    Ok(line[i - 1..j].trim())
}

fn checked_f64(line: &str, i: usize, j: usize, what: &str) -> Result<f64> {
    let s = checked_col(line, i, j)?;
    s.parse::<f64>()
        .map_err(|_| PdbqtError::at(format!("{what} \"{s}\" is not valid."), line))
}

/// True if columns `[i, j]` are all whitespace; errors if the line is too short.
fn substring_is_blank(line: &str, i: usize, j: usize) -> Result<bool> {
    require_len(line, j)?;
    Ok(line[i - 1..j].bytes().all(|b| b.is_ascii_whitespace()))
}

/// Parse a single `ATOM`/`HETATM` line into a [`ParsedAtom`].
fn parse_pdbqt_atom_string(line: &str) -> Result<ParsedAtom> {
    let number_str = checked_col(line, 7, 11)?;
    let number = number_str
        .parse::<u32>()
        .map_err(|_| PdbqtError::at(format!("Atom number \"{number_str}\" is not valid."), line))?;
    let x = checked_f64(line, 31, 38, "Coordinate")?;
    let y = checked_f64(line, 39, 46, "Coordinate")?;
    let z = checked_f64(line, 47, 54, "Coordinate")?;
    // Charge: columns 69-76, only if not blank (the blankness check comes first,
    // and both the blank check and the conversion require the line to be long
    // enough).
    let charge = if substring_is_blank(line, 69, 76)? {
        0.0
    } else {
        checked_f64(line, 69, 76, "Charge")?
    };
    // AutoDock type: column 78 to end of line, trimmed (the read extends to
    // end-of-string rather than a fixed column).
    let name = line.get(77..).unwrap_or("").trim();

    let ad = string_to_ad_type(name);
    let mut ty = AtomType {
        ad,
        ..Default::default()
    };
    if is_non_ad_metal_name(name) {
        ty.xs = XS_TYPE_MET_D;
    }
    if !ty.acceptable_type() {
        return Err(PdbqtError::at(
            format!(
                "Atom type {name} is not a valid AutoDock type (atom types are case-sensitive)."
            ),
            line,
        ));
    }
    let atom = Atom {
        ty,
        charge,
        coords: Vec3::new(x, y, z),
        bonds: Vec::new(),
    };
    Ok(ParsedAtom { atom, number })
}

// Note on integer parsing: a C-style `istringstream >> int` reads a leading
// integer and stops at the first non-digit (so e.g. "2x" parses as 2). We
// intentionally parse the whole whitespace-delimited token and reject trailing
// garbage. This is strictly more conservative and identical for well-formed PDBQT
// (whose BRANCH/TORSDOF operands are pure integers); it differs only on malformed
// input, where rejecting is the safer behavior.
fn parse_one_unsigned(line: &str, start: &str) -> Result<u32> {
    let rest = &line[start.len()..];
    rest.split_whitespace()
        .next()
        .and_then(|t| t.parse::<i64>().ok())
        .filter(|&v| v >= 0)
        .map(|v| v as u32)
        .ok_or_else(|| PdbqtError::at("Syntax error.", line))
}

fn parse_two_unsigneds(line: &str, start: &str) -> Result<(u32, u32)> {
    let rest = &line[start.len()..];
    let mut it = rest.split_whitespace();
    let a = it.next().and_then(|t| t.parse::<i64>().ok());
    let b = it.next().and_then(|t| t.parse::<i64>().ok());
    match (a, b) {
        (Some(a), Some(b)) if a >= 0 && b >= 0 => Ok((a as u32, b as u32)),
        _ => Err(PdbqtError::at("Syntax error.", line)),
    }
}

// ---------------------------------------------------------------------------
// parsing_struct (the in-progress torsion tree)
// ---------------------------------------------------------------------------

/// A node in the parse tree: an atom, the index of its context line, and the
/// child branches that emanate from it.
struct PNode {
    context_index: usize,
    a: ParsedAtom,
    ps: Vec<ParsingStruct>,
    /// The movable-atom index this node was assigned during postprocess (used to
    /// reference branch axis atoms when building the mobility matrix).
    atom_index: Option<usize>,
}

/// An in-progress (sub)tree.
#[derive(Default)]
struct ParsingStruct {
    atoms: Vec<PNode>,
    /// Index in `atoms` of this branch's connecting atom (the `to` atom of its
    /// `BRANCH`), if any.
    immobile_atom: Option<usize>,
}

impl ParsingStruct {
    /// Record an atom and the index of its (already-pushed) line.
    fn add(&mut self, a: ParsedAtom, context_len: usize) {
        self.atoms.push(PNode {
            context_index: context_len - 1,
            a,
            ps: Vec::new(),
            atom_index: None,
        });
    }

    /// True if this subtree holds only the immobile atom and no sub-branches.
    fn essentially_empty(&self) -> bool {
        for (i, nd) in self.atoms.iter().enumerate() {
            if let Some(im) = self.immobile_atom {
                if im != i {
                    return false;
                }
            }
            if !nd.ps.is_empty() {
                return false;
            }
        }
        true
    }
}

// ---------------------------------------------------------------------------
// Line cursor + recursive-descent parse
// ---------------------------------------------------------------------------

struct Cursor {
    lines: Vec<String>,
    pos: usize,
}

impl Cursor {
    fn new(text: &str) -> Self {
        // Split on '\n'; a trailing '\n' does not produce a final empty record.
        // Strip '\r' for CRLF inputs.
        let mut lines: Vec<String> = text
            .split('\n')
            .map(|l| l.trim_end_matches('\r').to_string())
            .collect();
        if lines.last().is_some_and(|l| l.is_empty()) {
            lines.pop();
        }
        Cursor { lines, pos: 0 }
    }

    fn next(&mut self) -> Option<String> {
        if self.pos < self.lines.len() {
            let s = self.lines[self.pos].clone();
            self.pos += 1;
            Some(s)
        } else {
            None
        }
    }
}

fn add_context(c: &mut Context, s: &str) {
    c.push((s.to_string(), None));
}

const MODEL_ERR: &str = "Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. \
                         Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.";

/// Consume up to and including `ROOT`, then parse the root atoms.
fn parse_pdbqt_root(cur: &mut Cursor, p: &mut ParsingStruct, c: &mut Context) -> Result<()> {
    while let Some(str) = cur.next() {
        add_context(c, &str);
        if str.is_empty() || str.starts_with("WARNING") || str.starts_with("REMARK") {
            // ignore
        } else if str.starts_with("ROOT") {
            parse_pdbqt_root_aux(cur, p, c)?;
            break;
        } else if str.starts_with("MODEL") {
            return Err(PdbqtError::new(MODEL_ERR));
        } else {
            return Err(PdbqtError::at(
                "Unknown or inappropriate tag found in flex residue or ligand.",
                &str,
            ));
        }
    }
    Ok(())
}

/// Parse root atoms until `ENDROOT`.
fn parse_pdbqt_root_aux(cur: &mut Cursor, p: &mut ParsingStruct, c: &mut Context) -> Result<()> {
    while let Some(str) = cur.next() {
        add_context(c, &str);
        if str.is_empty() || str.starts_with("WARNING") || str.starts_with("REMARK") {
            // ignore
        } else if str.starts_with("ATOM  ") || str.starts_with("HETATM") {
            let a = parse_pdbqt_atom_string(&str)?;
            p.add(a, c.len());
        } else if str.starts_with("ENDROOT") {
            return Ok(());
        } else if str.starts_with("MODEL") {
            return Err(PdbqtError::new(MODEL_ERR));
        } else {
            return Err(PdbqtError::at(
                "Unknown or inappropriate tag found in flex residue or ligand.",
                &str,
            ));
        }
    }
    Ok(())
}

/// Parse the root, then BRANCH/TORSDOF until end (ligand: `residue=false`).
fn parse_pdbqt_aux(
    cur: &mut Cursor,
    p: &mut ParsingStruct,
    c: &mut Context,
    torsdof: &mut Option<u32>,
) -> Result<()> {
    parse_pdbqt_root(cur, p, c)?;
    while let Some(str) = cur.next() {
        add_context(c, &str);
        if str.is_empty() || str.starts_with("WARNING") || str.starts_with("REMARK") {
            // ignore
        } else if str.starts_with("BRANCH") {
            parse_pdbqt_branch_aux(cur, &str, p, c)?;
        } else if str.starts_with("TORSDOF") {
            if torsdof.is_some() {
                return Err(PdbqtError::new("TORSDOF keyword can be defined only once."));
            }
            *torsdof = Some(parse_one_unsigned(&str, "TORSDOF")?);
        } else if str.starts_with("MODEL") {
            return Err(PdbqtError::new(MODEL_ERR));
        } else {
            return Err(PdbqtError::at(
                "Unknown or inappropriate tag found in flex residue or ligand.",
                &str,
            ));
        }
    }
    Ok(())
}

/// Find the parent atom `first`, then recurse into a child branch.
fn parse_pdbqt_branch_aux(
    cur: &mut Cursor,
    str: &str,
    p: &mut ParsingStruct,
    c: &mut Context,
) -> Result<()> {
    let (first, second) = parse_two_unsigneds(str, "BRANCH")?;
    match p.atoms.iter().position(|n| n.a.number == first) {
        Some(i) => {
            p.atoms[i].ps.push(ParsingStruct::default());
            let child = p.atoms[i].ps.last_mut().unwrap();
            parse_pdbqt_branch(cur, child, c, first, second)
        }
        None => Err(PdbqtError::at(
            format!("Atom number {first} is missing in this branch."),
            str,
        )),
    }
}

/// Parse the atoms of one branch until the matching `ENDBRANCH`.
fn parse_pdbqt_branch(
    cur: &mut Cursor,
    p: &mut ParsingStruct,
    c: &mut Context,
    from: u32,
    to: u32,
) -> Result<()> {
    while let Some(str) = cur.next() {
        add_context(c, &str);
        if str.is_empty() || str.starts_with("WARNING") || str.starts_with("REMARK") {
            // ignore
        } else if str.starts_with("BRANCH") {
            parse_pdbqt_branch_aux(cur, &str, p, c)?;
        } else if str.starts_with("ENDBRANCH") {
            let (first, second) = parse_two_unsigneds(&str, "ENDBRANCH")?;
            if first != from || second != to {
                return Err(PdbqtError::at("Inconsistent branch numbers.", &str));
            }
            if p.immobile_atom.is_none() {
                return Err(PdbqtError::at(
                    format!("Atom {to} has not been found in this branch."),
                    &str,
                ));
            }
            return Ok(());
        } else if str.starts_with("ATOM  ") || str.starts_with("HETATM") {
            let a = parse_pdbqt_atom_string(&str)?;
            if a.number == to {
                p.immobile_atom = Some(p.atoms.len());
            }
            p.add(a, c.len());
        } else if str.starts_with("MODEL") {
            return Err(PdbqtError::new(MODEL_ERR));
        } else {
            return Err(PdbqtError::at(
                "Unknown or inappropriate tag found in flex residue or ligand.",
                &str,
            ));
        }
    }
    Err(PdbqtError::new("Missing ENDBRANCH."))
}

// ---------------------------------------------------------------------------
// postprocess: parsing_struct -> movable atoms + torsion tree
// ---------------------------------------------------------------------------

/// A movable atom carries both its absolute coordinates (in `atom.coords`) and
/// its frame-relative coordinates.
struct MovableAtom {
    atom: Atom,
    relative_coords: Vec3,
}

/// Collects movable atoms and the atom-atom mobility relationships while walking
/// the parse tree.
#[derive(Default)]
struct Postprocess {
    nr_atoms: Vec<MovableAtom>,
    /// Per-frame `[begin, end)` ranges → all within-frame pairs are FIXED.
    fixed_ranges: Vec<(usize, usize)>,
    /// `(axis_atom, begin, end)` → axis atom is FIXED to each atom in the range.
    axis_fixed: Vec<(usize, usize, usize)>,
    /// `(from, to)` rotatable-bond pairs → ROTOR (applied after all FIXED).
    rotors: Vec<(usize, usize)>,
}

impl Postprocess {
    /// Append `p_node`'s atom to the current frame, record its frame-relative
    /// coordinates, link its context line, and remember its index.
    fn insert(&mut self, c: &mut Context, frame_origin: Vec3, p_node: &mut PNode) {
        let idx = self.nr_atoms.len();
        let relative = p_node.a.atom.coords - frame_origin;
        c[p_node.context_index].1 = Some(idx);
        self.nr_atoms.push(MovableAtom {
            atom: p_node.a.atom.clone(),
            relative_coords: relative,
        });
        p_node.atom_index = Some(idx);
    }

    /// Fill the frame's atom range, record mobility, then recurse into child
    /// branches. Returns `(begin, end, children)`.
    fn postprocess_branch(
        &mut self,
        c: &mut Context,
        p: &mut ParsingStruct,
        node_frame: &Frame,
    ) -> (usize, usize, Vec<Branch>) {
        let begin = self.nr_atoms.len();
        let n = p.atoms.len();
        for i in 0..n {
            // Skip this branch's own immobile atom — already inserted by the parent.
            if p.immobile_atom != Some(i) {
                // Borrow p.atoms[i] without aliasing the rest of `p`.
                let p_node = &mut p.atoms[i];
                self.insert(c, node_frame.origin, p_node);
            }
            // Each child's connecting (immobile) atom goes here.
            let m = p.atoms[i].ps.len();
            for j in 0..m {
                if !p.atoms[i].ps[j].atoms.is_empty() {
                    let im = p.atoms[i].ps[j].immobile_atom.expect("immobile atom set");
                    let p_im = &mut p.atoms[i].ps[j].atoms[im];
                    self.insert(c, node_frame.origin, p_im);
                }
            }
        }
        let end = self.nr_atoms.len();
        self.fixed_ranges.push((begin, end));

        let mut children = Vec::new();
        for i in 0..n {
            let axis_begin = p.atoms[i].atom_index; // the parent ("from") atom
            let from_coords = p.atoms[i].a.atom.coords;
            let m = p.atoms[i].ps.len();
            for j in 0..m {
                if !p.atoms[i].ps[j].essentially_empty() {
                    let im = p.atoms[i].ps[j].immobile_atom.expect("immobile atom set");
                    let origin = p.atoms[i].ps[j].atoms[im].a.atom.coords;
                    let axis_end = p.atoms[i].ps[j].atoms[im].atom_index; // the "to" atom
                    let seg = Segment::new(origin, 0, 0, from_coords, node_frame);
                    let child_frame = seg.frame;
                    let (cb, ce, cchildren) =
                        self.postprocess_branch(c, &mut p.atoms[i].ps[j], &child_frame);
                    if let (Some(ab), Some(ae)) = (axis_begin, axis_end) {
                        self.axis_fixed.push((ab, cb, ce));
                        self.axis_fixed.push((ae, cb, ce));
                        self.rotors.push((ab, ae));
                    }
                    let mut seg = seg;
                    seg.begin = cb;
                    seg.end = ce;
                    children.push(Branch {
                        node: seg,
                        children: cchildren,
                    });
                }
            }
        }
        (begin, end, children)
    }

    /// Build the ligand's flexible body from the root.
    fn postprocess_ligand(
        &mut self,
        c: &mut Context,
        p: &mut ParsingStruct,
        torsdof: u32,
    ) -> Ligand {
        let root_origin = p.atoms[0].a.atom.coords;
        let root_frame = Frame::new(root_origin);
        let (begin, end, children) = self.postprocess_branch(c, p, &root_frame);
        let node = RigidBody::new(root_origin, begin, end);
        let mut body = FlexibleBody::new(node);
        body.children = children;
        let mut lig = Ligand {
            body,
            begin: 0,
            end: 0,
            degrees_of_freedom: torsdof,
            cont: Vec::new(),
            pairs: Vec::new(),
        };
        lig.set_range();
        lig
    }

    /// Build the atom-atom mobility matrix from the recorded relationships:
    /// every pair starts `Variable`; within-frame and axis pairs become `Fixed`;
    /// the rotatable-bond pairs become `Rotor` (applied last so they win).
    fn build_mobility(&self) -> DistanceTypeMatrix {
        let n = self.nr_atoms.len();
        let mut mob = DistanceTypeMatrix::new(n, DistanceType::Variable);
        for &(begin, end) in &self.fixed_ranges {
            for i in begin..end {
                for j in (i + 1)..end {
                    mob.set(i, j, DistanceType::Fixed);
                }
            }
        }
        for &(axis, begin, end) in &self.axis_fixed {
            for i in begin..end {
                let (a, b) = if axis < i { (axis, i) } else { (i, axis) };
                mob.set(a, b, DistanceType::Fixed);
            }
        }
        for &(from, to) in &self.rotors {
            let (a, b) = if from < to { (from, to) } else { (to, from) };
            mob.set(a, b, DistanceType::Rotor);
        }
        mob
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parse a prepared ligand PDBQT (from text) into a [`Model`].
pub fn parse_ligand_pdbqt_from_string(text: &str, atom_typing: AtomTyping) -> Result<Model> {
    let mut cur = Cursor::new(text);
    let mut p = ParsingStruct::default();
    let mut c: Context = Vec::new();
    let mut torsdof: Option<u32> = None;

    parse_pdbqt_aux(&mut cur, &mut p, &mut c, &mut torsdof)?;

    if p.atoms.is_empty() {
        return Err(PdbqtError::new("No atoms in this ligand."));
    }
    let torsdof = torsdof.ok_or_else(|| PdbqtError::new("Missing TORSDOF keyword."))?;

    let mut pp = Postprocess::default();
    let mut lig = pp.postprocess_ligand(&mut c, &mut p, torsdof);
    lig.cont = c;
    let mobility = pp.build_mobility();

    // Split absolute (coords) and relative (atom.coords) coordinates.
    let mut m = Model::new(atom_typing);
    m.atoms.reserve(pp.nr_atoms.len());
    m.coords.reserve(pp.nr_atoms.len());
    for ma in &pp.nr_atoms {
        let mut model_atom = ma.atom.clone();
        model_atom.coords = ma.relative_coords; // model atoms store relative coords
        m.atoms.push(model_atom);
        m.coords.push(ma.atom.coords); // absolute coords
    }
    m.set_num_movable_atoms(pp.nr_atoms.len());
    m.ligands.push(lig);

    // model::initialize — bond detection, XS typing, intramolecular pairs.
    m.initialize(&mobility);
    Ok(m)
}

/// Parse a prepared ligand PDBQT file.
pub fn parse_ligand_pdbqt_from_file(
    path: impl AsRef<std::path::Path>,
    atom_typing: AtomTyping,
) -> Result<Model> {
    let text = std::fs::read_to_string(path.as_ref())
        .map_err(|e| PdbqtError::new(format!("cannot read {}: {e}", path.as_ref().display())))?;
    parse_ligand_pdbqt_from_string(&text, atom_typing)
}

/// Read a rigid receptor's atoms (no torsion tree).
fn parse_pdbqt_rigid(text: &str) -> Result<Vec<Atom>> {
    let mut atoms = Vec::new();
    for str in text.split('\n') {
        let str = str.trim_end_matches('\r');
        if str.is_empty()
            || str.starts_with("TER")
            || str.starts_with("END")
            || str.starts_with("WARNING")
            || str.starts_with("REMARK")
        {
            // ignore
        } else if str.starts_with("ATOM  ") || str.starts_with("HETATM") {
            atoms.push(parse_pdbqt_atom_string(str)?.atom);
        } else if str.starts_with("MODEL") {
            return Err(PdbqtError::new(
                "Unexpected multi-MODEL tag found in rigid receptor. \
                 Only one model can be used for the rigid receptor.",
            ));
        } else {
            return Err(PdbqtError::at(
                "Unknown or inappropriate tag found in rigid receptor.",
                str,
            ));
        }
    }
    Ok(atoms)
}

/// Parse a rigid receptor PDBQT (from text) into a [`Model`] whose `grid_atoms`
/// hold the receptor (no flex residues).
pub fn parse_receptor_pdbqt_from_string(text: &str, atom_typing: AtomTyping) -> Result<Model> {
    let atoms = parse_pdbqt_rigid(text)?;
    let mut m = Model::new(atom_typing);
    m.grid_atoms = atoms;
    // Receptor atoms are mutually fixed; initialize with an empty movable-atom
    // mobility matrix so bonds/XS types are assigned to the grid atoms.
    let mobility = DistanceTypeMatrix::new(0, DistanceType::Variable);
    m.initialize(&mobility);
    Ok(m)
}

/// Parse a rigid receptor PDBQT file.
pub fn parse_receptor_pdbqt_from_file(
    path: impl AsRef<std::path::Path>,
    atom_typing: AtomTyping,
) -> Result<Model> {
    let text = std::fs::read_to_string(path.as_ref())
        .map_err(|e| PdbqtError::new(format!("cannot read {}: {e}", path.as_ref().display())))?;
    parse_receptor_pdbqt_from_string(&text, atom_typing)
}

#[cfg(test)]
mod tests {
    use super::*;

    // A well-formed 1iep ligand ATOM line (79 columns).
    const ATOM_LINE: &str =
        "ATOM      1  N   UNL     1      16.600  51.810  14.798  1.00  0.00    -0.322 N ";

    #[test]
    fn parses_a_wellformed_atom_line() {
        let pa = parse_pdbqt_atom_string(ATOM_LINE).expect("valid line");
        assert_eq!(pa.number, 1);
        assert_eq!(pa.atom.coords, Vec3::new(16.600, 51.810, 14.798));
        assert!((pa.atom.charge - (-0.322)).abs() < 1e-12);
        assert_eq!(pa.atom.ty.ad, crate::atom::constants::AD_TYPE_N);
    }

    #[test]
    fn too_short_line_errors_instead_of_truncating() {
        // Cut mid-coordinate: the y/z fields' end columns now exceed the length.
        // This must raise "This line is too short." rather than silently reading
        // a truncated coordinate.
        let truncated = &ATOM_LINE[..35];
        let err = parse_pdbqt_atom_string(truncated).expect_err("should reject short line");
        assert!(err.message.contains("too short"), "got: {}", err.message);
    }

    #[test]
    fn blank_charge_field_defaults_to_zero() {
        // Same columns, but blank the charge field (69-76) while keeping the line
        // long enough; the type stays at column 78.
        let mut chars: Vec<char> = ATOM_LINE.chars().collect();
        for c in chars.iter_mut().take(76).skip(68) {
            *c = ' ';
        }
        let line: String = chars.into_iter().collect();
        let pa = parse_pdbqt_atom_string(&line).expect("valid line with blank charge");
        assert_eq!(pa.atom.charge, 0.0);
    }

    #[test]
    fn invalid_coordinate_reports_which_field() {
        // Put letters in the x-coordinate field (columns 31-38).
        let mut bytes = ATOM_LINE.as_bytes().to_vec();
        for b in bytes.iter_mut().take(38).skip(30) {
            *b = b'z';
        }
        let line = String::from_utf8(bytes).unwrap();
        let err = parse_pdbqt_atom_string(&line).expect_err("bad coordinate");
        assert!(err.message.contains("Coordinate"), "got: {}", err.message);
    }
}
