// SPDX-License-Identifier: Apache-2.0
//! Torsion tree and coordinate generation.
//!
//! The ligand is a [`FlexibleBody`]: a rigid-body root with a tree of
//! [`Branch`] children, one torsion per branch. `set_conf` walks the tree,
//! placing each frame from the conformation and writing absolute coordinates
//! for that frame's atoms.
//!
//! Atoms store frame-relative coordinates (set at parse time); `set_coords`
//! maps them to the lab frame via `origin + orientation * local`.

use crate::atom::Atom;
use crate::math::{angle_to_quaternion_axis, Mat3, Quat, Vec3, EPSILON_FL, IDENTITY, ZERO};

use super::conf::{LigandChange, LigandConf, RigidConf};

/// A coordinate frame: origin + orientation (as both quaternion and the cached
/// rotation matrix).
#[derive(Debug, Clone, Copy)]
pub struct Frame {
    pub origin: Vec3,
    orientation_q: Quat,
    orientation_m: Mat3,
}

impl Frame {
    pub fn new(origin: Vec3) -> Self {
        Frame {
            origin,
            orientation_q: IDENTITY,
            orientation_m: IDENTITY.to_r3(),
        }
    }

    /// `local_to_lab(local) = origin + orientation_m * local`.
    #[inline]
    pub fn local_to_lab(&self, local: &Vec3) -> Vec3 {
        self.origin + self.orientation_m.mul_vec(local)
    }

    /// `local_to_lab_direction(local) = orientation_m * local`.
    #[inline]
    pub fn local_to_lab_direction(&self, local: &Vec3) -> Vec3 {
        self.orientation_m.mul_vec(local)
    }

    #[inline]
    pub fn orientation(&self) -> Quat {
        self.orientation_q
    }

    /// `set_orientation(q)` — does **not** normalize.
    #[inline]
    fn set_orientation(&mut self, q: Quat) {
        self.orientation_q = q;
        self.orientation_m = q.to_r3();
    }
}

/// Write absolute coordinates for atoms in `[begin, end)` from their
/// frame-relative coordinates.
fn set_coords(frame: &Frame, begin: usize, end: usize, atoms: &[Atom], coords: &mut [Vec3]) {
    for i in begin..end {
        coords[i] = frame.local_to_lab(&atoms[i].coords);
    }
}

/// Net force and torque (about `origin`) of the atoms in `[begin, end)`.
fn sum_force_and_torque(
    origin: Vec3,
    begin: usize,
    end: usize,
    coords: &[Vec3],
    forces: &[Vec3],
) -> (Vec3, Vec3) {
    let mut force = ZERO;
    let mut torque = ZERO;
    for i in begin..end {
        force += forces[i];
        torque += (coords[i] - origin).cross(&forces[i]);
    }
    (force, torque)
}

/// The ligand root.
#[derive(Debug, Clone)]
pub struct RigidBody {
    pub frame: Frame,
    pub begin: usize,
    pub end: usize,
}

impl RigidBody {
    pub fn new(origin: Vec3, begin: usize, end: usize) -> Self {
        RigidBody {
            frame: Frame::new(origin),
            begin,
            end,
        }
    }

    /// Place origin/orientation from the rigid conf, then write this frame's
    /// atom coordinates.
    pub fn set_conf(&mut self, atoms: &[Atom], coords: &mut [Vec3], c: &RigidConf) {
        self.frame.origin = c.position;
        self.frame.set_orientation(c.orientation);
        set_coords(&self.frame, self.begin, self.end, atoms, coords);
    }
}

/// A non-root branch node. Holds the rotation axis plus the parent-relative
/// axis/origin captured at build time (when every frame's orientation was
/// identity).
#[derive(Debug, Clone)]
pub struct Segment {
    pub frame: Frame,
    pub begin: usize,
    pub end: usize,
    axis: Vec3,
    relative_axis: Vec3,
    relative_origin: Vec3,
}

impl Segment {
    /// `parent` must have identity orientation (true at build time). The axis
    /// points from `axis_root` (parent connecting atom) to `origin` (this
    /// branch's connecting atom).
    pub fn new(origin: Vec3, begin: usize, end: usize, axis_root: Vec3, parent: &Frame) -> Self {
        debug_assert!(parent.orientation() == IDENTITY);
        let diff = origin - axis_root;
        let nrm = diff.norm();
        debug_assert!(nrm >= EPSILON_FL);
        let axis = (1.0 / nrm) * diff;
        Segment {
            frame: Frame::new(origin),
            begin,
            end,
            axis,
            relative_axis: axis,
            relative_origin: origin - parent.origin,
        }
    }

    /// Consume one torsion, re-derive origin/axis from the parent frame,
    /// compose the torsion rotation with the parent orientation,
    /// normalize-approx, then write coordinates.
    fn set_conf(
        &mut self,
        parent: &Frame,
        atoms: &[Atom],
        coords: &mut [Vec3],
        torsions: &[f64],
        idx: &mut usize,
    ) {
        let torsion = torsions[*idx];
        *idx += 1;
        self.frame.origin = parent.local_to_lab(&self.relative_origin);
        self.axis = parent.local_to_lab_direction(&self.relative_axis);
        let mut tmp = angle_to_quaternion_axis(&self.axis, torsion).mul(&parent.orientation());
        tmp.normalize_approx();
        self.frame.set_orientation(tmp);
        set_coords(&self.frame, self.begin, self.end, atoms, coords);
    }
}

/// A torsion-tree branch.
#[derive(Debug, Clone)]
pub struct Branch {
    pub node: Segment,
    pub children: Vec<Branch>,
}

impl Branch {
    pub fn new(node: Segment) -> Self {
        Branch {
            node,
            children: Vec::new(),
        }
    }

    fn set_conf(
        &mut self,
        parent: &Frame,
        atoms: &[Atom],
        coords: &mut [Vec3],
        torsions: &[f64],
        idx: &mut usize,
    ) {
        self.node.set_conf(parent, atoms, coords, torsions, idx);
        let node_frame = self.node.frame;
        for child in &mut self.children {
            child.set_conf(&node_frame, atoms, coords, torsions, idx);
        }
    }

    /// `count_torsions` — this segment contributes one, plus children.
    fn count_torsions(&self, s: &mut usize) {
        *s += 1;
        for child in &self.children {
            child.count_torsions(s);
        }
    }
}

fn branches_set_conf(
    branches: &mut [Branch],
    parent: &Frame,
    atoms: &[Atom],
    coords: &mut [Vec3],
    torsions: &[f64],
    idx: &mut usize,
) {
    for b in branches {
        b.set_conf(parent, atoms, coords, torsions, idx);
    }
}

impl Branch {
    /// Accumulate this segment's force/torque (including children) and write
    /// its torsion derivative.
    fn derivative(
        &self,
        coords: &[Vec3],
        forces: &[Vec3],
        torsions: &mut [f64],
        idx: &mut usize,
    ) -> (Vec3, Vec3) {
        let mut ft = sum_force_and_torque(
            self.node.frame.origin,
            self.node.begin,
            self.node.end,
            coords,
            forces,
        );
        let d_idx = *idx;
        *idx += 1;
        branches_derivative(
            &self.children,
            self.node.frame.origin,
            coords,
            forces,
            &mut ft,
            torsions,
            idx,
        );
        // torsion deriv = torque . axis.
        torsions[d_idx] = ft.1.dot(&self.node.axis);
        ft
    }
}

/// `branches_derivative` — fold each child's force/torque into `out` with the
/// parent lever-arm cross product.
fn branches_derivative(
    branches: &[Branch],
    origin: Vec3,
    coords: &[Vec3],
    forces: &[Vec3],
    out: &mut (Vec3, Vec3),
    torsions: &mut [f64],
    idx: &mut usize,
) {
    for b in branches {
        let ft = b.derivative(coords, forces, torsions, idx);
        out.0 += ft.0;
        let r = b.node.frame.origin - origin;
        out.1 += r.cross(&ft.0) + ft.1;
    }
}

/// The ligand tree.
#[derive(Debug, Clone)]
pub struct FlexibleBody {
    pub node: RigidBody,
    pub children: Vec<Branch>,
}

impl FlexibleBody {
    pub fn new(node: RigidBody) -> Self {
        FlexibleBody {
            node,
            children: Vec::new(),
        }
    }

    /// Place the whole tree from a ligand conformation, writing absolute
    /// coordinates for every frame's atoms.
    pub fn set_conf(&mut self, atoms: &[Atom], coords: &mut [Vec3], c: &LigandConf) {
        self.node.set_conf(atoms, coords, &c.rigid);
        let node_frame = self.node.frame;
        let mut idx = 0usize;
        branches_set_conf(
            &mut self.children,
            &node_frame,
            atoms,
            coords,
            &c.torsions,
            &mut idx,
        );
        debug_assert_eq!(idx, c.torsions.len(), "torsion count mismatch in set_conf");
    }

    /// `count_torsions()` — total rotatable bonds in this tree.
    pub fn count_torsions(&self) -> usize {
        let mut s = 0;
        for child in &self.children {
            child.count_torsions(&mut s);
        }
        s
    }

    /// Convert per-atom forces into the conformation-space gradient.
    pub fn derivative(&self, coords: &[Vec3], forces: &[Vec3], c: &mut LigandChange) {
        let mut ft = sum_force_and_torque(
            self.node.frame.origin,
            self.node.begin,
            self.node.end,
            coords,
            forces,
        );
        let mut idx = 0usize;
        branches_derivative(
            &self.children,
            self.node.frame.origin,
            coords,
            forces,
            &mut ft,
            &mut c.torsions,
            &mut idx,
        );
        debug_assert_eq!(idx, c.torsions.len());
        c.rigid.position = ft.0;
        c.rigid.orientation = ft.1;
    }

    /// Union of all node atom ranges (min begin, max end).
    pub fn atom_range(&self) -> (usize, usize) {
        let mut begin = self.node.begin;
        let mut end = self.node.end;
        fn visit(b: &Branch, begin: &mut usize, end: &mut usize) {
            if *begin > b.node.begin {
                *begin = b.node.begin;
            }
            if *end < b.node.end {
                *end = b.node.end;
            }
            for c in &b.children {
                visit(c, begin, end);
            }
        }
        for child in &self.children {
            visit(child, &mut begin, &mut end);
        }
        (begin, end)
    }
}
