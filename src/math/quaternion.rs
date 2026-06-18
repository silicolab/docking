// SPDX-License-Identifier: Apache-2.0
//! Quaternion (`boost::math::quaternion<double>`).
//!
//! Component order is `(w, x, y, z)` with the real/scalar part first, matching
//! Boost's `R_component_1..4`. Multiplication is the Hamilton product in the
//! exact term order Boost uses. The rotation-matrix formula and the
//! `normalize_approx` early-out (tolerance `1e-6`) are preserved exactly
//! because they affect coordinate reproducibility.

use super::mat3::Mat3;
use super::vec3::Vec3;
use super::{normalize_angle, EPSILON_FL};

/// A quaternion stored as `(w, x, y, z)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quat {
    pub w: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// The identity quaternion `(1, 0, 0, 0)`.
pub const IDENTITY: Quat = Quat {
    w: 1.0,
    x: 0.0,
    y: 0.0,
    z: 0.0,
};

impl Quat {
    #[inline]
    pub const fn new(w: f64, x: f64, y: f64, z: f64) -> Self {
        Quat { w, x, y, z }
    }

    /// Squared norm `w^2 + x^2 + y^2 + z^2`.
    #[inline]
    pub fn norm_sqr(&self) -> f64 {
        self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Hamilton product `self * rhs`, in Boost's exact term order.
    #[inline]
    pub fn mul(&self, rhs: &Quat) -> Quat {
        let (a, b, c, d) = (self.w, self.x, self.y, self.z);
        let (e, f, g, h) = (rhs.w, rhs.x, rhs.y, rhs.z);
        Quat {
            w: a * e - b * f - c * g - d * h,
            x: a * f + b * e + c * h - d * g,
            y: a * g - b * h + c * e + d * f,
            z: a * h + b * g - c * f + d * e,
        }
    }

    /// Normalizes `q` with tolerance `1e-6`: if the squared norm is within
    /// `tolerance` of 1, leave `q` unchanged (the common path); otherwise
    /// divide by the norm (via reciprocal multiply).
    #[inline]
    pub fn normalize_approx(&mut self) {
        const TOLERANCE: f64 = 1e-6;
        let s = self.norm_sqr();
        if (s - 1.0).abs() < TOLERANCE {
            // most likely scenario — do nothing
        } else {
            let a = s.sqrt();
            debug_assert!(a > EPSILON_FL);
            let inv = 1.0 / a;
            self.w *= inv;
            self.x *= inv;
            self.y *= inv;
            self.z *= inv;
        }
    }

    /// The rotation matrix for a (assumed normalized) quaternion. Formula and
    /// term order are preserved exactly for coordinate reproducibility.
    pub fn to_r3(&self) -> Mat3 {
        let (a, b, c, d) = (self.w, self.x, self.y, self.z);
        let aa = a * a;
        let ab = a * b;
        let ac = a * c;
        let ad = a * d;
        let bb = b * b;
        let bc = b * c;
        let bd = b * d;
        let cc = c * c;
        let cd = c * d;
        let dd = d * d;
        Mat3::from_rows(
            aa + bb - cc - dd,
            2.0 * (-ad + bc),
            2.0 * (ac + bd),
            2.0 * (ad + bc),
            aa - bb + cc - dd,
            2.0 * (-ab + cd),
            2.0 * (-ac + bd),
            2.0 * (ab + cd),
            aa - bb - cc + dd,
        )
    }
}

/// Quaternion for a rotation of `angle` about `axis` (assumed a unit vector).
/// The angle is normalized to `[-pi, pi]` first.
pub fn angle_to_quaternion_axis(axis: &Vec3, mut angle: f64) -> Quat {
    debug_assert!((axis.norm() - 1.0).abs() < 1e-3);
    normalize_angle(&mut angle);
    let c = (angle / 2.0).cos();
    let s = (angle / 2.0).sin();
    Quat::new(c, s * axis[0], s * axis[1], s * axis[2])
}

/// Quaternion for a rotation vector `rotation == angle * axis`.
pub fn angle_to_quaternion(rotation: &Vec3) -> Quat {
    let angle = rotation.norm();
    if angle > EPSILON_FL {
        let axis = (1.0 / angle) * *rotation;
        angle_to_quaternion_axis(&axis, angle)
    } else {
        IDENTITY
    }
}

/// Rotate `q` by the rotation vector, then `normalize_approx` (the
/// normalization step was added in Vina 1.1.2).
pub fn quaternion_increment(q: &mut Quat, rotation: &Vec3) {
    *q = angle_to_quaternion(rotation).mul(q);
    q.normalize_approx();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::PI;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn identity_multiplication() {
        let q = Quat::new(0.5, 0.5, 0.5, 0.5);
        assert_eq!(IDENTITY.mul(&q), q);
        assert_eq!(q.mul(&IDENTITY), q);
    }

    #[test]
    fn identity_rotation_matrix() {
        let m = IDENTITY.to_r3();
        let v = Vec3::new(1.0, 2.0, 3.0);
        assert_eq!(m.mul_vec(&v), v);
    }

    #[test]
    fn rotation_about_z_90_degrees() {
        // 90 deg about z: x -> y, y -> -x.
        let q = angle_to_quaternion_axis(&Vec3::new(0.0, 0.0, 1.0), PI / 2.0);
        let m = q.to_r3();
        let rx = m.mul_vec(&Vec3::new(1.0, 0.0, 0.0));
        let ry = m.mul_vec(&Vec3::new(0.0, 1.0, 0.0));
        assert!(approx_eq(rx[0], 0.0, 1e-12) && approx_eq(rx[1], 1.0, 1e-12));
        assert!(approx_eq(ry[0], -1.0, 1e-12) && approx_eq(ry[1], 0.0, 1e-12));
    }

    #[test]
    fn rotation_matrix_is_orthonormal_and_preserves_length() {
        let axis = Vec3::new(1.0, 2.0, -2.0);
        let axis = (1.0 / axis.norm()) * axis;
        let q = angle_to_quaternion_axis(&axis, 0.7);
        let m = q.to_r3();
        let v = Vec3::new(3.0, -1.0, 2.0);
        let rv = m.mul_vec(&v);
        // A rotation preserves vector length.
        assert!(approx_eq(rv.norm(), v.norm(), 1e-12));
        // Columns are orthonormal.
        let col = |j: usize| Vec3::new(m.at(0, j), m.at(1, j), m.at(2, j));
        assert!(approx_eq(col(0).norm(), 1.0, 1e-12));
        assert!(approx_eq(col(0).dot(&col(1)), 0.0, 1e-12));
        assert!(approx_eq(col(1).dot(&col(2)), 0.0, 1e-12));
    }

    #[test]
    fn angle_axis_composition_adds_angles() {
        // Two 0.3 rad rotations about the same axis == one 0.6 rad rotation.
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let q1 = angle_to_quaternion_axis(&axis, 0.3);
        let combined = q1.mul(&q1);
        let direct = angle_to_quaternion_axis(&axis, 0.6);
        for (a, b) in [
            (combined.w, direct.w),
            (combined.x, direct.x),
            (combined.y, direct.y),
            (combined.z, direct.z),
        ] {
            assert!(approx_eq(a, b, 1e-12), "{a} vs {b}");
        }
    }

    #[test]
    fn normalize_approx_leaves_near_unit_untouched() {
        let mut q = Quat::new(1.0, 1e-4, 0.0, 0.0); // norm_sqr ~ 1 + 1e-8 < 1e-6 off
        let before = q;
        q.normalize_approx();
        assert_eq!(q, before, "near-unit quaternion must be left unchanged");
    }

    #[test]
    fn normalize_approx_rescales_when_far() {
        let mut q = Quat::new(2.0, 0.0, 0.0, 0.0);
        q.normalize_approx();
        assert!(approx_eq(q.norm_sqr(), 1.0, 1e-15));
    }
}
