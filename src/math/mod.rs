// SPDX-License-Identifier: Apache-2.0
//! Numeric primitives: [`Vec3`], [`Mat3`], [`Quat`], and the scalar
//! constants/helpers used throughout the engine.
//!
//! All floating point is `f64`.

mod mat3;
mod quaternion;
mod vec3;

pub use mat3::Mat3;
pub use quaternion::{
    angle_to_quaternion, angle_to_quaternion_axis, quaternion_increment, Quat, IDENTITY,
};
pub use vec3::{distance_sqr, Vec3, ZERO};

/// `pi` — the nearest `f64` to π, bit-identical to [`std::f64::consts::PI`].
pub const PI: f64 = std::f64::consts::PI;

/// Machine epsilon for `f64`.
pub const EPSILON_FL: f64 = f64::EPSILON;

/// `x^2`.
#[inline]
pub fn sqr(x: f64) -> f64 {
    x * x
}

/// `normalize_angle(x)`: add/subtract enough `2*pi` to bring `x` into
/// `[-pi, pi]`, recursing for very large/small inputs.
pub fn normalize_angle(x: &mut f64) {
    if *x > 3.0 * PI {
        // very large
        let n = (*x - PI) / (2.0 * PI); // how many 2*pi's to subtract
        *x -= 2.0 * PI * n.ceil();
        normalize_angle(x);
    } else if *x < -3.0 * PI {
        // very small
        let n = (-*x - PI) / (2.0 * PI); // how many 2*pi's to add
        *x += 2.0 * PI * n.ceil();
        normalize_angle(x);
    } else if *x > PI {
        // in (pi, 3*pi]
        *x -= 2.0 * PI;
    } else if *x < -PI {
        // in [-3*pi, -pi)
        *x += 2.0 * PI;
    }
    debug_assert!(*x >= -PI && *x <= PI);
}

/// `normalized_angle(x)` — value-returning form.
#[inline]
pub fn normalized_angle(mut x: f64) -> f64 {
    normalize_angle(&mut x);
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pi_is_the_ieee754_nearest_double_to_pi() {
        // The nearest f64 to π has bit pattern 0x400921FB54442D18. Confirm ours is identical.
        assert_eq!(PI.to_bits(), 0x4009_21FB_5444_2D18);
    }

    #[test]
    fn normalize_angle_wraps_into_range() {
        for &x in &[0.0, 3.0, -3.0, PI, -PI, 4.0, -4.0, 100.0, -100.0, 7.0 * PI] {
            let n = normalized_angle(x);
            assert!((-PI..=PI).contains(&n), "{x} -> {n} out of range");
            // Differs from x by an integer multiple of 2*pi.
            let k = ((x - n) / (2.0 * PI)).round();
            assert!((x - n - k * 2.0 * PI).abs() < 1e-9, "{x} -> {n}");
        }
    }
}
