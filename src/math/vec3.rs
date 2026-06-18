// SPDX-License-Identifier: Apache-2.0
//! 3-vector.
//!
//! Stores three `f64`s as `[0]=x, [1]=y, [2]=z`. Dot and `norm_sqr`
//! accumulate x, then y, then z (left-associative) because floating-point
//! summation order affects reproducibility.

use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Sub, SubAssign};

/// A 3-component `f64` vector.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3 {
    pub data: [f64; 3],
}

/// The zero vector.
pub const ZERO: Vec3 = Vec3 {
    data: [0.0, 0.0, 0.0],
};

impl Vec3 {
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { data: [x, y, z] }
    }

    #[inline]
    pub const fn x(&self) -> f64 {
        self.data[0]
    }
    #[inline]
    pub const fn y(&self) -> f64 {
        self.data[1]
    }
    #[inline]
    pub const fn z(&self) -> f64 {
        self.data[2]
    }

    /// `data[0]^2 + data[1]^2 + data[2]^2`.
    #[inline]
    pub fn norm_sqr(&self) -> f64 {
        self.data[0] * self.data[0] + self.data[1] * self.data[1] + self.data[2] * self.data[2]
    }

    /// `sqrt(norm_sqr)`.
    #[inline]
    pub fn norm(&self) -> f64 {
        self.norm_sqr().sqrt()
    }

    /// Dot product. Order: x, then y, then z.
    #[inline]
    pub fn dot(&self, v: &Vec3) -> f64 {
        self.data[0] * v.data[0] + self.data[1] * v.data[1] + self.data[2] * v.data[2]
    }

    /// Cross product.
    #[inline]
    pub fn cross(&self, b: &Vec3) -> Vec3 {
        Vec3::new(
            self.data[1] * b.data[2] - self.data[2] * b.data[1],
            self.data[2] * b.data[0] - self.data[0] * b.data[2],
            self.data[0] * b.data[1] - self.data[1] * b.data[0],
        )
    }
}

impl Index<usize> for Vec3 {
    type Output = f64;
    #[inline]
    fn index(&self, i: usize) -> &f64 {
        &self.data[i]
    }
}

impl IndexMut<usize> for Vec3 {
    #[inline]
    fn index_mut(&mut self, i: usize) -> &mut f64 {
        &mut self.data[i]
    }
}

impl Add for Vec3 {
    type Output = Vec3;
    #[inline]
    fn add(self, v: Vec3) -> Vec3 {
        Vec3::new(
            self.data[0] + v.data[0],
            self.data[1] + v.data[1],
            self.data[2] + v.data[2],
        )
    }
}

impl Sub for Vec3 {
    type Output = Vec3;
    #[inline]
    fn sub(self, v: Vec3) -> Vec3 {
        Vec3::new(
            self.data[0] - v.data[0],
            self.data[1] - v.data[1],
            self.data[2] - v.data[2],
        )
    }
}

impl AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, v: Vec3) {
        self.data[0] += v.data[0];
        self.data[1] += v.data[1];
        self.data[2] += v.data[2];
    }
}

impl SubAssign for Vec3 {
    #[inline]
    fn sub_assign(&mut self, v: Vec3) {
        self.data[0] -= v.data[0];
        self.data[1] -= v.data[1];
        self.data[2] -= v.data[2];
    }
}

/// Scalar * vector.
impl Mul<Vec3> for f64 {
    type Output = Vec3;
    #[inline]
    fn mul(self, v: Vec3) -> Vec3 {
        Vec3::new(self * v.data[0], self * v.data[1], self * v.data[2])
    }
}

/// Squared Euclidean distance between two points.
#[inline]
pub fn distance_sqr(a: &Vec3, b: &Vec3) -> f64 {
    let dx = a.data[0] - b.data[0];
    let dy = a.data[1] - b.data[1];
    let dz = a.data[2] - b.data[2];
    dx * dx + dy * dy + dz * dz
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dot_and_norm() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, -5.0, 6.0);
        assert_eq!(a.dot(&b), 1.0 * 4.0 + 2.0 * -5.0 + 3.0 * 6.0);
        assert_eq!(a.norm_sqr(), 14.0);
        assert!((a.norm() - 14.0_f64.sqrt()).abs() < 1e-15);
    }

    #[test]
    fn cross_is_right_handed() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        assert_eq!(x.cross(&y), Vec3::new(0.0, 0.0, 1.0));
        assert_eq!(y.cross(&x), Vec3::new(0.0, 0.0, -1.0));
    }

    #[test]
    fn ops_and_distance() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(1.0, 0.0, -1.0);
        assert_eq!(a + b, Vec3::new(2.0, 2.0, 2.0));
        assert_eq!(a - b, Vec3::new(0.0, 2.0, 4.0));
        assert_eq!(2.0 * a, Vec3::new(2.0, 4.0, 6.0));
        assert_eq!(distance_sqr(&a, &b), 0.0 + 4.0 + 16.0);
    }
}
