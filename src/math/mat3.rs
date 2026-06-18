// SPDX-License-Identifier: Apache-2.0
//! 3x3 matrix.
//!
//! Storage is **column-major**: linear index = `row + 3*col`, so
//! `(i, j)` maps to `data[i + 3*j]`. The element constructor takes
//! arguments in row-major reading order but stores them column-major.

use super::vec3::Vec3;

/// A 3x3 `f64` matrix in column-major storage.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mat3 {
    pub data: [f64; 9],
}

impl Mat3 {
    /// Construct from elements given in row-major reading order
    /// (`xx, xy, xz, yx, yy, yz, zx, zy, zz`), stored column-major.
    #[inline]
    #[allow(clippy::too_many_arguments)]
    pub const fn from_rows(
        xx: f64,
        xy: f64,
        xz: f64,
        yx: f64,
        yy: f64,
        yz: f64,
        zx: f64,
        zy: f64,
        zz: f64,
    ) -> Self {
        let mut data = [0.0; 9];
        data[0] = xx;
        data[3] = xy;
        data[6] = xz;
        data[1] = yx;
        data[4] = yy;
        data[7] = yz;
        data[2] = zx;
        data[5] = zy;
        data[8] = zz;
        Mat3 { data }
    }

    /// `(i, j)` accessor: `data[i + 3*j]` (column-major).
    #[inline]
    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.data[i + 3 * j]
    }

    /// Matrix * vector. Term summation order is preserved for reproducibility.
    #[inline]
    pub fn mul_vec(&self, v: &Vec3) -> Vec3 {
        let d = &self.data;
        Vec3::new(
            d[0] * v[0] + d[3] * v[1] + d[6] * v[2],
            d[1] * v[0] + d[4] * v[1] + d[7] * v[2],
            d[2] * v[0] + d[5] * v[1] + d[8] * v[2],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn column_major_layout() {
        // from_rows args are row-major; storage is column-major (row + 3*col).
        let m = Mat3::from_rows(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
        assert_eq!(m.at(0, 0), 1.0);
        assert_eq!(m.at(0, 1), 2.0);
        assert_eq!(m.at(0, 2), 3.0);
        assert_eq!(m.at(1, 0), 4.0);
        assert_eq!(m.at(2, 2), 9.0);
        // Column-major: data[0..3] is the first column (rows 0,1,2 of col 0).
        assert_eq!(m.data[0], 1.0);
        assert_eq!(m.data[1], 4.0);
        assert_eq!(m.data[2], 7.0);
    }

    #[test]
    fn identity_times_vec() {
        let id = Mat3::from_rows(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        let v = Vec3::new(3.0, -1.0, 2.0);
        assert_eq!(id.mul_vec(&v), v);
    }

    #[test]
    fn mul_vec_matches_formula() {
        let m = Mat3::from_rows(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
        let v = Vec3::new(1.0, 0.0, -1.0);
        // row i dot v
        assert_eq!(m.mul_vec(&v), Vec3::new(1.0 - 3.0, 4.0 - 6.0, 7.0 - 9.0));
    }
}
