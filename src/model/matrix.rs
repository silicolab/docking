// SPDX-License-Identifier: Apache-2.0
//! Packed triangular matrices (only the parts the engine uses).

/// Packed index for `i <= j`: upper-triangular packing **including** the
/// diagonal. `i + j*(j+1)/2`.
#[inline]
pub fn triangular_index(i: usize, j: usize) -> usize {
    debug_assert!(i <= j);
    i + j * (j + 1) / 2
}

/// Order-independent packed index.
#[inline]
pub fn triangular_index_permissive(i: usize, j: usize) -> usize {
    if i <= j {
        triangular_index(i, j)
    } else {
        triangular_index(j, i)
    }
}

/// Upper-triangular matrix including the diagonal, `n*(n+1)/2` elements.
/// Used for the precalculated potential tables.
#[derive(Debug, Clone)]
pub struct TriangularMatrix<T> {
    data: Vec<T>,
    dim: usize,
}

impl<T: Clone> TriangularMatrix<T> {
    pub fn new(n: usize, filler: T) -> Self {
        TriangularMatrix {
            data: vec![filler; n * (n + 1) / 2],
            dim: n,
        }
    }

    #[inline]
    pub fn dim(&self) -> usize {
        self.dim
    }

    /// Element `(i, j)` with `i <= j`.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> &T {
        &self.data[triangular_index(i, j)]
    }

    #[inline]
    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut T {
        &mut self.data[triangular_index(i, j)]
    }

    /// Element by precomputed packed index.
    #[inline]
    pub fn at(&self, index: usize) -> &T {
        &self.data[index]
    }

    /// Symmetric access: `(i, j)` for any ordering (`triangular_index_permissive`).
    #[inline]
    pub fn get_sym(&self, i: usize, j: usize) -> &T {
        &self.data[triangular_index_permissive(i, j)]
    }
}

/// Strictly-upper-triangular matrix, excluding the diagonal, `n*(n-1)/2`
/// elements. Index for `i < j` is `i + j*(j-1)/2`. Used for the atom-atom
/// mobility (distance-type) matrix.
#[derive(Debug, Clone)]
pub struct StrictlyTriangularMatrix<T> {
    data: Vec<T>,
    dim: usize,
}

impl<T: Clone> StrictlyTriangularMatrix<T> {
    pub fn new(n: usize, filler: T) -> Self {
        StrictlyTriangularMatrix {
            data: vec![filler; n * (n.saturating_sub(1)) / 2],
            dim: n,
        }
    }

    #[inline]
    pub fn dim(&self) -> usize {
        self.dim
    }

    #[inline]
    fn index(i: usize, j: usize) -> usize {
        debug_assert!(i < j);
        i + j * (j - 1) / 2
    }

    /// Element `(i, j)` with `i < j`.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> &T {
        &self.data[Self::index(i, j)]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, value: T) {
        let idx = Self::index(i, j);
        self.data[idx] = value;
    }
}

/// Relative mobility of an atom pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistanceType {
    Fixed,
    Rotor,
    Variable,
}

/// Atom-atom relative mobility matrix.
pub type DistanceTypeMatrix = StrictlyTriangularMatrix<DistanceType>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangular_indices_are_packed() {
        assert_eq!(triangular_index(0, 0), 0);
        assert_eq!(triangular_index(0, 1), 1);
        assert_eq!(triangular_index(1, 1), 2);
        assert_eq!(triangular_index(0, 2), 3);
        assert_eq!(triangular_index(2, 2), 5);
        assert_eq!(triangular_index_permissive(2, 1), triangular_index(1, 2));
    }

    #[test]
    fn triangular_matrix_roundtrips() {
        let mut m = TriangularMatrix::new(4, 0i32);
        for j in 0..4 {
            for i in 0..=j {
                *m.get_mut(i, j) = (i * 10 + j) as i32;
            }
        }
        for j in 0..4 {
            for i in 0..=j {
                assert_eq!(*m.get(i, j), (i * 10 + j) as i32);
            }
        }
    }

    #[test]
    fn strictly_triangular_mobility() {
        let mut m: DistanceTypeMatrix = StrictlyTriangularMatrix::new(5, DistanceType::Variable);
        assert_eq!(*m.get(0, 4), DistanceType::Variable);
        m.set(1, 3, DistanceType::Rotor);
        assert_eq!(*m.get(1, 3), DistanceType::Rotor);
        assert_eq!(*m.get(0, 1), DistanceType::Variable);
    }
}
