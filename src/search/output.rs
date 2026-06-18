// SPDX-License-Identifier: Apache-2.0
//! The ranked pose container used by the search.

use crate::math::{distance_sqr, Vec3};
use crate::model::conf::Conf;

/// A candidate pose: its conformation, energy, and heavy-atom coordinates (used
/// for RMSD clustering).
#[derive(Debug, Clone)]
pub struct OutputType {
    pub conf: Conf,
    pub e: f64,
    pub coords: Vec<Vec3>,
}

impl OutputType {
    pub fn new(conf: Conf, e: f64) -> Self {
        OutputType {
            conf,
            e,
            coords: Vec::new(),
        }
    }
}

/// No-superposition RMSD over equal-length, same-order coordinate sets.
pub fn rmsd_upper_bound(a: &[Vec3], b: &[Vec3]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    if a.is_empty() {
        return 0.0;
    }
    let acc: f64 = a.iter().zip(b).map(|(x, y)| distance_sqr(x, y)).sum();
    (acc / a.len() as f64).sqrt()
}

/// Vina's symmetry-aware lower-bound RMSD: each atom is matched to the nearest
/// atom of the same element in the other pose, taking the larger of the two
/// asymmetric directions. `els[k]` is the element type of `a[k]` and `b[k]`
/// (both poses share the same atom order, so a single table applies).
pub fn rmsd_lower_bound(els: &[usize], a: &[Vec3], b: &[Vec3]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), els.len());
    let asymmetric = |from: &[Vec3], to: &[Vec3]| -> f64 {
        if from.is_empty() {
            return 0.0;
        }
        let mut acc = 0.0;
        for (i, p) in from.iter().enumerate() {
            let mut closest = f64::MAX;
            for (j, q) in to.iter().enumerate() {
                if els[i] == els[j] {
                    closest = closest.min(distance_sqr(p, q));
                }
            }
            acc += closest;
        }
        (acc / from.len() as f64).sqrt()
    };
    asymmetric(a, b).max(asymmetric(b, a))
}

/// `(index, rmsd)` of the container pose closest to `coords`.
fn find_closest(coords: &[Vec3], out: &[OutputType]) -> (usize, f64) {
    let mut best = (out.len(), f64::MAX);
    for (i, o) in out.iter().enumerate() {
        let r = rmsd_upper_bound(coords, &o.coords);
        if i == 0 || r < best.1 {
            best = (i, r);
        }
    }
    best
}

/// A sorted (by energy, ascending) container of poses.
#[derive(Debug, Default, Clone)]
pub struct OutputContainer {
    pub poses: Vec<OutputType>,
}

impl OutputContainer {
    pub fn new() -> Self {
        OutputContainer { poses: Vec::new() }
    }

    pub fn len(&self) -> usize {
        self.poses.len()
    }
    pub fn is_empty(&self) -> bool {
        self.poses.is_empty()
    }

    /// Sort by energy ascending.
    pub fn sort(&mut self) {
        self.poses
            .sort_by(|a, b| a.e.partial_cmp(&b.e).unwrap_or(std::cmp::Ordering::Equal));
    }

    /// Insert `t`, merging with a close-by pose (keep the lower energy) or
    /// evicting the worst when full; keep sorted.
    pub fn add(&mut self, t: OutputType, min_rmsd: f64, max_size: usize) {
        let (idx, rmsd) = find_closest(&t.coords, &self.poses);
        if idx < self.poses.len() && rmsd < min_rmsd {
            // A very similar pose already exists; keep the better one.
            if t.e < self.poses[idx].e {
                self.poses[idx] = t;
            }
        } else if self.poses.len() < max_size {
            self.poses.push(t);
        } else if let Some(last) = self.poses.last() {
            if t.e < last.e {
                *self.poses.last_mut().unwrap() = t;
            }
        }
        self.sort();
    }
}
