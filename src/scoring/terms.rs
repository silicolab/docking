// SPDX-License-Identifier: Apache-2.0
//! Vina scoring potential terms (the Vinardo/AD4 variants are not yet implemented).
//!
//! Each term evaluates on an XS atom-type pair `(t1, t2)` at separation `r`.

use crate::atom::constants::{
    is_glue_type, xs_h_bond_possible, xs_is_hydrophobic, xs_radius, XS_TYPE_SIZE,
};
use crate::math::sqr;

/// `slope_step(x_bad, x_good, x)` — a clamped linear ramp from 0 at `x_bad` to 1
/// at `x_good` (handles either ordering of bad/good).
#[inline]
pub fn slope_step(x_bad: f64, x_good: f64, x: f64) -> f64 {
    if x_bad < x_good {
        if x <= x_bad {
            return 0.0;
        }
        if x >= x_good {
            return 1.0;
        }
    } else {
        if x >= x_bad {
            return 0.0;
        }
        if x <= x_good {
            return 1.0;
        }
    }
    (x - x_bad) / (x_good - x_bad)
}

/// `optimal_distance(t1, t2)` — sum of XS vdW radii, or 0 for glue types.
#[inline]
pub fn optimal_distance(t1: usize, t2: usize) -> f64 {
    if is_glue_type(t1) || is_glue_type(t2) {
        0.0
    } else {
        xs_radius(t1) + xs_radius(t2)
    }
}

/// The Vina scoring terms, in the fixed order used by `SF_VINA`.
#[derive(Debug, Clone, Copy)]
pub enum VinaTerm {
    /// `vina_gaussian(offset, width, cutoff)`.
    Gaussian {
        offset: f64,
        width: f64,
        cutoff: f64,
    },
    /// `vina_repulsion(offset, cutoff)`.
    Repulsion { offset: f64, cutoff: f64 },
    /// `vina_hydrophobic(good, bad, cutoff)`.
    Hydrophobic { good: f64, bad: f64, cutoff: f64 },
    /// `vina_non_dir_h_bond(good, bad, cutoff)`.
    NonDirHBond { good: f64, bad: f64, cutoff: f64 },
    /// `linearattraction(cutoff)` — macrocycle glue.
    LinearAttraction { cutoff: f64 },
}

impl VinaTerm {
    pub fn cutoff(&self) -> f64 {
        match *self {
            VinaTerm::Gaussian { cutoff, .. }
            | VinaTerm::Repulsion { cutoff, .. }
            | VinaTerm::Hydrophobic { cutoff, .. }
            | VinaTerm::NonDirHBond { cutoff, .. }
            | VinaTerm::LinearAttraction { cutoff } => cutoff,
        }
    }

    /// Evaluate the term for XS type pair `(t1, t2)` at separation `r`.
    pub fn eval(&self, t1: usize, t2: usize, r: f64) -> f64 {
        match *self {
            VinaTerm::Gaussian {
                offset,
                width,
                cutoff,
            } => {
                if r >= cutoff {
                    return 0.0;
                }
                gauss(r - (optimal_distance(t1, t2) + offset), width)
            }
            VinaTerm::Repulsion { offset, cutoff } => {
                if r >= cutoff {
                    return 0.0;
                }
                let d = r - (optimal_distance(t1, t2) + offset);
                if d > 0.0 {
                    0.0
                } else {
                    d * d
                }
            }
            VinaTerm::Hydrophobic { good, bad, cutoff } => {
                if r >= cutoff {
                    return 0.0;
                }
                if xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2) {
                    slope_step(bad, good, r - optimal_distance(t1, t2))
                } else {
                    0.0
                }
            }
            VinaTerm::NonDirHBond { good, bad, cutoff } => {
                if r >= cutoff {
                    return 0.0;
                }
                if xs_h_bond_possible(t1, t2) {
                    slope_step(bad, good, r - optimal_distance(t1, t2))
                } else {
                    0.0
                }
            }
            VinaTerm::LinearAttraction { cutoff } => {
                if r >= cutoff {
                    return 0.0;
                }
                if crate::atom::constants::is_glued(t1, t2) {
                    r
                } else {
                    0.0
                }
            }
        }
    }
}

/// `gauss(x) = exp(-(x/width)^2)`.
#[inline]
fn gauss(x: f64, width: f64) -> f64 {
    (-sqr(x / width)).exp()
}

/// Exact derivative of [`slope_step`] with respect to `x`: the ramp slope
/// `1/(x_good - x_bad)` inside the ramp, zero outside.
#[inline]
pub fn slope_step_deriv(x_bad: f64, x_good: f64, x: f64) -> f64 {
    if x_bad < x_good {
        if x <= x_bad || x >= x_good {
            return 0.0;
        }
    } else if x >= x_bad || x <= x_good {
        return 0.0;
    }
    1.0 / (x_good - x_bad)
}

impl VinaTerm {
    /// Exact `(value, d value / d r)` for the term — the analytic derivative of
    /// the continuous potential (used to validate the engine's interpolated
    /// derivative path via finite differences). Discontinuities at the cutoff and
    /// at `slope_step` knees are measure-zero and irrelevant for generic poses.
    pub fn eval_deriv_r(&self, t1: usize, t2: usize, r: f64) -> (f64, f64) {
        if r >= self.cutoff() {
            return (0.0, 0.0);
        }
        match *self {
            VinaTerm::Gaussian { offset, width, .. } => {
                let x = r - (optimal_distance(t1, t2) + offset);
                let g = gauss(x, width);
                (g, g * (-2.0 * x / (width * width)))
            }
            VinaTerm::Repulsion { offset, .. } => {
                let d = r - (optimal_distance(t1, t2) + offset);
                if d > 0.0 {
                    (0.0, 0.0)
                } else {
                    (d * d, 2.0 * d)
                }
            }
            VinaTerm::Hydrophobic { good, bad, .. } => {
                if xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2) {
                    let s = r - optimal_distance(t1, t2);
                    (slope_step(bad, good, s), slope_step_deriv(bad, good, s))
                } else {
                    (0.0, 0.0)
                }
            }
            VinaTerm::NonDirHBond { good, bad, .. } => {
                if xs_h_bond_possible(t1, t2) {
                    let s = r - optimal_distance(t1, t2);
                    (slope_step(bad, good, s), slope_step_deriv(bad, good, s))
                } else {
                    (0.0, 0.0)
                }
            }
            VinaTerm::LinearAttraction { .. } => {
                if crate::atom::constants::is_glued(t1, t2) {
                    (r, 1.0)
                } else {
                    (0.0, 0.0)
                }
            }
        }
    }
}

/// The `SF_VINA` potential set, in order, with their fixed parameters.
/// Hydrogens (XS type `>= XS_TYPE_SIZE`) are handled by the caller, which skips
/// type pairs out of range.
pub const VINA_TERMS: [VinaTerm; 6] = [
    VinaTerm::Gaussian {
        offset: 0.0,
        width: 0.5,
        cutoff: 8.0,
    },
    VinaTerm::Gaussian {
        offset: 3.0,
        width: 2.0,
        cutoff: 8.0,
    },
    VinaTerm::Repulsion {
        offset: 0.0,
        cutoff: 8.0,
    },
    VinaTerm::Hydrophobic {
        good: 0.5,
        bad: 1.5,
        cutoff: 8.0,
    },
    VinaTerm::NonDirHBond {
        good: -0.7,
        bad: 0.0,
        cutoff: 8.0,
    },
    VinaTerm::LinearAttraction { cutoff: 20.0 },
];

/// Number of XS atom types (table dimension for the per-type precalculate).
pub const NUM_XS_TYPES: usize = XS_TYPE_SIZE;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::constants::*;

    #[test]
    fn slope_step_ramps() {
        // bad=1.5, good=0.5 (bad > good): 1 below good, 0 above bad, linear between.
        assert_eq!(slope_step(1.5, 0.5, 0.0), 1.0);
        assert_eq!(slope_step(1.5, 0.5, 2.0), 0.0);
        assert!((slope_step(1.5, 0.5, 1.0) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn optimal_distance_is_sum_of_radii() {
        // C_H radius 1.9, O_A radius 1.7 -> 3.6.
        assert!((optimal_distance(XS_TYPE_C_H, XS_TYPE_O_A) - 3.6).abs() < 1e-12);
        // Glue types -> 0.
        assert_eq!(optimal_distance(XS_TYPE_G0, XS_TYPE_C_H), 0.0);
    }

    #[test]
    fn gaussian_peaks_at_optimal() {
        let g = VinaTerm::Gaussian {
            offset: 0.0,
            width: 0.5,
            cutoff: 8.0,
        };
        let opt = optimal_distance(XS_TYPE_C_H, XS_TYPE_C_H);
        // At r == optimal, gauss(0) == 1.
        assert!((g.eval(XS_TYPE_C_H, XS_TYPE_C_H, opt) - 1.0).abs() < 1e-12);
        // Beyond cutoff -> 0.
        assert_eq!(g.eval(XS_TYPE_C_H, XS_TYPE_C_H, 8.0), 0.0);
    }

    #[test]
    fn repulsion_only_when_closer_than_optimal() {
        let rep = VinaTerm::Repulsion {
            offset: 0.0,
            cutoff: 8.0,
        };
        let opt = optimal_distance(XS_TYPE_C_H, XS_TYPE_C_H); // 3.8
                                                              // Closer than optimal: positive d^2.
        let r = opt - 1.0;
        assert!((rep.eval(XS_TYPE_C_H, XS_TYPE_C_H, r) - 1.0).abs() < 1e-12);
        // Farther than optimal: zero.
        assert_eq!(rep.eval(XS_TYPE_C_H, XS_TYPE_C_H, opt + 1.0), 0.0);
    }

    #[test]
    fn hydrophobic_requires_two_hydrophobic_atoms() {
        let hyd = VinaTerm::Hydrophobic {
            good: 0.5,
            bad: 1.5,
            cutoff: 8.0,
        };
        // C_H is hydrophobic; N_P is not.
        assert_eq!(hyd.eval(XS_TYPE_C_H, XS_TYPE_N_P, 3.0), 0.0);
        // Two hydrophobic, at contact (d <= good) -> 1.
        let opt = optimal_distance(XS_TYPE_C_H, XS_TYPE_C_H);
        assert_eq!(hyd.eval(XS_TYPE_C_H, XS_TYPE_C_H, opt + 0.5), 1.0);
    }

    #[test]
    fn h_bond_requires_donor_acceptor() {
        let hb = VinaTerm::NonDirHBond {
            good: -0.7,
            bad: 0.0,
            cutoff: 8.0,
        };
        // N_D (donor) + O_A (acceptor): possible.
        let opt = optimal_distance(XS_TYPE_N_D, XS_TYPE_O_A);
        assert_eq!(hb.eval(XS_TYPE_N_D, XS_TYPE_O_A, opt - 0.7), 1.0);
        // Two carbons: not possible.
        assert_eq!(hb.eval(XS_TYPE_C_H, XS_TYPE_C_H, opt), 0.0);
    }
}
