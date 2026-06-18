// SPDX-License-Identifier: Apache-2.0
//! The Vina scoring function: weighted potential terms, the precalculate
//! "fast" energy lookup, the `curl` energy cap, and the conformation-independent
//! torsional term.

pub mod eval;
pub mod grid;
pub mod terms;

use crate::atom::constants::XS_TYPE_SIZE;
use crate::math::{sqr, EPSILON_FL};
use terms::{VinaTerm, VINA_TERMS};

/// True when `x` is comfortably below the `f64::MAX` sentinel Vina uses to mean
/// "no cap".
#[inline]
fn not_max(x: f64) -> bool {
    x < 0.1 * f64::MAX
}

/// `curl(e, v)` — the soft-curl energy cap. Only positive energies are
/// scaled, by `v / (v + e)`; attractive energies pass through unchanged.
#[inline]
pub fn curl(e: &mut f64, v: f64) {
    if *e > 0.0 && not_max(v) {
        let tmp = if v < EPSILON_FL { 0.0 } else { v / (v + *e) };
        *e *= tmp;
    }
}

/// `curl(e, deriv, v)` — the vector form: scales a positive energy by
/// `v/(v+e)` and its derivative by the square of that factor.
#[inline]
pub fn curl_deriv(e: &mut f64, deriv: &mut crate::math::Vec3, v: f64) {
    if *e > 0.0 && not_max(v) {
        let tmp = if v < EPSILON_FL { 0.0 } else { v / (v + *e) };
        *e *= tmp;
        *deriv = (tmp * tmp) * *deriv;
    }
}

/// `x / y`, but treats near-zero operands specially: a near-zero numerator
/// yields 0, and a near-zero denominator saturates to `±f64::MAX` (signed by
/// the operands) instead of producing an infinity.
#[inline]
fn conf_smooth_div(x: f64, y: f64) -> f64 {
    if x.abs() < EPSILON_FL {
        0.0
    } else if y.abs() < EPSILON_FL {
        if x * y > 0.0 {
            f64::MAX
        } else {
            -f64::MAX
        }
    } else {
        x / y
    }
}

/// The precalculate radial sampling factor (`factor = 32`).
pub const FACTOR: f64 = 32.0;

/// Vina default weights for the six potentials plus the conf-independent term,
/// in order: gauss1, gauss2, repulsion, hydrophobic, h-bond, glue, conf_indep.
///
/// The conf-independent weight is `5 * weight_rot / 0.1 - 1` with
/// `weight_rot = 0.05846`, matching Vina's weights.
pub fn vina_default_weights() -> [f64; 7] {
    let weight_rot = 0.05846;
    [
        -0.035579, // gauss1
        -0.005156, // gauss2
        0.840245,  // repulsion
        -0.035069, // hydrophobic
        -0.587439, // h-bond
        50.0,      // glue (linearattraction)
        5.0 * weight_rot / 0.1 - 1.0,
    ]
}

/// The Vina scoring function (`ScoringFunction` for `SF_VINA`).
#[derive(Debug, Clone)]
pub struct ScoringFunction {
    terms: [VinaTerm; 6],
    /// Length 7: one weight per potential, then the conf-independent weight.
    weights: [f64; 7],
    cutoff: f64,
    max_cutoff: f64,
}

impl Default for ScoringFunction {
    fn default() -> Self {
        Self::vina()
    }
}

impl ScoringFunction {
    /// The standard Vina scoring function with default weights.
    pub fn vina() -> Self {
        ScoringFunction {
            terms: VINA_TERMS,
            weights: vina_default_weights(),
            cutoff: 8.0,
            max_cutoff: 20.0,
        }
    }

    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }
    pub fn cutoff_sqr(&self) -> f64 {
        sqr(self.cutoff)
    }
    pub fn max_cutoff(&self) -> f64 {
        self.max_cutoff
    }

    /// The raw weighted sum of potentials (no cutoff check; each term self-zeros
    /// beyond its cutoff). Type indices must be valid XS types.
    #[inline]
    pub fn eval_raw(&self, t1: usize, t2: usize, r: f64) -> f64 {
        let mut acc = 0.0;
        for i in 0..6 {
            acc += self.weights[i] * self.terms[i].eval(t1, t2, r);
        }
        acc
    }

    /// The "fast" binned energy lookup computed on the fly: the binned
    /// average `(eval_raw(rs[i]) + eval_raw(rs[i+1])) / 2` with `i = ⌊factor·r²⌋`
    /// and `rs[k] = √(k/factor)`. Bit-identical to the precomputed [`grid`]
    /// table value, which is what makes the on-the-fly and tabulated paths
    /// interchangeable.
    #[inline]
    pub fn eval_fast(&self, t1: usize, t2: usize, r2: f64) -> f64 {
        if t1 >= XS_TYPE_SIZE || t2 >= XS_TYPE_SIZE {
            return 0.0;
        }
        let i = (FACTOR * r2) as usize;
        let rs_i = (i as f64 / FACTOR).sqrt();
        let rs_i1 = ((i + 1) as f64 / FACTOR).sqrt();
        (self.eval_raw(t1, t2, rs_i) + self.eval_raw(t1, t2, rs_i1)) / 2.0
    }

    /// The weighted energy at radial sample `k` (`rs[k] = sqrt(k/factor)`).
    #[inline]
    fn smooth_first(&self, t1: usize, t2: usize, k: usize) -> f64 {
        let r = (k as f64 / FACTOR).sqrt();
        self.eval_raw(t1, t2, r)
    }

    /// `dor` (derivative-over-r) at radial sample `k`: central difference of
    /// the smoothed energy, scaled; zero at `k == 0` (and at the top boundary,
    /// which is unreachable for `r2 < cutoff`).
    #[inline]
    fn dor(&self, t1: usize, t2: usize, k: usize) -> f64 {
        if k == 0 {
            return 0.0;
        }
        let rk = (k as f64 / FACTOR).sqrt();
        let r_next = ((k + 1) as f64 / FACTOR).sqrt();
        let r_prev = ((k - 1) as f64 / FACTOR).sqrt();
        let delta = r_next - r_prev;
        (self.smooth_first(t1, t2, k + 1) - self.smooth_first(t1, t2, k - 1)) / (delta * rk)
    }

    /// Linear interpolation of `(energy, derivative-over-r)` between radial
    /// samples `i1 = floor(factor·r²)` and `i1+1`. The energy here is the
    /// **interpolated** `smooth` value (used by the optimizer), distinct from the
    /// binned [`Self::eval_fast`] used for scoring.
    #[inline]
    pub fn eval_deriv(&self, t1: usize, t2: usize, r2: f64) -> (f64, f64) {
        if t1 >= XS_TYPE_SIZE || t2 >= XS_TYPE_SIZE {
            return (0.0, 0.0);
        }
        let r2_factored = FACTOR * r2;
        let i1 = r2_factored as usize;
        let i2 = i1 + 1;
        let rem = r2_factored - i1 as f64;
        let e1 = self.smooth_first(t1, t2, i1);
        let e2 = self.smooth_first(t1, t2, i2);
        let d1 = self.dor(t1, t2, i1);
        let d2 = self.dor(t1, t2, i2);
        let e = e1 + rem * (e2 - e1);
        let dor = d1 + rem * (d2 - d1);
        (e, dor)
    }

    /// Exact continuous `(energy, d energy / d r)` — the weighted analytic
    /// potential and its derivative (no binning/interpolation). Used to validate
    /// the interpolated derivative path via finite differences.
    #[inline]
    pub fn eval_exact_deriv(&self, t1: usize, t2: usize, r: f64) -> (f64, f64) {
        if t1 >= XS_TYPE_SIZE || t2 >= XS_TYPE_SIZE {
            return (0.0, 0.0);
        }
        let mut e = 0.0;
        let mut de = 0.0;
        for i in 0..6 {
            let (v, d) = self.terms[i].eval_deriv_r(t1, t2, r);
            e += self.weights[i] * v;
            de += self.weights[i] * d;
        }
        (e, de)
    }

    /// `conf_independent` for Vina: the `num_tors_div` term applied to energy
    /// `e`. `weight = 0.1 * (w + 1)`; returns `conf_smooth_div(e, 1 + weight *
    /// num_tors / 5)`.
    pub fn conf_independent(&self, num_tors: f64, e: f64) -> f64 {
        let weight = 0.1 * (self.weights[6] + 1.0);
        conf_smooth_div(e, 1.0 + weight * num_tors / 5.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn curl_leaves_attractive_unchanged_and_softens_repulsive() {
        let mut e = -3.0;
        curl(&mut e, 1000.0);
        assert_eq!(e, -3.0);

        let mut e = 2.0;
        curl(&mut e, 1000.0);
        // 2 * 1000/1002
        assert!((e - 2.0 * 1000.0 / 1002.0).abs() < 1e-12);
    }

    #[test]
    fn conf_independent_matches_known_1iep_relation() {
        // For 1iep: inter = -17.634, num_tors = 7 -> EFE = -12.513.
        let sf = ScoringFunction::vina();
        let efe = sf.conf_independent(7.0, -17.634);
        assert!((efe - (-12.513)).abs() < 1e-3, "EFE = {efe}");
        // Torsional free energy (3) = EFE - inter.
        let torsional = efe - (-17.634);
        assert!((torsional - 5.121).abs() < 1e-3, "torsional = {torsional}");
    }

    #[test]
    fn conf_independent_weight_is_1_923() {
        let w = vina_default_weights();
        assert!((w[6] - 1.923).abs() < 1e-9);
    }
}
