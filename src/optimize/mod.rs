// SPDX-License-Identifier: Apache-2.0
//! BFGS local optimization.
//!
//! Minimizes the energy over conformation space using a BFGS approximation of
//! the inverse Hessian with a backtracking (Armijo) line search. The objective
//! `f(conf, &mut grad) -> energy` evaluates the energy at a conformation and
//! fills the conformation-space gradient (see [`crate::scoring::eval`]).

use crate::math::EPSILON_FL;
use crate::model::conf::{Change, Conf};
use crate::model::matrix::TriangularMatrix;

type FlMat = TriangularMatrix<f64>;

fn scalar_product(a: &Change, b: &Change, n: usize) -> f64 {
    let mut tmp = 0.0;
    for i in 0..n {
        tmp += a.get(i) * b.get(i);
    }
    tmp
}

/// `out = -m * in` (m symmetric).
fn minus_mat_vec_product(m: &FlMat, input: &Change, out: &mut Change) {
    let n = m.dim();
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..n {
            sum += *m.get_sym(i, j) * input.get(j);
        }
        *out.get_mut(i) = -sum;
    }
}

fn set_diagonal(m: &mut FlMat, x: f64) {
    for i in 0..m.dim() {
        *m.get_mut(i, i) = x;
    }
}

/// `b -= a`.
fn subtract_change(b: &mut Change, a: &Change, n: usize) {
    for i in 0..n {
        *b.get_mut(i) -= a.get(i);
    }
}

/// Rank-2 BFGS update of the inverse-Hessian approximation `h`.
fn bfgs_update(h: &mut FlMat, p: &Change, y: &Change, alpha: f64) -> bool {
    let n = h.dim();
    let yp = scalar_product(y, p, n);
    if alpha * yp < EPSILON_FL {
        return false;
    }
    let mut minus_hy = y.clone();
    minus_mat_vec_product(h, y, &mut minus_hy);
    let yhy = -scalar_product(y, &minus_hy, n);
    let r = 1.0 / (alpha * yp);
    for i in 0..n {
        for j in i..n {
            *h.get_mut(i, j) +=
                alpha * r * (minus_hy.get(i) * p.get(j) + minus_hy.get(j) * p.get(i))
                    + alpha * alpha * (r * r * yhy + r) * p.get(i) * p.get(j);
        }
    }
    true
}

/// Result of the backtracking line search.
struct LineSearch {
    alpha: f64,
    x_new: Conf,
    g_new: Change,
    f1: f64,
}

/// Backtracking Armijo line search along direction `p`.
#[allow(clippy::too_many_arguments)]
fn line_search<F>(
    f: &mut F,
    n: usize,
    x: &Conf,
    g: &Change,
    f0: f64,
    p: &Change,
    g_new_template: &Change,
    evalcount: &mut i64,
) -> LineSearch
where
    F: FnMut(&Conf, &mut Change) -> f64,
{
    const C0: f64 = 0.0001;
    const MAX_TRIALS: usize = 10;
    const MULTIPLIER: f64 = 0.5;
    let mut alpha = 1.0;

    let pg = scalar_product(p, g, n);

    let mut x_new = x.clone();
    let mut g_new = g_new_template.clone();
    let mut f1 = f0;
    for _ in 0..MAX_TRIALS {
        x_new = x.clone();
        x_new.increment(p, alpha);
        f1 = f(&x_new, &mut g_new);
        *evalcount += 1;
        if f1 - f0 < C0 * alpha * pg {
            break;
        }
        alpha *= MULTIPLIER;
    }
    LineSearch {
        alpha,
        x_new,
        g_new,
        f1,
    }
}

/// Minimize `f` from `x` (updated in place), with `g` holding the working
/// gradient. Returns the final energy.
// The `!(a >= b)` / `!(a <= b)` comparisons are deliberate: they also fire for
// NaN (where `a < b` would not), matching the convergence/restore guards.
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn bfgs<F>(
    f: &mut F,
    x: &mut Conf,
    g: &mut Change,
    max_steps: usize,
    evalcount: &mut i64,
) -> f64
where
    F: FnMut(&Conf, &mut Change) -> f64,
{
    let n = g.num_floats();
    let mut h = FlMat::new(n, 0.0);
    set_diagonal(&mut h, 1.0);

    let mut f0 = f(x, g);
    *evalcount += 1;

    let f_orig = f0;
    let g_orig = g.clone();
    let x_orig = x.clone();

    let mut p = g.clone();

    for step in 0..max_steps {
        minus_mat_vec_product(&h, g, &mut p); // p = -h * g
        let ls = line_search(f, n, x, g, f0, &p, g, evalcount);
        let alpha = ls.alpha;
        let mut y = ls.g_new.clone();
        subtract_change(&mut y, g, n); // y = g_new - g

        f0 = ls.f1;
        *x = ls.x_new;
        // Convergence test uses the pre-step gradient; breaks for NaNs too.
        if !(scalar_product(g, g, n).sqrt() >= 1e-5) {
            break;
        }
        *g = ls.g_new;

        if step == 0 {
            let yy = scalar_product(&y, &y, n);
            if yy.abs() > EPSILON_FL {
                set_diagonal(&mut h, alpha * scalar_product(&y, &p, n) / yy);
            }
        }

        bfgs_update(&mut h, &p, &y, alpha);
    }

    // If we somehow made things worse (or hit NaN), restore the start.
    if !(f0 <= f_orig) {
        f0 = f_orig;
        *x = x_orig;
        *g = g_orig;
    }
    f0
}
