// SPDX-License-Identifier: Apache-2.0
//! Precomputed affinity grids (the "cache").
//!
//! For each ligand XS atom type, a 3D grid samples the intermolecular energy at
//! `granularity`-spaced points (each sample is the exact pairwise sum over
//! receptor atoms). At evaluation, a ligand atom's energy and force come from
//! trilinear interpolation of its grid (the force is the interpolation's spatial
//! gradient, so energy and force are consistent). This is the field the BFGS
//! optimizer and Monte-Carlo search follow.
//!
//! The trilinear gradient formula and the small 3-element coordinate loops keep
//! their explicit form, so two clippy lints are allowed module-wide.
#![allow(clippy::neg_multiply, clippy::needless_range_loop)]

use crate::atom::constants::XS_TYPE_SIZE;
use crate::math::{Vec3, ZERO};
use crate::model::matrix::triangular_index_permissive;
use crate::model::Model;

use super::eval::SearchBox;
use super::{curl, curl_deriv, ScoringFunction, FACTOR};

/// Map the CG variants onto their base carbon grid and reject glue dummies.
/// Returns `None` if the atom has no grid (glue or out-of-range XS type).
fn grid_type(xs: usize) -> Option<usize> {
    use crate::atom::constants::*;
    if xs >= XS_TYPE_SIZE || is_glue_type(xs) {
        return None;
    }
    Some(match xs {
        XS_TYPE_C_H_CG0 | XS_TYPE_C_H_CG1 | XS_TYPE_C_H_CG2 | XS_TYPE_C_H_CG3 => XS_TYPE_C_H,
        XS_TYPE_C_P_CG0 | XS_TYPE_C_P_CG1 | XS_TYPE_C_P_CG2 | XS_TYPE_C_P_CG3 => XS_TYPE_C_P,
        other => other,
    })
}

/// Tabulated fast-energy lookup for the grid build: binned-average energy
/// per XS type pair over `r2 in [0, cutoff_sqr]`.
struct FastTable {
    /// `fast[type_pair_index][bin]`.
    fast: Vec<Vec<f64>>,
    max_bin: usize,
}

impl FastTable {
    fn build(sf: &ScoringFunction) -> Self {
        let cutoff_sqr = sf.cutoff_sqr();
        let max_bin = (FACTOR * cutoff_sqr) as usize + 1;
        let n_types = XS_TYPE_SIZE;
        let n_pairs = n_types * (n_types + 1) / 2;
        let mut fast = vec![Vec::new(); n_pairs];
        let smooth =
            |t1: usize, t2: usize, k: usize| sf.eval_raw(t1, t2, (k as f64 / FACTOR).sqrt());
        for t2 in 0..n_types {
            for t1 in 0..=t2 {
                let idx = triangular_index_permissive(t1, t2);
                let mut col = vec![0.0; max_bin + 1];
                for (bin, slot) in col.iter_mut().enumerate() {
                    *slot = (smooth(t1, t2, bin) + smooth(t1, t2, bin + 1)) / 2.0;
                }
                fast[idx] = col;
            }
        }
        FastTable { fast, max_bin }
    }

    #[inline]
    fn eval_fast(&self, t1: usize, t2: usize, r2: f64) -> f64 {
        let bin = (FACTOR * r2) as usize;
        if bin > self.max_bin {
            return 0.0;
        }
        self.fast[triangular_index_permissive(t1, t2)][bin]
    }
}

/// A single 3D affinity grid.
#[derive(Debug, Clone)]
struct Grid {
    init: [f64; 3],
    factor: [f64; 3],
    factor_inv: [f64; 3],
    dim_fl_minus_1: [f64; 3],
    dim: [usize; 3], // sample points per axis (= n_voxels + 1)
    data: Vec<f64>,
}

impl Grid {
    fn new(gbox: &SearchBox) -> Self {
        let dim = [
            gbox.n_voxels[0] + 1,
            gbox.n_voxels[1] + 1,
            gbox.n_voxels[2] + 1,
        ];
        let mut factor = [0.0; 3];
        let mut factor_inv = [0.0; 3];
        let mut dim_fl_minus_1 = [0.0; 3];
        for i in 0..3 {
            let range = gbox.span(i);
            dim_fl_minus_1[i] = (dim[i] - 1) as f64;
            factor[i] = dim_fl_minus_1[i] / range;
            factor_inv[i] = 1.0 / factor[i];
        }
        Grid {
            init: gbox.begin,
            factor,
            factor_inv,
            dim_fl_minus_1,
            dim,
            data: vec![0.0; dim[0] * dim[1] * dim[2]],
        }
    }

    #[inline]
    fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        x + self.dim[0] * (y + self.dim[1] * z)
    }

    fn index_to_argument(&self, x: usize, y: usize, z: usize) -> Vec3 {
        Vec3::new(
            self.init[0] + self.factor_inv[0] * x as f64,
            self.init[1] + self.factor_inv[1] * y as f64,
            self.init[2] + self.factor_inv[2] * z as f64,
        )
    }

    /// Trilinear interpolation with out-of-box penalty and (optionally) the
    /// spatial gradient.
    fn evaluate_aux(&self, location: &Vec3, slope: f64, v: f64, deriv: Option<&mut Vec3>) -> f64 {
        let mut s = [0.0; 3];
        for i in 0..3 {
            s[i] = (location[i] - self.init[i]) * self.factor[i];
        }
        let mut miss = [0.0; 3];
        let mut region = [0i32; 3];
        let mut a = [0usize; 3];
        for i in 0..3 {
            if s[i] < 0.0 {
                miss[i] = -s[i];
                region[i] = -1;
                a[i] = 0;
                s[i] = 0.0;
            } else if s[i] >= self.dim_fl_minus_1[i] {
                miss[i] = s[i] - self.dim_fl_minus_1[i];
                region[i] = 1;
                a[i] = self.dim[i] - 2;
                s[i] = 1.0;
            } else {
                region[i] = 0;
                a[i] = s[i] as usize;
                s[i] -= a[i] as f64;
            }
        }
        let penalty = slope
            * (miss[0] * self.factor_inv[0]
                + miss[1] * self.factor_inv[1]
                + miss[2] * self.factor_inv[2]);

        let (x0, y0, z0) = (a[0], a[1], a[2]);
        let (x1, y1, z1) = (x0 + 1, y0 + 1, z0 + 1);
        let f000 = self.data[self.idx(x0, y0, z0)];
        let f100 = self.data[self.idx(x1, y0, z0)];
        let f010 = self.data[self.idx(x0, y1, z0)];
        let f110 = self.data[self.idx(x1, y1, z0)];
        let f001 = self.data[self.idx(x0, y0, z1)];
        let f101 = self.data[self.idx(x1, y0, z1)];
        let f011 = self.data[self.idx(x0, y1, z1)];
        let f111 = self.data[self.idx(x1, y1, z1)];

        let (x, y, z) = (s[0], s[1], s[2]);
        let (mx, my, mz) = (1.0 - x, 1.0 - y, 1.0 - z);

        let mut f = f000 * mx * my * mz
            + f100 * x * my * mz
            + f010 * mx * y * mz
            + f110 * x * y * mz
            + f001 * mx * my * z
            + f101 * x * my * z
            + f011 * mx * y * z
            + f111 * x * y * z;

        match deriv {
            Some(d) => {
                let x_g = f000 * (-1.0) * my * mz
                    + f100 * my * mz
                    + f010 * (-1.0) * y * mz
                    + f110 * y * mz
                    + f001 * (-1.0) * my * z
                    + f101 * my * z
                    + f011 * (-1.0) * y * z
                    + f111 * y * z;
                let y_g = f000 * mx * (-1.0) * mz
                    + f100 * x * (-1.0) * mz
                    + f010 * mx * mz
                    + f110 * x * mz
                    + f001 * mx * (-1.0) * z
                    + f101 * x * (-1.0) * z
                    + f011 * mx * z
                    + f111 * x * z;
                let z_g = f000 * mx * my * (-1.0)
                    + f100 * x * my * (-1.0)
                    + f010 * mx * y * (-1.0)
                    + f110 * x * y * (-1.0)
                    + f001 * mx * my
                    + f101 * x * my
                    + f011 * mx * y
                    + f111 * x * y;
                let mut gradient = Vec3::new(x_g, y_g, z_g);
                curl_deriv(&mut f, &mut gradient, v);
                for i in 0..3 {
                    let gradient_everywhere = if region[i] == 0 { gradient[i] } else { 0.0 };
                    d[i] = self.factor[i] * gradient_everywhere + slope * region[i] as f64;
                }
                f + penalty
            }
            None => {
                curl(&mut f, v);
                f + penalty
            }
        }
    }
}

/// The precomputed affinity grids for all needed ligand atom types.
pub struct Cache {
    grids: Vec<Option<Grid>>, // indexed by XS type
    slope: f64,
}

impl Cache {
    /// Build the grids needed by the model's movable atoms.
    pub fn populate(model: &Model, sf: &ScoringFunction, gbox: &SearchBox) -> Self {
        let table = FastTable::build(sf);
        let cutoff_sqr = sf.cutoff_sqr();
        let mut grids: Vec<Option<Grid>> = (0..XS_TYPE_SIZE).map(|_| None).collect();

        // Needed grid types: distinct grid_type() of the movable atoms.
        let mut needed: Vec<usize> = Vec::new();
        for i in 0..model.num_movable_atoms() {
            if let Some(t) = grid_type(model.atoms[i].xs()) {
                if !needed.contains(&t) {
                    needed.push(t);
                    grids[t] = Some(Grid::new(gbox));
                }
            }
        }
        if needed.is_empty() {
            return Cache {
                grids,
                slope: gbox.slope,
            };
        }

        // Spatial bins over receptor atoms (cell size = cutoff) to keep the
        // per-sample-point receptor scan local.
        let cutoff = cutoff_sqr.sqrt();
        let bins = ReceptorBins::build(model, cutoff);

        let dim = grids[needed[0]].as_ref().unwrap().dim;
        let proto = grids[needed[0]].as_ref().unwrap().clone();
        let mut affinities = vec![0.0; needed.len()];
        for x in 0..dim[0] {
            for y in 0..dim[1] {
                for z in 0..dim[2] {
                    for a in affinities.iter_mut() {
                        *a = 0.0;
                    }
                    let probe = proto.index_to_argument(x, y, z);
                    bins.for_each_near(&probe, |b| {
                        let t1 = model.grid_atoms[b].xs();
                        if t1 >= XS_TYPE_SIZE {
                            return;
                        }
                        let r2 = crate::math::distance_sqr(&model.grid_atoms[b].coords, &probe);
                        if r2 <= cutoff_sqr {
                            for (j, &t2) in needed.iter().enumerate() {
                                affinities[j] += table.eval_fast(t1, t2, r2);
                            }
                        }
                    });
                    for (j, &t) in needed.iter().enumerate() {
                        let g = grids[t].as_mut().unwrap();
                        let i = g.idx(x, y, z);
                        g.data[i] = affinities[j];
                    }
                }
            }
        }
        Cache {
            grids,
            slope: gbox.slope,
        }
    }

    /// Intermolecular energy by grid interpolation.
    pub fn eval(&self, model: &Model, v: f64) -> f64 {
        let mut e = 0.0;
        for i in 0..model.num_movable_atoms() {
            if let Some(t) = grid_type(model.atoms[i].xs()) {
                let g = self.grids[t].as_ref().expect("grid for needed type");
                e += g.evaluate_aux(&model.coords[i], self.slope, v, None);
            }
        }
        e
    }

    /// Intermolecular energy + per-atom forces by grid interpolation.
    /// Writes `forces[i]` for every movable atom.
    pub fn eval_deriv(&self, model: &Model, v: f64, forces: &mut [Vec3]) -> f64 {
        let mut e = 0.0;
        for i in 0..model.num_movable_atoms() {
            match grid_type(model.atoms[i].xs()) {
                Some(t) => {
                    let g = self.grids[t].as_ref().expect("grid for needed type");
                    let mut deriv = ZERO;
                    e += g.evaluate_aux(&model.coords[i], self.slope, v, Some(&mut deriv));
                    forces[i] = deriv;
                }
                None => forces[i] = ZERO,
            }
        }
        e
    }
}

/// A uniform spatial hash of receptor atoms, for the grid build.
struct ReceptorBins {
    cell: f64,
    origin: [f64; 3],
    dims: [usize; 3],
    cells: Vec<Vec<usize>>, // atom indices per cell
}

impl ReceptorBins {
    fn build(model: &Model, cell: f64) -> Self {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for a in &model.grid_atoms {
            for j in 0..3 {
                lo[j] = lo[j].min(a.coords[j]);
                hi[j] = hi[j].max(a.coords[j]);
            }
        }
        if model.grid_atoms.is_empty() {
            lo = [0.0; 3];
            hi = [0.0; 3];
        }
        let mut dims = [1usize; 3];
        for j in 0..3 {
            dims[j] = (((hi[j] - lo[j]) / cell).floor() as usize) + 1;
        }
        let mut cells = vec![Vec::new(); dims[0] * dims[1] * dims[2]];
        let cell_index = |c: [usize; 3]| c[0] + dims[0] * (c[1] + dims[1] * c[2]);
        for (idx, a) in model.grid_atoms.iter().enumerate() {
            let mut c = [0usize; 3];
            for j in 0..3 {
                c[j] = (((a.coords[j] - lo[j]) / cell).floor() as usize).min(dims[j] - 1);
            }
            cells[cell_index(c)].push(idx);
        }
        ReceptorBins {
            cell,
            origin: lo,
            dims,
            cells,
        }
    }

    /// Visit every receptor atom in the cell containing `p` and its 26 neighbours
    /// (a superset of all atoms within one cell width of `p`).
    fn for_each_near(&self, p: &Vec3, mut f: impl FnMut(usize)) {
        let mut base = [0i64; 3];
        for j in 0..3 {
            base[j] = ((p[j] - self.origin[j]) / self.cell).floor() as i64;
        }
        for dz in -1..=1i64 {
            for dy in -1..=1i64 {
                for dx in -1..=1i64 {
                    let c = [base[0] + dx, base[1] + dy, base[2] + dz];
                    if c.iter()
                        .enumerate()
                        .any(|(j, &v)| v < 0 || v as usize >= self.dims[j])
                    {
                        continue;
                    }
                    let ci = c[0] as usize
                        + self.dims[0] * (c[1] as usize + self.dims[1] * c[2] as usize);
                    for &atom in &self.cells[ci] {
                        f(atom);
                    }
                }
            }
        }
    }
}
