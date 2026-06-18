// SPDX-License-Identifier: Apache-2.0
//! Seeded pseudo-random number generation for the Monte-Carlo search.
//!
//! Reproducing the exact distribution byte stream of a particular Boost build is
//! brittle and version-dependent, and is explicitly **not** a goal: the
//! Monte-Carlo trajectory is validated by statistical equivalence and
//! convergence, not bit-for-bit.
//!
//! We therefore use a clean, well-defined, seedable generator: the standard
//! 32-bit Mersenne Twister (the same engine family Vina uses) with textbook
//! distribution transforms. Given a seed the search is fully deterministic.

use crate::math::{Quat, Vec3, PI};

/// Standard MT19937 (32-bit Mersenne Twister).
#[derive(Debug, Clone)]
struct Mt19937 {
    mt: [u32; 624],
    index: usize,
}

const N: usize = 624;
const M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;

impl Mt19937 {
    fn new(seed: u32) -> Self {
        let mut mt = [0u32; N];
        mt[0] = seed;
        for i in 1..N {
            mt[i] = 1_812_433_253u32
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Mt19937 { mt, index: N }
    }

    fn generate(&mut self) {
        for i in 0..N {
            let y = (self.mt[i] & UPPER_MASK) | (self.mt[(i + 1) % N] & LOWER_MASK);
            let mut next = self.mt[(i + M) % N] ^ (y >> 1);
            if y & 1 != 0 {
                next ^= MATRIX_A;
            }
            self.mt[i] = next;
        }
        self.index = 0;
    }

    fn next_u32(&mut self) -> u32 {
        if self.index >= N {
            self.generate();
        }
        let mut y = self.mt[self.index];
        self.index += 1;
        // Tempering.
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }
}

/// The engine's RNG: an MT19937 plus a cached normal deviate (Box-Muller
/// produces two at a time).
#[derive(Debug, Clone)]
pub struct Rng {
    engine: Mt19937,
    normal_cache: Option<f64>,
}

impl Rng {
    /// Seed the generator.
    pub fn seed(seed: u32) -> Self {
        Rng {
            engine: Mt19937::new(seed),
            normal_cache: None,
        }
    }

    /// Uniform `f64` in `[0, 1)`.
    #[inline]
    fn unit(&mut self) -> f64 {
        // 32-bit resolution in [0, 1).
        self.engine.next_u32() as f64 / (u32::MAX as f64 + 1.0)
    }

    /// Uniform in `[a, b)` (a < b).
    pub fn random_fl(&mut self, a: f64, b: f64) -> f64 {
        debug_assert!(a < b);
        a + (b - a) * self.unit()
    }

    /// Uniform integer in `[a, b]` (a <= b).
    pub fn random_int(&mut self, a: i32, b: i32) -> i32 {
        debug_assert!(a <= b);
        // Width is computed in i64 so the full i32 range cannot overflow.
        let span = (b as i64 - a as i64 + 1) as u64;
        (a as i64 + (self.engine.next_u32() as u64 % span) as i64) as i32
    }

    /// Uniform `usize` in `[a, b]`.
    pub fn random_sz(&mut self, a: usize, b: usize) -> usize {
        self.random_int(a as i32, b as i32) as usize
    }

    /// Normal deviate via Box-Muller (caching the second deviate).
    pub fn random_normal(&mut self, mean: f64, sigma: f64) -> f64 {
        if let Some(z) = self.normal_cache.take() {
            return mean + sigma * z;
        }
        // Draw two uniforms in (0, 1].
        let mut u1 = self.unit();
        if u1 < 1e-300 {
            u1 = 1e-300;
        }
        let u2 = self.unit();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * PI * u2;
        self.normal_cache = Some(r * theta.sin());
        mean + sigma * (r * theta.cos())
    }

    /// A point uniformly inside the unit sphere (rejection sampling).
    pub fn random_inside_sphere(&mut self) -> Vec3 {
        loop {
            let v = Vec3::new(
                self.random_fl(-1.0, 1.0),
                self.random_fl(-1.0, 1.0),
                self.random_fl(-1.0, 1.0),
            );
            if v.norm_sqr() < 1.0 {
                return v;
            }
        }
    }

    /// Uniform inside the box.
    pub fn random_in_box(&mut self, corner1: &Vec3, corner2: &Vec3) -> Vec3 {
        Vec3::new(
            self.random_fl(corner1[0], corner2[0]),
            self.random_fl(corner1[1], corner2[1]),
            self.random_fl(corner1[2], corner2[2]),
        )
    }

    /// A uniformly random unit quaternion (four normal deviates, normalized).
    pub fn random_orientation(&mut self) -> Quat {
        loop {
            let q = Quat::new(
                self.random_normal(0.0, 1.0),
                self.random_normal(0.0, 1.0),
                self.random_normal(0.0, 1.0),
                self.random_normal(0.0, 1.0),
            );
            let nrm = q.norm_sqr().sqrt();
            if nrm > crate::math::EPSILON_FL {
                let inv = 1.0 / nrm;
                return Quat::new(q.w * inv, q.x * inv, q.y * inv, q.z * inv);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mt19937_known_vector() {
        // Reference output of MT19937 seeded with 5489 (the standard default):
        // the first generated 32-bit word is 3499211612.
        let mut e = Mt19937::new(5489);
        assert_eq!(e.next_u32(), 3_499_211_612);
        assert_eq!(e.next_u32(), 581_869_302);
    }

    #[test]
    fn deterministic_for_same_seed() {
        let mut a = Rng::seed(42);
        let mut b = Rng::seed(42);
        for _ in 0..100 {
            assert_eq!(a.engine.next_u32(), b.engine.next_u32());
        }
    }

    #[test]
    fn random_fl_in_range_and_roughly_uniform() {
        let mut r = Rng::seed(7);
        let mut sum = 0.0;
        let n = 100_000;
        for _ in 0..n {
            let x = r.random_fl(-2.0, 5.0);
            assert!((-2.0..5.0).contains(&x));
            sum += x;
        }
        let mean = sum / n as f64;
        // Expected mean 1.5; allow sampling slack.
        assert!((mean - 1.5).abs() < 0.05, "mean {mean}");
    }

    #[test]
    fn random_int_covers_endpoints() {
        let mut r = Rng::seed(123);
        let (mut lo, mut hi) = (false, false);
        for _ in 0..1000 {
            let x = r.random_int(0, 3);
            assert!((0..=3).contains(&x));
            if x == 0 {
                lo = true;
            }
            if x == 3 {
                hi = true;
            }
        }
        assert!(lo && hi);
    }

    #[test]
    fn random_normal_stats() {
        let mut r = Rng::seed(99);
        let n = 200_000;
        let mut sum = 0.0;
        let mut sumsq = 0.0;
        for _ in 0..n {
            let z = r.random_normal(1.0, 2.0);
            sum += z;
            sumsq += z * z;
        }
        let mean = sum / n as f64;
        let var = sumsq / n as f64 - mean * mean;
        assert!((mean - 1.0).abs() < 0.03, "mean {mean}");
        assert!((var - 4.0).abs() < 0.1, "var {var}");
    }

    #[test]
    fn orientation_is_unit() {
        let mut r = Rng::seed(3);
        for _ in 0..1000 {
            let q = r.random_orientation();
            assert!((q.norm_sqr() - 1.0).abs() < 1e-12);
        }
    }
}
