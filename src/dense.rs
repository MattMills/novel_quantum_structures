//! Exact dense state-vector simulation over an arbitrary mixed-radix qudit
//! chain. Exponential in chain length — used as ground truth to verify the
//! tensor-network backend on small systems.
//!
//! Basis convention: site 0 is the most-significant digit, so the amplitude
//! index of `|x_0, x_1, …, x_{n-1}⟩` is `((x_0·d_1 + x_1)·d_2 + x_2)…`.

use crate::c64::C64;
use crate::mat::Mat;

/// Dense state vector over sites with per-site dimensions `dims`.
#[derive(Clone)]
pub struct DenseState {
    pub dims: Vec<usize>,
    pub amps: Vec<C64>,
}

/// Product of all local dimensions (total Hilbert-space dimension).
pub fn total_dim(dims: &[usize]) -> usize {
    dims.iter().product()
}

/// Stride of each site in the flattened index (site 0 most significant).
pub fn strides(dims: &[usize]) -> Vec<usize> {
    let n = dims.len();
    let mut s = vec![1usize; n];
    for i in (0..n.saturating_sub(1)).rev() {
        s[i] = s[i + 1] * dims[i + 1];
    }
    s
}

impl DenseState {
    /// Computational basis state `|digits⟩`.
    pub fn basis_state(dims: &[usize], digits: &[usize]) -> DenseState {
        assert_eq!(dims.len(), digits.len());
        for (d, x) in dims.iter().zip(digits) {
            assert!(x < d, "digit {} out of range for dimension {}", x, d);
        }
        let n = total_dim(dims);
        let mut amps = vec![C64::ZERO; n];
        let idx: usize = digits
            .iter()
            .zip(strides(dims).iter())
            .map(|(x, s)| x * s)
            .sum();
        amps[idx] = C64::ONE;
        DenseState {
            dims: dims.to_vec(),
            amps,
        }
    }

    /// All-zeros basis state.
    pub fn zero_state(dims: &[usize]) -> DenseState {
        DenseState::basis_state(dims, &vec![0; dims.len()])
    }

    pub fn num_sites(&self) -> usize {
        self.dims.len()
    }

    pub fn total_dim(&self) -> usize {
        self.amps.len()
    }

    pub fn norm(&self) -> f64 {
        self.amps.iter().map(|a| a.abs2()).sum::<f64>().sqrt()
    }

    pub fn normalize(&mut self) {
        let n = self.norm();
        assert!(n > 0.0);
        let inv = 1.0 / n;
        for a in &mut self.amps {
            *a = a.scale(inv);
        }
    }

    /// `⟨self|other⟩`.
    pub fn inner(&self, other: &DenseState) -> C64 {
        assert_eq!(self.dims, other.dims);
        let mut acc = C64::ZERO;
        for (a, b) in self.amps.iter().zip(other.amps.iter()) {
            acc += a.conj() * *b;
        }
        acc
    }

    /// `|⟨self|other⟩|²` (assumes both are normalized).
    pub fn fidelity(&self, other: &DenseState) -> f64 {
        self.inner(other).abs2()
    }

    /// Apply a single-site gate.
    pub fn apply1(&mut self, site: usize, g: &Mat) {
        let d = self.dims[site];
        assert_eq!(g.rows, d);
        assert_eq!(g.cols, d);
        let stride = strides(&self.dims)[site];
        let n = self.amps.len();
        let mut buf = vec![C64::ZERO; d];
        let mut idx = 0;
        while idx < n {
            if (idx / stride).is_multiple_of(d) {
                for (x, b) in buf.iter_mut().enumerate() {
                    *b = self.amps[idx + x * stride];
                }
                for (x, _) in buf.iter().enumerate() {
                    let mut acc = C64::ZERO;
                    for (y, b) in buf.iter().enumerate() {
                        acc += g.at(x, y) * *b;
                    }
                    self.amps[idx + x * stride] = acc;
                }
            }
            idx += 1;
        }
    }

    /// Apply a two-site gate on the ordered pair `(a, b)` (any positions,
    /// `a != b`). The gate is `(da·db)×(da·db)` with basis index
    /// `x_a·db + x_b`.
    pub fn apply2(&mut self, a: usize, b: usize, g: &Mat) {
        assert!(a != b);
        let (da, db) = (self.dims[a], self.dims[b]);
        assert_eq!(g.rows, da * db);
        let st = strides(&self.dims);
        let (sa, sb) = (st[a], st[b]);
        let n = self.amps.len();
        let block = da * db;
        let mut buf = vec![C64::ZERO; block];
        for idx in 0..n {
            let xa = (idx / sa) % da;
            let xb = (idx / sb) % db;
            if xa != 0 || xb != 0 {
                continue;
            }
            for x in 0..da {
                for y in 0..db {
                    buf[x * db + y] = self.amps[idx + x * sa + y * sb];
                }
            }
            for x in 0..da {
                for y in 0..db {
                    let mut acc = C64::ZERO;
                    let row = x * db + y;
                    for (col, v) in buf.iter().enumerate() {
                        acc += g.at(row, col) * *v;
                    }
                    self.amps[idx + x * sa + y * sb] = acc;
                }
            }
        }
    }

    /// Probability of measuring basis digit `x` on `site` (marginal).
    pub fn site_probability(&self, site: usize, x: usize) -> f64 {
        let d = self.dims[site];
        let stride = strides(&self.dims)[site];
        self.amps
            .iter()
            .enumerate()
            .filter(|(idx, _)| (idx / stride) % d == x)
            .map(|(_, a)| a.abs2())
            .sum()
    }

    /// Digits of a flattened basis index.
    pub fn digits_of(&self, mut idx: usize) -> Vec<usize> {
        let mut out = vec![0; self.dims.len()];
        for i in (0..self.dims.len()).rev() {
            out[i] = idx % self.dims[i];
            idx /= self.dims[i];
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gates;
    use crate::mat::Rng;

    /// Build the full operator for a 1-site gate by Kronecker products.
    fn full_op_1(dims: &[usize], site: usize, g: &Mat) -> Mat {
        let mut m = Mat::identity(1);
        for (i, &d) in dims.iter().enumerate() {
            let f = if i == site {
                g.clone()
            } else {
                Mat::identity(d)
            };
            m = m.kron(&f);
        }
        m
    }

    /// Build the full operator for a 2-site gate elementwise.
    fn full_op_2(dims: &[usize], a: usize, b: usize, g: &Mat) -> Mat {
        let n = total_dim(dims);
        let st = strides(dims);
        let (da, db) = (dims[a], dims[b]);
        let mut m = Mat::zeros(n, n);
        for col in 0..n {
            // digits of col
            let mut rem = col;
            let mut digits = vec![0usize; dims.len()];
            for i in (0..dims.len()).rev() {
                digits[i] = rem % dims[i];
                rem /= dims[i];
            }
            let (xa, xb) = (digits[a], digits[b]);
            let base = col - xa * st[a] - xb * st[b];
            for ya in 0..da {
                for yb in 0..db {
                    let row = base + ya * st[a] + yb * st[b];
                    m.add_at(row, col, g.at(ya * db + yb, xa * db + xb));
                }
            }
        }
        m
    }

    fn apply_full(m: &Mat, v: &[C64]) -> Vec<C64> {
        let mut out = vec![C64::ZERO; m.rows];
        for (r, o) in out.iter_mut().enumerate() {
            for (c, vc) in v.iter().enumerate() {
                *o += m.at(r, c) * *vc;
            }
        }
        out
    }

    fn random_state(dims: &[usize], rng: &mut Rng) -> DenseState {
        let mut s = DenseState::zero_state(dims);
        for a in &mut s.amps {
            *a = rng.c_gaussian();
        }
        s.normalize();
        s
    }

    fn max_diff(a: &[C64], b: &[C64]) -> f64 {
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (*x - *y).abs())
            .fold(0.0, f64::max)
    }

    #[test]
    fn apply1_matches_kron() {
        let dims = [2, 3, 4];
        let mut rng = Rng::new(21);
        for site in 0..3 {
            let g = gates::random(dims[site], &mut rng);
            let mut s = random_state(&dims, &mut rng);
            let expect = apply_full(&full_op_1(&dims, site, &g), &s.amps);
            s.apply1(site, &g);
            assert!(max_diff(&s.amps, &expect) < 1e-12, "site {}", site);
        }
    }

    #[test]
    fn apply2_matches_full_operator_all_pairs() {
        let dims = [2, 3, 4];
        let mut rng = Rng::new(22);
        for a in 0..3 {
            for b in 0..3 {
                if a == b {
                    continue;
                }
                let g = gates::random(dims[a] * dims[b], &mut rng);
                let mut s = random_state(&dims, &mut rng);
                let expect = apply_full(&full_op_2(&dims, a, b, &g), &s.amps);
                s.apply2(a, b, &g);
                assert!(max_diff(&s.amps, &expect) < 1e-11, "pair {} {}", a, b);
            }
        }
    }

    #[test]
    fn cshift_entangles_to_generalized_bell() {
        // fourier on site 0 then cshift(0→1) on a (3,3) pair yields
        // Σ_j |j,j⟩/√3.
        let mut s = DenseState::zero_state(&[3, 3]);
        s.apply1(0, &gates::fourier(3));
        s.apply2(0, 1, &gates::cshift(3, 3));
        let r = 1.0 / 3.0_f64.sqrt();
        for j in 0..3 {
            let idx = j * 3 + j;
            assert!((s.amps[idx].abs() - r).abs() < 1e-12);
        }
        assert!((s.norm() - 1.0).abs() < 1e-12);
        for j in 0..3 {
            assert!((s.site_probability(0, j) - 1.0 / 3.0).abs() < 1e-12);
        }
    }

    #[test]
    fn unitarity_preserves_norm() {
        let dims = [2, 3, 4, 3];
        let mut rng = Rng::new(23);
        let mut s = random_state(&dims, &mut rng);
        for _ in 0..10 {
            let a = rng.below(4);
            let mut b = rng.below(4);
            while b == a {
                b = rng.below(4);
            }
            s.apply2(a, b, &gates::random(dims[a] * dims[b], &mut rng));
        }
        assert!((s.norm() - 1.0).abs() < 1e-10);
    }
}
