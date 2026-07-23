//! Dense complex matrices, a deterministic RNG, Haar-random unitaries, and a
//! one-sided Jacobi singular value decomposition.
//!
//! The SVD is the workhorse of the whole library: every bond truncation,
//! canonicalization sweep, operator-Schmidt decomposition and site split goes
//! through [`svd`] / [`svd_trunc`]. One-sided Jacobi was chosen because it is
//! short, numerically robust, and easy to verify — matrix sizes in this
//! library are small (hundreds at most), so asymptotic speed is irrelevant.

use crate::c64::C64;

/// Row-major dense complex matrix.
#[derive(Clone, PartialEq)]
pub struct Mat {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<C64>,
}

impl Mat {
    pub fn zeros(rows: usize, cols: usize) -> Mat {
        Mat {
            rows,
            cols,
            data: vec![C64::ZERO; rows * cols],
        }
    }

    pub fn identity(n: usize) -> Mat {
        let mut m = Mat::zeros(n, n);
        for i in 0..n {
            m.set(i, i, C64::ONE);
        }
        m
    }

    pub fn from_fn(rows: usize, cols: usize, mut f: impl FnMut(usize, usize) -> C64) -> Mat {
        let mut m = Mat::zeros(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                m.set(r, c, f(r, c));
            }
        }
        m
    }

    #[inline]
    pub fn at(&self, r: usize, c: usize) -> C64 {
        debug_assert!(r < self.rows && c < self.cols);
        self.data[r * self.cols + c]
    }

    #[inline]
    pub fn set(&mut self, r: usize, c: usize, v: C64) {
        debug_assert!(r < self.rows && c < self.cols);
        self.data[r * self.cols + c] = v;
    }

    #[inline]
    pub fn add_at(&mut self, r: usize, c: usize, v: C64) {
        self.data[r * self.cols + c] += v;
    }

    /// Matrix product `self * other`.
    pub fn mul(&self, other: &Mat) -> Mat {
        assert_eq!(
            self.cols, other.rows,
            "matmul dimension mismatch: {}x{} * {}x{}",
            self.rows, self.cols, other.rows, other.cols
        );
        let mut out = Mat::zeros(self.rows, other.cols);
        for i in 0..self.rows {
            for k in 0..self.cols {
                let aik = self.at(i, k);
                if aik.re == 0.0 && aik.im == 0.0 {
                    continue;
                }
                let orow = i * out.cols;
                let brow = k * other.cols;
                for j in 0..other.cols {
                    out.data[orow + j] += aik * other.data[brow + j];
                }
            }
        }
        out
    }

    /// Conjugate transpose.
    pub fn adjoint(&self) -> Mat {
        let mut out = Mat::zeros(self.cols, self.rows);
        for r in 0..self.rows {
            for c in 0..self.cols {
                out.set(c, r, self.at(r, c).conj());
            }
        }
        out
    }

    /// Kronecker product `self ⊗ other` (self is the most-significant factor).
    pub fn kron(&self, other: &Mat) -> Mat {
        let mut out = Mat::zeros(self.rows * other.rows, self.cols * other.cols);
        for a in 0..self.rows {
            for b in 0..self.cols {
                let sab = self.at(a, b);
                for c in 0..other.rows {
                    for d in 0..other.cols {
                        out.set(a * other.rows + c, b * other.cols + d, sab * other.at(c, d));
                    }
                }
            }
        }
        out
    }

    pub fn scale(&self, s: C64) -> Mat {
        Mat {
            rows: self.rows,
            cols: self.cols,
            data: self.data.iter().map(|z| *z * s).collect(),
        }
    }

    pub fn sub(&self, other: &Mat) -> Mat {
        assert!(self.rows == other.rows && self.cols == other.cols);
        Mat {
            rows: self.rows,
            cols: self.cols,
            data: self
                .data
                .iter()
                .zip(other.data.iter())
                .map(|(a, b)| *a - *b)
                .collect(),
        }
    }

    pub fn frobenius(&self) -> f64 {
        self.data.iter().map(|z| z.abs2()).sum::<f64>().sqrt()
    }

    pub fn max_abs(&self) -> f64 {
        self.data.iter().map(|z| z.abs()).fold(0.0, f64::max)
    }

    /// Largest entry-wise deviation from another matrix.
    pub fn max_abs_diff(&self, other: &Mat) -> f64 {
        self.sub(other).max_abs()
    }

    /// True when `self† self == I` within `tol` (columns orthonormal).
    pub fn is_isometry(&self, tol: f64) -> bool {
        let g = self.adjoint().mul(self);
        g.max_abs_diff(&Mat::identity(self.cols)) < tol
    }

    /// True when the matrix is square and unitary within `tol`.
    pub fn is_unitary(&self, tol: f64) -> bool {
        self.rows == self.cols && self.is_isometry(tol)
    }

    /// Matrix power by repeated multiplication (small exponents only).
    pub fn pow(&self, mut n: usize) -> Mat {
        assert_eq!(self.rows, self.cols);
        let mut acc = Mat::identity(self.rows);
        let mut base = self.clone();
        while n > 0 {
            if n & 1 == 1 {
                acc = acc.mul(&base);
            }
            base = base.mul(&base);
            n >>= 1;
        }
        acc
    }
}

impl std::fmt::Debug for Mat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Mat {}x{} [", self.rows, self.cols)?;
        for r in 0..self.rows {
            write!(f, "  ")?;
            for c in 0..self.cols {
                write!(f, "{} ", self.at(r, c))?;
            }
            writeln!(f)?;
        }
        write!(f, "]")
    }
}

// ---------------------------------------------------------------------------
// Deterministic RNG (SplitMix64) — reproducible experiments, no `rand` crate.
// ---------------------------------------------------------------------------

/// SplitMix64 pseudo-random generator with Box–Muller Gaussian sampling.
pub struct Rng {
    state: u64,
    spare_gauss: Option<f64>,
}

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng {
            state: seed,
            spare_gauss: None,
        }
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }

    /// Uniform in `[0, 1)`.
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform integer in `[0, n)`.
    pub fn below(&mut self, n: usize) -> usize {
        (self.uniform() * n as f64) as usize % n
    }

    /// Standard normal sample (Box–Muller).
    pub fn gaussian(&mut self) -> f64 {
        if let Some(g) = self.spare_gauss.take() {
            return g;
        }
        // Avoid log(0).
        let u1 = loop {
            let u = self.uniform();
            if u > 1e-300 {
                break u;
            }
        };
        let u2 = self.uniform();
        let r = (-2.0 * u1.ln()).sqrt();
        let th = 2.0 * std::f64::consts::PI * u2;
        self.spare_gauss = Some(r * th.sin());
        r * th.cos()
    }

    /// Standard complex Gaussian entry.
    pub fn c_gaussian(&mut self) -> C64 {
        C64::new(self.gaussian(), self.gaussian())
    }
}

/// Haar-distributed random unitary via Ginibre + (modified Gram-Schmidt) QR
/// with the Mezzadri phase fix.
pub fn random_unitary(n: usize, rng: &mut Rng) -> Mat {
    let g = Mat::from_fn(n, n, |_, _| rng.c_gaussian());
    // Columns of g, orthonormalized.
    let mut q = g;
    for j in 0..n {
        // Two rounds of MGS projection for numerical orthogonality.
        for _round in 0..2 {
            for i in 0..j {
                let mut dot = C64::ZERO;
                for r in 0..n {
                    dot += q.at(r, i).conj() * q.at(r, j);
                }
                for r in 0..n {
                    let v = q.at(r, j) - q.at(r, i) * dot;
                    q.set(r, j, v);
                }
            }
        }
        let mut nrm2 = 0.0;
        for r in 0..n {
            nrm2 += q.at(r, j).abs2();
        }
        let nrm = nrm2.sqrt();
        assert!(nrm > 1e-100, "degenerate Ginibre sample");
        // Phase of the R diagonal entry (before normalization the top
        // component is arbitrary; use the phase of the pivot to fix Haar).
        let phase_fix = q.at(j, j).phase().conj();
        for r in 0..n {
            let v = q.at(r, j).scale(1.0 / nrm) * phase_fix;
            q.set(r, j, v);
        }
    }
    q
}

// ---------------------------------------------------------------------------
// One-sided Jacobi SVD for complex matrices.
// ---------------------------------------------------------------------------

/// Result of a singular value decomposition `a = u * diag(s) * vh` with
/// `u: m×k`, `s: k` (descending, non-negative), `vh: k×n`, `k = min(m, n)`.
pub struct Svd {
    pub u: Mat,
    pub s: Vec<f64>,
    pub vh: Mat,
}

/// Full (economy) SVD via one-sided Jacobi rotations.
pub fn svd(a: &Mat) -> Svd {
    if a.rows >= a.cols {
        jacobi_svd_tall(a)
    } else {
        // A = (A†)† ; A† is tall. If A† = U' Σ V'† then A = V' Σ U'†.
        let inner = jacobi_svd_tall(&a.adjoint());
        Svd {
            u: inner.vh.adjoint(),
            s: inner.s,
            vh: inner.u.adjoint(),
        }
    }
}

/// One-sided Jacobi for `m >= n`: rotates column pairs of a working copy `w`
/// until all columns are mutually orthogonal; then `σ_j = ‖w_j‖`,
/// `u_j = w_j / σ_j`, and the accumulated rotations form `V`.
fn jacobi_svd_tall(a: &Mat) -> Svd {
    let (m, n) = (a.rows, a.cols);
    debug_assert!(m >= n);
    let mut w = a.clone();
    let mut v = Mat::identity(n);
    let tol = 1e-14;
    let max_sweeps = 100;
    // Columns whose squared norm falls below this are numerically zero
    // (their singular value is beneath roundoff resolution): freeze them so
    // rotations never grind them into the denormal range, where Gram
    // entries underflow and rotation parameters turn to garbage.
    let fnorm2: f64 = a.data.iter().map(|z| z.abs2()).sum();
    let zero2 = fnorm2 * 1e-32;

    for _sweep in 0..max_sweeps {
        let mut changed = false;
        for p in 0..n.saturating_sub(1) {
            for q in (p + 1)..n {
                // 2x2 Gram block of columns p, q.
                let mut alpha = 0.0;
                let mut beta = 0.0;
                let mut gamma = C64::ZERO;
                for i in 0..m {
                    let wp = w.at(i, p);
                    let wq = w.at(i, q);
                    alpha += wp.abs2();
                    beta += wq.abs2();
                    gamma += wp.conj() * wq;
                }
                if alpha <= zero2 || beta <= zero2 {
                    continue;
                }
                let g_abs = gamma.abs();
                // Separate square roots: alpha*beta can underflow.
                if g_abs <= tol * alpha.sqrt() * beta.sqrt() || g_abs == 0.0 {
                    continue;
                }
                changed = true;
                // Absorb the phase of gamma into column q so the 2x2 problem
                // becomes real, then apply the classic Jacobi rotation.
                let ph = gamma.phase().conj(); // e^{-i arg(gamma)}
                let tau = (beta - alpha) / (2.0 * g_abs);
                let t = if tau >= 0.0 {
                    1.0 / (tau + (1.0 + tau * tau).sqrt())
                } else {
                    -1.0 / (-tau + (1.0 + tau * tau).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                for i in 0..m {
                    let wp = w.at(i, p);
                    let wq = w.at(i, q) * ph;
                    w.set(i, p, wp.scale(c) - wq.scale(s));
                    w.set(i, q, wp.scale(s) + wq.scale(c));
                }
                for i in 0..n {
                    let vp = v.at(i, p);
                    let vq = v.at(i, q) * ph;
                    v.set(i, p, vp.scale(c) - vq.scale(s));
                    v.set(i, q, vp.scale(s) + vq.scale(c));
                }
            }
        }
        if !changed {
            break;
        }
    }

    // Column norms are the singular values.
    let mut order: Vec<usize> = (0..n).collect();
    let norms: Vec<f64> = (0..n)
        .map(|j| (0..m).map(|i| w.at(i, j).abs2()).sum::<f64>().sqrt())
        .collect();
    order.sort_by(|&x, &y| norms[y].partial_cmp(&norms[x]).unwrap());

    let mut u = Mat::zeros(m, n);
    let mut vh = Mat::zeros(n, n);
    let mut s = vec![0.0; n];
    let smax = norms.iter().cloned().fold(0.0, f64::max);
    for (out_j, &j) in order.iter().enumerate() {
        // A column norm at or below roundoff resolution is a numerical zero:
        // report σ = 0 and leave the U column empty rather than normalizing
        // rotation debris into a fake direction. (`svd_trunc` drops these.)
        if norms[j] > smax * 1e-15 {
            s[out_j] = norms[j];
            let inv = 1.0 / norms[j];
            for i in 0..m {
                u.set(i, out_j, w.at(i, j).scale(inv));
            }
        }
        // vh row out_j = column j of v, conjugated.
        for i in 0..n {
            vh.set(out_j, i, v.at(i, j).conj());
        }
    }
    Svd { u, s, vh }
}

/// Truncation policy for SVD-based bond compression.
#[derive(Clone, Copy, Debug)]
pub struct TruncSpec {
    /// Hard cap on the number of retained singular values.
    pub max_rank: usize,
    /// Relative discarded-weight tolerance: the smallest singular values are
    /// dropped while their cumulative squared weight stays below
    /// `cutoff * Σ s_i^2`.
    pub cutoff: f64,
}

impl TruncSpec {
    pub fn exact() -> TruncSpec {
        TruncSpec {
            max_rank: usize::MAX,
            cutoff: 0.0,
        }
    }

    pub fn new(max_rank: usize, cutoff: f64) -> TruncSpec {
        TruncSpec { max_rank, cutoff }
    }
}

/// Result of a truncated SVD: `a ≈ u * diag(s) * vh` plus the relative
/// squared weight that was discarded.
pub struct TruncatedSvd {
    pub u: Mat,
    pub s: Vec<f64>,
    pub vh: Mat,
    /// `Σ_dropped s_i^2 / Σ_all s_i^2` — 0 for exact splits.
    pub discarded_weight: f64,
}

/// Orthonormalize the columns of `a` by modified Gram–Schmidt (two
/// projection rounds), dropping columns that collapse to numerical zero.
fn orth_columns(a: &Mat) -> Mat {
    let (m, n) = (a.rows, a.cols);
    let scale = (0..n)
        .map(|j| (0..m).map(|i| a.at(i, j).abs2()).sum::<f64>().sqrt())
        .fold(0.0, f64::max)
        .max(1e-300);
    let mut cols: Vec<Vec<C64>> = Vec::with_capacity(n);
    for j in 0..n {
        let mut v: Vec<C64> = (0..m).map(|i| a.at(i, j)).collect();
        for _round in 0..2 {
            for q in &cols {
                let mut dot = C64::ZERO;
                for i in 0..m {
                    dot += q[i].conj() * v[i];
                }
                for i in 0..m {
                    v[i] -= q[i] * dot;
                }
            }
        }
        let nrm = v.iter().map(|z| z.abs2()).sum::<f64>().sqrt();
        if nrm > 1e-13 * scale {
            let inv = 1.0 / nrm;
            for z in &mut v {
                *z = z.scale(inv);
            }
            cols.push(v);
        }
    }
    let k = cols.len().max(1);
    Mat::from_fn(m, k, |i, j| {
        if j < cols.len() {
            cols[j][i]
        } else {
            C64::ZERO
        }
    })
}

/// Randomized range-finder SVD (Halko–Martinsson–Tropp with one power
/// iteration): returns a rank-≤`k` factorization capturing the dominant
/// subspace. The tail beyond the sampled subspace is *not* represented —
/// [`svd_trunc`] accounts for it against the true Frobenius weight.
fn randomized_svd(a: &Mat, k: usize, rng_seed: u64) -> Svd {
    let (m, n) = (a.rows, a.cols);
    let k = k.min(m).min(n);
    let mut rng = Rng::new(rng_seed);
    let omega = Mat::from_fn(n, k, |_, _| rng.c_gaussian());
    let q0 = orth_columns(&a.mul(&omega));
    // One power iteration sharpens the captured subspace.
    let q1 = orth_columns(&a.adjoint().mul(&q0));
    let q = orth_columns(&a.mul(&q1));
    let b = q.adjoint().mul(a); // k' × n, k' small
    let inner = svd(&b);
    Svd {
        u: q.mul(&inner.u),
        s: inner.s,
        vh: inner.vh,
    }
}

const RAND_OVERSAMPLE: usize = 8;

/// SVD followed by rank truncation. Exact zeros (relative to `s[0]`) are
/// always dropped; at least one singular value is always kept.
///
/// When the rank cap is far below the matrix's smaller side, a randomized
/// range-finder is used instead of the full Jacobi decomposition; any
/// weight outside the sampled subspace is charged to `discarded_weight`
/// (measured against the true Frobenius norm), so the accounting stays
/// honest either way.
pub fn svd_trunc(a: &Mat, spec: TruncSpec) -> TruncatedSvd {
    let min_side = a.rows.min(a.cols);
    let frob2: f64 = a.data.iter().map(|z| z.abs2()).sum();
    let use_randomized = spec.max_rank != usize::MAX
        && min_side >= 128
        && spec.max_rank + RAND_OVERSAMPLE <= min_side / 2;
    let full = if use_randomized {
        let seed = 0x9E3779B97F4A7C15 ^ ((a.rows as u64) << 32) ^ a.cols as u64;
        randomized_svd(a, spec.max_rank + RAND_OVERSAMPLE, seed)
    } else {
        svd(a)
    };
    let total: f64 = frob2;
    if total == 0.0 {
        // Zero matrix: keep a single null direction to preserve shapes.
        return TruncatedSvd {
            u: submat_cols(&full.u, 1),
            s: vec![0.0],
            vh: submat_rows(&full.vh, 1),
            discarded_weight: 0.0,
        };
    }
    let smax = full.s[0];
    let k_all = full.s.len();
    // Largest k satisfying the cutoff, scanning tail weight from the back.
    let mut keep = k_all.min(spec.max_rank);
    // Drop numerically-zero values regardless of cutoff.
    while keep > 1 && full.s[keep - 1] <= smax * 1e-14 {
        keep -= 1;
    }
    if spec.cutoff > 0.0 {
        let budget = spec.cutoff * total;
        // Weight already outside the returned factors (randomized-path
        // projection loss) counts against the budget from the start.
        let represented: f64 = full.s.iter().map(|x| x * x).sum();
        let mut tail: f64 =
            (total - represented).max(0.0) + full.s[keep..].iter().map(|x| x * x).sum::<f64>();
        while keep > 1 {
            let w = full.s[keep - 1] * full.s[keep - 1];
            if tail + w <= budget {
                tail += w;
                keep -= 1;
            } else {
                break;
            }
        }
    }
    // Discarded weight = values actually dropped, plus (randomized path
    // only) the weight left outside the sampled subspace. Computed from the
    // dropped values directly so exact splits report an exact zero.
    let dropped: f64 = full.s[keep..].iter().map(|x| x * x).sum();
    let projection_loss = if use_randomized {
        let represented: f64 = full.s.iter().map(|x| x * x).sum();
        (total - represented).max(0.0)
    } else {
        0.0
    };
    TruncatedSvd {
        u: submat_cols(&full.u, keep),
        s: full.s[..keep].to_vec(),
        vh: submat_rows(&full.vh, keep),
        discarded_weight: (dropped + projection_loss) / total,
    }
}

fn submat_cols(a: &Mat, k: usize) -> Mat {
    Mat::from_fn(a.rows, k, |r, c| a.at(r, c))
}

fn submat_rows(a: &Mat, k: usize) -> Mat {
    Mat::from_fn(k, a.cols, |r, c| a.at(r, c))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reconstruct(d: &Svd) -> Mat {
        let k = d.s.len();
        let mut us = d.u.clone();
        for j in 0..k {
            for i in 0..us.rows {
                let v = us.at(i, j).scale(d.s[j]);
                us.set(i, j, v);
            }
        }
        us.mul(&d.vh)
    }

    #[test]
    fn svd_reconstructs_random_shapes() {
        let mut rng = Rng::new(7);
        for &(m, n) in &[
            (1, 1),
            (2, 2),
            (3, 5),
            (5, 3),
            (7, 7),
            (20, 6),
            (6, 20),
            (30, 30),
        ] {
            let a = Mat::from_fn(m, n, |_, _| rng.c_gaussian());
            let d = svd(&a);
            let err = reconstruct(&d).max_abs_diff(&a);
            assert!(err < 1e-10, "shape {}x{} err {}", m, n, err);
            // Descending, non-negative.
            for w in d.s.windows(2) {
                assert!(w[0] >= w[1] - 1e-12);
            }
            assert!(d.s.iter().all(|&x| x >= 0.0));
            assert!(d.u.is_isometry(1e-10), "U not isometric {}x{}", m, n);
            assert!(d.vh.adjoint().is_isometry(1e-10), "V not isometric");
        }
    }

    #[test]
    fn svd_handles_rank_deficiency() {
        let mut rng = Rng::new(11);
        // Rank-2 6x5 matrix.
        let b = Mat::from_fn(6, 2, |_, _| rng.c_gaussian());
        let c = Mat::from_fn(2, 5, |_, _| rng.c_gaussian());
        let a = b.mul(&c);
        let d = svd(&a);
        assert!(reconstruct(&d).max_abs_diff(&a) < 1e-10);
        assert!(d.s[2] < 1e-10 * d.s[0]);
        let t = svd_trunc(&a, TruncSpec::exact());
        assert_eq!(t.s.len(), 2, "exact truncation should drop null space");
        assert!(t.discarded_weight < 1e-20);
    }

    #[test]
    fn svd_trunc_respects_max_rank_and_reports_weight() {
        let mut rng = Rng::new(13);
        let a = Mat::from_fn(8, 8, |_, _| rng.c_gaussian());
        let t = svd_trunc(&a, TruncSpec::new(3, 0.0));
        assert_eq!(t.s.len(), 3);
        let full = svd(&a);
        let total: f64 = full.s.iter().map(|x| x * x).sum();
        let dropped: f64 = full.s[3..].iter().map(|x| x * x).sum();
        assert!((t.discarded_weight - dropped / total).abs() < 1e-12);
    }

    #[test]
    fn randomized_path_matches_full_svd() {
        // A 140×200 matrix with a rank-18 dominant part and a tiny tail:
        // truncation at max_rank 30 must go through the randomized path and
        // agree with the full Jacobi decomposition.
        let mut rng = Rng::new(314);
        let b = Mat::from_fn(140, 18, |_, _| rng.c_gaussian());
        let c = Mat::from_fn(18, 200, |_, _| rng.c_gaussian());
        let mut a = b.mul(&c);
        for z in &mut a.data {
            *z += rng.c_gaussian().scale(1e-10);
        }
        let spec = TruncSpec::new(30, 1e-12);
        assert!(
            spec.max_rank + super::RAND_OVERSAMPLE <= 140 / 2,
            "test must exercise the randomized path"
        );
        let t = svd_trunc(&a, spec);
        let full = svd(&a);
        assert!(t.s.len() >= 18);
        for j in 0..18 {
            assert!(
                (t.s[j] - full.s[j]).abs() < 1e-8 * full.s[0],
                "s[{}]: {} vs {}",
                j,
                t.s[j],
                full.s[j]
            );
        }
        // Reconstruction error at the noise floor.
        let mut us = t.u.clone();
        for j in 0..t.s.len() {
            for i in 0..us.rows {
                let v = us.at(i, j).scale(t.s[j]);
                us.set(i, j, v);
            }
        }
        let err = us.mul(&t.vh).max_abs_diff(&a);
        assert!(err < 1e-7, "reconstruction err {}", err);
        assert!(t.u.is_isometry(1e-9));
    }

    #[test]
    fn svd_survives_null_space_grinding() {
        // Regression: matrices whose null-space columns get rotated down into
        // the denormal range used to corrupt the significant subspace
        // (reconstruction error ~1e-6). Structured low-rank inputs with
        // wide dynamic range must reconstruct to machine precision.
        let mut rng = Rng::new(2718);
        for &(m, n, rank) in &[(15, 9, 3), (9, 15, 3), (12, 12, 2), (45, 9, 4)] {
            let b = Mat::from_fn(m, rank, |_, _| rng.c_gaussian());
            let mut c = Mat::from_fn(rank, n, |_, _| rng.c_gaussian());
            // Widely spread scales across the rank-one components.
            for (r, scale) in (0..rank).zip([1.0, 1e-4, 1e-8, 1e-12]) {
                for j in 0..n {
                    let v = c.at(r, j).scale(scale);
                    c.set(r, j, v);
                }
            }
            let a = b.mul(&c);
            let d = svd(&a);
            let err = reconstruct(&d).max_abs_diff(&a);
            let smax = d.s[0];
            assert!(
                err < 1e-12 * smax,
                "{}x{} rank {}: err {} smax {}",
                m,
                n,
                rank,
                err,
                smax
            );
            assert!(d.u.is_isometry(1e-10) || d.s.contains(&0.0));
        }
    }

    #[test]
    fn random_unitary_is_unitary() {
        let mut rng = Rng::new(42);
        for n in [1, 2, 3, 5, 9] {
            let u = random_unitary(n, &mut rng);
            assert!(u.is_unitary(1e-10), "n={}", n);
        }
    }

    #[test]
    fn kron_and_pow() {
        let x = Mat::from_fn(2, 2, |r, c| if r != c { C64::ONE } else { C64::ZERO });
        let k = x.kron(&Mat::identity(3));
        assert_eq!(k.rows, 6);
        // X^2 = I
        assert!(x.pow(2).max_abs_diff(&Mat::identity(2)) < 1e-14);
    }
}
