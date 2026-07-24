//! Heterogeneous-dimension **matrix product state** (MPS) simulator.
//!
//! This is the engine that makes zigzag qudit chains *classically
//! computable*: a state over sites with local dimensions `d_0, …, d_{n-1}` is
//! stored as a chain of 3-index tensors `T_i[l, p, r]` (left bond, physical,
//! right bond) whose total size is linear in `n` for bounded bond dimension,
//! instead of the product `Π d_i` a dense vector needs.
//!
//! Unusual features relative to a garden-variety MPS code:
//!
//! * **Per-site local dimensions** — nothing assumes a uniform `d`, so
//!   dimension waves like `2,3,4,5,4,3,2` are first-class.
//! * **Long-range two-site gates** between *any* pair of sites (the
//!   cross-scale couplings of the zigzag structure), implemented by an
//!   operator-Schmidt decomposition of the gate into an MPO threaded through
//!   the intervening sites, followed by a two-pass compression sweep.
//! * **Scale morphing** — [`Mps::merge_sites`] fuses two neighbouring sites
//!   into one coarser qudit, [`Mps::split_site`] factors a coarse qudit into
//!   two finer ones, and [`Mps::promote_site`] / [`Mps::demote_site`] grow or
//!   shrink a local dimension in place. The *representation* expands and
//!   collapses while the state is preserved (exactly, up to reported
//!   truncation), so a single gate at a coarse scale acts as a whole
//!   subcircuit at the fine scale.
//!
//! ## Canonical-form invariant
//!
//! Tensors strictly left of `center` are left-canonical isometries, tensors
//! strictly right of it are right-canonical, and the norm of the state lives
//! in `tensors[center]`. Every public operation restores this invariant, and
//! truncations are only ever performed with a canonical environment (which is
//! what makes each local SVD truncation globally optimal).

use crate::c64::C64;
use crate::dense::DenseState;
use crate::mat::{svd, svd_trunc, Mat, TruncSpec};

/// One MPS site tensor `T[l, p, r]`, stored row-major as `(l·d + p)·dr + r`.
#[derive(Clone)]
pub struct SiteTensor {
    pub dl: usize,
    pub d: usize,
    pub dr: usize,
    pub data: Vec<C64>,
}

impl SiteTensor {
    pub fn zeros(dl: usize, d: usize, dr: usize) -> SiteTensor {
        SiteTensor {
            dl,
            d,
            dr,
            data: vec![C64::ZERO; dl * d * dr],
        }
    }

    #[inline]
    pub fn at(&self, l: usize, p: usize, r: usize) -> C64 {
        self.data[(l * self.d + p) * self.dr + r]
    }

    #[inline]
    pub fn set(&mut self, l: usize, p: usize, r: usize, v: C64) {
        self.data[(l * self.d + p) * self.dr + r] = v;
    }

    pub fn frobenius(&self) -> f64 {
        self.data.iter().map(|z| z.abs2()).sum::<f64>().sqrt()
    }

    /// Matrix slice at fixed physical index: `(dl × dr)`.
    fn phys_slice(&self, p: usize) -> Mat {
        Mat::from_fn(self.dl, self.dr, |l, r| self.at(l, p, r))
    }
}

/// Matrix product state over a heterogeneous qudit chain.
#[derive(Clone)]
pub struct Mps {
    pub dims: Vec<usize>,
    pub tensors: Vec<SiteTensor>,
    /// Orthogonality center (see module docs).
    pub center: usize,
    /// Truncation policy applied by entangling operations.
    pub trunc: TruncSpec,
    /// Accumulated relative squared Schmidt weight discarded by truncations —
    /// 0 for exact simulation. (A tally of per-truncation losses, not an
    /// exact global error; measure fidelity against a reference when
    /// exactness matters.)
    pub discarded_weight: f64,
}

impl Mps {
    /// Computational basis product state `|digits⟩`.
    pub fn basis_state(dims: &[usize], digits: &[usize], trunc: TruncSpec) -> Mps {
        assert_eq!(dims.len(), digits.len());
        let tensors = dims
            .iter()
            .zip(digits)
            .map(|(&d, &x)| {
                assert!(x < d);
                let mut t = SiteTensor::zeros(1, d, 1);
                t.set(0, x, 0, C64::ONE);
                t
            })
            .collect();
        Mps {
            dims: dims.to_vec(),
            tensors,
            center: 0,
            trunc,
            discarded_weight: 0.0,
        }
    }

    /// All-zeros basis state.
    pub fn zero_state(dims: &[usize], trunc: TruncSpec) -> Mps {
        Mps::basis_state(dims, &vec![0; dims.len()], trunc)
    }

    /// Build an MPS from a dense state by successive SVD splits
    /// (left-to-right), truncating with `trunc`.
    pub fn from_dense(state: &DenseState, trunc: TruncSpec) -> Mps {
        let dims = state.dims.clone();
        let n = dims.len();
        assert!(n >= 1);
        let mut tensors: Vec<SiteTensor> = Vec::with_capacity(n);
        let mut discarded = 0.0;
        // carry: (bond × remaining) matrix.
        let mut carry = Mat {
            rows: 1,
            cols: state.amps.len(),
            data: state.amps.clone(),
        };
        for (i, &d) in dims.iter().enumerate().take(n - 1) {
            let bond = carry.rows;
            let rest = carry.cols / d;
            // reshape to (bond·d) × rest
            let m = Mat::from_fn(bond * d, rest, |row, col| {
                let (l, p) = (row / d, row % d);
                carry.at(l, p * rest + col)
            });
            let ts = svd_trunc(&m, trunc);
            discarded += ts.discarded_weight;
            let k = ts.s.len();
            let mut t = SiteTensor::zeros(bond, d, k);
            for l in 0..bond {
                for p in 0..d {
                    for r in 0..k {
                        t.set(l, p, r, ts.u.at(l * d + p, r));
                    }
                }
            }
            tensors.push(t);
            let renorm = if ts.discarded_weight > 0.0 {
                1.0 / (1.0 - ts.discarded_weight).sqrt()
            } else {
                1.0
            };
            carry = Mat::from_fn(k, rest, |r, c| ts.vh.at(r, c).scale(ts.s[r] * renorm));
            let _ = i;
        }
        // Last site: carry is (bond × d_last).
        let d_last = dims[n - 1];
        assert_eq!(carry.cols, d_last);
        let mut t = SiteTensor::zeros(carry.rows, d_last, 1);
        for l in 0..carry.rows {
            for p in 0..d_last {
                t.set(l, p, 0, carry.at(l, p));
            }
        }
        tensors.push(t);
        Mps {
            dims,
            tensors,
            center: n - 1,
            trunc,
            discarded_weight: discarded,
        }
    }

    pub fn num_sites(&self) -> usize {
        self.dims.len()
    }

    /// Bond dimensions between consecutive sites (`n-1` entries).
    pub fn bond_dims(&self) -> Vec<usize> {
        self.tensors[..self.tensors.len() - 1]
            .iter()
            .map(|t| t.dr)
            .collect()
    }

    pub fn max_bond_dim(&self) -> usize {
        self.bond_dims().into_iter().max().unwrap_or(1)
    }

    /// Total number of stored complex parameters.
    pub fn param_count(&self) -> usize {
        self.tensors.iter().map(|t| t.data.len()).sum()
    }

    // -- canonical-form machinery -----------------------------------------

    /// SVD `tensors[m]` as `[(l·d) × r]`, keep the isometry at `m`, absorb
    /// `S·V†` into `m+1`. Renormalizes to preserve the global norm when
    /// truncation discards weight.
    fn push_right(&mut self, m: usize, spec: TruncSpec) {
        let t = &self.tensors[m];
        let (dl, d, dr) = (t.dl, t.d, t.dr);
        let mm = Mat::from_fn(dl * d, dr, |row, r| t.at(row / d, row % d, r));
        let ts = svd_trunc(&mm, spec);
        let k = ts.s.len();
        let renorm = self.absorb_discard(ts.discarded_weight);
        let mut u = SiteTensor::zeros(dl, d, k);
        for l in 0..dl {
            for p in 0..d {
                for r in 0..k {
                    u.set(l, p, r, ts.u.at(l * d + p, r));
                }
            }
        }
        self.tensors[m] = u;
        let carry = Mat::from_fn(k, dr, |a, b| ts.vh.at(a, b).scale(ts.s[a] * renorm));
        let t2 = &self.tensors[m + 1];
        let (d2, dr2) = (t2.d, t2.dr);
        let mut new2 = SiteTensor::zeros(k, d2, dr2);
        for a in 0..k {
            for b in 0..dr {
                let cab = carry.at(a, b);
                if cab.re == 0.0 && cab.im == 0.0 {
                    continue;
                }
                for p in 0..d2 {
                    for r in 0..dr2 {
                        let v = new2.at(a, p, r) + cab * t2.at(b, p, r);
                        new2.set(a, p, r, v);
                    }
                }
            }
        }
        self.tensors[m + 1] = new2;
    }

    /// Mirror of [`Mps::push_right`]: SVD `tensors[m]` as `[l × (d·r)]`, keep
    /// the isometry at `m`, absorb `U·S` into `m-1`.
    fn push_left(&mut self, m: usize, spec: TruncSpec) {
        let t = &self.tensors[m];
        let (dl, d, dr) = (t.dl, t.d, t.dr);
        let mm = Mat::from_fn(dl, d * dr, |l, col| t.at(l, col / dr, col % dr));
        let ts = svd_trunc(&mm, spec);
        let k = ts.s.len();
        let renorm = self.absorb_discard(ts.discarded_weight);
        let mut v = SiteTensor::zeros(k, d, dr);
        for a in 0..k {
            for p in 0..d {
                for r in 0..dr {
                    v.set(a, p, r, ts.vh.at(a, p * dr + r));
                }
            }
        }
        self.tensors[m] = v;
        let carry = Mat::from_fn(dl, k, |l, a| ts.u.at(l, a).scale(ts.s[a] * renorm));
        let t0 = &self.tensors[m - 1];
        let (dl0, d0) = (t0.dl, t0.d);
        let mut new0 = SiteTensor::zeros(dl0, d0, k);
        for l in 0..dl0 {
            for p in 0..d0 {
                for b in 0..k {
                    let mut acc = C64::ZERO;
                    for a in 0..dl {
                        acc += t0.at(l, p, a) * carry.at(a, b);
                    }
                    new0.set(l, p, b, acc);
                }
            }
        }
        self.tensors[m - 1] = new0;
    }

    fn absorb_discard(&mut self, dw: f64) -> f64 {
        if dw > 0.0 {
            self.discarded_weight += dw;
            1.0 / (1.0 - dw).sqrt()
        } else {
            1.0
        }
    }

    /// Rebuild the canonical form from scratch and truncate every bond with
    /// the state's policy: an exact left→right sweep followed by a
    /// truncating right→left sweep. Valid regardless of the current
    /// canonical state (used after manual tensor surgery, e.g. MPO
    /// application); leaves the center at site 0.
    pub fn recompress(&mut self) {
        let n = self.num_sites();
        if n == 1 {
            self.center = 0;
            return;
        }
        // Truncating against a non-canonical environment is unsafe: raw
        // (e.g. freshly-contracted MPO-product) bond coordinates distort the
        // Schmidt weights, so a rank cap there can slice through genuinely
        // needed directions. Every lossy step below therefore happens with a
        // canonical environment; the first two passes are exact
        // (rank-revealing only).
        //
        // Cost control is by sweep choreography ("meet in the middle"):
        // pass 0 shrinks bonds inward from the right edge while the right
        // side is the smaller side of the matricization, so pass 1's
        // left-to-right exact sweep always sees its right bond already
        // reduced to true-rank scale — no SVD ever runs at raw-product
        // width on both sides.
        for m in (1..n).rev() {
            let t = &self.tensors[m];
            if t.d * t.dr > t.dl {
                break;
            }
            self.push_left(m, TruncSpec::exact());
        }
        for m in 0..n - 1 {
            self.push_right(m, TruncSpec::exact());
        }
        self.center = n - 1;
        let spec = self.trunc;
        for m in (1..n).rev() {
            self.push_left(m, spec);
        }
        self.center = 0;
    }

    /// Move the orthogonality center to `target` (exact — only numerically
    /// null Schmidt directions are dropped).
    pub fn move_center_to(&mut self, target: usize) {
        assert!(target < self.num_sites());
        while self.center < target {
            self.push_right(self.center, TruncSpec::exact());
            self.center += 1;
        }
        while self.center > target {
            self.push_left(self.center, TruncSpec::exact());
            self.center -= 1;
        }
    }

    // -- gate application --------------------------------------------------

    /// Apply a single-site gate.
    pub fn apply1(&mut self, site: usize, g: &Mat) {
        let d = self.dims[site];
        assert_eq!(g.rows, d);
        assert_eq!(g.cols, d);
        self.move_center_to(site);
        let t = &self.tensors[site];
        let mut new = SiteTensor::zeros(t.dl, d, t.dr);
        for l in 0..t.dl {
            for p in 0..d {
                for r in 0..t.dr {
                    let mut acc = C64::ZERO;
                    for q in 0..d {
                        acc += g.at(p, q) * t.at(l, q, r);
                    }
                    new.set(l, p, r, acc);
                }
            }
        }
        self.tensors[site] = new;
    }

    /// Apply a two-site gate on the ordered pair `(a, b)` — any two distinct
    /// sites, adjacent or not. The gate is `(d_a·d_b)×(d_a·d_b)` over basis
    /// index `x_a·d_b + x_b` (matching [`crate::dense::DenseState::apply2`]).
    pub fn apply2(&mut self, a: usize, b: usize, g: &Mat) {
        assert!(a != b);
        let (da, db) = (self.dims[a], self.dims[b]);
        assert_eq!(g.rows, da * db, "gate size mismatch");
        assert_eq!(g.cols, da * db);
        if a > b {
            // Reorder the gate to the (b, a) orientation and recurse.
            let g2 = Mat::from_fn(da * db, da * db, |row, col| {
                let (y_out, x_out) = (row / da, row % da);
                let (y_in, x_in) = (col / da, col % da);
                g.at(x_out * db + y_out, x_in * db + y_in)
            });
            return self.apply2(b, a, &g2);
        }
        if b == a + 1 {
            self.apply2_adjacent(a, g);
        } else {
            self.apply2_long_range(a, b, g);
        }
    }

    /// Two-site gate on neighbouring sites `(i, i+1)` with SVD re-splitting.
    fn apply2_adjacent(&mut self, i: usize, g: &Mat) {
        self.move_center_to(i);
        let (t1, t2) = (&self.tensors[i], &self.tensors[i + 1]);
        let (dl, d1, bond) = (t1.dl, t1.d, t1.dr);
        let (d2, dr) = (t2.d, t2.dr);
        // theta[(p1·d2+p2), (l·dr+r)]
        let mut theta = Mat::zeros(d1 * d2, dl * dr);
        for l in 0..dl {
            for p1 in 0..d1 {
                for bb in 0..bond {
                    let v1 = t1.at(l, p1, bb);
                    if v1.re == 0.0 && v1.im == 0.0 {
                        continue;
                    }
                    for p2 in 0..d2 {
                        for r in 0..dr {
                            theta.add_at(p1 * d2 + p2, l * dr + r, v1 * t2.at(bb, p2, r));
                        }
                    }
                }
            }
        }
        let theta2 = g.mul(&theta);
        // M[(l·d1+q1), (q2·dr+r)] = theta2[(q1·d2+q2), (l·dr+r)]
        let m = Mat::from_fn(dl * d1, d2 * dr, |row, col| {
            let (l, q1) = (row / d1, row % d1);
            let (q2, r) = (col / dr, col % dr);
            theta2.at(q1 * d2 + q2, l * dr + r)
        });
        let ts = svd_trunc(&m, self.trunc);
        let k = ts.s.len();
        let renorm = self.absorb_discard(ts.discarded_weight);
        let mut u = SiteTensor::zeros(dl, d1, k);
        for l in 0..dl {
            for p in 0..d1 {
                for r in 0..k {
                    u.set(l, p, r, ts.u.at(l * d1 + p, r));
                }
            }
        }
        let mut sv = SiteTensor::zeros(k, d2, dr);
        for a in 0..k {
            for p in 0..d2 {
                for r in 0..dr {
                    sv.set(a, p, r, ts.vh.at(a, p * dr + r).scale(ts.s[a] * renorm));
                }
            }
        }
        self.tensors[i] = u;
        self.tensors[i + 1] = sv;
        self.center = i + 1;
    }

    /// Long-range gate on `(a, b)` with `a < b - 1`: operator-Schmidt
    /// decomposition `G = Σ_κ A_κ ⊗ B_κ` threaded through the intervening
    /// sites as an MPO, followed by a two-pass (exact →, truncating ←)
    /// compression sweep over the affected window.
    fn apply2_long_range(&mut self, a: usize, b: usize, g: &Mat) {
        let (da, db) = (self.dims[a], self.dims[b]);
        self.move_center_to(a);
        // Operator-Schmidt across the (a | b) operator bipartition:
        // M[(x_out·da + x_in), (y_out·db + y_in)] = G[(x_out·db+y_out), (x_in·db+y_in)]
        let m = Mat::from_fn(da * da, db * db, |row, col| {
            let (x_out, x_in) = (row / da, row % da);
            let (y_out, y_in) = (col / db, col % db);
            g.at(x_out * db + y_out, x_in * db + y_in)
        });
        let os = svd(&m);
        let smax = os.s[0].max(1e-300);
        let kk =
            os.s.iter()
                .take_while(|&&s| s > 1e-14 * smax)
                .count()
                .max(1);
        // A_κ[x, x'] = U[(x·da+x'), κ]·√s_κ ; B_κ[y, y'] = √s_κ·Vh[κ, (y·db+y')]
        let a_ops: Vec<Mat> = (0..kk)
            .map(|k| {
                let sq = os.s[k].sqrt();
                Mat::from_fn(da, da, |x, xp| os.u.at(x * da + xp, k).scale(sq))
            })
            .collect();
        let b_ops: Vec<Mat> = (0..kk)
            .map(|k| {
                let sq = os.s[k].sqrt();
                Mat::from_fn(db, db, |y, yp| os.vh.at(k, y * db + yp).scale(sq))
            })
            .collect();

        // Site a: right bond gains the MPO index κ (major).
        {
            let t = &self.tensors[a];
            let (dl, dr) = (t.dl, t.dr);
            let mut new = SiteTensor::zeros(dl, da, kk * dr);
            for l in 0..dl {
                for x in 0..da {
                    for (k, ak) in a_ops.iter().enumerate() {
                        for r in 0..dr {
                            let mut acc = C64::ZERO;
                            for xp in 0..da {
                                acc += ak.at(x, xp) * t.at(l, xp, r);
                            }
                            new.set(l, x, k * dr + r, acc);
                        }
                    }
                }
            }
            self.tensors[a] = new;
        }
        // Intervening sites: κ passes through diagonally.
        for m_site in (a + 1)..b {
            let t = &self.tensors[m_site];
            let (dl, d, dr) = (t.dl, t.d, t.dr);
            let mut new = SiteTensor::zeros(kk * dl, d, kk * dr);
            for k in 0..kk {
                for l in 0..dl {
                    for p in 0..d {
                        for r in 0..dr {
                            new.set(k * dl + l, p, k * dr + r, t.at(l, p, r));
                        }
                    }
                }
            }
            self.tensors[m_site] = new;
        }
        // Site b: κ terminates against B_κ.
        {
            let t = &self.tensors[b];
            let (dl, dr) = (t.dl, t.dr);
            let mut new = SiteTensor::zeros(kk * dl, db, dr);
            for (k, bk) in b_ops.iter().enumerate() {
                for l in 0..dl {
                    for y in 0..db {
                        for r in 0..dr {
                            let mut acc = C64::ZERO;
                            for yp in 0..db {
                                acc += bk.at(y, yp) * t.at(l, yp, r);
                            }
                            new.set(k * dl + l, y, r, acc);
                        }
                    }
                }
            }
            self.tensors[b] = new;
        }
        // Compression: exact left→right pass, then truncating right→left.
        for m_site in a..b {
            self.push_right(m_site, TruncSpec::exact());
        }
        let spec = self.trunc;
        for m_site in ((a + 1)..=b).rev() {
            self.push_left(m_site, spec);
        }
        self.center = a;
    }

    // -- inner products & diagnostics --------------------------------------

    /// `⟨self|other⟩` by transfer-matrix contraction (cost `O(n·χ³·d)`).
    pub fn inner(&self, other: &Mps) -> C64 {
        assert_eq!(self.dims, other.dims);
        let mut env = Mat::identity(1);
        for (ta, tb) in self.tensors.iter().zip(other.tensors.iter()) {
            let mut new_env = Mat::zeros(ta.dr, tb.dr);
            for p in 0..ta.d {
                let sa = ta.phys_slice(p); // dl_a × dr_a
                let sb = tb.phys_slice(p); // dl_b × dr_b
                let t = sa.adjoint().mul(&env).mul(&sb); // dr_a × dr_b
                for r in 0..new_env.rows {
                    for c in 0..new_env.cols {
                        new_env.add_at(r, c, t.at(r, c));
                    }
                }
            }
            env = new_env;
        }
        env.at(0, 0)
    }

    pub fn norm(&self) -> f64 {
        self.inner(self).re.max(0.0).sqrt()
    }

    pub fn normalize(&mut self) {
        let n = self.norm();
        assert!(n > 0.0, "cannot normalize the zero state");
        let inv = C64::real(1.0 / n);
        for v in &mut self.tensors[self.center].data {
            *v *= inv;
        }
    }

    /// `|⟨self|other⟩|²` (assumes both normalized).
    pub fn fidelity(&self, other: &Mps) -> f64 {
        self.inner(other).abs2()
    }

    /// Contract the network into a dense state vector (exponential — only
    /// for small chains / verification).
    pub fn to_dense(&self) -> DenseState {
        let mut acc = Mat::identity(1); // (prefix-dim × bond)
        for t in &self.tensors {
            let (dl, d, dr) = (t.dl, t.d, t.dr);
            debug_assert_eq!(acc.cols, dl);
            let mut next = Mat::zeros(acc.rows * d, dr);
            for gi in 0..acc.rows {
                for l in 0..dl {
                    let av = acc.at(gi, l);
                    if av.re == 0.0 && av.im == 0.0 {
                        continue;
                    }
                    for p in 0..d {
                        for r in 0..dr {
                            next.add_at(gi * d + p, r, av * t.at(l, p, r));
                        }
                    }
                }
            }
            acc = next;
        }
        debug_assert_eq!(acc.cols, 1);
        DenseState {
            dims: self.dims.clone(),
            amps: acc.data,
        }
    }

    /// The sum `self + other` as a state vector (not normalized), built by
    /// the standard direct-sum construction — bond dimensions add — and
    /// recompressed under `trunc`. The workhorse behind operator linear
    /// combinations (Kraus sums, traced boundaries, rank-one corrections).
    pub fn add(&self, other: &Mps, trunc: TruncSpec) -> Mps {
        assert_eq!(
            self.dims, other.dims,
            "cannot add states on different chains"
        );
        let n = self.num_sites();
        if n == 1 {
            let (a, b) = (&self.tensors[0], &other.tensors[0]);
            let mut t = SiteTensor::zeros(1, a.d, 1);
            for p in 0..a.d {
                t.set(0, p, 0, a.at(0, p, 0) + b.at(0, p, 0));
            }
            return Mps {
                dims: self.dims.clone(),
                tensors: vec![t],
                center: 0,
                trunc,
                discarded_weight: 0.0,
            };
        }
        let mut tensors = Vec::with_capacity(n);
        for (i, (a, b)) in self.tensors.iter().zip(other.tensors.iter()).enumerate() {
            let first = i == 0;
            let last = i == n - 1;
            let dl = if first { 1 } else { a.dl + b.dl };
            let dr = if last { 1 } else { a.dr + b.dr };
            let mut t = SiteTensor::zeros(dl, a.d, dr);
            for p in 0..a.d {
                for l in 0..a.dl {
                    for r in 0..a.dr {
                        t.set(l, p, r, a.at(l, p, r));
                    }
                }
                let (lo, ro) = (if first { 0 } else { a.dl }, if last { 0 } else { a.dr });
                for l in 0..b.dl {
                    for r in 0..b.dr {
                        let cur = t.at(lo + l, p, ro + r);
                        t.set(lo + l, p, ro + r, cur + b.at(l, p, r));
                    }
                }
            }
            tensors.push(t);
        }
        let mut out = Mps {
            dims: self.dims.clone(),
            tensors,
            center: 0,
            trunc,
            discarded_weight: 0.0,
        };
        out.recompress();
        out
    }

    /// Marginal probabilities of the basis levels of `site` (assumes the
    /// state is normalized). Moves the orthogonality center.
    pub fn site_probabilities(&mut self, site: usize) -> Vec<f64> {
        self.move_center_to(site);
        let t = &self.tensors[site];
        (0..t.d)
            .map(|p| {
                let mut acc = 0.0;
                for l in 0..t.dl {
                    for r in 0..t.dr {
                        acc += t.at(l, p, r).abs2();
                    }
                }
                acc
            })
            .collect()
    }

    /// Normalized Schmidt spectra across every bond (`n-1` entries). Leaves
    /// the state unchanged (up to canonical re-gauging); the center ends at
    /// the last site.
    pub fn schmidt_spectra(&mut self) -> Vec<Vec<f64>> {
        let n = self.num_sites();
        self.move_center_to(0);
        let mut out = Vec::with_capacity(n - 1);
        for m in 0..n - 1 {
            // SVD of the center matricized as [(l·d) × r]: its singular
            // values are the Schmidt coefficients across bond (m, m+1).
            self.push_right(m, TruncSpec::exact());
            self.center = m + 1;
            // The Schmidt values were absorbed into m+1; recover them from
            // the new center by an SVD of its left matricization.
            let t = &self.tensors[m + 1];
            let mm = Mat::from_fn(t.dl, t.d * t.dr, |l, col| t.at(l, col / t.dr, col % t.dr));
            let dec = svd(&mm);
            let total: f64 = dec.s.iter().map(|s| s * s).sum();
            let norm = total.sqrt().max(1e-300);
            out.push(dec.s.iter().map(|s| s / norm).collect());
        }
        out
    }

    /// Von Neumann entanglement entropy in **bits** across a *single* bond
    /// `m`, by canonicalizing to it and taking one SVD — `O(χ³)`, versus
    /// the `O(n·χ³)` full sweep of [`Mps::bond_entropies_bits`]. Use this
    /// when only one cut's entropy is wanted (e.g. an inter-strand A|B
    /// bond, [`crate::crossing`]). Leaves the center at `m`.
    pub fn bond_entropy_bits_at(&mut self, m: usize) -> f64 {
        assert!(m + 1 < self.num_sites(), "bond index out of range");
        self.move_center_to(m);
        let t = &self.tensors[m];
        let mm = Mat::from_fn(t.dl * t.d, t.dr, |row, r| t.at(row / t.d, row % t.d, r));
        let dec = svd(&mm);
        let total: f64 = dec.s.iter().map(|s| s * s).sum();
        let norm = total.sqrt().max(1e-300);
        dec.s
            .iter()
            .map(|s| {
                let p = (s / norm) * (s / norm);
                if p > 1e-300 {
                    -p * p.log2()
                } else {
                    0.0
                }
            })
            .sum::<f64>()
            .max(0.0)
    }

    /// Von Neumann entanglement entropy in **bits** across every bond.
    pub fn bond_entropies_bits(&mut self) -> Vec<f64> {
        self.schmidt_spectra()
            .iter()
            .map(|sv| {
                sv.iter()
                    .map(|s| {
                        let p = s * s;
                        if p > 1e-300 {
                            -p * p.log2()
                        } else {
                            0.0
                        }
                    })
                    .sum::<f64>()
                    .max(0.0)
            })
            .collect()
    }

    // -- scale morphing ----------------------------------------------------

    /// Fuse sites `i` and `i+1` into a single site of dimension
    /// `d_i · d_{i+1}` (exact). The coarse basis index is `p_i·d_{i+1} + p_{i+1}`.
    pub fn merge_sites(&mut self, i: usize) {
        assert!(i + 1 < self.num_sites());
        self.move_center_to(i);
        let (t1, t2) = (&self.tensors[i], &self.tensors[i + 1]);
        let (dl, d1, bond) = (t1.dl, t1.d, t1.dr);
        let (d2, dr) = (t2.d, t2.dr);
        let mut merged = SiteTensor::zeros(dl, d1 * d2, dr);
        for l in 0..dl {
            for p1 in 0..d1 {
                for bb in 0..bond {
                    let v1 = t1.at(l, p1, bb);
                    if v1.re == 0.0 && v1.im == 0.0 {
                        continue;
                    }
                    for p2 in 0..d2 {
                        for r in 0..dr {
                            let v = merged.at(l, p1 * d2 + p2, r) + v1 * t2.at(bb, p2, r);
                            merged.set(l, p1 * d2 + p2, r, v);
                        }
                    }
                }
            }
        }
        self.tensors[i] = merged;
        self.tensors.remove(i + 1);
        self.dims[i] *= d2;
        self.dims.remove(i + 1);
        // center stays at i: contraction of the center with a
        // right-canonical neighbour is again a valid center.
    }

    /// Factor site `i` (dimension `d1·d2`) into two sites of dimensions
    /// `(d1, d2)`, truncating the new bond with the state's policy. Inverse
    /// of [`Mps::merge_sites`] up to truncation.
    pub fn split_site(&mut self, i: usize, d1: usize, d2: usize) {
        assert_eq!(
            self.dims[i],
            d1 * d2,
            "split dimensions must factor the site"
        );
        self.move_center_to(i);
        let t = &self.tensors[i];
        let (dl, dr) = (t.dl, t.dr);
        let m = Mat::from_fn(dl * d1, d2 * dr, |row, col| {
            let (l, p1) = (row / d1, row % d1);
            let (p2, r) = (col / dr, col % dr);
            t.at(l, p1 * d2 + p2, r)
        });
        let ts = svd_trunc(&m, self.trunc);
        let k = ts.s.len();
        let renorm = self.absorb_discard(ts.discarded_weight);
        let mut u = SiteTensor::zeros(dl, d1, k);
        for l in 0..dl {
            for p in 0..d1 {
                for r in 0..k {
                    u.set(l, p, r, ts.u.at(l * d1 + p, r));
                }
            }
        }
        let mut sv = SiteTensor::zeros(k, d2, dr);
        for a in 0..k {
            for p in 0..d2 {
                for r in 0..dr {
                    sv.set(a, p, r, ts.vh.at(a, p * dr + r).scale(ts.s[a] * renorm));
                }
            }
        }
        self.tensors[i] = u;
        self.tensors.insert(i + 1, sv);
        self.dims[i] = d1;
        self.dims.insert(i + 1, d2);
        self.center = i + 1;
    }

    /// Grow the local dimension of `site` to `d_new ≥ d` by isometric
    /// embedding (new levels start unpopulated). Exact; the canonical
    /// structure is untouched.
    pub fn promote_site(&mut self, site: usize, d_new: usize) {
        let d = self.dims[site];
        assert!(d_new >= d);
        if d_new == d {
            return;
        }
        let t = &self.tensors[site];
        let mut new = SiteTensor::zeros(t.dl, d_new, t.dr);
        for l in 0..t.dl {
            for p in 0..d {
                for r in 0..t.dr {
                    new.set(l, p, r, t.at(l, p, r));
                }
            }
        }
        self.tensors[site] = new;
        self.dims[site] = d_new;
    }

    /// Shrink the local dimension of `site` to `d_new ≤ d` by projecting out
    /// the upper levels, renormalizing, and returning the projected-out
    /// probability weight (0 when the state had no support there).
    pub fn demote_site(&mut self, site: usize, d_new: usize) -> f64 {
        let d = self.dims[site];
        assert!(d_new <= d && d_new > 0);
        if d_new == d {
            return 0.0;
        }
        self.move_center_to(site);
        let before = self.tensors[site].frobenius();
        let t = &self.tensors[site];
        let mut new = SiteTensor::zeros(t.dl, d_new, t.dr);
        for l in 0..t.dl {
            for p in 0..d_new {
                for r in 0..t.dr {
                    new.set(l, p, r, t.at(l, p, r));
                }
            }
        }
        let after = new.frobenius();
        assert!(
            after > 1e-150 * before.max(1e-150),
            "state has no support on the retained levels"
        );
        let rescale = before / after;
        for v in &mut new.data {
            *v = v.scale(rescale);
        }
        self.tensors[site] = new;
        self.dims[site] = d_new;
        1.0 - (after / before) * (after / before)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dense::DenseState;
    use crate::gates;
    use crate::mat::Rng;

    fn random_dense(dims: &[usize], rng: &mut Rng) -> DenseState {
        let mut s = DenseState::zero_state(dims);
        for a in &mut s.amps {
            *a = rng.c_gaussian();
        }
        s.normalize();
        s
    }

    #[test]
    fn basis_state_roundtrip() {
        let dims = [2, 3, 4, 5];
        let digits = [1, 2, 0, 4];
        let m = Mps::basis_state(&dims, &digits, TruncSpec::exact());
        let d = m.to_dense();
        let expect = DenseState::basis_state(&dims, &digits);
        assert!((d.fidelity(&expect) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn from_dense_to_dense_roundtrip() {
        let mut rng = Rng::new(5);
        let dims = [2, 3, 4, 3];
        let s = random_dense(&dims, &mut rng);
        let m = Mps::from_dense(&s, TruncSpec::exact());
        assert!((m.norm() - 1.0).abs() < 1e-10);
        assert!((m.to_dense().fidelity(&s) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn random_circuit_matches_dense_exactly() {
        let dims = [2, 3, 4, 3, 2];
        let mut rng = Rng::new(31);
        let mut mps = Mps::zero_state(&dims, TruncSpec::exact());
        let mut dense = DenseState::zero_state(&dims);
        for step in 0..25 {
            match step % 3 {
                0 => {
                    let s = rng.below(dims.len());
                    let g = gates::random(dims[s], &mut rng);
                    mps.apply1(s, &g);
                    dense.apply1(s, &g);
                }
                1 => {
                    let s = rng.below(dims.len() - 1);
                    let g = gates::random(dims[s] * dims[s + 1], &mut rng);
                    mps.apply2(s, s + 1, &g);
                    dense.apply2(s, s + 1, &g);
                }
                _ => {
                    // Long-range pair, both orientations.
                    let a = rng.below(dims.len());
                    let mut b = rng.below(dims.len());
                    while b == a {
                        b = rng.below(dims.len());
                    }
                    let g = gates::random(dims[a] * dims[b], &mut rng);
                    mps.apply2(a, b, &g);
                    dense.apply2(a, b, &g);
                }
            }
        }
        assert!((mps.norm() - 1.0).abs() < 1e-10);
        let f = mps.to_dense().fidelity(&dense);
        assert!((f - 1.0).abs() < 1e-10, "fidelity {}", f);
        assert!(mps.discarded_weight < 1e-18);
    }

    #[test]
    fn move_center_preserves_state() {
        let dims = [2, 3, 4];
        let mut rng = Rng::new(77);
        let s = random_dense(&dims, &mut rng);
        let mut m = Mps::from_dense(&s, TruncSpec::exact());
        for &t in &[0usize, 2, 1, 0, 2] {
            m.move_center_to(t);
            assert!((m.to_dense().fidelity(&s) - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn single_bond_entropy_matches_full_sweep() {
        let dims = [2, 3, 4, 3, 2];
        let mut rng = Rng::new(88);
        let s = random_dense(&dims, &mut rng);
        let full = {
            let mut m = Mps::from_dense(&s, TruncSpec::exact());
            m.bond_entropies_bits()
        };
        for bond in 0..dims.len() - 1 {
            let mut m = Mps::from_dense(&s, TruncSpec::exact());
            let one = m.bond_entropy_bits_at(bond);
            assert!(
                (one - full[bond]).abs() < 1e-9,
                "bond {}: {} vs {}",
                bond,
                one,
                full[bond]
            );
        }
    }

    #[test]
    fn bell_state_entropy_is_one_bit() {
        let dims = [2, 2];
        let mut m = Mps::zero_state(&dims, TruncSpec::exact());
        m.apply1(0, &gates::hadamard());
        m.apply2(0, 1, &gates::cshift(2, 2));
        let e = m.bond_entropies_bits();
        assert!((e[0] - 1.0).abs() < 1e-10, "entropy {}", e[0]);
    }

    #[test]
    fn generalized_bell_entropy_is_log2_d() {
        for d in [3usize, 5] {
            let dims = [d, d];
            let mut m = Mps::zero_state(&dims, TruncSpec::exact());
            m.apply1(0, &gates::fourier(d));
            m.apply2(0, 1, &gates::cshift(d, d));
            let e = m.bond_entropies_bits();
            assert!(
                (e[0] - (d as f64).log2()).abs() < 1e-10,
                "d={} entropy {}",
                d,
                e[0]
            );
        }
    }

    #[test]
    fn truncation_caps_bond_and_reports_weight() {
        let dims = [2, 2];
        let mut m = Mps::zero_state(&dims, TruncSpec::new(1, 0.0));
        m.apply1(0, &gates::hadamard());
        m.apply2(0, 1, &gates::cshift(2, 2)); // would create a Bell pair
        assert_eq!(m.max_bond_dim(), 1);
        assert!((m.discarded_weight - 0.5).abs() < 1e-10);
        assert!((m.norm() - 1.0).abs() < 1e-10, "renormalized after cut");
    }

    #[test]
    fn merge_split_roundtrip() {
        let dims = [2, 3, 4, 3];
        let mut rng = Rng::new(9);
        let s = random_dense(&dims, &mut rng);
        let mut m = Mps::from_dense(&s, TruncSpec::exact());
        m.merge_sites(1); // dims -> [2, 12, 3]
        assert_eq!(m.dims, vec![2, 12, 3]);
        m.split_site(1, 3, 4); // back to [2, 3, 4, 3]
        assert_eq!(m.dims, vec![2, 3, 4, 3]);
        assert!((m.to_dense().fidelity(&s) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn coarse_gate_equals_fine_gate_pair() {
        // A gate applied to a merged (2,3) site must equal the same gate
        // applied as a two-site gate before merging.
        let dims = [2, 3, 2];
        let mut rng = Rng::new(15);
        let s = random_dense(&dims, &mut rng);
        let g = gates::random(6, &mut rng);

        let mut fine = Mps::from_dense(&s, TruncSpec::exact());
        fine.apply2(0, 1, &g);

        let mut coarse = Mps::from_dense(&s, TruncSpec::exact());
        coarse.merge_sites(0); // dims [6, 2]
        coarse.apply1(0, &g);
        coarse.split_site(0, 2, 3);

        assert!((coarse.fidelity(&fine) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn promote_demote_roundtrip() {
        let dims = [3, 4];
        let mut rng = Rng::new(19);
        let s = random_dense(&dims, &mut rng);
        let mut m = Mps::from_dense(&s, TruncSpec::exact());
        m.promote_site(0, 5);
        assert_eq!(m.dims, vec![5, 4]);
        // State untouched in the embedded subspace:
        let leak = m.demote_site(0, 3);
        assert!(leak < 1e-14);
        assert!((m.to_dense().fidelity(&s) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn demote_reports_leakage() {
        // Uniform superposition on a qutrit: cutting to a qubit loses 1/3.
        let mut m = Mps::zero_state(&[3, 2], TruncSpec::exact());
        m.apply1(0, &gates::fourier(3));
        let leak = m.demote_site(0, 2);
        assert!((leak - 1.0 / 3.0).abs() < 1e-10, "leak {}", leak);
        assert!((m.norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn long_range_gate_is_norm_preserving_with_null_bonds() {
        // Regression: an MPO-applied long-range gate over a low-rank state
        // produces matricizations with large null spaces; a Jacobi-SVD
        // pathology used to inflate the norm by ~1e-5 here.
        let mut rng = Rng::new(123);
        for dims in [vec![3usize, 5, 4], vec![3, 4, 5, 4], vec![2, 3, 4, 5, 4]] {
            let n = dims.len();
            let g = gates::random(dims[0] * dims[n - 1], &mut rng);
            let mut psi = Mps::zero_state(&dims, TruncSpec::exact());
            for (i, &d) in dims.iter().enumerate() {
                psi.apply1(i, &gates::fourier(d));
            }
            let mut dense = psi.to_dense();
            psi.apply2(0, n - 1, &g);
            dense.apply2(0, n - 1, &g);
            assert!((psi.norm() - 1.0).abs() < 1e-12, "dims {:?}", dims);
            let f = psi.to_dense().fidelity(&dense);
            assert!((f - 1.0).abs() < 1e-12, "dims {:?} fidelity {}", dims, f);
        }
    }

    #[test]
    fn inner_matches_dense_inner() {
        let dims = [2, 3, 3];
        let mut rng = Rng::new(23);
        let s1 = random_dense(&dims, &mut rng);
        let s2 = random_dense(&dims, &mut rng);
        let m1 = Mps::from_dense(&s1, TruncSpec::exact());
        let m2 = Mps::from_dense(&s2, TruncSpec::exact());
        let a = m1.inner(&m2);
        let b = s1.inner(&s2);
        assert!((a - b).abs() < 1e-10);
    }
}
