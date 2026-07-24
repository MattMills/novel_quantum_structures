//! Matrix product operators: circuit blocks as first-class objects — the
//! "circuit of circuits" layer.
//!
//! An [`Mpo`] represents an operator on the qudit chain, stored — via the
//! Choi/vectorization isomorphism — as an ordinary [`Mps`] over the *doubled*
//! chain whose site `i` has dimension `d_i²` (flattened index
//! `p_out·d + p_in`). Read physically: each site of the carrier state holds a
//! (future, past) leg pair, so the operator is literally **a state stretched
//! between the past boundary and the future boundary of a circuit block**,
//! both boundaries being full `n`-site slices rather than single states.
//!
//! Consequences of that reading, all measurable here:
//!
//! * The carrier's bond spectrum across a cut is the **operator
//!   entanglement** of the block — how wide the "pipeline of expectations"
//!   between past and future must be at that point in space
//!   ([`Mpo::operator_entanglement_bits`]).
//! * Applying the block to a state ([`Mpo::apply_to`]) and composing two
//!   blocks ([`Mpo::compose_after`]) are both bond-contractions followed by
//!   recompression — gate composition lifted to *circuit-block composition*.
//!   Recompression is where the interesting physics happens: a pipeline of
//!   blocks can collapse to something far simpler than its parts (see the
//!   Fourier-space adder in `examples/circuit_of_circuits.rs`).
//! * A whole compiled block executes against a state in one MPO application:
//!   the block *is* a gate at the next scale up.
//!
//! Compilation ([`Mpo::from_circuit`]) reuses the full heterogeneous MPS
//! machinery (long-range gates included) by lifting each gate `G` to its
//! left-multiplication superoperator `G ⊗ I` on the doubled chain.
//!
//! Note on truncation semantics: bond truncation on the carrier preserves
//! the carrier's Frobenius norm (inherited from [`Mps`]); for the near-exact
//! compilations used here the rescaling is negligible.

use crate::c64::C64;
use crate::circuit::{Circuit, Op};
use crate::mat::{Mat, TruncSpec};
use crate::mps::{Mps, SiteTensor};

/// A matrix product operator over a heterogeneous qudit chain.
#[derive(Clone)]
pub struct Mpo {
    /// Physical dimensions of the chain the operator acts on.
    pub dims: Vec<usize>,
    /// Choi carrier: an MPS over dimensions `d_i²` (index `p_out·d + p_in`).
    pub carrier: Mps,
}

/// Lift a one-site gate `g` to its left-multiplication superoperator
/// `g ⊗ I` on the doubled index `p_out·d + p_in`.
fn lift1(g: &Mat, d: usize) -> Mat {
    Mat::from_fn(d * d, d * d, |row, col| {
        let (po, pi) = (row / d, row % d);
        let (qo, qi) = (col / d, col % d);
        if pi == qi {
            g.at(po, qo)
        } else {
            C64::ZERO
        }
    })
}

/// Lift a two-site gate on `(da, db)` to the doubled pair `(da², db²)`.
fn lift2(g: &Mat, da: usize, db: usize) -> Mat {
    let n = da * da * db * db;
    Mat::from_fn(n, n, |row, col| {
        // row = (xa·da + ia)·db² + (xb·db + ib)
        let (ra, rb) = (row / (db * db), row % (db * db));
        let (xa, ia) = (ra / da, ra % da);
        let (xb, ib) = (rb / db, rb % db);
        let (ca, cb) = (col / (db * db), col % (db * db));
        let (ya, ja) = (ca / da, ca % da);
        let (yb, jb) = (cb / db, cb % db);
        if ia == ja && ib == jb {
            g.at(xa * db + xb, ya * db + yb)
        } else {
            C64::ZERO
        }
    })
}

impl Mpo {
    /// The identity operator (carrier bond dimension 1).
    pub fn identity(dims: &[usize], trunc: TruncSpec) -> Mpo {
        let carrier_dims: Vec<usize> = dims.iter().map(|&d| d * d).collect();
        let tensors: Vec<SiteTensor> = dims
            .iter()
            .map(|&d| {
                let mut t = SiteTensor::zeros(1, d * d, 1);
                for p in 0..d {
                    t.set(0, p * d + p, 0, C64::ONE);
                }
                t
            })
            .collect();
        Mpo {
            dims: dims.to_vec(),
            carrier: Mps {
                dims: carrier_dims,
                tensors,
                center: 0,
                trunc,
                discarded_weight: 0.0,
            },
        }
    }

    /// Compile a circuit into an operator by absorbing its gates in order.
    pub fn from_circuit(c: &Circuit, trunc: TruncSpec) -> Mpo {
        let mut m = Mpo::identity(&c.dims, trunc);
        for op in &c.ops {
            m.absorb_after(op);
        }
        m
    }

    /// The exact MPO of an arbitrary permutation `x → π(x)` of the chain's
    /// ring, built through its dense Choi carrier and compressed exactly —
    /// so its bond dimensions *are* the permutation's operator Schmidt
    /// ranks (pinned against [`crate::width::perm_cut_rank`] by the
    /// tests). Verification-grade: cost and memory are `O(N²)` in the ring
    /// size `N = Π dᵢ`, so this is for small chains; structured
    /// permutations at scale go through [`crate::cascade::Transducer`].
    pub fn from_permutation(
        dims: &[usize],
        perm: impl Fn(usize) -> usize,
        trunc: TruncSpec,
    ) -> Mpo {
        use crate::dense::DenseState;
        let n_sites = dims.len();
        let doubled: Vec<usize> = dims.iter().map(|&d| d * d).collect();
        let mut st = DenseState::zero_state(&doubled);
        st.amps[0] = C64::ZERO;
        let n: usize = dims.iter().product();
        let mut hit = vec![false; n];
        for x in 0..n {
            let y = perm(x);
            assert!(y < n && !hit[y], "not a permutation of the chain's ring");
            hit[y] = true;
            let (mut xv, mut yv) = (x, y);
            let mut digs = vec![0usize; n_sites];
            for i in (0..n_sites).rev() {
                let d = dims[i];
                digs[i] = (yv % d) * d + (xv % d);
                xv /= d;
                yv /= d;
            }
            let mut idx = 0usize;
            for (i, &dd) in doubled.iter().enumerate() {
                idx = idx * dd + digs[i];
            }
            st.amps[idx] = C64::ONE;
        }
        Mpo {
            dims: dims.to_vec(),
            carrier: Mps::from_dense(&st, trunc),
        }
    }

    /// Absorb one more gate applied *after* the current operator
    /// (`self ← G ∘ self`).
    pub fn absorb_after(&mut self, op: &Op) {
        match op {
            Op::One(s, g) => {
                let d = self.dims[*s];
                self.carrier.apply1(*s, &lift1(g, d));
            }
            Op::Two(a, b, g) => {
                let (da, db) = (self.dims[*a], self.dims[*b]);
                self.carrier.apply2(*a, *b, &lift2(g, da, db));
            }
        }
    }

    pub fn num_sites(&self) -> usize {
        self.dims.len()
    }

    pub fn bond_dims(&self) -> Vec<usize> {
        self.carrier.bond_dims()
    }

    pub fn max_bond_dim(&self) -> usize {
        self.carrier.max_bond_dim()
    }

    pub fn param_count(&self) -> usize {
        self.carrier.param_count()
    }

    /// Apply the operator to a state: site-wise bond contraction followed by
    /// recompression under the *state's* truncation policy.
    pub fn apply_to(&self, psi: &Mps) -> Mps {
        assert_eq!(self.dims, psi.dims, "operator/state chain mismatch");
        let n = self.num_sites();
        let mut tensors = Vec::with_capacity(n);
        for (w, a) in self.carrier.tensors.iter().zip(psi.tensors.iter()) {
            let d = (w.d as f64).sqrt() as usize;
            debug_assert_eq!(d * d, w.d);
            debug_assert_eq!(a.d, d);
            let (lw, rw) = (w.dl, w.dr);
            let (la, ra) = (a.dl, a.dr);
            let mut t = SiteTensor::zeros(lw * la, d, rw * ra);
            for i_lw in 0..lw {
                for i_la in 0..la {
                    let l = i_lw * la + i_la;
                    for po in 0..d {
                        for i_rw in 0..rw {
                            for i_ra in 0..ra {
                                let mut acc = C64::ZERO;
                                for pi in 0..d {
                                    acc += w.at(i_lw, po * d + pi, i_rw) * a.at(i_la, pi, i_ra);
                                }
                                t.set(l, po, i_rw * ra + i_ra, acc);
                            }
                        }
                    }
                }
            }
            tensors.push(t);
        }
        let mut out = Mps {
            dims: psi.dims.clone(),
            tensors,
            center: 0,
            trunc: psi.trunc,
            discarded_weight: psi.discarded_weight,
        };
        out.recompress();
        out
    }

    /// Operator composition `self ∘ inner` (apply `inner` first), with the
    /// result recompressed under `trunc` — circuit-block composition.
    pub fn compose_after(&self, inner: &Mpo, trunc: TruncSpec) -> Mpo {
        assert_eq!(self.dims, inner.dims, "operator chain mismatch");
        let n = self.num_sites();
        let mut tensors = Vec::with_capacity(n);
        for ((wo, wi), &d) in self
            .carrier
            .tensors
            .iter()
            .zip(inner.carrier.tensors.iter())
            .zip(self.dims.iter())
        {
            let (lo, ro) = (wo.dl, wo.dr);
            let (li, ri) = (wi.dl, wi.dr);
            let mut t = SiteTensor::zeros(lo * li, d * d, ro * ri);
            for a_lo in 0..lo {
                for a_li in 0..li {
                    let l = a_lo * li + a_li;
                    for po in 0..d {
                        for pi in 0..d {
                            for a_ro in 0..ro {
                                for a_ri in 0..ri {
                                    let mut acc = C64::ZERO;
                                    for m in 0..d {
                                        acc += wo.at(a_lo, po * d + m, a_ro)
                                            * wi.at(a_li, m * d + pi, a_ri);
                                    }
                                    t.set(l, po * d + pi, a_ro * ri + a_ri, acc);
                                }
                            }
                        }
                    }
                }
            }
            tensors.push(t);
        }
        let mut carrier = Mps {
            dims: self.carrier.dims.clone(),
            tensors,
            center: 0,
            trunc,
            discarded_weight: 0.0,
        };
        carrier.recompress();
        Mpo {
            dims: self.dims.clone(),
            carrier,
        }
    }

    /// Operator sum `self + other`, recompressed under `trunc`. Together
    /// with [`Mpo::scale`] this makes operator *linear combinations*
    /// first-class: Kraus sums, rank-one corrections, traced boundaries.
    pub fn add(&self, other: &Mpo, trunc: TruncSpec) -> Mpo {
        assert_eq!(self.dims, other.dims, "operator chain mismatch");
        Mpo {
            dims: self.dims.clone(),
            carrier: self.carrier.add(&other.carrier, trunc),
        }
    }

    /// Scalar multiple `c · self`.
    pub fn scale(&self, c: C64) -> Mpo {
        let mut carrier = self.carrier.clone();
        for v in &mut carrier.tensors[0].data {
            *v *= c;
        }
        Mpo {
            dims: self.dims.clone(),
            carrier,
        }
    }

    /// The rank-one basis transfer operator `|to⟩⟨from|` — a bond-1 product
    /// MPO (each site contributes `|to_i⟩⟨from_i|`). With [`Mpo::add`] and
    /// [`Mpo::scale`] this builds rank-one corrections such as damped
    /// projectors `I − (1−γ)|x⟩⟨x|`.
    pub fn basis_transfer(dims: &[usize], to: &[usize], from: &[usize], trunc: TruncSpec) -> Mpo {
        assert_eq!(dims.len(), to.len());
        assert_eq!(dims.len(), from.len());
        let tensors: Vec<SiteTensor> = dims
            .iter()
            .zip(to.iter().zip(from.iter()))
            .map(|(&d, (&t, &f))| {
                assert!(t < d && f < d);
                let mut w = SiteTensor::zeros(1, d * d, 1);
                w.set(0, t * d + f, 0, C64::ONE);
                w
            })
            .collect();
        Mpo {
            dims: dims.to_vec(),
            carrier: Mps {
                dims: dims.iter().map(|&d| d * d).collect(),
                tensors,
                center: 0,
                trunc,
                discarded_weight: 0.0,
            },
        }
    }

    /// The control-select operator `Σ_a |a⟩⟨a|_ctrl ∘ branches[a]`: acts
    /// as `branches[a]` on states whose `ctrl` site holds digit `a`.
    /// `branches.len()` must equal `dims[ctrl]`, and each branch should
    /// act as the identity on the control site (whatever it does there is
    /// composed with the projector). For a qubit control and branches
    /// `[I, U]` this is the controlled-`U`; with a qudit control and
    /// branches `[I, U, U², …]` it is the controlled power
    /// `|a⟩|x⟩ → |a⟩ U^a|x⟩` — the kernel of phase estimation. Combined
    /// with [`crate::cascade::Transducer::mult_skipping`] it makes
    /// controlled modular arithmetic first-class (see
    /// `examples/shor_kernel.rs`).
    pub fn select_on(ctrl: usize, branches: &[Mpo], trunc: TruncSpec) -> Mpo {
        assert!(!branches.is_empty());
        let dims = branches[0].dims.clone();
        assert!(ctrl < dims.len());
        let d = dims[ctrl];
        assert_eq!(branches.len(), d, "one branch per control level");
        let mut acc: Option<Mpo> = None;
        for (a, branch) in branches.iter().enumerate() {
            assert_eq!(branch.dims, dims, "branch chain mismatch");
            let proj = Mat::from_fn(d, d, |r, c| {
                if r == a && c == a {
                    C64::ONE
                } else {
                    C64::ZERO
                }
            });
            let mut term = branch.clone();
            term.absorb_after(&Op::One(ctrl, proj));
            acc = Some(match acc {
                None => term,
                Some(s) => s.add(&term, trunc),
            });
        }
        acc.expect("at least one branch")
    }

    /// Hermitian adjoint.
    pub fn adjoint(&self) -> Mpo {
        let mut carrier = self.carrier.clone();
        for (t, &d) in carrier.tensors.iter_mut().zip(self.dims.iter()) {
            let mut new = SiteTensor::zeros(t.dl, t.d, t.dr);
            for l in 0..t.dl {
                for po in 0..d {
                    for pi in 0..d {
                        for r in 0..t.dr {
                            new.set(l, po * d + pi, r, t.at(l, pi * d + po, r).conj());
                        }
                    }
                }
            }
            *t = new;
        }
        Mpo {
            dims: self.dims.clone(),
            carrier,
        }
    }

    /// Operator entanglement (in bits) across every bond: the entanglement
    /// entropy of the normalized Choi state — the width of the
    /// past↔future correlation pipeline at each cut.
    pub fn operator_entanglement_bits(&self) -> Vec<f64> {
        let mut c = self.carrier.clone();
        c.normalize();
        c.bond_entropies_bits()
    }

    /// Hilbert–Schmidt inner product `tr(self† · other)` via the carriers.
    pub fn hs_inner(&self, other: &Mpo) -> C64 {
        self.carrier.inner(&other.carrier)
    }

    /// Normalized Hilbert–Schmidt overlap
    /// `|tr(self†·other)|² / (‖self‖² ‖other‖²)` — 1 iff the operators are
    /// equal up to a scalar.
    pub fn hs_fidelity(&self, other: &Mpo) -> f64 {
        let num = self.hs_inner(other).abs2();
        let na = self.carrier.inner(&self.carrier).re;
        let nb = other.carrier.inner(&other.carrier).re;
        num / (na * nb)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dense::DenseState;
    use crate::gates;
    use crate::mat::Rng;
    use crate::zigzag;

    fn random_mps(dims: &[usize], rng: &mut Rng) -> Mps {
        let mut s = DenseState::zero_state(dims);
        for a in &mut s.amps {
            *a = rng.c_gaussian();
        }
        s.normalize();
        Mps::from_dense(&s, TruncSpec::exact())
    }

    fn random_circuit(dims: &[usize], ops: usize, rng: &mut Rng) -> Circuit {
        let mut c = Circuit::new(dims.to_vec());
        for step in 0..ops {
            if step % 2 == 0 {
                let s = rng.below(dims.len());
                c.one(s, gates::random(dims[s], rng));
            } else {
                let a = rng.below(dims.len());
                let mut b = rng.below(dims.len());
                while b == a {
                    b = rng.below(dims.len());
                }
                c.two(a, b, gates::random(dims[a] * dims[b], rng));
            }
        }
        c
    }

    #[test]
    fn identity_mpo_is_identity() {
        let dims = [2usize, 3, 4];
        let mut rng = Rng::new(41);
        let id = Mpo::identity(&dims, TruncSpec::exact());
        assert_eq!(id.max_bond_dim(), 1);
        let psi = random_mps(&dims, &mut rng);
        let out = id.apply_to(&psi);
        assert!((out.fidelity(&psi) - 1.0).abs() < 1e-10);
        // Identity has zero operator entanglement everywhere.
        for e in id.operator_entanglement_bits() {
            assert!(e.abs() < 1e-10);
        }
    }

    #[test]
    fn compiled_circuit_acts_like_the_circuit() {
        let dims = [2usize, 3, 4, 3];
        let mut rng = Rng::new(43);
        let c = random_circuit(&dims, 12, &mut rng);
        let m = Mpo::from_circuit(&c, TruncSpec::exact());

        let psi = random_mps(&dims, &mut rng);
        // Route 1: gate-by-gate on the state.
        let mut direct = psi.clone();
        c.run_mps(&mut direct);
        // Route 2: one MPO application.
        let via_mpo = m.apply_to(&psi);
        assert!(
            (via_mpo.fidelity(&direct) - 1.0).abs() < 1e-8,
            "fidelity {}",
            via_mpo.fidelity(&direct)
        );
    }

    #[test]
    fn composition_matches_concatenation() {
        let dims = [2usize, 3, 3];
        let mut rng = Rng::new(47);
        let c1 = random_circuit(&dims, 6, &mut rng);
        let c2 = random_circuit(&dims, 6, &mut rng);
        let m1 = Mpo::from_circuit(&c1, TruncSpec::exact());
        let m2 = Mpo::from_circuit(&c2, TruncSpec::exact());
        // m2 ∘ m1 (c1 first) vs compiling c1 then c2.
        let composed = m2.compose_after(&m1, TruncSpec::exact());
        let mut concat = Circuit::new(dims.to_vec());
        for op in c1.ops.into_iter().chain(c2.ops.into_iter()) {
            concat.ops.push(op);
        }
        let reference = Mpo::from_circuit(&concat, TruncSpec::exact());
        assert!(
            (composed.hs_fidelity(&reference) - 1.0).abs() < 1e-8,
            "hs fidelity {}",
            composed.hs_fidelity(&reference)
        );
    }

    #[test]
    fn adjoint_inverts_unitary_blocks() {
        let dims = [2usize, 3, 2];
        let mut rng = Rng::new(53);
        let c = random_circuit(&dims, 8, &mut rng);
        let m = Mpo::from_circuit(&c, TruncSpec::exact());
        let inv = m.adjoint();
        let id = inv.compose_after(&m, TruncSpec::exact());
        let reference = Mpo::identity(&dims, TruncSpec::exact());
        assert!(
            (id.hs_fidelity(&reference) - 1.0).abs() < 1e-8,
            "U†U != I: hs fidelity {}",
            id.hs_fidelity(&reference)
        );
    }

    #[test]
    fn select_on_is_controlled_power() {
        // Qutrit control at site 0 over the sub-ring Z_15 (sites 1, 2):
        // |a, x⟩ → |a, 2^a·x mod 15⟩ with branches [I, ×2, ×4].
        use crate::cascade::Transducer;
        let dims = vec![3usize, 3, 5];
        let spec = TruncSpec::exact();
        let id = Mpo::identity(&dims, spec);
        let m2 = Transducer::mult_skipping(&dims, 0, 2).to_mpo(spec);
        let m4 = m2.compose_after(&m2, spec);
        let cm = Mpo::select_on(0, &[id, m2, m4], spec);
        for a in 0..3usize {
            for x in 0..15usize {
                let digits_of = |mut v: usize| -> Vec<usize> {
                    let mut d = vec![0usize; dims.len()];
                    for (i, &dim) in dims.iter().enumerate().rev() {
                        d[i] = v % dim;
                        v /= dim;
                    }
                    d
                };
                let input = Mps::basis_state(&dims, &digits_of(a * 15 + x), spec);
                let target =
                    Mps::basis_state(&dims, &digits_of(a * 15 + (1 << a) * x % 15), spec);
                let f = cm.apply_to(&input).fidelity(&target);
                assert!((f - 1.0).abs() < 1e-9, "a={} x={}: {}", a, x, f);
            }
        }
    }

    #[test]
    fn phase_kickback_is_bond_free() {
        // Control |+⟩, register in an eigenstate of ×7 on Z_12
        // (7² ≡ 1, orbit {1, 7}): |u⟩ = (|1⟩ − |7⟩)/√2 has U|u⟩ = −|u⟩,
        // so C-U sends |+⟩⊗|u⟩ to |−⟩⊗|u⟩ — a product state. The phase
        // lands on the control without a single unit of bond dimension
        // crossing the control cut.
        use crate::cascade::Transducer;
        use crate::gates;
        let dims = vec![2usize, 3, 4];
        let spec = TruncSpec::exact();
        let id = Mpo::identity(&dims, spec);
        let m7 = Transducer::mult_skipping(&dims, 0, 7).to_mpo(spec);
        let cm = Mpo::select_on(0, &[id, m7], spec);

        let sqrt_half = C64::real(1.0 / 2.0_f64.sqrt());
        let mut u = DenseState::zero_state(&dims);
        u.amps[0] = C64::ZERO;
        u.amps[1] = sqrt_half;
        u.amps[7] = C64::ZERO - sqrt_half;
        let mut plus_u = u.clone();
        plus_u.apply1(0, &gates::hadamard()); // |+⟩ ⊗ |u⟩
        let mut minus_u = plus_u.clone();
        minus_u.apply1(0, &gates::phase_diag(&[0.0, std::f64::consts::PI])); // |−⟩ ⊗ |u⟩

        let out = cm.apply_to(&Mps::from_dense(&plus_u, spec));
        let f = out.to_dense().fidelity(&minus_u);
        assert!((f - 1.0).abs() < 1e-9, "kickback fidelity {}", f);
        assert_eq!(out.bond_dims()[0], 1, "control bond must stay trivial");
    }

    #[test]
    fn controlled_width_formula() {
        // Across an interior cut with register split (L | S), the
        // controlled multiplier has width w + 1 when k ≢ 1 (mod S) and
        // exactly w when k ≡ 1 (mod S) — the control is free on resonant
        // powers. On [2,3,4] (sub-ring Z_12, cut (3|4), w = 3 for both):
        // ×7 (7 mod 4 = 3): bonds [2, 4]; ×5 (5 mod 4 = 1): bonds [2, 3].
        use crate::cascade::Transducer;
        let dims = vec![2usize, 3, 4];
        let spec = TruncSpec::exact();
        for (k, expect) in [(7usize, vec![2, 4]), (5, vec![2, 3])] {
            let id = Mpo::identity(&dims, spec);
            let m = Transducer::mult_skipping(&dims, 0, k).to_mpo(spec);
            assert_eq!(m.bond_dims(), vec![1, 3], "k={}", k);
            let cm = Mpo::select_on(0, &[id, m], spec);
            assert_eq!(cm.bond_dims(), expect, "k={}", k);
        }
    }

    #[test]
    fn from_permutation_matches_the_cascade() {
        // The dense-Choi route and the transducer route must produce the
        // same operator with the same exact bond profile.
        use crate::cascade::Transducer;
        let dims = vec![2usize, 3, 4]; // N = 24
        let spec = TruncSpec::exact();
        let via_dense = Mpo::from_permutation(&dims, |x| 7 * x % 24, spec);
        let via_cascade = Transducer::mult(&dims, 7).to_mpo(spec);
        assert_eq!(via_dense.bond_dims(), via_cascade.bond_dims());
        assert!((via_dense.hs_fidelity(&via_cascade) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn bowtie_block_operator_entanglement_is_finite_and_symmetric() {
        let profile = zigzag::diamond(2, 4);
        let m = Mpo::from_circuit(&zigzag::bowtie(&profile), TruncSpec::exact());
        let e = m.operator_entanglement_bits();
        let n = e.len();
        for i in 0..n / 2 {
            assert!(
                (e[i] - e[n - 1 - i]).abs() < 1e-8,
                "mirror symmetry broken: {:?}",
                e
            );
        }
        assert!(e.iter().all(|&x| x.is_finite()));
    }
}
