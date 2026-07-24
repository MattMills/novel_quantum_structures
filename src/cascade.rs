//! Stepwise cascade operators: finite-state **transducers** lifted to
//! quantum operators on the chain.
//!
//! A [`Transducer`] is a classical machine that sweeps the chain once,
//! carrying a message of dimension `m` on the moving front: at each site it
//! reads `(message_in, digit_in)`, writes `(message_out, digit_out)` by a
//! local reversible rule, and moves on. Lifted to an MPO, **the message
//! rides the virtual bond** — the operator's bond dimension *is* the width
//! of the classical information flowing through the cascade. The
//! carry-propagation adder ([`crate::radix::adder_mpo_exact`]) is the
//! `m = 2` member of this family; this module generalizes it.
//!
//! The flagship new member: **modular multiplication**. Schoolbook
//! multiplication by a constant `k` processes digits from least significant
//! to most with a multiplication carry in `[0, k)`:
//!
//! ```text
//!   v = k·digit_in + carry_in ;   digit_out = v mod d ;   carry_out = v div d
//! ```
//!
//! and `(digit_in, carry_in) ↔ v ↔ (carry_out, digit_out)` is a bijection at
//! every site, so the cascade is an exact permutation transducer with
//! message dimension `k`. Dropping the final carry is reduction mod `N`;
//! the resulting operator `|x⟩ → |k·x mod N⟩` is unitary **iff
//! `gcd(k, N) = 1`** — number theory surfacing as an operator property
//! (a violated gcd shows up as a norm anomaly, see the tests and
//! `examples/stepwise_cascade.rs`).
//!
//! Cascades compose recursively — `mult(k₁) ∘ mult(k₂) = mult(k₁k₂ mod N)`
//! — so modular exponentiation by repeated squaring becomes iterated
//! cascade composition, with the *measured true operator width* of each
//! product the honest cost meter (message dimension is only an upper
//! bound: the chain's geometric budget caps every bond).
//!
//! Direction convention: arithmetic carries flow from the least significant
//! site (`n-1`) toward the most significant (`0`) — [`Direction::RightToLeft`].

use crate::c64::C64;
use crate::mat::TruncSpec;
use crate::mpo::Mpo;
use crate::mps::{Mps, SiteTensor};

/// Which way the message front sweeps the chain.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Direction {
    /// Message enters at site 0, exits (dropped) after site `n-1`.
    LeftToRight,
    /// Message enters at site `n-1`, exits (dropped) after site 0 — the
    /// direction of arithmetic carries in this crate's digit convention.
    RightToLeft,
}

/// Boundary condition for the message front at either end of the sweep.
///
/// A classical machine starts in one state and its final state is read out
/// or discarded. An MPO boundary can do more: **enter in a superposition of
/// every state and postselect the exit** — the contraction keeps exactly the
/// branches whose entry guess turns out consistent. This is what lets the
/// geometrically opposed machine of [`Transducer::div`] run: the division
/// remainder machine guesses the wrap multiple `m` at entry (`SumAll`) and
/// keeps only exact divisions at exit (`Fixed(0)`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Boundary {
    /// The front enters/exits in a definite state.
    Fixed(usize),
    /// Sum over all front states at this end.
    SumAll,
}

/// A classical reversible transducer over the chain: per-site rules
/// `(message_in, digit_in) → (message_out, digit_out)`.
pub struct Transducer {
    pub profile: Vec<usize>,
    pub msg_dim: usize,
    /// `rules[i][msg_in * d_i + digit_in] = (msg_out, digit_out)`.
    pub rules: Vec<Vec<(usize, usize)>>,
    pub direction: Direction,
}

impl Transducer {
    /// The `+c mod N` adder: message = the additive carry bit (`m = 2`).
    pub fn adder(profile: &[usize], c_add: u128) -> Transducer {
        let n = profile.len();
        let big_n: u128 = profile.iter().map(|&d| d as u128).product();
        let mut rem = c_add % big_n;
        let mut c_digits = vec![0usize; n];
        for i in (0..n).rev() {
            c_digits[i] = (rem % profile[i] as u128) as usize;
            rem /= profile[i] as u128;
        }
        let rules = profile
            .iter()
            .enumerate()
            .map(|(i, &d)| {
                let mut r = vec![(0usize, 0usize); 2 * d];
                for msg_in in 0..2 {
                    for p_in in 0..d {
                        let v = p_in + c_digits[i] + msg_in;
                        r[msg_in * d + p_in] = (v / d, v % d);
                    }
                }
                r
            })
            .collect();
        Transducer {
            profile: profile.to_vec(),
            msg_dim: 2,
            rules,
            direction: Direction::RightToLeft,
        }
    }

    /// The `×k mod N` multiplier: message = the multiplication carry
    /// (`m = k`). Unitary iff `gcd(k, N) = 1`. Direct construction is
    /// limited to modest `k` (site tensors hold `k²·d²` entries); obtain
    /// larger multipliers by *composing* cascades.
    pub fn mult(profile: &[usize], k: usize) -> Transducer {
        Transducer::mult_directed(profile, k, Direction::RightToLeft)
    }

    /// [`Transducer::mult`] with an explicit sweep direction.
    /// `RightToLeft` multiplies the standard reading of the chain (LSB at
    /// site `n-1`); `LeftToRight` multiplies the reversed-place reading
    /// (LSB at site 0) — the encoding produced by
    /// [`crate::radix::mixed_radix_qft_reversed`], where the Fourier
    /// scaling theorem lands (see the tests: conjugating a multiplier by
    /// the frame inverts `k` *and reverses the cascade direction*).
    pub fn mult_directed(profile: &[usize], k: usize, direction: Direction) -> Transducer {
        assert!(k >= 1, "k must be positive");
        assert!(
            k <= 512,
            "direct construction is for modest k; compose cascades for larger multipliers"
        );
        let rules = profile
            .iter()
            .map(|&d| {
                let mut r = vec![(0usize, 0usize); k * d];
                for msg_in in 0..k {
                    for p_in in 0..d {
                        let v = k * p_in + msg_in;
                        r[msg_in * d + p_in] = (v / d, v % d);
                    }
                }
                r
            })
            .collect();
        Transducer {
            profile: profile.to_vec(),
            msg_dim: k,
            rules,
            direction,
        }
    }

    /// Modular **division**: the geometrically opposed machine to
    /// [`Transducer::mult`], computing the same inverse operation with a
    /// dual state set running the other way.
    ///
    /// `×k⁻¹ mod N` as a forward (LSB→MSB) carry machine needs message
    /// dimension `k⁻¹ mod N` — typically enormous. But long division's
    /// remainders flow **MSB→LSB**: with `y = k⁻¹x mod N` there is a unique
    /// wrap multiple `m ∈ [0, k)` with `k·y = x + m·N`, and dividing
    /// `x + m·N` by `k` most-significant-first needs only the remainder
    /// state `r ∈ [0, k)`:
    ///
    /// ```text
    ///   t = r·d + digit_in ;   digit_out = t div k ;   r' = t mod k
    /// ```
    ///
    /// The unknown `m` is exactly the *initial* remainder — so lift with
    /// [`Boundary::SumAll`] at entry (guess every `m`) and
    /// [`Boundary::Fixed`]`(0)` at exit (keep only exact divisions): for each
    /// input one branch survives, and the MPO is the exact permutation
    /// `|x⟩ → |k⁻¹·x mod N⟩` at message width `k` instead of `k⁻¹ mod N`.
    pub fn div(profile: &[usize], k: usize) -> Transducer {
        assert!(k >= 1, "k must be positive");
        assert!(
            k <= 512,
            "direct construction is for modest k; compose cascades for larger divisors"
        );
        let rules = profile
            .iter()
            .map(|&d| {
                let mut r = vec![(0usize, 0usize); k * d];
                for rem_in in 0..k {
                    for p_in in 0..d {
                        let t = rem_in * d + p_in;
                        r[rem_in * d + p_in] = (t % k, t / k);
                    }
                }
                r
            })
            .collect();
        Transducer {
            profile: profile.to_vec(),
            msg_dim: k,
            rules,
            direction: Direction::LeftToRight,
        }
    }

    /// The MPO of [`Transducer::div`] with its dual boundaries
    /// (`SumAll` entry, `Fixed(0)` exit): the exact operator
    /// `|x⟩ → |k⁻¹·x mod N⟩`, requiring `gcd(k, N) = 1`.
    pub fn div_mpo(profile: &[usize], k: usize, trunc: TruncSpec) -> Mpo {
        Transducer::div(profile, k).to_mpo_with(Boundary::SumAll, Boundary::Fixed(0), trunc)
    }

    /// The **looped boundary**: the exiting message is fed back into the
    /// entry — the machine eats its own tail, `M = Σ_m ⟨m|machine|m⟩`.
    ///
    /// Closing the loop changes the arithmetic: the traced adder is
    /// end-around carry (ones'-complement addition), computing
    /// `mod N-1` instead of `mod N`; the traced multiplier likewise becomes
    /// `×k mod N-1`. On the diamond this turns the composite ring `Z_2880`
    /// into the *prime field* `Z_2879` — dissolving every gcd obstruction —
    /// at the price of a **seam**: the two representatives of zero
    /// (`0` and `N-1`) make the traced identity the defective operator
    /// `I + |0⟩⟨N-1|`, whose iteration self-stabilizes the double-zero
    /// (see [`crate::stabilize`]).
    pub fn to_mpo_looped(&self, trunc: TruncSpec) -> Mpo {
        let id: Vec<usize> = (0..self.msg_dim).collect();
        self.to_mpo_looped_twisted(&id, trunc)
    }

    /// The looped boundary **twisted** by a message permutation `σ`: the
    /// exiting message is fed back into the entry *through* `σ`,
    /// `M = Σ_m ⟨σ(m)| machine |m⟩` ([`Transducer::to_mpo_looped`] is the
    /// identity twist).
    ///
    /// The twist selects the ring a second time. Two clean cases (both
    /// proved in THEORY.md and pinned by the tests):
    ///
    /// * **σ = identity** — end-around carry: arithmetic mod `N−1`, with a
    ///   *double zero* (`0` and `N−1` both represent it).
    /// * **σ = reversal** (`m ↦ M−1−m`) — **diminished-one arithmetic
    ///   mod `N+1`**, the encoding of Fermat-number-transform hardware:
    ///   chain value `x` represents `v = x+1 ∈ [1, N] ⊂ Z_{N+1}`. The
    ///   twisted adder `+c` realizes `v → v + (c+1) mod N+1` and the
    ///   twisted multiplier `×k` realizes `v → k·v mod N+1` exactly. The
    ///   unrepresentable zero of `Z_{N+1}` appears as an **annihilated
    ///   branch** — a hole, the precise dual of the straight loop's double
    ///   zero (a seam).
    ///
    /// One machine body thus computes in a family of rings selected purely
    /// at the boundary: `Z_N` (open), `Z_{N−1}` (looped), `Z_{N+1}`
    /// (reversal-twisted) — with unitarity governed by `gcd` against
    /// whichever ring the boundary chose.
    ///
    /// Note on truncation: the branch sum is accumulated under `trunc`,
    /// and for large message dimensions the *partial* sums can transiently
    /// need bond dimension near the profile's geometric cap — supply a
    /// rank cap that accommodates it, or study large-message twists
    /// branch-by-branch (`examples/boundary_twists.rs` does the latter
    /// for `×43` on the diamond).
    pub fn to_mpo_looped_twisted(&self, sigma: &[usize], trunc: TruncSpec) -> Mpo {
        assert_eq!(sigma.len(), self.msg_dim, "twist must permute the message set");
        let mut seen = vec![false; self.msg_dim];
        for &v in sigma {
            assert!(v < self.msg_dim && !seen[v], "twist must be a permutation");
            seen[v] = true;
        }
        let mut acc: Option<Mpo> = None;
        for m in 0..self.msg_dim {
            let branch = self.to_mpo_with(Boundary::Fixed(m), Boundary::Fixed(sigma[m]), trunc);
            acc = Some(match acc {
                None => branch,
                Some(a) => a.add(&branch, trunc),
            });
        }
        acc.expect("msg_dim >= 1")
    }

    /// Lift the transducer to an MPO with the standard machine boundaries:
    /// the front enters in state 0 and the exiting state is dropped (for
    /// `mod N` arithmetic, the wrap).
    pub fn to_mpo(&self, trunc: TruncSpec) -> Mpo {
        self.to_mpo_with(Boundary::Fixed(0), Boundary::SumAll, trunc)
    }

    /// Lift the transducer to an MPO with explicit boundary conditions on
    /// the message front. The message legs become virtual bonds; the
    /// carrier is canonicalized under an **exact** policy — cascade
    /// branches can carry arbitrarily small Frobenius weight (a deep carry
    /// has weight `~1/N`) while being free to keep in rank.
    pub fn to_mpo_with(&self, entry: Boundary, exit: Boundary, trunc: TruncSpec) -> Mpo {
        let n = self.profile.len();
        let m = self.msg_dim;
        let (entry_site, exit_site) = match self.direction {
            Direction::RightToLeft => (n - 1, 0),
            Direction::LeftToRight => (0, n - 1),
        };
        let tensors: Vec<SiteTensor> = (0..n)
            .map(|i| {
                let d = self.profile[i];
                // Bond holding message flowing *out* of this site / *into*
                // this site, by direction.
                let (dl, dr) = match self.direction {
                    Direction::RightToLeft => (
                        if i == exit_site { 1 } else { m },
                        if i == entry_site { 1 } else { m },
                    ),
                    Direction::LeftToRight => (
                        if i == entry_site { 1 } else { m },
                        if i == exit_site { 1 } else { m },
                    ),
                };
                let mut t = SiteTensor::zeros(dl, d * d, dr);
                let entering: Vec<usize> = if i == entry_site {
                    match entry {
                        Boundary::Fixed(s) => vec![s],
                        Boundary::SumAll => (0..m).collect(),
                    }
                } else {
                    (0..m).collect()
                };
                for &msg_in in &entering {
                    for p_in in 0..d {
                        let (msg_out, p_out) = self.rules[i][msg_in * d + p_in];
                        if i == exit_site {
                            if let Boundary::Fixed(s) = exit {
                                if msg_out != s {
                                    continue;
                                }
                            }
                        }
                        let bond_in = if i == entry_site { 0 } else { msg_in };
                        let bond_out = if i == exit_site { 0 } else { msg_out };
                        let (l, r) = match self.direction {
                            Direction::RightToLeft => (bond_out, bond_in),
                            Direction::LeftToRight => (bond_in, bond_out),
                        };
                        let cur = t.at(l, p_out * d + p_in, r);
                        t.set(l, p_out * d + p_in, r, cur + C64::ONE);
                    }
                }
                t
            })
            .collect();
        let mut carrier = Mps {
            dims: self.profile.iter().map(|&d| d * d).collect(),
            tensors,
            center: 0,
            trunc: TruncSpec::exact(),
            discarded_weight: 0.0,
        };
        carrier.recompress();
        carrier.trunc = trunc;
        Mpo {
            dims: self.profile.clone(),
            carrier,
        }
    }
}

/// `1 − hs_fidelity(A†∘A, I)`: zero for unitary operators, positive when
/// the cascade's global map fails injectivity (e.g. `gcd(k, N) ≠ 1`).
pub fn unitarity_defect(a: &Mpo, trunc: TruncSpec) -> f64 {
    let gram = a.adjoint().compose_after(a, trunc);
    let id = Mpo::identity(&a.dims, trunc);
    1.0 - gram.hs_fidelity(&id)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dense::total_dim;
    use crate::mps::Mps;
    use crate::radix;
    use crate::zigzag;

    const SPEC: TruncSpec = TruncSpec {
        max_rank: 96,
        cutoff: 1e-12,
    };

    fn basis(profile: &[usize], mut x: u128) -> Mps {
        let mut d = vec![0usize; profile.len()];
        for (i, &dim) in profile.iter().enumerate().rev() {
            d[i] = (x % dim as u128) as usize;
            x /= dim as u128;
        }
        Mps::basis_state(profile, &d, SPEC)
    }

    #[test]
    fn transducer_adder_matches_carry_mpo() {
        for profile in [vec![2usize, 3, 4], zigzag::diamond(2, 4)] {
            let a = Transducer::adder(&profile, 41).to_mpo(SPEC);
            let b = radix::adder_mpo_exact(&profile, 41, SPEC);
            let f = a.hs_fidelity(&b);
            assert!((f - 1.0).abs() < 1e-10, "profile {:?}: {}", profile, f);
            assert!(a.max_bond_dim() <= 2);
        }
    }

    #[test]
    fn mult_cascade_is_modular_multiplication() {
        let profile = vec![2usize, 3, 4]; // N = 24
        let n_total = total_dim(&profile) as u128;
        let k = 5; // gcd(5, 24) = 1
        let m = Transducer::mult(&profile, k).to_mpo(SPEC);
        assert!(m.max_bond_dim() <= k, "bond {} > k {}", m.max_bond_dim(), k);
        for x in 0..n_total {
            let out = m.apply_to(&basis(&profile, x)).to_dense();
            let target = (k as u128 * x % n_total) as usize;
            assert!(
                (out.amps[target].abs() - 1.0).abs() < 1e-9,
                "x={}: weight {}",
                x,
                out.amps[target].abs()
            );
        }
        assert!(unitarity_defect(&m, SPEC) < 1e-8);
    }

    #[test]
    fn mult_cascade_handles_deep_carries_at_scale() {
        // ×7 on the 25-site wave (N ≈ 8.6e12), worst-case input N-1.
        let wave = zigzag::wave(2, 5, 4);
        let big_n: u128 = wave.iter().map(|&d| d as u128).product();
        let m = Transducer::mult(&wave, 7).to_mpo(SPEC);
        for x in [1u128, big_n - 1, big_n / 3] {
            let out = m.apply_to(&basis(&wave, x));
            let expect = basis(&wave, 7 * x % big_n);
            let f = out.fidelity(&expect);
            assert!((f - 1.0).abs() < 1e-9, "x={}: fidelity {}", x, f);
        }
    }

    #[test]
    fn cascades_compose_recursively() {
        let profile = vec![2usize, 3, 4]; // N = 24
        let m5 = Transducer::mult(&profile, 5).to_mpo(SPEC);
        let m7 = Transducer::mult(&profile, 7).to_mpo(SPEC);
        // 5·7 = 35 ≡ 11 (mod 24)
        let m11 = Transducer::mult(&profile, 11).to_mpo(SPEC);
        let composed = m5.compose_after(&m7, SPEC);
        let f = composed.hs_fidelity(&m11);
        assert!((f - 1.0).abs() < 1e-8, "m5∘m7 vs m11: {}", f);
        // Adders and multipliers interleave: x → 5(x+3) = 5x + 15.
        let a3 = Transducer::adder(&profile, 3).to_mpo(SPEC);
        let a15 = Transducer::adder(&profile, 15).to_mpo(SPEC);
        let lhs = m5.compose_after(&a3, SPEC);
        let rhs = a15.compose_after(&m5, SPEC);
        let f2 = lhs.hs_fidelity(&rhs);
        assert!((f2 - 1.0).abs() < 1e-8, "affine identity: {}", f2);
    }

    #[test]
    fn gcd_obstruction_breaks_unitarity() {
        let profile = vec![2usize, 3, 4]; // N = 24
        let m2 = Transducer::mult(&profile, 2).to_mpo(SPEC); // gcd(2,24)=2
        let defect = unitarity_defect(&m2, SPEC);
        assert!(
            defect > 0.1,
            "gcd≠1 should break unitarity: defect {}",
            defect
        );
        // The collision is visible on states: |0⟩ and |12⟩ both map to |0⟩.
        let psi_a = m2.apply_to(&basis(&profile, 0));
        let psi_b = m2.apply_to(&basis(&profile, 12));
        let overlap = psi_a.fidelity(&psi_b);
        assert!(
            (overlap - 1.0).abs() < 1e-9,
            "colliding images: {}",
            overlap
        );
        // A unitary comparison case for contrast.
        let m5 = Transducer::mult(&profile, 5).to_mpo(SPEC);
        assert!(unitarity_defect(&m5, SPEC) < 1e-8);
    }

    #[test]
    fn division_is_the_opposed_inverse_machine() {
        // ×k⁻¹ via the MSB-first remainder machine at width k, checked
        // three ways on Z_72: against the adjoint of ×k, against a
        // directly built ×k⁻¹, and on every basis state.
        let profile = vec![2usize, 3, 4, 3]; // N = 72
        let n_total = total_dim(&profile) as u128;
        let k = 5;
        let k_inv = 29u128; // 5·29 = 145 ≡ 1 (mod 72)
        let dv = Transducer::div_mpo(&profile, k, SPEC);
        assert!(
            dv.max_bond_dim() <= k,
            "division width {} > k",
            dv.max_bond_dim()
        );

        let mk = Transducer::mult(&profile, k).to_mpo(SPEC);
        let f_adj = dv.hs_fidelity(&mk.adjoint());
        assert!((f_adj - 1.0).abs() < 1e-9, "÷5 vs (×5)†: {}", f_adj);

        let mk_inv = Transducer::mult(&profile, k_inv as usize).to_mpo(SPEC);
        let f_dir = dv.hs_fidelity(&mk_inv);
        assert!((f_dir - 1.0).abs() < 1e-9, "÷5 vs ×29: {}", f_dir);

        for x in 0..n_total {
            let out = dv.apply_to(&basis(&profile, x)).to_dense();
            let target = (k_inv * x % n_total) as usize;
            assert!(
                (out.amps[target].abs() - 1.0).abs() < 1e-9,
                "x={}: weight {}",
                x,
                out.amps[target].abs()
            );
        }
    }

    #[test]
    fn opposed_fronts_annihilate() {
        // ÷k ∘ ×k = identity: the right-moving carry front and the
        // left-moving remainder front cancel exactly.
        let profile = zigzag::diamond(2, 4);
        let mk = Transducer::mult(&profile, 7).to_mpo(SPEC);
        let dv = Transducer::div_mpo(&profile, 7, SPEC);
        let id = Mpo::identity(&profile, SPEC);
        let f1 = dv.compose_after(&mk, SPEC).hs_fidelity(&id);
        let f2 = mk.compose_after(&dv, SPEC).hs_fidelity(&id);
        assert!((f1 - 1.0).abs() < 1e-9, "÷7∘×7: {}", f1);
        assert!((f2 - 1.0).abs() < 1e-9, "×7∘÷7: {}", f2);
    }

    #[test]
    fn dual_front_ratio_is_the_modular_quotient() {
        // ×5 ∘ ÷7 = ×(5·7⁻¹ mod 72) = ×11: two opposed narrow fronts
        // realize a multiplier whose single-direction machine would need
        // message width 11 — and for larger rings, far worse.
        let profile = vec![2usize, 3, 4, 3]; // N = 72; 7⁻¹ = 31; 5·31 ≡ 11
        let mk5 = Transducer::mult(&profile, 5).to_mpo(SPEC);
        let dv7 = Transducer::div_mpo(&profile, 7, SPEC);
        let ratio = mk5.compose_after(&dv7, SPEC);
        let mk11 = Transducer::mult(&profile, 11).to_mpo(SPEC);
        let f = ratio.hs_fidelity(&mk11);
        assert!((f - 1.0).abs() < 1e-9, "×5∘÷7 vs ×11: {}", f);
    }

    #[test]
    fn reflection_duality_negation() {
        // The other geometric opposition: ring reflection. x → -x mod N is
        // complement (width 1, digit-local) then +1 (width 2), matching the
        // directly built ×(N-1) at total machine width 2 instead of N-1.
        use crate::circuit::Circuit;
        use crate::gates;
        let profile = vec![2usize, 3, 4]; // N = 24
        let mut comp_circuit = Circuit::new(profile.clone());
        for (i, &d) in profile.iter().enumerate() {
            comp_circuit.one(i, gates::complement(d));
        }
        let comp = Mpo::from_circuit(&comp_circuit, SPEC);
        assert_eq!(comp.max_bond_dim(), 1);
        let plus1 = Transducer::adder(&profile, 1).to_mpo(SPEC);
        let neg = plus1.compose_after(&comp, SPEC);
        let m23 = Transducer::mult(&profile, 23).to_mpo(SPEC);
        let f = neg.hs_fidelity(&m23);
        assert!((f - 1.0).abs() < 1e-9, "(+1)∘complement vs ×23: {}", f);
        assert!(neg.max_bond_dim() <= 2);
    }

    #[test]
    fn looped_adder_is_ones_complement_arithmetic() {
        // Feeding the carry back into the entry (end-around carry) changes
        // the ring: the traced adder computes mod N-1, with a two-branch
        // seam where both representatives of zero appear.
        let profile = vec![2usize, 3, 4]; // N = 24 → traced ring Z_23
        let n_total = total_dim(&profile) as u128;
        let c = 5u128;
        let m = Transducer::adder(&profile, c).to_mpo_looped(SPEC);
        for x in [0u128, 3, 10, 17, 20, 22] {
            if (x + c).is_multiple_of(n_total - 1) && x + c != 0 {
                continue; // seam handled below
            }
            let target = (x + c) % (n_total - 1);
            let out = m.apply_to(&basis(&profile, x)).to_dense();
            assert!(
                (out.amps[target as usize].abs() - 1.0).abs() < 1e-9,
                "x={}: weight {}",
                x,
                out.amps[target as usize].abs()
            );
        }
        // The seam: x + c = N-1 lands on the double zero — BOTH
        // representatives appear, in equal superposition.
        let x_seam = n_total - 1 - c;
        let out = m.apply_to(&basis(&profile, x_seam)).to_dense();
        let a_hi = out.amps[(n_total - 1) as usize].abs();
        let a_lo = out.amps[0].abs();
        assert!(
            (a_hi - 1.0).abs() < 1e-9 && (a_lo - 1.0).abs() < 1e-9,
            "seam branches: |N-1⟩ {} |0⟩ {}",
            a_hi,
            a_lo
        );
    }

    #[test]
    fn traced_identity_is_the_seam_jordan_block() {
        // The looped +0 machine is exactly I + |0⟩⟨N-1|.
        let profile = vec![2usize, 3, 4];
        let m = Transducer::adder(&profile, 0).to_mpo_looped(SPEC);
        let seam_digits: Vec<usize> = profile.iter().map(|&d| d - 1).collect();
        let expect = Mpo::identity(&profile, SPEC).add(
            &Mpo::basis_transfer(&profile, &[0, 0, 0], &seam_digits, SPEC),
            SPEC,
        );
        let f = m.hs_fidelity(&expect);
        assert!((f - 1.0).abs() < 1e-10, "I + |0⟩⟨N-1|: {}", f);
    }

    #[test]
    fn looped_boundary_heals_the_gcd_obstruction() {
        // ×6 on Z_2880 is 6-to-1 (defect 5/6). Close the loop: the traced
        // multiplier acts on Z_2879 — a prime field — where gcd(6, 2879)=1
        // and the machine is exactly unitary again.
        let profile = zigzag::diamond(2, 5); // N = 2880, N-1 = 2879 prime
        let n_total = total_dim(&profile) as u128;
        let open = Transducer::mult(&profile, 6).to_mpo(SPEC);
        let looped = Transducer::mult(&profile, 6).to_mpo_looped(SPEC);
        let d_open = unitarity_defect(&open, SPEC);
        let d_loop = unitarity_defect(&looped, SPEC);
        assert!(d_open > 0.8, "open defect {}", d_open);
        assert!(d_loop < 1e-8, "looped defect {}", d_loop);
        // And it really is ×6 mod 2879 (0 and N-1 are its fixed double-zero).
        for x in [1u128, 500, 2879] {
            let target = if x == n_total - 1 {
                n_total - 1
            } else {
                6 * x % (n_total - 1)
            };
            let f = looped
                .apply_to(&basis(&profile, x))
                .fidelity(&basis(&profile, target));
            assert!((f - 1.0).abs() < 1e-8, "x={}", x);
        }
    }

    #[test]
    fn mps_add_matches_dense_addition() {
        use crate::dense::DenseState;
        use crate::mat::Rng;
        let dims = vec![2usize, 3, 4];
        let mut rng = Rng::new(77);
        let mut a = DenseState::zero_state(&dims);
        let mut b = DenseState::zero_state(&dims);
        for v in &mut a.amps {
            *v = rng.c_gaussian();
        }
        for v in &mut b.amps {
            *v = rng.c_gaussian();
        }
        let ma = Mps::from_dense(&a, SPEC);
        let mb = Mps::from_dense(&b, SPEC);
        let sum = ma.add(&mb, SPEC).to_dense();
        let diff: f64 = sum
            .amps
            .iter()
            .zip(a.amps.iter().zip(b.amps.iter()))
            .map(|(s, (x, y))| (*s - (*x + *y)).abs())
            .fold(0.0, f64::max);
        assert!(diff < 1e-10, "max diff {}", diff);
    }

    #[test]
    fn reversal_twisted_adder_is_diminished_one_mod_n_plus_one() {
        // Twisting the carry loop by the swap turns the +c machine into
        // the diminished-one adder of Fermat-transform hardware:
        // value v = x+1 ∈ [1, N] represents Z_{N+1}, the map is
        // v → v + (c+1) mod N+1, and the unrepresentable zero of Z_{N+1}
        // is an annihilated branch (a hole) at x = N−1−c.
        let profile = vec![2usize, 3, 4]; // N = 24, twist ring Z_25
        let n = total_dim(&profile) as u128;
        let c = 5u128;
        let m = Transducer::adder(&profile, c).to_mpo_looped_twisted(&[1, 0], SPEC);
        for x in [0u128, 3, 17, 20, 23] {
            let v_out = (x + 1 + c + 1) % (n + 1);
            assert_ne!(v_out, 0, "sample hit the hole");
            let out = m.apply_to(&basis(&profile, x));
            let f = out.fidelity(&basis(&profile, v_out - 1));
            assert!((f - 1.0).abs() < 1e-9, "x={}: fidelity {}", x, f);
        }
        // The hole: x = N−1−c lands on the missing zero and is annihilated.
        let hole = m.apply_to(&basis(&profile, n - 1 - c));
        assert!(hole.norm() < 1e-9, "hole norm {}", hole.norm());
    }

    #[test]
    fn reversal_twisted_multiplier_computes_mod_n_plus_one() {
        // N = 72 sits between the twin primes 71 and 73: BOTH closed
        // boundary rings are fields. ×6 (gcd(6, 72) = 6, defect 5/6 open)
        // is healed by the straight loop (mod 71) AND by the reversal
        // twist (mod 73) — where it is exactly v → 6v mod 73 on the
        // diminished-one encoding, hole-free since gcd(6, 73) = 1.
        let profile = vec![2usize, 3, 4, 3]; // N = 72
        let n = total_dim(&profile) as u128;
        let k = 6usize;
        let sigma: Vec<usize> = (0..k).rev().collect();
        let open = Transducer::mult(&profile, k).to_mpo(SPEC);
        let looped = Transducer::mult(&profile, k).to_mpo_looped(SPEC);
        let twisted = Transducer::mult(&profile, k).to_mpo_looped_twisted(&sigma, SPEC);
        assert!(unitarity_defect(&open, SPEC) > 0.8);
        assert!(unitarity_defect(&looped, SPEC) < 1e-8);
        assert!(unitarity_defect(&twisted, SPEC) < 1e-8);
        for x in [0u128, 35, 50, 71] {
            let v_out = (k as u128 * (x + 1)) % (n + 1);
            let out = twisted.apply_to(&basis(&profile, x));
            let f = out.fidelity(&basis(&profile, v_out - 1));
            assert!((f - 1.0).abs() < 1e-8, "x={}: fidelity {}", x, f);
        }
    }

    #[test]
    fn twist_ring_gcd_governs_unitarity_both_ways() {
        // The twist can BREAK as well as heal: ×5 on Z_24 is unitary open
        // (gcd(5, 24) = 1) but gcd(5, 25) = 5 on the twist ring — holes
        // appear at v ≡ 0 (mod 5) and the twisted operator is defective.
        let profile = vec![2usize, 3, 4]; // N = 24, twist ring Z_25 = 5²
        let m5_open = Transducer::mult(&profile, 5).to_mpo(SPEC);
        let m5_tw = Transducer::mult(&profile, 5).to_mpo_looped_twisted(&[4, 3, 2, 1, 0], SPEC);
        assert!(unitarity_defect(&m5_open, SPEC) < 1e-8);
        assert!(unitarity_defect(&m5_tw, SPEC) > 0.1);
        // v = 5 (x = 4): 5·5 = 25 ≡ 0 (mod 25) — annihilated.
        let hole = m5_tw.apply_to(&basis(&profile, 4));
        assert!(hole.norm() < 1e-9, "hole norm {}", hole.norm());
    }

    #[test]
    fn twisted_traced_identity_is_the_unilateral_shift() {
        // The swap-twisted +0 machine is the unilateral shift
        // Σ_{x<N−1} |x+1⟩⟨x| = A_1 − |0⟩⟨N−1|: everything conveys toward
        // the missing zero and drains — the dual of the straight-traced
        // identity I + |0⟩⟨N−1|, which pumps the double zero.
        use crate::c64::C64;
        let profile = vec![2usize, 3, 4]; // N = 24
        let m = Transducer::adder(&profile, 0).to_mpo_looped_twisted(&[1, 0], SPEC);
        let a1 = Transducer::adder(&profile, 1).to_mpo(SPEC);
        let seam_digits: Vec<usize> = profile.iter().map(|&d| d - 1).collect();
        let expected = a1.add(
            &Mpo::basis_transfer(&profile, &[0, 0, 0], &seam_digits, SPEC)
                .scale(C64::real(-1.0)),
            SPEC,
        );
        let f = m.hs_fidelity(&expected);
        assert!((f - 1.0).abs() < 1e-10, "A_1 − |0⟩⟨N−1|: {}", f);
        assert!(m.apply_to(&basis(&profile, 23)).norm() < 1e-9);
    }

    #[test]
    fn fourier_conjugation_inverts_the_multiplier_and_reverses_the_cascade() {
        // The scaling theorem F ∘ M_k ∘ F† = M_{k⁻¹ mod N}, transported
        // through the reversed-digit frame: conjugating by V inverts the
        // multiplier AND flips the carry direction (the reversed encoding
        // has its least significant digit at site 0).
        let profile = vec![2usize, 3, 4, 3]; // N = 72, non-palindromic
        let k = 5;
        let k_inv = 29; // 5·29 = 145 = 2·72 + 1
        let v = Mpo::from_circuit(&radix::mixed_radix_qft_reversed(&profile, 0.0), SPEC);
        let mk = Transducer::mult(&profile, k).to_mpo(SPEC);
        let conj = v.compose_after(&mk.compose_after(&v.adjoint(), SPEC), SPEC);
        let mk_inv_rev =
            Transducer::mult_directed(&profile, k_inv, Direction::LeftToRight).to_mpo(SPEC);
        let f = conj.hs_fidelity(&mk_inv_rev);
        assert!((f - 1.0).abs() < 1e-7, "V∘M5∘V† vs reversed M29: {}", f);
    }
}
