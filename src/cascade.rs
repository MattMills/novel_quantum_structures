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
