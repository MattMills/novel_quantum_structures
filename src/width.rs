//! The **a-priori width calculus** of streamed permutations: the operator
//! Schmidt rank of modular multiplication across a chain cut, computed from
//! pure number theory — no tensors involved.
//!
//! For a cut splitting the chain into a left ring of size `L` and a right
//! block of size `S` (`N = L·S`), the multiplier `M_k : |x⟩ → |k·x mod N⟩`
//! factors branch-wise over the *crossing message* — the multiplication
//! carry entering the cut:
//!
//! ```text
//!   x = a·S + b   ⇒   M_k |a, b⟩ = |(k·a + c(b)) mod L,  k·b mod S⟩,
//!   c(b) = ⌊k·b / S⌋,
//! ```
//!
//! and the operator Schmidt rank across the cut is **exactly**
//! `|{c(b) mod L : b ∈ [0, S)}|` (the cut-rank theorem, THEORY.md §8).
//! The bond dimension of the exactly-recompressed cascade MPO is therefore
//! a number computable *in advance*; the tests pin the two against each
//! other bond-for-bond, including through operator composition.
//!
//! Two regimes:
//!
//! * `k ≤ S` — `b ↦ ⌊kb/S⌋` steps by 0 or 1 and reaches `k−1`, so it is a
//!   surjection onto `[0, k)`: the width is `min(k, L)` in closed form.
//! * `k > S` — the sequence skips and number theory takes over (this is
//!   where the resonances of README finding 16 live: `×2401` thinner than
//!   `×7`); the width is found by enumerating the `S` crossing values.
//!
//! `k` may always be reduced mod `N`: replacing `k` by `k + LS` changes
//! `c(b)` by `L·b ≡ 0 (mod L)`.
//!
//! The counting theorem is a property of *ring-native* multiplication —
//! the modulus and the chain agreeing on `N`. For any other permutation of
//! `Z_N` — foreign-modulus multiplication `×k mod M` with `M ≠ N`
//! ([`foreign_mult`]) chief among them — no closed form is known, and
//! [`perm_cut_rank`] fills the gap: the exact operator Schmidt rank of an
//! *arbitrary* permutation across a cut, computed from a Gram matrix on
//! the smaller side of the cut, still with no tensors. This is the
//! instrument behind the **reduction penalty** measurements
//! (`examples/foreign_rings.rs`, THEORY.md §8.13): native multiplication
//! costs its resonant width `w`, while foreign reduction costs `≈ (k+1)²`
//! capped by the squared geometry — and at the boundary rings `N ∓ 1` it
//! costs exactly the *square of the minimal opposed-front width*
//! `(a + b − 1)²` over representations `k ≡ a·b⁻¹`.

use crate::c64::C64;
use crate::mat::{svd, Mat};
use std::collections::{HashMap, HashSet};

/// Largest right-block size the enumerating branch will scan. Mid-chain
/// cuts of the 25-site wave (`S ≈ 4.1·10⁶`) stay inside, so every bond of
/// the wave is exactly computable from one side or the other.
pub const ENUM_LIMIT: u128 = 8_000_000;

/// A per-cut width statement: exact rank, or only the a-priori bound.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Width {
    /// The exact operator Schmidt rank across the cut.
    Exact(u128),
    /// Only the bound `min(k mod N, L, S)`: the crossing set was too large
    /// to enumerate and no closed form applies (`k > S > ENUM_LIMIT`).
    UpperBound(u128),
}

impl Width {
    pub fn value(self) -> u128 {
        match self {
            Width::Exact(v) | Width::UpperBound(v) => v,
        }
    }

    pub fn is_exact(self) -> bool {
        matches!(self, Width::Exact(_))
    }
}

impl std::fmt::Display for Width {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Width::Exact(v) => write!(f, "{}", v),
            Width::UpperBound(v) => write!(f, "≤{}", v),
        }
    }
}

/// Operator width of `×k` across a cut with left ring `l` and right block
/// `s`: the count `|{⌊k·b/s⌋ mod l : b ∈ [0, s)}|`. Equals the operator
/// Schmidt rank of the multiplier across the cut when `gcd(k, l·s) = 1`
/// (the count is well-defined regardless).
pub fn mult_cut_width(k: u128, l: u128, s: u128) -> Width {
    assert!(l >= 1 && s >= 1, "cut sizes must be positive");
    let k = k % (l * s);
    if k == 0 {
        return Width::Exact(1);
    }
    if k <= s {
        return Width::Exact(k.min(l));
    }
    if s <= ENUM_LIMIT {
        // Count distinct residues with a bitset when the left ring is
        // small enough; fall back to hashing for huge L (where the count
        // is bounded by s anyway).
        const BITSET_LIMIT: u128 = 1 << 25;
        let count = if l <= BITSET_LIMIT {
            let mut seen = vec![false; l as usize];
            let mut count = 0u128;
            for b in 0..s {
                let c = ((k * b / s) % l) as usize;
                if !seen[c] {
                    seen[c] = true;
                    count += 1;
                }
            }
            count
        } else {
            let mut seen: HashSet<u128> = HashSet::new();
            for b in 0..s {
                seen.insert((k * b / s) % l);
            }
            seen.len() as u128
        };
        return Width::Exact(count);
    }
    Width::UpperBound(k.min(l).min(s))
}

/// Per-bond width profile of `×k` on a dimension profile: entry `m` is the
/// width across the bond between sites `m` and `m+1` (left ring
/// `Π_{i≤m} d_i`, right block `Π_{i>m} d_i`).
pub fn mult_width_profile(profile: &[usize], k: u128) -> Vec<Width> {
    let n: u128 = profile.iter().map(|&d| d as u128).product();
    let mut left = 1u128;
    profile[..profile.len() - 1]
        .iter()
        .map(|&d| {
            left *= d as u128;
            mult_cut_width(k, left, n / left)
        })
        .collect()
}

/// Exact operator Schmidt rank of an **arbitrary permutation** of
/// `Z_{l·s}` across the cut with left ring `l` and right block `s` — the
/// rank of the reshuffled 0/1 matrix `R[(a',a),(b',b)] = [π(a·s+b) =
/// a'·s+b']`, whose closed form (`mult_cut_width`) exists only for
/// ring-native multiplication.
///
/// Computed with no tensors: the `N = l·s` ones of `R` are grouped by
/// their key on the larger side of the cut, the Gram matrix `G = R·Rᵀ` is
/// accumulated over the smaller side's distinct keys (`≤ min(N, l², s²)`
/// of them), and the rank is read from `G`'s spectrum (`G` is an integer
/// PSD matrix, so its null space is exact and the spectral gap at zero is
/// clean). Cost: `O(N)` to bin, `O(Σ group²) ≤ O(N·min(l,s))` to
/// accumulate, one SVD of the small Gram. Feasible wherever `N` itself is
/// enumerable — this is the instrument for permutations *outside* the
/// native arithmetic family, where no a-priori counting formula is known.
/// The tests pin it against `mult_width_profile` on the native family and
/// against exactly-recompressed dense-Choi MPOs off it.
pub fn perm_cut_rank<F: Fn(usize) -> usize>(l: usize, s: usize, perm: &F) -> usize {
    assert!(l >= 1 && s >= 1, "cut sizes must be positive");
    let n = l * s;
    let mut groups: HashMap<u64, Vec<u32>> = HashMap::new();
    let mut small_ids: HashMap<u64, u32> = HashMap::new();
    let mut hit = vec![false; n];
    for x in 0..n {
        let y = perm(x);
        assert!(y < n && !hit[y], "not a permutation of Z_{{l·s}}");
        hit[y] = true;
        let (a, b) = (x / s, x % s);
        let (ap, bp) = (y / s, y % s);
        let (small_key, large_key) = if l <= s {
            ((ap * l + a) as u64, (bp * s + b) as u64)
        } else {
            ((bp * s + b) as u64, (ap * l + a) as u64)
        };
        let next = small_ids.len() as u32;
        let id = *small_ids.entry(small_key).or_insert(next);
        groups.entry(large_key).or_default().push(id);
    }
    let m = small_ids.len();
    let mut g = Mat::zeros(m, m);
    for ids in groups.values() {
        for &r1 in ids {
            for &r2 in ids {
                let cur = g.at(r1 as usize, r2 as usize);
                g.set(r1 as usize, r2 as usize, cur + C64::ONE);
            }
        }
    }
    let sv = svd(&g);
    if sv.s.is_empty() || sv.s[0] <= 0.0 {
        return 0;
    }
    let tol = sv.s[0] * (m as f64) * 1e-12;
    sv.s.iter().filter(|&&v| v > tol).count()
}

/// Per-bond operator Schmidt ranks of a permutation on a dimension
/// profile — [`perm_cut_rank`] across every bond.
pub fn perm_width_profile<F: Fn(usize) -> usize>(profile: &[usize], perm: &F) -> Vec<usize> {
    let n: usize = profile.iter().product();
    let mut left = 1usize;
    profile[..profile.len() - 1]
        .iter()
        .map(|&d| {
            left *= d;
            perm_cut_rank(left, n / left, perm)
        })
        .collect()
}

/// The **foreign multiplier**: `x → k·x mod m` on `[0, m)`, identity on
/// the tail `[m, n)` — modular multiplication in a ring the chain does
/// *not* natively carry, extended to a permutation of the chain's own
/// `Z_n`. Requires `gcd(k, m) = 1` (else it is not a permutation) and
/// `m ≤ n`. This is the object whose cut ranks measure the *reduction
/// penalty* of foreign-modulus arithmetic (`examples/foreign_rings.rs`).
pub fn foreign_mult(n: usize, k: usize, m: usize) -> impl Fn(usize) -> usize {
    assert!(m >= 1 && m <= n, "foreign modulus must satisfy 1 ≤ m ≤ n");
    assert_eq!(gcd(k as u128, m as u128), 1, "gcd(k, m) must be 1");
    move |x| if x < m { k * x % m } else { x }
}

fn gcd(a: u128, b: u128) -> u128 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cascade::Transducer;
    use crate::mat::TruncSpec;
    use crate::zigzag;

    /// Direct enumeration, used to pin the closed-form branch.
    fn enumerated(k: u128, l: u128, s: u128) -> u128 {
        let k = k % (l * s);
        let mut seen: HashSet<u128> = HashSet::new();
        for b in 0..s {
            seen.insert((k * b / s) % l);
        }
        seen.len() as u128
    }

    fn mod_inverse(k: u128, n: u128) -> u128 {
        let (mut old_r, mut r) = (k as i128, n as i128);
        let (mut old_s, mut s) = (1i128, 0i128);
        while r != 0 {
            let q = old_r / r;
            (old_r, r) = (r, old_r - q * r);
            (old_s, s) = (s, old_s - q * s);
        }
        assert_eq!(old_r, 1, "gcd(k, n) must be 1");
        old_s.rem_euclid(n as i128) as u128
    }

    #[test]
    fn closed_form_matches_enumeration_when_k_is_small() {
        for l in [2u128, 6, 24, 120] {
            for s in [6u128, 24, 120] {
                for k in [1u128, 2, 3, 5, s - 1, s] {
                    assert_eq!(
                        mult_cut_width(k, l, s).value(),
                        enumerated(k, l, s),
                        "k={} l={} s={}",
                        k,
                        l,
                        s
                    );
                }
            }
        }
    }

    #[test]
    fn known_resonances_at_the_diamond_waist() {
        // (L, S) = (24, 120): the waist cut of diamond(2, 5), README
        // finding 16 measured widths 7 → 24 → 6 → 3 along the squaring
        // orbit of 7, and 7 for the inverse multiplier 823.
        for (k, expect) in [(7u128, 7u128), (49, 24), (2401, 6), (1921, 3), (823, 7)] {
            let w = mult_cut_width(k, 24, 120);
            assert!(w.is_exact());
            assert_eq!(w.value(), expect, "k={}", k);
        }
    }

    #[test]
    fn width_is_symmetric_under_inversion() {
        // rank(A) = rank(A†) and M_k† = M_{k⁻¹}: the count must agree.
        let n = 2880u128;
        for k in [7u128, 11, 49, 121, 823, 2401] {
            let ki = mod_inverse(k, n);
            assert_eq!(
                mult_cut_width(k, 24, 120).value(),
                mult_cut_width(ki, 24, 120).value(),
                "k={} k⁻¹={}",
                k,
                ki
            );
        }
    }

    #[test]
    fn atlas_matches_recompressed_cascade_bonds() {
        // The theorem, executed: exactly-recompressed cascade MPO bond
        // dimensions equal the number-theoretic prediction, bond for bond.
        let cases: Vec<(Vec<usize>, Vec<usize>)> = vec![
            (vec![2, 3, 4], vec![5, 7, 11]),
            (vec![2, 3, 4, 3], vec![5, 7, 29]),
            (zigzag::diamond(2, 4), vec![7, 11, 49]),
        ];
        for (profile, ks) in cases {
            for k in ks {
                let m = Transducer::mult(&profile, k).to_mpo(TruncSpec::exact());
                let predicted: Vec<u128> = mult_width_profile(&profile, k as u128)
                    .into_iter()
                    .map(|w| {
                        assert!(w.is_exact());
                        w.value()
                    })
                    .collect();
                let measured: Vec<u128> =
                    m.bond_dims().iter().map(|&b| b as u128).collect();
                assert_eq!(measured, predicted, "profile {:?} k {}", profile, k);
            }
        }
    }

    #[test]
    fn composition_realizes_the_predicted_resonance() {
        // ×49 = ×7 ∘ ×7 on diamond(2, 4): the atlas predicts [2,3,3,2] —
        // narrower than ×7 itself ([2,6,6,2]) and far below the message
        // bound 49. Exact recompression of the composed operator must land
        // exactly there, and match the directly built ×49.
        let profile = zigzag::diamond(2, 4);
        let exact = TruncSpec::exact();
        let m7 = Transducer::mult(&profile, 7).to_mpo(exact);
        assert_eq!(m7.bond_dims(), vec![2, 6, 6, 2]);
        let m49 = m7.compose_after(&m7, exact);
        let predicted: Vec<usize> = mult_width_profile(&profile, 49)
            .into_iter()
            .map(|w| w.value() as usize)
            .collect();
        assert_eq!(predicted, vec![2, 3, 3, 2]);
        assert_eq!(m49.bond_dims(), predicted);
        let direct = Transducer::mult(&profile, 49).to_mpo(exact);
        assert!((m49.hs_fidelity(&direct) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn reduction_mod_n_is_sound() {
        // c(b) changes by L·b ≡ 0 (mod L) when k gains a multiple of N.
        for k in [5u128, 49, 823] {
            assert_eq!(
                mult_cut_width(k, 24, 120).value(),
                mult_cut_width(k + 2880, 24, 120).value()
            );
        }
    }

    #[test]
    fn perm_rank_reproduces_the_native_counting_theorem() {
        // On ring-native multiplication the Gram-rank instrument must
        // agree with Theorem 8.5's count, bond for bond.
        let cases: Vec<(Vec<usize>, Vec<usize>)> = vec![
            (vec![2, 3, 4], vec![5, 7, 11]),
            (vec![2, 3, 4, 3], vec![5, 7, 29]),
            (zigzag::diamond(2, 4), vec![7, 11, 49]),
        ];
        for (profile, ks) in cases {
            let n: usize = profile.iter().product();
            for k in ks {
                let counted: Vec<usize> = mult_width_profile(&profile, k as u128)
                    .into_iter()
                    .map(|w| w.value() as usize)
                    .collect();
                let grammed = perm_width_profile(&profile, &foreign_mult(n, k, n));
                assert_eq!(grammed, counted, "profile {:?} k {}", profile, k);
            }
        }
    }

    #[test]
    fn perm_rank_matches_dense_choi_bonds_off_the_native_family() {
        // Off the native family no counting theorem exists; the Gram rank
        // must still equal the bond profile of the exactly-compressed
        // dense-Choi MPO — the definition of operator Schmidt rank.
        use crate::mpo::Mpo;
        let profile = vec![2usize, 3, 4, 3]; // N = 72
        let n: usize = profile.iter().product();
        for (k, m) in [(7usize, 71usize), (7, 69), (5, 61), (7, 72)] {
            let grammed = perm_width_profile(&profile, &foreign_mult(n, k, m));
            let mpo = Mpo::from_permutation(&profile, foreign_mult(n, k, m), TruncSpec::exact());
            assert_eq!(grammed, mpo.bond_dims(), "×{} mod {}", k, m);
        }
    }

    #[test]
    fn small_foreign_rings_are_screened_by_the_cut() {
        // A modulus that fits inside the right block never crosses the
        // cut: the operator is (perm inside the block) ⊕ identity, rank
        // 1 + 1 regardless of k. Foreign cost comes from straddling.
        let r = perm_cut_rank(24, 120, &foreign_mult(2880, 7, 97));
        assert_eq!(r, 2);
    }
}
