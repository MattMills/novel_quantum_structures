//! Twisted zigzag qudit structures: dimension-wave profiles and the
//! cross-scale circuit families that live on them.
//!
//! The basic object is a *dimension profile* — the list of local dimensions
//! along the chain. A [`diamond`] expands from a low dimension up to a waist
//! and collapses back (`2,3,4,5,4,3,2`); a [`wave`] repeats the diamond,
//! sharing the low-dimensional valleys:
//!
//! ```text
//! 5        ●                    ●
//! 4      ●   ●                ●   ●
//! 3    ●       ●            ●       ●
//! 2  ●           ●        ●           ●
//!    ·  ·  ·  ·  ·  (wave continues)  ·
//! ```
//!
//! Palindromic profiles admit **mirror pairs** `(i, n-1-i)` — equal-dimension
//! sites on opposite flanks of the wave. Because their dimensions match,
//! they support the full [`crate::gates::xswap`] exchange and generalized
//! Bell pairing: these are the "bidirectional pairs" that couple the
//! structure across scales without touching the intermediate dimensions.
//! The [`bowtie`] circuit Bell-pairs every mirror pair, producing a state
//! that is *maximally entangled across the waist* yet has tiny bond
//! dimension — the concrete form of the state-space compression this crate
//! explores.

use crate::circuit::Circuit;
use crate::gates;
use crate::mat::Rng;

/// Expand-and-collapse profile `lo, lo+1, …, hi, …, lo+1, lo`.
pub fn diamond(lo: usize, hi: usize) -> Vec<usize> {
    assert!(lo >= 2 && hi >= lo);
    let mut p: Vec<usize> = (lo..=hi).collect();
    p.extend((lo..hi).rev());
    p
}

/// A repeated diamond sharing its valleys: `periods` copies of
/// `lo..hi..lo` glued at the `lo` sites. Length `2·periods·(hi-lo) + 1`.
pub fn wave(lo: usize, hi: usize, periods: usize) -> Vec<usize> {
    assert!(periods >= 1);
    let d = diamond(lo, hi);
    let mut p = d.clone();
    for _ in 1..periods {
        p.extend(d.iter().skip(1));
    }
    p
}

/// Mirror pairs `(i, n-1-i)` of a profile (excluding a middle fixed point).
/// For palindromic profiles (diamonds and waves) the paired sites have equal
/// dimension.
pub fn mirror_pairs(profile: &[usize]) -> Vec<(usize, usize)> {
    let n = profile.len();
    (0..n / 2).map(|i| (i, n - 1 - i)).collect()
}

/// Sites at the maximum dimension of the profile (the waists).
pub fn waist_positions(profile: &[usize]) -> Vec<usize> {
    let hi = *profile.iter().max().unwrap();
    profile
        .iter()
        .enumerate()
        .filter(|(_, &d)| d == hi)
        .map(|(i, _)| i)
        .collect()
}

/// Sites at the minimum dimension of the profile (the valleys — the binary
/// edges when `lo == 2`).
pub fn valley_positions(profile: &[usize]) -> Vec<usize> {
    let lo = *profile.iter().min().unwrap();
    profile
        .iter()
        .enumerate()
        .filter(|(_, &d)| d == lo)
        .map(|(i, _)| i)
        .collect()
}

/// Total dense Hilbert-space dimension of a profile, as `f64` (chains beyond
/// a few sites overflow `usize` quickly — that is rather the point).
pub fn dense_dimension(profile: &[usize]) -> f64 {
    profile.iter().map(|&d| d as f64).product()
}

/// ASCII sketch of the dimension wave.
pub fn ascii_profile(profile: &[usize]) -> String {
    let hi = *profile.iter().max().unwrap();
    let lo = *profile.iter().min().unwrap();
    let mut out = String::new();
    for level in (lo..=hi).rev() {
        out.push_str(&format!("{:>3} ", level));
        for &d in profile {
            out.push_str(if d == level { " ● " } else { "   " });
        }
        out.push('\n');
    }
    out.push_str("    ");
    for i in 0..profile.len() {
        out.push_str(&format!("{:^3}", i % 100));
    }
    out.push('\n');
    out
}

// ---------------------------------------------------------------------------
// Named circuit families
// ---------------------------------------------------------------------------

/// Bell-pair every mirror pair of a palindromic profile:
/// `fourier(d)` on the left partner, then a controlled shift onto the right
/// partner, yielding the generalized Bell state `Σ_j |j,j⟩/√d` per pair.
/// The result is the "bowtie" state — maximal entanglement across the waist
/// at minimal bond dimension.
pub fn bowtie(profile: &[usize]) -> Circuit {
    let mut c = Circuit::new(profile.to_vec());
    for (i, j) in mirror_pairs(profile) {
        assert_eq!(profile[i], profile[j], "bowtie needs a palindromic profile");
        let d = profile[i];
        c.one(i, gates::fourier(d));
        c.two(i, j, gates::cshift(d, d));
    }
    c
}

/// `layers` of brickwork Haar-random nearest-neighbour gates (even bonds,
/// then odd bonds, per layer).
pub fn brickwork_random(profile: &[usize], layers: usize, rng: &mut Rng) -> Circuit {
    let mut c = Circuit::new(profile.to_vec());
    let n = profile.len();
    for _ in 0..layers {
        for parity in 0..2 {
            let mut i = parity;
            while i + 1 < n {
                c.two(i, i + 1, gates::random(profile[i] * profile[i + 1], rng));
                i += 2;
            }
        }
    }
    c
}

/// One deterministic cross-scale entangling round:
/// 1. `fourier` on every site (spread each qudit over all its levels),
/// 2. nearest-neighbour controlled shifts left→right (a wavefront running up
///    and down the dimension wave),
/// 3. controlled phases on every mirror pair (the cross couplings).
pub fn crosswave_round(profile: &[usize]) -> Circuit {
    let n = profile.len();
    let mut c = Circuit::new(profile.to_vec());
    for (i, &d) in profile.iter().enumerate() {
        c.one(i, gates::fourier(d));
    }
    for i in 0..n - 1 {
        c.two(i, i + 1, gates::cshift(profile[i], profile[i + 1]));
    }
    for (i, j) in mirror_pairs(profile) {
        c.two(i, j, gates::cphase(profile[i], profile[j]));
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diamond_and_wave_shapes() {
        assert_eq!(diamond(2, 5), vec![2, 3, 4, 5, 4, 3, 2]);
        assert_eq!(diamond(2, 2), vec![2]);
        let w = wave(2, 5, 3);
        assert_eq!(w.len(), 3 * 6 + 1);
        assert_eq!(w[0], 2);
        assert_eq!(w[w.len() - 1], 2);
        // Palindromic:
        let rev: Vec<usize> = w.iter().rev().cloned().collect();
        assert_eq!(w, rev);
        assert_eq!(waist_positions(&w), vec![3, 9, 15]);
        assert_eq!(valley_positions(&w), vec![0, 6, 12, 18]);
    }

    #[test]
    fn mirror_pairs_have_equal_dims_on_palindromes() {
        for profile in [diamond(2, 5), wave(2, 4, 2), wave(3, 6, 2)] {
            for (i, j) in mirror_pairs(&profile) {
                assert_eq!(profile[i], profile[j]);
            }
        }
    }

    #[test]
    fn dense_dimension_of_diamond() {
        assert_eq!(dense_dimension(&diamond(2, 5)), 2880.0);
    }

    #[test]
    fn circuits_build() {
        let p = diamond(2, 5);
        assert_eq!(bowtie(&p).len(), 6); // 3 pairs × (fourier + cshift)
        let mut rng = Rng::new(1);
        assert!(!brickwork_random(&p, 2, &mut rng).is_empty());
        assert!(!crosswave_round(&p).is_empty());
    }
}
