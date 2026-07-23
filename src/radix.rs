//! Mixed-radix Fourier machinery: the **structured tail-radix phase web**.
//!
//! A qudit chain with dimensions `d_0, …, d_{n-1}` (site 0 most significant)
//! encodes the ring `Z_N`, `N = Π d_k`. The Fourier transform over `Z_N`
//! factorizes over the chain with a remarkably rigid structure. Writing the
//! output digit of site `j` in *reversed place value* (place `Q_j = Π_{k<j}
//! d_k`), the kernel splits as
//!
//! ```text
//!   exp(2πi·x·y/N) = Π_j Π_{i≥j} exp(2πi · x_i · y'_j / D_{j,i}),
//!   D_{j,i} = d_j · d_{j+1} ··· d_i          (the radix segment product)
//! ```
//!
//! — every coupling with `i < j` carries an integer phase and vanishes.
//! So output digit `j` is phase-correlated **only with the tail of the
//! register** (sites `i ≥ j`), through [`crate::gates::cp_frac`] gates whose
//! angles `2π/D_{j,i}` are set by the radix segment between the two sites.
//! Because `D_{j,i}` grows as the product of every dimension in between, the
//! couplings decay super-exponentially with distance: the web is effectively
//! banded, which is exactly why these operators compress well as MPOs.
//!
//! The digit reversal that finishes the transform maps site `i` to site
//! `n-1-i` — on palindromic profiles (diamonds, waves) this is the **mirror
//! pair network**, implemented by `xswap` gates: the zigzag's bidirectional
//! pairs are load-bearing in the Fourier transform over its own ring.
//!
//! On top of the QFT sits Draper-style Fourier arithmetic:
//! [`fourier_phase_ramp`] is the diagonal block `D_c = diag(ω^{c·y})`
//! (a bond-dimension-1 product of single-site gates), and
//! `QFT† ∘ D_c ∘ QFT = |x⟩ → |x + c mod N⟩` — three circuit blocks whose
//! composition collapses to a permutation. See
//! `examples/circuit_of_circuits.rs`.

use crate::c64::C64;
use crate::circuit::Circuit;
use crate::gates;
use crate::mat::TruncSpec;
use crate::mpo::Mpo;
use crate::mps::{Mps, SiteTensor};
use std::f64::consts::PI;

/// Is the profile its own mirror image?
pub fn is_palindromic(profile: &[usize]) -> bool {
    let n = profile.len();
    (0..n / 2).all(|i| profile[i] == profile[n - 1 - i])
}

/// The mixed-radix QFT over `Z_N` with output digits in **reversed place
/// value** (site `j` ends holding the digit of place `Q_j = Π_{k<j} d_k`).
/// Works for any profile. Couplings with angle below `min_angle` radians are
/// dropped (pass 0.0 for the exact transform).
pub fn mixed_radix_qft_reversed(profile: &[usize], min_angle: f64) -> Circuit {
    let n = profile.len();
    let mut c = Circuit::new(profile.to_vec());
    for j in 0..n {
        c.one(j, gates::fourier(profile[j]));
        let mut seg = profile[j] as f64;
        for i in (j + 1)..n {
            seg *= profile[i] as f64;
            let angle = 2.0 * PI / seg;
            if angle < min_angle {
                break; // segments only grow — every later coupling is smaller
            }
            c.two(j, i, gates::cp_frac(profile[j], profile[i], seg));
        }
    }
    c
}

/// The mixed-radix QFT over `Z_N` in standard digit order — the tail-radix
/// phase web followed by the mirror-pair swap network. Requires a
/// palindromic profile (the digit reversal permutes site contents across the
/// mirror, which is only dimension-consistent when `d_i = d_{n-1-i}`).
pub fn mixed_radix_qft(profile: &[usize], min_angle: f64) -> Circuit {
    assert!(
        is_palindromic(profile),
        "standard-order mixed-radix QFT needs a palindromic profile; \
         use mixed_radix_qft_reversed for general chains"
    );
    let mut c = mixed_radix_qft_reversed(profile, min_angle);
    for (i, j) in crate::zigzag::mirror_pairs(profile) {
        c.two(i, j, gates::xswap(profile[i], profile[j]));
    }
    c
}

/// Draper phase ramp `D_c = diag(exp(2πi·c·y/N))` as a product of
/// single-site diagonal gates (in standard digit order):
/// `D_c = ⊗_j diag_y(exp(2πi·c·y/P_j))` with `P_j = Π_{k≤j} d_k`.
pub fn fourier_phase_ramp(profile: &[usize], c_add: usize) -> Circuit {
    let mut c = Circuit::new(profile.to_vec());
    let mut prefix = 1.0f64;
    for (j, &d) in profile.iter().enumerate() {
        prefix *= d as f64;
        let phases: Vec<f64> = (0..d)
            .map(|y| 2.0 * PI * (c_add as f64) * (y as f64) / prefix)
            .collect();
        c.one(j, gates::phase_diag(&phases));
    }
    c
}

/// The same ramp expressed on **reversed-place** digits (site `j` holding
/// the digit of place `Q_j = Π_{k<j} d_k`):
/// `D_c = ⊗_j diag_y(exp(2πi·c·y·Q_j/N)) = ⊗_j diag_y(exp(2πi·c·y/T_j))`
/// with `T_j = Π_{k≥j} d_k`.
///
/// This is the ramp that composes with [`mixed_radix_qft_reversed`]: in
/// `V† ∘ D_c ∘ V` the digit reversal cancels, so the pipeline is the exact
/// adder `|x⟩ → |x+c mod N⟩` in *standard* digits on **any** profile — no
/// mirror-swap stage, and every block stays MPO-compact.
pub fn fourier_phase_ramp_reversed(profile: &[usize], c_add: usize) -> Circuit {
    fourier_phase_ramp_reversed_real(profile, c_add as f64)
}

/// [`fourier_phase_ramp_reversed`] for a **real-valued** shift amount: the
/// factorization `exp(2πi·c·y/N) = Π_j exp(2πi·c·y_j·Q_j/N)` is exact for
/// any real `c`, so this is the generator of the *continuous* translation
/// flow on the ring — the heart of fractional operator powers
/// (see [`crate::flow`]).
pub fn fourier_phase_ramp_reversed_real(profile: &[usize], c_add: f64) -> Circuit {
    let n_total: f64 = profile.iter().map(|&d| d as f64).product();
    let mut c = Circuit::new(profile.to_vec());
    let mut q_place = 1.0f64;
    for (j, &d) in profile.iter().enumerate() {
        let phases: Vec<f64> = (0..d)
            .map(|y| 2.0 * PI * c_add * (y as f64) * q_place / n_total)
            .collect();
        c.one(j, gates::phase_diag(&phases));
        q_place *= d as f64;
    }
    c
}

/// Reversal-free Fourier-space adder: `V† ∘ D_c ∘ V` with
/// `V = mixed_radix_qft_reversed`. Valid on any profile; input and output
/// both in standard digit order.
pub fn adder_circuit_reversal_free(profile: &[usize], c_add: usize) -> Circuit {
    let mut c = mixed_radix_qft_reversed(profile, 0.0);
    c.extend(fourier_phase_ramp_reversed(profile, c_add));
    c.extend(mixed_radix_qft_reversed(profile, 0.0).adjoint());
    c
}

/// The modular adder `|x⟩ → |x + c mod N⟩` built *directly* as a bond-2
/// carry-propagation MPO — the collapsed form that Fourier-pipeline
/// composition converges toward.
///
/// Each bond carries one classical bit: the carry flowing from the least
/// significant site (`n-1`) toward the most significant (`0`, where the
/// overflow is dropped — that is the mod-`N` wrap). This form is *exact for
/// any ring size*: unlike the Fourier route, whose deepest-carry information
/// lives in tail-radix phases of order `2π/N` and falls below `f64`
/// resolution for `N ≳ 10^12`, the carry MPO is a 0/1 permutation tensor.
pub fn adder_mpo_exact(profile: &[usize], c_add: u128, trunc: TruncSpec) -> Mpo {
    let n = profile.len();
    let big_n: u128 = profile.iter().map(|&d| d as u128).product();
    let mut rem = c_add % big_n;
    // Digits of the constant, least significant first (site n-1 backwards).
    let mut c_digits = vec![0usize; n];
    for i in (0..n).rev() {
        c_digits[i] = (rem % profile[i] as u128) as usize;
        rem /= profile[i] as u128;
    }
    let tensors: Vec<SiteTensor> = (0..n)
        .map(|i| {
            let d = profile[i];
            // Left bond = carry toward site 0 (dropped at the edge);
            // right bond = carry arriving from less-significant sites.
            let dl = if i == 0 { 1 } else { 2 };
            let dr = if i == n - 1 { 1 } else { 2 };
            let mut t = SiteTensor::zeros(dl, d * d, dr);
            for p_in in 0..d {
                for carry_in in 0..dr {
                    let sum = p_in + c_digits[i] + carry_in;
                    let p_out = sum % d;
                    let carry_out = sum / d;
                    let l = if i == 0 { 0 } else { carry_out };
                    if i == 0 || carry_out < 2 {
                        t.set(l, p_out * d + p_in, carry_in, C64::ONE);
                    }
                }
            }
            t
        })
        .collect();
    let mut carrier = Mps {
        dims: profile.iter().map(|&d| d * d).collect(),
        tensors,
        center: 0,
        // Canonicalize under an *exact* policy: the deep-carry branch has
        // relative Frobenius weight ~1/N, so any weight-based cutoff would
        // amputate it — even though keeping it is free in rank. (Rare but
        // semantically critical operator branches must be compressed by
        // rank, not by weight.)
        trunc: TruncSpec::exact(),
        discarded_weight: 0.0,
    };
    carrier.recompress();
    carrier.trunc = trunc;
    Mpo {
        dims: profile.to_vec(),
        carrier,
    }
}

/// The full Fourier-space adder circuit `|x⟩ → |x + c mod N⟩`:
/// `QFT† ∘ D_c ∘ QFT`, concatenated at the gate level. (The interesting
/// route builds the same operator by *composing the three blocks as MPOs* —
/// see `examples/circuit_of_circuits.rs`.)
pub fn adder_circuit(profile: &[usize], c_add: usize) -> Circuit {
    let mut c = mixed_radix_qft(profile, 0.0);
    c.extend(fourier_phase_ramp(profile, c_add));
    c.extend(mixed_radix_qft(profile, 0.0).adjoint());
    c
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::c64::C64;
    use crate::dense::{total_dim, DenseState};
    use crate::mat::Rng;
    use crate::zigzag;

    /// Direct O(N²) DFT of a dense chain state: (Fψ)[y] = Σ_x ω^{xy} ψ[x]/√N.
    fn dft_dense(psi: &DenseState) -> DenseState {
        let n = psi.amps.len();
        let mut out = vec![C64::ZERO; n];
        for (y, o) in out.iter_mut().enumerate() {
            for (x, a) in psi.amps.iter().enumerate() {
                let phase = 2.0 * PI * ((x as u64 * y as u64) % n as u64) as f64 / n as f64;
                *o += C64::cis(phase) * *a;
            }
            *o = o.scale(1.0 / (n as f64).sqrt());
        }
        DenseState {
            dims: psi.dims.clone(),
            amps: out,
        }
    }

    #[test]
    fn reversed_qft_matches_dft_kernel() {
        // For each basis input x, the amplitude on output digits (y'_j) must
        // be ω^{x·y}/√N with y = Σ_j y'_j Q_j.
        for profile in [vec![2usize, 3], vec![2, 3, 4], vec![3, 2, 5]] {
            let n_total = total_dim(&profile);
            // Reversed place values Q_j.
            let mut q = vec![1usize; profile.len()];
            for j in 1..profile.len() {
                q[j] = q[j - 1] * profile[j - 1];
            }
            let circuit = mixed_radix_qft_reversed(&profile, 0.0);
            for x in 0..n_total {
                let mut psi = DenseState::zero_state(&profile);
                psi.amps[0] = C64::ZERO;
                psi.amps[x] = C64::ONE;
                circuit.run_dense(&mut psi);
                for (idx, amp) in psi.amps.iter().enumerate() {
                    // digits of idx are the y'_j values.
                    let digits = psi.digits_of(idx);
                    let y: usize = digits.iter().zip(q.iter()).map(|(d, qq)| d * qq).sum();
                    let expect = C64::cis(2.0 * PI * ((x * y) % n_total) as f64 / n_total as f64)
                        .scale(1.0 / (n_total as f64).sqrt());
                    assert!(
                        (*amp - expect).abs() < 1e-10,
                        "profile {:?} x={} idx={} amp {:?} expect {:?}",
                        profile,
                        x,
                        idx,
                        amp,
                        expect
                    );
                }
            }
        }
    }

    #[test]
    fn standard_qft_on_palindromes_is_the_dft() {
        let mut rng = Rng::new(97);
        for profile in [vec![2usize, 2], vec![2, 3, 2], zigzag::diamond(2, 4)] {
            let mut psi = DenseState::zero_state(&profile);
            for a in &mut psi.amps {
                *a = rng.c_gaussian();
            }
            psi.normalize();
            let expect = dft_dense(&psi);
            let mut got = psi.clone();
            mixed_radix_qft(&profile, 0.0).run_dense(&mut got);
            let diff: f64 = got
                .amps
                .iter()
                .zip(expect.amps.iter())
                .map(|(a, b)| (*a - *b).abs())
                .fold(0.0, f64::max);
            assert!(diff < 1e-9, "profile {:?}: max diff {}", profile, diff);
        }
    }

    #[test]
    fn adder_circuit_shifts_the_ring() {
        let profile = vec![2usize, 3, 2]; // N = 12, palindromic
        let n_total = total_dim(&profile);
        let c_add = 5;
        let circuit = adder_circuit(&profile, c_add);
        for x in [0usize, 3, 7, 11] {
            let mut psi = DenseState::zero_state(&profile);
            psi.amps[0] = C64::ZERO;
            psi.amps[x] = C64::ONE;
            circuit.run_dense(&mut psi);
            let target = (x + c_add) % n_total;
            assert!(
                (psi.amps[target].abs() - 1.0).abs() < 1e-9,
                "x={} → weight at target {}",
                x,
                psi.amps[target].abs()
            );
        }
    }

    #[test]
    fn reversal_free_adder_works_on_any_profile() {
        // The reversed-convention pipeline needs no palindrome: check a
        // non-palindromic chain end to end.
        for profile in [vec![2usize, 3, 4], vec![3, 5, 2, 4]] {
            let n_total = total_dim(&profile);
            let c_add = 7;
            let circuit = adder_circuit_reversal_free(&profile, c_add);
            for x in [0usize, 1, n_total / 2, n_total - 1] {
                let mut psi = DenseState::zero_state(&profile);
                psi.amps[0] = C64::ZERO;
                psi.amps[x] = C64::ONE;
                circuit.run_dense(&mut psi);
                let target = (x + c_add) % n_total;
                assert!(
                    (psi.amps[target].abs() - 1.0).abs() < 1e-9,
                    "profile {:?} x={}: weight {}",
                    profile,
                    x,
                    psi.amps[target].abs()
                );
            }
        }
    }

    #[test]
    fn exact_carry_mpo_is_the_adder() {
        use crate::mpo::Mpo;
        let spec = TruncSpec::new(64, 1e-12);
        // Small chains: against the dense cyclic shift.
        for profile in [vec![2usize, 3, 4], vec![2, 3, 2]] {
            let n_total = total_dim(&profile);
            let c_add = 5u128;
            let adder = adder_mpo_exact(&profile, c_add, spec);
            assert!(adder.max_bond_dim() <= 2, "carry MPO must be bond ≤ 2");
            for x in 0..n_total {
                let mut basis = DenseState::zero_state(&profile);
                basis.amps[0] = C64::ZERO;
                basis.amps[x] = C64::ONE;
                let psi = crate::mps::Mps::from_dense(&basis, spec);
                let out = adder.apply_to(&psi).to_dense();
                let target = (x + 5) % n_total;
                assert!(
                    (out.amps[target].abs() - 1.0).abs() < 1e-10,
                    "profile {:?} x={}",
                    profile,
                    x
                );
            }
            // And against the Fourier pipeline compilation.
            let pipeline = Mpo::from_circuit(&adder_circuit_reversal_free(&profile, 5), spec);
            let f = adder.hs_fidelity(&pipeline);
            assert!((f - 1.0).abs() < 1e-8, "hs fidelity {}", f);
        }
    }

    #[test]
    fn exact_carry_mpo_handles_the_deepest_carry_at_scale() {
        // On the 25-site wave (N ≈ 8.6e12) the +1 carry from |N-1⟩ wraps
        // through every site. The carry MPO does it exactly — this is the
        // regime where the Fourier route falls below f64 resolution.
        let wave = crate::zigzag::wave(2, 5, 4);
        let spec = TruncSpec::new(64, 1e-12);
        let big_n: u128 = wave.iter().map(|&d| d as u128).product();
        let adder = adder_mpo_exact(&wave, 1, spec);
        let all_max: Vec<usize> = wave.iter().map(|&d| d - 1).collect();
        let input = crate::mps::Mps::basis_state(&wave, &all_max, spec);
        let expect = crate::mps::Mps::basis_state(&wave, &vec![0; wave.len()], spec);
        let out = adder.apply_to(&input);
        let f = out.fidelity(&expect);
        assert!(
            (f - 1.0).abs() < 1e-10,
            "deep carry fidelity {} (N={})",
            f,
            big_n
        );
    }

    #[test]
    fn banded_qft_converges_with_cutoff() {
        // Dropping small-angle couplings gives the approximate QFT; the
        // error must shrink as the cutoff tightens.
        let profile = zigzag::diamond(2, 4);
        let mut rng = Rng::new(101);
        let mut psi = DenseState::zero_state(&profile);
        for a in &mut psi.amps {
            *a = rng.c_gaussian();
        }
        psi.normalize();
        let mut exact = psi.clone();
        mixed_radix_qft(&profile, 0.0).run_dense(&mut exact);
        let mut last = f64::INFINITY;
        for cutoff in [0.5, 0.1, 0.01, 0.0] {
            let mut approx = psi.clone();
            mixed_radix_qft(&profile, cutoff).run_dense(&mut approx);
            let infid = 1.0 - approx.fidelity(&exact);
            assert!(
                infid <= last + 1e-12,
                "cutoff {}: infidelity {} worse than {}",
                cutoff,
                infid,
                last
            );
            last = infid;
        }
        assert!(last.abs() < 1e-12);
    }
}
