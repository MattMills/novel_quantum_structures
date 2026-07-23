//! Stepwise cascades: operators that compute by sweeping a message along
//! the chain.
//!
//! A cascade (a [`Transducer`]) is a finite-state machine lifted to a
//! quantum operator: at each site it reads (message, digit), writes
//! (message, digit) by a reversible local rule, and moves on — the message
//! rides the MPO bond, so **bond dimension is literally the width of the
//! classical information front**. The carry adder is the m = 2 member;
//! this example explores the family:
//!
//! 1. modular multiplication `x → kx mod N` as a message-`k` cascade,
//!    exact on the diamond's ring and on the 8.6-trillion-element wave;
//! 2. recursive composition: multipliers compose into multipliers,
//!    adders and multipliers braid into affine maps;
//! 3. **repeated squaring toward modular exponentiation**, measuring the
//!    *true operator width* of `mult(7^(2^j) mod N)` at each squaring —
//!    where cascade width meets the chain's geometric budget;
//! 4. the Fourier scaling theorem as cascade geometry: conjugating a
//!    multiplier by the frame inverts `k` and *reverses the sweep
//!    direction*;
//! 5. the gcd obstruction: `gcd(k, N) ≠ 1` breaks unitarity, visible as a
//!    norm anomaly and colliding images.
//!
//! Run with: `cargo run --release --example stepwise_cascade`
//! (the big squaring step takes a minute or two)

use novel_quantum_structures::cascade::{unitarity_defect, Direction, Transducer};
use novel_quantum_structures::dense::total_dim;
use novel_quantum_structures::mpo::Mpo;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{radix, zigzag, TruncSpec};
use std::time::Instant;

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

fn main() {
    let profile = zigzag::diamond(2, 5);
    let n_total = total_dim(&profile) as u128;
    println!(
        "diamond {:?}, ring Z_{} — cascades sweep right-to-left (carry direction)\n",
        profile, n_total
    );

    // ---- 1. Multiplication as a message-k cascade --------------------------
    println!("=== 1. Modular multiplication as a stepwise cascade ===\n");
    println!("adder (+c):    message dim 2 (the carry bit)");
    println!("mult (×k):     message dim k (the multiplication carry)\n");
    for k in [7usize, 11] {
        let m = Transducer::mult(&profile, k).to_mpo(SPEC);
        println!(
            "×{:>2}: true bond profile {:?}  (message bound {}, geometry caps [4,36,576,576,36,4])",
            k,
            m.bond_dims(),
            k
        );
        // Spot-verify against u128 arithmetic, including the worst carry.
        for x in [1u128, 411, n_total - 1] {
            let f = m
                .apply_to(&basis(&profile, x))
                .fidelity(&basis(&profile, k as u128 * x % n_total));
            assert!((f - 1.0).abs() < 1e-8);
        }
        println!(
            "     verified on basis states; unitarity defect {:.1e}",
            unitarity_defect(&m, SPEC)
        );
    }
    // At scale: the wave.
    let wave = zigzag::wave(2, 5, 4);
    let wave_n: u128 = wave.iter().map(|&d| d as u128).product();
    let t0 = Instant::now();
    let m7w = Transducer::mult(&wave, 7).to_mpo(SPEC);
    let f = m7w
        .apply_to(&basis(&wave, wave_n - 1))
        .fidelity(&basis(&wave, 7 * (wave_n - 1) % wave_n));
    println!(
        "\n×7 on the 25-site wave (N = {}): built in {:.0?}, worst-case carry fidelity {:.9}\n",
        wave_n,
        t0.elapsed(),
        f
    );

    // ---- 2. Recursive composition ------------------------------------------
    println!("=== 2. Cascades compose recursively ===\n");
    let m7 = Transducer::mult(&profile, 7).to_mpo(SPEC);
    let m11 = Transducer::mult(&profile, 11).to_mpo(SPEC);
    let m77 = Transducer::mult(&profile, 77).to_mpo(SPEC);
    let composed = m7.compose_after(&m11, SPEC);
    println!(
        "×7 ∘ ×11 = ×77:  hs fidelity {:.9}   (composed bonds {:?}, direct-77 bonds {:?})",
        composed.hs_fidelity(&m77),
        composed.bond_dims(),
        m77.bond_dims()
    );
    let a5 = Transducer::adder(&profile, 5).to_mpo(SPEC);
    let a35 = Transducer::adder(&profile, 35).to_mpo(SPEC);
    let lhs = m7.compose_after(&a5, SPEC);
    let rhs = a35.compose_after(&m7, SPEC);
    println!(
        "×7 ∘ (+5) = (+35) ∘ ×7:  hs fidelity {:.9}   (affine maps braid)\n",
        lhs.hs_fidelity(&rhs)
    );

    // ---- 3. Repeated squaring: the cost curve of exponentiation ------------
    println!("=== 3. Repeated squaring toward modular exponentiation ===\n");
    println!("squaring ×7 cascades: k → k² mod 2880, composing operators exactly");
    println!("(exact recompression reveals the TRUE operator width of each power)\n");
    let exact = TruncSpec::exact();
    let mut acc = Transducer::mult(&profile, 7).to_mpo(exact);
    let mut k_val: u128 = 7;
    println!(
        "{:>12} {:>8} {:>26}",
        "multiplier", "max chi", "bond profile"
    );
    println!(
        "{:>12} {:>8} {:>26?}",
        k_val,
        acc.max_bond_dim(),
        acc.bond_dims()
    );
    for _ in 0..3 {
        let t = Instant::now();
        acc = acc.compose_after(&acc, exact);
        k_val = k_val * k_val % n_total;
        // Verify against u128 arithmetic on sampled basis states.
        for x in [1u128, 1439, n_total - 1] {
            let f = acc
                .apply_to(&basis(&profile, x))
                .fidelity(&basis(&profile, k_val * x % n_total));
            assert!((f - 1.0).abs() < 1e-7, "x={}: {}", x, f);
        }
        println!(
            "{:>12} {:>8} {:>26?}   ({:.1?})",
            k_val,
            acc.max_bond_dim(),
            acc.bond_dims(),
            t.elapsed()
        );
    }
    println!("\n7 → 49 → 2401 → 1921 (mod 2880), each verified on basis states.");
    println!("Message dimension is only an upper bound: the measured operator width");
    println!("is capped by the chain's geometric budget — cascade width and zigzag");
    println!("geometry meet in the same number.\n");

    // ---- 4. Fourier conjugation reverses the cascade -----------------------
    println!("=== 4. The Fourier frame inverts k and reverses the sweep ===\n");
    let small = vec![2usize, 3, 4, 3]; // N = 72, non-palindromic
    let v = Mpo::from_circuit(&radix::mixed_radix_qft_reversed(&small, 0.0), SPEC);
    let mk = Transducer::mult(&small, 5).to_mpo(SPEC);
    let conj = v.compose_after(&mk.compose_after(&v.adjoint(), SPEC), SPEC);
    let mk_inv_rev = Transducer::mult_directed(&small, 29, Direction::LeftToRight).to_mpo(SPEC);
    println!("on Z_72:  V ∘ (×5) ∘ V†  vs  (×29 with carries sweeping the other way):");
    println!(
        "hs fidelity {:.9}   (5·29 ≡ 1 mod 72 — inverted multiplier, reversed direction)\n",
        conj.hs_fidelity(&mk_inv_rev)
    );

    // ---- 5. The gcd obstruction --------------------------------------------
    println!("=== 5. Number theory as an operator property ===\n");
    let m6 = Transducer::mult(&profile, 6).to_mpo(SPEC); // gcd(6, 2880) = 6
    println!(
        "×6 on Z_2880 (gcd = 6): unitarity defect {:.3}   (×7: {:.1e})",
        unitarity_defect(&m6, SPEC),
        unitarity_defect(&m7, SPEC)
    );
    let img_a = m6.apply_to(&basis(&profile, 0));
    let img_b = m6.apply_to(&basis(&profile, 480)); // 6·480 = 2880 ≡ 0
    println!(
        "collision: |0⟩ and |480⟩ both map to |0⟩ — image overlap {:.9}",
        img_a.fidelity(&img_b)
    );
    println!("\nA cascade is unitary exactly when its arithmetic is invertible mod N:");
    println!("gcd(k, N) = 1. The obstruction is not a numerical artifact — it is the");
    println!("operator remembering number theory.");
}
