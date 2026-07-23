//! Geometrically opposed state sets: doing the operation with the dual
//! machine.
//!
//! Multiplication's carries flow LSB → MSB. Division's remainders flow
//! **MSB → LSB** — the opposed front. Both compute with a state set of size
//! `k`, but they compute *inverse* operations: so any `×k⁻¹ mod N`, whose
//! forward carry machine would need `k⁻¹ mod N` states (astronomical), runs
//! as the opposed remainder machine with just `k` states.
//!
//! The ingredient that makes the opposed machine runnable is quantum
//! boundary conditions: the remainder front *enters in a superposition of
//! every initial state* (the unknown wrap multiple `m`) and is
//! *postselected at exit* (`remainder = 0`, exact division) — the MPO
//! contraction keeps, for each input, exactly the one branch whose guess
//! was consistent. A classical FSM cannot do this; a tensor-network
//! operator does it for free.
//!
//! Also here: opposed fronts annihilate (`÷k ∘ ×k = 1`), two narrow
//! opposed fronts realize modular ratios whose single-direction machine
//! would be enormous, and the *other* geometric opposition — ring
//! reflection via digit complement — gives negation at machine width 2.
//!
//! Run with: `cargo run --release --example opposed_fronts`

use novel_quantum_structures::cascade::Transducer;
use novel_quantum_structures::circuit::Circuit;
use novel_quantum_structures::dense::total_dim;
use novel_quantum_structures::mpo::Mpo;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{gates, zigzag, TruncSpec};

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

/// Modular inverse by extended Euclid (requires gcd(k, n) = 1).
fn mod_inverse(k: u128, n: u128) -> u128 {
    let (mut old_r, mut r) = (k as i128, n as i128);
    let (mut old_s, mut s) = (1i128, 0i128);
    while r != 0 {
        let q = old_r / r;
        (old_r, r) = (r, old_r - q * r);
        (old_s, s) = (s, old_s - q * s);
    }
    assert_eq!(old_r, 1, "gcd(k, N) must be 1");
    old_s.rem_euclid(n as i128) as u128
}

fn main() {
    let profile = zigzag::diamond(2, 5);
    let n_total = total_dim(&profile) as u128;
    println!("diamond {:?}, ring Z_{}\n", profile, n_total);
    println!("×k: carry front sweeps LSB → MSB, state set [0, k)");
    println!("÷k: remainder front sweeps MSB → LSB, state set [0, k) — the opposed dual");
    println!("boundaries: enter in superposition over the wrap guess, postselect exact division\n");

    // ---- 1. The inverse-width collapse -------------------------------------
    println!("=== 1. Computing ×k⁻¹ with the opposed machine ===\n");
    let k = 7u128;
    let k_inv = mod_inverse(k, n_total);
    println!("on Z_{}: 7⁻¹ = {}", n_total, k_inv);
    println!(
        "  forward carry machine for ×{}: message width {}",
        k_inv, k_inv
    );
    println!("  opposed remainder machine ÷7:  message width {}\n", k);
    let dv7 = Transducer::div_mpo(&profile, 7, SPEC);
    let m7 = Transducer::mult(&profile, 7).to_mpo(SPEC);
    println!("÷7 bond profile {:?}", dv7.bond_dims());
    println!(
        "÷7 vs (×7)†:  hs fidelity {:.9}   (the adjoint bound, realized as a machine)",
        dv7.hs_fidelity(&m7.adjoint())
    );
    for x in [1u128, 411, n_total - 1] {
        let f = dv7
            .apply_to(&basis(&profile, x))
            .fidelity(&basis(&profile, k_inv * x % n_total));
        assert!((f - 1.0).abs() < 1e-8);
    }
    println!(
        "|x⟩ → |{}·x mod {}⟩ verified on basis states",
        k_inv, n_total
    );

    // At wave scale the collapse is astronomical.
    let wave = zigzag::wave(2, 5, 4);
    let wave_n: u128 = wave.iter().map(|&d| d as u128).product();
    let wave_inv = mod_inverse(7, wave_n);
    let dvw = Transducer::div_mpo(&wave, 7, SPEC);
    let mw = Transducer::mult(&wave, 7).to_mpo(SPEC);
    let x0 = wave_n - 1;
    let round_trip = dvw
        .apply_to(&mw.apply_to(&basis(&wave, x0)))
        .fidelity(&basis(&wave, x0));
    println!("\non the 25-site wave (N = {}):", wave_n);
    println!(
        "  7⁻¹ = {} — forward machine width ≈ {:.1e}",
        wave_inv, wave_inv as f64
    );
    println!(
        "  opposed machine width 7; ÷7∘×7 round trip on |N-1⟩: fidelity {:.9}\n",
        round_trip
    );

    // ---- 2. Opposed fronts annihilate --------------------------------------
    println!("=== 2. Opposed fronts annihilate ===\n");
    let id = Mpo::identity(&profile, SPEC);
    let ann = dv7.compose_after(&m7, SPEC);
    println!(
        "÷7 ∘ ×7 = 1:  hs fidelity {:.9}, bond profile {:?}\n",
        ann.hs_fidelity(&id),
        ann.bond_dims()
    );

    // ---- 3. Dual narrow fronts realize a huge single-front machine ---------
    println!("=== 3. Two opposed fronts vs one enormous front ===\n");
    let inv11 = mod_inverse(11, n_total);
    let ratio_k = 7 * inv11 % n_total;
    println!(
        "×(7/11) on Z_{} = ×{}   (11⁻¹ = {})",
        n_total, ratio_k, inv11
    );
    println!("  single-direction machine: message width {}", ratio_k);
    println!("  dual opposed fronts:      ×7 (width 7) then ÷11 (width 11)");
    let dv11 = Transducer::div_mpo(&profile, 11, SPEC);
    let ratio = dv11.compose_after(&m7, SPEC);
    for x in [1u128, 1000, n_total - 1] {
        let f = ratio
            .apply_to(&basis(&profile, x))
            .fidelity(&basis(&profile, ratio_k * x % n_total));
        assert!((f - 1.0).abs() < 1e-8, "x={}", x);
    }
    println!(
        "  composed ÷11∘×7: bond profile {:?}, verified as ×{} on basis states\n",
        ratio.bond_dims(),
        ratio_k
    );

    // ---- 4. The other opposition: ring reflection --------------------------
    println!("=== 4. Reflection duality: negation at machine width 2 ===\n");
    let mut comp_circuit = Circuit::new(profile.clone());
    for (i, &d) in profile.iter().enumerate() {
        comp_circuit.one(i, gates::complement(d));
    }
    let comp = Mpo::from_circuit(&comp_circuit, SPEC);
    let plus1 = Transducer::adder(&profile, 1).to_mpo(SPEC);
    let neg = plus1.compose_after(&comp, SPEC);
    println!(
        "x → -x  =  (+1) ∘ complement:  bond profile {:?}   (naive ×{} machine: width {})",
        neg.bond_dims(),
        n_total - 1,
        n_total - 1
    );
    // ×(N-7) = ×(-7) = neg ∘ ×7: width max(2, 7) instead of N-7.
    let m_minus7 = neg.compose_after(&m7, SPEC);
    let target_k = n_total - 7;
    for x in [1u128, 500] {
        let f = m_minus7
            .apply_to(&basis(&profile, x))
            .fidelity(&basis(&profile, target_k * x % n_total));
        assert!((f - 1.0).abs() < 1e-8);
    }
    println!(
        "×{} = neg ∘ ×7:  bond profile {:?}, verified on basis states",
        target_k,
        m_minus7.bond_dims()
    );

    println!("\nEvery operation has a forward machine and a geometrically opposed one —");
    println!("opposite sweep direction (carry/remainder), opposite ring orientation");
    println!("(k / N-k). Their widths differ enormously; the operator is only as wide");
    println!("as the *narrowest* machine that computes it, and quantum boundary");
    println!("conditions (enter in superposition, postselect the exit) are what let");
    println!("the opposed machine actually run.");
}
