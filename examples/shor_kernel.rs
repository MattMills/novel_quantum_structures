//! The Shor kernel, end-to-end on the diamond: controlled modular
//! exponentiation, phase estimation with one reused valley qubit, and
//! order recovery by continued fractions — with every operator width
//! predicted by number theory before construction.
//!
//! Layout: the diamond `[2,3,4,5,4,3,2]` splits into a **control** — the
//! edge valley qubit at site 0, a literal qubit — and an **arithmetic
//! register** — sites 1..6, encoding the sub-ring `Z_1440`. The two Shor
//! registers live on one dimension wave at different scales.
//!
//! The controlled multiplier is boundary-machine algebra:
//! `Transducer::mult_skipping` builds `I ⊗ (×k mod 1440)` (the carry
//! front tunnels through the control site), and `Mpo::select_on` composes
//! projector branches into `C-U = |0⟩⟨0|⊗I + |1⟩⟨1|⊗(×k)`. Powers
//! `C-U^{2^j}` come from operator squaring — and their exact bond
//! profiles obey the **controlled-width formula** (THEORY.md §8.11):
//! across a register cut `(L|S)`, `χ(C-U) = w + 1` if `k ≢ 1 (mod S)`,
//! and `χ(C-U) = w` — *control for free* — if `k ≡ 1 (mod S)`. The deep
//! QPE powers are resonant on every cut, so their controlled forms cost
//! exactly their uncontrolled widths.
//!
//! Measured here: (1) the atlas prices all ten `C-U^{2^j}` a priori,
//! bond-for-bond; (2) phase kickback is *bond-free* — on an eigenstate the
//! control cut stays at χ = 1 while the phase lands; (3) semiclassical
//! (Kitaev) phase estimation with the single reused control qubit, plus
//! continued fractions, recovers `ord(7 mod 1440) = 12`.
//!
//! Run with: `cargo run --release --example shor_kernel`

use novel_quantum_structures::cascade::Transducer;
use novel_quantum_structures::dense::DenseState;
use novel_quantum_structures::mpo::Mpo;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::width::mult_width_profile;
use novel_quantum_structures::{gates, zigzag, Rng, TruncSpec, C64};
use std::f64::consts::PI;

const SPEC: TruncSpec = TruncSpec {
    max_rank: 96,
    cutoff: 1e-12,
};

fn digits_of(profile: &[usize], mut v: u128) -> Vec<usize> {
    let mut d = vec![0usize; profile.len()];
    for (i, &dim) in profile.iter().enumerate().rev() {
        d[i] = (v % dim as u128) as usize;
        v /= dim as u128;
    }
    d
}

/// Last continued-fraction convergent of `x` with denominator ≤ `q_max`.
fn cf_convergent(mut x: f64, q_max: u64) -> (u64, u64) {
    let (mut p0, mut q0): (u64, u64) = (1, 0);
    let (mut p1, mut q1): (u64, u64) = (x.floor() as u64, 1);
    for _ in 0..30 {
        let f = x - x.floor();
        if f < 1e-9 {
            break;
        }
        x = 1.0 / f;
        let a = x.floor() as u64;
        let p2 = a.saturating_mul(p1).saturating_add(p0);
        let q2 = a.saturating_mul(q1).saturating_add(q0);
        if q2 > q_max {
            break;
        }
        (p0, q0, p1, q1) = (p1, q1, p2, q2);
    }
    (p1, q1)
}

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// ⟨X⟩ and ⟨Y⟩ of the control qubit at site 0.
fn control_xy(psi: &Mps) -> (f64, f64) {
    let mut a = psi.clone();
    a.apply1(0, &gates::hadamard());
    let px = a.site_probabilities(0)[0];
    let mut b = psi.clone();
    b.apply1(0, &gates::phase_diag(&[0.0, -PI / 2.0]));
    b.apply1(0, &gates::hadamard());
    let py = b.site_probabilities(0)[0];
    (2.0 * px - 1.0, 2.0 * py - 1.0)
}

fn main() {
    let profile = zigzag::diamond(2, 5);
    let sub_profile: Vec<usize> = profile[1..].to_vec();
    let n_sub: u128 = sub_profile.iter().map(|&d| d as u128).product();
    let k = 7u128;
    // ord(7) mod 1440, classically (the answer the quantum protocol must find).
    let mut ord = 1u128;
    let mut acc = k;
    while acc != 1 {
        acc = acc * k % n_sub;
        ord += 1;
    }
    println!("diamond {:?}: control = valley qubit (site 0),", profile);
    println!(
        "register = sites 1..6 = Z_{};  ord({} mod {}) = {} (to be recovered)\n",
        n_sub, k, n_sub, ord
    );

    // ---- 1. The atlas prices the controlled powers -------------------------
    println!("=== 1. Pricing C-U^(2^j) before building it ===\n");
    let exact = TruncSpec::exact();
    let id = Mpo::identity(&profile, exact);
    // Sub-cut right-block sizes S for the controlled-width indicator.
    let mut s_right = vec![0u128; sub_profile.len() - 1];
    {
        let mut s = 1u128;
        for i in (1..sub_profile.len()).rev() {
            s *= sub_profile[i] as u128;
            s_right[i - 1] = s;
        }
    }
    // Build M' powers by operator squaring; price and verify each distinct one.
    let m_bits = 10usize;
    let mut powers: Vec<Mpo> = Vec::with_capacity(m_bits);
    let mut k_vals: Vec<u128> = Vec::with_capacity(m_bits);
    let mut kk = k;
    for j in 0..m_bits {
        let mp = if j == 0 {
            Transducer::mult_skipping(&profile, 0, k as usize).to_mpo(exact)
        } else {
            powers[j - 1].compose_after(&powers[j - 1], exact)
        };
        powers.push(mp);
        k_vals.push(kk);
        kk = kk * kk % n_sub;
    }
    println!(
        "{:>6} {:>28} {:>28}",
        "k", "predicted C-U bonds", "measured"
    );
    let mut seen: Vec<u128> = Vec::new();
    for j in 0..m_bits {
        if seen.contains(&k_vals[j]) {
            continue;
        }
        seen.push(k_vals[j]);
        let w = mult_width_profile(&sub_profile, k_vals[j]);
        let mut predicted = vec![2usize];
        for (i, wi) in w.iter().enumerate() {
            let free = k_vals[j] % s_right[i] == 1;
            predicted.push(wi.value() as usize + if free { 0 } else { 1 });
        }
        let cu = Mpo::select_on(0, &[id.clone(), powers[j].clone()], exact);
        let measured = cu.bond_dims();
        println!(
            "{:>6} {:>28} {:>28}",
            k_vals[j],
            format!("{:?}", predicted),
            format!("{:?}", measured)
        );
        assert_eq!(measured, predicted, "pricing failed for ×{}", k_vals[j]);
    }
    println!("\nEvery bond agrees with the controlled-width formula w + [k ≢ 1 mod S].");
    println!("The deep powers (961, 481) are resonant on every cut — k ≡ 1 mod S");
    println!("everywhere — so their CONTROL IS FREE: the controlled operator is no");
    println!("wider than the bare multiplier. The whole protocol is priced by number");
    println!("theory before a single tensor exists.\n");

    // The controlled powers used by the protocol.
    let cu: Vec<Mpo> = powers
        .iter()
        .map(|p| Mpo::select_on(0, &[id.clone(), p.clone()], exact))
        .collect();

    // ---- 2. Phase kickback is bond-free ------------------------------------
    println!("=== 2. Phase kickback without bond dimension ===\n");
    let r = ord as usize;
    let s_eig = 1usize;
    let mut u1 = DenseState::zero_state(&profile);
    u1.amps[0] = C64::ZERO;
    let mut v = 1usize;
    let amp = 1.0 / (r as f64).sqrt();
    for j in 0..r {
        let phase = -2.0 * PI * (j * s_eig) as f64 / r as f64;
        u1.amps[v] = C64::cis(phase).scale(amp);
        v = v * k as usize % n_sub as usize;
    }
    u1.apply1(0, &gates::hadamard()); // control |+⟩
    let base = Mps::from_dense(&u1, SPEC);
    println!("register in the eigenstate |u_1⟩ of ×7 (eigenvalue e^(2πi/12)):");
    println!(
        "{:>8} {:>12} {:>12} {:>14} {:>12}",
        "power", "θ measured", "θ predicted", "control bond", ""
    );
    for t in 0..4 {
        let out = cu[t].apply_to(&base);
        let (x, y) = control_xy(&out);
        let theta = y.atan2(x).rem_euclid(2.0 * PI);
        let theta_pred =
            2.0 * PI * ((s_eig as u128 * (1u128 << t)) % ord) as f64 / ord as f64;
        assert!(
            (theta - theta_pred).abs() < 1e-6,
            "t={}: {} vs {}",
            t,
            theta,
            theta_pred
        );
        assert_eq!(out.bond_dims()[0], 1, "kickback must be bond-free");
        println!(
            "{:>8} {:>11.4}° {:>11.4}° {:>14} {:>12}",
            format!("2^{}", t),
            theta.to_degrees(),
            theta_pred.to_degrees(),
            out.bond_dims()[0],
            "χ = 1 ✓"
        );
    }
    // Contrast: a non-eigenstate entangles the control.
    let mut one = Mps::basis_state(&profile, &digits_of(&profile, 1), SPEC);
    one.apply1(0, &gates::hadamard());
    let ent = cu[0].apply_to(&one);
    println!(
        "\ncontrast — register |1⟩ (a superposition of all 12 eigenstates):"
    );
    println!(
        "after C-U the control bond is {} — kickback is bond-free ONLY on",
        ent.bond_dims()[0]
    );
    println!("eigenstates; that is the resource statement of phase estimation.\n");

    // ---- 3. Semiclassical order finding ------------------------------------
    println!("=== 3. Kitaev phase estimation with one reused control qubit ===\n");
    println!(
        "register |1⟩, {} bits per run, measurement by projection, feedback",
        m_bits
    );
    println!("rotations conditioned on earlier bits; continued fractions per run:\n");
    let proj = |b: usize| {
        novel_quantum_structures::Mat::from_fn(2, 2, |rr, cc| {
            if rr == b && cc == b {
                C64::ONE
            } else {
                C64::ZERO
            }
        })
    };
    let mut r_lcm = 1u64;
    for run in 0..6 {
        let mut rng = Rng::new(20260724 + run);
        let mut psi = Mps::basis_state(&profile, &digits_of(&profile, 1), SPEC);
        let mut bits = vec![0u8; m_bits + 1]; // bits[t], t = 1..=m_bits
        let mut last = 0usize;
        for t in (1..=m_bits).rev() {
            if last == 1 {
                psi.apply1(0, &gates::shift(2)); // X: reset control to |0⟩
            }
            psi.apply1(0, &gates::hadamard());
            psi = cu[t - 1].apply_to(&psi);
            let omega: f64 = (t + 1..=m_bits)
                .map(|l| bits[l] as f64 / 2f64.powi((l - t + 1) as i32))
                .sum();
            psi.apply1(0, &gates::phase_diag(&[0.0, -2.0 * PI * omega]));
            psi.apply1(0, &gates::hadamard());
            let p = psi.site_probabilities(0);
            let b = if rng.uniform() < p[0] { 0 } else { 1 };
            psi.apply1(0, &proj(b));
            psi.normalize();
            bits[t] = b as u8;
            last = b;
        }
        let phi: f64 = (1..=m_bits).map(|t| bits[t] as f64 / 2f64.powi(t as i32)).sum();
        // Convergents are already in lowest terms; the denominator divides
        // r. Any unit's order divides the Carmichael exponent
        // λ(1440) = lcm(λ(32), λ(9), λ(5)) = 24 — the classical bound that
        // plays the role of Shor's r < N when picking the convergent.
        let (s_frac, q_frac) = cf_convergent(phi, 24);
        let bit_str: String = (1..=m_bits).map(|t| char::from(b'0' + bits[t])).collect();
        println!(
            "  run {}: bits 0.{}  φ = {:.5}  ≈ {}/{}",
            run, bit_str, phi, s_frac, q_frac
        );
        r_lcm = r_lcm / gcd(r_lcm, q_frac.max(1)) * q_frac.max(1);
    }
    // The lcm of denominators divides the order; find the smallest verified multiple.
    let mut recovered = 0u64;
    for mult in 1..=8 {
        let cand = r_lcm * mult;
        let mut acc = 1u128;
        for _ in 0..cand {
            acc = acc * k % n_sub;
        }
        if acc == 1 {
            recovered = cand;
            break;
        }
    }
    println!("\nlcm of denominators = {}; smallest verified multiple with", r_lcm);
    println!(
        "7^r ≡ 1 (mod 1440):  r = {}   (classical check: ord = {})",
        recovered, ord
    );
    assert_eq!(recovered as u128, ord, "order recovery failed");

    println!("\nThe complete Shor kernel — controlled modular exponentiation, phase");
    println!("estimation, order recovery — executed as boundary machines on one");
    println!("dimension wave, with the control register living in the wave's valleys");
    println!("and every operator width predicted by the atlas in advance.");
}
