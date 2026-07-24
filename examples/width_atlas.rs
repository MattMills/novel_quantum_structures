//! The width atlas: the cost landscape of modular multiplication computed
//! **before any tensor is built** — and pinned against the real operators.
//!
//! The cut-rank theorem (THEORY.md §8) says the operator Schmidt rank of
//! `×k` across a chain cut with left ring `L` and right block `S` is
//! exactly `|{⌊k·b/S⌋ mod L : b < S}|` — pure number theory. This example
//! uses it as an instrument:
//!
//! 1. the **full atlas** of `Z_2880`: waist widths of every invertible
//!    multiplier, as a histogram — how common are resonances?
//! 2. **squaring chains** toward modular exponentiation, with per-bond
//!    width profiles predicted a priori and cross-checked against exactly
//!    recompressed cascade MPOs, bond for bond;
//! 3. **intrinsic width vs machine alphabets**: multipliers whose true
//!    width sits far below every single-front machine that could build
//!    them — narrowness reachable only through composition + recompression;
//! 4. the same calculus **at scale** on the 25-site wave (`N ≈ 8.6·10¹²`),
//!    where the atlas still answers exactly (using width(k) = width(k⁻¹)
//!    to enumerate from whichever side is small).
//!
//! Run with: `cargo run --release --example width_atlas`

use novel_quantum_structures::cascade::Transducer;
use novel_quantum_structures::width::{mult_cut_width, mult_width_profile, Width};
use novel_quantum_structures::{zigzag, TruncSpec};

fn gcd(a: u128, b: u128) -> u128 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn mod_inverse(k: u128, n: u128) -> u128 {
    let (mut old_r, mut r) = (k as i128, n as i128);
    let (mut old_s, mut s) = (1i128, 0i128);
    while r != 0 {
        let q = old_r / r;
        (old_r, r) = (r, old_r - q * r);
        (old_s, s) = (s, old_s - q * s);
    }
    assert_eq!(old_r, 1);
    old_s.rem_euclid(n as i128) as u128
}

/// Exact width wherever either `k` or `k⁻¹` admits it (rank(M) = rank(M†)).
fn best_width(k: u128, k_inv: u128, l: u128, s: u128) -> Width {
    let w = mult_cut_width(k, l, s);
    if w.is_exact() {
        return w;
    }
    let wi = mult_cut_width(k_inv, l, s);
    if wi.is_exact() {
        return wi;
    }
    Width::UpperBound(w.value().min(wi.value()))
}

fn main() {
    let profile = zigzag::diamond(2, 5);
    let n: u128 = profile.iter().map(|&d| d as u128).product();
    println!(
        "diamond {:?}, ring Z_{} — waist cut (L, S) = (24, 120), geometric cap 24\n",
        profile, n
    );

    // ---- 1. The full atlas of Z_2880 ---------------------------------------
    println!("=== 1. Waist widths of every invertible multiplier (no tensors) ===\n");
    let mut hist = std::collections::BTreeMap::new();
    let mut coprime = 0u32;
    for k in 1..n {
        if gcd(k, n) != 1 {
            continue;
        }
        coprime += 1;
        let w = mult_cut_width(k, 24, 120).value();
        *hist.entry(w).or_insert(0u32) += 1;
    }
    println!("{:>8} {:>8} {:>10}", "width", "count", "");
    for (w, c) in &hist {
        let bar = "#".repeat((*c as usize) / 8);
        println!("{:>8} {:>8}   {}", w, c, bar);
    }
    let resonant: u32 = hist.iter().filter(|(w, _)| **w < 24).map(|(_, c)| c).sum();
    println!(
        "\nφ({}) = {} multipliers; {} ({:.0}%) saturate the cap, {} ({:.0}%) are resonant (< 24).",
        n,
        coprime,
        coprime - resonant,
        100.0 * (coprime - resonant) as f64 / coprime as f64,
        resonant,
        100.0 * resonant as f64 / coprime as f64
    );
    println!("The thinnest operators and their widths:");
    let mut thin: Vec<(u128, u128)> = (1..n)
        .filter(|&k| gcd(k, n) == 1)
        .map(|k| (mult_cut_width(k, 24, 120).value(), k))
        .filter(|&(w, _)| w <= 3)
        .collect();
    thin.sort();
    for (w, k) in thin.iter().take(10) {
        println!("  ×{:<6} width {}", k, w);
    }

    // ---- 2. Squaring chains, predicted then measured -----------------------
    println!("\n=== 2. Squaring chains: bond profiles predicted a priori ===\n");
    for base in [7u128, 11] {
        let mut k = base;
        println!("orbit of {} under squaring mod {}:", base, n);
        for _ in 0..4 {
            let prof: Vec<String> = mult_width_profile(&profile, k)
                .iter()
                .map(|w| w.to_string())
                .collect();
            println!("  ×{:<6} predicted bonds [{}]", k, prof.join(", "));
            k = k * k % n;
        }
    }
    println!("\ncross-check against the real operators (exact recompression):");
    let exact = TruncSpec::exact();
    let p4 = zigzag::diamond(2, 4);
    let m7 = Transducer::mult(&p4, 7).to_mpo(exact);
    let m49 = m7.compose_after(&m7, exact);
    let m2401 = m49.compose_after(&m49, exact);
    for (k, m) in [(7u128, &m7), (49, &m49), (2401, &m2401)] {
        let predicted: Vec<usize> = mult_width_profile(&p4, k)
            .iter()
            .map(|w| w.value() as usize)
            .collect();
        let measured = m.bond_dims();
        println!(
            "  diamond(2,4) ×{:<5} predicted {:?}  measured {:?}  match {}",
            k,
            predicted,
            measured,
            predicted == measured
        );
        assert_eq!(predicted, measured);
    }

    // ---- 3. Intrinsic width vs every single-front machine ------------------
    println!("\n=== 3. Operators thinner than any machine that builds them ===\n");
    println!("A single-front machine for ×k has alphabet min(k, k⁻¹, N−k, (N−k)⁻¹)");
    println!("(forward carries, opposed remainders, and their reflections). The");
    println!("intrinsic width can sit far below all four:\n");
    println!(
        "{:>8} {:>10} {:>16} {:>10}",
        "k", "intrinsic", "best alphabet", "ratio"
    );
    let mut divergent: Vec<(u128, u128, u128)> = (2..n)
        .filter(|&k| gcd(k, n) == 1)
        .map(|k| {
            let ki = mod_inverse(k, n);
            let alpha = k.min(ki).min(n - k).min(n - ki);
            (mult_cut_width(k, 24, 120).value(), alpha, k)
        })
        .filter(|&(w, alpha, _)| alpha > 8 * w)
        .map(|(w, alpha, k)| (alpha / w, w, k))
        .collect();
    divergent.sort_by(|a, b| b.0.cmp(&a.0));
    for &(ratio, w, k) in divergent.iter().take(5) {
        let ki = mod_inverse(k, n);
        let alpha = k.min(ki).min(n - k).min(n - ki);
        println!("{:>8} {:>10} {:>16} {:>9}×", k, w, alpha, ratio);
    }
    println!(
        "\n{} of {} multipliers have intrinsic width below every single-front",
        divergent.len(),
        coprime
    );
    println!("alphabet by a factor > 8: their narrow presentations are reachable only");
    println!("by composition + recompression — no sweep of the chain carries so little.");

    // ---- 4. At scale: the 25-site wave -------------------------------------
    println!("\n=== 4. The wave (N ≈ 8.6·10¹²): the atlas still answers exactly ===\n");
    let wave = zigzag::wave(2, 5, 4);
    let wn: u128 = wave.iter().map(|&d| d as u128).product();
    println!("squaring orbit of 7 mod {}, waist-cut and max widths:", wn);
    // Waist cut of the central diamond: left = product of first 12 sites.
    let l_waist: u128 = wave[..12].iter().map(|&d| d as u128).product();
    let s_waist = wn / l_waist;
    println!(
        "  (central cut: L = {}, S = {}; widths marked ≤ are bounds)",
        l_waist, s_waist
    );
    let mut k = 7u128;
    println!("{:>28} {:>12} {:>12}", "k", "waist width", "max width");
    for _ in 0..6 {
        let ki = mod_inverse(k, wn);
        let waist = best_width(k, ki, l_waist, s_waist);
        let mut left = 1u128;
        let mut maxw: u128 = 0;
        let mut all_exact = true;
        for &d in &wave[..wave.len() - 1] {
            left *= d as u128;
            let w = best_width(k, ki, left, wn / left);
            all_exact &= w.is_exact();
            maxw = maxw.max(w.value());
        }
        println!(
            "{:>28} {:>12} {:>11}{}",
            k,
            waist.to_string(),
            maxw,
            if all_exact { " " } else { "≤" }
        );
        k = k * k % wn;
    }
    println!("\nverification at scale: ×7 cascade on the wave vs the atlas —");
    let m7w = Transducer::mult(&wave, 7).to_mpo(exact);
    let predicted: Vec<usize> = mult_width_profile(&wave, 7)
        .iter()
        .map(|w| w.value() as usize)
        .collect();
    let measured = m7w.bond_dims();
    println!("  predicted {:?}", predicted);
    println!("  measured  {:?}", measured);
    assert_eq!(predicted, measured);
    println!("  bond-for-bond agreement across all {} cuts.", measured.len());

    println!("\nThe cost curve of modular exponentiation on a chain is a number-");
    println!("theoretic object, computable in advance: the operator only executes it.");
}
