//! The Fourier frame as a flow: fractional powers of the QFT itself,
//! built from four operators the crate already owns.
//!
//! `F⁴ = I` on any ring `Z_N`, so the spectral projectors of the QFT are
//! polynomials in `F` and the matrix-power fractional Fourier transform is
//! the exact four-term combination
//!
//! ```text
//!   F^t = c₀(t)·I + c₁(t)·F + c₂(t)·Π + c₃(t)·F†,     Π: x → −x,
//! ```
//!
//! assembled by operator linear algebra ([`Mpo::add`]/[`Mpo::scale`]) —
//! no eigensolver. This answers an open direction of THEORY.md §12 and
//! asks the next question: does "width sees the integers" (README
//! finding 12) generalize from the *arithmetic* flow to the *frame* flow?
//!
//! Measured answer: **operator entanglement sees every integer** — the
//! opEE curve is pinned to the pure-power value at each whole `t`, an
//! extremum each time (minimum at even `t`: 0 and 1 bit; maximum at odd
//! `t`: 5.170 bits) — but **bond dimension collapses only at
//! `t ≡ 0, 2 (mod 4)`**, because `F` itself saturates the geometric cap
//! while `I` (χ = 1) and the reflection Π (χ = 2, a width-2 machine) sit
//! far below it. Width sees the frame exactly where the frame is narrower
//! than the geometry.
//!
//! Run with: `cargo run --release --example fractional_fourier`

use novel_quantum_structures::dense::total_dim;
use novel_quantum_structures::flow::FrameFlow;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{zigzag, TruncSpec};

const SPEC: TruncSpec = TruncSpec {
    max_rank: 96,
    cutoff: 1e-12,
};

fn digits(profile: &[usize], mut x: usize) -> Vec<usize> {
    let mut d = vec![0usize; profile.len()];
    for (i, &dim) in profile.iter().enumerate().rev() {
        d[i] = x % dim;
        x /= dim;
    }
    d
}

fn main() {
    let profile = zigzag::diamond(2, 4); // [2,3,4,3,2], N = 144
    let n_total = total_dim(&profile);
    println!(
        "diamond {:?}, ring Z_{} — the flow F^t through the Fourier frame\n",
        profile, n_total
    );
    let flow = FrameFlow::qft(&profile, SPEC);

    println!("the four pure powers and the geometry:");
    let names = ["I", "F", "Π = x→−x", "F†"];
    for (m, name) in names.iter().enumerate() {
        println!(
            "  F^{}  = {:<10} χ = {:>2}   opEE(max) = {:.3} bits",
            m,
            name,
            flow.powers[m].max_bond_dim(),
            flow.powers[m]
                .operator_entanglement_bits()
                .iter()
                .cloned()
                .fold(0.0, f64::max)
        );
    }
    println!("  operator-width caps of this profile: [4, 36, 36, 4]\n");

    // ---- 1. Scanning the frame flow ----------------------------------------
    println!("=== 1. χ and operator entanglement along F^t, t ∈ [0, 4] ===\n");
    let x0 = 100usize;
    let start = Mps::basis_state(&profile, &digits(&profile, x0), SPEC);
    println!(
        "{:>6} {:>5} {:>10} {:>15}  width",
        "t", "chi", "opEE(max)", "participation"
    );
    for k in 0..=16 {
        let t = k as f64 / 4.0;
        let u = flow.at(t);
        let chi = u.max_bond_dim();
        let opee = u
            .operator_entanglement_bits()
            .iter()
            .cloned()
            .fold(0.0, f64::max);
        let out = u.apply_to(&start).to_dense();
        let participation: f64 = out.amps.iter().map(|a| a.abs2() * a.abs2()).sum();
        let bar = "#".repeat(chi.min(60));
        println!(
            "{:>6.2} {:>5} {:>10.3} {:>15.4}  {}",
            t, chi, opee, participation, bar
        );
    }
    println!("\nopEE is pinned to the pure-power value at EVERY integer — an extremum");
    println!("each time (minima 0 and 1 bit at even t, maxima 5.170 at odd t) — but");
    println!("χ collapses only at t ≡ 0, 2 (mod 4), where the pure power (I, Π) is");
    println!("narrower than the geometric cap that F itself saturates. Participation");
    println!("of F^t|x₀⟩ breathes with period 2: localized at even t (identity /");
    println!("reflection), maximally flat at odd t (Fourier / inverse Fourier).\n");

    // ---- 2. The group law, measured ----------------------------------------
    println!("=== 2. An exact one-parameter group of period 4 ===\n");
    let half = flow.at(0.5);
    let f1 = half.compose_after(&half, SPEC).hs_fidelity(&flow.powers[1]);
    println!("F^½ ∘ F^½  = F    :  hs fidelity {:.9}", f1);
    let f2 = flow
        .at(1.3)
        .compose_after(&flow.at(0.7), SPEC)
        .hs_fidelity(&flow.powers[2]);
    println!("F^0.7 ∘ F^1.3 = Π :  hs fidelity {:.9}", f2);
    let f3 = flow
        .at(3.5)
        .compose_after(&flow.at(0.5), SPEC)
        .hs_fidelity(&flow.powers[0]);
    println!("F^3.5 ∘ F^0.5 = I :  hs fidelity {:.9}", f3);
    let ff = flow.powers[1]
        .compose_after(&flow.powers[1], SPEC)
        .hs_fidelity(&flow.powers[2]);
    println!("F ∘ F = Π         :  hs fidelity {:.9}   (the reflection, χ = 2)", ff);
    let out = flow.at(0.7).apply_to(&start);
    println!(
        "unitarity of a fractional snapshot: ‖F^0.7|x₀⟩‖ = {:.9}",
        out.norm()
    );

    println!("\nThe QFT is not one gate but a curve, and the curve is spanned by just");
    println!("four operators — the frame's own algebra (F⁴ = I) is what makes its");
    println!("fractional powers first-class objects on the chain.");
}
