//! Operators as state-components over increasing operation width.
//!
//! A circuit block is not one object but a **flow** `t ↦ U^t` — a geodesic
//! through operator space from the identity to the full operation — and
//! every snapshot is a first-class MPO whose bond dimension is its
//! *operation width*. Four experiments on the diamond's ring `Z_2880`:
//!
//! 1. **Operator width sees the integers**: χ(U^t) scanned over t ∈ [0, 2]
//!    dips sharply at whole shifts (crisp permutations) and widens at
//!    fractional ones (delocalized kernels).
//! 2. **Geometric imputation is not unique**: the same coarse block splits
//!    into geodesic halves (U^½·U^½) or causal halves (through the Fourier
//!    frame) — same product, very different intermediate widths.
//! 3. **The width cursor**: a pipeline refined again and again at its
//!    leading edge — coarse past, increasingly fine present — total
//!    invariant at every zoom level.
//! 4. **Walking the ring continuously**: a basis state stepped through the
//!    refined pipeline drifts linearly with an exact "lattice wobble"
//!    (amplitude 1/2π, the discreteness of the ring pushing back), breathes
//!    in participation, and relocalizes exactly at integer width.
//!
//! Run with: `cargo run --release --example operator_width_cursor`

use novel_quantum_structures::c64::C64;
use novel_quantum_structures::dense::total_dim;
use novel_quantum_structures::flow::{FourierFlow, WidthCursor};
use novel_quantum_structures::mpo::Mpo;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{radix, zigzag, TruncSpec};
use std::f64::consts::PI;

const SPEC: TruncSpec = TruncSpec {
    max_rank: 96,
    cutoff: 1e-12,
};

fn main() {
    let profile = zigzag::diamond(2, 5);
    let n_total = total_dim(&profile);
    println!(
        "diamond {:?}, ring Z_{} — flow U^t = continuous shift by t ring units\n",
        profile, n_total
    );
    let flow = FourierFlow::adder(&profile, 1.0, SPEC);

    // ---- 1. Operator width sees the integers ------------------------------
    println!("=== 1. Operator width along the flow: chi(U^t) for t in [0, 2] ===\n");
    println!("{:>7} {:>5} {:>10}  width", "t", "chi", "opEE(max)");
    for k in 0..=16 {
        let t = k as f64 / 8.0;
        let u = flow.at(t);
        let chi = u.max_bond_dim();
        let opee = u
            .operator_entanglement_bits()
            .iter()
            .cloned()
            .fold(0.0, f64::max);
        let bar = "#".repeat(chi);
        println!("{:>7.3} {:>5} {:>10.3}  {}", t, chi, opee, bar);
    }
    println!("\nThe flow collapses to a narrow permutation exactly at integer widths");
    println!("and delocalizes in between: quantization, read off bond dimension.\n");

    // ---- 2. Two imputations of the same coarse block -----------------------
    println!("=== 2. Geometric imputation of a coarse block into halves ===\n");
    let whole = flow.at(1.0);
    let geo_half = flow.at(0.5);
    let v = Mpo::from_circuit(&radix::mixed_radix_qft_reversed(&profile, 0.0), SPEC);
    let half_ramp = Mpo::from_circuit(
        &radix::fourier_phase_ramp_reversed_real(&profile, 0.5),
        SPEC,
    );
    // Causal halves: W1 = D^{1/2} ∘ V (enter the frame, half the phases),
    //                W2 = V† ∘ D^{1/2} (other half, leave the frame).
    let w1 = half_ramp.compose_after(&v, SPEC);
    let w2 = v.adjoint().compose_after(&half_ramp, SPEC);

    let geo_prod = geo_half.compose_after(&geo_half, SPEC);
    let causal_prod = w2.compose_after(&w1, SPEC);
    println!(
        "coarse block A = U^1:          chi {:>3}",
        whole.max_bond_dim()
    );
    println!(
        "geodesic halves  U^0.5 · U^0.5:  chi {:>3} each   (product vs A: hs fidelity {:.9})",
        geo_half.max_bond_dim(),
        geo_prod.hs_fidelity(&whole)
    );
    println!(
        "causal halves    W2 · W1:        chi {:>3},{:>2} each (product vs A: hs fidelity {:.9})",
        w2.max_bond_dim(),
        w1.max_bond_dim(),
        causal_prod.hs_fidelity(&whole)
    );
    println!(
        "the two half-points differ:      |⟨W1, U^0.5⟩|² (normalized) = {:.6}",
        w1.hs_fidelity(&geo_half)
    );
    println!("\nSame coarse operator, two curves through operator space: the geodesic");
    println!("midpoint is a wide delocalized kernel; the causal path detours through");
    println!("the Fourier frame and stays narrow. Imputation has geometry, and the");
    println!("width cursor can choose its route.\n");

    // ---- 3. The width cursor: fine grain at the leading edge ---------------
    println!("=== 3. Dyadic zoom: increasingly fine grain toward the present ===\n");
    let mut cursor = WidthCursor::whole(&flow, 1.0);
    for level in 0..=4 {
        if level > 0 {
            let last = cursor.segments.len() - 1;
            cursor.refine(last);
        }
        let total_fid = cursor.total().hs_fidelity(&whole);
        let desc: Vec<String> = cursor
            .segments
            .iter()
            .map(|s| format!("[{:.4},{:.4})χ{}", s.t0, s.t1, s.op.max_bond_dim()))
            .collect();
        println!(
            "zoom {}: {} segments, total ≡ A to {:.9}",
            level,
            cursor.segments.len(),
            total_fid
        );
        println!("        {}", desc.join(" "));
    }
    println!("\nEvery refinement is imputed from the flow; the composition is invariant.");
    println!("The past stays coarse (one wide block), the present resolves finer and");
    println!("finer — the operation-width cursor in action.\n");

    // ---- 4. Walking the ring through the refined pipeline ------------------
    println!("=== 4. A state stepped through the refined pipeline ===\n");
    let mut walk = WidthCursor::whole(&flow, 1.0);
    walk.refine(0);
    walk.refine(0);
    walk.refine(2); // quarters
    let x0 = 1000usize;
    let start = Mps::basis_state(&profile, &digits(&profile, x0), SPEC);
    let traj = walk.trajectory(&start);

    println!(
        "{:>6} {:>12} {:>12} {:>10} {:>14}",
        "t", "circ. pos", "wobble", "|⟨ω^x⟩|", "participation"
    );
    let tau = 2.0 * PI;
    let nf = n_total as f64;
    for (k, state) in traj.iter().enumerate() {
        let t = k as f64 * 0.25;
        let d = state.to_dense();
        let mut moment = C64::ZERO;
        let mut part = 0.0;
        for (x, a) in d.amps.iter().enumerate() {
            moment += C64::cis(tau * x as f64 / nf).scale(a.abs2());
            part += a.abs2() * a.abs2();
        }
        let pos = (moment.arg() / tau * nf).rem_euclid(nf);
        let wobble_pred = -(tau * t).sin() / tau;
        println!(
            "{:>6.2} {:>12.4} {:>12.4} {:>10.6} {:>14.4}",
            t,
            pos,
            wobble_pred,
            moment.abs(),
            part
        );
    }
    println!("\nThe circular position follows x0 + t + wobble exactly, with the wobble");
    println!("-sin(2πt)/2π the ring's discreteness pushing back on the continuous");
    println!("flow — vanishing at integer t, where the state relocalizes: the cursor's");
    println!("finest observable grain agrees with arithmetic exactly when the");
    println!("operation width is whole.");
}

fn digits(profile: &[usize], mut x: usize) -> Vec<usize> {
    let mut d = vec![0usize; profile.len()];
    for (i, &dim) in profile.iter().enumerate().rev() {
        d[i] = x % dim;
        x /= dim;
    }
    d
}
