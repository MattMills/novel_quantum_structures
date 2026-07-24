//! The **operative window**: how much computation — gates, circuits,
//! time-steps — is folded into one operator and operated on at once, and how
//! the cost scales as that window grows.
//!
//! Finding 28 measured `χ` as a resource across *space* (operator
//! entanglement across a cut). This measures the dual axis — *time*: the
//! window `W` is the temporal aperture (how much computation the operator
//! represents simultaneously), and its cost is again operator entanglement,
//! now the width of the past↔future pipeline the whole window must sustain.
//!
//! The result is the structure/scrambling boundary, transposed onto time:
//!
//! * **Arithmetic windows are unbounded.** A reversible/classical
//!   computation — adders, multipliers, any permutation — folds into one
//!   operator of bounded `χ` *no matter how deep*. You can hold an
//!   arbitrarily long arithmetic computation in one narrow operator.
//! * **Generic quantum windows are capped.** A generic circuit's operator
//!   entanglement grows with depth and saturates the geometric budget after
//!   only a layer or two — the aperture limit.
//! * **Operating on the whole window at once pays.** Folding `W` gates into
//!   one operator and applying it in a single MPO contraction beats running
//!   the `W` gates sequentially, and amortizes over every reuse.
//!
//! And the two-sided reading (dual to finding 28): the aperture limit *is*
//! the classical-simulability boundary on the time axis — as long as the
//! window stays narrow, the whole temporal block is classically held and
//! applied; when operator entanglement saturates, it is not.
//!
//! Run with: `cargo run --release --example operative_window`

use novel_quantum_structures::circuit::Op;
use novel_quantum_structures::dense::DenseState;
use novel_quantum_structures::mpo::Mpo;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::cascade::Transducer;
use novel_quantum_structures::{gates, radix, zigzag, Rng, TruncSpec};
use std::time::Instant;

const SPEC: TruncSpec = TruncSpec {
    max_rank: 256,
    cutoff: 1e-12,
};

fn opee_max(m: &Mpo) -> f64 {
    m.operator_entanglement_bits()
        .iter()
        .cloned()
        .fold(0.0, f64::max)
}

/// Fold `layers` of Haar-random brickwork into one operator, in place.
fn fold_random_layer(op: &mut Mpo, profile: &[usize], rng: &mut Rng) {
    for i in (0..profile.len() - 1).step_by(2) {
        op.absorb_after(&Op::Two(i, i + 1, gates::random(profile[i] * profile[i + 1], rng)));
    }
    for i in (1..profile.len() - 1).step_by(2) {
        op.absorb_after(&Op::Two(i, i + 1, gates::random(profile[i] * profile[i + 1], rng)));
    }
}

fn main() {
    let p = zigzag::diamond(2, 4); // N = 144
    println!("diamond {:?} — the operative window folded into one operator\n", p);

    // ---- 1. The two regimes -----------------------------------------------
    println!("=== 1. Two regimes: arithmetic vs generic, folded into one operator ===\n");
    let v = Mpo::from_circuit(&radix::mixed_radix_qft_reversed(&p, 0.0), SPEC);
    let adder = |c: usize| {
        v.adjoint().compose_after(
            &Mpo::from_circuit(&radix::fourier_phase_ramp_reversed(&p, c), SPEC).compose_after(&v, SPEC),
            SPEC,
        )
    };

    println!("ARITHMETIC window — compose K adders (a K-deep additive computation):");
    println!("{:>6} {:>10} {:>18}", "K adds", "folded χ", "op-entangle (bits)");
    let mut acc = Mpo::identity(&p, SPEC);
    for k in 1..=16 {
        acc = adder(7).compose_after(&acc, SPEC);
        if [1usize, 2, 4, 8, 16].contains(&k) {
            println!("{:>6} {:>10} {:>18.2}", k, acc.max_bond_dim(), opee_max(&acc));
        }
    }
    println!("→ χ flat at 2 however deep: an arbitrary-length arithmetic window is FREE.\n");

    println!("MULTIPLICATIVE window — compose K ×7 (deep modular exponentiation ×7^K):");
    println!("{:>6} {:>10} {:>18}", "K mults", "folded χ", "op-entangle (bits)");
    let m7 = Transducer::mult(&p, 7).to_mpo(SPEC);
    let mut mk = Mpo::identity(&p, SPEC);
    for k in 1..=8 {
        mk = m7.compose_after(&mk, SPEC);
        if [1usize, 2, 4, 8].contains(&k) {
            println!("{:>6} {:>10} {:>18.2}", k, mk.max_bond_dim(), opee_max(&mk));
        }
    }
    println!("→ bounded by geometry (the resonant width), never the depth.\n");

    println!("GENERIC quantum window — fold D layers of Haar-random brickwork:");
    println!("{:>6} {:>10} {:>18}", "depth D", "folded χ", "op-entangle (bits)");
    let mut rng = Rng::new(3);
    let mut op = Mpo::identity(&p, SPEC);
    for d in 1..=8 {
        fold_random_layer(&mut op, &p, &mut rng);
        if [1usize, 2, 3, 4, 8].contains(&d) {
            println!("{:>6} {:>10} {:>18.2}", d, op.max_bond_dim(), opee_max(&op));
        }
    }
    println!("→ operator entanglement grows with depth and SATURATES: the aperture limit.\n");

    // ---- 2. Window capacity ------------------------------------------------
    println!("=== 2. Window capacity — depth before saturation, by geometry ===\n");
    println!("How many generic layers fit before operator entanglement fills the budget:\n");
    println!("{:>14} {:>12} {:>16} {:>14}", "profile", "χ cap", "opEE cap (bits)", "capacity (D)");
    for hi in [3usize, 4] {
        let prof = zigzag::diamond(2, hi);
        let mut rng = Rng::new(11);
        let mut o = Mpo::identity(&prof, SPEC);
        let mut prev = 0.0;
        let mut capacity = 0;
        let mut cap_chi = 1;
        let mut cap_ee = 0.0;
        for d in 1..=10 {
            fold_random_layer(&mut o, &prof, &mut rng);
            let ee = opee_max(&o);
            if ee - prev < 0.05 {
                capacity = d;
                cap_chi = o.max_bond_dim();
                cap_ee = ee;
                break;
            }
            prev = ee;
            capacity = d;
            cap_chi = o.max_bond_dim();
            cap_ee = ee;
        }
        println!(
            "{:>14} {:>12} {:>16.2} {:>14}",
            format!("diamond(2,{})", hi),
            cap_chi,
            cap_ee,
            capacity
        );
    }
    println!("\n→ a wider geometry holds a deeper generic window; an arithmetic window is");
    println!("  unbounded on any of them (χ never grows with depth).\n");

    // ---- 3. Operating on the whole window at once --------------------------
    println!("=== 3. The payoff — one contraction over the whole window vs sequential ===\n");
    let window = adder(50); // a folded QFT–ramp–QFT window
    let circuit = radix::adder_circuit_reversal_free(&p, 50);
    let w_gates = circuit.len();
    let mut dense = DenseState::zero_state(&p);
    for a in &mut dense.amps {
        *a = rng.c_gaussian();
    }
    dense.normalize();
    let psi = Mps::from_dense(&dense, SPEC);

    for reps in [1usize, 20, 100] {
        let t0 = Instant::now();
        for _ in 0..reps {
            let _ = window.apply_to(&psi);
        }
        let folded = t0.elapsed();
        let t1 = Instant::now();
        for _ in 0..reps {
            let mut s = psi.clone();
            circuit.run_mps(&mut s);
        }
        let seq = t1.elapsed();
        println!(
            "  window={} gates, applied {:>3}×:  folded {:>8.1?}  sequential {:>8.1?}  ({:.0}× faster)",
            w_gates,
            reps,
            folded,
            seq,
            seq.as_secs_f64() / folded.as_secs_f64()
        );
    }
    println!("\n→ operating on the whole window at once (one MPO contraction) is an order");
    println!("  of magnitude faster than replaying its gates, and amortizes over reuse.\n");

    // ---- 4. Synthesis ------------------------------------------------------
    println!("=== The operative window is a temporal aperture ===\n");
    println!("Its cost is operator entanglement — the width of the past↔future pipeline");
    println!("the whole folded window must sustain. Arithmetic (reversible) computation");
    println!("holds an UNBOUNDED aperture in one narrow operator; generic quantum");
    println!("computation is aperture-limited, its window capped by operator-entanglement");
    println!("saturation. That cap is the classical-simulability boundary read on the time");
    println!("axis — dual to finding 28's space reading: as long as the aperture stays");
    println!("narrow, an arbitrarily long block of computation is held, composed, and");
    println!("applied as a single object. χ prices the window in time exactly as it prices");
    println!("the operator in space.");
}
