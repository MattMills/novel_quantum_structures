//! Circuit-of-circuits: operations reified as pipeline states between the
//! past and future boundaries of the chain, composed and recursed as
//! first-class objects.
//!
//! The cast, all over the diamond's ring `Z_2880`:
//!
//! * `V`  — the mixed-radix QFT (reversed-digit convention): the structured
//!   **tail-radix phase web** compiled into an MPO;
//! * `D_c` — the Draper phase ramp, a bond-dimension-1 diagonal block;
//! * `A_c = V† ∘ D_c ∘ V` — three blocks composed *as operators*, which
//!   collapses into the modular adder `|x⟩ → |x+c mod N⟩`.
//!
//! Then recursion: adders composed with adders stay adders (a closed compact
//! family), and conjugating back through `V` recovers the diagonal block.
//! Finale: the same pipeline on a 25-site wave — exact `+c` arithmetic on a
//! ring of 8.6 trillion elements, the operator held in a few dozen KB.
//!
//! Run with: `cargo run --release --example circuit_of_circuits`

use novel_quantum_structures::dense::{total_dim, DenseState};
use novel_quantum_structures::mpo::Mpo;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{radix, zigzag, Rng, TruncSpec};
use std::time::Instant;

const SPEC: TruncSpec = TruncSpec {
    max_rank: 96,
    cutoff: 1e-12,
};

fn ent_rounded(e: &[f64]) -> Vec<f64> {
    e.iter().map(|x| (x * 100.0).round() / 100.0).collect()
}

fn main() {
    let profile = zigzag::diamond(2, 5);
    let n_total = total_dim(&profile);
    println!("=== I. The structured tail-radix phase web ===\n");
    println!("chain {:?} encodes Z_{}", profile, n_total);
    println!("output digit j couples ONLY to input digits i ≥ j (the tail),");
    println!("with angle 2π/D where D = d_j·d_(j+1)···d_i (radix segment product):\n");
    println!(
        "{:>6} {:>6} {:>10} {:>12}",
        "j", "i", "segment D", "angle (rad)"
    );
    for j in [0usize, 2, 3] {
        let mut seg = 1.0f64;
        for (i, &di) in profile.iter().enumerate().skip(j) {
            seg *= di as f64;
            println!(
                "{:>6} {:>6} {:>10} {:>12.2e}",
                j,
                i,
                seg,
                2.0 * std::f64::consts::PI / seg
            );
        }
    }
    println!("\nAngles decay super-exponentially with distance: the web is effectively");
    println!("banded, which is exactly why these blocks compress as operators.\n");

    // ---- II. Blocks as pipeline states -----------------------------------
    println!("=== II. Circuit blocks as pipeline states (past ⊗ future) ===\n");
    let t0 = Instant::now();
    let qft_circuit = radix::mixed_radix_qft_reversed(&profile, 0.0);
    let v = Mpo::from_circuit(&qft_circuit, SPEC);
    println!(
        "V = QFT over Z_2880: {} gates compiled in {:.0?} into ONE operator object:",
        qft_circuit.len(),
        t0.elapsed()
    );
    println!("  MPO bond dims          {:?}", v.bond_dims());
    println!(
        "  operator entanglement  {:?} bits per cut",
        ent_rounded(&v.operator_entanglement_bits())
    );
    println!(
        "  compile discard        {:.1e}   parameters {}",
        v.carrier.discarded_weight,
        v.param_count()
    );

    let c_add = 1234;
    let ramp = Mpo::from_circuit(&radix::fourier_phase_ramp_reversed(&profile, c_add), SPEC);
    println!(
        "\nD_{} ramp (diagonal block): bond dims {:?} — a product operator",
        c_add,
        ramp.bond_dims()
    );
    println!("\nEach block is literally a state on the doubled chain: site i carries a");
    println!("(future, past) leg pair, and the bond spectrum is the width of the");
    println!("pipeline of expectations between the two boundaries at that cut.\n");

    // ---- III. Composition: the pipeline collapses ------------------------
    println!("=== III. Composing blocks: gate composition at the next scale up ===\n");
    let t1 = Instant::now();
    let adder = v
        .adjoint()
        .compose_after(&ramp.compose_after(&v, SPEC), SPEC);
    println!(
        "A_{} = V† ∘ D_{} ∘ V   (two MPO compositions, {:.0?})",
        c_add,
        c_add,
        t1.elapsed()
    );
    println!("  QFT block bonds   {:?}", v.bond_dims());
    println!(
        "  adder bonds       {:?}  ← the pipeline NARROWS under composition",
        adder.bond_dims()
    );
    println!(
        "  adder opEE        {:?} bits (carry logic, ≤ 2 bits per cut)",
        ent_rounded(&adder.operator_entanglement_bits())
    );

    // Verify against the exact cyclic shift, dense.
    let mut rng = Rng::new(2026);
    let mut dense = DenseState::zero_state(&profile);
    for a in &mut dense.amps {
        *a = rng.c_gaussian();
    }
    dense.normalize();
    let psi = Mps::from_dense(&dense, SPEC);
    let shifted = adder.apply_to(&psi).to_dense();
    let mut expect = DenseState::zero_state(&profile);
    for x in 0..n_total {
        expect.amps[(x + c_add) % n_total] = dense.amps[x];
    }
    println!(
        "  verification      A_{}|ψ⟩ vs exact shift on a random state: fidelity {:.12}",
        c_add,
        shifted.fidelity(&expect)
    );

    // The contrast: standard-order QFT pays for its digit reversal.
    let p4 = zigzag::diamond(2, 4);
    let v4_rev = Mpo::from_circuit(&radix::mixed_radix_qft_reversed(&p4, 0.0), SPEC);
    let t2 = Instant::now();
    let v4_std = Mpo::from_circuit(&radix::mixed_radix_qft(&p4, 0.0), SPEC);
    println!("\nWhere does operator entanglement actually live? On diamond(2,4):");
    println!(
        "  QFT core (reversed digits): max chi {:>3}",
        v4_rev.max_bond_dim()
    );
    println!(
        "  QFT + mirror-swap reversal: max chi {:>3}  ({:.1?})",
        v4_std.max_bond_dim(),
        t2.elapsed()
    );
    println!("The Fourier core is cheap; the digit reversal — the bowtie permutation");
    println!("across the waist — carries the entanglement. In the adder pipeline the");
    println!("reversal cancels between V and V†, so it is never paid at all.\n");

    // ---- IV. Recursion ----------------------------------------------------
    println!("=== IV. Recursion: computing between operations ===\n");
    let a3 = {
        let r3 = Mpo::from_circuit(&radix::fourier_phase_ramp_reversed(&profile, 3), SPEC);
        v.adjoint().compose_after(&r3.compose_after(&v, SPEC), SPEC)
    };
    let mut acc = Mpo::identity(&profile, SPEC);
    print!("A_3 self-composed: chi per step:");
    for k in 1..=6 {
        acc = a3.compose_after(&acc, SPEC);
        print!(" {}", acc.max_bond_dim());
        let direct = {
            let rk = Mpo::from_circuit(&radix::fourier_phase_ramp_reversed(&profile, 3 * k), SPEC);
            v.adjoint().compose_after(&rk.compose_after(&v, SPEC), SPEC)
        };
        assert!((acc.hs_fidelity(&direct) - 1.0).abs() < 1e-8);
    }
    println!("  (each step verified = A_(3k): a closed compact family)");

    let back = v.compose_after(&adder.compose_after(&v.adjoint(), SPEC), SPEC);
    println!(
        "V ∘ A_{} ∘ V† bonds: {:?} — conjugating back through the Fourier frame",
        c_add,
        back.bond_dims()
    );
    println!(
        "recovers the diagonal ramp (hs fidelity {:.10}): recursion runs both ways.\n",
        back.hs_fidelity(&ramp)
    );

    // ---- V. Scale, and the second budget ----------------------------------
    println!("=== V. An 8.6-trillion-element ring, and the two budgets ===\n");
    let wave = zigzag::wave(2, 5, 4);
    let big_n: u64 = wave.iter().map(|&d| d as u64).product();
    println!(
        "wave(2,5,4): {} sites, N = {} ≈ {:.2e}",
        wave.len(),
        big_n,
        big_n as f64
    );
    let t3 = Instant::now();
    let vw = Mpo::from_circuit(&radix::mixed_radix_qft_reversed(&wave, 0.0), SPEC);
    let rw = Mpo::from_circuit(&radix::fourier_phase_ramp_reversed(&wave, 1), SPEC);
    let aw = vw
        .adjoint()
        .compose_after(&rw.compose_after(&vw, SPEC), SPEC);
    println!(
        "\nFourier-pipeline A_1 over Z_{} built in {:.1?}: max chi {}, ≈ {:.0} KB",
        big_n,
        t3.elapsed(),
        aw.max_bond_dim(),
        aw.param_count() as f64 * 16.0 / 1e3
    );
    // The collapsed form, built directly: a bond-2 carry-propagation MPO.
    let exact = radix::adder_mpo_exact(&wave, 1, SPEC);
    println!(
        "exact carry MPO for the same operator: max chi {}, ≈ {:.1} KB",
        exact.max_bond_dim(),
        exact.param_count() as f64 * 16.0 / 1e3
    );
    println!(
        "operator-level agreement: hs fidelity {:.9}",
        aw.hs_fidelity(&exact)
    );

    let digits_of = |mut x: u64| -> Vec<usize> {
        let mut d = vec![0usize; wave.len()];
        for (i, &dim) in wave.iter().enumerate().rev() {
            d[i] = (x % dim as u64) as usize;
            x /= dim as u64;
        }
        d
    };
    println!(
        "\n{:>24} {:>18} {:>14}",
        "input", "Fourier pipeline", "carry MPO"
    );
    for &x in &[0u64, 123_456_789_012, big_n - 1] {
        let input = Mps::basis_state(&wave, &digits_of(x), SPEC);
        let expect = Mps::basis_state(&wave, &digits_of((x + 1) % big_n), SPEC);
        let f_pipeline = aw.apply_to(&input).fidelity(&expect);
        let f_exact = exact.apply_to(&input).fidelity(&expect);
        println!(
            "{:>24} {:>18.9} {:>14.9}",
            format!("|{}⟩", x),
            f_pipeline,
            f_exact
        );
    }
    println!("\nThe last row is the finding. Composition is governed by TWO budgets:");
    println!("entanglement (bond dimension) and spectral weight (precision). The");
    println!("|N-1⟩ → |0⟩ carry rides tail-radix couplings of angle 2π/N ≈ 7e-13 and");
    println!("occupies relative operator weight ~1/N ≈ 1e-13 — below any f64 cutoff —");
    println!("so the Fourier-composed pipeline loses exactly that one column of the");
    println!("operator while remaining machine-accurate on typical inputs. The");
    println!("collapsed carry form needs no weight at all: bond dimension 2, exact at");
    println!("any ring size, provided it is compressed by RANK, never by weight.");
    println!("Average-case operator fidelity is not worst-case operator fidelity.");
}
