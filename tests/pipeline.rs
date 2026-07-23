//! Circuit-of-circuits verification: QFT, phase-ramp, and inverse-QFT blocks
//! composed as MPOs must collapse into the Fourier-space adder over the
//! diamond's ring, and adders must compose recursively into adders.

use novel_quantum_structures::c64::C64;
use novel_quantum_structures::dense::{total_dim, DenseState};
use novel_quantum_structures::mpo::Mpo;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{radix, zigzag, Rng, TruncSpec};

// Rank cap 64 never bites on this geometry (the diamond(2,4) carrier is
// capped at 36 by the profile itself) but keeps composition products in the
// randomized-SVD regime.
const SPEC: TruncSpec = TruncSpec {
    max_rank: 64,
    cutoff: 1e-12,
};

fn adder_pipeline(profile: &[usize], c_add: usize) -> Mpo {
    let qft = Mpo::from_circuit(&radix::mixed_radix_qft(profile, 0.0), SPEC);
    let ramp = Mpo::from_circuit(&radix::fourier_phase_ramp(profile, c_add), SPEC);
    // QFT† ∘ (D_c ∘ QFT): blocks composed as operators, not concatenated gates.
    qft.adjoint()
        .compose_after(&ramp.compose_after(&qft, SPEC), SPEC)
}

#[test]
fn composed_blocks_form_the_adder() {
    let profile = zigzag::diamond(2, 4); // [2,3,4,3,2], N = 144
    let n_total = total_dim(&profile);
    let c_add = 41;
    let adder = adder_pipeline(&profile, c_add);

    // 1. Against the gate-level compilation of the same pipeline.
    let reference = Mpo::from_circuit(&radix::adder_circuit(&profile, c_add), SPEC);
    let f = adder.hs_fidelity(&reference);
    assert!(
        (f - 1.0).abs() < 1e-8,
        "hs fidelity vs compiled circuit {}",
        f
    );

    // 2. Against the exact permutation on basis states.
    for x in [0usize, 1, 17, 99, 143] {
        let mut basis = DenseState::zero_state(&profile);
        basis.amps[0] = C64::ZERO;
        basis.amps[x] = C64::ONE;
        let psi = Mps::from_dense(&basis, SPEC);
        let shifted = adder.apply_to(&psi).to_dense();
        let target = (x + c_add) % n_total;
        assert!(
            (shifted.amps[target].abs() - 1.0).abs() < 1e-8,
            "x={}: weight at x+c {}",
            x,
            shifted.amps[target].abs()
        );
    }

    // 3. The composed permutation is a *simpler* operator than its Fourier
    //    factors: the adder's operator entanglement collapses to carry
    //    logic. (Measured: QFT block chi is an order of magnitude larger.)
    let qft = Mpo::from_circuit(&radix::mixed_radix_qft(&profile, 0.0), SPEC);
    assert!(
        adder.max_bond_dim() < qft.max_bond_dim(),
        "adder chi {} !< qft chi {}",
        adder.max_bond_dim(),
        qft.max_bond_dim()
    );
}

#[test]
fn adders_compose_recursively_into_adders() {
    let profile = zigzag::diamond(2, 4);
    let n_total = total_dim(&profile);
    let a3 = adder_pipeline(&profile, 3);
    let a7 = adder_pipeline(&profile, 7);
    let a10 = adder_pipeline(&profile, 10);

    // Block recursion: A7 ∘ A3 = A10, as operators.
    let composed = a7.compose_after(&a3, SPEC);
    let f = composed.hs_fidelity(&a10);
    assert!((f - 1.0).abs() < 1e-8, "A7∘A3 vs A10: hs fidelity {}", f);
    // The composite stays in the compact family (carry-logic bond dims).
    assert!(
        composed.max_bond_dim() <= a10.max_bond_dim() + 2,
        "composition left the compact family: {} vs {}",
        composed.max_bond_dim(),
        a10.max_bond_dim()
    );

    // Iterated self-composition: A3^k = A_{3k mod N}.
    let mut acc = Mpo::identity(&profile, SPEC);
    for _ in 0..5 {
        acc = a3.compose_after(&acc, SPEC);
    }
    let direct = adder_pipeline(&profile, 15 % n_total);
    let f2 = acc.hs_fidelity(&direct);
    assert!((f2 - 1.0).abs() < 1e-8, "A3^5 vs A15: hs fidelity {}", f2);
}

#[test]
fn conjugating_back_recovers_the_diagonal_block() {
    // QFT ∘ adder ∘ QFT† must return to the bond-dimension-1 phase ramp:
    // recursion runs in both directions through the Fourier frame.
    let profile = zigzag::diamond(2, 4);
    let c_add = 29;
    let adder = adder_pipeline(&profile, c_add);
    let qft = Mpo::from_circuit(&radix::mixed_radix_qft(&profile, 0.0), SPEC);
    let back = qft.compose_after(&adder.compose_after(&qft.adjoint(), SPEC), SPEC);
    let ramp = Mpo::from_circuit(&radix::fourier_phase_ramp(&profile, c_add), SPEC);
    let f = back.hs_fidelity(&ramp);
    assert!((f - 1.0).abs() < 1e-8, "recovered ramp: hs fidelity {}", f);
    assert_eq!(ramp.max_bond_dim(), 1, "ramp must be a product operator");
}

#[test]
fn one_mpo_application_equals_the_whole_block() {
    // The compiled QFT block acts on states as one operation, identical to
    // running its gates one by one.
    let profile = zigzag::diamond(2, 4);
    let circuit = radix::mixed_radix_qft(&profile, 0.0);
    let block = Mpo::from_circuit(&circuit, SPEC);

    let mut rng = Rng::new(7);
    let mut dense = DenseState::zero_state(&profile);
    for a in &mut dense.amps {
        *a = rng.c_gaussian();
    }
    dense.normalize();
    let psi = Mps::from_dense(&dense, SPEC);

    let via_block = block.apply_to(&psi);
    let mut via_gates = psi.clone();
    circuit.run_mps(&mut via_gates);
    let f = via_block.fidelity(&via_gates);
    assert!((f - 1.0).abs() < 1e-8, "block vs gates fidelity {}", f);
}
