//! Cross-validation of the MPS backend against exact dense simulation on the
//! flagship structure: the 2→3→4→5→4→3→2 diamond.

use novel_quantum_structures::dense::DenseState;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{gates, zigzag, Rng, TruncSpec};

#[test]
fn bowtie_entropies_match_theory() {
    let profile = zigzag::diamond(2, 5);
    let mut psi = Mps::zero_state(&profile, TruncSpec::exact());
    zigzag::bowtie(&profile).run_mps(&mut psi);

    // Mirror Bell pairs: (0,6) with d=2, (1,5) with d=3, (2,4) with d=4;
    // waist site 3 stays unentangled. Bond k carries every pair it cuts.
    let log3 = 3f64.log2();
    let expected = [
        1.0,              // pair (0,6)
        1.0 + log3,       // + (1,5)
        1.0 + log3 + 2.0, // + (2,4) = log2 24
        1.0 + log3 + 2.0,
        1.0 + log3,
        1.0,
    ];
    let got = psi.bond_entropies_bits();
    for (b, (g, e)) in got.iter().zip(expected.iter()).enumerate() {
        assert!((g - e).abs() < 1e-9, "bond {}: got {} expected {}", b, g, e);
    }
    // Maximal possible Schmidt rank across the central bond is
    // min(2·3·4, 5·4·3·2) = 24 — the bowtie state saturates it.
    assert_eq!(psi.max_bond_dim(), 24);

    // And the whole thing matches the dense ground truth.
    let mut reference = DenseState::zero_state(&profile);
    zigzag::bowtie(&profile).run_dense(&mut reference);
    let f = psi.to_dense().fidelity(&reference);
    assert!((f - 1.0).abs() < 1e-9, "fidelity {}", f);
}

#[test]
fn crosswave_rounds_agree_with_dense() {
    let profile = zigzag::diamond(2, 5);
    let mut psi = Mps::zero_state(&profile, TruncSpec::exact());
    let mut reference = DenseState::zero_state(&profile);
    for _ in 0..2 {
        let c = zigzag::crosswave_round(&profile);
        c.run_mps(&mut psi);
        c.run_dense(&mut reference);
    }
    let f = psi.to_dense().fidelity(&reference);
    assert!((f - 1.0).abs() < 1e-8, "fidelity {}", f);
    assert!(psi.discarded_weight < 1e-16);
}

#[test]
fn random_cross_scale_circuit_agrees_with_dense() {
    let profile = zigzag::diamond(2, 5);
    let mut rng = Rng::new(2026);
    let mut psi = Mps::zero_state(&profile, TruncSpec::exact());
    let mut reference = DenseState::zero_state(&profile);

    // Alternate brickwork layers with random long-range mirror/waist gates.
    let mut c = zigzag::brickwork_random(&profile, 1, &mut rng);
    for (i, j) in zigzag::mirror_pairs(&profile) {
        c.two(i, j, gates::random(profile[i] * profile[j], &mut rng));
    }
    let waist = zigzag::waist_positions(&profile)[0];
    for v in zigzag::valley_positions(&profile) {
        c.two(
            v,
            waist,
            gates::random(profile[v] * profile[waist], &mut rng),
        );
    }
    c.run_mps(&mut psi);
    c.run_dense(&mut reference);

    assert!((psi.norm() - 1.0).abs() < 1e-8);
    let f = psi.to_dense().fidelity(&reference);
    assert!((f - 1.0).abs() < 1e-7, "fidelity {}", f);
}

#[test]
fn truncated_simulation_degrades_gracefully() {
    let profile = zigzag::diamond(2, 5);
    let mut rng = Rng::new(7);
    let c = zigzag::brickwork_random(&profile, 3, &mut rng);

    let mut reference = DenseState::zero_state(&profile);
    c.run_dense(&mut reference);

    let mut last_infidelity = f64::INFINITY;
    for max_bond in [4usize, 8, 16, 24] {
        let mut psi = Mps::zero_state(&profile, TruncSpec::new(max_bond, 0.0));
        c.run_mps(&mut psi);
        assert!((psi.norm() - 1.0).abs() < 1e-8);
        let infidelity = 1.0 - psi.to_dense().fidelity(&reference);
        // More bond dimension must not make things significantly worse.
        assert!(
            infidelity <= last_infidelity + 1e-9,
            "chi {}: infidelity {} worse than previous {}",
            max_bond,
            infidelity,
            last_infidelity
        );
        last_infidelity = infidelity;
    }
    // At chi = 24 the diamond is exactly representable.
    assert!(
        last_infidelity < 1e-7,
        "infidelity at full rank {}",
        last_infidelity
    );
}
