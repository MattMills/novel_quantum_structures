//! Sub-simulation of binary systems inside the zigzag chain: a 10-qubit
//! circuit compiled onto the 7-site diamond must reproduce the pure-qubit
//! reference exactly, on both backends.

use novel_quantum_structures::dense::DenseState;
use novel_quantum_structures::embed::{qubit_bit_swap, HostRegister, QubitOp};
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{gates, zigzag, Rng, TruncSpec};

fn random_qubit_circuit(num_qubits: usize, ops: usize, rng: &mut Rng) -> Vec<QubitOp> {
    let mut out = Vec::new();
    for step in 0..ops {
        if step % 2 == 0 {
            out.push(QubitOp::One(rng.below(num_qubits), gates::random(2, rng)));
        } else {
            let a = rng.below(num_qubits);
            let mut b = rng.below(num_qubits);
            while b == a {
                b = rng.below(num_qubits);
            }
            out.push(QubitOp::Two(a, b, gates::random(4, rng)));
        }
    }
    out
}

#[test]
fn diamond_hosts_ten_qubits_exactly() {
    let profile = zigzag::diamond(2, 5); // capacities 1,1,2,2,2,1,1
    let reg = HostRegister::new(&profile);
    assert_eq!(reg.num_qubits(), 10);

    let mut rng = Rng::new(99);
    let qops = random_qubit_circuit(10, 40, &mut rng);

    // Reference: pure qubit chain.
    let mut reference = DenseState::zero_state(&[2; 10]);
    reg.compile_reference(&qops).run_dense(&mut reference);

    // Hosted: compiled onto the diamond, dense backend.
    let mut hosted = DenseState::zero_state(&profile);
    reg.compile(&qops).run_dense(&mut hosted);

    let (qv, leakage) = reg.extract_qubit_state(&hosted);
    assert!(leakage < 1e-18, "leakage {}", leakage);
    let overlap: f64 = qv
        .iter()
        .zip(reference.amps.iter())
        .map(|(a, b)| (a.conj() * *b).re)
        .sum();
    // Global phase is preserved by construction, so amplitudes match 1:1.
    let diff: f64 = qv
        .iter()
        .zip(reference.amps.iter())
        .map(|(a, b)| (*a - *b).abs())
        .fold(0.0, f64::max);
    assert!(
        diff < 1e-9,
        "max amplitude diff {} (overlap {})",
        diff,
        overlap
    );
}

#[test]
fn hosted_circuit_runs_on_mps_backend_too() {
    let profile = zigzag::diamond(2, 5);
    let reg = HostRegister::new(&profile);
    let mut rng = Rng::new(123);
    let qops = random_qubit_circuit(10, 25, &mut rng);

    let mut reference = DenseState::zero_state(&[2; 10]);
    reg.compile_reference(&qops).run_dense(&mut reference);

    let mut psi = Mps::zero_state(&profile, TruncSpec::exact());
    reg.compile(&qops).run_mps(&mut psi);

    let (qv, leakage) = reg.extract_qubit_state(&psi.to_dense());
    assert!(leakage < 1e-14, "leakage {}", leakage);
    let diff: f64 = qv
        .iter()
        .zip(reference.amps.iter())
        .map(|(a, b)| (*a - *b).abs())
        .fold(0.0, f64::max);
    assert!(diff < 1e-7, "max amplitude diff {}", diff);
}

#[test]
fn edge_qubit_shuttles_into_the_waist_and_back() {
    // Prepare a random qubit state on edge site 0, swap it into bit 0 of the
    // waist's 2-qubit register, act there, swap back.
    let profile = zigzag::diamond(2, 5);
    let waist = zigzag::waist_positions(&profile)[0];
    let mut rng = Rng::new(55);
    let u = gates::random(2, &mut rng);

    let mut psi = DenseState::zero_state(&profile);
    psi.apply1(0, &u); // edge qubit now in a random state

    let swap_in = qubit_bit_swap(profile[waist], 0); // (qubit, host) pair order
    psi.apply2(0, waist, &swap_in);

    // The edge site must now be |0⟩ ...
    assert!((psi.site_probability(0, 0) - 1.0).abs() < 1e-12);
    // ... and the waist register bit 0 must carry the state: levels {2,3}
    // (bit0 = 1) hold |u_1|^2 of the weight.
    let p_bit1 = psi.site_probability(waist, 2) + psi.site_probability(waist, 3);
    assert!((p_bit1 - u.at(1, 0).abs2()).abs() < 1e-12);

    // Round trip restores the original state exactly.
    psi.apply2(0, waist, &swap_in);
    let mut expect = DenseState::zero_state(&profile);
    expect.apply1(0, &u);
    assert!((psi.fidelity(&expect) - 1.0).abs() < 1e-12);
}

#[test]
fn coarse_fourier_on_merged_pair_is_the_two_qubit_qft() {
    // F_4 on a 4-level site whose levels encode two big-endian bits equals
    // the textbook 2-qubit QFT circuit: H(q0), CS, H(q1), SWAP.
    let f4 = gates::fourier(4);
    let h = gates::hadamard();
    let cs = gates::phase_diag(&[0.0, 0.0, 0.0, std::f64::consts::FRAC_PI_2]);
    let swap = gates::swap(2);
    let id = gates::identity(2);
    // Circuit order: first applied is rightmost in the product.
    let qft = swap.mul(&id.kron(&h)).mul(&cs).mul(&h.kron(&id));
    assert!(
        f4.max_abs_diff(&qft) < 1e-12,
        "F_4 != 2-qubit QFT circuit (diff {})",
        f4.max_abs_diff(&qft)
    );
}
