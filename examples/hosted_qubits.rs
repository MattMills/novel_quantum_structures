//! Sub-simulating a binary system inside the zigzag chain.
//!
//! The 7-site diamond `2,3,4,5,4,3,2` hosts 10 logical qubits
//! (capacities 1,1,2,2,2,1,1): the valley sites *are* qubits — the binary
//! set on either side — and the interior sites carry 1–2 logical qubits each
//! plus slack. A 10-qubit GHZ circuit compiled onto the diamond reproduces
//! the pure-qubit reference exactly, and an edge qubit can shuttle *into*
//! the waist's register and back — the bidirectional pair in action.
//!
//! Run with: `cargo run --release --example hosted_qubits`

use novel_quantum_structures::dense::DenseState;
use novel_quantum_structures::embed::{capacity, qubit_bit_swap, HostRegister, QubitOp};
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{gates, zigzag, TruncSpec};

fn main() {
    let profile = zigzag::diamond(2, 5);
    let reg = HostRegister::new(&profile);

    println!("=== Hosting a 10-qubit register in the 7-site diamond ===\n");
    println!(
        "{:>6} {:>6} {:>10} {:>8}",
        "site", "dim", "capacity", "qubits"
    );
    let mut q0 = 0;
    for (site, &d) in profile.iter().enumerate() {
        let k = capacity(d);
        let names: Vec<String> = (q0..q0 + k).map(|q| format!("q{}", q)).collect();
        println!("{:>6} {:>6} {:>10} {:>8}", site, d, k, names.join(","));
        q0 += k;
    }
    println!(
        "\ntotal: {} qubits in a {}-dimensional chain (2^10 = 1024; overhead ×{:.2})\n",
        reg.num_qubits(),
        zigzag::dense_dimension(&profile),
        zigzag::dense_dimension(&profile) / 1024.0
    );

    // ---- GHZ-10 through the diamond --------------------------------------
    let mut qops = vec![QubitOp::One(0, gates::hadamard())];
    for q in 0..9 {
        qops.push(QubitOp::Two(q, q + 1, gates::cshift(2, 2)));
    }

    let mut reference = DenseState::zero_state(&[2; 10]);
    reg.compile_reference(&qops).run_dense(&mut reference);

    // On the diamond — MPS backend.
    let mut psi = Mps::zero_state(&profile, TruncSpec::exact());
    reg.compile(&qops).run_mps(&mut psi);

    let (qv, leakage) = reg.extract_qubit_state(&psi.to_dense());
    let diff = qv
        .iter()
        .zip(reference.amps.iter())
        .map(|(a, b)| (*a - *b).abs())
        .fold(0.0, f64::max);
    println!("=== GHZ-10 compiled onto the diamond (MPS backend) ===\n");
    println!(
        "|amp(|0…0⟩)| = {:.6}, |amp(|1…1⟩)| = {:.6}",
        qv[0].abs(),
        qv[1023].abs()
    );
    println!("max |amplitude - reference| : {:.3e}", diff);
    println!("slack-subspace leakage      : {:.3e}", leakage);
    println!("MPS bond dims on the diamond: {:?}", psi.bond_dims());
    println!("(a 10-qubit GHZ costs chi = 2 — the wave hosts it with room to spare)\n");

    // ---- Shuttle: edge qubit → waist register → back ---------------------
    println!("=== Bidirectional pair: edge qubit shuttling through the waist ===\n");
    let waist = zigzag::waist_positions(&profile)[0];
    let d_waist = profile[waist];

    // Random state on the two edge qubits (sites 0 and 6).
    let mut rng = novel_quantum_structures::Rng::new(42);
    let ua = gates::random(2, &mut rng);
    let ub = gates::random(2, &mut rng);
    let mut chain = DenseState::zero_state(&profile);
    chain.apply1(0, &ua);
    chain.apply1(6, &ub);

    // Swap both edge qubits into the waist's 2-bit register: site 0 → bit 0,
    // site 6 → bit 1. The 5-dit now holds both edge qubits; its slack level
    // rides along untouched.
    chain.apply2(0, waist, &qubit_bit_swap(d_waist, 0));
    chain.apply2(6, waist, &qubit_bit_swap(d_waist, 1));
    println!(
        "after shuttling in: P(site0 = |0⟩) = {:.12}",
        chain.site_probability(0, 0)
    );
    println!(
        "                    P(site6 = |0⟩) = {:.12}",
        chain.site_probability(6, 0)
    );

    // One coarse gate on the waist = a joint two-qubit operation: F_4 on the
    // register is the full 2-qubit QFT (verified in the test suite).
    let f4_lifted =
        novel_quantum_structures::embed::lift_register_gate(&gates::fourier(4), d_waist);
    chain.apply1(waist, &f4_lifted);

    // Shuttle back out.
    chain.apply2(6, waist, &qubit_bit_swap(d_waist, 1));
    chain.apply2(0, waist, &qubit_bit_swap(d_waist, 0));

    // Reference: apply the 2-qubit QFT directly on the edge pair (0,6).
    let reg2 = HostRegister::new(&profile);
    let mut expect = DenseState::zero_state(&profile);
    expect.apply1(0, &ua);
    expect.apply1(6, &ub);
    let qft_ops = vec![QubitOp::Two(0, 9, gates::fourier(4))]; // q0 = site0, q9 = site6
    reg2.compile(&qft_ops).run_dense(&mut expect);

    println!(
        "\nround trip (in → coarse F_4 on the waist → out) vs direct 2-qubit QFT on the edges:"
    );
    println!("fidelity = {:.15}", chain.fidelity(&expect));
    println!("\nThe edge qubits and the waist register are bidirectional pairs: binary");
    println!("states move up into the high-dimensional site, get processed by single");
    println!("coarse gates, and return — 'acting at multiple scales' made literal.");
}
