//! The expand/collapse machinery: the *representation* changes scale while
//! the state is preserved.
//!
//! An 8-qubit chain is entangled, then merged pairwise into a 4-dit chain
//! and again into a 16-dit pair — the same state at three different scales,
//! with identical entanglement profiles. At the coarse scale, one gate
//! (`F_4`) executes what takes a multi-gate subcircuit at the fine scale;
//! splitting back down recovers the fine-grained chain exactly.
//!
//! Run with: `cargo run --release --example scale_morphing`

use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{gates, zigzag, Rng, TruncSpec};

fn entropies_rounded(psi: &mut Mps) -> Vec<f64> {
    psi.bond_entropies_bits()
        .iter()
        .map(|e| (e * 1000.0).round() / 1000.0)
        .collect()
}

fn main() {
    let mut rng = Rng::new(7);

    // ---- One state, three scales -----------------------------------------
    println!("=== One state, three scales ===\n");
    let fine_dims = vec![2usize; 8];
    let mut psi = Mps::zero_state(&fine_dims, TruncSpec::exact());
    zigzag::brickwork_random(&fine_dims, 2, &mut rng).run_mps(&mut psi);
    let original = psi.clone();

    println!("scale 0 (8 qubits)   dims {:?}", psi.dims);
    println!("  bond entropies {:?}", entropies_rounded(&mut psi));

    for i in 0..4 {
        psi.merge_sites(i); // merge (i, i+1) repeatedly: [2,2,2,2,2,2,2,2] → [4,4,4,4]
    }
    println!("scale 1 (4 ququarts) dims {:?}", psi.dims);
    println!("  bond entropies {:?}", entropies_rounded(&mut psi));

    psi.merge_sites(0);
    psi.merge_sites(1);
    println!("scale 2 (2 sixteen-dits) dims {:?}", psi.dims);
    println!("  bond entropies {:?}", entropies_rounded(&mut psi));
    println!("\nMerging is exact: interior entanglement becomes intra-site structure;");
    println!("the entropies that survive are the same numbers from the fine scale.\n");

    // Split all the way back down and verify.
    psi.split_site(0, 4, 4);
    psi.split_site(2, 4, 4);
    for i in (0..4).rev() {
        psi.split_site(i, 2, 2);
    }
    println!(
        "round trip back to 8 qubits: fidelity vs original = {:.15}\n",
        psi.fidelity(&original)
    );

    // ---- One coarse gate = a fine subcircuit ------------------------------
    println!("=== One coarse gate = a fine-scale subcircuit ===\n");
    // Fine route: the textbook 2-qubit QFT circuit on qubits (0,1):
    // H(q0), CS(q0,q1), H(q1), SWAP — four gates.
    let dims2 = vec![2usize, 2];
    let mut fine = Mps::zero_state(&dims2, TruncSpec::exact());
    zigzag::brickwork_random(&dims2, 1, &mut rng).run_mps(&mut fine);
    let mut coarse = fine.clone();

    let cs = gates::phase_diag(&[0.0, 0.0, 0.0, std::f64::consts::FRAC_PI_2]);
    fine.apply1(0, &gates::hadamard());
    fine.apply2(0, 1, &cs);
    fine.apply1(1, &gates::hadamard());
    fine.apply2(0, 1, &gates::swap(2));

    // Coarse route: collapse the pair into one 4-level site, apply the
    // single gate F_4, expand again — one entangling operation.
    coarse.merge_sites(0);
    coarse.apply1(0, &gates::fourier(4));
    coarse.split_site(0, 2, 2);

    println!("4-gate QFT circuit at qubit scale  vs  merge → F_4 → split:");
    println!("fidelity = {:.15}", coarse.fidelity(&fine));

    // ---- Climbing the dimension ladder ------------------------------------
    println!("\n=== Climbing the zigzag dimension ladder in place ===\n");
    // promote: qutrit → 5-dit (exact embedding); demote back (leakage 0 if
    // the extra levels were never populated).
    let mut ladder = Mps::zero_state(&[3, 3], TruncSpec::exact());
    ladder.apply1(0, &gates::fourier(3));
    ladder.apply2(0, 1, &gates::cshift(3, 3));
    let before = ladder.clone();
    ladder.promote_site(0, 5);
    ladder.promote_site(1, 4);
    println!(
        "promoted dims {:?} — gates of the larger algebra now apply",
        ladder.dims
    );
    // Act with a 5-dit unitary that preserves the populated subspace: the
    // clock gate only adds phases per level.
    ladder.apply1(0, &gates::clock(5));
    ladder.apply1(0, &gates::clock(5).adjoint()); // undo
    let leak0 = ladder.demote_site(0, 3);
    let leak1 = ladder.demote_site(1, 3);
    println!(
        "demoted back to {:?}: leakage ({:.1e}, {:.1e}), fidelity vs original {:.15}",
        ladder.dims,
        leak0,
        leak1,
        ladder.fidelity(&before)
    );
}
