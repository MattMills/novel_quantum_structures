//! The flagship structure: the 2→3→4→5→4→3→2 diamond with its mirror pairs
//! Bell-coupled across the waist (the "bowtie" state), cross-validated
//! against exact dense simulation, plus a bond-dimension budget study.
//!
//! Run with: `cargo run --release --example diamond_bowtie`

use novel_quantum_structures::dense::DenseState;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{zigzag, Rng, TruncSpec};

fn main() {
    let profile = zigzag::diamond(2, 5);
    println!("The dimension diamond:\n");
    println!("{}", zigzag::ascii_profile(&profile));
    println!(
        "mirror pairs (bidirectional, equal dimension): {:?}",
        zigzag::mirror_pairs(&profile)
    );
    println!(
        "waist: {:?}   valleys (binary edges): {:?}\n",
        zigzag::waist_positions(&profile),
        zigzag::valley_positions(&profile)
    );

    // ---- The bowtie state -------------------------------------------------
    println!("=== Bowtie: generalized Bell pair on every mirror pair ===\n");
    let mut psi = Mps::zero_state(&profile, TruncSpec::exact());
    zigzag::bowtie(&profile).run_mps(&mut psi);

    let mut reference = DenseState::zero_state(&profile);
    zigzag::bowtie(&profile).run_dense(&mut reference);
    let fidelity = psi.to_dense().fidelity(&reference);

    println!("bond dimensions:      {:?}", psi.bond_dims());
    let entropies = psi.bond_entropies_bits();
    println!(
        "bond entropies (bits): {:?}",
        entropies
            .iter()
            .map(|e| (e * 1000.0).round() / 1000.0)
            .collect::<Vec<_>>()
    );
    println!("  (theory: 1, 1+log2 3 = 2.585, log2 24 = 4.585, 4.585, 2.585, 1)");
    println!("max Schmidt rank across the waist bond: {} (the maximum possible for this profile: min(2·3·4, 5·4·3·2) = 24 — saturated)", psi.max_bond_dim());
    println!("fidelity vs dense ground truth: {:.15}", fidelity);
    println!(
        "storage: dense = {} amplitudes; MPS = {} parameters (a single diamond is small — the tensor-network win appears on extended waves, see the long_wave example)\n",
        reference.total_dim(),
        psi.param_count()
    );

    // ---- How much bond dimension does diamond dynamics need? -------------
    println!("=== Bond-dimension budget for deep cross-scale dynamics ===\n");
    let mut rng = Rng::new(20260723);
    let mut circuit = zigzag::brickwork_random(&profile, 3, &mut rng);
    for _ in 0..2 {
        circuit.extend(zigzag::crosswave_round(&profile));
    }
    let mut exact = DenseState::zero_state(&profile);
    circuit.run_dense(&mut exact);

    println!(
        "circuit: 3 brickwork layers + 2 cross-scale waves ({} gates)",
        circuit.len()
    );
    println!(
        "{:>8} {:>14} {:>16} {:>12}",
        "chi", "infidelity", "discarded tally", "params"
    );
    for max_bond in [2usize, 4, 6, 8, 12, 16, 24] {
        let mut trial = Mps::zero_state(&profile, TruncSpec::new(max_bond, 0.0));
        circuit.run_mps(&mut trial);
        let infid = (1.0 - trial.to_dense().fidelity(&exact)).max(0.0);
        println!(
            "{:>8} {:>14.3e} {:>16.3e} {:>12}",
            max_bond,
            infid,
            trial.discarded_weight,
            trial.param_count()
        );
    }
    println!("\nchi = 24 is exact for this geometry: the profile itself caps the reachable");
    println!("Schmidt rank — the dimension wave is a built-in entanglement budget.");
}
