//! A four-diamond dimension wave: 25 sites, dense dimension ≈ 8.6·10¹²
//! (≈ 138 TB of amplitudes — no dense simulator will ever hold it).
//!
//! Two contrasting experiments, both impossible dense:
//!
//! * **Structured cross-scale dynamics** — Fourier spreading, a controlled-
//!   shift staircase riding the wave, per-diamond mirror ("bowtie")
//!   couplings, valley↔valley and edge↔edge Bell couplings spanning up to
//!   24 sites. This runs *numerically exactly* (zero truncation) at modest
//!   bond dimension: the zigzag's classical computability, quantified.
//! * **Generic scrambling** — Haar-random brickwork. Entanglement grows to
//!   whatever the geometry allows and truncation becomes severe: compression
//!   is a property of *structured* multi-scale dynamics, not of the lattice
//!   by itself (the same lesson MERA teaches).
//!
//! Run with: `cargo run --release --example long_wave`

use novel_quantum_structures::circuit::Circuit;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{gates, zigzag, Rng, TruncSpec};
use std::time::Instant;

/// The structured cross-scale circuit on a wave profile: a heterogeneous
/// cousin of cluster-state preparation. Fourier layers make every site
/// uniform; controlled-phase couplings (which *do* entangle uniform states)
/// wire up mirror pairs, valleys, and the outermost edges across scales;
/// a final controlled-shift staircase propagates the phase structure along
/// the wave.
fn structured_crosswave(profile: &[usize]) -> Circuit {
    let n = profile.len();
    let valleys = zigzag::valley_positions(profile);
    let mut c = Circuit::new(profile.to_vec());
    // Spread every site over all its levels.
    for (i, &d) in profile.iter().enumerate() {
        c.one(i, gates::fourier(d));
    }
    // Within-diamond mirror couplings (the bowtie flanks of each diamond).
    for w in valleys.windows(2) {
        let (lo, hi) = (w[0], w[1]);
        let span = hi - lo;
        for k in 1..span / 2 {
            let (i, j) = (lo + k, hi - k);
            c.two(i, j, gates::cphase(profile[i], profile[j]));
        }
    }
    // Valley↔valley couplings: each tunnels through a whole diamond.
    for w in valleys.windows(2) {
        c.two(w[0], w[1], gates::cphase(2, 2));
    }
    // Edge↔edge coupling across the entire wave (span n-1).
    c.two(0, n - 1, gates::cphase(2, 2));
    // Controlled-shift staircase up and down the dimension wave.
    for i in 0..n - 1 {
        c.two(i, i + 1, gates::cshift(profile[i], profile[i + 1]));
    }
    c
}

fn main() {
    let profile = zigzag::wave(2, 5, 4);
    let n = profile.len();
    println!("=== Four-period dimension wave ===\n");
    println!("{}", zigzag::ascii_profile(&profile));
    let dense_dim = zigzag::dense_dimension(&profile);
    println!(
        "sites: {}   dense dimension: {:.3e} amplitudes (≈ {:.0} TB — dense simulation is impossible)",
        n,
        dense_dim,
        dense_dim * 16.0 / 1e12
    );
    println!(
        "valleys (qubit edges): {:?}   waists: {:?}\n",
        zigzag::valley_positions(&profile),
        zigzag::waist_positions(&profile)
    );

    // ---- Part A: structured cross-scale dynamics, exactly ----------------
    println!("=== A. Structured cross-scale circuit — exact at modest chi ===\n");
    let circuit = structured_crosswave(&profile);
    let t0 = Instant::now();
    let mut psi = Mps::zero_state(&profile, TruncSpec::new(256, 0.0));
    circuit.run_mps(&mut psi);
    let elapsed = t0.elapsed();
    println!(
        "gates: {} (incl. mirror couplings spanning up to {} sites)",
        circuit.len(),
        n - 1
    );
    println!("max bond dimension reached: {}", psi.max_bond_dim());
    println!(
        "truncation:                 {:.1e} (numerically exact)",
        psi.discarded_weight
    );
    println!(
        "storage: {} parameters ≈ {:.2} MB   (dense: {:.0} TB)   wall time {:.1?}",
        psi.param_count(),
        psi.param_count() as f64 * 16.0 / 1e6,
        dense_dim * 16.0 / 1e12,
        elapsed
    );

    let ent = psi.bond_entropies_bits();
    println!("\nentanglement profile (bits per bond) — the wave shapes where correlations live:");
    let maxe = ent.iter().cloned().fold(0.0, f64::max).max(1e-9);
    for (b, e) in ent.iter().enumerate() {
        let bar = "#".repeat(((e / maxe) * 40.0).round() as usize);
        println!(
            "  bond {:2} ({}|{})  {:6.3}  {}",
            b,
            profile[b],
            profile[b + 1],
            e,
            bar
        );
    }

    // Convergence in chi: a capped rerun agrees with the uncapped run.
    let chi_cap = psi.max_bond_dim();
    let mut capped = Mps::zero_state(&profile, TruncSpec::new(chi_cap, 0.0));
    circuit.run_mps(&mut capped);
    println!(
        "\nrerun with chi capped at {}: fidelity vs uncapped = {:.12}",
        chi_cap,
        capped.fidelity(&psi)
    );

    // ---- Part B: generic scrambling saturates any budget ------------------
    println!("\n=== B. Haar-random brickwork — scrambling defeats fixed budgets ===\n");
    let mut rng = Rng::new(11);
    let scramble = zigzag::brickwork_random(&profile, 2, &mut rng);
    let t1 = Instant::now();
    let mut chaos = Mps::zero_state(&profile, TruncSpec::new(32, 0.0));
    // Start from the structured state to make the contrast direct.
    circuit.run_mps(&mut chaos);
    scramble.run_mps(&mut chaos);
    println!("2 Haar brickwork layers on top of the structured state, chi capped at 32:");
    println!(
        "max chi {}   truncation tally {:.3} (severe)   wall time {:.1?}",
        chaos.max_bond_dim(),
        chaos.discarded_weight,
        t1.elapsed()
    );
    println!("\nGeneric dynamics fills the geometry's entanglement budget; the zigzag's");
    println!("compression is a property of structured multi-scale circuits — the same");
    println!("boundary MERA draws between renormalizable and volume-law dynamics.");
}
