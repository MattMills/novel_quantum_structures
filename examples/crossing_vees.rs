//! Crossing dimension-wave strands: two peaks laid antiparallel into an X,
//! interacting at a chosen dimension level.
//!
//! The structure the user asked for: `[1,2,3,4,5,4,3,2,1]` crossing a
//! second copy running the other way, coupled where their dimensions
//! match — at level 5, 4, 3, 2, or 1 — with cycles that couple level by
//! level. Six experiments on one shared MPS:
//!
//! 1. the geometry: the X, its crossing pairs, the predicted A|B budget;
//! 2. single-level crossings: measured A|B entropy and rank vs the
//!    closed-form budget, cross-validated against dense — and the
//!    **multiplicity-vs-dimension tradeoff** (the shoulder d=4 injects
//!    more than the peak d=5; the d=1 pinch injects nothing);
//! 3. cycles: couple level by level, entanglement accumulating additively;
//! 4. **layout is the cost** — the same crossing is χ ≤ hi in interleaved
//!    layout and Π d in block layout; choosing the ordering is a ~1000×
//!    performance decision, and the gap is itself a measurement;
//! 5. at scale: a 50-site twin wave, full crossing, in microseconds;
//! 6. the valley dual: wide ends paired, pinch unique.
//!
//! Run with: `cargo run --release --example crossing_vees`

use novel_quantum_structures::crossing::{Crossing, Layout};
use novel_quantum_structures::dense::DenseState;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{gates, zigzag, Rng, TruncSpec};
use std::time::Instant;

const SPEC: TruncSpec = TruncSpec {
    max_rank: 4096,
    cutoff: 1e-14,
};

fn main() {
    let strand = zigzag::diamond(1, 5); // [1,2,3,4,5,4,3,2,1]
    let x = Crossing::twin(&strand, true);

    // ---- 1. The geometry ---------------------------------------------------
    println!("=== 1. Two peaks crossed into an X ===\n");
    println!("strand A:  {}", fmt_row(&x.a));
    println!("           {}", crossing_art(&x));
    println!(
        "strand B:  {}   (antiparallel: A[i] ~ B[n−1−i])\n",
        fmt_row(&x.b)
    );
    println!(
        "shared MPS (block layout): sites 0..{} are A, {}..{} are B; A|B bond = {}.",
        x.n(),
        x.n(),
        2 * x.n(),
        x.ab_cut()
    );
    println!("crossing pairs by level (the interaction knob):\n");
    println!(
        "{:>7} {:>7} {:>26} {:>14} {:>8}",
        "level", "#pairs", "A-sites ~ B-sites", "entropy (bits)", "rank"
    );
    for level in [5usize, 4, 3, 2, 1] {
        let pairs = x.pairs_at(&[level]);
        let desc: Vec<String> = pairs
            .iter()
            .map(|&(ga, gb, _)| format!("{}~{}", ga, gb))
            .collect();
        println!(
            "{:>7} {:>7} {:>26} {:>14.4} {:>8}",
            level,
            pairs.len(),
            desc.join(" "),
            x.bell_entropy_bits(&[level]),
            x.bell_rank(&[level])
        );
    }
    println!("\nThe peak (d=5) is UNIQUE; every lower level is PAIRED. Multiplicity is");
    println!("about to beat dimension.\n");

    // ---- 2. Single-level crossings: budget, tradeoff, dense check ----------
    println!("=== 2. Interacting at one level: the entanglement budget ===\n");
    let dims = x.dims();
    let small_strand = zigzag::diamond(1, 4); // peak 4, dim 144 — dense-checkable
    let xs = Crossing::twin(&small_strand, true);
    let sdims = xs.dims();

    println!(
        "{:>7} {:>16} {:>16} {:>12} {:>14}",
        "level", "measured S(A|B)", "predicted", "max χ", "dense fidelity"
    );
    let maxe = 4.0_f64;
    for level in [5usize, 4, 3, 2, 1] {
        let mut m = Mps::zero_state(&dims, SPEC);
        x.bell(&[level]).run_mps(&mut m);
        let measured = x.ab_entropy_bits(&mut m);
        let predicted = x.bell_entropy_bits(&[level]);
        let fid = if small_strand.contains(&level) {
            let mut ms = Mps::zero_state(&sdims, SPEC);
            xs.bell(&[level]).run_mps(&mut ms);
            let mut d = DenseState::zero_state(&sdims);
            xs.bell(&[level]).run_dense(&mut d);
            Some(ms.to_dense().fidelity(&d))
        } else {
            None
        };
        let bar = "#".repeat(((measured / maxe) * 40.0).round() as usize);
        println!(
            "{:>7} {:>16.4} {:>16.4} {:>12} {:>14}",
            level,
            measured,
            predicted,
            m.max_bond_dim(),
            fid.map(|f| format!("{:.9}", f)).unwrap_or_else(|| "—".into())
        );
        println!("        {}", bar);
    }
    println!("\nS(A|B) = Σ log2 d over the ACTIVE pairs. The maximum is the SHOULDER");
    println!("(d=4, two pairs, 4.00 bits), not the peak (d=5, one pair, 2.32 bits):");
    println!("a palindrome's peak is a bottleneck of multiplicity, not of dimension.");
    println!("Level 1 is the pinch — cshift(1,1) = I — a crossing that does nothing.\n");

    // ---- 3. Cycles: couple level by level ----------------------------------
    println!("=== 3. Cycles: interacting at a different level each round ===\n");
    println!("from |0…0⟩ on the peak-4 twin, one level activated per cycle:\n");
    println!(
        "{:>7} {:>10} {:>18} {:>18} {:>10}",
        "cycle", "level", "cumulative S(A|B)", "predicted (Σ)", "max χ"
    );
    let mut psi = Mps::zero_state(&sdims, SPEC);
    let mut active: Vec<usize> = Vec::new();
    for (cycle, level) in [4usize, 3, 2].iter().enumerate() {
        xs.bell(&[*level]).run_mps(&mut psi);
        active.push(*level);
        println!(
            "{:>7} {:>10} {:>18.4} {:>18.4} {:>10}",
            cycle + 1,
            level,
            xs.ab_entropy_bits(&mut psi),
            xs.bell_entropy_bits(&active),
            psi.max_bond_dim()
        );
    }
    println!("\nDisjoint levels commute, so the cycles ADD: 2.00 + 3.17 + 2.00 = 7.17 bits,");
    println!("the full peak-4 crossing budget reached one level at a time.\n");

    // ---- 4. Layout is the cost --------------------------------------------
    println!("=== 4. The same crossing, two layouts: a ~1000× cost gap ===\n");
    println!("A Bell crossing is a product of LOCAL pairs. Block layout forces all of");
    println!("it through one contiguous cut; interleaved layout keeps each pair local.\n");
    println!(
        "{:>13} {:>10} {:>10} {:>12} {:>12}",
        "layout", "max χ", "params", "time", "A|B bits"
    );
    for (name, layout) in [("block", Layout::Block), ("interleaved", Layout::Interleaved)] {
        let xc = Crossing::twin(&small_strand, true).with_layout(layout);
        let t = Instant::now();
        let mut m = Mps::zero_state(&xc.dims(), SPEC);
        xc.bell(&[1, 2, 3, 4]).run_mps(&mut m);
        let elapsed = t.elapsed();
        // A|B entropy: a single bond in block; the closed-form budget in
        // interleaved (where A|B is a comb, not one cut).
        let ab = match layout {
            Layout::Block => xc.ab_entropy_bits(&mut m),
            Layout::Interleaved => xc.bell_entropy_bits(&[1, 2, 3, 4]),
        };
        println!(
            "{:>13} {:>10} {:>10} {:>12.1?} {:>12.3}",
            name,
            m.max_bond_dim(),
            m.param_count(),
            elapsed,
            ab
        );
    }
    println!("\nSame state, same 7.17 bits of inter-strand entanglement — but block layout");
    println!("pays Π d = 144 in bond dimension and interleaved pays hi = 4. The A|B");
    println!("entanglement is real; whether it is EXPENSIVE is a choice of site order.");
    println!("Interleave to compute; use the block bond to read the entanglement off.\n");

    // ---- 5. At scale -------------------------------------------------------
    println!("=== 5. At scale: a 50-site twin wave, fully crossed ===\n");
    let big = zigzag::wave(1, 5, 3); // three peaks per strand
    let strand_dim: f64 = big.iter().map(|&d| d as f64).product();
    let xb = Crossing::twin(&big, true).with_layout(Layout::Interleaved);
    let t = Instant::now();
    let mut m = Mps::zero_state(&xb.dims(), SPEC);
    xb.bell(&[1, 2, 3, 4, 5]).run_mps(&mut m);
    println!(
        "strands = wave(1,5,3) [{} sites/strand → {}-site MPS; dim {:.2e} each, dense ≈ {:.1e}]",
        big.len(),
        2 * big.len(),
        strand_dim,
        strand_dim * strand_dim
    );
    println!(
        "full crossing (all 5 levels, {} pairs): max χ = {}, {} params, {:.0?}",
        xb.pairs().len(),
        m.max_bond_dim(),
        m.param_count(),
        t.elapsed()
    );
    println!(
        "predicted A|B entanglement: {:.2} bits — held in {:.1} KB.",
        xb.bell_entropy_bits(&[1, 2, 3, 4, 5]),
        m.param_count() as f64 * 16.0 / 1e3
    );
    // Scrambling contrast on the small twin.
    let mut rng = Rng::new(2026);
    let mut scrambled = Mps::zero_state(&sdims, TruncSpec::new(64, 1e-12));
    xs.spread(&[2, 3, 4]).run_mps(&mut scrambled);
    for (ga, gb, d) in xs.pairs_at(&[2, 3, 4]) {
        scrambled.apply2(ga, gb, &gates::random(d * d, &mut rng));
    }
    println!(
        "\ncontrast (peak-4 twin): Haar-random cross-coupling saturates any budget —"
    );
    println!(
        "max χ = {} (capped 64), truncation {:.2e}: structure, not the geometry, is",
        scrambled.max_bond_dim(),
        scrambled.discarded_weight
    );
    println!("what keeps a crossing compressible.\n");

    // ---- 6. The valley dual -----------------------------------------------
    println!("=== 6. The valley dual: multiplicity reflected ===\n");
    let vee = zigzag::valley(1, 5); // [5,4,3,2,1,2,3,4,5]
    let xv = Crossing::twin(&vee, true);
    println!(
        "strand:  {}   (a valley — wide ends, narrow pinch)",
        fmt_row(&vee)
    );
    println!("{:>7} {:>7} {:>14} {:>8}", "level", "#pairs", "entropy", "rank");
    for level in [5usize, 4, 3, 2, 1] {
        println!(
            "{:>7} {:>7} {:>14.4} {:>8}",
            level,
            xv.pairs_at(&[level]).len(),
            xv.bell_entropy_bits(&[level]),
            xv.bell_rank(&[level])
        );
    }
    println!("\nNow the WIDE ENDS (d=5) are paired and the PINCH (d=1) is unique: the");
    println!("valley crosses most richly at its rims (2·log2 5 = 4.64 bits), the");
    println!("diamond at its shoulders. The same five dimensions, reflected");
    println!("multiplicities — where a crossing can entangle two waves is a property");
    println!("of the SHAPE, not just the sizes.");
}

fn fmt_row(profile: &[usize]) -> String {
    profile
        .iter()
        .map(|d| format!("{:>2}", d))
        .collect::<Vec<_>>()
        .join(" ")
}

/// A crude ASCII of the antiparallel crossing arms.
fn crossing_art(x: &Crossing) -> String {
    let n = x.n();
    (0..n)
        .map(|i| {
            let c = if i < n / 2 {
                '\\'
            } else if i > n / 2 {
                '/'
            } else {
                'X'
            };
            format!("{:>2}", c)
        })
        .collect::<Vec<_>>()
        .join(" ")
}
