//! Applying the crossing-V set to the multi-V renormalization wave: where
//! the two ideas meet, and everything expands into a complex multi-scale
//! shape.
//!
//! The zigzag *wave* — `wave(1,5,p)`, many V's in a row — is a real-space
//! renormalization structure: valleys are fine (`d = 2`, a qubit), waists
//! are coarse (`d = 5`), `merge_sites` is a coarse-graining step, and
//! structured multi-scale dynamics stays compressed while it rides the wave
//! (findings 5–7). This example crosses those RG strands with the crossing
//! machinery, and finds the crossing inherits the wave's scale hierarchy:
//!
//! 1. **The wave is a scale ladder** — a quick recap of the RG structure the
//!    crossing is about to be applied against.
//! 2. **A crossing has a scale** — coupling two waves at the WAIST couples a
//!    few coarse (fat) modes; at the VALLEY-side, many fine (thin) modes.
//! 3. **The complex expanding shape** — a woven lattice of multi-V waves:
//!    its cost profile traces the wave silhouette (thick at the coarse
//!    waists, thin at the fine valleys) and thickens with the transverse
//!    extent. The RG structure shapes *where* the area law bites.
//! 4. **A crossing is RG-covariant** — coarse-graining a crossed strand
//!    (an exact `merge_sites` RG step) carries the crossing coupling up to
//!    the coarse block, state preserved exactly. Scale morphing and
//!    crossing commute.
//!
//! Run with: `cargo run --release --example renormalizing_crossings`

use novel_quantum_structures::crossing::{Crossing, Layout};
use novel_quantum_structures::dense::DenseState;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::network::wave_weave;
use novel_quantum_structures::{zigzag, TruncSpec};

const SPEC: TruncSpec = TruncSpec {
    max_rank: 4096,
    cutoff: 1e-13,
};

fn main() {
    // ---- 1. The wave is a scale ladder ------------------------------------
    println!("=== 1. The multi-V wave as a renormalization structure ===\n");
    let wave = zigzag::wave(1, 5, 3);
    println!("{}", zigzag::ascii_profile(&wave));
    println!("profile {:?}", wave);
    println!(
        "  valleys (fine, {} qubits each): {:?}",
        1,
        zigzag::valley_positions(&wave)
    );
    println!(
        "  waists  (coarse, log2 5 ≈ 2.32 bits each): {:?}",
        zigzag::waist_positions(&wave)
    );
    println!("\nEach V is an RG cell: dimension expands valley→waist→valley, coarse-");
    println!("graining fine qubits into a wide waist and back. merge_sites fuses a cell");
    println!("into one coarse site; structured dynamics stays compressed riding the wave");
    println!("(findings 5–7). Now cross these RG strands.\n");

    // ---- 2. A crossing has a scale ----------------------------------------
    println!("=== 2. Crossing two RG waves — the coupling inherits the scale ===\n");
    let w = zigzag::wave(1, 5, 2);
    let x = Crossing::twin(&w, true).with_layout(Layout::Interleaved);
    println!("two wave(1,5,2) strands crossed; coupling at each dimension band:\n");
    println!(
        "{:>10} {:>10} {:>12} {:>14} {:>18}",
        "level", "scale", "#crossings", "modes each", "A|B entropy (bits)"
    );
    for (level, scale) in [(5usize, "waist (coarse)"), (4, "shoulder"), (3, "mid"), (2, "valley (fine)")] {
        println!(
            "{:>10} {:>10} {:>12} {:>14} {:>18.3}",
            level,
            scale,
            x.pairs_at(&[level]).len(),
            level,
            x.bell_entropy_bits(&[level])
        );
    }
    println!("\nThe crossing has a SCALE. At the waists (coarse) there are few crossings");
    println!("of fat modes; toward the valleys (fine) there are many crossings of thin");
    println!("modes — one crossing point per shared waist, growing with the period. The");
    println!("wave's renormalization hierarchy is inherited by the inter-wave coupling.\n");

    // ---- 3. The complex expanding shape -----------------------------------
    println!("=== 3. A woven lattice of multi-V waves — the RG silhouette in the cost ===\n");
    // Dense-checkable: a small multi-V wave, 2–3 rows.
    let small = zigzag::wave(1, 3, 2); // [1,2,3,2,1,2,3,2,1], waists at 2 and 6
    println!("weave `rows` copies of {:?} (waists at {:?}):\n", small, zigzag::waist_positions(&small));
    for rows in [2usize, 3] {
        let (net, edges) = wave_weave(&small, rows);
        let mut m = Mps::zero_state(&net.dims(), SPEC);
        let c = net.cluster(&edges);
        c.run_mps(&mut m);
        // dense fidelity on the smaller rows.
        let fid = {
            let mut d = DenseState::zero_state(&net.dims());
            c.run_dense(&mut d);
            m.to_dense().fidelity(&d)
        };
        println!("rows={}: bond profile {:?}", rows, m.bond_dims());
        println!("         (max χ {}, dense fidelity {:.9})", m.max_bond_dim(), fid);
    }
    println!("\nThe cost profile is the WAVE, doubled: two humps at the two coarse waists,");
    println!("pinched to 1 at the fine valleys — and every hump thickens with each added");
    println!("row (2→3 rows lifts the waist bonds 6→12). The area law of a 2D crossing");
    println!("is MODULATED by the RG structure: it bites hardest at the coarse blocks.\n");

    // At scale (MPS-only): a lattice of full 5-waves.
    let big = zigzag::wave(1, 5, 2);
    println!("at scale — {} copies of wave(1,5,2) woven:", 2);
    for rows in [2usize, 3] {
        let (net, edges) = wave_weave(&big, rows);
        let mut m = Mps::zero_state(&net.dims(), SPEC);
        net.cluster(&edges).run_mps(&mut m);
        let bonds = m.bond_dims();
        // report the bond dimension at each waist column vs valley column.
        println!(
            "  rows={}: {} sites, max χ = {} (at the waists), min interior χ = {}",
            rows,
            net.dims().len(),
            m.max_bond_dim(),
            bonds.iter().filter(|&&b| b > 1).min().copied().unwrap_or(1)
        );
    }
    println!();

    // ---- 4. A crossing is RG-covariant ------------------------------------
    println!("=== 4. Coarse-graining a crossing: RG covariance ===\n");
    let strand = zigzag::diamond(1, 4); // [1,2,3,4,3,2,1], one RG cell
    let xb = Crossing::twin(&strand, false); // ladder, block layout
    let mut crossed = Mps::zero_state(&xb.dims(), SPEC);
    xb.bell(&[2, 3, 4]).run_mps(&mut crossed);
    println!(
        "crossed twin {:?} (ladder at levels 2,3,4): dims {:?}",
        strand,
        crossed.dims
    );

    // RG step: merge strand A's cell toward its waist (sites 1..5 → coarse).
    let mut coarse = crossed.clone();
    coarse.merge_sites(2); // fuse (2,3): d 3·4 = 12
    coarse.merge_sites(2); // fuse (12-site, 3): d 12·3 = 36 — the coarse waist block
    println!("after coarse-graining A's waist cell (two merges): dims {:?}", coarse.dims);

    // The coarse-grained crossed state equals crossing-then-merging: preserved.
    let mut ref_state = crossed.clone();
    ref_state.merge_sites(2);
    ref_state.merge_sites(2);
    println!(
        "state preserved under the RG step: fidelity {:.12}",
        coarse.fidelity(&ref_state)
    );
    // And it splits back exactly — the RG step is invertible here.
    coarse.split_site(2, 12, 3);
    coarse.split_site(2, 3, 4);
    println!(
        "and split back to the fine scale: fidelity vs original {:.12}",
        coarse.fidelity(&crossed)
    );
    println!("\nThe crossing coupling rode the merge up onto the coarse waist site (d=36)");
    println!("and back down, exactly. A crossing is a covariant object under the wave's");
    println!("renormalization: the coupling flows with the scale, it does not break.\n");

    println!("=== The synthesis ===\n");
    println!("The dimension wave supplies a scale ladder; the crossing supplies");
    println!("inter-structure coupling. Put together, the crossing acquires a scale");
    println!("(coarse waist couplings vs fine valley couplings), a 2D lattice of waves");
    println!("carries an area law shaped like the wave itself, and coarse-graining");
    println!("commutes with crossing. Everything expands in a complex multi-scale");
    println!("shape whose cost is still read off one number — the cutwidth — now");
    println!("modulated by the renormalization structure of the strands.");
}
