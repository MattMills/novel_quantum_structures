//! Networks of crossing strands: from bundles to weaves, and the cutwidth
//! that separates the cheap from the expensive.
//!
//! One X of two strands was finding 24. This is the general geometry — many
//! strands, coupled in a graph — organized by one law:
//!
//!   **max bond dimension = d^(edges crossing the worst contiguous cut),
//!   minimized over site orderings by the coupling graph's CUTWIDTH.**
//!
//! Five configurations, each a data point on that law:
//!
//! 1. **bundle** — K crossing V's coupled in ONE direction (a GHZ rung per
//!    site): cutwidth 1, χ = hi for any K;
//! 2. **overlay** — two strands crossed BOTH ways (X + ladder): cutwidth 2,
//!    χ = (hi−1)² (the shoulder squared);
//! 3. **weave** — strands in TWO transverse directions (a 2D lattice):
//!    cutwidth min(rows,cols), χ = d^rows — an area law, exponential in the
//!    transverse extent;
//! 4. **wave weave** — the same lattice built from dimension waves: the area
//!    law with a wave-shaped base;
//! 5. **multi-period crossing** — two VVVVV waves crossed at every shared
//!    peak: multiple crossing points, still cutwidth-1 cheap.
//!
//! Run with: `cargo run --release --example crossing_networks`

use novel_quantum_structures::crossing::{Crossing, Layout};
use novel_quantum_structures::dense::DenseState;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::network::{bundle, overlay, wave_weave, weave, Network};
use novel_quantum_structures::circuit::Circuit;
use novel_quantum_structures::{zigzag, TruncSpec};

const SPEC: TruncSpec = TruncSpec {
    max_rank: 4096,
    cutoff: 1e-13,
};

fn built(net: &Network, c: &Circuit) -> Mps {
    let mut m = Mps::zero_state(&net.dims(), SPEC);
    c.run_mps(&mut m);
    m
}

fn dense_fid(net: &Network, c: &Circuit) -> f64 {
    let m = built(net, c);
    let mut d = DenseState::zero_state(&net.dims());
    c.run_dense(&mut d);
    m.to_dense().fidelity(&d)
}

fn main() {
    println!("=== The cutwidth law ===\n");
    println!("A crossing network is a coupling graph on dimension-wave strands, laid");
    println!("along one MPS. Its max bond dimension is d^(edges crossing the worst cut),");
    println!("minimized over orderings by the graph's CUTWIDTH. Everything below is a");
    println!("reading of that one number.\n");

    // ---- 1. Bundle: many V's, one direction -------------------------------
    println!("=== 1. Bundle — K crossing V's, ONE direction (cutwidth 1) ===\n");
    println!("K copies of [1,2,3,4,5,4,3,2,1], GHZ-coupled across the bundle at each");
    println!("site (Σ|k…k⟩ rungs). Column-major, every rung is local:\n");
    println!("{:>5} {:>10} {:>12}", "K", "max χ", "params");
    let profile = zigzag::diamond(1, 5);
    for k in [2usize, 3, 5, 8] {
        let (net, groups) = bundle(&profile, k);
        let m = built(&net, &net.ghz(&groups));
        println!("{:>5} {:>10} {:>12}", k, m.max_bond_dim(), m.param_count());
    }
    // Correctness: a small bundle against dense.
    let (sn, sg) = bundle(&zigzag::diamond(1, 3), 3);
    println!(
        "\n(3 strands of [1,2,3,2,1] vs dense: fidelity {:.9})",
        dense_fid(&sn, &sn.ghz(&sg))
    );
    println!("max χ = hi = 5 no matter how many strands: one direction of crossing is");
    println!("cheap however wide the bundle. GHZ across K strands is a single mode.\n");

    // ---- 2. Overlay: two directions on one pair ---------------------------
    println!("=== 2. Overlay — one pair crossed BOTH ways (cutwidth 2) ===\n");
    println!("An X (A[i]~B[n−1−i]) AND a ladder (A[i]~B[i]) at once: a union of two");
    println!("matchings = disjoint 4-cycles. Ordering each cycle contiguously:\n");
    println!("{:>10} {:>10} {:>16} {:>16}", "strand", "max χ", "(hi−1)²", "dense fidelity");
    for hi in [3usize, 4, 5] {
        let p = zigzag::diamond(1, hi);
        let (net, edges) = overlay(&p);
        let c = net.cluster(&edges);
        let m = built(&net, &c);
        println!(
            "{:>10} {:>10} {:>16} {:>16.9}",
            format!("(1,{})", hi),
            m.max_bond_dim(),
            (hi - 1) * (hi - 1),
            dense_fid(&net, &c)
        );
    }
    println!("\nmax χ = (hi−1)²: overlaying a second direction adds ONE to the cutwidth.");
    println!("It is the SHOULDER squared, not the peak — the unique peak carries only");
    println!("the ladder edge (cutwidth 1); the paired shoulders reach cutwidth 2.\n");

    // ---- 3. Weave: two transverse directions ------------------------------
    println!("=== 3. Weave — strands in TWO directions, a 2D lattice (area law) ===\n");
    println!("rows × cols qudits (d=2), coupled on every horizontal AND vertical edge");
    println!("(a 2D cluster state). Column-major: the cut crosses `rows` bonds.\n");
    println!("{:>7} {:>7} {:>10} {:>14} {:>14}", "rows", "cols", "max χ", "d^rows", "dense fid");
    for &(rows, cols) in &[(2usize, 6usize), (3, 6), (4, 6), (5, 6), (3, 3), (3, 9)] {
        let (net, edges) = weave(2, rows, cols);
        let c = net.cluster(&edges);
        let m = built(&net, &c);
        let fid = if rows * cols <= 9 {
            format!("{:.6}", dense_fid(&net, &c))
        } else {
            "—".into()
        };
        println!(
            "{:>7} {:>7} {:>10} {:>14} {:>14}",
            rows,
            cols,
            m.max_bond_dim(),
            2usize.pow(rows as u32),
            fid
        );
    }
    println!("\nmax χ = 2^rows — EXPONENTIAL in the transverse extent, and flat in the");
    println!("length (rows=3 is χ=8 whether cols=3, 6, or 9). Crossing in two directions");
    println!("is the area law: the boundary where a 1D tensor network stops being");
    println!("efficient — the same line that separates MPS from PEPS.\n");

    // ---- 4. Wave weave: the area law with a wave-shaped base ---------------
    println!("=== 4. Wave weave — a lattice of dimension waves ===\n");
    println!("`rows` copies of [1,2,3,2,1] woven: vertical bonds at profile[c], horizontal");
    println!("bonds at min(profile[c],profile[c+1]). The area-law base MODULATES:\n");
    for rows in [2usize, 3] {
        let (net, edges) = wave_weave(&zigzag::diamond(1, 3), rows);
        let c = net.cluster(&edges);
        let m = built(&net, &c);
        println!(
            "rows={}: bond profile {:?}  (dense fidelity {:.9})",
            rows,
            m.bond_dims(),
            dense_fid(&net, &c)
        );
    }
    println!("\nThe bond profile is wave-shaped — thin at the dimension-1 rims, thick at");
    println!("the shared peak — and the whole profile thickens with each added row. The");
    println!("area law inherits the wave's silhouette.\n");

    // ---- 5. Multi-period crossing: many crossing points -------------------
    println!("=== 5. Multi-period waves crossing at every shared peak ===\n");
    println!("Two VVVVV waves wave(1,5,p) crossed antiparallel meet at p peak-regions —");
    println!("p crossing points — yet a single X is still cutwidth 1 (interleaved):\n");
    println!(
        "{:>8} {:>7} {:>16} {:>18} {:>10}",
        "periods", "sites", "level-5 crossings", "A|B budget (bits)", "max χ"
    );
    for periods in [1usize, 2, 3, 4] {
        let w = zigzag::wave(1, 5, periods);
        let xi = Crossing::twin(&w, true).with_layout(Layout::Interleaved);
        let mut m = Mps::zero_state(&xi.dims(), SPEC);
        xi.bell(&[1, 2, 3, 4, 5]).run_mps(&mut m);
        println!(
            "{:>8} {:>7} {:>16} {:>18.2} {:>10}",
            periods,
            xi.dims().len(),
            xi.pairs_at(&[5]).len(),
            xi.bell_entropy_bits(&[1, 2, 3, 4, 5]),
            m.max_bond_dim()
        );
    }
    println!("\nThe number of peak-to-peak crossings grows with the period, so the total");
    println!("inter-strand entanglement grows too — but each crossing stays local, so");
    println!("the bond dimension is pinned at hi = 5. Many crossing points, one direction:");
    println!("still cheap. The cost is the cutwidth, and the cutwidth is the geometry.\n");

    // ---- Summary ----------------------------------------------------------
    println!("=== The cutwidth ladder ===\n");
    println!("{:>22} {:>10} {:>18}", "geometry", "cutwidth", "max χ");
    println!("{:>22} {:>10} {:>18}", "single X (2 strands)", 1, "hi");
    println!("{:>22} {:>10} {:>18}", "bundle (K strands)", 1, "hi");
    println!("{:>22} {:>10} {:>18}", "multi-period crossing", 1, "hi");
    println!("{:>22} {:>10} {:>18}", "overlay (X + ladder)", 2, "(hi−1)²");
    println!("{:>22} {:>10} {:>18}", "weave (2D lattice)", "min(r,c)", "d^min(r,c)");
    println!("\nOne direction of crossing — however many strands, however many crossing");
    println!("points — is a constant cutwidth and stays classically cheap. A SECOND");
    println!("transverse direction is an area law. That threshold, not the number of");
    println!("strands, is what governs whether a crossing network is simulable.");
}
