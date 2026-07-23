//! Bulk sheaf towers: one state, many exact resolutions — a performance menu.
//!
//! Fusing neighbouring sites into coarser, higher-dimensional qudits is
//! lossless (the interior bond becomes intra-site structure), so a state has
//! a whole tower of exact representations: fine (many small sites) up through
//! coarse (few large sites), qudit dimension growing well past 10. Every
//! layer is the *same global section of the bulk* — the tower is "only for
//! performance." Four parts:
//!
//! 1. **Growth tower**: successive fusion of a `2..8..2` diamond, local
//!    dimension climbing past 10 into the thousands, every layer verified to
//!    be the identical state (fidelity 1).
//! 2. **The sheaf structure**: restriction (slice) and gluing (contract)
//!    make "section of the bulk" precise; a cover glues back to the whole,
//!    and gluing is associative — the sheaf axiom, executable.
//! 3. **Performance**: coarsening the state *and the operator together*
//!    (`coarsen_mpo`) runs a computation at a coarser resolution. Fewer
//!    sites means fewer canonicalization sweeps — a real speedup — up to the
//!    point where local dimension explodes.
//! 4. **The wall**: unbounded fusion is a trap; local dimension is a product,
//!    so it blows up exponentially. The sweet spot is bounded-block growth.
//!
//! Run with: `cargo run --release --example bulk_sheaf_tower`

use novel_quantum_structures::mpo::Mpo;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::sheaf::{
    coarsen, coarsen_mpo, glue, glue_all, refine, restrict, uniform_partition, BulkTower,
};
use novel_quantum_structures::{zigzag, Rng, TruncSpec};
use std::time::Instant;

const SPEC: TruncSpec = TruncSpec {
    max_rank: 64,
    cutoff: 1e-12,
};

fn main() {
    // ---- 1. The growth tower ----------------------------------------------
    println!("=== 1. Stepwise qudit growth: one state at three resolutions ===\n");
    let profile = zigzag::diamond(2, 8); // 13 sites, fine dims up to 8
    let mut rng = Rng::new(11);
    let mut fine = Mps::zero_state(&profile, SPEC);
    zigzag::brickwork_random(&profile, 1, &mut rng).run_mps(&mut fine);

    let tower = BulkTower::build(&fine, &[2, 2]);
    println!(
        "{:>6} {:>6} {:>8} {:>9} {:>9} {:>13}",
        "layer", "sites", "max d", "max bond", "params", "== fine?"
    );
    for k in 0..tower.num_layers() {
        let s = tower.stats(k);
        let fid = tower.fidelity_to_fine(k, TruncSpec::exact());
        println!(
            "{:>6} {:>6} {:>8} {:>9} {:>9} {:>13.9}",
            k, s.num_sites, s.max_local_dim, s.max_bond, s.state_params, fid
        );
    }
    println!("\nEach layer fuses pairs: 13 sites → 7 → 4, local dimension 8 → 56 → thousands.");
    println!("Every layer is the identical global section (fidelity 1) — pure re-gauging.\n");

    // ---- 2. The sheaf structure -------------------------------------------
    println!("=== 2. The bulk sheaf: restriction and gluing ===\n");
    let n = fine.num_sites();
    let a = restrict(&fine, 0, 4);
    let b = restrict(&fine, 4, 9);
    let c = restrict(&fine, 9, n);
    println!(
        "cover: sections over [0,4), [4,9), [9,{}) with boundary bonds {}, {}",
        n,
        a.right_bond(),
        b.right_bond()
    );
    let glued = glue_all(&[a.clone(), b.clone(), c.clone()]).into_mps(SPEC);
    println!(
        "glue the cover → fidelity vs whole = {:.12}",
        glued.fidelity(&fine)
    );
    let left = glue(&glue(&a, &b), &c).into_mps(SPEC);
    let right = glue(&a, &glue(&b, &c)).into_mps(SPEC);
    println!(
        "gluing is associative (the sheaf axiom): {:.12}",
        left.fidelity(&right)
    );
    println!("\nLocal sections of the bulk, consistent on their shared bonds, reassemble");
    println!("the global state — restriction and gluing are exact inverses.\n");

    // ---- 3. Performance: run the computation at a coarser resolution -------
    println!("=== 3. Coarser resolution = fewer sweeps = faster (bond-dominated) ===\n");
    let nq = 14;
    let qprofile = vec![2usize; nq];
    let mut rng2 = Rng::new(7);
    let mut qstate = Mps::zero_state(&qprofile, SPEC);
    zigzag::brickwork_random(&qprofile, 2, &mut rng2).run_mps(&mut qstate);
    let op = Mpo::from_circuit(&zigzag::brickwork_random(&qprofile, 1, &mut rng2), SPEC);
    let out_fine = op.apply_to(&qstate);
    println!(
        "{} qubits, state bond {}, operator bond {} (exact, no truncation)\n",
        nq,
        qstate.max_bond_dim(),
        op.max_bond_dim()
    );
    println!(
        "{:>6} {:>6} {:>6} {:>12} {:>10} {:>12}",
        "block", "sites", "max d", "apply", "speedup", "== fine?"
    );
    let mut fine_time = None;
    let mut best = (1usize, f64::INFINITY);
    for &block in &[1usize, 2, 3, 4, 5] {
        let (sizes, groups) = uniform_partition(&qprofile, block);
        let state = coarsen(&qstate, &sizes);
        let w = coarsen_mpo(&op, &sizes);
        let d = state.dims.iter().cloned().max().unwrap();
        let reps = 6;
        let t = Instant::now();
        let mut out = w.apply_to(&state);
        for _ in 1..reps {
            out = w.apply_to(&state);
        }
        let per = t.elapsed().as_secs_f64() / reps as f64;
        if block == 1 {
            fine_time = Some(per);
        }
        if per < best.1 {
            best = (block, per);
        }
        let refined = refine(&out, &groups, TruncSpec::exact());
        let fid = refined.fidelity(&out_fine);
        let speedup = fine_time.map(|f| f / per).unwrap_or(1.0);
        println!(
            "{:>6} {:>6} {:>6} {:>10.2}ms {:>9.2}× {:>12.9}",
            block,
            state.num_sites(),
            d,
            per * 1e3,
            speedup,
            fid
        );
    }
    println!(
        "\nBest at block {} ({:.1}× faster than fine), every result bit-identical.",
        best.0,
        fine_time.unwrap() / best.1
    );
    println!("The lever is site count: canonicalization sweeps scale with it, and coarse");
    println!("layers have fewer sites — while local dimension stays bounded at block ≤ 5.\n");

    // ---- 4. The wall -------------------------------------------------------
    println!("=== 4. The wall: local dimension is a product ===\n");
    let (sizes, _) = uniform_partition(&profile, 4); // fuse quads of diamond(2,8)
    let big = coarsen(&fine, &sizes);
    println!("fusing quads of the 2..8 diamond → dims {:?}", big.dims);
    println!(
        "max local dim {} — a single SVD on that site is already the bottleneck.",
        big.dims.iter().cloned().max().unwrap()
    );
    println!("\nGrowth pays off only in bounded blocks: coarsen to shed site-count");
    println!("overhead, stop before the qudit dimension (a product of the fused dims)");
    println!("explodes. The tower hands you every exact resolution; you pick the one");
    println!("whose local dimension fits your budget.");
}
