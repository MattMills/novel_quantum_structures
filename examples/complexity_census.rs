//! A complexity census: the **work** (build time, stored size) and
//! **output** (bond dimension = width, entanglement, semantics) of every
//! geometric object in the crate, measured — the empirical companion to
//! COMPLEXITY.md.
//!
//! The organizing question: an MPS over `n` sites of dimension `d` with
//! bond dimension `χ` costs `O(n · χ³ · d³)` to evolve and `O(n · χ² · d)`
//! to store, against the `O(dⁿ)` of a dense vector. So **χ is the entire
//! complexity story**, and every object here is classified by how its χ
//! grows:
//!
//!   A. constant width  — χ = O(1),   linear time & space   (cheap)
//!   B. geometry-bounded — χ ≤ budget, polynomial            (structured)
//!   C. exponential      — χ = d^extent, area/volume law     (hard)
//!
//! plus two cross-cutting phenomena: **predict-vs-pay** (the width atlas
//! computes χ in `O(S)` arithmetic without building the operator) and
//! **layout-dependence** (the same state is χ = Πd or χ = hi depending on
//! site order).
//!
//! Run with: `cargo run --release --example complexity_census`

use novel_quantum_structures::cascade::{unitarity_defect, Transducer};
use novel_quantum_structures::crossing::{Crossing, Layout};
use novel_quantum_structures::flow::{FourierFlow, FrameFlow};
use novel_quantum_structures::mpo::Mpo;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::network::{bundle, weave};
use novel_quantum_structures::stabilize::iterate;
use novel_quantum_structures::width::mult_width_profile;
use novel_quantum_structures::{radix, zigzag, C64, Rng, TruncSpec};
use std::time::{Duration, Instant};

const SPEC: TruncSpec = TruncSpec {
    max_rank: 4096,
    cutoff: 1e-13,
};

fn time<T>(f: impl FnOnce() -> T) -> (T, Duration) {
    let t = Instant::now();
    let r = f();
    (r, t.elapsed())
}

/// Kilobytes of a param count (16 bytes / complex).
fn kb(params: usize) -> f64 {
    params as f64 * 16.0 / 1e3
}

fn basis(profile: &[usize], mut x: u128) -> Mps {
    let mut d = vec![0usize; profile.len()];
    for (i, &dim) in profile.iter().enumerate().rev() {
        d[i] = (x % dim as u128) as usize;
        x /= dim as u128;
    }
    Mps::basis_state(profile, &d, SPEC)
}

fn main() {
    println!("COMPLEXITY CENSUS — work (time, size) and output (width χ) of each object");
    println!("cost model: evolve O(n·χ³·d³), store O(n·χ²·d), vs dense O(dⁿ)\n");

    // ======================================================================
    // CLASS A — constant width χ = O(1): linear time and space
    // ======================================================================
    println!("========== CLASS A: constant width (χ = O(1), linear cost) ==========\n");

    // The carry adder MPO: χ = 2 for ANY ring size.
    println!("-- carry adder MPO  |x⟩→|x+1⟩  (χ = 2, exact at any N) --");
    println!(
        "{:>14} {:>6} {:>6} {:>10} {:>12} {:>12}",
        "profile", "n", "χ", "params", "size", "build"
    );
    for p in [1usize, 2, 3, 4] {
        let wave = zigzag::wave(2, 5, p);
        let n = wave.len();
        let dense = zigzag::dense_dimension(&wave);
        let (adder, dt) = time(|| radix::adder_mpo_exact(&wave, 1, SPEC));
        println!(
            "{:>14} {:>6} {:>6} {:>10} {:>9.1} KB {:>12.0?}   (dense {:.1e})",
            format!("wave(2,5,{})", p),
            n,
            adder.max_bond_dim(),
            adder.param_count(),
            kb(adder.param_count()),
            dt,
            dense
        );
    }
    println!("→ χ flat at 2, params & time LINEAR in n; dense space is exponential.\n");

    // Bundle network: χ = hi for any number of strands K.
    println!("-- bundle of K crossing strands  (χ = hi, GHZ rungs) --");
    println!("{:>6} {:>6} {:>6} {:>10} {:>12}", "K", "sites", "χ", "params", "build");
    let strand = zigzag::diamond(2, 5);
    for k in [2usize, 4, 8, 16] {
        let (net, groups) = bundle(&strand, k);
        let ((m, _), dt) = time(|| {
            let mut m = Mps::zero_state(&net.dims(), SPEC);
            net.ghz(&groups).run_mps(&mut m);
            (m, ())
        });
        println!(
            "{:>6} {:>6} {:>6} {:>10} {:>12.0?}",
            k,
            net.dims().len(),
            m.max_bond_dim(),
            m.param_count(),
            dt
        );
    }
    println!("→ χ flat at hi = 5 however many strands; cost linear in K.\n");

    // ======================================================================
    // CLASS B — geometry-bounded width: polynomial, χ ≤ budget
    // ======================================================================
    println!("========== CLASS B: geometry-bounded width (χ ≤ budget) ==========\n");

    // Bowtie state: χ = the waist entanglement budget.
    println!("-- bowtie state on diamonds  (χ = waist budget = ∏ half-dims) --");
    println!(
        "{:>14} {:>6} {:>8} {:>10} {:>12} {:>10}",
        "profile", "n", "χ=budget", "waist bits", "params", "build"
    );
    for hi in [3usize, 4, 5, 6] {
        let profile = zigzag::diamond(2, hi);
        let ((mut m, _), dt) = time(|| {
            let mut m = Mps::zero_state(&profile, TruncSpec::exact());
            zigzag::bowtie(&profile).run_mps(&mut m);
            (m, ())
        });
        let ent = m
            .bond_entropies_bits()
            .iter()
            .cloned()
            .fold(0.0, f64::max);
        println!(
            "{:>14} {:>6} {:>8} {:>10.3} {:>12} {:>10.0?}",
            format!("diamond(2,{})", hi),
            profile.len(),
            m.max_bond_dim(),
            ent,
            m.param_count(),
            dt
        );
    }
    println!("→ χ set by geometry (the waist budget), independent of chain length.\n");

    // Modular multiplier: χ = the number-theoretic true width.
    println!("-- modular multiplier ×k on Z_2880  (χ = resonant true width ≤ k) --");
    println!(
        "{:>6} {:>10} {:>10} {:>12} {:>10}",
        "k", "atlas χ", "measured", "defect", "build"
    );
    let dia = zigzag::diamond(2, 5);
    for k in [7usize, 11, 49] {
        // Direct construction is capped at modest k; larger come from
        // composition (below).
        let atlas = mult_width_profile(&dia, k as u128)
            .iter()
            .map(|w| w.value() as usize)
            .max()
            .unwrap();
        let (m, dt) = time(|| Transducer::mult(&dia, k).to_mpo(SPEC));
        assert_eq!(atlas, m.max_bond_dim());
        println!(
            "{:>6} {:>10} {:>10} {:>12.1e} {:>10.0?}",
            k,
            atlas,
            m.max_bond_dim(),
            unitarity_defect(&m, SPEC),
            dt
        );
    }
    print!("  atlas-only (deep powers exceed the direct cap; built by composition):\n  squaring orbit 7 → 49 → 2401 → 1921 has widths ");
    for k in [7u128, 49, 2401, 1921] {
        print!(
            "{} ",
            mult_width_profile(&dia, k).iter().map(|w| w.value()).max().unwrap()
        );
    }
    println!("\n→ χ is the cut-rank |{{⌊kb/S⌋ mod L}}|, non-monotone in k (×1921 < ×7).\n");

    // QFT MPO: reversed core cheap, digit reversal expensive.
    println!("-- QFT operator: the tail-radix phase web  (χ = past↔future width) --");
    for hi in [3usize, 4] {
        let p = zigzag::diamond(2, hi);
        let (vr, dtr) = time(|| Mpo::from_circuit(&radix::mixed_radix_qft_reversed(&p, 0.0), SPEC));
        let (vs, dts) = time(|| Mpo::from_circuit(&radix::mixed_radix_qft(&p, 0.0), SPEC));
        println!(
            "  diamond(2,{}): reversed core χ={:<3} ({:.0?}),  +digit-reversal χ={:<3} ({:.0?})",
            hi,
            vr.max_bond_dim(),
            dtr,
            vs.max_bond_dim(),
            dts
        );
    }
    println!("→ the Fourier core is O(1) wide; the bowtie digit-reversal carries the χ.\n");

    // ======================================================================
    // CLASS C — exponential width: area / volume law
    // ======================================================================
    println!("========== CLASS C: exponential width (area / volume law) ==========\n");

    // Weave: two-direction crossing, χ = d^rows.
    println!("-- weave (2D lattice, d=2): χ = 2^rows, EXPONENTIAL in the transverse extent --");
    println!(
        "{:>6} {:>6} {:>8} {:>10} {:>12} {:>10}",
        "rows", "sites", "χ", "2^rows", "params", "build"
    );
    for rows in [2usize, 3, 4, 5] {
        let (net, edges) = weave(2, rows, 6);
        let ((m, _), dt) = time(|| {
            let mut m = Mps::zero_state(&net.dims(), SPEC);
            net.cluster(&edges).run_mps(&mut m);
            (m, ())
        });
        println!(
            "{:>6} {:>6} {:>8} {:>10} {:>12} {:>10.0?}",
            rows,
            net.dims().len(),
            m.max_bond_dim(),
            1usize << rows,
            m.param_count(),
            dt
        );
    }
    println!("→ each added transverse row DOUBLES χ (and roughly the work): the area law.\n");

    // Haar scrambling: volume law, truncation blows up.
    println!("-- Haar-random brickwork on a wave  (volume law, χ saturates) --");
    println!("{:>10} {:>8} {:>16}", "χ cap", "reached", "discarded weight");
    let sw = zigzag::wave(2, 4, 2);
    let mut rng = Rng::new(7);
    let scramble = zigzag::brickwork_random(&sw, 3, &mut rng);
    for cap in [8usize, 16, 32, 64] {
        let mut m = Mps::zero_state(&sw, TruncSpec::new(cap, 0.0));
        scramble.run_mps(&mut m);
        println!(
            "{:>10} {:>8} {:>16.3e}",
            cap,
            m.max_bond_dim(),
            m.discarded_weight
        );
    }
    println!("→ generic dynamics fills any budget: truncation error stays large at every cap.\n");

    // ======================================================================
    // CROSS-CUTTING I — predict vs pay (the atlas)
    // ======================================================================
    println!("========== predict-vs-pay: the width atlas ==========\n");
    println!("The width of ×k is a number-theoretic object — computable WITHOUT the tensor.");
    let wave = zigzag::wave(2, 5, 4); // N ≈ 8.6e12
    println!("{:>10} {:>18} {:>18} {:>10}", "k", "predict (arithmetic)", "pay (build MPO)", "χ");
    for k in [7usize, 11] {
        let (pred, dt_pred) = time(|| {
            mult_width_profile(&wave, k as u128)
                .iter()
                .map(|w| w.value() as usize)
                .max()
                .unwrap()
        });
        let (m, dt_pay) = time(|| Transducer::mult(&wave, k).to_mpo(SPEC));
        assert_eq!(pred, m.max_bond_dim());
        println!(
            "{:>10} {:>18.0?} {:>18.0?} {:>10}",
            k, dt_pred, dt_pay, pred
        );
    }
    println!("→ predicting the output width is orders faster than producing it, and exact.\n");

    // ======================================================================
    // CROSS-CUTTING II — layout is the cost
    // ======================================================================
    println!("========== layout-dependence: the SAME state, two costs ==========\n");
    let cs = zigzag::diamond(1, 4);
    println!("{:>14} {:>8} {:>10} {:>12}", "layout", "χ", "params", "build");
    for (name, layout) in [("block", Layout::Block), ("interleaved", Layout::Interleaved)] {
        let x = Crossing::twin(&cs, true).with_layout(layout);
        let ((m, _), dt) = time(|| {
            let mut m = Mps::zero_state(&x.dims(), SPEC);
            x.bell(&[1, 2, 3, 4]).run_mps(&mut m);
            (m, ())
        });
        println!(
            "{:>14} {:>8} {:>10} {:>12.0?}",
            name,
            m.max_bond_dim(),
            m.param_count(),
            dt
        );
    }
    println!("→ identical 7.17 bits of entanglement; the cost is the site ordering (cutwidth).\n");

    // ======================================================================
    // FLOWS — width along a one-parameter family
    // ======================================================================
    println!("========== flows: width as a function of the flow parameter ==========\n");
    let fp = zigzag::diamond(2, 5);
    let flow = FourierFlow::adder(&fp, 1.0, SPEC);
    println!("FourierFlow U^t on Z_2880 — χ sees the integers:");
    print!("  t:    ");
    for k in 0..=8 {
        print!("{:>5.2}", k as f64 / 4.0);
    }
    println!();
    print!("  χ:    ");
    for k in 0..=8 {
        print!("{:>5}", flow.at(k as f64 / 4.0).max_bond_dim());
    }
    println!("\n→ χ collapses at whole shifts (crisp permutations), widens between.\n");

    let small = vec![2usize, 3, 2];
    let (frame, dt_f) = time(|| FrameFlow::qft(&small, SPEC));
    println!(
        "FrameFlow F^t on Z_12 (built in {:.0?}) — χ of the pure powers:",
        dt_f
    );
    println!(
        "  F^0=I: χ={}   F^1=F: χ={}   F^2=Π: χ={}   F^3=F†: χ={}",
        frame.powers[0].max_bond_dim(),
        frame.powers[1].max_bond_dim(),
        frame.powers[2].max_bond_dim(),
        frame.powers[3].max_bond_dim()
    );
    println!("→ the QFT flow is a 4-term combination; width sees F⁴ = I.\n");

    // ======================================================================
    // BOUNDARIES — convergence complexity of a dynamical system
    // ======================================================================
    println!("========== boundaries: iteration count to self-stabilize ==========\n");
    let bp = vec![2usize, 3, 4];
    let n_total: u128 = bp.iter().map(|&d| d as u128).product();
    let seam: Vec<usize> = bp.iter().map(|&d| d - 1).collect();
    let pump = Transducer::adder(&bp, 0).to_mpo_looped(SPEC);
    let bare = iterate(&pump, &basis(&bp, n_total - 1), 3000, 1e-12);
    println!(
        "  bare seam (Jordan block, polynomial 1/k²): {} steps to 1e-12",
        bare.converged_at.map(|s| s.to_string()).unwrap_or_else(|| ">3000".into())
    );
    for gamma in [0.5f64, 0.2, 0.05] {
        let damper = Mpo::identity(&bp, SPEC).add(
            &Mpo::basis_transfer(&bp, &seam, &seam, SPEC).scale(C64::real(-(1.0 - gamma))),
            SPEC,
        );
        let m = pump.compose_after(&damper, SPEC);
        let run = iterate(&m, &basis(&bp, n_total - 1), 500, 1e-12);
        println!(
            "  damped seam γ={:<4} (eigenvalue, exponential γ^k): {} steps",
            gamma,
            run.converged_at.map(|s| s.to_string()).unwrap_or_else(|| ">500".into())
        );
    }
    println!("→ the boundary sets the convergence CLASS: polynomial vs exponential.\n");

    // ======================================================================
    // The dense wall — what χ buys
    // ======================================================================
    println!("========== the headline: what bounded χ buys ==========\n");
    println!("{:>16} {:>8} {:>16} {:>14} {:>10}", "object", "sites", "dense dim", "MPS/MPO size", "ratio");
    let w = zigzag::wave(2, 5, 4);
    let dense = zigzag::dense_dimension(&w);
    let adder = radix::adder_mpo_exact(&w, 1, SPEC);
    println!(
        "{:>16} {:>8} {:>16.2e} {:>11.1} KB {:>10.0e}",
        "+1 adder",
        w.len(),
        dense * dense, // operator on N-dim space
        kb(adder.param_count()),
        (dense * dense) / (adder.param_count() as f64)
    );
    println!(
        "\nThe 25-site wave's adder acts on a {:.1e}-dimensional operator space and is",
        dense * dense
    );
    println!("stored in {:.1} KB — the entire point of a bounded-χ representation.", kb(adder.param_count()));
}
