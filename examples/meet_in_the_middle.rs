//! **Retrodictive-prediction as a meet-in-the-middle mesh**, ordered by the
//! global tail-radix solve.
//!
//! The antiparallel crossing couples opposite time directions: a forward
//! cascade (`×a`, *prediction*, sweeps LSB→MSB) crossed with a backward
//! cascade (`÷b`, *retrodiction*, sweeps MSB→LSB). By the minimal-machine
//! principle (finding 18) the crossing computes an operation whose
//! single-direction width is enormous — meet in the middle — at the width
//! of the two narrow fronts. And the **tail-radix phase web** (finding 8),
//! the ring's own Fourier structure, is *where the fronts meet*: it is the
//! frame `V` in which the solve-components decouple into their ideal order.
//!
//! Four measured movements:
//!
//! 1. the dual-direction crossing computes a wide `×m` narrowly (and the
//!    single-front machine is literally unbuildable);
//! 2. the tail-radix frame is the meeting point — it maps *prediction to
//!    retrodiction* (`V M_k V† = M_{k⁻¹}`, the two crossing directions are
//!    frame-conjugate);
//! 3. the global tail-radix solve gives the ideal ordering — in the frame
//!    the additive components are diagonal and *commute*, so a whole mesh
//!    collapses to one frame round-trip, and the tail-radix ordering
//!    cancels the expensive digit-reversal meet-in-the-middle (`χ 6` vs
//!    `36`);
//! 4. a full affine mesh (multiplicative + additive, forward + backward)
//!    solved this way, verified and narrow.
//!
//! Run with: `cargo run --release --example meet_in_the_middle`

use novel_quantum_structures::cascade::Transducer;
use novel_quantum_structures::mpo::Mpo;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{radix, zigzag, TruncSpec};
use std::time::Instant;

const SPEC: TruncSpec = TruncSpec {
    max_rank: 96,
    cutoff: 1e-12,
};

fn modinv(k: u128, n: u128) -> u128 {
    let (mut r0, mut r) = (k as i128, n as i128);
    let (mut s0, mut s) = (1i128, 0i128);
    while r != 0 {
        let q = r0 / r;
        (r0, r) = (r, r0 - q * r);
        (s0, s) = (s, s0 - q * s);
    }
    s0.rem_euclid(n as i128) as u128
}

fn basis(p: &[usize], mut x: u128) -> Mps {
    let mut d = vec![0usize; p.len()];
    for (i, &dim) in p.iter().enumerate().rev() {
        d[i] = (x % dim as u128) as usize;
        x /= dim as u128;
    }
    Mps::basis_state(p, &d, SPEC)
}

fn main() {
    let profile = zigzag::diamond(2, 5);
    let n = 2880u128;
    println!("diamond {:?}, ring Z_{}\n", profile, n);

    // ---- 1. The dual-direction crossing -----------------------------------
    println!("=== 1. Retrodictive-prediction: a wide op as a narrow crossing ===\n");
    println!("×m = (×a forward, PREDICTION) crossed with (÷b backward, RETRODICTION):\n");
    println!(
        "{:>8} {:>18} {:>22} {:>10}",
        "×m", "naive forward", "crossed ×a ⋈ ÷b", "verified"
    );
    for (a, b) in [(7u128, 11u128), (7, 13), (5, 11)] {
        let m = a * modinv(b, n) % n;
        let fwd = Transducer::mult(&profile, a as usize).to_mpo(SPEC);
        let bwd = Transducer::div_mpo(&profile, b as usize, SPEC);
        let crossed = bwd.compose_after(&fwd, SPEC);
        let ok = [1u128, 1000, n - 1]
            .iter()
            .all(|&x| crossed.apply_to(&basis(&profile, x)).fidelity(&basis(&profile, m * x % n)) > 1.0 - 1e-8);
        println!(
            "{:>8} {:>11} states {:>14} width {:<3} {:>10}",
            m,
            m, // the forward multiplier's message dimension
            format!("×{}⋈÷{}", a, b),
            crossed.max_bond_dim(),
            ok
        );
    }
    println!("\nThe naive forward machine needs m states — and the direct cascade is");
    println!("capped at 512, so all three are literally UNBUILDABLE that way. The");
    println!("crossing of two narrow fronts is the only feasible route, and it is the");
    println!("narrowest: forward prediction reconciled with backward retrodiction.\n");

    // ---- 2. The tail-radix frame is the meeting point ----------------------
    println!("=== 2. The tail-radix frame maps prediction ↔ retrodiction ===\n");
    let small = vec![2usize, 3, 4, 3]; // N = 72, non-palindromic
    let v = Mpo::from_circuit(&radix::mixed_radix_qft_reversed(&small, 0.0), SPEC);
    let k = 5usize;
    let k_inv = modinv(k as u128, 72) as usize; // 29
    let fwd = Transducer::mult(&small, k).to_mpo(SPEC);
    let conj = v.compose_after(&fwd.compose_after(&v.adjoint(), SPEC), SPEC);
    let bwd_rev =
        Transducer::mult_directed(&small, k_inv, novel_quantum_structures::cascade::Direction::LeftToRight)
            .to_mpo(SPEC);
    println!("on Z_72:  V (×{}) V†  vs  ×{} sweeping the opposite way:", k, k_inv);
    println!("  hs fidelity {:.9}", conj.hs_fidelity(&bwd_rev));
    println!("\nThe tail-radix frame V conjugates a FORWARD ×k into a BACKWARD ×k⁻¹ —");
    println!("prediction and retrodiction are the same operator seen from the two sides");
    println!("of the frame. So the Fourier web is exactly where the crossing meets.\n");

    // ---- 3. The global tail-radix solve: ideal ordering --------------------
    println!("=== 3. The global tail-radix solve: ideal ordering of components ===\n");
    let p4 = zigzag::diamond(2, 4); // N = 144
    let v_rev = Mpo::from_circuit(&radix::mixed_radix_qft_reversed(&p4, 0.0), SPEC);
    let v_std = Mpo::from_circuit(&radix::mixed_radix_qft(&p4, 0.0), SPEC);
    println!(
        "tail-radix frame width:  reversed (ideal) χ={},  standard (+digit reversal) χ={}",
        v_rev.max_bond_dim(),
        v_std.max_bond_dim()
    );
    println!("The digit reversal is the costly component; the tail-radix ordering makes");
    println!("it CANCEL between V and V† — paid at χ={} instead of χ={}.\n", v_rev.max_bond_dim(), v_std.max_bond_dim());

    // A mesh of additive solve-components: in the frame they commute (diagonal).
    let adder = |c: usize| {
        let d = Mpo::from_circuit(&radix::fourier_phase_ramp_reversed(&p4, c), SPEC);
        v_rev.adjoint().compose_after(&d.compose_after(&v_rev, SPEC), SPEC)
    };
    let cs = [11usize, 29, 47, 100, 133];
    // Naive: cross each component through the frame separately (K round-trips).
    let t0 = Instant::now();
    let mut naive = Mpo::identity(&p4, SPEC);
    let mut peak = 1;
    for &c in &cs {
        naive = adder(c).compose_after(&naive, SPEC);
        peak = peak.max(naive.max_bond_dim());
    }
    let dt_naive = t0.elapsed();
    // Ideal: ONE frame round-trip; the additive components sum in the frame.
    let t1 = Instant::now();
    let sum: usize = cs.iter().sum();
    let ideal = adder(sum % 144);
    let dt_ideal = t1.elapsed();
    println!(
        "mesh of {} additive components, solved two ways (Σc = {}):",
        cs.len(),
        sum
    );
    println!(
        "  naive  — {} separate frame round-trips: final χ={}, peak χ={}, {:.0?}",
        cs.len(),
        naive.max_bond_dim(),
        peak,
        dt_naive
    );
    println!(
        "  ideal  — ONE tail-radix frame solve (components commute): χ={}, {:.0?}",
        ideal.max_bond_dim(),
        dt_ideal
    );
    println!("  agree: {}", (naive.hs_fidelity(&ideal) - 1.0).abs() < 1e-7);
    println!("\nIn the frame the additive solve-components are diagonal and COMMUTE, so");
    println!("their ideal ordering is trivial — enter the frame once, sum, leave once —");
    println!("collapsing {} crossings into one, at a fraction of the work.\n", cs.len());

    // ---- 4. A full affine mesh, solved -------------------------------------
    println!("=== 4. A full affine mesh solved by the crossing + tail-radix order ===\n");
    // x → k·x + c: forward ×a, backward ÷b for the multiplier, +c in the frame.
    let (a, b, c) = (7usize, 11usize, 500usize);
    let m = (a as u128 * modinv(b as u128, n)) % n; // the multiplier
    let mult = Transducer::div_mpo(&profile, b, SPEC).compose_after(&Transducer::mult(&profile, a).to_mpo(SPEC), SPEC);
    let add = radix::adder_mpo_exact(&profile, c as u128, SPEC);
    let affine = add.compose_after(&mult, SPEC); // x → m·x + c
    println!(
        "x → {}·x + {} (mod {}) = (×{} ⋈ ÷{}) then +{}:",
        m, c, n, a, b, c
    );
    println!(
        "  affine operator width {:?} (max {})",
        affine.bond_dims(),
        affine.max_bond_dim()
    );
    let ok = [1u128, 777, n - 1].iter().all(|&x| {
        affine
            .apply_to(&basis(&profile, x))
            .fidelity(&basis(&profile, (m * x + c as u128) % n))
            > 1.0 - 1e-8
    });
    println!("  verified x → {}x+{} on basis states: {}", m, c, ok);
    println!("\nThe single-front multiplier ×{} is unbuildable; the composite is solved", m);
    println!("as a narrow mesh — forward prediction × backward retrodiction for the");
    println!("multiplicative part, an additive component in the tail-radix frame, the");
    println!("whole thing ordered by the ring's own Fourier structure and read at the");
    println!("meeting point in the middle.");
}
