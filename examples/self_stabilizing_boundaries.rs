//! Self-stabilizing boundary systems: closing the message loop turns a
//! machine into a dynamical system with a designed attractor.
//!
//! Everything so far kept the message front's boundaries *open* (fixed
//! entry, dropped or postselected exit). Here we **close the loop** — feed
//! the exiting message back into the entry — and study the resulting
//! (non-unitary) operator by iterating it. A unitary map can't
//! self-stabilize; the non-unitarity a boundary introduces is exactly the
//! resource that makes an attractor.
//!
//! Four findings, on the diamond's ring unless noted:
//!
//! 1. **Closing the loop changes the arithmetic.** The traced adder is
//!    end-around carry — ones'-complement — computing mod `N-1`, so the
//!    diamond's composite `Z_2880` becomes the *prime field* `Z_2879`.
//! 2. **A prime field has no gcd obstruction.** `×6` is 6-to-1 (non-unitary)
//!    on `Z_2880`; the looped machine is exactly unitary on `Z_2879`. The
//!    boundary *heals* the operator.
//! 3. **The seam self-stabilizes.** The two representatives of zero
//!    (`0` and `N-1`) make the traced identity the defective operator
//!    `I + |0⟩⟨N-1|`; iterating it pumps negative-zero to canonical zero
//!    with error exactly `1/(k²+1)` — polynomial, from a Jordan block.
//! 4. **The convergence law is designable.** Composing a damped seam
//!    `I − (1−γ)|N-1⟩⟨N-1|` turns the Jordan block into a genuine
//!    eigenvalue `γ`, and stabilization becomes exponential at a rate you
//!    choose.
//!
//! Run with: `cargo run --release --example self_stabilizing_boundaries`

use novel_quantum_structures::c64::C64;
use novel_quantum_structures::cascade::Transducer;
use novel_quantum_structures::dense::total_dim;
use novel_quantum_structures::mpo::Mpo;
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::stabilize::{fixed_point_residual, iterate};
use novel_quantum_structures::{zigzag, TruncSpec};

const SPEC: TruncSpec = TruncSpec {
    max_rank: 128,
    cutoff: 1e-13,
};

fn basis(profile: &[usize], mut x: u128) -> Mps {
    let mut d = vec![0usize; profile.len()];
    for (i, &dim) in profile.iter().enumerate().rev() {
        d[i] = (x % dim as u128) as usize;
        x /= dim as u128;
    }
    Mps::basis_state(profile, &d, SPEC)
}

fn main() {
    let profile = zigzag::diamond(2, 5);
    let n = total_dim(&profile) as u128;
    println!(
        "diamond {:?}: open ring Z_{}, looped ring Z_{} ({} = 2879 is prime)\n",
        profile,
        n,
        n - 1,
        n - 1
    );

    // ---- 1. The loop changes the arithmetic -------------------------------
    println!("=== 1. Closing the loop = end-around carry = mod N-1 ===\n");
    let add5_open = Transducer::adder(&profile, 5).to_mpo(SPEC);
    let add5_loop = Transducer::adder(&profile, 5).to_mpo_looped(SPEC);
    for x in [10u128, 2873, 2878] {
        let open = (x + 5) % n;
        let loop_t = (x + 5) % (n - 1);
        let f_open = add5_open
            .apply_to(&basis(&profile, x))
            .fidelity(&basis(&profile, open));
        let f_loop = add5_loop
            .apply_to(&basis(&profile, x))
            .fidelity(&basis(&profile, loop_t));
        println!(
            "  x={:>4}:  open → {:>4} (fid {:.6}),  looped → {:>4} mod 2879 (fid {:.6})",
            x, open, f_open, loop_t, f_loop
        );
    }
    println!(
        "  looped bond profile {:?} (vs open {:?})\n",
        add5_loop.bond_dims(),
        add5_open.bond_dims()
    );

    // ---- 2. The prime field heals the gcd obstruction ---------------------
    println!("=== 2. The boundary heals a non-unitary operator ===\n");
    let m6_open = Transducer::mult(&profile, 6).to_mpo(SPEC);
    let m6_loop = Transducer::mult(&profile, 6).to_mpo_looped(SPEC);
    let defect = |m: &Mpo| {
        let g = m.adjoint().compose_after(m, SPEC);
        1.0 - g.hs_fidelity(&Mpo::identity(&profile, SPEC))
    };
    println!(
        "  ×6 on Z_2880 (gcd=6):   unitarity defect {:.4}  (6-to-1, not invertible)",
        defect(&m6_open)
    );
    println!(
        "  ×6 on Z_2879 (prime):   unitarity defect {:.1e}  (bijective — healed)",
        defect(&m6_loop)
    );
    // Confirm it is genuinely ×6 mod 2879.
    let ok = [1u128, 800, 2000].iter().all(|&x| {
        m6_loop
            .apply_to(&basis(&profile, x))
            .fidelity(&basis(&profile, 6 * x % (n - 1)))
            > 1.0 - 1e-7
    });
    println!("  verified |x⟩ → |6x mod 2879⟩ on samples: {}\n", ok);

    // ---- 3. The seam self-stabilizes (polynomial) -------------------------
    println!("=== 3. The double-zero seam self-stabilizes ===\n");
    let pump = Transducer::adder(&profile, 0).to_mpo_looped(SPEC);
    println!("  traced identity = I + |0⟩⟨N-1|  (both 0 and 2879 represent zero)");
    println!(
        "  ordinary state |1234⟩ is already a fixed point: residual {:.1e}",
        fixed_point_residual(&pump, &basis(&profile, 1234))
    );
    println!("  iterating from negative-zero |2879⟩ toward canonical |0⟩:\n");
    println!("  {:>4} {:>16} {:>16}", "step", "error", "1/(k²+1) law");
    let mut psi = basis(&profile, n - 1);
    psi.normalize();
    for k in 1..=8 {
        psi = pump.apply_to(&psi);
        psi.normalize();
        let err = 1.0 - psi.fidelity(&basis(&profile, 0));
        let law = 1.0 / ((k * k) as f64 + 1.0);
        println!("  {:>4} {:>16.6e} {:>16.6e}", k, err, law);
    }
    println!("  → error = 1/(k²+1) exactly: a Jordan block, polynomial stabilization\n");

    // ---- 4. Designable convergence with a damped seam ---------------------
    println!("=== 4. A damped seam makes stabilization exponential ===\n");
    let seam: Vec<usize> = profile.iter().map(|&d| d - 1).collect();
    println!("  damper = I − (1−γ)|N-1⟩⟨N-1|  (bond-1 product operator)");
    println!("  {:>8} {:>16} {:>20}", "γ", "steps to 1e-12", "regime");
    let bare = iterate(&pump, &basis(&profile, n - 1), 5000, 1e-12);
    println!(
        "  {:>8} {:>16} {:>20}",
        "1.0 (none)",
        bare.converged_at
            .map(|s| s.to_string())
            .unwrap_or_else(|| ">5000".into()),
        "polynomial"
    );
    for gamma in [0.5f64, 0.2, 0.05] {
        let damper = Mpo::identity(&profile, SPEC).add(
            &Mpo::basis_transfer(&profile, &seam, &seam, SPEC).scale(C64::real(-(1.0 - gamma))),
            SPEC,
        );
        let m = pump.compose_after(&damper, SPEC);
        let run = iterate(&m, &basis(&profile, n - 1), 500, 1e-12);
        println!(
            "  {:>8} {:>16} {:>20}",
            gamma,
            run.converged_at
                .map(|s| s.to_string())
                .unwrap_or_else(|| ">500".into()),
            "exponential"
        );
    }

    println!("\nA closed boundary loop turns a machine into a dynamical system. The");
    println!("loop chooses the ring (composite → prime, healing obstructions), the");
    println!("seam is an attractor that repairs the number representation, and a");
    println!("rank-one damper sets the convergence rate. Self-stabilization is a");
    println!("property you engineer at the boundary, not inside the operator.");
}
