//! Boundary twists: one machine body, a family of rings selected purely at
//! the boundary.
//!
//! A cascade's message loop can be closed straight (`m → m`) or **twisted**
//! through a permutation of the message set. The two clean twists of the
//! carry machines complete a ring family:
//!
//! ```text
//!   open  (enter 0, drop exit)      →  arithmetic mod N
//!   loop  (σ = id, end-around)      →  arithmetic mod N−1   (double zero)
//!   twist (σ = reversal)            →  arithmetic mod N+1   (missing zero)
//! ```
//!
//! The reversal twist realizes **diminished-one arithmetic** — the
//! encoding of Fermat-number-transform hardware: chain value `x`
//! represents `v = x+1 ∈ [1, N] ⊂ Z_{N+1}`, the twisted `+c` machine maps
//! `v → v + (c+1) mod N+1`, the twisted `×k` machine maps `v → k·v mod
//! N+1` — and the unrepresentable zero of `Z_{N+1}` appears as an
//! *annihilated branch*: a hole, the exact dual of the straight loop's
//! double zero. Unitarity on each closed ring is governed by gcd against
//! *that* ring — so a twist can heal an obstruction or create one, and
//! primality of `N∓1` is the design criterion.
//!
//! Run with: `cargo run --release --example boundary_twists`

use novel_quantum_structures::cascade::{unitarity_defect, Transducer};
use novel_quantum_structures::dense::{total_dim, DenseState};
use novel_quantum_structures::mps::Mps;
use novel_quantum_structures::{zigzag, TruncSpec};

const SPEC: TruncSpec = TruncSpec {
    max_rank: 96,
    cutoff: 1e-12,
};

fn basis(profile: &[usize], mut x: u128) -> Mps {
    let mut d = vec![0usize; profile.len()];
    for (i, &dim) in profile.iter().enumerate().rev() {
        d[i] = (x % dim as u128) as usize;
        x /= dim as u128;
    }
    Mps::basis_state(profile, &d, SPEC)
}

fn reversal(m: usize) -> Vec<usize> {
    (0..m).rev().collect()
}

fn main() {
    let profile = vec![2usize, 3, 4, 3];
    let n = total_dim(&profile) as u128; // 72
    println!(
        "chain {:?}, N = {} — flanked by the twin primes {} and {}\n",
        profile,
        n,
        n - 1,
        n + 1
    );

    // ---- 1. One machine, three rings ---------------------------------------
    println!("=== 1. The +5 machine under its three boundary conditions ===\n");
    let c = 5u128;
    let open = Transducer::adder(&profile, c).to_mpo(SPEC);
    let looped = Transducer::adder(&profile, c).to_mpo_looped(SPEC);
    let twisted = Transducer::adder(&profile, c).to_mpo_looped_twisted(&[1, 0], SPEC);
    println!("x → x+5 open (mod 72), looped (mod 71), twisted (diminished-one mod 73):");
    println!(
        "{:>6} {:>12} {:>12} {:>22}",
        "x", "open", "looped", "twisted (v = x+1)"
    );
    for x in [10u128, 65, 70, 71] {
        let f_of = |m: &novel_quantum_structures::mpo::Mpo, target: u128| {
            m.apply_to(&basis(&profile, x)).fidelity(&basis(&profile, target))
        };
        let t_open = (x + c) % n;
        let t_loop = (x + c) % (n - 1);
        let v_out = (x + 1 + c + 1) % (n + 1);
        println!(
            "{:>6} {:>4} ({:.4}) {:>4} ({:.4}) {:>7} ({:.4})",
            x,
            t_open,
            f_of(&open, t_open),
            t_loop,
            f_of(&looped, t_loop),
            v_out - 1,
            f_of(&twisted, v_out - 1)
        );
    }
    // The defect input is the SAME x for both closures: x + c = N−1 means
    // the loop lands on its double zero and the twist on its missing zero.
    let x_defect = n - 1 - c;
    let seam_out = looped.apply_to(&basis(&profile, x_defect));
    let hole_out = twisted.apply_to(&basis(&profile, x_defect));
    println!("\nthe aligned defect, x = {} (x + 5 = 71):", x_defect);
    println!(
        "  looped:  BOTH zeros appear — ⟨0|M|x⟩ = ⟨71|M|x⟩ = 1, ‖M|x⟩‖ = {:.4} (the seam)",
        seam_out.norm()
    );
    println!(
        "  twisted: the branch is ANNIHILATED — ‖M|x⟩‖ = {:.1e} (the hole)",
        hole_out.norm()
    );
    println!("A double zero and a missing zero are the two ways a boundary can");
    println!("misalign the N-state chain with its chosen (N∓1)-element ring.\n");

    // ---- 2. Twin primes: healing on both sides -----------------------------
    println!("=== 2. gcd healing on both flanks ===\n");
    let k = 6usize; // gcd(6, 72) = 6: badly non-invertible on the open ring
    let m_open = Transducer::mult(&profile, k).to_mpo(SPEC);
    let m_loop = Transducer::mult(&profile, k).to_mpo_looped(SPEC);
    let m_twist = Transducer::mult(&profile, k).to_mpo_looped_twisted(&reversal(k), SPEC);
    println!("×6:  open  (Z_72, gcd 6):  unitarity defect {:.4}", unitarity_defect(&m_open, SPEC));
    println!("     loop  (Z_71, prime):  unitarity defect {:.1e}", unitarity_defect(&m_loop, SPEC));
    println!("     twist (Z_73, prime):  unitarity defect {:.1e}", unitarity_defect(&m_twist, SPEC));
    for x in [0u128, 35, 71] {
        let v_out = (k as u128 * (x + 1)) % (n + 1);
        let f = m_twist
            .apply_to(&basis(&profile, x))
            .fidelity(&basis(&profile, v_out - 1));
        assert!((f - 1.0).abs() < 1e-8);
    }
    println!("     twist verified as v → 6v mod 73 on samples.");
    println!("\nN = 72 sits between twin primes: BOTH closed boundary rings are fields,");
    println!("so every multiplier the open chain breaks is healed by either closure.\n");

    // ---- 3. The twist can also break ---------------------------------------
    println!("=== 3. Each ring has its own obstructions ===\n");
    let p3 = vec![2usize, 3, 4]; // N = 24: twist ring 25 = 5²
    let m5_open = Transducer::mult(&p3, 5).to_mpo(SPEC);
    let m5_twist = Transducer::mult(&p3, 5).to_mpo_looped_twisted(&reversal(5), SPEC);
    println!("×5 on {:?} (N = 24, twist ring 25 = 5²):", p3);
    println!("  open:  defect {:.1e}   (gcd(5, 24) = 1 — unitary)", unitarity_defect(&m5_open, SPEC));
    println!("  twist: defect {:.3}   (gcd(5, 25) = 5 — holes at v ≡ 0 mod 5)", unitarity_defect(&m5_twist, SPEC));
    for x in [4u128, 9, 14, 19] {
        println!(
            "    |{}⟩ (v = {:>2}): ‖M|x⟩‖ = {:.1e}",
            x,
            x + 1,
            m5_twist.apply_to(&basis(&p3, x)).norm()
        );
    }
    let diamond = zigzag::diamond(2, 5);
    let dn = total_dim(&diamond) as u128;
    println!(
        "\nOn the diamond (N = {}): loop ring {} is prime, twist ring {} = 43·67.",
        dn,
        dn - 1,
        dn + 1
    );
    // The twisted ×43 operator is Σ_m ⟨σ(m)|T|m⟩; its defects show up
    // branch-by-branch (for each input at most one branch survives, and
    // this is exact — no summed operator to truncate):
    use novel_quantum_structures::cascade::Boundary;
    let t43 = Transducer::mult(&diamond, 43);
    let sigma43 = reversal(43);
    println!("twisted ×43, branch by branch:");
    for x in [0u128, 67, 66] {
        let mut survivor: Option<(usize, u128)> = None;
        for (m, &sm) in sigma43.iter().enumerate() {
            let branch =
                t43.to_mpo_with(Boundary::Fixed(m), Boundary::Fixed(sm), SPEC);
            let res = branch.apply_to(&basis(&diamond, x));
            if res.norm() > 0.5 {
                let dense = res.to_dense();
                let img = dense
                    .amps
                    .iter()
                    .position(|a| a.abs() > 0.5)
                    .expect("basis image") as u128;
                survivor = Some((m, img));
                break;
            }
        }
        match survivor {
            Some((m, img)) => println!(
                "  |{:>2}⟩ (v = {:>2}): branch m = {:>2} survives → |{}⟩",
                x,
                x + 1,
                m,
                img
            ),
            None => println!(
                "  |{:>2}⟩ (v = {:>2}): NO branch survives — the hole (43·67 ≡ 0 mod 2881)",
                x,
                x + 1
            ),
        }
    }
    println!("|0⟩ and |67⟩ collide on |42⟩ (43·1 ≡ 43·68 ≡ 43 mod 2881), and v = 67");
    println!("is annihilated: the twist ring's factorization surfaces as collisions");
    println!("plus holes.\n");

    // ---- 4. Seam vs hole: pump vs drain ------------------------------------
    println!("=== 4. The two traced identities: a pump and a drain ===\n");
    let drain = Transducer::adder(&p3, 0).to_mpo_looped_twisted(&[1, 0], SPEC);
    println!("loop  of +0: I + |0⟩⟨N−1|      — merges the double zero (a pump)");
    println!("twist of +0: Σ_x<N−1 |x+1⟩⟨x|  — the unilateral shift (a drain):");
    println!("             everything conveys toward the missing zero and annihilates.\n");
    let n3 = total_dim(&p3) as f64;
    let mut uniform = DenseState::zero_state(&p3);
    let amp = novel_quantum_structures::C64::real(1.0 / n3.sqrt());
    for a in &mut uniform.amps {
        *a = amp;
    }
    let mut psi = Mps::from_dense(&uniform, SPEC);
    println!("drain applied to the uniform state — measured vs the exact law (N−k)/N:");
    println!("{:>6} {:>14} {:>14}", "k", "‖M^k u‖²", "(N−k)/N");
    for k_step in 1..=6 {
        psi = drain.apply_to(&psi);
        let nn = psi.norm();
        println!(
            "{:>6} {:>14.6} {:>14.6}",
            k_step,
            nn * nn,
            (n3 - k_step as f64) / n3
        );
    }
    println!("\nA seam repairs the representation (polynomial or damped-exponential");
    println!("stabilization, README finding 19); a hole removes what the ring cannot");
    println!("represent. Boundary design chooses which — and on which ring — without");
    println!("touching a single local rule of the machine.");
}
