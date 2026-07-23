//! Self-stabilization dynamics of boundary systems: iterate a
//! (generally non-unitary) operator and watch it converge.
//!
//! Unitary machines cannot self-stabilize — they preserve distances, so
//! nothing attracts. The non-unitarity that boundary systems introduce
//! (postselected exits, looped message fronts, damped rank-one seams) is
//! exactly the resource stabilization needs: iterating
//! `|ψ⟩ → M|ψ⟩ / ‖M|ψ⟩‖` flows toward the dominant eigenspace of `M` — the
//! **attractor selected by the boundary conditions**.
//!
//! [`iterate`] runs the flow and records its diagnostics: the per-step norm
//! ratio (the contraction the boundary exerts) and the overlap with the
//! previous iterate (convergence of the direction). Two convergence
//! regimes appear in `examples/self_stabilizing_boundaries.rs`:
//!
//! * a **defective** (Jordan-block) boundary seam — e.g. the traced
//!   identity `I + |0⟩⟨N-1|` — stabilizes *polynomially* (error `~ 1/k²`);
//! * composing a damped seam `I − (1−γ)|N-1⟩⟨N-1|` makes the bad direction
//!   a genuine eigenvalue `γ < 1` and convergence *exponential* at rate γ.
//!
//! The boundary system is a design surface: it selects the attractor and
//! sets the convergence law.

use crate::mpo::Mpo;
use crate::mps::Mps;

/// Diagnostics of one stabilization run.
pub struct StabilizeRun {
    /// Normalized state after the last step taken.
    pub final_state: Mps,
    /// `‖M ψ_k‖ / ‖ψ_k‖` per step — the boundary's contraction/amplification.
    pub norm_ratio: Vec<f64>,
    /// `|⟨ψ_k | ψ_{k+1}⟩|²` per step — directional convergence.
    pub overlap_prev: Vec<f64>,
    /// First step at which `1 − overlap_prev < tol`, if reached.
    pub converged_at: Option<usize>,
}

/// Iterate `|ψ⟩ → M|ψ⟩/‖M|ψ⟩‖` from `psi0` for up to `max_steps`, stopping
/// early once the direction settles to within `tol`.
pub fn iterate(m: &Mpo, psi0: &Mps, max_steps: usize, tol: f64) -> StabilizeRun {
    let mut psi = psi0.clone();
    psi.normalize();
    let mut norm_ratio = Vec::new();
    let mut overlap_prev = Vec::new();
    let mut converged_at = None;
    for step in 0..max_steps {
        let mut next = m.apply_to(&psi);
        let ratio = next.norm();
        assert!(ratio > 1e-300, "state annihilated by the boundary system");
        next.normalize();
        let ov = next.fidelity(&psi);
        norm_ratio.push(ratio);
        overlap_prev.push(ov);
        psi = next;
        if converged_at.is_none() && 1.0 - ov < tol {
            converged_at = Some(step + 1);
            break;
        }
    }
    StabilizeRun {
        final_state: psi,
        norm_ratio,
        overlap_prev,
        converged_at,
    }
}

/// `1 − fidelity(M ψ/‖M ψ‖, ψ)`: zero iff `ψ` is a fixed direction of `M`.
pub fn fixed_point_residual(m: &Mpo, psi: &Mps) -> f64 {
    let mut out = m.apply_to(psi);
    out.normalize();
    let mut p = psi.clone();
    p.normalize();
    (1.0 - out.fidelity(&p)).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cascade::Transducer;
    use crate::dense::total_dim;
    use crate::mat::TruncSpec;
    use crate::mpo::Mpo;

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

    #[test]
    fn traced_identity_pumps_the_double_zero() {
        // M = I + |0⟩⟨N-1|: iterating from |N-1⟩ converges to |0⟩ —
        // the boundary seam self-stabilizes the number representation.
        let profile = vec![2usize, 3, 4]; // N = 24
        let n_total = total_dim(&profile) as u128;
        let m = Transducer::adder(&profile, 0).to_mpo_looped(SPEC);
        let run = iterate(&m, &basis(&profile, n_total - 1), 200, 1e-10);
        let f0 = run.final_state.fidelity(&basis(&profile, 0));
        assert!(f0 > 0.999, "converged to |0⟩: fidelity {}", f0);
        // Ordinary states are already fixed points.
        assert!(fixed_point_residual(&m, &basis(&profile, 7)) < 1e-10);
    }

    #[test]
    fn damped_seam_converges_exponentially() {
        use crate::c64::C64;
        let profile = vec![2usize, 3, 4]; // N = 24
        let n_total = total_dim(&profile);
        let seam_digits: Vec<usize> = profile.iter().map(|&d| d - 1).collect();
        let gamma = 0.5;
        let damper = Mpo::identity(&profile, SPEC).add(
            &Mpo::basis_transfer(&profile, &seam_digits, &seam_digits, SPEC)
                .scale(C64::real(-(1.0 - gamma))),
            SPEC,
        );
        let pump = Transducer::adder(&profile, 0).to_mpo_looped(SPEC);
        let m = pump.compose_after(&damper, SPEC);

        let run = iterate(&m, &basis(&profile, n_total as u128 - 1), 60, 1e-12);
        let f0 = run.final_state.fidelity(&basis(&profile, 0));
        assert!(f0 > 1.0 - 1e-9, "converged to |0⟩: fidelity {}", f0);
        // Exponential regime: it must settle far faster than the bare pump
        // does at the same tolerance.
        let bare = iterate(&pump, &basis(&profile, n_total as u128 - 1), 2000, 1e-12);
        let fast = run.converged_at.expect("damped run converged");
        // If the bare pump converged at all (polynomial → slow), the damped
        // run must be far faster; if it didn't converge in 2000 steps, the
        // contrast is already established.
        if let Some(slow) = bare.converged_at {
            assert!(fast * 4 < slow, "damped {} vs bare {}", fast, slow);
        }
    }
}
