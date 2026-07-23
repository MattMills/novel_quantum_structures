//! Operator flows and the **operation-width cursor**: operators as
//! state-components over a continuously divisible extent of computation.
//!
//! A circuit block usually exists at one granularity: the whole unitary `U`.
//! This module treats a block instead as a **one-parameter flow**
//! `t ↦ U^t` — the geodesic through operator space from the identity
//! (`t = 0`) to the full operation (`t = 1`) and beyond — with every
//! snapshot a first-class [`Mpo`]. Fractional powers are computed exactly in
//! a diagonalizing Fourier frame:
//!
//! ```text
//!   U^t = V† · D(t) · V,      D(t) a bond-1 diagonal, linear in t
//! ```
//!
//! so the flow is an exact one-parameter group (`U^s ∘ U^t = U^{s+t}`) at
//! MPO precision. For the ring translation (the modular adder) this is the
//! continuous shift `x → x + t·c`: at integer `t·c` a crisp permutation, at
//! fractional values a delocalized sinc kernel — and the operator's bond
//! dimension *measures* that difference (see
//! `examples/operator_width_cursor.rs`: operator width sees the integers).
//!
//! On top of the flow sits the [`WidthCursor`]: a pipeline of segments
//! `[t_0,t_1), [t_1,t_2), …`, each holding the collapsed operator for its
//! span. Segments **refine** (split at their midpoint, both halves imputed
//! from the flow — the geometric interpolation consistent with the coarse
//! block) and **coarsen** (compose back into a wider block), while the total
//! composition is invariant. Refining repeatedly at a focus point yields an
//! adaptive temporal resolution — coarse far from the cursor, increasingly
//! fine at it.

use crate::mat::TruncSpec;
use crate::mpo::Mpo;
use crate::mps::Mps;
use crate::radix;

/// A one-parameter operator flow `t ↦ V† · D(t·c) · V` over a qudit chain's
/// ring, realized in the reversed-digit Fourier frame. `at(1.0)` is the
/// modular adder by `c`; `at(t)` for fractional `t` is the continuous
/// translation by `t·c` ring units.
pub struct FourierFlow {
    pub profile: Vec<usize>,
    /// Shift per unit flow parameter, in ring units.
    pub c: f64,
    pub trunc: TruncSpec,
    v: Mpo,
    v_dag: Mpo,
}

impl FourierFlow {
    /// Build the flow of the modular adder by `c` (compiles the Fourier
    /// frame once; every subsequent snapshot costs two compositions).
    pub fn adder(profile: &[usize], c: f64, trunc: TruncSpec) -> FourierFlow {
        let v = Mpo::from_circuit(&radix::mixed_radix_qft_reversed(profile, 0.0), trunc);
        let v_dag = v.adjoint();
        FourierFlow {
            profile: profile.to_vec(),
            c,
            trunc,
            v,
            v_dag,
        }
    }

    /// The diagonal generator block `D(t·c)` (always bond dimension 1).
    pub fn generator_at(&self, t: f64) -> Mpo {
        Mpo::from_circuit(
            &radix::fourier_phase_ramp_reversed_real(&self.profile, t * self.c),
            self.trunc,
        )
    }

    /// The flow snapshot `U^t` as a collapsed operator.
    pub fn at(&self, t: f64) -> Mpo {
        let d = self.generator_at(t);
        self.v_dag
            .compose_after(&d.compose_after(&self.v, self.trunc), self.trunc)
    }

    /// The operator for a span `[t0, t1)` of the flow (`= U^{t1-t0}`).
    pub fn span(&self, t0: f64, t1: f64) -> Mpo {
        self.at(t1 - t0)
    }

    /// Apply the snapshot `U^t` to a state without materializing the
    /// collapsed operator (three block applications through the frame).
    pub fn apply_at(&self, t: f64, psi: &Mps) -> Mps {
        let a = self.v.apply_to(psi);
        let b = self.generator_at(t).apply_to(&a);
        self.v_dag.apply_to(&b)
    }
}

/// One segment of a [`WidthCursor`] pipeline: the collapsed operator for the
/// flow span `[t0, t1)`.
pub struct Segment {
    pub t0: f64,
    pub t1: f64,
    pub op: Mpo,
}

impl Segment {
    pub fn width(&self) -> f64 {
        self.t1 - self.t0
    }
}

/// An adaptive-granularity pipeline over an operator flow: consecutive
/// segments covering `[0, T)`, each held at its own temporal resolution.
/// The composition of all segments is invariant under [`WidthCursor::refine`]
/// and [`WidthCursor::coarsen`] (up to the truncation policy).
pub struct WidthCursor<'a> {
    pub flow: &'a FourierFlow,
    pub segments: Vec<Segment>,
}

impl<'a> WidthCursor<'a> {
    /// Start with a single segment covering `[0, t_total)`.
    pub fn whole(flow: &'a FourierFlow, t_total: f64) -> WidthCursor<'a> {
        WidthCursor {
            flow,
            segments: vec![Segment {
                t0: 0.0,
                t1: t_total,
                op: flow.span(0.0, t_total),
            }],
        }
    }

    /// Split segment `i` at its midpoint: both halves are imputed from the
    /// flow (the geometric interpolation consistent with the coarse block —
    /// their composition reproduces the parent by the group law).
    pub fn refine(&mut self, i: usize) {
        let (t0, t1) = (self.segments[i].t0, self.segments[i].t1);
        let tm = 0.5 * (t0 + t1);
        let first = Segment {
            t0,
            t1: tm,
            op: self.flow.span(t0, tm),
        };
        let second = Segment {
            t0: tm,
            t1,
            op: self.flow.span(tm, t1),
        };
        self.segments[i] = first;
        self.segments.insert(i + 1, second);
    }

    /// Merge segments `i` and `i+1` into one wider block by composition.
    pub fn coarsen(&mut self, i: usize) {
        assert!(i + 1 < self.segments.len());
        let second = self.segments.remove(i + 1);
        let first = &self.segments[i];
        let op = second.op.compose_after(&first.op, self.flow.trunc);
        self.segments[i] = Segment {
            t0: first.t0,
            t1: second.t1,
            op,
        };
    }

    /// Compose the whole pipeline into one operator (later segments applied
    /// after earlier ones).
    pub fn total(&self) -> Mpo {
        let mut acc = Mpo::identity(&self.flow.profile, self.flow.trunc);
        for seg in &self.segments {
            acc = seg.op.compose_after(&acc, self.flow.trunc);
        }
        acc
    }

    /// Run a state through the pipeline segment by segment, returning every
    /// intermediate state (the trajectory of the computation at the
    /// pipeline's current resolution). `out[k]` is the state after the
    /// first `k` segments; `out[0]` is the input.
    pub fn trajectory(&self, psi: &Mps) -> Vec<Mps> {
        let mut out = vec![psi.clone()];
        for seg in &self.segments {
            let next = seg.op.apply_to(out.last().unwrap());
            out.push(next);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::c64::C64;
    use crate::dense::{total_dim, DenseState};
    use crate::zigzag;

    const SPEC: TruncSpec = TruncSpec {
        max_rank: 96,
        cutoff: 1e-12,
    };

    #[test]
    fn flow_is_a_one_parameter_group() {
        let profile = zigzag::diamond(2, 4); // N = 144
        let flow = FourierFlow::adder(&profile, 5.0, SPEC);
        // U^0 = identity.
        let id = Mpo::identity(&profile, SPEC);
        let f0 = flow.at(0.0).hs_fidelity(&id);
        assert!((f0 - 1.0).abs() < 1e-8, "U^0 vs I: {}", f0);
        // Group law with fractional pieces.
        let a = flow.at(0.3);
        let b = flow.at(0.7);
        let ab = b.compose_after(&a, SPEC);
        let whole = flow.at(1.0);
        let f = ab.hs_fidelity(&whole);
        assert!((f - 1.0).abs() < 1e-8, "U^.3∘U^.7 vs U^1: {}", f);
        // U^1 is the exact adder.
        let exact = radix::adder_mpo_exact(&profile, 5, SPEC);
        let f1 = whole.hs_fidelity(&exact);
        assert!((f1 - 1.0).abs() < 1e-8, "U^1 vs exact adder: {}", f1);
    }

    #[test]
    fn integer_snapshots_are_narrow_fractional_are_wide() {
        // Operator width sees the integers: at whole shifts the flow
        // collapses to a permutation (tiny chi); at half-integer shifts it
        // is a delocalized kernel (much larger chi).
        let profile = zigzag::diamond(2, 4);
        let flow = FourierFlow::adder(&profile, 1.0, SPEC);
        let chi_int = flow.at(1.0).max_bond_dim();
        let chi_frac = flow.at(0.5).max_bond_dim();
        assert!(chi_int <= 4, "integer snapshot chi {}", chi_int);
        assert!(
            chi_frac >= 2 * chi_int,
            "fractional chi {} not ≫ integer chi {}",
            chi_frac,
            chi_int
        );
    }

    #[test]
    fn fractional_halves_compose_to_the_full_step_on_states() {
        let profile = zigzag::diamond(2, 4);
        let n_total = total_dim(&profile);
        let flow = FourierFlow::adder(&profile, 7.0, SPEC);
        let half = flow.at(0.5);
        for x in [0usize, 40, 143] {
            let mut basis = DenseState::zero_state(&profile);
            basis.amps[0] = C64::ZERO;
            basis.amps[x] = C64::ONE;
            let psi = Mps::from_dense(&basis, SPEC);
            let out = half.apply_to(&half.apply_to(&psi)).to_dense();
            let target = (x + 7) % n_total;
            assert!(
                (out.amps[target].abs() - 1.0).abs() < 1e-7,
                "x={}: weight {}",
                x,
                out.amps[target].abs()
            );
        }
    }

    #[test]
    fn cursor_refine_and_coarsen_preserve_the_total() {
        let profile = zigzag::diamond(2, 4);
        let flow = FourierFlow::adder(&profile, 3.0, SPEC);
        let reference = flow.at(1.0);

        let mut cursor = WidthCursor::whole(&flow, 1.0);
        // Dyadic zoom on the leading edge: [0,.5],[.5,.75],[.75,.875],[.875,1].
        cursor.refine(0);
        cursor.refine(1);
        cursor.refine(2);
        assert_eq!(cursor.segments.len(), 4);
        let widths: Vec<f64> = cursor.segments.iter().map(|s| s.width()).collect();
        assert!((widths.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        let f = cursor.total().hs_fidelity(&reference);
        assert!((f - 1.0).abs() < 1e-7, "refined total vs whole: {}", f);

        // Coarsen everything back into one block.
        while cursor.segments.len() > 1 {
            cursor.coarsen(0);
        }
        let f2 = cursor.total().hs_fidelity(&reference);
        assert!((f2 - 1.0).abs() < 1e-7, "coarsened total vs whole: {}", f2);
    }

    #[test]
    fn trajectory_walks_the_ring_continuously() {
        // Stepping |x⟩ through quarter-steps of a unit shift. The circular
        // first moment ⟨ω^x⟩ advances *exactly* linearly (Weyl covariance:
        // A^t† W A^t = ω^{t·c} W), staying at modulus 1, while lattice-basis
        // participation breathes: delocalized at fractional t, relocalized
        // at t = 1.
        let profile = zigzag::diamond(2, 4);
        let n_total = total_dim(&profile);
        let flow = FourierFlow::adder(&profile, 1.0, SPEC);
        let mut cursor = WidthCursor::whole(&flow, 1.0);
        cursor.refine(0);
        cursor.refine(0);
        cursor.refine(2); // widths [1/4, 1/4, 1/4, 1/4]

        let x0 = 70usize;
        let mut basis = DenseState::zero_state(&profile);
        basis.amps[0] = C64::ZERO;
        basis.amps[x0] = C64::ONE;
        let traj = cursor.trajectory(&Mps::from_dense(&basis, SPEC));
        assert_eq!(traj.len(), 5);

        let circular = |s: &Mps| -> C64 {
            let d = s.to_dense();
            let mut acc = C64::ZERO;
            for (x, a) in d.amps.iter().enumerate() {
                acc += C64::cis(2.0 * std::f64::consts::PI * x as f64 / n_total as f64)
                    .scale(a.abs2());
            }
            acc
        };
        let participation = |s: &Mps| -> f64 {
            let d = s.to_dense();
            d.amps.iter().map(|a| a.abs2() * a.abs2()).sum()
        };

        // Exact first-moment law of the fractional shift on a discrete ring:
        // ⟨ω^x⟩ = ω^{x0+tc} · ((N-1) + e^{-2πi·tc}) / N — linear drift plus
        // a "lattice wobble" from the wrap term, vanishing at integer tc.
        let tau = 2.0 * std::f64::consts::PI;
        for (k, state) in traj.iter().enumerate() {
            let tc = k as f64 * 0.25;
            let nf = n_total as f64;
            let pred = C64::cis(tau * (x0 as f64 + tc) / nf)
                * (C64::real(nf - 1.0) + C64::cis(-tau * tc)).scale(1.0 / nf);
            let got = circular(state);
            assert!(
                (got - pred).abs() < 1e-5,
                "step {}: moment {:?} predicted {:?}",
                k,
                got,
                pred
            );
        }
        // Participation breathes: fully localized at the integer endpoints,
        // spread out in between.
        assert!((participation(&traj[0]) - 1.0).abs() < 1e-9);
        assert!((participation(&traj[4]) - 1.0).abs() < 1e-6);
        for state in &traj[1..4] {
            assert!(
                participation(state) < 0.7,
                "fractional step not delocalized"
            );
        }
        // Endpoint fully localized on x0 + 1.
        let final_dense = traj[4].to_dense();
        assert!((final_dense.amps[(x0 + 1) % n_total].abs() - 1.0).abs() < 1e-7);
    }
}
