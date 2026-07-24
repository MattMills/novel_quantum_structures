//! **Crossing dimension-wave strands**: X-shaped and networked geometries
//! built from two (or more) waves that interact across scales.
//!
//! Everything else in this crate lives on a single chain. This module
//! studies *two* chains — two strands `A` and `B`, each a dimension wave —
//! **crossed** so that they interact at sites of matching dimension. The
//! canonical picture is two peaks `[1,2,3,4,5,4,3,2,1]` laid antiparallel,
//! forming an X whose strands meet at the center:
//!
//! ```text
//!   A:  1  2  3  4  5  4  3  2  1
//!        ╲  ╲  ╲  ╲ ╳ ╱  ╱  ╱  ╱      antiparallel pairing A[i] ~ B[n−1−i]
//!   B:  1  2  3  4  5  4  3  2  1      (every pair has equal dimension)
//! ```
//!
//! The two strands share one MPS: global sites `0..n` are strand `A`,
//! `n..2n` are strand `B`, so the **A|B bipartition is the single bond
//! `n−1`** and its entropy is exactly the inter-strand entanglement
//! ([`Crossing::ab_entropy_bits`]). Crossing couplers are then ordinary
//! long-range two-site gates (the [`crate::mps`] engine threads them across
//! the boundary), and the *interaction level* — which dimension band `5`,
//! `4`, `3`, `2`, or `1` carries the coupling — is a knob with a precise,
//! measured consequence:
//!
//! * **The crossing entanglement budget.** Bell-coupling the pairs at a
//!   level set `L` injects exactly `Σ log2 d` bits across A|B (Schmidt
//!   rank `Π d`), summed over the *active pairs* — so multiplicity matters.
//! * **Multiplicity beats dimension.** On a peak the waist (`d = hi`) is
//!   *unique* but every lower level is *paired*, so the most entanglement
//!   enters at the **shoulder** `d = hi−1` (two pairs), not the peak — for
//!   `[1..5..1]`, `2·log2 4 = 4` bits at the shoulder against `log2 5 =
//!   2.32` at the peak. A valley flips this: its wide ends are paired.
//! * **The dimension-1 pinch is a decoupled crossing.** `cshift(1,1)` and
//!   `xswap(1,·)` are the identity, so coupling "at level 1" does nothing;
//!   a `d = 1` site carries no entanglement of its own yet still routes
//!   bond correlation past it — a spacer, not a wall.
//! * **Cycles add.** Couplings at disjoint levels commute and their A|B
//!   entropies sum, so "interact at `d` in cycle `k`" accumulates linearly
//!   up to the total budget `Σ_i log2 P_i`.
//!
//! See `examples/crossing_vees.rs` and THEORY.md §11.

use crate::circuit::Circuit;
use crate::gates;
use crate::mps::Mps;

/// How the two strands' sites are ordered along the shared MPS.
///
/// This is the crossing's central performance knob, and it is *purely a
/// representation choice* — the physical state is identical either way.
///
/// * [`Layout::Block`] lays strand `A` then strand `B`, so the A|B
///   bipartition is the single bond `n−1` and its entropy is read
///   directly ([`Crossing::ab_entropy_bits`]). But every crossing pair is
///   then long-range, and a fully-coupled crossing forces *all* the
///   inter-strand entanglement through one contiguous cut — bond
///   dimension up to `Π d` (2880 for twin `[1..5..1]`). Clean to measure,
///   costly to compute.
/// * [`Layout::Interleaved`] places each crossing pair on adjacent sites,
///   so every coupler is nearest-neighbour and a Bell crossing is a
///   near-product of *local* pairs: every contiguous cut sees at most one
///   pair, so bond dimension stays `≤ hi` however many levels are coupled.
///   The A|B entanglement is then a non-contiguous (comb) partition rather
///   than one bond, but the *computation* is orders of magnitude cheaper.
///
/// The rule of thumb: **interleave to compute, block to read off the
/// inter-strand entanglement** — and the gap between their bond
/// dimensions is itself a measurement (see `examples/crossing_vees.rs`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Layout {
    Block,
    Interleaved,
}

/// Two dimension-wave strands sharing one MPS, crossed pairwise.
pub struct Crossing {
    /// Strand `A` profile.
    pub a: Vec<usize>,
    /// Strand `B` profile.
    pub b: Vec<usize>,
    /// Antiparallel pairing `A[i] ~ B[n−1−i]` (the X); else parallel
    /// `A[i] ~ B[i]` (a ladder).
    pub antiparallel: bool,
    /// Site ordering along the shared MPS.
    pub layout: Layout,
}

impl Crossing {
    pub fn new(a: &[usize], b: &[usize], antiparallel: bool) -> Crossing {
        assert_eq!(
            a.len(),
            b.len(),
            "strands must have equal length to cross pairwise"
        );
        Crossing {
            a: a.to_vec(),
            b: b.to_vec(),
            antiparallel,
            layout: Layout::Block,
        }
    }

    /// Two identical strands (the usual X), in [`Layout::Block`].
    pub fn twin(strand: &[usize], antiparallel: bool) -> Crossing {
        Crossing::new(strand, strand, antiparallel)
    }

    /// Same crossing, re-laid in `layout` (consumes and returns `self`).
    pub fn with_layout(mut self, layout: Layout) -> Crossing {
        self.layout = layout;
        self
    }

    /// Combined chain dimensions, in this crossing's layout order.
    pub fn dims(&self) -> Vec<usize> {
        match self.layout {
            Layout::Block => {
                let mut d = self.a.clone();
                d.extend_from_slice(&self.b);
                d
            }
            Layout::Interleaved => {
                let n = self.n();
                let mut d = vec![0usize; 2 * n];
                for i in 0..n {
                    d[self.ga(i)] = self.a[i];
                    d[self.gb(self.partner(i))] = self.b[self.partner(i)];
                }
                d
            }
        }
    }

    /// Sites per strand.
    pub fn n(&self) -> usize {
        self.a.len()
    }

    /// Global (MPS) index of strand-A site `i`.
    pub fn ga(&self, i: usize) -> usize {
        match self.layout {
            Layout::Block => i,
            Layout::Interleaved => 2 * i,
        }
    }

    /// Global (MPS) index of strand-B site `j`.
    pub fn gb(&self, j: usize) -> usize {
        match self.layout {
            Layout::Block => self.n() + j,
            // B[j] is the partner of A[i] with partner(i) = j; that pair
            // sits at MPS positions (2i, 2i+1).
            Layout::Interleaved => {
                let i = if self.antiparallel { self.n() - 1 - j } else { j };
                2 * i + 1
            }
        }
    }

    /// The bond whose entropy is the A|B entanglement — only meaningful in
    /// [`Layout::Block`], where A and B are contiguous.
    pub fn ab_cut(&self) -> usize {
        assert_eq!(
            self.layout,
            Layout::Block,
            "the A|B cut is a single bond only in Block layout"
        );
        self.a.len() - 1
    }

    fn partner(&self, i: usize) -> usize {
        if self.antiparallel {
            self.n() - 1 - i
        } else {
            i
        }
    }

    /// Crossing pairs `(global A index, global B index, shared dim)` — only
    /// where the paired sites have equal dimension.
    pub fn pairs(&self) -> Vec<(usize, usize, usize)> {
        (0..self.n())
            .filter_map(|i| {
                let j = self.partner(i);
                if self.a[i] == self.b[j] {
                    Some((self.ga(i), self.gb(j), self.a[i]))
                } else {
                    None
                }
            })
            .collect()
    }

    /// Crossing pairs whose shared dimension is in `levels`.
    pub fn pairs_at(&self, levels: &[usize]) -> Vec<(usize, usize, usize)> {
        self.pairs()
            .into_iter()
            .filter(|&(_, _, d)| levels.contains(&d))
            .collect()
    }

    /// Generalized-Bell coupling on every crossing pair at `levels`:
    /// `fourier(d)` on the A partner, `cshift(d, d)` onto the B partner —
    /// the state `Σ_k |k,k⟩/√d` per pair, maximal A|B entanglement.
    pub fn bell(&self, levels: &[usize]) -> Circuit {
        let mut c = Circuit::new(self.dims());
        for (ga, gb, d) in self.pairs_at(levels) {
            c.one(ga, gates::fourier(d));
            c.two(ga, gb, gates::cshift(d, d));
        }
        c
    }

    /// Controlled-phase coupling on every crossing pair at `levels`
    /// (entangles only *pre-spread* strands — compose after [`spread`]).
    ///
    /// [`spread`]: Crossing::spread
    pub fn phase(&self, levels: &[usize]) -> Circuit {
        let mut c = Circuit::new(self.dims());
        for (ga, gb, d) in self.pairs_at(levels) {
            c.two(ga, gb, gates::cphase(d, d));
        }
        c
    }

    /// Fourier-spread both partners of every crossing pair at `levels`.
    pub fn spread(&self, levels: &[usize]) -> Circuit {
        let mut c = Circuit::new(self.dims());
        for (ga, gb, d) in self.pairs_at(levels) {
            c.one(ga, gates::fourier(d));
            c.one(gb, gates::fourier(d));
        }
        c
    }

    /// Predicted A|B entanglement (bits) of Bell-coupling `levels`:
    /// `Σ log2 d` over the active pairs.
    pub fn bell_entropy_bits(&self, levels: &[usize]) -> f64 {
        self.pairs_at(levels)
            .iter()
            .map(|&(_, _, d)| (d as f64).log2())
            .sum()
    }

    /// Predicted A|B Schmidt rank of Bell-coupling `levels`: `Π d`.
    pub fn bell_rank(&self, levels: &[usize]) -> usize {
        self.pairs_at(levels).iter().map(|&(_, _, d)| d).product()
    }

    /// Measured A|B entanglement (bits) of a state on this crossing —
    /// [`Layout::Block`] only, where it is the single bond [`ab_cut`].
    /// Uses the targeted single-bond entropy (`O(χ³)`), not a full sweep.
    ///
    /// [`ab_cut`]: Crossing::ab_cut
    pub fn ab_entropy_bits(&self, psi: &mut Mps) -> f64 {
        psi.bond_entropy_bits_at(self.ab_cut())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dense::DenseState;
    use crate::mat::Rng;
    use crate::mps::Mps;
    use crate::zigzag;
    use crate::TruncSpec;

    fn run_both(dims: &[usize], c: &Circuit) -> (Mps, DenseState) {
        let mut m = Mps::zero_state(dims, TruncSpec::exact());
        c.run_mps(&mut m);
        let mut d = DenseState::zero_state(dims);
        c.run_dense(&mut d);
        (m, d)
    }

    #[test]
    fn bell_crossing_matches_dense() {
        // Full X coupling on twin peaks, cross-validated against dense.
        let strand = zigzag::diamond(1, 4); // [1,2,3,4,3,2,1], dim 144
        let x = Crossing::twin(&strand, true);
        let dims = x.dims();
        let (m, d) = run_both(&dims, &x.bell(&[1, 2, 3, 4]));
        assert!((m.norm() - 1.0).abs() < 1e-12);
        assert!((m.to_dense().fidelity(&d) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn ab_entropy_matches_the_budget() {
        // Each single level's measured A|B entropy = Σ log2 d over its
        // active pairs, and rank = Π d, verified on the peak [1..5..1].
        let strand = zigzag::diamond(1, 5);
        let x = Crossing::twin(&strand, true);
        let dims = x.dims();
        for level in [5usize, 4, 3, 2, 1] {
            let mut m = Mps::zero_state(&dims, TruncSpec::exact());
            x.bell(&[level]).run_mps(&mut m);
            let measured = x.ab_entropy_bits(&mut m);
            let predicted = x.bell_entropy_bits(&[level]);
            assert!(
                (measured - predicted).abs() < 1e-9,
                "level {}: measured {} predicted {}",
                level,
                measured,
                predicted
            );
            // Rank across the A|B cut.
            m.move_center_to(x.ab_cut());
            assert_eq!(m.bond_dims()[x.ab_cut()], x.bell_rank(&[level]));
        }
    }

    #[test]
    fn shoulder_beats_the_peak_and_pinch_is_trivial() {
        // Multiplicity vs dimension: the shoulder (d=4, two pairs) injects
        // more A|B entanglement than the unique peak (d=5); the d=1 pinch
        // injects none.
        let x = Crossing::twin(&zigzag::diamond(1, 5), true);
        let peak = x.bell_entropy_bits(&[5]);
        let shoulder = x.bell_entropy_bits(&[4]);
        assert!((peak - 5f64.log2()).abs() < 1e-12);
        assert!((shoulder - 4.0).abs() < 1e-12);
        assert!(shoulder > peak, "shoulder {} !> peak {}", shoulder, peak);
        assert_eq!(x.bell_entropy_bits(&[1]), 0.0);
        assert_eq!(x.bell_rank(&[1]), 1);
        assert_eq!(x.pairs_at(&[5]).len(), 1); // unique peak
        assert_eq!(x.pairs_at(&[4]).len(), 2); // paired shoulder
    }

    #[test]
    fn valley_flips_the_multiplicity() {
        // A valley's wide ends are paired and its pinch is unique — the
        // exact reflection of the diamond.
        let x = Crossing::twin(&zigzag::valley(1, 5), true);
        assert_eq!(x.pairs_at(&[5]).len(), 2); // paired wide ends
        assert_eq!(x.pairs_at(&[1]).len(), 1); // unique pinch
        assert!((x.bell_entropy_bits(&[5]) - 2.0 * 5f64.log2()).abs() < 1e-12);
        assert_eq!(x.bell_entropy_bits(&[1]), 0.0);
    }

    #[test]
    fn cycles_add_disjoint_levels() {
        // Coupling at level 4 then level 3 in a second cycle: the A|B
        // entropies add (disjoint sites commute), and the composite equals
        // one round at both levels — verified exactly against dense.
        let strand = zigzag::diamond(1, 4);
        let x = Crossing::twin(&strand, true);
        let dims = x.dims();
        let mut two_cycle = Mps::zero_state(&dims, TruncSpec::exact());
        x.bell(&[4]).run_mps(&mut two_cycle);
        x.bell(&[3]).run_mps(&mut two_cycle);
        let one_round = {
            let mut m = Mps::zero_state(&dims, TruncSpec::exact());
            x.bell(&[3, 4]).run_mps(&mut m);
            m
        };
        assert!((two_cycle.fidelity(&one_round) - 1.0).abs() < 1e-10);
        let measured = x.ab_entropy_bits(&mut two_cycle);
        let predicted = x.bell_entropy_bits(&[4]) + x.bell_entropy_bits(&[3]);
        assert!(
            (measured - predicted).abs() < 1e-9,
            "measured {} predicted {}",
            measured,
            predicted
        );
    }

    #[test]
    fn parallel_and_antiparallel_share_the_budget() {
        // The X (antiparallel) and the ladder (parallel) couple the same
        // dimension multiset, so their A|B entropy budgets match — the
        // difference is which physical sites pair, not how much entangles.
        let strand = zigzag::diamond(1, 4);
        for level in [4usize, 3, 2] {
            let xa = Crossing::twin(&strand, true);
            let xp = Crossing::twin(&strand, false);
            assert!(
                (xa.bell_entropy_bits(&[level]) - xp.bell_entropy_bits(&[level])).abs() < 1e-12
            );
        }
    }

    #[test]
    fn interleaved_layout_keeps_the_full_crossing_cheap() {
        // The pivotal performance fact: in interleaved layout a full Bell
        // crossing is a near-product of local pairs, so every contiguous
        // bond stays ≤ hi — no matter how many levels are coupled — while
        // block layout drives the same state to Π d = 2880.
        let strand = zigzag::diamond(1, 5);
        let hi = 5;
        let xi = Crossing::twin(&strand, true).with_layout(Layout::Interleaved);
        let mut m = Mps::zero_state(&xi.dims(), TruncSpec::exact());
        xi.bell(&[1, 2, 3, 4, 5]).run_mps(&mut m);
        assert!(
            m.max_bond_dim() <= hi,
            "interleaved full crossing max χ {} > hi {}",
            m.max_bond_dim(),
            hi
        );
        // Block layout, same physics, is forced through one contiguous cut.
        let xb = Crossing::twin(&zigzag::diamond(1, 4), true); // smaller: Π d = 144
        let mut mb = Mps::zero_state(&xb.dims(), TruncSpec::exact());
        xb.bell(&[1, 2, 3, 4]).run_mps(&mut mb);
        assert_eq!(mb.max_bond_dim(), xb.bell_rank(&[1, 2, 3, 4])); // = 144
        assert!(mb.max_bond_dim() > hi);
    }

    #[test]
    fn interleaved_couplings_are_nearest_neighbour() {
        // Every crossing pair is adjacent in interleaved layout — the
        // reason it is cheap (no long-range threading).
        let xi = Crossing::twin(&zigzag::diamond(1, 5), true).with_layout(Layout::Interleaved);
        for (ga, gb, _) in xi.pairs() {
            assert_eq!(ga.abs_diff(gb), 1, "pair {}~{} not adjacent", ga, gb);
        }
    }

    #[test]
    fn interleaved_matches_dense() {
        // Correctness of the interleaved build against dense in the same
        // site order.
        let xi = Crossing::twin(&zigzag::diamond(1, 4), true).with_layout(Layout::Interleaved);
        let dims = xi.dims();
        let circuit = xi.bell(&[1, 2, 3, 4]);
        let mut m = Mps::zero_state(&dims, TruncSpec::exact());
        circuit.run_mps(&mut m);
        let mut d = DenseState::zero_state(&dims);
        circuit.run_dense(&mut d);
        assert!((m.to_dense().fidelity(&d) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn coarse_graining_preserves_a_crossing() {
        // RG covariance: merging a crossed strand toward its waist (an
        // exact coarse-graining step) carries the crossing coupling up to
        // the coarse block and splits back, state preserved exactly.
        let strand = zigzag::diamond(1, 4); // [1,2,3,4,3,2,1]
        let x = Crossing::twin(&strand, false); // ladder, block layout
        let mut crossed = Mps::zero_state(&x.dims(), TruncSpec::exact());
        x.bell(&[2, 3, 4]).run_mps(&mut crossed);

        let mut coarse = crossed.clone();
        coarse.merge_sites(2); // fuse (2,3): d = 12
        coarse.merge_sites(2); // fuse (·,3): d = 36 — the coarse waist block
        assert_eq!(coarse.dims[2], 36);
        coarse.split_site(2, 12, 3);
        coarse.split_site(2, 3, 4);
        assert!(
            (coarse.fidelity(&crossed) - 1.0).abs() < 1e-10,
            "RG step broke the crossing: {}",
            coarse.fidelity(&crossed)
        );
    }

    #[test]
    fn crossing_has_a_scale() {
        // On a wave, coupling at the waist (coarse) is few crossings of fat
        // modes; toward the valley (fine), many crossings of thin modes.
        let w = zigzag::wave(1, 5, 2);
        let x = Crossing::twin(&w, true);
        // Two periods → two shared waists → two level-5 crossings.
        assert_eq!(x.pairs_at(&[5]).len(), 2);
        // The fine valley band (d=2) has more crossings than the coarse waist.
        assert!(x.pairs_at(&[2]).len() > x.pairs_at(&[5]).len());
        // Coarse crossings carry more entropy per crossing (fat modes).
        let per_waist = x.bell_entropy_bits(&[5]) / x.pairs_at(&[5]).len() as f64;
        let per_valley = x.bell_entropy_bits(&[2]) / x.pairs_at(&[2]).len() as f64;
        assert!(per_waist > per_valley, "{} !> {}", per_waist, per_valley);
    }

    #[test]
    fn structured_crossing_beats_scrambling() {
        // A structured single-level crossing stays at its budget rank; a
        // random cross coupling on the same pairs demands far more.
        let strand = zigzag::diamond(1, 4);
        let x = Crossing::twin(&strand, true);
        let dims = x.dims();
        let mut structured = Mps::zero_state(&dims, TruncSpec::exact());
        x.bell(&[3]).run_mps(&mut structured);
        structured.move_center_to(x.ab_cut());
        let struct_rank = structured.bond_dims()[x.ab_cut()];

        let mut rng = Rng::new(7);
        let mut scrambled = Mps::zero_state(&dims, TruncSpec::exact());
        // Spread, then random 2-site gates on the same crossing pairs.
        x.spread(&[3]).run_mps(&mut scrambled);
        for (ga, gb, d) in x.pairs_at(&[3]) {
            scrambled.apply2(ga, gb, &gates::random(d * d, &mut rng));
        }
        scrambled.move_center_to(x.ab_cut());
        let scrambled_rank = scrambled.bond_dims()[x.ab_cut()];
        assert_eq!(struct_rank, 9, "structured rank");
        assert!(
            scrambled_rank >= struct_rank,
            "scrambled {} vs structured {}",
            scrambled_rank,
            struct_rank
        );
    }
}
