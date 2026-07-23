//! Bulk sheaf towers: stepwise qudit growth as a multiresolution cache.
//!
//! A single state has many exact tensor-network representations. Fusing two
//! neighbouring sites into one coarser, higher-dimensional qudit
//! ([`crate::mps::Mps::merge_sites`]) is lossless — the interior bond is
//! absorbed into the enlarged local space — so repeatedly fusing builds a
//! **tower** of representations of *one invariant global state*, each layer
//! coarser (fewer sites, larger local dimension `d`) than the last. Growth
//! carries the profile well past `d = 10`: a `2,3,…,12,…,3,2` diamond fuses
//! to sites of dimension in the hundreds.
//!
//! Read as geometry, each layer is a **section of the bulk** at a coarser
//! resolution: the tower is a *cellular sheaf* whose edge stalks are the
//! Schmidt (bond) spaces, whose restriction map is canonicalize-and-slice
//! ([`restrict`]), and whose gluing is contraction ([`glue`]). The gluing
//! axiom holds exactly — local sections over a cover reassemble the global
//! section — which is what makes the whole tower "only for performance": no
//! layer is more or less true than another; they differ only in cost.
//!
//! And cost is what you choose a layer for. Coarsening shrinks site count
//! and absorbs interior bonds (cheaper for long-range operations and
//! bond-dominated states) while growing local dimension (dearer for
//! local operators, whose cost scales with `d²`). There is an optimal layer
//! per task, and [`BulkTower`] is the menu; `examples/bulk_sheaf_tower.rs`
//! measures the trade.

use crate::mps::{Mps, SiteTensor};
use crate::TruncSpec;

/// Fuse the sites of `m` in contiguous groups: `group_sizes` partitions the
/// chain (must sum to the site count), and each group becomes one site of
/// dimension equal to the product of its fine dimensions. Exact (merging
/// never truncates).
pub fn coarsen(m: &Mps, group_sizes: &[usize]) -> Mps {
    assert_eq!(
        group_sizes.iter().sum::<usize>(),
        m.num_sites(),
        "group sizes must partition the chain"
    );
    let mut w = m.clone();
    // After processing group `pos`, it occupies site index `pos` in the
    // shrinking chain, so the next group's fine sites begin at `pos + 1`.
    for (pos, &g) in group_sizes.iter().enumerate() {
        assert!(g >= 1, "group size must be positive");
        for _ in 1..g {
            w.merge_sites(pos);
        }
    }
    w
}

/// Inverse of [`coarsen`]: split each coarse site back into the finer
/// dimensions recorded in `groups` (`groups[j]` = the ordered fine
/// dimensions fused into coarse site `j`). Splits under `trunc`; pass
/// [`TruncSpec::exact`] to invert a fusion losslessly.
pub fn refine(m: &Mps, groups: &[Vec<usize>], trunc: TruncSpec) -> Mps {
    assert_eq!(groups.len(), m.num_sites(), "one group per coarse site");
    let mut w = m.clone();
    w.trunc = trunc;
    let mut pos = 0;
    for grp in groups {
        let prod: usize = grp.iter().product();
        assert_eq!(
            prod, w.dims[pos],
            "group product must equal the coarse dimension"
        );
        for j in 0..grp.len().saturating_sub(1) {
            let left = grp[j];
            let rest: usize = grp[j + 1..].iter().product();
            w.split_site(pos, left, rest);
            pos += 1;
        }
        pos += 1;
    }
    w
}

/// Partition a dimension list into contiguous groups of at most `block`
/// sites, returning the group sizes and the grouped dimension values.
pub fn uniform_partition(dims: &[usize], block: usize) -> (Vec<usize>, Vec<Vec<usize>>) {
    assert!(block >= 1);
    let mut sizes = Vec::new();
    let mut groups = Vec::new();
    let mut i = 0;
    while i < dims.len() {
        let end = (i + block).min(dims.len());
        sizes.push(end - i);
        groups.push(dims[i..end].to_vec());
        i = end;
    }
    (sizes, groups)
}

/// Decode a fused Choi-carrier physical index into its `(out, in)` halves.
/// Fusing carrier sites interleaves the doubled legs as
/// `(out₀,in₀,out₁,in₁,…)`; the coarse operator needs them regrouped as
/// `(out₀…out_{g-1}, in₀…in_{g-1})`, i.e. `out·D + in` with `D = ∏ dᵢ`.
fn deinterleave_choi(mut idx: usize, group: &[usize]) -> (usize, usize) {
    let g = group.len();
    let mut radices = Vec::with_capacity(2 * g);
    for &d in group {
        radices.push(d);
        radices.push(d);
    }
    let mut digits = vec![0usize; 2 * g];
    for pos in (0..2 * g).rev() {
        digits[pos] = idx % radices[pos];
        idx /= radices[pos];
    }
    let mut out = 0;
    let mut inn = 0;
    for (k, &d) in group.iter().enumerate() {
        out = out * d + digits[2 * k];
        inn = inn * d + digits[2 * k + 1];
    }
    (out, inn)
}

/// Coarsen an operator the same way as a state: fuse the Choi carrier's
/// sites in contiguous groups, so a fused site's doubled dimension becomes
/// `(∏ dᵢ)²`. Applying a coarsened operator to the identically coarsened
/// state reproduces the fine result — the whole computation lives at a
/// coarser resolution. Fewer sites carry the bond-dimension cost, which is
/// the performance lever when entanglement (not local dimension) dominates.
pub fn coarsen_mpo(w: &crate::mpo::Mpo, group_sizes: &[usize]) -> crate::mpo::Mpo {
    assert_eq!(
        group_sizes.iter().sum::<usize>(),
        w.dims.len(),
        "group sizes must partition the operator's chain"
    );
    let fused = coarsen(&w.carrier, group_sizes);
    // Regroup each fused site's interleaved doubled legs into out/in halves.
    let mut groups: Vec<Vec<usize>> = Vec::with_capacity(group_sizes.len());
    let mut pos = 0;
    for &g in group_sizes {
        groups.push(w.dims[pos..pos + g].to_vec());
        pos += g;
    }
    let mut tensors = Vec::with_capacity(groups.len());
    for (t, group) in fused.tensors.iter().zip(groups.iter()) {
        let dbig: usize = group.iter().product();
        let p = dbig * dbig;
        debug_assert_eq!(p, t.d);
        let mut nt = SiteTensor::zeros(t.dl, p, t.dr);
        for idx in 0..p {
            let (out, inn) = deinterleave_choi(idx, group);
            let q = out * dbig + inn;
            for l in 0..t.dl {
                for r in 0..t.dr {
                    nt.set(l, q, r, t.at(l, idx, r));
                }
            }
        }
        tensors.push(nt);
    }
    let dims: Vec<usize> = groups.iter().map(|g| g.iter().product()).collect();
    let carrier = Mps {
        dims: dims.iter().map(|&d| d * d).collect(),
        tensors,
        center: 0,
        trunc: w.carrier.trunc,
        discarded_weight: 0.0,
    };
    crate::mpo::Mpo { dims, carrier }
}

/// Per-layer cost summary.
#[derive(Clone, Debug)]
pub struct LayerStats {
    pub num_sites: usize,
    pub dims: Vec<usize>,
    pub max_local_dim: usize,
    pub max_bond: usize,
    /// Total stored complex parameters (state size).
    pub state_params: usize,
}

/// A tower of exact representations of one state at increasing coarseness.
pub struct BulkTower {
    /// `layers[0]` is the finest; each later layer fuses the previous one.
    pub layers: Vec<Mps>,
    /// `groups[k]` refines `layers[k+1]` back to `layers[k]`
    /// (per coarse site, the fine dimensions of `layers[k]` it fused).
    pub groups: Vec<Vec<Vec<usize>>>,
}

impl BulkTower {
    /// Build a tower from `fine` by fusing with a uniform block size at each
    /// layer transition (`schedule[k]` is the block size applied to layer
    /// `k` to make layer `k+1`).
    pub fn build(fine: &Mps, schedule: &[usize]) -> BulkTower {
        let mut layers = vec![fine.clone()];
        let mut groups = Vec::new();
        for &block in schedule {
            let cur = layers.last().unwrap();
            let (sizes, grp) = uniform_partition(&cur.dims, block);
            let coarse = coarsen(cur, &sizes);
            groups.push(grp);
            layers.push(coarse);
        }
        BulkTower { layers, groups }
    }

    pub fn num_layers(&self) -> usize {
        self.layers.len()
    }

    pub fn stats(&self, k: usize) -> LayerStats {
        let m = &self.layers[k];
        LayerStats {
            num_sites: m.num_sites(),
            dims: m.dims.clone(),
            max_local_dim: m.dims.iter().cloned().max().unwrap_or(1),
            max_bond: m.max_bond_dim(),
            state_params: m.param_count(),
        }
    }

    /// Refine `layers[k]` all the way back to the finest layout and return
    /// its fidelity with the stored fine state — 1 confirms the layer is the
    /// same global section (the "only for performance" invariant).
    pub fn fidelity_to_fine(&self, k: usize, trunc: TruncSpec) -> f64 {
        let mut w = self.layers[k].clone();
        for level in (0..k).rev() {
            w = refine(&w, &self.groups[level], trunc);
        }
        w.fidelity(&self.layers[0])
    }

    /// Refine `layers[k]` down to the finest layout (for applying a
    /// fine-scale observable or comparing results across layers).
    pub fn to_fine(&self, k: usize, trunc: TruncSpec) -> Mps {
        let mut w = self.layers[k].clone();
        for level in (0..k).rev() {
            w = refine(&w, &self.groups[level], trunc);
        }
        w
    }
}

// ---------------------------------------------------------------------------
// The cellular sheaf: sections, restriction, gluing.
// ---------------------------------------------------------------------------

/// A section of the bulk over a contiguous interval: an open MPS fragment
/// with dangling left/right boundary bonds (the edge stalks it must be
/// glued along).
#[derive(Clone)]
pub struct Section {
    pub dims: Vec<usize>,
    pub tensors: Vec<SiteTensor>,
}

impl Section {
    pub fn left_bond(&self) -> usize {
        self.tensors.first().map_or(1, |t| t.dl)
    }

    pub fn right_bond(&self) -> usize {
        self.tensors.last().map_or(1, |t| t.dr)
    }

    /// A closed section (both boundary bonds trivial) is a whole state.
    pub fn is_closed(&self) -> bool {
        self.left_bond() == 1 && self.right_bond() == 1
    }

    /// Interpret a closed section as an MPS.
    pub fn into_mps(self, trunc: TruncSpec) -> Mps {
        assert!(self.is_closed(), "only a closed section is a state");
        Mps {
            dims: self.dims,
            tensors: self.tensors,
            center: 0,
            trunc,
            discarded_weight: 0.0,
        }
    }
}

/// Restrict a state to the interval `[start, end)` — the sheaf's restriction
/// map. Slices out the site tensors in the state's current gauge; sections
/// taken from the same state share a gauge, so gluing a cover reproduces it
/// exactly.
pub fn restrict(m: &Mps, start: usize, end: usize) -> Section {
    assert!(start < end && end <= m.num_sites());
    Section {
        dims: m.dims[start..end].to_vec(),
        tensors: m.tensors[start..end].to_vec(),
    }
}

/// Glue two adjacent sections along their shared boundary bond — the sheaf's
/// gluing. Requires the boundary-bond dimensions to match.
pub fn glue(a: &Section, b: &Section) -> Section {
    assert_eq!(
        a.right_bond(),
        b.left_bond(),
        "sections do not agree on the shared boundary bond"
    );
    let mut dims = a.dims.clone();
    dims.extend_from_slice(&b.dims);
    let mut tensors = a.tensors.clone();
    tensors.extend_from_slice(&b.tensors);
    Section { dims, tensors }
}

/// Glue an ordered cover of sections into a single section (left to right).
pub fn glue_all(sections: &[Section]) -> Section {
    assert!(!sections.is_empty());
    let mut acc = sections[0].clone();
    for s in &sections[1..] {
        acc = glue(&acc, s);
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dense::DenseState;
    use crate::mat::Rng;
    use crate::zigzag;

    const SPEC: TruncSpec = TruncSpec {
        max_rank: 256,
        cutoff: 1e-13,
    };

    fn random_state(profile: &[usize], rng: &mut Rng) -> Mps {
        let mut s = DenseState::zero_state(profile);
        for a in &mut s.amps {
            *a = rng.c_gaussian();
        }
        s.normalize();
        Mps::from_dense(&s, SPEC)
    }

    #[test]
    fn coarsen_refine_round_trip_grows_past_ten() {
        let profile = zigzag::diamond(2, 5); // [2,3,4,5,4,3,2]
        let mut rng = Rng::new(1);
        let m = random_state(&profile, &mut rng);
        let (sizes, groups) = uniform_partition(&profile, 2);
        let coarse = coarsen(&m, &sizes);
        assert_eq!(coarse.dims, vec![6, 20, 12, 2]);
        assert!(
            coarse.dims.iter().cloned().max().unwrap() > 10,
            "grew past d=10"
        );
        let back = refine(&coarse, &groups, TruncSpec::exact());
        assert_eq!(back.dims, profile);
        assert!((back.fidelity(&m) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn every_tower_layer_is_the_same_section() {
        let profile = zigzag::diamond(2, 6); // max fine dim 6
        let mut rng = Rng::new(2);
        let m = random_state(&profile, &mut rng);
        let tower = BulkTower::build(&m, &[2, 2]);
        // Local dimension climbs well past ten as we grow.
        assert!(tower.stats(1).max_local_dim > 10);
        assert!(tower.stats(2).max_local_dim > tower.stats(1).max_local_dim);
        // Site count strictly falls.
        assert!(tower.stats(2).num_sites < tower.stats(1).num_sites);
        assert!(tower.stats(1).num_sites < tower.stats(0).num_sites);
        for k in 0..tower.num_layers() {
            let f = tower.fidelity_to_fine(k, TruncSpec::exact());
            assert!((f - 1.0).abs() < 1e-9, "layer {} fidelity {}", k, f);
        }
    }

    #[test]
    fn gluing_axiom_reassembles_the_global_section() {
        let profile = zigzag::diamond(2, 5);
        let mut rng = Rng::new(3);
        let m = random_state(&profile, &mut rng);
        let n = m.num_sites();
        // Whole restriction round-trips.
        let whole = restrict(&m, 0, n).into_mps(SPEC);
        assert!((whole.fidelity(&m) - 1.0).abs() < 1e-10);
        // A three-piece cover glues back to the state...
        let a = restrict(&m, 0, 2);
        let b = restrict(&m, 2, 5);
        let c = restrict(&m, 5, n);
        let glued = glue_all(&[a.clone(), b.clone(), c.clone()]).into_mps(SPEC);
        assert!((glued.fidelity(&m) - 1.0).abs() < 1e-10);
        // ...and gluing is associative (the sheaf axiom).
        let left_assoc = glue(&glue(&a, &b), &c).into_mps(SPEC);
        let right_assoc = glue(&a, &glue(&b, &c)).into_mps(SPEC);
        assert!((left_assoc.fidelity(&right_assoc) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn restriction_composes() {
        // Restricting an already-restricted section equals restricting the
        // global state directly (presheaf composition).
        let profile = zigzag::diamond(2, 5);
        let mut rng = Rng::new(4);
        let m = random_state(&profile, &mut rng);
        let outer = restrict(&m, 1, 6); // sites 1..6
                                        // Sub-restrict [2,4) of the *global* — same tensors as slicing outer.
        let direct = restrict(&m, 2, 4);
        // Outer holds sites 1..6, so its local indices 1..3 are global 2..4.
        let via_outer = Section {
            dims: outer.dims[1..3].to_vec(),
            tensors: outer.tensors[1..3].to_vec(),
        };
        assert_eq!(direct.dims, via_outer.dims);
        for (x, y) in direct.tensors.iter().zip(via_outer.tensors.iter()) {
            assert_eq!(x.data.len(), y.data.len());
            let diff = x
                .data
                .iter()
                .zip(y.data.iter())
                .map(|(a, b)| (*a - *b).abs())
                .fold(0.0, f64::max);
            assert!(diff < 1e-14, "restriction did not compose: {}", diff);
        }
    }

    #[test]
    fn coarsened_operator_reproduces_fine_application() {
        // A coarsened operator applied to the identically coarsened state,
        // refined back, must equal the fine application — exactly. Exercises
        // the Choi-leg de-interleaving across block sizes 2 and 3.
        use crate::cascade::Transducer;
        let profile = zigzag::diamond(2, 5);
        let n_total: u128 = profile.iter().map(|&d| d as u128).product();
        let x0 = 987u128;
        let mut digits = vec![0usize; profile.len()];
        let mut t = x0;
        for i in (0..profile.len()).rev() {
            digits[i] = (t % profile[i] as u128) as usize;
            t /= profile[i] as u128;
        }
        let fine = Mps::basis_state(&profile, &digits, SPEC);
        let w_fine = Transducer::adder(&profile, 250).to_mpo(SPEC);
        let out_fine = w_fine.apply_to(&fine);
        let target = (x0 + 250) % n_total;
        let mut tgt = vec![0usize; profile.len()];
        let mut tt = target;
        for i in (0..profile.len()).rev() {
            tgt[i] = (tt % profile[i] as u128) as usize;
            tt /= profile[i] as u128;
        }
        let expect = Mps::basis_state(&profile, &tgt, SPEC);

        for &block in &[2usize, 3] {
            let (sizes, groups) = uniform_partition(&profile, block);
            let state = coarsen(&fine, &sizes);
            let w = coarsen_mpo(&w_fine, &sizes);
            assert_eq!(w.dims, state.dims, "operator/state dims must match");
            let out = w.apply_to(&state);
            let refined = refine(&out, &groups, TruncSpec::exact());
            assert!(
                (refined.fidelity(&expect) - 1.0).abs() < 1e-9,
                "block {}: vs target {}",
                block,
                refined.fidelity(&expect)
            );
            assert!(
                (refined.fidelity(&out_fine) - 1.0).abs() < 1e-9,
                "block {}: vs fine",
                block
            );
        }
    }

    #[test]
    fn tower_preserves_state_under_a_shared_operation() {
        // Applying the same cascade at any layer, then refining, gives the
        // same result — the tower is a performance choice, not a physics one.
        use crate::cascade::Transducer;
        let profile = zigzag::diamond(2, 5);
        let n_total: u128 = profile.iter().map(|&d| d as u128).product();
        let x0 = 1234u128;
        let mut digits = vec![0usize; profile.len()];
        let mut t = x0;
        for i in (0..profile.len()).rev() {
            digits[i] = (t % profile[i] as u128) as usize;
            t /= profile[i] as u128;
        }
        let fine = Mps::basis_state(&profile, &digits, SPEC);
        let tower = BulkTower::build(&fine, &[2]);
        let target = (x0 + 100) % n_total;

        for k in 0..tower.num_layers() {
            let adder = Transducer::adder(&tower.layers[k].dims, 100).to_mpo(SPEC);
            let out = adder.apply_to(&tower.layers[k]);
            // Refine result to the fine layout and check the integer.
            let refined = if k == 0 {
                out
            } else {
                let mut w = out;
                for level in (0..k).rev() {
                    w = refine(&w, &tower.groups[level], TruncSpec::exact());
                }
                w
            };
            let mut tgt = vec![0usize; profile.len()];
            let mut tt = target;
            for i in (0..profile.len()).rev() {
                tgt[i] = (tt % profile[i] as u128) as usize;
                tt /= profile[i] as u128;
            }
            let f = refined.fidelity(&Mps::basis_state(&profile, &tgt, SPEC));
            assert!((f - 1.0).abs() < 1e-9, "layer {}: fidelity {}", k, f);
        }
    }
}
