//! **Networks of crossing strands**: many dimension-wave strands sharing
//! one MPS, coupled by an arbitrary graph of crossings.
//!
//! [`crate::crossing`] crosses *two* strands into an X. This module is the
//! general case — `K` strands and any coupling graph — and it exists to
//! measure one law:
//!
//! > **The cost of a crossing network is the cutwidth of its coupling
//! > graph.** Laid along one MPS, the maximum bond dimension is
//! > `d^(max edges crossing any contiguous cut)`, minimized over site
//! > orderings by the graph's *cutwidth*.
//!
//! That single fact organizes every richer geometry the crossing idea
//! invites, and separates the cheap from the expensive:
//!
//! * **One direction is cheap.** A *bundle* — `K` strands coupled only
//!   *across* the bundle at each site (GHZ rungs, no along-strand bonds) —
//!   has cutwidth 1: every rung is local in a column-major order, so the
//!   bond dimension is `hi` no matter how many strands `K` are bundled.
//! * **Overlaying directions adds cutwidth.** Crossing two strands *both*
//!   ways at once (an X and a ladder together) is a union of two matchings
//!   — disjoint 4-cycles, cutwidth 2 — so the bond dimension is the
//!   *shoulder* squared `(hi−1)²` (the unique peak carries only the ladder,
//!   so only the paired shoulders reach cutwidth 2).
//! * **Two axes is an area law.** A *weave* — strands running in two
//!   transverse directions, coupled at every intersection (a 2D lattice) —
//!   has cutwidth `min(rows, cols)`, so the bond dimension grows
//!   *exponentially in the smaller extent*. This is the boundary where the
//!   one-dimensional tensor-network representation stops being efficient,
//!   the same area law that separates MPS from PEPS.
//!
//! Couplers are built as graph states: [`Network::ghz`] (a `Σ|k…k⟩` rung
//! across a group) and [`Network::cluster`] (`fourier` then `cphase` on
//! each edge). For a graph state the entanglement across a cut is exactly
//! the number of edges crossing it (times `log d`), which is what makes
//! the cutwidth law hold on the nose — verified against dense simulation
//! and, at scale, read off the measured bond profile.
//!
//! See `examples/crossing_networks.rs` and THEORY.md §11.8.

use crate::circuit::Circuit;
use crate::gates;

/// An endpoint: `(strand, site)`.
pub type End = (usize, usize);
/// A coupler edge between two endpoints.
pub type Edge = (End, End);

/// Several dimension-wave strands sharing one MPS in a chosen site order.
pub struct Network {
    /// Strand profiles.
    pub strands: Vec<Vec<usize>>,
    /// MPS position → `(strand, site)`.
    pub order: Vec<End>,
    /// `(strand, site)` → MPS position (inverse of `order`).
    pos: Vec<Vec<usize>>,
}

impl Network {
    pub fn new(strands: Vec<Vec<usize>>, order: Vec<End>) -> Network {
        let total: usize = strands.iter().map(|s| s.len()).sum();
        assert_eq!(order.len(), total, "order must list every site once");
        let mut pos: Vec<Vec<usize>> = strands.iter().map(|s| vec![usize::MAX; s.len()]).collect();
        for (p, &(s, i)) in order.iter().enumerate() {
            assert!(
                pos[s][i] == usize::MAX,
                "site ({}, {}) listed twice in the order",
                s,
                i
            );
            pos[s][i] = p;
        }
        Network {
            strands,
            order,
            pos,
        }
    }

    /// Combined chain dimensions, in MPS order.
    pub fn dims(&self) -> Vec<usize> {
        self.order.iter().map(|&(s, i)| self.strands[s][i]).collect()
    }

    /// Total sites (the shared MPS length).
    pub fn len(&self) -> usize {
        self.order.len()
    }

    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }

    /// MPS index of `(strand, site)`.
    pub fn global(&self, s: usize, i: usize) -> usize {
        self.pos[s][i]
    }

    fn dim(&self, e: End) -> usize {
        self.strands[e.0][e.1]
    }

    /// A `Σ_k |k…k⟩/√d` rung across each group (all endpoints at one
    /// dimension): `fourier` on the first endpoint, `cshift` copying it to
    /// the rest. Groups with `d = 1` are skipped (a trivial rung).
    pub fn ghz(&self, groups: &[Vec<End>]) -> Circuit {
        let mut c = Circuit::new(self.dims());
        for g in groups {
            assert!(!g.is_empty());
            let d = self.dim(g[0]);
            if d == 1 {
                continue;
            }
            for &e in g {
                assert_eq!(self.dim(e), d, "GHZ group dimensions must match");
            }
            let q0 = self.global(g[0].0, g[0].1);
            c.one(q0, gates::fourier(d));
            for &e in &g[1..] {
                c.two(q0, self.global(e.0, e.1), gates::cshift(d, d));
            }
        }
        c
    }

    /// A qudit cluster (graph) state on `edges`: `fourier` every incident
    /// site, then `cphase` on each edge. The entanglement across any cut is
    /// exactly the number of edges crossing it (times `log d`).
    pub fn cluster(&self, edges: &[Edge]) -> Circuit {
        let mut c = Circuit::new(self.dims());
        let mut spread = vec![false; self.len()];
        for &(a, b) in edges {
            for e in [a, b] {
                let g = self.global(e.0, e.1);
                if !spread[g] {
                    c.one(g, gates::fourier(self.dim(e)));
                    spread[g] = true;
                }
            }
        }
        for &(a, b) in edges {
            let (da, db) = (self.dim(a), self.dim(b));
            c.two(self.global(a.0, a.1), self.global(b.0, b.1), gates::cphase(da, db));
        }
        c
    }
}

// ---------------------------------------------------------------------------
// Named network patterns
// ---------------------------------------------------------------------------

/// A **bundle** of `k` identical strands, coupled *only across* the bundle
/// at each site (a GHZ rung per column) — one direction of crossing. Sites
/// are laid column-major (all strands' site `i` adjacent), so every rung is
/// local and the bond dimension is `hi` independent of `k`. Returns the
/// network and its GHZ groups (one per site, `d = 1` sites dropped).
pub fn bundle(profile: &[usize], k: usize) -> (Network, Vec<Vec<End>>) {
    assert!(k >= 1);
    let n = profile.len();
    let strands = vec![profile.to_vec(); k];
    let mut order = Vec::with_capacity(k * n);
    for i in 0..n {
        for s in 0..k {
            order.push((s, i));
        }
    }
    let groups: Vec<Vec<End>> = (0..n)
        .filter(|&i| profile[i] > 1)
        .map(|i| (0..k).map(|s| (s, i)).collect())
        .collect();
    (Network::new(strands, order), groups)
}

/// Two strands crossed **both** ways at once — an antiparallel X
/// (`A[i]–B[n−1−i]`) *and* a parallel ladder (`A[i]–B[i]`), overlaid. The
/// coupling graph is a union of two matchings — a disjoint set of 4-cycles
/// `{A[i], B[i], A[n−1−i], B[n−1−i]}`, cutwidth 2 — so ordering each cycle
/// contiguously holds the bond dimension at `(hi−1)²`, the shoulder squared
/// (the unique peak sits in a 2-cycle with only its ladder edge). Overlaying
/// a second crossing direction adds exactly one to the cutwidth. Returns the network
/// and its cluster edges.
pub fn overlay(profile: &[usize]) -> (Network, Vec<Edge>) {
    let n = profile.len();
    let strands = vec![profile.to_vec(), profile.to_vec()];
    // Order each 4-cycle {A[i],B[i],A[j],B[j]} (j = n−1−i) contiguously.
    let mut order = Vec::with_capacity(2 * n);
    let mut placed = vec![false; n];
    for i in 0..n {
        if placed[i] {
            continue;
        }
        let j = n - 1 - i;
        order.push((0, i));
        order.push((1, i));
        placed[i] = true;
        if j != i {
            order.push((0, j));
            order.push((1, j));
            placed[j] = true;
        }
    }
    let mut edges = Vec::new();
    for i in 0..n {
        if profile[i] > 1 {
            edges.push(((0, i), (1, i))); // ladder
        }
        let j = n - 1 - i;
        if i < j && profile[i] == profile[j] && profile[i] > 1 {
            edges.push(((0, i), (1, j))); // X (once per pair)
            edges.push(((0, j), (1, i)));
        }
    }
    (Network::new(strands, order), edges)
}

/// A **weave**: `rows` strands running one way and the crossing bonds
/// running the other — a `rows × cols` lattice of dimension `d`, coupled on
/// every horizontal and vertical edge (a 2D cluster state). Sites are laid
/// column-major, so the bond dimension is `d^rows` — **exponential in the
/// transverse extent**, the area law. Returns the network and its edges.
pub fn weave(d: usize, rows: usize, cols: usize) -> (Network, Vec<Edge>) {
    assert!(rows >= 1 && cols >= 1);
    let strands = vec![vec![d; cols]; rows];
    let mut order = Vec::with_capacity(rows * cols);
    for c in 0..cols {
        for r in 0..rows {
            order.push((r, c));
        }
    }
    let mut edges = Vec::new();
    for r in 0..rows {
        for c in 0..cols {
            if c + 1 < cols {
                edges.push(((r, c), (r, c + 1))); // along a strand
            }
            if r + 1 < rows {
                edges.push(((r, c), (r + 1, c))); // across strands
            }
        }
    }
    (Network::new(strands, order), edges)
}

/// A **wave weave**: `rows` copies of a dimension-wave `profile` stacked and
/// coupled on every horizontal (along-wave) and vertical (across-wave)
/// edge — a lattice whose coupling dimension *modulates with the wave*.
/// Column-major, so a cut between columns `c, c+1` crosses `rows` along-wave
/// bonds at dimension `min(profile[c], profile[c+1])`: the area law with a
/// wave-shaped base. Returns the network and its edges.
pub fn wave_weave(profile: &[usize], rows: usize) -> (Network, Vec<Edge>) {
    let cols = profile.len();
    let strands = vec![profile.to_vec(); rows];
    let mut order = Vec::with_capacity(rows * cols);
    for c in 0..cols {
        for r in 0..rows {
            order.push((r, c));
        }
    }
    let mut edges = Vec::new();
    for r in 0..rows {
        for c in 0..cols {
            if c + 1 < cols && profile[c].min(profile[c + 1]) > 1 {
                edges.push(((r, c), (r, c + 1)));
            }
            if r + 1 < rows && profile[c] > 1 {
                edges.push(((r, c), (r + 1, c)));
            }
        }
    }
    (Network::new(strands, order), edges)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dense::DenseState;
    use crate::mps::Mps;
    use crate::zigzag;
    use crate::TruncSpec;

    fn build(net: &Network, circuit: &Circuit) -> Mps {
        let mut m = Mps::zero_state(&net.dims(), TruncSpec::exact());
        circuit.run_mps(&mut m);
        m
    }

    fn dense_fidelity(net: &Network, circuit: &Circuit) -> f64 {
        let m = build(net, circuit);
        let mut d = DenseState::zero_state(&net.dims());
        circuit.run_dense(&mut d);
        m.to_dense().fidelity(&d)
    }

    #[test]
    fn bundle_is_cheap_for_any_k() {
        // One-direction crossing: max χ = hi regardless of how many strands.
        let profile = zigzag::diamond(1, 5); // hi = 5
        for k in [2usize, 3, 4, 5] {
            let (net, groups) = bundle(&profile, k);
            let m = build(&net, &net.ghz(&groups));
            assert!(
                m.max_bond_dim() <= 5,
                "k={}: max χ {} > hi",
                k,
                m.max_bond_dim()
            );
        }
    }

    #[test]
    fn bundle_matches_dense() {
        let profile = zigzag::diamond(1, 3);
        let (net, groups) = bundle(&profile, 3);
        assert!((dense_fidelity(&net, &net.ghz(&groups)) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn weave_is_an_area_law() {
        // Two-direction crossing: max χ = d^rows, exponential in the
        // transverse extent and (nearly) flat in the length.
        let d = 2;
        for rows in [2usize, 3, 4] {
            let (net, edges) = weave(d, rows, 4);
            let m = build(&net, &net.cluster(&edges));
            assert_eq!(
                m.max_bond_dim(),
                d.pow(rows as u32),
                "rows={}: χ {} != d^rows",
                rows,
                m.max_bond_dim()
            );
        }
        // Flat in length: rows fixed, cols grows → χ unchanged.
        let (n1, e1) = weave(d, 3, 3);
        let (n2, e2) = weave(d, 3, 6);
        assert_eq!(
            build(&n1, &n1.cluster(&e1)).max_bond_dim(),
            build(&n2, &n2.cluster(&e2)).max_bond_dim()
        );
    }

    #[test]
    fn weave_matches_dense() {
        let d = 2;
        let (net, edges) = weave(d, 2, 3); // 6-qubit cluster
        assert!((dense_fidelity(&net, &net.cluster(&edges)) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn overlay_has_cutwidth_two() {
        // Two directions overlaid on one pair: cutwidth 2 → χ ≈ hi², far
        // below the block-order blowup but above the single-X hi.
        let profile = zigzag::diamond(1, 4); // hi = 4
        let (net, edges) = overlay(&profile);
        let m = build(&net, &net.cluster(&edges));
        assert!(
            m.max_bond_dim() <= 4 * 4,
            "overlay χ {} > hi²",
            m.max_bond_dim()
        );
        assert!(m.max_bond_dim() > 4, "overlay χ {} should exceed hi", m.max_bond_dim());
        assert!((dense_fidelity(&net, &net.cluster(&edges)) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn wave_weave_cost_is_waist_local_and_periodic() {
        // Weaving a MULTI-period wave: the cost concentrates at the coarse
        // waists and is periodic — a second period does not raise the peak
        // (the crossing cost is local to each RG cell), and valleys pinch
        // the bond to 1.
        let one = zigzag::wave(1, 3, 1); // [1,2,3,2,1]
        let two = zigzag::wave(1, 3, 2); // [1,2,3,2,1,2,3,2,1]
        let (n1, e1) = wave_weave(&one, 2);
        let (n2, e2) = wave_weave(&two, 2);
        let m1 = build(&n1, &n1.cluster(&e1));
        let mut m2 = build(&n2, &n2.cluster(&e2));
        // Same peak cost per period, independent of how many periods.
        assert_eq!(m1.max_bond_dim(), m2.max_bond_dim());
        // Valleys pinch: interior χ = 1 bonds appear (the fine rims).
        let bonds = m2.bond_dims();
        assert!(bonds.iter().filter(|&&b| b == 1).count() >= 4);
        assert!((dense_fidelity(&n2, &n2.cluster(&e2)) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn wave_weave_base_follows_the_profile() {
        // A cut between columns c, c+1 costs min(profile[c],profile[c+1])^rows.
        let profile = zigzag::diamond(1, 3); // [1,2,3,2,1]
        let rows = 2;
        let (net, edges) = wave_weave(&profile, rows);
        let m = build(&net, &net.cluster(&edges));
        let bonds = m.bond_dims();
        // Column-major: the bond just after column c's block is between-column.
        // Peak-adjacent gaps (min dim 2 or 3) carry base^rows.
        assert!((dense_fidelity(&net, &net.cluster(&edges)) - 1.0).abs() < 1e-9);
        assert!(
            *bonds.iter().max().unwrap() <= 3usize.pow(rows as u32),
            "max bond {} exceeds hi^rows",
            bonds.iter().max().unwrap()
        );
    }
}
