# novel_quantum_structures

**Twisted zigzag qudit chains: a zero-dependency Rust laboratory for
heterogeneous-dimension quantum structures.**

This library explores a question: what happens if a qudit chain's local
dimension *expands and collapses* along the chain — qubit → qutrit → 4-dit →
5-dit → back down — with **cross-scale couplings** that pair sites across the
dimension wave as bidirectional partners? Is such a structure classically
computable, and can it sub-simulate ordinary binary quantum systems on either
side of its high-dimensional waist?

The short answers, all verified executable in this repository: **yes, via
heterogeneous-dimension tensor networks; yes, exactly; and the interesting
boundary is which circuit families keep the representation compressed.**

```text
  5           ●                 ●
  4         ●   ●             ●   ●
  3       ●       ●         ●       ●
  2     ●           ●     ●           ●        … the dimension wave
        └───────x─────────┘
            └───x───┘        mirror pairs: equal-dimension sites on opposite
              └─x─┘          flanks — full-SWAP / Bell-couplable partners
```

## The idea, mapped to known ground

The zigzag structure is a cousin of ideas from tensor networks and
renormalization:

* **Classical computability.** A chain state is stored as a *matrix product
  state* (MPS) — one 3-index tensor per site — whose cost is linear in chain
  length for bounded bond dimension χ, instead of the product of all local
  dimensions. Nothing in the MPS formalism requires uniform local dimension,
  so dimension waves are first-class citizens here; the profile
  `[2,3,4,5,4,3,2]` is just another chain.
* **The profile is an entanglement budget.** The maximal Schmidt rank across
  any bond is `min(∏ dims left, ∏ dims right)` — for the diamond, 24 at the
  central bonds. The dimension wave *caps* how much entanglement can cross
  each cut: geometry doubles as resource bound, the same intuition that makes
  MERA-style hierarchical ansätze efficient.
* **Cross-scale couplings are MPO applications.** A two-site gate between
  *any* pair of sites (adjacent or 24 sites apart, equal or unequal
  dimensions) is decomposed by an operator-Schmidt SVD into an MPO threaded
  through the intervening tensors, then recompressed. The "1↔5, 5↔1"
  bidirectional couplings are ordinary operations, not special cases.
* **Sub-simulating binary systems is subspace embedding.** A dimension-`d`
  site hosts `floor(log2 d)` logical qubits with the remaining levels as
  slack. Lifted gates never mix register and slack, so hosted qubit circuits
  are *exact*, not approximate.

## What is in the box

| module | contents |
|---|---|
| `c64`, `mat` | complex arithmetic, deterministic RNG (SplitMix64), Haar unitaries, and a one-sided Jacobi complex SVD — all in-crate, **zero dependencies** |
| `dense` | exact mixed-radix state-vector simulator (ground truth) |
| `mps` | heterogeneous-dimension MPS: canonical forms, adjacent + long-range gates (operator-Schmidt MPO), Schmidt spectra / bond entropies, truncation with discarded-weight accounting |
| `mps` (scale morphing) | `merge_sites`, `split_site`, `promote_site`, `demote_site` — the representation itself expands and collapses while preserving the state |
| `gates` | qudit gates incl. heterogeneous couplers: `cshift`, `cphase`, and `xswap` (subspace exchange — the bidirectional-pair gate) |
| `embed` | hosted-qubit registers: gate lifting, compilation of qubit circuits onto qudit chains, exact extraction, `qubit_bit_swap` shuttling |
| `zigzag` | diamond/wave profiles, mirror pairs, valleys & waists, named circuit families (`bowtie`, `crosswave_round`, `brickwork_random`) |
| `circuit` | backend-agnostic gate lists so every experiment cross-validates dense vs MPS |

## Findings (all reproducible from `examples/`)

**1. The bowtie state saturates the diamond's waist — cheaply.**
Bell-pairing the three mirror pairs of `[2,3,4,5,4,3,2]` produces bond
entropies `[1, 2.585, 4.585, 4.585, 2.585, 1]` bits — exactly
`[log2 2, log2 6, log2 24, …]` — and Schmidt rank 24 across the waist, the
*maximum this geometry admits*, at MPS bond dimension 24. Verified against
dense simulation to fidelity 1 − 10⁻¹⁵. (`diamond_bowtie`)

**2. Deep scrambling needs exactly the geometric budget, no more.**
A 50-gate random cross-scale circuit on the diamond is exact at χ = 24 and
degrades monotonically below it — the dimension wave is a hard entanglement
ceiling. (`diamond_bowtie`)

**3. The diamond hosts 10 exact logical qubits (2.8× overhead).**
Capacities `1,1,2,2,2,1,1` across the seven sites: the valley sites *are*
qubits — the binary sets on either side — and random 40-gate qubit circuits
compiled onto the chain reproduce a pure-qubit reference to machine
precision with zero slack leakage, on both backends. GHZ-10 costs χ = 2.
(`hosted_qubits`)

**4. Edge qubits and the waist register are bidirectional pairs.**
`qubit_bit_swap` shuttles both edge qubits *into* the 5-dit waist's 2-bit
register; one coarse gate there (`F₄` = the whole 2-qubit QFT) processes
them jointly; shuttling back reproduces the direct computation to fidelity 1.
"Acting at multiple scales" made literal. (`hosted_qubits`)

**5. Scale morphing is exact and gates transfer across scales.**
An entangled 8-qubit chain merges to `[4,4,4,4]`, then `[16,16]`, with the
surviving bond entropies numerically identical at every scale; splitting back
recovers the original at fidelity 1 − 10⁻¹³. Merge → `F₄` → split equals the
four-gate qubit-scale QFT circuit. (`scale_morphing`)

**6. A 25-site wave (dense: 8.6·10¹² amplitudes ≈ 138 TB) runs exactly in
0.76 MB.** A 62-gate structured cross-scale circuit — Fourier layers,
per-diamond mirror couplings, valley↔valley and edge↔edge couplings spanning
up to 24 sites, a controlled-shift staircase — executes with **zero
truncation** at χ = 48 in ~60 ms, with a periodic entanglement profile
riding the dimension wave. (`long_wave`)

**7. The compression boundary is structure, not lattice.** Haar-random
brickwork on the same wave blows through any fixed χ (severe truncation at
χ = 32). The zigzag's classical computability is a property of *structured
multi-scale dynamics* — the same line MERA draws between renormalizable and
volume-law circuits. (`long_wave`)

## Quick start

```rust
use novel_quantum_structures::{mps::Mps, zigzag, TruncSpec};

// The 2→3→4→5→4→3→2 diamond, Bell-paired across its mirror flanks.
let profile = zigzag::diamond(2, 5);
let mut psi = Mps::zero_state(&profile, TruncSpec::new(64, 1e-12));
zigzag::bowtie(&profile).run_mps(&mut psi);

assert_eq!(psi.max_bond_dim(), 24);          // waist-saturating, yet tiny
let bits = psi.bond_entropies_bits();        // [1, 2.585, 4.585, ...]
```

Long-range cross-scale couplings and scale morphing:

```rust
use novel_quantum_structures::gates;

psi.apply2(0, 6, &gates::xswap(2, 2));       // exchange the edge qubits
psi.merge_sites(2);                          // 4·5 → one 20-dit site
psi.split_site(2, 4, 5);                     // and back
```

## Running

```sh
cargo test                                   # 49 tests, dense-vs-MPS cross-validation
cargo run --release --example diamond_bowtie
cargo run --release --example hosted_qubits
cargo run --release --example scale_morphing
cargo run --release --example long_wave
```

No dependencies; builds with any reasonably recent stable Rust.

## Design notes

* **Verification-first.** Every mechanism is cross-checked against the exact
  dense simulator on small chains (49 tests), including randomized circuits
  over both orientations of long-range gates, canonical-form invariance, and
  analytic entropy values.
* **In-crate numerics.** The SVD is a one-sided Jacobi with two
  hardening details worth knowing: numerically-null columns are frozen
  before rotations grind them into the denormal range (a pathology that
  corrupts the significant subspace), and singular values below roundoff
  resolution are reported as exact zeros rather than normalized into fake
  directions.
* **Truncation is honest.** Every discarded Schmidt weight is tallied on the
  state (`discarded_weight`); exact runs report exactly 0.

## Open directions

* **Which circuit families stay compressed?** The structured/scrambling
  boundary measured in `long_wave` invites a systematic study: Clifford-like
  qudit circuits, cross-scale couplings only, dimension-commensurate gates
  (`gcd`-respecting `cshift`/`cphase` webs)…
* **`xswap` as a disentangler.** For mirror-symmetric correlations the
  subspace exchange can relocate entanglement toward the waist before
  truncation — a MERA-style disentangler adapted to the wave. The
  diagnostics (`schmidt_spectra`) are in place to measure whether it earns
  its keep.
* **Dynamic profiles.** `promote`/`demote`/`merge`/`split` allow the wave
  itself to evolve during a computation — an adaptive-geometry simulator
  where the dimension profile tracks where entanglement wants to live.
* **Beyond chains.** Mirror pairs hint at a ladder; the natural next
  structure is a tree or bowtie *network* of dimension waves.
