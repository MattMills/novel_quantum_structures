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

A second question follows: can the *operations themselves* become pipeline
objects — states stretched between the past and the future boundary of a
circuit block, composed and recursed as first-class values? Also yes: the
`mpo`/`radix` layers reify circuit blocks as operators, and the structured
**tail-radix phase web** (the mixed-radix Fourier transform over the chain's
own ring) gives them something real to compute: `QFT† ∘ ramp ∘ QFT`
composed as three blocks collapses into a modular adder, adders compose
recursively into adders, and the whole scheme runs — with one measured,
instructive precision boundary — on rings of 8.6 trillion elements
(findings 8–11 below).

And a third: is an operation one object, or a *flow* — divisible into ever
finer temporal grain? The `flow` layer extends every block to an exact
one-parameter group `t ↦ U^t` and drives an **operation-width cursor** over
it: a pipeline of flow segments refined at the focus and coarsened behind
it, with the total invariant. The headline measurements: bond dimension
along the flow *sees the integers* (χ = 3 at whole shifts, 9–11 between),
the same coarse block admits Hilbert–Schmidt-*orthogonal* half-imputations
of different widths, and states walking the refined pipeline follow an
exact drift-plus-lattice-wobble law (findings 12–14 below).

Finally, the operators that *compute* here are **stepwise cascades**:
finite-state transducers whose message rides the MPO bond. Modular
multiplication joins the adder as an exact cascade, the family composes
recursively toward modular exponentiation, and the measured operator width
of `×k` turns out to be number-theoretically resonant — non-monotone in k
(findings 15–17 below).

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
| `mpo` | **circuit-of-circuits**: circuit blocks as matrix product operators — Choi-carrier states on the doubled chain, block composition (`compose_after`), one-shot application (`apply_to`), operator-entanglement diagnostics |
| `radix` | the **tail-radix phase web**: mixed-radix QFT over the chain's ring `Z_N` (standard and reversal-free), Draper phase ramps, the exact bond-2 carry adder MPO |
| `flow` | **operator flows and the operation-width cursor**: exact fractional powers `U^t` via the Fourier frame (`FourierFlow`), and `WidthCursor` — a pipeline of flow segments at adaptive temporal resolution (refine/coarsen with invariant total) |
| `cascade` | **stepwise cascade operators**: finite-state transducers lifted to MPOs with the message riding the bond — the carry adder (m = 2), modular multiplication `x → kx mod N` (m = k, unitary iff gcd(k, N) = 1), recursive composition |
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

**8. The chain's ring has a structured tail-radix phase web.** The Fourier
transform over `Z_N` (`N = Π d_i`) factorizes exactly so that output digit
`j` couples *only to the tail* of the register (sites `i ≥ j`), through
fractional controlled phases of angle `2π/(d_j·d_{j+1}···d_i)` — the radix
segment product. Angles decay super-exponentially with distance, so the web
is effectively banded. Verified digit-by-digit against the dense DFT kernel.
(`radix`, `circuit_of_circuits`)

**9. Circuit blocks are pipeline states, and composition collapses them.**
The QFT over `Z_2880` compiles to an MPO with bond dims `[2,6,8,7,6,2]`
(≈ 0.8 bits of past↔future correlation per cut). Composing three blocks as
operators — `QFT† ∘ ramp ∘ QFT` — *collapses* the pipeline into the modular
adder at bonds `[2,4,4,4,2,1]`: the composite is simpler than its factors,
and recompression finds that automatically. Verified against the exact
cyclic shift. A measured aside: the QFT *core* is operator-cheap (χ = 6–8);
it is the digit-reversal stage — the bowtie mirror permutation — that costs
χ = 36, and the adder pipeline cancels it between `V` and `V†`.
(`circuit_of_circuits`, `tests/pipeline.rs`)

**10. Operator recursion is closed and runs both ways.** Adders composed
with adders stay adders (`A_7 ∘ A_3 = A_10`, χ flat at 3 across six
self-compositions, each verified), and conjugating back through the Fourier
frame recovers the bond-1 diagonal ramp. Gate composition has become
circuit-of-circuits composition. (`circuit_of_circuits`)

**11. Two budgets govern operator composition: entanglement *and*
precision.** On the 25-site wave the same pipeline builds the adder over
`Z_8599633920000` as a χ = 6 operator (~84 KB) that is machine-accurate on
typical inputs — but the deepest carry `|N-1⟩ → |0⟩` rides tail-radix
couplings of angle `2π/N ≈ 7·10⁻¹³` whose relative operator weight `~1/N`
falls below any `f64` cutoff, and that single column of the operator is
lost (Hilbert–Schmidt fidelity vs the exact adder is still 1 − 10⁻⁹:
average-case operator fidelity is not worst-case). The collapsed form —
a hand-built bond-2 carry-propagation MPO, `radix::adder_mpo_exact` — is
exact at any ring size, *provided it is compressed by rank, never by
weight*: its deep-carry branch has relative Frobenius weight `~1/N`, free
to keep in rank but fatal to cut. (`circuit_of_circuits`,
`radix::exact_carry_mpo_handles_the_deepest_carry_at_scale`)

**12. Operators are flows, and width sees the integers.** Every block
extends to an exact one-parameter group `U^t = V†·D(t)·V` (fractional
powers are free in the Fourier frame; `U^s ∘ U^t = U^{s+t}` verified at
10⁻⁸). Scanning the continuous shift over `Z_2880` for `t ∈ [0, 2]`: bond
dimension collapses to 3 exactly at integer widths (crisp permutations) and
sits at 9–11 in between (delocalized sinc kernels), with operator
entanglement peaking mid-flow — **arithmetic quantization read directly off
operation width**. (`flow`, `operator_width_cursor`)

**13. Imputation has geometry.** The same coarse block factors into
geodesic halves `U^½·U^½` (χ = 11 each) or causal halves through the
Fourier frame (χ = 8 each); both products reproduce the block at fidelity
1.000000000, yet the two midpoints are *Hilbert–Schmidt orthogonal*
(overlap 0.000000) — genuinely different curves through operator space
with identical endpoints. Decomposing a coarse operation into finer ones is
a geometric *choice*, with measurable width costs per route.
(`operator_width_cursor`)

**14. The operation-width cursor works.** A pipeline over the flow refines
dyadically at its leading edge — coarse past, increasingly fine present —
with the total composition invariant (fidelity 1.000000000 at every zoom
level). A state stepped through the refined pipeline drifts along the ring
following the exact first-moment law `⟨ω^x⟩ = ω^{x₀+t}·((N−1)+e^{−2πit})/N`
— linear motion plus a **lattice wobble** of amplitude `1/2π` (the ring's
discreteness pushing back on the continuous flow), participation breathing
`1 → ⅔ → ⅓ → ⅔ → 1`, and exact relocalization at whole width.
(`flow::WidthCursor`, `operator_width_cursor`)

**15. Stepwise cascades are a first-class operator family.** A transducer —
a finite-state machine sweeping the chain with its message riding the MPO
bond — lifts to an operator whose bond dimension *is* the width of the
classical information front. The carry adder is the `m = 2` member; modular
multiplication `x → kx mod N` is the `m = k` member (the digit-carry rule
`(digit, carry) ↔ k·digit + carry` is a bijection at every site). `×7` on
the 8.6-trillion-element wave builds in under a millisecond with the
worst-case carry exact; the family composes recursively (`×7 ∘ ×11 = ×77`
with the composed and directly-built operators matching bond-for-bond, and
adders braiding with multipliers into affine maps). (`cascade`,
`stepwise_cascade`)

**16. The true width of modular multiplication is number-theoretic — and
non-monotone.** Repeated squaring `7 → 49 → 2401 → 1921 (mod 2880)` gives
measured operator widths `7 → 24 → 6 → 3`. The message bound `k` collapses
across each cut: with right-side size `S` and left ring `L = N/S`, the
carry entering the cut is `c = ⌊kb/S⌋` and only `c mod L` matters — at the
diamond's waist (`S = 120, L = 24`), `×2401` has `c = 20b`, taking
`24/gcd(20,24) = 6` values, and `×1921` has `c = 16b`, taking
`24/gcd(16,24) = 3`. **Multiplying by 1921 is a thinner operator than
multiplying by 7.** The cost curve of modular exponentiation on a chain is
resonant, not monotone. (`stepwise_cascade`)

**17. Fourier conjugation reverses the cascade; gcd breaks unitarity.**
`V ∘ (×k) ∘ V†` equals `×k⁻¹` *with the carry sweep running in the opposite
direction* (the scaling theorem transported through the reversed-digit
frame; verified at fidelity 1 on a non-palindromic chain). And a cascade is
unitary exactly when its arithmetic is invertible: `×6` on `Z_2880`
(gcd = 6) shows unitarity defect 0.833 = 1 − 1/6 with an explicit image
collision `|0⟩, |480⟩ → |0⟩` — number theory surfacing as an operator
property. (`cascade`, `stepwise_cascade`)

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

Circuit-of-circuits — blocks composed as operators:

```rust
use novel_quantum_structures::{mpo::Mpo, radix};

let spec = TruncSpec::new(96, 1e-12);
let v = Mpo::from_circuit(&radix::mixed_radix_qft_reversed(&profile, 0.0), spec);
let d = Mpo::from_circuit(&radix::fourier_phase_ramp_reversed(&profile, 1234), spec);
let adder = v.adjoint().compose_after(&d.compose_after(&v, spec), spec);
// adder.bond_dims() == [2, 4, 4, 4, 2, 1] — the pipeline collapsed.
let shifted = adder.apply_to(&psi);          // |x⟩ → |x + 1234 mod 2880⟩
```

## Running

```sh
cargo test                                   # 78 tests, dense-vs-MPS cross-validation
cargo run --release --example diamond_bowtie
cargo run --release --example hosted_qubits
cargo run --release --example scale_morphing
cargo run --release --example long_wave
cargo run --release --example circuit_of_circuits
cargo run --release --example operator_width_cursor
cargo run --release --example stepwise_cascade
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
  directions. Strongly rank-capped truncations of large matrices go through
  a randomized range-finder (Halko-style, one power iteration) whose
  unsampled weight is charged honestly to the discard tally.
* **Truncate only against canonical environments.** Raw (freshly
  contracted) bond coordinates distort Schmidt weights — a rank cap applied
  there can slice through degenerate subspaces that carry real weight, a
  bug this library hit and now tests against. `Mps::recompress` therefore
  shrinks bonds exactly by a meet-in-the-middle sweep choreography (no SVD
  ever runs at raw-product width on both sides) and applies the lossy pass
  only in canonical form.
* **Truncation is honest.** Every discarded Schmidt weight is tallied on the
  state (`discarded_weight`); exact runs report exactly 0.

## Open directions

* **Which circuit families stay compressed?** The structured/scrambling
  boundary measured in `long_wave` invites a systematic study: Clifford-like
  qudit circuits, cross-scale couplings only, dimension-commensurate gates
  (`gcd`-respecting `cshift`/`cphase` webs)…
* **Fourier-space arithmetic beyond addition.** The adder pipeline
  (finding 9) begs for multiplication: `x → k·x mod N` as a composed block,
  and from there modular exponentiation — the Shor kernel — as a circuit of
  circuits over the zigzag's ring, with operator entanglement as the cost
  meter.
* **`xswap` as a disentangler.** For mirror-symmetric correlations the
  subspace exchange can relocate entanglement toward the waist before
  truncation — a MERA-style disentangler adapted to the wave. The
  diagnostics (`schmidt_spectra`) are in place to measure whether it earns
  its keep.
* **Temporal pipelines.** The width cursor is the first rung: flow
  segments at adaptive grain. The process-tensor generalization —
  boundaries spanning *several time slices*, contracted along the space
  axis — is the next.
* **Flows beyond the shift family.** `FourierFlow` diagonalizes the adder
  family; fractional powers of the QFT itself (the fractional Fourier
  transform over `Z_N` as an MPO flow) need `V`'s own eigenframe — open,
  and the natural test of whether "width sees structure" generalizes.
* **Width-optimal imputation.** Finding 13 shows factorization routes have
  different width costs; searching decomposition space for minimal-width
  pipelines (the causal detour beat the geodesic here) is an optimization
  problem this library can now pose concretely.
* **Higher-precision carriers.** Finding 11 shows `f64` spectral weight is
  a real resource boundary at `N ≳ 10¹²`; a `f128`/double-double carrier
  would push Fourier-composed arithmetic several orders further.
* **Dynamic profiles.** `promote`/`demote`/`merge`/`split` allow the wave
  itself to evolve during a computation — an adaptive-geometry simulator
  where the dimension profile tracks where entanglement wants to live.
* **Beyond chains.** Mirror pairs hint at a ladder; the natural next
  structure is a tree or bowtie *network* of dimension waves.
