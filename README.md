# novel_quantum_structures

**Twisted zigzag qudit chains: a zero-dependency Rust laboratory for
heterogeneous-dimension quantum structures.**

> **[THEORY.md](THEORY.md)** derives the mathematics behind every finding
> below — definitions, propositions, and proofs, each tied to the module
> that implements it and the example that measures it.

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
— and every machine has a *geometrically opposed dual* sweeping the other
way, often exponentially narrower (findings 15–18 below).

And a fifth: close the boundary loop, and the machine becomes a
*self-stabilizing dynamical system*. The traced adder computes in a prime
field that heals gcd obstructions, the double-zero seam is an attractor
that repairs the number representation with an exact `1/(k²+1)` law, and a
rank-one damper makes the repair exponential at a rate you choose
(finding 19 below).

Three further instruments close the loop between the theory document and
the laboratory: the **width atlas** computes operator widths from pure
number theory before any tensor exists (finding 20), **boundary twists**
re-wire a machine's message loop so one body computes mod `N−1`, `N`, or
`N+1` (finding 21), and the **frame flow** takes fractional powers of the
QFT itself out of its own operator algebra, `F⁴ = I` (finding 22).
Controlled cascades then run the **full Shor kernel** — controlled
modular exponentiation, phase estimation on one reused valley qubit,
order recovery by continued fractions — priced bond-for-bond in advance
(finding 23).

And the first steps past a single chain: **crossing strands** and
**networks**. Two dimension waves `[1,2,3,4,5,4,3,2,1]` laid antiparallel
into an X, coupled where their dimensions match — and *which* level (5, 4,
3, 2, 1) carries the coupling selects the entanglement geometry: the
shoulder injects more than the peak, the pinch nothing, and the cost is a
choice of site ordering, not a property of the physics (finding 24). Push
to many strands and arbitrary coupling graphs and one law organizes all of
it — the cost is the graph's **cutwidth**: one crossing direction is flat
(bundles, multi-period crossings stay `χ = hi` however many strands or
crossing points), a second is an area law (a woven lattice is
`χ = d^min(rows,cols)`) — the MPS/PEPS boundary, measured (finding 25).

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
| `mpo` | **circuit-of-circuits**: circuit blocks as matrix product operators — Choi-carrier states on the doubled chain, block composition (`compose_after`), one-shot application (`apply_to`), operator-entanglement diagnostics; `select_on` — control-selection `Σ_a \|a⟩⟨a\| ∘ U^a`, making controlled blocks (the `C-U^{2^j}` of phase estimation) first-class |
| `radix` | the **tail-radix phase web**: mixed-radix QFT over the chain's ring `Z_N` (standard and reversal-free), Draper phase ramps, the exact bond-2 carry adder MPO |
| `flow` | **operator flows and the operation-width cursor**: exact fractional powers `U^t` via the Fourier frame (`FourierFlow`), and `WidthCursor` — a pipeline of flow segments at adaptive temporal resolution (refine/coarsen with invariant total); `FrameFlow` — fractional powers of the **QFT itself** as the four-term projector combination `F^t = Σ c_m(t)·F^m` (`F⁴ = I`), no eigensolver needed |
| `cascade` | **stepwise cascade operators**: finite-state transducers lifted to MPOs with the message riding the bond — the carry adder (m = 2), modular multiplication `x → kx mod N` (m = k, unitary iff gcd(k, N) = 1), recursive composition, and the **geometrically opposed dual machine** (`div`: MSB-first remainders with superposed-entry/postselected-exit boundaries, computing `×k⁻¹` at width k) |
| `stabilize` | **self-stabilizing boundary systems**: looped message boundaries (`to_mpo_looped`) and **twisted loops** (`to_mpo_looped_twisted` — the reversal twist selects diminished-one arithmetic mod `N+1`, completing the ring family `N−1 / N / N+1`), iteration dynamics with an attractor, plus operator linear combinations (`Mpo::add`/`scale`/`basis_transfer`) |
| `width` | the **a-priori width calculus**: the cut-rank theorem `χ_cut(×k) = \|{⌊kb/S⌋ mod L}\|` as executable number theory — per-bond width profiles of modular multiplication computed with no tensors, pinned bond-for-bond against recompressed cascade MPOs |
| `crossing` | **crossing dimension-wave strands**: two waves sharing one MPS, crossed pairwise into an X (`antiparallel`) or ladder (`parallel`) and coupled at a chosen dimension level; `Layout::{Block, Interleaved}` — the A|B entanglement is one bond in Block, a `≤ hi` near-product in Interleaved (a ~1000× cost knob for the same state) |
| `network` | **networks of crossing strands**: `K` strands and any coupling graph, with the law that MPS cost = `d^cutwidth`; `bundle` (one direction, `χ = hi` for any `K`), `overlay` (two directions on a pair, `(hi−1)²`), `weave`/`wave_weave` (a 2D lattice — area law `d^min(rows,cols)`); GHZ and cluster graph-state couplers |
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

**18. Every machine has a geometrically opposed dual — and it can be
exponentially narrower.** Multiplication's carries sweep LSB→MSB;
division's remainders sweep MSB→LSB with the *same* state-set size `k`.
Since `÷k = ×k⁻¹`, the inverse multiplier — whose forward carry machine
needs `k⁻¹ mod N` states — runs as the opposed remainder machine at width
`k`: on `Z_2880`, `×823` at width 7 instead of 823; on the 25-site wave,
width 7 instead of ≈ 6·10¹². The enabling ingredient is quantum boundary
conditions (`Boundary::SumAll` entry, `Boundary::Fixed(0)` exit): the
remainder front *enters in superposition over the unknown wrap multiple*
and is *postselected on exact division* — for each input exactly one
branch survives, something no classical FSM can do. Corollaries, all
verified at fidelity 1: opposed fronts annihilate (`÷7 ∘ ×7` collapses to
the all-ones bond profile), two narrow opposed fronts realize enormous
single-front machines (`×2357 = ÷11 ∘ ×7`, measured true width
`[2,6,17,13,6,2]`), and the second geometric opposition — ring reflection
by digit complement — gives negation `x → −x` at machine width 2 and
`×(N−7)` at width 8 instead of `N−7`. The operator is only as wide as the
*narrowest* machine that computes it. (`cascade::Transducer::div`,
`opposed_fronts`)

**19. A closed boundary loop is a self-stabilizing dynamical system.**
Feeding a cascade's exiting message back into its entry (`to_mpo_looped`)
turns the machine non-unitary — and non-unitarity is what an attractor
needs. Four exact results on the diamond: (a) the traced adder is
end-around carry, computing `mod N-1`, so the composite ring `Z_2880`
becomes the **prime field `Z_2879`**; (b) in a prime field there is no gcd
obstruction — `×6`, which is 6-to-1 (unitarity defect 0.833) on `Z_2880`,
is exactly unitary (defect 10⁻¹⁶) when looped onto `Z_2879`: the boundary
*heals* the operator; (c) the two representatives of zero make the traced
identity the defective operator `I + |0⟩⟨N-1|`, whose iteration pumps
negative-zero to canonical zero with error exactly `1/(k²+1)` — a Jordan
block, polynomial; (d) composing a rank-one damped seam
`I − (1−γ)|N-1⟩⟨N-1|` turns the Jordan block into a genuine eigenvalue and
makes stabilization exponential — 1001 iterations at γ = 1 collapse to
20 / 11 / 7 at γ = 0.5 / 0.2 / 0.05. Self-stabilization is engineered at
the boundary, not inside the operator. (`stabilize`,
`self_stabilizing_boundaries`)

**20. The width of a multiplication operator is computable before the
operator exists — and resonance survives at trillion scale.** The
cut-rank theorem (THEORY.md §8) reduces the operator Schmidt rank of `×k`
across any cut to counting `|{⌊kb/S⌋ mod L}|` — executable number theory,
no tensors (`width`). The full atlas of `Z_2880`: 668 of the 768
invertible multipliers (87%) saturate the waist cap 24; 100 are resonant,
down to width 2 (`×1441`, `×2879 = −1`) and width 3 (`×961`, `×1439`,
`×1921`). Predictions match exactly-recompressed cascade MPOs
bond-for-bond — `×7`, `×49 = ×7∘×7`, `×2401 = ×49∘×49` on diamond(2,4),
and all 24 bonds of `×7` on the wave. The sharpened minimal-machine
statement: **584 of 768 multipliers have intrinsic width below every
single-front machine alphabet** `min(k, k⁻¹, N−k, (N−k)⁻¹)` by a factor
> 8 (extreme: `×1441`, width 2 vs best alphabet 1439 — a 719× gap), so
their narrow presentations are reachable only through composition +
recompression: no single sweep of the chain carries so little. At scale,
the squaring orbit of 7 mod `8.6·10¹²` has central-cut widths
`7 → 49 → 2401 → 2 073 600 (cap) → 1 077 111 → 1 928 737` — the cost
curve of modular exponentiation stays number-theoretically resonant on a
trillion-element ring, and it is computed a priori. (`width`,
`width_atlas`)

**21. Twisting the loop selects a third ring: one machine computes mod
N−1, N, and N+1.** Feeding the exiting message back through the
*reversal* permutation (`to_mpo_looped_twisted`) turns the carry machines
into **diminished-one arithmetic mod `N+1`** — the encoding of
Fermat-number-transform hardware: chain value `x` represents `v = x+1`,
the twisted `+c` maps `v → v+(c+1) mod N+1`, the twisted `×k` maps
`v → k·v mod N+1` exactly, and the unrepresentable zero of `Z_{N+1}` is
an *annihilated branch* — a **hole**, the dual of the straight loop's
double-zero **seam** (the same input `x = N−1−c` hits both defects). On
`[2,3,4,3]` (`N = 72`, flanked by the twin primes 71 and 73) both closed
rings are fields: `×6` (defect 0.833 open) is healed by either closure
(defect ≤ 10⁻¹⁶). The twist can also *break*: `×5` is unitary on `Z_24`
but defective on the twist ring `25 = 5²` (defect 0.833, four annihilated
holes), and the diamond's twist ring `2881 = 43·67` shows the
factorization branch-by-branch (`|0⟩` and `|67⟩` collide on `|42⟩`;
`v = 67` annihilates). The twisted traced identity is the unilateral
shift `Σ|x+1⟩⟨x|` — a *drain* into the missing zero obeying the exact law
`‖M^k·uniform‖² = (N−k)/N` — where the straight trace was a pump.
Primality of `N∓1` is the design criterion, selectable by profile.
(`cascade::Transducer::to_mpo_looped_twisted`, `boundary_twists`)

**22. The QFT is itself a flow, spanned by four operators.** `F⁴ = I` on
any ring, so the QFT's spectral projectors are polynomials in `F` and the
fractional Fourier transform is the exact four-term combination
`F^t = c₀(t)·I + c₁(t)·F + c₂(t)·Π + c₃(t)·F†` (`Π` = the ring
reflection `x → −x`, a width-2 machine) — assembled by operator linear
algebra, no eigensolver (`flow::FrameFlow`). Measured on diamond(2,4): an
exact one-parameter group of period 4 (`F^½ ∘ F^½ = F`,
`F^0.7 ∘ F^1.3 = Π`, `F^3.5 ∘ F^0.5 = I`, each at fidelity
1.000000000), unitary at every `t`. The frame-quantization signature
refines finding 12: operator entanglement is *pinned* to the pure-power
value at every integer (extrema 0 / 5.170 / 1.000 / 5.170 bits), but bond
dimension collapses only at `t ≡ 0, 2 (mod 4)` (χ = 1, 2 against the
flat cap 36 elsewhere) — `F` itself saturates the geometry, so **width
sees the frame exactly where the frame is narrower than the geometry**.
Participation of `F^t|x₀⟩` breathes with period 2: localized at even `t`,
maximally flat (1/N) at odd `t`. (`flow::FrameFlow`,
`fractional_fourier`)

**23. The Shor kernel runs end-to-end — controlled, priced, and
bond-free where it counts.** `Transducer::mult_skipping` (the carry front
tunnels through a control site, leaving `I ⊗ ×k` at bond 1 across it) and
`Mpo::select_on` (projector-branch composition) make controlled modular
arithmetic first-class. The exact cost obeys the **controlled-width
formula** `χ = w + [k ≢ 1 mod S]` (THEORY.md Prop 8.11) — and the deep
QPE powers are resonant on every cut, so **their control is free**: the
`C-U^{2^j}` family for 7 mod 1440 on the diamond measures
`[2,4,8,8,6,2] / [2,4,13,24,6,2] / [2,3,3,3,3,2] / [2,3,3,3,3,2]`, every
bond as predicted, with `961 = 7⁴` and `481 = 7⁸` satisfying
`k ≡ 1 (mod S)` at every cut. Phase kickback is **bond-free**: on an
eigenstate the control cut stays at `χ = 1` while the phase lands
(30°/60°/120°/240° measured exactly); a non-eigenstate register costs
`χ = 2` — the resource statement of phase estimation, read off a bond
dimension. Semiclassical Kitaev phase estimation — the single valley
qubit at site 0 reused for all 10 bits, measurement by projection,
feedback rotations conditioned on earlier bits — plus continued fractions
bounded by the Carmichael exponent `λ(1440) = 24` recovers
**`ord(7 mod 1440) = 12`** from six runs (measured fractions
`1/3, 5/12, 1/4, 1/3, 2/3, 5/6`; the register is never measured — each
run collapses onto an eigenstate through control backaction alone). Both
Shor registers live on one dimension wave — control in the valleys,
arithmetic in the interior — and the whole protocol runs in a fifth of a
second. (`cascade::Transducer::mult_skipping`, `mpo::Mpo::select_on`,
`shor_kernel`)

**24. Two crossing waves entangle by shape, and their cost is a layout
choice.** Two strands `[1,2,3,4,5,4,3,2,1]` laid antiparallel into an X
(`A[i] ~ B[n−1−i]`) share one MPS and couple where their dimensions match.
Bell-coupling the pairs at a level set injects exactly `Σ log2 d` bits
across the A|B bipartition (Schmidt rank `Π d`) — verified against the
closed form and dense simulation for every level. The content is in the
*multiplicity*: a peak's waist (`d=5`) is **unique** but every lower level
is **paired**, so the most inter-strand entanglement enters at the
**shoulder** (`d=4`, two pairs, `4.00` bits) — not the peak (`d=5`, one
pair, `2.32` bits); the `d=1` pinch is a **decoupled crossing**
(`cshift(1,1) = I`). Cycles add: coupling at 5 then 4 then 3 then 2
accumulates `2.32+4.00+3.17+2.00 = 11.49` bits, the full budget one level
per cycle. A **valley** dual reflects the multiplicities (paired wide
ends, unique pinch), crossing richest at its rims. And the profiling
headline: the same crossing is `χ = Π d` (144, ~70 ms) in **Block** layout
but `χ ≤ hi` (4, ~30 µs — a **~1000× gap**) in **Interleaved** layout,
because a Bell crossing is a product of *local* pairs; a 50-site twin wave
fully crossed holds `34.5` bits of inter-strand entanglement in 8 KB at
`χ = 5`. The entanglement is physics; whether it is *expensive* is a
choice of site order — interleave to compute, read it off the block bond
(`Mps::bond_entropy_bits_at`, an `O(χ³)` local diagnostic).
(`crossing`, `crossing_vees`)

**25. Crossing networks: the cost is the coupling graph's cutwidth, and
the cutwidth is the dimensionality of the crossing.** Generalizing the X to
`K` strands and any coupling graph, the MPS bond dimension is exactly
`d^(edges crossing the worst cut)`, minimized over site orderings by the
graph's **cutwidth** (proved for graph-state couplers, dense-validated).
Reading that one integer sorts every richer geometry: a **bundle** of `K`
strands coupled in one direction (a GHZ rung per site) is cutwidth 1 —
`χ = hi` for `K = 2, 3, 5, 8`, *any* `K`; two multi-period `wave(1,5,p)`
strands cross at `p` peak-regions, and the inter-strand entanglement climbs
(`11.49 → 45.97` bits for `p = 1…4`) while `χ` stays pinned at `hi = 5`
(**many crossing points, one direction, still cheap**). **Overlaying** a
second direction on the same pair (X + ladder) unions two matchings into
4-cycles — cutwidth 2, `χ = (hi−1)²` (`4, 9, 16` for `hi = 3, 4, 5`; the
*shoulder* squared, since the unique peak carries only the ladder). A
**weave** — strands in two transverse directions, a 2D lattice — is
cutwidth `min(rows,cols)`, so `χ = d^rows` (`4, 8, 16, 32` for `rows = 2…5`,
*flat* in the length): an **area law**, the boundary where the 1D tensor
network stops being efficient (the MPS/PEPS line), and a **wave weave**
inherits it with a wave-shaped bond profile `[1,1,2,4,6,4,2,1,1]`. One
direction of crossing — however many strands or crossing points — is
constant cutwidth and classically cheap; a second transverse direction is
exponential. That threshold, not the strand count, governs simulability.
(`network`, `crossing_networks`)

**26. Crossing the renormalization wave: the coupling acquires a scale.**
The zigzag *wave* (many V's) is a real-space renormalization structure —
valleys fine, waists coarse, `merge_sites` a coarse-graining step,
structured dynamics riding the scale hierarchy (findings 5–7). Applying the
crossing to *these* strands, the crossing inherits the hierarchy, three
measured ways. **A crossing has a scale**: coupling two `wave(1,5,2)`
strands at the coarse waists is few crossings of *fat* modes (2 crossings,
5 modes), at the fine valleys many of *thin* modes (4 crossings, 2 modes).
**A woven lattice of waves carries an RG-shaped area law**: the bond profile
of a wave weave *is* the wave — `[1,1,2,4,6,4,2,1,1,1,2,4,6,4,2,1,1]`, two
humps at the two waists pinched to 1 at the valleys — thickening with each
row (`6 → 12` at 2 → 3 rows) and *periodic* (a second period doesn't raise
the peak; the cost is local to each RG cell). **A crossing is
RG-covariant**: coarse-graining a crossed strand toward its waist (an exact
`merge_sites`) carries the coupling up onto the coarse `d = 36` block and
splits back at fidelity `1.000000000000` — scale morphing and crossing
commute. The wave supplies a scale ladder, the crossing supplies coupling,
and everything expands into a complex multi-scale shape still costed by one
number — the cutwidth, now modulated by the renormalization structure.
(`crossing`, `network`, `renormalizing_crossings`)

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
cargo test                                   # 125 tests, dense-vs-MPS cross-validation
cargo run --release --example diamond_bowtie
cargo run --release --example hosted_qubits
cargo run --release --example scale_morphing
cargo run --release --example long_wave
cargo run --release --example circuit_of_circuits
cargo run --release --example operator_width_cursor
cargo run --release --example stepwise_cascade
cargo run --release --example opposed_fronts
cargo run --release --example self_stabilizing_boundaries
cargo run --release --example width_atlas
cargo run --release --example boundary_twists
cargo run --release --example fractional_fourier
cargo run --release --example shor_kernel
cargo run --release --example crossing_vees
cargo run --release --example crossing_networks
cargo run --release --example renormalizing_crossings
```

No dependencies; builds with any reasonably recent stable Rust.

## Design notes

* **Verification-first.** Every mechanism is cross-checked against the exact
  dense simulator on small chains (125 tests), including randomized circuits
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
* **The Shor kernel at scale.** Multiplication, recursive exponentiation,
  controlled powers, and order finding all run (findings 15, 20, 23). Next:
  coherent readout — a full QFT over a *hosted control register* (the
  wave's valley qubits, finding 3) instead of the semiclassical single
  qubit — and order finding on the 25-site wave, where the atlas already
  prices the `C-U` family over `Z_{8.6·10¹²}`.
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
* **Beyond chains.** The first step is taken — `crossing` puts two waves on
  one MPS as an X (finding 24), and shows the layout, not the entanglement,
  is the cost. The next structures are genuine *networks*: a ladder
  (parallel crossing), then trees and necklaces of waves, where no single
  interleaving localizes every coupler and choosing the site order becomes
  a combinatorial layout-optimization problem (the crossing pattern's
  tree-width). Do lattices of crossed waves obey an area law, and can
  `xswap` relocate inter-strand entanglement toward chosen cuts?
