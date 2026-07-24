# Theory notes

**The mathematics of twisted zigzag qudit chains — definitions, propositions,
and proofs behind the measurements in the [README](README.md).**

The README reports *what this library measures*; this document derives *why
those numbers are what they are*. It is written as a short monograph:
numbered propositions with proofs, each tied to the module that implements
the object and the example or test that measures the claim. Where a result
is stated without proof it is marked as **measured** — an executable fact of
the repository, reproducible from `examples/` — and where a derivation here
sharpens or corrects an idealized claim in the README, that is said
explicitly (see e.g. Proposition 9.6 on "orthogonal" imputations).

Nothing here claims quantum advantage. The subject is a classical
*representation theory of structured operators* — which states, circuit
blocks, machines, and flows admit narrow tensor-network presentations, what
the width of a presentation means, and what new dynamics appear when the
presentations' boundary conditions are opened, dualized, or closed into
loops. The quantum formalism (unitarity, entanglement, superposed boundary
vectors) is the *language* in which those questions become sharp.

### Reading map

| § | subject | modules | README findings |
|---|---------|---------|-----------------|
| 1 | chains, rings, digit frames, zigzag geometry | `zigzag`, `dense` | — |
| 2 | the entanglement budget | `mps`, `zigzag` | 1, 2, 6, 7 |
| 3 | heterogeneous MPS mechanics | `mps`, `mat` | — |
| 4 | hosted qubits | `embed` | 3, 4 |
| 5 | scale morphing | `mps` | 5 |
| 6 | operators as states | `mpo` | 9 (setup) |
| 7 | the tail-radix phase web | `radix` | 8–11 |
| 8 | stepwise cascades, duals, the atlas, the Shor kernel | `cascade`, `width`, `mpo` | 15–18, 20, 23 |
| 9 | operator flows, the width cursor, the frame flow | `flow` | 12–14, 22 |
| 10 | self-stabilizing boundary systems and twisted loops | `stabilize`, `cascade` | 19, 21 |
| 11 | crossing strands, networks, the cutwidth law, and the RG wave | `crossing`, `network` | 24, 25, 26 |
| 12 | numerical foundations | `mat`, `c64` | design notes |
| 13 | dictionary and open problems | — | open directions |

Throughout: sites are indexed `0 … n−1`, local dimensions `d_0 … d_{n−1}`,
`N = Π_k d_k`, and `ω = ω_N = e^{2πi/N}`. "Bond `m`" is the cut between
sites `m` and `m+1`. All logarithms in entropy statements are base 2.


## 1. Kinematics: chains, rings, and digit frames

**Definition 1.1 (chain, profile).** A *heterogeneous qudit chain* is a
Hilbert space `H = C^{d_0} ⊗ ··· ⊗ C^{d_{n−1}}` together with its
*dimension profile* `(d_0, …, d_{n−1})`, `d_i ≥ 2`. Nothing requires the
`d_i` equal; profiles that rise and fall are the objects of study.

**Definition 1.2 (digit frames).** The chain's computational basis is
identified with the ring `Z_N` in two ways.

* **Standard place value** — site 0 most significant:

  ```text
  P_i = Π_{k>i} d_k ,        x = Σ_i x_i · P_i .
  ```

* **Reversed place value** — site 0 least significant:

  ```text
  Q_j = Π_{k<j} d_k ,        y = Σ_j y'_j · Q_j .
  ```

The *digit reversal* `ρ : Z_N → Z_N` sends the number with standard digits
`(x_0, …, x_{n−1})` to the number with standard digits
`(x_{n−1}, …, x_0)`. It is an involution, and it is realizable by site-local
swaps exactly when the profile is palindromic (`d_i = d_{n−1−i}`), since a
digit moving from site `i` to site `n−1−i` must find a slot of its own
dimension there. Both frames are used constantly: the mixed-radix Fourier
transform (§7) naturally *produces* reversed-place digits, and choosing to
stay in that convention is what makes Fourier arithmetic profile-agnostic.

**Definition 1.3 (zigzag geometry, `zigzag`).** The *diamond*
`diamond(lo, hi)` is the profile `lo, lo+1, …, hi, …, lo+1, lo`; the *wave*
`wave(lo, hi, p)` glues `p` diamonds at their shared `lo`-valleys. For a
palindromic profile the *mirror pairs* are `(i, n−1−i)`; the *waists* are
the maximal-dimension sites, the *valleys* the minimal ones. On a
palindrome, mirror pairs have equal dimension — so they admit the full
subspace exchange `xswap(d, d) = SWAP` and generalized Bell pairing. These
are the "bidirectional pairs": couplings that jump across the dimension
wave without touching the intermediate scales.

The flagship instance is `diamond(2, 5) = [2,3,4,5,4,3,2]`, whose ring is
`Z_2880`; the flagship large instance is `wave(2,5,4)`, 25 sites, ring size
`N = 8 599 633 920 000 ≈ 8.6·10¹²` (a dense state vector would occupy
≈ 138 TB).


## 2. The entanglement budget

The first structural fact is that a dimension profile is not merely a list
of storage sizes — it is a *hard ceiling on entanglement*, bond by bond.

**Proposition 2.1 (rank ceiling).** Let `|ψ⟩ ∈ H` and let
`Λ_m = Π_{i≤m} d_i` be the left volume at bond `m`. The Schmidt rank of
`|ψ⟩` across bond `m` satisfies

```text
rank_m(ψ) ≤ B_m := min( Λ_m , N/Λ_m ),
```

and consequently the bond entropy satisfies `S_m ≤ log B_m`.

*Proof.* The Schmidt rank across a cut is the matrix rank of `ψ` reshaped
to a `Λ_m × (N/Λ_m)` matrix, and matrix rank is at most the smaller side
length. ∎

For `diamond(2,5)` the budget profile is

```text
B = [ 2, 6, 24, 24, 6, 2 ],
```

peaking at 24 across the two bonds flanking the waist. The profile is
literally an entanglement resource specification: geometry doubles as
budget, the same intuition that makes hierarchical (MERA-style) ansätze
efficient. On a wave the budget rises and falls periodically, which is why
the measured entanglement profile of structured dynamics "rides the wave"
(finding 6).

**Proposition 2.2 (bowtie spectrum).** Let the profile be a palindrome of
odd length with strictly increasing flank `d_0 < d_1 < ··· < d_c` (a
diamond). The *bowtie state* — the generalized Bell state
`|Φ_d⟩ = Σ_{j<d} |j,j⟩/√d` on every mirror pair, `|0⟩` on the center — has,
across bond `m` (take `m < c` by mirror symmetry):

```text
rank_m = Π_{i≤m} d_i = Λ_m = B_m ,      S_m = Σ_{i≤m} log d_i  bits,
```

with a flat Schmidt spectrum. In particular the bowtie **saturates the
budget of Proposition 2.1 at every bond**.

*Proof.* The state is a tensor product over pairs. A pair `(i, n−1−i)` with
`i ≤ m` has one leg on each side of the cut and contributes a factor with
Schmidt rank `d_i` and flat spectrum `1/√d_i`; a pair with `i > m` lies
entirely on the right of the cut (both `i > m` and `n−1−i > m` hold, the
latter because `m < c ≤ n−1−i`) and contributes rank 1, as does the center
site. Schmidt ranks and spectra multiply across tensor products, giving
rank `Π_{i≤m} d_i` and entropy `Σ_{i≤m} log d_i`. On the increasing flank
`Λ_m ≤ N/Λ_m`, so this equals `B_m`; bonds right of the center follow by
mirror symmetry. ∎

For `diamond(2,5)`: entropies `[1, log 6, log 24, log 24, log 6, 1] =
[1, 2.585, 4.585, 4.585, 2.585, 1]` bits and waist rank 24 — exactly the
measured values (finding 1), at MPS bond dimension 24 versus dense
dimension 2880. The state that is *maximally* entangled across the waist is
*cheap*: maximal entanglement under a geometric budget is small-rank by
definition of the budget.

**Necessity and the structure boundary (measured).** The budget is
necessary as well as sufficient for the diamond's dynamics: a deep random
cross-scale circuit is exact at `χ = 24` and degrades monotonically below
it (finding 2). But the budget alone does not make long chains cheap — on
the 25-site wave the mid-chain budget is astronomically large. What is
measured there (findings 6, 7) is the classic tensor-network dichotomy:

* a *structured* 62-gate cross-scale circuit (Fourier layers, mirror
  couplings, valley and edge couplings spanning up to 24 sites, a
  controlled-shift staircase) runs with zero truncation at `χ = 48`,
  0.76 MB of parameters;
* Haar-random brickwork on the same geometry blows through any fixed `χ`.

Compressibility is a property of *structured multi-scale dynamics* on the
geometry, not of the geometry alone — the boundary MERA draws between
renormalizable and volume-law circuits, and the reason polynomially-bounded
bond dimension implies classical simulability at all [1, 2, 3].


## 3. Heterogeneous MPS mechanics

This section fixes what the engine (`mps`) actually guarantees; readers
familiar with MPS folklore can skim, but two points (the canonical-gauge
truncation rule and the operator-Schmidt threading) carry real content.

**Definition 3.1.** An MPS stores `|ψ⟩` as site tensors `T_i[l, p, r]`
(left bond, physical, right bond) with per-site physical dimension `d_i`;
`|ψ(x_0…x_{n−1})⟩`'s amplitude is the product of the matrices
`T_i[·, x_i, ·]`. The *canonical-form invariant*: tensors strictly left of
the orthogonality `center` are left-isometries, tensors strictly right are
right-isometries, and the state's norm lives in the center tensor. Every
public operation restores this invariant.

**Why the gauge matters.** With a canonical environment, the singular
values of the center's matricization *are* the Schmidt coefficients of the
global state across that bond, so a local truncation is a globally optimal
projection (Eckart–Young through an isometric embedding) [1]. Without it,
freshly-contracted "raw" bond coordinates distort Schmidt weights, and a
rank cap can slice through a degenerate subspace that carries real weight —
a bug this library hit and now regression-tests. Hence the two operational
rules implemented by `Mps::recompress`:

1. **lossy truncation only against canonical environments**, and
2. **cost control by sweep choreography** ("meet in the middle": shrink
   bonds inward exactly from the cheap side first) rather than by early
   lossy cuts, so no SVD ever runs at raw-product width on both sides.

**Long-range gates (operator-Schmidt threading).** A two-site gate `G` on
arbitrary sites `a < b` is decomposed across the operator bipartition
`(a | b)`:

```text
G = Σ_{κ=1}^{K} A_κ ⊗ B_κ        (operator-Schmidt / SVD of the reshuffled G),
```

then applied as a width-`K` MPO: site `a` absorbs `A_κ` and opens a bond
index `κ`, intervening sites pass `κ` through diagonally, site `b` closes
it against `B_κ`; a two-pass (exact →, truncating ←) sweep over the window
recompresses. This makes the zigzag's cross-scale couplings — pairs 24
sites apart, unequal dimensions — ordinary operations. Correctness is
cross-validated against the dense backend over random circuits including
both orientations of long-range gates (`tests/backend_agreement.rs`).

**Truncation accounting.** A `TruncSpec` is a pair (hard rank cap, relative
discarded-weight cutoff). Every discarded Schmidt weight is tallied on the
state (`discarded_weight`), with the retained part renormalized; exact runs
report exactly 0. The tally is a per-truncation sum, not a global error
bound — where exactness matters the library measures fidelity against a
reference instead of trusting the tally.

**Linear combinations.** `Mps::add` implements the direct-sum construction
(bond dimensions add, then recompress). It looks like a utility; it is
load-bearing: operator sums, traced boundaries, and rank-one dampers (§10)
are all built on it.


## 4. Hosted qubits: exact sub-simulation

**Definition 4.1.** A site of dimension `d` *hosts* `κ(d) = ⌊log2 d⌋`
qubits: levels `0 … 2^{κ}−1` are read as a big-endian `κ`-bit *register*,
levels `≥ 2^κ` are *slack*. A qubit gate is *lifted* to the site by acting
as itself on the register block and as the identity on slack
(`embed::lift_register_gate`); two-qubit gates lift likewise, within a site
or across two sites.

**Proposition 4.2 (exactness).** Lifted gates are block-diagonal with
respect to the (register ⊕ slack) decomposition of every site, blocks are
preserved under products, and the restriction of a lifted circuit to the
global register subspace equals the reference qubit circuit. Hence a state
prepared in the register subspace stays in it (zero leakage), and the
hosted computation is *exactly* recoverable — the binary system is a view
into the qudit chain, not an approximation.

*Proof.* Block-diagonality is by construction of the lift; products of
block-diagonal operators are block-diagonal with multiplied blocks; the
register blocks are, by construction, exactly the reference gates. ∎

On `diamond(2,5)` the capacities are `1,1,2,2,2,1,1` — ten qubits in a
2880-dimensional space, overhead `2880/1024 = 2.8125` (finding 3). Two
readings are worth separating:

* the valley sites *are* qubits — the binary sets on either flank of the
  wave are literal, not encoded;
* interior sites host one to two logical qubits each *plus slack*, and the
  slack is not waste: it is the room into which the same site's algebra
  extends when the chain is used as a qudit machine (§7–§8).

**Bidirectional pairs.** Two gates realize the "pairing" language of §1:

* `xswap(da, db)` exchanges the shared `min(da,db)`-dimensional subspaces
  of two sites and fixes the rest — a permutation and an involution on the
  shared subspace. For equal dimensions it is the full SWAP; for `(2, 5)`
  it exchanges a physical qubit with the qubit-view of a 5-dit.
* `qubit_bit_swap(d, t)` swaps a physical qubit with *bit `t`* of a host
  register — a full SWAP on the embedded pair, identity on slack, and an
  involution.

The measured shuttle (finding 4) is then conjugation, nothing more: swap
both edge qubits into the waist register, apply one coarse gate there
(`F_4`, which on the register *is* the two-qubit QFT), swap back — equal to
the direct two-qubit computation at fidelity 1. Its significance is
semantic: a single `d = 5` site processes, in one native gate, a joint
operation on two qubits that live 6 sites apart — "acting at multiple
scales" is implemented by an exact isometric dictionary, not analogy.


## 5. Scale morphing

The representation itself can expand and collapse.

**Proposition 5.1.** (i) `merge_sites(i)` — contracting neighbouring
tensors into one site of dimension `d_i·d_{i+1}` with coarse index
`p_i·d_{i+1} + p_{i+1}` — is exact: it is the associativity isomorphism
`C^{d_i} ⊗ C^{d_{i+1}} ≅ C^{d_i d_{i+1}}` applied to the state.
(ii) `split_site(i, d1, d2)` inverts it by an SVD refactoring, exact up to
the state's truncation policy (a genuinely new bond appears, carrying the
entanglement that was intra-site). (iii) Applying a coarse gate `G` between
a merge and a split equals applying `G` as a two-site gate at the fine
scale. (iv) Bond entropies at surviving cuts are invariant under merging —
interior entanglement becomes intra-site structure and returns intact.

*Proof.* (i), (iii), (iv) are immediate from the associativity reading: the
state is unchanged as a vector in `H`; only the tensor factorization
recorded in `dims` changes, and a two-site gate *is* a one-site gate on the
fused factor. (ii) is the standard MPS split; exactness up to reported
truncation is the SVD. ∎

`promote_site` / `demote_site` move a single site up or down the dimension
ladder in place — an isometric embedding (new levels unpopulated), and a
projection that *reports its leakage* (the projected-out probability),
respectively. Together with (i)–(ii) these give the library its
"adaptive-geometry" surface: the profile is a dynamical variable, not a
type parameter (finding 5: an entangled 8-qubit chain morphs to `[4,4,4,4]`
to `[16,16]` and back at fidelity `1 − 10⁻¹³`, entropies identical at every
scale; merge → `F_4` → split equals the four-gate qubit-scale QFT).


## 6. Operators as states: the circuit-of-circuits layer

**Definition 6.1 (Choi carrier, `mpo`).** An operator `U` on the chain is
stored via vectorization (the Choi isomorphism [4]) as an ordinary MPS —
the *carrier* — over the doubled chain whose site `i` has dimension `d_i²`,
flattened index `p_out·d_i + p_in`. Read physically: each carrier site
holds a (future, past) leg pair, so the operator is literally **a state
stretched between the past boundary and the future boundary of a circuit
block**, both boundaries being full `n`-site slices.

Three consequences, all load-bearing:

* **Operator entanglement.** The entanglement entropy of the *normalized*
  carrier across bond `m` (`Mpo::operator_entanglement_bits`) is the
  operator entanglement of `U` at that cut [5] — the width of the
  correlation pipeline between past and future. For a unitary,
  `‖U‖_F = √N`, and zero operator entanglement everywhere characterizes
  product operators (e.g. the identity, or any diagonal phase ramp).
* **Everything is MPS algebra.** Compilation absorbs each gate `G` as its
  left-multiplication superoperator `G ⊗ I` on the doubled chain — so the
  full heterogeneous machinery of §3, long-range gates included, compiles
  circuits into blocks. Application to a state and composition of two
  blocks are site-wise bond contractions followed by `recompress`; the
  adjoint is a site-local leg swap plus conjugation; `add`, `scale`, and
  the rank-one `basis_transfer` `|to⟩⟨from|` make the space of operators a
  vector space with a computable basis of product operators.
* **Scalar-invariant comparison.** `hs_fidelity(A, B) =
  |tr(A†B)|² / (‖A‖²‖B‖²)` equals 1 iff `A = cB` — the workhorse equality
  test for operators held as carriers.

**The collapse principle.** Composition contracts two carriers and then
recompresses; recompression computes the *intrinsic* bond profile of the
product. Nothing forces a product to be as wide as its factors — and the
central measured phenomenon of the `mpo`/`radix` layers (finding 9) is that
a pipeline of blocks can collapse to something far simpler than its parts:

```text
QFT (χ up to 8)  ∘  ramp (χ = 1)  ∘  QFT† (χ up to 8)   →   adder (χ ≤ 4).
```

Composition-then-recompression is how the representation *does algebra*:
the simplification `V† D_c V = X^c` is discovered numerically, not
programmed. This is what "circuit of circuits" means operationally — a
compiled block is a gate at the next scale up, and the algebra of blocks is
executed by the same contraction machinery that executes gates.


## 7. The tail-radix phase web

The zigzag chain's ring `Z_N` carries a Fourier transform, and its
factorization over the chain has a rigid — and, for compression, very
fortunate — structure.

**Theorem 7.1 (tail-radix factorization).** Write inputs in standard place
value (`x = Σ x_i P_i`) and outputs in reversed place value
(`y = Σ y'_j Q_j`). Let `D_{j,i} = d_j · d_{j+1} ··· d_i` for `j ≤ i` (the
*radix segment product*). Then

```text
                     ┌ 1 / D_{j,i}          if j ≤ i
  P_i · Q_j / N  =   │
                     └ Π_{i<k<j} d_k  ∈ Z   if j > i ,
```

and therefore, modulo 1,

```text
  x·y/N  ≡  Σ_{j ≤ i}  x_i · y'_j / D_{j,i} ,
```

i.e. the DFT kernel `ω^{xy}` splits into pair couplings in which output
digit `j` is phase-correlated **only with the tail of the register**
(input sites `i ≥ j`), through the angle `2π/D_{j,i}`.

*Proof.* For `j ≤ i` the index sets `{k < j}`, `{j ≤ k ≤ i}`, `{k > i}`
partition `{0,…,n−1}`, so `N = Q_j · D_{j,i} · P_i`, giving the first case.
For `j > i` the sets `{k > i}` and `{k < j}` cover everything and overlap
exactly on `{i < k < j}`, so `P_i Q_j = N · Π_{i<k<j} d_k` — an integer
multiple of `N`, hence a vanishing phase. Summing `x_i y'_j P_i Q_j / N`
over all `(i, j)` and dropping integer terms gives the congruence. ∎

**Circuit realization** (`radix::mixed_radix_qft_reversed`). One
left-to-right pass: at site `j`, the local Fourier gate `F_{d_j}`
(realizing the diagonal `j = i` term and creating the `y'_j` superposition),
then fractional controlled phases `cp_frac` of angle `2π/D_{j,i}` from the
now-output site `j` to each still-input site `i > j`. The pass structure
works because couplings for output digit `j` involve only *unprocessed*
sites. For uniform `d = 2` this is the textbook QFT circuit with
bit-reversed output; Theorem 7.1 is its mixed-radix generalization, with
the *radix segment product* replacing the power of two.

**Corollary 7.2 (banding).** `D_{j,i} ≥ 2^{i−j+1}`, so coupling angles
decay at least exponentially in site separation — a factor `d_k ≥ 2` per
intervening site (on the wave, a factor 2–5). Dropping couplings below a
threshold angle gives the approximate transform, with monotonically
improving fidelity as the threshold tightens (tested); this is the
mixed-radix analogue of the banded approximate QFT [6], and it is *why the
web compresses*: the operator's long-range content is exponentially weak.

**Where the width lives (measured).** Compiled over `Z_2880`, the reversed
QFT block has bonds `[2,6,8,7,6,2]` — operator-cheap. The *digit reversal*
`ρ` needed for standard output order is the expensive part: it is exactly
the mirror-pair permutation of the zigzag (realized by the `xswap` network
on a palindrome — the bowtie geometry is load-bearing in the chain's own
Fourier transform), and it costs `χ = 36` on `diamond(2,4)`. The
reversal-free convention keeps every block compact, and in the arithmetic
pipeline below the reversal cancels between `V` and `V†` — it is never
paid at all (finding 9, aside).

**Fourier arithmetic.** The Draper ramp [7] is the diagonal block

```text
D_c = diag_y( e^{2πi·c·y/N} ) = ⊗_j diag_{y'_j}( e^{2πi·c·y'_j·Q_j/N} ),
```

a bond-1 product operator — and the factorization holds verbatim for
**real** `c`, which is what makes the translation *flow* of §9 exact. The
shift theorem `F† D_c F = X^c` (`X` the cyclic shift `|x⟩ → |x+1 mod N⟩`)
then gives:

**Proposition 7.3 (reversal-free adder).** With `V` the reversed-digit QFT
and `D_c` the reversed-place ramp, `V† ∘ D_c ∘ V` is the exact modular
adder `|x⟩ → |x + c mod N⟩` in standard digits, on **any** profile —
palindromic or not.

Composing the three blocks as MPOs collapses the pipeline (finding 9):
adder bonds `[2,4,4,4,2,1]`, simpler than either Fourier factor. The
family is closed and runs both ways (finding 10): `A_a ∘ A_b = A_{a+b}` —
trivial as permutation algebra, but the measured content is
*representational*: the bond profile stays flat (`χ ≈ 3`) across six
self-compositions, i.e. the compact presentation is closed under the
algebra, and conjugating back through the frame (`V A_c V†`) recovers the
bond-1 diagonal ramp.

### 7.4 The two budgets: rank and weight

Finding 11 is the most consequential numerical lesson in the repository,
and it deserves its precise statement. The adder is the permutation
operator `A_c = Σ_x |x+c mod N⟩⟨x|`, with `‖A_c‖_F = √N`. Consider the
single deepest-carry column, `|N−1⟩ → |c−1⟩` (for `c = 1`:
`|N−1⟩ → |0⟩`), the one input whose carry ripples through every site.

* **As spectral weight, it is invisible.** That column is one matrix
  element out of `N`: its relative squared Frobenius weight is `1/N` —
  `≈ 1.2·10⁻¹³` on the 25-site wave. Any *relative weight cutoff*
  `ε > 1/N` (the standard MPS truncation policy, e.g. `ε = 10⁻¹²`)
  licenses discarding it. In the Fourier route the same information rides
  tail-radix couplings of angle `2π/N ≈ 7·10⁻¹³`, at the edge of what
  accumulated `f64` roundoff can represent coherently. The measured
  consequence: the Fourier-composed adder over `Z_{8.6·10¹²}` is
  machine-accurate on typical inputs, has Hilbert–Schmidt fidelity
  `1 − 10⁻⁹` with the exact adder — and maps `|N−1⟩` wrongly. *Average-case
  operator fidelity is not worst-case operator fidelity*: HS distance is an
  average over columns, and one wrong column costs only `~1/N` of it.

* **As rank, it is free.** The hand-built carry MPO
  (`radix::adder_mpo_exact`) holds the same operator at bond dimension 2 —
  one classical carry bit riding the bond, 0/1 permutation tensors, exact
  at any ring size. The deep-carry branch costs *no extra rank*.

Hence the rule stamped into the code (`adder_mpo_exact`,
`Transducer::to_mpo_with`, both of which canonicalize under an *exact*
policy): **operators with semantically critical low-weight branches must be
compressed by rank, never by weight.** Two independent budgets govern
operator composition — entanglement (bond dimension) and spectral weight
(floating-point precision) — and they fail in different places: rank fails
on scrambling, weight fails on rare-but-exact arithmetic branches at
`N ≳ 1/ε`.


## 8. Stepwise cascades: machines on the bond

**Definition 8.1 (transducer, `cascade`).** A *transducer* is a classical
machine that sweeps the chain once in a fixed direction, carrying a
*message* `m ∈ [0, M)` on the moving front; at site `i` it applies a
deterministic local rule

```text
(message_in, digit_in) → (message_out, digit_out).
```

Lifted to an MPO, the message rides the virtual bond: the site tensor is
the 0/1 tensor of the rule, and boundary vectors at the entry and exit ends
implement a *boundary condition* `Boundary ∈ {Fixed(s), SumAll}` (a basis
vector or the all-ones vector on the message bond).

**Proposition 8.2 (lift).** (i) The lifted MPO has bond dimension `≤ M`
at every interior cut. (ii) With entry `Fixed(s₀)` and exit `SumAll`
("start in `s₀`, drop the final state"), the operator is the graph
`Σ_x |f(x)⟩⟨x|` of the machine's global function `f`, and it is unitary iff
`f` is a bijection. (iii) *(Measured refinement)* recompression generally
finds bond dimensions strictly below `M` — the message alphabet is only an
upper bound on the true width (Theorem 8.5 computes the true value).

*Proof.* (i) by construction; (ii) determinism means each input basis state
selects exactly one message path, so each column of the operator has a
single 1; a 0/1 matrix with one 1 per column is `Σ|f(x)⟩⟨x|`, unitary iff
`f` is injective. ∎

This is the point where the crate's slogan becomes a definition: **bond
dimension is the width of the classical information front.** (The reading
of low-bond MPOs as weighted finite automata is classical [8]; the cascade
layer runs the correspondence in reverse, *designing* the automaton and
inheriting the operator.)

Two arithmetic members generate everything else:

* **Adder** (`m = 2`): rule `v = digit + c_i + carry`, output
  `(v div d, v mod d)` — the bond-2 carry MPO of §7.4, rebuilt as a machine.
* **Multiplier ×k** (`m = k`): rule `v = k·digit + carry`, output
  `(v div d, v mod d)`. Here `(digit, carry) ↔ v ↔ (carry', digit')` is a
  *bijection* `[0,d)×[0,k) ↔ [0,k)×[0,d)` at every site (both sides
  enumerate `[0, kd)`), so the cascade is exactly the schoolbook
  multiplication automaton; dropping the exiting carry is reduction
  mod `N`, and the global map is `x ↦ kx mod N`.

**Proposition 8.3 (the gcd obstruction, quantitatively).** Let
`g = gcd(k, N)` and `A` the mult-cascade MPO. Then
`unitarity_defect(A) := 1 − hs_fidelity(A†A, I) = 1 − 1/g`.

*Proof.* `A = Σ_x |kx mod N⟩⟨x|`, so `A†A = Σ_{x,x'} [kx ≡ kx'] |x'⟩⟨x|` —
the block-all-ones matrix over the fibers of `x ↦ kx`. Each nonempty fiber
has exactly `g` elements (the kernel of ×k on `Z_N` has order `g`), and
there are `N/g` fibers. Then `tr(A†A) = N` (diagonal pairs), and
`‖A†A‖_F² = Σ_fibers g² = (N/g)·g² = Ng`. So
`hs_fidelity(A†A, I) = N² / (Ng · N) = 1/g`. ∎

Measured: `×6` on `Z_2880` (`g = 6`) has defect `0.833 = 1 − 1/6`, with the
explicit collision `|0⟩, |480⟩ → |0⟩` (finding 17). Number theory is not an
analogy here — it *is* an operator property, with the gcd read off a
Frobenius norm.

**Composition (measured, finding 15).** Cascades compose recursively —
`×k₁ ∘ ×k₂ = ×(k₁k₂ mod N)` with composed and directly-built operators
matching bond-for-bond, and adders braid with multipliers by the affine
relation `×k ∘ (+c) = (+kc) ∘ ×k`. The cascade family therefore realizes
the affine group `{x ↦ ax + b : gcd(a, N) = 1}` of the ring in compact
presentations, and modular exponentiation (the Shor kernel [9]) is iterated
cascade composition, with the *measured true operator width* the honest
cost meter.

### 8.4 The true width of a streamed permutation

The message bound `M` and the geometric budget `B_m` are both upper bounds.
The intrinsic width has an exact formula.

**Theorem 8.5 (cut rank).** Fix a cut splitting the chain into left/right
value groups: `x = a·S + b`, `a ∈ Z_L`, `b ∈ Z_S`, `N = L·S`. Let
`gcd(k, N) = 1` and define the *crossing message*

```text
c(b) = ⌊ k·b / S ⌋ ,        Γ = { c(b) mod L : b ∈ [0, S) } .
```

Then the operator Schmidt rank of `M_k : |x⟩ → |kx mod N⟩` across the cut
is exactly `|Γ|`. In particular `width ≤ min(k, L, S)`.

*Proof.* Writing `kb = c(b)·S + (kb mod S)`,

```text
M_k |a, b⟩ = | (ka + c(b)) mod L ,  kb mod S ⟩ ,
```

so `M_k = Σ_{γ∈Γ} (X_L^γ · M_k^{(L)}) ⊗ (M_k^{(S)} · Π_γ)`, where `X_L` is
the cyclic shift on the left group, `M_k^{(·)}` the local multipliers
(bijections since `gcd(k, L) = gcd(k, S) = 1`), and `Π_γ` projects onto
`{b : c(b) ≡ γ (mod L)}`. The left factors are linearly independent
(`Σ_γ λ_γ X_L^γ M_k^{(L)} = 0` implies `Σ λ_γ X_L^γ = 0`, and the powers of
`X_L` have disjoint supports); the right factors are nonzero with disjoint
column supports, hence independent. The rank of `Σ_γ A_γ ⊗ B_γ` with both
families independent is the number of terms. ∎

**Corollary 8.6 (resonance — finding 16).** At the diamond's waist cut
(`L = 24`, `S = 120`), for `k ≡ 1 (mod S)` — i.e. `k = aS + 1` — the
crossing message is `c(b) = a·b`, so the width is `L / gcd(a, L)`. The
repeated-squaring chain `7 → 49 → 2401 → 1921 (mod 2880)` gives:

```text
k        k mod S   crossing message         width
7        7         ⌊7b/120⌋ ∈ {0..6}        7
49       49        ⌊49b/120⌋ ∈ {0..48}      24   (= L: budget-capped)
2401     1         20·b  mod 24             24/gcd(20,24) = 6
1921     1         16·b  mod 24             24/gcd(16,24) = 3
```

— exactly the measured operator widths. **Multiplying by 1921 is a thinner
operator than multiplying by 7.** The cost curve of modular exponentiation
on a chain is number-theoretically resonant, non-monotone in `k`, and
computable in advance from `gcd` data; the operator merely *remembers* it.

The theorem also certifies narrowness where the message bound is
astronomical: for `k' = 823 = 7⁻¹ (mod 2880)`, `823·7 = 48·120 + 1` makes
`c(b)` periodic with period 7 up to multiples of 24, so `|Γ| = 7` — the
inverse multiplier is intrinsically 7 wide at the waist, three orders below
its own message bound. Which raises the question the dual machine answers:
what *machine* attains that width?

### 8.7 The geometrically opposed dual

**Theorem 8.8 (division machine).** Let `gcd(k, N) = 1`. The transducer
sweeping **MSB → LSB** with remainder state `r ∈ [0, k)` and site rule

```text
t = r·d + digit_in ;    digit_out = t div k ;    r' = t mod k
```

(a bijection `[0,k)×[0,d) ↔ [0,d)×[0,k)` via `t ∈ [0, kd)`), lifted with
boundary conditions **entry `SumAll`, exit `Fixed(0)`**, is exactly the
unitary `|x⟩ → |k⁻¹·x mod N⟩`, at message width `k`.

*Proof.* Fix an entry remainder `r₀`. Sweeping most-significant-first, the
machine performs schoolbook long division of the (mixed-radix) number
`r₀·N + x` by `k`: the invariant after processing a prefix of value `X_p`
with place-product `D_p` is `r = (r₀·D_p + X_p) mod k`, with the emitted
digits encoding the partial quotient (each emitted digit is `< d` at its
site since `t < kd`). At exit, the output value is `⌊(r₀N + x)/k⌋` and the
exit remainder is `(r₀N + x) mod k`. The exit condition `Fixed(0)` keeps
the branch iff `k | r₀N + x`; since `gcd(k, N) = 1`, for each `x` exactly
one `r₀ ∈ [0, k)` satisfies `r₀ ≡ −x·N⁻¹ (mod k)`. The surviving branch
outputs `y = (r₀N + x)/k < N` with `k·y ≡ x (mod N)`, i.e.
`y = k⁻¹x mod N`. Formally, each fixed-boundary branch
`T_{r₀} = ⟨0|_exit machine |r₀⟩_entry` is a partial isometry whose domain is
the residue class `{x ≡ −r₀N (mod k)}` and whose image is the band
`[r₀N/k, (r₀+1)N/k)`; the `k` domains and the `k` images each partition
`Z_N`, so `Σ_{r₀} T_{r₀}` is the total bijection. ∎

Three remarks make this the conceptual center of the cascade layer:

* **The quantum boundary is the enabling ingredient — and it is virtual.**
  A classical FSM must *choose* its initial state; this machine *enters in
  superposition over the unknown wrap multiple* (`SumAll` = the all-ones
  boundary vector) and is *postselected on exact division* (`Fixed(0)`).
  No runtime postselection occurs: the assembled operator is exactly
  unitary, because for every input precisely one branch survives the
  contraction. The boundary conditions carve a unitary out of a sum of
  partial isometries — linear algebra, not measurement.
* **The width collapse is exponential.** The forward carry machine for
  `×k⁻¹` needs message width `k⁻¹ mod N` — `823` on the diamond, `≈ 6·10¹²`
  on the 25-site wave. The opposed remainder machine needs `k` — `7`. The
  measured `÷7` bond profile realizes the intrinsic width that Theorem 8.5
  certified (finding 18).
* **Duality is geometric.** Multiplication's carries flow LSB → MSB;
  division's remainders flow MSB → LSB, with the *same* state-set size. The
  two fronts are opposed in direction and inverse in semantics.

**Corollaries (all measured at fidelity 1, finding 18).**

* *Annihilation:* `÷k ∘ ×k = 1` — the composed opposed fronts collapse to
  the all-ones bond profile.
* *Ratio machines:* `×(k₁k₂⁻¹ mod N) = ÷k₂ ∘ ×k₁` — two narrow opposed
  fronts realize a single-front machine of astronomical width
  (`×2357 = ÷11 ∘ ×7`, measured true width `[2,6,17,13,6,2]`).
* *Reflection duality:* the digit complement `x_i ↦ d_i−1−x_i` maps
  `x ↦ N−1−x` (since `Σ_i (d_i−1)P_i = Σ_i (P_{i−1} − P_i) = N−1`,
  telescoping with `P_{−1} := N`), so negation is
  `(+1) ∘ complement` — machine width 2 — and `×(N−k) = neg ∘ ×k` at width
  `max(2, k)` instead of `N−k`.

**The minimal-machine principle.** Operator Schmidt rank across a cut is
intrinsic (Theorem 8.5 computes it for this family); every implementing
machine upper-bounds it by its crossing-message alphabet; and recompression
is the oracle that finds it. Every measured case in this family attains its
intrinsic width through one of the geometric duals — sweep reversal or ring
reflection. **The operator is only as wide as the narrowest machine that
computes it**, and "find the dual machine" is the discrete analogue of
choosing a better factorization route in §9.5.

**Proposition 8.9 (Fourier conjugation reverses the cascade).** For
`gcd(k, N) = 1`, `F M_k F† = M_{k⁻¹}` (the scaling theorem: from
`M_k F |x⟩ = Σ_y ω^{xy}|ky⟩/√N = F|k⁻¹x⟩`). In the crate's reversed-digit
frame `V = ρ∘F`, conjugation transports the multiplier to
`ρ M_{k⁻¹} ρ` — the inverse multiplier *acting on the reversed encoding*,
i.e. with its carry sweep running in the opposite direction. Verified on a
non-palindromic chain (finding 17): the Fourier frame implements, in one
stroke, both the arithmetic inversion and the geometric front reversal that
§8.7 realizes combinatorially.

### 8.10 The executable atlas

Theorem 8.5 is constructive: `|{⌊kb/S⌋ mod L}|` is elementary arithmetic,
so the width of a multiplication operator is known *before the operator
exists*. The `width` module implements the calculus (closed form
`min(k, L)` when `k ≤ S` — the crossing sequence steps by 0 or 1 and is
onto `[0, k)` — enumeration otherwise, with `k` always reducible mod `N`),
and the test suite pins it against exactly-recompressed cascade MPOs
bond-for-bond, including through composition: `×49 = ×7 ∘ ×7` on
`diamond(2,4)` lands exactly on the predicted `[2,3,3,2]` — *narrower than
its own factor* `×7` at `[2,6,6,2]`.

Measured atlas of `Z_2880` (finding 20, `examples/width_atlas.rs`):
668 of the 768 invertible multipliers saturate the waist cap 24; the
resonant tail reaches width 2 (`×1441` — an involution, `1441² ≡ 1` — and
`×2879 = −1`) and width 3 (`×961`, `×1439`, `×1921`). On the wave
(`N ≈ 8.6·10¹²`) the squaring orbit of 7 has central-cut widths
`7 → 49 → 2401 → 2 073 600 (cap) → 1 077 111 → 1 928 737`: resonance is
not a small-ring artifact, and it is computed a priori.

The atlas also sharpens §8.7's minimal-machine principle into a *negative*
result about construction. Let the best single-front alphabet be
`B(k) = min(k, k⁻¹, N−k, (N−k)⁻¹)` — the smallest message set any of the
four dual sweeps (forward/opposed, direct/reflected) must carry. Measured:
**584 of 768 multipliers have intrinsic waist width below `B(k)` by a
factor > 8** (extreme: `×1441`, width 2 against `B = 1439` — a 719× gap).
Theorem 8.5's `Γ` *is* the forward machine's realized crossing set, so
per-cut recompression always attains the intrinsic rank; but no single
sweep's **alphabet** comes close — the narrow presentation exists only
after composition + recompression. Machines are upper-bound certificates
for construction, and the intrinsic width is in general strictly finer
than every one of them (see the revised open problem 2 in §13).

### 8.11 Controlled cascades: the priced Shor kernel

Phase estimation needs controlled powers
`C-U^{2^j} = |0⟩⟨0| ⊗ I + |1⟩⟨1| ⊗ U^{2^j}`. On the chain these are
boundary-machine algebra: `Transducer::mult_skipping` builds
`I_{d_c} ⊗ (×k mod N/d_c)` — the carry front *tunnels through* the
control site, and the lifted operator is a pure tensor product across it
(bond 1) — and `Mpo::select_on` composes projector branches into the
controlled operator. Its exact cost:

**Proposition 8.11 (controlled cut rank).** Let
`CM = P₀ ⊗ I + P₁ ⊗ M_k` with the control left of a register cut
`(L | S)`, and let `w` be the width of `M_k` at that cut (Theorem 8.5).
Then

```text
  χ_cut(CM) = w + 1   if k ≢ 1 (mod S) ,
  χ_cut(CM) = w       if k ≡ 1 (mod S)    — the control is free,
```

and across the control's own cut `χ = 2` (for `M_k ≠ I`).

*Proof.* Decompose `M_k = Σ_γ (X_L^γ M_k^{(L)}) ⊗ (M_k^{(S)} Π_γ)` as in
Theorem 8.5. The left family `{P₀⊗I_L} ∪ {P₁⊗X^γ M^{(L)}}` is always
independent (the projectors have disjoint supports), so the rank is
`1 + w` unless the right family `{I_S} ∪ {M^{(S)}Π_γ}` is dependent.
`I_S = Σ_γ c_γ M^{(S)}Π_γ` forces `Σ_γ c_γ Π_γ = (M^{(S)})⁻¹`, whose left
side is diagonal — possible iff `M^{(S)} = I`, i.e. `k ≡ 1 (mod S)`, where
`c_γ = 1` works (`Σ_γ Π_γ = I`). In that case
`CM = Σ_γ (P₀⊗I_L + P₁⊗X^γ M^{(L)}) ⊗ Π_γ` is a `w`-term decomposition
with both families independent. ∎

Measured, the formula prices the entire order-finding protocol on the
diamond (finding 23, `examples/shor_kernel.rs`): the `C-U^{2^j}` family
for `k = 7` on `Z_1440` has bond profiles `[2,4,8,8,6,2]`,
`[2,4,13,24,6,2]`, `[2,3,3,3,3,2]`, `[2,3,3,3,3,2]` — every bond as
predicted, and the deep powers `961 = 7⁴` and `481 = 7⁸` satisfy
`k ≡ 1 (mod S)` on *every* cut (`961 ≡ 1` mod 480, 120, 24, 6, 2), so
their control is free: the controlled operator is exactly as wide as the
bare multiplier. Resonance and control-freeness are the same phenomenon
seen from two sides — deep in a squaring orbit, the multiplier acts
almost-locally, and there is nothing for the control to correlate with.

Two further measured facts complete the kernel:

* **Phase kickback is bond-free.** On an eigenstate `|u_s⟩` of `×k`, the
  controlled operator sends `|+⟩⊗|u_s⟩` to a *product*
  `((|0⟩ + e^{2πi·s·2^j/r}|1⟩)/√2) ⊗ |u_s⟩`: the control cut stays at
  `χ = 1` while the phase lands (measured phases 30°/60°/120°/240° exact
  at tolerance `10⁻⁶`; a non-eigenstate register gives `χ = 2`). The
  resource statement of phase estimation, read directly off a bond
  dimension.
* **Order finding closes end-to-end.** Semiclassical (Kitaev) phase
  estimation with one reused valley qubit — measurement by projection,
  feedback rotations conditioned on earlier bits — plus continued
  fractions bounded by the Carmichael exponent `λ(1440) = 24` recovers
  `ord(7 mod 1440) = 12` from six 10-bit runs (measured fractions
  `1/3, 5/12, 1/4, 1/3, 2/3, 5/6`, lcm of denominators 12). The register
  is never measured; each run collapses onto an eigenstate through the
  control's backaction alone. Both Shor registers live on one dimension
  wave: control in a valley, arithmetic in the interior.


## 9. Operator flows and the operation-width cursor

Is an operation one object, or a divisible extent? The `flow` layer makes
the question quantitative by extending blocks to one-parameter groups.

**Definition 9.1 (`FourierFlow`).** For the adder family, define

```text
U^t := V† · D(t·c) · V ,      D(s) = diag_y( e^{2πi·s·y/N} )
```

(the real-`c` ramp of §7 makes `D(s)` exact for all real `s`). Since
`D(s)D(t) = D(s+t)` exactly, `t ↦ U^t` is an exact one-parameter group,
`U^0 = I`, `U^1 = A_c`; the group law and endpoints are verified at MPO
precision (`U^{0.3} ∘ U^{0.7} = U^1` at `10⁻⁸`). `span(t₀,t₁) = U^{t₁−t₀}`
— the flow is autonomous.

A branch choice is being made here: `U^t = exp(it·c·H)` with
`H = V†(2π ŷ/N)V` — the logarithm generated by the frame's number ramp.
Other logarithms of `A_c` exist (rewind any phase winding by `2π`); this
one is canonical for the frame and linear in `t`, which is what makes every
snapshot cheap (two compositions through a fixed compiled `V`).

**Proposition 9.2 (snapshot kernel).** `U^s` (shift extent `s = t·c`) is
the circulant

```text
U^s = Σ_{j ∈ Z_N} w_j(s) · X^j ,     w_j(s) = (1/N) Σ_y e^{2πi(s−j)y/N}
                                            = (1/N) · (e^{2πi(s−j)} − 1)/(e^{2πi(s−j)/N} − 1),
```

the periodic Dirichlet (sinc) kernel centred at `j = s`. At integer `s` it
degenerates to the single shift `X^s`; between integers it is a genuinely
delocalized superposition of shifts with `1/|j−s|` tails.

*Proof.* `U^s = F†D(s)F` is diagonalized by the Fourier basis with symbol
`e^{2πisy/N}`; its circulant coefficients are the inverse DFT of the
symbol — a geometric sum. ∎

This is the mechanism behind **finding 12 — width sees the integers**: at
integer `s` the operator is a permutation whose exact presentation is the
bond-2 carry MPO (measured `χ = 3`, the exact form plus compile residue at
the working cutoff); at fractional `s` it is a coherent mixture of many
shifts, and the measured width sits at `χ = 9–11` with operator
entanglement peaking mid-flow. *Bond dimension along a flow is an
arithmetic quantization detector*: `χ(U^t)` dips precisely where the
operation extent is whole.

**Proposition 9.3 (first-moment law — the lattice wobble).** Let
`|ψ_t⟩ = U^t|x₀⟩` and `⟨ω^x⟩_t := Σ_x |ψ_t(x)|² ω^x` the circular first
moment. Then, exactly,

```text
⟨ω^x⟩_t = ω^{x₀ + s} · ( (N−1) + e^{−2πis} ) / N ,        s = t·c .
```

*Proof.* Let `W = Σ_x ω^x |x⟩⟨x|` (clock) and `φ = F|x₀⟩`, i.e.
`φ(y) = ω^{x₀y}/√N`. A direct computation gives `F W F† = S†` (the lowering
shift on frequency space), and conjugating by the ramp,

```text
D(s)† S† D(s) = e^{2πis/N} ( S† + (e^{−2πis} − 1)|N−1⟩⟨0| ) :
```

the ramp phases cancel to a constant `e^{2πis/N}` on every step
`y → y−1` *except at the wrap* `0 → N−1`, which acquires the anomalous
factor `e^{−2πis}`. Sandwiching in `φ` (where `⟨φ|S†|φ⟩ = ω^{x₀}` and
`⟨φ|N−1⟩⟨0|φ⟩ = ω^{x₀}/N`) yields the formula. ∎

Consequences, all measured (finding 14): the extracted ring position is

```text
pos(t) = x₀ + s − sin(2πs)/2π + O(1/N) ,
```

— **linear drift plus a lattice wobble** `−sin(2πs)/2π`, of amplitude
`1/2π` ring units *independent of `N`*: the ring's discreteness pushes back
on the continuous flow through exactly one rank-one wrap term, vanishing at
integer extents. Exact Weyl covariance would be wobble-free; the wobble
*is* the failure of the ramp to be a character of the cyclic group at
fractional `s`.

**Proposition 9.4 (participation breathing).** With `τ` the fractional part
of `s`, the participation `P(t) = Σ_x |ψ_t(x)|⁴` satisfies

```text
P(t)  →  1 − (2/3)·sin²(πτ)        (N → ∞),
```

so quarter-, half-, and whole-extent snapshots breathe as
`1 → 2/3 → 1/3 → 2/3 → 1` — the measured values.

*Proof sketch.* `|ψ_t(x)|² = sin²(πτ) / (N² sin²(πδ/N))` with
`δ = x₀+s−x ≡ τ (mod 1)`. Summing fourth powers and using
`Σ_{m∈Z} (m+τ)⁻⁴ = π⁴(3 − 2sin²πτ)/(3 sin⁴πτ)` (differentiate the
cotangent series twice) gives the limit. ∎

### 9.5 Imputation has geometry

Splitting a coarse block into finer ones is not unique. Two measured routes
factor the unit block `A = U^1` into halves (finding 13):

* **geodesic**: `A = U^{1/2} · U^{1/2}` — midpoint on the flow, measured
  `χ = 11` each (a delocalized Dirichlet kernel);
* **causal**: `A = W₂ · W₁` with `W₁ = D(½)V`, `W₂ = V†D(½)` — detour
  through the Fourier frame, measured `χ = 8` each (the frame's own core
  width: multiplying by a diagonal cannot change bond dimension, so
  `χ(W₁) = χ(V)`).

Both products reproduce `A` at fidelity `1 − 10⁻⁹`, yet the two midpoints
are measured as Hilbert–Schmidt orthogonal to six decimals. The precise
statement is prettier than the idealized one:

**Proposition 9.6 (Gauss-sum orthogonality of midpoints).** The normalized
overlap of the two half-way operators is exactly

```text
hs_fidelity( U^{1/2} , W₁ )  =  |tr V|² / N² ,
```

independent of the ramp. For the plain DFT frame `tr F` is the quadratic
Gauss sum (`|tr F|² = 2` for `N ≡ 0 mod 4`); for the reversed-digit frame
`tr V = N^{−1/2} Σ_x ω^{x·ρ(x)}` is a twisted Gauss-type sum, measured
`|tr V|² ≈ 0.0112` on `Z_2880`. Either way the overlap is `1/N²`-scale —
the identity is verified numerically to six significant figures on
`Z_2880`, where both sides equal `1.351·10⁻⁹`, printing as `0.000000`:
the midpoints are *asymptotically*, not identically, orthogonal.

*Proof.* `tr(U^{1/2†} W₁) = tr(V† D(−½) V · D(½) V)`. Cycling the trace and
using `VV† = I` collapses it to `tr(D(−½) V D(½))`; the two diagonal ramps
cancel pointwise on the diagonal, leaving `tr V`. Normalize by
`‖U^{1/2}‖_F ‖W₁‖_F = N`. ∎

The moral survives the refinement: decomposing a coarse operation into
finer ones is a **choice of path through operator space** — routes with
identical endpoints pass through essentially disjoint midpoints, at
route-dependent width cost (here the causal detour is strictly narrower
than the geodesic). "Width-optimal imputation" is thereby a well-posed
optimization problem, posed concretely by this library (§13).

### 9.7 The width cursor

The `WidthCursor` holds a pipeline of flow segments `[t_0,t_1), [t_1,t_2), …`
each as a collapsed MPO. `refine(i)` splits a segment at its midpoint,
*imputing* both halves from the flow (consistent with the parent by the
group law); `coarsen(i)` composes neighbours back into a wider block. The
total composition is invariant under both — measured at fidelity
`1.000000000` through four dyadic zoom levels (finding 14) — so the cursor
maintains **adaptive temporal resolution**: a coarse past, an increasingly
fine present, with the semantics pinned. `trajectory` then walks a state
through the pipeline at its current grain; Propositions 9.3–9.4 are the
laws that walk obeys, with exact relocalization at whole extent — the
cursor's finest observable grain agrees with arithmetic exactly when the
operation width is whole.

### 8.12 Meet-in-the-middle: retrodictive-prediction and the tail-radix solve

The minimal-machine principle (§8.7) says an operator is only as wide as
its narrowest machine, and §8.7's dual machines realize it. Read
geometrically, that dual is a *crossing of opposite time directions*, and
the tail-radix frame (§7) is where the two directions meet — the synthesis
measured in finding 27.

**Proposition 8.12 (the dual-direction crossing).** For
`gcd(a, N) = gcd(b, N) = 1`, the composite

```text
  ÷b ∘ ×a  =  ×(a·b⁻¹ mod N)
```

crosses a *forward* carry front (`×a`, sweeping LSB→MSB — **prediction**)
with a *backward* remainder front (`÷b`, sweeping MSB→LSB —
**retrodiction**, §8.7). Its operator width is bounded by
`max(width(×a), width(÷b))`, whereas the naive forward machine for
`m = a·b⁻¹` carries a multiplication carry in `[0, m)` and needs message
dimension `m`.

*Proof.* Composition of the two opposed cascades; the width bound is
Theorem 8.5 applied to each narrow factor, and the product's cut rank is at
most the product of the factors' — here dominated by the wider of the two
narrow fronts. ∎

This is the geometric content of **retrodictive-prediction**: a forward
prediction reconciled with a backward retrodiction at the crossing,
computing an operation neither front reaches alone within the width budget.
Measured (finding 27): on `Z_2880`, `×2357 = ×7 ⋈ ÷11` at width 17,
`×2659 = ×7 ⋈ ÷13` at 19, `×2095 = ×5 ⋈ ÷11` at 15 — each verified, and
each *unbuildable* as a single forward cascade (`m > 512`, the direct-
construction cap). The crossing is not merely narrower; it is the only
feasible route.

**Proposition 8.13 (the tail-radix frame is the meeting point).** The
reversed-digit QFT `V` (§7) conjugates a forward multiplier into a backward
one, `V M_k V† = ρ M_{k⁻¹} ρ` (Prop. 8.9): prediction and retrodiction are
the *same operator seen from the two sides of the frame* (measured hs
fidelity 1.0). And `V` diagonalizes the additive family (the shift
theorem, §7), so a **mesh of additive solve-components commutes in the
frame** — its ideal ordering is trivial: enter the frame once, sum the
diagonal ramps, leave once. The tail-radix banding (Cor. 7.2) makes that
frame narrow — `χ = 6` reversed against `χ = 36` with the digit reversal,
which **cancels meet-in-the-middle** between `V` and `V†` (§7.3, finding 9).

The **global tail-radix solve** (finding 27) is then: to compute a wide
affine/modular operation, cross a forward prediction with a backward
retrodiction for the multiplicative part (Prop. 8.12), sum the additive
components in the tail-radix frame (Prop. 8.13), and order the whole solve
by the ring's own banded Fourier web — reading the answer at the meeting
point in the middle. Measured (`examples/meet_in_the_middle.rs`): a mesh of
5 additive components collapses from 5 frame round-trips (~4 ms) to one
(~0.7 ms), agreeing exactly; and the full affine map
`x → 2357·x + 500 mod 2880` — whose single-front multiplier is unbuildable
— solves as a width-17 mesh, verified on basis states. The efficiency is
the minimal-machine width of §8.7 realized *geometrically*, plus the
frame's decoupling of the additive mesh; it is a width/representation
efficiency (bond dimension), not a quantum speedup.

**Proposition 8.14 (bond dimension is a two-sided resource).** Let `U` have
operator Schmidt rank `χ` across a cut. Then (i) simulating `U` costs
`O(n·χ³·d³)` time and `O(n·χ²·d²)` space (§12), and (ii) any circuit
implementing `U` must send `Ω(log₂ χ)` two-qudit gates across that cut —
each gate raises operator entanglement by `O(1)`, so `log₂ χ` bits of it
require that many crossing gates. `χ` therefore prices *both* simulation and
implementation, and the constructions of §8.12 minimize it on both sides at
once.

The measured consequences (finding 28, `examples/resource_efficiency.rs`),
stated with their honest scope:

* **Ancilla.** A single-front `×m` cascade carries a `⌈log₂ m⌉`-qubit carry
  register to realize an operator whose true operator entanglement is only
  `log₂(width)`. On `Z_2880`, `×2357` spends 12 carry qubits for a
  `3.9`-bit operator; the crossing `×7 ⋈ ÷11` (Prop. 8.12) achieves the
  intrinsic width at 4 — the machine no longer carries more than the
  operator contains. (Real, but a statement about the cascade realization,
  not a bound over all circuits.)
* **Gate count.** The reversal-free ordering (Prop. 8.13) drops the
  digit-reversal swap network (`Θ(n)` gates), and additive components batch
  in the frame: `K` additions cost one QFT round-trip plus `K` cheap phase
  layers rather than `K` round-trips — a measured `3.1×…7.7×` for
  `K = 4…32`, approaching `(2·|QFT|+|ramp|)/|ramp| ≈ 10×`. This *is* Draper's
  Fourier arithmetic [7]; the only novelty is that it runs on the
  heterogeneous mixed-radix wave and composes with the crossing.
* **What is not claimed.** No asymptotic advantage over the best known
  modular-arithmetic circuits. Indeed the framing cuts the other way:
  representation efficiency *is* the quantum-advantage boundary — a
  bounded-`χ` computation is classically simulable (§12) and so carries no
  quantum advantage. These structured computations sit on the simulable side
  by construction; the value is the unified meter `χ` and a geometry that
  minimizes it for simulation and implementation together.

### 9.8 The frame flow: fractional powers of the QFT

The first edition of this document (§13, problem 3) supposed fractional
powers of the QFT needed `V`'s eigenframe. For the standard-order frame
they need no eigensolver at all:

**Proposition 9.7 (projector flow of the QFT — finding 22).** On any ring
`Z_N`, `F² = Π` (the reflection `x ↦ −x mod N`) and `F⁴ = I`. Hence the
spectral projectors of `F` are polynomials in `F`,
`P_j = ¼ Σ_{m<4} i^{−jm} F^m`, and the matrix-power fractional Fourier
transform [12]

```text
  F^t := Σ_j i^{j·t} P_j = Σ_{m=0}^{3} c_m(t) · F^m ,
  c_m(t) = ¼ Σ_{j=0}^{3} e^{iπ·j·(t−m)/2} ,
```

is an exact one-parameter group of period 4, unitary at every `t`,
assembled from four MPOs the crate already owns — `I` (χ = 1), `F`,
`Π` (the width-2 reflection machine of §8.7), and `F†` — by operator
linear combination (`flow::FrameFlow`).

*Proof.* `F²|x⟩ = (1/N) Σ_{y,z} ω^{xz+zy} |y⟩ = Σ_y δ_{x+y≡0} |y⟩ =
|−x mod N⟩`, so `F⁴ = Π² = I`. The projector algebra
(`P_j P_k = δ_{jk} P_j`, `Σ_j P_j = I`) reduces to `F⁴ = I`; the group law
and unitarity follow from the eigenvalue reading `F^t = Σ_j i^{jt} P_j`. ∎

Measured on `diamond(2,4)` (`examples/fractional_fourier.rs`): the group
law and endpoint identities hold at fidelity `1.000000000`
(`F^½ ∘ F^½ = F`, `F^0.7 ∘ F^1.3 = Π`, `F^3.5 ∘ F^0.5 = I`), with unit
norm on states at every `t`. The width signature refines finding 12
instructively. Operator entanglement is *pinned* to the pure-power value
at every integer — the curve passes through `0 → 5.170 → 1.000 → 5.170 →
0` bits, an extremum each time — but bond dimension collapses only at
`t ≡ 0, 2 (mod 4)` (χ = 1 and 2 against a flat 36 elsewhere), because `F`
itself *saturates* the profile's operator-width cap: at odd integers,
rank has nothing to collapse to. **Width sees the frame exactly where the
frame is narrower than the geometry** — quantization detection by bond
dimension requires headroom between the operator and the geometric
budget. Participation of `F^t|x₀⟩` breathes with period 2 (localized at
even `t`, maximally flat at odd `t`), the frame-flow counterpart of
Proposition 9.4.


## 10. Self-stabilizing boundary systems

Sections 7–9 kept the message boundaries *open* (fixed entry, dropped or
postselected exit). Closing the loop — feeding the exiting message back
into the entry — is the third boundary regime, and it changes the category:
the operator becomes non-unitary, and non-unitarity is precisely the
resource an attractor needs (a unitary map preserves distances; nothing
contracts, so nothing stabilizes).

**Definition 10.1.** For a transducer `T` with message dimension `M`, the
*looped* operator is `T_loop = Σ_{m<M} ⟨m|T|m⟩` — the partial trace over
the message bond (`to_mpo_looped`).

**Theorem 10.2 (closing the loop changes the arithmetic).** The looped
adder is *end-around carry* — ones'-complement addition — computing
`mod N−1`:

* for `x + c ≤ N−2`: exactly the branch `m = 0` is consistent
  (`m_out = 0`), giving `x + c`;
* for `N ≤ x + c ≤ 2N−3`: exactly the branch `m = 1` is consistent,
  giving `x + c + 1 − N = (x+c) mod (N−1)`;
* for `x + c = N−1` (and the wrap image `2N−2`): **both** branches are
  consistent, and the column contains the two representatives of zero —
  `|N−1⟩` and `|0⟩` — in equal superposition (the *seam*).

Likewise the looped multiplier computes `×k mod (N−1)` (from
`kx + m = mN + r ⇒ r ≡ kx (mod N−1)`, with a unique consistent `m` except
on the seam). On `diamond(2,5)` this replaces the highly composite ring
`Z_2880` by the **prime field** `Z_2879`.

*Proof.* Case analysis on the carry: with entry carry `m ∈ {0,1}` the open
adder computes `x + c + m` with exit carry `m' = ⌊(x+c+m)/N⌋`; the loop
keeps branches with `m = m'`. If `x+c ≤ N−2`, `m = 0` gives no overflow
(`m' = 0` ✓) while `m = 1` gives `x+c+1 ≤ N−1 < N` (`m' = 0 ≠ 1` ✗). If
`x+c ≥ N`, `m = 0` overflows (`m' = 1 ≠ 0` ✗) while `m = 1` yields
`x+c+1−N ∈ [1, N−1]` with `m' = 1` ✓. At `x+c = N−1` both checks pass:
`m = 0` gives `|N−1⟩` without overflow, `m = 1` gives `x+c+1 = N ⇒ |0⟩`
with carry out 1. ∎

This is a piece of classical computer arithmetic — the end-around-carry
adder of ones'-complement machines — resurfacing as the trace of a quantum
boundary, seam and all: ones' complement famously has "negative zero", and
here it appears as an *operator defect* rather than an engineering nuisance.

**Corollary 10.3 (the boundary heals the operator — finding 19b).** `×6` on
`Z_2880` has `gcd = 6`, unitarity defect `1 − 1/6 = 0.833` (Prop. 8.3). The
*same transducer looped* acts on `Z_2879`; since 2879 is prime,
`gcd(6, 2879) = 1` and the traced operator is exactly unitary (measured
defect `10⁻¹⁶`) away from the seam. No property of the local rules changed
— the boundary alone selected a ring in which the obstruction does not
exist.

**Proposition 10.4 (the seam is a Jordan block — finding 19c).** The looped
`+0` machine is exactly

```text
M = I + |0⟩⟨N−1| ,
```

a defective operator: eigenvalue 1 with a rank-one nilpotent part. Then
`M^k = I + k·|0⟩⟨N−1|`, so iterating from the negative zero,

```text
M^k |N−1⟩ ∝ |N−1⟩ + k|0⟩ ,      1 − |⟨0|ψ_k⟩|² = 1/(k²+1)  exactly,
```

— polynomial self-stabilization of the number representation, with the
measured error matching `1/(k²+1)` digit for digit. (Ordinary states are
already fixed points; only the seam flows.) The normalized iteration
`ψ → Mψ/‖Mψ‖` is power iteration, and a defective dominant eigenvalue is
precisely the case where power iteration converges polynomially. A detail
that explains the measured run lengths: the driver's stopping rule watches
the *successive* overlap, and `1 − |⟨ψ_k|ψ_{k+1}⟩|² ≈ k⁻⁴` for this flow —
so tolerance `10⁻¹²` halts near `k ≈ 10³` (measured: 1001) while the true
error is then only `~10⁻⁶`.

**Proposition 10.5 (the convergence law is designable — finding 19d).**
Compose a *damped seam* — the bond-1 correction
`I − (1−γ)|N−1⟩⟨N−1|`, built from `identity + basis_transfer·scale` — with
the pump:

```text
M_γ = (I + |0⟩⟨N−1|)(I − (1−γ)|N−1⟩⟨N−1|) = I − (1−γ)|N−1⟩⟨N−1| + γ|0⟩⟨N−1| .
```

On the seam plane `span{|N−1⟩, |0⟩}` the spectrum is `{γ, 1}` — the Jordan
block has been split into genuine eigenvalues — and in closed form

```text
M_γ^k |N−1⟩ = γ^k |N−1⟩ + γ(1−γ^k)/(1−γ) |0⟩ ,
error ≈ ((1−γ)/γ)² · γ^{2k}      (exponential at amplitude rate γ).
```

Predicted steps to `10⁻¹²`: `≈ 20 / 10 / 6` at `γ = 0.5 / 0.2 / 0.05`;
measured: `20 / 11 / 7` (the stopping rule detects settling one step after
it happens). Self-stabilization is thus **engineered at the boundary, not
inside the operator**: the loop chooses the ring (composite → prime,
dissolving gcd obstructions), the seam supplies the attractor that repairs
the representation, and a rank-one damper sets the convergence law —
polynomial by default, exponential at a chosen rate.

### 10.6 Twisted loops: the mod N±1 family

The loop of Definition 10.1 admits a twist: feed the exit back through a
permutation `σ` of the message set, `T_σ = Σ_m ⟨σ(m)| T |m⟩`
(`Transducer::to_mpo_looped_twisted`; Definition 10.1 is `σ = id`). The
reversal twist completes a ring family around the open machine:

**Theorem 10.7 (reversal twist = diminished-one mod N+1 — finding 21).**
For the carry machines, twisting the loop by the message reversal
`σ(m) = M−1−m` yields exactly the **diminished-one arithmetic** of
`Z_{N+1}` — the encoding of Fermat-number-transform hardware [13] — in
which chain value `x` represents `v = x+1 ∈ [1, N]`:

* the twisted adder `+c` maps `v ↦ v + (c+1) mod (N+1)`, with the branch
  **annihilated** when the result is the unrepresentable zero (a *hole*
  at `x = N−1−c`);
* the twisted multiplier `×k` maps `v ↦ k·v mod (N+1)` exactly, hole-free
  and unitary iff `gcd(k, N+1) = 1`.

*Proof.* Adder: with entry carry `m ∈ {0,1}` the open machine computes
`x+c+m` with exit carry `m′ = ⌊(x+c+m)/N⌋`; the twist keeps branches with
`m′ = 1−m`. For `x+c ≤ N−2` only `m = 1` is consistent (no overflow),
giving `x+c+1`; for `x+c ≥ N` only `m = 0` is (overflow), giving `x+c−N`;
at `x+c = N−1` neither is — the hole. In value coordinates the two live
branches are `v ↦ v+(c+1)` and `v ↦ v+(c+1)−(N+1)`. Multiplier: entry
`m ∈ [0,k)` computes `kx+m` with exit `m′ = ⌊(kx+m)/N⌋`; imposing
`m′ = k−1−m` and reducing mod `N+1` (where `N ≡ −1`) forces
`out ≡ kx+k−1`, i.e. `v_out ≡ k·v (mod N+1)`; the linear relation
`m(N+1) = (k−1)N + out − kx` then has exactly one solution `m ∈ [0, k)`
with `out ∈ [0, N)` for each `x` — except when `k·v ≡ 0 (mod N+1)`, the
holes. ∎

The family is symmetric around the open ring, and the two closures fail
in dual ways: **mod N−1 has a double zero** (the seam — Proposition
10.4's Jordan pump), **mod N+1 has a missing zero** (the hole — its
traced identity is the unilateral shift `Σ_{x<N−1} |x+1⟩⟨x|`, a *drain*
under which the uniform state obeys the exact law
`‖T_σ^k · u‖² = (N−k)/N`, measured to six digits). The same input
`x = N−1−c` hits both defects. Unitarity on each closed ring is governed
by gcd against *that* ring: on `[2,3,4,3]` — `N = 72`, flanked by the
twin primes 71 and 73 — both closures are fields and heal every open
obstruction (`×6`: defect `0.833 → ≤ 10⁻¹⁶` both ways); `×5` is unitary
open but breaks on the twist ring `25 = 5²`; and the diamond's twist ring
`2881 = 43·67` exhibits its factorization operationally as collisions
plus holes, branch by branch (`examples/boundary_twists.rs`). Primality
of `N ∓ 1` is a design criterion *selectable by profile*.

The taxonomy of boundary conditions, assembled:

```text
entry        exit         semantics
Fixed(0)     SumAll       mod-N arithmetic (drop the wrap)          §8
SumAll       Fixed(0)     postselected exact division = ×k⁻¹        §8.7
looped (σ = id)           mod N−1 arithmetic + seam attractor       §10.2–5
looped (σ = reversal)     mod N+1 diminished-one + hole drain       §10.6
```

One machine body, an operator category per boundary — unitary arithmetic
in three different rings, unitary inverse arithmetic by dual sweep, and
non-unitary dynamical systems with designed defects (a seam that repairs,
a hole that drains) — all selected by vectors, and one permutation, on a
bond of dimension `k`.


## 11. Crossing strands: dimension-wave networks

Every structure so far lives on one chain. The `crossing` module takes the
first step past that — the README's "beyond chains" direction — with *two*
dimension-wave strands sharing one MPS, crossed so they interact across
scales. The canonical object is two peaks `[1,2,3,4,5,4,3,2,1]` laid
antiparallel into an X.

**Definition 11.1 (crossing).** Given strands `A`, `B` of equal length `n`,
the *antiparallel* pairing couples `A[i]` to `B[n−1−i]` (the X); the
*parallel* pairing couples `A[i]` to `B[i]` (a ladder). A pairing is
*dimension-consistent* where `dim A[i] = dim B[partner(i)]`; for twin
palindromic strands every pair matches. The *interaction level* selects
which dimension band `d` carries a coupler. The two strands share one MPS
under a **layout** (Definition 11.6).

**Proposition 11.2 (crossing budget).** Bell-couple (fourier on the A
partner, `cshift(d,d)` onto the B partner) every dimension-consistent pair
whose dimension lies in a level set `L`. Then the entanglement across the
A|B bipartition is exactly

```text
  S(A|B) = Σ_{pairs p active}  log2 d_p       bits,
  rank(A|B) = Π_{pairs p active} d_p ,
```

with a flat spectrum, summed over the *active pairs* (dimension in `L`).

*Proof.* Each active pair becomes the generalized Bell state
`Σ_k |k⟩_{A[i]}|k⟩_{B[j]}/√d` — one leg in `A`, one in `B` — of Schmidt
rank `d` and flat spectrum across the A|B cut. Inactive sites stay in
`|0⟩`. The global state is a tensor product of these pairs with a product
remainder, and Schmidt ranks/entropies multiply/add across tensor
products. ∎

The budget is measured to machine precision against the closed form and
cross-validated against dense simulation for every single level
(`examples/crossing_vees.rs`, finding 24). Its content is in the
*multiplicity*:

**Corollary 11.3 (multiplicity beats dimension).** On a peak
`diamond(1, hi)` with antiparallel pairing, the single-level A|B entropy is
`m(d)·log2 d` where `m(d)` is the number of pairs at dimension `d`. The
peak `d = hi` is *unique* (`m = 1`); every lower level is *paired*
(`m = 2`). So the maximum is at the **shoulder** `d = hi−1`, not the peak,
whenever `2 log2(hi−1) > log2 hi`, i.e. `(hi−1)² > hi` — true for all
`hi ≥ 3`. For `[1..5..1]`: `4` bits at the shoulder against `2.32` at the
peak. A palindrome's peak is a bottleneck of *multiplicity*, not of
dimension.

**Corollary 11.4 (the dimension-1 pinch).** A coupler at `d = 1` is the
identity (`cshift(1,1) = [1]`, `xswap(1,·) = I`), so level-1 coupling
injects nothing — the pinch is a *decoupled* crossing. A `d = 1` site
carries no entanglement of its own (its physical space is one-dimensional)
yet still routes bond correlation past it: a spacer, not a wall. (Verified
end-to-end: the engine handles `d = 1` sites and long-range gates threaded
through them at fidelity 1.)

**Proposition 11.5 (cycles add).** Bell couplings at disjoint level sets
act on disjoint sites, hence commute, and their A|B entropies sum. So
"interact at level `d` in cycle `k`" accumulates linearly and
order-independently, up to the total budget `Σ_i log2 P_i` reached when
every level is active. Measured: `2.32 + 4 + 3.17 + 2 = 11.49` bits on the
peak-5 twin, one level per cycle. This is the precise sense of the user's
"interact at 5, 4, 3, 2, or 1, in a second cycle."

**The valley dual.** Since `valley[i] = hi+lo − diamond[i]`
(Def. `zigzag::valley`), the multiplicities reflect: a valley's *wide ends*
are paired and its *pinch* is unique. So a valley crosses most richly at
its rims (`2 log2 hi`, e.g. `4.64` bits at `d = 5`), a diamond at its
shoulders. Where a crossing can entangle two waves is a property of the
*shape*, not just the dimensions present.

### 11.6 Layout is the cost, not the entanglement

The performance lesson — and the direct answer to "can this be more
efficient?" — is that the inter-strand entanglement is *physics* but its
*cost* is a representation choice.

**Definition 11.6 (layouts).** A crossing orders its `2n` sites along the
shared MPS either as **Block** (`A` then `B`, so the A|B bipartition is the
single bond `n−1`) or **Interleaved** (each crossing pair on adjacent
sites).

**Proposition 11.7 (layout cost gap).** For a Bell crossing:

* **Block** layout carries the full A|B entanglement on one contiguous
  bond, so its bond dimension there is `rank(A|B) = Π d` (up to 2880 for
  twin `[1..5..1]`), and every coupler is long-range — threaded across the
  A|B boundary through `Θ(n)` intervening tensors.
* **Interleaved** layout makes every coupler nearest-neighbour, and every
  contiguous cut severs *at most one* Bell pair, so bond dimension stays
  `≤ hi` regardless of how many levels are coupled.

*Proof of the interleaved bound.* Place pair `k` (partners `A[i]`,
`B[partner(i)]`) at MPS positions `(2k, 2k+1)`. A contiguous cut falls
either inside a pair (severing that one Bell pair, rank `≤ hi`) or between
pairs `k` and `k+1` (all pairs `≤ k` fully left, all pairs `> k` fully
right — a product, rank 1). ∎

The measured gap is decisive (finding 24): the full crossing of the peak-4
twin is `χ = 144`, 94 546 parameters, ~70 ms in Block; the *same physical
state* is `χ = 4`, 88 parameters, ~30 µs Interleaved — a ~1000× cost gap
with identical `S(A|B) = 7.17` bits. A 50-site twin wave, fully crossed,
is `χ = 5` in 8 KB. The block layout's expense is not the entanglement (a
genuine 7.17 bits, and `χ = 144` is the minimum for a *contiguous* cut
carrying it) — it is the choice to make A|B contiguous. **Interleave to
compute; read the entanglement off the block bond.** The single-bond
diagnostic `Mps::bond_entropy_bits_at` (an `O(χ³)` local SVD, not the
`O(nχ³)` full sweep) makes that read cheap.

Compressibility, as everywhere in this crate, is structural: a Haar-random
cross-coupling on the same pairs saturates any bond budget. And the
interleaving that makes one X cheap is the first case of a general law —
the subject of §11.8.

### 11.8 Crossing networks and the cutwidth law

The `network` module generalizes the crossing from two strands to `K`, with
an arbitrary coupling graph, and settles what the layout lesson of §11.6
was a special case of.

**Definition 11.8 (crossing network).** A set of dimension-wave strands
with a *coupling graph* `G` whose vertices are sites `(strand, i)` and whose
edges are couplers, laid along one MPS by an ordering `π` of all sites.
Couplers are realized as graph states — a GHZ rung `Σ_k|k…k⟩` across a
group ([`Network::ghz`]), or a `cphase` on each edge ([`Network::cluster`]).

**Theorem 11.9 (the cutwidth law).** For a `cphase` graph state of uniform
dimension `d` on `G`, the MPS bond dimension across the cut after MPS
position `p` is exactly `d^{c_π(p)}`, where `c_π(p)` is the number of edges
of `G` with one endpoint in `π[0..p]` and one in `π[p..]`. Hence

```text
  max_p χ_p = d^{ cutwidth_π(G) } ,     min_π  = d^{ cutwidth(G) },
```

the minimum over orderings being `d` raised to the graph's **cutwidth**.

*Proof.* A `cphase` graph state `|G⟩ = Π_{(u,v)∈E} CZ_{uv} |+⟩^{⊗V}` is a
stabilizer (graph) state; the reduced state on a contiguous block `L = π[0..p]`
has entanglement rank `d^{r}` with `r` the `Z_d`-rank of the biadjacency
matrix `A_{L,\bar L}` of edges crossing the cut. For a simple graph that
matrix is a 0/1 incidence block; its rank equals the number of crossing
edges whenever they are `Z_d`-independent, which holds for the planar
lattices and matchings here (each crossing edge meets a distinct boundary
vertex). Maximizing over `p` gives the ordering's cutwidth; minimizing over
`π` is the graph-theoretic cutwidth. ∎

Every crossing geometry is then a reading of one integer — measured against
the theorem to the exact bond, and dense-validated (finding 25):

| geometry | coupling graph | cutwidth | max χ |
|---|---|---|---|
| single X (§11) | perfect matching | 1 (pair) | `hi` |
| bundle of `K` | `K`-stars per site | 1 | `hi` |
| multi-period crossing | matching, `p` peak-pairs | 1 | `hi` |
| overlay (X + ladder) | union of 2 matchings = 4-cycles | 2 | `(hi−1)²` |
| weave (`r×c` lattice) | 2D grid | `min(r,c)` | `d^{min(r,c)}` |

Three consequences carry the physics:

* **One direction is flat.** A bundle of `K` strands GHZ-coupled at each
  site (Corollary of Thm 11.9: `K`-stars have cutwidth 1 in column-major
  order) stays at `χ = hi` for every `K` — measured `2, 3, 5, 8` strands,
  all `χ = 5`. A GHZ across `K` parties is *one* Schmidt mode across any
  cut; widening the bundle adds parties, not width. Likewise two
  multi-period waves cross at `p` peak-regions, and the inter-strand
  entanglement grows with `p` (`11.49 → 45.97` bits for `p = 1…4`) while
  `χ` stays pinned at `hi = 5`: **many crossing points, one direction,
  still cheap.**
* **Overlaying adds one.** Crossing the *same* pair both ways
  (X and ladder) unions two matchings into disjoint 4-cycles, cutwidth 2,
  so `χ = (hi−1)²` — measured `4, 9, 16` for `hi = 3, 4, 5`. It is the
  *shoulder* squared, not the peak's: the unique peak sits in a 2-cycle
  carrying only its ladder edge (cutwidth 1), and only the paired shoulders
  reach cutwidth 2 — the multiplicity theme of Cor. 11.3, resurfacing as a
  cutwidth statement.
* **Two directions is an area law.** A weave — strands running two
  transverse ways, coupled at every intersection — is a 2D lattice, whose
  cutwidth is `min(rows, cols)`. So `χ = d^{min(rows,cols)}`: measured
  `2^rows` exactly (`4, 8, 16, 32` for `rows = 2…5`), and *flat in the
  length* (`rows = 3` is `χ = 8` at `cols = 3, 6, 9` alike). This is the
  boundary where the one-dimensional representation stops being efficient —
  precisely the area law that separates MPS from PEPS [14]. A *wave* weave
  (a lattice of `[1,2,3,2,1]` strands) inherits the law with a wave-shaped
  base: its bond profile `[1,1,2,4,6,4,2,1,1]` is thin at the `d=1` rims and
  thick at the shared peak, thickening with each added row.

The organizing statement: **the cost of a crossing network is the cutwidth
of its coupling graph, and the cutwidth is set by the *dimensionality* of
the crossing, not the number of strands.** One transverse direction —
however many strands, however many crossing points — is constant cutwidth
and classically cheap; a second direction is an area law. Finding a
minimal-cutwidth ordering is the layout-optimization problem, `NP`-hard in
general but exactly solvable for the structured graphs here, and the honest
successor to "interleave to compute" (§11.6).

### 11.9 Crossing the renormalization wave

The strands crossed so far were flat. The zigzag *wave* `wave(lo,hi,p)` is
not — it is a real-space renormalization structure (§2, findings 5–7):
valleys are fine (`d = lo`), waists coarse (`d = hi`), `merge_sites` is a
coarse-graining step, and structured multi-scale dynamics stays compressed
riding the wave. Applying the crossing to *these* strands, the crossing
inherits the scale hierarchy — three measured effects (finding 26).

**Proposition 11.10 (a crossing has a scale).** Cross two `wave(lo,hi,p)`
strands antiparallel. The dimension band `d` has multiplicity `m(d)` among
the crossing pairs, and Bell-coupling it injects `m(d)·log d` bits in
`m(d)` Schmidt modes of rank `d` each. The waist band `d = hi` has
`m = p` (one crossing per shared waist — few, coarse, *fat* modes); the
fine bands have `m = 2p` (many, *thin* modes). So the coupling's mode
content is coarse at the waists and fine at the valleys — the wave's RG
hierarchy, transported onto the inter-strand coupling.

*Proof.* Multiplicity counting on the palindromic wave: `hi` appears once
per period (`p` times), each interior band twice per period (`2p`); the
crossing budget of Prop. 11.2 then reads off `m(d)·log d` in `m(d)` modes.
∎ (Measured on `wave(1,5,2)`: level 5 → 2 crossings of 5 modes, level 2 →
4 crossings of 2 modes.)

**Proposition 11.11 (the woven wave lattice — RG-shaped area law).** Weave
`rows` copies of a wave `P` (vertical bonds at `P[c]`, horizontal at
`min(P[c],P[c+1])`, column-major). The bond dimension across the cut after
column `c` is `min(P[c],P[c+1])^{rows}` — the area law of Theorem 11.9 with
a *wave-shaped base*: thick `hi^{rows}` at the coarse waists, pinched to 1
at the `d = 1` rims, and **periodic** — a second period does not raise the
peak, because the crossing cost is local to each RG cell.

*Proof.* Column-major order makes the vertical bonds (within a column of
`rows` sites) local and the `rows` horizontal bonds between columns `c` and
`c+1` the crossing set of a cut there; those form a matching at dimension
`min(P[c],P[c+1])`, so Theorem 11.9 gives rank `min(P[c],P[c+1])^{rows}`.
The profile of column-gaps is the wave's, repeated per period. ∎

Measured: weaving `wave(1,3,2)` at 2 and 3 rows gives bond profiles
`[1,1,2,4,6,4,2,1,1,1,2,4,6,4,2,1,1]` and `[…,8,12,12,8,…]` — two humps at
the two waists, thickening `6 → 12` with the added row, valleys pinched to
1, dense-validated. **The RG structure decides where the area law bites:
hardest at the coarse blocks.**

**Proposition 11.12 (RG-covariance of a crossing).** Let `|Ψ⟩` be a crossed
state and `R` an exact coarse-graining (`merge_sites`) on a window of one
strand disjoint from the cut being read. Then `R|Ψ⟩` is again a crossed
state with the coupling carried onto the coarse block, and `R` is
invertible (`split_site`) with `split ∘ merge = 1`: coarse-graining and
crossing commute.

*Proof.* `merge_sites` is the associativity isomorphism
`C^{d_i} ⊗ C^{d_{i+1}} ≅ C^{d_i d_{i+1}}` (Prop. 5.1) applied to `|Ψ⟩` as a
vector; it changes the tensor factorization, not the state, so any coupling
supported on the merged sites now acts on a sublevel of the coarse site,
and `split_site` inverts it exactly. ∎

Measured: crossing twin `[1,2,3,4,3,2,1]`, then merging strand A's waist
cell into one `d = 36` site and splitting back, returns the original at
fidelity `1.000000000000`. The crossing coupling rides the RG step onto the
coarse waist and back — a covariant object under the wave's
renormalization.

Together: **the wave supplies a scale ladder, the crossing supplies
coupling, and the two compose** — the coupling acquires a scale, a 2D
lattice of waves carries an area law shaped like the wave, and the RG
transformation commutes with the crossing. The complex, expanding,
multi-scale geometry is still costed by one number, the cutwidth, now
modulated by the renormalization structure of the strands. What remains
open is the genuine RG *flow*: iterating coarse-graining on a woven wave
lattice and asking whether the effective crossing coupling runs to a fixed
point — the natural meeting of this crate's MERA-like geometry with its
crossing calculus (§13, problem 7).


## 12. Numerical foundations

The entire stack — complex arithmetic, RNG, QR, SVD — is in-crate and
dependency-free, so every numerical claim above rests on ~900 audited lines
(`c64`, `mat`). The design decisions that carry theory-level weight:

* **One-sided Jacobi SVD**, chosen for verifiable robustness at the small
  matrix sizes this library produces, with two hardening rules that matter
  for tensor networks specifically:
  1. *numerically null columns are frozen* before rotations grind them into
     the denormal range (where Gram entries underflow and rotation
     parameters turn to garbage, corrupting the *significant* subspace) —
     the pathology arises structurally here, because low-rank states
     threaded by MPOs produce matricizations with large null spaces
     (regression: `long_range_gate_is_norm_preserving_with_null_bonds`);
  2. *singular values below roundoff resolution are reported as exact
     zeros* rather than normalizing rotation debris into fake directions.
* **Truncation policy** (`TruncSpec`): a hard rank cap plus a *relative*
  discarded-weight cutoff, with every discard tallied and the retained part
  renormalized. Strongly rank-capped truncations of large matrices go
  through a randomized range-finder (Halko-style, one power iteration [10])
  whose unsampled weight is *charged to the discard tally* — approximation
  is never silent.
* **Canonical-gauge discipline** (§3): lossy operations only against
  canonical environments; exact rank-revealing sweeps arranged
  meet-in-the-middle so no SVD runs at raw-product width on both sides.
* **Determinism**: SplitMix64 with fixed seeds; Haar unitaries by Ginibre +
  two-round MGS QR with the phase fix [11]. Every experiment is
  reproducible bit-for-bit.

The epistemology of the crate follows from §7.4's lesson: *rank-side*
guarantees (exact canonicalization, discard tallies) and *weight-side*
guarantees (f64 phase resolution, relative cutoffs) are different
promises, and the test suite exercises both — 125 tests, with every
structural mechanism cross-validated against the dense ground-truth
simulator and, where possible, against closed-form laws (Schmidt spectra,
entropy values, moment laws, convergence rates) rather than against
snapshots of its own output.


## 13. Dictionary, and open problems

The unifying observation of the library is that **one representation — an
MPS over a heterogeneous chain — supports five semantic layers**, and that
its bond dimension means something different, and true, in each:

| layer | object | carrier | bond dimension reads as | boundary reads as |
|---|---|---|---|---|
| state | `\|ψ⟩` | MPS | entanglement across the cut (§2) | — |
| block | circuit block `U` | Choi MPS on `d²` | operator entanglement: past↔future pipeline width (§6) | — |
| machine | transducer | 0/1 MPO | width of the classical message front (§8) | drop = mod N; postselect = exact division; trace = mod N−1 |
| flow | `t ↦ U^t` | MPO family | operation width — integer-quantized (§9) | — |
| network | `K` crossed strands | MPS on all sites | `d^cutwidth` of the coupling graph — layout-dependent, one direction flat, two an area law (§11) | — |

and, orthogonally, that **number theory surfaces as operator properties**:

```text
gcd(k, N) = g            ↔   unitarity defect 1 − 1/g               (Prop 8.3)
k mod S, gcd with L      ↔   operator width across the cut          (Thm 8.5)
N ∓ 1 prime              ↔   closed (looped/twisted) machine unitary (Cor 10.3, Thm 10.7)
quadratic Gauss sums     ↔   near-orthogonality of imputations      (Prop 9.6)
ones'-complement seam    ↔   Jordan block, 1/(k²+1) stabilization   (Prop 10.4)
F⁴ = I                   ↔   the QFT's own flow, four operators deep (Prop 9.7)
```

with two independent budgets — geometric rank (Prop 2.1) and spectral
weight (§7.4) — governing what any of it costs.

Problems this document sharpens beyond the README's open directions:

1. **The width function of multiplication.** Theorem 8.5 reduces
   `k ↦ χ_cut(M_k)` to the counting problem `|{⌊kb/S⌋ mod L}|`, and the
   `width` module now computes it wholesale — the full atlas of `Z_2880`
   and the wave's squaring orbits (finding 20, §8.10). What remains open
   is the *closed form*: the `k ≤ S` and `k ≡ 1 (mod S)` regimes are
   solved (`min(k, L)` and `L/gcd(a, L)`), but the general `k > S` count
   is still enumeration, and the atlas's measured structure (87%
   cap-saturation, the width-2/3 resonant tail) awaits a theorem — as does
   the atlas for streamed permutations beyond the affine family.
2. **The alphabet gap.** The first edition conjectured that some dual
   sweep always "attains the cut rank" — trivially true as stated, since
   Theorem 8.5's `Γ` *is* the forward machine's realized crossing set.
   The atlas reframes the real question and answers it negatively: 584 of
   768 multipliers on `Z_2880` have intrinsic width below every
   single-front *construction alphabet* `min(k, k⁻¹, N−k, (N−k)⁻¹)` by a
   factor > 8 — extreme case `×1441`, width 2 against alphabet 1439
   (§8.10). Open: a machine model whose construction alphabet meets the
   intrinsic width — multi-front sweeps, branching messages, or a
   composition calculus with certified intermediate widths.
3. **Flows beyond finite-order frames.** For the standard-order QFT the
   problem dissolved: `F⁴ = I` makes `F^t` an exact four-term operator
   combination, no eigenframe needed (Proposition 9.7, finding 22). The
   genuinely open cases: frames without small order — the reversed-digit
   frame `V = ρ∘F` obeys no low power identity in general, so its
   fractional flow needs an actual eigenframe — and frame flows *at wave
   scale*, where the standard-order frame is unavailable (the reversal
   stage saturates χ) and a reversal-free fractional transform is the
   missing object.
4. **Width-optimal imputation.** Proposition 9.6 shows factorization routes
   with identical endpoints and `O(1/N²)` mutual overlap; the causal route
   beat the geodesic by `χ = 8` vs `11`. Characterize the minimal-width
   path between `I` and a given block — a discrete geodesic problem in
   operator space with bond dimension as the metric.
5. **Exact finite-`N` breathing law.** Proposition 9.4 is asymptotic; the
   finite-`N` participation of the Dirichlet kernel should admit a closed
   form (the measured values match the limit to 4 digits already at
   `N = 2880`).
6. **The boundary-design calculus.** §10 now exhibits *four* boundary
   regimes, and the twist axis is partly mapped: the identity and reversal
   twists yield mod `N−1` and mod `N+1` arithmetic (Theorem 10.7,
   finding 21). Open: the remaining twists (a general `σ` mixes affine
   branches over different domains — what algebra do they generate?),
   weighted and partial traces interpolating open ↔ closed, and *coupled*
   loops feeding one machine's exit into another's entry. Which rings,
   attractors, and convergence laws are reachable by boundary engineering
   alone, for a fixed transducer body, remains wide open — but the
   instrument set now includes the twist.
7. **Networks of strands.** §11.8 settles the two-strand question into a
   law — the cost of a crossing network is its coupling graph's cutwidth
   (Theorem 11.9, finding 25) — and confirms the area law for a woven
   lattice directly. What remains open is the *optimization*: cutwidth
   minimization is `NP`-hard in general, and the structured graphs here
   (bundles, weaves, overlays) are the easy cases; a crossing pattern of
   many waves at *arbitrary* angles poses a real layout search, and the
   `xswap`-as-disentangler question (can subspace exchange relocate
   inter-strand entanglement toward chosen cuts, lowering the realized
   cutwidth below the naive graph value?) is now concretely measurable. And
   §11.9 opens the RG axis: crossing the multi-V renormalization wave gives
   the coupling a scale and shapes the area law like the wave (finding 26),
   but the genuine renormalization *flow* — iterate `merge_sites` on a woven
   wave lattice and watch the effective crossing coupling run, toward a
   fixed point or away — is untouched, and is where this crate's MERA-like
   geometry and its crossing calculus would truly meet. The physics beyond:
   fermionic or frustrated couplers on a weave, and whether any crossing
   network with a genuinely two-dimensional coupling graph can dodge the
   area law through the heterogeneity of the dimension wave.


## References

[1] U. Schollwöck, *The density-matrix renormalization group in the age of
matrix product states*, Ann. Phys. **326**, 96 (2011).

[2] G. Vidal, *Efficient classical simulation of slightly entangled quantum
computations*, Phys. Rev. Lett. **91**, 147902 (2003).

[3] G. Vidal, *Entanglement renormalization*, Phys. Rev. Lett. **99**,
220405 (2007); I. Markov and Y. Shi, *Simulating quantum computation by
contracting tensor networks*, SIAM J. Comput. **38**, 963 (2008).

[4] M.-D. Choi, *Completely positive linear maps on complex matrices*,
Linear Algebra Appl. **10**, 285 (1975).

[5] P. Zanardi, *Entanglement of quantum evolutions*, Phys. Rev. A **63**,
040304(R) (2001).

[6] D. Coppersmith, *An approximate Fourier transform useful in quantum
factoring*, IBM Research Report RC 19642 (1994), arXiv:quant-ph/0201067.

[7] T. G. Draper, *Addition on a quantum computer*,
arXiv:quant-ph/0008033 (2000).

[8] G. M. Crosswhite and D. Bacon, *Finite automata for caching in matrix
product algorithms*, Phys. Rev. A **78**, 012356 (2008).

[9] P. W. Shor, *Polynomial-time algorithms for prime factorization and
discrete logarithms on a quantum computer*, SIAM J. Comput. **26**, 1484
(1997); V. Vedral, A. Barenco, and A. Ekert, *Quantum networks for
elementary arithmetic operations*, Phys. Rev. A **54**, 147 (1996).

[10] N. Halko, P.-G. Martinsson, and J. A. Tropp, *Finding structure with
randomness*, SIAM Review **53**, 217 (2011).

[11] F. Mezzadri, *How to generate random matrices from the classical
compact groups*, Notices Amer. Math. Soc. **54**, 592 (2007).

[12] J. H. McClellan and T. W. Parks, *Eigenvalue and eigenvector
decomposition of the discrete Fourier transform*, IEEE Trans. Audio
Electroacoust. **20**, 66 (1972); Ç. Candan, M. A. Kutay, and
H. M. Ozaktas, *The discrete fractional Fourier transform*, IEEE Trans.
Signal Process. **48**, 1329 (2000).

[13] R. C. Agarwal and C. S. Burrus, *Fast convolution using Fermat number
transforms with applications to digital filtering*, IEEE Trans. Acoust.
Speech Signal Process. **22**, 87 (1974); L. M. Leibowitz, *A simplified
binary arithmetic for the Fermat number transform*, IEEE Trans. Acoust.
Speech Signal Process. **24**, 356 (1976).

[14] F. Verstraete and J. I. Cirac, *Renormalization algorithms for
quantum-many body systems in two and higher dimensions*,
arXiv:cond-mat/0407066 (2004); J. Eisert, M. Cramer, and M. B. Plenio,
*Area laws for the entanglement entropy*, Rev. Mod. Phys. **82**, 277
(2010).
