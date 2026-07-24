# Complexity notes

**The work and the output of every geometric object in the crate — a
cost-model overview, a measured census, and the classification that falls
out of it.**

Companion to [THEORY.md](THEORY.md) (which derives *what* each object is)
and reproducible from `examples/complexity_census.rs` (which *measures*
every number quoted here, in ~2 seconds). Where THEORY.md asks "is it
correct?", this asks "what does it cost, and what does it produce?"

The one-sentence answer: **for every object, the bond dimension `χ` is the
complexity.** An MPS of `n` sites, local dimension `d`, bond dimension `χ`
costs `O(n·χ³·d³)` to evolve and `O(n·χ²·d)` to store, against the `O(dⁿ)`
of the dense vector it represents. So the entire census reduces to one
question per object — *how does its `χ` grow?* — and the objects sort into
three classes by the answer.


## 1. The cost model

Everything is built from a handful of primitives (`mat`, `mps`, `mpo`),
whose costs compose into the object-level figures.

| primitive | cost | where |
|---|---|---|
| Jacobi SVD of `m×n` (`m≥n`) | `O(m·n²)` (constant sweeps) | `mat::svd` |
| single-site gate | `O(χ²·d²)` | `Mps::apply1` |
| adjacent two-site gate | `O(χ³·d³)` (an SVD at bond `χd`) | `Mps::apply2_adjacent` |
| long-range gate over span `ℓ` | `O(ℓ·κ³·χ³·d)`, operator-Schmidt `κ≤d²` | `Mps::apply2_long_range` |
| canonicalize / recompress | `O(n·χ³·d)` | `Mps::recompress` |
| inner product / fidelity | `O(n·χ³·d)` | `Mps::inner` |
| single-bond entropy | `O(χ³·d)` (one local SVD) | `Mps::bond_entropy_bits_at` |
| MPO apply to state | `O(n·(χ_U·χ_ψ)³·d)` | `Mpo::apply_to` |
| MPO compose | `O(n·(χ_1·χ_2)³·d²)` | `Mpo::compose_after` |

Two master formulas follow, and they are the whole story:

```text
  store an MPS/MPO:   O(n · χ² · d)      (MPO carrier: O(n · χ_U² · d²))
  evolve g gates:     O(g · χ³ · d³)     (+ span factors for long-range)
  dense equivalent:   O(dⁿ) space, O(d^{2n}) per gate
```

The saving is `dⁿ / (n·χ²·d)` — exponential in `n` **whenever `χ` stays
bounded**, and nothing whenever it does not. Hence the classification.


## 2. Three complexity classes

Every object sits in one of three classes by how `χ` behaves. All numbers
below are measured (`complexity_census.rs`).

### Class A — constant width `χ = O(1)`: linear time and space

The cheap regime. `χ` does not grow with the chain, so cost is *linear* in
`n` and the dense saving is total.

**The carry adder MPO** `|x⟩ → |x+1 mod N⟩` — `χ = 2` for *any* ring size,
because one classical carry bit rides the bond:

```text
   profile      n    χ    params      size      build
wave(2,5,1)     7    2       316     5.1 KB      69µs      (dense 2.9e3)
wave(2,5,2)    13    2       632    10.1 KB      41µs      (dense 4.1e6)
wave(2,5,3)    19    2       948    15.2 KB      54µs      (dense 6.0e9)
wave(2,5,4)    25    2      1264    20.2 KB      79µs      (dense 8.6e12)
```

`χ` flat, params and time strictly linear in `n`, while the dense dimension
runs to `8.6·10¹²`. Also in Class A: the **diagonal phase ramp** (`χ = 1`,
a product operator), **GHZ / cluster states** on a bounded-degree graph,
and the **bundle** network — `K` crossing strands GHZ-coupled in one
direction stay at `χ = hi` however many strands (measured `χ = 5` for
`K = 2, 4, 8, 16`, cost linear in `K`).

### Class B — geometry-bounded width `χ ≤ budget`: polynomial

The structured regime. `χ` is capped by the geometry or by number theory —
bounded, but larger than `O(1)`, and often the object's headline quantity.

**The bowtie state** saturates the diamond's entanglement budget — `χ` is
the waist rank `∏` of half the dimensions, set by geometry and independent
of chain length:

```text
   profile      n    χ=budget    waist bits    params    build
diamond(2,3)    3       2           1.000          20      8µs
diamond(2,4)    5       6           2.585         224     33µs
diamond(2,5)    7      24           4.585        4112    763µs
diamond(2,6)    9     120           6.907      116432     96ms
```

**The modular multiplier `×k`** — `χ` is the *number-theoretic true width*
(the cut-rank `|{⌊kb/S⌋ mod L}|`, THEORY §8.5), non-monotone in `k`:

```text
   k    atlas χ   measured    defect
   7        7         7        0        (unitary)
  11       11        11        ~1e-16
  49       24        24        ~1e-16
  squaring orbit 7 → 49 → 2401 → 1921 has widths  7 → 24 → 6 → 3
```

`×1921` is a *thinner* operator than `×7`. Also Class B: the **mixed-radix
QFT** (reversed core `χ = 6`, but the digit-reversal / bowtie permutation
carries `χ = 36` — measured `diamond(2,4)`: `6` vs `36`), the **single X
crossing** (`χ ≤ hi` interleaved), and the **Fourier adder pipeline**,
whose *build* touches `χ = 36` but whose *output collapses to `χ = 4`* — an
object that pays Class-B work to produce a Class-A result.

### Class C — exponential width: area law and volume law

The hard regime. `χ` grows exponentially in some extent, and the
one-dimensional representation stops being efficient.

**The weave** — strands crossing in two transverse directions (a 2D
lattice) — has `χ = d^{rows}`, exponential in the transverse extent:

```text
  rows   sites     χ    2^rows    params    build
    2      12       4       4        296     200µs
    3      18       8       8       1704     495µs
    4      24      16      16       8872       3ms
    5      30      32      32      43688      21ms
```

Each added transverse row *doubles* `χ` (and roughly the work) — the **area
law** (THEORY §11.9), the MPS/PEPS boundary. **Haar-random scrambling** is
worse still — a *volume* law that fills any budget: brickwork on a wave has
truncation error `0.41 / 0.12 / 8e-3` at caps `8 / 16 / 32`, reaching the
full rank `64` before it is exact. Compressibility is a property of
*structure*, not of the geometry.


## 3. Two cross-cutting phenomena

Beyond the three classes, two facts about *where the cost lives* recur.

### Predict vs pay — the width atlas

An object's output width is usually a number-theoretic invariant, and the
`width` module computes it **without building the operator**: the cut-rank
`|{⌊kb/S⌋ mod L}|` is `O(S)` arithmetic (or `O(1)` closed form), against
`O(n·w³)` to build and recompress the tensor.

```text
   k    predict (arithmetic)    pay (build MPO)    χ
   7           5µs                    794µs         7
  11           2µs                      4ms        11
```

Predicting the cost is orders of magnitude cheaper than paying it, and
exact (the census `assert!`s equality). This is a genuine complexity
separation: **the cost of a modular-arithmetic operator is knowable in time
independent of the cost of realizing it.** It is what lets the Shor kernel
(`examples/shor_kernel.rs`) price all ten controlled powers a priori.

### Layout is the cost — cutwidth

The *same physical state* can be Class A or Class C depending only on the
site ordering. A fully Bell-crossed twin `[1,2,3,4,3,2,1]`:

```text
      layout        χ    params    build
       block      144     94546     68ms
 interleaved        4        88     42µs
```

Identical `7.17` bits of inter-strand entanglement, a `~1000×` cost gap.
The complexity of a crossing network is the **cutwidth** of its coupling
graph (THEORY §11.9): one crossing direction is constant cutwidth (cheap),
a second is an area law — and finding a minimal-cutwidth ordering is the
`NP`-hard layout-optimization problem underneath "interleave to compute."


## 4. What `χ` *means* — the width dictionary

`χ` is one number with a different meaning in each layer, and each meaning
is a complexity statement.

| layer | object | `χ` reads as | measured range |
|---|---|---|---|
| state | `\|ψ⟩` | entanglement rank across a cut | `2` (GHZ) … budget `120` |
| block | operator `U` | operator entanglement (past↔future width) | `1` (ramp) … `36` (QFT reversal) |
| machine | transducer | width of the classical message front | `2` (adder) … resonant `w ≤ k` |
| flow | `t ↦ U^t` | operation width, integer-quantized | `1–2` (whole) … `9–11` (fractional) |
| network | crossed strands | `d^{cutwidth}` of the coupling graph | `hi` (bundle) … `d^{rows}` (weave) |

Measured flow signatures (the width *seeing* arithmetic and frame
structure):

```text
FourierFlow U^t, Z_2880:   t   0.00 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00
                           χ      1    9    9    9    2    9    9    9    2
FrameFlow F^t, Z_12:       F^0=I χ=1   F^1=F χ=4   F^2=Π χ=2   F^3=F† χ=4   (F⁴=I)
```


## 5. A second complexity axis — convergence

The boundary systems (`stabilize`) are not costed by width but by
**iteration count**: a non-unitary looped operator is iterated to its
attractor, and the boundary design sets the *convergence class*.

```text
  bare seam (Jordan block, error ~ 1/k²):   1001 steps to 1e-12   (polynomial)
  damped seam γ=0.5 (eigenvalue, error ~ γ^{2k}):  20 steps       (exponential)
  damped seam γ=0.2:                               11 steps
  damped seam γ=0.05:                               7 steps
```

Same fixed point, same per-step cost `O(n·χ³·d)`; a rank-one damper moves
the whole computation from polynomial to exponential convergence — the
convergence *rate* is engineered at the boundary (THEORY §10).


## 6. Census summary

Every geometric object, its work to build, and its output — `n` sites,
`d` max local dimension, `N = ∏dᵢ`, `χ` bond dimension, `w` cascade width,
`ℓ` gate span, `g` gate count, `K` strands, `r×c` weave.

| object | build work | output `χ` | output size | class |
|---|---|---|---|---|
| carry adder MPO | `O(n·d³)` | `2` | `O(n)` — 20 KB @ n=25 | A |
| phase ramp `D_c` | `O(n·d)` | `1` | `O(n)` | A |
| GHZ / NN cluster | `O(n·d³)` | `O(1)` | `O(n)` | A |
| bundle of `K` strands | `O(K·n·hi³)` | `hi` | `O(K·n)` | A |
| bowtie state | `O(n·ℓ·budget³·d)` | budget | `O(n·budget²)` | B |
| modular multiplier `×k` | `O(n·k·d²)+O(n·w³·d)` | `w ≤ k` (resonant) | `O(n·w²)` | B |
| QFT (reversed core) | `O(n²·d⁶)` | `6–8` | `O(n·χ²·d²)` | B |
| QFT (standard, +reversal) | `O(n²·36³·d²)` | `36` | ″ | B |
| single X (interleaved) | `O(n·hi³)` | `hi` | `O(n)` | B |
| Fourier adder pipeline | `O(n·36³·d²)` | `4` (collapsed) | `O(n)` | B→A |
| `FourierFlow` `U^t` | `O(n·χ_V³·d²)` per snapshot | `2`/`9–11` | — | B |
| `FrameFlow` `F^t` | `O(n·χ_F³)` per snapshot | `1`/`2`/`36` | — | B |
| controlled `C-U^{2ʲ}` (Shor) | `O(n·(w+1)³·d²)` | `w+[k≢1]` | `O(n·χ²·d²)` | B |
| looped / twisted boundary | `O(m·n·χ³·d²)` build | modest | `O(n·χ²·d²)` | B |
| single X (block layout) | `O(n·(∏d)³)` | `∏d` | `O(n·(∏d)²)` | C* |
| weave `r×c` (2 directions) | `O(c·d^{3r})` | `d^{min(r,c)}` | `O(c·d^{2r})` | C |
| Haar scrambling | `O(g·χ³·d³)` → saturates | full rank | volume | C |

`*` the block-layout X is Class C only by its layout; interleaved it is
Class B — the cutwidth point of §3.

Reproduce every row: `cargo run --release --example complexity_census`.
The headline it ends on: the 25-site wave's `+1` adder acts on a
`7.4·10²⁵`-dimensional operator space and is stored in **20.2 KB** — a
`6·10²²×` compression, which is Class A doing what Class A does.
