#![forbid(unsafe_code)]
//! # Novel quantum structures: twisted zigzag qudit chains
//!
//! A zero-dependency Rust library for exploring a family of quantum
//! structures built from *heterogeneous-dimension qudit chains* whose local
//! dimension expands and collapses along the chain — e.g. `2→3→4→5→4→3→2` —
//! with *cross-scale couplings* pairing sites across the dimension wave:
//!
//! ```text
//!   5           ●
//!   4         ●   ●
//!   3       ●       ●
//!   2     ●           ●
//!         └─────x─────┘        mirror pairs: (0,6), (1,5), (2,4)
//!           └───x───┘          — equal dimension, full xswap /
//!             └─x─┘              generalized-Bell couplable
//! ```
//!
//! The library provides two interchangeable simulation backends:
//!
//! * [`dense`] — an exact dense state-vector simulator over arbitrary mixed
//!   radix (ground truth for verification), and
//! * [`mps`] — a heterogeneous-dimension **matrix product state** simulator
//!   whose cost is polynomial in chain length for bounded bond dimension.
//!   This is what makes long zigzag chains *classically computable*: the raw
//!   Hilbert space of a many-diamond dimension wave is astronomically large,
//!   but its tensor-network representation stays linear in the number of
//!   sites. Long-range (cross-scale) two-site gates are applied through an
//!   operator-Schmidt MPO decomposition, so the "1↔5" couplings are
//!   first-class operations, not an afterthought.
//!
//! On top of the simulators sit the structural ideas this crate exists to
//! explore:
//!
//! * [`zigzag`] — dimension-wave profiles (diamonds `2,3,4,5,4,3,2` and
//!   multi-period waves), their mirror pairs and valley/waist geometry, and
//!   named cross-scale circuit families ([`zigzag::bowtie`],
//!   [`zigzag::crosswave_round`], …).
//! * [`embed`] — *hosted qubits*: a dimension-`d` site carries
//!   `floor(log2 d)` binary degrees of freedom, so qubit circuits compile
//!   onto a zigzag chain and are recovered exactly — the binary sets on
//!   either side of the wave are literal qubits, and the waist hosts two
//!   logical qubits per site. [`embed::qubit_bit_swap`] makes an edge qubit
//!   and a logical bit inside a higher site a *bidirectional pair*.
//! * [`mps`] scale morphing — [`mps::Mps::merge_sites`],
//!   [`mps::Mps::split_site`], [`mps::Mps::promote_site`],
//!   [`mps::Mps::demote_site`]: the representation itself expands and
//!   collapses while preserving the state, so a single coarse gate acts as a
//!   whole fine-scale subcircuit.
//! * [`mpo`] — **circuit-of-circuits**: whole circuit blocks reified as
//!   matrix product operators, stored (via the Choi isomorphism) as states
//!   on the doubled chain — pipeline states stretched between the past and
//!   future boundaries of the block. Blocks compose
//!   ([`mpo::Mpo::compose_after`]) and apply ([`mpo::Mpo::apply_to`]) as
//!   first-class objects, and their bond spectra measure operator
//!   entanglement — the width of the past↔future correlation pipeline.
//!   Control-selection ([`mpo::Mpo::select_on`]) makes controlled blocks
//!   first-class — with [`cascade::Transducer::mult_skipping`] it builds
//!   the `C-U^{2^j}` of phase estimation, and the full Shor kernel runs
//!   end-to-end (`examples/shor_kernel.rs`).
//! * [`radix`] — the **structured tail-radix phase web**: the mixed-radix
//!   Fourier transform over the chain's ring `Z_N`, whose output digit `j`
//!   couples only to input digits `i ≥ j` with angle `2π/(d_j···d_i)`;
//!   Draper phase ramps; and the exact bond-2 carry adder MPO. Together
//!   with [`mpo`] these realize recursive Fourier-space arithmetic on the
//!   zigzag chain (`QFT† ∘ ramp ∘ QFT` collapsing to a modular adder).
//! * [`cascade`] — **stepwise cascade operators**: finite-state transducers
//!   lifted to MPOs, with the classical message riding the virtual bond —
//!   bond dimension *is* the width of the information front. The carry
//!   adder is the `m = 2` member; modular multiplication `x → kx mod N` is
//!   the `m = k` member (unitary iff `gcd(k, N) = 1` — number theory as an
//!   operator property), and the family is closed under composition, so
//!   modular exponentiation is iterated cascade composition. Every machine
//!   has a **geometrically opposed dual** ([`cascade::Transducer::div`]):
//!   division's remainders sweep the opposite direction from
//!   multiplication's carries, computing `×k⁻¹` at width `k` instead of
//!   `k⁻¹ mod N` — enabled by quantum boundary conditions
//!   ([`cascade::Boundary`]: enter in superposition, postselect the exit).
//! * [`flow`] — **operators as state-components over operation width**: a
//!   block becomes a one-parameter flow `t ↦ U^t` (exact fractional powers
//!   in a diagonalizing Fourier frame), and the [`flow::WidthCursor`] holds
//!   a pipeline of flow segments at adaptive temporal resolution — refine
//!   at the focus, coarsen behind it, total composition invariant. Bond
//!   dimension along the flow *measures* operation width: it collapses at
//!   integer shifts and widens at fractional ones. [`flow::FrameFlow`]
//!   extends the idea to the Fourier frame itself: `F⁴ = I` makes the
//!   fractional QFT `F^t` a four-term combination of `I`, `F`, the ring
//!   reflection, and `F†` — no eigensolver needed.
//! * [`width`] — the **a-priori width calculus**: the cut-rank theorem
//!   `χ_cut(×k) = |{⌊k·b/S⌋ mod L}|` as executable number theory, so the
//!   bond profile of modular multiplication is computed *before* any
//!   tensor is built — and pinned against recompressed cascade MPOs by
//!   the tests.
//! * [`stabilize`] — **self-stabilizing boundary systems**: close a cascade's
//!   message loop ([`cascade::Transducer::to_mpo_looped`]) and it becomes a
//!   non-unitary dynamical system with a designed attractor. The traced
//!   adder is ones'-complement (mod `N-1`), turning the diamond's composite
//!   ring into a prime field that *heals* gcd obstructions; the double-zero
//!   seam is a Jordan block whose iteration self-stabilizes the number
//!   representation, at a convergence rate a rank-one damper sets.
//!   *Twisting* the loop by the message reversal
//!   ([`cascade::Transducer::to_mpo_looped_twisted`]) selects the third
//!   ring of the family — diminished-one arithmetic mod `N+1`, with the
//!   missing zero as an annihilating hole dual to the seam. Operator
//!   linear combinations ([`mpo::Mpo::add`], [`mpo::Mpo::scale`],
//!   [`mpo::Mpo::basis_transfer`]) make these boundary systems first-class.
//! * [`circuit`] — a backend-agnostic gate list so every experiment can be
//!   cross-validated dense-vs-MPS.
//!
//! Everything (complex arithmetic, RNG, QR, a one-sided Jacobi SVD) is
//! implemented in-crate; there are no dependencies.
//!
//! ## Quick start
//!
//! ```
//! use novel_quantum_structures::{mps::Mps, zigzag, TruncSpec};
//!
//! // The 2→3→4→5→4→3→2 diamond, Bell-paired across its mirror flanks.
//! let profile = zigzag::diamond(2, 5);
//! let mut psi = Mps::zero_state(&profile, TruncSpec::new(64, 1e-12));
//! zigzag::bowtie(&profile).run_mps(&mut psi);
//!
//! // Maximally entangled across the waist — yet the bond dimension is 24,
//! // not the dense 2880.
//! assert_eq!(psi.max_bond_dim(), 24);
//! let bits = psi.bond_entropies_bits();
//! assert!((bits[2] - 24f64.log2()).abs() < 1e-9);
//! ```

pub mod c64;
pub mod cascade;
pub mod circuit;
pub mod dense;
pub mod embed;
pub mod flow;
pub mod gates;
pub mod mat;
pub mod mpo;
pub mod mps;
pub mod radix;
pub mod stabilize;
pub mod width;
pub mod zigzag;

pub use c64::C64;
pub use mat::{Mat, Rng, TruncSpec};
