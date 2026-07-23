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
pub mod circuit;
pub mod dense;
pub mod embed;
pub mod gates;
pub mod mat;
pub mod mps;
pub mod zigzag;

pub use c64::C64;
pub use mat::{Mat, Rng, TruncSpec};
