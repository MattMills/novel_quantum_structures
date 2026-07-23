//! Backend-agnostic circuits over a heterogeneous qudit chain.
//!
//! A [`Circuit`] is a straight-line list of one- and two-site unitaries tied
//! to a dimension profile. The same circuit runs on the exact dense backend
//! and on the MPS backend, which is how every structural claim in this crate
//! gets verified.

use crate::dense::DenseState;
use crate::mat::Mat;
use crate::mps::Mps;

/// A single circuit operation.
pub enum Op {
    /// One-site gate at `(site)`.
    One(usize, Mat),
    /// Two-site gate on the ordered pair `(a, b)` (adjacent or long-range).
    Two(usize, usize, Mat),
}

/// A gate list bound to a dimension profile.
pub struct Circuit {
    pub dims: Vec<usize>,
    pub ops: Vec<Op>,
}

impl Circuit {
    pub fn new(dims: Vec<usize>) -> Circuit {
        Circuit {
            dims,
            ops: Vec::new(),
        }
    }

    /// Append a one-site gate.
    pub fn one(&mut self, site: usize, g: Mat) -> &mut Self {
        assert_eq!(g.rows, self.dims[site], "gate/site dimension mismatch");
        assert_eq!(g.rows, g.cols);
        self.ops.push(Op::One(site, g));
        self
    }

    /// Append a two-site gate on the ordered pair `(a, b)`.
    pub fn two(&mut self, a: usize, b: usize, g: Mat) -> &mut Self {
        assert!(a != b);
        assert_eq!(g.rows, self.dims[a] * self.dims[b], "gate size mismatch");
        assert_eq!(g.rows, g.cols);
        self.ops.push(Op::Two(a, b, g));
        self
    }

    pub fn extend(&mut self, other: Circuit) -> &mut Self {
        assert_eq!(self.dims, other.dims);
        self.ops.extend(other.ops);
        self
    }

    pub fn len(&self) -> usize {
        self.ops.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ops.is_empty()
    }

    /// Run on the exact dense backend.
    pub fn run_dense(&self, state: &mut DenseState) {
        assert_eq!(state.dims, self.dims);
        for op in &self.ops {
            match op {
                Op::One(s, g) => state.apply1(*s, g),
                Op::Two(a, b, g) => state.apply2(*a, *b, g),
            }
        }
    }

    /// Run on the MPS backend (truncation policy comes from the state).
    pub fn run_mps(&self, state: &mut Mps) {
        assert_eq!(state.dims, self.dims);
        for op in &self.ops {
            match op {
                Op::One(s, g) => state.apply1(*s, g),
                Op::Two(a, b, g) => state.apply2(*a, *b, g),
            }
        }
    }
}
