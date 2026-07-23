//! Hosted qubits: sub-simulating binary systems inside a qudit chain.
//!
//! A dimension-`d` site contains `k = floor(log2 d)` binary degrees of
//! freedom: levels `0 .. 2^k` are read as a big-endian `k`-bit register and
//! levels `≥ 2^k` are *slack*. A [`HostRegister`] assigns a global qubit
//! index to every hosted bit of a dimension profile, and lifts qubit gates to
//! qudit gates that act on the embedded register and leave slack levels
//! untouched. Because lifted gates never mix register and slack levels, a
//! state prepared in the register subspace stays there, and the hosted qubit
//! computation is *exactly* recoverable — the binary system is a view into
//! the qudit chain, not an approximation.
//!
//! The diamond profile `2,3,4,5,4,3,2` hosts `1+1+2+2+2+1+1 = 10` qubits in
//! a 2880-dimensional space (overhead ×2.8 over `2^10 = 1024`): the edge
//! sites are *literal* qubits (the binary set on either side), while the
//! waist hosts two logical qubits per site.
//!
//! [`qubit_bit_swap`] realizes the "bidirectional pair" between a physical
//! qubit site and one logical bit inside a larger site: a full SWAP of the
//! qubit with bit `t` of the host register (slack untouched, involution on
//! the register subspace). Combined with lifted gates this lets edge qubits
//! shuttle *into* the high-dimensional waist, be acted on at coarse scale,
//! and shuttle back.

use crate::c64::C64;
use crate::circuit::{Circuit, Op};
use crate::dense::DenseState;
use crate::mat::Mat;

/// Number of qubits a dimension-`d` site hosts: `floor(log2 d)`.
pub fn capacity(d: usize) -> usize {
    assert!(d >= 1);
    (usize::BITS - 1 - d.leading_zeros()) as usize
}

#[inline]
fn bit_of(level: usize, t: usize, k: usize) -> usize {
    (level >> (k - 1 - t)) & 1
}

#[inline]
fn with_bit(level: usize, t: usize, k: usize, b: usize) -> usize {
    let mask = 1usize << (k - 1 - t);
    (level & !mask) | (b * mask)
}

/// Lift a gate on the full `k`-bit register of a site (a `2^k × 2^k` matrix)
/// to the site's `d × d` space: acts on the register subspace, identity on
/// slack levels.
pub fn lift_register_gate(g: &Mat, d: usize) -> Mat {
    let k = capacity(d);
    let reg = 1usize << k;
    assert_eq!(g.rows, reg, "gate must act on the full {}-bit register", k);
    Mat::from_fn(d, d, |row, col| {
        if row < reg && col < reg {
            g.at(row, col)
        } else if row == col {
            C64::ONE
        } else {
            C64::ZERO
        }
    })
}

/// A 2×2 gate on bit `t` of a `k`-bit register, as a `2^k × 2^k` matrix.
pub fn one_bit_gate(g2: &Mat, k: usize, t: usize) -> Mat {
    assert!(t < k);
    let left = Mat::identity(1 << t);
    let right = Mat::identity(1 << (k - 1 - t));
    left.kron(g2).kron(&right)
}

/// A 4×4 gate on bits `(ta, tb)` of a `k`-bit register (gate basis index is
/// `bit_a·2 + bit_b`), as a `2^k × 2^k` matrix.
pub fn two_bit_gate(g4: &Mat, k: usize, ta: usize, tb: usize) -> Mat {
    assert!(ta < k && tb < k && ta != tb);
    assert_eq!(g4.rows, 4);
    let n = 1usize << k;
    let mut m = Mat::zeros(n, n);
    for col in 0..n {
        let (u, v) = (bit_of(col, ta, k), bit_of(col, tb, k));
        for up in 0..2 {
            for vp in 0..2 {
                let row = with_bit(with_bit(col, ta, k, up), tb, k, vp);
                m.add_at(row, col, g4.at(up * 2 + vp, u * 2 + v));
            }
        }
    }
    m
}

/// Full SWAP between a physical qubit and bit `t` of the register of a
/// dimension-`d_host` site, as a gate on the ordered pair
/// `(qubit site, host site)`. Slack levels of the host are untouched.
pub fn qubit_bit_swap(d_host: usize, t: usize) -> Mat {
    let k = capacity(d_host);
    assert!(t < k);
    let reg = 1usize << k;
    let n = 2 * d_host;
    Mat::from_fn(n, n, |row, col| {
        let (a, l) = (col / d_host, col % d_host);
        let out = if l < reg {
            let b = bit_of(l, t, k);
            b * d_host + with_bit(l, t, k, a)
        } else {
            col
        };
        if row == out {
            C64::ONE
        } else {
            C64::ZERO
        }
    })
}

/// Where a lifted two-qubit gate lands.
pub enum Lifted2 {
    /// Both qubits live in the same site: a one-site qudit gate.
    Same(usize, Mat),
    /// Qubits in different sites: a two-site gate on the ordered site pair.
    Cross(usize, usize, Mat),
}

/// Assignment of global qubit indices to the hosted bits of a dimension
/// profile. Global order: sites ascending, then bit offset ascending
/// (offset 0 is the most significant bit of a site's register).
pub struct HostRegister {
    pub dims: Vec<usize>,
    pub caps: Vec<usize>,
    /// `qubit_map[q] = (site, bit offset)`.
    pub qubit_map: Vec<(usize, usize)>,
}

impl HostRegister {
    pub fn new(dims: &[usize]) -> HostRegister {
        let caps: Vec<usize> = dims.iter().map(|&d| capacity(d)).collect();
        let mut qubit_map = Vec::new();
        for (site, &k) in caps.iter().enumerate() {
            for t in 0..k {
                qubit_map.push((site, t));
            }
        }
        HostRegister {
            dims: dims.to_vec(),
            caps,
            qubit_map,
        }
    }

    pub fn num_qubits(&self) -> usize {
        self.qubit_map.len()
    }

    /// Lift a single-qubit gate on global qubit `q` to a one-site qudit gate.
    pub fn lift1(&self, g2: &Mat, q: usize) -> (usize, Mat) {
        assert_eq!(g2.rows, 2);
        let (site, t) = self.qubit_map[q];
        let k = self.caps[site];
        let reg_gate = one_bit_gate(g2, k, t);
        (site, lift_register_gate(&reg_gate, self.dims[site]))
    }

    /// Lift a two-qubit gate on global qubits `(qa, qb)` (gate basis index
    /// `bit_a·2 + bit_b`).
    pub fn lift2(&self, g4: &Mat, qa: usize, qb: usize) -> Lifted2 {
        assert_eq!(g4.rows, 4);
        assert!(qa != qb);
        let (sa, ta) = self.qubit_map[qa];
        let (sb, tb) = self.qubit_map[qb];
        if sa == sb {
            let k = self.caps[sa];
            let reg_gate = two_bit_gate(g4, k, ta, tb);
            return Lifted2::Same(sa, lift_register_gate(&reg_gate, self.dims[sa]));
        }
        let (da, db) = (self.dims[sa], self.dims[sb]);
        let (ka, kb) = (self.caps[sa], self.caps[sb]);
        let (rega, regb) = (1usize << ka, 1usize << kb);
        let n = da * db;
        let mut m = Mat::zeros(n, n);
        for la in 0..da {
            for lb in 0..db {
                let col = la * db + lb;
                if la < rega && lb < regb {
                    let (u, v) = (bit_of(la, ta, ka), bit_of(lb, tb, kb));
                    for up in 0..2 {
                        for vp in 0..2 {
                            let row = with_bit(la, ta, ka, up) * db + with_bit(lb, tb, kb, vp);
                            m.add_at(row, col, g4.at(up * 2 + vp, u * 2 + v));
                        }
                    }
                } else {
                    m.add_at(col, col, C64::ONE);
                }
            }
        }
        Lifted2::Cross(sa, sb, m)
    }

    /// The chain basis digits encoding the qubit basis state `|bits⟩`
    /// (slack levels unused).
    pub fn qubit_basis_to_digits(&self, bits: &[usize]) -> Vec<usize> {
        assert_eq!(bits.len(), self.num_qubits());
        let mut digits = vec![0usize; self.dims.len()];
        for (q, &b) in bits.iter().enumerate() {
            assert!(b < 2);
            let (site, t) = self.qubit_map[q];
            let k = self.caps[site];
            digits[site] = with_bit(digits[site], t, k, b);
        }
        digits
    }

    /// Read a chain state as a hosted-qubit state: returns the `2^Q` qubit
    /// amplitude vector plus the probability weight found outside the
    /// register subspace (0 when only lifted gates were applied). The qubit
    /// vector is *not* renormalized.
    pub fn extract_qubit_state(&self, state: &DenseState) -> (Vec<C64>, f64) {
        assert_eq!(state.dims, self.dims);
        let q = self.num_qubits();
        let mut out = vec![C64::ZERO; 1 << q];
        let mut leakage = 0.0;
        for (idx, &amp) in state.amps.iter().enumerate() {
            if amp.re == 0.0 && amp.im == 0.0 {
                continue;
            }
            let digits = state.digits_of(idx);
            let mut in_register = true;
            let mut qidx = 0usize;
            for (site, (&d, &k)) in digits.iter().zip(self.caps.iter()).enumerate() {
                let _ = site;
                if d >= (1usize << k) {
                    in_register = false;
                    break;
                }
                qidx = (qidx << k) | d;
            }
            if in_register {
                out[qidx] = amp;
            } else {
                leakage += amp.abs2();
            }
        }
        (out, leakage)
    }

    /// Compile a qubit circuit (over global qubit indices) into a chain
    /// circuit of lifted gates.
    pub fn compile(&self, qubit_ops: &[QubitOp]) -> Circuit {
        let mut c = Circuit::new(self.dims.clone());
        for op in qubit_ops {
            match op {
                QubitOp::One(q, g) => {
                    let (site, lifted) = self.lift1(g, *q);
                    c.ops.push(Op::One(site, lifted));
                }
                QubitOp::Two(qa, qb, g) => match self.lift2(g, *qa, *qb) {
                    Lifted2::Same(site, lifted) => c.ops.push(Op::One(site, lifted)),
                    Lifted2::Cross(sa, sb, lifted) => c.ops.push(Op::Two(sa, sb, lifted)),
                },
            }
        }
        c
    }

    /// Compile the same qubit circuit for a reference all-qubit chain
    /// (dims `[2; Q]`), for verification.
    pub fn compile_reference(&self, qubit_ops: &[QubitOp]) -> Circuit {
        let mut c = Circuit::new(vec![2; self.num_qubits()]);
        for op in qubit_ops {
            match op {
                QubitOp::One(q, g) => {
                    c.ops.push(Op::One(*q, g.clone()));
                }
                QubitOp::Two(qa, qb, g) => {
                    c.ops.push(Op::Two(*qa, *qb, g.clone()));
                }
            }
        }
        c
    }
}

/// A gate in a qubit-level circuit, addressed by global qubit index.
pub enum QubitOp {
    One(usize, Mat),
    Two(usize, usize, Mat),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gates;
    use crate::mat::Rng;

    #[test]
    fn capacities() {
        assert_eq!(capacity(2), 1);
        assert_eq!(capacity(3), 1);
        assert_eq!(capacity(4), 2);
        assert_eq!(capacity(5), 2);
        assert_eq!(capacity(8), 3);
        assert_eq!(capacity(9), 3);
    }

    #[test]
    fn lifted_gates_are_unitary() {
        let mut rng = Rng::new(4);
        let reg = HostRegister::new(&[2, 3, 4, 5, 4, 3, 2]);
        assert_eq!(reg.num_qubits(), 10);
        for q in 0..reg.num_qubits() {
            let (_, g) = reg.lift1(&gates::random(2, &mut rng), q);
            assert!(g.is_unitary(1e-10), "lift1 q={}", q);
        }
        for qa in 0..reg.num_qubits() {
            for qb in 0..reg.num_qubits() {
                if qa == qb {
                    continue;
                }
                match reg.lift2(&gates::random(4, &mut rng), qa, qb) {
                    Lifted2::Same(_, g) => assert!(g.is_unitary(1e-10)),
                    Lifted2::Cross(_, _, g) => assert!(g.is_unitary(1e-10)),
                }
            }
        }
    }

    #[test]
    fn qubit_bit_swap_is_unitary_involution() {
        for d in [4usize, 5, 6, 8] {
            for t in 0..capacity(d) {
                let g = qubit_bit_swap(d, t);
                assert!(g.is_unitary(1e-12));
                assert!(g.mul(&g).max_abs_diff(&Mat::identity(2 * d)) < 1e-12);
            }
        }
    }

    #[test]
    fn qubit_bit_swap_moves_a_bit() {
        // Host d=5 (k=2). Swap qubit with bit t=0 (MSB). |a=1, l=0⟩ → |0, l=2⟩.
        let g = qubit_bit_swap(5, 0);
        let col = 5; // |a=1, l=0⟩ = 1·5 + 0
        let row = 2; // |a=0, l=2⟩ = 0·5 + 2
        assert_eq!(g.at(row, col), C64::ONE);
        // Slack level 4 untouched: |1, 4⟩ fixed.
        assert_eq!(g.at(9, 9), C64::ONE);
    }

    #[test]
    fn basis_encoding_roundtrip() {
        let reg = HostRegister::new(&[2, 3, 4, 5]);
        // Qubits: site0 bit0, site1 bit0, site2 bits(0,1), site3 bits(0,1) → 6 qubits.
        assert_eq!(reg.num_qubits(), 6);
        let bits = [1, 0, 1, 1, 0, 1];
        let digits = reg.qubit_basis_to_digits(&bits);
        // site2 register (1,1) → level 3; site3 (0,1) → level 1.
        assert_eq!(digits, vec![1, 0, 3, 1]);
        let dense = DenseState::basis_state(&reg.dims, &digits);
        let (qv, leak) = reg.extract_qubit_state(&dense);
        assert!(leak < 1e-15);
        let qidx = bits.iter().fold(0usize, |acc, &b| (acc << 1) | b);
        assert!((qv[qidx].abs() - 1.0).abs() < 1e-12);
    }
}
