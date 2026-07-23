//! Generalized qudit gates, including gates that couple sites of *different*
//! local dimension — the raw material for cross-scale interactions.
//!
//! ## Conventions
//!
//! * A single-site gate on a dimension-`d` site is a `d×d` unitary.
//! * A two-site gate on the ordered pair `(a, b)` with dimensions `(da, db)`
//!   is a `(da·db)×(da·db)` unitary over the basis `|x_a, x_b⟩` with index
//!   `x_a·db + x_b` (first site is the most-significant digit).
//! * Gates act as `|out⟩ = G |in⟩`, i.e. columns index the input basis state.
//!
//! ## Cross-dimension couplers
//!
//! Between sites of different dimension there is no canonical SWAP, but there
//! are natural unitaries:
//!
//! * [`cshift`] — the heterogeneous controlled-shift `|a,b⟩ → |a, b+a mod db⟩`
//!   (a generalized CNOT that is well-defined for any dimension pair);
//! * [`cphase`] — `|a,b⟩ → exp(2πi·ab/L)|a,b⟩` with `L = lcm(da, db)`;
//! * [`xswap`] — the *subspace exchange*: it swaps the shared
//!   `min(da,db)`-dimensional subspaces of the two sites and leaves the slack
//!   levels alone. `xswap(2,5)` literally exchanges a qubit with the qubit
//!   view of a 5-dit: applied twice it is the identity, which is what makes a
//!   small site and a large site a **bidirectional pair**.

use crate::c64::C64;
use crate::mat::{random_unitary, Mat, Rng};
use std::f64::consts::PI;

/// `d×d` identity.
pub fn identity(d: usize) -> Mat {
    Mat::identity(d)
}

/// Quantum Fourier transform on a `d`-level system:
/// `F[j,k] = exp(2πi·jk/d)/√d`. `fourier(2)` is the Hadamard gate.
pub fn fourier(d: usize) -> Mat {
    let norm = 1.0 / (d as f64).sqrt();
    Mat::from_fn(d, d, |j, k| {
        C64::cis(2.0 * PI * (j * k) as f64 / d as f64).scale(norm)
    })
}

/// Generalized Pauli-X (cyclic shift): `X|k⟩ = |k+1 mod d⟩`.
pub fn shift(d: usize) -> Mat {
    Mat::from_fn(d, d, |j, k| {
        if j == (k + 1) % d {
            C64::ONE
        } else {
            C64::ZERO
        }
    })
}

/// Generalized Pauli-Z (clock): `Z|k⟩ = exp(2πi·k/d)|k⟩`.
pub fn clock(d: usize) -> Mat {
    Mat::from_fn(d, d, |j, k| {
        if j == k {
            C64::cis(2.0 * PI * k as f64 / d as f64)
        } else {
            C64::ZERO
        }
    })
}

/// Hadamard on a qubit (= `fourier(2)`).
pub fn hadamard() -> Mat {
    fourier(2)
}

/// Diagonal phase gate `|k⟩ → exp(i·phases[k])|k⟩`.
pub fn phase_diag(phases: &[f64]) -> Mat {
    let d = phases.len();
    Mat::from_fn(d, d, |j, k| {
        if j == k {
            C64::cis(phases[k])
        } else {
            C64::ZERO
        }
    })
}

/// Digit complement `|a⟩ → |d-1-a⟩` — the reflection of a site's levels.
/// Applied to every site of a chain it maps the ring value `x` to
/// `N-1-x`: composing with a `+1` carry cascade gives negation
/// `x → -x mod N` at total machine width 2, however large `N` is.
pub fn complement(d: usize) -> Mat {
    Mat::from_fn(
        d,
        d,
        |j, k| {
            if j == d - 1 - k {
                C64::ONE
            } else {
                C64::ZERO
            }
        },
    )
}

/// Haar-random `d×d` unitary.
pub fn random(d: usize, rng: &mut Rng) -> Mat {
    random_unitary(d, rng)
}

/// Heterogeneous controlled-shift `|a,b⟩ → |a, (b + a) mod db⟩` on a
/// `(dc, dt)` pair — the generalized CNOT. Well-defined for any dimensions.
pub fn cshift(dc: usize, dt: usize) -> Mat {
    let n = dc * dt;
    Mat::from_fn(n, n, |row, col| {
        let (a, b) = (col / dt, col % dt);
        let out = a * dt + (b + a) % dt;
        if row == out {
            C64::ONE
        } else {
            C64::ZERO
        }
    })
}

/// Inverse of [`cshift`]: `|a,b⟩ → |a, (b - a) mod db⟩`.
pub fn cshift_inv(dc: usize, dt: usize) -> Mat {
    cshift(dc, dt).adjoint()
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// Heterogeneous controlled-phase `|a,b⟩ → exp(2πi·ab/L)|a,b⟩` with
/// `L = lcm(da, db)`. For `da == db == d` this is the standard qudit CZ.
pub fn cphase(da: usize, db: usize) -> Mat {
    let l = da / gcd(da, db) * db;
    let n = da * db;
    Mat::from_fn(n, n, |row, col| {
        if row != col {
            return C64::ZERO;
        }
        let (a, b) = (col / db, col % db);
        C64::cis(2.0 * PI * (a * b) as f64 / l as f64)
    })
}

/// Fractional controlled phase `|a,b⟩ → exp(2πi·ab/denom)|a,b⟩` for an
/// arbitrary real denominator. This is the **tail-radix coupling**: in the
/// mixed-radix Fourier transform over `Z_N`, output digit `j` couples to
/// input digit `i ≥ j` with exactly this gate, `denom` being the product of
/// the radix segment `d_j·d_{j+1}···d_i` (see [`crate::radix`]). Diagonal,
/// hence unitary for any `denom`.
pub fn cp_frac(da: usize, db: usize, denom: f64) -> Mat {
    let n = da * db;
    Mat::from_fn(n, n, |row, col| {
        if row != col {
            return C64::ZERO;
        }
        let (a, b) = (col / db, col % db);
        C64::cis(2.0 * PI * (a * b) as f64 / denom)
    })
}

/// Subspace exchange between a `da`-dim and a `db`-dim site: with
/// `m = min(da, db)`, maps `|a,b⟩ → |b,a⟩` whenever both `a < m` and `b < m`,
/// and acts as the identity otherwise. A permutation, hence unitary; it is an
/// involution on any state supported in the shared subspace.
///
/// For equal dimensions this is the full SWAP.
pub fn xswap(da: usize, db: usize) -> Mat {
    let m = da.min(db);
    let n = da * db;
    Mat::from_fn(n, n, |row, col| {
        let (a, b) = (col / db, col % db);
        let out = if a < m && b < m { b * db + a } else { col };
        if row == out {
            C64::ONE
        } else {
            C64::ZERO
        }
    })
}

/// Full SWAP between two equal-dimension sites.
pub fn swap(d: usize) -> Mat {
    xswap(d, d)
}

/// Controlled unitary with qudit control: `|a,b⟩ → |a⟩ (U^a)|b⟩`.
/// For `dc = 2` this is the ordinary controlled-U.
pub fn controlled_pow(dc: usize, u: &Mat) -> Mat {
    let dt = u.rows;
    assert_eq!(u.rows, u.cols);
    let n = dc * dt;
    let mut out = Mat::zeros(n, n);
    for a in 0..dc {
        let ua = u.pow(a);
        for r in 0..dt {
            for c in 0..dt {
                out.set(a * dt + r, a * dt + c, ua.at(r, c));
            }
        }
    }
    out
}

/// Tensor product of two single-site gates viewed as a two-site gate on the
/// ordered pair (first ⊗ second).
pub fn kron2(a: &Mat, b: &Mat) -> Mat {
    a.kron(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_generators_are_unitary() {
        let mut rng = Rng::new(3);
        for d in 2..=6 {
            assert!(fourier(d).is_unitary(1e-12), "fourier {}", d);
            assert!(shift(d).is_unitary(1e-12), "shift {}", d);
            assert!(clock(d).is_unitary(1e-12), "clock {}", d);
            assert!(random(d, &mut rng).is_unitary(1e-10), "random {}", d);
        }
        for &(a, b) in &[(2, 2), (2, 3), (3, 2), (2, 5), (5, 2), (4, 5), (3, 4)] {
            assert!(cshift(a, b).is_unitary(1e-12), "cshift {} {}", a, b);
            assert!(cphase(a, b).is_unitary(1e-12), "cphase {} {}", a, b);
            assert!(xswap(a, b).is_unitary(1e-12), "xswap {} {}", a, b);
            assert!(
                controlled_pow(a, &fourier(b)).is_unitary(1e-12),
                "cpow {} {}",
                a,
                b
            );
        }
    }

    #[test]
    fn fourier2_is_hadamard() {
        let h = hadamard();
        let r = 1.0 / 2.0_f64.sqrt();
        assert!((h.at(0, 0).re - r).abs() < 1e-14);
        assert!((h.at(1, 1).re + r).abs() < 1e-14);
        assert!(h.at(1, 0).im.abs() < 1e-14);
    }

    #[test]
    fn cshift_acts_on_basis() {
        // |a=2, b=1⟩ on (3,4) → |2, 3⟩
        let g = cshift(3, 4);
        let col = 2 * 4 + 1;
        let row = 2 * 4 + 3;
        assert_eq!(g.at(row, col), C64::ONE);
        assert_eq!(g.at(col, col), C64::ZERO);
    }

    #[test]
    fn xswap_is_involution_and_exchanges_subspace() {
        for &(a, b) in &[(2, 5), (3, 4), (4, 4)] {
            let g = xswap(a, b);
            assert!(
                g.mul(&g).max_abs_diff(&Mat::identity(a * b)) < 1e-14,
                "involution {} {}",
                a,
                b
            );
        }
        // (2,5): |1, 0⟩ → |0, 1⟩ ; slack |0, 4⟩ untouched.
        let g = xswap(2, 5);
        assert_eq!(g.at(1, 5), C64::ONE); // col |1,0⟩=5 → row |0,1⟩=1
        assert_eq!(g.at(4, 4), C64::ONE); // |0,4⟩ fixed
    }

    #[test]
    fn shift_clock_commutation() {
        // Z X = ω X Z with ω = exp(2πi/d).
        for d in 2..=5 {
            let zx = clock(d).mul(&shift(d));
            let xz = shift(d).mul(&clock(d)).scale(C64::cis(2.0 * PI / d as f64));
            assert!(zx.max_abs_diff(&xz) < 1e-12, "d={}", d);
        }
    }
}
