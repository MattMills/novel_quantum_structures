//! Minimal complex double-precision arithmetic.
//!
//! Implemented in-crate (rather than pulling `num-complex`) so the library is
//! fully self-contained.

use std::fmt;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

/// A complex number with `f64` real and imaginary parts.
#[derive(Clone, Copy, PartialEq, Default)]
pub struct C64 {
    pub re: f64,
    pub im: f64,
}

impl C64 {
    pub const ZERO: C64 = C64 { re: 0.0, im: 0.0 };
    pub const ONE: C64 = C64 { re: 1.0, im: 0.0 };
    pub const I: C64 = C64 { re: 0.0, im: 1.0 };

    #[inline]
    pub fn new(re: f64, im: f64) -> C64 {
        C64 { re, im }
    }

    /// Purely real complex number.
    #[inline]
    pub fn real(re: f64) -> C64 {
        C64 { re, im: 0.0 }
    }

    /// `exp(i * theta)` — the unit phase at angle `theta`.
    #[inline]
    pub fn cis(theta: f64) -> C64 {
        C64 {
            re: theta.cos(),
            im: theta.sin(),
        }
    }

    #[inline]
    pub fn conj(self) -> C64 {
        C64 {
            re: self.re,
            im: -self.im,
        }
    }

    /// Squared modulus `|z|^2`.
    #[inline]
    pub fn abs2(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    /// Modulus `|z|`.
    #[inline]
    pub fn abs(self) -> f64 {
        self.abs2().sqrt()
    }

    /// Argument (phase angle) of `z`.
    #[inline]
    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }

    /// Multiply by a real scalar.
    #[inline]
    pub fn scale(self, s: f64) -> C64 {
        C64 {
            re: self.re * s,
            im: self.im * s,
        }
    }

    /// The unit-modulus phase `z / |z|`; returns 1 for `z == 0`.
    #[inline]
    pub fn phase(self) -> C64 {
        let a = self.abs();
        if a == 0.0 {
            C64::ONE
        } else {
            self.scale(1.0 / a)
        }
    }
}

impl Add for C64 {
    type Output = C64;
    #[inline]
    fn add(self, o: C64) -> C64 {
        C64::new(self.re + o.re, self.im + o.im)
    }
}

impl Sub for C64 {
    type Output = C64;
    #[inline]
    fn sub(self, o: C64) -> C64 {
        C64::new(self.re - o.re, self.im - o.im)
    }
}

impl Mul for C64 {
    type Output = C64;
    #[inline]
    fn mul(self, o: C64) -> C64 {
        C64::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}

impl Div for C64 {
    type Output = C64;
    #[inline]
    fn div(self, o: C64) -> C64 {
        let d = o.abs2();
        C64::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
}

impl Neg for C64 {
    type Output = C64;
    #[inline]
    fn neg(self) -> C64 {
        C64::new(-self.re, -self.im)
    }
}

impl AddAssign for C64 {
    #[inline]
    fn add_assign(&mut self, o: C64) {
        self.re += o.re;
        self.im += o.im;
    }
}

impl SubAssign for C64 {
    #[inline]
    fn sub_assign(&mut self, o: C64) {
        self.re -= o.re;
        self.im -= o.im;
    }
}

impl MulAssign for C64 {
    #[inline]
    fn mul_assign(&mut self, o: C64) {
        *self = *self * o;
    }
}

impl fmt::Debug for C64 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:+.6}{:+.6}i", self.re, self.im)
    }
}

impl fmt::Display for C64 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:+.4}{:+.4}i", self.re, self.im)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic() {
        let a = C64::new(1.0, 2.0);
        let b = C64::new(-3.0, 0.5);
        assert_eq!(a + b, C64::new(-2.0, 2.5));
        assert_eq!(a - b, C64::new(4.0, 1.5));
        // (1+2i)(-3+0.5i) = -3 + 0.5i - 6i + i^2 = -4 - 5.5i
        assert_eq!(a * b, C64::new(-4.0, -5.5));
        let q = (a / b) * b - a;
        assert!(q.abs() < 1e-12);
    }

    #[test]
    fn cis_and_phase() {
        let z = C64::cis(0.7);
        assert!((z.abs() - 1.0).abs() < 1e-15);
        assert!((z.arg() - 0.7).abs() < 1e-15);
        let w = C64::new(3.0, -4.0);
        assert!((w.phase().abs() - 1.0).abs() < 1e-15);
        assert!((w.phase().scale(w.abs()) - w).abs() < 1e-12);
    }

    #[test]
    fn conj_mul_gives_abs2() {
        let z = C64::new(2.5, -1.5);
        let p = z.conj() * z;
        assert!((p.re - z.abs2()).abs() < 1e-12);
        assert!(p.im.abs() < 1e-12);
    }
}
