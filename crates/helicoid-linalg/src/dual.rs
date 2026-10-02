//! Forward-mode dual numbers (`docs/PHASE2.md` §3, decision `0003`).
//!
//! `Dual<S, N>` carries a value and `N` directional derivatives and is itself a [`Real`], so code
//! generic over `S: Real` differentiates itself, Taylor branches included.
//!
//! The value part runs the same `S` operation on the same operands as plain `S` evaluation, so it
//! is bitwise the plain result (`dual_value_is_plain_value`), except the sign and payload of a NaN
//! that arithmetic or a `libm` call produces: Rust leaves them unspecified and a release build may
//! commute operands. Comparisons and `select` read the value part only. Every derivative rule is
//! written with `S` operations, so `S` may itself be a `Dual` and second derivatives fall out of
//! nesting.
//!
//! | value | derivative component `i` |
//! |---|---|
//! | `a * b` | `a.d[i] * b.v + a.v * b.d[i]` |
//! | `a / b` | `(a.d[i] - q * b.d[i]) / b.v`, `q = a.v / b.v`: no `b.v^2` to overflow |
//! | `sqrt(x)` | `x.d[i] / (2 sqrt x.v)` |
//! | `cbrt(x)` | `x.d[i] / (3 c^2)`, `c = cbrt x.v` |
//! | `sin_cos(x)` | `(x.d[i] cos x.v, -x.d[i] sin x.v)` |
//! | `atan2(y, x)` | `(x.v y.d[i] - y.v x.d[i]) / (x.v^2 + y.v^2)` |
//! | `abs(x)` | `sgn(x.v) x.d[i]`, `sgn(+-0) = +1` |
//! | `copysign(x, s)` | `sgn(x.v) sgn(s) x.d[i]`, `sgn(+-0) = +1` for `x`, the sign bit for `s` |
//!
//! The quotient derivative is NaN once `q` overflows (the value is already infinite) and loses
//! relative accuracy while `q` is subnormal. The other derivative domains are on the methods. No
//! derivative domain is asserted: a check would panic on inputs where plain `S` evaluation does not.

use crate::real::{Blend, Mask, Precision, Real};
use core::array;
use core::ops::{Add, Div, Mul, Neg, Sub};

/// A value and its `N` directional derivatives.
///
/// There is no `PartialEq` or `PartialOrd`: a comparison of duals is a comparison of their value
/// parts, and goes through [`Real::lt`]/[`Real::le`].
#[derive(Clone, Copy, Debug)]
pub struct Dual<S, const N: usize> {
    /// The value.
    pub v: S,
    /// The derivative of the value along each of the `N` independent variables.
    pub d: [S; N],
}

impl<S: Real, const N: usize> Dual<S, N> {
    /// A constant: every derivative is zero.
    #[inline]
    pub fn constant(v: S) -> Self {
        Self {
            v,
            d: [S::zero(); N],
        }
    }

    /// The `i`-th independent variable: `d` is the `i`-th unit vector.
    ///
    /// # Domain
    ///
    /// `i < N`; checked by `debug_assert!`. A release build with `i >= N` returns a constant.
    #[inline]
    pub fn variable(v: S, i: usize) -> Self {
        debug_assert!(i < N, "Dual::variable: index out of range");
        Self {
            v,
            d: array::from_fn(|k| if k == i { S::one() } else { S::zero() }),
        }
    }
}

// The product and quotient rules live in inherent methods so that clippy's
// `suspicious_arithmetic_impl`, which reads a `+` inside `Mul::mul` as a typo, needs no local
// `allow`.
impl<S: Real, const N: usize> Dual<S, N> {
    #[inline]
    fn product(self, o: Self) -> Self {
        Self {
            v: self.v * o.v,
            d: array::from_fn(|i| self.d[i] * o.v + self.v * o.d[i]),
        }
    }

    #[inline]
    fn quotient(self, o: Self) -> Self {
        let q = self.v / o.v;
        Self {
            v: q,
            d: array::from_fn(|i| (self.d[i] - q * o.d[i]) / o.v),
        }
    }
}

impl<S: Real, const N: usize> Add for Dual<S, N> {
    type Output = Self;
    #[inline]
    fn add(self, o: Self) -> Self {
        Self {
            v: self.v + o.v,
            d: array::from_fn(|i| self.d[i] + o.d[i]),
        }
    }
}

impl<S: Real, const N: usize> Sub for Dual<S, N> {
    type Output = Self;
    #[inline]
    fn sub(self, o: Self) -> Self {
        Self {
            v: self.v - o.v,
            d: array::from_fn(|i| self.d[i] - o.d[i]),
        }
    }
}

impl<S: Real, const N: usize> Mul for Dual<S, N> {
    type Output = Self;
    #[inline]
    fn mul(self, o: Self) -> Self {
        self.product(o)
    }
}

impl<S: Real, const N: usize> Div for Dual<S, N> {
    type Output = Self;
    #[inline]
    fn div(self, o: Self) -> Self {
        self.quotient(o)
    }
}

impl<S: Real, const N: usize> Neg for Dual<S, N> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self {
            v: -self.v,
            d: self.d.map(Neg::neg),
        }
    }
}

/// Selects the value and every derivative under the same mask, so a branch taken on the value
/// part carries the derivative of the arm taken.
impl<S: Real, const N: usize> Blend<S> for Dual<S, N> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            v: S::select(m, t.v, f.v),
            d: array::from_fn(|i| S::select(m, t.d[i], f.d[i])),
        }
    }
}

impl<S: Real, const N: usize> Real for Dual<S, N> {
    type Mask = S::Mask;
    const PRECISION: Precision = S::PRECISION;

    #[inline]
    fn lit(x: f64) -> Self {
        Self::constant(S::lit(x))
    }
    #[inline]
    fn zero() -> Self {
        Self::constant(S::zero())
    }
    #[inline]
    fn one() -> Self {
        Self::constant(S::one())
    }
    #[inline]
    fn lt(self, rhs: Self) -> S::Mask {
        self.v.lt(rhs.v)
    }
    #[inline]
    fn le(self, rhs: Self) -> S::Mask {
        self.v.le(rhs.v)
    }
    #[inline]
    fn select(m: S::Mask, t: Self, f: Self) -> Self {
        <Self as Blend<S>>::blend(m, t, f)
    }

    /// The value is `S::sqrt` of the value.
    ///
    /// # Domain
    ///
    /// `self.v >= 0` or NaN, as for `S::sqrt`, which checks it in debug builds. The derivative
    /// `d / (2 sqrt v)` needs `self.v > 0`: at `v = 0` it is `+-inf` where `d` is not zero and NaN
    /// (`0 / 0`) where it is, untouched components included; the value is unaffected. Unchecked,
    /// since the value at 0 is legal. The safe-argument pattern keeps a zero out of any selected
    /// arm.
    #[inline]
    fn sqrt(self) -> Self {
        let r = self.v.sqrt();
        let two_r = r + r;
        Self {
            v: r,
            d: self.d.map(|d| d / two_r),
        }
    }

    /// The value is `S::cbrt` of the value, defined for every argument.
    ///
    /// # Domain
    ///
    /// The derivative `d / (3 c^2)`, `c = cbrt v`, needs `self.v != 0`: at `v = +-0` it is `+-inf`
    /// (the sign of `d`, since `c^2` is `+0` for both zeros) where `d` is not zero and NaN
    /// (`0 / 0`) where it is, untouched components included; the value is unaffected. Unchecked,
    /// since the value at 0 is legal. The safe-argument pattern keeps a zero out of any selected
    /// arm. Elsewhere it is accurate: `c^2` is a normal number for every nonzero finite `v` at
    /// both precisions. At `v = +-inf` it is `0` for finite `d`.
    #[inline]
    fn cbrt(self) -> Self {
        let c = self.v.cbrt();
        let three_c2 = S::lit(3.0) * (c * c);
        Self {
            v: c,
            d: self.d.map(|d| d / three_c2),
        }
    }

    #[inline]
    fn sin_cos(self) -> (Self, Self) {
        let (s, c) = self.v.sin_cos();
        (
            Self {
                v: s,
                d: self.d.map(|d| c * d),
            },
            Self {
                v: c,
                d: self.d.map(|d| -(s * d)),
            },
        )
    }

    /// `self` is `y`. The value is `S::atan2` and has no restricted domain.
    ///
    /// # Domain
    ///
    /// The derivative divides by `x.v^2 + y.v^2` and is accurate only where that sum is a normal
    /// number, i.e. the larger of `|x.v|` and `|y.v|` lies in about `1e-154 .. 1e154` (`f64`) or
    /// `1e-19 .. 1e19` (`f32`). Below it the sum is subnormal or zero: the derivative loses digits,
    /// then is `+-inf`, and NaN at the origin (`0 / 0`). Above it the sum is infinite and the
    /// derivative `0`. The value is unaffected; the derivative is unchecked, since the value is
    /// legal everywhere.
    #[inline]
    fn atan2(self, x: Self) -> Self {
        let r2 = x.v * x.v + self.v * self.v;
        Self {
            v: self.v.atan2(x.v),
            d: array::from_fn(|i| (x.v * self.d[i] - self.v * x.d[i]) / r2),
        }
    }

    #[inline]
    fn abs(self) -> Self {
        let neg = self.v.lt(S::zero());
        Self {
            v: self.v.abs(),
            d: self.d.map(|d| S::select(neg, -d, d)),
        }
    }

    /// The derivative is with respect to `self` only; `s` contributes none. The sign taken from
    /// `s` is its sign bit, as in the value: `copysign(x, -0.0)` is `-|x|` and differentiates as
    /// such.
    #[inline]
    fn copysign(self, s: Self) -> Self {
        let neg_x = self.v.lt(S::zero());
        let neg_s = S::one().copysign(s.v).lt(S::zero());
        let flip = neg_x.or(neg_s).and(neg_x.and(neg_s).not());
        Self {
            v: self.v.copysign(s.v),
            d: self.d.map(|d| S::select(flip, -d, d)),
        }
    }

    #[inline]
    fn value_f64(self) -> f64 {
        self.v.value_f64()
    }
}
