//! [`Libm`], the shipped scalar with `atan2` alone `libm`'s: the arrangement `0105` replaced.
//!
//! Self-contained (nothing from `crate::`) so `xtask/benches/atan2.rs` compiles the same file.

use helicoid_linalg::{Precision, Real};

/// The shipped scalar `T` with `atan2` alone `libm`'s at `T`'s precision.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub(crate) struct Libm<T>(pub(crate) T);

macro_rules! libm_scalar {
    ($t:ty, $atan2:path) => {
        impl core::ops::Add for Libm<$t> {
            type Output = Self;
            fn add(self, o: Self) -> Self {
                Libm(self.0 + o.0)
            }
        }
        impl core::ops::Sub for Libm<$t> {
            type Output = Self;
            fn sub(self, o: Self) -> Self {
                Libm(self.0 - o.0)
            }
        }
        impl core::ops::Mul for Libm<$t> {
            type Output = Self;
            fn mul(self, o: Self) -> Self {
                Libm(self.0 * o.0)
            }
        }
        impl core::ops::Div for Libm<$t> {
            type Output = Self;
            fn div(self, o: Self) -> Self {
                Libm(self.0 / o.0)
            }
        }
        impl core::ops::Neg for Libm<$t> {
            type Output = Self;
            fn neg(self) -> Self {
                Libm(-self.0)
            }
        }
        impl From<Libm<$t>> for f64 {
            fn from(x: Libm<$t>) -> f64 {
                f64::from(x.0)
            }
        }

        // Every method but `atan2` is the shipped scalar's, so the twin changes one variable.
        impl Real for Libm<$t> {
            type Mask = bool;
            const PRECISION: Precision = <$t as Real>::PRECISION;
            fn lit(x: f64) -> Self {
                Libm(<$t as Real>::lit(x))
            }
            fn round_from_f64(x: f64) -> Self {
                Libm(<$t as Real>::round_from_f64(x))
            }
            fn zero() -> Self {
                Libm(<$t as Real>::zero())
            }
            fn one() -> Self {
                Libm(<$t as Real>::one())
            }
            fn lt(self, rhs: Self) -> bool {
                <$t as Real>::lt(self.0, rhs.0)
            }
            fn le(self, rhs: Self) -> bool {
                <$t as Real>::le(self.0, rhs.0)
            }
            fn select(m: bool, t: Self, f: Self) -> Self {
                Libm(<$t as Real>::select(m, t.0, f.0))
            }
            fn sqrt(self) -> Self {
                Libm(Real::sqrt(self.0))
            }
            fn cbrt(self) -> Self {
                Libm(Real::cbrt(self.0))
            }
            fn sin_cos(self) -> (Self, Self) {
                let (s, c) = Real::sin_cos(self.0);
                (Libm(s), Libm(c))
            }
            fn sin(self) -> Self {
                Libm(Real::sin(self.0))
            }
            fn cos(self) -> Self {
                Libm(Real::cos(self.0))
            }
            fn acos(self) -> Self {
                Libm(Real::acos(self.0))
            }
            fn exp(self) -> Self {
                Libm(Real::exp(self.0))
            }
            fn exp_m1(self) -> Self {
                Libm(Real::exp_m1(self.0))
            }
            /// The variable: `libm`'s.
            fn atan2(self, x: Self) -> Self {
                Libm($atan2(self.0, x.0))
            }
            fn abs(self) -> Self {
                Libm(Real::abs(self.0))
            }
            fn copysign(self, sign: Self) -> Self {
                Libm(Real::copysign(self.0, sign.0))
            }
            fn value_f64(self) -> f64 {
                Real::value_f64(self.0)
            }
        }
    };
}

libm_scalar!(f64, libm::atan2);
libm_scalar!(f32, libm::atan2f);

impl From<f64> for Libm<f64> {
    fn from(x: f64) -> Self {
        Libm(x)
    }
}
