//! `Mask` for `bool` and `Real` for `f64`/`f32`, every transcendental through `libm` (D16).

use crate::real::{Mask, Precision, Real};

impl Mask for bool {
    #[inline]
    fn and(self, o: Self) -> Self {
        self & o
    }
    #[inline]
    fn or(self, o: Self) -> Self {
        self | o
    }
    #[inline]
    fn not(self) -> Self {
        !self
    }
    #[inline]
    fn all(self) -> bool {
        self
    }
    #[inline]
    fn any(self) -> bool {
        self
    }
    #[inline]
    fn decide<T>(
        self,
        t: impl FnOnce() -> T,
        f: impl FnOnce() -> T,
        _blend: impl FnOnce(Self, T, T) -> T,
    ) -> T {
        if self {
            t()
        } else {
            f()
        }
    }
}

macro_rules! impl_real {
    ($t:ty, $prec:expr, $sqrt:ident, $cbrt:ident, $sincos:ident, $sin:ident, $cos:ident, $acos:ident, $atan2:ident, $fabs:ident, $copysign:ident, $lit:ident) => {
        impl Real for $t {
            type Mask = bool;
            const PRECISION: Precision = $prec;
            #[inline]
            fn lit(x: f64) -> Self {
                $lit(x)
            }
            #[inline]
            fn zero() -> Self {
                0.0
            }
            #[inline]
            fn one() -> Self {
                1.0
            }
            #[inline]
            fn lt(self, rhs: Self) -> bool {
                self < rhs
            }
            #[inline]
            fn le(self, rhs: Self) -> bool {
                self <= rhs
            }
            #[inline]
            fn select(m: bool, t: Self, f: Self) -> Self {
                m.decide(|| t, || f, |_, a, _| a)
            }
            #[inline]
            fn sqrt(self) -> Self {
                debug_assert!(!Real::lt(self, 0.0), "Real::sqrt: negative argument");
                libm::$sqrt(self)
            }
            #[inline]
            fn cbrt(self) -> Self {
                libm::$cbrt(self)
            }
            #[inline]
            fn sin_cos(self) -> (Self, Self) {
                libm::$sincos(self)
            }
            #[inline]
            fn sin(self) -> Self {
                libm::$sin(self)
            }
            #[inline]
            fn cos(self) -> Self {
                libm::$cos(self)
            }
            #[inline]
            fn acos(self) -> Self {
                libm::$acos(self)
            }
            #[inline]
            fn atan2(self, x: Self) -> Self {
                libm::$atan2(self, x)
            }
            #[inline]
            fn abs(self) -> Self {
                libm::$fabs(self)
            }
            #[inline]
            fn copysign(self, sign: Self) -> Self {
                libm::$copysign(self, sign)
            }
            #[inline]
            fn value_f64(self) -> f64 {
                f64::from(self)
            }
        }
    };
}

fn lit_f64(x: f64) -> f64 {
    x
}

fn lit_f32(x: f64) -> f32 {
    let y = x as f32;
    debug_assert!(
        f64::from(y).to_bits() == x.to_bits(),
        "Real::lit takes exactly representable constants"
    );
    y
}

impl_real!(
    f64,
    Precision::F64,
    sqrt,
    cbrt,
    sincos,
    sin,
    cos,
    acos,
    atan2,
    fabs,
    copysign,
    lit_f64
);
impl_real!(
    f32,
    Precision::F32,
    sqrtf,
    cbrtf,
    sincosf,
    sinf,
    cosf,
    acosf,
    atan2f,
    fabsf,
    copysignf,
    lit_f32
);
