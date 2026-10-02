//! The scalar model (`docs/PHASE2.md` §2, decision `0003`).

use core::ops::{Add, Div, Mul, Neg, Sub};

/// A truth value that may be one bool or one per SIMD lane.
///
/// `impl Mask for bool` is the only place in the workspace where a mask reaches an `if`; scalar
/// `Real::select` goes through [`Mask::decide`].
pub trait Mask: Copy {
    /// Lane-wise conjunction.
    fn and(self, o: Self) -> Self;
    /// Lane-wise disjunction.
    fn or(self, o: Self) -> Self;
    /// Lane-wise negation.
    fn not(self) -> Self;
    /// True iff every lane is set. Tests and `debug_assert!` only.
    fn all(self) -> bool;
    /// True iff some lane is set. Tests and `debug_assert!` only.
    fn any(self) -> bool;
    /// The branching policy belongs to the mask: `bool` evaluates one arm; a lane mask
    /// evaluates both and blends.
    fn decide<T>(
        self,
        t: impl FnOnce() -> T,
        f: impl FnOnce() -> T,
        blend: impl FnOnce(Self, T, T) -> T,
    ) -> T;
}

/// The floating-point width a [`Real`] evaluates in; selects generated constants at
/// monomorphization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Precision {
    /// Single precision.
    F32,
    /// Double precision.
    F64,
}

/// A scalar that cannot say "less than": comparisons return a [`Mask`], never a `bool`.
///
/// `Real` has neither `PartialOrd` nor `PartialEq`, so generic code cannot branch on a float. Both
/// of these fail to compile:
///
/// ```compile_fail,E0369
/// use helicoid_linalg::Real;
/// fn less<S: Real>(a: S, b: S) -> bool { a < b }
/// ```
///
/// ```compile_fail,E0369
/// use helicoid_linalg::Real;
/// fn same<S: Real>(a: S, b: S) -> bool { a == b }
/// ```
///
/// Stable rustdoc does not check the error code, so this positive control proves that the two
/// failures above come from the missing bound and not from a broken import:
///
/// ```
/// use helicoid_linalg::Real;
/// fn less<S: Real + PartialOrd>(a: S, b: S) -> bool { a < b }
/// fn same<S: Real + PartialEq>(a: S, b: S) -> bool { a == b }
/// assert!(less(1.0_f64, 2.0) && same(1.0_f32, 1.0));
/// ```
///
/// Use [`Real::lt`]/[`Real::le`] with [`Real::select`] or [`Real::branch`]; the exact arm of a
/// branch evaluates at a safe argument, `S::select(small, S::one(), x)`, so it is finite in every
/// lane.
pub trait Real:
    Copy
    + 'static
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
    /// The mask type comparisons produce.
    type Mask: Mask;
    /// The precision of this scalar.
    const PRECISION: Precision;
    /// Exactly representable constants only (0.5, 2, 3, ...). Anything that must be correctly
    /// rounded per precision is generated (`coeffs::generated`).
    ///
    /// # Domain
    ///
    /// `x` is exactly representable at this precision; the `f32` impl checks it in debug builds.
    /// Release builds round.
    fn lit(x: f64) -> Self;
    /// Additive identity, `+0`.
    fn zero() -> Self;
    /// Multiplicative identity.
    fn one() -> Self;
    /// `self < rhs`.
    fn lt(self, rhs: Self) -> Self::Mask;
    /// `self <= rhs`.
    fn le(self, rhs: Self) -> Self::Mask;
    /// `t` where `m` is set, else `f`.
    fn select(m: Self::Mask, t: Self, f: Self) -> Self;
    /// Evaluate `t` or `f` under the mask's policy: one arm for `bool`, both and a blend for lanes.
    fn branch<T: Blend<Self>>(m: Self::Mask, t: impl FnOnce() -> T, f: impl FnOnce() -> T) -> T {
        m.decide(t, f, |m, a, b| T::blend(m, a, b))
    }
    /// Square root.
    ///
    /// # Domain
    ///
    /// `self >= 0` or NaN; the scalar impls check it in debug builds. Release builds return NaN
    /// for a negative argument. Lane code evaluates both arms, so it passes a safe argument.
    fn sqrt(self) -> Self;
    /// The real cube root: odd, `cbrt(-x) = -cbrt(x)`, so a zero and an infinity keep their sign.
    ///
    /// Every argument is in the domain. The derivative of a `Dual` is singular at 0
    /// ([`Dual::cbrt`](crate::Dual)).
    fn cbrt(self) -> Self;
    /// `(sin self, cos self)`.
    fn sin_cos(self) -> (Self, Self);
    /// The angle of the point `(x, self)`, i.e. `atan2(y = self, x)`.
    fn atan2(self, x: Self) -> Self;
    /// Absolute value.
    fn abs(self) -> Self;
    /// `self` with the sign bit of `sign`. The sign bit of a NaN `sign` is the target's (`0018`).
    fn copysign(self, sign: Self) -> Self;
    /// The value part as `f64`. Tests and `debug_assert!` only (`docs/API.md` R4).
    fn value_f64(self) -> f64;
}

/// Set iff `x` is finite: `x * 0` is `+-0` for a finite `x` and NaN for NaN and `+-inf`.
#[inline]
pub(crate) fn is_finite<S: Real>(x: S) -> S::Mask {
    (x * S::zero()).le(S::zero())
}

/// A value type that a mask can select between, so [`Real::branch`] can return tuples of
/// coefficients and whole group elements.
pub trait Blend<S: Real>: Sized {
    /// `t` where `m` is set, else `f`.
    fn blend(m: S::Mask, t: Self, f: Self) -> Self;
}

impl<S: Real> Blend<S> for S {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        S::select(m, t, f)
    }
}

macro_rules! blend_tuple {
    ($($T:ident $t:ident $f:ident),+) => {
        impl<S: Real, $($T: Blend<S>),+> Blend<S> for ($($T,)+) {
            #[inline]
            fn blend(m: S::Mask, ($($t,)+): Self, ($($f,)+): Self) -> Self {
                ($($T::blend(m, $t, $f),)+)
            }
        }
    };
}

blend_tuple!(A a0 b0);
blend_tuple!(A a0 b0, B a1 b1);
blend_tuple!(A a0 b0, B a1 b1, C a2 b2);
blend_tuple!(A a0 b0, B a1 b1, C a2 b2, D a3 b3);
blend_tuple!(A a0 b0, B a1 b1, C a2 b2, D a3 b3, E a4 b4);
blend_tuple!(A a0 b0, B a1 b1, C a2 b2, D a3 b3, E a4 b4, F a5 b5);
blend_tuple!(A a0 b0, B a1 b1, C a2 b2, D a3 b3, E a4 b4, F a5 b5, G a6 b6);
blend_tuple!(A a0 b0, B a1 b1, C a2 b2, D a3 b3, E a4 b4, F a5 b5, G a6 b6, H a7 b7);

impl<S: Real, T: Blend<S>, const N: usize> Blend<S> for [T; N] {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        let mut f = f.into_iter();
        t.map(|a| match f.next() {
            Some(b) => T::blend(m, a, b),
            None => {
                // `f` has exactly `N` items; release returns `a` rather than panicking (D11).
                debug_assert!(false, "Blend for [T; N]: arrays differ in length");
                a
            }
        })
    }
}
