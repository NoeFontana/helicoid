//! Fixed-size column vectors (`docs/PHASE2.md` §4).
//!
//! Every reduction runs left to right from the first term, so a result is a fixed sequence of
//! `S` operations (D16) and a lone `-0` stays `-0`. Nothing here compares a scalar.

use crate::real::{Blend, Real};
use core::array;
use core::ops::{Add, Neg, Sub};

/// A vector of `N` scalars. Layout is not a contract (D2).
///
/// There is no `PartialEq`: a float comparison is not the equality of two values. This does not
/// compile:
///
/// ```compile_fail,E0369
/// use helicoid_linalg::Vector;
/// fn same(a: Vector<f64, 3>, b: Vector<f64, 3>) -> bool { a == b }
/// ```
///
/// Positive control, so the failure above comes from the missing impl and not from an import:
///
/// ```
/// use helicoid_linalg::Vector;
/// fn sum(a: Vector<f64, 3>, b: Vector<f64, 3>) -> Vector<f64, 3> { a + b }
/// assert_eq!(sum(Vector([1.0; 3]), Vector([2.0; 3])).0[0].to_bits(), 3.0_f64.to_bits());
/// ```
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct Vector<S, const N: usize>(pub [S; N]);

/// A 2-vector.
pub type Vec2<S> = Vector<S, 2>;
/// A 3-vector.
pub type Vec3<S> = Vector<S, 3>;

/// The sum of `terms` added left to right starting from the first, `+0` when there are none.
#[inline]
pub(crate) fn sum<S: Real>(mut terms: impl Iterator<Item = S>) -> S {
    match terms.next() {
        Some(first) => terms.fold(first, |acc, t| acc + t),
        None => S::zero(),
    }
}

impl<S: Real, const N: usize> Vector<S, N> {
    /// `s * self`, componentwise.
    #[inline]
    pub fn scale(self, s: S) -> Self {
        Self(self.0.map(|x| x * s))
    }

    /// `sum_i self_i o_i`, accumulated in index order.
    #[inline]
    pub fn dot(self, o: Self) -> S {
        sum(self.0.into_iter().zip(o.0).map(|(a, b)| a * b))
    }

    /// `self . self`, the form every squared norm in the workspace takes (`NUMERICS.md` §2.1).
    #[inline]
    pub fn norm_sq(self) -> S {
        self.dot(self)
    }

    /// `sqrt(self . self)`; `S::sqrt` only.
    ///
    /// # Domain
    ///
    /// `norm_sq` must be a normal finite number: it is not scaled, so a finite vector whose
    /// entries are beyond about `1e154` (`f64`) or `1e19` (`f32`) returns `inf`, and one whose
    /// largest entry is below about `1e-154` or `1e-19` returns `0` or loses digits. Over `Dual`
    /// the derivative at the zero vector is NaN, the `sqrt` rule at 0 (`Dual::sqrt`, decided in
    /// `0020`); a caller that can meet it uses `norm_sq` behind the safe-argument pattern. The
    /// guarded form is silent below the domain: for entries under about `1e-154` it returns
    /// `d = 0` where the true derivative is a unit vector.
    ///
    /// ```
    /// use helicoid_linalg::{Dual, Real, Vector};
    ///
    /// fn guarded_norm<S: Real, const N: usize>(v: Vector<S, N>) -> S {
    ///     let n2 = v.norm_sq();
    ///     let zero = n2.le(S::zero());
    ///     S::select(zero, S::zero(), S::select(zero, S::one(), n2).sqrt())
    /// }
    /// let var = |x: [f64; 3]| {
    ///     Vector::<Dual<f64, 3>, 3>(core::array::from_fn(|i| Dual::variable(x[i], i)))
    /// };
    /// assert!(var([0.0; 3]).norm().d[0].is_nan());
    /// assert_eq!(guarded_norm(var([0.0; 3])).d, [0.0; 3]);
    /// let n = guarded_norm(var([3.0, 4.0, 0.0]));
    /// assert_eq!((n.v, n.d), (5.0, [0.6, 0.8, 0.0]));
    /// ```
    #[inline]
    pub fn norm(self) -> S {
        self.norm_sq().sqrt()
    }
}

impl<S: Real> Vector<S, 3> {
    /// `self x o`: `(a1 b2 - a2 b1, a2 b0 - a0 b2, a0 b1 - a1 b0)`, right-handed.
    #[inline]
    pub fn cross(self, o: Self) -> Self {
        let ([a0, a1, a2], [b0, b1, b2]) = (self.0, o.0);
        Self([a1 * b2 - a2 * b1, a2 * b0 - a0 * b2, a0 * b1 - a1 * b0])
    }
}

impl<S: Real, const N: usize> Add for Vector<S, N> {
    type Output = Self;
    #[inline]
    fn add(self, o: Self) -> Self {
        Self(array::from_fn(|i| self.0[i] + o.0[i]))
    }
}

impl<S: Real, const N: usize> Sub for Vector<S, N> {
    type Output = Self;
    #[inline]
    fn sub(self, o: Self) -> Self {
        Self(array::from_fn(|i| self.0[i] - o.0[i]))
    }
}

impl<S: Real, const N: usize> Neg for Vector<S, N> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self(self.0.map(Neg::neg))
    }
}

impl<S: Real, const N: usize> Blend<S> for Vector<S, N> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self(<[S; N]>::blend(m, t.0, f.0))
    }
}
