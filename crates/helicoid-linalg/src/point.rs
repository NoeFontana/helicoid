//! Points: an affine space over [`Vector`] (`docs/PHASE2.md` §4).
//!
//! `Point - Point` is a [`Vector`] and `Point + Vector` is a `Point`. There is no `Point + Point`
//! and no scaling: neither is meaningful on an affine space.

use crate::real::{Blend, Real};
use crate::vector::Vector;
use core::array;
use core::ops::{Add, Sub};

/// A position in `N`-space. Layout is not a contract (D2).
///
/// There is no `PartialEq` and no `Point + Point`. Neither of these compiles:
///
/// ```compile_fail,E0369
/// use helicoid_linalg::Point;
/// fn same(a: Point<f64, 3>, b: Point<f64, 3>) -> bool { a == b }
/// ```
///
/// ```compile_fail,E0369
/// use helicoid_linalg::Point;
/// fn sum(a: Point<f64, 3>, b: Point<f64, 3>) -> Point<f64, 3> { a + b }
/// ```
///
/// Positive control: the affine operations that do exist.
///
/// ```
/// use helicoid_linalg::{Point, Vector};
/// fn shifted(a: Point<f64, 3>, b: Point<f64, 3>) -> Point<f64, 3> { a + (b - a) }
/// let q = shifted(Point([1.0; 3]), Point([2.0; 3]));
/// assert_eq!(q.0[2].to_bits(), 2.0_f64.to_bits());
/// let _: Vector<f64, 3> = Point([2.0; 3]) - Point([1.0; 3]);
/// ```
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct Point<S, const N: usize>(pub [S; N]);

/// A 2-D point.
pub type Point2<S> = Point<S, 2>;
/// A 3-D point.
pub type Point3<S> = Point<S, 3>;

impl<S: Real, const N: usize> Sub for Point<S, N> {
    type Output = Vector<S, N>;
    #[inline]
    fn sub(self, o: Self) -> Vector<S, N> {
        Vector(array::from_fn(|i| self.0[i] - o.0[i]))
    }
}

impl<S: Real, const N: usize> Add<Vector<S, N>> for Point<S, N> {
    type Output = Self;
    #[inline]
    fn add(self, v: Vector<S, N>) -> Self {
        Self(array::from_fn(|i| self.0[i] + v.0[i]))
    }
}

impl<S: Real, const N: usize> Blend<S> for Point<S, N> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self(<[S; N]>::blend(m, t.0, f.0))
    }
}
