//! `hat` and `vee` between `Vec3` and skew `Mat3` (`docs/NUMERICS.md` §1: `varphi^ = [varphi]_x`).

use crate::matrix::Mat3;
use crate::real::Real;
use crate::vector::{Vec3, Vector};

/// The cross-product matrix `[v]_x`, so that `hat(a) * b == a.cross(b)` for finite `b`
/// (the structural zeros of `hat` give `0 * inf = NaN` where the cross product has no such term).
///
/// ```text
///          [  0  -z   y ]
/// [v]_x =  [  z   0  -x ]
///          [ -y   x   0 ]
/// ```
#[inline]
pub fn hat<S: Real>(v: Vec3<S>) -> Mat3<S> {
    let [x, y, z] = v.0;
    let o = S::zero();
    Mat3::from_rows([Vector([o, -z, y]), Vector([z, o, -x]), Vector([-y, x, o])])
}

/// The inverse of [`hat`]: reads the entries `(2, 1)`, `(0, 2)`, `(1, 0)`.
///
/// It does not check that `m` is skew and does not symmetrize: for a general `m` the result is
/// those three entries, not the skew part `(m - m^T)/2`.
#[inline]
pub fn vee<S: Real>(m: Mat3<S>) -> Vec3<S> {
    Vector([m.get(2, 1), m.get(0, 2), m.get(1, 0)])
}
