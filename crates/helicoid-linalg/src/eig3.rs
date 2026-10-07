//! Eigendecomposition of a symmetric 3x3 matrix (`docs/PHASE2.md` §6, decision `0017`).
//!
//! Symmetric 3x3 eigendecomposition generic over `S: Real`: Smith's closed form
//! for eigenvalues and cross products of rows of `A - lambda I` for eigenvectors.
//! Every decision is a mask; square roots and divisions use safe arguments. The frame
//! is always orthonormal.

use crate::cubic::pi;
use crate::matrix::Mat3;
use crate::real::{is_finite, Blend, Mask, Real};
use crate::vector::{Vec3, Vector};

/// The unit vector along `v`, or `fallback` where `v` is zero or has an entry that is not finite.
/// `v` is divided by its largest entry first, so no square under- or overflows and the result has
/// unit length at every magnitude.
#[inline]
fn unit_or<S: Real>(v: Vec3<S>, fallback: Vec3<S>) -> Vec3<S> {
    let [x, y, z] = v.0.map(Real::abs);
    let big = S::select(x.lt(y), y, x);
    let big = S::select(big.lt(z), z, big);
    let ok = S::zero()
        .lt(big)
        .and(is_finite(x))
        .and(is_finite(y))
        .and(is_finite(z));
    let big = S::select(ok, big, S::one());
    let w = Vector(v.0.map(|c| c / big));
    let n = w.norm();
    Vec3::blend(ok, Vector(w.0.map(|c| c / n)), fallback)
}

/// A unit vector orthogonal to the unit vector `v`: `v x e`, `e` the axis of its smallest
/// component, so `|v x e|^2 = 1 - v_e^2 >= 2/3`. The first axis on a tie.
#[inline]
fn any_orthogonal<S: Real>(v: Vec3<S>) -> Vec3<S> {
    let [x, y, z] = v.0.map(Real::abs);
    let (zero, one) = (S::zero(), S::one());
    let e = Vec3::blend(
        x.le(y).and(x.le(z)),
        Vector([one, zero, zero]),
        Vec3::blend(
            y.le(z),
            Vector([zero, one, zero]),
            Vector([zero, zero, one]),
        ),
    );
    unit_or(v.cross(e), e)
}

/// The longest cross product of two rows of `A - lambda I`; `m` is `[xx, yy, zz, xy, xz, yz]`. For
/// a simple eigenvalue `lambda` the rows span a plane, and the cross product of two of them is an
/// eigenvector.
#[inline]
fn null_vec<S: Real>(m: [S; 6], lambda: S) -> Vec3<S> {
    let [xx, yy, zz, xy, xz, yz] = m;
    let r0 = Vector([xx - lambda, xy, xz]);
    let r1 = Vector([xy, yy - lambda, yz]);
    let r2 = Vector([xz, yz, zz - lambda]);
    let (c01, c02, c12) = (r0.cross(r1), r0.cross(r2), r1.cross(r2));
    let (n01, n02, n12) = (c01.norm_sq(), c02.norm_sq(), c12.norm_sq());
    let first = n02.le(n01).and(n12.le(n01));
    Vec3::blend(first, c01, Vec3::blend(n12.le(n02), c02, c12))
}

/// The eigenvalues, ascending, and an orthonormal basis of eigenvectors of the symmetric matrix `a`.
///
/// Returns `(values, vectors)`: `values.0[0] <= values.0[1] <= values.0[2]`, column `k` of `vectors`
/// an eigenvector of `values.0[k]`, `vectors` orthonormal with `det = +1`. The sign of a column is
/// what the recovery produces; only the parity of the three is fixed.
///
/// With `q = trace / 3`, `B = A - q I` and `p = sqrt(|B|_F^2 / 6)` the eigenvalues are Smith's
/// `q + 2 p cos((acos(r) + 2 pi k) / 3)`, `k = 0, 1, 2`, `r = det(B) / (2 p^3)` clamped to
/// `[-1, 1]`, sorted by a compare-exchange network (`acos` and `pi` as in
/// [`solve_cubic`](crate::solve_cubic)). For the eigenvalues `l0 <= l1 <= l2`, the vector of `l2` is
/// the longest cross product of two rows of `A - l2 I`, normalised (`v2`); that of `l1` is the same
/// for `l1`, projected orthogonal to `v2` and normalised (`v1`); the third is `v1 x v2`. Where a
/// cross product is zero, `v2` is `e_z` and `v1` a fixed unit vector orthogonal to `v2`: a multiple
/// of the identity gives the identity.
///
/// # Domain
///
/// Every input is legal and nothing is asserted (`docs/API.md` R6, `0023` (draft)): only the lower
/// triangle is read, and the caller owns symmetry. The frame is finite and orthonormal to a few `u`
/// for every input; the eigenvalues are not finite (`NaN` or infinite) for a non-finite entry of the
/// lower triangle, never a plausible number. For entries of magnitude `m`, `p^3` must be a normal
/// number: about `1e-100 < m < 1e100` (`f64`) and `1e-12 < m < 1e12` (`f32`). Digits are lost while
/// it is subnormal (`1e-8 |A|` at `2^-350`, `f64`), and the eigenvalues are not finite where it
/// overflows or underflows to 0, unless the spread of `A - q I` is below the rounding of `q`, when
/// they are `q`. The vectors
/// need the squared length of a cross product, `|A|^4`, normal: `1e-75 < m < 1e75` and
/// `1e-9 < m < 1e9`.
///
/// **Near a double eigenvalue the closed form is not backward stable.** `r` sits near `+-1`, where
/// `acos` has an infinite slope: an error `u` in `r` is `sqrt(u)` in the angle, and the two close
/// eigenvalues are off by up to `sqrt(u) |A|` (`6e-9 |A|` on a rank-1 matrix at `f64`). A cross
/// product is a difference of products of size `|A|^2` that cancels to `|A| gap`, so its rounding
/// error tilts the vector toward the third by `u |A| / gap`. With `gap` the smaller gap between
/// adjacent eigenvalues and `c = 1 + min(|A| / gap, u^(-1/2))`, the recorded bounds (`eig3_tests`,
/// a fit to the measured errors) are `(4 c + 24) u |A|` for the eigenvalue error and the residual
/// `|A v - l v|`, and `(4 c + 24) u |A| / gap_k` for the angle of `v_k`, `gap_k` its own gap. The
/// bars of `docs/PHASE2.md` §6 name no constant, so none is claimed met; Kopp's hybrid is owed.
///
/// **Where the two largest eigenvalues are within about `2 sqrt(u) |A|` of each other, equal ones
/// included, no column is reliable.** The eigenvalues are right (the pair to `sqrt(u) |A|`), but
/// `v2` is rounding noise or the fallback `e_z`, and `v1` and `v0` are built from it: the vector of
/// the isolated `l0` has a residual of the order of `|A|` however wide its gap (`0.8 |A|` for
/// `diag(5, 1, 5)`).
///
/// A `Dual` result differentiates the arm taken: the derivative is `NaN` where `acos` sees
/// `r = +-1` (a double eigenvalue), and for a multiple of the identity the arms are constants.
///
/// # Example
///
/// ```
/// use helicoid_linalg::{eig3, Mat3, Vector};
///
/// // [[2, 1, 0], [1, 2, 0], [0, 0, 5]] has the eigenvalues 1, 3, 5.
/// let a = Mat3::from_rows([
///     Vector([2.0_f64, 1.0, 0.0]),
///     Vector([1.0, 2.0, 0.0]),
///     Vector([0.0, 0.0, 5.0]),
/// ]);
/// let (values, vectors) = eig3(&a);
/// assert!(values.0.iter().zip([1.0, 3.0, 5.0]).all(|(l, w)| (l - w).abs() < 1e-12));
/// assert!((vectors.get(2, 2).abs() - 1.0).abs() < 1e-12); // the largest: the z axis
/// ```
#[inline]
pub fn eig3<S: Real>(a: &Mat3<S>) -> (Vec3<S>, Mat3<S>) {
    let (zero, one) = (S::zero(), S::one());
    let (two, three) = (S::lit(2.0), S::lit(3.0));
    let m = [
        a.get(0, 0),
        a.get(1, 1),
        a.get(2, 2),
        a.get(1, 0),
        a.get(2, 0),
        a.get(2, 1),
    ];
    let [xx, yy, zz, xy, xz, yz] = m;

    let q = (xx + yy + zz) / three;
    let (bxx, byy, bzz) = (xx - q, yy - q, zz - q);
    let sum_sq = bxx * bxx + byy * byy + bzz * bzz + two * (xy * xy + xz * xz + yz * yz);
    let p_sq = sum_sq / S::lit(6.0);
    // Each degenerate arm is taken on a test that `NaN` fails, so `NaN` reaches the general arm and
    // propagates. `sqrt` is singular at 0 for a `Dual`; a multiple of the identity has `p = 0` as a
    // constant.
    let p = S::select(p_sq.le(zero), zero, p_sq.sqrt());
    let p_cubed = p_sq * p;
    let det_b = bxx * (byy * bzz - yz * yz) - xy * (xy * bzz - yz * xz) + xz * (xy * yz - byy * xz);
    // `p^3 = 0` is a multiple of the identity (`r = 0`, `l = q`) or `p^2` or `p^3` underflowed. Then
    // `r = 0` is right only where the spread of `B` is below the rounding of `q`; otherwise `q`
    // three times would be a wrong answer that looks right, so it is `NaN`.
    let no_p3 = p_cubed.le(zero);
    let r = S::branch(
        no_p3,
        || {
            let spread = bxx.abs() + byy.abs() + bzz.abs() + xy.abs() + xz.abs() + yz.abs();
            let rounded = q.abs();
            S::select((rounded + spread).le(rounded), zero, S::zero() / S::zero())
        },
        || det_b / (two * S::select(no_p3, one, p_cubed)),
    );
    let r = S::select(r.lt(-one), -one, r);
    let r = S::select(one.lt(r), one, r);

    let theta = r.acos();
    let two_pi = two * pi::<S>();
    let lambda = |k: usize| {
        let phi = (theta + two_pi * S::lit(k as f64)) / three;
        q + two * p * phi.cos()
    };
    let order = |x: S, y: S| {
        let swap = y.lt(x);
        (S::select(swap, y, x), S::select(swap, x, y))
    };
    let (l0, l1) = order(lambda(0), lambda(1));
    let (l1, l2) = order(l1, lambda(2));
    let (l0, l1) = order(l0, l1);

    let e_z = Vector([zero, zero, one]);
    let v2 = unit_or(null_vec(m, l2), e_z);
    let w = any_orthogonal(v2);
    let x = v2.cross(w);
    let u1 = unit_or(null_vec(m, l1), w);
    let in_plane = w.scale(w.dot(u1)) + x.scale(x.dot(u1));
    let v1 = unit_or(in_plane, w);
    let v0 = unit_or(v1.cross(v2), x);
    (Vector([l0, l1, l2]), Mat3::from_cols([v0, v1, v2]))
}
