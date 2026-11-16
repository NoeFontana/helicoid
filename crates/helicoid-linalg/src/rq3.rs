//! The RQ and LQ factorizations of a 3x3 matrix (`docs/PHASE2.md` §6, `docs/NUMERICS.md` §22,
//! decision `0117`).
//!
//! Gram–Schmidt from the last row, the projection taken twice under Kahan's test, the first row
//! the cross product; each row scaled by its own power of two. Generic over `S: Real`: every
//! decision is a mask, every square root and division has a safe argument. `lq3` is `rq3` of the
//! reversed matrix, reversed.

use crate::eig3::{into_window, smallest_normal, unit, window};
use crate::matrix::Mat3;
use crate::real::{Blend, Mask, Real};
use crate::tri::{LowerTri3, UpperTri3};
use crate::vector::{Vec3, Vector};

/// `v / d`, entry by entry: a division, not a product by `1 / d`, so that `d / d = 1`
/// (`maths/rq3.md` RQ.8).
#[inline(always)]
fn over<S: Real>(v: Vec3<S>, d: S) -> Vec3<S> {
    Vector([v.0[0] / d, v.0[1] / d, v.0[2] / d])
}

fn axis<S: Real>(i: usize) -> Vec3<S> {
    let mut e = [S::zero(); 3];
    e[i] = S::one();
    Vector(e)
}

/// A unit vector orthogonal to the unit `r`: Gram–Schmidt of the axis of its smallest
/// `|component|`, the **last** of equals, whose residual has squared norm at least `2/3`. `e_2` for
/// `r = e_3`, so the zero matrix factors with `R = I` (`0117` decision 2).
#[cold]
#[inline(never)]
fn completion<S: Real>(r0: S, r1: S, r2: S) -> Vec3<S> {
    let r = Vector([r0, r1, r2]);
    let [x, y, z] = r.0.map(|c| c.abs());
    // `z` unless another is strictly smaller; between `x` and `y`, `y` unless `x` is.
    let not_z = x.lt(z).or(y.lt(z));
    let e = Vec3::blend(not_z, Vec3::blend(x.lt(y), axis(0), axis(1)), axis(2));
    let w = e - r.scale(r.dot(e));
    unit(over(w, w.norm()))
}

/// Row `v` brought into `eig3`'s window by its own power of two: `(row, finite, back)`, the row
/// zero where an entry is not finite, `back` the factor that undoes the scaling.
#[inline(always)]
fn scaled<S: Real>(v: Vec3<S>) -> (Vec3<S>, S::Mask, S) {
    let (e, finite, _, back) = into_window(v.0, window::<S>());
    (Vector(e), finite, back)
}

/// `(U, R, valid)` with `m = U R`: `U` upper triangular, `R` a rotation.
///
/// Gram–Schmidt from the last row (`docs/NUMERICS.md` §22, `0117`):
/// - `u33 = |m3|`, `r3 = m3 / u33`;
/// - `r2` the residual of `m2` against `r3`, the projection taken twice, of norm `u22`; where the
///   second projection keeps less than half the first's squared norm, `m2` lies in the span of
///   `r3` to working precision (Kahan's test, `maths/rq3.md` RQ.4) and `u22 = 0`;
/// - `r1 = r2 x r3`, so `R` is a rotation by construction;
/// - `u1j = m1 . rj` and `u23 = m2 . r3`.
///
/// Each row is first scaled by its own exact power of two into a window (`[2^-128, 2^128)` at
/// `f64`, `[2^-32, 2^32)` at `f32`), which leaves `R` unchanged, and the rows of `U` are scaled
/// back.
///
/// `u22` and `u33` are never negative and `u11` has the sign of `det m`. Where `m2` and `m3` are
/// independent this is the one factorization with `R` in SO(3) and `u22, u33 > 0`, whatever the
/// sign of `det m`. `valid` is `u11 > 0 && u22 > 0 && u33 > 0`, read on the computed values, as
/// `chol`'s mask is; in exact arithmetic it is `det m > 0`, where the factorization has a positive
/// diagonal and is unique. Under a clear mask the factorization is still returned.
///
/// # Exactness
///
/// An upper triangular `m` with a positive diagonal, `+0` below it and no `-0` returns `U = m`
/// and `R = I` bit for bit where each row's entries lie within `2^894` (`f64`), `2^94` (`f32`) of
/// its largest and `m22` within `2^383`, `2^31` of its row's largest: a diagonal one at every
/// finite magnitude, subnormal entries included (`maths/rq3.md` RQ.8, RQ.9). A `-0` off the
/// diagonal returns the same values, the zeros of `U` and `R` of either sign.
///
/// # Accuracy
///
/// Backward stable row by row: `m_i = (U R)_i + e_i` with `|e_i|` a few `u |m_i|`, and `R`
/// orthonormal to a few `u`, however close `m2` and `m3` are to parallel (`maths/rq3.md` RQ.4,
/// RQ.5).
///
/// # Domain
///
/// Every input is legal and nothing is asserted (`docs/API.md` R6). For every finite input the
/// outputs are finite, unless an entry of `U` exceeds the largest finite number. Where `m2` and
/// `m3` are parallel, or `m3` is zero, `u22` (and `u33`) is `0` and `R` is completed by a fixed
/// rule; the zero matrix gives `U = 0`, `R = I`. A non-finite entry gives `U = I`, `R = I` and a
/// clear mask, never a plausible factorization.
///
/// A `Dual` result differentiates the construction: wherever the factorization is unique it is the
/// derivative of the factorization (`maths/rq3.md` RQ.6), triangular input included. Where `m2`
/// and `m3` are dependent the lanes are finite and in general no directional derivative.
///
/// # Example
///
/// ```
/// use helicoid_linalg::{rq3, Mat3, Vector};
///
/// // An upper triangle with a positive diagonal is its own factor, bit for bit.
/// let m = Mat3::from_rows([
///     Vector([2.0_f64, 0.5, -1.0]),
///     Vector([0.0, 3.0, 0.25]),
///     Vector([0.0, 0.0, 0.7]),
/// ]);
/// let (u, r, valid) = rq3(&m);
/// assert!(valid);
/// assert_eq!(u.to_row_major(), [2.0, 0.5, -1.0, 3.0, 0.25, 0.7]);
/// assert_eq!(r.get(0, 0), 1.0);
/// assert_eq!(r.get(1, 0), 0.0);
/// ```
#[inline]
pub fn rq3<S: Real>(m: &Mat3<S>) -> (UpperTri3<S>, Mat3<S>, S::Mask) {
    let (zero, one) = (S::zero(), S::one());
    let small = smallest_normal::<S>();
    let (m1, f1, b1) = scaled(m.row(0));
    let (m2, f2, b2) = scaled(m.row(1));
    let (m3, f3, b3) = scaled(m.row(2));
    let finite = f1.and(f2).and(f3);

    // Row 3. A row below the smallest normal after scaling is the zero row.
    let a3 = m3.norm_sq();
    let null3 = a3.lt(small);
    let s3 = S::select(null3, one, a3).sqrt();
    let r3 = Vec3::blend(null3, axis(2), unit(over(m3, s3)));
    let u33 = S::select(null3, zero, s3);

    // Row 2: the projection twice, kept where the second leaves at least half of the first's
    // squared norm (`maths/rq3.md` RQ.4); otherwise `m2` is in `r3`'s span to working precision.
    let u23 = m2.dot(r3);
    let w1 = m2 - r3.scale(u23);
    let w2 = w1 - r3.scale(w1.dot(r3));
    let a2 = w2.norm_sq();
    let dep = a2.lt(small).or(a2.lt(S::lit(0.5) * w1.norm_sq()));
    let s2 = S::select(dep, one, a2).sqrt();
    let r2 = S::branch(
        dep,
        move || completion(r3.0[0], r3.0[1], r3.0[2]),
        move || unit(over(w2, s2)),
    );
    let u22 = S::select(dep, zero, s2);

    // Row 1: a rotation by construction.
    let r1 = unit(r2.cross(r3));
    let (u11, u22, u33) = (m1.dot(r1) * b1, u22 * b2, u33 * b3);
    // On the returned diagonal, not the scaled one: a diagonal entry that underflows on the way
    // back clears the mask rather than contradict it.
    let valid = finite.and(zero.lt(u11)).and(zero.lt(u22)).and(zero.lt(u33));
    let u = UpperTri3::from_row_major([u11, m1.dot(r2) * b1, m1.dot(r3) * b1, u22, u23 * b2, u33]);
    let r = Mat3::from_rows([r1, r2, r3]);
    let (u, r) =
        <(UpperTri3<S>, Mat3<S>)>::blend(finite, (u, r), (UpperTri3::identity(), Mat3::identity()));
    (u, r, valid)
}

/// `J m J`, `J` the reversal: entry `(i, j)` from `(2 - i, 2 - j)`.
#[inline(always)]
fn flip<S: Real>(m: &Mat3<S>) -> Mat3<S> {
    Mat3::from_rows([2, 1, 0].map(|i| Vector([m.get(i, 2), m.get(i, 1), m.get(i, 0)])))
}

/// `(L, Q, valid)` with `m = L Q`: `L` lower triangular, `Q` a rotation.
///
/// `rq3` of `J m J` reversed: `L = J U J`, `Q = J R J` (`J` the reversal, `maths/rq3.md` RQ.2), bit
/// for bit, so Gram–Schmidt runs from the first row. `l11` and `l22` are never negative and `l33`
/// has the sign of `det m`; `valid` is `l11 > 0 && l22 > 0 && l33 > 0`, in exact arithmetic
/// `det m > 0`. Every statement of `rq3` holds with the rows reversed: a lower triangle with a
/// positive diagonal returns `L = m`, `Q = I` bit for bit.
///
/// # Domain
///
/// Every input is legal and nothing is asserted (`docs/API.md` R6), as `rq3`.
///
/// # Example
///
/// ```
/// use helicoid_linalg::{lq3, Mat3, Vector};
///
/// let m = Mat3::from_rows([
///     Vector([1.0_f64, 2.0, 0.0]),
///     Vector([0.0, 1.0, 1.0]),
///     Vector([1.0, 0.0, 1.0]),
/// ]);
/// let (l, q, valid) = lq3(&m);
/// assert!(valid);
/// let back = l.to_matrix() * q;
/// assert!((0..9).all(|k| (back.get(k / 3, k % 3) - m.get(k / 3, k % 3)).abs() < 1e-14));
/// ```
#[inline]
pub fn lq3<S: Real>(m: &Mat3<S>) -> (LowerTri3<S>, Mat3<S>, S::Mask) {
    let (u, r, valid) = rq3(&flip(m));
    let [a, b, c, d, e, f] = u.to_row_major();
    (
        LowerTri3::from_row_major([f, e, d, c, b, a]),
        flip(&r),
        valid,
    )
}
