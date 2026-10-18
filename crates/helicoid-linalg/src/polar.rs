//! The rotation nearest to a 3x3 matrix (`docs/NUMERICS.md` §19, decision `0102`, draft).
//!
//! Newton's polar iteration where the nearest rotation is the orthogonal polar factor and the
//! determinant resolves it, a reflection through `svd3`'s last right singular vector where the
//! determinant is negative, and `svd3`'s `U V^T` between; generic over `S: Real`.

use crate::eig3::{eig3, ladder, largest};
use crate::matrix::Mat3;
use crate::real::{is_finite, pow2, Mask, Precision, Real};
use crate::svd3::{svd3, WINDOW_F32, WINDOW_F64};
use crate::vector::{Vec3, Vector};

/// At most eight Newton steps: from `|det M| >= tau |M|_F^3` the scaled iteration converged in at
/// most seven on every family of `0102`'s study.
const STEPS: usize = 8;

/// The columns of the cofactor matrix, `X^-T = C / det X`, and `det X`.
#[inline(always)]
fn cofactors<S: Real>([x0, x1, x2]: [Vec3<S>; 3]) -> ([Vec3<S>; 3], S) {
    let c = [x1.cross(x2), x2.cross(x0), x0.cross(x1)];
    let det = x0.dot(c[0]);
    (c, det)
}

fn norm_sq<S: Real>(x: &[Vec3<S>; 3]) -> S {
    x[0].norm_sq() + x[1].norm_sq() + x[2].norm_sq()
}

/// One step of the scaled Newton iteration, `X <- (mu X + C / (mu det)) / 2`, `mu = (|C|_F / (|det|
/// |X|_F))^(1/2)` (Higham's Frobenius scaling, `maths/svd3.md` PO.1), and the squared length of the
/// step. `det X` must not be zero.
#[inline(always)]
fn step<S: Real>(x: [Vec3<S>; 3]) -> ([Vec3<S>; 3], S) {
    let half = S::lit(0.5);
    let (c, det) = cofactors(x);
    // `(|C| / |X|)^(1/2)` before the division by `|det|`: the squares of `C` and `det` would not
    // both stay finite inside the `f32` window.
    let mu = ((norm_sq(&c) / norm_sq(&x)).sqrt() / det.abs()).sqrt();
    let k = S::one() / (mu * det);
    let next = [0, 1, 2].map(|j| (x[j].scale(mu) + c[j].scale(k)).scale(half));
    let d = [0, 1, 2].map(|j| next[j] - x[j]);
    (next, norm_sq(&d))
}

/// Newton's iteration from `x`, `det x > 0`, until a step is at most `4u sqrt(3)` or `STEPS` have
/// run; the step that detects convergence is applied, the one the tangents need (`0102` decision
/// 3).
#[inline]
fn newton<S: Real>(x: [Vec3<S>; 3]) -> [Vec3<S>; 3] {
    let tol = S::lit(match S::PRECISION {
        Precision::F64 => const { 48.0 * pow2(-106) },
        Precision::F32 => const { 48.0 * pow2(-48) },
    });
    let (mut x, mut d) = step(x);
    let mut k = 1;
    while k < STEPS {
        let calm = d.le(tol);
        (x, d) = S::branch(calm, || (x, d), || step(x));
        k = calm.decide(|| STEPS, || k + 1, |_, _, next| next);
    }
    x
}

/// `svd3`'s `U V^T`, then one Newton step, which squares the frame's error of the product of two
/// frames (7.8-8.6 `u` without it) and projects the tangent onto the one it already has.
#[inline]
fn through_svd3<S: Real>(a: &Mat3<S>) -> [Vec3<S>; 3] {
    let (u, _, v) = svd3(a);
    let vt = v.transpose();
    step([0, 1, 2].map(|j| u * vt.col(j))).0
}

/// The reflection arm, `det x <= -tau |x|_F^3` (`0102` decision 4): `Q`, the orthogonal polar
/// factor of `x` (`det Q = -1`), by Newton; `H = sym(Q^T x)`, whose eigenvalues are the singular
/// values; its smallest eigenvalue `l` from `eig3`, refined by one Newton step on `det(H - l I)`;
/// the projector `P = adj(H - l I) / tr adj(H - l I)` onto its eigenvector; `R = Q (I - 2P)`
/// (`maths/svd3.md` PO.2). Every step after `eig3` is rational in `x` and `l`, and the refined `l`
/// has the root's derivative whatever `eig3`'s lanes are, so the tangent is exact wherever the
/// smallest singular value is simple, the two largest tied or not. Where `tr adj(H - l I)` is
/// below `16 u |x|_F^2` the smallest singular value is tied to rounding and the nearest rotation
/// is not unique: `svd3`'s `U V^T`, one of the nearest. A last Newton step polishes the frame.
#[inline]
fn reflect<S: Real>(x: [Vec3<S>; 3], a: &Mat3<S>, n2: S) -> [Vec3<S>; 3] {
    let (half, two) = (S::lit(0.5), S::lit(2.0));
    let q = Mat3::from_cols(newton(x));
    let h0 = q.transpose() * *a;
    let h = (h0 + h0.transpose()).scale(half);
    let (l, _) = eig3(&h);
    let shifted = |l: S| {
        let mut y = h;
        for k in 0..3 {
            y.set(k, k, h.get(k, k) - l);
        }
        let (c, det) = cofactors([y.col(0), y.col(1), y.col(2)]);
        (c, det, c[0].0[0] + c[1].0[1] + c[2].0[2])
    };
    let (_, det, tr) = shifted(l.0[0]);
    let floor = S::lit(match S::PRECISION {
        Precision::F64 => const { 16.0 * pow2(-53) },
        Precision::F32 => const { 16.0 * pow2(-24) },
    }) * n2;
    S::branch(
        tr.le(floor),
        || through_svd3(a),
        || {
            // `d/dl det(H - l I) = -tr adj(H - l I)`.
            let (c, _, tr) = shifted(l.0[0] + det / tr);
            let k = two / tr;
            // `adj(H - l I)` is symmetric, so its cofactor columns are its columns.
            // One more step of the iteration squares the frame's error left by `P`'s rounding
            // (16-20 `u` without it); it projects the tangent onto the one it already has.
            step([0, 1, 2].map(|j| q.col(j) - (q * c[j]).scale(k))).0
        },
    )
}

/// The rotation nearest to `m` in the Frobenius norm.
///
/// Returns `R` in SO(3), orthonormal to a few `u` with `det R > 0`, for every finite input: the
/// nearest rotation `U diag(1, 1, det U V^T) V^T` where it is unique, one of the nearest where it is
/// not (`maths/so3.md` SO.13(c)). With `svd3`'s signed convention this is `U V^T`; this routine is
/// that projection, faster and with the right tangents near SO(3) (`docs/NUMERICS.md` §19, `0102`,
/// draft):
///
/// - where `det m >= tau |m|_F^3` (`tau = 2^-26` at `f64`, `2^-12` at `f32`): Newton's polar
///   iteration with Frobenius scaling, `X <- (mu X + X^-T / mu) / 2`, `X^-T` from the cofactors,
///   two or three steps near SO(3);
/// - where `det m <= -tau |m|_F^3`: the same iteration gives the orthogonal polar factor `Q`, `det
///   Q = -1`, and the nearest rotation is `Q (I - 2 P)`, `P` the projector onto the eigenvector of
///   the smallest eigenvalue of `sym(Q^T m)` (`maths/svd3.md` PO.2);
/// - between: `svd3`'s `U V^T`.
///
/// The matrix is first scaled by an exact power of two, as `svd3`'s.
///
/// # Accuracy
///
/// Within a few `u |m|_F / (sigma_2 + d sigma_3)` of the nearest rotation, the conditioning of
/// the projection (`maths/svd3.md` PO.3), `sigma` its singular values and `d` the sign of `det m`.
///
/// # Domain
///
/// Every input is legal and nothing is asserted (`docs/API.md` R6). The zero matrix gives the
/// identity. A non-finite entry gives a matrix of `NaN`, never a plausible rotation.
///
/// A `Dual` result differentiates the iteration: in the Newton and reflection arms it is the
/// derivative of the nearest rotation wherever that is unique, ties among the singular values
/// included (`maths/svd3.md` PO.3), to the convergence of the exit. Between them it is `svd3`'s:
/// exact where the singular values are simple, finite and in general no directional derivative at
/// a tie.
///
/// # Example
///
/// ```
/// use helicoid_linalg::{nearest_rotation, Mat3, Vector};
///
/// // A rotation by 90 degrees about z, scaled and perturbed.
/// let m = Mat3::from_rows([
///     Vector([0.01_f64, -2.0, 0.0]),
///     Vector([2.0, 0.0, 0.02]),
///     Vector([0.0, -0.01, 2.0]),
/// ]);
/// let r = nearest_rotation(&m);
/// assert!((r.get(1, 0) - 1.0).abs() < 1e-3 && (r.get(0, 1) + 1.0).abs() < 1e-3);
/// ```
#[inline]
pub fn nearest_rotation<S: Real>(m: &Mat3<S>) -> Mat3<S> {
    let (zero, one) = (S::zero(), S::one());
    let e: [S; 9] = core::array::from_fn(|k| m.get(k % 3, k / 3));
    let finite = e[1..]
        .iter()
        .fold(is_finite(e[0]), |ok, &x| ok.and(is_finite(x)));
    let e = e.map(|x| S::select(finite, x, zero));
    let (lo, hi, steps) = match S::PRECISION {
        Precision::F64 => WINDOW_F64,
        Precision::F32 => WINDOW_F32,
    };
    let mx = largest(&e);
    let inside = S::lit(lo).le(mx).and(mx.lt(S::lit(hi))).or(mx.le(zero));
    let (e, _, _) = S::branch(inside, || (e, mx, one), || ladder(&steps, (e, mx, one)));
    let x = [
        Vector([e[0], e[1], e[2]]),
        Vector([e[3], e[4], e[5]]),
        Vector([e[6], e[7], e[8]]),
    ];
    let a = Mat3::from_cols(x);

    let (_, det) = cofactors(x);
    let n2 = norm_sq(&x);
    let tau = S::lit(match S::PRECISION {
        Precision::F64 => const { pow2(-26) },
        Precision::F32 => const { pow2(-12) },
    });
    let th = tau * n2 * n2.sqrt();
    let r = S::branch(
        th.lt(det),
        || newton(x),
        || S::branch(det.lt(-th), || reflect(x, &a, n2), || through_svd3(&a)),
    );
    let nan = S::zero() / S::zero();
    Mat3::from_cols(r.map(|c| Vector(c.0.map(|y| S::select(finite, y, nan)))))
}
