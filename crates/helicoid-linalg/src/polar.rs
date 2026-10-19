//! The rotation nearest to a 3x3 matrix (`docs/NUMERICS.md` §19, decision `0102`).
//!
//! Newton's polar iteration where the nearest rotation is the orthogonal polar factor and the
//! determinant resolves it, a reflection through `svd3`'s last right singular vector where the
//! determinant is negative, and `svd3`'s `U V^T` between; generic over `S: Real`.

use crate::eig3::{eig3_values, into_window};
use crate::matrix::Mat3;
use crate::real::{pow2, Mask, Precision, Real};
use crate::svd3::{svd3, WINDOW_F32, WINDOW_F64};
use crate::vector::{Vec3, Vector};

/// At most eight Newton steps, one more than `0102`'s study needed past `tau`.
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

/// One step of the scaled Newton iteration from `x`, its cofactors `c`, `det` and `n2 = |x|_F^2`:
/// `X <- (mu X + C / (mu det)) / 2`, `mu = (|C|_F / (|det| |X|_F))^(1/2)` (Higham's Frobenius
/// scaling, `maths/svd3.md` PO.1), and the squared length of the step. `det` must not be zero.
#[inline(always)]
fn step_from<S: Real>(x: [Vec3<S>; 3], c: [Vec3<S>; 3], det: S, n2: S) -> ([Vec3<S>; 3], S) {
    let half = S::lit(0.5);
    // `(|C| / |X|)^(1/2)` before the division by `|det|`: the squares of `C` and `det` would not
    // both stay finite inside the `f32` window.
    let mu = ((norm_sq(&c) / n2).sqrt() / det.abs()).sqrt();
    // The `1/2` in both scalars, exactly: the bits of halving the sum.
    let (a, k) = (half * mu, half / (mu * det));
    let next = [0, 1, 2].map(|j| x[j].scale(a) + c[j].scale(k));
    let d = [0, 1, 2].map(|j| next[j] - x[j]);
    (next, norm_sq(&d))
}

#[inline(always)]
fn step<S: Real>(x: [Vec3<S>; 3]) -> ([Vec3<S>; 3], S) {
    let (c, det) = cofactors(x);
    step_from(x, c, det, norm_sq(&x))
}

/// Newton's iteration after its first step `(x, d)`, until a step is at most `4u sqrt(3)` or
/// `STEPS` have run; the step that detects convergence is applied, the one the tangents need
/// (`0102` decision 3).
#[inline(always)]
fn newton<S: Real>((mut x, mut d): ([Vec3<S>; 3], S)) -> [Vec3<S>; 3] {
    let tol = S::lit(match S::PRECISION {
        Precision::F64 => const { 48.0 * pow2(-106) },
        Precision::F32 => const { 48.0 * pow2(-48) },
    });
    let mut k = 1;
    while k < STEPS {
        let calm = d.le(tol);
        (x, d) = S::branch(calm, || (x, d), || step(x));
        k = calm.decide(|| STEPS, || k + 1, |_, _, next| next);
    }
    x
}

/// The last Newton steps of an arm whose `x` is orthogonal to a small `d`: a step leaves an error of
/// about the square of its length, so a second runs only where the first was longer than `sqrt(u)`
/// (`0102` decision 5). It projects the tangent onto the one it already has.
#[inline]
fn polish<S: Real>(x: [Vec3<S>; 3]) -> [Vec3<S>; 3] {
    let u = S::lit(match S::PRECISION {
        Precision::F64 => const { pow2(-53) },
        Precision::F32 => const { pow2(-24) },
    });
    let (x, d) = step(x);
    S::branch(d.le(u), || x, || step(x).0)
}

/// `svd3`'s `U V^T`, then `polish`, which squares the error of the product of two frames.
#[cold]
#[inline(never)]
fn through_svd3<S: Real>(x: [Vec3<S>; 3]) -> [Vec3<S>; 3] {
    let (u, _, v) = svd3(&Mat3::from_cols(x));
    let vt = v.transpose();
    polish([0, 1, 2].map(|j| u * vt.col(j)))
}

/// The reflection arm, `det x < -tau |x|_F^3` (`0102` decision 4): `Q`, the orthogonal polar
/// factor of `x` (`det Q = -1`), by Newton; `H = sym(Q^T x)`, whose eigenvalues are the singular
/// values; its smallest eigenvalue `l` from `eig3`, refined by one Newton step on `det(H - l I)`;
/// the projector `P = adj(H - l I) / tr adj(H - l I)` onto its eigenvector; `R = Q (I - 2P)`
/// (`maths/svd3.md` PO.2), then `polish`. Every step after `eig3` is rational in `x` and `l`, and
/// the refined `l` has the root's derivative whatever `eig3`'s lanes are, so the tangent is exact
/// wherever the smallest singular value is simple, the two largest tied or not.
///
/// The margin returned is at most zero where the gap from `l` to the next eigenvalue is at
/// rounding, `e = 16 u |x|_F`: `tr adj(H - l I)`, the product of the gaps to the other two, at most
/// `e (g + e)`, `g` the spread of the three. The nearest rotation is not unique there, and the
/// caller takes `svd3`'s, one of the nearest; the frame returned is `x`.
#[cold]
#[inline(never)]
fn reflect<S: Real>(x: [Vec3<S>; 3], n2: S) -> ([Vec3<S>; 3], S) {
    let two = S::lit(2.0);
    let q = Mat3::from_cols(newton(step(x)));
    let h0 = q.transpose() * Mat3::from_cols(x);
    let h = (h0 + h0.transpose()).scale(S::lit(0.5));
    let l = eig3_values(&h);
    let shifted = |l: S| {
        let mut y = h;
        for k in 0..3 {
            y.set(k, k, h.get(k, k) - l);
        }
        let (c, det) = cofactors([y.col(0), y.col(1), y.col(2)]);
        (c, det, c[0].0[0] + c[1].0[1] + c[2].0[2])
    };
    let (_, det, tr) = shifted(l.0[0]);
    let e = S::lit(match S::PRECISION {
        Precision::F64 => const { 16.0 * pow2(-53) },
        Precision::F32 => const { 16.0 * pow2(-24) },
    }) * n2.sqrt();
    let margin = tr - e * ((l.0[2] - l.0[0]) + e);
    let r = S::branch(
        margin.le(S::zero()),
        || x,
        || {
            // `d/dl det(H - l I) = -tr adj(H - l I)`.
            let (c, _, tr) = shifted(l.0[0] + det / tr);
            let k = two / tr;
            // `adj(H - l I)` is symmetric, so its cofactor columns are its columns.
            polish([0, 1, 2].map(|j| q.col(j) - (q * c[j]).scale(k)))
        },
    );
    (r, margin)
}

/// The rotation nearest to `m` in the Frobenius norm.
///
/// Returns `R` in SO(3), orthonormal to a few `u` with `det R > 0`, for every finite input: the
/// nearest rotation `U diag(1, 1, det U V^T) V^T` where it is unique, one of the nearest where it is
/// not (`maths/so3.md` SO.13(c)). With `svd3`'s signed convention this is `U V^T`; this routine is
/// that projection, faster and with the right tangents near SO(3) (`docs/NUMERICS.md` §19,
/// `0102`):
///
/// - where `det m > tau |m|_F^3` (`tau = 2^-26` at `f64`, `2^-12` at `f32`): Newton's polar
///   iteration with Frobenius scaling, `X <- (mu X + X^-T / mu) / 2`, `X^-T` from the cofactors,
///   two or three steps near SO(3);
/// - where `det m < -tau |m|_F^3`: the same iteration gives the orthogonal polar factor `Q`, `det
///   Q = -1`, and the nearest rotation is `Q (I - 2 P)`, `P` the projector onto the eigenvector of
///   the smallest eigenvalue of `sym(Q^T m)` (`maths/svd3.md` PO.2);
/// - between, and where the smallest singular value of a reflection is tied to rounding: `svd3`'s
///   `U V^T`.
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
    // A non-finite entry iterates on the zero matrix and reports a matrix of NaN; the rotation
    // needs no scaling back.
    let window = match S::PRECISION {
        Precision::F64 => &WINDOW_F64,
        Precision::F32 => &WINDOW_F32,
    };
    let e: [S; 9] = core::array::from_fn(|k| m.get(k % 3, k / 3));
    let (e, finite, _, _) = into_window(e, window);
    let x = [
        Vector([e[0], e[1], e[2]]),
        Vector([e[3], e[4], e[5]]),
        Vector([e[6], e[7], e[8]]),
    ];
    let (c, det) = cofactors(x);
    let n2 = norm_sq(&x);
    let tau = S::lit(match S::PRECISION {
        Precision::F64 => const { pow2(-26) },
        Precision::F32 => const { pow2(-12) },
    });
    let th = tau * n2 * n2.sqrt();
    let r = S::branch(
        th.lt(det),
        || newton(step_from(x, c, det, n2)),
        || {
            // Between the arms the margin is zero: `svd3`'s, shared with a tied reflection.
            let (r, margin) = S::branch(det.lt(-th), || reflect(x, n2), || (x, S::zero()));
            S::branch(margin.le(S::zero()), || through_svd3(x), || r)
        },
    );
    let nan = S::zero() / S::zero();
    Mat3::from_cols(r.map(|c| Vector(c.0.map(|y| S::select(finite, y, nan)))))
}
