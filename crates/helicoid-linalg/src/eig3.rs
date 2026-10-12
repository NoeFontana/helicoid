//! Eigendecomposition of a symmetric 3x3 matrix (`docs/PHASE2.md` §6, `docs/NUMERICS.md` §17,
//! decision `0081`).
//!
//! Cyclic Jacobi with a fixed sweep count, generic over `S: Real`: every decision is a mask, every
//! square root and division has a safe argument, and the only operations are `+ - * /`, `sqrt`,
//! `abs` and `copysign`.

use crate::cubic::pow2;
use crate::matrix::Mat3;
use crate::real::{is_finite, Blend, Mask, Precision, Real};
use crate::vector::{Vec3, Vector};

/// Sweeps over the pairs `(0, 1), (0, 2), (1, 2)`. Four leave every off-diagonal entry below
/// `3e-4 u |A|` at `f64` and `3e-13 u |A|` at `f32`, three leave up to `3e-5 |A|` (`0081`).
const SWEEPS: usize = 4;

/// The window `[lo, hi)` of the largest `|entry|` inside which the matrix is not scaled, and the
/// two-sided power-of-two steps that bring any other finite matrix into it (`0081` decision 3).
/// Inside it `d^2 + 4 a_pq^2` cannot overflow and every square that matters is normal.
const WINDOW_F64: (f64, f64, [f64; 4]) = (
    pow2(-128),
    pow2(128),
    [pow2(512), pow2(256), pow2(128), pow2(64)],
);
const WINDOW_F32: (f64, f64, [f64; 4]) =
    (pow2(-32), pow2(32), [pow2(64), pow2(32), pow2(16), pow2(8)]);

/// The lower triangle `[xx, yy, zz, xy, xz, yz]` and the power of two that undoes its scaling.
type Scaled<S> = ([S; 6], S);

/// [`eig3`]'s scaling steps, out of line: rarely taken, as `cubic`'s ladders. At `2^k` the entries
/// are divided by it where the largest is at least `2^k` and multiplied where it is below
/// `2^(1-k)`; every product is exact.
#[cold]
#[inline(never)]
fn ladder<S: Real>(steps: &[f64], start: Scaled<S>) -> Scaled<S> {
    let (one, two) = (S::one(), S::lit(2.0));
    steps.iter().fold(start, |(e, back), &k| {
        let (big, small) = (S::lit(k), S::lit(1.0 / k));
        let m = largest(&e);
        let (down, up) = (big.le(m), (m * big).lt(two));
        let f = S::select(down, small, S::select(up, big, one));
        let b = S::select(down, big, S::select(up, small, one));
        (e.map(|x| x * f), back * b)
    })
}

/// The largest `|entry|`.
#[inline]
fn largest<S: Real>(e: &[S; 6]) -> S {
    e.iter().fold(S::zero(), |m, &x| {
        let x = x.abs();
        S::select(m.lt(x), x, m)
    })
}

/// One rotation in the `(p, q)` plane, `r` the third index: `a_pq` is annihilated, `(a_pp, a_qq,
/// a_rp, a_rq)` returned updated and the columns `v_p`, `v_q` rotated (`NUMERICS.md` §17).
///
/// With `d = a_qq - a_pp`, `w = 2 a_pq sign(d)`, `den = |d| + sqrt(d^2 + w^2)` and `g = sqrt(den^2 +
/// w^2)`: `t = w / den`, `c = den / g`, `s = w / g`. Below the smallest normal `d^2 + w^2` is `1`
/// under the root (the safe argument), which leaves an orthogonal rotation by about `w` of a block
/// below `2^-383 |A|` (`f64`), `2^-31 |A|` (`f32`) inside the window: `c` and `s` share `g`, so
/// `c^2 + s^2 = 1` to rounding only while `g` is normal.
#[inline(always)]
#[allow(clippy::type_complexity)]
fn rotate<S: Real>(
    (app, aqq, apq, arp, arq): (S, S, S, S, S),
    (vp, vq): (Vec3<S>, Vec3<S>),
) -> ((S, S, S, S), (Vec3<S>, Vec3<S>)) {
    let (one, two) = (S::one(), S::lit(2.0));
    let tiny = S::lit(match S::PRECISION {
        Precision::F64 => pow2(-1022),
        Precision::F32 => pow2(-126),
    });
    let d = aqq - app;
    let w = two * apq * one.copysign(d);
    let h2 = d * d + w * w;
    let den = d.abs() + S::select(h2.lt(tiny), one, h2).sqrt();
    let t = w / den;
    let g = (den * den + w * w).sqrt();
    let (c, s) = (den / g, w / g);
    let ta = t * apq;
    (
        (app - ta, aqq + ta, c * arp - s * arq, s * arp + c * arq),
        (vp.scale(c) - vq.scale(s), vp.scale(s) + vq.scale(c)),
    )
}

/// The eigenvalues, ascending, and an orthonormal basis of eigenvectors of the symmetric matrix `a`.
///
/// Returns `(values, vectors)`: `values.0[0] <= values.0[1] <= values.0[2]`, column `k` of `vectors`
/// an eigenvector of `values.0[k]`, `vectors` orthonormal with `det = +1`. The sign of a column is
/// what the rotations produce; only the parity of the three is fixed.
///
/// Cyclic Jacobi (`docs/NUMERICS.md` §17, `0081`): four sweeps of the rotations `(0, 1)`, `(0, 2)`,
/// `(1, 2)`, each annihilating its entry; the eigenvalues are the diagonal that is left and the
/// vectors the product of the rotations, sorted ascending by a compare-exchange network whose
/// exchanges flip the first column's sign to keep `det = +1`, each column normalised once at the
/// end (the rotations' `c^2 + s^2 - 1` is a length error: `19 u` of frame error without it, `7.7 u`
/// with it). The count is fixed: no convergence
/// test, no data-dependent loop. The matrix is first scaled by an exact power of two into a window
/// (`[2^-128, 2^128)` at `f64`, `[2^-32, 2^32)` at `f32`), and the eigenvalues scaled back.
///
/// # Accuracy
///
/// Backward stable: the eigenvalues are within a few `u |A|` of the exact ones and column `k` within
/// a few `u |A| / gap_k` radians, `gap_k` its distance to the nearest other eigenvalue, at every
/// gap, ties and a multiple of the identity included. The corpus id `eig3` reads at most 4.76 `u`
/// (`f64`) and 3.24 `u` (`f32`) over every stratum (`0081`); the bounds `eig3_tests` asserts are its
/// header's.
///
/// # Domain
///
/// Every input is legal and nothing is asserted (`docs/API.md` R6, `0081`): only the lower triangle
/// is read, and the caller owns symmetry. For every finite input the eigenvalues are finite unless
/// one exceeds the largest finite number, and the frame is orthonormal to a few `u`. A matrix whose
/// largest entry is subnormal keeps only the digits it has. A non-finite entry gives `NaN`
/// eigenvalues, never a plausible number, and the identity frame.
///
/// A `Dual` result differentiates the rotations taken. Where the eigenvalues are simple that is
/// first-order perturbation theory, to the convergence of the last sweep; at a repeated eigenvalue
/// it is one of the one-sided derivatives (a multiple of the identity: each eigenvalue has the
/// derivative of one diagonal entry, and the vectors none).
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
    let e = [
        a.get(0, 0),
        a.get(1, 1),
        a.get(2, 2),
        a.get(1, 0),
        a.get(2, 0),
        a.get(2, 1),
    ];
    // A non-finite entry iterates on the zero matrix, the safe argument, and reports through the
    // eigenvalues.
    let finite = e
        .iter()
        .fold(is_finite(e[0]), |ok, &x| ok.and(is_finite(x)));
    let e = e.map(|x| S::select(finite, x, zero));
    let (lo, hi, steps) = match S::PRECISION {
        Precision::F64 => WINDOW_F64,
        Precision::F32 => WINDOW_F32,
    };
    let m = largest(&e);
    let inside = S::lit(lo).le(m).and(m.lt(S::lit(hi))).or(m.le(zero));
    let ([mut d0, mut d1, mut d2, mut o01, mut o02, mut o12], back) =
        S::branch(inside, || (e, one), || ladder(&steps, (e, one)));

    let (mut v0, mut v1, mut v2) = (
        Vector([one, zero, zero]),
        Vector([zero, one, zero]),
        Vector([zero, zero, one]),
    );
    for _ in 0..SWEEPS {
        let ((p, q, rp, rq), (vp, vq)) = rotate((d0, d1, o01, o02, o12), (v0, v1));
        (d0, d1, o02, o12, o01, v0, v1) = (p, q, rp, rq, zero, vp, vq);
        let ((p, q, rp, rq), (vp, vq)) = rotate((d0, d2, o02, o01, o12), (v0, v2));
        (d0, d2, o01, o12, o02, v0, v2) = (p, q, rp, rq, zero, vp, vq);
        let ((p, q, rp, rq), (vp, vq)) = rotate((d1, d2, o12, o01, o02), (v1, v2));
        (d1, d2, o01, o02, o12, v1, v2) = (p, q, rp, rq, zero, vp, vq);
    }

    let mut sign = one;
    let mut exchange = |(x, vx): (S, Vec3<S>), (y, vy): (S, Vec3<S>)| {
        let swap = y.lt(x);
        sign = S::select(swap, -sign, sign);
        (
            (S::select(swap, y, x), Vec3::blend(swap, vy, vx)),
            (S::select(swap, x, y), Vec3::blend(swap, vx, vy)),
        )
    };
    let (p0, p1) = exchange((d0, v0), (d1, v1));
    let (p1, p2) = exchange(p1, (d2, v2));
    let (p0, p1) = exchange(p0, p1);
    let nan = S::zero() / S::zero();
    let value = |l: S| S::select(finite, l * back, nan);
    let unit = |v: Vec3<S>| {
        let n = v.norm();
        Vector(v.0.map(|x| x / n))
    };
    (
        Vector([value(p0.0), value(p1.0), value(p2.0)]),
        Mat3::from_cols([unit(p0.1.scale(sign)), unit(p1.1), unit(p2.1)]),
    )
}
