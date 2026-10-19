//! The singular value decomposition of a 3x3 matrix (`docs/PHASE2.md` §6, `docs/NUMERICS.md` §18,
//! decision `0101`).
//!
//! One-sided Jacobi on `eig3`'s rotation, with `eig3`'s exit and polish, generic over `S: Real`:
//! every decision is a mask, every square root and division has a safe argument, and the only
//! operations are `+ - * /`, `sqrt`, `abs` and `copysign`.

use crate::eig3::{angle, into_window, largest, smallest_normal, unit, Window};
use crate::matrix::Mat3;
use crate::real::{pow2, Blend, Mask, Precision, Real};
use crate::vector::{Vec3, Vector};

/// At most twelve rotations, `eig3`'s bound: a scratch study of the exit below reached at most
/// twelve on every family but rank-deficient `f32` input, whose extra rotations move nothing scored
/// (`0101`).
const ROTATIONS: usize = 12;

/// The window `[lo, hi)` of the largest `|entry|` inside which the matrix is not scaled, and its
/// steps. The rotation's radicand is of degree four in `A`: inside the window it cannot overflow,
/// and a pair whose radicand is below the smallest normal, which `angle`'s safe arm leaves
/// unrotated, costs the finish at most `2^-256` (`f64`), `2^-32` (`f32`), that is `2^-128 |A|` and
/// `u |A| / 2` at the floor (`maths/svd3.md` SV.6). The `f32` window reaches `2^28` so that a
/// covariance in mm² takes no ladder; its steps are `eig3`'s with a second `2^64`, which carries
/// the smallest subnormal inside.
const WINDOW_F64: Window = Window {
    lo: pow2(-128),
    hi: pow2(128),
    steps: &[pow2(512), pow2(256), pow2(128), pow2(64)],
};
const WINDOW_F32: Window = Window {
    lo: pow2(-7),
    hi: pow2(28),
    steps: &[pow2(64), pow2(64), pow2(32), pow2(16), pow2(8)],
};
const _: () = assert!(WINDOW_F64.reaches(Precision::F64) && WINDOW_F32.reaches(Precision::F32));

/// The columns of `B = A V` and of `V` in the frame of the next rotation's roles `(p, q, r)`, and
/// their Gram entries `[a_p, a_q, a_r, g_pq, g_qr, g_rp]`, `a = |b|^2`, `g = b . b'`.
type State<S> = ([Vec3<S>; 3], [Vec3<S>; 3], [S; 6]);

/// The Gram entries of three columns, in [`State`]'s order.
#[inline(always)]
fn gram<S: Real>([b0, b1, b2]: [Vec3<S>; 3]) -> [S; 6] {
    [
        b0.norm_sq(),
        b1.norm_sq(),
        b2.norm_sq(),
        b0.dot(b1),
        b1.dot(b2),
        b2.dot(b0),
    ]
}

/// `(x, y) -> (c x - s y, s x + c y)`.
#[inline(always)]
fn rot<S: Real>(c: S, s: S, x: Vec3<S>, y: Vec3<S>) -> (Vec3<S>, Vec3<S>) {
    (x.scale(c) - y.scale(s), x.scale(s) + y.scale(c))
}

/// The rotation of a pair of columns that annihilates their inner product `g` (`NUMERICS.md` §18,
/// `maths/svd3.md` SV.1): `eig3`'s on the Gram block `[[a_p, g], [g, a_q]]`. `(c, s)`.
#[inline(always)]
fn pair<S: Real>(ap: S, aq: S, g: S) -> (S, S) {
    let w = S::lit(2.0) * g;
    let (_, c, s, _, _) = angle(aq - ap, g, w * w);
    (c, s)
}

/// One rotation of the pair `(p, q)`; the next roles are `(q, r, p)`. `a_r` is carried: `b_r` did
/// not move.
#[inline(always)]
fn turn<S: Real>(([bp, bq, br], [vp, vq, vr], g): State<S>) -> State<S> {
    let (c, s) = pair(g[0], g[1], g[3]);
    let (bp, bq) = rot(c, s, bp, bq);
    let (vp, vq) = rot(c, s, vp, vq);
    let g = [
        bq.norm_sq(),
        g[2],
        bp.norm_sq(),
        bq.dot(br),
        br.dot(bp),
        bp.dot(bq),
    ];
    ([bq, br, bp], [vq, vr, vp], g)
}

/// The last step, from a state within the exit: the rotations of `(p, q)` and `(q, r)`, both
/// computed from that state, applied in turn; the pair `(r, p)` is the one the last rotation
/// annihilated. Their values move by at most the entries; their `Dual` lanes annihilate the
/// entries' lanes to first order, which no value test can see (`0084` decision 4).
///
/// The columns leave in the order `(r, p, q)`: index order for an exit after the first rotation,
/// whose roles are `(1, 2, 0)`, so a diagonal matrix keeps its axes before the sort and the zero
/// matrix keeps `V = I`. A later exit leaves another cyclic order, which the sort decides.
#[inline(always)]
fn polish<S: Real>(([b0, b1, b2], [v0, v1, v2], g): State<S>) -> ([Vec3<S>; 3], [Vec3<S>; 3]) {
    let (ca, sa) = pair(g[0], g[1], g[3]);
    let (cb, sb) = pair(g[1], g[2], g[4]);
    let ((b0, b1), (v0, v1)) = (rot(ca, sa, b0, b1), rot(ca, sa, v0, v1));
    let ((b1, b2), (v1, v2)) = (rot(cb, sb, b1, b2), rot(cb, sb, v1, v2));
    ([b2, b0, b1], [v2, v0, v1])
}

/// A unit vector orthogonal to the unit `u`: Gram–Schmidt of the axis of its smallest `|component|`
/// (the first of equals), whose residual has squared norm at least `2/3`. `e_2` for `u = e_1`, so
/// the zero matrix keeps the identity frame.
#[inline]
fn orthogonal<S: Real>(u: Vec3<S>) -> Vec3<S> {
    let (zero, one) = (S::zero(), S::one());
    let [x, y, z] = u.0.map(|c| c.abs());
    // The smallest of `(x, y, z)`, the first of equals: `x` unless another is strictly smaller.
    let not_x = y.lt(x).or(z.lt(x));
    let z_not_y = z.lt(y);
    let axis = |i: usize| -> Vec3<S> {
        let mut e = [zero; 3];
        e[i] = one;
        Vector(e)
    };
    let e = Vec3::blend(not_x, Vec3::blend(z_not_y, axis(2), axis(1)), axis(0));
    let w = e - u.scale(u.dot(e));
    w.scale(one / w.norm())
}

/// `(U, sigma, V)` with `a = U diag(sigma) V^T`: `U` and `V` rotations, `sigma_1 >= sigma_2 >=
/// |sigma_3|` and `sigma_3` of the sign of `det a`.
///
/// One-sided Jacobi with a deterministic exit (`docs/NUMERICS.md` §18, `0101`): the columns of
/// `B = a V` are rotated in pairs `(0, 1), (1, 2), (2, 0), ...` by `eig3`'s rotation of their Gram
/// block, each annihilating its pair's inner product, until every pair is orthogonal to `4u`
/// relative or holds a null column, at most `u |B|`, or twelve have run; then one step that rotates
/// the two remaining pairs, which moves the values by at most those inner products and gives a
/// `Dual` its first-order lanes. The columns are sorted by length with exchanges that are quarter
/// turns, so `V` stays a rotation; `U` is Gram–Schmidt on them, the projection taken twice, with
/// `u_3 = u_1 x u_2`, so it is a rotation by construction, converged or not, and
/// `sigma_3 = u_3 . b_3` carries the sign of `det a`; one masked pass orders rounding-level ties.
/// The exit compares values that every target computes to the bit (D16). The matrix is first scaled
/// by an exact power of two into a window (`[2^-128, 2^128)` at `f64`, `[2^-7, 2^28)` at `f32`),
/// and `sigma` scaled back.
///
/// With this ordering `U V^T` is the rotation nearest to `a` in the Frobenius norm wherever that is
/// unique (`maths/svd3.md` SV.4).
///
/// # Accuracy
///
/// Backward stable column by column: `a = U diag(sigma) V^T + E` with `|e_j|` a few `u |a_j|`, so
/// the singular values are within a few `u |a|` of the exact ones and each singular pair within a
/// few `u |a| / gap_i` radians, `gap_i` its distance to the nearest other `|sigma_j|` (or to
/// `-|sigma_j|`). The relative accuracy one-sided Jacobi has on a matrix whose columns are graded
/// is measured by the corpus and not promised (`0101` decision 1).
///
/// # Domain
///
/// Every input is legal and nothing is asserted (`docs/API.md` R6): for every finite input the
/// singular values are finite unless one exceeds the largest finite number, and `U`, `V` are
/// rotations to a few `u`. Every finite matrix, a subnormal one included, is scaled into the window
/// exactly, so it loses only the rounding of `sigma` scaled back. The zero matrix gives `sigma = 0`
/// and `U = V = I`. A non-finite entry gives `NaN` singular values, never a plausible number, and
/// `U = V = I`.
///
/// A `Dual` result differentiates the rotations taken, the last step's included: where `sigma_i`
/// is simple, `d sigma_i = u_i^T (dA) v_i`, and where no other `|sigma_j|` equals `|sigma_i|`, the
/// singular pair's lanes are first-order perturbation theory, to the convergence of the exit. At a
/// repeated singular value the lanes are finite and in general no directional derivative.
///
/// # Example
///
/// ```
/// use helicoid_linalg::{svd3, Mat3, Vector};
///
/// // A reflection scaled by (3, 2, 1): det < 0, so the last singular value is -1.
/// let a = Mat3::from_rows([
///     Vector([0.0_f64, 2.0, 0.0]),
///     Vector([3.0, 0.0, 0.0]),
///     Vector([0.0, 0.0, 1.0]),
/// ]);
/// let (_, sigma, _) = svd3(&a);
/// assert!(sigma.0.iter().zip([3.0, 2.0, -1.0]).all(|(s, w)| (s - w).abs() < 1e-12));
/// ```
#[inline]
pub fn svd3<S: Real>(a: &Mat3<S>) -> (Mat3<S>, Vec3<S>, Mat3<S>) {
    let (zero, one) = (S::zero(), S::one());
    // A non-finite entry iterates on the zero matrix and reports through the singular values.
    let window = match S::PRECISION {
        Precision::F64 => &WINDOW_F64,
        Precision::F32 => &WINDOW_F32,
    };
    let e: [S; 9] = core::array::from_fn(|k| a.get(k % 3, k / 3));
    let (e, finite, _, back) = into_window(e, window);

    let b = [
        Vector([e[0], e[1], e[2]]),
        Vector([e[3], e[4], e[5]]),
        Vector([e[6], e[7], e[8]]),
    ];
    let ident = [
        Vector([one, zero, zero]),
        Vector([zero, one, zero]),
        Vector([zero, zero, one]),
    ];
    let mut st = turn((b, ident, gram(b)));
    // The exit (`0101` decision 3): every pair has `g^2 <= 16 u^2 a_i a_j`, squared so that it
    // takes no root, or holds a null column, `a <= u^2 max a`, whose direction the finish does not
    // need: it costs at most its own norm, `u |B|`.
    let (u2, tol) = match S::PRECISION {
        Precision::F64 => (const { pow2(-106) }, const { pow2(-102) }),
        Precision::F32 => (const { pow2(-48) }, const { pow2(-44) }),
    };
    let (u2, tol) = (S::lit(u2), S::lit(tol));
    let mut k = 1;
    while k < ROTATIONS {
        let g = st.2;
        let f = u2 * largest(&[g[0], g[1], g[2]]);
        let [np, nq, nr] = [g[0], g[1], g[2]].map(|x| x.le(f));
        let ok = |gij: S, ai: S, aj: S| (gij * gij).le(tol * ai * aj);
        let calm = (ok(g[3], g[0], g[1]).or(np).or(nq))
            .and(ok(g[4], g[1], g[2]).or(nq).or(nr))
            .and(ok(g[5], g[2], g[0]).or(nr).or(np));
        st = S::branch(calm, || st, || turn(st));
        // A scalar leaves at its exit; lanes run on, each one held by the branch above.
        k = calm.decide(|| ROTATIONS, || k + 1, |_, _, next| next);
    }
    let (b, v) = polish(st);

    // Sort by length, descending; each exchange is a quarter turn, `(x, y) -> (y, -x)`, on `B`
    // and `V` alike, so `det V` stays +1.
    let mut cols = [0, 1, 2].map(|i| (b[i].norm_sq(), b[i], v[i]));
    for (p, q) in [(0, 1), (1, 2), (0, 1)] {
        let ((ap, bp, vp), (aq, bq, vq)) = (cols[p], cols[q]);
        let swap = ap.lt(aq);
        cols[p] = (
            S::select(swap, aq, ap),
            Vec3::blend(swap, bq, bp),
            Vec3::blend(swap, vq, vp),
        );
        cols[q] = (
            S::select(swap, ap, aq),
            Vec3::blend(swap, -bp, bq),
            Vec3::blend(swap, -vp, vq),
        );
    }
    let [(a0, b0, v0), (_, b1, v1), (_, b2, _)] = cols;

    // `U` by Gram–Schmidt with `u_3 = u_1 x u_2`: a rotation by construction, and the residuals
    // `u_i . b_j` are the converged inner products (`maths/svd3.md` SV.3).
    let null1 = a0.le(zero);
    let s1 = S::select(null1, one, a0).sqrt();
    let u1 = Vec3::blend(null1, ident[0], unit(b0.scale(one / s1)));
    let s1 = S::select(null1, zero, s1);
    // The projection twice: once leaves `u1 . w` at `u |b1|`, not `u |w|`, where `b1` leaves the
    // exit (null, or at the cap) far from orthogonal to `b0`; twice is enough. A residual whose
    // square is below the smallest normal keeps too few digits in its root to be normalised: it is
    // below `2^-511` (`f64`), `2^-63` (`f32`), at most `2^-383 |A|`, `2^-56 |A|` inside the window,
    // and counts as zero.
    let w = b1 - u1.scale(u1.dot(b1));
    let w = w - u1.scale(u1.dot(w));
    let aw = w.norm_sq();
    let null2 = aw.lt(smallest_normal());
    let s2 = S::select(null2, one, aw).sqrt();
    let u2 = S::branch(null2, || orthogonal(u1), || unit(w.scale(one / s2)));
    let s2 = S::select(null2, zero, s2);
    let u3 = unit(u1.cross(u2));
    let s3 = u3.dot(b2);
    // `V` is a product of rotations: re-orthonormalised without a root.
    let v0 = unit(v0);
    let v1 = unit(v1 - v0.scale(v0.dot(v1)));
    let v2 = unit(v0.cross(v1));

    // Rounding-level inversions at ties (`maths/svd3.md` SV.5), each keeping `U Sigma V^T` and
    // both determinants.
    type Pair<S> = (S, S, Vec3<S>, Vec3<S>, Vec3<S>, Vec3<S>);
    let first = |(s1, s2, u1, u2, v1, v2): Pair<S>| -> Pair<S> {
        let swap = s1.lt(s2);
        (
            S::select(swap, s2, s1),
            S::select(swap, s1, s2),
            Vec3::blend(swap, u2, u1),
            Vec3::blend(swap, -u1, u2),
            Vec3::blend(swap, v2, v1),
            Vec3::blend(swap, -v1, v2),
        )
    };
    let (s1, s2, u1, u2, v0, v1) = first((s1, s2, u1, u2, v0, v1));
    let swap = s2.lt(s3.abs());
    let sg = one.copysign(s3);
    let (s2, s3, u2, u3, v1, v2) = (
        S::select(swap, sg * s3, s2),
        S::select(swap, sg * s2, s3),
        Vec3::blend(swap, u3.scale(sg), u2),
        Vec3::blend(swap, u2.scale(-sg), u3),
        Vec3::blend(swap, v2, v1),
        Vec3::blend(swap, -v1, v2),
    );
    let (s1, s2, u1, u2, v0, v1) = first((s1, s2, u1, u2, v0, v1));

    let nan = S::zero() / S::zero();
    let value = |s: S| S::select(finite, s * back, nan);
    (
        Mat3::from_cols([u1, u2, u3]),
        Vector([value(s1), value(s2), value(s3)]),
        Mat3::from_cols([v0, v1, v2]),
    )
}
