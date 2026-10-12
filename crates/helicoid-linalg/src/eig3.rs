//! Eigendecomposition of a symmetric 3x3 matrix (`docs/PHASE2.md` §6, `docs/NUMERICS.md` §17,
//! decisions `0081`, `0084`, draft).
//!
//! Cyclic Jacobi with a deterministic exit, generic over `S: Real`: every decision is a mask, every
//! square root and division has a safe argument, and the only operations are `+ - * /`, `sqrt`,
//! `abs` and `copysign`.

use crate::matrix::Mat3;
use crate::real::{is_finite, pow2, pow2_step, Blend, Mask, Precision, Real};
use crate::vector::{Vec3, Vector};

/// At most twelve rotations, four sweeps: they leave every off-diagonal entry below `3e-4 u |A|` at
/// `f64` and `3e-13 u |A|` at `f32`, and three sweeps up to `3e-5 |A|` (`0081`).
const ROTATIONS: usize = 12;

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

/// The lower triangle `[xx, yy, zz, xy, xz, yz]`, its largest `|entry|` and the power of two that
/// undoes its scaling.
type Scaled<S> = ([S; 6], S, S);

/// [`eig3`]'s scaling steps, out of line: rarely taken, as `cubic`'s ladders. At `2^k` the entries
/// are divided by it where the largest is at least `2^k` and multiplied where it is below
/// `2^(1-k)` (`real::pow2_step`). The largest entry stays normal, so it scales exactly and is carried
/// as `m f`; an entry that falls below the smallest normal rounds, by less than `2^-1074` of `2^k`,
/// far below `u |A|`.
#[cold]
#[inline(never)]
fn ladder<S: Real>(steps: &[f64], start: Scaled<S>) -> Scaled<S> {
    steps.iter().fold(start, |(e, m, back), &k| {
        let (f, b) = pow2_step(k, m);
        (e.map(|x| x * f), m * f, back * b)
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

/// The state before a rotation, in the frame of its roles `(p, q, r)`: the diagonal `[a_pp, a_qq,
/// a_rr]`, `d = a_qq - a_pp`, the pair's entry `a_pq`, `w2 = (2 a_pq)^2` and `y = a_rq`, `a_rp`
/// being zero (the previous rotation annihilated it); then the columns `[v_p, v_q, v_r]`.
type Roles<S> = ([S; 7], [Vec3<S>; 3]);

/// The rotation annihilating `a_pq` (`NUMERICS.md` §17): `(ta, c, s, c^2, t)`.
///
/// `h = sqrt(d^2 + w2)`, `den = |d| + h`, `t = w / den`, `ta = t a_pq = sign(d) (w2 / 2) / den`,
/// `c^2 = den / 2h`, `s = c t`. `ta` is read from `w2`, the square the root saw, not from `a_pq`:
/// the diagonal then moves by exactly the block's eigenvalue shift for that root (`0084`, draft).
/// Below the smallest normal `d^2 + w2` is `1` under the root and `c = 1` (the safe argument,
/// `0081`): that arm rotates by about `w`, a block below `2^-383 |A|` (`f64`), `2^-31 |A|` (`f32`)
/// inside the window.
#[inline(always)]
fn angle<S: Real>(d: S, apq: S, w2: S) -> (S, S, S, S, S) {
    let (one, two, half) = (S::one(), S::lit(2.0), S::lit(0.5));
    let tiny = S::lit(match S::PRECISION {
        Precision::F64 => const { pow2(-1022) },
        Precision::F32 => const { pow2(-126) },
    });
    let sgn = one.copysign(d);
    let h2 = d * d + w2;
    let small = h2.lt(tiny);
    let h = S::select(small, one, h2).sqrt();
    let den = d.abs() + h;
    let t = two * apq * sgn / den;
    let c2 = den / (h + h);
    // `c2` is at least 1/2 in both arms, so its root is safe. The select sits after the root:
    // `select(..).sqrt()` here crashes rustc 1.98.1's LLVM (`0084`).
    let c = S::select(small, one, c2.sqrt());
    (
        sgn * (half * w2) / den,
        c,
        c * t,
        S::select(small, one, c2),
        t,
    )
}

/// One rotation, `x = a_rp` (`FIRST` only; zero after it): the state of the next rotation, whose
/// roles are `(q, r, p)`. The next pair's entry `c z`, `z = t x + y`, reaches the next rotation as
/// `w2 = 4 z^2 c^2`, `c^2` a square root earlier than `c` (`0084`, draft); `d` is recomputed from
/// the stored diagonal.
#[inline(always)]
fn turn<S: Real, const FIRST: bool>(
    ([dp, dq, dr, d, apq, w2, y], [vp, vq, vr]): Roles<S>,
    x: S,
) -> Roles<S> {
    let four = S::lit(4.0);
    let (ta, c, s, c2, t) = angle(d, apq, w2);
    let (z, rp) = if FIRST {
        (t * x + y, c * x - s * y)
    } else {
        (y, -(s * y))
    };
    (
        [
            dq + ta,
            dr,
            dp - ta,
            dr - (dq + ta),
            c * z,
            four * (z * z) * c2,
            rp,
        ],
        [vp.scale(s) + vq.scale(c), vr, vp.scale(c) - vq.scale(s)],
    )
}

/// The last step, from a state whose entries are at most `u |A|`: the rotations of `(p, q)` and
/// `(q, r)`, both computed from that state, in index order (`0084`, draft). Their values move by at
/// most the entries; their `Dual` lanes annihilate the entries' lanes to first order, which no
/// value test can see.
///
/// The index order is that of an exit after the first rotation, whose roles are `(1, 2, 0)`: a
/// multiple of the identity and a diagonal matrix keep their axes in order. A later exit leaves the
/// three in another cyclic order, which the sort makes ascending; it decides only between equal
/// eigenvalues.
#[inline(always)]
fn polish<S: Real>(([dp, dq, dr, d, apq, w2, y], [vp, vq, vr]): Roles<S>) -> [(S, Vec3<S>); 3] {
    let wb = S::lit(2.0) * y;
    let (ta, ca, sa, _, _) = angle(d, apq, w2);
    let (tb, cb, sb, _, _) = angle(dr - dq, y, wb * wb);
    let (vp, vq) = (vp.scale(ca) - vq.scale(sa), vp.scale(sa) + vq.scale(ca));
    let (vq, vr) = (vq.scale(cb) - vr.scale(sb), vq.scale(sb) + vr.scale(cb));
    [(dr + tb, vr), (dp - ta, vp), ((dq + ta) - tb, vq)]
}

/// The eigenvalues, ascending, and an orthonormal basis of eigenvectors of the symmetric matrix `a`.
///
/// Returns `(values, vectors)`: `values.0[0] <= values.0[1] <= values.0[2]`, column `k` of `vectors`
/// an eigenvector of `values.0[k]`, `vectors` orthonormal with `det = +1`. The sign of a column is
/// what the rotations produce; only the parity of the three is fixed.
///
/// Cyclic Jacobi with a deterministic exit (`docs/NUMERICS.md` §17, `0081`, `0084`, draft):
/// rotations of `(0, 1)`, `(1, 2)`, `(2, 0)`, ... each annihilating its entry, until both entries
/// left are at most `u |A|` or twelve have run; then one step that rotates the two remaining pairs,
/// which moves the values by at most those entries and gives a `Dual` its first-order lanes. The
/// eigenvalues are the diagonal that is left and the vectors the product of the rotations, sorted
/// ascending by a compare-exchange network whose exchanges flip the first column's sign to keep
/// `det = +1`, each column normalised once at the end by a Newton step. The exit compares values
/// that every target computes to the bit (D16). The matrix is first scaled by an exact power of two
/// into a window (`[2^-128, 2^128)` at `f64`, `[2^-32, 2^32)` at `f32`), and the eigenvalues scaled
/// back.
///
/// # Accuracy
///
/// Backward stable: the eigenvalues are within a few `u |A|` of the exact ones and column `k` within
/// a few `u |A| / gap_k` radians, `gap_k` its distance to the nearest other eigenvalue, at every
/// gap, ties and a multiple of the identity included. The corpus id `eig3` reads at most 3.06 `u`
/// (`f64`) and 2.67 `u` (`f32`) over every stratum (`0084`, draft); the bounds `eig3_tests`
/// asserts are its header's.
///
/// # Domain
///
/// Every input is legal and nothing is asserted (`docs/API.md` R6, `0081`): only the lower triangle
/// is read, and the caller owns symmetry. For every finite input the eigenvalues are finite unless
/// one exceeds the largest finite number, and the frame is orthonormal to a few `u`. A matrix whose
/// largest entry is subnormal keeps only the digits it has. A non-finite entry gives `NaN`
/// eigenvalues, never a plausible number, and the identity frame.
///
/// A `Dual` result differentiates the rotations taken, the last step's included: the exit reads
/// values only, and that step annihilates the lanes of the entries it leaves (a diagonal input with
/// an off-diagonal lane gets its vectors' derivative there). Where the eigenvalues are simple that
/// is first-order perturbation theory, to the convergence of the exit. At a repeated eigenvalue
/// the eigenvalues are not differentiable, and the lanes are finite but in general no directional
/// derivative: at a multiple of the identity, a diagonal lane gives each eigenvalue one diagonal
/// entry's derivative, and an off-diagonal lane gives the eigenvalues none and the vectors a
/// rotation.
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
    let finite = e[1..]
        .iter()
        .fold(is_finite(e[0]), |ok, &x| ok.and(is_finite(x)));
    let e = e.map(|x| S::select(finite, x, zero));
    let (lo, hi, steps) = match S::PRECISION {
        Precision::F64 => WINDOW_F64,
        Precision::F32 => WINDOW_F32,
    };
    let m = largest(&e);
    let inside = S::lit(lo).le(m).and(m.lt(S::lit(hi))).or(m.le(zero));
    let ([d0, d1, d2, o01, o02, o12], m, back) =
        S::branch(inside, || (e, m, one), || ladder(&steps, (e, m, one)));

    let two = S::lit(2.0);
    let ident = [
        Vector([one, zero, zero]),
        Vector([zero, one, zero]),
        Vector([zero, zero, one]),
    ];
    // The first rotation, of `(0, 1)`, is general: `a_rp = a_02` is not yet zero.
    let w = two * o01;
    let mut st = turn::<S, true>(([d0, d1, d2, d1 - d0, o01, w * w, o12], ident), o02);
    // The exit (`0084`, draft): both entries left at most `u |A|`, `|A|` the largest entry.
    let th = S::lit(match S::PRECISION {
        Precision::F64 => const { pow2(-53) },
        Precision::F32 => const { pow2(-24) },
    }) * m;
    let mut k = 1;
    while k < ROTATIONS {
        let calm = st.0[4].abs().le(th).and(st.0[6].abs().le(th));
        st = S::branch(calm, || st, || turn::<S, false>(st, zero));
        // A scalar leaves at its exit; lanes run on, each one held by the branch above.
        k = calm.decide(|| ROTATIONS, || k + 1, |_, _, next| next);
    }
    let [(d0, v0), (d1, v1), (d2, v2)] = polish(st);

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
    // Each column is within about `17 u` of unit length (measured), so one Newton step on
    // `1 / sqrt(n^2)`, `(3 - n^2) / 2`, normalises it to rounding with no root and no division
    // (`0084`, draft).
    let (three_halves, half) = (S::lit(1.5), S::lit(0.5));
    let unit = |v: Vec3<S>| v.scale(three_halves - half * v.norm_sq());
    (
        Vector([value(p0.0), value(p1.0), value(p2.0)]),
        Mat3::from_cols([unit(p0.1.scale(sign)), unit(p1.1), unit(p2.1)]),
    )
}
