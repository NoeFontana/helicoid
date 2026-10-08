//! Real roots of a cubic (`docs/PHASE2.md` §6, decision `0017`).
//!
//! Depressed-cubic solver generic over `S: Real`: Cardano for one real root, Viète's
//! trigonometric form for three, and explicit repeated-root formulas. Branches evaluate
//! at safe arguments using `Real::branch` and mask decisions. Valid roots report via mask.

use crate::real::{is_finite, Mask, Precision, Real};
use crate::vector::{Vec3, Vector};

/// `2^e` as an `f64`.
///
/// # Domain
///
/// `-1022 <= e <= 1023`, asserted; exact at both precisions while `-126 <= e <= 127`. Outside it the
/// shift would discard the sign-extended high bits and return an arbitrary finite float, so the
/// assert stands rather than a `debug_assert!` (D11 governs release *checks*, and every caller
/// evaluates this in a `const` context, where the assert is a compile error and costs nothing at
/// run time).
const fn pow2(e: i32) -> f64 {
    assert!(
        e >= -1022 && e <= 1023,
        "pow2: exponent outside the normal range"
    );
    f64::from_bits(((1023 + e) as u64) << 52)
}

/// `2^ULPS_LOG2` unit roundoffs: `2^(ULPS_LOG2 - 53)` for `f64`, `2^(ULPS_LOG2 - 24)` for `f32`.
///
/// `ULPS_LOG2` is a const parameter and each arm is a `const` block, so the value is a literal in
/// the compiled code and `pow2`'s domain is checked when this instantiates, not per call.
#[inline]
fn tol<S: Real, const ULPS_LOG2: i32>() -> S {
    match S::PRECISION {
        Precision::F64 => S::lit(const { pow2(ULPS_LOG2 - 53) }),
        Precision::F32 => S::lit(const { pow2(ULPS_LOG2 - 24) }),
    }
}

/// `2^13 u`, the relative band of the discriminant.
const DISCRIMINANT_BAND: i32 = 13;

/// `max(x, y)` for finite arguments.
#[inline]
fn max<S: Real>(x: S, y: S) -> S {
    S::select(x.lt(y), y, x)
}

/// [`homogenise`]'s window edges `L`, `H` and their squares and cubes: `s >= L` is `|B| >= L` or
/// `|C| >= L^2` or `|D| >= L^3`, and `s < H` the same with `<` and `and`.
const WINDOW_F64: ([f64; 3], [f64; 3]) = (
    [pow2(-128), pow2(-256), pow2(-384)],
    [pow2(128), pow2(256), pow2(384)],
);
const WINDOW_F32: ([f64; 3], [f64; 3]) = (
    [pow2(-7), pow2(-14), pow2(-21)],
    [pow2(16), pow2(32), pow2(48)],
);

/// `2^k` for the two-sided steps of [`homogenise`], down to the edge of its window: after them the
/// scale lies in `[2^-127, 2^128)` (`f64`) or `[2^-7, 2^8)` (`f32`).
const STEPS_F64: [f64; 3] = [pow2(512), pow2(256), pow2(128)];
const STEPS_F32: [f64; 4] = [pow2(64), pow2(32), pow2(16), pow2(8)];

/// One two-sided step of [`homogenise`] at `K = 2^k`: where the scale `s = max(|B|, sqrt|C|,
/// cbrt|D|)` is at least `K`, divide the cubic by it (`B/K`, `C/K^2`, `D/K^3`); where it is below
/// `2/K`, multiply. Each comparison is on an exact power-of-two multiple of a coefficient, so
/// nothing is rounded and no `sqrt` or `cbrt` is taken; a product that underflows or overflows
/// compares on the side it belongs to.
#[inline]
fn step<S: Real>(k: f64, (b, c, d, m): (S, S, S, S)) -> (S, S, S, S) {
    let (one, two) = (S::one(), S::lit(2.0));
    let (big, small) = (S::lit(k), S::lit(1.0 / k));
    let (ab, ac, ad) = (b.abs(), c.abs(), d.abs());
    let down = big
        .le(ab)
        .or(big.le(ac * small))
        .or(big.le(ad * small * small));
    let up = (ab * big)
        .lt(two)
        .and((ac * big * big).lt(two * two))
        .and((ad * big * big * big).lt(two * two * two));
    let f = S::select(down, small, S::select(up, big, one));
    let g = S::select(down, big, S::select(up, small, one));
    (b * f, c * f * f, d * f * f * f, m * g)
}

/// The cubic `y^3 + B y^2 + C y + D` homogenised to `x = m y` (`0031` L2): `(B/m, C/m^2, D/m^3, m)`
/// with `m` a power of two, every value exact.
///
/// What the scaling owes is that neither summand of `disc` over- or underflows: `p` is of the
/// order `s^2` or, cancelling, at least `u s^2`, and `q` of `s^3`, so a scale `s` in `[2^-128,
/// 2^128)` (`f64`) or `[2^-7, 2^16)` (`f32`) keeps `p^3` and `q^2` normal, and the one tolerance of
/// the arms, the band, compares quantities of one degree in `s`. Inside that window `m = 1` and the cubic is
/// untouched; outside it the ladder brings `s` in, as a lazy branch, so a scalar cubic in the window
/// pays three comparisons per coefficient. Below the smallest normal, and at `B = C = D = 0`, the
/// scale stays below the window (the ladders sum to `2^896`, `2^120`), which leaves the arithmetic
/// exact.
// `always`: left to LLVM it is a call returning four values through memory, measured 10-33 ns of
// a 14-45 ns solver; inline, the window test is a few comparisons.
#[inline(always)]
fn homogenise<S: Real>(b: S, c: S, d: S) -> (S, S, S, S) {
    let (steps, (lo, hi)): (&[f64], ([f64; 3], [f64; 3])) = match S::PRECISION {
        Precision::F64 => (&STEPS_F64, WINDOW_F64),
        Precision::F32 => (&STEPS_F32, WINDOW_F32),
    };
    let ([lo1, lo2, lo3], [hi1, hi2, hi3]) = (lo.map(S::lit), hi.map(S::lit));
    let (ab, ac, ad) = (b.abs(), c.abs(), d.abs());
    let inside = ab
        .lt(hi1)
        .and(ac.lt(hi2))
        .and(ad.lt(hi3))
        .and(lo1.le(ab).or(lo2.le(ac)).or(lo3.le(ad)));
    let start = (b, c, d, S::one());
    S::branch(inside, || start, || ladder(steps, start))
}

/// [`homogenise`]'s steps, out of line: rarely taken, and a call is an arm LLVM cannot speculate
/// into the common path, which it does with the arithmetic inline.
#[cold]
#[inline(never)]
fn ladder<S: Real>(steps: &[f64], start: (S, S, S, S)) -> (S, S, S, S) {
    steps.iter().fold(start, |x, &k| step(k, x))
}

/// `|a|`'s window for `1/a`: inside it `1/a` is finite and normal, outside it all four coefficients
/// are scaled by one power of two (`0031` L3), which changes no root.
const LEADING_F64: (f64, f64, [f64; 2]) = (pow2(-1020), pow2(1020), [pow2(512), pow2(256)]);
const LEADING_F32: (f64, f64, [f64; 2]) = (pow2(-124), pow2(124), [pow2(64), pow2(32)]);

/// `(a, b, c, d)` times a power of two that brings `|a|` into the window of `1/a`, exact: the same
/// cubic, so neither a subnormal `a` nor one near the largest finite number loses a cubic to `1/a`
/// overflowing or going subnormal. A lazy branch around a cold ladder, as [`homogenise`]'s.
#[inline(always)]
fn lead<S: Real>(a: S, b: S, c: S, d: S) -> (S, S, S, S) {
    let (lo, hi, steps) = match S::PRECISION {
        Precision::F64 => LEADING_F64,
        Precision::F32 => LEADING_F32,
    };
    let aa = a.abs();
    let inside = S::lit(lo).le(aa).and(aa.lt(S::lit(hi)));
    S::branch(
        inside,
        || (a, b, c, d),
        || lead_ladder(&steps, (a, b, c, d)),
    )
}

/// [`lead`]'s steps, two-sided on `|a|`, out of line for [`ladder`]'s reason.
#[cold]
#[inline(never)]
fn lead_ladder<S: Real>(steps: &[f64], start: (S, S, S, S)) -> (S, S, S, S) {
    let (one, two) = (S::one(), S::lit(2.0));
    steps.iter().fold(start, |(a, b, c, d), &k| {
        let (big, small) = (S::lit(k), S::lit(1.0 / k));
        let aa = a.abs();
        let f = S::select(big.le(aa), small, S::select((aa * big).lt(two), big, one));
        (a * f, b * f, c * f, d * f)
    })
}

/// `pi`, correctly rounded at this precision; a `Dual` sees a constant.
///
/// The literal per precision, not `atan2(+0, -1)`: that spelling is bit-equal at both precisions
/// (`pi_and_acos_are_within_their_ulps` asserts it against `core::f64::consts::PI` and
/// `core::f32::consts::PI`) but is a `libm` call for a compile-time constant, and `libm::atan2` is
/// neither generic nor `#[inline]`, so without workspace LTO it cannot be folded away. `f32`'s
/// value is widened exactly into `lit`'s `f64`, never the `f64` constant rounded down to `f32`,
/// which would be a double rounding.
#[inline]
pub(crate) fn pi<S: Real>() -> S {
    match S::PRECISION {
        Precision::F64 => S::lit(core::f64::consts::PI),
        Precision::F32 => S::lit(core::f32::consts::PI as f64),
    }
}

/// The real roots of `a x^3 + b x^2 + c x + d`, and which slots of the result are roots.
///
/// Returns `(roots, valid)`: slot `k` is reported as a root iff `valid[k]` is set (`# Domain` says
/// where a reported root is not one). A slot whose mask is clear holds exactly `+0` (never NaN or
/// infinity) and is not a root, even when it reads like one. A root that is exactly zero is `+0`,
/// never `-0`. The order of the slots is that of the arm that produced them; it is not sorted:
/// three roots come as the phases `k = 0, 1, 2`, largest first up to rounding.
///
/// With `B = b/a`, `C = c/a`, `D = d/a` (`1/a` and three products, after `a`, `b`, `c`, `d` are
/// scaled by one power of two into `1/a`'s range), the cubic is homogenised (`0031` L2): `x = m y`,
/// `m` a power of two, gives the monic cubic in `y` with coefficients `B/m`, `C/m^2`, `D/m^3`, formed
/// exactly, whose scale `max(|B|, sqrt|C|, cbrt|D|)` keeps `p^3` and `q^2` normal (`m = 1` for a
/// scale in `[2^-128, 2^128)`, `[2^-7, 2^16)` at `f32`). On those, `y = t - B/(3m)` gives
/// `t^3 + p t + q`, `p = C - B^2/3`, `q = 2 B^3/27 - B C/3 + D` (scaled), and
/// `disc = q^2/4 + p^3/27` (positive for one real root). With
/// `band = max(q^2/4, |p^3/27|) 2^-40` (`f64`), the arms are, in this order:
///
/// | arm | taken when | valid slots and their values (each `- B/(3m)`, then times `m`) |
/// |---|---|---|
/// | one real root | `disc > band` | 0: `w + v` where `p <= 0`, `-q / (w^2 + p/3 + v^2)` where `p > 0`; `w = cbrt(-q/2 + copysign(sqrt disc, -q/2))`, `v = -(p/3)/w` |
/// | three real roots | `disc < -band` or `p < 0` | 0, 1, 2: `2 r cos(acos(3q / (2 p r)) / 3 - 2 pi k / 3)`, `r = sqrt(-p/3)`, the `acos` argument clamped to `[-1, 1]` |
/// | triple root | else | 0: `cbrt(-q)` |
///
/// The one-real-root arm pairs its two cube roots through `w v = -p/3` (`0031` L1): `w` is the one
/// whose radicand is a sum of like signs, so it never cancels, and `v = -(p/3)/w`. Where `p > 0` the
/// two have opposite signs and their sum cancels for a root small against the scale, so the root
/// is taken from `w^3 + v^3 = -q` instead, over `w^2 - w v + v^2`, a sum of positive terms. Both
/// are formed and one selected, without a branch.
///
/// The triple-root arm is the only other case: with `p >= 0` both summands of `disc` are
/// non-negative, so `disc` lies within `band` of zero only when both are `+0`, which on the
/// homogenised cubic is `p` and `q` zero or below about `2^-340` (`f64`), a triple root of the
/// scale. A slot whose value is not finite is not valid: a root beyond the largest finite number
/// after the `m` scaling.
///
/// # Tolerances
///
/// The one threshold is a power of two, exact at both precisions, a fixed multiple of the unit
/// roundoff `u` (`2^-53`, `2^-24`) (`0017`), relative to the two summands of `disc` (`0031`):
///
/// | threshold | multiple | `f64` | `f32` |
/// |---|---|---|---|
/// | relative discriminant band | `2^13 u` | `2^-40` | `2^-11` |
///
/// Changing it is a changelog line naming this function.
///
/// # Domain
///
/// Every input is legal and nothing is asserted (`docs/API.md` R4, R6): the mask is the report. A
/// non-finite coefficient, `a = 0`, and one of `b/a`, `c/a`, `d/a` overflowing are not a cubic and
/// give no valid slot; there is no fallback to a quadratic, as `a -> 0` sends a root to infinity
/// (`0031` L3). Scaling every coefficient by the same power of two changes no root and no mask, at
/// both precisions, for every finite nonzero `a` and down to a scale `max(|B|, sqrt|C|, cbrt|D|)` at
/// the smallest normal.
///
/// The formulas are not backward stable near a multiple root. A double root is found to about
/// `sqrt(u)` and a triple root to `u^(1/3)` of the scale, and a complex pair whose imaginary part is
/// below about `sqrt(band)` of the scale is reported as a real double root (the price of never
/// dropping a repeated root).
///
/// The band is relative to the two summands of `disc`, not to `B`, `C`, `D`. Where `p` and `q`
/// cancel (roots close together, far from the origin) the rounding error of `disc` can exceed it, a
/// pair of real roots reads as positive and the arm has one root: `(x - 1.09375)^2 (x - 0.921875)`
/// at `f32` loses its double root (`0031` L4).
///
/// A `Dual` root differentiates the arm taken, the mask reading the value part. The derivative is
/// infinite or NaN where the true one is (a multiple root) and where the `acos` argument is `+-1`
/// (`sqrt` at 0). In the one-real-root arm it is finite at every simple root.
///
/// # Example
///
/// ```
/// use helicoid_linalg::solve_cubic;
///
/// // (x - 1)(x - 2)(x - 3), and x^3 + x^2 + x + 1 = (x + 1)(x^2 + 1), whose only real root is -1.
/// let (roots, valid) = solve_cubic(1.0_f64, -6.0, 11.0, -6.0);
/// assert!(valid.iter().all(|&v| v));
/// let mut sorted = roots.0;
/// sorted.sort_by(f64::total_cmp);
/// assert!(sorted.iter().zip([1.0, 2.0, 3.0]).all(|(r, w)| (r - w).abs() < 1e-12));
///
/// let (roots, valid) = solve_cubic(1.0_f64, 1.0, 1.0, 1.0);
/// assert_eq!(valid, [true, false, false]);
/// assert!((roots.0[0] + 1.0).abs() < 1e-12 && roots.0[1] == 0.0);
/// ```
#[inline]
pub fn solve_cubic<S: Real>(a: S, b: S, c: S, d: S) -> (Vec3<S>, [S::Mask; 3]) {
    let (zero, one) = (S::zero(), S::one());
    let (two, three) = (S::lit(2.0), S::lit(3.0));

    // `1/a` and three products, not three divisions (`0031` L5, measured against), after `a` is
    // brought into `1/a`'s range.
    let usable = is_finite(a).and(zero.lt(a.abs()));
    let (a, b, c, d) = lead(S::select(usable, a, one), b, c, d);
    let inv_a = one / a;
    let (big_b, big_c, big_d) = (b * inv_a, c * inv_a, d * inv_a);
    let is_cubic = usable
        .and(is_finite(big_b))
        .and(is_finite(big_c))
        .and(is_finite(big_d));
    // Where there is no cubic every lane evaluates `x^3 = 0`, so no arm sees a value it cannot take.
    let (big_b, big_c, big_d) = (
        S::select(is_cubic, big_b, zero),
        S::select(is_cubic, big_c, zero),
        S::select(is_cubic, big_d, zero),
    );
    let (big_b, big_c, big_d, m) = homogenise(big_b, big_c, big_d);

    let shift = big_b / three;
    let p = big_c - big_b * big_b / three;
    let q = two * big_b * big_b * big_b / S::lit(27.0) - big_b * big_c / three + big_d;
    let term_q = q * q * S::lit(0.25);
    let term_p = p * p * p / S::lit(27.0);
    let disc = term_q + term_p;
    let band = max(term_q, term_p.abs()) * tol::<S, DISCRIMINANT_BAND>();

    let one_root = band.lt(disc);
    let three_roots = one_root.not().and(disc.lt(-band).or(p.lt(zero)));

    let t: [S; 3] = S::branch(
        one_root,
        || {
            // `w v = -p/3` (`0031` L1). An inactive lane has `disc = 1`, so `|w| >= 1`; an active
            // one has `disc > band >= 0`, so `w != 0` either way and `v` is safe. Where `p > 0`, `w`
            // and `v` have opposite signs and `w + v` cancels when the root is small against the
            // scale, so the root is `-q / (w^2 - w v + v^2)`, a sum of positive terms: the same
            // root, from `w^3 + v^3 = -q`. Its denominator is at least `|p|/3` in every lane, so
            // both are formed and one selected: a `p` of either sign costs no misprediction.
            let root_disc = S::select(one_root, disc, one).sqrt();
            let half = S::lit(-0.5) * q;
            let third = p / three;
            let w = (half + root_disc.copysign(half)).cbrt();
            let v = -third / w;
            let quotient = -q / (w * w + third + v * v);
            [S::select(zero.lt(p), quotient, w + v), zero, zero]
        },
        || {
            S::branch(
                three_roots,
                || {
                    let p = S::select(three_roots, p, -three);
                    let r = (-p / three).sqrt();
                    let arg = (three * q) / (two * p * r);
                    let arg = S::select(arg.lt(-one), -one, arg);
                    let arg = S::select(one.lt(arg), one, arg);
                    // `theta / three`, `2 pi / 3` and `2 r` are loop invariants. Hoisting them is
                    // bit-identical: `2 pi k / 3` for `k = 0, 1` is unchanged, and for `k = 2` the
                    // two spellings differ only by a factor of two, which rounding commutes with.
                    let third = arg.acos() / three;
                    let step = two * pi::<S>() / three;
                    let two_r = two * r;
                    // `Real::cos` and not `sin_cos().1`: the sine is discarded, and the two are
                    // bit-identical (`0052`, `0022`), so this drops three kernels and no bit.
                    core::array::from_fn(|k| two_r * (third - step * S::lit(k as f64)).cos())
                },
                || [(-q).cbrt(), zero, zero],
            )
        },
    );

    let active = [
        is_cubic,
        is_cubic.and(three_roots),
        is_cubic.and(three_roots),
    ];
    let slot = |k: usize| {
        // Times `m`, exact: the root of the cubic as given.
        let alpha = (t[k] - shift) * m;
        let valid = active[k].and(is_finite(alpha));
        // `+ 0` turns `-0` into `+0` and leaves every other value, and its derivative, alone.
        (S::select(valid, alpha + zero, zero), valid)
    };
    let ((r0, v0), (r1, v1), (r2, v2)) = (slot(0), slot(1), slot(2));
    (Vector([r0, r1, r2]), [v0, v1, v2])
}
