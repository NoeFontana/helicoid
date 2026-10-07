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
/// `-1022 <= e <= 0`, asserted; exact at both precisions while `e >= -126`. Outside it the shift
/// would discard the sign-extended high bits and return an arbitrary finite float, so the assert
/// stands rather than a `debug_assert!` (D11 governs release *checks*, and every caller evaluates
/// this in a `const` block, where the assert is a compile error and costs nothing at run time).
const fn pow2(e: i32) -> f64 {
    assert!(
        e >= -1022 && e <= 0,
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

/// `2^7 u`, the floor on the leading coefficient and the triple-root test.
const NEAR_ZERO: i32 = 7;
/// `2^13 u`, the relative band of the discriminant.
const DISCRIMINANT_BAND: i32 = 13;

/// `max(x, y)` for finite arguments.
#[inline]
fn max<S: Real>(x: S, y: S) -> S {
    S::select(x.lt(y), y, x)
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
/// With `B = b/a`, `C = c/a`, `D = d/a` the substitution `x = t - B/3` gives `t^3 + p t + q`,
/// `p = C - B^2/3`, `q = 2 B^3/27 - B C/3 + D`, and `disc = q^2/4 + p^3/27` (positive for one real
/// root). With `band = max(q^2/4, |p^3/27|) 2^-40` (`f64`), the arms are, in this order:
///
/// | arm | taken when | valid slots and their values (each `- B/3`) |
/// |---|---|---|
/// | one real root | `disc > band` | 0: `cbrt(-q/2 + sqrt disc) + cbrt(-q/2 - sqrt disc)` |
/// | three real roots | `disc < -band` or `p < 0` | 0, 1, 2: `2 r cos(acos(3q / (2 p r)) / 3 - 2 pi k / 3)`, `r = sqrt(-p/3)`, the `acos` argument clamped to `[-1, 1]` |
/// | triple root | else, and `\|p\| <= max(\|p\|, \|q\|, 1) 2^-46` | 0: `cbrt(-q)` |
/// | single and double root | else | 0: `3q/p`; 1: `-3q/(2p)`, the double root, once |
///
/// The last arm is unreachable for coefficients that pass the floor: with `p >= 0` both summands of
/// `disc` are non-negative, so it lies within `band` of zero only when both underflow to `+0`, and
/// the triple-root test takes that. A slot whose value is not finite is not valid: `2 p r`
/// underflows for `x^3 - 1e-250 x`, the `acos` argument is NaN, and no root is reported although
/// there are three.
///
/// # Tolerances
///
/// Three thresholds are powers of two, exact at both precisions, a fixed multiple of the unit
/// roundoff `u` (`2^-53`, `2^-24`) (`0017`):
///
/// | threshold | multiple | `f64` | `f32` |
/// |---|---|---|---|
/// | leading coefficient floor: `max(scale, 1) tol` | `2^7 u` | `2^-46` | `2^-17` |
/// | relative discriminant band | `2^13 u` | `2^-40` | `2^-11` |
/// | relative triple-root test | `2^7 u` | `2^-46` | `2^-17` |
///
/// Changing one is a changelog line naming this function.
///
/// # Domain
///
/// Every input is legal and nothing is asserted (`docs/API.md` R4, R6): the mask is the report. A
/// non-finite coefficient, `a = 0` and a leading coefficient below its floor are not a cubic and
/// give no valid slot; there is no fallback to a quadratic, as `a -> 0` sends a root to infinity.
/// The floor also bounds each of `|B|`, `|C|`, `|D|` (`b/a`, `c/a`, `d/a`) by `1/tol`, `2^17` at
/// `f32` and `2^46` at `f64`: at `f32`, `x^3 - 10^5 x` is solved and `x^3 - 10^6 x` (roots `0`,
/// `+-1000`) and `(x - 100)(x - 101)(x + 99)` (`D = 999900`) are not.
///
/// The formulas are not backward stable. Near a multiple root they cancel: a double root is found
/// to about `sqrt(u)` and a triple root to `u^(1/3)` of the scale, and a complex pair whose
/// imaginary part is below about `sqrt(band)` of the scale is reported as a real double root (the
/// price of never dropping a repeated root). In the one-real-root arm the root is a difference of
/// two cube roots, one of them of a `disc` that cancels when `p^3` is small against `q^2`: for
/// `x^3 + p x - 1` the root is off by up to `4e-6` at `p` near `1e-5` (`f64`), `4e-3` at `p` near
/// `1e-2` (`f32`), and `p = 0` is exact.
///
/// The band is relative to the two summands of `disc`, not to `B`, `C`, `D`. Where `p` and `q`
/// cancel (roots close together, far from the origin) the rounding error of `disc` can exceed it, a
/// pair of real roots reads as positive and the arm has one root: `(x - 1.09375)^2 (x - 0.921875)`
/// at `f32` loses its double root. Where `p^3` and `q^2` underflow (roots of the size `sqrt|p|`
/// below about `4e-8` at `f32`, `2e-54` at `f64`) `disc` reads 0 and `p < 0` sends every cubic to
/// the trigonometric arm: with one real root the three valid slots are not roots
/// (`x^3 - 10^-16 x + 10^-24` at `f32`).
///
/// A `Dual` root differentiates the arm taken, the mask reading the value part. The derivative is
/// infinite or NaN where the true one is (a multiple root) and where the `acos` argument is `+-1`
/// (`sqrt` at 0). In the one-real-root arm it is `-inf` wherever `h - s` cancels to exactly 0,
/// although the root is simple: at `p = 0` (`x^3 + q`), where `cbrt` is evaluated at 0
/// (`Real::cbrt`), and for `x^3 + p x - 1` at `f64` for every `p` below about `1.3e-5`, where the
/// value is off by up to `4e-6`; above that it loses accuracy faster than the value.
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

    let abs_a = a.abs();
    let scale = max(max(abs_a, b.abs()), max(c.abs(), d.abs()));
    let is_cubic = is_finite(a)
        .and(is_finite(b))
        .and(is_finite(c))
        .and(is_finite(d))
        .and(abs_a.le(max(scale, one) * tol::<S, NEAR_ZERO>()).not());
    // Where there is no cubic every lane evaluates `x^3 = 0`, so no arm sees a value it cannot take.
    let (a, b, c, d) = (
        S::select(is_cubic, a, one),
        S::select(is_cubic, b, zero),
        S::select(is_cubic, c, zero),
        S::select(is_cubic, d, zero),
    );

    let inv_a = one / a;
    let (big_b, big_c, big_d) = (b * inv_a, c * inv_a, d * inv_a);
    let shift = big_b / three;
    let p = big_c - big_b * big_b / three;
    let q = two * big_b * big_b * big_b / S::lit(27.0) - big_b * big_c / three + big_d;
    let term_q = q * q * S::lit(0.25);
    let term_p = p * p * p / S::lit(27.0);
    let disc = term_q + term_p;
    let band = max(term_q, term_p.abs()) * tol::<S, DISCRIMINANT_BAND>();
    let p_scale = max(max(p.abs(), q.abs()), one);

    let one_root = band.lt(disc);
    let three_roots = one_root.not().and(disc.lt(-band).or(p.lt(zero)));
    let repeated = one_root.not().and(three_roots.not());
    let triple = p.abs().le(p_scale * tol::<S, NEAR_ZERO>());
    // The fourth arm is unreachable for coefficients that pass the floor, and this is where that
    // claim is executable rather than prose: `repeated` forces `p >= 0`, so both summands of `disc`
    // are non-negative and `fl(term_q + term_p) >= max(term_q, term_p)`, so `disc <= band` only
    // when both are `+0`, which `triple` then takes. It holds *because* the floor bounds `|B|`,
    // `|C|`, `|D|`, which keeps either summand from overflowing to `+inf` (there `band` is `+inf`
    // too and nothing is `one_root` or `three_roots`). Lowering that floor is the change that would
    // make the arm live again, and `0031` (draft) L3 recommends exactly that, so this fires first.
    debug_assert!(
        !repeated.and(triple.not()).any(),
        "solve_cubic: the single-and-double-root arm became reachable; see `0031` (draft) L3"
    );

    let t: [S; 3] = S::branch(
        one_root,
        || {
            let root_disc = S::select(one_root, disc, one).sqrt();
            let half = S::lit(-0.5) * q;
            [
                (half + root_disc).cbrt() + (half - root_disc).cbrt(),
                zero,
                zero,
            ]
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
                || {
                    S::branch(
                        triple,
                        || [(-q).cbrt(), zero, zero],
                        || {
                            let p = S::select(triple, one, p);
                            [three * q / p, S::lit(-1.5) * q / p, zero]
                        },
                    )
                },
            )
        },
    );

    let active = [
        is_cubic,
        is_cubic.and(three_roots.or(repeated.and(triple.not()))),
        is_cubic.and(three_roots),
    ];
    let slot = |k: usize| {
        let alpha = t[k] - shift;
        let valid = active[k].and(is_finite(alpha));
        // `+ 0` turns `-0` into `+0` and leaves every other value, and its derivative, alone.
        (S::select(valid, alpha + zero, zero), valid)
    };
    let ((r0, v0), (r1, v1), (r2, v2)) = (slot(0), slot(1), slot(2));
    (Vector([r0, r1, r2]), [v0, v1, v2])
}
