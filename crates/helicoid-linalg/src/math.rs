//! Portable elementary functions: each is one fixed sequence of IEEE 754 basic operations, so its
//! bits are the same on every target, with or without a fused multiply-add (D16, `0105`).
//!
//! [`atan2`] and [`atan2f`] are what [`Real::atan2`](crate::Real::atan2) computes at `f64` and
//! `f32`. Their polynomials and tables are `generated.rs`, written by `just atan-fit` and checked
//! by `just atan-fit-check`; nothing here is a typed coefficient.
//!
//! **binary64.** `t = n/d ∈ [0, 1]`, `n = min`, `d = max` of `|x|, |y|`, is reduced about the
//! nearest of the `CENTRES + 1` centres `c = j/C`: `atan(t) = atan(c) + atan(v)`,
//! `v = (n − c·d)/(d + c·n)`, `|v| ≤ 1/(2C)`. `c` has at most `log₂ C + 1` significant bits, so
//! `c·d` and `c·n` are each two exact products once the low bits of `d` and `n` are split off;
//! `n − c·d` is then exact (Sterbenz, then a representable difference) and `d + c·n` is an exact
//! pair. The quotient's residual comes from Dekker's halves, without a fused multiply-add. The tail
//! `v·s·P(s)` is at most `v²/3 ≤ 2⁻⁹·⁶` of the answer, so it needs no compensation.
//!
//! **binary32.** Evaluated in binary64: folded at `t = ½` by `(n − d)/(n + d)`, whose numerator
//! and denominator are exact for binary32 operands, then one rounding of the quotient and a
//! polynomial good to `2⁻⁴⁰·⁶`, then one rounding to binary32.

mod generated;

use generated::{ATAN_HI, ATAN_LO, CENTRES, P32, P64, PI4_HI, PI4_LO};

/// `atan2(y, x)`: the angle of `(x, y)` from the positive `x` axis, in `[−π, π]`.
///
/// Within 0.504 ulp of the exact value on every draw read (`0105`: at most 0.5032 ulp over
/// `3 × 2²²` draws with `|y/x|` within `2^±4`, `2^±60` and `2^±1000`), and IEEE 754's limits,
/// correctly rounded, where an argument is zero or infinite. Odd in `y`, and unchanged when both
/// arguments are scaled by a power of two that keeps them normal. The same bits on every target.
///
/// # Domain
///
/// Every input; a NaN argument gives a NaN.
///
/// ```
/// use helicoid_linalg::math::atan2;
///
/// assert_eq!(atan2(1.0, 1.0).to_bits(), core::f64::consts::FRAC_PI_4.to_bits());
/// assert_eq!(atan2(-0.0, -1.0).to_bits(), (-core::f64::consts::PI).to_bits());
/// ```
#[must_use]
pub fn atan2(y: f64, x: f64) -> f64 {
    let (ax, ay) = (x.abs(), y.abs());
    let swap = ax < ay;
    let n = if ay < ax { ay } else { ax };
    let d = if ax > ay { ax } else { ay };
    let (n, d) = if WINDOW.0 <= d && d <= WINDOW.1 {
        (n, d)
    } else {
        rescale(n, d)
    };
    let (j, num, den, den_lo) = centre(n, d);
    // Two divisions, this and the residual's `1/den`: `r = 1/den`, `v = num·r` reads the same
    // accuracy but puts a multiply on the chain, 28.0 against 27.3 ns, and the divider is not the
    // throughput bound (`0105`).
    let v = num / den;

    // The quadrant: `q·π/4 ± (atan(c) + atan(v))`, `−` where exactly one of `x < 0`, `|y| > |x|`.
    let neg_x = x.is_sign_negative();
    let sign = if neg_x != swap { -1.0 } else { 1.0 };
    let q = if swap {
        2.0
    } else if neg_x {
        4.0
    } else {
        0.0
    };
    // `q·π/4_hi` is exact; `|q·π/4| ≥ atan(c)` and `|q·π/4 ± atan(c)| ≥ |v|` unless zero.
    let (h0, l0) = fast_two_sum(q * PI4_HI, sign * ATAN_HI[j]);
    let us = sign * v;
    let (h1, l1) = fast_two_sum(h0, us);
    let s = v * v;
    let [k0, k1, k2, k3, k4] = P64;
    let s2 = s * s;
    let p = (k0 + k1 * s) + s2 * ((k2 + k3 * s) + s2 * k4);
    // A branch, not a select: in one basic block LLVM's SLP vectorizer packs the residual's
    // products into lanes and lengthens the chain by a third (`0105`).
    let dv = if RESIDUAL <= num.abs() {
        // `num − v·(den + den_lo)`: `vh·bh` and `vh·bl` are exact and `num − vh·bh` is exact
        // (Sterbenz); the rest rounds at `2⁻²⁵` of the residual.
        let bc = SPLIT * den;
        let bh = bc - (bc - den);
        let bl = den - bh;
        let vh = f64::from_bits(v.to_bits() & !((1 << 27) - 1));
        let vl = v - vh;
        let rem = ((num - vh * bh) - vh * bl) - (vl * den + v * den_lo);
        // `atan(v + dv) = atan(v) + dv/(1 + s)`, `1/(1 + s)` to `s² ≤ 2⁻¹⁶`.
        sign * (rem * ((1.0 / den) * (1.0 - s)))
    } else {
        0.0
    };
    let early = (l1 + l0) + (q * PI4_LO + sign * ATAN_LO[j]);
    let mag = h1 + ((early + dv) + (us * s) * p);
    if x.is_nan() || y.is_nan() {
        x + y
    } else {
        mag.copysign(y)
    }
}

/// `atan2(y, x)` at binary32: the angle of `(x, y)` from the positive `x` axis, in `[−π, π]`.
///
/// Evaluated in binary64 and rounded once; correctly rounded on every draw read (`0105`: `2²¹`
/// draws with `|y/x|` within `2^±4` and `2^±120`), and IEEE 754's limits, correctly rounded, where
/// an argument is zero or infinite. Odd in `y`. The same bits on every target.
///
/// # Domain
///
/// Every input; a NaN argument gives a NaN.
///
/// ```
/// use helicoid_linalg::math::atan2f;
///
/// assert_eq!(atan2f(1.0, 0.0).to_bits(), core::f32::consts::FRAC_PI_2.to_bits());
/// ```
#[must_use]
pub fn atan2f(y: f32, x: f32) -> f32 {
    let (ax, ay) = (f64::from(x.abs()), f64::from(y.abs()));
    let swap = ax < ay;
    let n = if ay < ax { ay } else { ax };
    let d = if ax > ay { ax } else { ay };
    let (n, d) = if 0.0 < d && d < f64::INFINITY {
        (n, d)
    } else {
        limits(n, d)
    };
    // `n − d` and `n + d` are exact: binary32 operands within a factor of two.
    let fold = 0.5 * d < n;
    let u = (n - pick(fold, d, 0.0)) / (d + pick(fold, n, 0.0));
    let neg_x = x.is_sign_negative();
    let sign = pick(neg_x != swap, -1.0, 1.0);
    let k = pick(swap, 2.0, pick(neg_x, 4.0, 0.0)) + pick(fold, sign, 0.0);
    let us = sign * u;
    let s = u * u;
    let [c0, c1, c2, c3, c4, c5, c6, c7] = P32;
    let (s2, s4) = (s * s, (s * s) * (s * s));
    let p = ((c0 + c1 * s) + s2 * (c2 + c3 * s)) + s4 * ((c4 + c5 * s) + s2 * (c6 + c7 * s));
    let mag = k * PI4_HI + (us + ((us * s) * p + k * PI4_LO));
    if x.is_nan() || y.is_nan() {
        x + y
    } else {
        (mag as f32).copysign(y)
    }
}

/// `2^e`, for `e` a normal exponent.
const fn pow2(e: i32) -> f64 {
    f64::from_bits(((1023 + e) as u64) << 52)
}

/// The window on `d` the compensations need: `(2²⁷ + 1)·2d` finite below `2^995`.
const WINDOW: (f64, f64) = (pow2(-460), pow2(990));
/// Below this `|num|` the residual's products are not all normal; the quotient is then exact enough.
const RESIDUAL: f64 = pow2(-969);
/// Veltkamp's splitter for a 53-bit significand.
const SPLIT: f64 = pow2(27) + 1.0;

/// `(s, e)` with `s + e = a + b` exactly, for `|a| ≥ |b|` or `a = 0`.
fn fast_two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    (s, b - (s - a))
}

/// `a` where `m`, else `b`, by masks rather than a branch.
fn pick(m: bool, a: f64, b: f64) -> f64 {
    let mask = u64::from(m).wrapping_neg();
    f64::from_bits((a.to_bits() & mask) | (b.to_bits() & !mask))
}

/// `(n, d)` with `d` brought into [`WINDOW`] by an exact power of two, or, where `d` is zero or
/// infinite, IEEE 754's limits: `(0, 1)`, or each infinite argument `1` and each finite one `0`.
#[cold]
#[inline(never)]
fn rescale(n: f64, d: f64) -> (f64, f64) {
    if d == 0.0 || d == f64::INFINITY {
        limits(n, d)
    } else {
        let scale = if d > 1.0 { pow2(-512) } else { pow2(512) };
        (n * scale, d * scale)
    }
}

/// IEEE 754's limits where `d` is zero or infinite: `(0, 1)`, or each infinite argument `1` and
/// each finite one `0`. A NaN passes through.
#[cold]
#[inline(never)]
fn limits(n: f64, d: f64) -> (f64, f64) {
    if d == 0.0 {
        (0.0, 1.0)
    } else if d == f64::INFINITY {
        (if n == f64::INFINITY { 1.0 } else { 0.0 }, 1.0)
    } else {
        (n, d)
    }
}

/// `(j, num, den, den_lo)` for `0 ≤ n ≤ d` in [`WINDOW`]: `c = j/C` the centre nearest `n/d`,
/// `num = n − c·d` exactly, `den + den_lo = d + c·n` exactly.
pub(crate) fn centre(n: f64, d: f64) -> (usize, f64, f64, f64) {
    let w = 1.0 / CENTRES as f64;
    // `c` sums, as a tree, `w` for each midpoint passed; the midpoints round, which moves `|v|`'s
    // bound by an ulp. Counting them as an integer `j` puts a popcount and an integer-to-float
    // conversion ahead of the division: 39.7 against 27.3 ns at `f64` (`0105`).
    let mut parts = [0.0; CENTRES];
    for (i, part) in parts.iter_mut().enumerate() {
        let mid = (2 * i + 1) as f64 * (0.5 * w);
        *part = if d * mid < n { w } else { 0.0 };
    }
    let mut len = CENTRES;
    while len > 1 {
        len /= 2;
        for i in 0..len {
            parts[i] = parts[2 * i] + parts[2 * i + 1];
        }
    }
    let c = parts[0];
    // Exact: `c = j/C`. Off the critical path; only the tables read it.
    let j = ((c * CENTRES as f64) as usize).min(CENTRES);
    // `c` has at most `log₂ C + 1` significant bits: clearing that many low bits of `d` and `n`
    // leaves halves whose products by `c` are exact.
    let keep = !(2 * CENTRES as u64 - 1);
    let (dh, nh) = (
        f64::from_bits(d.to_bits() & keep),
        f64::from_bits(n.to_bits() & keep),
    );
    // `n − c·dh` is exact by Sterbenz (`n/d` is within `1/(2C)` of `c`), and `n − c·d` is
    // representable, so the second subtraction is exact too.
    let num = (n - c * dh) - c * (d - dh);
    let (den, e1) = fast_two_sum(d, c * nh);
    let (den, e2) = fast_two_sum(den, c * (n - nh));
    (j, num, den, e1 + e2)
}

#[cfg(test)]
pub(crate) mod fit {
    //! The generated tables, for the tests that check their shape.
    pub(crate) use super::generated::{ATAN_HI, ATAN_LO, CENTRES, PI4_HI, PI4_LO};
}
