//! The portable `atan2` candidate of `0105` (draft): its fit, its two kernels, and [`Portable`], the
//! shipped scalar with `atan2` alone replaced by them.
//!
//! Self-contained (nothing from `crate::`) so `xtask/benches/atan2.rs` compiles the same file.
//!
//! **binary64.** `t = n/d ∈ [0, 1]`, `n = min`, `d = max` of `|x|, |y|`, is reduced about the
//! nearest of the [`CENTRES`]` + 1` centres `c = j/C`: `atan(t) = atan(c) + atan(v)`,
//! `v = (n − c·d)/(d + c·n)`, `|v| ≤ 1/(2C)`. `c` has at most `log₂ C + 1` significant bits, so
//! `c·d` and `c·n` are each two exact products once the low bits of `d` and `n` are split off;
//! `n − c·d` is then exact (Sterbenz, then a representable difference) and `d + c·n` is kept as a
//! pair. The quotient's residual comes from Dekker's halves without a fused multiply-add (D16).
//! The tail `v·s·P(s)` is at most `v²/3 ≤ 2⁻⁹·⁶` of the answer, so it needs no compensation.
//!
//! **binary32.** Evaluated in binary64: folded at `t = ½` by `(n − d)/(n + d)`, whose numerator
//! and denominator are exact for binary32 operands, then one rounding of the quotient and a
//! polynomial good to `2⁻⁴⁰·⁶`, then one rounding to binary32.
//!
//! `P`, the centres' `atan` and the split of `π/4` are `conformance/fits/atan2.json`, written by
//! `just atan-fit`; nothing here is a typed coefficient.

use std::sync::OnceLock;

use helicoid_linalg::{Precision, Real};
use serde::Deserialize;

/// The binary64 kernel's centres, a power of two, and `P`'s terms at each precision: the shape the
/// code below is written for, checked against the committed fit.
pub(crate) const CENTRES: usize = 8;
pub(crate) const TERMS64: usize = 5;
pub(crate) const TERMS32: usize = 8;

/// The committed fit, as the generator wrote it.
const FIT: &str = include_str!("../../../../conformance/fits/atan2.json");

/// The binary64 kernel's polynomial, the centres' `atan` as `hi + lo`, and `π/4 = hi + lo`.
#[derive(Debug)]
pub(crate) struct Fit64 {
    pub(crate) coeffs: [f64; TERMS64],
    pub(crate) atan_hi: [f64; CENTRES + 1],
    pub(crate) atan_lo: [f64; CENTRES + 1],
    pub(crate) pi4_hi: f64,
    pub(crate) pi4_lo: f64,
}

/// The binary32 kernel's polynomial, at binary64.
#[derive(Debug)]
pub(crate) struct Fit32 {
    pub(crate) coeffs: [f64; TERMS32],
}

#[derive(Debug)]
pub(crate) struct Fits {
    pub(crate) f64: Fit64,
    pub(crate) f32: Fit32,
}

#[derive(Deserialize)]
struct FitFile {
    f64: Json64,
    f32: Json32,
}

#[derive(Deserialize)]
struct Json64 {
    centres: usize,
    coeffs: Vec<String>,
    atan_hi: Vec<String>,
    atan_lo: Vec<String>,
    pi4_hi: String,
    pi4_lo: String,
}

#[derive(Deserialize)]
struct Json32 {
    evaluated_at: String,
    fold: String,
    coeffs: Vec<String>,
}

fn bits(hex: &str) -> Result<f64, String> {
    let digits = hex.strip_prefix("0x").ok_or(format!("{hex}: no 0x"))?;
    u64::from_str_radix(digits, 16)
        .map(f64::from_bits)
        .map_err(|e| format!("{hex}: {e}"))
}

/// `values` as exactly `K` binary64s; another length is a different kernel.
fn array<const K: usize>(what: &str, values: &[String]) -> Result<[f64; K], String> {
    let parsed = values
        .iter()
        .map(|v| bits(v))
        .collect::<Result<Vec<_>, _>>()?;
    parsed
        .try_into()
        .map_err(|v: Vec<f64>| format!("{what}: {} values, the kernel takes {K}", v.len()))
}

/// The fit, parsed once.
pub(crate) fn fits() -> Result<&'static Fits, &'static str> {
    static FITS: OnceLock<Result<Fits, String>> = OnceLock::new();
    FITS.get_or_init(|| {
        let file: FitFile = serde_json::from_str(FIT).map_err(|e| format!("atan2.json: {e}"))?;
        let (a, b) = (file.f64, file.f32);
        if a.centres != CENTRES || b.evaluated_at != "f64" || b.fold != "1/2" {
            return Err(format!(
                "a fit for {} centres, binary32 at {} folded at {}; the kernels are written for \
                 {CENTRES}, f64 and 1/2",
                a.centres, b.evaluated_at, b.fold
            ));
        }
        Ok(Fits {
            f64: Fit64 {
                coeffs: array("f64 coeffs", &a.coeffs)?,
                atan_hi: array("atan_hi", &a.atan_hi)?,
                atan_lo: array("atan_lo", &a.atan_lo)?,
                pi4_hi: bits(&a.pi4_hi)?,
                pi4_lo: bits(&a.pi4_lo)?,
            },
            f32: Fit32 {
                coeffs: array("f32 coeffs", &b.coeffs)?,
            },
        })
    })
    .as_ref()
    .map_err(String::as_str)
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

/// `(c, num, den, den_lo)` for `0 ≤ n ≤ d` in [`WINDOW`]: `c = j/C` the centre nearest `n/d`,
/// `num = n − c·d` exactly, `den + den_lo = d + c·n` to `2⁻¹⁰⁶`.
pub(crate) fn centre(n: f64, d: f64) -> (f64, f64, f64, f64) {
    let w = 1.0 / CENTRES as f64;
    // `c` sums the centres passed; the midpoints round, which moves `|v|`'s bound by an ulp.
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
    (c, num, den, e1 + e2)
}

/// `atan2(y, x)` by the binary64 kernel.
pub(crate) fn atan2_f64(y: f64, x: f64, fit: &Fit64) -> f64 {
    let (ax, ay) = (x.abs(), y.abs());
    let swap = ax < ay;
    let n = if ay < ax { ay } else { ax };
    let d = if ax > ay { ax } else { ay };
    let (n, d) = if WINDOW.0 <= d && d <= WINDOW.1 {
        (n, d)
    } else {
        rescale(n, d)
    };
    let (c, num, den, den_lo) = centre(n, d);
    let v = num / den;
    let inv = 1.0 / den;
    // `num − v·(den + den_lo)`: `vh·bh` and `vh·bl` are exact and `num − vh·bh` is exact
    // (Sterbenz); the rest rounds at `2⁻²⁵` of the residual.
    let bc = SPLIT * den;
    let bh = bc - (bc - den);
    let bl = den - bh;
    let vh = f64::from_bits(v.to_bits() & !((1 << 27) - 1));
    let vl = v - vh;
    let rem = ((num - vh * bh) - vh * bl) - (vl * den + v * den_lo);

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
    let j = ((c * CENTRES as f64) as usize).min(CENTRES);
    // `q·π/4_hi` is exact; `|q·π/4| ≥ atan(c)` and `|q·π/4 ± atan(c)| ≥ |v|` unless zero.
    let (h0, l0) = fast_two_sum(q * fit.pi4_hi, sign * fit.atan_hi[j]);
    let us = sign * v;
    let (h1, l1) = fast_two_sum(h0, us);
    let s = v * v;
    let k = &fit.coeffs;
    let s2 = s * s;
    let p = (k[0] + k[1] * s) + s2 * ((k[2] + k[3] * s) + s2 * k[4]);
    // `atan(v + dv) = atan(v) + dv/(1 + s)`, `1/(1 + s)` to `s² ≤ 2⁻¹⁶`.
    let ok = if RESIDUAL <= num.abs() { sign } else { 0.0 };
    let early = (l1 + l0) + (q * fit.pi4_lo + sign * fit.atan_lo[j]);
    let mag = h1 + ((early + rem * (ok * (inv * (1.0 - s)))) + (us * s) * p);
    if x.is_nan() || y.is_nan() {
        x + y
    } else {
        mag.copysign(y)
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

/// `atan2(y, x)` by the binary32 kernel, evaluated in binary64.
pub(crate) fn atan2_f32(y: f32, x: f32, fit: &Fit32, pi4: (f64, f64)) -> f32 {
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
    let c = &fit.coeffs;
    let (s2, s4) = (s * s, (s * s) * (s * s));
    let p = ((c[0] + c[1] * s) + s2 * (c[2] + c[3] * s))
        + s4 * ((c[4] + c[5] * s) + s2 * (c[6] + c[7] * s));
    let mag = k * pi4.0 + (us + ((us * s) * p + k * pi4.1));
    if x.is_nan() || y.is_nan() {
        x + y
    } else {
        (mag as f32).copysign(y)
    }
}

/// The shipped scalar `T` with `atan2` alone replaced by the kernel at `T`'s precision.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub(crate) struct Portable<T>(pub(crate) T);

fn kernel_f64(y: f64, x: f64, fits: &Fits) -> f64 {
    atan2_f64(y, x, &fits.f64)
}

fn kernel_f32(y: f32, x: f32, fits: &Fits) -> f32 {
    atan2_f32(y, x, &fits.f32, (fits.f64.pi4_hi, fits.f64.pi4_lo))
}

macro_rules! portable {
    ($t:ty, $kernel:ident) => {
        impl core::ops::Add for Portable<$t> {
            type Output = Self;
            fn add(self, o: Self) -> Self {
                Portable(self.0 + o.0)
            }
        }
        impl core::ops::Sub for Portable<$t> {
            type Output = Self;
            fn sub(self, o: Self) -> Self {
                Portable(self.0 - o.0)
            }
        }
        impl core::ops::Mul for Portable<$t> {
            type Output = Self;
            fn mul(self, o: Self) -> Self {
                Portable(self.0 * o.0)
            }
        }
        impl core::ops::Div for Portable<$t> {
            type Output = Self;
            fn div(self, o: Self) -> Self {
                Portable(self.0 / o.0)
            }
        }
        impl core::ops::Neg for Portable<$t> {
            type Output = Self;
            fn neg(self) -> Self {
                Portable(-self.0)
            }
        }
        impl From<Portable<$t>> for f64 {
            fn from(x: Portable<$t>) -> f64 {
                f64::from(x.0)
            }
        }

        // Every method but `atan2` is the shipped scalar's, so the twin changes one variable.
        impl Real for Portable<$t> {
            type Mask = bool;
            const PRECISION: Precision = <$t as Real>::PRECISION;
            fn lit(x: f64) -> Self {
                Portable(<$t as Real>::lit(x))
            }
            fn round_from_f64(x: f64) -> Self {
                Portable(<$t as Real>::round_from_f64(x))
            }
            fn zero() -> Self {
                Portable(<$t as Real>::zero())
            }
            fn one() -> Self {
                Portable(<$t as Real>::one())
            }
            fn lt(self, rhs: Self) -> bool {
                <$t as Real>::lt(self.0, rhs.0)
            }
            fn le(self, rhs: Self) -> bool {
                <$t as Real>::le(self.0, rhs.0)
            }
            fn select(m: bool, t: Self, f: Self) -> Self {
                Portable(<$t as Real>::select(m, t.0, f.0))
            }
            fn sqrt(self) -> Self {
                Portable(Real::sqrt(self.0))
            }
            fn cbrt(self) -> Self {
                Portable(Real::cbrt(self.0))
            }
            fn sin_cos(self) -> (Self, Self) {
                let (s, c) = Real::sin_cos(self.0);
                (Portable(s), Portable(c))
            }
            fn sin(self) -> Self {
                Portable(Real::sin(self.0))
            }
            fn cos(self) -> Self {
                Portable(Real::cos(self.0))
            }
            fn acos(self) -> Self {
                Portable(Real::acos(self.0))
            }
            fn exp(self) -> Self {
                Portable(Real::exp(self.0))
            }
            fn exp_m1(self) -> Self {
                Portable(Real::exp_m1(self.0))
            }
            /// The candidate. An unparsed fit answers NaN; the twin answers nothing first.
            fn atan2(self, x: Self) -> Self {
                Portable(match fits() {
                    Ok(fits) => $kernel(self.0, x.0, fits),
                    Err(_) => <$t>::NAN,
                })
            }
            fn abs(self) -> Self {
                Portable(Real::abs(self.0))
            }
            fn copysign(self, sign: Self) -> Self {
                Portable(Real::copysign(self.0, sign.0))
            }
            fn value_f64(self) -> f64 {
                Real::value_f64(self.0)
            }
        }
    };
}

portable!(f64, kernel_f64);
portable!(f32, kernel_f32);

impl From<f64> for Portable<f64> {
    fn from(x: f64) -> Self {
        Portable(x)
    }
}
