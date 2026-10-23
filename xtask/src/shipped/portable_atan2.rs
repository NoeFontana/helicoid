//! `helicoid:portable-atan2`: the shipped program at [`Portable`], whose `atan2` is the portable
//! candidate of `0105` (draft) and whose every other operation is the shipped scalar's.
//!
//! The kernel is basic operations only, in one fixed sequence per precision:
//! - `n = min`, `d = max` of `|x|, |y|`, rescaled by an exact power of two where `d` leaves the
//!   window the compensations need;
//! - the fold `n − d` over `n + d` where `d/2 < n`, so `|u| ≤ ½` and the numerator is exact
//!   (Sterbenz), the denominator's rounding kept by `fast_two_sum`;
//! - the quotient's residual by Dekker's product, with no fused multiply-add (D16);
//! - `|atan2| = k·π/4 + atan(u)` with `k·π/4_hi` exact for `k ≤ 4`, `atan(u) = u + u s P(s)`.
//!
//! `P` and the split of `π/4` are `conformance/fits/atan2.json`, written by
//! `python -m gen atan-fit`; nothing here is a typed coefficient. The window and the splitter are
//! the formats' own: `2^⌈p/2⌉ + 1` splits a `p`-bit significand, and the window keeps
//! `c·(n + d)` finite and the residual's products normal.
//!
//! Planted, as `helicoid:host-std` is: a twin is not a candidate, and no bar reads its rows. It
//! answers both precisions, since every operation but `atan2` is the shipped scalar's own.

use std::sync::OnceLock;

use helicoid_linalg::{Mask, Precision, Real};
use serde::Deserialize;

use super::{program, Helicoid};
use crate::conformance::corpus::Record;
use crate::conformance::subject::{Output, Registered, Subject};

/// The subject's name, and the value of `--subject`.
pub(crate) const NAME: &str = "helicoid:portable-atan2";

/// The committed fit, as the generator wrote it.
const FIT: &str = include_str!("../../../conformance/fits/atan2.json");

/// One precision's polynomial and `π/4` split.
#[derive(Debug)]
pub(crate) struct Fit<S> {
    /// `P`'s coefficients, constant term first.
    pub(crate) coeffs: Vec<S>,
    pub(crate) pi4_hi: S,
    pub(crate) pi4_lo: S,
}

#[derive(Deserialize)]
struct FitFile {
    f64: FitJson,
    f32: FitJson,
}

#[derive(Deserialize)]
struct FitJson {
    umax: String,
    terms: usize,
    coeffs: Vec<String>,
    pi4_hi: String,
    pi4_lo: String,
}

fn bits(hex: &str) -> Result<f64, String> {
    let digits = hex.strip_prefix("0x").ok_or(format!("{hex}: no 0x"))?;
    u64::from_str_radix(digits, 16)
        .map(f64::from_bits)
        .map_err(|e| format!("{hex}: {e}"))
}

impl FitJson {
    /// The fit at `S`, every value exactly representable there.
    fn at<S: Real>(&self) -> Result<Fit<S>, String> {
        // The kernel's fold is at `d/2 < n`, so `|u| ≤ ½`; a fit on another interval is a
        // different kernel.
        if self.umax != "1/2" || self.coeffs.len() != self.terms {
            return Err(format!(
                "fit on |u| <= {} with {} of {} terms; the kernel folds to 1/2",
                self.umax,
                self.coeffs.len(),
                self.terms
            ));
        }
        let exact = |hex: &str| -> Result<S, String> {
            let x = bits(hex)?;
            let at = S::round_from_f64(x);
            (at.value_f64().to_bits() == x.to_bits())
                .then_some(at)
                .ok_or(format!("{hex} is not exact at {:?}", S::PRECISION))
        };
        Ok(Fit {
            coeffs: self
                .coeffs
                .iter()
                .map(|c| exact(c))
                .collect::<Result<_, _>>()?,
            pi4_hi: exact(&self.pi4_hi)?,
            pi4_lo: exact(&self.pi4_lo)?,
        })
    }
}

/// The binary64 fit, then the binary32 one.
pub(crate) type Fits = (Fit<f64>, Fit<f32>);

/// The two fits, parsed once.
pub(crate) fn fits() -> Result<&'static Fits, &'static str> {
    static FITS: OnceLock<Result<Fits, String>> = OnceLock::new();
    FITS.get_or_init(|| {
        let file: FitFile = serde_json::from_str(FIT).map_err(|e| format!("atan2.json: {e}"))?;
        Ok((file.f64.at()?, file.f32.at()?))
    })
    .as_ref()
    .map_err(String::as_str)
}

/// The format's constants: the window on `d`, the powers of two that bring `d` back into it, the
/// splitter, and the numerator below which the residual's products are not all normal.
struct Format<S> {
    lo: S,
    hi: S,
    up: S,
    down: S,
    split: S,
    residual: S,
}

impl<S: Real> Format<S> {
    fn of() -> Self {
        let pow2 = |e: i32| S::lit(libm::ldexp(1.0, e));
        match S::PRECISION {
            // `(2^27 + 1)·2d` stays below `2^1024` for `d < 2^995`; `num·2^-53` is normal for
            // `|num| ≥ 2^-969`.
            Precision::F64 => Self {
                lo: pow2(-460),
                hi: pow2(990),
                up: pow2(512),
                down: pow2(-512),
                split: pow2(27) + S::one(),
                residual: pow2(-969),
            },
            // `(2^12 + 1)·2d` stays below `2^128` for `d < 2^114`; `num·2^-24` is normal for
            // `|num| ≥ 2^-102`.
            Precision::F32 => Self {
                lo: pow2(-50),
                hi: pow2(110),
                up: pow2(64),
                down: pow2(-64),
                split: pow2(12) + S::one(),
                residual: pow2(-102),
            },
        }
    }
}

/// `(s, e)` with `s + e = a + b` exactly, for `|a| ≥ |b|` or `a = 0`.
fn fast_two_sum<S: Real>(a: S, b: S) -> (S, S) {
    let s = a + b;
    (s, b - (s - a))
}

/// `(h, l)` with `h + l = a·b` exactly where no product underflows (Dekker), without `mul_add`.
fn two_product<S: Real>(a: S, b: S, split: S) -> (S, S) {
    let halves = |x: S| {
        let c = split * x;
        let hi = c - (c - x);
        (hi, x - hi)
    };
    let h = a * b;
    let ((ah, al), (bh, bl)) = (halves(a), halves(b));
    (h, ((ah * bh - h) + ah * bl + al * bh) + al * bl)
}

/// `(n, d)` with `d` brought into the window by an exact power of two, or, where `d` is infinite,
/// IEEE 754's limits: each infinite argument `1`, each finite one `0`. A NaN passes through.
#[cold]
#[inline(never)]
fn rescale<S: Real>(n: S, d: S, f: &Format<S>) -> (S, S) {
    let inf = S::lit(f64::INFINITY);
    let (n_inf, d_inf) = (inf.le(n), inf.le(d));
    let big = f.hi.lt(d);
    let scale = S::select(big, f.down, f.up);
    let n = S::select(d_inf, S::select(n_inf, S::one(), S::zero()), n * scale);
    let d = S::select(d_inf, S::one(), d * scale);
    (n, d)
}

/// `atan2(y, x)` by the portable kernel of `0105` (draft), at `fit`'s precision.
pub(crate) fn atan2<S: Real>(y: S, x: S, fit: &Fit<S>) -> S {
    let (zero, one) = (S::zero(), S::one());
    let f = Format::<S>::of();
    let (ax, ay) = (x.abs(), y.abs());
    let swap = ax.lt(ay);
    let (n, d) = (S::select(swap, ax, ay), S::select(swap, ay, ax));
    let inside = f.lo.le(d).and(d.le(f.hi));
    let fr = &f;
    let (n, d) = S::branch(inside, move || (n, d), move || rescale(n, d, fr));

    // `|u| ≤ ½`: the fold's numerator is exact by Sterbenz (`d/2 < n ≤ d`), its denominator's
    // rounding is `den_lo`.
    let fold = (S::lit(0.5) * d).lt(n);
    let (sum, sum_lo) = fast_two_sum(d, n);
    let num = S::select(fold, n - d, n);
    let empty = d.le(zero);
    let den = S::select(empty, one, S::select(fold, sum, d));
    let den_lo = S::select(fold, sum_lo, zero);
    let u = num / den;
    // `num − ph` is exact (Sterbenz: `ph` is `num` to a rounding).
    let (ph, pl) = two_product(u, den, f.split);
    let rem = ((num - ph) - pl) - u * den_lo;
    let normal = f.residual.le(num.abs());
    let du = S::select(normal.and(empty.not()), rem / den, zero);

    // The octant: `k·π/4 ± atan(u)`, the sign flipped where exactly one of `x < 0`, `|y| > |x|`.
    let neg_x = one.copysign(x).lt(zero);
    let flip = neg_x.and(swap.not()).or(swap.and(neg_x.not()));
    let sign = S::select(flip, -one, one);
    let q = S::select(swap, S::lit(2.0), S::select(neg_x, S::lit(4.0), zero));
    let k = q + S::select(fold, sign, zero);
    let (us, dus) = (sign * u, sign * du);

    let s = us * us;
    let mut p = zero;
    for &c in fit.coeffs.iter().rev() {
        p = p * s + c;
    }
    // `atan(u + du) = atan(u) + du/(1 + s)`, `1/(1 + s)` to first order.
    let tail = (us * s) * p + dus * (one - s);
    let (hi, lo) = fast_two_sum(k * fit.pi4_hi, us);
    let mag = hi + (lo + (k * fit.pi4_lo + tail));
    let nan = x.le(x).and(y.le(y)).not();
    S::select(nan, x + y, mag.copysign(y))
}

/// The shipped scalar `T` with `atan2` alone replaced by [`atan2`].
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub(crate) struct Portable<T>(pub(crate) T);

macro_rules! portable {
    ($t:ty, $fit:tt) => {
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
            /// The candidate. An unparsed fit answers NaN; [`Twin::eval`] answers nothing first.
            fn atan2(self, x: Self) -> Self {
                Portable(match fits() {
                    Ok(fits) => atan2(self.0, x.0, &fits.$fit),
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

portable!(f64, 0);
portable!(f32, 1);

impl From<f64> for Portable<f64> {
    fn from(x: f64) -> Self {
        Portable(x)
    }
}

/// The shipped program at [`Portable`].
pub(crate) struct Twin;

impl Subject for Twin {
    fn name(&self) -> &str {
        NAME
    }

    fn supports(&self, fn_id: &str) -> bool {
        Helicoid.supports(fn_id)
    }

    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        if fits().is_err() {
            return Output::new();
        }
        program::<Portable<f64>, Portable<f32>>(fn_id, record, precision)
    }
}

pub(crate) fn registered() -> Registered {
    let version = format!("{}@portable-atan2", env!("CARGO_PKG_VERSION"));
    Registered {
        version: version.clone(),
        version_f32: version,
        no_f32: None,
        planted: true,
        subject: Box::new(Twin),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// splitmix64, so a draw is the same on every host.
    struct Draw(u64);

    impl Draw {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }

        /// `±m·2^e`, `m` uniform in `[1, 2)`, `e` uniform in `[lo, hi]`.
        fn wide(&mut self, lo: i32, hi: i32) -> f64 {
            let r = self.next();
            let span = u64::try_from(hi - lo + 1).unwrap_or(1);
            let e = lo + i32::try_from(r % span).unwrap_or(0);
            let m = 1.0 + (r >> 11) as f64 / (1u64 << 53) as f64;
            let x = libm::ldexp(m, e);
            if r & 1 == 0 {
                x
            } else {
                -x
            }
        }
    }

    /// A double-double reference, about `2^-100` relative: enough to read an error in binary64
    /// ulps to many digits. `π/4` is its own `atan(1)`, so nothing here is a typed constant.
    mod dd {
        use super::super::{fast_two_sum, two_product};

        #[derive(Clone, Copy, Debug)]
        pub(super) struct Dd(pub(super) f64, pub(super) f64);

        const SPLIT: f64 = 134_217_729.0;

        fn two_sum(a: f64, b: f64) -> (f64, f64) {
            let s = a + b;
            let bb = s - a;
            (s, (a - (s - bb)) + (b - bb))
        }

        impl Dd {
            pub(super) fn of(x: f64) -> Self {
                Dd(x, 0.0)
            }
            pub(super) fn add(self, o: Self) -> Self {
                let (s, e) = two_sum(self.0, o.0);
                let (t, f) = two_sum(self.1, o.1);
                let (s, e) = fast_two_sum(s, e + t);
                let (s, e) = fast_two_sum(s, e + f);
                Dd(s, e)
            }
            pub(super) fn neg(self) -> Self {
                Dd(-self.0, -self.1)
            }
            pub(super) fn mul(self, o: Self) -> Self {
                let (p, e) = two_product(self.0, o.0, SPLIT);
                let (p, e) = fast_two_sum(p, e + (self.0 * o.1 + self.1 * o.0));
                Dd(p, e)
            }
            pub(super) fn div(self, o: Self) -> Self {
                let q1 = self.0 / o.0;
                let r = self.add(o.mul(Dd::of(q1)).neg());
                let q2 = r.0 / o.0;
                let r = r.add(o.mul(Dd::of(q2)).neg());
                let q3 = r.0 / o.0;
                Dd::of(q1).add(Dd::of(q2)).add(Dd::of(q3))
            }
            pub(super) fn sqrt(self) -> Self {
                let s = self.0.sqrt();
                let r = self.add(Dd::of(s).mul(Dd::of(s)).neg());
                Dd::of(s).add(Dd::of(r.0 / (2.0 * s)))
            }
        }

        /// `atan(t)` for `0 ≤ t ≤ 1`: four half-angle steps, then the series.
        pub(super) fn atan(t: Dd) -> Dd {
            let one = Dd::of(1.0);
            let mut z = t;
            for _ in 0..4 {
                z = z.div(one.add(one.add(z.mul(z)).sqrt()));
            }
            let z2 = z.mul(z);
            let (mut sum, mut pow) = (Dd::of(0.0), z);
            for j in 0..20 {
                let term = pow.div(Dd::of(f64::from(2 * j + 1)));
                sum = if j % 2 == 0 {
                    sum.add(term)
                } else {
                    sum.add(term.neg())
                };
                pow = pow.mul(z2);
            }
            sum.mul(Dd::of(16.0))
        }

        /// `atan2(y, x)` for finite, non-zero arguments.
        pub(super) fn atan2(y: f64, x: f64) -> Dd {
            let (ax, ay) = (x.abs(), y.abs());
            let (n, d) = if ay > ax { (ax, ay) } else { (ay, ax) };
            let a = atan(Dd::of(n).div(Dd::of(d)));
            let pi4 = atan(Dd::of(1.0));
            let pi2 = pi4.mul(Dd::of(2.0));
            let pi = pi4.mul(Dd::of(4.0));
            let a = match (x < 0.0, ay > ax) {
                (false, false) => a,
                (false, true) => pi2.add(a.neg()),
                (true, true) => pi2.add(a),
                (true, false) => pi.add(a.neg()),
            };
            if y < 0.0 {
                a.neg()
            } else {
                a
            }
        }
    }

    /// `|got − want|` in ulps of `want` at `bits` significand bits.
    fn ulps(got: f64, want: dd::Dd, bits: i32) -> f64 {
        let e = libm::ilogb(want.0).max(if bits == 53 { -1022 } else { -126 });
        ((got - want.0) - want.1).abs() / libm::ldexp(1.0, e - (bits - 1))
    }

    /// The worst error of `subject` at binary64 over `draws` pairs, and where.
    fn f64_error(
        draws: usize,
        seed: u64,
        (lo, hi): (i32, i32),
        subject: impl Fn(f64, f64) -> Result<f64, String>,
    ) -> Result<(f64, f64, f64), String> {
        let mut draw = Draw(seed);
        let mut worst = (0.0, 0.0, 0.0);
        for _ in 0..draws {
            let (y, x) = (draw.wide(lo, hi), draw.wide(lo, hi));
            let err = ulps(subject(y, x)?, dd::atan2(y, x), 53);
            if err > worst.0 {
                worst = (err, y, x);
            }
        }
        Ok(worst)
    }

    /// The tuning readout: the candidate's, `libm`'s and the host's worst binary64 error against
    /// the double-double reference, near and far ratios.
    ///
    /// `cargo nextest run -p xtask --run-ignored only --no-capture -- portable_atan2`.
    #[test]
    #[ignore = "a measurement, printed for 0105 (draft) to cite"]
    #[allow(clippy::print_stdout)]
    fn the_binary64_kernel_reads_its_error() -> Result<(), String> {
        for (range, label) in [((-2, 2), "|y/x| in 2^±4"), ((-30, 30), "|y/x| in 2^±60")] {
            let (err, y, x) = f64_error(1 << 20, 5, range, at64)?;
            println!("f64 {label}: max {err:.4} ulp at atan2({y:e}, {x:e})");
            let (err, y, x) = f64_error(1 << 20, 5, range, |y, x| Ok(libm::atan2(y, x)))?;
            println!("libm::atan2 {label}: max {err:.4} ulp at atan2({y:e}, {x:e})");
            let (err, y, x) = f64_error(1 << 20, 5, range, |y, x| Ok(y.atan2(x)))?;
            println!("host std {label}: max {err:.4} ulp at atan2({y:e}, {x:e})");
        }
        Ok(())
    }

    /// The candidate's worst binary64 error per arm: `k` and whether the quotient folded.
    #[test]
    #[ignore = "a measurement, printed for 0105 (draft) to cite"]
    #[allow(clippy::print_stdout)]
    fn the_binary64_kernel_reads_its_error_per_arm() -> Result<(), String> {
        let mut draw = Draw(6);
        let mut worst = std::collections::BTreeMap::<(u8, bool), (f64, f64)>::new();
        for _ in 0..1 << 20 {
            let (y, x) = (draw.wide(-2, 2), draw.wide(-2, 2));
            let (n, d) = (y.abs().min(x.abs()), y.abs().max(x.abs()));
            let (fold, swap, neg_x) = (0.5 * d < n, y.abs() > x.abs(), x < 0.0);
            // The kernel's `k = q ± fold`, the sign `−` where exactly one of `neg_x`, `swap`.
            let q: u8 = if swap {
                2
            } else if neg_x {
                4
            } else {
                0
            };
            let k = match (fold, neg_x != swap) {
                (false, _) => q,
                (true, false) => q + 1,
                (true, true) => q - 1,
            };
            let err = ulps(at64(y, x)?, dd::atan2(y, x), 53);
            let cell = worst.entry((k, fold)).or_insert((0.0, 0.0));
            if err > cell.0 {
                *cell = (err, n / d);
            }
        }
        for ((k, fold), (err, t)) in &worst {
            println!("k {k} fold {fold:5}: max {err:.4} ulp at n/d = {t:.4}");
        }
        Ok(())
    }

    /// The reference against itself where the answer is known: `atan2(1, 1) = π/4` and
    /// `atan2(1, 0+) ≈ π/2` both round to `core`'s constants.
    #[test]
    fn the_reference_rounds_to_the_known_constants() {
        let pi4 = dd::atan2(1.0, 1.0);
        assert_eq!(pi4.0.to_bits(), core::f64::consts::FRAC_PI_4.to_bits());
        assert!(ulps(core::f64::consts::FRAC_PI_4, pi4, 53) < 0.5);
        let pi = dd::atan2(1e-300, -1.0);
        assert_eq!(pi.0.to_bits(), core::f64::consts::PI.to_bits());
    }

    fn at64(y: f64, x: f64) -> Result<f64, String> {
        Ok(atan2(y, x, &fits()?.0))
    }

    fn at32(y: f32, x: f32) -> Result<f32, String> {
        Ok(atan2(y, x, &fits()?.1))
    }

    fn same64(a: f64, b: f64) -> bool {
        a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan())
    }

    fn same32(a: f32, b: f32) -> bool {
        a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan())
    }

    #[test]
    fn the_committed_fit_parses() -> Result<(), String> {
        let (f64_fit, f32_fit) = fits()?;
        assert!(!f64_fit.coeffs.is_empty() && !f32_fit.coeffs.is_empty());
        Ok(())
    }

    /// `k·π/4_hi` is exact for `k ≤ 4`: `π/4_hi` has two spare low bits at each precision.
    #[test]
    fn the_high_part_of_pi_over_4_has_two_spare_bits() -> Result<(), String> {
        let (a, b) = fits()?;
        assert!(a.pi4_hi.to_bits().trailing_zeros() >= 2);
        assert!(b.pi4_hi.to_bits().trailing_zeros() >= 2);
        Ok(())
    }

    /// Wherever an argument is zero, infinite or NaN, the answer is IEEE 754's, correctly rounded:
    /// `±0`, `±π/2`, `±π`, `±π/4`, `±3π/4`, or NaN. At `f64` that is `libm`'s to the bit. At `f32`
    /// it is `libm`'s binary64 answer rounded, since `libm::atan2f` returns `π/2` one ulp low where
    /// `y` is infinite or dwarfs `x` (`0x3fc90fda` against the nearest, `0x3fc90fdb`); those
    /// binary64 constants are nowhere near a binary32 midpoint.
    #[test]
    fn special_values_are_ieee_s() -> Result<(), String> {
        let mut atan2f_off = Vec::new();
        let edge = [0.0, -0.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN];
        let any = [
            0.0,
            -0.0,
            1.0,
            -1.0,
            3.5e-310,
            -2.0e300,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
        ];
        for &a in &edge {
            for &b in &any {
                for (y, x) in [(a, b), (b, a)] {
                    let (got, want) = (at64(y, x)?, libm::atan2(y, x));
                    assert!(
                        same64(got, want),
                        "f64 atan2({y:e}, {x:e}) = {got:e}, libm {want:e}"
                    );
                    let (y, x) = (y as f32, x as f32);
                    let want = libm::atan2(f64::from(y), f64::from(x)) as f32;
                    let got = at32(y, x)?;
                    assert!(
                        same32(got, want),
                        "f32 atan2({y:e}, {x:e}) = {got:e}, want {want:e}"
                    );
                    if !same32(libm::atan2f(y, x), want) {
                        atan2f_off.push((y, x));
                    }
                }
            }
        }
        let inf = f32::INFINITY;
        assert_eq!(atan2f_off, [(inf, 1.0), (-inf, 1.0)]);
        Ok(())
    }

    /// `atan2(−y, x) = −atan2(y, x)` bit for bit: the kernel works on `|y|` and restores its sign.
    #[test]
    fn atan2_is_odd_in_y() -> Result<(), String> {
        let mut draw = Draw(1);
        for _ in 0..100_000 {
            let (y, x) = (draw.wide(-1074, 1023), draw.wide(-1074, 1023));
            assert!(same64(at64(-y, x)?, -at64(y, x)?), "({y:e}, {x:e})");
            let (y, x) = (y as f32, x as f32);
            assert!(same32(at32(-y, x)?, -at32(y, x)?), "({y:e}, {x:e})");
        }
        Ok(())
    }

    /// The rescale arm is exact: scaling both arguments by `2^k` changes no bit, through the window's
    /// edges and the cold arm on either side, while neither scaled argument is subnormal.
    #[test]
    fn scaling_both_arguments_by_a_power_of_two_changes_no_bit() -> Result<(), String> {
        let mut draw = Draw(2);
        for _ in 0..2_000 {
            let (y, x) = (draw.wide(-4, 4), draw.wide(-4, 4));
            let (want64, want32) = (at64(y, x)?, at32(y as f32, x as f32)?);
            for k in (-1015..=1015).step_by(7) {
                let s = libm::ldexp(1.0, k);
                assert!(
                    same64(at64(y * s, x * s)?, want64),
                    "({y:e}, {x:e}) at 2^{k}"
                );
            }
            for k in -118..=120 {
                let s = libm::ldexpf(1.0, k);
                let (y, x) = (y as f32, x as f32);
                assert!(
                    same32(at32(y * s, x * s)?, want32),
                    "({y:e}, {x:e}) at 2^{k}"
                );
            }
        }
        Ok(())
    }

    /// Monotone in `y` at `x = 1` through the fold (`y = ½`), the swap (`y = 1`) and the swapped
    /// fold (`y = 2`), a thousand ulps either side of each.
    #[test]
    fn monotone_through_the_fold_and_the_swap() -> Result<(), String> {
        for at in [0.5f64, 1.0, 2.0] {
            let mut y = at;
            for _ in 0..1000 {
                y = f64::from_bits(y.to_bits() - 1);
            }
            let mut last = at64(y, 1.0)?;
            for _ in 0..2000 {
                y = f64::from_bits(y.to_bits() + 1);
                let next = at64(y, 1.0)?;
                assert!(last <= next, "f64 at y = {y:e}: {last:e} then {next:e}");
                last = next;
            }
            let mut y = at as f32;
            for _ in 0..1000 {
                y = f32::from_bits(y.to_bits() - 1);
            }
            let mut last = at32(y, 1.0)?;
            for _ in 0..2000 {
                y = f32::from_bits(y.to_bits() + 1);
                let next = at32(y, 1.0)?;
                assert!(last <= next, "f32 at y = {y:e}: {last:e} then {next:e}");
                last = next;
            }
        }
        Ok(())
    }

    /// The binary32 kernel's error in its own ulps, against `libm`'s binary64 `atan2` (whose error is
    /// `2^-29` of a binary32 ulp), over every octant, both folds and arguments `2^±60` apart.
    fn f32_error(
        draws: usize,
        seed: u64,
        subject: impl Fn(f32, f32) -> Result<f32, String>,
    ) -> Result<(f64, f32, f32), String> {
        let mut draw = Draw(seed);
        let mut worst = (0.0, 0.0f32, 0.0f32);
        for i in 0..draws {
            let (y, x) = if i % 2 == 0 {
                (draw.wide(-60, 60), draw.wide(-60, 60))
            } else {
                (draw.wide(-2, 2), draw.wide(-2, 2))
            };
            let (y, x) = (y as f32, x as f32);
            let want = libm::atan2(f64::from(y), f64::from(x));
            let got = f64::from(subject(y, x)?);
            let ulp = libm::ldexp(1.0, libm::ilogb(want).max(-126) - 23);
            let err = (got - want).abs() / ulp;
            if err > worst.0 {
                worst = (err, y, x);
            }
        }
        Ok(worst)
    }

    #[test]
    #[ignore = "a measurement, printed for 0105 (draft) to cite"]
    #[allow(clippy::print_stdout)]
    fn the_binary32_kernel_reads_its_error() -> Result<(), String> {
        let (err, y, x) = f32_error(1 << 20, 3, at32)?;
        println!("f32 max {err:.4} ulp at atan2({y:e}, {x:e})");
        let (err, y, x) = f32_error(1 << 20, 3, |y, x| Ok(libm::atan2f(y, x)))?;
        println!("libm::atan2f max {err:.4} ulp at atan2({y:e}, {x:e})");
        Ok(())
    }

    /// Binary64 against `libm`: a smoke bound, since `libm` is itself up to an ulp off; the measured
    /// readings are the corpus's (`helicoid:portable-atan2` over `real_atan2`).
    #[test]
    fn binary64_is_within_two_ulps_of_libm() -> Result<(), String> {
        let mut draw = Draw(4);
        let mut worst = 0u64;
        for _ in 0..1_000_000 {
            let (y, x) = (draw.wide(-30, 30), draw.wide(-30, 30));
            let (got, want) = (at64(y, x)?, libm::atan2(y, x));
            worst = worst.max(got.to_bits().abs_diff(want.to_bits()));
        }
        assert!(worst <= 2, "{worst} ulps");
        Ok(())
    }

    /// The guard every run keeps: within an ulp of the reference at both precisions, on a sample
    /// the ignored readouts above take sixteen times over.
    #[test]
    fn the_kernel_is_within_an_ulp() -> Result<(), String> {
        let (err, y, x) = f64_error(1 << 16, 7, (-30, 30), at64)?;
        assert!(err < 1.0, "f64 {err} ulp at atan2({y:e}, {x:e})");
        let (err, y, x) = f32_error(1 << 16, 7, at32)?;
        assert!(err < 1.0, "f32 {err} ulp at atan2({y:e}, {x:e})");
        Ok(())
    }
}
