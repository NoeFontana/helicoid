//! The corpus's number formats, exactly: hex-float inputs (`float.hex()`), decimal references, and
//! the correctly rounded conversion of an integer ratio, and of its square root, to binary64.
//! Nothing here goes through `str::parse::<f64>` or a float rounding, so an error formed from
//! these is exact (`docs/maths/error-analysis.md` EA.19(b)).

use std::str::FromStr;

use num_bigint::BigUint;

/// `±mant · 2^exp`, the exact value of a finite binary64.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Dyadic {
    pub(crate) neg: bool,
    pub(crate) mant: u64,
    pub(crate) exp: i32,
}

impl Dyadic {
    /// `x` must be finite: the harness counts non-finite outputs before it forms an error.
    pub(crate) fn of(x: f64) -> Self {
        let bits = x.to_bits();
        let biased = ((bits >> 52) & 0x7ff) as i32;
        let frac = bits & ((1 << 52) - 1);
        let (mant, exp) = if biased == 0 {
            (frac, -1074)
        } else {
            (frac | (1 << 52), biased - 1075)
        };
        Self {
            neg: bits >> 63 == 1,
            mant,
            exp,
        }
    }
}

/// `±mant · 10^exp10`, the exact value of a reference string. A zero has `mant == 0`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Decimal {
    pub(crate) neg: bool,
    pub(crate) mant: u128,
    pub(crate) exp10: i32,
}

impl FromStr for Decimal {
    type Err = String;

    /// `[+-]digits[.digits][e[+-]digits]`; at most 38 digits (30 are written).
    fn from_str(s: &str) -> Result<Self, String> {
        let bad = || format!("bad decimal `{s}`");
        let (neg, rest) = split_sign(s);
        let (num, exp) = match rest.split_once(['e', 'E']) {
            Some((n, e)) => (n, e.parse::<i32>().map_err(|_| bad())?),
            None => (rest, 0),
        };
        let (int, frac) = num.split_once('.').unwrap_or((num, ""));
        if int.is_empty() && frac.is_empty() {
            return Err(bad());
        }
        let mut mant = 0u128;
        for c in int.chars().chain(frac.chars()) {
            let d = c.to_digit(10).ok_or_else(bad)?;
            mant = mant
                .checked_mul(10)
                .and_then(|m| m.checked_add(u128::from(d)))
                .ok_or_else(|| format!("decimal `{s}` has more than 38 digits"))?;
        }
        let exp10 = i32::try_from(frac.len())
            .ok()
            .and_then(|f| exp.checked_sub(f))
            .ok_or_else(bad)?;
        Ok(Self { neg, mant, exp10 })
    }
}

fn split_sign(s: &str) -> (bool, &str) {
    match s.strip_prefix('-') {
        Some(r) => (true, r),
        None => (false, s.strip_prefix('+').unwrap_or(s)),
    }
}

/// Parses `[-]0xH.HHHp[+-]D` (Python's `float.hex()`) to the binary64 it denotes. A value that is
/// not exactly a binary64 is an error: an input the subject cannot receive is a corpus defect.
pub(crate) fn parse_hex(s: &str) -> Result<f64, String> {
    let bad = || format!("bad hex float `{s}`");
    let (neg, rest) = split_sign(s);
    let (mant_txt, exp_txt) = rest
        .strip_prefix("0x")
        .and_then(|r| r.split_once('p'))
        .ok_or_else(bad)?;
    let (int, frac) = mant_txt.split_once('.').unwrap_or((mant_txt, ""));
    if int.is_empty() && frac.is_empty() {
        return Err(bad());
    }
    let mut mant = 0u128;
    for c in int.chars().chain(frac.chars()) {
        let d = c.to_digit(16).ok_or_else(bad)?;
        mant = mant
            .checked_mul(16)
            .and_then(|m| m.checked_add(u128::from(d)))
            .ok_or_else(bad)?;
    }
    let exp = i32::try_from(frac.len())
        .ok()
        .and_then(|f| exp_txt.parse::<i32>().ok()?.checked_sub(4 * f))
        .ok_or_else(bad)?;
    binary64(neg, mant, exp).ok_or_else(|| format!("hex float `{s}` is not a binary64"))
}

/// `±mant · 2^exp` as a binary64 when it is exactly one (subnormals included).
fn binary64(neg: bool, mut mant: u128, mut exp: i32) -> Option<f64> {
    let sign = u64::from(neg) << 63;
    if mant == 0 {
        return Some(f64::from_bits(sign));
    }
    let tz = mant.trailing_zeros();
    mant >>= tz;
    exp = exp.checked_add(i32::try_from(tz).ok()?)?;
    let width = i32::try_from(128 - mant.leading_zeros()).ok()?;
    let top = exp.checked_add(width - 1)?;
    if width > 53 || top > 1023 {
        return None;
    }
    let mant = u64::try_from(mant).ok()?;
    if top >= -1022 {
        let frac = (mant << (53 - width)) & ((1 << 52) - 1);
        Some(f64::from_bits(
            sign | (u64::try_from(top + 1023).ok()? << 52) | frac,
        ))
    } else if exp >= -1074 {
        Some(f64::from_bits(sign | (mant << (exp + 1074))))
    } else {
        None
    }
}

/// `2^k` for `-1074 <= k <= 1023`.
fn pow2(k: i64) -> f64 {
    if k >= -1022 {
        f64::from_bits(((k + 1023) as u64) << 52)
    } else {
        f64::from_bits(1 << (k + 1074))
    }
}

/// `num / den` correctly rounded to binary64 (nearest, ties to even, subnormals and overflow
/// included). `den` must be non-zero. The quotient is cut at the last bit the result keeps, so
/// there is one rounding and no double rounding through a wider format. It rounds the seeded
/// kernels' series constants and is the test oracles' reference rounding.
pub(crate) fn ratio_to_f64(num: &BigUint, den: &BigUint) -> f64 {
    ratio_to_binary(num, den, Format::F64)
}

/// [`ratio_to_f64`] for binary32: the same single rounding, cut at the last bit binary32 keeps.
pub(crate) fn ratio_to_f32(num: &BigUint, den: &BigUint) -> f32 {
    // A binary32 value or infinity, exact in binary64, so the cast rounds nothing.
    ratio_to_binary(num, den, Format::F32) as f32
}

/// A binary format by its fraction bits and its normal exponent range.
#[derive(Clone, Copy)]
struct Format {
    fraction: i64,
    emin: i64,
    emax: i64,
}

impl Format {
    const F64: Self = Self {
        fraction: 52,
        emin: -1022,
        emax: 1023,
    };
    const F32: Self = Self {
        fraction: 23,
        emin: -126,
        emax: 127,
    };
}

/// `num / den` rounded once to `format`, returned as a binary64 (exact for either format).
fn ratio_to_binary(num: &BigUint, den: &BigUint, format: Format) -> f64 {
    if num.bits() == 0 {
        return 0.0;
    }
    let mut e = num.bits() as i64 - den.bits() as i64;
    let below = if e >= 0 {
        *num < (den << e as usize)
    } else {
        (num << (-e) as usize) < *den
    };
    e -= i64::from(below);
    if e > format.emax {
        return f64::INFINITY;
    }
    let quantum = e.max(format.emin) - format.fraction;
    let (n, d) = if quantum >= 0 {
        (num.clone(), den << quantum as usize)
    } else {
        (num << (-quantum) as usize, den.clone())
    };
    let mut q = &n / &d;
    let twice = (&n - &q * &d) << 1usize;
    if twice > d || (twice == d && q.bit(0)) {
        q += 1u32;
    }
    // `q <= 2^(fraction + 1)`: exact as an `f64`, and `q * 2^quantum` is representable by
    // construction or overflows to infinity.
    let q = q.iter_u64_digits().next().unwrap_or(0);
    q as f64 * pow2(quantum)
}

/// `sqrt(num / den)` correctly rounded to binary64 (nearest, ties to even), for ratios far outside
/// the `f64` range: an error of `1e-300` in units of `u` squares to nothing. `den` must be non-zero.
/// The ratio is scaled by `2^(2k)` until its integer root has at least 64 bits, the root is taken
/// exactly with a sticky bit for whatever the division and the root cut, and the result is rounded
/// once, at the last bit it keeps: nothing passes through a wider or narrower float first.
pub(crate) fn sqrt_ratio(num: &BigUint, den: &BigUint) -> f64 {
    if num.bits() == 0 {
        return 0.0;
    }
    // num/den is in (2^(e-1), 2^(e+1)): the root is at least 2^1024 above 2049, below 2^-1075
    // under -2154.
    let e = num.bits() as i64 - den.bits() as i64;
    if e > 2049 {
        return f64::INFINITY;
    }
    if e < -2154 {
        return 0.0;
    }
    // 2k >= 127 - e: num * 2^(2k) / den >= 2^126, so the root has at least 64 bits.
    let k = (128 - e).div_euclid(2);
    let (n, d) = if k >= 0 {
        (num << (2 * k) as usize, den.clone())
    } else {
        (num.clone(), den << (-2 * k) as usize)
    };
    let scaled = &n / &d;
    let root = scaled.sqrt();
    // The true root, times 2^k, is in [root, root + 1), and equal to `root` only when neither the
    // division nor the root left a remainder.
    let sticky = &scaled * &d != n || &root * &root != scaled;
    let top = root.bits() as i64 - 1 - k;
    if top > 1023 {
        return f64::INFINITY;
    }
    // The result keeps the bits of `root` from `cut` up; `cut >= 11` in the normal range, and
    // larger below it, so the half below the cut is never empty.
    let quantum = top.max(-1022) - 52;
    let cut = usize::try_from(k + quantum).unwrap_or(1);
    let mut kept = &root >> cut;
    let rest = &root - (&kept << cut);
    let half = BigUint::from(1u32) << (cut - 1);
    if rest > half || (rest == half && (sticky || kept.bit(0))) {
        kept += 1u32;
    }
    // `kept <= 2^53`: exact as an `f64`, and `kept * 2^quantum` is exact or overflows to infinity.
    kept.iter_u64_digits().next().unwrap_or(0) as f64 * pow2(quantum)
}

#[cfg(test)]
impl Decimal {
    /// The reference as an exact integer ratio, for the roundings below.
    fn ratio(self) -> (BigUint, BigUint) {
        let mant = BigUint::from(self.mant);
        let ten = BigUint::from(10u32);
        let exp = self.exp10.unsigned_abs();
        if self.exp10 >= 0 {
            (mant * ten.pow(exp), BigUint::from(1u32))
        } else {
            (mant, ten.pow(exp))
        }
    }

    /// The correctly rounded binary64 of the reference: what a perfect subject returns.
    pub(crate) fn to_f64(self) -> f64 {
        let (num, den) = self.ratio();
        let x = ratio_to_f64(&num, &den);
        if self.neg {
            -x
        } else {
            x
        }
    }

    /// The correctly rounded binary32, rounded once from the ratio and not through binary64: what
    /// a perfect `f32` subject returns.
    pub(crate) fn to_f32(self) -> f32 {
        let (num, den) = self.ratio();
        let x = ratio_to_f32(&num, &den);
        if self.neg {
            -x
        } else {
            x
        }
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::conformance::testkit::splitmix;

    /// Python's `float.hex()`, written independently of `parse_hex`.
    fn hex(x: f64) -> String {
        let bits = x.to_bits();
        let sign = if bits >> 63 == 1 { "-" } else { "" };
        let biased = (bits >> 52) & 0x7ff;
        let frac = bits & ((1 << 52) - 1);
        match (biased, frac) {
            (0, 0) => format!("{sign}0x0.0p+0"),
            (0, f) => format!("{sign}0x0.{f:013x}p-1022"),
            (b, f) => format!("{sign}0x1.{f:013x}p{:+}", b as i64 - 1023),
        }
    }

    #[test]
    fn hex_floats_parse_exactly() -> Result<(), String> {
        let cases = [
            ("0x1.0000000000000p+0", 1.0),
            ("0x1.8p+1", 3.0),
            ("-0x1.4p-2", -0.3125),
            ("0x1p-1074", f64::from_bits(1)),
            ("0x0.0000000000001p-1022", f64::from_bits(1)),
            ("0x1.fffffffffffffp+1023", f64::MAX),
            ("0x1.0000000000001p+0", 1.0 + f64::EPSILON),
        ];
        for (s, want) in cases {
            assert_eq!(parse_hex(s)?.to_bits(), want.to_bits(), "{s}");
        }
        assert_eq!(parse_hex("-0x0.0p+0")?.to_bits(), (-0.0f64).to_bits());
        assert_eq!(parse_hex("0x0.0p+0")?.to_bits(), 0);
        for bad in [
            "0x1p-1075",
            "0x1p+1024",
            "0x1.00000000000001p+0",
            "1.0",
            "0x.p+0",
            "0x1.8",
            "0x1.gp+0",
            "inf",
            "nan",
        ] {
            assert!(parse_hex(bad).is_err(), "{bad}");
        }
        Ok(())
    }

    #[test]
    fn hex_round_trips_random_bit_patterns_including_subnormals() -> Result<(), String> {
        let mut state = 7;
        for i in 0..20_000 {
            let mut bits = splitmix(&mut state);
            if i % 4 == 0 {
                bits &= !(0x7ff << 52);
            }
            let x = f64::from_bits(bits);
            if x.is_finite() {
                assert_eq!(parse_hex(&hex(x))?.to_bits(), bits, "{}", hex(x));
                let d = Dyadic::of(x);
                assert_eq!(d.mant as f64 * pow2(i64::from(d.exp)), x.abs());
                assert_eq!(d.neg, x.is_sign_negative());
            }
        }
        Ok(())
    }

    #[test]
    fn decimals_parse_as_exact_integers_and_exponents() -> Result<(), String> {
        let d = |s: &str| s.parse::<Decimal>();
        let want = |neg, mant, exp10| Decimal { neg, mant, exp10 };
        assert_eq!(
            d("5.00000000000000000000000000000e-1")?,
            want(false, 5 * 10u128.pow(29), -30)
        );
        assert_eq!(d("-1.5e3")?, want(true, 15, 2));
        assert_eq!(d("1e5")?, want(false, 1, 5));
        assert_eq!(d("0.00000000000000000000000000000e0")?.mant, 0);
        assert_eq!(d("2.44699470068014011431711870931e-12")?.exp10, -12 - 29);
        assert_eq!(d("+7")?, want(false, 7, 0));
        for bad in ["", "e5", "1.2.3", "1e", "abc", "1e+x", &"9".repeat(39)] {
            assert!(d(bad).is_err(), "{bad}");
        }
        Ok(())
    }

    #[test]
    fn ratio_rounds_to_nearest_even_once() {
        let r = |n: u128, d: u128| ratio_to_f64(&BigUint::from(n), &BigUint::from(d));
        assert_eq!(r(1, 3), 1.0 / 3.0);
        assert_eq!(r(0, 5), 0.0);
        let two53 = 1u128 << 53;
        // Ties: 2^53 + 1 is halfway between 2^53 and 2^53 + 2 (even mantissa: down); 2^53 + 3 up.
        assert_eq!(r(two53 + 1, 1), two53 as f64);
        assert_eq!(r(two53 + 3, 1), (two53 + 4) as f64);
        // Off the tie: 2^53 + 1.5 is nearer 2^53 + 2 than 2^53.
        assert_eq!(r(2 * two53 + 3, 2), (two53 + 2) as f64);
        assert_eq!(r((1 << 120) * 3, 1 << 119), 6.0);
        // Subnormal: 1.5 * 2^-1074 ties to the even 2 * 2^-1074; 0.5 * 2^-1074 ties to 0.
        let tiny = |k: usize, d: usize| {
            ratio_to_f64(&BigUint::from(k), &(BigUint::from(1u32) << (1074 + d)))
        };
        assert_eq!(tiny(3, 1), f64::from_bits(2));
        assert_eq!(tiny(1, 1), 0.0);
        assert_eq!(tiny(3, 0), f64::from_bits(3));
        assert_eq!(tiny(5, 1), f64::from_bits(2));
        // The largest finite is 2^1024 - 2^971; the tie above it rounds to the even 2^1024: overflow.
        let one = BigUint::from(1u32);
        let top = &one << 1024usize;
        assert_eq!(ratio_to_f64(&(&top - (&one << 971usize)), &one), f64::MAX);
        assert_eq!(
            ratio_to_f64(&(&top - (&one << 970usize) - &one), &one),
            f64::MAX
        );
        assert_eq!(
            ratio_to_f64(&(&top - (&one << 970usize)), &one),
            f64::INFINITY
        );
    }

    #[test]
    fn ratio_rounds_to_binary32_once() {
        let r = |n: u128, d: u128| ratio_to_f32(&BigUint::from(n), &BigUint::from(d));
        let (two24, one) = (1u128 << 24, BigUint::from(1u32));
        // Ties go to the even neighbour.
        assert_eq!(r(two24 + 1, 1), 16_777_216.0);
        assert_eq!(r(two24 + 3, 1), 16_777_220.0);
        // Just above the tie 2^24 + 1 by 2^-60: binary64 holds it as the tie itself, and rounding
        // that to binary32 again would give 2^24. One rounding gives 2^24 + 2.
        assert_eq!(r(((two24 + 1) << 60) + 1, 1 << 60), 16_777_218.0);
        assert_eq!(r(1, 3).to_bits(), (1.0f32 / 3.0).to_bits());
        // Subnormal: 1.5 * 2^-149 ties to the even 2 * 2^-149; 0.5 * 2^-149 ties to 0.
        let tiny =
            |k: u32, d: usize| ratio_to_f32(&BigUint::from(k), &(BigUint::from(1u32) << (149 + d)));
        assert_eq!(tiny(3, 1), f32::from_bits(2));
        assert_eq!(tiny(1, 1), 0.0);
        assert_eq!(tiny(3, 0), f32::from_bits(3));
        // The largest finite is 2^128 - 2^104; the tie above it rounds to the even 2^128: overflow.
        let big = |k: u32| BigUint::from(1u32) << k;
        let (max, tie) = (big(128) - big(104), big(128) - big(103));
        assert_eq!(ratio_to_f32(&max, &one), f32::MAX);
        assert_eq!(ratio_to_f32(&(&tie - &one), &one), f32::MAX);
        assert_eq!(ratio_to_f32(&tie, &one), f32::INFINITY);
        // Integers a binary32 holds: the quotient of two binary32 values is correctly rounded.
        let mut state = 5;
        for _ in 0..2000 {
            let n = crate::conformance::testkit::splitmix(&mut state) >> 40;
            let d = (crate::conformance::testkit::splitmix(&mut state) >> 40) | 1;
            let want = n as f32 / d as f32;
            assert_eq!(r(n.into(), d.into()).to_bits(), want.to_bits(), "{n}/{d}");
        }
    }

    #[test]
    fn roots_reach_below_the_square_underflow() {
        let big = |k: usize| BigUint::from(1u32) << k;
        let one = BigUint::from(1u32);
        assert_eq!(sqrt_ratio(&one, &BigUint::from(4u32)), 0.5);
        assert_eq!(sqrt_ratio(&BigUint::default(), &one), 0.0);
        assert_eq!(sqrt_ratio(&big(2000), &one), pow2(1000));
        assert_eq!(sqrt_ratio(&big(2001), &one), pow2(1000) * 2f64.sqrt());
        assert_eq!(sqrt_ratio(&one, &big(2140)), f64::from_bits(1 << 4));
        assert_eq!(sqrt_ratio(&one, &big(2200)), 0.0);
        // 2^-1075 is the tie between 0 and the smallest subnormal (even: 0); sqrt(2) 2^-1075 is above it.
        assert_eq!(sqrt_ratio(&one, &big(2150)), 0.0);
        assert_eq!(
            sqrt_ratio(&BigUint::from(2u32), &big(2150)),
            f64::from_bits(1)
        );
        assert_eq!(sqrt_ratio(&big(3000), &one), f64::INFINITY);
        // The overflow tie: 2^1024 - 2^970 is halfway between MAX (odd) and 2^1024; a hair below it
        // is MAX.
        let mid = big(1024) - big(970);
        assert_eq!(sqrt_ratio(&(&mid * &mid), &one), f64::INFINITY);
        let below = mid - 1u32;
        assert_eq!(sqrt_ratio(&(&below * &below), &one), f64::MAX);
        assert_eq!(sqrt_ratio(&big(3), &big(1)), 2.0);
    }

    /// Whether `r` is `sqrt(num / den)` rounded to nearest: `num / den` lies between the squares of
    /// the midpoints to `r`'s neighbours, in integers over `2^2150` (`r` finite, positive, not the
    /// largest finite).
    fn is_nearest_root(r: f64, num: &BigUint, den: &BigUint) -> bool {
        let scaled = |x: f64| {
            let d = Dyadic::of(x);
            BigUint::from(d.mant) << (d.exp + 1074) as usize
        };
        let (below, here, above) = (scaled(r.next_down()), scaled(r), scaled(r.next_up()));
        let square = |x: BigUint| &x * &x * den;
        let target = num << 2150usize;
        square(&below + &here) <= target && target <= square(&here + &above)
    }

    #[test]
    fn roots_are_correctly_rounded() {
        // A 1 ulp miss of the old float route: 15661141135.162043.
        let big = |s: &str| BigUint::parse_bytes(s.as_bytes(), 10).unwrap_or_default();
        let (n, d) = (
            big("131610346136799443493804746066869303032"),
            big("536590802857326633"),
        );
        assert_eq!(sqrt_ratio(&n, &d).to_bits(), 0x420d_2bd0_6479_4bdd);
        // Ties: (2^53 + 1)^2 / 2^106 is halfway between 1 and 1 + 2^-52 (even: 1); with 3 it is
        // halfway between 1 + 2^-52 and 1 + 2^-51 (even: the latter).
        let tie = |k: u128| {
            sqrt_ratio(
                &BigUint::from(((1u128 << 53) + k).pow(2)),
                &(BigUint::from(1u32) << 106usize),
            )
        };
        assert_eq!(tie(1), 1.0);
        assert_eq!(tie(3), 1.0 + 2.0 * f64::EPSILON);
        assert_eq!(tie(2), 1.0 + f64::EPSILON);
        let mut state = 5;
        let (mut subnormal, mut checked) = (0, 0);
        for i in 0..20_000 {
            let wide = (u128::from(splitmix(&mut state)) << 64) | u128::from(splitmix(&mut state));
            let num = BigUint::from(wide) << (splitmix(&mut state) % 3000) as usize;
            let den =
                BigUint::from(splitmix(&mut state) | 1) << (splitmix(&mut state) % 3000) as usize;
            let r = sqrt_ratio(&num, &den);
            if r.is_finite() && r > 0.0 && r < f64::MAX {
                assert!(is_nearest_root(r, &num, &den), "{i}: {num} / {den} = {r:e}");
                checked += 1;
                subnormal += usize::from(r < f64::MIN_POSITIVE);
            }
        }
        assert!(checked > 10_000 && subnormal > 100, "{checked} {subnormal}");
    }

    #[test]
    fn decimals_round_like_rust_parse_on_random_strings() -> Result<(), String> {
        let mut state = 11;
        for _ in 0..5_000 {
            let digits = splitmix(&mut state) % 10u64.pow(17);
            let exp = (splitmix(&mut state) % 640) as i32 - 330;
            let neg = if splitmix(&mut state) & 1 == 1 {
                "-"
            } else {
                ""
            };
            let s = format!(
                "{neg}{}.{:016}e{exp}",
                digits / 10u64.pow(16),
                digits % 10u64.pow(16)
            );
            let want: f64 = s.parse().map_err(|e| format!("{s}: {e}"))?;
            assert_eq!(
                s.parse::<Decimal>()?.to_f64().to_bits(),
                want.to_bits(),
                "{s}"
            );
            // The same decimal at binary32: Rust's parser rounds once, as `to_f32` must (a
            // binary64 in between would round twice).
            let want: f32 = s.parse().map_err(|e| format!("{s}: {e}"))?;
            let got = s.parse::<Decimal>()?.to_f32();
            assert_eq!(got.to_bits(), want.to_bits(), "{s}");
        }
        // 1 + 2^-24 + 2^-60 is above the binary32 tie 1 + 2^-24, which binary64 holds as the tie.
        let above = "1.00000005960464477625798673799e0".parse::<Decimal>()?;
        assert_eq!(above.to_f32(), 1.0 + f32::EPSILON);
        assert_eq!(above.to_f64() as f32, 1.0);
        Ok(())
    }
}
