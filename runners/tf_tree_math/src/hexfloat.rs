//! Python's `float.hex()` and its inverse: the format of the corpus's inputs and of this runner's
//! outputs. One spelling per value, so a round trip is byte-exact and no rounding can enter.
//! `nan`, `inf` and `-inf` carry a non-finite output to the harness, which counts it.

/// `[-]0x1.<13 hex digits>p<±exp>` for a normal number, `[-]0x0.<13 digits>p-1022` for a
/// subnormal, `[-]0x0.0p+0` for a zero, `inf`, `-inf`, `nan`.
pub(crate) fn format(x: f64) -> String {
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf" } else { "inf" }.into();
    }
    let bits = x.to_bits();
    let sign = if bits >> 63 == 1 { "-" } else { "" };
    let frac = bits & ((1 << 52) - 1);
    let biased = i32::try_from((bits >> 52) & 0x7ff).unwrap_or(0);
    match (biased, frac) {
        (0, 0) => format!("{sign}0x0.0p+0"),
        (0, f) => format!("{sign}0x0.{f:013x}p-1022"),
        (e, f) => format!("{sign}0x1.{f:013x}p{:+}", e - 1023),
    }
}

/// The value [`format()`] spells `s` from; any other spelling is an error.
pub(crate) fn parse(s: &str) -> Result<f64, String> {
    match s {
        "nan" => return Ok(f64::NAN),
        "inf" => return Ok(f64::INFINITY),
        "-inf" => return Ok(f64::NEG_INFINITY),
        _ => {}
    }
    let bad = || format!("bad hex float `{s}`");
    let (neg, rest) = s.strip_prefix('-').map_or((false, s), |r| (true, r));
    let (lead, rest) = match (rest.strip_prefix("0x1."), rest.strip_prefix("0x0.")) {
        (Some(r), _) => (1, r),
        (None, Some(r)) => (0, r),
        (None, None) => return Err(bad()),
    };
    let (digits, exp) = rest.split_once('p').ok_or_else(bad)?;
    let exp: i32 = exp.parse().map_err(|_| bad())?;
    let frac = |d: &str| match d.len() == 13 && d.bytes().all(|b| b.is_ascii_hexdigit()) {
        true => u64::from_str_radix(d, 16).map_err(|_| bad()),
        false => Err(bad()),
    };
    let magnitude = match (lead, digits, exp) {
        (0, "0", 0) => 0,
        (0, d, -1022) => frac(d)?,
        (1, d, -1022..=1023) => (u64::try_from(exp + 1023).map_err(|_| bad())? << 52) | frac(d)?,
        _ => return Err(bad()),
    };
    // A subnormal spelled with no digit set is a zero spelled the long way.
    match (lead, magnitude) {
        (0, 0) if digits != "0" => Err(bad()),
        _ => Ok(f64::from_bits((u64::from(neg) << 63) | magnitude)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Each pair is `float.hex()` of the value, worked out by hand from its bits.
    const CASES: [(u64, &str); 8] = [
        (0x3ff0_0000_0000_0000, "0x1.0000000000000p+0"),
        (0x4008_0000_0000_0000, "0x1.8000000000000p+1"),
        (0xc004_0000_0000_0000, "-0x1.4000000000000p+1"),
        (0x0000_0000_0000_0000, "0x0.0p+0"),
        (0x8000_0000_0000_0000, "-0x0.0p+0"),
        (0x0000_0000_0000_0001, "0x0.0000000000001p-1022"),
        (0x0010_0000_0000_0000, "0x1.0000000000000p-1022"),
        (0x7fef_ffff_ffff_ffff, "0x1.fffffffffffffp+1023"),
    ];

    #[test]
    fn the_spelling_is_pythons_and_reads_back_bit_for_bit() -> Result<(), String> {
        for (bits, text) in CASES {
            assert_eq!(format(f64::from_bits(bits)), text);
            assert_eq!(parse(text)?.to_bits(), bits, "{text}");
        }
        Ok(())
    }

    #[test]
    fn non_finite_values_have_their_own_words() -> Result<(), String> {
        assert_eq!(format(f64::NAN), "nan");
        assert_eq!(format(f64::INFINITY), "inf");
        assert_eq!(format(f64::NEG_INFINITY), "-inf");
        assert!(parse("nan")?.is_nan());
        assert_eq!(parse("-inf")?.to_bits(), f64::NEG_INFINITY.to_bits());
        Ok(())
    }

    #[test]
    fn any_other_spelling_is_refused() {
        for bad in [
            "1.0",
            "0x1p+0",
            "0x1.8p+1",
            "0x1.0000000000000p+1024",
            "0x1.0000000000000p-1023",
            "0x0.0000000000000p-1022",
            "0x0.0p-1022",
            "0x1.000000000000gp+0",
            "0x1.0000000000000",
            "",
        ] {
            assert!(parse(bad).is_err(), "{bad}");
        }
    }
}
