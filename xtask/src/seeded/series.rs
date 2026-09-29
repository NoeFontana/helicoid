//! The series constants: `conformance/corpus/coeff_series.jsonl` read at run time, each exact
//! `"num/den"` rounded once, in integers, at the precision it is used at (0003 item 6: no
//! `f64 -> f32` double rounding). Nothing is typed.

use std::path::Path;

use helicoid_linalg::{Precision, Real};
use num_bigint::BigUint;
use serde::Deserialize;

use super::kernel::Coeff;
use crate::conformance::number::{ratio_to_f32, ratio_to_f64};

/// The series file, in the corpus directory.
pub(crate) const FILE: &str = "coeff_series.jsonl";

/// `"num/den"` (`"-1/48"`, a denominator is required) as sign, numerator and denominator.
fn parse_rational(text: &str) -> Result<(bool, BigUint, BigUint), String> {
    let bad = || format!("bad rational `{text}`");
    let (num, den) = text.split_once('/').ok_or_else(bad)?;
    let (neg, digits) = match num.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, num),
    };
    let big = |s: &str| {
        let digits_only = !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
        digits_only
            .then(|| BigUint::parse_bytes(s.as_bytes(), 10))
            .flatten()
    };
    let (num, den) = (big(digits).ok_or_else(bad)?, big(den).ok_or_else(bad)?);
    if den.bits() == 0 {
        return Err(bad());
    }
    Ok((neg, num, den))
}

/// The rational `text` correctly rounded at `S`, once, in integers.
pub(crate) fn rational<S: Real>(text: &str) -> Result<S, String> {
    let (neg, num, den) = parse_rational(text)?;
    let x = match S::PRECISION {
        Precision::F64 => ratio_to_f64(&num, &den),
        // Exact in binary64: `Real::lit` receives a value `S` holds.
        Precision::F32 => f64::from(ratio_to_f32(&num, &den)),
    };
    Ok(S::lit(if neg { -x } else { x }))
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Line {
    id: u64,
    coeff: String,
    branch: String,
    prefactor: String,
    series: Vec<String>,
}

/// The series of `k, a, b, c, d, e` in the branch variable `θ²`, rounded at `S`; one non-empty
/// series per coefficient, all of one length.
#[derive(Clone, Debug)]
pub(crate) struct Series<S>(Vec<Vec<S>>);

/// The committed series file, compiled in: until the `f32` sweep generates its own constants
/// (`docs/decisions/0016` item 3) the seeded kernels at `f32` read the exact rationals from it.
const COMMITTED: &str = include_str!("../../../conformance/corpus/coeff_series.jsonl");

impl<S: Real> Series<S> {
    pub(crate) fn load(corpus: &Path) -> Result<Self, String> {
        let path = corpus.join(FILE);
        let text =
            std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        Self::parse(&text, &path.display().to_string())
    }

    /// The compiled-in series, rounded at `S`.
    pub(crate) fn committed() -> Result<Self, String> {
        Self::parse(COMMITTED, FILE)
    }

    /// The series of `text`, a `coeff_series.jsonl` read from `origin`.
    fn parse(text: &str, origin: &str) -> Result<Self, String> {
        let mut found: [Option<Vec<S>>; 6] = Default::default();
        for (i, line) in text.lines().enumerate() {
            let at = |e: String| format!("{origin}:{}: {e}", i + 1);
            let l: Line = serde_json::from_str(line).map_err(|e| at(e.to_string()))?;
            if l.id != i as u64 {
                return Err(at(format!("id {} on line {}", l.id, i + 1)));
            }
            // `r` has its own branch variable and no seeded kernel.
            let Some(c) = Coeff::of_name(&l.coeff) else {
                continue;
            };
            if l.branch != "theta^2" || l.prefactor != "1" {
                return Err(at(format!(
                    "`{}` is not a plain series in theta^2",
                    l.coeff
                )));
            }
            if found[c.index()].is_some() {
                return Err(at(format!("`{}` appears twice", l.coeff)));
            }
            let terms = l.series.iter().map(|s| rational::<S>(s));
            found[c.index()] = Some(terms.collect::<Result<_, _>>().map_err(at)?);
        }
        let mut all = Vec::new();
        for c in Coeff::ALL {
            let series = found[c.index()]
                .take()
                .ok_or_else(|| format!("{origin}: no series for `{}`", c.name()))?;
            if series.is_empty() || series.len() != all.first().map_or(series.len(), Vec::len) {
                return Err(format!(
                    "{origin}: `{}` has {} terms",
                    c.name(),
                    series.len()
                ));
            }
            all.push(series);
        }
        Ok(Self(all))
    }

    /// The terms of `c`, lowest power first.
    pub(crate) fn of(&self, c: Coeff) -> &[S] {
        &self.0[c.index()]
    }

    /// The number of terms every series has.
    pub(crate) fn terms(&self) -> usize {
        self.of(Coeff::K).len()
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::conformance::testkit::{splitmix, Scratch};

    fn committed() -> Result<String, String> {
        let path = crate::conformance::corpus_dir()?.join("coeff_series.jsonl");
        std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))
    }

    /// A random integer of exactly `bits` bits, odd.
    fn big(bits: usize, state: &mut u64) -> BigUint {
        let mut n = BigUint::default();
        while (n.bits() as usize) < bits {
            n = (n << 64usize) | BigUint::from(splitmix(state));
        }
        let excess = n.bits() as usize - bits;
        (n >> excess) | BigUint::from(1u32)
    }

    /// `num / den` through the decimal expansion to 80 digits and Rust's own (correctly rounded)
    /// parser: nothing shared with `ratio_to_f64`.
    fn via_decimal<T: std::str::FromStr>(num: &BigUint, den: &BigUint) -> Result<T, String> {
        let digits = |n: &BigUint| n.to_string().len() as i64;
        let shift = 80 + digits(den) - digits(num);
        let ten = |k: i64| BigUint::from(10u32).pow(k.unsigned_abs() as u32);
        let q = if shift >= 0 {
            num * ten(shift) / den
        } else {
            num / (den * ten(shift))
        };
        format!("{q}e{}", -shift)
            .parse()
            .map_err(|_| format!("{q}e{}", -shift))
    }

    #[test]
    fn a_rational_rounds_like_division_and_like_the_decimal_parser() -> Result<(), String> {
        let mut state = 3;
        for _ in 0..3000 {
            let (n, d) = (splitmix(&mut state) >> 11, (splitmix(&mut state) >> 11) | 1);
            let x = rational::<f64>(&format!("{n}/{d}"))?;
            assert_eq!(x.to_bits(), (n as f64 / d as f64).to_bits(), "{n}/{d}");
            assert_eq!(rational::<f64>(&format!("-{n}/{d}"))?, -x);
            // Both operands are binary32 values: the quotient is correctly rounded there too.
            let (n, d) = (n >> 29, (d >> 29) | 1);
            let y = rational::<f32>(&format!("{n}/{d}"))?;
            assert_eq!(y.to_bits(), (n as f32 / d as f32).to_bits(), "{n}/{d}");
        }
        let (mut overflow, mut subnormal) = (0, 0);
        let (mut overflow32, mut subnormal32) = (0, 0);
        for i in 0..3000 {
            let (num, den) = (
                big(1 + splitmix(&mut state) as usize % 1300, &mut state),
                big(1 + splitmix(&mut state) as usize % 1300, &mut state),
            );
            let text = format!("{num}/{den}");
            let want = via_decimal::<f64>(&num, &den)?;
            assert_eq!(
                rational::<f64>(&text)?.to_bits(),
                want.to_bits(),
                "{i}: {text}"
            );
            overflow += usize::from(want.is_infinite());
            subnormal += usize::from(want != 0.0 && want < f64::MIN_POSITIVE);
            let want = via_decimal::<f32>(&num, &den)?;
            assert_eq!(
                rational::<f32>(&text)?.to_bits(),
                want.to_bits(),
                "{i}: {text}"
            );
            overflow32 += usize::from(want.is_infinite());
            subnormal32 += usize::from(want != 0.0 && want < f32::MIN_POSITIVE);
        }
        assert!(overflow > 20 && subnormal > 5, "{overflow} {subnormal}");
        assert!(
            overflow32 > 100 && subnormal32 > 5,
            "{overflow32} {subnormal32}"
        );
        Ok(())
    }

    /// `(2^70 + k 2^46 ± 1) / 2^70` is `1 + k 2^-24` off by `2^-70`, under half a binary64 ulp: it
    /// rounds to that binary32 tie in binary64, and the tie goes to even. Rounded once at binary32
    /// it is on its own side, one binary32 ulp from where the double rounding lands.
    #[test]
    fn a_rational_next_to_a_binary32_tie_rounds_once_at_binary32() -> Result<(), String> {
        let one = BigUint::from(1u32);
        let den = &one << 70usize;
        let text = |k: u32, above: bool| {
            let num = &den + (BigUint::from(k) << 46usize);
            format!("{}/{den}", if above { num + &one } else { num - &one })
        };
        let ulp = f32::EPSILON;
        // Above the tie between 1 and 1 + 2^-23: up, where double rounding ties to the even 1.
        let (above, below) = (text(1, true), text(3, false));
        assert_eq!(rational::<f32>(&above)?, 1.0 + ulp);
        assert_eq!(rational::<f64>(&above)? as f32, 1.0);
        // Below the tie between 1 + 2^-23 and 1 + 2^-22: down, where double rounding ties to even up.
        assert_eq!(rational::<f32>(&below)?, 1.0 + ulp);
        assert_eq!(rational::<f64>(&below)? as f32, 1.0 + 2.0 * ulp);
        assert_eq!(rational::<f32>(&format!("-{above}"))?, -(1.0 + ulp));
        Ok(())
    }

    #[test]
    fn a_parse_error_names_its_origin_and_line() {
        // The `f32` kernels read the compiled-in text: a broken one is reported as this, not lost.
        let e = Series::<f32>::parse("{}\n", FILE).err().unwrap_or_default();
        assert!(e.starts_with("coeff_series.jsonl:1: "), "{e}");
        assert!(Series::<f32>::committed().is_ok());
    }

    #[test]
    fn a_rational_needs_a_sign_digits_and_a_nonzero_denominator() -> Result<(), String> {
        assert_eq!(rational::<f64>("1/1")?, 1.0);
        assert_eq!(rational::<f64>("-1/48")?, -1.0 / 48.0);
        assert_eq!(rational::<f64>("0/7")?.to_bits(), 0);
        assert_eq!(rational::<f32>("-1/48")?, -1.0 / 48.0);
        for bad in [
            "", "1", "1/", "/1", "1/0", "+1/2", "1/-2", "1/2/3", "a/2", "1.5/2", " 1/2",
        ] {
            assert!(rational::<f64>(bad).is_err(), "{bad}");
        }
        Ok(())
    }

    /// Exact term `j` of the series of a coefficient in `z = θ²` (`docs/maths/coefficients.md`
    /// CO.2), as `(negative, num, den)`.
    fn general_term(c: Coeff, j: u32) -> (bool, BigUint, BigUint) {
        let fact = |n: u32| (1..=u64::from(n)).fold(BigUint::from(1u32), |a, i| a * i);
        let (num, den) = match c {
            Coeff::K => (
                1,
                (BigUint::from(1u32) << (2 * j as usize + 1)) * fact(2 * j + 1),
            ),
            Coeff::A => (1, fact(2 * j + 2)),
            Coeff::B => (1, fact(2 * j + 3)),
            Coeff::D => (1, fact(2 * j + 4)),
            Coeff::E => (j + 1, fact(2 * j + 5)),
            Coeff::C => (0, BigUint::default()),
        };
        (j % 2 == 1, BigUint::from(num), den)
    }

    #[test]
    fn the_committed_series_are_the_general_terms() -> Result<(), String> {
        let text = committed()?;
        let mut checked = 0;
        for line in text.lines() {
            let l: Line = serde_json::from_str(line).map_err(|e| e.to_string())?;
            let Some(c) = Coeff::of_name(&l.coeff) else {
                continue;
            };
            assert_eq!(l.series.len(), 16, "{}", l.coeff);
            for (j, term) in l.series.iter().enumerate() {
                let (neg, num, den) = parse_rational(term)?;
                if c == Coeff::C {
                    // The leading four of `NUMERICS.md` §4; the tail is Bernoulli numbers.
                    let lead = [(1u32, 12u32), (1, 720), (1, 30240), (1, 1_209_600)];
                    if let Some(&(n, d)) = lead.get(j) {
                        assert!(!neg && num * d == BigUint::from(n) * &den, "c[{j}]");
                    }
                    continue;
                }
                let (want_neg, want_num, want_den) = general_term(c, j as u32);
                let same = num * &want_den == want_num * &den;
                assert!(neg == want_neg && same, "{}[{j}] = {term}", l.coeff);
                checked += 1;
            }
        }
        assert_eq!(checked, 5 * 16);
        Ok(())
    }

    #[test]
    fn the_loader_refuses_a_missing_series_a_wrong_id_and_a_series_not_plain() -> Result<(), String>
    {
        let dir = Scratch::new("series");
        std::fs::create_dir_all(&dir.0).map_err(|e| e.to_string())?;
        let good = committed()?;
        let load = |text: &str| -> Result<String, String> {
            std::fs::write(dir.0.join("coeff_series.jsonl"), text).map_err(|e| e.to_string())?;
            Ok(Series::<f64>::load(&dir.0).err().unwrap_or_default())
        };
        assert_eq!(load(&good)?, "");
        let renamed = good.replace("\"coeff\":\"e\"", "\"coeff\":\"z\"");
        assert!(load(&renamed)?.contains("no series for `e`"));
        let plain = "\"prefactor\":\"1\"";
        assert!(load(&good.replace(plain, "\"prefactor\":\"2\""))?.contains("plain series"));
        assert!(load(&good.replacen("\"id\":1", "\"id\":9", 1))?.contains("id 9 on line 2"));
        assert!(load("")?.contains("no series for `k`"));
        // One series longer than the rest, and six empty ones: no kernel could run on either.
        let longer = good.replacen("\"series\":[", "\"series\":[\"1/1\",", 1);
        assert!(load(&longer)?.contains("`a` has 16 terms"));
        let empty = |i: usize, c: &str| {
            let line = "\"branch\":\"theta^2\",\"prefactor\":\"1\",\"series\":[]";
            format!("{{\"id\":{i},\"coeff\":\"{c}\",{line}}}\n")
        };
        let none: String = Coeff::ALL
            .iter()
            .enumerate()
            .map(|(i, c)| empty(i, c.name()))
            .collect();
        assert!(load(&none)?.contains("`k` has 0 terms"));
        Ok(())
    }
}
