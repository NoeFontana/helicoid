//! The forward error of `docs/NUMERICS.md` §11, per record, in units of `u`, formed exactly.
//!
//! For each output field, `F = ‖ŷ − y‖ / (max(‖y‖, floor) · u)`, the norm over the whole field
//! (vectors: 2-norm, matrices: Frobenius, which is the same sum of squares). A record's score is
//! the largest `F` over its fields. `ŷ` is the subject's binary64 output, `y` the 30-digit
//! reference: every entry is scaled to an integer over one `2^p · 5^q`, the sums of squares are
//! integers, and only the root of their ratio is rounded, once, to binary64
//! (`number::sqrt_ratio`). Parsing `y` to `f64` first would quantize `F` to an ulp
//! (`docs/maths/error-analysis.md` EA.19(b)).
//!
//! Readings where `NUMERICS.md` §11 is silent, each in [`TABLE`]; 0014 (draft) proposes them and
//! settles none:
//!
//! - **Floor.** §11 names 1 for rotations and Jacobians and the ρ-scale for translations, and none
//!   for a tangent output or a coefficient. Those get the smallest normal, `2^-1022`: a relative
//!   error, except that an exact zero is not divided by and a subnormal output is judged by its
//!   quantization (a half-spacing rounding reads at most 1 per entry). A floor of 1 would make
//!   `theta:1e-k` measure absolute error, `EA.4(b)`. The ρ-scale is the 2-norm of the translation
//!   entries of the input tangent.
//! - **Sign.** A quaternion output is compared after sign alignment, which is the smaller of
//!   `‖ŷ − y‖` and `‖ŷ + y‖`. `so3_log`'s tangent uses the same form when the input quaternion has
//!   `w = +0` (bit pattern zero), and only there. The SE_N(3) tangent does not: `−τ` is not the
//!   logarithm the input's sign selects, its translation block differs.
//! - **Fields.** No aggregate says which field of a record scored; the maximum is taken.
//! - **Tangents** are one vector, rotation block and translation block together, as §11 says: at
//!   `rho:1e4` the translation block dominates `‖y‖` and a rotation error that is many `u` of `φ`
//!   can read under 1.
//! - **Overflow.** An `F` that does not fit in binary64 (a gross error against a zero or tiny
//!   reference) is [`Score::NonFinite`]: a finite output that wrong fails the run like a NaN.
//! - **`f32`** is refused: the corpus inputs are binary64 and `docs/PHASE1.md` §4.4 does not say
//!   whether an `f32` subject receives them rounded.

use helicoid_linalg::Precision;
use num_bigint::{BigInt, BigUint, Sign};

use super::corpus::{Record, Tensor};
use super::number::{sqrt_ratio, Decimal, Dyadic};
use super::subject::Output;

/// The scale `‖y‖` is floored at.
pub(crate) enum Floor {
    /// 1: rotations, quaternions, Jacobians, adjoints, rotated points.
    Unit,
    /// `2^-1022`: tangents and coefficients, whose scale is their own.
    Tiny,
    /// The 2-norm of `input[skip..]`: the translation scale of an `Exp` (rotation-first tangent).
    Scale { input: &'static str, skip: usize },
}

/// How the sign of the output enters the error.
pub(crate) enum SignRule {
    Fixed,
    /// `min(‖ŷ − y‖, ‖ŷ + y‖)`.
    Align,
    /// `Align` when the input `q` has `w = +0`, `Fixed` otherwise.
    AlignAtW0,
}

pub(crate) struct FieldRule {
    pub(crate) field: &'static str,
    floor: Floor,
    sign: SignRule,
}

pub(crate) enum Rule {
    Forward(&'static [FieldRule]),
    /// `NUMERICS.md` §11: a matrix input is scored by backward error only, which is owed; until
    /// then its outputs are checked for shape and finiteness (`finiteness`).
    BackwardOnly,
}

impl Rule {
    /// The score of `out` against the reference of `rec` under this rule.
    pub(crate) fn score(
        &self,
        rec: &Record,
        out: &Output,
        precision: Precision,
    ) -> Result<Score, String> {
        match self {
            Rule::Forward(fields) => score(fields, rec, out, precision),
            Rule::BackwardOnly => finiteness(rec, out),
        }
    }
}

const fn field(field: &'static str, floor: Floor, sign: SignRule) -> FieldRule {
    FieldRule { field, floor, sign }
}

const COEFF: &[FieldRule] = &[
    field("value", Floor::Tiny, SignRule::Fixed),
    field("d_branch", Floor::Tiny, SignRule::Fixed),
];
/// One output of a coefficient on its own: the curve of a defect is a field's, not the record's
/// maximum. `b`'s value errs by `θ^-2` and its derivative through `Dual` by `θ^-4`
/// (`docs/maths/coefficients.md` CO.6, CO.14).
pub(crate) static COEFF_VALUE: Rule = Forward(&[field("value", Floor::Tiny, SignRule::Fixed)]);
pub(crate) static COEFF_D_BRANCH: Rule =
    Forward(&[field("d_branch", Floor::Tiny, SignRule::Fixed)]);
const JAC: &[FieldRule] = &[field("J", Floor::Unit, SignRule::Fixed)];
const AD: &[FieldRule] = &[field("Ad", Floor::Unit, SignRule::Fixed)];
const TANGENT: &[FieldRule] = &[field("tau", Floor::Tiny, SignRule::Fixed)];
const SO3_EXP: &[FieldRule] = &[field("q", Floor::Unit, SignRule::Align)];
const SO3_LOG: &[FieldRule] = &[field("phi", Floor::Tiny, SignRule::AlignAtW0)];
const SO3_ACT: &[FieldRule] = &[field("Rp", Floor::Unit, SignRule::Fixed)];
const SEN3_EXP: &[FieldRule] = &[
    field("q", Floor::Unit, SignRule::Align),
    field(
        "x",
        Floor::Scale {
            input: "tau",
            skip: 3,
        },
        SignRule::Fixed,
    ),
];
const SO2_EXP: &[FieldRule] = &[field("z", Floor::Unit, SignRule::Fixed)];
const SO2_LOG: &[FieldRule] = &[field("theta", Floor::Tiny, SignRule::Fixed)];
const SE2_EXP: &[FieldRule] = &[
    field("z", Floor::Unit, SignRule::Fixed),
    field(
        "t",
        Floor::Scale {
            input: "tau",
            skip: 1,
        },
        SignRule::Fixed,
    ),
];

use Rule::{BackwardOnly, Forward};

/// The rule of every corpus family; `sen3_*` share one row over `_n1`, `_n2`, `_n3`.
const TABLE: &[(&str, Rule)] = &[
    ("coeff_k", Forward(COEFF)),
    ("coeff_a", Forward(COEFF)),
    ("coeff_b", Forward(COEFF)),
    ("coeff_c", Forward(COEFF)),
    ("coeff_d", Forward(COEFF)),
    ("coeff_e", Forward(COEFF)),
    ("coeff_r", Forward(COEFF)),
    ("so3_exp", Forward(SO3_EXP)),
    ("so3_log", Forward(SO3_LOG)),
    ("so3_act", Forward(SO3_ACT)),
    ("so3_from_matrix", BackwardOnly),
    ("so3_jr", Forward(JAC)),
    ("so3_jl", Forward(JAC)),
    ("so3_jr_inv", Forward(JAC)),
    ("so3_jl_inv", Forward(JAC)),
    ("sen3_exp", Forward(SEN3_EXP)),
    ("sen3_log", Forward(TANGENT)),
    ("sen3_ad", Forward(AD)),
    ("sen3_jr", Forward(JAC)),
    ("sen3_jl", Forward(JAC)),
    ("sen3_jr_inv", Forward(JAC)),
    ("sen3_jl_inv", Forward(JAC)),
    ("so2_exp", Forward(SO2_EXP)),
    ("so2_log", Forward(SO2_LOG)),
    ("se2_exp", Forward(SE2_EXP)),
    ("se2_log", Forward(TANGENT)),
    ("se2_ad", Forward(AD)),
    ("se2_jr", Forward(JAC)),
    ("se2_jl", Forward(JAC)),
    ("se2_jr_inv", Forward(JAC)),
    ("se2_jl_inv", Forward(JAC)),
];

/// The rule of function id `fn_id`, if the table has one.
pub(crate) fn rule(fn_id: &str) -> Option<&'static Rule> {
    let family = match fn_id.rsplit_once("_n") {
        Some((f, n)) if f.starts_with("sen3_") && matches!(n, "1" | "2" | "3") => f,
        _ => fn_id,
    };
    TABLE.iter().find(|(f, _)| *f == family).map(|(_, r)| r)
}

/// The score of one record.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Score {
    /// The record's error, in units of `u`.
    Finite(f64),
    /// Some output entry is NaN or infinite, or the error of a finite output overflows binary64;
    /// no error is reported.
    NonFinite,
    /// A finite output of a function scored by backward error only, which is owed: shape and
    /// finiteness are all that was checked.
    Unscored,
}

/// `log2(1/u)`: `f32` is refused (see the module docs).
pub(crate) fn unit_bits(precision: Precision) -> Result<usize, String> {
    match precision {
        Precision::F64 => Ok(53),
        Precision::F32 => Err(
            "f32 is not supported: the corpus inputs are binary64 and PHASE1 §4.4 \
             does not say whether an f32 subject receives them rounded"
                .into(),
        ),
    }
}

/// The reference tensor of `field` and the subject's values for it, of the same length.
fn pair<'a>(
    rec: &'a Record,
    out: &'a Output,
    field: &str,
) -> Result<(&'a Tensor<Decimal>, &'a [f64]), String> {
    let y = rec
        .reference
        .get(field)
        .ok_or_else(|| format!("the reference has no `{field}`"))?;
    let yh = out
        .get(field)
        .ok_or_else(|| format!("the subject returned no `{field}`"))?;
    if yh.len() != y.data.len() {
        return Err(format!(
            "`{field}`: the subject returned {} values, the reference has {} (shape {:?})",
            yh.len(),
            y.data.len(),
            y.shape
        ));
    }
    Ok((y, yh))
}

/// Scores `out` against the reference of `rec` under `rules`.
pub(crate) fn score(
    rules: &[FieldRule],
    rec: &Record,
    out: &Output,
    precision: Precision,
) -> Result<Score, String> {
    let bits = unit_bits(precision)?;
    let mut worst = 0.0f64;
    for r in rules {
        let (y, yh) = pair(rec, out, r.field)?;
        if !yh.iter().all(|x| x.is_finite()) {
            return Ok(Score::NonFinite);
        }
        let floor = match r.floor {
            Floor::Unit => vec![Dyadic {
                neg: false,
                mant: 1,
                exp: 0,
            }],
            Floor::Tiny => vec![Dyadic {
                neg: false,
                mant: 1,
                exp: -1022,
            }],
            Floor::Scale { input, skip } => rec
                .input(input)
                .and_then(|v| v.get(skip..))
                .ok_or_else(|| format!("the record has no input `{input}[{skip}..]`"))?
                .iter()
                .map(|&x| Dyadic::of(x))
                .collect(),
        };
        let aligned = match r.sign {
            SignRule::Fixed => false,
            SignRule::Align => true,
            SignRule::AlignAtW0 => rec
                .input("q")
                .and_then(<[f64]>::first)
                .is_some_and(|w| w.to_bits() == 0),
        };
        worst = worst.max(forward_u(&y.data, yh, &floor, aligned, bits));
    }
    Ok(if worst.is_finite() {
        Score::Finite(worst)
    } else {
        Score::NonFinite
    })
}

/// The score of a function scored by backward error only (`Rule::BackwardOnly`): every reference
/// field must be answered with its length, and every value must be finite.
pub(crate) fn finiteness(rec: &Record, out: &Output) -> Result<Score, String> {
    let mut finite = true;
    for field in rec.reference.keys() {
        finite &= pair(rec, out, field)?.1.iter().all(|x| x.is_finite());
    }
    Ok(if finite {
        Score::Unscored
    } else {
        Score::NonFinite
    })
}

/// A shift or power that is non-negative by construction.
fn nonneg(n: i64) -> usize {
    debug_assert!(n >= 0);
    usize::try_from(n).unwrap_or(0)
}

/// `‖ŷ − y‖ / (max(‖y‖, ‖floor‖) · 2^-bits)`, all of `y`, `ŷ` and `floor` exact; `aligned` takes
/// the smaller of `‖ŷ − y‖` and `‖ŷ + y‖`. `0/0` is 0, `x/0` is infinite, and so is an `F` past binary64.
fn forward_u(y: &[Decimal], yh: &[f64], floor: &[Dyadic], aligned: bool, bits: usize) -> f64 {
    let yh: Vec<Dyadic> = yh.iter().map(|&x| Dyadic::of(x)).collect();
    // Times S = 2^p · 5^q every value is an integer: A·10^c is A·5^(c+q)·2^(c+p), m·2^k is
    // m·2^(k+p)·5^q. A zero has no exponent to speak of and is 0.
    let q = y
        .iter()
        .filter(|d| d.mant != 0)
        .map(|d| -i64::from(d.exp10))
        .fold(0, i64::max);
    let low = |v: &[Dyadic]| {
        let exps = v.iter().filter(|d| d.mant != 0).map(|d| -i64::from(d.exp));
        exps.fold(0, i64::max)
    };
    let p = q.max(low(&yh)).max(low(floor));
    let five = |n: i64| BigUint::from(5u32).pow(nonneg(n) as u32);
    let s5 = five(q);
    let signed = |neg, mag| BigInt::from_biguint(if neg { Sign::Minus } else { Sign::Plus }, mag);
    let dyadic = |d: &Dyadic| {
        if d.mant == 0 {
            return BigInt::default();
        }
        signed(
            d.neg,
            (BigUint::from(d.mant) << nonneg(i64::from(d.exp) + p)) * &s5,
        )
    };
    let decimal = |d: &Decimal| {
        if d.mant == 0 {
            return BigInt::default();
        }
        let c = i64::from(d.exp10);
        signed(
            d.neg,
            (BigUint::from(d.mant) * five(c + q)) << nonneg(c + p),
        )
    };
    let sy: Vec<BigInt> = y.iter().map(decimal).collect();
    let square = |b: BigInt| b.magnitude() * b.magnitude();
    let (mut minus, mut plus) = (BigUint::default(), BigUint::default());
    for (h, b) in yh.iter().map(dyadic).zip(&sy) {
        minus += square(&h - b);
        if aligned {
            plus += square(&h + b);
        }
    }
    let err = if aligned { minus.min(plus) } else { minus };
    let sum_y: BigUint = sy.into_iter().map(square).sum();
    let sum_floor: BigUint = floor.iter().map(|d| square(dyadic(d))).sum();
    let scale = sum_y.max(sum_floor);
    if scale.bits() == 0 {
        return if err.bits() == 0 { 0.0 } else { f64::INFINITY };
    }
    sqrt_ratio(&(err << (2 * bits)), &scale)
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::conformance::testkit::record;

    const U: f64 = 1.0 / (1u64 << 53) as f64;

    fn field_score(
        y: &[&str],
        yh: &[f64],
        floor: Floor,
        sign: SignRule,
        inputs: &[(&str, &[f64])],
    ) -> Result<f64, String> {
        let rec = record(inputs, &[("v", y)])?;
        let out = Output::from([("v".to_string(), yh.to_vec())]);
        finite(score(
            &[field("v", floor, sign)],
            &rec,
            &out,
            Precision::F64,
        )?)
    }

    fn finite(score: Score) -> Result<f64, String> {
        match score {
            Score::Finite(f) => Ok(f),
            other => Err(format!("{other:?}")),
        }
    }

    fn plain(y: &[&str], yh: &[f64]) -> Result<f64, String> {
        field_score(y, yh, Floor::Unit, SignRule::Fixed, &[])
    }

    fn near(a: f64, b: f64) {
        assert!((a - b).abs() <= 1e-15 * b.abs(), "{a} vs {b}");
    }

    #[test]
    fn one_ulp_above_one_is_two_u() -> Result<(), String> {
        assert_eq!(plain(&["1e0"], &[1.0 + f64::EPSILON])?, 2.0);
        assert_eq!(plain(&["1e0"], &[1.0 - f64::EPSILON / 2.0])?, 1.0);
        assert_eq!(plain(&["1e0"], &[1.0])?, 0.0);
        Ok(())
    }

    #[test]
    fn the_error_is_norm_wise_over_the_whole_field() -> Result<(), String> {
        // ulp(3) = 2^-51, ulp(4) = 2^-50 (4 sits on a binade edge, so the step up is 2^-50).
        near(
            plain(&["3e0", "4e0"], &[3.0 + 2.0 * f64::EPSILON, 4.0])?,
            0.8,
        );
        let both = plain(
            &["3e0", "4e0"],
            &[3.0 + 2.0 * f64::EPSILON, 4.0 + 4.0 * f64::EPSILON],
        )?;
        near(both, 4.0 * 5.0f64.sqrt() / 5.0);
        Ok(())
    }

    #[test]
    fn the_reference_is_never_rounded_first() -> Result<(), String> {
        // 0.1f64 = 3602879701896397 / 2^55; minus 1/10 is 2^-55 / 5, and (2^-55 / 5) / (0.1 u) = 1/2.
        let relative = |y: &str, yh| field_score(&[y], &[yh], Floor::Tiny, SignRule::Fixed, &[]);
        assert_eq!(relative("1e-1", 0.1)?, 0.5);
        // 5e2 has a positive decimal exponent; ulp(500) = 2^-44: 2^-44 / (500 u) = 1.024.
        near(plain(&["5e2"], &[500.0 + 2.0f64.powi(-44)])?, 1.024);
        // A reference off every double by 1e-30 relative still resolves 1e-14 of u.
        let f = plain(&["1.00000000000000000000000000001e0"], &[1.0])?;
        assert!((f * U - 1e-29).abs() < 1e-40, "{f}");
        Ok(())
    }

    #[test]
    fn the_floor_keeps_a_zero_reference_finite_and_sets_the_scale() -> Result<(), String> {
        let unit = |y: &str, yh: f64| plain(&[y], &[yh]);
        assert_eq!(unit("0e0", 0.0)?, 0.0);
        assert_eq!(unit("0e0", -0.0)?, 0.0);
        assert_eq!(unit("0e0", U)?, 1.0);
        let tiny = |y: &str, yh: f64| field_score(&[y], &[yh], Floor::Tiny, SignRule::Fixed, &[]);
        assert_eq!(tiny("0e0", 0.0)?, 0.0);
        assert_eq!(tiny("0e0", f64::from_bits(1))?, 2.0);
        // A total loss of 1e-20: absolute under floor 1, the full 1/u under the smallest normal.
        near(unit("1e-20", 0.0)?, 1e-20 / U);
        assert_eq!(tiny("1e-20", 0.0)?, 9_007_199_254_740_992.0);
        Ok(())
    }

    #[test]
    fn the_translation_floor_is_the_norm_of_the_inputs_translation_entries() -> Result<(), String> {
        let scale = Floor::Scale {
            input: "tau",
            skip: 1,
        };
        let tau: &[(&str, &[f64])] = &[("tau", &[0.5, 3000.0, 4000.0])];
        // ‖y‖ = 10 sits under the scale ‖ρ‖ = 5000; ulp(10) = 2^-49: 2^-49 / (5000 u) = 16 / 5000.
        let f = field_score(
            &["10e0", "0e0"],
            &[10.0 + 2.0f64.powi(-49), 0.0],
            scale,
            SignRule::Fixed,
            tau,
        )?;
        near(f, 16.0 / 5000.0);
        let scale = Floor::Scale {
            input: "tau",
            skip: 1,
        };
        assert!(field_score(&["1e0"], &[1.0], scale, SignRule::Fixed, &[]).is_err());
        Ok(())
    }

    #[test]
    fn signs_align_for_quaternions_and_at_w_plus_zero_only() -> Result<(), String> {
        let q = ["1e0", "0e0", "0e0", "0e0"];
        let flipped = [-1.0, 0.0, 0.0, 0.0];
        let f = |sign, w: f64| {
            field_score(
                &q,
                &flipped,
                Floor::Unit,
                sign,
                &[("q", &[w, 0.0, 0.0, 0.0])],
            )
        };
        assert_eq!(f(SignRule::Fixed, 1.0)?, 2.0 / U);
        assert_eq!(f(SignRule::Align, 1.0)?, 0.0);
        assert_eq!(f(SignRule::AlignAtW0, 0.0)?, 0.0);
        assert_eq!(f(SignRule::AlignAtW0, -0.0)?, 2.0 / U);
        assert_eq!(f(SignRule::AlignAtW0, 1e-3)?, 2.0 / U);
        // Alignment picks the nearer sign, not always the flip: one ulp off after flipping is 2 u.
        let off = [-(1.0 + f64::EPSILON), 0.0, 0.0, 0.0];
        assert_eq!(
            field_score(&q, &off, Floor::Unit, SignRule::Align, &[])?,
            2.0
        );
        Ok(())
    }

    #[test]
    fn non_finite_outputs_are_not_scored() -> Result<(), String> {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let rec = record(&[], &[("v", &["1e0", "2e0"])])?;
            let out = Output::from([("v".to_string(), vec![1.0, bad])]);
            let s = score(
                &[field("v", Floor::Unit, SignRule::Fixed)],
                &rec,
                &out,
                Precision::F64,
            )?;
            assert_eq!(s, Score::NonFinite);
        }
        Ok(())
    }

    #[test]
    fn subnormals_are_judged_by_their_quantization() -> Result<(), String> {
        let tiny = |y: &str, yh: f64| field_score(&[y], &[yh], Floor::Tiny, SignRule::Fixed, &[]);
        let y = "1e-311";
        let rounded = y.parse::<Decimal>()?.to_f64();
        assert!(rounded > 0.0 && rounded < f64::MIN_POSITIVE);
        assert!(tiny(y, rounded)? <= 1.0);
        let up = tiny(y, rounded.next_up())?;
        assert!((1.0..=3.0).contains(&up), "{up}");
        // 2^-1074 to 30 digits, met by the smallest subnormal: the reference's own rounding only.
        let f = tiny("4.94065645841246544176568792868e-324", f64::from_bits(1))?;
        assert!(f > 0.0 && f < 1e-29, "{f}");
        // Under a floor of 1 the same total loss is invisible, which is why tangents do not use it.
        assert!(plain(&[y], &[0.0])? < 1e-290);
        Ok(())
    }

    #[test]
    fn a_record_scores_its_worst_field() -> Result<(), String> {
        let rec = record(&[], &[("a", &["1e0"]), ("b", &["1e0"])])?;
        let out = Output::from([
            ("a".to_string(), vec![1.0 + f64::EPSILON]),
            ("b".to_string(), vec![1.0]),
        ]);
        let rules = [
            field("a", Floor::Unit, SignRule::Fixed),
            field("b", Floor::Unit, SignRule::Fixed),
        ];
        assert_eq!(
            score(&rules, &rec, &out, Precision::F64)?,
            Score::Finite(2.0)
        );
        let short = Output::from([("a".to_string(), vec![1.0])]);
        assert!(score(&rules, &rec, &short, Precision::F64).is_err());
        let long = Output::from([
            ("a".to_string(), vec![1.0, 2.0]),
            ("b".to_string(), vec![1.0]),
        ]);
        assert!(score(&rules, &rec, &long, Precision::F64).is_err());
        Ok(())
    }

    #[test]
    fn a_backward_only_record_is_checked_for_shape_and_finiteness() -> Result<(), String> {
        let rec = record(&[], &[("q", &["1e0", "0e0"])])?;
        let out = |q: Vec<f64>| Output::from([("q".to_string(), q)]);
        assert_eq!(finiteness(&rec, &out(vec![9.0, 9.0]))?, Score::Unscored);
        assert_eq!(
            finiteness(&rec, &out(vec![f64::NAN, 0.0]))?,
            Score::NonFinite
        );
        assert!(finiteness(&rec, &out(vec![1.0])).is_err());
        assert!(finiteness(&rec, &Output::new()).is_err());
        Ok(())
    }

    #[test]
    fn f32_is_an_error() -> Result<(), String> {
        let rec = record(&[], &[("v", &["1e0"])])?;
        let out = Output::from([("v".to_string(), vec![1.0])]);
        let e = score(
            &[field("v", Floor::Unit, SignRule::Fixed)],
            &rec,
            &out,
            Precision::F32,
        );
        assert!(e.is_err_and(|e| e.contains("f32 is not supported")));
        assert!(unit_bits(Precision::F32).is_err());
        Ok(())
    }

    #[test]
    fn the_table_finds_a_family_through_its_n_suffix() {
        let fields = |id| match rule(id) {
            Some(Rule::Forward(f)) => f.len(),
            _ => 0,
        };
        assert_eq!(
            (fields("sen3_exp_n3"), fields("sen3_exp"), fields("se2_exp")),
            (2, 2, 2)
        );
        assert!(matches!(rule("so3_from_matrix"), Some(Rule::BackwardOnly)));
        assert!(rule("sen3_exp_n4").is_none() && rule("nope").is_none());
    }

    /// The score of `id`'s rule on a record whose every rule field is the one value `y`, answered
    /// by its correct rounding everywhere but `field`, which is answered `yh`.
    fn probe(
        id: &str,
        field: &str,
        y: &str,
        yh: f64,
        inputs: &[(&str, &[f64])],
    ) -> Result<Score, String> {
        let Some(Rule::Forward(rules)) = rule(id) else {
            return Err(format!("`{id}` has no forward rule"));
        };
        let exact = y.parse::<Decimal>()?.to_f64();
        let one = [y];
        let reference: Vec<(&str, &[&str])> = rules.iter().map(|r| (r.field, &one[..])).collect();
        let out: Output = rules
            .iter()
            .map(|r| {
                (
                    r.field.to_string(),
                    vec![if r.field == field { yh } else { exact }],
                )
            })
            .collect();
        score(rules, &record(inputs, &reference)?, &out, Precision::F64)
    }

    /// The floor a row of `TABLE` has, as the tests below state it and `TABLE` does not.
    #[derive(Clone, Copy)]
    enum Want {
        Tiny,
        Unit,
        /// The 2-norm of the translation entries of the input `tau`.
        Scale,
    }

    /// Every scored field of every forward id: `(id, field, floor, quaternion sign alignment)`.
    fn wanted() -> Vec<(String, &'static str, Want, bool)> {
        let mut rows = Vec::new();
        let mut add = |names: &[&str], field, want, aligned| {
            for name in names {
                let ids = if name.starts_with("sen3_") {
                    (1..=3).map(|n| format!("{name}_n{n}")).collect()
                } else {
                    vec![name.to_string()]
                };
                rows.extend(ids.into_iter().map(|id| (id, field, want, aligned)));
            }
        };
        let coeffs = [
            "coeff_k", "coeff_a", "coeff_b", "coeff_c", "coeff_d", "coeff_e", "coeff_r",
        ];
        add(&coeffs, "value", Want::Tiny, false);
        add(&coeffs, "d_branch", Want::Tiny, false);
        add(&["so3_exp", "sen3_exp"], "q", Want::Unit, true);
        add(&["sen3_exp"], "x", Want::Scale, false);
        add(&["so3_log"], "phi", Want::Tiny, false);
        add(&["sen3_log", "se2_log"], "tau", Want::Tiny, false);
        add(&["so3_act"], "Rp", Want::Unit, false);
        let jac = ["jr", "jl", "jr_inv", "jl_inv"];
        for group in ["so3", "sen3", "se2"] {
            let names: Vec<String> = jac.iter().map(|j| format!("{group}_{j}")).collect();
            add(
                &names.iter().map(String::as_str).collect::<Vec<_>>(),
                "J",
                Want::Unit,
                false,
            );
        }
        add(&["sen3_ad", "se2_ad"], "Ad", Want::Unit, false);
        add(&["so2_exp", "se2_exp"], "z", Want::Unit, false);
        add(&["so2_log"], "theta", Want::Tiny, false);
        add(&["se2_exp"], "t", Want::Scale, false);
        rows
    }

    /// An input tangent whose translation entries have norm 5, the whole vector 13, and the
    /// entries after the first translation entry 4: a floor over the wrong slice reads otherwise.
    fn tau_of(id: &str) -> Vec<f64> {
        if id.starts_with("se2_") {
            return vec![12.0, 3.0, 4.0];
        }
        let n = id.chars().last().and_then(|c| c.to_digit(10)).unwrap_or(1) as usize;
        let mut tau = vec![0.0, 0.0, 12.0, 3.0, 4.0, 0.0];
        tau.resize(3 + 3 * n, 0.0);
        tau
    }

    #[test]
    fn the_table_scores_exactly_the_fields_this_test_lists() {
        let wanted = wanted();
        for (id, ..) in &wanted {
            let mut want: Vec<&str> = wanted.iter().filter(|w| w.0 == *id).map(|w| w.1).collect();
            let mut got: Vec<&str> = match rule(id) {
                Some(Rule::Forward(rules)) => rules.iter().map(|r| r.field).collect(),
                _ => Vec::new(),
            };
            want.sort_unstable();
            got.sort_unstable();
            assert_eq!(got, want, "{id}");
        }
        for (family, r) in TABLE {
            let listed =
                |w: &(String, &str, Want, bool)| w.0 == *family || w.0 == format!("{family}_n1");
            assert!(
                matches!(r, BackwardOnly) || wanted.iter().any(listed),
                "{family}"
            );
        }
    }

    #[test]
    fn each_row_has_the_floor_this_test_lists() -> Result<(), String> {
        // A total loss of a 1e-20 reference: 1/u under the smallest normal, 1e-20/u under a
        // floor of 1, 1e-20/(5 u) under the translation scale.
        let unit = 1e-20 * (1u64 << 53) as f64;
        for (id, field, want, _) in wanted() {
            let tau = tau_of(&id);
            let f = finite(probe(&id, field, "1e-20", 0.0, &[("tau", &tau)])?)?;
            match want {
                Want::Tiny => assert_eq!(f, 9_007_199_254_740_992.0, "{id} {field}"),
                Want::Unit => near(f, unit),
                Want::Scale => near(f, unit / 5.0),
            }
        }
        Ok(())
    }

    #[test]
    fn only_quaternions_align_and_so3_log_at_w_plus_zero() -> Result<(), String> {
        // A field answered with its own negative: 0 when aligned, else 2/u over the floor.
        let two_over_u = 2.0 * (1u64 << 53) as f64;
        for (id, field, want, aligned) in wanted() {
            let tau = tau_of(&id);
            let inputs: [(&str, &[f64]); 2] = [("tau", &tau), ("q", &[0.5, 0.5, 0.5, 0.5])];
            let f = finite(probe(&id, field, "1e0", -1.0, &inputs)?)?;
            match (aligned, want) {
                (true, _) => assert_eq!(f, 0.0, "{id} {field}"),
                (false, Want::Scale) => near(f, two_over_u / 5.0),
                (false, _) => assert_eq!(f, two_over_u, "{id} {field}"),
            }
        }
        let at = |id, field, w: f64| {
            finite(probe(
                id,
                field,
                "1e0",
                -1.0,
                &[("q", &[w, 1.0, 0.0, 0.0])],
            )?)
        };
        assert_eq!(at("so3_log", "phi", 0.0)?, 0.0);
        assert_eq!(at("so3_log", "phi", -0.0)?, two_over_u);
        assert_eq!(at("so3_log", "phi", 0.5)?, two_over_u);
        assert_eq!(at("sen3_log_n1", "tau", 0.0)?, two_over_u);
        Ok(())
    }

    #[test]
    fn a_quaternion_aligns_as_a_whole_not_per_component() -> Result<(), String> {
        let Some(Rule::Forward(rules)) = rule("so3_exp") else {
            return Err("so3_exp has no forward rule".into());
        };
        let rec = record(&[], &[("q", &["5e-1", "5e-1", "5e-1", "5e-1"])])?;
        let f = |q: [f64; 4]| {
            let out = Output::from([("q".to_string(), q.to_vec())]);
            score(rules, &rec, &out, Precision::F64)
        };
        assert_eq!(f([-0.5; 4])?, Score::Finite(0.0));
        // One component flipped: 1/u either way, the whole-vector minimum is not per-component.
        assert_eq!(
            f([-0.5, 0.5, 0.5, 0.5])?,
            Score::Finite(9_007_199_254_740_992.0)
        );
        Ok(())
    }

    #[test]
    fn an_error_past_binary64_is_non_finite_not_infinite() -> Result<(), String> {
        // 1 against a zero reference under the smallest normal: 2^1075 in u.
        assert_eq!(probe("so3_log", "phi", "0e0", 1.0, &[])?, Score::NonFinite);
        // A zero reference under a zero translation scale.
        let tau: &[(&str, &[f64])] = &[("tau", &[0.0, 0.0, 0.0])];
        assert_eq!(probe("se2_exp", "t", "0e0", 1e-300, tau)?, Score::NonFinite);
        assert_eq!(probe("se2_exp", "t", "0e0", 0.0, tau)?, Score::Finite(0.0));
        Ok(())
    }
}
