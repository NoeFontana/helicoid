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
//! - **`f32`** (`docs/decisions/0016` item 2) is scored in units of `u = 2^-24` on the `@f32` strata,
//!   whose inputs are exactly binary32 and whose references are the function at those inputs. The
//!   smallest-normal floor is `2^-126`, `EA.4(b)`'s; an output that is not exactly a binary32 is an
//!   error, not a score: a subject that computed in binary64 has not been measured at `f32`.

use helicoid_linalg::Precision;
use num_bigint::{BigInt, BigUint, Sign};

use super::corpus::{exact_f32, Record, Tensor};
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
    /// `eig3` (`0056` decision 3): `lambda` by forward error, `V` column by column, sign-aligned and
    /// weighted by the reference gap, [`eigen`].
    Eigen,
    /// `solve_cubic`: the distance between the valid slots and the reference roots, [`roots`].
    Roots,
    /// `chol`: the reported `mask` must be the reference's; `fields` are scored where it is set,
    /// [`masked`].
    Masked {
        mask: &'static str,
        fields: &'static [FieldRule],
    },
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
            Rule::Eigen => eigen(rec, out, precision),
            Rule::Roots => roots(rec, out, precision),
            Rule::Masked { mask, fields } => masked(mask, fields, rec, out, precision),
        }
    }

    /// The reference fields the rule reads, sorted; `None` for [`Rule::BackwardOnly`], which reads
    /// whatever the reference holds.
    #[cfg(test)]
    pub(crate) fn reference_fields(&self) -> Option<Vec<&'static str>> {
        let mut fields = match self {
            Rule::Forward(fields) => fields.iter().map(|f| f.field).collect(),
            Rule::BackwardOnly => return None,
            Rule::Eigen => vec!["lambda", "V"],
            Rule::Roots => vec!["re", "im"],
            Rule::Masked { mask, fields } => {
                let mut v: Vec<&str> = fields.iter().map(|f| f.field).collect();
                v.push(mask);
                v
            }
        };
        fields.sort_unstable();
        Some(fields)
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
/// The geodesic's pose (`PHASE4.md` §4). Its quaternion is `Align`, as `so3_exp`'s and
/// `sen3_exp`'s are: the output is a **pose**, and `−q` is the same rotation, where `sen3_log`'s
/// tangent is `Fixed` because `−τ` is a different logarithm.
///
/// The translation's floor is `‖x₀‖`, the scale of the problem rather than of the answer: §11
/// names "the ρ-scale" for a translation and this id has no input tangent to read it from. It
/// binds only where the answer cancels against the poses — between two nearby poses the answer is
/// their own size either way — which is the reading `0014` (draft) question 1 leaves open and the
/// one `sen3_exp`'s `Scale` already takes.
const SO3_GEODESIC: &[FieldRule] = &[field("q", Floor::Unit, SignRule::Align)];
const SE3_GEODESIC: &[FieldRule] = &[
    field("q", Floor::Unit, SignRule::Align),
    field(
        "x",
        Floor::Scale {
            input: "x0",
            skip: 0,
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

/// `0056` decision 3: relative, as a tangent's is (`Floor::Tiny`), for the factor, the solution,
/// the projected quaternion and every value and derivative of `Dual`.
const CHOL: &[FieldRule] = &[field("L", Floor::Tiny, SignRule::Fixed)];
const CHOL_SOLVE: &[FieldRule] = &[field("x", Floor::Tiny, SignRule::Fixed)];
const RENORMALIZE: &[FieldRule] = &[field("q", Floor::Tiny, SignRule::Fixed)];
const REAL_UNARY: &[FieldRule] = &[
    field("value", Floor::Tiny, SignRule::Fixed),
    field("d", Floor::Tiny, SignRule::Fixed),
];
const REAL_SIN_COS: &[FieldRule] = &[
    field("sin", Floor::Tiny, SignRule::Fixed),
    field("cos", Floor::Tiny, SignRule::Fixed),
    field("d_sin", Floor::Tiny, SignRule::Fixed),
    field("d_cos", Floor::Tiny, SignRule::Fixed),
];
const REAL_ATAN2: &[FieldRule] = &[
    field("value", Floor::Tiny, SignRule::Fixed),
    field("d_y", Floor::Tiny, SignRule::Fixed),
    field("d_x", Floor::Tiny, SignRule::Fixed),
];
const REAL_DIV: &[FieldRule] = &[
    field("value", Floor::Tiny, SignRule::Fixed),
    field("d_n", Floor::Tiny, SignRule::Fixed),
    field("d_d", Floor::Tiny, SignRule::Fixed),
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
    ("coeff_cos_half", Forward(COEFF)),
    ("coeff_r", Forward(COEFF)),
    ("so3_exp", Forward(SO3_EXP)),
    ("so3_log", Forward(SO3_LOG)),
    ("so3_act", Forward(SO3_ACT)),
    ("so3_from_matrix", BackwardOnly),
    ("so3_jr", Forward(JAC)),
    ("so3_jl", Forward(JAC)),
    ("so3_jr_inv", Forward(JAC)),
    ("so3_jl_inv", Forward(JAC)),
    ("so3_geodesic", Forward(SO3_GEODESIC)),
    ("se3_geodesic", Forward(SE3_GEODESIC)),
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
    ("solve_cubic", Rule::Roots),
    ("eig3", Rule::Eigen),
    (
        "chol",
        Rule::Masked {
            mask: "valid",
            fields: CHOL,
        },
    ),
    ("chol_solve", Forward(CHOL_SOLVE)),
    ("quat_renormalize", Forward(RENORMALIZE)),
    ("real_sqrt", Forward(REAL_UNARY)),
    ("real_cbrt", Forward(REAL_UNARY)),
    ("real_acos", Forward(REAL_UNARY)),
    ("real_sin_cos", Forward(REAL_SIN_COS)),
    ("real_atan2", Forward(REAL_ATAN2)),
    ("real_div", Forward(REAL_DIV)),
];

/// The rule of function id `fn_id`, if the table has one.
pub(crate) fn rule(fn_id: &str) -> Option<&'static Rule> {
    let family = match fn_id.rsplit_once("_n") {
        Some((f, n)) if f.starts_with("sen3_") && matches!(n, "1" | "2" | "3") => f,
        Some((f, n)) if matches!(f, "chol" | "chol_solve") && matches!(n, "3" | "6") => f,
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

/// `log2(1/u)`.
pub(crate) fn unit_bits(precision: Precision) -> usize {
    match precision {
        Precision::F64 => 53,
        Precision::F32 => 24,
    }
}

/// The exponent of the smallest normal number.
fn min_exp(precision: Precision) -> i32 {
    match precision {
        Precision::F64 => -1022,
        Precision::F32 => -126,
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
    let bits = unit_bits(precision);
    let mut worst = 0.0f64;
    for r in rules {
        let (y, yh) = pair(rec, out, r.field)?;
        if !yh.iter().all(|x| x.is_finite()) {
            return Ok(Score::NonFinite);
        }
        if precision == Precision::F32 && yh.iter().any(|&x| exact_f32(x).is_none()) {
            return Err(format!("`{}` is not exactly binary32", r.field));
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
                exp: min_exp(precision),
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

/// Every value of one score over one denominator `S = 2^p · 5^q`: times `S`, a decimal `A·10^c` is
/// `A·5^(c+q)·2^(c+p)` and a dyadic `m·2^k` is `m·2^(k+p)·5^q`, both integers. A zero has no
/// exponent to speak of and is 0.
struct Exact {
    p: i64,
    q: i64,
    s5: BigUint,
}

impl Exact {
    fn of<'a>(
        decimals: impl IntoIterator<Item = &'a Decimal>,
        dyadics: impl IntoIterator<Item = &'a Dyadic>,
    ) -> Self {
        let q = decimals
            .into_iter()
            .filter(|d| d.mant != 0)
            .map(|d| -i64::from(d.exp10))
            .fold(0, i64::max);
        let p = dyadics
            .into_iter()
            .filter(|d| d.mant != 0)
            .map(|d| -i64::from(d.exp))
            .fold(q, i64::max);
        Self { p, q, s5: five(q) }
    }

    fn dyadic(&self, d: &Dyadic) -> BigInt {
        if d.mant == 0 {
            return BigInt::default();
        }
        signed(
            d.neg,
            (BigUint::from(d.mant) << nonneg(i64::from(d.exp) + self.p)) * &self.s5,
        )
    }

    fn decimal(&self, d: &Decimal) -> BigInt {
        if d.mant == 0 {
            return BigInt::default();
        }
        let c = i64::from(d.exp10);
        signed(
            d.neg,
            (BigUint::from(d.mant) * five(c + self.q)) << nonneg(c + self.p),
        )
    }

    /// `S²`.
    fn square(&self) -> BigUint {
        (&self.s5 * &self.s5) << nonneg(2 * self.p)
    }
}

fn five(n: i64) -> BigUint {
    BigUint::from(5u32).pow(nonneg(n) as u32)
}

fn signed(neg: bool, mag: BigUint) -> BigInt {
    BigInt::from_biguint(if neg { Sign::Minus } else { Sign::Plus }, mag)
}

fn square(b: &BigInt) -> BigUint {
    b.magnitude() * b.magnitude()
}

/// `sqrt(num / den)` in `u`, given `num` before its `2^(2·bits)`: `0/0` is 0, `x/0` is infinite.
fn ratio_u(num: &BigUint, den: &BigUint, bits: usize) -> f64 {
    if den.bits() == 0 {
        return if num.bits() == 0 { 0.0 } else { f64::INFINITY };
    }
    sqrt_ratio(&(num << (2 * bits)), den)
}

/// `‖ŷ − y‖ / (max(‖y‖, ‖floor‖) · 2^-bits)`, all of `y`, `ŷ` and `floor` exact; `aligned` takes
/// the smaller of `‖ŷ − y‖` and `‖ŷ + y‖`. `0/0` is 0, `x/0` is infinite, and so is an `F` past binary64.
fn forward_u(y: &[Decimal], yh: &[f64], floor: &[Dyadic], aligned: bool, bits: usize) -> f64 {
    let yh: Vec<Dyadic> = yh.iter().map(|&x| Dyadic::of(x)).collect();
    let ex = Exact::of(y, yh.iter().chain(floor));
    let sy: Vec<BigInt> = y.iter().map(|d| ex.decimal(d)).collect();
    let (mut minus, mut plus) = (BigUint::default(), BigUint::default());
    for (h, b) in yh.iter().map(|d| ex.dyadic(d)).zip(&sy) {
        minus += square(&(&h - b));
        if aligned {
            plus += square(&(&h + b));
        }
    }
    let err = if aligned { minus.min(plus) } else { minus };
    let sum_y: BigUint = sy.iter().map(square).sum();
    let sum_floor: BigUint = floor.iter().map(|d| square(&ex.dyadic(d))).sum();
    ratio_u(&err, &sum_y.max(sum_floor), bits)
}

/// The subject's values of `field`, `n` of them, finite and of the precision: `Ok(None)` when one
/// is not finite.
fn values<'a>(
    out: &'a Output,
    field: &str,
    n: usize,
    precision: Precision,
) -> Result<Option<&'a [f64]>, String> {
    let v = out
        .get(field)
        .ok_or_else(|| format!("the subject returned no `{field}`"))?;
    if v.len() != n {
        return Err(format!(
            "`{field}`: the subject returned {} values, not {n}",
            v.len()
        ));
    }
    if !v.iter().all(|x| x.is_finite()) {
        return Ok(None);
    }
    if precision == Precision::F32 && v.iter().any(|&x| exact_f32(x).is_none()) {
        return Err(format!("`{field}` is not exactly binary32"));
    }
    Ok(Some(v))
}

/// The reference values of `field`, `n` of them.
fn reference<'a>(rec: &'a Record, field: &str, n: usize) -> Result<&'a [Decimal], String> {
    let y = rec
        .reference
        .get(field)
        .ok_or_else(|| format!("the reference has no `{field}`"))?;
    if y.data.len() != n {
        return Err(format!(
            "the reference's `{field}` has {} values, not {n}",
            y.data.len()
        ));
    }
    Ok(&y.data)
}

/// A mask the subject reports as numbers: 1 set, 0 clear, anything else an error.
fn mask(values: &[f64], field: &str) -> Result<Vec<bool>, String> {
    values
        .iter()
        .map(|&v| match v {
            1.0 => Ok(true),
            0.0 => Ok(false),
            _ => Err(format!("`{field}` holds {v}, not a mask")),
        })
        .collect()
}

fn finite_or_nonfinite(worst: f64) -> Score {
    if worst.is_finite() {
        Score::Finite(worst)
    } else {
        Score::NonFinite
    }
}

/// `NUMERICS.md` §11's eigenvector error (`0056` decision 3), with `lambda`'s forward error beside
/// it: column `i` of `V` scores `min(‖v̂ᵢ − vᵢ‖, ‖v̂ᵢ + vᵢ‖) · gapᵢ / (‖λ‖₂ u)`, `gapᵢ` the
/// reference's `min_{j≠i} |λᵢ − λⱼ|`, so a repeated eigenvalue's column weighs 0.
fn eigen(rec: &Record, out: &Output, precision: Precision) -> Result<Score, String> {
    let bits = unit_bits(precision);
    let (lambda, v) = (reference(rec, "lambda", 3)?, reference(rec, "V", 9)?);
    let (Some(lh), Some(vh)) = (
        values(out, "lambda", 3, precision)?,
        values(out, "V", 9, precision)?,
    ) else {
        return Ok(Score::NonFinite);
    };
    let tiny = [Dyadic {
        neg: false,
        mant: 1,
        exp: min_exp(precision),
    }];
    let mut worst = forward_u(lambda, lh, &tiny, false, bits);
    let vh: Vec<Dyadic> = vh.iter().map(|&x| Dyadic::of(x)).collect();
    let ex = Exact::of(lambda.iter().chain(v), &vh);
    let l: Vec<BigInt> = lambda.iter().map(|d| ex.decimal(d)).collect();
    let norm: BigUint = l.iter().map(square).sum();
    for i in 0..3 {
        let gap = (0..3)
            .filter(|&j| j != i)
            .map(|j| square(&(&l[i] - &l[j])))
            .min()
            .unwrap_or_default();
        let (mut minus, mut plus) = (BigUint::default(), BigUint::default());
        for r in 0..3 {
            let (h, y) = (ex.dyadic(&vh[3 * i + r]), ex.decimal(&v[3 * i + r]));
            minus += square(&(&h - &y));
            plus += square(&(&h + &y));
        }
        // (d S)² (gap S)² / ((‖λ‖ S)² S²): every factor scaled by the one `S`.
        let num = minus.min(plus) * gap;
        worst = worst.max(ratio_u(&num, &(&norm * ex.square()), bits));
    }
    Ok(finite_or_nonfinite(worst))
}

/// `NUMERICS.md` §11's root-set distance (`0056` decision 3): the Hausdorff distance between the
/// valid slots and the three reference roots `re + i·im`, each real root to its nearest valid slot
/// and each valid slot to its nearest root, over `max(‖Z‖₂, smallest normal) · u`. A real root no
/// valid slot answers counts `‖Z‖₂`, so it reads `1/u`: finite, a score and not a failed run.
fn roots(rec: &Record, out: &Output, precision: Precision) -> Result<Score, String> {
    let bits = unit_bits(precision);
    let (re, im) = (reference(rec, "re", 3)?, reference(rec, "im", 3)?);
    let valid = mask(
        out.get("valid").ok_or("the subject returned no `valid`")?,
        "valid",
    )?;
    let slots = out.get("roots").ok_or("the subject returned no `roots`")?;
    if valid.len() != 3 || slots.len() != 3 {
        return Err("`roots` and `valid` hold three slots each".into());
    }
    let answered: Vec<f64> = slots
        .iter()
        .zip(&valid)
        .filter_map(|(&r, &ok)| ok.then_some(r))
        .collect();
    if !answered.iter().all(|x| x.is_finite()) {
        return Ok(Score::NonFinite);
    }
    if precision == Precision::F32 && answered.iter().any(|&x| exact_f32(x).is_none()) {
        return Err("`roots` is not exactly binary32".into());
    }
    let rh: Vec<Dyadic> = answered.iter().map(|&x| Dyadic::of(x)).collect();
    let tiny = Dyadic {
        neg: false,
        mant: 1,
        exp: min_exp(precision),
    };
    let ex = Exact::of(re.iter().chain(im), rh.iter().chain([&tiny]));
    let z: Vec<(BigInt, BigInt)> = re
        .iter()
        .zip(im)
        .map(|(a, b)| (ex.decimal(a), ex.decimal(b)))
        .collect();
    let rh: Vec<BigInt> = rh.iter().map(|d| ex.dyadic(d)).collect();
    let norm: BigUint = z.iter().map(|(a, b)| square(a) + square(b)).sum();
    let dist = |r: &BigInt, (a, b): &(BigInt, BigInt)| square(&(r - a)) + square(b);
    let mut far = BigUint::default();
    for zi in z.iter().filter(|(_, b)| b.sign() == Sign::NoSign) {
        let near = rh.iter().map(|r| dist(r, zi)).min();
        far = far.max(near.unwrap_or_else(|| norm.clone()));
    }
    for r in &rh {
        far = far.max(z.iter().map(|zi| dist(r, zi)).min().unwrap_or_default());
    }
    let floor = square(&ex.dyadic(&tiny));
    Ok(finite_or_nonfinite(ratio_u(&far, &norm.max(floor), bits)))
}

/// `NUMERICS.md` §11's mask rule (`0056` decision 3): a reported mask that is not the reference's
/// reads `1/u`; where both are clear the values beside it are not scored, where both are set they
/// are, under `fields`.
fn masked(
    field: &str,
    fields: &[FieldRule],
    rec: &Record,
    out: &Output,
    precision: Precision,
) -> Result<Score, String> {
    let want = reference(rec, field, 1)?[0].mant != 0;
    let got = mask(
        out.get(field)
            .ok_or_else(|| format!("the subject returned no `{field}`"))?,
        field,
    )?;
    let &[got] = got.as_slice() else {
        return Err(format!("`{field}` holds one mask"));
    };
    if got != want {
        return Ok(Score::Finite(2f64.powi(unit_bits(precision) as i32)));
    }
    if !want {
        return Ok(Score::Finite(0.0));
    }
    score(fields, rec, out, precision)
}

/// A shift or power that is non-negative by construction.
fn nonneg(n: i64) -> usize {
    debug_assert!(n >= 0);
    usize::try_from(n).unwrap_or(0)
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

    /// `field_score` at `f32`, with the floor of a coefficient.
    fn f32_score(y: &str, yh: f64) -> Result<Score, String> {
        let rec = record(&[], &[("v", &[y])])?;
        let out = Output::from([("v".to_string(), vec![yh])]);
        let rule = [field("v", Floor::Tiny, SignRule::Fixed)];
        score(&rule, &rec, &out, Precision::F32)
    }

    #[test]
    fn f32_is_scored_in_units_of_two_to_the_minus_24() -> Result<(), String> {
        let bits = (unit_bits(Precision::F64), unit_bits(Precision::F32));
        assert_eq!(bits, (53, 24));
        // One binary32 ulp above 1 is 2^-23 = 2 u.
        let up = f64::from(1.0f32 + f32::EPSILON);
        assert_eq!(f32_score("1e0", up)?, Score::Finite(2.0));
        assert_eq!(f32_score("1e0", 1.0)?, Score::Finite(0.0));
        // 0.1f32 is 0.100000001490116119384765625: its own reference scores 0, and against 1e-1
        // (never rounded first) the difference 2^-26/10 is a quarter of 0.1 u.
        let tenth = "1.00000001490116119384765625e-1";
        assert_eq!(f32_score(tenth, f64::from(0.1f32))?, Score::Finite(0.0));
        assert_eq!(f32_score("1e-1", f64::from(0.1f32))?, Score::Finite(0.25));
        Ok(())
    }

    /// Every value of every field is a binary32, not the first of each.
    #[test]
    fn every_value_of_every_output_field_must_be_a_binary32_at_f32() -> Result<(), String> {
        let rec = record(&[], &[("a", &["1e0", "2e0"]), ("b", &["1e0"])])?;
        let rules = [
            field("a", Floor::Tiny, SignRule::Fixed),
            field("b", Floor::Tiny, SignRule::Fixed),
        ];
        let out = |a: [f64; 2], b: f64| {
            Output::from([("a".to_string(), a.to_vec()), ("b".to_string(), vec![b])])
        };
        let at = |o: Output| score(&rules, &rec, &o, Precision::F32);
        assert_eq!(at(out([1.0, 2.0], 1.0))?, Score::Finite(0.0));
        let second_value = out([1.0, 2.0 + 2.0 * f64::EPSILON], 1.0);
        let second_field = out([1.0, 2.0], 1.0 + f64::EPSILON);
        for (bad, field) in [(second_value, "`a`"), (second_field, "`b`")] {
            let e = at(bad).err().unwrap_or_default();
            assert!(
                e.starts_with(field) && e.contains("not exactly binary32"),
                "{e}"
            );
        }
        Ok(())
    }

    #[test]
    fn an_f32_output_is_a_binary32_and_the_floor_is_the_smallest_normal_of_binary32(
    ) -> Result<(), String> {
        let e = f32_score("1e0", 1.0 + f64::EPSILON)
            .err()
            .unwrap_or_default();
        assert!(e.contains("not exactly binary32"), "{e}");
        assert_eq!(f32_score("1e0", f64::NAN)?, Score::NonFinite);
        // The smallest binary32 subnormal against 0 is 2^-149 over 2^-126 u: 2 (a half-spacing
        // rounding reads 1); the floor 2^-1022 of binary64 would read 2^897.
        let tiny = f64::from(f32::from_bits(1));
        assert_eq!(f32_score("0e0", tiny)?, Score::Finite(2.0));
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
        /// The 2-norm of the id's scale input: the translation entries of `tau`, or the whole of
        /// `x0` for the geodesic. [`SCALE_X0`] and [`tau_of`] both make that norm 5, so one
        /// expectation covers either reading.
        Scale,
    }

    /// An `x0` of norm 5, the geodesic's scale input (`Floor::Scale { input: "x0", skip: 0 }`).
    const SCALE_X0: [f64; 3] = [3.0, 4.0, 0.0];

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
            "coeff_k",
            "coeff_a",
            "coeff_b",
            "coeff_c",
            "coeff_d",
            "coeff_e",
            "coeff_cos_half",
            "coeff_r",
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
        add(&["so3_geodesic", "se3_geodesic"], "q", Want::Unit, true);
        add(&["se3_geodesic"], "x", Want::Scale, false);
        add(&["chol_solve_n3", "chol_solve_n6"], "x", Want::Tiny, false);
        add(&["quat_renormalize"], "q", Want::Tiny, false);
        for name in ["real_sqrt", "real_cbrt", "real_acos"] {
            add(&[name], "value", Want::Tiny, false);
            add(&[name], "d", Want::Tiny, false);
        }
        for f in ["sin", "cos", "d_sin", "d_cos"] {
            add(&["real_sin_cos"], f, Want::Tiny, false);
        }
        for f in ["value", "d_y", "d_x"] {
            add(&["real_atan2"], f, Want::Tiny, false);
        }
        for f in ["value", "d_n", "d_d"] {
            add(&["real_div"], f, Want::Tiny, false);
        }
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
            let listed = |w: &(String, &str, Want, bool)| {
                w.0 == *family || ["1", "3"].iter().any(|n| w.0 == format!("{family}_n{n}"))
            };
            let own_test = matches!(
                r,
                BackwardOnly | Rule::Eigen | Rule::Roots | Rule::Masked { .. }
            );
            assert!(own_test || wanted.iter().any(listed), "{family}");
        }
    }

    #[test]
    fn each_row_has_the_floor_this_test_lists() -> Result<(), String> {
        // A total loss of a 1e-20 reference: 1/u under the smallest normal, 1e-20/u under a
        // floor of 1, 1e-20/(5 u) under the translation scale.
        let unit = 1e-20 * (1u64 << 53) as f64;
        for (id, field, want, _) in wanted() {
            let tau = tau_of(&id);
            let inputs: [(&str, &[f64]); 2] = [("tau", &tau), ("x0", &SCALE_X0)];
            let f = finite(probe(&id, field, "1e-20", 0.0, &inputs)?)?;
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
            let inputs: [(&str, &[f64]); 3] = [
                ("tau", &tau),
                ("q", &[0.5, 0.5, 0.5, 0.5]),
                ("x0", &SCALE_X0),
            ];
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

    /// `rule`'s score of `out` against a record of `reference`.
    fn rule_score(
        r: &Rule,
        reference: &[(&str, &[&str])],
        out: &[(&str, &[f64])],
    ) -> Result<Score, String> {
        let out: Output = out
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_vec()))
            .collect();
        r.score(&record(&[], reference)?, &out, Precision::F64)
    }

    const LAMBDA: [&str; 3] = ["1e0", "2e0", "4e0"];
    const IDENTITY: [&str; 9] = [
        "1e0", "0e0", "0e0", "0e0", "1e0", "0e0", "0e0", "0e0", "1e0",
    ];
    const EYE: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

    fn eig(lambda: [&str; 3], v: [f64; 9]) -> Result<Score, String> {
        let l: Vec<f64> = lambda
            .iter()
            .map(|s| s.parse::<Decimal>().map(|d| d.to_f64()))
            .collect::<Result<_, _>>()?;
        rule_score(
            &Rule::Eigen,
            &[("lambda", &lambda), ("V", &IDENTITY)],
            &[("lambda", &l), ("V", &v)],
        )
    }

    #[test]
    fn an_eigenvector_scores_its_angle_times_its_gap_and_each_column_aligns_alone(
    ) -> Result<(), String> {
        assert_eq!(eig(LAMBDA, EYE)?, Score::Finite(0.0));
        let mut flipped = EYE;
        flipped[4] = -1.0;
        assert_eq!(eig(LAMBDA, flipped)?, Score::Finite(0.0));
        // Column 0 tilted by e = 2^-40 into axis 1: ‖v̂ − v‖ = e (to first order the length change
        // is e²/2, below the score's resolution here), gap₀ = 1, ‖λ‖ = √21.
        let e = 2f64.powi(-40);
        let mut tilted = EYE;
        tilted[1] = e;
        let want = e / (21f64.sqrt() * U);
        let Score::Finite(f) = eig(LAMBDA, tilted)? else {
            return Err("non-finite".into());
        };
        near(f, want);
        // The same tilt between two equal eigenvalues weighs nothing: their vectors are not unique.
        let Score::Finite(f) = eig(["1e0", "1e0", "4e0"], tilted)? else {
            return Err("non-finite".into());
        };
        assert_eq!(f, 0.0);
        Ok(())
    }

    fn cubic(
        re: [&str; 3],
        im: [&str; 3],
        roots: [f64; 3],
        valid: [f64; 3],
    ) -> Result<Score, String> {
        rule_score(
            &Rule::Roots,
            &[("re", &re), ("im", &im)],
            &[("roots", &roots), ("valid", &valid)],
        )
    }

    #[test]
    fn the_root_set_distance_needs_no_order_and_reads_a_lost_root_as_one_over_u(
    ) -> Result<(), String> {
        let (re, im) = (["1e0", "2e0", "2e0"], ["0e0"; 3]);
        assert_eq!(
            cubic(re, im, [2.0, 1.0, 2.0], [1.0; 3])?,
            Score::Finite(0.0)
        );
        // A double root reported once is still every real root answered.
        assert_eq!(
            cubic(re, im, [2.0, 1.0, 7.0], [1.0, 1.0, 0.0])?,
            Score::Finite(0.0)
        );
        // The root 1 answered only by the slots at 2: 1 / (‖Z‖ u), ‖Z‖ = 3.
        let Score::Finite(f) = cubic(re, im, [2.0, 2.0, 2.0], [1.0; 3])? else {
            return Err("non-finite".into());
        };
        near(f, 1.0 / (3.0 * U));
        // No valid slot at all: every real root is unanswered, ‖Z‖ / (‖Z‖ u).
        assert_eq!(
            cubic(re, im, [2.0, 2.0, 2.0], [0.0; 3])?,
            Score::Finite(1.0 / U)
        );
        // A slot off by one ulp of 2: 2^-51 / (3 u) = 4/3.
        let Score::Finite(f) = cubic(re, im, [1.0, 2.0 + 4.0 * f64::EPSILON / 2.0, 2.0], [1.0; 3])?
        else {
            return Err("non-finite".into());
        };
        near(f, 4.0 / 3.0);
        Ok(())
    }

    #[test]
    fn a_complex_pair_reported_as_a_real_double_root_costs_its_imaginary_part() -> Result<(), String>
    {
        // Roots -1 and 1 ± 2^-30 i; the slots 1, 1 stand for the pair. ‖Z‖² = 3 + 2·2^-60.
        let (re, im) = (
            ["-1e0", "1e0", "1e0"],
            [
                "0e0",
                "-9.31322574615478515625e-10",
                "9.31322574615478515625e-10",
            ],
        );
        let Score::Finite(f) = cubic(re, im, [-1.0, 1.0, 1.0], [1.0; 3])? else {
            return Err("non-finite".into());
        };
        near(f, 2f64.powi(-30) / (3f64.sqrt() * U));
        // Clear slots are not read, whatever they hold.
        assert_eq!(
            cubic(re, im, [-1.0, f64::NAN, 9.0], [1.0, 0.0, 0.0])?,
            Score::Finite(0.0)
        );
        assert!(cubic(re, im, [-1.0, 0.0, 0.0], [1.0, 0.5, 0.0]).is_err());
        Ok(())
    }

    #[test]
    fn a_mask_that_disagrees_reads_one_over_u_and_a_clear_one_hides_its_values(
    ) -> Result<(), String> {
        let rule = Rule::Masked {
            mask: "valid",
            fields: CHOL,
        };
        let at = |want: &str, got: f64, l: f64| {
            rule_score(
                &rule,
                &[("valid", &[want]), ("L", &["2e0"])],
                &[("valid", &[got]), ("L", &[l])],
            )
        };
        assert_eq!(at("1e0", 1.0, 2.0)?, Score::Finite(0.0));
        assert_eq!(
            at("1e0", 1.0, 2.0 + 2.0 * f64::EPSILON)?,
            Score::Finite(2.0)
        );
        assert_eq!(at("1e0", 0.0, 2.0)?, Score::Finite(1.0 / U));
        assert_eq!(at("0e0", 1.0, 2.0)?, Score::Finite(1.0 / U));
        assert_eq!(at("0e0", 0.0, 99.0)?, Score::Finite(0.0));
        Ok(())
    }

    #[test]
    fn the_new_rules_read_the_fields_they_name() {
        let fields = |id| rule(id).and_then(Rule::reference_fields);
        assert_eq!(fields("eig3"), Some(vec!["V", "lambda"]));
        assert_eq!(fields("solve_cubic"), Some(vec!["im", "re"]));
        assert_eq!(fields("chol_n6"), Some(vec!["L", "valid"]));
        assert_eq!(fields("chol_solve_n3"), Some(vec!["x"]));
        assert!(rule("chol_n4").is_none() && fields("so3_from_matrix").is_none());
    }
}
