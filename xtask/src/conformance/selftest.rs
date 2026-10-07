//! `cargo xtask conformance --self-test` (`docs/PHASE1.md` §10): the instrument shown to detect what
//! it claims to. It runs the correct seeded kernel, which no mechanism may fire on, and every
//! planted defect (`crate::seeded`), which must fire its named mechanism, over the `coeff_*` corpus
//! ids (the `Log` defects over `so3_log`: `selftest_so3`; the SE(3) defects over `sen3_*`:
//! `selftest_se3`; the envelope half of the planted `c`: `selftest_envelope`), and fails when either
//! does not hold.
//! Nothing is written to `conformance/results/`.
//!
//! Readings where §10 is silent, each the smallest:
//!
//! - **The curve** is the value field of `coeff_b` (`θ^-2`; through `Dual` the exact arm's
//!   derivative errs by `θ^-4`, so a record's maximum over both fields would fit `p = 4`); §10
//!   says "error curve" and 0014 (draft) asks which field (open question 7). It is
//!   the per-stratum `max_u` against `θ = 10^-k`, the stratum label, for `k = 8, ..., 2`, by least
//!   squares of `ln max_u` on `ln θ`. The window is inclusive and the fit quality is reported, not
//!   gated. The stratum `theta:1e-8` reaches below `sqrt(6u)`, where `b` evaluates to 0 and the
//!   error saturates at `1/u`; that biases `p` low (`docs/maths/coefficients.md`, "Checked
//!   (CO.5-CO.7)").
//! - **At `f32`** (`docs/decisions/0016` item 2) the same mechanisms run on the `@f32` strata of the
//!   coefficient ids, in units of `2^-24`, by the same subject through the same adapter: the curve
//!   is `theta:1e-3@f32` to `theta:1e-1@f32` (`DECADES_F32`; a reading, 0014 (draft) question 26),
//!   in the same window, the non-finite count `theta:exact0@f32`. The correct kernel there runs the
//!   `f32` constants of `generated.rs` (`0016` item 3), whose errors are printed and equal the
//!   CSV's `f32` rows. The sweep's ranking of the planted `c` is binary64's, so `c` is not run at
//!   `f32` (`Registered::no_f32` refuses it); the `Log` and SE(3) defects have no `@f32` stratum to
//!   run on.
//! - **Non-finite** counts the (record, field) pairs of `theta:exact0` in every coefficient id, so a
//!   NaN in `value` or in `d_branch` both count.
//! - **The sweep's ranking** is for `c` with switch `1e-8` (as `θ`: 0014 (draft) question 14) and
//!   two terms (`seeded::C_PLANTED`): the sweep of `coeff_c` (`crate::thresholds`) scores that
//!   candidate like any it sweeps, and it is detected when its objective exceeds the chosen
//!   candidate's by more than [`DOMINATED_BY`], a factor of `10^6` (question 15), "strictly worse
//!   by a stated margin". The reading printed is the objective, the chosen one and the candidate's
//!   rank among the swept ones; §10's other half, the envelope, is `selftest_envelope`. Every
//!   subject is measured with its own `c`, so the correct kernel, whose `c` is the chosen
//!   candidate, must be silent. The rank is of the candidate a subject declares, so the
//!   run also fails when the objective the harness measures for the subject's `c` is not the ranked
//!   one, to the bit: a subject cannot change what it computes and keep its declaration.
//! - **The correct kernel** runs the generated switches (`seeded::generated`, from
//!   `conformance/sweeps/thresholds-seeded.csv`), not an optimum in general: the sweep's choice over
//!   this corpus and grid. Its errors are printed, and equal the CSV's per-field maxima.
//! - **Only the named mechanism is gated** for a defect. The other mechanism's reading is printed
//!   and not gated: `b` by its definition is `0/0` at `z = 0`, so it fires `nonfinite` too
//!   (`the_b_defect_also_fires_the_nonfinite_mechanism_and_the_k_defect_not_the_curve`).

use std::fmt::Write;
use std::path::Path;

use helicoid_linalg::Precision;

use super::metric::{COEFF_D_BRANCH, COEFF_VALUE};
use super::report::Row;
use super::subject::Registered;
use super::{corpus, evaluate_by, selftest_envelope, selftest_linalg, selftest_se3, selftest_so3};
use crate::seeded::{Coeff, Defect, Seeded};
use crate::thresholds::{Ranker, Ranking};

/// The exponent `p` of `θ^-p` that the `b` defect's curve must fit, inclusive (§10).
const WINDOW: (f64, f64) = (1.8, 2.2);
/// The strata `theta:1e-8` to `theta:1e-2` (§10).
const DECADES: std::ops::RangeInclusive<i32> = 2..=8;
/// At `f32`, `theta:1e-3@f32` to `theta:1e-1@f32`: the strata wholly above `√(6u)`, `6.0e-4`,
/// below which `b` by its definition is 0 and its error a plateau of `1/u`. `DECADES` keeps the
/// one stratum that holds `2.6e-8`, its plateau one point of seven (`p = 1.937`); here it would be
/// one of four and fit `p = 1.64`, outside §10's window, which does not depend on `u`. The top
/// stratum, `[0.1, 1)`, is above D12's switch, where the defect and the D12 kernel this reading was
/// made against ran one exact arm; the `f32` sweep leaves the range as it is, and the two strata
/// below it fit `p = 2.188` (0014 (draft) question 26).
const DECADES_F32: std::ops::RangeInclusive<i32> = 1..=3;
/// A candidate is dominated when its objective is more than this factor above the chosen one's.
const DOMINATED_BY: f64 = 1e6;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mechanism {
    ErrorCurve,
    NonfiniteAtZero,
    SweepRank,
}

impl Mechanism {
    const ALL: [Mechanism; 3] = [
        Mechanism::ErrorCurve,
        Mechanism::NonfiniteAtZero,
        Mechanism::SweepRank,
    ];

    /// The mechanisms meaningful at `precision`: the sweep that ranks `c` is binary64's, and no
    /// planted `c` is a candidate of the `f32` sweep.
    fn at(precision: Precision) -> &'static [Mechanism] {
        match precision {
            Precision::F64 => &Self::ALL,
            Precision::F32 => &Self::ALL[..2],
        }
    }

    fn label(self) -> &'static str {
        match self {
            Mechanism::ErrorCurve => "b value curve fits theta^-p",
            Mechanism::NonfiniteAtZero => "nonfinite in theta:exact0",
            Mechanism::SweepRank => "sweep ranks c dominated",
        }
    }

    /// The mechanism that must detect the coefficient `defect` (§10); the `Log` defects are
    /// judged by `selftest_so3` and the SE(3) ones by `selftest_se3`, and have none here.
    fn of(defect: Defect) -> Option<Self> {
        match defect {
            Defect::BNoSeries => Some(Mechanism::ErrorCurve),
            Defect::KSqrtUnsafe => Some(Mechanism::NonfiniteAtZero),
            Defect::CTwoTermsEarly => Some(Mechanism::SweepRank),
            Defect::LogAcos
            | Defect::LogNoFlip
            | Defect::Se3ExpTranslationFirst
            | Defect::QMinusHalf => None,
        }
    }
}

/// A subject, the mechanism it must fire (`None` for the correct kernel, on which none may) and
/// the `(terms, switch)` it runs `c` with, which the sweep ranks.
struct Case {
    subject: Registered,
    must_fire: Option<Mechanism>,
    c: (usize, f64),
}

impl Case {
    fn new(seeded: Seeded, must_fire: Option<Mechanism>) -> Self {
        let c = seeded.candidate(Coeff::C);
        Self {
            subject: seeded.registered(),
            must_fire,
            c: (c.terms, c.switch_z.v),
        }
    }
}

fn cases() -> Vec<Case> {
    let mut all = vec![Case::new(Seeded::generated(), None)];
    for d in Defect::COEFFICIENT {
        all.push(Case::new(Seeded::planted(d), Mechanism::of(d)));
    }
    all
}

/// The cases whose mechanism is meaningful at `f32`: the correct kernel and the two defects the
/// error curve and the non-finite count detect. `c`'s is the sweep's ranking, binary64's for now.
fn cases_f32() -> Vec<Case> {
    let mut all = cases();
    all.retain(|c| c.must_fire != Some(Mechanism::SweepRank));
    all
}

/// The decades of the curve's strata, `theta:1e-k`.
fn decades(precision: Precision) -> std::ops::RangeInclusive<i32> {
    match precision {
        Precision::F64 => DECADES,
        Precision::F32 => DECADES_F32,
    }
}

/// The suffix of the strata a precision is scored on.
fn suffix(precision: Precision) -> &'static str {
    match precision {
        Precision::F64 => "",
        Precision::F32 => corpus::F32_SUFFIX,
    }
}

/// One subject's rows over the coefficient ids, one output field at a time, and where the sweep
/// ranks its `c`.
struct Rows {
    value: Vec<Row>,
    d_branch: Vec<Row>,
    /// Only at binary64 (the sweep).
    ranking: Option<Ranking>,
}

/// Whether a candidate's objective is beyond [`DOMINATED_BY`] times the chosen one's.
fn dominated(r: &Ranking) -> bool {
    r.objective > DOMINATED_BY * r.chosen
}

/// The objective the harness measures for the subject's `c`: the larger of the two fields' maxima
/// over the `theta:*` strata of `coeff_c`, `inf` where an output is not finite, as the sweep scores
/// a candidate; `None` without such a row.
fn measured_c(rows: &Rows) -> Option<f64> {
    let of_c = |r: &&Row| r.fn_id == "coeff_c" && r.stratum.starts_with("theta:");
    let mut seen = rows
        .value
        .iter()
        .chain(&rows.d_branch)
        .filter(of_c)
        .peekable();
    seen.peek()?;
    Some(seen.fold(0.0, |worst, r| match r.nonfinite {
        0 => f64::max(worst, r.max_u),
        _ => f64::INFINITY,
    }))
}

/// `(p, r²)` of the least-squares line `ln y = c - p ln x`; `None` unless there are two points, all
/// finite and positive.
fn power_law(points: &[(f64, f64)]) -> Option<(f64, f64)> {
    let usable = |&(x, y): &(f64, f64)| x.is_finite() && y.is_finite() && x > 0.0 && y > 0.0;
    if points.len() < 2 || !points.iter().all(usable) {
        return None;
    }
    let logs: Vec<(f64, f64)> = points
        .iter()
        .map(|&(x, y)| (libm::log(x), libm::log(y)))
        .collect();
    let n = logs.len() as f64;
    let (mx, my) = (
        logs.iter().map(|p| p.0).sum::<f64>() / n,
        logs.iter().map(|p| p.1).sum::<f64>() / n,
    );
    let sum = |f: &dyn Fn(&(f64, f64)) -> f64| logs.iter().map(f).sum::<f64>();
    let sxx = sum(&|p| (p.0 - mx) * (p.0 - mx));
    let sxy = sum(&|p| (p.0 - mx) * (p.1 - my));
    let syy = sum(&|p| (p.1 - my) * (p.1 - my));
    let r2 = sxy * sxy / (sxx * syy);
    Some((-sxy / sxx, r2))
}

/// The `max_u` of `coeff_b`'s value over the strata of `DECADES` (`DECADES_F32`), as
/// `(θ, max_u)`, from the lowest θ.
fn curve(rows: &Rows, precision: Precision) -> Result<Vec<(f64, f64)>, String> {
    let mut points = Vec::new();
    for k in decades(precision).rev() {
        let stratum = format!("theta:1e-{k}{}", suffix(precision));
        let row = rows
            .value
            .iter()
            .find(|r| r.fn_id == "coeff_b" && r.stratum == stratum)
            .ok_or_else(|| format!("no `coeff_b` row for `{stratum}`"))?;
        let theta = format!("1e-{k}")
            .parse::<f64>()
            .map_err(|e| e.to_string())?;
        points.push((theta, row.max_u));
    }
    Ok(points)
}

/// Whether `mechanism` fires on `rows`, and what it saw.
fn observe(
    mechanism: Mechanism,
    rows: &Rows,
    precision: Precision,
    window: (f64, f64),
) -> Result<(bool, String), String> {
    match mechanism {
        Mechanism::ErrorCurve => Ok(match power_law(&curve(rows, precision)?) {
            Some((p, r2)) => (
                (window.0..=window.1).contains(&p),
                format!("p = {p:.3}, r2 = {r2:.3}"),
            ),
            None => (false, "no curve: a zero or non-finite max".to_string()),
        }),
        Mechanism::NonfiniteAtZero => {
            let fields = [("value", &rows.value), ("d_branch", &rows.d_branch)];
            let exact0 = format!("theta:exact0{}", suffix(precision));
            let hits: Vec<String> = fields
                .iter()
                .flat_map(|&(field, rows)| {
                    let zero = rows.iter().filter(|r| r.stratum == exact0);
                    zero.filter(|r| r.nonfinite > 0)
                        .map(move |r| format!("{} {field}: {}", r.fn_id, r.nonfinite))
                })
                .collect();
            let seen = if hits.is_empty() {
                "nonfinite = 0".to_string()
            } else {
                hits.join(", ")
            };
            Ok((!hits.is_empty(), seen))
        }
        Mechanism::SweepRank => {
            let r = rows
                .ranking
                .as_ref()
                .ok_or("no sweep ranks a subject at f32")?;
            let seen = format!(
                "c {:.3e} u against {:.3e} u chosen (margin {DOMINATED_BY:.0e}), rank {} of {}",
                r.objective, r.chosen, r.rank, r.of
            );
            Ok((dominated(r), seen))
        }
    }
}

/// The largest `max_u` of `fn_id`, with its stratum.
pub(super) fn worst<'a>(rows: &'a [Row], fn_id: &str) -> Option<&'a Row> {
    let mut of_fn = rows
        .iter()
        .filter(|r| r.fn_id == fn_id && !r.max_u.is_nan());
    of_fn
        .next()
        .map(|first| of_fn.fold(first, |w, r| if r.max_u > w.max_u { r } else { w }))
}

pub(super) struct Report {
    pub(super) text: String,
    pub(super) failures: Vec<String>,
}

/// Runs `cases` over the coefficient ids of `dir`, at `precision` (its strata alone), and judges
/// every mechanism meaningful there against `window`.
fn check(
    dir: &Path,
    cases: Vec<Case>,
    precision: Precision,
    window: (f64, f64),
) -> Result<Report, String> {
    let mut entries = corpus::manifest(dir)?;
    let judged = |e: &corpus::Entry| Coeff::of_fn(&e.fn_id).is_some();
    entries.retain(|e| judged(e) && cases.iter().any(|c| c.subject.subject.supports(&e.fn_id)));
    // One sweep of `c` ranks every subject's `c`.
    let ranker = match precision {
        Precision::F64 => Some(Ranker::new(dir, Coeff::C)?),
        Precision::F32 => None,
    };
    let mut ranked = Vec::new();
    let (mut subjects, mut expected) = (Vec::new(), Vec::new());
    for case in cases {
        ranked.push(
            ranker
                .as_ref()
                .map(|r| r.rank(case.c.0, case.c.1))
                .transpose()?,
        );
        subjects.push(case.subject);
        expected.push(case.must_fire);
    }
    let by = |rule| evaluate_by(dir, &entries, &subjects, precision, &|_| Some(rule));
    let (value, d_branch) = (by(&COEFF_VALUE)?, by(&COEFF_D_BRANCH)?);

    let mut text = String::new();
    let mut failures = Vec::new();
    let mut verdicts = String::new();
    let mut curves = String::new();
    let per_subject = subjects.iter().zip(value).zip(d_branch).zip(ranked);
    for ((((s, value), d_branch), ranking), must_fire) in per_subject.zip(expected) {
        let name = s.subject.name();
        let rows = Rows {
            value,
            d_branch,
            ranking,
        };
        let measured = measured_c(&rows);
        if let Some(ranking) = &rows.ranking {
            if measured.map(f64::to_bits) != Some(ranking.objective.to_bits()) {
                failures.push(format!(
                    "{name}: the sweep ranks `c` at {:e} u, but the subject's `c` measures {} u: \
                     its declared candidate is not what it runs",
                    ranking.objective,
                    measured.map_or("no row".to_string(), |m| format!("{m:e}"))
                ));
            }
        }
        if must_fire.is_none() {
            let version = s.version_at(precision);
            let _ = writeln!(text, "{name} {version}: the max over strata, in u");
            for c in Coeff::ALL {
                let fn_id = format!("coeff_{}", c.name());
                let cell = |rows: &[Row]| {
                    worst(rows, &fn_id).map_or("-".to_string(), |r| {
                        format!("{:.3e} ({})", r.max_u, r.stratum)
                    })
                };
                let _ = writeln!(
                    text,
                    "  {fn_id}  value {:<28} d_branch {}",
                    cell(&rows.value),
                    cell(&rows.d_branch)
                );
            }
        }
        let shown: Vec<String> = curve(&rows, precision)?
            .iter()
            .map(|p| format!("{:.2e}", p.1))
            .collect();
        let _ = writeln!(curves, "  {name:<24} {}", shown.join("  "));
        for &m in Mechanism::at(precision) {
            let (fired, seen) = observe(m, &rows, precision, window)?;
            let (expects, verdict, ok) = match must_fire {
                None => ("silent", if fired { "FAIL" } else { "ok" }, !fired),
                Some(f) if f == m => ("detected", if fired { "ok" } else { "FAIL" }, fired),
                Some(_) => ("-", if fired { "fires" } else { "quiet" }, true),
            };
            let _ = writeln!(
                verdicts,
                "  {name:<24} -> {:<30} -> {expects:<8} {verdict:<5} {seen}",
                m.label()
            );
            if !ok {
                failures.push(format!("{name}: `{}` is not {expects}: {seen}", m.label()));
            }
        }
    }
    let (k, end) = (decades(precision), suffix(precision));
    let (first, last) = (k.end(), k.start());
    let _ = writeln!(
        text,
        "\nvalue max_u of coeff_b at theta:1e-{first}{end} .. theta:1e-{last}{end} (the curve, p in \
         [{}, {}]):\n{curves}\nsubject -> mechanism -> verdict (`detected`: must fire; `silent`: \
         must not; `-`: the other defect's mechanism, printed and not gated):\n{verdicts}",
        window.0, window.1
    );
    Ok(Report { text, failures })
}

/// The halves over `dir` (the coefficients, SO(3), SE_N(3), the envelope's and `0056`'s at
/// binary64, then the coefficients at `f32`): their reports joined, and every failure among them.
fn halves(
    dir: &Path,
    coefficients: Vec<Case>,
    so3: Vec<selftest_so3::Case>,
    se3: Vec<selftest_se3::Case>,
    envelope: Seeded,
) -> Result<Report, String> {
    let coefficients = check(dir, coefficients, Precision::F64, WINDOW)?;
    let so3 = selftest_so3::check(dir, so3, selftest_so3::BAR)?;
    let se3 = selftest_se3::check(dir, se3, selftest_se3::BAR)?;
    let envelope = selftest_envelope::check(dir, envelope)?;
    let linalg = selftest_linalg::check(dir, &selftest_linalg::Defect::ALL)?;
    let f32 = check(dir, cases_f32(), Precision::F32, WINDOW)?;
    let f32_title = "f32 (the coefficients' @f32 strata, u = 2^-24):";
    let text = format!(
        "{}\n{}\n{}\n{}\n{}\n{f32_title}\n{}",
        coefficients.text, so3.text, se3.text, envelope.text, linalg.text, f32.text
    );
    let mut failures = coefficients.failures;
    failures.extend(so3.failures);
    failures.extend(se3.failures);
    failures.extend(envelope.failures);
    failures.extend(linalg.failures);
    failures.extend(f32.failures.into_iter().map(|f| format!("f32: {f}")));
    Ok(Report { text, failures })
}

#[allow(clippy::print_stdout)]
pub(crate) fn run(dir: &Path) -> Result<(), String> {
    let cases = (cases(), selftest_so3::cases(), selftest_se3::cases());
    let planted = selftest_envelope::planted();
    let report = halves(dir, cases.0, cases.1, cases.2, planted)?;
    print!("{}", report.text);
    match report.failures.as_slice() {
        [] => Ok(()),
        failed => Err(format!("self-test: {}", failed.join("; "))),
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;
    use crate::seeded::{d12, Candidate, Series, C_PLANTED, D1};

    /// The verdict lines of a report that name `subject`.
    fn failed(report: &Report, subject: &str) -> bool {
        report.failures.iter().any(|f| f.starts_with(subject))
    }

    fn series() -> Result<Series<D1>, String> {
        Series::load(&corpus_dir()?)
    }

    /// The correct kernels at the `D12` prior, all six coefficients: the candidate the sweep
    /// scores as its prior, not the generated one.
    fn prior(series: &Series<D1>, must_fire: Option<Mechanism>) -> Result<Case, String> {
        Ok(Case::new(Seeded::uniform(series, d12(series)?), must_fire))
    }

    fn run_cases(cases: Vec<Case>, window: (f64, f64)) -> Result<Report, String> {
        check(&corpus_dir()?, cases, Precision::F64, window)
    }

    fn run_f32(cases: Vec<Case>, window: (f64, f64)) -> Result<Report, String> {
        check(&corpus_dir()?, cases, Precision::F32, window)
    }

    /// The lines of a report that begin with `prefix`, runs of blanks collapsed to one.
    fn lines(report: &Report, prefix: &str) -> Vec<String> {
        let squeeze = |l: &str| l.split_whitespace().collect::<Vec<_>>().join(" ");
        let all = report.text.lines().map(squeeze);
        all.filter(|l| l.starts_with(prefix)).collect()
    }

    #[test]
    fn the_self_test_runs_every_half_over_the_committed_corpus() -> Result<(), String> {
        run(&corpus_dir()?)
    }

    #[test]
    fn a_failure_in_the_se3_half_fails_the_self_test() -> Result<(), String> {
        let dir = corpus_dir()?;
        let (so3, se3) = (selftest_so3::cases(), selftest_se3::defect_not_planted());
        let planted = selftest_envelope::planted();
        let report = halves(&dir, cases(), so3, se3, planted)?;
        let want = "seeded:correct: `every stratum of sen3_jr, sen3_jl fails` is not detected";
        assert_eq!(report.failures.len(), 1, "{:?}", report.failures);
        assert!(
            report.failures[0].starts_with(want),
            "{:?}",
            report.failures
        );
        Ok(())
    }

    #[test]
    fn a_failure_in_the_envelope_half_fails_the_self_test() -> Result<(), String> {
        let dir = corpus_dir()?;
        let (so3, se3) = (selftest_so3::cases(), selftest_se3::cases());
        // The correct kernel planted as the defect: nothing for the envelope to detect.
        let report = halves(&dir, cases(), so3, se3, Seeded::generated())?;
        let want = "seeded:correct: `envelope fails against seeded:correct` is not detected";
        assert_eq!(report.failures.len(), 1, "{:?}", report.failures);
        assert!(
            report.failures[0].starts_with(want),
            "{:?}",
            report.failures
        );
        Ok(())
    }

    #[test]
    fn a_power_law_fits_its_exponent_and_a_broken_curve_does_not_fit() {
        let law = |p: i32| -> Vec<(f64, f64)> {
            (2..=8)
                .map(|k| (10f64.powi(-k), 3.0 * 10f64.powi(p * k)))
                .collect()
        };
        for p in [-1, 0, 2, 4] {
            let (fit, r2) = power_law(&law(p)).unwrap_or((f64::NAN, f64::NAN));
            assert!((fit - f64::from(p)).abs() < 1e-12, "{p}: {fit}");
            // A constant curve has no variance to explain; its r2 is not meaningful.
            assert!(p == 0 || (r2 - 1.0).abs() < 1e-12, "{p}: {r2}");
        }
        assert_eq!(power_law(&[(1.0, 1.0)]), None);
        assert_eq!(power_law(&[(1.0, 1.0), (0.1, 0.0)]), None);
        assert_eq!(power_law(&[(1.0, 1.0), (0.1, f64::NAN)]), None);
        assert_eq!(power_law(&[(0.0, 1.0), (0.1, 1.0)]), None);
    }

    #[test]
    fn the_committed_corpus_passes_and_the_report_says_how() -> Result<(), String> {
        let report = run_cases(cases(), WINDOW)?;
        assert_eq!(report.failures, Vec::<String>::new());
        let t = &report.text;
        assert!(
            t.contains("seeded:correct generated") && t.contains("coeff_e  value"),
            "{t}"
        );
        assert!(
            t.contains("seeded:b-no-series") && t.contains("-> detected ok    p = 1."),
            "{t}"
        );
        assert!(t.contains("coeff_k d_branch: 1"), "{t}");
        let c = "seeded:c-two-terms-1e-8 -> sweep ranks c dominated -> detected ok c ";
        assert_eq!(lines(&report, c).len(), 1, "{t}");
        Ok(())
    }

    #[test]
    fn the_window_the_decades_and_the_measured_exponent_are_pinned() -> Result<(), String> {
        // §10: `theta:1e-8` to `theta:1e-2`, p in [1.8, 2.2]; seven points, the saturated first
        // one included, `p` = 1.937 (`docs/maths/index.md`, the open item on the `b` window).
        assert_eq!(WINDOW, (1.8, 2.2));
        assert_eq!(DECADES, 2..=8);
        let report = run_cases(cases(), WINDOW)?;
        let curve = lines(&report, "seeded:b-no-series 9.01e15");
        let points = curve.first().map(|l| l.split(' ').count());
        assert_eq!(points, Some(1 + 7), "{curve:?}");
        let fit = "seeded:b-no-series -> b value curve fits theta^-p -> detected ok p = 1.937";
        assert_eq!(lines(&report, fit).len(), 1, "{}", report.text);
        Ok(())
    }

    #[test]
    fn the_correct_kernels_errors_are_the_sweeps_per_field_maxima() -> Result<(), String> {
        // The generated kernel's max over strata, per output field, is the CSV's, digit for digit
        // as printed (`{:.3e}`): the harness and the sweep score the same kernel by two routes, at
        // both precisions (at `f32` the `f32` rows and the `@f32` strata).
        let csv_path = crate::conformance::root()?.join("conformance/sweeps/thresholds-seeded.csv");
        let csv = std::fs::read_to_string(&csv_path).map_err(|e| e.to_string())?;
        let rows: Vec<Vec<&str>> = csv
            .lines()
            .skip(1)
            .map(|r| r.split(',').collect())
            .collect();
        assert_eq!(rows.len(), 16);
        let generated = || vec![Case::new(Seeded::generated(), None)];
        for (precision, report) in [
            ("f64", run_cases(generated(), WINDOW)?),
            ("f32", run_f32(generated(), WINDOW)?),
        ] {
            for c in Coeff::ALL {
                let row = rows.iter().find(|r| r[0] == c.name() && r[1] == precision);
                let row = row.ok_or_else(|| format!("no `{precision}` row for {c:?}"))?;
                let shown = |x: &str| {
                    x.parse::<f64>()
                        .map(|x| format!("{x:.3e}"))
                        .map_err(|e| e.to_string())
                };
                // Columns 7 and 8 of the CSV: `value_max_u`, `deriv_max_u`.
                let want = format!(
                    "coeff_{} value {} d_branch {}",
                    c.name(),
                    shown(row[7])?,
                    shown(row[8])?
                );
                let got = lines(&report, &format!("coeff_{}", c.name()));
                let got: Vec<String> = got
                    .iter()
                    .map(|l| {
                        let words: Vec<&str> = l.split(' ').collect();
                        // `coeff_x value V (stratum) d_branch D (stratum)`
                        format!(
                            "{} {} {} {} {}",
                            words[0], words[1], words[2], words[4], words[5]
                        )
                    })
                    .collect();
                assert_eq!(got, [want], "{c:?} {precision}");
            }
        }
        Ok(())
    }

    #[test]
    fn the_d12_priors_errors_are_pinned_per_coefficient() -> Result<(), String> {
        // Bit-deterministic (D16): the max over strata, in u, of the D12 candidate, per output
        // field. A changed series length, switch, operand or exact-arm form moves a cell.
        let report = run_cases(vec![prior(&series()?, None)?], WINDOW)?;
        let want = [
            "coeff_k value 1.315e0 (theta:dense) d_branch 8.689e3 (theta:dense)",
            "coeff_a value 4.523e1 (theta:dense) d_branch 2.225e5 (theta:dense)",
            "coeff_b value 3.660e2 (theta:1e-1) d_branch 1.159e6 (theta:dense)",
            "coeff_c value 2.144e3 (theta:dense) d_branch 2.750e7 (theta:dense)",
            "coeff_d value 2.826e3 (theta:dense) d_branch 1.166e7 (theta:dense)",
            "coeff_e value 2.131e6 (theta:dense) d_branch 9.221e9 (theta:dense)",
        ];
        assert_eq!(lines(&report, "coeff_"), want);
        assert!(report.text.contains("seeded:correct 4terms-z0.01"));
        Ok(())
    }

    #[test]
    fn a_defect_that_is_not_planted_is_a_failure_and_so_is_a_noisy_correct_kernel(
    ) -> Result<(), String> {
        // Each defect's slot holds the correct kernel: nothing fires, so all three are missed.
        let correct = || Case::new(Seeded::generated(), None);
        let expecting = |m| Case::new(Seeded::generated(), Some(m));
        let missed = run_cases(
            vec![
                correct(),
                expecting(Mechanism::ErrorCurve),
                expecting(Mechanism::NonfiniteAtZero),
                expecting(Mechanism::SweepRank),
            ],
            WINDOW,
        )?;
        assert_eq!(missed.failures.len(), 3, "{:?}", missed.failures);
        assert!(missed
            .failures
            .iter()
            .all(|f| f.contains("is not detected")));
        // A defect where the correct kernel belongs: its own mechanism breaks the silence.
        for d in Defect::COEFFICIENT {
            let loud = run_cases(vec![Case::new(Seeded::planted(d), None)], WINDOW)?;
            assert!(failed(&loud, &format!("seeded:{}", d.name())), "{d:?}");
        }
        Ok(())
    }

    #[test]
    fn the_b_defect_also_fires_the_nonfinite_mechanism_and_the_k_defect_not_the_curve(
    ) -> Result<(), String> {
        // `b` by its definition is 0/0 at z = 0. The other mechanism is printed, not gated.
        let report = run_cases(cases(), WINDOW)?;
        assert_eq!(report.failures, Vec::<String>::new());
        let b = "seeded:b-no-series -> nonfinite in theta:exact0 -> - fires coeff_b value: 1";
        let k = "seeded:k-sqrt-unsafe -> b value curve fits theta^-p -> - quiet";
        for cell in [b, k] {
            assert_eq!(lines(&report, cell).len(), 1, "{cell}\n{}", report.text);
        }
        Ok(())
    }

    #[test]
    fn moving_the_window_breaks_the_self_test_in_both_directions() -> Result<(), String> {
        // Too high: the defect's fitted p is below it, the defect is missed.
        let high = run_cases(cases(), (2.5, 3.0))?;
        assert!(failed(&high, "seeded:b-no-series"), "{:?}", high.failures);
        // Wide open: the correct kernel's curve (p < 0) fits it, the silence is broken.
        let open = run_cases(cases(), (-1.0, 3.0))?;
        assert!(failed(&open, "seeded:correct"), "{:?}", open.failures);
        // The named mechanism is a different one for `k`: the window does not touch it.
        assert!(!failed(&high, "seeded:k-sqrt-unsafe"));
        Ok(())
    }

    #[test]
    fn a_candidate_that_switches_below_the_cancellation_is_the_b_defect_in_disguise(
    ) -> Result<(), String> {
        // The correct kernel with its switch at z = 1e-16 (theta = 1e-8): `b` is the exact arm
        // across the window, and the curve fits theta^-2 like the planted defect.
        let series = series()?;
        let late = Candidate::new(&series, 4, "1e-16")?;
        let case = Case::new(Seeded::uniform(&series, late), Some(Mechanism::ErrorCurve));
        let report = run_cases(vec![case], WINDOW)?;
        assert_eq!(report.failures, Vec::<String>::new());
        Ok(())
    }

    /// The swept candidates: `TERMS` series lengths times the admissible grid, which `0039`'s
    /// domain rule ends at the last point below `π²`. Derived rather than typed, so that lifting
    /// either limit again moves this with it.
    const SWEPT: usize = 17408;

    #[test]
    fn the_planted_c_is_detected_by_the_sweeps_ranking_and_only_by_it() -> Result<(), String> {
        let report = run_cases(cases(), WINDOW)?;
        assert_eq!(report.failures, Vec::<String>::new());
        // The correct kernel's `c` is the chosen candidate: first of the swept candidates,
        // `TERMS` lengths times the admissible grid (`0039`), ratio 1.
        let silent = lines(
            &report,
            "seeded:correct -> sweep ranks c dominated -> silent ok c ",
        );
        assert!(
            silent.len() == 1 && silent[0].ends_with(&format!("rank 1 of {SWEPT}")),
            "{silent:?}"
        );
        // The planted one is far beyond the margin, among the last of them; the other defects
        // keep the chosen `c` and are quiet.
        let planted = "seeded:c-two-terms-1e-8 -> sweep ranks c dominated -> detected ok c ";
        let planted = lines(&report, planted);
        assert_eq!(planted.len(), 1, "{}", report.text);
        for other in ["b-no-series", "k-sqrt-unsafe"] {
            let quiet = format!("seeded:{other} -> sweep ranks c dominated -> - quiet c ");
            assert_eq!(lines(&report, &quiet).len(), 1, "{other}\n{}", report.text);
        }
        // The planted defect fires no other mechanism.
        let (curve, zero) = (
            "b value curve fits theta^-p -> - quiet",
            "nonfinite in theta:exact0 -> - quiet",
        );
        for m in [curve, zero] {
            let line = format!("seeded:c-two-terms-1e-8 -> {m}");
            assert_eq!(lines(&report, &line).len(), 1, "{m}\n{}", report.text);
        }
        Ok(())
    }

    #[test]
    fn a_subject_that_is_not_what_its_declared_candidate_says_fails_the_run() -> Result<(), String>
    {
        let mismatched = |case: Case| -> Result<bool, String> {
            let report = run_cases(vec![case], WINDOW)?;
            let says = "declared candidate is not what it runs";
            Ok(report.failures.iter().any(|f| f.contains(says)))
        };
        let chosen = Case::new(Seeded::generated(), None).c;
        // The correct kernel declaring the planted `c`: the sweep ranks a `c` it does not run.
        let subject = || Seeded::generated().registered();
        assert!(mismatched(Case {
            subject: subject(),
            must_fire: None,
            c: C_PLANTED
        })?);
        // The planted `c` declaring the chosen one: the ranking is silent and would go unseen
        // were the declaration not checked against the output.
        assert!(mismatched(Case {
            subject: Seeded::planted(Defect::CTwoTermsEarly).registered(),
            must_fire: Some(Mechanism::SweepRank),
            c: chosen
        })?);
        // Where the declaration is the subject's, the two agree to the bit.
        assert!(!mismatched(Case::new(Seeded::generated(), None))?);
        let planted = Seeded::planted(Defect::CTwoTermsEarly);
        assert!(!mismatched(Case::new(planted, Some(Mechanism::SweepRank)))?);
        // Also where `c` is not finite: the exact arm alone is `0/0` at `z = 0`, `inf` to the
        // sweep and, its records counted as non-finite, to the harness.
        let series = series()?;
        let exact = Candidate::new(&series, 1, "0")?;
        assert!(!mismatched(Case::new(
            Seeded::uniform(&series, exact),
            None
        ))?);
        Ok(())
    }

    #[test]
    fn a_candidate_is_dominated_only_beyond_the_stated_margin() -> Result<(), String> {
        let at = |objective| Ranking {
            objective,
            chosen: 2.0,
            rank: 1,
            of: SWEPT,
        };
        assert_eq!(DOMINATED_BY, 1e6);
        assert!(!dominated(&at(2.0)) && !dominated(&at(2.0 * DOMINATED_BY)));
        assert!(dominated(&at(2.0 * DOMINATED_BY * 1.000_001)) && dominated(&at(f64::INFINITY)));
        assert!(!dominated(&at(f64::NAN)));
        // The `D12` prior's `c` is worse than the chosen one, by tens and not by 10^6: not a
        // defect, and a self-test expecting the sweep to rank it dominated fails.
        let series = series()?;
        let report = run_cases(vec![prior(&series, Some(Mechanism::SweepRank))?], WINDOW)?;
        assert!(failed(&report, "seeded:correct"), "{:?}", report.failures);
        // The same kernel with `c` switched at z = 1e-16 and cut to two terms is.
        let planted = Candidate::new(&series, 2, "1e-16")?;
        let case = Case::new(
            Seeded::uniform(&series, planted),
            Some(Mechanism::SweepRank),
        );
        assert_eq!(
            run_cases(vec![case], WINDOW)?.failures,
            Vec::<String>::new()
        );
        Ok(())
    }

    #[test]
    fn the_f32_half_detects_its_defects_and_the_correct_kernel_is_silent() -> Result<(), String> {
        // Fitted `p`: binary64 1.937 (seven strata, the plateau one of them), binary32 2.134
        // (three, above the plateau); the correct kernel's curve does not fit at either.
        let report = run_f32(cases_f32(), WINDOW)?;
        assert_eq!(report.failures, Vec::<String>::new());
        let (curve, zero) = ("b value curve fits theta^-p", "nonfinite in theta:exact0");
        let want = [
            format!("seeded:b-no-series -> {curve} -> detected ok p = 2.134, r2 = 1.000"),
            format!("seeded:correct -> {curve} -> silent ok p = -0.021, r2 = 0.920"),
            format!("seeded:correct -> {zero} -> silent ok nonfinite = 0"),
            format!("seeded:k-sqrt-unsafe -> {zero} -> detected ok coeff_k d_branch: 1"),
        ];
        for line in &want {
            assert_eq!(lines(&report, line).len(), 1, "{line}\n{}", report.text);
        }
        // The planted `c` is the sweep's, and the sweep is binary64's: no `f32` subject is ranked.
        assert!(!report.text.contains("sweep ranks"), "{}", report.text);
        Ok(())
    }

    #[test]
    fn a_moved_f32_window_or_a_missed_f32_defect_fails_the_f32_half() -> Result<(), String> {
        let high = run_f32(cases_f32(), (2.5, 3.0))?;
        assert!(failed(&high, "seeded:b-no-series"), "{:?}", high.failures);
        let expecting = Case::new(Seeded::generated(), Some(Mechanism::ErrorCurve));
        let missed = run_f32(vec![expecting], WINDOW)?;
        assert_eq!(missed.failures.len(), 1, "{:?}", missed.failures);
        Ok(())
    }

    #[test]
    fn the_plateau_stratum_would_pull_the_f32_fit_below_the_window() -> Result<(), String> {
        // `theta:1e-4@f32` holds `sqrt(6u) = 6e-4`: `b` is 0 there, the error `1/u` = 2^24, and
        // with it the four strata fit `p = 1.64`; the three above it fit 2.134.
        let dir = corpus_dir()?;
        let b = Seeded::planted(Defect::BNoSeries).registered();
        let entries: Vec<_> = corpus::manifest(&dir)?
            .into_iter()
            .filter(|e| e.fn_id == "coeff_b")
            .collect();
        let rows = evaluate_by(&dir, &entries, &[b], Precision::F32, &|_| {
            Some(&COEFF_VALUE)
        })?;
        let max = |k: i32| {
            let name = format!("theta:1e-{k}@f32");
            rows[0].iter().find(|r| r.stratum == name).map(|r| r.max_u)
        };
        assert_eq!(max(4), Some(16_777_216.0));
        let fit = |ks: &[i32]| {
            let points: Option<Vec<(f64, f64)>> = ks
                .iter()
                .map(|&k| Some((10f64.powi(-k), max(k)?)))
                .collect();
            points
                .and_then(|p| power_law(&p))
                .map(|(p, _)| (p * 1e3).round() / 1e3)
        };
        assert_eq!(fit(&[4, 3, 2, 1]), Some(1.644));
        assert_eq!(fit(&[3, 2, 1]), Some(2.134));
        // `theta:1e-1@f32` is above D12's switch: the two strata below it fit 0.012 from the
        // window's edge (0014 (draft) question 26).
        assert_eq!(fit(&[3, 2]), Some(2.188));
        Ok(())
    }
}
