//! `cargo xtask conformance --self-test` (`docs/PHASE1.md` §10): the instrument shown to detect what
//! it claims to. It runs the correct seeded kernel, which no mechanism may fire on, and every
//! planted defect (`crate::seeded`), which must fire its named mechanism, over the `coeff_*` corpus
//! ids (the `Log` defects over `so3_log`: `selftest_so3`), and fails when either does not hold.
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
//! - **Non-finite** counts the (record, field) pairs of `theta:exact0` in every coefficient id, so a
//!   NaN in `value` or in `d_branch` both count.
//! - **The sweep's ranking** is for `c` with switch `1e-8` (as `θ`: 0014 (draft) question 14) and
//!   two terms (`seeded::C_PLANTED`): the sweep of `coeff_c` (`crate::thresholds`) scores that
//!   candidate like any it sweeps, and it is detected when its objective exceeds the chosen
//!   candidate's by more than [`DOMINATED_BY`], a factor of `10^6` (question 15), "strictly worse
//!   by a stated margin"; §10's other half, the envelope, is not started. The reading printed is the objective, the chosen one and the candidate's rank among
//!   the swept ones. Every subject is measured with its own `c`, so the correct kernel, whose `c` is
//!   the chosen candidate, must be silent. The rank is of the candidate a subject declares, so the
//!   run also fails when the objective the harness measures for the subject's `c` is not the ranked
//!   one, to the bit: a subject cannot change what it computes and keep its declaration.
//! - **The correct kernel** runs the generated switches (`seeded::generated`, from
//!   `conformance/sweeps/thresholds.csv`), not an optimum in general: the sweep's choice over this
//!   corpus and grid. Its errors are printed, and equal the CSV's per-field maxima.
//! - **Only the named mechanism is gated** for a defect. The other mechanism's reading is printed
//!   and not gated: `b` by its definition is `0/0` at `z = 0`, so it fires `nonfinite` too
//!   (`the_b_defect_also_fires_the_nonfinite_mechanism_and_the_k_defect_not_the_curve`).

use std::fmt::Write;
use std::path::Path;

use helicoid_linalg::Precision;

use super::metric::{COEFF_D_BRANCH, COEFF_VALUE};
use super::report::Row;
use super::subject::Registered;
use super::{corpus, evaluate_by, selftest_so3};
use crate::seeded::{Coeff, Defect, Seeded};
use crate::thresholds::{Ranker, Ranking};

/// The exponent `p` of `θ^-p` that the `b` defect's curve must fit, inclusive (§10).
const WINDOW: (f64, f64) = (1.8, 2.2);
/// The strata `theta:1e-8` to `theta:1e-2` (§10).
const DECADES: std::ops::RangeInclusive<i32> = 2..=8;
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

    fn label(self) -> &'static str {
        match self {
            Mechanism::ErrorCurve => "b value curve fits theta^-p",
            Mechanism::NonfiniteAtZero => "nonfinite in theta:exact0",
            Mechanism::SweepRank => "sweep ranks c dominated",
        }
    }

    /// The mechanism that must detect the coefficient `defect` (§10); the `Log` defects are
    /// judged by `so3` and have none here.
    fn of(defect: Defect) -> Option<Self> {
        match defect {
            Defect::BNoSeries => Some(Mechanism::ErrorCurve),
            Defect::KSqrtUnsafe => Some(Mechanism::NonfiniteAtZero),
            Defect::CTwoTermsEarly => Some(Mechanism::SweepRank),
            Defect::LogAcos | Defect::LogNoFlip => None,
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

/// One subject's rows over the coefficient ids, one output field at a time, and where the sweep
/// ranks its `c`.
struct Rows {
    value: Vec<Row>,
    d_branch: Vec<Row>,
    ranking: Ranking,
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

/// The `max_u` of `coeff_b`'s value at `theta:1e-8` to `theta:1e-2`, as `(θ, max_u)`.
fn curve(rows: &Rows) -> Result<Vec<(f64, f64)>, String> {
    let mut points = Vec::new();
    for k in DECADES.rev() {
        let stratum = format!("theta:1e-{k}");
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
    window: (f64, f64),
) -> Result<(bool, String), String> {
    match mechanism {
        Mechanism::ErrorCurve => Ok(match power_law(&curve(rows)?) {
            Some((p, r2)) => (
                (window.0..=window.1).contains(&p),
                format!("p = {p:.3}, r2 = {r2:.3}"),
            ),
            None => (false, "no curve: a zero or non-finite max".to_string()),
        }),
        Mechanism::NonfiniteAtZero => {
            let fields = [("value", &rows.value), ("d_branch", &rows.d_branch)];
            let hits: Vec<String> = fields
                .iter()
                .flat_map(|&(field, rows)| {
                    let zero = rows.iter().filter(|r| r.stratum == "theta:exact0");
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
            let r = &rows.ranking;
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

/// Runs `cases` over the coefficient ids of `dir` and judges every mechanism against `window`.
fn check(dir: &Path, cases: Vec<Case>, window: (f64, f64)) -> Result<Report, String> {
    let mut entries = corpus::manifest(dir)?;
    let judged = |e: &corpus::Entry| Coeff::of_fn(&e.fn_id).is_some();
    entries.retain(|e| judged(e) && cases.iter().any(|c| c.subject.subject.supports(&e.fn_id)));
    // One sweep of `c` ranks every subject's `c`.
    let ranker = Ranker::new(dir, Coeff::C)?;
    let mut ranked = Vec::new();
    let (mut subjects, mut expected) = (Vec::new(), Vec::new());
    for case in cases {
        ranked.push(ranker.rank(case.c.0, case.c.1)?);
        subjects.push(case.subject);
        expected.push(case.must_fire);
    }
    let by = |rule| evaluate_by(dir, &entries, &subjects, Precision::F64, &|_| Some(rule));
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
        if measured.map(f64::to_bits) != Some(rows.ranking.objective.to_bits()) {
            failures.push(format!(
                "{name}: the sweep ranks `c` at {:e} u, but the subject's `c` measures {} u: its \
                 declared candidate is not what it runs",
                rows.ranking.objective,
                measured.map_or("no row".to_string(), |m| format!("{m:e}"))
            ));
        }
        if must_fire.is_none() {
            let _ = writeln!(text, "{name} {}: the max over strata, in u", s.version);
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
        let shown: Vec<String> = curve(&rows)?
            .iter()
            .map(|p| format!("{:.2e}", p.1))
            .collect();
        let _ = writeln!(curves, "  {name:<24} {}", shown.join("  "));
        for m in Mechanism::ALL {
            let (fired, seen) = observe(m, &rows, window)?;
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
    let _ = writeln!(
        text,
        "\nvalue max_u of coeff_b at theta:1e-8 .. theta:1e-2 (the curve, p in [{}, {}]):\n{curves}\
         \nsubject -> mechanism -> verdict (`detected`: must fire; `silent`: must not; `-`: the other \
         defect's mechanism, printed and not gated):\n{verdicts}",
        window.0, window.1
    );
    Ok(Report { text, failures })
}

#[allow(clippy::print_stdout)]
pub(crate) fn run(dir: &Path) -> Result<(), String> {
    let coefficients = check(dir, cases(), WINDOW)?;
    let so3 = selftest_so3::check(dir, selftest_so3::cases(), selftest_so3::BAR)?;
    print!("{}\n{}", coefficients.text, so3.text);
    let mut failures = coefficients.failures;
    failures.extend(so3.failures);
    match failures.as_slice() {
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
        check(&corpus_dir()?, cases, window)
    }

    /// The lines of a report that begin with `prefix`, runs of blanks collapsed to one.
    fn lines(report: &Report, prefix: &str) -> Vec<String> {
        let squeeze = |l: &str| l.split_whitespace().collect::<Vec<_>>().join(" ");
        let all = report.text.lines().map(squeeze);
        all.filter(|l| l.starts_with(prefix)).collect()
    }

    #[test]
    fn the_self_test_runs_both_halves_over_the_committed_corpus() -> Result<(), String> {
        run(&corpus_dir()?)
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
        // as printed (`{:.3e}`): the harness and the sweep score the same kernel by two routes.
        let report = run_cases(vec![Case::new(Seeded::generated(), None)], WINDOW)?;
        let csv_path = crate::conformance::root()?.join("conformance/sweeps/thresholds.csv");
        let csv = std::fs::read_to_string(&csv_path).map_err(|e| e.to_string())?;
        let rows: Vec<Vec<&str>> = csv
            .lines()
            .skip(1)
            .map(|r| r.split(',').collect())
            .collect();
        assert_eq!(rows.len(), 6);
        for (c, row) in Coeff::ALL.into_iter().zip(rows) {
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
            assert_eq!(got, [want], "{c:?}");
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

    #[test]
    fn the_planted_c_is_detected_by_the_sweeps_ranking_and_only_by_it() -> Result<(), String> {
        let report = run_cases(cases(), WINDOW)?;
        assert_eq!(report.failures, Vec::<String>::new());
        // The correct kernel's `c` is the chosen candidate: first of the 8200 swept, ratio 1.
        let silent = lines(
            &report,
            "seeded:correct -> sweep ranks c dominated -> silent ok c ",
        );
        assert!(
            silent.len() == 1 && silent[0].ends_with("rank 1 of 8200"),
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
            of: 8200,
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
}
