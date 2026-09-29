//! `cargo xtask conformance --self-test` (`docs/PHASE1.md` §10): the instrument shown to detect what
//! it claims to. It runs the correct seeded kernel, which no mechanism may fire on, and every
//! planted defect (`crate::seeded`), which must fire its named mechanism, over the `coeff_*` corpus
//! ids, and fails when either does not hold. Nothing is written to `conformance/results/`.
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
//! - **The correct kernel** is the `tf_tree` D12 prior for every coefficient, an evaluation
//!   candidate, not an optimum; its errors are printed. `NUMERICS.md` §4 lists the prior for `a`,
//!   `b`, `c` only, so the rows of `k`, `d`, `e` are printed as D12 applied to a coefficient it was
//!   not defined for.
//! - **Only the named mechanism is gated** for a defect. The other mechanism's reading is printed
//!   and not gated: `b` by its definition is `0/0` at `z = 0`, so it fires `nonfinite` too
//!   (`the_b_defect_also_fires_the_nonfinite_mechanism_and_the_k_defect_not_the_curve`).

use std::fmt::Write;
use std::path::Path;

use helicoid_linalg::Precision;

use super::metric::{COEFF_D_BRANCH, COEFF_VALUE};
use super::report::Row;
use super::subject::Registered;
use super::{corpus, evaluate_by};
use crate::seeded::{d12, Defect, Seeded, Series, D1};

/// The exponent `p` of `θ^-p` that the `b` defect's curve must fit, inclusive (§10).
const WINDOW: (f64, f64) = (1.8, 2.2);
/// The strata `theta:1e-8` to `theta:1e-2` (§10).
const DECADES: std::ops::RangeInclusive<i32> = 2..=8;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mechanism {
    ErrorCurve,
    NonfiniteAtZero,
}

impl Mechanism {
    const ALL: [Mechanism; 2] = [Mechanism::ErrorCurve, Mechanism::NonfiniteAtZero];

    fn label(self) -> &'static str {
        match self {
            Mechanism::ErrorCurve => "b value curve fits theta^-p",
            Mechanism::NonfiniteAtZero => "nonfinite in theta:exact0",
        }
    }

    /// The mechanism that must detect `defect` (§10).
    fn of(defect: Defect) -> Self {
        match defect {
            Defect::BNoSeries => Mechanism::ErrorCurve,
            Defect::KSqrtUnsafe => Mechanism::NonfiniteAtZero,
        }
    }
}

/// A subject and the mechanism it must fire; `None` for the correct kernel, on which none may.
struct Case {
    subject: Registered,
    must_fire: Option<Mechanism>,
}

fn cases(series: &Series<D1>) -> Result<Vec<Case>, String> {
    let correct = Seeded::correct(series.clone(), d12(series)?);
    let mut all = vec![Case {
        subject: correct.registered(),
        must_fire: None,
    }];
    for d in Defect::ALL {
        all.push(Case {
            subject: Seeded::planted(series.clone(), d)?.registered(),
            must_fire: Some(Mechanism::of(d)),
        });
    }
    Ok(all)
}

/// The coefficient ids in print order, and whether `NUMERICS.md` §4 lists the D12 prior for them.
const COEFFS: [(&str, bool); 6] = [
    ("coeff_a", true),
    ("coeff_b", true),
    ("coeff_c", true),
    ("coeff_k", false),
    ("coeff_d", false),
    ("coeff_e", false),
];

/// One subject's rows over the coefficient ids, one output field at a time.
struct Rows {
    value: Vec<Row>,
    d_branch: Vec<Row>,
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
    }
}

/// The largest `max_u` of `fn_id`, with its stratum.
fn worst<'a>(rows: &'a [Row], fn_id: &str) -> Option<&'a Row> {
    let mut of_fn = rows
        .iter()
        .filter(|r| r.fn_id == fn_id && !r.max_u.is_nan());
    of_fn
        .next()
        .map(|first| of_fn.fold(first, |w, r| if r.max_u > w.max_u { r } else { w }))
}

struct Report {
    text: String,
    failures: Vec<String>,
}

/// Runs `cases` over the coefficient ids of `dir` and judges every mechanism against `window`.
fn check(dir: &Path, cases: Vec<Case>, window: (f64, f64)) -> Result<Report, String> {
    let mut entries = corpus::manifest(dir)?;
    entries.retain(|e| cases.iter().any(|c| c.subject.subject.supports(&e.fn_id)));
    let (subjects, expected): (Vec<Registered>, Vec<Option<Mechanism>>) =
        cases.into_iter().map(|c| (c.subject, c.must_fire)).unzip();
    let by = |rule| evaluate_by(dir, &entries, &subjects, Precision::F64, &|_| Some(rule));
    let (value, d_branch) = (by(&COEFF_VALUE)?, by(&COEFF_D_BRANCH)?);

    let mut text = String::new();
    let mut failures = Vec::new();
    let mut verdicts = String::new();
    let mut curves = String::new();
    let per_subject = subjects.iter().zip(value).zip(d_branch).zip(expected);
    for (((s, value), d_branch), must_fire) in per_subject {
        let name = s.subject.name();
        let rows = Rows { value, d_branch };
        if must_fire.is_none() {
            let _ = writeln!(text, "{name} {}: the max over strata, in u", s.version);
            let mut group = None;
            for (fn_id, prior) in COEFFS {
                if group != Some(prior) {
                    group = Some(prior);
                    let _ = writeln!(
                        text,
                        "  {}",
                        if prior {
                            "the D12 prior (NUMERICS.md §4):"
                        } else {
                            "D12 applied to a coefficient it was not defined for:"
                        }
                    );
                }
                let cell = |rows: &[Row]| {
                    worst(rows, fn_id).map_or("-".to_string(), |r| {
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
        let _ = writeln!(curves, "  {name:<22} {}", shown.join("  "));
        for m in Mechanism::ALL {
            let (fired, seen) = observe(m, &rows, window)?;
            let (expects, verdict, ok) = match must_fire {
                None => ("silent", if fired { "FAIL" } else { "ok" }, !fired),
                Some(f) if f == m => ("detected", if fired { "ok" } else { "FAIL" }, fired),
                Some(_) => ("-", if fired { "fires" } else { "quiet" }, true),
            };
            let _ = writeln!(
                verdicts,
                "  {name:<22} -> {:<30} -> {expects:<8} {verdict:<5} {seen}",
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
    let report = check(dir, cases(&Series::load(dir)?)?, WINDOW)?;
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
    use crate::seeded::Candidate;

    /// The verdict lines of a report that name `subject`.
    fn failed(report: &Report, subject: &str) -> bool {
        report.failures.iter().any(|f| f.starts_with(subject))
    }

    fn series() -> Result<Series<D1>, String> {
        Series::load(&corpus_dir()?)
    }

    fn correct(series: &Series<D1>, must_fire: Option<Mechanism>) -> Result<Case, String> {
        Ok(Case {
            subject: Seeded::correct(series.clone(), d12(series)?).registered(),
            must_fire,
        })
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
        let report = run_cases(cases(&series()?)?, WINDOW)?;
        assert_eq!(report.failures, Vec::<String>::new());
        let t = &report.text;
        assert!(
            t.contains("seeded:correct 4terms-z0.01") && t.contains("coeff_e  value"),
            "{t}"
        );
        assert!(
            t.contains("seeded:b-no-series") && t.contains("-> detected ok    p = 1."),
            "{t}"
        );
        assert!(t.contains("coeff_k d_branch: 1"), "{t}");
        Ok(())
    }

    #[test]
    fn the_window_the_decades_and_the_measured_exponent_are_pinned() -> Result<(), String> {
        // §10: `theta:1e-8` to `theta:1e-2`, p in [1.8, 2.2]; seven points, the saturated first
        // one included, `p` = 1.937 (`docs/maths/index.md`, the open item on the `b` window).
        assert_eq!(WINDOW, (1.8, 2.2));
        assert_eq!(DECADES, 2..=8);
        let report = run_cases(cases(&series()?)?, WINDOW)?;
        let curve = lines(&report, "seeded:b-no-series 9.01e15");
        let points = curve.first().map(|l| l.split(' ').count());
        assert_eq!(points, Some(1 + 7), "{curve:?}");
        let fit = "seeded:b-no-series -> b value curve fits theta^-p -> detected ok p = 1.937";
        assert_eq!(lines(&report, fit).len(), 1, "{}", report.text);
        Ok(())
    }

    #[test]
    fn the_correct_kernels_errors_are_pinned_per_coefficient() -> Result<(), String> {
        // Bit-deterministic (D16): the max over strata, in u, of the D12 candidate, per output
        // field. A changed series length, switch, operand or exact-arm form moves a cell.
        let report = run_cases(vec![correct(&series()?, None)?], WINDOW)?;
        let want = [
            "coeff_a value 4.523e1 (theta:dense) d_branch 2.225e5 (theta:dense)",
            "coeff_b value 3.660e2 (theta:1e-1) d_branch 1.159e6 (theta:dense)",
            "coeff_c value 2.144e3 (theta:dense) d_branch 2.750e7 (theta:dense)",
            "coeff_k value 1.315e0 (theta:dense) d_branch 8.689e3 (theta:dense)",
            "coeff_d value 2.826e3 (theta:dense) d_branch 1.166e7 (theta:dense)",
            "coeff_e value 2.131e6 (theta:dense) d_branch 9.221e9 (theta:dense)",
        ];
        assert_eq!(lines(&report, "coeff_"), want);
        // The rows of `k`, `d`, `e` say the prior was not defined for them.
        let t = &report.text;
        let (prior, other) = (
            t.find("the D12 prior"),
            t.find("D12 applied to a coefficient"),
        );
        let (c, k) = (t.find("coeff_c  value"), t.find("coeff_k  value"));
        assert!(prior < c && c < other && other < k, "{t}");
        Ok(())
    }

    #[test]
    fn a_defect_that_is_not_planted_is_a_failure_and_so_is_a_noisy_correct_kernel(
    ) -> Result<(), String> {
        let series = series()?;
        // Each defect's slot holds the correct kernel: nothing fires, so both are missed.
        let missed = run_cases(
            vec![
                correct(&series, None)?,
                correct(&series, Some(Mechanism::ErrorCurve))?,
                correct(&series, Some(Mechanism::NonfiniteAtZero))?,
            ],
            WINDOW,
        )?;
        assert_eq!(missed.failures.len(), 2, "{:?}", missed.failures);
        assert!(missed
            .failures
            .iter()
            .all(|f| f.contains("is not detected")));
        // A defect where the correct kernel belongs: both mechanisms' silence is broken.
        for d in Defect::ALL {
            let case = Case {
                subject: Seeded::planted(series.clone(), d)?.registered(),
                must_fire: None,
            };
            let loud = run_cases(vec![case], WINDOW)?;
            assert!(failed(&loud, &format!("seeded:{}", d.name())), "{d:?}");
        }
        Ok(())
    }

    #[test]
    fn the_b_defect_also_fires_the_nonfinite_mechanism_and_the_k_defect_not_the_curve(
    ) -> Result<(), String> {
        // `b` by its definition is 0/0 at z = 0. The other mechanism is printed, not gated.
        let report = run_cases(cases(&series()?)?, WINDOW)?;
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
        let series = series()?;
        // Too high: the defect's fitted p is below it, the defect is missed.
        let high = run_cases(cases(&series)?, (2.5, 3.0))?;
        assert!(failed(&high, "seeded:b-no-series"), "{:?}", high.failures);
        // Wide open: the correct kernel's curve (p < 0) fits it, the silence is broken.
        let open = run_cases(cases(&series)?, (-1.0, 3.0))?;
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
        let case = Case {
            subject: Seeded::correct(series, late).registered(),
            must_fire: Some(Mechanism::ErrorCurve),
        };
        let report = run_cases(vec![case], WINDOW)?;
        assert_eq!(report.failures, Vec::<String>::new());
        Ok(())
    }
}
