//! The `0116` half of `--self-test` (`docs/PHASE1.md` §10): decision 5's six spline defects,
//! planted in the seeded evaluator (`seeded::spline`), against that evaluator clean.
//!
//! - Four defects must read more than [`FACTOR`] times the clean maximum (at least 1, as
//!   `selftest_gamma`'s) on every stratum of their ids but decision 5's exceptions, where the
//!   defect's size, about `‖d‖` or `‖d‖²`, falls under `u`: those are printed and not judged.
//! - `λ̈` without its `h⁻²` must read that far above on the `h:*` strata, and every other stratum,
//!   all at `h = 1`, must read the clean maximum to the bit.
//! - The Gram matrix without its `h^{1−2r}` must miss its 0-ulp bar at `h ≠ 1` and meet it at
//!   `h = 1`, where the clean evaluator meets it everywhere.
//!
//! The clean evaluator's own maxima, over every spline id at both precisions, are printed first:
//! they must be finite, and the Gram's 0.

use std::fmt::Write;
use std::path::Path;

use helicoid_linalg::Precision;

use super::corpus::{self, F32_SUFFIX};
use super::evaluate;
use super::report::Row;
use super::selftest::Report;
use super::subject::{Registered, Subject};
use crate::seeded::spline::{Planted, Spline, IDS};

/// How far above the clean maximum a planted one must read, as `selftest_gamma`'s.
const FACTOR: f64 = 1024.0;

/// What a defect must do, decision 5's table.
enum Criterion {
    /// Above `FACTOR` times the clean max on every stratum but these (their names without `@f32`
    /// apply at both precisions, with it at `f32` alone).
    Exceeds(&'static [&'static str]),
    /// Above on the `h:*` strata, the clean max to the bit on every other.
    OnlyAtH,
    /// Not 0 at `h ≠ 1`, 0 at `h = 1`.
    GramOffOne,
}

fn criterion(p: Planted) -> Criterion {
    match p {
        Planted::NonCumulative | Planted::JrInvForJlInv => Criterion::Exceeds(&["inc:tiny@f32"]),
        Planted::AdSignFlipped | Planted::AdNotInverted => {
            Criterion::Exceeds(&["inc:tiny", "inc:mixed"])
        }
        Planted::NoInvHSquared => Criterion::OnlyAtH,
        Planted::GramUnscaled => Criterion::GramOffOne,
    }
}

fn excepted(except: &[&str], stratum: &str) -> bool {
    except
        .iter()
        .any(|e| *e == stratum || *e == stratum.trim_end_matches(F32_SUFFIX))
}

/// The rows of `subject` over `fn_id` at `precision`.
fn rows(
    dir: &Path,
    fn_id: &str,
    subject: Box<dyn Subject>,
    precision: Precision,
) -> Result<Vec<Row>, String> {
    let entries: Vec<_> = corpus::manifest(dir)?
        .into_iter()
        .filter(|e| e.fn_id == fn_id)
        .collect();
    if entries.is_empty() {
        return Err(format!("no corpus file for `{fn_id}`"));
    }
    let registered = [Registered {
        version: "self-test".into(),
        version_f32: "self-test".into(),
        no_f32: None,
        subject,
        planted: true,
    }];
    Ok(evaluate(dir, &entries, &registered, precision)?.remove(0))
}

/// The record's `inv_h` is the stratum's: `h:1` is the only Gram stratum at `h = 1`.
fn at_unit_h(stratum: &str) -> bool {
    stratum.trim_end_matches(F32_SUFFIX) == "h:1"
}

pub(super) fn check(dir: &Path) -> Result<Report, String> {
    let precisions = [Precision::F64, Precision::F32];
    let mut text = String::from("0116: the clean seeded spline, max u per stratum\n");
    let mut failures = Vec::new();
    for fn_id in IDS {
        for precision in precisions {
            let clean = rows(dir, fn_id, Box::new(Spline(None)), precision)?;
            let cells: Vec<String> = clean
                .iter()
                .map(|r| format!("{} {:.3}", r.stratum, r.max_u))
                .collect();
            let _ = writeln!(text, "  {fn_id} ({precision:?}): {}", cells.join(", "));
            for r in &clean {
                let gram_fails = fn_id == "rn_spline_gram" && r.max_u != 0.0;
                if !r.max_u.is_finite() || r.nonfinite > 0 || gram_fails {
                    failures.push(format!(
                        "seeded:spline: {fn_id} {} reads {} ({} non-finite)",
                        r.stratum, r.max_u, r.nonfinite
                    ));
                }
            }
        }
    }
    let _ = writeln!(
        text,
        "0116: the planted spline defects, > {FACTOR} x the clean max but decision 5's exceptions"
    );
    for planted in Planted::ALL {
        for &fn_id in planted.ids() {
            for precision in precisions {
                let clean = rows(dir, fn_id, Box::new(Spline(None)), precision)?;
                let got = rows(dir, fn_id, Box::new(Spline(Some(planted))), precision)?;
                for (c, p) in clean.iter().zip(&got) {
                    let (stratum, cm, pm) = (&c.stratum, c.max_u, p.max_u);
                    let above = !pm.is_finite() || p.nonfinite > 0 || pm > FACTOR * cm.max(1.0);
                    let (judged, ok) = match criterion(planted) {
                        Criterion::Exceeds(except) if excepted(except, stratum) => (false, true),
                        Criterion::Exceeds(_) => (true, above),
                        Criterion::OnlyAtH if stratum.starts_with("h:") => (true, above),
                        Criterion::OnlyAtH => (true, pm.to_bits() == cm.to_bits()),
                        Criterion::GramOffOne if at_unit_h(stratum) => (true, pm == 0.0),
                        Criterion::GramOffOne => (true, pm > 0.0),
                    };
                    let _ =
                        writeln!(
                        text,
                        "  {} {fn_id} {stratum} ({precision:?}): clean {cm:.3}, planted {pm:.3e}{}",
                        planted.name(),
                        if judged { "" } else { " (an exception, not judged)" }
                    );
                    if !ok {
                        failures.push(format!(
                            "{}: not as decision 5 states on {fn_id} {stratum} ({pm} against {cm})",
                            planted.name()
                        ));
                    }
                }
            }
        }
    }
    Ok(Report { text, failures })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;

    #[test]
    fn every_spline_defect_is_detected_as_decision_5_states() -> Result<(), String> {
        let report = check(&corpus_dir()?)?;
        assert!(
            report.failures.is_empty(),
            "{:?}\n{}",
            report.failures,
            report.text
        );
        Ok(())
    }

    #[test]
    fn an_exception_names_a_stratum_at_both_precisions_or_its_twin_alone() {
        assert!(excepted(&["inc:tiny"], "inc:tiny") && excepted(&["inc:tiny"], "inc:tiny@f32"));
        assert!(excepted(&["inc:tiny@f32"], "inc:tiny@f32"));
        assert!(!excepted(&["inc:tiny@f32"], "inc:tiny") && !excepted(&["inc:tiny"], "inc:large"));
    }
}
