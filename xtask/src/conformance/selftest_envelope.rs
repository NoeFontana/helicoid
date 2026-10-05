//! The envelope half of `--self-test` (`docs/PHASE1.md` §10): the planted `c` "fails against the
//! correct seeded kernel". Both run over `coeff_c` with the harness's rule, a record scored by the
//! larger of its fields, as `just conformance` writes them, and the domination bar of
//! `crate::envelope` reads the correct kernel as the oracle: the planted subject must lose on some
//! stratum, and the correct one on none against itself. The other half of the defect, the sweep's
//! ranking, is `selftest`'s.

use std::fmt::Write;
use std::path::Path;

use helicoid_linalg::Precision;

use super::report::Row;
use super::selftest::Report;
use super::{corpus, evaluate};
use crate::envelope::beaten_strata;
use crate::seeded::{Defect, Seeded};

/// The subject this half plants.
pub(super) fn planted() -> Seeded {
    Seeded::planted(Defect::CTwoTermsEarly)
}

pub(super) fn check(dir: &Path, planted: Seeded) -> Result<Report, String> {
    let entries: Vec<_> = corpus::manifest(dir)?
        .into_iter()
        .filter(|e| e.fn_id == "coeff_c")
        .collect();
    let subjects = vec![Seeded::generated().registered(), planted.registered()];
    let names: Vec<String> = subjects
        .iter()
        .map(|s| s.subject.name().to_string())
        .collect();
    let per_subject = evaluate(dir, &entries, &subjects, Precision::F64)?;
    let [correct, defect] = <[Vec<Row>; 2]>::try_from(per_subject)
        .map_err(|_| "the envelope half runs two subjects".to_string())?;
    let (own, beaten) = (
        beaten_strata(&correct, &correct),
        beaten_strata(&defect, &correct),
    );
    let mut text = format!(
        "envelope: domination over {} on coeff_c, {} strata\n",
        names[0],
        correct.len()
    );
    let mut failures = Vec::new();
    let _ = writeln!(
        text,
        "  {} loses on {} strata, {} on {}",
        names[1],
        beaten.len(),
        names[0],
        own.len()
    );
    if let Some(first) = beaten.first() {
        let _ = writeln!(text, "  first: {first}");
    }
    if beaten.is_empty() {
        failures.push(format!(
            "{}: `envelope fails against {}` is not detected: no stratum of coeff_c is lost",
            names[1], names[0]
        ));
    }
    if !own.is_empty() {
        failures.push(format!(
            "{}: the envelope fails the correct kernel against itself on {} strata",
            names[0],
            own.len()
        ));
    }
    Ok(Report { text, failures })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;

    #[test]
    fn the_planted_c_loses_on_strata_and_the_correct_kernel_on_none() -> Result<(), String> {
        let report = check(&corpus_dir()?, planted())?;
        assert!(report.failures.is_empty(), "{:?}", report.failures);
        assert!(report
            .text
            .contains("seeded:c-two-terms-1e-8 loses on 10 strata"));
        assert!(
            report.text.contains("seeded:correct on 0"),
            "{}",
            report.text
        );
        Ok(())
    }

    #[test]
    fn a_defect_that_changes_nothing_is_not_detected() -> Result<(), String> {
        let report = check(&corpus_dir()?, Seeded::generated())?;
        let want = "seeded:correct: `envelope fails against seeded:correct` is not detected";
        assert_eq!(report.failures.len(), 1, "{:?}", report.failures);
        assert!(
            report.failures[0].starts_with(want),
            "{:?}",
            report.failures
        );
        Ok(())
    }
}
