//! The SE_N(3) half of `--self-test` (`docs/PHASE1.md` §10): the correct seeded `Exp`, `J_r` and
//! `J_l` must stay under [`BAR`] on every stratum of `sen3_{exp,jr,jl}_n{1,2,3}` and fire no
//! mechanism, and each planted defect must fire its own, over the committed corpus.
//!
//! Readings where §10 is silent, each the smallest (0014 (draft) questions 21 to 23):
//!
//! - **The bar** is [`BAR`] `u` per stratum on the max, one for all nine ids as for SO(3): 1.35
//!   times the measured 8.90 (`J_l`, `rho:1e4/theta=pi-1e-6`; `Exp` at most 3.47), of the
//!   per-coefficient kernel (question 24). §10 states none (question 21).
//! - **"Fails"** is [`FAILS_AT`] `u`, the line SO(3)'s defects use, a million times the bar. A
//!   non-finite output fails; a `NaN` or unscored stratum does not, so it cannot fire a mechanism
//!   (question 22).
//! - **"Every `rho:*` stratum"** is the 25 cells `rho:<scale>/theta=<θ>` (§4.4) of each
//!   `sen3_exp_n<N>` the subject answers: 75 for the translation-first defect, which is planted at
//!   every `N`. **"`sen3_jr*` fails"** is every stratum (52) of every `sen3_jr_n<N>` and
//!   `sen3_jl_n<N>` the subject answers: the defect is in `Q`, which both use. A file with fewer
//!   strata than that fires nothing. A mechanism reads nothing the subject does not answer
//!   (question 23).
//! - **`Dual` comparison** (§10, `Q`'s second mechanism) needs forward-mode differentiation of
//!   the shipped `Exp`, which is the `helicoid` subject (Phase 3); it is not run here.

use std::collections::BTreeMap;
use std::fmt::Write;
use std::path::Path;

use helicoid_linalg::Precision;

use super::report::Row;
use super::selftest::{worst, Report};
use super::selftest_so3::{over_the_bar, FAILS_AT};
use super::subject::Registered;
use super::{corpus, evaluate};
use crate::seeded::se3::{self, Op};
use crate::seeded::{Defect, Seeded};

/// The largest `max_u` the correct kernel may show on a stratum.
pub(super) const BAR: f64 = 12.0;
/// The `rho:*` cells of a `sen3_exp` file (§4.4: 5 scales times 5 angles).
const RHO_CELLS: usize = 25;
/// The strata of a `sen3_jr` or `sen3_jl` file: 27 `theta:*` (`theta:dense` excluded) and the cells.
const STRATA: usize = 52;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mechanism {
    RhoStrata,
    Jacobians,
}

impl Mechanism {
    const ALL: [Mechanism; 2] = [Mechanism::RhoStrata, Mechanism::Jacobians];

    fn label(self) -> &'static str {
        match self {
            Mechanism::RhoStrata => "every rho:* stratum of sen3_exp fails",
            Mechanism::Jacobians => "every stratum of sen3_jr, sen3_jl fails",
        }
    }

    fn of(defect: Defect) -> Option<Self> {
        match defect {
            Defect::Se3ExpTranslationFirst => Some(Mechanism::RhoStrata),
            Defect::QMinusHalf => Some(Mechanism::Jacobians),
            _ => None,
        }
    }

    /// Whether the mechanism reads `row`.
    fn reads(self, row: &Row) -> bool {
        match (self, se3::parse(&row.fn_id)) {
            (Mechanism::RhoStrata, Some((Op::Exp, _))) => row.stratum.starts_with("rho:"),
            (Mechanism::Jacobians, Some((Op::Jr | Op::Jl, _))) => true,
            _ => false,
        }
    }

    /// The strata a file must have for the mechanism to fire on it.
    fn cells(self) -> usize {
        match self {
            Mechanism::RhoStrata => RHO_CELLS,
            Mechanism::Jacobians => STRATA,
        }
    }
}

/// A subject and the mechanism it must fire (`None` for the correct kernel, on which none may).
pub(super) struct Case {
    subject: Registered,
    must_fire: Option<Mechanism>,
}

impl Case {
    fn new(seeded: Seeded, must_fire: Option<Mechanism>) -> Self {
        Self {
            subject: seeded.registered(),
            must_fire,
        }
    }
}

pub(super) fn cases() -> Vec<Case> {
    let mut all = vec![Case::new(Seeded::generated(), None)];
    for d in Defect::SE3 {
        all.push(Case::new(Seeded::planted(d), Mechanism::of(d)));
    }
    all
}

/// The correct kernel and a second one expected to fire the `sen3_jr` mechanism, which it cannot:
/// the cases of a run whose defect is not planted.
#[cfg(test)]
pub(super) fn defect_not_planted() -> Vec<Case> {
    vec![
        Case::new(Seeded::generated(), None),
        Case::new(Seeded::generated(), Some(Mechanism::Jacobians)),
    ]
}

/// Whether `mechanism` fires on `rows` (every file it reads has all its strata, all failing), and
/// what it saw.
fn observe(mechanism: Mechanism, rows: &[Row]) -> (bool, String) {
    let mut files: BTreeMap<&str, Vec<&Row>> = BTreeMap::new();
    for r in rows.iter().filter(|r| mechanism.reads(r)) {
        files.entry(&r.fn_id).or_default().push(r);
    }
    if files.is_empty() {
        return (false, "the subject answers none of them".to_string());
    }
    let read: Vec<&&Row> = files.values().flatten().collect();
    let failing = read
        .iter()
        .filter(|r| r.nonfinite > 0 || r.max_u >= FAILS_AT)
        .count();
    let scored = read.iter().filter(|r| r.nonfinite == 0).map(|r| r.max_u);
    let (least, most) = scored.fold((f64::INFINITY, 0.0), |(l, m), u| (l.min(u), f64::max(m, u)));
    let short: Vec<String> = files
        .iter()
        .filter(|(_, strata)| strata.len() != mechanism.cells())
        .map(|(fn_id, strata)| format!("{fn_id} has {}", strata.len()))
        .collect();
    let mut seen = format!(
        "{failing} of {} strata fail over {} file(s), {least:.1e} to {most:.1e} u",
        read.len(),
        files.len()
    );
    if !short.is_empty() {
        seen += &format!("; {} (a file has {})", short.join(", "), mechanism.cells());
    }
    (short.is_empty() && failing == read.len(), seen)
}

/// Runs `cases` over the `sen3_{exp,jr,jl}` ids of `dir`, the correct kernel against `bar`.
pub(super) fn check(dir: &Path, cases: Vec<Case>, bar: f64) -> Result<Report, String> {
    let mut entries = corpus::manifest(dir)?;
    entries.retain(|e| se3::parse(&e.fn_id).is_some());
    let (subjects, expected): (Vec<Registered>, Vec<Option<Mechanism>>) =
        cases.into_iter().map(|c| (c.subject, c.must_fire)).unzip();
    if !expected.contains(&None) {
        return Err("no correct kernel among the subjects".into());
    }
    let per_subject = evaluate(dir, &entries, &subjects, Precision::F64)?;

    let (mut text, mut failures, mut verdicts) = (String::new(), Vec::new(), String::new());
    for ((s, rows), must_fire) in subjects.iter().zip(&per_subject).zip(expected) {
        let name = s.subject.name();
        if must_fire.is_none() {
            let _ = writeln!(
                text,
                "{name} {}: the max over strata, in u (bar {bar})",
                s.version
            );
            for e in &entries {
                match worst(rows, &e.fn_id) {
                    Some(r) => {
                        let _ = writeln!(text, "  {}  {:.3e} ({})", e.fn_id, r.max_u, r.stratum);
                    }
                    None => failures.push(format!("{name}: no score for {}", e.fn_id)),
                }
            }
            failures.extend(over_the_bar(name, rows, bar));
        }
        for m in Mechanism::ALL {
            let (fired, seen) = observe(m, rows);
            let (expects, verdict, ok) = match must_fire {
                None => ("silent", if fired { "FAIL" } else { "ok" }, !fired),
                Some(f) if f == m => ("detected", if fired { "ok" } else { "FAIL" }, fired),
                Some(_) => ("-", if fired { "fires" } else { "quiet" }, true),
            };
            let _ = writeln!(
                verdicts,
                "  {name:<32} -> {:<40} -> {expects:<8} {verdict:<5} {seen}",
                m.label()
            );
            if !ok {
                failures.push(format!("{name}: `{}` is not {expects}: {seen}", m.label()));
            }
        }
    }
    let _ = writeln!(
        text,
        "\nsubject -> mechanism -> verdict (`detected`: must fire; `silent`: must not; `-`: the other \
         defect's mechanism, printed and not gated):\n{verdicts}"
    );
    Ok(Report { text, failures })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;

    fn run(cases: Vec<Case>, bar: f64) -> Result<Report, String> {
        check(&corpus_dir()?, cases, bar)
    }

    fn has(report: &Report, needle: &str) -> bool {
        let squeeze = |t: &str| t.split_whitespace().collect::<Vec<_>>().join(" ");
        squeeze(&report.text).contains(needle)
    }

    fn row(fn_id: &str, stratum: &str, max_u: f64, nonfinite: usize) -> Row {
        Row {
            fn_id: fn_id.to_string(),
            stratum: stratum.to_string(),
            precision: Precision::F64,
            subject: "s".to_string(),
            version: "v".to_string(),
            n: 1,
            max_u,
            p99_u: max_u,
            argmax_id: 0,
            nonfinite,
        }
    }

    /// The strata of `fn_id`: `rho` cells and `theta` strata, each scored `u`.
    fn file(fn_id: &str, rho: usize, theta: usize, u: f64) -> Vec<Row> {
        let strata = (0..rho)
            .map(|i| format!("rho:{i}"))
            .chain((0..theta).map(|i| format!("theta:{i}")));
        strata.map(|s| row(fn_id, &s, u, 0)).collect()
    }

    #[test]
    fn the_committed_corpus_passes_and_the_report_says_how() -> Result<(), String> {
        let report = run(cases(), BAR)?;
        assert_eq!(report.failures, Vec::<String>::new());
        // The correct kernel's max per id, in u (bit-deterministic, D16), and each defect's own
        // mechanism; the other defect's is not read.
        for line in [
            "sen3_exp_n1 3.468e0 (rho:1e3/theta=pi-1e-6)",
            "sen3_jl_n1 8.903e0 (rho:1e4/theta=pi-1e-6)",
            "seeded:correct -> every rho:* stratum of sen3_exp fails -> silent ok 0 of 75 strata fail \
             over 3 file(s), 4.5e-1 to 3.5e0 u",
            "seeded:se3-exp-translation-first -> every rho:* stratum of sen3_exp fails -> detected ok \
             75 of 75 strata fail over 3 file(s), 8.9e15 to 2.8e22 u",
            "seeded:q-minus-half -> every stratum of sen3_jr, sen3_jl fails -> detected ok 312 of 312 \
             strata fail over 6 file(s), 5.2e9 to 2.4e16 u",
            "seeded:q-minus-half -> every rho:* stratum of sen3_exp fails -> - quiet the subject \
             answers none of them",
        ] {
            assert!(has(&report, line), "{line}\n{}", report.text);
        }
        let ids = report.text.lines().filter(|l| l.starts_with("  sen3_"));
        assert_eq!(ids.count(), 9);
        assert_eq!((BAR, FAILS_AT, RHO_CELLS, STRATA), (12.0, 1e7, 25, 52));
        Ok(())
    }

    #[test]
    fn a_defect_that_is_not_planted_is_missed_and_one_in_the_correct_slot_is_loud(
    ) -> Result<(), String> {
        let missed = run(
            vec![
                Case::new(Seeded::generated(), None),
                Case::new(Seeded::generated(), Some(Mechanism::RhoStrata)),
                Case::new(Seeded::generated(), Some(Mechanism::Jacobians)),
            ],
            BAR,
        )?;
        assert_eq!(missed.failures.len(), 2, "{:?}", missed.failures);
        assert!(missed
            .failures
            .iter()
            .all(|f| f.contains("is not detected")));
        for d in Defect::SE3 {
            let loud = run(vec![Case::new(Seeded::planted(d), None)], BAR)?;
            let name = format!("seeded:{}", d.name());
            assert!(loud.failures.iter().any(|f| f.starts_with(&name)), "{d:?}");
        }
        let alone = run(
            vec![Case::new(Seeded::planted(Defect::QMinusHalf), None)],
            BAR,
        );
        assert!(alone.is_ok());
        let none = run(
            vec![Case::new(
                Seeded::planted(Defect::QMinusHalf),
                Some(Mechanism::Jacobians),
            )],
            BAR,
        );
        assert!(none.err().unwrap_or_default().contains("no correct kernel"));
        Ok(())
    }

    #[test]
    fn a_bar_under_the_measured_maximum_fails_the_correct_kernel_on_every_id() -> Result<(), String>
    {
        let report = run(cases(), 2.0)?;
        for op in ["exp", "jr", "jl"] {
            for n in 1..=3 {
                let id = format!("seeded:correct: sen3_{op}_n{n} ");
                let over = |f: &String| f.starts_with(&id) && f.contains("over the bar of 2 u");
                assert!(report.failures.iter().any(over), "{id}");
            }
        }
        assert!(report
            .failures
            .iter()
            .all(|f| f.starts_with("seeded:correct")));
        Ok(())
    }

    #[test]
    fn a_mechanism_fires_only_where_every_stratum_of_a_whole_file_fails() {
        let rho = |rows: Vec<Row>| observe(Mechanism::RhoStrata, &rows).0;
        // The `theta:*` strata of the file pass and are not read; another file's rows are not either.
        let all = || {
            let mut rows = file("sen3_exp_n1", RHO_CELLS, 0, 1e9);
            rows.extend(file("sen3_exp_n1", 0, 27, 1.0));
            rows
        };
        assert!(rho(all()));
        let mut quiet = all();
        quiet.extend(file("sen3_exp_n2", RHO_CELLS, 0, 1.0));
        assert!(!rho(quiet));
        // One stratum under the line, one missing, a `NaN`, and a non-finite output.
        let mut under = all();
        under[3].max_u = 9.9e6;
        assert!(!rho(under));
        let mut short = all();
        short.remove(0);
        assert!(!rho(short));
        let mut nan = all();
        (nan[0].max_u, nan[0].nonfinite) = (f64::NAN, 0);
        assert!(!rho(nan));
        let mut lost = all();
        (lost[0].max_u, lost[0].nonfinite) = (f64::NAN, 2);
        assert!(rho(lost));
        // `Jacobians` reads both sides and every stratum, and only those.
        let jac = |rows: Vec<Row>| observe(Mechanism::Jacobians, &rows);
        let mut rows = file("sen3_jr_n1", 25, 27, 1e9);
        rows.extend(file("sen3_jl_n3", 25, 27, 1e9));
        rows.extend(file("sen3_exp_n1", 25, 27, 1.0));
        assert!(jac(rows.clone()).0);
        rows[60].max_u = 1e6;
        let (fired, seen) = jac(rows);
        assert!(
            !fired && seen.starts_with("103 of 104 strata fail over 2 file(s)"),
            "{seen}"
        );
        assert!(!jac(Vec::new()).0);
    }
}
