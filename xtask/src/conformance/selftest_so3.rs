//! The SO(3) half of `--self-test` (`docs/PHASE1.md` §10): the correct seeded `Exp` and `Log` must
//! stay under [`BAR`] on every stratum of `so3_exp` and `so3_log` and fire no mechanism, and each
//! planted `Log` defect must fire its own, over the committed corpus.
//!
//! Readings where §10 is silent, each the smallest (0014 (draft) questions 16 to 18):
//!
//! - **The bar** is [`BAR`] `u` per stratum, on the max: a few roundings (`docs/maths/so3.md`
//!   SO.6 measures `2 atan2(n, w)` at 2.2 `u`, `docs/maths/coefficients.md` CO.16 the `r` exact arm
//!   at 1.9 `u`), 1.4 times the measured 2.91 `u`. The spec states none for the correct kernel.
//! - **`acos`** is detected when every stratum of `theta:1e-k` from `k =` [`ACOS_FROM`]`.0` and of
//!   `theta:pi-1e-k` from `k =` `.1` reaches [`FAILS_AT`] `u`, §10's `10^7`. §10 cannot mean every
//!   stratum: `θ` read from `χ = cos θ` errs by `≈ u/(2θ²)` (`docs/maths/so3.md` SO.6(a)), which
//!   is `10^7 u` only below `θ ≈ 2·10^-4` and, near `π`, from `π − θ ≈ 10^-7` (the diagonal form
//!   loses more than `χ` alone). The two `k` are those boundaries, measured (8.2e7 `u` at
//!   `theta:1e-4`, 1.5e7 at `theta:pi-1e-7`); the family maxima are printed.
//! - **The negated half** is the records whose input has `w < 0`: the corpus holds each quaternion
//!   and then its negative in every `so3_log` stratum but `q:w0` (`conformance/generate/README.md`).
//!   A stratum's negated half fails when its max is at least [`FAILS_AT`] `u`; the missing flip is
//!   `2π` in `φ`, so at least `2/u`. A stratum has nothing to detect where every negated record
//!   takes the series arm of `r` (`seeded::takes_series_arm`, `n² = 0`: `theta:exact0`,
//!   `theta:subnormal`), because `2/w` is analytic in the sign of `w`, and the defect answers it
//!   bit for bit as the correct kernel. §10's "every `so3_log` stratum" holds for the rest.
//!   Detected when every stratum fails or has nothing to detect, and one fails.
//! - **A gap is a failure**: a stratum of §10's two families with no score, a non-finite output
//!   and a `NaN` are read as failing the bar, never as passing it.

use std::fmt::Write;
use std::path::Path;

use helicoid_linalg::Precision;

use super::corpus::{self, Record};
use super::evaluate;
use super::metric::{self, Score};
use super::report::Row;
use super::selftest::{worst, Report};
use super::subject::{Output, Registered};
use crate::seeded::{takes_series_arm, Defect, Seeded};

/// The largest `max_u` the correct kernel may show on a stratum.
pub(super) const BAR: f64 = 4.0;
/// §10's `10^7 u` for `acos`, and the negated half's failure line.
pub(super) const FAILS_AT: f64 = 1e7;
/// The strata `theta:1e-k` and `theta:pi-1e-k` of §10, `k = 1, ..., 12`.
const DECADES: std::ops::RangeInclusive<u32> = 1..=12;
/// The two families of §10, by the prefix of their stratum names.
const STEMS: [&str; 2] = ["theta:", "theta:pi-"];
/// The `k` from which every stratum of each family must reach [`FAILS_AT`] under `acos`.
const ACOS_FROM: (usize, usize) = (4, 7);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mechanism {
    AcosLoss,
    NegatedHalf,
}

impl Mechanism {
    const ALL: [Mechanism; 2] = [Mechanism::AcosLoss, Mechanism::NegatedHalf];

    fn label(self) -> &'static str {
        match self {
            Mechanism::AcosLoss => "theta:1e-k, pi-1e-k >= 1e7 u",
            Mechanism::NegatedHalf => "negated half of so3_log fails",
        }
    }

    fn of(defect: Defect) -> Option<Self> {
        match defect {
            Defect::LogAcos => Some(Mechanism::AcosLoss),
            Defect::LogNoFlip => Some(Mechanism::NegatedHalf),
            _ => None,
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
    for d in Defect::LOG {
        all.push(Case::new(Seeded::planted(d), Mechanism::of(d)));
    }
    all
}

/// One stratum's negated half of one subject.
struct Half {
    stratum: String,
    /// The largest score over the negated records; infinite if one is not finite.
    worst: f64,
    /// Every negated record is answered bit for bit as the correct kernel answers it.
    identical: bool,
    /// Every negated record takes the series arm of `r`, where the flip changes nothing.
    hidden: bool,
}

impl Half {
    /// A stratum where the missing flip cannot show, and the subject answers as the correct kernel.
    fn nothing_to_detect(&self) -> bool {
        self.hidden && self.identical && self.worst < FAILS_AT
    }
}

fn same_bits(a: &Output, b: &Output) -> bool {
    let bits = |o: &Output| -> Vec<(String, Vec<u64>)> {
        let field =
            |(k, v): (&String, &Vec<f64>)| (k.clone(), v.iter().map(|x| x.to_bits()).collect());
        o.iter().map(field).collect()
    };
    bits(a) == bits(b)
}

/// The negated half of every stratum that has one, in corpus order.
fn halves(records: &[Record], out: &[Output], correct: &[Output]) -> Result<Vec<Half>, String> {
    let rule = metric::rule("so3_log").ok_or("no metric rule for `so3_log`")?;
    let mut all: Vec<Half> = Vec::new();
    for ((rec, o), c) in records.iter().zip(out).zip(correct) {
        let Some(&q) = rec.input("q").and_then(|q| q.first_chunk::<4>()) else {
            continue;
        };
        if q[0] >= 0.0 {
            continue;
        }
        let score = match rule.score(rec, o, Precision::F64)? {
            Score::Finite(f) => f,
            _ => f64::INFINITY,
        };
        let at = match all.iter().position(|h| h.stratum == rec.stratum) {
            Some(at) => at,
            None => {
                all.push(Half {
                    stratum: rec.stratum.clone(),
                    worst: 0.0,
                    identical: true,
                    hidden: true,
                });
                all.len() - 1
            }
        };
        let h = &mut all[at];
        h.worst = h.worst.max(score);
        h.identical &= same_bits(o, c);
        h.hidden &= takes_series_arm(q);
    }
    Ok(all)
}

/// The `max_u` of `so3_log` at `<stem>1e-k` for `k = 1, ..., 12`: infinite where a record is not
/// finite, `None` where the stratum has no row or no score.
fn tail(rows: &[Row], stem: &str) -> Vec<Option<f64>> {
    let at = |k: u32| {
        let name = format!("{stem}1e-{k}");
        let row = rows
            .iter()
            .find(|r| r.fn_id == "so3_log" && r.stratum == name)?;
        match (row.nonfinite, row.max_u) {
            (0, u) if u.is_nan() => None,
            (0, u) => Some(u),
            _ => Some(f64::INFINITY),
        }
    };
    DECADES.map(at).collect()
}

/// The strata of §10's two families with no score in `rows`.
fn gaps(rows: &[Row]) -> Vec<String> {
    let mut all = Vec::new();
    for stem in STEMS {
        let missing = DECADES.zip(tail(rows, stem)).filter(|(_, u)| u.is_none());
        all.extend(missing.map(|(k, _)| format!("{stem}1e-{k}")));
    }
    all
}

/// The `k` from which every stratum of `tail` reaches [`FAILS_AT`], if the last does.
fn reaches_from(tail: &[Option<f64>]) -> Option<usize> {
    let reaches = |u: &&Option<f64>| u.is_some_and(|u| u >= FAILS_AT);
    let run = tail.iter().rev().take_while(reaches).count();
    (run > 0).then(|| tail.len() + 1 - run)
}

fn observe(mechanism: Mechanism, rows: &[Row], negated: &[Half]) -> (bool, String) {
    match mechanism {
        Mechanism::AcosLoss => {
            let (small, near_pi) = (tail(rows, STEMS[0]), tail(rows, STEMS[1]));
            let peak = |t: &[Option<f64>]| t.iter().flatten().copied().fold(0.0, f64::max);
            let from =
                |t: &[Option<f64>]| reaches_from(t).map_or("none".to_string(), |k| k.to_string());
            let seen = format!(
                "theta:1e-k max {:.1e} u, every stratum from k = {} (gate {}); theta:pi-1e-k max \
                 {:.1e} u, from k = {} (gate {})",
                peak(&small),
                from(&small),
                ACOS_FROM.0,
                peak(&near_pi),
                from(&near_pi),
                ACOS_FROM.1
            );
            let by = |t: &[Option<f64>], gate: usize| reaches_from(t).is_some_and(|k| k <= gate);
            (by(&small, ACOS_FROM.0) && by(&near_pi, ACOS_FROM.1), seen)
        }
        Mechanism::NegatedHalf => {
            let fails = negated.iter().filter(|h| h.worst >= FAILS_AT).count();
            let idle: Vec<&str> = negated
                .iter()
                .filter(|h| h.nothing_to_detect())
                .map(|h| h.stratum.as_str())
                .collect();
            let undecided = negated.len() - fails - idle.len();
            let hidden = match idle.len() {
                0 => "none".to_string(),
                n if n == negated.len() => "all".to_string(),
                _ => idle.join(", "),
            };
            let seen = format!(
                "{fails} of {} strata fail, {undecided} undecided (nothing to detect where n² = 0 \
                 and answered as the correct kernel: {hidden})",
                negated.len()
            );
            (fails > 0 && undecided == 0, seen)
        }
    }
}

/// The rows that fail `bar`: a non-finite output, a `NaN` and a `max_u` above it.
pub(super) fn over_the_bar(name: &str, rows: &[Row], bar: f64) -> Vec<String> {
    let mut failures = Vec::new();
    for r in rows {
        if r.nonfinite > 0 {
            failures.push(format!(
                "{name}: {} {} has {} non-finite outputs",
                r.fn_id, r.stratum, r.nonfinite
            ));
        } else if r.max_u.is_nan() || r.max_u > bar {
            failures.push(format!(
                "{name}: {} {} is {} u, over the bar of {bar} u",
                r.fn_id, r.stratum, r.max_u
            ));
        }
    }
    failures
}

/// Runs `cases` over `so3_exp` and `so3_log` of `dir`, the correct kernel against `bar`.
pub(super) fn check(dir: &Path, cases: Vec<Case>, bar: f64) -> Result<Report, String> {
    let mut entries = corpus::manifest(dir)?;
    let ids = ["so3_exp", "so3_log"];
    entries.retain(|e| ids.contains(&e.fn_id.as_str()));
    let (subjects, expected): (Vec<Registered>, Vec<Option<Mechanism>>) =
        cases.into_iter().map(|c| (c.subject, c.must_fire)).unzip();
    let correct = expected
        .iter()
        .position(Option::is_none)
        .ok_or("no correct kernel among the subjects")?;
    let per_subject = evaluate(dir, &entries, &subjects, Precision::F64)?;
    let log = entries
        .iter()
        .find(|e| e.fn_id == "so3_log")
        .ok_or("no corpus file for `so3_log`")?;
    let records = corpus::read(dir, log)?;
    let answers = |s: &Registered| -> Vec<Output> {
        let eval = |r: &Record| s.subject.eval("so3_log", r, Precision::F64);
        records.iter().map(eval).collect()
    };
    let reference = answers(&subjects[correct]);

    let (mut text, mut failures) = (String::new(), Vec::new());
    let mut verdicts = String::new();
    for ((s, rows), must_fire) in subjects.iter().zip(&per_subject).zip(expected) {
        let name = s.subject.name();
        let unscored = gaps(rows);
        if !unscored.is_empty() {
            failures.push(format!("{name}: no score in {}", unscored.join(", ")));
        }
        if must_fire.is_none() {
            let _ = writeln!(
                text,
                "{name} {}: the max over strata, in u (bar {bar})",
                s.version
            );
            for fn_id in ids {
                let cell = worst(rows, fn_id).map_or("-".to_string(), |r| {
                    format!("{:.3e} ({})", r.max_u, r.stratum)
                });
                let _ = writeln!(text, "  {fn_id}  {cell}");
                if worst(rows, fn_id).is_none() {
                    failures.push(format!("{name}: no score for {fn_id}"));
                }
            }
            failures.extend(over_the_bar(name, rows, bar));
        }
        let negated = if s.subject.supports("so3_log") {
            halves(&records, &answers(s), &reference)?
        } else {
            Vec::new()
        };
        for m in Mechanism::ALL {
            let (fired, seen) = observe(m, rows, &negated);
            let (expects, verdict, ok) = match must_fire {
                None => ("silent", if fired { "FAIL" } else { "ok" }, !fired),
                Some(f) if f == m => ("detected", if fired { "ok" } else { "FAIL" }, fired),
                Some(_) => ("-", if fired { "fires" } else { "quiet" }, true),
            };
            let _ = writeln!(
                verdicts,
                "  {name:<20} -> {:<30} -> {expects:<8} {verdict:<5} {seen}",
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
    use crate::conformance::subject::Subject;
    use crate::conformance::testkit::Fixed;

    fn run(cases: Vec<Case>, bar: f64) -> Result<Report, String> {
        check(&corpus_dir()?, cases, bar)
    }

    fn correct() -> Case {
        Case::new(Seeded::generated(), None)
    }

    fn planted(d: Defect, must_fire: Option<Mechanism>) -> Case {
        Case::new(Seeded::planted(d), must_fire)
    }

    fn has(report: &Report, needle: &str) -> bool {
        report
            .text
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(" ")
            .contains(needle)
    }

    fn so3_log_records() -> Result<Vec<Record>, String> {
        let dir = corpus_dir()?;
        let entry = corpus::manifest(&dir)?
            .into_iter()
            .find(|e| e.fn_id == "so3_log")
            .ok_or("no so3_log")?;
        corpus::read(&dir, &entry)
    }

    fn ask(records: &[Record], s: Seeded) -> Vec<Output> {
        let s = s.registered();
        let eval = |r: &Record| s.subject.eval("so3_log", r, Precision::F64);
        records.iter().map(eval).collect()
    }

    fn negated(r: &Record) -> bool {
        r.input("q").is_some_and(|q| q[0] < 0.0)
    }

    #[test]
    fn the_committed_corpus_passes_and_the_report_says_how() -> Result<(), String> {
        let report = run(cases(), BAR)?;
        assert_eq!(report.failures, Vec::<String>::new());
        // The correct kernel's max over strata, in u, per id (bit-deterministic, D16).
        for line in [
            "so3_exp 2.743e0 (theta:pi-1e-9)",
            "so3_log 2.901e0 (theta:dense)",
            "seeded:log-acos -> theta:1e-k, pi-1e-k >= 1e7 u -> detected ok theta:1e-k max 9.0e15 u, \
             every stratum from k = 4 (gate 4); theta:pi-1e-k max 8.5e7 u, from k = 7 (gate 7)",
            "seeded:log-no-flip -> negated half of so3_log fails -> detected ok 27 of 29 strata fail, \
             0 undecided (nothing to detect where n² = 0 and answered as the correct kernel: \
             theta:exact0, theta:subnormal)",
            // Each defect also fires the other's line where its errors are that large, and not where not.
            "seeded:log-no-flip -> theta:1e-k, pi-1e-k >= 1e7 u -> - fires",
            "seeded:log-acos -> negated half of so3_log fails -> - quiet 17 of 29 strata fail, 12 undecided",
        ] {
            assert!(has(&report, line), "{line}\n{}", report.text);
        }
        assert_eq!(
            (BAR, FAILS_AT, DECADES, ACOS_FROM),
            (4.0, 1e7, 1..=12, (4, 7))
        );
        Ok(())
    }

    #[test]
    fn the_flip_is_all_the_defect_removes() -> Result<(), String> {
        // Where w >= 0 the defect answers as the correct kernel does, to the bit; where w < 0 it
        // differs on every stratum but `theta:exact0` and `theta:subnormal`, and `q:w0` has no such
        // record.
        let records = so3_log_records()?;
        let (good, bad) = (
            ask(&records, Seeded::generated()),
            ask(&records, Seeded::planted(Defect::LogNoFlip)),
        );
        let mut differing: Vec<&str> = Vec::new();
        for ((r, g), b) in records.iter().zip(&good).zip(&bad) {
            assert_eq!(
                same_bits(g, b),
                !negated(r) || r.stratum == "theta:exact0" || r.stratum == "theta:subnormal",
                "{}",
                r.id
            );
            if !same_bits(g, b) && !differing.contains(&r.stratum.as_str()) {
                differing.push(&r.stratum);
            }
        }
        assert_eq!(differing.len(), 27);
        assert!(!differing.contains(&"q:w0"));
        let halves = halves(&records, &bad, &good)?;
        assert_eq!(halves.len(), 29);
        assert!(halves.iter().all(|h| h.stratum != "q:w0"));
        // Only those two strata take the series arm of `r` for every negated record, derived from
        // the records and not listed.
        let hidden: Vec<&str> = halves
            .iter()
            .filter(|h| h.hidden)
            .map(|h| h.stratum.as_str())
            .collect();
        assert_eq!(hidden, ["theta:exact0", "theta:subnormal"]);
        assert!(halves
            .iter()
            .all(|h| h.nothing_to_detect() == h.hidden && h.identical == h.hidden));
        Ok(())
    }

    #[test]
    fn a_defect_that_is_not_planted_is_missed_and_one_in_the_correct_slot_is_loud(
    ) -> Result<(), String> {
        let missed = run(
            vec![
                correct(),
                Case::new(Seeded::generated(), Some(Mechanism::AcosLoss)),
                Case::new(Seeded::generated(), Some(Mechanism::NegatedHalf)),
            ],
            BAR,
        )?;
        assert_eq!(missed.failures.len(), 2, "{:?}", missed.failures);
        assert!(missed
            .failures
            .iter()
            .all(|f| f.contains("is not detected")));
        for d in Defect::LOG {
            let loud = run(vec![planted(d, None)], BAR)?;
            let name = format!("seeded:{}", d.name());
            assert!(loud.failures.iter().any(|f| f.starts_with(&name)), "{d:?}");
        }
        let none = run(
            vec![planted(Defect::LogAcos, Mechanism::of(Defect::LogAcos))],
            BAR,
        );
        assert!(none.err().unwrap_or_default().contains("no correct kernel"));
        Ok(())
    }

    #[test]
    fn a_bar_below_the_measured_maximum_fails_the_correct_kernel_on_each_id() -> Result<(), String>
    {
        let report = run(cases(), 2.5)?;
        // Each id's own maximum (2.743, 2.901) is over 2.5, so the bar is held on both.
        for id in ["so3_exp", "so3_log"] {
            let over = |f: &String| {
                f.starts_with(&format!("seeded:correct: {id} "))
                    && f.contains("over the bar of 2.5 u")
            };
            assert!(
                report.failures.iter().any(over),
                "{id}: {:?}",
                report.failures
            );
        }
        // The defects are not held to the bar.
        assert!(report
            .failures
            .iter()
            .all(|f| f.starts_with("seeded:correct")));
        Ok(())
    }

    /// The correct kernel, answering `NaN` in the first entry of every record `poison` names, and
    /// only the id `only` if there is one.
    struct Poisoned {
        inner: Registered,
        poison: fn(&Record) -> bool,
        only: Option<&'static str>,
    }

    impl Subject for Poisoned {
        fn name(&self) -> &str {
            "poisoned"
        }

        fn supports(&self, fn_id: &str) -> bool {
            self.inner.subject.supports(fn_id) && self.only.is_none_or(|o| o == fn_id)
        }

        fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
            let mut out = self.inner.subject.eval(fn_id, record, precision);
            if (self.poison)(record) {
                out.values_mut().for_each(|v| v[0] = f64::NAN);
            }
            out
        }
    }

    #[test]
    fn a_non_finite_answer_fails_the_correct_kernel_whatever_the_bar() -> Result<(), String> {
        let poisoned = Poisoned {
            inner: Seeded::generated().registered(),
            poison: |r| r.stratum == "theta:dense",
            only: None,
        };
        let case = Case {
            subject: Registered::new("v", Box::new(poisoned)),
            must_fire: None,
        };
        let report = run(vec![case], 1e300)?;
        let hit = |f: &String| f.starts_with("poisoned: so3_log theta:dense has ");
        assert!(report.failures.iter().any(hit), "{:?}", report.failures);
        assert!(report
            .failures
            .iter()
            .any(|f| f.starts_with("poisoned: so3_exp")));
        Ok(())
    }

    #[test]
    fn a_correct_kernel_that_answers_only_one_id_is_not_silent() -> Result<(), String> {
        let case = |only| Case {
            subject: Registered::new(
                "v",
                Box::new(Poisoned {
                    inner: Seeded::generated().registered(),
                    poison: |_| false,
                    only: Some(only),
                }),
            ),
            must_fire: None,
        };
        let no_exp = run(vec![case("so3_log")], BAR)?.failures;
        assert!(
            no_exp.contains(&"poisoned: no score for so3_exp".to_string()),
            "{no_exp:?}"
        );
        let no_log = run(vec![case("so3_exp")], BAR)?.failures;
        assert!(
            no_log.contains(&"poisoned: no score for so3_log".to_string()),
            "{no_log:?}"
        );
        let gap = "poisoned: no score in theta:1e-1, theta:1e-2, theta:1e-3";
        assert!(no_log.iter().any(|f| f.starts_with(gap)), "{no_log:?}");
        Ok(())
    }

    #[test]
    fn the_bar_fails_a_non_finite_output_a_nan_and_a_max_above_it_and_passes_the_bar_itself() {
        let rows = [
            row("a", 1.0, 0),
            row("b", f64::NAN, 0),
            row("c", 5.0, 0),
            row("d", 1.0, 2),
            row("e", 4.0, 0),
        ];
        let failed = over_the_bar("s", &rows, 4.0);
        let names: Vec<&str> = failed
            .iter()
            .filter_map(|f| f.split_whitespace().nth(2))
            .collect();
        assert_eq!(names, ["b", "c", "d"], "{failed:?}");
        assert!(failed[2].ends_with("has 2 non-finite outputs"));
    }

    #[test]
    fn a_non_finite_negated_record_is_a_failure_and_not_a_score_of_zero() -> Result<(), String> {
        let records = so3_log_records()?;
        let good = ask(&records, Seeded::generated());
        let mut nan = good.clone();
        for (r, o) in records.iter().zip(&mut nan) {
            if negated(r) && r.stratum == "theta:dense" {
                o.values_mut().for_each(|v| v[0] = f64::NAN);
            }
        }
        let halves = halves(&records, &nan, &good)?;
        let dense = halves.iter().find(|h| h.stratum == "theta:dense");
        let dense = dense.ok_or("no theta:dense half")?;
        assert!(dense.worst.is_infinite() && !dense.identical);
        Ok(())
    }

    #[test]
    fn a_defect_that_flips_only_where_it_is_easy_is_not_detected() -> Result<(), String> {
        // The correct `Log` where `n² < 0.01` and the flip missing elsewhere: it answers most
        // strata as the correct kernel does without `n² = 0` to excuse it, so the gate must not.
        let (good, bad) = (
            Seeded::generated().registered(),
            Seeded::planted(Defect::LogNoFlip).registered(),
        );
        let near = Fixed::new("flips-near", move |r| {
            let n2: f64 = r
                .input("q")
                .map_or(0.0, |q| q[1..].iter().map(|c| c * c).sum());
            let s = if n2 < 0.01 { &good } else { &bad };
            s.subject.eval("so3_log", r, Precision::F64)
        })
        .only("so3_log");
        let case = Case {
            subject: Registered::new("v", Box::new(near)),
            must_fire: Some(Mechanism::NegatedHalf),
        };
        let report = run(vec![correct(), case], BAR)?;
        assert_eq!(report.failures.len(), 1, "{:?}", report.failures);
        let f = &report.failures[0];
        assert!(
            f.starts_with("flips-near: `negated half of so3_log fails` is not detected"),
            "{f}"
        );
        Ok(())
    }

    fn row(stratum: &str, max_u: f64, nonfinite: usize) -> Row {
        Row {
            fn_id: "so3_log".to_string(),
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

    /// The 24 strata of §10's two families, each scored `(max_u, nonfinite)` by its `k`.
    fn family(
        small: impl Fn(u32) -> (f64, usize),
        near_pi: impl Fn(u32) -> (f64, usize),
    ) -> Vec<Row> {
        let mut rows = Vec::new();
        for k in DECADES {
            let ((a, an), (b, bn)) = (small(k), near_pi(k));
            rows.push(row(&format!("theta:1e-{k}"), a, an));
            rows.push(row(&format!("theta:pi-1e-{k}"), b, bn));
        }
        rows
    }

    /// `1e8` from `k = from`, `1` before.
    fn from(from: u32) -> impl Fn(u32) -> (f64, usize) {
        move |k| (if k >= from { 1e8 } else { 1.0 }, 0)
    }

    fn acos_fires(rows: &[Row]) -> bool {
        observe(Mechanism::AcosLoss, rows, &[]).0
    }

    #[test]
    fn acos_needs_both_families_from_their_own_k_and_reads_a_non_finite_as_failing() {
        assert!(acos_fires(&family(from(4), from(7))));
        assert!(acos_fires(&family(from(1), from(1))));
        // One family alone, either one.
        assert!(!acos_fires(&family(from(4), from(99))));
        assert!(!acos_fires(&family(from(99), from(7))));
        // A decade late, in either family, and a single stratum.
        assert!(!acos_fires(&family(from(5), from(7))));
        assert!(!acos_fires(&family(from(4), from(8))));
        let last = |k| (if k == 12 { 1e8 } else { 1.0 }, 0);
        assert!(!acos_fires(&family(last, last)));
        // A non-finite record is an infinite error, whatever `max_u` says.
        let lost = |k| (1.0, usize::from(k >= 4));
        let lost_pi = |k| (1.0, usize::from(k >= 7));
        assert!(acos_fires(&family(lost, lost_pi)));
    }

    #[test]
    fn a_stratum_with_no_score_is_a_gap_not_a_narrower_family() {
        assert_eq!(gaps(&family(from(4), from(7))), Vec::<String>::new());
        let mut rows = family(from(4), from(7));
        rows.retain(|r| r.stratum != "theta:pi-1e-5");
        rows.iter_mut()
            .filter(|r| r.stratum == "theta:1e-9")
            .for_each(|r| r.max_u = f64::NAN);
        assert_eq!(gaps(&rows), ["theta:1e-9", "theta:pi-1e-5"]);
        // A stratum whose records are all non-finite has no `max_u` and is scored, infinite.
        let mut rows = family(from(4), from(7));
        rows.iter_mut()
            .filter(|r| r.stratum == "theta:1e-9")
            .for_each(|r| (r.max_u, r.nonfinite) = (f64::NAN, 3));
        assert_eq!(gaps(&rows), Vec::<String>::new());
    }

    #[test]
    fn a_stratum_is_excused_only_where_the_flip_cannot_show_and_the_subject_is_the_correct_one() {
        let half = |worst, identical, hidden| Half {
            stratum: "s".to_string(),
            worst,
            identical,
            hidden,
        };
        let fires = |negated: &[Half]| observe(Mechanism::NegatedHalf, &[], negated).0;
        let fails = || half(1e16, false, false);
        assert!(fires(&[fails(), half(0.5, true, true)]));
        assert!(fires(&[fails(), fails()]));
        // Answered as the correct kernel where the flip shows, or hidden but answered otherwise.
        assert!(!fires(&[fails(), half(0.5, true, false)]));
        assert!(!fires(&[fails(), half(0.5, false, true)]));
        // Under the line, and nothing failing.
        assert!(!fires(&[fails(), half(5.0, false, false)]));
        assert!(!fires(&[half(0.5, true, true)]));
        assert!(!fires(&[]));
    }

    #[test]
    fn the_acos_family_is_read_from_the_end_and_a_gap_stops_the_run() {
        let some = |t: [f64; 12]| t.map(Some);
        let mut t = some([1.0; 12]);
        assert_eq!(reaches_from(&t), None);
        t[11] = Some(1e8);
        assert_eq!(reaches_from(&t), Some(12));
        t[10] = Some(1e7);
        t[8] = Some(1e9);
        assert_eq!(reaches_from(&t), Some(11));
        assert_eq!(reaches_from(&some([1e7; 12])), Some(1));
        // A stratum with no score does not reach, and a run stops there.
        t[10] = None;
        assert_eq!(reaches_from(&t), Some(12));
        assert_eq!(reaches_from(&[None; 12]), None);
    }
}
