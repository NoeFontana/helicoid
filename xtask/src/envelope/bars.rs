//! The bars (`docs/PHASE1.md` §8, `docs/maths/error-analysis.md` EA.12), per `(fn, stratum,
//! precision)` and on the max, never a mean. A pure function of rows and the baseline: nothing is
//! read or written here.
//!
//! - **Domination**: the candidate's `max_u` at most the smallest oracle `max_u`, ties passing; a
//!   stratum with no scored oracle row is unpaired, not failed. An oracle row that is non-finite
//!   or unscored has no maximum, which is `+inf`: it is never the best, so a stratum where every
//!   oracle is bad does not lower the bar (§7). Its row stays in the oracle's own CSV.
//! - **No-regress**: the candidate's `max_u` at most the baseline's, compared exactly. A baseline
//!   with no row for a scored stratum, a baseline row no scored stratum answers (a stratum is
//!   never deleted or narrowed, D7) and a different `n` (the corpus changed under a maximum over
//!   fewer or more records) each fail.
//! - **Non-finite**: any non-finite candidate output fails (`NUMERICS.md` §11, EA.11(e)).
//!
//! A candidate row with nothing scored is judged by none when its function is scored by backward
//! error only (`metric::Rule::BackwardOnly`: `so3_from_matrix`, backward error owed); on any other
//! function it fails, since nothing but a defect in the subject or the harness leaves one.

use std::collections::{BTreeMap, BTreeSet};

use crate::conformance::metric::{self, Rule};
use crate::conformance::report::Row;

use super::rows::{Baseline, Key, Subject};

/// What a failure is about.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Bar {
    Domination,
    /// A baseline that is missing, a scored row it lacks, a maximum over its, or a wider `n`: the
    /// baseline is what `--bless` rewrites, so it moves these.
    NoRegress,
    /// A baseline stratum no scored row answers, or one over more records than the candidate's: a
    /// stratum is never deleted or narrowed (D7), so `--bless` refuses these.
    Shrunk,
    NonFinite,
    /// A row with nothing scored on a function that is scored.
    Unscored,
    /// An oracle scored on another number of records.
    Corpus,
    Coverage,
    /// An in-process candidate that wrote no rows.
    Candidate,
    Drift,
}

impl Bar {
    /// The word a `FAIL` line and the docs use.
    pub(super) fn name(self) -> &'static str {
        match self {
            Bar::Domination => "domination",
            Bar::NoRegress | Bar::Shrunk => "no-regress",
            Bar::NonFinite => "non-finite",
            Bar::Unscored => "unscored",
            Bar::Corpus => "corpus",
            Bar::Coverage => "coverage",
            Bar::Candidate => "candidate",
            Bar::Drift => "drift",
        }
    }

    /// Whether a bless moves the failure into the baseline instead of refusing it.
    pub(super) fn moves_baseline(self) -> bool {
        self == Bar::NoRegress
    }
}

pub(super) struct Failure {
    pub(super) bar: Bar,
    pub(super) text: String,
}

impl Failure {
    pub(super) fn new(bar: Bar, text: String) -> Self {
        Self { bar, text }
    }
}

#[derive(Default)]
pub(super) struct Verdict {
    pub(super) failures: Vec<Failure>,
    /// Candidate rows with a maximum.
    pub(super) scored: usize,
    /// Of those, the ones with a scored oracle row, and the ones without.
    pub(super) paired: usize,
    pub(super) unpaired: usize,
    /// Candidate rows with nothing scored on a function that has no forward error yet.
    pub(super) unscored: usize,
    /// Scored rows strictly under their baseline: the baseline is stale, not the candidate wrong.
    pub(super) improved: usize,
    /// Domination failures where the winning oracle shares our transcendental backend, so the
    /// difference is in the program and not in `libm` (`0036`, draft).
    pub(super) dominated_same_backend: usize,
    /// Domination failures where only a host-`std` oracle won, which D16 may account for.
    pub(super) dominated_host_std: usize,
}

enum Status {
    Scored(f64),
    Unscored,
    NonFinite,
}

/// Whether `r` has a maximum: the rows a baseline holds.
pub(super) fn is_scored(r: &Row) -> bool {
    matches!(status(r), Status::Scored(_))
}

fn status(r: &Row) -> Status {
    match (r.nonfinite, r.max_u.is_finite()) {
        (0, true) => Status::Scored(r.max_u),
        (0, false) => Status::Unscored,
        _ => Status::NonFinite,
    }
}

/// The oracles' rows by key.
pub(super) struct Index<'a>(Vec<(&'a Subject, BTreeMap<Key, &'a Row>)>);

impl<'a> Index<'a> {
    pub(super) fn new(oracles: &'a [Subject]) -> Self {
        let by_key = |o: &'a Subject| o.rows.iter().map(|r| (Key::of(r), r)).collect();
        Self(oracles.iter().map(|o| (o, by_key(o))).collect())
    }

    /// The oracle with the smallest maximum among those that scored `key` on `n` records, the
    /// first of equal ones; `other` gets each oracle that scored it on another number, which is
    /// no comparison.
    pub(super) fn best(
        &self,
        key: &Key,
        n: usize,
        other: &mut Vec<(&'a str, usize)>,
    ) -> Option<(&'a Subject, f64)> {
        let mut best: Option<(&Subject, f64)> = None;
        for (o, rows) in &self.0 {
            let Some(theirs) = rows.get(key) else {
                continue;
            };
            if theirs.n != n {
                other.push((&o.name, theirs.n));
            } else if let Status::Scored(m) = status(theirs) {
                if best.is_none_or(|(_, b)| m < b) {
                    best = Some((o, m));
                }
            }
        }
        best
    }

    /// Whether any oracle on the `libm` crate scored `key` on `n` records and beat `max`.
    ///
    /// Not "is the best oracle on the `libm` crate": if a same-backend oracle beats the candidate at
    /// all, the gap is in the program whatever a glibc-backed oracle does, and the two questions
    /// disagree — on this corpus, 30 strata against 26 (`0036`, draft).
    pub(super) fn beaten_by_libm_crate(&self, key: &Key, n: usize, max: f64) -> bool {
        self.0.iter().any(|(o, rows)| {
            crate::conformance::oracle_backend(&o.name)
                == Some(crate::conformance::Backend::LibmCrate)
                && rows.get(key).is_some_and(|theirs| {
                    theirs.n == n && matches!(status(theirs), Status::Scored(m) if m < max)
                })
        })
    }
}

/// `candidate` against `oracles` (in the order given: the first of equal bests is named) and, for
/// no-regress, `baseline`; `None` is a baseline file that does not exist, one failure and no row
/// read against it.
pub(super) fn judge(
    candidate: &Subject,
    oracles: &[Subject],
    baseline: Option<&Baseline>,
) -> Verdict {
    let index = Index::new(oracles);
    let mut v = Verdict::default();
    let mut answered = BTreeSet::new();
    if baseline.is_none() {
        v.failures.push(Failure::new(
            Bar::NoRegress,
            format!(
                "no baseline for `{}`; `just envelope --bless` writes it",
                candidate.name
            ),
        ));
    }
    for r in &candidate.rows {
        let key = Key::of(r);
        let max = match status(r) {
            Status::Scored(max) => max,
            Status::Unscored if matches!(metric::rule(&r.fn_id), Some(Rule::BackwardOnly)) => {
                v.unscored += 1;
                continue;
            }
            Status::Unscored => {
                let text = format!(
                    "{key}: nothing is scored in {} records, none non-finite",
                    r.n
                );
                v.failures.push(Failure::new(Bar::Unscored, text));
                continue;
            }
            Status::NonFinite => {
                let text = format!("{key}: {} of {} records are non-finite", r.nonfinite, r.n);
                v.failures.push(Failure::new(Bar::NonFinite, text));
                continue;
            }
        };
        v.scored += 1;
        answered.insert(key.clone());
        let mut other = Vec::new();
        let best = index.best(&key, r.n, &mut other);
        for (name, n) in other {
            let text = format!("{key}: `{name}` scored {n} records, the candidate {}", r.n);
            v.failures.push(Failure::new(Bar::Corpus, text));
        }
        match best {
            None => v.unpaired += 1,
            Some((o, m)) => {
                v.paired += 1;
                if max > m {
                    // Which oracle won decides what the failure can mean, so it is said here and
                    // not left to a reader with a script (`0036`, draft).
                    // An oracle whose backend is not declared is left unclassified: saying
                    // "only a host-`std` oracle beats this" of it would be a claim with nothing
                    // behind it.
                    let why = if index.beaten_by_libm_crate(&key, r.n, max) {
                        v.dominated_same_backend += 1;
                        " [a `libm`-crate oracle also beats this: not D16's cost]"
                    } else if crate::conformance::oracle_backend(&o.name)
                        == Some(crate::conformance::Backend::HostStd)
                    {
                        v.dominated_host_std += 1;
                        " [only a host-`std` oracle beats this: D16 is a candidate, `0032` draft]"
                    } else {
                        ""
                    };
                    let text = format!(
                        "{key}: {max:e} u, over `{}` ({}) at {m:e} u{why}",
                        o.name, o.version
                    );
                    v.failures.push(Failure::new(Bar::Domination, text));
                }
            }
        }
        let Some(baseline) = baseline else { continue };
        match baseline.get(&key) {
            None => {
                let text = format!("{key}: no baseline row; `just envelope --bless` writes it");
                v.failures.push(Failure::new(Bar::NoRegress, text));
            }
            Some(&(n, _)) if n != r.n => {
                let text = format!("{key}: {} records, the baseline's maximum is over {n}", r.n);
                let bar = if r.n < n { Bar::Shrunk } else { Bar::NoRegress };
                v.failures.push(Failure::new(bar, text));
            }
            Some(&(_, base)) if max > base => {
                let text = format!("{key}: {max:e} u, over the baseline's {base:e} u");
                v.failures.push(Failure::new(Bar::NoRegress, text));
            }
            Some(&(_, base)) => v.improved += usize::from(max < base),
        }
    }
    for key in baseline.into_iter().flat_map(Baseline::keys) {
        if !answered.contains(key) {
            let text = format!("{key}: the baseline has it and no scored candidate row answers");
            v.failures.push(Failure::new(Bar::Shrunk, text));
        }
    }
    v
}

#[cfg(test)]
mod tests {
    use super::super::rows::testing::row;
    use super::*;

    const S: &str = "theta:1e-3";

    fn subject(name: &str, rows: Vec<Row>) -> Subject {
        Subject {
            name: name.into(),
            version: format!("{name}@1"),
            rows,
        }
    }

    fn one(name: &str, max: f64) -> Subject {
        subject(name, vec![row(name, "f", S, max)])
    }

    fn baseline(max: f64) -> Baseline {
        Baseline::from([(Key::of(&row("c", "f", S, 0.0)), (64, max))])
    }

    fn bars(v: &Verdict) -> Vec<&str> {
        v.failures.iter().map(|f| f.bar.name()).collect()
    }

    /// A failure is classified by whether **any** `libm`-crate oracle beats the candidate, not by
    /// whether the best one does: a same-backend oracle winning at all rules out D16's cost, and the
    /// two questions disagree on the real corpus by 30 against 26 (`0036`, draft).
    #[test]
    fn a_domination_failure_says_whether_a_same_backend_oracle_also_won() {
        let (c, base) = (one("c", 3.0), baseline(9.0));
        // sophus-rs is best, but `tf_tree_math` is on the `libm` crate and also beats us.
        let v = judge(
            &c,
            &[one("tf_tree_math", 2.0), one("sophus_rs", 1.0)],
            Some(&base),
        );
        assert_eq!(bars(&v), ["domination"]);
        assert_eq!(
            (v.dominated_same_backend, v.dominated_host_std),
            (1, 0),
            "{:?}",
            v.failures.iter().map(|f| &f.text).collect::<Vec<_>>()
        );
        assert!(
            v.failures[0].text.contains("not D16's cost"),
            "{:?}",
            v.failures[0].text
        );

        // The same best oracle, but now the `libm`-crate one does not beat us.
        let v = judge(
            &c,
            &[one("tf_tree_math", 4.0), one("sophus_rs", 1.0)],
            Some(&base),
        );
        assert_eq!((v.dominated_same_backend, v.dominated_host_std), (0, 1));
        assert!(
            v.failures[0].text.contains("D16 is a candidate"),
            "{:?}",
            v.failures[0].text
        );

        // An oracle with no declared backend is classified as neither, and says nothing.
        let v = judge(&c, &[one("someone_else", 1.0)], Some(&base));
        assert_eq!((v.dominated_same_backend, v.dominated_host_std), (0, 0));
        assert!(
            !v.failures[0].text.contains('['),
            "{:?}",
            v.failures[0].text
        );

        // No failure, no classification.
        let v = judge(&one("c", 0.5), &[one("tf_tree_math", 2.0)], Some(&base));
        assert!(v.failures.is_empty());
        assert_eq!((v.dominated_same_backend, v.dominated_host_std), (0, 0));
    }

    #[test]
    fn domination_passes_under_and_at_the_best_oracle_and_fails_one_ulp_over() {
        let dom = |max: f64| judge(&one("c", max), &[one("o", 2.0)], Some(&baseline(9.0)));
        for max in [1.0, 2.0] {
            let v = dom(max);
            assert!(v.failures.is_empty(), "{max}: {:?}", bars(&v));
            assert_eq!((v.scored, v.paired, v.unpaired), (1, 1, 0));
        }
        let v = dom(2.0f64.next_up());
        assert_eq!(bars(&v), ["domination"]);
        assert!(v.failures[0]
            .text
            .starts_with("f/theta:1e-3/f64: 2.0000000000000004e0 u, over `o` (o@1) at 2e0 u"));
    }

    #[test]
    fn the_bar_is_the_smallest_oracle_and_a_bad_oracle_is_no_bar() {
        let b = baseline(9.0);
        let weak_and_strong = [one("weak", 5.0), one("strong", 1.0)];
        let v = judge(&one("c", 3.0), &weak_and_strong, Some(&b));
        assert!(v.failures.len() == 1 && v.failures[0].text.contains("`strong`"));
        // A non-finite or unscored oracle row has no maximum, whatever its finite records scored:
        // every oracle bad, nothing to beat.
        let unscored = Row {
            max_u: f64::NAN,
            ..row("p", "f", S, 0.0)
        };
        for max in [f64::NAN, 0.5] {
            let bad = Row {
                nonfinite: 1,
                max_u: max,
                ..row("o", "f", S, 0.0)
            };
            let bad = [
                subject("o", vec![bad]),
                subject("p", vec![unscored.clone()]),
            ];
            let v = judge(&one("c", 1e30), &bad, Some(&baseline(1e30)));
            assert!(v.failures.is_empty() && (v.paired, v.unpaired) == (0, 1));
            // Beside a good one it leaves the good one's bar.
            let [bad, _] = bad;
            let v = judge(&one("c", 3.0), &[bad, one("g", 1.0)], Some(&b));
            assert_eq!(bars(&v), ["domination"]);
            assert!(v.failures[0].text.contains("`g`"));
        }
    }

    #[test]
    fn equal_best_oracles_name_the_first_in_the_order_given() {
        let b = baseline(9.0);
        let names = |order: [&str; 2]| {
            let oracles = order.map(|n| one(n, 1.0));
            let v = judge(&one("c", 3.0), &oracles, Some(&b));
            v.failures.into_iter().map(|f| f.text).collect::<Vec<_>>()
        };
        assert!(names(["a", "b"])[0].contains("over `a` (a@1)"));
        assert!(names(["b", "a"])[0].contains("over `b` (b@1)"));
    }

    #[test]
    fn a_stratum_no_oracle_answers_is_unpaired_and_only_the_baseline_judges_it() {
        let v = judge(&one("c", 7.0), &[], Some(&baseline(7.0)));
        assert!(v.failures.is_empty() && (v.paired, v.unpaired) == (0, 1));
        let other = subject(
            "o",
            vec![row("o", "g", S, 0.0), row("o", "f", "other", 0.0)],
        );
        let v = judge(&one("c", 8.0), &[other], Some(&baseline(7.0)));
        assert_eq!((bars(&v), v.unpaired), (vec!["no-regress"], 1));
    }

    #[test]
    fn no_regress_is_exact_equality_passes_and_one_ulp_worse_fails() {
        let no_regress = |max: f64| judge(&one("c", max), &[], Some(&baseline(3.0)));
        assert!(no_regress(3.0).failures.is_empty());
        let better = no_regress(3.0f64.next_down());
        assert!(better.failures.is_empty() && better.improved == 1);
        let worse = no_regress(3.0f64.next_up());
        assert_eq!(bars(&worse), ["no-regress"]);
        assert!(worse.failures[0]
            .text
            .contains("3.0000000000000004e0 u, over the baseline's 3e0 u"));
    }

    #[test]
    fn a_missing_baseline_file_row_stratum_or_size_fails() {
        let v = judge(&one("c", 1.0), &[], None);
        assert_eq!(bars(&v), ["no-regress"]);
        assert!(v.failures[0].text.contains("no baseline for `c`"));
        // A row the baseline lacks, and a baseline row nothing answers.
        let two = subject("c", vec![row("c", "f", S, 1.0), row("c", "f", "new", 1.0)]);
        let v = judge(&two, &[], Some(&baseline(1.0)));
        assert!(
            bars(&v) == ["no-regress"]
                && v.failures[0].text.starts_with("f/new/f64: no baseline row")
        );
        let v = judge(
            &one("c", 1.0),
            &[],
            Some(&Baseline::from([(
                Key::of(&row("c", "f", "gone", 0.0)),
                (64, 1.0),
            )])),
        );
        assert_eq!(bars(&v), ["no-regress", "no-regress"]);
        assert!(v.failures[1]
            .text
            .starts_with("f/gone/f64: the baseline has it"));
        // A dropped stratum is never the baseline's to move (D7); a new one is.
        assert_eq!(
            v.failures.iter().map(|f| f.bar).collect::<Vec<_>>(),
            [Bar::NoRegress, Bar::Shrunk]
        );
        // The same maximum over another number of records is another stratum: over more records
        // than the baseline's it moves the baseline, over fewer it narrows it.
        let over = |n: usize| {
            let mut b = baseline(1.0);
            b.values_mut().for_each(|e| e.0 = n);
            judge(&one("c", 1.0), &[], Some(&b))
        };
        let (wider, narrower) = (over(63), over(65));
        assert_eq!(bars(&wider), ["no-regress"]);
        assert!(wider.failures[0]
            .text
            .contains("64 records, the baseline's maximum is over 63"));
        assert_eq!(wider.failures[0].bar, Bar::NoRegress);
        assert_eq!(narrower.failures[0].bar, Bar::Shrunk);
        // An oracle scored on another corpus is no comparison.
        let other = Row {
            n: 65,
            ..row("o", "f", S, 1.0)
        };
        let v = judge(
            &one("c", 1.0),
            &[subject("o", vec![other])],
            Some(&baseline(1.0)),
        );
        assert_eq!((bars(&v), v.unpaired), (vec!["corpus"], 1));
    }

    #[test]
    fn a_non_finite_candidate_fails_and_only_a_backward_error_id_may_be_unscored() {
        // The harness leaves the finite records' maximum in a row with some non-finite ones.
        for max in [f64::NAN, 3.0] {
            let bad = Row {
                nonfinite: 3,
                max_u: max,
                ..row("c", "f", S, 0.0)
            };
            let v = judge(&subject("c", vec![bad]), &[], Some(&Baseline::new()));
            assert_eq!(bars(&v), ["non-finite"]);
            let text = &v.failures[0].text;
            assert!(text.contains("3 of 64 records are non-finite") && v.scored == 0);
        }
        let nothing = |fn_id: &str| Row {
            max_u: f64::NAN,
            ..row("c", fn_id, S, 0.0)
        };
        let judged = |fn_id: &str| {
            judge(
                &subject("c", vec![nothing(fn_id)]),
                &[],
                Some(&Baseline::new()),
            )
        };
        let v = judged("so3_from_matrix");
        assert!(v.failures.is_empty() && v.unscored == 1);
        // Any other function scored nothing only through a defect.
        let v = judged("so3_exp");
        assert_eq!((bars(&v), v.unscored), (vec!["unscored"], 0));
    }
}
