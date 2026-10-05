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
    /// A row of the exception table that no longer excepts anything: the stratum is dominated
    /// after all, or it is not scored, or the candidate has no row for it
    /// ([`0046`](../../../../docs/decisions/0046-explained-by-record-needs-a-record-to-point-at.md)
    /// item 4). An exception is a debt with a test attached and cannot outlive the defect.
    Exception,
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
            Bar::Exception => "exception",
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
    /// Domination failures where only a host-`std` oracle won and the one-variable twin closes
    /// the gap: D16's cost, measured and not inferred (`0037`, draft).
    pub(super) dominated_libm_bound: usize,
    /// Domination failures where only a host-`std` oracle won and the twin does **not** close the
    /// gap, so D16 does not account for it either. 37 of 45 on the real corpus (`0037`, draft).
    pub(super) dominated_unexplained: usize,
    /// Domination failures where only a host-`std` oracle won and the twin has no row to attribute
    /// them with.
    pub(super) dominated_host_std: usize,
    /// Strata an oracle won that the exception table explains (`0046`). They are **not** failures
    /// and **not** `paired` wins: the maximum still goes to the baseline and still has to not
    /// regress, so an excepted stratum is watched more closely than a dominated one, not less.
    pub(super) excepted: Vec<Excepted>,
}

/// One stratum the exception table explains, as the evidence page prints it.
pub(super) struct Excepted {
    pub(super) key: String,
    /// The candidate's maximum and the best oracle's, in `u`, and that oracle's name.
    pub(super) candidate: f64,
    pub(super) oracle: f64,
    pub(super) oracle_name: String,
    pub(super) record: String,
    pub(super) reason: String,
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

/// Whether the one-variable twin closes a domination failure: `Some(true)` when the candidate's
/// own program on the host's transcendentals reaches the best oracle's maximum, so the gap is
/// D16's; `Some(false)` when it does not, so D16 does not account for it; `None` when the twin
/// has no comparable row (`0037`, draft).
fn swap_closes(twin: Option<&Subject>, key: &Key, n: usize, best: f64) -> Option<bool> {
    let row = twin?.rows.iter().find(|r| Key::of(r) == *key && r.n == n)?;
    match status(row) {
        Status::Scored(m) => Some(m <= best),
        _ => None,
    }
}

/// `candidate` against `oracles` (in the order given: the first of equal bests is named) and, for
/// no-regress, `baseline`; `None` is a baseline file that does not exist, one failure and no row
/// read against it. With no twin, a failure only a host-`std` oracle wins is left unattributed.
pub(super) fn judge(
    candidate: &Subject,
    oracles: &[Subject],
    baseline: Option<&Baseline>,
) -> Verdict {
    judge_with(candidate, oracles, baseline, None)
}

/// [`judge`] with `twin`, the candidate's own program with one variable changed, which attributes
/// a failure no oracle on our backend wins. The twin is never judged and never a bar.
pub(super) fn judge_with(
    candidate: &Subject,
    oracles: &[Subject],
    baseline: Option<&Baseline>,
    twin: Option<&Subject>,
) -> Verdict {
    judge_excepting(
        candidate,
        oracles,
        baseline,
        twin,
        &super::exceptions::Exceptions::default(),
    )
}

/// [`judge_with`] with the exception table (`0046`): a listed `(fn, stratum, precision)` the
/// oracle wins is recorded rather than failed, and a listed one it does **not** win fails as a
/// stale row. It excepts domination and nothing else -- no-regress, non-finite, unscored, corpus
/// and coverage are untouched by it, which is item 2 and the reason the table is read here and
/// not at the call site.
pub(super) fn judge_excepting(
    candidate: &Subject,
    oracles: &[Subject],
    baseline: Option<&Baseline>,
    twin: Option<&Subject>,
    exceptions: &super::exceptions::Exceptions,
) -> Verdict {
    let index = Index::new(oracles);
    let mut used: BTreeSet<String> = BTreeSet::new();
    // Owned, so the `None` arm can name the twin the run actually read -- the candidate's own
    // where one exists, a stand-in's otherwise -- instead of a fixed name that may blame a subject
    // the run never read.
    let no_row = match twin {
        Some(t) => format!(
            " [only a host-`std` oracle beats this; `{}` has no row to attribute it]",
            t.name
        ),
        None => {
            " [only a host-`std` oracle beats this; no twin was read to attribute it]".to_string()
        }
    };
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
                    if let Some(e) = exceptions.get(&key.to_string()) {
                        used.insert(e.key());
                        v.excepted.push(Excepted {
                            key: e.key(),
                            candidate: max,
                            oracle: m,
                            oracle_name: o.name.clone(),
                            record: e.record.clone(),
                            reason: e.reason.clone(),
                        });
                        let Some(baseline) = baseline else { continue };
                        no_regress(&mut v, &key, r, max, baseline);
                        continue;
                    }
                    // Which oracle won decides what the failure can mean, so it is said here and
                    // not left to a reader with a script (`0036`, draft).
                    // An oracle whose backend is not declared is left unclassified: saying
                    // "only a host-`std` oracle beats this" of it would be a claim with nothing
                    // behind it.
                    let host_std = crate::conformance::oracle_backend(&o.name)
                        == Some(crate::conformance::Backend::HostStd);
                    let why = if index.beaten_by_libm_crate(&key, r.n, max) {
                        v.dominated_same_backend += 1;
                        " [a `libm`-crate oracle also beats this: not D16's cost]"
                    } else if !host_std {
                        ""
                    } else {
                        // Only a glibc-backed oracle won, so D16 is a candidate. Whether it is the
                        // cause is the twin's to say, and saying it without the twin was a claim
                        // with nothing behind it on 37 of 45 strata (`0037`, draft).
                        match swap_closes(twin, &key, r.n, m) {
                            Some(true) => {
                                v.dominated_libm_bound += 1;
                                " [the swap to the host's transcendentals closes this: D16's cost, \
                                 `0032` draft]"
                            }
                            Some(false) => {
                                v.dominated_unexplained += 1;
                                " [the swap to the host's transcendentals does not close this: not \
                                 D16's cost either, `0037` draft]"
                            }
                            None => {
                                v.dominated_host_std += 1;
                                // The twin that has no row is named by the caller: it is the
                                // candidate's own where one exists and a stand-in's otherwise, so
                                // a fixed name here would blame a subject the run never read.
                                &no_row
                            }
                        }
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
        no_regress(&mut v, &key, r, max, baseline);
    }
    for key in baseline.into_iter().flat_map(Baseline::keys) {
        if !answered.contains(key) {
            let text = format!("{key}: the baseline has it and no scored candidate row answers");
            v.failures.push(Failure::new(Bar::Shrunk, text));
        }
    }
    // Every row of the table has to be earning its keep: one that excepted nothing this run is
    // stale, and the message names it for deletion (`0046` item 4). This catches the three ways a
    // row stops being needed -- the stratum is dominated after all, it is no longer scored, and the
    // candidate has no row for it at all -- without having to tell them apart, which is the point:
    // an exception that explains nothing is to be removed whatever the reason.
    for e in exceptions.iter() {
        if !used.contains(&e.key()) {
            let text = format!(
                "`{}` excepts `{}` and no oracle wins it: delete the row (it cites `{}`)",
                super::exceptions::PATH,
                e.key(),
                e.record
            );
            v.failures.push(Failure::new(Bar::Exception, text));
        }
    }
    v
}

/// The no-regress bar for one scored candidate row, which an excepted stratum pays like any other
/// (`0046` item 2: the exception touches domination and nothing else).
fn no_regress(v: &mut Verdict, key: &Key, r: &Row, max: f64, baseline: &Baseline) {
    match baseline.get(key) {
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

    /// The exception table (`0046`), on the five cases item 2's plan names.
    ///
    /// A row excepts **domination and nothing else**: the fourth case is the one that matters, an
    /// excepted stratum that also regresses, which must still fail. If an exception ever silenced
    /// no-regress, the table would become the waiver list every waiver list becomes.
    #[test]
    fn an_exception_moves_domination_only_and_a_stale_row_fails() -> Result<(), String> {
        let table = |fn_id: &str, stratum: &str| {
            super::super::exceptions::Exceptions::parse(&format!(
                "[[exception]]\nfn = \"{fn_id}\"\nstratum = \"{stratum}\"\n\
                 precision = \"f64\"\nrecord = \"0046\"\nreason = \"a tie\"\n"
            ))
        };
        let oracles = [one("tf_tree_math", 1.0)];
        let judged = |cand: f64, base: f64, t: &super::super::exceptions::Exceptions| {
            let c = subject("c", vec![row("c", "f", S, cand)]);
            judge_excepting(&c, &oracles, Some(&baseline(base)), None, t)
        };

        // 1. An excepted failure: no domination failure, one excepted row, and the reason reaches
        //    the evidence page with both maxima.
        let t = table("f", S)?;
        let v = judged(2.0, 2.0, &t);
        assert_eq!(bars(&v), Vec::<&str>::new());
        assert_eq!(v.excepted.len(), 1);
        let e = v.excepted.first().ok_or("no excepted row")?;
        assert_eq!(
            (e.key.as_str(), e.record.as_str()),
            ("f/theta:1e-3/f64", "0046")
        );
        assert_eq!(
            (e.candidate, e.oracle, e.oracle_name.as_str()),
            (2.0, 1.0, "tf_tree_math")
        );
        // The excepted stratum is not counted as a win, and it is still paired.
        assert_eq!((v.paired, v.unpaired), (1, 0));

        // 2. A stale row: the stratum is dominated after all, so the row fails and is named.
        let v = judged(0.5, 0.5, &t);
        assert_eq!(bars(&v), ["exception"]);
        let text = &v.failures.first().ok_or("no failure")?.text;
        assert!(
            text.contains("f/theta:1e-3/f64") && text.contains("delete the row"),
            "{text}"
        );

        // 3. A row for a stratum the candidate does not answer is stale too: an exception that
        //    explains nothing goes, whatever the reason it stopped explaining.
        let v = judged(2.0, 2.0, &table("other_fn", S)?);
        assert_eq!(bars(&v), ["domination", "exception"]);

        // 4. **An exception does not touch no-regress.** The stratum is excepted *and* over its
        //    baseline: domination is silent, no-regress is not.
        let v = judged(2.0, 1.5, &t);
        assert_eq!(bars(&v), ["no-regress"]);
        assert_eq!(v.excepted.len(), 1);

        // 5. Nor the other bars: a non-finite candidate row is excepted by nothing. The
        //    `no-regress` here is the pre-existing `Shrunk` bar — the row is not scored, so the
        //    baseline's entry for it goes unanswered — and the exception is stale because a
        //    stratum with no maximum cannot be dominated.
        let mut bad = row("c", "f", S, f64::NAN);
        bad.nonfinite = 1;
        let c = subject("c", vec![bad]);
        let v = judge_excepting(&c, &oracles, Some(&baseline(2.0)), None, &t);
        assert_eq!(bars(&v), ["non-finite", "no-regress", "exception"]);
        Ok(())
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

        // The same best oracle, but now the `libm`-crate one does not beat us. With no twin the
        // failure is unattributed: D16 is a candidate and nothing here says it is the cause.
        let only_glibc = [one("tf_tree_math", 4.0), one("sophus_rs", 1.0)];
        let v = judge(&c, &only_glibc, Some(&base));
        assert_eq!((v.dominated_same_backend, v.dominated_host_std), (0, 1));
        assert!(
            v.failures[0].text.ends_with(
                "[only a host-`std` oracle beats this; no twin was read to attribute it]"
            ),
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

    /// What the twin adds: a failure only a glibc-backed oracle wins is D16's cost when the
    /// candidate's own program on the host's transcendentals reaches that oracle, and is **not**
    /// D16's when it does not. Saying "D16 is a candidate" without asking was wrong on 37 of 45
    /// strata of the real corpus (`0037`, draft).
    #[test]
    fn the_twin_says_whether_the_swap_closes_a_failure_only_glibc_won() {
        let (c, base) = (one("c", 3.0), baseline(9.0));
        let oracles = [one("tf_tree_math", 4.0), one("sophus_rs", 1.0)];
        let twin = |max: f64| one(crate::seeded::TWIN, max);

        // The swap reaches the best oracle: D16's cost, measured.
        let v = judge_with(&c, &oracles, Some(&base), Some(&twin(1.0)));
        assert_eq!(
            (
                v.dominated_libm_bound,
                v.dominated_unexplained,
                v.dominated_host_std
            ),
            (1, 0, 0)
        );
        assert!(v.failures[0].text.contains("closes this: D16's cost"));

        // The swap changes nothing: not D16's either, which is the 37.
        let v = judge_with(&c, &oracles, Some(&base), Some(&twin(3.0)));
        assert_eq!(
            (
                v.dominated_libm_bound,
                v.dominated_unexplained,
                v.dominated_host_std
            ),
            (0, 1, 0)
        );
        assert!(v.failures[0].text.contains("not D16's cost either"));

        // A twin whose row is for another stratum attributes nothing, as no twin does.
        let mut elsewhere = twin(1.0);
        elsewhere.rows[0].stratum = "another".to_string();
        let v = judge_with(&c, &oracles, Some(&base), Some(&elsewhere));
        assert_eq!(
            (
                v.dominated_libm_bound,
                v.dominated_unexplained,
                v.dominated_host_std
            ),
            (0, 0, 1)
        );

        // A same-backend oracle winning answers the question first: the twin is not consulted.
        let also_tf = [one("tf_tree_math", 2.0), one("sophus_rs", 1.0)];
        let v = judge_with(&c, &also_tf, Some(&base), Some(&twin(1.0)));
        assert_eq!((v.dominated_same_backend, v.dominated_libm_bound), (1, 0));

        // The twin is never a bar: it is not an oracle and cannot make a failure.
        let v = judge_with(
            &one("c", 0.5),
            &[one("sophus_rs", 1.0)],
            Some(&base),
            Some(&twin(0.1)),
        );
        assert!(v.failures.is_empty());
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
