//! `cargo xtask envelope [--candidate NAME] [--oracle NAME]... [--check | --bless [--dry-run]]`
//! (`docs/PHASE1.md` §8, `0006`): merges the result CSVs `conformance` writes under
//! `conformance/results/` and judges the candidate's rows by the bars of [`bars`], then
//! [`coverage`]. `--bless` writes `conformance/baseline/<candidate>.csv`, the candidate's scored
//! result rows with `git_rev` blank (the no-regress bar, and a diff of the numbers alone), and
//! `docs/evidence/ENVELOPE.md` ([`evidence`]); `--dry-run` says what it would
//! write and writes nothing; `--check` fails when either file is not what `--bless` would write,
//! as `thresholds --check` does, so an improved maximum is blessed in its own PR.
//!
//! **Binary64 only.** `f32` scoring and the `@f32` strata (`0016`) landed after this module
//! and it reads the binary64 result files alone, although `PHASE1.md` §8 states the bars per
//! precision. The `@f32` half of the envelope is owed.
//!
//! Where the specs are silent, the smallest reading was taken:
//!
//! - **The candidate** is [`CANDIDATE`], the in-process subject `helicoid` (Phase 3), read from
//!   its own CSV, and the oracles are the runners of `conformance::oracle_names`, each of which
//!   must have written its CSV. While the candidate is not an in-process subject that
//!   `conformance` runs, having no rows is the state before Phase 3: nothing is judged, the run
//!   says so and checks coverage only, and a baseline with no candidate rows fails. Once it is
//!   registered, no rows fail the run, so the gate switches on with the subject and not with a
//!   code change. `--candidate` names another registered subject (`seeded:correct` is the
//!   stand-in) for a read-only run; any other name is refused.
//! - **`--bless`** takes the candidate's scored rows as they are and reports, rather than fails
//!   on, what only moves the baseline (a worse max, a stratum the baseline lacks, more records);
//!   the diff is the review (§8). It writes nothing while any other failure stands: domination
//!   (`PHASE3.md` §10 blesses a dominated envelope), non-finite or unscored rows, coverage, an
//!   oracle over another corpus, and a baseline stratum no row answers or one over more records
//!   than the candidate's, which D7 forbids.
//! - **`--check`** compares the page without the columns the oracle runners fill in
//!   ([`evidence::candidate_view`]), and the baseline exactly. A runner's numbers depend on its own
//!   resolved dependencies, and a runner on the host's `std` on the host's system library too, which
//!   D16 does not cover (`docs/maths/error-analysis.md` EA.13(d)). Not every runner is host-bound:
//!   `tf_tree_math` routes through the `libm` crate as we do, which is what makes it the control in
//!   `0032` (draft) and what the domination split reports (`conformance::Backend`).
//! - **The CSVs read** are `results_path`'s for the named subjects, so a stale `--fn` run's file
//!   or a planted defect's is never merged.

mod bars;
mod coverage;
mod evidence;
mod rows;

use std::fmt::Write;
use std::path::Path;

use crate::conformance::report::{csv, Row};
use crate::conformance::{corpus, oracle_names, results_path, root, subject};
use bars::{Bar, Failure};
use helicoid_linalg::Precision;
use rows::Subject;

const USAGE: &str = "usage: cargo xtask envelope [--candidate NAME] [--oracle NAME]... \
[--check | --bless [--dry-run]]";

/// The subject the bars judge (`docs/PHASE3.md` §10). It joins `subject::registry` in Phase 3;
/// nothing writes rows under this name before.
const CANDIDATE: &str = "helicoid";

const CORPUS: &str = "conformance/corpus";
const RESULTS: &str = "conformance/results";
const BASELINE: &str = "conformance/baseline";

#[derive(Clone)]
struct Options {
    candidate: String,
    /// The oracles by name; the runners when not given.
    oracles: Option<Vec<String>>,
    check: bool,
    bless: bool,
    dry_run: bool,
    /// The in-process subjects `conformance` runs: a candidate among them owes rows.
    registered: Vec<String>,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            candidate: CANDIDATE.to_string(),
            oracles: None,
            check: false,
            bless: false,
            dry_run: false,
            registered: Vec::new(),
        }
    }
}

fn parse_args(args: &[String]) -> Result<Options, String> {
    let mut o = Options {
        registered: subject::registry()
            .iter()
            .map(|s| s.subject.name().to_string())
            .collect(),
        ..Options::default()
    };
    let mut rest = args.iter();
    while let Some(arg) = rest.next() {
        let (flag, inline) = match arg.split_once('=') {
            Some((f, v)) => (f, Some(v.to_string())),
            None => (arg.as_str(), None),
        };
        match (flag, &inline) {
            ("--check", None) => o.check = true,
            ("--bless", None) => o.bless = true,
            ("--dry-run", None) => o.dry_run = true,
            ("--candidate" | "--oracle", _) => {
                let value = inline.or_else(|| rest.next().cloned());
                let value = value.ok_or_else(|| format!("`{flag}` needs a value; {USAGE}"))?;
                match flag {
                    "--candidate" => o.candidate = value,
                    _ => o.oracles.get_or_insert_with(Vec::new).push(value),
                }
            }
            _ => return Err(format!("unknown argument `{arg}`; {USAGE}")),
        }
    }
    let known = o.candidate == CANDIDATE || o.registered.contains(&o.candidate);
    match (
        o.check && o.bless,
        o.dry_run && !o.bless,
        o.bless && o.candidate != CANDIDATE,
        known,
    ) {
        (true, ..) => Err(format!(
            "`--check` and `--bless` exclude each other; {USAGE}"
        )),
        (_, true, ..) => Err(format!("`--dry-run` is for `--bless`; {USAGE}")),
        (_, _, true, _) => Err(format!(
            "`--bless` writes the baseline of `{CANDIDATE}` only, not a stand-in's; {USAGE}"
        )),
        (.., false) => Err(format!(
            "`{}` is neither `{CANDIDATE}` nor an in-process subject ({}); {USAGE}",
            o.candidate,
            o.registered.join(", ")
        )),
        _ => Ok(o),
    }
}

/// What a run says, and what failed.
struct Report {
    text: String,
    failures: Vec<Failure>,
}

/// `parse` of the file at `root/rel`; `None` when there is none.
fn read_optional<T>(
    root: &Path,
    rel: &Path,
    parse: impl FnOnce(&str) -> Result<T, String>,
) -> Result<Option<T>, String> {
    let path = root.join(rel);
    match std::fs::read_to_string(&path) {
        Ok(text) => parse(&text)
            .map(Some)
            .map_err(|e| format!("{}: {e}", path.display())),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(e) => Err(format!("{}: {e}", path.display())),
    }
}

fn read_subject(root: &Path, name: &str) -> Result<Option<Subject>, String> {
    let rel = results_path(Path::new(RESULTS), name, None, Precision::F64);
    read_optional(root, &rel, |text| rows::parse_results(name, text))
}

fn write_file(root: &Path, rel: &str, text: &str) -> Result<(), String> {
    let path = root.join(rel);
    let dir = path.parent().ok_or("no parent dir")?;
    std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    std::fs::write(&path, text).map_err(|e| format!("{}: {e}", path.display()))
}

/// The oracles' rows: the named ones, else the runners'; each must have written its CSV.
fn read_oracles(root: &Path, options: &Options) -> Result<Vec<Subject>, String> {
    let names = options.oracles.clone();
    let names = names.unwrap_or_else(|| oracle_names().into_iter().map(String::from).collect());
    let mut oracles: Vec<Subject> = Vec::new();
    for name in names {
        if name == options.candidate || oracles.iter().any(|o| o.name == name) {
            return Err(format!("`{name}` is the candidate or named twice"));
        }
        let missing = || format!("no result rows for oracle `{name}`: run its oracle recipe");
        oracles.push(read_subject(root, &name)?.ok_or_else(missing)?);
    }
    Ok(oracles)
}

/// The run's opening lines: what was read and what each bar covered.
fn summary(
    c: &Subject,
    oracles: &[Subject],
    v: &bars::Verdict,
    has_baseline: bool,
    twin: Option<&str>,
) -> String {
    let name = &c.name;
    let versions: Vec<String> = oracles
        .iter()
        .map(|o| format!("{} ({})", o.name, o.version))
        .collect();
    let no_regress = match has_baseline {
        true => format!(
            "{} scored rows, {} under their baseline",
            v.scored, v.improved
        ),
        false => format!("no baseline: {} scored rows are not judged", v.scored),
    };
    // The split by backend is the first thing to know about a domination failure: an oracle on the
    // `libm` crate computes the same transcendentals we do, so a stratum it wins is the program's
    // (`0036`, draft).
    let classified = v.dominated_same_backend
        + v.dominated_libm_bound
        + v.dominated_unexplained
        + v.dominated_host_std;
    let split = match classified {
        0 => String::new(),
        _ => format!(
            "\n  domination failures: {} the program's (an oracle on the `libm` crate won), \
             {} D16's (`{}` closes the gap), {} neither (it does not), {} unattributed (no twin \
             row)",
            v.dominated_same_backend,
            v.dominated_libm_bound,
            twin.unwrap_or("no twin"),
            v.dominated_unexplained,
            v.dominated_host_std
        ),
    };
    format!(
        "envelope: candidate `{name}` ({}), {} rows; oracles: {}\n  \
         domination: {} strata paired with an oracle, {} with none{split}\n  \
         no-regress: {no_regress}\n  \
         not scored: {} rows\n",
        c.version,
        c.rows.len(),
        if versions.is_empty() {
            "none".to_string()
        } else {
            versions.join(", ")
        },
        v.paired,
        v.unpaired,
        v.unscored
    )
}

/// The strata, as their failure lines, on which `oracle` beats `candidate` by the domination bar of
/// [`bars`] alone: `--self-test` reads a planted defect against the correct kernel (`PHASE1.md`
/// §10) with it.
pub(crate) fn beaten_strata(candidate: &[Row], oracle: &[Row]) -> Vec<String> {
    let subject = |rows: &[Row]| Subject {
        name: rows.first().map_or_else(String::new, |r| r.subject.clone()),
        version: rows.first().map_or_else(String::new, |r| r.version.clone()),
        rows: rows.to_vec(),
    };
    let verdict = bars::judge(&subject(candidate), &[subject(oracle)], None);
    let beaten = verdict.failures.into_iter();
    beaten
        .filter(|f| f.bar == Bar::Domination)
        .map(|f| f.text)
        .collect()
}

/// The function ids of the corpus whose file exists and has records: `MANIFEST.json` lists a
/// file, it does not make one.
fn corpus_ids(root: &Path) -> Result<Vec<String>, String> {
    let dir = root.join(CORPUS);
    let entries = corpus::manifest(&dir)?;
    let present =
        |e: &corpus::Entry| e.records > 0 && dir.join(format!("{}.jsonl", e.fn_id)).is_file();
    Ok(entries
        .into_iter()
        .filter(present)
        .map(|e| e.fn_id)
        .collect())
}

/// One run over the checkout at `root`: everything it reads and writes is under it, and its
/// text names paths relative to it, so two runs on the same files say the same bytes.
fn execute(root: &Path, options: &Options) -> Result<Report, String> {
    let name = options.candidate.as_str();
    let cover = coverage::check(&corpus_ids(root)?, &|phase, row| {
        coverage::landed_in(root, phase, row)
    });
    let candidate = read_subject(root, name)?.filter(|s| !s.rows.is_empty());
    let baseline_rel = results_path(Path::new(BASELINE), name, None, Precision::F64);
    let baseline = read_optional(root, &baseline_rel, |text| {
        rows::baseline_of(&rows::parse_results(name, text)?)
    })?;
    let baseline_rel = baseline_rel.to_string_lossy().into_owned();
    // Oracle rows are read only when there is something to judge against them.
    let oracles = match candidate {
        Some(_) => read_oracles(root, options)?,
        None => Vec::new(),
    };
    // The twin is optional and never an oracle: it attributes a failure, it does not score one
    // (`0037`, draft). Absent, the failures it would explain are reported unattributed.
    // A candidate's own program with one variable changed is `<candidate>:host-std` where that
    // subject exists; `seeded:host-std` is the fallback, which is the stand-in's program and so
    // can only attribute the ids that program answers. Preferring the candidate's own is what
    // `bars::judge_with`'s contract says the twin is.
    let twin = match &candidate {
        Some(c) => {
            let own = format!("{}{}", c.name, crate::seeded::TWIN_SUFFIX);
            match read_subject(root, &own)?.filter(|s| !s.rows.is_empty()) {
                Some(s) => Some(s),
                None => read_subject(root, crate::seeded::TWIN)?.filter(|s| !s.rows.is_empty()),
            }
        }
        None => None,
    };

    let mut failures = cover.failures;
    let verdict = match &candidate {
        Some(c) => bars::judge_with(c, &oracles, baseline.as_ref(), twin.as_ref()),
        None => bars::Verdict::default(),
    };
    let registered = options.registered.iter().any(|r| r == name);
    let mut text = match &candidate {
        Some(c) => summary(
            c,
            &oracles,
            &verdict,
            baseline.is_some(),
            twin.as_ref().map(|t| t.name.as_str()),
        ),
        None if registered => format!(
            "envelope: no result rows for the candidate `{name}`: nothing is judged, and it is \
             owed\n"
        ),
        None => format!(
            "envelope: no result rows for the candidate `{name}`: domination and no-regress are \
             NOT evaluated (its in-process subject is Phase 3), only coverage is\n"
        ),
    };
    failures.extend(verdict.failures);
    if candidate.is_none() {
        if registered {
            let owed = format!(
                "`{name}` is an in-process subject and wrote no rows: run \
                 `cargo xtask conformance --subject {name}`"
            );
            failures.push(Failure::new(Bar::Candidate, owed));
        }
        if baseline.is_some() {
            let stale = format!("`{baseline_rel}` exists and the candidate wrote no rows");
            failures.push(Failure::new(Bar::Shrunk, stale));
        }
    }
    let _ = writeln!(
        text,
        "  coverage: {} required ids have corpus files, {} owed ids excused until their phase's row is Done",
        cover.required, cover.excused
    );

    // What `--bless` writes, and what `--check` holds the committed files to.
    let mut files = vec![(
        evidence::PATH.to_string(),
        evidence::page(name, candidate.as_ref(), &oracles),
    )];
    if let Some(c) = &candidate {
        let scored: Vec<Row> = c
            .rows
            .iter()
            .filter(|r| bars::is_scored(r))
            .cloned()
            .collect();
        files.push((baseline_rel, csv(&scored, "")));
    }
    if options.bless {
        // The new baseline replaces the old: what it only moves is shown, not failed.
        let (moved, hard): (Vec<_>, Vec<_>) =
            failures.into_iter().partition(|f| f.bar.moves_baseline());
        for f in &moved {
            let _ = writeln!(text, "  changes the baseline: {}", f.text);
        }
        failures = hard;
        if failures.is_empty() {
            let verb = if options.dry_run {
                "would write"
            } else {
                "wrote"
            };
            for (rel, contents) in &files {
                let _ = writeln!(text, "  {verb} {rel} ({} lines)", contents.lines().count());
                if !options.dry_run {
                    write_file(root, rel, contents)?;
                }
            }
        } else {
            let _ = writeln!(
                text,
                "  wrote nothing: the failures below are not the baseline's to move"
            );
        }
    }
    if options.check {
        for (rel, fresh) in &files {
            let view = |t: &str| match rel.as_str() {
                evidence::PATH => evidence::candidate_view(t),
                _ => t.to_string(),
            };
            let drift = match std::fs::read_to_string(root.join(rel)) {
                Ok(committed) => evidence::first_difference(&view(&committed), &view(fresh))
                    .map(|at| format!("`{rel}` is not what `--bless` writes ({at})")),
                Err(_) => Some(format!("`{rel}` does not exist")),
            };
            if let Some(d) = drift {
                let text = format!("{d}; run `just envelope --bless` and review the diff");
                failures.push(Failure::new(Bar::Drift, text));
            }
        }
    }
    for f in &failures {
        let _ = writeln!(text, "FAIL {}: {}", f.bar.name(), f.text);
    }
    Ok(Report { text, failures })
}

#[allow(clippy::print_stdout)]
pub(crate) fn run(args: &[String]) -> Result<(), String> {
    let options = parse_args(args)?;
    let report = execute(&root()?, &options)?;
    print!("{}", report.text);
    match report.failures.len() {
        0 => Ok(()),
        n => Err(format!("{n} failure(s)")),
    }
}

#[cfg(test)]
mod tests {
    use super::rows::testing::row;
    use super::*;
    use crate::conformance::report;
    use crate::conformance::testkit::Scratch;

    fn args(a: &[&str]) -> Vec<String> {
        a.iter().map(ToString::to_string).collect()
    }

    #[test]
    fn arguments_take_both_spellings_and_refuse_the_rest() -> Result<(), String> {
        let a = [
            "--candidate=seeded:correct",
            "--oracle",
            "a",
            "--oracle=b",
            "--check",
        ];
        let o = parse_args(&args(&a))?;
        assert_eq!((o.candidate.as_str(), o.check), ("seeded:correct", true));
        assert_eq!(o.oracles, Some(args(&["a", "b"])));
        assert!(o.registered.contains(&o.candidate));
        let d = parse_args(&[])?;
        assert_eq!(
            (d.candidate.as_str(), d.oracles.is_none()),
            (CANDIDATE, true)
        );
        assert!(parse_args(&args(&["--bless", "--dry-run"]))?.dry_run);
        for bad in [
            &["--check", "--bless"][..],
            &["--dry-run"],
            &["--dry-run", "--check"],
            &["--bless", "--candidate", "seeded:correct"],
            &["--candidate"],
            &["--candidate", "nosuch"],
            &["--oracle"],
            &["--bless=x"],
            &["--nope"],
            &["x"],
        ] {
            assert!(parse_args(&args(bad)).is_err(), "{bad:?}");
        }
        let e = parse_args(&args(&["--candidate", "nosuch"])).err();
        assert!(e.is_some_and(|e| e.contains("neither `helicoid` nor an in-process subject")));
        Ok(())
    }

    /// The manifest of a checkout that lists `ids`, each with one record.
    fn manifest(ids: &[String]) -> String {
        let files: Vec<String> = ids
            .iter()
            .map(|id| format!("\"{id}.jsonl\":{{\"kind\":\"corpus\",\"records\":1}}"))
            .collect();
        format!("{{\"files\":{{{}}}}}", files.join(","))
    }

    /// A checkout with the corpus files the specs require (a stub line each: nothing reads a
    /// record), the docs coverage reads with no phase landed, and no results.
    fn world(tag: &str) -> Result<Scratch, String> {
        let scratch = Scratch::new(tag);
        let dir = scratch.0.join(CORPUS);
        std::fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
        let ids = coverage::ids_for_tests();
        for id in &ids {
            std::fs::write(dir.join(format!("{id}.jsonl")), "{}\n").map_err(|e| e.to_string())?;
        }
        std::fs::write(dir.join("MANIFEST.json"), manifest(&ids)).map_err(|e| e.to_string())?;
        for (path, text) in coverage::docs_for_tests(&[]) {
            write_file(&scratch.0, &path, &text)?;
        }
        Ok(scratch)
    }

    fn put(scratch: &Scratch, subject: &str, rows: &[Row]) -> Result<(), String> {
        let rel = results_path(Path::new(RESULTS), subject, None, Precision::F64);
        write_file(&scratch.0, &rel.to_string_lossy(), &csv(rows, "rev"))
    }

    fn read(scratch: &Scratch, rel: &str) -> Result<String, String> {
        std::fs::read_to_string(scratch.0.join(rel)).map_err(|e| e.to_string())
    }

    fn options(candidate: &str, oracles: &[&str], mode: &str) -> Options {
        Options {
            candidate: candidate.into(),
            oracles: Some(oracles.iter().map(ToString::to_string).collect()),
            check: mode == "check",
            bless: mode == "bless" || mode == "dry",
            dry_run: mode == "dry",
            registered: Vec::new(),
        }
    }

    /// `candidate` against `oracles`; `mode` is `""`, `check`, `bless` or `dry` (a dry bless).
    fn go(
        scratch: &Scratch,
        candidate: &str,
        oracles: &[&str],
        mode: &str,
    ) -> Result<Report, String> {
        execute(&scratch.0, &options(candidate, oracles, mode))
    }

    fn bars(r: &Report) -> Vec<&str> {
        r.failures.iter().map(|f| f.bar.name()).collect()
    }

    const BASELINE_OF_CAND: &str = "conformance/baseline/cand.csv";

    fn strata(subject: &str, maxima: &[f64]) -> Vec<Row> {
        let all = ["theta:1e-3", "theta:1e-2", "theta:1e-1"];
        let rows = maxima.iter().zip(all);
        rows.map(|(&m, s)| row(subject, "so3_exp", s, m)).collect()
    }

    #[test]
    fn until_the_candidate_has_rows_only_coverage_runs_and_the_run_says_so() -> Result<(), String> {
        let scratch = world("no-candidate")?;
        let run = |mode: &str| go(&scratch, CANDIDATE, &[], mode);
        // The page does not exist yet: that is drift, and the only failure.
        let r = run("check")?;
        assert_eq!(bars(&r), ["drift"]);
        assert!(r.text.contains("NOT evaluated") && r.text.contains("46 required ids"));
        let r = run("bless")?;
        assert!(r.failures.is_empty() && r.text.contains("wrote docs/evidence/ENVELOPE.md"));
        assert!(run("check")?.failures.is_empty());
        // A baseline that nothing answers fails, and is not blessed away.
        let stale = csv(&strata(CANDIDATE, &[1.0]), "");
        write_file(&scratch.0, "conformance/baseline/helicoid.csv", &stale)?;
        assert_eq!(bars(&run("")?), ["no-regress"]);
        let r = run("bless")?;
        assert!(r.text.contains("wrote nothing"), "{}", r.text);
        assert_eq!(bars(&r), ["no-regress"]);
        assert_eq!(read(&scratch, "conformance/baseline/helicoid.csv")?, stale);
        Ok(())
    }

    #[test]
    fn a_registered_candidate_with_no_rows_fails_and_the_gate_needs_no_code_change(
    ) -> Result<(), String> {
        let scratch = world("registered")?;
        let registered = |mode: &str| Options {
            registered: args(&["cand"]),
            ..options("cand", &[], mode)
        };
        let owed = "FAIL candidate: `cand` is an in-process subject and wrote no rows: run \
                    `cargo xtask conformance --subject cand`";
        // No file, and a file with a header only, are the same.
        for header_only in [false, true] {
            if header_only {
                put(&scratch, "cand", &[])?;
            }
            let r = execute(&scratch.0, &registered(""))?;
            assert_eq!(bars(&r), ["candidate"]);
            assert!(r.text.contains(owed), "{}", r.text);
            assert!(!r.text.contains("NOT evaluated"), "{}", r.text);
            // A bless writes nothing, not even the page.
            let r = execute(&scratch.0, &registered("bless"))?;
            assert_eq!(bars(&r), ["candidate"]);
            assert!(!scratch.0.join(evidence::PATH).exists());
        }
        // Not registered (the state before Phase 3), the same files are the pass-with-message.
        assert!(go(&scratch, "cand", &[], "")?.failures.is_empty());
        // With rows it is judged, and nothing is owed.
        put(&scratch, "cand", &strata("cand", &[1.0]))?;
        let r = execute(&scratch.0, &registered(""))?;
        assert_eq!(bars(&r), ["no-regress"]);
        Ok(())
    }

    /// The bars that failed, and whether `so2_exp`'s missing file is named.
    fn so2_exp_missing(r: &Report) -> (Vec<&str>, bool) {
        let named = r
            .text
            .contains("FAIL coverage: `so2_exp` has no corpus file");
        (bars(r), named)
    }

    #[test]
    fn a_required_id_whose_file_is_missing_empty_or_unlisted_fails_the_run() -> Result<(), String> {
        let scratch = world("no-id")?;
        let dir = scratch.0.join(CORPUS);
        let file = dir.join("so2_exp.jsonl");
        let want = (vec!["coverage"], true);
        // Listed, and deleted: the manifest alone is no file.
        std::fs::remove_file(&file).map_err(|e| e.to_string())?;
        assert_eq!(so2_exp_missing(&go(&scratch, CANDIDATE, &[], "")?), want);
        // Listed, present, and with no records.
        std::fs::write(&file, "").map_err(|e| e.to_string())?;
        let manifest = read(&scratch, "conformance/corpus/MANIFEST.json")?;
        let empty = manifest.replace(
            "\"so2_exp.jsonl\":{\"kind\":\"corpus\",\"records\":1}",
            "\"so2_exp.jsonl\":{\"kind\":\"corpus\",\"records\":0}",
        );
        assert_ne!(empty, manifest);
        std::fs::write(dir.join("MANIFEST.json"), empty).map_err(|e| e.to_string())?;
        assert_eq!(so2_exp_missing(&go(&scratch, CANDIDATE, &[], "")?), want);
        // Not listed: the manifest reader refuses the unlisted file.
        let cut = manifest.replace("\"so2_exp.jsonl\":{\"kind\":\"corpus\",\"records\":1},", "");
        std::fs::write(dir.join("MANIFEST.json"), cut).map_err(|e| e.to_string())?;
        std::fs::remove_file(&file).map_err(|e| e.to_string())?;
        assert_eq!(so2_exp_missing(&go(&scratch, CANDIDATE, &[], "")?), want);
        Ok(())
    }

    #[test]
    fn an_owed_group_is_due_when_its_row_is_done_whatever_the_rest_of_the_table_says(
    ) -> Result<(), String> {
        let scratch = world("owed")?;
        for (path, text) in coverage::docs_for_tests(&["corpus ids (§6)"]) {
            write_file(&scratch.0, &path, &text)?;
        }
        let r = go(&scratch, CANDIDATE, &[], "")?;
        assert_eq!(bars(&r), ["coverage", "coverage", "coverage"]);
        assert!(r.text.contains(
            "FAIL coverage: phase 2's row `corpus ids (§6)` is Done and `eig3` (PHASE2.md §6) has no corpus file"
        ));
        // A checkout whose docs lack the row cannot tell, and does not excuse.
        std::fs::remove_file(scratch.0.join("docs/PHASE6.md")).map_err(|e| e.to_string())?;
        let r = go(&scratch, CANDIDATE, &[], "")?;
        assert_eq!(r.failures.len(), 4);
        assert!(r.text.contains("docs/PHASE6.md: "), "{}", r.text);
        Ok(())
    }

    /// The rows every candidate below has: three strata of `so3_exp`, and one of
    /// `so3_from_matrix` with nothing scored (backward error is owed).
    fn cand_rows(maxima: &[f64]) -> Vec<Row> {
        let mut rows = strata("cand", maxima);
        rows.push(Row {
            max_u: f64::NAN,
            p99_u: f64::NAN,
            ..row("cand", "so3_from_matrix", "q:w0", 0.0)
        });
        rows
    }

    #[test]
    fn a_candidate_is_judged_blessed_and_checked_and_the_report_is_a_function_of_the_files(
    ) -> Result<(), String> {
        let scratch = world("e2e")?;
        put(&scratch, "cand", &cand_rows(&[1.0, 2.0, 3.0]))?;
        put(&scratch, "o1", &strata("o1", &[1.0, 2.5, 9.0]))?;
        put(&scratch, "o2", &strata("o2", &[4.0, 2.0, 3.5]))?;
        let run = |mode: &str| go(&scratch, "cand", &["o1", "o2"], mode);
        // Ties (1.0 against o1, 2.0 against o2) and a win pass; there is no baseline yet.
        let r = run("")?;
        assert_eq!(bars(&r), ["no-regress"]);
        assert!(r
            .text
            .contains("3 strata paired with an oracle, 0 with none"));
        assert!(r.text.contains("no baseline for `cand`"));
        assert!(r.text.contains("no baseline: 3 scored rows are not judged"));
        // A dry bless writes nothing, neither file; a bless writes the same bytes twice.
        let dry = run("dry")?;
        assert!(dry
            .text
            .contains("would write conformance/baseline/cand.csv (4 lines)"));
        assert!(!scratch.0.join(BASELINE).exists());
        assert!(!scratch.0.join(evidence::PATH).exists());
        assert!(run("bless")?.failures.is_empty());
        let (baseline, page) = (
            read(&scratch, BASELINE_OF_CAND)?,
            read(&scratch, evidence::PATH)?,
        );
        // The unscored row is no baseline row (the next run would abort on it), and p99 is not max.
        let mut want = String::from(report::HEADER) + "\n";
        let p99 = ["5e-1", "1e0", "1.5e0"];
        for ((m, p), k) in ["1e0", "2e0", "3e0"]
            .iter()
            .zip(p99)
            .zip(["1e-3", "1e-2", "1e-1"])
        {
            want += &format!("so3_exp,theta:{k},f64,cand,cand@1,64,{m},{p},0,0,\n");
        }
        assert_eq!(baseline, want);
        assert!(page.contains("| `so3_exp` | `theta:1e-2` | f64 | 64 | 2.0000e0 | 1.0000e0 | `o2` | o2@1 | 2.0000e0 |"));
        assert!(page.contains("| `so3_from_matrix` | `q:w0` | f64 | 64 | — | — | — | — | — |"));
        assert!(run("bless")?.failures.is_empty());
        assert_eq!(read(&scratch, BASELINE_OF_CAND)?, baseline);
        assert_eq!(read(&scratch, evidence::PATH)?, page);
        let (a, b) = (run("check")?, run("check")?);
        assert!(a.failures.is_empty() && a.text == b.text, "{}", a.text);
        assert!(a.text.contains("3 scored rows, 0 under their baseline"));
        assert!(a.text.contains("not scored: 1 rows"));
        // The oracles' columns are the runners' on their own libm: a different best oracle and
        // fifth digit is not drift.
        put(&scratch, "o2", &strata("o2", &[4.0, 2.0625, 3.5]))?;
        put(&scratch, "o1", &strata("o1", &[1.0, 2.03125, 9.0]))?;
        assert!(run("check")?.failures.is_empty());
        assert_ne!(read(&scratch, evidence::PATH)?, {
            let fresh = run("bless")?;
            assert!(fresh.failures.is_empty());
            read(&scratch, evidence::PATH)?
        });
        Ok(())
    }

    #[test]
    fn a_planted_regression_and_a_planted_dominated_case_fail_and_an_improvement_is_drift(
    ) -> Result<(), String> {
        let scratch = world("planted")?;
        put(&scratch, "o1", &strata("o1", &[1.0, 2.5, 9.0]))?;
        put(&scratch, "cand", &strata("cand", &[1.0, 2.0, 3.0]))?;
        let run = |mode: &str| go(&scratch, "cand", &["o1"], mode);
        assert!(run("bless")?.failures.is_empty());
        // One ulp worse on one stratum fails no-regress there and nowhere else.
        let mut rows = strata("cand", &[1.0, 2.0, 3.0]);
        rows[1].max_u = 2.0f64.next_up();
        put(&scratch, "cand", &rows)?;
        let r = run("")?;
        assert_eq!(bars(&r), ["no-regress"]);
        let worse = "so3_exp/theta:1e-2/f64: 2.0000000000000004e0 u, over the baseline's 2e0 u";
        assert!(
            r.text.contains(&format!("FAIL no-regress: {worse}")),
            "{}",
            r.text
        );
        // Blessing the worse maximum passes, and says what it moved.
        let r = run("bless")?;
        assert!(
            r.failures.is_empty()
                && r.text
                    .contains("changes the baseline: so3_exp/theta:1e-2/f64")
        );
        // One ulp over the oracle fails domination, and a bless then writes nothing.
        rows[0].max_u = 1.0f64.next_up();
        put(&scratch, "cand", &rows)?;
        let (baseline, page) = (
            read(&scratch, BASELINE_OF_CAND)?,
            read(&scratch, evidence::PATH)?,
        );
        for mode in ["", "bless"] {
            let r = run(mode)?;
            let dominated: Vec<String> = r
                .failures
                .into_iter()
                .filter(|f| f.bar == Bar::Domination)
                .map(|f| f.text)
                .collect();
            let want = "so3_exp/theta:1e-3/f64: 1.0000000000000002e0 u, over `o1` (o1@1) at 1e0 u";
            assert_eq!(dominated, [want], "{mode}");
        }
        assert!(run("bless")?.text.contains("wrote nothing"));
        assert_eq!(read(&scratch, BASELINE_OF_CAND)?, baseline);
        assert_eq!(read(&scratch, evidence::PATH)?, page);
        // An improvement passes the bars and fails the drift check until it is blessed.
        put(&scratch, "cand", &strata("cand", &[0.5, 1.0, 1.0]))?;
        assert!(run("")?.failures.is_empty());
        let r = run("check")?;
        assert!(r.failures.iter().all(|f| f.bar == Bar::Drift) && !r.failures.is_empty());
        assert!(r
            .text
            .contains("conformance/baseline/cand.csv` is not what `--bless` writes (line 2"));
        assert!(run("bless")?.failures.is_empty() && run("check")?.failures.is_empty());
        Ok(())
    }

    #[test]
    fn a_bless_never_deletes_or_narrows_a_stratum_and_writes_nothing_beside_a_failure(
    ) -> Result<(), String> {
        let scratch = world("shrink")?;
        put(&scratch, "cand", &strata("cand", &[1.0, 2.0, 3.0]))?;
        let run = |mode: &str| go(&scratch, "cand", &[], mode);
        assert!(run("bless")?.failures.is_empty());
        let (baseline, page) = (
            read(&scratch, BASELINE_OF_CAND)?,
            read(&scratch, evidence::PATH)?,
        );
        let refused = |r: &Report, why: &str| -> Result<(), String> {
            assert_eq!(bars(r), ["no-regress"], "{why}\n{}", r.text);
            assert!(r.text.contains("wrote nothing"), "{why}");
            assert_eq!(read(&scratch, BASELINE_OF_CAND)?, baseline, "{why}");
            assert_eq!(read(&scratch, evidence::PATH)?, page, "{why}");
            Ok(())
        };
        // A stratum the candidate no longer answers.
        put(&scratch, "cand", &strata("cand", &[1.0, 2.0]))?;
        refused(&run("bless")?, "dropped")?;
        // A stratum over fewer records.
        let mut rows = strata("cand", &[1.0, 2.0, 3.0]);
        rows[0].n = 32;
        put(&scratch, "cand", &rows)?;
        refused(&run("bless")?, "narrowed")?;
        // Over more records it moves the baseline, as a worse maximum does.
        rows[0].n = 128;
        put(&scratch, "cand", &rows)?;
        let r = run("bless")?;
        assert!(
            r.failures.is_empty()
                && r.text
                    .contains("changes the baseline: so3_exp/theta:1e-3/f64: 128 records")
        );
        assert!(read(&scratch, BASELINE_OF_CAND)?.contains("f64,cand,cand@1,128,1e0,"));
        // A coverage failure, beside a candidate that is fine.
        let dir = scratch.0.join(CORPUS);
        let mut ids = coverage::ids_for_tests();
        let gone = ids.remove(0);
        std::fs::remove_file(dir.join(format!("{gone}.jsonl"))).map_err(|e| e.to_string())?;
        std::fs::write(dir.join("MANIFEST.json"), manifest(&ids)).map_err(|e| e.to_string())?;
        let (baseline, page) = (
            read(&scratch, BASELINE_OF_CAND)?,
            read(&scratch, evidence::PATH)?,
        );
        rows[1].max_u = 1.0;
        put(&scratch, "cand", &rows)?;
        let r = run("bless")?;
        assert_eq!(bars(&r), ["coverage"]);
        assert!(r.text.contains("wrote nothing"));
        assert_eq!(read(&scratch, BASELINE_OF_CAND)?, baseline);
        assert_eq!(read(&scratch, evidence::PATH)?, page);
        Ok(())
    }

    #[test]
    fn a_hand_edit_of_the_baseline_or_the_candidate_columns_fails_the_check_and_writes_nothing(
    ) -> Result<(), String> {
        let scratch = world("hand-edit")?;
        put(&scratch, "cand", &strata("cand", &[1.0, 2.0, 3.0]))?;
        put(&scratch, "o1", &strata("o1", &[9.0, 9.0, 9.0]))?;
        let run = |mode: &str| go(&scratch, "cand", &["o1"], mode);
        assert!(run("bless")?.failures.is_empty());
        for rel in [evidence::PATH, BASELINE_OF_CAND] {
            let good = read(&scratch, rel)?;
            let edited = good.replacen("2.0", "2.1", 1).replacen("2e0", "2e1", 1);
            assert_ne!(edited, good);
            std::fs::write(scratch.0.join(rel), &edited).map_err(|e| e.to_string())?;
            let r = run("check")?;
            let named = |f: &Failure| {
                f.bar == Bar::Drift && f.text.contains(rel) && f.text.contains("line ")
            };
            assert!(r.failures.iter().any(named), "{}", r.text);
            assert_eq!(read(&scratch, rel)?, edited);
            std::fs::write(scratch.0.join(rel), good).map_err(|e| e.to_string())?;
        }
        assert!(run("check")?.failures.is_empty());
        Ok(())
    }

    #[test]
    fn a_missing_oracle_file_and_a_non_finite_candidate_stop_the_run() -> Result<(), String> {
        let scratch = world("stops")?;
        put(&scratch, "cand", &strata("cand", &[1.0, 2.0, 3.0]))?;
        let e = go(&scratch, "cand", &["o1"], "").err().unwrap_or_default();
        assert!(e.contains("no result rows for oracle `o1`"), "{e}");
        let e = go(&scratch, "cand", &["cand"], "")
            .err()
            .unwrap_or_default();
        assert!(e.contains("is the candidate or named twice"), "{e}");
        // The harness leaves the finite records' maximum in a row with a non-finite one.
        let mut rows = strata("cand", &[1.0, 2.0, 3.0]);
        rows[0].nonfinite = 1;
        put(&scratch, "cand", &rows)?;
        put(&scratch, "o1", &strata("o1", &[9.0, 9.0, 9.0]))?;
        assert_eq!(
            bars(&go(&scratch, "cand", &["o1"], "")?),
            ["no-regress", "non-finite"]
        );
        let r = go(&scratch, "cand", &["o1"], "bless")?;
        assert!(r.text.contains("wrote nothing") && bars(&r) == ["non-finite"]);
        assert!(!scratch.0.join(BASELINE).exists() && !scratch.0.join(evidence::PATH).exists());
        // A row with nothing scored on a function that is scored is a defect, not a skip.
        rows[0].nonfinite = 0;
        rows[0].max_u = f64::NAN;
        put(&scratch, "cand", &rows)?;
        let r = go(&scratch, "cand", &["o1"], "")?;
        assert_eq!(bars(&r), ["no-regress", "unscored"]);
        Ok(())
    }

    #[test]
    fn a_planted_c_fails_domination_against_the_correct_seeded_kernel_on_the_real_corpus(
    ) -> Result<(), String> {
        use crate::conformance::{corpus_dir, evaluate};
        use crate::seeded::{Defect, Seeded};
        use helicoid_linalg::Precision;
        let dir = corpus_dir()?;
        let entries: Vec<_> = corpus::manifest(&dir)?
            .into_iter()
            .filter(|e| e.fn_id == "coeff_c")
            .collect();
        let subjects = vec![
            Seeded::generated().registered(),
            Seeded::planted(Defect::CTwoTermsEarly).registered(),
        ];
        let per_subject = evaluate(&dir, &entries, &subjects, Precision::F64)?;
        let scratch = world("planted-c")?;
        for (s, rows) in subjects.iter().zip(&per_subject) {
            put(&scratch, s.subject.name(), rows)?;
        }
        let renamed = |r: &Row| Row {
            subject: "twin".into(),
            ..r.clone()
        };
        put(
            &scratch,
            "twin",
            &per_subject[0].iter().map(renamed).collect::<Vec<_>>(),
        )?;
        // The correct kernel ties its own rows everywhere; only the baseline is missing.
        let r = go(&scratch, "seeded:correct", &["twin"], "")?;
        let tied = per_subject[0].len();
        assert_eq!(bars(&r), ["no-regress"]);
        assert!(r
            .text
            .contains(&format!("{tied} strata paired with an oracle")));
        // The planted `c` ties where its series is the correct kernel's and is beaten from
        // `theta:1e-8`, whose `z` is the first at its switch, to `theta:1e-1`, and on `theta:dense`.
        let r = go(&scratch, "seeded:c-two-terms-1e-8", &["seeded:correct"], "")?;
        let beaten: Vec<&str> = r
            .text
            .lines()
            .filter_map(|l| l.strip_prefix("FAIL domination: coeff_c/"))
            .filter_map(|l| l.split_once("/f64: ").map(|(stratum, _)| stratum))
            .collect();
        let mut want: Vec<String> = (1..=8).rev().map(|k| format!("theta:1e-{k}")).collect();
        want.push("theta:dense".into());
        assert_eq!(beaten, want, "{}", r.text);
        assert!(beaten.len() < tied);
        Ok(())
    }

    #[test]
    fn a_subject_one_ulp_worse_fails_no_regress_where_its_max_rose_on_the_real_corpus(
    ) -> Result<(), String> {
        use crate::conformance::subject::Registered;
        use crate::conformance::testkit::Perfect;
        use crate::conformance::{corpus_dir, evaluate};
        use helicoid_linalg::Precision;
        let dir = corpus_dir()?;
        let entries: Vec<_> = corpus::manifest(&dir)?
            .into_iter()
            .filter(|e| e.fn_id == "so2_exp")
            .collect();
        let subjects = vec![
            Registered::new("1", Box::new(Perfect::exact())),
            Registered::new("1", Box::new(Perfect::next_up())),
        ];
        let [exact, worse] =
            <[Vec<Row>; 2]>::try_from(evaluate(&dir, &entries, &subjects, Precision::F64)?)
                .map_err(|_| "two subjects")?;
        // The same subject as it was blessed, then one ulp worse in every answer.
        let as_cand = |rows: &[Row]| -> Vec<Row> {
            let rename = |r: &Row| Row {
                subject: "cand".into(),
                ..r.clone()
            };
            rows.iter().map(rename).collect()
        };
        let scratch = world("one-ulp")?;
        put(&scratch, "cand", &as_cand(&exact))?;
        assert!(go(&scratch, "cand", &[], "bless")?.failures.is_empty());
        assert!(go(&scratch, "cand", &[], "check")?.failures.is_empty());
        // Real rows reach the baseline and the page: a line per stratum in each.
        let baseline = read(&scratch, BASELINE_OF_CAND)?;
        assert_eq!(baseline.lines().count(), exact.len() + 1);
        let page = read(&scratch, evidence::PATH)?;
        let on_page = page.lines().filter(|l| l.starts_with("| `so2_exp` | `"));
        assert_eq!(on_page.count(), exact.len());
        put(&scratch, "cand", &as_cand(&worse))?;
        let r = go(&scratch, "cand", &[], "")?;
        let rose = worse
            .iter()
            .zip(&exact)
            .filter(|(w, e)| w.max_u > e.max_u)
            .count();
        assert!(
            rose > 0 && r.failures.len() == rose,
            "{rose} of {}\n{}",
            exact.len(),
            r.text
        );
        assert!(r.failures.iter().all(|f| f.bar == Bar::NoRegress));
        Ok(())
    }
}
