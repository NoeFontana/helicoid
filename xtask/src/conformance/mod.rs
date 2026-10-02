//! `cargo xtask conformance [--subject NAME] [--fn ID]` (`docs/PHASE1.md` §5): runs the in-process
//! subjects over the committed corpus, scores each record with the forward error of
//! `docs/NUMERICS.md` §11 in units of `u`, aggregates per `(fn, stratum, precision, subject)`, writes
//! `conformance/results/<subject>.csv` and prints the table by `max_u` descending. It fails on any
//! non-finite output, and when nothing was scored. A run reads the whole corpus (each file parses,
//! ids are line numbers, record counts match `MANIFEST.json`); a run that scored nothing is not a
//! pass. `--self-test` runs the seeded kernels and defects instead (`selftest`).
//!
//! Exact by construction: inputs are hex floats, references are 30-digit decimals, and the error
//! is formed in integers ([`metric`]), never through `f64` parsing of a reference. The output is a
//! function of the corpus, the subject and `git_rev`: no hash order, no threads, no timestamps.
//!
//! Where the specs are silent, the smallest reading was taken (details in [`metric`]):
//!
//! - **p99** is nearest-rank over the records with finite outputs, `⌈0.99·m⌉`-th smallest.
//! - **Non-finite** records count in `n` and `nonfinite`, not in `max_u`, `p99_u` or `argmax_id`;
//!   so does a finite output whose error overflows binary64.
//! - **Backward-only ids** (`so3_from_matrix`) are still evaluated: shape and finiteness are
//!   checked, `n` and `nonfinite` are counted, `max_u` and `p99_u` are `NaN` until backward error
//!   exists.
//! - **A record's score** is the largest error over its output fields; a `--fn` run writes
//!   `<subject>--<fn>.csv` so it never replaces a full result.
//! - **`git_rev`** is `HEAD`, plus `-dirty` when the tree has uncommitted changes.
//! - **Not implemented**: backward error (`Log` near π, `from_matrix`), `f32`, oracle runners, the
//!   envelope, the `helicoid` subject, and the `Dual` comparison of the planted `Q` defect.

pub(crate) mod corpus;
pub(crate) mod metric;
pub(crate) mod number;
mod report;
mod selftest;
mod selftest_se3;
mod selftest_so3;
pub(crate) mod subject;

#[cfg(test)]
mod sanity;
#[cfg(test)]
pub(crate) mod testkit;

use std::path::{Path, PathBuf};
use std::process::Command;

use helicoid_linalg::Precision;

use metric::Rule;
use report::{Aggregate, Row};
use subject::Registered;

const USAGE: &str = "usage: cargo xtask conformance [--subject NAME] [--fn ID] | --self-test";

#[derive(Default)]
struct Options {
    subject: Option<String>,
    fn_id: Option<String>,
    self_test: bool,
}

fn parse_args(args: &[String]) -> Result<Options, String> {
    let mut options = Options::default();
    let mut rest = args.iter();
    while let Some(arg) = rest.next() {
        let (flag, inline) = match arg.split_once('=') {
            Some((f, v)) => (f, Some(v.to_string())),
            None => (arg.as_str(), None),
        };
        let slot = match (flag, &inline) {
            ("--subject", _) => &mut options.subject,
            ("--fn", _) => &mut options.fn_id,
            ("--self-test", None) => {
                options.self_test = true;
                continue;
            }
            _ => return Err(format!("unknown argument `{arg}`; {USAGE}")),
        };
        let value = inline.or_else(|| rest.next().cloned());
        *slot = Some(value.ok_or_else(|| format!("`{flag}` needs a value; {USAGE}"))?);
    }
    if options.self_test && (options.subject.is_some() || options.fn_id.is_some()) {
        return Err(format!("`--self-test` takes no other argument; {USAGE}"));
    }
    Ok(options)
}

pub(crate) fn root() -> Result<PathBuf, String> {
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    manifest
        .parent()
        .map(Path::to_path_buf)
        .ok_or_else(|| "xtask has no parent dir".into())
}

pub(crate) fn corpus_dir() -> Result<PathBuf, String> {
    Ok(root()?.join("conformance/corpus"))
}

/// Every subject over every entry; one row list per subject, in `subjects` order. Each corpus
/// file is read (and so validated) once, whether or not a subject supports it.
fn evaluate(
    dir: &Path,
    entries: &[corpus::Entry],
    subjects: &[Registered],
    precision: Precision,
) -> Result<Vec<Vec<Row>>, String> {
    evaluate_by(dir, entries, subjects, precision, &metric::rule)
}

/// [`evaluate`] with the rule of each function id chosen by `rule_of`.
fn evaluate_by(
    dir: &Path,
    entries: &[corpus::Entry],
    subjects: &[Registered],
    precision: Precision,
    rule_of: &dyn Fn(&str) -> Option<&'static Rule>,
) -> Result<Vec<Vec<Row>>, String> {
    metric::unit_bits(precision)?;
    let mut rows: Vec<Vec<Row>> = subjects.iter().map(|_| Vec::new()).collect();
    for entry in entries {
        let fn_id = entry.fn_id.as_str();
        let rule = rule_of(fn_id).ok_or_else(|| {
            format!("no metric rule for `{fn_id}`: add its row to `metric::TABLE`")
        })?;
        let records = corpus::read(dir, entry)?;
        for (s, out) in subjects.iter().zip(&mut rows) {
            if !s.subject.supports(fn_id) {
                continue;
            }
            let mut aggregate = Aggregate::default();
            for record in &records {
                let output = s.subject.eval(fn_id, record, precision);
                let score = rule.score(record, &output, precision).map_err(|e| {
                    format!("{} on {fn_id} record {}: {e}", s.subject.name(), record.id)
                })?;
                aggregate.add(record.id, &record.stratum, score);
            }
            out.extend(aggregate.rows(fn_id, precision, s.subject.name(), &s.version));
        }
    }
    Ok(rows)
}

fn git_rev(root: &Path) -> String {
    let git = |args: &[&str]| {
        let out = Command::new("git")
            .arg("-C")
            .arg(root)
            .args(args)
            .output()
            .ok();
        out.filter(|o| o.status.success())
            .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
    };
    match (git(&["rev-parse", "HEAD"]), git(&["status", "--porcelain"])) {
        (Some(rev), Some(status)) if status.is_empty() => rev,
        (Some(rev), _) => format!("{rev}-dirty"),
        (None, _) => "unknown".into(),
    }
}

/// `<subject>.csv`, with a `--fn` run apart: `<subject>--<fn>.csv`.
fn results_path(results: &Path, subject: &str, fn_id: Option<&str>) -> PathBuf {
    let clean = |s: &str| {
        s.replace(
            |c: char| !(c.is_ascii_alphanumeric() || "_.-".contains(c)),
            "-",
        )
    };
    results.join(match fn_id {
        Some(f) => format!("{}--{}.csv", clean(subject), clean(f)),
        None => format!("{}.csv", clean(subject)),
    })
}

/// Selects by `options`, evaluates, writes one CSV per subject under `results` and prints its
/// table. The files are written before the verdict, so a failing run leaves its rows to read.
#[allow(clippy::print_stdout)]
fn run_with(
    corpus_dir: &Path,
    results: &Path,
    rev: &str,
    options: &Options,
    mut subjects: Vec<Registered>,
) -> Result<(), String> {
    let mut entries = corpus::manifest(corpus_dir)?;
    if let Some(id) = &options.fn_id {
        entries.retain(|e| e.fn_id == *id);
        if entries.is_empty() {
            return Err(format!("no corpus file for `--fn {id}`"));
        }
    }
    match &options.subject {
        Some(name) => {
            subjects.retain(|s| s.subject.name() == name);
            if subjects.is_empty() {
                return Err(format!("no in-process subject `{name}` is registered"));
            }
        }
        None => subjects.retain(|s| !s.planted),
    }
    let per_subject = evaluate(corpus_dir, &entries, &subjects, Precision::F64)?;
    if per_subject.iter().all(Vec::is_empty) {
        return Err(if subjects.is_empty() {
            format!(
                "nothing was scored: no in-process subject is registered ({} corpus files read, \
                 record counts match MANIFEST.json)",
                entries.len()
            )
        } else {
            "nothing was scored: no selected subject supports the selected function".into()
        });
    }
    let mut failing = 0;
    for (s, rows) in subjects.iter().zip(&per_subject) {
        let path = results_path(results, s.subject.name(), options.fn_id.as_deref());
        std::fs::create_dir_all(results).map_err(|e| format!("{}: {e}", results.display()))?;
        std::fs::write(&path, report::csv(rows, rev))
            .map_err(|e| format!("{}: {e}", path.display()))?;
        println!(
            "{} {}\n{}",
            s.subject.name(),
            s.version,
            report::table(rows)
        );
        eprintln!("conformance: wrote {}", path.display());
        failing += rows.iter().map(|r| r.nonfinite).sum::<usize>();
    }
    match failing {
        0 => Ok(()),
        n => Err(format!("{n} record(s) with a non-finite output")),
    }
}

pub(crate) fn run(args: &[String]) -> Result<(), String> {
    let root = root()?;
    let results = root.join("conformance/results");
    let options = parse_args(args)?;
    if options.self_test {
        return selftest::run(&corpus_dir()?);
    }
    let rev = git_rev(&root);
    run_with(
        &corpus_dir()?,
        &results,
        &rev,
        &options,
        subject::registry(),
    )
}

#[cfg(test)]
mod tests {
    use super::subject::{Output, Subject};
    use super::*;
    use testkit::{Fixed, Perfect, Scratch};

    fn args(a: &[&str]) -> Vec<String> {
        a.iter().map(ToString::to_string).collect()
    }

    #[test]
    fn arguments_take_both_spellings_and_refuse_the_rest() -> Result<(), String> {
        let o = parse_args(&args(&["--subject", "x", "--fn=so2_exp"]))?;
        assert_eq!(
            (o.subject.as_deref(), o.fn_id.as_deref()),
            (Some("x"), Some("so2_exp"))
        );
        for bad in [
            &["--subject"][..],
            &["--nope"],
            &["--self-test=x"],
            &["--self-test", "--fn", "so2_exp"],
            &["--subject=x", "--self-test"],
            &["x"],
        ] {
            assert!(parse_args(&args(bad)).is_err(), "{bad:?}");
        }
        assert!(parse_args(&args(&["--self-test"]))?.self_test);
        Ok(())
    }

    #[test]
    fn a_full_run_and_a_filtered_one_write_different_files() {
        let results = Path::new("/r");
        assert_eq!(
            results_path(results, "seeded:b", None),
            Path::new("/r/seeded-b.csv")
        );
        let p = results_path(results, "helicoid", Some("so3_exp"));
        assert_eq!(p, Path::new("/r/helicoid--so3_exp.csv"));
    }

    #[test]
    fn f32_is_refused_even_with_nothing_to_run() -> Result<(), String> {
        let dir = corpus_dir()?;
        let e = evaluate(&dir, &[], &[], Precision::F32)
            .err()
            .unwrap_or_default();
        assert!(e.contains("f32 is not supported"), "{e}");
        Ok(())
    }

    fn registered(subject: impl Subject + 'static) -> Registered {
        Registered::new("1", Box::new(subject))
    }

    fn only(fn_id: &str) -> Result<Vec<corpus::Entry>, String> {
        let all = corpus::manifest(&corpus_dir()?)?;
        Ok(all.into_iter().filter(|e| e.fn_id == fn_id).collect())
    }

    #[test]
    fn a_subject_that_omits_a_field_or_a_value_is_an_error_not_a_score() -> Result<(), String> {
        let (dir, entries) = (corpus_dir()?, only("so2_exp")?);
        let run = |subject| evaluate(&dir, &entries, &[registered(subject)], Precision::F64);
        let e = run(Fixed::new("m", |_| Output::new()))
            .err()
            .unwrap_or_default();
        assert!(e.contains("no `z`"), "{e}");
        let short = Fixed::new("s", |_| [("z".to_string(), vec![1.0])].into());
        let e = run(short).err().unwrap_or_default();
        assert!(e.contains("returned 1 values, the reference has 2"), "{e}");
        Ok(())
    }

    #[test]
    fn a_subject_is_only_asked_for_the_ids_it_supports() -> Result<(), String> {
        let (dir, entries) = (corpus_dir()?, only("so2_exp")?);
        let picky = Fixed::new("p", |_| Output::new()).only("so2_log");
        assert!(evaluate(&dir, &entries, &[registered(picky)], Precision::F64)?[0].is_empty());
        Ok(())
    }

    fn run_in(
        scratch: &Scratch,
        args: &[&str],
        subject: impl Subject + 'static,
    ) -> Result<(), String> {
        let options = parse_args(&args.iter().map(ToString::to_string).collect::<Vec<_>>())?;
        run_with(
            &corpus_dir()?,
            &scratch.0,
            "rev",
            &options,
            vec![registered(subject)],
        )
    }

    #[test]
    fn two_runs_write_identical_bytes() -> Result<(), String> {
        let (a, b) = (Scratch::new("a"), Scratch::new("b"));
        let fns = ["so2_exp", "se2_ad", "sen3_exp_n1", "so3_log"];
        for scratch in [&a, &b] {
            for f in fns {
                run_in(scratch, &["--fn", f], Perfect::next_up())?;
            }
        }
        for f in fns {
            let name = format!("next-up--{f}.csv");
            let read = |s: &Scratch| std::fs::read(s.0.join(&name)).map_err(|e| e.to_string());
            let (x, y) = (read(&a)?, read(&b)?);
            assert!(
                x.starts_with(report::HEADER.as_bytes()) && x.len() > 300,
                "{name}"
            );
            assert_eq!(x, y, "{name}");
        }
        Ok(())
    }

    #[test]
    fn a_non_finite_output_fails_the_run_after_the_rows_are_written() -> Result<(), String> {
        let scratch = Scratch::new("nan");
        let nan = Fixed::new("nan", |r| {
            let mut out = Perfect::exact().eval("", r, Precision::F64);
            if r.id == 7 {
                out.insert("theta".into(), vec![f64::NAN]);
            }
            out
        });
        let e = run_in(&scratch, &["--fn=so2_log"], nan)
            .err()
            .unwrap_or_default();
        assert!(e.contains("1 record(s) with a non-finite output"), "{e}");
        let csv = std::fs::read_to_string(scratch.0.join("nan--so2_log.csv"))
            .map_err(|e| e.to_string())?;
        // Record 7 is in the first stratum: counted (n = 128, nonfinite = 1), not in the max.
        let first = csv.lines().nth(1).unwrap_or_default();
        assert!(
            first.starts_with("so2_log,theta:1e-12,f64,nan,1,128,"),
            "{first}"
        );
        assert!(
            first.ends_with(",1,rev") && !first.contains(",7,1,rev"),
            "{first}"
        );
        let ok = run_in(&scratch, &["--fn", "so2_log"], Perfect::exact());
        assert!(ok.is_ok(), "{ok:?}");
        let unknown = run_in(&scratch, &["--fn", "nope"], Perfect::exact())
            .err()
            .unwrap_or_default();
        assert!(unknown.contains("no corpus file"), "{unknown}");
        let other = run_in(&scratch, &["--subject", "other"], Perfect::exact())
            .err()
            .unwrap_or_default();
        assert!(other.contains("no in-process subject"), "{other}");
        Ok(())
    }

    #[test]
    fn a_run_that_scored_nothing_fails_after_reading_the_corpus() -> Result<(), String> {
        let scratch = Scratch::new("empty");
        let options = parse_args(&[])?;
        let e = run_with(&corpus_dir()?, &scratch.0, "rev", &options, Vec::new())
            .err()
            .unwrap_or_default();
        assert!(
            e.contains("nothing was scored: no in-process subject"),
            "{e}"
        );
        assert!(e.contains("record counts match"), "{e}");
        assert!(!scratch.0.exists());
        let e = run_with(&scratch.0, &scratch.0, "rev", &options, Vec::new())
            .err()
            .unwrap_or_default();
        assert!(e.contains("MANIFEST.json"), "{e}");
        // A subject that supports nothing selected leaves nothing scored either.
        let picky = Fixed::new("p", |_| Output::new()).only("so2_log");
        let e = run_in(&scratch, &["--fn", "so2_exp"], picky)
            .err()
            .unwrap_or_default();
        assert!(e.contains("no selected subject supports"), "{e}");
        assert!(!scratch.0.exists());
        Ok(())
    }

    #[test]
    fn a_plain_run_skips_planted_subjects_and_naming_one_runs_it() -> Result<(), String> {
        let scratch = Scratch::new("planted");
        let run = |args: &[&str]| {
            let args: Vec<String> = args.iter().map(ToString::to_string).collect();
            let all = subject::registry();
            run_with(&corpus_dir()?, &scratch.0, "rev", &parse_args(&args)?, all)
        };
        let ok = run(&["--fn", "coeff_k"]);
        assert!(ok.is_ok(), "{ok:?}");
        let (correct, planted) = (
            "seeded-correct--coeff_k.csv",
            "seeded-k-sqrt-unsafe--coeff_k.csv",
        );
        assert!(scratch.0.join(correct).exists() && !scratch.0.join(planted).exists());
        // The planted `k` has a NaN derivative at `theta:exact0`: named, it fails like any subject.
        let e = run(&["--subject", "seeded:k-sqrt-unsafe", "--fn", "coeff_k"])
            .err()
            .unwrap_or_default();
        assert!(e.contains("non-finite output"), "{e}");
        assert!(scratch.0.join(planted).exists());
        Ok(())
    }

    #[test]
    fn a_backward_only_id_is_evaluated_and_a_non_finite_output_fails_it() -> Result<(), String> {
        let scratch = Scratch::new("from-matrix");
        let ok = run_in(&scratch, &["--fn", "so3_from_matrix"], Perfect::exact());
        assert!(ok.is_ok(), "{ok:?}");
        let csv = std::fs::read_to_string(scratch.0.join("perfect--so3_from_matrix.csv"))
            .map_err(|e| e.to_string())?;
        assert!(
            csv.lines().count() > 5 && csv.contains(",NaN,NaN,"),
            "{csv}"
        );
        let nan = Fixed::new("nan", |r| {
            let mut out = Perfect::exact().eval("", r, Precision::F64);
            out.insert("q".into(), vec![f64::NAN; 4]);
            out
        });
        let e = run_in(&scratch, &["--fn", "so3_from_matrix"], nan)
            .err()
            .unwrap_or_default();
        assert!(e.contains("non-finite output"), "{e}");
        Ok(())
    }

    #[test]
    fn git_rev_is_head_and_marks_a_dirty_tree() -> Result<(), String> {
        let scratch = Scratch::new("git");
        let io = |e: std::io::Error| e.to_string();
        std::fs::create_dir_all(&scratch.0).map_err(io)?;
        assert_eq!(git_rev(&scratch.0), "unknown");
        let git = |args: &[&str]| {
            let config = ["-c", "user.name=t", "-c", "user.email=t@t"];
            let out = Command::new("git")
                .arg("-C")
                .arg(&scratch.0)
                .args(config)
                .args(["-c", "commit.gpgsign=false"])
                .args(args)
                .output()
                .map_err(io)?;
            match out.status.success() {
                true => Ok(String::from_utf8_lossy(&out.stdout).trim().to_string()),
                false => Err(format!("git {args:?} failed")),
            }
        };
        git(&["init", "-q"])?;
        std::fs::write(scratch.0.join("a"), "a").map_err(io)?;
        git(&["add", "a"])?;
        git(&["commit", "-q", "-m", "m"])?;
        let head = git(&["rev-parse", "HEAD"])?;
        assert_eq!(git_rev(&scratch.0), head);
        std::fs::write(scratch.0.join("b"), "b").map_err(io)?;
        assert_eq!(git_rev(&scratch.0), format!("{head}-dirty"));
        Ok(())
    }
}
