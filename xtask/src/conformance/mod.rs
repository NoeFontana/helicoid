//! `cargo xtask conformance [--subject NAME] [--fn ID] [--precision f64|f32]` (`docs/PHASE1.md` §5): runs the in-process
//! subjects over the committed corpus, scores each record with the forward error of
//! `docs/NUMERICS.md` §11 in units of `u`, aggregates per `(fn, stratum, precision, subject)`, writes
//! `conformance/results/<subject>.csv` and prints the table by `max_u` descending. It fails on any
//! non-finite output (an oracle's is recorded, not failed: `docs/PHASE1.md` §7), and when nothing
//! was scored. A run reads the whole corpus (each file parses,
//! ids are line numbers, record counts match `MANIFEST.json`); a run that scored nothing is not a
//! pass. `--self-test` runs the seeded kernels and defects instead (`selftest`); `--oracle NAME`
//! runs an out-of-process oracle runner and scores its answer files (`oracle`).
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
//! - **`@f32` strata** (`docs/decisions/0016`) are `f32`'s and no other precision's: an `f32` run
//!   scores them alone, with `u = 2^-24` ([`metric`]), and a binary64 run skips them. The inputs an
//!   `f32` subject receives are asserted exactly binary32 (a lossless cast), its outputs likewise.
//!   A subject asked for `f32` on an id with no `@f32` stratum is an error, and so is one with no
//!   `f32` kernel (`Registered::no_f32`). `--precision f32` without `--fn` runs the ids that have
//!   such strata and names the others that a selected subject supports (left out, not asked); with
//!   `--fn` it runs the id named or fails. Its result file is `<subject>[--<fn>]--f32.csv`, the
//!   `precision` column `f32`.
//! - **`helicoid`** (`crate::shipped`) is a plain subject over the eight `coeff_*` ids at both
//!   precisions; a plain run of any of them prints its rows beside `seeded:correct`'s.
//! - **Not implemented**: backward error (`Log` near π, `from_matrix`), the container oracle
//!   runners, `helicoid` on any other id, and the `Dual` comparison of the planted `Q` defect. The
//!   envelope reads the CSVs this module writes (`crate::envelope`).

pub(crate) mod corpus;
pub(crate) mod metric;
pub(crate) mod number;
mod oracle;
pub(crate) mod report;
mod selftest;
mod selftest_envelope;
mod selftest_linalg;
mod selftest_se3;
mod selftest_so3;
pub(crate) mod subject;

#[cfg(test)]
mod committed_linalg;
#[cfg(test)]
mod measure_geodesic;
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

const USAGE: &str =
    "usage: cargo xtask conformance [--subject NAME | --oracle NAME] [--fn ID] [--precision f64|f32] | --self-test";

#[derive(Default)]
struct Options {
    subject: Option<String>,
    oracle: Option<String>,
    fn_id: Option<String>,
    /// `f64` when absent.
    precision: Option<Precision>,
    self_test: bool,
}

fn parse_precision(name: &str) -> Result<Precision, String> {
    match name {
        "f64" => Ok(Precision::F64),
        "f32" => Ok(Precision::F32),
        _ => Err(format!("unknown precision `{name}`; {USAGE}")),
    }
}

fn parse_args(args: &[String]) -> Result<Options, String> {
    let mut options = Options::default();
    let mut rest = args.iter();
    while let Some(arg) = rest.next() {
        let (flag, inline) = match arg.split_once('=') {
            Some((f, v)) => (f, Some(v.to_string())),
            None => (arg.as_str(), None),
        };
        if flag == "--self-test" && inline.is_none() {
            options.self_test = true;
            continue;
        }
        if !matches!(flag, "--subject" | "--oracle" | "--fn" | "--precision") {
            return Err(format!("unknown argument `{arg}`; {USAGE}"));
        }
        let value = inline.or_else(|| rest.next().cloned());
        let value = value.ok_or_else(|| format!("`{flag}` needs a value; {USAGE}"))?;
        match flag {
            "--subject" => options.subject = Some(value),
            "--oracle" => options.oracle = Some(value),
            "--fn" => options.fn_id = Some(value),
            _ => options.precision = Some(parse_precision(&value)?),
        }
    }
    let other = options.subject.is_some()
        || options.oracle.is_some()
        || options.fn_id.is_some()
        || options.precision.is_some();
    if options.self_test && other {
        return Err(format!("`--self-test` takes no other argument; {USAGE}"));
    }
    if options.subject.is_some() && options.oracle.is_some() {
        return Err(format!(
            "`--subject` and `--oracle` exclude each other; {USAGE}"
        ));
    }
    Ok(options)
}

/// The names of the oracle runners, the subjects the envelope's domination bar is over.
pub(crate) fn oracle_names() -> Vec<&'static str> {
    oracle::names()
}

pub(crate) use oracle::{backend as oracle_backend, Backend};

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
/// file is read (and so validated) once, whether or not a subject supports it. At `f32` only the
/// `@f32` strata are scored, and a supporting subject is an error on an entry that has none, and
/// where it has no `f32` kernel.
pub(crate) fn evaluate(
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
            let f32 = precision == Precision::F32;
            if let (true, Some(why)) = (f32, &s.no_f32) {
                return Err(format!("{} has no f32 kernel: {why}", s.subject.name()));
            }
            let scored = records.iter().filter(|r| r.is_f32_stratum() == f32);
            let mut scored = scored.peekable();
            if f32 && scored.peek().is_none() {
                return Err(format!(
                    "{} was asked for f32 on `{fn_id}`, which has no `@f32` stratum: only the \
                     scalar coefficient ids have any (docs/decisions/0016 item 2)",
                    s.subject.name()
                ));
            }
            let mut aggregate = Aggregate::default();
            for record in scored {
                if f32 {
                    record
                        .require_binary32()
                        .map_err(|e| format!("{fn_id} record {}: {e}", record.id))?;
                }
                let output = s.subject.eval(fn_id, record, precision);
                let score = rule.score(record, &output, precision).map_err(|e| {
                    format!("{} on {fn_id} record {}: {e}", s.subject.name(), record.id)
                })?;
                aggregate.add(record.id, &record.stratum, score);
            }
            let version = s.version_at(precision);
            out.extend(aggregate.rows(fn_id, precision, s.subject.name(), version));
        }
    }
    Ok(rows)
}

/// The entries an `f32` run without `--fn` scores, those that have `@f32` strata, and the ids of the
/// others that a subject supports: left out, and named, so the run does not skip them in silence.
fn split_f32(
    dir: &Path,
    entries: Vec<corpus::Entry>,
    subjects: &[Registered],
) -> Result<(Vec<corpus::Entry>, Vec<String>), String> {
    let (mut with, mut left_out) = (Vec::new(), Vec::new());
    for e in entries {
        if corpus::mentions_f32(dir, &e)? {
            with.push(e);
        } else if subjects.iter().any(|s| s.subject.supports(&e.fn_id)) {
            left_out.push(e.fn_id);
        }
    }
    Ok((with, left_out))
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

/// `<subject>.csv`, with a `--fn` run apart, `<subject>--<fn>.csv`, and an `f32` run apart again,
/// `<subject>[--<fn>]--f32.csv`.
pub(crate) fn results_path(
    results: &Path,
    subject: &str,
    fn_id: Option<&str>,
    precision: Precision,
) -> PathBuf {
    let clean = |s: &str| {
        s.replace(
            |c: char| !(c.is_ascii_alphanumeric() || "_.-".contains(c)),
            "-",
        )
    };
    let name = match fn_id {
        Some(f) => format!("{}--{}", clean(subject), clean(f)),
        None => clean(subject),
    };
    results.join(match precision {
        Precision::F64 => format!("{name}.csv"),
        Precision::F32 => format!("{name}--f32.csv"),
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
    let precision = options.precision.unwrap_or(Precision::F64);
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
    if precision == Precision::F32 && options.fn_id.is_none() {
        let (with, left_out) = split_f32(corpus_dir, entries, &subjects)?;
        eprintln!(
            "conformance: f32 runs the {} ids that have @f32 strata",
            with.len()
        );
        if !left_out.is_empty() {
            eprintln!(
                "conformance: f32 leaves out (no @f32 stratum, supported by a selected subject): {}",
                left_out.join(", ")
            );
        }
        entries = with;
    }
    let per_subject = evaluate(corpus_dir, &entries, &subjects, precision)?;
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
        let fn_id = options.fn_id.as_deref();
        let path = results_path(results, s.subject.name(), fn_id, precision);
        std::fs::create_dir_all(results).map_err(|e| format!("{}: {e}", results.display()))?;
        std::fs::write(&path, report::csv(rows, rev))
            .map_err(|e| format!("{}: {e}", path.display()))?;
        println!(
            "{} {}\n{}",
            s.subject.name(),
            s.version_at(precision),
            report::table(rows)
        );
        eprintln!("conformance: wrote {}", path.display());
        let nonfinite = rows.iter().map(|r| r.nonfinite).sum::<usize>();
        match (nonfinite, options.oracle.is_some()) {
            (0, _) => {}
            // Oracles may be wrong: the count is in the row, the run still passes (§7).
            (n, true) => eprintln!(
                "conformance: {} has {n} record(s) with a non-finite output, recorded",
                s.subject.name()
            ),
            (n, false) => failing += n,
        }
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
    let corpus = corpus_dir()?;
    let subjects = match &options.oracle {
        Some(name) => vec![oracle::run(
            name,
            &root,
            &corpus,
            &results,
            options.fn_id.as_deref(),
        )?],
        None => subject::registry(),
    };
    run_with(&corpus, &results, &rev, &options, subjects)
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
        let p = |a: &[&str]| parse_args(&args(a)).map(|o| o.precision);
        assert_eq!(p(&["--precision", "f32"]), Ok(Some(Precision::F32)));
        assert_eq!(p(&["--precision=f64"]), Ok(Some(Precision::F64)));
        assert_eq!(p(&[]), Ok(None));
        let o = parse_args(&args(&["--oracle=tf_tree_math", "--fn", "so2_exp"]))?;
        assert_eq!(o.oracle.as_deref(), Some("tf_tree_math"));
        for bad in [
            &["--subject"][..],
            &["--nope"],
            &["--self-test=x"],
            &["--self-test", "--fn", "so2_exp"],
            &["--subject=x", "--self-test"],
            &["--precision", "f32", "--self-test"],
            &["--precision"],
            &["--precision=f16"],
            &["--oracle", "x", "--self-test"],
            &["--oracle", "x", "--subject", "y"],
            &["--oracle"],
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
        let at = |f, p| results_path(results, "seeded:b", f, p);
        assert_eq!(at(None, Precision::F64), Path::new("/r/seeded-b.csv"));
        let p = results_path(results, "helicoid", Some("so3_exp"), Precision::F64);
        assert_eq!(p, Path::new("/r/helicoid--so3_exp.csv"));
        assert_eq!(at(None, Precision::F32), Path::new("/r/seeded-b--f32.csv"));
        let p = at(Some("coeff_k"), Precision::F32);
        assert_eq!(p, Path::new("/r/seeded-b--coeff_k--f32.csv"));
    }

    #[test]
    fn a_subject_asked_for_f32_on_an_id_without_f32_strata_is_an_error() -> Result<(), String> {
        let (dir, entries) = (corpus_dir()?, only("so3_exp")?);
        let asked = |s: Registered| evaluate(&dir, &entries, &[s], Precision::F32);
        let e = asked(registered(Perfect::exact()))
            .err()
            .unwrap_or_default();
        assert!(e.contains("asked for f32 on `so3_exp`"), "{e}");
        assert!(e.contains("no `@f32` stratum"), "{e}");
        // Not asked: a subject that does not support the id, and nothing to run.
        let picky = Fixed::new("p", |_| Output::new()).only("so2_log");
        assert!(asked(registered(picky))?[0].is_empty());
        assert!(evaluate(&dir, &[], &[], Precision::F32)?.is_empty());
        Ok(())
    }

    #[test]
    fn an_f32_input_that_is_not_a_binary32_is_an_error_not_a_rounding() -> Result<(), String> {
        let scratch = Scratch::new("f32-input");
        let io = |e: std::io::Error| e.to_string();
        std::fs::create_dir_all(&scratch.0).map_err(io)?;
        let manifest = r#"{"files":{"coeff_k.jsonl":{"kind":"corpus","records":1}}}"#;
        std::fs::write(scratch.0.join("MANIFEST.json"), manifest).map_err(io)?;
        // 1 + 2^-52 is a binary64 and no binary32; 1 + 2^-23 is both.
        for (theta, lossless) in [("0x1.0000000000001p+0", false), ("0x1.00000p+0", true)] {
            let record = format!(
                r#"{{"id":0,"in":{{"theta":["{theta}"]}},"out":{{"d_branch":["1e0"],"value":["1e0"]}},"stratum":"s@f32"}}"#
            );
            std::fs::write(scratch.0.join("coeff_k.jsonl"), record).map_err(io)?;
            let entries = corpus::manifest(&scratch.0)?;
            let subject = [registered(Perfect::exact())];
            let run = evaluate(&scratch.0, &entries, &subject, Precision::F32);
            assert_eq!(run.is_ok(), lossless, "{theta}");
            let e = run.err().unwrap_or_default();
            let named = "coeff_k record 0: input `theta` holds";
            assert_eq!(e.contains(named), !lossless, "{e}");
        }
        Ok(())
    }

    #[test]
    fn an_f32_run_scores_the_f32_strata_alone_and_writes_its_own_file() -> Result<(), String> {
        let scratch = Scratch::new("f32-run");
        run_in(&scratch, &["--precision=f32"], Perfect::exact())?;
        let path = scratch.0.join("perfect--f32.csv");
        let csv = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
        let rows: Vec<Vec<&str>> = csv
            .lines()
            .skip(1)
            .map(|l| l.split(',').collect())
            .collect();
        // The 9 coefficient ids' 28 each (`coeff_alpha`, `0062`, the ninth) and `coeff_r`'s
        // `q:w0@f32`, then `0056`'s: `solve_cubic`
        // 12, `eig3` 19, `chol_n*` 6 each, `chol_solve_n*` 5 each, `quat_renormalize` 5, `real_sqrt`
        // and `real_cbrt` 8 each, `real_sin_cos` 5, `real_acos` 4, `real_atan2` 7, `real_div` 2;
        // `0061`'s `mat2_inverse_adj` 3.
        let new = 12 + 19 + 2 * 6 + 2 * 5 + 5 + 2 * 8 + 5 + 4 + 7 + 2 + 3;
        assert_eq!(rows.len(), 9 * 28 + 1 + new);
        assert!(rows.iter().all(|c| c[1].ends_with("@f32") && c[2] == "f32"));
        // An id named without any is an error, not an empty run.
        let named = ["--precision", "f32", "--fn", "so3_exp"];
        let e = run_in(&scratch, &named, Perfect::exact()).err();
        assert!(e.unwrap_or_default().contains("asked for f32 on `so3_exp`"));
        Ok(())
    }

    #[test]
    fn a_subject_with_no_f32_kernel_is_an_error_at_f32_and_a_subject_at_f64() -> Result<(), String>
    {
        let (dir, entries) = (corpus_dir()?, only("coeff_k")?);
        let bare = || Registered {
            no_f32: Some("no kernel".into()),
            ..registered(Perfect::exact())
        };
        let e = evaluate(&dir, &entries, &[bare()], Precision::F32);
        assert_eq!(
            e.err().unwrap_or_default(),
            "perfect has no f32 kernel: no kernel"
        );
        assert!(evaluate(&dir, &entries, &[bare()], Precision::F64).is_ok());
        // Asked only for what it supports: an id it does not is not an f32 request of it.
        let elsewhere = only("so2_exp")?;
        let picky = Registered {
            no_f32: Some("no kernel".into()),
            ..registered(Fixed::new("p", |_| Output::new()).only("coeff_k"))
        };
        assert!(evaluate(&dir, &elsewhere, &[picky], Precision::F32)?[0].is_empty());
        Ok(())
    }

    #[test]
    fn a_planted_subject_the_f32_kernels_do_not_model_is_an_error_at_f32_and_writes_nothing(
    ) -> Result<(), String> {
        let scratch = Scratch::new("planted-c-f32");
        let run = |name: &str, precision: &str, fn_id: Option<&str>| {
            let mut args = vec!["--subject", name, "--precision", precision];
            args.extend(fn_id.into_iter().flat_map(|f| ["--fn", f]));
            let args: Vec<String> = args.iter().map(ToString::to_string).collect();
            let all = subject::registry();
            run_with(&corpus_dir()?, &scratch.0, "rev", &parse_args(&args)?, all)
        };
        // The planted `c` is a candidate of the binary64 sweep: at `f32` it would answer as the D12
        // kernel under its own name, on one coefficient id or on every one that has `@f32` strata.
        let c = "seeded:c-two-terms-1e-8";
        for fn_id in [Some("coeff_c"), Some("coeff_b"), None] {
            let e = run(c, "f32", fn_id).err().unwrap_or_default();
            assert!(e.starts_with(&format!("{c} has no f32 kernel")), "{e}");
        }
        assert!(!scratch.0.exists());
        // At binary64 it runs; a defect the `f32` kernels do model gets past the refusal, and fails
        // on the non-finite output its own mechanism detects.
        let ok = run(c, "f64", Some("coeff_c"));
        assert!(ok.is_ok(), "{ok:?}");
        let e = run("seeded:b-no-series", "f32", Some("coeff_b"));
        assert!(e.err().unwrap_or_default().contains("non-finite output"));
        assert!(scratch
            .0
            .join("seeded-b-no-series--coeff_b--f32.csv")
            .exists());
        Ok(())
    }

    #[test]
    fn an_f32_run_without_fn_names_the_ids_it_leaves_out() -> Result<(), String> {
        let dir = corpus_dir()?;
        let split = |subjects: &[Registered]| split_f32(&dir, corpus::manifest(&dir)?, subjects);
        let (with, left_out) = split(&[registered(Perfect::exact())])?;
        // The 9 coefficient ids (`0016`, `coeff_alpha` of `0062`), `0056`'s 13 and `0061`'s one.
        assert_eq!(with.len(), 9 + 13 + 1);
        let of_0056 = |id: &str| {
            ["solve_cubic", "eig3", "chol", "quat_", "real_", "mat2_"]
                .iter()
                .any(|p| id.starts_with(p))
        };
        assert!(with
            .iter()
            .all(|e| e.fn_id.starts_with("coeff_") || of_0056(&e.fn_id)));
        let total = corpus::manifest(&dir)?.len();
        assert_eq!(left_out.len(), total - 9 - 13 - 1);
        assert!(["so3_exp", "sen3_exp_n1", "so2_exp"]
            .iter()
            .all(|id| left_out.iter().any(|l| l == id)));
        // Only what a selected subject supports is named.
        let of = |id: &str| registered(Fixed::new("p", |_| Output::new()).only(id));
        assert_eq!(split(&[of("so3_log")])?.1, ["so3_log"]);
        assert!(split(&[of("coeff_k")])?.1.is_empty());
        // The correct seeded subject and `helicoid` answer `so3_*` and `sen3_*`, which have no
        // `@f32` stratum. The `ad`, `log` and `*_inv` widths are `helicoid`'s alone: the seeded
        // kernel is `PHASE1.md` §10's `exp`, `jr` and `jl`.
        let plain: Vec<Registered> = subject::registry()
            .into_iter()
            .filter(|r| !r.planted)
            .collect();
        let want = [
            // The SE(3) chart ids (`0060` decision 8), binary64 only like every vector id.
            "se3_decoupled_local",
            "se3_decoupled_retract",
            // The geodesic ids, `PHASE4.md` §4's, which the `helicoid` subject answers at `f64`.
            "se3_geodesic",
            "se3_screw_local",
            "se3_screw_retract",
            "se3_world_local",
            "se3_world_retract",
            "sen3_ad_n1",
            "sen3_ad_n2",
            "sen3_ad_n3",
            "sen3_exp_n1",
            "sen3_exp_n2",
            "sen3_exp_n3",
            "sen3_jl_inv_n1",
            "sen3_jl_inv_n2",
            "sen3_jl_inv_n3",
            "sen3_jl_n1",
            "sen3_jl_n2",
            "sen3_jl_n3",
            "sen3_jr_inv_n1",
            "sen3_jr_inv_n2",
            "sen3_jr_inv_n3",
            "sen3_jr_n1",
            "sen3_jr_n2",
            "sen3_jr_n3",
            "sen3_log_n1",
            "sen3_log_n2",
            "sen3_log_n3",
            // Every `so3_*` id the `helicoid` subject answers. None has an `@f32` stratum until a
            // record extends `0016`, so a plain `f32` run names them all and runs none.
            "so3_act",
            "so3_exp",
            "so3_from_matrix",
            "so3_geodesic",
            "so3_jl",
            "so3_jl_inv",
            "so3_jr",
            "so3_jr_inv",
            "so3_log",
        ];
        assert_eq!(split(&plain)?.1, want);
        Ok(())
    }

    #[test]
    fn a_rows_subject_version_is_the_one_of_the_precision_it_ran_at() -> Result<(), String> {
        let scratch = Scratch::new("version");
        let versions = || Registered {
            version: "v64".into(),
            version_f32: "v32".into(),
            ..registered(Perfect::exact())
        };
        let cells = |file: &str| -> Result<Vec<String>, String> {
            let csv = std::fs::read_to_string(scratch.0.join(file)).map_err(|e| e.to_string())?;
            let cell = |l: &str| l.split(',').nth(4).unwrap_or_default().to_string();
            Ok(csv.lines().skip(1).map(cell).collect())
        };
        let runs = [
            (&["--fn", "coeff_k"][..], "perfect--coeff_k.csv", "v64"),
            (
                &["--precision", "f32", "--fn", "coeff_k"],
                "perfect--coeff_k--f32.csv",
                "v32",
            ),
        ];
        for (args, file, want) in runs {
            let args: Vec<String> = args.iter().map(ToString::to_string).collect();
            let options = parse_args(&args)?;
            run_with(
                &corpus_dir()?,
                &scratch.0,
                "rev",
                &options,
                vec![versions()],
            )?;
            let cells = cells(file)?;
            assert!(
                cells.len() > 20 && cells.iter().all(|c| c == want),
                "{file}: {cells:?}"
            );
        }
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

    /// `Perfect::exact` but for a NaN `theta` at record 7, in the first stratum of `so2_log`.
    fn nan_at_record_7() -> Fixed {
        Fixed::new("nan", |r| {
            let mut out = Perfect::exact().eval("", r, Precision::F64);
            if r.id == 7 {
                out.insert("theta".into(), vec![f64::NAN]);
            }
            out
        })
    }

    #[test]
    fn a_non_finite_oracle_answer_is_recorded_and_the_run_passes() -> Result<(), String> {
        let scratch = Scratch::new("oracle-nan");
        run_in(
            &scratch,
            &["--oracle=nan", "--fn=so2_log"],
            nan_at_record_7(),
        )?;
        let csv = std::fs::read_to_string(scratch.0.join("nan--so2_log.csv"))
            .map_err(|e| e.to_string())?;
        let first = csv.lines().nth(1).unwrap_or_default();
        assert!(
            first.starts_with("so2_log,theta:1e-12,f64,nan,1,128,"),
            "{first}"
        );
        Ok(())
    }

    #[test]
    fn a_non_finite_output_fails_the_run_after_the_rows_are_written() -> Result<(), String> {
        let scratch = Scratch::new("nan");
        let nan = nan_at_record_7();
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
