//! `cargo xtask conformance --oracle NAME` (`docs/PHASE1.md` §7): runs the oracle runner
//! `runners/<NAME>` (a workspace-excluded crate with its own lockfile, run by
//! `cargo run --release --locked`) over the corpus and scores what it wrote with the exact metric
//! of an in-process subject, as a [`FileSubject`]. The file protocol, the smallest reading of §7
//! (0014 (draft) question 26):
//!
//! - `<runner> --version` prints the pin on one line: the `subject_version` of every row.
//! - `<runner> --out DIR FILE.jsonl…` writes `DIR/<fn>.jsonl` for each corpus file whose function
//!   it answers, and nothing for the others.
//! - A line is `{"id":N,"out":{"<field>":["<hex float>",…]}}`: `N` the corpus line's `id`, a field
//!   flat in the reference's order (a matrix column-major, no `shape`), a value a hex float or
//!   `nan`, `inf`, `-inf`; a scalar field may be the bare string.
//! - An answer file holds every record of its corpus file, ids in order. A runner that skipped a
//!   record, answered a function the corpus does not have or wrote a malformed line is an error,
//!   not a score; so is a function in [`Runner::answers`] with no file, since a missing file
//!   otherwise reads as "not answered". A wrong answer is a score, and a non-finite one a recorded
//!   row (`nonfinite`), not a failed run: oracles may be wrong.
//!
//! `DIR` is `conformance/results/oracle/<NAME>` (`<NAME>--<fn>` for a `--fn` run, as its CSV),
//! emptied first and left for inspection.

use std::collections::BTreeMap;
use std::ffi::OsString;
use std::path::Path;
use std::process::{Command, ExitStatus, Stdio};

use helicoid_linalg::Precision;
use serde::Deserialize;
use serde_json::Value;

use super::corpus::{self, Entry, Object, Record};
use super::number::parse_hex;
use super::subject::{Output, Registered, Subject};

/// A runner crate under `runners/`; `name` is the subject's name in the result rows.
struct Runner {
    name: &'static str,
    dir: &'static str,
    /// The function ids it must answer: the list its own `--out` advertises, held again here so
    /// that an id it stops answering fails the run instead of dropping out of the results.
    answers: &'static [&'static str],
}

const RUNNERS: &[Runner] = &[
    Runner {
        name: "tf_tree_math",
        dir: "runners/tf_tree_math",
        answers: &["so3_exp", "so3_log", "sen3_exp_n1", "sen3_log_n1"],
    },
    Runner {
        name: "sophus_rs",
        dir: "runners/sophus_rs",
        answers: &[
            "so3_exp",
            "so3_log",
            "so3_jl",
            "so3_jr",
            "so3_jl_inv",
            "so3_jr_inv",
            "sen3_exp_n1",
            "sen3_log_n1",
            "sen3_jl_n1",
            "sen3_jr_n1",
            "sen3_jl_inv_n1",
            "sen3_jr_inv_n1",
        ],
    },
];

pub(super) fn names() -> Vec<&'static str> {
    RUNNERS.iter().map(|r| r.name).collect()
}

/// A runner's answers, read back from its files.
pub(crate) struct FileSubject {
    name: String,
    /// Per function id, the answer to record `i` at index `i` (a record's `id` is its line number).
    answers: BTreeMap<String, Vec<Output>>,
}

impl Subject for FileSubject {
    fn name(&self) -> &str {
        &self.name
    }

    fn supports(&self, fn_id: &str) -> bool {
        self.answers.contains_key(fn_id)
    }

    /// The runner's answer, or none, which the metric reports as a missing field: ingestion has
    /// already checked that no record is missing. Its answers are binary64, so at another
    /// precision there is none.
    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        if precision != Precision::F64 {
            return Output::new();
        }
        let at = usize::try_from(record.id).ok();
        let answer = at.and_then(|i| self.answers.get(fn_id)?.get(i));
        answer.cloned().unwrap_or_default()
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Answer {
    id: u64,
    out: Object,
}

fn float(v: &Value) -> Result<f64, String> {
    match v.as_str().ok_or_else(|| format!("not a string: {v}"))? {
        "nan" => Ok(f64::NAN),
        "inf" => Ok(f64::INFINITY),
        "-inf" => Ok(f64::NEG_INFINITY),
        hex => parse_hex(hex),
    }
}

/// One answer line: the record's `id` and the output fields.
fn parse_answer(line: &str) -> Result<(u64, Output), String> {
    let answer: Answer = serde_json::from_str(line).map_err(|e| e.to_string())?;
    let mut out = Output::new();
    for (key, value) in answer.out.0 {
        let values = match &value {
            Value::Array(items) => items.iter().map(float).collect(),
            v => float(v).map(|x| vec![x]),
        };
        out.insert(key.clone(), values.map_err(|e| format!("`{key}`: {e}"))?);
    }
    Ok((answer.id, out))
}

impl FileSubject {
    /// The answers in `dir`, one `<fn>.jsonl` per function of `entries` the runner answers.
    pub(crate) fn ingest(name: &str, dir: &Path, entries: &[Entry]) -> Result<Self, String> {
        let io = |p: &Path, e: std::io::Error| format!("{}: {e}", p.display());
        let file_of = |e: &Entry| format!("{}.jsonl", e.fn_id);
        for file in std::fs::read_dir(dir).map_err(|e| io(dir, e))? {
            let file = file.map_err(|e| io(dir, e))?.file_name();
            if !entries.iter().any(|e| *file == *file_of(e)) {
                let file = file.to_string_lossy();
                return Err(format!(
                    "{}: `{file}` is not a corpus function",
                    dir.display()
                ));
            }
        }
        let mut answers = BTreeMap::new();
        for entry in entries {
            let path = dir.join(file_of(entry));
            let text = match std::fs::read_to_string(&path) {
                Ok(text) => text,
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => continue,
                Err(e) => return Err(io(&path, e)),
            };
            let mut file = Vec::with_capacity(entry.records);
            for (i, line) in text.lines().enumerate() {
                let at = |e: String| format!("{}:{}: {e}", path.display(), i + 1);
                let (id, out) = parse_answer(line).map_err(at)?;
                if id != i as u64 {
                    return Err(at(format!("id {id} on line {}", i + 1)));
                }
                file.push(out);
            }
            if file.len() != entry.records {
                return Err(format!(
                    "{}: {} answers, the corpus has {} records",
                    path.display(),
                    file.len(),
                    entry.records
                ));
            }
            answers.insert(entry.fn_id.clone(), file);
        }
        Ok(Self {
            name: name.to_string(),
            answers,
        })
    }
}

/// `cargo run --release --locked` on the runner, arguments after `--` still to add.
fn runner_command(root: &Path, runner: &Runner, quiet: bool) -> Command {
    let cargo = std::env::var_os("CARGO").unwrap_or_else(|| OsString::from("cargo"));
    let mut command = Command::new(cargo);
    command
        .current_dir(root)
        .args(["run", "--release", "--locked"])
        .args(quiet.then_some("--quiet"))
        .arg("--manifest-path")
        .arg(root.join(runner.dir).join("Cargo.toml"))
        .arg("--");
    command
}

fn succeeded(what: &str, status: ExitStatus) -> Result<(), String> {
    match status.success() {
        true => Ok(()),
        false => Err(format!("{what} failed: {status}")),
    }
}

/// The pin a runner printed: exactly one non-empty line.
fn version_line(stdout: &[u8]) -> Result<String, String> {
    let text = std::str::from_utf8(stdout).map_err(|e| format!("`--version` is not UTF-8: {e}"))?;
    let line = text.strip_suffix('\n').unwrap_or(text);
    match line.is_empty() || line.contains(['\n', '\r']) {
        true => Err(format!("`--version` must print one line, not {text:?}")),
        false => Ok(line.to_string()),
    }
}

/// Builds and runs the runner `name` over the corpus (only `fn_id`'s file when given) and reads
/// its answers: the subject, its `subject_version` the runner's pin.
pub(crate) fn run(
    name: &str,
    root: &Path,
    corpus_dir: &Path,
    results: &Path,
    fn_id: Option<&str>,
) -> Result<Registered, String> {
    let runner = RUNNERS.iter().find(|r| r.name == name).ok_or_else(|| {
        let known: Vec<&str> = RUNNERS.iter().map(|r| r.name).collect();
        format!("no oracle runner `{name}`; have {known:?}")
    })?;
    let launch = |quiet| runner_command(root, runner, quiet);
    execute(runner, &launch, corpus_dir, results, fn_id)
}

/// `dir`, empty: whatever an earlier run left there is not an answer of this one.
fn fresh_dir(dir: &Path) -> Result<(), String> {
    if let Err(e) = std::fs::remove_dir_all(dir) {
        if e.kind() != std::io::ErrorKind::NotFound {
            return Err(format!("{}: {e}", dir.display()));
        }
    }
    std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))
}

/// [`run`] with the way to start the runner given: `launch(quiet)` is the command up to the
/// arguments of the protocol.
fn execute(
    runner: &Runner,
    launch: &dyn Fn(bool) -> Command,
    corpus_dir: &Path,
    results: &Path,
    fn_id: Option<&str>,
) -> Result<Registered, String> {
    let mut entries = corpus::manifest(corpus_dir)?;
    if let Some(id) = fn_id {
        entries.retain(|e| e.fn_id == id);
        if entries.is_empty() {
            return Err(format!("no corpus file for `--fn {id}`"));
        }
    }
    let version = {
        let out = launch(false)
            .arg("--version")
            .stderr(Stdio::inherit())
            .output();
        let out = out.map_err(|e| format!("{}: {e}", runner.name))?;
        succeeded(&format!("`{} --version`", runner.name), out.status)?;
        version_line(&out.stdout)?
    };
    let out_dir = results.join("oracle").join(match fn_id {
        Some(id) => format!("{}--{id}", runner.name),
        None => runner.name.to_string(),
    });
    fresh_dir(&out_dir)?;
    let files = entries
        .iter()
        .map(|e| corpus_dir.join(format!("{}.jsonl", e.fn_id)));
    let status = launch(true)
        .arg("--out")
        .arg(&out_dir)
        .args(files)
        .status()
        .map_err(|e| format!("{}: {e}", runner.name))?;
    succeeded(&format!("runner `{}`", runner.name), status)?;
    let subject = FileSubject::ingest(runner.name, &out_dir, &entries)?;
    let unanswered = runner
        .answers
        .iter()
        .find(|id| entries.iter().any(|e| e.fn_id == **id) && !subject.supports(id));
    if let Some(id) = unanswered {
        return Err(format!(
            "runner `{}` wrote no answers for `{id}`",
            runner.name
        ));
    }
    Ok(Registered {
        version: version.clone(),
        version_f32: version,
        // `tf_tree_math` is binary64 only (`tf_tree` D6), so a run at `f32` is an error on every id
        // it answers rather than another kernel's answer under its name (`0016` item 2).
        no_f32: Some("tf_tree_math has no f32 kernel".to_string()),
        subject: Box::new(subject),
        planted: false,
    })
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use serde_json::json;

    use super::super::testkit::{self, Perfect, Scratch};
    use super::super::{corpus_dir, evaluate};
    use super::*;

    /// Python's `float.hex()`, the corpus's spelling, for what a subject can return.
    fn hex(x: f64) -> String {
        let bits = x.to_bits();
        let sign = if bits >> 63 == 1 { "-" } else { "" };
        let (biased, frac) = ((bits >> 52) & 0x7ff, bits & ((1 << 52) - 1));
        match (biased, frac) {
            (0, 0) => format!("{sign}0x0.0p+0"),
            (0, f) => format!("{sign}0x0.{f:013x}p-1022"),
            (e, f) => format!("{sign}0x1.{f:013x}p{:+}", e as i64 - 1023),
        }
    }

    fn write_answers(dir: &Path, fn_id: &str, answers: &[Output]) -> Result<(), String> {
        std::fs::create_dir_all(dir).map_err(|e| e.to_string())?;
        let mut text = String::new();
        for (id, out) in answers.iter().enumerate() {
            let fields: BTreeMap<&String, Vec<String>> = out
                .iter()
                .map(|(k, v)| (k, v.iter().map(|&x| hex(x)).collect()))
                .collect();
            text += &format!("{}\n", json!({ "id": id, "out": fields }));
        }
        std::fs::write(dir.join(format!("{fn_id}.jsonl")), text).map_err(|e| e.to_string())
    }

    fn entry(fn_id: &str, records: usize) -> Entry {
        Entry {
            fn_id: fn_id.to_string(),
            records,
        }
    }

    #[test]
    fn an_ingested_answer_scores_exactly_like_the_subject_that_wrote_it() -> Result<(), String> {
        let dir = corpus_dir()?;
        let entries: Vec<Entry> = corpus::manifest(&dir)?
            .into_iter()
            .filter(|e| e.fn_id == "so2_exp")
            .collect();
        let records = corpus::read(&dir, &entries[0])?;
        let scratch = Scratch::new("ingest");
        let writer = Perfect::next_up();
        let answers: Vec<Output> = records
            .iter()
            .map(|r| writer.eval("so2_exp", r, Precision::F64))
            .collect();
        write_answers(&scratch.0, "so2_exp", &answers)?;
        let file = FileSubject::ingest("next-up", &scratch.0, &entries)?;
        assert!(file.supports("so2_exp") && !file.supports("so2_log"));
        let score = |s: Box<dyn Subject>| {
            evaluate(&dir, &entries, &[Registered::new("1", s)], Precision::F64)
        };
        let (from_file, direct) = (score(Box::new(file))?, score(Box::new(writer))?);
        let strata: std::collections::BTreeSet<&str> =
            records.iter().map(|r| r.stratum.as_str()).collect();
        assert_eq!(from_file[0].len(), strata.len());
        assert!(from_file[0].iter().any(|r| r.max_u > 0.0));
        assert_eq!(from_file, direct);
        Ok(())
    }

    #[test]
    fn a_line_holds_an_id_and_hex_floats_bare_strings_and_non_finite_words() -> Result<(), String> {
        let line =
            r#"{"id":3,"out":{"a":"0x1.0000000000000p+0","z":["nan","-inf","inf","-0x0.0p+0"]}}"#;
        let (id, out) = parse_answer(line)?;
        assert_eq!((id, out["a"].clone()), (3, vec![1.0]));
        let z = &out["z"];
        assert!(
            z[0].is_nan() && z[1].is_infinite() && z[1] < 0.0 && z[2].is_infinite() && z[2] > 0.0
        );
        assert_eq!(z[3].to_bits(), (-0.0f64).to_bits());
        Ok(())
    }

    #[test]
    fn a_malformed_line_is_an_error() {
        let bad = |line: &str, needle: &str| {
            let e = parse_answer(line).err().unwrap_or_default();
            assert!(e.contains(needle), "{line}: {e}");
        };
        bad(r#"{"id":0,"out":{"z":["1.5"]}}"#, "bad hex float");
        bad(r#"{"id":0,"out":{"z":[1]}}"#, "not a string");
        bad(
            r#"{"id":0,"out":{"z":["nan"],"z":["nan"]}}"#,
            "duplicate key",
        );
        bad(r#"{"id":0,"out":{},"stratum":"s"}"#, "unknown field");
        bad(r#"{"out":{}}"#, "missing field");
        bad(r#"{"id":0,"out":{"z":["NaN"]}}"#, "bad hex float");
    }

    /// The error of ingesting `lines` (as `f.jsonl`) against a corpus file of `records` records.
    fn ingest_error(tag: &str, records: usize, lines: &[&str]) -> Result<String, String> {
        let scratch = Scratch::new(tag);
        std::fs::create_dir_all(&scratch.0).map_err(|e| e.to_string())?;
        let text: String = lines.iter().map(|l| format!("{l}\n")).collect();
        std::fs::write(scratch.0.join("f.jsonl"), text).map_err(|e| e.to_string())?;
        let e = FileSubject::ingest("x", &scratch.0, &[entry("f", records)]).err();
        Ok(e.unwrap_or_default())
    }

    #[test]
    fn a_skipped_record_a_wrong_id_and_a_stray_file_are_errors_not_scores() -> Result<(), String> {
        let one = r#"{"id":0,"out":{"z":["0x1.0000000000000p+0"]}}"#;
        let two = r#"{"id":1,"out":{"z":["0x1.0000000000000p+0"]}}"#;
        let e = ingest_error("short", 2, &[one])?;
        assert!(e.contains("1 answers, the corpus has 2 records"), "{e}");
        let e = ingest_error("id", 2, &[one, one])?;
        assert!(e.contains("f.jsonl:2: id 0 on line 2"), "{e}");
        assert!(ingest_error("ok", 2, &[one, two])?.is_empty());
        let scratch = Scratch::new("stray");
        std::fs::create_dir_all(&scratch.0).map_err(|e| e.to_string())?;
        std::fs::write(scratch.0.join("g.jsonl"), "").map_err(|e| e.to_string())?;
        let e = FileSubject::ingest("x", &scratch.0, &[entry("f", 0)]).err();
        assert!(e
            .unwrap_or_default()
            .contains("`g.jsonl` is not a corpus function"));
        Ok(())
    }

    #[test]
    fn the_pin_is_one_line() -> Result<(), String> {
        assert_eq!(version_line(b"tf_tree@abc\n")?, "tf_tree@abc");
        assert_eq!(version_line(b"tf_tree@abc")?, "tf_tree@abc");
        for bad in [&b""[..], b"\n", b"a\nb\n", b"a\r\n", b"\xff"] {
            assert!(version_line(bad).is_err(), "{bad:?}");
        }
        Ok(())
    }

    #[test]
    fn an_unknown_runner_and_an_unknown_function_stop_before_anything_runs() -> Result<(), String> {
        let (root, dir) = (super::super::root()?, corpus_dir()?);
        let scratch = Scratch::new("unknown");
        let go = |name, fn_id| {
            run(name, &root, &dir, &scratch.0, fn_id)
                .err()
                .unwrap_or_default()
        };
        assert!(go("nope", None)
            .contains("no oracle runner `nope`; have [\"tf_tree_math\", \"sophus_rs\"]"));
        assert!(go("tf_tree_math", Some("nope")).contains("no corpus file for `--fn nope`"));
        assert!(!scratch.0.exists());
        Ok(())
    }

    /// A misspelled id would otherwise never be looked for: `execute` checks only the owed ids
    /// that the corpus has.
    #[test]
    fn every_owed_function_is_a_corpus_function_of_a_runner_that_exists() -> Result<(), String> {
        let (root, entries) = (super::super::root()?, corpus::manifest(&corpus_dir()?)?);
        for runner in RUNNERS {
            assert!(
                root.join(runner.dir).join("Cargo.toml").is_file(),
                "{}",
                runner.name
            );
            for id in runner.answers {
                assert!(
                    entries.iter().any(|e| e.fn_id == *id),
                    "{}: {id}",
                    runner.name
                );
            }
        }
        Ok(())
    }

    #[test]
    fn a_runner_that_writes_more_than_the_corpus_has_is_an_error() -> Result<(), String> {
        let one = r#"{"id":0,"out":{"z":["0x1.0000000000000p+0"]}}"#;
        let two = r#"{"id":1,"out":{"z":["0x1.0000000000000p+0"]}}"#;
        let e = ingest_error("long", 1, &[one, two])?;
        assert!(e.contains("2 answers, the corpus has 1 records"), "{e}");
        Ok(())
    }

    #[test]
    fn a_function_the_runner_does_not_answer_is_not_supported_and_binary64_only(
    ) -> Result<(), String> {
        let scratch = Scratch::new("partial");
        std::fs::create_dir_all(&scratch.0).map_err(|e| e.to_string())?;
        let line = r#"{"id":0,"out":{"z":["0x1.0000000000000p+0"]}}"#;
        std::fs::write(scratch.0.join("f.jsonl"), format!("{line}\n"))
            .map_err(|e| e.to_string())?;
        let entries = [entry("f", 1), entry("g", 1)];
        let file = FileSubject::ingest("x", &scratch.0, &entries)?;
        assert!(file.supports("f") && !file.supports("g"));
        let record = testkit::record(&[], &[])?;
        assert_eq!(file.eval("f", &record, Precision::F64)["z"], [1.0]);
        assert!(file.eval("f", &record, Precision::F32).is_empty());
        Ok(())
    }

    #[test]
    fn the_runner_is_a_locked_release_run_of_its_own_manifest() {
        let command = runner_command(Path::new("/r"), &RUNNERS[0], true);
        let args: Vec<_> = command.get_args().map(|a| a.to_string_lossy()).collect();
        assert_eq!(
            args,
            [
                "run",
                "--release",
                "--locked",
                "--quiet",
                "--manifest-path",
                "/r/runners/tf_tree_math/Cargo.toml",
                "--"
            ]
        );
        assert_eq!(command.get_current_dir(), Some(Path::new("/r")));
    }

    /// The runner `stub`, answering `so2_exp`; `script` is the `sh` that stands in for it.
    const STUB: Runner = Runner {
        name: "stub",
        dir: "",
        answers: &["so2_exp"],
    };

    /// Answers `--version` with `stub@1` and `--out DIR …` with a copy of the files in `$ANSWERS`.
    const ANSWERING: &str = "case \"$1\" in
  --version) echo stub@1 ;;
  --out) for f in \"$ANSWERS\"/*.jsonl; do [ -e \"$f\" ] && cp \"$f\" \"$2\"; done ;;
esac
exit 0
";

    /// A scratch directory with the answers of `so2_exp` the harness's own next-up subject gives.
    struct Stub {
        scratch: Scratch,
        script: PathBuf,
    }

    impl Stub {
        fn new(tag: &str, script: &str, answers: bool) -> Result<Self, String> {
            let scratch = Scratch::new(tag);
            let dir = scratch.0.join("answers");
            std::fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
            if answers {
                let corpus = corpus_dir()?;
                let entries: Vec<Entry> = corpus::manifest(&corpus)?
                    .into_iter()
                    .filter(|e| e.fn_id == "so2_exp")
                    .collect();
                let records = corpus::read(&corpus, &entries[0])?;
                let writer = Perfect::next_up();
                let out: Vec<Output> = records
                    .iter()
                    .map(|r| writer.eval("so2_exp", r, Precision::F64))
                    .collect();
                write_answers(&dir, "so2_exp", &out)?;
            }
            let path = scratch.0.join("stub.sh");
            std::fs::write(&path, script).map_err(|e| e.to_string())?;
            Ok(Self {
                scratch,
                script: path,
            })
        }

        fn execute(&self, fn_id: Option<&str>) -> Result<Registered, String> {
            let launch = |_| {
                let mut command = Command::new("sh");
                command.arg(&self.script);
                command.env("ANSWERS", self.scratch.0.join("answers"));
                command
            };
            let results = self.scratch.0.join("results");
            execute(&STUB, &launch, &corpus_dir()?, &results, fn_id)
        }

        fn out(&self, dir: &str) -> PathBuf {
            self.scratch.0.join("results/oracle").join(dir)
        }
    }

    #[test]
    fn the_runner_is_registered_under_the_pin_it_prints() -> Result<(), String> {
        let stub = Stub::new("registered", ANSWERING, true)?;
        let registered = stub.execute(None)?;
        assert_eq!(registered.version, "stub@1");
        assert!(!registered.planted);
        assert_eq!(registered.subject.name(), "stub");
        assert!(registered.subject.supports("so2_exp") && !registered.subject.supports("so2_log"));
        assert!(stub.out("stub").join("so2_exp.jsonl").exists());
        Ok(())
    }

    #[test]
    fn what_an_earlier_run_left_is_removed_and_a_fn_run_has_a_directory_of_its_own(
    ) -> Result<(), String> {
        let stub = Stub::new("stale", ANSWERING, true)?;
        let (full, only) = (stub.out("stub"), stub.out("stub--so2_exp"));
        std::fs::create_dir_all(&full).map_err(|e| e.to_string())?;
        std::fs::write(full.join("so2_log.jsonl"), "stale\n").map_err(|e| e.to_string())?;
        stub.execute(None)?;
        assert!(!full.join("so2_log.jsonl").exists() && full.join("so2_exp.jsonl").exists());
        std::fs::create_dir_all(&only).map_err(|e| e.to_string())?;
        std::fs::write(only.join("stale.jsonl"), "").map_err(|e| e.to_string())?;
        stub.execute(Some("so2_exp"))?;
        assert!(!only.join("stale.jsonl").exists() && only.join("so2_exp.jsonl").exists());
        assert!(
            full.join("so2_exp.jsonl").exists(),
            "a full run's answers stay"
        );
        Ok(())
    }

    #[test]
    fn a_runner_that_fails_or_prints_no_pin_is_an_error() -> Result<(), String> {
        let error = |tag: &str, script: &str| -> Result<String, String> {
            Ok(Stub::new(tag, script, false)?
                .execute(None)
                .err()
                .unwrap_or_default())
        };
        let e = error("no-version", "exit 3\n")?;
        assert!(e.contains("`stub --version` failed"), "{e}");
        let e = error(
            "no-out",
            "case \"$1\" in --version) echo stub@1 ;; *) exit 1 ;; esac\n",
        )?;
        assert!(e.contains("runner `stub` failed"), "{e}");
        let e = error("two-lines", "echo a; echo b\n")?;
        assert!(e.contains("must print one line"), "{e}");
        Ok(())
    }

    #[test]
    fn a_runner_that_stops_answering_an_id_is_an_error_when_that_id_is_selected(
    ) -> Result<(), String> {
        let stub = Stub::new("silent", ANSWERING, false)?;
        for fn_id in [None, Some("so2_exp")] {
            let e = stub.execute(fn_id).err().unwrap_or_default();
            assert!(
                e.contains("runner `stub` wrote no answers for `so2_exp`"),
                "{e}"
            );
        }
        // Not selected, not owed: the run is empty, and the harness says nothing was scored.
        let registered = stub.execute(Some("so2_log"))?;
        assert!(!registered.subject.supports("so2_exp"));
        Ok(())
    }
}
