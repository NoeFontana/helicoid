//! `cargo xtask conformance --oracle NAME` (`docs/PHASE1.md` §7): runs the oracle runner `NAME`
//! over the corpus and scores what it wrote with the exact metric of an in-process subject, as a
//! [`FileSubject`]. A runner is a workspace-excluded crate under `runners/` with its own lockfile,
//! run by `cargo run --release --locked`, or an image built from `docker/oracles/<NAME>/Dockerfile`
//! and run offline ([`Launch`], `0072`). The file protocol, the smallest reading of §7
//! (0014 (draft) question 26):
//!
//! - `<runner> --version` prints the pin on one line: the `subject_version` of every row, at
//!   either precision.
//! - `<runner> --out DIR FILE.jsonl…` writes `DIR/<fn>.jsonl` for each corpus file whose function
//!   it answers, and nothing for the others, in binary64.
//! - `<runner> --precision f32 --out DIR FILE.jsonl…` does the same in binary32 for the `@f32`
//!   records alone (`0077`): it narrows each input losslessly, computes in binary32 and writes
//!   binary32 values. A runner with no [`Runner::answers_f32`] is never asked, and refuses the
//!   flag as any unknown argument.
//! - A line is `{"id":N,"out":{"<field>":["<hex float>",…]}}`: `N` the corpus line's `id`, a field
//!   flat in the reference's order (a matrix column-major, no `shape`), a value a hex float or
//!   `nan`, `inf`, `-inf`; a scalar field may be the bare string.
//! - An answer file holds every record of its corpus file at its precision, ids in order: every
//!   record in binary64, the `@f32` ones in binary32. A runner that skipped a
//!   record, answered a function the corpus does not have or wrote a malformed line is an error,
//!   not a score; so is a function in [`Runner::answers`] with no file, since a missing file
//!   otherwise reads as "not answered". A wrong answer is a score, and a non-finite one a recorded
//!   row (`nonfinite`), not a failed run: oracles may be wrong.
//!
//! `DIR` is `conformance/results/oracle/<NAME>` (`<NAME>--<fn>` for a `--fn` run, then `--f32` at
//! binary32, as its CSV), emptied first and left for inspection.

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

/// An oracle runner; `name` is the subject's name in the result rows.
struct Runner {
    name: &'static str,
    launch: Launch,
    /// The function ids it must answer: the list its own `--out` advertises, held again here so
    /// that an id it stops answering fails the run instead of dropping out of the results.
    answers: &'static [&'static str],
    /// Those it must answer in binary32 too, each with `@f32` strata (`0077`). Empty: it answers
    /// binary64 only, and a run at `f32` launches nothing.
    answers_f32: &'static [&'static str],
    /// Where its transcendentals come from: see [`Backend`].
    backend: Backend,
    /// The library has no binary32 kernel, so it answers the binary64 strata of an id that also
    /// has `@f32` ones; every other runner answers every such id at `f32` too.
    binary64_only: bool,
}

/// Where a runner's transcendentals come from, which decides what a domination failure can mean.
///
/// A runner on the **host's `std`** gets glibc on Linux, whose `atan2` is nearly correctly rounded
/// where `libm`'s is 1.96 u, so a stratum such a runner wins may cost nothing but D16
/// (`error-analysis.md` EA.13(d)); the one-variable twin says whether it does (`0037`, draft).
/// Every runner is on it: the one on the `libm` crate, `tf_tree_math`, is a consumer and not an
/// oracle ([`0107`](../../../docs/decisions/0107-the-consumer-is-not-an-oracle.md)).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Backend {
    /// The host's `std`, so its system library.
    HostStd,
}

/// How a runner is built and started; either way the protocol's arguments follow.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Launch {
    /// A workspace-excluded crate in this directory: `cargo run --release --locked`.
    Crate(&'static str),
    /// The image of this directory's `Dockerfile`, run with no network and a read-only root, the
    /// corpus mounted read-only and the answers' directory writable, each at its host path so the
    /// protocol's paths need no translation, as the owner of that directory (`0072`).
    Container(&'static str),
}

#[cfg(test)]
impl Launch {
    /// The file whose presence says the runner exists.
    fn manifest(self) -> String {
        match self {
            Launch::Crate(dir) => format!("{dir}/Cargo.toml"),
            Launch::Container(dir) => format!("{dir}/Dockerfile"),
        }
    }
}

const RUNNERS: &[Runner] = &[
    Runner {
        name: "sophus_rs",
        launch: Launch::Crate("runners/sophus_rs"),
        // sophus-rs is Rust `std`, so the host's library.
        backend: Backend::HostStd,
        binary64_only: true,
        // `sophus_autodiff` implements its scalar trait for `f64` alone.
        answers_f32: &[],
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
    Runner {
        name: "nalgebra",
        launch: Launch::Crate("runners/nalgebra"),
        // nalgebra is Rust `std`, so the host's library (`0056` decision 4).
        backend: Backend::HostStd,
        binary64_only: false,
        // Generic over `RealField`, and every id it answers has `@f32` strata.
        answers_f32: &[
            "eig3",
            "eig_sym_n4",
            "eig_sym_n9",
            "eig_sym_n16",
            "svd3",
            "nearest_rotation",
            "chol_n3",
            "chol_n6",
            "chol_solve_n3",
            "chol_solve_n6",
            "solve_cubic",
            "quat_renormalize",
            "mat2_inverse_adj",
        ],
        answers: &[
            "eig3",
            "eig_sym_n4",
            "eig_sym_n9",
            "eig_sym_n16",
            "svd3",
            "nearest_rotation",
            "chol_n3",
            "chol_n6",
            "chol_solve_n3",
            "chol_solve_n6",
            "solve_cubic",
            "quat_renormalize",
            "mat2_inverse_adj",
        ],
    },
    Runner {
        name: "sophus",
        launch: Launch::Container("docker/oracles/sophus"),
        // C++ on the image's glibc (`0072`).
        backend: Backend::HostStd,
        binary64_only: false,
        // `Sophus::Sim3<float>`, `Sophus::SO3<float>` of a matrix (`0085`) and its `fitToSO3`
        // (`0102`), Eigen's `SelfAdjointEigenSolver<Matrix3f>` (`0083`) and
        // `JacobiSVD<Matrix3f>` (`0101`); no other id it answers has `@f32` strata.
        answers_f32: &[
            "so3_from_matrix",
            "sim3_exp",
            "sim3_log",
            "sim3_ad",
            "sim3_act",
            "eig3",
            "eig_sym_n4",
            "eig_sym_n9",
            "eig_sym_n16",
            "svd3",
            "nearest_rotation",
        ],
        answers: &[
            "so2_exp",
            "so2_log",
            "se2_exp",
            "se2_log",
            "se2_ad",
            "so3_exp",
            "so3_log",
            "so3_act",
            "so3_from_matrix",
            "so3_jl",
            "so3_jr",
            "so3_jl_inv",
            "so3_jr_inv",
            "so3_geodesic",
            "se3_geodesic",
            "sen3_exp_n1",
            "sen3_log_n1",
            "sen3_ad_n1",
            "sen3_jl_n1",
            "sen3_jr_n1",
            "sen3_jl_inv_n1",
            "sen3_jr_inv_n1",
            "sim3_exp",
            "sim3_log",
            "sim3_ad",
            "sim3_act",
            // The image's Eigen, not Sophus (`0083`; `0101`; `0111`).
            "eig3",
            "eig_sym_n4",
            "eig_sym_n9",
            "eig_sym_n16",
            "svd3",
            // Sophus's `fitToSO3` (`0102`).
            "nearest_rotation",
        ],
    },
    Runner {
        name: "gtsam",
        launch: Launch::Container("docker/oracles/gtsam"),
        // C++ on the image's glibc (`0074`).
        backend: Backend::HostStd,
        binary64_only: true,
        // GTSAM is `double` throughout.
        answers_f32: &[],
        answers: &[
            "so2_exp",
            "so2_log",
            "se2_exp",
            "se2_log",
            "se2_ad",
            "se2_jl",
            "se2_jr",
            "se2_jl_inv",
            "se2_jr_inv",
            "so3_exp",
            "so3_log",
            "so3_act",
            "so3_inv_act",
            "se3_inv_act",
            "so3_from_matrix",
            "so3_jl",
            "so3_jr",
            "so3_jl_inv",
            "so3_jr_inv",
            "so3_gamma2",
            "so3_geodesic",
            "se3_geodesic",
            "se3_screw_retract",
            "se3_screw_local",
            "sen3_exp_n1",
            "sen3_log_n1",
            "sen3_ad_n1",
            "sen3_jl_n1",
            "sen3_jr_n1",
            "sen3_jl_inv_n1",
            "sen3_jr_inv_n1",
            "sen3_exp_n2",
            "sen3_log_n2",
            "sen3_ad_n2",
            "sen3_jl_n2",
            "sen3_jr_n2",
            "sen3_jl_inv_n2",
            "sen3_jr_inv_n2",
            "sen3_exp_n3",
            "sen3_log_n3",
            "sen3_ad_n3",
            "sen3_jl_n3",
            "sen3_jr_n3",
            "sen3_jl_inv_n3",
            "sen3_jr_inv_n3",
            "sim3_exp",
            "sim3_log",
            "sim3_ad",
            "sim3_act",
        ],
    },
    Runner {
        name: "manif",
        launch: Launch::Container("docker/oracles/manif"),
        // C++ on the image's glibc (`0075`).
        backend: Backend::HostStd,
        binary64_only: false,
        // `se3_inv_act` (`0103`) is the one id it answers with `@f32` strata, in
        // `manif::SE3<float>`.
        answers_f32: &["se3_inv_act"],
        answers: &[
            "so2_exp",
            "so2_log",
            "se2_exp",
            "se2_log",
            "se2_ad",
            "se2_jl",
            "se2_jr",
            "se2_jl_inv",
            "se2_jr_inv",
            "so3_exp",
            "so3_log",
            "so3_act",
            "so3_inv_act",
            "se3_inv_act",
            "so3_jl",
            "so3_jr",
            "so3_jl_inv",
            "so3_jr_inv",
            "so3_geodesic",
            "se3_geodesic",
            "se3_screw_retract",
            "se3_screw_local",
            "se3_world_retract",
            "se3_world_local",
            "sen3_exp_n1",
            "sen3_log_n1",
            "sen3_ad_n1",
            "sen3_jl_n1",
            "sen3_jr_n1",
            "sen3_jl_inv_n1",
            "sen3_jr_inv_n1",
            "sen3_exp_n2",
            "sen3_log_n2",
            "sen3_ad_n2",
            "sen3_jl_n2",
            "sen3_jr_n2",
            "sen3_jl_inv_n2",
            "sen3_jr_inv_n2",
        ],
    },
];

pub(super) fn names() -> Vec<&'static str> {
    RUNNERS.iter().map(|r| r.name).collect()
}

/// Where `name`'s transcendentals come from, or `None` if it is not a registered runner.
pub(crate) fn backend(name: &str) -> Option<Backend> {
    RUNNERS.iter().find(|r| r.name == name).map(|r| r.backend)
}

/// Whether the runner `name`'s library has no binary32 kernel, as against one that has no id with
/// `@f32` strata to answer: what a run at `f32` that launches nothing says.
pub(crate) fn binary64_only(name: &str) -> bool {
    RUNNERS.iter().any(|r| r.name == name && r.binary64_only)
}

/// Whether the runner `name` answers in binary32, and so owes `f32` result rows.
pub(crate) fn answers_f32(name: &str) -> bool {
    RUNNERS
        .iter()
        .any(|r| r.name == name && !r.answers_f32.is_empty())
}

/// A runner's answers at one precision, read back from its files.
pub(crate) struct FileSubject {
    name: String,
    precision: Precision,
    /// Per function id, the answer to each record of that precision, by record `id`.
    answers: BTreeMap<String, BTreeMap<u64, Output>>,
}

impl Subject for FileSubject {
    fn name(&self) -> &str {
        &self.name
    }

    fn supports(&self, fn_id: &str) -> bool {
        self.answers.contains_key(fn_id)
    }

    /// The runner's answer, or none, which the metric reports as a missing field: ingestion has
    /// already checked that no record is missing. At another precision than its run's there is
    /// none.
    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        if precision != self.precision {
            return Output::new();
        }
        let answer = self.answers.get(fn_id).and_then(|a| a.get(&record.id));
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

/// The ids of an answer file of `entry` at `precision`, in order: every record's in binary64, the
/// `@f32` records' in binary32, which only the corpus in `corpus_dir` can say.
fn answered_ids(
    corpus_dir: &Path,
    entry: &Entry,
    precision: Precision,
) -> Result<Vec<u64>, String> {
    match precision {
        Precision::F64 => Ok((0..entry.records as u64).collect()),
        Precision::F32 => Ok(corpus::read(corpus_dir, entry)?
            .iter()
            .filter(|r| r.is_f32_stratum())
            .map(|r| r.id)
            .collect()),
    }
}

impl FileSubject {
    /// The answers in `dir` at `precision`, one `<fn>.jsonl` per function of `entries` the runner
    /// answers.
    pub(crate) fn ingest(
        name: &str,
        dir: &Path,
        entries: &[Entry],
        precision: Precision,
        corpus_dir: &Path,
    ) -> Result<Self, String> {
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
            let ids = answered_ids(corpus_dir, entry, precision)?;
            let (mut file, mut lines) = (BTreeMap::new(), 0);
            for (i, line) in text.lines().enumerate() {
                let at = |e: String| format!("{}:{}: {e}", path.display(), i + 1);
                let (id, out) = parse_answer(line).map_err(at)?;
                if ids.get(i).is_some_and(|&want| want != id) {
                    return Err(at(format!("id {id} on line {}", i + 1)));
                }
                file.insert(id, out);
                lines += 1;
            }
            if lines != ids.len() {
                let of = match precision {
                    Precision::F64 => "",
                    Precision::F32 => "@f32 ",
                };
                return Err(format!(
                    "{}: {lines} answers, the corpus has {} {of}records",
                    path.display(),
                    ids.len()
                ));
            }
            answers.insert(entry.fn_id.clone(), file);
        }
        Ok(Self {
            name: name.to_string(),
            precision,
            answers,
        })
    }
}

/// `cargo run --release --locked` on the crate in `dir`, arguments after `--` still to add.
fn crate_command(root: &Path, dir: &str, quiet: bool) -> Command {
    let cargo = std::env::var_os("CARGO").unwrap_or_else(|| OsString::from("cargo"));
    let mut command = Command::new(cargo);
    command
        .current_dir(root)
        .args(["run", "--release", "--locked"])
        .args(quiet.then_some("--quiet"))
        .arg("--manifest-path")
        .arg(root.join(dir).join("Cargo.toml"))
        .arg("--");
    command
}

/// The tag a container runner's image is built under, for a reader of `docker images`. It is run
/// by the id its build printed, which a later build cannot move.
fn image_tag(name: &str) -> String {
    format!("helicoid-oracle-{name}")
}

/// `docker build --quiet` of `dir`'s `Dockerfile`, which prints the image id alone.
fn build_command(root: &Path, name: &str, dir: &str) -> Command {
    let mut command = Command::new("docker");
    command
        .args(["build", "--quiet", "--tag", &image_tag(name)])
        .arg(root.join(dir));
    command
}

/// The image id `docker build --quiet` printed: its last line, `sha256:` and hex digits.
fn image_id(stdout: &[u8]) -> Result<String, String> {
    let text = std::str::from_utf8(stdout).map_err(|e| format!("image id is not UTF-8: {e}"))?;
    let id = text.lines().map(str::trim).rfind(|l| !l.is_empty());
    let id = id.unwrap_or_default();
    let hex = id.strip_prefix("sha256:").unwrap_or_default();
    match !hex.is_empty() && hex.bytes().all(|b| b.is_ascii_hexdigit()) {
        true => Ok(id.to_string()),
        false => Err(format!("`docker build` printed no image id: {text:?}")),
    }
}

/// A bind mount of `path` at the same path, which `--mount`'s comma-separated syntax can spell
/// only for an absolute path without a comma.
fn bind(path: &Path, readonly: bool) -> Result<String, String> {
    let text = path
        .to_str()
        .filter(|t| path.is_absolute() && !t.contains(','))
        .ok_or_else(|| format!("{}: not an absolute path without a comma", path.display()))?;
    let ro = if readonly { ",readonly" } else { "" };
    Ok(format!("--mount=type=bind,src={text},dst={text}{ro}"))
}

/// `docker run` of `image` with `mounts` (from [`bind`]), as `user` (`uid:gid`) when given; the
/// protocol's arguments still to add.
fn container_command(image: &str, mounts: &[String], user: Option<(u32, u32)>) -> Command {
    let mut command = Command::new("docker");
    command
        .args(["run", "--rm", "--network=none", "--read-only"])
        .args(user.map(|(uid, gid)| format!("--user={uid}:{gid}")))
        .args(mounts)
        .arg(image);
    command
}

/// The owner of `dir`, whom a container must run as for the files it writes there to be theirs.
#[cfg(unix)]
fn owner(dir: &Path) -> Result<Option<(u32, u32)>, String> {
    use std::os::unix::fs::MetadataExt;
    let meta = std::fs::metadata(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    Ok(Some((meta.uid(), meta.gid())))
}

#[cfg(not(unix))]
fn owner(_dir: &Path) -> Result<Option<(u32, u32)>, String> {
    Ok(None)
}

/// Builds a container runner's image, which docker caches after the first build, and returns its
/// id.
fn build_image(root: &Path, name: &str, dir: &str) -> Result<String, String> {
    eprintln!("conformance: building {dir}/Dockerfile");
    let out = build_command(root, name, dir)
        .stderr(Stdio::inherit())
        .output()
        .map_err(|e| format!("`docker build` of {dir}: {e}"))?;
    succeeded(&format!("`docker build` of {dir}"), out.status)?;
    image_id(&out.stdout)
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

/// Builds and runs the runner `name` over the corpus (only `fn_id`'s file when given) at
/// `precision` and reads its answers: the subject, its `subject_version` the runner's pin. `None`
/// at `f32` for a runner that answers binary64 only, which is not launched.
pub(crate) fn run(
    name: &str,
    root: &Path,
    corpus_dir: &Path,
    results: &Path,
    fn_id: Option<&str>,
    precision: Precision,
) -> Result<Option<Registered>, String> {
    let runner = RUNNERS.iter().find(|r| r.name == name).ok_or_else(|| {
        let known: Vec<&str> = RUNNERS.iter().map(|r| r.name).collect();
        format!("no oracle runner `{name}`; have {known:?}")
    })?;
    if precision == Precision::F32 && runner.answers_f32.is_empty() {
        return Ok(None);
    }
    let at = Run {
        corpus_dir,
        results,
        fn_id,
        precision,
    };
    let registered = match runner.launch {
        Launch::Crate(dir) => {
            let launch = |quiet| crate_command(root, dir, quiet);
            execute(runner, &launch, &at)
        }
        Launch::Container(dir) => {
            // Every answer directory is under `oracle/`, so it is the one writable mount. It must
            // exist first: docker would create a missing one as root.
            let out = results.join("oracle");
            std::fs::create_dir_all(&out).map_err(|e| format!("{}: {e}", out.display()))?;
            let mounts = [bind(corpus_dir, true)?, bind(&out, false)?];
            let user = owner(&out)?;
            let image = build_image(root, runner.name, dir)?;
            let launch = |_| container_command(&image, &mounts, user);
            execute(runner, &launch, &at)
        }
    };
    registered.map(Some)
}

/// What one run of a runner reads, writes and asks.
struct Run<'a> {
    corpus_dir: &'a Path,
    results: &'a Path,
    fn_id: Option<&'a str>,
    precision: Precision,
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
    at: &Run,
) -> Result<Registered, String> {
    let Run {
        corpus_dir,
        results,
        fn_id,
        precision,
    } = *at;
    let mut entries = corpus::manifest(corpus_dir)?;
    if let Some(id) = fn_id {
        entries.retain(|e| e.fn_id == id);
        if entries.is_empty() {
            return Err(format!("no corpus file for `--fn {id}`"));
        }
    }
    let (owed, precision_args) = match precision {
        Precision::F64 => (runner.answers, &[][..]),
        Precision::F32 => (runner.answers_f32, &["--precision", "f32"][..]),
    };
    if precision == Precision::F32 {
        // A file with no `@f32` record has nothing to answer in binary32.
        let mut with = Vec::new();
        for e in entries {
            if corpus::mentions_f32(corpus_dir, &e)? {
                with.push(e);
            }
        }
        entries = with;
        if let (Some(id), true) = (fn_id, entries.is_empty()) {
            return Err(format!(
                "`{id}` has no `@f32` stratum to answer in binary32"
            ));
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
    let mut dir = match fn_id {
        Some(id) => format!("{}--{id}", runner.name),
        None => runner.name.to_string(),
    };
    if precision == Precision::F32 {
        dir += "--f32";
    }
    let out_dir = results.join("oracle").join(dir);
    fresh_dir(&out_dir)?;
    let files = entries
        .iter()
        .map(|e| corpus_dir.join(format!("{}.jsonl", e.fn_id)));
    let status = launch(true)
        .args(precision_args)
        .arg("--out")
        .arg(&out_dir)
        .args(files)
        .status()
        .map_err(|e| format!("{}: {e}", runner.name))?;
    succeeded(&format!("runner `{}`", runner.name), status)?;
    let subject = FileSubject::ingest(runner.name, &out_dir, &entries, precision, corpus_dir)?;
    let unanswered = owed
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
        // One pin at both precisions: the row's `precision` column names the instantiation.
        version_f32: version,
        // The run's answers are of its precision alone, so a score at the other is an error rather
        // than another kernel's answer under this name (`0016` item 2).
        no_f32: match precision {
            Precision::F64 => Some(format!("the `{}` runner was run in binary64", runner.name)),
            Precision::F32 => None,
        },
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

    fn write_answers(dir: &Path, fn_id: &str, answers: &[(u64, Output)]) -> Result<(), String> {
        std::fs::create_dir_all(dir).map_err(|e| e.to_string())?;
        let mut text = String::new();
        for (id, out) in answers {
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

    /// An answer file's lines: each record `id` with its answer.
    type Answers = Vec<(u64, Output)>;

    /// `fn_id`'s corpus entry and the harness's next-up answers to its records at `precision`: every
    /// record in binary64, the `@f32` ones in binary32, as a runner writes them.
    fn next_up_answers(fn_id: &str, precision: Precision) -> Result<(Vec<Entry>, Answers), String> {
        let dir = corpus_dir()?;
        let entries: Vec<Entry> = corpus::manifest(&dir)?
            .into_iter()
            .filter(|e| e.fn_id == fn_id)
            .collect();
        let records = corpus::read(&dir, &entries[0])?;
        let writer = Perfect::next_up();
        let answers = records
            .iter()
            .filter(|r| precision == Precision::F64 || r.is_f32_stratum())
            .map(|r| (r.id, writer.eval(fn_id, r, precision)))
            .collect();
        Ok((entries, answers))
    }

    /// At each precision, a file of the next-up subject's answers scores as that subject does;
    /// at `f32` the file holds the `@f32` records alone and its values are binary32.
    #[test]
    fn an_ingested_answer_scores_exactly_like_the_subject_that_wrote_it() -> Result<(), String> {
        let dir = corpus_dir()?;
        for (fn_id, precision) in [
            ("so2_exp", Precision::F64),
            ("quat_renormalize", Precision::F64),
            ("quat_renormalize", Precision::F32),
        ] {
            let (entries, answers) = next_up_answers(fn_id, precision)?;
            let scratch = Scratch::new("ingest");
            write_answers(&scratch.0, fn_id, &answers)?;
            let file = FileSubject::ingest("next-up", &scratch.0, &entries, precision, &dir)?;
            assert!(file.supports(fn_id) && !file.supports("so2_log"));
            let score = |s: Box<dyn Subject>| {
                evaluate(&dir, &entries, &[Registered::new("1", s)], precision)
            };
            let from_file = score(Box::new(file))?;
            let direct = score(Box::new(Perfect::next_up()))?;
            let records = corpus::read(&dir, &entries[0])?;
            let strata: std::collections::BTreeSet<&str> = records
                .iter()
                .filter(|r| r.is_f32_stratum() == (precision == Precision::F32))
                .map(|r| r.stratum.as_str())
                .collect();
            assert_eq!(from_file[0].len(), strata.len(), "{fn_id}");
            assert!(from_file[0].iter().any(|r| r.max_u > 0.0));
            assert!(from_file[0].iter().all(|r| r.precision == precision));
            assert_eq!(from_file, direct, "{fn_id}");
        }
        Ok(())
    }

    /// A binary32 file is the `@f32` records' answers: one answering every record, or skipping
    /// the first `@f32` one, is an error naming what the corpus has.
    #[test]
    fn a_binary32_file_holds_the_f32_records_and_no_other() -> Result<(), String> {
        let dir = corpus_dir()?;
        let (entries, all) = next_up_answers("mat2_inverse_adj", Precision::F64)?;
        let (_, f32s) = next_up_answers("mat2_inverse_adj", Precision::F32)?;
        assert!(!f32s.is_empty() && f32s.len() < all.len());
        let error = |answers: &[(u64, Output)]| -> Result<String, String> {
            let scratch = Scratch::new("f32-ids");
            write_answers(&scratch.0, "mat2_inverse_adj", answers)?;
            let e = FileSubject::ingest("x", &scratch.0, &entries, Precision::F32, &dir).err();
            Ok(e.unwrap_or_default())
        };
        let e = error(&all)?;
        assert!(
            e.contains("mat2_inverse_adj.jsonl:1: id 0 on line 1"),
            "{e}"
        );
        let e = error(&f32s[1..])?;
        assert!(e.contains(&format!("id {} on line 1", f32s[1].0)), "{e}");
        let e = error(&f32s[..1])?;
        let want = format!("1 answers, the corpus has {} @f32 records", f32s.len());
        assert!(e.contains(&want), "{e}");
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

    /// Binary64 answers in `dir` to `entries`, which are not in the corpus: at binary64 the ids
    /// are the line numbers, so no corpus file is read.
    fn ingest64(dir: &Path, entries: &[Entry]) -> Result<FileSubject, String> {
        FileSubject::ingest("x", dir, entries, Precision::F64, Path::new("/no-corpus"))
    }

    /// The error of ingesting `lines` (as `f.jsonl`) against a corpus file of `records` records.
    fn ingest_error(tag: &str, records: usize, lines: &[&str]) -> Result<String, String> {
        let scratch = Scratch::new(tag);
        std::fs::create_dir_all(&scratch.0).map_err(|e| e.to_string())?;
        let text: String = lines.iter().map(|l| format!("{l}\n")).collect();
        std::fs::write(scratch.0.join("f.jsonl"), text).map_err(|e| e.to_string())?;
        let e = ingest64(&scratch.0, &[entry("f", records)]).err();
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
        let e = ingest64(&scratch.0, &[entry("f", 0)]).err();
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
        let go = |name, fn_id, precision| {
            run(name, &root, &dir, &scratch.0, fn_id, precision)
                .err()
                .unwrap_or_default()
        };
        assert!(go("nope", None, Precision::F64).contains(
            "no oracle runner `nope`; have [\"sophus_rs\", \"nalgebra\", \"sophus\", \"gtsam\", \"manif\"]"
        ));
        let e = go("sophus_rs", Some("nope"), Precision::F64);
        assert!(e.contains("no corpus file for `--fn nope`"));
        assert!(!scratch.0.exists());
        Ok(())
    }

    /// A runner with no binary32 ids is not launched at `f32` (no cargo, no docker), and the run is
    /// not a failure: an `f32` pass over every runner skips it (`0077`).
    #[test]
    fn a_binary64_runner_is_not_launched_at_f32() -> Result<(), String> {
        let (root, dir) = (super::super::root()?, corpus_dir()?);
        let scratch = Scratch::new("f64-only");
        for name in ["sophus_rs", "gtsam"] {
            let none = run(name, &root, &dir, &scratch.0, None, Precision::F32)?;
            assert!(none.is_none(), "{name}");
            assert!(!answers_f32(name), "{name}");
        }
        assert!(answers_f32("nalgebra") && answers_f32("sophus") && answers_f32("manif"));
        assert!(!answers_f32("nope"));
        assert!(!scratch.0.exists());
        Ok(())
    }

    /// A misspelled id would otherwise never be looked for: `execute` checks only the owed ids
    /// that the corpus has. A binary32 id is one the runner answers whose file has `@f32` strata,
    /// and every id it answers with such strata is one: answering it in binary64 only would leave
    /// those strata unpaired with no word said.
    #[test]
    fn every_owed_function_is_a_corpus_function_of_a_runner_that_exists() -> Result<(), String> {
        let dir = corpus_dir()?;
        let (root, entries) = (super::super::root()?, corpus::manifest(&dir)?);
        for runner in RUNNERS {
            assert!(
                root.join(runner.launch.manifest()).is_file(),
                "{}",
                runner.name
            );
            for id in runner.answers {
                let entry = entries.iter().find(|e| e.fn_id == *id);
                let entry = entry.ok_or_else(|| format!("{}: {id}", runner.name))?;
                let f32 = runner.answers_f32.contains(id);
                let has_f32 = corpus::mentions_f32(&dir, entry)?;
                assert_eq!(
                    f32,
                    has_f32 && !runner.binary64_only,
                    "{}: {id}",
                    runner.name
                );
            }
            for id in runner.answers_f32 {
                assert!(runner.answers.contains(id), "{}: {id}", runner.name);
            }
        }
        Ok(())
    }

    /// An image's build context is its own directory, so the container runners each hold a copy
    /// of the hex-float codec (`0075`); a fix to one copy that missed another fails here.
    #[test]
    fn the_container_runners_share_one_hex_float_codec() -> Result<(), String> {
        let root = super::super::root()?;
        let mut copies = RUNNERS.iter().filter_map(|r| match r.launch {
            Launch::Container(dir) => Some(root.join(dir).join("hexfloat.hpp")),
            Launch::Crate(_) => None,
        });
        let first = copies.next().ok_or("no container runner")?;
        let read = |p: &Path| std::fs::read(p).map_err(|e| format!("{}: {e}", p.display()));
        let codec = read(&first)?;
        for copy in copies {
            assert!(
                read(&copy)? == codec,
                "{} differs from {}",
                copy.display(),
                first.display()
            );
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

    /// A file subject answers at the precision of its run alone.
    #[test]
    fn a_function_the_runner_does_not_answer_is_not_supported_and_a_run_has_one_precision(
    ) -> Result<(), String> {
        let scratch = Scratch::new("partial");
        std::fs::create_dir_all(&scratch.0).map_err(|e| e.to_string())?;
        let line = r#"{"id":0,"out":{"z":["0x1.0000000000000p+0"]}}"#;
        std::fs::write(scratch.0.join("f.jsonl"), format!("{line}\n"))
            .map_err(|e| e.to_string())?;
        let entries = [entry("f", 1), entry("g", 1)];
        let file = ingest64(&scratch.0, &entries)?;
        assert!(file.supports("f") && !file.supports("g"));
        let record = testkit::record(&[], &[])?;
        assert_eq!(file.eval("f", &record, Precision::F64)["z"], [1.0]);
        assert!(file.eval("f", &record, Precision::F32).is_empty());
        Ok(())
    }

    #[test]
    fn the_runner_is_a_locked_release_run_of_its_own_manifest() {
        let command = crate_command(Path::new("/r"), "runners/sophus_rs", true);
        let args: Vec<_> = command.get_args().map(|a| a.to_string_lossy()).collect();
        assert_eq!(
            args,
            [
                "run",
                "--release",
                "--locked",
                "--quiet",
                "--manifest-path",
                "/r/runners/sophus_rs/Cargo.toml",
                "--"
            ]
        );
        assert_eq!(command.get_current_dir(), Some(Path::new("/r")));
    }

    fn args_of(command: &Command) -> Vec<String> {
        command
            .get_args()
            .map(|a| a.to_string_lossy().into_owned())
            .collect()
    }

    #[test]
    fn a_container_runner_is_built_from_its_dockerfile_and_run_offline_by_image_id(
    ) -> Result<(), String> {
        let build = build_command(Path::new("/r"), "sophus", "docker/oracles/sophus");
        assert_eq!(build.get_program(), "docker");
        assert_eq!(
            args_of(&build),
            [
                "build",
                "--quiet",
                "--tag",
                "helicoid-oracle-sophus",
                "/r/docker/oracles/sophus"
            ]
        );
        let mounts = [
            bind(Path::new("/r/corpus"), true)?,
            bind(Path::new("/r/results/oracle"), false)?,
        ];
        let run = container_command("sha256:ab12", &mounts, Some((1000, 100)));
        assert_eq!(run.get_program(), "docker");
        assert_eq!(
            args_of(&run),
            [
                "run",
                "--rm",
                "--network=none",
                "--read-only",
                "--user=1000:100",
                "--mount=type=bind,src=/r/corpus,dst=/r/corpus,readonly",
                "--mount=type=bind,src=/r/results/oracle,dst=/r/results/oracle",
                "sha256:ab12"
            ]
        );
        let anyone = container_command("sha256:ab12", &mounts, None);
        assert!(!args_of(&anyone).iter().any(|a| a.starts_with("--user")));
        Ok(())
    }

    #[test]
    fn a_mount_is_an_absolute_path_without_a_comma() {
        for bad in ["relative/dir", "/a,b"] {
            let e = bind(Path::new(bad), true).err().unwrap_or_default();
            assert!(
                e.contains("not an absolute path without a comma"),
                "{bad}: {e}"
            );
        }
    }

    #[test]
    fn the_image_id_is_the_last_line_docker_build_prints() -> Result<(), String> {
        let hex = "0123456789abcdef".repeat(4);
        let id = format!("sha256:{hex}");
        assert_eq!(image_id(format!("{id}\n").as_bytes())?, id);
        assert_eq!(image_id(format!("noise\n{id}\n\n").as_bytes())?, id);
        for bad in [
            &b""[..],
            b"\n",
            b"sha256:\n",
            b"sha256:xyz\n",
            b"Successfully built 1\n",
        ] {
            assert!(image_id(bad).is_err(), "{bad:?}");
        }
        Ok(())
    }

    /// The owner of a directory this process created is this process's user, whom the container
    /// runs as.
    #[cfg(unix)]
    #[test]
    fn the_container_runs_as_the_owner_of_the_answers_directory() -> Result<(), String> {
        use std::os::unix::fs::MetadataExt;
        let scratch = Scratch::new("owner");
        std::fs::create_dir_all(&scratch.0).map_err(|e| e.to_string())?;
        let meta = std::fs::metadata(&scratch.0).map_err(|e| e.to_string())?;
        assert_eq!(owner(&scratch.0)?, Some((meta.uid(), meta.gid())));
        assert!(owner(&scratch.0.join("missing")).is_err());
        Ok(())
    }

    /// Each registered runner's backend, and an unregistered name answering `None`.
    ///
    /// sophus-rs and nalgebra are Rust `std`, and the Sophus, GTSAM and manif images are C++ on
    /// their glibc. The domination split reads this, so a wrong entry would silently mislabel a
    /// failure (`0036`, draft).
    #[test]
    fn every_runner_declares_where_its_transcendentals_come_from() {
        assert_eq!(backend("sophus_rs"), Some(Backend::HostStd));
        assert_eq!(backend("nalgebra"), Some(Backend::HostStd));
        assert_eq!(backend("sophus"), Some(Backend::HostStd));
        assert_eq!(backend("gtsam"), Some(Backend::HostStd));
        assert_eq!(backend("manif"), Some(Backend::HostStd));
        assert_eq!(backend("not_a_runner"), None);
        assert_eq!(names().len(), RUNNERS.len());
    }

    /// The runner `stub`, answering `so2_exp` and, at both precisions, `quat_renormalize`;
    /// `script` is the `sh` that stands in for it.
    const STUB: Runner = Runner {
        name: "stub",
        launch: Launch::Crate(""),
        backend: Backend::HostStd,
        binary64_only: false,
        answers: &["so2_exp", "quat_renormalize"],
        answers_f32: &["quat_renormalize"],
    };

    /// Answers `--version` with `stub@1`, and `[--precision P] --out DIR FILE…` with a copy of each
    /// file's answers in `$ANSWERS/P` (`f64` without the flag); logs the arguments of an answering
    /// run to `$ANSWERS/args`.
    const ANSWERING: &str = "p=f64
if [ \"$1\" = --precision ]; then p=$2; shift 2; fi
case \"$1\" in
  --version) echo stub@1 ;;
  --out) echo \"$p $*\" >> \"$ANSWERS/args\"; d=$2; shift 2
    for f in \"$@\"; do a=\"$ANSWERS/$p/${f##*/}\"; [ -e \"$a\" ] && cp \"$a\" \"$d\"; done ;;
esac
exit 0
";

    /// A scratch directory with the answers the harness's own next-up subject gives, to `so2_exp`
    /// and `quat_renormalize` in binary64 and to `quat_renormalize`'s `@f32` records in binary32.
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
                for (fn_id, precision, sub) in [
                    ("so2_exp", Precision::F64, "f64"),
                    ("quat_renormalize", Precision::F64, "f64"),
                    ("quat_renormalize", Precision::F32, "f32"),
                ] {
                    let (_, out) = next_up_answers(fn_id, precision)?;
                    write_answers(&dir.join(sub), fn_id, &out)?;
                }
            }
            let path = scratch.0.join("stub.sh");
            std::fs::write(&path, script).map_err(|e| e.to_string())?;
            Ok(Self {
                scratch,
                script: path,
            })
        }

        fn execute(&self, fn_id: Option<&str>) -> Result<Registered, String> {
            self.execute_at(fn_id, Precision::F64)
        }

        fn execute_at(
            &self,
            fn_id: Option<&str>,
            precision: Precision,
        ) -> Result<Registered, String> {
            let launch = |_| {
                let mut command = Command::new("sh");
                command.arg(&self.script);
                command.env("ANSWERS", self.scratch.0.join("answers"));
                command
            };
            let (corpus, results) = (corpus_dir()?, self.scratch.0.join("results"));
            let at = Run {
                corpus_dir: &corpus,
                results: &results,
                fn_id,
                precision,
            };
            execute(&STUB, &launch, &at)
        }

        /// The arguments of each answering run, after its precision.
        fn args(&self) -> Result<String, String> {
            let path = self.scratch.0.join("answers/args");
            std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))
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
        let e = stub
            .execute_at(None, Precision::F32)
            .err()
            .unwrap_or_default();
        assert!(
            e.contains("runner `stub` wrote no answers for `quat_renormalize`"),
            "{e}"
        );
        Ok(())
    }

    /// At `f32` the runner is asked `--precision f32 --out DIR` with the files that have `@f32`
    /// strata alone, into a directory of its own, and its answers score at `f32` only and as the
    /// subject that wrote them.
    #[test]
    fn a_binary32_run_asks_for_the_f32_files_and_scores_at_f32() -> Result<(), String> {
        let stub = Stub::new("f32", ANSWERING, true)?;
        let registered = stub.execute_at(None, Precision::F32)?;
        assert!(registered.no_f32.is_none() && stub.execute(None)?.no_f32.is_some());
        assert_eq!(registered.version_at(Precision::F32), "stub@1");
        let args = stub.args()?;
        let f32_run = args
            .lines()
            .find(|l| l.starts_with("f32 "))
            .unwrap_or_default();
        let out = stub.out("stub--f32");
        assert!(f32_run.starts_with(&format!("f32 --out {} ", out.display())));
        let corpus = corpus_dir()?;
        let files: Vec<&str> = f32_run.split(' ').skip(3).collect();
        assert!(files.contains(
            &corpus
                .join("quat_renormalize.jsonl")
                .to_string_lossy()
                .as_ref()
        ));
        for file in &files {
            assert!(std::fs::read_to_string(file)
                .map_err(|e| e.to_string())?
                .contains("@f32\""));
        }
        assert!(
            !files.iter().any(|f| f.ends_with("/so2_exp.jsonl")),
            "{f32_run}"
        );
        assert!(out.join("quat_renormalize.jsonl").exists() && stub.out("stub").exists());
        let subject = registered.subject;
        assert!(subject.supports("quat_renormalize") && !subject.supports("so2_exp"));
        let entries: Vec<Entry> = corpus::manifest(&corpus)?
            .into_iter()
            .filter(|e| e.fn_id == "quat_renormalize")
            .collect();
        let scored = evaluate(
            &corpus,
            &entries,
            &[Registered::new("1", subject)],
            Precision::F32,
        )?;
        let direct = evaluate(
            &corpus,
            &entries,
            &[Registered::new("1", Box::new(Perfect::next_up()))],
            Precision::F32,
        )?;
        assert!(!scored[0].is_empty());
        assert_eq!(
            scored[0].iter().map(|r| r.max_u).collect::<Vec<_>>(),
            direct[0].iter().map(|r| r.max_u).collect::<Vec<_>>()
        );
        // A selected id with no `@f32` stratum has nothing to answer in binary32.
        let e = stub
            .execute_at(Some("so2_exp"), Precision::F32)
            .err()
            .unwrap_or_default();
        assert!(e.contains("`so2_exp` has no `@f32` stratum"), "{e}");
        let only = stub.execute_at(Some("quat_renormalize"), Precision::F32)?;
        assert!(only.subject.supports("quat_renormalize"));
        assert!(stub.out("stub--quat_renormalize--f32").exists());
        Ok(())
    }
}
