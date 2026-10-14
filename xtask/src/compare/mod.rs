//! `cargo xtask compare --a NAME --b NAME [--fn ID] [--precision f64|f32]`: a paired comparison of
//! two in-process subjects over the committed corpus (`0038` decision 4). It reads each corpus file
//! once, scores both subjects on every record with the conformance metric (`NUMERICS.md` §11), and
//! prints, per `(fn, stratum, precision)`, then per id with its strata pooled, then over every id
//! pooled, the paired direction, an exact sign test, the geometric mean of the per-record ratio
//! `b/a` and its extremes in both directions, and the ratio of the two per-stratum maxima
//! (decisions 5 and 6; [`stats`]).
//!
//! A reporting task and not a bar: it writes no result CSV, touches no baseline, and is read by no
//! bar, so it passes whatever it finds and fails only when it cannot run. A non-finite output is
//! counted and compared with nothing, where `conformance` would fail the run. Both subjects may be
//! planted; neither may be an out-of-process oracle. The readings `0038` leaves open are `0078`'s.

mod stats;

use std::collections::BTreeMap;
use std::fmt::Write;
use std::path::Path;

use helicoid_linalg::Precision;

use crate::conformance::subject::{self, Registered};
use crate::conformance::{corpus, corpus_dir, metric, report};
use stats::Tally;

const USAGE: &str = "usage: cargo xtask compare --a NAME --b NAME [--fn ID] [--precision f64|f32]";

#[derive(Debug, Default, PartialEq)]
struct Options {
    a: String,
    b: String,
    fn_id: Option<String>,
    precision: Option<Precision>,
}

fn parse_args(args: &[String]) -> Result<Options, String> {
    let mut o = Options::default();
    let mut rest = args.iter();
    while let Some(arg) = rest.next() {
        let (flag, inline) = match arg.split_once('=') {
            Some((f, v)) => (f, Some(v.to_string())),
            None => (arg.as_str(), None),
        };
        if !matches!(flag, "--a" | "--b" | "--fn" | "--precision") {
            return Err(format!("unknown argument `{arg}`; {USAGE}"));
        }
        let value = inline.or_else(|| rest.next().cloned());
        let value = value.ok_or_else(|| format!("`{flag}` needs a value; {USAGE}"))?;
        match flag {
            "--a" => o.a = value,
            "--b" => o.b = value,
            "--fn" => o.fn_id = Some(value),
            _ => {
                o.precision = Some(match value.as_str() {
                    "f64" => Precision::F64,
                    "f32" => Precision::F32,
                    _ => return Err(format!("unknown precision `{value}`; {USAGE}")),
                })
            }
        }
    }
    if o.a.is_empty() || o.b.is_empty() {
        return Err(format!("`--a` and `--b` are both required; {USAGE}"));
    }
    Ok(o)
}

/// One id's strata, in order of first appearance, each with its tally.
struct IdResult {
    fn_id: String,
    strata: Vec<(String, Tally)>,
}

impl IdResult {
    /// The strata pooled (`0038` decision 5: per id first).
    fn pooled(&self) -> Tally {
        let mut t = Tally::default();
        for (_, s) in &self.strata {
            t.merge(s);
        }
        t
    }
}

/// `subject` registered under `name`, or the names there are.
fn pick(all: &mut Vec<Registered>, name: &str) -> Result<Registered, String> {
    match all.iter().position(|s| s.subject.name() == name) {
        Some(i) => Ok(all.remove(i)),
        None => {
            let names: Vec<&str> = all.iter().map(|s| s.subject.name()).collect();
            Err(format!(
                "no in-process subject `{name}` is registered; have: {}",
                names.join(", ")
            ))
        }
    }
}

/// Both subjects over `entries`: each file read once, every record scored on both sides. An id
/// either subject does not support is left out and named in the second list. At `f32` only the
/// `@f32` strata are scored, as in `conformance`, and a subject without an `f32` kernel is an
/// error.
fn compare(
    dir: &Path,
    entries: &[corpus::Entry],
    a: &Registered,
    b: &Registered,
    precision: Precision,
) -> Result<(Vec<IdResult>, Vec<String>), String> {
    let f32 = precision == Precision::F32;
    if f32 {
        for s in [a, b] {
            if let Some(why) = &s.no_f32 {
                return Err(format!("{} has no f32 kernel: {why}", s.subject.name()));
            }
        }
    }
    let (mut out, mut left_out) = (Vec::new(), Vec::new());
    for entry in entries {
        let fn_id = entry.fn_id.as_str();
        if !(a.subject.supports(fn_id) && b.subject.supports(fn_id)) {
            if a.subject.supports(fn_id) || b.subject.supports(fn_id) {
                left_out.push(entry.fn_id.clone());
            }
            continue;
        }
        let rule = metric::rule(fn_id).ok_or_else(|| {
            format!("no metric rule for `{fn_id}`: add its row to `metric::TABLE`")
        })?;
        let mut strata: Vec<(String, Tally)> = Vec::new();
        let mut index: BTreeMap<String, usize> = BTreeMap::new();
        for record in corpus::read(dir, entry)?
            .iter()
            .filter(|r| r.is_f32_stratum() == f32)
        {
            if f32 && !corpus::converts(fn_id) {
                record
                    .require_binary32()
                    .map_err(|e| format!("{fn_id} record {}: {e}", record.id))?;
            }
            let score = |s: &Registered| {
                let output = s.subject.eval(fn_id, record, precision);
                rule.score(record, &output, precision).map_err(|e| {
                    format!("{} on {fn_id} record {}: {e}", s.subject.name(), record.id)
                })
            };
            let (sa, sb) = (score(a)?, score(b)?);
            let at = *index.entry(record.stratum.clone()).or_insert_with(|| {
                strata.push((record.stratum.clone(), Tally::default()));
                strata.len() - 1
            });
            strata[at].1.add(sa, sb);
        }
        if strata.is_empty() {
            return Err(format!(
                "`{fn_id}` has no record at {}: only the ids with `@f32` strata compare at f32 \
                 (docs/decisions/0016 item 2)",
                report::precision_name(precision)
            ));
        }
        for (_, t) in &mut strata {
            t.close_stratum();
        }
        out.push(IdResult {
            fn_id: entry.fn_id.clone(),
            strata,
        });
    }
    Ok((out, left_out))
}

/// `b`, `a` or `=`: which subject is better on more records.
fn direction(t: &Tally) -> &'static str {
    match t.better.cmp(&t.worse) {
        std::cmp::Ordering::Greater => "b",
        std::cmp::Ordering::Less => "a",
        std::cmp::Ordering::Equal => "=",
    }
}

fn ratio(x: Option<f64>) -> String {
    x.map_or_else(|| "-".to_string(), |r| format!("{r:.3}"))
}

fn max_u(x: Option<f64>) -> String {
    x.map_or_else(|| "-".to_string(), |r| format!("{r:.4e}"))
}

/// Three decimals down to 0.01, two significant digits below; `0` is a `p` past binary64's range.
fn p_value(p: f64) -> String {
    match p {
        0.0 => "<1e-307".to_string(),
        p if p >= 0.01 => format!("{p:.3}"),
        p => format!("{p:.2e}"),
    }
}

const HEADER: [&str; 15] = [
    "fn",
    "stratum",
    "n",
    "b<a",
    "b>a",
    "tie",
    "dir",
    "p_sign",
    "gmean_b/a",
    "worst_b/a",
    "worst_a/b",
    "max_a_u",
    "max_b_u",
    "max_b/max_a",
    "other",
];

/// One row of the table: a stratum, an id pooled, or the comparison pooled.
fn cells(fn_id: &str, stratum: &str, t: &Tally) -> [String; 15] {
    let mut other = Vec::new();
    for (count, what) in [
        (t.zero_sided, "zero-sided"),
        (t.nonfinite_a, "nonfinite-a"),
        (t.nonfinite_b, "nonfinite-b"),
    ] {
        if count > 0 {
            other.push(format!("{count} {what}"));
        }
    }
    [
        fn_id.to_string(),
        stratum.to_string(),
        t.n.to_string(),
        t.better.to_string(),
        t.worse.to_string(),
        t.ties.to_string(),
        direction(t).to_string(),
        p_value(t.sign_p()),
        ratio(t.geomean()),
        ratio(t.worst_b),
        ratio(t.worst_a),
        max_u(t.max_a),
        max_u(t.max_b),
        ratio(t.max_ratio()),
        other.join(", "),
    ]
}

/// The strata tally of a pooled row: how many strata have `max_b / max_a` below, above and at 1.
fn strata_line(t: &Tally) -> String {
    format!(
        "strata by max_b/max_a: {} b better, {} b worse, {} equal; max_b/max_a from {} to {}",
        t.strata_better,
        t.strata_worse,
        t.strata_tied,
        ratio(t.strata_ratio.map(|r| r.0)),
        ratio(t.strata_ratio.map(|r| r.1)),
    )
}

/// The report: a header naming the two subjects, one row per stratum, a pooled row per id with
/// its strata tally, and the pooled row over every id beside them, never instead (decision 5).
fn render(a: &Registered, b: &Registered, precision: Precision, ids: &[IdResult]) -> String {
    let p = report::precision_name(precision);
    let mut rows: Vec<[String; 15]> = Vec::new();
    let mut tallies: Vec<String> = Vec::new();
    let mut all = Tally::default();
    for id in ids {
        for (stratum, t) in &id.strata {
            rows.push(cells(&id.fn_id, stratum, t));
        }
        let pooled = id.pooled();
        let label = format!("(id: {} strata)", id.strata.len());
        rows.push(cells(&id.fn_id, &label, &pooled));
        tallies.push(format!("{}: {}", id.fn_id, strata_line(&pooled)));
        all.merge(&pooled);
    }
    let label = format!("(all: {} ids pooled)", ids.len());
    rows.push(cells("*", &label, &all));
    tallies.push(format!("*: {}", strata_line(&all)));

    let width: Vec<usize> = (0..HEADER.len())
        .map(|c| {
            rows.iter()
                .map(|r| r[c].len())
                .fold(HEADER[c].len(), usize::max)
        })
        .collect();
    let mut out = format!(
        "compare at {p}: a = {} {}, b = {} {}\n\
         b<a: records where b's error is smaller; gmean_b/a over the records that differ, both \
         errors positive; worst_b/a and worst_a/b the largest ratio each way; ids pooled share \
         kernels, so the last row is not independent evidence (0038 decision 5)\n",
        a.subject.name(),
        a.version_at(precision),
        b.subject.name(),
        b.version_at(precision),
    );
    let mut line = |cols: &mut dyn Iterator<Item = &str>| {
        for (c, text) in cols.enumerate() {
            let w = width[c];
            let _ = if c < 2 || c + 1 == HEADER.len() {
                write!(out, "{text:<w$}  ")
            } else {
                write!(out, "{text:>w$}  ")
            };
        }
        out.truncate(out.trim_end().len());
        out.push('\n');
    };
    line(&mut HEADER.into_iter());
    for r in &rows {
        line(&mut r.iter().map(String::as_str));
    }
    out.push('\n');
    for t in &tallies {
        out.push_str(t);
        out.push('\n');
    }
    out
}

/// Selects by `options`, compares and returns the report.
fn run_with(dir: &Path, options: &Options, mut all: Vec<Registered>) -> Result<String, String> {
    let precision = options.precision.unwrap_or(Precision::F64);
    let mut entries = corpus::manifest(dir)?;
    if let Some(id) = &options.fn_id {
        entries.retain(|e| e.fn_id == *id);
        if entries.is_empty() {
            return Err(format!("no corpus file for `--fn {id}`"));
        }
    }
    let a = pick(&mut all, &options.a)?;
    // `--a X --b X` is a control: every record ties. The registry is rebuilt for the second pick.
    let b = match options.a == options.b {
        true => pick(&mut subject::registry(), &options.b)?,
        false => pick(&mut all, &options.b)?,
    };
    if precision == Precision::F32 && options.fn_id.is_none() {
        let mut with = Vec::new();
        for e in entries {
            if corpus::mentions_f32(dir, &e)? {
                with.push(e);
            }
        }
        entries = with;
    }
    let (ids, left_out) = compare(dir, &entries, &a, &b, precision)?;
    if !left_out.is_empty() {
        eprintln!(
            "compare: answered by one subject only, left out: {}",
            left_out.join(", ")
        );
    }
    if ids.is_empty() {
        return Err(format!(
            "nothing to compare: no corpus id at {} is answered by both `{}` and `{}`",
            report::precision_name(precision),
            options.a,
            options.b
        ));
    }
    Ok(render(&a, &b, precision, &ids))
}

#[allow(clippy::print_stdout)]
pub(crate) fn run(args: &[String]) -> Result<(), String> {
    let options = parse_args(args)?;
    let report = run_with(&corpus_dir()?, &options, subject::registry())?;
    print!("{report}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::metric::Score;

    fn args(a: &[&str]) -> Vec<String> {
        a.iter().map(ToString::to_string).collect()
    }

    #[test]
    fn arguments_take_both_spellings_and_refuse_the_rest() -> Result<(), String> {
        let o = parse_args(&args(&[
            "--a",
            "x",
            "--b=y",
            "--fn",
            "f",
            "--precision=f32",
        ]))?;
        let want = Options {
            a: "x".into(),
            b: "y".into(),
            fn_id: Some("f".into()),
            precision: Some(Precision::F32),
        };
        assert_eq!(o, want);
        for bad in [
            &["--a", "x"][..],
            &["--b", "y"],
            &["--a", "x", "--b"],
            &["--a", "x", "--b", "y", "--subject", "z"],
            &["--a", "x", "--b", "y", "--precision", "f16"],
            &["x"],
            &[],
        ] {
            assert!(parse_args(&args(bad)).is_err(), "{bad:?}");
        }
        Ok(())
    }

    #[test]
    fn a_subject_against_itself_ties_on_every_record() -> Result<(), String> {
        let o = parse_args(&args(&[
            "--a", "helicoid", "--b", "helicoid", "--fn", "so3_exp",
        ]))?;
        let report = run_with(&corpus_dir()?, &o, subject::registry())?;
        let pooled = report
            .lines()
            .find(|l| l.starts_with("so3_exp") && l.contains("(id:"))
            .unwrap_or_default();
        let cols: Vec<&str> = pooled.split_whitespace().collect();
        // n, b<a, b>a, tie after the four-word label.
        let n = cols.get(4).copied().unwrap_or_default();
        assert_eq!(cols.get(5..8), Some(&["0", "0", n][..]), "{pooled}");
        assert!(report.contains("so3_exp: strata by max_b/max_a: 0 b better, 0 b worse"));
        Ok(())
    }

    #[test]
    fn an_unknown_subject_and_an_id_neither_answers_are_errors() -> Result<(), String> {
        let dir = corpus_dir()?;
        let run = |a: &[&str]| run_with(&dir, &parse_args(&args(a))?, subject::registry());
        let e = run(&["--a", "helicoid", "--b", "nope"])
            .err()
            .unwrap_or_default();
        assert!(e.contains("no in-process subject `nope`"), "{e}");
        let e = run(&[
            "--a",
            "helicoid",
            "--b",
            "helicoid:unit-r",
            "--fn",
            "so3_exp",
        ]);
        assert!(e.err().unwrap_or_default().contains("nothing to compare"));
        let e = run(&[
            "--a",
            "helicoid",
            "--b",
            "helicoid:host-std",
            "--precision",
            "f32",
        ]);
        assert!(e.err().unwrap_or_default().contains("has no f32 kernel"));
        let e = run(&[
            "--a",
            "helicoid",
            "--b",
            "helicoid",
            "--fn",
            "so3_exp",
            "--precision=f32",
        ]);
        assert!(e.err().unwrap_or_default().contains("no record at f32"));
        Ok(())
    }

    #[test]
    fn the_report_pools_per_id_and_then_beside_the_ids() {
        let tally = |pairs: &[(f64, f64)]| {
            let mut t = Tally::default();
            for &(a, b) in pairs {
                t.add(Score::Finite(a), Score::Finite(b));
            }
            t.close_stratum();
            t
        };
        let ids = [
            IdResult {
                fn_id: "f".into(),
                strata: vec![
                    ("s1".into(), tally(&[(2.0, 1.0), (2.0, 1.0)])),
                    ("s2".into(), tally(&[(1.0, 1.0)])),
                ],
            },
            IdResult {
                fn_id: "g".into(),
                strata: vec![("s1".into(), tally(&[(1.0, 3.0)]))],
            },
        ];
        let reg = || subject::registry().remove(0);
        let text = render(&reg(), &reg(), Precision::F64, &ids);
        let row = |fn_id: &str, stratum: &str| {
            text.lines()
                .map(|l| l.split_whitespace().collect::<Vec<_>>())
                .find(|c| c.starts_with(&[fn_id, stratum]))
                .unwrap_or_default()
        };
        assert_eq!(row("f", "s1")[2..7], ["2", "2", "0", "0", "b"]);
        assert_eq!(row("f", "(id:")[4..9], ["3", "2", "0", "1", "b"]);
        assert_eq!(row("*", "(all:")[5..10], ["4", "2", "1", "1", "b"]);
        assert!(text.contains("f: strata by max_b/max_a: 1 b better, 0 b worse, 1 equal"));
        assert!(text.contains("g: strata by max_b/max_a: 0 b better, 1 b worse, 0 equal"));
        assert!(text.contains("*: strata by max_b/max_a: 1 b better, 1 b worse, 1 equal"));
    }

    /// The task's own control (`0078` decision 5): `0038`'s tally of the `W²` form twin by the
    /// maximum over the nine `sen3_*` ids, reproduced from the per-stratum rows. `0038` published 21
    /// strata better and 104 worse; since `0089` (draft) moved the seeded kernels' closed arms the
    /// same comparison reads 36 and 94, every id still net worse, and the `sen3_exp` ratio floors
    /// published as 0.73 and 0.64 read 0.77, 0.96 and 0.90.
    #[test]
    fn the_task_reproduces_the_tally_of_0038() -> Result<(), String> {
        let dir = corpus_dir()?;
        let mut all = subject::registry();
        let a = pick(&mut all, "seeded:correct")?;
        let b = pick(&mut all, "seeded:w2-identity")?;
        let entries: Vec<corpus::Entry> = corpus::manifest(&dir)?
            .into_iter()
            .filter(|e| e.fn_id.starts_with("sen3_"))
            .collect();
        let (ids, _) = compare(&dir, &entries, &a, &b, Precision::F64)?;
        assert_eq!(ids.len(), 9);
        let mut pooled = Tally::default();
        for id in &ids {
            let t = id.pooled();
            // "Every one of the nine ids is net worse on its own."
            assert!(t.worse > t.better, "{}", id.fn_id);
            pooled.merge(&t);
        }
        assert_eq!((pooled.strata_better, pooled.strata_worse), (36, 94));
        let low: Vec<String> = ids
            .iter()
            .filter(|i| i.fn_id.starts_with("sen3_exp"))
            .filter_map(|i| i.pooled().strata_ratio)
            .map(|(lo, _)| format!("{lo:.2}"))
            .collect();
        assert_eq!(low, ["0.77", "0.96", "0.90"]);
        Ok(())
    }

    /// The twin's positive control (`0078` decision 3): `0076` decision 1's hand tally, made by
    /// editing `Quat::to_matrix` and rerunning `conformance`, reproduced by this task to the
    /// stratum and to the two decimals of its ratio column.
    #[test]
    fn the_unit_r_twin_reproduces_the_tally_of_0076() -> Result<(), String> {
        let dir = corpus_dir()?;
        let mut all = subject::registry();
        let (a, b) = (
            pick(&mut all, "helicoid")?,
            pick(&mut all, "helicoid:unit-r")?,
        );
        let entries: Vec<corpus::Entry> = corpus::manifest(&dir)?
            .into_iter()
            .filter(|e| b.subject.supports(&e.fn_id))
            .collect();
        // (precision, ids pooled, strata better, worse, ratio range) as `0076` prints them.
        let table = [
            (Precision::F64, "sen3_ad_n", 22, 20, "0.28-2.12"),
            (Precision::F64, "sim3_ad", 36, 51, "0.28-2.02"),
            (Precision::F32, "sim3_ad", 43, 32, "0.45-1.94"),
        ];
        for (precision, prefix, better, worse, range) in table {
            let f32 = precision == Precision::F32;
            let mut chosen = Vec::new();
            for e in &entries {
                if !f32 || corpus::mentions_f32(&dir, e)? {
                    chosen.push(corpus::Entry {
                        fn_id: e.fn_id.clone(),
                        records: e.records,
                    });
                }
            }
            let (ids, _) = compare(&dir, &chosen, &a, &b, precision)?;
            let mut pooled = Tally::default();
            for id in ids.iter().filter(|i| i.fn_id.starts_with(prefix)) {
                let t = id.pooled();
                // "22 each": every width alone, not only their sum.
                assert_eq!(
                    (t.strata_better, t.strata_worse),
                    (better, worse),
                    "{}",
                    id.fn_id
                );
                pooled.merge(&t);
            }
            let (lo, hi) = pooled.strata_ratio.unwrap_or_default();
            assert_eq!(format!("{lo:.2}-{hi:.2}"), range, "{prefix} {precision:?}");
        }
        Ok(())
    }
}
