//! `cargo xtask determinism` (`docs/PHASE6.md` §3, D16): the shipped subject's outputs over the
//! whole corpus, digested per function id and precision, and the digests of several targets
//! compared.
//!
//! `--write FILE` evaluates `helicoid` (`crate::shipped`) on every record it supports, the plain
//! strata at `f64` and the `@f32` strata at `f32`, and writes one line per (id, precision): the id,
//! the precision, the record count and the SHA-256 of the outputs. Each output value enters the hash
//! as its binary64 bit pattern in hexadecimal, so the digest sees every bit; a NaN enters as `nan`,
//! because D16 leaves NaN sign and payload apart. `--compare A B …` fails on any line that differs,
//! naming its id, and on an id one file has and another lacks.
//!
//! Nothing here is scored: a digest is evidence of identity, never of accuracy (`just envelope`).

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::Path;

use helicoid_linalg::Precision;
use sha2::{Digest, Sha256};

use crate::conformance::corpus::{self, Record};
use crate::conformance::corpus_dir;
use crate::conformance::subject::Output;

const USAGE: &str = "usage: cargo xtask determinism --write FILE | --compare FILE FILE [FILE …]";

pub(crate) fn run(args: &[String]) -> Result<(), String> {
    match args {
        [flag, file] if flag == "--write" => {
            let table = digests(&corpus_dir()?)?;
            std::fs::write(file, table).map_err(|e| format!("{file}: {e}"))?;
            eprintln!("determinism: wrote {file}");
            Ok(())
        }
        [flag, files @ ..] if flag == "--compare" && files.len() >= 2 => {
            let tables = files
                .iter()
                .map(|f| {
                    let text = std::fs::read_to_string(f).map_err(|e| format!("{f}: {e}"))?;
                    Ok((f.as_str(), text))
                })
                .collect::<Result<Vec<_>, String>>()?;
            compare(&tables)
        }
        _ => Err(USAGE.into()),
    }
}

/// The digest table of the shipped subject over the corpus in `dir`.
fn digests(dir: &Path) -> Result<String, String> {
    let registered = crate::shipped::registered();
    let subject = registered.subject.as_ref();
    let mut table = format!(
        "# helicoid {} on {}: fn precision records sha256\n",
        registered.version,
        env!("XTASK_TARGET"),
    );
    for entry in corpus::manifest(dir)? {
        let fn_id = entry.fn_id.as_str();
        if !subject.supports(fn_id) {
            continue;
        }
        let records = corpus::read(dir, &entry)?;
        for precision in [Precision::F64, Precision::F32] {
            let f32 = precision == Precision::F32;
            let mut hash = Sha256::new();
            let mut n = 0_usize;
            for record in records.iter().filter(|r| r.is_f32_stratum() == f32) {
                if f32 {
                    record
                        .require_binary32()
                        .map_err(|e| format!("{fn_id} record {}: {e}", record.id))?;
                }
                hash.update(line(record, &subject.eval(fn_id, record, precision)));
                n += 1;
            }
            if n > 0 {
                let hex: String = hash.finalize().iter().map(|b| format!("{b:02x}")).collect();
                let name = if f32 { "f32" } else { "f64" };
                let _ = writeln!(table, "{fn_id} {name} {n} {hex}");
            }
        }
    }
    Ok(table)
}

/// One record's outputs: its id, then each field and its values' bit patterns, `nan` for any NaN.
fn line(record: &Record, output: &Output) -> String {
    let mut s = record.id.to_string();
    for (field, values) in output {
        let _ = write!(s, " {field}");
        for v in values {
            if v.is_nan() {
                s.push_str(" nan");
            } else {
                let _ = write!(s, " {:016x}", v.to_bits());
            }
        }
    }
    s.push('\n');
    s
}

/// Fails on any digest line on which the tables disagree, naming the id and the files.
fn compare(tables: &[(&str, String)]) -> Result<(), String> {
    // (fn, precision) -> file -> "records sha256"
    let mut seen: BTreeMap<(String, String), BTreeMap<&str, String>> = BTreeMap::new();
    for (file, text) in tables {
        for row in text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let cells: Vec<&str> = row.split_whitespace().collect();
            let [fn_id, precision, n, sha] = cells[..] else {
                return Err(format!("{file}: malformed line `{row}`"));
            };
            let key = (fn_id.to_string(), precision.to_string());
            seen.entry(key)
                .or_default()
                .insert(file, format!("{n} {sha}"));
        }
    }
    let mut differ = Vec::new();
    for ((fn_id, precision), by_file) in &seen {
        let mut values = by_file.values();
        let first = values.next();
        if by_file.len() != tables.len() || values.any(|v| Some(v) != first) {
            let files: Vec<&str> = tables
                .iter()
                .map(|(f, _)| *f)
                .filter(|f| !by_file.contains_key(f))
                .collect();
            let missing = match files.is_empty() {
                true => String::new(),
                false => format!(" (absent from {})", files.join(", ")),
            };
            differ.push(format!("  {fn_id} at {precision}{missing}"));
        }
    }
    for (file, text) in tables {
        let header = text.lines().next().unwrap_or_default();
        eprintln!("determinism: {file}: {header}");
    }
    match differ.is_empty() {
        true => {
            eprintln!(
                "determinism: {} digests identical across {} files",
                seen.len(),
                tables.len()
            );
            Ok(())
        }
        false => Err(format!(
            "{} of {} digests differ (D16):\n{}",
            differ.len(),
            seen.len(),
            differ.join("\n")
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(rows: &[&str]) -> String {
        format!("# header\n{}\n", rows.join("\n"))
    }

    #[test]
    fn identical_tables_pass_and_any_difference_is_named() {
        let a = table(&["so3_exp f64 10 aa", "so3_exp f32 4 bb"]);
        assert_eq!(compare(&[("a", a.clone()), ("b", a.clone())]), Ok(()));
        let b = table(&["so3_exp f64 10 aa", "so3_exp f32 4 bc"]);
        let err = compare(&[("a", a.clone()), ("b", b)])
            .err()
            .unwrap_or_default();
        assert!(
            err.contains("so3_exp at f32") && !err.contains("at f64"),
            "{err}"
        );
        let c = table(&["so3_exp f64 10 aa"]);
        let err = compare(&[("a", a), ("c", c)]).err().unwrap_or_default();
        assert!(err.contains("so3_exp at f32 (absent from c)"), "{err}");
        assert!(compare(&[("a", table(&["one two"])), ("b", table(&[]))]).is_err());
    }

    #[test]
    fn a_nan_is_one_value_and_every_other_bit_counts() -> Result<(), String> {
        let record = corpus::parse_line(
            r#"{"id":7,"in":{"x":["0x1p+0"]},"out":{"y":["1"]},"stratum":"s"}"#,
        )?;
        let out = |v: f64| Output::from([("y".to_string(), vec![v])]);
        assert_eq!(
            line(&record, &out(f64::NAN)),
            line(&record, &out(-f64::NAN))
        );
        assert_eq!(line(&record, &out(1.0)), "7 y 3ff0000000000000\n");
        assert_ne!(line(&record, &out(0.0)), line(&record, &out(-0.0)));
        Ok(())
    }

    #[test]
    fn usage_is_refused() {
        let args = |a: &[&str]| a.iter().map(ToString::to_string).collect::<Vec<_>>();
        assert!(run(&args(&[])).is_err());
        assert!(run(&args(&["--compare", "one"])).is_err());
        assert!(run(&args(&["--write"])).is_err());
    }
}
