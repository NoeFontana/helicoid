//! The nalgebra oracle (`docs/PHASE1.md` §7, `docs/decisions/0056` decision 4): reads
//! corpus JSONL and writes one answer file per function id it supports. It never reads a reference
//! and never scores; `cargo xtask conformance --oracle nalgebra` does, with the harness's exact
//! metric.
//!
//! ```text
//! nalgebra_runner --version               the pin, one line: the result rows' subject_version
//! nalgebra_runner --out DIR FILE.jsonl…   DIR/<id>.jsonl for each FILE whose id is supported
//! ```
//!
//! An answer line is `{"id":N,"out":{"<field>":["<hex float>",…]}}`, `N` the record's `id`, fields
//! flat (matrices column-major, no `shape`), values as in [`hexfloat`].

mod convert;
mod hexfloat;

use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;
use std::process::ExitCode;

use serde_json::{json, Map, Value};

use convert::{answer, supported, Fields};

/// The oracle, as `Cargo.toml` pins it.
const PIN: &str = "nalgebra@0.35.0";

/// The `in` object of a corpus record: every entry is a hex-float string or an array of them,
/// except a matrix's sibling `shape`, two integers, which is kept as two numbers.
fn inputs(object: &Map<String, Value>) -> Result<Fields, String> {
    let mut fields = BTreeMap::new();
    for (key, value) in object {
        let one = |v: &Value| match (key.as_str(), v) {
            ("shape", v) => v
                .as_u64()
                .and_then(|n| u32::try_from(n).ok())
                .map(f64::from)
                .ok_or_else(|| format!("`shape`: not a dimension: {v}")),
            (_, Value::String(s)) => hexfloat::parse(s),
            (_, v) => Err(format!("`{key}`: not a string: {v}")),
        };
        let values = match value {
            Value::Array(items) => items.iter().map(one).collect::<Result<_, _>>()?,
            v => vec![one(v)?],
        };
        fields.insert(key.clone(), values);
    }
    Ok(fields)
}

/// The answer file of `fn_id` for the corpus file `text`: one line per record, in order.
fn answer_file(fn_id: &str, text: &str) -> Result<String, String> {
    let mut out = String::with_capacity(text.len() / 2);
    for (i, line) in text.lines().enumerate() {
        let at = |e: String| format!("line {}: {e}", i + 1);
        let record: Value = serde_json::from_str(line).map_err(|e| at(e.to_string()))?;
        let id = record
            .get("id")
            .and_then(Value::as_u64)
            .ok_or_else(|| at("no `id`".into()))?;
        let object = record
            .get("in")
            .and_then(Value::as_object)
            .ok_or_else(|| at("no `in` object".into()))?;
        let given = inputs(object).map_err(at)?;
        let answered = answer(fn_id, &given)
            .ok_or_else(|| format!("`{fn_id}` is not supported"))?
            .map_err(at)?;
        let out_fields: Map<String, Value> = answered
            .into_iter()
            .map(|(k, v)| {
                (
                    k,
                    v.iter()
                        .map(|&x| Value::from(hexfloat::format(x)))
                        .collect(),
                )
            })
            .collect();
        out.push_str(&json!({ "id": id, "out": out_fields }).to_string());
        out.push('\n');
    }
    Ok(out)
}

fn run(args: &[String], stdout: &mut impl Write, stderr: &mut impl Write) -> Result<(), String> {
    let io = |e: std::io::Error| e.to_string();
    match args {
        [flag] if flag == "--version" => writeln!(stdout, "{PIN}").map_err(io),
        [flag, dir, files @ ..] if flag == "--out" && !files.is_empty() => {
            let dir = Path::new(dir);
            std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
            let (mut answered, mut skipped) = (Vec::new(), 0usize);
            for file in files {
                let path = Path::new(file);
                let fn_id = path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or_default();
                if !supported().any(|id| id == fn_id) {
                    skipped += 1;
                    continue;
                }
                let text = std::fs::read_to_string(path).map_err(|e| format!("{file}: {e}"))?;
                let lines = answer_file(fn_id, &text).map_err(|e| format!("{file}: {e}"))?;
                let target = dir.join(format!("{fn_id}.jsonl"));
                std::fs::write(&target, lines).map_err(|e| format!("{}: {e}", target.display()))?;
                answered.push(fn_id);
            }
            writeln!(
                stderr,
                "answered {answered:?}; {skipped} files of other ids skipped"
            )
            .map_err(io)
        }
        _ => Err("usage: nalgebra_runner --version | --out DIR FILE.jsonl…".into()),
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match run(&args, &mut std::io::stdout(), &mut std::io::stderr()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            // Nothing else to do with a failed write to stderr.
            let _ = writeln!(std::io::stderr(), "nalgebra_runner: {e}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CORPUS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../conformance/corpus");

    /// A corpus line as the generator writes it, its reference elided to what the runner ignores.
    const LINE: &str = r#"{"id":0,"in":{"A":["0x1.0000000000000p+2","0x1.0000000000000p+1","0x1.0000000000000p+1","0x1.0000000000000p+1","0x1.4000000000000p+2","0x1.8000000000000p+1","0x1.0000000000000p+1","0x1.8000000000000p+1","0x1.8000000000000p+2"],"shape":[3,3]},"out":{"valid":["1e0"]},"stratum":"chol:spd"}"#;

    #[test]
    fn a_record_is_answered_with_its_id_and_hex_floats() -> Result<(), String> {
        let out = answer_file("chol_n3", LINE)?;
        let line: Value = serde_json::from_str(out.trim_end()).map_err(|e| e.to_string())?;
        assert_eq!(line["id"], 0);
        let floats = |key: &str| -> Result<Vec<f64>, String> {
            line["out"][key]
                .as_array()
                .ok_or(format!("no {key}"))?
                .iter()
                .map(|v| hexfloat::parse(v.as_str().unwrap_or_default()))
                .collect()
        };
        // `[[4, 2, 2], [2, 5, 3], [2, 3, 6]]`, whose factor is exact.
        let l = [2.0, 1.0, 1.0, 0.0, 2.0, 1.0, 0.0, 0.0, 2.0];
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&floats("L")?), bits(&l));
        assert_eq!(bits(&floats("valid")?), bits(&[1.0]));
        assert!(line["out"].get("shape").is_none());
        assert_eq!(out.matches('\n').count(), 1);
        Ok(())
    }

    #[test]
    fn a_bad_line_or_an_unsupported_id_is_an_error_naming_the_line() {
        let e = answer_file("chol_n3", &format!("{LINE}\n{{}}\n"))
            .err()
            .unwrap_or_default();
        assert!(e.contains("line 2"), "{e}");
        let bad = LINE.replacen("0x1.0000000000000p+2", "4.0", 1);
        assert!(answer_file("chol_n3", &bad)
            .err()
            .unwrap_or_default()
            .contains("bad hex float"));
        let bad = LINE.replace("[3,3]", "[3,\"3\"]");
        assert!(answer_file("chol_n3", &bad)
            .err()
            .unwrap_or_default()
            .contains("not a dimension"));
        assert!(answer_file("chol_n6", LINE).is_err());
        assert!(answer_file("so3_exp", LINE).is_err());
    }

    /// Every input of every supported corpus file is `float.hex()` to the byte, every matrix of an
    /// id nalgebra answers from one triangle exactly symmetric (`mat2_inverse_adj`'s are general),
    /// and every record answered without a panic. A non-finite answer is the harness's to record (§7), not this test's to refuse. The
    /// corpus is committed, so this reads it as it is.
    #[test]
    fn the_whole_supported_corpus_is_answered() -> Result<(), String> {
        for fn_id in supported() {
            let path = format!("{CORPUS}/{fn_id}.jsonl");
            let text = std::fs::read_to_string(&path).map_err(|e| format!("{path}: {e}"))?;
            for line in text.lines() {
                let record: Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
                let object = record["in"].as_object().ok_or("in")?;
                for (key, value) in object.iter().filter(|(k, _)| *k != "shape") {
                    let items = match value {
                        Value::Array(items) => items.clone(),
                        v => vec![v.clone()],
                    };
                    for s in &items {
                        let s = s.as_str().ok_or("string")?;
                        assert_eq!(hexfloat::format(hexfloat::parse(s)?), s, "{fn_id} {key}");
                    }
                }
                let given = inputs(object)?;
                if let Some(a) = given.get("A").filter(|_| fn_id != "mat2_inverse_adj") {
                    let n = (a.len() as f64).sqrt() as usize;
                    for (i, j) in (0..n).flat_map(|i| (0..n).map(move |j| (i, j))) {
                        assert_eq!(a[i + n * j].to_bits(), a[j + n * i].to_bits(), "{fn_id}");
                    }
                }
            }
            let answers = answer_file(fn_id, &text)?;
            assert_eq!(answers.lines().count(), text.lines().count(), "{fn_id}");
        }
        Ok(())
    }

    #[test]
    fn the_pin_is_the_version_of_the_manifest_and_of_the_lockfile() -> Result<(), String> {
        let version = PIN.strip_prefix("nalgebra@").ok_or("PIN")?;
        let pinned = format!("nalgebra = {{ version = \"={version}\"");
        assert!(include_str!("../Cargo.toml").contains(&pinned));
        let locked = format!("name = \"nalgebra\"\nversion = \"{version}\"\nsource = \"registry+");
        assert!(include_str!("../Cargo.lock").contains(&locked));
        Ok(())
    }

    #[test]
    fn the_command_line_is_version_or_out_with_files() {
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let mut go = |args: &[&str]| {
            let args: Vec<String> = args.iter().map(ToString::to_string).collect();
            run(&args, &mut out, &mut err)
        };
        assert!(go(&["--version"]).is_ok());
        for bad in [&[][..], &["--out", "d"], &["--nope"], &["x.jsonl"]] {
            assert!(go(bad).is_err(), "{bad:?}");
        }
        assert_eq!(String::from_utf8_lossy(&out), format!("{PIN}\n"));
    }
}
