//! The `tf_tree_math` oracle (`docs/PHASE1.md` §7, `docs/decisions/0010`): reads corpus JSONL and
//! writes one answer file per function id it supports. It never reads a reference and never
//! scores; `cargo xtask conformance --oracle tf_tree_math` does, with the harness's exact metric.
//!
//! ```text
//! tf_tree_math_runner --version               the pin, one line: the result rows' subject_version
//! tf_tree_math_runner --out DIR FILE.jsonl…   DIR/<id>.jsonl for each FILE whose id is supported
//! ```
//!
//! An answer line is `{"id":N,"out":{"<field>":["<hex float>",…]}}`, `N` the record's `id`, fields
//! flat in the reference's order (matrices column-major, no `shape`), values as in [`hexfloat`].

mod convert;
mod hexfloat;

use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;
use std::process::ExitCode;

use serde_json::{json, Map, Value};

use convert::{answer, supported, Fields};

/// The oracle, as `Cargo.toml` pins it: the `tf_tree` commit that `tf_tree_math` is built from.
const PIN: &str = "tf_tree@20bc5a0518ec791318777d2105ab3296047d84ec";

/// The `in` object of a corpus record: every entry is a hex-float string or an array of them,
/// except `shape`.
///
/// `shape` is the sibling a matrix field carries (`docs/PHASE1.md` §4.3) and holds two JSON
/// integers, not values; a record holds at most one matrix, so the dimensions are the field's own
/// and the conversion that reads it names the layout (`convert::helicoid_to_tf_tree_rot3`).
/// Parsing it as a hex float is what made `so3_from_matrix` fail with `string` before it was
/// skipped here.
const SHAPE: &str = "shape";

fn inputs(object: &Map<String, Value>) -> Result<Fields, String> {
    let mut fields = BTreeMap::new();
    for (key, value) in object.iter().filter(|(k, _)| *k != SHAPE) {
        let one = |v: &Value| match v.as_str() {
            Some(s) => hexfloat::parse(s),
            None => Err(format!("`{key}`: not a string: {v}")),
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
        _ => Err("usage: tf_tree_math_runner --version | --out DIR FILE.jsonl…".into()),
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match run(&args, &mut std::io::stdout(), &mut std::io::stderr()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            // Nothing else to do with a failed write to stderr.
            let _ = writeln!(std::io::stderr(), "tf_tree_math_runner: {e}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CORPUS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../conformance/corpus");

    /// A corpus line as the generator writes it, its reference elided to what the runner ignores.
    const LINE: &str = r#"{"id":0,"in":{"phi":["0x0.0p+0","0x0.0p+0","0x1.921fb54442d18p+0"]},"out":{"q":["1e0"]},"stratum":"s"}"#;

    #[test]
    fn a_record_is_answered_with_its_id_and_hex_floats() -> Result<(), String> {
        let out = answer_file("so3_exp", LINE)?;
        let line: Value = serde_json::from_str(out.trim_end()).map_err(|e| e.to_string())?;
        assert_eq!(line["id"], 0);
        let q: Vec<f64> = line["out"]["q"]
            .as_array()
            .ok_or("no q")?
            .iter()
            .map(|v| hexfloat::parse(v.as_str().unwrap_or_default()))
            .collect::<Result<_, _>>()?;
        assert_eq!(q.len(), 4);
        // The angle is `0x1.921fb54442d18p+0` = 1.5707963267948966, so `q = (√½, 0, 0, √½)`.
        assert!((q[0] - core::f64::consts::FRAC_1_SQRT_2).abs() < 4e-16 && q[1].to_bits() == 0);
        assert_eq!(out.matches('\n').count(), 1);
        Ok(())
    }

    #[test]
    fn a_bad_line_or_an_unsupported_id_is_an_error_naming_the_line() {
        let e = answer_file("so3_exp", &format!("{LINE}\n{{}}\n"))
            .err()
            .unwrap_or_default();
        assert!(e.contains("line 2"), "{e}");
        let bad = LINE.replace("0x1.921fb54442d18p+0", "1.5");
        assert!(answer_file("so3_exp", &bad)
            .err()
            .unwrap_or_default()
            .contains("bad hex float"));
        assert!(answer_file("so3_jr", LINE).is_err());
    }

    /// Every input of every supported corpus file is `float.hex()` to the byte, every record is
    /// answered, and every output is finite. The corpus is committed, so this reads it as it is.
    ///
    /// [`SHAPE`] is the one entry that is not a value: it is checked as the dimensions it is —
    /// integers whose product is the length of the matrix field beside it — so the exception
    /// [`inputs`] makes is asserted here rather than only skipped.
    #[test]
    fn the_whole_supported_corpus_is_answered() -> Result<(), String> {
        for fn_id in supported() {
            let path = format!("{CORPUS}/{fn_id}.jsonl");
            let text = std::fs::read_to_string(&path).map_err(|e| format!("{path}: {e}"))?;
            for line in text.lines() {
                let record: Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
                let object = record["in"].as_object().ok_or("in")?;
                for (key, value) in object.iter().filter(|(k, _)| *k != SHAPE) {
                    for s in value.as_array().ok_or(format!("{key}: array"))? {
                        let s = s.as_str().ok_or(format!("{key}: string"))?;
                        assert_eq!(hexfloat::format(hexfloat::parse(s)?), s);
                    }
                }
                if let Some(shape) = object.get(SHAPE) {
                    let dims: Vec<u64> = shape
                        .as_array()
                        .ok_or("shape: array")?
                        .iter()
                        .map(|d| d.as_u64().ok_or("shape: integer".to_string()))
                        .collect::<Result<_, _>>()?;
                    let entries = dims.iter().product::<u64>() as usize;
                    let of = |v: &Value| v.as_array().map_or(0, Vec::len);
                    assert!(
                        object.values().any(|v| of(v) == entries),
                        "{fn_id}: shape {dims:?} fits no field"
                    );
                }
            }
            let answers = answer_file(fn_id, &text)?;
            assert_eq!(answers.lines().count(), text.lines().count(), "{fn_id}");
            assert!(
                !answers.contains("nan") && !answers.contains("inf"),
                "{fn_id}"
            );
        }
        Ok(())
    }

    #[test]
    fn the_pin_is_the_commit_of_the_manifest_and_of_the_lockfile() -> Result<(), String> {
        let rev = PIN.strip_prefix("tf_tree@").ok_or("PIN")?;
        assert_eq!(rev.len(), 40);
        assert!(include_str!("../Cargo.toml").contains(&format!("rev = \"{rev}\"")));
        let source = format!("git+https://github.com/NoeFontana/tf_tree?rev={rev}#{rev}");
        assert!(include_str!("../Cargo.lock").contains(&source));
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
