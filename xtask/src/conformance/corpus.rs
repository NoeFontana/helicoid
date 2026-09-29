//! The committed corpus (`conformance/generate/README.md`, "Record schema"): one JSONL file per
//! function id, listed with its `kind` in `MANIFEST.json`. Only `kind: "corpus"` files are function
//! ids; `coeff_series` is skipped (`docs/PHASE1.md` §4.3). A `.jsonl` the manifest does not list,
//! a `kind` other than `corpus` and `series`, a repeated key, a wrong `id` and a wrong record count
//! are errors. The manifest's `sha256` is not read: `just corpus-check` compares bytes.

use std::collections::BTreeMap;
use std::fmt;
use std::path::Path;

use serde::de::{Deserializer, Error as _, MapAccess, Visitor};
use serde::Deserialize;
use serde_json::Value;

use super::number::{parse_hex, Decimal};

/// One named value of a record: a scalar, a vector, or a column-major matrix (`shape` is
/// `Some((rows, cols))`).
#[derive(Debug)]
pub(crate) struct Tensor<T> {
    pub(crate) shape: Option<(usize, usize)>,
    pub(crate) data: Vec<T>,
}

type Fields<T> = BTreeMap<String, Tensor<T>>;

/// One corpus line: exact binary64 inputs and the exact decimal reference.
#[derive(Debug)]
pub(crate) struct Record {
    pub(crate) id: u64,
    pub(crate) stratum: String,
    pub(crate) inputs: Fields<f64>,
    pub(crate) reference: Fields<Decimal>,
}

impl Record {
    /// The flat, column-major values of input `key` (a scalar is one value).
    pub(crate) fn input(&self, key: &str) -> Option<&[f64]> {
        self.inputs.get(key).map(|t| t.data.as_slice())
    }

    /// An `@f32` stratum (`docs/decisions/0016`): exact binary32 inputs, scored at `f32` only.
    pub(crate) fn is_f32_stratum(&self) -> bool {
        self.stratum.ends_with(F32_SUFFIX)
    }

    /// Every input is exactly a binary32, so an `f32` subject receives it by a lossless cast.
    pub(crate) fn require_binary32(&self) -> Result<(), String> {
        for (key, t) in &self.inputs {
            if let Some(x) = t.data.iter().find(|&&x| exact_f32(x).is_none()) {
                return Err(format!("input `{key}` holds {x:e}, not a binary32"));
            }
        }
        Ok(())
    }
}

/// The suffix of the `f32`-exact strata's names (`docs/decisions/0016`).
pub(crate) const F32_SUFFIX: &str = "@f32";

/// `x` as a binary32, when that is lossless: the cast an `f32` subject makes of an `@f32` input.
pub(crate) fn exact_f32(x: f64) -> Option<f32> {
    let y = x as f32;
    (f64::from(y).to_bits() == x.to_bits()).then_some(y)
}

/// A JSON object that refuses a repeated key (`BTreeMap`'s own `Deserialize` keeps the last).
pub(crate) struct Object(pub(crate) BTreeMap<String, Value>);

impl<'de> Deserialize<'de> for Object {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct Unique;
        impl<'de> Visitor<'de> for Unique {
            type Value = Object;

            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("an object")
            }

            fn visit_map<A: MapAccess<'de>>(self, mut entries: A) -> Result<Object, A::Error> {
                let mut map = BTreeMap::new();
                while let Some((key, value)) = entries.next_entry::<String, Value>()? {
                    if map.contains_key(&key) {
                        return Err(A::Error::custom(format!("duplicate key `{key}`")));
                    }
                    map.insert(key, value);
                }
                Ok(Object(map))
            }
        }
        d.deserialize_map(Unique)
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Line {
    id: u64,
    stratum: String,
    #[serde(rename = "in")]
    input: Object,
    out: Object,
}

/// Splits `shape` off `object` and reads every other entry with `elem`. A `shape` sits beside
/// exactly one array, the matrix's flat data.
fn fields<T>(
    mut object: BTreeMap<String, Value>,
    elem: impl Fn(&str) -> Result<T, String>,
) -> Result<Fields<T>, String> {
    let shape = match object.remove("shape") {
        None => None,
        Some(v) => match serde_json::from_value::<[usize; 2]>(v) {
            Ok([r, c]) => Some((r, c)),
            Err(e) => return Err(format!("bad `shape`: {e}")),
        },
    };
    fn text(v: &Value) -> Result<&str, String> {
        v.as_str().ok_or_else(|| format!("not a string: {v}"))
    }
    let mut out = Fields::new();
    for (key, value) in object {
        let data = match &value {
            Value::Array(items) => items.iter().map(|v| elem(text(v)?)).collect(),
            v => elem(text(v)?).map(|x| vec![x]),
        }
        .map_err(|e| format!("`{key}`: {e}"))?;
        let shape = shape.filter(|_| value.is_array());
        out.insert(key, Tensor { shape, data });
    }
    if let Some((r, c)) = shape {
        let mut arrays = out.iter().filter(|(_, t)| t.shape.is_some());
        match (arrays.next(), arrays.next()) {
            (Some((_, t)), None) if t.data.len() == r * c => {}
            (Some((k, t)), None) => {
                return Err(format!(
                    "`{k}` has {} values for shape {r}x{c}",
                    t.data.len()
                ))
            }
            _ => return Err("`shape` needs exactly one array beside it".into()),
        }
    }
    Ok(out)
}

/// Parses one line of a function file.
pub(crate) fn parse_line(line: &str) -> Result<Record, String> {
    let l: Line = serde_json::from_str(line).map_err(|e| e.to_string())?;
    Ok(Record {
        id: l.id,
        stratum: l.stratum,
        inputs: fields(l.input.0, parse_hex).map_err(|e| format!("in: {e}"))?,
        reference: fields(l.out.0, |s| s.parse::<Decimal>()).map_err(|e| format!("out: {e}"))?,
    })
}

/// A function file the manifest lists as `kind: "corpus"`.
pub(crate) struct Entry {
    pub(crate) fn_id: String,
    pub(crate) records: usize,
}

#[derive(Deserialize)]
struct ManifestFile {
    kind: String,
    records: usize,
}

#[derive(Deserialize)]
struct Manifest {
    files: BTreeMap<String, ManifestFile>,
}

/// The function ids of `dir/MANIFEST.json`, sorted by file name. Every `.jsonl` in `dir` must be
/// listed, with a `kind` of `corpus` or `series`.
pub(crate) fn manifest(dir: &Path) -> Result<Vec<Entry>, String> {
    let path = dir.join("MANIFEST.json");
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let m: Manifest =
        serde_json::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))?;
    let unlisted = |e: std::io::Result<std::fs::DirEntry>| {
        let name = e.ok()?.file_name().into_string().ok()?;
        (name.ends_with(".jsonl") && !m.files.contains_key(&name)).then_some(name)
    };
    let dir_err = |e: std::io::Error| format!("{}: {e}", dir.display());
    let mut extra: Vec<String> = std::fs::read_dir(dir)
        .map_err(dir_err)?
        .filter_map(unlisted)
        .collect();
    extra.sort();
    if let Some(name) = extra.first() {
        return Err(format!("{}: `{name}` is not listed", path.display()));
    }
    let mut entries = Vec::new();
    for (name, f) in m.files {
        match f.kind.as_str() {
            "corpus" => entries.push(Entry {
                fn_id: name.strip_suffix(".jsonl").unwrap_or(&name).to_string(),
                records: f.records,
            }),
            "series" => {}
            kind => return Err(format!("{}: `{name}` has kind `{kind}`", path.display())),
        }
    }
    Ok(entries)
}

/// Whether `entry` may have `@f32` strata: no false negative (a stratum name ends the way no input
/// or reference does), and a false positive is caught where the file is read.
pub(crate) fn mentions_f32(dir: &Path, entry: &Entry) -> Result<bool, String> {
    let path = dir.join(format!("{}.jsonl", entry.fn_id));
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(text.contains(&format!("{F32_SUFFIX}\"")))
}

/// Every record of `entry`, in file order; the count must be the manifest's and each `id` its
/// 0-based line number.
pub(crate) fn read(dir: &Path, entry: &Entry) -> Result<Vec<Record>, String> {
    let path = dir.join(format!("{}.jsonl", entry.fn_id));
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut records = Vec::with_capacity(entry.records);
    for (i, line) in text.lines().enumerate() {
        let at = |e: String| format!("{}:{}: {e}", path.display(), i + 1);
        let r = parse_line(line).map_err(at)?;
        if r.id != i as u64 {
            return Err(at(format!("id {} on line {}", r.id, i + 1)));
        }
        records.push(r);
    }
    if records.len() != entry.records {
        return Err(format!(
            "{}: {} records, MANIFEST.json says {}",
            path.display(),
            records.len(),
            entry.records
        ));
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;
    use crate::conformance::testkit::Scratch;

    const VECTOR: &str = r#"{"id":3,"in":{"phi":["0x1.8p+1","-0x0.0p+0","0x1.0000000000001p-1022"]},"out":{"q":["9.99999999999999999999999251527e-1","-1.06240087782203108924896608142e-13"]},"stratum":"theta:1e-12"}"#;
    const MATRIX: &str = r#"{"id":0,"in":{"R":["0x1p+0","0x0.0p+0","0x0.0p+0","0x0.0p+0","0x1p+0","0x0.0p+0"],"shape":[3,2]},"out":{"J":["1e0","0e0","0e0","0e0","1e0","2.5e-311"],"shape":[3,2],"s":"1e0"},"stratum":"s"}"#;

    #[test]
    fn vectors_scalars_and_matrices_parse_exactly() -> Result<(), String> {
        let r = parse_line(VECTOR)?;
        assert_eq!((r.id, r.stratum.as_str()), (3, "theta:1e-12"));
        assert_eq!(r.input("phi").map(<[f64]>::len), Some(3));
        let phi = r.input("phi").ok_or("phi")?;
        assert_eq!((phi[0], phi[1].to_bits()), (3.0, (-0.0f64).to_bits()));
        assert_eq!(phi[2].to_bits(), (1 << 52) + 1);
        let q = &r.reference["q"];
        assert_eq!(q.shape, None);
        assert!(q.data[1].neg && q.data[1].exp10 == -13 - 29);

        let m = parse_line(MATRIX)?;
        assert_eq!(m.inputs["R"].shape, Some((3, 2)));
        assert_eq!(m.reference["J"].shape, Some((3, 2)));
        assert_eq!(m.reference["J"].data[5].exp10, -311 - 1);
        assert_eq!(m.reference["s"].shape, None);
        assert_eq!(m.reference["s"].data.len(), 1);
        Ok(())
    }

    #[test]
    fn schema_violations_are_named() {
        let bad = |line: &str, needle: &str| {
            let e = parse_line(line).err().unwrap_or_default();
            assert!(e.contains(needle), "{line}: {e}");
        };
        bad(&VECTOR.replace("0x1.8p+1", "3.0"), "bad hex float");
        bad(
            &VECTOR.replace("9.99999999999999999999999251527e-1", "x"),
            "bad decimal",
        );
        bad(&MATRIX.replace("[3,2]", "[3,3]"), "6 values for shape 3x3");
        bad(
            &MATRIX.replace(r#""s":"1e0""#, r#""s":["1e0"]"#),
            "exactly one array",
        );
        bad(&MATRIX.replace("[3,2],\"s\"", "[3],\"s\""), "bad `shape`");
        bad(&VECTOR.replace("\"id\":3,", ""), "missing field");
        bad(
            &VECTOR.replace("\"stratum\"", "\"extra\":1,\"stratum\""),
            "unknown field",
        );
        bad(&VECTOR.replace("\"-0x0.0p+0\"", "1"), "not a string");
    }

    #[test]
    fn the_committed_manifest_lists_every_function_file_and_it_reads() -> Result<(), String> {
        let dir = super::super::corpus_dir()?;
        let entries = manifest(&dir)?;
        assert!(entries.iter().all(|e| e.fn_id != "coeff_series"));
        assert!(entries.len() >= 45, "{} ids", entries.len());
        let small = entries
            .iter()
            .find(|e| e.fn_id == "so2_exp")
            .ok_or("so2_exp")?;
        let records = read(&dir, small)?;
        assert_eq!(records.len(), 3419);
        assert_eq!(records[0].stratum, "theta:1e-12");
        Ok(())
    }

    /// `docs/decisions/0016`: an `@f32` stratum holds its binary64 twin's inputs rounded to
    /// binary32 (Rust's cast rounds to nearest even), so an `f32` subject receives them by a
    /// lossless cast. `theta:subnormal@f32` has a decade of its own, below the smallest normal.
    #[test]
    fn the_f32_strata_are_the_binary64_ones_rounded_to_binary32() -> Result<(), String> {
        let round = |x: f64| f64::from(x as f32).to_bits();
        let dir = super::super::corpus_dir()?;
        let mut twins = 0;
        for e in manifest(&dir)?
            .iter()
            .filter(|e| e.fn_id.starts_with("coeff_"))
        {
            let records = read(&dir, e)?;
            let of = |name: &str| -> Vec<&Record> {
                records.iter().filter(|r| r.stratum == name).collect()
            };
            let names: BTreeSet<&str> = records.iter().map(|r| r.stratum.as_str()).collect();
            for name in names.iter().filter_map(|s| s.strip_suffix("@f32")) {
                let twin = format!("{name}@f32");
                assert_eq!(of(name).len(), of(&twin).len(), "{} {name}", e.fn_id);
                for (w, n) in of(name).into_iter().zip(of(&twin)) {
                    for (key, t) in &n.inputs {
                        for (a, b) in w.inputs[key].data.iter().zip(&t.data) {
                            assert_eq!(round(*b), b.to_bits(), "{} {twin} {key}", e.fn_id);
                            match (name == "theta:subnormal", key.as_str()) {
                                (false, _) => assert_eq!(round(*a), b.to_bits(), "{twin}"),
                                (true, "w") => assert_eq!(b.to_bits(), 1f64.to_bits()),
                                (true, _) => assert!(*b < f64::from(f32::MIN_POSITIVE), "{twin}"),
                            }
                        }
                    }
                }
                twins += 1;
            }
        }
        assert_eq!(
            twins,
            8 * 28 + 1,
            "28 theta strata in each of 8 ids, and `q:w0@f32`"
        );
        Ok(())
    }

    /// Every value of every input is a binary32, not the first of each: a vector input and
    /// `coeff_r`'s two keys reach an `f32` subject through the same cast.
    #[test]
    fn every_value_of_every_input_must_be_a_binary32() -> Result<(), String> {
        use crate::conformance::testkit::record;
        let (bad, good) = (0.1, 0.5);
        let named = |r: Record, want: &str| {
            let e = r.require_binary32().err().unwrap_or_default();
            assert!(e.starts_with(want), "{e}");
        };
        assert_eq!(
            record(&[("n", &[good, 0.25]), ("w", &[1.0, 2.0])], &[])?.require_binary32(),
            Ok(())
        );
        named(
            record(&[("phi", &[good, bad])], &[])?,
            "input `phi` holds 1e-1",
        );
        named(
            record(&[("n", &[good]), ("w", &[bad])], &[])?,
            "input `w` holds 1e-1",
        );
        named(
            record(&[("n", &[good, bad]), ("w", &[good])], &[])?,
            "input `n` holds 1e-1",
        );
        Ok(())
    }

    /// A corpus directory holding `MANIFEST.json` and `files`.
    fn corpus_dir(tag: &str, manifest: &str, files: &[(&str, &str)]) -> Result<Scratch, String> {
        let scratch = Scratch::new(tag);
        let io = |e: std::io::Error| e.to_string();
        std::fs::create_dir_all(&scratch.0).map_err(io)?;
        std::fs::write(scratch.0.join("MANIFEST.json"), manifest).map_err(io)?;
        for (name, text) in files {
            std::fs::write(scratch.0.join(name), text).map_err(io)?;
        }
        Ok(scratch)
    }

    fn line(id: u64) -> String {
        VECTOR.replace("\"id\":3", &format!("\"id\":{id}"))
    }

    fn manifest_of(records: usize) -> String {
        format!(r#"{{"files":{{"f.jsonl":{{"kind":"corpus","records":{records}}}}}}}"#)
    }

    /// The error of reading `f.jsonl` as `text` under a manifest that says `records`.
    fn read_error(tag: &str, records: usize, text: &str) -> Result<String, String> {
        let dir = corpus_dir(tag, &manifest_of(records), &[("f.jsonl", text)])?;
        let entries = manifest(&dir.0)?;
        Ok(read(&dir.0, &entries[0]).err().unwrap_or_default())
    }

    #[test]
    fn a_wrong_id_a_wrong_count_a_blank_line_and_a_repeated_key_are_errors() -> Result<(), String> {
        let good = format!("{}\n{}\n", line(0), line(1));
        let dir = corpus_dir("ok", &manifest_of(2), &[("f.jsonl", &good)])?;
        let entries = manifest(&dir.0)?;
        assert_eq!(read(&dir.0, &entries[0])?.len(), 2);
        let skipped = format!("{}\n{}\n", line(0), line(2));
        assert!(read_error("id", 2, &skipped)?.contains("id 2 on line 2"));
        assert!(read_error("count", 3, &good)?.contains("2 records, MANIFEST.json says 3"));
        let blank = format!("{}\n\n{}\n", line(0), line(1));
        assert!(read_error("blank", 2, &blank)?.contains("f.jsonl:2:"));
        let twice = VECTOR.replace("\"phi\":[", "\"phi\":[\"0x1p+0\"],\"phi\":[");
        let twice = twice.replace("\"id\":3", "\"id\":0");
        assert!(read_error("twice", 1, &twice)?.contains("duplicate key `phi`"));
        Ok(())
    }

    #[test]
    fn the_manifest_lists_every_jsonl_and_knows_every_kind() -> Result<(), String> {
        let empty = ("f.jsonl", "");
        let series = r#"{"files":{"f.jsonl":{"kind":"corpus","records":0},"s.jsonl":{"kind":"series","records":0}}}"#;
        let dir = corpus_dir("series", series, &[empty, ("s.jsonl", "")])?;
        assert_eq!(manifest(&dir.0)?.len(), 1);
        let dir = corpus_dir("extra", &manifest_of(0), &[empty, ("g.jsonl", "")])?;
        let e = manifest(&dir.0).err().unwrap_or_default();
        assert!(e.contains("`g.jsonl` is not listed"), "{e}");
        let odd = r#"{"files":{"f.jsonl":{"kind":"corps","records":0}}}"#;
        let dir = corpus_dir("kind", odd, &[empty])?;
        let e = manifest(&dir.0).err().unwrap_or_default();
        assert!(e.contains("kind `corps`"), "{e}");
        Ok(())
    }
}
