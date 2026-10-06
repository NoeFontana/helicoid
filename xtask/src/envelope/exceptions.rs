//! The domination exception table, `conformance/baseline/exceptions.toml`
//! ([`0046`](../../../../docs/decisions/0046-explained-by-record-needs-a-record-to-point-at.md)).
//!
//! `PHASE3.md` §10 and `PHASE4.md` §5.2 both say a stratum an oracle wins blocks the phase "until
//! fixed or **explained by record**", and until this table there was nowhere to record an
//! explanation: the escape hatch the specs presuppose did not exist, so `--bless` could never write
//! a first baseline. A row here excepts **domination and nothing else** (item 2): no-regress, a
//! non-finite output, an unscored row, an oracle over another record count and coverage stay
//! unconditional, and the excepted stratum's `max_u` still goes to the baseline and still has to
//! not regress — so an excepted stratum is watched *more* closely than a dominated one.
//!
//! An exception is a citation, and it cannot outlive the defect it describes: `cargo xtask lint`
//! fails a row whose `record` is missing or is a `draft` (item 3, and `0040` — "a draft authorises
//! nothing"), and a row whose stratum turns out to be dominated after all fails the run and is
//! named for deletion (item 4).
//!
//! The format is read here by hand rather than through a TOML crate: it is a flat array of tables
//! with five string keys and nothing else is accepted, so a strict forty-line reader is both the
//! whole of the grammar and the whole of the validation — an unknown key, a missing field, a
//! repeated field or a duplicate `(fn, stratum, precision)` is an error with a line number.

use std::collections::BTreeMap;
use std::path::Path;

/// One row: what it excepts, the record that explains it, and the reason in one line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Exception {
    pub(crate) fn_id: String,
    pub(crate) stratum: String,
    pub(crate) precision: String,
    /// The record number as written, e.g. `0046`; the lint resolves it to a file and a status.
    pub(crate) record: String,
    pub(crate) reason: String,
    /// The 1-based line of the `[[exception]]` header, for an error that names the row.
    pub(crate) line: usize,
}

impl Exception {
    /// `fn/stratum/precision`, the spelling a `FAIL` line and the evidence page use.
    pub(crate) fn key(&self) -> String {
        format!("{}/{}/{}", self.fn_id, self.stratum, self.precision)
    }
}

/// The table, indexed by `fn/stratum/precision`, in file order.
#[derive(Default)]
pub(crate) struct Exceptions {
    by_key: BTreeMap<String, Exception>,
    order: Vec<String>,
}

/// The file the table lives in, relative to the repository root.
pub(crate) const PATH: &str = "conformance/baseline/exceptions.toml";

const FIELDS: [&str; 5] = ["fn", "stratum", "precision", "record", "reason"];

impl Exceptions {
    /// The table at `root`, or an empty one when the file does not exist: a repository with no
    /// exceptions is the normal state, not an error.
    pub(crate) fn read(root: &Path) -> Result<Self, String> {
        let path = root.join(PATH);
        match std::fs::read_to_string(&path) {
            Ok(text) => Self::parse(&text).map_err(|e| format!("{PATH}: {e}")),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(Self::default()),
            Err(e) => Err(format!("{}: {e}", path.display())),
        }
    }

    /// `[[exception]]` tables of exactly [`FIELDS`], each a basic string on one line.
    pub(crate) fn parse(text: &str) -> Result<Self, String> {
        let mut out = Self::default();
        let mut current: Option<(usize, BTreeMap<&str, String>)> = None;
        for (i, raw) in text.lines().enumerate() {
            let at = i + 1;
            let line = strip_comment(raw).trim();
            if line.is_empty() {
                continue;
            }
            if line == "[[exception]]" {
                if let Some((line_no, fields)) = current.take() {
                    out.push(finish(line_no, fields)?)?;
                }
                current = Some((at, BTreeMap::new()));
                continue;
            }
            let Some((key, value)) = line.split_once('=') else {
                return Err(format!("line {at}: `{line}` is not `key = \"value\"`"));
            };
            let key = key.trim();
            let Some((_, fields)) = current.as_mut() else {
                return Err(format!("line {at}: `{key}` before any `[[exception]]`"));
            };
            if !FIELDS.contains(&key) {
                return Err(format!(
                    "line {at}: unknown key `{key}`; the fields are {}",
                    FIELDS.join(", ")
                ));
            }
            let value = value.trim();
            let text = value
                .strip_prefix('"')
                .and_then(|v| v.strip_suffix('"'))
                .ok_or_else(|| format!("line {at}: `{key}`'s value is not a quoted string"))?;
            if text.contains('"') || text.contains('\\') {
                return Err(format!("line {at}: `{key}`'s value is not a basic string"));
            }
            if fields.insert(key, text.to_string()).is_some() {
                return Err(format!("line {at}: `{key}` is given twice"));
            }
        }
        if let Some((line_no, fields)) = current {
            out.push(finish(line_no, fields)?)?;
        }
        Ok(out)
    }

    fn push(&mut self, e: Exception) -> Result<(), String> {
        let key = e.key();
        if self.by_key.contains_key(&key) {
            return Err(format!("line {}: `{key}` is excepted twice", e.line));
        }
        self.order.push(key.clone());
        self.by_key.insert(key, e);
        Ok(())
    }

    /// The row for `fn/stratum/precision`, if any.
    pub(crate) fn get(&self, key: &str) -> Option<&Exception> {
        self.by_key.get(key)
    }

    /// The rows in file order.
    pub(crate) fn iter(&self) -> impl Iterator<Item = &Exception> {
        self.order.iter().filter_map(|k| self.by_key.get(k))
    }

    #[cfg(test)]
    pub(crate) fn is_empty(&self) -> bool {
        self.order.is_empty()
    }

    /// The file this table would be written as, for a test that round-trips it.
    #[cfg(test)]
    pub(crate) fn render(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::new();
        for e in self.iter() {
            let _ = write!(
                out,
                "[[exception]]\nfn = \"{}\"\nstratum = \"{}\"\nprecision = \"{}\"\n\
                 record = \"{}\"\nreason = \"{}\"\n\n",
                e.fn_id, e.stratum, e.precision, e.record, e.reason
            );
        }
        out
    }
}

/// `raw` without a trailing comment: the first `#` **outside** a quoted span, since a `reason` is
/// free text and may hold one. Splitting at the first `#` would have truncated such a value into a
/// string with no closing quote -- an error rather than silence, but a baffling one.
fn strip_comment(raw: &str) -> &str {
    let mut quoted = false;
    for (i, c) in raw.char_indices() {
        match c {
            '"' => quoted = !quoted,
            '#' if !quoted => return &raw[..i],
            _ => {}
        }
    }
    raw
}

fn finish(line: usize, mut fields: BTreeMap<&str, String>) -> Result<Exception, String> {
    let take = |fields: &mut BTreeMap<&str, String>, k: &str| {
        fields
            .remove(k)
            .ok_or_else(|| format!("line {line}: no `{k}`"))
    };
    let fn_id = take(&mut fields, "fn")?;
    let stratum = take(&mut fields, "stratum")?;
    let precision = take(&mut fields, "precision")?;
    let record = take(&mut fields, "record")?;
    let reason = take(&mut fields, "reason")?;
    if !matches!(precision.as_str(), "f64" | "f32") {
        return Err(format!(
            "line {line}: `precision` is `{precision}`, not `f64` or `f32`"
        ));
    }
    if reason.trim().is_empty() {
        return Err(format!("line {line}: `reason` is empty"));
    }
    Ok(Exception {
        fn_id,
        stratum,
        precision,
        record,
        reason,
        line,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const ONE: &str = "\
[[exception]]
fn = \"so3_exp\"
stratum = \"theta:1e-5\"
precision = \"f64\"
record = \"0046\"
reason = \"ten-digit tie: both programs under 1/2 u, margin 5e-11 relative\"
";

    #[test]
    fn a_row_round_trips_and_is_found_by_its_key() -> Result<(), String> {
        let t = Exceptions::parse(ONE)?;
        assert_eq!(t.iter().count(), 1);
        let e = t.get("so3_exp/theta:1e-5/f64").ok_or("not found")?;
        assert_eq!((e.record.as_str(), e.line), ("0046", 1));
        assert!(e.reason.starts_with("ten-digit tie"));
        assert!(t.get("so3_exp/theta:1e-4/f64").is_none());
        // `render` is the input again, so the page and the file cannot drift apart.
        assert_eq!(Exceptions::parse(&t.render())?.iter().count(), 1);
        Ok(())
    }

    #[test]
    fn an_absent_file_is_an_empty_table_and_a_malformed_one_is_an_error() -> Result<(), String> {
        let t = Exceptions::read(Path::new("/nonexistent-root-for-this-test"))?;
        assert!(t.is_empty() && t.iter().next().is_none());
        let err = |text: &str| Exceptions::parse(text).err().unwrap_or_default();
        assert!(err("fn = \"a\"").contains("before any `[[exception]]`"));
        assert!(err("[[exception]]\nwho = \"a\"").contains("unknown key `who`"));
        assert!(err("[[exception]]\nfn = a").contains("not a quoted string"));
        assert!(err("[[exception]]\nfn = \"a\"").contains("no `stratum`"));
        assert!(err(&ONE.replace("\"f64\"", "\"f16\"")).contains("not `f64` or `f32`"));
        let blank = ONE
            .lines()
            .map(|l| {
                if l.starts_with("reason") {
                    "reason = \"  \""
                } else {
                    l
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        assert!(err(&blank).contains("is empty"));
        // A `#` inside the value is part of it, not the start of a comment.
        let hashed = ONE.replace("ten-digit tie", "ten-digit tie #1");
        assert!(Exceptions::parse(&hashed)?
            .iter()
            .any(|e| e.reason.contains("#1")));
        assert!(err(&format!("{ONE}{ONE}")).contains("is excepted twice"));
        assert!(err("[[exception]]\nfn = \"a\"\nfn = \"b\"").contains("`fn` is given twice"));
        assert!(err("[[exception]]\nfn\n").contains("is not `key = \"value\"`"));
        Ok(())
    }

    /// A comment and blank lines are skipped, and a trailing comment on a value is not part of it.
    #[test]
    fn comments_and_blank_lines_are_skipped() -> Result<(), String> {
        let text = format!("# the table\n\n{ONE}\n# trailing\n");
        let t = Exceptions::parse(&text)?;
        assert_eq!(t.iter().count(), 1);
        Ok(())
    }
}
