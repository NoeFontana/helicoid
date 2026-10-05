//! An exception is a citation, and this is the check
//! ([`0046`](../../../docs/decisions/0046-explained-by-record-needs-a-record-to-point-at.md)
//! item 3).
//!
//! Every row of `conformance/baseline/exceptions.toml` names the record that explains it. The row
//! is only as good as that citation, so `cargo xtask lint` fails one whose record does not exist
//! or is not `ready`: a `draft` cannot except a bar, because
//! [`0040`](../../../docs/decisions/0040-a-draft-is-not-a-parking-space.md) says a draft
//! authorises nothing. The status reader is `drafts`', so the `**Status:**` convention is spelled
//! once.
//!
//! A malformed table is one violation here rather than a panic in the envelope: `just lint` runs
//! on every change and `just envelope` does not, so this is where a typo should surface.

use super::{drafts, File, Violation};
use crate::envelope::exceptions::Exceptions;

pub(crate) fn check(files: &[File]) -> Vec<Violation> {
    let path = crate::envelope::exceptions::PATH;
    let Some(table) = files.iter().find(|f| f.path == path) else {
        return Vec::new();
    };
    let say = |line: usize, message: String| Violation {
        path: path.to_string(),
        line,
        message: format!("[exceptions] {message}"),
    };
    let parsed = match Exceptions::parse(&table.text) {
        Ok(t) => t,
        Err(e) => return vec![say(1, e)],
    };
    let mut out = Vec::new();
    for e in parsed.iter() {
        let prefix = format!("docs/decisions/{}-", e.record);
        let record = files.iter().find(|f| f.path.starts_with(&prefix));
        match record.and_then(|f| drafts::status(&f.text)) {
            Some(s) if s == "ready" => {}
            Some(s) => out.push(say(
                e.line,
                format!(
                    "`{}` cites `{}`, whose status is `{s}`: only a `ready` record may except a \
                     bar (`0040`)",
                    e.key(),
                    e.record
                ),
            )),
            None if record.is_some() => out.push(say(
                e.line,
                format!(
                    "`{}` cites `{}`, which has no `**Status:**` line",
                    e.key(),
                    e.record
                ),
            )),
            None => out.push(say(
                e.line,
                format!(
                    "`{}` cites `{}`, and no `docs/decisions/{}-*.md` exists",
                    e.key(),
                    e.record,
                    e.record
                ),
            )),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(record: &str) -> String {
        format!(
            "[[exception]]\nfn = \"so3_exp\"\nstratum = \"theta:1e-5\"\nprecision = \"f64\"\n\
             record = \"{record}\"\nreason = \"a tie\"\n"
        )
    }

    fn files(record: &str, status: &str) -> Vec<File> {
        let mut out = vec![File::new(crate::envelope::exceptions::PATH, &row(record))];
        if !status.is_empty() {
            out.push(File::new(
                "docs/decisions/0046-x.md",
                &format!("# 0046\n\n**Status:** {status}\n"),
            ));
        }
        out
    }

    fn messages(files: &[File]) -> Vec<String> {
        check(files).iter().map(ToString::to_string).collect()
    }

    #[test]
    fn a_ready_record_passes_and_a_draft_a_missing_one_and_a_typo_do_not() {
        assert!(messages(&files("0046", "ready")).is_empty());
        let draft = messages(&files("0046", "draft"));
        assert_eq!(draft.len(), 1);
        assert!(draft[0].contains("whose status is `draft`"), "{draft:?}");
        let gone = messages(&files("0046", ""));
        assert_eq!(gone.len(), 1);
        assert!(
            gone[0].contains("no `docs/decisions/0046-*.md` exists"),
            "{gone:?}"
        );
        // A record whose file exists without a status line is its own message.
        let blank = vec![
            File::new(crate::envelope::exceptions::PATH, &row("0046")),
            File::new("docs/decisions/0046-x.md", "# 0046\n"),
        ];
        assert!(messages(&blank)[0].contains("no `**Status:**` line"));
        // A malformed table is one violation at line 1, not a panic later.
        let bad = vec![File::new(
            crate::envelope::exceptions::PATH,
            "[[exception]]\nwho = \"a\"\n",
        )];
        let m = messages(&bad);
        assert_eq!(m.len(), 1);
        assert!(m[0].contains("unknown key `who`"), "{m:?}");
        // No table is no violation: a checkout with nothing excepted is the normal state.
        assert!(check(&[]).is_empty());
    }

    /// The row's line number is reported, so a table of many says which one.
    #[test]
    fn the_violation_names_the_rows_line() {
        let two = format!("{}\n{}", row("0046"), row("0099").replace("1e-5", "1e-4"));
        let files = vec![File::new(crate::envelope::exceptions::PATH, &two)];
        let m = messages(&files);
        assert_eq!(m.len(), 2);
        assert!(m[0].starts_with(&format!("{}:1:", crate::envelope::exceptions::PATH)));
        assert!(m[1].contains(":8:"), "{m:?}");
    }
}
