//! Check 2: a `draft` decision record authorises nothing (`docs/decisions/README.md`, Lifecycle).
//! Where a draft is cited in a spec status table (`## 0.0` or `## 0.`), a Rust comment under
//! `crates/` or `xtask/`, or an amendment banner (a `>` block that says "amend"), the word
//! `draft` must sit within three words of the citation. A line that names `tf_tree` cites that
//! project's own record numbers and is skipped. No other context is inspected: prose elsewhere
//! and `docs/PROJECT.md` §5.1 are out of scope for now.

use super::{comments, File, Violation};

/// `Some(number)` for `docs/decisions/NNNN-slug.md`.
fn record_number(path: &str) -> Option<&str> {
    let name = path.strip_prefix("docs/decisions/")?.strip_suffix(".md")?;
    let number = name
        .get(..4)
        .filter(|n| n.bytes().all(|c| c.is_ascii_digit()))?;
    name[4..].starts_with('-').then_some(number)
}

/// First word of the value of the record's `**Status:**` line, lowercased.
fn status(text: &str) -> Option<String> {
    let value = text.lines().find_map(|l| l.strip_prefix("**Status:**"))?;
    let word = value.split_whitespace().next()?;
    Some(
        word.trim_matches(|c: char| !c.is_alphanumeric())
            .to_lowercase(),
    )
}

/// Byte ranges in `line` where `number` appears as a token (not part of a longer alphanumeric
/// run, a version or a decimal); a link to `NNNN-slug.md` counts through its number.
fn citations(line: &str, number: &str) -> Vec<(usize, usize)> {
    let b = line.as_bytes();
    let glued = |c: u8| c.is_ascii_alphanumeric() || c == b'_';
    line.match_indices(number)
        .filter(|&(at, _)| {
            let before = at.checked_sub(1).map(|i| b[i]);
            let after = b.get(at + 4).copied();
            !before.is_some_and(|c| glued(c) || c == b'.') && !after.is_some_and(glued)
        })
        .map(|(at, _)| (at, at + 4))
        .collect()
}

/// Whether the word `draft` is within three words of `line[s..e]` on either side, without
/// crossing another four-digit record number.
fn says_draft(line: &str, (s, e): (usize, usize)) -> bool {
    let word = |w: &str| {
        w.trim_matches(|c: char| !c.is_alphanumeric())
            .to_lowercase()
    };
    let is_number = |w: &str| w.len() == 4 && w.bytes().all(|c| c.is_ascii_digit());
    let near = |words: &mut dyn Iterator<Item = &str>| {
        words
            .map(word)
            .take(3)
            .take_while(|w| !is_number(w))
            .any(|w| w == "draft")
    };
    near(&mut line[e..].split_whitespace()) || near(&mut line[..s].split_whitespace().rev())
}

pub(crate) fn check(files: &[File]) -> Vec<Violation> {
    let mut out = Vec::new();
    let mut drafts = Vec::new();
    for f in files {
        let Some(number) = record_number(&f.path) else {
            continue;
        };
        match status(&f.text).as_deref() {
            Some("draft") => drafts.push(number),
            Some(_) => {}
            None => out.push(Violation::new(
                &f.path,
                1,
                "[drafts] record has no `**Status:**` line",
            )),
        }
    }
    if drafts.is_empty() {
        return out;
    }
    for f in files {
        for (line, text) in contexts(f) {
            if text.contains("tf_tree") {
                continue;
            }
            for number in &drafts {
                if citations(text, number)
                    .into_iter()
                    .any(|c| !says_draft(text, c))
                {
                    out.push(Violation::new(
                        &f.path,
                        line,
                        format!(
                            "[drafts] cites draft record {} without saying `draft` next to it",
                            number
                        ),
                    ));
                }
            }
        }
    }
    out
}

/// The `(line, text)` pairs of `f` where citing a draft is inspected.
fn contexts(f: &File) -> Vec<(usize, &str)> {
    let in_code =
        (f.path.starts_with("crates/") || f.path.starts_with("xtask/")) && f.path.ends_with(".rs");
    if in_code {
        return comments::of(&f.path, &f.text)
            .iter()
            .map(|c| (c.line, c.text))
            .collect();
    }
    if !f.path.ends_with(".md") {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut status_level: Option<usize> = None;
    let (mut fenced, mut block, mut amends) = (false, Vec::new(), false);
    for (i, l) in f.text.lines().chain(std::iter::once("")).enumerate() {
        let t = l.trim_start();
        // A banner is the whole contiguous `>` block that says "amend" anywhere.
        if t.starts_with('>') && !fenced {
            block.push((i + 1, l));
            amends |= t.to_lowercase().contains("amend");
            continue;
        }
        if amends {
            out.append(&mut block);
        }
        (block, amends) = (Vec::new(), false);
        if t.starts_with("```") || t.starts_with("~~~") {
            fenced = !fenced;
            continue;
        }
        if fenced {
            continue;
        }
        let level = t.bytes().take_while(|&c| c == b'#').count();
        if (1..=6).contains(&level) && t[level..].starts_with(' ') {
            if status_level.is_some_and(|s| level <= s) {
                status_level = None;
            }
            if status_level.is_none() && is_status_heading(t[level..].trim()) {
                status_level = Some(level);
            }
        } else if status_level.is_some() && t.starts_with('|') {
            out.push((i + 1, l));
        }
    }
    out
}

/// `0.0 ...` or `0. ...` (a spec's status section), as a whole leading token.
fn is_status_heading(text: &str) -> bool {
    let rest = text
        .strip_prefix("0.")
        .map(|r| r.strip_prefix('0').unwrap_or(r));
    rest.is_some_and(|r| r.is_empty() || r.starts_with(char::is_whitespace))
}

#[cfg(test)]
mod tests {
    use super::*;

    const D13: &str = "docs/decisions/0013-simd.md";

    fn tree(extra: &[(&str, &str)]) -> Vec<File> {
        let mut v = vec![
            File::new(D13, "# 0013: x\n\n**Status:** draft\n"),
            File::new("docs/decisions/0007-b.md", "# 0007\n\n**Status:** ready\n"),
        ];
        v.extend(extra.iter().map(|(p, t)| File::new(p, t)));
        v
    }

    fn lines(extra: &[(&str, &str)]) -> Vec<String> {
        check(&tree(extra))
            .iter()
            .map(ToString::to_string)
            .collect()
    }

    #[test]
    fn status_table_rows() {
        let spec = "## 0.0 Implementation status\n\n| Area | Status |\n|---|---|\n| SIMD | Done per 0013 |\n| SIMD | see 0013 (draft) |\n| B | 0007 |\n\n## 1. Next\n\n| x | 0013 |\n";
        let v = lines(&[("docs/PHASE9.md", spec)]);
        assert_eq!(v.len(), 1, "{v:?}");
        assert!(v[0].starts_with("docs/PHASE9.md:5: "), "{v:?}");
    }

    #[test]
    fn link_and_file_name_forms() {
        let row = "## 0.0 S\n| a | [x](./decisions/0013-simd.md) |\n";
        assert_eq!(lines(&[("docs/P.md", row)]).len(), 1);
    }

    #[test]
    fn rust_comments_only() {
        let src = "// settled by 0013\nlet s = \"0013\"; // 0013 (draft)\n/// per 0007\n";
        let v = lines(&[("crates/a/src/lib.rs", src)]);
        assert_eq!(v.len(), 1, "{v:?}");
        assert!(v[0].starts_with("crates/a/src/lib.rs:1: "), "{v:?}");
        assert!(lines(&[("runners/a/src/lib.rs", "// 0013\n")]).is_empty());
        assert_eq!(lines(&[("xtask/src/x.rs", "// 0013\n")]).len(), 1);
    }

    #[test]
    fn amendment_banners_are_whole_blocks() {
        let bad = "> **Amended** by 0013.\n> Companions: 0013\n\nprose 0013\n> quote 0013\n";
        let v = lines(&[(
            "docs/decisions/0002-c.md",
            &format!("**Status:** ready\n{bad}"),
        )]);
        assert_eq!(v.len(), 2, "{v:?}");
        assert!(v[0].contains(":2: ") && v[1].contains(":3: "), "{v:?}");
        assert!(lines(&[("docs/A.md", "> amended by 0013, a draft\n")]).is_empty());
        // A wrapped banner whose citation sits on the continuation line, at end of file too.
        let wrapped = "> **Amended** by the following record:\n> 0013 settles it.";
        assert_eq!(lines(&[("docs/A.md", wrapped)]).len(), 1);
    }

    #[test]
    fn status_scope_ignores_fences_and_sub_headings() {
        let spec = "## 0.0 Status\n| a | 0013 |\n### 0.0.1 sub\n| b | 0013 |\n### Other\n| c | 0013 |\n```sh\n# comment\n| f | 0013 |\n```\n| d | 0013 |\n\n## 1. Next\n| e | 0013 |\n";
        let at: Vec<String> = lines(&[("docs/E.md", spec)])
            .iter()
            .map(|v| v.split(": ").next().unwrap_or_default().to_string())
            .collect();
        assert_eq!(
            at,
            ["docs/E.md:2", "docs/E.md:4", "docs/E.md:6", "docs/E.md:11"]
        );
        for other in [
            "## 10.0 X\n| a | 0013 |\n",
            "## 1.0.0 Release\n| a | 0013 |\n",
            "# T\n```sh\n# bump to 0.0.1\n```\n## 2. Other\n| a | 0013 |\n",
        ] {
            assert!(lines(&[("docs/P.md", other)]).is_empty(), "{other}");
        }
        assert_eq!(
            lines(&[("docs/N.md", "## 0. Status\n| a | 0013 |\n")]).len(),
            1
        );
    }

    #[test]
    fn draft_must_be_a_word_next_to_the_citation() {
        let row = |cell: &str| format!("## 0.0 S\n| a | {cell} |\n");
        for ok in [
            "0013 (draft)",
            "Draft 0013",
            "the draft record 0013",
            "0013 is a draft",
        ] {
            assert!(lines(&[("docs/P.md", &row(ok))]).is_empty(), "{ok}");
        }
        for bad in [
            "0013 is the redrafted plan",
            "0013 settled, 0014 (draft)",
            "0013 settled the SIMD question in a draft",
        ] {
            assert_eq!(lines(&[("docs/P.md", &row(bad))]).len(), 1, "{bad}");
        }
    }

    #[test]
    fn number_is_a_token_and_other_contexts_are_ignored() {
        let spec = "## 0.0 S\n| a | 00131 and 20013 and 0.0013 and 0x0013 and id_0013 |\n";
        assert!(lines(&[("docs/P.md", spec)]).is_empty());
        assert!(lines(&[("docs/P.md", "prose citing 0013 freely\n| a | 0013 |\n")]).is_empty());
        assert!(lines(&[("crates/a/src/lib.rs", "// mask 0x0013_ffff\n")]).is_empty());
        assert!(lines(&[("docs/P.md", "## 0.0 S\n| a | tf_tree 0013 rows |\n")]).is_empty());
    }

    #[test]
    fn only_drafts_are_policed_and_status_is_required() {
        let spec = "## 0.0 S\n| a | 0007 |\n";
        assert!(lines(&[("docs/P.md", spec)]).is_empty());
        let files = [File::new("docs/decisions/0099-x.md", "# no status\n")];
        assert_eq!(check(&files).len(), 1);
        assert!(check(&[File::new("docs/decisions/template.md", "# t\n")]).is_empty());
    }

    #[test]
    fn non_draft_statuses_do_not_trigger() {
        let files = [
            File::new(
                "docs/decisions/0050-x.md",
                "**Status:** superseded by 0051\n",
            ),
            File::new("docs/P.md", "## 0.0 S\n| a | 0050 |\n"),
        ];
        assert!(check(&files).is_empty());
    }
}
