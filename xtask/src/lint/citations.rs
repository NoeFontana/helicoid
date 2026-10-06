//! Check 1: no `<name>.rs:<digits>` or `<name>.rs#L<digits>` line citation in markdown or in the
//! comments of code and config files (`docs/PHASE1.md` §3). Cite a symbol instead.
//!
//! Check 2: a cited symbol **resolves**. "Cite a symbol, never a line number" is only half a rule
//! while a renamed symbol passes silently -- the whole point of the convention is that a citation
//! a rename breaks says so, and a stale one is worse than a line number because `grep` returns
//! nothing and the reader cannot tell whether the claim moved or was deleted. So a
//! `<module>_tests::<name>` citation is resolved against that module's file: the module has to be
//! one of the tracked files, and `name` has to be defined in it. Two such citations shipped stale
//! in the PR that added this check (`0048`).

use super::{comments, File, Violation};

pub(crate) fn check(files: &[File]) -> Vec<Violation> {
    let mut out = Vec::new();
    for f in files {
        for c in comments::of(&f.path, &f.text) {
            if let Some(cite) = find(c.text) {
                out.push(Violation::new(
                    &f.path,
                    c.line,
                    format!(
                        "[citations] line citation `{cite}`; cite a symbol, never a line number"
                    ),
                ));
            }
        }
    }
    out
}

/// Check 2: every `<module>_tests::<name>` citation resolves to a definition in that module.
///
/// A module this tree does not hold is skipped, not failed: a comment may name another crate's
/// test module, and the lint reads only tracked files.
pub(crate) fn check_symbols(files: &[File]) -> Vec<Violation> {
    let modules: Vec<(String, &str)> = files
        .iter()
        .filter_map(|f| {
            let stem = f.path.rsplit('/').next()?.strip_suffix(".rs")?;
            stem.ends_with("_tests")
                .then(|| (stem.to_string(), f.text.as_str()))
        })
        .collect();
    let mut out = Vec::new();
    for f in files {
        for c in comments::of(&f.path, &f.text) {
            for (module, name) in cited(c.text) {
                let Some((_, text)) = modules.iter().find(|(m, _)| *m == module) else {
                    continue;
                };
                if !defines(text, name) {
                    out.push(Violation::new(
                        &f.path,
                        c.line,
                        format!(
                            "[symbols] `{module}::{name}` is not defined in `{module}.rs`; \
                             a citation a rename breaks has to break loudly"
                        ),
                    ));
                }
            }
        }
    }
    out
}

/// Every `<ident>_tests::<ident>[::<ident>]*` in `s`, as one `(module, segment)` per segment.
///
/// Each segment is resolved, not just the first: a citation that goes through an inner module
/// would otherwise be checked only as far as that module, and the test name -- the part a rename
/// moves -- would go unread.
fn cited(s: &str) -> Vec<(String, &str)> {
    let ident = |c: char| c.is_ascii_alphanumeric() || c == '_';
    let mut out = Vec::new();
    for (at, _) in s.match_indices("_tests::") {
        let start = s[..at].rfind(|c: char| !ident(c)).map_or(0, |p| p + 1);
        // A qualified path names the same module, so a leading `crate::` or a longer prefix
        // changes nothing: only the identifier immediately before `_tests` is the module.
        let module = &s[start..at + "_tests".len()];
        if module.starts_with(|c: char| c.is_ascii_digit()) {
            continue;
        }
        let mut rest = &s[at + "_tests::".len()..];
        loop {
            let end = rest.find(|c: char| !ident(c)).unwrap_or(rest.len());
            if end == 0 {
                break;
            }
            out.push((module.to_string(), &rest[..end]));
            match rest[end..].strip_prefix("::") {
                Some(r) => rest = r,
                None => break,
            }
        }
    }
    out
}

/// Whether `text` defines `name` as an item: a cited test is a `fn`, a cited group of them a
/// `mod`, and the rest are what else a citation may reasonably point at.
fn defines(text: &str, name: &str) -> bool {
    const KINDS: [&str; 7] = [
        "fn ", "mod ", "struct ", "enum ", "const ", "type ", "trait ",
    ];
    KINDS.iter().any(|k| {
        text.match_indices(k).any(|(at, _)| {
            let tail = &text[at + k.len()..];
            tail.strip_prefix(name)
                .is_some_and(|r| !r.starts_with(|c: char| c.is_ascii_alphanumeric() || c == '_'))
        })
    })
}

/// The first `<ident>.rs:<digit>` / `<ident>.rs#L<digit>` in `s`, with its name and separator.
fn find(s: &str) -> Option<&str> {
    let b = s.as_bytes();
    let ident = |c: u8| c.is_ascii_alphanumeric() || c == b'_';
    for (at, _) in s.match_indices(".rs") {
        // A closing backtick may sit between the name and the separator.
        let tail = b[at + 3..].strip_prefix(b"`").unwrap_or(&b[at + 3..]);
        let ticks = b.len() - at - 3 - tail.len();
        let sep = match tail {
            [b':', d, ..] if d.is_ascii_digit() => 1,
            [b'#', b'L', d, ..] if d.is_ascii_digit() => 2,
            _ => continue,
        };
        if at == 0 || !ident(b[at - 1]) {
            continue;
        }
        let start = b[..at]
            .iter()
            .rposition(|&c| !(ident(c) || c == b'/' || c == b'-'));
        let end = at
            + 3
            + ticks
            + sep
            + tail[sep..]
                .iter()
                .take_while(|c| c.is_ascii_digit())
                .count();
        return Some(&s[start.map_or(0, |p| p + 1)..end]);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::super::tests::line_cite;
    use super::*;

    fn hits(path: &str, text: &str) -> usize {
        check(&[File::new(path, text)]).len()
    }

    #[test]
    fn a_cited_test_has_to_exist_and_an_unknown_module_is_skipped() {
        let m = File::new(
            "crates/h/src/quat_tests.rs",
            "mod inner { fn kept() {} }\nfn here() {}",
        );
        let cite = |s: &str| {
            let files = vec![
                File::new("crates/h/src/quat_tests.rs", &m.text),
                File::new("crates/h/src/quat.rs", &format!("/// {s}\nfn q() {{}}")),
            ];
            check_symbols(&files).len()
        };
        assert_eq!(cite("see `quat_tests::here`"), 0);
        // A module, and a module-qualified path: only the segment after `_tests::` is resolved.
        assert_eq!(cite("see `quat_tests::inner`"), 0);
        assert_eq!(cite("see `quat_tests::inner::kept`"), 0);
        // Every segment is resolved, so a rename of the last one is caught.
        assert_eq!(cite("see `quat_tests::inner::gone`"), 1);
        assert_eq!(cite("see `crate::quat_tests::here`"), 0);
        assert_eq!(cite("see `quat_tests::gone`"), 1);
        // A prefix of a real name is not that name.
        assert_eq!(cite("see `quat_tests::her`"), 1);
        // A module this tree does not hold is skipped, not failed.
        assert_eq!(cite("see `foreign_tests::whatever`"), 0);
        // Two in one comment are two violations.
        assert_eq!(cite("`quat_tests::a` and `quat_tests::b`"), 2);
    }

    #[test]
    fn flags_colon_and_anchor_forms() {
        assert_eq!(find(&line_cite("src/lib", 42)), Some("src/lib.rs:42"));
        assert_eq!(find(&format!("a/b.rs#L{}-L{}", 7, 9)), Some("a/b.rs#L7"));
        assert_eq!(find(&format!("x.rs:{}:{}", 1, 2)), Some("x.rs:1"));
    }

    #[test]
    fn digits_in_names_later_matches_and_backticks() {
        assert_eq!(find("see v1.rs:12"), Some("v1.rs:12"));
        assert_eq!(find("see mod.rs and lib.rs:12"), Some("lib.rs:12"));
        assert_eq!(find("see `lib.rs`:12"), Some("lib.rs`:12"));
        assert_eq!(find("see `lib.rs` for 12"), None);
    }

    #[test]
    fn placeholders_and_prose_pass() {
        assert_eq!(find("a `path.rs:NNN` citation"), None);
        assert_eq!(find("see .rs:12 or the mod.rs file"), None);
        assert_eq!(find("`crate::foo` in lib.rs, then rs:5"), None);
    }

    #[test]
    fn markdown_flags_every_line() {
        assert_eq!(
            hits("docs/A.md", &format!("ok\n{}", line_cite("foo", 3))),
            1
        );
    }

    #[test]
    fn code_flags_comments_only() {
        let src = format!(
            "// {}\nlet s = \"{}\";\n",
            line_cite("a", 1),
            line_cite("b", 2)
        );
        assert_eq!(hits("crates/x/src/lib.rs", &src), 1);
        let cfg = format!("# {}key = \"x\"\n", line_cite("c", 3));
        assert_eq!(hits("Cargo.toml", &cfg), 1);
        assert_eq!(
            hits("Cargo.toml", &format!("key = \"{}\"\n", line_cite("c", 3))),
            0
        );
    }

    #[test]
    fn reports_the_line() {
        let v = check(&[File::new("a.md", &format!("x\ny\n{}", line_cite("f", 5)))]);
        assert_eq!(v.iter().map(|v| v.line).collect::<Vec<_>>(), [3]);
    }
}
