//! Check 1: no `<name>.rs:<digits>` or `<name>.rs#L<digits>` line citation in markdown or in the
//! comments of code and config files (`docs/PHASE1.md` §3). Cite a symbol instead.

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
