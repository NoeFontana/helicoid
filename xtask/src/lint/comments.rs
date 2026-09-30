//! Comment extraction: what a "comment" is for each file kind, one entry per source line.

/// A comment fragment on one line (1-based).
pub(crate) struct Comment<'a> {
    pub(crate) line: usize,
    pub(crate) text: &'a str,
}

enum Kind {
    /// Prose: every line counts.
    Markdown,
    /// `//` and nested `/* */`, skipping string and char literals.
    Rust,
    /// `//` and `/* */` without literal tracking (C-family runners).
    Slash,
    /// `#` at a line start or after whitespace, outside double quotes (TOML, YAML, shell, Python,
    /// `justfile`, Make, ignore and attribute files). Whitespace-only opening keeps `${#x}`, `$#`
    /// and `it's` from hiding or faking a comment.
    Hash,
}

fn kind(path: &str) -> Option<Kind> {
    let name = path.rsplit('/').next().unwrap_or(path);
    let ext = name.rsplit_once('.').map(|(_, e)| e);
    match (name, ext) {
        (_, Some("md" | "markdown")) => Some(Kind::Markdown),
        (_, Some("rs")) => Some(Kind::Rust),
        (_, Some("c" | "cc" | "cpp" | "cxx" | "h" | "hh" | "hpp" | "js" | "ts" | "cu")) => {
            Some(Kind::Slash)
        }
        (
            "justfile" | "Justfile" | "Makefile" | "CMakeLists.txt" | "Dockerfile" | "CODEOWNERS"
            | ".gitattributes" | ".dockerignore",
            _,
        ) => Some(Kind::Hash),
        _ if name.starts_with("Dockerfile.") => Some(Kind::Hash),
        (
            _,
            Some(
                "toml" | "yml" | "yaml" | "py" | "sh" | "bash" | "ini" | "cmake" | "mk"
                | "dockerfile" | "gitignore" | "editorconfig",
            ),
        ) => Some(Kind::Hash),
        _ => None,
    }
}

/// Comments of `text`; empty for a file kind with no comment syntax we know.
pub(crate) fn of<'a>(path: &str, text: &'a str) -> Vec<Comment<'a>> {
    match kind(path) {
        Some(Kind::Markdown) => text
            .lines()
            .enumerate()
            .map(|(i, text)| Comment { line: i + 1, text })
            .collect(),
        Some(Kind::Rust) => slash(text, true),
        Some(Kind::Slash) => slash(text, false),
        Some(Kind::Hash) => hash(text),
        None => Vec::new(),
    }
}

fn hash(text: &str) -> Vec<Comment<'_>> {
    let mut out = Vec::new();
    for (i, l) in text.lines().enumerate() {
        let (mut quoted, mut escaped, mut prev) = (false, false, ' ');
        for (j, c) in l.char_indices() {
            match c {
                _ if escaped => escaped = false,
                '\\' if quoted => escaped = true,
                '"' => quoted = !quoted,
                '#' if !quoted && prev.is_whitespace() => {
                    out.push(Comment {
                        line: i + 1,
                        text: &l[j..],
                    });
                    break;
                }
                _ => {}
            }
            prev = c;
        }
    }
    out
}

/// Lexer for `//` and `/* */`; with `literals`, Rust strings (incl. raw) and chars are skipped so
/// a fixture string is never mistaken for a comment.
fn slash(text: &str, literals: bool) -> Vec<Comment<'_>> {
    let b = text.as_bytes();
    let (mut i, mut line, mut out) = (0, 1, Vec::new());
    while i < b.len() {
        match b[i] {
            b'\n' => {
                line += 1;
                i += 1;
            }
            b'/' if b.get(i + 1) == Some(&b'/') => {
                let end = text[i..].find('\n').map_or(b.len(), |n| i + n);
                out.push(Comment {
                    line,
                    text: &text[i..end],
                });
                i = end;
            }
            b'/' if b.get(i + 1) == Some(&b'*') => {
                let (mut depth, mut seg, mut j) = (1, i, i + 2);
                while j < b.len() && depth > 0 {
                    match (b[j], b.get(j + 1)) {
                        (b'/', Some(b'*')) if literals => (depth, j) = (depth + 1, j + 2),
                        (b'*', Some(b'/')) => (depth, j) = (depth - 1, j + 2),
                        (b'\n', _) => {
                            out.push(Comment {
                                line,
                                text: &text[seg..j],
                            });
                            (line, seg, j) = (line + 1, j + 1, j + 1);
                        }
                        _ => j += 1,
                    }
                }
                out.push(Comment {
                    line,
                    text: &text[seg..j.min(b.len())],
                });
                i = j;
            }
            b'"' if literals => {
                let (end, nl) = skip_string(b, i + 1);
                (i, line) = (end, line + nl);
            }
            b'r' if literals && raw_hashes(b, i).is_some() => {
                let hashes = raw_hashes(b, i).unwrap_or(0);
                let close = format!("\"{}", "#".repeat(hashes));
                let start = i + 2 + hashes;
                let end = text[start..]
                    .find(&close)
                    .map_or(b.len(), |n| start + n + close.len());
                line += b[i..end].iter().filter(|&&c| c == b'\n').count();
                i = end;
            }
            b'\'' if literals => i = skip_char_or_lifetime(text, i),
            _ => i += 1,
        }
    }
    out
}

/// `Some(n)` when `b[i]` starts a raw string `r#..#"` with `n` hashes (also `br`, `cr`).
fn raw_hashes(b: &[u8], i: usize) -> Option<usize> {
    let base = if i > 0 && matches!(b[i - 1], b'b' | b'c') {
        i - 1
    } else {
        i
    };
    if base > 0 && (b[base - 1].is_ascii_alphanumeric() || b[base - 1] == b'_') {
        return None;
    }
    let n = b[i + 1..].iter().take_while(|&&c| c == b'#').count();
    (b.get(i + 1 + n) == Some(&b'"')).then_some(n)
}

/// Index after the closing quote of a normal string starting at `i`, and newlines crossed.
fn skip_string(b: &[u8], mut i: usize) -> (usize, usize) {
    let mut nl = 0;
    while i < b.len() {
        match b[i] {
            b'\\' => {
                nl += usize::from(b.get(i + 1) == Some(&b'\n'));
                i += 2;
            }
            b'"' => return (i + 1, nl),
            c => {
                nl += usize::from(c == b'\n');
                i += 1;
            }
        }
    }
    (b.len(), nl)
}

fn skip_char_or_lifetime(text: &str, i: usize) -> usize {
    let rest = &text[i + 1..];
    let mut chars = rest.chars();
    match chars.next() {
        Some('\\') => {
            // Skip the escaped character first: it may itself be the quote (`'\''`).
            let body = 1 + chars.next().map_or(0, char::len_utf8);
            rest.get(body..)
                .and_then(|r| r.find('\''))
                .map_or(i + 1, |n| i + 1 + body + n + 1)
        }
        Some(c) if chars.next() == Some('\'') => i + 1 + c.len_utf8() + 1,
        _ => i + 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn texts(path: &str, src: &str) -> Vec<(usize, String)> {
        of(path, src)
            .iter()
            .map(|c| (c.line, c.text.to_string()))
            .collect()
    }

    #[test]
    fn rust_strings_are_not_comments() {
        let src = "let a = \"// no\"; // yes\nlet b = r#\"/* no \"# ;\nlet c = '\"'; let d: &'static str = \"x\"; // yes2\n";
        assert_eq!(
            texts("a.rs", src),
            [(1, "// yes".to_string()), (3, "// yes2".to_string())]
        );
    }

    #[test]
    fn rust_block_comments_nest_and_split_by_line() {
        let src = "/* a /* b */\n c */ x // d\n";
        assert_eq!(
            texts("a.rs", src),
            [
                (1, "/* a /* b */".to_string()),
                (2, " c */".to_string()),
                (2, "// d".to_string())
            ]
        );
    }

    #[test]
    fn hash_comments_respect_quotes() {
        let src = "k = \"a # b\" # c\n# d\n";
        assert_eq!(
            texts("Cargo.toml", src),
            [(1, "# c".to_string()), (2, "# d".to_string())]
        );
    }

    #[test]
    fn hash_comment_needs_leading_whitespace_and_ignores_apostrophes() {
        let src = "run: echo it's # a\necho ${#x} $# # b\nk = \"q\\\" # no\" # c\n";
        assert_eq!(
            texts("ci/a.yml", src),
            [
                (1, "# a".to_string()),
                (2, "# b".to_string()),
                (3, "# c".to_string())
            ]
        );
    }

    #[test]
    fn every_known_hash_and_slash_kind_has_comments() {
        for path in [
            "a.yml",
            "a.yaml",
            "a.sh",
            "a.bash",
            "a.py",
            "a.ini",
            "a.editorconfig",
            "a.gitignore",
            "justfile",
            "Makefile",
            "CMakeLists.txt",
            "CODEOWNERS",
            "Dockerfile",
            "Dockerfile.oracle",
            "x.dockerfile",
            ".gitattributes",
            ".dockerignore",
            "d/Cargo.toml",
        ] {
            assert_eq!(texts(path, "# c\n"), [(1, "# c".to_string())], "{path}");
        }
        for path in ["a.cpp", "a.cc", "a.cxx", "a.c", "a.h", "a.hpp", "a.js"] {
            assert_eq!(texts(path, "x // c\n"), [(1, "// c".to_string())], "{path}");
        }
        assert_eq!(texts("a.markdown", "p\n").len(), 1);
    }

    #[test]
    fn rust_literals_do_not_desynchronise_the_lexer() {
        let cases: [(&str, &[(usize, &str)]); 6] = [
            // Raw string whose interior holds a quote and a comment leader.
            ("let a = r#\"\"q\" // no\"#; // yes\n", &[(1, "// yes")]),
            // The hash count matters: `\"#` does not close an `r##` string.
            ("let a = r##\"x\"# // no\"##; // yes\n", &[(1, "// yes")]),
            // An escaped quote does not end a string; an escaped backslash does not escape it.
            (
                "let a = \"\\\" // no\"; // y1\nlet b = \"\\\\\"; // y2\n",
                &[(1, "// y1"), (2, "// y2")],
            ),
            // Newlines inside plain and raw strings advance the line counter.
            ("let a = \"x\ny\nz\"; // yes\n", &[(3, "// yes")]),
            ("let a = r\"x\ny\"; // yes\n", &[(2, "// yes")]),
            // Escaped-quote char literals.
            (
                "matches!(c, '\\''|'\"') // y1\nlet s = \"// no\"; '\\u{41}' // y2\n",
                &[(1, "// y1"), (2, "// y2")],
            ),
        ];
        for (src, want) in cases {
            let want: Vec<(usize, String)> = want.iter().map(|&(n, t)| (n, t.into())).collect();
            assert_eq!(texts("a.rs", src), want, "{src}");
        }
    }

    #[test]
    fn unknown_kind_has_no_comments() {
        assert!(texts("LICENSE-MIT", "# x\n// y\n").is_empty());
    }
}
