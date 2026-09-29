//! `cargo xtask lint` (`docs/PHASE1.md` §3). Each check is a pure function from a virtual file
//! set to violations; only [`load_tree`] touches the disk, so every check tests on fixtures.

mod citations;
mod comments;
mod drafts;
mod generated;

use std::fmt;
use std::io::ErrorKind;
use std::path::Path;
use std::process::Command;

/// One tracked text file: repo-relative `/` path and content.
pub(crate) struct File {
    pub(crate) path: String,
    pub(crate) text: String,
}

impl File {
    #[cfg(test)]
    fn new(path: &str, text: &str) -> Self {
        Self {
            path: path.into(),
            text: text.into(),
        }
    }
}

/// `path:line: message`, 1-based line.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Violation {
    pub(crate) path: String,
    pub(crate) line: usize,
    pub(crate) message: String,
}

impl Violation {
    pub(crate) fn new(path: &str, line: usize, message: impl Into<String>) -> Self {
        Self {
            path: path.into(),
            line,
            message: message.into(),
        }
    }
}

impl fmt::Display for Violation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}: {}", self.path, self.line, self.message)
    }
}

type Check = fn(&[File]) -> Vec<Violation>;

/// Independent checks; add a line here and nothing else.
const CHECKS: &[Check] = &[citations::check, drafts::check, generated::check_registry];

/// Runs every check over `files`; sorted by path then line.
pub(crate) fn check_all(files: &[File]) -> Vec<Violation> {
    let mut all: Vec<Violation> = CHECKS.iter().flat_map(|c| c(files)).collect();
    all.sort();
    all
}

/// Lints the checkout that contains this crate.
pub(crate) fn run() -> Result<Vec<Violation>, String> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .ok_or("xtask has no parent dir")?;
    Ok(check_all(&load_tree(root)?))
}

/// Tracked regular files as text (`git add` is the precondition, so untracked scratch files never
/// fail the lint). Symlinks, submodules, binary and deleted files are skipped.
pub(crate) fn load_tree(root: &Path) -> Result<Vec<File>, String> {
    let out = Command::new("git")
        .arg("-C")
        .arg(root)
        .args(["ls-files", "-z", "--cached"])
        .output()
        .map_err(|e| format!("cannot run git: {e}"))?;
    if !out.status.success() {
        return Err(format!(
            "git ls-files failed: {}",
            String::from_utf8_lossy(&out.stderr)
        ));
    }
    let names = String::from_utf8(out.stdout).map_err(|e| format!("non-UTF-8 path: {e}"))?;
    let mut files = Vec::new();
    for path in names.split('\0').filter(|p| !p.is_empty()) {
        let full = root.join(path);
        match std::fs::symlink_metadata(&full) {
            Ok(m) if m.is_file() => {}
            Ok(_) => continue,
            Err(e) if e.kind() == ErrorKind::NotFound => continue,
            Err(e) => return Err(format!("{path}: {e}")),
        }
        match std::fs::read_to_string(&full) {
            Ok(text) => files.push(File {
                path: path.into(),
                text,
            }),
            Err(e) if matches!(e.kind(), ErrorKind::NotFound | ErrorKind::InvalidData) => {}
            Err(e) => return Err(format!("{path}: {e}")),
        }
    }
    Ok(files)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn real_tree_is_clean() -> Result<(), String> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .ok_or("no parent")?;
        let files = load_tree(root)?;
        assert!(files.len() > 20, "the tree loaded suspiciously few files");
        let violations = check_all(&files);
        let report: Vec<String> = violations.iter().map(ToString::to_string).collect();
        assert!(
            violations.is_empty(),
            "lint violations:\n{}",
            report.join("\n")
        );
        Ok(())
    }

    /// One violation of each kind, so a check dropped from `CHECKS` fails here.
    #[test]
    fn every_check_runs() {
        let files = [
            File::new("docs/decisions/0013-x.md", "**Status:** draft\n"),
            File::new("docs/A.md", &line_cite("foo", 3)),
            File::new("crates/a/src/lib.rs", "// settled by 0013\n"),
            File::new("docs/G.md", "<!-- @generated -->\n"),
        ];
        let out: Vec<String> = check_all(&files).iter().map(ToString::to_string).collect();
        for tag in ["[citations]", "[drafts]", "[generated]"] {
            assert_eq!(
                out.iter().filter(|l| l.contains(tag)).count(),
                1,
                "{tag}: {out:?}"
            );
        }
        assert_eq!(out.len(), 3, "{out:?}");
    }

    #[cfg(unix)]
    #[test]
    fn load_tree_reads_tracked_regular_files_only() -> Result<(), String> {
        let io = |e: std::io::Error| e.to_string();
        let root = std::env::temp_dir().join(format!("xtask-lint-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(root.join("d")).map_err(io)?;
        for (name, bytes) in [
            ("d/a.md", &b"text\n"[..]),
            ("bin.dat", &[0xff, 0xfe][..]),
            ("gone.md", b"x"),
            ("untracked.md", b"scratch"),
        ] {
            std::fs::write(root.join(name), bytes).map_err(io)?;
        }
        std::os::unix::fs::symlink("d", root.join("dirlink")).map_err(io)?;
        let git = |args: &[&str]| {
            let ok = Command::new("git").arg("-C").arg(&root).args(args).output();
            ok.map_err(io).and_then(|o| {
                o.status
                    .success()
                    .then_some(())
                    .ok_or(format!("git {args:?}"))
            })
        };
        git(&["init", "-q"])?;
        git(&["add", "d/a.md", "bin.dat", "gone.md", "dirlink"])?;
        std::fs::remove_file(root.join("gone.md")).map_err(io)?;
        let files = load_tree(&root);
        let _ = std::fs::remove_dir_all(&root);
        let names: Vec<String> = files?.into_iter().map(|f| f.path).collect();
        assert_eq!(names, ["d/a.md"]);
        Ok(())
    }

    #[test]
    fn violations_format_and_sort() {
        let files = [
            File::new("b.md", &line_cite("x", 3)),
            File::new("a.md", &line_cite("y", 9)),
        ];
        let out: Vec<String> = check_all(&files).iter().map(ToString::to_string).collect();
        assert_eq!(out.len(), 2);
        assert!(out[0].starts_with("a.md:1: "), "{out:?}");
        assert!(out[1].starts_with("b.md:1: "), "{out:?}");
    }

    /// A line citation built at runtime, so this file never contains one.
    pub(super) fn line_cite(stem: &str, n: u32) -> String {
        format!("see {stem}.rs:{n}\n")
    }
}
