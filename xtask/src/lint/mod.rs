//! `cargo xtask lint` (`docs/PHASE1.md` §3). Each check is a pure function from a virtual file
//! set (or parsed `cargo metadata`) to violations; only [`load_tree`] and `metadata::load` touch
//! the disk, so every check tests on fixtures.

mod citations;
mod closure;
mod comments;
mod drafts;
mod generated;
mod metadata;
mod sweep;

use std::fmt;
use std::io::ErrorKind;
use std::path::Path;
use std::process::Command;

use metadata::Metadata;

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

/// Checks over the parsed `cargo metadata` of the workspace.
type ManifestCheck = fn(&Metadata) -> Vec<Violation>;

const MANIFEST_CHECKS: &[ManifestCheck] = &[closure::check_metadata, sweep::check_metadata];

/// Runs every check over `files`; sorted by path then line.
pub(crate) fn check_all(files: &[File]) -> Vec<Violation> {
    let mut all: Vec<Violation> = CHECKS.iter().flat_map(|c| c(files)).collect();
    all.sort();
    all
}

/// Runs every manifest check over `meta`; sorted by path then line.
pub(crate) fn check_manifests(meta: &Metadata) -> Vec<Violation> {
    let mut all: Vec<Violation> = MANIFEST_CHECKS.iter().flat_map(|c| c(meta)).collect();
    all.sort();
    all
}

/// Lints the checkout that contains this crate.
pub(crate) fn run() -> Result<Vec<Violation>, String> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .ok_or("xtask has no parent dir")?;
    check_root(root)
}

/// Every check over the workspace rooted at `root`, files and manifests together.
fn check_root(root: &Path) -> Result<Vec<Violation>, String> {
    let mut all = check_all(&load_tree(root)?);
    all.extend(check_manifests(&metadata::load(&root.join("Cargo.toml"))?));
    all.sort();
    Ok(all)
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
        assert!(
            load_tree(root)?.len() > 20,
            "the tree loaded suspiciously few files"
        );
        let violations = check_root(root)?;
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
        let mut files = vec![
            File::new("docs/decisions/0013-x.md", "**Status:** draft\n"),
            File::new("docs/A.md", &line_cite("foo", 3)),
            File::new("crates/a/src/lib.rs", "// settled by 0013\n"),
            File::new("docs/G.md", "<!-- @generated -->\n"),
        ];
        files.extend(generated::stubs());
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

    /// One violation of each manifest check, so one dropped from `MANIFEST_CHECKS` fails here.
    #[test]
    fn every_manifest_check_runs() {
        use metadata::fixture::Fx;
        let m = Metadata::fixture(&[
            Fx::member("helicoid").deps(&[("nalgebra", None, false, &["__sweep"])]),
            Fx::dep("nalgebra"),
        ]);
        let out: Vec<String> = check_manifests(&m)
            .iter()
            .map(ToString::to_string)
            .collect();
        for tag in ["[closure]", "[sweep]"] {
            assert!(out.iter().any(|l| l.contains(tag)), "{tag}: {out:?}");
        }
    }

    /// `check_root`, the gate itself, runs the manifest checks: a planted `nalgebra` in a tracked
    /// temporary workspace comes back from it.
    #[test]
    fn check_root_runs_the_manifest_checks() -> Result<(), String> {
        let ok = "[dependencies]\nlibm = { path = \"../libm\" }\n";
        let helicoid = format!(
            "{ok}helicoid-linalg = {{ path = \"../helicoid-linalg\" }}\nnalgebra = {{ path = \"../nalgebra\" }}\n"
        );
        let crates = [
            ("libm", ""),
            ("nalgebra", ""),
            ("helicoid-linalg", ok),
            ("helicoid", helicoid.as_str()),
        ];
        let manifest = metadata::tests::workspace("root", &crates)?;
        let root = manifest.parent().ok_or("no parent")?;
        let git = |args: &[&str]| {
            let out = Command::new("git").arg("-C").arg(root).args(args).output();
            out.map_err(|e| e.to_string()).and_then(|o| {
                o.status
                    .success()
                    .then_some(())
                    .ok_or(format!("git {args:?}"))
            })
        };
        for stub in generated::stubs() {
            let path = root.join(&stub.path);
            let dir = path.parent().ok_or("no parent")?;
            std::fs::create_dir_all(dir).map_err(|e| e.to_string())?;
            std::fs::write(&path, &stub.text).map_err(|e| e.to_string())?;
        }
        let out = git(&["init", "-q"])
            .and_then(|()| git(&["add", "."]))
            .and_then(|()| check_root(root));
        let _ = std::fs::remove_dir_all(root);
        let out: Vec<String> = out?.iter().map(ToString::to_string).collect();
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(out[0].contains("[closure] `nalgebra`"), "{out:?}");
        Ok(())
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
        let mut files = vec![
            File::new("b.md", &line_cite("x", 3)),
            File::new("a.md", &line_cite("y", 9)),
        ];
        files.extend(generated::stubs());
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
