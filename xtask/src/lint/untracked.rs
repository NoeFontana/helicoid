//! Check 4: a file git does not have is a file no other check reads.
//!
//! Every other check runs over `load_tree`'s output, which is `git ls-files --cached`: tracked
//! files, by design, so that scratch work in the checkout cannot fail the lint. The cost of that
//! design is silence in the one case that matters — a file written and not yet added is invisible,
//! so the checks that would have read it pass, and the commit does not contain it either.
//!
//! Measured on this repository: a citation of `0033` (draft) in `docs/PHASE1.md` §0.0 was accepted
//! by `cargo xtask lint` and by `just lint`, and became a `drafts` violation the moment `git add`
//! made the record visible — after the commit that claimed the lint was clean. Every new decision
//! record is untracked at the moment it is written, so this is the normal path and not an unlucky
//! one.
//!
//! `--exclude-standard` applies `.gitignore`, so a build artefact or an ignored baseline is not
//! reported. What is reported is a file inside the checkout that git has never been told about.

use std::path::{Path, PathBuf};
use std::process::Command;

use super::Violation;

/// Untracked, non-ignored paths, as `git` reports them from `root`.
pub(crate) fn load(root: &Path) -> Result<Vec<PathBuf>, String> {
    let out = Command::new("git")
        .arg("-C")
        .arg(root)
        .args(["ls-files", "-z", "--others", "--exclude-standard"])
        .output()
        .map_err(|e| format!("cannot run git: {e}"))?;
    if !out.status.success() {
        return Err(format!(
            "git ls-files --others failed: {}",
            String::from_utf8_lossy(&out.stderr)
        ));
    }
    let names = String::from_utf8(out.stdout).map_err(|e| format!("non-UTF-8 path: {e}"))?;
    Ok(names
        .split('\0')
        .filter(|p| !p.is_empty())
        .map(PathBuf::from)
        .collect())
}

/// One violation per untracked path, at line 1: there is no line to blame, only a missing `git add`.
pub(crate) fn check(paths: &[PathBuf]) -> Vec<Violation> {
    paths
        .iter()
        .map(|p| {
            Violation::new(
                &p.display().to_string(),
                1,
                "[untracked] every other check reads tracked files only, so this file is \
                 unchecked and is not in the next commit either; `git add` it, or ignore it",
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Each untracked path is reported once, at line 1, naming the remedy.
    #[test]
    fn every_untracked_path_is_reported() {
        let paths = [
            PathBuf::from("docs/decisions/0034-a-new-record.md"),
            PathBuf::from("crates/helicoid/src/new.rs"),
        ];
        let out: Vec<String> = check(&paths).iter().map(ToString::to_string).collect();
        assert_eq!(out.len(), 2, "{out:?}");
        assert!(
            out[0].starts_with("docs/decisions/0034-a-new-record.md:1: [untracked]"),
            "{out:?}"
        );
        assert!(out[0].contains("git add"), "{out:?}");
        assert!(
            out[1].starts_with("crates/helicoid/src/new.rs:1: [untracked]"),
            "{out:?}"
        );
    }

    /// A clean checkout reports nothing, which is the state the lint is normally run in.
    #[test]
    fn a_clean_checkout_is_silent() {
        assert!(check(&[]).is_empty());
    }

    /// `load` reads the real `git`: a planted file is reported, and an ignored one is not.
    #[test]
    fn git_reports_the_planted_file_and_not_the_ignored_one() -> Result<(), String> {
        let root = std::env::temp_dir().join("helicoid-lint-untracked");
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(root.join("docs")).map_err(|e| e.to_string())?;
        let git = |args: &[&str]| -> Result<(), String> {
            let ok = Command::new("git")
                .arg("-C")
                .arg(&root)
                .args(args)
                .output()
                .map_err(|e| format!("git {args:?}: {e}"))?;
            ok.status
                .success()
                .then_some(())
                .ok_or(format!("git {args:?}"))
        };
        git(&["init", "-q"])?;
        std::fs::write(root.join(".gitignore"), "ignored.md\n").map_err(|e| e.to_string())?;
        std::fs::write(root.join("docs/new.md"), "x\n").map_err(|e| e.to_string())?;
        std::fs::write(root.join("ignored.md"), "x\n").map_err(|e| e.to_string())?;
        git(&["add", ".gitignore"])?;
        let paths = load(&root)?;
        let shown: Vec<String> = paths.iter().map(|p| p.display().to_string()).collect();
        assert!(shown.contains(&"docs/new.md".to_string()), "{shown:?}");
        assert!(!shown.iter().any(|p| p.contains("ignored")), "{shown:?}");
        Ok(())
    }
}
