//! `cargo xtask release-notes VERSION [--out FILE]` (`0100`): the notes of the GitHub release
//! `release.yml` cuts for the tag `vVERSION`, from what the repository already keeps.
//!
//! - **The changelog's section** `## [VERSION]`, verbatim, its relative links made absolute at the
//!   tag. It is the curated half, and the only one a human writes.
//! - **The decision records** added or amended since the previous `v*` tag (`git diff`), each with
//!   its title and its `**Status:**` line at the tag. A `draft` record's change ships, but the
//!   record is not yet the contract, and the table says so.
//! - **The crates and the MSRV**, from the workspace manifest.
//!
//! It fails, and the release with it, where the tag is not the workspace's version or the changelog
//! has no non-empty section for it: a release is cut from a stamped changelog (`CONTRIBUTING.md`).
//! GitHub's own list of the merged pull requests is appended by `gh release create --generate-notes`.

use std::path::Path;
use std::process::Command;

const USAGE: &str = "usage: cargo xtask release-notes VERSION [--out FILE]";
const CRATES: [&str; 2] = ["helicoid", "helicoid-linalg"];

pub(crate) fn run(args: &[String]) -> Result<(), String> {
    let (version, out) = match args {
        [v] => (v.as_str(), None),
        [v, flag, file] if flag == "--out" => (v.as_str(), Some(file.as_str())),
        _ => return Err(USAGE.into()),
    };
    let root = crate::conformance::root()?;
    let manifest = read(&root.join("Cargo.toml"))?;
    let ws = Workspace::parse(&manifest)?;
    if ws.version != version {
        return Err(format!(
            "the tag is v{version} and the workspace is {}: bump Cargo.toml, or tag v{}",
            ws.version, ws.version
        ));
    }
    let section = changelog_section(&read(&root.join("CHANGELOG.md"))?, version)?;
    let tag = format!("v{version}");
    let since = previous_tag(&git(&root, &["tag", "--list", "v*"])?, version);
    let changed = match &since {
        Some(prev) => git(
            &root,
            &[
                "diff",
                "--name-status",
                "--find-renames",
                &format!("{prev}..HEAD"),
                "--",
                "docs/decisions",
            ],
        )?,
        None => git(&root, &["ls-files", "docs/decisions"])?
            .lines()
            .map(|l| format!("A\t{l}"))
            .collect::<Vec<_>>()
            .join("\n"),
    };
    let mut records = Vec::new();
    for (path, how) in decision_paths(&changed) {
        let text = read(&root.join(&path))?;
        records.push(Record::parse(&path, &text, how)?);
    }
    let notes = render(&ws, &tag, since.as_deref(), &section, &records);
    match out {
        Some(file) => std::fs::write(file, notes).map_err(|e| format!("{file}: {e}"))?,
        None => std::io::Write::write_all(&mut std::io::stdout(), notes.as_bytes())
            .map_err(|e| format!("stdout: {e}"))?,
    }
    Ok(())
}

fn read(path: &Path) -> Result<String, String> {
    std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))
}

fn git(root: &Path, args: &[&str]) -> Result<String, String> {
    let out = Command::new("git")
        .current_dir(root)
        .args(args)
        .output()
        .map_err(|e| format!("git {}: {e}", args.join(" ")))?;
    if !out.status.success() {
        return Err(format!(
            "git {}: {}",
            args.join(" "),
            String::from_utf8_lossy(&out.stderr).trim()
        ));
    }
    String::from_utf8(out.stdout).map_err(|e| format!("git {}: {e}", args.join(" ")))
}

/// The workspace manifest's `[workspace.package]` fields the notes print.
#[derive(Debug, PartialEq)]
struct Workspace {
    version: String,
    rust_version: String,
    repository: String,
}

impl Workspace {
    fn parse(manifest: &str) -> Result<Self, String> {
        let mut in_package = false;
        let (mut version, mut rust_version, mut repository) = (None, None, None);
        for line in manifest.lines() {
            let line = line.trim();
            if line.starts_with('[') {
                in_package = line == "[workspace.package]";
                continue;
            }
            if !in_package {
                continue;
            }
            let Some((key, value)) = line.split_once('=') else {
                continue;
            };
            let value = value.trim().trim_matches('"').to_string();
            match key.trim() {
                "version" => version = Some(value),
                "rust-version" => rust_version = Some(value),
                "repository" => repository = Some(value),
                _ => {}
            }
        }
        let need =
            |v: Option<String>, k: &str| v.ok_or(format!("[workspace.package] has no `{k}`"));
        Ok(Self {
            version: need(version, "version")?,
            rust_version: need(rust_version, "rust-version")?,
            repository: need(repository, "repository")?,
        })
    }
}

/// The body of `## [VERSION] - date`: its date line kept, up to the next `## ` heading.
fn changelog_section(changelog: &str, version: &str) -> Result<String, String> {
    let head = format!("## [{version}]");
    let mut lines = changelog.lines().skip_while(|l| !l.starts_with(&head));
    let Some(first) = lines.next() else {
        return Err(format!(
            "CHANGELOG.md has no `{head}` section: stamp it before tagging"
        ));
    };
    let date = first[head.len()..].trim().trim_start_matches('-').trim();
    let body: Vec<&str> = lines.take_while(|l| !l.starts_with("## ")).collect();
    let body = body.join("\n").trim().to_string();
    if body.is_empty() {
        return Err(format!("CHANGELOG.md's `{head}` section is empty"));
    }
    Ok(match date.is_empty() {
        true => body,
        false => format!("_Released {date}._\n\n{body}"),
    })
}

/// `text` with its repository-relative links (`](./x)`, `](docs/x)`, `](CHANGELOG.md)`) made
/// absolute at `tag`, so they resolve on the release page. Absolute links and anchors are left.
fn absolute_links(text: &str, repository: &str, tag: &str) -> String {
    let base = format!("{repository}/blob/{tag}/");
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(i) = rest.find("](") {
        out.push_str(&rest[..i + 2]);
        rest = &rest[i + 2..];
        let target = rest.split(')').next().unwrap_or("");
        let relative = !(target.contains("://") || target.starts_with('#') || target.is_empty());
        if relative {
            out.push_str(&base);
            rest = rest.strip_prefix("./").unwrap_or(rest);
        }
    }
    out.push_str(rest);
    out
}

/// The highest `vX.Y.Z` tag below `version`, by numeric comparison.
fn previous_tag(tags: &str, version: &str) -> Option<String> {
    let parse = |v: &str| -> Option<Vec<u64>> { v.split('.').map(|p| p.parse().ok()).collect() };
    let current = parse(version)?;
    tags.lines()
        .filter_map(|t| {
            t.trim()
                .strip_prefix('v')
                .and_then(|v| Some((parse(v)?, t.trim())))
        })
        .filter(|(v, _)| *v < current)
        .max_by(|a, b| a.0.cmp(&b.0))
        .map(|(_, t)| t.to_string())
}

/// `git diff --name-status` lines naming a decision record (`NNNN-*.md`): its path at the tag, and
/// whether it was added or amended. A deleted record has no text to cite and is left out.
fn decision_paths(name_status: &str) -> Vec<(String, How)> {
    let mut out: Vec<(String, How)> = name_status
        .lines()
        .filter_map(|l| {
            let mut cols = l.split('\t');
            let status = cols.next()?;
            let path = cols.next_back()?.to_string();
            let file = path.rsplit('/').next()?;
            let numbered = file.len() > 5 && file[..4].bytes().all(|b| b.is_ascii_digit());
            if !numbered || !file.ends_with(".md") {
                return None;
            }
            match status.chars().next()? {
                'A' => Some((path, How::Added)),
                'M' | 'R' => Some((path, How::Amended)),
                _ => None,
            }
        })
        .collect();
    out.sort();
    out
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum How {
    Added,
    Amended,
}

/// A record as the notes cite it.
#[derive(Debug, PartialEq)]
struct Record {
    path: String,
    number: String,
    title: String,
    status: String,
    how: How,
}

impl Record {
    fn parse(path: &str, text: &str, how: How) -> Result<Self, String> {
        let heading = text
            .lines()
            .find(|l| l.starts_with("# "))
            .ok_or(format!("{path}: no `# NNNN: title` heading"))?;
        let (number, title) = heading[2..]
            .split_once(':')
            .ok_or(format!("{path}: the heading is not `# NNNN: title`"))?;
        let status = text
            .lines()
            .find_map(|l| l.strip_prefix("**Status:**"))
            .ok_or(format!("{path}: no `**Status:**` line"))?;
        Ok(Self {
            path: path.to_string(),
            number: number.trim().to_string(),
            title: title.trim().to_string(),
            status: status.trim().to_string(),
            how,
        })
    }
}

fn render(
    ws: &Workspace,
    tag: &str,
    since: Option<&str>,
    section: &str,
    records: &[Record],
) -> String {
    let mut out = absolute_links(section, &ws.repository, tag);
    out.push_str("\n\n## Decision records\n\n");
    let span = match since {
        Some(prev) => format!("since `{prev}`"),
        None => "in this first release".to_string(),
    };
    if records.is_empty() {
        out.push_str(&format!("No record was added or amended {span}.\n"));
    } else {
        let drafts = records
            .iter()
            .filter(|r| r.status.starts_with("draft"))
            .count();
        out.push_str(&format!(
            "{} added or amended {span}. A `draft` record's change ships in this release, but the \
             record is not yet the contract ({drafts} here).\n\n| Record | Status | |\n|---|---|---|\n",
            records.len()
        ));
        for r in records {
            let how = match r.how {
                How::Added => "added",
                How::Amended => "amended",
            };
            out.push_str(&format!(
                "| [`{}`]({}/blob/{tag}/{}) {} | {} | {how} |\n",
                r.number, ws.repository, r.path, r.title, r.status
            ));
        }
    }
    let version = tag.trim_start_matches('v');
    out.push_str(&format!(
        "\n## Crates\n\n{}, MSRV {}. `0.0.x` is not ordinary semver: every release may break every \
         other, so pin it exactly:\n\n```toml\nhelicoid = \"={version}\"\n```\n",
        CRATES
            .iter()
            .map(|c| format!("[`{c} {version}`](https://crates.io/crates/{c}/{version})"))
            .collect::<Vec<_>>()
            .join(" and "),
        ws.rust_version
    ));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const MANIFEST: &str = "[workspace]\nmembers = [\"a\"]\n\n[workspace.package]\nversion = \
        \"0.0.5\"\nedition = \"2021\"\nrust-version = \"1.85\"\nrepository = \
        \"https://github.com/o/r\"\n\n[workspace.dependencies]\nversion = \"9\"\n";

    #[test]
    fn the_workspace_fields_are_read_from_their_table_only() -> Result<(), String> {
        let ws = Workspace::parse(MANIFEST)?;
        assert_eq!(
            ws,
            Workspace {
                version: "0.0.5".into(),
                rust_version: "1.85".into(),
                repository: "https://github.com/o/r".into(),
            }
        );
        assert!(Workspace::parse("[workspace.package]\nversion = \"1\"\n").is_err());
        Ok(())
    }

    #[test]
    fn a_section_is_its_body_to_the_next_heading_and_must_exist() -> Result<(), String> {
        let log = "# Changelog\n\n## [Unreleased]\n\n## [0.0.5] - 2026-10-16\n\nWhat.\n\n### \
                   Changed\n\n- x\n\n## [0.0.4] - 2026-10-14\n\nOld.\n";
        assert_eq!(
            changelog_section(log, "0.0.5")?,
            "_Released 2026-10-16._\n\nWhat.\n\n### Changed\n\n- x"
        );
        assert!(changelog_section(log, "0.0.6").is_err());
        assert!(changelog_section("## [0.0.7]\n\n## [0.0.6]\n", "0.0.7").is_err());
        Ok(())
    }

    #[test]
    fn the_committed_changelog_has_a_section_for_the_workspace_version() -> Result<(), String> {
        let root = crate::conformance::root()?;
        let ws = Workspace::parse(&read(&root.join("Cargo.toml"))?)?;
        changelog_section(&read(&root.join("CHANGELOG.md"))?, &ws.version).map(|_| ())
    }

    #[test]
    fn relative_links_point_at_the_tag() {
        let text = "[a](./docs/x.md) [b](CHANGELOG.md) [c](https://e.org/y) [d](#h) [`0099`](./docs/decisions/0099-z.md)";
        assert_eq!(
            absolute_links(text, "https://github.com/o/r", "v0.0.5"),
            "[a](https://github.com/o/r/blob/v0.0.5/docs/x.md) \
             [b](https://github.com/o/r/blob/v0.0.5/CHANGELOG.md) [c](https://e.org/y) [d](#h) \
             [`0099`](https://github.com/o/r/blob/v0.0.5/docs/decisions/0099-z.md)"
        );
    }

    #[test]
    fn the_previous_tag_is_the_highest_below_by_number() {
        let tags = "v0.0.1\nv0.0.10\nv0.0.2\nv0.0.9\nnot-a-tag\nv0.1.0\n";
        assert_eq!(previous_tag(tags, "0.0.10").as_deref(), Some("v0.0.9"));
        assert_eq!(previous_tag(tags, "0.1.0").as_deref(), Some("v0.0.10"));
        assert_eq!(previous_tag(tags, "0.0.1"), None);
    }

    #[test]
    fn records_are_the_numbered_files_added_or_amended() -> Result<(), String> {
        let diff = "A\tdocs/decisions/0099-a.md\nM\tdocs/decisions/0089-b.md\nM\tdocs/decisions/\
                    README.md\nD\tdocs/decisions/0001-c.md\nR100\tdocs/decisions/0002-old.md\t\
                    docs/decisions/0002-new.md\nA\tdocs/decisions/template.md\n";
        assert_eq!(
            decision_paths(diff),
            [
                ("docs/decisions/0002-new.md".to_string(), How::Amended),
                ("docs/decisions/0089-b.md".to_string(), How::Amended),
                ("docs/decisions/0099-a.md".to_string(), How::Added),
            ]
        );
        let r = Record::parse(
            "docs/decisions/0099-a.md",
            "# 0099: The series sums smallest first\n\n**Status:** draft\n",
            How::Added,
        )?;
        assert_eq!(
            (r.number.as_str(), r.title.as_str(), r.status.as_str()),
            ("0099", "The series sums smallest first", "draft")
        );
        Ok(())
    }

    #[test]
    fn the_notes_carry_the_section_the_records_and_the_pin() -> Result<(), String> {
        let ws = Workspace::parse(MANIFEST)?;
        let r = Record::parse(
            "docs/decisions/0099-a.md",
            "# 0099: T\n**Status:** draft\n",
            How::Added,
        )?;
        let notes = render(
            &ws,
            "v0.0.5",
            Some("v0.0.4"),
            "Body [x](./docs/y.md).",
            &[r],
        );
        assert!(notes.starts_with("Body [x](https://github.com/o/r/blob/v0.0.5/docs/y.md)."));
        assert!(notes.contains("1 added or amended since `v0.0.4`"));
        assert!(notes.contains("| [`0099`](https://github.com/o/r/blob/v0.0.5/docs/decisions/0099-a.md) T | draft | added |"));
        assert!(notes.contains("helicoid = \"=0.0.5\""));
        assert!(notes.contains("MSRV 1.85"));
        Ok(())
    }
}
