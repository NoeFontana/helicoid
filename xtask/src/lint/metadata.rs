//! The slice of `cargo metadata --format-version 1` that the manifest checks read. Only
//! [`load`] runs cargo; every check is a pure function of a parsed [`Metadata`], so it tests on
//! fixtures built by `Metadata::fixture` and, for the parsing itself, on a real temporary
//! workspace (`tests` below).
//!
//! The resolve graph is the whole workspace's, features unified across every member: a feature
//! that `xtask` enables on a dependency it shares with a library shows up in that library's
//! closure. The lints fail closed on that (a false positive names the crate, never a false
//! negative), and cargo offers no per-package resolve to avoid it.

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use std::process::Command;

use serde::Deserialize;

use super::Violation;

#[derive(Debug, Deserialize)]
pub(crate) struct Metadata {
    pub(crate) packages: Vec<Package>,
    pub(crate) workspace_members: Vec<String>,
    pub(crate) workspace_root: String,
    pub(crate) resolve: Resolve,
}

#[derive(Debug, Deserialize)]
pub(crate) struct Package {
    pub(crate) name: String,
    pub(crate) id: String,
    pub(crate) manifest_path: String,
    /// The `[features]` table: feature name to its members (`"feat"`, `"dep:x"`, `"x/feat"`).
    pub(crate) features: BTreeMap<String, Vec<String>>,
    pub(crate) dependencies: Vec<Dependency>,
}

/// A declared dependency; `kind` is `None` for a normal one, else `"dev"` or `"build"`.
#[derive(Debug, Deserialize)]
pub(crate) struct Dependency {
    pub(crate) name: String,
    /// The `package = ".."` rename, which is also the name of the implicit feature.
    pub(crate) rename: Option<String>,
    pub(crate) kind: Option<String>,
    pub(crate) optional: bool,
    /// Features this manifest requests on the dependency.
    pub(crate) features: Vec<String>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct Resolve {
    pub(crate) nodes: Vec<Node>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct Node {
    pub(crate) id: String,
    pub(crate) deps: Vec<NodeDep>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct NodeDep {
    pub(crate) pkg: String,
    pub(crate) dep_kinds: Vec<DepKind>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct DepKind {
    pub(crate) kind: Option<String>,
}

impl Dependency {
    /// The key the manifest declares it under: the rename if any, else the package name.
    pub(crate) fn key(&self) -> &str {
        self.rename.as_deref().unwrap_or(&self.name)
    }
}

impl Package {
    /// Whether `default` reaches, through `[features]`, a feature name or member (`"feat"`,
    /// `"dep:x"`, `"x/feat"`) for which `hit` holds.
    pub(crate) fn default_reaches(&self, hit: impl Fn(&str) -> bool) -> bool {
        let mut seen = BTreeSet::new();
        let mut stack = vec!["default"];
        while let Some(feature) = stack.pop() {
            if hit(feature) {
                return true;
            }
            if seen.insert(feature) {
                let members = self.features.get(feature).map_or(&[][..], Vec::as_slice);
                if members.iter().any(|m| hit(m)) {
                    return true;
                }
                stack.extend(
                    members
                        .iter()
                        .map(String::as_str)
                        .filter(|s| !s.contains(['/', ':'])),
                );
            }
        }
        false
    }
}

impl Metadata {
    pub(crate) fn package(&self, id: &str) -> Option<&Package> {
        self.packages.iter().find(|p| p.id == id)
    }

    /// The workspace members, `xtask` excluded: the crates whose manifests the checks bind.
    pub(crate) fn libraries(&self) -> impl Iterator<Item = &Package> {
        self.workspace_members
            .iter()
            .filter_map(|id| self.package(id))
            .filter(|p| p.name != "xtask")
    }

    /// Repo-relative `/` manifest path for `p`, for a [`Violation`].
    pub(crate) fn violation(&self, p: &Package, message: impl Into<String>) -> Violation {
        let rel = p
            .manifest_path
            .strip_prefix(&self.workspace_root)
            .map_or(p.manifest_path.as_str(), |r| r.trim_start_matches('/'));
        Violation::new(rel, 1, message)
    }
}

/// Runs cargo over the whole workspace: all features, every target platform, `--locked`. Needs
/// the registry index and crate manifests in cargo's cache (a cold cache needs the network).
pub(crate) fn load(manifest: &Path) -> Result<Metadata, String> {
    let cargo = std::env::var_os("CARGO").unwrap_or_else(|| "cargo".into());
    let out = Command::new(cargo)
        .args([
            "metadata",
            "--format-version",
            "1",
            "--locked",
            "--all-features",
        ])
        .arg("--manifest-path")
        .arg(manifest)
        .output()
        .map_err(|e| format!("cannot run cargo metadata: {e}"))?;
    if !out.status.success() {
        return Err(format!(
            "cargo metadata failed: {}",
            String::from_utf8_lossy(&out.stderr)
        ));
    }
    serde_json::from_slice(&out.stdout).map_err(|e| format!("cargo metadata output: {e}"))
}

#[cfg(test)]
pub(super) mod fixture {
    use super::*;

    /// A declared dependency: `(name, kind, optional, requested features)`.
    pub(crate) type D<'a> = (&'a str, Option<&'a str>, bool, &'a [&'a str]);

    /// A fixture package.
    #[derive(Clone, Copy)]
    pub(crate) struct Fx<'a> {
        pub(crate) name: &'a str,
        pub(crate) member: bool,
        pub(crate) deps: &'a [D<'a>],
        pub(crate) features: &'a [(&'a str, &'a [&'a str])],
    }

    impl<'a> Fx<'a> {
        pub(crate) fn member(name: &'a str) -> Self {
            Self {
                name,
                member: true,
                deps: &[],
                features: &[],
            }
        }
        pub(crate) fn dep(name: &'a str) -> Self {
            Self {
                member: false,
                ..Self::member(name)
            }
        }
        pub(crate) fn deps(self, deps: &'a [D<'a>]) -> Self {
            Self { deps, ..self }
        }
        pub(crate) fn features(self, features: &'a [(&'a str, &'a [&'a str])]) -> Self {
            Self { features, ..self }
        }
    }

    impl Metadata {
        /// Ids are package names; the resolve graph mirrors every declared dependency, so a
        /// fixture states each edge once.
        pub(crate) fn fixture(pkgs: &[Fx<'_>]) -> Self {
            let owned = |s: &str| s.to_owned();
            Self {
                workspace_root: "/ws".into(),
                workspace_members: pkgs
                    .iter()
                    .filter(|p| p.member)
                    .map(|p| owned(p.name))
                    .collect(),
                packages: pkgs
                    .iter()
                    .map(|p| Package {
                        name: owned(p.name),
                        id: owned(p.name),
                        manifest_path: format!("/ws/{}/Cargo.toml", p.name),
                        features: p
                            .features
                            .iter()
                            .map(|(f, m)| (owned(f), m.iter().map(|s| owned(s)).collect()))
                            .collect(),
                        dependencies: p
                            .deps
                            .iter()
                            .map(|&(name, kind, optional, feats)| Dependency {
                                name: owned(name),
                                rename: None,
                                kind: kind.map(owned),
                                optional,
                                features: feats.iter().map(|s| owned(s)).collect(),
                            })
                            .collect(),
                    })
                    .collect(),
                resolve: Resolve {
                    nodes: pkgs
                        .iter()
                        .map(|p| Node {
                            id: owned(p.name),
                            deps: p
                                .deps
                                .iter()
                                .map(|&(name, kind, ..)| NodeDep {
                                    pkg: owned(name),
                                    dep_kinds: vec![DepKind {
                                        kind: kind.map(owned),
                                    }],
                                })
                                .collect(),
                        })
                        .collect(),
                },
            }
        }
    }
}

/// The closure and `__sweep` checks over a real temporary workspace, through a real
/// `cargo metadata`: the parse and the `dep_kinds` semantics are cargo's, not the fixtures'.
#[cfg(test)]
pub(super) mod tests {
    use std::path::PathBuf;

    use super::super::check_manifests;
    use super::*;

    const MEMBERS: [&str; 3] = ["helicoid-linalg", "helicoid", "xtask"];

    /// Path-only stub crates (`(name, dependency/feature lines)`); returns the manifest path.
    /// Cargo.lock is generated offline, so no registry is touched.
    pub(crate) fn workspace(tag: &str, crates: &[(&str, &str)]) -> Result<PathBuf, String> {
        let io = |e: std::io::Error| e.to_string();
        let root = std::env::temp_dir().join(format!("xtask-ws-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        // The budget's own crates are members; everything else is a stub, as libm is in the real
        // tree (a path dependency inside the root would otherwise become a member).
        let quoted = |member: bool| -> String {
            let names: Vec<String> = crates
                .iter()
                .filter(|(n, _)| MEMBERS.contains(n) == member)
                .map(|(n, _)| format!("\"{n}\""))
                .collect();
            names.join(", ")
        };
        std::fs::create_dir_all(&root).map_err(io)?;
        std::fs::write(
            root.join("Cargo.toml"),
            format!(
                "[workspace]\nresolver = \"2\"\nmembers = [{}]\nexclude = [{}]\n",
                quoted(true),
                quoted(false)
            ),
        )
        .map_err(io)?;
        for (name, body) in crates {
            std::fs::create_dir_all(root.join(name).join("src")).map_err(io)?;
            std::fs::write(root.join(name).join("src/lib.rs"), "").map_err(io)?;
            let manifest = format!(
                "[package]\nname = \"{name}\"\nversion = \"0.0.0\"\nedition = \"2021\"\n{body}"
            );
            std::fs::write(root.join(name).join("Cargo.toml"), manifest).map_err(io)?;
        }
        let manifest = root.join("Cargo.toml");
        let cargo = std::env::var_os("CARGO").unwrap_or_else(|| "cargo".into());
        let lock = Command::new(cargo)
            .args(["generate-lockfile", "--offline", "--manifest-path"])
            .arg(&manifest)
            .output()
            .map_err(io)?;
        if lock.status.success() {
            Ok(manifest)
        } else {
            let _ = std::fs::remove_dir_all(&root);
            Err(String::from_utf8_lossy(&lock.stderr).into_owned())
        }
    }

    fn remove(manifest: &Path) {
        let _ = std::fs::remove_dir_all(manifest.parent().unwrap_or(manifest));
    }

    const OK: &str = "[dependencies]\nhelicoid-linalg = { path = \"../helicoid-linalg\" }\nlibm = { path = \"../libm\" }\n";
    const LEAF: &str = "";

    /// Manifest violations of the budget's crates (`libm`, `mint`, `helicoid-linalg`, `helicoid`)
    /// with the given `libm` and `helicoid` bodies, plus `extra` crates.
    fn lint(
        tag: &str,
        libm: &str,
        helicoid: &str,
        extra: &[(&str, &str)],
    ) -> Result<Vec<String>, String> {
        let linalg = "[features]\n__sweep = []\n[dependencies]\nlibm = { path = \"../libm\" }\nmint = { path = \"../mint\", optional = true }\n";
        let mut crates = vec![
            ("libm", libm),
            ("mint", LEAF),
            ("helicoid-linalg", linalg),
            ("helicoid", helicoid),
        ];
        crates.extend_from_slice(extra);
        let manifest = workspace(tag, &crates)?;
        let meta = load(&manifest);
        remove(&manifest);
        Ok(check_manifests(&meta?)
            .iter()
            .map(ToString::to_string)
            .collect())
    }

    #[test]
    fn the_budget_passes_with_dev_edges_optional_mint_and_an_unrestricted_xtask(
    ) -> Result<(), String> {
        let helicoid = format!("{OK}mint = {{ path = \"../mint\", optional = true }}\n[dev-dependencies]\nnalgebra = {{ path = \"../nalgebra\" }}\n");
        let xtask = "[dependencies]\nnalgebra = { path = \"../nalgebra\" }\nhelicoid = { path = \"../helicoid\" }\nhelicoid-linalg = { path = \"../helicoid-linalg\", features = [\"__sweep\"] }\n";
        let out = lint(
            "ok",
            LEAF,
            &helicoid,
            &[("nalgebra", LEAF), ("xtask", xtask)],
        )?;
        assert_eq!(out, Vec::<String>::new());
        Ok(())
    }

    #[test]
    fn a_planted_nalgebra_is_rejected_with_its_path() -> Result<(), String> {
        let helicoid = format!("{OK}nalgebra = {{ path = \"../nalgebra\" }}\n");
        let out = lint("nalgebra", LEAF, &helicoid, &[("nalgebra", LEAF)])?;
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].starts_with("helicoid/Cargo.toml:1: [closure] `nalgebra`"),
            "{out:?}"
        );
        assert!(out[0].ends_with("helicoid -> nalgebra"), "{out:?}");
        Ok(())
    }

    #[test]
    fn a_planted_serde_behind_libm_is_rejected_for_both_crates() -> Result<(), String> {
        let libm = "[dependencies]\nserde = { path = \"../serde\" }\n";
        let out = lint("serde", libm, OK, &[("serde", LEAF)])?;
        assert_eq!(out.len(), 2, "{out:?}");
        assert!(
            out[0].ends_with("helicoid-linalg -> libm -> serde"),
            "{out:?}"
        );
        assert!(out[1].ends_with("helicoid -> libm -> serde"), "{out:?}");
        Ok(())
    }

    #[test]
    fn a_library_enabling_sweep_on_a_dependency_is_rejected() -> Result<(), String> {
        let helicoid = OK.replace(
            "{ path = \"../helicoid-linalg\" }",
            "{ path = \"../helicoid-linalg\", features = [\"__sweep\"] }",
        );
        let out = lint("sweep", LEAF, &helicoid, &[])?;
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].contains("[sweep] `helicoid` enables `__sweep` on dependency `helicoid-linalg`"),
            "{out:?}"
        );
        Ok(())
    }
    #[test]
    fn an_optional_offender_behind_a_feature_is_seen_because_all_features_is_on(
    ) -> Result<(), String> {
        let helicoid = format!(
            "[features]\nextra = [\"dep:nalgebra\"]\n{OK}nalgebra = {{ path = \"../nalgebra\", optional = true }}\n"
        );
        let out = lint("optional", LEAF, &helicoid, &[("nalgebra", LEAF)])?;
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(out[0].ends_with("helicoid -> nalgebra"), "{out:?}");
        Ok(())
    }

    #[test]
    fn a_normal_dependency_that_is_also_a_dev_dependency_is_still_normal() -> Result<(), String> {
        let helicoid = format!(
            "{OK}nalgebra = {{ path = \"../nalgebra\" }}\n[dev-dependencies]\nnalgebra = {{ path = \"../nalgebra\" }}\n"
        );
        let out = lint("both", LEAF, &helicoid, &[("nalgebra", LEAF)])?;
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(out[0].ends_with("helicoid -> nalgebra"), "{out:?}");
        Ok(())
    }

    #[test]
    fn default_reaching_sweep_or_mint_is_rejected() -> Result<(), String> {
        let helicoid = format!("[features]\n__sweep = []\ndefault = [\"__sweep\"]\n{OK}");
        let out = lint("dsweep", LEAF, &helicoid, &[])?;
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].contains("`helicoid` enables `__sweep` by default"),
            "{out:?}"
        );
        let helicoid = format!(
            "[features]\ndefault = [\"mint\"]\n{OK}mint = {{ path = \"../mint\", optional = true }}\n"
        );
        let out = lint("dmint", LEAF, &helicoid, &[])?;
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].contains("`helicoid` enables `mint` by default"),
            "{out:?}"
        );
        Ok(())
    }

    #[test]
    fn a_stale_lockfile_is_an_error_because_load_is_locked() -> Result<(), String> {
        let crates = [("helicoid", "[dependencies]\n"), ("nalgebra", LEAF)];
        let manifest = workspace("stale", &crates)?;
        let member = manifest.with_file_name("helicoid/Cargo.toml");
        let edited = std::fs::read_to_string(&member)
            .map(|t| format!("{t}nalgebra = {{ path = \"../nalgebra\" }}\n"));
        let written = edited
            .and_then(|t| std::fs::write(&member, t))
            .map_err(|e| e.to_string());
        let loaded = load(&manifest);
        remove(&manifest);
        written?;
        assert!(loaded.is_err(), "a stale Cargo.lock must not be refreshed");
        Ok(())
    }
}
