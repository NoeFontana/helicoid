//! Check 4: the normal-dependency closure of each library crate equals its allowed set
//! (`0007`, `docs/PHASE1.md` §3). Normal edges only, all features, every target platform; dev and
//! build edges are not part of what a consumer links. A crate outside the set is reported once,
//! at the edge where the closure leaves the allowed set, with the path that brings it and the
//! crates hidden behind it, so adding `nalgebra` is one line and not forty.
//!
//! Membership is by package name: a path or git crate that is named `libm` passes here, and
//! `cargo deny check sources` (`just audit`) is what rejects a non-crates.io registry or git
//! source. The resolve graph is workspace-wide (see `metadata`), so the report can name a crate
//! that only `xtask`'s features bring in.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use super::metadata::{Metadata, Package};
use super::Violation;

/// `(library crate, required closure)`: `0007` item 1. A workspace member other than `xtask`
/// that is absent here is a violation (fail closed), so a new crate registers its budget, and the
/// record that authorises it, in its own PR.
const BUDGET: &[(&str, &[&str])] = &[
    ("helicoid-linalg", &["libm"]),
    ("helicoid", &["helicoid-linalg", "libm"]),
];

/// Allowed in the closure only as an optional dependency declared by a `BUDGET` crate.
const OPTIONAL: &[&str] = &["mint"];

/// How many hidden crates a violation names before saying "and N more".
const NAMED: usize = 6;

pub(crate) fn check_metadata(m: &Metadata) -> Vec<Violation> {
    let mut out = Vec::new();
    for lib in m.libraries() {
        match BUDGET.iter().find(|(n, _)| *n == lib.name) {
            Some((_, required)) => check_crate(m, lib, required, &mut out),
            None => out.push(m.violation(
                lib,
                format!(
                    "[closure] `{}` has no dependency budget; register it in `closure::BUDGET` in the PR that adds the crate",
                    lib.name
                ),
            )),
        }
    }
    out
}

fn check_crate<'m>(m: &'m Metadata, lib: &'m Package, required: &[&str], out: &mut Vec<Violation>) {
    for d in lib
        .dependencies
        .iter()
        .filter(|d| d.kind.is_none() && OPTIONAL.contains(&d.name.as_str()))
    {
        if !d.optional {
            out.push(m.violation(
                lib,
                format!(
                    "[closure] `{}` must be an optional dependency of `{}`",
                    d.name, lib.name
                ),
            ));
        } else if lib.default_reaches(|f| {
            let key = d.key();
            f == key
                || f.strip_prefix("dep:") == Some(key)
                || f.strip_prefix(key).is_some_and(|r| r.starts_with('/'))
        }) {
            out.push(m.violation(
                lib,
                format!(
                    "[closure] `{}` enables `{}` by default; it is a feature only",
                    lib.name, d.name
                ),
            ));
        }
    }
    let name_of = |id: &'m str| m.package(id).map_or(id, |p| p.name.as_str());
    let mut parent: BTreeMap<&str, &str> = BTreeMap::new();
    let mut queue = VecDeque::from([lib.id.as_str()]);
    let mut order: Vec<(&str, &str)> = Vec::new();
    while let Some(id) = queue.pop_front() {
        let Some(node) = m.resolve.nodes.iter().find(|n| n.id == id) else {
            continue;
        };
        let mut deps: Vec<&str> = node
            .deps
            .iter()
            .filter(|d| d.dep_kinds.iter().any(|k| k.kind.is_none()))
            .map(|d| d.pkg.as_str())
            .collect();
        deps.sort_by_key(|d| name_of(d));
        for dep in deps {
            if dep != lib.id && !parent.contains_key(dep) {
                parent.insert(dep, id);
                order.push((dep, id));
                queue.push_back(dep);
            }
        }
    }
    let allowed = |id: &'m str, from: &'m str| {
        let (n, from) = (name_of(id), name_of(from));
        required.contains(&n) || (OPTIONAL.contains(&n) && BUDGET.iter().any(|(b, _)| *b == from))
    };
    let offenders: BTreeSet<&str> = order
        .iter()
        .filter(|&&(id, from)| !allowed(id, from))
        .map(|&(id, _)| id)
        .collect();
    let mut hidden: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for &id in &offenders {
        let mut top = id;
        while let Some(&p) = parent.get(top) {
            if !offenders.contains(p) {
                break;
            }
            top = p;
        }
        if top != id {
            hidden.entry(top).or_default().push(name_of(id));
        } else {
            hidden.entry(top).or_default();
        }
    }
    for (top, behind) in hidden {
        let mut path = vec![name_of(top)];
        let mut at = top;
        while let Some(&p) = parent.get(at) {
            path.push(name_of(p));
            at = p;
        }
        path.reverse();
        let mut msg = format!(
            "[closure] `{}` is outside the allowed dependency set (`0007`): {}",
            name_of(top),
            path.join(" -> ")
        );
        if !behind.is_empty() {
            let more = behind.len().saturating_sub(NAMED);
            msg.push_str(&format!(
                "; brings in {}",
                behind[..behind.len() - more].join(", ")
            ));
            if more > 0 {
                msg.push_str(&format!(" and {more} more"));
            }
        }
        out.push(m.violation(lib, msg));
    }
    let present: BTreeSet<&str> = order.iter().map(|&(id, _)| name_of(id)).collect();
    for r in required.iter().filter(|r| !present.contains(*r)) {
        out.push(m.violation(
            lib,
            format!("[closure] `{}` no longer depends on `{r}`, which its closure must contain (`0007`)", lib.name),
        ));
    }
}

#[cfg(test)]
mod tests {
    use super::super::metadata::fixture::{Fx, D};
    use super::*;

    const fn n(name: &str) -> D<'_> {
        (name, None, false, &[])
    }
    const fn opt(name: &str) -> D<'_> {
        (name, None, true, &[])
    }
    const fn kind<'a>(name: &'a str, kind: &'a str) -> D<'a> {
        (name, Some(kind), false, &[])
    }

    fn report(m: &Metadata) -> Vec<String> {
        check_metadata(m).iter().map(ToString::to_string).collect()
    }

    /// The workspace of `0007` with the given dependency lists, plus `others` (leaf packages,
    /// extra members).
    fn tree(linalg: &[D<'_>], helicoid: &[D<'_>], others: &[Fx<'_>]) -> Metadata {
        let mut pkgs = vec![
            Fx::member("helicoid-linalg").deps(linalg),
            Fx::member("helicoid").deps(helicoid),
            Fx::dep("libm"),
            Fx::dep("mint"),
        ];
        pkgs.extend_from_slice(others);
        Metadata::fixture(&pkgs)
    }

    const LINALG: [D<'_>; 2] = [n("libm"), opt("mint")];
    const HELICOID: [D<'_>; 3] = [n("helicoid-linalg"), n("libm"), opt("mint")];

    #[test]
    fn the_budget_passes_and_xtask_is_unrestricted() {
        let deps = [n("serde"), n("nalgebra")];
        let xtask = Fx::member("xtask").deps(&deps);
        let m = tree(
            &LINALG,
            &HELICOID,
            &[xtask, Fx::dep("serde"), Fx::dep("nalgebra")],
        );
        assert_eq!(report(&m), Vec::<String>::new());
    }

    #[test]
    fn planted_nalgebra_is_rejected_with_its_path() {
        let helicoid = [n("helicoid-linalg"), n("libm"), n("nalgebra")];
        let out = report(&tree(&LINALG, &helicoid, &[Fx::dep("nalgebra")]));
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].starts_with("helicoid/Cargo.toml:1: [closure] `nalgebra`"),
            "{out:?}"
        );
        assert!(out[0].ends_with("helicoid -> nalgebra"), "{out:?}");
    }

    #[test]
    fn a_dependency_of_a_dependency_names_the_whole_path() {
        let linalg = [n("libm"), n("nalgebra")];
        let out = report(&tree(&linalg, &HELICOID, &[Fx::dep("nalgebra")]));
        assert_eq!(out.len(), 2, "{out:?}");
        assert!(out[0].ends_with("helicoid-linalg -> nalgebra"), "{out:?}");
        assert!(
            out[1].ends_with("helicoid -> helicoid-linalg -> nalgebra"),
            "{out:?}"
        );
    }

    #[test]
    fn planted_serde_through_a_transitive_dependency() {
        let serde = [n("serde")];
        let foo = Fx::dep("foo").deps(&serde);
        let helicoid = [n("helicoid-linalg"), n("libm"), n("foo")];
        let out = report(&tree(&LINALG, &helicoid, &[foo, Fx::dep("serde")]));
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].contains("`foo`") && out[0].contains("helicoid -> foo"),
            "{out:?}"
        );
        assert!(out[0].ends_with("; brings in serde"), "{out:?}");
        // Behind an allowed crate the offender is itself the frontier.
        let libm = Fx::dep("libm").deps(&serde);
        let out = report(&Metadata::fixture(&[
            Fx::member("helicoid-linalg").deps(&LINALG),
            Fx::member("helicoid").deps(&HELICOID),
            libm,
            Fx::dep("mint"),
            Fx::dep("serde"),
        ]));
        assert_eq!(out.len(), 2, "{out:?}");
        assert!(
            out[0].ends_with("helicoid-linalg -> libm -> serde"),
            "{out:?}"
        );
        assert!(out[1].ends_with("helicoid -> libm -> serde"), "{out:?}");
    }

    #[test]
    fn hidden_crates_are_summarised() {
        let names = ["a", "b", "c", "d", "e", "f", "g", "h"];
        let deps: Vec<D<'_>> = names.iter().map(|x| n(x)).collect();
        let big = Fx::dep("big").deps(&deps);
        let helicoid = [n("helicoid-linalg"), n("libm"), n("big")];
        let leaves: Vec<Fx<'_>> = names.iter().map(|x| Fx::dep(x)).collect();
        let mut others = vec![big];
        others.extend(leaves);
        let out = report(&tree(&LINALG, &helicoid, &others));
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].ends_with("; brings in a, b, c, d, e, f and 2 more"),
            "{out:?}"
        );
    }

    #[test]
    fn dev_and_build_edges_are_not_part_of_the_closure() {
        let helicoid = [
            n("helicoid-linalg"),
            n("libm"),
            kind("nalgebra", "dev"),
            kind("serde", "build"),
        ];
        let others = [Fx::dep("nalgebra"), Fx::dep("serde")];
        assert_eq!(
            report(&tree(&LINALG, &helicoid, &others)),
            Vec::<String>::new()
        );
    }

    #[test]
    fn mint_is_allowed_only_as_an_optional_direct_dependency() {
        let helicoid = [n("helicoid-linalg"), n("libm"), n("mint")];
        let out = report(&tree(&LINALG, &helicoid, &[]));
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].contains("`mint` must be an optional dependency of `helicoid`"),
            "{out:?}"
        );
        let mint = [n("mint")];
        let libm = Fx::dep("libm").deps(&mint);
        let out = report(&Metadata::fixture(&[
            Fx::member("helicoid-linalg").deps(&[n("libm")]),
            Fx::member("helicoid").deps(&[n("helicoid-linalg"), n("libm")]),
            libm,
            Fx::dep("mint"),
        ]));
        assert_eq!(out.len(), 2, "{out:?}");
        assert!(
            out[0].ends_with("helicoid-linalg -> libm -> mint"),
            "{out:?}"
        );
        assert!(out[1].ends_with("helicoid -> libm -> mint"), "{out:?}");
    }

    #[test]
    fn a_chain_is_reported_once_at_its_frontier() {
        let (b, c) = ([n("b")], [n("c")]);
        let others = [Fx::dep("a").deps(&b), Fx::dep("b").deps(&c), Fx::dep("c")];
        let helicoid = [n("helicoid-linalg"), n("libm"), n("a")];
        let out = report(&tree(&LINALG, &helicoid, &others));
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(out[0].ends_with("helicoid -> a; brings in b, c"), "{out:?}");
    }

    #[test]
    fn the_reported_path_is_the_shortest() {
        // `off` is reachable in two hops through `a` and in three through `z`; a search that
        // follows the longer branch first would file it under `z`.
        let (off, zq, q) = ([n("off")], [n("q")], [n("off")]);
        let others = [
            Fx::dep("a").deps(&off),
            Fx::dep("z").deps(&zq),
            Fx::dep("q").deps(&q),
            Fx::dep("off"),
        ];
        let helicoid = [n("helicoid-linalg"), n("libm"), n("a"), n("z")];
        let out = report(&tree(&LINALG, &helicoid, &others));
        assert_eq!(out.len(), 2, "{out:?}");
        assert!(out[0].ends_with("helicoid -> a; brings in off"), "{out:?}");
        assert!(out[1].ends_with("helicoid -> z; brings in q"), "{out:?}");
    }

    #[test]
    fn mint_is_a_feature_never_a_default() {
        let features: &[(&str, &[&str])] = &[("default", &["mint"])];
        let m = Metadata::fixture(&[
            Fx::member("helicoid-linalg").deps(&LINALG),
            Fx::member("helicoid").deps(&HELICOID).features(features),
            Fx::dep("libm"),
            Fx::dep("mint"),
        ]);
        let out = report(&m);
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].contains("`helicoid` enables `mint` by default"),
            "{out:?}"
        );
        let features: &[(&str, &[&str])] = &[("default", &["std"]), ("std", &["dep:mint"])];
        let m = Metadata::fixture(&[
            Fx::member("helicoid-linalg")
                .deps(&LINALG)
                .features(features),
            Fx::member("helicoid").deps(&HELICOID),
            Fx::dep("libm"),
            Fx::dep("mint"),
        ]);
        let out = report(&m);
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].contains("`helicoid-linalg` enables `mint` by default"),
            "{out:?}"
        );
        // A weak forward does not enable the dependency.
        let features: &[(&str, &[&str])] = &[("default", &["mint?/serde"])];
        let m = Metadata::fixture(&[
            Fx::member("helicoid-linalg").deps(&LINALG),
            Fx::member("helicoid").deps(&HELICOID).features(features),
            Fx::dep("libm"),
            Fx::dep("mint"),
        ]);
        assert_eq!(report(&m), Vec::<String>::new());
    }

    #[test]
    fn mint_must_be_optional_in_both_crates_but_a_dev_dependency_is_free() {
        let linalg = [n("libm"), n("mint")];
        let out = report(&tree(&linalg, &HELICOID, &[]));
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].contains("`mint` must be an optional dependency of `helicoid-linalg`"),
            "{out:?}"
        );
        let helicoid = [n("helicoid-linalg"), n("libm"), kind("mint", "dev")];
        assert_eq!(report(&tree(&LINALG, &helicoid, &[])), Vec::<String>::new());
    }

    #[test]
    fn a_missing_required_dependency_and_an_unregistered_member_are_violations() {
        let out = report(&tree(&[opt("mint")], &HELICOID, &[]));
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].contains("helicoid-linalg` no longer depends on `libm`"),
            "{out:?}"
        );
        let out = report(&tree(&LINALG, &HELICOID, &[Fx::member("helicoid-spline")]));
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].contains("`helicoid-spline` has no dependency budget"),
            "{out:?}"
        );
    }
}
