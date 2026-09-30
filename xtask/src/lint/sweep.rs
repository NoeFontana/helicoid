//! Check 5: no crate other than `xtask` enables the private `__sweep` feature (`0004`,
//! `docs/PHASE1.md` §3, §6), so the sweep's evaluators can never reach a consumer's build. A
//! workspace member breaks the rule by requesting it on a dependency, by naming `dep/__sweep` in
//! its own `[features]`, or by reaching its own `__sweep` from `default`.

use super::metadata::Metadata;
use super::Violation;

const SWEEP: &str = "__sweep";

pub(crate) fn check_metadata(m: &Metadata) -> Vec<Violation> {
    let mut out = Vec::new();
    for lib in m.libraries() {
        for d in lib
            .dependencies
            .iter()
            .filter(|d| d.features.iter().any(|f| f == SWEEP))
        {
            out.push(m.violation(
                lib,
                format!(
                    "[sweep] `{}` enables `{SWEEP}` on dependency `{}`; only `xtask` may",
                    lib.name, d.name
                ),
            ));
        }
        for (feature, members) in &lib.features {
            for member in members.iter().filter(|s| s.ends_with(&format!("/{SWEEP}"))) {
                out.push(m.violation(
                    lib,
                    format!(
                        "[sweep] feature `{feature}` of `{}` enables `{member}`; only `xtask` may",
                        lib.name
                    ),
                ));
            }
        }
        if lib.default_reaches(|f| f == SWEEP) {
            out.push(m.violation(
                lib,
                format!("[sweep] `{}` enables `{SWEEP}` by default", lib.name),
            ));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::super::metadata::fixture::Fx;
    use super::*;

    fn report(m: &Metadata) -> Vec<String> {
        check_metadata(m).iter().map(ToString::to_string).collect()
    }

    #[test]
    fn xtask_may_enable_it_and_a_library_may_define_it() {
        let m = Metadata::fixture(&[
            Fx::member("helicoid").features(&[
                ("__sweep", &[]),
                ("default", &["std"]),
                ("std", &[]),
            ]),
            Fx::member("xtask").deps(&[("helicoid", None, false, &["__sweep"])]),
        ]);
        assert_eq!(report(&m), Vec::<String>::new());
    }

    #[test]
    fn requesting_it_on_a_dependency_is_rejected() {
        let m = Metadata::fixture(&[
            Fx::member("helicoid").deps(&[("helicoid-linalg", None, false, &["__sweep"])]),
            Fx::dep("helicoid-linalg"),
        ]);
        let out = report(&m);
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(out[0].starts_with("helicoid/Cargo.toml:1: [sweep] `helicoid` enables `__sweep` on dependency `helicoid-linalg`"), "{out:?}");
        // A dev-dependency too: it would still build the measured code into a test binary.
        let m = Metadata::fixture(&[Fx::member("helicoid").deps(&[(
            "helicoid-linalg",
            Some("dev"),
            false,
            &["__sweep"],
        )])]);
        assert_eq!(report(&m).len(), 1);
    }

    #[test]
    fn forwarding_it_through_a_feature_is_rejected() {
        let m = Metadata::fixture(&[Fx::member("helicoid").features(&[
            ("fast", &["helicoid-linalg?/__sweep"]),
            ("mint", &["dep:mint"]),
        ])]);
        let out = report(&m);
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].contains("feature `fast` of `helicoid` enables `helicoid-linalg?/__sweep`"),
            "{out:?}"
        );
    }

    #[test]
    fn default_reaching_it_directly_or_transitively_is_rejected() {
        for default in [&["__sweep"][..], &["a"][..]] {
            let m = Metadata::fixture(&[Fx::member("helicoid").features(&[
                ("default", default),
                ("a", &["b"]),
                ("b", &["__sweep"]),
                ("__sweep", &[]),
            ])]);
            let out = report(&m);
            assert_eq!(out.len(), 1, "{out:?}");
            assert!(
                out[0].ends_with("`helicoid` enables `__sweep` by default"),
                "{out:?}"
            );
        }
    }

    #[test]
    fn a_feature_cycle_terminates() {
        let m = Metadata::fixture(&[
            Fx::member("helicoid").features(&[("default", &["a"]), ("a", &["default"])])
        ]);
        assert_eq!(report(&m), Vec::<String>::new());
    }
}
