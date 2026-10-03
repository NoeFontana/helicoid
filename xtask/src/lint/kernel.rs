//! Check 5: `libm::` appears only in the one private kernel module (D16).
//!
//! D16 routes every transcendental through the `libm` crate so that an output is a function of the
//! input bits alone. That holds today by construction: `helicoid-linalg`'s `float.rs` is one macro
//! `impl Real for f32/f64` and holds every call — `sqrt`, `cbrt`, `sincos`, `atan2`, `fabs`,
//! `copysign` — while `Dual` inherits them through `Real` and has none of its own.
//!
//! `no_std` enforces half of that: with no `std`, `f64::sin` as an inherent method does not exist,
//! so a second `libm` call site cannot appear by *forgetting* the trait. Nothing stopped one being
//! added on purpose, and a transcendental reached outside the kernel is a bit-identity claim nobody
//! measured — the envelope would not catch it, because the envelope compares against a reference
//! that routes the same way.
//!
//! Test files are exempt: a test that pins `Real::atan2` against `libm::atan2` directly is the twin
//! comparison D6 asks for, and the whole measurement in `0032` (draft) was swapping one call for
//! another.
//! A doc comment that names `libm::atan2` in prose is not a call, so only code lines are read.

use std::collections::BTreeMap;

use super::{comments, File, Violation};

/// The one module that may call `libm`, and the crate-relative paths that are tests.
const KERNEL: &str = "crates/helicoid-linalg/src/float.rs";

/// Whether `path` is a test file, by the repository's own `*_tests.rs` / `tests.rs` convention.
fn is_test(path: &str) -> bool {
    path.ends_with("_tests.rs") || path.ends_with("/tests.rs")
}

pub(crate) fn check(files: &[File]) -> Vec<Violation> {
    let mut out = Vec::new();
    for f in files {
        if !f.path.ends_with(".rs") || f.path == KERNEL || is_test(&f.path) {
            continue;
        }
        if !f.path.starts_with("crates/") {
            continue;
        }
        // A line's comment fragments are removed before the test, so prose that names `libm::sin`
        // is not a call. `cubic.rs` has two such lines today.
        let mut prose: BTreeMap<usize, Vec<&str>> = BTreeMap::new();
        for c in comments::of(&f.path, &f.text) {
            prose.entry(c.line).or_default().push(c.text);
        }
        for (i, line) in f.text.lines().enumerate() {
            if !line.contains("libm::") {
                continue;
            }
            let mut code = line.to_string();
            for fragment in prose.get(&(i + 1)).into_iter().flatten() {
                if let Some(at) = code.find(*fragment) {
                    code.replace_range(at..at + fragment.len(), "");
                }
            }
            if code.contains("libm::") {
                out.push(Violation::new(
                    &f.path,
                    i + 1,
                    format!(
                        "[kernel] `libm::` outside `{KERNEL}`; D16 routes every transcendental \
                         through one private kernel, reached by `Real`"
                    ),
                ));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A call outside the kernel is a violation; the kernel itself and a test are not.
    #[test]
    fn only_the_kernel_may_call_libm() {
        let files = [
            File::new(KERNEL, "fn sqrt(self) -> Self { libm::sqrt(self) }\n"),
            File::new(
                "crates/helicoid/src/quat.rs",
                "fn bad(x: f64) -> f64 {\n    libm::atan2(x, 1.0)\n}\n",
            ),
            File::new(
                "crates/helicoid/src/quat_tests.rs",
                "assert_eq!(a, libm::atan2(x, y));\n",
            ),
            File::new(
                "crates/helicoid-linalg/src/tests.rs",
                "let want = libm::sin(x);\n",
            ),
            File::new("xtask/src/oracle.rs", "let y = libm::sin(x);\n"),
        ];
        let out: Vec<String> = check(&files).iter().map(ToString::to_string).collect();
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(
            out[0].starts_with("crates/helicoid/src/quat.rs:2: [kernel]"),
            "{out:?}"
        );
    }

    /// Prose that names `libm::atan2` is not a call: `cubic.rs` has two such lines today.
    #[test]
    fn a_doc_comment_naming_libm_is_not_a_call() {
        let files = [File::new(
            "crates/helicoid-linalg/src/cubic.rs",
            "/// and `libm::atan2` is the one that costs.\n// libm::sin(x) would be a call\nfn f() {}\n",
        )];
        assert!(check(&files).is_empty(), "{:?}", check(&files));
    }
}
