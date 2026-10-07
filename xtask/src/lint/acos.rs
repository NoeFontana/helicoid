//! Check 8: `acos` is not reached from the group crate (D5).
//!
//! D5 is "`Log` goes through the quaternion `atan2`; no `acos` of a trace anywhere in the
//! workspace", and until [`0022`] it held by **absence**: `Real` had no `acos`, so a group routine
//! that wanted one had to spell `atan2(sqrt(..), x)` and a reviewer would see it. `0022` adds
//! `Real::acos` for `solve_cubic` and `eig3`, which is where it belongs — the argument there is a
//! normalized determinant ratio, not a trace — and in doing so opens the hole D5 names.
//!
//! So the rule is made explicit instead: **`crates/helicoid` does not name `acos`.** That is the
//! crate where rotations live, so every `acos` it could want is one of D5's — of a trace, of a
//! quaternion dot, of anything whose slope is infinite exactly where rotations are interesting.
//! `helicoid-linalg` is exempt: `cubic` and `eig3` are its callers and `float.rs` is where the
//! `libm` routing lives.
//!
//! Test files are exempt for check 5's reason: a test that pins `Real::acos` against `libm::acos`
//! is the twin comparison D6 asks for. A doc comment that names `acos` in prose is not a call, and
//! several do — `so3.rs` explains at length why slerp's `acos` is not taken — so only code is read.
//!
//! [`0022`]: ../../../docs/decisions/0022-real-owes-acos-and-cos.md

use std::collections::BTreeMap;

use super::{comments, is_test, File, Violation};

/// The crate this check guards: D5's subject matter is rotations, which live here.
const GUARDED: &str = "crates/helicoid/";

pub(crate) fn check(files: &[File]) -> Vec<Violation> {
    let mut out = Vec::new();
    for f in files {
        if !f.path.ends_with(".rs") || !f.path.starts_with(GUARDED) || is_test(&f.path) {
            continue;
        }
        let mut prose: BTreeMap<usize, Vec<&str>> = BTreeMap::new();
        for c in comments::of(&f.path, &f.text) {
            prose.entry(c.line).or_default().push(c.text);
        }
        for (i, line) in f.text.lines().enumerate() {
            if !line.contains("acos") {
                continue;
            }
            let mut code = line.to_string();
            for fragment in prose.get(&(i + 1)).into_iter().flatten() {
                code = code.replace(fragment, "");
            }
            if code.contains("acos") {
                out.push(Violation::new(
                    &f.path,
                    i + 1,
                    "[acos] `acos` in the group crate; D5 routes `Log` through the quaternion \
                     `atan2` and `acos`'s slope is infinite where rotations are interesting. \
                     `solve_cubic` and `eig3` are its callers (`0022`)"
                        .to_string(),
                ));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_is_flagged_and_prose_is_not() {
        let hits = |path: &str, text: &str| check(&[File::new(path, text)]).len();
        // A call in the group crate.
        assert_eq!(hits("crates/helicoid/src/so3.rs", "let a = x.acos();"), 1);
        // The same word in prose, which `so3.rs` really does carry.
        assert_eq!(
            hits(
                "crates/helicoid/src/so3.rs",
                "/// `atan2` where slerp has `acos`, whose slope is infinite.\nlet a = 1;"
            ),
            0
        );
        // A code line with a trailing comment that also names it: the code half still counts.
        assert_eq!(
            hits(
                "crates/helicoid/src/so3.rs",
                "let a = x.acos(); // acos here"
            ),
            1
        );
        // `helicoid-linalg` is exempt -- `cubic` and `eig3` are the sanctioned callers.
        assert_eq!(hits("crates/helicoid-linalg/src/cubic.rs", "r.acos()"), 0);
        // Tests are exempt, as check 5's are.
        assert_eq!(hits("crates/helicoid/src/so3_tests.rs", "x.acos()"), 0);
    }
}
