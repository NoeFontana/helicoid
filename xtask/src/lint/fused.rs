//! Check 6: no fused multiply-add in a library crate (D16).
//!
//! D16 forbids a fused multiply-add so that an output is a function of the input bits alone.
//! Nothing enforced it before this check.
//!
//! **The clause is not decorative.** `a * b + c` is *not* contracted into a fused op by rustc at
//! `-O` on either architecture — verified in the instruction stream, `mulsd`/`addsd` on x86_64 and
//! `fmul`/`fadd` on aarch64 — and a scan of the compiled code found no fused op in any of the 4068
//! functions of `helicoid`, `helicoid-linalg` and the `libm` crate built for
//! `aarch64-unknown-linux-gnu`, the architecture where such a scan has teeth, since FMA is in
//! ARMv8's baseline and LLVM emits it wherever the source asks. Asking costs a rounding where the
//! catalogue's error analysis counts two (`docs/maths/coefficients.md` CO.6, CO.9).
//!
//! **Three routes to asking are already closed, and one is not.** `f64::mul_add` is a `std`
//! inherent method and the library crates are `no_std`, so it does not exist there (the mechanism
//! [`super::kernel`] notes for `libm`); `core::intrinsics` and `#[target_feature]` need `unsafe`,
//! which every library root forbids (`0007`); and `Real` declares no fused operation, so generic
//! numeric code cannot call one whatever the scalar is (`0003`) — the compiler refuses it. What is
//! open is **`libm::fma`**, which `libm` 0.2.16 exports and which compiles inside the one module
//! `kernel`'s check lets call `libm::`. That is the route this check exists for; `mul_add` is
//! checked too, because the first library crate to gain `std` would reopen it silently.
//!
//! A `rustflags` or `target-cpu` in a cargo config is not Rust and is not read here;
//! `.cargo/config.toml` says it carries none, and checking that is owed.
//!
//! Tests are exempt, as in [`super::kernel`]: a test that pins a fused op against two roundings is
//! the comparison that shows why the rule exists. A comment naming `mul_add` is not a call, so only
//! code lines are read — the two lines in the tree today are both of that kind.

use std::collections::BTreeMap;

use super::{comments, is_test, File, Violation};

/// How a fused multiply-add can be spelled in a library crate: the `libm` crate's `fma`, `fmaf`
/// and their wider siblings, and `std`'s inherent `mul_add` should a crate ever gain `std`.
const FUSED: [&str; 2] = ["libm::fma", "mul_add"];

pub(crate) fn check(files: &[File]) -> Vec<Violation> {
    let mut out = Vec::new();
    for f in files {
        if !f.path.ends_with(".rs") || !f.path.starts_with("crates/") || is_test(&f.path) {
            continue;
        }
        let mut prose: BTreeMap<usize, Vec<&str>> = BTreeMap::new();
        for c in comments::of(&f.path, &f.text) {
            prose.entry(c.line).or_default().push(c.text);
        }
        for (i, line) in f.text.lines().enumerate() {
            if !FUSED.iter().any(|w| line.contains(w)) {
                continue;
            }
            let mut code = line.to_string();
            for fragment in prose.get(&(i + 1)).into_iter().flatten() {
                if let Some(at) = code.find(*fragment) {
                    code.replace_range(at..at + fragment.len(), "");
                }
            }
            for spelling in FUSED.iter().filter(|w| code.contains(**w)) {
                out.push(Violation::new(
                    &f.path,
                    i + 1,
                    format!(
                        "[fused] `{spelling}` in a library crate; D16 forbids a fused \
                         multiply-add, and `a * b + c` is two roundings on every target"
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

    /// A call in a library crate is a violation; a test, an `xtask` and a comment are not.
    #[test]
    fn a_library_crate_may_not_fuse_a_multiply_and_an_add() {
        let files = [
            // The one route that compiles today: `libm::fma` in the module `kernel`'s check
            // lets call `libm::`.
            File::new(
                "crates/helicoid-linalg/src/float.rs",
                "fn f(a: f64, b: f64, c: f64) -> f64 {\n    libm::fma(a, b, c)\n}\n",
            ),
            File::new(
                "crates/helicoid/src/coeffs/kernel.rs",
                "fn horner(s: f64, z: f64, p: f64) -> f64 {\n    z.mul_add(p, s)\n}\n",
            ),
            File::new(
                "crates/helicoid/src/coeffs/kernel_tests.rs",
                "assert_ne!(z.mul_add(p, s), z * p + s);\n",
            ),
            File::new("xtask/src/seeded/kernel.rs", "let h = z.mul_add(p, s);\n"),
        ];
        let out: Vec<String> = check(&files).iter().map(ToString::to_string).collect();
        assert_eq!(out.len(), 2, "{out:?}");
        assert!(
            out[0].starts_with("crates/helicoid-linalg/src/float.rs:2: [fused] `libm::fma`"),
            "{out:?}"
        );
        assert!(
            out[1].starts_with("crates/helicoid/src/coeffs/kernel.rs:2: [fused] `mul_add`"),
            "{out:?}"
        );
    }

    /// The two lines in the tree today both say the rule rather than break it.
    #[test]
    fn a_comment_that_names_the_rule_is_not_a_call() {
        let files = [
            File::new(
                "crates/helicoid/src/coeffs/kernel.rs",
                "/// `p_j = t_j + z p_{j+1}` (CO.9); no `mul_add`. Under `Dual` the same.\nfn f() {}\n",
            ),
            File::new(
                "crates/helicoid-linalg/src/float.rs",
                "// no mul_add here\nfn g() {}\n",
            ),
        ];
        assert!(check(&files).is_empty(), "{:?}", check(&files));
    }
}
