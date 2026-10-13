//! The resolved `libm` version, read from the workspace lock file into `HELICOID_LIBM_VERSION`,
//! and the target triple into `XTASK_TARGET`.
//!
//! D16's bit identity holds per resolved `libm` version (`docs/decisions/0034`, draft), and `libm`
//! is a caret dependency, so a patch release may move bits. The rows the shipped subject writes
//! and the files the sweep generates name the version they were computed under, so a bump reads as
//! "`libm` moved" rather than as unexplained drift.

use std::path::Path;

fn main() -> Result<(), String> {
    let lock = Path::new(env!("CARGO_MANIFEST_DIR")).join("../Cargo.lock");
    println!("cargo:rerun-if-changed={}", lock.display());
    let text = std::fs::read_to_string(&lock).map_err(|e| format!("{}: {e}", lock.display()))?;
    let version = libm_version(&text)?;
    println!("cargo:rustc-env=HELICOID_LIBM_VERSION={version}");
    // `cargo xtask determinism` names the target its digests were computed on (`PHASE6.md` §3).
    let target = std::env::var("TARGET").map_err(|e| format!("TARGET: {e}"))?;
    println!("cargo:rustc-env=XTASK_TARGET={target}");
    Ok(())
}

/// The version of the one `libm` package in `lock`; two are an error, since the version the
/// kernel linked could then not be named.
fn libm_version(lock: &str) -> Result<String, String> {
    let mut found = Vec::new();
    let mut lines = lock.lines();
    while let Some(line) = lines.next() {
        if line == "name = \"libm\"" {
            let version = lines
                .next()
                .and_then(|l| l.strip_prefix("version = \""))
                .and_then(|l| l.strip_suffix('"'))
                .ok_or("Cargo.lock: `libm` has no version line after its name")?;
            found.push(version.to_string());
        }
    }
    match found.as_slice() {
        [one] => Ok(one.clone()),
        [] => Err("Cargo.lock has no `libm` package".to_string()),
        many => Err(format!(
            "Cargo.lock resolves `libm` more than once: {many:?}"
        )),
    }
}
