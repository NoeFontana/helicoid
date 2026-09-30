//! `cargo xtask <task>`. Tasks land with the `docs/PHASE1.md` steps that own them; an unknown
//! task name is refused rather than silently succeeding.

#![allow(clippy::print_stderr)]

mod lint;

use std::process::ExitCode;

/// Every implemented task; the usage line and the unknown-task error read this list.
const TASKS: &[&str] = &["lint"];

fn main() -> ExitCode {
    let task = std::env::args().nth(1);
    if dispatch(task.as_deref(), lint::run) {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

/// Runs `task`, reporting on stderr; `true` only when it ran and found nothing to complain about.
/// `lint` is injected so the exit policy tests without touching the disk.
fn dispatch(
    task: Option<&str>,
    lint: impl FnOnce() -> Result<Vec<lint::Violation>, String>,
) -> bool {
    match task {
        Some("lint") => match lint() {
            Ok(violations) if violations.is_empty() => true,
            Ok(violations) => {
                for v in &violations {
                    eprintln!("{v}");
                }
                eprintln!("xtask lint: {} violation(s)", violations.len());
                false
            }
            Err(e) => {
                eprintln!("xtask lint: {e}");
                false
            }
        },
        Some(task) => {
            eprintln!("xtask: `{task}` is not implemented (docs/PHASE1.md §0.0); have: {TASKS:?}");
            false
        }
        None => {
            eprintln!("usage: cargo xtask <task>; tasks: {TASKS:?}");
            false
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exit_policy() {
        let violation = || lint::Violation::new("a.md", 1, "bad");
        assert!(dispatch(Some("lint"), || Ok(vec![])));
        assert!(!dispatch(Some("lint"), || Ok(vec![violation()])));
        assert!(!dispatch(Some("lint"), || Err("boom".into())));
        assert!(!dispatch(Some("nope"), || Ok(vec![])));
        assert!(!dispatch(None, || Ok(vec![])));
    }
}
