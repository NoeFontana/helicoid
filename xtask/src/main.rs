//! `cargo xtask <task>`. Tasks land with the `docs/PHASE1.md` steps that own them; an unknown
//! task name is refused rather than silently succeeding.

#![allow(clippy::print_stderr)]

mod conformance;
mod envelope;
mod lint;
mod seeded;
mod shipped;
mod thresholds;

use std::process::ExitCode;

/// Every implemented task; the usage line and the unknown-task error read this list.
const TASKS: &[&str] = &["lint", "conformance", "thresholds", "envelope"];

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if dispatch(
        &args,
        lint::run,
        conformance::run,
        thresholds::run,
        envelope::run,
    ) {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

/// Runs the task named by `args[0]` with the rest as its arguments, reporting on stderr; `true`
/// only when it ran and found nothing to complain about. The tasks are injected so the exit policy
/// tests without touching the disk.
fn dispatch(
    args: &[String],
    lint: impl FnOnce() -> Result<Vec<lint::Violation>, String>,
    conformance: impl FnOnce(&[String]) -> Result<(), String>,
    thresholds: impl FnOnce(&[String]) -> Result<(), String>,
    envelope: impl FnOnce(&[String]) -> Result<(), String>,
) -> bool {
    match args.first().map(String::as_str) {
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
        Some("conformance") => match conformance(&args[1..]) {
            Ok(()) => true,
            Err(e) => {
                eprintln!("xtask conformance: {e}");
                false
            }
        },
        Some("thresholds") => match thresholds(&args[1..]) {
            Ok(()) => true,
            Err(e) => {
                eprintln!("xtask thresholds: {e}");
                false
            }
        },
        Some("envelope") => match envelope(&args[1..]) {
            Ok(()) => true,
            Err(e) => {
                eprintln!("xtask envelope: {e}");
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
        let run = |args: &[&str],
                   lint: Result<Vec<lint::Violation>, String>,
                   conf: Result<(), String>| {
            let args: Vec<String> = args.iter().map(ToString::to_string).collect();
            dispatch(
                &args,
                || lint,
                |_| conf,
                |_| Err("unused".into()),
                |_| Err("unused".into()),
            )
        };
        assert!(run(&["lint"], Ok(vec![]), Ok(())));
        assert!(!run(&["lint"], Ok(vec![violation()]), Ok(())));
        assert!(!run(&["lint"], Err("boom".into()), Ok(())));
        assert!(!run(&["nope"], Ok(vec![]), Ok(())));
        assert!(!run(&[], Ok(vec![]), Ok(())));
        assert!(run(
            &["conformance", "--fn", "x"],
            Err("unused".into()),
            Ok(())
        ));
        assert!(!run(&["conformance"], Ok(vec![]), Err("nonfinite".into())));
        // `thresholds` gets the arguments after its name, and its `Err` fails the run.
        let thr = |args: &[&str]| {
            let args: Vec<String> = args.iter().map(ToString::to_string).collect();
            let check = |rest: &[String]| match *rest == ["--check"] {
                true => Ok(()),
                false => Err("usage".to_string()),
            };
            dispatch(
                &args,
                || Ok(vec![]),
                |_| Ok(()),
                check,
                |_| Err("unused".into()),
            )
        };
        assert!(thr(&["thresholds", "--check"]) && !thr(&["thresholds"]));
        // `envelope` gets the arguments after its name, and its `Err` fails the run.
        let env = |args: &[&str], out: Result<(), String>| {
            let args: Vec<String> = args.iter().map(ToString::to_string).collect();
            let unused = |_: &[String]| Err("unused".to_string());
            dispatch(
                &args,
                || Ok(vec![]),
                unused,
                unused,
                |rest| {
                    assert_eq!(rest, ["--check"]);
                    out
                },
            )
        };
        assert!(env(&["envelope", "--check"], Ok(())));
        assert!(!env(&["envelope", "--check"], Err("1 failure(s)".into())));
    }

    #[test]
    fn conformance_receives_the_arguments_after_its_name() {
        let args: Vec<String> = ["conformance", "--fn", "so2_exp"]
            .map(String::from)
            .to_vec();
        let seen = std::cell::RefCell::new(Vec::new());
        let ok = dispatch(
            &args,
            || Ok(vec![]),
            |rest| {
                seen.borrow_mut().extend(rest.iter().cloned());
                Ok(())
            },
            |_| Err("unused".into()),
            |_| Err("unused".into()),
        );
        assert!(ok && *seen.borrow() == ["--fn", "so2_exp"]);
    }
}
