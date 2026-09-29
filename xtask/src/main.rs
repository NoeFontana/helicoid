//! `cargo xtask <task>`. Tasks land with the `docs/PHASE1.md` steps that own them; until then a
//! task name is refused rather than silently succeeding.

#![allow(clippy::print_stderr)]

use std::process::ExitCode;

fn main() -> ExitCode {
    match std::env::args().nth(1).as_deref() {
        Some(task) => eprintln!("xtask: `{task}` is not implemented (docs/PHASE1.md §0.0)"),
        None => eprintln!("usage: cargo xtask <task>; no tasks are implemented yet"),
    }
    ExitCode::FAILURE
}
