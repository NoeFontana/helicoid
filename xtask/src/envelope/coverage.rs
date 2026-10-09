//! Coverage (`docs/PHASE1.md` §8): every function id the specs define has a corpus file once its
//! phase has landed. `NUMERICS.md` and `API.md` §3 name routines, not corpus ids; the ids are
//! `PHASE1.md` §4.3's (`required`, one per routine of `NUMERICS.md` §3 to §6, a test equates them
//! with the table) and those the later phases name (`OWED`). An owed id is due when the §0.0 row
//! of its phase that owns it says `Done`, not when the whole table does: those tables carry rows
//! (migrations in other repositories, the 1.0 criteria) that no corpus file waits for.
//!
//! An owed id with a corpus file fails too: it is then required, and stays so, since the corpus
//! never loses a function file (D7). A file counts when the manifest lists it with records and it
//! exists: `MANIFEST.json` alone does not make one.

use super::bars::{Bar, Failure};

/// `(phase, text of the §0.0 row that owns them, where named, ids)`: what later phases add to the
/// corpus. The row text must match exactly one row of its phase's table
/// (a test reads the docs), so a renamed row fails a test, not coverage silently.
const OWED: &[(u8, &str, &str, &[&str])] = &[
    (2, "corpus ids (§6)", "PHASE2.md §6", &["svd3"]),
    (
        5,
        "S² and its chart",
        "PHASE5.md §2",
        &["s2_retract", "s2_local"],
    ),
    (
        5,
        "Sim(3)",
        "PHASE5.md §3",
        &["sim3_exp", "sim3_log", "sim3_jr", "sim3_jr_inv", "sim3_ad"],
    ),
    (5, "Γ₁, Γ₂", "PHASE5.md §4", &["so3_gamma2"]),
    (
        6,
        "ambient Jacobians",
        "PHASE6.md §1",
        &["so3_plus_jacobian", "se3_plus_jacobian"],
    ),
];

/// The ids of `docs/PHASE1.md` §4.3 that corpus v1 holds.
pub(super) fn required() -> Vec<String> {
    let jacobians = ["jr", "jl", "jr_inv", "jl_inv"];
    let group = |prefix: &str, fns: &[&str], suffix: &str| -> Vec<String> {
        let all = fns.iter().chain(&jacobians);
        all.map(|f| format!("{prefix}_{f}{suffix}")).collect()
    };
    // `cos_half` is §4.3's eighth coefficient id (`NUMERICS.md` §3.1), added after this module
    // was written, and `alpha` (`0062`) its ninth; `the_required_ids_are_the_ids_of_the_definitions_table`
    // is what notices.
    let mut ids: Vec<String> = ["k", "a", "b", "c", "d", "e", "r", "cos_half", "alpha"]
        .iter()
        .map(|c| format!("coeff_{c}"))
        .collect();
    ids.extend(group("so3", &["exp", "log", "act", "from_matrix"], ""));
    // `PHASE4.md` §4's two, required since their corpus files landed (`0045` plan step 3): an
    // owed id with a file fails coverage, so the file, this move and the §4.3 rows are one step.
    ids.extend(["so3_geodesic".to_string(), "se3_geodesic".to_string()]);
    // `0060` decision 8: the three SE(3) charts' `retract` and `local`.
    for chart in ["screw", "decoupled", "world"] {
        ids.extend(["retract", "local"].map(|op| format!("se3_{chart}_{op}")));
    }
    for n in 1..=3 {
        ids.extend(group("sen3", &["exp", "log", "ad"], &format!("_n{n}")));
    }
    ids.extend(["so2_exp".to_string(), "so2_log".to_string()]);
    ids.extend(group("se2", &["exp", "log", "ad"], ""));
    // `0056`'s, the routines D7 did not reach: `PHASE2.md` §6's two that exist, `chol` and
    // `chol_solve` at the widths a consumer factors, `renormalize`, and `Dual` as `real_*`.
    ids.extend(
        [
            "solve_cubic",
            "eig3",
            "chol_n3",
            "chol_n6",
            "chol_solve_n3",
            "chol_solve_n6",
            "quat_renormalize",
            "real_sqrt",
            "real_cbrt",
            "real_sin_cos",
            "real_acos",
            "real_atan2",
            "real_div",
        ]
        .map(String::from),
    );
    // `0061`'s, with its corpus file in the same step.
    ids.push("mat2_inverse_adj".to_string());
    ids
}

/// The ids a corpus must hold, for a test that builds one.
#[cfg(test)]
pub(super) fn ids_for_tests() -> Vec<String> {
    required()
}

pub(super) struct Coverage {
    pub(super) failures: Vec<Failure>,
    pub(super) required: usize,
    pub(super) excused: usize,
}

/// `ids` are the function files of the corpus; `landed(phase, row)` is whether that row of the
/// phase's §0.0 says `Done`, or why it cannot be told.
pub(super) fn check(ids: &[String], landed: &dyn Fn(u8, &str) -> Result<bool, String>) -> Coverage {
    let has = |id: &str| ids.iter().any(|i| i == id);
    let fail = |text: String| Failure::new(Bar::Coverage, text);
    let mut c = Coverage {
        failures: Vec::new(),
        required: required().len(),
        excused: 0,
    };
    for id in required().iter().filter(|id| !has(id)) {
        c.failures
            .push(fail(format!("`{id}` has no corpus file (PHASE1.md §4.3)")));
    }
    for &(phase, row, source, owed) in OWED {
        let due = landed(phase, row);
        if let Err(why) = &due {
            c.failures.push(fail(why.clone()));
        }
        for id in owed {
            match (has(id), &due) {
                (true, _) => c.failures.push(fail(format!(
                    "`{id}` has a corpus file and is listed as owed: move it to the required ids"
                ))),
                (false, Ok(true)) => c.failures.push(fail(format!(
                    "phase {phase}'s row `{row}` is Done and `{id}` ({source}) has no corpus file"
                ))),
                (false, Ok(false)) => c.excused += 1,
                (false, Err(_)) => {}
            }
        }
    }
    c
}

/// `(area, status)` of each data row of the table under `## 0.0`. The first cell is the area and
/// the rest of the row is the status, since a status may hold a `|` in code.
fn status_rows(text: &str) -> Vec<(String, String)> {
    let (mut rows, mut inside) = (Vec::new(), false);
    for line in text.lines() {
        if line.starts_with("## ") {
            inside = line.starts_with("## 0.0 ");
        } else if let (true, Some((area, status))) = (
            inside,
            line.strip_prefix('|').and_then(|l| l.split_once('|')),
        ) {
            let status = status.trim_end_matches('|').trim().trim_start_matches('*');
            if !(status.starts_with("Status") || status.starts_with(['-', ':'])) {
                rows.push((area.trim().to_string(), status.to_string()));
            }
        }
    }
    rows
}

/// Whether the one row of `docs/PHASE<phase>.md`'s §0.0 whose area contains `row` says `Done`.
/// No such file, no such row and several are errors: a coverage that cannot tell must not excuse.
pub(super) fn landed_in(root: &std::path::Path, phase: u8, row: &str) -> Result<bool, String> {
    let file = format!("docs/PHASE{phase}.md");
    let text = std::fs::read_to_string(root.join(&file)).map_err(|e| format!("{file}: {e}"))?;
    let rows = status_rows(&text);
    let mut hits = rows.iter().filter(|(area, _)| area.contains(row));
    match (hits.next(), hits.next()) {
        (Some((_, status)), None) => Ok(status.starts_with("Done")),
        (None, _) => Err(format!(
            "{file} §0.0 has no row containing `{row}`: coverage cannot tell whether its ids are due"
        )),
        (Some(_), Some(_)) => Err(format!(
            "{file} §0.0 has several rows containing `{row}`: coverage cannot tell which owns its ids"
        )),
    }
}

/// The docs a checkout needs for coverage to read them: per phase a §0.0 table with one row per
/// owed group, `Done` for those named in `done` and `Not started` for the rest.
#[cfg(test)]
pub(super) fn docs_for_tests(done: &[&str]) -> Vec<(String, String)> {
    let mut phases: Vec<u8> = OWED.iter().map(|o| o.0).collect();
    phases.dedup();
    let table = |phase: u8| {
        let rows = OWED.iter().filter(|o| o.0 == phase).map(|o| {
            let status = if done.contains(&o.1) {
                "Done"
            } else {
                "Not started"
            };
            format!("| {} (§1) | {status} |\n", o.1)
        });
        format!(
            "## 0.0 Implementation status\n\n| Area | Status |\n|---|---|\n{}",
            rows.collect::<String>()
        )
    };
    phases
        .into_iter()
        .map(|p| (format!("docs/PHASE{p}.md"), table(p)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus;
    use crate::conformance::testkit::Scratch;

    fn messages(c: &Coverage) -> Vec<&str> {
        c.failures.iter().map(|f| f.text.as_str()).collect()
    }

    /// `landed` for a check that says the rows named in `done` are Done and the others are not.
    fn rows(done: &'static [&'static str]) -> impl Fn(u8, &str) -> Result<bool, String> {
        move |_, row| Ok(done.contains(&row))
    }

    #[test]
    fn the_committed_corpus_covers_what_is_required_and_holds_nothing_owed() -> Result<(), String> {
        let dir = crate::conformance::corpus_dir()?;
        let have: Vec<String> = corpus::manifest(&dir)?
            .into_iter()
            .map(|e| e.fn_id)
            .collect();
        let c = check(&have, &rows(&[]));
        assert!(c.failures.is_empty(), "{:?}", messages(&c));
        assert_eq!((c.required, c.excused), (required().len(), 11));
        Ok(())
    }

    #[test]
    fn the_required_ids_are_the_ids_of_the_definitions_table() -> Result<(), String> {
        let path = crate::conformance::root()?.join("docs/PHASE1.md");
        let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
        let (mut table, mut inside) = (Vec::new(), false);
        for line in text.lines() {
            if line.starts_with('#') {
                inside = line.starts_with("### 4.3 ");
            } else if inside && line.starts_with("| `") {
                let first = line.split('|').nth(1).unwrap_or_default();
                for span in first.split('`').skip(1).step_by(2) {
                    match span.split_once('{') {
                        Some((stem, n)) => {
                            let n = n.trim_end_matches('}').split(',');
                            table.extend(n.map(|n| format!("{stem}{n}")));
                        }
                        // The series file is not a corpus file (its manifest `kind`).
                        None if span == "coeff_series" => {}
                        None => table.push(span.to_string()),
                    }
                }
            }
        }
        table.sort();
        let mut have = required();
        have.sort();
        assert_eq!(table, have);
        Ok(())
    }

    #[test]
    fn a_missing_required_id_fails_and_an_owed_id_is_excused_until_its_row_is_done() {
        let mut have = required();
        have.retain(|i| i != "sen3_jr_inv_n2");
        let c = check(&have, &rows(&[]));
        assert_eq!(
            messages(&c),
            ["`sen3_jr_inv_n2` has no corpus file (PHASE1.md §4.3)"]
        );
        // One row Done: its ids fail, and no other group's.
        let c = check(&required(), &rows(&["corpus ids (§6)"]));
        assert_eq!(c.failures.len(), 1);
        assert!(messages(&c)[0].starts_with(
            "phase 2's row `corpus ids (§6)` is Done and `svd3` (PHASE2.md §6) has no corpus file"
        ));
        assert_eq!(c.excused, 10);
    }

    #[test]
    fn every_owed_group_is_due_by_its_own_row_and_by_no_other() {
        for &(phase, row, _, ids) in OWED {
            let only = |p: u8, r: &str| -> Result<bool, String> { Ok(p == phase && r == row) };
            let c = check(&required(), &only);
            assert_eq!(c.failures.len(), ids.len(), "{row}");
            assert_eq!(c.excused, 11 - ids.len(), "{row}");
        }
    }

    #[test]
    fn a_row_that_cannot_be_read_fails_once_per_group_instead_of_excusing_it() {
        let unreadable = |_: u8, row: &str| Err(format!("no row `{row}`"));
        let c = check(&required(), &unreadable);
        assert_eq!(c.failures.len(), OWED.len());
        assert_eq!(c.excused, 0);
    }

    #[test]
    fn an_owed_id_that_has_a_file_fails() {
        let mut have = required();
        have.push("s2_retract".into());
        let c = check(&have, &rows(&[]));
        assert_eq!(c.failures.len(), 1);
        assert!(c.failures[0]
            .text
            .contains("`s2_retract` has a corpus file and is listed as owed"));
        have.push("s2_local".into());
        assert_eq!(check(&have, &rows(&[])).failures.len(), 2);
    }

    #[test]
    fn a_row_is_done_when_its_own_status_starts_with_done_whatever_the_rest_says(
    ) -> Result<(), String> {
        let doc = |rows: &str| {
            format!("# t\n\n## 0.0 Implementation status\n\n| Area | Status |\n|---|---|\n{rows}\n## 1. Next\n\n| x | Not started |\n")
        };
        let read = |rows: &str| status_rows(&doc(rows));
        let got = read("| a (§1) | Done |\n| b | **Done** |\n| c | Partial: `x | y` |\n");
        let status: Vec<&str> = got.iter().map(|r| r.1.as_str()).collect();
        assert_eq!(status, ["Done", "Done**", "Partial: `x | y`"]);
        assert_eq!(got[0].0, "a (§1)");
        assert!(read("").is_empty() && status_rows("| a | Done |\n").is_empty());

        let scratch = Scratch::new("landed");
        let write = |name: &str, text: String| {
            let dir = scratch.0.join("docs");
            std::fs::create_dir_all(&dir).and_then(|()| std::fs::write(dir.join(name), text))
        };
        let table = "| own row (§1) | Done |\n| unrelated | Not started |\n| another | Not Done yet |\n| word | Partial: Done in part |\n";
        assert!(write("PHASE2.md", doc(table)).is_ok());
        let of = |row: &str| landed_in(&scratch.0, 2, row);
        // An unrelated row that is not Done does not hold the phase back, and "Done" inside a
        // status that does not start with it is not Done.
        assert_eq!(of("own row"), Ok(true));
        assert_eq!(of("unrelated"), Ok(false));
        assert_eq!((of("another"), of("word")), (Ok(false), Ok(false)));
        let no_row = of("nothing").err().unwrap_or_default();
        assert!(no_row.contains("no row containing `nothing`"), "{no_row}");
        // `o` is in three areas.
        assert!(of("o").err().unwrap_or_default().contains("several rows"));
        let no_file = landed_in(&scratch.0, 3, "own row")
            .err()
            .unwrap_or_default();
        assert!(no_file.starts_with("docs/PHASE3.md"), "{no_file}");
        Ok(())
    }

    #[test]
    fn every_owed_row_is_one_row_of_its_phases_table_in_this_checkout() -> Result<(), String> {
        let root = crate::conformance::root()?;
        for &(phase, row, _, _) in OWED {
            landed_in(&root, phase, row)?;
        }
        Ok(())
    }
}
