//! `cargo xtask thresholds [--check]` (`docs/PHASE1.md` §6, `docs/decisions/0004`, `0016` item 3):
//! for each of `k, a, b, c, d, e`, `cos θ/2` and `r` at `f64` and at `f32`, the series terms and
//! switch that minimise the maximum, over every `theta:*` record of `coeff_<c>` (the `@f32`
//! strata at `f32`), of the larger of the value error and the `d/dz` error through `Dual<S, 1>`,
//! both exact and in units of `u` of the precision (`conformance::metric`, the `coeff_*` rule of
//! one field each). Writes `conformance/sweeps/thresholds-seeded.csv` and, from it,
//! `xtask/src/seeded/generated.rs` (`emit`); `--check` compares a fresh sweep with both, byte for
//! byte, and writes nothing. The sweep measures the seeded kernel at a candidate it is handed,
//! never at `generated.rs`, so that file is an output and not an input; the seeded correct kernel
//! then runs it, which is why a hand edit that does not compile stops the tool building (`emit`).
//! Phase 3 points the tool at `helicoid::coeffs`, whose own files are `conformance/sweeps/
//! thresholds.csv` and `crates/helicoid/src/coeffs/generated.rs`.
//!
//! **Candidates.** `1..=8` terms times `grid`: 64 points per decade of the branch variable
//! `z = θ²` from `1e-16` (`θ = 1e-8`) to `1` (`θ = 1`), ends included, 8200 candidates, at each
//! precision the points correctly rounded there (`grid`); the per-decade-of-`z` reading is
//! 0014 (draft) question 8. A candidate uses the series arm where `z < switch` (`z = fl(θ·θ)` at the
//! precision, as `Seeded::eval` forms it) and the exact arm elsewhere; each arm's error is formed
//! once per record and a candidate only selects. Ties, on the exact `f64` objective, go to fewer
//! terms, then to the larger switch. The prior, `tf_tree` D12 (four terms, `θ < 0.1`: `z < 0.01`,
//! on the grid at `f64`, rounded to binary32 at `f32`), is scored as a named candidate, for every
//! id but `r` though `NUMERICS.md` §4 lists it for `a`, `b`, `c` only. `r` has no D12 value in its
//! branch variable: its prior is the same (four terms, `0.01`) taken in `s`, which is `θ = 0.2`
//! where D12's `θ < 0.1` is `s < 2.5e-3`, so its `prior_*` cells and the figure printed beside it
//! are that candidate's, not D12's. An error that is not finite is `inf`, so a candidate selecting
//! it never wins.
//!
//! **`r` is a reading, pending the maintainer, not a spec** (0014 (draft) question 29, which
//! holds the reading and its differences from 0015 (draft) NU.6): the branch variable is
//! `s = n²/w²`, formed by a division, the grid is over `s`, only a record with `w > 0` enters the
//! objective, and a candidate takes the series arm iff `w > 0` and `s < switch` (`w = 0` selects
//! the exact arm since `s` is infinite, and `w < 0`, which only S² charts produce, always does).
//! The `Dual` is seeded on `n²` at fixed `w`, the derivative the corpus stores.
//!
//! **Where a row is not a measured optimum.** `terms` stops at §6's `m ≤ 8`. The grid stops at
//! `θ = 1` and `z < switch` is strict, so the corpus record at `θ = 1` exactly is always on the
//! exact arm and no switch above it is swept: a top of `nextUp(1)` gives `e` 8 terms and an
//! objective 3.3 times lower at `f64` (0014 (draft) question 9). `theta:dense` (200 points per
//! decade of `θ`) meets the `f64` grid at every 16th index, where a record's arm is the last bit
//! of `fl(θ·θ)` against the rounded grid point, so quote no switch beyond `below_objective` and
//! `above_objective`.
//!
//! **Columns**, one row per coefficient and precision (`f64` rows first), a function of the corpus
//! and the kernel:
//!
//! | column | meaning |
//! |---|---|
//! | `coeff`, `precision` | `k, a, b, c, d, e, cos_half, r`; `f64` or `f32` |
//! | `terms` | series terms of the chosen candidate |
//! | `switch_bits`, `switch_z`, `switch_theta` | its switch: the bit pattern `generated.rs` emits (16 hex digits at `f64`, 8 at `f32`), its shortest decimal at the precision, and its square root (`n/w` for `r`) |
//! | `grid_index` | its place in the grid, `0..=1024`; `1024` is the top, `θ = 1` |
//! | `value_max_u`, `deriv_max_u` | its maxima over the records of the value error and of the `d/dz` error |
//! | `objective` | their larger |
//! | `argmax_field`, `argmax_stratum`, `argmax_id`, `argmax_z` | what attains `objective`: `value` or `deriv` (`value` on a tie), and the first such record's stratum, id and branch variable |
//! | `below_objective`, `above_objective` | the chosen terms one grid point below and above the chosen switch; empty at an end |
//! | `top_objective` | the best objective among the candidates at the top of the grid |
//! | `prior_objective` | the same for the prior: D12, and for `r` the candidate above |
//! | `prior_rank` | 1 + the grid candidates whose objective is smaller than the prior's |
//! | `tied` | grid candidates with exactly the chosen objective, among which the tie-break chose |
//! | `next_objective` | the smallest grid objective above the chosen one; empty when none |
//!
//! Numbers are shortest round-trip decimals (`{:e}`), those of an `f32` row's switch and branch
//! variable at `f32`. Not swept: SE(2)'s `α`, `β` (no corpus id), a switch shared by a call-site
//! group (CO.18), and `crates/helicoid/src/coeffs/generated.rs` (Phase 3).

mod emit;
mod grid;
mod measure;
mod search;

use std::path::Path;

use helicoid_linalg::{Dual, Precision, Real};

use crate::conformance::root;
use crate::seeded::{d12, Coeff, Series, Swept, D1, SERIES_FILE};
use search::{Field, Sweep};

const USAGE: &str = "usage: cargo xtask thresholds [--check]";
/// The committed sweep of the seeded kernels, relative to the repository root. The name leaves
/// `conformance/sweeps/thresholds.csv` to `helicoid::coeffs` (Phase 3).
const CSV: &str = "conformance/sweeps/thresholds-seeded.csv";
const HEADER: &str = "coeff,precision,terms,switch_bits,switch_z,switch_theta,grid_index,\
value_max_u,deriv_max_u,objective,argmax_field,argmax_stratum,argmax_id,argmax_z,below_objective,\
above_objective,top_objective,prior_objective,prior_rank,tied,next_objective";

/// One coefficient's sweep at one precision and the record that attains its objective.
struct Row {
    id: Swept,
    precision: Precision,
    sweep: Sweep,
    field: Field,
    stratum: String,
    record: u64,
    z: f64,
}

/// The sweep of `ids` over the corpus in `dir` at the precision of `S`.
fn sweep_at<S: Real + Into<f64>>(dir: &Path, ids: &[Swept]) -> Result<Vec<Row>, String> {
    let series = Series::<Dual<S, 1>>::load(dir)?;
    let prior = d12(&series)?;
    let grid = grid::grid(S::PRECISION);
    // The grid is the precision's own. On this corpus a binary64 grid swept at `f32` selects the
    // same arms, so only this refusal shows the wiring.
    if let Some(g) = grid.iter().find(|&&g| !grid::holds(S::PRECISION, g)) {
        let p = precision_name(S::PRECISION);
        return Err(format!("grid point {g:e} is not a value {p} holds"));
    }
    ids.iter()
        .map(|&id| {
            let m = measure::samples(dir, &series, id)?;
            let prior = (prior.terms, prior.switch_z.v.into());
            let sweep = search::search(&m.samples, &grid, prior)?;
            let (field, at) = sweep.chosen.argmax();
            let (stratum, record) = m.records[at].clone();
            let z = m.samples[at].z;
            Ok(Row {
                id,
                precision: S::PRECISION,
                sweep,
                field,
                stratum,
                record,
                z,
            })
        })
        .collect()
}

/// The sweep of `ids` at `f64`, then at `f32`.
fn sweep(dir: &Path, ids: &[Swept]) -> Result<Vec<Row>, String> {
    let mut rows = sweep_at::<f64>(dir, ids)?;
    rows.extend(sweep_at::<f32>(dir, ids)?);
    Ok(rows)
}

/// The name of `precision` in the CSV.
fn precision_name(precision: Precision) -> &'static str {
    match precision {
        Precision::F64 => "f64",
        Precision::F32 => "f32",
    }
}

/// `x`, a value `precision` holds, as the shortest decimal that reads back at `precision`.
fn shown(precision: Precision, x: f64) -> String {
    match precision {
        Precision::F64 => format!("{x:e}"),
        Precision::F32 => format!("{:e}", x as f32),
    }
}

/// The square root of `z` at `precision`.
fn root_of(precision: Precision, z: f64) -> f64 {
    match precision {
        Precision::F64 => libm::sqrt(z),
        Precision::F32 => f64::from(libm::sqrtf(z as f32)),
    }
}

/// The CSV of `rows`.
fn render(rows: &[Row]) -> String {
    let mut out = format!("{HEADER}\n");
    let opt = |x: Option<f64>| x.map_or(String::new(), |x| format!("{x:e}"));
    for r in rows {
        let (s, z, p) = (&r.sweep, r.sweep.chosen.switch, r.precision);
        let bits = match p {
            Precision::F64 => format!("0x{:016x}", z.to_bits()),
            Precision::F32 => format!("0x{:08x}", (z as f32).to_bits()),
        };
        let cells = [
            r.id.name().to_string(),
            precision_name(p).to_string(),
            s.chosen.terms.to_string(),
            bits,
            shown(p, z),
            shown(p, root_of(p, z)),
            s.index.to_string(),
            format!("{:e}", s.chosen.value),
            format!("{:e}", s.chosen.deriv),
            format!("{:e}", s.chosen.objective()),
            r.field.name().to_string(),
            r.stratum.clone(),
            r.record.to_string(),
            shown(p, r.z),
            opt(s.below),
            opt(s.above),
            format!("{:e}", s.top),
            format!("{:e}", s.prior.objective()),
            s.prior_rank.to_string(),
            s.tied.to_string(),
            opt(s.next),
        ];
        out.push_str(&cells.join(","));
        out.push('\n');
    }
    out
}

/// The one coefficient's sweep at `f64`, kept to rank a candidate against its choice: what
/// `conformance --self-test` asks of the planted `c` (`docs/PHASE1.md` §10).
pub(crate) struct Ranker {
    samples: Vec<search::Sample>,
    grid: Vec<f64>,
    chosen: f64,
}

/// Where a candidate stands among the candidates swept for its coefficient.
pub(crate) struct Ranking {
    pub(crate) objective: f64,
    /// The chosen candidate's objective.
    pub(crate) chosen: f64,
    /// 1 + the swept candidates whose objective is strictly smaller.
    pub(crate) rank: usize,
    /// The candidates swept: `TERMS` times the grid.
    pub(crate) of: usize,
}

impl Ranker {
    pub(crate) fn new(dir: &Path, c: Coeff) -> Result<Self, String> {
        let series = Series::<D1>::load(dir)?;
        let prior = d12(&series)?;
        let grid = grid::grid(Precision::F64);
        let samples = measure::samples(dir, &series, Swept::Coeff(c))?.samples;
        let sweep = search::search(&samples, &grid, (prior.terms, prior.switch_z.v))?;
        Ok(Self {
            samples,
            grid,
            chosen: sweep.chosen.objective(),
        })
    }

    /// `terms` series terms below `switch` (`z`), scored like a swept candidate; `switch` need not
    /// be on the grid.
    pub(crate) fn rank(&self, terms: usize, switch: f64) -> Result<Ranking, String> {
        if !(1..=search::TERMS).contains(&terms) {
            return Err(format!(
                "{terms} terms are not in the swept 1..={}",
                search::TERMS
            ));
        }
        let objective = search::score(&self.samples, terms, switch).objective();
        Ok(Ranking {
            objective,
            chosen: self.chosen,
            rank: search::rank_of(&self.samples, &self.grid, objective),
            of: search::TERMS * self.grid.len(),
        })
    }
}

/// The files a sweep writes, as `(path relative to the repository root, text)`: the CSV and the
/// generated file, both computed before either is written.
fn files(dir: &Path, rows: &[Row]) -> Result<Vec<(&'static str, String)>, String> {
    let csv = render(rows);
    let path = dir.join(SERIES_FILE);
    let series_file = std::fs::read(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let (wide, narrow) = (Series::<f64>::load(dir)?, Series::<f32>::load(dir)?);
    let generated = emit::render(&csv, &series_file, &wide, &narrow)?;
    Ok(vec![(CSV, csv), (emit::PATH, generated)])
}

/// `Ok` when `fresh` is what `path` (repository-relative) holds, else the first line that differs.
fn compare(path: &str, committed: &str, fresh: &str) -> Result<(), String> {
    if committed == fresh {
        return Ok(());
    }
    let mut lines = committed.lines().zip(fresh.lines()).enumerate();
    let at = lines.find(|(_, (a, b))| a != b).map_or_else(
        || "one is a prefix of the other".to_string(),
        |(i, (a, b))| format!("line {}: committed `{a}`, fresh `{b}`", i + 1),
    );
    Err(format!(
        "{path} is not what a fresh sweep writes ({at}); run `just thresholds` and review the diff"
    ))
}

/// With `check`, `Ok` when `root/path` holds `text` and nothing is written; otherwise writes `text`.
fn finish(root: &Path, path: &str, text: &str, check: bool) -> Result<(), String> {
    let full = root.join(path);
    if check {
        let committed =
            std::fs::read_to_string(&full).map_err(|e| format!("{}: {e}", full.display()))?;
        return compare(path, &committed, text);
    }
    let dir = full.parent().ok_or("no parent dir")?;
    std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    std::fs::write(&full, text).map_err(|e| format!("{}: {e}", full.display()))
}

/// [`finish`] for every file; a check names every file that differs, not the first.
fn finish_all(root: &Path, files: &[(&str, String)], check: bool) -> Result<(), String> {
    let failed: Vec<String> = files
        .iter()
        .filter_map(|(path, text)| finish(root, path, text, check).err())
        .collect();
    match failed.is_empty() {
        true => Ok(()),
        false => Err(failed.join("; ")),
    }
}

pub(crate) fn run(args: &[String]) -> Result<(), String> {
    let root = root()?;
    run_at(&root, &root.join("conformance/corpus"), args)
}

/// [`run`] with the files under `root` and the corpus in `corpus`.
#[allow(clippy::print_stdout)]
fn run_at(root: &Path, corpus: &Path, args: &[String]) -> Result<(), String> {
    let check = match args {
        [] => false,
        [flag] if flag == "--check" => true,
        _ => return Err(USAGE.to_string()),
    };
    let rows = sweep(corpus, &Swept::ALL)?;
    let files = files(corpus, &rows)?;
    finish_all(root, &files, check)?;
    if !check {
        for r in &rows {
            println!(
                "{} {}: {} terms, switch root {} (grid {}), objective {:e} u; prior {:e} u",
                r.id.name(),
                precision_name(r.precision),
                r.sweep.chosen.terms,
                shown(r.precision, root_of(r.precision, r.sweep.chosen.switch)),
                r.sweep.index,
                r.sweep.chosen.objective(),
                r.sweep.prior.objective()
            );
        }
        for (path, _) in &files {
            eprintln!("thresholds: wrote {}", root.join(path).display());
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;
    use crate::conformance::testkit::Scratch;
    use crate::seeded::C_PLANTED;

    fn committed() -> Result<String, String> {
        let path = root()?.join(CSV);
        std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))
    }

    /// The files of a fresh sweep of every id at both precisions.
    fn fresh() -> Result<Vec<(&'static str, String)>, String> {
        let dir = corpus_dir()?;
        files(&dir, &sweep(&dir, &Swept::ALL)?)
    }

    #[test]
    fn arguments_are_refused_and_a_difference_names_its_file_and_line() {
        assert!(run(&["--nope".to_string()]).is_err() && run(&["a".into(), "b".into()]).is_err());
        assert!(compare("p", "a\nb\n", "a\nb\n").is_ok());
        let e = compare("x/p.rs", "a\nb\n", "a\nc\n")
            .err()
            .unwrap_or_default();
        assert!(
            e.starts_with("x/p.rs is not what a fresh sweep writes"),
            "{e}"
        );
        assert!(e.contains("line 2: committed `b`, fresh `c`"), "{e}");
        assert!(compare("p", "a\n", "a\nb\n").is_err());
    }

    #[test]
    fn the_committed_files_are_what_a_fresh_sweep_writes() -> Result<(), String> {
        let fresh = fresh()?;
        let paths: Vec<&str> = fresh.iter().map(|f| f.0).collect();
        assert_eq!(
            paths,
            [
                "conformance/sweeps/thresholds-seeded.csv",
                "xtask/src/seeded/generated.rs"
            ]
        );
        finish_all(&root()?, &fresh, true)
    }

    #[test]
    fn two_sweeps_write_the_same_bytes() -> Result<(), String> {
        let two = [Swept::Coeff(Coeff::B), Swept::Coeff(Coeff::E)];
        let run = || sweep(&corpus_dir()?, &two).map(|s| render(&s));
        let (first, second) = (run()?, run()?);
        assert_eq!(first, second);
        // The header, then `b`, `e` at `f64`, then at `f32`.
        assert_eq!(first.lines().count(), 5);
        Ok(())
    }

    #[test]
    fn the_committed_rows_are_the_documented_columns() -> Result<(), String> {
        let text = committed()?;
        let mut rows = text.lines();
        assert_eq!(rows.next(), Some(HEADER));
        let columns: Vec<&str> = HEADER.split(',').collect();
        let mut names = Vec::new();
        for row in rows {
            let cells: Vec<&str> = row.split(',').collect();
            assert_eq!(cells.len(), columns.len(), "{row}");
            let cell = |name: &str| {
                let i = columns.iter().position(|c| *c == name).ok_or(name)?;
                Ok::<_, String>(cells[i])
            };
            let num = |name: &str| {
                let x = cell(name)?;
                x.parse::<f64>()
                    .map_err(|e| format!("{row}: {name} `{x}`: {e}"))
            };
            // The objective is the larger maximum, and no grid candidate beats it, D12 included.
            let objective = num("objective")?;
            assert_eq!(
                num("value_max_u")?.max(num("deriv_max_u")?).to_bits(),
                objective.to_bits()
            );
            let (top, prior) = (num("top_objective")?, num("prior_objective")?);
            assert!(objective <= top && objective <= prior && objective < num("next_objective")?);
            // The named field is the one that attains it, and its record is a `theta:*` one.
            let field = format!("{}_max_u", cell("argmax_field")?);
            assert_eq!(num(&field)?.to_bits(), objective.to_bits(), "{row}");
            assert!(cell("argmax_stratum")?.starts_with("theta:"), "{row}");
            // A neighbour is no better than the chosen switch, and is empty only at an end.
            let end = cell("grid_index")? == "1024";
            assert_eq!(cell("above_objective")?.is_empty(), end, "{row}");
            for n in ["below_objective", "above_objective"] {
                assert!(cell(n)?.is_empty() || num(n)? >= objective, "{row}");
            }
            // The switch is 16 hex digits of bits at `f64` and 8 at `f32`, the row's own precision.
            let width = if cell("precision")? == "f64" { 18 } else { 10 };
            assert_eq!(cell("switch_bits")?.len(), width, "{row}");
            names.push((cells[0], cells[1]));
        }
        let ids = ["k", "a", "b", "c", "d", "e", "cos_half", "r"];
        let want = ["f64", "f32"].map(|p| ids.map(|c| (c, p)));
        assert_eq!(names, want.concat());
        Ok(())
    }

    #[test]
    fn a_check_refuses_a_drifted_file_and_writes_nothing_while_a_sweep_writes_it(
    ) -> Result<(), String> {
        let scratch = Scratch::new("finish");
        let path = scratch.0.join("sweeps/thresholds.csv");
        let read = || std::fs::read_to_string(&path).map_err(|e| e.to_string());
        let (root, rel) = (scratch.0.as_path(), "sweeps/thresholds.csv");
        // No file: a check is an error and creates nothing.
        assert!(finish(root, rel, "a\nb\n", true).is_err() && !path.exists());
        // A write creates the directory and the file; a check of the same text passes.
        finish(root, rel, "a\nb\n", false)?;
        assert_eq!(read()?, "a\nb\n");
        finish(root, rel, "a\nb\n", true)?;
        // Drift, including a missing trailing newline, is an error naming the line, and the
        // committed file is left as it was.
        let e = finish(root, rel, "a\nc\n", true).err().unwrap_or_default();
        assert!(e.contains("line 2: committed `b`, fresh `c`"), "{e}");
        assert!(finish(root, rel, "a\nb", true).is_err());
        assert_eq!(read()?, "a\nb\n");
        // A write replaces it.
        finish(root, rel, "a\nc\n", false)?;
        assert_eq!(read()?, "a\nc\n");
        Ok(())
    }

    /// `text` with the character at `at` replaced by a different one.
    fn one_character_off(text: &str, at: usize) -> String {
        let swap = |c: char| match c {
            'x' => 'y',
            _ => 'x',
        };
        let mut chars: Vec<char> = text.chars().collect();
        chars[at] = swap(chars[at]);
        chars.into_iter().collect()
    }

    #[test]
    fn a_one_character_hand_edit_of_a_generated_file_fails_the_check() -> Result<(), String> {
        // `docs/PHASE1.md` §11: the sweep regenerates the seeded `generated.rs` byte for byte, and
        // `thresholds --check` fails on a hand edit. Every 7th character of each file, the first
        // and the last, and a hex digit of a switch and a digit of a literal, which still parse.
        let fresh = fresh()?;
        let scratch = Scratch::new("hand-edit");
        finish_all(&scratch.0, &fresh, false)?;
        finish_all(&scratch.0, &fresh, true)?;
        for (rel, text) in &fresh {
            let full = scratch.0.join(rel);
            let n = text.chars().count();
            let mut edits: Vec<String> = (0..n)
                .step_by(7)
                .map(|i| one_character_off(text, i))
                .collect();
            edits.push(one_character_off(text, n - 1));
            edits.push(text.replacen("0x3f", "0x3e", 1));
            edits.push(text.replacen("e-1,", "e-2,", 1));
            for edited in edits.iter().filter(|e| *e != text) {
                std::fs::write(&full, edited).map_err(|e| e.to_string())?;
                // The edit differs in one character, or in one digit through `replacen`.
                let e = finish_all(&scratch.0, &fresh, true).err();
                let e = e.ok_or_else(|| format!("{rel}: an edit went unseen"))?;
                assert!(e.contains(rel) && e.contains("line "), "{e}");
                // A check writes nothing: the edit is still there.
                let kept = std::fs::read_to_string(&full).map_err(|e| e.to_string())?;
                assert_eq!(&kept, edited);
            }
            assert!(edits.len() > 100, "{rel}: {} edits", edits.len());
            std::fs::write(&full, text).map_err(|e| e.to_string())?;
            finish_all(&scratch.0, &fresh, true)?;
        }
        Ok(())
    }

    #[test]
    fn the_generated_file_is_a_function_of_the_sweep_it_is_checked_with() -> Result<(), String> {
        // Only the CSV drifts: the check names it and not the generated file, and the other way.
        let fresh = fresh()?;
        let scratch = Scratch::new("which-file");
        finish_all(&scratch.0, &fresh, false)?;
        for (rel, text) in &fresh {
            std::fs::write(scratch.0.join(rel), text.replacen('e', "E", 1))
                .map_err(|e| e.to_string())?;
            let e = finish_all(&scratch.0, &fresh, true)
                .err()
                .unwrap_or_default();
            // The one file that drifted is named; the other's path is not the subject of an error.
            assert!(e.starts_with(rel), "{e}");
            assert_eq!(
                e.matches("is not what a fresh sweep writes").count(),
                1,
                "{e}"
            );
            std::fs::write(scratch.0.join(rel), text).map_err(|e| e.to_string())?;
        }
        Ok(())
    }

    #[test]
    fn a_run_writes_both_files_and_its_check_fails_a_hand_edit_and_writes_nothing(
    ) -> Result<(), String> {
        // `run_at` is what `just thresholds` and `just thresholds-check` reach, under a scratch root.
        let (corpus, scratch) = (corpus_dir()?, Scratch::new("run-at"));
        let (at, check) = (scratch.0.as_path(), ["--check".to_string()]);
        let read = |dir: &Path, path: &str| {
            std::fs::read_to_string(dir.join(path)).map_err(|e| format!("{path}: {e}"))
        };
        run_at(at, &corpus, &[])?;
        for path in [CSV, emit::PATH] {
            assert_eq!(read(at, path)?, read(&root()?, path)?, "{path}");
        }
        run_at(at, &corpus, &check)?;
        // A hand edit fails the check, which names the file and the line and leaves the edit.
        let edited = read(at, emit::PATH)?.replacen("0x3f", "0x3e", 1);
        std::fs::write(at.join(emit::PATH), &edited).map_err(|e| e.to_string())?;
        let e = run_at(at, &corpus, &check).err();
        let e = e.ok_or("the check passed a hand edit")?;
        assert!(e.starts_with(emit::PATH) && e.contains("line "), "{e}");
        assert_eq!(read(at, emit::PATH)?, edited);
        Ok(())
    }

    #[test]
    fn the_planted_c_of_section_10_is_dominated() -> Result<(), String> {
        // `c` with switch 1e-8 (the first grid point) and two terms: the exact arm's `12u/θ²`
        // above it, so the sweep ranks it far below the chosen candidate.
        let ranker = Ranker::new(&corpus_dir()?, Coeff::C)?;
        let (terms, switch) = C_PLANTED;
        assert_eq!(switch.to_bits(), grid::grid(Precision::F64)[0].to_bits());
        let r = ranker.rank(terms, switch)?;
        assert!(
            r.objective > 1e6 * r.chosen,
            "{:e} vs {:e}",
            r.objective,
            r.chosen
        );
        assert_eq!(r.of, 8200);
        // Only candidates that also leave nearly every record to the exact arm rank with it.
        assert!(r.rank * 10 > r.of * 9, "{} of {}", r.rank, r.of);
        Ok(())
    }

    #[test]
    fn a_candidate_is_scored_at_the_terms_it_is_given() -> Result<(), String> {
        // Below a switch on the series arm every record is `terms` terms of it: fewer terms score
        // worse and rank lower. (At the planted `c`'s switch of `1e-16` they do not differ.)
        let ranker = Ranker::new(&corpus_dir()?, Coeff::C)?;
        let at = |terms| ranker.rank(terms, 1e-2);
        let (one, two, eight) = (at(1)?, at(2)?, at(8)?);
        assert!(one.objective > two.objective && two.objective > eight.objective);
        assert!(one.rank > two.rank && two.rank > eight.rank);
        Ok(())
    }

    #[test]
    fn the_chosen_candidate_ranks_first_and_a_candidate_off_the_grid_is_scored(
    ) -> Result<(), String> {
        let dir = corpus_dir()?;
        let ranker = Ranker::new(&dir, Coeff::A)?;
        let chosen = committed()?;
        let row = chosen
            .lines()
            .find(|l| l.starts_with("a,f64,"))
            .ok_or("no row for a")?;
        let cells: Vec<&str> = row.split(',').collect();
        let bits = u64::from_str_radix(cells[3].trim_start_matches("0x"), 16)
            .map_err(|e| e.to_string())?;
        let terms: usize = cells[2].parse().map_err(|e| format!("{e}"))?;
        let r = ranker.rank(terms, f64::from_bits(bits))?;
        assert_eq!((r.rank, r.objective.to_bits()), (1, r.chosen.to_bits()));
        // A switch off the grid is scored like any other, and 0 or 9 terms are refused.
        let off = ranker.rank(4, 0.01)?;
        assert!(off.objective > r.chosen);
        assert!(ranker.rank(0, 0.01).is_err() && ranker.rank(9, 0.01).is_err());
        Ok(())
    }
}
