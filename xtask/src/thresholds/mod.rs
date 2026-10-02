//! `cargo xtask thresholds [--check]` (`docs/PHASE1.md` §6, `docs/decisions/0004`): for each of
//! `k, a, b, c, d, e` at `f64`, the series terms and switch that minimise the maximum, over every
//! `theta:*` record of `coeff_<c>`, of the larger of the value error and the `d/dz` error through
//! `Dual<f64, 1>`, both exact and in units of `u` (`conformance::metric`, the `coeff_*` rule of one
//! field each). Writes `conformance/sweeps/thresholds.csv`; `--check` compares a fresh sweep with
//! it. The subject is the seeded kernel; Phase 3 points the tool at `helicoid::coeffs`.
//!
//! **Candidates.** `1..=8` terms times `grid`: 64 points per decade of the branch variable
//! `z = θ²` from `1e-16` (`θ = 1e-8`) to `1` (`θ = 1`), ends included, 8200 candidates; the
//! per-decade-of-`z` reading is 0014 (draft) question 8. A candidate uses the series arm where
//! `z < switch` (`z = fl(θ·θ)`, as `Seeded::eval` forms it) and the exact arm elsewhere; each
//! arm's error is formed once per record and a candidate only selects. Ties, on the exact `f64`
//! objective, go to fewer terms, then to the larger switch. The prior, `tf_tree` D12 (four terms,
//! `θ < 0.1`: `z < 0.01`, on the grid), is scored as a named candidate, for all six though
//! `NUMERICS.md` §4 lists it for `a`, `b`, `c`. An error that is not finite is `inf`, so a
//! candidate selecting it never wins.
//!
//! **Where a row is not a measured optimum.** `terms` stops at §6's `m ≤ 8`. The grid stops at
//! `θ = 1` and `z < switch` is strict, so the corpus record at `θ = 1` exactly is always on the
//! exact arm and no switch above it is swept: a top of `nextUp(1)` gives `e` 8 terms and an
//! objective 3.3 times lower (0014 (draft) question 9). `theta:dense` (200 points per decade of `θ`) meets the
//! grid at every 16th index, where a record's arm is the last bit of `fl(θ·θ)` against the
//! rounded grid point; `c`'s switch is one of them, so quote no switch beyond `below_objective`
//! and `above_objective`.
//!
//! **Columns**, one row per coefficient, a function of the corpus and the kernel:
//!
//! | column | meaning |
//! |---|---|
//! | `coeff`, `precision` | `k, a, b, c, d, e`; `f64` |
//! | `terms` | series terms of the chosen candidate |
//! | `switch_bits`, `switch_z`, `switch_theta` | its switch: the bit pattern `generated.rs` will emit, its shortest decimal, and `sqrt(z)` |
//! | `grid_index` | its place in the grid, `0..=1024`; `1024` is the top, `θ = 1` |
//! | `value_max_u`, `deriv_max_u` | its maxima over the records of the value error and of the `d/dz` error |
//! | `objective` | their larger |
//! | `argmax_field`, `argmax_stratum`, `argmax_id`, `argmax_z` | what attains `objective`: `value` or `deriv` (`value` on a tie), and the first such record's stratum, id and `z` |
//! | `below_objective`, `above_objective` | the chosen terms one grid point below and above the chosen switch; empty at an end |
//! | `top_objective` | the best objective among the candidates at the top of the grid |
//! | `prior_objective` | the same for the D12 prior |
//! | `prior_rank` | 1 + the grid candidates whose objective is smaller than the prior's |
//! | `tied` | grid candidates with exactly the chosen objective, among which the tie-break chose |
//! | `next_objective` | the smallest grid objective above the chosen one; empty when none |
//!
//! Numbers are shortest round-trip decimals (`{:e}`). Not swept: `r` (its branch variable
//! `n²/w²` and `w ≤ 0` domain are open, `docs/maths/coefficients.md` CO.16), `f32` (no
//! `f32`-exact corpus inputs), SE(2)'s `α`, `β` and `cos θ/2` (no corpus id), a switch shared by a
//! call-site group (CO.18), and `generated.rs`.

mod grid;
mod measure;
mod search;

use std::path::Path;

use crate::conformance::root;
use crate::seeded::{d12, Coeff, Series, D1};
use search::{Field, Sweep};

const USAGE: &str = "usage: cargo xtask thresholds [--check]";
/// The committed sweep, relative to the repository root.
const CSV: &str = "conformance/sweeps/thresholds.csv";
const HEADER: &str = "coeff,precision,terms,switch_bits,switch_z,switch_theta,grid_index,\
value_max_u,deriv_max_u,objective,argmax_field,argmax_stratum,argmax_id,argmax_z,below_objective,\
above_objective,top_objective,prior_objective,prior_rank,tied,next_objective";

/// One coefficient's sweep and the record that attains its objective.
struct Row {
    coeff: Coeff,
    sweep: Sweep,
    field: Field,
    stratum: String,
    id: u64,
    z: f64,
}

/// The sweep of `coeffs` over the corpus in `dir`.
fn sweep(dir: &Path, coeffs: &[Coeff]) -> Result<Vec<Row>, String> {
    let series = Series::<D1>::load(dir)?;
    let prior = d12(&series)?;
    let grid = grid::grid();
    coeffs
        .iter()
        .map(|&coeff| {
            let m = measure::samples(dir, &series, coeff)?;
            let sweep = search::search(&m.samples, &grid, (prior.terms, prior.switch_z.v))?;
            let (field, at) = sweep.chosen.argmax();
            let (stratum, id) = m.records[at].clone();
            let z = m.samples[at].z;
            Ok(Row {
                coeff,
                sweep,
                field,
                stratum,
                id,
                z,
            })
        })
        .collect()
}

/// The CSV of `rows`.
fn render(rows: &[Row]) -> String {
    let mut out = format!("{HEADER}\n");
    let opt = |x: Option<f64>| x.map_or(String::new(), |x| format!("{x:e}"));
    for r in rows {
        let (s, z) = (&r.sweep, r.sweep.chosen.switch);
        let cells = [
            r.coeff.name().to_string(),
            "f64".to_string(),
            s.chosen.terms.to_string(),
            format!("0x{:016x}", z.to_bits()),
            format!("{z:e}"),
            format!("{:e}", libm::sqrt(z)),
            s.index.to_string(),
            format!("{:e}", s.chosen.value),
            format!("{:e}", s.chosen.deriv),
            format!("{:e}", s.chosen.objective()),
            r.field.name().to_string(),
            r.stratum.clone(),
            r.id.to_string(),
            format!("{:e}", r.z),
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

/// `Ok` when `fresh` is `committed`, else the first line that differs.
fn compare(committed: &str, fresh: &str) -> Result<(), String> {
    if committed == fresh {
        return Ok(());
    }
    let mut lines = committed.lines().zip(fresh.lines()).enumerate();
    let at = lines.find(|(_, (a, b))| a != b).map_or_else(
        || "one is a prefix of the other".to_string(),
        |(i, (a, b))| format!("line {}: committed `{a}`, fresh `{b}`", i + 1),
    );
    Err(format!(
        "{CSV} is not what a fresh sweep writes ({at}); run `just thresholds` and review the diff"
    ))
}

/// With `check`, `Ok` when `path` holds `text` and nothing is written; otherwise writes `text`.
fn finish(path: &Path, text: &str, check: bool) -> Result<(), String> {
    if check {
        let committed =
            std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        return compare(&committed, text);
    }
    let dir = path.parent().ok_or("no parent dir")?;
    std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    std::fs::write(path, text).map_err(|e| format!("{}: {e}", path.display()))
}

#[allow(clippy::print_stdout)]
pub(crate) fn run(args: &[String]) -> Result<(), String> {
    let check = match args {
        [] => false,
        [flag] if flag == "--check" => true,
        _ => return Err(USAGE.to_string()),
    };
    let root = root()?;
    let rows = sweep(&root.join("conformance/corpus"), &Coeff::ALL)?;
    let path = root.join(CSV);
    finish(&path, &render(&rows), check)?;
    if !check {
        for r in &rows {
            println!(
                "{}: {} terms, switch theta {:e} (grid {}), objective {:e} u; D12 prior {:e} u",
                r.coeff.name(),
                r.sweep.chosen.terms,
                libm::sqrt(r.sweep.chosen.switch),
                r.sweep.index,
                r.sweep.chosen.objective(),
                r.sweep.prior.objective()
            );
        }
        eprintln!("thresholds: wrote {}", path.display());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;
    use crate::conformance::testkit::Scratch;

    fn committed() -> Result<String, String> {
        let path = root()?.join(CSV);
        std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))
    }

    #[test]
    fn arguments_are_refused_and_a_difference_names_its_line() {
        assert!(run(&["--nope".to_string()]).is_err() && run(&["a".into(), "b".into()]).is_err());
        assert!(compare("a\nb\n", "a\nb\n").is_ok());
        let e = compare("a\nb\n", "a\nc\n").err().unwrap_or_default();
        assert!(e.contains("line 2: committed `b`, fresh `c`"), "{e}");
        assert!(compare("a\n", "a\nb\n").is_err());
    }

    #[test]
    fn the_committed_sweep_is_what_a_fresh_one_writes() -> Result<(), String> {
        compare(&committed()?, &render(&sweep(&corpus_dir()?, &Coeff::ALL)?))
    }

    #[test]
    fn two_sweeps_write_the_same_bytes() -> Result<(), String> {
        let run = || sweep(&corpus_dir()?, &[Coeff::B, Coeff::E]).map(|s| render(&s));
        let (first, second) = (run()?, run()?);
        assert_eq!(first, second);
        assert_eq!(first.lines().count(), 3);
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
            names.push(cells[0]);
        }
        assert_eq!(names, ["k", "a", "b", "c", "d", "e"]);
        Ok(())
    }

    #[test]
    fn a_check_refuses_a_drifted_file_and_writes_nothing_while_a_sweep_writes_it(
    ) -> Result<(), String> {
        let scratch = Scratch::new("finish");
        let path = scratch.0.join("sweeps/thresholds.csv");
        let read = || std::fs::read_to_string(&path).map_err(|e| e.to_string());
        // No file: a check is an error and creates nothing.
        assert!(finish(&path, "a\nb\n", true).is_err() && !path.exists());
        // A write creates the directory and the file; a check of the same text passes.
        finish(&path, "a\nb\n", false)?;
        assert_eq!(read()?, "a\nb\n");
        finish(&path, "a\nb\n", true)?;
        // Drift, including a missing trailing newline, is an error naming the line, and the
        // committed file is left as it was.
        let e = finish(&path, "a\nc\n", true).err().unwrap_or_default();
        assert!(e.contains("line 2: committed `b`, fresh `c`"), "{e}");
        assert!(finish(&path, "a\nb", true).is_err());
        assert_eq!(read()?, "a\nb\n");
        // A write replaces it.
        finish(&path, "a\nc\n", false)?;
        assert_eq!(read()?, "a\nc\n");
        Ok(())
    }

    #[test]
    fn the_planted_c_of_section_10_is_dominated() -> Result<(), String> {
        // `c` with switch 1e-8 (the first grid point) and two terms: the exact arm's `12u/θ²`
        // above it, so the sweep ranks it far below the chosen candidate.
        let dir = corpus_dir()?;
        let series = Series::<D1>::load(&dir)?;
        let (grid, prior) = (grid::grid(), d12(&series)?);
        let samples = measure::samples(&dir, &series, Coeff::C)?.samples;
        let best = search::search(&samples, &grid, (prior.terms, prior.switch_z.v))?;
        let planted = search::score(&samples, 2, grid[0]).objective();
        assert!(planted > 1e6 * best.chosen.objective(), "{planted:e}");
        Ok(())
    }
}
