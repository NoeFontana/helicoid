//! Reading criterion's raw samples (`docs/PHASE1.md` §9).
//!
//! `target/criterion/<group>/<bench>/<slot>/sample.json` holds the run's `iters` and `times`, one
//! entry per sample: `times[i]` nanoseconds for `iters[i]` iterations. The per-iteration time is
//! their quotient — criterion's own point estimate is not read, because the gate needs the
//! distribution and not a summary.
//!
//! The estimates file is deliberately ignored: it carries criterion's own CI, computed by its own
//! unpaired bootstrap over one run, which answers a different question from §9's.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::Deserialize;
use serde::Serialize;

/// One benchmark's samples, in criterion's order.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub(super) struct Samples {
    /// `times[i] / iters[i]`, nanoseconds per iteration.
    pub(super) per_iter: Vec<f64>,
}

#[derive(Deserialize)]
struct RawSample {
    iters: Vec<f64>,
    times: Vec<f64>,
}

/// Every benchmark under `root`, keyed by `<group>/<bench>` with `/` restored from criterion's
/// flattened directory names, reading slot `slot` (`new` or `base`).
pub(super) fn read_all(root: &Path, slot: &str) -> Result<BTreeMap<String, Samples>, String> {
    let mut out = BTreeMap::new();
    let groups = std::fs::read_dir(root).map_err(|e| format!("{}: {e}", root.display()))?;
    for group in groups {
        let group = group.map_err(|e| e.to_string())?.path();
        if !group.is_dir() {
            continue;
        }
        let benches = std::fs::read_dir(&group).map_err(|e| format!("{}: {e}", group.display()))?;
        for bench in benches {
            let bench = bench.map_err(|e| e.to_string())?.path();
            let file = bench.join(slot).join("sample.json");
            if !file.is_file() {
                continue;
            }
            let name = key(&group, &bench)?;
            out.insert(name, read_one(&file)?);
        }
    }
    if out.is_empty() {
        return Err(format!(
            "no `{slot}/sample.json` under {}; run the benches first",
            root.display()
        ));
    }
    Ok(out)
}

/// `<group>/<bench>` from the two directory names. Criterion writes the names with `/` replaced by
/// `_`, so the key is the directory text and not the `bench_function` id: the two agree as long as
/// no id contains `_`, and the key is only ever compared with another run's key, never parsed.
fn key(group: &Path, bench: &Path) -> Result<String, String> {
    let part = |p: &Path| {
        p.file_name()
            .and_then(|s| s.to_str())
            .map(str::to_string)
            .ok_or_else(|| format!("{}: not a readable directory name", p.display()))
    };
    Ok(format!("{}/{}", part(group)?, part(bench)?))
}

/// One criterion `sample.json`, as `read_all` reads it.
///
/// Exposed so a recording can be replayed through exactly this parser: a recording is a byte copy
/// of criterion's own file, so the division and the validation below are replayed too, and not
/// merely trusted (`0035`, draft).
pub(super) fn read_file(file: &Path) -> Result<Samples, String> {
    read_one(&file.to_path_buf())
}

fn read_one(file: &PathBuf) -> Result<Samples, String> {
    let text = std::fs::read_to_string(file).map_err(|e| format!("{}: {e}", file.display()))?;
    let raw: RawSample =
        serde_json::from_str(&text).map_err(|e| format!("{}: {e}", file.display()))?;
    if raw.iters.len() != raw.times.len() || raw.iters.is_empty() {
        return Err(format!("{}: iters and times disagree", file.display()));
    }
    let mut per_iter = Vec::with_capacity(raw.iters.len());
    for (i, (&n, &t)) in raw.iters.iter().zip(&raw.times).enumerate() {
        // A NaN or non-positive iteration count would divide; spelled out rather than as a
        // negated comparison, which clippy reads as a mistake and which hides the NaN case.
        if n.is_nan() || n <= 0.0 || !t.is_finite() || t < 0.0 {
            return Err(format!(
                "{}: sample {i} is not a measurement",
                file.display()
            ));
        }
        per_iter.push(t / n);
    }
    Ok(Samples { per_iter })
}

/// Write the candidate's samples as the committed baseline, one file per benchmark.
///
/// Per-iteration times, not criterion's raw `iters`/`times`: the gate only ever reads the quotient,
/// and a file of quotients cannot be mistaken for a criterion slot. The directory is emptied first,
/// so a benchmark that no longer exists does not linger as a baseline nothing answers.
pub(super) fn write_all(dir: &Path, all: &BTreeMap<String, Samples>) -> Result<(), String> {
    if dir.exists() {
        std::fs::remove_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    }
    std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    for (name, s) in all {
        let file = dir.join(format!("{}.json", name.replace('/', "__")));
        let body = serde_json::to_string(s).map_err(|e| format!("{name}: {e}"))?;
        std::fs::write(&file, body).map_err(|e| format!("{}: {e}", file.display()))?;
    }
    Ok(())
}

/// Read the committed baseline `write_all` wrote.
pub(super) fn read_committed(dir: &Path) -> Result<BTreeMap<String, Samples>, String> {
    let entries = std::fs::read_dir(dir).map_err(|_| {
        format!(
            "no committed baseline at {}; `cargo xtask bench-gate --bless` writes one, and it is \
             this host's: timings do not travel",
            dir.display()
        )
    })?;
    let mut out = BTreeMap::new();
    for e in entries {
        let path = e.map_err(|e| e.to_string())?.path();
        let Some(stem) = path.file_stem().and_then(|s| s.to_str()) else {
            continue;
        };
        if path.extension().and_then(|s| s.to_str()) != Some("json") {
            continue;
        }
        let text =
            std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        let s: Samples =
            serde_json::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))?;
        out.insert(stem.replace("__", "/"), s);
    }
    if out.is_empty() {
        return Err(format!("{}: no baseline files", dir.display()));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn scratch(name: &str) -> Result<PathBuf, String> {
        let dir = std::env::temp_dir().join(format!("helicoid-bench-{name}"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        Ok(dir)
    }

    fn write(
        root: &Path,
        group: &str,
        bench: &str,
        slot: &str,
        iters: &[f64],
        times: &[f64],
    ) -> Result<(), String> {
        let dir = root.join(group).join(bench).join(slot);
        fs::create_dir_all(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        let body =
            format!("{{\"sampling_mode\":\"Linear\",\"iters\":{iters:?},\"times\":{times:?}}}");
        fs::write(dir.join("sample.json"), body).map_err(|e| format!("{}: {e}", dir.display()))
    }

    #[test]
    fn a_run_is_read_as_per_iteration_times() -> Result<(), String> {
        let root = scratch("read")?;
        write(&root, "g", "b", "new", &[10.0, 20.0], &[100.0, 220.0])?;
        let all = read_all(&root, "new")?;
        assert_eq!(all.len(), 1);
        let got = &all["g/b"].per_iter;
        assert_eq!(got.len(), 2);
        assert_eq!(got[0].to_bits(), 10.0_f64.to_bits());
        assert_eq!(got[1].to_bits(), 11.0_f64.to_bits());
        Ok(())
    }

    /// A slot that does not exist is not an empty run: the benches were not executed.
    #[test]
    fn a_missing_slot_is_an_error_not_an_empty_map() -> Result<(), String> {
        let root = scratch("missing")?;
        write(&root, "g", "b", "new", &[1.0], &[1.0])?;
        assert!(read_all(&root, "base").is_err());
        Ok(())
    }

    /// Nonsense in a sample file is an error, never a timing: a zero iteration count would divide.
    #[test]
    fn a_malformed_sample_is_an_error() -> Result<(), String> {
        let root = scratch("bad")?;
        write(&root, "g", "b", "new", &[0.0], &[1.0])?;
        assert!(read_all(&root, "new").is_err());
        let root = scratch("bad2")?;
        write(&root, "g", "b", "new", &[1.0, 2.0], &[1.0])?;
        assert!(read_all(&root, "new").is_err());
        Ok(())
    }

    /// A blessed baseline round-trips, and the keys survive the `/`-to-`__` filename escape.
    #[test]
    fn a_blessed_baseline_round_trips() -> Result<(), String> {
        let root = scratch("bless")?;
        let mut all = BTreeMap::new();
        all.insert(
            "coeffs_f64/exp_near".to_string(),
            Samples {
                per_iter: vec![1.5, 2.5],
            },
        );
        all.insert(
            "coeffs_f32/jr_near".to_string(),
            Samples {
                per_iter: vec![3.0],
            },
        );
        write_all(&root, &all)?;
        assert_eq!(read_committed(&root)?, all);
        // Blessing again drops a benchmark that no longer exists.
        let mut fewer = BTreeMap::new();
        fewer.insert(
            "coeffs_f64/exp_near".to_string(),
            Samples {
                per_iter: vec![1.5, 2.5],
            },
        );
        write_all(&root, &fewer)?;
        assert_eq!(read_committed(&root)?, fewer);
        Ok(())
    }

    /// No baseline is an error that says how to make one, never an empty comparison.
    #[test]
    fn a_missing_baseline_is_an_error() -> Result<(), String> {
        let root = scratch("nobase")?;
        let e = read_committed(&root.join("absent"))
            .err()
            .unwrap_or_default();
        assert!(e.contains("--bless"), "{e}");
        fs::create_dir_all(root.join("empty")).map_err(|e| e.to_string())?;
        assert!(read_committed(&root.join("empty")).is_err());
        Ok(())
    }

    /// Every group and bench is found, and the key joins the two directory names.
    #[test]
    fn every_group_and_bench_is_keyed_by_its_two_directories() -> Result<(), String> {
        let root = scratch("keys")?;
        for (g, b) in [
            ("coeffs_f64", "exp_near"),
            ("coeffs_f64", "jr_near"),
            ("coeffs_f32", "exp_near"),
        ] {
            write(&root, g, b, "new", &[2.0], &[4.0])?;
        }
        let all = read_all(&root, "new")?;
        assert_eq!(
            all.keys().map(String::as_str).collect::<Vec<_>>(),
            [
                "coeffs_f32/exp_near",
                "coeffs_f64/exp_near",
                "coeffs_f64/jr_near"
            ]
        );
        Ok(())
    }
}
