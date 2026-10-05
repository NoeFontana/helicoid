//! `cargo xtask bench-gate [--aa] [--bless] [--dry-run]` (`docs/PHASE1.md` §9).
//!
//! Runs `crates/helicoid`'s criterion benches twice, interleaved by criterion's own slots (`base`
//! then `new`, one `cargo bench` each), pairs the two runs' samples and reports the **paired
//! bootstrap 95% CI of the ratio** per benchmark. A run fails when the whole CI lies above
//! `1 + δ`, with `δ` the host's A/A noise floor from `--aa`, recorded in `baseline/HOST.md`.
//!
//! **A pair is adjacent in time, not just in index.** Each benchmark is run twice back to back and
//! those two runs are its pair; the whole set is not run twice. Measured: one benchmark run four
//! times in a row spreads 0.6%, while the same benchmark's two measurements taken ~2 minutes apart
//! (one pass of 60, then another) spread up to 90% — thermal and frequency state, and whatever else
//! the machine was doing, drift over a pass and do not cancel. Index pairing cannot recover that;
//! temporal adjacency can, and it costs the same total time.
//!
//! **`δ` is measured beside the comparison, not stored.** §9 records the A/A floor in
//! `baseline/HOST.md`, and on this host that file cannot be the allowance: the same protocol, the
//! same binary against itself, gave a worst `δ` of 0.0161 in one session and **0.6864** in another
//! a few hours later — 23 of 60 benchmarks above 2% where the first run had none. A floor from the
//! quiet session fails identical code; one from the noisy session would pass a 50% regression. So
//! `--against` runs its own A/A control for every benchmark, in the same invocation, bracketing the
//! candidate (`b₁ c b₂`), and judges the candidate against that. `HOST.md` stays a log of what the
//! machine has shown, which is what tells a human whether it is worth running at all.
//!
//! Two allowances that look more frugal were measured and refuted first. **One floor per
//! benchmark**, stored: across two quiet sessions the distribution over the 60 reproduces (median
//! 0.0029 then 0.0030) while single benchmarks move tenfold (0.0019 → 0.0190), so the benchmarks
//! that measured quiet get an allowance nothing can meet and 8 of 60 failed identical code.
//! **An accumulated maximum** over A/A runs: refuted by the first run after it, which pinned the
//! allowance at 68% permanently. Sixty draws estimate an extreme; one does not; and no number
//! measured last Tuesday describes this minute.
//!
//! **A regression has to appear in every pair, and the pairs are spread over the whole run.** Each
//! benchmark yields `REPLICATES` × 2 candidate pairs — both brackets of each triplet — and the gate
//! fails only when every one puts its whole CI above `1 + floor`. §9's rule is the per-pair test;
//! requiring all of them is what a *short* disturbance cannot satisfy. It is not enough on its own:
//! with the replicate loop innermost, a benchmark's six pairs were taken inside about 30 s and so
//! were not independent, and 2 of 60 identical-code benchmarks failed with a quiet control. The
//! replicate loop is therefore the outer one, which puts a benchmark's three triplets a pass apart
//! — about 20 min here — at no extra cost, while keeping the three windows of a triplet contiguous.
//!
//! A middle-window bias was the first explanation and the data refused it: over 60 benchmarks the
//! median candidate-to-baseline ratio is 0.9994 and 0.9989 under two `criterion` versions, with 24
//! and 18 of 60 above 1.0. The candidate always runs second of three, so a position effect would
//! have shown as a systematic offset, and there is none.
//!
//! **`--against <binary>` gates; the committed-baseline comparison reports and never gates.** §9's
//! "interleaved baseline/candidate runs" needs the two *binaries* present at once, so that a
//! benchmark's baseline and candidate measurements are adjacent in time the way `--aa`'s pair is.
//! `--against` does exactly that and fails. Comparing a *committed* set of samples against a later
//! run cannot be adjacent: the baseline is from another point in time by construction, and measured
//! on this host that costs far more than any signal — identical code against samples blessed
//! earlier the same day read **0.49× to 2.04×**, 27 of 60 benchmarks past their own floor, against
//! an adjacent floor of at most **1.90%**. Part of that is the host's own CPU allocation changing
//! between runs (`available_parallelism` reported 1 at one bless and 8 later), which no statistic
//! repairs. So the default mode prints ratios beside each benchmark's floor and fails nothing: a
//! gate whose floor does not cover its own comparison is an arbitrary percentage wearing a
//! measurement's clothes, which is what §9 refuses.
//!
//! **A failure says the two binaries differ there, not that the change caused it.** Measured by
//! negative control: slowing `coeffs::kernel::exact_a` failed the five `jr_coeffs` benchmarks that
//! evaluate it, at 1.38x to 1.57x, *and* four `exp_coeffs` benchmarks at 1.016x to 1.041x, which
//! have no path to `exact_a` at all — one caller. Between the two binaries only 2 of 3192 text
//! symbols changed size, `.text` grew 64 bytes, and 2997 symbols kept an identical size at a
//! different address: byte-identical code, relocated, and alignment is worth percent. Source
//! locality does not survive to machine level, so a failing row that cannot route through the
//! change is a layout finding, confirmed with a symbol diff (`nm --defined-only -S`, size against
//! address) and not with another benchmark run.
//!
//! Oracle rows are reported, never gated (§9), and none exist: no oracle answers a latency question.

mod bootstrap;
mod samples;

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use bootstrap::{ratio_ci, Interval, RESAMPLES};

/// The bootstrap's seed. Fixed so a verdict is a function of the samples alone.
const SEED: u64 = 0x5EED_B3C0_DEAD_BEEF;

/// Where the measured floor lives.
const HOST: &str = "baseline/HOST.md";

/// The bench target `--bench` selects when it is not given. `coeffs` is the coefficient kernel's
/// 60 benchmarks; `groups` is `PHASE3.md` §11's, in a binary of its own so a change to one does not
/// pay for the other's windows.
const BENCH: &str = "coeffs";

/// The two named slots each benchmark's adjacent pair writes.
///
/// **Not `base`/`new`.** Criterion copies its own `new` into `base` *within* a run, so after any
/// number of runs the two are byte-identical and every ratio is exactly 1 — an A/A floor of
/// `0.0000` across every benchmark, which is what a gate reading those slots would report, and
/// which would pass anything. `--save-baseline NAME` writes a slot that persists, so two passes
/// leave two measurements.
const AA_SLOTS: (&str, &str) = ("aa-1", "aa-2");

/// The third slot `--against` needs: a candidate window is **bracketed** by two baseline windows.
///
/// `b₁ c b₂` in time order, so the control pair `(b₁, b₂)` spans *more* elapsed time than either
/// candidate pair `(b₁, c)` and `(b₂, c)`. A floor measured over a longer span than the comparison
/// it licenses can only be conservative, and the reverse order (`b₁ b₂ c`) would measure the floor
/// over a shorter span than the comparison, which is the asymmetry that licenses a false verdict.
const BRACKET_SLOT: &str = "aa-3";

/// Where `--bless` keeps the committed samples the gate compares against.
const BASELINE: &str = "baseline/bench";

/// How many adjacent pairs a benchmark measures. In `--aa` the **least disturbed** one carries the
/// verdict; in `--against` every triplet's pairs must agree, and the triplets are a pass apart.
///
/// One pair is not enough, and the reason is measured. The median removes contamination *within* a
/// window, but it cannot tell a uniformly slow window from slower code — the interquartile spread
/// is tight either way. Those whole-window shifts are rare and large: the same binary against
/// itself put 6 of 60 verdicts past their own floor, one at 1.41x, on a benchmark that measured
/// 10.81 ns in six consecutive runs before and after.
///
/// Pooling the three does not fix it either. A disturbance long enough to cover two consecutive
/// runs — and one benchmark's three pairs take about 15 s here — puts two of a side's three windows
/// over, which is past a median's 50% breakdown point: pooling gave one benchmark a floor of 58.84%
/// where ten consecutive windows of it spread 2.38%.
///
/// Contamination is one-sided at **both** levels: a disturbed sample is slow and never fast, and so
/// is a disturbed window. So the window level takes the same treatment as the sample level — the
/// pair with the lowest total cost is the one least disturbed, and it alone is reported. Three
/// survives two bad pairs. The statistic stays §9's: one paired bootstrap CI of a ratio over one
/// pair of adjacent windows.
const REPLICATES: usize = 3;

struct Options {
    aa: bool,
    bless: bool,
    dry_run: bool,
    /// A prebuilt baseline bench binary to alternate with this tree's, per benchmark.
    against: Option<PathBuf>,
    /// Where to persist every window `--against` measures, so the rule can be replayed from it.
    record: Option<PathBuf>,
    /// A recording to run the rule over, executing no benchmark.
    replay: Option<PathBuf>,
    /// Which `[[bench]]` target to build and run ([`BENCH`] by default).
    bench: String,
    /// Run only the benchmarks whose id contains one of these, instead of every id the target
    /// lists. Repeat the flag to keep several.
    ///
    /// The rule is per benchmark, so a subset is the same measurement over fewer rows -- what it
    /// cannot do is report the suite's worst, so a *gate* run claims nothing here. It is for
    /// deciding one change against the rows that change can move, with the A/A control measured
    /// over the same subset.
    only: Vec<String>,
}

const USAGE: &str = "usage: cargo xtask bench-gate [--aa] [--bless] [--dry-run] \
                     [--bench <target>] [--only <substring>] \
                     [--against <bench-binary> [--record <dir>]] [--replay <dir>]";

fn parse_args(args: &[String]) -> Result<Options, String> {
    let mut o = Options {
        aa: false,
        bless: false,
        dry_run: false,
        against: None,
        record: None,
        replay: None,
        bench: BENCH.to_string(),
        only: Vec::new(),
    };
    let mut rest = args.iter();
    while let Some(a) = rest.next() {
        match a.as_str() {
            "--aa" => o.aa = true,
            "--bless" => o.bless = true,
            "--dry-run" => o.dry_run = true,
            "--against" => {
                let path = rest
                    .next()
                    .ok_or_else(|| format!("`--against` takes a path; {USAGE}"))?;
                o.against = Some(PathBuf::from(path));
            }
            "--record" => {
                let path = rest
                    .next()
                    .ok_or_else(|| format!("`--record` takes a directory; {USAGE}"))?;
                o.record = Some(PathBuf::from(path));
            }
            "--replay" => {
                let path = rest
                    .next()
                    .ok_or_else(|| format!("`--replay` takes a directory; {USAGE}"))?;
                o.replay = Some(PathBuf::from(path));
            }
            "--bench" => {
                let name = rest
                    .next()
                    .ok_or_else(|| format!("`--bench` takes a target name; {USAGE}"))?;
                o.bench = name.clone();
            }
            "--only" => {
                let pat = rest
                    .next()
                    .ok_or_else(|| format!("`--only` takes a substring; {USAGE}"))?;
                o.only.push(pat.clone());
            }
            _ => return Err(format!("unknown argument `{a}`; {USAGE}")),
        }
    }
    if o.aa && o.against.is_some() {
        return Err(
            "`--aa` measures one binary against itself; `--against` is the other mode".to_string(),
        );
    }
    if o.replay.is_some() && (o.aa || o.against.is_some() || o.bless) {
        return Err(
            "`--replay` runs the rule over a recording and measures nothing, so it takes no other \
             mode"
                .to_string(),
        );
    }
    if o.record.is_some() && o.against.is_none() {
        return Err(
            "`--record` persists what `--against` measures; pass `--against` too".to_string(),
        );
    }
    Ok(o)
}

pub(crate) fn run(args: &[String]) -> Result<(), String> {
    let o = parse_args(args)?;
    let root = crate::conformance::root()?;
    let criterion = root.join("target").join("criterion");
    // Before `list` and `build_bench`: a replay runs no benchmark, so it must not need one built.
    if let Some(dir) = &o.replay {
        let mode = format!("--replay {}", dir.display());
        return verdict(&read_recording(dir)?, Some(&root), &mode);
    }
    let all = list(&root, &o.bench)?;
    let ids: Vec<String> = match o.only.is_empty() {
        true => all,
        false => {
            let kept: Vec<String> = all
                .into_iter()
                .filter(|id| o.only.iter().any(|p| id.contains(p)))
                .collect();
            if kept.is_empty() {
                return Err(format!(
                    "`--only` matched no benchmark of `{}`: {:?}",
                    o.bench, o.only
                ));
            }
            kept
        }
    };
    let mine = build_bench(&root, &o.bench)?;
    eprintln!("bench-gate: {} benchmarks", ids.len());

    if o.aa {
        let mut rows = Vec::with_capacity(ids.len());
        for (i, id) in ids.iter().enumerate() {
            eprintln!("bench-gate --aa: {}/{} {id}", i + 1, ids.len());
            rows.push(pooled_pair(&root, &criterion, id, &mine, &mine)?);
        }
        let report = Report { rows };
        eprintln!("{}", report.text(None));
        let worst = report.worst_delta();
        eprintln!(
            "bench-gate --aa: δ per benchmark, worst {:.4} ({:.2}%), median {:.4}; {RESAMPLES} resamples",
            worst,
            worst * 100.0,
            report.median_delta()
        );
        if worst <= 0.0 {
            return Err(
                "every δ is zero: the two runs measured the same bytes, so this is not an \
                        A/A floor"
                    .to_string(),
            );
        }
        if o.bless && !o.dry_run {
            let path = root.join(HOST);
            let previous = read_host(&path).unwrap_or_default();
            let runs = write_host(&path, &report, &previous)?;
            eprintln!(
                "bench-gate: wrote {} — floor {:.4} over {} recorded A/A run(s)",
                path.display(),
                last_max(&runs),
                runs.len()
            );
        } else if o.bless {
            eprintln!("bench-gate --dry-run: would write {HOST}");
        }
        return Ok(());
    }

    if let Some(baseline) = o.against {
        return against(
            &root,
            &criterion,
            &ids,
            &baseline,
            &mine,
            o.record.as_deref(),
        );
    }

    bench(&root, &mine, AA_SLOTS.1, None)?;
    let candidate = samples::read_all(&criterion, AA_SLOTS.1)?;

    if o.bless {
        let dir = root.join(BASELINE);
        if o.dry_run {
            eprintln!(
                "bench-gate --dry-run: would write {} benchmark(s) under {BASELINE}",
                candidate.len()
            );
            return Ok(());
        }
        samples::write_all(&dir, &candidate)?;
        eprintln!(
            "bench-gate --bless: wrote {} benchmark(s) under {BASELINE} ({})",
            candidate.len(),
            host_line()
        );
        return Ok(());
    }

    let floor = last_max(&read_host(&root.join(HOST))?);
    let baseline = samples::read_committed(&root.join(BASELINE))?;
    let report = compare(&baseline, &candidate)?;
    eprintln!("{}", report.text(Some(floor)));
    // Reported, not gated: see the module doc. The counts are printed so the drift is visible and
    // so the owed `--against` mode has a number to beat.
    let over = report.rows.iter().filter(|r| r.fails(floor)).count();
    let worst = report
        .rows
        .iter()
        .map(|r| r.ci.point)
        .fold(0.0_f64, f64::max);
    eprintln!(
        "bench-gate: {} benchmarks against the blessed samples; {over} would fail 1 + δ \
         (δ {floor:.4}), worst ratio {worst:.4}.",
        report.rows.len()
    );
    eprintln!(
        "bench-gate: NOT a verdict. A committed baseline is not adjacent in time to this run, and \
         on this host that drift is about 60x the adjacent floor (`baseline/HOST.md`). §9's gate \
         needs `--against <ref>`, which is owed."
    );
    Ok(())
}

/// §9's gate: the candidate **bracketed** by the baseline, with the floor measured in the same run.
///
/// This is the mode a stored floor cannot be. Each benchmark runs `b1 c b2` — baseline, candidate,
/// baseline — `REPLICATES` times. `(b1, b2)` is an A/A control measured *in this invocation, on
/// this benchmark, over a longer span than the comparison*, and it is the floor the candidate's two
/// pairs are judged against. A benchmark fails only when **every** pair, over every replicate, puts
/// its whole CI above `1 + floor`.
///
/// A floor read from `HOST.md` is not valid here, and the measurement says so: the same `--aa`
/// protocol on this host gave a worst `δ` of 0.0161 in one session and **0.6864** in another a few
/// hours later, 18 of 60 benchmarks above 5% where the first had none above 2%. Nothing about the
/// protocol changed; the machine's state did, and a shared VM does not announce it. The quiet
/// session's floor fails the identical binary against itself; the noisy one's would pass a 50%
/// regression. Only a control measured beside the comparison cancels that, and it costs one extra
/// window per replicate.
///
/// The baseline is a prebuilt criterion bench binary: in another tree,
/// `cargo bench -p helicoid --bench coeffs --features __sweep --no-run` prints the path to pass.
/// `--against <ref>`, which would build it from a git worktree, is owed (`0033`, draft).
fn against(
    root: &Path,
    criterion: &Path,
    ids: &[String],
    baseline: &Path,
    mine: &Path,
    recording: Option<&Path>,
) -> Result<(), String> {
    if !baseline.is_file() {
        return Err(format!("{}: not a bench binary", baseline.display()));
    }
    // **The replicate loop is the outer one.** Measured: with it inside, a benchmark's six pairs are
    // all taken inside about 30 s, so they are not independent and a disturbance outlasting that span
    // defeats the all-pairs rule — 2 of 60 identical-code benchmarks failed that way, both with a
    // quiet control. Outside, a benchmark's three triplets are a whole pass apart, about 20 min here,
    // so such a disturbance can spoil one triplet of three and never all six pairs. The windows
    // within a triplet stay contiguous, which is the adjacency the pairing needs.
    let mut windows = Vec::with_capacity(REPLICATES * ids.len());
    for rep in 0..REPLICATES {
        for (i, id) in ids.iter().enumerate() {
            eprintln!(
                "bench-gate --against: replicate {}/{REPLICATES}, {}/{} {id}",
                rep + 1,
                i + 1,
                ids.len()
            );
            let w = measure(root, criterion, id, baseline, mine)?;
            if let Some(dir) = recording {
                record(dir, criterion, rep + 1, &w.name)?;
            }
            windows.push(w);
        }
    }
    if let Some(dir) = recording {
        eprintln!(
            "bench-gate --against: recorded {} triplets under {}; `--replay` re-runs the rule over \
             them and measures nothing",
            windows.len(),
            dir.display()
        );
    }
    verdict(&windows, Some(root), "--against")
}

/// One benchmark's `--against` verdict: the floor measured beside it, and every candidate pair.
#[derive(Debug, Clone)]
struct Verdict {
    name: String,
    /// The noisiest control pair of the `REPLICATES` replicates, as a swap-invariant `δ`.
    ///
    /// The **maximum**, not the median: the floor has to cover the span the candidate's own pairs
    /// were measured over, and a replicate that was disturbed at all was disturbed for the
    /// candidate too.
    floor: f64,
    /// `REPLICATES` × 2 CIs of candidate-over-baseline: each replicate's two brackets.
    pairs: Vec<Interval>,
}

impl Verdict {
    /// §9's rule per pair, and the regression must hold in **every** pair.
    ///
    /// One pair is §9's own test. Requiring all of them is what stops an isolated disturbance from
    /// producing a verdict: a real slowdown is in every window, and the measurement of what is not
    /// is in this module's doc.
    fn fails(&self) -> bool {
        !self.pairs.is_empty() && self.pairs.iter().all(|c| c.low > 1.0 + self.floor)
    }

    /// The least any pair claims: the slowdown every pair agrees on at least.
    fn claim(&self) -> f64 {
        self.pairs
            .iter()
            .map(|c| c.low)
            .fold(f64::INFINITY, f64::min)
    }

    /// The middle pair's ratio, for ordering the report.
    fn point(&self) -> f64 {
        let points: Vec<f64> = self.pairs.iter().map(|c| c.point).collect();
        bootstrap::median(&points)
    }
}

fn table(verdicts: &[Verdict]) -> String {
    let mut out =
        String::from("  benchmark                                  point  least  floor\n");
    let mut rows: Vec<&Verdict> = verdicts.iter().collect();
    rows.sort_by(|a, b| b.point().total_cmp(&a.point()));
    for v in rows {
        out.push_str(&format!(
            "  {:<40} {:6.4} {:6.4} {:6.4}{}\n",
            v.name,
            v.point(),
            v.claim(),
            v.floor,
            if v.fails() { "  FAIL" } else { "" }
        ));
    }
    out
}

/// One `b1 c b2` triplet: the control's `δ` and the candidate's two pairs, adjacent in time.
///
/// One triplet, not `REPLICATES` of them, because the replicate loop is the **outer** one — see
/// `against`. Within a triplet the three windows are contiguous, which is what makes the pairing
/// cancel drift at all.
struct Triplet {
    name: String,
    floor: f64,
    pairs: [Interval; 2],
}

/// One triplet's three windows, before any statistic touches them.
///
/// **This is the seam.** Measuring is slow and needs a quiet machine; the rule that turns these
/// three vectors into a verdict is a pure, seeded function and needs neither. Keeping them apart is
/// what lets `--replay` answer in milliseconds what `--against` answers in twenty minutes
/// (`0035`, draft).
#[derive(Debug, Clone, PartialEq)]
struct Windows {
    name: String,
    /// The baseline window before the candidate.
    first: Vec<f64>,
    candidate: Vec<f64>,
    /// The baseline window after it, which is what makes the control span the comparison.
    second: Vec<f64>,
}

/// The rule, over one triplet. Pure: no process, no clock, no filesystem.
fn triplet_of(w: &Windows) -> Result<Triplet, String> {
    let control = ratio_ci(&w.first, &w.second, SEED)
        .ok_or_else(|| format!("{}: the control windows do not pair", w.name))?;
    let pair = |base: &[f64]| {
        ratio_ci(base, &w.candidate, SEED)
            .ok_or_else(|| format!("{}: the candidate windows do not pair", w.name))
    };
    Ok(Triplet {
        name: w.name.clone(),
        floor: delta(&control),
        pairs: [pair(&w.first)?, pair(&w.second)?],
    })
}

/// Measure one triplet: `b1 c b2`, contiguous.
fn measure(
    root: &Path,
    criterion: &Path,
    id: &str,
    baseline: &Path,
    mine: &Path,
) -> Result<Windows, String> {
    let mut name = String::new();
    // The slots persist, and `read_all` reads every benchmark that has one, so without this each
    // triplet would re-read the previous one's samples as extra rows.
    clear_slots(criterion)?;
    bench(root, baseline, AA_SLOTS.0, Some(id))?;
    bench(root, mine, AA_SLOTS.1, Some(id))?;
    bench(root, baseline, BRACKET_SLOT, Some(id))?;
    let first = one(criterion, AA_SLOTS.0, id, &mut name)?;
    let candidate = one(criterion, AA_SLOTS.1, id, &mut name)?;
    let second = one(criterion, BRACKET_SLOT, id, &mut name)?;
    Ok(Windows {
        name,
        first,
        candidate,
        second,
    })
}

/// Accumulate triplets into one verdict per benchmark: the floor is the noisiest control, and every
/// pair of every triplet must agree before the gate fails.
fn accumulate(windows: &[Windows]) -> Result<Vec<Verdict>, String> {
    let mut acc: BTreeMap<String, Verdict> = BTreeMap::new();
    for w in windows {
        let t = triplet_of(w)?;
        let v = acc.entry(t.name.clone()).or_insert_with(|| Verdict {
            name: t.name.clone(),
            floor: 0.0,
            pairs: Vec::with_capacity(2 * REPLICATES),
        });
        v.floor = v.floor.max(t.floor);
        v.pairs.extend_from_slice(&t.pairs);
    }
    Ok(acc.into_values().collect())
}

/// The report and the exit verdict, shared by `--against` and `--replay` so that a replay cannot
/// drift from the gate it stands in for.
fn verdict(windows: &[Windows], root: Option<&Path>, mode: &str) -> Result<(), String> {
    let verdicts = accumulate(windows)?;
    if verdicts.is_empty() {
        return Err(format!("bench-gate {mode}: no windows"));
    }
    eprintln!("{}", table(&verdicts));

    // What this run could resolve at all, which is a property of the machine it ran on and not of
    // the candidate. Printed whatever the verdict: a pass under a loose floor and a pass under a
    // tight one are not the same statement.
    let floors: Vec<f64> = verdicts.iter().map(|v| v.floor).collect();
    let worst = floors.iter().copied().fold(0.0_f64, f64::max);
    let median = bootstrap::median(&floors);
    eprintln!(
        "bench-gate {mode}: concurrent A/A floor, median {:.4} ({:.2}%), worst {:.4} ({:.2}%) over \
         {} benchmarks. Nothing smaller than a benchmark's own floor is resolvable by this run.",
        median,
        median * 100.0,
        worst,
        worst * 100.0,
        verdicts.len()
    );
    if let Some(root) = root {
        if let Ok(recorded) = read_host(&root.join(HOST)) {
            eprintln!(
                "bench-gate {mode}: {HOST} records {:.4} over {} A/A run(s); it is not the \
                 allowance — this run's own control is.",
                last_max(&recorded),
                recorded.len()
            );
        }
    }

    let failed: Vec<&Verdict> = verdicts.iter().filter(|v| v.fails()).collect();
    if failed.is_empty() {
        eprintln!(
            "bench-gate {mode}: {} benchmarks, none above its own concurrent floor in all {} pairs.",
            verdicts.len(),
            2 * REPLICATES
        );
        return Ok(());
    }
    Err(format!(
        "bench-gate {mode}: {} of {} benchmark(s) regressed past the floor measured beside \
         them:\n{}",
        failed.len(),
        verdicts.len(),
        failed
            .iter()
            .map(|v| format!(
                "  {} at least {:.4} in every one of {} pairs (floor {:.4})",
                v.name,
                v.claim(),
                v.pairs.len(),
                v.floor
            ))
            .collect::<Vec<_>>()
            .join("\n")
    ))
}

/// The one benchmark a `--exact` filtered run wrote to `slot`, and its name.
fn one(criterion: &Path, slot: &str, id: &str, name: &mut String) -> Result<Vec<f64>, String> {
    let all = samples::read_all(criterion, slot)?;
    if all.len() != 1 {
        return Err(format!(
            "{id}: the filter matched {} benchmarks in `{slot}`, not 1",
            all.len()
        ));
    }
    let mut out = Vec::new();
    for (key, s) in &all {
        *name = key.clone();
        out.extend_from_slice(&s.per_iter);
    }
    Ok(out)
}

/// Criterion's filename for a run's raw samples.
const SAMPLE: &str = "sample.json";

/// The slots a triplet writes, in the order they are measured.
const TRIPLET_SLOTS: [&str; 3] = [AA_SLOTS.0, AA_SLOTS.1, BRACKET_SLOT];

/// Persist one triplet, so a slow run stops deleting its evidence.
///
/// **The recording is criterion's own `sample.json`, copied byte for byte.** That is the right
/// format for this data and the reasoning is fidelity, not size: a copy is lossless by
/// construction, it needs no second parser, it can be diffed against a live `target/criterion`
/// tree, and replaying it goes through `samples::read_file` — so the `times`/`iters` division and
/// its validation are replayed rather than assumed. Storing the derived per-iteration quotients
/// instead, which was the first design, would bake today's derivation into the recording and make
/// exactly that step unreplayable. `iters` is redundant under criterion's `Linear` mode, an
/// arithmetic sequence expressible in three numbers, and is kept anyway: dropping it would end the
/// byte-copy property for a few kilobytes (`0035`, draft).
fn record(dir: &Path, criterion: &Path, rep: usize, name: &str) -> Result<(), String> {
    for slot in TRIPLET_SLOTS {
        let from = criterion.join(name).join(slot).join(SAMPLE);
        let into = dir.join(name).join(format!("rep{rep}")).join(slot);
        std::fs::create_dir_all(&into).map_err(|e| format!("{}: {e}", into.display()))?;
        let to = into.join(SAMPLE);
        std::fs::copy(&from, &to)
            .map_err(|e| format!("{} -> {}: {e}", from.display(), to.display()))?;
    }
    Ok(())
}

/// Read a recording: `<group>/<bench>/rep<N>/<slot>/sample.json`, in path order, so that a replay
/// is a function of the directory alone.
fn read_recording(dir: &Path) -> Result<Vec<Windows>, String> {
    let mut found: BTreeMap<(String, String), PathBuf> = BTreeMap::new();
    collect(dir, dir, &mut found)?;
    let mut out = Vec::with_capacity(found.len());
    for ((name, _rep), at) in &found {
        let window = |slot: &str| -> Result<Vec<f64>, String> {
            Ok(samples::read_file(&at.join(slot).join(SAMPLE))?.per_iter)
        };
        out.push(Windows {
            name: name.clone(),
            first: window(TRIPLET_SLOTS[0])?,
            candidate: window(TRIPLET_SLOTS[1])?,
            second: window(TRIPLET_SLOTS[2])?,
        });
    }
    if out.is_empty() {
        return Err(format!(
            "{}: no recorded triplets; `--against <binary> --record <dir>` writes them",
            dir.display()
        ));
    }
    Ok(out)
}

/// Every `rep<N>` directory under `root`, keyed by the `<group>/<bench>` path above it and by the
/// replicate, so the map's order is the path order.
fn collect(
    root: &Path,
    at: &Path,
    out: &mut BTreeMap<(String, String), PathBuf>,
) -> Result<(), String> {
    for e in std::fs::read_dir(at).map_err(|e| format!("{}: {e}", at.display()))? {
        let path = e.map_err(|e| e.to_string())?.path();
        if !path.is_dir() {
            continue;
        }
        let leaf = path
            .file_name()
            .and_then(|s| s.to_str())
            .ok_or_else(|| format!("{}: not a readable directory name", path.display()))?;
        if let Some(rep) = leaf.strip_prefix("rep") {
            let relative = path
                .parent()
                .and_then(|p| p.strip_prefix(root).ok())
                .ok_or_else(|| format!("{}: not under the recording root", path.display()))?;
            let name = relative
                .to_str()
                .ok_or_else(|| format!("{}: not a readable path", relative.display()))?
                .replace(std::path::MAIN_SEPARATOR, "/");
            out.insert((name, rep.to_string()), path.clone());
        } else {
            collect(root, &path, out)?;
        }
    }
    Ok(())
}

/// One benchmark's verdict: `REPLICATES` alternations, the least disturbed pair reported.
///
/// The two binaries are the same file for `--aa`. They alternate rather than running in blocks, so
/// each pair stays adjacent in time, and the pair selected is the one whose two windows cost least
/// in total — see `REPLICATES`.
fn pooled_pair(
    root: &Path,
    criterion: &Path,
    id: &str,
    first_binary: &Path,
    second_binary: &Path,
) -> Result<Row, String> {
    let mut name = String::new();
    let mut best: Option<(f64, Vec<f64>, Vec<f64>)> = None;
    for _ in 0..REPLICATES {
        // The slots persist, and `read_all` reads every benchmark that has one, so without this
        // each replicate would re-read the previous one's samples as extra rows.
        clear_slots(criterion)?;
        bench(root, first_binary, AA_SLOTS.0, Some(id))?;
        bench(root, second_binary, AA_SLOTS.1, Some(id))?;
        let first = samples::read_all(criterion, AA_SLOTS.0)?;
        let second = samples::read_all(criterion, AA_SLOTS.1)?;
        if first.len() != 1 || second.len() != 1 {
            return Err(format!(
                "{id}: the filter matched {} and {} benchmarks, not 1",
                first.len(),
                second.len()
            ));
        }
        let mut a = Vec::new();
        let mut b = Vec::new();
        for (key, s) in &first {
            name = key.clone();
            a.extend_from_slice(&s.per_iter);
        }
        for s in second.values() {
            b.extend_from_slice(&s.per_iter);
        }
        let cost = bootstrap::median(&a) + bootstrap::median(&b);
        if best
            .as_ref()
            .is_none_or(|&(best_cost, _, _)| cost < best_cost)
        {
            best = Some((cost, a, b));
        }
    }
    let (_, a, b) = best.ok_or_else(|| format!("{id}: no pair was measured"))?;
    let ci = ratio_ci(&a, &b, SEED)
        .ok_or_else(|| format!("{id}: {} and {} samples do not pair", a.len(), b.len()))?;
    Ok(Row { name, ci })
}

/// Remove both A/A slots wherever criterion wrote them, so the next pair's read sees only what
/// that pair measured. Criterion's other slots (`new`, a blessed baseline) are left alone.
fn clear_slots(criterion: &Path) -> Result<(), String> {
    if !criterion.is_dir() {
        return Ok(());
    }
    for group in std::fs::read_dir(criterion).map_err(|e| e.to_string())? {
        let group = group.map_err(|e| e.to_string())?.path();
        if !group.is_dir() {
            continue;
        }
        for bench in std::fs::read_dir(&group).map_err(|e| e.to_string())? {
            let bench = bench.map_err(|e| e.to_string())?.path();
            for slot in [AA_SLOTS.0, AA_SLOTS.1, BRACKET_SLOT] {
                let dir = bench.join(slot);
                if dir.is_dir() {
                    std::fs::remove_dir_all(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
                }
            }
        }
    }
    Ok(())
}

/// The benchmark ids, from the harness itself: `--list` prints `<id>: benchmark`.
fn list(root: &Path, bench: &str) -> Result<Vec<String>, String> {
    let out = Command::new("cargo")
        .current_dir(root)
        .args([
            "bench",
            "-p",
            "helicoid",
            "--bench",
            bench,
            "--features",
            "__sweep",
            "--",
            "--list",
        ])
        .output()
        .map_err(|e| format!("cargo bench --list: {e}"))?;
    if !out.status.success() {
        return Err(format!("cargo bench --list: {}", out.status));
    }
    let text = String::from_utf8_lossy(&out.stdout);
    let ids: Vec<String> = text
        .lines()
        .filter_map(|l| l.strip_suffix(": benchmark"))
        .map(str::to_string)
        .collect();
    if ids.is_empty() {
        return Err("cargo bench --list named no benchmarks".to_string());
    }
    Ok(ids)
}

/// Build this tree's bench binary and return its path, with `__sweep` on the command line and
/// never in a manifest (`cargo xtask lint` check 5: only `xtask` may request it).
///
/// Built once and then executed directly, so that **every** measurement in every mode goes through
/// one invocation path. Running one side under `cargo` and the other direct is not symmetric —
/// cargo spends its own CPU immediately before the run it wraps — and a floor measured with cargo
/// on both sides does not cover a comparison with cargo on one. Measured: that asymmetry alone put
/// 3 of 60 identical-code verdicts past their floor, at ratios of 1.007 to 1.018.
fn build_bench(root: &Path, bench: &str) -> Result<PathBuf, String> {
    let out = Command::new("cargo")
        .current_dir(root)
        .args([
            "bench",
            "-p",
            "helicoid",
            "--bench",
            bench,
            "--features",
            "__sweep",
            "--no-run",
            "--message-format=json",
        ])
        .output()
        .map_err(|e| format!("cargo bench --no-run: {e}"))?;
    if !out.status.success() {
        return Err(format!(
            "cargo bench --no-run: {}\n{}",
            out.status,
            String::from_utf8_lossy(&out.stderr)
        ));
    }
    let text = String::from_utf8_lossy(&out.stdout);
    let mut found = None;
    for line in text.lines() {
        let Ok(v) = serde_json::from_str::<serde_json::Value>(line) else {
            continue;
        };
        if v.get("reason").and_then(|r| r.as_str()) != Some("compiler-artifact") {
            continue;
        }
        if v.pointer("/target/name").and_then(|n| n.as_str()) != Some(bench) {
            continue;
        }
        if let Some(exe) = v.get("executable").and_then(|e| e.as_str()) {
            found = Some(PathBuf::from(exe));
        }
    }
    found.ok_or_else(|| format!("cargo bench --no-run named no `{bench}` executable"))
}

/// One run of a criterion bench binary over one benchmark, into one slot.
///
/// `--bench` is what cargo passes to put criterion in bench mode; the cwd is the workspace root so
/// that every binary writes the same `target/criterion` tree.
fn bench(root: &Path, binary: &Path, slot: &str, only: Option<&str>) -> Result<(), String> {
    let mut cmd = Command::new(binary);
    cmd.current_dir(root)
        .args(["--bench", "--save-baseline", slot]);
    if let Some(id) = only {
        // criterion's filter is a regex, so an id's `.` would match any character; `--exact` is
        // what makes the id literal.
        cmd.args(["--exact", id]);
    }
    let status = cmd
        .status()
        .map_err(|e| format!("{}: {e}", binary.display()))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("{}: {status}", binary.display()))
    }
}

struct Row {
    name: String,
    ci: Interval,
}

impl Row {
    /// §9: fail when the **whole** CI lies above `1 + δ`.
    fn fails(&self, floor: f64) -> bool {
        self.ci.low > 1.0 + floor
    }

    /// The allowance this row's interval justifies, **invariant under swapping the two runs**.
    ///
    /// Which of a pair is "baseline" is a slot name, so a floor read off `high` alone charges the
    /// pair for drift in one direction and nothing for the same drift in the other: on this host
    /// four A/A pairs came out 22% to 43% *faster*, and `high - 1` gave each of them the tightest
    /// floor in the run. Relabelling replaces the interval by its reciprocal, so the largest
    /// slowdown the pair exhibits under either labelling is `max(high - 1, 1/low - 1)`.
    fn delta(&self) -> f64 {
        delta(&self.ci)
    }
}

struct Report {
    rows: Vec<Row>,
}

impl Report {
    fn worst_delta(&self) -> f64 {
        self.rows.iter().map(Row::delta).fold(0.0, f64::max)
    }

    /// The median per-benchmark δ, for the summary line: the worst alone says nothing about the
    /// host's typical behaviour, and on a machine with one bad benchmark the two differ a lot.
    fn median_delta(&self) -> f64 {
        let mut d: Vec<f64> = self.rows.iter().map(Row::delta).collect();
        d.sort_by(f64::total_cmp);
        d.get(d.len() / 2).copied().unwrap_or(0.0)
    }

    fn text(&self, floor: Option<f64>) -> String {
        let mut out =
            String::from("  benchmark                                    low  point   high\n");
        let mut rows: Vec<&Row> = self.rows.iter().collect();
        rows.sort_by(|a, b| b.ci.point.total_cmp(&a.ci.point));
        for r in rows {
            let flag = match floor {
                Some(f) if r.fails(f) => "  FAIL",
                _ => "",
            };
            out.push_str(&format!(
                "  {:<40} {:6.4} {:6.4} {:6.4}{flag}\n",
                r.name, r.ci.low, r.ci.point, r.ci.high
            ));
        }
        out
    }
}

/// The CI per benchmark the two runs share. A benchmark in one run and not the other is an error:
/// the two passes did not measure the same set, so no verdict covers it.
fn compare(
    base: &BTreeMap<String, samples::Samples>,
    new: &BTreeMap<String, samples::Samples>,
) -> Result<Report, String> {
    if base.keys().ne(new.keys()) {
        return Err(
            "the two passes benched different sets; re-run with a clean target/criterion"
                .to_string(),
        );
    }
    let mut rows = Vec::with_capacity(base.len());
    for (name, b) in base {
        let n = &new[name];
        let ci = ratio_ci(&b.per_iter, &n.per_iter, SEED).ok_or_else(|| {
            format!(
                "{name}: {} and {} samples do not pair",
                b.per_iter.len(),
                n.per_iter.len()
            )
        })?;
        rows.push(Row {
            name: name.clone(),
            ci,
        });
    }
    Ok(Report { rows })
}

/// One recorded `--aa` run. The floor is the maximum over every run `HOST.md` holds.
#[derive(Debug, Clone, PartialEq)]
struct Aa {
    benchmarks: usize,
    median: f64,
    max: f64,
    cpu: String,
}

/// What this host showed when it was last measured: the most recent run's worst `δ`.
///
/// **Not an allowance, and not a maximum over the runs.** Accumulating the maximum was tried and
/// the measurement refuted it in one run: a session whose worst `δ` was 0.6864 would pin the floor
/// there forever, and a gate with a 68% allowance passes a 50% regression. The most recent run is
/// the only one that describes the host now, and even that describes it only for as long as the
/// machine's state holds — which is why `--against` measures its own control (see `against`).
fn last_max(runs: &[Aa]) -> f64 {
    runs.last().map_or(0.0, |r| r.max)
}

/// The allowance an interval justifies, **invariant under swapping the two runs**.
///
/// Which of a pair is "baseline" is a slot name, so a floor read off `high` alone charges the pair
/// for drift one way and nothing for the same drift the other way. Relabelling replaces the
/// interval by its reciprocal, so the larger of the two one-sided reaches is the allowance, and on
/// this host taking `high` alone gave the four noisiest pairs — 22% to 43% apart — the tightest
/// floors in the run.
fn delta(ci: &Interval) -> f64 {
    let reciprocal = if ci.low > 0.0 {
        1.0 / ci.low - 1.0
    } else {
        0.0
    };
    (ci.high - 1.0).max(reciprocal).max(0.0)
}

/// `baseline/HOST.md`: the floor, the runs it is the maximum over, and the host they ran on.
///
/// `previous` is what `read_host` found, so a bless appends. Runs measured on a different CPU are
/// dropped: their floor is another machine's.
fn write_host(path: &PathBuf, report: &Report, previous: &[Aa]) -> Result<Vec<Aa>, String> {
    let parent = path.parent().ok_or("HOST.md has no parent")?;
    std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    let cpu = host_line();
    let mut runs: Vec<Aa> = previous.iter().filter(|r| r.cpu == cpu).cloned().collect();
    runs.push(Aa {
        benchmarks: report.rows.len(),
        median: report.median_delta(),
        max: report.worst_delta(),
        cpu,
    });
    let floor = last_max(&runs);

    let mut history =
        String::from("| run | benchmarks | median `δ` | max `δ` | CPU |\n|---|---|---|---|---|\n");
    for (i, r) in runs.iter().enumerate() {
        history.push_str(&format!(
            "| {} | {} | {:.4} | {:.4} | {} |\n",
            i + 1,
            r.benchmarks,
            r.median,
            r.max,
            r.cpu
        ));
    }
    let mut per_bench = String::from("| benchmark | `δ` |\n|---|---|\n");
    let mut rows: Vec<&Row> = report.rows.iter().collect();
    rows.sort_by(|a, b| a.name.cmp(&b.name));
    for r in rows {
        per_bench.push_str(&format!("| `{}` | {:.4} |\n", r.name, r.delta()));
    }

    let body = format!(
        "<!-- @generated by `cargo xtask bench-gate` — do not edit; `--aa --bless` writes it. -->\n\
         # What this host's A/A noise has measured\n\n\
         `δ` is an **A/A** spread (`docs/PHASE1.md` §9): the largest slowdown the paired-bootstrap\n\
         95% CI reaches when the *same* code is run against itself, so it is this machine and never\n\
         a change. Measured, never asserted.\n\n\
         **This file is a log, not the gate's allowance.** §9 records `δ` here and `--against`\n\
         deliberately does not read it, because on this host the number does not survive the walk\n\
         from one session to the next — see the table. The allowance `--against` uses is a control\n\
         measured beside the comparison, in the same invocation, on the same benchmark\n\
         (`0033`, draft). What this log is for: knowing whether the machine is quiet enough to bother.\n\n\
         Last measured: **{floor:.4} ({pct:.2}%)** worst, over {n} recorded run(s).\n\n{history}\n\
         **A stored floor does not reproduce, and the spread is not subtle.** Two sessions a few\n\
         hours apart, same protocol, same binary against itself: worst `δ` 0.0161, then 0.6864; p75\n\
         0.0036, then 0.0832; 0 of 60 benchmarks above 2%, then 23. A floor from the quiet session\n\
         fails the identical binary against itself, and one from the noisy session would pass a 50%\n\
         regression. Nothing in the protocol changed between them; the machine's state did, and a\n\
         shared VM does not announce that (`available_parallelism` read 1 at one bless and 8 at\n\
         another).\n\n\
         **Neither one floor per benchmark nor an accumulated maximum survives that.** Per\n\
         benchmark was measured first: across two quiet sessions the *distribution* over the 60\n\
         reproduces (min 0.0008 then 0.0009, median 0.0029 then 0.0030) while single benchmarks do\n\
         not — `exp_coeffs_straddle-0.9` read 0.0019, then 0.0190 — and 8 of 60 then failed the\n\
         identical binary. Accumulating the maximum over runs instead was refuted by the first run\n\
         after it: one noisy session pins the allowance at 68% permanently. The per-benchmark table\n\
         below is evidence of where the noise is, and is not a set of floors.\n\n\
         **The statistic is the ratio of medians, from the least disturbed of {reps} alternations.**\n\
         Contamination here is one-sided at both levels — a disturbed sample is slow and never\n\
         fast, and so is a disturbed window — so both take a low order statistic. Measured on this\n\
         host: one benchmark over eight consecutive windows spread 2.97% by window mean and 1.04%\n\
         by window median, while *inside* a window the interquartile spread was 1.1% to 2.7% and\n\
         the extremes reached 96%. Taking the mean, then letting a bootstrap resample those extremes\n\
         with replacement, gave a worst `δ` of 71.72%; the median gives 1.90% on the same protocol.\n\
         Pooling three pairs did not finish the job — a disturbance covering two consecutive runs\n\
         puts two of a side's three windows over, past a median's breakdown point, and left one\n\
         benchmark at 58.84% where ten consecutive windows of it spread 2.38%. Reporting the\n\
         cheapest pair instead survives two disturbed pairs of the three.\n\n\
         **`δ` is `max(high - 1, 1/low - 1)`, not `high - 1`.** Which run of a pair is called the\n\
         baseline is a slot name; relabelling replaces the interval by its reciprocal. Taking only\n\
         `high` charges a pair for drift one way and nothing for the same drift the other way, and\n\
         on this host that gave the four noisiest pairs — 22% to 43% apart — the tightest floors\n\
         in the run.\n\n\
         **The two runs of a pair are adjacent in time.** Measured on this host: one benchmark run\n\
         four times in a row spreads 0.6%, while the same benchmark measured once per pass of all\n\
         60 — the two measurements about two minutes apart — spreads up to 90%. Thermal and\n\
         frequency state drift over a pass and do not cancel, and no amount of resampling recovers\n\
         it. A floor measured the way the gate measures is the only one that means anything.\n\n\
         Every measurement in every mode runs the bench **binary** directly, built once. One side\n\
         under `cargo` and the other direct is not symmetric, and a floor measured both-under-cargo\n\
         does not cover it: that asymmetry alone put 3 of 60 identical-code verdicts past their\n\
         floor, at 1.007 to 1.018. {RESAMPLES} resamples. `taskset` is the caller's: the floor is\n\
         only as quiet as the machine was.\n\n\
         ## Per benchmark, this run\n\n{per_bench}\n\
         ## The A/A run that set them\n\n```text\n{table}```\n",
        floor = floor,
        pct = floor * 100.0,
        n = runs.len(),
        history = history,
        reps = REPLICATES,
        per_bench = per_bench,
        table = report.text(None),
    );
    std::fs::write(path, body).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(runs)
}

/// The recorded A/A runs, in the order they were blessed.
fn read_host(path: &PathBuf) -> Result<Vec<Aa>, String> {
    let text = std::fs::read_to_string(path).map_err(|_| {
        format!(
            "no {}; measure this host's floor with `cargo xtask bench-gate --aa --bless`",
            path.display()
        )
    })?;
    let mut out = Vec::new();
    for line in text.lines() {
        // `| <run> | <benchmarks> | <median> | <max> | <cpu> |`, the history table and no other
        // row: the header, the separator and the two-column tables all fail to parse.
        let cells: Vec<&str> = line.trim().split('|').collect();
        let [_, run, benchmarks, median, max, cpu, _] = cells.as_slice() else {
            continue;
        };
        let (Ok(_), Ok(benchmarks), Ok(median), Ok(max)) = (
            run.trim().parse::<usize>(),
            benchmarks.trim().parse::<usize>(),
            median.trim().parse::<f64>(),
            max.trim().parse::<f64>(),
        ) else {
            continue;
        };
        out.push(Aa {
            benchmarks,
            median,
            max,
            cpu: cpu.trim().to_string(),
        });
    }
    if out.is_empty() {
        return Err(format!(
            "{}: no recorded A/A run; re-measure with `cargo xtask bench-gate --aa --bless`",
            path.display()
        ));
    }
    Ok(out)
}

/// The host, for `HOST.md` and for saying which machine a baseline was blessed on: a committed
/// timing baseline is not portable, and the message is the only place that says so at bless time.
fn host_line() -> String {
    let cpu = read_first("/proc/cpuinfo", "model name").unwrap_or_else(|| "unknown".into());
    let cores = std::thread::available_parallelism()
        .map(|n| n.get().to_string())
        .unwrap_or_else(|_| "unknown".into());
    format!("{cpu}, {cores} logical CPUs")
}

fn read_first(file: &str, key: &str) -> Option<String> {
    let text = std::fs::read_to_string(file).ok()?;
    text.lines()
        .find(|l| l.starts_with(key))
        .and_then(|l| l.split_once(':'))
        .map(|(_, v)| v.trim().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(a: &[&str]) -> Vec<String> {
        a.iter().map(ToString::to_string).collect()
    }

    /// The three flags, `--against <path>`, and nothing else. `--bless` means two different
    /// artifacts by mode — the floor with `--aa`, the committed samples without — so it is valid
    /// either way.
    #[test]
    fn arguments_are_refused_unless_the_mode_combination_is_meaningful() -> Result<(), String> {
        let o = parse_args(&args(&["--aa", "--bless", "--dry-run"]))?;
        assert!(o.aa && o.bless && o.dry_run && o.against.is_none());
        let o = parse_args(&args(&["--bless"]))?;
        assert!(o.bless && !o.aa && !o.dry_run);
        let o = parse_args(&args(&["--against", "/tmp/coeffs"]))?;
        assert_eq!(o.against, Some(PathBuf::from("/tmp/coeffs")));
        assert!(parse_args(&args(&[])).is_ok());
        assert!(parse_args(&args(&["--nope"])).is_err());
        assert!(parse_args(&args(&["--aa=1"])).is_err());
        // A path is required, and the two comparison modes are exclusive: `--aa` runs one binary
        // against itself, so a second binary has no role in it.
        assert!(parse_args(&args(&["--against"])).is_err());
        assert!(parse_args(&args(&["--aa", "--against", "/tmp/c"])).is_err());
        // `--record` persists what a measurement produced, so it needs one; `--replay` measures
        // nothing, so it takes no other mode.
        let o = parse_args(&args(&["--against", "/tmp/c", "--record", "/tmp/r"]))?;
        assert_eq!(o.record, Some(PathBuf::from("/tmp/r")));
        let o = parse_args(&args(&["--replay", "/tmp/r"]))?;
        assert_eq!(o.replay, Some(PathBuf::from("/tmp/r")));
        assert!(parse_args(&args(&["--record", "/tmp/r"])).is_err());
        assert!(parse_args(&args(&["--record"])).is_err());
        assert!(parse_args(&args(&["--replay"])).is_err());
        assert!(parse_args(&args(&["--replay", "/tmp/r", "--aa"])).is_err());
        assert!(parse_args(&args(&["--replay", "/tmp/r", "--bless"])).is_err());
        assert!(parse_args(&args(&["--replay", "/tmp/r", "--against", "/tmp/c"])).is_err());
        Ok(())
    }

    fn row(name: &str, low: f64, point: f64, high: f64) -> Row {
        Row {
            name: name.to_string(),
            ci: Interval { low, point, high },
        }
    }

    /// §9's rule, exactly: the **whole** CI must lie above `1 + δ`. A CI that straddles it passes,
    /// however high its point estimate, because the measurement does not separate it from noise.
    #[test]
    fn only_a_ci_wholly_above_the_floor_fails() {
        let floor = 0.05;
        assert!(row("slow", 1.20, 1.30, 1.40).fails(floor));
        assert!(!row("straddles", 1.00, 1.30, 1.60).fails(floor));
        assert!(!row("at the floor", 1.05, 1.08, 1.10).fails(floor));
        assert!(!row("faster", 0.80, 0.85, 0.90).fails(floor));
    }

    /// The median δ summarises the host where the worst alone would not: one bad benchmark moves
    /// the worst and not the middle.
    #[test]
    fn the_median_delta_is_not_moved_by_one_bad_benchmark() {
        let r = Report {
            rows: vec![
                row("a", 0.99, 1.00, 1.01),
                row("b", 0.99, 1.00, 1.02),
                row("c", 1.80, 1.85, 1.90),
            ],
        };
        assert!(
            (r.median_delta() - 0.02).abs() < 1e-12,
            "{}",
            r.median_delta()
        );
        assert!((r.worst_delta() - 0.90).abs() < 1e-12);
    }

    /// The floor does not depend on which run of the pair is called the baseline: the two
    /// labellings give reciprocal intervals and the same δ.
    #[test]
    fn the_floor_is_the_same_under_either_labelling_of_the_pair() {
        let slower = row("a", 1.05, 1.20, 1.40);
        let faster = row(
            "a",
            1.0 / slower.ci.high,
            1.0 / slower.ci.point,
            1.0 / slower.ci.low,
        );
        assert!(
            (slower.delta() - faster.delta()).abs() < 1e-12,
            "{} vs {}",
            slower.delta(),
            faster.delta()
        );
        // And it is the larger of the two deviations, which here is the interval's own upper end.
        assert!((slower.delta() - 0.40).abs() < 1e-12, "{}", slower.delta());
        // A pair that came out 29% faster is 1.41x apart the other way round, not noiseless.
        assert!(
            (row("b", 0.71, 0.74, 0.78).delta() - 0.4085).abs() < 1e-4,
            "{}",
            row("b", 0.71, 0.74, 0.78).delta()
        );
    }

    /// The floor is the largest per-benchmark δ, and never negative: a degenerate interval leaves
    /// no allowance rather than a negative one.
    #[test]
    fn the_floor_is_the_largest_delta_and_never_negative() {
        let r = Report {
            rows: vec![row("a", 0.9, 1.0, 1.02), row("b", 0.95, 1.01, 1.07)],
        };
        // `a`'s low of 0.9 is 1.111x the other way round, which is the larger deviation in the run.
        assert!(
            (r.worst_delta() - 0.1111).abs() < 1e-4,
            "{}",
            r.worst_delta()
        );
        let unit = Report {
            rows: vec![row("a", 1.0, 1.0, 1.0)],
        };
        assert_eq!(unit.worst_delta().to_bits(), 0.0_f64.to_bits());
        // A non-positive bound cannot be inverted; it leaves no allowance either.
        assert_eq!(row("a", 0.0, 0.0, 0.0).delta().to_bits(), 0.0_f64.to_bits());
    }

    /// A HOST.md round trip: the floor the gate reads is the one `--bless` wrote, and a second
    /// bless appends rather than replacing, so the floor is the maximum over both runs.
    #[test]
    fn the_floor_round_trips_through_host_md() -> Result<(), String> {
        let dir = std::env::temp_dir().join("helicoid-bench-host");
        let _ = std::fs::remove_dir_all(&dir);
        let path = dir.join("HOST.md");
        let loud = Report {
            rows: vec![
                row("coeffs_f64/exp", 0.98, 1.01, 1.0423),
                // Wholly below 1: the allowance is what the pair shows when relabelled.
                row("coeffs_f64/jr", 0.90, 0.95, 0.9800),
            ],
        };
        let quiet = Report {
            rows: vec![row("coeffs_f64/exp", 0.999, 1.001, 1.0020)],
        };
        let runs = write_host(&path, &loud, &[])?;
        assert_eq!(runs.len(), 1);
        // Read back at the file's four decimals, which is the precision the gate gates on.
        let back = read_host(&path)?;
        assert_eq!(back.len(), 1, "{back:?}");
        assert_eq!(back[0].cpu, runs[0].cpu);
        // The loudest of the two benchmarks sets the run's max: 1/0.90 - 1, not 1.0423 - 1.
        assert!((last_max(&back) - 0.1111).abs() < 1e-4, "{back:?}");
        assert_eq!(back[0].benchmarks, 2);

        // A quieter second run appends, and the log reports the newest run: a number measured in
        // an earlier session does not describe the machine now, which is why nothing gates on it.
        let runs = write_host(&path, &quiet, &back)?;
        assert_eq!(runs.len(), 2);
        let back = read_host(&path)?;
        assert_eq!(back.len(), 2, "{back:?}");
        assert!((last_max(&back) - 0.0020).abs() < 1e-4, "{back:?}");
        assert!((back[0].max - 0.1111).abs() < 1e-4, "{back:?}");

        // A run from another CPU is another machine's floor: the table starts over.
        let elsewhere = vec![Aa {
            benchmarks: 60,
            median: 0.5,
            max: 0.9,
            cpu: "a different machine".to_string(),
        }];
        let runs = write_host(&path, &quiet, &elsewhere)?;
        assert_eq!(runs.len(), 1, "{runs:?}");
        assert!((last_max(&runs) - 0.0020).abs() < 1e-4, "{runs:?}");

        // A file with no recorded run, and a missing file, are both errors rather than a default
        // floor: a gate with no measured allowance is the asserted percentage §9 refuses.
        std::fs::write(&path, "# nothing here\n").map_err(|e| e.to_string())?;
        assert!(read_host(&path).is_err());
        assert!(read_host(&dir.join("absent.md")).is_err());
        Ok(())
    }

    /// The log reports its newest run and never a maximum over them: accumulating the maximum
    /// lets one noisy session speak for the host forever, which a 0.6864 run demonstrated.
    #[test]
    fn the_log_reports_its_newest_run() {
        let aa = |max: f64| Aa {
            benchmarks: 60,
            median: 0.003,
            max,
            cpu: "x".to_string(),
        };
        assert_eq!(last_max(&[]).to_bits(), 0.0_f64.to_bits());
        assert!((last_max(&[aa(0.0161), aa(0.6864), aa(0.0150)]) - 0.0150).abs() < 1e-12);
    }

    fn ci(low: f64, point: f64, high: f64) -> Interval {
        Interval { low, point, high }
    }

    /// A benchmark fails only when **every** pair of every replicate is wholly above the floor
    /// measured beside it. One disturbed window is what this rule exists to absorb.
    #[test]
    fn a_regression_must_appear_in_every_pair() {
        let v = |floor: f64, pairs: Vec<Interval>| Verdict {
            name: "coeffs_f64/jr_coeffs_near-pi".to_string(),
            floor,
            pairs,
        };
        // Every pair above the floor: a regression.
        let all = v(
            0.01,
            vec![
                ci(1.04, 1.05, 1.06),
                ci(1.03, 1.05, 1.07),
                ci(1.05, 1.06, 1.07),
            ],
        );
        assert!(all.fails());
        assert!((all.claim() - 1.03).abs() < 1e-12);
        // One disturbed window and two clean ones: not a verdict.
        let one = v(
            0.01,
            vec![
                ci(1.40, 1.50, 1.60),
                ci(0.99, 1.00, 1.01),
                ci(0.98, 1.00, 1.02),
            ],
        );
        assert!(!one.fails());
        // The same intervals under a floor that covers them: the machine, not the change.
        assert!(!v(0.10, vec![ci(1.04, 1.05, 1.06), ci(1.03, 1.05, 1.07)]).fails());
        // A CI that straddles its floor is not wholly above it (§9).
        assert!(!v(0.01, vec![ci(1.00, 1.30, 1.60), ci(1.02, 1.30, 1.60)]).fails());
        // Faster is never a failure, and no pair is never a verdict.
        assert!(!v(0.01, vec![ci(0.80, 0.85, 0.90)]).fails());
        assert!(!v(0.01, vec![]).fails());
    }

    /// A clean A/A triplet: deterministic jitter, the candidate scaled by `scale`.
    fn synthetic(name: &str, seed: usize, scale: f64) -> Windows {
        let jitter = |k: usize| 10.0 + f64::from(u8::try_from((k + seed) % 7).unwrap_or(0)) / 100.0;
        let window = |o: usize| (0..100).map(|i| jitter(i + o)).collect::<Vec<f64>>();
        Windows {
            name: name.to_string(),
            first: window(0),
            candidate: window(1).iter().map(|x| x * scale).collect(),
            second: window(2),
        }
    }

    /// `REPLICATES` triplets of one benchmark, each replicate's candidate scaled by its own factor.
    fn run_of(scales: &[f64]) -> Vec<Windows> {
        scales
            .iter()
            .enumerate()
            .map(|(r, &s)| synthetic("coeffs_f64/exp_coeffs_generic-1", r, s))
            .collect()
    }

    /// The seam: a recording in criterion's own format replays to the verdict a measurement gave.
    ///
    /// This is what makes `0035`'s (draft) first layer possible — the rule is a function of these
    /// vectors, so it is exercised without a quiet machine and without a twenty-minute run. The
    /// fixture is written the way criterion writes it, so `samples::read_file`'s division and
    /// validation are on the replay path and not merely trusted.
    #[test]
    fn a_recording_round_trips_and_replays_to_the_same_verdict() -> Result<(), String> {
        let dir = std::env::temp_dir().join("helicoid-bench-recording");
        let _ = std::fs::remove_dir_all(&dir);
        let windows = run_of(&[1.0, 1.0, 1.0]);
        for (r, w) in windows.iter().enumerate() {
            let slots = [&w.first, &w.candidate, &w.second];
            for (slot, per_iter) in TRIPLET_SLOTS.iter().zip(slots) {
                let into = dir.join(&w.name).join(format!("rep{}", r + 1)).join(slot);
                std::fs::create_dir_all(&into).map_err(|e| e.to_string())?;
                // `iters` integral and `times` their product, so the parser's quotient returns the
                // per-iteration time this triplet stands for.
                let iters: Vec<f64> = (1..=per_iter.len())
                    .map(|i| f64::from(u32::try_from(i * 1000).unwrap_or(u32::MAX)))
                    .collect();
                let times: Vec<f64> = iters.iter().zip(per_iter).map(|(n, p)| n * p).collect();
                let body = format!(
                    "{{\"sampling_mode\":\"Linear\",\"iters\":{iters:?},\"times\":{times:?}}}"
                );
                std::fs::write(into.join(SAMPLE), body).map_err(|e| e.to_string())?;
            }
        }
        let back = read_recording(&dir)?;
        assert_eq!(back.len(), windows.len(), "{back:?}");
        let replayed = accumulate(&back)?;
        assert_eq!(replayed.len(), 1);
        assert_eq!(replayed[0].pairs.len(), 2 * REPLICATES);
        assert_eq!(replayed[0].name, "coeffs_f64/exp_coeffs_generic-1");
        assert!(!replayed[0].fails(), "identical code must not fail");
        // An empty directory is an error that says how to make a recording.
        let empty = std::env::temp_dir().join("helicoid-bench-recording-empty");
        let _ = std::fs::remove_dir_all(&empty);
        std::fs::create_dir_all(&empty).map_err(|e| e.to_string())?;
        let e = read_recording(&empty).err().unwrap_or_default();
        assert!(e.contains("--record"), "{e}");
        Ok(())
    }

    /// One disturbed window of nine is not a verdict: the other triplets disagree.
    ///
    /// `0033` (draft) paid a twenty-minute run to find this class of thing. It is now arithmetic.
    #[test]
    fn one_disturbed_window_is_not_a_verdict() -> Result<(), String> {
        let mut windows = run_of(&[1.0, 1.0, 1.0]);
        for x in windows[0].candidate.iter_mut() {
            *x *= 1.60;
        }
        let v = accumulate(&windows)?;
        assert!(!v[0].fails(), "{:?} floor {}", v[0].pairs, v[0].floor);
        Ok(())
    }

    /// A disturbance covering one whole triplet is not a verdict either, and this is the case the
    /// replicate loop's order exists for: when a benchmark's triplets ran consecutively, one such
    /// span covered all six of its pairs and 2 of 60 identical-code benchmarks failed.
    #[test]
    fn a_disturbance_spanning_one_triplet_is_not_a_verdict() -> Result<(), String> {
        let mut windows = run_of(&[1.0, 1.0, 1.0]);
        for x in windows[0].candidate.iter_mut() {
            *x *= 1.30;
        }
        let v = accumulate(&windows)?;
        assert!(!v[0].fails(), "{:?} floor {}", v[0].pairs, v[0].floor);
        Ok(())
    }

    /// Drift is separated from a regression, which is what bracketing buys: every window of every
    /// triplet slow by 30% is a slow machine, and the ratios stay at 1.
    #[test]
    fn a_uniformly_slow_run_is_drift_and_not_a_regression() -> Result<(), String> {
        let mut windows = run_of(&[1.0, 1.0, 1.0]);
        for w in windows.iter_mut() {
            for v in [&mut w.first, &mut w.candidate, &mut w.second] {
                for x in v.iter_mut() {
                    *x *= 1.30;
                }
            }
        }
        let v = accumulate(&windows)?;
        assert!(!v[0].fails(), "{:?} floor {}", v[0].pairs, v[0].floor);
        Ok(())
    }

    /// A slowdown present in every triplet is a verdict: that is the gate's whole purpose.
    ///
    /// It is also the rule's limit, stated here because nothing in the data separates the two: a
    /// disturbance confined to the candidate's windows in *every* triplet has this same shape. The
    /// floor's own size is the only report of that, which is why `--against` prints it every run.
    #[test]
    fn a_slowdown_in_every_triplet_is_a_verdict() -> Result<(), String> {
        let v = accumulate(&run_of(&[1.05, 1.05, 1.05]))?;
        assert!(v[0].fails(), "{:?} floor {}", v[0].pairs, v[0].floor);
        assert!(v[0].claim() > 1.0);
        Ok(())
    }

    /// The report orders by the middle pair's ratio, so the noisiest single window cannot put a
    /// benchmark at the top of the table.
    #[test]
    fn the_report_orders_by_the_middle_pair() {
        let v = Verdict {
            name: "g/b".to_string(),
            floor: 0.01,
            pairs: vec![
                ci(1.90, 2.00, 2.10),
                ci(0.99, 1.00, 1.01),
                ci(1.00, 1.01, 1.02),
            ],
        };
        assert!((v.point() - 1.01).abs() < 1e-12, "{}", v.point());
        let text = table(&[v]);
        assert!(text.contains("g/b"), "{text}");
        assert!(!text.contains("FAIL"), "{text}");
    }

    /// Two passes that benched different sets have no verdict: the pairing is by name.
    #[test]
    fn a_mismatched_pair_of_passes_is_an_error() {
        let s = |v: Vec<f64>| samples::Samples { per_iter: v };
        let mut base = BTreeMap::new();
        base.insert("g/a".to_string(), s(vec![1.0, 2.0]));
        let mut new = BTreeMap::new();
        new.insert("g/b".to_string(), s(vec![1.0, 2.0]));
        assert!(compare(&base, &new).is_err());
        // Same names, unequal sample counts: not pairs either.
        let mut new = BTreeMap::new();
        new.insert("g/a".to_string(), s(vec![1.0]));
        assert!(compare(&base, &new).is_err());
    }
}
