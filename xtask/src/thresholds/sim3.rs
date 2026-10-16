//! The Sim(3) coefficient group's sweep (`docs/decisions/0070` decision 4): `K` and `r_×²` of the
//! total-degree series, `Δ_b` of `v₁`'s form b, `σ_w²` and the terms of `w₀`'s series, at each
//! precision, from `conformance/corpus/sim3_ref.jsonl` (the eight scalars on rays in `(σ, θ)`) and
//! the arms `helicoid::coeffs` ships, through `__sweep`. Part of the `helicoid` target: it writes
//! `conformance/sweeps/thresholds-sim3.csv` and `crates/helicoid/src/coeffs/generated_sim3.rs`.
//!
//! **The objective is the consumer's.** At a point, a scalar's error `|ĉ - c|` is weighted by the
//! power of `θ` that multiplies it in `V`, `Q` and `Γ₂(M)ρ` — `θ⁰` for `v₀, v₁, w₀` (`v₁` is
//! `Q`'s `-v₁Ξ`, which nothing suppresses), `θ¹` for `v₂, w₁`, `θ²` for `w₂, ∂_z v₁`, `θ³` for
//! `∂_z v₂` (`docs/maths/sim3.md` SM.7(c), SM.11(e)) — over `u` times the block's scale
//! `max(1, v₀)`; a point scores the worst of its eight. A candidate's objective is the maximum over
//! the points of the arm it selects there.
//!
//! **The search.** `w₀` first, on its own error: `terms` in `1..=W0_TERMS_MAX` and `σ_w²` on the
//! switch grid, the series below. Then each point's closed-arm score for each form of `v₁` and each
//! depth's series score; for each `r_×²` on the grid the best `Δ_b ≥ r_×²` by prefix maxima over
//! the points sorted by `Δ`, and the pair `(K, r_×²)` minimizing the larger of the two arms.
//! Candidates within [`SLACK`] of the best are equal, and the cheapest wins: fewer terms, then the
//! larger series region for `(K, r_×²)` and the smaller one for `w₀` (its exact arm reuses
//! `expm1 σ`, which the closed arm has already formed). Last, three shallower tiers of the series
//! arm, each from the one below it (`short_arm`, `0099`, draft).

use std::fmt::Write as _;
use std::path::Path;

use helicoid::__sweep as k;
use helicoid_linalg::{Precision, Real};
use num_bigint::{BigInt, BigUint, Sign};
use serde_json::Value;

use super::grid::grid;
use crate::conformance::number::{parse_hex, ratio_to_f64, Decimal};

pub(super) const REF_FILE: &str = "sim3_ref.jsonl";
pub(super) const SERIES_FILE: &str = "sim3_series.jsonl";
pub(super) const CSV: &str = "conformance/sweeps/thresholds-sim3.csv";
pub(super) const GENERATED: &str = "crates/helicoid/src/coeffs/generated_sim3.rs";

const NAMES: [&str; 8] = ["v0", "v1", "v2", "w0", "w1", "w2", "dv1", "dv2"];
/// The power of `θ` each scalar is multiplied by in the blocks.
const WEIGHT: [i32; 8] = [0, 0, 1, 0, 1, 2, 2, 3];
const K_RANGE: core::ops::RangeInclusive<usize> = 4..=24;
const W0_TERMS_MAX: usize = 20;
/// Two candidates whose objectives differ by less than this, in `u`, are equal: a quarter of a
/// unit, or 2% of the best, whichever is larger.
const SLACK: (f64, f64) = (0.25, 0.02);
/// How much worse than depth `K`, in `u`, a shorter series arm may score at any one point.
const TIER_SLACK: f64 = 0.0625;

/// A reference split exactly into `hi + lo`, two binary64s.
#[derive(Clone, Copy)]
struct Split {
    hi: f64,
    lo: f64,
}

struct Point {
    sigma: f64,
    z: f64,
    theta: f64,
    refs: [Split; 8],
    scale: f64,
}

/// The chosen parameters at one precision and what they score.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct Choice {
    pub(super) series_terms: usize,
    pub(super) series_below: f64,
    pub(super) form_b_below: f64,
    pub(super) w0_below: f64,
    pub(super) w0_terms: usize,
    pub(super) objective: f64,
    pub(super) series_objective: f64,
    pub(super) closed_objective: f64,
    pub(super) w0_objective: f64,
    pub(super) points: usize,
    pub(super) short_terms: usize,
    pub(super) short_below: f64,
    pub(super) mid_terms: usize,
    pub(super) mid_below: f64,
    pub(super) upper_terms: usize,
    pub(super) upper_below: f64,
}

fn signed(neg: bool, m: BigUint) -> BigInt {
    BigInt::from_biguint(if neg { Sign::Minus } else { Sign::Plus }, m)
}

/// `d` as the exact rational `num / den`.
fn rational(d: &Decimal) -> (BigInt, BigUint) {
    let ten = BigUint::from(10u32);
    let m = BigUint::from(d.mant);
    if d.exp10 >= 0 {
        (
            signed(d.neg, m * ten.pow(d.exp10.unsigned_abs())),
            BigUint::from(1u32),
        )
    } else {
        (signed(d.neg, m), ten.pow(d.exp10.unsigned_abs()))
    }
}

fn round(num: &BigInt, den: &BigUint) -> f64 {
    let x = ratio_to_f64(num.magnitude(), den);
    if num.sign() == Sign::Minus {
        -x
    } else {
        x
    }
}

/// `d = hi + lo` with `hi` its binary64 rounding and `lo` the rounding of the exact remainder: the
/// error of an output `y` is then `|(y - hi) - lo|`, exact to far below the error it measures.
fn split(d: &Decimal) -> Split {
    let (num, den) = rational(d);
    let hi = round(&num, &den);
    if hi == 0.0 {
        return Split { hi, lo: 0.0 };
    }
    let bits = hi.to_bits();
    let (exp, frac) = (((bits >> 52) & 0x7ff) as i64, bits & ((1 << 52) - 1));
    let (m, e) = if exp == 0 {
        (frac, -1074)
    } else {
        (frac | (1 << 52), exp - 1075)
    };
    let m = signed(hi < 0.0, BigUint::from(m));
    // num/den - m 2^e over one denominator.
    let (n2, d2) = if e >= 0 {
        (num - m * (BigInt::from(1) << e), den)
    } else {
        let p = BigUint::from(1u32) << (-e) as usize;
        (
            num * BigInt::from(p.clone()) - m * BigInt::from(den.clone()),
            den * p,
        )
    };
    Split {
        hi,
        lo: round(&n2, &d2),
    }
}

/// The reference points of `precision`.
fn points(dir: &Path, precision: Precision) -> Result<Vec<Point>, String> {
    let path = dir.join(REF_FILE);
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let want = match precision {
        Precision::F64 => "f64",
        Precision::F32 => "f32",
    };
    let mut out = Vec::new();
    for (n, line) in text.lines().enumerate() {
        let bad = |what: &str| format!("{}:{}: {what}", path.display(), n + 1);
        let v: Value = serde_json::from_str(line).map_err(|e| bad(&e.to_string()))?;
        if v["precision"] != want {
            continue;
        }
        let hex = |key: &str| v[key].as_str().ok_or_else(|| bad(key)).and_then(parse_hex);
        let (sigma, z) = (hex("sigma")?, hex("z")?);
        let mut refs = [Split { hi: 0.0, lo: 0.0 }; 8];
        for (r, name) in refs.iter_mut().zip(NAMES) {
            let s = v[name].as_str().ok_or_else(|| bad(name))?;
            *r = split(&s.parse::<Decimal>()?);
        }
        out.push(Point {
            sigma,
            z,
            theta: z.sqrt(),
            scale: refs[0].hi.abs().max(1.0),
            refs,
        });
    }
    Ok(out)
}

fn unit(precision: Precision) -> f64 {
    match precision {
        Precision::F64 => (-53f64).exp2(),
        Precision::F32 => (-24f64).exp2(),
    }
}

/// The weighted error of scalar `j` of `y` at `p`, in `u`; `inf` when not finite.
fn scalar(p: &Point, j: usize, y: f64, u: f64) -> f64 {
    if !y.is_finite() {
        return f64::INFINITY;
    }
    let r = p.refs[j];
    ((y - r.hi) - r.lo).abs() * p.theta.powi(WEIGHT[j]) / (u * p.scale)
}

fn score(p: &Point, y: &[f64; 8], u: f64) -> f64 {
    (0..8).fold(0.0, |w, j| w.max(scalar(p, j, y[j], u)))
}

fn within(x: f64, best: f64) -> bool {
    x <= best + SLACK.0.max(SLACK.1 * best)
}

fn eval<S: Real<Mask = bool> + Into<f64>>(
    f: impl Fn(S, S) -> [S; 8],
    p: &Point,
    to: impl Fn(f64) -> S,
) -> [f64; 8] {
    f(to(p.sigma), to(p.z)).map(Into::into)
}

/// The sweep at `S`, whose `from_f64` is exact on the points (an `@f32` point holds a binary32).
pub(super) fn sweep<S: Real<Mask = bool> + Into<f64>>(
    dir: &Path,
    precision: Precision,
    to: impl Fn(f64) -> S + Copy,
) -> Result<Choice, String> {
    let pts = points(dir, precision)?;
    let u = unit(precision);
    let grid = grid(precision);
    // `Δ` as the shipped mask forms it, at `S`.
    let delta: Vec<f64> = pts
        .iter()
        .map(|p| {
            let (s, z) = (to(p.sigma), to(p.z));
            (s * s + z).into()
        })
        .collect();

    // w₀: its own error, the series below `σ_w²`, the exact form above.
    let w0_err = |p: &Point, terms: usize| -> f64 {
        let y = eval::<S>(|s, z| k::sim3_closed(s, z, false, terms), p, to);
        scalar(p, 3, y[3], u)
    };
    let exact: Vec<f64> = pts.iter().map(|p| w0_err(p, 0)).collect();
    let sig2: Vec<f64> = pts.iter().map(|p| p.sigma * p.sigma).collect();
    let mut w0_best = (f64::INFINITY, 0, 0.0);
    let mut w0_cands = Vec::new();
    for terms in 1..=W0_TERMS_MAX {
        let ser: Vec<f64> = pts.iter().map(|p| w0_err(p, terms)).collect();
        for &s in &grid {
            let obj = (0..pts.len()).fold(0.0f64, |w, i| {
                w.max(if sig2[i] < s { ser[i] } else { exact[i] })
            });
            w0_cands.push((obj, terms, s));
            if obj < w0_best.0 {
                w0_best = (obj, terms, s);
            }
        }
    }
    // Cheapest within the slack: fewer terms, then the smaller series region.
    let (w0_objective, w0_terms, w0_below) = w0_cands
        .iter()
        .copied()
        .filter(|c| within(c.0, w0_best.0))
        .min_by(|a, b| a.1.cmp(&b.1).then(a.2.total_cmp(&b.2)))
        .unwrap_or(w0_best);

    // Each point's score on each arm.
    let order: Vec<usize> = {
        let mut o: Vec<usize> = (0..pts.len()).collect();
        o.sort_by(|&a, &b| delta[a].total_cmp(&delta[b]));
        o
    };
    let sorted_delta: Vec<f64> = order.iter().map(|&i| delta[i]).collect();
    let closed = |form_b: bool| -> Vec<f64> {
        order
            .iter()
            .map(|&i| {
                let p = &pts[i];
                if delta[i] <= 0.0 {
                    return f64::INFINITY;
                }
                let terms = if p.sigma * p.sigma < w0_below {
                    w0_terms
                } else {
                    0
                };
                score(
                    p,
                    &eval::<S>(|s, z| k::sim3_closed(s, z, form_b, terms), p, to),
                    u,
                )
            })
            .collect()
    };
    let (ca, cb) = (closed(false), closed(true));
    let n = pts.len();
    // Suffix maxima of form b's score: `cb_suffix[i]` is the worst over points `i..`.
    let mut cb_suffix = vec![0.0f64; n + 1];
    for i in (0..n).rev() {
        cb_suffix[i] = cb_suffix[i + 1].max(cb[i]);
    }
    let first_at = |x: f64| sorted_delta.partition_point(|&d| d < x);
    // For each candidate `r_×²`: the best closed objective over `Δ_b ≥ r_×²`, and that `Δ_b`.
    let mut best_closed = Vec::with_capacity(grid.len());
    for (gi, &s) in grid.iter().enumerate() {
        let i0 = first_at(s);
        let (mut run_a, mut best) = (0.0f64, (cb_suffix[i0], s));
        let mut i = i0;
        for &b in &grid[gi + 1..] {
            let ib = first_at(b);
            while i < ib {
                run_a = run_a.max(ca[i]);
                i += 1;
            }
            let obj = run_a.max(cb_suffix[ib]);
            if obj < best.0 {
                best = (obj, b);
            }
        }
        best_closed.push(best);
    }
    let mut cands = Vec::new();
    let mut best = f64::INFINITY;
    let mut sers = Vec::new();
    for k in K_RANGE {
        let ser: Vec<f64> = order
            .iter()
            .map(|&i| {
                score(
                    &pts[i],
                    &eval::<S>(|s, z| k::sim3_series(s, z, k), &pts[i], to),
                    u,
                )
            })
            .collect();
        let mut run = 0.0f64;
        let mut i = 0;
        for (gi, &s) in grid.iter().enumerate() {
            let is = first_at(s);
            while i < is {
                run = run.max(ser[i]);
                i += 1;
            }
            let (c, b) = best_closed[gi];
            let obj = run.max(c);
            cands.push((obj, k, s, b, run, c));
            best = best.min(obj);
        }
        sers.push((k, ser));
    }
    // Cheapest within the slack: fewer terms, then the larger series region.
    let chosen = cands
        .iter()
        .copied()
        .filter(|c| within(c.0, best))
        .min_by(|a, b| a.1.cmp(&b.1).then(b.2.total_cmp(&a.2)))
        .ok_or("sim3: no candidate")?;
    let (short_terms, short_below) =
        short_arm(&grid, &sorted_delta, &sers, chosen.1, (0.0, chosen.2));
    let (mid_terms, mid_below) = short_arm(
        &grid,
        &sorted_delta,
        &sers,
        chosen.1,
        (short_below, chosen.2),
    );
    let (upper_terms, upper_below) =
        short_arm(&grid, &sorted_delta, &sers, chosen.1, (mid_below, chosen.2));
    Ok(Choice {
        series_terms: chosen.1,
        series_below: chosen.2,
        form_b_below: chosen.3,
        w0_below,
        w0_terms,
        objective: chosen.0.max(w0_objective),
        series_objective: chosen.4,
        closed_objective: chosen.5,
        w0_objective,
        points: n,
        short_terms,
        short_below,
        mid_terms,
        mid_below,
        upper_terms,
        upper_below,
    })
}

/// A shorter series arm over `lo ≤ Δ < s` (`0099`, draft), `s` on the grid below `r_×²`: a depth
/// `K_s < K` that scores within [`TIER_SLACK`] of depth `K` at **every** point there, its
/// truncation negligible point by point and not only at the maximum (within [`SLACK`] of the
/// maximum, the 4- and 6-term `f32` arms read 1.04× `main`'s total `sim3_log` error). Of those, the fewest terms
/// summed over the points in `[lo, r_×²)`, then the larger region; `(K, lo)`, no arm, where none is
/// cheaper. Called with `lo = 0` for the short arm, then above each arm for the next, twice.
/// `sorted_delta` and each depth's scores are in the same `Δ` order.
fn short_arm(
    grid: &[f64],
    sorted_delta: &[f64],
    sers: &[(usize, Vec<f64>)],
    k_full: usize,
    (lo, series_below): (f64, f64),
) -> (usize, f64) {
    let first_at = |x: f64| sorted_delta.partition_point(|&d| d < x);
    let Some((_, full)) = sers.iter().find(|(k, _)| *k == k_full) else {
        return (k_full, lo);
    };
    let (start, end) = (first_at(lo), first_at(series_below));
    let mut best = (k_full, lo, (end - start) * k_full);
    for (k, ser) in sers.iter().filter(|(k, _)| *k < k_full) {
        let (mut worse, mut i) = (0.0f64, start);
        for &s in grid.iter().filter(|&&s| lo < s && s < series_below) {
            let is = first_at(s);
            while i < is {
                worse = worse.max(ser[i] - full[i]);
                i += 1;
            }
            let cost = (is - start) * k + (end - is) * k_full;
            if worse <= TIER_SLACK && (cost, -s) < (best.2, -best.1) {
                best = (*k, s, cost);
            }
        }
    }
    (best.0, best.1)
}

/// `1/n!` from `sim3_series.jsonl`, at binary64 and binary32.
pub(super) fn series(dir: &Path) -> Result<(Vec<f64>, Vec<f32>), String> {
    let path = dir.join(SERIES_FILE);
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let (mut wide, mut narrow) = (Vec::new(), Vec::new());
    for (n, line) in text.lines().enumerate() {
        let v: Value =
            serde_json::from_str(line).map_err(|e| format!("{SERIES_FILE}:{}: {e}", n + 1))?;
        let hex = |key: &str| {
            v[key]
                .as_str()
                .ok_or_else(|| format!("{SERIES_FILE}:{}: no `{key}`", n + 1))
                .and_then(parse_hex)
        };
        wide.push(hex("f64")?);
        narrow.push(hex("f32")? as f32);
    }
    Ok((wide, narrow))
}

/// Whether the compiled-in `1/n!` are the corpus's.
pub(super) fn series_is_current(dir: &Path) -> Result<bool, String> {
    let (wide, narrow) = series(dir)?;
    let (w, n) = k::sim3_recip_fact();
    let bits = |a: &[f64], b: &[f64]| {
        a.iter()
            .map(|x| x.to_bits())
            .eq(b.iter().map(|x| x.to_bits()))
    };
    let narrow64: Vec<f64> = narrow.iter().map(|&x| f64::from(x)).collect();
    let n64: Vec<f64> = n.iter().map(|&x| f64::from(x)).collect();
    Ok(bits(&wide, w) && bits(&narrow64, &n64))
}

fn lit64(x: f64) -> String {
    format!("f64::from_bits(0x{:016x})", x.to_bits())
}

fn lit32(x: f32) -> String {
    format!("f32::from_bits(0x{:08x})", x.to_bits())
}

/// `(path, text)` of each file a run writes.
type Files = Vec<(&'static str, String)>;

/// The two files: the CSV and `generated_sim3.rs`. `choices` is `None` for the placeholder a stale
/// series writes.
pub(super) fn render(dir: &Path, choices: Option<(&Choice, &Choice)>) -> Result<Files, String> {
    let (wide, narrow) = series(dir)?;
    let series_text = std::fs::read(dir.join(SERIES_FILE)).map_err(|e| e.to_string())?;
    let series_sha = super::emit::sha256_hex(&series_text);
    let mut csv = String::from(
        "precision,series_terms,series_below,form_b_below,w0_terms,w0_below,objective,series_objective,closed_objective,w0_objective,points,short_terms,short_below,mid_terms,mid_below,upper_terms,upper_below\n",
    );
    let mut rs = format!(
        "// @generated by `cargo xtask thresholds` from {CSV} — do not edit.\n\
         // Source series: conformance/corpus/{SERIES_FILE}, sha256 {series_sha}.\n\
         {}\n\
         use super::sim3::Sim3Switch;\n\n\
         /// `1/n!` for `n = 0..={}`, correctly rounded (`sim3_series.jsonl`).\n\
         pub(crate) const RECIP_FACT_F64: [f64; {}] = [\n",
        super::emit::libm_line(),
        wide.len() - 1,
        wide.len()
    );
    for x in &wide {
        let _ = writeln!(rs, "    {},", lit64(*x));
    }
    let _ = write!(rs, "];\n\n/// [`RECIP_FACT_F64`] at binary32.\npub(crate) const RECIP_FACT_F32: [f32; {}] = [\n", narrow.len());
    for x in &narrow {
        let _ = writeln!(rs, "    {},", lit32(*x));
    }
    rs.push_str("];\n");
    let placeholder = Choice {
        series_terms: 16,
        series_below: 0.5,
        form_b_below: 1.0,
        w0_below: 0.25,
        w0_terms: 16,
        objective: f64::NAN,
        series_objective: f64::NAN,
        closed_objective: f64::NAN,
        w0_objective: f64::NAN,
        points: 0,
        short_terms: 16,
        short_below: 0.0,
        mid_terms: 16,
        mid_below: 0.0,
        upper_terms: 16,
        upper_below: 0.0,
    };
    let (c64, c32) = choices.unwrap_or((&placeholder, &placeholder));
    for (name, c, wide) in [("f64", c64, true), ("f32", c32, false)] {
        let _ = writeln!(
            csv,
            "{name},{},{:e},{:e},{},{:e},{:e},{:e},{:e},{:e},{},{},{:e},{},{:e},{},{:e}",
            c.series_terms,
            c.series_below,
            c.form_b_below,
            c.w0_terms,
            c.w0_below,
            c.objective,
            c.series_objective,
            c.closed_objective,
            c.w0_objective,
            c.points,
            c.short_terms,
            c.short_below,
            c.mid_terms,
            c.mid_below,
            c.upper_terms,
            c.upper_below
        );
        let lit = |x: f64| if wide { lit64(x) } else { lit32(x as f32) };
        let ty = if wide { "f64" } else { "f32" };
        let _ = write!(
            rs,
            "\n// Objective (max u, the blocks' weighted error over {points} points): {obj:e}; \
             series {sobj:e}, closed {cobj:e}, w0 {wobj:e}.\n\
             pub(crate) const SIM3_{upper}: Sim3Switch<{ty}> = Sim3Switch {{\n    \
             // Δ < {sb:e}\n    series_below: {sbl},\n    series_terms: {k},\n    \
             // Δ < {shb:e}\n    short_below: {shbl},\n    short_terms: {sk},\n    \
             // Δ < {mb:e}\n    mid_below: {mbl},\n    mid_terms: {mk},\n    \
             // Δ < {ub:e}\n    upper_below: {ubl},\n    upper_terms: {uk},\n    \
             // Δ ≥ {fb:e}\n    form_b_below: {fbl},\n    \
             // σ² < {wb:e}\n    w0_below: {wbl},\n    w0_terms: {wt},\n}};\n",
            points = c.points,
            obj = c.objective,
            sobj = c.series_objective,
            cobj = c.closed_objective,
            wobj = c.w0_objective,
            upper = name.to_uppercase(),
            sb = c.series_below,
            sbl = lit(c.series_below),
            k = c.series_terms,
            shb = c.short_below,
            shbl = lit(c.short_below),
            sk = c.short_terms,
            mb = c.mid_below,
            mbl = lit(c.mid_below),
            mk = c.mid_terms,
            ub = c.upper_below,
            ubl = lit(c.upper_below),
            uk = c.upper_terms,
            fb = c.form_b_below,
            fbl = lit(c.form_b_below),
            wb = c.w0_below,
            wbl = lit(c.w0_below),
            wt = c.w0_terms,
        );
    }
    Ok(vec![(CSV, csv), (GENERATED, rs)])
}

/// The sweep at both precisions and its two files.
pub(super) fn files(dir: &Path) -> Result<(Files, [Choice; 2]), String> {
    let c64 = sweep::<f64>(dir, Precision::F64, |x| x)?;
    let c32 = sweep::<f32>(dir, Precision::F32, |x| x as f32)?;
    let files = render(dir, Some((&c64, &c32)))?;
    Ok((files, [c64, c32]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_split_is_exact() -> Result<(), String> {
        let d: Decimal = "1.00000000000000000000000000001e-1".parse()?;
        let s = split(&d);
        assert_eq!(s.hi.to_bits(), 0.1f64.to_bits());
        // 0.1f64 = 0.1000000000000000055511151231257827..., so lo ≈ -5.55e-18 + 1e-31.
        assert!((s.lo + 5.551_115_123_125_783e-18).abs() < 1e-30, "{}", s.lo);
        Ok(())
    }
}
