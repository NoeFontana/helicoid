//! The switch grid of `docs/PHASE1.md` §6: 64 points per decade of the branch variable `z = θ²`
//! (`s = n²/w²` for `r`), from `z = 1e-16` (`θ = 1e-8`) to `z = 10`, both ends included.
//!
//! Point `i` is `10^((i - 1024)/64)` correctly rounded to the precision it is swept at (binary32
//! points are held as the binary64 they equal), from the 64th root of a power of ten in integers:
//! the grid does not depend on the quality of any `pow`, and a binary32 point is never a binary64
//! one rounded again (`docs/decisions/0016` item 3).
//!
//! **The grid spans the domain and the selection rule stops below it**
//! (`docs/decisions/0039` item 7). `NUMERICS.md` §12 allows `z ∈ [0, π²]`, and a grid that stopped
//! at `z = 1` excluded the optimum for three of eight binary64 coefficients and six of seven at
//! binary32; a grid that stopped *above* `π²` let the sweep retire the exact arm, which measured
//! worse. So the grid is a geometric object that runs past the domain, and [`DOMAIN`] is a
//! semantic constraint the search applies — which keeps
//! `the_grid_has_64_points_per_decade_and_both_ends` a statement about the grid alone, and means a
//! later domain change does not reshape it.

use helicoid_linalg::Precision;
use num_bigint::BigUint;

use crate::conformance::number::{ratio_to_f32, ratio_to_f64};

/// Grid points per decade of `z`.
const PER_DECADE: usize = 64;
/// Decades of `z` below `1`: `θ` from `1e-8` to `1`.
const DECADES_BELOW: usize = 16;
/// Decades of `z` above `1`. One is enough: it contains `π²`, the top of `NUMERICS.md` §12's
/// domain, and [`DOMAIN`] stops the search strictly below that.
const DECADES_ABOVE: usize = 1;

/// The top of the branch variable's domain, `π²` (`NUMERICS.md` §12): the search refuses a switch
/// at or above it, because a series arm selected there would serve inputs the domain does not
/// contain while retiring the exact arm on inputs it does.
///
/// Typed, and not a `0004` quantity: `0004` governs a switch *point*, which is an accuracy choice
/// the sweep makes. This is the edge of the domain the coefficients are defined on, which no
/// measurement moves.
pub(super) const DOMAIN: f64 = core::f64::consts::PI * core::f64::consts::PI;

/// `10^(e/64)` rounded to nearest at `precision`. `m = ⌊2^P · 10^(e/64)⌋` is the integer 64th
/// root of `2^(64 P) · 10^e` (the floor of a root is the root of the floor); `(2m + 1) / 2^(P + 1)`
/// lies strictly between `m` and `m + 1` in units of `2^-P`, so the one rounding is never a tie: it
/// is the nearest value to the root, irrational or a power of ten (no dyadic).
///
/// The same construction on both sides of `1`, so a point above it is no less exact than one below
/// and a binary32 point is still binary32's own (`0016` item 3).
fn pow10_over_64(e: i32, precision: Precision) -> f64 {
    const P: usize = 128;
    let unit = BigUint::from(1u32) << (64 * P);
    let ten = |k: u32| BigUint::from(10u32).pow(k);
    let scaled = match e >= 0 {
        true => unit * ten(e.unsigned_abs()),
        false => unit / ten(e.unsigned_abs()),
    };
    let m = scaled.nth_root(64);
    rounded(
        &((m << 1usize) + 1u32),
        &(BigUint::from(1u32) << (P + 1)),
        precision,
    )
}

/// `num/den` rounded once to nearest at `precision`, held as a binary64.
fn rounded(num: &BigUint, den: &BigUint, precision: Precision) -> f64 {
    match precision {
        Precision::F64 => ratio_to_f64(num, den),
        Precision::F32 => f64::from(ratio_to_f32(num, den)),
    }
}

/// Whether `x` is a value `precision` holds, as every point of [`grid`] is at its own precision.
pub(super) fn holds(precision: Precision, x: f64) -> bool {
    match precision {
        Precision::F64 => true,
        Precision::F32 => f64::from(x as f32).to_bits() == x.to_bits(),
    }
}

/// The switch points, increasing: `1089` values from `1e-16` to `10`.
pub(super) fn grid(precision: Precision) -> Vec<f64> {
    let (lo, hi) = (
        (DECADES_BELOW * PER_DECADE) as i32,
        (DECADES_ABOVE * PER_DECADE) as i32,
    );
    (-lo..=hi).map(|e| pow10_over_64(e, precision)).collect()
}

/// The index of the largest grid point strictly below [`DOMAIN`]: the last switch the search may
/// choose, and the one it reports as the domain binding when it does
/// (`docs/decisions/0039` item 8).
pub(super) fn last_admissible(grid: &[f64]) -> usize {
    grid.iter().rposition(|&z| z < DOMAIN).unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::number::Dyadic;

    const BOTH: [Precision; 2] = [Precision::F64, Precision::F32];

    #[test]
    fn the_grid_has_64_points_per_decade_and_both_ends() -> Result<(), String> {
        for p in BOTH {
            let g = grid(p);
            assert_eq!(g.len(), 17 * 64 + 1);
            assert!(g.windows(2).all(|w| w[0] < w[1]));
            // Each decade starts at its power of ten, which Rust's parser rounds independently; the
            // ends, `1e-16` and `1e1`, are two of them.
            for j in -1i32..=16 {
                let text = format!("1e{}", -j);
                let ten = match p {
                    Precision::F64 => text.parse::<f64>().map_err(|e| format!("{e}"))?,
                    Precision::F32 => f64::from(text.parse::<f32>().map_err(|e| format!("{e}"))?),
                };
                let i = ((16 - j) * 64) as usize;
                assert_eq!(g[i].to_bits(), ten.to_bits(), "{p:?} 1e{}", -j);
            }
            // The domain bound falls inside the top decade, so a point above it exists and the
            // selection rule has something to refuse (`0039` item 7).
            let last = last_admissible(&g);
            assert!(g[last] < DOMAIN && g[last + 1] >= DOMAIN, "{p:?}");
            assert!(
                last < g.len() - 1,
                "{p:?}: the domain bound is not inside the grid"
            );
        }
        Ok(())
    }

    /// `10^(e/64)` lies between the midpoints from `x` to its neighbours `lo` and `hi`: `mid^64`
    /// against `10^e`, in integers over `2^1075` (`x` normal and at least `1e-16`). A negative `e`
    /// multiplies the left side, a positive one the right, so the comparison is exact on both
    /// sides of `1`.
    fn is_nearest(x: f64, (lo, hi): (f64, f64), e: i32) -> bool {
        let scaled = |y: f64| {
            let d = Dyadic::of(y);
            BigUint::from(d.mant) << (d.exp + 1074) as usize
        };
        let ten = |k: u32| BigUint::from(10u32).pow(k);
        let (up, down) = match e >= 0 {
            true => (ten(e.unsigned_abs()), ten(0)),
            false => (ten(0), ten(e.unsigned_abs())),
        };
        let power = |mid: BigUint| mid.pow(64) * &down;
        let one = (BigUint::from(1u32) << (64 * 1075)) * &up;
        let here = scaled(x);
        power(&scaled(lo) + &here) <= one && one <= power(&here + &scaled(hi))
    }

    /// Every grid point, at both precisions, has a reference record for every coefficient
    /// (`conformance/corpus/coeff_switch_ref.jsonl`, `docs/decisions/0039` plan step 0).
    ///
    /// This is what makes the artifact's key sound: it is built by a second implementation of this
    /// grid, in Python, and keyed by the point's binary64 bit pattern. A drift between the two
    /// constructions would otherwise show up only at whichever switch the sweep happened to
    /// choose — one missing key out of 1089, found by luck. Here it is 17 424 keys, found at once.
    #[test]
    fn every_grid_point_has_a_reference_record_for_every_coefficient() -> Result<(), String> {
        let path = crate::conformance::corpus_dir()?.join("coeff_switch_ref.jsonl");
        let text =
            std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        let mut keys = std::collections::BTreeSet::new();
        for line in text.lines() {
            let of = |name: &str| -> Option<String> {
                let at = line.find(&format!("\"{name}\":\""))?;
                let rest = &line[at + name.len() + 4..];
                Some(rest[..rest.find('"')?].to_string())
            };
            let (bits, coeff, precision) = (of("bits"), of("coeff"), of("precision"));
            match (bits, coeff, precision) {
                (Some(b), Some(c), Some(p)) => keys.insert((b, c, p)),
                _ => return Err(format!("a record with no key: `{line}`")),
            };
        }
        let names = ["k", "a", "b", "c", "d", "e", "cos_half", "r"];
        let mut want = 0;
        for (p, name) in [(Precision::F64, "f64"), (Precision::F32, "f32")] {
            for z in grid(p) {
                for coeff in names {
                    let key = (format!("{:016x}", z.to_bits()), coeff.into(), name.into());
                    assert!(keys.contains(&key), "no reference for {key:?} at z = {z:e}");
                    want += 1;
                }
            }
        }
        // And nothing else: a record whose point this grid does not hold is a drift too.
        assert_eq!(keys.len(), want);
        Ok(())
    }

    #[test]
    fn every_point_is_the_nearest_value_to_its_power_of_ten() {
        let f32s = |x: f64| {
            let y = x as f32;
            (f64::from(y.next_down()), f64::from(y.next_up()))
        };
        for (i, &x) in grid(Precision::F64).iter().enumerate() {
            let near = (x.next_down(), x.next_up());
            assert!(is_nearest(x, near, i as i32 - 1024), "f64 {i}: {x:e}");
        }
        // Binary32 points are on binary32's own spacing: not a binary64 point cast down.
        for (i, &x) in grid(Precision::F32).iter().enumerate() {
            assert!(holds(Precision::F32, x), "point {i} is not a binary32");
            assert!(is_nearest(x, f32s(x), i as i32 - 1024), "f32 {i}: {x:e}");
        }
    }

    #[test]
    fn a_binary32_point_is_rounded_once_from_its_rational_and_not_from_a_binary64_one() {
        // `1 + 2^-24 + 2^-70` is above the binary32 tie between 1 and `1 + 2^-23`. Binary64 holds
        // it as the tie itself, which a second rounding sends to the even 1. No point of the grid
        // is such a value, so its outputs cannot show a double rounding; this pins the rounding.
        let den = BigUint::from(1u32) << 70usize;
        let num = &den + (BigUint::from(1u32) << 46usize) + 1u32;
        let (once, wide) = (
            rounded(&num, &den, Precision::F32),
            rounded(&num, &den, Precision::F64),
        );
        assert_eq!(once.to_bits(), (1.0 + f64::from(f32::EPSILON)).to_bits());
        assert_eq!(
            wide.to_bits(),
            (1.0 + f64::from(f32::EPSILON) / 2.0).to_bits()
        );
        assert_eq!(f64::from(wide as f32).to_bits(), 1.0f64.to_bits());
    }

    #[test]
    fn a_grid_point_is_held_at_its_precision_and_a_binary64_one_is_not_held_at_binary32() {
        for p in BOTH {
            assert!(grid(p).iter().all(|&x| holds(p, x)), "{p:?}");
        }
        // A binary64 grid is not a binary32 one, so `sweep_at::<f32>` refuses it.
        let wide = grid(Precision::F64);
        assert!(wide.iter().any(|&x| !holds(Precision::F32, x)));
    }
}
