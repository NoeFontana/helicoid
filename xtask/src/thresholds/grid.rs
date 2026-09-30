//! The switch grid of `docs/PHASE1.md` §6: 64 points per decade of the branch variable `z = θ²`
//! (`s = n²/w²` for `r`), from `z = 1e-16` (`θ = 1e-8`) to `z = 1` (`θ = 1`), both ends included.
//!
//! Point `i` is `10^((i - 1024)/64)` correctly rounded to the precision it is swept at (binary32
//! points are held as the binary64 they equal), from the 64th root of a power of ten in integers:
//! the grid does not depend on the quality of any `pow`, and a binary32 point is never a binary64
//! one rounded again (`docs/decisions/0016` item 3).

use helicoid_linalg::Precision;
use num_bigint::BigUint;

use crate::conformance::number::{ratio_to_f32, ratio_to_f64};

/// Grid points per decade of `z`.
const PER_DECADE: usize = 64;
/// Decades of `z` spanned: `θ` from `1e-8` to `1`.
const DECADES: usize = 16;

/// `10^(-n/64)` rounded to nearest at `precision`. `m = ⌊2^P · 10^(-n/64)⌋` is the integer 64th
/// root of `2^(64 P) / 10^n` (the floor of a root is the root of the floor); `(2m + 1) / 2^(P + 1)`
/// lies strictly between `m` and `m + 1` in units of `2^-P`, so the one rounding is never a tie: it
/// is the nearest value to the root, irrational or a power of ten (no dyadic).
fn pow10_neg(n: u32, precision: Precision) -> f64 {
    const P: usize = 128;
    let scaled = (BigUint::from(1u32) << (64 * P)) / BigUint::from(10u32).pow(n);
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

/// The switch points, increasing: `1025` values from `1e-16` to `1`.
pub(super) fn grid(precision: Precision) -> Vec<f64> {
    let top = DECADES * PER_DECADE;
    (0..=top)
        .map(|i| pow10_neg((top - i) as u32, precision))
        .collect()
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
            assert_eq!(g.len(), 16 * 64 + 1);
            assert!(g.windows(2).all(|w| w[0] < w[1]));
            // Each decade starts at its power of ten, which Rust's parser rounds independently; the
            // ends, `1e-16` and `1`, are two of them.
            for j in 0..=16 {
                let text = format!("1e-{j}");
                let ten = match p {
                    Precision::F64 => text.parse::<f64>().map_err(|e| format!("{e}"))?,
                    Precision::F32 => f64::from(text.parse::<f32>().map_err(|e| format!("{e}"))?),
                };
                assert_eq!(g[(16 - j) * 64].to_bits(), ten.to_bits(), "{p:?} 1e-{j}");
            }
        }
        Ok(())
    }

    /// `10^(-n/64)` lies between the midpoints from `x` to its neighbours `lo` and `hi`: `mid^64 ·
    /// 10^n` against 1, in integers over `2^1075` (`x` normal and at least `1e-16`).
    fn is_nearest(x: f64, (lo, hi): (f64, f64), n: u32) -> bool {
        let scaled = |y: f64| {
            let d = Dyadic::of(y);
            BigUint::from(d.mant) << (d.exp + 1074) as usize
        };
        let power = |mid: BigUint| mid.pow(64) * BigUint::from(10u32).pow(n);
        let one = BigUint::from(1u32) << (64 * 1075);
        let here = scaled(x);
        power(&scaled(lo) + &here) <= one && one <= power(&here + &scaled(hi))
    }

    #[test]
    fn every_point_is_the_nearest_value_to_its_power_of_ten() {
        let f32s = |x: f64| {
            let y = x as f32;
            (f64::from(y.next_down()), f64::from(y.next_up()))
        };
        for (i, &x) in grid(Precision::F64).iter().enumerate() {
            let near = (x.next_down(), x.next_up());
            assert!(is_nearest(x, near, (1024 - i) as u32), "f64 {i}: {x:e}");
        }
        // Binary32 points are on binary32's own spacing: not a binary64 point cast down.
        for (i, &x) in grid(Precision::F32).iter().enumerate() {
            assert!(holds(Precision::F32, x), "point {i} is not a binary32");
            assert!(is_nearest(x, f32s(x), (1024 - i) as u32), "f32 {i}: {x:e}");
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
