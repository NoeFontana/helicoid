//! The switch grid of `docs/PHASE1.md` §6: 64 points per decade of the branch variable `z = θ²`,
//! from `z = 1e-16` (`θ = 1e-8`) to `z = 1` (`θ = 1`), both ends included.
//!
//! Point `i` is `10^((i - 1024)/64)` correctly rounded to binary64, from the 64th root of a power
//! of ten in integers: the grid does not depend on the quality of any `pow`.

use num_bigint::BigUint;

use crate::conformance::number::ratio_to_f64;

/// Grid points per decade of `z`.
const PER_DECADE: usize = 64;
/// Decades of `z` spanned: `θ` from `1e-8` to `1`.
const DECADES: usize = 16;

/// `10^(-n/64)` rounded to nearest. `m = ⌊2^P · 10^(-n/64)⌋` is the integer 64th root of
/// `2^(64 P) / 10^n` (the floor of a root is the root of the floor); `(2m + 1) / 2^(P + 1)` lies
/// strictly between `m` and `m + 1` in units of `2^-P`, so the one rounding of `ratio_to_f64` is
/// never a tie: it is the nearest binary64 to the root, irrational or a power of ten (no dyadic).
fn pow10_neg(n: u32) -> f64 {
    const P: usize = 128;
    let scaled = (BigUint::from(1u32) << (64 * P)) / BigUint::from(10u32).pow(n);
    let m = scaled.nth_root(64);
    ratio_to_f64(&((m << 1usize) + 1u32), &(BigUint::from(1u32) << (P + 1)))
}

/// The switch points, increasing: `1025` values from `1e-16` to `1`.
pub(super) fn grid() -> Vec<f64> {
    let top = DECADES * PER_DECADE;
    (0..=top).map(|i| pow10_neg((top - i) as u32)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::number::Dyadic;

    #[test]
    fn the_grid_has_64_points_per_decade_and_both_ends() -> Result<(), String> {
        let g = grid();
        assert_eq!(g.len(), 16 * 64 + 1);
        assert!(g.windows(2).all(|w| w[0] < w[1]));
        // Each decade starts at its power of ten, which Rust's parser rounds independently.
        for j in 0..=16 {
            let ten: f64 = format!("1e-{j}").parse().map_err(|e| format!("{e}"))?;
            assert_eq!(g[(16 - j) * 64].to_bits(), ten.to_bits(), "1e-{j}");
        }
        assert_eq!(
            (g[0].to_bits(), g[1024].to_bits()),
            (1e-16f64.to_bits(), 1.0f64.to_bits())
        );
        Ok(())
    }

    /// `10^(-n/64)` lies between the midpoints from `x` to its neighbours: `mid^64 · 10^n` against
    /// 1, in integers over `2^1075` (`x` normal and at least `1e-16`).
    fn is_nearest(x: f64, n: u32) -> bool {
        let scaled = |y: f64| {
            let d = Dyadic::of(y);
            BigUint::from(d.mant) << (d.exp + 1074) as usize
        };
        let power = |mid: BigUint| mid.pow(64) * BigUint::from(10u32).pow(n);
        let one = BigUint::from(1u32) << (64 * 1075);
        let (lo, here, hi) = (scaled(x.next_down()), scaled(x), scaled(x.next_up()));
        power(&lo + &here) <= one && one <= power(&here + &hi)
    }

    #[test]
    fn every_point_is_the_nearest_binary64_to_its_power_of_ten() {
        let g = grid();
        for (i, &x) in g.iter().enumerate() {
            assert!(is_nearest(x, (1024 - i) as u32), "point {i}: {x:e}");
        }
    }
}
