//! The cumulative weights `λ`, `λ̇`, `λ̈` at one segment parameter (`NUMERICS.md` §21).

use helicoid_linalg::{Mask, Real};

use crate::order::{Order, Tables, ValidOrder};

/// The cumulative weights of a segment at `u`, up to derivative order `D`.
///
/// `D` (0, 1 or 2) is the highest derivative formed. Rows above it read 0 and cost nothing, and
/// every evaluation takes its `D` from its weights, so a value-only caller forms `λ` alone
/// (`0114` decision 5). Slot 0 is `λ_0 = 1`, `λ̇_0 = λ̈_0 = 0`, which the cumulative form never
/// multiplies by anything.
#[derive(Clone, Copy, Debug)]
pub struct Weights<S, const K: usize, const D: usize> {
    /// `λ_j(u)`.
    pub lambda: [S; K],
    /// `λ̇_j = h⁻¹ dλ_j/du`; zero when `D < 1`.
    pub dot: [S; K],
    /// `λ̈_j = h⁻² d²λ_j/du²`; zero when `D < 2`.
    pub ddot: [S; K],
}

impl<S: Real, const K: usize, const D: usize> Weights<S, K, D>
where
    Order<K>: ValidOrder,
{
    /// The weights at `u`, for a knot spacing `h = 1 / inv_h`.
    ///
    /// Each row is Horner on the exact numerators of `0115`, then one division by `(K-1)!`, then
    /// the `h⁻¹` product: the order `NUMERICS.md` §21 states.
    ///
    /// # Domain
    ///
    /// `u` is finite; outside `[0, 1]` the polynomial extrapolates (R6). `inv_h` is finite and
    /// positive, checked by `debug_assert!`; it is not read when `D = 0`.
    #[inline]
    pub fn at(u: S, inv_h: S) -> Self {
        const { assert!(D <= 2, "Weights: D is 0, 1 or 2") };
        let () = Tables::<K>::EXACT;
        let table = &Tables::<K>::CUMULATIVE;
        let f = S::lit(Tables::<K>::FACTORIAL as f64);
        let mut lambda = [S::zero(); K];
        let mut dot = [S::zero(); K];
        let mut ddot = [S::zero(); K];
        lambda[0] = S::one();
        for j in 1..K {
            lambda[j] = horner(u, &table[j], 0) / f;
        }
        if D >= 1 {
            debug_assert!(S::zero().lt(inv_h).all(), "Weights::at: inv_h is positive");
            for j in 1..K {
                dot[j] = inv_h * (horner(u, &table[j], 1) / f);
            }
        }
        if D >= 2 {
            let inv_h2 = inv_h * inv_h;
            for j in 1..K {
                ddot[j] = inv_h2 * (horner(u, &table[j], 2) / f);
            }
        }
        Self { lambda, dot, ddot }
    }
}

/// `Σ_{k ≥ r} k!/(k-r)! · c_k · u^{k-r}` by Horner: the `r`-th derivative of the row's
/// polynomial. Each coefficient is an integer below `2^24` (`Tables::EXACT`), so `S::lit` is exact.
#[inline]
fn horner<S: Real, const K: usize>(u: S, c: &[i64; K], r: usize) -> S {
    let coeff = |k: usize| {
        let mut f = 1;
        for i in 0..r {
            f *= (k - i) as i64;
        }
        S::lit((f * c[k]) as f64)
    };
    let mut acc = coeff(K - 1);
    for k in (r..K - 1).rev() {
        acc = acc * u + coeff(k);
    }
    acc
}
