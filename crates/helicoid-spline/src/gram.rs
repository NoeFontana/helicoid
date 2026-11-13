//! Derivative Gram matrices of an Rⁿ spline (`NUMERICS.md` §21, `0115` decision 5).

use helicoid_linalg::{Mask, Matrix, Real};

use crate::order::{Order, Tables, ValidOrder};

/// `G^(R)`: the energy `∫ ‖p^(R)(t)‖² dt` over one segment is `cᵀ (G^(R) ⊗ I) c` in its `K`
/// controls. A weighted energy `Q⁻¹` is the caller's Kronecker product.
///
/// Each entry is `S::lit(num) / S::lit(den)`, one correctly rounded division, times `h⁻¹` raised to
/// `2R - 1` by repeated products, so a power-of-two `inv_h` scales exactly. `R` is 1 or 2, asserted
/// at monomorphization; `G^(2)` at `K = 2` is zero.
///
/// # Domain
///
/// `inv_h` is finite and positive, checked by `debug_assert!`.
pub fn derivative_gram<S: Real, const K: usize, const R: usize>(inv_h: S) -> Matrix<S, K, K>
where
    Order<K>: ValidOrder,
{
    const { assert!(R == 1 || R == 2, "derivative_gram: R is 1 or 2") };
    let () = Tables::<K>::EXACT;
    debug_assert!(
        S::zero().lt(inv_h).all(),
        "derivative_gram: inv_h is positive"
    );
    let (num, den) = &Tables::<K>::GRAM[R - 1];
    let mut scale = inv_h;
    for _ in 1..2 * R - 1 {
        scale = scale * inv_h;
    }
    let den = S::lit(*den as f64);
    let mut g = Matrix::<S, K, K>::identity().scale(S::zero());
    for (s, row) in num.iter().enumerate() {
        for (t, &n) in row.iter().enumerate() {
            g.set(s, t, S::lit(n as f64) / den * scale);
        }
    }
    g
}
