//! The integrated exponentials `Γ_m(φ) = Σ Wⁿ/(n+m)!` of `NUMERICS.md` §7 and their directional
//! Jacobians (`docs/PHASE5.md` §4, `0064`; derivations in `docs/maths/gamma-gaussian.md` GG.1–GG.6).
//!
//! `Γ₀ = Exp`, `Γ₁ = J_l` and `Γ₂ = ½I + bW + dW²` (GG.2(a)). With body rate `ω` and specific force
//! `f` constant over `Δt`, `φ = ωΔt`: `ΔR = Exp φ`, `Δv = Δt Γ₁(φ) f`, `Δp = Δt² Γ₂(φ) f` (GG.4).
//! Preintegration itself is not `helicoid`'s (`0009`). `Γ` is entire, so nothing here has a
//! domain beyond the overflow of `θ²` that `Exp` has (GG.2, CO.15(c)).

use core::array;

use helicoid_linalg::{hat, Dual, Mat3, Matrix, Real, Vec3, Vector};

use super::{mul_hat, norm_sq, SO3Tangent, SO3};
use crate::coeffs::{gamma2_coeffs, jr_coeffs};
use crate::traits::LieGroup;

/// `Γ₁(φ) = Σ Wⁿ/(n+1)! = J_l(φ)` (GG.2(a)): [`SO3::jl`] itself, to the bit, so the velocity
/// increment and `SE_N(3)`'s `Exp` translation column (SE.3(b)) are one evaluation.
#[inline]
pub fn gamma1<S: Real>(phi: &SO3Tangent<S>) -> Mat3<S> {
    SO3::jl(phi)
}

/// `Γ₂(φ) = Σ Wⁿ/(n+2)! = ½I + bW + dW²` (`NUMERICS.md` §7, GG.2(a)), `(b, d)` from one grouped
/// `branch` on `θ²` with `b`'s and `d`'s swept switches (`0064`).
///
/// Assembled as `SO3::jr` is, `W²` through the structured product and not `φφᵀ − θ²I`, so the
/// three `Γ` share one evaluation order.
#[inline]
pub fn gamma2<S: Real>(phi: &SO3Tangent<S>) -> Mat3<S> {
    let w = hat(phi.phi);
    let (b, d) = gamma2_coeffs(norm_sq(phi.phi));
    (Matrix::identity().scale(S::lit(0.5)) + w.scale(b)) + mul_hat(&w, phi.phi).scale(d)
}

/// `(Γ_M(φ) v, ∂(Γ_M(φ) v)/∂φ)` for `M` 1 or 2 (`PHASE5.md` §4, `0064`): GG.2(d)'s
/// `v/M! + σ φ×v + σ' φ×(φ×v)` evaluated once on `Dual<S, 3>` seeded at `φ`, so the value is the
/// real part and the Jacobian, row `i` column `j` `∂(Γ_M v)_i/∂φ_j`, the lanes.
///
/// `M` is const and checked at monomorphization; `M = 3` has coefficients but loses `2u/θ²` on its
/// exact arm (GG.3(c)) and `M ≥ 4` needs a coefficient outside `NUMERICS.md` §4 (GG.2), so neither
/// is offered.
///
/// The Jacobian is the derivative of the computed function (GG.6(a)): finite at `φ = 0`, subnormal
/// included, where it is `−v^/(M+1)!` (GG.5(b)); on the series arm worse than the value by about
/// `2K/θ` for `K` terms (GG.6(c)); on the exact arm about `u/θ` (`M = 1`) or `u/θ²` (`M = 2`)
/// (GG.6(d)). The value is the vector form, not `gamma2(φ) · v`: a different rounding of the same
/// number.
///
/// ```
/// use helicoid::{so3, SO3Tangent};
/// use helicoid_linalg::Vector;
///
/// let phi = SO3Tangent { phi: Vector([0.0, 0.0, 0.0_f64]) };
/// let (gv, j) = so3::gamma_apply_jacobian::<2, _>(&phi, Vector([1.0, 2.0, 3.0]));
/// assert_eq!(gv.0, [0.5, 1.0, 1.5]);
/// assert_eq!(j.get(0, 1), 3.0 / 6.0); // −v^/3! at φ = 0
/// ```
///
/// `M = 3` does not build:
///
/// ```compile_fail,E0080
/// use helicoid::{so3, SO3Tangent};
/// use helicoid_linalg::Vector;
///
/// let phi = SO3Tangent { phi: Vector([0.1, 0.0, 0.0_f64]) };
/// let _ = so3::gamma_apply_jacobian::<3, _>(&phi, Vector([1.0, 0.0, 0.0]));
/// ```
pub fn gamma_apply_jacobian<const M: usize, S: Real>(
    phi: &SO3Tangent<S>,
    v: Vec3<S>,
) -> (Vec3<S>, Mat3<S>) {
    const { assert!(M == 1 || M == 2, "gamma_apply_jacobian: M is 1 or 2") };
    let p: Vec3<Dual<S, 3>> = Vector(array::from_fn(|i| Dual::variable(phi.phi.0[i], i)));
    let v: Vec3<Dual<S, 3>> = Vector(v.0.map(Dual::constant));
    let z = norm_sq(p);
    // `M` is a const: one arm survives monomorphization.
    let (head, s1, s2) = if M == 1 {
        let (a, b) = jr_coeffs(z);
        (v, a, b)
    } else {
        let (b, d) = gamma2_coeffs(z);
        (v.scale(Dual::lit(0.5)), b, d)
    };
    let pv = p.cross(v);
    let g = (head + pv.scale(s1)) + p.cross(pv).scale(s2);
    (
        Vector(g.0.map(|x| x.v)),
        Matrix::from_rows(g.0.map(|x| Vector(x.d))),
    )
}
