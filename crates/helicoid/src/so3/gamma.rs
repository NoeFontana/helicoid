//! The integrated exponentials `Γ_m(φ) = Σ Wⁿ/(n+m)!` of `NUMERICS.md` §7 and their directional
//! Jacobians (`docs/PHASE5.md` §4, `0064`; derivations in `docs/maths/gamma-gaussian.md` GG.1–GG.6).
//!
//! `Γ₀ = Exp`, `Γ₁ = J_l` and `Γ₂ = ½I + bW + dW²` (GG.2(a)). With body rate `ω` and specific force
//! `f` constant over `Δt`, `φ = ωΔt`: `ΔR = Exp φ`, `Δv = Δt Γ₁(φ) f`, `Δp = Δt² Γ₂(φ) f` (GG.4).
//! Preintegration itself is not `helicoid`'s (`0009`). `Γ` is entire, so nothing here has a
//! domain beyond the overflow of `θ²` that `Exp` has (GG.2, CO.15(c)).

use helicoid_linalg::{Dual, Mat3, Matrix, Real, Vec3, Vector};

use super::{norm_sq, HatSq, SO3Tangent, SO3};
use crate::coeffs::{gamma2_coeffs, jr_coeffs};
use crate::traits::LieGroup;

/// `Γ₁(φ) = Σ Wⁿ/(n+1)! = J_l(φ)` (GG.2(a)): [`SO3::jl`] itself, to the bit, so the velocity
/// increment and `SE_N(3)`'s `Exp` translation column (SE.3(b)) are one evaluation.
#[inline]
pub fn gamma1<S: Real>(phi: &SO3Tangent<S>) -> Mat3<S> {
    SO3::jl(phi)
}

/// `Γ₂(φ) = Σ Wⁿ/(n+2)! = ½I + bW + dW²` (`NUMERICS.md` §7, GG.2(a)), `(b, d)` from one grouped
/// `branch` on `θ²` with `b`'s and `d`'s swept switches (`0064`), assembled as [`SO3::jr`] is
/// (`0066`).
#[inline]
pub fn gamma2<S: Real>(phi: &SO3Tangent<S>) -> Mat3<S> {
    let w2 = HatSq::of(phi.phi);
    let (b, d) = gamma2_coeffs(norm_sq(phi.phi));
    w2.poly(S::lit(0.5), b, d)
}

/// `(Γ_M(φ) v, ∂(Γ_M(φ) v)/∂φ, ∂(Γ_M(φ) v)/∂v = Γ_M(φ))` for `M` 1 or 2 (`PHASE5.md` §4, `0064`,
/// `0066`): the Jacobians in argument order, as every `*_jacobians` returns them. A preintegration
/// step needs all three (`Δp = Δt² Γ₂ f`: the increment, its gyro and its accelerometer
/// Jacobian), and they share one coefficient evaluation.
///
/// GG.2(d)'s vector form `g = v/M! + σ₁ φ×v + σ₂ φ×(φ×v)`, with `(σ₁, σ₂)` `jr_coeffs`' `(a, b)`
/// (`M = 1`) or `gamma2_coeffs`' `(b, d)` (`M = 2`), evaluated once on `Dual<S, 1>` seeded at
/// `z = θ²`, and the chain rule
///
/// `∂g/∂φ = −σ₁ v^ + σ₂ ((φ·v) I + φvᵀ − 2vφᵀ) + 2 (σ₁′ φ×v + σ₂′ φ×(φ×v)) φᵀ`.
///
/// The value and the third output are the plain evaluation's bits (`0003`): the third is
/// [`gamma1`] or [`gamma2`] to the bit. The value is the vector form, not `Γ_M(φ) · v`: a different
/// rounding of the same number.
///
/// The Jacobian is the derivative of the computed coefficients (GG.6(a)): finite at `φ = 0`,
/// subnormal included, where it is `−v^/(M+1)!` (GG.5(b)); worse than the value by about `2K/θ`
/// on a `K`-term series arm (GG.6(c)); and on the exact arm, which runs from `θ = 2.89` (`M = 1`)
/// or `3.11` (`M = 2`) at `f64` and from `2.79` at `f32`, about `u/θ` (`M = 1`) or `u/θ²` (`M = 2`)
/// (GG.6(d)).
///
/// `M` is const and checked at monomorphization; `M = 3` has coefficients but loses `2u/θ²` on its
/// exact arm (GG.3(c)) and `M ≥ 4` needs a coefficient outside `NUMERICS.md` §4 (GG.2), so neither
/// is offered.
///
/// ```
/// use helicoid::{so3, SO3Tangent};
/// use helicoid_linalg::Vector;
///
/// let phi = SO3Tangent { phi: Vector([0.0, 0.0, 0.0_f64]) };
/// let (gv, j, g) = so3::gamma_apply_jacobians::<2, _>(&phi, Vector([1.0, 2.0, 3.0]));
/// assert_eq!(gv.0, [0.5, 1.0, 1.5]);
/// assert_eq!(j.get(0, 1), 3.0 / 6.0); // −v^/3! at φ = 0
/// assert_eq!(g.get(0, 0), 0.5); // Γ₂(0) = ½I
/// ```
///
/// `M = 3` does not build:
///
/// ```compile_fail,E0080
/// use helicoid::{so3, SO3Tangent};
/// use helicoid_linalg::Vector;
///
/// let phi = SO3Tangent { phi: Vector([0.1, 0.0, 0.0_f64]) };
/// let _ = so3::gamma_apply_jacobians::<3, _>(&phi, Vector([1.0, 0.0, 0.0]));
/// ```
#[inline]
pub fn gamma_apply_jacobians<const M: usize, S: Real>(
    phi: &SO3Tangent<S>,
    v: Vec3<S>,
) -> (Vec3<S>, Mat3<S>, Mat3<S>) {
    const { assert!(M == 1 || M == 2, "gamma_apply_jacobians: M is 1 or 2") };
    let p = phi.phi;
    // Everything that does not read the coefficients goes before them (`HatSq`'s reason).
    let w2 = HatSq::of(p);
    let pv = p.cross(v);
    let ppv = p.cross(pv);
    let pdotv = p.dot(v);
    let z = Dual::<S, 1>::variable(norm_sq(p), 0);
    // `M` is a const: one arm survives monomorphization.
    let (s0, s1, s2) = if M == 1 {
        let (a, b) = jr_coeffs(z);
        (S::one(), a, b)
    } else {
        let (b, d) = gamma2_coeffs(z);
        (S::lit(0.5), b, d)
    };
    let (ds1, ds2, s1, s2) = (s1.d[0], s2.d[0], s1.v, s2.v);
    let value = (v.scale(s0) + pv.scale(s1)) + ppv.scale(s2);
    // `J = s2 (φ·v) I − s1 v^ + ((s2 φ) vᵀ + k φᵀ) − (2 s2 v) φᵀ`, `k = 2 (σ₁′ φ×v + σ₂′ φ×(φ×v))`:
    // the two `φᵀ` terms kept apart, the grouping `0066` measured best near `π`.
    let diag = s2 * pdotv;
    let a = p.scale(s2);
    let k = (pv.scale(ds1) + ppv.scale(ds2)).scale(S::lit(2.0));
    let m2 = v.scale(S::lit(2.0) * s2);
    let [x, y, w] = v.0.map(|e| s1 * e);
    let skew = [[S::zero(), w, -y], [-w, S::zero(), x], [y, -x, S::zero()]];
    let jac = Matrix::from_rows(core::array::from_fn(|i| {
        Vector(core::array::from_fn(|j| {
            let outer = (a.0[i] * v.0[j] + k.0[i] * p.0[j]) - m2.0[i] * p.0[j];
            if i == j {
                diag + outer
            } else {
                skew[i][j] + outer
            }
        }))
    }));
    (value, jac, w2.poly(s0, s1, s2))
}
