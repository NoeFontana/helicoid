//! The integrated exponentials `Γ_m(φ) = Σ Wⁿ/(n+m)!` of `NUMERICS.md` §7 and their directional
//! Jacobians (`docs/PHASE5.md` §4, `0064`; derivations in `docs/maths/gamma-gaussian.md` GG.1–GG.6).
//!
//! `Γ₀ = Exp`, `Γ₁ = J_l` and `Γ₂ = ½I + bW + dW²` (GG.2(a)). With body rate `ω` and specific force
//! `f` constant over `Δt`, `φ = ωΔt`: `ΔR = Exp φ`, `Δv = Δt Γ₁(φ) f`, `Δp = Δt² Γ₂(φ) f` (GG.4).
//! Preintegration itself is not `helicoid`'s (`0009`). `Γ` is entire, so nothing here has a
//! domain beyond the overflow of `θ²` that `Exp` has (GG.2, CO.15(c)).

use helicoid_linalg::{Dual, Mat3, Matrix, Real, Vec3, Vector};

use super::{exp_from, norm_sq, HatSq, SO3Tangent, SO3};
use crate::coeffs::{exp_gamma_coeffs, exp_gamma_coeffs_dz, gamma2_coeffs, jr_coeffs_sinc};
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

/// `(Γ_M(φ) v, ∂(Γ_M(φ) v)/∂φ, ∂(Γ_M(φ) v)/∂v = Γ_M(φ))`: what [`gamma_apply_jacobians`] returns,
/// and [`exp_gamma_apply_jacobians`] for each `M`.
pub type GammaJacobians<S> = (Vec3<S>, Mat3<S>, Mat3<S>);

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
/// on a `K`-term series arm (GG.6(c)); and on the exact arm, which runs from `θ = 2.41` (`M = 1`)
/// or `2.79` (`M = 2`) at `f64` and from `3.00` at `f32`, about `u/θ` (`M = 1`) or `u/θ²` (`M = 2`)
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
) -> GammaJacobians<S> {
    const { assert!(M == 1 || M == 2, "gamma_apply_jacobians: M is 1 or 2") };
    let o = Operands::of(phi.phi, v);
    let z = Dual::<S, 1>::variable(norm_sq(phi.phi), 0);
    // `M` is a const: one arm survives monomorphization.
    // `M = 1` takes `SO3::jl`'s diagonal, so the third output is `gamma1` to the bit
    // (`0089`); `M = 2` has no such rewrite.
    let (s0, s1, s2, far) = if M == 1 {
        let (a, b, sinc, far) = jr_coeffs_sinc(z);
        (S::one(), a, b, Some((sinc.v, far)))
    } else {
        let (b, d) = gamma2_coeffs(z);
        (S::lit(0.5), b, d, None)
    };
    apply_jacobians_from(&o, s0, s1, s2, far)
}

/// What `Γ_M v` and its Jacobians read besides the coefficients: `φ`, `v`, `W²`, `φ×v`, `φ×(φ×v)`
/// and `φ·v`. Formed before the coefficients (`HatSq`'s reason), and once for both `M` in
/// [`exp_gamma_apply_jacobians`].
#[derive(Clone, Copy)]
struct Operands<S: Real> {
    p: Vec3<S>,
    v: Vec3<S>,
    w2: HatSq<S>,
    pv: Vec3<S>,
    ppv: Vec3<S>,
    pdotv: S,
}

impl<S: Real> Operands<S> {
    #[inline(always)]
    fn of(p: Vec3<S>, v: Vec3<S>) -> Self {
        let w2 = HatSq::of(p);
        let pv = p.cross(v);
        let ppv = p.cross(pv);
        let pdotv = p.dot(v);
        Self {
            p,
            v,
            w2,
            pv,
            ppv,
            pdotv,
        }
    }
}

/// [`gamma_apply_jacobians`] after its coefficients, `σ₁` and `σ₂` given on `Dual`'s lane: the one
/// body [`exp_gamma_apply_jacobians`] runs for each `M`, so each of its triples is that function's
/// to the bit (`0119`). `far` is `Γ₁`'s `0089` diagonal. Always inlined: out of line, its triple
/// returns through memory and `gamma_apply_jacobians` took 105–155 instructions more per call.
#[inline(always)]
fn apply_jacobians_from<S: Real>(
    o: &Operands<S>,
    s0: S,
    s1: Dual<S, 1>,
    s2: Dual<S, 1>,
    far: Option<(S, S::Mask)>,
) -> GammaJacobians<S> {
    let Operands {
        p,
        v,
        w2,
        pv,
        ppv,
        pdotv,
    } = *o;
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
    let g = match far {
        Some((sinc, far)) => w2.poly_far(s0, s1, s2, sinc, far),
        None => w2.poly(s0, s1, s2),
    };
    (value, jac, g)
}

/// `(Exp φ, Γ₁(φ) v, Γ₂(φ) v)` from one coefficient evaluation (`0119`): [`SO3::exp`],
/// `gamma1(φ) * v` and `gamma2(φ) * v`, each **to the bit**, for a step that integrates a constant
/// rate and force, `ΔR = Exp φ`, `Δv = Δt Γ₁ f`, `Δp = Δt² Γ₂ f` (GG.4).
///
/// One `θ²`, one `W²` and one coefficient evaluation serve all three, each coefficient on the arm
/// its own call takes. `Γ_M v` is `Mat3 * Vec3`, not [`gamma_apply_jacobians`]' vector form: a
/// different rounding of the same number.
///
/// ```
/// use helicoid::{so3, LieGroup, SO3Tangent, SO3};
/// use helicoid_linalg::Vector;
///
/// let phi = SO3Tangent { phi: Vector([0.3, -0.2, 0.1_f64]) };
/// let v = Vector([1.0, 2.0, 3.0]);
/// let (r, g1v, g2v) = so3::exp_gamma_apply(&phi, v);
/// let (q, e) = (r.quat(), SO3::exp(&phi).quat());
/// assert_eq!([q.w, q.x, q.y, q.z], [e.w, e.x, e.y, e.z]);
/// assert_eq!(g1v.0, (so3::gamma1(&phi) * v).0);
/// assert_eq!(g2v.0, (so3::gamma2(&phi) * v).0);
/// ```
///
/// # Domain
///
/// [`SO3::exp`]'s: every `φ` whose `θ²` does not overflow (GG.2, CO.15(c)).
#[inline]
pub fn exp_gamma_apply<S: Real>(phi: &SO3Tangent<S>, v: Vec3<S>) -> (SO3<S>, Vec3<S>, Vec3<S>) {
    let w2 = HatSq::of(phi.phi);
    let c = exp_gamma_coeffs(norm_sq(phi.phi));
    let r = exp_from(phi.phi, c.k, c.cos_half);
    // `Γ₁ = J_l = J_r(−φ)`: `HatSq::of(−φ)` is `HatSq::of(φ)` and `(−a)(−x)` is `ax`, exactly
    // (`gamma1_is_jl_to_the_bit`).
    let g = c.gamma;
    let g1 = w2.poly_far(S::one(), g.a, g.b1, g.sinc, g.far);
    let g2 = w2.poly(S::lit(0.5), g.b2, g.d);
    (r, g1 * v, g2 * v)
}

/// [`exp_gamma_apply`]'s rotation and [`gamma_apply_jacobians`]' triple for `M = 1` and for `M = 2`,
/// `(Exp φ, (Γ₁ v, ∂(Γ₁ v)/∂φ, Γ₁), (Γ₂ v, ∂(Γ₂ v)/∂φ, Γ₂))`, each **to the bit**, from one
/// coefficient evaluation on `Dual<S, 1>` (`0119`): a preintegration step's increment and both of
/// its Jacobians.
///
/// `b` is evaluated once for both triples except between the groups' second switches, where each
/// group's own switch picks its arm; `W²`, `φ×v`, `φ×(φ×v)` and `φ·v` are formed once. The
/// rotation reads the value lanes of `k` and `cos θ/2`, which are `exp_coeffs`' plain bits (`0003`
/// item 5). `Γ₁`'s third output keeps `SO3::jl`'s diagonal (`0089`). The values are the vector
/// form, as in [`gamma_apply_jacobians`], not [`exp_gamma_apply`]'s `Γ_M * v`.
///
/// ```
/// use helicoid::{so3, LieGroup, SO3Tangent, SO3};
/// use helicoid_linalg::Vector;
///
/// let phi = SO3Tangent { phi: Vector([0.3, -0.2, 0.1_f64]) };
/// let v = Vector([1.0, 2.0, 3.0]);
/// let (r, g1, g2) = so3::exp_gamma_apply_jacobians(&phi, v);
/// let (q, e) = (r.quat(), SO3::exp(&phi).quat());
/// assert_eq!([q.w, q.x, q.y, q.z], [e.w, e.x, e.y, e.z]);
/// let (w1, w2) = (
///     so3::gamma_apply_jacobians::<1, _>(&phi, v),
///     so3::gamma_apply_jacobians::<2, _>(&phi, v),
/// );
/// assert_eq!((g1.0 .0, g1.1.get(0, 1), g1.2.get(2, 2)), (w1.0 .0, w1.1.get(0, 1), w1.2.get(2, 2)));
/// assert_eq!((g2.0 .0, g2.1.get(0, 1), g2.2.get(2, 2)), (w2.0 .0, w2.1.get(0, 1), w2.2.get(2, 2)));
/// ```
///
/// # Domain
///
/// [`SO3::exp`]'s.
#[inline]
pub fn exp_gamma_apply_jacobians<S: Real>(
    phi: &SO3Tangent<S>,
    v: Vec3<S>,
) -> (SO3<S>, GammaJacobians<S>, GammaJacobians<S>) {
    let o = Operands::of(phi.phi, v);
    let c = exp_gamma_coeffs_dz(norm_sq(phi.phi));
    let r = exp_from(phi.phi, c.k, c.cos_half);
    let g = c.gamma;
    let g1 = apply_jacobians_from(&o, S::one(), g.a, g.b1, Some((g.sinc.v, g.far)));
    let g2 = apply_jacobians_from(&o, S::lit(0.5), g.b2, g.d, None);
    (r, g1, g2)
}
