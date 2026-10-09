# helicoid — Phase 5 Implementation Specification: Extended Geometry and Retraction Migrations

> **Companions:** `docs/NUMERICS.md` §7–§9, [`0012`](./decisions/0012-a-retraction-is-a-chart.md),
> [`0005`](./decisions/0005-the-jacobian-is-a-dual-matrix.md).

**Deliverable:** charts as types; S²; Sim(3); Γ₁/Γ₂ with directional Jacobians; `Gaussian` with a
stated side; and consumer migrations onto `helicoid` without unrecorded regressions. Sections marked
**NORMATIVE** are requirements.

## 0.0 Implementation status

| Area | Status |
|---|---|
| `Chart`, `Manifold`, `WithChart`; group charts (§1) | Done ([`0060`](./decisions/0060-the-charts-go-first-and-name-their-frame.md), #119): `RightChart`, `LeftChart` on every group, one `Manifold` impl per group, `Lifted`; the chart laws, `Dual` Jacobians and CH.6 on SO(3), SE(3), SE₂(3), SE₃(3), ℝ³, SO(3)×ℝ³ |
| SE(3) charts `Screw`, `Decoupled`, `WorldTranslation` (§1.3) | Done (`0060`, #120, #123): `TwistBlockJac`, `Se3Chart`, `SE3::chart_transition` (all nine pairs against `Dual`); corpus ids `se3_{screw,decoupled,world}_{retract,local}` |
| S² and its chart (§2) | Not started |
| Sim(3) (§3) | Not started |
| Γ₁, Γ₂, directional Jacobians (§4) | Done ([`0064`](./decisions/0064-the-integrated-exponentials-reuse-the-swept-switches.md), #127, #128): `so3::{gamma1, gamma2, gamma_apply_jacobians::<M>}`, `M` 1 or 2 (one coefficient evaluation, [`0066`](./decisions/0066-gamma-shares-its-coefficients-and-a-gaussian-factors-once.md)); the twin `reference::gamma_apply`; corpus id `so3_gamma2` at both precisions; the seeded defect Γ₂ from Γ₁'s coefficients |
| `Gaussian<S, G, Sd, D>` (§5) | Done ([`0065`](./decisions/0065-a-gaussian-names-its-side-and-is-stored-symmetric.md), #129, #130): `to_left`, `to_right`, `propagate`, `mahalanobis_sq`, and `whitener` for many-point gating ([`0066`](./decisions/0066-gamma-shares-its-coefficients-and-a-gaussian-factors-once.md)); the round-trip bound and the twins on SE(3), SE₂(3), SO(3), SO(3)×ℝ³; corpus ids `gaussian_mahalanobis_{se3,se23}` |
| Consumer migrations (§6) | Not started |

## 0. Non-goals and guardrails — read first

**NORMATIVE.** No preintegration (Γ is the primitive; integrators live downstream), no
sampling or RNG in `Gaussian`, no solver-facing ambient Jacobians (Phase 6). Sim(3) implementation
awaits completion of `NUMERICS.md` §9. §1, §2, §4 and §5 do not wait for Phase 4 or for
`PHASE3.md` §6, and §1 comes first ([`0060`](./decisions/0060-the-charts-go-first-and-name-their-frame.md)).

## 1. Charts

**NORMATIVE** ([`0012`](./decisions/0012-a-retraction-is-a-chart.md)).

### 1.1 Traits

```rust
pub trait Chart<S: Real, M>: Copy {
    type Tangent: Tangent<S>;
    type Jac: Jac<S, Self::Tangent>;
    fn at(base: &M) -> Self;
    fn base(&self) -> M;
    fn retract(&self, delta: &Self::Tangent) -> M;
    fn local(&self, other: &M) -> Self::Tangent;
    fn retract_jacobian(&self, delta: &Self::Tangent) -> Self::Jac;
    fn local_jacobian(&self, other: &M) -> Self::Jac;
    /// Provided: `(local(other), local_jacobian(other))`, overridden to share the `Log` (`0060`).
    fn local_with_jacobian(&self, other: &M) -> (Self::Tangent, Self::Jac);
}

pub trait Manifold<S: Real>: Copy + Blend<S> {
    type Tangent: Tangent<S>;
    type Chart: Chart<S, Self, Tangent = Self::Tangent>;
    const DOF: usize;
}

#[repr(transparent)]
pub struct WithChart<M, C>(pub M, core::marker::PhantomData<C>);
```

Every `LieGroup` is a `Manifold` with `Chart = RightChart<G>` (retract $= X\,\mathrm{Exp}(\delta)$,
local $= Y \ominus_R X$). `LeftChart<G>` is the left twin. Each group has its own `Manifold` impl: a
blanket over `LieGroup` beside `WithChart`'s is E0119. A chart holds its frozen base, so it carries
`S`. `WithChart<M, C>`'s chart is `Lifted<C>`, which delegates to `C`
([`0060`](./decisions/0060-the-charts-go-first-and-name-their-frame.md)).

### 1.2 Rationale

Linearization uses `M::Chart` by default. Alternative retractions are explicit in types:
`WithChart<SE3<f64>, Decoupled<f64>>`.

### 1.3 SE(3) charts

| Chart | `retract(δ = [φ; ρ])` | `local(Y)` | Target |
|---|---|---|---|
| `Screw` (`RightChart<SE3>`) | $X\,\mathrm{Exp}(\delta)$ | $\mathrm{Log}(X^{-1}Y)$ | Default Lie exponential |
| `Decoupled` | $(R\,\mathrm{Exp}(\varphi),\ t + R\rho)$ | $(\mathrm{Log}(R^\top R_Y),\ R^\top(t_Y - t))$ | Decoupled rotation/translation |
| `WorldTranslation` | $(R\,\mathrm{Exp}(\varphi),\ t + \rho)$ | $(\mathrm{Log}(R^\top R_Y),\ t_Y - t)$ | `Product<SO3, R3>` chart |

`Screw` and `Decoupled` agree to first order at $\delta = 0$. `WorldTranslation` equals `Decoupled`
after $\mathrm{diag}(I, R)$ on the translation tangent, exactly, and agrees with the other two to
first order only if $R = I$ (`docs/maths/charts.md` CH.4(b), CH.5(c); `0060`).
`SE3::chart_transition::<From, To>()` returns that first-order transition. All three take the tangent
`Twist<S>`. Corpus ids: `se3_{screw,decoupled,world}_{retract,local}` over the `theta:*` and
`rho:*` strata, with a generic base (`0060` decision 8).

### 1.4 Jacobian types

`Screw` uses `SEn3Jac<S, 1>`. `Decoupled` and `WorldTranslation` use `TwistBlockJac<S>`, a newtype
over `ProductJac<Mat3<S>, Mat3<S>>` that implements `Jac<S, Twist<S>>` only, since their diagonal
blocks differ (a second `Jac` impl on `ProductJac` itself is E0283 at every concrete call; `0060`). Charts without group structure use dense matrices implementing
`Jac` (S²: `Mat2<S>`). Validated against `Dual` dual numbers.

## 2. S²

**NORMATIVE.** `S2<S>(Vec3<S>)`, `S2Chart<S> { base, b1, b2 }`, formulas `NUMERICS.md` §8.
Householder basis computed once in `Chart::at`, with $\varsigma = +1$ unless $n_z < 0$. `S2` is held
unit: `from_vec_unchecked`, `from_vec_normalized` and `renormalize` with the quaternion's bounds
(`NUMERICS.md` §12). Domain: `local(m)` requires $m \ne -n$.
Corpus ids: `s2_retract`, `s2_local`. Strata: `s2:nz0` (records at $n_z = +0$ and $-0$),
`s2:near-antipode`, `s2:generic`. `s2_local`'s reference is the geometric logarithm (Rodrigues, then
the quaternion `atan2`), never `mp.logm`. `s2_retract` is scored by `NUMERICS.md` §11's direction
metric ([`0063`](./decisions/0063-the-sphere-reads-its-sign-by-comparison-and-is-held-unit.md)).

## 3. Sim(3)

Awaits `NUMERICS.md` §9. `Sim3<S> { q, t, sigma }`, tangent `Sim3Tangent { phi, rho, sigma }`,
dense `Mat<7>` Jacobians ([`0005`](./decisions/0005-the-jacobian-is-a-dual-matrix.md)).
Corpus ids: `sim3_exp`, `sim3_log`, `sim3_jr`, `sim3_jr_inv`, `sim3_ad`.

## 4. Integrated exponentials

**NORMATIVE.** `so3::gamma1(φ) = jl(φ)`, `so3::gamma2(φ)` via `gamma2_coeffs` → $(b, d)$
(`NUMERICS.md` §7); corpus id `so3_gamma2`. Directional Jacobian computes
$\partial(\Gamma_m(\varphi)v)/\partial\varphi$ from `gamma_m`'s coefficients on `Dual<S, 1>` and the
chain rule ([`0066`](./decisions/0066-gamma-shares-its-coefficients-and-a-gaussian-factors-once.md)).

[`0064`](./decisions/0064-the-integrated-exponentials-reuse-the-swept-switches.md) fixes the shape:
- `gamma_apply_jacobians::<M, S>(&φ, v) -> (Γ_M v, ∂/∂φ, ∂/∂v = Γ_M)` for `M` 1 or 2, asserted
  at monomorphization ([`0066`](./decisions/0066-gamma-shares-its-coefficients-and-a-gaussian-factors-once.md) adds `Γ_M` and renames it).
- `gamma2_coeffs` takes no switch of its own.
- `so3_gamma2` has `so3_jl`'s strata and their `@f32` twins.
- The seeded defect "Γ₂ from Γ₁'s coefficients" is in `PHASE1.md` §10.

## 5. `Gaussian`

**NORMATIVE.**

```rust
pub struct Gaussian<S: Real, G: LieGroup<S>, Sd: Side, const D: usize> {
    pub mean: G,
    cov: Matrix<S, D, D>, // read through `cov()`; exactly symmetric (`0066`)
    _side: core::marker::PhantomData<Sd>,
}
```

- Enforces `const { assert!(D == G::DOF) }`.
- `to_left()` / `to_right()`: $\Sigma_L = \mathrm{Ad}_\mu\,\Sigma_R\,\mathrm{Ad}_\mu^\top$ via `Jac::sandwich`.
- `propagate(&self, j: &G::Jac, mean: G)`: $J\Sigma J^\top$.
- `mahalanobis_sq(&self, x: &G) -> (S, S::Mask)`: $\|L^{-1}(x \ominus_{Sd} \mu)\|^2$ with `chol` mask.
- `whitener(&self) -> (Whitener, S::Mask)`: $L$ once, then `Whitener::mahalanobis_sq(x) -> S` per
  point, `mahalanobis_sq`'s value to the bit ([`0066`](./decisions/0066-gamma-shares-its-coefficients-and-a-gaussian-factors-once.md)).
- Side is encoded in the type for invariant safety.

[`0065`](./decisions/0065-a-gaussian-names-its-side-and-is-stored-symmetric.md) fixes what this leaves open:
- `new` is the only constructor. Only `cov`'s lower triangle is read, and every value a method
  returns has its upper triangle copied from its lower one.
- `to_right` uses the adjoint of the computed inverse.
- `propagate`'s `j` is documented as the `Sd`-side Jacobian, which the type cannot check.
- The round trip `to_right ∘ to_left` stays within `γ_{4D+42}(|B||A|)|Σ|(|B||A|)ᵀ` componentwise.
- Corpus ids `gaussian_mahalanobis_se3` and `_se23`, at both precisions.

## 6. Consumer migrations

**NORMATIVE gates.** Tracked in downstream consumer repositories.
- **Vision/PnP pipelines:** Retractions become charts. Gate: discrete outputs and residual evaluations identical.
- **Tag/Pose estimators:** `Pose::retract` adopts `WithChart<SE3<f64>, Decoupled<f64>>`. Matrix representations convert cleanly.
- **Fusion estimators:** Direct consumption of `SE23`, $\Gamma$, both sides, and `Gaussian`.

## 7. Definition of done

- [ ] §1, §2, §4, §5 implemented with corpus, envelope, and dual-number twins.
- [ ] `NUMERICS.md` §9 completed by record; §3 implemented and dominated against oracles.
- [ ] Downstream migrations validated against regression gates.

## Appendix: suggested implementation order

1. `Chart`, `Manifold`, `WithChart`, group charts.
2. SE(3) charts and Jacobians.
3. Γ₂ and `Dual` directional Jacobian.
4. `Gaussian`.
5. S².
6. `NUMERICS.md` §9 record, then Sim(3).
7. Downstream migrations.
