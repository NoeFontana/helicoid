# helicoid — Phase 5 Implementation Specification: Extended Geometry and Retraction Migrations

> **Companions:** `docs/NUMERICS.md` §7–§9, [`0012`](./decisions/0012-a-retraction-is-a-chart.md),
> [`0005`](./decisions/0005-the-jacobian-is-a-dual-matrix.md).

**Deliverable:** charts as types; S²; Sim(3); Γ₁/Γ₂ with directional Jacobians; `Gaussian` with a
stated side; and consumer migrations onto `helicoid` without unrecorded regressions. Sections marked
**NORMATIVE** are requirements.

## 0.0 Implementation status

| Area | Status |
|---|---|
| `Chart`, `Manifold`, `WithChart`; group charts (§1) | Not started |
| SE(3) charts `Screw`, `Decoupled`, `WorldTranslation` (§1.3) | Not started |
| S² and its chart (§2) | Not started |
| Sim(3) (§3) | Not started |
| Γ₁, Γ₂, directional Jacobians (§4) | Not started |
| `Gaussian<S, G, Sd, D>` (§5) | Not started |
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
$\partial(\Gamma_m(\varphi)v)/\partial\varphi$ via `gamma_m` on `Dual<S, 3>`.

## 5. `Gaussian`

**NORMATIVE.**

```rust
pub struct Gaussian<S: Real, G: LieGroup<S>, Sd: Side, const D: usize> {
    pub mean: G,
    pub cov: Matrix<S, D, D>,
    _side: core::marker::PhantomData<Sd>,
}
```

- Enforces `const { assert!(D == G::DOF) }`.
- `to_left()` / `to_right()`: $\Sigma_L = \mathrm{Ad}_\mu\,\Sigma_R\,\mathrm{Ad}_\mu^\top$ via `Jac::sandwich`.
- `propagate(&self, j: &G::Jac, mean: G)`: $J\Sigma J^\top$.
- `mahalanobis_sq(&self, x: &G) -> (S, S::Mask)`: $\|L^{-1}(x \ominus_{Sd} \mu)\|^2$ with `chol` mask.
- Side is encoded in the type for invariant safety.

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
