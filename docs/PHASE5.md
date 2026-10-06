# helicoid — Phase 5 Implementation Specification: Extended Geometry and the Retraction Migrations

> **Companions:** `docs/NUMERICS.md` §7–§9, [`0012`](./decisions/0012-a-retraction-is-a-chart.md),
> [`0005`](./decisions/0005-the-jacobian-is-a-dual-matrix.md).

**Deliverable:** charts as types; S²; Sim(3); Γ₁/Γ₂ with directional Jacobians; `Gaussian` with a
stated side; and omnisac, locus-tag and locus_fusion moved onto `helicoid` without an unrecorded
change to any of their outputs. Sections marked **NORMATIVE** are requirements.

## 0.0 Implementation status

| Area | Status |
|---|---|
| `Chart`, `Manifold`, `WithChart`; group charts (§1) | Not started |
| SE(3) charts `Screw`, `Decoupled`, `WorldTranslation` (§1.3) | Not started |
| S² and its chart (§2) | Not started |
| Sim(3) — **blocked on `NUMERICS.md` §9** (§3) | Not started |
| Γ₁, Γ₂, directional Jacobians (§4) | Not started |
| `Gaussian<S, G, Sd, D>` (§5) | Not started |
| omnisac, locus-tag, locus_fusion migrations (§6) | Not started |

## 0. Non-goals and guardrails — read first

**NORMATIVE.** No preintegration (Γ is the primitive; the integrator is locus_fusion's), no
sampling or RNG in `Gaussian`, no solver-facing ambient Jacobians (Phase 6). **Sim(3) code does not
start until `NUMERICS.md` §9 is completed by record.**

## 1. Charts

**NORMATIVE** ([`0012`](./decisions/0012-a-retraction-is-a-chart.md)).

### 1.1 Traits

```rust
pub trait Chart<S: Real, M>: Copy {
    type Tangent: Tangent<S>;
    /// §1.4: the group's `Jac` where the tangent is a group's, a fixed dense matrix otherwise.
    type Jac: Jac<S, Self::Tangent>;
    /// Freeze everything the chart needs at `base` (e.g. S²'s basis).
    fn at(base: &M) -> Self;
    fn base(&self) -> M;
    fn retract(&self, delta: &Self::Tangent) -> M;
    fn local(&self, other: &M) -> Self::Tangent;
    /// ∂ retract(δ)/∂δ, expressed in the chart frozen at retract(δ).
    fn retract_jacobian(&self, delta: &Self::Tangent) -> Self::Jac;
    /// ∂ local(m)/∂m, with m perturbed in the chart frozen at m.
    fn local_jacobian(&self, other: &M) -> Self::Jac;
}

pub trait Manifold<S: Real>: Copy + Blend<S> {
    type Tangent: Tangent<S>;
    type Chart: Chart<S, Self, Tangent = Self::Tangent>;
    const DOF: usize;
}

/// `M` with chart `C` instead of `M::Chart`. Zero cost.
#[repr(transparent)]
pub struct WithChart<M, C>(pub M, core::marker::PhantomData<C>);
```

Every `LieGroup` is a `Manifold` with `Chart = RightChart<G>` (retract $= X\,\mathrm{Exp}(\delta)$,
local $= Y \ominus_R X$). `LeftChart<G>` is the left twin.

### 1.2 Why charts are types

A solver generic over `Manifold` linearizes in `M::Chart`; a consumer that wants a different
retraction spells it `WithChart<SE3<f64>, Decoupled>`, visible in the type and the review.

### 1.3 SE(3) charts

| Chart | `retract(δ = [φ; ρ])` | `local(Y)` | Used by |
|---|---|---|---|
| `Screw` (= `RightChart<SE3>`) | $X\,\mathrm{Exp}(\delta)$ | $\mathrm{Log}(X^{-1}Y)$ | default; locus_fusion, locus-calib |
| `Decoupled` | $(R\,\mathrm{Exp}(\varphi),\ t + R\rho)$ | $(\mathrm{Log}(R^\top R_Y),\ R^\top(t_Y - t))$ | locus-tag `Pose::retract` |
| `WorldTranslation` | $(R\,\mathrm{Exp}(\varphi),\ t + \rho)$ | $(\mathrm{Log}(R^\top R_Y),\ t_Y - t)$ | `Product<SO3, R3>`'s chart |

They agree to first order at $\delta = 0$ and differ at second order, which is exactly why a silent
swap changes LM iterates.

### 1.4 Jacobian types for charts

`Screw` is `RightChart<SE3>` and uses `SEn3Jac<S, 1>`. `Decoupled` and `WorldTranslation` do
**not** have dual-matrix Jacobians — their diagonal blocks differ (e.g. `Decoupled`'s
`retract_jacobian` is $\mathrm{diag}(J_r(\varphi),\ \mathrm{Exp}(-\varphi))$) — so they use the
block-diagonal `ProductJac<Mat3<S>, Mat3<S>>`. Charts without group structure use a fixed dense
matrix implementing `Jac` (S²: `Mat2<S>`). Each chart's Jacobians are checked against `Dual`
through `retract`/`local` (`chart_jacobians_match_dual_*`).

## 2. S²

**NORMATIVE.** `S2<S>(Vec3<S>)`, `S2Chart<S> { base, b1, b2 }`, formulas `NUMERICS.md` §8; the
Householder basis is computed once in `Chart::at`. Domain: `local(m)` requires $m \ne -n$.
Corpus ids `s2_retract`, `s2_local` (reference: minimal rotation via `mp.expm`/`mp.logm`); strata
`s2:nz0` ($n_z = \pm 0$, the basis discontinuity), `s2:near-antipode`, `s2:generic`, and the
`theta:*` decades for $\|\delta\|$.

## 3. Sim(3)

**Blocked on `NUMERICS.md` §9.** Then: `Sim3<S> { q, t, sigma }`, tangent `Sim3Tangent { phi, rho,
sigma }`, dense `Mat<7>` Jacobians (0005); corpus ids `sim3_exp`, `sim3_log`, `sim3_jr`,
`sim3_jr_inv`, `sim3_ad` (reference: `mp.expm`/`mp.logm` of the 4×4 similarity matrix, series of
$\mathrm{ad}$); strata on the grid $\sigma \in \{0, \pm10^{-12}, \dots, \pm10^{-1}, \pm1, \pm3\}$ ×
the `theta:*` decades, **including the joint limit** $\sigma^2 + \theta^2 \to 0$ along several rays.

## 4. Integrated exponentials

**NORMATIVE.** `so3::gamma1(φ) = jl(φ)`, `so3::gamma2(φ)` via `gamma2_coeffs` → $(b, d)$
(`NUMERICS.md` §7); corpus id `so3_gamma2` (series $\sum W^n/(n+2)!$).
`so3::gamma_apply_jacobian::<M>(φ, v) -> Mat3` computes $\partial(\Gamma_m(\varphi)v)/\partial\varphi$ by
evaluating `gamma_m` on `Dual<S, 3>` — **this is the implementation, not a test aid**, until the
ledger row in `PROJECT.md` §5.1 is opened by a bench.

## 5. `Gaussian`

**NORMATIVE.**

```rust
pub struct Gaussian<S: Real, G: LieGroup<S>, Sd: Side, const D: usize> {
    pub mean: G,
    pub cov: Matrix<S, D, D>,
    _side: core::marker::PhantomData<Sd>,
}
```

- Construction asserts `const { assert!(D == G::DOF) }` (a post-monomorphization error, stable).
- `to_left()` / `to_right()`: $\Sigma_L = \mathrm{Ad}_\mu\,\Sigma_R\,\mathrm{Ad}_\mu^\top$ (and inverse),
  via `Jac::sandwich`; reference twin dense.
- `propagate(&self, j: &G::Jac, mean: G)`: $J\Sigma J^\top$.
- `mahalanobis_sq(&self, x: &G) -> (S, S::Mask)`: $\|L^{-1}(x \ominus_{Sd} \mu)\|^2$ with
  `chol`'s positive-definiteness mask.
- The side is part of the type: a `Gaussian<_, SE3<f64>, Right, 6>` cannot be passed where `Left`
  is expected. This is the calibration → VIO handoff's whole safety argument.

## 6. Consumer migrations

**NORMATIVE gates.** Each migration is a PR in the consumer, citing this section.

- **omnisac** (PnP, essential, rigid refinement): retractions become charts. Gate: RunRecords'
  discrete outputs and `residual_evals` identical; values within the conformance bars. A chart
  change (as opposed to a port) is an omnisac record.
- **locus-tag**: `Pose::retract` becomes `WithChart<SE3<f64>, Decoupled>`; its `[v, ω]` LM delta
  crosses into `helicoid` through `Twist::from_translation_first` at the LM boundary only;
  `quat_from_so3` becomes `SO3::from_matrix` and its Müller-hang regression test moves with it;
  nalgebra types convert through `mint`. Gate: the paired-CI comparison under `tools/bench/`
  (`compare/`, `metrics.py`, `strata.py`), driven through `tools/cli.py bench` — no regression on
  every reported metric, **rotation p95 included**. Bit identity is not expected (nalgebra's
  `Rotation3::new` vs `SO3::exp` differ in rounding), so the gate is statistical by design.
  Measured scope, as of [`0049`](./decisions/0049-the-boundary-is-what-removes-a-way-to-be-wrong.md) decision 6: the whole *production* Lie surface is five
  items — `Pose::{retract, inverse, adjoint}`, the two covariance reframes and `quat_from_so3` (two
  non-test callers, both in `locus-py`) — and everything else matching `retract`/`Rotation3` under
  `crates/locus-core/src/` is a test. `Pose::adjoint` exists **only** to sandwich a 6x6, so
  `Jac::sandwich` replaces it and the dense `Matrix6` never materializes (`0005`). The unrecorded
  constraint is representation: `Pose`'s `rotation` field is a public `Matrix3<f64>` where
  `helicoid` is quaternion-backed, so either the retract converts per LM iteration — up to twenty
  per pose, each paying `from_matrix`'s pivot selection — or `Pose`'s representation changes, which
  is a locus-tag published-API change. Which one is a measurement on that gate, not a reading.
- **locus_fusion**: if the estimator is Rust, it uses `SE23`, Γ, both sides and `Gaussian`
  directly. If it stays C++, its test suite consumes the corpus JSONL through a Sophus-order
  converter; a failing C++ stratum is a locus_fusion bug report, not a `helicoid` change.

## 7. Definition of done

- [ ] §1, §2, §4, §5 implemented with corpus, envelope and twins; baseline blessed.
- [ ] `NUMERICS.md` §9 completed by record; §3 implemented and dominated where oracles exist
      (Sophus `Sim3`).
- [ ] The three migrations merged in their repositories with their gates quoted.

## Appendix: suggested implementation order

1. `Chart`, `Manifold`, `WithChart`, group charts.
2. SE(3) charts and their Jacobians (locus-tag's migration can start here).
3. Γ₂ and the `Dual` directional Jacobian.
4. `Gaussian`.
5. S².
6. `NUMERICS.md` §9 record, then Sim(3).
7. Consumer migrations.
