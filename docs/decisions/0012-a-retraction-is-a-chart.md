# 0012: A retraction is a chart

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

> **Amended by [`0060`](./0060-the-charts-go-first-and-name-their-frame.md), 2026-10-09.** "They
> agree to first order" holds for `Screw` and `Decoupled` only. `WorldTranslation` agrees with
> `Decoupled` after the linear map `diag(I, R)` on the translation tangent, exactly, and not to
> first order unless `R = I` (`docs/maths/charts.md` CH.4(b), CH.5(c)). The Decision is unchanged.

The stack already uses three SE(3) retractions: `tf_tree_math`'s $X\,\mathrm{Exp}(\xi)$, locus-tag's
`Pose::retract` $(R\,\mathrm{Exp}(\omega),\ t + R v)$ — chosen deliberately, its rustdoc cites a
`JᵀWJ` conditioning improvement — and the product retraction $(R\,\mathrm{Exp}(\omega),\ t + v)$
implied by `Product<SO3, R3>`. They agree to first order and differ at second, so replacing one with
another changes LM iterates, RunRecords and bench tails without changing any test that checks the
fixed point. S² adds a second reason: by the hairy ball theorem no continuous global tangent basis
exists, so its chart must be frozen per linearization. InEKF-style filters need left charts; Ceres-
style solvers need ambient Jacobians per chart (Phase 6).

## Decision

1. A chart is a type implementing `Chart<S, M>` (`PHASE5.md` §1.1): `at(base)` freezes it;
   `retract`, `local` and their Jacobians use only what was frozen.
2. `Manifold<S>` carries a default `Chart`; every `LieGroup` is a `Manifold` with
   `RightChart<G>`. `LeftChart<G>` exists for every group.
3. `WithChart<M, C>` is a zero-cost newtype selecting a non-default chart, so the choice is visible
   in the type: locus-tag's LM variable is `WithChart<SE3<f64>, Decoupled>`.
4. SE(3) has exactly three named charts: `Screw` (the default, `= RightChart<SE3>`), `Decoupled`,
   `WorldTranslation` (`PHASE5.md` §1.3). A fourth is a record.
5. S²'s chart freezes the Householder basis at `at` (`NUMERICS.md` §8).

## Rationale

Making the retraction a type turns a silent semantic change into a reviewable diff, and lets a
consumer migrate onto `helicoid` without changing its numerics (locus-tag keeps `Decoupled`). The
alternatives lost: a runtime enum (a branch on every retract, and the choice invisible in
signatures); one blessed retraction per group (forces behaviour changes on migration).

## Consequences

- Solvers are generic over `Manifold` and linearize in `M::Chart`.
- Chart Jacobians are `Dual`-checked like group Jacobians.

## Implementation plan

`PHASE5.md` appendix steps 1–2 — verified by `chart_jacobians_match_dual_*` and the corpus ids of
`PHASE5.md` §2.

## Open questions

None.
