# helicoid — Phase 7 Requirements: Continuous Time (`helicoid-spline`)

> **GATED by [`0011`](./decisions/0011-continuous-time-waits-for-a-consumer.md), not scheduled.**
> Requirements artifact, not an implementation authorization. Questions are answered by consumer records.

## 0.0 Gate

| Gate | Status |
|---|---|
| A consumer record exists naming the problem | **Unmet** |
| That record specifies: groups, knots, order, derivatives, precision, problem size, latency budget | **Unmet** |
| That record resolves §3.2 (baseline parity parameterization) | **Unmet** |
| The core owes nothing new (§1 delivered) | Met on Phase 3 completion |

Triggers: downstream camera–IMU calibration, continuous-time estimation, or `tf_tree` cumulative B-splines.

## 1. What the core already owes

`Ad`, `ad`, `Jr`, `Jr⁻¹`, `Jl`, `Jl⁻¹` on every group (Phase 3); geodesics and Jacobians (Phase 4);
`Dual` through all operations (Phase 2). No spline-shaped hook is added to core in advance.

## 2. Requirements sketch

- Crate `helicoid-spline`, same budget as `helicoid` ([`0007`](./decisions/0007-the-budget-a-foundation-can-afford.md)).
- Uniform cumulative B-splines of order `K` (const generic), basis matrix computed at compile time.
- Generic over Lie groups: SE(3) and `Product<SO3, R3>` share implementation.
- Value, body velocity, body acceleration by Sommer et al. recursions; $O(K)$ control-point Jacobians.
- Verification: `Dual` through splines, mpmath corpus fixtures, oracle benches.

## 3. Open design questions

1. **Non-uniform knots:** Uniform vs non-uniform spacing needs.
2. **Baseline parity:** Euclidean splines on rotation vectors vs Lie-group splines.
3. **`tf_tree` arena integration:** Stored control points vs sampled evaluation.
4. **Precision:** `f32` requirements.
5. **SE(3) parameterization:** Cumulative SE(3) vs split `SO3 × R3`.
6. **Gaussian-process priors:** STEAM-style alternative formulations.
