# helicoid — Phase 7 Implementation Specification: Continuous Time (`helicoid-spline`)

> **Opened by [`0113`](./decisions/0113-phase-7-opens-for-a-calibration-consumer.md)**, which answers
> [`0011`](./decisions/0011-continuous-time-waits-for-a-consumer.md)'s gate for a calibration
> consumer. The formulas are `NUMERICS.md` §21; the records are `0112`–`0116`.

## 0.0 Implementation status

| Area | Status |
|---|---|
| Gate: a consumer record names the problem | Met ([`0113`](./decisions/0113-phase-7-opens-for-a-calibration-consumer.md)) |
| Gate: groups, knots, order, derivatives, precision, problem size, latency budget | Met (`0113` decision 1) |
| Gate: §3.2, baseline parity | Met (`0113` decision 2: no parity spline) |
| Gate: the core owes nothing new | Met, but one primitive: `Jac::add`, `Jac::scale` ([`0112`](./decisions/0112-a-jacobian-sums-by-name.md)) |
| `Jac::add`, `Jac::scale` on every impl (`0112`) | Not started |
| Corpus ids, references, seeded defects ([`0116`](./decisions/0116-the-spline-is-measured-before-it-ships.md)) | Not started |
| Crate, `Order`, `Weights`, tables, `derivative_gram`, twins ([`0114`](./decisions/0114-a-spline-evaluates-the-window-it-is-handed.md) step 1, [`0115`](./decisions/0115-the-spline-tables-are-exact-integers.md)) | Not started |
| `Segment`: value and kinematics (`0114` step 2) | Not started |
| `JacSegment`: control Jacobians, both sides (`0114` step 3) | Not started |
| Shipped subject answers the spline ids (`0114` step 4) | Not started |
| basalt-headers oracle (`0116` step 3) | Not started |
| Bench rows and baseline (`0116` step 4) | Not started |

## 1. What the core already owes

`Ad`, `ad`, `Jr`, `Jr⁻¹`, `Jl`, `Jl⁻¹` on every group (Phase 3); geodesics and Jacobians (Phase 4);
`Dual` through all operations (Phase 2). The one addition is `0112`'s named Jacobian sum and scale.
No spline-shaped hook is added to the core.

## 2. Scope

**NORMATIVE.** As `0114`:
- the crate `helicoid-spline`, under `0007`'s budget: `helicoid` and `helicoid-linalg` only;
- uniform cumulative B-splines of order `K = 2..6`, with exact basis tables computed at compile time
  (`0115`);
- generic over `LieGroup`, so that SE(3) and `Product<SO3, R3>` share one implementation;
- value, body velocity and body acceleration by Sommer et al.'s recursions, with `O(K)` control
  Jacobians on either side.

Verification is `0116`'s: the corpus ids, mpmath references, basalt-headers, seeded defects and bench
rows.

## 3. Design questions

All six are answered in `0113` decision 2. Non-uniform knots, third derivatives and
Gaussian-process priors stay gated by `0011`, each until a consumer asks.
