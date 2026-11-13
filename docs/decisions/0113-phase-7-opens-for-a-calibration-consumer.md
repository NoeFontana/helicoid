# 0113: Phase 7 opens for a calibration consumer

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** —

## Context

[`0011`](./0011-continuous-time-waits-for-a-consumer.md) gates `helicoid-spline`. The gate opens on a
consumer record that answers `PHASE7.md` §0.0's rows: groups, knots, order, derivatives, precision,
problem size and latency budget, and §3.2 (baseline parity).

A downstream camera–IMU calibration consumer now needs a continuous-time trajectory and bias
splines. Its records are private and are not cited from here. This record states its answers
itself, and the owner vouches for them as the consumer's.

What the consumer needs that ships today: `Exp`, `Log`, `Ad`, `ad`, `J_r`, `J_r⁻¹`, `J_l⁻¹` on SO(3),
SE(3), Rⁿ and products (Phase 3); `Dual` (Phase 2); `svd3`, `eig3`, `eig_sym`, `chol`; `S2` and
`S2Chart`; `Gaussian`. The spline itself does not ship.

## Decision

1. **The gate's rows, answered:**

   | Row | Answer |
   |---|---|
   | Groups | `Product<SO3, R3>` (the trajectory) and `Rⁿ` (biases). SE(3) comes from the same generic code (open question 2). |
   | Knots | Uniform. The consumer chooses the segment and `u` from its own integer time; gaps are separate splines. |
   | Order | `K = 4` to `6` for the groups (open question 3 for `K = 2, 3`). |
   | Derivatives | Value, velocity and acceleration, right-trivialized (body rate and body angular acceleration on SO(3)), with their Jacobians on every control. No third derivative. |
   | Precision | The consumer runs binary64 (open question 1 for binary32). |
   | Problem size | About `3 × 10⁵` evaluations per Levenberg–Marquardt iteration, windows of at most 6 controls, many samples per segment (an IMU stream). |
   | Latency budget | **None stated.** The bar is the bench gate's no-regress on committed rows, and the kernels take the fastest form the measurements support (`O(K)` backward accumulation, per-segment caching). No absolute number is set until a consumer profile shows evaluation dominating its iteration. |

2. **`PHASE7.md` §3's questions:**
   1. *Non-uniform knots*: no.
   2. *Baseline parity*: no parity spline. The consumer scores calibrated outputs against truth
      generated outside every spline space, so Kalibr's rotation-vector spline is part of Kalibr's
      method, not a confound. `helicoid` ships no rotation-vector spline.
   3. *`tf_tree` arena integration*: not decided here; it belongs to a `tf_tree` record. The window API
      ([`0114`](./0114-a-spline-evaluates-the-window-it-is-handed.md)) serves stored control points and
      sampled evaluation alike.
   4. *Precision*: generic over `Real`, as every `helicoid` routine is.
   5. *SE(3) or split*: generic over `LieGroup`. The consumer uses the split `Product<SO3, R3>`.
   6. *Gaussian-process priors*: not requested; they stay gated by `0011`.
3. **The work is four more records:** [`0112`](./0112-a-jacobian-sums-by-name.md) (the one core
   primitive), [`0114`](./0114-a-spline-evaluates-the-window-it-is-handed.md) (crate and API),
   [`0115`](./0115-the-spline-tables-are-exact-integers.md) (basis and Gram tables) and
   [`0116`](./0116-the-spline-is-measured-before-it-ships.md) (conformance). The formulas go into a new
   `NUMERICS.md` section with the first of them to become `ready`.

## Rationale

- **Restating, not citing.** `0011` asks for a consumer record so that the answers are written down
  before the API is. Here they are, without exposing a private repository.
- **No latency number.** A number the consumer cannot defend would be a guess, which `0011` exists
  to avoid. No-regress plus `O(K)` kernels is what "as fast as possible" can be checked against.

## Consequences

- When this is `ready`: `PHASE7.md` §0.0's three consumer rows read Met, citing it; `PHASE7.md` drops
  its GATED banner; `0011`'s *Implementation* line cites this record; `PROJECT.md` §5.1's Phase 7 row
  moves.
- Third derivatives, non-uniform knots and GP priors stay gated, each until a consumer asks.

## Implementation plan

1. The spec text of *Consequences*, in the PR that makes this `ready`. Verified by `just lint`.

## Open questions

1. **Binary32 spline ids.** Every `helicoid` id is scored at both precisions; this consumer is
   binary64 only. Keep the `@f32` strata (about half the spline corpus, `0116`), or record an
   exception and score binary64 only?
2. **SE(3) ids now.** The generic code makes SE(3) nearly free to build, but its corpus rows are the
   largest per record (`0116`), and `tf_tree`, the consumer that would use them, has not asked.
3. **Orders 2 and 3.** Generate and score them for Rⁿ and the Gram matrices, or only `K = 4..6`?
