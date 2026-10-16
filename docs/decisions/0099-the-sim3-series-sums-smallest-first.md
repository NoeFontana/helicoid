# 0099: The Sim(3) series sums smallest first

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

`Sim3::log` was the largest remaining gap to the oracles. On the deep tier it is above the best
oracle on 24 of 119 `f64` strata and 25 of 104 `f32` ones, and on 4-record committed cells up to
3.2×.
- **Where the error is.** A per-field split of `Sim3::log` against the corpus reproduces the corpus
  scores exactly, and on the worst records it is in `ρ = V⁻¹t`, not `φ`.
  - At `sigma:0.3/theta:1e-8`, 4.80 `u`, all three components of `ρ` carry the same relative error,
    `+4.5 u`. That is a scale, `1/v₀`.
- **Why `v₀`.** It is on the group's series arm (`Δ = 0.09 < r_×² = 0.237`), where
  `coeffs::sim3::series_at` sums `v₀ = Σ σʲ/(j+1)!` and its siblings **upward**, largest term first.
  Each of the `K = 13` additions rounds at the sum's scale, about `1.17`.
- **What depends on it.** `v₀` is `V`'s eigenvalue on the axis, so `Exp`, `Log` and the four
  Jacobians all carry it.

## Decision

1. **`v₀ = g₁(σ)` and `w₀ = g₂(σ)` are summed by Horner** on the series arm (`g_series`, the closed
   arm's own). They are polynomials in `σ` alone. The recurrence for `ξ`, `ζ`, `ξ'`, `ζ'` and their
   sums is unchanged.
2. **The sweep regenerates `K` and `r_×²`** (`cargo xtask thresholds`, `0004`). With `v₀` accurate,
   the series beats the closed arm much further out:

   | | `K` | `r_×²` | objective (blocks' weighted max) |
   |---|---|---|---|
   | `f64` | 13 → 19 | 0.237 → 1.98 | 4.72 → **2.32** `u` |
   | `f32` | 9 → 11 | 0.835 → 2.37 | 2.93 → **2.26** `u` |

   The old switch sat at `0.237` because the forward-summed series lost to the closed arm beyond
   it, and the closed arm near small `Δ` read about 4 `u`.
3. **The series arm has four tiers.** Below three generated switches in `Δ` it runs at shallower
   depths, chosen by the sweep (`short_arm`).
   - **The criterion:** a tier's depth must score within `1/16 u` of depth `K` at **every** sweep
     point under it, so its truncation is negligible point by point.
   - **The choice:** of the valid candidates, the fewest terms summed over the arm's sweep points,
     each tier searched above the one below it.
   - **The result** (switches in `Δ`, depth below them):

     | | `≤ 5.4·10⁻⁵` | `≤ 0.042` | `≤ 0.26` | `≤ 1.98` |
     |---|---|---|---|---|
     | `f64` | 6 | 11 | 14 | 19 |

     | | `≤ 9.6·10⁻⁴` | `≤ 0.13` | `≤ 0.63` | `≤ 2.37` |
     |---|---|---|---|---|
     | `f32` | 4 | 7 | 9 | 11 |

   - **In code:** `Sim3Switch` gains `short`, `mid` and `upper` pairs. `coeffs::sim3::tiered` takes
     one branch per tier, so each depth is a constant its loop unrolls by.
4. **`NUMERICS.md` §9** states the Horner sums and the tiers. `docs/maths/sim3.md`'s SM.6–SM.7
   measurements describe the upward sum and are not re-measured here.

## Rationale

**Accuracy.** Every subject was rescored on both tiers. `cargo xtask determinism` against `main` at
`2b5a8a6` moves exactly `sim3_exp`, `sim3_log`, `sim3_jr`, `sim3_jl`, `sim3_jr_inv` and
`sim3_jl_inv`, at both precisions. `sim3_ad` and `sim3_act` read no coefficient.

| deep tier | records better / worse | total error | strata above best oracle |
|---|---|---|---|
| `sim3_log` `f64` | 1009 / 518 | 0.898× | 24 → 18 of 119 |
| `sim3_log` `f32` | 868 / 494 | 0.919× | 25 → 20 of 104 |
| `sim3_exp` `f64` | 1023 / 483 | 0.922× | 19 → 13 of 122 |
| `sim3_exp` `f32` | 872 / 414 | 0.939× | 12 → 8 of 116 |
| the four Jacobians | 1.8–2.2 : 1 | 0.95–0.98× | 0 → 0 |

- **`sim3_log`'s worst deep stratum** falls from 5.14 to 4.31 `u` (`f64`). The committed tier
  reads the same way, for example `sim3_log` `f32` at 0.860× total error with strata above the best
  oracle 26 → 21.
- **Envelope.** Deep tier: 1 domination failure, the one `0098` excepts, before and after, with
  unresolved 60 → 51 and unresolvable 140 → 138. Committed tier: 0 failures and 23 unresolved,
  unresolvable 498 → 485.

**Latency.** `Sim3::exp` and `Sim3::log` in a consumer loop, 4096 records per band, best of 40,
one pinned core, `σ ~ σ̄·U[0.5, 1.5]`, `main` → this (ns, `Exp`/`Log`):

| `f64` | `θ ~ 10⁻⁶..10⁻²` | `0.05..0.2` | `0.5..2` | `2.8..3` |
|---|---|---|---|---|
| `σ̄ = 10⁻³` | 34.1/50.5 → 21.0/34.2 | 37.0/67.0 → 36.2/65.5 | 41.7/80.7 → 57.3/94.5 | 64.0/93.0 → 64.2/93.0 |
| `σ̄ = 0.1` | 34.1/50.5 → 32.6/49.1 | 37.0/67.2 → 38.3/67.9 | 41.8/81.0 → 57.5/94.5 | 64.0/93.2 → 64.0/92.6 |
| `σ̄ = 1` | 48.3/60.5 → 56.7/72.1 | 55.3/80.5 → 60.9/89.5 | 55.3/87.8 → 62.1/93.6 | 77.8/106.4 → 78.0/106.4 |

| `f32` | `θ ~ 10⁻⁶..10⁻²` | `0.05..0.2` | `0.5..2` | `2.8..3` |
|---|---|---|---|---|
| `σ̄ = 10⁻³` | 25.8/38.4 → 18.5/25.6 | 26.0/54.7 → 23.6/51.6 | 33.4/64.8 → 40.3/70.4 | 42.7/64.9 → 42.7/64.6 |
| `σ̄ = 0.1` | 25.8/38.3 → 23.5/35.7 | 26.0/54.6 → 23.5/51.5 | 33.3/65.0 → 38.1/68.7 | 42.7/64.8 → 42.6/64.7 |
| `σ̄ = 1` | 38.8/55.5 → 35.0/50.5 | 42.5/71.7 → 36.5/64.4 | 44.7/79.5 → 44.3/76.0 | 58.6/83.0 → 58.8/83.3 |

- **Near the identity it is faster:** 38% at `f64`, 28% at `f32`.
- **It costs where `Δ ∈ [0.26, 1.98)`:** the 19-term series that replaced the closed arm, 37% on
  `f64`'s `0.5..2` band at small `σ` and 17% at `σ̄ = 1`. That is where the closed arm read 4 `u`.
- **The trade, for the owner:** the regenerated `r_×²` buys the objective's halving, and the
  `f64` series beyond `Δ ≈ 0.26` is the price.

**Rejected:**
- **Every sequence summed backward** (store the levels, sum down). Within `0.002×` of this record's
  accuracy at the old switches, but the 24 × 8 buffer spills: the series arm went 34 → 65 ns.
- **`v₀` by Horner in two chains** (`E(σ²) + σO(σ²)`). Half the latency of one chain, but half the
  gain (`sim3_log` deep 0.951× against 0.916× at the old switches).
- **One depth `K = 19` for the whole series arm.** That is the sweep's own choice, but +65% near
  the identity, where most calls are.
- **A selected depth instead of a branch per tier.** The loop bound stops being a constant, and
  near the identity it read 5 ns slower.
- **Tiers within `SLACK` of `K`'s maximum.** `f32` `sim3_log` read 1.04× `main`'s total error: a
  shallow tier raised typical errors while keeping the worst.
- **Keeping the old switches with the Horner sums.** `just thresholds-check` fails it: the switches
  are the sweep's to choose (`0004`).

## Consequences

- **`generated_sim3.rs` and `thresholds-sim3.csv` are regenerated,** and the CSV gains six
  columns: `short_terms`, `short_below`, `mid_terms`, `mid_below`, `upper_terms`, `upper_below`.
- **`0096`'s note that no reader of `χΔ` takes the series arm still holds,** since `b`'s exact arm
  starts at `θ² = 7.77` (`f64`) and `8.98` (`f32`), above both `r_×²`.

## Implementation plan

1. The Horner sums, `Sim3Switch`'s tiers, `tiered`, the sweep's `short_arm`, the regenerated
   files, and `NUMERICS.md` §9. Verified by:
   - `just thresholds-check`;
   - `just test`, including `the_rays_meet_the_recorded_objective`;
   - `just lint`, `just doc`, `just msrv`, `just no-std`;
   - `cargo xtask determinism` against `main`;
   - `cargo xtask envelope` with and without `--deep`.

## Open questions

1. **Whether the `f64` series beyond `Δ ≈ 0.26` is worth its latency.** It costs 17–37% there, and
   buys the halved objective. That is the owner's call. A cheaper deep series, a bivariate polynomial
   in `(σ, z)` with generated coefficients, would be a record of its own.
