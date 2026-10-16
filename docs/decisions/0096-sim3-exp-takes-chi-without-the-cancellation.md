# 0096: `Sim3::exp` takes `χ` without the cancellation

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

`Sim3::exp` applies `V = v₀I + v₁W + v₂W²` by two cross products (`NUMERICS.md` §9):
`v₀ρ + v₁ (φ × ρ) + v₂ (φ × (φ × ρ))`.
- **The cancellation.** Its `ρ` part is `χρ`, with `χ = v₀ − zv₂ = Re g₁(σ + iθ)` (SM.5(c)). At
  `σ = 0` that is `sin θ/θ`, formed out of `v₀ = 1` and `zv₂ → 1` near `π`: the cancellation of
  [`0095`](./0095-exp-translation-takes-the-far-diagonal.md) (draft) in `SEn3::exp`, carried by a
  coefficient instead of a diagonal. Elsewhere the two terms part by less, but `zv₂` is the larger
  wherever `|σ|` is not small against `θ`.
- **The standing.** On the deep tier, `helicoid` is above the best oracle on 25 of `sim3_exp`'s 122
  `f64` strata. Its worst is `sigma:10/theta:3`, at 10.60 `u` against Sophus's 5.86.

## Decision

1. **The group's coefficients carry `χΔ`** (`Sim3Coeffs::chi_delta`), on the closed arm as
   `σ²v₀ + e^σ z sinc θ − σe^σ za`.
   - **Where it comes from:** SM.5(c)'s `Δ Re g₁ = σ(e^σ cos θ − 1) + θe^σ sin θ`, with
     `e^σ cos θ − 1 = m₁ − e^σ za`. It is scaled as every other closed-arm scalar is.
   - **Order:** the two non-negative terms are added first, then the one of `−σ`'s sign.
   - **On the series arm** it is `(v₀ − zv₂)Δ`, the definition. No reader takes that arm.
   - **No new transcendental or switch:** every input is one the arm already holds.
2. **`Sim3::exp` takes the same far arm as `SEn3::exp`.** Where `b` is on its exact arm, `Vρ` is
   `χρ + v₁ (φ × ρ) + v₂ (φ·ρ) φ`, with `χ = χΔ/Δ` divided inside that arm.
   - **The mask:** `coeffs::b_takes_exact_arm`, a predicate on `b`'s switch, as `0051`'s
     `log_ratio_takes_short_arm` is on `r`'s. It is the mask `jr_coeffs_sinc` returns, so no switch
     is typed.
   - **Below the switch,** the bits are unchanged.
   - **`exp_and_jl`'s `Exp`** takes the same arm, so it stays `exp`'s to the bit.
3. **`sim3_v` returns `(v, χΔ)`.** `NUMERICS.md` §9 states the far arm.
4. **Unchanged:** `Log` and the inverse Jacobians keep `v₀ − zv₂` (see *Rejected*).

## Rationale

**Accuracy.** Every subject was rescored on both tiers, against `main` at `3b031d2`.
`cargo xtask determinism` moves exactly `sim3_exp`, at both precisions.

| `sim3_exp` | records better / worse | strata maxima better / worse | total error | above best oracle |
|---|---|---|---|---|
| deep, `f64` | 864 / 492 | 22 / 4 | 0.954× | 25 → 19 of 122 |
| deep, `f32` | 763 / 432 | 15 / 10 | 0.960× | 13 → 12 of 116 |
| committed, `f64` | 57 / 32 | 14 / 8 | 0.968× | 33 → 29 of 122 |
| committed, `f32` | 47 / 26 | 12 / 8 | 0.948× | 9 → 10 of 116 |

- **The worst strata come in,** on the deep tier at `f64`:
  - `sigma:10/theta:3`: 10.60 → 6.87 `u`;
  - `sigma:100/theta:3`: 10.56 → 7.88;
  - `sigma:3/theta:pi-1e-6`: 8.59 → 5.87;
  - `sigma:3/theta:3`: 7.23 → 5.75;
  - `sigma:1e-2/theta:3`: 4.80 → 2.73.
- **The committed tier's worst moves the other way** (`sigma:30/theta:3`, 5.68 → 7.13 `u`), on
  4 records. On the deep tier's 64 the same cell reads 7.66 → 7.13.
- **The envelope.** Deep tier: 2 domination failures before and after, with unresolved 65 → 63
  and unresolvable 151 → 147. Committed tier: 0 failures and 23 unresolved either way,
  unresolvable 505 → 502.

**Latency.** `Sim3::exp` in a consumer loop, 4096 tangents per band, best of 40, one pinned core,
`ρ ~ U[−2, 2]³`, `σ ~ σ̄·U[0.5, 1.5]`, `main` → this (ns):

| | `θ ~ 10⁻⁶..10⁻²` | `0.5..2` | `2.8..3` | `3.11..3.14` |
|---|---|---|---|---|
| `f64`, `σ̄ = 10⁻³` | 36.1 → 37.1 | 39.0 → 41.2 | 60.6 → 63.2 | 55.5 → 58.2 |
| `f64`, `σ̄ = 1` | 45.0 → 47.3 | 52.0 → 54.7 | 75.0 → 77.5 | 69.7 → 71.8 |
| `f32`, `σ̄ = 10⁻³` | 25.1 → 27.1 | 31.5 → 32.8 | 38.6 → 41.0 | 33.7 → 35.9 |
| `f32`, `σ̄ = 1` | 35.0 → 37.5 | 38.9 → 41.2 | 54.0 → 54.6 | 51.1 → 49.1 |

That is +3 to +8%, and −4% at `f32` near `π`. `χΔ` in the kernel costs nothing measurable: with
`Exp`'s old assembly the two builds read the same. The cost is the branch itself, as `0095`'s was
before its quaternion moved, and no ordering found here removed it. **This is the trade the owner
decides:** about 4.5% of `sim3_exp`'s total error and its worst strata, for 1–3 ns.

**Rejected:**
- **`χ` inside the kernel's division.** Same accuracy, but 1 ns slower on the far bands, where the
  third division contends with `v₂`'s and `v₁`'s.
- **Branch-free per-arm scalars** (`0095`'s five-term formula): +4 ns at `f64` and up to +17 ns at
  `f32`.
- **`t` formed before `q`:** no better.
- **`χ` in `Log`'s `V⁻¹`** (`Λ = χ² + zv₁²`, `i₂`), for `sim3_log`, `sim3_jr_inv` and `sim3_jl_inv`.
  Their means moved by ±1% and their records split evenly. `V⁻¹`'s conditioning, not `χ`'s
  rounding, sets their error, as `0091` found for `SEn3::log`'s `β` form.

## Consequences

- **`Sim3Coeffs` carries a ninth scalar.** Every arm forms it. `sim3_v`'s callers read `(v, χΔ)`.
- **`b_takes_exact_arm` is a second exported predicate** on a generated switch.
  `the_rays_meet_the_recorded_objective`'s scan now also holds `sim3_v`'s `χΔ` to the full group's
  bits.

## Implementation plan

1. `Sim3Coeffs::chi_delta`, `coeffs::b_takes_exact_arm`, `Sim3::exp` and `exp_and_jl` through
   `v_apply`, `NUMERICS.md` §9. Verified by:
   - `just test`, including `chi_delta_is_v0_minus_z_v2_without_the_cancellation` (`sin θ/θ` to
     `4u` at `σ = 0` near `π`, where the subtraction reads about `10⁶ u`);
   - `just lint`, `just doc`, `just msrv`, `just no-std`;
   - `cargo xtask determinism` against `main`;
   - `cargo xtask envelope` with and without `--deep`.

## Open questions

1. **Whether 1–3 ns of `Sim3::exp` buys 4.5% of its error and its worst strata.** That is the
   owner's call.
