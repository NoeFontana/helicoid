# 0064: The integrated exponentials reuse the swept switches

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** (filled in as work lands)

## Context

`PHASE5.md` §4 asks for `so3::gamma1(φ) = jl(φ)`, `so3::gamma2(φ)` from `gamma2_coeffs → (b, d)`,
a directional Jacobian `∂(Γ_m(φ) v)/∂φ` through `Dual<S, 3>`, and the corpus id `so3_gamma2`. A
fusion estimator's IMU increments are `Δv = Δt Γ₁(ωΔt) f` and `Δp = Δt² Γ₂(ωΔt) f` (GG.4), at
`θ = ‖ω‖Δt` between `3.5e-3` and `1.7e-2` for 200 °/s at 200 Hz to 1 kHz, which is where the exact
arm of `b` loses digits in `Γ₂` (GG.3(b)). Every formula is derived and checked in
`docs/maths/gamma-gaussian.md` GG.1–GG.6. Four questions were open (`docs/maths/index.md`):

1. Which `m` does `gamma_apply_jacobian` support?
2. What length and tolerance does the §14 twin's "dense sum" use?
3. Does `NUMERICS.md` §12 owe `Γ` a row?
4. `NUMERICS.md` §7 calls the `Dual` Jacobian "exact to rounding, Taylor branches included", but
   GG.6 measures the derivative of the computed function instead: worse than the value by `≈ 2K/θ`
   on a `K`-term series arm, and by `u/θ` (`m = 1`) or `u/θ²` (`m = 2`) on an exact arm.

`PHASE5.md` §4 also names no strata for `so3_gamma2`. `0016` keeps every vector id at binary64
"until a record asks for one".

**The switches are per coefficient, not per group.** `coeffs/generated.rs` holds one
`Switch<S, M>` per coefficient and per precision (`B_F64`, `D_F64`, …), and a call-site group is
`kernel::grouped` over its members' arms. A group `(b, d)` therefore needs no sweep and adds no
`generated.rs` row. Its short prefix, the smallest second switch among its members, is `d`'s at
both precisions (`2.64e-3` at `f64`, `1.38e-3` at `f32`), and that is also `q_coeffs`' `(b, d, e)`,
since `e`'s is larger. So the group's members are `q_coeffs`' `b` and `d` to the bit.

### Measured

**The twin.** 20 000 draws: `θ` log-uniform on `[1e-8, π]`, a uniform axis, `v` uniform in
`[−1, 1)³`, plus `φ = 0` and a subnormal `φ`. The error is `‖ŷ − y‖ / max(‖y‖, 1)` in `u`, the
worst per band. An `f32` subject is compared with the twin at `f64`, from the same `f32` inputs. An
`f64` subject is compared with the `f64` twin (D6's comparison, so it carries the twin's own
rounding).

| `gamma_apply_jacobian::<M>` | `θ < 1e-2` | `< 1e-1` | `< 1` | `≤ π` |
|---|---|---|---|---|
| `M = 1`, `f64`: value / Jacobian | 2.44 / 1.94 | 3.28 / 2.30 | 5.49 / 3.63 | 7.83 / 4.75 |
| `M = 2`, `f64` | 1.44 / 1.03 | 2.08 / 1.15 | 3.08 / 1.37 | 4.01 / 1.42 |
| `M = 1`, `f32` | 1.19 / 1.28 | 1.26 / 1.16 | 1.63 / 1.35 | 3.67 / 4.93 |
| `M = 2`, `f32` | 0.72 / 0.76 | 0.77 / 0.70 | 1.00 / 0.72 | 1.56 / 1.86 |

The `f64` value rows at `M = 1` are the twin's rounding: it sums terms up to `θ²/2` against a
result near 1. The `f32` rows are errors. The generated switches of `b` and `d` lie at `θ² ≈ 9.6`
(`f64`) and `7.8`, `9.6` (`f32`), so on `θ ≤ π` the code that runs is the series arm, and GG.6's
`620 u` at an illustrative `0.12` switch does not arise. `Γ₂`'s matrix against the twin's
columns: at most 2.95 `u` (`f64`) and 1.84 `u` (`f32`).

**The corpus.** `so3_gamma2`, over `so3_jl`'s strata and their `@f32` twins, scores at most
**1.46 `u` at `f64`** (`theta:pi-1e-3`) and **1.68 `u` at `f32`** (`theta:pi-1e-4@f32`). Its
reference is the series `Σ Wⁿ/(n+2)!` at 120 digits, cross-checked by GG.2(b)
(`I + Γ₂W = Γ₁`), which pins `Γ₂` off the kernel of `W`, and by `Γ₂φ = φ/2`, which pins it on that
kernel.

## Decision

1. **`gamma2_coeffs(z) -> (b, d)`** is a `pub(crate)` group in `coeffs`, built as
   `grouped([B.arm(), D.arm()], exact_b_d, z)`. Its exact arm shares one `θ` and computes `sin θ`
   and `sin(θ/2)`, without a cosine (`0052`). It adds no switch, no series and no `generated.rs`
   row, so `cargo xtask thresholds` writes the same bytes. Its members are `q_coeffs`' to the bit,
   value and derivative, at both precisions.
2. **`pub mod so3`** holds three functions. `SO3` and `SO3Tangent` stay at the crate root too.
   - `gamma1<S: Real>(phi: &SO3Tangent<S>) -> Mat3<S>` is `SO3::jl`, to the bit.
   - `gamma2<S: Real>(phi: &SO3Tangent<S>) -> Mat3<S>` is `½I + bW + dW²`. It is assembled as
     `SO3::jr` is: `W²` through the structured product, not `φφᵀ − θ²I`.
   - `gamma_apply_jacobian<const M: usize, S: Real>(phi, v: Vec3<S>) -> (Vec3<S>, Mat3<S>)`
     returns `(Γ_M(φ) v, ∂(Γ_M(φ) v)/∂φ)`.
     - It evaluates GG.2(d)'s `v/M! + σ φ×v + σ' φ×(φ×v)` once on `Dual<S, 3>` seeded at `φ`. The
       coefficients come from `jr_coeffs` (`M = 1`) or `gamma2_coeffs` (`M = 2`), selected by the
       const `M`.
     - `M ∈ {1, 2}` is asserted at monomorphization. `M = 3` has coefficients but loses `2u/θ²` on
       its exact arm (GG.3(c)), and `M ≥ 4` needs `σ₆`, which is outside `NUMERICS.md` §4.
     - The value is the vector form, a different rounding of `gamma_M(φ) · v`.
3. **The twin** (`NUMERICS.md` §14) is `reference::gamma_apply(m, φ, v)`, the first
   `GAMMA_TERMS = 40` terms of `Σ Wⁿ v/(n+m)!` summed left to right, with `Wⁿ v` as repeated cross
   products.
   - Its domain is `θ ≤ π`, where the last term is below `10⁻²⁷` of `v`.
   - Its proptests (`gamma_apply_jacobian_matches_reference_m{1,2}_{f64,f32}`) bound value and
     Jacobian per `θ` band at twice the table above.
4. **`NUMERICS.md`.**
   - §7's sentence becomes GG.6's: the `Dual` Jacobian is the derivative of the computed function,
     finite at `φ = 0`, and as accurate as the arm that runs.
   - §12 states that `Γ` has no row: it is entire, and its only limit is the overflow of `θ²` that
     `Exp` has (CO.15(c)).
5. **Corpus id `so3_gamma2`.**
   - Strata: `so3_jl`'s and each one's `@f32` twin, the binary64 draw rounded to binary32 per
     component. This extends `0016` to this one vector id. Every other vector id stays at
     binary64.
   - Reference: `Σ Wⁿ/(n+2)!` at 120 digits.
   - Scoring: output `G` (column-major), scored as a Jacobian is (floor 1).
   - Answered by `helicoid` at both precisions and by `helicoid:host-std` at binary64.
   - It moves from `coverage.rs`'s owed ids to `required()` with its corpus file.
6. **Seeded defect** (`PHASE1.md` §10): `Γ₂` built from `Γ₁`'s coefficients, `½I + aW + bW²` (the
   off-by-one in `m` that GG.2(a)'s `σ_{m+1}`, `σ_{m+2}` makes easy). It is planted as
   `gamma1(φ) − ½I`.
   - Its error is `≈ θ/3`, so it is required to read more than `1024×` the clean maximum on every
     stratum whose smallest `θ` is at least `12 288 u` (`θ/3 = 4096 u`), at both precisions.
   - Below that, a metric with a floor of 1 cannot see it: `theta:exact0`, where it is exact, and
     `theta:subnormal`.

## Rationale

- **No new switch.** "`gamma2_coeffs` through the sweep" is satisfied by the swept switches of its
  members, because the sweep is per coefficient. A group-level switch would be a second number
  for a coefficient that already has one, and `0004` forbids typing it.
- **The vector form under `Dual`, not the matrix.**
  - It is what GG.6 measured, and two cross products cost less than assembling `Γ` on `Dual` lanes.
  - The closed form `𝒥₁ = Q − [J_l ρ]_× J_l` (GG.5(c)) needs two groups and a product, and stays
    gated on a consumer bench (`PROJECT.md` §5.1).
  - `𝒥₂`'s closed form needs `τ₆`, outside §4.
- **`@f32` strata for this id alone.** "Both precisions" was asked for, and `Γ₂` at binary32 is a
  real consumer's case (an MCU attitude-and-heading filter). The `@f32` inputs are exact binary32s,
  as `0016` requires, so the score is exact.
- **A detection threshold, not "every stratum".** The planted error shrinks with `θ`. Requiring
  every stratum would require the metric to see `10⁻¹²/3` against `u = 6e-8`, which no metric
  with a floor of 1 can.

## Consequences

- `helicoid::so3` is a public path. Its module documentation is now public.
- `b` and `d` now have three call-site groups between them. The kernel's
  `groups_are_their_members_own_arm` and `gamma2_coeffs_are_q_coeffs_to_the_bit` hold them
  together.
- `so3_gamma2` is the first vector id an `f32` conformance run scores.

## Implementation plan

1. This record, with the `NUMERICS.md` §7, §12, §14, `PHASE5.md` §4, `API.md` §2 and
   `docs/maths/index.md` edits. Verified by `just lint` and `just doc`.
2. `gamma2_coeffs`, `so3::{gamma1, gamma2, gamma_apply_jacobian}`, `reference::gamma_apply` and
   `gamma_tests.rs`. Verified by `just test`, `just thresholds-check` (unchanged bytes),
   `just no-std` and `just msrv`.
3. `so3_gamma2`: generator, corpus, subject, metric, coverage, `PHASE1.md` §4.3 and §10, and
   `selftest_gamma`. Verified by `just corpus-check`, `just conformance` at both precisions,
   `just conformance --self-test` and `just envelope`.

## Open questions

None.

## Further work

- The closed-form `𝒥₁` (GG.5(c)), if a consumer's bench shows the `Dual` path is its bottleneck
  (`PROJECT.md` §5.1).
- `Γ₃`, if a consumer needs it; it owes the accuracy of `σ₅ = (d − 2e)/3` on its exact arm
  (GG.3(c)).
