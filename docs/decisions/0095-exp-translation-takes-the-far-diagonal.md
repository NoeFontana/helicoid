# 0095: `Exp`'s translation takes `J_r`'s far diagonal

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

`SEn3::exp` applies `J_l(φ)` to each `ρᵢ` by two cross products (`NUMERICS.md` §5.1):
`ρ + a (φ × ρ) + b (φ × (φ × ρ))`.
- **The cancellation.** `φ × (φ × ρ) = (φ·ρ) φ − θ² ρ`, and `bθ² → 1` at `π`. So the `ρ` part of
  the column is `(1 − bθ²) ρ = (sin θ/θ) ρ`, formed out of two terms of opposite sign. It is
  `J_r`'s diagonal, which [`0089`](./0089-the-closed-arms-take-theta-squared-as-given.md)
  decision 3 rewrote for the matrix.
- **The standing.** `sen3_exp` is the SE_N(3) family furthest behind the oracles. On `main` it is
  at or below the best oracle on 100 of 156 deep strata. Its worst deep stratum is
  `sen3_exp_n1/theta:pi-1e-2`, at 3.44 `u` against manif's 1.88.
- **The quaternion is not the target.** The maxima that sit at the quaternion's floor near `π`
  (`sen3_exp_n*/theta:pi-1e-9`, 2.743 `u` at every `N`) are `SO3::exp`'s. A cheap rewrite of it
  measured mixed, and the uniform one costs about 50% of `SO3::exp`.

## Decision

1. **Where `b` is on its exact arm, the column is `(sin θ/θ) ρ + a (φ × ρ) + b (φ·ρ) φ`,** with the
   `sin θ/θ` that arm formed `b` from (`jr_coeffs_sinc`, `0089`). Below the switch it stays
   `ρ + a (φ × ρ) + b (φ × (φ × ρ))`, to the bit.
2. **One `branch` on `b`'s own mask.** It is the predicate `jr_coeffs_sinc`'s group has already
   branched on, so there is no new switch. The quaternion is formed before the branch.
3. **`NUMERICS.md` §5.1 states both forms.**

## Rationale

**Accuracy.** Every subject was rescored on both tiers, against `main` at `3b031d2`.
`cargo xtask determinism` moves exactly `sen3_exp_n1..3` and `se3_screw_retract` at `f64`.

| | `main` | this |
|---|---|---|
| `sen3_exp`, deep: strata better / worse | | 38 / 1 |
| `sen3_exp`, deep: worst stratum | 3.64 `u` | 2.79 |
| `sen3_exp`, deep: strata above the best oracle | 46 | 31 |
| `sen3_exp`, committed: strata better / worse | | 22 / 10 |
| deep tier: domination failures | 2 | **1** |
| deep tier: unresolved / unresolvable | 65 / 151 | 62 / 144 |
| committed tier: failures / unresolved / unresolvable | 0 / 23 / 505 | 0 / 23 / 501 |

- **The failure removed** is `sen3_exp_n2/theta:pi-1e-9`. Its maximum is still the quaternion's
  2.743 `u`, but the stratum no longer resolves against GTSAM.
- **Per record near `π`:** on `sen3_exp_n1` the mean falls in all 12 `theta:pi-1e-k` strata, for
  example from 1.42 to 1.11 `u` at `pi-1e-6` and from 1.41 to 1.14 at `pi-1e-12`.
- **`se3_screw_retract`** composes this `Exp`.
  - **Committed tier:** 9 strata better and 2 worse.
  - **Deep tier:** the maximum moves both ways (4 better, 8 worse; worst 4.44 → 5.85 `u` at
    `theta:pi-1e-1`).
  - **Per record:** more records improve than worsen in every one of the 12 near-`π` strata, and
    the mean falls in 11.
  - **On the bar:** none of its strata fails.

**Latency.** `SEn3::exp` in a consumer loop, 4096 tangents per band, best of 40, one pinned core,
`ρᵢ ~ U[−2, 2]`, `main` → this (ns):

| | `θ ~ 10⁻⁶..10⁻²` | `0.5..2` | `2.8..3` | `3.11..3.14` |
|---|---|---|---|---|
| `f64`, `N = 1` | 6.55 → 6.33 | 13.2 → 12.9 | 32.5 → 33.1 | 29.8 → 30.7 |
| `f64`, `N = 3` | 10.2 → 11.1 | 17.3 → 18.2 | 35.7 → 37.6 | 33.4 → 35.3 |
| `f32`, `N = 1` | 4.77 → 5.26 | 9.06 → 9.50 | 14.5 → 14.8 | 27.7 → 19.4 |
| `f32`, `N = 3` | 8.86 → 8.55 | 13.5 → 12.6 | 18.4 → 17.9 | 29.5 → 24.2 |

The arithmetic is the same size on both arms. The differences are code generation, within about
±10% except at `f32` near `π`, where the far arm is faster.

**Rejected:**
- **One formula everywhere, `(1 − bz) ρ + a (φ × ρ) + b (φ·ρ) φ`.** Below the switch it reads up
  to 1.9× `main` (`theta:1e-8`: 0.93 → 1.81 `u`); 65 deep strata were worse.
- **Both columns formed and selected per component.** Same bits as this record, but +4.2 ns at
  `N = 3` near the identity, where most calls are.
- **Per-arm scalars into one five-term formula** (`(1, a, b, 0)` and `(sin θ/θ, a, 0, b)`). Same
  bits and no branch, but +1.5 and +3.7 ns near the identity.
- **The branch with the quaternion formed after it.** Same bits, and 4–7 ns slower at `f64` on
  the far arm.
- **The matrix `J_l ρ`, GTSAM's.** `0036` measured it worse away from `π`.

## Consequences

- **`tf_tree_math::exp_se3` parity** holds below `b`'s switch (`θ² < 7.77` at `f64`) and not above
  it, where `helicoid` is the more accurate.
- **`Exp`'s quaternion near `π` is the remaining loss** (`sim3_exp/sigma:-3/theta:pi-1e-6@f32`,
  the deep tier's one failure). It is `SO3::exp`'s, not this record's.

## Implementation plan

1. `SEn3::exp`; `NUMERICS.md` §5.1. Verified by:
   - `just test`, including `exp_translation_is_jl_times_rho_on_both_arms` (both arms against
     `SO3::jl(φ) ρ`, through the switch and up to `π`) and the group laws;
   - `just lint`, `just doc`;
   - `cargo xtask determinism` against `main`;
   - `cargo xtask envelope` with and without `--deep`, every subject's rows at one rev.

## Open questions

None of substance; it is a `draft` until the owner accepts the new evaluation.

## Further work

- **`Sim3::exp`'s `V`** has the same `χ = v₀ − zv₂` cancellation near `σ = 0`, `θ → π`.
- **One `sin(θ/2)` for `SEn3::exp`.** `exp_coeffs` and `jr_coeffs_sinc` each evaluate `sin(θ/2)` on
  their exact arms. A group holding `k`, `cos(θ/2)`, `a` and `b` would evaluate it once with the
  same bits, which is one `libm` call of the far arm's three.
