# 0097: The Jacobian's column is five words

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

`SEn3::jr` evaluates `NUMERICS.md` §5.3's `Q(−ρᵢ, −φ)` as written.
- **The cost.** Seven hat products per column, each `18` multiplications, and an eighth for `W²`.
  `jl` is the provided `J_r(−τ)`, and `rplus_jacobians` and `lplus_jacobians` call them.
- **The arrangement is one of many.** `docs/maths/se3.md` SE.8 remark (ii) reduces the same block to
  five words, `Q = aΞ + b(WΞ + ΞW) − (φ·ρ)[(b − 2d)W + 2eW²]` with `a = ½ − zd`. It did not adopt
  them: another arrangement is a §5.3 edit and a record.
- **The precedent.** [`0093`](./0093-the-inverse-jacobian-column-is-four-words.md) (draft) took the
  inverse's column to four words and split it into a symmetric part, shared by the two sides, and
  a skew part, whose sign they differ by. The forward column splits the same way:
  - `S(ρ) = b(WΞ + ΞW) − 2e(φ·ρ)W²` is unchanged when `φ` and `ρ` are both negated;
  - `[aρ − (φ·ρ)(b − 2d)φ]^` changes sign.
- **The standing.** `sen3_jr_n*` and `sen3_jl_n*` were above the best oracle on 49 of 312 committed
  strata and 30 of 312 deep ones. Their worst stratum read 7.05 `u` on the deep tier.

## Decision

1. **`J(τ)`'s column is evaluated in five words** (`NUMERICS.md` §5.3, SE.8 remark (ii)): `J_l`'s is
   `S(ρ) + h^` and `J_r`'s is `S(ρ) − h^`, with `h = aρ − (φ·ρ)(b − 2d)φ`.
   - **`S` is formed as `0093`'s is** (`Cols`, which serves both): each diagonal entry is
     `−(2b(φⱼρⱼ + φₖρₖ) − 2e(φ·ρ)(φⱼ² + φₖ²))`, so no `(φ·ρ)` or `θ²` that would cancel is formed.
2. **The coefficients come from `coeffs::q_five_coeffs`.** It returns `(a, b, b − 2d, e)` and
   `J_r(φ)`'s `(sin θ/θ, far)`, from `jr_coeffs_sinc` and `q_coeffs`, each member its own group's to
   the bit.
   - **`a` is SO(3)'s `(1 − cos θ)/θ²`** (CO.4(a)), not `½ − zd`.
   - **`b − 2d` is `(2a − sin θ/θ)/z` on `b`'s exact arm**, two terms of one sign near `π`, where
     `b − 2d` loses a factor of about 2.5. Below that arm it is `b − 2d`.
   - **There is no new switch:** the choice is `b`'s own mask, as `0089`'s `J_r` diagonal is.
3. **The diagonal is `SO3::jr`'s,** formed from the same call (`HatSq::poly_far`). `J_l`'s is its
   transpose, `SO3::jl` to the bit.
4. **`jl` is overridden** with the same code at the opposite sign of `h` and the transposed
   diagonal. It is the provided `J_r(−τ)` to the bit, and
   `the_fused_inverses_are_the_separate_ones_to_the_bit` compares the two.

## Rationale

**Accuracy.** Every subject was rescored on both tiers, against `main` at `f78086d`.
`cargo xtask determinism` moves exactly `sen3_jr_n1..3` and `sen3_jl_n1..3` at `f64`.

| tier | records better / worse | total error | near-`π` strata | strata maxima better / worse | worst stratum | above best oracle |
|---|---|---|---|---|---|---|
| committed | 900 / 562 | 0.929× | 0.876× | 147 / 67 | 4.08 → 2.68 `u` | 49 → 48 of 312 |
| deep | 9194 / 6333 | 0.932× | 0.880× | 113 / 103 | 7.05 → 4.96 `u` | 30 → 28 of 312 |

- **The coefficients' part.** With `a = ½ − zd` and `b − 2d` as differences, the deep tier read
  0.946× total error and 5.38 `u` worst. Decision 2 moved 68 deep strata better and 10 worse
  beyond that.
- **The envelope:** both tiers read as on `main`. The deep tier has 1 domination failure, `SO3::exp`'s
  quaternion near `π`, with unresolved 62 → 60 and unresolvable 144 → 141. The committed tier has
  0 failures and 23 unresolved, with unresolvable 501 → 500.
- **Near `π` the gap to the best oracle remains** (`theta:pi-1e-3`: 2.31 `u` against Sophus's 1.34),
  unresolved. It is a single record per stratum, and SO(3)'s own `J_r` there is mid-field among the
  oracles (2.36 against 1.67–3.17).

**Latency.** Consumer loop, 4096 tangents per band, best of 40, one pinned core, `ρᵢ ~ U[−2, 2]`,
`f64`, `main` → this (ns):

| call | `N` | `θ ~ 10⁻⁶..10⁻²` | `0.5..2` | `2.8..3` | `3.11..3.14` |
|---|---|---|---|---|---|
| `jr` | 1 | 40.3 → 16.4 | 53.3 → 27.1 | 87.2 → 50.8 | 81.7 → 47.1 |
| `jr` | 3 | 109.1 → 33.1 | 119.7 → 42.8 | 152.6 → 66.3 | 147.3 → 62.7 |
| `jl` | 3 | 112.7 → 28.2 | 122.8 → 39.3 | 155.9 → 65.2 | 150.5 → 61.9 |
| `rplus_jacobians` | 1 | 52.0 → 26.0 | 73.8 → 39.3 | 121.9 → 85.2 | 113.1 → 76.7 |
| `rplus_jacobians` | 3 | 149.2 → 66.1 | 163.4 → 80.3 | 217.9 → 131.9 | 209.1 → 125.1 |

At `f32`, `jr` at `N = 3` goes from 92.8, 97.6, 98.3 and 114.5 ns to 27.5, 35.1, 35.6 and 51.3.

**Rejected:**
- **§5.3's seven words, as written.** They are slower, and less accurate by the table above.
- **The five words with `a = ½ − zd` and `b − 2d` formed in `sen3.rs`.** That is the coefficients'
  part above, and it would evaluate a cancelling coefficient outside `helicoid::coeffs`.
- **`Q` as a matrix product, GTSAM's `H` of `applyLeft`.** It cannot be cheaper than five words,
  for `0093`'s reason.

## Consequences

- **`q_block` is removed.** `jr` and `jl` go through `forward`, and `0093`'s `InvCols` becomes `Cols`,
  which serves both directions. `inverse`'s and `inverses`' bits do not move.
- **`q_coeffs` is read by the library only through `q_five_coeffs`.** Its re-export is kept for the
  sweep and the tests.
- **SE.8 remark (ii) is the evaluation and §5.3's nine-term `Q` the definition,** as SE.17 and SE.12
  are for the inverse.

## Implementation plan

1. `coeffs::q_five_coeffs`; `Cols`; `SEn3::{jr, jl}` through `forward`; `NUMERICS.md` §5.3;
   `docs/maths/se3.md` SE.8 remark (ii). Verified by:
   - `just test`, including `q_five_coeffs_forms_b_minus_2d_without_the_cancellation`,
     `the_fused_inverses_are_the_separate_ones_to_the_bit` (now also `jl` against `J_r(−τ)`, 8004
     blocks) and the group laws' `jacobians_match_dual`;
   - `just lint`, `just doc`;
   - `cargo xtask determinism` against `main`;
   - `cargo xtask envelope` with and without `--deep`.

## Open questions

None of substance; it is a `draft` until the owner accepts the new evaluation.

## Further work

- **The near-`π` maxima of `sen3_jr`/`jl`** (1.7× the best oracle on one record per stratum, at
  `pi-1e-3`): not the columns' coefficients, and not resolved.
