# 0086: `from_matrix` keeps Shepperd's quotients on the sphere

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record's PR (#162)

## Context

locus-sensors builds its frame conventions (FLU, FRD, RDF, RUB, ENU, NED) with `SO3::from_matrix`
on the 24 signed permutation matrices of determinant `+1`. They are exact rotations, and their
quaternion components lie in `{0, ±½, ±√½, ±1}`. Since `0079`, `from_matrix` normalizes the pivot's
column of `K + I`. On a quarter turn that column is `(2, 0, 0, 2)`, its norm is `fl(√8) =
2 fl(√2)`, and both nonzero components round to `fl(1/fl(√2))`. That is one ulp below `fl(√½)`, so
`‖q‖² − 1 = −1.6 u` exactly and `fl(‖q‖²) = 1 − 2u`. `SO3::act` normalizes nothing and carries
`2|η|‖v‖` (`0058`), so every later `act` inherits the deficit. Two equal components can never give
`fl(‖q‖²) = 1` on a quarter turn, so only a formula that rounds them apart can.

[`0085`](./0085-from-matrix-is-measured-on-exact-rotations.md) scores `from_matrix` by the angle to
the polar factor and the drift off the sphere, on a `matrix:exact` stratum and at both precisions.
Its first run put `0079` above the best oracle on 18 of 31 `f64` strata. `tf_tree_math`, whose
`quat_from_rot3` is Shepperd's quotients with no normalization on the same `libm`, reads 0.49 `u`
on the `θ` decades where `0079` reads 0.9 to 1.4.
[`0088`](./0088-the-domination-bar-must-tell-worse-from-equal.md) (draft) resolves 11 of those 18
as real losses and attributes all 11 to the program, not to D16.

## Decision

1. **Shepperd's quotients, in `tf_tree_math`'s arithmetic.** `NUMERICS.md` §3.4 now reads: pick
   the largest `c_j` as before. Each axis candidate adds its own diagonal first, `c_1 = ((1 + R00)
   − R11) − R22` and its two rotations. Form `q_j = ½√c_j` and the other three as the column's
   entries over `2√c_j`.
2. **Kept within `ε = 16u` of the sphere.** When `|fl(‖q‖²) − 1| ≤ 16u` the quotients are returned.
   Otherwise the matrix is scaled, non-orthogonal or NaN, and `0079`'s normalized column is returned,
   its NaN guard included. A matrix off SO(3) whose quotients land within `16u` keeps them, so its
   output is within `16u` of the sphere rather than `0079`'s `3.5u`.
3. **Tests.**
   - `from_matrix_keeps_exact_rotations_on_the_sphere`: all 24 rotations at both precisions,
     `fl(‖q‖²) = 1` exactly at `f64` and within `u` at `f32`, and `act` on each axis returns the
     column to `4u`.
   - `from_matrix_rounds_a_quarter_turn_both_ways`: the bits.
   - `from_matrix_keeps_the_quotients_on_rotations`: 20 000 rotations per precision, uniform, near
     `π` and near the identity, equal to the bare quotients to the bit, so none falls back.
   - `from_matrix_falls_back_on_every_pivot`: a rotation scaled by `1 + 2^-16`, one per pivot,
     comes back unit and denotes the rotation, so the fallback's column tree stays the quotients'.
     Swapping one arm's `sx` and `sy` fails it.

## Rationale

**Why `16u`.** Rotations from `SO3::exp` then `to_matrix`, rounded to each precision, 10⁷ per
family (uniform `θ`, `π − 10^-k` for `k` in 1, 3, 6, 9, 12, and `θ ~ 10^-k`), never put the
computed `|fl(‖q‖²) − 1|` of the quotients above `8u`, at either precision. The first draft's `2u`
sent 2.5% of uniform rotations and 8 to 15% of those near `π` to the column. `16u` is twice the
measured ceiling. It is a test of whether the input is a rotation, not a series switch, so `0004`
does not apply.

**Under `0085`'s instrument**, from `cargo xtask conformance --fn so3_from_matrix` at both
precisions. "Above the best oracle" is a stratum's maximum against the smallest of `gtsam`,
`sophus` and `tf_tree_math` (`sophus` alone at `f32`):

| | `0079` | first draft (`2u`, old `c_j` order) | this |
|---|---|---|---|
| `f64` strata above the best oracle | 18 of 31 | 13 of 31 | **8 of 31** |
| `f32` strata above the best oracle | 14 of 31 | 11 of 31 | **9 of 31** |
| `f64` max, `q:nonunit` aside | 1.567 | 1.659 | 1.714 |
| `f32` max, `q:nonunit` aside | 1.834 | 1.544 | 1.525 |
| `matrix:exact`, `f64` / `f32` | 0.798 / 0.287 | 0.707 / 0.287 | 0.707 / 0.287 |
| `q:w0`, `f64` | 1.289 | 1.659 | 0.593 |
| `q:nonunit`, `f64` | 89.41 | 89.41 | 89.41 |

- **At `f64` this is `tf_tree_math`** on every stratum's maximum but `q:nonunit`, where the fallback
  reads 89.41 `u` against its 162.99.
- **Against `0079`,** 6 `f64` and 5 `f32` strata near `π` read higher, by up to 0.50 `u`
  (`theta:pi-1e-1`, 1.218 → 1.714). On those strata the oracles disagree with each other by as much:
  at `theta:pi-1e-1` `tf_tree_math` reads 1.714 and Sophus and GTSAM 1.442, and at `theta:pi-1e-2`
  Sophus and GTSAM read 2.161. Near `π` every candidate here and every oracle takes the same axis
  pivot, so their exact Shepperd values coincide and EA.9(d)'s definitional gap to the polar
  factor (up to `7/4 ε` per component) is common to all; what differs there is rounding. Where
  the trace pivot and an axis pivot compete, the gap is not common, and a stratum maximum ranks
  conventions as well as rounding. `0079` has no committed baseline under `0085`'s metric, so no
  no-regress bar applies. Whether any of them is a real loss is for `0088` (draft) to say, at its
  64-record tier.
- **Rejected candidates:**
  - **The first draft's `c_j` order,** `((1 − R00) + R11) − R22` and `((1 − R00) − R11) + R22`:
    13 strata above the oracles, and `q:w0` at 1.659.
  - **`2u`:** as above. At `16u` every stratum's maximum equals the first draft's, so neither
    threshold sends a corpus record to the column. The difference is the in-domain inputs the
    corpus does not hold.
  - **One `sqrt` ahead of the pivot tree,** the tree returning both the quotients and the column.
    It is bit-identical (114 of 114 `cargo xtask determinism` digests) and has no stack spill, but
    the gate reads `f32` `mixed-64` at 2.77× 0.0.3 against 1.62×, and `f64` at 1.12–1.28× against
    1.03–1.08×. In the branch's arms a predicted pivot starts its `sqrt` at once; hoisted, the
    `sqrt` waits on the comparisons.
  - **The pivot written as selects,** all four divisions and the pivot's slot replaced: 2.1× (`f32`)
    and 1.55× (`f64`).
  - **One reciprocal and four multiplies:** the same bits on exact rotations, more strata above
    the oracles.
  - **The column first,** swapping to the quotients only when its own drift reaches `2u`: still 10
    `f64` strata above `0079` near `π`.
  - **The bare quotients,** with no fallback, are `tf_tree_math`: 162.99 `u` at `q:nonunit` and
    `1e9 u` on a matrix perturbed by `1e-6`.
  - **A Newton step** on the quotients (corpus `f64` maximum 1.759, no convergence at `1e-6`) or on
    the column (its factor rounds to `1` on a quarter turn).

## Consequences

- **Bits.** Of the 114 `cargo xtask determinism` digests, only `so3_from_matrix`'s change.
- **tf_tree.** With its `helicoid` feature, `tf_tree_math::quat_from_rot3` delegates here (its
  `0063`). Under `0079` it was up to 2.85× less accurate than its own native body on the `θ`
  decades. Under this record the two share their arithmetic wherever they pick the same pivot.
  The native body takes the trace pivot whenever `tr R > 0`, and this the largest candidate.
- **`Dual`.** The derivative is the arm's. Along SO(3) both arms agree. Off it, the quotients'
  derivative has a component normal to the sphere, which `0079`'s normalized column did not. An
  ambient Jacobian of `from_matrix` (`PHASE6.md` §1) has to say which one it states.
- **Cost.** One `sqrt` and three divisions, branch-selected, where `0079` had one `sqrt` and four
  divisions that `f32` packs into one `divps`. On a quiet host the gate against 0.0.3 reads
  1.63–1.85× at `f32` (three of four rows resolved, floors under 1%) and 1.03× at `f64` (`near-pi`
  resolved; the other rows inside their floors). The `f32` cost is Shepperd's: `tf_tree_math`'s
  native body does the same divisions.

## Implementation plan

1. Decisions 1–3, `NUMERICS.md` §3.4 and the rustdoc. Verified by `just test`, by
   `cargo xtask conformance --fn so3_from_matrix` at both precisions against `0085`'s oracle rows,
   by `cargo xtask determinism` against `main`, and by `cargo xtask bench-gate --bench groups --only
   so3/from_matrix --against <0.0.3's bench binary>`.

## Open questions

1. **The `π` strata against `0079`.** They wait on the 64-record tier of `0088` (draft).

Answered by the owner, on merging this record's PR: the `f32` latency, 1.63–1.85× 0.0.3, is
accepted for 9 strata above the oracles instead of 14 and exact rotations on the sphere.

## Further work

- Correctly rounded components (scipy's tables) read 5.41 `u` on `act` at `f64` per the consumer's
  measurement. Shepperd's mixed rounding beats them, and nothing here pursues correct rounding.
