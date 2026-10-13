# 0079: `from_matrix` normalizes a column of `K + I`

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record's PR (#153)

## Context

`SO3::from_matrix` is Shepperd's method (`NUMERICS.md` §3.4). It takes `d = √c_j` of the largest
pivot candidate, divides the three other numerators by `2d`, then normalizes. That is two `sqrt`
and seven divisions. locus-tag's tag-pose path converts every chosen IPPE branch through it, and
measured it at about 2× nalgebra's Shepperd, which does not normalize (2026-10-12, downstream
validation on locus-tag's `feat/helicoid-pose` branch).

`docs/maths/so3.md` SO.12(b) gives a cheaper form of the same quaternion: column `j` of
`K(R) + I` is `4 q̂_j q̂`, so normalizing that column returns `q̂` with its pivot component
positive, exactly the current output's sign rule. The column's entries are already formed (the
pivot candidate `c_j` and the six sums and differences), so the extraction needs no division
before the normalization: **one `sqrt`, four divisions**. The selected column's pivot entry is at
least `1` for every input (SO.11(c)), so the norm is at least `1` and the `sqrt` needs no safe
argument.

## Measured

Protocol:
- **Accuracy:** every record of `so3_from_matrix.jsonl` (2594 records, 30 strata, read
  column-major as the harness reads it), through both kernels at `f64` and on the matrix rounded
  to `f32`. Each output is scored exactly in Python against Shepperd's formula evaluated in exact
  rational arithmetic on the same input matrix (pivot by exact comparison, 80-digit `sqrt`). That
  isolates each kernel's own rounding. The error is the componentwise max of `|q − s q_exact|/u`,
  with `s = ±1` the better sign. The id itself is unscored: the harness reports matrix inputs as
  backward error only (`NUMERICS.md` §11), so `just conformance` cannot adjudicate this change.
- **Speed:** `cargo xtask bench-gate --bench groups --only so3/from_matrix --against <base>`, the
  base binary built clean from `main` with this PR's two bench rows added: one per stratum, and
  `mixed-64`, 64 rotations whose pivot changes unpredictably from call to call. `taskset -c 2`,
  EPYC-Milan. Ratios are candidate/base, over two runs.

| | `f64` | `f32` |
|---|---|---|
| Time, three strata (identity, generic, near-`π`) | **0.54–0.60** (floor 0.3–0.8%) | **0.70** (floor 0.3–0.6%) |
| Time, `mixed-64` | **0.57** | **0.70** |
| Max error over all records, old → new | 1.557 → **1.507** `u` | 2.013 → **1.465** `u` |
| Mean error, old → new | 0.362 → 0.416 `u` | 0.329 → 0.355 `u` |
| Records new better / worse | 359 / 542 | 284 / 396 |
| Strata whose max grows | 11 of 30, by at most +0.32 `u` (`theta:pi-1e-7`, 1.150 → 1.472) | 4 of 30, by at most +0.26 `u` (`theta:1e0`, 1.208 → 1.465) |

The mean grows because the old denominator `2√c_j` is one rounded `sqrt` of one number, and its
final normalization rescales a vector already unit to `O(u)`. The new denominator
`√(Σ col²)` carries the rounding of four squares and three sums. The worst case shrinks because the
old form's error compounds through two divisions where the new one has one.

## Decision

1. **`NUMERICS.md` §3.4:** the extraction is "column `j` of `K(R) + I` (SO.12(b)), normalized",
   with the trace pivot's column `(1 + tr R, R₂₁ − R₁₂, R₀₂ − R₂₀, R₁₀ − R₀₁)` written out. The
   pivot rule, the nested `branch`es and the sign rule are unchanged.
2. **`SO3::from_matrix`** forms the selected column inside the branches and normalizes once
   outside them, with the four divisions of `Quat::from_wxyz_normalized`.
3. **The arms return scalars, captured by `move`.** Arms that borrow their candidates let LLVM
   merge the four arms' reads into one load through a selected stack address: every candidate is
   spilled, and the pivot compare feeds the load's address through `cmov` instead of a predicted
   branch. That kernel, the PR's first, read 0.88 (`f64`) and 1.05–1.10 (`f32`) of the old time.
4. **An overflowed norm gives NaN.** `n = sqrt(n2) + n2 · 0`: the added term is `+0` for a
   finite `n2`, so `n` and every in-domain output are unchanged to the bit, and NaN once `n2`
   overflows. That is an entry beyond about `2^511` (`f64`) or `2^63` (`f32`), where `q / ‖q‖` would
   otherwise be a finite zero quaternion. It costs one multiplication beside the `sqrt` and one
   addition after it. Rescaling instead, through a cold branch, kept the old kernel's range and
   cost the gain: 0.98 (`f64`) and 1.18–1.21 (`f32`).

## Rationale

The alternative that keeps the old accuracy and drops the second `sqrt` is to skip the
normalization: `q = col/(2√c_j)`. It breaks §3.4's "then normalize" and the rustdoc's promise of a
unit quaternion for a non-orthogonal input (`|η|` up to `6ε` by SO.12(d)), and `SO3`'s invariants
assume `|η|` near `u`. Rejected.

## Consequences

- Output bits change on about half the records (1218 of 2594 at `f64`), each within the bounds
  above. No other function's bits change.
- Out of domain, the range shrinks. The old kernel returned a finite unit quaternion for entries
  up to about `2^1022`; this one returns NaN beyond `2^511` (`f64`) or `2^63` (`f32`).
  `from_matrix_is_nan_where_the_column_norm_overflows` pins it. A rotation's entries are at most `1`.
- Every matrix with finite entries below that overflow gives a unit quaternion, the zero matrix
  included: the selected column's pivot entry is at least `1`.
- Branch-free pivot selection measured slower at every row, `mixed-64` included: blending whole
  columns costs more than the branch it replaces, mispredictions and all. With the rescaling
  guard it read 1.42 (`f64`) and 1.78 (`f32`) on the strata, and 1.13 and
  1.44 on `mixed-64`.
- The no-regress bar (`0006`) is per stratum, on the max. If `so3_from_matrix` were scored,
  11 `f64` strata would regress, by up to 0.32 `u`, while its overall max improves. The
  bar does not read this id today. The owner accepted that trade for the speed (2026-10-13).

## Implementation plan

1. One PR: the kernel, the §3.4 edit, the `so3/from_matrix` bench rows, the CHANGELOG. Verified by
   `just test`, the measurement above and the bench gate above.

## Open questions

None.

## Further work

- A forward-error score for matrix-input ids against the exact-arithmetic Shepperd, as this
  record's protocol computes it, so the instrument can read `from_matrix` changes itself.
