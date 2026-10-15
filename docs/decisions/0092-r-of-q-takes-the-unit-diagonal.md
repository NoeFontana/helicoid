# 0092: `R(q)` takes the unit diagonal

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

[`0076`](./0076-the-oracles-first-verdict-r-of-q-stays-and-drift-is-excepted.md) decision 1 kept
`Quat::to_matrix`'s homogeneous diagonal, `(w² − ‖u‖²) + 2uᵢ²`, over Eigen's unit diagonal
`1 − 2(uⱼ² + uₖ²)`:
- **Its reason.** On a rounded unit `q`, `‖q‖² = 1 + η` with `|η| ≈ u`, and every reference
  evaluates `R(q/‖q‖)`. The homogeneous form errs by `ηR` at every angle, the unit form by
  `η(R − I)`: nothing at the identity, `2√2|η|` at `π`.
- **Its measurement.** Under the strict bar the unit form was better on 22 strata of each `sen3_ad`
  id and worse on 20. Each form won one end of `θ`.

Two things have changed since:
- **The bar resolves a loss.** [`0088`](./0088-the-domination-bar-must-tell-worse-from-equal.md)
  fails a stratum only where the swap test resolves it.
- **The deep tier has run end to end.** On its 64-draw strata, judged against the best oracle with
  `q`/`−q` twins as one draw, `Ad`'s `R` block is the largest family of resolved losses, 67 of
  ~100. All sit near the identity, at `θ` from `1e-8` to `1e-1`, where the homogeneous form reads
  1.1–1.9 `u` and GTSAM, `tf_tree_math`, manif and Sophus all read 0.5–0.6. 15 lose to
  `tf_tree_math`, on the `libm` crate, so none is D16's.
- **Near `π` the unit form loses nothing.** It reads 2.7–3.4 `u` against the homogeneous 1.9–2.1,
  and GTSAM, the best oracle there, *is* the unit form, so it ties. A tie is no failure.

## Decision

1. **`Quat::to_matrix` evaluates `R(q) = (1 − 2‖u‖²) I + 2uuᵀ + 2w[u]×`**: each diagonal entry
   `1 − (2uⱼ² + 2uₖ²)` as written, the off-diagonal entries and the nine products unchanged
   (`NUMERICS.md` §1).
   - **Unit `q`:** the same matrix.
   - **Any other `q`:** `‖q‖² R(q/‖q‖) − (‖q‖² − 1) I`, which is `(1 + η)R − ηI`.
2. **`NUMERICS.md` §12's carried contract for `to_matrix` becomes `(1 + η)R − ηI`.** Applied to
   `v` that is `Rv + η(Rv − v)`, `act`'s own term. Today `to_matrix` and `act` disagree on a drifted
   quaternion; after this they agree.
3. **The twin of [`0078`](./0078-a-comparison-drops-ties-and-prints-both-maxima.md) is inverted.**
   `helicoid:unit-r` becomes `helicoid:homogeneous-r`, the form this record replaces, so that
   `just compare --a helicoid:homogeneous-r --b helicoid` re-measures the choice. It reproduces
   `0076`'s table to the stratum (`the_homogeneous_r_twin_reproduces_the_tally_of_0076`).

## Rationale

The deep tier, `sen3_ad_n1..3` and `sim3_ad` at both precisions, with twins as one draw and
`0088`'s `α = 0.001`:

| `R(q)` diagonal | resolved losses | strata above the best oracle | summed maxima (`u`) | `to_matrix` latency |
|---|---|---|---|---|
| homogeneous (main) | 67 | 136 | 1272 | 2.59 ns `f64`, 3.57 `f32` |
| **unit (this)** | **5** | 23 | 1846 | 2.77 (+7%), 2.77 (−22%) |
| unit near the identity, homogeneous past `\|w\| = ½` | 5 | 20 | 1464 | 3.88 (+50%), 4.06 (+14%) |
| homogeneous divided by `‖q‖²` | 74 | 154 | 552 | not timed |
| first-order normalized, `t + tη` per entry | 4 | 35 | 422 | 5.33 (+106%), 5.81 (+63%) |

Latency is a consumer loop: 4 096 rotations per band, best of 40, one pinned core. The SE(3)
adjoint goes from 5.13 to 5.00 ns at `f64` and from 6.51 to 7.35 at `f32`.

- **The unit form takes 62 of the 67.** The 5 left are `sim3_ad`, whose excess is the scale.
- **Summed maxima rise from 1272 to 1846 `u`, all in `q:nonunit`.** That stratum is deliberately
  far from unit: 257 → 434 `u` at `sen3_ad_n1`, which is already excepted (`0076` decision 2) and
  stays a gross loss either way.
- **The select by `|w|` keeps the `π` lead, but not at this cost.** LLVM turns two of its three
  selects into branches, and `to_matrix` is in every consumer's hot path.
- **Normalizing fixes `q:nonunit`** (257 → 1.6 `u`), which `0058` declined for carried
  quaternions. At twice the latency it is not worth it, for a stratum the exception table already
  explains.

## Consequences

- **Bits.** `cargo xtask determinism` against `main` moves exactly the `Ad` digests:
  `sen3_ad_n1..3` at `f64` and `sim3_ad` at both precisions. No other id the corpus scores reads
  `R(q)` as a matrix.
- **`0076` decision 1 is superseded** when this record is `ready`; its decision 2 (the two
  exception rows) stands. `Quat::to_matrix`'s rustdoc now states this record, not `0076`.
- **A public behaviour change for a non-unit `q`.** `Quat::to_matrix` is no longer a scaled
  rotation there: `q = 2i` gives `diag(1, −7, −7)`, where it gave `4·diag(1, −1, −1)`. A consumer
  that recovers `R` as `to_matrix() / ‖q‖²` must normalize `q` instead. On `0.0.x` the release
  notes say so.
- **`to_matrix` and `act` are one map.** The unit form is exactly the matrix of §3.3's vector form,
  since `[u]×² = uuᵀ − ‖u‖²I`. `docs/maths/so3.md` SO.2(d) and SO.14(c) say so, and `act`'s and
  `to_matrix`'s *Domain* give the shared error, `η(Rv − v)`.
- **The `f32` `Ad` cost is codegen.** The disassembly is branch-free with no stack traffic either
  way. LLVM's vectorizer packs the new diagonal with 45 `shufps` against 33, and regrouping the
  sums in the source does not change its choice.
- **`quat_tests::matrix_hand_cases`** pins the non-unit contract: `2·1 ↦ I`, `2i ↦ diag(1, −7, −7)`.
  **`matrix_sums_group_as_written`** pins the new grouping on a quaternion where `(1 − 2y²) − 2z²` and
  `(1 − 2z²) − 2y²` each move all three diagonal entries by an ulp.

## Implementation plan

1. `Quat::to_matrix`, `NUMERICS.md` §1 and §12, the tests above, the drift-band test reading
   `act`'s term from `to_matrix`, and the inverted twin. Verified by:
   - `just test`, `just lint`;
   - `cargo xtask determinism` against `main`;
   - `cargo xtask envelope` with and without `--deep`.

## Open questions

1. **Accepting `π`'s tie for the identity's win.** Near `π` the unit form gives back up to 1.3 `u`
   and ties GTSAM. That is the owner's call.
