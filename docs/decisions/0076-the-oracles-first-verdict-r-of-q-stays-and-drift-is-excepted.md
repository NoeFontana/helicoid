# 0076: the oracles' first verdict: `R(q)` stays, and carried drift is excepted

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #145

## Context

The three container oracles ([`0072`](./0072-a-container-oracle-is-an-image-on-the-same-protocol.md),
[`0074`](./0074-the-gtsam-oracle-is-built-from-source-under-our-flags.md),
[`0075`](./0075-the-manif-oracle-reads-se23-as-position-then-velocity.md)) take `just envelope` from
148 domination failures to 709, 38 of them D16's. Most are gaps of 1–3× on maxima of 4–6
records, which [`0037`](./0037-six-records-cannot-answer-a-one-variable-question.md) (draft)
measured as unresolved at that count and [`0038`](./0038-a-program-comparison-is-not-a-bar.md)
keeps on the max: they stay failures, and none of them is a formula to chase. Two findings are
not of that kind, and each asks for a decision.

1. **`Ad`'s `R` block.** Every agent that replayed a large `Ad` gap — `sim3_ad/sigma:-3/theta:0.3`
   (1.350 against Sophus's 0.426 `u`), `sen3_ad_n2/theta:1e-2` (1.565 against manif's 0.445) —
   found it in `Quat::to_matrix`, and Eigen's `1 - 2(y² + z²)` diagonal, which assumes `‖q‖ = 1`,
   reproducing the oracle's number. `sen3_ad_n1/theta:1e-2` already lost 3.3× to `tf_tree_math`,
   which is on the `libm` crate, so it is not D16.
2. **`q:nonunit`.** `so3_act` reads 510.3 `u` against Sophus's 7.867 and `sen3_ad_n1` 256.1
   against 2.132. Sophus normalizes `q` when it builds the element; `helicoid` carries it
   ([`0058`](./0058-a-drifted-quaternion-is-carried-not-vouched-for.md)).

## Decision

1. **`R(q)` stays `(w² - ‖u‖²) I + 2uuᵀ + 2w[u]×`** (`NUMERICS.md` §1). Measured by making
   `to_matrix` the unit form `(1 - 2‖u‖²) I + 2uuᵀ + 2w[u]×` and rerunning `just conformance` at
   both precisions; only the `Ad` ids move:

   | id | strata better | worse | ratio, new / old |
   |---|---|---|---|
   | `sen3_ad_n1`, `_n2`, `_n3` | 22 each | 20 each | 0.28–2.12 |
   | `sim3_ad` `f64` | 36 | 51 | 0.28–2.02 |
   | `sim3_ad` `f32` | 43 | 32 | 0.45–1.94 |

   The unit form is better at `θ ≤ 1` and worse at `θ ≥ 3`, per record and pooled per id
   (`cargo xtask compare`, `p` ≤ 1.1e-6), though not on every stratum's maximum. The cause is not
   the rounding of `w²`: a corpus quaternion is a rounded unit one, `‖q‖² = 1 + η` with
   `|η| ≈ u`, and every reference evaluates `R(q/‖q‖)`. The homogeneous form errs by `ηR̂`,
   norm `√3|η|` at every angle; the unit form by `η(R̂ - I)`, zero at the identity and `2√2|η|`
   at `π`. Each form wins one end of `θ`, and an oracle that normalizes wins both, for a `sqrt` and
   a division per element that `0058` declined. A per-entry choice between the forms would halve the
   diagonal's term and leave the off-diagonal `ηR̂ᵢⱼ` of both; it is not worth a new formula. The
   homogeneous form also keeps §12's `to_matrix` contract, `(1 + η)R`.
2. **`so3_act/q:nonunit/f64` and `sen3_ad_n1/q:nonunit/f64` are excepted**
   ([`0046`](./0046-explained-by-record-needs-a-record-to-point-at.md)), citing this record: the
   reference is the rotation of `q/‖q‖`, Sophus's answer is that rotation, and `helicoid`'s is
   `0058`'s stated first-order term, `Rv + η(Rv - v)` for `act` and `(1 + η)R` in `Ad`'s blocks,
   within `0058`'s measured bounds. `sen3_log_n1/q:nonunit` (1.717 against 1.646) is **not**:
   `Log` is scale-invariant (§12), so drift does not explain it.

## Consequences

- `just envelope` reads 707 domination failures, the two rows above excepted; no output moves.
- `Quat::to_matrix`'s rustdoc states decision 1, so the unit form is not proposed again on a
  replay of one stratum.

## Implementation plan

1. One PR (#145) on top of #144: this record, the two rows of `conformance/baseline/exceptions.toml`,
   the rustdoc. Verified by `just lint` (the rows' record is `ready`) and `just envelope`.

## Open questions

None.

## Further work

- ~~`cargo xtask compare` (`0038` decision 4) would have answered decision 1 as a direction and a
  magnitude per stratum, with a twin of the unit form, instead of this record's hand tally.~~
  **Done** ([`0078`](./0078-a-comparison-drops-ties-and-prints-both-maxima.md)):
  `just compare --a helicoid --b helicoid:unit-r` reproduces decision 1's table exactly, every
  strata count and the ratio column to two decimals, at both precisions. Per record it confirms the
  direction at each end of `θ` (`p` from 1.1e-6 to 8.0e-15), but not on every stratum: by the
  maximum, `sim3_ad` at `f64` is worse on 20 of its 81 strata at `θ ≤ 1`. `0078` has the rows.
