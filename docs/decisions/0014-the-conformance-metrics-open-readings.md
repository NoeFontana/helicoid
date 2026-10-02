# 0014: The conformance metric's open readings

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** none; a draft authorises nothing.

## Context

`NUMERICS.md` §11 defines the forward error and is silent on several things a harness must decide
before it can score a record. `xtask/src/conformance/metric.rs` (`TABLE` and its module docs) took
the smallest reading of each so that the harness core runs (`PHASE1.md` §5). The readings fix the
bar every envelope baseline is blessed against (`0006`); none is blessed yet (`PHASE1.md` §8 is not
started), so each is still free to change. `xtask/src/thresholds` (`PHASE1.md` §6) did the same for
the threshold sweep (questions 8 to 11); no `generated.rs` is written from it yet.

## Decision

None yet. The record lists the readings so they are argued before the first baseline; each one
that survives is a `NUMERICS.md` §11 edit.

## Rationale

—

## Consequences

Until `ready`: the readings are provisional, `just envelope --bless` does not exist, and no baseline
is committed. A changed reading changes every `max_u` under it.

## Implementation plan

None until the open questions are resolved.

## Open questions

1. **Floor of a tangent or a coefficient.** The code uses `2^-1022`, so `theta:1e-k` reads a
   relative error (`EA.4(b)`, which argues it for `Log`'s tangent and for subnormals only). Is one
   floor right for the coefficients too?
2. **SE_N(3) tangent at `w = +0`.** The code gives it no sign forgiveness: `−τ` is not the logarithm
   the input's sign selects. `(−φ, J_l(−φ)⁻¹x)` is also a logarithm of the same `X`. Should §11's
   sign invariance cover it?
3. **One score per record.** The code takes the maximum over a record's output fields. Per-field
   rows would say which field scored.
4. **Whole-tangent norm.** §11 says vectors take the 2-norm, so a tangent's rotation and translation
   blocks share one norm: at `rho:1e4` a rotation error many `u` of `φ` reads under 1. Score the
   blocks apart?
5. **The ρ-scale.** The code takes the 2-norm of the input tangent's translation entries (all `N`
   blocks together) for `sen3_exp.x` and `se2_exp.t`. §11 says "‖ρ‖-scale".
6. **Overflow.** An `F` past binary64 (a gross error against a zero reference) counts as non-finite
   and fails the run. Is that the intended gate, or a failing score of its own?
7. **The field of the `b` defect's curve** (`PHASE1.md` §10, not §11). The self-test fits the `value`
   field's per-stratum `max_u` over `theta:1e-8`…`theta:1e-2` (`p` = 1.937). Through `Dual` the
   exact arm's derivative errs by θ⁻⁴, so a record's maximum over `value` and `d_branch`, which is
   what a results row holds, fits `p` = 3.98 and misses the window. Score the field, or widen the
   window and score the record?
8. **Per decade of what** (`PHASE1.md` §6). "A log grid of 64 points per decade" of the branch
   variable, "spanning θ ∈ [1e-8, 1]". The sweep takes decades of `z = θ²`: 1025 points,
   `10^((i−1024)/64)`, 128 per decade of θ. Decades of θ are 513 points, the even indices; they move
   `a` to θ = 0.835 (37.8 u against 25.8), `b` to the top (64.4 u against 63.5) and `e`'s switch to
   0.965 (the same objective), and leave `k`, `c`, `d`. Which is meant?
9. **The top of the grid.** It ends at `z = 1` and a candidate takes the series arm where
   `z < switch`, so the corpus record at θ = 1 exactly is always on the exact arm and no switch
   above θ = 1 is tried. That record is `e`'s objective (8.69e3 u, 7 terms); a top of `nextUp(1.0)`
   gives `e` 8 terms and 2.60e3 u, and leaves `k`, `a`, `b`, `c`, `d`. Is the range closed at 1, or
   does the grid extend to the first switch that puts that record on the series arm?
10. **The D12 prior for all six.** §6 describes it as one switch for `a`, `b`, `c` and evaluates it
    as a named candidate; the sweep scores it for all six coefficients, as `seeded:correct` runs it.
    Report it for `k`, `d`, `e`?
11. **Ties.** Exact `f64` equality of the objective (`tied` counts them). A tolerance would make
    near-equal candidates tie, and the cheaper one win. Exact, or a tolerance?
