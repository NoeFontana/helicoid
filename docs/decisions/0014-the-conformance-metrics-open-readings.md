# 0014: The conformance metric's open readings

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** none; a draft authorises nothing.

## Context

`NUMERICS.md` §11 defines the forward error and is silent on several things a harness must decide
before it can score a record. `xtask/src/conformance/metric.rs` (`TABLE` and its module docs) took
the smallest reading of each so that the harness core runs (`PHASE1.md` §5). The readings fix the
bar every envelope baseline is blessed against (`0006`); none is blessed yet (`PHASE1.md` §8 is not
started), so each is still free to change.

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
