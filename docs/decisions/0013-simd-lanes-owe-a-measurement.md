# 0013: SIMD lanes owe a measurement

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** none; a draft authorises nothing.

## Context

0003 makes a lane `Real` possible without touching any algorithm: `Mask::decide` blends instead of
branching. Batched consumers exist in principle (`tf_tree`'s `at_many`, batched residuals in BA and
calibration). `tf_tree` 0016 is the evidence against doing it speculatively: `pulp` for ~12% on one
batch path cost nine dependencies and was rejected; `-C target-cpu=x86-64-v3` made the path slower;
autovectorization over `[f64; 4]` was a pessimization.

## Decision

None yet. The record exists so the request has a place to be argued.

## Rationale

—

## Consequences

Until `ready`: no lane type, no SIMD dependency, no `target-cpu`, and no algorithm written "for
lanes" beyond what 0003 already requires.

## Implementation plan

None until the open questions are resolved.

## Open questions

1. **Which consumer, which call, how large a batch?** A named gate in that consumer, not a
   microbenchmark here.
2. **What bar?** A win on that consumer's gate that pays for the dependency cost, at least as
   demanding as 0016's (a ~12% win did not pay for nine crates).
3. **Which mechanism?** `core::simd` is nightly (D17 excludes it); `wide` and `pulp` cost
   dependencies (0007); a hand-written lane type on stable must beat 0016's autovectorization
   result.
4. **Determinism:** lane results must be bit-identical to scalar ones (D16) — do the candidate
   mechanisms guarantee lane-wise `libm` semantics?
