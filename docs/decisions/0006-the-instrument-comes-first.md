# 0006: The instrument comes first

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

"State of the art" is a claim per regime: near identity, near π, large translation, joint
limits. Finite differences cannot see a near-identity Jacobian error smaller than their own step
error; hand-typed high-precision tables (`tf_tree_math`'s sweep) do not scale to Jacobians, SE_N(3)
or Sim(3); and a mean hides the stratum where a library is wrong. The stack's standing rule
applies: instrument before algorithm — a measurement built after the code tends to measure what the
code already does.

## Decision

1. **Corpus from definitions.** A Python + pinned mpmath generator at 120 digits computes every
   reference from its definition (matrix exp/log, defining series, high-precision inverses), never
   from `NUMERICS.md`'s closed forms; inputs are exact binary64; the corpus is committed and
   regenerates byte-identically (`PHASE1.md` §2, §4).
2. **Strata** are named, versioned input families; every bar is per function, per stratum, per
   precision, on the max (p99 reported).
3. **Oracle envelope.** `tf_tree_math`, sophus-rs, Sophus, manif, GTSAM through one file protocol;
   pinned versions; their errors recorded.
4. **Bars.** Domination: `helicoid`'s max ≤ the best oracle's max. No-regress: ≤ the committed
   baseline, compared exactly (D16). Coverage: every public numeric id has a corpus file.
5. **Seeded defects** prove the instrument detects what it claims to (`PHASE1.md` §10).
6. **Evidence** is generated into `docs/evidence/ENVELOPE.md` on every bless.

## Rationale

Harness independence is what makes the corpus a specification rather than a second
implementation. Exact comparison against the baseline is only possible because outputs are
deterministic; that is why D16 and this record are one design. Alternatives lost: tolerance-based
unit tests (assert constants, hide regressions below the tolerance); FD Jacobian checks (blind
near identity); a single oracle (its bugs become ours).

## Consequences

- Phase 1 ships no group code; every later phase starts by adding corpus ids.
- Oracles run nightly in containers; the in-process harness runs on every PR.
- A stratum is never deleted or narrowed because it fails.

## Implementation plan

`PHASE1.md` appendix, steps 2–7 — verified by `PHASE1.md` §11.

## Open questions

None.
