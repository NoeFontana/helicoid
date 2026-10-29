# 0109: Determinism is required before Phase 6

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #210

## Context

D16 promises bit-identical outputs on x86_64, aarch64 and wasm32 for a given `helicoid` and a given
resolved `libm`. `PHASE6.md` §3 makes `just determinism` its check. CI already runs it:
- `determinism`, a matrix over x86_64, x86_64 at `ndebug`, aarch64 (`ubuntu-24.04-arm`) and
  `wasm32-wasip1`, each writing one digest table;
- `determinism (compare)`, which compares the four tables.

`PHASE1.md` §3 makes both required "from Phase 6", and neither is a required check today. A change
that moves one target's bits can therefore merge on a red `determinism (compare)`.

locus-calib:0007 requires `result.json` to be bit-identical on x86_64 and aarch64. `helicoid`'s share
of that promise is D16, and this job is the only evidence for it. The first consumer that depends on
D16 for its own contract has arrived before Phase 6.

**Measured** (`gh run list --branch main --workflow ci.yml`, the 30 most recent runs):
- `determinism (compare)` succeeded on all 30.
- On the latest run, each digest leg took 36–58 s and the comparison finished 1 min 40 s after the
  start, inside the 3 min 49 s – 4 min 47 s that the two required `test` legs took. Requiring it
  adds no wall-clock time to a pull request.

**A skipped job reports success.** GitHub reports a job skipped because one of its `needs` failed
as successful, and a required check that is skipped does not block a merge. As written,
`determinism (compare)` has `needs: determinism`. If the aarch64 leg fails, the comparison is
skipped, and requiring the comparison alone would pass.

## Decision

1. **`determinism (compare)` is a required check now**, beside the nine of the ruleset. The digest
   legs are not named in the ruleset: the comparison stands for them (decision 2).
2. **The comparison cannot be skipped into a pass.**
   - The job runs `if: ${{ !cancelled() }}`.
   - Its first step fails unless `needs.determinism.result` is `success`, so a failed or skipped
     leg makes it red.
   - A cancelled run still reports no verdict, as every other job does.
3. **Spec text.**
   - `PHASE1.md` §3's "required from Phase 6" becomes required, citing this record.
   - `PHASE6.md` §0.0's row says the job is required.
   - `PHASE6.md` §7's checkbox "Determinism required and green" and its appendix item 3 stay
     Phase 6's: green on every target is still Phase 6's criterion, and the requirement is met early.
   - `ci.yml`'s comment above the job says the same.

## Rationale

- **Requiring the legs too.** Five required names for one property. A renamed matrix entry
  silently stops being required. The comparison, made unskippable, already fails on any leg.
- **Waiting for Phase 6.** Phase 6 is gated on work locus-calib does not need (`PHASE6.md` §1, §2).
  Meanwhile every merge to `main` could break the one property that consumer's contract names.
- **Amending `0007`.** `0007` is the dependency budget. D16 is the decision this check enforces,
  and `PHASE1.md` §3 is the text that schedules it. A short record keeps the reason findable.

## Consequences

- A pull request that moves any target's bits is blocked until the digests agree again.
- The ruleset change is the owner's: CI code cannot make a check required.
- Every new corpus id enters the digests as soon as the shipped subject supports it
  (`determinism::digests` walks the manifest), so later ids are covered without a CI change.

## Implementation plan

1. The `ci.yml` change and the spec text (decisions 2 and 3). Verified by the workflow passing on
   the pull request.
2. The ruleset (decision 1), applied by the owner. Verified by a throwaway draft pull request whose
   aarch64 output differs by one ulp: `determinism (compare)` is red and the merge is blocked.

## Open questions

None.
