# 0035: Validation has three layers, and one of them needs a quiet machine

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** proposes `bench-gate --record`/`--replay`, a committed replay fixture, and the
two controls as ordinary tests. No change to the gate's statistic or its verdict.

## Context

`PHASE1.md` §9's gate works and is validated
([`0033`](./0033-a-latency-floor-is-measured-beside-the-comparison.md), draft, both controls). It is
also unusable as a development loop. Measured on this host:

| unit | windows | cost |
|---|---|---|
| one window (0.5 s warmup, 1.5 s collect, process start, criterion's analysis) | 1 | ~2.2 s |
| one `b₁ c b₂` triplet | 3 | ~6.6 s |
| one benchmark, `REPLICATES = 3` | 9 | ~20 s |
| the 60-benchmark set | 540 | ~20 min |
| one validation cycle: positive **and** negative control | 1080 | **~40 min** |

Six cycles were spent in one day reaching `0033`'s conclusions — roughly four hours — and every one
of them was answering a question about a **pure function**. The gate's decision rule is
`median` → `ratio_ci` → `delta` → "every pair above `1 + floor`": a deterministic, seeded function
from sample vectors to a verdict. None of it needs a CPU to be quiet. What needs a quiet CPU is
*obtaining* sample vectors.

The two are currently welded together, and `bench::clear_slots` is where the damage is done: it
deletes every window's samples immediately after a single read. A 20-minute run produces 540 vectors
of 100 per-iteration times — 54 000 measurements — and discards all of them, so the next iteration
re-measures from scratch. Each of `0033`'s findings (the mean's 0.7172 floor, pooling's 0.5884, the
asymmetric `δ`, the refuted position bias, the 30-second independence failure) was paid for at full
price, and none of them left a regression test behind.

## Decision

Proposed, not settled. **Three layers, distinguished by what they can possibly be wrong about.**

**L1 — the rule, in milliseconds, in `just test`.** The decision rule is validated by replaying
*recorded* sample vectors, not by running benchmarks.
- `bench-gate --record <dir>` persists every window `--against` measures, as
  `<group>/<bench>/rep<N>/<slot>/sample.json`. The slow run stops throwing its evidence away.
- **The stored format is criterion's own `sample.json`, copied byte for byte**, and the reason is
  fidelity rather than size. A copy is lossless by construction; it needs no second parser; a
  recording can be diffed against a live `target/criterion` tree; and a replay goes through
  `samples::read_file`, so the `times`/`iters` division and its validation are *replayed* instead of
  assumed. Storing the derived per-iteration quotients — the first design — would bake today's
  derivation into the recording and make precisely that step unreplayable. `iters` is redundant under
  criterion's `Linear` mode (an arithmetic sequence expressible in three numbers) and is kept anyway,
  because dropping it would end the byte-copy property to save a few kilobytes.
- `bench-gate --replay <dir>` runs the whole rule over a recording and prints the same report with
  the same exit status, executing no benchmark.
- A **committed replay fixture** — a stratified subset, both precisions × three θ strata — in two
  recordings: one A/A (the same binary both sides) and one with a known effect (the `exact_a`
  slowdown). The two controls then become ordinary `#[test]`s, deterministic because `SEED` is
  fixed, and `just test` runs them in milliseconds.
- **Synthetic contamination**, one property per test: a single disturbed window (no verdict); a
  disturbance covering one whole triplet (no verdict — the shape that defeated the rule when a
  benchmark's triplets were consecutive); every window of every triplet slow together (drift, no
  verdict — what bracketing buys); and a candidate slower in every triplet (a verdict). The last is
  also the rule's **limit**, because a disturbance confined to the candidate's windows across the
  whole run has that same shape and nothing in the data separates them; the floor's own size is the
  only report of it, which is why `--against` prints it every run.

**L2 — the host, in a minute, before trusting any run.** A stratified `--aa` over a handful of
benchmarks answers "what can this session resolve at all?", which is a property of the machine and
not of the code. ~8 benchmarks is ~2 min. It reports; it never gates.

**L3 — certification, rarely.** The full 60-benchmark `--against` for a real candidate, or to
re-characterise the host after a toolchain, `criterion` or CPU change. This is the only layer that
costs 20 minutes, and it is the only one that has to.

The fixture is **test input, not a baseline.** `baseline/bench/` stays `.gitignore`d because a
timing is one core under one load and does not travel (D16's reasoning applies to accuracy, and its
converse applies here). A replay fixture makes no claim about any host: it is a committed vector of
numbers fed to a pure function, exactly as `conformance/`'s corpus is, and it is reviewed for shape
rather than for value.

## Rationale

The layering follows from what each thing can be wrong about. A statistic can be wrong about
arithmetic — that is checkable without a machine, and so it must be checked without one, because a
check that costs 40 minutes is a check that gets skipped. A protocol can be wrong about *what data
it collects* — adjacency, bracket order, loop order — and that can only be caught against a real
machine, but it changes rarely. A host can be wrong about being quiet, and that is true or false per
session, so it is a precondition and not a test.

The alternative is to keep the controls as end-to-end runs only. It is simpler and it is what exists;
its cost is measured above, and its real cost is worse than the wall clock: the 30-second
independence failure was found by *accident*, on a session that happened to be noisy enough to
expose it, three estimator revisions after the mechanism had been introduced. Under L1 it is a
four-line test that fails deterministically.

A second alternative — shorten the windows — is refused. Window length and sample count determine
the noise the rule must survive, so a cheaper window measures a different problem, and L3's job is
to measure *this* one.

## Consequences

- Estimator iteration goes from ~40 min to under a second, and every finding in `0033` becomes a
  permanent test rather than institutional memory.
- `just test` gains the two controls, so a change to the rule that breaks either one fails in CI
  (when CI returns) rather than in a human's afternoon.
- A recording is an artifact with provenance: host, `criterion` version, the two binaries' hashes.
  Without that it is a pile of numbers whose meaning decays.
- L1 cannot catch a protocol error, by construction. Anything about *ordering in time* — adjacency,
  bracketing, loop order — still needs L3, and the record should say so rather than let a green
  `just test` imply the gate is sound.
- The fixture must be small enough to review and large enough to be realistic. A stratified subset
  (6 benchmarks × 9 windows × 100 samples) is ~5 400 numbers.

## Implementation plan

1. `--record <dir>`, copying every window; `--replay <dir>`, running the rule over a recording —
   verified by a round-trip test over a synthetic recording in criterion's own shape. `--aa`
   recording is owed; only `--against` records today.
2. The committed replay fixture, from one L3 run of each kind, with its provenance header.
3. The two controls as `#[test]`s over the fixture: A/A must produce no failure; the `exact_a`
   recording must produce exactly its five rows.
4. The synthetic contamination tests, one per property, including the two that `0033` paid 40
   minutes each to discover.
5. `--stratum <n>` or an equivalent subset selector for L2, and a `just bench-precheck` recipe.
6. Only then, the owed `--against <ref>`.

## Open questions

1. **Does a cheaper L3 exist at all?** Two separate processes cannot interleave finer than a window,
   ~2 s. One binary that contains *both* versions could interleave at the sample level, ~15 ms —
   Cargo permits depending on the same crate twice under different keys, one by `path` and one by
   `git`/`rev`. That pairs 100× finer, likely removes the need for the bracket *and* for replicates,
   and would cut the full set from ~20 min to ~2–4 min while subsuming `--against <ref>`. It also
   puts both versions in one layout, so the relative-layout confound `0033` measured (a 64-byte
   `.text` shift moving four untouched benchmarks by 1.6–4.1%) changes character rather than
   disappearing. Worth measuring before L3's shape is fixed; a Cargo-surface decision, so not this
   record's to take alone.
2. Should the replay fixture live under `xtask/src/bench/fixtures/` and be `include_str!`'d, making
   the tests hermetic, or under `conformance/` beside the other committed corpus? The first is
   simpler; the second is where a reader would look.
3. Does L2's subset need to be *measured* as representative, or is "both precisions × the fastest,
   a middling and the slowest stratum" enough? The honest answer is that a subset's floor is not the
   set's floor, and L2 only claims an order of magnitude.
