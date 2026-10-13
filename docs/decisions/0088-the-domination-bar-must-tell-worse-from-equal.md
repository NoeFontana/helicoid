# 0088: The domination bar must tell worse from equal

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** proposes a `PHASE1.md` §8 and `0006` decision 4 change, a deep corpus tier, and
a fix queue; no code change.

## Context

`PHASE1.md` §8's domination bar fails a `(fn, stratum, precision)` whenever `helicoid`'s maximum
exceeds the best oracle's, by any amount. On `main` at `cb27cee`, with all six oracles rerun on the
current corpus, `just envelope` reads **2656 strata paired with an oracle and 794 domination
failures**, plus the 3 that `conformance/baseline/exceptions.toml` excepts. `--bless` writes nothing
while one stands, so **no baseline has ever been blessed**, and the no-regress bar judges none of
the 4417 scored rows. Every numerics change since Phase 3 has therefore been argued by a study of
its own, outside the envelope, and some of those studies were wrong (`0086`, a draft in #162,
printed a table the rows did not reproduce).

[`0046`](./0046-explained-by-record-needs-a-record-to-point-at.md) gave the bar a hatch: a row in
the exception table that cites a `ready` record. Its item 6 refuses to use the hatch for ordinary
losses. [`0038`](./0038-a-program-comparison-is-not-a-bar.md) decision 1 kept the bar on the
maximum. `0037` (draft) measured that a 6-record stratum cannot resolve the median gap, and
proposed more records. This record measures which of the 794 failures are real and why more
records alone cannot clear the others.

## Measurement

**Protocol.** On branch `docs/envelope-statistics` at `cb27cee`:
- `just oracles` reran all six oracles on the current corpus. The in-process candidate and its twin
  came from `just envelope`.
- A scratch hook in `Aggregate::rows` wrote every subject's per-record errors. It is not committed.
  `max_u` is `NUMERICS.md` §11's metric, unchanged.
- Each stratum was paired with its best oracle: the one with the smallest maximum, as the bar
  chooses it. The pairs `(c_i, o_i)` run over the records both scored.

**The statistic.** `T = max c − max o`, which is the bar's own quantity. Under the null hypothesis
that the two programs are equally good on this stratum, each pair is exchangeable, so swapping
`c_i` and `o_i` changes nothing. The p-value is `P(T' ≥ T)` over the swaps. It is enumerated
exactly up to `n = 16` and drawn from 20 000 seeded swaps above that. The test is paired, so the
correlation between two programs on one input is kept, and nothing is averaged.

### 1. The bar fails a program that is exactly as good as its oracle

Summed over the 2656 paired strata, `P(T' > 0)` under that null is **1098**. That is the number of
failures to expect from a candidate whose errors are exchangeable with its best oracle's. The
library has 797. The number does not depend on `n`: for exchangeable pairs, `P(max c > max o) =
½` at every sample size. More records sharpen each maximum, but they never decide the order of two
maxima that converge on one supremum. Selecting the *best* of up to six oracles makes it worse.

The control is the oracles themselves, judged by the bar against one another on the strata they
share:

| subject | against | shared strata | bit-identical | the bar fails | resolved at `p ≤ 0.001` |
|---|---|---|---|---|---|
| Sophus | GTSAM | 1251 | 224 | **587** | 38 |
| GTSAM | manif | 1466 | 303 | 413 | 25 |
| Sophus | sophus-rs | 484 | 36 | 152 | 1 |
| Sophus | manif | 786 | 254 | 260 | 9 |
| `helicoid` | `helicoid:host-std` | 3053 | 2884 | 80 | **0** |

Sophus fails the bar against GTSAM on 47% of their shared strata. The last row is `helicoid`
against itself on the host's transcendentals: the twin of `0037` (draft). It fails 80 strata and
resolves none.

### 2. Three in four paired strata cannot resolve anything

A swap test on `n` records cannot reach a p-value below `2⁻ⁿ`. At `n = 4`, every Sim(3) cell, that
floor is 1/16. At `n = 6`, every SE_N(3) and SE(3)-chart cell, it is 1/64.

- **2083 of the 2656 paired strata (78%) have `n < 10`**, so they cannot reach `p ≤ 0.001` whatever
  the errors are.
- **609 of the 797 failures sit there.**

That is `0037`'s finding generalised, and it is the bar's resolution, not a property of any kernel.

### 3. The triage

Every failure, by `0037`'s attribution (draft), the envelope's own words, and by whether the swap
test resolves it at `p ≤ 0.001`:

| attribution | resolved | unresolved, `n ≥ 10` | unresolvable, `n < 10` | total |
|---|---|---|---|---|
| program: a `libm`-crate oracle also wins | 14 | 14 | 28 | 56 |
| D16: the host-`std` twin closes it | 0 | 6 | 11 | 17 |
| neither: the twin does not close it | 16 | 114 | 495 | 625 |
| host-`std` oracle, no twin row (all `@f32`) | 1 | 21 | 74 | 96 |
| excepted (`0046`, `0076`) | 1 | 1 | 1 | 3 |
| **total** | **32** | **156** | **609** | **797** |

By family, the resolved failures are SO(3) (26 at `f64`, 1 at `f32`) and SE_N(3) (5, all
`sen3_log` at 12 records). Sim(3), the SE(3) charts, `chol` and `quat_renormalize` resolve none.
The resolved failures have a median ratio of maxima of **1.73**; the unresolved ones 1.26, p90
2.01.

Not one D16 failure resolves. `0032` (draft) proposed a third verdict for exactly that class, but
the test already leaves it unresolved.

### 4. What `n` buys, measured

The 122 failing strata with `n ≥ 64` were subsampled to `m` records, 100 draws each, with the test
rerun on each subsample. "Real" means the 32 strata resolved at the full `n`.

| `m` | real resolved at `0.001` | real resolved at `0.01` | other resolved at `0.001` |
|---|---|---|---|
| 6 | 0.0% | 0.2% | 0.0% |
| 12 | 5.9% | 21.8% | 0.0% |
| 24 | 31.7% | 57.8% | 0.8% |
| 32 | 49.7% | 69.1% | 1.2% |
| 48 | 64.4% | 84.2% | 0.8% |
| 64 | **76.7%** | **94.9%** | 0.1% |

### 5. The budget

The corpus is 62.3 MiB against [`0071`](./0071-the-corpus-budget-is-64-mib.md)'s 64.
- Writing an exactly-zero entry as `0` (`0037` item 3, draft) recovers **9.4 MiB**: 309 523 entries
  of 33 characters each.
- Raising every stratum below 64 records to 64, with the zeros already written as `0`, comes to
  **309 MiB**: Sim(3) 141 MiB, SE_N(3) 115 MiB, the SE(3) charts 12.5, SE(2) 8.6.
- With xz the current corpus is 9.1 MB, so the size is digits, not information.

64 records will not fit in the committed corpus at any plausible budget.

## Decision

Proposed, not settled.

1. **Domination fails only when the loss is resolved.** A paired stratum fails when two things hold:
   - `helicoid`'s maximum exceeds the best oracle's, as today;
   - the paired swap test on `T = max c − max o` rejects exchangeability at **`α = 0.001`**.

   **What the bar does not do:**
   - It averages nothing; the statistic is a difference of two maxima.
   - It does not pass a candidate for being typically fine. A stratum lost on its worst records
     still fails, and the 32 resolved failures in finding 3 are exactly those.

   **Making it deterministic:**
   - The test enumerates exactly to `n = 16`. Above that it uses 2¹⁶ swaps drawn from a fixed
     `splitmix64` seed, the generator's own.
   - The verdict is therefore a function of the rows alone, like everything else under D16. A false
     failure is a stable one, which a reviewer reads once.

   **What the envelope prints:** it still lists every unresolved excess, with its `p`, as a report
   line rather than a failure, so nothing disappears from view.
2. **A stratum that cannot resolve is not a domination stratum.** Where `2⁻ⁿ > α` (`n < 10`), the
   envelope prints "unresolvable at `n`" and reads the stratum from the deep tier (decision 3). It
   is neither passed nor failed.
3. **A deep tier, for the bar that is already local.** `conformance/corpus-deep/` is generated by
   `just corpus-deep` and git-ignored.
   - **Contents:** each id whose strata hold fewer than 64 records, raised to 64.
   - **Nested:** the committed records are a byte-identical prefix, as `0037` found the generator's
     draw already is.
   - **Checked:** a committed `MANIFEST-deep.json` holds each file's SHA-256, and `just
     corpus-deep-check` compares against it.
   - **Who reads it:** the oracles and `envelope` read the deep tier where a file exists. No-regress
     stays on the committed corpus, where it is exact and runs in CI.
   - **No CI change:** the oracles and the envelope run only locally by the owner's rule, so the
     deep tier adds nothing to CI and nothing to the clone.
4. **No-regress is unchanged.** It stays exact against the committed baseline. D16 makes that
   meaningful, and it is the bar this record exists to switch on.

## Rationale

**Larger strata alone** fail by finding 1. At any `n`, an equally good candidate fails about half
its paired strata, and more than half against the best of several oracles. Finding 2 makes it
necessary to raise `n`, but raising `n` is not enough on its own.

**Keeping the bar and excepting the rest** (`0046`'s hatch) would take 763 rows today, each
needing a `ready` record. That list would have to be regenerated every time an oracle or a kernel
moves. `0046` item 6 refuses this for good reason: it is a waiver list.

**A tolerance in `u`** fails both ways:
- It passes a resolved 1.1× defect on a large stratum.
- It still fails a coin flip at 1.3× on four records. The unresolved failures reach 2.01× at p90,
  because the maximum of four draws is that noisy.

**A false-discovery rate across the run** (Benjamini–Hochberg) makes one stratum's verdict depend on
how many others failed. A bar has to be local.

**A per-record majority or unanimity rule** (`0037`, draft; `0038`) answers a different question,
"is the candidate usually worse?". `0038` rejected it for a gate, rightly.

The swap test is the comparison `0006` asks for, the worst case against the worst case, made to
distinguish a loss from a coin flip. `0033` already gave latency the same shape: a paired statistic
and a floor measured beside the comparison.

**`α = 0.001` over `0.01`:**
- About 2000 paired strata are in no position to be worse, so `0.01` would plant around 20 stable
  false failures, each an investigation. `0.001` plants about 2.
- At 64 records, `0.001` still finds 77% of the real losses, against 95% for `0.01`.
- One kernel defect shows on several strata. `so3_from_matrix`'s resolves on 11 of 31, so a missed
  stratum rarely hides a defect.

## Consequences

- **Under this bar today, domination fails 31 strata, not 794.** The 32nd resolved one,
  `so3_act/q:nonunit`, is already excepted by `0076`. Nothing in this record moves a maximum.
- **The fix queue** is the resolved failures not already excepted, worst ratio first within each id:

  | id | strata | worst stratum: `helicoid` against the best oracle | attribution |
  |---|---|---|---|
  | `so3_from_matrix` `f64` | 11 | `theta:1e-3`: 1.408 against 0.494 `u` (`tf_tree_math`), 2.85× | program, all 11 |
  | `so3_jr`, `so3_jl` `f64` | 4 + 4 | `theta:pi-1e-8`: 4.097 against 2.363 (sophus-rs), 1.73× | neither |
  | `sen3_log_n2`, `sen3_log_n1` `f64` | 3 + 2 | `rho:1e-6/theta=1e-1`: 1.995 against 0.888 (manif), 2.25×, all 12 records worse | neither / program |
  | `so3_gamma2` `f64` | 2 | `theta:pi-1e-3`: 1.459 against 1.091 (GTSAM) | neither |
  | `so3_log` `f64` | 2 | `theta:1e-3`: 2.376 against 1.631 (sophus-rs) | `0032` (draft) measured `atan2` |
  | `so3_jr_inv`, `so3_jl_inv` `f64` | 1 + 1 | `theta:pi-1e-6`: 1.283 against 0.892 (Sophus) | neither |
  | `so3_from_matrix` `f32` | 1 | `theta:dense@f32`: 1.370 against 0.825 (Sophus) | no twin row |

  `so3_from_matrix` is 12 of the 31. `tf_tree_math`, on the same `libm`, reads 0.49 `u` where
  `helicoid` reads 0.9–1.4 on six decade strata, so this is the program and not D16. That is
  evidence for `0086` (draft, #162).
- **D16 costs no resolved stratum.** `so3_log`'s two are the only candidates, and `0032` (draft)
  already attributes them to `atan2`. `0032` option 1 becomes unnecessary for the gate.
- **Under the deep tier,** Sim(3), SE_N(3), the SE(3) charts and SE(2) are judged for the first
  time. Today the bar there is a coin flip with no power. Expect real losses to appear: at 64
  records, finding 4 says the test finds three in four of them.
- `0006` decision 4 and `PHASE1.md` §8's first bullet change wording, and `0038` decision 1 is
  reopened on evidence that record did not have (findings 1 and 2). That is the owner's call
  (open question 1).
- `ENVELOPE.md` gains a section of unresolved excesses with their `p`, beside the exceptions
  section. That keeps every number `0046` item 5 wanted visible.

## Implementation plan

1. **Per-record errors as a harness output:** `cargo xtask conformance --records`, one line per
   record and subject, beside the result CSV. They are read by the envelope only and never by
   no-regress. Verified by a test that their per-stratum maximum is the row's `max_u`, bit for bit.
2. **The swap test in `envelope/bars.rs`,** with its seeded draw. Verified by:
   - a unit test that two exchangeable subjects never fail;
   - a planted defect, `seeded:correct` with one stratum's errors doubled, that fails at
     `n = 64` and is "unresolvable" at `n = 4`;
   - reproducing this record's 32 / 156 / 609 split on today's rows.
3. **The deep tier:** `just corpus-deep`, `corpus-deep-check`, the manifest, and the oracles and
   envelope reading it. Verified by the nested prefix being byte-identical to the committed corpus,
   and by `corpus-deep-check` on a second run.
4. **The fix queue,** one record per id where the fix is a formula, starting with `so3_from_matrix`
   (`0086`, draft), then `so3_jr`/`so3_jl` near π.
5. **Bless,** owner only: the baseline and `ENVELOPE.md`, when the bar and coverage pass.

## Open questions

1. **Does the owner accept a resolved-loss domination bar** (decision 1)? It changes `0006`
   decision 4 and `PHASE1.md` §8. The hard rule's words, "on the max — never a mean", hold, and the
   bar no longer fails a coin flip. Rejecting it leaves the hatch, 763 exception rows, as the
   only route to a blessed baseline.
2. **`α = 0.001` or `0.01`?** That is about 2 or about 20 stable false failures, against 77% or 95%
   power at 64 records.
3. **The deep tier at 64 records, git-ignored with committed hashes,** or a larger committed budget
   than 64 MiB? The latter costs 245 MiB more of clone and `corpus-check` time for a bar that runs
   only locally.

## Further work

- `helicoid:host-std` has no `f32` kernel, so the 96 `@f32` failures carry no D16 attribution. A
  binary32 twin would give them one.
- Writing zeros as `0` in the committed corpus (`0037` item 3, draft) frees 9.4 MiB. That would
  lift SE_N(3) to 12 records, which is resolvable at `0.001` (`2⁻¹² = 2.4 × 10⁻⁴`), with no deep
  tier needed for it.
- The swap test is also the per-stratum statistic `0038` decision 2's `cargo xtask compare` lacks:
  a paired verdict on the maxima beside its sign test on the records.
