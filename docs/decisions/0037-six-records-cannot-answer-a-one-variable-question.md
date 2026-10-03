# 0037: Six records cannot answer a one-variable question

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** the twin and its two controls; the rest proposes `PHASE1.md` §4.2, §4.4 and §5 edits.

## Context

[`0036`](./0036-the-silence-in-5-1-costs-twenty-domination-failures.md) (draft) plan step 2 owes one
measurement: the one-variable swap for `sen3_jl_n1` and `sen3_jr_n1`, 30 of the 45 domination
failures where only a glibc-backed oracle beats the candidate and "libm-bound" is a hypothesis.
[`0032`](./0032-domination-charges-helicoid-for-d16.md) (draft) ran that swap for `so3_log` in a
script that was not kept, and its open question 1 asks whether D16's price is a per-function
constant from `error-analysis.md` EA.13(d) or a per-stratum measurement.

One subject answers all three, so it was built as a subject rather than a script.

## Measurement

### The instrument, and its positive control

`seeded:host-std` is `seeded:correct`'s program at a scalar whose six transcendentals come from
Rust `std` — the host's libm, glibc 2.39 on x86_64 — where every other subject here routes them
through the `libm` crate (D16). The switches, the series, the associations, the application form and
the generic code are shared by construction; six method bodies differ.

Three of the six are not a variable. Over 100 000 samples, `sqrt`, `abs` and `copysign` agree with
the `libm` crate bit for bit, as [`0018`](./0018-libm-arch-is-bit-identical-for-exact-operations.md)
says exact operations must. `sin_cos` disagrees on 5 844 and `atan2` on 16 836. `cbrt` disagrees on
**0** over signed arguments spanning $2^{-64}$ to $2^{64}$, so D16 costs
[`0017`](./0017-cbrt-and-mask-valued-roots.md)'s `cbrt` nothing on this host
(`the_scalar_differs_from_the_libm_crate_in_the_transcendentals_only`).

**The twin reproduces `0032`'s host-glibc column on 8 of 8 `so3_log` strata**, to the four decimals
that record prints, including the one it flags as *not* matching sophus-rs (`theta:1e-7`, 1.8492
against sophus-rs's 1.9396). A retired script is now a test
(`the_twin_reproduces_the_glibc_column_of_0032`), and the whole corpus takes 3.6 s.

### What the swap changes

751 rows, binary64, every id the correct kernel runs. The twin differs from `seeded:correct` on
**15**:

| id | strata that move | direction |
|---|---|---|
| `so3_log` | 11 of 30 | 8 better (the 8 failures), `q:nonunit`, `theta:1e-8`, `theta:1e0` worse |
| `coeff_r` | 2 of 29 | both better |
| `coeff_d` | 1 of 28 | **worse**: `theta:1e0` 830.73 u against the host's 907.92 u |
| `so3_exp` | 1 of 28 | the same `max_u`, a different `p99_u` |

**Every `sen3_*` row is bit-identical**, on all nine ids and all 52 strata each. So of the 75
domination failures: `so3_log`'s 8 close under the swap, as `0032` measured; the other 37 in
`0032`'s class do not move at all.

### Why that is not a refutation

The swap had almost no opportunity. `sen3_jl_n1` is 52 strata of **6 records**; in **1 of its 312
records** do the two libraries disagree on `sin` or `cos` at any argument the kernel reaches, and
that record is in `theta:1e0`. Across the 17 strata at $\theta$ near $\pi$ — 102 records, and 12 of
this id's 15 failures — **no record reaches a disagreement at all**
(`the_swap_has_no_power_on_the_strata_near_pi`). `sen3_jr_n1` and `sen3_exp_n1` read the same.

The rate is not the reason; the sample size is. At $\theta \approx 1$ the two libraries disagree on
**15.8%** of arguments, over $[0.1, 1]$ on 11.1%, and over $[\pi - 0.2, \pi]$ on 2.1% (1 000 000
samples each). Six draws per stratum is 0.95 expected divergent records where the rate is highest.

So `0036` step 2 is not answered by this corpus: for 37 of the 45, "libm-bound" is neither confirmed
nor refuted, and D16 is measured to cost 8 strata, not 45.

### What six records costs the bar

`PHASE1.md` §4.4 states the count and calls it a trade-off: *"Each SE_N(3) stratum is 6 records, not
64: five dense matrices per record at 64 would be 190 MB. 6 is a budget trade-off ... not the largest
count that fits."* It does not say what 6 resolves.

- **66 of the 75 domination failures are on a `sen3_*` id**, so they are a comparison of two maxima
  over 6 samples, under a bar `0006` puts on the max and never a mean.
- `p99_u` equals `max_u` on **741 of 751** rows. That is arithmetic, not chance, and `report.rs`
  already says so — nearest-rank p99 is the max below 100 records — so the schema's second aggregate
  separates from the first on `theta:dense` only.

### Where the records are

The corpus is **32.98 MB** against §4.4's 50 MB cap; the SE_N(3) family is 16.71 MB. Sixty-four
records would be 186 MB, which is §4.4's 190 MB estimate. But most of what is stored is not
information:

- Every record of the twelve `sen3_{jl,jr}{,_inv}_n{1,2,3}` ids and the three `sen3_ad_n*` stores a
  dense matrix whose diagonal 3×3 blocks are **string-identical** — 1 distinct of 4 on 312 of 312
  records at $N = 3$ — and whose entries off the dual-matrix structure are exactly 0: 81 of 144 at
  $N = 3$, 36 of 81 at $N = 2$, 9 of 36 at $N = 1$. Only 36 of 144 entries carry information at
  $N = 3$, a **4.00×** redundancy; 8.56 MB over the twelve Jacobian ids.
- A zero is written `0.00000000000000000000000000000e0`, 33 characters. Across those fifteen ids
  that is **6.87 MB of 15.72 MB — 21% of the whole corpus — spent on the digits of zero.**

The two savings are not equally cheap. **Writing an exactly-zero entry as `0` changes no metric**:
`NUMERICS.md` §11's norms are over the same values, and the subject writes exact zeros in those
positions too (`se3::jacobian` fills `S::zero()`), so each contributes 0 to $\lVert \hat y - y\rVert$
and to $\lVert y \rVert$ either way. It recovers ≈6.2 MB, enough to take every SE_N(3) stratum from
6 records to about 12 and leave the corpus near today's size. **Dropping the repeated diagonal
blocks does change the field the metric norms**, so it is a §11 question and not free.

### What six draws resolve, measured

Open question 2 asked what count a per-stratum maximum needs. It is answerable without changing
anything: regenerate one id at a higher count and compare.

`sen3_jl_n1` was regenerated at **64 records per stratum** with the committed generator, one
constant changed (`gen/strata.py`'s `SEN3_SAMPLES`), 2 min 40 s, 4.95 MB against 0.46 MB. The draw
is **nested**: record 0 of the 64 is byte-identical to the committed record 0, reference included,
so the committed 6 are a prefix of the 64 and the two maxima are comparable rather than two
designs.

**At 64 records the maximum rises on 47 of 52 strata** — median **1.266×**, p90 1.624×, up to
**2.008×** (`theta:1e-1`, 0.4728 u to 0.9491 u). The other 5 tie, which is what a nested draw
means. So a 6-record stratum understates its own maximum by about a quarter, and the understatement
is not uniform.

The variance matters more than the bias, because domination is a comparison. Resampling the 64
records of each stratum, the maximum over `m` draws:

| `m` | median 64-record max / median `m`-draw max | p95/p5 spread of the `m`-draw max | P(two designs differ by > 1.13×) |
|---|---|---|---|
| **6** | 1.230 | **1.59×** | **0.525** |
| 12 | 1.135 | 1.37× | 0.391 |
| 16 | 1.095 | 1.31× | 0.342 |
| 20 | 1.061 | 1.28× | 0.305 |
| 24 | 1.061 | 1.21× | 0.239 |
| 32 | 1.000 | 1.16× | 0.127 |
| 48 | 1.000 | 1.10× | 0.018 |

1.13× is the **median observed domination gap** of the 45 failures in this class (p90 1.746, max
2.361; the 30 algorithmic ones are median 1.306, max 2.304). So **at six records, two designs of
the same stratum, the same subject and the same library disagree by more than the median measured
signal in 52.5% of trials.** Every one of the 75 gaps is under 2.86×, and 2.86× is the p95 ratio
between two 6-draw designs on the 64-record strata the whole corpus offers for this check.

That does not make the gaps false — the comparison is paired, both sides see the same six records,
and a correlated design cancels part of the noise. It makes them **unresolved at this count**, which
is a different and worse thing than wrong: no amount of reasoning about `libm` or about assembly
forms can be checked against a number whose own design variance is larger than the effect.

### What the budget can buy

Exact byte accounting over the fifteen block ids, three schemas, same records:

| reference schema | SE_N(3) family | corpus | records/stratum within 50 MB | with a 15% reserve |
|---|---|---|---|---|
| today | 17.50 MB | 32.98 MB | 11.8 | 10.1 |
| an exactly-zero entry written `0` | 10.83 MB | 26.31 MB | 19.1 | 16.2 |
| a dual matrix written as one | **7.66 MB** | 23.14 MB | **27.1** | 23.0 |

Sixty-four records would be 178 MB, which is §4.4's 190 MB estimate confirmed. **Nothing inside the
50 MB cap reaches the count at which the design noise falls below the signal**: writing a dual
matrix as a dual matrix buys 24 records, where two designs still disagree by more than 1.13× in 24%
of trials. The count that would settle these 37 failures is near 48, and it does not fit.

So the resolution question is not only a byte question. At a fixed budget, the instrument's power is
a property of the **statistic** as much as of the sample size, and a maximum over six records uses
one record and discards five. Under exchangeable paired errors, "the candidate is worse on every one
of the six records" has a false-positive rate of $2^{-6} = 1.6\%$ against the 52.5% above — a 33×
reduction at the same records, and no new bytes. `0006` chose the max over the mean to stop a mean
hiding a bad case, which it does; a paired rule is neither, and `0033` already reached for exactly
this shape (a paired statistic, swap-invariant) for latency while accuracy stayed at a bare max.

## Decision

Proposed, not settled.

1. **A one-variable twin is a subject, not a script.** `PHASE1.md` §5 lists in-process subjects as
   `helicoid` and `seeded:<defect>`. A twin is a third kind: the candidate's own program with
   exactly one thing changed, registered `planted` so no plain run and no bar reads its rows, and
   **validated by reproducing a published measurement** before it is used to make a new one. §5's
   list should say so. A twin is not a candidate; the envelope comparing `helicoid` with a
   differently-rounded copy of itself would measure nothing.
2. **A stratum's record count is part of the instrument, so §4.4 should state what each count
   resolves**, not only what it costs. Six records is six draws of a maximum, and at six draws the
   design noise exceeds the median measured signal more often than not. Where the bar is per
   stratum on the max, the stratum size *is* the bar's resolution.
   **And raising the count alone cannot fix it inside the cap**, so the bar's statistic is in scope
   too: a paired rule over the records of a stratum costs no bytes and cuts the false-positive rate
   by 33× at six records. That is `0006`'s to decide, and this record asks it (open question 6).
3. **Spend the digits of zero on records.** §4.2 should write an exactly-zero reference entry as
   `0`, and §4.4 should spend what that recovers on the SE_N(3) record count.

## Rationale

A script that is not kept cannot be a control. `0032`'s swap was correct and its numbers stood up —
the twin reproduces all eight — but nothing in the tree could re-run it, so `0036` inherited a
hypothesis for 37 strata it had no way to test, and the envelope began annotating them with a claim
about D16 that this record now shows is unsupported for most of them. The cost of making it a
subject was one scalar and one `Subject` impl, because
[`0003`](./0003-the-scalar-that-cannot-say-less-than.md) already made every kernel generic over
`Real` for an unrelated reason.

The positive control is the load-bearing part. A twin that reproduces nothing published is a second
unverified number; one that reproduces eight is an instrument. It also costs almost nothing to keep
honest: the two controls are pinned rather than bounded, so a `libm` or glibc version that moves
these bits fails a test, which is the signal
[`0034`](./0034-a-second-transcendental-backend-owes-a-measurement.md) (draft) asks for.

Reporting the corpus finding here rather than in a record of its own is deliberate: it was not
sought, it is what the measurement found when it came back empty, and the chain from "no change" to
"six draws" to "21% of the bytes are zeros" is one argument.

## Consequences

- **`0036`'s 45 must not be read as 45 libm-bound.** Eight are measured so; 37 are unexplained. The
  envelope's own annotation says *"only a host-`std` oracle beats this: D16 is a candidate, `0032`
  draft"* on all 45, which overstates what is known for 37 of them. That text should say what it
  knows, and the twin's rows are what let it.
- **`0032` open question 1 is answerable per stratum now**, for any id, by this subject rather than
  by EA.13(d)'s constants — and the corpus does not support those constants as a charge: `sin`/`cos`
  move 1 of 751 rows, and that row moves against the host (`coeff_d`, `theta:1e0`, 830.73 u to
  907.92 u). `0032`'s proposed `libm-bound` verdict applies to 8 strata.
- **Raising a stratum's record count re-blesses every baseline for that id**, and the maxima in
  `0032`, `0036` and this record are not comparable across such a change. That is a cost of the
  proposal, not an argument against it: a maximum over 6 draws was never the same quantity as a
  maximum over 12.
- A twin is a diagnostic and stays out of the plain run by `Registered::planted`, so `just
  conformance` and `just envelope` are unchanged until step 2 lands.
- `0017`'s `cbrt` and the exact operations have a measured cost of zero on this host. D16's price,
  as measured, is `atan2` on `Log`-shaped strata and nothing else.

## Implementation plan

1. The twin, the positive control against `0032`'s column, and the power test — verified by
   `the_twin_reproduces_the_glibc_column_of_0032`,
   `the_swap_has_no_power_on_the_strata_near_pi` and
   `the_scalar_differs_from_the_libm_crate_in_the_transcendentals_only`. **Landed.**
2. The envelope reads the twin's rows and each domination failure says which of the three it is:
   beaten by a same-backend oracle, closed by the swap, or unexplained — verified by
   `the_twin_says_whether_the_swap_closes_a_failure_only_glibc_won`, and by the run itself, which
   now reports **30 the program's, 8 D16's, 37 neither**, the partition this record measured.
   `just envelope` reruns `conformance-twin` so the attribution is never read from a stale twin;
   with no twin rows a failure is reported unattributed and no claim about D16 is made.
   **Landed.**
3. §4.2's exactly-zero entry, with `just corpus` regenerating and `just corpus-check` the byte
   comparison. **Owed, blocked on open question 1.**
4. §4.4's SE_N(3) record count, once step 3 has paid for it, with every affected baseline re-blessed
   in the same PR. **Owed, blocked on open questions 1 and 2.**

## Open questions

1. Is `0` still "a decimal string with 30 significant digits" (§4.2), or does the sentence need the
   exception written out? The generator, the parser and the metric all read it as the same value;
   only the prose is at issue.
2. ~~What count does a per-stratum maximum need?~~ **Measured above.** Below about 48 records two
   designs of one stratum disagree by more than the median observed signal too often to read a
   1.1–2.4× gap, and 48 does not fit the cap under any of the three schemas. What is *not* settled
   is the count to adopt given that — 24 (the dual-matrix schema's ceiling with a reserve) buys a
   2.1× reduction in false positives for a full corpus regeneration, and does not make the bar
   sound on its own.
3. **Should a dual matrix's reference be written as a dual matrix?** It is the difference between 19
   and 27 records per stratum, so it is now load-bearing rather than an extra. It needs
   `NUMERICS.md` §11 to say what field the metric norms; `0014` (draft) question 25's sibling, and
   not asked there. Note the zeros alone are metric-neutral and this is not.
4. Should a twin ever be an *oracle* rather than a diagnostic? It is the only "oracle" that is
   `libm`-neutral by construction, which is `0032` open question 3 in a sharper form; it is also the
   candidate's own program, so domination against it is not a comparison with anyone else's
   algorithm.
5. Is there a second twin worth keeping — the three application forms `0036` measured, as subjects
   rather than as a branch that was reverted? Same mechanism, same validation shape.
6. **Does `0006`'s bar need a statistic that uses more than one record of a small stratum?** The
   max is 1 of 6; a paired rule is 6 of 6 and costs nothing. `0006` says "on the max, never a mean",
   and a paired per-record rule is neither — `0033` took that route for latency. This record raises
   it; `0006`'s to answer, and until it does, the 37 failures in this class should be read as
   unresolved and not as an algorithm to chase.
