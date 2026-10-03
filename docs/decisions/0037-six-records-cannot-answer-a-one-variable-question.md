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

## Decision

Proposed, not settled.

1. **A one-variable twin is a subject, not a script.** `PHASE1.md` §5 lists in-process subjects as
   `helicoid` and `seeded:<defect>`. A twin is a third kind: the candidate's own program with
   exactly one thing changed, registered `planted` so no plain run and no bar reads its rows, and
   **validated by reproducing a published measurement** before it is used to make a new one. §5's
   list should say so. A twin is not a candidate; the envelope comparing `helicoid` with a
   differently-rounded copy of itself would measure nothing.
2. **A stratum's record count is part of the instrument, so §4.4 should state what each count
   resolves.** Six records is six draws of a maximum, and §4.4 justifies the number in bytes alone.
   Where the bar is per stratum on the max, the stratum size *is* the bar's resolution.
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
2. **What count does a per-stratum maximum need?** Twelve is what the zeros pay for, not an answer.
   The question is what resolution `0006`'s bar requires of a maximum, and it has never been asked
   of any stratum family — `coeff_*` at 64 and `so3_log` at 128 were not chosen by it either.
3. Should the repeated diagonal blocks go too? It is 8.56 MB more and a further doubling, but it
   changes the field `NUMERICS.md` §11 norms, so it needs §11 to say what a dual matrix's field is.
   `0014` (draft) question 25's sibling, and not asked there.
4. Should a twin ever be an *oracle* rather than a diagnostic? It is the only "oracle" that is
   `libm`-neutral by construction, which is `0032` open question 3 in a sharper form; it is also the
   candidate's own program, so domination against it is not a comparison with anyone else's
   algorithm.
5. Is there a second twin worth keeping — the three application forms `0036` measured, as subjects
   rather than as a branch that was reverted? Same mechanism, same validation shape.
