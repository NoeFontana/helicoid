# 0038: A program comparison is not a bar

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** the `W²` form twin; the rest proposes how two programs are compared, and corrects `0037`.

## Context

[`0037`](./0037-six-records-cannot-answer-a-one-variable-question.md) (draft) found that 66 of the 75
domination failures are maxima over six records and asked (its open question 6) whether
[`0006`](./0006-the-instrument-comes-first.md)'s bar needs a statistic that uses more than one
record of a small stratum. It proposed a paired *unanimity* rule and quoted its false-positive rate.

Measuring it changed the question. The bar is not the problem; using the bar to **choose between two
programs** is, and this repo has been doing that — `0036` (draft) ranked three application forms by
their domination-failure counts, and the first thing I did with the twin below was judge it by "0
failures fixed, 7 newly broken".

## Measurement

### A second twin, and a negative result

`seeded:w2-identity` is `seeded:correct` with one variable changed: the `b W²` term of
`J_l(φ) = I + aW + bW²` formed as `φφᵀ − θ²I` instead of as the matrix product `W·W`. `Q`'s own
`W·W` stays a product under both, so the variable is that one term. The two are the same matrix in
real arithmetic, which is the twin's control, since a form twin has no published column to reproduce
(`the_two_w2_forms_are_one_matrix_to_rounding`: they agree to a few ulp of `θ²` and differ in the
bits on 1 in 3 random directions). This answers `0037` open question 5: a form twin is worth keeping,
and this is the first.

The identity form is **worse**. Over the committed corpus it moves **125 of 751** rows, across all
nine `sen3_*` ids, and the direction is one-sided:

| | rows moved | identity worse | identity better |
|---|---|---|---|
| all nine ids | 125 | **104** | 21 |

An exact two-sided sign test over the 125 moved rows gives **p = 2.1 × 10⁻¹⁴**. Every one of the nine
ids is net worse on its own. The ratio of maxima is median 1.027, geometric mean 1.061, p90 1.254.
Where the identity wins it is near `θ = π` at large `ρ` and by small margins, with two exceptions
(`sen3_exp_n*` at `rho:1e4/theta=pi-1e-6`, 0.73× and 0.64×).

Why: the identity's diagonal entry is `φᵢ² − θ²`, one subtraction of like-sized quantities, where
`W·W`'s is `(0 − φⱼ²) + (−φₖ²)` — a sum of terms of one sign, with no cancellation. `0036` found the
identity *better* near π, but that was the **application** `J_l ρ` as
`r + a(φ×r) + b(φ(φ·r) − θ²r)`, where no matrix is formed and no diagonal exists. Different
operation, and the two results do not conflict.

### What the bar said about it

The domination bar's report of that same change: **0 failures fixed, 7 newly broken**, 75 to 82. In
isolation — 7 moved verdicts out of 751 rows — that reads as noise, and by `0037`'s resolution
measurement it *is* within the design noise of a six-record stratum. The bar could not have reached
p = 2.1 × 10⁻¹⁴, because it throws away the two things that produce it: the per-record pairing and
the pooling across strata and ids.

### Why the bar cannot be fixed by adding records

`sen3_jl_n1` at 64 records per stratum (`0037`'s regeneration), per-record errors for both forms.
Twenty-seven of the 52 strata are bit-identical on all 64 records — the form is not a variable there
at all. Twenty have at least 8 of 64 records differing and a net direction; on those:

| rule on `m` records | m = 6 | m = 12 | m = 24 |
|---|---|---|---|
| **per-stratum max** (today's bar) | 1.7 : 1 | 1.4 : 1 | **1.2 : 1** |
| per-record majority (paired sign) | 2.8 : 1 | 3.7 : 1 | **5.6 : 1** |
| max of the per-record ratios | 2.0 : 1 | 2.0 : 1 | 1.8 : 1 |
| worse on every record (unanimity) | silent 99.7% | silent 100% | silent 100% |

Odds of recovering the 64-record per-record direction, right : wrong, ignoring the draws it declines
to call.

**The paired majority converges and the max does not.** That is not noise: at the **full** 64
records, the max direction and the per-record direction **agree on 9 of 20 strata**. They are
different quantities. `theta:pi-1e-7` is worse under the identity on 29 of 64 records against 7
better, and its maximum is *lower* under the identity; both statements are true. The median stratum
has 5% of its 64 records within 10% of its maximum, so a maximum is an estimate from about three
records however many are drawn.

**So `0037`'s open question 6 was the wrong question**, and its arithmetic was half of one: a
unanimity rule does have a 1.6% false-positive rate, and it is silent on 99.7% of six-record draws,
so the 33× I quoted there was bought with all of the power. Corrected below.

## Decision

Proposed, not settled.

1. **`0006`'s domination bar stays on the max.** A maximum is the right *gate*: a consumer feels the
   worst case, and a bar that passes a program because it is typically fine is the failure mode
   `0006` was written against. Nothing here argues for changing it, and `0037` open question 6
   should be read as answered *no*.
2. **A program comparison is a separate instrument, and it is paired and pooled.** Choosing between
   two readings of a formula — which `0014` (draft) questions 19, 24 and 25 and `0036` all require —
   is a question about two programs on the same inputs, not about either program's worst case. It
   should be answered by the per-record differences of two twins over every stratum and id both
   answer, reported as a direction with a sign test, and never by counting domination failures.
3. **Counting domination failures to rank programs is retired.** `0036`'s ranking of the three
   application forms (36 / 16 / 19 failures) is a *bar* reading of a *comparison* question. It
   should be re-reported as a paired comparison before it decides anything, and `0036`'s open
   question 1 should not be settled on the failure counts as they stand.

## Rationale

The two questions differ in what they are statistics *of*. A domination failure is a property of one
stratum's worst record; a program choice is a property of the whole corpus. Pooling is legitimate for
the second and illegitimate for the first — `0006` forbids a mean precisely so that a bad stratum
cannot be averaged away, and that reasoning does not transfer to "which of these two assemblies
should ship", where every record is evidence.

The repo already contains the right pattern in the other half of the instrument. `0033` chose a
paired statistic for latency — bracketed control, paired bootstrap of a ratio, swap-invariant — and
gave the reason: a difference measured on the same machine at the same moment is the only one that
is about the code. The accuracy side kept a bare max for both jobs, and so has been answering
program questions with gate machinery.

This record's own negative result is the demonstration. The identity form is worse at
p = 2.1 × 10⁻¹⁴, which is not close; the bar's view of the same change was seven rows, which is not
readable. Had the bar been the only instrument, the honest conclusion would have been "no effect",
and a worse program would have been as defensible as the current one.

## Consequences

- **`seeded:w2-identity` is a retired hypothesis, recorded so it is not retried.** `W² = φφᵀ − θ²I`
  in `J_l`'s formation does not explain the 30 `sen3_j{l,r}_n1` domination failures; it makes the
  library measurably worse. The twin stays in the tree as the evidence.
- **The 30 `sen3_j{l,r}_n1` failures remain unexplained**, and two candidate causes are now
  eliminated by measurement rather than by argument: the transcendental library (`0037`) and the
  `W²` form (here). What has never been varied is `Q`'s word association (`NUMERICS.md` §5.3's
  `d_words = (W²X + XW²) − 3·WXW`, a difference of like-sized matrices) and the per-coefficient
  versus grouped switch (`0014` question 24).
- **`0037` open question 6 is answered, and its unanimity arithmetic is corrected here**, not in
  that record: a draft is not edited to match a later finding.
- A comparison instrument needs per-record errors of two subjects, which nothing in the tree
  reports today — the result CSV is one row per `(fn, stratum, precision, subject)` by `PHASE1.md`
  §5. That is the implementation cost of item 2 and the subject of open question 1.

## Implementation plan

1. The `W²` form twin, its control and its measurement — verified by
   `the_two_w2_forms_are_one_matrix_to_rounding` and by `just conformance --subject
   seeded:w2-identity`. **Landed.**
2. A paired comparison of two registered subjects: per-record differences over every shared
   stratum, a direction and a sign test per id and pooled, printed and never a bar. **Owed, blocked
   on open question 1.**
3. `0036`'s three application forms re-reported through it, replacing the failure counts.
   **Owed, after 2.**
4. `Q`'s word association as the next one-variable twin, since it is what is left.
   **Owed, after 2.**

## Open questions

1. **Where does a paired comparison live?** It needs per-record errors of two subjects, which
   `PHASE1.md` §5's result CSV does not carry. A second file, a flag on `conformance` that compares
   instead of scoring, or an `xtask compare` of its own — and whichever it is, §5 gains a sentence.
2. Should a comparison be allowed to pool across *ids*, or only across strata of one id? The
   p = 2.1 × 10⁻¹⁴ above pools nine ids that share a kernel, which is a real dependency between the
   rows it counts. Pooling per id and then combining is the conservative reading and gives the same
   direction nine times out of nine.
3. Does a comparison owe a *magnitude*, not only a direction? "Worse at p = 2e-14 by a geometric
   mean of 1.061" is a complete answer; a sign test alone is not, and a formula change should be
   worth stating in ULP.
4. Is there a stratum where the form choice should differ, as `0036` found for the application
   form? Two `sen3_exp_n*` strata improve by 0.64× to 0.73× under the identity, both at
   `rho:1e4/theta=pi-1e-6`. A per-regime choice is `0004`'s territory (a switch that is generated,
   never typed) and would need the sweep extended to assembly forms, which `0036` open question 1
   already weighs.
