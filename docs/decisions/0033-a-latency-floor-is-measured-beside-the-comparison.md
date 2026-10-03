# 0033: A latency floor is measured beside the comparison

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** `cargo xtask bench-gate` (`xtask/src/bench/`) — the estimator choices below are
implemented; the bar change `0006` would need, and `--against <ref>`, are not.

## Context

[`PHASE1.md`](../PHASE1.md) §9 fixes the gate's *shape* and deliberately not its estimator:
"interleaved baseline/candidate runs, paired bootstrap 95% CI of the ratio; fail when the whole CI
lies above $1 + \delta$, with $\delta$ the host's A/A noise floor measured by `bench-gate --aa` and
recorded in `baseline/HOST.md`. A measured floor, not an asserted percentage."

It names the interval, the pairing, the bootstrap and the comparison. It does not name the location
statistic, how a pair is formed, how many pairs a verdict takes, or whether $\delta$ is one number
or sixty. Every one of those turned out to decide whether the gate works at all, and each was
settled by a measurement that refuted the first guess. This record is where those measurements
live, because the code can only carry the conclusions.

[`0006`](./0006-the-instrument-comes-first.md) does name a statistic, for the accuracy bars:
domination and no-regress are judged "per function, per stratum, per precision, on the **max** —
never a mean". The question this record opens is whether that sentence governs §9 too.

## Measurement

One host, `AMD EPYC-Milan Processor`, 8 logical CPUs, a shared VM whose effective CPU allocation
changed between runs (`available_parallelism` read 1 at one bless and 8 at another). 60 benchmarks
(`crates/helicoid/benches/coeffs.rs`, four coefficient functions × five θ strata × two precisions),
criterion 0.5.1, 100 samples per window, 1024 bootstrap resamples, a fixed seed. $\delta$ per
benchmark is $\max(\text{high} - 1,\ 1/\text{low} - 1)$ over the paired 95% CI of the ratio.

| estimator question | first guess | what it measured | what the measurement forced |
|---|---|---|---|
| location statistic | mean of per-iteration times | worst $\delta$ **0.7172** | median: **0.0190** on the same windows |
| pair formation | index-paired across two passes of 60 | two measurements ~2 min apart spread to **90%** | adjacent in time: four runs in a row spread **0.6%** |
| pairs per verdict | one | identical bytes put **6 of 60** past their floor, one at **1.41×** | three alternations, the least disturbed reported |
| combining the three | pool all three windows | one benchmark read **0.5884** where ten consecutive windows of it spread **2.38%** | select, do not pool: **0.0022** |
| which way $\delta$ points | $\text{high} - 1$ | the four noisiest pairs (22–43% apart) got the run's **tightest** floors | swap-invariant $\max(\text{high}-1, 1/\text{low}-1)$ |
| invocation path | baseline direct, candidate under `cargo` | **3 of 60** identical-code failures at 1.007–1.018 | one prebuilt binary, exec'd directly on both sides |
| $\delta$ per benchmark or per host | per benchmark | **8 of 60** identical-code failures; one benchmark's $\delta$ moved 0.0019 → 0.0190 | one floor for the host |
| where $\delta$ comes from | `HOST.md`, as §9 says | worst $\delta$ 0.0161 one session, **0.6864** the next | a control measured beside the comparison |
| how many pairs must agree | one, as §9 says | every false verdict in this module's history was one disturbed window | all `REPLICATES` × 2 of them |

The last row is the one that bears on `0006`. Two A/A runs of all 60 benchmarks, same protocol,
same binary against itself:

| | run 1 | run 2 |
|---|---|---|
| min $\delta$ | 0.0008 | 0.0009 |
| median $\delta$ | 0.0029 | 0.0030 |
| max $\delta$ | 0.0161 | 0.0201 |

The *distribution* reproduces to a digit. Individual benchmarks do not: `exp_coeffs_straddle-0.9`
read 0.0019 and then 0.0190, `jr_coeffs_series-0.5` 0.0022 then 0.0187. The max over 60 — the only
summary the gate's own question ("does any benchmark fail?") asks for — moved by a quarter of
itself.

Contamination explains the asymmetry. A disturbed sample is slow and never fast, and so is a
disturbed window, so the noise is one-sided at both levels and $\delta$'s upper tail is set by rare
large events. A median over 100 samples and a selection over 3 windows together suppress those
events but cannot abolish them; what survives lands on whichever benchmarks were unlucky, and never
the same ones twice.

### What the positive control exposed, after the first fix

The bracketed control passed its positive control once (60 benchmarks, 0 failures) and then failed
it: **2 of 60** on a later session, both `log_ratio` at binary32, each claiming at least 1.0109 and
1.0155 in *every one of six pairs* while its own concurrent control was quiet (floors 0.0040 and
0.0033). Identical binary on both sides.

The first hypothesis was a position effect — the candidate always runs second of three, so anything
specific to the middle window is invisible to a control measured between windows one and three. The
data refuses it. Over 60 benchmarks the median candidate-to-baseline ratio is **0.9994** under
`criterion` 0.5.1 and **0.9989** under 0.8.2, with 24 and 18 of 60 above 1.0: no systematic offset,
in either direction, under either version.

What remains is an independence failure. With the replicate loop *inside* the benchmark loop, a
benchmark's three triplets ran consecutively, so all six of its pairs were taken inside about 30 s.
The all-pairs rule survives one disturbed window; it does not survive a disturbance that outlasts the
whole span, and such a disturbance leaves the control quiet because it covers the control's two
windows equally. So the rule as first stated — "an isolated disturbance cannot satisfy it" — was
true of a window and false of a span.

The fix costs nothing: run the replicate loop **outermost**, a full pass of all 60 benchmarks per
replicate. A benchmark's triplets are then about 20 min apart here, the windows within a triplet stay
contiguous (which the adjacency measurement above requires), and a 30-second disturbance can spoil at
most one triplet of three. The same total windows are measured.

### The measurement that decides the floor

Four A/A runs of the same protocol, same binary, same machine, same flags, over one day:

| | run 1 | run 2 | run 3 | run 4 |
|---|---|---|---|---|
| min $\delta$ | 0.0008 | 0.0009 | 0.0012 | — |
| median $\delta$ | 0.0029 | 0.0030 | 0.0089 | 0.0039 |
| p75 $\delta$ | 0.0036 | — | 0.0832 | — |
| p95 $\delta$ | — | — | 0.5776 | — |
| max $\delta$ | 0.0161 | 0.0201 | **0.6864** | 0.0650 |
| above 2% | 0 of 60 | 0 of 60 | 23 of 60 | — |
| above 5% | 0 of 60 | 0 of 60 | 18 of 60 | — |

Run 3 is not one outlier: the whole distribution moved, by a factor of 40 at the top. Nothing in
the protocol changed between any of them, so what moved is the machine's state, and a shared VM
does not announce it (`available_parallelism` read 1 at one bless and 8 at another, which is the
same fact from a different angle). Run 4, measured after the gate was rebuilt, is a fourth value
again — 0.0650, four times run 2 and a tenth of run 3. The median is the stable thing across all
four; the maximum, which is the only summary the gate's question asks for, is not.

This refutes the stored floor itself, not a choice of statistic for it. **A floor from the quiet
session fails the identical binary against itself; a floor from the noisy session would pass a 50%
regression.** No summary of run 1 and 2 predicts run 3, and the accumulating maximum proposed to
cover that — the first thing tried, and landed for an hour — is refuted by run 3 in the other
direction: one noisy session pins the allowance at 68% permanently, for every candidate, forever.

The gate's time budget is what makes a concurrent control affordable: a benchmark's pair takes about
5 s, so a third window per replicate costs 50% more wall clock and buys a floor measured in the
same minute as the thing it licenses.

## Decision

Proposed, not settled; items 1 and 2 are implemented because a gate that reads a stale floor is
worse than no gate, and the §9 reading in item 2 is the one question this record most needs
answered.

1. **`0006`'s "on the max, never a mean" governs accuracy and not latency.** For an accuracy bar
   the max is the right statistic *because the measurement is deterministic*: D16 makes a stratum's
   error a function of the input bits, the corpus is committed, and re-running reproduces the max
   exactly. A max over reproducible quantities is a bound. A max over random ones is an extreme
   value whose sampling error grows with the number of draws — which is why 60 latency benchmarks
   produce a max that drifts while their median does not, and why 180 draws of it drift further.
   `0006`'s sentence should name the property it depends on (reproducibility), not the bar.
2. **$\delta$ is measured beside the comparison, in the same invocation, and `HOST.md` becomes a
   log.** `--against` runs each benchmark as `b₁ c b₂` — baseline, candidate, baseline — three
   times. The control pair $(b_1, b_2)$ spans *more* elapsed time than either candidate pair, so
   the floor it yields is conservative for the comparison it licenses; the candidate's two pairs
   are judged against the noisiest control of the three replicates. This is a **deviation from
   §9's letter**, which says $\delta$ is "recorded in `baseline/HOST.md`" — the file is still
   written by `--aa` and is still how a human learns whether the machine is worth using, but the
   gate does not read it. §9's *substance* — "a measured floor, not an asserted percentage" — is
   what the concurrent control serves, and the stored floor as literally specified is the one thing
   the measurement rules out.
3. **A regression must appear in every pair, and the pairs must be spread over the run.** Each
   benchmark yields `REPLICATES` × 2 candidate pairs, and the gate fails only when every one puts
   its whole CI above $1 + \text{floor}$. §9's rule is the per-pair test; requiring all of them is
   what a *short* disturbance cannot satisfy. That is not sufficient on its own, and the correction
   is in the Measurement below: the **replicate loop is the outer one**, so a benchmark's three
   triplets are a whole pass apart rather than consecutive, at no extra cost.
4. **The gate reports its own resolution, every run.** A pass under a 0.3% floor and a pass under a
   30% floor are not the same statement, so `--against` prints the median and worst concurrent
   floor beside the verdict. An effect smaller than a benchmark's own floor is reported as "below
   the floor of this run", never as "no regression".

## Rationale

The alternative to one floor is the per-benchmark floor, and it is the better instrument in
principle: `jr_coeffs` at near-identity is genuinely noisier than `q_coeffs` at θ = 1, and one
scalar lets the worst benchmark buy allowance for the other 59. It was implemented first for that
reason. It fails for a reason that is not about principle: a single pair estimates one benchmark's
floor with a sampling error the size of the estimate, so 60 stored floors are 60 under-powered
guesses and the verdict is dominated by which ones guessed low. The bracketed control keeps the
per-benchmark resolution that made it attractive and drops the storage that made it wrong.

The alternative to a concurrent control is a margin — the stored floor times some factor. That
factor would be an asserted percentage wearing a measurement's clothes, which is the exact failure
§9 names, and run 3 shows no finite factor works: 0.0161 and 0.6864 are both honest measurements of
the same protocol.

The alternative to selecting the least-disturbed pair is pooling, which is what a bootstrap over all
the data would normally want. It is refused by measurement, not taste: a disturbance long enough to
cover two of a side's three windows — about 15 s per benchmark here — is past a median's 50%
breakdown point, and pooling then reports the disturbance as the signal.

The alternative to bracketing (`b₁ c b₂`) is `b₁ b₂ c`, which is cheaper to reason about and wrong:
it measures the floor over a *shorter* span than the comparison it licenses, so the drift the
candidate suffers is not in the control. The order is the whole point.

## Consequences

- `0006` §"the bars" needs the reproducibility qualifier, or it reads as mandating a statistic for
  §9 that §9's own measurements refute. That edit is this record's to propose and not to make.
- §9's `baseline/HOST.md` row changes meaning: it is a log of what the host has shown, not the
  gate's allowance. A spec edit is owed if this record goes `ready`.
- `--against` costs 9 windows per benchmark instead of 6: about 20 minutes for the 60 coefficient
  benchmarks on this host. `PHASE3.md` §11's group benches will add to that, and the budget is the
  reason `REPLICATES` is 3 and not 5.
- The gate's sensitivity is a property of the run, not of the repository, and is printed. On a noisy
  session the honest verdict is that nothing below (say) 8% is resolvable — not that the candidate
  passed.
- **The gate detects reproducible latency differences between two binaries; it does not attribute
  them.** The negative control measured 1.6% to 4.1% on four benchmarks whose machine code is
  byte-identical across the two binaries and merely relocated by a 64-byte `.text` shift. Source
  locality does not survive to machine level, so a failure on a benchmark that cannot route through
  the change is evidence of layout, and the way to confirm it is a symbol diff (`nm --defined-only
  -S`, compare size against address). Reporting such a row as "the change made `exp_coeffs` slower"
  would be false; reporting it as "these two binaries differ here" is what was measured.
- A `#[repr(align)]`-style pinning of the benchmark loops, or a layout-randomising harness that
  measures the *distribution* over link orders, is the only way to separate the two. Neither is
  specified and neither is cheap; `0013`'s "owes a measurement" posture applies.

## Implementation plan

1. The estimators of the table above — median location statistic, adjacent pairing, three
   replicates, least-disturbed selection, swap-invariant $\delta$, one prebuilt binary exec'd
   directly — verified by `xtask/src/bench/`'s tests and by `HOST.md`'s recorded runs. **Landed.**
2. The bracketed concurrent control and the every-pair rule — verified by
   `a_regression_must_appear_in_every_pair` and `the_report_orders_by_the_middle_pair`. **Landed.**
2a. The replicate loop outermost, so a benchmark's triplets are a pass apart. **Landed**; its own
   re-validation is step 3's, re-run.
3. Positive control: the identical binary as its own baseline, through the whole `--against` path.
   **Landed: 60 benchmarks, 0 failures**, on a session whose worst concurrent floor was 0.5909 —
   the condition under which every stored floor produced false verdicts.
4. Negative control: four extra `sqrt` round trips in `coeffs::kernel::exact_a`. **Landed**, and it
   corrected its own prediction twice. The prediction was six `jr_coeffs` rows; it is **five**,
   because `A_F32` and `B_F32` both switch at θ² < 1 so `straddle-0.9` is not a straddle at
   binary32 and no exact arm runs there. All five fired, claiming at least 1.38× to 1.57× in every
   one of six pairs. Four *more* rows fired, 1.016× to 1.041×, all `exp_coeffs`, which cannot reach
   `exact_a` — one caller, `jr_coeffs`. The cause is measured, not guessed: between the two
   binaries only **2 of 3192** text symbols changed size (the two `Real::branch` monomorphisations
   that inline `exact_a`, +32 bytes each), `.text` grew 64 bytes, and **2997 symbols kept an
   identical size at a different address**. Byte-identical code, relocated. The remaining 51 rows
   did not fire.
5. `--against <ref>`: build the baseline binary from a git worktree at a ref, so `just bench-check`
   can name a commit instead of a path. Today `--against` takes a prebuilt binary. **Owed.**
6. `0006`'s bar sentence gains the reproducibility qualifier, and `PHASE1.md` §9's `HOST.md` row
   says log rather than allowance — both verified by this record moving to `ready`. **Owed.**

## Open questions

1. **Does §9 accept a concurrent floor?** §9 says $\delta$ is "measured by `bench-gate --aa` and
   recorded in `baseline/HOST.md`", and the gate as built measures its own. The measurement says
   the stored floor cannot work on this host; the spec says where $\delta$ lives. This is the
   question that blocks `ready`, and it is the owner's.
2. Does `0006`'s max survive anywhere in §9, or is the median the only defensible summary of a
   latency distribution? The gate needs no location statistic *across* benchmarks — it fails per
   benchmark — so the question is only about what gets reported and read.
3. Is there a floor at all on a machine this unstable, or is the right answer a precondition: refuse
   to render a verdict when the concurrent control is above some measured level? That level would
   be a percentage, and asserting one is what §9 refuses, so the honest form is "report the
   resolution and let the reader decide" — which is what item 4 does. Whether `just bench-check`
   may therefore *pass* on a session that could not have detected anything is unresolved.
4. Does the floor belong per *stratum family* rather than per benchmark? Five θ strata × four
   functions is 20 cells per precision, and 3 pairs each is 180 draws per cell per run, which may
   estimate a cell's extreme where it cannot estimate a benchmark's. Not measured.
5. §9's `--aa` measures *this* set of 60. `PHASE3.md` §11 adds the group benches; whether their
   noise is the same is unmeasured, and the log makes that visible rather than assumed.
