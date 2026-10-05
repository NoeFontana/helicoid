# 0047: The second arm is admitted by agreement, not by the objective

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this PR — `coeffs::Switch` carries the second arm, `PHASE1.md` §6 states stage 2,
and both generated files are regenerated.

## Context

[`0039`](./0039-the-sweeps-grid-stops-below-its-own-optimum.md) is `ready` and its items 5, 6 and
10 are the **second switch**: a series arm's term count is set by the largest branch variable it
serves and paid by the smallest, so one arm makes every input pay the hardest one's price. Item 6
fixes the search's shape — stage 1 is the objective, unchanged; stage 2 is a cost search run
afterwards, because a shorter arm underneath a maximum the long arm already attains cannot lower
that maximum. Item 10 fixes the cost model: one term, one unit, weighted by the corpus's own
records.

Item 6 also fixes stage 2's **feasibility rule**, in four words: *"the cheapest arm that **holds**
stage 1's objective"*. Implementing it raised a question the record does not answer, because the
record never asks it.

**The objective is a maximum and a bar is not.** `Score::objective` is the largest error over every
`theta:*` record of one coefficient at one precision.
[`0006`](./0006-the-instrument-comes-first.md) reads domination and no-regress **per function, per
stratum, per precision**, on that stratum's own maximum. A prefix that holds the objective is free
to raise any record whose error is below it — which is every record but one. So the rule as written
admits a second arm that leaves the objective byte-identical and still moves, or breaks, a stratum
that the objective never looked at.

### How much it would move, measured

Both rules run on the committed corpus through `helicoid::__sweep`, at the switch and term count
stage 1 chose. The strict rule admits a prefix only where it agrees with the whole arm **to the
bit** — at every grid point and at every corpus record below the second switch, in the value and in
$\mathrm{d}/\mathrm{d}z$ alike. The loose rule is item 6's, literally. `worst ratio` is the largest
factor by which the admitted prefix raises a record's error above the whole arm's at that same
record, over the records it serves:

| coefficient | strict $m_0$, $z_0$ | loose $m_0$, $z_0$ | loose worst ratio |
|---|---|---|---|
| `k` f64 | 5, 6.98e-3 | 4, 5.83e-4 | 6.8x |
| `a` f64 | 5, 2.74e-3 | 5, 5.83e-3 | 43x |
| `b` f64 | 5, 2.84e-3 | 5, 8.66e-3 | 40x |
| `c` f64 | 6, 1.24e-2 | 5, 8.66e-3 | **13370x** |
| `d` f64 | 5, 2.64e-3 | 5, 1.24e-2 | 146x |
| `e` f64 | 5, 3.52e-3 | 4, 6.49e-4 | **5084x** |
| `cos θ/2` f64 | 5, 5.23e-3 | 4, 3.65e-4 | 171x |
| `r` f64 | 5, 5.05e-5 | 5, 2.29e-4 | 1709x |
| `k` f32 | 3, 3.52e-3 | 3, 1.78e-2 | 4.2x |
| `a` f32 | 3, 1.00e-3 | 3, 8.98e-3 | 31x |
| `b` f32 | 3, 1.38e-3 | 3, 2.29e-2 | 53x |
| `c` f32 | 4, 4.53e-2 | 4, 1.60e-1 | 36x |
| `d` f32 | 3, 1.38e-3 | 3, 3.79e-2 | 145x |
| `e` f32 | 3, 2.46e-3 | 3, 3.40e-2 | **50224x** |
| `cos θ/2` f32 | 3, 1.91e-3 | 3, 8.98e-3 | 145x |
| `r` f32 | 4, 1.49e-3 | 4, 7.77e-3 | 33x |

And what the two cost, in corpus-weighted term units over the records the series arm serves:

| | one arm | strict | loose |
|---|---|---|---|
| binary64 catalogue | 172427 | 92264 (**1.87x**) | 87196 (1.98x) |
| binary32 catalogue | 119547 | 60082 (**1.99x**) | 54990 (2.17x) |

**The loose rule is 5.8% cheaper at binary64 and 9.3% at binary32, and it multiplies individual
corpus rows by up to four orders of magnitude.** `e` at binary32 is the extreme: a prefix that
holds a 15.6 `u` objective is 50224x worse than the whole arm at some record below its switch,
because the objective is attained near $\pi$ where the *exact* arm runs and says nothing whatever
about the series arm at $z = 10^{-2}$.

### What the strict rule costs, which is the question the ratios above do not answer

Reach, not terms. The two rules pick the same $m_0$ on eleven of sixteen rows; where they differ it
is by one term. What differs is $z_0$, and by less than a decade on every row: a prefix's truncation
is $O(z^{m_0})$, so admitting it up to a tolerance of 94 `u` instead of half an ulp moves $z_0$ by
$94^{1/m_0}$ — a factor of 1.3 at $m_0 = 16$, 1.6 at $m_0 = 8$ — and the records are log-uniform
across seventeen decades (`error-analysis.md` EA.21). A fraction of a decade of a seventeen-decade
sample is the whole of the 5.8%.

## Decision

1. **Stage 2 admits a prefix where it agrees with the whole arm to the bit, not where it holds the
   objective.** `0039` item 6's feasibility rule is amended: the second switch is the largest grid
   point such that at **every grid point and every corpus record below it** the prefix and the whole
   arm give the same bits, in the value and in $\mathrm{d}/\mathrm{d}z$. Item 6's rule is implied —
   equal bits have equal errors — so this is strictly inside the feasible set item 6 describes, and
   it is the part of that set in which the second arm is a latency change and nothing else.
2. **The cost model is item 10's, unchanged**: one term, one unit, weighted by the corpus's own
   records, minimised over the prefix length, `O(m · n)` over a prefix maximum. A *timed* model
   would make `generated.rs` a function of the host and break D16.
3. **A second arm is a switch and a shorter prefix, or neither.** The emitter refuses a row with one
   and not the other. "Neither" is a second switch of `0`, which no branch variable is below, so a
   coefficient with no cheaper prefix needs no special case in the kernel.
4. **The agreement is checked at the grid, in the library, beside the constants it is about.**
   `coeffs::tests::second_arm_is_the_first` reads the grid out of
   `conformance/corpus/coeff_switch_ref.jsonl` — one record per grid point, so its keys *are* the
   grid — and compares the two arms at every point below the second switch. It needs no reference
   and no corpus run: a prefix that reaches further than its agreement does fails it.
5. **Between those points the two arms can part by one place, and that is measured, not assumed.**
   The rounding of each Horner step jitters the crossing, so a sample a few grid steps below the
   second switch can differ where the grid points around it do not. Worst over a dense log-uniform
   sample: **1.875 `u` at binary64** (`b`) and **3.000 `u` at binary32** (`k`), against a catalogue
   whose gentlest objective is 1.747 `u` and whose worst is 94.5. `coeffs::tests` carries it as a
   no-regress bound at the measurement, not twice it, because a prefix reaching past its agreement
   is the one regression the bound exists to catch.
6. **The second arm is `helicoid::coeffs`' and not the seeded kernel's.** `xtask/src/seeded` exists
   to validate the search (`PHASE1.md` §10's planted defects) and its generated file keeps one arm.
   A seeded twin of a prefix that is bit-identical by construction would validate nothing, and
   `0004` item 4's "the code measured is the code shipped" is about the `helicoid` target, which is
   swept through `__sweep` as before.
7. **A call-site group reads *one* second switch, its members' smallest.** Measured, not assumed:
   a mask and a branch per member costs one comparison and one branch on **every** call, including
   every call above the second switch where no member takes a short arm, and at near-identity
   `theta` that overhead **exceeded the terms it saved** on all three multi-member groups at
   binary32 — `exp_coeffs` 1.17x, `jr_coeffs` 1.15x and `q_coeffs` 1.15x *slower* than the one-arm
   kernel, each past its own concurrently measured A/A floor in all six bootstrap pairs — while the
   two one-member groups won 3.8x and 1.6x. Above the second switch it measured 1.20x to 1.28x
   slower, scaling with the group's size. Reading one switch for the group puts the short arms one
   comparison from the entry point and leaves the rest of `grouped` byte-identical to what it was.

   This does **not** cost accuracy: below the group's switch every member is below its own, so each
   takes the prefix the sweep admitted; between the group's and a member's own the member takes the
   *whole* arm, which is the arm the prefix was admitted to agree with there. Both readings give the
   same bits at every grid point and every corpus record, which is why the byte comparison in
   *Consequences* holds for either. The groups' smallest and largest second switches are within a
   factor of 1.33 (`exp_coeffs`), 1.04 (`jr_coeffs`) and 1.33 (`q_coeffs`) of one another, so the
   reach given up is a rounding of a seventeen-decade sample.

   `groups_are_their_members_own_arm` is therefore a statement about the group and not the
   coefficient for the second switch, and about the coefficient for the first: a member is its own
   exact arm above its own switch whatever it is grouped with, and its own *prefix* below its
   group's. The test carries both, with the group memberships of `PHASE3.md` §3 spelled out.

   The two one-member call sites, `jr_inv_coeff` and `log_ratio`, are written the same way for the
   same reason, and `log_ratio` is where the shape mattered most: with the prefix selected *inside*
   the series arm it was **1.16x slower** at `theta = 0.5` — where `s = 0.065` takes the whole arm —
   at both precisions, past a floor of 0.0029 in all six pairs, even though that spelling also saved
   a division. One comparison on the way in buys the prefix; a comparison under the series arm buys
   it and pays for two Horner bodies in one block. `SO3::log` is the consumer, and `r`'s two switches in
   `s = tan^2(theta/2)` are `5.05e-5` and `0.0898`, so that was a regression over the whole band
   `theta ∈ [0.0142, 0.583]` — every rotation between a degree and a third of a radian.

## Rationale

The alternative was to take item 6 as written and then check what it broke — which is the loop this
project has a rule against (D7, `0006`): the instrument before the algorithm. The instrument here is
free, because bit-identity needs no reference: two prefixes of one Horner compared at one argument.
Measuring both rules before shipping either is what turned a four-word phrase into a 5.8%-against-
13370x trade, and that trade is not close.

Loosening the bars to fit the loose rule was never a candidate:
[`0038`](./0038-a-program-comparison-is-not-a-bar.md) decided that the maximum stays the gate, and
[`0046`](./0046-explained-by-record-needs-a-record-to-point-at.md) decided that an exception is a
row with a record behind it, not a tolerance. A second arm that needed either would be paying for
latency with the project's only accuracy guarantee.

Choosing the measurement rather than twice it in item 5 departs from `laws::Bounds`' convention
deliberately. A `Bounds` row guards an algorithm whose error is a continuous function of its input,
where a factor of two absorbs a different host's `libm`. This bound guards a *selection rule*, where
the quantity measured is the rounding of one Horner step and a regression means the generated
`short_terms` reaches past the agreement the sweep verified. Two is slack enough to miss that.

## What it bought, measured

`cargo xtask bench-gate --bench coeffs --only near-identity-7.5e-8 --against`, the baseline built in
this tree from `main` (a build directory moves a row 2-3%, `PHASE3.md` §0.0), three replicates and a
paired bootstrap, point ratio with each benchmark's concurrently measured A/A floor:

| group | `f64` | floor | `f32` | floor |
|---|---|---|---|---|
| `jr_inv_coeff` | **0.2605** | 0.0038 | **0.2353** | 0.0068 |
| `jr_coeffs` | 0.7476 | 0.0462 | **0.5282** | 0.0131 |
| `exp_coeffs` | 0.8760 | 0.0050 | **0.5525** | 0.0537 |
| `log_ratio` | 0.6365 | 0.0043 | 0.5975 | 0.0052 |
| `q_coeffs` | 0.6585 | 0.0081 | 0.7826 | 0.0389 |

Ten of ten faster — 1.14x to **4.25x** — every effect outside its own floor. `jr_inv_coeff` is the
extreme because `c` is the one coefficient that takes the term cap: sixteen terms down to six at
binary64 and thirteen down to four at binary32.

### What it costs on the arms that do not take it

A second switch is one comparison a kernel did not make, and every call that does **not** reach the
short arm pays it. Same baseline, same protocol, at `near-pi` — where `theta = pi - 1e-6` puts every
one of the five groups on its **exact** arm, so this is the extra comparison and nothing else:

| group | `f64` | floor | `f32` | floor |
|---|---|---|---|---|
| `jr_inv_coeff` | 1.0014 | 0.0041 | 0.9997 | 0.0218 |
| `q_coeffs` | 1.0190 | 0.0068 | 1.0081 | 0.0069 |
| `jr_coeffs` | 1.0204 | 0.0021 | 1.0350 | 0.0101 |
| `exp_coeffs` | 1.0275 | 0.0066 | 1.0269 | 0.0116 |
| `log_ratio` | 1.1288 | 0.0785 | **1.3036** | 0.0512 |

**Four of the five cost 1% to 3.5%**, which is one comparison against a `sqrt` and a `sin_cos`, and
is the price of the 1.14x to 4.25x above. `jr_inv_coeff` costs nothing measurable.

**`log_ratio` costs 1.37 ns, and on its exact arm that is 12% at binary64 and 30% at binary32.**
The mechanism is not code size: adding `#[inline]` to all five entry points and to `grouped` was
measured and changed **nothing** (`log_ratio` f64 10.97 -> 12.34 ns with the attribute, 1.125,
against 1.1288 without), so the attribute is not shipped. What is left is the shape: `r`'s masks are
`w > 0 and s < switch`, so a second switch adds a comparison *and* an `and` to a chain that already
waits on the division forming `s`, in front of a body that is one `sqrt` and one `atan2`. 1.37 ns is
about four cycles, which is what that costs.

The alternative is to test the **first** switch first and the second inside the series arm, which
makes the exact arm byte-identical to the one-arm kernel and moves the 1.37 ns onto the whole series
arm instead. That was measured once, in the shape this record's item 7 rejected, at **1.16x** — on a
15-term Horner the same absolute cost is a larger ratio than on a `sqrt` plus an `atan2`. So the
ordering shipped puts the cost on the bigger body, which is the better trade at binary64 and is
*not* obviously the better one at binary32, where the exact arm is cheaper. Measuring that swap at
binary32 is **Further work**, with the numbers above as its baseline.

What the ordering buys, against `r`'s three regions in `theta` (its switches are `s = 5.05e-5` and
`s = 0.0898` in `s = tan^2(theta/2)`): the short arm serves `theta < 0.0142`, the whole arm
`[0.0142, 0.582)`, the exact arm the remaining 81.5% of `[0, pi]`. A uniform `theta` mostly reaches
the exact arm; `tf_tree`'s kilohertz traffic between consecutive poses is in the first region, which
is the one `PHASE3.md` §11 names first and where this measures **1.67x faster**.

## Consequences

- **Every conformance result is unchanged, byte for byte.** Verified directly: all nineteen files
  under `conformance/results/` are identical to the committed ones apart from the git revision each
  run stamps. So no envelope verdict, no baseline row, no `branch_continuity` bound and no
  `at_switch_*` column moves, and `0046`'s exception table is not touched.
- `conformance/sweeps/thresholds.csv` gains six columns — `short_terms`, `short_switch_bits`,
  `short_switch_z`, `short_grid_index`, `cost_one`, `cost_two` — and **not one of its existing cells
  changes**, which is stage 1's independence from stage 2, checked rather than asserted.
- `coeffs::Switch` gains `short_below` and `short_terms`, and `Switch::first` two arguments. The
  struct stays one object per coefficient per precision (`0039` item 10): nothing about the shape of
  `generated.rs` is hand-written.
- A group's branch structure gains one level, and **grouping stays a cost and never a value**: a
  member takes the arm its own two switches select, whichever members share a branch, which
  `coeffs::tests::groups_are_their_members_own_arm` pins at both precisions.
- `PHASE1.md` §6 is NORMATIVE and is edited: the search is two stages, and the second one's rule is
  item 1 above.

## Implementation plan

1. `coeffs::Switch` and `coeffs::kernel` carry the second arm; `search::second` and
   `measure::grid_arms` choose it; the CSV and both emitters report it; `PHASE1.md` §6 states it —
   verified by `second_arm_is_the_first` at both precisions, `the_two_series_arms_part_by_at_most_an_ulp`,
   `the_second_arm_is_the_cheapest_prefix_that_agrees_to_the_bit`,
   `a_second_arm_emits_both_switches_and_half_of_one_is_refused`, `just thresholds-check`, and a
   byte comparison of `conformance/results/` across the change.
2. The latency it was for — `cargo xtask bench-gate --bench coeffs --against` against a baseline
   built in this tree, with the ratio per fixture recorded in the PR and in `PHASE3.md` §11. This
   step is what found decision 7: the first reading shipped a 1.15x to 1.17x regression at the one
   stratum `PHASE3.md` §11 names first, and no amount of reasoning about term counts would have
   said so.

## Further work

- **The groups' own composition is measurable without risking a bit, and is probably not worth it.**
  `SEn3::jr` evaluates `b` twice — once in `jr_coeffs` through `SO3::jr`, once in `q_coeffs` — with a
  second `sqrt` and a second `sin_cos` beside it, and `q_coeffs` evaluates `b` and `d`'s Horners
  *and* its exact arm for $z \in [8.98, 9.65)$. A member's value is `small_i ? horner_i : exact_i`
  whatever it is grouped with (decision 7's test), so a merged `(a, b, d, e)` group is **bit-identical**
  — which corrects `PHASE3.md` §0.0's bench row, where the same proposal was retired partly for
  "moving the arm boundaries". What that row measured still stands and is the reason to leave it:
  `se3/jr`'s near-identity row, with every coefficient on a series arm and no `sin_cos` at all, is
  only 30 ns under its generic row, so the transcendentals are **~7% of that row** and Barfoot's `Q`
  is the cost. A `sqrt` and a `sin_cos` saved is worth about 7% of one method.
- **A switch shared by a call-site group** — one swept switch per group rather than per coefficient —
  is `0015` (draft) NU.5 option (a) and *would* move values. It is also the standing lead for the
  four `theta ~ 1` regressions `0039` left (CO.18). Separate record.
- **`log_ratio`'s switch ordering at binary32.** Testing the first switch first and the second
  inside the series arm makes the exact arm byte-identical to the one-arm kernel and moves 1.37 ns
  onto the whole series arm. At binary64 that is the worse trade (a 15-term Horner is a smaller body
  than a `sqrt` plus an `atan2`); at binary32, where the exact arm measured 30% sensitive to 1.37 ns,
  it may be the better one. One `bench-gate --only` run each at `near-pi` and at the series fixture
  decides it. Not swapped on a guess.
- **`cargo xtask bench-gate --only` is unreliable on a repeated run.** Three of six invocations in
  this PR's measurement failed with `no 'aa-1/sample.json' under target/criterion` or
  `the filter matched 2 benchmarks in 'aa-3', not 1`, on a cleaned `target/criterion` and an
  argument that had just worked. The full-suite run does not fail that way. That is why the series
  fixture's numbers here are the pre-item-7 shape's and not the shipped one's, and it is a defect in
  the gate, not in the subject.
- **The Horner's dependency chain is what is left.** A five-term arm is five dependent
  multiply-adds, about 40 cycles of latency that only overlap because consecutive calls are
  independent; a shorter chain for the same terms is an *association* change, so it is a
  `NUMERICS.md` §4 edit (CO.9 states Horner), a new sweep and its own record — the error constants
  are the association's, not the polynomial's. Not speculated on here; noted so the next person
  measuring the kernel knows where the remaining cost is.
