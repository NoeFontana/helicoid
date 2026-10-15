# 0091: A twin is one draw to the swap test

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

[`0088`](./0088-the-domination-bar-must-tell-worse-from-equal.md)'s swap test assumes each record's
pair of errors is exchangeable on its own. The corpus breaks that by design:
- **Twins.** `gen/so3.py` and `gen/sen3.py` write every quaternion and then its negative, the half
  a `Log` without its flip gets wrong. `q` and `−q` are one rotation, and every program answers
  both through one arithmetic once it canonicalizes, so the two records carry one error twice.
- **Where.** Twins occur in `so3_log` (29 of 30 strata) and `sen3_log_n1`, `_n2`, `_n3` (53 of
  54 each), and nowhere else. The generator also writes some inputs more than once: in
  `real_acos`, `x:near±1` holds 37 distinct inputs in 64 records, 22 and 20 at `@f32`; in
  `real_sin_cos`, `x:near-k-pi/2` holds 61. A deterministic program answers a repeat as it answers
  the first, so a repeat is no more evidence than a twin.
- **What it does to the test.** Swapping one twin and not the other is no arrangement the null
  hypothesis allows. Counting both roughly squares the p-value, and lets a 6-draw stratum look as
  if it had 12.

Main's envelope, after `0089` and `0090`, failed one stratum:
- **`sen3_log_n2/theta:pi-1e-2` at `f64`:** 1.747 `u` against GTSAM's 0.794, resolved at
  `p = 9.8e-4`.
- **Its 12 records are 6 rotations, each twinned.** With the twins swapped together, `p = 3.1e-2`,
  and 6 draws cannot reach `α = 0.001` at all (the floor is `2⁻⁶`).
- **On 64 rotations the loss reverses.** On the deep tier, `helicoid` reads 2.192 `u` against
  GTSAM's 2.696, and at `_n3` 2.152 against 2.901. GTSAM's 0.794 was the luck of 6 rotations.
- **Recomputing every excess with twins swapped together**, it is the only verdict that changes on
  the committed tier.

Two arithmetic alternatives for `SEn3::log`'s translation were measured on the deep tier and lost:
- **`ρ = βx + c(φ·x)φ − ½φ × x`**, with `β = (θ/2)cot(θ/2)`, which does not cancel near `π`: 10
  resolved losses against main's 7 over the three ids, summed maxima 242 against 220.
- **GTSAM's matrix `J_l⁻¹ x`:** 12 and 245.

`SEn3::log` keeps its arithmetic.

## Decision

1. **The swap test's unit is a draw, not a record.** Two records of a stratum are one draw when
   every input is equal, a quaternion input up to its sign.
   - **Which inputs are quaternions:** they are declared once, as `corpus::QUATERNION_INPUTS`
     (`q`, `q0`, `q1`), and a test holds the list to the corpus.
   - **Pairing:** each draw is one pair, each side's maximum over the draw's records. The draws
     swap whole, so the maxima are what a swap moves.
2. **The conformance run writes the draw.** It reads the corpus's validated inputs and writes the
   draw beside each record's error, as a sixth cell of the records file: the smallest id of the
   records one draw with it.
   - **No stale grouping:** the grouping is the one the errors were scored under, so a corpus drawn
     again after a run cannot regroup that run's records.
   - **The envelope reads the corpus nowhere:** it parses no input and pays nothing.
   - **Old files read as before:** a records file in `0088`'s five-cell format reads each record
     as its own draw, which is `0088` exactly.
3. **Resolvability counts the draws the test was over.** `0088` decision 2's floor is computed
   over the draws both subjects scored, so a twinned stratum of 12 records is 6 draws. Under 10
   draws it is unresolvable and read from the deep tier, as a stratum under 10 records was. The
   excess line prints that count. Where no records were read it prints the row's records, and says
   so.

## Rationale

`cargo xtask envelope`, the candidate rescored at one rev, both tiers:

| | `0088` | this |
|---|---|---|
| committed tier: domination failures | 1 | **0** |
| committed tier: unresolved excesses | 33 | 23 |
| committed tier: unresolvable excesses | 726 | 691 |
| deep tier: domination failures | 100 | 97 |

- **Every committed-tier move is a `Log` stratum.** None crosses `α` the other way, so no loss
  appears. The repeats in `real_acos` and `real_sin_cos` move no verdict.
- **The deep tier keeps its real losses.** 64 draws resolve them as before. The 97 are the next
  fix queue; 66 are `Ad`'s `R` block near the identity, which `0076` decision 1 kept.
- **Rejected: twins as a separate stratum.** That moves corpus bytes, every baseline row and the
  sign coverage the twins exist for, to fix a statistics question.
- **Rejected: the envelope reading the corpus (this record's first form).**
  - **Staleness:** a corpus drawn again after the run regrouped the run's records with no error.
  - **Two parsers:** it parsed the corpus a second time, without `corpus::read`'s checks.
  - **Cost:** the envelope took 1.6 s where it had taken 0.16, before the deep tier's 382 MB.
- **Rejected: grouping by equal errors.** It reads the answer the test is about. Two independent
  draws may tie by chance, and two twins of a program that does not canonicalize may differ.

## Consequences

- **`0088` decisions 1 and 2 read "draw" for "record"** when this record is `ready`.
  `PHASE1.md` §8's "per-record" pairing gains the same sentence.
- **The records file gains a `draw` cell.** Every subject's file written from now on carries it.
  Only the candidate's draws are read, so an oracle's older file still pairs.
- **The excess line counts draws**, so a twinned stratum's report shows its real `n`.
- **`MANIFEST-deep.json` is regenerated** in the same PR. It was stale (`0088`'s owed item), and
  the deep evidence above was drawn against it.

## Implementation plan

1. `corpus::draws`, the records file's `draw` cell, `envelope::units` and the grouping in
   `bars::swap_of`:
   - tests `a_quaternion_and_its_negative_are_one_draw`, `quaternion_inputs_are_four_vectors`,
     `the_draw_column_is_read_and_0088_files_are_one_draw_per_record` and
     `a_twin_pair_is_one_draw`;
   - the regenerated `MANIFEST-deep.json`.

   Verified by `cargo nextest run -p xtask`, `just lint`, `just corpus-deep-check`, and
   `cargo xtask envelope` with and without `--deep`.

## Further work

- **The generator's repeats.** `real_acos` writes 27 of 64 `x:near±1` inputs twice or more.
  Drawing distinct inputs would give those strata their 64 draws back.

## Open questions

1. **Whether a draw is the right unit for the no-regress bar too.** It is exact per stratum and
   needs no test, so this record leaves it alone.
