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
- **How many.** 53 of `sen3_log_n2`'s 54 strata are twin pairs, and so are most strata of every id
  with a quaternion input.
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

1. **The swap test's unit is a draw, not a record.** Two records of a stratum are one unit when
   every input is equal, a quaternion input (`q`, `q0`, `q1`: four components) up to its sign.
   - **Reading.** The units are read from the corpus inputs alone, never from an error, so the
     grouping cannot lean on the answer it tests.
   - **Pairing.** Each unit is one pair, each side's maximum over the unit's records. The units
     swap whole, so the maxima are what a swap moves.
2. **Resolvability counts units.** `0088` decision 2's floor is computed over the units, so a
   twinned stratum of 12 records is 6 draws. A stratum under 10 draws is unresolvable and is read
   from the deep tier, as one under 10 records was.
3. **Each tier reads its own corpus.** The committed tier groups by `conformance/corpus`, the deep
   tier by `conformance/corpus-deep`. Where a tier has no manifest, each record is its own unit,
   which is `0088` exactly.

## Rationale

`cargo xtask envelope`, the candidate rescored at one rev, both tiers:

| | `0088` | this |
|---|---|---|
| committed tier: domination failures | 1 | **0** |
| committed tier: unresolved excesses | 33 | 23 |
| committed tier: unresolvable excesses | 726 | 691 |
| deep tier: domination failures | 100 | 97 |

- **Every committed-tier move is a twinned stratum.** None crosses `α` the other way, so no loss
  appears.
- **The deep tier keeps its real losses.** 64 draws resolve them as before. The 97 are the next
  fix queue; 66 are `Ad`'s `R` block near the identity, which `0076` decision 1 kept.
- **Rejected: twins as a separate stratum.** That moves corpus bytes, every baseline row and the
  sign coverage the twins exist for, to fix a statistics question.
- **Rejected: grouping by equal errors.** It reads the answer the test is about. Two independent
  draws may tie by chance, and two twins of a program that does not canonicalize may differ.

## Consequences

- **`0088` decisions 1 and 2 read "draw" for "record"** when this record is `ready`.
  `PHASE1.md` §8's "per-record" pairing gains the same sentence.
- **The envelope reads the corpus inputs once per run**: 63 MB committed, 382 MB deep, in seconds.
- **The excess line says "independent draws"**, so a twinned stratum's report shows its real `n`.
- **`MANIFEST-deep.json` is regenerated** in the same PR. It was stale (`0088`'s owed item), and
  the deep evidence above was drawn against it.

## Implementation plan

1. `envelope::units` and the grouping in `bars::swap_of`:
   - tests `a_quaternion_and_its_negative_are_one_unit_and_nothing_else_is`,
     `only_q_named_four_vectors_are_quaternions` and `a_twin_pair_is_one_draw`;
   - the regenerated `MANIFEST-deep.json`.

   Verified by `cargo nextest run -p xtask`, `just lint`, `just corpus-deep-check`, and
   `cargo xtask envelope` with and without `--deep`.

## Open questions

1. **Whether a draw is the right unit for the no-regress bar too.** It is exact per stratum and
   needs no test, so this record leaves it alone.
