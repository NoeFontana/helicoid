# 0078: A comparison drops ties and prints both maxima

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #147

## Context

[`0038`](./0038-a-program-comparison-is-not-a-bar.md) decisions 4 to 6 fix the shape of
`cargo xtask compare`: two in-process subjects, the corpus read once, both evaluated on every
record, and per `(fn, stratum, precision)` a paired direction, a sign test, the geometric mean of
the per-record error ratio and "the per-stratum maximum ratio in both directions", pooled per id
and then combined. Implementing it met six readings that record does not fix: what a tie is, what a
zero error or a non-finite output does to each statistic, which of two quantities "maximum ratio"
names, what is pooled and how, and what a run's scope and exit status are. The owner's instruction
was to decide them conservatively and record them here.

[`0076`](./0076-the-oracles-first-verdict-r-of-q-stays-and-drift-is-excepted.md) decision 1 is
the first comparison the task can re-measure: a hand tally of `R(q)`'s unit form against the
homogeneous one, made by editing `Quat::to_matrix` and rerunning `conformance`. Its *Further work*
asks for exactly this.

## Measurement

### `helicoid` against `helicoid:unit-r`, both precisions

`just compare --a helicoid --b helicoid:unit-r [--precision f32]`, release build, 0.07 s each.
The pooled rows and the per-id tally of strata by `max_b / max_a`:

| id | prec | records `b<a` / `b>a` / tie | sign `p` | geomean `b/a` | strata better / worse / equal | `max_b/max_a` |
|---|---|---|---|---|---|---|
| `sen3_ad_n1` | f64 | 111 / 96 / 117 | 0.331 | 0.936 | **22 / 20** / 12 | 0.303–1.951 |
| `sen3_ad_n2` | f64 | 107 / 100 / 117 | 0.677 | 0.933 | **22 / 20** / 12 | 0.284–2.116 |
| `sen3_ad_n3` | f64 | 107 / 100 / 117 | 0.677 | 0.942 | **22 / 20** / 12 | 0.361–1.884 |
| `sim3_ad` | f64 | 163 / 141 / 212 | 0.228 | 0.913 | **36 / 51** / 35 | 0.276–2.020 |
| `sim3_ad` | f32 | 151 / 103 / 238 | 3.11e-3 | 0.937 | **43 / 32** / 41 | 0.452–1.940 |

**`0076`'s tally is reproduced exactly**: the strata counts of every id at both precisions, and the
ratio column to its two decimals (0.28–2.12 over the three widths, 0.28–2.02, 0.45–1.94). So are the
two replays that found the cause: `sim3_ad/sigma:-3/theta:0.3` goes from 1.350 to 0.373 `u`
(Sophus 0.426) and `sen3_ad_n2/theta:1e-2` from 1.565 to 0.445 `u` (manif 0.445).

The pooled rows are not a verdict, and should not be: each id pools two ends of `θ` that go opposite
ways. Summing the task's per-stratum rows by the angle the stratum names (`θ ≤ 1`, with `exact0` and
`subnormal`; `θ ≥ 3`, with `q:w0`) gives the direction `0076` states, per record:

| id | prec | `θ ≤ 1`: records `b<a` / `b>a`, `p` | `θ ≥ 3`: records `b<a` / `b>a`, `p` |
|---|---|---|---|
| `sen3_ad_n1` | f64 | 91 / 19, 1.9e-12 | 19 / 72, 2.0e-8 |
| `sen3_ad_n2` | f64 | 87 / 23, 6.0e-10 | 19 / 72, 2.0e-8 |
| `sen3_ad_n3` | f64 | 88 / 22, 1.6e-10 | 18 / 73, 5.0e-9 |
| `sim3_ad` | f64 | 122 / 49, 2.3e-8 | 23 / 84, 2.4e-9 |
| `sim3_ad` | f32 | 108 / 22, 8.0e-15 | 26 / 75, 1.1e-6 |

It does **not** hold stratum by stratum, as `0076`'s prose has it ("wins every stratum at `θ ≤ 1`
and loses every one at `θ ≥ 3`"). By the maximum, at `θ ≤ 1` each SE_N(3) width has 22 better, 11
equal and 2 worse (`rho:1e3/theta=1e-4` and `rho:1e4/theta=1`, 1.009–1.177×), and `sim3_ad` has 30
better, 31 equal and **20 worse** at `f64` (up to 1.492×), 36 / 37 / 4 at `f32`. At `θ ≥ 3` each
SE_N(3) width has 17 worse and one equal (`theta:pi-1e-8`); `sim3_ad` 30 worse, 2 better, 2 equal at
`f64` and 28 / 3 / 1 at `f32`. The counts in `0076`'s table were right and its sentence was
stronger than they are; on four to six records a maximum is one record (`0038`), and the
exceptions are within that. `0076` decision 1 does not rest on "every" and is unchanged.

### Which unit form

The variable is Eigen's `toRotationMatrix`, whose diagonal is `1 − (2y² + 2z²)`. Spelled as `0076`
writes the identity, `(1 − 2‖u‖²) + 2x²` from the same products, it is a different program — a
subtraction of like-sized terms that Eigen's form avoids — and it does not reproduce the tally:
8 / 24 / 22 strata on `sen3_ad_n1`, 17 / 42 / 63 on `sim3_ad` at `f64`. The tally is
discriminating, and the twin is the program `0076` measured.

### The task against `0038`'s own number

`just compare --a seeded:correct --b seeded:w2-identity` reproduces `0038`'s table: over the nine
`sen3_*` ids, **21 strata better and 104 worse** by the maximum, each id net worse on its own, and
the two `sen3_exp_n*` exceptions at 0.73× and 0.64×. Per record, the nine pool to 172 better and
362 worse. `so3_jl` and `so3_jr` also move (9 better, 5 worse each; 159 / 330 records,
`p = 8.0e-15`), which `0038`'s table, over the `sen3_*` ids, does not list.

## Decision

1. **A record is a pair of conformance scores**, `metric::Rule::score` at the precision (the
   largest `NUMERICS.md` §11 error over the output fields). Read conservatively:
   - **a tie** is the two errors equal to the bit; it is counted, and it is out of the sign test and
     of every ratio;
   - **exactly one zero error** is a direction, so it is in the sign test; it is out of the
     geometric mean, whose logarithm is infinite there, and it reads `inf` in the extremes;
   - **a non-finite output** on either side, and **an unscored record** (a backward-only id), are
     counted per side and are out of every statistic. A non-finite output is a conformance failure,
     not a magnitude, and the comparison does not turn it into one. It does not fail the run.
2. **The statistics**, `b` read against `a`:
   - the **direction** is which subject is better on more records, and the **sign test** is the exact
     two-sided binomial over the records that differ, ties dropped;
   - the **geometric mean** of `b/a` is over the records that differ with both errors positive, so
     the strata the variable never reaches do not dilute it toward 1 (the number of such records is
     the row's `b<a` plus `b>a`, less its zero-sided ones);
   - "the per-stratum maximum ratio in both directions" has two readings, and **both are printed**:
     `worst_b/a` and `worst_a/b`, the largest per-record ratio each way, which is decision 6's
     "better on average and worse somewhere"; and `max_b/max_a`, the bar's two maxima as one ratio,
     with each id's count of strata below, above and at 1 and its range, which is the form of every
     hand tally so far (`0036`, `0038`, `0076`), so that a published tally can be reproduced.
3. **`helicoid:unit-r` is a one-variable twin**: `Ad` of SE_N(3) and Sim(3) with `R(q)` in Eigen's
   unit form, assembled in `xtask` from the public API, because `SEn3Jac` and `Sim3Jac`'s blocks are
   crate-private ([`0028`](./0028-two-ready-specs-disagree-on-the-blocks-of-sen3jac.md)). Its control
   is that the assembly with `Quat::to_matrix` is `helicoid`'s answer to the bit on every record of
   the four ids (`the_assembly_with_the_shipped_r_is_the_shipped_ad`); its validation is
   reproducing `0076`'s tally (`the_unit_r_twin_reproduces_the_tally_of_0076`). It is planted, as
   `helicoid:host-std` is. That a twin is validated by a published measurement before it makes a
   new one is [`0037`](./0037-six-records-cannot-answer-a-one-variable-question.md) (draft) decision
   1; it is adopted here for this twin, and decides nothing for that record.
4. **Pooling and scope.** Per id is every record of the id's strata at the run's precision, one
   sample; the combination is every id's records, printed after the ids and labelled as not
   independent. No Fisher or Stouffer combination: each assumes the independence `0038` decision 5
   says the ids do not have. A run is one precision (`--precision`, `f64` by default; at `f32` the
   `@f32` strata alone, as `conformance` scores them), over the ids both subjects answer; an id one
   of them answers is named on stderr. Either subject may be any registered in-process subject,
   planted included, and the two may be the same, which is a control: every record ties. The task
   fails only when it cannot run (an unknown subject, nothing both answer, an `f32` run of a subject
   without an `f32` kernel, a metric error); a finding never fails it.
5. **The task is validated as a twin is**: by reproducing `0038`'s published `W²` tally
   (`the_task_reproduces_the_tally_of_0038`), and its statistics by unit tests that read no corpus
   (`compare::stats`).

## Rationale

Every reading above is the one that hides least. Dropping ties is the textbook sign test and the
only reading under which a stratum the variable does not touch says nothing; counting them would
let 27 bit-identical strata (`0038`'s `sen3_jl_n1` at 64 records) argue for "no effect". A zero
error is the strongest evidence of direction there is and the weakest of magnitude, so it is in
the one and out of the other. Printing both readings of the maximum costs two columns; picking one
would have made either `0076`'s tally irreproducible or decision 6's "worse somewhere" invisible.

The `θ` split above is the case for keeping the pooled row beside the ids and never instead: the
pooled `sen3_ad_n1` row says `p = 0.33`, which is two strong and opposite effects cancelling, not
"no effect". A stratum filter would let the task print the split itself; it is *Further work*,
because summing the printed rows answers it today and a filter is a new surface to argue about.

## Consequences

- `PHASE1.md` §5 gains one sentence and its schema is unchanged; the §0.0 row names the task.
- `0038` decision 4 is implemented; its plan steps 3 (`0036`'s application forms re-reported) and 4
  (`Q`'s word association as a twin) are now unblocked.
- A hand tally made by editing the library is now a twin and a test. A change that moves the
  `Ad` assembly, or the corpus of the four ids, fails `the_unit_r_twin_reproduces_the_tally_of_0076`
  and has to say why.
- `0076`'s "every stratum" sentence is corrected here, not there; its *Further work* bullet points
  at this record.

## Implementation plan

1. `cargo xtask compare`, `just compare`, the twin and the two controls, on one branch — verified by
   the tests named in decisions 3 and 5, `compare::stats`' unit tests, and `just lint`, `just test`,
   `just doc`.

## Open questions

None.

## Further work

1. A `--stratum` filter (a substring or a list), so a per-regime direction such as `θ ≤ 1` against
   `θ ≥ 3` is the task's own pooled row rather than a sum of its rows.
2. A per-record listing (`--records`), for the replay of one stratum that found `0076`'s cause.
