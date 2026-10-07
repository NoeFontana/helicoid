# 0053: `acos` is better, the roots are not necessarily

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record and `0022`'s two steps land together

## Context

[`0022`](./0022-real-owes-acos-and-cos.md) is `ready` and unimplemented. It adds `Real::acos` and
`Real::cos`, deletes `cubic.rs`'s private `atan2(sqrt((1 - x)(1 + x)), x)` spelling and takes the
sine `eig3` and `cubic` discard, on three stated reasons: **it is faster**, **it is more accurate**,
and **it restores bit-identity with omnisac**. Implementing it as written means re-measuring, and
the re-measurement does not support two of the three. `0022` is `ready`, so the index's rules allow
it to be amended only by a signed banner from a later record; this is that record.

Every figure below was taken on this host, with the protocol stated, against a reference this
repository can reproduce.

**`acos` the function: `0022` is right, independently.** Scored against `mpmath` at 60 digits, in
ulps of the exact value:

| region | `atan2` form | `libm::acos` |
|---|---|---|
| grid over `[-1, 1]`, 4095 points | 1.659 | **0.852** |
| `x -> 1-` (`1 - 2^-k`, `k = 1..52`) — the double root | 1.148 | **0.489** |
| `x -> -1+` | 0.520 | 0.520 |
| `x` near 0 (`+-2^-k`) | 0.855 | **0.520** |
| exact `+-1`, `0` | 0.276 | 0.276 |

Better or equal everywhere, by about 1.9x on the grid and 2.3x in the region `eig3` cares about.
`0022`'s own table says the same with a different ulp convention.

**Cost: `0022` is right.** 20 M calls over the region, three alternating passes: `libm::acos`
**6.3 ns** against the `atan2` form's **16.6 ns**, a factor of **2.63** against its stated 2.61;
`libm::cos` **7.7 ns** against `sin_cos().1`'s **9.6 ns**, so about 1.9 ns per phase and three
phases per solve.

**`cos` changes no value at all**, which `0022` assumes and does not establish. `libm::cos` and
`libm::sincos(..).1` are identical **by construction** — the same `rem_pio2`, the same kernels, the
same octant table and, unlike `sin`
([`0052`](./0052-real-owes-sin-and-the-corpus-does-not-move.md)), the same small-argument cut — and
measured: 0 disagreements over 800 010 binary64 arguments and over **all `2^32`** binary32 patterns.

**Reason 2 does not survive at the routine level.** A more accurate `acos` is not a more accurate
root: the root is `2 r cos(theta/3 - 2 pi k / 3)`, and the rounding downstream of `theta` is not
monotone in `theta`'s error. Measured, on the three `GOLDEN` rows that move, against `mpmath`'s
exact roots:

| row | slot | `atan2` form | `libm::acos` | |
|---|---|---|---|---|
| 0 | 0 | **0.522 u** | 7.478 u | the old spelling is closer |
| 1 | 1 | 10.339 u | **3.339 u** | the new one is |
| 2 | 1 | 4.068 u | **1.068 u** | the new one is |

And `eig3`'s worst `value / bound` over its fixture moves **0.150 -> 0.206** on `ev` and `res` —
`ang` does not move — which is a 1.37x degradation of the ratio, five times inside the bound.
`0022` says "more accurate, by a little, everywhere measured"; that is true of `acos` and false of
both routines that call it.

**Reason 3 is not verifiable in this repository.** `0022` says the trigonometric arm's values move
"back toward omnisac". omnisac is not a dependency, no committed artefact records its per-row
output, and `bits_are_omnisacs` **passed** before this change — so the port reproduced all 36 stored
triples, and after it three of them differ. Whether the new three are omnisac's cannot be told from
here.

**And nothing in the corpus can see any of this.** The committed corpus has 40 ids and **none** of
them is `cubic_*`, `eig3_*` or `real_*`: `PHASE2.md` §8 owes them and the Phase 1 generator has not
produced them. Verified rather than assumed — `just conformance` before and after moves **zero** of
1565 numeric rows. So `0006`'s domination and no-regress bars do not reach `solve_cubic` or `eig3`
at all, and the whole accuracy argument for this change rests on the two test modules' own bounds.

## Decision

1. **`0022`'s steps 1 and 2 are implemented as written.** `Real` gains `acos` and `cos`, the private
   spelling goes, `cubic` and `eig3` take `Real::cos` where they discarded a sine, and the moved
   rows are regenerated. The decision stands on reason 1 plus the function-level half of reason 2:
   **2.63x on `acos`, 20% on `cos`, and a strictly better `acos`**, with every recorded bound of
   both routines still holding by a factor of five.

2. **`0022`'s reason 2 is narrowed and its reason 3 withdrawn**, in the words of the Context above.
   The defensible claim is "the function is strictly better and the routines' recorded bounds still
   hold", not "the roots are more accurate" and not "bit-identity with omnisac is restored".

3. **The figures that moved are re-recorded where they live**, not only here: `cubic_tests`'s header
   loses the 5 `u` / 11 `u` deviation note and gains what is actually known about the three
   regenerated rows, and `eig3_tests`'s header moves 0.15 to 0.206 and says why.

4. **`bits_are_omnisacs` keeps its name and its doc says what it now pins.** Renaming it would lose
   the fact that the 36 rows *were* omnisac's when they were drawn, which is the only link to the
   reference implementation this repository still has. What it no longer is, is a live comparison.

5. **D5 becomes a lint instead of an absence.** D5 is "no `acos` of a trace anywhere in the
   workspace", and until now `Real` had no `acos`, so a group routine wanting one had to spell it
   out where a reviewer would see it. Adding `Real::acos` opens that hole, so `cargo xtask lint`
   gains check 8: **`crates/helicoid` does not name `acos`** in code. `helicoid-linalg` is exempt —
   `cubic` and `eig3` are the sanctioned callers, and their argument is a normalized determinant
   ratio, not a trace.

6. **`dual_value_is_plain_value`'s probe is corrected.** `PHASE2.md` §0.0 claims it runs "over every
   `Real` method"; the probe is a hand-written list of 26 expressions and `0052` added `Real::sin`
   without extending it, so the claim was already false. It becomes 29 with `sin`, `cos` and `acos`.

7. **`Dual::acos`'s derivative feeds `sqrt` a safe argument.** `-d / sqrt((1 - v)(1 + v))` is the
   rule `0022` states, and the product is negative outside `[-1, 1]` where `Real::sqrt` carries a
   `debug_assert!` — so a `select` substitutes `0`, giving `-+inf`, which is the convention
   `Real::atan2`'s rule already sets for a degenerate derivative. Not a clamp for the caller's
   benefit: it is the only spelling that neither panics in debug nor calls `libm::sqrt` directly,
   which lint check 5 forbids there.

## Rationale

The alternative is to implement `0022` and leave its reasons standing. That is the failure mode
`0022` itself diagnoses one level down — "a cost comparison against the function being replaced was
never run, and three documents reasoned about the trade without one" — and repeating it with an
accuracy comparison instead of a cost one would be worse, because the record would then assert
something three measurements contradict.

The alternative to implementing it at all is to keep the private spelling on the strength of the
routine-level figures. That loses 2.63x on a function `eig3` calls once and `cubic` once per solve,
to protect a 7 `u` movement on one polynomial of thirty-six and a bound ratio that stays at a fifth
of its bound. `0006` ranks accuracy first where a bar exists; here no bar exists, both routines'
own bounds hold, and the function being called is unambiguously the better one.

**Why the mixed routine-level result is not a reason to stop:** the three `GOLDEN` rows are a sample
of three, chosen by which rows moved rather than by what they represent, and the one that worsened
is at 7.5 `u` against a bound of `16 u ((W + Wd)/P + R)` that it does not approach. The measurement
that would settle the question is a corpus id for `solve_cubic` and `eig3`, which `PHASE2.md` §8
owes and *Further work* 1 names.

## Consequences

- `Real` is two methods wider; every implementor gains both, including the four test doubles and
  `xtask`'s two seeded lanes. An out-of-tree `impl Real` breaks, which `0.0.x` permits.
- The trigonometric arm of `solve_cubic` is one `atan2`, one `sqrt` and two multiplies lighter per
  call, and three `sin` kernels lighter per solve; `eig3` likewise.
- `cubic.rs` keeps `pi`: it is a compile-time constant and `0022`'s Context settles it separately.
- `bits_are_omnisacs` pins this port's bits on 36 rows that were omnisac's when drawn. The link is
  historical from here on, and decision 4 records that rather than hiding it behind a rename.
- D5 is enforced for the first time. Check 8 would have flagged nothing before this change, because
  there was nothing to flag.
- `PHASE2.md` §2's `libm` list, §3's `Dual` rules and §6's "`Real` has no `acos`" sentence move, as
  do `API.md` §2's row and `0017` decision 1's list of required methods.

## Implementation plan

1. This record; `Real::{acos, cos}` with their `Dual` rules and every impl; the two consumers; the
   three regenerated `GOLDEN` roots; both headers' figures; lint check 8; the probe; the doc edits —
   verified by `just lint`, `just test`, `just doc`, `just no-std`, `just msrv`, and
   `just envelope` holding at **113** (it cannot move: no corpus id reaches either routine).
2. `cargo xtask bench-gate --bench groups --against <baseline built in this tree>` — `eig3` and
   `solve_cubic` have no `benches/groups.rs` row, so this measures only that nothing else moved;
   the per-call figures of the Context are what stand for the win, and *Further work* 2 owes the
   rows.

## Open questions

None.

## Further work

1. **`solve_cubic` and `eig3` have no corpus id**, so `0006`'s bars do not reach them and this
   record's accuracy argument had to be made out of their own test modules. `PHASE2.md` §8 owes
   `real_*` ids for `sqrt`, `cbrt`, `sin_cos` and `atan2` through the Phase 1 generator; the two
   routines above want their own. Until then a change like this one is judged by bounds that are
   fits, not by a maximum over a stratum.
2. **Neither routine has a bench row.** `benches/groups.rs` covers the group crate only, so the
   2.63x and 20% here are per-call figures from a standalone harness rather than a gated
   measurement. A `linalg` bench target is owed and would have priced this change end to end.
3. **`cubic_tests`'s bound ratios are not reproducible from this repository.** Its header quotes
   "worst `err / bound` 0.88" from a scratch harness that is not committed, where `eig3_tests` has
   `EIG3_REPORT`. The same hook for `cubic_tests` would make the figure citable under `0040` item 4;
   this record did not add it because the ratios it quotes did not change.
