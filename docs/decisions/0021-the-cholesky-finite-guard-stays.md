# 0021: The Cholesky finite guard stays

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #36, records the decision only; no code change.

## Context

`chol` (`NUMERICS.md` §15) stores an entry `x = (a_ij - s) / L_jj` only if `ok_j` and `fin(x) = (x *
0 <= 0)`, else `+0`, and folds `fin` into the mask. Per entry, `N(N-1)/2` times, the guard costs a
multiply and a compare (`fin`), four mask operations (`or`, `not` and `and` for the fold, `and` for
the store) and a select. The question is whether a cheaper guard gives a bit-identical `(L, mask)`.
Bit identity is the bar because §15.2 and `PHASE2.md` §0.0 pin the contract: mask `0 < d_j` and
every entry finite, `L_jj = 1` and a zero column after a failed pivot, an overflowing entry stored
as `+0`, the summation order to the bit.

**What each guard protects** (`ok_j` set means `d_j` finite and positive, so `L_jj = sqrt(d_j)` is
finite and positive):

- The entry guard is needed. `x` is non-finite for three reasons, not one: `a_ij` is NaN or `inf`
  (read from `a`); `s` overflows (the stored `L` entries are finite, their products are not) or is
  NaN (`inf - inf`); or the numerator is finite and the quotient overflows (`L_jj < 1`). Since
  `L_jj` is finite and positive, a non-finite numerator gives a non-finite `x`, so testing `x` alone
  covers all three. Dropping it stores `inf` or NaN, which feeds later columns; flushing at the end
  instead changes every later entry that would have read the `+0`.
- When `ok_j` is clear the entry guard is redundant in value and mask (the entry is `+0`, the mask
  term `fin or not ok` is true). Lanes evaluate every arm, so it is computed anyway.
- The pivot test needs `fin(d)` for exactly one input class: `d = +inf` (`a_jj = +inf`, `s` finite).
  `0 < d` already rejects NaN, `-inf` and every `d <= 0`. Without `fin`, `sqrt(inf)` is stored on
  the diagonal. It is per column, `N` per call, not the cost centre.
- A NaN in `a_jj` is caught by the mask (`0 < NaN` is false); a NaN in `a_ij` by the entry guard
  when `ok_j`, and is never read into a stored value when not. The upper triangle is not read.
- For `Dual` the guard tests the value lane; derivative lanes follow `sqrt` and the quotient (§12).
- **Why a sticky sum is a sound finiteness test** (the `fast` variant below). `acc + x` over the
  stored entries is non-finite iff some `x` is (`inf - inf` is NaN, NaN and `inf` are absorbing),
  unless finite terms overflow the sum. They cannot when every pivot is ok: `0 < d_j` forces `L_jk^2
  <= a_jj`, so every stored entry has `|L_jk| <= sqrt(a_jj) <= sqrt(MAX)` and at most `N(N-1)/2` of
  them are added. When a pivot fails the mask is false either way and `L` comes from the fall-back.

## Measurements

Scratch crate outside the repository (harness not committed; the variants are in the Appendix), path
dependency on `helicoid-linalg`, `--release`, fat LTO, `codegen-units = 1`, one AMD EPYC-Milan
machine (x86_64), pinned with `taskset` to a single core, two different cores tried. 2048
pre-generated matrices per batch, `black_box` on input and output, batches interleaved across
variants, 41 repetitions, **min** ns per factorisation as the statistic (on quiet cores the median
is within 1% of it). An earlier sweep on a loaded machine (load 12 on 8 cores) reported medians up
to 28% away from the min; it is discarded, its outliers were contention. Variants, all bit-identical
to `chol` (below):

- **hoist**: the mask fold once per column.
- **fast**: unguarded run, one sticky accumulator `acc + x`, and the current code as fall-back
  through `Real::branch`. The mask is the fast run's verdict (analysis above).
- **rows**: `L` accumulated row-major, so `L_ik` and `L_jk` are contiguous in `k`.
- **fast_rows**: both.
- **noguard**: no entry guard at all. Not a candidate (it breaks the contract); the guard's cost
  when its work is deleted, for reference only (see below).

Min change against the current `chol`, dense SPD input (mask set), range over three runs on each of
two cores:

| scalar, `N` | hoist | fast | rows | fast_rows | noguard (reference) |
|---|---|---|---|---|---|
| `f64`, 3 | -0% to +1% | -5% to -3% | 0% | -13% to -1% | -5% to -4% |
| `f64`, 6 | -2% to 0% | **+5% to +18%** | -2% | -5% to -2% | -6% to -3% |
| `f32`, 3 | 0% | -6% to -1% | 0% to +1% | -6% to -1% | -7% to -5% |
| `f32`, 6 | -4% to -2% | **+10% to +18%** | -5% to -2% | -7% to -4% | -7% to -5% |
| `Dual<f64, 2>`, 6 | 0% | +2% to +7% | **+37% to +52%** | **+38% to +41%** | -1% to **+11%** |

Failing input (the fall-back runs), `f64`, `N = 6`, one core: `fast` +40% at 37% failures, +67%
at 77%, +112% at 100%; `fast_rows` +34%, +60%, +104%. Across `f32` and `N = 3` the range is +10% to
+123%. `hoist` and `rows` stay within about 10% of `current` there.

Reading the table:

- The guard is off the critical path: the sums are fixed left to right and the division is serial,
  so the multiply, compare and select overlap the dependency chain.
- Code layout dominates below about 10%. The same binary gives `f64`, `N = 3`, `fast_rows` -1% in
  one run and -13% in another, an `f64`, `N = 3` conditioned-near-1 row shows `noguard` at +18%, and
  on `Dual` deleting the guard is 11% slower on one core and 1% faster on the other. `noguard` is
  therefore one observed codegen, not a bound on what a variant can win; `fast_rows` overlaps it at
  `f32`, `N = 6` and reaches it at `f64`, `N = 3`.
- Row-major storage loses on `Dual` (+37% to +52%); the cause was not investigated.
- `fast` loses on the success path at `N = 6` and loses more on failure; the cause was not
  investigated. This differs from an earlier, noisy sweep that showed it at about 0%.

Differential test, bitwise on `(L, mask)`, 10^6 cases per (variant, scalar, `N`) for `f64`, `f32`
and `Dual<f64, 2>` at `N = 3` and `6`, 24 combinations, **0 mismatches**, plus
`helicoid_linalg::chol` against the scratch copy of the reference. A reviewer's independent re-run
at 5 * 10^4 cases per combination also found 0. Ten generator modes: pure specials (`+-0`, `5e-324`,
`MIN_POSITIVE`, `+-1e300`, `MAX`, `+-inf`, NaN), symmetric specials, SPD with one poked special, SPD
scaled over 10^-150 to 10^150, rank-deficient, indefinite, non-symmetric dense, near-singular
pivots, subnormal, and a tiny pivot under near-overflow off-diagonals. 24% to 26% of `f64` cases and
15% of `f32` cases set the mask. A variant with the fall-back removed (a seeded defect, `fast`
without `branch`) fails on more than half of the cases, so the corpus discriminates. The generator
is not committed; the modes above are its specification.

## Decision

Keep the per-entry guard and the current `chol` as it is. No code change, no new proptest, no
change to `NUMERICS.md` §15 or the `PHASE2.md` status row.

Three standing rules follow from the evidence, so the question is not reopened by a benchmark alone:

1. **The bar for rewriting this kernel.** A gain of at least 15% at `N = 6` on a pinned core, on the
   success path, with no regression on failing input or on `Dual`, and `(L, mask)` bit-identical on
   every input (NaN sign and payload aside, D16 and `0018`). The best bit-identical variant measured
   at most 7% at `N = 6` (`f32`, `fast_rows`) and lost 38% to 41% on `Dual`. The number is this
   record's choice for this kernel, in the spirit of `tf_tree` 0016 (a ~12% win did not pay for nine
   crates), as `0013` cites it; it is not a repo-wide rule.
2. **`L` stays finite for every input.** It is a public guarantee (the `chol` rustdoc, `NUMERICS.md`
   §15.2 and §12, the `PHASE2.md` row): a caller that solves or samples from `L` never has to test
   it for `inf` or NaN. A guard that only cleared the mask would push that check onto every caller,
   so this record proposes no spec change to §15.2 or §12.
3. **Re-measure on a trigger, not on a hunch.** When `just bench` exists (`PHASE1.md`) and can time
   aarch64, wasm32 or `thumbv7em`, or when a lane `Real` exists (`0013`), the variants in the
   Appendix are the starting point.

## Rationale

`NUMERICS.md`, `0006` and `0013` define no gain a scalar kernel rewrite must reach, and `chol` has
no bench in `just bench` or the baseline, so no repo gate is at stake. Of the four bit-identical
variants tried, `fast` and `fast_rows` add a data-dependent second pass that costs +10% to +123% on
failing input, `fast` also loses on the success path at `N = 6`, `rows` loses 37% to 52% on `Dual`,
and `hoist` is within noise. The best success-path result at `N = 6` is `fast_rows`, `f64` -2% to
-5% (`noguard` -3% to -6%) and `f32` -4% to -7% (`noguard` -5% to -7%). That is a few percent on one
machine for a rewrite of a function whose summation order is pinned to the bit: complexity the
contract does not need. Decision rule 1 states the bar this decision is measured against. Under
`0013` a lane mask evaluates both arms, so `fast` is expected to lose there; that is reasoned, not
measured.

## Consequences

- `chol` stays the single implementation; no reference twin is added (`NUMERICS.md` §14 has no
  `chol` row).
- A future SIMD `Real` (`0013`) or an `f32` embedded target owes a re-measurement (`0013` points
  here); this record is evidence for one x86_64 machine only.

## Implementation plan

None.

## Open questions

None. The unmeasured targets and the run-to-run spread (about 10% at `N = 3` from code layout
alone, two cores of one machine) do not change a "no change" decision; they are the re-measurement
trigger above.

## Appendix: variants

`sum` is the fold used by the crate (first term, then `acc + t`); `is_finite(x) = (x * 0 <= 0)`. The
reference `chol_ref` is the crate's `chol`. `fast_rows` is `fast` with the row-major storage of
`rows`. The bodies are given as deltas against `chol_ref`.

```rust
// chol_ref: per column j
//   d   = a[j][j] - sum(L[k][j]^2, k < j)
//   ok  = 0 < d && is_finite(d);   pd &= ok
//   ljj = select(ok, d, 1).sqrt(); L[j][j] = ljj
//   per i > j:  s = sum(L[k][i] * L[k][j], k < j);  x = (a[i][j] - s) / ljj
//               fin = is_finite(x);  pd &= fin || !ok;  L[j][i] = select(ok && fin, x, 0)

// hoist: replace the per-entry `pd &= fin || !ok` by, per column,
//   colfin &= fin  (per entry);   pd = pd & ok & (colfin || !ok)  (after the i loop)

// fast: per entry  acc = acc + x;  L[j][i] = x            (no fin, no select)
//   verdict = pd && is_finite(acc)
//   L = branch(verdict, || L, || chol_ref(a).0);  mask = verdict

// rows: store L row-major (rows[i][j]); the sums read rows[i][..j] and rows[j][..j]

// noguard: per entry  L[j][i] = x  and no fold into pd (contract-breaking; reference only)

// seeded defect: fast without the branch fall-back
```

