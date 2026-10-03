# 0026: `write_dense` pays for zeros, not for checks

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

`0025` *Consequences* left one item open: `Jac::write_dense` reaches a [`StridedMut`] only through
`set`, which asserts `r < rows && c < cols` and recomputes a saturating index on every store, even
though the impl's `debug_assert!` validated the shape once. At `DOF = 9` the emitted release asm is
**546 instructions, 81 branches, 48 compares and 81 panic calls for 81 stores**, on the path a
solver crosses per residual block. The obvious reading is that the checks are the cost and a
checked-once bulk writer on `StridedMut` is the fix.

**That reading is wrong, and the measurements say so.** Static instruction count is a bad proxy
here: the 81 panic calls sit on never-taken branches that the predictor handles for free and the
layout places cold, and the loop is unrolled at a `const` `N`, so `r == c` folds away and LLVM
already common-subexpression-eliminates a third of the index arithmetic (48 compares, not 162).

Four implementations of "write `k I` into a view", all in one binary so that inlining is shared,
`RnJac<f64, N>`, minimum of 7 repetitions of 200 000 iterations each, release, `black_box` on both
the buffer and the scalar, AMD EPYC-Milan. `set` is today's per-entry loop; `write_col` materializes
a `[S; N]` column and hands it over; `fill` passes a closure and lets the view walk itself with
induction variables and one bounds check; `zero+struct` fills the whole view in bulk and then writes
only the structural non-zeros. *Ceiling* is the same result written straight into a `&mut [f64]`
with `slice::fill` and `N` diagonal stores, with no view at all — the best any API could do.

| destination | DOF | `set` | `write_col` | `fill` | `zero+struct` | ceiling |
|---|---|---|---|---|---|---|
| column-major, packed | 3 | **2.62 ns** | +326% | +221% | 5.83 ns (+123%) | 1.62 ns |
| column-major, packed | 6 | 13.11 ns | +155% | +81% | **8.38 ns (−36%)** | 3.72 ns |
| column-major, packed | 9 | 29.68 ns | +53% | +69% | **10.94 ns (−63%)** | 6.20 ns |
| column-major, `cs = 16` | 9 | 29.67 ns | +63% | +69% | **22.27 ns (−25%)** | — |
| row-major (Ceres) | 9 | 29.69 ns | +84% | +62% | 38.15 ns (+28%) | — |
| block in a 64-row matrix | 9 | 29.73 ns | +63% | +60% | **25.25 ns (−15%)** | — |

Three things follow.

**The checks are not the cost.** Both check-removing designs lose everywhere, by 53% to 326%.
`write_col` loses because materializing a column doubles the write traffic: each cell is written
once to the stack and once to the destination, and for `DOF = 3` the `memcpy` setup dwarfs 9 stores.
`fill` loses because the closure defeats the unrolling — `rows` and `cols` are runtime values, so
the loop cannot unroll and each cell costs a call and an induction step instead of one store from a
straight line. Today's loop runs at 0.33 ns per entry, about one store per cycle, which is the
store-port limit for scalar stores: there is no headroom to win by removing arithmetic.

**The headroom is the structural zeros, and it is large.** The ceiling is 4.8x faster than today at
`DOF = 9`. A structured Jacobian is sparse by construction — that is what `0005` is for — so most
of what `write_dense` stores is zero: 72 of 81 entries for `k I`. Storing them one at a time is 72
stores; filling 648 contiguous bytes vectorizes to about 11. `zero+struct` captures 2.7x of the
available 4.8x through the existing `set` for the non-zeros.

**The win is conditional on size and on layout, and `write_dense` cannot see either.** It loses at
`DOF = 3` (+123%), where a bulk fill cannot amortize, and it loses on a row-major destination
(+28%), where a column is strided, nothing fills in bulk, and the second pass over the non-zeros is
paid for nothing. A `Jac` impl can read `rows()` and `cols()` but not the strides, so it cannot
make that choice; and a crossover in `DOF` is a switch point, which `0004` requires be measured and
generated rather than typed.

## Decision

**1. The per-entry `set` loop stays, and this record says why so that it is not "fixed" again.**
`Jac::write_dense` is not changed and no bulk writer is added to `StridedMut` by this record. The
rustdoc of `Jac::write_dense` gains one line: the per-entry checks are measured to be free, and the
cost is the structural zeros.

**2. `write_col` and `fill` are rejected on measurement**, in the shapes given in *Context*. Neither
is to be revisited without a destination layout that makes a column a contiguous source the impl
already holds — which a structured Jacobian, by definition, does not.

**3. The subject is instrumented before it is optimized.** A criterion bench for
`Jac::write_dense` lands under the `PHASE3.md` §11 benches, over the three destination layouts
(packed column-major, padded column-major, row-major) and every implemented `DOF`, gated by
`just bench-check`'s paired bootstrap against `baseline/`. `0006` is literal here: there is no
bench for this path today, so no variant of it may ship.

**4. The two-phase form is specified but deferred to its real subject.** `RnJac` is the *most*
sparse Jacobian the workspace will ever have — `N` non-zeros in `N²` — so the table above is its
best case and not a prediction for the Jacobians a solver actually calls. The decision whether to
adopt it belongs to the PR that lands `SEn3Jac` (`PHASE3.md` §5, blocks, `DOF` 6 and 9) and
`ProductJac` (§7, block-diagonal), measured on their own sparsity with the bench of decision 3. If
it is adopted, the shape is:

```rust
fn write_dense(&self, out: &mut StridedMut<'_, S>) {
    debug_assert!(out.rows() == N && out.cols() == N, "...");
    out.fill_value(S::zero());          // one bounds check; bulk where the layout allows
    for i in 0..N {
        out.set(i, i, self.k);          // the structural non-zeros only
    }
}
```

with `StridedMut::fill_value(&mut self, v: S)` checking `fits` once and then filling the whole view
as one span when `rs == 1 && cs == rows`, each column as a span when `rs == 1`, and strided
otherwise. Values are unchanged and no arithmetic is reordered, so D16 bit-identity is untouched and
the existing bounds hold with no re-measurement.

## Rationale

The alternative to decision 1 was to add the bulk writer now on the strength of the asm. It was
measured instead, and it lost: a 4.6x headroom exists but points at the zeros, not the checks, so
the API addition that the instruction count seemed to justify would have made every call slower.

Decision 3 over "adopt the two-phase form now for `RnJac`": `Rn` is not a hot group — a translation
Jacobian is `±I` and a consumer that needs one usually needs no dense form at all — so adopting it
there would buy a benchmark result rather than a solver improvement, and would set the crossover
from the least representative Jacobian in the workspace. `0021` is the precedent for recording a
measured non-change: four `chol` variants were measured and none cleared its stated bar.

The ceiling column is kept in the table because it bounds what any future design can claim. An
`unsafe` unchecked store is not among the alternatives (`#![forbid(unsafe_code)]`), and it would not
help: the ceiling is reached with safe `slice::fill`.

## Consequences

- One rustdoc line on `Jac::write_dense`, and this record, stand between the next reader and a
  plausible optimization that is a regression. That is the deliverable.
- A bench and a `baseline/` entry for `write_dense` exist from decision 3 on, so the §11 bench gate
  covers this path and a future regression is caught rather than argued.
- `StridedMut` gains no public item now; `PHASE2.md` §5's surface is unchanged, and the `fill_value`
  sketch in decision 4 is a *proposal*, not a contract. A draft authorises nothing.
- If decision 4 is later adopted, the crossover in `DOF` and the layout predicate are measured
  inputs, not typed constants, per `0004`.

## Implementation plan

1. The rustdoc line on `Jac::write_dense` and this record — verified by `just doc` and
   `cargo xtask lint`.
2. The criterion bench of decision 3 and its `baseline/` entry, over the three layouts and every
   implemented `DOF` — verified by `just bench-check` passing on an unchanged tree, and by the
   numbers of *Context* reproducing within the paired-bootstrap interval.
3. Deferred to the `SEn3Jac`/`ProductJac` PRs: measure the two-phase form on their sparsity and
   either adopt it behind a measured layout predicate or record that it lost there too — verified by
   `just bench-check` and the existing `jac_dense_order` law, which pins the written values and so
   makes the change observable only in time.

## Open questions

- **How does a `Jac` impl learn that a bulk fill will pay?** It can read `rows()` and `cols()` but
  not the strides, and D2 says layout is not a contract, so exposing `rs`/`cs` is not on. The
  candidates are a predicate (`fn row_contiguous(&self) -> bool`), moving the whole two-phase write
  inside `StridedMut` behind a structure descriptor, or having `fill_value` be the only primitive
  and accepting the row-major regression. Not answerable before decision 3's bench exists.
- **Where is the `DOF` crossover, and is it one number or one per layout?** Measured at 3 (loses)
  and 6 (wins) for the sparsest Jacobian only. `0004` governs how it is recorded if it is needed.
- **Is `RnJac` worth exempting even if the two-phase form is adopted?** It is the case the
  crossover hurts most at `DOF = 3`, and `Rn<S, 3>` is a real instantiation.
