# 0019: A Cholesky solve without the transpose

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #34.

## Context

`solve_upper`'s documentation says that for a Cholesky factor `L` you pass `l.transpose()`, so the
way to solve `A x = b` from `chol`'s factor is

```rust
solve_upper(&l.transpose(), solve_lower(&l, b))
```

Every consumer that solves against a factor has to know to write that line, and it builds an
`N x N` copy whose only reader is a solve that could read `L`'s column `i` instead: in
`solve_upper`, `u_ik = l_ki`, and for `k > i` that is column `i` of `L` below its diagonal, in
storage order (`Matrix` is column-major).

**Measurement** (advisory, not a gate: a scratch crate outside the repository, not committed, since
no document prescribes a benchmark for `helicoid-linalg`; `chol_solve` as in the Decision against the
composition above). AMD EPYC-Milan on a shared host (load average 4 to 10 on 8 cores), `rustc`
1.98.1, `f64`, pinned to one core with `taskset`. Data: 64 `(L, b)` pairs (`L = chol(B B^T + I)`,
`B` entries in `[-0.5, 0.5)`); each pass calls the solve on every pair through `black_box` and sums
one output lane; 20,000 passes per repetition, the order of the two variants alternated between
repetitions, the two outputs' bits checked equal first. Two call shapes: inlined into the loop, and
behind an `#[inline(never)]` wrapper. The scratch crate is its own workspace with a path dependency
on `helicoid-linalg`; `[profile.release]` is absent for the first row, which is cargo's default
(`lto` unset, thin-local LTO, 16 codegen units), and the other rows set the named key. The ratio is
`chol_solve` over the composition, on the minimum over repetitions, range over the runs; absolute
times drift between runs on this host (2x), so only ratios within a run are compared.

| build | `N` | shape | ratio |
|---|---|---|---|
| `opt-level = 3`, cargo's default | 3 | inlined | 0.98 to 1.00 |
| | 3 | not inlined | 0.97 to 1.00 |
| | 6 | inlined | 0.96 to 1.01 |
| | 6 | not inlined | 0.98 to 1.00 |
| `opt-level = 3`, `lto = "fat"`, 1 codegen unit | 3 | inlined | 0.99 to 1.00 |
| | 3 | not inlined | 0.99 to 1.01 |
| | 6 | inlined | 0.99 to 1.03 |
| | 6 | not inlined | 0.98 to 1.00 |
| `opt-level = 3`, `lto = "off"` | 3 | inlined | 0.90 to 0.91 |
| | 3 | not inlined | 0.72 to 0.89 |
| | 6 | inlined | 0.71 to 0.76 |
| | 6 | not inlined | 0.83 to 0.84 |
| `opt-level = 1` | 3 | inlined, not inlined | 1.00 to 1.01 |
| | 6 | inlined, not inlined | 0.68 to 0.70 |
| `opt-level = 0` | 3 | inlined | 0.68 to 0.69 |
| | 6 | inlined | 0.64 to 0.66 |

Every row is three runs of 51 repetitions, except `opt-level = 0` (two runs of 15). **Whether the
copy costs anything depends on the optimizer removing it.** In cargo's default release profile and
with fat LTO it does, and the two are equal to within the run-to-run spread of this host (up to 5%
at `N = 6`, in either direction). Without cross-crate inlining (`lto = "off"`), or without
optimization, it does not: `chol_solve` is 10% to 30% faster at `opt-level = 3`, and 30% to 35%
faster at `opt-level` 0, and at `N = 6` and `opt-level` 1. The `opt-level = 0` not-inlined `N = 3`
cell had one run at 1.00 (median ratio 0.81 in that run, 0.57 in the other), which is the noise of
this host and not a result.

## Decision

One public item in `helicoid-linalg`, beside `chol`, `solve_lower`, `solve_upper`:

```rust
pub fn chol_solve<S: Real, const N: usize>(l: &Matrix<S, N, N>, b: Vector<S, N>) -> Vector<S, N>;
```

It solves `A x = b` from the Cholesky factor `l` of `A`: `solve_lower(l, b)` (unchanged), then a
back substitution against `l^T` that reads `l` directly (`NUMERICS.md` §15.6):

```text
x_i = (y_i - sum_{k>i} l_ki x_k) / l_ii,   i = N-1 down to 0,  k ascending from i+1
```

- **Signature.** The conventions of `solve_lower`/`solve_upper`: the matrix by reference, `b` by
  value, the vector back, generic over `S: Real` and `const N`. It takes the factor, not `A`, so a
  factor is reused across right-hand sides and the `chol` mask stays the caller's (R4). No mask is
  returned: a solve has no comparison to report.
- **Contract.** The output equals, to the bit in every lane of every `S`,
  `solve_upper(&l.transpose(), solve_lower(&l, b))`, for every input (`NUMERICS.md` §15.6 has the
  proof: `transpose` computes nothing, and the two run the same operations on the same operands in
  the same order). The one exception is the sign and payload of a NaN produced by arithmetic, which
  Rust leaves unspecified (`PHASE2.md` §3, the same exception as `dual_value_is_plain_value`); a
  NaN is a NaN in both. `-0`, `+-inf`, subnormals and the extremes of the range agree to the bit.
- **Summation order** (D16), stated in the rustdoc and pinned by
  `the_transposed_solve_sums_left_to_right`: each sum in increasing `k` from its first term with
  `sum` (`+0` when empty), the result subtracted from `y_i` once, then divided by `l_ii`; `i`
  descends.
- **Domain** (R6, `NUMERICS.md` §12, §15.3): every `l_ii` nonzero and not NaN, `debug_assert!`ed;
  release never panics and a zero diagonal gives `+-inf` or NaN. The forward solve runs first and
  sees the same diagonal, so a debug build panics with `solve_lower: zero diagonal`, as the
  composition does. The strictly upper triangle of `l` is never read.
- **The column solve is not public.** `solve_lower_transposed(l, y)` (`pub(crate)`, the same
  loop with the same `debug_assert!` under its own message) exists so `chol_solve` is two named
  steps and so it is tested on its own against `solve_upper(&l.transpose(), y)`.
- **Reference twin** (D6, `NUMERICS.md` §14): the composition above stays as the definition of
  correct. The twin is written out of two public functions, so there is no `reference` item (§14's
  prose names the exception); the row of §14 spells it, and `chol_solve_matches_reference` (`f64`,
  `f32`, `Dual<f64, 2>`, `N` = 1..=6, lower triangle and `b` drawn from `+-0`, `+-inf`, `+-NaN`,
  `f64::MAX`, the smallest normal and the smallest subnormal as well as plain values; the diagonal
  in the domain; junk above the diagonal) and `chol_solve_matches_reference_on_plain_values`
  compare with every lane's `to_bits`, all NaNs as one value.
- `solve_lower` and `solve_upper` are unchanged, except that `solve_upper`'s rustdoc points to
  `chol_solve` instead of telling the reader to pass `l.transpose()`.

## Rationale

**Why ship it, and what it does not buy.** Not speed in a release build: in cargo's default profile
and with fat LTO the two are equal to within the noise (Context), so a consumer that builds either
way loses nothing by writing the composition. What it buys is one call whose equality with the
documented recipe is tested, so a consumer never re-derives it or writes the forward solve against
the wrong triangle; and a result that does not depend on the optimizer removing an `N x N` copy,
which builds with `lto = "off"`, at `opt-level` 0 or 1, or behind a call the inliner declines do not
get (10% to 35% slower, about 30% at `N = 6`). The cost is about 30 lines and one function whose
result is defined to be the old one, so it can be adopted or ignored without changing any number.
The alternative "do nothing and document the composition" is coherent and is what a default-profile
or fat-LTO consumer measures as equal; this record rejects it for the tested recipe and the
optimizer independence, not for a release-profile speed-up, and the owner may still prefer it.

**No consumer is named.** No document in `docs/` cites an `A^{-1} b` from a factor:
`Gaussian::mahalanobis_sq` (`PHASE5.md`) uses `solve_lower` alone. `chol_solve` is a convenience
for the recipe `solve_upper` documented, added on the owner's request, with no caller in the
workspace to migrate. It is a small, speculative surface; the column solve stays private because it
would be a second one, with no use even in the recipe.

**Why one public item and not two.** A public `solve_lower_transposed` (solve `L^T x = y`) would
make the trio `solve_lower`, `solve_upper`, `solve_lower_transposed`, and then `solve_upper` with
a transposed argument for symmetry: a family with no consumer. `A^{-1} b` is the recipe `solve_upper`
already documented and the one solve a factor is made for; `L^{-T} z` alone (sampling from an
information matrix) is a different use, with no named caller. When one appears it is a record that
promotes the `pub(crate)` function, with the contract already written and tested here.

**Alternatives rejected.**

- *`solve_upper` with a `transposed: bool`.* A boolean mode flag, and a second summation order
  behind one name; the domain and `debug_assert!` messages stop being one function's.
- *`chol_solve(a, b)` that factors first.* It re-factors per right-hand side and has to return a
  mask; the caller already holds `(L, mask)` from `chol` (R4).
- *A multiple-right-hand-side form.* No consumer; `b` is a `Vector`, a loop is the caller's.
- *Reordering the sum (row by row, or a column-oriented axpy update) for the column reads.* Not
  bit-identical to the composition, so the twin could not be the definition; it would need its own
  error bound and corpus stratum ([`0006`](./0006-the-instrument-comes-first.md)) for no measured
  gain. The column read alone gives the contiguous access and keeps the sums.

## Consequences

- One new public item in `helicoid-linalg` (`API.md` §2, `PHASE2.md` §4, `CHANGELOG.md`);
  `0.0.x`, nothing breaks.
- A new invariant to maintain: `chol_solve` equals its twin to the bit (NaN sign and payload aside)
  for every `S`. Any future change to `solve_lower` or `solve_upper` (a different sum order, say)
  changes `chol_solve` in step or fails `chol_solve_matches_reference`.
- No corpus stratum, no bar change: `chol_solve` is a composite of `solve_lower` and `solve_upper`,
  which are covered by the Higham-bound proptests and have no corpus stratum yet (`PHASE2.md` §8),
  so `chol_solve` meets [`0006`](./0006-the-instrument-comes-first.md) no better and no worse than
  its parts and adds no stratum of its own. The error bound is that of the composition (Thm 10.4,
  §15.4), which `cholesky_solve_leaves_a_small_residual` now checks through `chol_solve`.
- The row in `NUMERICS.md` §14 has no `reference` item (`helicoid-linalg` has no such module); its
  prose names the exception, and the twin-table lint owed by `PHASE1.md` §0.0 has to accept a row
  whose twin is a composition of public functions.
- Budget (0007) unchanged: no dependency, no `alloc`, no `unsafe`.

## Implementation plan

1. `chol_solve`, `solve_lower_transposed`, the tests above, the `NUMERICS.md` §12, §14, §15.3,
   §15.6 edits, the `PHASE2.md` §0.0 row and §4 text, the `API.md` row, `CHANGELOG.md`, this
   record's index rows; one PR — verified by `just lint`, `just test`, `just doc`, `just no-std`,
   `just msrv`, and `PROPTEST_CASES=1000000 PROPTEST_RNG_SEED=1 cargo nextest run --release -p
   helicoid-linalg chol_solve` (both twin proptests, 10^6 cases each).

## Open questions

None.
