# 0057: `eig3` anchors the isolated end

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** the PR that lands this record

## Context

`eig3`'s eigenvector frame was anchored on the largest eigenvalue, as omnisac's is.
- `v2` is the longest cross product of the rows of `A − l2 I`;
- `v1` is that of `A − l1 I`, projected orthogonal to `v2`;
- `v0 = v1 × v2`.

Where the top pair is a tie, `v2` is rounding noise, and every column built from it is wrong. That
includes the isolated `l0`'s column, whose own gap is wide. [`0023`](./0023-eig3-departs-from-omnisac-and-its-limits.md)
(draft) records this as a limit and asks, in its open question 2, "anchor by mask, or build the
hybrid?".

[`0056`](./0056-the-routines-d7-does-not-reach.md)'s corpus id `eig3` measured the limit for the
first time. On the gap-weighted eigenvector scale, which reads `O(1)` at the conditioning, every
`eig:gap-1e-k/top` stratum with `k ≥ 8` read between `3.9e15` and `1.4e16 u`, and so did
`eig:rank1`, a double eigenvalue at 0. The same gaps at the bottom of the spectrum read `5e7`–`1e8`.
nalgebra's `SymmetricEigen` reads 5–13 on all of them.

## Decision

1. **The frame is anchored on the more isolated end**, chosen by a mask on the two computed gaps:
   `l2` where `l2 − l1 ≥ l1 − l0`, else `l0`.
   - The anchor's vector is the longest cross product of the rows of `A − l I`.
   - `v1` is that of `A − l1 I`, projected orthogonal to the anchor's.
   - The third closes the frame with `det = +1`: `v0 = v1 × v2` when anchored on `l2`, and
     `v2 = v0 × v1` when anchored on `l0`. It is one cross product, signed by the mask.
2. Where `l2 − l1 ≥ l1 − l0` the result is **bit-identical** to the frame anchored on `l2`.
3. This answers `0023` open question 2 and nothing else of it. Its questions 1 (a mask beside the
   result) and 4 (scaling) stay open. Its question 3 is answered by `0056`.

## Rationale

The cross product of two rows of `A − l I` is accurate to `u |A| / gap` in angle, `gap` being
`l`'s own. Anchoring on the end with the wider gap makes the anchor the best-conditioned vector of
the three. Once the anchor is right, the plane orthogonal to it is the invariant subspace of the
other two, so a tied pair's vector is right whatever it is.

Kopp's hybrid, the alternative, iterates where its accuracy estimate fails. It would also fix what
anchoring cannot, the eigenvalue error at a double root, but it is a data-dependent loop and a much
larger change. Anchoring costs a comparison and a few selects.

## Consequences

- The test that pinned the failure becomes
  `eig3_tests::a_tie_of_the_top_pair_anchors_on_the_isolated_end`. Residuals that were above
  `0.01 |A|` (and `0.4 |A|` on axis-aligned ties) are now below `64 sqrt(u) |A|`.
- What remains is the eigenvalue error of a close pair, up to `sqrt(u) |A|` from `acos` at
  `r = ±1`. It now dominates every gap stratum, at both ends alike: `5e7`–`1e8 u` at `f64` for gaps
  of `10^-8` and below. That is where `0023`'s hybrid, or an exact `2 × 2` deflation of the pair,
  is owed.
- `0056`'s committed rows for `eig3` are re-recorded.

## Implementation plan

1. Decision 1, the rustdoc and the test above, with the committed rows re-recorded — verified by
   `the_unverified_routines_reproduce_their_committed_rows` and `benches/linalg.rs`.

## Open questions

None.

## Measured

`max_u` per stratum, the shipped subject against the committed corpus. Only the strata that moved
are listed; every other stratum is equal or lower.

| stratum | before | after | nalgebra |
|---|---|---|---|
| `eig:gap-1e-8/top` | 3.9e15 | **1.1e8** | 6.7 |
| `eig:gap-1e-9/top` | 1.1e16 | **6.8e7** | 6.3 |
| `eig:gap-1e-10/top` | 1.2e16 | **7.0e7** | 7.5 |
| `eig:gap-1e-11/top` | 1.4e16 | **5.2e7** | 12.8 |
| `eig:gap-1e-12/top` | 7.6e15 | **6.3e7** | 12.8 |
| `eig:rank1` | 1.2e16 | **6.3e7** | 13.3 |
| `eig:scale-up` | 29.8 | **22.5** | 8.0 |
| `eig:random@f32` | 11.1 | **6.27** | — |
| `eig:gap-1e-4/top@f32` | 2.8e5 | **3.7e3** | — |
| `eig:gap-1e-5/top@f32` | 1.2e6 | **2.3e3** | — |
| `eig:gap-1e-6/top@f32` | 6.2e5 | **3.1e3** | — |
| `eig:rank1@f32` | 2.1e7 | **3.5e3** | — |

Latency, `benches/linalg.rs`, `bench-gate --against` a baseline built in this tree, core 5 pinned,
no other process above 20% CPU (concurrent A/A floor median 0.28%):
- `f64`: 1.017–1.018 on every row, about 2.7 ns of 164.
- `f32`: within its floor.

## Further work

1. The double root's eigenvalue: `0023`'s hybrid, or an exact `2 × 2` deflation once the anchor is
   accurate. Either is a formula and a record of its own.
