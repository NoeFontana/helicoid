# 0081: `eig3` is four Jacobi sweeps

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** the PR that lands this record

## Context

`eig3` was Smith's closed form: `q = trace / 3`, `p = sqrt(|A - q I|_F^2 / 6)`, the eigenvalues
`q + 2 p cos((acos(r) + 2 pi k) / 3)` with `r = det(A - q I) / (2 p^3)`, and the vectors cross
products of rows of `A - l I`, anchored on the more isolated end
([`0057`](./0057-eig3-anchors-the-isolated-end.md)). Against nalgebra's `SymmetricEigen`
([`0056`](./0056-the-routines-d7-does-not-reach.md) decision 4, `0077`) it failed domination on
27 of 31 strata at binary64 and 16 of 19 at binary32: to `1.07e8 u` at `eig:gap-1e-8/top`
against 6.66, and 3683 against 7.66 at `eig:gap-1e-4/top@f32`.

**The loss is the eigenvalues, and it is `acos`.** Scored field by field (the `eigen` metric with
`lambda` and `V` printed apart; a scratch print, not committed), `lambda` alone is the maximum of
every stratum except `eig:gap-1e-8/{bottom,top}`, where `V` reads `1.0e8`–`1.1e8` against
`lambda`'s `6.4e7`, the eigenvalue error carried into the null vector. `lambda`'s maximum grows as
`1 / gap` and stops at `sqrt(u)`: binary64 `48–56`, `327–392`, `2.9e3–4.3e3`, `2.9e4–3.9e4`,
`3.8e5–5.4e5`, `3.7e6–5.9e6 u` at relative gaps `10^-2` to `10^-7`, then `5.2e7–8.4e7`
(`0.5`–`0.9 / sqrt(u)`) from `10^-8` down and at `eig:rank1`. That is the slope of `acos` at
`r = ±1`: `1 - r^2` is of the order `(gap / |A|)^2`, so the rounding of `r` costs `u |A|^2 / gap`
on the close pair, capped at `sqrt(u) |A|` where the slope is infinite. `eig:random`'s 22.7 and
`eig:scale-{up,down}`'s 22.5, 15.5 are the same record (`id` 44 and its scalings), a pair at
relative gap `0.012`. `eig:triple` reads 2.26 and 2.38 (`@f32`), the closed form's own roundings
(`trace / 3` among them), against an oracle that is exact there (`4e-14` and 0).

`0023` (draft) open questions 1 (accept "every input has an answer, non-finite eigenvalues are the
report" and "the frame is orthonormal for every input", or take a mask) and 4 (scale the matrix,
or keep the magnitude cliff, `p^3` normal: `1e-100 < m < 1e100` at `f64`, `1e-12 < m < 1e12` at
`f32`) are unanswered.

## Decision

1. **`eig3` is cyclic Jacobi, four sweeps.** `NUMERICS.md` §17 states it: each sweep rotates the
   pairs `(0, 1), (0, 2), (1, 2)` in that order, each rotation annihilating its entry; the count is
   fixed, with no convergence test. The eigenvalues are the diagonal left, the vectors the product
   of the rotations, each column normalised once at the end, sorted ascending by a compare-exchange
   network whose exchanges flip the first column's sign so that `det = +1`. `PHASE2.md` §6's
   "Smith's closed form" becomes this.
2. **The rotation** (Golub and Van Loan's `sym.schur2`, rearranged so nothing divides by `a_pq`):
   `d = a_qq - a_pp`, `w = 2 a_pq copysign(1, d)`, `den = |d| + sqrt(d^2 + w^2)`, `g = sqrt(den^2 +
   w^2)`, `t = w / den`, `c = den / g`, `s = w / g`; `a_pp -= t a_pq`, `a_qq += t a_pq`, `a_pq = 0`,
   `(a_rp, a_rq) = (c a_rp - s a_rq, s a_rp + c a_rq)`, `(v_p, v_q) = (c v_p - s v_q, s v_p + c
   v_q)`. Where `d^2 + w^2` is below the smallest normal the root takes `1` (the safe argument):
   `c` and `s` share `g`, which must be normal for `c^2 + s^2 = 1` to rounding.
3. **Scaling** (`0023` question 4): the lower triangle is multiplied by a power of two, exactly,
   that brings its largest entry into `[2^-128, 2^128)` (`f64`) or `[2^-32, 2^32)` (`f32`), by
   two-sided steps `2^512, 2^256, 2^128, 2^64` (`2^64, 2^32, 2^16, 2^8`) in `cubic`'s pattern, a
   lazy branch around a cold ladder; the eigenvalues are scaled back. Inside the window
   `d^2 + w^2` cannot overflow and a block whose square is subnormal is below `2^-383 |A|` (`f64`),
   `2^-31 |A|` (`f32`), so the safe argument of decision 2 costs nothing that is scored.
4. **The contract** (`0023` question 1): every input has an answer, no mask. For every finite
   input the eigenvalues are finite unless one exceeds the largest finite number, and the frame is
   orthonormal; a non-finite entry iterates on the zero matrix and gives `NaN` eigenvalues and the
   identity frame. `0023` decisions 1 and 2 are kept in substance; their limits are gone.
5. `0023`'s remaining questions are answered here (1 and 4), and its decision 3 and the limits of
   its decision 4 no longer describe the code: it stays a draft whose history this record closes.

## Rationale

**Measured against the alternatives.** Rows are `cargo xtask conformance --fn eig3 [--precision
f32]`; nalgebra's are `cargo xtask conformance --oracle nalgebra --fn eig3 [--precision f32]`
(nalgebra 0.35.0).

| | binary64 max (stratum) | binary32 max | dominates nalgebra |
|---|---|---|---|
| closed form (before) | 1.07e8 (`gap-1e-8/top`) | 3683 (`gap-1e-4/top`) | 4 of 31, 3 of 19 |
| Jacobi, 3 sweeps | 8.3e8 (`random`) | 3.24 | — |
| 4 sweeps, `c = 1 / sqrt(1 + t^2)`, `s = t c` | 6.43 (`random`) | 4.84 | 31, 19 |
| 4 sweeps, decision 2, no final normalisation | 4.76 (`gap-1e-2/bottom`) | 3.40 | not tabulated |
| **4 sweeps, decision 2, normalised (this)** | **4.76** (`gap-1e-2/bottom`) | **3.24** | **31 of 31, 19 of 19** |
| 5 sweeps, as this | 4.76 | 3.24 | not tabulated |

- **Four sweeps, not three or five.** A convergence study of the rotation above in plain `f64` and
  `f32` (a scratch harness, not committed: `2e5` Haar rotations of each of 25 spectra, gaps
  `10^0` to `10^-16` at either end, near-triples, rank 1 and a graded one; `3e6` of each of four
  random families, uniform and log-uniform gaps and Gaussian entries; and at `f64` a hill climb on
  the entries, 400 restarts) puts the largest off-diagonal entry after three sweeps at
  `2.8e11 u |A|` (`3e-5 |A|`) at `f64` and `330 u |A|` at `f32`, and after four at `2.8e-4 u |A|`
  and `3e-13 u |A|`; the fifth sweep changes no corpus maximum. Three sweeps happen to pass the
  binary32 corpus (64 random records) but not that study, so the count is four at both
  precisions.
- **`c = den / g`, not `1 / sqrt(1 + t^2)`.** One rounding each instead of three; it lowers the
  binary64 maximum from 6.43 to 4.76. Its cost is the normal `g` of decision 2: before the safe
  argument moved from `d^2 + w^2 = 0` to "below the smallest normal", a subnormal block left the
  frame `1.8e6 u` from orthonormal at `f32` (a planted `gap-1e-10/top` spectrum).
- **The final normalisation** takes the frame from `19 u` to `7.7 u` (`|V^T V - I|`, `|det V - 1|`,
  the worst over `10^6` planted cases): the rotations' `c^2 + s^2 - 1` is a length error, which
  orthogonality inherits, and the binary32 maximum moves from 3.40 to 3.24.
- **Kopp's hybrid, not taken.** It keeps the closed form where an error estimate allows, and the
  closed form reads `48–56 u` at a relative gap of `10^-2` and 22.7 on `eig:random`, against
  nalgebra's 5.9–8.7 and 8.2: domination would route nearly every matrix to the fallback, and
  `eig:triple`'s exact oracle is out of reach of `trace / 3`. An `atan2` form with a
  cancellation-free discriminant would mend the `acos` slope but not `trace / 3`; it was not built.
- **The magnitude range comes free.** The rotations need no cube and no `|A|^4`; an exact
  power-of-two scale (decision 3) costs no rounding, where `0023` question 4 priced the division
  by the largest entry at one rounding per entry.

**Cost.** A rotation is two square roots and three divisions on one dependency chain, twelve per
call. Interleaved scratch microbenchmark on `benches/linalg.rs`'s three fixtures, the minimum of
41 runs of `2e5` calls each, one pinned core of an AMD EPYC-Milan VM under load (two runs): `f64`
146 → 389–394 ns (random), 154 → 368 (`gap-1e-8`), 148 → 338 (triple), nalgebra 410, 252, 122;
`f32` 122 → 289–290, 102 → 325, 199 → 289, nalgebra 288–289, 112, 70. That is 1.5–3.2× the closed
form, within 5% of nalgebra on a random spectrum, and slower than nalgebra where its iteration
stops early (a near pair, a multiple of the identity). The final normalisation is 10–27 ns of it.
`cargo xtask bench-gate --bench linalg --only eig3 --against <the closed form's bench binary>`
(six interleaved pairs, one pinned core) reads, as point ratios, `f64` 2.13 (random), 2.43
(`gap-1e-8`), 2.28 (triple) and `f32` 2.62, 2.79, 1.38, the A/A floor beside them 0.3–18%: the
gate fails all six, and is meant to. This record buys the accuracy with that latency.

## Consequences

- `PHASE2.md` §6 and its §0.0 row, `NUMERICS.md` §12's `eig3` row, §13 and the new §17,
  `API.md` R6 (unchanged in substance) follow this record.
- `eig3` calls no transcendental: `solve_cubic` is `Real::acos`'s one library caller
  (`xtask/src/lint/acos.rs`'s exemption stays, its text names one caller).
- The committed rows of `0056` decision 5 are re-recorded: every `eig3` stratum is lower or equal,
  none higher (table below), at both precisions.
- `eig3_tests` records one bound at every gap: `ev`, `res <= 32 u |A|`, `ang_k <= 32 u |A| /
  gap_k`, the frame within `16 u`, with the measured worst ratios in its header; the closed form's
  pinned limits (the tie, the clamp of `r`, the magnitude cliff) are replaced by tests that a
  repeated eigenvalue costs no digits and that every finite magnitude is in the domain.
- A `Dual` derivative is that of the rotations taken. At a multiple of the identity each
  eigenvalue has the derivative of one diagonal entry, where the closed form gave `1/3` of each:
  both are one-sided derivatives of a non-differentiable point.

## Implementation plan

1. Decisions 1–4, `NUMERICS.md` §17 and the documents above, the tests, the re-recorded rows —
   verified by `the_unverified_routines_reproduce_their_committed_rows`, `eig3_tests`,
   `cargo xtask conformance --fn eig3` at both precisions against nalgebra's rows; the cost read
   by `cargo xtask bench-gate --bench linalg --only eig3 --against <the closed form's bench
   binary>`, which fails by the ratios above.

## Open questions

None.

## Measured

`max_u` per stratum, the shipped subject against the committed corpus, before and after, beside
nalgebra's. Domination holds on all 50; `eig:triple` is a tie (`3.99e-14`, 0), which passes.

| stratum | before | after | nalgebra |
|---|---|---|---|
| `eig:random` | 22.7 | **3.53** | 8.19 |
| `eig:gap-1e-{0,1}/{bottom,top}` | 4.60–5.63 | **2.30–2.80** | 6.76–8.78 |
| `eig:gap-1e-2/bottom`, `/top` | 47.9, 55.8 | **4.76**, **2.01** | 5.88, 8.67 |
| `eig:gap-1e-{3,…,7}/{bottom,top}` | 327–5.93e6 | **1.74–2.84** | 4.90–14.9 |
| `eig:gap-1e-{8,…,12}/{bottom,top}` | 5.21e7–1.07e8 | **1.08–1.93** | 6.28–12.8 |
| `eig:triple` | 2.26 | **3.99e-14** | 3.99e-14 |
| `eig:rank1` | 6.33e7 | **1.66** | 13.3 |
| `eig:scale-up`, `-down` | 22.5, 15.5 | **3.07**, **2.30** | 7.95, 7.57 |
| `eig:random@f32` | 6.27 | **2.36** | 6.15 |
| `eig:gap-1e-{0,1}/{bottom,top}@f32` | 5.25–6.69 | **1.79–2.11** | 6.27–6.79 |
| `eig:gap-1e-{2,…,6}/{bottom,top}@f32` | 40.7–3683 | **1.39–1.96** | 6.29–10.4 |
| `eig:triple@f32` | 2.38 | **0** | 0 |
| `eig:rank1@f32` | 3524 | **1.99** | 12.8 |
| `eig:scale-up@f32`, `-down@f32` | 12.0, 51.1 | **1.99**, **3.24** | 8.00, 8.46 |

## Further work

1. A shorter chain per rotation: the next rotation's entry is `c (a_rp - t a_rq)`, and `c^2 =
   1 / (1 + t^2)` is available a square root earlier; worth a measurement if a consumer's budget
   needs the 200 ns back.
2. `0023` can be marked superseded by this record when its owner next edits the index.
