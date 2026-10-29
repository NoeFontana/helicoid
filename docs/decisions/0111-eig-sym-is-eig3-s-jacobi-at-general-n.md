# 0111: `eig_sym` is `eig3`'s Jacobi at general `N`

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #212–#216, the implementation plan below

## Context

locus-calib's M5 D4 is a degeneracy diagnostic. It diagonalises the equilibrated `N_g x N_g` Schur
matrix of the global calibration parameters and flags `λ_min < 1e-5` or `κ > 1e5`; it reports the
eigenvectors of a flagged cluster as the parameters that are poorly determined together.
`helicoid-linalg` has one symmetric eigensolver, `eig3`, at `N = 3`. Nothing else in the workspace
diagonalises a matrix, and `0009` puts the fixed-size linear algebra a consumer needs here.
The owner has set the range at **`4 <= N <= 16`**: `N_g` is at most 9 for one camera, and 16 leaves room
for globals shared across cameras.

**Why Jacobi, again.** `0084` measured the alternative at `N = 3`. Eigen's Householder
tridiagonalisation and QL port read 3.0–11.4 `u` at p99.9 where `eig3`'s Jacobi reads at most 3.9.
At general `N`, Jacobi also has the property tridiagonal QR lacks. For a positive definite
`A = D Â D` (`D = diag(sqrt(a_ii))`, `Â` unit-diagonal), Jacobi with a relative threshold computes
every eigenvalue to `c(N) u κ(Â) λ_k` (Demmel and Veselić, *Jacobi's method is more accurate than
QR*, SIAM J. Matrix Anal. Appl. 13(4), 1992). QR computes only `c(N) u |A|`. The two differ where
`A` is graded: the small eigenvalues of a Schur matrix that was never equilibrated.

**What D4 itself needs.** On the equilibrated matrix, `|Â| <= N`. At `f64` the absolute bound is about
`1e-13`, far below the `1e-5` threshold, so D4's decisions are not limited by rounding even without
the relative property. At `f32` the same bound is about `1e-6`, the threshold's own order: D4 runs
at `f64`. The relative property protects a consumer that does not equilibrate, and the corpus is
built to see it.

**The instrument does not see it today.** `NUMERICS.md` §11's eigenvector metric scores `λ` in
`u |λ|_2`. That is an absolute reading, under which an absolute rotation threshold is as backward
stable as a relative one. A planted absolute threshold would pass every stratum.

## Decision

1. **The function.**
   `helicoid_linalg::eig_sym<S: Real, const N: usize>(a: &Matrix<S, N, N>) -> (Vector<S, N>,
   Matrix<S, N, N>)`, re-exported at the crate root beside `eig3`.
   - `const { assert!(4 <= N && N <= 16) }`. `N = 3` is `eig3`, the one answer at that size. The
     bounds are pinned by two `compile_fail,E0080` doctests.
   - The contract is `eig3`'s: eigenvalues ascending; column `k` a unit eigenvector of value `k`;
     the frame orthonormal with `det = +1`; only the lower triangle read; nothing asserted (R6); a
     non-finite entry gives `NaN` values and the identity frame.
2. **The algorithm** is `NUMERICS.md` §20:
   - `eig3`'s window, ladder and safe arm.
   - Cyclic-by-row sweeps over `p < q`, each rotation `eig3`'s `angle`.
   - A rotation where `|a_pq| > max(u r_p r_q, tiny)`, `r_i = sqrt(|a_ii|)` cached. Elsewhere the
     pair is kept: its entries are selected as they were, not rotated by `(1, 0)`.
   - A small pair, `|2 a_pq| < sqrt(tiny)`, is not `angle`'s alone. Where `|d|` is below
     `sqrt(tiny)` too, `d` and `a_pq` are scaled by `2^512` (`2^64` at `f32`) before `angle` and
     `ta` scaled back. Where `w^2` still underflows, `ta = t a_pq`.
   - The exit: after a sweep that rotates nothing (`Mask::decide`), or at `CAP[N]`.
   - One polish sweep that rotates every pair (`0084` decision 4 at `N`).
   - A stable odd-even transposition sort whose exchanges are counted into column 0's sign.
   - `eig3`'s Newton normalisation of each column.
   - `eig3`'s `angle`, `WINDOW_F64`/`WINDOW_F32` and `unit` are shared, `unit` generic over `N`,
     and so is the scaling rule, `window_scale`, which `into_window` now calls; `eig3`'s bits do
     not move.
3. **The cap**, by `N`: `CAP = [8, 10, 11, 11, 12, 12, 12, 13, 13, 14, 14, 14, 15]` for
   `N = 4..=16`. Each entry is one more than the most sweeps any matrix of the convergence study
   (*Measured*) needed to leave at that `N`, the sweep that rotates nothing counted. The cap is one
   table for both precisions: `f32` never needed more than `f64`.
4. **The corpus**: ids `eig_sym_n{4,9,16}` (`PHASE1.md` §4.3, §4.4), 8.6 MB. That is 0.4 MB past
   `0110`'s estimate of 8.2 MB, the two steep `graded` strata, taken from the 4.4 MiB it leaves.
   - 25 binary64 strata and 20 `@f32` twins: `random`, `gap-1e-k/{bottom,top}`,
     `cluster-{2,3}`, `rank-def-{1,3}`, `graded-1e-k`, `unit-diag-cond-1e-k`, `diagonal`,
     `identity-multiple`, `scale-{up,down}`.
   - 32, 16 and 4 records per stratum at `N = 4, 9, 16`, and 64 in the deep tier.
   - `graded-1e-k` reaches `k = 200` (`30` at `f32`), below `sqrt(tiny)`, where a pair is small.
   - The reference is `mp.eigsy` at 120 digits plus the digits the diagonal spans, so a graded
     `A`'s smallest eigenvalue is as accurate as its largest; plus `kappa = κ_2(D^-1 A D^-1)` where
     the stored `A` is exactly positive definite, else 0.
   - An exact null space of two or more dimensions takes the rational RREF basis, orthonormalised,
     so the reference is one basis at every precision.
5. **The metric** is `Rule::EigenRelative` (`NUMERICS.md` §11's new relative reading). It applies
   where `kappa > 0` and `kappa u <= 2^-10`:
   - per eigenvalue, `|λ̂_k - λ_k| / (κ λ_k u)`;
   - per vector, `min |v̂_k ∓ v_k| min(relgap_k, 1) / (κ u)`;
   - the score is the larger of this reading and `Rule::Eigen`'s.
   - A relative gap above 1 is capped. It would claim an angle below the `u` to which a unit
     vector's components are rounded.
6. **The oracles**: nalgebra's `SymmetricEigen` (`runners/nalgebra`) and Eigen's
   `SelfAdjointEigenSolver` (the Sophus image), at each `N` and both precisions. The bar is `0088`'s
   resolved domination, read on the deep tier, where the committed strata are too small to resolve.
7. **Self-test** (`selftest_linalg`):
   - `rotates-an-eigenvector-at-n` at `eig_sym_n16/eig:random`;
   - `loses-the-small-eigenvalue`, `u |λ|_2` added to the smallest value at
     `eig_sym_n9/eig:graded-1e-16`, one `u` to the absolute reading and `1e16/κ` to the relative
     one.
8. **Tests and rows.**
   - `eig_sym_tests`: exact diagonal and identity cases, power-of-two scaling bit for bit,
     non-finite input, the upper triangle ignored, permutation, the `Dual` value path and lanes,
     first-order perturbation theory with the polish,
     `an_absolute_threshold_loses_the_graded_eigenvalues`,
     `a_block_below_the_safe_arm_keeps_its_digits` and `a_held_lane_keeps_the_sign_of_zero`.
   - `measure_eig_sym` and `hill_climb_eig_sym`, both ignored.
   - `committed_linalg` rows, and bench rows `eig_sym/{4,9,16}/{random,graded,cluster-3}` at `f64`
     and `f32`.

## Rationale

- **Tridiagonal QR** (Eigen's or nalgebra's, ported) is faster: nalgebra reads 0.60 / 3.40 /
  11.4 µs at `N = 4 / 9 / 16` against 0.69 / 6.2 / 28 here (*Measured*). It loses `0084`'s
  accuracy at `N = 3`. On the graded strata it reads up to `5e6 u` against at most 5.4 here: the
  property that answers D4's unequilibrated callers. A call made once per solve buys accuracy
  with latency, as `0081` did.
- **One-sided Jacobi on a Cholesky factor** (Veselić and Hari; Drmač's `dgejsv` preconditioning)
  converges faster on a positive definite input. It needs one, and the `rank-def` strata (the
  gauge, the DLT null space) are not.
- **An absolute threshold**, `|a_pq| <= u |A|`, is `eig3`'s exit. It leaves a graded matrix's small
  eigenvalues to `u |A|`: `an_absolute_threshold_loses_the_graded_eigenvalues` reads `>= 1e6`
  against at most 9 for the relative one.
- **The squared test**, `a_pq² <= u² |a_pp a_qq|`, underflows at `f32` on `graded`: a diagonal of
  `2^-91` squares to `2^-182`. The test would then never mask, and only the cap would end the loop.
  The cached roots cost two square roots per applied rotation, not per test.
- **The floor `tiny`** keeps two exactly zero diagonals (`rank-def`) from rotating on an entry far
  below anything scored. It also bounds a rotating `a_pq` below, which is what lets one exact
  scaling bring every small pair's `w^2` into the normal range.
- **A small pair.** `angle`'s safe arm sets in at `d^2 + w^2 < tiny`, about `2^-511 |A|` at
  `f64` and `2^-63` at `f32`, far above the floor. There it leaves the pair unrotated. Where only
  `w^2` underflows, its `ta` is lost, and `ta` can be as large as the small diagonal. Either is
  harmless under `eig3`'s absolute exit and fatal to the relative one. Before the fix,
  `graded-1e-30` at `f32` read `1.9e6 u` relative and a `2^-70` block `5.6e6`.
  - The scaling is exact, and `angle` is invariant under it.
  - It sits behind a branch on `|2 a_pq|`, out of line. A scalar pays only the comparison, and
    the latency was measured unchanged.
  - `ta = t a_pq` everywhere would change every rotation's bits from `eig3`'s, against `0084`'s
    reason for reading `ta` from `w^2`.
- **A kept pair is selected, not rotated by `(1, 0)`.** `1 x - 0 y` is `+0` where `x = -0` and
  `y < 0`, so a lane held beside a rotating one would part from its scalar run in the sign of a
  zero, which a later rotation's sign can carry.
- **The mirror in a second pass.** One pass that also stores the strided row measured 5–7%
  slower at `N = 9, 16`: it keeps LLVM from vectorising the column.
- **A cap per `N`**: the worst case grows from 7 sweeps at `N = 4` to 14 at 16. One constant would
  be 15 everywhere. The cap only bounds lanes and a pathological input; a scalar leaves at
  convergence.
- **The polish** is `0084`'s reason at `N`. Without it, a diagonal input with an off-diagonal lane
  returns vectors with no derivative (`derivative_is_first_order_perturbation_theory`'s third
  case).
- **A stable transposition network**, not Batcher's (fewer comparators). Stability is what returns a
  multiple of the identity as `(cI, I)` and a diagonal input in its axes' order, bit for bit. At
  `N = 16` its 120 exchanges cost little beside the sweeps.
- **Newton–Schulz on the frame**, `V(3I - VᵀV)/2`, would bring the frame to `~u`, against 20.7 `u`
  measured at `N = 15`. The corpus metric reads the columns and finds no stratum above an oracle,
  so it is further work, not taken.
- **Round-robin orderings** (disjoint pairs per step, for instruction-level parallelism) change the
  bits and need their own measurement. Nothing D4 does is latency-bound.
- **`N <= 16`, not 20**: the owner's choice. A record that widens it owes the convergence study and
  the corpus at the new `N`.

## Consequences

- `NUMERICS.md` gains §20 and §11's relative reading, and §0 and §12 rows. `API.md` §2 and R6 name
  `eig_sym`. `PHASE2.md` §6 and its §0.0 row, and `PHASE1.md` §4.3 and §4.4, list it.
- The three ids enter the determinism digests and `committed_linalg`'s no-regress rows.
- `eig3`'s committed rows and digests are unchanged by the shared helpers. That is the bar for the
  refactor. `eig3`, `svd3` and `nearest_rotation` latency is unchanged: `window_scale` takes the
  arm that scales the entries as a closure, which only that arm runs.
- D16: `+ - * /`, `sqrt` and `abs` only, the exit a comparison of values every target computes to
  the bit, and no fused multiply-add. LLVM may vectorise the column updates, which reassociates
  nothing.

## Implementation plan

1. This record and its spec text. Verified by `just lint` and `just doc`.
2. The instrument (`0006`):
   - the generator's `eig_sym` module, its checks and tests;
   - the corpus and its deep tier;
   - `Rule::EigenRelative`, `coverage::required`, the oracle rows.
   Verified by `just corpus-check`, `just corpus-test`, `just corpus-deep-check`,
   `cargo nextest run -p xtask`, `just oracle-nalgebra` and `just oracle-sophus`.
3. The kernel: `eig_sym.rs`, `eig_sym_tests.rs`, the shared helpers in `eig3.rs`, the shipped
   subject, the self-test defects, the committed rows and the bench rows. Verified by:
   - `just test`;
   - `just conformance` at both precisions and `just conformance --self-test`;
   - `just conformance-deep`, then `just envelope --deep`;
   - `just no-std`, `just msrv`, `just wasm` and `just determinism`.

## Open questions

None.

## Further work

1. Newton–Schulz on the frame (above), if a consumer reads `VᵀV`.
2. A values-only entry, as `eig3_values` (`0104`), if D4's profile asks: the frame is about half of
   the work at `N = 16`.
3. Round-robin orderings, if a consumer puts `eig_sym` on a hot path.

## Measured

**Convergence and accuracy** (`measure_eig_sym`; `2e4` matrices of each of 20 families at every
`N`, the first 300 scored against the double-double relative Jacobi reference; `--release`, seed
1). The worst over the families. `rel` is read where `κ u <= 2^-10`:

| N | sweeps f64 / f32 | abs f64 / f32 | rel f64 / f32 | frame f64 / f32 |
|---|---|---|---|---|
| 4 | 7 / 6 | 4.71 / 3.99 | 1.80 / 1.54 | 7.50 / 5.36 |
| 5 | 9 / 8 | 5.29 / 4.86 | 2.08 / 2.09 | 8.50 / 6.81 |
| 6 | 10 / 9 | 5.25 / 5.67 | 3.28 / 1.84 | 10.20 / 8.27 |
| 7 | 10 / 9 | 6.74 / 5.47 | 3.13 / 1.93 | 10.88 / 9.54 |
| 8 | 11 / 9 | 6.00 / 5.67 | 4.25 / 1.61 | 12.00 / 12.33 |
| 9 | 11 / 10 | 7.86 / 6.40 | 4.20 / 2.81 | 12.19 / 11.42 |
| 10 | 11 / 10 | 7.89 / 6.99 | 2.50 / 1.67 | 16.12 / 14.44 |
| 11 | 12 / 10 | 9.52 / 7.48 | 2.78 / 2.26 | 14.53 / 13.00 |
| 12 | 12 / 11 | 9.11 / 6.49 | 3.03 / 1.92 | 14.00 / 12.69 |
| 13 | 13 / 11 | 10.35 / 8.03 | 3.24 / 2.30 | 18.42 / 13.57 |
| 14 | 13 / 12 | 9.68 / 7.71 | 2.74 / 2.15 | 17.00 / 14.63 |
| 15 | 13 / 12 | 9.68 / 7.95 | 3.29 / 2.49 | 20.66 / 14.41 |
| 16 | 14 / 12 | 12.63 / 8.96 | 3.15 / 2.56 | 19.50 / 17.30 |

- The slowest families are `rank-def-3` and `unit-diag-1e-12`, a near-singular coupling the
  relative threshold resolves to `u r_p r_q`.
- The hill climb (`hill_climb_eig_sym`, 150 restarts of 300 shrinking perturbations from each
  family) finds at most 7, 10 and 13 sweeps at `N = 4, 9, 16`, none above the random study.
- With the cap lifted, no matrix needs more sweeps than the table's.
- The families with small pairs, `graded-1e-30` at `f32` and `graded-1e-200` at `f64`, read at
  most 2.21 `rel`. Before the small-pair rotation they read `2.6e6` and `1.2e15`. Every other
  row is the same to the printed digit, except `diag-tie+1e-17`'s `f32` frame, within 0.5 `u`:
  its polish now rotates small pairs.
- `eig_sym_tests`' recorded bounds are twice these worsts: `abs <= 25.3`, `rel <= 9`,
  `frame <= 42`.

**Corpus, committed tier** (`cargo xtask conformance`, the largest `max_u` over each id's strata):

| id | helicoid f64 / f32 | nalgebra 0.35 f64 / f32 | Eigen 3.4 f64 / f32 |
|---|---|---|---|
| `eig_sym_n4` | 3.66 / 3.01 | 2.48e19 / 1.50e19 | 2.48e19 / 1.86e10 |
| `eig_sym_n9` | 5.35 / 5.70 | 5.10e15 / 3.32e7 | 9.26e19 / 1.23e10 |
| `eig_sym_n16` | 7.83 / 5.66 | 9.57e16 / `NaN` | 4.11e17 / 8.22e7 |

- The oracles' maxima are on `graded-1e-200` and `graded-1e-30@f32`, where `helicoid` reads at
  most 1.30 and 1.50. nalgebra returns `NaN` on all four `eig_sym_n16/graded-1e-30@f32` records.
- On every one of the 135 (id, stratum, precision) cells, `helicoid`'s maximum is at or below the
  best oracle's (median ratio 0.37). The ties are `diagonal` and `identity-multiple`, which every
  program solves exactly.
- **`just envelope --deep`**, every oracle rerun on both tiers: no `eig_sym` stratum is a
  domination failure, and none is printed as an excess, resolved or not, at either tier. The run's
  one failure is `main`'s standing "no baseline for `helicoid`".

**Latency** (`benches/linalg.rs`, `f64`, one unpinned run on a shared host; nalgebra's
`SymmetricEigen::new` on the same fixtures from a scratch binary, best of 31 × 2000 calls):

| N | `random` | `graded` | `cluster-3` | nalgebra `random` / `graded` |
|---|---|---|---|---|
| 4 | 0.69 µs | 0.45 | 0.39 | 0.60 / 0.46 |
| 9 | 6.2 | 3.9 | 6.8 | 3.40 / 2.57 |
| 16 | 28.0 | 17.6 | 28.2 | 11.4 / 8.9 |

- The small-pair branch, the selected kept pair and `r` left out of the polish cost nothing on
  the scalar path. Against the first implementation, paired, one binary per variant, pinned:
  2–9% faster at `f32` and about 4% at `f64`, at every `N`.
