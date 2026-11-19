# 0118: `eig_sym` reaches `N = 64`

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** the implementation plan below

> Amends [`0111`](./0111-eig-sym-is-eig3-s-jacobi-at-general-n.md) (the bound, the cap, the corpus
> ids) and [`0110`](./0110-the-corpus-budget-is-80-mib-and-checked.md) (the cap is 100 MiB).

## Context

`0111` set `eig_sym`'s range at `4 <= N <= 16` for a calibration consumer's degeneracy diagnostic,
and said what widening it owes: "the convergence study and the corpus at the new `N`".

The same diagnostic now runs on larger Schur matrices. The matrix of global parameters is 9 x 9 in
the consumer's base model and 24 x 24 in its next one. It reaches about 60 once time offset, skew
and the camera intrinsics are estimated globally.

Nothing in `eig_sym` depends on `N` but the sweep cap. The rotation, the relative threshold, the
small-pair scaling, the polish, the sort and the normalisation are generic (`NUMERICS.md` §20). D1
asks for one implementation, so the route is the same Jacobi with a longer cap table, not a second
solver above 16.

**The corpus cannot hold it under `0110`.** The tracked corpus is 82 363 319 bytes (78.55 MiB) on
`main` at `f71a88d`, which leaves 1.45 MiB. `0110`'s size formula, about `60 N^2 + 38 N + 138`
bytes a record, matches `eig_sym_n16`'s 2 904 844 bytes exactly. At 24, 40 and 64 a record is
35.6, 97.7 and 248 KB. A single record of each of the 45 strata at `N = 64` is 11.2 MB.

## Decision

1. **The bound** is `4 <= N <= 64`: `const { assert!(4 <= N && N <= 64) }`. The two
   `compile_fail,E0080` doctests are `N = 3` and `N = 65`. Every other clause of `0111` decision 1
   stands.
2. **The algorithm** is `0111` decision 2, unchanged. At `N <= 16` every bit is `0111`'s.
3. **The cap** extends to `N = 64`: `CAP` has 61 entries, the first 13 `0111`'s. Each entry is one
   more than the most sweeps any matrix of the study (*Measured*) needed at that `N` **or below**,
   the sweep that rotates nothing counted: 15 at 17, 18 at 24, 21 at 40, 24 at 64.
   - "Or below": a matrix bordered by a zero row and column is a matrix at `N + 1` that takes the
     same rotations. A zero `a_pq` is never rotated (the threshold is at least `tiny`), and neither
     the scaling nor the order of the other pairs moves. The worst case cannot fall as `N` grows,
     so a dip in the measured maxima (16 sweeps at 23, after 17 at 22) is the sample's.
   - At `N <= 16` the rule gives `0111`'s table, entry for entry.
   - `f32` needs at least one sweep fewer than `f64` at every `N`, so the table is `f64`'s.
4. **The corpus**: ids `eig_sym_n{24,40,64}`, `0111` decision 4's 45 strata and generators at
   each, with the same `GAP_K`, `GRADED_K` and `UNIT_DIAG_K`.
   - **2, 1 and 1 records per stratum**, about 3.2, 4.4 and 11.2 MB, 18.8 MB in all.
   - In the deep tier, **64, 32 and 16**, not `0088`'s flat 64. A swap test on `n` records reaches
     `2^-n`, so 16 resolves `p <= 0.001`. A flat 64 at `N = 64` would be 715 MB of git-ignored
     data, and hours of `mp.eigsy`.
   - `strata.deep(n, m)` takes the deep count per id; `m` defaults to 64, so no other id moves.
5. **The budget** (amends `0110` decision 1): the committed corpus stays under **100 MiB**
   (104 857 600 bytes), checked by `lint::budget` as before. With this record it is about
   96.5 MiB.
6. **The oracles**: nalgebra's `SymmetricEigen` and Eigen's `SelfAdjointEigenSolver` at the three
   `N` and both precisions. The bar is `0088`'s resolved domination on the deep tier, as `0111`
   decision 6.
7. **Tests.**
   - `eig_sym_tests`' recorded bounds take a second band for `17 <= N <= 64`: twice the worst
     the study reads there, `[abs, rel, frame] = [49.2, 12.4, 128.0]` from 24.58 `u` at `N = 57`,
     6.18 at 55 (`f32`) and 64.0 at 60.
   - `every_family_is_within_its_bounds_*` at `N = 24, 40, 64` and a sparse set between them.
   - `measure_eig_sym` covers `4..=64` on a thread pool. `hill_climb_eig_sym` runs at 4, 9, 16, 24,
     32, 40, 48 and 64.
   - `committed_linalg` rows and bench rows `eig_sym/{24,40,64}/random` at `f64` and `f32`.
8. **The rustdoc** states the cost (*Measured*): latency, idle and beside a busy hardware thread,
   and the stack, about `4 N^2 size_of::<S>()`, four times that unoptimised.
9. **Row `p` is mirrored once per row of pairs**, after the last pair `(p, q)`, not after each.
   No rotation of the row reads it: column `q`'s entry at row `p` only feeds turned values that
   each rotation then assigns (`a_pp`, `a_qq`, `a_pq`, `a_qp`). Every entry is computed from the
   same operands in the same order, so every output keeps its bits at every `N`; row `q` is still
   mirrored per rotation, as later pairs read it. It halves the strided stores: 3% faster at
   `N = 9`, 9 to 16% from 12 to 64, and slower nowhere (*Measured*).

## Rationale

- **One solver, not a switch.** A tridiagonal QR above 16 would be faster: nalgebra's is 4.3
  times faster at 64 on an idle core, up to 11 times beside a busy hardware thread (*Measured*). It
  is a second implementation (D1), and it loses the relative accuracy that `0111` argued for:
  nalgebra reads up to `5e6 u` on `graded` at `N <= 16`. The diagnostic runs once per solve, so its
  latency, 1.3 to 4 ms at 64, is not on a hot path.
- **The gap is the sweeps, not the code.** Idle, a rotation costs about 0.7 cycles for each entry
  it updates, vectorised; padding the columns against cache-set aliasing changed nothing, and
  dropping the frame's accumulation entirely saves at most 22% idle, 27% contended. Jacobi's 20-odd
  sweeps of `N^3` work against QR's few are the price of the relative accuracy.
- **Only the stores that keep bits, and only the simplest.** Decision 9 removes stores nothing
  reads. Going further, mirroring row `q` only into the columns the sweep still reads and
  restoring the rest once per sweep, keeps every bit too and is a third faster at `N = 64`. But
  its loops lose their constant trip count: it is 15% slower at `N = 6` and ties decision 9 from
  16 to 32. Taking it above a measured `N` would add a second scheme and a switch point to buy 2 to
  17% on a diagnostic run once per solve (*Further work*).
- **A cap per `N`, measured at every `N`.** The cap bounds lanes and pathological inputs. A cap
  below a matrix's need returns it unconverged, so each entry is measured, not fitted; the running
  maximum (decision 3) is the measurement read with what the kernel proves, not a fit.
- **64, not 60.** 60 is the consumer's estimate today. 64 leaves room for a global parameter more,
  and costs nothing beyond the study and the ids.
- **24, 40 and 64 as the ids.** 24 and 64 are the consumer's sizes. 40 is a midpoint, so a defect
  specific to one monomorphisation, such as a vectorisation choice LLVM makes at one width, has two
  chances to show between 16 and 64. The study covers every `N` between them.
- **All 45 strata at `N = 64`, at one record each.** One record is the least a stratum can hold,
  and every stratum exercises a path: `graded-1e-200` the small pair, `rank-def` the floor,
  `scale-{up,down}` the window. The determinism digests are per monomorphisation, and these are
  new ones.
- **Raise the cap, not cut the strata.** `0110`'s reasons for the cap are clone size and
  `corpus-check`'s time. The new records are full-width digits, which compress less than the rest:
  `xz -9` takes their 18.8 MB to 4.9 MB, a quarter. On an 8-vCPU host under other load,
  `corpus-check` takes under 2 minutes warm, as pull requests run it, and 17 fresh, as `main` and
  the nightly run it. Leaving `n40` and `n64` to the deep tier alone would ship two
  monomorphisations that CI never runs (`0006`).

## Consequences

- At `N <= 16`, `eig_sym`'s bits, its committed rows and its determinism digests are unchanged.
  That is the bar for the change.
- `NUMERICS.md` §20 and its §12 row, `API.md` §2, `PHASE2.md` §6 and its §0.0 row, and `PHASE1.md`
  §4.3 and §4.4 state `4 <= N <= 64` and the new ids. `PHASE1.md` §4.4's cap sentence says
  100 MiB.
- `0111` and `0110` take an amendment banner that cites this record.
- `corpus-deep` takes about 2.4 times as long: 178 minutes on an 8-vCPU host under other load,
  from about 75.

## Implementation plan

1. This record and its spec text. Verified by `just lint` and `just doc`.
2. The instrument (`0006`):
   - the generator's ids, `strata.deep(n, m)` and their tests;
   - the corpus and its deep tier;
   - `lint::budget` at 100 MiB;
   - the oracle, coverage and metric rows;
   - nalgebra's and the Sophus image's instantiations.
   Verified by `just corpus-check`, `just corpus-test`, `just corpus-deep-check`,
   `cargo nextest run -p xtask`, `just oracle-nalgebra` and `just oracle-sophus`.
3. The kernel: the bound, `CAP`, the row mirror (decision 9), the rustdoc, `eig_sym_tests`, the
   shipped subject, the committed rows and the bench rows. Verified by:
   - the determinism digests of the kernel before decision 9 and after, identical for every id;
   - `just test`;
   - `just conformance` at both precisions and `just conformance --self-test`;
   - `just conformance-deep`, then `just envelope --deep`;
   - `just no-std`, `just msrv`, `just wasm` and `just determinism`.

## Open questions

None.

## Further work

1. `0111`'s further work, which grows with `N`: Newton–Schulz on the frame if a consumer reads
   `VᵀV`, and a values-only entry if the diagnostic's profile asks.
2. From `N = 32`, `a` and `v` (`2 N^2` entries) outgrow half a core's L1, and beside a busy
   hardware thread the time grows to 1.4 times idle at 40 and 2.3 at 64. Two schedules keep every
   bit and are owed a consumer whose profile shows it: row `q` mirrored only into the columns the
   sweep still reads, with one restoring pass per sweep, above `N = 20` (measured: 0.67 of the
   kernel before decision 9 at 64, against decision 9's 0.84; L1 write misses down 85%); and `v`
   turned by blocks of rows after each row of pairs.

## Measured

**The convergence study**: `measure_eig_sym --release`, `2e4` matrices of each of the 20 families
at every `N` in `4..=64` and both precisions, 300 of them scored against the double-double
reference, seed 1. **The hill climb**: `hill_climb_eig_sym --release`, 8 restarts of 300 steps
per family, `f64`, at 4, 9, 16, 24, 32, 40, 48 and 64. Sweeps are maxima over the families, the
sweep that rotates nothing counted; abs, rel and frame are maxima in `u`. Every `N <= 16` row
reproduces `0111`'s. The most sweeps came from `unit-diag-1e-12` at every `N >= 11`.

| `N` | sweeps `f64` | `f32` | hill climb | `CAP` | abs `f64` / `f32` | rel | frame |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 7 | 6 | 7 | 8 | 4.71 / 3.99 | 1.80 / 1.54 | 7.5 / 5.4 |
| 9 | 11 | 10 | 9 | 12 | 7.86 / 6.40 | 4.20 / 2.81 | 12.2 / 11.4 |
| 16 | 14 | 12 | 12 | 15 | 12.63 / 8.96 | 3.15 / 2.56 | 19.5 / 17.3 |
| 17 | 14 | 13 |  | 15 | 9.94 / 9.64 | 3.16 / 2.55 | 22.0 / 17.2 |
| 20 | 15 | 14 |  | 16 | 10.92 / 9.35 | 3.30 / 3.17 | 21.0 / 18.4 |
| 22 | 17 | 13 |  | 18 | 15.66 / 10.11 | 2.95 / 2.69 | 23.4 / 18.1 |
| 24 | 17 | 14 | 14 | 18 | 13.36 / 10.44 | 3.27 / 3.40 | 26.0 / 17.4 |
| 28 | 19 | 14 |  | 20 | 13.05 / 10.45 | 3.48 / 3.46 | 27.0 / 20.5 |
| 32 | 18 | 15 | 16 | 20 | 15.50 / 12.52 | 3.29 / 3.79 | 32.0 / 22.8 |
| 36 | 20 | 17 |  | 21 | 14.13 / 11.77 | 3.68 / 3.89 | 30.0 / 23.0 |
| 40 | 19 | 16 | 17 | 21 | 16.61 / 14.50 | 4.40 / 4.49 | 34.0 / 28.2 |
| 48 | 20 | 16 | 18 | 21 | 16.19 / 14.57 | 5.34 / 4.92 | 40.0 / 30.2 |
| 51 | 21 | 16 |  | 22 | 15.71 / 16.01 | 4.13 / 5.10 | 42.0 / 34.9 |
| 56 | 21 | 16 |  | 22 | 14.59 / 18.66 | 4.36 / 5.25 | 52.0 / 34.5 |
| 57 | 21 | 16 |  | 22 | 24.58 / 16.72 | 4.70 / 4.83 | 48.0 / 37.0 |
| 60 | 21 | 17 |  | 22 | 16.53 / 17.53 | 4.54 / 5.97 | 64.0 / 36.8 |
| 62 | 22 | 16 |  | 23 | 17.30 / 15.32 | 5.08 / 5.11 | 53.0 / 39.2 |
| 64 | 23 | 17 | 19 | 24 | 18.21 / 17.58 | 5.41 / 5.17 | 56.0 / 42.4 |

The `f64` sweeps at every `N` in `17..=64`: 14 15 15 15 16 17 16 17 17 17 17 19 18 18 18 18 18 19
19 20 19 19 19 19 20 19 20 20 20 20 20 20 20 20 21 20 21 21 21 21 21 21 21 21 21 22 21 23. With
the cap lifted to 30, no matrix of the study or the climb needed more than the table's.

**Latency**, of the kernel before decision 9: `f64`, a random symmetric matrix, the least of 200 × 5 calls, one pinned core of an
8-vCPU EPYC Milan (Zen 3), interleaved with nalgebra 0.35's `SymmetricEigen` on the same matrix.
"Idle" is the host at load 6, "busy" at load 17 with the core's other hardware thread loaded.

| `N` | idle | nalgebra | ratio | busy | nalgebra | ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 30.8 µs | 10.2 µs | 3.0 | 35.5 µs | 11.8 µs | 3.0 |
| 24 | 94.5 µs | 24.2 µs | 3.9 | 113 µs | 29.3 µs | 3.8 |
| 40 | 389 µs | 99.7 µs | 3.9 | 538 µs | 105 µs | 5.1 |
| 64 | 1.52 ms | 345 µs | 4.4 | 3.4–3.9 ms | 344 µs | 9.3–11.3 |

Idle, the time per `N^3` is flat, 5.8–7.5 ns from 16 to 64. Without the frame's accumulation it
falls by 3% at 16, 13% at 32 and 22% at 64 idle, 27% at 64 busy.

With decision 9, at load 3 and the same protocol: 27.5 µs, 88.5 µs, 358 µs and 1.29 ms at
`N = 16, 24, 40, 64`, against nalgebra's 10.1, 24.3, 84.1 and 302 µs: 2.7, 3.6, 4.3 and 4.3 times.

**Decision 9 and the further schedule**, against the kernel before decision 9: the median of 25
interleaved rounds of the least of 3 × 6 to 60 calls, `f64`, one pinned core at load 12 to 15,
outputs equal to the bit. Callgrind's simulated L1 write misses at `N = 64`: 2.32 M for three calls
before, 1.44 M with decision 9, 0.35 M with the further schedule.

| `N` | 4 | 6 | 9 | 12 | 16 | 20 | 24 | 32 | 40 | 64 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| decision 9 | 1.00 | 1.00 | 0.97 | 0.88 | 0.88 | 0.87 | 0.89 | 0.89 | 0.86 | 0.84 |
| further | 1.02 | 1.15 | 1.03 | 0.94 | 0.89 | 0.87 | 0.87 | 0.88 | 0.84 | 0.67 |

**Stack**, measured on a spawned thread: `f64` at `N = 64` needs 160 KiB optimised, 640 KiB
unoptimised; `Dual<f64, 6>` 1 MiB and 4 MiB.

**The corpus** against the oracles. On the deep tier, `eig_sym` is at or below the better of
nalgebra and Eigen on every stratum at `N = 24, 40, 64`, at both precisions, and `just envelope
--deep` fails none. Its maxima there are 13.0, 15.2 and 13.5 `u` in binary64 and 8.26, 9.37 and
13.6 in binary32. The closest stratum is `unit-diag-cond-1e-12` at 64: 13.5 `u` against
nalgebra's 14.0.

On the committed tier no stratum fails. Four are unresolvable at one record, and the deep tier
reverses every one, because the oracles' maxima rise with the records:

| `N`, stratum | committed, `eig_sym` / nalgebra | deep, `eig_sym` / nalgebra / Eigen |
|---|---:|---:|
| 40, `unit-diag-cond-1e-12` | 9.71 / 3.16 | 11.2 / 18.1 / 30.7 |
| 64, `unit-diag-cond-1e-12` | 11.3 / 8.78 | 13.5 / 14.0 / 47.1 |
| 40, `gap-1e-1/bottom@f32` | 6.12 / 5.31 | 6.77 / 13.7 / 29.6 |
| 40, `scale-down@f32` | 5.75 / 5.10 | 7.62 / 15.1 / 36.5 |
