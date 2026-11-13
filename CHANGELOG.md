# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Released
entries live in [`docs/changelog/`](./docs/changelog/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [Unreleased]

### Breaking

- **`Jac` gains two required methods, `add` and `scale`**, entry by entry and named, never an
  operator ([`0112`](docs/decisions/0112-a-jacobian-sums-by-name.md)). Every impl in the crate has
  them and stays structured; an external `Jac` impl must add both. `NUMERICS.md` §14 gains the
  dense twins `reference::{sen3jac, sim3jac}_{add, scale}`. Every existing output keeps its bits.

### Added

- **`helicoid-spline`**, Phase 7's crate ([`0113`](docs/decisions/0113-phase-7-opens-for-a-calibration-consumer.md),
  [`0114`](docs/decisions/0114-a-spline-evaluates-the-window-it-is-handed.md),
  [`0115`](docs/decisions/0115-the-spline-tables-are-exact-integers.md)): uniform cumulative
  B-splines of order `K = 2..=6` on any `LieGroup`, under `helicoid`'s budget. `Segment` evaluates
  the value, body velocity and body acceleration by Sommer et al.'s recursion; `JacSegment` their
  Jacobians on every control, either side, in `O(K)`; `derivative_gram` the energy matrices. Every
  table is exact integers at compile time. On the 16 new corpus ids it reads at most 9.69 `u`
  (binary64) and 8.50 `u` (binary32), and bit-identical across x86_64 and `wasm32-wasip1`.

### Instrument

Not shipped in the crates.
- **The spline's corpus ids, references and seeded defects** ([`0116`](docs/decisions/0116-the-spline-is-measured-before-it-ships.md)):
  `{so3, rn, se3}_spline_*` and `rn_spline_gram`, 1.75 MiB of the 80 MiB budget, with mpmath
  references and an independent `NUMERICS.md` §21 evaluator (`seeded:spline`) that all six planted
  defects are measured against. The shipped subject answers every id; `bench --bench spline` adds
  the criterion rows.
- **`bench-gate` reads criterion's samples where criterion writes them**: `$CRITERION_HOME`, else
  `$CARGO_TARGET_DIR/criterion`, else cargo's target directory, where it read
  `<root>/target/criterion` and found no samples under a moved target directory. A relative
  variable is refused, since the bench binary and cargo would resolve it from different
  directories. Found in a consumer that copied the gate.
- **The bench tests write to per-process scratch directories** removed on drop, where they shared
  fixed `/tmp/helicoid-bench-*` paths, so concurrent test runs in two worktrees deleted each
  other's fixtures; so does the `untracked` lint's git test (`/tmp/helicoid-lint-untracked`).
  Nothing breaks.

## [0.0.9] - 2026-10-29

The general symmetric eigensolver release. `helicoid_linalg::eig_sym::<S, N>` is `eig3`'s Jacobi at `4 <= N <= 16` with a relative threshold: a positive definite input's eigenvalues are accurate relative to each one, down to a diagonal spanning `1e-200`, and on every committed stratum it reads at or below nalgebra's `SymmetricEigen` and Eigen's `SelfAdjointEigenSolver` (at most 7.8 `u`), at 1.2–2.5× nalgebra's latency. The corpus gains `eig_sym_n{4,9,16}` and the `real_atan2` strata where `math::atan2`'s formulas meet; `corpus-check` proves each thing once, `determinism (compare)` is a required check, and the corpus budget is 80 MiB. **No API break**, and every existing output keeps its bits.

### Instrument

Not shipped in the crates.
- **`real_atan2` strata where `math::atan2`'s formulas meet** (`0105`): the centre midpoints, the
  binary32 fold, the diagonal, the axes, the window's low edge and both ends of the derivative's
  domain, 896 records. `math::atan2`'s value reads at most 0.939 `u` on them at `f64` (`libm`
  1.579) and 0.981 at `f32` (1.679).
- **`corpus-check` proves each thing once** (`0108`): `just corpus-verify` checks the manifest
  and every file's hash in seconds; CI skips the regeneration on a pull request whose generator,
  interpreter and manifest a full check already passed, and otherwise replays unchanged tasks from
  a memo keyed on the code each ran, 2% of them run again. `main` and a nightly job run every task.
  The generator's tests run a test class per core in their own job. Locally `just corpus-check`
  takes about a minute with the memo warm, against 17.
- **`determinism (compare)` is a required check** (`0109`): the comparison of the x86_64,
  aarch64 and `wasm32-wasip1` digest tables now blocks a merge, and it runs and fails when a digest
  leg did not succeed, where a skipped job would have reported success.
- **The corpus budget is 80 MiB and `cargo xtask lint` checks it** (`0110`, amending `0071`):
  the committed corpus had reached 67.8 MiB, past `0071`'s 64, with nothing failing.
- **`eig_sym_n{4,9,16}`** (`0111`): 25 binary64 and 20 binary32 strata per id, among them
  graded positive definite matrices down to a diagonal spanning `1e-200` (`1e-30` at `f32`),
  unit-diagonal ones, exact null spaces and exact clusters; 32, 16 and 4 records per stratum
  (8.6 MB), 64 in the deep tier. The reference adds the digits the diagonal spans to its
  precision and carries `kappa`, and `NUMERICS.md` §11 gains a relative reading on positive
  definite records, which an absolute one cannot see: nalgebra and Eigen read up to `9e19 u` on
  `graded-1e-200` under it. Both answer the ids as oracles.

### Added

- **`helicoid_linalg::eig_sym::<S, N>`** (`0111`): the eigendecomposition of a symmetric `N x N`
  matrix, `4 <= N <= 16`, with `eig3`'s contract. `eig3`'s Jacobi with a relative threshold, so a
  positive definite input's eigenvalues are accurate relative to each one, a measured sweep cap per
  `N` and a polish sweep for `Dual`; a pair below the square root of the smallest normal is
  scaled by an exact power of two first, so the relative accuracy holds down to a diagonal
  spanning `1e-200`. On every committed stratum of `eig_sym_n{4,9,16}` it reads at or below
  nalgebra's `SymmetricEigen` and Eigen's `SelfAdjointEigenSolver` (at most 7.8 `u`), at 1.2–2.5×
  nalgebra's latency (0.69, 6.2 and 28 µs at `N = 4, 9, 16`). `eig3`'s bits and latency are
  unchanged.

## Released versions

Frozen history, newest first, in [`docs/changelog/`](./docs/changelog/); leave it
out of searches and sweeps.

- [`0.0.8`](./docs/changelog/0.0.8.md) — 2026-10-24 (the portable `atan2`, `eig3_values`)
- [`0.0.7`](./docs/changelog/0.0.7.md) — 2026-10-22 (`nearest_rotation`, `inv_act` and its Jacobians)
- [`0.0.6`](./docs/changelog/0.0.6.md) — 2026-10-19 (`svd3` in the signed convention)
- [`0.0.5`](./docs/changelog/0.0.5.md) — 2026-10-16 (closed-form SE_N(3) Jacobians, `Exp` near `π`, the deep tier)
- [`0.0.4`](./docs/changelog/0.0.4.md) — 2026-10-14 (`from_matrix` on the unit sphere, `round_from_f64`, `just determinism`)
- [`0.0.3`](./docs/changelog/0.0.3.md) — 2026-10-13 (Phase 5, `eig3` as Jacobi, three C++ oracles)
- [`0.0.2`](./docs/changelog/0.0.2.md) — 2026-10-09 (the carried path and a faster geodesic)
- [`0.0.1`](./docs/changelog/0.0.1.md) — 2026-10-08 (first release)

[0.0.1]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.1
[0.0.2]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.2
[0.0.3]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.3
[0.0.4]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.4
[0.0.5]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.5
[0.0.6]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.6
[0.0.7]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.7
[0.0.8]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.8
[0.0.9]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.9
