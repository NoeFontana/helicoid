# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Released
entries live in [`docs/changelog/`](./docs/changelog/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [Unreleased]

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

## [0.0.8] - 2026-10-24

The portable `atan2` release. `Real::atan2` is now `helicoid`'s own kernel at both precisions, public as `helicoid_linalg::math::{atan2, atan2f}`: within 0.504 ulp at `f64` and correctly rounded at `f32` on every draw read, where `libm` reads 1.40 and 1.47, with the same bits on every target whatever `libm` a consumer resolves. `SO3::log` rounds each component once on its exact arm, and with the new `atan2` its maximum falls from 2.317 to 1.405 `u`; a composed `Log` stays inlinable. `eig3_values` returns `eig3`'s eigenvalues without the frame in up to half the time. **No API break**; every `atan2` and `Log` caller's bits move.

### Changed

- **`Real::atan2` is `helicoid`'s own portable kernel at both precisions** (`0105`), where it was
  the `libm` crate's. Binary64 reads within 0.504 ulp of the exact value on every draw read
  (`libm` 1.40), binary32 is correctly rounded on every draw read (`atan2f` 1.47); IEEE 754's
  limits at zeros and infinities, and the same bits on every target independent of `libm`. It
  costs 1.1–1.25× `libm`'s latency, 1.07× on `SO3::log`'s exact arm and nothing on its series
  arms. Every `atan2` caller's bits move: `so3_log`'s maximum falls from 1.680 to 1.405 `u`, and
  `Log` of SO(2), SE(2), SE_N(3) and Sim(3), `S2::local` and `coeffs`' `r` move with it. The
  regenerated sweep moves `r`'s series switch at both precisions (`f32`'s series is one term
  shorter).

- **`SO3::log` rounds each component once on its exact arm** (`0106`): `r = 2 atan2(n, w)/n`
  is carried as two terms and `n` with its own error, where `r · u` rounded `r` just above `2` and
  then each product. `so3_log`'s maximum falls from 2.317 to 1.680 `u`, `theta:dense` to 1.578 and
  `theta:pi-1e-k` from 1.24–1.79 to 0.68–0.88; `sen3_log`, `sim3_log`, the SE(3) charts' `local`
  and `gaussian_mahalanobis_se3`/`_se23` move with it. Every caller's bits move on that arm.
  The compensated arm is a function of its own, so `log` stays small enough to inline: inline, it
  pushed `log` out of line at callers with several call sites, where reading back a just-written
  quaternion stalls and a composed `Log` took 25 ns near the identity instead of 10.

### Added

- **`helicoid_linalg::math::{atan2, atan2f}`** (`0105`): the kernels behind `Real::atan2`, public,
  with `libm`'s names and argument order. Their polynomials and tables are generated by
  `just atan-fit` and drift-checked in CI (`atan-fit-check`).

- **`helicoid_linalg::eig3_values`** (`0104`): `eig3`'s eigenvalues without its frame, bit for bit
  at every `S: Real` (`eig3_values_matches_reference`; the composition `eig3(a).0` is its twin).
  It takes 0.46–0.83× `eig3`'s time, the least gain on a random spectrum and the most at a triple
  eigenvalue. The two share one sort network, and `eig3`, `eig3_values` and `svd3` read one
  scale-back factor that is `NaN` for a non-finite input, one product per value where a product
  and a select were; every latency is within noise of its time before. No existing output moves.

### Removed

- **The `tf_tree_math` oracle** (`0107`): `tf_tree` is migrating onto `helicoid`, so its math is a
  consumer, not an independent oracle. `runners/tf_tree_math`, `just oracle-tf-tree-math` and the
  envelope's `libm`-crate class of domination failures are gone; the one-variable twin attributes
  every failure. No stratum loses its only oracle (2797 paired before and after). Domination
  loosens on the 42 strata `tf_tree_math` was best on; no-regress, once a baseline is blessed,
  holds them exactly. Library crates are unchanged.

## Released versions

Frozen history, newest first, in [`docs/changelog/`](./docs/changelog/); leave it
out of searches and sweeps.

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
