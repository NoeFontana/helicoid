# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Released
entries live in [`docs/changelog/`](./docs/changelog/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [Unreleased]

### Added

- **`helicoid_linalg::nearest_rotation`** (`0102`): the Frobenius-nearest rotation of a 3×3,
  `U diag(1, 1, det UVᵀ) Vᵀ`, one of the nearest where it is not unique. It runs Newton's polar
  iteration with Frobenius scaling where `det M > τ‖M‖³`. Where `det M < −τ‖M‖³` it takes the
  orthogonal polar factor and a projector reflection, and between them `svd3`'s `UVᵀ`. Its `Dual`
  tangent is the derivative of the projection wherever that is unique off the `svd3` arm, ties
  among the singular values included. It reads at most 5.08 `u` over the corpus at both
  precisions, and dominates nalgebra's `from_matrix_eps` and Sophus's `fitToSO3` on every stratum
  at both precisions, by at least 4.24× (`f64`) and 3.24× (`f32`). It runs in 99.5 ns on a
  rotation `1e-8` off SO(3) at `f64` and 52.6 at `f32`, 286–603 ns and 227–486 elsewhere.

### Instrument

Not shipped in the crates.
- **The corpus id `nearest_rotation`** (`0102`): 17 binary64 strata (296 records) and 15 binary32
  (264), against `mp.svd_r` at 120 digits. The metric `Rule::Polar` weights the forward error by
  the projection's conditioning, the weight exact in integers. A self-test pair goes with it, and
  nalgebra's `Rotation3::from_matrix_eps` (its input at unit scale) and Sophus's `SO3::fitToSO3`
  are oracles.
- **Test profiles.** `just test`'s second pass runs at the new `ndebug` profile, the dev profile with
  `debug_assert!` and overflow checks off, instead of `--release`. The dev profile optimizes the
  hot third-party crates and `xtask`, and keeps line tables only. After an edit to
  `helicoid-linalg`, the two passes take 125 s instead of 400 s (392 CPU-s instead of 1252), and
  the dev target directory is 30% smaller. `just determinism` and CI also digest the corpus at
  `ndebug`, so an opt-level that moved an output bit fails the comparison.

## [0.0.6] - 2026-10-19

The 3×3 SVD release. `helicoid_linalg::svd3` uses the signed convention, so `UVᵀ` is the nearest rotation where that is unique. It is backward stable at every finite input and dominates nalgebra's `SVD` and Eigen's `JacobiSVD` on every stratum at both precisions. `eig3`'s scaling ladder gets faster with bit-identical outputs. **No API break.**

### Added

- **`helicoid_linalg::svd3`** (`0101`): the 3×3 SVD in the signed convention, `U, V ∈ SO(3)`,
  `σ₁ ≥ σ₂ ≥ |σ₃|`, `sign σ₃ = sign det A`, so `UVᵀ` is the nearest rotation where that is unique.
  One-sided Jacobi on `eig3`'s rotation, with `eig3`'s exit and polish, then Gram–Schmidt with the
  projection taken twice. It is backward stable column by column at every finite input, subnormal
  ones included: 6.03 `u` at most over the corpus at `f64`, 5.08 at `f32`. It dominates nalgebra's
  `SVD` (McAdams, `4.2×10⁷ u` at `graded-1e-8`) and Eigen's `JacobiSVD` on every stratum at both
  precisions. It runs in 265–412 ns at `f64` and 167–364 ns at `f32`.

### Changed

- **`eig3`'s ladder** (`0101`) passes only scalars and takes no product on its chain. Its `f32`
  arm runs 0.81× as long; every output is bit-identical.

### Instrument

Not shipped in the crates.
- **The corpus id `svd3`** (`0101`): 34 binary64 strata (580 records) and 27 binary32 strata
  (468), against `mp.svd_r` at 120 digits. The metric `Rule::Svd` and a self-test pair go with it,
  and nalgebra's `SVD` and Eigen's `JacobiSVD` (in the Sophus image) are oracles. A record part may
  key its `shape` by field.

## Released versions

Frozen history, newest first, in [`docs/changelog/`](./docs/changelog/); leave it
out of searches and sweeps.

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
