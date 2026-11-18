# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Released
entries live in [`docs/changelog/`](./docs/changelog/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [Unreleased]

### Added

- **`so3::exp_gamma_apply` and `so3::exp_gamma_apply_jacobians`**
  ([`0119`](docs/decisions/0119-exp-and-gamma-share-one-evaluation.md)): a preintegration step's
  `Exp φ`, `Γ₁ v` and `Γ₂ v`, and both directional Jacobians, from one coefficient evaluation. Every
  output is `SO3::exp`, `gamma1(φ) * v`, `gamma2(φ) * v` or `gamma_apply_jacobians::<1 | 2>` to the
  bit, `Dual` included, so moving from the separate calls changes no result. They take 9–44%
  (binary64) and 15–37% (binary32) fewer instructions for the values, and 17–36% and 5–25% fewer
  with the Jacobians. The alias `so3::GammaJacobians<S>` names `gamma_apply_jacobians`'
  `(Vec3, Mat3, Mat3)`. Every existing output keeps its bits.

## [0.0.11] - 2026-11-15

The triangular factorization release. `helicoid_linalg::rq3` and `lq3` factor a 3×3 matrix as a triangle times a rotation, `M = U R` and `M = L Q`, with `UpperTri3` and `LowerTri3`: classical Gram–Schmidt twice under Kahan's test on rows each scaled by its own power of two, a triangle with a positive diagonal returned with the identity bit for bit, and at most 4.33 `u` (binary64) and 4.00 `u` (binary32) on the new corpus id `rq3`, below nalgebra's and Eigen's Householder QR on every stratum (`0117`). **No API break**, and every existing output keeps its bits.

### Added

- **`rq3` and `lq3`**, the 3×3 RQ and LQ factorizations, with the triangle types `UpperTri3` and
  `LowerTri3` ([`0117`](docs/decisions/0117-rq-and-lq-are-gram-schmidt-twice-on-scaled-rows.md)).
  `M = U R` and `M = L Q` with a rotation and a triangle whose diagonal is positive where
  `det M > 0`, reported by a mask; under a clear mask the factorization is kept. Classical
  Gram–Schmidt twice under Kahan's test on rows each scaled by its own power of two: backward stable
  row by row and orthogonal to a few `u` however close two rows are to parallel; `lq3` is `rq3`
  reversed, bit for bit. A triangle with a positive diagonal returns itself and the identity
  rotation bit for bit, at both precisions and in `Dual` values. On the corpus it
  reads at most 4.33 `u` (binary64) and 4.00 `u` (binary32): below nalgebra's and Eigen's
  Householder QR on every stratum, by at least 1.5× wherever it reads above rounding. Every existing
  output keeps its bits.

### Instrument

Not shipped in the crates.
- **The corpus id `rq3`** ([`0117`](docs/decisions/0117-rq-and-lq-are-gram-schmidt-twice-on-scaled-rows.md)
  decision 9): 24 binary64 and 23 binary32 strata against RQ.1's construction at 120 digits, the
  metric `Rule::Rq` (rowwise backward error, orthogonality, forward errors over the reference's κ),
  two seeded defects, and nalgebra's `QR` and Eigen's `HouseholderQR` as oracles.

## Released versions

Frozen history, newest first, in [`docs/changelog/`](./docs/changelog/); leave it
out of searches and sweeps.

- [`0.0.10`](./docs/changelog/0.0.10.md) — 2026-11-13 (the continuous-time spline, `helicoid-spline`)
- [`0.0.9`](./docs/changelog/0.0.9.md) — 2026-10-29 (general symmetric eigensolver `eig_sym`)
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
[0.0.10]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.10
[0.0.11]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.11
