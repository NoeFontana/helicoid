# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Released
entries live in [`docs/changelog/`](./docs/changelog/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [0.0.12] - 2026-11-18

The fused preintegration step release. `so3::exp_gamma_apply` returns `Exp φ`, `Γ₁(φ) v` and `Γ₂(φ) v`, and `so3::exp_gamma_apply_jacobians` adds both directional Jacobians, from one coefficient evaluation: one root, one `sin_cos` and at most one `sin` where the separate calls take three, one and four. Every output is the separate calls' to the bit, `Dual` included, at 5–43% fewer instructions per call (`0119`). **No API break**, and every existing output keeps its bits.

### Added

- **`so3::exp_gamma_apply` and `so3::exp_gamma_apply_jacobians`**
  ([`0119`](docs/decisions/0119-exp-and-gamma-share-one-evaluation.md)): a preintegration step's
  `Exp φ`, `Γ₁ v` and `Γ₂ v`, and both directional Jacobians, from one coefficient evaluation. Every
  output is `SO3::exp`, `gamma1(φ) * v`, `gamma2(φ) * v` or `gamma_apply_jacobians::<1 | 2>` to the
  bit, `Dual` included, so moving from the separate calls changes no result. They take 10–43%
  (binary64) and 15–37% (binary32) fewer instructions for the values, and 17–36% and 5–25% fewer
  with the Jacobians. The alias `so3::GammaJacobians<S>` names `gamma_apply_jacobians`'
  `(Vec3, Mat3, Mat3)`. Every existing output keeps its bits.

## Released versions

Frozen history, newest first, in [`docs/changelog/`](./docs/changelog/); leave it
out of searches and sweeps.

- [`0.0.11`](./docs/changelog/0.0.11.md) — 2026-11-15 (the triangular factorizations `rq3` and `lq3`)
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
[0.0.12]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.12
