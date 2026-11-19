# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Released
entries live in [`docs/changelog/`](./docs/changelog/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [Unreleased]

## [0.0.13] - 2026-11-19

The `eig_sym` to 64 release. `helicoid_linalg::eig_sym` accepts `4 <= N <= 64`, from `4 <= N <= 16`: the same Jacobi with its sweep cap measured at every `N` (`0118`). It is also 3 to 16% faster from `N = 9`. **No API break**, and every existing output keeps its bits.

### Changed

- **`eig_sym` accepts `4 <= N <= 64`**, from `4 <= N <= 16`
  ([`0118`](docs/decisions/0118-eig-sym-reaches-n-64.md)). It is the same Jacobi with its sweep cap
  measured at every `N` to 64. At `N <= 16` every output keeps its bits.
- **`eig_sym` is 3 to 16% faster from `N = 9`**, every output's bits unchanged: row `p` of the
  working matrix is mirrored once per row of pairs, not after each rotation
  ([`0118`](docs/decisions/0118-eig-sym-reaches-n-64.md) decision 9).

### Instrument

Not shipped in the crates.
- **The corpus ids `eig_sym_n{24,40,64}`**
  ([`0118`](docs/decisions/0118-eig-sym-reaches-n-64.md)): `0111`'s 45 strata at 2, 1 and 1 records
  each (64, 32 and 16 in the deep tier), against nalgebra and Eigen. The committed corpus cap is
  100 MiB, up from 80.

## Released versions

Frozen history, newest first, in [`docs/changelog/`](./docs/changelog/); leave it
out of searches and sweeps.

- [`0.0.12`](./docs/changelog/0.0.12.md) — 2026-11-18 (the fused preintegration step `exp_gamma_apply`)
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
[0.0.13]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.13
