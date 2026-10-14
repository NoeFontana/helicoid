# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Released
entries live in [`docs/changelog/`](./docs/changelog/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [Unreleased]

## [0.0.4] - 2026-10-14

A patch release for locus-sensors' requests: `SO3::from_matrix` keeps exact rotations on the unit sphere and matches `tf_tree_math`'s accuracy at `f64` (`0086`), `Real::round_from_f64` rounds a runtime `f64` (`0087`), and `just determinism` checks D16 across x86_64, aarch64 and wasm32. No API break: `round_from_f64` is a provided method.

### Added

- **`just determinism`** (`PHASE6.md` §3, D16): `cargo xtask determinism` digests the shipped subject's outputs over the whole corpus, one SHA-256 per id and precision, NaN canonicalized. CI digests x86_64, aarch64 and `wasm32-wasip1` and fails on any difference. At `0.0.3` all 114 digests agree on x86_64 and `wasm32-wasip1`.

- **`Real::round_from_f64(x: f64) -> Self`** (`0087`): a runtime `f64` correctly rounded to nearest-even at the scalar's precision, the identity at `f64`, `as f32` at `f32` (IEEE's conversion, the same bits on every target), a constant for `Dual`. Every `f64` is in the domain: at `f32` it overflows to `±inf` from `(2 − 2⁻²⁴)·2¹²⁷` and reaches `±0` at `2⁻¹⁵⁰`. A provided method, so an `impl Real` outside the workspace keeps compiling. Corpus id `real_round_from_f64`, with binary32 ties, near-ties and both edges: every stratum reads only the reference's 30-digit rounding.

### Changed

- **`SO3::from_matrix` returns Shepperd's quotients** (`0086`), in `tf_tree_math`'s arithmetic, whenever they land within `16u` of the unit sphere, which every rotation rounded to the precision does. A scaled, non-orthogonal or NaN matrix still gets 0.0.3's normalized column. Exact rotations, such as quarter turns and frame conventions' signed permutations, come back with `fl(‖q‖²) = 1` at `f64` (0.0.3: `1 − 2u`, which every later `act` carried). Under `0085`'s score, strata above the best oracle fall from 18 to 8 of 31 at `f64` and from 14 to 9 at `f32`, and `f64` now equals `tf_tree_math` on every stratum's maximum but `q:nonunit` (89.4 against 163.0 `u`). Six `f64` strata near `π` read up to 0.50 `u` higher than 0.0.3.

## Released versions

Frozen history, newest first, in [`docs/changelog/`](./docs/changelog/); leave it
out of searches and sweeps.

- [`0.0.3`](./docs/changelog/0.0.3.md) — 2026-10-13 (Phase 5, `eig3` as Jacobi, three C++ oracles)
- [`0.0.2`](./docs/changelog/0.0.2.md) — 2026-10-09 (the carried path and a faster geodesic)
- [`0.0.1`](./docs/changelog/0.0.1.md) — 2026-10-08 (first release)

[0.0.1]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.1
[0.0.2]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.2
[0.0.3]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.3
[0.0.4]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.4
