# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Released
entries live in [`docs/changelog/`](./docs/changelog/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [0.0.10] - 2026-11-13

The continuous-time spline release (`helicoid-spline`). Adds uniform cumulative B-splines of order `K = 2..=6` on any `LieGroup` with `O(K)` control Jacobians and compile-time exact integer tables (`0113`–`0115`), `Jac::add` and `Jac::scale` on `helicoid` (`0112`), and the 16 spline corpus ids with `basalt-headers` oracle (`0116`). **Breaking:** external `Jac` impls must add `add` and `scale`; every existing output keeps its bits.

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
- **`just oracle-basalt`**: basalt-headers (Sommer and Usenko's implementation of the same
  recursions) as a pinned container oracle for the spline ids, binary64 only. On the committed tier
  `helicoid` reads 3.9–4.2 `u` on `so3_spline_value_jac` against its 249–274 `u`, and 0 `u` on the
  Gram matrix against 752 `u`.
- **`bench-gate` reads criterion's samples where criterion writes them**: `$CRITERION_HOME`, else
  `$CARGO_TARGET_DIR/criterion`, else cargo's target directory, where it read
  `<root>/target/criterion` and found no samples under a moved target directory. A relative
  variable is refused, since the bench binary and cargo would resolve it from different
  directories. Found in a consumer that copied the gate.
- **The bench tests write to per-process scratch directories** removed on drop, where they shared
  fixed `/tmp/helicoid-bench-*` paths, so concurrent test runs in two worktrees deleted each
  other's fixtures; so does the `untracked` lint's git test (`/tmp/helicoid-lint-untracked`).
  Nothing breaks.

## Released versions

Frozen history, newest first, in [`docs/changelog/`](./docs/changelog/); leave it
out of searches and sweeps.

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
