# 0072: a container oracle is an image on the same protocol

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #142

## Context

`PHASE1.md` §7 owes three container runners, `docker/oracles/{sophus,manif,gtsam}`, and §0.0
lists them pending. The three runners that exist are workspace-excluded crates started by
`cargo run --release --locked` (`xtask/src/conformance/oracle.rs`, `Runner`), and none of their
oracles has Sim(3): `PHASE5.md` §7 leaves "§3 dominated against oracles" open for want of one
([`0070`](./0070-sim3-a-structured-jacobian-and-a-series-on-delta.md) *Consequences*). Sophus C++
has `Sim3`, header-only, so it is the first container.

What the file protocol needs from a launcher is a command prefix to which `--version` or
`--out DIR FILE…` is appended, the paths being the host's. Nothing in §7 says how an image is
built, which paths it sees, who owns what it writes, or what a runner does with an input its
library refuses.

Sophus 1.24.6 (`d0b7315`, its latest release) exports, for ids the corpus has: `SO2`/`SE2`
`exp`, `log`, `Adj`; `SO3` `exp`, `log`, the action, the quaternion of a matrix, `leftJacobian`,
`leftJacobianInverse`; `SE3` the same Jacobians, `exp`, `log`, `Adj`; `interpolate`, which is
`X₀ Exp(t Log(X₀⁻¹ X₁))` (`NUMERICS.md` §10); `Sim3` `exp`, `log`, `Adj` and the action. It has
no Jacobian of `SO2`, `SE2` or `Sim3`. A failed `SOPHUS_ENSURE` aborts unless
`SOPHUS_ENABLE_ENSURE_HANDLER` names a handler; `RxSO3`'s constructors ensure a scale in
`[10⁻¹⁰, 10¹⁰)`, and `RxSO3::exp` clamps to that range without one.

## Decision

1. **The launcher.** `oracle.rs`'s `Runner` gains `launch: Launch`, `Crate(dir)` (unchanged) or
   `Container(dir)`. A container runner is built by `docker build --quiet --tag
   helicoid-oracle-<name> docker/oracles/<name>` and run by the image id that build printed, so
   a concurrent rebuild cannot change what a run executes. Every launch is `docker run --rm
   --network=none --read-only --user=<uid>:<gid>` with two bind mounts, each at its host path so
   the protocol's arguments pass unchanged: the corpus directory read-only, and
   `conformance/results/oracle` (created first, so docker does not create it as root) writable.
   `uid:gid` is that directory's owner. A path that is relative or holds a comma is an error
   before docker runs. The file protocol, `FileSubject` and every check of `execute` are
   unchanged.
2. **The image** (`docker/oracles/sophus/Dockerfile`): `ubuntu:24.04` by digest; `g++`,
   `libeigen3-dev`, `nlohmann-json3-dev` from the Ubuntu archive snapshot `20261001T000000Z`
   (gcc 13.3.0, Eigen 3.4.0, nlohmann-json 3.11.3, glibc 2.39); Sophus by commit and the SHA-256 of
   its archive. Compiled `-std=c++17 -O2 -ffp-contract=off -DEIGEN_DONT_VECTORIZE` with warnings
   denied, so the oracle is the same scalar program on every host. The build compiles and runs
   `convert_test`, then the runner. `--version` prints
   `sophus@<version>+<commit> eigen@… gcc@… glibc@…`.
3. **The ids** (`Runner::answers`, 26): `so2_{exp,log}`, `se2_{exp,log,ad}`,
   `so3_{exp,log,act,from_matrix,jl,jr,jl_inv,jr_inv,geodesic}`, `se3_geodesic`,
   `sen3_{exp,log,ad,jl,jr,jl_inv,jr_inv}_n1`, `sim3_{exp,log,ad,act}`. The right Jacobians and
   their inverses are the left ones at `−τ` (`NUMERICS.md` §1), as sophus-rs's runner answers them.
   Backend `HostStd`: C++ on the image's glibc.
4. **The conversions** (`convert.hpp`), named after the foreign convention, each tested by hand
   in `convert_test.cpp`:
   - quaternion: helicoid `[w, x, y, z]` is `Eigen::Quaterniond(w, x, y, z)`, which stores
     `[x, y, z, w]`;
   - tangent: helicoid's component `i` sits at Sophus's `at[i]`, with `at` = `[2, 0, 1]` for
     SE(2) (`[υ; θ]`), `[3, 4, 5, 0, 1, 2]` for SE(3) (`[υ; ω]`) and `[3, 4, 5, 0, 1, 2, 6]` for
     Sim(3) (`[υ; ω; σ]`); a Jacobian or adjoint's helicoid `(i, j)` is Sophus's
     `(at[i], at[j])`, written column-major;
   - elements: `SO2(z₀, z₁)`, `SO3(q)` (both normalize), `SE3(SO3(q), x)`, and
     `Sim3(RxSO3(e^σ, SO3(q)), x)`; back, `q` is `rxso3().so3()`'s quaternion and
     `σ = log(scale())`, `scale()` being the stored `√s · q`'s squared norm, as `RxSO3::log`
     reads it.
5. **A refused record is answered `nan`.** `refused.cpp` makes a failed `SOPHUS_ENSURE` throw;
   the runner answers that record with `nan` in every output field and counts it on stderr. It is
   answered, so no stratum is the runner's choice (§7), and non-finite, so its row is unscored
   and never the best (`xtask/src/envelope/bars.rs`). A clamped answer is not a refusal: it is
   scored.
6. **Binary64 only**, as every runner (`execute`'s `no_f32`): the protocol has no precision, and
   `Sophus::Sim3<float>` rows are left to *Further work*.

## Rationale

- **Same-path mounts** keep the protocol the one `execute` already checks; translating paths
  would put a second argument grammar between the harness and the runner. `--mount` over `-v`
  because `-v` creates a missing source as root and splits on `:`.
- **The image id, not the tag**, because the tag is mutable and the build is cached anyway.
- **`--network=none`, `--read-only`, `--user`**: everything the runner computes with is fixed at
  build time, it writes only the answers' directory, and those files are the user's to delete.
- **The archive snapshot** pins the compiler, Eigen and glibc together; pinning three package
  versions would leave gcc's and glibc's point releases to the day of the build.
- **`nan` over failing the run** (sophus-rs's runner fails, but no corpus record reaches its
  refusal): 216 Sim(3) records are outside `RxSO3`'s domain, and failing would leave Sophus
  unusable on every id. **`nan` over disabling the ensures**: the result would be a guess the
  library itself refuses.
- **No contraction, no SIMD**: whether gcc fuses `a·b + c` depends on the dialect and the target
  (aarch64 has the instruction in its base ISA), and Eigen's packet sums reassociate, so without
  the flags one pin could name two programs.

## Consequences

Measured on `a6d3ced` (`perf/sim3-bench`) with the four oracles' rows, `just envelope` (no
baseline in the tree, so no no-regress verdict):

- **Pairing.** 434 binary64 strata newly have an oracle, all Sim(3): `sim3_exp` 122, `sim3_log`,
  `sim3_ad`, `sim3_act` 104 each. The other 18 strata of each of the last three
  (`|σ| ∈ {30, 100, 700}`) are refused (`RxSO3`'s scale), so unpaired. `sim3_exp` there is
  clamped and scored at `2·10¹⁵`–`9·10¹⁵ u`. `so2_*`/`se2_*` rows wait for SO(2)/SE(2) (`PHASE3.md`
  §0.0).
- **Domination failures**: 148 → 301, the 153 new ones being `sim3_*` 87, `so3_*` 34 and
  `sen3_*_n1` 32; 9 of them close under `helicoid:host-std` (D16's cost). Sophus becomes the best
  oracle of 35 failures that already failed. Two are not near ties: `so3_act/q:nonunit`
  (510.3 against 7.867 `u`) and `sen3_ad_n1/q:nonunit` (256.1 against 2.132), where Sophus
  normalizes the drifted `q` at construction and `helicoid` carries it
  ([`0058`](./0058-a-drifted-quaternion-is-carried-not-vouched-for.md)). The other 151 are
  within 3.17× on maxima of at most 10.3 `u`.
- **Sim(3)** fails 87 of its 434 paired strata (`sim3_exp` 15, `sim3_log` 17, `sim3_ad` 33,
  `sim3_act` 22), by 1.00×–3.17×, the worst `sim3_ad/sigma:-3/theta:0.3` (1.350 against
  0.4263 `u`). `PHASE5.md` §7's Sim(3) line stays open. That stratum's loss is `Ad`'s `R` block
  (1.0–1.3 `u` of its 1.35, by replaying its four records), not `e^σ` (`helicoid:host-std` reads
  the same 1.35): `Quat::to_matrix` forms the diagonal as `(w² - |u|²) + 2x²`, the homogeneous
  `R(q)` of `NUMERICS.md` §1 that the reference also evaluates, and carries the rounding of
  `w² ≈ 1` near the identity; Eigen's `1 - 2(y² + z²)` does not, because it assumes `‖q‖ = 1`.
- `just oracle-sophus` needs docker; the first build fetches from the snapshot and GitHub, every
  later run is offline and takes seconds.

## Implementation plan

1. This record, the launcher and its tests, `docker/oracles/sophus`, `just oracle-sophus`, the
   status rows. Verified by `cargo nextest run -p xtask -- oracle`, the image build (it runs
   `convert_test`), `just oracle-sophus`, `just envelope`, `just lint`, `just test`, `just doc`.

## Open questions

None.

## Further work

- manif and GTSAM on the same `Launch::Container`.
- `f32` oracle rows: a protocol precision and `Sophus::Sim3<float>` would let Sim(3) be judged at
  both precisions.
- Whether the two `q:nonunit` failures are an exception under `0058` (`0046`'s registry).
- Whether `R(q)` should drop `w²` where `q` is held unit, which would win `sim3_ad`'s and
  `so3_act`'s near-identity strata at the cost of erring by `‖q‖² - 1` on a drifted `q`: a
  `NUMERICS.md` §1 change, so a record of its own.
- The nightly `oracles` CI job (`PHASE1.md` §3); the `envelope` job that would run
  `just oracle-sophus` is disabled.
