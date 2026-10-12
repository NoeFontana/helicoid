# 0074: the GTSAM oracle is built from source under our flags

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #143

## Context

`PHASE1.md` §7 owes `docker/oracles/gtsam` ("`Rot3`/`Pose3` `Expmap`/`Logmap` and their
derivatives"); [`0072`](./0072-a-container-oracle-is-an-image-on-the-same-protocol.md) built the
launcher (`Launch::Container`) and the first image, Sophus. Sophus answers four Sim(3) ids and fails
`helicoid` on 87 of the 434 strata it pairs; nothing answers `sen3_*_n2`, `sen3_*_n3`, `so3_gamma2`
or the SE(3) charts.

GTSAM 4.3.0 (`71a25ca`, 2026-09-18, its latest release) exports, for ids the corpus has, read from
its source at that commit:

- `Rot2` `Expmap`/`Logmap`; `Pose2` `Expmap`, `Logmap`, `AdjointMap`, `ExpmapDerivative`,
  `LogmapDerivative` of a tangent, tangent `[v; ω]`.
- `Rot3` `Expmap`, `Logmap` (Shepperd on the matrix, then `atan`), `rotate`, `toQuaternion`,
  `ExpmapDerivative`, `LogmapDerivative` of a tangent, `slerp`; `so3::DexpFunctor::Gamma()`, the
  kernel `½I ± CΩ + EΩ²`.
- `ExtendedPose3<K>`, SE_K(3) with tangent `[ω; ρ₁; …; ρ_K]` and element `(R, x₁, …, x_K)`, the
  columns of a `3 × K` block: `Expmap`, `Logmap`, `AdjointMap`, `ExpmapDerivative`,
  `LogmapDerivative` of a tangent. `Pose3` is `ExtendedPose3<1, Pose3>` with its own `Expmap`;
  `NavState` is `ExtendedPose3<2, NavState>` with columns `(t, v)`; `Se23` is `ExtendedPose3<2>`.
  `Pose3::retract`/`localCoordinates` go through `ChartAtOrigin`, which is `Expmap`/`Logmap` under
  the default `GTSAM_POSE3_EXPMAP`; `slerp` is `interpolate`, `X Expmap(t Logmap(X⁻¹Y))`.
- `Similarity3`, stored `(R, t, s)` as `[R t; 0 s⁻¹]`, acting as `s(Rp + t)`, tangent `[ω; ρ; λ]`
  with `s = e^λ`: `Expmap`, `Logmap`, `AdjointMap`, `transformFrom`. Neither map takes a Jacobian
  argument and no other function returns one, so `sim3_j*` has no GTSAM row.
- `ExpmapDerivative` is the right Jacobian (`Expmap(ξ + δ) ≈ Expmap(ξ) Expmap(J δ)`) and
  `LogmapDerivative` of a tangent its inverse; no left one is exported but `DexpFunctor`'s SO(3)
  pair.
- `Unit3`'s chart has its own basis, not `NUMERICS.md` §8's, so `s2_*` is another function; the
  decoupled SE(3) chart is a compile-time alternative to the default (`GTSAM_POSE3_EXPMAP=OFF`).

A pinned build can be GTSAM's C++ source or the `gtsam` wheel. Both run GTSAM's C++, but the wheel
is a binary compiled with flags we do not choose: whether it contracts `a·b + c` and whether Eigen
vectorizes is per wheel, so one pin would name a program per architecture (0072's *Rationale*).

## Decision

1. **The image** (`docker/oracles/gtsam/Dockerfile`), on 0072's launcher unchanged: `ubuntu:24.04`
   by 0072's digest; `g++`, `cmake`, `make`, `nlohmann-json3-dev` from the same archive snapshot
   `20261001T000000Z` (gcc 13.3.0, cmake 3.28.3, nlohmann-json 3.11.3, glibc 2.39); GTSAM 4.3.0 by
   commit `71a25ca36c084cbad1f872e812d6d97fbadfdb05` and the SHA-256 of its archive, with the Eigen
   it bundles (3.4.0). GTSAM is **built in the image**, a static `Release` library, and the runner
   is built against its installed CMake package, both under `-ffp-contract=off
   -DEIGEN_DONT_VECTORIZE` and with `GTSAM_BUILD_WITH_MARCH_NATIVE=OFF`. Every option that reaches
   the numerics keeps GTSAM's default (`Rot3` a matrix, `GTSAM_POSE3_EXPMAP`,
   `GTSAM_ROT3_EXPMAP`); the others drop what the runner never calls (Boost, TBB, METIS, tests,
   examples, `gtsam_unstable`) and GTSAM's own `-Werror`, which gcc 13 trips in `Gal3.cpp` without
   Eigen SIMD. The runner and `convert_test` are compiled `-Wall -Wextra -Wpedantic -Werror`; the
   build runs `convert_test` and stops on a failure. `--version` prints
   `gtsam@4.3.0+<commit> eigen@3.4.0 gcc@… glibc@…`.
2. **The ids** (`Runner::answers`, 47): `so2_{exp,log}`,
   `se2_{exp,log,ad,jl,jr,jl_inv,jr_inv}`,
   `so3_{exp,log,act,from_matrix,jl,jr,jl_inv,jr_inv,gamma2,geodesic}`, `se3_geodesic`,
   `se3_screw_{retract,local}`, `sen3_{exp,log,ad,jl,jr,jl_inv,jr_inv}_n{1,2,3}` (`Pose3` at
   `N = 1`, `ExtendedPose3<N>` above), `sim3_{exp,log,ad,act}`. Every right Jacobian is
   `ExpmapDerivative`, every right inverse `LogmapDerivative` of the tangent, and the left ones the
   same at `−τ` (`NUMERICS.md` §1), for SO(3), SE(2) and SE_N(3) alike. `so3_from_matrix` is
   `Rot3(R).toQuaternion()`, Eigen's quaternion of the matrix as given. Backend `HostStd`.
3. **The conversions** (`convert.hpp`), named after GTSAM, each tested by hand in
   `convert_test.cpp`:
   - a quaternion is `Rot3::Quaternion(w, x, y, z)`, the matrix Eigen's `toRotationMatrix` forms,
     **not normalized**; back, `toQuaternion()`'s `(w, x, y, z)`;
   - SO(3), SE_N(3) and Sim(3) tangents are in helicoid's order; SE(2)'s component `i` sits at
     GTSAM's `[2, 0, 1][i]`, a matrix's `(i, j)` at `(at[i], at[j])`; matrices are column-major;
   - SE_N(3)'s `[x₁; …; x_N]` are `ExtendedPose3<N>`'s columns in order;
   - SO(2) is `Rot2::fromCosSin(c, s)` (which normalizes past `10⁻¹⁰`), SE(2) `Pose2(Rot2, t)`;
   - Sim(3) `(q, x, σ)` is `Similarity3(Rot3(q), x / e^σ, e^σ)`; back, `x = s · t`, `σ = log s`.
   The derivatives' sides are tested too: `Rot3::ExpmapDerivative` at `φ = (0, 0, π/2)` against
   `J_r`'s `2/π`, `Pose3`'s and `ExtendedPose3<2>`'s at `φ = 0` against `[I 0; −½ρ^ I]`.
4. **No refusal path.** GTSAM refuses no corpus record, so an exception fails the run. A
   non-finite answer is GTSAM's (`sim3_log` at `σ = −700`, where `V⁻¹`'s cofactors overflow) and is
   recorded, unscored.
5. **Binary64 only**, as every runner.

## Rationale

- **Source over the wheel**: the flags are the oracle. A wheel pins bytes, not a program we can
  state; built here, the library and the runner are the scalar program 0072 requires of Sophus.
  The cost is a first build of minutes (the library's 237 translation units); every later run is
  offline and cached.
- **`ExtendedPose3<N>` over `NavState`** at `N = 2`: the same code (`NavState` adds no Lie
  operation), and columns in order need no names. `NavState`'s `(t, v)` and helicoid's `SE23`
  `(v, p)` disagree on the names; SE_N(3) does not depend on them, so either column order answers
  the same function, but only the generic one also serves `N = 3`.
- **`Pose3` at `N = 1`**: it is the class GTSAM's users call, and its `Expmap` is its own.
- **Unnormalized quaternions**: `Rot3::Quaternion` is GTSAM's entry point and does not normalize;
  normalizing in the runner would answer for a step GTSAM does not take. GTSAM's `q:nonunit` rows
  carry the drift (`so3_act` 510.8 `u`) and win nothing there.
- **Right Jacobians at `−τ`** rather than `DexpFunctor`'s SO(3) left pair: one rule for every
  group, and the negation is exact, so the two are the same numbers.

## Consequences

Measured on this branch (`d4856a5` + this record) with all five oracles' rows, `just envelope`
(no baseline in the tree):

- **Pairing.** 2001 binary64 strata have an oracle, 919 more than with four: `sen3_*_n2` and
  `sen3_*_n3` 736 (their first oracle), `se3_screw_{retract,local}` 104, `so3_gamma2` 28, and the
  51 Sim(3) strata Sophus refuses (`|σ| ∈ {30, 100, 700}`; `sim3_log` at `σ = −700` is GTSAM's
  `nan`). Sim(3) is paired on 485 strata.
- **Domination failures**: 301 → 634. The 333 new ones are `sen3_*_n2/n3` 169, `sen3_*_n1` 31,
  `se3_screw_*` 43, `sim3_*` 64, `so3_*` 26 (`so3_act` 12, `so3_gamma2` 9); 10 of them close under
  `helicoid:host-std` (D16's cost), the rest do not. GTSAM becomes the best oracle of 49 failures
  that already failed. None is gone.
- **Ratios.** Every new failure is within 2.5× except `sen3_ad_n2/theta:1e-2` (1.565 against 0.4450
  `u`, 3.52×) and `sen3_ad_n3` there (2.77×); the largest new maximum is 7.80 `u`
  (`sen3_jr_n2`). Replaying `sen3_ad_n2`'s argmax: helicoid's `R` block errs by 1.37–1.74 `u` on the
  diagonal and Eigen's by at most 0.49, 0072's diagnosis again — `Quat::to_matrix`'s homogeneous
  `(w² − |u|²) + 2x²` carries the rounding of `w² ≈ 1`.
- **Sim(3)** fails 151 of 485 paired strata (`sim3_exp` 30, `sim3_log` 27, `sim3_ad` 48,
  `sim3_act` 46); GTSAM beats `helicoid` on 23, 17, 38 and 38 of them. Three are past 3×, all
  already failing against Sophus: `sim3_act/sigma:1/theta:pi-1e-6` (6.660 against GTSAM's 1.221
  `u`), `sim3_ad/sigma:-3/theta:0.3` (1.350 against 0.3730, 0072's `R` block) and
  `sim3_exp/sigma:1/theta:3` (5.815 against 1.924). The two at `σ = 1` are half D16:
  `helicoid:host-std` reads 3.549 and 3.254. Replaying `sim3_act`'s argmax (record 277): the
  `libm` crate's `exp(1)` is one ulp above the correctly rounded `e` (6.660 → 3.549 with glibc's),
  and the rest is the quaternion sandwich's 1.41 `u` on `Rp`, which `s R p + t` amplifies (`|sRp|
  = e` against `|Xp| = 1.56`); with `Rp` rounded exactly the same tail gives 1.109. `PHASE5.md`
  §7's Sim(3) line stays open.
- GTSAM's own losses, for the record: SE(2)'s Jacobians cancel below `θ = 10⁻¹` (`se2_jl_inv`
  8.2·10¹³ `u` at `θ = 10⁻⁵`), `Log` through the matrix loses `theta:subnormal` and the quaternion
  sign at `q:w0` (`sen3_log` 1.8·10¹⁶ `u`, the other branch at exactly π), `Similarity3::Expmap` is
  1.4·10⁵ `u` at `σ = 700`. Oracles may be wrong (§7); the best one sets the bar.
- `just oracle-gtsam` needs docker; the first build fetches from the snapshot and GitHub and
  compiles GTSAM.

## Implementation plan

1. This record, `docker/oracles/gtsam`, the `RUNNERS` entry, `just oracle-gtsam`, the status rows.
   Verified by the image build (it runs `convert_test`), `just oracle-gtsam`, `just envelope`,
   `cargo nextest run -p xtask -- oracle`, `just lint`, `just test`, `just doc`.

## Open questions

None.

## Further work

- manif on the same launcher (`PHASE1.md` §7).
- The `R(q)` diagonal of 0072's *Further work* now explains failures on `sen3_ad_n2/n3` and
  `sim3_ad`; the sandwich-versus-matrix rounding in `act` on `so3_act` and `sim3_act`.
- `f32` rows, and a Sim(3) Jacobian oracle: GTSAM 4.3.0 exports none.
