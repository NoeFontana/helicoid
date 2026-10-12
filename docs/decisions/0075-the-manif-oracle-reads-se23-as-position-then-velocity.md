# 0075: the manif oracle reads SE₂(3) as position, then rotation, then velocity

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #144

## Context

`PHASE1.md` §7 owes `docker/oracles/manif` ("header-only, small C++ runner"), and
[`0072`](./0072-a-container-oracle-is-an-image-on-the-same-protocol.md) built the launcher it
runs on (`Launch::Container`). No runner answers an SE₂(3) id: the `sen3_*_n2` strata, and the
SE(3) charts' `se3_{screw,world}_{retract,local}`, have no oracle.

manif 0.0.5 (`c49fec8`, its latest release) is header-only and vendors `tl::optional`. For ids the
corpus has, it exports `SO2`, `SE2`, `SO3`, `SE3` and `SE_2_3`, each with `exp`, `log`, `adj`, and
for every tangent `ljac`, `rjac`, `ljacinv`, `rjacinv`, each its own formula (SO(3)'s right ones
are the transposes); `SO3::act`; `rplus`/`rminus` on every group, including
`Bundle<SO3, R3>`; and `interpolate_slerp`, `X₀.rplus(X₁.rminus(X₀) · t)`, which is
`NUMERICS.md` §10's `γ`. It has no Sim(3), no SE₃(3), no quaternion of a matrix and nothing
that computes the decoupled chart's `x + Rρ`.

Its conventions, read from the source and the hand tests below: tangents translation-first,
SE(2) and SE(3) `[ρ; θ]`; **SE₂(3) `[ρ; θ; ν]`**, the rotation between the two translations
(`SE_2_3Base::log` writes `[J_l⁻¹ t; θ; J_l⁻¹ v]`), its element `(t, q, v)` with `t` the position
and `v` the velocity, composed as helicoid's (`(R₁R₂, R₁t₂ + t₁, R₁v₂ + v₁)`). helicoid's SE₂(3)
is `(R, v, p)` with tangent `[φ; ρ_v; ρ_p]` (`docs/maths/se3.md` SE.1). Quaternions are stored
`[x, y, z, w]`. Every constructor checks the unit norm, `|‖q‖ − 1| < 100 ε`, and throws
`manif::invalid_argument` when it fails; it never normalizes. The check is `MANIF_ASSERT`, which
`NDEBUG` or `MANIF_NO_DEBUG` compiles out.

## Decision

1. **The image** (`docker/oracles/manif/Dockerfile`): `0072`'s base digest and archive snapshot
   (gcc 13.3.0, Eigen 3.4.0, nlohmann-json 3.11.3, glibc 2.39); manif 0.0.5 by commit
   `c49fec8150e62821d9be7e2a3ac49db75bf7bb92` and the SHA-256 of its archive, keeping `include/`
   and `external/tl`. Compiled `-std=c++17 -O2 -ffp-contract=off -DEIGEN_DONT_VECTORIZE`, warnings
   denied, **without `NDEBUG`** (`convert.hpp` refuses to compile with it). The build compiles and
   runs `convert_test`, then the runner. `--version` prints
   `manif@0.0.5+<commit> eigen@… gcc@… glibc@…`.
2. **The ids** (`Runner::answers`, 36), every one manif computes and the corpus holds:
   `so2_{exp,log}`, `se2_{exp,log,ad,jl,jr,jl_inv,jr_inv}`,
   `so3_{exp,log,act,jl,jr,jl_inv,jr_inv,geodesic}`, `se3_geodesic`,
   `se3_{screw,world}_{retract,local}`, and `sen3_{exp,log,ad,jl,jr,jl_inv,jr_inv}` at `_n1`
   (`SE3`) and `_n2` (`SE_2_3`). Each Jacobian is manif's own function for it, not another at
   `−τ`. The screw chart is `SE3::rplus`/`rminus`; the world chart is `Bundle<SO3, R3>`'s, whose
   tangent is `[φ; ρ]`. The geodesics are `interpolate_slerp`. Backend `HostStd`.
3. **The conversions** (`convert.hpp`), named after the foreign convention, each tested by hand in
   `convert_test.cpp`:
   - quaternion: helicoid `[w, x, y, z]` is `Eigen::Quaterniond(w, x, y, z)`;
   - tangent: helicoid's component `i` sits at manif's `at[i]`, with `at` = `[2, 0, 1]` for SE(2),
     `[3, 4, 5, 0, 1, 2]` for SE(3) and **`[3, 4, 5, 6, 7, 8, 0, 1, 2]` for SE₂(3)**
     (`φ → θ`, `ρ_v → ν`, `ρ_p → ρ`); a Jacobian or adjoint's helicoid `(i, j)` is manif's
     `(at[i], at[j])`, written column-major;
   - elements: `SO2(z₀, z₁)`, `SE2(t₀, t₁, z₀, z₁)`, `SE3(x, q)`, and helicoid's `(q, [v; p])` as
     `SE_2_3(p, q, v)`; back, `x = [linearVelocity(); translation()]`.
4. **A refused record is answered `nan`**, as `0072` decision 5: the runner catches
   `manif::invalid_argument` and `manif::runtime_error` (`interpolate_slerp` checks
   `t ∈ [0, 1]`).
5. **Shared helpers are copied.** An image's build context is its own directory (`0072`
   decision 1), so `hexfloat.hpp` is a byte-identical copy of Sophus's, held so by the xtask test
   `the_container_runners_share_one_hex_float_codec`, and the runner's file half (`inputs`,
   `answer`, `answer_file`, `run`) is Sophus's with the names changed.
6. **Binary64 only**, as every runner.

## Rationale

- **The SE₂(3) map follows the labels, not the slots.** manif's group law and formulas treat `t`
  and `v` alike, so mapping helicoid's `x₁` to manif's `t` would score the same; but helicoid's
  `x₁` is the velocity by convention (SE.1) and manif's `v` is its velocity, and a hand test that
  fixes the labels (`se23: manif's translation is helicoid's p = x₂`) fails on the other map,
  where a symmetric test could not.
- **The checks stay on.** With `NDEBUG` manif carries a non-unit quaternion silently; with the
  checks, a record outside the domain manif states is a refusal, as Sophus's ensures are in
  `0072`. No corpus record trips one: the `q:nonunit` strata drift less than `100 ε`.
- **The world chart through `Bundle`**: it is manif's own product manifold, whose `rplus` is
  `(q ∘ Exp φ, x + ρ)` with manif's `compose`; spelling it as an `SO3` and a vector add in the
  runner would be our formula, not manif's.
- **A copy over a shared context**: widening the build context to `docker/oracles` would change
  `0072`'s launcher and every image's cache key for one 30-line header; the test makes the copy
  as safe as a share.

## Consequences

Measured on `d4856a5` (`feat/oracle-sophus`) with the five oracles' rows, `just envelope` (no
baseline in the tree, so no no-regress verdict):

- **Pairing.** 1082 → 1658 binary64 strata have an oracle; the 576 new are `sen3_*_n2` (368) and
  `se3_{screw,world}_*` (208). `so2_*`/`se2_*` rows wait for SO(2)/SE(2) (`PHASE3.md` §0.0).
  No record is refused.
- **Domination failures**: 301 → 468, all 167 new ones over manif: `sen3_*_n2` 71 (`ad` 22,
  `exp` 14, `jl` 12, `jr` 12, `log` 8, `jl_inv` 2, `jr_inv` 1), the charts 75
  (`se3_world_local` 26, `se3_screw_retract` 22, `se3_world_retract` 18, `se3_screw_local` 9),
  and in strata Sophus or `tf_tree_math` already paired, `so3_act` 12, `so3_exp` 3,
  `sen3_jl_inv_n1` 2, `so3_log` 2, `sen3_exp_n1` 1, `sen3_log_n1` 1. 3 close under
  `helicoid:host-std` (`sen3_log_n2` at `rho:{1e-6,1e-3,1e3}/theta=1`). manif becomes the best
  oracle of 20 failures that already failed. SE₂(3) fails 71 of its 368 strata and the charts 75
  of 208.
- **One is not near a tie**: `sen3_ad_n2/theta:1e-2`, 1.565 against 0.4450 `u` (3.52×). Replaying
  its argmax record, `Ad` with the homogeneous `R(q)` of `NUMERICS.md` §1 reads 1.565 `u`, 1.289
  of it in the three `R` blocks, and with Eigen's `1 − 2(y² + z²)` diagonal 0.373; manif's argmax
  record replays to its 0.4450 with Eigen's. It is `0072`'s near-identity `R(q)` loss, and
  `sen3_ad_n1/theta:1e-2` already fails by 3.30× against `tf_tree_math`, Sophus and manif alike.
  The other 166 are within 2.79×, on maxima of at most 7.80 `u`
  (`sen3_jr_n2/rho:1e4/theta=pi-1e-6`, against 3.80).
- **manif is wrong where it is not best, as §7 allows**: its small-angle `SO3::log` ignores the
  sign of `w`, so `so3_log/theta:1e-12` reads `1.8·10¹⁶ u` on the `w < 0` records; its `J_l⁻¹`
  reads `4·10⁷ u` at `θ = π − 10⁻⁸` and SE(3)'s `J_l` `9·10⁶ u` at `θ = 10⁻⁴` (`(1 − cos θ)/θ²`
  and its relatives in closed form down to `θ² = 100 ε`). Those rows are recorded and set no bar.

## Implementation plan

1. This record, `docker/oracles/manif`, the `Runner` entry and the hex-codec test,
   `just oracle-manif`, the status rows. Verified by the image build (it runs `convert_test`),
   `cargo nextest run -p xtask -- oracle`, `just oracle-manif`, `just envelope`, `just lint`,
   `just test`, `just doc`.

## Open questions

None.

## Further work

- The GTSAM container on the same launcher; then `just oracles`.
- `sen3_*_n2`'s 71 and the charts' 75 failures, on maxima under 8 `u`: which are `R(q)`'s
  (`0072` *Further work*) and which are not, stratum by stratum.
- The decoupled chart still has no oracle.
- `f32` rows, as `0072`: manif is templated on the scalar.
