# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Released
entries live in [`docs/changelog/`](./docs/changelog/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [Unreleased]

### Added

- **`just determinism`** (`PHASE6.md` §3, D16): `cargo xtask determinism` digests the shipped subject's outputs over the whole corpus, one SHA-256 per id and precision, NaN canonicalized. CI digests x86_64, aarch64 and `wasm32-wasip1` and fails on any difference. At `0.0.3` all 114 digests agree on x86_64 and `wasm32-wasip1`.

- **`Real::round_from_f64(x: f64) -> Self`** (`0087`): a runtime `f64` correctly rounded to nearest-even at the scalar's precision, the identity at `f64`, `as f32` at `f32` (IEEE's conversion, the same bits on every target), a constant for `Dual`. Every `f64` is in the domain: at `f32` it overflows to `±inf` from `(2 − 2⁻²⁴)·2¹²⁷` and reaches `±0` at `2⁻¹⁵⁰`. A provided method, so an `impl Real` outside the workspace keeps compiling. Corpus id `real_round_from_f64`, with binary32 ties, near-ties and both edges: every stratum reads only the reference's 30-digit rounding.

### Changed

- **`SO3::from_matrix` returns Shepperd's quotients** (`0086`), in `tf_tree_math`'s arithmetic, whenever they land within `16u` of the unit sphere, which every rotation rounded to the precision does. A scaled, non-orthogonal or NaN matrix still gets 0.0.3's normalized column. Exact rotations, such as quarter turns and frame conventions' signed permutations, come back with `fl(‖q‖²) = 1` at `f64` (0.0.3: `1 − 2u`, which every later `act` carried). Under `0085`'s score, strata above the best oracle fall from 18 to 8 of 31 at `f64` and from 14 to 9 at `f32`, and `f64` now equals `tf_tree_math` on every stratum's maximum but `q:nonunit` (89.4 against 163.0 `u`). Six `f64` strata near `π` read up to 0.50 `u` higher than 0.0.3.

## [0.0.3] - 2026-10-13

Phase 5: charts, S², Sim(3), the integrated exponentials and `Gaussian`. `eig3` is Jacobi and dominates nalgebra and Eigen; `chol` and `SO3::from_matrix` are faster; three C++ container oracles. MSRV 1.85.

### Breaking

- **`Real` gains two required methods, `exp` and `exp_m1`** (`0068`): an implementor of `Real` outside this workspace adds them. No public item was removed or changed.

### Added

- **Charts** (`PHASE5.md` §1, `0012`, `0060`): `Chart<S, M>`, `Manifold<S>`, `WithChart<M, C>` (its chart `Lifted<C>`), `RightChart<G>` and `LeftChart<G>` on every group. The SE(3) charts `Screw<S>`, `Decoupled<S>` (locus-tag's `Pose::retract`) and `WorldTranslation<S>` share the tangent `Twist<S>`; the last two have `TwistBlockJac<S>` Jacobians, and `SE3::chart_transition::<From, To>()` gives `DΦ(0)` (`I`, `diag(I, R)`, `diag(I, Rᵀ)`) for moving a covariance between them. `Chart::local_with_jacobian` returns a residual and its Jacobian with the `Log` taken once, bit-identical to the two calls. Every chart Jacobian is checked against `Dual`, with both CH.6 identities. Corpus ids `se3_{screw,decoupled,world}_{retract,local}`: at most 5.94 `u`.

- **`helicoid::sinc(θ²) -> (sin θ/θ, d/dθ²)`** and **`helicoid::sinc_value(θ²)`** (`0062`): the value `2k cos(θ/2)` of `exp_coeffs`, the derivative the `Dual` lane of the same evaluation; `sinc_value` is the value alone on plain `S`, bit-identical, without the `Dual` cost. No new switch, `coeffs/generated.rs` unchanged; accuracy stated for `θ ≤ π`, the range the corpus scores. Corpus id `coeff_alpha` (the `theta:*` strata at both precisions): at most 2.59 `u` (`f64`) and 2.30 `u` (`f32`).

- **`Mat2::inverse_adj`** (`0061`): `(adj/det, det)` with `Mat3::inverse_adj`'s contract; the adjugate is exact, so `det` is the only rounding. Corpus id `mat2_inverse_adj` (`cond:1`, `1e4`, `1e8`, an exactly singular stratum; `cond:1e8` at binary64 only), scored with `NUMERICS.md` §11's expected-non-finite rule; nalgebra's `Matrix2::try_inverse` answers it in the nalgebra runner, equal on every stratum.

- **Integrated exponentials** (`PHASE5.md` §4, `0064`):
  - `helicoid::so3::{gamma1, gamma2, gamma_apply_jacobians::<M>}`. `Γ₁` is `SO3::jl` to the bit; `Γ₂ = ½I + bW + dW²` comes from a new group `gamma2_coeffs` over the swept switches of `b` and `d`, with no new switch.
  - `gamma_apply_jacobians` returns `(Γ_M v, ∂(Γ_M v)/∂φ, Γ_M)` for `M` 1 or 2 from one coefficient evaluation on `Dual<S, 1>` and the chain rule (`0066`). Its twin is `reference::gamma_apply`, whose `m` is bounded so that `m!` is exact.
  - Corpus id `so3_gamma2` at both precisions: at most 1.46 `u` (`f64`) and 1.68 `u` (`f32`). A planted Γ₂-from-Γ₁ defect is detected.

- **`Gaussian<S, G, Sd, D>`** (`PHASE5.md` §5, `0065`):
  - The side is in the type; `new` is the only constructor, and it asserts `D == G::DOF`.
  - `new` reads `Σ` by its lower triangle; the field is private and `cov()` returns it exactly symmetric (`0066`).
  - `whitener()` factors `Σ` once for many-point gating: `Whitener::mahalanobis_sq` is `mahalanobis_sq` to the bit (`0066`).
  - `to_left`/`to_right` stay within `γ_{4D+42}(|B||A|)|Σ|(|B||A|)ᵀ` on the round trip.
  - `propagate` documents its side contract; `mahalanobis_sq` returns `chol`'s mask.
  - Corpus ids `gaussian_mahalanobis_{se3,se23}`, each stratum within 0.44 of its componentwise condition number.

- **S² and its chart** (`PHASE5.md` §2, `0063`, `0067`, `0082`): `S2<S>`, held unit by `from_vec_unchecked`, `from_vec_normalized` and `renormalize`, and `S2Chart<S>`, its Householder basis frozen per linearization and oriented (`b₂ = ς H e_y`, so `b₁ × b₂ = n` on both hemispheres: a bearing behind a camera has the handedness of one in front), exposed as `rotation_basis()` (`B`, `δ`'s rotation vector `Bδ`) and `displacement_basis()` (`T = [-b₂ b₁] = ∂retract/∂δ` at `0`, `Tᵀ = ∂local/∂m` at `n`), both `Matrix<S, 3, 2>`, on `RnTangent<S, 2>` coordinates with `Mat2<S>` Jacobians (`Mat2<S>: Jac<S, RnTangent<S, 2>>`); `S2` is a `Manifold`. `ς` is read by comparison, so `n_z = ±0` give one basis. `local` takes `α/s` from `log_ratio`, and its Jacobian forms `Exp` of the rotation vector rather than dividing by `1 + n·m`. The twin `reference::s2_local`. Corpus ids `s2_retract` and `s2_local` over `s2:nz0`, `s2:near-antipode` and `s2:generic` at both precisions, `s2_retract` scored by `NUMERICS.md` §11's new direction metric: at most 6.03 `u` and 1.89 `u` off the antipode.

- **`Real::exp` and `Real::exp_m1`** (`0068`), through `libm`, with `Dual` rules (`exp_m1`'s derivative from `exp`, so it keeps its digits past `-37.4`). Corpus ids `real_exp` and `real_expm1` at both precisions: at most 1.22 `u`. Sim(3)'s first prerequisite.

- **Sim(3)** (`PHASE5.md` §3, `0070`): `Sim3<S>`, `Sim3Tangent<S>` and the structured `Sim3Jac<S>` (`[A 0 0; B C f; 0 0 e]`, closed under product and inverse, `117` multiplications against `343` dense), with dense twins `reference::sim3jac_*`. `Exp`, `Log` (no `mp.logm`), `Ad`, `ad`, the Jacobians and their inverses by blocks, finite to `σ = ±700` by the scaled arrangement. Its coefficients are one group: the total-degree series below a generated `r_×²`, closed forms above, both switched by a new ray sweep in `cargo xtask thresholds` (`f64` 4.72 `u`, `f32` 2.93 `u` on the Jacobian blocks). Corpus ids `sim3_{exp, log, jr, jl, jr_inv, jl_inv, ad, act}` at both precisions: the Jacobians within 3.64 `u` on the joint-limit rays to `r_Δ = 10⁻¹²`, every value finite to `σ = ±700`. The corpus budget is 64 MiB (`0071`).

- **The Sophus C++ oracle** (`PHASE1.md` §7, `0072`), the first container runner: `cargo xtask conformance --oracle sophus` (`just oracle-sophus`) builds `docker/oracles/sophus`, pinned by base digest, Ubuntu archive snapshot and Sophus commit, and runs it offline on the existing file protocol. It answers 26 ids at binary64 across SO(2), SE(2), SO(3), SE(3) and Sim(3), its convention conversions tested by hand at image build. Sim(3) has its first oracle: `helicoid` is not dominated on 87 of the 434 paired strata, so `PHASE5.md` §7 stays open.

- **The GTSAM C++ oracle** (`PHASE1.md` §7, `0074`): `just oracle-gtsam` builds `docker/oracles/gtsam`, GTSAM 4.3.0 pinned by commit and checksum and compiled in the image under the runner's flags (no contraction, no Eigen SIMD), and runs it offline. It answers 47 ids at binary64: SO(2), SE(2), SO(3) with `Γ₂`, SE_N(3) for `N` = 1, 2, 3, the SE(3) screw chart and Sim(3) `exp`/`log`/`ad`/`act`. `sen3_*_n2`, `sen3_*_n3`, `so3_gamma2` and the screw chart get their first oracle (919 more strata paired); Sim(3)'s domination failures go from 87 of 434 to 151 of 485, so `PHASE5.md` §7 stays open.
- **The manif oracle** (`PHASE1.md` §7, `0075`), the second container: `just oracle-manif` builds `docker/oracles/manif`, manif 0.0.5 pinned by commit and checksum on the Sophus image's base and snapshot, and answers 36 ids at binary64 across SO(2), SE(2), SO(3), SE(3), SE₂(3) and the SE(3) screw and world charts, every Jacobian by manif's own formula. SE₂(3) (`sen3_*_n2`) and the charts have their first oracle: `helicoid` is not dominated on 146 of the 576 newly paired strata, all within 3.52×.
- **Oracles answer at the precision they are asked** (`0077`): the runner protocol takes `--precision f32`, and an `f32` run scores binary32 answers as an in-process `f32` subject is scored. nalgebra answers its ids and Sophus `sim3_{exp, log, ad, act}` in binary32; GTSAM, manif, `tf_tree_math` and sophus-rs stay binary64. Sim(3) at `f32` fails 74 of 428 paired strata, by at most 2.53×.
- **Eigen judges `eig3`** (`0083`): the Sophus image answers `eig3` with Eigen 3.4.0's `SelfAdjointEigenSolver<Matrix3<T>>`, at both precisions, a second oracle on all 50 strata beside nalgebra; `computeDirect` was measured once and is not a subject.
- **`cargo xtask compare`** (`just compare`, `0038` decision 4, `0078`): a paired comparison of two in-process subjects per id and stratum (an exact sign test with ties dropped, the geometric mean and both extremes of the per-record error ratio, and the ratio of the two maxima), printed only. Its first twin is `helicoid:unit-r`, `Ad` with `R(q)` in the unit form.
- **The envelope with every oracle**: domination failures go from 148 (0.0.2) to 805, 707 at binary64 and 98 at binary32, with 3 exceptions (`0076`). Most of the rise is strata that had no oracle before, so it counts newly measured gaps more than regressions. The oracles run locally (`just oracles`), not in CI.

### Changed

- **`chol` runs unrolled for `N` from 4 to 9** (`0080`): macros write the loop's recurrence out per column, the same operations in the same order, so `(L, mask)` is bit-identical on every input (tested against the loop, every special value, all lanes). LLVM did not unroll the loop nest past `N = 3`. Measured by `bench-gate --against`: `f64` `N = 6` 57 ns (was 97), `N = 8` 0.60×, `N = 9` 0.68×; `f32` 0.73–0.74×; `Dual<f64, 2>` `N = 6` 0.74×; failing input 0.63–0.83×. `N <= 3` and `N > 9` run the loop unchanged. New bench rows `linalg/*/chol/*`.
- **`SO3::from_matrix` normalizes a column of `K + I`** (`0079`): one `sqrt` and four divisions where it took two and seven; 0.54–0.60 of its time at `f64`, 0.70 at `f32`, also over 64 rotations with an unpredictable pivot. An entry beyond about `2^511` (`f64`) or `2^63` (`f32`) now gives NaN, not a finite quaternion. Output bits change on about half the corpus records: the max error over the corpus falls (`f64` 1.557 → 1.507 `u`, `f32` 2.013 → 1.465), the mean rises 0.03–0.05 `u`.
- **`eig3` leaves Jacobi when it has converged** (`0084`): after each rotation it stops once both entries left are at most `u |A|`, a deterministic test (D16), then one polishing step of the two remaining pairs that keeps `Dual` lanes first-order (a diagonal input with an off-diagonal lane gets its vectors' derivative). Each rotation's next entry reaches the next rotation squared, from `c²` a square root before `c`; columns are normalised by one Newton step. Latency 0.72× / 0.52× `#154`'s on a random spectrum (`f64` / `f32`), 0.40× / 0.25× at a pair `1e-8` apart, 0.28× / 0.21× at a triple. Corpus max 3.06 `u` (`f64`, was 4.76) and 2.67 `u` (`f32`, was 3.24), dominating nalgebra and Eigen on every stratum; `eig3_tests` records `36 u |A|` (was 32) and gains the ignored `measure_eig3`, a double-double large-sample table. Ties may come out in another order among equal eigenvalues; a diagonal matrix keeps its axes.
- **`eig3` is four cyclic Jacobi sweeps** (`0081`), replacing Smith's closed form, whose `acos` lost `u |A|² / gap` on a close pair, up to `sqrt(u) |A|`. Every `eig3` stratum now dominates nalgebra's `SymmetricEigen` (was 4 of 31 at binary64, 3 of 19 at binary32): at most 4.76 `u` (`f64`, was 1.07e8) and 3.24 `u` (`f32`, was 3683); `eig:triple` ties nalgebra's exact answer. Every finite matrix is in the domain (an exact power-of-two scaling; was `1e-100 < m < 1e100`, `1e-12 < m < 1e12`), and `eig3` calls no transcendental. It costs 2.1–2.8× the closed form's latency (1.4× at an `f32` triple eigenvalue), about nalgebra's on a random spectrum. A `Dual` at a multiple of the identity now gives each eigenvalue the derivative of one diagonal entry, where it gave a third of each.
- **The resolved `libm` version is recorded** (`0034` plan step 2, in part): helicoid's conformance rows carry `subject_version` `0.0.2+libm.0.2.16`, and the files `cargo xtask thresholds` generates name it in their header. `libm` stays a caret dependency: an exact pin would force every consumer's graph onto one patch release (`0034`). No output moves.
- **MSRV 1.85** (D17): omnisac, a consumer, is on 1.85, and the workspace already built there unchanged. `just msrv` and the CI job read the new floor.
- **`R(q)` stays homogeneous** (`0076`): the unit form `(1 − 2‖u‖²)I + …` wins `Ad` near the identity and loses it near `π`, and it breaks §12's `to_matrix` contract `(1 + η)R`, so `Quat::to_matrix` is unchanged and its rustdoc says why. `so3_act/q:nonunit/f64` and `sen3_ad_n1/q:nonunit/f64` are excepted as `0058`'s stated drift term.
- **Sim(3)'s `⊖` and `⊕` Jacobians share their coefficient evaluations** (`0073`): the group at `σ` and `-σ` from one branch and one set of transcendentals, each to the bit, so `rminus_jacobians` and `lminus_jacobians` run one pair where they ran three evaluations, and `rplus_jacobians`, `lplus_jacobians` one where they ran two. No output bit moves. `rminus_jacobians` takes 0.65–0.84 of its time (`f64` 270–294 ns, was 387–468), `rplus_jacobians` 0.73–0.88 and `jr_inv` 0.74–0.84.
- **`Matrix × Matrix` is a plain loop** (`0069`): the nested `array::from_fn` built the product through the stack, 87 stores for a `3 x 3`; the same sums in the same order, so every output is bit-identical. `Jac::mul` on SO(3), SE(3) and SE₂(3) takes 0.09–0.11 of its former time, SE_N(3) `jr_inv` and `rminus_jacobians` 0.40–0.60.
- **`Vector::dot` is a plain loop** (`0067`): the same left-to-right sum, so every output is bit-identical, but `f32` callers no longer reload each lane through an `array::IntoIter` on the stack. `f32` `S2Chart::local` takes 0.13–0.21 of its former time and the SE₂(3) Gaussian gate 0.87; a few `f64` rows move by ±5% with register allocation.

## Released versions

Frozen history, newest first, in [`docs/changelog/`](./docs/changelog/); leave it
out of searches and sweeps.

- [`0.0.2`](./docs/changelog/0.0.2.md) — 2026-10-09 (the carried path and a faster geodesic)
- [`0.0.1`](./docs/changelog/0.0.1.md) — 2026-10-08 (first release)

[0.0.1]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.1
[0.0.2]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.2
[0.0.3]: https://github.com/NoeFontana/helicoid/releases/tag/v0.0.3
