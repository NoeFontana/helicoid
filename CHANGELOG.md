# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [Unreleased]

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
- **`cargo xtask compare`** (`just compare`, `0038` decision 4, `0078`): a paired comparison of two in-process subjects per id and stratum (an exact sign test with ties dropped, the geometric mean and both extremes of the per-record error ratio, and the ratio of the two maxima), printed only. Its first twin is `helicoid:unit-r`, `Ad` with `R(q)` in the unit form.
- **The envelope with every oracle**: domination failures go from 148 (0.0.2) to 805, 707 at binary64 and 98 at binary32, with 3 exceptions (`0076`). Most of the rise is strata that had no oracle before, so it counts newly measured gaps more than regressions. The oracles run locally (`just oracles`), not in CI.

### Changed

- **`chol` runs unrolled for `N` from 4 to 9** (`0080`): macros write the loop's recurrence out per column, the same operations in the same order, so `(L, mask)` is bit-identical on every input (tested against the loop, every special value, all lanes). LLVM did not unroll the loop nest past `N = 3`. Measured by `bench-gate --against`: `f64` `N = 6` 57 ns (was 97), `N = 8` 0.60×, `N = 9` 0.68×; `f32` 0.73–0.74×; `Dual<f64, 2>` `N = 6` 0.74×; failing input 0.63–0.83×. `N <= 3` and `N > 9` run the loop unchanged. New bench rows `linalg/*/chol/*`.
- **`SO3::from_matrix` normalizes a column of `K + I`** (`0079`): one `sqrt` and four divisions where it took two and seven; 0.54–0.60 of its time at `f64`, 0.70 at `f32`, also over 64 rotations with an unpredictable pivot. An entry beyond about `2^511` (`f64`) or `2^63` (`f32`) now gives NaN, not a finite quaternion. Output bits change on about half the corpus records: the max error over the corpus falls (`f64` 1.557 → 1.507 `u`, `f32` 2.013 → 1.465), the mean rises 0.03–0.05 `u`.
- **`eig3` leaves Jacobi when it has converged** (`0084`, draft): after each rotation it stops once both entries left are at most `u |A|`, a deterministic test (D16), then one polishing step of the two remaining pairs that keeps `Dual` lanes first-order (a diagonal input with an off-diagonal lane gets its vectors' derivative). Each rotation's next entry reaches the next rotation squared, from `c²` a square root before `c`; columns are normalised by one Newton step. Latency 0.72× / 0.52× `#154`'s on a random spectrum (`f64` / `f32`), 0.40× / 0.25× at a pair `1e-8` apart, 0.28× / 0.21× at a triple. Corpus max 3.06 `u` (`f64`, was 4.76) and 2.67 `u` (`f32`, was 3.24), dominating nalgebra and Eigen on every stratum; `eig3_tests` records `36 u |A|` (was 32) and gains the ignored `measure_eig3`, a double-double large-sample table. Ties may come out in another order among equal eigenvalues; a diagonal matrix keeps its axes.
- **`eig3` is four cyclic Jacobi sweeps** (`0081`), replacing Smith's closed form, whose `acos` lost `u |A|² / gap` on a close pair, up to `sqrt(u) |A|`. Every `eig3` stratum now dominates nalgebra's `SymmetricEigen` (was 4 of 31 at binary64, 3 of 19 at binary32): at most 4.76 `u` (`f64`, was 1.07e8) and 3.24 `u` (`f32`, was 3683); `eig:triple` ties nalgebra's exact answer. Every finite matrix is in the domain (an exact power-of-two scaling; was `1e-100 < m < 1e100`, `1e-12 < m < 1e12`), and `eig3` calls no transcendental. It costs 2.1–2.8× the closed form's latency (1.4× at an `f32` triple eigenvalue), about nalgebra's on a random spectrum. A `Dual` at a multiple of the identity now gives each eigenvalue the derivative of one diagonal entry, where it gave a third of each.
- **The resolved `libm` version is recorded** (`0034` plan step 2, in part): helicoid's conformance rows carry `subject_version` `0.0.2+libm.0.2.16`, and the files `cargo xtask thresholds` generates name it in their header. `libm` stays a caret dependency: an exact pin would force every consumer's graph onto one patch release (`0034`). No output moves.
- **MSRV 1.85** (D17): omnisac, a consumer, is on 1.85, and the workspace already built there unchanged. `just msrv` and the CI job read the new floor.
- **`R(q)` stays homogeneous** (`0076`): the unit form `(1 − 2‖u‖²)I + …` wins `Ad` near the identity and loses it near `π`, and it breaks §12's `to_matrix` contract `(1 + η)R`, so `Quat::to_matrix` is unchanged and its rustdoc says why. `so3_act/q:nonunit/f64` and `sen3_ad_n1/q:nonunit/f64` are excepted as `0058`'s stated drift term.
- **Sim(3)'s `⊖` and `⊕` Jacobians share their coefficient evaluations** (`0073`): the group at `σ` and `-σ` from one branch and one set of transcendentals, each to the bit, so `rminus_jacobians` and `lminus_jacobians` run one pair where they ran three evaluations, and `rplus_jacobians`, `lplus_jacobians` one where they ran two. No output bit moves. `rminus_jacobians` takes 0.65–0.84 of its time (`f64` 270–294 ns, was 387–468), `rplus_jacobians` 0.73–0.88 and `jr_inv` 0.74–0.84.
- **`Matrix × Matrix` is a plain loop** (`0069`): the nested `array::from_fn` built the product through the stack, 87 stores for a `3 x 3`; the same sums in the same order, so every output is bit-identical. `Jac::mul` on SO(3), SE(3) and SE₂(3) takes 0.09–0.11 of its former time, SE_N(3) `jr_inv` and `rminus_jacobians` 0.40–0.60.
- **`Vector::dot` is a plain loop** (`0067`): the same left-to-right sum, so every output is bit-identical, but `f32` callers no longer reload each lane through an `array::IntoIter` on the stack. `f32` `S2Chart::local` takes 0.13–0.21 of its former time and the SE₂(3) Gaussian gate 0.87; a few `f64` rows move by ±5% with register allocation.

## [0.0.2] - 2026-10-09

A patch release: the geodesic below `r`'s second switch is faster and more accurate (`0059`), and the carried quaternion's domain is stated (`0058`). No API change.

### Changed

- **`SO3::from_quat_unchecked`'s stated domain is the drift band** (`0058`): a struct literal moved into it is the *carried* path, `|‖q‖² − 1| ≤ 2^-26.29`, nothing asserted. `NUMERICS.md` §12 states each SO(3) operation's first-order error there, and §3.6 states the drift budget. No code changes.
- **`SO3::geodesic` and `SE3::geodesic` below `r`'s second switch are faster and more accurate** (`0059`): the provided body is written out so that its two series no longer wait on its products, and `horner` no longer starts from `0 · z`. `so3_geodesic` and `se3_geodesic` at `geo:consecutive` go from 1.572 to 0.993 `u`; every other conformance row is byte-identical.

## [0.0.1] - 2026-10-08

First release: `helicoid-linalg` and `helicoid` at `0.0.1` (`0049`).

### Added

- **Corpus ids for the routines D7 did not reach** (`0056`): `solve_cubic`, `eig3`, `chol_n{3,6}`, `chol_solve_n{3,6}`, `quat_renormalize` and `real_*` (`Dual`'s derivatives), at both precisions; root-set, gap-weighted eigenvector and mask metrics; a nalgebra oracle runner (`just oracle-nalgebra`); 198 committed per-stratum maxima reproduced bit for bit by the shipped subject.

- **`Real` transcendental primitives (`sin`, `cos`, `acos`)**: added to `Real` via `libm` to avoid discarding unneeded components from `sin_cos` (`0022`, `0052`, `0053`).
- **`SEn3` constructors for general $N$**: `from_parts` and `parts` by value through `SO3` to preserve unit invariants across all orders (`0042`).
- **`Product` constructors**: `Product::from_parts` and `parts` to enable composition without `Exp`/`Log` roundtrips (`0029`, `0040`).
- **`SEn3` drift repair and relative transforms**: `SEn3::renormalize` for explicit drift repair, and `SEn3::{mul_inv, inv_mul}` to optimize relative tree compositions (`0044`, `0048`).
- **Geodesic API and velocity**: `LieGroup::geodesic` and `geodesic_velocity` with structured Jacobians (`PHASE4.md` §1.1, `NUMERICS.md` §10).
- **Domination exceptions harness**: `conformance/baseline/exceptions.toml` allows excepting ten-digit numerical ties from domination failures (`0046`).
- **Symbol resolution linter**: `cargo xtask lint` verifies that test symbols cited in comments resolve to definitions (`0048`).
- **Geodesic conformance instruments**: `so3_geodesic` and `se3_geodesic` corpus datasets (180 records each) for conformance validation (`PHASE4.md` §4, `0045`).
- **Second series arm in coefficient kernels**: short polynomial prefixes below a second switch reduce Horner evaluations for small angles while preserving exact bit agreement (`0039`, `0047`).
- **Switch reference points (`coeff_switch_ref.jsonl`)**: exact reference values generated at grid points to verify continuity bounds at switch thresholds (`0039`).
- **Group benchmarks**: criterion benchmark suites for group operations in `crates/helicoid/benches/groups.rs` (`PHASE3.md` §11).
- **Dual-matrix Jacobian verification**: proptest suite verifying every group Jacobian row against `Dual<S, DOF>` automatic differentiation (`PHASE3.md` §8).
- **Linear algebra routines in `helicoid-linalg`**: `Real::cbrt`, `solve_cubic`, and symmetric `eig3` (`PHASE2.md`, `0017`).
- **Cholesky factor solve**: `chol_solve(&L, b)` solving from factors with column reads (`0019`).
- **Conformance and oracle runner harnesses**: runner harnesses for `tf_tree_math` and `sophus-rs` (`PHASE1.md` §7, `0010`), plus `cargo xtask envelope` judging domination and regression bars (`PHASE1.md` §8, `0006`).

### Changed

- **`eig3` anchors its frame on the more isolated end** (`0057`): a tie of the top pair no longer leaves every column wrong (`eig:gap-1e-12/top` and `eig:rank1`: `1e16` → `6e7 u`; residual `0.8 |A|` → below `64 sqrt(u) |A|` on `diag(5, 1, 5)`); the closing vector is no longer renormalised; `f64` 0.95× latency.
- **`solve_cubic` stops inheriting omnisac's limits** (`0031`): the one-real-root arm pairs its cube roots and takes `-q/(w² + p/3 + v²)` where `p > 0` (`x³ + px − 1` at small `p`: 1.8e10 → 1.7 `u`); the monic cubic is homogenised by exact powers of two, so a scaled cubic keeps its roots (`1e-100`-scaled: no root → 51 `u`) and `p³`, `q²` never underflow into a wrong arm; the leading-coefficient floor is gone, a cubic is `a ≠ 0` with `1/a`, `b/a`, `c/a`, `d/a` finite. `GOLDEN` regenerated; `benches/linalg.rs` added (latency 0.97–1.13× of before on every row it already solved; a cubic it used to reject as out of range, about 100 ns).

- **`Real::branch` is inlined; `abs`/`copysign` are `core`'s**: the coefficient kernels run 0.55× and the group operations 0.915× (geomean) of before, `so3/exp` 0.57×, results bit-identical (`0055`).
- **`SE3::geodesic` is the screw twin**: `SO3::geodesic`'s rotation, and the translation in the world frame — the definition's closed form below `r`'s second switch, GE.12's dual-quaternion power above it. `se3_geodesic` dominates `tf_tree_math`'s `ScLerp` on all three strata, 0.59–0.92× the provided body's latency, 4.6–6.1% faster than `screw_pow` inside `tf_tree`'s lookup bench; through `Dual` it differentiates like the provided body below the switch. `SE₂(3)` keeps the provided body (`PHASE4.md` §1.2, `0054`).
- **`SO3::geodesic` two-arm dispatch**: provided body below `log_ratio` short switch and blend above, preserving near-identity performance without accuracy loss (`0051`, superseding `0050`).
- **`SO3::geodesic` blend denominator**: denominator recomputed as $\sin(\alpha)$ rather than $\|\mathbf{v}\|$ to ensure endpoint exactness and dominate oracle baselines (`0050`).
- **Geodesic Jacobian cancellation-free form**: evaluated via $(1-t) J_l((1-t)d) J_l^{-1}(d)$ with `scale(k)`, eliminating $1/(1-t)$ numerical cancellation (`0043`).
- **Trigonometric solver backends**: `solve_cubic` and `eig3` route through `libm::acos` and `libm::cos` for improved latency and accuracy (`0022`, `0053`).
- **`SEn3Jac` internal block representation**: fields scoped to `pub(crate)` to keep representation private while supporting internal operations (`0028`, `0040`).
- **Threshold sweep domain expansion**: search grid lifted across the entire domain up to $\pi^2$, improving accuracy across all coefficient series (`0039`).
- **Benchmark gate loop order**: replicate loop moved outermost to minimize thermal and scheduling drift across comparison pairs (`0033`).
- **Structural zero skipping**: hat-structured matrix multiplications skip structural zeros to optimize throughput.
- **Relative transform Jacobian fusion**: `SEn3` `rminus_jacobians` and `lminus_jacobians` fuse inverse operations to avoid redundant allocations.
- **Determinism through `libm` `arch` feature**: target-optimized hardware instructions enabled for exactly rounded primitives (`0018`).
- **Threaded `dot_acc`**: `Tangent::dot` implemented via `dot_acc` to ensure flat sum semantics across nested products (`0025`).

### Fixed

- **Phase 4 geodesic reference validation near $\pi$**: reference computation switches to geometric quaternion `Log` where `mp.logm` branches into complex preimages (`0045`).
- **Product right-invariance check**: corrected to verify SE(3) cross-coupling rather than abelian product invariance (`0045`).
- **Scale-relative error measurement**: `laws::e_at` computes relative errors against normalized denominators.
- **`Dual::acos` domain protection**: safe argument selection prevents singular derivatives at boundaries.

### Removed

- **Unused `Quat` methods (`dot`, `norm`)**: removed in favor of `rminus` and `norm_sq` to prevent unvalidated domain assumptions (`0048`).
- **Static benchmark noise floors**: stored host noise floor replaced by active control measurements beside comparisons (`0033`).
