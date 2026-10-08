# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [Unreleased]

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

- **`eig3` anchors its frame on the more isolated end** (`0057`): a tie of the top pair no longer leaves every column wrong (`eig:gap-1e-12/top` and `eig:rank1`: `1e16` → `6e7 u`; residual `0.8 |A|` → below `64 sqrt(u) |A|` on `diag(5, 1, 5)`); bit-identical where the top is the isolated end; `f64` +1.7%.
- **`solve_cubic` stops inheriting omnisac's limits** (`0031`): the one-real-root arm pairs its cube roots and takes `-q/(w² + p/3 + v²)` where `p > 0` (`x³ + px − 1` at small `p`: 1.8e10 → 1.7 `u`); the monic cubic is homogenised by exact powers of two, so a scaled cubic keeps its roots (`1e-100`-scaled: no root → 51 `u`) and `p³`, `q²` never underflow into a wrong arm; the leading-coefficient floor is gone, a cubic is `a ≠ 0` with `1/a`, `b/a`, `c/a`, `d/a` finite. `GOLDEN` regenerated; `benches/linalg.rs` added (`f64` 1.01–1.30×, `f32` 1.24–1.37× latency).

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
