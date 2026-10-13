# Architectural decisions

A significant architectural change starts as a decision record in this folder, not as a PR.
Skeleton: [`template.md`](./template.md).

Cite [`docs/PROJECT.md`](../PROJECT.md) (D1–D18 in §5) and [`docs/NUMERICS.md`](../NUMERICS.md);
[`docs/API.md`](../API.md) is what a public-surface decision is checked against. A `path.rs:NNN`
citation is unreliable: cite the symbol.

The table does not restate status; read it from the record:

```sh
grep -m1 -H '^\*\*Status:' docs/decisions/0*.md
```

| Record | Decided |
|---|---|
| [`0001`](./0001-record-architectural-decisions.md) | the meta-decision that this folder exists. |
| [`0002`](./0002-one-convention-for-a-stack-that-already-disagrees.md) | Hamilton, `w` first, active, rotation-first tangents, right perturbation by default, sides always spelled; converters for every other order. |
| [`0003`](./0003-the-scalar-that-cannot-say-less-than.md) | `Real` without `PartialOrd`; masks, `branch`, `select`, the safe argument; `Dual` inherits the value path bitwise. |
| [`0004`](./0004-switch-points-are-generated-not-typed.md) | switch points and series lengths generated per coefficient and precision from a committed sweep; drift-gated. |
| [`0005`](./0005-the-jacobian-is-a-dual-matrix.md) | structured Jacobians: SE_N(3) as dual matrices, `Jr⁻¹` by algebra, dense only as a write; stable-Rust const-generic workarounds. |
| [`0006`](./0006-the-instrument-comes-first.md) | mpmath corpus from definitions, strata, the oracle envelope, domination and exact no-regress bars, seeded defects. |
| [`0007`](./0007-the-budget-a-foundation-can-afford.md) | `libm` (+ optional `mint`); `no_std`, no `alloc`, `forbid(unsafe_code)`; MSRV rule; determinism through `libm`. |
| [`0008`](./0008-the-name-helicoid.md) | the name `helicoid`; crate layout; crate boundaries. |
| [`0009`](./0009-what-helicoid-does-not-own.md) | the ownership test and what it excludes: sensors, solvers, preintegration, storage, `tf2` quirks. |
| [`0010`](./0010-seeded-from-tf-tree-math.md) | SO(3)/SE(3)/ScLERP ported from `tf_tree_math`; `tf_tree_math` is oracle #1 and must be dominated before `tf_tree` migrates. |
| [`0011`](./0011-continuous-time-waits-for-a-consumer.md) | `helicoid-spline` is gated on a consumer record; the core owes only `Ad`, `ad`, `Jr`, `Jr⁻¹`. |
| [`0012`](./0012-a-retraction-is-a-chart.md) | retractions are chart types; SE(3) has `Screw`, `Decoupled`, `WorldTranslation`; S² charts are frozen. |
| [`0013`](./0013-simd-lanes-owe-a-measurement.md) | SIMD `Real` lanes: the question and the bar; nothing decided. |
| [`0014`](./0014-the-conformance-metrics-open-readings.md) | conformance metric specifications, threshold sweep grids, and oracle protocols; nothing decided. |
| [`0015`](./0015-specification-gaps-found-while-building-the-instrument.md) | specification gaps found during instrument construction; nothing decided. |
| [`0016`](./0016-f32-exact-strata-for-the-scalar-coefficient-ids.md) | `f32` is swept and scored on `@f32` strata with exact binary32 inputs. |
| [`0017`](./0017-cbrt-and-mask-valued-roots.md) | `Real::cbrt`; `solve_cubic` returns roots with an `S::Mask` validity array. |
| [`0018`](./0018-libm-arch-is-bit-identical-for-exact-operations.md) | `libm` `arch` feature: exactly-rounded operations guarantee bit-identical results. |
| [`0019`](./0019-a-cholesky-solve-without-the-transpose.md) | `chol_solve(&L, b)` solves from factor with column reads, bit-identical to transposed copy. |
| [`0020`](./0020-dual-sqrt-at-zero-keeps-its-derivative.md) | `Dual::sqrt` at 0 retains `d / (2 sqrt v)` to report singularity. |
| [`0021`](./0021-the-cholesky-finite-guard-stays.md) | `chol` per-entry finite guard retained after benchmarking variants. |
| [`0022`](./0022-real-owes-acos-and-cos.md) | `Real` gains `acos` and `cos` via `libm` for improved performance and accuracy in `solve_cubic` and `eig3`. |
| [`0023`](./0023-eig3-departs-from-omnisac-and-its-limits.md) | `eig3` departures from baseline (orthonormal frame construction, non-finite handling); nothing decided. |
| [`0025`](./0025-a-structured-jacobian-and-a-sealed-side.md) | structured Jacobian traits: Rⁿ Jacobian structure, sealed `Side`, and threaded `dot_acc`. |
| [`0026`](./0026-write-dense-pays-for-zeros-not-for-checks.md) | `Jac::write_dense` zeroing costs vs bounds checks; decision left to PRs. |
| [`0027`](./0027-a-normalizing-constructor-normalizes.md) | `*_normalized` divides by norm; `renormalize` retains Newton step for drift repair. |
| [`0028`](./0028-two-ready-specs-disagree-on-the-blocks-of-sen3jac.md) | `SEn3Jac` internal block representation made `pub(crate)`. |
| [`0029`](./0029-a-product-needs-a-way-to-be-built-from-its-factors.md) | `Product::from_parts`/`parts` constructor added alongside factor types. |
| [`0030`](./0030-a-consumer-comparison-is-a-trigger-0021-does-not-list.md) | **superseded by `0080`.** Cholesky against external subjects. `0080` answered its questions 2 and 5, and moved 3 and 4 to its *Further work*. |
| [`0031`](./0031-what-the-cubic-port-inherits-from-omnisac.md) | **ready.** `solve_cubic` pairs its cube roots (and takes a quotient where `p > 0`), homogenises the monic cubic by powers of two, and drops the leading-coefficient floor; `1/a` stays (L5 measured against). Dominates nalgebra on every binary64 stratum. |
| [`0032`](./0032-domination-charges-helicoid-for-d16.md) | attribution of `so3_log` domination failures to `libm::atan2` vs host glibc; nothing decided. |
| [`0033`](./0033-a-latency-floor-is-measured-beside-the-comparison.md) | measuring benchmark latency floors beside comparison rather than stored; nothing decided. |
| [`0034`](./0034-a-second-transcendental-backend-owes-a-measurement.md) | **ready.** D16 is a property: bit-identical across targets for a given resolved `libm`, through one private kernel (lint-enforced), its three prohibitions literal. A bit change is a minor release from 1.0. A second backend needs a consumer benchmark that clears tf_tree `0016`'s bar end to end. |
| [`0035`](./0035-validation-has-three-layers-and-one-needs-a-quiet-machine.md) | multi-layered validation protocol for benchmark gating; nothing decided. |
| [`0036`](./0036-the-silence-in-5-1-costs-twenty-domination-failures.md) | analysis of `sen3_exp_n1` domination failures across Jacobian evaluation forms; nothing decided. |
| [`0037`](./0037-six-records-cannot-answer-a-one-variable-question.md) | evaluation of host-std kernel swaps and stratum sample counts; nothing decided. |
| [`0038`](./0038-a-program-comparison-is-not-a-bar.md) | paired pooled comparisons for algorithm selection vs domination max gating; nothing decided. |
| [`0039`](./0039-the-sweeps-grid-stops-below-its-own-optimum.md) | extending coefficient sweep grid across full domain and adding second series arm. Settled (`0040`). |
| [`0040`](./0040-a-draft-is-not-a-parking-space.md) | clarification of `ready` vs `draft` criteria to unblock settled records. |
| [`0042`](./0042-a-general-sen3-needs-a-way-to-be-built-from-its-parts.md) | **ready.** General `SEn3::from_parts`/`parts` constructor preserving unit invariants. |
| [`0041`](./0041-the-integration-is-an-adapter-at-tf-tree-math.md) | **ready.** `tf_tree` integration as an adapter layer at `tf_tree_math` behind a cargo feature. |
| [`0043`](./0043-the-geodesic-jacobian-ships-the-cancellation-free-form.md) | **ready.** Geodesic Jacobian uses cancellation-free form `(1-t) J_l((1-t)d) J_l⁻¹(d)`; adds `scale(k)`. |
| [`0044`](./0044-four-primitives-the-first-consumer-names-and-no-spec-does.md) | **ready.** `SEn3::renormalize` and `mul_inv` primitives for consumer composition. |
| [`0045`](./0045-two-phase-4-checks-cannot-be-taken-as-written.md) | **ready.** Uses geometric quaternion `Log` for Phase 4 geodesic references near $\pi$. |
| [`0046`](./0046-explained-by-record-needs-a-record-to-point-at.md) | **ready.** Domination exceptions registry in `conformance/baseline/exceptions.toml`. |
| [`0047`](./0047-the-second-arm-is-admitted-by-agreement-not-by-the-objective.md) | **ready.** Strict bit-agreement criterion for admitting second series arms below switch. |
| [`0048`](./0048-the-relative-transform-pair-earns-the-surface-dot-and-norm-do-not.md) | **ready.** Removed redundant `Quat::{dot, norm}`; added `SEn3::{mul_inv, inv_mul}` and scale-relative errors. |
| [`0049`](./0049-the-boundary-is-what-removes-a-way-to-be-wrong.md) | **ready.** Consumer boundary principle: retain operations that remove domain/order error modes. |
| [`0050`](./0050-the-geodesic-s-denominator-is-the-whole-domination-gap.md) | **ready.** `SO3::geodesic` blend denominator recomputed as $\sin(\alpha)$ to eliminate domination gap. |
| [`0051`](./0051-two-arms-on-the-switch-the-sweep-already-chose.md) | **ready.** Two-arm geodesic dispatch using `log_ratio` short arm to preserve near-identity performance. |
| [`0052`](./0052-real-owes-sin-and-the-corpus-does-not-move.md) | **ready.** `Real::sin` added, bit-identical to `sin_cos().0`, avoiding unused cosine calculations. |
| [`0053`](./0053-acos-is-better-the-roots-are-not-necessarily.md) | **ready.** Implements `Real::acos` and `cos` via `libm`; verified bounds across cubic and eigen solvers. |
| [`0054`](./0054-the-screw-twin-takes-two-arms-in-the-world-frame.md) | **ready.** `SE3::geodesic`'s screw twin: `SO3::geodesic`'s rotation, the world-frame translation in two arms; `se3_geodesic` dominates `ScLerp`, 4.6–6.1% faster inside `tf_tree`'s lookup bench. |
| [`0055`](./0055-the-branch-is-inlined-and-the-sign-is-a-bit.md) | **ready.** `Real::branch` is `#[inline]` and `abs`/`copysign` are `core`'s: coefficient kernels 0.55×, groups 0.915× geomean, bit-identical. |
| [`0056`](./0056-the-routines-d7-does-not-reach.md) | **ready.** Corpus ids for `solve_cubic`, `eig3`, `chol`/`chol_solve`, `quat_renormalize` and `real_*` at both precisions; root-set, gap-weighted eigenvector and mask metrics; a nalgebra oracle; committed rows reproduced bit for bit. |
| [`0057`](./0057-eig3-anchors-the-isolated-end.md) | **ready.** `eig3` anchors its frame on the more isolated end of the spectrum (`0023` question 2): top-pair ties and rank-1 go from `1e16 u` to `1e8`, 0.95× latency at `f64`. |
| [`0058`](./0058-a-drifted-quaternion-is-carried-not-vouched-for.md) | **ready.** Two entry points named: the *vouched* `from_wxyz_unchecked` (`2^-40`, asserted, NaN fails) and the *carried* struct literal into `SO3::from_quat_unchecked` (the drift band `2^-26.29`, nothing asserted). §12 states each operation's first-order error on a carried `q`, and §3.6 states the drift budget. No code change.
| [`0059`](./0059-the-short-arm-is-written-out.md) | **ready.** `SO3::geodesic`'s arm below the switch is written out by two identities exact in `R`, so `r`'s and `Exp`'s series stop waiting on its products; `horner` drops its `0 · z` step. `tf_tree`'s delegated `lerpslerp` goes from +15 % to parity, and `geo:consecutive` from 1.572 to 0.993 `u`.
| [`0060`](./0060-the-charts-go-first-and-name-their-frame.md) | **ready.** `PHASE5.md` §1, §2, §4, §5 proceed beside Phase 4 and SO(2)/SE(2), charts first. Corrects `0012`'s "agree to first order": `WorldTranslation` agrees after `diag(I, R)` (CH.5(c)). Chart types carry `S`; one `Manifold` impl per group (a blanket is E0119); `Lifted<C>` for `WithChart`; `TwistBlockJac` instead of a second `Jac` impl on `ProductJac` (E0283); `SE3::chart_transition`; six `se3_*_{retract,local}` ids.
| [`0061`](./0061-mat2-keeps-its-adjugate.md) | **ready.** `Mat2::inverse_adj -> (Mat2, S)` with `Mat3`'s contract (the adjugate is exact; `det` is the only rounding); id `mat2_inverse_adj` by conditioning, an exactly singular stratum, nalgebra as oracle. Answers `0015` (draft) PH.2.
| [`0062`](./0062-sin-theta-over-theta-is-public-and-differentiates-its-branch.md) | **ready.** `helicoid::sinc(θ²) -> (α, dα/dz)`: the value `2k cos(θ/2)` from `exp_coeffs`, the derivative by `Dual` through the same branch (measured 2.3× better than `(b − a)/2` on `[1, π]`); no new switch; id `coeff_alpha`.
| [`0063`](./0063-the-sphere-reads-its-sign-by-comparison-and-is-held-unit.md) | **ready.** S²: `ς = +1` unless `n_z < 0`; `S2` vouched, normalized and renormalized entry points with the quaternion's bounds (band measured on 3-vectors); `α/s = ½ r(s², w)`; `s2_local`'s reference without `mp.logm`; §11's direction metric.
| [`0064`](./0064-the-integrated-exponentials-reuse-the-swept-switches.md) | **ready.** `so3::{gamma1, gamma2, gamma_apply_jacobian::<M>}` (`M` 1 or 2): `gamma2_coeffs` is a group over the swept switches of `b` and `d` (no new switch, `q_coeffs`' to the bit); the Jacobian is `Dual<S, 3>` through the vector form, twin `reference::gamma_apply` (40 terms, `θ ≤ π`); `so3_gamma2` at both precisions (worst 1.46 / 1.68 `u`); the seeded defect Γ₂ from Γ₁'s coefficients.
| [`0065`](./0065-a-gaussian-names-its-side-and-is-stored-symmetric.md) | **ready.** `Gaussian<S, G, Sd, D>`: `new` the only constructor, `cov` read by its lower triangle and returned exactly symmetric; `to_left`/`to_right` (the adjoint of the computed inverse); `propagate`'s side contract documented; `mahalanobis_sq` returns `chol`'s mask; the round trip within `γ_{4D+42}(|B||A|)|Σ|(|B||A|)ᵀ` (worst 0.27); ids `gaussian_mahalanobis_{se3,se23}`, every maximum within 0.44 of its componentwise condition number.
| [`0066`](./0066-gamma-shares-its-coefficients-and-a-gaussian-factors-once.md) | **ready.** Corrects `0064`: the switches are on `θ²`, so Γ's exact arm runs from `θ = 2.79`–`3.11`. `so3::gamma_apply_jacobians -> (Γ_M v, ∂/∂φ, Γ_M)` from one `Dual<S, 1>` coefficient evaluation and the chain rule (0.46–0.52× the time); `HatSq` assembles `J_r`, `J_r⁻¹`, `Γ₂` entrywise (bit-identical bar zero signs); `Gaussian::cov` private, `Gaussian::whitener` / `Whitener` for many-point gating; `reference::gamma_apply` bounds `m`; the lower-triangle sandwich refused by count.
| [`0067`](./0067-the-sphere-takes-plane-coordinates-and-both-precisions.md) | **ready.** S²'s implementation readings: `RnTangent<S, 2>` coordinates with `Mat2<S>: Jac`, `S2: Manifold`; `renormalize(&mut self)`; `retract` from one `exp_coeffs`, `retract_jacobian` in the frame `retract` returns, `local_jacobian` through `Exp` of the rotation vector (never `1/(1 + n·m)`); `local(-n) = 0` by a safe argument and §12's far-hemisphere bound; corpus `s2_retract`/`s2_local` at both precisions, two plantings; `Vector::dot` by index (`f32` `local` 6x faster, bit-identical).
| [`0068`](./0068-real-owes-exp-and-exp-m1.md) | **ready.** `Real::{exp, exp_m1}` through `libm` (D16); `Dual::exp_m1`'s derivative from `exp`, not `1 + exp_m1`; corpus ids `real_exp`, `real_expm1` at both precisions; `PHASE2.md` §2's listing reconciled with `0022`. Sim(3)'s first prerequisite.
| [`0069`](./0069-the-sphere-review-and-the-matrix-product.md) | **ready.** The S² review: edge tests near the antipode and on the equator; no fused coefficient groups (below `θ ≈ 3.1` both are Horner); `Matrix × Matrix` as a loop — the nested `array::from_fn` built the output through the stack, so `jac_mul` takes 0.09–0.11 of its time and SE_N(3) `jr_inv` 0.40–0.60, bit-identical; one S² answer for the clean and planted subjects.
| [`0070`](./0070-sim3-a-structured-jacobian-and-a-series-on-delta.md) | **ready.** `NUMERICS.md` §9 stated from `docs/maths/sim3.md`: `Sim3`, `Sim3Tangent`, the structured `Sim3Jac` (amends `0005`); `Exp`, `Log` (no `mp.logm`), `Ad`, `ad`, the Jacobians by blocks; the eight coefficients as one group — limits at `Δ = 0`, the total-degree series below `r_×²`, closed forms above, scaled for `σ > 0`; switches generated by a ray sweep on the blocks' error; §12's rows to `ln MAX`; eight corpus ids at both precisions; a four-PR plan.
| [`0071`](./0071-the-corpus-budget-is-64-mib.md) | **ready.** Corpus v1 stays under 64 MiB (was 50 MB, reached by `sim3_ref.jsonl`), and Sim(3)'s eight ids take four records per cell (17 MB at six); a family over budget still reduces its sample count. |
| [`0072`](./0072-a-container-oracle-is-an-image-on-the-same-protocol.md) | **ready.** Container oracles: `Launch::Container` builds `docker/oracles/<name>` and runs the image by id offline, read-only, as the answers' owner, with same-path mounts so the file protocol is unchanged. The Sophus image (snapshot-pinned toolchain, Sophus 1.24.6 by commit and checksum, no contraction, no SIMD) answers 26 ids at binary64, a refused record `nan`. Sim(3) gets its first oracle: 87 of 434 strata not dominated. |
| [`0073`](./0073-sim3-two-sides-share-one-coefficient-evaluation.md) | **ready.** The Sim(3) coefficient group at `σ` and `-σ` from one branch, each to the bit (`Δ` and the transcendentals are even in `σ`, the series' recurrence is homogeneous); `⊖`'s Jacobians from one pair and one SO(3) pair, `⊕`'s from one evaluation; the §2.3 rows law stays exact. |
| [`0074`](./0074-the-gtsam-oracle-is-built-from-source-under-our-flags.md) | **ready.** The GTSAM container on `0072`'s launcher: GTSAM 4.3.0 by commit and checksum, built from source in the image with the runner under one set of flags (no contraction, no SIMD), not a wheel. 47 ids at binary64 — SO(2), SE(2), SO(3), `Γ₂`, SE_N(3) for `N` = 1–3 (`ExtendedPose3`), the screw chart, Sim(3) `exp`/`log`/`ad`/`act` — right Jacobians from `ExpmapDerivative`, left at `−τ`; `Similarity3`'s `x = s·t`. 919 more strata paired; domination failures 301 → 634; Sim(3) fails 151 of 485. |
| [`0075`](./0075-the-manif-oracle-reads-se23-as-position-then-velocity.md) | **ready.** The manif container on `0072`'s launcher: manif 0.0.5 by commit and checksum on the Sophus image's toolchain, its unit-norm checks kept as refusals. 36 ids at binary64, every Jacobian by manif's own formula; SE₂(3) `[ρ_p; θ; ν]` maps to `[φ; ρ_v; ρ_p]`; the screw and world charts are `SE3` and `Bundle<SO3, R3>` `rplus`/`rminus`. SE₂(3) and the charts get their first oracle: 146 of 576 strata not dominated. |
| [`0076`](./0076-the-oracles-first-verdict-r-of-q-stays-and-drift-is-excepted.md) | **ready.** The container oracles' verdict: `R(q)` keeps its homogeneous form (Eigen's unit form wins `θ ≤ 1` and loses `θ ≥ 3`; the gap is the rounded input's `η`, which only normalizing removes); `so3_act` and `sen3_ad_n1` at `q:nonunit` excepted as `0058`'s carried drift. |
| [`0077`](./0077-an-oracle-answers-at-the-precision-it-is-asked.md) | **ready.** The oracle protocol gains `--precision f32`: the `@f32` records alone, inputs narrowed losslessly, computed and written in binary32; `--version` unchanged; `Runner::answers_f32` declares the binary32 ids, and a binary64-only runner (GTSAM, manif, `tf_tree_math`, sophus-rs) is not launched at `f32` and passes. nalgebra (all eight ids) and Sophus (`Sim3<float>`) answer; the envelope reads both precisions. 489 binary32 strata paired; failures 707 → 805; Sim(3) fails 74 of 428 at `f32`. |
| [`0078`](./0078-a-comparison-drops-ties-and-prints-both-maxima.md) | **ready.** `cargo xtask compare`'s readings of `0038`: ties dropped, a lone zero a direction and not a ratio, non-finite counted and compared with nothing; both "maximum ratios" printed (per-record extremes and `max_b/max_a` with a strata tally); per id, then pooled beside. The twin `helicoid:unit-r` reproduces `0076`'s tally exactly; the task reproduces `0038`'s 21 / 104. |
| [`0079`](./0079-from-matrix-normalizes-a-column-of-k.md) | **ready.** `SO3::from_matrix` normalizes a column of `K + I` (SO.12(b)): one `sqrt`, 0.54–0.60 of the time at `f64`, 0.70 at `f32`; overall max error down, mean up 0.03–0.05 `u`, 11 strata's max up to +0.32 `u`, a trade the owner accepted. |
| [`0080`](./0080-chol-unrolls-what-llvm-will-not.md) | **ready.** `chol` runs fully unrolled for `4 <= N <= 9` by two `macro_rules!`, the loop's operations in its order, so `(L, mask)` is bit-identical (two proptests over every special, all lanes); the loop stays for `N <= 3` (LLVM unrolls and unswitches it) and `N > 9`. Meets `0021`'s bar under the gate: `f64` `N = 6` 0.588×, `f32` 0.726×, `Dual` 0.743×, failing input 0.63–0.83×. Answers `0030` question 2. |
| [`0081`](./0081-eig3-is-four-jacobi-sweeps.md) | **ready.** `eig3` is four cyclic Jacobi sweeps (`NUMERICS.md` §17) with an exact power-of-two scaling, replacing Smith's closed form, whose `acos` lost up to `sqrt(u) |A|` at a close pair: every `eig3` stratum dominates nalgebra at both precisions (max 4.76 `u` `f64`, 3.24 `u` `f32`), every finite matrix is in the domain, at 2.1–2.8× the latency. Answers `0023` (draft) questions 1 and 4. |
| [`0082`](./0082-the-sphere-chart-is-oriented-and-shows-its-frames.md) | **ready.** S²'s frame is oriented (`b₂ = ς He_y`, `b₁ × b₂ = n` on both hemispheres) and public: `S2Chart::rotation_basis` (`B`) and `displacement_basis` (`T = [-b₂ b₁] = ∂retract/∂δ` at `0`). |
| [`0083`](./0083-eig3-is-also-judged-against-eigen.md) | **ready.** The Sophus image answers `eig3` with its Eigen 3.4.0 `SelfAdjointEigenSolver<Matrix3<T>>::compute`, at both precisions: a second oracle on all 50 strata, which `helicoid` dominates (ratio of maxima 0.19–0.59); `computeDirect` measured once and not a subject. |
| [`0084`](./0084-eig3-leaves-jacobi-when-it-has-converged.md) | **ready.** `eig3` leaves its Jacobi loop when both entries left are at most `u |A|` (deterministic, D16), then one polishing step keeps `Dual` lanes first-order; a rotation whose next entry reaches the next rotation squared, a square root earlier. 0.21–0.72× `#154`'s latency, faster than Eigen on a random spectrum; corpus max 3.06 / 2.67 `u`, dominating nalgebra and Eigen; overall and frame maxima below `#154`'s; per-family figures recorded, not a bar. |
| [`0085`](./0085-from-matrix-is-measured-on-exact-rotations.md) | **ready.** `so3_from_matrix` is scored forward by §11's direction metric in ℝ⁴, sign-free, against the polar factor; new stratum `matrix:exact` (the 24 signed permutation rotations, from integers) and `@f32` twins of every stratum; `Rule::BackwardOnly` retired. The `0079` kernel is dominated on 18 + 14 strata, `matrix:exact` at `f64` among them (0.7985 against 0.7071); the kernel fix is owed. |
| [`0086`](./0086-from-matrix-keeps-shepperd-s-quotients-on-the-sphere.md) | **ready.** `SO3::from_matrix` returns Shepperd's quotients in `tf_tree_math`'s arithmetic (each axis `c_j` adds its own diagonal first) when they are within `16u` of the sphere, twice the 8u any rounded rotation reaches, and `0079`'s normalized column otherwise. Exact rotations are unit (a quarter turn: `fl(‖q‖²) = 1`, was `1 − 2u`). Strata above the best oracle: 8 of 31 at `f64` (was 18), 9 at `f32` (was 14); `f64` equals `tf_tree_math` but on `q:nonunit`. `f32` costs 1.63–1.85×, accepted. Open: the `π` strata under `0088`. |
| [`0088`](./0088-the-domination-bar-must-tell-worse-from-equal.md) | **draft.** A strict max-over-max bar fails a program as good as its oracle on about half its strata, at any `n` (expected 1098 of 2656; Sophus fails it against GTSAM on 47%). Proposes: domination fails only when a paired swap test on `max c − max o` resolves the loss at `α = 0.001`; strata under 10 records (78% of the paired ones) read a git-ignored 64-record deep tier. Today that leaves 31 resolved failures, the fix queue (12 `so3_from_matrix`), out of 794. Owner's call: it reopens `0038` decision 1. |

## Lifecycle

Four statuses on one document type; no folder moves.

- **draft**: open questions present. **A draft authorises nothing**: no spec §0.0 row, amendment
  banner, code comment or status table may cite one as settled (`just lint` checks this).
- **draft (awaiting a decision)**: the evidence is in and a choice is owed. A record whose
  *Decision* states options and a recommendation is in this state, not the one above, and says so
  (`0040` item 3).
- **ready**: every open question the *Decision* depends on is resolved and the implementation plan
  is concrete. **A question the Decision does not depend on does not block**: it goes in *Further
  work*, which a `ready` record may have and an *Open questions* section may not duplicate
  (`0040` items 1 and 2). A measurement whose protocol another session can rerun is citable as
  measured from `ready` on, whatever the record proposes (item 4).
- **implemented**: code shipped, PRs linked, document frozen. A record is not edited to match the
  code; a new record supersedes it. The only edits are the status line and an amendment banner
  signed and dated by a later record.
- **superseded by NNNN**: replaced; the document stays in place.

A `ready` record is the contract: the *Decision* is implemented as stated, each plan step lands as
one PR in order with its listed verification, and an open question found mid-implementation means
stop and ask.

Numbers are sequential four-digit and never renumbered. Filenames are
`NNNN-kebab-case-noun-phrase.md`.
