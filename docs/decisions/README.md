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
| [`0030`](./0030-a-consumer-comparison-is-a-trigger-0021-does-not-list.md) | Cholesky benchmarking against external subjects; nothing decided. |
| [`0031`](./0031-what-the-cubic-port-inherits-from-omnisac.md) | **ready.** `solve_cubic` pairs its cube roots (and takes a quotient where `p > 0`), homogenises the monic cubic by powers of two, and drops the leading-coefficient floor; `1/a` stays (L5 measured against). Dominates nalgebra on every binary64 stratum. |
| [`0032`](./0032-domination-charges-helicoid-for-d16.md) | attribution of `so3_log` domination failures to `libm::atan2` vs host glibc; nothing decided. |
| [`0033`](./0033-a-latency-floor-is-measured-beside-the-comparison.md) | measuring benchmark latency floors beside comparison rather than stored; nothing decided. |
| [`0034`](./0034-a-second-transcendental-backend-owes-a-measurement.md) | determinism bounds across `libm` versions and criteria for alternate backends; nothing decided. |
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
