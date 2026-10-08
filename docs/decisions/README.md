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
