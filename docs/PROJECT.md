# helicoid — Project Overview

> **Read this before any phase spec.** This document explains *what* and *why*; `NUMERICS.md`
> states *which formula*; the phase specs explain *how*. When a spec does not answer a question,
> consult the decision log in §5.

## 1. What this is

`helicoid` is the Lie-group layer of one robotics perception stack: SO(2), SO(3), SE(2),
SE_N(3) (SE(3), SE₂(3), …), Sim(3), S², their products, and the Jacobians, charts and geodesics
that estimation needs. It is `no_std`, allocation-free, `unsafe`-free, generic over the scalar
(`f64`, `f32`, forward-mode dual numbers), and bit-identical across x86_64, aarch64 and wasm32.

It exists because the stack implements the same mathematics several times, with two tangent orders
and three retractions:

| Where | Representation | Tangent order | Retraction |
|---|---|---|---|
| `tf_tree_math` (published) | Hamilton, `w` first, `Iso3 { q, t }` | `[ω, v]` | `T̂ · exp(ξ^)` (SE(3) exponential) |
| locus-tag `Pose::retract` | nalgebra `Rotation3` + `Vector3` | `[v, ω]` (`delta[0..3]` is translation) | `(R·Exp(ω), R·v + t)` — decoupled, not SE(3) Exp |
| omnisac refinement (PnP, E, rigid) | private | private | private |
| locus_fusion | Sophus (Eigen) | `[υ; ω]` | Sophus right-⊕ |
| locus-calib (planned) | — | — | would be a fifth copy |

`helicoid` replaces all of them, and is the only place a small-angle series is allowed to exist.

**Non-goals.** No sensor models (`locus-camera`, `locus-imu`, `locus-lidar`), no solver, no IMU
preintegration, no storage formats, no splines until a consumer opens the gate
([`0009`](./decisions/0009-what-helicoid-does-not-own.md),
[`0011`](./decisions/0011-continuous-time-waits-for-a-consumer.md)). `helicoid` computes;
consumers store, estimate and decide.

## 2. The problems being solved

| Problem | Consequence | Our answer |
|---|---|---|
| Lie code reimplemented per consumer | Divergent conventions (`[ω, v]` vs `[v, ω]`, `w`-first vs `xyzw`); a covariance handed between libraries is silently wrong by a permutation or an `Ad` | One layer; conventions in one record ([`0002`](./decisions/0002-one-convention-for-a-stack-that-already-disagrees.md)); `Gaussian<G, Side>` |
| Small-angle switch asserted (`1e-8` in most libraries) | `(θ − sin θ)/θ³` keeps ~7 digits at `θ = 1e-4`; Barfoot's `Q` coefficient keeps ~5 at `1e-2` | Switch points and series lengths generated per coefficient and precision from a committed sweep ([`0004`](./decisions/0004-switch-points-are-generated-not-typed.md)) |
| `log` via `acos((tr R − 1)/2)` | ~8 digits lost near 0 and near π | Quaternion `atan2` form only (D5) |
| Jacobians checked by finite differences | FD cannot see a near-identity error below its own step error | Forward-mode `Dual` through the shipped code; an mpmath corpus computed from definitions ([`0006`](./decisions/0006-the-instrument-comes-first.md)) |
| Dense 6×6 Jacobian algebra | 216 multiplications per product; one more closed form (`Jr⁻¹`) to get wrong | Dual-matrix structure: 81 multiplications; `Jr⁻¹` derived by algebra ([`0005`](./decisions/0005-the-jacobian-is-a-dual-matrix.md)) |
| No SE₂(3), S², Sim(3), Γ functions in the stack | VIO and calibration re-derive them | Phases 3 and 5 |
| One retraction baked into each LM | Changing it silently changes iterates, RunRecords and bench tails | Charts are named types ([`0012`](./decisions/0012-a-retraction-is-a-chart.md)) |
| Branchy numeric code | Cannot run on SIMD lanes; unsafe under reverse-mode AD | `Real` has no `PartialOrd`; masks and `branch` ([`0003`](./decisions/0003-the-scalar-that-cannot-say-less-than.md)) |
| Platform `libm` differences | Native and wasm32 outputs differ by an ulp; Fuse's log-rebuildable stores cannot be bit-proven | `libm` crate on every target (D16) |
| "State of the art" asserted | No evidence which strata a library is worse on | Domination bar over an oracle envelope, per stratum, on the max (D8) |

## 3. Architecture in one page

**Crates.** `helicoid-linalg` is the leaf: `Real`/`Mask`/`Blend`/`Dual`, fixed-size
`Vector`/`Matrix`/`Point`, strided views, `eig3`/`svd3`/`solve_cubic`, fixed-size Cholesky.
`helicoid` holds the coefficient kernel (`coeffs`, the only home of cancelling expressions), the
groups, sides, structured Jacobians, charts, geodesics, `Gaussian`, and the `reference` twins. Both
are `no_std`, no `alloc`, `forbid(unsafe_code)`, and depend on `libm` alone (plus optional `mint`).

**The verification spine** is built before any group (Phase 1) and never leaves:

1. `conformance/generate/` — Python + pinned mpmath at 120 digits, computing every quantity from
   its *definition* (matrix exponential and logarithm, defining series), never from the closed
   forms under test. Inputs are exact binary64; the corpus is committed and regenerates
   byte-identically.
2. `xtask conformance` — runs in-process subjects (helicoid, seeded defects) over the corpus;
   per-stratum max and p99 error in units of `u`.
3. `runners/` and `docker/oracles/` — `tf_tree_math`, sophus-rs, Sophus, manif, GTSAM through one
   file protocol.
4. `xtask envelope` — domination over the best oracle and no-regress against the committed
   baseline, per function, per stratum, per precision.
5. `Dual` — every Jacobian is compared with forward-mode differentiation *of the shipped code*,
   Taylor branches included.
6. `reference` — composite routines keep an obvious, slow twin, proptested against the fast one.
7. `xtask thresholds` — the sweep that generates every switch point and series length.

**Consumers** sit above: `tf_tree` (composition, geodesics, `Ad`), omnisac and locus-tag
(retractions, charts, `from_matrix`), locus-calib (charts, ambient Jacobians, `Gaussian`),
locus_fusion (SE₂(3), Γ, both sides), fuse-geometry (wasm32, determinism). The solver crate and the
`locus-*` sensor-model crates depend on `helicoid`; `helicoid` depends on none of them.

**The load-bearing consequence:** the corpus specifies *correct*, the envelope specifies *state of
the art*, and determinism (D16) makes the no-regress bar exact — equal-or-better, bit for bit, on
each stratum's max.

## 4. Roadmap

- **Phase 1 — skeleton and the instrument.** Workspace, lints, CI on four targets, the reference
  generator and corpus v1, the conformance harness, the threshold sweep, oracle runners, the
  envelope, the bench gate; all validated against planted defects. **No group code.**
  `docs/PHASE1.md`.
- **Phase 2 — `helicoid-linalg`.** The scalar model (`Real`, `Mask`, `Blend`), `Dual`, fixed-size
  types, strided views, `eig3`/`svd3`/`solve_cubic` migrated from omnisac. `docs/PHASE2.md`.
- **Phase 3 — core groups.** The coefficient kernel and its generated thresholds; SO(2), SO(3),
  SE(2), SE_N(3), Rⁿ, products; sides; every Jacobian; reference twins; the first envelope.
  `docs/PHASE3.md`.
- **Phase 4 — geodesics and the `tf_tree` migration.** Geodesic, its Jacobians, the dual-quaternion
  fast twin, invariance tests; the parity table, envelope and bench evidence a `tf_tree` record
  needs. `docs/PHASE4.md`.
- **Phase 5 — extended geometry and the retraction migrations.** Charts, S², Sim(3), Γ₁/Γ₂,
  `Gaussian`; omnisac, locus-tag and locus_fusion onto `helicoid`. `docs/PHASE5.md`.
- **Phase 6 — interop, determinism, locus-calib, 1.0.** Ambient (Ceres-style) Jacobians, `mint`,
  the cross-target bit-identity gate, locus-calib's first commit, the 1.0 criteria.
  `docs/PHASE6.md`.
- **Phase 7 — continuous time, gated.** `helicoid-spline`. `docs/PHASE7.md` is a **requirements
  artifact, not an implementation authorization**
  ([`0011`](./decisions/0011-continuous-time-waits-for-a-consumer.md)).

Status per phase: the §0.0 table heading each spec is authoritative.

## 5. Decision log

- **D1 — One Lie layer for the stack.** `helicoid` replaces `tf_tree_math`'s SE(3), omnisac's and
  locus-tag's retractions, and locus_fusion's Sophus usage where it is Rust; locus-calib uses it
  from its first commit. *Do not* add a second implementation of `Exp`, `Log` or a Jacobian anywhere
  in the stack; a consumer needing a variant gets a chart here (D12), not a local fork.
- **D2 — Values, not storage.** Every type is a `Copy` value. `#[repr(C)]` makes layout predictable
  but **layout is not a semver contract**: consumers own their storage formats (`tf_tree`'s arena
  records stay `tf_tree`'s, its D4). *Do not* add `Pod`, `Zeroable` or `serde` to a `helicoid`
  type.
- **D3 — One convention record.** Hamilton, `w` first, active, `a * b = T_a_x · T_x_b`,
  rotation-first tangents, right perturbation by default, sides always spelled
  ([`0002`](./decisions/0002-one-convention-for-a-stack-that-already-disagrees.md)). *Do not* add
  an `Add`/`Sub` impl for ⊕/⊖, a `From<[S; 6]>`, or an unlabeled dense order.
- **D4 — Generic scalar; comparisons return masks**
  ([`0003`](./decisions/0003-the-scalar-that-cannot-say-less-than.md)). `Real` has no
  `PartialOrd`. *Do not* add `PartialOrd`, `PartialEq` or a `bool`-returning comparison to `Real`.
- **D5 — Numerics are measured, not assumed** (inherits `tf_tree` D12). `Log` goes through the
  quaternion `atan2` form; cancellation-free rewrites are used where they exist (`NUMERICS.md` §4);
  switch points and series lengths are generated per coefficient and per precision
  ([`0004`](./decisions/0004-switch-points-are-generated-not-typed.md)). *Do not* hand-edit
  `coeffs/generated.rs`; *do not* "simplify" a switch to `1e-8`.
- **D6 — Reference plus fast, forever** (inherits `tf_tree` D13). Every row of `NUMERICS.md` §14
  keeps an obvious, slow twin in `helicoid::reference` and a `*_matches_reference` proptest.
  Primitives (`Exp`, `Log`, coefficients) are checked against the corpus instead.
- **D7 — Instrument before algorithm** ([`0006`](./decisions/0006-the-instrument-comes-first.md)).
  A public numeric routine without a corpus stratum is unverified and may not ship in a release.
  *Do not* delete or narrow a stratum because it fails.
- **D8 — Domination and no-regress, never means.** Per function, stratum and precision: max error
  ≤ the best oracle's max error, and ≤ the committed baseline. Performance is reported; only
  self-regression and consumer never-regress gates are enforced.
- **D9 — Jacobians are structured** ([`0005`](./decisions/0005-the-jacobian-is-a-dual-matrix.md)).
  SE_N(3)'s Jacobians are dual matrices closed under product and inverse; dense only through
  `write_dense` and `sandwich`.
- **D10 — Dependency and unsafe budget**
  ([`0007`](./decisions/0007-the-budget-a-foundation-can-afford.md)). `libm` (+ optional `mint`);
  `no_std`, no `alloc`, `forbid(unsafe_code)`. `xtask` and runners are unrestricted and never a
  normal dependency of a library crate.
- **D11 — No panic on valid input.** A routine's domain is a `# Domain` rustdoc section and a
  `debug_assert!`; release builds never check it. Out-of-domain release behaviour is unspecified
  but never a panic and never UB. The one documented panic class is out-of-bounds strided access (`get`, `set`, `block`).
- **D12 — A retraction is a chart** ([`0012`](./decisions/0012-a-retraction-is-a-chart.md)).
  SE(3) has three named charts (`Screw`, `Decoupled`, `WorldTranslation`); S² charts are frozen
  per linearization. *Do not* write a retraction inline in a consumer.
- **D13 — What `helicoid` does not own**
  ([`0009`](./decisions/0009-what-helicoid-does-not-own.md)). Ownership test: a feature belongs
  here iff it needs no knowledge of a sensor, a solver, or a consumer's state or storage layout.
- **D14 — Continuous time is gated**
  ([`0011`](./decisions/0011-continuous-time-waits-for-a-consumer.md)). `helicoid-spline` opens on
  a consumer's record, not on schedule. What splines will consume (`Ad`, `ad`, `Jr`, `Jr⁻¹`) ships
  in Phase 3 regardless.
- **D15 — Seeded from `tf_tree_math`, and measured against it**
  ([`0010`](./decisions/0010-seeded-from-tf-tree-math.md)). `tf_tree_math` is oracle #1;
  `helicoid` must dominate it on every paired stratum before `tf_tree` migrates.
- **D16 — Determinism by construction.** Every transcendental goes through the `libm` crate on
  every target; no `mul_add`, no `target-cpu`, no fast-math. Outputs are bit-identical across
  x86_64, aarch64 and wasm32 (NaN sign and payload apart); `just determinism` checks it. This is what
  makes D8's no-regress bar exact. `libm`'s `arch` feature is on: it routes only exactly rounded
  operations to the target's instruction ([`0018`](./decisions/0018-libm-arch-is-bit-identical-for-exact-operations.md)).
- **D17 — MSRV never exceeds a consumer's.** Today **1.87** (`tf_tree`). Raising it is a record
  that names every consumer's MSRV.
- **D18 — Apache-2.0 / MIT dual license**, matching `tf_tree` D20.

### 5.1 Gated-work ledger

**Nothing in this table is authorised.** A row is where the argument goes when the work is
scheduled; adding a row is not a decision, removing one is.

| Entry | Source | Gate |
|---|---|---|
| `helicoid-spline` | [`PHASE7.md`](./PHASE7.md), [`0011`](./decisions/0011-continuous-time-waits-for-a-consumer.md) | A consumer record: locus-calib camera–IMU, locus_fusion continuous time, or `tf_tree` Phase 6 |
| SIMD `Real` impl | [`0013`](./decisions/0013-simd-lanes-owe-a-measurement.md) (draft) | A consumer gate win that pays for its dependency cost; `tf_tree` 0016's evidence is the bar to beat |
| Closed-form Γ directional Jacobians | [`PHASE5.md`](./PHASE5.md) §4 | A bench showing the `Dual` path is a consumer bottleneck |
| `helicoid-py` | no record | locus-calib's Python tooling asks for batched NumPy access |
| The small dense LM (omnisac's LM core) | [`0009`](./decisions/0009-what-helicoid-does-not-own.md) | Not here: solver tier. Listed so the request has somewhere to be refused |

## 6. Design smells — stop if you catch yourself doing these

- Writing `if` on a float, adding `PartialOrd` to `Real`, or returning `bool` from a numeric
  comparison (D4)
- Typing a series coefficient or a switch point; evaluating `(θ − sin θ)/θ³`, `c`, `d` or `e`
  outside `coeffs`; `acos` near a rotation (D5)
- Returning `[[S; 6]; 6]` from a group method; forming a dense Jacobian to multiply two structured
  ones (D9)
- Adding a dependency, `unsafe`, `std`, `alloc` or a feature to a library crate (D10); `mul_add` or
  a `target-cpu` flag (D16)
- A camera, IMU or LiDAR model; a solver; preintegration; a spline (D13, D14)
- A retraction written inline in a consumer (D12); a covariance without its side
- `+`/`-` for ⊕/⊖; an unlabeled `[S; 6]` crossing an API (D3)
- Quoting a mean error, or one number aggregated over strata (D8)
- A fast routine without a reference twin (D6); deleting a stratum because it fails (D7)
- `Pod`, `serde`, or a layout test that a consumer's storage format depends on (D2)

## 7. Glossary

| Term | Meaning |
|---|---|
| Side | Which end a perturbation multiplies: `Right` (`X·Exp(τ)`, default) or `Left` (`Exp(τ)·X`). Each side's Jacobians are expressed in that side's convention. |
| Chart | A retraction and its local inverse, frozen at a base point. A type, not a flag ([`0012`](./decisions/0012-a-retraction-is-a-chart.md)). |
| Dense order | The order a tangent's components take in a flat array: rotation first (`NUMERICS.md` §1). Visible only through `write_dense`/`read_dense` and named converters. |
| Mask | The result of comparing two `Real`s: `bool` for scalars and duals, a lane mask for SIMD lanes. |
| Safe argument | The exact arm of a branch evaluates at `select(small, 1, θ)`, so no arm ever produces a non-finite value. |
| Coefficient | A scalar function of `θ²` from `NUMERICS.md` §4 with a closed form and a series. |
| Switch point | The branch-variable value below which a coefficient uses its series; generated, per precision. |
| Stratum | A named input family of the corpus (`theta:1e-8`, `theta:pi-1e-6`, `rho:1e4`, …). Every bar is per stratum. |
| Corpus | Committed JSONL of exact inputs and 30-digit references, computed from definitions by mpmath. |
| Subject / oracle | Anything evaluated over the corpus; an oracle is a subject that is not `helicoid`. |
| Envelope | Per (function, stratum, precision), the best oracle max error. |
| Baseline | The committed per-stratum max errors of `helicoid`; the no-regress bar. |
| `u` | Unit roundoff: `2⁻⁵³` for `f64` (and `Dual<f64, N>`), `2⁻²⁴` for `f32`. |
| Reference twin | The obvious, slow implementation a fast composite routine is proptested against (D6). |

## 8. Document map

`PHASE1.md`–`PHASE6.md` are the normative phase specs; `PHASE7.md` is **gated by
[`0011`](./decisions/0011-continuous-time-waits-for-a-consumer.md)**. `NUMERICS.md` is not a phase:
it is the formula contract every phase implements, read before touching a formula. `API.md` is not
a phase: the six rules (§1), the surfaces (§2–§4), stability (§5), the new-surface check (§6), read
before adding public API. `decisions/` holds the records; `CLAUDE.md` at the root is the agent
entry point.
