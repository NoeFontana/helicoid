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
| [`0008`](./0008-the-name-helicoid.md) | the name `helicoid`; crate layout; the `helicoid` / `locus-*` boundary. |
| [`0009`](./0009-what-helicoid-does-not-own.md) | the ownership test and what it excludes: sensors, solvers, preintegration, storage, `tf2` quirks. |
| [`0010`](./0010-seeded-from-tf-tree-math.md) | SO(3)/SE(3)/ScLERP ported from `tf_tree_math`; `tf_tree_math` is oracle #1 and must be dominated before `tf_tree` migrates. |
| [`0011`](./0011-continuous-time-waits-for-a-consumer.md) | `helicoid-spline` is gated on a consumer record; the core owes only `Ad`, `ad`, `Jr`, `Jr⁻¹`. |
| [`0012`](./0012-a-retraction-is-a-chart.md) | retractions are chart types; SE(3) has `Screw`, `Decoupled`, `WorldTranslation`; S² charts are frozen. |
| [`0013`](./0013-simd-lanes-owe-a-measurement.md) | SIMD `Real` lanes: the question and the bar; nothing decided. |
| [`0014`](./0014-the-conformance-metrics-open-readings.md) | the readings the conformance metric takes where `NUMERICS.md` §11 is silent, the field of `PHASE1.md` §10's `b` curve, and the threshold sweep's grid, top, prior and ties; nothing decided. |
| [`0015`](./0015-specification-gaps-found-while-building-the-instrument.md) | the gaps the generator, `helicoid-linalg`, the lint and the derivation pages found in the specs: passage, evidence, options and a recommendation each; nothing decided. |
| [`0018`](./0018-libm-arch-is-bit-identical-for-exact-operations.md) | `libm` `arch` on: the routed `sqrt`/`fma`/`rint` are exactly rounded, so every non-NaN output is bit-identical; per-target dispatch read; NaN bits stay outside D16; a pinned-digest `sqrt` test. |
| [`0019`](./0019-a-cholesky-solve-without-the-transpose.md) | `chol_solve(&L, b)`: `A x = b` from the factor, `L^T` read by column, bit-identical to the transposed-copy composition; the column solve stays private. |
| [`0020`](./0020-dual-sqrt-at-zero-keeps-its-derivative.md) | `Dual::sqrt` at 0 keeps `d / (2 sqrt v)`: the NaN is the report, the §10 row needs it; a zero-safe norm is a later, additive record. |
| [`0021`](./0021-the-cholesky-finite-guard-stays.md) | the `chol` per-entry finite guard stays: four bit-identical variants measured, none worth a rewrite; the 15% bar and the finiteness of `L` are stated. |
| [`0025`](./0025-a-structured-jacobian-and-a-sealed-side.md) | the trait layer as built: Rⁿ's Jacobian is structured, correcting `PHASE3.md` §7 rather than `0005`; `Side` is sealed while the side selector stays with the SO(3) PR; `Tangent` gains a threaded `dot_acc` so `Product` can meet the bound-`0` dense-order law; `read_dense` poisons a missing entry; a group exposes an `Add`-carrying field only where it is abelian. |
| [`0027`](./0027-a-normalizing-constructor-normalizes.md) | `*_normalized` divides by the norm and `renormalize` keeps the Newton step, settling `API.md` R6 against `NUMERICS.md` §3.6; the step's domain contains `from_wxyz_unchecked`'s, so it stays drift repair. |

## Lifecycle

Four statuses on one document type; no folder moves.

- **draft**: open questions present. **A draft authorises nothing**: no spec §0.0 row, amendment
  banner, code comment or status table may cite one as settled (`just lint` checks this).
- **ready**: every open question the *Decision* depends on is resolved and the implementation plan
  is concrete.
- **implemented**: code shipped, PRs linked, document frozen. A record is not edited to match the
  code; a new record supersedes it. The only edits are the status line and an amendment banner
  signed and dated by a later record.
- **superseded by NNNN**: replaced; the document stays in place.

A `ready` record is the contract: the *Decision* is implemented as stated, each plan step lands as
one PR in order with its listed verification, and an open question found mid-implementation means
stop and ask.

Numbers are sequential four-digit and never renumbered. Filenames are
`NNNN-kebab-case-noun-phrase.md`.
