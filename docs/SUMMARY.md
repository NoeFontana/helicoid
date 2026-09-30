# Summary

[Overview](README.md)

# Start here
- [Project overview and roadmap](PROJECT.md)

# Reference
- [Numerics: formulas, series, domains](NUMERICS.md)
- [API contract](API.md)

# Derivations (non-normative)
- [Maths: notation, conventions, how results are checked](maths/index.md)
- [Lie groups, adjoints and the Jacobians of Exp](maths/lie-groups.md)
- [SO(3): quaternion, Exp, Log, Jacobians, from_matrix](maths/so3.md)
- [The coefficient catalogue: series, cancellation, switch points](maths/coefficients.md)
- [SE_N(3): Exp, Log, adjoints, the Q block, the dual-matrix algebra](maths/se3.md)
- [SO(2) and SE(2): Exp, Log, adjoints, the rotation-first Jacobians](maths/so2-se2.md)
- [Error analysis, conditioning, forward-mode AD, the reference corpus](maths/error-analysis.md)

# Specifications
- [Phase 1: skeleton and the instrument](PHASE1.md)
- [Phase 2: helicoid-linalg](PHASE2.md)
- [Phase 3: core groups](PHASE3.md)
- [Phase 4: geodesics and the tf_tree migration](PHASE4.md)
- [Phase 5: extended geometry and the retraction migrations](PHASE5.md)
- [Phase 6: interop, determinism, locus-calib, 1.0](PHASE6.md)
- [Phase 7: continuous time (gated)](PHASE7.md)

# Decisions
- [Index and lifecycle](decisions/README.md)
- [0001: Record architectural decisions in `docs/decisions/`](decisions/0001-record-architectural-decisions.md)
- [0002: One convention, for a stack that already disagrees](decisions/0002-one-convention-for-a-stack-that-already-disagrees.md)
- [0003: The scalar that cannot say "less than"](decisions/0003-the-scalar-that-cannot-say-less-than.md)
- [0004: Switch points are generated, not typed](decisions/0004-switch-points-are-generated-not-typed.md)
- [0005: The Jacobian is a dual matrix](decisions/0005-the-jacobian-is-a-dual-matrix.md)
- [0006: The instrument comes first](decisions/0006-the-instrument-comes-first.md)
- [0007: The budget a foundation can afford](decisions/0007-the-budget-a-foundation-can-afford.md)
- [0008: The name `helicoid`](decisions/0008-the-name-helicoid.md)
- [0009: What `helicoid` does not own](decisions/0009-what-helicoid-does-not-own.md)
- [0010: Seeded from `tf_tree_math`, and measured against it](decisions/0010-seeded-from-tf-tree-math.md)
- [0011: Continuous time waits for a consumer](decisions/0011-continuous-time-waits-for-a-consumer.md)
- [0012: A retraction is a chart](decisions/0012-a-retraction-is-a-chart.md)
- [0013: SIMD lanes owe a measurement (draft)](decisions/0013-simd-lanes-owe-a-measurement.md)
- [0018: libm arch is bit-identical for exact operations](decisions/0018-libm-arch-is-bit-identical-for-exact-operations.md)
- [0019: A Cholesky solve without the transpose](decisions/0019-a-cholesky-solve-without-the-transpose.md)
