# helicoid

Lie groups for robotics perception: SO(2), SO(3), SE(2), SE_N(3) (SE(3), SE₂(3), …), Sim(3), S²
and their products, with structured Jacobians, charts and geodesics — `no_std`, no allocation, no
`unsafe`, generic over the scalar (`f64`, `f32`, forward-mode dual numbers), bit-identical across
x86_64, aarch64 and wasm32.

Every switch point is generated from a committed sweep; every routine is measured per stratum
against an mpmath corpus and an envelope of oracles (Sophus, manif, GTSAM, sophus-rs,
`tf_tree_math`), and must dominate it.

**Status:** pre-implementation. The design is in [`docs/`](./docs/); start at
[`docs/PROJECT.md`](./docs/PROJECT.md). Agents: [`CLAUDE.md`](./CLAUDE.md).

Crates: `helicoid-linalg` (the leaf: scalar model, fixed-size linear algebra), `helicoid` (groups,
Jacobians, charts, geodesics). `helicoid-spline` is gated
([`0011`](./docs/decisions/0011-continuous-time-waits-for-a-consumer.md)).

License: MIT OR Apache-2.0.
