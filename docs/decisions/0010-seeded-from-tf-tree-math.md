# 0010: Seeded from `tf_tree_math`, and measured against it

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

`tf_tree_math` is the best-tested Lie code in the stack: `log_so3` through the quaternion `atan2`,
a measured small-angle switch (its D12), reference twins with a 10⁵-pair differential proptest at
1e-14 for ScLERP (its D13), left/right invariance tests including a positive failure for
`LerpSlerp`, and a hardcoded high-precision sweep over $[10^{-12}, \pi]$. It is published, `no_std`,
`forbid(unsafe_code)`, and its conventions are the ones 0002 adopts. locus-tag contributes one
regression worth keeping: nalgebra's iterative matrix-to-quaternion loop that never converges on
degenerate input, which its `quat_from_so3` avoids with Shepperd's closed form.

## Decision

1. Phase 3's `SO3`/`SE3` and Phase 4's geodesic and dual-quaternion fast twin are **ports** of
   `tf_tree_math`, generalized to `S: Real`; the porting PRs cite the `tf_tree` commit.
2. `tf_tree_math`'s tests are carried over in spirit and superseded in form: its high-precision
   table by the corpus, its differential proptest by `se3_geodesic_matches_reference` and
   `tf_tree_math_differential`, its invariance tests verbatim.
3. `tf_tree_math` is **oracle #1** (`runners/tf_tree_math`). `helicoid` must dominate it on every
   paired stratum before `tf_tree` migrates (`PHASE4.md` §5.2).
4. The migration itself, the amendment of `tf_tree` D12/D13/D14, and `tf_tree_math`'s fate are
   **`tf_tree`'s record**, not this one.
5. locus-tag's degenerate matrices are `from_matrix_never_iterates` fixtures.

## Rationale

Porting the tested implementation, then measuring the port against the original, is cheaper and
safer than a clean-room implementation that must rediscover the same traps. Keeping the original
as an oracle turns "we did not make it worse" into a bar rather than a hope.

## Consequences

- If `tf_tree_math` wins a stratum, the phase stops until fixed or explained by record.
- `tf_tree` keeps `LerpSlerp` and its `Interp` trait (0009).

## Implementation plan

1. `runners/tf_tree_math` (Phase 1) — verified by envelope rows.
2. SO(3)/SE(3) port (Phase 3) — verified by domination over the runner.
3. Geodesic and fast twin port (Phase 4) — verified by the invariance tests and domination.

## Open questions

None.
