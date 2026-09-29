# 0007: The budget a foundation can afford

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

`helicoid` sits under crates with strict budgets: `tf_tree` D14 (`tf_tree_core` = `libm` +
`bytemuck` + `blake3`; `tf_tree_math` = `libm` + `bytemuck`), omnisac's leaf invariants (`no_std`,
faer-free, `forbid(unsafe_code)`), fuse-geometry on wasm32. `tf_tree` 0016 is the measured warning:
`pulp` would have taken `tf_tree_math` from 2 dependencies to 11 for ~12% on one batch path, and
`-C target-cpu=x86-64-v3` made that path slower.

## Decision

1. **Dependencies.** `helicoid-linalg` = `libm`. `helicoid` = `helicoid-linalg` + `libm`.
   Optional feature `mint` on both. Nothing else: no `nalgebra`, `faer`, `num-traits`, `bytemuck`,
   `serde`, SIMD crate, or logging.
2. **`no_std`, no `alloc`**, `#![forbid(unsafe_code)]` on every library root.
3. **`libm` for every transcendental on every target** (D16); no `mul_add`; no `target-cpu` in any
   config or CI job.
4. **MSRV never exceeds the lowest consumer MSRV** (D17): 1.87 today.
5. `xtask` and `runners/` are unrestricted, never a normal dependency of a library crate; runners
   are excluded from the workspace.
6. `cargo xtask lint` asserts each library crate's normal-dependency closure equals item 1's set;
   `deny.toml` bans the list above for library crates.

## Rationale

With this budget, adding `helicoid` to `tf_tree` adds **no new third-party crate** (`libm` is
already there), which is what makes D1 affordable for every consumer. `mint` is the interop
standard nalgebra, glam and cgmath already implement, so no per-library feature is needed.
Alternatives lost: a `nalgebra` feature (feature unification would pull it into consumers that
forbid it); `std` math (platform-dependent bits, D16).

## Consequences

- SIMD lanes need their own record and must beat 0016's bar (0013).
- `bytemuck` stays in consumers that store `helicoid` values in their own `Pod` records (D2).

## Implementation plan

1. Workspace, `deny.toml`, the closure lint — verified by the planted-`nalgebra` test
   (`PHASE1.md` §11).

## Open questions

None.
