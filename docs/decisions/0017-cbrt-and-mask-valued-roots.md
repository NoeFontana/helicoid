# 0017: `Real::cbrt` and mask-valued roots for `solve_cubic`

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

`PHASE2.md` §6 asks for omnisac's `solve_cubic` "with its current signature, unchanged" so that its
migration is a pure move. That signature returns `ArrayVec<f64, 3>`, which cannot be kept:
`arrayvec` is outside the `0007` budget, and a variable root count derived from float comparisons
cannot be a `usize` under the scalar model (`0003`, `API.md` R4: an answer that depends on a
comparison is an `S::Mask`). Its one-real-root arm needs a cube root, which `Real` lacks
(`sqrt`, `sin_cos`, `atan2`, `abs`, `copysign` only). Its tolerances are `f64` literals with no
`f32` analogue, and `Real::lit` accepts only exactly representable constants.

## Decision

1. **`Real::cbrt`**, routed through `libm::cbrt` / `cbrtf` (D16). `Dual`: `d/(3·c²)` with
   `c = cbrt(v)`; like `sqrt` it is infinite at `v = 0` where `d ≠ 0`, and the safe-argument pattern
   applies. `NUMERICS.md` §12 gains a `cbrt` note (no domain restriction; derivative singular at 0);
   `PHASE2.md` §2 (signature), §3 (the `Dual` rule and its zero) and §8 (its derivative check) and
   `API.md` §2 list it; `dual_value_is_plain_value` covers it.
2. **`solve_cubic` returns its roots with a validity mask**: a `Vec3<S>` of roots and a
   `[S::Mask; 3]` in which slot `k` is a root iff its mask is set (a thin named wrapper struct is
   permitted; no `ArrayVec`, no `usize` count). The algorithm is otherwise omnisac's, ported as is
   with every branch through `S::branch`/`S::select`. omnisac's migration therefore needs a small
   adapter (masks to `ArrayVec`), not a pure move; `PHASE2.md` §6 and §9 are amended accordingly.
3. **Tolerances are per-precision power-of-two literals**, exactly representable in both
   precisions, stated in rustdoc with the omnisac value they replace and the reasoning
   (a fixed multiple of the unit roundoff). They are not coefficients of `NUMERICS.md` §4 and are
   not generated; changing one is a changelog line naming the function.
4. `eig3` and `svd3` are ported the same way: omnisac's algorithm, generic over `S: Real`.

## Rationale

`cbrt` is a primitive with a correctly-defined derivative and a `libm` implementation on every
target, cheaper and more accurate than composing it from `exp`/`ln` (which `Real` lacks). A mask
array keeps `solve_cubic` lane-correct and AD-correct. Alternatives lost: amending `0007` for
`arrayvec` (a dependency for one signature); a `usize` count restricted to `Real<Mask = bool>`
(splits the scalar model).

## Consequences

- `Real` gains a method: every `Real` impl (`f64`, `f32`, `Dual`) implements it.
- omnisac's migration adapts the root container at its call sites (`PHASE2.md` §9).

## Implementation plan

1. `Real::cbrt`, the `Dual` rule, the `NUMERICS.md` §12 note and docs edits — verified by
   `dual_value_is_plain_value` and a derivative check against mpmath.
2. `solve_cubic` — verified by the strata of `PHASE2.md` §6 against mpmath and a differential
   against omnisac.
3. `eig3`, then `svd3` — verified likewise (`PHASE2.md` §6, §8).

## Open questions

None.
