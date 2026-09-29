# 0003: The scalar that cannot say "less than"

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

The same numeric code must run on `f64`, `f32`, `Dual<S, N>` (for Jacobians of the shipped code),
and possibly SIMD lanes later (0013). The usual pattern, `if theta < EPS { series } else { exact }`,
is wrong for three of those: a lane vector has no single truth value; reverse-mode AD through a
branch taken on a traced value is unsound; and the exact arm evaluated at `θ = 0` (as lanes must)
produces `NaN`, which a later blend does not clean up in every AD system. `tf_tree` 0053 adds the
codegen evidence: a mask-select written to be branchless compiled to a branch anyway — so a model
cannot promise branchless code, only lane- and AD-correct code.

## Decision

1. `Real` (`PHASE2.md` §2) has **no `PartialOrd` and no `PartialEq`**. Generic code cannot write
   `<`; comparisons are `lt`/`le` returning `S::Mask`.
2. `Mask::decide` owns the branching policy: `bool` evaluates one arm, a lane mask evaluates both
   and blends. `Real::branch` is provided on top of it; `Real::select` blends values.
3. **The safe argument**: the exact arm evaluates at `select(small, 1, θ²)`, so no arm is ever
   non-finite.
4. `Blend` is implemented for every value type, so `branch` can return tuples of coefficients and
   whole group elements.
5. `Dual<S, N>` compares on its value part and **its value path is bitwise identical** to plain
   `S` evaluation.
6. `PRECISION` is an associated const; generated constants are selected by it at monomorphization.
   `lit` is for exactly representable constants only.
7. The one `if` on a float comparison in the workspace is inside `impl Mask for bool`.

## Rationale

A lint for "no `if` on floats" is syntactic and leaky; removing `PartialOrd` makes the wrong code
not compile. Putting the policy in `Mask` keeps scalar code lazy (one arm) while making lane code
correct by construction. Alternatives lost: `num-traits::Float` (brings `PartialOrd`, a dependency,
and no mask concept); a `simd` feature flag switching implementations (two code paths to verify).

## Consequences

- Every numeric routine is generic over `S: Real`; concrete `f64` code in a library crate is a
  review failure outside the `Real` impls.
- Public API returns `S::Mask` where an answer depends on a comparison (API R4).
- A SIMD impl (0013) needs no change to any algorithm.

## Implementation plan

1. `Mask`, `Real`, `Blend`, `Precision`, `f64`/`f32` impls — verified by the `compile_fail`
   doctest `real_has_no_partial_ord` and `just no-std`.
2. `Dual` — verified by `dual_value_is_plain_value`.
3. The seeded `NaN` defect — verified by `just conformance --self-test`.

## Open questions

None.
