# 0114: A spline evaluates the window it is handed

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #226, #228, #231

## Context

[`0113`](./0113-phase-7-opens-for-a-calibration-consumer.md) opens Phase 7: a uniform cumulative
B-spline of order `K` on any `LieGroup`, with value, velocity, acceleration and their Jacobians on
every control. The caller owns its control points and its time. `helicoid` has no time type, and
generic code cannot compare a float (D4).

The formulas are `NUMERICS.md` §21, with derivations in `docs/maths/splines.md`. Every one was
checked against central differences at 50 digits on SO(3) and SE(3), at `K = 6`, on both sides.
The worst Jacobian entry differs by `2e-25`, and the basis and Gram values were checked in exact
rationals.

## Decision

1. **A new crate, `helicoid-spline`.** It is `no_std`, has no `alloc`, has
   `#![forbid(unsafe_code)]`, and depends on `helicoid` and `helicoid-linalg` only. It calls no
   transcendental itself, so it needs no `libm`: every `sin`, `atan2` and `sqrt` reaches D16's
   kernel through the groups. The `cargo xtask lint` dependency closure gains the crate.
2. **The caller hands over a window.** The segment index and `u` come from the caller's own time,
   so no float comparison is made here. The controls stay in the caller's storage. Arrays have
   length `K` with slot 0 unused, because `K − 1` as a length needs `generic_const_exprs`.
   ```rust
   pub struct Order<const K: usize>;
   pub trait ValidOrder: sealed::Sealed {}          // Order<K> for the orders 0115 generates

   /// `D` (0, 1 or 2, asserted at monomorphization) is the highest derivative formed; the rows
   /// above it are zero and cost nothing. Every evaluation below takes its `D` from here.
   pub struct Weights<S, const K: usize, const D: usize> {
       pub lambda: [S; K], pub dot: [S; K], pub ddot: [S; K],
   }
   impl<S: Real, const K: usize, const D: usize> Weights<S, K, D> where Order<K>: ValidOrder {
       /// λ, λ̇ = h⁻¹ dλ/du, λ̈ = h⁻² d²λ/du², by Horner on 0115's rows.
       pub fn at(u: S, inv_h: S) -> Self;
   }

   /// `x0` and `d_j = Log(X_{j-1}⁻¹ X_j)`, once per segment.
   pub struct Segment<S, G, const K: usize> { /* private */ }
   /// A `Segment` plus `J_r⁻¹(d_j)` and `J_l⁻¹(d_j)`, for Jacobian evaluation.
   pub struct JacSegment<S, G, const K: usize> { /* private */ }

   impl<S: Real, G: LieGroup<S>, const K: usize> Segment<S, G, K> where Order<K>: ValidOrder {
       pub fn new(window: &[G; K]) -> Self;
       pub fn value<const D: usize>(&self, w: &Weights<S, K, D>) -> G;
       /// Orders above `D` are zero.
       pub fn kinematics<const D: usize>(&self, w: &Weights<S, K, D>) -> Kinematics<S, G>;
   }
   impl<S: Real, G: LieGroup<S>, const K: usize> JacSegment<S, G, K> where Order<K>: ValidOrder {
       pub fn new(window: &[G; K]) -> Self;
       pub fn segment(&self) -> &Segment<S, G, K>;
       pub fn kinematics_jacobians<Sd: Side, const D: usize>(&self, w: &Weights<S, K, D>)
           -> (Kinematics<S, G>, ControlJacobians<S, G, K>);
   }

   pub struct Kinematics<S, G: LieGroup<S>> { pub value: G, pub vel: G::Tangent, pub acc: G::Tangent }
   /// Row `m` is the Jacobian with respect to control `m`. Orders above `D` are zero.
   pub struct ControlJacobians<S, G: LieGroup<S>, const K: usize> {
       pub value: [G::Jac; K], pub vel: [G::Jac; K], pub acc: [G::Jac; K],
   }
   ```
3. **The side is named (R1).** `kinematics_jacobians::<Right, _>` perturbs every control as
   `X_m Exp(δ_m)` and the value as `X Exp(ε)`. That is the backward accumulation `NUMERICS.md`
   states. `::<Left, _>` is derived from it: `J_L = Ad_X J_R Ad_{X_m⁻¹}` for the value, and
   `J_R Ad_{X_m⁻¹}` for the velocity and acceleration, which are tangents with no output side.
   This follows `act_jacobians::<Sd>` and `PHASE4.md`'s `geodesic_jacobians::<Sd>` (owed).
4. **One generic recursion, no per-group code.**
   - Products and Rⁿ need no override. `Product`'s operations are componentwise and `RnJac` is
     `k I`, so the generic code on `Product<SO3, R3>` equals its two components to the bit.
   - On Rⁿ, `Ad = I`, `ad = 0` and `J = I` reduce the recursion to the non-cumulative weights
     `b_m = λ_m − λ_{m+1}`.
   - Bitwise equality with the components is a proptest. Agreement with `Σ b_m c_m` is a proptest
     within a stated bound. It is not bitwise: the cumulative form forms `c_j − c_{j−1}` first.
   - The one core change the recursion needs is [`0112`](./0112-a-jacobian-sums-by-name.md).
5. **Speed is designed in, not tuned later** (`0113`: no budget, as fast as measured).
   - `Log`, `J_r⁻¹` and `J_l⁻¹` run once per segment.
   - Value and kinematics never compute the inverse Jacobians, because only `JacSegment` holds them.
   - Jacobians are `O(K)` by backward accumulation.
   - Above the requested `D`, nothing is formed, `Weights::at` included: a value-only caller pays
     one Horner pass per weight.
   - A further kernel change is adjudicated in a consumer-shaped loop (`.claude/rules/numerics.md`),
     against the bench rows of [`0116`](./0116-the-spline-is-measured-before-it-ships.md).
6. **Domains (R6, D11).**
   - `u` is any finite value: outside `[0, 1]` the polynomial extrapolates.
   - `inv_h` is finite and positive, `debug_assert!`ed.
   - Consecutive controls whose relative rotation nears π inherit `Log`'s documented behaviour.
     There is no failure value, because `Log` is total.
7. **Reference twins (D6).** `reference::{spline_value, spline_kinematics, spline_jacobians}`
   accumulate forward per control, in `O(K²)`, with dense Jacobians. They are rows of
   `NUMERICS.md` §14.

## Rationale

- **A trait over the groups with an Rⁿ override** would compute `Σ b_m c_m` directly, but a product
  could not use it. `ProductJac`'s blocks are private (`0025`), so an outside crate cannot assemble
  one from component Jacobians. That would mean a second core change, and the product-equals-components
  property would hold only by care rather than by construction.
- **One `Segment` holding the inverse Jacobians** (the first sketch) charges value-only callers two
  `3 × 3` series evaluations per increment for nothing.
- **Right only** would leave R1's side implicit in a method name. The `Left` conversion costs
  `K + 1` adjoints and up to `4K` Jacobian products (two per value Jacobian, one per velocity and
  acceleration Jacobian), and no new formula beyond `Ad`.

## Consequences

- API.md gains the crate's surface, which answers §6's checklist:
  1. R1 holds: the side is a type parameter. R3 holds: tangents and Jacobians are the groups' own
     types, and the only arrays are indexed by control, as the rustdoc states.
  2. The kernels have corpus ids (`0116`) and twins.
  3. Domains are stated and `debug_assert!`ed.
  4. No `bool`, no dense Jacobian, no unlabeled array.
  5. No second spelling: `geodesic` is the `K = 2` linear case, and stays the way to write it.
  6. Ownership holds (`0009`, `0011`).
  7. The crate is `no_std`, allocation-free and D16-bit-identical.
  8. It removes a sign and a side in each Jacobian that a consumer would otherwise hand-derive.
- `just dev` gains a `spline` target.

## Implementation plan

1. The crate skeleton, `Order`, `Weights`, `0115`'s tables and the three reference twins of
   decision 7, so that `0116`'s seeded defects run before any kernel. Verified by
   `weights_match_exact_rationals`, the twins against `0116`'s corpus, `just no-std`, `just msrv`,
   and the lint's dependency closure rejecting a planted `alloc`.
2. `Segment`, `value` and `kinematics`. Verified by
   `spline_*_matches_reference`, `Dual` through `value` for `vel` and `acc`, continuity across a knot
   to order `min(K − 2, 2)`, and a product equal to its components bitwise.
3. `JacSegment` and `kinematics_jacobians` at both sides. Verified by `spline_jacobians_matches_reference`
   and `Dual` along `X_m Exp(ε e_k)`.
4. The shipped conformance subject answers `0116`'s ids. Verified by `just conformance` at both
   precisions.

## Open questions

None. `Kinematics` and `ControlJacobians` read zero above `D` rather than splitting the entry points
by order. `Weights` carries `D`, so `inv_h` is unused at `D = 0` and the caller passes any finite
positive value.

## Further work

- `J_l⁻¹(d) = J_r⁻¹(d) − ad_d` exactly on every Lie group: the two series differ only in the odd
  Bernoulli terms, and every odd one past `B₁` vanishes. `JacSegment::new` could form `J_l⁻¹` from
  `J_r⁻¹` with one `ad` and one `add`, not a second series. That saves `K − 1` series evaluations
  per segment, but it rounds differently near the switch points. It needs a `NUMERICS.md` identity
  and a measurement against the separate evaluation.
