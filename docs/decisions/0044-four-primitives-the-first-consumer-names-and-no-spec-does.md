# 0044: Four primitives the first consumer names and no spec does

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** landed. `Quat::{dot, norm}`, `SEn3::renormalize`, `SEn3<S, 1>::mul_inv`, `API.md` §3's two rows, `NUMERICS.md` §14's `mul_inv` row and its twin test. Item 5's debt stands: no corpus id for any of the four, owed before `0.1.0`

## Context

[`0041`](./0041-the-integration-is-an-adapter-at-tf-tree-math.md) enumerates, function by function,
what `tf_tree_math`'s bodies delegate to. Four of its names have no counterpart here and no row in
any spec — not in `PHASE3.md` §4 or §5's NORMATIVE lists, not in `API.md` §3, not in
`PHASE4.md` §5.1's parity table:

| `0041` | wave | `helicoid` today |
|---|---|---|
| `Quat::dot` | 1 | nothing |
| `Quat::norm` | 1 | `Quat::norm_sq` only |
| `Iso3::normalized` | 2 | `SO3::renormalize` and `Quat::renormalize`; **nothing on `SEn3`** |
| `Iso3::mul_inv` | 2 | nothing; `rminus`/`lminus` go through `Log`, a different operation |

Three of the four are arithmetic with no decision in them. The interesting ones are the last two.

**`SEn3` cannot be renormalized from outside at all.** Its `q` field is private
(`0028` option A narrowed `SEn3Jac`, and `SEn3`'s fields were never public), `Mul`'s own rustdoc
says "the product quaternion is not renormalized … `SO3::renormalize` is the caller's step
(`0027`)", and the only route a caller has is `parts()` → `SO3::renormalize` → `from_parts`. So the
type documents a caller obligation it gives the caller no way to discharge in one step. That is a
hole, not a missing convenience.

**`mul_inv` is not a spelling convenience.** `tf_tree_math`'s is documented as "compute
`self · rhs⁻¹` directly, without materializing `rhs⁻¹`": rotation `q_a q_b*`, translation
`t_a − (q_a q_b*)·t_b`. The composition `a * b.inverse()` instead computes
`t_a + R_a·(−R_b^{\mathsf T}t_b)` — **two** rotations of a vector where `mul_inv` does one. The two
are equal in exact arithmetic and **not bit-identical**, because `R_a(R_b^{\mathsf T}t_b)` rounds
twice where `R(q_aq_b^*)\,t_b` rounds once. So it is a routine with its own rounding, which under
D6 means its own twin, and under D8 its own scored rows if it ever gets a corpus id.

A tree walk computes relative transforms constantly, so the saved rotation is on the hot path of
the first consumer.

## Decision

1. **`Quat` gains `dot` and `norm`**, and `API.md` §3's `Quat` row is corrected while we are in it:
   it lists five constructors and omits the shipped `conjugate`, `norm_sq`, `to_matrix` and
   `renormalize`, so it is already stale and `PHASE3.md` §4's list is the fuller one.
   - `pub fn dot(&self, o: &Self) -> S` — `w w' + x x' + y y' + z z'`, summed left to right in that
     order, which is `Quat::norm_sq`'s order at `o == self` and is asserted to be so by a test.
     No domain.
   - `pub fn norm(&self) -> S` — `norm_sq().sqrt()`. Two roundings, not one; the rustdoc says so,
     because the alternative (a scaled hypot) is a different routine and nothing in the stack has
     asked for it. **Domain:** `norm_sq()` finite — the squaring overflows for
     $\lVert q\rVert \gtrsim 1.3\times10^{154}$, which is a documented `# Domain` and no release
     check (D11).
2. **`SEn3` gains `renormalize`**, `&mut self`, delegating to `SO3::renormalize` on the rotation
   and leaving the columns alone, exactly as `SO3::renormalize` delegates to `Quat`'s. `API.md` R2
   makes `&mut self` the shape for an in-place operation, and `0027`'s domain carries unchanged:
   one Newton step is a normalization only within $\lvert\eta\rvert \le 2^{-26.29}$ (`f64`) or
   $2^{-11.79}$ (`f32`), which is what a chain of compositions drifts into and not what an
   arbitrary quaternion sits in.
   There is **no** `normalized(self) -> Self`: `tf_tree_math`'s `Iso3::normalized` is a by-value
   wrapper the adapter writes in one line over this, and `API.md` R3 keeps one spelling per
   operation.
3. **`SEn3<S, 1>` gains `mul_inv`**, by `0041`'s reading, not `tf_tree`'s spelling:

   ```rust
   /// `self · other⁻¹`, without forming `other⁻¹`.
   pub fn mul_inv(&self, other: &Self) -> Self
   ```

   Defined for `N = 1` only, where `SEn3`'s single column is a translation and the saving is one
   `act`; for general `N` the same algebra holds column by column and the row is added when a
   consumer asks (`PROJECT.md` §5.1 is where that argument goes, not here).
4. **It is a `NUMERICS.md` §14 row with a reference twin**, the twin being `a * b.inverse()` —
   the composition it is *not* bit-identical to, which is exactly what makes the twin worth having.
   `mul_inv_matches_reference` is the proptest, its bound measured and then no-regress, following
   `sen3_jr_inv_matches_reference_n1`'s shape.
5. **No corpus id for any of the four.** `PHASE3.md` §0.0 already owes "a `renormalize` corpus id
   and any `Quat` stratum of its own"; these four join that debt rather than inventing a stratum
   family each, which §4.4 makes a record in its own right. D7's bar is that a *public numeric
   routine* without a stratum is unverified and **may not ship in a release** — `helicoid` is
   unpublished, and [`0041`](./0041-the-integration-is-an-adapter-at-tf-tree-math.md) decision 7's
   `0.0.1` is the release that makes this due, so the debt is named with a deadline rather than
   left open: the `quat_*` stratum family is owed before `0.1.0`, not before `0.0.1`.

## Rationale

The alternative for each is "let the adapter write it", and for two of the four that is right —
which is why `normalized` is not on this list and `rotate`, `norm_squared` and `normalize` map onto
shipped items under other names.

It is wrong for `mul_inv`, because D1 forbids a second implementation of a group operation in a
consumer and this one has its own rounding: an adapter writing `q_a q_b*` and one `act` would be
implementing an SE(3) composition in `tf_tree_math`, which is the thing this whole migration
exists to stop. It is wrong for `renormalize`, because the data a caller needs is behind a private
field, so "the adapter writes it" means a `parts`/`from_parts` round trip on a hot path to work
around an encapsulation we chose. And `dot`/`norm` are genuinely trivial, but `Quat` is a public
type whose arithmetic is specified to the bit elsewhere (`norm_sq` has committed left-to-right
goldens), so two more members with stated summation order cost less than two consumers guessing.

## Consequences

- Four public items, one `NUMERICS.md` §14 row, one `API.md` §3 row corrected.
- The `quat_*` corpus debt gains a deadline (`0.1.0`) and three more routines.
- `SEn3::renormalize` makes `Mul`'s existing rustdoc discharge-able as written; no behaviour of
  `Mul` changes, and nothing is renormalized implicitly (`0027`).
- `mul_inv` is **not** a drop-in for `a * b.inverse()` in a bit-identity test. Any existing
  assertion comparing the two forms must name which it means; the twin test is the one place they
  are compared, with a tolerance.

## Implementation plan

1. This record and `API.md` §3's corrected `Quat` row — verified by `just lint`.
2. `Quat::{dot, norm}` with the summation-order test against `norm_sq` — verified by `just test` at
   `f64`, `f32` and `Dual<f64, 4>`, as `quat_tests.rs` already covers the rest of the type.
3. `SEn3::renormalize` — verified by a test that a drifted product returns to unit norm with
   `0027`'s quadratic convergence, mirroring `so3_tests`' step test.
4. `SEn3::mul_inv`, its §14 row and `mul_inv_matches_reference` — verified by the proptest with its
   measured bound, and by a test that it is **not** bit-identical to `a * b.inverse()` on at least
   one pair, so the claim in this record stays true.

## Open questions

None.

## Further work

1. The `quat_*` stratum family — `norm`, `dot`, `renormalize` and `SEn3::renormalize` — before
   `0.1.0`. `PHASE3.md` §0.0 owes it already; this record dates it.
2. `mul_inv` for general `N`, when a consumer has one.
3. Whether `Quat` owes `inverse` for a non-unit quaternion (`conjugate`'s rustdoc says it *is* the
   inverse of a unit quaternion and nothing covers the rest). No consumer asks; `API.md` §6's
   question when one does.
