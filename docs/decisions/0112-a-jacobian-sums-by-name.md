# 0112: A Jacobian sums by name

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

`Jac` has `identity`, `mul`, `inverse`, `neg`, `apply`, `apply_transpose`, `write_dense` and
`sandwich`. It cannot add two Jacobians or multiply one by a scalar, and nothing builds one from
dense entries. Generic code over `G::Jac` therefore cannot write a sum.

The multivariate chain rule is a sum. A cumulative B-spline's control Jacobians
([`0114`](./0114-a-spline-evaluates-the-window-it-is-handed.md)) are, per control `m`:

- the value: `Σ_j Ad_{P_j⁻¹} λ_j J_r(λ_j d_j) ∂d_j/∂δ_m`, two terms when both `d_m` and `d_{m+1}`
  depend on `δ_m`;
- one recursion step of the velocity: `λ_j ad(Ad_{A_j⁻¹} ω⁽ʲ⁻¹⁾) J_r(λ_j d_j) + λ̇_j I`,
  `A_j = Exp(λ_j d_j)`;
- on Rⁿ: `b_m I`, a scalar times the identity.

None of these can be written as `[G::Jac; K]` today. [`0011`](./0011-continuous-time-waits-for-a-consumer.md)
decision 3 lets the core change for a spline only when a record proves a primitive missing. This is
that record.

**What stands in the way.** `SEn3Jac`'s rustdoc says "There is no `PartialEq`, `Add` or `Sub`".
[`0028`](./0028-two-ready-specs-disagree-on-the-blocks-of-sen3jac.md) item 1 counts "an unlabelled
Jacobian sum is reachable" as a harm of public fields. Its reason is `0002`'s: no composition happens
without naming its side.

## Decision

1. **`Jac` gains two required methods**, named, like `Tangent::add` and `Tangent::scale`:
   ```rust
   /// `self + o`, entry by entry. Both are in the same side's convention; the sum does not check it.
   fn add(&self, o: &Self) -> Self;
   /// `k * self`, entry by entry.
   fn scale(&self, k: S) -> Self;
   ```
   No Jacobian type implements `Add`, `Sub` or `Mul<S>`. API.md R1's operator rule is unchanged.
2. **Every impl stays structured.** Each representation is a linear subspace of the dense matrices,
   so it is closed under both operations:

   | Impl | `add` and `scale` |
   |---|---|
   | `Mat3` (SO(3)), `Mat2` (S²) | entry by entry |
   | `RnJac` (`k I`) | `k + k'`, `s k` |
   | `SEn3Jac` (`diag`, `col`) | `diag` and each `col[i]` |
   | `Sim3Jac` (`[A 0 0; B C f; 0 0 e]`) | `A`, `B`, `C`, `f`, `e` |
   | `ProductJac` | per component |
   | `TwistBlockJac` | delegates to its `ProductJac` |

   Each stored entry is one rounding of the dense operation's entry, so `write_dense(a.add(b))`
   equals `write_dense(a) + write_dense(b)` bit for bit on every stored entry. A structural zero
   is written `+0`, whereas the dense `k · 0` is `−0` for a negative `k` and NaN for a non-finite
   one, so structural zeros are compared as zeros, not as bits.
3. **The rustdocs that deny a sum** change: `SEn3Jac`'s and `ProductJac`'s "There is no
   `PartialEq`, `Add` or `Sub`" read "There is no `PartialEq`, and no `Add` or `Sub` operator: a
   sum is `Jac::add`", and `Jac`'s "closed under `mul`, `inverse` and `neg`" gains `add` and
   `scale`. `0028`'s decision (crate-private blocks) stands. Its item 1 is answered here rather
   than reopened: the fields stay private, and the one sum that exists is named.
4. **`NUMERICS.md` §14** gains rows `SEn3Jac::add`, `scale` and `Sim3Jac::add`, `scale`, each against
   the dense entrywise operation, bit for bit. The `Mat3`, `Mat2`, `RnJac` and product impls are
   the entrywise operation itself, and have no twin.

## Rationale

- **A sum carries no side, as `mul` carries none.** `0028`'s worry is a right Jacobian added to a left
  one. `mul` already composes them without a check, and the guard against mislabelling has always
  been the `*_jacobians` method names (R1), not the missing sum. Leaving out `add` does not stop the
  bug. It only stops the chain rule.
- **A spline-local Jacobian type** would keep the core unchanged. But `0005` and API.md R5 make
  `G::Jac` the Jacobian everywhere, a consumer composing a spline Jacobian with a group's would need
  a conversion, and the spline would re-derive each group's structure.
- **`from_dense` instead of `add`** would let any sum be built densely. That puts back the dense
  path R5 closes, and it loses the structure.
- **One `axpy`-style method** (`self + k o`) is `scale` then `add`, two roundings either way, since
  D16 forbids fused multiply-add. It would be a second spelling (API.md §6 item 5).

## Consequences

- `Jac` is not sealed, so an external impl must add both methods. The workspace is `0.0.x`.
- Every new `Jac` impl owes `add` and `scale` and keeps the subspace property. A representation not
  closed under sums is not a `Jac`.
- `Dual` flows through both with no extra code.

## Implementation plan

1. `Jac::add` and `Jac::scale` on the seven impls and on `heis_tests`' `HJac`; the rustdocs of
   decision 3; the `NUMERICS.md` §14
   rows and `reference::{sen3jac, sim3jac}_{add, scale}`. Verified by
   `jac_add_and_scale_match_dense` (a proptest per impl at `f64`, `f32` and `Dual`, bit for bit)
   and `*_matches_reference` for the two §14 rows; `just dev helicoid`, `just test`, `just no-std`.

## Open questions

None. A sum carries no side marker, as `mul` carries none (Rationale).
