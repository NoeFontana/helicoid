# 0005: The Jacobian is a dual matrix

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

Every SE_N(3) Jacobian that matters — $\mathrm{Ad}$, $\mathrm{ad}$, $J_r$, $J_l$ and their
inverses — has the rotation-first shape $\begin{bmatrix} A & 0 \\ B_i & A \end{bmatrix}$
(`NUMERICS.md` §2.2, §5). Dense 6×6 algebra ignores it: 216 multiplications per product, a
separate closed form for $J_r^{-1}$ to derive and get wrong, and for SE₂(3) 729 multiplications.
Stable Rust cannot express `Matrix<S, {3 + 3N}, {3 + 3N}>` (`generic_const_exprs` is unstable).

## Decision

1. `G::Jac` is a structured associated type. SE_N(3): `SEn3Jac<S, N> { diag: Mat3, col: [Mat3; N] }`
   with $(A + \epsilon B)(C + \epsilon D) = AC + \epsilon(B_iC + AD_i)$ and
   $(A + \epsilon B)^{-1} = A^{-1} - \epsilon A^{-1}B_iA^{-1}$. SO(3): `Mat3`. SO(2): scalar.
   SE(2): dense `Mat3`. Rⁿ: identity-structured. `Product`: block-diagonal `ProductJac`.
   Sim(3): dense 7×7 (the scale column breaks the dual structure).
2. $J_r^{-1}(\tau)$ and $J_l^{-1}(\tau)$ for SE_N(3) are **the dual-matrix inverse** of $J_r$,
   $J_l$, using SO(3)'s closed-form inverse for $A$. No second closed form exists in the code.
3. Dense matrices appear only through `write_dense(&mut StridedMut)` and `sandwich` into a
   fixed-size covariance (API R5).
4. Const-generic workarounds: tangents hold `[Vec3; N]`; `DOF` is an associated const on both
   `Tangent` and `LieGroup`, tied by a `const` assertion; `Gaussian<S, G, Sd, const D>` and
   `sandwich<const D>` carry `D` as a separate parameter with `const { assert!(D == G::DOF) }`.

## Rationale

Product cost $27 + 54N$ vs $(3 + 3N)^3$: 81 vs 216 for SE(3), 135 vs 729 for SE₂(3). The inverse
comes from algebra, so one formula ($Q$) covers four Jacobians. The dense twin stays as the
reference (`NUMERICS.md` §14), so nothing is taken on faith. Alternatives lost: dense
`[[S; 6]; 6]` (cost, extra formulas); nightly `generic_const_exprs` (MSRV, D17).

## Consequences

- Solvers receive dense blocks by writing, never by return value.
- Adding a group means stating its Jacobian structure in `NUMERICS.md` first.

## Implementation plan

1. `SEn3Jac` algebra with dense reference twins — verified by `sen3jac_mul_matches_reference`,
   `sen3jac_inverse_matches_reference`.
2. `jr_inv` via the inverse — verified by the `sen3_jr_inv_n*` corpus ids.

## Open questions

None.
