# 0043: The geodesic Jacobian ships the cancellation-free form, and `Jac` gains one method

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** none yet; `PHASE4.md` §2 is the deliverable

## Context

`NUMERICS.md` §10 states the right geodesic Jacobians as

$$\partial X(t)/\partial X_1 = t\,J_r(td)\,J_r^{-1}(d),\qquad
  \partial X(t)/\partial X_0 = \mathrm{Ad}_{\mathrm{Exp}(-td)} - t\,J_r(td)\,J_r^{-1}(d)\,\mathrm{Ad}_{\Delta^{-1}}.$$

`docs/maths/geodesics.md` GE.7(a) proves an **equal** form for the second,
$J_0 = (1-t)\,J_l\big((1-t)d\big)\,J_l^{-1}(d)$ — $J_1$ of the swapped pair $(X_1, X_0)$ at
$1-t$ — and GE.7(d) derives the left pair, which §10 does not state at all although `PHASE4.md`
§2 takes a side.

The two forms are not interchangeable in floating point. §10's is a difference of two $O(1)$
matrices that tends to zero like $1-t$, so its absolute error is $O(u\lVert\mathrm{Ad}\rVert)$ and
its **relative** error grows like $1/(1-t)$. Measured (SE(3), $\theta(d) = 1$, both operands
rounded to `f64`, one sample, mpmath 1.3.0 at 110 digits for the exact value; protocol in GE.7's
*Checked* block, scripts not committed):

| $1-t$ | $10^{-1}$ | $10^{-3}$ | $10^{-5}$ | $10^{-9}$ | $10^{-12}$ |
|---|---|---|---|---|---|
| §10's difference | $19\,u$ | $820\,u$ | $1.5\times10^5\,u$ | $9.7\times10^8\,u$ | $1.3\times10^{12}\,u$ |
| GE.7(a) | $\le 0.72\,u$ | $\le 0.72\,u$ | $\le 0.72\,u$ | $\le 0.72\,u$ | $\le 0.72\,u$ |

GE.7(b) also asks for $(J_0, J_1) = (0, I)$ at $t = 1$. GE.7(a) gives it **exactly**, by
$x\cdot(+0)$; the difference form leaves a residue of the order of its own error. And `PHASE4.md`
§4's own strata sample $t = 1 - 10^{-9}$, where the table reads $9.7\times10^8\,u$ — so this is not
a corner a bar would miss, it is a stratum the spec already requires.

The arrangement is not free to adopt: `Jac` (`crates/helicoid/src/traits.rs`) has `identity`,
`mul`, `inverse`, `neg`, `apply`, `apply_transpose`, `write_dense` and `sandwich`, and **no scalar
multiple**. `0005` and `0025` sealed that surface deliberately, so either form needs an addition.

## Decision

1. **`NUMERICS.md` §10 ships GE.7(a)**, with the §10 difference form kept beside it as the
   definition it is equal to, and the measurement above as the reason the implementation takes the
   second. §10 also gains the **left** pair from GE.7(d), collapsed onto the left tangent
   $d_L = X_1 \ominus_L X_0 = \mathrm{Ad}_{X_0}d$:

   | | right, $d = X_1 \ominus_R X_0$ | left, $d_L = X_1 \ominus_L X_0$ |
   |---|---|---|
   | $\partial X(t)/\partial X_1$ | $t\,J_r(td)\,J_r^{-1}(d)$ | $t\,J_l(t\,d_L)\,J_l^{-1}(d_L)$ |
   | $\partial X(t)/\partial X_0$ | $(1-t)\,J_l((1-t)d)\,J_l^{-1}(d)$ | $(1-t)\,J_r((1-t)d_L)\,J_r^{-1}(d_L)$ |

   The left column is GE.7(d)'s $\mathrm{Ad}_{X_i}$-conjugated forms with the conjugations
   collapsed, which is legitimate because every factor is a power series in $\mathrm{ad}_d$ and
   $\mathrm{Ad}_{\Delta} = \mathrm{Ad}_{\mathrm{Exp}\,d}$ is one too, so they commute; that step is
   the load-bearing one and the PR verifies it numerically rather than on the page.
2. **`Jac` gains exactly one method**, mirroring `Tangent::scale`:

   ```rust
   fn scale(&self, k: S) -> Self;
   ```

   **No `add`, no `sub`, and no fused `k·self + l·other`.** This is not minimalism for its own
   sake: `Ad - J\,Ad` is precisely the shape an FMA contracts, and D16 forbids the contraction —
   so *not having the method* is the cheapest way to guarantee nobody writes one. A fused
   `affine` would additionally need $0\cdot x$, which is NaN for a non-finite $x$ and so poisons
   the degenerate uses the way `Matrix::mul` already may not fold.
3. **`geodesic_jacobians::<Sd: Side>(x0, x1, t)` is a provided body**, as `API.md` §3 lists it,
   selecting the column of decision 1's table by `match Sd::IS_RIGHT` — the existing precedent of
   `SO3::compose_jacobians` and `SEn3::compose_jacobians`. It resolves at monomorphization and
   neither arm survives. Both rows are one expression, `J::mul` then `J::scale`, because
   $J_0$ *is* $J_1$ of the swapped pair at $1-t$: swapping negates $d$, and $J_l(-x) = J_r(x)$
   exchanges the two coefficient pairs. That sentence is the whole derivation and belongs in the
   rustdoc.
4. **`mul` and `scale` stay two rounded steps.** `A.mul(&B).scale(t)` must not become a
   `mul_scaled` forming $t\sum_k a_{ik}b_{kj}$ in one expression: that is reassociable and
   FMA-contractible, and it would break `SEn3Jac::mul`'s bit-identity with
   `reference::sen3jac_mul`, which a §14 twin test pins.
5. **Two groups override, for a saving that already exists.** `SEn3::geodesic_jacobians` goes
   through `sen3.rs`'s `inverses`, which returns both $\ominus$ inverses from **one** evaluation of
   `q_words`' seven matrix products — the saving `rminus_jacobians` already banks. `SO3` uses
   $J_l^{-1}(\varphi) = J_r^{-1}(\varphi)^{\mathsf T}$ bit for bit, which `SO3::rminus_jacobians`
   already relies on. Only the evaluations **at $d$** fuse; those at $td$ and $(1-t)d$ have
   different arguments and do not.
6. **`Product` overrides `geodesic`, `geodesic_velocity` and `geodesic_jacobians`** as its factors'
   bodies side by side. The provided bodies are already componentwise, because `Product` overrides
   `rplus`/`rminus`; the override exists so that a factor which overrides `geodesic` with a fast
   twin is the one the product uses, which is what `Product`'s own header promises.
7. **`geodesic_velocity` takes no side.** GE.8 makes it the body twist $d = X_1 \ominus_R X_0$,
   whose rows are `rminus_jacobians`; `API.md` §3 lists it without a parameter. The spatial twist
   is $\mathrm{Ad}_{X_0}d = X_1 \ominus_L X_0$, which is `lminus` and needs no method. The rustdoc
   says so.

## Rationale

The alternative arrangements all cost more than one method:

- **`scale` + `sub`, implementing §10 literally** — two methods, the worse numerics by twelve
  orders of magnitude at a stratum the spec samples, an extra `Ad_{Δ⁻¹}` to form (a group inverse,
  an `adjoint` and a third `mul`), and it writes down the one FMA-contractible shape in the phase.
- **A fused `affine`/`fma` method** — the contractible shape by construction, plus the $0\cdot x$
  problem.
- **`sub_scaled`** — still needs `scale` for $J_1$, or a `Jac::zero` the trait does not have.
- **A per-group required `geodesic_jacobians`** — contradicts `API.md` §3, which lists it as
  *provided*, and forces one hand derivation of §10 per group, every future group included, where
  one provided body serves.
- **Routing through `Side`'s four hooks** (`Sd::plus_jacobians` ∘ `Sd::minus_jacobians` by the
  chain rule) is side-generic for free and reuses each group's closed forms — and reproduces
  *exactly* §10's cancelling difference, because `plus_jacobians`' first row is
  $\mathrm{Ad}_{\mathrm{Exp}\sigma}^{-1}$ and `minus_jacobians`' second is $-J_l^{-1}(d)$. GE.7(a)
  is not expressible through the four hooks. This is the alternative that looks best and is worst.

## Consequences

- `Jac` is a seven-method trait plus `scale`; every implementor gains one line
  (`SEn3Jac`: `diag` and each `col` scaled; `RnJac`: `k·s`; `Mat3`: delegate to the inherent
  `Matrix::scale` so the two can never differ in the bits; `ProductJac`: per block; the test-only
  `HJac`: per block).
- `NUMERICS.md` §10 states both arrangements and which one ships, so a reader who derives the
  difference form is not surprised by the code.
- `S::one() - t` is exact for $t \in [0.5, 1]$ by Sterbenz, which is the whole regime the
  cancellation argument covers, `PHASE4.md` §4's $t = 1 - 10^{-9}$ included; below $0.5$ it rounds
  by at most $u/2$ where $J_0 \to I$. One comment, not a branch.
- `Mat3`'s inherent `Matrix::scale` shadows the trait method on a concrete `Mat3`, as
  `Matrix::identity` already shadows `Jac::identity`. A test asserts the two agree bit for bit.

## Implementation plan

1. `NUMERICS.md` §10's edit and this record — verified by `just lint` and the `decisions/README.md`
   row.
2. `Jac::scale` and its five bodies, plus the `Mat3` shadowing test — verified by `just test` and
   by `laws::jac_dense_order` still exactly `0` on every instantiated group.
3. The provided `geodesic_jacobians`, the `SEn3`/`SO3`/`Product` overrides, and decision 1's
   commuting step checked numerically — verified by `geodesic_jacobians_match_reference` (`Dual`
   through `reference::geodesic`) and `jacobians_match_dual_geodesic` (`PHASE4.md` §2), each bound
   twice the measured worst as `laws::Bounds` requires.
4. The `t = 0` reading — verified by a test: GE.7(a) computes $J_l(d)J_l^{-1}(d)$ there, a matrix
   times its own closed-form inverse, so $J_0 = I + O(\kappa u)$ rather than exactly $I$. The
   bound comes from the measurement, and if it is large near $\pi$ the fallback is a `branch` on
   $t$, which would reinstate a `Jac` difference — so this is measured before the bound is written,
   not after.

## Open questions

None.

## Further work

1. `PHASE4.md` §4 has no corpus id for the geodesic **Jacobians**, only for the geodesic. Whether
   the `Dual` twin plus the laws is enough, or they owe a stratum of their own, is a `PHASE4.md`
   §4 question this record does not depend on.
2. GE.7(c) notes `Product`'s $\mathbb R^3$ blocks are exactly $(1-t)I$ and $tI$. Whether `RnJac`
   should short-circuit the two `mul`s is a measurement, not a correctness question.
3. For $t \notin [0, 1]$ the formulas hold in exact arithmetic (extrapolation along the same screw)
   and their error grows with $\lvert t\rvert$ (GE.13(a)). Nothing in the API restricts $t$, and
   whether it should is `API.md` §6's question, not this record's.
