# helicoid — Phase 2 Implementation Specification: `helicoid-linalg`

> **Companions:** `docs/PROJECT.md` (D4, D10, D16), `docs/API.md` §2,
> [`0003`](./decisions/0003-the-scalar-that-cannot-say-less-than.md) (the scalar model),
> [`0007`](./decisions/0007-the-budget-a-foundation-can-afford.md) (the budget).

**Deliverable:** the leaf crate every other piece stands on — the scalar model, forward-mode dual
numbers, fixed-size linear algebra, strided views, and the three 3×3 decompositions omnisac already
depends on — verified by the Phase 1 instrument. Sections marked **NORMATIVE** are requirements;
signatures in code blocks are normative.

## 0.0 Implementation status

| Area | Status |
|---|---|
| `Mask`, `Real`, `Blend`, `Precision`; `f64`/`f32` impls (§2) | Partial: all four traits and both scalar impls done, `compile_fail` doctests for `<` and `==` in place; `Real::sqrt` (`f64`, `f32`) is checked bit for bit against an integer square root over `2 x 10^6` cases with a pinned digest, with `libm`'s `arch` feature on (`0018`); checked by hand on x86_64, aarch64 (`qemu-user`) and wasm32 (node WASI), run by CI on x86_64 and aarch64; `Dual` is a `Real` and a `Blend<S>` (§3); `Blend` for `Vector`/`Matrix`/`Point` done (§4); the `Blend` impls for `helicoid` value types land with their types |
| `Dual<S, N>` (§3) | Partial: type, `constant`/`variable`, `Real` with every §3 rule, `Blend<S>`, `dual_value_is_plain_value` over every `Real` method (`f64`, `f32`, nested, poisoned derivatives; bitwise up to NaN sign and payload of arithmetic outputs), second order through nesting on `sin_cos` and `atan2`, `sqrt` at 0 tested; `copysign` takes `sgn(s)` from the sign bit; the `atan2`, `sqrt` and quotient derivative domains are documented and pinned; `dual_matches_mpmath_derivative` runs on an inline mpmath fixture (`sqrt`, `sin_cos`, `atan2`, quotient, product), not yet on the corpus ids `real_*` (§8, needs the Phase 1 generator) |
| `Vector`, `Matrix`, `Point`, `hat`/`vee`, `Mat3::inverse_adj`, `chol` (§4) | Partial: `Vector`, `Point`, `Matrix`, the aliases, every listed operation, `Blend`, `hat`/`vee`, `Mat3::inverse_adj` done; algebra proptests under `f64`, `f32`, `Dual<f64, 2>` with derived and measured bounds (10^6 cases, seeded recipe in the test header); summation order and `-0` pinned to the bit; no `PartialEq` and no `Point + Point` pinned by `compile_fail` doctests; the magnitude range of `norm` and `inverse_adj` documented and pinned; `chol`, `solve_lower`, `solve_upper` done (Cholesky-Crout, `NUMERICS.md` §15; mask `0 < pivot` and every entry finite; a failed pivot gives `L_jj = 1` and a zero column, an overflowing entry is stored as `+0`, so `L` is finite for every input; summation order pinned to the bit; proptests at `N = 1..=6` under `f64`, `f32`, `Dual<f64, 2>` against Higham's Thm 8.5, 10.3, 10.4 bounds, rank-deficient, near-singular and wide-dynamic-range inputs, the mask at a zero pivot, the solves' `debug_assert!` under `should_panic`; no corpus stratum or reference twin yet, §8); `Mat2` adjugate not started |
| `Strided`, `StridedMut` (§5) | Done: `col_major`, `row_major`, `with_strides`, `block`, `get`, `rows`, `cols`, `set`; out-of-bounds access panics in release (saturating index, never wraps); the constructors `debug_assert!` the fit; tested against the index formula, faer/Ceres layouts and an exact `u128` model over both view types; `write_dense` (Phase 3) is the first consumer |
| `eig3`, `svd3`, `solve_cubic` + corpus ids (§6) | Not started |
| `mint` feature (§7) | Done: optional feature `mint` (`mint` >= 0.5.7, no default features), `From`/`Into` both ways for `Vector<S, 2..=4>`, `Point<S, 2..=3>` and `Matrix<S, N, N>` at `N` = 2..=4 (`ColumnMatrixN`, field `x` is column 0), `S = f32, f64`; bitwise round trips (`-0`, infinities, subnormals, signalling and payload NaNs pinned in every slot, plus random bit patterns) and component and column order tested; `mint` has no `Point4`, so `Point` stops at 3; `just lint test` cover the feature, `just msrv no-std wasm` build it |
| omnisac migration (§9) | Not started |

## 0. Non-goals and guardrails — read first

**NORMATIVE.** No quaternion, no group, no Lie anything (those are `helicoid`). No heap, no `dyn`,
no `std`. **No `PartialOrd`, `PartialEq` or `bool`-returning comparison on `Real`.** No `mul_add`.
No generic matrix inverse beyond `Mat2`/`Mat3` adjugates and fixed Cholesky. Dependencies: `libm`,
optional `mint`.

## 1. Crate layout

`real.rs` (`Mask`, `Real`, `Blend`, `Precision`), `float.rs` (`f64`/`f32` impls), `dual.rs`,
`vector.rs`, `matrix.rs`, `point.rs`, `skew.rs` (`hat`, `vee`), `strided.rs`, `chol.rs`, `eig3.rs`,
`svd3.rs`, `cubic.rs`, `mint.rs` (feature-gated).

## 2. The scalar model

**NORMATIVE** ([`0003`](./decisions/0003-the-scalar-that-cannot-say-less-than.md)).

```rust
pub trait Mask: Copy {
    fn and(self, o: Self) -> Self;
    fn or(self, o: Self) -> Self;
    fn not(self) -> Self;
    /// Tests and `debug_assert!` only.
    fn all(self) -> bool;
    fn any(self) -> bool;
    /// The branching policy belongs to the mask: `bool` evaluates one arm; a lane mask
    /// evaluates both and blends.
    fn decide<T>(self, t: impl FnOnce() -> T, f: impl FnOnce() -> T,
                 blend: impl FnOnce(Self, T, T) -> T) -> T;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Precision { F32, F64 }

pub trait Real:
    Copy + 'static
    + core::ops::Add<Output = Self> + core::ops::Sub<Output = Self>
    + core::ops::Mul<Output = Self> + core::ops::Div<Output = Self>
    + core::ops::Neg<Output = Self>
{
    type Mask: Mask;
    const PRECISION: Precision;
    /// Exactly representable constants only (0.5, 2, 3, …). Anything that must be correctly
    /// rounded per precision is generated (`coeffs::generated`).
    fn lit(x: f64) -> Self;
    fn zero() -> Self;
    fn one() -> Self;
    fn lt(self, rhs: Self) -> Self::Mask;
    fn le(self, rhs: Self) -> Self::Mask;
    fn select(m: Self::Mask, t: Self, f: Self) -> Self;
    fn branch<T: Blend<Self>>(m: Self::Mask, t: impl FnOnce() -> T, f: impl FnOnce() -> T) -> T {
        m.decide(t, f, |m, a, b| T::blend(m, a, b))
    }
    fn sqrt(self) -> Self;
    fn sin_cos(self) -> (Self, Self);
    fn atan2(self, x: Self) -> Self;
    fn abs(self) -> Self;
    fn copysign(self, sign: Self) -> Self;
    /// The value part as f64. Tests and `debug_assert!` only (API.md R4).
    fn value_f64(self) -> f64;
}

pub trait Blend<S: Real>: Sized {
    fn blend(m: S::Mask, t: Self, f: Self) -> Self;
}
```

- **`f64`/`f32` impls** route every transcendental through `libm` (`sqrt`, `sincos`, `atan2`,
  `fabs`, `copysign`), never `std` (D16). `Mask = bool`. **`impl Mask for bool` is the only place in
  the workspace where a float comparison's result reaches an `if`.**
- **`Blend`** is implemented for `S`, tuples up to arity 8, `[T; N]` where `T: Blend<S>`, `Vector`,
  `Matrix`, `Point`, and (in `helicoid`) every value type.
- **The safe-argument pattern** is how the exact arm stays finite in every lane:
  `let t2s = S::select(small, S::one(), theta2); let theta = t2s.sqrt();`.
- NaN and infinity: inputs are the caller's responsibility; no routine produces a non-finite value
  from finite in-domain input (Phase 1's `nonfinite` column checks it).

## 3. Dual numbers

**NORMATIVE.**

```rust
#[derive(Clone, Copy, Debug)]
pub struct Dual<S, const N: usize> { pub v: S, pub d: [S; N] }
impl<S: Real, const N: usize> Dual<S, N> {
    pub fn constant(v: S) -> Self;
    pub fn variable(v: S, i: usize) -> Self;   // d = e_i; debug_assert!(i < N)
}
impl<S: Real, const N: usize> Real for Dual<S, N> { type Mask = S::Mask; /* … */ }
```

- `PRECISION = S::PRECISION`; comparisons and masks act on the value part only.
- Rules: `a / b` → $(a_d - q\,b_d)/b_v$ with $q = a_v/b_v$ (the quotient rule with no $b_v^2$ to
  overflow); `sqrt` → $d/(2\sqrt v)$; `sin_cos` → $(d\cos v, -d\sin v)$; `atan2(y, x)` →
  $(x_v y_d - y_v x_d)/(x_v^2 + y_v^2)$; `abs` → $\mathrm{sgn}(v)\,d$ with $\mathrm{sgn}(\pm 0) = +1$;
  `copysign(x, s)` → $\mathrm{sgn}(x_v)\,\mathrm{sgn}(s_v)\,x_d$, with $\mathrm{sgn}(\pm 0) = +1$
  for $x_v$ and $\mathrm{sgn}(s_v) = -1$ iff the **sign bit** of $s_v$ is set, so
  $\mathrm{copysign}(x, -0.0) = -\lvert x\rvert$ and differentiates as such: the derivative of
  the value returned.
- **The value path of `Dual<S, N>` is bitwise identical to the plain `S` evaluation** of the same
  generic code (`dual_value_is_plain_value`), except the sign and payload of a NaN produced by
  arithmetic or a `libm` call: Rust leaves them unspecified and a release build may commute
  operands. Sign-bit operations and `select` are compared exactly, NaN included. This is what
  makes "the derivative of the shipped code" true.
- `sqrt` at 0 has an infinite derivative by the rule above (NaN, $0/0$, in a component whose $d$
  is zero, untouched components included); the safe-argument pattern keeps it out of any selected
  arm, and Phase 1's seeded defect proves the harness notices when it does not.
- **Derivative domains** (`# Domain` on each method, `NUMERICS.md` §12), none asserted: a check
  would panic on inputs the plain value path accepts. `atan2` divides by $x_v^2 + y_v^2$ and is
  accurate only while that sum is normal (the larger argument in about $10^{\pm154}$ for `f64`,
  $10^{\pm19}$ for `f32`), and `a / b` is NaN once $q$ overflows and loses relative accuracy while
  $q$ is subnormal; the value is unaffected in every case. Two tests pin them
  (`dual_atan2_derivative_is_accurate_to_the_edge_of_its_domain_and_only_there`,
  `dual_quotient_derivative_is_nan_once_the_quotient_overflows`).
- Nesting (`Dual<Dual<f64, M>, N>`) is supported and tested to second order on `sin_cos`, `atan2`.

## 4. Fixed-size types

**NORMATIVE.**

```rust
#[repr(C)] #[derive(Clone, Copy, Debug)] pub struct Vector<S, const N: usize>(pub [S; N]);
#[repr(C)] #[derive(Clone, Copy, Debug)] pub struct Point<S, const N: usize>(pub [S; N]);
/// Column-major: `cols[c][r]`.
#[repr(C)] #[derive(Clone, Copy, Debug)] pub struct Matrix<S, const R: usize, const C: usize> { cols: [[S; R]; C] }
pub type Vec2<S> = Vector<S, 2>;  pub type Vec3<S> = Vector<S, 3>;
pub type Mat2<S> = Matrix<S, 2, 2>; pub type Mat3<S> = Matrix<S, 3, 3>;
pub type Point2<S> = Point<S, 2>; pub type Point3<S> = Point<S, 3>;
```

- Operations: `+ - neg`, scalar `scale`, `dot`, `cross` (3), `norm_sq`, `norm`, matrix×matrix
  (generic `R, K, C`), matrix×vector, `transpose`, `identity`, `from_cols`, `from_rows`, `col`,
  `row`, `get(r, c)`, `set(r, c, v)`. `Point − Point = Vector`, `Point + Vector = Point`; no
  `Point + Point`.
- `hat(v) -> Mat3`, `vee(m) -> Vec3` (skew part, no symmetrization check).
- `Mat3::inverse_adj(&self) -> (Mat3, S)` returns adjugate/det and det; the caller decides.
- `chol<S, const N: usize>(a: &Matrix<S, N, N>) -> (Matrix<S, N, N>, S::Mask)`: lower-triangular
  `L` and the positive-definiteness mask; `solve_lower`, `solve_upper` beside it.
- **No `PartialEq` on these types** (it would be a float comparison, R4); tests compare through
  `value_f64`.

## 5. Strided views

**NORMATIVE.**

```rust
pub struct StridedMut<'a, S> { data: &'a mut [S], rows: usize, cols: usize, rs: usize, cs: usize }
impl<'a, S: Copy> StridedMut<'a, S> {
    pub fn col_major(data: &'a mut [S], rows: usize, cols: usize) -> Self;
    pub fn row_major(data: &'a mut [S], rows: usize, cols: usize) -> Self;
    pub fn with_strides(data: &'a mut [S], rows: usize, cols: usize, rs: usize, cs: usize) -> Self;
    pub fn block(&mut self, r0: usize, c0: usize, rows: usize, cols: usize) -> StridedMut<'_, S>;
    pub fn set(&mut self, r: usize, c: usize, v: S);
}
```

`Strided<'a, S>` is the read-only twin. Bounds are checked by safe slice indexing, in release too, on
`get`, `set` and `block`: an out-of-bounds access panics and is documented under `# Panics` (D11's
one class). A view that does not fit its slice is a `debug_assert!` at construction and a panic on
the first access that reaches past the slice. Strides are non-negative `usize`. Solvers hand `helicoid` a block of
their Jacobian through these (faer `MatMut` and Ceres row-major buffers both map onto
`with_strides`).

## 6. `eig3`, `svd3`, `solve_cubic`

**NORMATIVE.**

- `eig3(a: &Mat3<S>) -> (Vec3<S>, Mat3<S>)`: symmetric input; eigenvalues ascending, orthonormal
  eigenvectors as columns with $\det = +1$. **Port omnisac's closed form (Smith) as is.** Bars, per
  stratum, against `mp.eigsy`: eigenvalue error in $u\|A\|$; eigenvector angular error in
  $u\|A\|/\mathrm{gap}$ (the Davis–Kahan scale). Strata: `eig:gap-1e-k` ($k = 0, \dots, 12$, two
  eigenvalues at relative gap $10^{-k}$), `eig:triple`, `eig:rank1`, `eig:random`. **If the
  eigenvector bar fails on a clustering stratum**, add Kopp's hybrid (analytic with an iterative
  fallback selected by mask) and record which strata switched; do not widen the bar.
- `svd3(a: &Mat3<S>) -> (Mat3<S>, Vec3<S>, Mat3<S>)`: McAdams et al., fixed Jacobi sweep count
  (no convergence loop), signed so $\det U = \det V = +1$ and only $\sigma_3$ may be negative
  (`nearest_rotation` is $U V^\top$ directly). Against `mp.svd_r`.
- `solve_cubic`: **omnisac's current signature, unchanged**, so its migration is a pure move.
  Against `mp.polyroots`; strata: distinct roots, near-double, triple, one real root.
- All three branch only through `S::branch`/`S::select`.

## 7. `mint`

`From`/`Into` between `Vector<S, N>`/`Point<S, N>` (`N` = 2, 3, 4) and `mint::VectorN`/`PointN`
(`mint` defines no `Point4`, so `Point` converts at `N` = 2, 3), and between `Matrix<S, N, N>` and
`mint::ColumnMatrixN`, for `S = f32, f64`. Nothing else.

## 8. Tests

- Proptests for every operation's algebra (associativity to a recorded bound, transpose/product
  identities), under `f64`, `f32` and `Dual<f64, 2>`.
- `dual_value_is_plain_value` over every `Real` method; `dual_matches_mpmath_derivative` for
  `sqrt`, `sin_cos`, `atan2` on corpus ids `real_*` (Phase 1 generator, added here).
- `eig3`/`svd3`/`solve_cubic` as conformance subjects over their strata; envelope rows include
  nalgebra's `SymmetricEigen`/`SVD` through a runner (reported, not gated).
- `just no-std` builds the crate with no `alloc`.

## 9. The omnisac migration

omnisac never ships its own numeric crate: its primitives depend on `helicoid-linalg`.
**Gate:** omnisac's RunRecords are bit-identical where omnisac already routed through `libm`;
otherwise discrete outputs (inlier sets, `residual_evals`, iteration counts) are identical and value
differences are recorded per stratum in the omnisac PR. omnisac's LM core stays in omnisac
([`0009`](./decisions/0009-what-helicoid-does-not-own.md)).

## 10. Definition of done

- [ ] §2–§7 implemented; `#![forbid(unsafe_code)]`, `no_std`, no `alloc`; dependency closure is
      `{libm}` (+ `mint` with the feature).
- [ ] `Real` has no `PartialOrd`/`PartialEq` — a `compile_fail` doctest asserts `a < b` does not
      compile for `S: Real`.
- [ ] `dual_value_is_plain_value` green for every method.
- [ ] `eig3`, `svd3`, `solve_cubic` pass their strata; baseline blessed.
- [ ] omnisac depends on `helicoid-linalg`; §9's gate met.

## Appendix: suggested implementation order

1. `Mask`, `Real`, `Blend`, `f64`/`f32` impls, the `compile_fail` doctest.
2. `Dual`, its rules, `dual_value_is_plain_value`.
3. `Vector`, `Point`, `Matrix`, `hat`/`vee`, adjugate inverses, `chol`.
4. Strided views.
5. `solve_cubic`, `eig3`, `svd3` ports and their corpus ids and strata.
6. `mint`.
7. The omnisac migration PR (in omnisac).
