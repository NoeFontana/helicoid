# helicoid — Phase 2 Implementation Specification: `helicoid-linalg`

> **Companions:** `docs/PROJECT.md` (D4, D10, D16), `docs/API.md` §2,
> [`0003`](./decisions/0003-the-scalar-that-cannot-say-less-than.md) (the scalar model),
> [`0007`](./decisions/0007-the-budget-a-foundation-can-afford.md) (the budget).

**Deliverable:** the leaf crate every other piece stands on — the scalar model, forward-mode dual
numbers, fixed-size linear algebra, strided views, and core 3×3 decompositions — verified by the Phase 1 instrument. Sections marked **NORMATIVE** are requirements;
signatures in code blocks are normative.

## 0.0 Implementation status

| Area | Status |
|---|---|
| `Mask`, `Real`, `Blend`, `Precision`; `f64`/`f32` impls (§2) | Partial: traits and scalar impls implemented, including `cbrt` and `sqrt` bit identity with `libm` `arch` (`0017`, `0018`). |
| `Dual<S, N>` (§3) | Partial: dual numbers generic over `Real` with automatic differentiation; corpus ids `real_sqrt`, `real_cbrt`, `real_sin_cos`, `real_acos`, `real_atan2`, `real_div` at both precisions, ≤ 3.65 `u` ([`0056`](./decisions/0056-the-routines-d7-does-not-reach.md)). |
| `Vector`, `Matrix`, `Point`, `hat`/`vee`, `Mat3::inverse_adj`, `chol` (§4) | Partial: fixed-size types, operations, and `chol_solve` implemented (`0019`); corpus ids `chol_n{3,6}`, `chol_solve_n{3,6}` (`0056`). `Mat2` adjugate pending. |
| `Strided`, `StridedMut` (§5) | Done: column-major, row-major, strided slice views. |
| `eig3`, `svd3`, `solve_cubic` + corpus ids (§6) | Partial: `solve_cubic` and `eig3` implemented with `libm` transcendentals (`0017`, `0022`, `0053`). Corpus ids `solve_cubic` and `eig3` at both precisions, `eig3` against nalgebra (`0056`). `svd3` and its id pending. |
| `mint` feature (§7) | Done: optional feature conversions for fixed-size types. |
| Downstream perception migration (§9) | Not started |

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
    fn cbrt(self) -> Self;
    fn sin_cos(self) -> (Self, Self);
    fn sin(self) -> Self;
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

- **`f64`/`f32` impls** route every transcendental through `libm` (`sqrt`, `cbrt`, `sincos`,
  `sin`, `atan2`, `fabs`, `copysign`), never `std` (D16). **`sin` is bit-identical to
  `sin_cos().0`** and that is a requirement on the impl, not an observation: [`0052`](./decisions/0052-real-owes-sin-and-the-corpus-does-not-move.md) adds it for
  the callers that discard the cosine, every recorded figure of a routine re-spelled from one to the
  other depends on the identity, and `libm` is a caret dependency whose two functions are *not* the
  same expression — `sin` and `sincos` cut their small-argument fast path at different thresholds,
  and at binary32 their octant arms are exact negations of each other.
  `sin_tests` enforces it, exhaustively over all `2^32` binary32 patterns. `Mask = bool`. **`impl Mask for bool` is the
  only place in the workspace where a float comparison's result reaches an `if`.**
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
  overflow); `sqrt` → $d/(2\sqrt v)$; `cbrt` → $d/(3c^2)$, $c = \mathrm{cbrt}(v)$; `sin_cos` →
  $(d\cos v, -d\sin v)$ and `sin` → $d\cos v$, taken *through* `sin_cos` because a dual's
  derivative needs the cosine, so the cheaper call saves nothing here and the value is identical by
  construction ([`0052`](./decisions/0052-real-owes-sin-and-the-corpus-does-not-move.md)); `atan2(y, x)` → $(x_v y_d - y_v x_d)/(x_v^2 + y_v^2)$; `abs` →
  $\mathrm{sgn}(v)\,d$ with $\mathrm{sgn}(\pm 0) = +1$;
  `copysign(x, s)` → $\mathrm{sgn}(x_v)\,\mathrm{sgn}(s_v)\,x_d$, with $\mathrm{sgn}(\pm 0) = +1$
  for $x_v$ and $\mathrm{sgn}(s_v) = -1$ iff the **sign bit** of $s_v$ is set, so
  $\mathrm{copysign}(x, -0.0) = -\lvert x\rvert$ and differentiates as such: the derivative of
  the value returned.
- **The value path of `Dual<S, N>` is bitwise identical to the plain `S` evaluation** of the same
  generic code (`dual_value_is_plain_value`), except the sign and payload of a NaN produced by
  arithmetic or a `libm` call: Rust leaves them unspecified and a release build may commute
  operands. Sign-bit operations and `select` are compared exactly, NaN included. This is what
  makes "the derivative of the shipped code" true.
- `sqrt` and `cbrt` at 0 have an infinite derivative by the rules above (NaN, $0/0$, in a
  component whose $d$ is zero, untouched components included); the safe-argument pattern keeps a
  zero out of any selected arm, and Phase 1's seeded defect, on `sqrt`, proves the harness notices
  when it does not ([`0020`](./decisions/0020-dual-sqrt-at-zero-keeps-its-derivative.md)).
- **Derivative domains** (`# Domain` on each method, `NUMERICS.md` §12), none asserted: a check
  would panic on inputs the plain value path accepts. `atan2` divides by $x_v^2 + y_v^2$ and is
  accurate only while that sum is normal (the larger argument in about $10^{\pm154}$ for `f64`,
  $10^{\pm19}$ for `f32`), and `a / b` is NaN once $q$ overflows and loses relative accuracy while
  $q$ is subnormal; the value is unaffected in every case. Two tests pin them
  (`dual_atan2_derivative_is_accurate_to_the_edge_of_its_domain_and_only_there`,
  `dual_quotient_derivative_is_nan_once_the_quotient_overflows`).
- Nesting (`Dual<Dual<f64, M>, N>`) is supported and tested to second order on `sin_cos`,
  `atan2`, `cbrt`.

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
- `chol_solve<S, const N: usize>(l: &Matrix<S, N, N>, b: Vector<S, N>) -> Vector<S, N>`: solves
  `A x = b` from the factor `l` of `A`; `solve_lower` then a back substitution against `l^T` read
  by column (`NUMERICS.md` §15.6), bit-identical to `solve_upper(&l.transpose(), solve_lower(&l,
  b))`, its reference twin. The column solve itself is not public
  ([`0019`](./decisions/0019-a-cholesky-solve-without-the-transpose.md)).
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
  eigenvectors as columns with $\det = +1$. Smith's closed form. Bars against `mp.eigsy`: eigenvalue error in $u\|A\|$; eigenvector angular error in
  $u\|A\|/\mathrm{gap}$ (`NUMERICS.md` §11), each the committed per-stratum maximum and no typed
  constant. Strata: `eig:gap-1e-k/{bottom,top}` ($k = 0, \dots, 12$; `@f32` $0, \dots, 6$),
  `eig:triple`, `eig:rank1`, `eig:random`, `eig:scale-{up,down}`
  ([`0056`](./decisions/0056-the-routines-d7-does-not-reach.md)).
- `svd3(a: &Mat3<S>) -> (Mat3<S>, Vec3<S>, Mat3<S>)`: McAdams et al., fixed Jacobi sweep count, signed so $\det U = \det V = +1$.
- `solve_cubic<S: Real>(a: S, b: S, c: S, d: S) -> (Vec3<S>, [S::Mask; 3])`: real roots of
  `a x³ + b x² + c x + d` with validity mask ([`0017`](./decisions/0017-cbrt-and-mask-valued-roots.md)); power-of-two tolerances; `Real::acos` and `Real::cos` via `libm` ([`0022`](./decisions/0022-real-owes-acos-and-cos.md), [`0053`](./decisions/0053-acos-is-better-the-roots-are-not-necessarily.md)). Against
  `mp.polyroots` of the stored coefficients, by root-set distance (`NUMERICS.md` §11). Strata:
  `cubic:distinct`, `cubic:double`, `cubic:triple`, `cubic:one-real` (planted, exact),
  `cubic:near-double-1e-k` ($k = 2, 4, 6, 8$), `cubic:one-real-p-small`,
  `cubic:coeff-scale-{up,down}` (`0056`).
- All three branch only through `S::branch`/`S::select`.

## 7. `mint`

`From`/`Into` conversions between fixed-size types and `mint` equivalents for `S = f32, f64`.

## 8. Tests

- Algebraic proptests for all operations under `f64`, `f32`, and `Dual<f64, 2>`.
- `dual_value_is_plain_value` over all `Real` methods; derivative comparisons against mpmath.
- `eig3`/`svd3`/`solve_cubic` conformance tests against reference strata.
- `#![forbid(unsafe_code)]`, `no_std`, zero allocations.

## 9. Downstream perception migration

Downstream estimation crates migrate from ad-hoc math routines to `helicoid-linalg` primitives, adapting call sites to the masked root and ascending eigenpair formats.

## 10. Definition of done

- [ ] §2–§7 implemented; `#![forbid(unsafe_code)]`, `no_std`, no `alloc`; dependencies: `libm`, optional `mint`.
- [ ] `Real` has no `PartialOrd`/`PartialEq` (`compile_fail` doctests enforced).
- [ ] `dual_value_is_plain_value` passes for all methods.
- [ ] `eig3`, `svd3`, `solve_cubic` meet conformance bars.
- [ ] Downstream crates depend on `helicoid-linalg`.

## Appendix: suggested implementation order

1. `Mask`, `Real`, `Blend`, `f64`/`f32` impls.
2. `Dual`, autodiff rules, `dual_value_is_plain_value`.
3. `Vector`, `Point`, `Matrix`, `hat`/`vee`, adjugates, `chol`.
4. Strided views.
5. `solve_cubic`, `eig3`, `svd3`.
6. `mint` integration.
7. Downstream consumer integration.
