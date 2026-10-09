# helicoid — the API contract

> **Companions:** [`PROJECT.md`](./PROJECT.md) (decision log D1–D18),
> [`NUMERICS.md`](./NUMERICS.md) (formulas), [`PHASE1.md`](./PHASE1.md)–[`PHASE6.md`](./PHASE6.md)
> (per-phase specs), [`PHASE7.md`](./PHASE7.md) (gated).

The rules that generate every public item. Sections marked **NORMATIVE** are requirements; item
lists are the target surface, delivered by the phase named in each row. **Status:** ready; this
document schedules nothing.

## 1. The six rules

A question none of these answers is a decision record, not an API choice.

### R1 — Spell the side

**NORMATIVE.** Every perturbation operation names its side: `rplus`, `lplus`, `rminus`, `lminus`,
and their `*_jacobians`. Generic code takes `Sd: Side` (`Left`, `Right`). **No group, tangent or
manifold type implements `Add`, `Sub`, `AddAssign` or `SubAssign`.** A group type exposes a field
whose type does (`Rn`'s `Vector`) only where the group is **abelian**, so that `⊕_R = ⊕_L` and the
reach-through can name no side; a tangent is a vector space and is exempt
([`0025`](./decisions/0025-a-structured-jacobian-and-a-sealed-side.md)). `Mul` is composition
(`a * b` is `T_a_x · T_x_b`) and action (`x * p`). Each side's Jacobians are expressed in that
side's perturbation convention (`NUMERICS.md` §2.3).

### R2 — Values in, values out

**NORMATIVE.** Groups, tangents, Jacobians, charts and `Gaussian` are `Copy`. No method allocates,
locks, reads global state or mutates hidden state; the only `&mut self` numeric method is
`renormalize`. Batch operations write into caller memory (`act_many(&self, &[Point3<S>],
&mut [Point3<S>])`); there is no `Vec` anywhere.

### R3 — Order is stated, never inferred

**NORMATIVE.** Tangents have named fields (`phi`, `rho`, `theta`, `sigma`). A flat array appears
only through `write_dense`/`read_dense`, whose order is `NUMERICS.md` §1's (rotation first), and a
foreign order only through a converter named after that order: `from_translation_first`,
`to_translation_first`, `Quat::from_xyzw`, `Quat::to_xyzw`, `Quat::from_jpl`. There is no
`From<[S; N]>`, no `Into<[S; N]>`, and no `Index` on a tangent.

### R4 — No boolean from a float

**NORMATIVE.** `Real` has no `PartialOrd`/`PartialEq`. Comparisons return `S::Mask`; control flow
goes through `S::branch`/`S::select`. A public function whose answer depends on a float comparison
returns `S::Mask` (e.g. `chol` returns `(L, S::Mask)` for positive-definiteness). `value_f64` exists
for tests and `debug_assert!` only.

### R5 — Jacobians are typed; dense is a write

**NORMATIVE.** A group's Jacobian is `G::Jac`, a structured type closed under `mul` and `inverse`
([`0005`](./decisions/0005-the-jacobian-is-a-dual-matrix.md)). A dense matrix is produced only by
writing into caller memory (`write_dense(&mut StridedMut)`) or by `sandwich` into a fixed-size
covariance. Ambient Jacobians (Phase 6) follow the same rule.

### R6 — The domain is part of the contract

**NORMATIVE.** Every function with a restricted domain carries a `# Domain` rustdoc section naming
it (`NUMERICS.md` §12) and a `debug_assert!` enforcing it. Release builds never check; out-of-domain
release behaviour is an unspecified value, never a panic (bar out-of-bounds strided access, D11),
never UB. The functions that assert nothing are `chol`, `solve_cubic` and `eig3`: each accepts every
input; the first two report through a mask (R4), `eig3` reads one triangle and reports by non-finite
eigenvalues (`0023` (draft)).
Unit-norm inputs are checked by `*_unchecked` constructors in debug only; `*_normalized`
constructors normalize.

## 2. `helicoid-linalg` — the leaf (Phase 2)

| Item | Kind | Notes |
|---|---|---|
| `Mask` | trait | `and`, `or`, `not`, `all`, `any`, `decide`. `impl Mask for bool`. |
| `Real` | trait | `PRECISION`, `lit`, `zero`, `one`, `lt`, `le`, `select`, `branch` (provided), `sqrt`, `cbrt`, `sin_cos`, `sin`, `cos`, `acos`, `atan2`, `abs`, `copysign`, `value_f64`. `impl` for `f64`, `f32`, `Dual<S, N>`. |
| `Blend<S>` | trait | Lane-wise select for tuples, arrays, `Vector`, `Matrix` and every `helicoid` value type. |
| `Precision` | enum | `F32`, `F64`; selects generated constants. |
| `Dual<S, const N: usize>` | struct | `{ v: S, d: [S; N] }`; `variable(v, i)`, `constant(v)`; nests. |
| `Vector<S, N>`, `Matrix<S, R, C>`, `Point<S, N>` | structs | `repr(C)`; `Matrix` column-major `[[S; R]; C]`; aliases `Vec2`, `Vec3`, `Mat2`, `Mat3`, `Point2`, `Point3`. |
| `hat`, `vee` | fns | `Vec3 ↔ Mat3` skew. |
| `Mat2::inverse_adj`, `Mat3::inverse_adj` | fns | `(adjugate/det, det)`; the caller decides what `det` means. `Mat2`'s adjugate is exact, so `det` is its only rounding ([`0061`](./decisions/0061-mat2-keeps-its-adjugate.md)). |
| `chol<S, N>`, `solve_lower`, `solve_upper`, `chol_solve` | fns | Fixed-size Cholesky; `(L, S::Mask)`. Triangular solves reading one triangle each; `chol_solve(&L, b)` solves `A x = b` from the factor without a transpose ([`0019`](./decisions/0019-a-cholesky-solve-without-the-transpose.md)); `NUMERICS.md` §15. |
| `Strided<'a, S>`, `StridedMut<'a, S>` | structs | `col_major`, `row_major`, `block`; bounds-checked (`# Panics`). |
| `eig3`, `svd3`, `solve_cubic` | fns | `NUMERICS.md` §13 references; signatures in `PHASE2.md` §6. |
| `mint` | feature | `From`/`Into` for `Vector`, `Matrix` at sizes 2–4 and `Point` at 2–3 (`mint` has no 4-D point). |

## 3. `helicoid` — groups and geometry

| Item | Kind | Phase | Notes |
|---|---|---|---|
| `Tangent<S>`, `Jac<S, T>`, `LieGroup<S>`, `Side`, `Left`, `Right` | traits/ZSTs | 3 | Signatures `PHASE3.md` §2. `Side` is sealed: `Left` and `Right` are the only impls. `Tangent::dot_acc` is required, `dot` provided (`0025`). |
| `Quat<S>` | struct | 3 | `{ w, x, y, z }`; `from_wxyz_unchecked`, `from_wxyz_normalized`, `from_xyzw`, `to_xyzw`, `from_jpl`, `conjugate`, `norm_sq`, `to_matrix`, `renormalize` ([`0027`](./decisions/0027-a-normalizing-constructor-normalizes.md)). No `inverse`: `conjugate` is it on the unit sphere and nothing covers the rest. No `norm` or `dot` ([`0048`](./decisions/0048-the-relative-transform-pair-earns-the-surface-dot-and-norm-do-not.md)): the fields are `pub` and `norm_sq` is public, so each is one line in a consumer, and the angle between two rotations is `rminus`'s tangent norm through D5's `atan2`, not a cosine. |
| `SO2<S>`, `SO3<S>`, `SE2<S>` | structs | 3 | |
| `SEn3<S, const N: usize>`; `SE3<S> = SEn3<S, 1>`, `SE23<S> = SEn3<S, 2>` | struct, aliases | 3 | `{ q: Quat<S>, x: [Vec3<S>; N] }` (private); `from_parts`, `parts`, `rotation`, `renormalize` ([`0044`](./decisions/0044-four-primitives-the-first-consumer-names-and-no-spec-does.md)), `mul_inv`, `inv_mul` (Phase 4, [`0048`](./decisions/0048-the-relative-transform-pair-earns-the-surface-dot-and-norm-do-not.md)), and at `N = 1` `translation`, `from_rt`, `from_quat_translation`, `act_many`, `act_jacobians` and `Mul<Point3>` (the single-point action has no named method: the operator *is* it, R1); at `N = 2` `velocity`, `position`. No `normalized(self) -> Self`: R3 keeps one spelling per operation and R2 makes `&mut self` the shape for the in-place one. The composition `a * b.inverse()` stays available by writing it — it is `mul_inv`'s §14 reference twin — and gets no named spelling of its own. |
| `SEn3Tangent<S, N>`; `Twist<S> = SEn3Tangent<S, 1>` | struct, alias | 3 | `{ phi, rho: [Vec3; N] }`; `Twist::omega()`, `Twist::v()`. |
| `SEn3Jac<S, N>` | struct | 3 | `{ diag: Mat3, col: [Mat3; N] }` ([`0005`](./decisions/0005-the-jacobian-is-a-dual-matrix.md)). |
| `Rn<S, N>`, `RnTangent<S, N>`, `RnJac<S, N>`, `Product<A, B>`, `ProductJac<A, B>` | structs | 3 | Block-diagonal Jacobians; `RnJac` is `k I`. |
| `act`, `act_many`, `act_jacobians` | methods | 3 | SO(2), SO(3), SE(2), SE(3), Sim(3) only. |
| `from_matrix` | method | 3 | SO(3): Shepperd, closed form, never iterative. |
| `LieGroup::geodesic`, `geodesic_jacobians`, `geodesic_velocity` | provided methods | 4 | Reference twin = definition. |
| `Chart<S, M>`, `Manifold<S>`, `WithChart<M, C>`, `Lifted<C>` | traits/structs | 5 | [`0012`](./decisions/0012-a-retraction-is-a-chart.md). One `Manifold` impl per group, no blanket; `WithChart`'s chart is `Lifted<C>` ([`0060`](./decisions/0060-the-charts-go-first-and-name-their-frame.md)). |
| `RightChart<G>`, `LeftChart<G>`; `Screw<S> = RightChart<SE3<S>>`, `Decoupled<S>`, `WorldTranslation<S>` | charts | 5 | Each holds its frozen base. The three SE(3) charts share the tangent `Twist<S>` (`0060`). |
| `TwistBlockJac<S>` | struct | 5 | `Decoupled`'s and `WorldTranslation`'s `Jac`: a newtype over `ProductJac<Mat3, Mat3>` implementing `Jac<S, Twist<S>>` only; `rotation_block`, `translation_block` (`0060`). |
| `Se3Chart<S>`, `SE3::chart_transition::<From, To>` | sealed trait, method | 5 | `DΦ(0)` between two SE(3) charts at a base: `I`, `diag(I, R)` or `diag(I, Rᵀ)` (`0060`). |
| `S2<S>`, `S2Chart<S>` | structs | 5 | Frozen Householder basis. `S2::{from_vec_unchecked, from_vec_normalized, renormalize, vec}` ([`0063`](./decisions/0063-the-sphere-reads-its-sign-by-comparison-and-is-held-unit.md)). |
| `sinc`, `sinc_value` | fns | 5 | `(sin θ/θ, its derivative in θ²)` at `θ²`, and the value alone on plain `S`, evaluated in `coeffs` ([`0062`](./decisions/0062-sin-theta-over-theta-is-public-and-differentiates-its-branch.md)). |
| `Sim3<S>` | struct | 5 | Formula block owed (`NUMERICS.md` §9). |
| `so3::gamma1`, `so3::gamma2`, `so3::gamma_apply_jacobian` | fns | 5 | `pub mod so3`. `Γ₁ = SO3::jl` to the bit; `Γ₂ = ½I + bW + dW²`; `gamma_apply_jacobian::<M, S>(&φ, v) -> (Γ_M v, ∂(Γ_M v)/∂φ)`, `M ∈ {1, 2}`, through `Dual<S, 3>` ([`0064`](./decisions/0064-the-integrated-exponentials-reuse-the-swept-switches.md)). |
| `Gaussian<S, G, Sd, const D: usize>` | struct | 5 | `D == G::DOF` asserted at compile time. |
| `AmbientChart` | trait | 6 | Ceres-style `PlusJacobian`/`MinusJacobian`. |
| `reference::*` | module | 3–5 | Public and documented: the twins are the definition of *correct* (D6). |

`coeffs` is `pub(crate)`, except `sinc` and `sinc_value` (`0062`). Its `__sweep` feature exposes evaluators to `xtask` only and is not part of
the API ([`0004`](./decisions/0004-switch-points-are-generated-not-typed.md)); `just lint` fails if
any crate other than `xtask` enables it.

## 4. Interop

- **`mint`** (optional, both crates): `Vector`, `Matrix`, `Point`; `Quat ↔ mint::Quaternion`
  (`{ v, s }` — field names, no order ambiguity). nalgebra, glam and cgmath convert through `mint`;
  there is no `nalgebra` or `faer` feature ([`0007`](./decisions/0007-the-budget-a-foundation-can-afford.md)).
- **Dense order converters** (R3): the only way a translation-first or `xyzw` array enters or
  leaves.
- **Solvers** consume `Manifold`, `Chart`, `Jac::write_dense` and `AmbientChart`; `helicoid` adds
  nothing solver-specific ([`0009`](./decisions/0009-what-helicoid-does-not-own.md)).

## 5. Stability

`0.0.x` until [`PHASE6.md`](./PHASE6.md) §6 holds: every release may break every other. From 1.0:
`cargo-semver-checks` gates every release. **A change to a generated switch point or series length
changes outputs without changing the API**: it ships in a minor release with a changelog line
naming the functions, strata and the old/new max error.

## 6. The check a new surface has to pass

A PR adding public API answers each line in its description:

1. Which rule of §1 could it violate, and why does it not?
2. Does it have a corpus stratum, or is it a composite with a reference twin (`NUMERICS.md` §14)?
3. Is its domain stated (`# Domain`) and `debug_assert!`ed?
4. Does it return a `bool` from a float, a dense Jacobian, or an unlabeled array?
5. Does it duplicate an existing path (a second spelling)? Document the one that exists instead.
6. Does it pass the ownership test ([`0009`](./decisions/0009-what-helicoid-does-not-own.md))?
7. Is it `no_std`, allocation-free, dependency-free, and bit-identical across targets?
8. Does it **remove a way to be wrong** that the consumer cannot see — a tangent order, a side, a
   chart, a switch point, a series, a branch cut, a quaternion sign — or does it only relocate code
   the consumer already has right? ([`0049`](./decisions/0049-the-boundary-is-what-removes-a-way-to-be-wrong.md); question 6 asks whether `helicoid` *may* own it, which is
   a different thing, and three items have now shipped past it.)
