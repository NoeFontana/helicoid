# helicoid — Phase 3 Implementation Specification: Core Groups

> **Companions:** `docs/NUMERICS.md` (every formula this phase implements), `docs/API.md` §1 and §3,
> [`0002`](./decisions/0002-one-convention-for-a-stack-that-already-disagrees.md),
> [`0004`](./decisions/0004-switch-points-are-generated-not-typed.md),
> [`0005`](./decisions/0005-the-jacobian-is-a-dual-matrix.md),
> [`0010`](./decisions/0010-seeded-from-tf-tree-math.md).

**Deliverable:** the coefficient kernel with its generated thresholds; SO(2), SO(3), SE(2),
SE_N(3), Rⁿ and products; both sides; every Jacobian of `NUMERICS.md` §2–§6; reference twins; the
first blessed envelope. SO(3)/SE(3) start as a port of `tf_tree_math` generalized to `S: Real`
(0010). Sections marked **NORMATIVE** are requirements; trait signatures are normative.

## 0.0 Implementation status

| Area | Status |
|---|---|
| Traits: `Tangent`, `Jac`, `LieGroup`, `Side` (§2) | Partial: the four traits and `Left`/`Right` as written in §2, `Side` delegating to `rplus`/`lplus`/`rminus`/`lminus` and their `*_jacobians`; the `DOF` tie (in each impl's `identity` and in every generic law) and `sandwich`'s `D == T::DOF` are `const` assertions, so a mismatch fails `cargo build`/`cargo test`, not `cargo check` (`sandwich`'s pinned by a `compile_fail` doctest); implemented by `Rn` and by a test-only non-abelian group (Heisenberg, `heis_tests.rs`; 2-step nilpotent, so no series term past the first is checked and its `Log` has no branch), both checked by the generic laws of `laws.rs` (axioms, `Exp`/`Log`, `Ad`, `ad` as a bracket with its series against `Ad`, `J_r`, `J_l`, `J_l = Ad J_r`, `⊕`/`⊖` round trips, `Side`, the rows of `NUMERICS.md` §2.3 against `Ad` and `J`, dense order, `sandwich`, `dual_value_is_plain_value`) under `f64`, `f32`, `Dual<f64, 3>`; `Side` has no first-class selector for the side row of `compose_jacobians`/`inverse_jacobian` (the test group compares `TypeId`s; whether `Side` gets a method or the groups use `TypeId` is for the SO(3) PR to decide); the rows are not compared with `Dual` differentiation: `jacobians_match_dual_*` and the twins (§8) not started; `Tangent` and `Jac` are also implemented by `SEn3Tangent` and `SEn3Jac` (§5), which the three-dimensional checks of `laws.rs` do not cover (their own tests do) |
| Coefficient kernel + `__sweep` + generated thresholds (§3) | Not started |
| `Quat`, `SO3` (§4) | Partial: `Quat` in `quat.rs`: `Blend`, Hamilton product as `Mul`, `conjugate`, `norm_sq`, `to_matrix` (`R(q)` of `NUMERICS.md` §1, a scaled rotation for a non-unit `q`), `from_wxyz_unchecked` (`debug_assert!` of the §3.6 bound, by `S::PRECISION`), `from_wxyz_normalized` and `renormalize` (the first-order Newton step, asserting nothing as §3.6 and §12 state no domain; its accuracy domain, `\|η\| <= 2^-26.29` (`f64`) or `2^-11.79` (`f32`), is documented, and outside it the step is not a normalization; `API.md` R6 reads `*_normalized` as normalizing, which §3.6 does not: for the maintainer to reconcile), `from_xyzw`, `to_xyzw`, `from_jpl` (array arguments, unchecked permutations); checked under `f64`, `f32`, `Dual<f64, 4>` (associativity, `q q* = 1`, `R` orthogonal with `det = 1` and equal to `q v q*`, `from_jpl` against the JPL matrix, the `from_wxyz_unchecked` domain on both sides, quadratic convergence of the step, the `Dual` value path bit-identical to `f64`); `f64` and `f32` only: the Hamilton identities, the bit-exact left-to-right goldens of the product and `norm_sq`, the hand cases; the `Dual` derivative lanes are checked for `Mul` and `norm_sq` only; not started: `SO3` and everything after `Quat` in the list of §4 and the `Quat` corpus strata |
| `SEn3<S, N>`, `SEn3Tangent`, `SEn3Jac`; `SE3`, `SE23` (§5) | Partial: `SEn3Tangent` (`sen3.rs`: `Tangent` with `DOF = 3 + 3N`, `dot` one flat sum in index order) with `Twist`, `omega()`, `v()` and `from_translation_first`/`to_translation_first` (the permutation of `docs/maths/se3.md` SE.14(a), exact); `SEn3Jac` (`dualmat.rs`) as a `Jac`: `mul` (27 + 54N multiplications, measured with a counting scalar), `inverse` (`A⁻¹` from `Mat3::inverse_adj`; the computed `det A` nonzero and finite is `debug_assert!`ed, and a release build returns non-finite entries at `det A = 0` where `NUMERICS.md` §12 promises `jr_inv` a finite value; that and `NUMERICS.md` §5.4 and `0005` taking `A⁻¹` from the SO(3) closed form, which `inverse(&self)` cannot receive: for the `jr_inv` PR to settle), `neg`, `apply`, `apply_transpose`, `write_dense`, `sandwich`, `identity`; checked under `f64`, `f32`, `Dual<f64, 3>` for `N = 1, 2, 3` against the twins (§9 row below), `write_dense` against the matrices of `NUMERICS.md` §2.2 (written out for `N = 1, 2`, entries for `N = 3`), `sandwich` and `apply` against dense, exact integer algebra, the `Dual` value path bit-identical to `f64`, `twist_translation_first_round_trip`; not started: `SEn3`, `SE3`, `SE23` and everything of §5 that needs the group (`exp`, `log`, `adjoint`, `ad`, `jr`, `jl`, `jr_inv`, `jl_inv`, the action) and the corpus ids |
| `SO2`, `SE2` (§6) | Not started |
| `Rn`, `Product` (§7) | Partial: `Rn<S, N>` with tangent `RnTangent { rho }` and `RnJac`, a scalar multiple `k I` (`I`, `-I` and `ad = 0`); every §2.3 row for both sides, all `±I` and symmetric; addition is `Mul` and `Rn + Rn` is pinned by a `compile_fail` doctest; `Product<A, B>` (`product.rs`): `LieGroup` with every method, the provided ones included, the factors' method componentwise (so a factor's override of a provided method is used; no factor overrides one yet, so that is untested until SO(3)), tangent the pair `(A::Tangent, B::Tangent)` (`Tangent` for a pair, dense order `A` then `B`, `dot` the sum of the factors' partial sums, not one flat sum, as `Tangent::dot` now states), `DOF = A::DOF + B::DOF` (the `identity` tie plus the factors' own), every §2.3 row per side as the factors' rows side by side, `Mul` and `Blend` componentwise; `ProductJac<JA, JB>` (its parameters are the factors' Jacobian types, `Product<A, B>::Jac = ProductJac<A::Jac, B::Jac>`) as a `Jac`: `mul`, `inverse`, `neg`, `apply`, `apply_transpose` per block, `write_dense` into the two diagonal blocks (the off-diagonal blocks `+0`), `sandwich` from the factors' `apply` alone (`M = J Σ`, then `(J Mᵀ)ᵀ`; the factors' own `sandwich` needs a block of size `A::DOF`, which a generic `D` cannot build); checked by the generic laws of `laws.rs` under `f64`, `f32`, `Dual<f64, D>` on `Product<Rn<3>, Rn<2>>`, `Product<Rn<2>, Heis>`, `Product<Heis, Rn<2>>` and the nested `Product<Product<Heis, Rn<2>>, Heis>` (`D = 8`), plus hand cases (block-diagonal dense layout on a padded view, `Jac::mul` operand order against the dense product on non-commuting blocks, every row against the factors' rows, `sandwich` cross blocks, `Blend`); the only non-abelian factor is the test-only Heisenberg group (2-step nilpotent), so the block structure is exercised on a non-abelian factor with no series term past the first, and a rotation factor (`Product<SO3, Rn<3>>`) only when SO(3) lands; not started: the twin-table row, `jacobians_match_dual_*` and the `Gaussian`/action uses of a product |
| Side Jacobians, action Jacobians (§8) | Not started |
| Reference twins and proptests (§9) | Partial: `reference::sen3jac_mul` and `sen3jac_inverse` (dense product; Gauss–Jordan with row exchange by `select`), writing into caller memory, and `sen3jac_mul_matches_reference`, `sen3jac_inverse_matches_reference` (the inverse against `κ u` of the dense matrix, on entries uniform on `[-1, 1)` with `κ u <= 1e-3`: ill-conditioned inputs are unchecked until the `sen3_jr_inv_n*` corpus ids); every other row of `NUMERICS.md` §14 and the twin-table lint not started |
| Envelope blessed; `docs/evidence/ENVELOPE.md` (§10) | Not started |
| Benches (§11) | Not started |

## 0. Non-goals and guardrails — read first

**NORMATIVE.** Not in this phase: geodesics (Phase 4); charts other than the implicit right/left
perturbation (Phase 5); S², Sim(3), Γ₂, `Gaussian` (Phase 5); ambient Jacobians (Phase 6). No
group method returns a dense Jacobian; no public item returns a `bool` from a float; no
`Add`/`Sub` on any group or tangent. Every formula comes from `NUMERICS.md`; if one is missing,
stop and ask — do not derive it in the PR.

## 1. Module layout

```
crates/helicoid/src/
├── lib.rs              # conventions doc (0002), re-exports
├── traits.rs           # Tangent, Jac, LieGroup
├── side.rs             # Side, Left, Right
├── coeffs/{mod.rs, kernel.rs, generated.rs, sweep.rs}   # sweep.rs behind `__sweep`
├── quat.rs  so2.rs  so3.rs  se2.rs  sen3.rs  rn.rs  product.rs
├── dualmat.rs          # SEn3Jac algebra
└── reference.rs        # twins (NUMERICS.md §14)
```

## 2. Traits

**NORMATIVE.**

```rust
pub trait Tangent<S: Real>: Copy + Blend<S> {
    const DOF: usize;
    fn zero() -> Self;
    fn add(&self, o: &Self) -> Self;
    fn sub(&self, o: &Self) -> Self;
    fn neg(&self) -> Self;
    fn scale(&self, k: S) -> Self;
    fn dot(&self, o: &Self) -> S;
    /// NUMERICS.md §1 order. debug_assert!(out.len() == DOF).
    fn write_dense(&self, out: &mut [S]);
    fn read_dense(src: &[S]) -> Self;
}

pub trait Jac<S: Real, T: Tangent<S>>: Copy + Blend<S> {
    fn identity() -> Self;
    fn mul(&self, o: &Self) -> Self;
    fn inverse(&self) -> Self;            // exact algebraic inverse (0005)
    fn neg(&self) -> Self;
    fn apply(&self, t: &T) -> T;
    fn apply_transpose(&self, t: &T) -> T;
    fn write_dense(&self, out: &mut StridedMut<'_, S>);
    /// `J Σ Jᵀ` for a dense covariance; `const { assert!(D == T::DOF) }`.
    fn sandwich<const D: usize>(&self, cov: &Matrix<S, D, D>) -> Matrix<S, D, D>;
}

pub trait LieGroup<S: Real>: Copy + Blend<S> + core::ops::Mul<Output = Self> {
    type Tangent: Tangent<S>;
    type Jac: Jac<S, Self::Tangent>;
    const DOF: usize;
    fn identity() -> Self;
    fn inverse(&self) -> Self;
    fn exp(tau: &Self::Tangent) -> Self;
    fn log(&self) -> Self::Tangent;
    fn adjoint(&self) -> Self::Jac;
    fn ad(tau: &Self::Tangent) -> Self::Jac;
    fn jr(tau: &Self::Tangent) -> Self::Jac;
    fn jr_inv(tau: &Self::Tangent) -> Self::Jac;

    // Provided; NUMERICS.md §2.3. Overridden only where a closed form is cheaper, and then the
    // provided body becomes the reference twin.
    fn jl(tau: &Self::Tangent) -> Self::Jac { Self::jr(&tau.neg()) }
    fn jl_inv(tau: &Self::Tangent) -> Self::Jac { Self::jr_inv(&tau.neg()) }
    fn rplus(&self, tau: &Self::Tangent) -> Self { *self * Self::exp(tau) }
    fn lplus(&self, tau: &Self::Tangent) -> Self { Self::exp(tau) * *self }
    /// `self ⊖_R base = Log(base⁻¹ · self)`.
    fn rminus(&self, base: &Self) -> Self::Tangent { (base.inverse() * *self).log() }
    /// `self ⊖_L base = Log(self · base⁻¹)`.
    fn lminus(&self, base: &Self) -> Self::Tangent { (*self * base.inverse()).log() }
    fn rplus_jacobians(&self, tau: &Self::Tangent) -> (Self::Jac, Self::Jac);
    fn lplus_jacobians(&self, tau: &Self::Tangent) -> (Self::Jac, Self::Jac);
    fn rminus_jacobians(&self, base: &Self) -> (Self::Jac, Self::Jac);
    fn lminus_jacobians(&self, base: &Self) -> (Self::Jac, Self::Jac);
    fn compose_jacobians<Sd: Side>(&self, rhs: &Self) -> (Self::Jac, Self::Jac);
    fn inverse_jacobian<Sd: Side>(&self) -> Self::Jac;
}

pub trait Side: Copy + 'static {
    fn plus<S: Real, G: LieGroup<S>>(x: &G, tau: &G::Tangent) -> G;
    fn minus<S: Real, G: LieGroup<S>>(y: &G, x: &G) -> G::Tangent;
    fn plus_jacobians<S: Real, G: LieGroup<S>>(x: &G, tau: &G::Tangent) -> (G::Jac, G::Jac);
    fn minus_jacobians<S: Real, G: LieGroup<S>>(y: &G, x: &G) -> (G::Jac, G::Jac);
}
#[derive(Clone, Copy, Debug)] pub struct Right;
#[derive(Clone, Copy, Debug)] pub struct Left;
```

Each `*_jacobians` returns `(∂/∂first, ∂/∂second)` in its side's convention, exactly the rows of
`NUMERICS.md` §2.3. **`DOF` is duplicated on `Tangent` and `LieGroup`** because stable Rust cannot
compute array lengths from associated consts; a `const` assertion ties them in every impl.

## 3. The coefficient kernel

**NORMATIVE** ([`0004`](./decisions/0004-switch-points-are-generated-not-typed.md)).

- `coeffs` is `pub(crate)` and **the only module that evaluates `NUMERICS.md` §4**. Call-site
  groups: `exp_coeffs(θ²) -> (k, cos_half)`, `jr_coeffs(θ²) -> (a, b)`,
  `jr_inv_coeff(θ²) -> c`, `q_coeffs(θ²) -> (b, d, e)`, `log_ratio(n², w) -> r`,
  `se2_coeffs(θ²) -> (α, β)` (`gamma2_coeffs` arrives in Phase 5). Each is **one** `S::branch` over a
  tuple; the exact arm uses the safe argument.
- `generated.rs` holds a `Switch` per coefficient per precision; kernels select by
  `S::PRECISION` (a `const` match, resolved at monomorphization).
- `sweep.rs` (feature `__sweep`) exposes `exact_<name>(θ²)` and `series_<name>(θ², terms)` for
  `xtask thresholds`. Not API; `just lint` fails if a crate other than `xtask` enables it.
- **Run the sweep** (Phase 1 §6) and commit `generated.rs` and `conformance/sweeps/thresholds.csv`
  in the PR that lands the kernel. The PR description quotes the generated switch points beside the
  `tf_tree` D12 prior.

## 4. `Quat` and SO(3)

**NORMATIVE.**

```rust
#[repr(C)] #[derive(Clone, Copy, Debug)]
pub struct Quat<S> { pub w: S, pub x: S, pub y: S, pub z: S }
#[repr(transparent)] #[derive(Clone, Copy, Debug)]
pub struct SO3<S>(Quat<S>);
```

- Constructors: `Quat::from_wxyz_unchecked` (debug-asserts `NUMERICS.md` §3.6's bound),
  `from_wxyz_normalized`, `from_xyzw`, `to_xyzw`, `from_jpl` (Hamilton $(w, -x, -y, -z)$ from a JPL
  $(x, y, z, w)$ — same matrix; Sommer et al. 2018; frame semantics stay the caller's).
  `SO3::from_quat_unchecked`, `from_quat_normalized`, `quat()`.
- `exp`, `log` (`NUMERICS.md` §3.1–§3.2, with the `copysign` flip), `act`, `act_many`,
  `to_matrix`, `from_matrix` (§3.4; nested `branch`, closed form), `renormalize`, Hamilton product
  as `Mul`, conjugate as `inverse`, `adjoint = R` (as `Mat3`), `ad = W`, `jr`, `jr_inv`, `jl`,
  `jl_inv` from `coeffs`. `SO3::Jac = Mat3<S>`.
- **`from_matrix_never_iterates`**: locus-tag's degenerate near-singular IPPE matrices (the inputs
  of its `quat_from_so3` regression) are fixtures; the test asserts completion and backward error.

## 5. SE_N(3)

**NORMATIVE.**

```rust
#[repr(C)] #[derive(Clone, Copy, Debug)]
pub struct SEn3<S, const N: usize> { q: Quat<S>, x: [Vec3<S>; N] }
pub type SE3<S> = SEn3<S, 1>;
pub type SE23<S> = SEn3<S, 2>;

#[derive(Clone, Copy, Debug)]
pub struct SEn3Tangent<S, const N: usize> { pub phi: Vec3<S>, pub rho: [Vec3<S>; N] }
pub type Twist<S> = SEn3Tangent<S, 1>;       // omega() = phi, v() = rho[0]

/// Rotation-first [[A, 0], [B_i, A]] (NUMERICS.md §2.2).
#[derive(Clone, Copy, Debug)]
pub struct SEn3Jac<S, const N: usize> { pub diag: Mat3<S>, pub col: [Mat3<S>; N] }
```

- `SE3` accessors: `rotation()`, `translation()`, `from_rt`, `from_quat_translation`; `SE23`:
  `rotation()`, `velocity()`, `position()` with the order $x_1 = v$, $x_2 = p$ stated in rustdoc
  (InEKF convention).
- `exp`, `log`, `adjoint`, `ad`, `jr`, `jl`, `jr_inv`, `jl_inv` per `NUMERICS.md` §5; `jr_inv` is
  the dual-matrix inverse of `jr` (§5.4), never a separate closed form.
- `SEn3Jac`: `mul`, `inverse`, `apply`, `apply_transpose`, `write_dense` per §2.2; cost 27 + 54N
  multiplications per product.
- Action for `SE3` only (`impl<S> Mul<Point3<S>> for SE3<S>`), `act_many`, `act_jacobians::<Sd>` →
  `(Matrix<S, 3, 6>, Mat3<S>)` per `NUMERICS.md` §2.4.
- Twist converters: `Twist::from_translation_first([S; 6])`, `to_translation_first` (locus-tag,
  Sophus, manif orders), per API R3.

## 6. SO(2) and SE(2)

`SO2<S> { c, s }`, `SE2<S> { rot: SO2<S>, t: Vec2<S> }`, tangent `SE2Tangent { theta, rho: Vec2 }`,
`SE2::Jac = Mat3<S>` (dense; three DOF do not pay for structure). Formulas `NUMERICS.md` §6; the
PR adds the rotation-first permuted SE(2) Jacobians to `NUMERICS.md` §6 before code.

## 7. Rⁿ and products

`Rn<S, N>(Vector<S, N>)`: addition as the group law, tangent `RnTangent { rho: Vector<S, N> }`,
`Jac = RnJac`, the matrices `k I` (`jr = I`, `Ad = I`, `ad = 0`, `-I` in the `⊖` rows; a `Mat<N>`
has no inverse at general `N`), so generic code needs no special case. `Product<A, B>(A, B)` with
tangent `(A::Tangent, B::Tangent)`, `ProductJac(A::Jac, B::Jac)` block-diagonal,
`DOF = A::DOF + B::DOF`, dense order A then B. `Product<SO3<S>, Rn<S, 3>>` is the tf2-semantics
pose.

## 8. Side and action Jacobians

Implement every row of `NUMERICS.md` §2.3 for both sides, and §2.4. Every one is checked two ways:
against `Dual<S, DOF>` differentiation of the implemented operation (`jacobians_match_dual_*`), and
against the chained-primitive reference twin (`*_matches_reference`).

## 9. Tests

**NORMATIVE.** Named proptests, each per group and per precision (`f64`, `f32`) unless noted; bounds
are recorded measurements, not asserted constants.

- `group_axioms_*`: associativity, identity, inverse.
- `exp_log_roundtrip_*` ($\theta < \pi$): backward error.
- `adjoint_identity_*`: $X\,\mathrm{Exp}(\tau)\,X^{-1} = \mathrm{Exp}(\mathrm{Ad}_X\tau)$.
- `jl_is_ad_jr_*`: $J_l = \mathrm{Ad}_{\mathrm{Exp}\tau} J_r$.
- `jacobians_match_dual_*`: every Jacobian vs `Dual` through the shipped operation.
- `dual_value_is_plain_value_*`: every group method under `Dual<f64, DOF>`.
- `branch_continuity_*`: at each generated switch point (`NUMERICS.md` §4).
- `so3_log_w_flip`, `so3_log_at_w0_is_a_function_of_the_sign`.
- `from_matrix_never_iterates` (§4).
- `tf_tree_math_differential` (`f64`, runner-backed in CI): helicoid vs `tf_tree_math` on 10⁵
  random pairs, near-identity and near-π included; the recorded max difference is quoted in the
  Phase 4 migration evidence.
- Corpus: every public numeric function id has its corpus file; `just envelope` green.

## 10. The first envelope

Run `just oracles` and `just envelope`; **every paired stratum must be dominated**. A stratum where
an oracle wins blocks the phase: fix the implementation or write a record explaining why the oracle
is measuring something else. Bless the baseline; commit `docs/evidence/ENVELOPE.md`.

## 11. Benches

`exp`, `log`, `compose`, `act`, `act_many`, `jr`, `jr_inv`, `adjoint`, `SEn3Jac::mul` for SO(3),
SE(3), SE₂(3); fixtures per stratum (near-identity, generic, near-π); `f64` and `f32`.
`tf_tree_math` and sophus-rs rows reported beside them. Blessed into `baseline/`.

## 12. Definition of done

- [ ] §2–§8 implemented; `NUMERICS.md` §6 carries the permuted SE(2) Jacobians.
- [ ] `generated.rs` committed from a sweep; `thresholds-check` green.
- [ ] §9 tests green on x86_64, aarch64, wasm32; `no-std` build green.
- [ ] Envelope: every paired stratum dominated; baseline and `ENVELOPE.md` committed.
- [ ] `cargo xtask lint` twin-table check green for every Phase 3 row of `NUMERICS.md` §14.
- [ ] Benches blessed.

## Appendix: suggested implementation order

1. Traits, `Side`, `Rn` (the trivial group proves the traits compile generically).
2. Coefficient kernel, `__sweep`, run the sweep, commit `generated.rs`.
3. `Quat`, `SO3` (port from `tf_tree_math`), its corpus ids, envelope rows.
4. `SEn3` with `N = 1`, then generic `N`; `SEn3Jac` algebra.
5. `SO2`, `SE2`.
6. `Product`.
7. Side and action Jacobians; `Dual` checks.
8. Reference twins; twin-table lint.
9. Envelope, bless; benches.
