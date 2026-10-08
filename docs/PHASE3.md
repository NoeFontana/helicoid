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
| Traits: `Tangent`, `Jac`, `LieGroup`, `Side` (§2) | Partial: core traits implemented and verified under `laws.rs` for `Rn`, `SO3`, `SEn3` (`0025`). Chained primitive twins pending. |
| Coefficient kernel + `__sweep` + generated thresholds (§3) | Partial: `helicoid::coeffs` implemented with two-stage sweep and generated switches (`0039`, `0047`). `se2_coeffs` pending. |
| `Quat`, `SO3` (§4) | Partial: `Quat` and `SO3` implemented, dominating oracles across all `so3_*` corpus ids. `renormalize` scored as `quat_renormalize`, its `renorm:*` strata inside `NUMERICS.md` §3.6's band at both precisions ([`0056`](./decisions/0056-the-routines-d7-does-not-reach.md)). |
| `SEn3<S, N>`, `SEn3Tangent`, `SEn3Jac`; `SE3`, `SE23` (§5) | Partial: `SEn3` implemented for general $N$ (`0042`); structural zeros skipped in products; verified across 21 `sen3_*` corpus ids. |
| `SO2`, `SE2` (§6) | Not started |
| `Rn`, `Product` (§7) | Partial: `Rn` and `Product<A, B>` implemented with structured `ProductJac` and constructors (`0025`, `0029`). |
| Side Jacobians, action Jacobians (§8) | Partial: `Dual<S, DOF>` automatic differentiation verified for all group Jacobian rows; action Jacobians implemented. |
| Reference twins and proptests (§9) | Partial: dense reference twins and proptests for `SEn3Jac`, `ProductJac`, and `SEn3::jr_inv` implemented. |
| Envelope blessed; `docs/evidence/ENVELOPE.md` (§10) | Not started: 1559 candidate rows evaluated; baseline blessing pending final Phase 3 criteria. |
| Benches (§11) | Partial: criterion group benchmarks in `benches/groups.rs` with replicate-outermost gating (`0033` draft). |

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
    /// The products of the dense components accumulated into `acc`, in index order.
    fn dot_acc(&self, o: &Self, acc: S) -> S;
    /// Provided. The `+0` seed is normative: `+0 + -0` is `+0`.
    fn dot(&self, o: &Self) -> S { self.dot_acc(o, S::zero()) }
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

// Sealed: `Right` and `Left` are the only implementations (0025). The SO(3) PR took the
// first-class selector over a `TypeId` comparison: `IS_RIGHT` is an associated const, so
// `match Sd::IS_RIGHT` resolves at monomorphization and neither arm survives into the emitted
// code. Sealing is what kept that choice internal.
pub trait Side: sealed::Sealed + Copy + 'static {
    const IS_RIGHT: bool;
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
compute array lengths from associated consts; a `const` assertion ties them in every impl. The same
limitation is why `dot_acc`, not `dot`, is the required operation: a composite tangent cannot size
a flattening buffer from `DOF`, so it threads the accumulator instead
([`0025`](./decisions/0025-a-structured-jacobian-and-a-sealed-side.md)).

## 3. The coefficient kernel

**NORMATIVE** ([`0004`](./decisions/0004-switch-points-are-generated-not-typed.md)).

- `coeffs` is `pub(crate)` and **the only module that evaluates `NUMERICS.md` §4**. Call-site
  groups: `exp_coeffs(θ²) -> (k, cos_half)`, `jr_coeffs(θ²) -> (a, b)`,
  `jr_inv_coeff(θ²) -> c`, `q_coeffs(θ²) -> (b, d, e)`, `log_ratio(n², w) -> r`,
  `se2_coeffs(θ²) -> (α, β)` (`gamma2_coeffs` arrives in Phase 5). Each is **one** `S::branch` on
  "every member is on its series arm"; the exact arm uses the safe argument and runs once for the
  group, and each member then selects its own arm by its own mask — so **grouping is a cost and
  never a value**, whichever members share a branch
  ([`0047`](./decisions/0047-the-second-arm-is-admitted-by-agreement-not-by-the-objective.md)).
- `generated.rs` holds a `Switch` per coefficient per precision; kernels select by
  `S::PRECISION` (a `const` match, resolved at monomorphization). A `Switch` carries **two** series
  arms, the second the first $m_0$ terms of the same series below a second switch (Phase 1 §6): the
  sweep admits it only where it agrees with the whole arm to the bit, so it is a latency degree of
  freedom and the member's value is unchanged.
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
  `from_wxyz_normalized` (divides by the norm, [`0027`](./decisions/0027-a-normalizing-constructor-normalizes.md);
  `renormalize` keeps the Newton step), `from_xyzw`, `to_xyzw`, `from_jpl` (Hamilton $(w, -x, -y, -z)$ from a JPL
  $(x, y, z, w)$ — same matrix; Sommer et al. 2018; frame semantics stay the caller's).
  `SO3::from_quat_unchecked`, `from_quat_normalized` (divides by the norm, `0027` decision 4),
  `quat()`.
- `exp`, `log` (`NUMERICS.md` §3.1–§3.2, with the `copysign` flip), `act`, `act_many`,
  `to_matrix`, `from_matrix` (§3.4; nested `branch`, closed form), `renormalize`, Hamilton product
  as `Mul`, conjugate as `inverse`, `adjoint = R` (as `Mat3`), `ad = W`, `jr`, `jr_inv`, `jl`,
  `jl_inv` from `coeffs`. `SO3::Jac = Mat3<S>`.
- **`from_matrix_never_iterates`**: degenerate near-singular IPPE matrices
  are fixtures; the test asserts completion and backward error.

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

/// Rotation-first [[A, 0], [B_i, A]] (NUMERICS.md §2.2). The blocks are a representation of the
/// dual matrix and not values of it, so they are crate-private and reached from outside through
/// `write_dense` and `apply` (0025 decision 5; 0028 option A corrected this declaration, which
/// spelled them `pub`).
#[derive(Clone, Copy, Debug)]
pub struct SEn3Jac<S, const N: usize> { pub(crate) diag: Mat3<S>, pub(crate) col: [Mat3<S>; N] }
```

- Accessors: `rotation()` and, at every `N`, the exact pair `from_parts(SO3<S>, [Vec3<S>; N])` /
  `parts()` ([`0042`](./decisions/0042-a-general-sen3-needs-a-way-to-be-built-from-its-parts.md));
  `SE3` adds `translation()`, `from_rt`, `from_quat_translation`; `SE23` adds `velocity()`,
  `position()` with the order $x_1 = v$, $x_2 = p$ stated in rustdoc (InEKF convention). The named
  forms are the ones to prefer where they exist, and the generic pair is what a width without names
  has — the corpus's `sen3_*_n3` records are an element's parts and nothing else can read them.
- `exp`, `log`, `adjoint`, `ad`, `jr`, `jl`, `jr_inv`, `jl_inv` per `NUMERICS.md` §5; `jr_inv` is
  the dual-matrix inverse of `jr` (§5.4), never a separate closed form.
- `SEn3Jac`: `mul`, `inverse`, `apply`, `apply_transpose`, `write_dense` per §2.2; cost 27 + 54N
  multiplications per product.
- Action for `SE3` only (`impl<S> Mul<Point3<S>> for SE3<S>`), `act_many`, `act_jacobians::<Sd>` →
  `(Matrix<S, 3, 6>, Mat3<S>)` per `NUMERICS.md` §2.4.
- Twist converters: `Twist::from_translation_first([S; 6])`, `to_translation_first` (translation-first,
  Sophus, manif orders), per API R3.

## 6. SO(2) and SE(2)

`SO2<S> { c, s }`, `SE2<S> { rot: SO2<S>, t: Vec2<S> }`, tangent `SE2Tangent { theta, rho: Vec2 }`,
`SE2::Jac = Mat3<S>` (dense; three DOF do not pay for structure). Formulas `NUMERICS.md` §6; the
PR adds the rotation-first permuted SE(2) Jacobians to `NUMERICS.md` §6 before code.

## 7. Rⁿ and products

`Rn<S, N>(pub Vector<S, N>)`: addition as the group law, tangent `RnTangent { rho: Vector<S, N> }`,
`Jac = RnJac`, the matrices `k I` (`jr = I`, `Ad = I`, `ad = 0`, `-I` in the `⊖` rows; a `Mat<N>`
has no inverse at general `N`), so generic code needs no special case. `Product<A, B>(A, B)` with
tangent `(A::Tangent, B::Tangent)`, `ProductJac(A::Jac, B::Jac)` block-diagonal,
`DOF = A::DOF + B::DOF`, dense order A then B, and `dot_acc` threaded as A then B so that the sum
is the flat one §2 specifies. `Product<SO3<S>, Rn<S, 3>>` is the tf2-semantics pose. Rⁿ's field is
`pub` because the group is abelian, so a reach-through to `Vector`'s `Add` can name no side; a
non-abelian group exposes no such field
([`0025`](./decisions/0025-a-structured-jacobian-and-a-sealed-side.md)).

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
is measuring something else — and that explanation is a row in
`conformance/baseline/exceptions.toml` naming the `ready` record, which excepts **domination only**
and fails the run once the stratum is dominated after all
([`0046`](./decisions/0046-explained-by-record-needs-a-record-to-point-at.md)). Bless the baseline;
commit `docs/evidence/ENVELOPE.md`.

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
