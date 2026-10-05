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
| Traits: `Tangent`, `Jac`, `LieGroup`, `Side` (§2) | Partial ([`0025`](./decisions/0025-a-structured-jacobian-and-a-sealed-side.md)): the four traits as written in §2, `Side` sealed to `Left`/`Right`, `dot_acc` the required operation; the `DOF` ties and `sandwich`'s `D == T::DOF` are `const` assertions, so a mismatch fails `cargo build`/`cargo test`, not `cargo check` (`sandwich`'s pinned by a `compile_fail` doctest); implemented by `Rn` and by a test-only non-abelian group (Heisenberg, `heis_tests.rs`), both under every generic law of `laws.rs` at `f64`, `f32` and `Dual<f64, 3>`. The side selector of `compose_jacobians`/`inverse_jacobian` is `Side::IS_RIGHT`, an associated const resolved at monomorphization (the SO(3) PR's call, taken over a `TypeId` comparison; `SO3` is the first group whose two sides differ, and `laws::side_delegation` holds at 1.118 `u` for it). Not started: comparing the §2.3 rows with `Dual` differentiation (`jacobians_match_dual_*` and the §8 twins); `Tangent` and `Jac` are also implemented by `SEn3Tangent` and `SEn3Jac` (§5), for which the generic laws of `laws.rs` are not instantiated (their own tests cover them) |
| Coefficient kernel + `__sweep` + generated thresholds (§3) | Partial: `helicoid::coeffs` (`pub(crate)`; SO(3) is now its first consumer, so a default build reaches `exp_coeffs`, `jr_coeffs`, `jr_inv_coeff` and `log_ratio` from `SO3`, and the `cfg` the module carried is gone. `q_coeffs` is SE_N(3)'s alone and `#[allow(dead_code)]` until §5, rather than `cfg`'d, which would split its exact arms across two builds; the generated catalogue is allowed whole, since a hand edit fails `just thresholds-check`. Its `nonnegative` assert now reads "no lane is negative" instead of "every lane is >= 0": the two differ only on NaN, which the assert is not there to reject, and every arm returns NaN for a NaN branch variable — the answer a value function owes, and what `laws::dual_value_is_plain_value` asks of `SO3::exp` over `f64::ANY`): `exp_coeffs`, `jr_coeffs`, `jr_inv_coeff`, `q_coeffs` and `log_ratio`, each one `S::branch` on "every member is on its series arm", off which the group's exact arms run once, at the safe argument, and each member selects by its own mask; the exact arms in the operand order of the seeded kernels (`docs/maths/coefficients.md` CO.6), the series by Horner in the branch variable, the tables chosen by `S::PRECISION`; `generated.rs` (a `Switch` per coefficient per precision, taking its terms from the eight swept series that follow it) and `conformance/sweeps/thresholds.csv`, written by `cargo xtask thresholds helicoid` from a sweep of the shipped `exact_*`/`series_*` through `__sweep`, a fixed point (a second run writes the same bytes; `thresholds-check` passes); the shipped arms are the seeded arms bit for bit on every record where the series arm is defined (`w > 0` for `r`), and the shipped groups score the recorded objective over the `theta:*` strata to the bit, and are finite on every stratum, at `f64` and `f32`; tests in `coeffs/tests.rs`: `branch_continuity_*` at every switch (the jump at most the sum of the two arms' errors at the two records that bracket it, a sample and not a bound over the interval; not §4's wording: 0015 (draft) NU.4), a group is its members' own arms, no exact arm below a group's smallest switch and each `sin_cos` once above it (a counting scalar), the exact arms' bits at two `z`, `dual_value_is_plain_value_*`, no non-finite operation in a lane that evaluates both arms (`θ² = 0`, subnormal, each switch, `r` down to `w = 1e-100`), the series at `0`. Readings, provisional: a group's members each have their own switch (0015 (draft) NU.5) so no exact arm runs below the group's smallest switch (`Exp`: `θ² < 5.6e-15` at `f64`, `2.6e-6` at `f32`, where `cos θ/2` switches, so from `θ = 7.5e-8`, `1.6e-3` it is on its exact arm, at one `sqrt` and one `sin_cos` as §3.1 states; `jr_coeffs`: `0.72` at `f64`) and above it each member selects by its own mask (a scalar mask runs the Horner of the series members only; one `sin_cos` of `jr_coeffs` for `θ² ∈ [0.72, 0.96)` and of `q_coeffs` for `[0.96, 1)` is evaluated for nothing at `f64`, a division of `exp_coeffs`; unmeasured, no bench yet); `r` takes its series iff `w > 0` and `n²/w² < switch`, the quotient formed by a division (0014 (draft) question 29, 0015 (draft) NU.6), so a lane that evaluates both arms sees a non-finite operation where `w > 0` and the quotient overflows (`w` below about `1e-154` at `f64`, `5e-20` at `f32`; the value is right); each recorded objective is a maximum over the corpus's records, not over the domain, and 8 of the 16 switches (`k`, `d` at `f64`; `k`, `a`…`e` at `f32`) are `θ² = 1`, the top of the sweep's grid, not a measured optimum (0014 (draft) question 9). Missing: `se2_coeffs` (0015 (draft) PH.3: `β` is odd in `θ`), `gamma2_coeffs` (Phase 5), the per-call-site bench; `just no-std` and `just wasm` build the kernel with `--features __sweep` |
| `Quat`, `SO3` (§4) | Partial: `Quat` in `quat.rs`: `Blend`, Hamilton product as `Mul`, `conjugate`, `norm_sq`, `to_matrix` (`R(q)` of `NUMERICS.md` §1, a scaled rotation for a non-unit `q`), `from_wxyz_unchecked` (`debug_assert!` of the §3.6 bound, by `S::PRECISION`), `from_wxyz_normalized` (divides by the norm, `0027`; domain `\|q\|²` normal, `debug_assert!` on the result) and `renormalize` (the first-order Newton step, asserting nothing as §3.6 states no domain for it; a normalization only within `\|η\| <= 2^-26.29` (`f64`) or `2^-11.79` (`f32`), which contains the `from_wxyz_unchecked` domain by 13.7 bits and 4.2 bits), `from_xyzw`, `to_xyzw`, `from_jpl` (array arguments, unchecked permutations); checked under `f64`, `f32`, `Dual<f64, 4>` (associativity, `q q* = 1`, `R` orthogonal with `det = 1` and equal to `q v q*`, `from_jpl` against the JPL matrix, the `from_wxyz_unchecked` domain on both sides, quadratic convergence of the step, `from_wxyz_normalized` unit from a non-unit `q`, the `Dual` value path bit-identical to `f64`); `f64` and `f32` only: the Hamilton identities, the bit-exact left-to-right goldens of the product and `norm_sq`, the hand cases; the `Dual` derivative lanes are checked for `Mul` and `norm_sq` only. **`SO3` done**: `from_quat_unchecked`, `from_quat_normalized`, `quat`, `exp`, `log`, `act`, `act_many`, `to_matrix`, `from_matrix`, `renormalize`, `Mul`, `inverse`, `adjoint = R`, `ad = W`, `jr`, `jr_inv` and every §2.3 row, with `Jac = Mat3<S>` implementing `Jac` (its `inverse` the algebraic one through `Mat3::inverse_adj`, whose quotient is already divided by `det`); `jl`/`jl_inv` are **not** overridden, the provided `jr(−τ)` being `I + aW + bW²` to the bit, so one code path and §14's twin stays the provided body. `Exp` builds its quaternion field by field rather than through `from_wxyz_unchecked`, whose assert is for a caller claiming unit and would turn a NaN tangent into a panic. Under every generic law of `laws.rs` at `f64`, `f32` and `Dual<f64, 3>`, `f64` and `Dual` agreeing to every digit, with `jacobian_rows` and `tangent_dense_order` at exactly **0** — so each §2.3 row is reproduced bit for bit — the worst law over 60 000 samples being `ad_consistency` at 7.693 `u` (`f64`, `Dual`) and `plus_minus` at 7.478 (`f32`); every bound is twice the larger measured maximum, as `Bounds` asks, and the module doc carries the table. The sample count is part of a bound here: a first pass set them from 600 samples and was wrong in both directions — `f32` `group_axioms` at 6.0 against a measured 4.29 that `proptest` then beat with 6.32, and `f64` `ad_consistency` and `plus_minus` left at 8.0 against maxima of 7.69 and 7.43. `Side` gained `IS_RIGHT`, the §2 selector left to this PR. Scored by the `helicoid` subject over **every `so3_*` corpus id**, `f64`, no non-finite output: `so3_exp` 2.743 (`theta:pi-1e-9`, equal to the seeded subject's), `so3_log` 2.627 (`theta:1e-4`; `theta:dense` 2.465 against the seeded reading's 2.901, the swept `r` switch on `s = n²/w²` beating the seeded `n² = 0` one), `so3_jr` and `so3_jl` 4.097 (`theta:pi-1e-8`), `so3_jr_inv` and `so3_jl_inv` 2.112 (`theta:dense`), `so3_act` 4.899 (`theta:pi-1e-5`), `so3_from_matrix` shape and finiteness only (its `max_u` is NaN by design, §4.3). `so3_act`'s `q:nonunit` row is 510.3, which is the **open reading** of §3.3 off the unit sphere and not a defect: it and the scaled form's 257.8 are `2·2^-45/u` and `2^-45/u`, so the corpus reference is the normalized one (maths index). `0029`'s `from_parts`/`parts` land here, on `Product`, since SO(3) is what makes `Product<SO3, Rn<3>>` constructible (§7's row). `act_many`'s §14 twin now has its `*_matches_reference` proptest (D6), and it found that **§14's `3 u` for that row is too tight**: measured 7.587 `u` near `π`, 4.873 at `θ ~ 1`, 3.379 at `θ ~ 1e-6`. Not an implementation choice — §3.3 prescribes forming `R(q)` once, so its nine entries round before any point is touched while `act` rounds a sandwich per point. Raising it is a `NUMERICS.md` §14 edit and a record, owed. Not started: a `renormalize` corpus id and any `Quat` stratum of its own (the `q:w0` and `q:nonunit` records `so3_log`, `so3_act` and `so3_from_matrix` carry are scored, 64 each), and `jacobians_match_dual_*` |
| `SEn3<S, N>`, `SEn3Tangent`, `SEn3Jac`; `SE3`, `SE23` (§5) | Partial: `SEn3Tangent` (`sen3.rs`: `Tangent` with `DOF = 3 + 3N`, `dot` one flat sum in index order) with `Twist`, `omega()`, `v()` and `from_translation_first`/`to_translation_first` (the permutation of `docs/maths/se3.md` SE.14(a), exact); `SEn3Jac` (`dualmat.rs`) as a `Jac`: `mul` (27 + 54N multiplications, measured with a counting scalar), `inverse` (`A⁻¹` from `Mat3::inverse_adj`; the computed `det A` nonzero and finite is `debug_assert!`ed, and a release build returns non-finite entries at `det A = 0` where `NUMERICS.md` §12 promises `jr_inv` a finite value; that and `NUMERICS.md` §5.4 and `0005` taking `A⁻¹` from the SO(3) closed form, which `inverse(&self)` cannot receive: for the `jr_inv` PR to settle), `neg`, `apply`, `apply_transpose`, `write_dense`, `identity`; `sandwich` forms `M = J Σ` by columns and `M Jᵀ` by rows (row `r` of `M Jᵀ` is `J` applied to row `r` of `M`) through `from_cols`/`from_rows`, so it holds no `D x D` scratch and reaches no entry through `Matrix::get`/`set`, and forms neither `Aᵀ` nor the `N` transposes of `col`: bit-identical to the block form it replaces, checked over 150 000 random `(J, Σ)` at every `N` and every scalar with not one differing entry, so the recorded bounds are unmoved rather than re-measured, and measured at 545 → 48 ns (`f64`, `D = 6`) and 8078 → 2328 ns (`Dual<f64, 3>`, `D = 12`); the blocks are `pub(crate)`, as `0025` decision 5 requires and `ProductJac` already read it: `0028` option A, which corrected this section's declaration from the `pub diag`/`pub col` it spelled; checked under `f64`, `f32`, `Dual<f64, 3>` for `N = 1, 2, 3` against the twins (§9 row below), `write_dense` against the matrices of `NUMERICS.md` §2.2 (written out for `N = 1, 2`, entries for `N = 3`), `sandwich`, `apply` and `apply_transpose` against their §14 twins, `neg` and `Blend` over every `Dual` lane, exact integer algebra, the `Dual` value path bit-identical to `f64`, `twist_translation_first_round_trip`; not started: `SEn3`, `SE3`, `SE23` and everything of §5 that needs the group (`exp`, `log`, `adjoint`, `ad`, `jr`, `jl`, `jr_inv`, `jl_inv`, the action) and the corpus ids |
| `SO2`, `SE2` (§6) | Not started |
| `Rn`, `Product` (§7) | Partial ([`0025`](./decisions/0025-a-structured-jacobian-and-a-sealed-side.md)): `Rn<S, N>` with tangent `RnTangent { rho }` and `RnJac`, a scalar multiple `k I` (`I`, `-I` and `ad = 0`), correcting §7's `Jac = Mat<N>`, which cannot implement `Jac::inverse`; every §2.3 row for both sides, all `±I` and symmetric; addition is `Mul`, with `Rn + Rn` and `RnTangent + RnTangent` pinned by `compile_fail` doctests; `Product<A, B>` (`product.rs`): `LieGroup` with every method, the provided ones included, the factors' method componentwise (so a factor's override of a provided method is what the product uses; no factor overrides one yet, so that is untested until SO(3)), tangent the pair `(A::Tangent, B::Tangent)` (`Tangent` for a pair, dense order `A` then `B`, `dot_acc` threaded `A` then `B`, so `dot` is the one flat index-order sum §2 specifies at any nesting depth and `tangent_dense_order` holds at a recorded `f64` bound of `0`), `DOF = A::DOF + B::DOF` (the `identity` tie plus the factors' own), every §2.3 row per side as the factors' rows side by side, `Mul` and `Blend` componentwise; `ProductJac<JA, JB>` (its parameters are the factors' Jacobian types, `Product<A, B>::Jac = ProductJac<A::Jac, B::Jac>`) as a `Jac`: `mul`, `inverse`, `neg`, `apply`, `apply_transpose` per block, `write_dense` into the two diagonal blocks (the off-diagonal blocks `+0`, as two rectangular loops, not one `DOF x DOF` sweep with a test per entry), `sandwich` from the factors' `apply` alone (`M = J Σ` by columns, then `M Jᵀ` by rows, since row `r` of `M Jᵀ` is `J` applied to row `r` of `M`; the factors' own `sandwich` needs a block of size `A::DOF`, which a generic `D` cannot build); the fields of neither type are `pub`, as §7 spells them and `0025` decision 5 lists, so outside the crate a `Product` is built and read only through `exp`/`log` and `write_dense`/`read_dense` (R3) — §7 names `Product<SO3, Rn<3>>` the tf2 pose but gives no way to build one from a rotation and a translation without an `Exp`/`Log` round trip, which for SO(3) is neither exact nor cheap and is worst-conditioned at a half turn: `0029` settles it: option A, `from_parts`/`parts` land with the SO(3) PR; checked by the generic laws of `laws.rs` under `f64`, `f32`, `Dual<f64, D>` on `Product<Rn<3>, Rn<2>>`, `Product<Rn<2>, Heis>`, `Product<Heis, Rn<2>>` and the nested `Product<Product<Heis, Rn<2>>, Heis>` (`D = 8`), plus hand cases (block-diagonal dense layout on a padded view, `Jac::mul` operand order against the dense product on non-commuting blocks, every row against the factors' rows, `sandwich` cross blocks, `Blend`, `dot` one flat sum, the out-of-domain behaviour of both profiles) and the §9 twin row below; the only non-abelian factor is the test-only Heisenberg group (2-step nilpotent), so the block structure is exercised on a non-abelian factor with no series term past the first, and on a rotation factor (`Product<SO3, Rn<3>>`) only when SO(3) lands; not started: `jacobians_match_dual_*` and the `Gaussian`/action uses of a product |
| Side Jacobians, action Jacobians (§8) | Not started |
| Reference twins and proptests (§9) | Partial: every `SEn3Jac` row of `NUMERICS.md` §14 — `reference::sen3jac_mul`, `sen3jac_inverse` (dense product; Gauss–Jordan with row exchange by `select`), `sen3jac_apply`, `sen3jac_apply_transpose` (dense `J x` and `Jᵀ x`) and `sen3jac_sandwich` (dense `J Σ Jᵀ`, which the Phase 5 `Gaussian` row reuses), all writing into caller memory, with `sen3jac_*_matches_reference` per scalar and `N = 1, 2, 3` (the inverse against `κ u` of the dense matrix, on entries uniform on `[-1, 1)` with `κ u <= 1e-3`: ill-conditioned inputs are unchecked until the `sen3_jr_inv_n*` corpus ids; `apply` and `sandwich` on the same entries, which are not exactly representable, so the association the structured form takes is what the bound covers); `reference::productjac_sandwich` (the dense `J Σ Jᵀ` of a block-diagonal `J`, sharing the `D⁴` body of `sen3jac_sandwich`), with `productjac_sandwich_matches_reference` over the four product groups per scalar, on the value part (the `Dual` lanes of `sandwich` wait for `jacobians_match_dual_*`, §8); every other row of §14 and the twin-table lint not started |
| Envelope blessed; `docs/evidence/ENVELOPE.md` (§10) | Not started, and no longer vacuous: `cargo xtask envelope` judges `helicoid`'s 455 rows — **170 strata paired with an oracle**, 255 with none — and fails on no baseline and **23 domination failures**, every one of them attributed (5 the program's, 2 D16's, 16 neither). The 16 are `so3_jr` and `so3_jl` on the eight `theta:pi-1e*` strata, 1.3–1.7x behind sophus-rs on 64 records each, and they are the assembly and not the backend by measurement: the `libm` crate and glibc do not disagree on one `sin`/`cos` argument those strata reach, and `seeded:host-std`'s rows there are `helicoid`'s to the bit (`host::the_swap_has_no_power_near_pi_on_the_so3_jacobians`; 0037 (draft)). The one loss to `tf_tree_math` is `so3_exp`/`theta:1e-5`, 0.4952463054990454 u against 0.4952463054734329 u at the same argmax record — the two agree to ten digits, both under ½ `u`, and `PHASE4.md` §5.2 reads it as a stratum the oracle wins (0038 (draft)). `--bless` writes nothing while a bar fails, so there is no baseline, the 425 scored rows carry no no-regress bar, and the page still reads that nothing is compared |
| Benches (§11) | Not started; owed first, `exp` at near-identity `θ` from `7.5e-8` to `1`, where the kernel runs one `sqrt` and one `sin_cos` (the kernel row above) |

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
  `from_wxyz_normalized` (divides by the norm, [`0027`](./decisions/0027-a-normalizing-constructor-normalizes.md);
  `renormalize` keeps the Newton step), `from_xyzw`, `to_xyzw`, `from_jpl` (Hamilton $(w, -x, -y, -z)$ from a JPL
  $(x, y, z, w)$ — same matrix; Sommer et al. 2018; frame semantics stay the caller's).
  `SO3::from_quat_unchecked`, `from_quat_normalized` (divides by the norm, `0027` decision 4),
  `quat()`.
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

/// Rotation-first [[A, 0], [B_i, A]] (NUMERICS.md §2.2). The blocks are a representation of the
/// dual matrix and not values of it, so they are crate-private and reached from outside through
/// `write_dense` and `apply` (0025 decision 5; 0028 option A corrected this declaration, which
/// spelled them `pub`).
#[derive(Clone, Copy, Debug)]
pub struct SEn3Jac<S, const N: usize> { pub(crate) diag: Mat3<S>, pub(crate) col: [Mat3<S>; N] }
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
