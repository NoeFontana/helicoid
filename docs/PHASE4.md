# helicoid — Phase 4 Implementation Specification: Geodesics and the `tf_tree` Migration

> **Companions:** `docs/NUMERICS.md` §10, [`0010`](./decisions/0010-seeded-from-tf-tree-math.md),
> and — in the `tf_tree` repository — `docs/PROJECT.md` D5, D12, D13, D14 and `docs/PHASE1.md` §3.

**Deliverable:** geodesics on every group, their Jacobians and velocity, the dual-quaternion fast
twin for SE(3), the invariance tests; and everything a `tf_tree` decision record needs to move
`tf_tree_math` onto `helicoid`: the parity table, the envelope, the bench evidence. Sections marked
**NORMATIVE** are requirements.

## 0.0 Implementation status

| Area | Status |
|---|---|
| `LieGroup::geodesic`, reference twin (§1) | **Done** (§1.1, §1.3). `SO3::geodesic` overridden with GE.14's two-arm switch ([`0050`](./decisions/0050-the-geodesic-s-denominator-is-the-whole-domination-gap.md), [`0051`](./decisions/0051-two-arms-on-the-switch-the-sweep-already-chose.md)). `Product` delegates per factor. Verified via `laws::geodesic`. |
| SE(3) dual-quaternion fast twin (§1.2) | **Done** ([`0054`](./decisions/0054-the-screw-twin-takes-two-arms-in-the-world-frame.md)). `SE3::geodesic` overridden: `SO3::geodesic`'s rotation; the translation in the world frame, the definition's closed form below `r`'s second switch and GE.12's power above it (GE.15). `se3_geodesic` 1.572 / 2.057 / 2.721 `u` against `ScLerp`'s 2.336 / 2.502 / 3.253, all three dominated (`envelope` 113 → 112); 0.59–0.92× the provided body by `bench-gate`; 4.6–6.1% faster than `screw_pow` inside `tf_tree`'s `lookup` (0054's A/B). `N ≥ 2` keeps the provided body. |
| `geodesic_jacobians`, `geodesic_velocity` (§2) | Partial: `LieGroup::geodesic_velocity` ships (`x₁ ⊖_R x₀`); `laws::geodesic_legs` passes. Owed: `geodesic_jacobians` and dual twin ([`0043`](./decisions/0043-the-geodesic-jacobian-ships-the-cancellation-free-form.md)). |
| Invariance tests (§3) | Partial: Left/right invariance tested on all groups via `laws::geodesic`. `Product<SO3, Rn<3>>` asserts exact right-invariance under its own law and expected positive failure under SE(3) ([`0045`](./decisions/0045-two-phase-4-checks-cannot-be-taken-as-written.md)). |
| Corpus id `se3_geodesic`, `so3_geodesic`; strata (§4) | **Done.** 180 records each in `conformance/corpus/{so3,se3}_geodesic.jsonl` across `geo:consecutive`, `geo:generic`, and `geo:near-pi` ([`0045`](./decisions/0045-two-phase-4-checks-cannot-be-taken-as-written.md)). |
| Parity table and evidence for `tf_tree` (§5) | Not started |

## 0. Non-goals and guardrails — read first

**NORMATIVE.** `helicoid` does not change `tf_tree`. The migration is a **`tf_tree` decision
record**, drafted in `tf_tree`, citing this spec; `helicoid`'s side is API parity, the envelope and
the bench evidence. `tf_tree`'s `LerpSlerp` is a `tf2`-compatibility emulation owned by `tf_tree`
([`0009`](./decisions/0009-what-helicoid-does-not-own.md)); `helicoid` does not reproduce its
small-angle LERP fallback.

## 1. Geodesic

**NORMATIVE.**

### 1.1 Definition

```rust
fn geodesic(x0: &Self, x1: &Self, t: S) -> Self;               // provided on LieGroup
```

$X(t) = X_0\,\mathrm{Exp}(t\,\mathrm{Log}(X_0^{-1} X_1))$ (`NUMERICS.md` §10). The provided body **is**
`reference::geodesic`; groups may override with a fast twin.

### 1.2 SE(3) fast twin

Port `tf_tree_math::dualquat::screw_pow` (one `atan2`, one `sin_cos`), generalized to `S: Real`.
Overrides `SE3::geodesic` only if `bench-gate` shows speedup over reference. Agreement tested via
`se3_geodesic_matches_reference` ($10^5$ random pairs).

Shipped as `NUMERICS.md` §10 states it ([`0054`](./decisions/0054-the-screw-twin-takes-two-arms-in-the-world-frame.md)): no switch
of its own and no series outside §4. Below `r`'s second switch, where `tf_tree_math`'s series arm
runs, the translation is the definition's closed form on §4's coefficients (GE.15(b)), so GE.12's
`0/0` at `‖v‖ = 0` and its `Dual` loss never arise.

### 1.3 `Product<SO3, R3>`

Its geodesic is slerp + lerp by the product structure — `tf2`'s semantics without `tf2`'s fallback.
Enables `tf_tree` to express `LerpSlerp` on `helicoid` types.

## 2. Jacobians and velocity

**NORMATIVE.** `geodesic_jacobians::<Right>(x0, x1, t) -> (J0, J1)` and
`geodesic_velocity(x0, x1) -> Tangent` per `NUMERICS.md` §10. The reference twin is `Dual` through
`reference::geodesic` (`geodesic_jacobians_match_reference`); `jacobians_match_dual_geodesic`
covers the fast twin.

## 3. Invariance tests

**NORMATIVE**, mirroring `tf_tree` `PHASE1.md` §3.4:

- **Left:** `geodesic(G·a, G·b, s) == G·geodesic(a, b, s)` for every group, both SE(3) implementations.
- **Right:** `geodesic(a·H, b·H, s) == geodesic(a, b, s)·H` for SE(3), both implementations.
- **`Product<SO3, R3>` under its own law is right-invariant** (`Ad_H = diag(R_H, I)` preserves `θ`).
- **`Product<SO3, R3>` fails right-invariance under SE(3) coupling**: fixed-seed test asserts
  `max_err > 1e-6`. Lerp is straight-line, while right-invariance forces screw coupling `ρ = J_l⁻¹(φ) t`
  (`docs/maths/geodesics.md` GE.5, [`0045`](./decisions/0045-two-phase-4-checks-cannot-be-taken-as-written.md)).

Left-invariance is measured at `‖t_G‖ ∈ {0, 1, 1e4}`.

## 4. Corpus additions

Function ids `so3_geodesic`, `se3_geodesic` (inputs $X_0, X_1, t$; reference `mp.expm` of $t\,d$
with geometric $\mathrm{Log}$). Group identities cross-check $\theta > 3.0$
([`0045`](./decisions/0045-two-phase-4-checks-cannot-be-taken-as-written.md)).
Strata: `geo:consecutive` (relative motion $\|d\| \in [10^{-9}, 10^{-3}]$, $\|t_0\|$ up to $10^4$),
`geo:generic`, `geo:near-pi` (relative rotation $\pi - 10^{-k}$);
$t \in \{0, 10^{-9}, 0.25, 0.5, 1 - 10^{-9}, 1\}$ plus uniform samples.
`tf_tree_math`'s `ScLerp` and `slerp` are the oracles.

## 5. What `helicoid` delivers to the `tf_tree` record

### 5.1 Parity table

| `tf_tree_math` | `helicoid` | Note |
|---|---|---|
| `Quat` `{w, x, y, z}` | `Quat<f64>` | Same field order ([`0002`](./decisions/0002-one-convention-for-a-stack-that-already-disagrees.md)) |
| `Iso3 { q, t }` | `SE3<f64>` | Same composition semantics; layout not a contract (D2) |
| `Vec3` | `Vec3<f64>` (`helicoid-linalg`) | |
| `exp_se3([f64; 6])` / `log_se3` | `SE3::exp(&Twist)` / `log` + `write_dense` | Both `[ω, v]` |
| `exp_so3` / `log_so3` | `SO3::exp` / `log` | Quaternion `atan2` form in both |
| `Twist { omega, v }` | `Twist<f64>` | `omega()`, `v()` |
| `Ad` in `twist.rs` | `SE3::adjoint()` (`SEn3Jac`) | |
| `ScLerp`, `reference::sclerp` | `SE3::geodesic`, `reference::geodesic` | |
| `dualquat::screw_pow` | fast twin of `SE3::geodesic` | §1.2 |
| `slerp` | `SO3::geodesic` | |
| `LerpSlerp` | stays in `tf_tree` | §0 |
| `quat_from_rot3` | `SO3::from_matrix` | Shepperd, closed form. `SO3::from_matrix` normalizes ([`0049`](./decisions/0049-the-boundary-is-what-removes-a-way-to-be-wrong.md)). |
| `Interp` trait | stays in `tf_tree` | `tf_tree`'s API, implemented over `helicoid` |

### 5.2 Envelope precondition

`helicoid` dominates `tf_tree_math` on every paired stratum in §5.1, or exceptions must be recorded
in `conformance/baseline/exceptions.toml` ([`0046`](./decisions/0046-explained-by-record-needs-a-record-to-point-at.md)).

### 5.3 Bench precondition

`helicoid`'s `exp`, `log`, `compose`, `geodesic` at `f64` are not slower than `tf_tree_math`'s by
`bench-gate`. Binding gate is `tf_tree`'s own `cargo xtask bench-gate` and `just interp-accuracy`.

### 5.4 What the `tf_tree` record must amend

- **D12** — thresholds become `helicoid`'s generated ones.
- **D13** — reference twins move with the code.
- **D14** — `helicoid` and `helicoid-linalg` add no new third-party crate.
- **`tf_tree_math`** — re-export shim for one `0.0.x` release, then removal.
- **Arena boundary** — `Iso3` conversion to `SE3<f64>` is a field move.

## 6. Definition of done

- [ ] §1–§4 implemented; invariance tests green, including the positive failure.
- [ ] Fast-twin decision recorded with its bench numbers.
- [ ] §5.2 and §5.3 evidence produced and linked from a `draft` record in `tf_tree`.

## Appendix: suggested implementation order

1. Provided `geodesic`, `reference::geodesic`, corpus ids and strata.
2. Invariance tests.
3. Jacobians and velocity with their twins.
4. SE(3) fast twin; its differential proptest; the bench decision.
5. Parity evidence; the `tf_tree` draft record.
