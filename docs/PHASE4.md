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
| `LieGroup::geodesic`, reference twin (§1) | Not started |
| SE(3) dual-quaternion fast twin (§1.2) | Not started |
| `geodesic_jacobians`, `geodesic_velocity` (§2) | Not started |
| Invariance tests (§3) | Not started |
| Corpus id `se3_geodesic`, `so3_geodesic`; strata (§4) | Not started |
| Parity table and evidence for `tf_tree` (§5) | Not started |

## 0. Non-goals and guardrails — read first

**NORMATIVE.** `helicoid` does not change `tf_tree`. The migration is a **`tf_tree` decision
record**, drafted in `tf_tree`, citing this spec; `helicoid`'s side is API parity, the envelope and
the bench evidence. `tf_tree`'s `LerpSlerp` is a `tf2`-compatibility emulation owned by `tf_tree`
([`0009`](./decisions/0009-what-helicoid-does-not-own.md): it needs information only `tf_tree`
has — `tf2`'s behaviour); `helicoid` does not reproduce its small-angle LERP fallback.

## 1. Geodesic

**NORMATIVE.**

### 1.1 Definition

```rust
fn geodesic(x0: &Self, x1: &Self, t: S) -> Self;               // provided on LieGroup
```

$X(t) = X_0\,\mathrm{Exp}(t\,\mathrm{Log}(X_0^{-1} X_1))$ (`NUMERICS.md` §10). The provided body **is**
`reference::geodesic`; groups may override with a fast twin.

### 1.2 SE(3) fast twin

Port `tf_tree_math::dualquat::screw_pow` (one `atan2`, one `sin_cos`), generalized to `S: Real`
with every branch through `S::branch` and the coefficients through `coeffs`. It becomes
`SE3::geodesic`'s implementation **only if** `bench-gate` shows it faster than the reference on the
Phase 3 fixtures (`tf_tree`'s measurement says it is; this re-measures under `helicoid`'s
kernels). `se3_geodesic_matches_reference`: $10^5$ random pairs including near-identity relative
motion with far-from-origin $X_0$ and near-π relative rotation; the recorded agreement is quoted
against `tf_tree`'s 1e-14 NORMATIVE.

### 1.3 `Product<SO3, R3>`

Its geodesic is slerp + lerp by the product structure — `tf2`'s semantics without `tf2`'s
fallback. It exists so `tf_tree` can express `LerpSlerp` on `helicoid` types if its record decides
to.

## 2. Jacobians and velocity

**NORMATIVE.** `geodesic_jacobians::<Right>(x0, x1, t) -> (J0, J1)` and
`geodesic_velocity(x0, x1) -> Tangent` per `NUMERICS.md` §10. The reference twin is `Dual` through
`reference::geodesic` (`geodesic_jacobians_match_reference`); `jacobians_match_dual_geodesic`
covers the fast twin.

## 3. Invariance tests

**NORMATIVE**, mirroring `tf_tree` `PHASE1.md` §3.4:

- **Left:** `geodesic(G·a, G·b, s) == G·geodesic(a, b, s)` for every group, both SE(3)
  implementations.
- **Right:** `geodesic(a·H, b·H, s) == geodesic(a, b, s)·H` for SE(3), both implementations.
- **`Product<SO3, R3>` under its own law is right-invariant**, exactly — `Ad_H = diag(R_H, I)`
  preserves `θ` — and a test asserts that it is.
- **`Product<SO3, R3>` fails right-invariance, positively, under the SE(3) reading of `(R, t)`**
  (`a·H = (R_a R_H, R_a t_H + t_a)`): a fixed-seed test asserts `max_err > 1e-6`. **Do not "fix"
  it** — the lerp is a world-frame straight line, while right-invariance forces the screw coupling
  `ρ = J_l⁻¹(φ) t`, and fixing it gives `ScLerp`, a different function. The gap is
  `|μ_s(θ)|·‖t_{H⊥}‖` with `|μ_s| = ½s(1−s)θ² + O(θ⁴)` and `t_{H⊥}` the part of `t_H` orthogonal
  to the axis, so it **vanishes** at `s ∈ {0, 1}`, at `θ = 0` and for `t_H` parallel to the axis:
  the fixture takes `s = ½` and `θ²‖t_{H⊥}‖ ≳ 1e-5`, and carries that inequality in its comment so
  a later edit cannot make the assertion vacuous a second way
  (`docs/maths/geodesics.md` GE.5, [`0045`](./decisions/0045-two-phase-4-checks-cannot-be-taken-as-written.md)).

Bounds are recorded measurements from the first green run, then no-regress. Left-invariance is
measured at `‖t_G‖ ∈ {0, 1, 1e4}` and its bound set from those: `geodesic(G·a, G·b, s)` forms
`a⁻¹G⁻¹Gb`, which equals `a⁻¹b` only to rounding, and §4's strata put `‖t₀‖` at `1e4`.

## 4. Corpus additions

Function ids `so3_geodesic`, `se3_geodesic` (inputs $X_0, X_1, t$; reference `mp.expm` of $t\,d$
with $d$ the **geometric** $\mathrm{Log}$ — the quaternion `atan2` form and a $\mathsf V$ solve —
**not `mp.logm`**, which returns a complex, non-principal logarithm from $\theta = 3.03$ and so is
wrong by $O(1)$ on exactly the `geo:near-pi` stratum below; it is the generation-time cross-check
for $\theta \le 3.0$, where the two agree to $9.2\times10^{-41}$, and the group identities are the
check above it:
[`0045`](./decisions/0045-two-phase-4-checks-cannot-be-taken-as-written.md)).
Strata: `geo:consecutive` (relative motion $\|d\| \in [10^{-9}, 10^{-3}]$ — `tf_tree`'s
kilohertz-edge regime — with $\|t_0\|$ up to $10^4$), `geo:generic`, `geo:near-pi` (relative
rotation $\pi - 10^{-k}$); $t \in \{0, 10^{-9}, 0.25, 0.5, 1 - 10^{-9}, 1\}$ plus uniform samples.
`tf_tree_math`'s `ScLerp` and `slerp` are the oracles. `geo:near-pi` stops strictly below $\pi$
and records its margin: at $\pi$ the two preimages give geodesics $O(1)$ apart, so the stratum
measures conditioning, not agreement. $t = 0$ and $t = 1$ stay and are exact rows — at $t = 0$ the
answer is $X_0$ with zero deviation.

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
| `quat_from_rot3` | `SO3::from_matrix` | Shepperd, closed form |
| `Interp` trait | stays in `tf_tree` | `tf_tree`'s API, implemented over `helicoid` |

### 5.2 Envelope precondition

`helicoid` dominates `tf_tree_math` on **every** paired stratum of every function in §5.1. A
stratum where `tf_tree_math` wins blocks the migration until fixed or explained by record — the
explanation being a row in `conformance/baseline/exceptions.toml` that names the `ready` record
([`0046`](./decisions/0046-explained-by-record-needs-a-record-to-point-at.md)), so an exception is
auditable and cannot outlive the defect it describes.

### 5.3 Bench precondition

`helicoid`'s `exp`, `log`, `compose`, `geodesic` at `f64` are not slower than `tf_tree_math`'s by
`bench-gate` (near-identity fixture included). This is evidence for the `tf_tree` record; the
**binding** gate is `tf_tree`'s own `cargo xtask bench-gate` (its `0013` rows) and
`just interp-accuracy`.

### 5.4 What the `tf_tree` record must amend

- **D12** — thresholds become `helicoid`'s generated ones; "θ < 0.1, four terms" becomes the
  recorded prior.
- **D13** — reference twins for the moved routines move with the code; `tf_tree` keeps its
  `tf2`-facing twins.
- **D14** — the dependency budget gains `helicoid` and `helicoid-linalg`, which add **no new
  third-party crate** (`libm` is already there); `bytemuck` stays `tf_tree`'s, for its arena
  records.
- **`tf_tree_math`'s fate** — suggested: a re-export shim for one `0.0.x` release, then removal.
  `tf_tree`'s call.
- **The arena boundary** — the arena keeps its own `Iso3` records (`tf_tree` D4, API R4); the
  conversion to `SE3<f64>` is a field move, costed by `tf_tree`'s bench.

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
