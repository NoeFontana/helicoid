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
| `LieGroup::geodesic`, reference twin (§1) | **Done** for §1.1 and §1.3; §1.2's fast twin is the row below. `LieGroup::{geodesic, geodesic_velocity}` are provided methods and the body **is** `reference::geodesic` — which they call, so the two are one expression and a group's override has something to be measured against that it cannot shadow. `Product` delegates per factor, so `Product<SO3, Rn<3>>` is slerp and lerp by construction (§1.3). **`SO3::geodesic` is overridden** with GE.14's blend ([`0050`](./decisions/0050-the-geodesic-s-denominator-is-the-whole-domination-gap.md), implemented): six spellings of the one function GE.14 equates, scored over the committed corpus, put the whole domination gap in the blend's *denominator* — `sin α` recomputed from `α`, never `‖v‖`, which is the same number for a unit quaternion and not the same floating-point value. Only the recomputed form makes **both** endpoints exact, each weight being one number divided by itself at its own endpoint, and every maximum the `‖v‖` spelling scored sat on or beside an endpoint. Measured: **1.644 / 1.738 / 1.642** `u`, against the provided body's 1.572 / 2.721 / 2.429 and `tf_tree_math::slerp`'s 2.187 / 1.834 / 1.642 — **oracle #1 dominated on all three** where it won two, and `just envelope` 115 → **113**, its only remaining geodesic failure `se3_geodesic/geo:generic`, which `SEn3`'s own path reads and `SO3::geodesic` does not reach. The cost is 1.046× at `geo:consecutive`, the one stratum the provided body still wins — and the sentence below is why, holding *there and only there*:  GE.14 proves the provided body is shortest-arc slerp with `Log`'s sign rule, in the better-conditioned spelling — `atan2` where slerp has `acos`, whose slope is infinite where the quaternions are nearly equal — and the corpus now says so, `helicoid` 1.572 `u` against `tf_tree_math::slerp`'s 2.187 at `geo:consecutive`. Measured against the oracle on the other two strata it is behind: `so3_geodesic` 2.721 / 2.429 against 1.834 / 1.642 at `geo:generic` / `geo:near-pi`, and `se3_geodesic` 1.572 / 2.721 / 2.429 against `ScLerp`'s 2.336 / 2.502 / 3.253 — **three** domination failures, both `so3_geodesic` strata above and `se3_geodesic` at `geo:generic`. `se3_geodesic` read 2.556 / 3.240 / 3.466 until [`0048`](./decisions/0048-the-relative-transform-pair-earns-the-surface-dot-and-norm-do-not.md) made `rminus` the one-rotation `inv_mul`, which this id reads through `geodesic`: its three maxima now equal `so3_geodesic`'s to every digit, so the translation block is no longer the scored maximum on any geodesic stratum and what is left is the shared rotation path. That is the motivation for §1.2 restated in numbers: the provided body is `Log` then `Exp` where `ScLerp` is one `atan2` and one `sin_cos`, so the reference pays for roundings the screw form does not have. Verification: `laws::geodesic` in `laws_for!` for every instantiated group (`Rn`, `SO3`, `SE3`, `SE_2(3)`, `Heis` and the four products) at `f64`, `f32` and `Dual<f64, D>`, **seven** legs, each bounded apart in `laws::GEODESIC_LEGS`'s order — both endpoints, GE.2(c)'s symmetry, the velocity in GE.2(b)'s finite form (the row below has its figures), GE.4's left **and** right invariance, which holds on every group and not only SE(3), and the `twin` leg against `reference::geodesic` — with bounds twice the worst of 10^6 draws of a stream of its own, folded into no single number: a folded bound lets the tight legs rot behind the loose one. `t = 0` is the left endpoint **bit for bit** on SO(3) and SE_2(3), and since [`0050`](./decisions/0050-the-geodesic-s-denominator-is-the-whole-domination-gap.md) `t = 1` is the right endpoint bit for bit on SO(3) too — including past `π`, where `Log`'s flip makes it `−q₁`, the same rotation on the short arc — which is the property the override was taken for and which took that law leg from 7.160 `u` to the 1.118 floor. The laws' `gerr` cannot say either (its floor for two bitwise equal quaternions is 1.118 `u`) and three dedicated tests do — with the one exception their rustdoc names, a `-0.0` component of `x₀` coming back `+0.0`. The corpus scores it: `xtask/src/shipped.rs`'s `Geodesic` family answers both ids |
| SE(3) dual-quaternion fast twin (§1.2) | Not started |
| `geodesic_jacobians`, `geodesic_velocity` (§2) | Partial: `LieGroup::geodesic_velocity` ships (`x₁ ⊖_R x₀`, constant along the curve) and is checked by `laws::geodesic_legs`'s `velocity` leg in GE.2(b)'s finite form — `Log(γ(t)⁻¹γ(t+⅓)) = ⅓·velocity`, which ties the returned tangent to the curve where comparing it with `rminus` would be a tautology (`rminus` *is* its body). Worst, in `u`: 8.353 at SE_2(3), 7.909 at SE(3), 6.505 at SO(3), 3.538 at `Rn`. **Owed:** `geodesic_jacobians` and §2's `Dual`-through-`reference::geodesic` twin, with [`0043`](./decisions/0043-the-geodesic-jacobian-ships-the-cancellation-free-form.md)'s cancellation-free `J₀` and `Jac::scale`; GE.13(c) says the translation block needs a per-stratum tolerance of order `10²α⁻¹u` and no single-coefficient series removes it |
| Invariance tests (§3) | Partial: left and right invariance hold on **every** group — `γ(h x₀, h x₁, t) = h γ(x₀, x₁, t)` and `γ(x₀ h, x₁ h, t) = γ(x₀, x₁, t) h`, the second from `Log(h⁻¹Δh) = Ad_{h⁻¹}d` with no commutativity used — and are two legs of `laws::geodesic`, so every instantiated group gets them at three scalars. Worst, in `u`, over 10^6 `laws::Rng::shaped` draws through each group's own `measure_geodesic`, taken across `f64`, `f32` and `Dual` — the protocol the bounds are recorded from, which the figures previously quoted here did not state and no run reproduced: left / right **9.276 / 10.979** at SO(3) ([`0050`](./decisions/0050-the-geodesic-s-denominator-is-the-whole-domination-gap.md)), 8.502 / 9.360 at SE(3), 10.731 / 8.316 at SE_2(3), 4.648 / 5.036 at `Heis`, 4.031 / 4.031 at `Rn`. **Done** for `Product<SO3, Rn<3>>`, the group §1.3 names, which the four generic product instantiations do not cover (they are abelian or Heisenberg): `product_tests::so3_r3` runs the generic legs on a quaternion factor, asserts §1.3's slerp and lerp **bit for bit** against the two factor curves, and asserts the positive failure under the **SE(3) reading** of `(R, t)` not merely past `1e-6` but **equal to GE.5(b)'s closed form** `R₀M_s t_H`, so it cannot pass on an unrelated defect (`0045` item 4). Its fixture carries GE.5(b)'s inequality: `s = ½`, `θ = 1` about `z`, `t_H = (1,2,3)`, so `θ²‖t_{H⊥}‖ ≈ 2.2` against the `1e-5` the bar needs. `0045` item 5's bound is measured per `‖t_G‖ ∈ {0, 1, 1e4}` — **6.585, 5.812, 3.800** `u` since [`0050`](./decisions/0050-the-geodesic-s-denominator-is-the-whole-domination-gap.md) (5.866, 6.811, 3.800 before it, so the largest moved from `1` to `0`), the **smallest at 1e4** either way, which is the worry `0045` item 5 raises: the absolute cancellation in `a⁻¹G⁻¹Gb` does grow with the scale, but §11's metric divides by `max(‖Log b‖, 1)`, which grows with it too |
| Corpus id `se3_geodesic`, `so3_geodesic`; strata (§4) | **Done.** 180 records each, `conformance/corpus/{so3,se3}_geodesic.jsonl` (37 MB of the 50 MB budget, from 34), through `conformance/generate/gen/geodesic.py`: the geodesic of the two rotations the records denote, `gamma = X_0 Exp(t Log Delta)` with the **geometric** `Log` (the `so3_log` quaternion route, and `sen3_log_n1`'s $\mathsf V$ solve for the translation block) and `sen3_exp_n1`'s matrix exponential. `mp.logm` appears nowhere in the forward path and is the generation-time cross-check below $\theta = 3$, where the two agree to **8.2e-113** over the whole corpus (6.8e-121 at `geo:consecutive`); `geo:near-pi` is above that limit at every record, so it is checked by the group identities alone, which is what [`0045`](./decisions/0045-two-phase-4-checks-cannot-be-taken-as-written.md) item 2 states. Those identities run at **every** record, each costing one further evaluation of the reference at a different argument: the two endpoints, the symmetry $\gamma(X_0,X_1,t) = \gamma(X_1,X_0,1-t)$, and $\mathrm{Log}(X(t)^{-1}X(t+\tfrac13)) = \tfrac13 d$ — constant body velocity with no difference quotient. Strata, three, `PHASE1.md` §4.4: `geo:consecutive`, `geo:generic`, `geo:near-pi`, 6 pose pairs crossed with 10 parameters (§4's six and four uniform), $\lVert x_0\rVert$ cycling $1, 10^2, 10^4$ in each. Measured and recorded there: at $\lVert x_0\rVert = 10^4$ with $\theta(d) = 10^{-9}$ the stored $x_1 - x_0$ keeps about three digits — the relative motion is **below the ulp of the poses it is between**, which is `tf_tree`'s kilohertz edge and the reason §4 asks for that scale. Verification: 13 generator tests (`tests/test_geodesic.py`) — four hand-computed interpolations, GE.4(a)'s left invariance and the matrix route as properties no cross-check uses, the stated margins and band, and 198 mutants of the cross-checks of which the only ones that pass are the 36 that negate the **whole** quaternion, which is the same rotation; `just corpus-check` byte-compares a regeneration. The oracle answers both ids (`tf_tree_math`'s `slerp` and `ScLerp`), so §5.1's parity rows are measurable; no candidate answers them until §1 lands, which is what the empty `helicoid` rows mean |
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
| `quat_from_rot3` | `SO3::from_matrix` | Shepperd, closed form — and **not a parity row**: `quat_from_rot3` returns an un-normalized quaternion by its own rustdoc where `SO3::from_matrix` normalizes, so this row records a behaviour change with a number on it ([`0049`](./decisions/0049-the-boundary-is-what-removes-a-way-to-be-wrong.md) decision 4). Nothing breaks: `tf_tree_c::layout::read` rejects `|det R - 1| > 1e-6` upstream, so every reachable argument is already a rotation to that tolerance |
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
