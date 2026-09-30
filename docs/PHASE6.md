# helicoid — Phase 6 Implementation Specification: Interop, Determinism, locus-calib, 1.0

> **Companions:** `docs/API.md` §4–§5, `docs/PROJECT.md` D16, [`0009`](./decisions/0009-what-helicoid-does-not-own.md).

**Deliverable:** what a solver needs beyond `Chart` (ambient Jacobians for over-parameterized
storage), `mint` on the group types, the cross-target bit-identity gate as a required check,
locus-calib built on `helicoid` from its first commit, and the criteria that end the `0.0.x` line.
Sections marked **NORMATIVE** are requirements.

## 0.0 Implementation status

| Area | Status |
|---|---|
| `AmbientChart`; quaternion and SE(3) ambient Jacobians (§1) | Not started |
| `mint` for `Quat` (§2) | Not started |
| `just determinism` required in CI (§3) | Not started |
| locus-calib depends on `helicoid` (§4) | Not started |
| 1.0 criteria (§6) | Not met |

## 0. Non-goals and guardrails — read first

**NORMATIVE.** No solver, no cost function, no loss, no `Problem` type — those are the solver
crate's ([`0009`](./decisions/0009-what-helicoid-does-not-own.md)). No `ceres` or `faer` feature.

## 1. Ambient Jacobians

**NORMATIVE.** For solvers that store the ambient parameterization (Ceres `Manifold`, `gtsam`
retract-on-storage): `PlusJacobian` $= \partial\,\mathrm{retract}(\delta)/\partial\delta|_{\delta=0}$ in
ambient coordinates, `MinusJacobian` its left inverse.

```rust
pub trait AmbientChart<S: Real>: Chart<S, Self::Point> {
    type Point;
    const AMBIENT: usize;
    fn write_plus_jacobian(&self, out: &mut StridedMut<'_, S>);   // AMBIENT × DOF
    fn write_minus_jacobian(&self, out: &mut StridedMut<'_, S>);  // DOF × AMBIENT
}
```

Quaternion $(w, x, y, z)$, right chart $q \otimes \mathrm{Exp}(\delta)$, with $v = (x, y, z)$:

$$
P = \tfrac12\begin{bmatrix} -v^\top \\ wI + [v]_\times \end{bmatrix}
= \tfrac12\begin{bmatrix} -x & -y & -z \\ w & -z & y \\ z & w & -x \\ -y & x & w \end{bmatrix},\qquad
M = 2\,\big[\,-v \;\big|\; wI - [v]_\times\,\big]
= 2\begin{bmatrix} -x & w & z & -y \\ -y & -z & w & x \\ -z & y & -x & w \end{bmatrix},
$$

$MP = I_3$ for unit $q$. SE(3) as $(q, t)$, 7 ambient, tangent $[\varphi; \rho]$, for `Screw` and
`Decoupled` alike (they agree to first order):
$P_{SE3} = \begin{bmatrix} P & 0 \\ 0 & R \end{bmatrix}$ ($7\times6$),
$M_{SE3} = \begin{bmatrix} M & 0 \\ 0 & R^\top \end{bmatrix}$ ($6\times7$). Corpus ids
`so3_plus_jacobian`, `se3_plus_jacobian` (reference: `mp.diff` of the retraction in ambient
coordinates).

## 2. `mint` for groups

`Quat<S> ↔ mint::Quaternion<S>` (`{ v, s }`, by field name). `SO3`/`SE3` convert through their
parts; no `mint` type for tangents (their order is `helicoid`'s, API R3).

## 3. Determinism

**NORMATIVE.** `just determinism` evaluates every in-process subject over the full corpus on
x86_64-linux, aarch64-linux and `wasm32-wasip1` (wasmtime), serializes outputs as hex floats in
corpus order, and compares SHA-256 digests. **Any difference fails.** Required in CI from this phase
(wired since Phase 1). A difference is a D16 violation: find the non-`libm` transcendental, the
`mul_add`, or the target flag. NaN sign and payload are not compared (`0018`); a NaN output already
fails the corpus.

## 4. locus-calib

locus-calib depends on `helicoid` from its first commit: charts for every state, ambient Jacobians
if it uses Ceres, `Gaussian<_, _, Right, _>` for published extrinsics. A `helicoid` change requested
by locus-calib is a `helicoid` record, like any other consumer's.

## 5. The solver crate

The separate solver crate (unnamed; faer-based; tiers per
[`0009`](./decisions/0009-what-helicoid-does-not-own.md)) consumes `Manifold`, `Chart`,
`Jac::write_dense`, `AmbientChart`. Anything it needs that is not in this list is a `helicoid`
record that must pass the ownership test.

## 6. The 1.0 criteria

**NORMATIVE.** All of:

1. Every Phase 1–6 §0.0 row **Done**; `NUMERICS.md` §0 all **Ready** (Sim(3) included).
2. Envelope: every paired stratum dominated; baseline no-regress green; determinism green.
3. At least two consumers migrated and shipping on `helicoid` (`tf_tree` and one of omnisac /
   locus-tag).
4. `cargo-semver-checks` baseline recorded; `API.md` reviewed against
   [`0002`](./decisions/0002-one-convention-for-a-stack-that-already-disagrees.md) item by item.
5. No `draft` record cited as settled anywhere (`just lint`).

## 7. Definition of done

- [ ] §1–§4 implemented; ambient Jacobians in the corpus and dominated where oracles exist
      (Ceres' `QuaternionManifold` through a container runner).
- [ ] Determinism required and green.
- [ ] §6 evaluated; either 1.0 is cut, or the unmet items are listed in §0.0.

## Appendix: suggested implementation order

1. `AmbientChart`, quaternion, SE(3); corpus ids.
2. `mint` for `Quat`.
3. Make `determinism` required.
4. locus-calib scaffolding (in locus-calib).
5. 1.0 review.
