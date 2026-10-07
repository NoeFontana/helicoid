# helicoid — Phase 6 Implementation Specification: Interop, Determinism, Calibration, 1.0

> **Companions:** `docs/API.md` §4–§5, `docs/PROJECT.md` D16, [`0009`](./decisions/0009-what-helicoid-does-not-own.md).

**Deliverable:** ambient Jacobians for over-parameterized storage, `mint` on group types,
cross-target bit-identity gate, calibration tooling support, and 1.0 release criteria.
Sections marked **NORMATIVE** are requirements.

## 0.0 Implementation status

| Area | Status |
|---|---|
| `AmbientChart`; quaternion and SE(3) ambient Jacobians (§1) | Not started |
| `mint` for `Quat` (§2) | Not started |
| `just determinism` required in CI (§3) | Not started |
| Calibration library depends on `helicoid` (§4) | Not started |
| 1.0 criteria (§6) | Not met |

## 0. Non-goals and guardrails — read first

**NORMATIVE.** No solver, cost function, loss, or `Problem` type — owned by solver crates
([`0009`](./decisions/0009-what-helicoid-does-not-own.md)). No `ceres` or `faer` feature.

## 1. Ambient Jacobians

**NORMATIVE.** For solvers storing ambient parameterizations (e.g. Ceres `Manifold`):
`PlusJacobian` $= \partial\,\mathrm{retract}(\delta)/\partial\delta|_{\delta=0}$, `MinusJacobian` its left inverse.

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
P = \tfrac12\begin{bmatrix} -v^\top \\ wI + [v]_\times \end{bmatrix},\qquad
M = 2\,\big[\,-v \;\big|\; wI - [v]_\times\,\big],
$$

$MP = I_3$ for unit $q$. SE(3) ($7\times6$ ambient $P_{SE3}$, $6\times7$ $M_{SE3}$).
Corpus ids: `so3_plus_jacobian`, `se3_plus_jacobian`.

## 2. `mint` for groups

`Quat<S> ↔ mint::Quaternion<S>`. `SO3`/`SE3` convert via components. Tangents keep internal order.

## 3. Determinism

**NORMATIVE.** `just determinism` evaluates in-process subjects across x86_64, aarch64, and
`wasm32-wasip1`, comparing SHA-256 digests of hex-float outputs. Any difference fails (D16).

## 4. Calibration integration

Downstream calibration tools consume `helicoid`: charts for states, ambient Jacobians,
and `Gaussian<_, _, Right, _>` for extrinsics.

## 5. The solver crate

Consumes `Manifold`, `Chart`, `Jac::write_dense`, `AmbientChart`.

## 6. The 1.0 criteria

**NORMATIVE.**
1. Every Phase 1–6 §0.0 row **Done**; `NUMERICS.md` §0 all **Ready**.
2. Envelope: every paired stratum dominated; baseline and determinism green.
3. At least two consumers migrated and shipping on `helicoid` (`tf_tree` and one additional consumer).
4. `cargo-semver-checks` baseline recorded; `API.md` reviewed against
   [`0002`](./decisions/0002-one-convention-for-a-stack-that-already-disagrees.md).
5. No `draft` records cited as settled.

## 7. Definition of done

- [ ] §1–§4 implemented; ambient Jacobians in corpus and dominated.
- [ ] Determinism required and green.
- [ ] §6 evaluated for 1.0 release.

## Appendix: suggested implementation order

1. `AmbientChart`, quaternion, SE(3); corpus ids.
2. `mint` for `Quat`.
3. Make `determinism` required.
4. Calibration scaffolding.
5. 1.0 review.
