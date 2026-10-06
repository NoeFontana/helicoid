# helicoid

Lie groups for robotics perception, with structured Jacobians and geodesics — `no_std`, no
allocation, `#![forbid(unsafe_code)]`, generic over the scalar (`f64`, `f32`, forward-mode dual
numbers), bit-identical across x86_64, aarch64 and wasm32.

```rust
use helicoid::{Jac, LieGroup, Right, Tangent, Twist, SE3};
use helicoid_linalg::Vector;

// Tangents are rotation-first: `phi`, then one `rho` per column of the group element.
let a = SE3::<f64>::exp(&Twist { phi: Vector([0.1, -0.2, 0.3]), rho: [Vector([1.0, 2.0, 3.0])] });
let b = a.rplus(&Twist { phi: Vector([0.0, 0.0, 0.05]), rho: [Vector([0.5, 0.0, 0.0])] });

// Every perturbation names its side: `b.rminus(&a)` is `Log(a⁻¹b)`, the body-frame delta.
let delta = b.rminus(&a);
// The geodesic is `a ⊕ t·delta`, one expression over those side-named primitives.
let mid = SE3::geodesic(&a, &b, 0.5);
// Jacobians are structured: `apply` pushes a tangent through without materializing a 6×6.
let (j_a, _j_b) = a.compose_jacobians::<Right>(&b);
let pushed = j_a.apply(&delta);

assert!(delta.dot(&delta) > 0.0 && pushed.dot(&pushed) > 0.0);
assert!(mid.rminus(&a).dot(&delta) > 0.0);
```

## What it is

- **One convention, stated and enforced.** Hamilton quaternions, `w` first, active rotations;
  `a * b` is `T_a_x · T_x_b`; tangents are rotation-first `[φ; ρ₁; …; ρ_N]`; right perturbation is
  the default and **every** perturbation names its side (`rplus`/`lplus`/`rminus`/`lminus`). There
  is no `Add`/`Sub` impl for ⊕/⊖, so a side is never implied.
- **Jacobians are structured, not dense.** A group's Jacobian is its own type (`G::Jac`) with
  `mul`, `inverse`, `apply`, `sandwich`; a dense matrix is produced only on request through
  `write_dense`. No method returns `[[S; 6]; 6]`.
- **Switch points are measured, not typed.** Every series/exact boundary and every series
  coefficient in the crate is generated from a committed sweep over an mpmath reference; a hand edit
  fails the repository's own check. No `if theta < 1e-8` appears anywhere.
- **Determinism is a property, not a hope.** Every transcendental goes through the `libm` crate; no
  `mul_add`, no `target-cpu`, no fast-math.
- **Measured per stratum, against oracles.** Each routine is scored in units of `u` against a
  110-digit mpmath corpus and against an envelope of independent implementations, per function, per
  stratum, per precision, on the **maximum** — never a mean.

## What ships today

`Quat`, `SO3`, `SEn3<S, N>` (with `SE3` and `SE23`), `Rn`, `Product`, the `Tangent`/`Jac`/
`LieGroup`/`Side` traits, the coefficient kernel, side and action Jacobians, geodesics and their
reference twins.

**Not yet:** `SO2`/`SE2`, `Sim(3)`, `S²`, charts as types, Γ functions, `Gaussian`. Each is
specified in [`docs/`](https://github.com/NoeFontana/helicoid/tree/main/docs) and none is
implemented; a spec's §0.0 table is the source of truth over any prose, this file included.

## Stability

`0.0.x`: **every release may break every other.** Nothing here is a stability promise yet; the
criteria for 1.0 are in `docs/PHASE6.md` §6.

MSRV 1.87. Dependencies: `helicoid-linalg` and `libm`, and nothing else.

License: MIT OR Apache-2.0.
