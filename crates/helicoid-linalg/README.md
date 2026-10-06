# helicoid-linalg

The scalar model and fixed-size linear algebra under [`helicoid`](https://crates.io/crates/helicoid)
— `no_std`, no allocation, `#![forbid(unsafe_code)]`, and deterministic across x86_64, aarch64 and
wasm32.

```rust
use helicoid_linalg::{Dual, Matrix, Real, Vector};

// Generic over the scalar: the same code runs at `f64`, `f32` or a forward-mode dual number.
fn scaled_norm_sq<S: Real>(v: &Vector<S, 3>, k: S) -> S {
    v.0.iter().fold(S::zero(), |acc, &x| acc + (x * k) * (x * k))
}

let v = Vector([1.0_f64, 2.0, 3.0]);
assert!((scaled_norm_sq(&v, 2.0) - 56.0).abs() < 1e-12);

// `Dual<f64, 1>` carries the derivative through the same function, with no second implementation.
let d = Vector([
    Dual::<f64, 1>::variable(1.0, 0),
    Dual::constant(2.0),
    Dual::constant(3.0),
]);
// d/dx of 4(x² + 4 + 9) at x = 1 is 8x = 8.
assert!((scaled_norm_sq(&d, Dual::constant(2.0)).d[0] - 8.0).abs() < 1e-12);

let _m: Matrix<f64, 3, 3> = Matrix::identity();
```

## Why `Real` has no `PartialOrd`

Deliberately. Generic numeric code cannot write `if theta < eps`, because at a dual number that
branch silently drops the derivative of the arm not taken, and at a SIMD lane it does not typecheck
at all. Comparisons return `S::Mask`; a branch goes through `S::branch`/`S::select`, with the
safe-argument pattern inside the exact arm. The type system carries the rule so review does not have
to.

## What ships today

`Real` and `Mask`; `Dual<S, D>` forward-mode duals; `Vector<S, N>`, `Point<S, N>`, `Matrix<S, R, C>`
and the `Vec3`/`Mat3` aliases; strided views; `chol`, `solve_cubic`, `eig3`; optional [`mint`]
interop behind the `mint` feature.

**Not yet:** `svd3`. Each spec's §0.0 table in
[`docs/`](https://github.com/NoeFontana/helicoid/tree/main/docs) is the source of truth over any
prose, this file included.

## Stability

`0.0.x`: **every release may break every other.** MSRV 1.87. Dependencies: `libm`, and `mint` only
when its feature is on.

License: MIT OR Apache-2.0.

[`mint`]: https://crates.io/crates/mint
