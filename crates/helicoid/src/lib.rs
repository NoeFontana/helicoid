//! Lie groups, Jacobians, charts and geodesics.
//!
//! Implemented so far (`docs/PHASE3.md` §0.0): the traits [`Tangent`], [`Jac`] and [`LieGroup`],
//! [`Side`] with [`Left`] and [`Right`], the translation group [`Rn`], which proves the traits
//! compile generically, the quaternion [`Quat`], the SE_N(3) tangent [`SEn3Tangent`] with its
//! [`Twist`] converters, and the dual-matrix Jacobian [`SEn3Jac`] with its dense twins in
//! [`mod@reference`]. Every group is written against the traits and generic over the scalar
//! `S: Real` of `helicoid-linalg`.
//!
//! # Conventions
//!
//! 1. **Quaternions:** Hamilton, stored `w` first (`Quat { w, x, y, z }`), active rotations.
//! 2. **Composition:** `a * b` is `T_a_x · T_x_b`; `X * p` maps a point from the child frame into the
//!    parent frame.
//! 3. **Tangents are rotation-first:** `[φ; ρ₁; …; ρ_N]` for SE_N(3), `[ω; v]` for a twist, `[θ; ρ]`
//!    for SE(2), `[φ; ρ; σ]` for Sim(3). Named fields everywhere; the order is visible only through
//!    `write_dense`/`read_dense`.
//! 4. **Right perturbation is the default; both sides are first-class and always spelled**
//!    (`rplus`, `lplus`, `rminus`, `lminus`, `Side`). No `Add`/`Sub` for ⊕/⊖.
//! 5. **`Log`'s canonical branch** is $\theta \in [0, \pi]$ via the `w ≥ 0` flip; at `w = +0` it is a
//!    function of the quaternion's sign (`NUMERICS.md` §3.2).
//! 6. **Converters for every other order**, named after it (API R3): `Quat::from_xyzw`, `to_xyzw`,
//!    `from_jpl` (Hamilton $(w, -x, -y, -z)$ from JPL $(x, y, z, w)$, same matrix — Sommer et al.
//!    2018), `Twist::from_translation_first`, `to_translation_first`.
//! 7. `NUMERICS.md` states every formula already permuted into these conventions.

#![no_std]
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod dualmat;
mod quat;
pub mod reference;
mod rn;
mod sen3;
mod side;
mod traits;

pub use dualmat::SEn3Jac;
pub use quat::Quat;
pub use rn::{Rn, RnJac, RnTangent};
pub use sen3::{SEn3Tangent, Twist};
pub use side::{Left, Right, Side};
pub use traits::{Jac, LieGroup, Tangent};

// `proptest` and the `Vec` the tests collect into need `std`; the library itself stays `no_std`.
#[cfg(test)]
extern crate std;

#[cfg(test)]
mod dualmat_tests;
#[cfg(test)]
mod heis_tests;
#[cfg(test)]
mod laws;
#[cfg(test)]
mod quat_tests;
#[cfg(test)]
mod rn_tests;
#[cfg(test)]
mod sen3_tests;
