//! Lie groups, Jacobians, charts and geodesics.
//!
//! Implemented so far (`docs/PHASE3.md` §0.0): the traits [`Tangent`], [`Jac`] and [`LieGroup`],
//! [`Side`] with [`Left`] and [`Right`], the translation group [`Rn`], which proves the traits
//! compile generically, the quaternion [`Quat`], the SE_N(3) tangent [`SEn3Tangent`] with its
//! [`Twist`] converters, and the dual-matrix Jacobian [`SEn3Jac`] with its dense twins in
//! [`mod@reference`], the rotation group [`SO3`] with its tangent [`SO3Tangent`] and `Mat3`
//! Jacobian, the products [`Product`] with the block-diagonal [`ProductJac`], the side-typed
//! [`Gaussian`], [`sinc`], the one coefficient of the kernel that is public
//! (`docs/decisions/0062`), and the integrated exponentials of [`so3`] (`docs/decisions/0064`). Every
//! group is written against the traits and generic over the scalar `S: Real` of
//! `helicoid-linalg`.
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

// SO(3) is the kernel's first consumer (`exp_coeffs`, `jr_coeffs`, `jr_inv_coeff`, `log_ratio`),
// so the `cfg` the module carried until this point is gone: a default build reaches every grouped
// entry point through `SO3`. `q_coeffs` is SE_N(3)'s and stays gated inside until §5 lands.
mod chart;
mod coeffs;
mod dualmat;
mod gaussian;
mod product;
mod quat;
pub mod reference;
mod rn;
mod se3_charts;
mod sen3;
mod side;
pub mod so3;
#[cfg(test)]
mod so3_tests;
mod traits;

pub use chart::{Chart, LeftChart, Lifted, Manifold, RightChart, WithChart};
#[cfg(feature = "__sweep")]
#[doc(hidden)]
pub use coeffs::sweep as __sweep;
pub use coeffs::{sinc, sinc_value};
pub use dualmat::SEn3Jac;
pub use gaussian::Gaussian;
pub use product::{Product, ProductJac};
pub use quat::Quat;
pub use rn::{Rn, RnJac, RnTangent};
pub use se3_charts::{Decoupled, Screw, Se3Chart, TwistBlockJac, WorldTranslation};
pub use sen3::{SEn3, SEn3Tangent, Twist, SE23, SE3};
pub use side::{Left, Right, Side};
pub use so3::{SO3Tangent, SO3};
pub use traits::{Jac, LieGroup, Tangent};

/// The crate README, compiled as a doctest so its example cannot rot.
///
/// `cfg(doctest)`, so the module exists only while `cargo test --doc` is collecting: it is not in
/// the rustdoc front page (which is this file's own header) and not in the built library.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
mod readme {}

// `proptest` and the `Vec` the tests collect into need `std`; the library itself stays `no_std`.
#[cfg(test)]
extern crate std;

#[cfg(test)]
mod chart_tests;
#[cfg(test)]
mod dualmat_tests;
#[cfg(test)]
mod gamma_tests;
#[cfg(test)]
mod gaussian_tests;
#[cfg(test)]
mod heis_tests;
#[cfg(test)]
mod laws;
#[cfg(test)]
mod product_tests;
#[cfg(test)]
mod quat_tests;
#[cfg(test)]
mod rn_tests;
#[cfg(test)]
mod se3_charts_tests;
#[cfg(test)]
mod sen3_tests;
