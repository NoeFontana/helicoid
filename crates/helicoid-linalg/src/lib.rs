//! The leaf of the `helicoid` workspace: the scalar model and fixed-size linear algebra.
//!
//! Implemented so far: the scalar model (`docs/PHASE2.md` §2), i.e. [`Mask`], [`Real`], [`Blend`]
//! and [`Precision`], with `f64` and `f32` as `Real` (mask `bool`, every transcendental through
//! `libm`), forward-mode [`Dual`] numbers (§3), also a `Real`, and the fixed-size [`Vector`],
//! [`Point`] and column-major [`Matrix`] with [`hat`]/[`vee`], [`Mat3::inverse_adj`] and the
//! Cholesky factorization [`chol`] with [`solve_lower`]/[`solve_upper`] (§4).
//! Numeric code is written once, generic over `S: Real`; a comparison yields `S::Mask` and control
//! flow goes through [`Real::branch`] and [`Real::select`].
//!
//! ```
//! use helicoid_linalg::Real;
//!
//! /// `sqrt(|x|)` with a branch that evaluates one arm for `f64`.
//! fn root<S: Real>(x: S) -> S {
//!     let neg = x.lt(S::zero());
//!     S::branch(neg, || (-x).sqrt(), || x.sqrt())
//! }
//! assert_eq!(root(-4.0_f64).value_f64().to_bits(), 2.0_f64.to_bits());
//! ```

#![no_std]
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod chol;
mod dual;
mod float;
mod matrix;
mod point;
mod real;
mod skew;
mod vector;

pub use chol::{chol, solve_lower, solve_upper};
pub use dual::Dual;
pub use matrix::{Mat2, Mat3, Matrix};
pub use point::{Point, Point2, Point3};
pub use real::{Blend, Mask, Precision, Real};
pub use skew::{hat, vee};
pub use vector::{Vec2, Vec3, Vector};

// `proptest` and the `format!` in its macros need `std`; the library itself stays `no_std`.
#[cfg(test)]
extern crate std;

#[cfg(test)]
mod chol_tests;
#[cfg(test)]
mod dual_tests;
#[cfg(test)]
mod linalg_tests;
#[cfg(test)]
mod tests;
