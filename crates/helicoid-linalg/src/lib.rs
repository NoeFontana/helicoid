//! The leaf of the `helicoid` workspace: the scalar model and fixed-size linear algebra.
//!
//! Implemented so far: the scalar model (`docs/PHASE2.md` §2), i.e. [`Mask`], [`Real`], [`Blend`]
//! and [`Precision`], with `f64` and `f32` as `Real` (mask `bool`, every transcendental through
//! `libm`). Numeric code is written once, generic over `S: Real`; a comparison yields `S::Mask`
//! and control flow goes through [`Real::branch`] and [`Real::select`].
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

mod float;
mod real;

pub use real::{Blend, Mask, Precision, Real};

#[cfg(test)]
mod tests;
