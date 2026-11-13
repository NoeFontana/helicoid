//! Uniform cumulative B-splines on Lie groups (`docs/PHASE7.md`, `docs/NUMERICS.md` §21).
//!
//! The caller hands over a window of `K` controls in its own storage, the segment parameter `u` and
//! `1/h`; it picks the segment from its own time, so nothing here compares a float or owns a time
//! type (`0114`). [`Segment`] forms the increments once per segment; [`Weights`] the basis at one
//! `u`; [`Segment::kinematics`] the value, body velocity and body acceleration by Sommer et al.'s
//! recursion, and [`JacSegment::kinematics_jacobians`] their `O(K)` Jacobians on every control,
//! either side. One generic recursion serves every group: on a product it runs componentwise, on Rⁿ
//! it reduces to the non-cumulative weights.
//!
//! Orders `K = 2..=6` ([`ValidOrder`]); every table is exact integers (`0115`).

#![no_std]
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod gram;
mod jacobian;
mod order;
pub mod reference;
mod segment;
mod weights;

#[cfg(test)]
mod tests;

pub use gram::derivative_gram;
pub use jacobian::{ControlJacobians, JacSegment};
pub use order::{Order, ValidOrder};
pub use segment::{Kinematics, Segment};
pub use weights::Weights;
