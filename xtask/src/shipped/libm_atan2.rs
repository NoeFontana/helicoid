//! `helicoid:libm-atan2`: the shipped program at [`Libm`], whose `atan2` is the `libm` crate's, the
//! arrangement `0105` replaced with `helicoid_linalg::math::atan2`, and whose every other operation
//! is the shipped scalar's. Kept as a subject so `cargo xtask compare` can re-measure the choice
//! (`0078`).
//!
//! Planted, as `helicoid:host-std` is: a twin is not a candidate, and no bar reads its rows. It
//! answers both precisions, since every operation but `atan2` is the shipped scalar's own.

pub(crate) mod scalar;

use helicoid_linalg::Precision;

use super::{program, Helicoid};
use crate::conformance::corpus::Record;
use crate::conformance::subject::{Output, Registered, Subject};
use scalar::Libm;

/// The subject's name, and the value of `--subject`.
pub(crate) const NAME: &str = "helicoid:libm-atan2";

/// The shipped program at [`Libm`].
pub(crate) struct Twin;

impl Subject for Twin {
    fn name(&self) -> &str {
        NAME
    }

    fn supports(&self, fn_id: &str) -> bool {
        Helicoid.supports(fn_id)
    }

    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        program::<Libm<f64>, Libm<f32>>(fn_id, record, precision)
    }
}

pub(crate) fn registered() -> Registered {
    let version = format!("{}@libm-atan2", env!("CARGO_PKG_VERSION"));
    Registered {
        version: version.clone(),
        version_f32: version,
        no_f32: None,
        planted: true,
        subject: Box::new(Twin),
    }
}
