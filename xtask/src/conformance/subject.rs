//! Subjects (`docs/PHASE1.md` §5): anything the harness evaluates over the corpus in process.

use std::collections::BTreeMap;

use helicoid_linalg::Precision;

use super::corpus::Record;

/// A subject's answer: output field name to its flat, column-major binary64 values, in the order
/// of the reference's field.
pub(crate) type Output = BTreeMap<String, Vec<f64>>;

pub(crate) trait Subject {
    fn name(&self) -> &str;
    fn supports(&self, fn_id: &str) -> bool;
    /// The record carries its reference; a subject reads `record.inputs` only.
    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output;
}

/// A subject with the `subject_version` its result rows carry (the trait has none:
/// `docs/PHASE1.md` §5 lists three methods).
pub(crate) struct Registered {
    pub(crate) version: String,
    pub(crate) subject: Box<dyn Subject>,
}

/// The in-process subjects. Empty: `helicoid` arrives with Phase 3 and `seeded:<defect>` with
/// `docs/PHASE1.md` §10; the harness's own oracles are test subjects and never rows a bar reads.
pub(crate) fn registry() -> Vec<Registered> {
    Vec::new()
}
