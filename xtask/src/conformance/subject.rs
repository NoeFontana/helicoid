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
    /// The version at `f32`, for a subject whose `f32` kernel is not its binary64 one.
    pub(crate) version_f32: String,
    /// Why the subject has no `f32` kernel, when it has none: a run at `f32` is then an error on
    /// every id it supports, never the answer of another kernel under this subject's name.
    pub(crate) no_f32: Option<String>,
    pub(crate) subject: Box<dyn Subject>,
    /// A seeded defect: evaluated only when a run names it (`--subject`) or in `--self-test`, so
    /// that a plain run's verdict and table are about real subjects.
    pub(crate) planted: bool,
}

impl Registered {
    /// The `subject_version` its rows at `precision` carry.
    pub(crate) fn version_at(&self, precision: Precision) -> &str {
        match precision {
            Precision::F64 => &self.version,
            Precision::F32 => &self.version_f32,
        }
    }
}

#[cfg(test)]
impl Registered {
    pub(crate) fn new(version: &str, subject: Box<dyn Subject>) -> Self {
        Self {
            version: version.to_string(),
            version_f32: version.to_string(),
            no_f32: None,
            subject,
            planted: false,
        }
    }
}

/// The in-process subjects: `helicoid` (`crate::shipped`) and the seeded kernels and defects
/// (`crate::seeded`). The harness's own oracles are test subjects and never rows a bar reads.
pub(crate) fn registry() -> Vec<Registered> {
    let mut all = vec![crate::shipped::registered()];
    all.extend(crate::seeded::registry());
    all
}
