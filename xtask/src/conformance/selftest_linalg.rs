//! The `0056` half of `--self-test` (`docs/PHASE1.md` §10): one defect per metric `0056` adds,
//! planted into the shipped program's answer. Each must raise its stratum's maximum by more than
//! [`FACTOR`] over the clean program's, which must itself be finite. A metric that cannot see a
//! dropped root, a rotated eigenvector, a flipped mask, half a Newton step or a doubled derivative
//! is not measuring what `NUMERICS.md` §11 says it does.

use std::fmt::Write;
use std::path::Path;

use helicoid_linalg::Precision;

use super::corpus::{self, Record};
use super::evaluate;
use super::selftest::Report;
use super::subject::{Output, Registered, Subject};
use crate::shipped::Helicoid;

/// How far above the clean maximum a planted one must read. Every defect below is gross, 2^24 u
/// or more where the clean program reads a few `u`: the bar is that a metric sees it, not how well.
const FACTOR: f64 = 1024.0;

/// A defect, planted into the shipped answer of one id.
#[derive(Clone, Copy, Debug)]
pub(super) enum Defect {
    /// `solve_cubic`: the first valid slot's mask cleared (the root-set distance).
    DropsARoot,
    /// `eig3`: the first two columns of `V` rotated by `1e-9` in their plane (the gap weight).
    RotatesAnEigenvector,
    /// `chol`: the mask negated (mask agreement).
    FlipsTheMask,
    /// `quat_renormalize`: half the Newton step, `(q + q')/2`.
    TakesHalfAStep,
    /// `real_sqrt`: the derivative doubled.
    DoublesTheDerivative,
}

impl Defect {
    pub(super) const ALL: [Defect; 5] = [
        Defect::DropsARoot,
        Defect::RotatesAnEigenvector,
        Defect::FlipsTheMask,
        Defect::TakesHalfAStep,
        Defect::DoublesTheDerivative,
    ];

    /// `(id, stratum)` the defect must fire on.
    fn site(self) -> (&'static str, &'static str) {
        match self {
            Defect::DropsARoot => ("solve_cubic", "cubic:distinct"),
            Defect::RotatesAnEigenvector => ("eig3", "eig:random"),
            Defect::FlipsTheMask => ("chol_n3", "chol:indefinite"),
            Defect::TakesHalfAStep => ("quat_renormalize", "renorm:eta-2^-27"),
            Defect::DoublesTheDerivative => ("real_sqrt", "x:1e0"),
        }
    }

    fn name(self) -> &'static str {
        match self {
            Defect::DropsARoot => "drops-a-root",
            Defect::RotatesAnEigenvector => "rotates-an-eigenvector",
            Defect::FlipsTheMask => "flips-the-mask",
            Defect::TakesHalfAStep => "takes-half-a-step",
            Defect::DoublesTheDerivative => "doubles-the-derivative",
        }
    }

    fn plant(self, record: &Record, mut out: Output) -> Output {
        let get = |out: &mut Output, k: &str| out.get_mut(k).map(std::mem::take);
        match self {
            Defect::DropsARoot => {
                if let Some(valid) = out.get_mut("valid") {
                    if let Some(v) = valid.iter_mut().find(|v| v.to_bits() == 1f64.to_bits()) {
                        *v = 0.0;
                    }
                }
            }
            Defect::RotatesAnEigenvector => {
                if let Some(mut v) = get(&mut out, "V").filter(|v| v.len() == 9) {
                    let (s, c) = 1e-9f64.sin_cos();
                    for r in 0..3 {
                        let (a, b) = (v[r], v[3 + r]);
                        (v[r], v[3 + r]) = (c * a + s * b, c * b - s * a);
                    }
                    out.insert("V".into(), v);
                }
            }
            Defect::FlipsTheMask => {
                if let Some(valid) = out.get_mut("valid") {
                    valid.iter_mut().for_each(|v| *v = 1.0 - *v);
                }
            }
            Defect::TakesHalfAStep => {
                if let (Some(q), Some(input)) = (out.get_mut("q"), record.input("q")) {
                    q.iter_mut()
                        .zip(input)
                        .for_each(|(o, i)| *o = 0.5 * (*o + i));
                }
            }
            Defect::DoublesTheDerivative => {
                if let Some(d) = out.get_mut("d") {
                    d.iter_mut().for_each(|d| *d *= 2.0);
                }
            }
        }
        out
    }
}

/// The shipped subject with `defect` planted into its answers.
struct Planted(Defect);

impl Subject for Planted {
    fn name(&self) -> &str {
        self.0.name()
    }

    fn supports(&self, fn_id: &str) -> bool {
        fn_id == self.0.site().0
    }

    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        self.0
            .plant(record, Helicoid.eval(fn_id, record, precision))
    }
}

/// The stratum's maximum under `subject`.
fn max_u(dir: &Path, subject: Box<dyn Subject>, fn_id: &str, stratum: &str) -> Result<f64, String> {
    let entries: Vec<_> = corpus::manifest(dir)?
        .into_iter()
        .filter(|e| e.fn_id == fn_id)
        .collect();
    let registered = [Registered {
        version: "self-test".into(),
        version_f32: "self-test".into(),
        no_f32: None,
        subject,
        planted: true,
    }];
    let rows = evaluate(dir, &entries, &registered, Precision::F64)?.remove(0);
    rows.iter()
        .find(|r| r.stratum == stratum)
        .map(|r| r.max_u)
        .ok_or_else(|| format!("{fn_id} has no stratum {stratum}"))
}

pub(super) fn check(dir: &Path, defects: &[Defect]) -> Result<Report, String> {
    let mut text =
        format!("0056: one planted defect per metric, each > {FACTOR} x the clean max\n");
    let mut failures = Vec::new();
    for &d in defects {
        let (fn_id, stratum) = d.site();
        let clean = max_u(dir, Box::new(Helicoid), fn_id, stratum)?;
        let planted = max_u(dir, Box::new(Planted(d)), fn_id, stratum)?;
        let _ = writeln!(
            text,
            "  {fn_id} {stratum}: clean {clean:.3}, {} {planted:.3e}",
            d.name()
        );
        if !clean.is_finite() {
            failures.push(format!(
                "{fn_id} {stratum}: the clean program reads {clean}"
            ));
        } else if planted.is_nan() || planted <= FACTOR * clean.max(1.0) {
            failures.push(format!(
                "{}: not detected on {fn_id} {stratum} ({planted} against {clean})",
                d.name()
            ));
        }
    }
    Ok(Report { text, failures })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;

    #[test]
    fn every_planted_defect_fires_and_the_clean_program_does_not() -> Result<(), String> {
        let report = check(&corpus_dir()?, &Defect::ALL)?;
        assert!(
            report.failures.is_empty(),
            "{:?}\n{}",
            report.failures,
            report.text
        );
        Ok(())
    }
}
