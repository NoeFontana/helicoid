//! The `0064` half of `--self-test` (`docs/PHASE1.md` §10): `Γ₂` built from `Γ₁`'s coefficients,
//! `½I + aW + bW²` for `½I + bW + dW²`, the off-by-one in `m` that GG.2(a)'s `σ_{m+1}`, `σ_{m+2}`
//! makes easy, planted as `so3::gamma1(φ) − ½I` over `so3_gamma2`.
//!
//! The planted answer is wrong by `(a − b)W + (b − d)W²`, about `θ/3` at small `θ`, so a metric
//! with a floor of 1 sees it exactly where `θ/3` is above rounding. A stratum must therefore fail by
//! more than [`FACTOR`] times the clean maximum when its smallest `θ` is at least [`VISIBLE`] `u`,
//! at each precision over its own strata; below that (`theta:exact0`, `theta:subnormal`, the
//! smallest decades) no metric of a floor of 1 can see it, and the clean and planted maxima are
//! printed and not judged.

use std::fmt::Write;
use std::path::Path;

use helicoid::{Jac, LieGroup, Tangent, SO3};
use helicoid_linalg::{Mat3, Precision, Real, StridedMut};

use super::corpus::{self, exact_f32, Record};
use super::evaluate;
use super::metric::unit_bits;
use super::selftest::Report;
use super::subject::{Output, Registered, Subject};
use crate::shipped::Helicoid;

/// How far above the clean maximum a planted one must read, as `selftest_linalg`'s.
const FACTOR: f64 = 1024.0;
/// The smallest `θ`, in `u`, at which the defect must be seen: `θ/3` is then `4096 u`, `FACTOR`
/// times a clean maximum of `4 u`.
const VISIBLE: f64 = 3.0 * FACTOR * 4.0;

/// `so3::gamma1(φ) − ½I` at `S`, `φ` read exactly; nothing for any other id.
fn planted<S: Real + Into<f64>>(record: &Record) -> Output {
    let Some(phi) = phi_of(record).and_then(|p| exact::<S>(&p)) else {
        return Output::new();
    };
    let tau = <SO3<S> as LieGroup<S>>::Tangent::read_dense(&phi);
    let g = helicoid::so3::gamma1(&tau) - Mat3::<S>::identity().scale(S::lit(0.5));
    let mut buf = [S::zero(); 9];
    Jac::<S, <SO3<S> as LieGroup<S>>::Tangent>::write_dense(
        &g,
        &mut StridedMut::col_major(&mut buf, 3, 3),
    );
    Output::from([("G".to_string(), buf.map(Into::into).to_vec())])
}

fn phi_of(record: &Record) -> Option<[f64; 3]> {
    record.input("phi")?.first_chunk::<3>().copied()
}

fn exact<S: Real>(phi: &[f64; 3]) -> Option<[S; 3]> {
    let one = |x: f64| match S::PRECISION {
        Precision::F64 => Some(S::lit(x)),
        Precision::F32 => exact_f32(x).map(|x| S::lit(f64::from(x))),
    };
    Some([one(phi[0])?, one(phi[1])?, one(phi[2])?])
}

/// The shipped subject with `Γ₂` answered from `Γ₁`'s coefficients.
struct FromGamma1;

impl Subject for FromGamma1 {
    fn name(&self) -> &str {
        "gamma2-from-gamma1"
    }

    fn supports(&self, fn_id: &str) -> bool {
        fn_id == "so3_gamma2"
    }

    fn eval(&self, _: &str, record: &Record, precision: Precision) -> Output {
        match precision {
            Precision::F64 => planted::<f64>(record),
            Precision::F32 => planted::<f32>(record),
        }
    }
}

/// `(stratum, max_u)` of `subject` over `so3_gamma2` at `precision`.
fn maxima(
    dir: &Path,
    subject: Box<dyn Subject>,
    precision: Precision,
) -> Result<Vec<(String, f64)>, String> {
    let entries: Vec<_> = corpus::manifest(dir)?
        .into_iter()
        .filter(|e| e.fn_id == "so3_gamma2")
        .collect();
    let registered = [Registered {
        version: "self-test".into(),
        version_f32: "self-test".into(),
        no_f32: None,
        subject,
        planted: true,
    }];
    let rows = evaluate(dir, &entries, &registered, precision)?.remove(0);
    Ok(rows.into_iter().map(|r| (r.stratum, r.max_u)).collect())
}

/// The smallest `θ` of each stratum of `so3_gamma2` at `precision`, in that precision's `u`.
fn smallest_theta(dir: &Path, precision: Precision) -> Result<Vec<(String, f64)>, String> {
    let entry = corpus::manifest(dir)?
        .into_iter()
        .find(|e| e.fn_id == "so3_gamma2")
        .ok_or("no so3_gamma2 corpus file")?;
    let u = (-(unit_bits(precision) as f64)).exp2();
    let mut out: Vec<(String, f64)> = Vec::new();
    for r in corpus::read(dir, &entry)? {
        if r.is_f32_stratum() != (precision == Precision::F32) {
            continue;
        }
        let theta = phi_of(&r).map_or(0.0, |p| p.iter().map(|x| x * x).sum::<f64>().sqrt()) / u;
        match out.last_mut() {
            Some((s, m)) if *s == r.stratum => *m = m.min(theta),
            _ => out.push((r.stratum.clone(), theta)),
        }
    }
    Ok(out)
}

pub(super) fn check(dir: &Path) -> Result<Report, String> {
    let mut text = format!(
        "0064: Γ₂ from Γ₁'s coefficients over so3_gamma2, > {FACTOR} x the clean max where θ ≥ {VISIBLE} u\n"
    );
    let mut failures = Vec::new();
    for precision in [Precision::F64, Precision::F32] {
        let clean = maxima(dir, Box::new(Helicoid), precision)?;
        let planted = maxima(dir, Box::new(FromGamma1), precision)?;
        let theta = smallest_theta(dir, precision)?;
        let mut judged = 0;
        for ((stratum, c), (_, p)) in clean.iter().zip(&planted) {
            let t = theta
                .iter()
                .find(|(s, _)| s == stratum)
                .map_or(0.0, |&(_, t)| t);
            let visible = t >= VISIBLE;
            let _ = writeln!(
                text,
                "  {stratum} ({precision:?}): clean {c:.3}, planted {p:.3e}{}",
                if visible {
                    ""
                } else {
                    " (below θ = VISIBLE u, not judged)"
                }
            );
            if !c.is_finite() {
                failures.push(format!("so3_gamma2 {stratum}: the clean program reads {c}"));
            } else if visible {
                judged += 1;
                if p.is_nan() || *p <= FACTOR * c.max(1.0) {
                    failures.push(format!(
                        "gamma2-from-gamma1: not detected on so3_gamma2 {stratum} ({p} against {c})"
                    ));
                }
            }
        }
        if judged == 0 {
            failures.push(format!("so3_gamma2 has no judged stratum at {precision:?}"));
        }
    }
    Ok(Report { text, failures })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;

    #[test]
    fn gamma2_from_gamma1_is_detected_and_the_clean_program_is_finite() -> Result<(), String> {
        let report = check(&corpus_dir()?)?;
        assert!(
            report.failures.is_empty(),
            "{:?}\n{}",
            report.failures,
            report.text
        );
        Ok(())
    }
}
