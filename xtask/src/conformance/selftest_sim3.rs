//! The `0070` half of `--self-test` (`docs/PHASE1.md` §10): two Sim(3) defects a convention slip
//! makes, planted through the public API.
//!
//! - **The scale ignored in `V`**: `Exp`'s translation as SE(3)'s `J_l(φ) ρ`. It errs by about
//!   `|σ| ‖ρ‖ / 2`, so a stratum of `sim3_exp` must fail where its smallest `|σ|` puts that
//!   [`FACTOR`] times above the clean maximum, and is printed, not judged, below.
//! - **`J_r`'s scale column with its sign flipped**, `+Γ₂(−M)ρ` read as `−`: wrong by `2Γ₂ρ`, of
//!   the size of `ρ`, so every `sim3_jr` stratum must fail.
//!
//! Every judged stratum must also be one where `FACTOR` times the clean maximum is below `1/u`, and
//! each id keeps a judged stratum at each precision.

use std::fmt::Write;
use std::path::Path;

use helicoid::{Jac, LieGroup, SEn3Tangent, Sim3, Sim3Tangent, Tangent, SE3, SO3};
use helicoid_linalg::{Precision, Real, StridedMut};

use super::corpus::{self, Record};
use super::evaluate;
use super::metric::unit_bits;
use super::selftest::Report;
use super::subject::{Output, Registered, Subject};
use crate::shipped::{scalars, Helicoid};

const FACTOR: f64 = 1024.0;

#[derive(Clone, Copy)]
enum Planted {
    ScaleIgnored,
    ScaleColumn,
}

impl Planted {
    fn name(self) -> &'static str {
        match self {
            Planted::ScaleIgnored => "sim3-scale-ignored",
            Planted::ScaleColumn => "sim3-scale-column",
        }
    }

    fn id(self) -> &'static str {
        match self {
            Planted::ScaleIgnored => "sim3_exp",
            Planted::ScaleColumn => "sim3_jr",
        }
    }

    fn answer<S: Real + Into<f64>>(self, record: &Record) -> Output {
        let Some(tau) = scalars::<S>(record, "tau").filter(|t| t.len() == 7) else {
            return Output::new();
        };
        let tau = Sim3Tangent::read_dense(&tau);
        match self {
            Planted::ScaleIgnored => {
                let q = SO3::exp(&helicoid::SO3Tangent { phi: tau.phi }).quat();
                let se = SEn3Tangent {
                    phi: tau.phi,
                    rho: [tau.rho],
                };
                let x = SE3::exp(&se).translation();
                Output::from([
                    (
                        "q".to_string(),
                        [q.w, q.x, q.y, q.z].map(Into::into).to_vec(),
                    ),
                    ("x".to_string(), x.0.map(Into::into).to_vec()),
                    ("sigma".to_string(), vec![tau.sigma.into()]),
                ])
            }
            Planted::ScaleColumn => {
                let mut buf = [S::zero(); 49];
                Sim3::<S>::jr(&tau).write_dense(&mut StridedMut::col_major(&mut buf, 7, 7));
                // Column 6, rows 3..6: the `ρ` rows of the scale column, column-major.
                for r in 3..6 {
                    buf[6 * 7 + r] = -buf[6 * 7 + r];
                }
                Output::from([("J".to_string(), buf.map(Into::into).to_vec())])
            }
        }
    }
}

struct Plant(Planted);

impl Subject for Plant {
    fn name(&self) -> &str {
        self.0.name()
    }

    fn supports(&self, fn_id: &str) -> bool {
        fn_id == self.0.id()
    }

    fn eval(&self, _: &str, record: &Record, precision: Precision) -> Output {
        match precision {
            Precision::F64 => self.0.answer::<f64>(record),
            Precision::F32 => self.0.answer::<f32>(record),
        }
    }
}

fn maxima(
    dir: &Path,
    fn_id: &str,
    subject: Box<dyn Subject>,
    precision: Precision,
) -> Result<Vec<(String, f64)>, String> {
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
    let rows = evaluate(dir, &entries, &registered, precision)?.remove(0);
    Ok(rows.into_iter().map(|r| (r.stratum, r.max_u)).collect())
}

/// The smallest `|σ|` of each stratum of `fn_id` at `precision`.
fn smallest_sigma(
    dir: &Path,
    fn_id: &str,
    precision: Precision,
) -> Result<Vec<(String, f64)>, String> {
    let entry = corpus::manifest(dir)?
        .into_iter()
        .find(|e| e.fn_id == fn_id)
        .ok_or_else(|| format!("no {fn_id} corpus file"))?;
    let mut out: Vec<(String, f64)> = Vec::new();
    for r in corpus::read(dir, &entry)? {
        if r.is_f32_stratum() != (precision == Precision::F32) {
            continue;
        }
        let s = r
            .input("tau")
            .and_then(|t| t.get(6))
            .map_or(0.0, |s| s.abs());
        match out.iter_mut().find(|(n, _)| *n == r.stratum) {
            Some((_, m)) => *m = m.min(s),
            None => out.push((r.stratum.clone(), s)),
        }
    }
    Ok(out)
}

pub(super) fn check(dir: &Path) -> Result<Report, String> {
    let mut text = format!(
        "0070: Sim(3)'s scale ignored in V and J_r's scale column flipped, > {FACTOR} x the clean max\n"
    );
    let mut failures = Vec::new();
    for planted in [Planted::ScaleIgnored, Planted::ScaleColumn] {
        let id = planted.id();
        for precision in [Precision::F64, Precision::F32] {
            let u = (-(unit_bits(precision) as f64)).exp2();
            let clean = maxima(dir, id, Box::new(Helicoid), precision)?;
            let got = maxima(dir, id, Box::new(Plant(planted)), precision)?;
            let sigma = smallest_sigma(dir, id, precision)?;
            let mut judged = 0;
            for ((stratum, c), (_, p)) in clean.iter().zip(&got) {
                let s = sigma
                    .iter()
                    .find(|(n, _)| n == stratum)
                    .map_or(0.0, |&(_, s)| s);
                let reach = FACTOR * c.max(1.0);
                let visible = reach * u < 1.0
                    && match planted {
                        // The error is `|σ| ‖ρ‖ / 2` over a scale of order one.
                        Planted::ScaleIgnored => s / 4.0 >= reach * u,
                        Planted::ScaleColumn => true,
                    };
                let _ = writeln!(
                    text,
                    "  {} {id} {stratum} ({precision:?}): clean {c:.3}, planted {p:.3e}{}",
                    planted.name(),
                    if visible { "" } else { " (not judged)" }
                );
                if !c.is_finite() {
                    failures.push(format!("{id} {stratum}: the clean program reads {c}"));
                } else if visible {
                    judged += 1;
                    if p.is_nan() || *p <= reach {
                        failures.push(format!(
                            "{}: not detected on {id} {stratum} ({p} against {c})",
                            planted.name()
                        ));
                    }
                }
            }
            if judged == 0 {
                failures.push(format!(
                    "{}: {id} has no judged stratum at {precision:?}",
                    planted.name()
                ));
            }
        }
    }
    Ok(Report { text, failures })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;

    #[test]
    fn both_sim3_defects_are_detected_and_the_clean_program_is_finite() -> Result<(), String> {
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
