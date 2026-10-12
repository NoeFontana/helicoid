//! The `0063` half of `--self-test` (`docs/PHASE1.md` §10): the two readings of `NUMERICS.md` §8
//! that `0063` settled, planted over `s2_retract` and `s2_local`.
//!
//! - **`ς` by the sign bit**, `copysign(1, n_z)` for decision 1's comparison: at `n_z = -0` the
//!   oriented basis is rotated by `T°` (`docs/maths/charts.md` CH.8(f), CH.11(d)), the identity only
//!   at azimuth `π/2 (mod π)`, so the `s2:nz0` strata, whose odd records hold `-0`, must fail
//!   and every other stratum must read what the clean program reads. Planted through the public
//!   API: the chart at `n_z = -2⁻¹⁴⁰` has `ν_z = -1` exactly, the sign-bit basis at `-0`.
//! - **§8's ratio as worded**, `r(s², w)` for `½ r(s², w)` (decision 3): `local` returns `2δ`.
//!
//! A stratum is judged where [`FACTOR`] times its clean maximum is below a total loss (`1/u`): near
//! the antipode at `f32`, `local` is at its conditioning (`u/(π - α)` is `0.1` at `π - α = 10⁻⁶`,
//! CH.11(e)), so no defect can read `FACTOR` times above it there. There the maxima are printed and
//! not judged, and each id must keep a judged stratum at each precision.

use std::fmt::Write;
use std::path::Path;

use helicoid::{Chart, S2Chart, S2};
use helicoid_linalg::{Precision, Real, Vector};

use super::corpus::{self, Record};
use super::evaluate;
use super::metric::unit_bits;
use super::selftest::Report;
use super::subject::{Output, Registered, Subject};
use crate::shipped::{s2_of, Helicoid, Sphere};

/// How far above the clean maximum a planted one must read, as `selftest_gamma`'s.
const FACTOR: f64 = 1024.0;

/// The two plantings.
#[derive(Clone, Copy)]
enum Planted {
    SignBit,
    FullRatio,
}

impl Planted {
    fn name(self) -> &'static str {
        match self {
            Planted::SignBit => "s2-sign-bit",
            Planted::FullRatio => "s2-full-ratio",
        }
    }

    fn ids(self) -> &'static [&'static str] {
        match self {
            Planted::SignBit => &["s2_retract", "s2_local"],
            Planted::FullRatio => &["s2_local"],
        }
    }

    /// Whether `stratum` must fail: `s2:nz0` for the sign bit, every judged one for the ratio.
    fn must_fail(self, stratum: &str) -> bool {
        match self {
            Planted::SignBit => stratum.starts_with("s2:nz0"),
            Planted::FullRatio => true,
        }
    }

    fn answer<S: Real + Into<f64>>(self, fn_id: &str, record: &Record) -> Output {
        let (Some(id), Some(n)) = (Sphere::of_fn(fn_id), s2_of::<S>(record, "n")) else {
            return Output::new();
        };
        match self {
            Planted::SignBit => id.answer_with(&S2Chart::at(&sign_bit(n)), record),
            Planted::FullRatio => {
                // `2δ`, exact at both precisions: the clean answer with §8's ratio as worded.
                let mut out = id.answer_with(&S2Chart::at(&n), record);
                if let Some(d) = out.get_mut("delta") {
                    d.iter_mut().for_each(|x| *x *= 2.0);
                }
                out
            }
        }
    }
}

/// `n` with a `-0` third component moved to `-2⁻¹⁴⁰`, where the comparison reads `ς = -1`: the
/// basis `copysign` would give at `-0`, since `-2⁻¹⁴⁰ - 1 = -1` at both precisions.
fn sign_bit<S: Real>(n: S2<S>) -> S2<S> {
    let [x, y, z] = n.vec().0;
    let negative_zero = z.value_f64().to_bits() == (-0.0_f64).to_bits();
    if negative_zero {
        S2::from_vec_unchecked(Vector([x, y, S::lit(-(2f64.powi(-140)))]))
    } else {
        n
    }
}

struct Plant(Planted);

impl Subject for Plant {
    fn name(&self) -> &str {
        self.0.name()
    }

    fn supports(&self, fn_id: &str) -> bool {
        self.0.ids().contains(&fn_id)
    }

    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        match precision {
            Precision::F64 => self.0.answer::<f64>(fn_id, record),
            Precision::F32 => self.0.answer::<f32>(fn_id, record),
        }
    }
}

/// `(stratum, max_u)` of `subject` over `fn_id` at `precision`.
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

pub(super) fn check(dir: &Path) -> Result<Report, String> {
    let mut text = format!(
        "0063: the sign bit and §8's ratio over s2_*, > {FACTOR} x the clean max where that is \
         below 1/u\n"
    );
    let mut failures = Vec::new();
    for planted in [Planted::SignBit, Planted::FullRatio] {
        for &fn_id in planted.ids() {
            for precision in [Precision::F64, Precision::F32] {
                let cap = (unit_bits(precision) as f64).exp2();
                let clean = maxima(dir, fn_id, Box::new(Helicoid), precision)?;
                let got = maxima(dir, fn_id, Box::new(Plant(planted)), precision)?;
                let mut judged = 0;
                for ((stratum, c), (_, p)) in clean.iter().zip(&got) {
                    let visible = FACTOR * c.max(1.0) < cap;
                    let fail = planted.must_fail(stratum);
                    let _ =
                        writeln!(
                        text,
                        "  {} {fn_id} {stratum} ({precision:?}): clean {c:.3}, planted {p:.3e}{}",
                        planted.name(),
                        if visible || !fail { "" } else { " (at its conditioning, not judged)" }
                    );
                    if !c.is_finite() {
                        failures.push(format!("{fn_id} {stratum}: the clean program reads {c}"));
                    } else if !fail {
                        if p.to_bits() != c.to_bits() {
                            failures.push(format!(
                                "{}: {fn_id} {stratum} moved ({p} against {c}), where n_z is not -0",
                                planted.name()
                            ));
                        }
                    } else if visible {
                        judged += 1;
                        if p.is_nan() || *p <= FACTOR * c.max(1.0) {
                            failures.push(format!(
                                "{}: not detected on {fn_id} {stratum} ({p} against {c})",
                                planted.name()
                            ));
                        }
                    }
                }
                if judged == 0 {
                    failures.push(format!(
                        "{}: {fn_id} has no judged stratum at {precision:?}",
                        planted.name()
                    ));
                }
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
    fn both_s2_readings_are_detected_and_the_clean_program_is_finite() -> Result<(), String> {
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
