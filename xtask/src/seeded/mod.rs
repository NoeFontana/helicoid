//! Seeded subjects (`docs/PHASE1.md` §10): the correct coefficient kernels and the planted
//! defects, run over the `coeff_*` corpus ids as in-process subjects. Only `k, a, b, c, d, e`
//! (`coeff_r` has no seeded kernel), and only the defects that need no sweep: `c` with a switch of
//! `1e-8` and two terms arrives with it. The subject receives `θ`, forms the branch variable
//! `z = fl(θ·θ)` itself (`conformance/generate/README.md`) and reports the value and `d/dz` of one
//! `Dual<f64, 1>` evaluation, so the value path is the plain one (`the_dual_value_path_is_the_plain_value`).

mod kernel;
mod series;

use std::path::Path;

use helicoid_linalg::Precision;

use crate::conformance::corpus::Record;
use crate::conformance::subject::{Output, Registered, Subject};
use kernel::{b_no_series, coefficient, k_sqrt_unsafe, Coeff};
pub(crate) use kernel::{d12, Candidate, D1};
pub(crate) use series::Series;

/// A planted defect (`docs/PHASE1.md` §10), applied to one coefficient; the others stay correct.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Defect {
    /// `b` by its definition, no series.
    BNoSeries,
    /// `k` with an unsafe `sqrt` of `θ²`, under `Dual`.
    KSqrtUnsafe,
}

impl Defect {
    pub(crate) const ALL: [Defect; 2] = [Defect::BNoSeries, Defect::KSqrtUnsafe];

    pub(crate) fn name(self) -> &'static str {
        match self {
            Defect::BNoSeries => "b-no-series",
            Defect::KSqrtUnsafe => "k-sqrt-unsafe",
        }
    }
}

pub(crate) struct Seeded {
    name: String,
    defect: Option<Defect>,
    candidate: Candidate<D1>,
    series: Series<D1>,
}

impl Seeded {
    /// The correct kernels at `candidate`.
    pub(crate) fn correct(series: Series<D1>, candidate: Candidate<D1>) -> Self {
        Self {
            name: "seeded:correct".to_string(),
            defect: None,
            candidate,
            series,
        }
    }

    /// The correct kernels at the `D12` prior, with `defect` planted.
    pub(crate) fn planted(series: Series<D1>, defect: Defect) -> Result<Self, String> {
        let candidate = d12(&series)?;
        Ok(Self {
            name: format!("seeded:{}", defect.name()),
            defect: Some(defect),
            candidate,
            series,
        })
    }

    fn kernel(&self, c: Coeff, z: D1) -> D1 {
        let series = self.series.of(c);
        match (self.defect, c) {
            (Some(Defect::BNoSeries), Coeff::B) => b_no_series(z),
            (Some(Defect::KSqrtUnsafe), Coeff::K) => k_sqrt_unsafe(z, self.candidate, series),
            _ => coefficient(c, z, self.candidate, series),
        }
    }

    pub(crate) fn registered(self) -> Registered {
        Registered {
            version: format!(
                "{}terms-z{}",
                self.candidate.terms, self.candidate.switch_z.v
            ),
            planted: self.defect.is_some(),
            subject: Box::new(self),
        }
    }
}

impl Subject for Seeded {
    fn name(&self) -> &str {
        &self.name
    }

    fn supports(&self, fn_id: &str) -> bool {
        Coeff::of_fn(fn_id).is_some()
    }

    /// `f32` is refused by the harness before any subject runs; an empty answer here is an error
    /// there, never a score.
    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        let (Some(c), Precision::F64) = (Coeff::of_fn(fn_id), precision) else {
            return Output::new();
        };
        let Some(&theta) = record.input("theta").and_then(<[f64]>::first) else {
            return Output::new();
        };
        let r = self.kernel(c, D1::variable(theta * theta, 0));
        Output::from([
            ("value".to_string(), vec![r.v]),
            ("d_branch".to_string(), r.d.to_vec()),
        ])
    }
}

/// The subjects `just conformance` knows: the correct kernels and every defect. The defects are
/// planted: a run that names no subject skips them.
pub(crate) fn registry(corpus: &Path) -> Result<Vec<Registered>, String> {
    let series = Series::load(corpus)?;
    let mut all = vec![Seeded::correct(series.clone(), d12(&series)?).registered()];
    for d in Defect::ALL {
        all.push(Seeded::planted(series.clone(), d)?.registered());
    }
    Ok(all)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;
    use crate::conformance::testkit::record;

    fn subject(defect: Option<Defect>) -> Result<Seeded, String> {
        let series = Series::load(&corpus_dir()?)?;
        match defect {
            None => Ok(Seeded::correct(series.clone(), d12(&series)?)),
            Some(d) => Seeded::planted(series, d),
        }
    }

    /// `d/dz` of the 16-term series, `Σ j s_j z^(j-1)`, in plain `f64`: independent of the kernels'
    /// exact arms and of `Dual`.
    fn series_derivative(s: &[f64], z: f64) -> f64 {
        let step = |acc: f64, (j, &s): (usize, &f64)| j as f64 * s + z * acc;
        s.iter().enumerate().skip(1).rev().fold(0.0, step)
    }

    #[test]
    fn a_subject_answers_the_value_and_the_derivative_in_the_branch_variable() -> Result<(), String>
    {
        let (s, theta) = (subject(None)?, 0.5);
        let rec = record(&[("theta", &[theta])], &[])?;
        let out = s.eval("coeff_a", &rec, Precision::F64);
        assert_eq!(out.keys().collect::<Vec<_>>(), ["d_branch", "value"]);
        // a = (1 - cos(theta)) / theta^2 and d/dz a = (theta sin(theta) - 2 (1 - cos(theta))) / (2 theta^4).
        let (sin, cos) = theta.sin_cos();
        let z = theta * theta;
        assert!((out["value"][0] - (1.0 - cos) / z).abs() < 1e-15);
        assert!(
            (out["d_branch"][0] - (theta * sin - 2.0 * (1.0 - cos)) / (2.0 * z * z)).abs() < 1e-13
        );
        Ok(())
    }

    #[test]
    fn the_derivative_is_in_z_not_in_theta_on_both_sides_of_the_switch() -> Result<(), String> {
        // At theta = 0.5, 2 theta = 1 and d/dz = d/dtheta: any other theta tells them apart.
        let (s, series) = (subject(None)?, Series::<f64>::load(&corpus_dir()?)?);
        for c in Coeff::ALL {
            let fn_id = format!("coeff_{}", c.name());
            for theta in [0.02, 0.06, 0.3, 0.7, 2.0, 3.0] {
                let rec = record(&[("theta", &[theta])], &[])?;
                let out = s.eval(&fn_id, &rec, Precision::F64);
                let (z, s16) = (theta * theta, series.of(c));
                let want = series_derivative(s16, z);
                let got = out["d_branch"][0];
                assert!(
                    (got - want).abs() <= 1e-7 * want.abs(),
                    "{fn_id} at {theta}: {got} vs {want}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn a_subject_supports_the_six_coefficients_and_answers_nothing_it_cannot() -> Result<(), String>
    {
        let s = subject(None)?;
        for id in [
            "coeff_k", "coeff_a", "coeff_b", "coeff_c", "coeff_d", "coeff_e",
        ] {
            assert!(s.supports(id), "{id}");
        }
        for id in ["coeff_r", "coeff_", "so3_exp", "coeff_series"] {
            assert!(!s.supports(id), "{id}");
        }
        let rec = record(&[("theta", &[0.5])], &[])?;
        assert!(s.eval("coeff_r", &rec, Precision::F64).is_empty());
        assert!(s.eval("coeff_k", &rec, Precision::F32).is_empty());
        assert!(s
            .eval("coeff_k", &record(&[], &[])?, Precision::F64)
            .is_empty());
        Ok(())
    }

    #[test]
    fn the_registry_holds_the_correct_kernel_and_every_defect_planted() -> Result<(), String> {
        let all = registry(&corpus_dir()?)?;
        let names: Vec<(&str, bool)> = all.iter().map(|r| (r.subject.name(), r.planted)).collect();
        let want = [
            ("seeded:correct", false),
            ("seeded:b-no-series", true),
            ("seeded:k-sqrt-unsafe", true),
        ];
        assert_eq!(names, want);
        assert!(all.iter().all(|r| r.version == "4terms-z0.01"));
        Ok(())
    }
}
