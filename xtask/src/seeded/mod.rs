//! Seeded subjects (`docs/PHASE1.md` §10): the correct coefficient kernels and the planted
//! defects, run over the `coeff_*` corpus ids as in-process subjects. Only `k, a, b, c, d, e`
//! (`coeff_r` has no seeded kernel).
//!
//! The correct kernel runs the **generated** switches of each coefficient: its series length,
//! switch and series terms are `generated.rs`'s, which `cargo xtask thresholds` writes
//! (`crate::thresholds`), and nothing here is typed. The sweep measures the same kernel through
//! [`coefficient`] with the candidate it scores, so the constants the subject runs are the
//! candidate the sweep chose (`the_generated_kernel_scores_the_objective_the_sweep_chose`).
//!
//! The planted defects are the correct kernels with one coefficient changed. `c` with a switch of
//! `1e-8` and two terms ([`C_PLANTED`]) is named by the sweep's ranking (`conformance::selftest`);
//! its other mechanism, the envelope, is not started. The subject receives `θ`, forms the branch
//! variable `z = fl(θ·θ)` itself (`conformance/generate/README.md`) and reports the value and
//! `d/dz` of one `Dual<f64, 1>` evaluation, so the value path is the plain one
//! (`the_dual_value_path_is_the_plain_value`).

mod generated;
mod kernel;
mod series;
mod switch;

use helicoid_linalg::Precision;

use crate::conformance::corpus::Record;
use crate::conformance::subject::{Output, Registered, Subject};
use generated::{A_F64, B_F64, C_F64, D_F64, E_F64, K_F64};
use kernel::{b_no_series, k_sqrt_unsafe};
pub(crate) use kernel::{coefficient, d12, Candidate, Coeff, D1};
pub(crate) use series::{Series, FILE as SERIES_FILE};

/// The planted `c` (`docs/PHASE1.md` §10): two terms below `z = 1e-16`, §10's `1e-8` read as `θ`
/// (0014 (draft) question 14), which is the first point of the sweep's grid.
pub(crate) const C_PLANTED: (usize, f64) = (2, 1e-16);

/// A planted defect (`docs/PHASE1.md` §10), applied to one coefficient; the others stay correct.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Defect {
    /// `b` by its definition, no series.
    BNoSeries,
    /// `k` with an unsafe `sqrt` of `θ²`, under `Dual`.
    KSqrtUnsafe,
    /// `c` with switch `1e-8` and two terms ([`C_PLANTED`]).
    CTwoTermsEarly,
}

impl Defect {
    pub(crate) const ALL: [Defect; 3] = [
        Defect::BNoSeries,
        Defect::KSqrtUnsafe,
        Defect::CTwoTermsEarly,
    ];

    pub(crate) fn name(self) -> &'static str {
        match self {
            Defect::BNoSeries => "b-no-series",
            Defect::KSqrtUnsafe => "k-sqrt-unsafe",
            Defect::CTwoTermsEarly => "c-two-terms-1e-8",
        }
    }
}

/// One coefficient's kernel: how many terms below which switch, and the series they are the first
/// terms of.
struct Arm {
    candidate: Candidate<D1>,
    series: Vec<D1>,
}

impl Arm {
    /// The generated switch `(below, series)`: every term of the series, below `below`.
    fn generated((below, series): (f64, &[f64])) -> Self {
        Self {
            candidate: Candidate {
                terms: series.len(),
                switch_z: D1::constant(below),
            },
            series: series.iter().map(|&x| D1::constant(x)).collect(),
        }
    }
}

pub(crate) struct Seeded {
    name: String,
    version: String,
    defect: Option<Defect>,
    /// By [`Coeff::index`].
    arms: [Arm; 6],
}

impl Seeded {
    /// The correct kernels at the generated switches, `k, a, b, c, d, e`.
    pub(crate) fn generated() -> Self {
        let switches = [
            K_F64.parts(),
            A_F64.parts(),
            B_F64.parts(),
            C_F64.parts(),
            D_F64.parts(),
            E_F64.parts(),
        ];
        Self {
            name: "seeded:correct".to_string(),
            version: "generated".to_string(),
            defect: None,
            arms: switches.map(Arm::generated),
        }
    }

    /// The correct kernels, all six at one `candidate` over the corpus's series: how a test drives
    /// the harness at a candidate the sweep also scores, the `D12` prior included.
    #[cfg(test)]
    pub(crate) fn uniform(series: &Series<D1>, candidate: Candidate<D1>) -> Self {
        Self {
            name: "seeded:correct".to_string(),
            version: format!("{}terms-z{}", candidate.terms, candidate.switch_z.v),
            defect: None,
            arms: Coeff::ALL.map(|c| Arm {
                candidate,
                series: series.of(c).to_vec(),
            }),
        }
    }

    /// The generated kernels with `defect` planted.
    pub(crate) fn planted(defect: Defect) -> Self {
        let mut seeded = Self::generated();
        seeded.name = format!("seeded:{}", defect.name());
        seeded.defect = Some(defect);
        if defect == Defect::CTwoTermsEarly {
            let (terms, z) = C_PLANTED;
            seeded.arms[Coeff::C.index()].candidate = Candidate {
                terms,
                switch_z: D1::constant(z),
            };
        }
        seeded
    }

    /// The series length and switch this subject runs `c` with.
    pub(crate) fn candidate(&self, c: Coeff) -> Candidate<D1> {
        self.arms[c.index()].candidate
    }

    /// The series this subject's `c` takes its terms from.
    #[cfg(test)]
    pub(crate) fn series(&self, c: Coeff) -> &[D1] {
        &self.arms[c.index()].series
    }

    fn kernel(&self, c: Coeff, z: D1) -> D1 {
        let Arm { candidate, series } = &self.arms[c.index()];
        match (self.defect, c) {
            (Some(Defect::BNoSeries), Coeff::B) => b_no_series(z),
            (Some(Defect::KSqrtUnsafe), Coeff::K) => k_sqrt_unsafe(z, *candidate, series),
            _ => coefficient(c, z, *candidate, series),
        }
    }

    pub(crate) fn registered(self) -> Registered {
        Registered {
            version: self.version.clone(),
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
pub(crate) fn registry() -> Vec<Registered> {
    let mut all = vec![Seeded::generated().registered()];
    all.extend(Defect::ALL.map(|d| Seeded::planted(d).registered()));
    all
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;
    use crate::conformance::testkit::record;

    fn subject(defect: Option<Defect>) -> Seeded {
        match defect {
            None => Seeded::generated(),
            Some(d) => Seeded::planted(d),
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
        let (s, theta) = (subject(None), 0.5);
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
        let (s, series) = (subject(None), Series::<f64>::load(&corpus_dir()?)?);
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
        let s = subject(None);
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
    fn the_registry_holds_the_correct_kernel_and_every_defect_planted() {
        let all = registry();
        let names: Vec<(&str, bool)> = all.iter().map(|r| (r.subject.name(), r.planted)).collect();
        let want = [
            ("seeded:correct", false),
            ("seeded:b-no-series", true),
            ("seeded:k-sqrt-unsafe", true),
            ("seeded:c-two-terms-1e-8", true),
        ];
        assert_eq!(names, want);
        assert!(all.iter().all(|r| r.version == "generated"));
    }

    #[test]
    fn a_defect_changes_one_coefficient_of_the_generated_kernel_and_no_other() {
        let generated = Seeded::generated();
        let planted = Seeded::planted(Defect::CTwoTermsEarly);
        for c in Coeff::ALL {
            let (a, b) = (generated.candidate(c), planted.candidate(c));
            let same = (a.terms, a.switch_z.v.to_bits()) == (b.terms, b.switch_z.v.to_bits());
            assert_eq!(same, c != Coeff::C, "{c:?}");
        }
        // Two terms below 1e-16: §10's 1e-8 as θ, and the grid's first point.
        let c = planted.candidate(Coeff::C);
        assert_eq!((c.terms, c.switch_z.v.to_bits()), (2, 1e-16f64.to_bits()));
        assert_eq!(C_PLANTED, (2, 1e-16));
        // Its series is the generated one, of which two terms are used.
        assert_eq!(planted.arms[Coeff::C.index()].series.len(), 8);
    }
}
