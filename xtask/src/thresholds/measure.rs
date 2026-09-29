//! The per-record errors of every arm, formed once from the seeded kernel and the committed
//! corpus (`docs/PHASE1.md` §6): the exact arm, and the arm of each of `1..=TERMS` series terms.

use std::path::Path;

use helicoid_linalg::Precision;

use super::search::{Errors, Sample, TERMS};
use crate::conformance::corpus::{self, Record};
use crate::conformance::metric::{Rule, Score, COEFF_D_BRANCH, COEFF_VALUE};
use crate::conformance::subject::Output;
use crate::seeded::{coefficient, Candidate, Coeff, Series, D1};

/// The kernel of `c` at `z` with `terms` terms below `switch`: `0` selects the exact arm at `z`
/// and `+inf` the series arm at every finite `z`. Both go through the kernel's own branch, so
/// each arm is what a candidate that selects it evaluates.
fn kernel(c: Coeff, series: &Series<D1>, z: f64, terms: usize, switch: f64) -> D1 {
    let candidate = Candidate {
        terms,
        switch_z: D1::constant(switch),
    };
    coefficient(c, D1::variable(z, 0), candidate, series.of(c))
}

/// The value and the `d/dz` of `r` against the reference of `rec`, exact, in units of `u`.
fn errors(rec: &Record, r: D1) -> Result<Errors, String> {
    let out = Output::from([
        ("value".to_string(), vec![r.v]),
        ("d_branch".to_string(), r.d.to_vec()),
    ]);
    let of = |rule: &Rule| match rule.score(rec, &out, Precision::F64)? {
        Score::Finite(u) => Ok::<_, String>(u),
        _ => Ok(f64::INFINITY),
    };
    Ok(Errors {
        value: of(&COEFF_VALUE)?,
        deriv: of(&COEFF_D_BRANCH)?,
    })
}

/// The samples of one coefficient and, at the same indices, the stratum and id of their records.
pub(super) struct Measured {
    pub(super) samples: Vec<Sample>,
    pub(super) records: Vec<(String, u64)>,
}

/// One sample per `theta:*` record of `coeff_<c>` (every record of it, today).
pub(super) fn samples(dir: &Path, series: &Series<D1>, c: Coeff) -> Result<Measured, String> {
    let fn_id = format!("coeff_{}", c.name());
    let entry = corpus::manifest(dir)?
        .into_iter()
        .find(|e| e.fn_id == fn_id)
        .ok_or_else(|| format!("no corpus file for `{fn_id}`"))?;
    let (mut out, mut records) = (Vec::new(), Vec::new());
    for rec in corpus::read(dir, &entry)? {
        if !rec.stratum.starts_with("theta:") {
            continue;
        }
        let theta = rec.input("theta").and_then(<[f64]>::first);
        // `z = fl(θ·θ)`, as `Seeded::eval` forms it: the same branch variable, the same selection.
        let z = theta
            .map(|t| t * t)
            .ok_or_else(|| format!("{fn_id} {}: no theta", rec.id))?;
        let arm = |terms, switch| errors(&rec, kernel(c, series, z, terms, switch));
        let exact = arm(1, 0.0)?;
        let mut s = Sample {
            z,
            exact,
            series: [exact; TERMS],
        };
        for terms in 1..=TERMS {
            s.series[terms - 1] = arm(terms, f64::INFINITY)?;
        }
        out.push(s);
        records.push((rec.stratum, rec.id));
    }
    match out.is_empty() {
        true => Err(format!("`{fn_id}` has no `theta:*` record")),
        false => Ok(Measured {
            samples: out,
            records,
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;
    use crate::conformance::metric::rule;
    use crate::conformance::subject::Subject;
    use crate::seeded::{d12, Seeded};
    use crate::thresholds::grid::grid;
    use crate::thresholds::search::score;

    /// The largest record score of the whole seeded kernel at `candidate`, through the harness's
    /// own subject and the `coeff_<c>` rule (the larger of the two fields, per record).
    fn whole_kernel(
        c: Coeff,
        series: &Series<D1>,
        candidate: Candidate<D1>,
    ) -> Result<f64, String> {
        let dir = corpus_dir()?;
        let fn_id = format!("coeff_{}", c.name());
        let entry = corpus::manifest(&dir)?
            .into_iter()
            .find(|e| e.fn_id == fn_id)
            .ok_or("no entry")?;
        let subject = Seeded::correct(series.clone(), candidate);
        let rule = rule(&fn_id).ok_or("no rule")?;
        let mut worst = 0.0f64;
        for rec in corpus::read(&dir, &entry)? {
            let out = subject.eval(&fn_id, &rec, Precision::F64);
            match rule.score(&rec, &out, Precision::F64)? {
                Score::Finite(u) => worst = worst.max(u),
                other => return Err(format!("{fn_id} {}: {other:?}", rec.id)),
            }
        }
        Ok(worst)
    }

    #[test]
    fn a_candidate_scores_what_the_whole_kernel_scores() -> Result<(), String> {
        let (dir, g) = (corpus_dir()?, grid());
        let series = Series::<D1>::load(&dir)?;
        let prior = d12(&series)?;
        // The prior, a switch inside the series-only decades, one at the top, and one in the middle:
        // the arms a candidate splices are the arms the whole kernel runs, bit for bit.
        let candidates = [
            (prior.terms, prior.switch_z.v),
            (1, g[500]),
            (5, g[900]),
            (8, g[1024]),
        ];
        for c in Coeff::ALL {
            let Measured { samples, records } = samples(&dir, &series, c)?;
            assert_eq!((samples.len(), records.len()), (1710, 1710), "{c:?}");
            for (terms, switch) in candidates {
                let cand = Candidate {
                    terms,
                    switch_z: D1::constant(switch),
                };
                let (spliced, whole) = (
                    score(&samples, terms, switch),
                    whole_kernel(c, &series, cand)?,
                );
                assert_eq!(
                    spliced.objective().to_bits(),
                    whole.to_bits(),
                    "{c:?} {terms} {switch:e}"
                );
            }
        }
        Ok(())
    }

    /// The first record of `coeff_<c>` in `stratum`.
    fn record(dir: &Path, c: Coeff, stratum: &str) -> Result<Record, String> {
        let fn_id = format!("coeff_{}", c.name());
        let entry = corpus::manifest(dir)?
            .into_iter()
            .find(|e| e.fn_id == fn_id)
            .ok_or("no entry")?;
        let mut all = corpus::read(dir, &entry)?.into_iter();
        all.find(|r| r.stratum == stratum)
            .ok_or_else(|| format!("{fn_id}: no `{stratum}` record"))
    }

    #[test]
    fn an_output_that_is_not_finite_is_an_infinite_error_in_its_own_field() -> Result<(), String> {
        let rec = record(&corpus_dir()?, Coeff::K, "theta:dense")?;
        let (nan, inf) = (f64::NAN, f64::INFINITY);
        let finite = errors(&rec, D1 { v: 0.5, d: [-0.1] })?;
        assert!(finite.value.is_finite() && finite.deriv.is_finite());
        // Neither zero nor NaN: a non-finite arm must lose every comparison it enters.
        for (v, d) in [(nan, -0.1), (inf, -0.1), (0.5, nan), (0.5, -inf)] {
            let e = errors(&rec, D1 { v, d: [d] })?;
            assert_eq!(e.value.is_infinite(), !v.is_finite(), "{v} {d}");
            assert_eq!(e.deriv.is_infinite(), !d.is_finite(), "{v} {d}");
        }
        Ok(())
    }

    #[test]
    fn the_exact_arm_of_b_at_zero_is_infinite_so_no_candidate_selects_it() -> Result<(), String> {
        // `b` is `0/0` at `z = 0`: the grid's smallest switch is 1e-16, so the series arm is
        // taken there, and only an infinite score keeps the exact arm from ever being chosen.
        let dir = corpus_dir()?;
        let series = Series::<D1>::load(&dir)?;
        let Measured { samples, records } = samples(&dir, &series, Coeff::B)?;
        let zero = samples.iter().position(|s| s.z == 0.0).ok_or("no z = 0")?;
        assert_eq!(records[zero].0, "theta:exact0");
        assert!(samples[zero].exact.value.is_infinite());
        assert!(samples[zero].series.iter().all(|e| e.value.is_finite()));
        Ok(())
    }
}
