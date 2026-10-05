//! The per-record errors of every arm, formed once from a kernel's arms ([`Arms`]: the seeded
//! kernel's or `helicoid::coeffs`' through the hidden `__sweep` feature) and the committed corpus
//! (`docs/PHASE1.md` §6): the exact arm, and the arm of each of `1..=TERMS` series terms. At `f32`
//! (`docs/decisions/0016` item 3) the records are the `@f32` strata's and the kernel runs on
//! `Dual<f32, 1>`; everything else, the objective included, is binary64's.

use std::path::Path;

use helicoid_linalg::{Dual, Precision, Real};

use super::search::{Errors, Sample, TERMS};
use crate::conformance::corpus::{self, Record};
use crate::conformance::metric::{Rule, Score, COEFF_D_BRANCH, COEFF_VALUE};
use crate::conformance::subject::Output;
use crate::seeded::{branch_variable, evaluate, input, Candidate, Coeff, Input, Series, Swept};

/// The two arms of one kernel at a record's arguments on `Dual<S, 1>`, seeded at the branch variable.
pub(super) trait Arms<S: Real> {
    /// The exact arm at `x`, with no safe argument: at `z = 0` it is not finite.
    fn exact(&self, id: Swept, x: Input<Dual<S, 1>>) -> Dual<S, 1>;
    /// The series arm of `terms` terms at `x`.
    fn series(&self, id: Swept, x: Input<Dual<S, 1>>, terms: usize) -> Dual<S, 1>;
}

/// The seeded kernel's arms: its own candidate with a switch that always or never selects the series.
pub(super) struct SeededArms<'a, S: Real>(pub(super) &'a Series<Dual<S, 1>>);

impl<S: Real> SeededArms<'_, S> {
    fn at(&self, id: Swept, x: Input<Dual<S, 1>>, terms: usize, switch: f64) -> Dual<S, 1> {
        let cand = Candidate {
            terms,
            switch_z: Dual::<S, 1>::lit(switch),
        };
        evaluate(id, x, cand, self.0.swept(id))
    }
}

impl<S: Real> Arms<S> for SeededArms<'_, S> {
    fn exact(&self, id: Swept, x: Input<Dual<S, 1>>) -> Dual<S, 1> {
        self.at(id, x, 1, 0.0)
    }

    fn series(&self, id: Swept, x: Input<Dual<S, 1>>, terms: usize) -> Dual<S, 1> {
        self.at(id, x, terms, f64::INFINITY)
    }
}

/// The arms `helicoid::coeffs` ships, through `helicoid::__sweep`.
pub(super) struct HelicoidArms;

impl<S: Real> Arms<S> for HelicoidArms {
    fn exact(&self, id: Swept, x: Input<Dual<S, 1>>) -> Dual<S, 1> {
        use helicoid::__sweep as k;
        match id {
            Swept::Coeff(Coeff::K) => k::exact_k(x.z),
            Swept::Coeff(Coeff::A) => k::exact_a(x.z),
            Swept::Coeff(Coeff::B) => k::exact_b(x.z),
            Swept::Coeff(Coeff::C) => k::exact_c(x.z),
            Swept::Coeff(Coeff::D) => k::exact_d(x.z),
            Swept::Coeff(Coeff::E) => k::exact_e(x.z),
            Swept::CosHalf => k::exact_cos_half(x.z),
            Swept::R => k::exact_r(x.z, x.w),
        }
    }

    fn series(&self, id: Swept, x: Input<Dual<S, 1>>, terms: usize) -> Dual<S, 1> {
        use helicoid::__sweep as k;
        match id {
            Swept::Coeff(Coeff::K) => k::series_k(x.z, terms),
            Swept::Coeff(Coeff::A) => k::series_a(x.z, terms),
            Swept::Coeff(Coeff::B) => k::series_b(x.z, terms),
            Swept::Coeff(Coeff::C) => k::series_c(x.z, terms),
            Swept::Coeff(Coeff::D) => k::series_d(x.z, terms),
            Swept::Coeff(Coeff::E) => k::series_e(x.z, terms),
            Swept::CosHalf => k::series_cos_half(x.z, terms),
            Swept::R => k::series_r(x.z, x.w, terms),
        }
    }
}

/// The value and the `d/dz` of `r` against the reference of `rec`, exact, in units of `u` of the
/// precision `r` is computed at.
fn errors<S: Real + Into<f64>>(rec: &Record, r: Dual<S, 1>) -> Result<Errors, String> {
    let out = Output::from([
        ("value".to_string(), vec![r.v.into()]),
        ("d_branch".to_string(), r.d.map(Into::into).to_vec()),
    ]);
    let of = |rule: &Rule| match rule.score(rec, &out, S::PRECISION)? {
        Score::Finite(u) => Ok::<_, String>(u),
        _ => Ok(f64::INFINITY),
    };
    Ok(Errors {
        value: of(&COEFF_VALUE)?,
        deriv: of(&COEFF_D_BRANCH)?,
    })
}

/// The `(value, d/dz)` of one arm, widened to binary64 and taken as bits: the comparison
/// [`super::search::second`] admits a prefix by. Widening is injective, so these agree exactly when the
/// arm's own bits do.
fn bits_of<S: Real + Into<f64>>(r: Dual<S, 1>) -> (u64, u64) {
    let v: f64 = r.v.into();
    let d: f64 = r.d[0].into();
    (v.to_bits(), d.to_bits())
}

/// Every series arm's bits at every grid point, through the same `Arms` the records go through:
/// the dense, reference-free half of stage 2's feasibility rule. The branch variable is the grid
/// point itself — `w = 1` makes `r`'s `n²/w²` the grid point and its `2/w` factor exact, so a
/// prefix agrees with the whole arm here exactly when it does inside `log_ratio`.
pub(super) fn grid_arms<S: Real + Into<f64>>(
    arms: &impl Arms<S>,
    id: Swept,
    grid: &[f64],
) -> Vec<[(u64, u64); TERMS]> {
    let at = |&z: &f64| {
        let x = Input {
            z: S::lit(z),
            w: S::one(),
        };
        core::array::from_fn(|i| bits_of(arms.series(id, x.seed(), i + 1)))
    };
    grid.iter().map(at).collect()
}

/// The samples of one coefficient and, at the same indices, the stratum and id of their records.
pub(super) struct Measured {
    pub(super) samples: Vec<Sample>,
    pub(super) records: Vec<(String, u64)>,
}

/// One sample per `theta:*` record of `coeff_<id>` at `S`'s precision: the binary64 strata's, or
/// at `f32` the `@f32` ones. That filter is what leaves only records with `w > 0` in the objective
/// of `r`, whose branch variable `n²/w²` is infinite at `w = 0` (an exact-arm point whatever the
/// switch): every `theta:*` record has `w > 0` and `q:w0` has none, which
/// `r_is_swept_over_the_records_with_w_above_zero_in_s_and_the_others_are_not` pins.
pub(super) fn samples<S: Real + Into<f64>>(
    dir: &Path,
    arms: &impl Arms<S>,
    id: Swept,
) -> Result<Measured, String> {
    let fn_id = format!("coeff_{}", id.name());
    let entry = corpus::manifest(dir)?
        .into_iter()
        .find(|e| e.fn_id == fn_id)
        .ok_or_else(|| format!("no corpus file for `{fn_id}`"))?;
    let f32 = S::PRECISION == Precision::F32;
    let (mut out, mut records) = (Vec::new(), Vec::new());
    for rec in corpus::read(dir, &entry)? {
        if !rec.stratum.starts_with("theta:") || rec.is_f32_stratum() != f32 {
            continue;
        }
        // `z = fl(θ·θ)` and `s = fl(n²/w²)`, as `Seeded::eval` forms them: the branch variable the
        // kernel compares, so a candidate selects as the kernel does.
        let x =
            input::<S>(id, &rec).ok_or_else(|| format!("{fn_id} {}: no usable input", rec.id))?;
        let exact = errors(&rec, arms.exact(id, x.seed()))?;
        let mut s = Sample {
            z: branch_variable(id, x).into(),
            exact,
            series: [exact; TERMS],
            bits: [(0, 0); TERMS],
        };
        for terms in 1..=TERMS {
            let arm = arms.series(id, x.seed(), terms);
            s.series[terms - 1] = errors(&rec, arm)?;
            s.bits[terms - 1] = bits_of(arm);
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
    use crate::seeded::{d12, Coeff, Defect, Seeded, D1};
    use crate::shipped::shipped;
    use crate::thresholds::grid::grid;
    use crate::thresholds::search::score;
    use crate::thresholds::{CSV_HELICOID, CSV_SEEDED};

    const BOTH: [Precision; 2] = [Precision::F64, Precision::F32];

    /// The largest record score, through the harness's own subject and the `coeff_<id>` rule (the
    /// larger of the two fields, per record), of `subject`, whatever it runs `id` with, over the
    /// `theta:*` records the sweep measures at `precision`.
    fn whole_subject(id: Swept, subject: &Seeded, precision: Precision) -> Result<f64, String> {
        let dir = corpus_dir()?;
        let fn_id = format!("coeff_{}", id.name());
        let entry = corpus::manifest(&dir)?
            .into_iter()
            .find(|e| e.fn_id == fn_id)
            .ok_or("no entry")?;
        let rule = rule(&fn_id).ok_or("no rule")?;
        let mut worst = 0.0f64;
        for rec in corpus::read(&dir, &entry)?.into_iter().filter(|r| {
            r.stratum.starts_with("theta:") && r.is_f32_stratum() == (precision == Precision::F32)
        }) {
            let out = subject.eval(&fn_id, &rec, precision);
            match rule.score(&rec, &out, precision)? {
                Score::Finite(u) => worst = worst.max(u),
                other => return Err(format!("{fn_id} {}: {other:?}", rec.id)),
            }
        }
        Ok(worst)
    }

    /// Every id at `S`'s precision: the objective a candidate scores spliced from the arms formed
    /// once is the one the whole kernel scores, bit for bit, and `subject` runs one candidate at
    /// every id.
    fn spliced_is_whole<S: Real + Into<f64>>(
        dir: &Path,
        subject: impl Fn(&Series<Dual<S, 1>>, Candidate<Dual<S, 1>>) -> Seeded,
    ) -> Result<(), String> {
        let (precision, g) = (S::PRECISION, grid(S::PRECISION));
        let series = Series::<Dual<S, 1>>::load(dir)?;
        let prior = d12(&series)?;
        // The prior, a switch inside the series-only decades, one at the top, and one in the middle:
        // the arms a candidate splices are the arms the whole kernel runs. Then three switches at
        // a record's own branch variable, where `z < switch` is false by one bit: the record's arm
        // is the one the strict comparison of the kernel picks, at the precision it is formed at.
        let fixed = [
            (prior.terms, prior.switch_z.v.into()),
            (1, g[500]),
            (5, g[900]),
            (8, g[1024]),
        ];
        for id in Swept::ALL {
            let Measured { samples, records } = samples(dir, &SeededArms(&series), id)?;
            assert_eq!(
                (samples.len(), records.len()),
                (1710, 1710),
                "{id:?} {precision:?}"
            );
            // Records above `1e-4` only: below it the exact arm is not finite, and a switch that
            // leaves one on it scores infinite spliced but is an error whole.
            let inner: Vec<f64> = samples
                .iter()
                .map(|s| s.z)
                .filter(|z| (1e-4..1.0).contains(z))
                .collect();
            assert!(inner.len() > 100, "{id:?} {precision:?}: {}", inner.len());
            let edges =
                [inner.len() / 4, inner.len() / 2, inner.len() * 3 / 4].map(|i| (3, inner[i]));
            for (terms, switch) in fixed.into_iter().chain(edges) {
                let cand = Candidate {
                    terms,
                    switch_z: Dual::<S, 1>::lit(switch),
                };
                let whole = whole_subject(id, &subject(&series, cand), precision)?;
                assert_eq!(
                    score(&samples, terms, switch).objective().to_bits(),
                    whole.to_bits(),
                    "{id:?} {precision:?} {terms} {switch:e}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn a_candidate_scores_what_the_whole_kernel_scores() -> Result<(), String> {
        let dir = corpus_dir()?;
        spliced_is_whole::<f64>(&dir, Seeded::uniform)?;
        spliced_is_whole::<f32>(&dir, Seeded::uniform_f32)
    }

    /// `(terms, switch_bits, value_max_u, deriv_max_u)` of the row of `id` at `at` in the committed
    /// sweep `csv`.
    fn committed_row(
        csv: &str,
        id: Swept,
        at: Precision,
    ) -> Result<(usize, u64, f64, f64), String> {
        let path = crate::conformance::root()?.join(csv);
        let text = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
        let header: Vec<&str> = super::super::HEADER.split(',').collect();
        let precision = format!("{at:?}").to_lowercase();
        let row = text
            .lines()
            .find(|l| l.starts_with(&format!("{},{precision},", id.name())))
            .ok_or("no row")?;
        let cells: Vec<&str> = row.split(',').collect();
        let cell = |name: &str| {
            let i = header.iter().position(|h| *h == name).ok_or(name)?;
            Ok::<_, String>(cells[i])
        };
        let float = |name: &str| cell(name)?.parse::<f64>().map_err(|e| e.to_string());
        let bits = u64::from_str_radix(cell("switch_bits")?.trim_start_matches("0x"), 16)
            .map_err(|e| e.to_string())?;
        let terms = cell("terms")?.parse().map_err(|e| format!("{e}"))?;
        Ok((terms, bits, float("value_max_u")?, float("deriv_max_u")?))
    }

    #[test]
    fn the_generated_kernel_scores_the_objective_the_sweep_chose() -> Result<(), String> {
        // The subject that runs `generated.rs` is the candidate the CSV records: the same terms and
        // switch bits, the same series (the largest score is the sweep's objective to the bit),
        // whatever order the constants are listed in, at both precisions.
        let generated = Seeded::generated();
        for at in BOTH {
            for id in Swept::ALL {
                let (terms, bits, value, deriv) = committed_row(CSV_SEEDED, id, at)?;
                let (n, switch, series) = generated.switch(id, at).ok_or("no kernel")?;
                let bits = match at {
                    Precision::F64 => f64::from_bits(bits),
                    Precision::F32 => f64::from(f32::from_bits(bits as u32)),
                };
                assert_eq!((n, series.len()), (terms, terms), "{id:?} {at:?}");
                assert_eq!(switch.to_bits(), bits.to_bits(), "{id:?} {at:?}");
                let worst = whole_subject(id, &generated, at)?;
                assert_eq!(worst.to_bits(), value.max(deriv).to_bits(), "{id:?} {at:?}");
            }
        }
        Ok(())
    }

    #[test]
    fn the_generated_series_is_the_corpus_series_rounded_once_at_its_precision(
    ) -> Result<(), String> {
        let (generated, dir) = (Seeded::generated(), corpus_dir()?);
        let (wide, narrow) = (Series::<f64>::load(&dir)?, Series::<f32>::load(&dir)?);
        for id in Swept::ALL {
            let want64: Vec<f64> = wide.swept(id).to_vec();
            let want32: Vec<f64> = narrow.swept(id).iter().map(|&x| f64::from(x)).collect();
            for (at, want) in [(Precision::F64, want64), (Precision::F32, want32)] {
                let got = generated.switch(id, at).ok_or("no kernel")?.2;
                let (got, want): (Vec<u64>, Vec<u64>) = (
                    got.iter().map(|x| x.to_bits()).collect(),
                    want.iter().map(|x| x.to_bits()).collect(),
                );
                assert_eq!(got, want[..got.len()], "{id:?} {at:?}");
            }
        }
        for c in Coeff::ALL {
            assert!(generated.series(c).iter().all(|x| x.d[0] == 0.0), "{c:?}");
        }
        Ok(())
    }

    #[test]
    fn the_planted_c_scores_what_the_sweep_ranks() -> Result<(), String> {
        // The subject that plants `c` scores, over the whole corpus, exactly the objective the
        // sweep's `score` gives its candidate: the ranking is of the subject that runs.
        let dir = corpus_dir()?;
        let series = Series::<D1>::load(&dir)?;
        let c = Swept::Coeff(Coeff::C);
        let samples = samples(&dir, &SeededArms(&series), c)?.samples;
        let (terms, switch) = crate::seeded::C_PLANTED;
        let planted = Seeded::planted(Defect::CTwoTermsEarly);
        assert_eq!(planted.candidate(Coeff::C).terms, terms);
        assert_eq!(
            planted.candidate(Coeff::C).switch_z.v.to_bits(),
            switch.to_bits()
        );
        let (whole, spliced) = (
            whole_subject(c, &planted, Precision::F64)?,
            score(&samples, terms, switch).objective(),
        );
        assert_eq!(whole.to_bits(), spliced.to_bits());
        // Every other coefficient of the planted subject is the generated one.
        let generated = Seeded::generated();
        for c in Coeff::ALL
            .into_iter()
            .filter(|&c| c != Coeff::C && c != Coeff::B && c != Coeff::K)
        {
            let id = Swept::Coeff(c);
            assert_eq!(
                whole_subject(id, &planted, Precision::F64)?.to_bits(),
                whole_subject(id, &generated, Precision::F64)?.to_bits(),
                "{c:?}"
            );
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
        let Measured { samples, records } =
            samples(&dir, &SeededArms(&series), Swept::Coeff(Coeff::B))?;
        let zero = samples.iter().position(|s| s.z == 0.0).ok_or("no z = 0")?;
        assert_eq!(records[zero].0, "theta:exact0");
        assert!(samples[zero].exact.value.is_infinite());
        assert!(samples[zero].series.iter().all(|e| e.value.is_finite()));
        Ok(())
    }

    #[test]
    fn r_is_swept_over_the_records_with_w_above_zero_in_s_and_the_others_are_not(
    ) -> Result<(), String> {
        // `q:w0` is `(n, +0)`: `s` is infinite, the exact arm whatever the switch, so it is not in
        // an objective that selects an arm. Its records are in the corpus, three each at `f64` and
        // `f32`, and none of the 1710 `theta:*` records has `w <= 0`.
        let dir = corpus_dir()?;
        let fn_id = "coeff_r";
        let entry = corpus::manifest(&dir)?
            .into_iter()
            .find(|e| e.fn_id == fn_id);
        let records = corpus::read(&dir, &entry.ok_or("no entry")?)?;
        let w0 = |r: &&Record| r.input("w").is_some_and(|w| w[0] <= 0.0);
        assert!(records
            .iter()
            .filter(w0)
            .all(|r| r.stratum.starts_with("q:w0")));
        assert_eq!(records.iter().filter(w0).count(), 6);
        let series = Series::<D1>::load(&dir)?;
        let m = samples(&dir, &SeededArms(&series), Swept::R)?;
        assert!(m.records.iter().all(|(s, _)| s.starts_with("theta:")));
        // `s = n²/w²` is `tan²(θ/2)`: from 0 at `theta:exact0` to `w ≈ 5e-13` at `theta:pi-1e-12`.
        let zs = m.samples.iter().map(|s| s.z);
        assert!(zs.clone().fold(f64::INFINITY, f64::min) == 0.0);
        assert!(zs.fold(0.0, f64::max) > 1e24);
        Ok(())
    }

    /// Equal to the bit, but for the payload of a NaN, which Rust leaves unspecified.
    fn same(a: f64, b: f64) -> bool {
        a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan())
    }

    /// The records of `coeff_<id>` at `S`'s precision: the strata without the suffix at binary64,
    /// the `@f32` ones at binary32.
    fn records_at<S: Real>(dir: &Path, id: Swept) -> Result<Vec<Record>, String> {
        let fn_id = format!("coeff_{}", id.name());
        let entry = corpus::manifest(dir)?
            .into_iter()
            .find(|e| e.fn_id == fn_id);
        let f32 = S::PRECISION == Precision::F32;
        let all = corpus::read(dir, &entry.ok_or("no entry")?)?;
        Ok(all
            .into_iter()
            .filter(|r| r.is_f32_stratum() == f32)
            .collect())
    }

    /// The shipped arms are the seeded arms bit for bit, value and derivative, exact and every
    /// series length, at every record of every stratum of every id, and there are that many.
    fn the_shipped_arms_are_the_seeded_arms<S: Real + Into<f64>>() -> Result<usize, String> {
        let dir = corpus_dir()?;
        let series = Series::<Dual<S, 1>>::load(&dir)?;
        let (seeded, mut n) = (SeededArms(&series), 0);
        for id in Swept::ALL {
            for rec in records_at::<S>(&dir, id)? {
                let x = input::<S>(id, &rec).ok_or("no usable input")?.seed();
                // `r` at `w = +0` (`q:w0`) has no series arm: the seeded kernel's mask sends every
                // candidate to the exact arm there, and the shipped series arm is `2/w` times one.
                let terms = if x.w.v.into() > 0.0 || id != Swept::R {
                    TERMS
                } else {
                    0
                };
                let arms = std::iter::once((seeded.exact(id, x), HelicoidArms.exact(id, x))).chain(
                    (1..=terms).map(|t| (seeded.series(id, x, t), HelicoidArms.series(id, x, t))),
                );
                for (arm, (a, b)) in arms.enumerate() {
                    let (v, d) = (
                        same(a.v.into(), b.v.into()),
                        same(a.d[0].into(), b.d[0].into()),
                    );
                    assert!(v && d, "{id:?} {} {} arm {arm}", rec.stratum, rec.id);
                    n += 1;
                }
            }
        }
        Ok(n)
    }

    #[test]
    fn the_shipped_arms_are_the_seeded_arms_bit_for_bit() -> Result<(), String> {
        // Eight ids, 3420 records per id split between the two precisions (`coeff_r`: 3426), and
        // `TERMS` series lengths beside the exact arm.
        let (wide, narrow) = (
            the_shipped_arms_are_the_seeded_arms::<f64>()?,
            the_shipped_arms_are_the_seeded_arms::<f32>()?,
        );
        // `TERMS + 1` arms of each record but the six `q:w0` records of `r`, which have the
        // exact one. The `+ 1` is the exact arm; the count scales with the cap `0039` lifted.
        assert_eq!(wide + narrow, (7 * 3420 + 3426) * (TERMS + 1) - 6 * TERMS);
        Ok(())
    }

    /// The shipped groups over the corpus: every record of every stratum finite in value and
    /// derivative (the safe argument at `θ = 0`, subnormal, `π`, `w = +0`), and over the `theta:*`
    /// strata the objective the sweep recorded for the switch it chose, to the bit.
    fn the_shipped_groups_score_the_objective_the_sweep_chose<S: Real + Into<f64>>(
    ) -> Result<(), String> {
        let dir = corpus_dir()?;
        let at = S::PRECISION;
        for id in Swept::ALL {
            let (mut value, mut deriv) = (0.0f64, 0.0f64);
            for rec in records_at::<S>(&dir, id)? {
                let x = input::<S>(id, &rec).ok_or("no usable input")?.seed();
                let e = errors(&rec, shipped(id, x))?;
                assert!(
                    e.value.is_finite() && e.deriv.is_finite(),
                    "{id:?} {}",
                    rec.stratum
                );
                if rec.stratum.starts_with("theta:") {
                    (value, deriv) = (value.max(e.value), deriv.max(e.deriv));
                }
            }
            let (_, _, want_value, want_deriv) = committed_row(CSV_HELICOID, id, at)?;
            let name = id.name();
            assert_eq!(
                value.to_bits(),
                want_value.to_bits(),
                "{name} {at:?}: value"
            );
            assert_eq!(
                deriv.to_bits(),
                want_deriv.to_bits(),
                "{name} {at:?}: derivative"
            );
        }
        Ok(())
    }

    #[test]
    fn the_shipped_groups_score_the_objective_the_sweep_chose_at_f64() -> Result<(), String> {
        the_shipped_groups_score_the_objective_the_sweep_chose::<f64>()
    }

    #[test]
    fn the_shipped_groups_score_the_objective_the_sweep_chose_at_f32() -> Result<(), String> {
        the_shipped_groups_score_the_objective_the_sweep_chose::<f32>()
    }

    fn b_is_the_same_in_both_groups<S: Real + Into<f64>>() -> Result<(), String> {
        use helicoid::__sweep as k;
        let dir = corpus_dir()?;
        let id = Swept::Coeff(Coeff::B);
        let records = records_at::<S>(&dir, id)?;
        assert!(!records.is_empty());
        for rec in records {
            let z = input::<S>(id, &rec).ok_or("no input")?.seed().z;
            let (jr, q) = (k::jr_coeffs(z).1, k::q_coeffs(z).0);
            assert!(
                same(jr.v.into(), q.v.into()) && same(jr.d[0].into(), q.d[0].into()),
                "{} {}",
                rec.stratum,
                rec.id
            );
        }
        Ok(())
    }

    #[test]
    fn b_is_the_same_in_both_groups_that_hold_it() -> Result<(), String> {
        b_is_the_same_in_both_groups::<f64>()
    }

    #[test]
    fn b_is_the_same_in_both_groups_that_hold_it_at_f32() -> Result<(), String> {
        b_is_the_same_in_both_groups::<f32>()
    }
}
