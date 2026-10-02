//! The `helicoid` conformance subject (`docs/PHASE1.md` §5) over the eight `coeff_*` ids: the
//! shipped groups of `helicoid::coeffs` through the hidden `__sweep` feature, so the code scored is
//! the code shipped (`0004` item 4), at both precisions on `Dual<S, 1>` seeded at the branch variable.
//!
//! The arguments are the seeded subject's (`seeded::input`: `z = fl(θ·θ)`, `n²` and `w` for `r`,
//! formed at the precision, `@f32` inputs cast losslessly) and the answer its fields, `value` and
//! `d_branch`. Each id reads the group that holds it at the call sites of `docs/PHASE3.md` §3:
//! `k` and `cos θ/2` from `exp_coeffs`, `a` and `b` from `jr_coeffs`, `c` from `jr_inv_coeff`, `d`
//! and `e` from `q_coeffs`, `r` from `log_ratio`. `b` also lives in `q_coeffs`; the two agree to
//! the bit at both precisions (`b_is_the_same_in_both_groups_that_hold_it` and its `_at_f32` twin,
//! sweep tests), so one is scored. No other id is supported: `helicoid` has no group code yet.
//!
//! What these tests do not cover: the boundary (`s = switch`, `w = 0`) and the safe-argument
//! pattern. A scalar mask makes `S::branch` lazy, so the exact arm never sees an unsafe argument
//! and a finite score is no evidence for it; the `helicoid` crate's lane and boundary tests are.
//!
//! `subject_version` is the workspace version: there is no other version of the shipped kernel.

use helicoid::__sweep as k;
use helicoid_linalg::{Dual, Precision, Real};

use crate::conformance::corpus::Record;
use crate::conformance::subject::{Output, Registered, Subject};
use crate::seeded::{input, Coeff, Input, Swept};

/// The group of `id` at `x`, as the shipped kernel evaluates it.
pub(crate) fn shipped<S: Real>(id: Swept, x: Input<Dual<S, 1>>) -> Dual<S, 1> {
    match id {
        Swept::Coeff(Coeff::K) => k::exp_coeffs(x.z).0,
        Swept::CosHalf => k::exp_coeffs(x.z).1,
        Swept::Coeff(Coeff::A) => k::jr_coeffs(x.z).0,
        Swept::Coeff(Coeff::B) => k::jr_coeffs(x.z).1,
        Swept::Coeff(Coeff::C) => k::jr_inv_coeff(x.z),
        Swept::Coeff(Coeff::D) => k::q_coeffs(x.z).1,
        Swept::Coeff(Coeff::E) => k::q_coeffs(x.z).2,
        Swept::R => k::log_ratio(x.z, x.w),
    }
}

/// The subject's answer at `S`: the value and `d/dz`, widened exactly to binary64; nothing when the
/// record holds no usable input.
fn answer<S: Real + Into<f64>>(id: Swept, record: &Record) -> Output {
    let Some(x) = input::<S>(id, record) else {
        return Output::new();
    };
    let r = shipped(id, x.seed());
    Output::from([
        ("value".to_string(), vec![r.v.into()]),
        ("d_branch".to_string(), r.d.map(Into::into).to_vec()),
    ])
}

pub(crate) struct Helicoid;

impl Subject for Helicoid {
    fn name(&self) -> &str {
        "helicoid"
    }

    fn supports(&self, fn_id: &str) -> bool {
        Swept::of_fn(fn_id).is_some()
    }

    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        let Some(id) = Swept::of_fn(fn_id) else {
            return Output::new();
        };
        match precision {
            Precision::F64 => answer::<f64>(id, record),
            Precision::F32 => answer::<f32>(id, record),
        }
    }
}

pub(crate) fn registered() -> Registered {
    let version = env!("CARGO_PKG_VERSION").to_string();
    Registered {
        version: version.clone(),
        version_f32: version,
        no_f32: None,
        subject: Box::new(Helicoid),
        planted: false,
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;
    use crate::conformance::corpus::{self, Record};
    use crate::conformance::corpus_dir;
    use crate::conformance::metric::{rule, Score};
    use crate::conformance::subject::registry;
    use crate::conformance::testkit::record;
    use crate::seeded::{d12, Seeded, Series, D1};

    type F32 = Dual<f32, 1>;

    /// Every record of `coeff_<id>` at `precision`: the strata without the suffix at binary64, the
    /// `@f32` ones at binary32.
    fn records(id: Swept, precision: Precision) -> Result<Vec<Record>, String> {
        let dir = corpus_dir()?;
        let fn_id = format!("coeff_{}", id.name());
        let entry = corpus::manifest(&dir)?
            .into_iter()
            .find(|e| e.fn_id == fn_id);
        let all = corpus::read(&dir, &entry.ok_or("no corpus file")?)?;
        let f32 = precision == Precision::F32;
        Ok(all
            .into_iter()
            .filter(|r| r.is_f32_stratum() == f32)
            .collect())
    }

    /// The largest score per stratum of `subject` on `id`, `None` for a non-finite record.
    fn maxima(
        subject: &dyn Subject,
        id: Swept,
        precision: Precision,
    ) -> Result<BTreeMap<String, Option<f64>>, String> {
        let fn_id = format!("coeff_{}", id.name());
        let rule = rule(&fn_id).ok_or("no rule")?;
        let mut out: BTreeMap<String, Option<f64>> = BTreeMap::new();
        for rec in records(id, precision)? {
            let answer = subject.eval(&fn_id, &rec, precision);
            let score = match rule.score(&rec, &answer, precision)? {
                Score::Finite(u) => Some(u),
                _ => None,
            };
            let slot = out.entry(rec.stratum.clone()).or_insert(Some(0.0));
            *slot = slot.zip(score).map(|(a, b)| a.max(b));
        }
        Ok(out)
    }

    /// Equal to the bit, but for the payload of a NaN, which Rust leaves unspecified.
    fn same(a: &Output, b: &Output) -> bool {
        let bits = |o: &Output| -> Vec<(String, Vec<u64>)> {
            let nan = f64::NAN.to_bits();
            let map = |x: &f64| if x.is_nan() { nan } else { x.to_bits() };
            o.iter()
                .map(|(k, v)| (k.clone(), v.iter().map(map).collect()))
                .collect()
        };
        bits(a) == bits(b)
    }

    #[test]
    fn the_registry_holds_helicoid_as_a_plain_subject_with_both_precisions() {
        let all = registry();
        let h = &all[0];
        assert_eq!((h.subject.name(), h.planted), ("helicoid", false));
        assert_eq!(h.version, env!("CARGO_PKG_VERSION"));
        assert_eq!(h.version_at(Precision::F32), h.version_at(Precision::F64));
        assert!(h.no_f32.is_none());
        assert_eq!(
            all.iter()
                .filter(|r| r.subject.name() == "helicoid")
                .count(),
            1
        );
    }

    #[test]
    fn helicoid_supports_the_eight_coefficient_ids_and_answers_nothing_it_cannot(
    ) -> Result<(), String> {
        for id in Swept::ALL {
            assert!(Helicoid.supports(&format!("coeff_{}", id.name())), "{id:?}");
        }
        for id in [
            "coeff_",
            "coeff_f",
            "coeff_series",
            "so3_exp",
            "sen3_exp_n1",
        ] {
            assert!(!Helicoid.supports(id), "{id}");
            let rec = record(&[("theta", &[0.5]), ("n", &[0.5]), ("w", &[1.0])], &[])?;
            assert!(Helicoid.eval(id, &rec, Precision::F64).is_empty(), "{id}");
        }
        // `r` reads `n` and `w`, the others `theta`; a record without them is answered with nothing.
        let theta = record(&[("theta", &[0.5])], &[])?;
        let nw = record(&[("n", &[0.5]), ("w", &[1.0])], &[])?;
        for precision in [Precision::F64, Precision::F32] {
            assert!(Helicoid.eval("coeff_r", &theta, precision).is_empty());
            assert!(Helicoid.eval("coeff_k", &nw, precision).is_empty());
            let out = Helicoid.eval("coeff_r", &nw, precision);
            assert_eq!(out.keys().collect::<Vec<_>>(), ["d_branch", "value"]);
            let out = Helicoid.eval("coeff_e", &theta, precision);
            assert_eq!(out.keys().collect::<Vec<_>>(), ["d_branch", "value"]);
        }
        // An `f32` subject receives a binary32 or nothing: 0.1 is no binary32, and is not rounded.
        let tenth = record(&[("theta", &[0.1])], &[])?;
        assert!(Helicoid.eval("coeff_k", &tenth, Precision::F32).is_empty());
        assert_eq!(Helicoid.eval("coeff_k", &tenth, Precision::F64).len(), 2);
        Ok(())
    }

    #[test]
    fn helicoid_is_the_seeded_correct_kernel_bit_for_bit_at_both_precisions() -> Result<(), String>
    {
        // The two kernels share every candidate (`the_shipped_arms_are_the_seeded_arms_bit_for_bit`
        // and the two committed sweeps), so they agree on every record, in value and derivative.
        let correct = Seeded::generated();
        let mut n = 0;
        for precision in [Precision::F64, Precision::F32] {
            for id in Swept::ALL {
                let fn_id = format!("coeff_{}", id.name());
                for rec in records(id, precision)? {
                    let (h, c) = (
                        Helicoid.eval(&fn_id, &rec, precision),
                        correct.eval(&fn_id, &rec, precision),
                    );
                    assert!(
                        !h.is_empty() && same(&h, &c),
                        "{fn_id} {precision:?} {}",
                        rec.stratum
                    );
                    n += 1;
                }
            }
        }
        // Seven ids of 3420 records and `coeff_r` of 3426, over both precisions.
        assert_eq!(n, 7 * 3420 + 3426);
        Ok(())
    }

    #[test]
    fn no_record_of_any_stratum_scores_non_finite_at_either_precision() -> Result<(), String> {
        for precision in [Precision::F64, Precision::F32] {
            for id in Swept::ALL {
                for (stratum, max) in maxima(&Helicoid, id, precision)? {
                    let max = max.ok_or(format!("{id:?} {precision:?} {stratum}: non-finite"))?;
                    assert!(max.is_finite() && max >= 0.0, "{id:?} {stratum}: {max}");
                }
            }
        }
        Ok(())
    }

    /// The strata on which the prior (four terms below `z = 0.01`) has a smaller maximum than
    /// `helicoid`, as `(id, precision, stratum)`. It is `tf_tree`'s D12 for `a`, `b`, `c`; for `k`,
    /// `d`, `e` and `cos θ/2` D12 applied to a coefficient it was not defined for, and for `r`
    /// `(4 terms, s < 0.01)`, `θ = 0.2` and not D12's (0014 (draft) questions 10 and 29).
    fn prior_wins() -> Result<Vec<(&'static str, Precision, String)>, String> {
        let dir = corpus_dir()?;
        let (s64, s32) = (Series::<D1>::load(&dir)?, Series::<F32>::load(&dir)?);
        let prior64 = Seeded::uniform(&s64, d12(&s64)?);
        let prior32 = Seeded::uniform_f32(&s32, d12(&s32)?);
        let mut wins = Vec::new();
        for precision in [Precision::F64, Precision::F32] {
            for id in Swept::ALL {
                let prior = match precision {
                    Precision::F64 => &prior64,
                    Precision::F32 => &prior32,
                };
                let (h, p) = (
                    maxima(&Helicoid, id, precision)?,
                    maxima(prior, id, precision)?,
                );
                for (stratum, h) in h {
                    if p[&stratum].zip(h).is_some_and(|(p, h)| p < h) {
                        wins.push((id.name(), precision, stratum));
                    }
                }
            }
        }
        Ok(wins)
    }

    #[test]
    fn the_prior_beats_helicoid_on_nine_cos_half_strata_and_on_no_other() -> Result<(), String> {
        // `cos θ/2`'s objective is `theta:1e0`'s, where the two arms are one expression, so every
        // switch whose arms stay under it ties, and the tie goes to fewer terms and the larger
        // switch (`z < 5.6e-15`, `θ < 7.5e-8`; `z < 2.6e-6` at `f32`). Above it the exact arm
        // costs 1.4 to 1.9 `u` where the prior's series costs less. A stratum's maximum is what
        // the bars compare (D8), the sweep's objective one maximum over all of them: open,
        // 0014 (draft) question 30.
        let want = [
            (Precision::F64, "theta:1e-3"),
            (Precision::F64, "theta:1e-4"),
            (Precision::F64, "theta:1e-5"),
            (Precision::F64, "theta:1e-6"),
            (Precision::F64, "theta:1e-7"),
            (Precision::F64, "theta:1e-8"),
            (Precision::F32, "theta:1e-2@f32"),
            (Precision::F32, "theta:1e-3@f32"),
            (Precision::F32, "theta:dense@f32"),
        ];
        let mut got = prior_wins()?;
        got.sort_by(|a, b| a.2.cmp(&b.2));
        let mut want = want.map(|(p, s)| ("cos_half", p, s.to_string())).to_vec();
        want.sort_by(|a, b| a.2.cmp(&b.2));
        assert_eq!(got, want);
        Ok(())
    }
}
