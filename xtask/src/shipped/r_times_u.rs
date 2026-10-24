//! `helicoid:r-times-u`: the shipped `SO3::log` with one variable changed, the exact arm of `r` as
//! `r · (sign u)`, `r` rounded and then each component (`0090` decision 1). It is the arrangement
//! `0106` replaces with `r` as two terms and one rounding per component, kept as a subject so
//! `cargo xtask compare` can re-measure the choice (`0078`).
//!
//! On a series arm the twin is the shipped `SO3::log`, by construction: it dispatches on the masks
//! `log_ratio_scale` itself branches on, and on the exact arm forms `r` with the kernel's own
//! `exact_r`. That the variable is taken and moves bits is pinned
//! (`the_exact_arm_is_taken_and_moves`).
//!
//! Planted, as `helicoid:host-std` is: a twin is not a candidate, and no bar reads its rows.

use helicoid::__sweep as k;
use helicoid::{LieGroup, Quat, Tangent, SO3};
use helicoid_linalg::Precision;

use super::{quat_of, Helicoid};
use crate::conformance::corpus::Record;
use crate::conformance::subject::{Output, Registered, Subject};

/// `Log q` with the exact arm `r · (sign u)`, and whether that arm was taken.
pub(crate) fn r_times_u(q: Quat<f64>) -> ([f64; 3], bool) {
    let sign = 1.0f64.copysign(q.w);
    let (x, y, z) = (q.x, q.y, q.z);
    let (n2, w) = ((x * x + y * y) + z * z, sign * q.w);
    if k::log_ratio_takes_exact_arm(n2, w) {
        let r = k::exact_r(n2, w);
        ([r * (sign * x), r * (sign * y), r * (sign * z)], true)
    } else {
        let mut phi = [0.0; 3];
        SO3::from_quat_unchecked(q).log().write_dense(&mut phi);
        (phi, false)
    }
}

pub(crate) struct RTimesU;

impl Subject for RTimesU {
    fn name(&self) -> &str {
        "helicoid:r-times-u"
    }

    fn supports(&self, fn_id: &str) -> bool {
        fn_id == "so3_log" && Helicoid.supports(fn_id)
    }

    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        match (fn_id, precision, quat_of::<f64>(record, "q")) {
            ("so3_log", Precision::F64, Some(q)) => {
                Output::from([("phi".to_string(), r_times_u(q).0.to_vec())])
            }
            _ => Output::new(),
        }
    }
}

pub(crate) fn registered() -> Registered {
    let version = format!("{}@r-times-u", env!("CARGO_PKG_VERSION"));
    Registered {
        version: version.clone(),
        version_f32: version,
        no_f32: Some("`so3_log` has no `@f32` strata (`0016`)".to_string()),
        subject: Box::new(RTimesU),
        planted: true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::{corpus, corpus_dir};

    fn bits(o: &Output) -> Vec<(String, Vec<u64>)> {
        o.iter()
            .map(|(k, v)| (k.clone(), v.iter().map(|x| x.to_bits()).collect()))
            .collect()
    }

    /// The control: the exact arm, the variable, is taken on a share of the records and moves
    /// about half of them (992 of 1988 when written), and both arms occur.
    #[test]
    fn the_exact_arm_is_taken_and_moves() -> Result<(), String> {
        let dir = corpus_dir()?;
        let (mut series, mut exact, mut moved) = (0, 0, 0);
        for entry in corpus::manifest(&dir)?
            .into_iter()
            .filter(|e| e.fn_id == "so3_log")
        {
            for rec in corpus::read(&dir, &entry)? {
                let shipped = Helicoid.eval(&entry.fn_id, &rec, Precision::F64);
                let twin = RTimesU.eval(&entry.fn_id, &rec, Precision::F64);
                assert!(!shipped.is_empty(), "{}", rec.id);
                let Some(q) = quat_of::<f64>(&rec, "q") else {
                    return Err(format!("{}: no q", rec.id));
                };
                if r_times_u(q).1 {
                    exact += 1;
                    moved += usize::from(bits(&twin) != bits(&shipped));
                } else {
                    series += 1;
                }
            }
        }
        assert!(series > 0 && exact > 0, "{series} series, {exact} exact");
        assert!(moved > exact / 3, "{moved} of {exact}");
        Ok(())
    }
}
