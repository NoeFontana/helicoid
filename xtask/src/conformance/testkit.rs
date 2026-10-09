//! Test-only: a seeded generator, hand-built records, and the subjects the harness checks itself
//! with. None of these is ever registered (`subject::registry`), so no envelope row can read them.

use std::collections::BTreeMap;
use std::path::PathBuf;

use helicoid_linalg::Precision;

use super::corpus::{Record, Tensor};
use super::metric::{self, Rule};
use super::number::Decimal;
use super::subject::{Output, Subject};

/// splitmix64, the generator's stream (`docs/maths/error-analysis.md` EA.23).
pub(crate) fn splitmix(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// A directory under the system temp dir, removed when the test ends.
pub(crate) struct Scratch(pub(crate) PathBuf);

impl Scratch {
    pub(crate) fn new(tag: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("xtask-conformance-{}-{tag}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        Self(dir)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A record with the given binary64 inputs and decimal reference strings.
pub(crate) fn record(
    inputs: &[(&str, &[f64])],
    reference: &[(&str, &[&str])],
) -> Result<Record, String> {
    fn tensor<T>(data: Vec<T>) -> Tensor<T> {
        Tensor { shape: None, data }
    }
    let mut record = Record {
        id: 0,
        stratum: "s".into(),
        inputs: BTreeMap::new(),
        reference: BTreeMap::new(),
    };
    for (key, values) in inputs {
        record.inputs.insert((*key).into(), tensor(values.to_vec()));
    }
    for (key, values) in reference {
        let data = values
            .iter()
            .map(|s| s.parse::<Decimal>())
            .collect::<Result<_, _>>()?;
        record.reference.insert((*key).into(), tensor(data));
    }
    Ok(record)
}

/// The harness's oracle: the reference correctly rounded at the precision asked, then moved `ulps`
/// places up there.
pub(crate) struct Perfect {
    name: &'static str,
    ulps: u32,
}

impl Perfect {
    pub(crate) fn exact() -> Self {
        Self {
            name: "perfect",
            ulps: 0,
        }
    }

    pub(crate) fn next_up() -> Self {
        Self {
            name: "next-up",
            ulps: 1,
        }
    }
}

impl Subject for Perfect {
    fn name(&self) -> &str {
        self.name
    }

    fn supports(&self, _: &str) -> bool {
        true
    }

    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        let up = |x: f32| (0..self.ulps).fold(x, |x, _| x.next_up());
        let round = |d: &Decimal| match precision {
            Precision::F64 => (0..self.ulps).fold(d.to_f64(), |x, _| x.next_up()),
            Precision::F32 => f64::from(up(d.to_f32())),
        };
        let field =
            |(k, t): (&String, &Tensor<Decimal>)| (k.clone(), t.data.iter().map(round).collect());
        let mut out: Output = record.reference.iter().map(field).collect();
        if matches!(metric::rule(fn_id), Some(Rule::Adjugate)) {
            // A singular reference's `inv` is zeros standing for "none": the right answer is a
            // non-finite inverse (`0061` decision 3).
            let singular = record
                .reference
                .get("det")
                .is_some_and(|t| t.data.iter().all(|d| d.mant == 0));
            if singular {
                out.insert("inv".to_string(), vec![f64::NAN; 4]);
            }
            return out;
        }
        if !matches!(metric::rule(fn_id), Some(Rule::Roots)) {
            return out;
        }
        // A root-finder's answer is slots and their mask: every real root in a valid slot, the
        // complex ones' slots clear (`0056` decision 3). Real is read off the reference, before
        // rounding moves an exact zero.
        let re = out.get("re").cloned().unwrap_or_default();
        let im = record
            .reference
            .get("im")
            .map(|t| t.data.as_slice())
            .unwrap_or_default();
        let real = |i: usize| im.get(i).is_some_and(|d| d.mant == 0);
        Output::from([
            (
                "roots".to_string(),
                (0..re.len())
                    .map(|i| if real(i) { re[i] } else { 0.0 })
                    .collect(),
            ),
            (
                "valid".to_string(),
                (0..re.len())
                    .map(|i| if real(i) { 1.0 } else { 0.0 })
                    .collect(),
            ),
        ])
    }
}

/// A subject given by a closure over the record.
pub(super) struct Fixed {
    name: String,
    answer: Box<dyn Fn(&Record) -> Output>,
    only: Option<String>,
}

impl Fixed {
    pub(super) fn new(name: &str, answer: impl Fn(&Record) -> Output + 'static) -> Self {
        Self {
            name: name.into(),
            answer: Box::new(answer),
            only: None,
        }
    }

    /// Supports one function id only.
    pub(super) fn only(mut self, fn_id: &str) -> Self {
        self.only = Some(fn_id.into());
        self
    }
}

impl Subject for Fixed {
    fn name(&self) -> &str {
        &self.name
    }

    fn supports(&self, fn_id: &str) -> bool {
        self.only.as_deref().is_none_or(|f| f == fn_id)
    }

    fn eval(&self, _: &str, record: &Record, _: Precision) -> Output {
        (self.answer)(record)
    }
}
