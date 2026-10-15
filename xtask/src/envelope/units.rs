//! The swap test's exchangeable units (`0091`, draft): the records of a stratum whose inputs are
//! one draw. The generator writes each quaternion and its negative, the half a `Log` without its
//! flip gets wrong (`gen/so3.py`, `gen/sen3.py`), so a pair of records can be one rotation, and a
//! program answers both in one arithmetic. Swapping one record of such a pair and not the other is
//! no arrangement the null hypothesis allows, and counting both as evidence squares the p-value.
//!
//! Two records are one unit when every input is equal, a quaternion input (`q`, `q0`, `q1`: four
//! components) up to its sign. Read from the inputs alone, never from an error, so the grouping
//! cannot lean on the answer it tests.

use std::collections::BTreeMap;
use std::path::Path;

use serde::Deserialize;
use serde_json::Value;

use crate::conformance::corpus;
use crate::conformance::number::parse_hex;

/// Per `(fn, stratum)`, each record's unit: the smallest id among the records one with it.
#[derive(Default)]
pub(super) struct Units(BTreeMap<(String, String), BTreeMap<u64, u64>>);

/// The part of a corpus line a unit is read from. A record's id is its 0-based line number, which
/// `corpus::read` checks where the records are scored.
#[derive(Default, Deserialize)]
#[serde(default)]
struct Line {
    stratum: String,
    #[serde(rename = "in")]
    input: BTreeMap<String, Value>,
}

/// Whether input `name` is a quaternion, whose negative is the same rotation.
fn is_quaternion(name: &str, values: &[u64]) -> bool {
    let digits = name
        .strip_prefix('q')
        .is_some_and(|d| d.chars().all(|c| c.is_ascii_digit()));
    digits && values.len() == 4
}

/// `v`'s binary64 bit patterns, depth first; a string that is not a number is an error.
fn flatten(v: &Value, out: &mut Vec<u64>) -> Result<(), String> {
    match v {
        Value::String(s) => out.push(parse_hex(s)?.to_bits()),
        Value::Number(n) => out.push(n.as_f64().ok_or("a number past binary64")?.to_bits()),
        Value::Array(items) => {
            for item in items {
                flatten(item, out)?;
            }
        }
        other => return Err(format!("an input that is not a number: {other}")),
    }
    Ok(())
}

/// A record's inputs by name, as bit patterns.
type Inputs<'a> = Vec<(&'a str, Vec<u64>)>;

/// The record's inputs, each quaternion with the sign that makes its first nonzero component
/// positive: equal for two records exactly when they are one unit.
fn canonical(line: &Line) -> Result<Inputs<'_>, String> {
    let mut key = Vec::with_capacity(line.input.len());
    for (name, v) in &line.input {
        let mut bits = Vec::new();
        flatten(v, &mut bits)?;
        if is_quaternion(name, &bits) {
            let first = bits.iter().map(|&b| f64::from_bits(b)).find(|x| *x != 0.0);
            if first.is_some_and(|x| x < 0.0) {
                bits = bits
                    .iter()
                    .map(|&b| (-f64::from_bits(b)).to_bits())
                    .collect();
            }
        }
        key.push((name.as_str(), bits));
    }
    Ok(key)
}

impl Units {
    /// Every corpus file of `dir` that exists; none, each record its own unit, without a manifest.
    pub(super) fn read(dir: &Path) -> Result<Self, String> {
        let mut out = BTreeMap::new();
        if !dir.join("MANIFEST.json").is_file() {
            return Ok(Self(out));
        }
        for entry in corpus::manifest(dir)? {
            let path = dir.join(format!("{}.jsonl", entry.fn_id));
            if !path.is_file() {
                continue;
            }
            let text =
                std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
            let lines = text
                .lines()
                .map(|l| serde_json::from_str::<Line>(l).map_err(|e| e.to_string()))
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| format!("{}: {e}", path.display()))?;
            let mut first: BTreeMap<(&str, Inputs), u64> = BTreeMap::new();
            for (id, line) in (0_u64..).zip(&lines) {
                let key = canonical(line).map_err(|e| format!("{}: {e}", path.display()))?;
                let unit = *first.entry((line.stratum.as_str(), key)).or_insert(id);
                out.entry((entry.fn_id.clone(), line.stratum.clone()))
                    .or_insert_with(BTreeMap::new)
                    .insert(id, unit);
            }
        }
        Ok(Self(out))
    }

    /// `id`'s unit in `(fn_id, stratum)`; itself where the corpus was not read.
    pub(super) fn of(&self, fn_id: &str, stratum: &str, id: u64) -> u64 {
        self.0
            .get(&(fn_id.to_string(), stratum.to_string()))
            .and_then(|m| m.get(&id))
            .copied()
            .unwrap_or(id)
    }

    /// `n` records of `(fn_id, stratum)` in twin pairs, `2k` and `2k + 1` one unit.
    #[cfg(test)]
    pub(super) fn twins(fn_id: &str, stratum: &str, n: u64) -> Self {
        let ids = (0..n).map(|id| (id, id - id % 2)).collect();
        Self(BTreeMap::from([(
            (fn_id.to_string(), stratum.to_string()),
            ids,
        )]))
    }

    /// The units among `ids` of `(fn_id, stratum)`.
    pub(super) fn count(
        &self,
        fn_id: &str,
        stratum: &str,
        ids: impl Iterator<Item = u64>,
    ) -> usize {
        let mut units: Vec<u64> = ids.map(|id| self.of(fn_id, stratum, id)).collect();
        units.sort_unstable();
        units.dedup();
        units.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(q: [f64; 4], x: f64) -> Line {
        let hex = |v: f64| serde_json::Number::from_f64(v).map_or(Value::Null, Value::Number);
        Line {
            stratum: "s".into(),
            input: BTreeMap::from([
                ("q".to_string(), Value::Array(q.map(hex).to_vec())),
                ("x".to_string(), Value::Array(vec![hex(x)])),
            ]),
        }
    }

    #[test]
    fn a_quaternion_and_its_negative_are_one_unit_and_nothing_else_is() -> Result<(), String> {
        let q = [0.5, -0.5, 0.5, 0.5];
        let neg = q.map(|c| -c);
        let w0 = [0.0, -0.6, 0.0, 0.8];
        let lines = [
            line(q, 1.0),
            line(neg, 1.0),
            line(neg, 2.0),
            line([0.5, 0.5, 0.5, 0.5], 1.0),
            line(w0, 1.0),
            line(w0.map(|c| -c), 1.0),
        ];
        let keys = lines.iter().map(canonical).collect::<Result<Vec<_>, _>>()?;
        assert_eq!(keys[0], keys[1]);
        // Another input differing, or a quaternion differing other than by sign, is another draw.
        assert_ne!(keys[0], keys[2]);
        assert_ne!(keys[0], keys[3]);
        // A zero leading component: the sign is read from the first nonzero one.
        assert_eq!(keys[4], keys[5]);
        Ok(())
    }

    #[test]
    fn only_q_named_four_vectors_are_quaternions() {
        assert!(is_quaternion("q", &[0; 4]));
        assert!(is_quaternion("q1", &[0; 4]));
        assert!(!is_quaternion("q", &[0; 3]));
        assert!(!is_quaternion("qx", &[0; 4]));
        assert!(!is_quaternion("A", &[0; 4]));
    }
}
