//! The swap test's exchangeable units (`0091`, draft): the records of a stratum that are one draw.
//! The conformance run writes each record's draw beside its error (`corpus::draws`), so the
//! grouping is the one the errors were scored under; a records file in `0088`'s format has none,
//! and each of its records is its own draw.

use std::collections::BTreeMap;

/// Per `(fn, stratum)`, each record's draw, the smallest id among the records one with it: kept
/// only where it is not the record's own id, which is every record's draw by default.
#[derive(Default)]
pub(super) struct Units(BTreeMap<(String, String), BTreeMap<u64, u64>>);

impl Units {
    /// Record `id` of `(fn_id, stratum)` is of draw `draw`.
    pub(super) fn insert(&mut self, fn_id: &str, stratum: &str, id: u64, draw: u64) {
        let key = (fn_id.to_string(), stratum.to_string());
        self.0.entry(key).or_default().insert(id, draw);
    }

    /// `(fn_id, stratum)`'s map from record to draw; `None` where no draw was read.
    pub(super) fn stratum(&self, fn_id: &str, stratum: &str) -> Option<&BTreeMap<u64, u64>> {
        self.0.get(&(fn_id.to_string(), stratum.to_string()))
    }

    /// `n` records of `(fn_id, stratum)` in twin pairs, `2k` and `2k + 1` one draw.
    #[cfg(test)]
    pub(super) fn twins(fn_id: &str, stratum: &str, n: u64) -> Self {
        let mut out = Self::default();
        for id in 0..n {
            out.insert(fn_id, stratum, id, id - id % 2);
        }
        out
    }
}
