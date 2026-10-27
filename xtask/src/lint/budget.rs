//! Check: the committed corpus stays under its budget (`0110`, amending `0071`): the bytes of
//! every tracked file under `conformance/corpus/`, the manifest included. The deep tier is
//! git-ignored and has none (`0088` decision 3).

use super::{File, Violation};

/// 80 MiB (`0110` decision 1).
pub(crate) const BUDGET: usize = 80 << 20;

const CORPUS: &str = "conformance/corpus/";

pub(crate) fn check(files: &[File]) -> Vec<Violation> {
    over(files, BUDGET)
}

fn over(files: &[File], budget: usize) -> Vec<Violation> {
    let total: usize = files
        .iter()
        .filter(|f| f.path.starts_with(CORPUS))
        .map(|f| f.text.len())
        .sum();
    if total <= budget {
        return Vec::new();
    }
    vec![Violation::new(
        "conformance/corpus/MANIFEST.json",
        1,
        format!(
            "[budget] the committed corpus is {total} bytes, over its {budget} (`0110`): a family \
             that would exceed it reduces its sample count (`PHASE1.md` §4.4)"
        ),
    )]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_cap_counts_only_the_committed_corpus_and_is_inclusive() {
        let files = [
            File::new("conformance/corpus/a.jsonl", "12345"),
            File::new("conformance/corpus/MANIFEST.json", "678"),
            File::new("conformance/corpus-deep/b.jsonl", "far over"),
            File::new("docs/c.md", "far over"),
        ];
        assert!(over(&files, 8).is_empty());
        let out = over(&files, 7);
        assert_eq!(out.len(), 1, "{out:?}");
        assert!(out[0].message.contains("8 bytes, over its 7"), "{out:?}");
    }
}
