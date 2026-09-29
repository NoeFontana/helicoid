# 0001: Record architectural decisions in `docs/decisions/`

**Status:** implemented
**Owner:** @NoeFontana
**Implementation:** initial commit

## Decision

Architectural Decision Records rooted at `docs/decisions/`, as in `tf_tree`:

- One document type, skeleton in [`template.md`](./template.md).
- Four statuses on the document itself: `draft → ready → implemented → superseded by NNNN`.
- Two gates: `draft → ready` (open questions resolved, plan concrete) and an `implemented`
  immutability lock.
- Sequential four-digit numbering, append-only.
- [`README.md`](./README.md) states the lifecycle and gates.

Anything touching the public API, a convention, a formula, a crate boundary, the dependency budget,
the corpus schema or the release process starts as a `draft` here.
