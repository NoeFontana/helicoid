# 0071: the corpus budget is 64 MiB

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #139

> **Amended by [`0110`](./0110-the-corpus-budget-is-80-mib-and-checked.md), 2026-10-27.** The cap is
> 80 MiB, and `cargo xtask lint` enforces it: the corpus had reached 67.8 MiB, past this record's 64,
> with nothing failing.

## Context

`PHASE1.md` §4.4 holds corpus v1 under 50 MB uncompressed, and a family that would exceed it
"reduces its sample count by an edit to this list, never silently". The corpus reached 52.2 MB
(49.8 MiB, so the budget reads as MiB) with `sim3_ref.jsonl`, the 2.9 MB of ray references the
Sim(3) sweep reads ([`0070`](./0070-sim3-a-structured-jacobian-and-a-series-on-delta.md), #137).
The eight Sim(3) ids at SE_N(3)'s density — six records per stratum, five of them `7 x 7` matrices,
both precisions — measured 17.3 MB. At the 0.2 MB left, Sim(3) would get about one record per
cell of its `σ × θ` grid, and its no-regress bar would read one draw per stratum.

## Decision

**Corpus v1 stays under 64 MiB uncompressed** (67 108 864 bytes). `PHASE1.md` §4.4's sentence
changes to say so; the rule that a family over budget reduces its sample count, by an edit to the
list, stands, and Sim(3) is its first use: at six records per cell the family measured 17.3 MB
and the corpus 69.5 MB, so it takes **four** per cell (eight per ray): the corpus is then 64.2 MB (61.2 MiB), 2.9 MB under.

## Rationale

- The budget guards clone size and `just corpus-check`'s time, not a property of the
  measurement. 14 MiB is a quarter more of both; one record per cell is a no-regress bar that
  cannot tell a defect from a draw.
- SE_N(3)'s six is the precedent `PHASE1.md` §4.4 argues as a budget trade-off; Sim(3)'s five
  `7 x 7` ids make the same trade one step lower.

## Consequences

- `PHASE1.md` §4.4 reads 64 MiB. Later phases' ids budget against what is left.

## Implementation plan

1. With the Sim(3) corpus ids (#139). Verified by the corpus's size after `just corpus`.

## Open questions

None.
