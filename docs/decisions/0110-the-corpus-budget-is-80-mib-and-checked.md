# 0110: The corpus budget is 80 MiB, and checked

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** the PR that lands this record

> Amends [`0071`](./0071-the-corpus-budget-is-64-mib.md): the cap is 80 MiB, and `cargo xtask lint`
> enforces it.

## Context

[`0071`](./0071-the-corpus-budget-is-64-mib.md) holds corpus v1 under 64 MiB uncompressed
(67 108 864 bytes). Its only verification is "the corpus's size after `just corpus`", which no
recipe or check reads.

**The corpus is over it.** The tracked files under `conformance/corpus/` total 71 108 287 bytes
(67.81 MiB) on `main` at `7ed33b6`. The total per commit that touched the corpus, from
`git ls-tree -r -l <commit> conformance/corpus`:

| commit | change | MiB |
|---|---|---|
| `f2e3851` | `0082` | 61.21 |
| `cb27cee` | `0085`, `so3_from_matrix` scored | 62.32 |
| `67b3806` | `0101`, `svd3` | 63.52 |
| `0e06ae9` | `0102`, `nearest_rotation` | 63.86 |
| `75897e1` | `0103`, `so3_inv_act` and `se3_inv_act` | **67.60** |
| `7ddc3c5` | `0105`, `real_atan2` strata | 67.81 |

`0103`'s 3.74 MiB crossed the cap, and nothing failed. `PHASE1.md` §4.4's rule, "a family that would
exceed it reduces its sample count by an edit to this list, never silently", was broken silently.

The next family is `eig_sym`, the symmetric eigensolver at 4 ≤ N ≤ 16 that locus-calib's D4 consumes:
three ids,
N ∈ {4, 9, 16}, at both precisions. At about 22 bytes per input entry and 38 per output (the hex and
30-digit formats), a record is about 22N² + 38(N² + N + 1) + 100 bytes: 1.25 KB at N = 4, 5.3 KB at
N = 9 and 16.1 KB at N = 16. Its 43 strata per id would hold 32, 16 and 4 records at the committed tier. The family is then about 8.2 MB, and the corpus about 75.6 MiB.

## Decision

1. **Corpus v1 stays under 80 MiB uncompressed** (83 886 080 bytes). That is the sum of every
   tracked file under `conformance/corpus/`, the manifest included.
   - `PHASE1.md` §4.4's sentence changes to say so and cites this record.
   - `0071`'s rule stands: a family that would exceed the cap reduces its sample count by an edit
     to that list.
2. **`cargo xtask lint` checks it** (`lint::budget`). A required check fails a pull request that
   takes the committed corpus over the cap, naming the total and the cap.
   - The check reads the tracked files `lint` already loads, so it adds no I/O.
   - It does not apply to the deep tier, which is git-ignored and has no budget
     ([`0088`](./0088-the-domination-bar-must-tell-worse-from-equal.md) decision 3).
3. **`eig_sym`'s allocation is 8.2 MB**, which leaves about 4.4 MiB for what follows. Its own
   record states the counts.

## Rationale

- **Why `lint`, not `gen verify`.** `corpus-verify` runs inside `corpus-check`, which is not a
  required check, and `0108` lets a pull request skip it. `lint` is required and already reads
  every tracked file as text.
- **Why 80, not 64 with fewer records.** `0071` names the cap's purpose: clone size and
  `corpus-check`'s time.
  - `0108` took a warm `corpus-check` to about a minute, from 17.
  - With `xz` the corpus is 9.1 MB (`0088` finding 5): the size is digits, not information.
  - Holding 64 MiB would take 3.8 MiB back out of `0103`'s or `0105`'s strata. It would also leave
    `eig_sym` nothing.
- **Why not zeros written as `0` first** (`0037` item 3, draft, about 9.4 MiB). It changes every
  id's file and the reader. It is worth its own record, and it does not depend on this one. A cap
  without a check would overrun again either way.
- **Why not a separate cap per family.** The precedent (`0071`'s Sim(3) cut) is one budget with
  allocations argued in the record that spends it. A per-family table would need an edit for every
  stratum.

## Consequences

- The cap now fails loudly: the next family that crosses it either cuts its records or writes a
  record that moves the cap.
- `0071` gets an amendment banner that cites this record. Its status is unchanged.

## Implementation plan

1. Decisions 1 and 2: `xtask/src/lint/budget.rs`, its line in `CHECKS`, its fixture in
   `lint::tests::every_check_runs`, `PHASE1.md` §4.4 and `0071`'s banner. Verified by
   `budget::tests`, `lint::tests::real_tree_is_clean` and `just lint`.

## Open questions

None.
