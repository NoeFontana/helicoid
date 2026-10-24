# 0108: The corpus check proves each thing once

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record's PR: `gen verify`, `gen key`, `gen/memo.py`, `tests/parallel.py`, the `corpus-check` and `corpus-test` CI jobs, `corpus-nightly.yml`.

## Context

`just corpus-check` regenerates the whole mpmath corpus, 106 546 records over 89 ids, into a
temporary directory and compares it byte for byte with `conformance/corpus/`
([`PHASE1.md`](../PHASE1.md) §4 item 4), then runs the generator's unit tests. It ran on every pull
request and every push to `main`.

- **CI.** On #205's run (`ubuntu-latest`, 4 vCPU), generation took 16.3 min and the 293 tests,
  serially, 3.2 min. The same job read 10.5 min on another run of the same day.
- **The work.** Profiled per task on an 8-core host under load: 4 985 tasks, one per
  (id, stratum), about 2 000 CPU-s. The longest task is 18 s (`so3_log theta:dense`), so the run is
  throughput-bound; the pool already spreads it. Sim(3)'s four Jacobians are 27% of it, SE_N(3)'s
  at `n = 3` 15%. A record's generation-time cross-check costs 1–8× its evaluation.
- **The files beside the registry's.** `coeff_series`, `coeff_switch_ref`, `sim3_series` and
  `sim3_ref` were built serially in the parent after the pool, about 37 s on this host.
- **gmpy2.** mpmath's gmpy2 backend wrote the same bytes on five ids' first six strata and ran
  0.85–1.9× as long: at 120 digits the cost is the interpreter's, not the integers'. Not taken.

The check proves one thing: the committed bytes are what this generator writes. That is a function
of the generator's sources, the interpreter and the corpus alone, and
[`manifest.py`](../../conformance/generate/gen/manifest.py) already hashes the first and the last.
A pull request that touches neither proved nothing new in its 17 minutes; #205, which touched only
`fit/`, outside the generator, was one.

## Decision

1. **Verify in seconds; regenerate only what has not been proven.**
   - `python -m gen verify` (`just corpus-verify`): the committed manifest names this generator,
     interpreter minor version, mpmath, precisions and seed, and lists exactly the files present,
     each with its sha256.
   - `python -m gen key` (`just corpus-key`): one sha256 of the generator identity, the exact
     interpreter (`sys.version`, machine) and the manifest's bytes, which hash every file.
   - CI's `corpus-check` verifies, then on a pull request looks the key up in the Actions cache.
     A hit means a full check of the same inputs passed: the regeneration is skipped. The marker is
     written only after a full check passes.
   - `main` always runs every task (`just corpus-regen fresh`), and `corpus-nightly.yml` runs
     `just corpus-check fresh` with no cache at all.
2. **A memo of tasks** ([`gen/memo.py`](../../conformance/generate/gen/memo.py)). Each (id,
   stratum) task, and each of the four files beside the registry's, is keyed on:
   - the files of `gen` in which a function started during the task, recorded by
     `sys.monitoring` (one event per code object per task, then disabled);
   - the static import closure of those files within `gen`, for module-level values read without
     running the module's code; an import of `registry`, `cli` or `__main__` adds nothing, since
     their only effect on a task is the spec below (a function of theirs that runs is a
     dependency as any other);
   - `describe` of the spec without its strata, and of the stratum: qualified names, code with the
     file name removed, closure values, defaults and dataclass fields. An object it cannot describe
     is not memoised;
   - Python's full version, mpmath's, the machine and `HELICOID_CORPUS_DEEP`.

   A task whose key and every dependency's hash are unchanged replays its records. With `--audit p`
   (2% in every recipe) a share of the replays, chosen by the run's seed, runs again, and a byte
   that differs fails the run. `--fresh` runs every task and rewrites its entry. The memo is
   `conformance/generate/.memo/`, git-ignored; CI restores the newest one `main` saved.
3. **The four files are tasks** of the same pool, started first.
4. **The generator's tests run a test class per core** (`python -m tests.parallel`,
   `just corpus-test`), and CI runs them as their own job, `corpus-test`, skipped on a pull request
   where they passed for the same corpus key, tests and `fit`.

## Rationale

**Skip by key, not by path filter.** A path filter on `conformance/generate/**` would have run
#205's check for an edit to `fit/`, and would skip a hand edit to the corpus that leaves the
generator alone unless the filter also names the corpus. The key is the generator's own identity
and the manifest: what the check is a function of, no more and no less. A cache entry can only
skip a proof already made on identical inputs; a pull request cannot write `main`'s entries.

**Dynamic dependencies, bounded by a static closure.** A file list declared per id would rot. What
a task ran is exact for code; what it read at module level is not, since `from .strata import
N_RANDOM` binds a value whose module may run nothing during the task. The static closure covers
that, at the price of some needless reruns (every task imports `corpus.py`'s imports). The registry
is the one hub left out, because the spec it builds is in the key; with it, adding an id would
rerun every task.

**The audit and the fresh runs are the soundness argument.** The memo is an optimisation of a check
that still runs in full on `main` and nightly. A dependency the rules miss, such as a file read at
run time, which no generator module does today, shows up as an audit failure on a pull request or
a mismatch on `main`, never as a corpus that passes everywhere and is wrong.

**Not taken.** gmpy2 (above). Weaker cross-checks: they are what the corpus's trust rests on, and
with the memo they run only for the tasks that change. Sharding the full run across jobs: with the
skip and the memo, a full run happens on `main` and nightly, off any pull request's path.

## Consequences

Measured on the 8-core host, under load 14:

| Run | Before | After |
|---|---|---|
| `just corpus` or `corpus-check`'s generation, every task | 15.4 min | the same (`fresh`) |
| the same, nothing changed (memo) | — | 32 s, 4 892 replayed, 97 audited |
| after an edit to `gen/real.py` | — | 15 s, 101 tasks run |
| `just corpus-test` | 190 s on CI, serial | 50 s, bounded by `GenerationTest` |
| `just corpus-check` | about 17 min | 1 min 10 s with the memo warm |

- On CI, a pull request that touches neither the generator nor the corpus spends the job's setup
  and `corpus-verify`; one that does replays from `main`'s memo.
- The memo is 87 MB uncompressed; CI keeps one per run in the Actions cache, newest restored.
- `main`'s `corpus-check` still takes the full time, off the merge path, and fills the memo.
- The generator's identity now includes `memo.py`: an edit to it reruns every task.
- `corpus-check` stays a check that is not required, as before.

## Implementation plan

1. This record, `gen verify` and `gen key`, the memo, the four files as tasks, the parallel test
   runner, the recipes and the CI jobs. Verified by `tests/test_memo.py`, `just corpus-check` (the
   corpus byte-identical; only the manifest's generator hash moves) and the first CI runs: a fresh
   run on `main`, then a pull request that skips.

## Further work

- Shard `main`'s fresh run across jobs by measured task cost, if its time comes to matter.
- `just corpus-deep` takes no memo; the deep tier draws rarely.
