# 0030: A consumer comparison is a trigger `0021` does not list

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** records a finding and five open questions; no code change.

## Context

[`0021`](./0021-the-cholesky-finite-guard-stays.md) closed the question of `chol`'s per-entry
finite guard against four bit-identical variants, set the bar for rewriting the kernel at **15% at
`N = 6`**, and in its standing rule 3 deferred re-measurement to a **trigger**: "when `just bench`
exists (`PHASE1.md`) and can time aarch64, wasm32 or `thumbv7em`, or when a lane `Real` exists".
Every variant it weighed was one of ours.

A trigger arrived that the rule does not list: **an external subject**. Evaluating `helicoid` as the
math backing for `locus-tag` put `chol` beside the two small symmetric-positive-definite solves that
consumer runs today — `nalgebra`'s `Cholesky` at `N = 6` in its pose Levenberg–Marquardt solvers, and
a hand-rolled `[f64; 64]` factor-and-solve at `N = 8` in its line-segment refinement. That consumer's
own micro-optimization protocol forbids regressing a pipeline stage, so "is `chol` not slower?" is a
gate on the migration, not a curiosity. `0021`'s bar answers "should we rewrite ours?" and cannot
answer "are we behind someone else's?".

No document prescribes a benchmark for `helicoid-linalg`: `PHASE1.md` §9 scopes `criterion` to
`crates/helicoid/benches` **from Phase 3**, and [`0019`](./0019-a-cholesky-solve-without-the-transpose.md)
already recorded its own measurement as advisory for that reason. The measurement below is advisory
on the same grounds.

## Measurement

Advisory, not a gate. Scratch crate outside the repository, not committed; path dependency on
`helicoid-linalg`, registry dependency on `nalgebra` 0.35 (`default-features = false`,
`features = ["std"]`). `--release`, `lto = "fat"`, `codegen-units = 1`, no `target-cpu` and no
`rustflags`, so `fma` is off for both subjects (D16, `0007`). One AMD EPYC-Milan host, 8 vCPU, KVM
guest (4 cores x 2 threads), `rustc` 1.98.1; pinned with `taskset`. `N = 6` over **ten runs on four
cores**, `N = 3` and `N = 8` over **six runs on two cores**. Another session held about one core
through the sweep (load average 2.5 to 4 on 8 cores, nothing like the load-12 sweep `0021`
discards); the **min** statistic is the defence, and the `nalgebra` rows' tightness is the evidence
it worked. 2048
pre-generated dense SPD matrices per batch (`A = BᵀB + N I`, column-major, `B` entries in `[-1, 1)`,
one fixed LCG seed per `N`), mask set on every one. `black_box` on the input **and on the whole
factor**; the statistic is the **min** ns per factorisation over 41 repetitions, as in `0021`.

**Keeping the whole factor live is load-bearing, and two earlier versions of this harness did not.**
Consuming one entry of `L` lets the optimizer delete most of the factorisation: that harness reported
`noguard` *slower* than the guarded kernel and a verbatim copy of `chol` 26% slower than `chol`, both
impossible. The rows below checksum every entry of `L` plus the mask. A separate row times the
checksum alone so the others can be read net of it.

Factorisation only, `f64`, range over the six runs:

| `N` | `chol` | `nalgebra` | ratio | checksum probe |
|---|---|---|---|---|
| 3 | 14.36–14.46 ns | 9.68–9.72 ns | 1.48–1.49× | 1.50–1.58 ns |
| 6 | 93.69–107.46 ns | 31.30–36.73 ns | **2.6–3.4×** | 6.53–8.46 ns |
| 8 | 165.95–174.38 ns | 121.32–126.02 ns | 1.32–1.44× | 14.78–15.29 ns |

`nalgebra` at `N = 6` is bimodal by core — 31.3 to 31.7 ns on two of the four, 36.5 to 36.7 ns on
the other two — which is an SMT-sibling effect the prescribed harness's per-host `δ` would capture
and this one only notices. `nalgebra`'s row builds a fresh lower-triangular
matrix to read the factor, work `chol` does not do, so the gap is if anything understated.

Factor **and** solve, separate harness, same host and protocol, `L` live through the solve:

| Solve | consumer's today | `chol` + `chol_solve` | ratio |
|---|---|---|---|
| `N = 6` factor + solve | 53.0 ns (`nalgebra`) | 118.7 ns | 2.24× |
| `N = 8` factor + solve | 137.5 ns (hand-rolled) | 193.5 ns | 1.41× |
| `N = 6` solve only, factor given | 27.2 ns (`nalgebra`) | 26.8 ns | **0.98×** |

Accuracy, max relative residual `‖Ax − b‖/‖b‖` over the same 256 systems, favours `chol`:
4.06e-16 against `nalgebra`'s 4.40e-16 at `N = 6`, and 4.16e-16 against the hand-rolled 5.11e-16 at
`N = 8`. The hand-rolled kernel also tests its pivot against an **absolute** `1e-12`, where `chol`'s
mask is scale-free.

### The attribution failed, and the harness is why

Two bit-identical candidates were measured against the shipped kernel at `N = 6`: reading `a` from a
plain column-major array instead of through `Matrix::get`, which is bounds-checked in release; and an
explicit indexed reduction in place of the iterator fold. Both were checked `(L, mask)`-bit-identical
to `chol` on every matrix in the batch at `N = 3, 6, 8`. A third — hoisting the pivot's reciprocal
so the entries multiply instead of divide — was rejected before timing: it differs from `chol` in
thousands of entries per batch, so D16 disqualifies it whatever it costs.

| Variant at `N = 6`, ten runs | vs shipped |
|---|---|
| verbatim copy of the shipped body (**control, must be 0%**) | **−0.4% to +11.3%** |
| `a` read without `Matrix::get` | −22.1% to +5.0% |
| plus an indexed reduction | −33.4% to −5.7% |
| guard deleted (`0021`'s `noguard`; breaks the contract, reference only) | −11.1% to −27.8% |

Two things to take from this table, in opposite directions.

**It cannot settle the question.** The shipped kernel's own absolute figure swings 93.69 to 107.46 ns
run to run, and the control — a verbatim copy of the body being measured — lands between −0.4% and
+11.3% of it. The control's deviation is also systematically **non-negative**, so it is not pure
noise but a real artifact of where the copy inlines; either way, a harness whose zero reads up to
+11% cannot adjudicate a 15% bar. `noguard` agrees with `0021` in sign and not in magnitude
(−11% to −28% here against its −3% to −6%), which is evidence about ad-hoc harnesses rather than
about the guard.

**But the indexed reduction may well clear the bar.** Its entire ten-run range is below zero, with a
midpoint near −20%, and it is bit-identical. That is a real candidate — and precisely why it must not
be taken on this evidence: `0021`'s bar has three clauses, and this sweep tested **one**. It measured
the `f64` success path only. Nothing here touches failing input or `Dual`, where `0021`'s own
row-major variants lost 38% to 41% — the case that carries every Jacobian. A rewrite blessed on the
success path alone could regress the arm the library exists to serve.

`chol_solve` is at parity with `nalgebra`, so whatever the cause, it is in the factorisation.

## Decision (proposed)

1. **Do not touch `chol` yet.** `0021`'s bar stands and no variant has been *shown* to clear it: the
   indexed reduction is a live candidate on the `f64` success path and untested on the two clauses
   that matter most. "Instrument before algorithm"
   ([`0006`](./0006-the-instrument-comes-first.md)) binds `helicoid`'s own latency, not only its
   accuracy: a perf change whose effect this repository cannot resolve does not ship, however
   promising the midpoint looks.
2. **An external subject is a trigger.** Amend `0021` rule 3 to list it alongside `just bench` and a
   lane `Real`: a consumer's incumbent kernel measured faster than ours reopens the measurement,
   though not the guard.
3. **Extend `PHASE1.md` §9's harness to `helicoid-linalg`** rather than inventing a second one:
   `criterion` (pinned) in `crates/helicoid-linalg/benches`, `cargo xtask bench-gate` with §9's
   paired bootstrap 95% CI of the ratio and the host's A/A noise floor `δ` from `bench-gate --aa`
   recorded in `baseline/HOST.md`. The control row above is what `δ` is for.
4. **Third-party comparison rows live in a workspace-excluded crate** under `runners/`, as the oracle
   runners do, so no comparison dependency enters the root `Cargo.lock` (`PHASE1.md` §0). Reported,
   never gated (§9).

## Open questions

1. **Does §9's "from Phase 3" bind?** The consumer evaluation needs a `helicoid-linalg` bench now,
   and Phase 3's group code does not exist. Amend §9, or carry the harness as Phase 2 work?
2. **Why is the gap worst at exactly `N = 6`** — 2.6–2.8×, against 1.5× at `N = 3` and 1.3–1.4× at
   `N = 8`? Non-monotonic in `N` points at a threshold effect, not at a constant per-entry cost, and
   `N = 6` is the size every SE(3) consumer uses.
3. **Is `nalgebra` a legitimate subject?** §9 names `tf_tree_math` and sophus-rs for bench rows, and
   [`0010`](./0010-seeded-from-tf-tree-math.md)'s runners are accuracy oracles. A latency-only
   subject that answers no corpus id may need its own concept.
4. **How much of the residual is D16?** The pinned left-to-right summation order forbids
   reassociation that `nalgebra` is free to take. Testable under the instrument as a reference row
   only — a variant free to reassociate is not a candidate, it is a measurement of what determinism
   costs.
5. **Do `0021`'s magnitudes survive?** Re-measure its Appendix variants under the instrument, or let
   the record stand as a decision whose numbers are no longer the best available?

## Consequences

- The `locus-tag` migration's first gate is **blocked on question 1**, not on `chol`: until a
  resolvable measurement exists, no claim that `chol` does not regress that consumer can be made.
  Its `N = 8` kernel is the narrower target (1.3–1.4×, and it buys a scale-free pivot and a better
  residual); its `N = 6` pose path is the one that would regress.
- Nothing here changes `NUMERICS.md` §15, `chol`'s rustdoc, or the `PHASE2.md` §0.0 row.

## Implementation plan

1. This record.
2. On `ready`: the §9 amendment and the `helicoid-linalg` bench harness with `bench-gate --aa` and
   `baseline/HOST.md`.
3. On `ready`: the excluded comparison runner.
4. On `ready`: re-measure `0021`'s variants and the two above under the instrument — **all three
   clauses of the bar**, success path, failing input and `Dual` — then decide on `chol`. The indexed
   reduction is the first candidate to put through it.
