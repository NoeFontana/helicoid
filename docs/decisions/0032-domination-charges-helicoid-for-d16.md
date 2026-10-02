# 0032: Domination charges `helicoid` for D16

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** records a measurement and proposes a bar change; no code change.

## Context

The envelope's first run over the real corpus reports **75 domination failures**, and `so3_log` is
8 of them — the one id where the candidate wins no stratum at all (`PHASE1.md` §8 row). Those 8 were
the obvious place to start looking for an algorithm to improve.

There is no algorithm to improve. `so3_log` is the same program in `helicoid` and in both oracles:
`Log` goes through the quaternion `atan2` (D5), `φ = r u` with `r = 2 atan2(n, w)/n`, `n² = (x² + y²) + z²`.
The entire gap is **which `atan2` ran**. `helicoid` routes every transcendental through the `libm`
crate so that outputs are a function of the input bits alone (D16,
[`0018`](./0018-libm-arch-is-bit-identical-for-exact-operations.md)); sophus-rs is Rust `std`, so its
`atan2` is the host's glibc.

`error-analysis.md` **EA.13(d) predicted exactly this**, including the magnitude and the conclusion:
`libm`'s $\nu$ is $1.96$ on `atan2` against glibc's $\le 1.00$, so "a `Log` stratum whose error is
dominated by its final `atan2` may start up to $\approx 1u$ behind an oracle with a nearly correctly
rounded one: **a domination failure there can be the library's and not the algorithm's**". What was
missing is the measurement on the corpus, and the count of strata it accounts for.

## Measurement

Advisory, as `0019`'s and `0022`'s were. The corpus's own 30-digit references,
`metric.rs`'s metric formed at 40 digits in `mpmath` ($F = \lVert \hat y - y\rVert / (\max(\lVert y
\rVert, 2^{-1022})\,u)$, 2-norm over the field), binary64. One variable changed: `atan2`, the `libm`
crate's against the host's (glibc 2.39, x86_64). Everything else — the flip, `n²`'s association,
`sqrt`, the division, the three multiplies — held fixed.

The two `atan2`s disagree on **692 of 4994** `so3_log` records (13.9%), and on 22 to 40 of the 128
records of each losing stratum (430 of 1602 on `theta:dense`).

| stratum | `libm` crate | host glibc | price | ratio |
|---|---|---|---|---|
| `theta:dense` | 2.9010 | 2.1605 | **+0.7405** | 1.34 |
| `theta:1e-1` | 2.2928 | 1.6643 | +0.6285 | 1.38 |
| `theta:1e-2` | 2.2160 | 1.9876 | +0.2284 | 1.11 |
| `theta:1e-3` | 2.3639 | 1.6312 | +0.7327 | **1.45** |
| `theta:1e-4` | 2.3651 | 2.0877 | +0.2774 | 1.13 |
| `theta:1e-5` | 1.7959 | 1.6585 | +0.1374 | 1.08 |
| `theta:1e-6` | 2.5411 | 1.9923 | +0.5488 | 1.28 |
| `theta:1e-7` | 2.1743 | 1.9396 | +0.2347 | 1.12 |

**The columns are the two subjects.** The `libm`-crate column reproduces the harness's
`seeded:correct` row on **8 of 8** strata; the glibc column reproduces `sophus_rs`'s on **7 of 8**
(`theta:1e-7`, where glibc gives 1.8492 against sophus-rs's 1.9396, so that stratum has a second
difference of its own and is not purely `atan2`-bound). `tf_tree_math`, which also uses the `libm`
crate, is **bit-identical to `helicoid` on every `so3_log` stratum** — the clean control: same
`atan2`, same answer.

So D16 costs `so3_log` at most **0.74 u**, inside EA.13(d)'s "up to about 1u", and it costs exactly
the 8 strata the bar fails.

## Decision (proposed)

**Domination must not charge the candidate for an oracle's `libm`.** The bar as written compares
`helicoid` on the `libm` crate with an oracle on its platform's, and `PHASE1.md` §8 states it as a
property of the algorithm. On an `atan2`- or `sin_cos`-dominated stratum it is partly a property of
the two `libm`s, and the candidate is the one that chose reproducibility.

Three shapes, narrowest first:

1. **Report, do not fail.** A third verdict beside pass and fail — `libm-bound` — for a stratum
   where the candidate is within the per-function `libm`-vs-glibc penalty of the best oracle
   (EA.13(d)'s $\nu$ ratio: $\approx 2$ for `atan2`, $\approx 1.1$ for `sin`/`cos`). The run still
   prints it; the gate does not fail on it. Reversible, and it keeps the number visible.
2. **Normalise the oracle.** Build the excluded runners against the `libm` crate so both sides share
   one implementation. Honest for `tf_tree_math` (already `libm`) and sophus-rs (Rust, so
   patchable); impossible for the three container oracles (`PHASE1.md` §7), which are C++ on their
   platform's libm. Partial by construction.
3. **State it and accept the failures.** Leave the bar and record in §8 that a `Log`-shaped stratum
   may fail for D16. Cheapest, and it leaves a standing 8-stratum red that no change can clear,
   which is the condition that teaches a team to ignore a gate.

This record proposes **1**, with the penalty per transcendental taken from EA.13(d) and named in
§8 rather than hidden in the envelope.

## Consequences

- 8 of the 75 domination failures are reclassified, not fixed. `so3_log` stops being "the id where
  we never win" and becomes "the id where we are `libm`-bound".
- **The other 67 are not explained by this.** `Exp`'s transcendental is `sin_cos`, where `libm`'s
  $\nu$ is $1.07$–$1.09$ against glibc's $\le 1.00$ — a $\approx1.1\times$ penalty, not $\approx2\times$ —
  so `sen3_exp_n1`'s 36 failures (ratios up to 2.30) cannot be `libm` alone. The same swap
  experiment on `sin_cos` is the next measurement, and `0031`'s open question 2 and `0014`'s
  question 24 are the live hypotheses for the remainder.
- `0014` question 40 already names the `libm` difference **for the page's oracle columns**; it does
  not reach the bar. This record is the bar's half of that question and should be read with it.
- `just determinism` (Phase 6) is what checks the D16 claim this record spends. Until it runs, D16
  is a design intent that costs 0.74 u and is not yet verified across targets.

## Open questions

1. Is the penalty a per-function constant from EA.13(d), or measured per stratum? A constant is
   simple and loose; per stratum needs the oracle to answer on both `libm`s, which option 2's
   partial normalisation would give for two oracles and never for the containers.
2. Does `no-regress` need the same treatment? EA.13(d) says no — it compares `helicoid` with itself
   — but a `libm` version bump changes the bits and so the baseline, which is the same hazard
   arriving by another route.
3. Should `tf_tree_math` be the preferred oracle for exactly this reason? It is on the `libm` crate,
   so it is the only oracle whose comparison with `helicoid` is `libm`-neutral today, and it is
   bit-identical to us on all 30 `so3_log` strata. That makes it a weak oracle for `Log` (it cannot
   distinguish us from itself) and a fair one for everything else.
4. Is option 1 a `PHASE1.md` §8 amendment or an `EA.13` one? The number comes from the maths page;
   the verdict belongs to the spec.
