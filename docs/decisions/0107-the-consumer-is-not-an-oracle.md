# 0107: The consumer is not an oracle

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

[`0010`](./0010-seeded-from-tf-tree-math.md) decision 3 made `tf_tree_math` oracle #1:
`runners/tf_tree_math`, a workspace-excluded crate on `tf_tree` at commit `20bc5a05` (`v0.0.6`),
answering nine ids at binary64 (`so3_exp`, `so3_log`, `so3_act`, `so3_from_matrix`,
`so3_geodesic`, `sen3_exp_n1`, `sen3_log_n1`, `sen3_ad_n1`, `se3_geodesic`). It was the only
runner on the `libm` crate (`Backend::LibmCrate`), which the envelope used to class a domination
failure as "the program's" (`0036`, draft).

An oracle is worth its pairing only while it is independent of `helicoid`. `tf_tree_math` no longer
is. `helicoid`'s SO(3), SE(3) and geodesics are ports of it (`0010` decision 1), and `tf_tree` is
migrating onto `helicoid` ([`0041`](./0041-the-integration-is-an-adapter-at-tf-tree-math.md);
`PHASE4.md` §5.4: a re-export shim, then removal). Pinned, it measures the port against its source;
once the shim lands, it measures `helicoid` against itself. The owner decided on 2026-10-24: "it's
our own library and we start integrating helicoid inside."

**The last run against it.** `just oracles`, then `just envelope`, at `a270325` on one host, each
oracle's per-stratum maxima read from its `conformance` stdout (five significant digits):

- `tf_tree_math` scored 285 strata. On none was it the only oracle. On 141 another oracle tied its
  maximum, another beat it on 102, and on 42 it was the strictly best. Those 42 are 23
  `sen3_log_n1`, 6 `so3_from_matrix`, 4 `sen3_exp_n1`, 4 `so3_act`, 2 `so3_geodesic`,
  2 `se3_geodesic` and 1 `so3_log`, listed under *Consequences*.
- On those 42, `helicoid`'s maximum equals `tf_tree_math`'s on 34 and is below it on 7. It is
  above it on 1: `sen3_exp_n1/theta:pi-1e-7`, 1.5457 against 1.5338 `u`, 6 draws, unresolvable.
- `envelope` before and after the removal: 2797 strata paired with an oracle and 1768 with none
  (both runs). No resolved domination failure, 24 unresolved and 485 unresolvable excesses, and
  3 excepted strata (both runs). The one failure is main's standing "no baseline for `helicoid`".
  Two excesses change their oracle:

| Stratum | `helicoid` | Over, before | Over, after |
|---|---|---|---|
| `sen3_exp_n1/theta:pi-1e-7/f64` | 1.54567 | `tf_tree_math` 1.53385 | `sophus` 1.54542 |
| `so3_exp/theta:dense/f64` | 0.681802 | `tf_tree_math` 0.659456 | `sophus_rs` 0.659456 |

## Decision

1. **An oracle is independent of `helicoid`.** A consumer migrating onto `helicoid` is not an
   oracle (`PHASE1.md` §7, *Independence*).
2. **`runners/tf_tree_math` is deleted**, along with its `RUNNERS` entry in
   `xtask/src/conformance/oracle.rs`, `just oracle-tf-tree-math` and its step in `just oracles` and
   `just conformance-deep`, its `exclude` entry in the workspace, its steps in the disabled
   `envelope` CI job, its dependabot entry, and `deny.toml`'s `allow-git` for `tf_tree`, the only git
   source. No tracked file was the runner's output: `conformance/results/` is git-ignored. Corpus,
   baseline and exception table are unchanged.
3. **The `libm`-crate class goes with its only runner.** `Backend::LibmCrate`,
   `Index::beaten_by_libm_crate` and `Verdict::dominated_same_backend` are deleted, and so is the
   envelope's "the program's" count. A domination failure won by a host-`std` oracle is attributed
   by the one-variable twin alone (`0037`, draft), as it already was whenever no `libm`-crate oracle
   won. `Backend` keeps `HostStd`, which every runner declares and the attribution reads.
4. **Spec text.** `PHASE4.md` §5.2's precondition (dominate `tf_tree_math` before `tf_tree`
   migrates) is withdrawn. The ids behind §5.1 stay under `PHASE1.md` §8's domination bar
   against the remaining oracles, which answer all nine. `PHASE3.md` §9 drops the planned
   `tf_tree_math_differential`, which was never written. D15 and `0010` (amendment banner) say
   oracle #1 ended here. Past records that cite `tf_tree_math` as an oracle stay as written.

## Rationale

The bar is independent evidence that `helicoid` is no worse than other careful implementations of
the same mathematics. A copy of `helicoid`'s own source adds none of that. Of the 42 strata where
`tf_tree_math` was best, `helicoid` matched its maximum on 34: the port reproducing its source.

- **Keep the pin at `v0.0.6`.** This freezes a measurement of the port against what it was ported
  from. The port is done (`0010` plan steps 2–3; `PHASE4.md` §0.0), and the last run above records
  the result. Keeping it would also keep a git dependency on a repository that is migrating away
  from the code it pins.
- **Keep it until `tf_tree_math` is a shim.** It would then answer `helicoid`'s own bits and tie
  every stratum, which the bar reads as passing. The comparison would be empty, and nothing in the
  envelope would show it.
- **Keep `Backend::LibmCrate` for a future runner.** No candidate exists. A variant that nothing
  constructs is dead code, and the twin already measures directly what the class inferred: whether
  D16 explains a failure.

## Consequences

- No stratum loses its only oracle: Sophus and GTSAM answer all nine ids, manif all but
  `so3_from_matrix`, and sophus-rs the four `exp`/`log` ones. The paired count does not move (2797).
- **The bar loosens on 42 strata**, from `tf_tree_math`'s maximum to the next oracle's. The
  maxima are binary64 `u`, five significant digits, sorted by the ratio of the new bar to the old:

| Stratum | `tf_tree_math` | Next oracle | Ratio | `helicoid` |
|---|---|---|---|---|
| `sen3_log_n1/theta:1e-4` | 0.39677 | manif 1.1768 | 2.97 | 0.39677 |
| `so3_from_matrix/q:w0` | 0.59283 | sophus 1.6589 | 2.80 | 0.59283 |
| `sen3_log_n1/rho:1e-3/theta=1e-4` | 0.63753 | gtsam 1.7262 | 2.71 | 0.63753 |
| `sen3_log_n1/rho:1e3/theta=1e-1` | 0.70524 | sophus_rs 1.4197 | 2.01 | 0.67403 |
| `se3_geodesic/geo:generic` | 2.5019 | gtsam 4.4303 | 1.77 | 2.0574 |
| `sen3_log_n1/rho:1e0/theta=1e-8` | 0.47205 | sophus_rs 0.76231 | 1.61 | 0.47205 |
| `so3_from_matrix/theta:pi-1e-2` | 1.3490 | sophus 2.1610 | 1.60 | 1.3490 |
| `sen3_log_n1/rho:1e3/theta=1e-8` | 0.49054 | sophus_rs 0.77555 | 1.58 | 0.49054 |
| `so3_from_matrix/theta:pi-1e-9` | 0.95415 | sophus 1.4964 | 1.57 | 0.95415 |
| `sen3_log_n1/theta:1e-6` | 0.72557 | gtsam 1.1057 | 1.52 | 0.72557 |
| `so3_geodesic/geo:near-pi` | 1.6417 | gtsam 2.3870 | 1.45 | 1.6417 |
| `so3_from_matrix/theta:pi-1e-11` | 0.86141 | sophus 1.2346 | 1.43 | 0.86141 |
| `sen3_log_n1/rho:1e3/theta=1e-4` | 0.57731 | sophus 0.80640 | 1.40 | 0.57731 |
| `sen3_log_n1/theta:1e-7` | 0.88027 | sophus_rs 1.2115 | 1.38 | 0.88027 |
| `sen3_exp_n1/theta:1e-4` | 0.60686 | gtsam 0.82056 | 1.35 | 0.60686 |
| `sen3_log_n1/theta:1e-3` | 0.82505 | manif 1.0928 | 1.32 | 0.82505 |
| `sen3_log_n1/theta:pi-1e-3` | 0.92156 | sophus 1.1915 | 1.29 | 0.92156 |
| `sen3_log_n1/rho:1e-3/theta=1e-8` | 0.70429 | sophus_rs 0.88297 | 1.25 | 0.70429 |
| `sen3_log_n1/rho:1e4/theta=1e-4` | 0.83351 | gtsam 1.0409 | 1.25 | 0.83351 |
| `sen3_log_n1/rho:1e0/theta=1` | 1.1628 | sophus 1.4336 | 1.23 | 0.83436 |
| `so3_act/theta:dense` | 1.8204 | gtsam 2.2411 | 1.23 | 1.8204 |
| `sen3_log_n1/theta:1e-8` | 0.64728 | sophus_rs 0.78210 | 1.21 | 0.64728 |
| `sen3_log_n1/theta:1e-9` | 0.48458 | sophus_rs 0.58360 | 1.20 | 0.48458 |
| `se3_geodesic/geo:near-pi` | 3.2528 | sophus 3.8301 | 1.18 | 2.7206 |
| `sen3_log_n1/rho:1e4/theta=1e-8` | 0.79681 | sophus_rs 0.93551 | 1.17 | 0.79681 |
| `sen3_exp_n1/rho:1e-6/theta=1` | 1.0671 | sophus_rs 1.2454 | 1.17 | 0.76842 |
| `sen3_log_n1/theta:1e-10` | 0.51203 | sophus_rs 0.59247 | 1.16 | 0.51203 |
| `so3_log/q:nonunit` | 1.8136 | sophus_rs 2.0900 | 1.15 | 1.8136 |
| `sen3_log_n1/theta:1e-2` | 0.92775 | sophus_rs 1.0565 | 1.14 | 0.92775 |
| `so3_act/theta:pi-1e-3` | 4.1980 | gtsam 4.7511 | 1.13 | 4.1980 |
| `so3_act/theta:pi-1e-1` | 3.7450 | gtsam 4.1499 | 1.11 | 3.7450 |
| `sen3_log_n1/rho:1e4/theta=pi-1e-6` | 2.1169 | sophus 2.3149 | 1.09 | 2.1169 |
| `sen3_log_n1/theta:pi-1e-10` | 1.5924 | sophus_rs 1.6753 | 1.05 | 1.5924 |
| `so3_from_matrix/theta:pi-1e-4` | 1.5186 | sophus 1.5892 | 1.05 | 1.5186 |
| `sen3_exp_n1/rho:1e-6/theta=1e-4` | 0.82325 | gtsam 0.85509 | 1.04 | 0.82325 |
| `so3_geodesic/geo:generic` | 1.8341 | sophus 1.8873 | 1.03 | 1.7382 |
| `so3_from_matrix/theta:pi-1e-7` | 1.4393 | sophus 1.4788 | 1.03 | 1.4393 |
| `sen3_log_n1/theta:1e-5` | 0.66551 | sophus_rs 0.68273 | 1.03 | 0.66551 |
| `so3_act/theta:pi-1e-12` | 3.1613 | gtsam 3.2056 | 1.01 | 3.1613 |
| `sen3_log_n1/rho:1e-3/theta=1e-1` | 1.1493 | sophus_rs 1.1637 | 1.01 | 0.73879 |
| `sen3_exp_n1/theta:pi-1e-7` | 1.5338 | sophus 1.5454 | 1.01 | 1.5457 |
| `sen3_log_n1/theta:pi-1e-9` | 0.82270 | sophus 0.82450 | 1.00 | 0.82270 |

  No-regress is untouched: once a baseline is blessed it holds each of these maxima exactly, so a
  regression on any of them fails whatever the oracles read.
- The 141 tied strata keep their bar (at the printed precision): the tying oracle holds it.
- Every runner now reads its transcendentals from a host's system library. A domination failure is
  either D16's cost or not by the twin's measurement (`0037`, draft), never by an oracle's
  backend.
- `tf_tree`'s migration evidence is `tf_tree`'s to keep (`0010` decision 4): its own bench gate and
  `just interp-accuracy` (`PHASE4.md` §5.3).

## Implementation plan

1. This record's PR: the deletions and spec edits of decisions 2–4 — verified by `just lint`,
   `just doc`, `just test`, and `just envelope` reading 2797 paired strata with only main's
   standing no-baseline failure.
