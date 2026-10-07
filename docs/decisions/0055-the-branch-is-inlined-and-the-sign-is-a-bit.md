# 0055: the branch is inlined, and the sign is a bit

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record and the two changes land together

## Context

[`0054`](./0054-the-screw-twin-takes-two-arms-in-the-world-frame.md) found two things by disassembly
while tuning the SE(3) screw twin, and listed both as *Further work*:

- **`Real::branch` has no `#[inline]`.** It is a provided method, generic over its two arm closures,
  so LLVM may inline it — and in the `groups` bench binary it did not: each `branch` was a separate
  function receiving its arms by value. The twin with two sequential `branch` calls read 1.18× the
  provided body until it was restructured to one. Every coefficient in `helicoid::coeffs` reaches
  its caller through `branch` (`grouped`'s two levels, `log_ratio`'s, `jr_inv_coeff`'s).
- **`Real::copysign` and `Real::abs` are calls to `libm::copysign` and `libm::fabs`**, out of line
  across the crate boundary on every `Log`'s sign flip, where each is one bit operation.

## Decision

1. `Real::branch` is `#[inline]`.
2. `Real::abs` and `Real::copysign` for `f64`/`f32` are `core`'s `f64::abs`/`copysign` (and the
   `f32` pair). Both are exact sign-bit operations: the same bits as `libm`'s for every input, NaN
   included. D16 routes every **transcendental** through `libm`; neither is one.
3. **Not** `#[inline]` on the coefficient kernel's entry points (`exp_coeffs`, `log_ratio`,
   `grouped`, …): measured below, it is a loss.

## Rationale

**Screen.** Five binaries of each bench target, each built in this tree — the base, A (decision 1),
B (decision 2), A+B, and C (`#[inline]` on the kernel) — run once over every row, core 5 pinned, in
the order base, A, B, A+B, C, base. Ratios of criterion medians to the base; the base's own
re-run is the drift, `p5`–`p95` 0.994–1.013 (`groups`) and 0.985–1.017 (`coeffs`):

| candidate | `groups` (188 rows) geomean / max | `coeffs` (60 rows) geomean / max |
|---|---|---|
| A, `branch` inlined | 0.923 / 1.017 | 0.583 / 1.005 |
| B, `core` sign ops | 0.995 / 1.044 | 1.003 / 1.041 |
| **A+B** | **0.918** / 1.035 | **0.554** / 1.008 |
| C, kernel inlined | 1.019 / **1.303** | 1.124 / **3.000** |

C makes `jr_inv_coeff` 2–3× slower: forcing the whole kernel into its callers defeats the
branch layout `0047` item 7 measured. B alone is small and confined where it should be — its best
rows are the ones that flip a sign (`so3/log` 0.89–0.93, `so3/geodesic` 0.89 at near-identity),
and its worst are rows that never call `copysign` (`se23/jr_inv`), so layout.

**Gate.** `cargo xtask bench-gate --against <base>` over every row of both targets, A+B built in
this tree, core 5 pinned: **248 benchmarks, none above its own concurrent floor in all 6 pairs**
(median floor 0.6 %). Point ratios:

| `coeffs` | | `groups` | |
|---|---|---|---|
| `q_coeffs` | 0.378 | `so3/exp` | 0.568 |
| `exp_coeffs` | 0.418 | `so3/jr` | 0.744 |
| `jr_coeffs` | 0.433 | `so3/log` | 0.768 |
| `log_ratio` | 0.742 | `se3/exp` | 0.780 |
| `jr_inv_coeff` | 1.000 | `so3/geodesic` | 0.818 |
| **all 60** | **0.551** | **all 188** | **0.915** |

(family geomeans over strata and precisions). 77 of the 188 `groups` rows are more than 5 %
faster. The largest point above 1 is 1.072, on an `adjoint` row whose own A/A floor is 44 %:
unresolvable, not a regression.

**Bits.** Inlining cannot change a result — Rust contracts no `mul_add` (D16) — and the sign
operations are exact. Checked: `cargo xtask conformance --subject helicoid` writes all 1566 rows
byte-identical to the base's, the git revision column aside.

**The consumer.** In a scratch copy of `tf_tree` with `ScLerp` delegating to the `0054` twin,
`lookup/depth3/sclerp` reads 254.8 / 256.4 / 256.7 ns with this record against 256.6 / 255.8 /
255.5 without, interleaved, and `screw_pow`'s 274.9 / 274.9 / 275.1: unchanged, because `0054`'s
restructure already left the twin inlined there. This record removes the trap for every caller
that has not been restructured around it.

## Consequences

- A routine may take as many `S::branch` calls as its arithmetic needs; `0054`'s single-branch
  restructure stays because it is also the clearer code, not because it is required.
- `libm::copysign` and `libm::fabs` are no longer reached; `cargo xtask lint`'s `kernel` check
  (`libm::` only in D16's kernel) still holds, `float.rs` being that kernel.
- Benchmarks recorded in earlier records (`0050`, `0051`, `0052`, `0054`) are ratios against
  baselines built before this one; they stay true of what they compared and are not re-run.

## Implementation plan

1. Decisions 1 and 2 — verified by the gate above, the byte-identical conformance rows, `just lint`,
   `just test`, `just doc`, `just msrv`, `just no-std`, `just wasm`.

## Open questions

None.

## Further work

1. **`Mask::decide`'s lane impls**, when a lane `Real` exists (`0013`): the same question, which
   only a lane bench can answer.
