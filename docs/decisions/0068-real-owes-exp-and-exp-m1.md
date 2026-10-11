# 0068: `Real` owes `exp` and `exp_m1`

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #134

## Context

Every Sim(3) routine needs `e^σ` and `expm1 σ` (`docs/maths/sim3.md` SM.4, SM.5(b)): the scale
`s = e^σ`, `v₀ = expm1(σ)/σ`, the functions `g_ℓ(σ)`, and the scaled arrangement's `e^{-σ}`,
`expm1(-σ)` and `e^{-σ/2}` (SM.7(g)). `Real` (`PHASE2.md` §2) has neither, and
`PHASE5.md` §3 waits on it. `docs/maths/index.md` lists this as an open item and names the shape of
its answer: a record like [`0022`](./0022-real-owes-acos-and-cos.md), with a `Dual` rule for each, a
`libm` route (D16), `real_*` strata (`PHASE2.md` §8) and `API.md` §6's answers. It also notes that
`PHASE2.md` §2's listing of `Real` was never updated for `0022`'s `cos` and `acos`.

`exp(σ) - 1` is not a substitute for `expm1`: it cancels for `|σ| ≪ 1`, exactly where `v₀` and
`g_ℓ` are evaluated near the joint limit, and at `σ = 10⁻¹⁰` keeps six digits of sixteen.

## Decision

1. **`Real` gains two required methods**, `exp` and `exp_m1`, routed through `libm::{exp, expm1}`
   at `f64` and `libm::{expf, expm1f}` at `f32` (D16), in `float.rs`'s one macro.
2. **`Dual`**: `exp` carries `d e^v` with the value's own `e^v`; `exp_m1` carries `d e^v` from
   `exp`, **not** from `1 + exp_m1(v)`, which rounds to `0` below `v = -37.4` (`f64`) and
   `-17.3` (`f32`) while `e^v` is still normal: that sum would lose every digit of the derivative.
3. **Domain** (`NUMERICS.md` §12): every argument. Above `ln MAX` (`709.78`, `88.72`) the value
   and the seeded lanes are `+∞`; `exp` is subnormal below `-708.40` (`-87.34`) and `+0` below
   `-745.13` (`-103.97`); `exp_m1` is `-1` exactly below `-37.4` (`-17.3`). Nothing is asserted.
4. **Corpus ids `real_exp` and `real_expm1`** (`PHASE1.md` §4.3), the value and the `Dual`
   derivative scored as the other `real_*` ids, 64 records per stratum at both precisions:
   `x:tiny` (`|x|` log-uniform in `[1e-300, 1e-8]`, `[1e-30, 1e-4]` at `f32`, both signs),
   `x:unit` (`[-1, 1]`), `x:moderate` (`|x| ∈ [1, 30]`, `[1, 15]`, both signs), `x:large` (to
   `ln MAX`); for `exp` also `x:negative` (a normal `e^x` far below 1) and `x:subnormal` (a
   subnormal `e^x`, judged by its quantization); for `expm1` also `x:saturated` (`expm1 = -1`,
   its derivative normal). References `mp.exp`, `mp.expm1`, cross-checked by `v e^{-x} = 1`,
   `d = v`, `d - v = 1` and `mp.diff`.
5. **`PHASE2.md` §2** lists `cos`, `acos`, `exp` and `exp_m1` beside the rest.

## Measured

`max_u` per stratum, the worse of the value and the `Dual` derivative (`just conformance` at both
precisions; no non-finite output):

| id | precision | `tiny` | `unit` | `moderate` | `large` | `negative` | `subnormal` | `saturated` |
|---|---|---|---|---|---|---|---|---|
| `real_exp` | `f64` | 0.051 | 0.925 | 0.993 | 0.900 | 0.912 | 0.985 | — |
| `real_exp` | `f32` | 0.935 | 0.885 | 0.984 | 0.978 | 1.003 | 1.218 | — |
| `real_expm1` | `f64` | 0.311 | 0.990 | 0.993 | 0.900 | — | — | 0.961 |
| `real_expm1` | `f32` | 0.935 | 1.183 | 0.984 | 0.978 | — | — | 0.879 |

Every other corpus file regenerates byte-identical.

## Rationale

- **Two methods, not `exp` alone.** `expm1` from `exp` cancels where Sim(3) needs it most; `exp`
  from `expm1` adds a rounding where it is not needed. `libm` has both.
- **The derivative from `exp`** costs one more call per `Dual` `exp_m1` and keeps the derivative's
  relative accuracy everywhere the value has it; `1 + exp_m1` does not, past `-37.4`.

## Consequences

- Every `impl Real` gains two methods: `f64`, `f32`, `Dual`, the test lanes (`L2`, `Lane`,
  `Counted`) and `xtask`'s `Host` and seeded `Lane`. An out-of-tree `impl Real` breaks, which
  `0.0.x` permits.
- `real_exp` and `real_expm1` are required by coverage.
- `docs/maths/index.md`'s open item (iv) on Sim(3) is answered.

## Implementation plan

1. One PR (#134): the methods, the `Dual` rules and tests (`dual_exp_and_exp_m1_match_mpmath`,
   `nesting_gives_second_derivatives_of_exp`), the two corpus ids, the docs. Verified by
   `just lint`, `just test`, `just corpus-check`, `just conformance` at both precisions.

## Open questions

None.
