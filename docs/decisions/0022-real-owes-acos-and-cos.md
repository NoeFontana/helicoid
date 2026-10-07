# 0022: `Real` owes `acos` and `cos`

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** `pi` as a literal has landed; steps 1 and 2 landed with [`0053`](./0053-acos-is-better-the-roots-are-not-necessarily.md).

> **Amended 2026-10-07 by [`0053`](./0053-acos-is-better-the-roots-are-not-necessarily.md), the signed edit the index's rules allow.** Both steps are
> implemented as written and the *Decision* stands, on reason 1 and the function-level half of
> reason 2. Two of the three reasons do not survive re-measurement:
>
> - **Reason 2 is narrowed.** `libm::acos` is better than the `atan2` form everywhere measured —
>   independently reproduced, 0.852 `u` against 1.659 on the grid and 0.489 against 1.148 as
>   `x -> 1-` — but a more accurate `acos` is **not** a more accurate root. On the three `GOLDEN`
>   rows that move, scored against `mpmath`, two roots get closer to the truth and one gets further,
>   0.522 `u` to 7.478; `eig3`'s worst `value / bound` over its fixture moves 0.150 to 0.206. "More
>   accurate, by a little, everywhere measured" is true of the function and false of both routines.
> - **Reason 3 is withdrawn.** "Back toward omnisac" is not verifiable in this repository: omnisac
>   is not a dependency, nothing committed records its per-row output, and `bits_are_omnisacs`
>   *passed* before the change, so the port already reproduced all 36 stored triples.
>
> Reason 1 is confirmed: 2.63x against this record's stated 2.61x for `acos`, and `libm::cos` at
> 7.7 ns against `sin_cos().1`'s 9.6. `0053` also records what this record assumed without
> establishing — that `cos` is bit-identical to `sin_cos().1`, which holds by construction and
> exhaustively over all `2^32` binary32 — and what neither record anticipated: **no corpus id
> reaches `solve_cubic` or `eig3`**, so `0006`'s bars were never the judge here.

> **Retitled and rescoped.** This record was *"`solve_cubic`'s `acos` and the limits it inherits"*
> and carried two open questions. The first is answered here, against the measurements its first
> draft did not have. The second — which of the four numerical limits the port inherits from omnisac
> get a fix — is a different subject with a larger plan and is now
> [`0031`](./0031-what-the-cubic-port-inherits-from-omnisac.md) (draft).

## Context

`0017` step 2 ports omnisac's `solve_cubic` "as is". omnisac's three-real-root arm calls
`f64::acos`, `f64::cos` and `core::f64::consts::PI`. `Real` has `sqrt`, `cbrt`, `sin_cos`, `atan2`,
`abs` and `copysign` only, and no constant but `lit`, so the port spelled all three in terms of what
it had:

```rust
pi  = atan2(+0, -1)
acos x = atan2(sqrt((1 - x)(1 + x)), x)
cos phi = phi.sin_cos().1
```

The first draft took these as obviously cheap and the factored `acos` as *more* accurate than
`libm`'s near a double root, on the reasoning that `sqrt((1 - x)(1 + x))` beats `sqrt(1 - x²)` where
`x -> ±1`. That reasoning is sound about the **subexpression** and says nothing about `libm::acos`,
which does not compute it that way. Nobody had compared the two.

## Measurement

Advisory, as `0019`'s was: no document prescribes a benchmark for `helicoid-linalg`
(`PHASE1.md` §9 scopes `criterion` to `crates/helicoid/benches` from Phase 3). Scratch crate
outside the repository, not committed, `--release`, `lto = "fat"`, `codegen-units = 1`, no
`target-cpu`, AMD EPYC-Milan (8 vCPU, KVM guest), `rustc` 1.98.1, `taskset`-pinned, min of 41 over
4096 pre-generated arguments, `black_box` both ends, two cores. Accuracy against `mpmath` at 60
digits, worst relative error in units of `2^-53`.

**Accuracy — the factored form is not better anywhere measured, including where it was supposed to
be:**

| region | `atan2` form | `libm::acos` |
|---|---|---|
| grid over `[-1, 1]`, 4095 points | 1.31 u | **0.98 u** |
| `x -> 1⁻` (the double root), `1 - 2^-k`, `k = 1..52` | 1.28 u | **0.92 u** |
| `x -> -1⁺`, `-1 + 2^-k` | 0.92 u | 0.92 u |
| `x` near 0, `±2^-k` | 0.99 u | **0.92 u** |

The two agree to within **1 ulp** over the same points, and both are bit-equal at `x = 0, ±1`.

**Cost:**

| call | ns | note |
|---|---|---|
| `libm::acos` | 4.37 | |
| the `atan2` form | 11.42 | **2.61×** the call it replaces, per solve |
| `libm::sincos` | 5.66 | what the arm calls, three times |
| `libm::cos` | 4.19 | **26%** cheaper; 4.4 ns per solve over three phases |

`pi` is the third of the three and is already settled: it is a compile-time constant, so it is now a
per-precision literal, bit-equal to the `atan2` spelling by
`pi_and_acos_are_within_their_ulps`. That landed as a behaviour-preserving change and needed no
decision.

## Decision

**`Real` gains two required methods, `acos` and `cos`**, routed through `libm` like every other
transcendental (D16), and `cubic.rs` loses its private `acos` and its three discarded sines.

Three reasons, and they do not conflict:

1. **It is faster.** 7.0 ns per solve from `acos`, 4.4 ns from `cos`: ~11 ns off an arm that is
   dominated by four or five transcendentals.
2. **It is more accurate**, by a little, everywhere measured — including the near-±1 region the
   factored form was chosen for.
3. **It restores bit-identity with the reference implementation.** The trigonometric arm is the one
   place this port is *not* bit-equal to omnisac: the `GOLDEN` header records up to 5 `u` of the
   largest root over the 4e4 rows the golden rows are drawn from and 11 `u` over 5e6 random
   polynomials. Routing `acos` and `cos` to the same `libm` functions omnisac calls removes that
   deviation, and `bits_are_omnisacs` becomes a statement about every arm. For a workspace whose
   verification strategy *is* the reference twin (D6), buying that back for one trait method is
   cheap.

4. **There is more than one caller.** `eig3` (`0017` step 3, on `linalg/eig3`) imports
   `cubic::{acos, pi}` and computes Smith's eigenvalues as `q + 2 p cos((acos(r) + 2 pi k) / 3)`,
   `k = 0, 1, 2` — the same shape as the trigonometric arm, one `acos` and three cosines per call.
   The first draft said "public surface for one caller" while `eig3` was still unwritten, so that
   was fair at the time; `eig3`'s own PR then **amended decision 1 of this record** to make the two
   helpers `pub(crate)` "(`eig3`, `0023` (draft), calls them too)" — recording the second caller
   without revisiting the conclusion it undermines. This record revisits it. `eig3` is also the
   consumer that cares *most* about accuracy near `±1`: its own rustdoc
   records that `r` sits near `±1` at a double eigenvalue, "where `acos` has an infinite slope", so
   an error `u` in `r` becomes `sqrt(u)` in the angle. That is precisely the region where the
   measurement above puts `libm::acos` ahead (0.92 u against 1.28 u).

`cbrt` is the precedent: `0017` added a required method to `Real` for one caller, on the same
reasoning. The cost here is the same — one more method on the four in-workspace impls (`f64`, `f32`,
`Dual`, the test lanes and `xtask`'s seeded lane), and an out-of-tree `impl Real` breaks, which
`0.0.x` permits — and it is paid for twice over.

## Rationale

The alternative is to keep the private spellings. Their argument was "public surface for one
caller", and it fails twice over: there are two callers, as `eig3`'s PR itself records, and the
argument holds only while the form is free and at least as good. It is neither — 2.6× the cost,
marginally worse, and it forfeits the bit-identity that makes the port checkable. The lesson is
narrower than the decision: a cost comparison against the function being replaced was never run,
and three documents reasoned about the trade without one.

Keeping `sin_cos` and discarding the sine is the smaller of the two: 4.4 ns per solve, and `cos`
is wanted by nothing else today. It goes in the same step because the two methods are one
`Real`-surface change, one set of `Dual` rules and one set of tests, and splitting them would pay
that overhead twice.

## Consequences

- `Dual` rules: `acos` is `-d / sqrt((1 - v)(1 + v))`, singular at `v = ±1` exactly as the current
  form is (`sqrt` at 0), so the `# Domain` row is unchanged in substance; `cos` is `-d sin v`, and
  the value lane must stay bit-equal to `sin_cos().1` so `dual_value_is_plain_value` holds.
- `PHASE2.md` §2's `libm` list, §3's `Dual` rules, §8's checks and §6's "`Real` has no `acos`"
  sentence all move, as does `API.md` §2 and `NUMERICS.md` §12's rows. `0017` decision 1's list of
  required methods gains two.
- **`eig3` moves with it.** Its eigenvalues shift by up to about 1 ulp of `acos`'s output, amplified
  by `sqrt(u)` near a double eigenvalue, so its golden rows and its fitted error bounds
  (`eig3_tests`) must be re-measured in the same step. `eig3` sits directly on top of `solve_cubic`
  in the review stack, so if it merges first this is a follow-up against it; the work is the same
  either way, it is only cheaper to pay once.
- `GOLDEN`'s header loses the 5 `u` / 11 `u` deviation note. The trigonometric arm's values move by
  up to 11 `u` **back toward** omnisac, so the golden rows for three-real-root polynomials must be
  regenerated from the new code and the mpmath strata re-measured; a polynomial sitting between the
  old and new values could in principle change root count, which `PHASE2.md` §9's second clause
  already covers.
- This does **not** touch any arm's formula. The four inherited numerical limits are `0031` (draft)
  and are independent of this record: `0031` can land before or after it.

## Implementation plan

1. `Real::acos` and `Real::cos`: the trait methods, the `f64`/`f32` `libm` routing, the `Dual`
   rules, every in-workspace impl, and the doc edits above. Verification to the bar `cbrt` set
   (`0017` step 1): `libm` routing bit for bit at both precisions, the IEEE corners, the derivative
   against mpmath at exact rows, second order through nesting, and the mutations of each rule.
2. `cubic.rs` drops the private `acos` and uses `Real::cos`; regenerate the `GOLDEN` three-root rows
   and re-measure the mpmath strata; update the header's deviation note and `PHASE2.md` §6.
