# 0059: the short arm is written out

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** the PR that lands this record

## Context

`tf_tree`'s first adapter wave (its record `0063`) delegates `tf_tree_math::slerp` to
`SO3::geodesic`. Its own `lookup/depth3/lerpslerp` bench read **207.9 ns native against 238.9 ns
delegated**, +15 %, over four interleaved runs. Every other `lookup` row was at parity. `0063`'s
decision 5 keeps the native body of a function that regresses on latency, so `slerp` would have
stayed native, giving up `0058`'s near-pair drift accuracy.

The fixture's keyframes are 1 to 20 ms apart, and their relative rotations sit below `r`'s second
switch (`s = tan²α` of order `1e-9` to `1e-6`). So the bench runs `0051`'s arm below the switch,
the provided body `q₀ Exp(t r v)`, and never the blend.

### Measured

**Protocol.** A scratch copy of `tf_tree` at its Wave 1 branch has `helicoid` as a path
dependency. One bench binary is built per variant. The binaries run interleaved, criterion's median
per run, pinned with `taskset -c 3`, one round being every binary once. The host is the one `0055`
used, and figures are averages over rounds.

**Where the time goes.** These are diagnostic builds, not candidates:

| Build | `lerpslerp` |
|---|---|
| native `slerp` | 208.3 ns |
| `helicoid` 0.0.1 | 239.6 ns |
| feature on, adapter calling the native body | 209.1 ns |
| `r` replaced by the constant `2` | 193.6 ns |
| `Exp`'s coefficients replaced by `(½, 1)` | 200.6 ns |
| the arm returns `q₀` once dispatched | 111.2 ns |

The whole gap is the geodesic. Removing either series saves more than the gap. The arm's
instructions are few, and the arm is latency-bound. Each step waits on the one before:
`d = q₀* q₁` → `n²` → `s = n²/w²` (a division) → `r`'s Horner → `φ = t r v` → `θ² = ‖φ‖²` →
`Exp`'s Horner → `q₀ Exp(φ)`. The depth-3 lookup runs one such chain per dynamic edge.

**Disassembly** of the delegated `LerpSlerp::eval` showed two more costs on that chain:

- `log_ratio` stays a call. The geodesic spills seven registers around it, and `log_ratio` then
  re-decides the arm that `log_ratio_takes_short_arm` has just decided.
- Every Horner starts its fold at zero, so the first step is `t + z · 0`. `z · 0` cannot be folded,
  because it is NaN at `z = ∞`.

## Decision

Four changes, each exact in `R`. The first three give the same bits.

1. **`horner` starts its fold at the last term.** Every selected argument is finite, and for those
   `t + z · 0` is `t` exactly. The first step is therefore the same bits, one multiply-add sooner,
   in every coefficient in `helicoid::coeffs`.
2. **`log_ratio_short(n2, w, short)`** is `log_ratio`'s short arm alone, at the safe argument
   `short` selects. The geodesic's arm below the switch, and the screw twin's `t = 1` translation,
   call it with the mask they already hold. Its expression is `log_ratio`'s short arm. The two
   predicates, `n² < T w²` and `n²/w² < T`, can part by one ulp exactly at `T`. There this arm is
   still the series of its own switch, and no corpus record lands there.
   It is `#[inline]`: one Horner, not the branch layout of the kernel's entry points that `0055`
   decision 3 keeps out of line.
3. **`n²` is formed from `d`, not from the flipped copy.** Because `flip` is `±1`, the bits are the
   same and the flip's multiply leaves the chain.
4. **The arm below the switch is written out** by two identities. In `geodesic_short`:
   - `θ² = ‖t r v‖² = r² (t² n²)`. `t² n²` waits on nothing, so `Exp`'s Horner starts two
     multiplies after `r` rather than five.
   - `q₀ (c, k t r v) = c q₀ + k (t r (q₀ (0, v)))`, by bilinearity. The product `q₀ (0, v)` waits
     on nothing either, so the end of the chain is one multiply-add rather than a Hamilton product.

   Rounding aside, this is the provided body. `θ` and the axis are `Exp`'s, so the arm stays
   scale-free on a carried quaternion (`0058`). `NUMERICS.md` §10 states the identities.
   `measure_geodesic` gains the arm as route F, an independent transcription through `__sweep`.
   Its tests now tie the shipped arm to F bit for bit, where they used to tie it to the provided
   body.

## Rationale

**The consumer.** The same protocol, three rounds per step, each step adding to the last:

| Build | `lerpslerp` |
|---|---|
| native | 208.3 ns |
| `helicoid` 0.0.1 | 239.6 ns |
| + decision 1 | 235.4 ns |
| + decision 2 | 234.6 ns |
| + decision 4 | 215.9 ns |
| + decision 3, and `t r` applied to `q₀ (0, v)` before `k` | **210.5 ns** |

The final A/B, five interleaved rounds: native **209.3 ns**, `helicoid` 0.0.1 **239.1 ns**, this
record **212.0 ns**. That is +1.3 % over native, down from +14.2 %. The `sclerp` rows are
unchanged within ±1 %.

**Accuracy.** `cargo xtask conformance --subject helicoid` writes 1 673 rows. 1 671 of them are
byte-identical to `main`'s, the git revision aside. The two that move improve:

| Row | Before | After | Oracle #1 |
|---|---|---|---|
| `so3_geodesic` `geo:consecutive` | 1.5721 `u` | **0.9929 `u`** | 2.187 (`slerp`) |
| `se3_geodesic` `geo:consecutive` | 1.5721 `u` | **0.9929 `u`** | 2.336 (`ScLerp`) |

`geo:generic` and `geo:near-pi` take the blend and do not move. Under `0058`'s carried draws, the
arm reads 5.363 `u` at every `η`, against 5.612 before. The blend's tilt ratio is unchanged at
1.0027. Every law and bound in `just test` holds as written.

**`helicoid`'s own gate.** `cargo xtask bench-gate --against <main>`, core 5 pinned, with
`main`'s binaries built in a worktree. On `coeffs`, **60 benchmarks, none above its own
concurrent floor in all 6 pairs** (median floor 0.68 %). The median ratio per row, as family
geomeans:

| `coeffs` | `f64` | `f32` |
|---|---|---|
| `exp_coeffs` | 0.956 | 0.938 |
| `jr_coeffs` | 0.979 | 1.000 |
| `jr_inv_coeff` | 0.946 | 0.899 |
| `log_ratio` | 0.936 | 0.904 |
| `q_coeffs` | 0.951 | 0.956 |
| **all 60** | **0.946** | |

The series-arm rows (`near-identity`, `small`) carry the gain, from 0.72 to 0.94. The exact-arm
rows read 1.00 ± 0.01, as decision 1 predicts: their Horner runs only on the members below
their own switch. The largest point is 1.008.

On `groups`, **186 of 188 benchmarks pass** (median floor 0.67 %, geomean 0.982). The geodesic
rows:

| `geodesic` | `f64` | `f32` |
|---|---|---|
| `so3` `near-identity-7.5e-8` | 0.633 | 0.629 |
| `se3` `near-identity-7.5e-8` | 0.915 | 0.875 |
| `se3` `consecutive-1e-3` | 0.910 | 0.879 |
| `so3` `generic-1` / `near-pi` | 0.978 / 0.986 | 0.966 / 0.985 |
| `se3` `generic-1` / `near-pi` | 0.986 / 0.983 | **1.034 / 1.028** |

**The two `f32` SE(3) rows above the blend's switch fail the gate, and this record accepts them.**
They run the long arm, which none of the four decisions touches. Its instructions are the same in
both binaries, up to spill slots. Each alternative below was gated on the 14 geodesic rows:

| Variant | `f32` `se3` `generic` / `near-pi` | Other rows past their floor |
|---|---|---|
| this record | 1.032 / 1.028 | none |
| every function aligned to 64 bytes (`-align-all-functions=6`), both sides | 1.031 / 1.026 | none |
| without decision 1 | 1.036 / 1.027 | none |
| the screw twin's `t = 1` call back on `log_ratio` | 1.033 / 1.028 | none |
| `θ² = (t r)² n²` | 1.029 / 1.023 | none |
| `log_ratio_short` not `#[inline]` | 1.032 / 1.029 | none |
| without decision 3 | 1.045 / 1.040 | `f32` `so3` long arm, 1.041 / 1.052 |
| without decision 4 | 0.962 / 0.972 | `f64` `se3` long arm 1.023 / 1.021, `f32` `so3` 1.008 / 1.012 |

So layout is not the cause, and neither is any one decision. The long arm shares one function with
the short arm, and LLVM allocates registers for both together. Any change to the short arm moves
the long-arm rows by up to 3 % in one direction or the other, and `main`'s figures are one draw of
that. Steering a register allocator from source would be fitted to one LLVM version and would be
code no reader could maintain. The short arm's gain is the depth of its dependency chain, which
holds under any allocation. The cost is 3 % on two `f32` rows, against 13–37 % on the rows below
the switch and parity for the `f64` consumer this was done for.

**Alternatives.**

- **Keep `slerp` native in `tf_tree`** (`0063` decision 5). This would cost the consumer
  `0058`'s near-pair accuracy (3 984 `u` against 6.5 `u` at the band's edge) and a stratum
  dominated by more than 2×.
- **`#[inline]` on the kernel.** `0055` decision 3 measured it as a loss across the kernel, and
  here the call was 1 ns of the 31. The chain is the cost, not the call.
- **Estrin's scheme, or a third, shorter series tier.** Either would halve the two Horners' depth,
  but each changes CO.9's evaluation form or the sweep's admission rule for every coefficient. That
  is a record of its own (*Further work*), and the consumer is at parity without it.

## Consequences

- `SO3::geodesic` and `SE3::geodesic` below `r`'s second switch are the written-out arm, not
  `reference::geodesic`'s expression. The `twin` leg of `laws::geodesic` already compares two
  different expressions there.
- `PHASE4.md` §0.0's `se3_geodesic` row reads 0.993 at `geo:consecutive`.
- `tf_tree`'s `0063` can keep `slerp` delegated, once a release carries this record.
- `groups/f32/se3/geodesic/{generic-1, near-pi}` read about 1.03× `0.0.1`'s, for the reason
  above. A later `bench-gate --against` this record's binary measures from the new figures.

## Implementation plan

1. Decisions 1–4, `NUMERICS.md` §10, `measure_geodesic`'s route F and its tests. Verified by the
   conformance diff above, `just test`, `just lint`, `just doc`, and `bench-gate --against` on
   `coeffs` (passes) and `groups` (passes apart from the two rows accepted above).

## Open questions

None.

## Further work

1. **Dispatch on `r`'s first switch.** `measure_geodesic`'s scan, rerun with route F below the
   threshold, reads 0.9929 / **1.4402** / 1.6417 `u` at `R_F64.below` (`θ ≈ 0.58`). It dominates
   all three strata, and `geo:generic` improves on the shipped 1.7382. It would also send far more
   pairs to the arm without a transcendental. `0051` rejected that threshold because the provided
   body lost there, at 2.5019. Raising it owes `f32` evidence that does not exist yet (`0051`
   *Further work* 2) and its own bench rows.
2. **Kernel latency.** Estrin's scheme, or a shorter tier for `s ≲ 1e-8`, measured against CO.9's
   Horner in a consumer's bench.
3. **The two-arm functions' register pressure.** Whether `geodesic_long` out of line
   (`#[inline(never)]`) separates the two arms' allocation for less than a call costs. Not
   measured here. It would be a trade on every long-arm geodesic, and this record does not need it.
