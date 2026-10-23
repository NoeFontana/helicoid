# 0105: A portable `atan2` owes the trigger and a generator

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** none; plan steps 2–5 are measurement, and no library code moves before `ready`.

## Context

A proposal replaces the `libm` crate's `atan2` behind `Real::atan2` with a portable kernel, generic
over `S: Real` and built from basic operations only. It has five parts:

- a fold before the division: `n − d` over `n + d` when `d/2 < n`, so the numerator is exact
  (Sterbenz);
- an exact `k · π/4` reconstruction, `k ∈ {0, …, 4}`, using a high part of `π/4` with trailing
  zero bits;
- a quotient compensated by Dekker's product and `fast_two_sum`, with no fused multiply-add;
- an odd minimax polynomial evaluated by Estrin;
- a power-of-two rescale for arguments near the ends of the exponent range.

It claims at most 0.52 ULP and 3–4× `libm`'s speed, with outputs bit-identical across x86_64,
aarch64 and wasm32.

That is ladder step (iii) of [`0034`](./0034-a-second-transcendental-backend-owes-a-measurement.md):
"faster **portable** kernels: pure Rust, basic operations only, a fixed operation sequence,
restricted to the Lie argument domains, verified by the corpus". `0034` puts the whole ladder
behind a trigger: a committed consumer benchmark where the measured transcendental share $f$ and
backend speedup $k$ clear tf_tree's `0016` bar **end to end**. No such benchmark exists. The one
profile in the stack, tf_tree's lookup, has no `libm` row.

**What a better `atan2` would buy in accuracy is already attributed.**
[`0090`](./0090-log-scales-its-vector-by-an-exact-two.md) places `so3_log`'s two resolved losses
on `r`'s series arms, where no `atan2` is called. D16's `atan2` costs only on the exact arm
($\theta \gtrsim 0.58$ at `f64`), for example 2.317 `u` against sophus-rs's 1.368 at $\theta = 1$.
The glibc twin of [`0037`](./0037-six-records-cannot-answer-a-one-variable-question.md) (draft)
moves 15 of 751 rows, and one of them (`coeff_d`, `theta:1e0`) gets **worse**. The `real_atan2`
rows are committed at 2.20, 2.05 and 1.58 `u` (`yx:generic`, `yx:ratio-1e-8`, `yx:ratio-1e-100`).
So the accuracy gain is bounded and specific: the exact arm of `r` and the `real_atan2` strata.

**Two of the proposal's parts already exist.** `Dual::atan2` takes its derivative analytically,
$(x_v y_d - y_v x_d)/(x_v^2 + y_v^2)$, with the domain `NUMERICS.md` §12 states. The fused
small-angle `log_ratio` is `coeffs::log_ratio`'s series arm in $s = n^2/w^2$, behind a switch that
`cargo xtask thresholds` generates. The proposal's own version uses a fixed `n ≤ w/2` switch and
an inline series, and both break hard rules: `0004` (switches are generated), and the rule that
only `coeffs` evaluates a cancelling coefficient.

## Measurement

**Polynomial length.** The proposal states that an odd polynomial of degree 11–13 on
$|u| \le \tfrac12$ reaches $2^{-61.4}$ relative error. Measured here:

- Protocol: the Chebyshev interpolant of $g(s) = \arctan(\sqrt s)/\sqrt s$ on $[0, S]$, in mpmath at
  60 digits; its maximum relative error on a 2001-point grid. Interpolation at Chebyshev points
  comes within a small factor of minimax, so the table puts a close lower bound on the length a
  minimax fit needs.
- Each row is $m$ terms in $s$, an odd polynomial of degree $2m - 1$ in $u$.

| terms in $s$ | degree in $u$ | $\lvert u\rvert \le \tfrac12$ (the proposal) | $\lvert u\rvert \le \tan\tfrac\pi8$ | $\lvert u\rvert \le \tfrac7{16}$ (FDLIBM) |
|---|---|---|---|---|
| 4 | 7 | $2^{-18.7}$ | — | $2^{-20.1}$ |
| 6 | 11 | $2^{-27.6}$ | $2^{-30.6}$ | $2^{-29.7}$ |
| 7 | 13 | $2^{-31.9}$ | $2^{-35.4}$ | $2^{-34.4}$ |
| 11 | 21 | $2^{-49.2}$ | $2^{-54.7}$ | $2^{-53.1}$ |
| 12 | 23 | $2^{-53.5}$ | $2^{-59.5}$ | $2^{-57.7}$ |
| 13 | 25 | $2^{-57.8}$ | $2^{-64.2}$ | $2^{-62.3}$ |

The proposal's six coefficients ($c_3 \ldots c_{13}$) give about $2^{-32}$, not $2^{-61.4}$. That
leaves about 21 bits of binary64 wrong, and the four-term binary32 kernel ($2^{-18.7}$) leaves
about 5 bits of binary32 wrong. Reaching $2^{-56}$ or better at binary64 on the proposal's interval
takes **13 terms**, against FDLIBM's 11 on its narrower one. So the two-interval fold buys no
fewer table entries and needs a longer polynomial. A symmetric fold at $\tan\frac\pi8$ would need
12 terms, but it loses Sterbenz on the numerator for $n/d \in (\sqrt2 - 1, \tfrac12)$.

**The rest of the proposal, read against the code.** Nothing below has been run.

1. **The coefficients are typed in, and they disagree with their own text.** The text says
   "degree-11, 5 terms" and the code has six. The binary32 `c5 = 0x3e4ccccd` is $0.2$ rounded, but
   its comment reads 0.19997856. Under `0004`, no coefficient is typed in any case.
2. **The window misroutes axis-aligned inputs.** It requires
   $2^{-480} \le \min(\lvert x\rvert, \lvert y\rvert)$, so `atan2(0, w)` takes the cold arm. That
   is `SO3::log` at the identity, and every `y/x` below $2^{-480}$. The cold arm then rescales by
   the maximum only, so for those inputs it changes nothing. The answer is still right; the "fast
   path for 99.999%" is not. The window should test the maximum alone.
3. **The ≤ 0.52 ULP bound is plausible, but it is not a bar.**
   - Given a polynomial correct to $\ll u$, the compensated quotient and the exact $k\,\pi/4_{hi}$
     make a sub-0.6 `u` maximum likely: the residual enters once, through the polynomial tail, and
     the variable `u3`, which also carries it, is never read.
   - The bar remains `0006`'s and `0088`'s measured maximum per stratum, not a derivation.
4. **`Mask` has no `xor`.** The octant sign is written with `and`, `or` and `not`. `S::branch` on a
   scalar is a real branch, so "branchless" holds per lane only. SIMD lanes are out of scope until
   [`0013`](./0013-simd-lanes-owe-a-measurement.md) (draft) has its measurement.
5. **The fused multiply-add route is already closed.** `cargo xtask lint` rejects `mul_add` and
   `libm::fma`, and Rust does not contract `a * b + c` on its own. The proposal's assembly check
   would guard what the lint already guards at the source.

## Decision

Proposed, not settled: **measure the kernel; ship nothing until `0034`'s trigger fires, or the
owner rules that the accuracy alone is enough (open question 1).**

If the kernel is built, it takes this shape:

1. **One private body behind `Real::atan2`** in `helicoid-linalg`'s kernel module, the D16
   chokepoint. It is pure Rust, uses basic operations only, and runs one fixed operation sequence
   per precision. `libm::atan2` stays callable as its twin (open question 5).
2. **The formulation, corrected:**
   - $n = \min$, $d = \max$ of $\lvert x\rvert, \lvert y\rvert$, with the fold at $d/2 < n$:
     Sterbenz-exact numerator, and the denominator's error recovered by `fast_two_sum(d, n)`.
   - $k\,\pi/4_{hi}$ exact for $k \le 4$, with the low part added last.
   - The quotient compensated by Dekker's product.
   - The rescale arm keyed on $\max(\lvert x\rvert, \lvert y\rvert)$ alone.
   - The interval (one fold or more) is open question 2.
3. **Generated coefficients.** A generator in `conformance/generate/` fits
   $(\arctan\sqrt s/\sqrt s - 1)/s$ by Remez in mpmath and rounds the coefficients to each
   precision. It commits the coefficients with a header naming its source and its error, and a
   `--check` recipe fails on drift. No coefficient appears in this record or by hand in any source
   file.
4. **Verified as a candidate before it replaces anything.**
   - A planted one-variable twin, `seeded:portable-atan2`, follows `seeded:host-std`'s pattern:
     `f64` with `atan2` alone swapped, validated by reproducing `seeded:correct`'s rows wherever
     `atan2` is not reached.
   - New `real_atan2` strata cover what the corpus lacks: the fold boundary $n = d/2 \pm k$ ulp,
     axis-aligned inputs, angles near $k\pi/4$, and the window's edges.

## Rationale

The ladder's order is `0034`'s, and the arithmetic still holds. At $f = 0.15$, even $k = 4$ gives
$1/(0.85 + 0.0375) = 1.127\times$ end to end, right at the bar tf_tree's `0016` set for rejecting a
dependency. $f$ for `atan2` alone is a fraction of that 0.15, because `Exp`, `act` and every
geodesic call `sin_cos`, not `atan2`. A latency case therefore needs a consumer measurement first.
The accuracy case is real but narrow, as `0090` bounds it.

Building the generator and the twin first costs nothing irreversible. Neither touches a library
crate, and both are the instruments the swap would need anyway. A kernel shipped first and
measured later would move every `atan2` consumer's bits (`so3_log`, `s2_local`, the SO(2)/SE(2)
`Log`, `sim3_log`, `coeff_r`'s exact arm) and the committed baselines, on a claim this record has
already found wrong in its main number.

The alternative of upstreaming the kernel to the `libm` crate keeps `0007` at "`libm` only", and
`0034` prefers it where possible. It is open question 3, not a rejection.

## Consequences

- Nothing changes in a library crate until this record is `ready`.
- If the kernel ships:
  - the bits of every `atan2` consumer move once. On `0.0.x` any release may move them; from 1.0
    that is a minor-version event (`0034` item 4);
  - `coeffs/generated.rs`'s `Kernel:` header names the new kernel, and `thresholds-check`,
    the conformance results and the baselines are regenerated in that PR;
  - D16's property, bit identity across targets for a given `helicoid` version, holds by
    construction for a kernel of basic operations, and `just determinism` checks it;
  - D16's resolved-`libm` qualifier no longer covers `atan2`.
- A new committed generator joins `corpus` and `thresholds` among the drift-gated artefacts.

## Implementation plan

1. This record, draft. Verified by `just lint` and `just doc`.
2. The coefficient generator and its committed output, with `--check` in `just lint`'s neighbours.
   Verified by the generator's own `--check`, and by a printed minimax error at or below the
   precision's target on the chosen interval.
3. `seeded:portable-atan2`, planted. Unit tests check: special values equal to `libm::atan2`
   (±0, ±∞, NaN, all four signed-zero quadrants); exact odd symmetry; monotone radial sweeps; the
   fold boundary; the window edges. Verified by
   `cargo nextest run -p xtask -- portable`, and by reproducing `seeded:correct` wherever `atan2`
   is unreached.
4. The new `real_atan2` strata, then the comparison: max and p99 `u` per stratum for `libm`, the
   candidate and `seeded:host-std`. Verified by `just corpus-check` and
   `just conformance --subject seeded:portable-atan2`.
5. Latency: an `atan2` group in `crates/helicoid/benches/linalg.rs` in three regimes (one octant,
   random octants, axis-aligned), `so3_log` end to end, and $f$ and $k$ from tf_tree's lookup
   bench. Verified by paired runs; the numbers go into this record.
6. On `ready` only: the swap in the kernel module, with the twin, the regenerated artefacts and a
   changelog entry. Verified by `just test`, `just conformance`, `just envelope` and
   `just determinism`.

## Open questions

1. Is the accuracy gain on `r`'s exact arm and on `real_atan2` enough to proceed without
   `0034`'s latency trigger?
2. Which interval: one fold to $\lvert u\rvert \le \tfrac12$ (13 terms), a symmetric fold at
   $\tan\frac\pi8$ (12 terms, numerator rounded and compensated), or FDLIBM's four intervals
   without a table (11 terms, more selects)? Step 5 measures the cost.
3. Should the kernel be upstreamed to the `libm` crate instead of living in `helicoid-linalg`?
4. Does binary32 get its own kernel, or stay on `libm::atan2f`? Under `0016` the `@f32` strata
   are the coefficient ids' alone, so a binary32 kernel would also owe strata.
5. Where does the twin live? `NUMERICS.md` §14's twins are composite routines. `libm::atan2` as
   a test-only twin needs the `kernel` containment lint to admit a test module.

## Further work

- `Dual::atan2`'s derivative overflows outside the larger argument's $10^{\pm154}$ (§12). A
  power-of-two rescale of $x_v^2 + y_v^2$ would widen it without moving the value. That is a
  separate record with no consumer yet.
