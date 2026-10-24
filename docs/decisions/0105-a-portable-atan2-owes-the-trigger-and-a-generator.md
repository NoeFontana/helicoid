# 0105: A portable `atan2` owes the trigger and a generator

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** none in a library crate; the generator, the twin and `xtask/benches/atan2.rs` are measurement.

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
on `r`'s series arms, where no `atan2` is called. D16's `atan2` is called only on the exact arm
($\theta \gtrsim 0.58$ at `f64`), which then read 2.317 `u` against sophus-rs's 1.368 at
$\theta = 1$; [`0106`](./0106-log-rounds-its-vector-once-on-the-exact-arm.md) has since rounded
that arm's vector once, and `so3_log` reads 1.680 `u`, below every oracle on every stratum.
The glibc twin of [`0037`](./0037-six-records-cannot-answer-a-one-variable-question.md) (draft)
moves 15 of 751 rows, and one of them (`coeff_d`, `theta:1e0`) gets **worse**. The `real_atan2`
rows are committed at 2.20, 2.05 and 1.58 `u` (`yx:generic`, `yx:ratio-1e-8`, `yx:ratio-1e-100`).
So the accuracy gain is bounded and specific: the exact arm of `r`, and `atan2`'s own value.

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
3. **The ≤ 0.52 ULP bound does not hold for this formulation.** The residual enters once,
   through the polynomial tail, and the variable `u3`, which also carries it, is never read. With
   a polynomial correct to $2^{-56.8}$ the maximum is still 0.78 ulp, all of it on one arm (below).
4. **`Mask` has no `xor`.** The octant sign is written with `and`, `or` and `not`. `S::branch` on a
   scalar is a real branch, so "branchless" holds per lane only. SIMD lanes are out of scope until
   [`0013`](./0013-simd-lanes-owe-a-measurement.md) (draft) has its measurement.
5. **The fused multiply-add route is already closed.** `cargo xtask lint` rejects `mul_add` and
   `libm::fma`, and Rust does not contract `a * b + c` on its own. The proposal's assembly check
   would guard what the lint already guards at the source.

### Tuning

A scratch harness, not committed, read about twenty variants against the double-double reference
($2^{21}$ draws, $\lvert y/x\rvert$ within $2^{\pm4}$). It also timed each one, with each call's result
fed into the next call's `y`, on an AMD EPYC Milan (Zen 3) at about 3.7 GHz. What it found:

1. **Rounding the coefficients one at a time keeps what naive rounding loses.** The generator now
   rounds each coefficient and then refits the rest. On $\lvert u\rvert \le \tfrac12$ with 13
   terms that keeps $2^{-62.0}$ where naive rounding kept $2^{-56.8}$, at no run-time cost.
2. **Only a smaller tail, or a compensated cubic term, closes the $k = 0$ gap.**

   | variant | max ulp |
   |---|---|
   | fold at $\tfrac12$, 13 terms (this PR's first candidate) | 0.780 |
   | the same, coefficients refit one at a time | 0.724 |
   | … and $s = u^2$ and $u s$ each compensated | 0.720 |
   | … and the cubic term $c_0 u^3$ compensated | 0.528 |
   | symmetric fold at $\tan\frac\pi8$, 12 terms, numerator by `fast_two_sum` | 0.628 |
   | … and the cubic term compensated | 0.515 |
   | 5 centres $j/4$, 7 terms (the reduction below) | 0.518 |
   | 9 centres $j/8$, 5 terms (the candidate) | 0.502 |

3. **The proposal's 3–4× speedup does not exist on this host.** Division plus polynomial with no
   compensation at all (12 terms, 1.49 ulp) already costs `libm`'s latency, about 23 ns. Without
   a fused multiply-add, every compensated variant above costs 30–37 ns.
4. **The residual belongs behind a branch.** In one basic block, LLVM's SLP vectorizer at the
   default `-O3` packs the residual's products into vector lanes, and the shuffles land on the
   critical path: 35 ns, against 25 ns under `-C llvm-args=-vectorize-slp=false`. Computing the
   residual under the branch that already skips it below the window (`if |num| ≥ 2⁻⁹⁶⁹`) puts it
   in a block of its own: the same bits on $2^{22}$ draws, at the SLP-off time. Skipping it
   where it is a small part of an ulp does not pay: at most 0.52 ulp needs it wherever
   $\lvert v\rvert$ exceeds $2^{-6}$ of the rest of the answer, and the test, late on the critical
   path, made the kernel slower (31–33 ns).
5. **Binary32 evaluated in binary64 needs no compensation.** Folded at $\tfrac12$ with 8 terms
   ($2^{-40.6}$), it misrounds **none** of $2^{21}$ draws. With 7 terms it misrounds 28, at most
   0.50016 ulp. `libm::atan2f` misrounds 10% of them, at most 1.47 ulp.
6. **Neither one division nor an integer centre index pays** (the committed bench, below). One
   division, $r = 1/(d + c\,n)$ and $v = (n - c\,d)\,r$ with the residual taken against that
   $v$, reads the same or lower maxima on $3 \times 2^{22}$ draws but costs 28.0 ns against
   27.3: the multiply lengthens the chain, and the divider is not the throughput bound (14.6
   against 14.8 ns). Counting the midpoints passed as an integer and converting it to $c$ puts a
   popcount and an integer-to-float conversion ahead of the division: 39.7 ns against 27.3. $c$
   is a tree sum of floating-point terms, and the table index is derived from it off the critical
   path.

### The candidate, measured

`helicoid:portable-atan2` is the shipped program at a scalar whose `atan2` alone is the candidate,
at both precisions. It is planted, so no bar reads it. Its kernels are in
`xtask/src/shipped/portable_atan2/kernel.rs`, and its fit is `conformance/fits/atan2.json`, read
at compile time, so a fit of another shape is a build error.

- **binary64.** $t = n/d$ is reduced about the nearest of 9 centres $c = j/8$:
  $\arctan t = \arctan c + \arctan v$, with $v = (n - c\,d)/(d + c\,n)$ and $\lvert v\rvert \le \tfrac1{16}$.
  - Because $c$ has at most 4 significant bits, $c\,d$ and $c\,n$ are each two exact products once
    the low 4 bits of $d$ and $n$ are split off. $n - c\,d$ is then exact (Sterbenz, then a
    representable difference), and $d + c\,n$ is an exact pair: both of its low parts are
    multiples of $\mathrm{ulp}(n)/8$ within $\mathrm{ulp}(d + c\,n)$, so their sum rounds nothing.
  - The quotient's residual takes Dekker's halves.
  - The tail $v s P(s)$ is at most $2^{-9.6}$ of the answer, so it is not compensated.
  - $P$ has 5 terms at $2^{-62.56}$, and $\arctan(j/8)$ is a table of `hi + lo` pairs.
- **binary32.** Evaluated in binary64 and folded at $\tfrac12$. The numerator and denominator are
  exact for binary32 operands. $P$ has 8 terms at $2^{-40.58}$, and the result is rounded to
  binary32 once.

To share the dispatch, the shipped subject's body became `shipped::program::<D, F>`. Its 123
`determinism` digests are unchanged.

**Accuracy against the double-double reference.** Each row is $2^{20}$ draws. The `f32` row reads
against `libm`'s binary64 `atan2`.

| | candidate | `libm` | host glibc 2.39 |
|---|---|---|---|
| `f64`, $\lvert y/x\rvert$ within $2^{\pm4}$ | 0.5013 ulp | 1.402 | 0.519 |
| `f64`, within $2^{\pm60}$ | 0.5008 | 1.377 | 0.515 |
| `f64`, worst arm ($j = 0$, $q = 0$, another seed) | 0.5031 | — | — |
| `f32`, within $2^{\pm4}$ and $2^{\pm120}$ | 0.5000 | 1.471 (`atan2f`) | — |

On $2^{22}$ draws per range the binary64 maximum is 0.5032 ulp within $2^{\pm4}$, 0.5025 within
$2^{\pm60}$ and 0.5010 within $2^{\pm1000}$; 268 of the $3 \times 2^{22}$ draws are above half an
ulp.

- **Exactness and monotonicity.** Unit tests check, in integer arithmetic, that $n - c\,d$ and
  $d + c\,n$ are exact, and that the kernels are monotone through every midpoint between centres,
  the fold and the swap.
- **Special values.** At both precisions the candidate returns IEEE 754's limits, correctly
  rounded. `libm::atan2f(±∞, 1)` returns $\pi/2$ one ulp low (`0x3fc90fda`); the candidate returns
  the nearest (`0x3fc90fdb`).

**Time per call** (`cargo bench -p xtask --bench atan2`, the same host, pinned to one core,
medians of 3 interleaved rounds, one binary per variant). Each row is 256 calls, either chained
through `y` (latency) or independent (throughput). `mixed` draws random octants; `axis` has
$y = 0$, where `libm` returns early. `so3_log` is `Log` with `0106`'s exact arm, $\theta$ uniform
in the range named, about random axes.

| | `libm` | candidate | host glibc |
|---|---|---|---|
| `f64` latency, mixed / one octant | 24.5 / 22.1 ns | 27.3 / 27.4 | 23.1 / 24.2 |
| `f64` throughput, mixed | 7.2 ns | 14.8 | 11.5 |
| `f32` latency, mixed / one octant | 18.0 / 15.8 ns | 21.9 / 21.4 | 23.6 / 21.0 |
| `f32` throughput, mixed | 5.9 ns | 6.0 | 7.8 |
| `so3_log`, $\theta \in [0.6, 3.1)$ (the exact arm) | 42.7 ns | 45.6 | — |
| `so3_log`, $\theta \in [0.01, 0.5)$ (a series arm, no `atan2`) | 11.9 ns | 11.9 | — |

The binary64 kernel costs 1.1–1.25× `libm`'s latency and 2× its throughput. On `so3_log`'s exact
arm, whose rows are independent calls, that is 2.9 ns of 42.7, 1.07×. An earlier reading here,
13.9 against 28.4 ns, came before `0106` and from a bench whose $\theta$ ran past each row's range
(37 of the series row's 256 rotations were on the exact arm) and whose candidate looked its fit
up on every call. The binary32 kernel matches `atan2f`'s throughput at 1.2–1.35× its latency.

**On the corpus** (`cargo xtask compare --a helicoid --b helicoid:portable-atan2`, with `0106`'s
`Log`):

- **`f64` overall:** 201 of 65 051 records move, 123 better and 78 worse (sign test
  $p = 1.8\times10^{-3}$). Of 3145 strata, 27 maxima fall, by down to 0.481×, and 6 rise, all in
  strata of 4–12 records: `sen3_log_n2`'s `rho:1e3/theta=1` by 1.449×, `se3_world_local`'s
  `rho:1e0/theta=1` by 1.286×, `sim3_log`'s `sigma:700/theta:1` by 1.268× (at 0.022 `u`).
- **`f32` overall:** 26 records move, 18 better and 8 worse ($p = 0.076$). Of 1420 strata, 3
  maxima fall and 2 rise, the largest `sim3_log`'s `sigma:1e-8/theta:3@f32` by 1.158×.
- **`real_atan2`:** at `f64`, 13 of 448 values move and all are better ($p = 2.4\times10^{-4}$);
  at `f32`, 6 move and all are better. But **no stratum maximum moves**: the committed maxima are
  read from fields other than the value.
- **`so3_log`:** 34 records better and 26 worse; 3 maxima fall and none rises. The overall
  maximum falls from 1.680 to 1.405 `u`, `theta:1e0` from 1.680 to 1.108, `theta:dense` from
  1.578 to 1.271 and `theta:1e-1` from 1.298 to 1.137. Before `0106` the same candidate moved
  `so3_log`'s maxima both ways (`theta:1e0` 1.741 to 2.175); with each component of `φ` rounded
  once, the exact arm's maxima now follow `atan2`'s error.

## Decision

Proposed, not settled: **do not swap binary64's `atan2`.** The candidate is within 0.503 ulp of
the reference, but it is slower than `libm` (1.1–1.25× latency, 2× throughput, 1.07× on
`so3_log`'s exact arm). So `0034`'s trigger, a speedup, cannot fire. Since `0106`, maxima move:
`so3_log`'s fall (1.680 to 1.405 `u` overall) and none of its strata rises, but 6 strata of other
ids rise, by up to 1.449×. The binary32 kernel is the one open trade: correctly rounded on every
draw read, at `atan2f`'s throughput (open question 4).

If the kernel is built, it takes this shape:

1. **One private body behind `Real::atan2`** in `helicoid-linalg`'s kernel module, the D16
   chokepoint. It is pure Rust, uses basic operations only, and runs one fixed operation sequence
   per precision. `libm::atan2` stays callable as its twin (open question 5).
2. **The formulation is the candidate's** (above): 9 centres at binary64 with an exact numerator,
   and binary32 evaluated in binary64. The rescale arm is keyed on
   $\max(\lvert x\rvert, \lvert y\rvert)$ alone.
3. **Generated coefficients and tables.** `conformance/generate/fit/` (`python -m fit atan2`) fits
   $(\arctan\sqrt s/\sqrt s - 1)/s$ by weighted Remez in mpmath, rounding the coefficients one at
   a time with a refit after each. It tabulates $\arctan(j/C)$ as pairs and splits $\pi/4$. It
   writes `conformance/fits/atan2.json` with the source and both errors, and
   `just atan-fit-check` fails on drift. The package sits beside `gen`, not in it,
   because `gen`'s sources are hashed into the corpus's provenance. No coefficient appears in this
   record or by hand in any source file.
4. **Verified as a candidate before it replaces anything.**
   - A planted one-variable twin, `helicoid:portable-atan2`, follows `helicoid:host-std`'s
     pattern: the shipped program at a scalar with `atan2` alone swapped, at both precisions.
   - New `real_atan2` strata cover what the corpus lacks: the fold boundary $n = d/2 \pm k$ ulp,
     axis-aligned inputs, angles near $k\pi/4$, and the window's edges.

## Rationale

The measurement settles the latency case: without a fused multiply-add, compensation is what an
`atan2` within 0.52 ulp costs, and it costs more than `libm`'s whole call. The ladder's order is
`0034`'s, and the arithmetic still holds. At $f = 0.15$, even $k = 4$ gives
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
2. The coefficient generator and its committed output. Verified by `just atan-fit-check`, a CI job,
   and `tests/test_fit_atan2.py`. **Done.**
3. `helicoid:portable-atan2`, planted. Unit tests check:
   - special values against IEEE 754's correctly rounded limits;
   - exact odd symmetry;
   - bit-exact invariance under $2^k$ scaling, through the window and the cold arm;
   - the exact numerator against the reference;
   - monotone sweeps through every place two formulas meet;
   - at most 0.51 ulp at binary64 and 0.5001 at binary32 against the reference.

   The readouts above are its ignored measurements. Verified by
   `cargo nextest run -p xtask -- portable_atan2`, and by the shipped subject's unchanged
   `determinism` digests. **Done.**
4. The new `real_atan2` strata, then the comparison: max and p99 `u` per stratum for `libm`, the
   candidate and `seeded:host-std`. Verified by `just corpus-check` and
   `just compare --a helicoid --b helicoid:portable-atan2`. **Owed** if a kernel goes forward.
5. Latency: `xtask/benches/atan2.rs`, in three regimes and with `so3_log` end to end, read above.
   **Done** on one x86_64 host; aarch64 and wasm32 are not timed. $f$ and $k$ from tf_tree's
   lookup bench are moot while $k < 1$.
6. On `ready` only: the swap in the kernel module, with the twin, the regenerated artefacts and a
   changelog entry. Verified by `just test`, `just conformance`, `just envelope` and
   `just determinism`.

## Open questions

1. Is a near-correctly-rounded `atan2` worth 1.07× on `so3_log`'s exact arm? Before `0106` it
   was not: most of that arm's error was `r` rounded, then `r · uᵢ` rounded. With `0106` each
   component rounds once, and the candidate lowers `so3_log`'s maximum from 1.680 to 1.405 `u`
   with no `so3_log` stratum rising, while 6 strata of other ids rise. This record still proposes
   no binary64 swap; the trade is the owner's to reweigh.
2. *(Answered by the tuning.)* The $k = 0$ gap closes when the tail shrinks: 9 centres $j/8$
   reach 0.503 ulp with a 5-term polynomial and no compensated tail.
3. Should the kernel be upstreamed to the `libm` crate instead of living in `helicoid-linalg`?
   A correctly rounded `atan2` is the form worth upstreaming. CORE-MATH has one for binary64
   (MIT; its fast phase's bound is empirical, so correct rounding is not proven), and `libm`
   0.2.16 already ships CORE-MATH's `cbrt`. Correct rounding makes the bits the same on every
   target with or without a fused multiply-add, but no fused multiply-add is reachable here
   anyway: `f64::mul_add` is not in `core`, the intrinsics need `unsafe` (`0007`), and
   `libm::fma` dispatches at run time on x86_64. A scratch port without one read correctly
   rounded on $4 \times 2^{22}$ draws at about 1.3× `libm`'s latency natively and 1.2× on
   wasm32. Upstream keeps D16 as worded.
4. Does binary32 get the binary64-evaluated kernel? It costs nothing in throughput and 1.2–1.35×
   in latency, and moves `atan2f` from 1.47 ulp to correctly rounded on every draw read, but no
   `@f32` maximum falls for it. Under `0016` the `@f32` strata are the coefficient ids' alone, so a
   binary32 kernel would also owe strata.
5. Where does the twin live? `NUMERICS.md` §14's twins are composite routines. `libm::atan2` as
   a test-only twin needs the `kernel` containment lint to admit a test module.

## Further work

- `Dual::atan2`'s derivative overflows outside the larger argument's $10^{\pm154}$ (§12). A
  power-of-two rescale of $x_v^2 + y_v^2$ would widen it without moving the value. That is a
  separate record with no consumer yet.
