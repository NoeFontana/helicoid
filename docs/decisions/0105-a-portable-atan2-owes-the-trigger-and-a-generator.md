# 0105: A portable `atan2` owes the trigger and a generator

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #205: `helicoid_linalg::math::{atan2, atan2f}` behind `Real::atan2`; `math/generated.rs` from `just atan-fit`; `helicoid:libm-atan2` the planted twin.

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

While this record was a draft, `helicoid:portable-atan2` was the shipped program at a scalar whose
`atan2` alone was the candidate, at both precisions, planted so that no bar read it; the readings
below are its. The kernels are now `helicoid_linalg::math`, and the twin is inverted:
`helicoid:libm-atan2` is the shipped program with `libm`'s `atan2`.

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

**At the swap** (`cargo xtask compare --a helicoid:libm-atan2 --b helicoid`, both at the
regenerated switches below):

- `f64`: 124 records better and 80 worse ($p = 2.5\times10^{-3}$); 27 stratum maxima fall and 7
  rise, by at most 1.449× (`sen3_log_n2`'s `rho:1e3/theta=1`, 0.820 to 1.188 `u`). `so3_log`'s
  maximum falls from 1.680 to 1.405 `u`, `coeff_r`'s from 36.64 to 28.11, `sim3_log`'s from 3.533
  to 2.796.
- `f32`: 20 better and 8 worse; 4 maxima fall and 1 rises (`sim3_log`'s `sigma:1e-8/theta:3@f32`,
  1.158×).
- **`r`'s switch moves.** `cargo xtask thresholds` reads the new exact arm and moves the switch at
  both precisions: at `f64` from $n^2/w^2 < 0.0898$ to $0.0866$, the derivative objective from
  36.64 to 28.11 `u`; at `f32` from $0.379$ to $0.328$ with one series term fewer (15 to 14), the
  derivative objective from 10.41 to 11.19 `u`.
- `envelope --check`: no domination failure. No `helicoid` baseline is blessed yet, so the
  no-regress bar does not run.

## Decision

Settled by the owner on 2026-10-24: **swap `atan2` at both precisions, and make the kernels public.**
The binary64 kernel is within 0.504 ulp of the reference on every draw read, against `libm`'s
1.40; the binary32 kernel is correctly rounded on every draw read, against `atan2f`'s 1.47 ulp. The
price is 1.1–1.25× `libm`'s latency, 1.07× on `so3_log`'s exact arm and nothing on its series
arms. `0034`'s trigger is a speedup and does not fire; the swap is taken for accuracy, which
`0034` (iii) admits.

1. **`helicoid_linalg::math`, public.** [`atan2`](../../crates/helicoid-linalg/src/math.rs) and
   `atan2f` are free functions, `libm`'s names and argument order, so a consumer migrating off
   `libm` changes a path. `Real::atan2` for `f64` and `f32` calls them; `Dual::atan2` inherits
   its value through `Real`. Both are pure Rust, basic operations only, one fixed operation
   sequence per precision.
2. **The formulation is the candidate's** (above): 9 centres at binary64 with an exact numerator and
   a compensated quotient; binary32 evaluated in binary64, folded at ½. The rescale arm is keyed on
   $\max(\lvert x\rvert, \lvert y\rvert)$ alone.
3. **Generated coefficients and tables.** `conformance/generate/fit/` (`python -m fit atan2`) fits
   $(\arctan\sqrt s/\sqrt s - 1)/s$ by weighted Remez in mpmath, rounding the coefficients one at
   a time with a refit after each. It tabulates $\arctan(j/C)$ as pairs, splits $\pi/4$, and
   writes `crates/helicoid-linalg/src/math/generated.rs`. `just atan-fit-check`, a CI job, fails
   on drift, and the generated-file lint holds the header. The package sits beside `gen`, not in
   it, because `gen`'s sources are hashed into the corpus's provenance. No coefficient appears in
   this record or by hand in any source file.
4. **`libm`'s `atan2` stays as the twin.** `helicoid:libm-atan2`, planted, is the shipped program
   with `libm`'s `atan2` at both precisions, so `cargo xtask compare` re-measures the choice
   (`0078`). `atan2` is a primitive, checked against the corpus, so it has no `NUMERICS.md` §14
   row; `math_tests` keeps a two-ulp smoke bound against `libm`.
5. **The strata `0006` asks for.** New `real_atan2` strata cover where the formulas meet: the
   centre midpoints, the binary32 fold, the diagonal, the axes, and the window's and the
   exponent range's edges. They follow the swap (Implementation plan, step 4).

## Rationale

The measurement settles the latency case: without a fused multiply-add, compensation is what an
`atan2` within 0.52 ulp costs, and it costs more than `libm`'s whole call. The accuracy case is what
changed. Before `0106`, most of `so3_log`'s exact-arm error was `r` rounded and then `r · uᵢ`
rounded, so a better `atan2` moved maxima both ways. With each component rounded once, the exact
arm's maxima follow `atan2`'s error: `so3_log` falls from 1.680 to 1.405 `u` overall, and none of
its strata rises. Six small strata of other ids rise, by up to 1.449×; the owner takes that for a
`Log` whose dominant term is now correctly rounded in all but 268 of $3 \times 2^{22}$ draws.

**A module, not a crate.** A `helicoid-math` crate was weighed. Every consumer already depends on
`helicoid-linalg`, so a crate buys no reach today and costs a publish, an MSRV, a `deny` surface and
a `0007` budget line. On `0.0.x` a later split is a move plus a re-export. The trigger for one is a
consumer that wants the kernels without the linear algebra, or the ladder growing past `atan2`.

**Not upstream, yet.** Upstreaming to the `libm` crate would keep D16 as worded. A correctly rounded
`atan2` is the form worth upstreaming (open question 3 records the scratch port); this kernel is
within 0.504 ulp, not correctly rounded, and waiting on a `libm` release would block the consumers.

## Consequences

- The bits of every `atan2` consumer move once: `so3_log`, `s2_local`, the SO(2)/SE(2) `Log`,
  `sim3_log`, `coeff_r`'s exact arm, `real_atan2`. On `0.0.x` any release may move them; from 1.0
  that is a minor-version event (`0034` item 4).
- **D16** reads: for a given `helicoid` version and resolved `libm`, outputs are bit-identical
  across targets; `atan2`'s bits are fixed by the `helicoid` version alone. `just determinism`
  checks it. The `kernel` lint still keeps every `libm::` call in `float.rs`; `math` has none.
- `coeffs/generated.rs`'s `Kernel:` header names `math`'s `atan2`, and the sweep is regenerated in
  the swap's PR; `r`'s switch moves (above).
- **`math` is public API** (`API.md` §2): its accuracy statement is a measured maximum, not a proof,
  and a tuning that changes a bit is a minor-version event from 1.0.
- A committed generator joins `corpus` and `thresholds` among the drift-gated artefacts, and the
  generated-file registry names a `just` recipe as well as `cargo xtask` tasks.

## Implementation plan

1. This record. Verified by `just lint` and `just doc`. **Done.**
2. The coefficient generator and its committed output. Verified by `just atan-fit-check`, a CI job,
   and `tests/test_fit_atan2.py`. **Done.**
3. The candidate as a planted twin, measured: special values, odd symmetry, `2^k` invariance, the
   exact numerator and denominator, monotone sweeps where formulas meet, the error readouts.
   **Done**; the tests moved with the kernels to `helicoid-linalg`'s `math_tests`.
4. The new `real_atan2` strata. Verified by `just corpus-check` and
   `just compare --a helicoid:libm-atan2 --b helicoid`. **Done**, after the swap: 8 binary64 and 6
   binary32 strata, 896 records, inside `Dual::atan2`'s derivative domain (§12), so the upper window
   edge and the exponent range's ends stay the kernel's unit tests'. Every stratum's maximum, read
   over all fields, is the derivative's and does not move; the value alone falls from 1.579 to
   0.939 `u` at `f64` and from 1.679 to 0.981 at `f32`, 69 and 53 records better, none worse.
5. Latency: `xtask/benches/atan2.rs`, in three regimes and with `so3_log` end to end. **Done** on one
   x86_64 host; aarch64 and wasm32 are not timed.
6. The swap: `math`, `Real::atan2` through it, `helicoid:libm-atan2`, the regenerated sweep,
   results and baselines, `API.md`, `NUMERICS.md`, D16 and a changelog entry. Verified by
   `just test`, `just conformance`, `just envelope --check` and `just determinism`. **Done.**

## Open questions

All answered on 2026-10-24.

1. *Is a near-correctly-rounded `atan2` worth 1.07× on `so3_log`'s exact arm?* Yes, at both
   precisions (Decision).
2. *(Answered by the tuning.)* The $k = 0$ gap closes when the tail shrinks: 9 centres $j/8$
   reach 0.503 ulp with a 5-term polynomial and no compensated tail.
3. *Upstream to the `libm` crate instead?* Not now: the kernels live in `helicoid-linalg`, public
   (Rationale). A correctly rounded `atan2` remains the form worth upstreaming. CORE-MATH has one
   for binary64 (MIT; its fast phase's bound is empirical, so correct rounding is not proven), and
   `libm` 0.2.16 already ships CORE-MATH's `cbrt`. A scratch port without a fused multiply-add read
   correctly rounded on $4 \times 2^{22}$ draws at about 1.3× `libm`'s latency natively and 1.2×
   on wasm32.
4. *Does binary32 get the binary64-evaluated kernel?* Yes. It costs nothing in throughput and
   1.2–1.35× in latency, and moves `atan2f` from 1.47 ulp to correctly rounded on every draw read.
5. *Where does the twin live?* As a planted subject, `helicoid:libm-atan2`, and as `math_tests`'
   smoke bound; `atan2` is a primitive, so no §14 row (Decision item 4).

## Further work

- `Dual::atan2`'s derivative overflows outside the larger argument's $10^{\pm154}$ (§12). A
  power-of-two rescale of $x_v^2 + y_v^2$ would widen it without moving the value. That is a
  separate record with no consumer yet.
