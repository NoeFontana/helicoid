# 0039: The sweep's grid stops below its own optimum

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** records a measurement and proposes a `PHASE1.md` §6 amendment; no code change here.

## Context

`PHASE1.md` §6 is **NORMATIVE** and fixes the sweep's search space: *"Series terms `m ∈ {1, …, 8}`;
switch point in the branch variable on a log grid of 64 points per decade, spanning
`θ ∈ [10⁻⁸, 1]`."* [`0004`](./0004-switch-points-are-generated-not-typed.md) then makes whatever
that search returns the shipped switch, and forbids typing one.

**Eight of the fifteen generated switches sit exactly on the grid's upper bound**, `θ² < 1e0`:
`K_F64` and `D_F64`, and six of the seven θ-based binary32 switches (`K`, `A`, `B`, `C`, `D`,
`E_F32`). `E_F64` sits one grid step below it, at `θ² < 0.9647`. A search whose answer is on the
boundary for more than half its subjects is not reporting an optimum; it is reporting the edge of
the box.

Nothing mathematical puts the box there. The coefficients are entire (`coefficients.md` CO.2), so
the series in `z = θ²` has no radius to respect, and `NUMERICS.md` §12's domain runs to `θ = π`,
that is `z ≤ π² ≈ 9.87`. The grid stops at `z = 1` because §6 says so.

## Measurement

One variable: the grid's span, extended by one decade to `z ≤ 10` — above the domain's own top —
with the same 64 points per decade, the same integer construction (`10^(n/64)` as the 64th root of
`2^(64P)·10ⁿ`, so a point is never a `pow` and never a binary64 rounded again,
[`0016`](./0016-f32-exact-strata-for-the-scalar-coefficient-ids.md) item 3), and the same
`m ≤ 8`. Both sweep targets rerun, the seeded kernels and the shipped `helicoid::coeffs` through
its `__sweep` feature. Reproducible from one constant plus a `pow10_pos` mirroring `pow10_neg`.

### The switches and the objective

Binary64, the sweep's own objective (max over `theta:dense` and every `theta:*` stratum), value and
derivative:

| coefficient | switch, `θ² <` | value, u | derivative, u |
|---|---|---|---|
| `k` | **1e0** → 1.9110 | 1.209 → 1.044 | 13.618 → **6.604** (2.06×) |
| `a` | 0.7234 → 0.7234 | 2.823 | 25.795 |
| `b` | 0.9647 → 0.9647 | 2.610 | 63.538 |
| `c` | 0.5623 → 0.5623 | 21.18 | 7395.5 |
| `d` | **1e0** → 1.5399 | 18.51 → 11.51 | 830.73 → **349.71** (2.38×) |
| `e` | 0.9647 → 1.6548 | 155.9 → **35.41** (4.40×) | 8687.8 → **1214.1** (7.15×) |
| `cos θ/2` | 5.623e-15 → same | 0.844 | 1.954 |
| `r` | unchanged | 2.283 | 260.26 |

The five that were not at the bound are **byte-identical**. The three that were move up, and
`coeff_e`'s derivative — **the largest error in the whole corpus** — falls by 7.15×.

### End to end, and in the direction `0038` asks for

The corpus rows of the shipped candidate confirm the objective, because the objective *is* the
corpus maximum (the sweep reads the same strata):

| row | committed grid | lifted grid | |
|---|---|---|---|
| `coeff_e` / `theta:dense` | 8687.8 u | **34.18 u** | 254× |
| `coeff_e` / `theta:1e-1` | 4273.7 u | 29.09 u | 147× |
| `coeff_d` / `theta:dense` | 112.19 u | 21.90 u | 5.1× |
| `sen3_jl_n1` / `rho:1e3/theta=1` | 4.8453 u | 1.0848 u | 4.5× |
| `sen3_jr_n1` / `rho:1e4/theta=1` | 3.5643 u | 1.4732 u | 2.4× |

Pooled over all 751 rows, which is how [`0038`](./0038-a-program-comparison-is-not-a-bar.md) (draft)
says a program comparison is read: the ceiling moves **43** rows — **37 better, 6 worse**, exact
two-sided sign test **p = 1.6 × 10⁻⁶**, geometric mean ratio **0.454**. The six that get worse are
all `sen3_j{l,r}_n*` at `theta:1e-1`, by 1.01× to 1.06×.

**Domination: 75 failures become 72. Three are fixed and none is broken**, and all three are in the
class `0037` (draft) leaves unexplained — `sen3_jl_n1` and `sen3_jr_n1` at `rho:1e3/theta=1`, and
`sen3_jr_n1` at `rho:1e4/theta=1`, each from above its best oracle to about a quarter of it.

**No group id's maximum moves.** Every `sen3_*` and `so3_*` maximum is identical to the bit, because
those maxima sit at `θ` near `π`, where the exact arm runs under either grid. The improvement is
real and it is in `θ ∈ [0.98, 1.29]`, which no group stratum maximises over. This is not a 7×
improvement in SE_N(3); it is a 7× improvement in a coefficient's derivative lane and a 2–5×
improvement on nineteen group rows that are not their stratum's worst.

### Why

`coefficients.md` CO.14 measured it before this record existed: through the **exact** arm the
relative error of the derivative is `Ĉ'u·θ^-(p+2)` with `p' = 2, 2, 4, 4, 4, 6` for `k, a, b, c, d,
e` and `Ĉ'` up to 19900 for `e`. The exact arm is therefore worst immediately **above** the switch,
the series arm worst immediately below it (CO.8 truncation), and the optimum is their crossing
(CO.10). For `e` that crossing is at `θ ≈ 1.29`. Capping the grid at `θ = 1` forces `e`'s exact arm
to run from `θ = 0.98`, where it is seven times worse than it needs to be, and the sweep had no
candidate that said otherwise.

The direction of the trade also explains why this costs nothing: §6 itself says *"the series arm is
the cheaper one"*, and moving a switch up moves inputs from the exact arm to the series arm.

### The term cap, and the two limits together

§6's other limit is `m ∈ {1, …, 8}` while the committed `coeff_series.jsonl` holds sixteen terms.
Varying both as a 2×2 — read-only, through `sweep_at`, which writes nothing — shows they are
**jointly** binding, and that measuring either alone understates the result by orders of magnitude.
Binary64, the sweep's derivative objective, in `u`:

| coefficient | `z ≤ 1`, `m ≤ 8` (committed) | `z ≤ 10`, `m ≤ 8` | `z ≤ 1`, `m ≤ 16` | `z ≤ 10`, `m ≤ 16` |
|---|---|---|---|---|
| `k` | 13.618 | 6.604 | 13.618 | **2.201** |
| `a` | 25.795 | 25.795 | 14.623 | **3.408** |
| `b` | 63.538 | 63.538 | 63.538 | **1.998** |
| `c` | 7395.5 | 7395.5 | 1846.5 | **94.478** |
| `d` | 830.73 | 349.71 | 830.73 | **1.860** |
| `e` | 8687.8 | 1214.1 | 8687.8 | **2.650** |
| `cos θ/2` | 1.954 | 1.954 | 1.954 | 1.747 |
| `r` | 260.26 | 260.26 | **36.637** | 36.637 |

Terms alone move `a`, `c` and `r` and nothing else; the ceiling alone moves `k`, `d` and `e` and
nothing else; together they move everything, `e` by **3278×**. More terms are useless if the switch
cannot rise to spend them, and a higher switch is useless if the series cannot reach it.

### Where the grid should stop: just below `π²`, not above it

`z ≤ 10` is above the domain (`π² = 9.8696`), so the sweep may put a switch where **no valid input
reaches the exact arm** and the series covers everything. Scored over the whole corpus — the
sweep's own choice, built as a subject and run through `conformance::evaluate`, writing no generated
file — that is not the best reading:

| | worst row | rows > 100 u | rows > 10 u | domination | fixed | broken |
|---|---|---|---|---|---|---|
| **A** committed, `z ≤ 1`, `m ≤ 8` | 8687.8 u | 10 | 26 | 75 | — | — |
| **B** `z ≤ 10`, `m ≤ 8` | 7395.5 u | — | — | 72 | 3 | **0** |
| **D** `z ≤ 10`, `m ≤ 16` | **94.5 u** | 0 | 6 | 75 | 15 | **15** |
| **E** `z < π²`, `m ≤ 16` | **94.5 u** | 0 | 12 | **70** | **5** | **0** |

Config **D** puts the switch above the domain and pays for it: 89 of 751 rows get worse, 92% of
them at `θ` near `π`, because the series arm replaces an exact arm that is **sub-ULP** there
(`coeff_d` at `theta:pi-1e-5` goes 0.1285 u → 1.0755 u, the worst ratio at 8.37× and the worst
absolute still under 2.7 u). Split by regime, D is a uniform **4.4×** win away from `π` (66 better,
7 worse, geometric mean 0.226) and a wash near it (112 better, 82 worse, geometric mean 0.763).

Config **E** keeps the last grid point **below** `π²` — `10^(63/64) = 9.6605` — so the exact arm
still serves `θ > 3.108`. It is strictly better than every other reading: **worst row 8687.8 u →
94.5 u (92×), no row above 100 u where the committed kernel has ten, domination 75 → 70 with five
fixed and none broken**, and pooled over all rows 67 better against 14 worse, `p = 1.9 × 10⁻⁹`,
geometric mean **0.274**. Per id:

| id | committed | E | |
|---|---|---|---|
| `coeff_e` | 8687.8 u | 19.87 u | **437×** |
| `coeff_d` | 830.73 u | 7.90 u | **105×** |
| `coeff_c` | 7395.5 u | 94.48 u | **78×** |
| `coeff_b` | 63.54 u | 5.67 u | 11.2× |
| `coeff_r` | 260.26 u | 36.64 u | 7.1× |
| `coeff_a` | 25.80 u | 3.41 u | 7.6× |
| `coeff_k` | 13.62 u | 2.69 u | 5.1× |
| every `sen3_*`, `so3_*` | unchanged | unchanged | 1.0× |

**What E does not do** is move a group id, for the reason `0037` gave: their maxima sit at `θ` near
`π`, where E leaves the exact arm in place. D moves them — `sen3_jr_n1` 8.6698 u → 3.3363 u, every
`sen3_j*` by 1.9–2.6× — and that is the whole of the case for D, against 15 broken verdicts near π.
Choosing between them is a decision and this record does not make it.

**A third limit now binds.** Under both D and E, `c` chooses **16 terms**, the cap, and its 94.5 u is
the corpus's worst row. Sixteen is also what `coeff_series.jsonl` holds (`PHASE1.md` §4.3), so the
next constraint is the corpus's own series length, not the sweep's.

### What it costs, which is the question the accuracy numbers above do not answer

Open question 6 asked whether retiring a coefficient's exact arm changes its cost. Measured with
criterion on this host, the same shape as `benches/coeffs.rs` (one call per iteration, `black_box`
on the input), on stand-in arms rather than the shipped groups:

| | ns | | | ns |
|---|---|---|---|---|
| series, 4 terms | 1.22 | | `sqrt` alone | 2.33 |
| series, 5 terms | 1.41 | | `sincos` alone | **7.02** |
| series, 8 terms | 2.41 | | one divide | 1.11 |
| series, 12 terms | 4.20 | | exact `a` (`sqrt`, `sincos`, 2 mul) | 8.07 |
| series, 16 terms | **6.43** | | exact `e` (`sqrt`, `sincos`, 8 flops, divide) | 11.46 |

**A sixteen-term Horner is cheaper than one `sincos`**, so on this host no exact arm can beat a
full-length series arm: 6.43 ns against a 7.02 ns floor. The accuracy the search space was hiding is
not paid for in latency — where it moves an input from the exact arm to the series.

Where it *doesn't* move the input, it is paid for. One switch against two, per regime:

| at `z` = | committed, 8 terms below 1 | config E, 16 terms below 9.66 | two switches, 5 then 16 |
|---|---|---|---|
| `5.6e-15` (near identity) | 2.39 ns | **6.70 ns** (2.81× slower) | **1.71 ns** (1.39× *faster*) |
| `0.81` | 2.38 ns | 6.65 ns (2.79× slower) | 6.75 ns |
| `9.0` (above every switch) | **11.64 ns** | **6.64 ns** (1.75× faster) | 6.76 ns |

So config E is **1.75× faster at large θ and 2.8× slower at small θ** — and small θ is exactly where
`PHASE3.md` §11 puts the priority ("`exp` at near-identity `θ` from `7.5e-8` to `1` first").

**That regression is avoidable for free.** A series arm's term count is set by the largest `z` it
serves and paid by the smallest; a second switch decouples them. Five terms below `z = 10⁻²` is
**1.39× faster than the committed kernel** at near-identity and truncates at about `0.004 u` there,
while the sixteen-term arm still serves the hard band. The extra branch costs 0.1 ns when taken.
Two switches are faster than today at both ends and never worse than config E.

**These are arms, not the shipped groups.** `helicoid::coeffs` evaluates a group inside one
`S::branch` and takes the exact closure as soon as one member passes its switch, so a group's cost
is which arm ran; the direction transfers and the magnitudes are the arms'. `benches/coeffs.rs`
measures the real thing and the real comparison needs the regenerated file.

### The second switch is free to search, because it is not part of the accuracy problem

Open question 6 asked what a two-switch search costs the sweep, and guessed `m₀ × m₁ × grid²`.
That is the wrong shape. Run as a joint minimisation over the sweep's own per-sample arm errors, the
two-switch search returns a **degenerate** answer on seven of eight coefficients — `m₀ = 1`,
`z₀ = 0`, the low arm serving nothing — for a structural reason:

**A second switch cannot improve the objective.** `Score::objective` is a maximum over the records,
and the long arm already attains it wherever it is selected; adding a shorter arm underneath can
only raise a record's error, never lower the maximum. So the second switch is not an accuracy
degree of freedom at all. It is a *cost* degree of freedom, and the search is two stages, not one:

1. minimise the objective with one switch — exactly today's search, unchanged;
2. then find the cheapest arm that **holds** that objective: the smallest `m₀` whose prefix maximum
   stays at or below it, and the largest `z₀` it reaches.

Stage 2 is `O(m · n)` over a prefix maximum per term count. There is no `grid²`.

### What stage 2 finds, at identical accuracy

For each coefficient, the term count that holds the one-switch objective over a given share of the
corpus's records, and the record-weighted cost of serving 90% from the short arm and the rest from
the long one (arm costs as measured above, interpolated for the lengths not benched):

| coefficient | objective | terms for 100% | 50% | 90% | 99% | two switches | one switch | |
|---|---|---|---|---|---|---|---|---|
| `k` | 2.201 u | 11 | 4 | **7** | 10 | 2.24 ns | 3.70 ns | 1.65× |
| `a` | 2.564 u | 14 | 4 | **7** | 12 | 2.40 ns | 5.34 ns | 2.22× |
| `b` | 1.998 u | 14 | 4 | **7** | 12 | 2.40 ns | 5.34 ns | 2.22× |
| `c` | 94.48 u | 16 | 3 | **7** | 14 | 2.51 ns | 6.43 ns | 2.56× |
| `d` | 1.873 u | 13 | 4 | **7** | 12 | 2.35 ns | 4.77 ns | 2.03× |
| `e` | 2.65 u | 14 | 4 | **7** | 12 | 2.40 ns | 5.34 ns | 2.22× |
| `cos θ/2` | 1.779 u | 10 | 4 | **7** | 10 | 2.19 ns | 3.20 ns | 1.46× |
| `r` | 36.64 u | 15 | 4 | 8 | 14 | 2.76 ns | 5.88 ns | 2.13× |
| **catalogue** | | | | | | **19.25 ns** | **40.00 ns** | **2.08×** |

**Seven terms hold the full objective for 90% of the records of every coefficient but `r`**, and the
last 10% is what drives the count to 10–16. One switch makes that 90% pay the last 10%'s price. The
ladder for `e`, at its 2.65 u:

| terms | holds 2.65 u for | share of records | ns |
|---|---|---|---|
| 2 | `z < 7.1e-15` | 22.2% | 0.80 |
| 3 | `z < 4.7e-7` | 46.7% | 1.01 |
| 4 | `z < 3.1e-4` | 68.5% | 1.22 |
| 5 | `z < 8.1e-3` | 79.3% | 1.41 |
| 8 | `z < 0.676` | 94.1% | 2.41 |
| 12 | `z < 7.862` | 99.0% | 4.20 |
| 14 | everywhere | 100.0% | 5.34 |

So `e`'s fourteenth term exists for **the last 1% of the domain**. A two-switch kernel at 8 and 14
costs 2.41 ns on 94% of records and 5.34 ns on 6% — **at parity with the committed kernel's
near-identity latency, and 3278× more accurate**.

**The weighting is the corpus's, not a consumer's.** Records are log-uniform across decades of `θ`
by design (`error-analysis.md` EA.21), so "90% of records" is 90% of a deliberately flat sample and
not a usage distribution; a consumer holding mostly near-identity poses gains more, one holding
mostly large rotations less. The term counts are exact; the shares and the nanoseconds are
corpus-weighted.

### Implementing it found the blocker, which is the corpus and not the sweep

Plan step 1 was built and **reverted**. The grid span, the `π²` refusal, the lowest-point failure
and the feasibility header all work, both targets regenerate, and the accuracy is what the
measurements above predict — at binary64 `e` 7.2×, `d` 2.4×, `k` 2.1×, the other five byte-identical.
At **binary32**, which the measurements above did not cover and where **six of seven** switches sat
on the old ceiling rather than two, it is far larger: `e` **283×** (4429 u → 15.65 u), `d` **80×**,
`c` **33×**, `b` 7.6×, `k` 6.1×, `a` 4.6×, and the catalogue's worst binary32 row goes
**4429 u → 49.01 u**. Per stratum, binary64 7 rows move, **7 better and 0 worse**; binary32 24 move,
**22 better and 2 worse** (`d`/`theta:pi-1e-1@f32` 5.307 → 10.40 u). Domination **75 → 72**, three
fixed and none broken.

**It does not ship, because `branch_continuity` fails and this corpus cannot say whether it should.**
`coefficients.md` CO.12 bounds the jump between the arms at a switch by the sum of their errors
there. At the new switches the measured jumps are `k` 1.084 u against 1.077 (0.70% over), `k` at
binary32 1.553 against 1.498 (3.69%), and **`d` 22.10 u against 6.30 u — 3.5×**. The bound's
right-hand side is sampled at the two corpus records that bracket the switch, which that test's own
doc already calls "a sample and not a bound", and the samples are now far from the switch:

**`k`, `d` and `e` move to `θ` = 1.38, 1.24 and 1.29, and `theta:dense` — the stratum whose stated
purpose in `PHASE1.md` §4.4 is "switch-point continuity" — spans `θ ∈ [10⁻⁴, 1]`.** Above it the
only records are `theta:1e0`'s 64 over `[1, π − 0.1)`: 21 of `coeff_d`'s 1710 binary64 records lie
in `θ ∈ [1.0, 1.4]`. Densifying that region is the obvious fix, and measuring it first showed it is
the wrong one.

### CO.12 holds at the lifted switch, and its right-hand side cannot be sampled

`d`'s two arms at the lifted switch against a 60-digit reference, each evaluated exactly as
`coeffs/kernel.rs` writes it:

| at `z` = | exact arm | series arm | sum | the jump |
|---|---|---|---|---|
| 1.5399 (lifted) | **19.178 u** | 2.926 u | **22.104 u** | **22.104 u** |
| 1.0 (committed) | 1.766 u | 0.215 u | 1.981 u | 1.551 u |

**The jump equals the sum of the arms' true errors to three decimals: CO.12 is satisfied exactly and
the lifted kernel is correct.** What failed is the test's estimate of the right-hand side, by 3.5×.

Sampling more finely cannot fix it, because the exact arm's error is a **sawtooth**. Over `θ` within
±0.4% of the lifted switch, 81 points, it swings from 0.098 u to 19.178 u — a factor of **195** —
and around the committed switch from 0.264 u to 23.913 u, a factor of **90**. The numerator
`θ² − 4 sin²(θ/2)` is a cancelling difference, so its rounding turns on where the operands' bits
fall and not smoothly on `θ`. A two-record sample of that is not a bound at 64 points per decade,
at 200, or at 2000.

**So `branch_continuity` passes at the committed switches by luck rather than by construction** —
`θ = 1` sits near the low end of a 90× range and the brackets happened to sample above it. That is
true today, before this record changes anything.

The blocker is therefore neither the corpus nor the maths: CO.12's right-hand side is a point sample
of an oscillating quantity. Loosening the test instead was considered and rejected — the slack needed
is 3.5×, which would leave it unable to fail anything CO.12 cares about.

**A per-stratum maximum was the next reading and it fails on both arms**, which measuring it
settled too:

| | at the switch | max over the stratum, the arm's own side |
|---|---|---|
| series arm | 2.926 u | **6.3 × 10⁶ u** (at `θ = 3.04`, where the 8-term series has diverged) |
| exact arm | 19.178 u | 19.722 u (over `θ ∈ [1, 1.241)`, a **2.8%** margin) |

The series arm's is vacuous: a bound of 6.3 million `u` cannot fail anything. The exact arm's covers
the switch by 2.8%, which is the same luck in a thinner form — the sawtooth peak near `θ = 1`
happens to sit just above the value at the switch, and nothing makes it so.

**What is sound is a reference at the switch**, and it is not circular. `conformance/sweeps/thresholds.csv`
already carries the chosen switch as a decimal, so the chain runs sweep → switch → reference →
test, and the sweep never reads the reference. A generator artifact holding each coefficient's value
and derivative at its committed switch, to 30 digits, the way `coeff_series` holds its series, lets
`branch_continuity` compute both arms there and compare each against the truth — no sampling, no
window, no slack. That is the measurement this section is made of: done by hand in mpmath at 60
digits, it gives 19.178 + 2.926 = 22.104 against a jump of 22.104.

## Decision

Settled under `0040` item 5; the measurements above are its basis.

1. **`PHASE1.md` §6's grid should span the domain and stop strictly below it**, not
   `[10⁻⁸, 1]`. The branch variable's range is what `NUMERICS.md` §12 allows, `z ∈ [0, π²]`; a grid
   that stops at `z = 1` excludes the optimum for three of eight binary64 coefficients and six of
   seven at binary32, and one that stops *above* `π²` lets the sweep retire the exact arm, which
   measured worse (config D). The last point should be the largest grid point below `π²`, and the
   construction extended with the same integer root so a point stays exact at both precisions
   (`0016` item 3).
2. **A boundary optimum should be an error, not a result.** A sweep that returns the grid's first or
   last point has not searched; `search.rs` already reports the chosen terms one grid point either
   side of the choice, so it knows. It should say so, and `just thresholds` should fail rather than
   write a switch that sits on an end.
3. **The `m ≤ 8` cap should go up with the grid, because the two are jointly binding.** Neither
   alone reaches below 1846 u; together they reach 94.5 u. The cap's natural value is the corpus's
   series length, sixteen (`PHASE1.md` §4.3), which is then the binding limit for `c` and the next
   thing to measure.
4. **Config E is the recommended reading of the two limits** — `z < π²`, `m ≤ 16` — because it is
   the only one that improves accuracy without a trade: five domination failures fixed, none broken,
   no row above 100 u, 67 rows better against 14 worse at `p = 1.9 × 10⁻⁹`. Config D's extra 1.9–2.6×
   on the SE_N(3) Jacobians costs 15 broken verdicts and is a separate decision.
5. **And §6 should search two switches, not one**, because a one-switch config E is 2.8× slower at
   near-identity `θ` — `PHASE3.md` §11's stated priority — while two switches are **faster than the
   committed kernel at both ends** (1.39× at near-identity, 1.72× above every switch) and cost one
   branch of 0.1 ns. The term count a series arm needs is set by the largest `z` it serves and paid
   by the smallest; one switch makes every input pay the hardest one's price, and at identical
   accuracy that price is **2.08× over the catalogue**. This is the third way the search space is
   too small, after the ceiling and the cap, and it is the one that changes the generated file's
   *shape* rather than its numbers.
6. **The second switch is searched after the objective, not with it**, because it cannot change the
   objective (the maximum is already attained by the long arm). Stage 1 is today's search unchanged;
   stage 2 is the cheapest arm that holds stage 1's objective, `O(m · n)` over a prefix maximum. §6
   should say so, since a joint search returns a degenerate answer and a reader would read that as
   "a second switch buys nothing".
7. **The grid keeps its decade structure and the domain bound moves into the selection rule.** The
   grid spans `z ∈ [10⁻¹⁶, 10]` by the same integer root at 64 points per decade, and the sweep
   **refuses** a switch at or above `π²`. A grid is a geometric object and a domain is a semantic
   constraint; putting the bound in the rule keeps
   `the_grid_has_64_points_per_decade_and_both_ends` intact and means a later domain change does not
   reshape the grid. Measured as "stop the grid at 9.66"; the same selection, better factored.
8. **A boundary optimum is an error at a grid end and a reported fact at the domain bound.** A
   choice on the grid's first or last point means the search space is too small and fails the run
   (item 2); a choice on the last point below `π²` is the domain binding, which `k`, `b` and `d` do
   under config E, and is printed rather than failed.
9. **`generated.rs` reports its feasibility, not only its objective**: the grid span, the term cap,
   and whether the choice touched either. The objective alone is what let eight switches sit on a
   wall for the file's whole existence.
10. **A second arm belongs to `Switch`, and stage 2 minimises the corpus-weighted term count.** The
    generated object stays one per coefficient, because `0004` makes `generated.rs` the single source
    and a hand-nested branch would be typed structure around generated numbers. The cost model is
    **one term, one unit** — linear, which the measured 0.33–0.50 ns per term supports — and the
    weighting is the corpus's own records. A *timed* cost model would make the generated file a
    function of the host and break D16; a term count keeps it a function of the corpus, as `0004`
    requires. The 90% share above is how the result was described, not the rule.
11. **The exact arm is not retired near `π`** (config D is not taken). It is the only reading that
    worsens rows — 89 of them — and the case for it rests on group *maxima*, which `0038` establishes
    must not be read from domination counts. Config E is strictly better on every measure available
    and is reversible; D waits for the paired comparison `0038` proposes.

## Rationale

The ceiling is not a numerics mistake, it is a search-space mistake, and those are invisible from
inside the result: `generated.rs` reports its objective and not its feasibility, so eight switches
have been sitting on a wall for as long as the file has existed, each one reported as *the*
generated choice. `0004` is right that a switch should be generated rather than typed, and that is
exactly why the generator's box has to be checked — a typed constant is at least visibly a
judgement, while a boundary optimum reads as a measurement.

The measurement is reported pooled as well as per stratum because `0038` (draft) was written two
hours earlier about this same mistake in a different place: 3 fixed domination failures out of 751
rows is the kind of number that reads as noise, and 37 better against 6 worse at p = 1.6 × 10⁻⁶ is
the same change seen properly.

## Consequences

- **The largest error in the shipped library falls by 7.15×** and the corpus's worst row moves from
  `coeff_e` to `coeff_c` (7395.5 u, untouched). The headline is per coefficient; the corpus-wide
  maximum improves only 1.17×.
- `coefficients.md` CO.14(ii) says a `jacobians_match_dual_*` tolerance (`PHASE3.md` §9) is
  `10²`–`10⁴ u` for `b, c, d, e`. Those tolerances have not been written yet, and `d` and `e`'s
  halve and fall 7× before they are.
- **At binary32 the ceiling still binds after one decade**: `d` and `e` choose the new top,
  `z = 10`, which is above `π²`. A switch above the domain means the series arm covers every input
  and the exact arm is unreachable at that precision — a legitimate outcome, and one §6 should say
  out loud rather than leave as an artifact of where the grid happens to stop.
- The `straddle-0.9` bench fixture (`crates/helicoid/benches/coeffs.rs`) straddles `a`'s and `b`'s
  switches, and neither moves, so it keeps its meaning. `k`, `d` and `e` move and no fixture names
  them.
- Regenerating is one PR that edits §6, extends the grid, and reruns both sweep targets; nothing is
  blessed, since no baseline exists. **Hazard:** changing the search space can make the tool see the
  series as stale, and it then writes placeholder switches and fails, by design — the second run is
  the real one (`justfile`, `thresholds`).
- D16 is unaffected: the grid points are integer-rooted and identical on every target, as the
  existing ones are. A coefficient whose series arm covers the domain calls no transcendental at
  all, which makes D16 trivially true for it — and is a cost question, not an accuracy one
  (open question 6).
- **The sweep's objective is a scalar max per coefficient and cannot see a per-stratum trade.**
  Config D improves that max by 3278× for `e` while making 89 rows worse, and the objective reports
  only the improvement. `0006`'s no-regress bar is exactly the instrument that would catch it, and
  it has no baseline yet (`PHASE3.md` §10), so nothing today stops a sweep from buying its objective
  with other strata. This is `0038`'s finding arriving by a third route.

The two limits are split because only one of them is free. **The grid lift alone is strictly
dominating at binary64**: an input between the old and the new switch moves from the exact arm
(8.07–11.46 ns) onto the 8-term series (2.41 ns) and every other input keeps its arm and term count,
so no input gets slower and three coefficients get 2.1–7.2× more accurate. **The term cap is not
free** — it is what costs 2.8× at near-identity — and the second arm is what pays for it, so they
land together or not at all. And step 0 comes first, because step 1 moves switches out of the
stratum that validates them.

0. **A reference at each switch, and `branch_continuity` reads it.**
   `conformance/generate` emits one record per `(coefficient, precision)` holding the value and the
   derivative at the switch `conformance/sweeps/thresholds.csv` names, at 30 digits — a derived
   artifact beside `coeff_series`, with its own `kind` in `MANIFEST.json` (`PHASE1.md` §4.2, §4.3).
   `branch_continuity` then evaluates both arms at the switch and compares each against the truth,
   so CO.12 is checked as stated: no sampling, no window, no slack. Verified by `branch_continuity`
   passing at the lifted switches with the recorded errors and no tolerance, and by the two numbers
   it prints matching this record's 19.178 u and 2.926 u. The four `at_switch_*` CSV columns stay,
   so the gap between a sample and the truth stays visible. **Owed, and it blocks step 1.**

   Two cheaper readings of this step were tried and refuted above, each by measuring it before
   building it: densifying `theta:dense` (the quantity oscillates 195× over 0.8% of `θ`, so no
   density makes a point sample a bound) and a per-stratum maximum (vacuous for the series arm at
   `6.3 × 10⁶ u`, and a 2.8% margin for the exact one).
1. **§6's grid span and selection rule, `m ≤ 8` unchanged**: the grid spans `z ∈ [10⁻¹⁶, 10]` by the
   same integer root (decision 7), a switch at or above `π²` is not a candidate, a choice on the
   grid's lowest point fails the run (decision 8), and `generated.rs`'s header carries the grid span
   and the term cap beside its objective (decision 9). Both targets regenerated. Verified by
   `the_grid_has_64_points_per_decade_and_both_ends` extended to the new end and the domain bound,
   `everything_tied_has_no_next_objective_and_an_empty_grid_is_an_error`,
   `the_committed_rows_are_the_documented_columns`, `a_constant_has_the_documented_shape`, the
   prior-comparison test and `just thresholds-check`. **Built once and reverted; owed after step 0.**
2. **The term cap to the corpus's series length, and `Switch`'s second arm, together** (decisions 6
   and 10): stage 2 after today's search, the arm chosen by corpus-weighted term count. This is the
   78–437× and the 2.08×, and it changes the generated file's shape and the shipped kernel's branch.
   **Owed, after step 1.**

## Open questions

None. Decisions 1–11 depend on no unresolved question; what the measurements suggest next is below
(`0040` item 2).

## Further work

0. **Should the per-switch reference be a corpus file or a sweep column?** Step 0 makes it a
   generator artifact, because the value is a function of the switch and the definition and nothing
   else. The alternative is the sweep writing it, which puts a 30-digit reference in a CSV the sweep
   also reads, and `0004` wants the generated file a function of the corpus alone.
1. **Does any other generated switch sit near a sawtooth peak?** The committed `θ = 1` sits near the
   low end of a 90× range, which is luck. Nothing checks it, and the sweep's objective cannot: it is
   a maximum over records, not a value at the switch.
1. **Does `PHASE1.md` §4.3's series length need to grow?** `c` takes all sixteen terms the corpus
   holds and is the corpus's worst row at 94.5 u; the catalogue's other seven are at 1.8–37 u. What
   the generator costs at, say, twenty-four terms is unmeasured. The Decision caps at the corpus's
   length, whatever it is, so it does not depend on this.
2. **Should the exact arm be retired near `π` after all?** Decision 11 says not now. The case is
   config D's 1.9–2.6× on the SE_N(3) Jacobian maxima, the only thing that moves a group maximum,
   against 15 broken verdicts all within 2.7 u absolute. Reopening it needs `0038`'s paired
   comparison, not another bar reading.
3. **Do the arm timings hold for the shipped groups?** `exp_coeffs` evaluates several coefficients in
   one `S::branch` and takes the exact closure as soon as one member passes its switch, so a group's
   switch is the smallest of its members' and a second switch changes which closure runs for a whole
   group. `benches/coeffs.rs` is the instrument and the comparison needs the regenerated file.
4. **What weighting should stage 2 use when a consumer supplies one?** Decision 10 takes the
   corpus's records, which `error-analysis.md` EA.21 made log-uniform across decades on purpose and
   which is a stand-in for a usage distribution, not one. locus-tag is the first consumer
   (`PHASE5.md` owns that migration) and
   a pose stream's `θ` distribution is nothing like flat.
