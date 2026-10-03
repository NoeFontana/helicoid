# 0039: The sweep's grid stops below its own optimum

**Status:** draft
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

## Decision

Proposed, not settled.

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
4. **Config E is the recommended reading** — `z < π²`, `m ≤ 16` — because it is the only one that
   improves without a trade: five domination failures fixed, none broken, no row above 100 u, 67
   rows better against 14 worse at `p = 1.9 × 10⁻⁹`. Config D's extra 1.9–2.6× on the SE_N(3)
   Jacobians costs 15 broken verdicts and is a separate decision.

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

## Implementation plan

1. §6's grid span, the extended construction, and the boundary check of item 2 — verified by the
   existing `the_grid_has_64_points_per_decade_and_both_ends` extended to the new end, and by a new
   test that a chosen end-point switch fails the run. **Owed, blocked on open question 1.**
2. Both sweep targets regenerated, with the objectives above as the expected result, and
   `just thresholds-check` green. **Owed, after 1.**
3. The `m ≤ 8` cap measured the same way, for `c`. **Owed.**

## Open questions

1. ~~Where exactly should the grid stop?~~ **Measured: just below `π²`.** Stopping above it
   (config D) retires the exact arm and makes 89 rows worse. Open: whether the bound is written as
   "the largest grid point below `π²`" or by extending the decade structure and refusing the points
   above — the first is what was measured, the second keeps the existing construction's test.
2. ~~Is a switch above `π²` allowed to stand?~~ **No, on this corpus.** But config D's SE_N(3)
   Jacobian improvement (1.9–2.6×) is real and is the only thing that moves a group maximum, so
   **should the exact arm be retired near `π` anyway?** That is the one live trade here: 15 broken
   domination verdicts, all within 2.7 u absolute, for a halving of the Jacobian maxima.
3. ~~Does the `m ≤ 8` cap bind for `c`?~~ **Yes, and it binds jointly with the ceiling.** Open: `c`
   now chooses all sixteen terms the corpus holds, so does `PHASE1.md` §4.3's series length need to
   grow, and what does the generator cost at, say, 24 terms?
4. **Is a boundary optimum against the *domain* an error?** Decision item 2 says a choice on a grid
   end has not searched, but under config E three of eight switches sit on the last point below
   `π²` — a bound that is physical, not arbitrary. The check must exempt it, which means it has to
   distinguish the two ends.
5. Should the sweep report its *feasibility* as well as its objective — the grid span, the term cap,
   and whether the choice touched either? `generated.rs` carries the objective alone today, which is
   what let this stand.
6. Does retiring a coefficient's exact arm change its **cost**? Under config E the series arm runs
   for almost every input and calls no transcendental, where the exact arm calls `sqrt` and
   `sin_cos`; the series is longer. §6 says the series arm is the cheaper one and nothing here
   measured it. `PHASE1.md` §9's gate is what would, and `coeffs.rs`'s `straddle-0.9` fixture
   straddles `a` and `b`, whose switches both move under E.
