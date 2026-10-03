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

## Decision

Proposed, not settled.

1. **`PHASE1.md` §6's grid should span the domain, not `[10⁻⁸, 1]`.** The branch variable's range is
   what `NUMERICS.md` §12 allows, `z ∈ [0, π²]`; a grid that stops at `z = 1` excludes the optimum
   for three of eight binary64 coefficients and six of seven at binary32. The upper end should be
   stated as the domain's, and the construction extended with the same integer root so a point
   stays exact at both precisions (`0016` item 3).
2. **A boundary optimum should be an error, not a result.** A sweep that returns the grid's first or
   last point has not searched; `search.rs` already reports the chosen terms one grid point either
   side of the choice, so it knows. It should say so, and `just thresholds` should fail rather than
   write a switch that sits on an end.
3. **The `m ≤ 8` cap is a separate question and is not answered here.** `coeff_c` is the one
   coefficient the ceiling does not help, it is now the largest error in the corpus at 7395.5 u, and
   it uses 8 of the 8 terms it is allowed. The committed series holds 16.

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
  existing ones are.

## Implementation plan

1. §6's grid span, the extended construction, and the boundary check of item 2 — verified by the
   existing `the_grid_has_64_points_per_decade_and_both_ends` extended to the new end, and by a new
   test that a chosen end-point switch fails the run. **Owed, blocked on open question 1.**
2. Both sweep targets regenerated, with the objectives above as the expected result, and
   `just thresholds-check` green. **Owed, after 1.**
3. The `m ≤ 8` cap measured the same way, for `c`. **Owed.**

## Open questions

1. **Where exactly should the grid stop?** `z = π²` is the domain, and a grid point is not generally
   a nice number there; one decade (`z ≤ 10`) covers it and keeps the decade structure the existing
   construction and its test depend on. `z ≤ 10` is the cheap answer and is slightly outside the
   domain, which item 2's boundary check would then flag for binary32's `d` and `e` — correctly, but
   it means the two decisions interact.
2. Is a switch above `π²` allowed to stand? It makes the exact arm unreachable for that coefficient
   and precision, so the reference twin (D6) becomes the only thing exercising it, and
   `0003` item 3's safe-argument pattern is then carried by lane code alone.
3. Does the `m ≤ 8` cap bind for `c`? Its 7395.5 u is now the corpus's worst and `coefficients.md`
   CO.14 attributes it to a `θ⁻⁴` cancellation with `Ĉ' ≈ 4450`, which more terms would push to
   larger `θ` exactly as the ceiling did for `e`. A 16-term run was started and not completed.
4. Should the sweep report its *feasibility* as well as its objective — the grid span, the term cap,
   and whether the choice touched either? `generated.rs` carries the objective alone today, which is
   what let this stand.
