# 0031: What the cubic port inherits from omnisac

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** the PR that lands its *Decision*, with `0056`'s corpus id as the instrument.
Carried over from [`0022`](./0022-real-owes-acos-and-cos.md)'s open question 2.

## Context

`PHASE2.md` §6 is **NORMATIVE** that `solve_cubic` is "omnisac's algorithm as is". Testing the port
against planted roots found four places where the algorithm is not backward stable or answers the
wrong question. All four are omnisac's too; all four are already stated in the rustdoc `# Domain`
and pinned by tests (`a_repeated_root_can_be_dropped`, `underflow_is_pinned_at_both_precisions`,
`thresholds_sit_at_their_table_values`), so none is a surprise and none is a regression. The
question this record exists to answer is **which of them we stop inheriting**, and in what order.

Fixing any of them amends a NORMATIVE "as is" clause and moves golden rows, which is why it is a
record and not a commit (`0022` plan step 2: *a `NUMERICS.md` edit and a record for the formula,
then code*).

Measurements below: `mpmath` at 60 digits as truth, binary64, the same arithmetic the port performs;
`mp.polyroots` with `extraprec = 300`. The error figures are relative errors of the root.

## L1 — the one-real-root arm cancels

`t = cbrt(h + s) + cbrt(h - s)` with `h = -q/2` and `s = sqrt(disc)`. When `|p|³ ≪ q²`, `s -> |h|`
and `h - s` cancels; at the limit it is exactly `+0` and the arm returns `cbrt(2h)`, the root of
`t³ + q` — the `p` term silently gone.

`solve_cubic(1, 0, 1e-5, -1)`: `h = 0.5`, `s = 0.5`, `h - s = 0` exactly, and the arm returns
**exactly `1.0`** for a root of `0.999996666666666679`. Relative error `3.33e-6`, about `1.5e10` u.

| `p` in `x³ + p x - 1` | current | paired (`u v = -p/3`) |
|---|---|---|
| `1e-1` | 1.97e-14 | 1.63e-16 |
| `1e-3` | 9.20e-12 | 3.04e-17 |
| `1e-5` | **3.33e-6** | 3.42e-17 |
| `1e-8` | 3.33e-9 | 5.38e-17 |
| `1e-12` | 3.33e-13 | 4.44e-17 |

The pairing is `u = cbrt(h + copysign(s, h))`, the sum that cannot cancel, then `t = u - p/(3u)`
from `u v = -p/3`. `Real` already has `copysign`.

**Over 16 010 random one-real-root cubics the paired form is not uniformly better**: worse on 2103
rows, better on 4401, identical on the rest, and 40.6% of values move at all. It **dominates on the
max**, 4.67e-10 → 2.63e-12, which is the bar this workspace uses
([`0006`](./0006-the-instrument-comes-first.md): per function, per stratum, on the max, never a
mean). It also removes a `Dual` derivative of `-inf` at a simple root — wherever `h - s` is exactly
0 today, the derivative is `-inf` although the true one is finite, so an optimizer differentiating
through the arm gets a non-finite Jacobian on a smooth input.

Safety under a lane mask is preserved: the arm already evaluates `S::select(one_root, disc, one)`, so
an inactive lane has `s >= 1`, hence `|h + copysign(s, h)| >= 1` and `u` is never 0. In an active
lane `disc > band >= 0` strictly, so `s > 0` and `u != 0`. No new unsafe argument.

**Recommendation: fix, first.** Largest error of the four, a wrong derivative as well as a wrong
value, a one-line formula, no mask change — so no root count moves, unlike every other item here.

## L2 — both summands of `disc` underflow, and `p < 0` then forces the wrong arm

`three_roots` is `disc < -band` **or** `p < 0`. Where `p³/27` and `q²/4` both underflow to `+0`,
`disc` and `band` read `+0`, `disc < -band` is false, and a negative `p` alone sends the cubic to
the trigonometric arm whatever its roots.

`solve_cubic(1f32, 0, -1e-16, 1e-24)` has **one** real root — exactly, `q²/4 + p³/27` is
`2.13e-49 > 0` — and comes back with **three valid slots**, none of them a root. That is a wrong
answer, not an inaccurate one, and it is the only item here that is. The threshold is a depressed
root magnitude of about `4e-8` at `f32` and `2e-54` at `f64`.

`underflow_is_pinned_at_both_precisions` asserts the current mask on exactly this row. It pins the
documented behaviour; the rustdoc `# Domain` says so in the same words. It is not a test that
blesses an undetected bug — but it is a test that will have to change.

Fix: homogenise. The depressed cubic's roots scale as `m` when `p` scales as `m²` and `q` as `m³`,
so with `m` a power of two near `max(sqrt|p|, cbrt|q|)`, `P = p/m²` and `Q = q/m³` are `O(1)`,
neither summand under- or overflows, the arm selection is made on representable quantities, and
`t = m T`. Powers of two keep it exact, so nothing is lost to the scaling itself.

**Recommendation: fix, second.** It is the only wrong *answer*. It is a larger change than L1 and it
moves the arm selection, so root counts move and every golden row must be regenerated.

## L3 — the leading-coefficient floor rejects well-conditioned cubics

`is_cubic` requires `|a| > max(scale, 1) · 2^7 u`, `scale` being the largest `|coefficient|`. The
floor exists only to stop `1/a` blowing up, but it is relative to the **largest coefficient**, so it
caps the span of the coefficients at `1/tol`: `2^17` at `f32`, `2^46` at `f64`.

`solve_cubic(1f32, -102, -9799, 999900)` is `(x - 100)(x - 101)(x + 99)`, three well-separated
roots; `scale = 999900`, the floor is `7.63`, `|a| = 1`, and the mask is `[false; 3]` — the empty set.
`x³ - 10⁶ x` (roots `0, ±1000`) likewise. `f64` accepts both, so the two precisions disagree on root
**count**, not only on ulps. The rustdoc `# Domain` names both examples.

Fix: guard the quantities the floor is a proxy for. `B`, `C`, `D` are `b/a`, `c/a`, `d/a`; requiring
each to be finite is scale-free, exact, and is what the floor was approximating. `0017` decision 3
settled the *tolerance values*, not the predicate they are used in, so this is not a reopening of
`0017`.

**Recommendation: fix, third**, and note it is the one item where `f32` and `f64` disagree
qualitatively, which matters for D16's cross-target story even though both are deterministic.

## L4 — the discriminant band is relative to the wrong thing

`band = max(q²/4, |p³/27|) · 2^13 u` is relative to `disc`'s own two summands, not to the rounding
error incurred in forming `p` and `q`. Where `p` and `q` cancel — roots close together and far from
the origin — the real error exceeds the band, `disc` reads positive, and the one-real-root arm
returns a single root.

`(x - 1.09375)² (x - 0.921875)` at `f32`, every coefficient exact: `disc` is mathematically 0, reads
positive, and the double root is dropped. `a_repeated_root_can_be_dropped` pins it.

The correct scale is already computed — in the test file. `cubic_tests::disc_noise` propagates the
rounding error of forming `p` and `q` (`u·pw`, `u·qw` with `pw`, `qw` the magnitudes of the terms)
and `regular` uses `noise <= band` to decide which planted rows the error model covers. Moving that
model into the solver is the fix, and it would also shrink the carve-out that currently excludes
whole families of planted rows from asserting anything.

Cost: the solver must track `pw` and `qw`, which is extra work in the prologue every call pays,
including the "not a cubic" lane. That trade needs a measurement, and no document prescribes a
benchmark for `helicoid-linalg` (`PHASE1.md` §9 scopes `criterion` to `crates/helicoid/benches` from
Phase 3), so the harness is owed before this item can be judged.

**Recommendation: fix last, and only with a bench.** It is an accuracy limit, not a wrong answer, and
it is the only one of the four whose fix is not obviously free.

## L5 — `1/a` then three multiplies

`inv_a = 1/a` and `B = b·inv_a` adds a rounding to each of `B`, `C`, `D` against `b/a`, to save two
divisions in a function that spends four or five transcendentals. Measured: the depressed
coefficients differ from the three-division form on **21.2%** of 16 010 random cubics, and `q`
amplifies the difference threefold through `2B³/27`.

The review that raised it also notes `1/a` is **subnormal** when `|a|` is within a factor of ~8 of
the format maximum, losing about three bits outright; confirmed, `a = 1.7e308` gives
`1/a = 5.88e-309` and `3/a != 3·(1/a)`. Whether such an `a` survives `is_cubic` depends on L3's
predicate, so the two interact.

**Recommendation: fold into whichever of L2 or L3 lands first**, since both regenerate every golden
row anyway and this change alone would not justify doing so. On its own the trade is defensible:
accuracy-over-speed is this crate's stated priority, but two divisions against 1 ulp per coefficient
is not the lopsided call L1 is.

## Consequences

- Every item but L1 moves `solve_cubic`'s **mask**, so `bits_are_omnisacs`, the mpmath strata, the
  planted-root model and the three pinning tests all change. `PHASE2.md` §9's second clause already
  licenses that ("a polynomial between the old and the new tolerance can change its root count").
- `PHASE2.md` §6's "omnisac's algorithm as is" needs amending for any of them, and §6's `# Domain`
  paragraph and the rustdoc shrink as each lands. `NUMERICS.md` gains the paired Cardano form for
  L1 and the scaling for L2.
- The port stays checkable only if the reference moves with it: once the arms differ from omnisac
  deliberately, `GOLDEN` stops being "omnisac's bits" for those arms and must be regenerated from
  this implementation, with mpmath as the accuracy oracle rather than omnisac. That is a change of
  kind, and it is the real cost of this record — it should be paid once, for all the items that
  land, not four times.
- `0022` is independent: it changes which `libm` functions the trigonometric arm calls, not any
  arm's formula, and can land before or after.

## Decision

`0056` gave `solve_cubic` a corpus id, so each item below is measured on it, per stratum, at both
precisions, against nalgebra's companion-matrix roots (`0056` decision 4). `NUMERICS.md` §16 states
the formulas.

1. **L1 lands, and goes one step further.** The one-real-root arm pairs its cube roots: `w` is the
   cube root whose radicand cannot cancel and `v = -(p/3)/w`. Where `p <= 0`, `w` and `v` share a
   sign and the root is `w + v`. **Where `p > 0`, `w + v` cancels** whenever the root is small
   against the cubic's scale (`x³ + x + q`, root near `-q`). That regime is not in L1 above, and
   the unpaired form cancels there too. There the root is `-q / (w² + p/3 + v²)`, the same root from
   `w³ + v³ = -q` over a sum of positive terms. Over 20 000 one-real-root cubics with `p` and `q`
   log-uniform in `[1e-8, 1e8]`, the pairing alone reads **8.7e19 `u`** on the `p > 0` side and
   the quotient 10.7; on the `p < 0` side the pairing reads 5.1 and the quotient 9.8. The arm takes
   each where it is the better, through a lazy branch on the sign of `p`. A stratum,
   `cubic:one-real-small-root`, now draws the regime.
2. **L2 lands as homogenisation of the monic cubic**, not of the depressed one. Once L3 removes the
   floor, `B` can be as large as the format allows and `B²` overflows before `p` exists. So the
   power-of-two scaling `x = m y` divides `B`, `C`, `D` by `m`, `m²`, `m³`. Only over- and
   underflow of `p³` and `q²` has to be prevented, not a unit scale: every tolerance compares
   quantities of one degree in the scale, and the triple-root test is reached only when both
   summands of `disc` are zero or underflow. So `m = 1` inside a window, `[2^-128, 2^128)` (`f64`)
   or `[2^-8, 2^16)` (`f32`), tested by comparisons, and a cold, out-of-line ladder of two-sided
   power-of-two steps brings the rest in. No `sqrt` or `cbrt` is spent on the scale.
3. **L3 lands as written.** A cubic is `a != 0` with `1/a`, `B`, `C`, `D` finite.
4. **L5 does not land**, measured against. The three divisions cost about 20 ns of a 20–60 ns
   solver once L2 is in, in a binary that times both forms side by side. On the corpus they win
   `cubic:triple@f32` (0 against `6.8e4` `u`: `1/a` is inexact for a planted `a`) and lose
   `cubic:near-double-1e-4` (`3.5e5` against `1.1e5`). `1/a` stays, and with it the subnormal `a`
   whose reciprocal overflows, which is "not a cubic" (`# Domain`).
5. **L4 stays a limit**, documented and pinned (`a_repeated_root_can_be_dropped`). It owes a bench
   to judge its prologue cost, and `benches/linalg.rs` is now that bench.

## Measured

`max_u` per stratum, the shipped subject against the committed corpus (`0056`), before and after,
with nalgebra's binary64 row:

| stratum | before | after | nalgebra |
|---|---|---|---|
| `cubic:one-real` | 1815 | **1.70** | 1551 |
| `cubic:one-real-p-small` | 1.8e10 | **1.71** | 3.03 |
| `cubic:one-real-small-root` | 0.65 | **0.0055** | 0.0056 |
| `cubic:coeff-scale-down` | 9.0e15 (no root) | **51.2** | 9.3e4 |
| `cubic:one-real-p-small@f32` | 3.0e4 | **1.39** | — |
| `cubic:coeff-scale-down@f32` | 1.7e7 (no root) | **77.9** | — |
| `cubic:one-real@f32` | 638 | **2.46** | — |
| every other stratum | — | equal | — |

`solve_cubic` now dominates nalgebra on every binary64 stratum. Before, it lost three. The
small-root stratum reads below 1 `u` before too, because the root-set metric is relative to `‖Z‖`,
which the complex pair dominates. The quotient's gain is in the small root's *relative* error, the
`8.7e19 u` above, which no stratum of `0056` scores.

Latency, `benches/linalg.rs`, `bench-gate --against` a baseline built in this tree. Core 5 pinned,
with no other process above 20% CPU during the run; concurrent A/A floor median 0.7%:

| row | `f64` | `f32` |
|---|---|---|
| `one-real` | 45.2 → 46.0 ns (1.014) | 27.7 → 34.6 ns (1.25) |
| `distinct` | 34.6 → 38.6 ns (1.11) | 22.0 → 27.2 ns (1.24) |
| `triple` | 14.3 → 19.3 ns (1.30) | 12.1 → 16.5 ns (1.37) |

Two traps cost most of the first attempt, both now in the code's comments:
- **`homogenise` left to `#[inline]`** was a call returning four values through memory, 10–33 ns.
- **A bench wrapper taking `[S; 4]` by value** measured a store-to-load forwarding stall, 10–20 ns,
  that moved with the routine's own load order.

## Open questions

None. Question 1 is `0056`'s corpus id. Question 2: the bar is max-domination per stratum, plus
`0056`'s committed rows, which are per-row no-regress in both directions. Question 3: the constants
keep their values, since on the homogenised cubic they compare quantities of one degree in the
scale. Question 4: two treatments; `eig3` is `0023`'s.

## Former open questions

1. **Does `solve_cubic` get a corpus id and a conformance subject?** `PHASE2.md` §6's status row
   says "no corpus id and no conformance subject yet". Every item above is an accuracy claim, and
   the workspace's instrument for accuracy claims is the corpus and the envelope
   ([`0006`](./0006-the-instrument-comes-first.md)). Fixing these behind inline mpmath fixtures
   rather than corpus strata is the thing `0006` exists to prevent. **This may be the first step,
   before any of L1–L5.**
2. L1's recommendation rests on domination on the max, with 2103 of 16 010 rows getting worse. Is
   max-domination the right bar for a root-finder, or does a per-row no-regress clause belong here?
3. L2's scaling changes the depressed cubic the solver actually factors, so every tolerance in the
   `# Tolerances` table is then relative to an `O(1)` quantity. Do the three constants keep their
   values, or does the sweep that `0017` decision 3 declined become necessary?
4. Which of these, if any, blocks `eig3` and `svd3`? **Not by sharing `solve_cubic`** — `eig3`
   (on `linalg/eig3`) has its own closed form and imports only `cubic::{acos, pi}`, so none of
   L1–L5 reaches it and L2's wrong root count does not become a wrong eigenvalue count. What they
   *do* share is the trigonometric shape and therefore `0022`'s two methods. The open question is
   narrower: `eig3`'s rustdoc records the same kind of limit as L4 (not backward stable near a
   double eigenvalue, `sqrt(u)` in the angle) and `0023` (draft) carries its own list. Do the two
   get one treatment or two?
