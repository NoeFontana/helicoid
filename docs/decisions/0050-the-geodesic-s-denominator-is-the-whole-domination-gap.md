# 0050: the geodesic's denominator is the whole domination gap

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #100 (the measurement and this record); steps 2 and 3 as noted below

> **Step 2 landed, and the predictions held.** `just envelope` **115 → 113**, the number this
> record's *Consequences* named, with `se3_geodesic/geo:generic` the only geodesic failure left.
> The corpus reads **1.6440 / 1.7382 / 1.6417** `u` at `geo:consecutive` / `geo:generic` /
> `geo:near-pi`, matching the measurement's figures to every digit it printed, and the `host-std`
> twin reads the same, so the override costs D16 nothing. Two law legs moved and they are one fact:
> SO(3)'s `t=1` leg fell from 7.160 `u` to **1.118**, `gerr`'s floor for two bitwise equal
> quaternions — the endpoint exactness, read by the laws without the corpus — and `twin` rose from
> 1.118 to **7.233**, because the leg now compares two genuinely different expressions instead of a
> call with itself. `symmetry` improved (8.951 → 6.013) and `velocity` and `right` degraded
> (6.505 → 7.107, 9.663 → 10.979); all are identity legs, and `0006` makes the corpus the bar.
> `Product<SO3, Rn<3>>` inherited the same pair of moves through its per-factor delegation, and its
> bit-for-bit slerp-and-lerp assertion still passes, which is *Further work* 3 answered.
>
> **Step 3 ran, and it contests decision 2.** The geodesic had **no bench row at all** — `PHASE4.md`
> §0.0 owed six and `benches/groups.rs` held none, which is why the first gate run could say only
> that nothing else regressed (168 benchmarks, none above its own floor). With the rows added and a
> baseline built in this tree from the pre-override `so3.rs`, `so3/geodesic` reads:
>
> | stratum | `f64` | `f32` |
> |---|---|---|
> | `near-identity-7.5e-8` | **1.5591x slower** (19.7 -> 30.7 ns) | 0.4886, **2.05x faster** |
> | `generic-1` | 0.8447, 1.18x faster | 0.8427, 1.19x faster |
> | `near-pi` | 0.9084, 1.10x faster | 0.9091, 1.10x faster |
>
> All twelve `se3`/`se23` rows sit at 0.99–1.01, which is the control: `SEn3` does not route through
> `SO3::geodesic`, so the rows measure what they claim, and **the SE(3) path where `tf_tree`'s
> 300 ns gate lives is untouched**. But `cargo xtask bench-gate` **exits 1** on the one row — at
> least 1.5525 in every one of 6 pairs against a floor of 0.0044 — so this is resolvable, not noise,
> and `PHASE1.md` §9's rule is per benchmark.
>
> Decision 2 weighed `0.072 u` of accuracy against a `0004` sweep and took one arm. It did not know
> about +11 ns at the stratum a transform tree lives in, because no row existed to say so. *Further
> work* 4 holds what the two-arm now looks like; the decision is not reversed here.

## Context

`so3_geodesic` loses two of its three strata to oracle #1 and `PHASE4.md` §0.0 records the figures:
`helicoid` **1.572** `u` at `geo:consecutive` against `tf_tree_math::slerp`'s 2.187, and **2.721**
/ **2.429** against **1.834** / **1.642** at `geo:generic` / `geo:near-pi`. Two of the project's
three remaining domination failures are those rows (the third is `se3_geodesic/geo:generic`). §0.0
also states why `SO3::geodesic` is *not* overridden: GE.14 proves the provided body is shortest-arc
slerp "in the better-conditioned spelling — `atan2` where slerp has `acos`, whose slope is infinite
where the quaternions are nearly equal".

Nothing in that was wrong, and it does not explain the losses. GE.14 proves three spellings
**equal**, so the split cannot be the function; and the conditioning claim is about *nearly equal*
quaternions, which is `geo:consecutive` — the one stratum `helicoid` **wins**. So the argument is
confirmed exactly where it applies and silent everywhere else, and
[`0049`](./0049-the-boundary-is-what-removes-a-way-to-be-wrong.md) *Further work* 1 misread it as a
prediction about `geo:near-pi`, where `acos`'s argument `q_0\cdot q_1` is near **zero** and `acos`
is at its best-conditioned point. That sentence of `0049` is corrected by this record.

So the question was measured instead. `xtask/src/conformance/measure_geodesic.rs` scores six
spellings of the one function GE.14 equates, record by record over the committed 180-record corpus,
through `metric::rule("so3_geodesic")` — the same exact-rational metric the conformance runner
uses, so the figures below are directly comparable with the committed rows and needed no conversion.

**Protocol.** `cargo test -p xtask -- --ignored --nocapture measure_the_geodesic_spellings`,
binary64, maximum per stratum with the worst record's id. Route A is the control: it has to
reproduce the committed `helicoid` rows, and
`measure_geodesic::tests::the_shipped_route_reproduces_the_committed_rows` asserts it does, to
`5e-4`, in `just test`. No route uses `acos`: `Real` has none (`0022`) and D5 bans it on this path,
so a spelling with one could not ship and measuring it would answer nothing.

| spelling | `geo:consecutive` | `geo:generic` | `geo:near-pi` |
|---|---|---|---|
| **A** shipped, `q₀Exp(t Log(q₀*q₁))` | **1.572** (#4) | 2.721 (#104) | 2.429 (#139) |
| **B** GE.14 grouped, `q₀(cos tα, ϖ_t v)` | 1.572 (#4) | 2.502 (#79) | 3.253 (#139) |
| **C** blend ÷ `‖v‖`, α from the product | 2.279 (#23) | 2.442 (#60) | 2.545 (#120) |
| **D** blend ÷ `‖v‖`, α from the 4-dot | 5.05e7 (#5) | 2.481 (#90) | **1.642** (#139) |
| **E** blend ÷ `sin α`, α from the product | 1.644 (#11) | **1.738** (#73) | **1.642** (#139) |
| **G** blend ÷ `sin α`, α from the chord | 1.644 (#11) | **1.738** (#73) | 1.726 (#153) |
| oracle `tf_tree_math::slerp` | 2.187 | 1.834 | 1.642 |

Four things that table settles, none of which was known before it:

1. **The gap is the denominator, and nothing else.** C and E differ in one token — `‖v‖` against
   `sin α` recomputed from the same `α` — and that token is worth **2.442 → 1.738** at
   `geo:generic` and **2.545 → 1.642** at `geo:near-pi`. Every other difference in the table is
   worth a few hundredths.
2. **Because it is the endpoints.** At `t = 0` the left weight is `sin α / den`, which is exactly
   one only when `den` *is* that same `sin α`; at `t = 1` the mirror. Dividing by `‖v‖` mixes a
   `sqrt`-derived denominator with `sin`-derived numerators, so the endpoint is off by that `ε`. The
   corpus says so in the record ids: every losing maximum sits on or beside an endpoint — C's
   `geo:generic` worst is #60 at `t = 0`, D's is #90 at `t = 0`, C's `geo:near-pi` worst is #120 at
   `t = 0`, and even **A's** `geo:generic` worst is #104 at `t = 1 − 10⁻⁹` — while E's worsts move
   to the interior (#73 at `t = 0.5`, #139 at `t = 0.727`). E is exact at **both** endpoints where
   the shipped body is exact only at `t = 0`, which `PHASE4.md` §0.0 already records one-sidedly,
   and `measure_geodesic::tests::the_recomputed_denominator_is_what_makes_both_endpoints_exact`
   pins both halves over every endpoint record, including that the shipped route still fails at
   `t = 1`.
3. **`PHASE4.md` §1.2's screw twin cannot fix these rows.** B is the rotation part of GE.12, and it
   is **bit-identical to A** at `geo:consecutive` (same `1.572`, same record #4) — as it must be,
   since `k_t·t·r·v` and `ϖ_t·v` are the same product re-associated — and at `geo:near-pi` it is
   **worse**, 3.253 against 2.429. So the twin's rotation path *is* the shipped path, and the
   remaining `so3_geodesic` failures were never its to close. `0049` *Further work* 1 said so from
   the algebra; this is the measurement.
4. **D identifies the oracle.** It reproduces `tf_tree_math::slerp`'s `geo:near-pi` row to the digit
   and on the same record (#139), and it explodes to `5.05e7 u` at `geo:consecutive` — which is the
   `1 − d²` cancellation, and therefore the reason that function carries a normalized-LERP
   fallback. Reading it confirms the identification: its first line is `h = ½‖q_a − q_b‖²`, the
   cancellation-free `1 − |d|`, and its large arc is `acos`/`sin` with `sin(angle)` as the
   denominator — the same choice E makes.

## Decision

1. **`SO3::geodesic` is overridden with route E**, which `NUMERICS.md` gains as GE.14's right-hand
   expression made explicit:

   $$\gamma(q_0, q_1, t) = \frac{\sin((1-t)\alpha)\,q_0 + \sin(t\alpha)\,q_1}{\sin\alpha},\qquad
   (w, v) = q_0^* q_1 \text{ with } w \ge 0,\quad \alpha = \mathrm{atan2}(\lVert v\rVert, w),$$

   with the denominator **`sin α` recomputed from `α`**, never `‖v‖`, and `‖v‖ = 0` returning `q₀`
   — the arc is a point and `q₀` is both its endpoints. This is already GE.14, so it is not a
   formula change; what the `NUMERICS.md` edit adds is that the denominator's *spelling* is
   normative, because the measurement says it is the entire difference between dominating oracle #1
   and losing to it.

2. **One arm, no switch.** A two-arm routine — A below a switch, E above — would read
   `1.572 / 1.738 / 1.642`, the best cell of every column, and it is **not taken**. The whole gain
   over single-arm E is `0.072 u` on one stratum; the cost is a `0004` sweep, a generated constant,
   a second arm under D6's twin obligation and a second spelling against `API.md` R3. E alone
   dominates oracle #1 on all three strata, which is the bar `0006` sets, so the switch buys
   nothing the bar reads.

3. **Route G is recorded and not taken.** It ties E at two strata, needs **no Hamilton product**
   (its `α` comes from the chord `h = ½‖q₀ − q₁‖²`), and is behind at `geo:near-pi` — 1.726 against
   1.642 — because Sterbenz buys nothing there: the two quaternions are not within a factor of two
   of each other, so `h` is near 1 and `h(2 − h)` carries its own cancellation. It is the candidate
   to revisit **if and only if** a bench shows E's product is the cost, and `0006` ranks accuracy
   first.

4. **`PHASE4.md` §0.0's "not overridden, on purpose" is superseded**, and its conditioning sentence
   kept: `atan2` against `acos` is still the right reading *at `geo:consecutive`*, which is the one
   stratum where the shipped body stays the better spelling and where E costs 1.046× (1.644 against
   1.572). The claim that stops being true is that the provided body is the better spelling
   everywhere.

5. **`0049` *Further work* 1's reading is corrected**: GE.14 predicts a `helicoid` win where the
   quaternions are *nearly equal*, which is `geo:consecutive`, and the corpus confirms it. There was
   never a paradox at `geo:near-pi`; `acos`'s slope is infinite at argument `±1`, not at `0`.

6. **`se3_geodesic/geo:generic` is not closed by this** and is not claimed to be. `SEn3`'s geodesic
   is the provided body over `SEn3::rminus`/`rplus`, which does not route through `SO3::geodesic`, so
   overriding SO(3) moves no `se3_geodesic` row. Its translation block is coupled through
   $\mathsf V(\varphi)$ and a blend of two poses is not its geodesic, so the SE(3) analogue is
   `PHASE4.md` §1.2's screw form and owes a measurement of its own before a line of it is written.

## Rationale

The alternative was to take §1.2's screw twin as the fix, which is what the approved plan's H4
implies and what `0049` *Further work* 1 doubted on algebraic grounds. The measurement settles it
against: the twin's rotation part is the shipped expression re-associated, bit-identical on the
stratum where both are best and worse where it matters. Writing it first would have spent the
expensive item and left both rows failing.

The alternative to overriding at all was to accept two domination failures as a property of the
`Exp`/`Log` route and except them under `0046`. That is what `0046` exists for and it is wrong here:
an exception needs a reason the program cannot help, and a one-token change that dominates every
stratum is the opposite of one.

**Not taking the two-arm** is the decision most likely to be revisited, so its ground is stated
plainly: `0006` reads domination and no-regress, E satisfies both on every stratum, and the only
thing the extra arm improves is a figure no bar reads. `0047` is the precedent — a second arm is
admitted on evidence, not on a hundredth of a `u`.

## Consequences

- `SO3::geodesic` stops being a provided method and becomes an override, so `laws::geodesic`'s seven
  legs are measured against a different implementation on SO(3) and their bounds are re-recorded at
  twice the new worst. The `twin` leg becomes load-bearing in a way it was not: `reference::geodesic`
  is the provided body, so the leg now compares two genuinely different expressions rather than one
  with itself, which is what D6 asks of every §14 row.
- `se3_geodesic` reads `SEn3`'s own path and does **not** move, so after step 2 the geodesic failures
  go from three to one and the remaining one is SE(3)'s alone (decision 6). `just envelope` is
  expected to go 115 → 113; if it does not, the override has reached a row this record did not
  predict and that is a finding, not a tolerance.
- `t = 1` becomes exact on SO(3), which no document currently claims for any group. `PHASE4.md`
  §0.0's bit-exactness sentence gains the right endpoint for SO(3) and keeps its `-0.0` exception,
  which E has for the same reason the provided body does: a weight of exactly zero added to a
  `-0.0` component gives `+0.0`.
- One more `libm` call site on the hot path (`sin α` as well as the two weight sines), against two
  fewer Hamilton products and both coefficient evaluations. The direction is not obvious from the
  count, so it is benched, and `0006` would keep E on accuracy even if it is slower.
- `measure_geodesic.rs` is a permanent instrument with two non-ignored tests, so the next spelling
  is compared rather than argued.

## Implementation plan

1. The measurement and this record, with `0049` *Further work* 1 corrected and `PHASE4.md` §0.0
   marked superseded for the override — verified by `just lint`, `just test` (the control and the
   endpoint test), and the row in `decisions/README.md`.
2. `SO3::geodesic` as route E, generic over `S: Real` with `S::Mask` for the `‖v‖ = 0` arm; the
   `NUMERICS.md` §10 and GE.14 edits decision 1 names; `laws::GEODESIC_LEGS`'s bounds re-recorded at
   twice the measured worst over `10⁶` draws — verified by `just test`, `just conformance` and
   `just envelope` reaching **113**.
3. `benches/groups.rs`'s three `so3/geodesic` rows re-run against a baseline binary built in this
   tree — verified by `cargo xtask bench-gate --bench groups --against <binary>`, with the number
   recorded here whichever way it goes (`0006`: accuracy decides, latency is reported).

## Open questions

None.

## Further work

1. **Decision 2 is contested by latency, and the cheap two-arm needs no sweep.** The banner's step-3
   table is the evidence. What makes this different from the two-arm decision 2 declined: route **B**
   is **bit-identical to route A** at `geo:consecutive`, and route B's coefficients are
   `r = log_ratio(..)` and `(k_t, cos t\alpha) = exp_coeffs(..)` — already-swept switches. So "route
   A below the switch, the blend above" is a dispatch on **a generated switch that already exists**,
   not a new one: no sweep, no typed constant, `0004` satisfied by reading the predicate `coeffs`
   already computes rather than adding a number. That would read the best accuracy cell of every
   column (1.572 / 1.738 / 1.642) **and** the best latency of every row. What it costs is a
   `pub(crate)` predicate out of `coeffs`, a second arm under D6's twin obligation, a second
   spelling against `API.md` R3, and re-measurement at the switch boundary where the two arms meet —
   which is exactly where `0047` found a prefix rule can be 50224x wrong, so the boundary is measured
   and not assumed. Not authorised here: it reverses a decision this record took, so it is its own
   record.
2. **The SE(3) analogue, which is the last geodesic failure.** Decision 6 says why this record does
   not reach it. The measurement to run first is the same shape: score candidate spellings of
   `SEn3::geodesic`'s translation block against `se3_geodesic`, with `measure_geodesic.rs` extended
   rather than copied, before §1.2's twin is written.
3. **Binary32.** Every figure here is binary64, which is where the three failures are. A spelling
   with a `sin` in the denominator has one more cancellation path at 24 bits than at 53, and §1.2's
   twin will want the column anyway.
4. **Whether `Rn` and the products want anything.** `Rn`'s geodesic is lerp and exact; `Product`
   delegates per factor, so `Product<SO3, Rn<3>>` inherits decision 1 for free and its bit-for-bit
   slerp-and-lerp assertion in `product_tests::so3_r3` has to be re-read against the new body.
