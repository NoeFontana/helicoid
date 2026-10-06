# 0051: two arms on the switch the sweep already chose

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record and the two arms land together

> **Step 2 ran. The accuracy claim holds exactly; the latency claim holds in part, and the part that
> does not is written here rather than rounded off.**
>
> `just envelope` **113**, with `so3_geodesic` at **1.5721 / 1.7382 / 1.6417** `u` — the best cell of
> every column, oracle #1 dominated on all three strata. Every figure in the scan reproduced.
>
> Latency, against the pre-`0050` provided body with both binaries built in this tree
> (`bench-gate --bench groups --only 'so3/geodesic' --against <baseline>`; the id separator is `/`,
> not the `_` criterion's directories use):
>
> | `so3/geodesic` | `f64` | `f32` |
> |---|---|---|
> | `near-identity-7.5e-8` | 1.1162 (least 1.0319, floor 0.0814) | **1.1251** (least 1.1124, floor 0.0271) FAIL |
> | `generic-1` | 0.9139, 1.09x faster | 0.8889, 1.12x faster |
> | `near-pi` | 0.9797, 1.02x faster | 0.9684, 1.03x faster |
>
> Against `0050`'s single blend, which read 1.5591 / 0.4886, 0.8447 / 0.8427, 0.9084 / 0.9091, three
> things are true and the third is the one to keep in view:
>
> 1. **The `f64` regression is cut from 1.56x to 1.12x**, 19.7 -> 30.7 ns becoming 19.7 -> 22.0 ns.
>    Removing the predicate's division — `n² < T w²` instead of `n²/w² < T`, which `log_ratio` would
>    otherwise compute twice — is worth about a third of what was left.
> 2. **The wins above the switch shrink**, 0.8447 -> 0.9139 and 0.9084 -> 0.9797, because the
>    dispatch is paid on the blend's side too.
> 3. **At `f32` the single blend was better on every row** — 0.4886 against this 1.1251 at
>    `near-identity` — so the two arms are a `f64` improvement and a `f32` cost. *Further work* 1
>    holds the one-line fix and why it is not taken.
>
> And a finding about the gate rather than the routine: the `f64` `near-identity` row **failed at a
> 0.0030 floor and passed at a 0.0814 one**, same code, two runs differing only in how many
> benchmarks shared the machine. Its point estimate is 1.1162 either way. So that row is not
> resolvable on a busy host and a verdict on it needs `PHASE1.md` §9's quiet machine; this record
> does not read the pass as a pass.

## Context

[`0050`](./0050-the-geodesic-s-denominator-is-the-whole-domination-gap.md) overrode `SO3::geodesic`
with GE.14's blend and closed two of the project's three remaining domination failures, `envelope`
115 → 113. Its step 3 then found what no row had been able to say, because the geodesic had **no
bench row at all** — `PHASE4.md` §0.0 owed six and `benches/groups.rs` held none:

| `so3/geodesic` | `f64` | `f32` |
|---|---|---|
| `near-identity-7.5e-8` | **1.5591, 1.56× slower** (19.7 → 30.7 ns) | 0.4886, 2.05× faster |
| `generic-1` | 0.8447, 1.18× faster | 0.8427, 1.19× faster |
| `near-pi` | 0.9084, 1.10× faster | 0.9091, 1.10× faster |

`cargo xtask bench-gate` exits 1 on that row — at least 1.5525 in every one of 6 pairs against a
floor of 0.0044, so resolvable and not noise — and `PHASE1.md` §9's rule is per benchmark. All
twelve `se3`/`se23` rows sit at 0.99–1.01, the control that says the rows measure what they claim
and that the SE(3) path carrying `tf_tree`'s 300 ns gate is untouched.

The cause is not subtle. Below the switch the provided body's `log_ratio` and `exp_coeffs` take
their **series** arms — polynomials, no transcendental — where the blend pays `atan2` and three
`sin` whatever the angle. At `f32` those transcendentals are cheap enough that the blend wins
anyway; at `f64` they are not.

`0050` *Further work* 1 proposed the two-arm and asserted it needed no sweep, because route B is
bit-identical to route A at `geo:consecutive` and rides switches that already exist. **That
assertion was half right and the half that was wrong is the whole design**, which is why this
record exists rather than an amendment. Dispatching on `log_ratio`'s series/exact switch —
`n²/w² < 8.977×10⁻²`, `θ ≈ 0.58` — reads:

| stratum | on `R_F64.below` | `0050`'s one arm | oracle #1 |
|---|---|---|---|
| `geo:consecutive` | **1.5721** | 1.6440 | 2.187 |
| `geo:generic` | **2.5019** | 1.7382 | 1.834 |
| `geo:near-pi` | 1.6417 | 1.6417 | 1.642 |

`envelope` 113 → **114**: the stratum `0050` bought is given back. **50 of `geo:generic`'s 60
records sit below that switch** and take the provided body, and the provided body is worse there.
The provided body stays the *faster* arm up to `θ ≈ 0.58` but stops being the more *accurate* one
almost immediately above identity, so the two crossovers are decades apart and `0006` makes the
accuracy one bind.

So the threshold was measured instead of assumed, over the committed corpus, by
`measure_geodesic::scan` — `cargo test -p xtask -- --ignored --nocapture
scan_the_dispatch_threshold`, binary64, per-stratum max of "the provided body where `s < T`, the
blend elsewhere":

| `T` on `s = tan²α` | `geo:consecutive` | `geo:generic` | `geo:near-pi` | |
|---|---|---|---|---|
| `R_F64.below`, `θ ≈ 5.8×10⁻¹` | 1.5721 | 2.5019 | 1.6417 | **loses** |
| **`R_F64.short`, `θ ≈ 1.4×10⁻²`** | **1.5721** | **1.7382** | **1.6417** | dominates |
| `10⁻⁶`, `θ ≈ 2.0×10⁻³` | 1.5721 | 1.7382 | 1.6417 | dominates |
| `10⁻⁸`, `θ ≈ 2.0×10⁻⁴` | 1.5721 | 1.7382 | 1.6417 | dominates |
| `10⁻¹²`, `θ ≈ 2.0×10⁻⁶` | 1.5721 | 1.7382 | 1.6417 | dominates |
| `0`, the blend alone | 1.6440 | 1.7382 | 1.6417 | dominates |

Two readings that decide the design. **`r`'s second switch is already in the right place** —
`0047`'s short arm, four decades of `s` below the series/exact one — and it reads the best cell of
every column. And **every threshold at or below it reads identically**, so the boundary has at least
five decades of slack: this is a choice, but not a delicate one.

## Decision

1. **`SO3::geodesic` is two arms of the one function GE.14 equates**, dispatched on
   `log_ratio`'s **second** generated switch:

   - `s = n²/w² <` `short_below`, and `w > 0`: the provided body,
     `q₀ Exp(t Log(q₀* q₁))`, which is also `reference::geodesic`.
   - otherwise: GE.14's blend with the denominator `sin α` recomputed (`0050` decision 1,
     unchanged).

   The relative quaternion `q₀* q₁` is formed **once** and both arms read it, so the shared prefix
   is the dispatch's cost and not an extra one.

2. **No switch of its own, and the one it borrows is named as borrowed.** `0004` forbids typing a
   switch point and this needs none, but `r`'s second switch is chosen by the sweep to minimise a
   corpus-weighted term count (`0047`) — not to place this boundary. So the coupling is made
   explicit rather than silent: `coeffs::log_ratio_takes_short_arm` is the single definition of the
   predicate, its rustdoc carries the scan's finding, and
   `measure_geodesic::tests::the_shipped_dispatch_threshold_still_dominates` **fails** if a later
   sweep moves that number into the region the table says loses.

3. **`t ≥ 1` takes the blend on both sides of the switch.** Only the blend is exact there — its
   right weight is `sin(1·α)/sin α`, one number over itself — and `0050` shipped that bit-exactness
   with `PHASE4.md` §0.0 recording it. Giving it up at `geo:consecutive` would be backwards:
   consecutive keyframes are exactly where a query *at* the later one happens. One `le` buys it, and
   `t > 1` is extrapolation, which GE.14 covers and the blend computes, so the same arm serves it.

   The predicate reads `t`, which is safe here and would not be in an early return. Both arms are
   implementations of the same function, so a `Dual` takes the selected arm's derivative and both
   are right; returning the constant `x₁` would zero it, which is the trap `PHASE4.md` warns of and
   this is not.

4. **The `t = 0` arm is not special-cased**, because it does not need to be: `scale(0)` gives `±0`,
   `Exp` of that is exactly the identity quaternion, and multiplying by it is exact. The provided
   body is already bit-exact at `t = 0`, which `so3_tests::geodesic_at_zero_is_the_left_endpoint_bit_for_bit`
   has asserted since before `0050`.

5. **`0050` decision 2 is superseded and its *Further work* 1 corrected.** Decision 2 took one arm
   on the ground that the two-arm's whole gain was `0.072 u` on one stratum; it was taken before any
   bench row existed, so it could not weigh `+11 ns` at the stratum a transform tree lives in.
   *Further work* 1's claim that the existing switch would serve is wrong by four decades, and the
   corrected claim is narrower and measured: the existing **second** switch serves.

## Rationale

The alternatives, each refuted by a number above rather than by preference:

- **One arm, the blend** (`0050` as shipped) fails `PHASE1.md` §9 on a resolvable row. `0006` ranks
  accuracy first, which is why the blend shipped at all, but §9 is a separate gate and a failing
  gate is a failing gate.
- **One arm, the provided body** (before `0050`) loses two strata to oracle #1.
- **Two arms on the series/exact switch** loses `geo:generic` and takes `envelope` back to 114.
- **A swept threshold of its own** is what *Further work* 1 feared and the scan makes unnecessary:
  the sweep's objective is a corpus-weighted term count over a coefficient's own records, not a
  per-stratum max over a group routine's, so a geodesic entry would be a new kind of sweep for a
  boundary that already has five decades of slack. If the guard test ever fires, this is the answer
  to reach for — and it will be reached for with a table rather than an argument.

**Why the accuracy crossover is where it is**, since the design rests on it: the provided body's
error is the `Log` → `scale` → `Exp` round trip's roundings, which grow with `θ` from essentially
nothing at identity; the blend's is flat, because its weights are `O(1)` and its endpoints exact.
So the provided body wins only in the band where the round trip has not yet accumulated anything,
which is narrower than the band where its series arms are still cheap. Two different mechanisms,
two different boundaries, and no reason for them to coincide.

## Consequences

- `coeffs` gains one `pub(crate)` predicate. It is not a new number and not a new arm: it is the
  comparison `log_ratio` already makes, named so a caller can make the same one.
- **`symmetry` is this record's one cost**, and it is structural. GE.2(c) asks
  `γ(x₀, x₁, t) = γ(x₁, x₀, 1−t)`, so at `t = 0` the swapped call is at `t = 1` and takes the blend
  where the original takes the provided body — two arms differing by about 7 `u` being asked to
  agree. SO(3)'s leg goes 6.013 → **7.670** `u` at `f32` and `Product<SO3, Rn<3>>`'s 5.4233 →
  6.4564; both bounds are re-recorded at twice the new worst. No other leg moves.
- `left_invariance_at_three_translation_scales` returns to its pre-`0050` figures to the digit
  (5.8658 / 6.8112 / 3.8000), because its draws are near-identity and take the provided body again.
  That the figures come back *exactly* is itself a check on the dispatch: the arm below the switch
  **is** what shipped before.
- Both endpoints stay bit-exact, so `t=0` and `t=1` both read `gerr`'s 1.118 floor.
- The `twin` leg stays off the floor (7.213 `u`), because one arm is not `reference::geodesic` — D6
  satisfied on the §14 row `0050` added.
- Four decades of `s` separate the shipped threshold from the one that loses, and a test now stands
  in the gap. That is a dependency on the sweep that did not exist before, and it is written down.

## Implementation plan

1. This record; `coeffs::log_ratio_takes_short_arm`; the two arms; the re-recorded `symmetry`
   bounds; `measure_geodesic`'s scan and its two guard tests — verified by `just lint`,
   `just test`, `just conformance` and `just envelope` holding at **113**.
2. `cargo xtask bench-gate --bench groups --only geodesic --against <baseline built in this tree>`
   showing the `near-identity` regression gone and the `generic` / `near-pi` wins kept — the figure
   recorded here whichever way it goes.

## Open questions

None.

## Further work

1. **The `f32` arms should probably be one, and `Real` already has the lever.** Measured, the single
   blend beats the two arms on **every** `f32` row — 0.4886 against 1.1251 at `near-identity`,
   0.8427 against 0.8889 at `generic`, 0.9091 against 0.9684 at `near-pi` — and costs `0.072 u` on a
   stratum that has no `f32` corpus row to score it. `Real::PRECISION` is a `const`, so
   "blend alone at binary32, two arms at binary64" is a compile-time choice with no runtime test.
   **Not taken**, for two reasons worth stating rather than assuming: no consumer is on `f32`
   (`tf_tree` and locus-tag are both `f64`), so the 3.4 ns is nobody's yet; and an algorithm that
   differs by scalar type is a real departure from `0003`'s uniform body, which should answer to a
   consumer and not to a bench row. One line when one exists.
2. **Binary32 has no corpus evidence.** Every figure in the scan is binary64, because
   `so3_geodesic` has no `f32` stratum. The `f32` threshold is `R_F32.short`, `s < 1.4855×10⁻³`
   (`θ ≈ 7.7×10⁻²`), and the laws' `as_f32` legs are all that cover it. An `@f32` geodesic stratum
   (`0016`) would close it, and §1.2's twin will want one anyway.
3. **The `symmetry` leg's cost is removable in principle** by dispatching on a quantity symmetric
   in `t ↔ 1−t`, which the `t ≥ 1` clause is not. Whether that is worth a third predicate is a
   question for a record with a measurement, not this one.
4. `0050` *Further work* 1's successor: **the SE(3) analogue** is still the last geodesic failure,
   and `SEn3::geodesic` reaches none of this.
