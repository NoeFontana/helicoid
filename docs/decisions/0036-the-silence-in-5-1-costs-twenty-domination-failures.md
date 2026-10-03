# 0036: The silence in §5.1 costs twenty domination failures

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** records a measurement and proposes a `NUMERICS.md` §5.1 addition; no code change.

## Context

The envelope's first run over the real corpus reports **75 domination failures** for the seeded
stand-in. [`0032`](./0032-domination-charges-helicoid-for-d16.md) (draft) accounted for 8 of them —
`so3_log`, where the whole gap is `libm::atan2` against the host's glibc, at most 0.74 u. The other
67 were unexplained, and the first useful step is not an experiment but a partition.

**Which oracle beat us is the discriminator.** `tf_tree_math` routes its transcendentals through the
`libm` *crate*, as we do, and `0032` measured it bit-identical to `helicoid` on all 30 `so3_log`
strata. sophus-rs uses Rust `std`, so glibc on Linux. So a stratum where `tf_tree_math` beats us
cannot be explained by D16's cost, and one where only sophus-rs beats us can be:

| class | count | reading |
|---|---|---|
| beaten by `tf_tree_math` (same `libm`) | **30** — `sen3_exp_n1` 29, `so3_exp` 1 | algorithmic; ours |
| beaten only by `sophus_rs` (glibc) | **45** — `sen3_jl_n1` 15, `sen3_jr_n1` 15, `so3_log` 8, `sen3_exp_n1` 7 | `0032`'s class; proven for `so3_log`'s 8 only |

One of the 30 is an artifact and is reported as a consequence below: `so3_exp` reads
$4.952463054990454 \times 10^{-1}$ u against `tf_tree_math`'s $4.952463054734329 \times 10^{-1}$ u —
the same number to ten significant figures, failing by $5 \times 10^{-11}$ relative.

That leaves **29 `sen3_exp_n1` strata where an oracle on our own `libm` is up to 2.30× more
accurate**, which is an algorithm question and the subject of this record.

## Measurement

`tf_tree_math`'s `exp_se3` never forms a matrix:

```rust
// V·v = v + c1·(ω × v) + c2·(ω × (ω × v)).
let wxv = w.cross(v);  let wxwxv = w.cross(wxv);
let t = v.add(wxv.scale(c1)).add(wxwxv.scale(c2));
```

The stand-in forms `J_l(φ) = I + aW + bW²` as a 3×3 matrix with `W²` as `W·W`, then applies it to
each `ρ_i`. The two are the same formula: `W r = φ × r` and `W² r = φ × (φ × r)`. They are not the
same program, and `NUMERICS.md` §5.1 does not say which to run — the stand-in's own module doc
records the choice as a reading where the spec is silent, and
[`0014`](./0014-the-conformance-metrics-open-readings.md) (draft) question 25 is the open question.

Three readings, one variable each, same corpus, same `libm`, same coefficients, binary64, 52
`sen3_exp_n1` strata. Domination failures for `sen3_exp_n1`, out of 75 total:

| reading | `sen3_exp_n1` failures | total |
|---|---|---|
| `I + aW + bW·W` as a matrix (today) | 36 | 75 |
| `r + a(φ×r) + b(φ×(φ×r))`, two cross products | **16** | 55 |
| `r + a(φ×r) + b(φ(φ·r) − θ²r)`, the `W² = φφᵀ − θ²I` identity | 19 | 58 |

Every other id is untouched by all three — `sen3_jl_n1` 15, `sen3_jr_n1` 15, `so3_log` 8, `so3_exp`
1, in every run — which is what makes this one variable.

**No reading dominates, and the structure is by θ regime.** On the worst max_u, which is `0006`'s
bar:

| regime | wins (of 3) | worst max_u: `W·W` / 2 cross / identity |
|---|---|---|
| θ away from π (35 strata) | 2 cross **23**, `W·W` 11, identity 1 | 2.1517 / **1.5087** / 1.5221 |
| θ near π (17 strata) | `W·W` **11**, identity 4, 2 cross 2 | 3.4678 / 4.6856 / **3.2855** |

The largest single improvements are at small and moderate θ — `rho:1e-6/theta=1e-1` goes 2.1517 u →
0.6000 u, a factor of **3.59** — and the largest regressions are at θ near π, where
`rho:1e-6/theta=pi-1e-6` goes 3.3739 u → 4.6856 u. Of 52 strata, the two-cross form improves 27,
worsens 17 and ties 8.

**All of this is measured on the stand-in, not on the shipped kernel.** The stand-in evaluates each
coefficient at its own generated switch where §4 groups them inside one `branch`, and its module doc
says so: "every error measured here is of this kernel, not of the grouped one". `PHASE3.md` §10 owes
`helicoid`'s own in-process subject for these ids, and until it exists no number here is a verdict on
the library. What transfers is the *question* and its price.

## Decision

Proposed, not settled.

1. **`NUMERICS.md` §5.1 should state how `J_l(φ) ρ` is applied**, because its silence is worth up to
   3.59× in ULP and 20 of 75 domination failures. A formula's accuracy is not a property of the
   formula; it is a property of the program, and a spec that fixes one and not the other has not
   fixed the result. This is a `NUMERICS.md` edit **and** this record (`CLAUDE.md`, decision
   workflow), and the edit is not made here.
2. **The evidence does not pick one form, so the decision is between two honest options**, and this
   record does not choose:
   - *Accept the trade*: take the two-cross form, best away from π and on 25 of 52 strata overall,
     and accept being worse near π, where the error is already 2–4 u and dominated by the
     coefficients.
   - *Choose by regime*: the two-cross form away from π, the identity form near it. That is a new
     switch, and [`0004`](./0004-switch-points-are-generated-not-typed.md) forbids typing one — it
     would have to come from `cargo xtask thresholds`, which today sweeps coefficient kernels and
     not assembly forms. That is a real extension of the sweep, not a line of code.
3. **Answering `0014` question 25 requires the grouped kernel.** The ranking above is the stand-in's.
   It should be re-measured once `PHASE3.md` §10's in-process subject exists, because the grouped
   `branch` changes which arm runs and therefore what the assembly is fed.

## Rationale

The partition came first for a reason: of 75 failures, 45 are in a class
[`0032`](./0032-domination-charges-helicoid-for-d16.md) already explains in kind, and chasing an
algorithm through them would have been chasing `libm`. Splitting by *which oracle won* costs one pass
over three CSVs and turns 67 unexplained failures into 30 worth an experiment and 45 worth a
confirmation.

Why the forms differ is visible in the arithmetic rather than asserted: `W·W`'s diagonal entries are
$\varphi_i^2 - \theta^2$, a difference of like-sized quantities, and the matrix path then spends
two three-term sums per output entry. The two-cross path spends two roundings per entry and forms no
matrix — which is why it wins where `a` and `b` are both significant, and loses near π where `a → 0`
leaves `b·(φ×(φ×r))` carrying the result and the second cross product's cancellation is exposed with
nothing to average it against.

## Consequences

- **Domination has no tolerance, and one of the 75 is an artifact of that.** `so3_exp` fails by
  $5 \times 10^{-11}$ relative — ten matching significant figures. `0006` says the bar is on the max
  and never a mean; it does not say the comparison is strict. A bar that reports a tie as a failure
  spends a reader's attention on nothing, and 1 of 75 here is nothing. Whether domination gains a
  tolerance (and what measures it) is a `0006` question this record raises and does not answer.
- The 45 in `0032`'s class are a confirmation, not a discovery: `sen3_jl_n1` and `sen3_jr_n1` are 30
  of them and have never had the one-variable test `so3_log` got. Until they do, "libm-bound" is a
  hypothesis for them and a measurement only for `so3_log`.
- If §5.1 fixes an application form, `helicoid`'s own `sen3_exp` must implement it and its reference
  twin (D6) must keep the other, since the twin's job is to disagree only where the fast path is
  wrong.

## Implementation plan

1. The partition above, as a reproducible analysis over `conformance/results/*.csv` — the three-CSV
   comparison is four lines of Python today and should be a `cargo xtask envelope` column so a
   reader sees *which* oracle won without writing a script. **Owed.**
2. The one-variable test for `sen3_jl_n1` and `sen3_jr_n1`, the shape `0032` used: swap only the
   transcendental and see whether the gap closes. **Owed.**
3. `NUMERICS.md` §5.1's application form, when this record is `ready`, with `0014` question 25
   closed by it. **Owed, blocked on the open questions.**
4. Re-measure the ranking against `PHASE3.md` §10's in-process subject once it exists. **Owed.**

## Open questions

1. **Accept the trade, or sweep the regime?** Option 2a takes the two-cross form and is worse near π;
   option 2b needs the threshold sweep extended from coefficient kernels to assembly forms, which is
   new machinery. The answer decides whether `0004` is involved at all.
2. Does the ranking survive the grouped kernel? The stand-in's readings are documented and differ
   from the shipped kernel's; `PHASE3.md` §10 is the only thing that can answer it.
3. Should domination carry a tolerance, and if so what measures it? Raised by the `so3_exp` artifact;
   `0006`'s to answer, not this record's.
4. Is there a fourth form? The three tested are the two `0014` q25 names plus the cross-product
   reading `tf_tree_math` happens to use. A Rodrigues-style assembly that shares `sin θ`/`cos θ` with
   the rotation was not tried, and `0014` question 19 (the `cos θ/2` reading) touches it.
