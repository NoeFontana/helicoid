# 0054: the screw twin is the oracle's translation under the shipped rotation

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record and the twin land together

## Context

`PHASE4.md` §1.2 owes `SE3::geodesic`'s dual-quaternion fast twin, a port of
`tf_tree_math::dualquat::screw_pow`, and [`0041`](./0041-the-integration-is-an-adapter-at-tf-tree-math.md)'s
Wave 3 cannot delegate `ScLerp` without it: `SE3::geodesic` is the provided body, `X₀ Exp(t Log Δ)`.
Three things were open.

**Accuracy.** `se3_geodesic/geo:generic` was the project's last geodesic domination failure, 2.721
against `ScLerp`'s 2.502 `u`, and `0048` had shown it to be the *rotation*.

**What the port needs.** `docs/maths/index.md` held four items against §1.2: the `0/0` guard at
`‖v‖ = 0` is `f64`'s literal `1e-290`, which `f32` cannot represent and no document chooses; the
transcendental-free small-angle arm `screw_pow` ships needs two series outside `NUMERICS.md` §4,
each a catalogue entry and a record; whether `geodesic` is exact at `t ∈ {0, 1}`; and the tolerance
of `se3_geodesic_matches_reference`.

**The latency premise, which was wrong.** `0041` (Context, constraint 2) put Wave 3's bar at "≈241
ns": `lookup/depth3/sclerp`'s committed 192.7 ns plus `PHASE1.md` §9's 25 %. Read in `tf_tree`, no
check applies that rule to `ScLerp` — 192.7 ns is prose in `tf_tree`'s `0013`, measured on a host
that record calls unfit; `cargo xtask bench-gate` prints latency `UNAVAILABLE`; `LATENCY_SLACK`
times `LerpSlerp`; and the rule is `tf_tree`'s `PHASE1.md` §11.3, not §9. What binds is
`NS_PER_STEP_ESTIMATE = 64` (`tf_tree_py`), published through `tf_tree`'s `API.md` §3.4: its
compile-time GIL-crossover assertion holds for 56..=66 only, so a depth-3 median above about
200 ns is a published-API change. **The budget is parity, not 25 %.** And the lookup fixture's
edges publish at 50 Hz – 1 kHz with angular rates near 1 rad/s, so its relative rotations are near
`1e-3` rad: the bench measures `screw_pow`'s series arm alone.

## Decision

1. **`SE3::geodesic` is GE.12's power, whose rotation is `SO3::geodesic`'s to the bit.**
   `so3::geodesic_parts` is `SO3::geodesic`'s body, unchanged in every bit it returns, and also
   returns the relative rotation `q_rᵗ = (cos tα, ϖ_t v)`, `ϖ_t` and the flipped `q₀* q₁`; the SE(3)
   twin forms each translation column from GE.12's dual part over them. `NUMERICS.md` §10 states the
   formulas.
2. **The small-angle arm is the catalogue's.** Below `r`'s second switch, `ϖ_t = k·t·r` and
   `cos tα = c` with `r` from `log_ratio` and `(k, c)` from `exp_coeffs` — both on their short
   arms, polynomials, no transcendental (GE.15(a)). No new coefficient, series or switch.
3. **Above it, `ϖ_t = sin(tα)/‖v‖`** (GE.12's own), not the blend's `sin(tα)/sin α`, with
   `cos tα` from the same `sin_cos`.
4. **Both products are Hamilton's** — `½(0, x_Δ) ⊗ (w, v)` and `q_dᵗ ⊗ (q_rᵗ)*` through `Quat`'s
   `Mul` — `screw_pow`'s association, not the dot/cross grouping of the same terms.
5. **The guard is `ϰ = 0` where `‖v‖² ≤ 0`**, with `0003`'s safe argument in the division; no range
   constant (GE.15(b)). This answers the index item on `1e-290`.
6. **`N = 1` only.** `SEn3::geodesic` dispatches on the const `N`; `N ≥ 2` keeps the provided body
   until an `se23_geodesic` stratum can verify the power there (`0006`).
7. **Endpoints.** `t = 0` is exact, as before. At `t = 1` the rotation is exact (the blend) and the
   translation is `LieGroup::geodesic`'s "to rounding" (GE.13(f)); no `select(t == 1, x₁, …)`, which
   would zero `∂/∂t` under `Dual`. A consumer whose contract is exact endpoints keeps its own early
   return — `tf_tree`'s `ScLerp` does.
8. **`se3_geodesic_matches_reference`**: `10⁵` pairs per precision over three regimes (consecutive,
   generic, near `π`), `‖x₀‖` to `1e4`, the translation read against `max(‖x₀‖, ‖x₁‖, 1)`. No band
   near `π` is excluded — both twins read the sign of one `q₀* q₁` — and `t ∉ [0, 1]` has its own
   bound, since both twins' errors grow with `|t|`.

## Rationale

**The measurement** (`measure_geodesic::se3`, extended from the SO(3) instrument as `0050` *Further
work* 2 asked; `cargo test -p xtask -- --ignored --nocapture measure_the_se3_geodesic_spellings`),
`se3_geodesic`, binary64, max `u` per stratum:

| stratum | A provided | B GE.12 (`screw_pow`) | C B's translation, SO(3)'s rotation | **shipped** | `ScLerp` |
|---|---|---|---|---|---|
| `geo:consecutive` | 1.572 | 1.572 | 1.572 | **1.572** | 2.336 |
| `geo:generic` | 2.721 | 2.502 | 1.738 | **1.738** | 2.502 |
| `geo:near-pi` | 2.429 | 3.253 | 3.112 | **3.112** | 3.253 |

B reproduces the oracle's rows to the record (#79, #139), so it is a faithful model of `screw_pow`'s
arithmetic. **B → C is the rotation alone**, and it is the whole of `geo:generic`; the shipped
twin is C. `just envelope`: **113 → 112**, the one removed being `se3_geodesic/geo:generic`; no
geodesic row fails. `so3_geodesic` is unchanged to the bit (1.5721 / 1.7382 / 1.6417), the
refactor's control.

**Why not the blend's weight for `ϖ_t` (decision 3)**, which read 2.433 / 1.765 on the two upper
strata — a better `geo:near-pi` and a worse `geo:generic` than shipped. Each stratum's maximum there
is one or two records (#64, #68 against #133), which is the six-record resolution problem at 60. So
the two were compared where the noise is not: 20 000 draws per regime against GE.12 at 60 digits from
the exact binary64 inputs, translation error relative to `max(1, ‖x₀‖, ‖x‖)`, the `f64` arithmetic
in the shipped order:

| regime | `sin α` (blend's): max / mean / wins | **`‖v‖`**: max / mean / wins |
|---|---|---|
| generic, `‖x₀‖ ~ 1` | 17.23 / 1.560 / 1743 | **15.79 / 1.281 / 6021** |
| generic, `‖x₀‖ ~ 100` | 1.39 / 0.344 / 67 | **1.00 / 0.343 / 212** |
| near `π` | 21.79 / 1.867 / 1958 | **16.66 / 1.396 / 7814** |
| `θ ∈ [1e-2, 0.3]` | 25.60 / 1.398 / 1634 | **21.60 / 1.245 / 4624** |

`‖v‖` wins every regime on every statistic. The reason is structural: `ϖ_t v` then has norm
`sin tα` whatever `‖q₀* q₁‖` is, as the short side's `k·t·r` does; `sin α` recomputed from
`atan2` is `‖v‖/‖q‖`, so the relative power carries the input's unit-norm defect into the translation.
The blend keeps `sin α` because only there is it one number over itself at `t = 1` (`0050`); the
translation has no such endpoint to protect.

**Why Hamilton's association (decision 4).** On the same 20 000 draws the dot/cross grouping and the
Hamilton product tie — maxima identical, means 1.281 / 1.283, 1.396 / 1.399, 1.245 / 1.251 — so the
association is noise and the tie goes to the spelling the oracle uses, which also reads the better
corpus maximum (`geo:generic` 1.738 against the grouping's 2.433, the same #64).

**Why `ϰ = 0` (decision 5)** rather than GE.13(b)'s second arm, the reference twin behind a range
constant: the terms `ϰ` multiplies are `h t(1−t²) α²/6 · n̂ + O(α⁴)`, with zero gradient at
`q_v = 0` (GE.15(b), checked at 110 digits). Dropping them is exact at every precision, costs one
`select`, and leaves no constant for `f32` to need. `the_screw_twin_s_guard_keeps_value_and_derivative`
measures it at `‖v‖ ∈ {0, 1e-150, 1e-160, 1e-170}` (`f64`, `Dual`) and `{0, 1e-15, 1e-20, 1e-25}`
(`f32`), and `the_screw_twin_s_guard_differentiates_like_the_reference` the derivative lanes.

**Latency.** `bench-gate --bench groups --only 'se3/geodesic' --against <provided body, built in
this tree>`, core 5 pinned, concurrent A/A floor at most 0.97 %:

| `se3/geodesic` | `f64` | `f32` |
|---|---|---|
| `consecutive-1e-3` (new row: `tf_tree`'s regime) | 0.863 | 0.848 |
| `near-identity-7.5e-8` | 0.862 | 0.847 |
| `generic-1` | 0.721 | 0.729 |
| `near-pi` | 0.705 | 0.750 |

So `PHASE4.md` §1.2's condition — override only if faster than the reference — holds on every row.

**Against `screw_pow`, in the consumer.** A same-host criterion comparison first read the twin at
62.5 ns against `ScLerp` at 44.4 ns, which would have failed Wave 3. It was codegen, not
arithmetic: adding unrelated bench functions to the same binary moved the unchanged `ScLerp` to
64.6 ns and the provided body from 79.5 to 59.5, and behind `#[inline(never)]` wrappers the order
reversed (twin 60.0, `ScLerp` 66.5 ns; 81.2 against 82.0 at `generic`). A routine of a few hundred
instructions is decided by whether LLVM inlines it, so the question was taken to where it is asked:
a scratch copy of `tf_tree` at `20bc5a0` with `ScLerp::eval`'s body delegating to this twin (its
early returns kept), `lookup` built both ways, A B A B A B, core 5 pinned:

| `tf_tree` `lookup` | `screw_pow` | this twin |
|---|---|---|
| `depth3/sclerp` | 275.35 / 272.15 / 273.94 ns | 272.57 / 273.22 / 273.27 ns |
| `depth1/sclerp` | 95.02 / 95.87 / 95.02 | 94.44 / 95.81 / 94.33 |
| `depth6/sclerp` | 188.66 / 187.94 / 188.72 | 187.85 / 186.80 / 187.64 |

**Parity**, inside the run-to-run spread (this host's absolute is not `tf_tree`'s 192.7 ns; the
ratio is what transfers). A second build without two speculative changes — `#[inline]` on the
coefficient kernel, and `core`'s `copysign`/`abs` for `libm`'s, which disassembly showed staying out
of line — read the same (272.94 / 275.39 / 273.47 against 272.57 / 274.16 / 272.11), so neither is
taken.

## Consequences

- `laws::geodesic`'s `twin` leg now compares two expressions at `N = 1`: 8.051 `u` at binary64,
  8.587 at binary32 (SO(3)'s reads 7.213 for the same reason). `sen3_tests`'s shared bounds are
  re-recorded from `10⁶` draws: symmetry 33 → 30, right 19 → 17, `twin` 3 → 18.
- `se3_geodesic_matches_reference` reads, over `10⁶` draws per regime, `[13.40, 14.77, 12.81]` `u`
  inside `[0, 1]` and `[28.77, 31.46, 30.31]` outside at `f64` (`f32`: `[7.61, 11.66, 13.24]`,
  `[27.80, 30.62, 36.01]`); the bounds are twice that. These compare two programs; the corpus is
  the arbiter (GE.13(e)).
- Under `Dual` the twin keeps GE.13(c)'s loss: the translation's derivative errs by about
  `10² α⁻¹ u` at small angle, inherent to the `ϰ` grouping. `geodesic_jacobians` (`0043`) is a
  closed form and is not affected; a `Dual` through the twin owes that tolerance.
- **Wave 3's bar is parity.** `0041` is amended to say so; its Wave 3 is measured as above, in
  `tf_tree`'s own bench, both bodies, interleaved — and `NS_PER_STEP_ESTIMATE` is not re-derived
  unless that A/B moves the median.
- `se3_geodesic`'s `geo:near-pi` margin is 3.112 against 3.253, 4.3 %, on one record (#133).

## Implementation plan

1. `so3::geodesic_parts`, `SEn3::geodesic`'s override, `measure_geodesic::se3`, the bench row, the
   tests, GE.15, the `NUMERICS.md` §10 and §14 rows and `PHASE4.md` §0.0 — verified by `just lint`,
   `just test`, `just doc`, `just msrv`, `just no-std`, `measure_geodesic::tests::
   the_shipped_se3_route_reproduces_its_rows_and_dominates`, and `just envelope` at 112 with no
   geodesic row failing.
2. `tf_tree`'s Wave 3, in its own record, after the publish and Waves 1–2 (`0049`'s order): the
   adapter bodies of `ScLerp`, `screw_pow`, `screw_twist`, `screw_pow_with_twist`, the tests run
   under both feature states, and the A/B above repeated on `tf_tree`'s host — verified by
   `lookup/depth3/sclerp` at parity and `NS_PER_STEP_ESTIMATE` unchanged.

## Open questions

None.

## Further work

1. **A regrouped coefficient** `ψ_t = (ϖ_t cos α − t cos tα)/sin² α` (→ `−t(1−t²)/3`) would write
   the `ϰ` terms as `q_d,w ψ_t v` with no `O(‖x‖/α)` part, which GE.13(c) says no single-coefficient
   series removes from the current grouping. It would retire the guard and the `Dual` loss; it is a
   new two-variable coefficient, so a §4 edit and a record.
2. **An `se23_geodesic` corpus id**, then decision 6's dispatch removed: the power is written for
   every `N` already.
3. **An `@f32` geodesic stratum** (`0051` *Further work* 2): the twin is measured at `f32` by its
   proptests only.
4. **`libm::copysign` and `libm::fabs` are calls.** Not needed here, but every `Real::copysign` and
   `abs` in the workspace pays one; `core`'s methods are the same bits inline. A measurement on the
   groups bench decides it.
