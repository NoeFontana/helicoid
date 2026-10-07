# 0054: the screw twin takes two arms, in the world frame

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record and the twin land together (#105)

> **Further work 4 and 5 taken by [`0055`](./0055-the-branch-is-inlined-and-the-sign-is-a-bit.md)** (2026-10-07): `Real::branch` is `#[inline]` and
> `abs`/`copysign` are `core`'s; coefficient kernels 0.55×, `groups` 0.915×, bit-identical.

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

1. **The rotation is `SO3::geodesic`'s, to the bit.** `so3::{geodesic_rel, geodesic_short,
   geodesic_long}` are that routine's prefix and two arms; `SO3::geodesic` is one `S::branch` over
   them, and so is the twin, which forms rotation and translation in the same arm.
2. **The translation is computed in the world frame**, on `δ = x₁ − x₀` and
   `v_w = vec(q₁ q₀*)` (Log's flip applied): the power conjugated by `R₀` (GE.15(a)), so no column
   is rotated in or out.
3. **Below `r`'s second switch** (exactly `t = 1` included), the translation is the definition,
   `t J_l(tφ) J_l⁻¹(φ) δ` with `φ = r v_w`, by two cross products each way as `log` and `exp` write
   them, its coefficients §4's `r`, `c`, `a`, `b` on their series arms. No transcendental, no new
   coefficient, and no cancelling one.
4. **Above it**, GE.12's dual part, `ϖ_t = sin(tα)/‖v‖` (not the blend's `sin(tα)/sin α`) and
   `cos tα` from one `sin_cos`, its products Hamilton's and `½(0, δ) ⊗ q` written without its zero
   lanes. The `select`s on `‖v‖ = 0` are `0003`'s safe arguments for a lane that is not selected;
   no range constant.
5. **`N = 1` only.** `SEn3::geodesic` dispatches on the const `N`; `N ≥ 2` keeps the provided body
   until an `se23_geodesic` stratum can verify the twin there (`0006`).
6. **Endpoints.** `t = 0` is exact, as before. At `t = 1` the rotation is exact (the blend) and the
   translation is `LieGroup::geodesic`'s "to rounding" (GE.13(f)); no `select(t == 1, x₁, …)`, which
   would zero `∂/∂t` under `Dual`. A consumer whose contract is exact endpoints keeps its own early
   return — `tf_tree`'s `ScLerp` does.
7. **Tests.** `se3_geodesic_matches_reference`: `10⁵` pairs per precision over three regimes
   (consecutive, generic, near `π`), `‖x₀‖` to `1e4`, the translation read against
   `max(‖x₀‖, ‖x₁‖, 1)`, no band near `π` (both twins read the sign of one `q₀ · q₁`), `t ∉ [0, 1]`
   its own bound. `the_screw_twin_differentiates_like_the_reference`: `Dual` lanes against the
   reference's, per angle.

## Rationale

**The first version was GE.12 at every angle, and its review found the defect that shapes this
one.** That version read the `ϰ` grouping below the switch too, with `ϖ_t = k·t·r` and
`cos tα = c` from `log_ratio` and `exp_coeffs`, and `ϰ = 0` at `‖v‖² = 0`. Its values were right.
Its derivatives were not: `ϰ` multiplies `ϖ_t w − t cos tα`, an `O(α²)` difference (GE.15(b)) —
a cancelling coefficient computed inline, which `CLAUDE.md` reserves for `helicoid::coeffs` — and
through `Dual` its `O(‖x‖/α)` terms cancel in the assembly (GE.13(c)). Twin against reference,
lanes on the relative tangent, worst over `t ∈ {¼, ⅓, 0.7314, 1}`:

| `θ` | `1e-12` | `1e-9` | `1e-6` | `1e-3` | `1e-2` | `0.1` | `1` |
|---|---|---|---|---|---|---|---|
| GE.12 at every angle | 1.7e12 | 2.0e9 | 5.4e6 | 2426 | 314 | 34 | 3.7 |
| **two arms** | **2.3** | **2.3** | **3.4** | **2.3** | **2.0** | 33.5 | 4.7 |

`ϰ = 0` was the limit only at exactly `‖v‖² = 0`; a normal `n²` of `1e-300` kept the loss. Below
the switch the definition has none — its coefficients are the catalogue's — and above it
`α ≥ α_s` bounds GE.12's.

**The measurement** (`measure_geodesic::se3`, extended from the SO(3) instrument as `0050` *Further
work* 2 asked; `cargo test -p xtask -- --ignored --nocapture measure_the_se3_geodesic_spellings`),
`se3_geodesic`, binary64, max `u` per stratum:

| stratum | A provided | B GE.12 (`screw_pow`) | C B's translation, SO(3)'s rotation | **shipped** | `ScLerp` |
|---|---|---|---|---|---|
| `geo:consecutive` | 1.572 | 1.572 | 1.572 | **1.572** | 2.336 |
| `geo:generic` | 2.721 | 2.502 | 1.738 | **2.057** | 2.502 |
| `geo:near-pi` | 2.429 | 3.253 | 3.112 | **2.721** | 3.253 |

B reproduces the oracle's rows to the record (#79, #139), so it is a faithful model of `screw_pow`'s
arithmetic, and B → C is the rotation alone, which is the whole of `geo:generic`'s loss. `just
envelope`: **113 → 112**, the one removed being `se3_geodesic/geo:generic`; no geodesic row fails.
`so3_geodesic` is unchanged to the bit (1.5721 / 1.7382 / 1.6417), the refactor's control.

The corpus cannot rank C against the shipped spelling: each stratum's maximum is one or two of its
60 records (#73, #75, #125, #133), the six-record problem again. So the choices were measured
where the noise is not: `2×10⁴` draws per regime against GE.12 at 60 digits from the exact binary64
inputs, the translation relative to `max(1, ‖x₀‖, ‖x‖)`, `f64` arithmetic in the shipped order.

**World frame (decision 2).** Against GE.12 in `x₀`'s frame, Hamilton products both:

| regime | body frame: max / mean / wins | **world frame** |
|---|---|---|
| generic, `‖x₀‖ ~ 1` | 15.79 / 1.283 / 3074 | **6.51 / 0.723 / 12018** |
| generic, `‖x₀‖ ~ 100` | 1.00 / 0.343 / 124 | **0.95 / 0.343 / 407** |
| near `π` | 16.66 / 1.399 / 3958 | **8.99 / 0.870 / 11762** |
| `θ ∈ [1e-2, 0.3]` | 21.60 / 1.251 / 2624 | **6.89 / 0.645 / 12348** |

Two rotations per column are two roundings of `δ`'s size; the world frame has none, and it is
cheaper.

**`ϖ_t = sin(tα)/‖v‖` (decision 4).** In the body frame, against the blend's `sin α`: generic
1.281 / 1.560 mean, 15.79 / 17.23 max, 6021 / 1743 wins; near `π` 1.396 / 1.867, 16.66 / 21.79,
7814 / 1958. `ϖ_t v` then has norm `sin tα` whatever `‖q₀* q₁‖` is; `sin α` recomputed from `atan2`
is `‖v‖/‖q‖`, which carries the input's unit-norm defect into the translation. The blend keeps
`sin α` because only there is it one number over itself at `t = 1` (`0050`).

**Hamilton's association (decision 4).** The dot/cross grouping of the same terms tied it on the
same draws (means 1.281 / 1.283, 1.396 / 1.399, maxima equal), so the tie went to `screw_pow`'s.

**Latency.** Each geodesic row now times the routine through an `#[inline(never)]` call
(`groups::geodesic_of`), because inlined into criterion's closure a routine this size is timed
with LLVM's choice there: an unchanged `screw_pow` read 44.4 and 64.6 ns in two builds differing
only in other bench functions. `bench-gate --bench groups --only geodesic --against <dispatch
reverted, built in this tree>`, core 5 pinned, concurrent A/A floor at most 2.98 %:

| `se3/geodesic` | `f64` | `f32` |
|---|---|---|
| `consecutive-1e-3` (new row: `tf_tree`'s regime) | 0.921 | 0.898 |
| `near-identity-7.5e-8` | 0.920 | 0.899 |
| `generic-1` | 0.615 | 0.594 |
| `near-pi` | 0.598 | 0.610 |

Every `so3` and `se23` geodesic row reads 1.00 ± 0.01. One branch on the hot path is load-bearing:
the same arithmetic with the rotation and the translation in **two** `S::branch` calls read
**1.18×** on the short rows, the bench binary keeping each `Real::branch` — which has no
`#[inline]` — out of line with its captures.

**Against `screw_pow`, in the consumer.** A scratch copy of `tf_tree` at `20bc5a0` with
`ScLerp::eval`'s body delegating to this twin (its early returns kept), `lookup` built both ways,
A B A B A B, core 5 pinned:

| `tf_tree` `lookup` (ns) | `screw_pow` | this twin | |
|---|---|---|---|
| `depth1/sclerp` | 94.41 / 94.45 / 95.17 | 90.41 / 90.01 / 90.35 | −4.6 % |
| `depth3/sclerp` | 272.66 / 272.13 / 275.23 | 256.77 / 255.28 / 255.45 | −6.1 % |
| `depth6/sclerp` | 186.60 / 186.26 / 188.43 | 176.83 / 176.59 / 176.46 | −5.3 % |

This host's absolute is not `tf_tree`'s 192.7 ns; the ratio is what transfers. Two speculative
changes — `#[inline]` on the coefficient kernel, `core`'s `copysign`/`abs` for `libm`'s, which
disassembly shows out of line — moved neither this A/B nor a pinned microbench, and are not taken.

## Consequences

- `laws::geodesic`'s `twin` leg now compares two expressions at `N = 1`: 8.051 `u` at binary64,
  8.579 at binary32 (SO(3)'s reads 7.213 for the same reason). `sen3_tests`'s shared bounds are
  re-recorded: symmetry 33 → 30, right 19 → 17, `twin` 3 → 18; `N = 1`'s symmetry fell 16.111 →
  11.415.
- `se3_geodesic_matches_reference` reads `[23.5, 24.1, 28.2]` `u` inside `[0, 1]` and
  `[54.9, 54.2, 54.3]` outside at `f64` (`f32`: `[35.9, 21.3, 24.2]`, `[56.9, 62.5, 86.5]`); the
  bounds are twice that. The difference is mostly the **reference's**, which rotates the
  translation out and back: at `t = 1`, where the truth is `x₁`, the draws differing by more than
  `12 u` read 2.6 `u` mean (11.1 max) for the twin and 12.6 (35.9) for the reference.
- Above the switch a `Dual` through the twin keeps GE.13(c)'s loss, bounded by the switch:
  33.5 `u` at `θ = 0.1`. `geodesic_jacobians` (`0043`) is a closed form and is not affected.
- **Wave 3's bar is parity**, and the twin clears it with room. `0041` is amended to say so; its
  Wave 3 is measured as above, both bodies, interleaved, on `tf_tree`'s host, and
  `NS_PER_STEP_ESTIMATE` is not re-derived unless that A/B moves the median.
- `se3_geodesic`'s tightest margin is `geo:generic`, 2.057 against 2.502.

## Implementation plan

1. `so3::{geodesic_rel, geodesic_short, geodesic_long}`, `SEn3::geodesic`'s override,
   `measure_geodesic::se3`, the `consecutive-1e-3` row and `geodesic_of`, the tests, GE.15, the
   `NUMERICS.md` §10 and §14 rows and `PHASE4.md` §0.0 — verified by `just lint`, `just test`,
   `just doc`, `just msrv`, `just no-std`, `just wasm`,
   `measure_geodesic::tests::the_shipped_se3_route_reproduces_its_rows_and_dominates`, and
   `just envelope` at 112 with no geodesic row failing.
2. `tf_tree`'s Wave 3, in its own record, after the publish and Waves 1–2 (`0049`'s order): the
   adapter bodies of `ScLerp`, `screw_pow`, `screw_twist`, `screw_pow_with_twist`, the tests run
   under both feature states, and the A/B above repeated on `tf_tree`'s host — verified by
   `lookup/depth3/sclerp` at parity or better and `NS_PER_STEP_ESTIMATE` unchanged.

## Open questions

None.

## Further work

1. **A regrouped coefficient** `ψ_t = (ϖ_t cos α − t cos tα)/sin² α` (→ `−t(1−t²)/3`) would let
   GE.12 run at every angle with no cancellation, retiring the short arm's three coefficients for
   one; it is a new two-variable coefficient, so a §4 edit and a record, and the short arm is
   already faster than `screw_pow` where it matters.
2. **An `se23_geodesic` corpus id**, then decision 5's dispatch removed: the twin is written for
   every `N` already.
3. **An `@f32` geodesic stratum** (`0051` *Further work* 2): the twin is measured at `f32` by its
   proptests only.
4. **`Real::branch` has no `#[inline]`.** Here it cost 18 % until the routine was restructured
   around it; every other two-`branch` routine in the workspace may be paying the same. The groups
   bench, re-gated with it inlined, decides it.
5. **`libm::copysign` and `libm::fabs` are calls**, where `core`'s methods are the same bits
   inline. Not needed here; the groups bench decides it.
