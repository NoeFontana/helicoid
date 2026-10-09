# 0063: the sphere reads its sign by comparison, and is held unit

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #118 (the record; its code is the S² PR, `PHASE5.md` §2)

## Context

`PHASE5.md` §2 and `NUMERICS.md` §8 specify S² and its frozen Householder chart. `docs/maths/charts.md`
derives and checks them (CH.8–CH.11), and in doing so finds four gaps that `docs/maths/index.md`
lists as open. Each one changes what an implementer would write.

1. **The sign of zero.** §8 sets `ν = n + sgn(n_z) e_z` with `sgn(0) = +1`, and §2 has a stratum
   at `n_z = ±0`. With `n_z = −0.0`, the reading `n_z ≥ 0` gives `ς = +1` and `copysign(1, n_z)`
   gives `−1`. Both bases are valid (CH.11(d)), and they differ by the reflection `T` of CH.11(a):
   `‖B₊ − B₋‖_F = 2`. So `s2_retract` and `s2_local` at that stratum differ by `O(1)` between the
   readings.
2. **The unit norm.** §12 has no row for `S2`, while the quaternion has two
   (`from_wxyz_unchecked`'s `2^-40`, the carried band of `0058`). Neither §8 nor §2 says whether
   `S2` renormalizes. With `‖n‖ = 1 + ε`, the `H` form of the basis keeps `BᵀB = I` but has
   `Bᵀn = O(ε)`, and the closed form CH.8(d) loses both by `2ε` (CH.8(e)).
3. **The near-antipode reference.** §2 defines `s2_local`'s reference as the minimal rotation "via
   `mp.expm`/`mp.logm`". `mp.logm` of a rotation is complex and non-principal from `α ≈ 3.03`:
   20 of 20 samples per angle at `3.03, 3.04, 3.05, 3.1, 3.14`, at 50 digits on mpmath 1.3.0 and
   1.4.1 (CH.9). So the `s2:near-antipode` stratum would be scored against a wrong answer.
4. **No metric for a direction.** §11 scores vectors by a relative 2-norm. A unit-vector output is
   a direction. Its error is the angle to the reference, which for small errors the 2-norm gives
   to first order but which is exact for any error. §11 also does not score drift off the sphere.

§8 also has a factor-of-two slip: "the ratio `α/‖n × m‖` is §4's `r`". By CH.9(c),
`α/s = ½ r(s², w) = r(s², 1 + w)`, with `r(n², w) = 2 atan2(n, w)/n` (`0015` (draft) NU.6).

### Measured

**The renormalizing step on 3-vectors.** `n ← n (3 − ‖n‖²)/2` is §3.6's step, here in three
components. numpy 2.5.3, 20 000 unit draws per row, `η = ‖n‖² − 1` set to `±2^ℓ × U(0.5, 1)`,
the step in the stated precision, and `|‖n′‖² − 1|` evaluated in extended precision. The worst
values were:

| binary64, `ℓ` | −40 | −30 | −27 | −26.29 | −26 | −24 |
|---|---|---|---|---|---|---|
| worst, `u` | 4 | 4 | 5 | 5 | 6 | 27 |

| binary32, `ℓ` | −16 | −13 | −12 | −11.79 | −11.5 | −10 |
|---|---|---|---|---|---|---|
| worst, `u` | 4.44 | 4.34 | 4.91 | 4.72 | 5.46 | 15.3 |

The quaternion's band, `2^-26.29` and `2^-11.79` (SO.15), is S²'s band too: inside it, one step
lands within `5u` of the sphere.

## Decision

1. **`ς = +1` unless `n_z < 0`.** In code that is `S::select(n_z.lt(S::zero()), −1, 1)`, so `+0`,
   `−0` and NaN all give `+1`. §8's "`sgn(0) = +1`" becomes this sentence, and CH.11(d)'s two
   valid choices become one. This does not touch §3.2, where the quaternion's `w` at `±0` is read
   by its sign bit (`0055`): there the bit picks between two logarithms of the same rotation
   (§3.2). Here it picks a frame, and a `−0` produced by arithmetic says nothing about
   geometry, so the chart should not depend on it.
2. **`S2` is held unit, with three entry points**, as the quaternion's (`0027`, `0058`). `S2<S>`
   wraps a `Vec3<S>` in a private field, read with `S2::vec`.
   - `S2::from_vec_unchecked(v)` vouches. It `debug_assert!`s `|‖v‖² − 1| ≤ 2^-40` (`f64`) or
     `2^-16` (`f32`), and a NaN fails.
   - `S2::from_vec_normalized(v)` divides by the norm (`0027`), with domain `‖v‖²` normal.
   - `S2::renormalize(&self)` is the one Newton step, a normalization for
     `|‖n‖² − 1| ≤ 2^-26.29` (`f64`) or `2^-11.79` (`f32`). It asserts nothing and is defined
     everywhere, as `renormalize` is for the quaternion.

   `S2Chart::at` and `local` take `S2` values, so their inputs are inside the vouched bound. §12
   gains three rows, one per entry point, copied in form from the quaternion's.
3. **§8's ratio.** `α/s = ½ r(s², w)`, taken from `log_ratio` in the first form, with the mask on
   `s²/w²` as CO.16(c) and (e) require. The second form, `r(s², 1 + w)`, is named and rejected:
   it is 3 to 5 times less accurate near the antipode (CH.9(c), CH.11(e)).
4. **`s2_local`'s reference** is the geometric logarithm, never `mp.logm`. At 60 digits it computes
   `E_nm = Exp(α m̂)` by Rodrigues, then `Log(E_nm)` through the quaternion's `atan2` (`0045`'s
   reference for the geodesics), projected by `Bᵀ`, with `B` the basis of decision 1 at the input
   `n`. This agrees with CH.9(b)'s closed form to `2.6e-104` (CH.9, *Checked*). `PHASE5.md` §2's
   sentence changes to say so. The `s2:nz0` stratum holds records at `n_z = +0` and at
   `n_z = −0`, and both references use `ς = +1`.
5. **§11's direction metric.** For an output `â` whose reference `a` is a unit vector
   (`s2_retract`, and any later bearing), the score is
   `max(atan2(‖â × a‖, â · a)/u, |‖â‖² − 1|/(2u))`: the angle between the two, and the drift off the
   sphere, which the angle cannot see. `atan2` keeps the angle accurate at both ends, where `acos`
   of the dot product loses half the digits near `0`. `s2_local`'s output is a tangent in `ℝ²` and
   keeps the 2-norm forward error.
6. **`docs/maths/index.md`.** The three S² open items (the sign of zero, the `mp.logm` reference,
   the unit-norm row) are marked answered by `0063`.

## Rationale

- **A comparison, not the bit.** The alternative, `copysign`, makes the frame of a chart depend
  on how a zero was produced. `−(+0)`, or a product with a negative factor that underflows, would
  move it by a reflection. Under decision 1 the chart at `n` is a function of the point.
  The quaternion's bit reading stays because there it does not choose a frame.
- **Three entry points, not a tolerance alone.** A tolerance tells the caller what is legal but
  gives a drifted `n` no way back. The quaternion's three-entry arrangement is already the house
  pattern, and its band was measured above to hold for three components.
- **The geometric reference.** It is the same answer the mathematics defines, without
  `mp.logm`'s branch failure. `0045` made the same change for the geodesics.
- **Angle plus drift.** The 2-norm agrees with the angle to first order, so this changes no small
  score. It is exact at large errors, and the drift term catches a retract that leaves the sphere
  while pointing the right way.

## Consequences

- `NUMERICS.md` §8, §11 and §12, and `PHASE5.md` §2, change as above. `API.md` §2 lists `S2`'s
  constructors and `renormalize`.
- The `S2` and `S2Chart` implementation (`PHASE5.md` §2) builds on this record. It is not in this
  record's plan.
- The direction metric is a new rule in `xtask/src/conformance/metric.rs`, owed when `s2_retract`
  lands.

## Implementation plan

1. This record and the edits of *Consequences*. Verified by `just lint` and `just doc`.
2. The S² PR (`PHASE5.md` §2) implements decisions 1, 2 and 3. It adds the generator's
   `s2_retract` and `s2_local` with decision 4's reference and the `s2:nz0` records at `±0`, and
   the decision 5 metric rule. Its tests include `basis_ignores_the_sign_of_zero` (bit-identical
   `B` at `±0`), `renormalize_lands_within_five_u_in_the_band` at both precisions, and the
   `from_vec_unchecked` assert. Verified by `just test`, `just corpus-check`, `just conformance`
   and `just envelope`.

## Open questions

None.

## Further work

- locus-tag's bearings (its camera models' unprojection) are directions too. If it adopts `S2`, its
  gates can use decision 5's metric.
