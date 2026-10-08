# 0058: a drifted quaternion is carried, not vouched for

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** the PR that lands this record

## Context

`NUMERICS.md` §3.6 says composition never normalizes, so a chain of products drifts. Repairing
the drift is the caller's explicit `renormalize`. §12 gives the domain of exactly one entry
point, `Quat::from_wxyz_unchecked`: `|η| ≤ 2^-40` (`f64`) or `2^-16` (`f32`), where
`η = ‖q‖² − 1`. A `debug_assert!` enforces it, and a NaN fails it. A quaternion reaches a group in
two ways today, and §12 describes only the first:

- **Vouched for.** The caller claims the quaternion is unit and the assert checks the claim:
  `Quat::from_wxyz_unchecked`, `from_wxyz_normalized`, and `SE3::from_quat_translation`, which
  goes through the first.
- **Carried.** A struct literal (`Quat`'s fields are public) or `Quat::from_xyzw`, then
  `SO3::from_quat_unchecked`, a move that checks nothing. `SE3::from_quat_translation`'s rustdoc
  already says so. `Exp`, `Mul` and `from_matrix` build their outputs this way.

`SO3::from_quat_unchecked`'s rustdoc still says the caller guarantees §3.6's bound. So a
consumer whose chains drift by design, such as `tf_tree`'s `Iso3`, has no stated contract once
`|η|` passes `2^-40`, although §3.6 says the chain will get there. `tf_tree`'s first adapter
wave (its draft record `0063`) ran into this. Its adapter built through `from_wxyz_unchecked`,
one test composing 50 000 times tripped the assert, and `Quat::rotate` was kept on the native
path because of it.

### Measured

All figures are in `u = 2^-53`, at `f64`. The protocols are the tests named here, all in
`so3_tests::carried`.

**Budget.** `a_repeated_step_drifts_linearly_and_one_step_repairs_it`: `10^4` products of one
repeated step drift by `0.608`, `0.730` and `0.000 u` per product for its three steps. Fresh
random steps drift at most `0.024 u` per product over `10^4` to `10^6` products (a scratch run
of the same loop). A repeated step therefore reaches `2^-40` after about `1.1 × 10^4` products,
and the edge of the band where one Newton step still normalizes, `2^-26.29`, after about
`1.6 × 10^8`. One `renormalize` from `2^-40` lands within `4 u`. `tf_tree`'s depth-limited
lookups (`MAX_DEPTH = 32`) measured a worst `|η|` of `4e-15` over 10 000 chains, about 230 times
under `2^-40`.

**Each operation, at a carried `λq`, `λ² = 1 + η`.** `each_operation_errs_by_its_first_order_term_across_the_drift_band`
and `the_geodesic_tilts_by_the_norm_difference_only_on_the_blend` take 4 000 draws at each of
`|η| ≤ 2^-40, 2^-36, 2^-32, 2^-28, 2^-26.29`. The worst error against each operation's
first-order term:

| Operation | First-order term at `λq` | Worst residual |
|---|---|---|
| `log` | the `Log` of `q/‖q‖`: §3.2's scale invariance | `3.959 u` relative |
| `act` | `Rv + η (Rv − v)`, so `|error| ≤ 2|η| ‖v‖` | `12.961 u ‖v‖` |
| `to_matrix` | `(1 + η) R`, §1's scaled rotation | `10.885 u ‖v‖` |
| `Mul` | `η_ab = η_a + η_b + η_a η_b`; norms multiply | `6.000 u` |
| `inverse` | `η` unchanged | bit-exact |
| `geodesic`, provided body (below `r`'s second switch) | the geodesic of `q₀/‖q₀‖` and `q₁/‖q₁‖`: scale-invariant | `5.612 u`, at every `η` |
| `geodesic`, the blend | tilted toward the longer endpoint by `sin((1−t)α) sin(tα)/sin α · |η₁ − η₀|` in rotation angle, with `α` the half-angle of `q₀*q₁` | at most `1.0027×` the term |
| `geodesic` at `t = 0` and `t = 1` | the carried endpoint | bit-exact |

On a carried quaternion, every operation's error is its first-order term plus the rounding
error it has on unit input. Only `act`, `to_matrix` and the blend have a nonzero term. The
blend's tilt follows from perturbing its weights. Scaling `q₁` by `1 + δ`, with `δ = (η₁ − η₀)/2`
to first order, moves the blended arc point by `sin((1−t)α) sin(tα)/sin α · δ` in quaternion
angle, which is twice that in rotation angle. The measured ratio stays at 1 to the third digit.

**The consumer it replaces.** A scratch run used the same draws through `tf_tree_math`'s native
`slerp` and `Quat::rotate`:
- `rotate` is bit-identical to `SO3::act` on every one of 100 000 drifted draws.
- `slerp`'s blend has the same tilt, within 0.2%.
- On near pairs (relative angle below `1.7e-4`), `slerp` loses `5.7 u` at `|η| ≤ 2^-40` and
  `3 984 u` at `2^-26.29`, because its chord-based angle reads the norm difference as rotation.
  `geodesic`'s provided body stays at `6.5 u` (those draws, in `u = 2^-53`).

So what a carrying consumer gives up by delegating is nothing. It gains the near-pair arm.

## Decision

1. **Two entry points, named.** The *vouched* path is `Quat::from_wxyz_unchecked` and the
   constructors that reach it. It is unchanged: `2^-40` or `2^-16`, a `debug_assert!`, and a NaN
   fails it. That bound is the claim the caller makes. The *carried* path is a struct literal or
   `Quat::from_xyzw`, moved into `SO3::from_quat_unchecked` (and `SE3::from_rt`). It checks
   nothing, in debug or release, and propagates NaN. Its domain is **the drift band**,
   `|η| ≤ 2^-26.29` (`f64`) or `2^-11.79` (`f32`), which is `renormalize`'s accuracy domain. One
   step brings any carried quaternion in the band back within the vouched bound. `SO3::from_quat_unchecked`'s
   rustdoc changes to say this. No code changes and no item is added.
2. **§12 states each operation's error on a carried quaternion**, as the table above does. These
   are first-order terms, not new formulas. Every operation computes exactly what it computed
   before. SE_N(3) inherits the rows: the rotation is SO(3)'s, and each translation column of a
   product passes through `act`, so it carries `2|η_a| ‖x_b‖`.
3. **§3.6 states the budget.** A repeated step drifts up to `0.73 u` per product. A consumer that
   composes more than about `10^4` times without a `renormalize` has left the vouched bound, and
   one that composes more than about `10^8` times has left the band. Where in the chain to
   renormalize is the consumer's decision (`0027`, `0049` decision 3).
4. **Unchanged:** the `2^-40` and `2^-16` bounds, the NaN-failing assert, `renormalize`, every
   formula, every corpus stratum, and the open question on which reading `act` takes at
   `q:nonunit` (`docs/maths/index.md`). At `|η| = 2^-45`, the measured `510.3 u` there is
   decision 2's `2|η|/u`.

## Rationale

Four alternatives were weighed.

- **Let NaN pass the vouching assert.** This was the first recommendation `tf_tree`'s integration
  got, and it is wrong. A NaN is the loudest false claim a caller can make, and that assert exists
  to catch false claims. The functions that must return NaN to their callers, `Exp` and
  `from_matrix`, already avoid it, each with its reason in its rustdoc. A consumer that is total
  in NaN takes the carried path, which never had an assert.
- **Widen `2^-40` to the drift band.** This would make every vouching caller's claim 13.7 bits
  weaker, to suit a consumer that is not making the claim. The bound is what the vouched path's
  accuracy rests on, and decision 2 states what the carried path costs without touching it.
- **Add an assertion-free constructor.** It already exists twice: the struct literal and
  `from_xyzw` followed by the move. A third spelling would break `API.md` R3.
- **Have the consumer renormalize in its hot path, or in the adapter.** This costs every call a
  Newton step that its lookups do not need: `4e-15` against `1.2e-8`. It changes bits between
  the two sides of the consumer's feature, hides a drift the consumer should see, and does
  nothing about NaN. The measured budget says where renormalizing pays, which is in the consumer's
  long accumulators. That is `0027`'s and `0049`'s position.

A debug assert at the band's edge on the carried path was also considered and rejected. It would
panic on NaN, which is exactly what a total consumer needs the carried path not to do.
`renormalize` asserts nothing for the same reason (§12).

## Consequences

- A consumer whose chains drift by design builds with the struct literal and
  `SO3::from_quat_unchecked`. Inside the band, it is inside a stated domain with a stated error.
  `tf_tree`'s adapter does this, and `Quat::rotate` delegates. Both are in `tf_tree`'s own
  record, not here.
- `log_is_scale_invariant_across_the_unchecked_domain`'s comment is corrected. Outside the vouched
  bound, a quaternion can be built in a debug build: by the carried path.
- Four tests in `so3_tests::carried` pin decision 2's rows, the NaN propagation and decision 3's
  budget, with bounds recorded as `Bounds` asks.
- No release is owed. The change is rustdoc, `NUMERICS.md` and tests.

## Implementation plan

1. This record, `NUMERICS.md` §3.6 and §12, `SO3::from_quat_unchecked`'s and
   `Quat::from_wxyz_unchecked`'s rustdoc, and the `so3_tests::carried` tests, verified by
   `just test`, `just lint` and `just doc`.

## Open questions

None.

## Further work

- The rows are measured at `f64`. The first-order terms do not depend on precision. An `f32`
  sweep would confirm the residuals and is owed before an `f32` consumer relies on them.
- The SE_N(3) inheritance in decision 2 is derived from `act`'s row and not separately sampled.
