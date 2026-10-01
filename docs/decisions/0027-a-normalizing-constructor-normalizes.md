# 0027: A normalizing constructor normalizes

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #22

## Context

Two normative documents disagreed about one function, and the code shipped one of the two readings.

`API.md` R6 says, normatively: "Unit-norm inputs are checked by `*_unchecked` constructors in debug
only; **`*_normalized` constructors normalize**." `NUMERICS.md` §3.6 says: "Construction from
external data and the explicit `renormalize` apply the first-order Newton step
$q \leftarrow q\,(3 - \|q\|^2)/2$." `Quat::from_wxyz_normalized` implemented §3.6, and
`PHASE3.md` §0.0 booked the contradiction as "for the maintainer to reconcile".

The step is not a normalization. With $\eta = \|q\|^2 - 1$ it leaves
$\|q'\|^2 - 1 = -\tfrac34\eta^2 + \tfrac14\eta^3$ exactly (`maths/so3.md` SO.14), below rounding
only for $\lvert\eta\rvert \le 2^{-26.29}$ (`f64`) or $2^{-11.79}$ (`f32`) (SO.15). Outside that it
is a worse guess the further out it starts, it returns **zero** at $\eta = 2$, and it reverses $q$
beyond. Measured consequences of the shipped behaviour, all silent in debug and in release:

| input | $\|q\|^2$ | returned | `to_matrix` of it |
|---|---|---|---|
| `(0, 1, 1, 1)` — 180° about an unnormalized $(1,1,1)$ | 3 | `(0, 0, 0, 0)` | the zero matrix |
| `(3, 0, 0, 0)` | 9 | `(-9, 0, 0, 0)` | $81\,R(q/\|q\|)$, reversed |
| `(0, 0, 0, 0)` | 0 | `(0, 0, 0, 0)` | the zero matrix |

A zero quaternion sends every acted point to the origin. Nothing reports it, because §3.6 and §12
state no domain for the step, so there is nothing for a `debug_assert!` to enforce.

The regime matters to the resolution. A unit quaternion drifts by at most about $16u$ per Hamilton
product (SO.15), so from unit norm it takes about $2^{22}$ compositions in `f64` to leave the step's
accuracy domain; and $2^{-40} \ll 2^{-26.29}$, $2^{-16} \ll 2^{-11.79}$, so **every** quaternion
`from_wxyz_unchecked` admits is inside that domain with 13.7 bits of margin in `f64` and 4.2 in
`f32`. The step is a true normalization for drift repair. It is wrong only where the norm was never
near one — which is exactly what a constructor taking external data receives.

## Decision

1. **`Quat::from_wxyz_normalized` divides by the norm**: one `Real::sqrt` of `norm_sq` and one
   division per component. `NUMERICS.md` §3.6 is amended to say so, and to keep the Newton step for
   `renormalize` alone. `API.md` R6 stands unchanged and is now satisfied as written.

2. **Its domain is `‖q‖²` normal** — nonzero, with the components within about $10^{\pm154}$
   (`f64`) or $10^{\pm19}$ (`f32`) — stated in a `# Domain` section and enforced by a
   `debug_assert!` **on the result**: that it is unit to the tolerance `from_wxyz_unchecked`
   accepts. One assertion covers the zero quaternion (quotients are NaN), the overflow (quotients
   are zero) and a NaN input, and it states the contract that matters: this constructor produces
   what that one admits. Release builds check nothing (D11); at `q = 0` the result is NaN, which
   reaches whatever the caller computes, as `read_dense` poisons a short buffer
   ([`0025`](./0025-a-structured-jacobian-and-a-sealed-side.md)).

3. **`renormalize` keeps the Newton step, its name and its silence.** It asserts nothing, is the
   only `&mut self` numeric method (R2), and its rustdoc points at `from_wxyz_normalized` for a
   quaternion that is not already near unit norm.

4. **The rule is general, not a patch to one function.** `*_normalized` divides by the norm;
   `renormalize` is the Newton step. This binds `SO3::from_quat_normalized` and the SO(2)
   normalization of `PHASE3.md` §4 and §6 before they are written.

5. **The `from_matrix` question is not settled here.**
   [`0015`](./0015-specification-gaps-found-while-building-the-instrument.md) NU.8 stays open: its
   recommendation (a) — divide by the norm in `from_matrix` — is consistent with this record and
   now has one less thing to decide, but the sign rule it also asks about is untouched. NU.8's
   parenthetical "keep Newton for `renormalize` and `from_wxyz_*`" is **superseded** for
   `from_wxyz_*` by decision 1; `0015` is a draft and authorises nothing, but a reader will find
   that line.

## Rationale

Each of the three operations gets one job, and the names say which:

| | vouches for the norm? | cost | algorithm |
|---|---|---|---|
| `from_wxyz_unchecked` | the caller does (`debug_assert!`) | none | none |
| `from_wxyz_normalized` | the caller cannot | one `sqrt`, four divisions | divide by `‖q‖` |
| `renormalize` | it was unit and drifted | nine multiplies | Newton step |

Two alternatives lost.

**Rename `from_wxyz_normalized` to `from_wxyz_renormalized`** and keep the step. This satisfies R6
by vacating the name rather than by meeting the contract, and leaves the method a spelling of
`let mut q = Quat { .. }; q.renormalize(); q` — three lines of caller code — while the need it was
added for, *make this unit, I cannot vouch for where it came from*, stays unserved. Renaming also
costs the same edits across `NUMERICS.md`, `API.md` §3, `PHASE3.md` §4 and the forward
`SO3::from_quat_normalized`, and buys nothing at the call site.

**Amend R6 to carve out an exception.** R6 is one sentence covering every `*_normalized`
constructor in the workspace. Weakening it so that one function can keep a cheaper algorithm
inverts the cost: every future reader of every such constructor has to go and check which kind it
is.

The cost of deciding this way is one `sqrt` and four divisions on a constructor. Composition is the
hot path, not construction, and a caller who knows the norm is already close has
`from_wxyz_unchecked` (free) or `renormalize` (no `sqrt`). Dividing four times rather than
multiplying by one reciprocal is deliberate: the quotient rounds once instead of twice, which is
worth three cycles on a path whose whole purpose is to produce a trustworthy unit quaternion.

Determinism is unaffected. `Real::sqrt` is exactly rounded on every target
([`0018`](./0018-libm-arch-is-bit-identical-for-exact-operations.md)) and division is an IEEE
operation, so the result is bit-identical on x86_64, aarch64 and wasm32 (D16).

## Consequences

- **A behaviour change on `0.0.x`.** `from_wxyz_normalized` returns different values for every
  input that is not already unit — including the ones in the table above, where it previously
  returned zero or a reversed quaternion. No consumer calls it yet; `SO3` is not written.
- **A new recorded bound.** `norm`, the worst $\lvert\eta\rvert$ of `from_wxyz_normalized` over
  $10^6$ seeded non-unit quaternions: **6.00 u** (`f64`) and **5.50 u** (`f32`), recorded as 12 and
  11 by the doubling rule of `quat_tests.rs`. Both are far inside the `debug_assert!`'s tolerance
  of $2^{13}u$ and $2^{8}u$, which is why the assertion on the result is affordable as the domain
  check.
- **`measure_worst_errors` gains a second RNG stream** for the new law's scale, so that the eight
  previously recorded figures are drawn from exactly the samples they were recorded from. They
  reproduce unchanged.
- **The step's own guarantees are now stated where they bind.** `renormalize` is documented as drift
  repair, and the proof that its domain contains `from_wxyz_unchecked`'s is in the Context above,
  not re-derived at each call site.
- **Still owed, and not changed here:** the `Quat` corpus stratum and reference twin
  ([`0006`](./0006-the-instrument-comes-first.md), `PHASE3.md` §0.0). The nine laws are proptests
  against recorded bounds, which is not domination over an oracle envelope, and this record does not
  pretend otherwise.

## Implementation plan

1. `from_wxyz_normalized` divides by the norm, with the `# Domain` section and the
   `debug_assert!`; `renormalize`'s rustdoc cross-reference is updated — verified by
   `from_wxyz_normalized_is_unit_far_from_unit_norm` (the three table rows and two exact cases),
   `newton_step_hand_case` (the step and the constructor now differ on one input),
   `out_of_domain::normalized_rejects_the_zero_quaternion`,
   `out_of_domain::normalized_rejects_an_overflowing_norm_in_f32`,
   `out_of_domain_does_not_panic_in_release` (NaN, not zero), and the `from_wxyz_normalized_is_unit`
   law under `f64`, `f32` and `Dual<f64, 4>`. Reverting the body to the step fails five of these in
   debug and six in release. **Landed in #22.**
2. `NUMERICS.md` §3.6 and the §12 row; `PHASE3.md` §0.0 and §4; this record. **Landed in #22.**
3. Apply decision 4 to `SO3::from_quat_normalized` when SO(3) lands, and to the SO(2)
   normalization if `0015` NU.10 is adopted — verified by the same pair of laws per type.
4. The `Quat` corpus stratum, which is what finally gates this surface
   ([`0006`](./0006-the-instrument-comes-first.md)) — owed by `PHASE3.md` §0.0, not by this record.

## Open questions

None. The `from_matrix` sign rule and the SO(2) normalization are `0015` NU.8 and NU.10, where they
already are.
