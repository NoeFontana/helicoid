# 0029: a `Product` needs a way to be built from its factors

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** —

## Context

`Product<A, B>` landed with both fields private, which
[`0025`](./0025-a-structured-jacobian-and-a-sealed-side.md) decision 5 requires and
[`PHASE3.md`](../PHASE3.md) §7 spells. It landed with no constructor either, and those two facts
together leave the type with no way in from its parts. Outside the crate the whole surface is
`LieGroup`: `exp`, `log`, `inverse`, `Mul`, `Blend`, the `*plus`/`*minus` family and the Jacobians.

`PHASE3.md` §7 names `Product<SO3<S>, Rn<S, 3>>` the tf2-semantics pose. A consumer holding a
rotation and a translation — which is what every pose source hands you — can therefore only write

```rust
let pose = Product::<SO3<S>, Rn<S, 3>>::exp(&(rot.log(), trans.log()));
```

and that composition is wrong in three separate ways for the job:

1. **It is not exact.** `Rn`'s `exp`/`log` are the identity, but SO(3)'s are not: the round trip
   spends an `atan2` (D5: `Log` goes through the quaternion `atan2`) and a `sin_cos`, and returns a
   rotation that differs from the input in the last bits. A pose that came in exact does not come out
   exact, and `tf_tree` composes poses along a chain, so the error accumulates over tree depth for no
   reason at all.
2. **It is worst-conditioned exactly where poses are ordinary.** `NUMERICS.md` §12 and
   `docs/maths/error-analysis.md` record `Log`'s conditioning growing as the rotation approaches a
   half turn. A half-turn pose is not a degenerate input to a transform tree — it is a camera looking
   backwards — but it is the worst input to this round trip.
3. **It costs transcendentals to copy six numbers.** The operation the consumer wants is a move.

The same need is already specified for the concrete group: `PHASE3.md` §5 plans `SE3` accessors
`rotation()`, `translation()`, `from_rt` and `from_quat_translation`. `Product` is the generic form of
exactly that shape and has none of them.

`PHASE3.md` §0.0's §7 row records this as the maintainer's to reconcile, which is where the question
has been parked; this record is that question written out. Reading is in the same state: there is no
way to recover the factors either, so a consumer that stores a pose as a `Product` cannot get the
rotation back out except through `log`.

## Decision

**Nothing is decided. A draft authorises nothing** — no constructor is added by this record.

A new public item is governed by [`API.md`](../API.md) §6, so the options are stated as surface:

**A. A named pair of accessors, fields staying private.**

```rust
impl<A, B> Product<A, B> {
    pub fn from_parts(a: A, b: B) -> Self;
    pub fn parts(&self) -> (&A, &B);
}
```

Additive, keeps `0025` decision 5 and `PHASE3.md` §7 both intact, and names the operation so that
`from_parts` reads as the move it is. It is also what `SE3::from_rt` will be, one level down, so the
two surfaces agree. R1 (every name says its frame or its side) is satisfied: there is no side to a
construction from parts, and `from_parts` is in the order the type parameters are.

**B. Public fields.** `Product<A, B>(pub A, pub B)`. One keyword each, and `exp`/`log` keep working.
It contradicts `0025` decision 5, which lists `Product` among the types that expose no field, and it
makes the tuple positions public surface rather than the names. It is also the irreversible
direction, which is why the implementation chose `pub(crate)` in the first place.

**C. Nothing until Phase 4 asks.** [`PHASE4.md`](../PHASE4.md) owns what `helicoid` owes the
`tf_tree` migration, and that migration is the first real consumer of a pose type. Deferring costs
nothing today because no consumer exists; it costs a release if `Product` ships `0.1` without it and
the idiom `exp(log(..))` has to stay supported afterwards.

The recommendation is **A**, and that it land with the `SO3` PR rather than before it: `SO3` is what
makes `Product<SO3, Rn<3>>` constructible at all, so it is the PR where the gap becomes real and
where the accessors can be tested against something that is not abelian.

## Rationale

**A over B** because `0025` decision 5 is `ready` and says `Product` exposes no field; a constructor
adds the capability without touching the rule, so no document has to change. Named accessors also
survive the type gaining a third factor through nesting, where tuple positions get worse:
`Product<Product<A, B>, C>::from_parts(ab, c)` still reads, while `.0.1` does not.

**A over C** because the cost of A is one additive item and the cost of C is a supported idiom. But
the *timing* argument for C is sound, which is why the recommendation adopts it: the record says what
the surface is so that the `SO3` PR implements rather than invents, and nothing lands until there is
a factor worth constructing.

`parts()` returning references rather than values is for a factor that is not `Copy`; every group
today is `Copy`, so this costs nothing and does not have to be revisited when one is not.

## Consequences

- If **A** lands, `Product` gains two public items and `API.md` §6's checklist applies to them:
  `from_parts` and `parts` are in §2's `helicoid` surface, and `PHASE3.md` §7 and its §0.0 row lose
  the "no way to build one" clause.
- The `exp`/`log` idiom keeps working either way, so nothing a consumer writes today breaks.
- A `Gaussian` over a `Product` (Phase 5) needs the same way in, so settling it here settles it
  there.
- Until it is settled, `PHASE3.md` §0.0's §7 row keeps carrying the gap, and the `tf_tree` migration
  of Phase 4 is where it stops being hypothetical.

## Implementation plan

1. The owner picks A, B or C and this record becomes `ready` — verified by
   `grep -m1 -H '^\*\*Status:' docs/decisions/0029-*.md`.
2. **A only**, with the `SO3` PR: `from_parts` and `parts`, the `PHASE3.md` §7 and §0.0 edits, and
   the `API.md` §2 surface entry — verified by a round-trip test that `from_parts(a, b)` then
   `parts()` returns both factors bit for bit, and by a test that `from_parts` on an `SO3` and an
   `Rn<3>` agrees with `exp(&(a.log(), b.log()))` to within the recorded `exp_log` bound while being
   bit-exact where the round trip is not.
3. **A only:** the generic laws gain an instance built through `from_parts` rather than `exp`, so the
   constructor is on the same footing as the rest of the surface — verified by `laws_for!` passing on
   it at `f64`, `f32` and `Dual`.

## Open questions

- **Does `parts()` want to be `into_parts()` as well?** Every group is `Copy` today, so a by-value
  form is free to add later and pointless now; if a factor ever holds a buffer, the answer changes.
- **Should the tuple struct become a named struct?** `Product { a, b }` would make `from_parts`
  redundant and the field names the surface. It is a larger change and `PHASE3.md` §7 spells the
  tuple form, so it is listed only so that the choice is visible.
- **Does `SE3` then need `Product` at all?** If `SE3` is its own type with `from_rt`, the tf2 pose may
  never be spelled as a `Product` in a consumer, and the gap matters only for products
  `helicoid` does not name. `PHASE4.md`'s migration answers this, and it is the strongest argument
  for option C's timing.
