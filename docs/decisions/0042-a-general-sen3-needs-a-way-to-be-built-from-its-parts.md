# 0042: a general `SEn3` needs a way to be built from its parts

**Status:** draft (awaiting a decision)
**Owner:** @NoeFontana
**Implementation:** `crates/helicoid/src/sen3.rs` (`SEn3::from_parts`, `SEn3::parts`), landed with
the SE_N(3) PR under the reading below; this record is the question written out.

## Context

[`0029`](./0029-a-product-needs-a-way-to-be-built-from-its-factors.md) asked this of `Product` and
answered it: private fields plus no constructor leave a type reachable only through an `Exp`/`Log`
round trip, which is inexact, worst-conditioned at a half turn, and spends two transcendentals to
copy numbers. It took option A — `from_parts`/`parts` — and noted that the concrete group's named
accessors, `PHASE3.md` §5's `SE3::from_rt`, `translation()` and the `SE23` pair, are the same shape
one level down.

[`PHASE3.md`](../PHASE3.md) §5 specifies those accessors **at `N = 1` and `N = 2` only**. The type
is generic in `N`, and the corpus is not: `PHASE1.md` §4.3 holds `sen3_{exp,log,ad,jr,jl,jr_inv,
jl_inv}_n{1,2,3}`, 21 ids, and seven of them at `N = 3`. Three of those seven — `log`, `ad` and
their widths — take a record whose input is `(q, x₁, x₂, x₃)`: a group element, handed to the
subject as its parts. With the §5 surface alone there is **no way to build it**, so `helicoid`
cannot answer its own corpus at `N = 3`, and the strata stay unscored rather than unpaired.

`Exp` is not a way round it: the element a `sen3_log_n3` record names is not `Exp` of anything the
record holds, and manufacturing one would score a different input than the reference's.

## Decision

**`SEn3<S, N>` gains `from_parts(SO3<S>, [Vec3<S>; N])` and `parts(&self) -> (SO3<S>, [Vec3<S>;
N])`**, `0029`'s names at `0029`'s level, with the fields staying private.

- By value, not by reference, where `0029` returns `(&A, &B)`: `SO3` and `[Vec3<S>; N]` are `Copy`
  and at most thirteen scalars, and `rotation()` already returns by value.
- `from_parts` takes an `SO3`, not a `Quat`, so the unit invariant is carried by the type that
  states it. That is a *typing* argument and not a checking one: `SO3::from_quat_unchecked` is a
  move with no assertion of its own (`Quat::from_wxyz_unchecked` carries `NUMERICS.md` §3.6's
  `debug_assert!`), so a constructor that states the §3.6 domain must reach that item — which
  `SE3::from_quat_translation` now does, and `from_parts` does not state a domain at all.
- The named accessors stay and stay preferred: `from_rt` and `from_quat_translation` read as what
  they are at the width where a reader knows what `t` means.
- `API.md` §6 governs a new public item; the pair is additive, names the operation as the move it
  is, and adds no trait, no dependency and no `Pod`-adjacent surface (D2).

## Rationale

The alternative is to let the conformance subject reach the fields — a `pub(crate)` constructor
plus a crate-internal subject — and that is the one thing `0006` forbids by construction: the
subject that scores the shipped code must be an *external* consumer of it, or the measurement is of
something other than what ships (`PHASE1.md` §5, `0004` item 4). Making the harness a privileged
reader would retire that guarantee for every id at once, to save two methods.

Keeping `N = 3` unanswered was the other option. It costs the corpus's whole `N = 3` third — the
only width where `SEn3Jac`'s `col` array has more than two blocks, so the only one where a
per-block bug cannot hide behind symmetry.

## Consequences

- Two public items on a type whose fields stay private, and a second way to build an `SE3` beside
  `from_rt` (`from_parts` of a one-element array). The named one is documented as preferred; no
  lint enforces it.
- `parts` is the first public reader of all `N` columns. A consumer storing an `SE23` can now
  recover `v` and `p` without `Log`, which is `0029`'s argument unchanged.

## Open questions

1. Should `from_rt` become `from_parts`'s alias rather than its own body? They differ only in the
   array literal today.
2. `SE23` has no `from_parts`-shaped named constructor (`from_rvp`?), which `PHASE3.md` §5 does not
   list either. Left unasked until a consumer holds a velocity and a position.

## Further work

1. `PHASE3.md` §5's accessor bullet is written per `N`; it gains the generic pair when this record
   is `ready`.
