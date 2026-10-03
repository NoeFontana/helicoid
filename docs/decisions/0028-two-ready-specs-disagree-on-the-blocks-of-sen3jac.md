# 0028: two ready specs disagree on the blocks of `SEn3Jac`

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

`SEn3Jac<S, N>` ships with both of its fields public:

```rust
pub struct SEn3Jac<S, const N: usize> {
    pub diag: Mat3<S>,
    pub col: [Mat3<S>; N],
}
```

Two `ready` documents say different things about that, and both are authoritative:

- [`PHASE3.md`](../PHASE3.md) §5 spells the type with `pub diag` and `pub col`. The implementation
  followed the spec, which is the rule.
- [`0025`](./0025-a-structured-jacobian-and-a-sealed-side.md) decision 5 states that a field which is
  a **representation choice** rather than the value stays private regardless, and names `RnJac`'s `k`
  as the case: `Rn`'s `Vector` is public because it *is* the value, `RnJac`'s `k` is private because
  `k I` is a storage decision (D2, values not storage).

`diag` and `col` are the second kind, and the type's own rustdoc says so in as many words — it
describes the dense matrix and then says it is "stored as `diag = A` and `col = [B₁, …, B_N]`". The
sibling `ProductJac`, landed under `0025`, took the other reading: its rustdoc states that the blocks
"are a representation choice and stay private, as `RnJac`'s `k` does". So the workspace currently
applies both rules, one per type, decided by which document the implementer read.

Three things follow from the fields being public, none of them intended by either document:

1. **An unlabelled Jacobian sum is reachable.** `Mat3` implements `Add` and `Sub`, so
   `SEn3Jac { diag: a.diag + b.diag, col: … }` type-checks. The type's rustdoc states the opposite
   — "There is no `PartialEq`, `Add` or `Sub`" — and `0002` forbids an `Add` impl for `⊕` precisely
   so that no composition happens without naming its side. The fields re-open by hand what the
   missing impl closes.
2. **A dense block escapes without `write_dense`.** `j.diag` hands out a `Mat3` directly, where
   [`API.md`](../API.md) R5 and `0005` have a structured Jacobian produce a dense form only through
   `write_dense` or `sandwich`.
3. **It is a permanent commitment.** D2 says layout is not a semver contract, and these fields *are*
   the layout. Widening later is additive; narrowing is breaking. The workspace is `0.0.x`, so the
   cost of narrowing is near zero today and rises with every consumer.

Nothing is measured here and nothing about the numerics changes either way: this is a surface
question only.

## Decision

**Option A: `SEn3Jac`'s block fields narrow to `pub(crate)`, and `PHASE3.md` §5's declaration is
corrected to match.** Decided under `0040` item 5, on reversibility: A is the only option that can
be undone without a breaking change, no external consumer reads `j.diag` today, and `0005`'s
successor Jacobians — `SO3Jac`, `Sim3Jac`, the Phase 6 ambient ones — are all unwritten, so the rule
is settled before they copy a precedent rather than after. `0025` decision 5 and its rustdoc's "no
`Add`" guarantee both stand unchanged, which B would have had to break.

The three options as they were weighed, with what each costs:

**A. Narrow to `pub(crate)` and correct `PHASE3.md` §5.** The reversible direction, and the one
`ProductJac` already took. `Gaussian` and the ambient Jacobians are inside `helicoid`, so
`pub(crate)` reaches every planned in-crate consumer. Cost: `PHASE3.md` §5's declaration is edited,
which is a `ready` spec changing because a later `ready` record overrode it — exactly what `0025`
did to §7's `Jac = Mat<N>`, so there is precedent for the mechanism. Any external consumer that
reads `j.diag` today breaks, and none exists.

**B. Keep `pub` and correct `0025` decision 5.** Decision 5 would have to say that a block of a dual
matrix is the value and not a representation, which is hard to write without also making `RnJac::k`
public, since `k I` is no more a representation than `A + εB` is. It would also leave the rustdoc's
"no `Add`" sentence false, so that sentence would have to go, and with it the guarantee.

**C. Keep `pub`, change neither document, and accept that the rule is per-type.** The cheapest now
and the most expensive later: the next structured Jacobian has no rule to follow, and the question
returns with `SO3Jac`, `Sim3Jac` and the Phase 6 ambient Jacobians.

**A is taken**, on reversibility alone: it is the only option that can be undone without a breaking
change, and `0005`'s successor Jacobians are all still unwritten, so the rule should be settled
before they copy a precedent rather than after.

## Rationale

The reason this is a record and not a PR is that both documents are `ready`. `CLAUDE.md` is explicit
that each spec's §0.0 status table is the source of truth and that an agent stops and asks when the
documents do not answer a question; here they answer it twice, differently. Narrowing a public field
on the strength of one of two conflicting specs would be choosing which spec wins by fiat, which is
what the decision folder exists to prevent.

`0021` is the precedent for a record whose content is "measured, and nothing changes"; this is the
surface analogue: argued, and nothing changes until the owner rules.

## Consequences

- Until this is settled, `SEn3Jac` and `ProductJac` demonstrate opposite conventions, and
  `PHASE3.md` §0.0's §5 row carries the conflict so that the next implementer meets it.
- If **A** is taken, the edit is two visibility keywords, one `PHASE3.md` §5 declaration, and the
  deletion of nothing else: no call site inside the workspace reads either field outside
  `dualmat.rs` and its tests, and `reference.rs`'s twins reach the blocks through `Jac`.
- If **B** is taken, the rustdoc's "no `PartialEq`, `Add` or `Sub`" sentence becomes false and must
  be replaced by what is actually guaranteed, and `0025` decision 5 needs a new boundary that keeps
  `RnJac::k` private for a reason that is not "it is stored".
- Either way, `SO3Jac`, `Sim3Jac` and the Phase 6 ambient Jacobians inherit a rule instead of a
  coin toss.

## Implementation plan

1. The owner picks A, B or C and this record becomes `ready` with that decision — verified by
   `grep -m1 -H '^\*\*Status:' docs/decisions/0028-*.md`.
2. **A only:** the two `pub` keywords, the `PHASE3.md` §5 declaration and the §0.0 row — verified by
   `just build`, `just test` and `just doc` on an unchanged test suite, since no test reads the
   fields through anything but `dualmat_tests.rs`, which is inside the crate.
3. **B only:** the `0025` decision-5 boundary, the `ProductJac` rustdoc sentence that cites it, and
   the `SEn3Jac` rustdoc's guarantee sentence — verified by `cargo xtask lint` and `just doc`.

## Open questions

- **Does any planned consumer need `A` alone?** `PHASE3.md` §0.0's §5 row already records that
  `SEn3Jac::inverse` wants `A⁻¹` from the SO(3) closed form, which `inverse(&self)` cannot receive;
  if the answer to that is a constructor or an argument rather than a field read, option A costs
  nothing at all.
- **Is the rule about `Add` or about layout?** If the objection to `pub` is only that `Mat3: Add`
  makes an unlabelled sum reachable, a newtype without `Add` would keep the fields public and close
  the hole; that is a third shape neither document considered.
- **Does `0025` decision 5 bind types it did not name?** It names `RnJac`, `Rn`, `SO3`, `SE2`,
  `SEn3`, `Sim3` and `Product`, and `SEn3Jac` is not in the list. Whether the decision is a rule or
  an enumeration decides whether `PHASE3.md` §5 is even in conflict.
