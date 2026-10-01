# 0025: a structured Jacobian and a sealed `Side`

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #21

## Context

The trait layer of `PHASE3.md` §2 and the group of §7 landed as one PR: `Tangent`, `Jac`,
`LieGroup`, `Side`/`Right`/`Left`, `Rn` with `RnTangent` and `RnJac`, and the generic law harness
`laws.rs` with its `Bounds` and `laws_for!`. Four things in it are not what §2 and §7 said, and the
PR edited those sections in the same change. `decisions/README.md` reserves that order for a
record: a new public item starts as a draft record, and the rationale lives in the record, not in a
§0.0 status cell.

This record was written **after** the code, during the review of #21, and says so rather than
implying it came first: decisions 1 and 5 describe what that PR already contained, decisions 2, 3
and 4 are changes the review added to it, and all five land in #21. The alternative — leaving the
two unbuilt decisions to a follow-up — was rejected so that the knowledge and the code land
together; nothing is published (`0.0.x`), so narrowing a surface now costs nothing.

**1. §7's `Jac = Mat<N>` cannot implement `Jac`.** §7 specified the Jacobian of Rⁿ as
`Jac = Mat<N>`, identity-structured. `Jac::inverse` is required and is the exact algebraic inverse,
not a numerical one (`0005`), and `Matrix<S, N, N>` has no inverse at general `N`: `PHASE2.md` §4
gives `Mat3::inverse_adj` and the Cholesky path, nothing else. §7 named a type that cannot satisfy
the trait it assigned, and it did so by making the Jacobian dense, which is what `0005` exists to
prevent.

**2. `Side` is public, unsealed, and carries no side discriminator.** A group selects its side row
by comparing `TypeId`s; the test group's `is_right` is the pattern every later group would copy.
`PHASE3.md` §0.0 defers *which* selector to the SO(3) PR. It does not address who may implement the
trait. A downstream `impl Side for MyRight` whose `plus` delegates to `rplus` gets the correct
`plus` and `minus` and the **Left** row from `compose_jacobians` and `inverse_jacobian` — a
silently wrong Jacobian handed to a solver, with no compile error and no runtime signal. Nothing
ties an impl's `plus` and `minus` to the same side either.

**3. `Tangent::dot`'s stated order is not achievable by `Product`.** §2 specifies `dot` as "the
Euclidean inner product of the dense components, summed in index order", and `laws::
tangent_dense_order` pins it bit-exactly: the recorded `f64` bound is `0` for both groups, against
a reference that accumulates `((0 + d₀e₀) + d₁e₁) + …` over the whole dense vector. §7 also
schedules `Product<A, B>`, whose natural `dot` is `self.0.dot(&o.0) + self.1.dot(&o.1)` — two
partial sums, not one index-order accumulation. It differs in the last bits and cannot meet a bound
of `0`. The bound is the point of that law: it is a *specification-conformance* statement about
dense order, not an accuracy statement, and consumers index that order (a solver's residual
blocks, `tf_tree`'s arena records).

**4. `Tangent::read_dense` reads a missing entry as `+0`.** §2 fixes only the signature; the fill
value is this PR's own rustdoc. `+0` is a valid component, so a wrongly sized caller buffer is
indistinguishable from a legitimately zero tangent, and the caller gets a plausible wrong update
instead of a diagnostic. D11 forbids the release check that would say so.

**5. `Rn` and `RnTangent` expose fields whose type implements `Add` and `Sub`.** `Rn<S, N>(pub
Vector<S, N>)` and `RnTangent { pub rho }`, and `Vector<S, N>` implements `Add`/`Sub` in
`helicoid-linalg`, so `Rn(a.0 + b.0)` is `a * b` and `RnTangent { rho: t.rho - u.rho }` is
`t.sub(&u)`. `API.md` R1 bans those impls on the group and tangent types, which holds; the two
`compile_fail,E0369` doctests pin only the outer types. `RnJac` keeps `k` private.

## Decision

**1. Rⁿ's Jacobian is structured, and §7 is corrected, not `0005`.** `RnTangent { rho:
Vector<S, N> }` is the tangent and `RnJac`, a scalar multiple `k I` with `k` private, is the
Jacobian: closed under `mul`, `inverse` and `neg`, holding exactly the values Rⁿ needs — `I` for
`Ad`, `J_r`, `J_l` and their inverses, `-I` for the second row of `⊖` and of `X⁻¹`, and `0` for
`ad`. `RnJac::inverse` states `k != 0` as its domain, since `ad` is `0` and has no inverse.
`PHASE3.md` §7 is edited to say so; `0005` is unchanged, and this closes a §7 violation of it.

**2. `Side` is sealed. The selector stays with the SO(3) PR.** A private supertrait makes `Right`
and `Left` the only inhabitants:

```rust
mod sealed {
    pub trait Sealed {}
}

pub trait Side: sealed::Sealed + Copy + 'static { /* unchanged */ }

impl sealed::Sealed for Right {}
impl sealed::Sealed for Left {}
```

No `unsafe`, no `alloc`, no runtime cost. **Whether a group reads the side from `TypeId` or from a
first-class selector on `Side` remains the SO(3) PR's decision**, as `PHASE3.md` §0.0 states; this
record only removes the third-party impl.

**3. The flat accumulation is the primitive.** `Tangent` gains a threaded accumulator as the
required method and `dot` becomes provided:

```rust
/// The Euclidean inner product of the dense components, accumulated into `acc` in index order.
fn dot_acc(&self, o: &Self, acc: S) -> S;

/// `dot_acc(o, +0)`.
fn dot(&self, o: &Self) -> S {
    self.dot_acc(o, S::zero())
}
```

A leaf tangent accumulates its own components in index order. `Product<A, B>` threads:

```rust
fn dot_acc(&self, o: &Self, acc: S) -> S {
    self.1.dot_acc(&o.1, self.0.dot_acc(&o.0, acc))
}
```

which is the flat left-to-right sum over the concatenated dense order at any nesting depth, with no
buffer. Seeding with `+0` is normative and is what the law's reference does: `0 + (-0)` is `+0`, so
seeding from the first product instead would differ on a signed zero. `dot`'s stated contract in §2
is unchanged — it is now achievable.

**4. `read_dense` poisons a missing entry.** A release build reads a missing entry as NaN, not
`+0`. The debug `debug_assert!(src.len() == DOF)` and the no-release-panic guarantee (D11) stand.
`write_dense` keeps writing `min(out.len(), DOF)` entries: it writes into the caller's buffer, so
there is nothing to poison.

**5. A field whose type carries `Add`/`Sub` is exposed only where no side can be mislabeled.**
`Rn`'s field and `RnTangent::rho` stay `pub`, and `PHASE3.md` §7 spells `pub`. The rule the later
groups inherit:

- A **group** type exposes such a field only when the group is **abelian**, so that `⊕_R = ⊕_L` and
  the reach-through cannot express an unlabeled side. `SO3`, `SE2`, `SEn3`, `Sim3` and `Product`
  expose none; their components are reached through `exp`/`log` and `write_dense`/`read_dense` (R3).
- A **tangent** type is exempt: a tangent is a vector space, `Tangent::add` and `sub` carry no side,
  and a field-wise sum cannot be mistaken for `⊕`.
- A field that is a **representation choice** rather than the value stays private regardless, which
  is why `RnJac::k` is private while `Rn`'s `Vector` is public (D2, values not storage).

## Rationale

Decision 1 is forced: §7 as written is unimplementable, and the structured form is what `0005`
already requires of every other group. Nothing was traded.

Decision 2 is about **preserving the SO(3) PR's freedom, not pre-empting it.** Sealing is
backward-compatible to add now and a breaking change to add later. Sealed, the exhaustiveness that
`is_right` already assumes becomes true, so a later move from `TypeId` to a `const IS_RIGHT: bool`
or a `pick` method is a non-breaking internal change. Unsealed, that freedom is spent the moment
the crate is published, and the cost of leaving it is a wrong Jacobian with no diagnostic — the
failure class `0005` and D11 are both written to avoid. Defensive alternatives (an `is_right` that
panics on an unknown type) were rejected: a sealed trait makes the case unreachable instead of
reporting it, and D11 reserves panics for the strided accesses.

Decision 3 chose the threaded accumulator over the two obvious alternatives. Narrowing the contract
to "in index order within each factor" would make the law's reference mirror the tangent's factor
structure, which destroys the genericity that makes `tangent_dense_order` worth having — it is
currently one flat reference for every group. Relaxing the bound from `0` to a few `u` would leave
the law unable to distinguish a correct dense order from a different one, which is the only thing
it tests. Re-flattening both tangents into a `DOF`-sized stack buffer inside `Product::dot` also
works and needs no `alloc`, but it copies on every call and the threaded form does not.

Decision 4 follows the harness's own practice: `laws::dt` and `laws::dj` already read through
NaN-poisoned buffers so that an unwritten entry fails a law. The same reasoning applied to the
public boundary gives NaN, not `+0`.

Decision 5 is narrower than it first looks, which is why `pub` survives. R1 exists so that `⊕` and
`⊖` always name their side; on an abelian group the sides coincide, so there is no side to
mislabel, and `a.0 + b.0` on Rⁿ is visibly vector addition rather than a disguised group law. On a
non-abelian group the same expression is not merely unlabeled but wrong — `SO3(q_a + q_b)` is not a
rotation — so exposure is barred there for a stronger reason than R1. Making the trivial group's
field private, with no accessor, would buy nothing real: an accessor returning `Vector` would
itself be a second dense path and so violate R3 while fixing R1.

## Consequences

- `PHASE3.md` §7 and `API.md` §3 state `RnTangent`/`RnJac` and `pub` on `Rn`'s field; `API.md` R1
  gains the abelian condition as a stated exception rather than leaving it to a reader.
- The §0.0 rows for §2 and §7 shrink to `Partial` plus this record's number. Their present contents
  — the `TypeId`-versus-selector question and the "rows are not compared with `Dual`
  differentiation" gap — move here and to §2/§7 prose; a status cell is not a record.
- `Side` becomes unimplementable downstream. This is a deliberate narrowing of a surface that is
  not yet published (the workspace is on `0.0.x`), and `0012` already routes a consumer that wants
  different behaviour to a `Chart` rather than to its own impl.
- Every `Tangent` impl owes `dot_acc`; `dot` is no longer overridden anywhere. The bound-`0`
  `tangent_dense_order` becomes a permanent guarantee for composite groups, not a property Rⁿ and
  the test group happen to have.
- `Tangent::read_dense`'s release behaviour changes, so `rn_tests`'
  `out_of_domain_does_not_panic_in_release` changes with it.
- **Out of scope, owed a record of its own:** `Jac::write_dense` reaches a `StridedMut` only
  through `set`, which range-checks and recomputes a saturating index per entry — 81 bounds
  branches and 81 panic calls for 81 stores at `DOF = 9`, measured from the emitted asm in release.
  Fixing it needs a checked-once bulk writer on `StridedMut`, a new public item in
  `helicoid-linalg` governed by `PHASE2.md` §5, not by this record.

## Implementation plan

Steps 1 to 4 land together in #21: the record is retrospective, so the usual one-step-one-PR
breakdown collapses into the PR being reviewed. Step 5 waits for `Product`.

1. Seal `Side` with the private supertrait; `PHASE3.md` §2 states that `Right` and `Left` are the
   only impls and that the selector is still the SO(3) PR's call — verified by a
   `compile_fail,E0277` doctest that implements `Side` for a local type with all four methods, so
   the only error is the unsatisfied supertrait. A doctest compiles as a downstream crate, which is
   exactly the case at issue; removing the seal makes the doctest fail, so it is not vacuous.
2. `Tangent::dot_acc` required, `dot` provided; `RnTangent` implements `dot_acc` — verified by
   `tangent_dense_order` holding at a recorded `f64` bound of `0` for both groups.
3. `read_dense` fills a missing entry with NaN; `traits.rs` and `rn_tests.rs` state it — verified by
   `out_of_domain_does_not_panic_in_release`.
4. `PHASE3.md` §2, §7 and the two §0.0 rows; `API.md` R1 and §3 — verified by `cargo xtask lint`
   (line and draft citations) and `just doc`.
5. `Product<A, B>` implements `dot_acc` by threading — verified by `tangent_dense_order` at bound
   `0` for the product group. Lands with §7's `Product`, not before.

## Open questions

None.

Two were resolved before this became `ready`. **Sealing lands in #21**, not a follow-up: amending
the PR advances the base branch of #22 without rewriting history, so the stack above it does not
rebase, and capturing the knowledge with the code was worth more than the open surface. **`dot_acc`
is not a second spelling under `API.md` §6 Q5**: it is the required operation and `dot` is its
unit-seeded specialization, the required/provided pattern `LieGroup` uses throughout, and only
`dot` appears in §2's surface list. The decisive argument is not taste — the flattening-buffer
alternative cannot be written on stable Rust at all, since `Product`'s buffer would have to be
sized from `DOF`, an associated const.
