# 0048: The relative-transform pair earns the surface; `dot` and `norm` do not

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** landed. `Quat::{dot, norm}` removed; `SEn3::mul_inv` at every `N` with its
mirror `SEn3::inv_mul`; `SEn3`'s `rminus`/`lminus` **are** those two; `laws::e_at`; `API.md` §3's
two rows, `NUMERICS.md` §14's two rows. Supersedes [`0044`](./0044-four-primitives-the-first-consumer-names-and-no-spec-does.md)
decisions 1 and 3; its decision 2 (`SEn3::renormalize`) and item 5 (the `quat_*` corpus debt,
now two routines lighter) stand

## Context

[`0044`](./0044-four-primitives-the-first-consumer-names-and-no-spec-does.md) added four public
items because [`0041`](./0041-the-integration-is-an-adapter-at-tf-tree-math.md) names them. A
review of the PR that landed them asked the question that record's own *Rationale* half-concedes —
"the alternative for each is *let the adapter write it*, and for two of the four that is right" —
and the measurements below say it picked the wrong two.

**`Quat::dot` and `Quat::norm` have no caller in either library crate.** `Quat`'s `w, x, y, z` are
`pub` and `norm_sq` is public, so each is one line in a consumer and neither reaches anything
private. Three further facts:

- **`dot`'s stated purpose is already shipped, better conditioned.** Its rustdoc offered
  `cos(α/2)` for "the angle `α` between the two rotations", for a consumer to read the sign of and
  pick the shorter arc. `x1.rminus(x0)` returns the tangent whose norm *is* that angle, through
  D5's quaternion `atan2` rather than a cosine whose slope vanishes where the angle does. `API.md`
  R3 keeps one spelling per operation.
- **Nothing here will ever read a dot's sign.** `docs/maths/geodesics.md` GE.14 proves the provided
  `geodesic` body already *is* shortest-arc slerp in the `atan2` spelling, so
  [`0045`](./0045-two-phase-4-checks-cannot-be-taken-as-written.md) took slerp off the geodesic
  path. The one caller the item was added for does not exist and is not coming.
- **Both shipped a `# Domain` that R6 cannot enforce.** R6 is NORMATIVE: a restricted domain is a
  `# Domain` section *naming it in `NUMERICS.md` §12* plus a `debug_assert!`. `norm` had neither —
  no §12 row to name, and no assert — and the domain it did state was one-sided: it named the
  overflow at `1.3e154` and not the underflow, so `q = (1e-170, 0, 0, 0)` returns `0.0` for a
  quaternion of norm `1e-170` while satisfying the stated contract word for word.

**`mul_inv` is the opposite case, and as landed it earned none of it.** It is a group operation
with its own rounding, so D1 forbids a consumer writing it, and `0044`'s algebra is right: one
rotation of a vector where `a * b.inverse()` does two. But

- **the library's own relative difference kept paying for two.** `LieGroup::lminus`'s shipped body
  is `(*self * base.inverse()).log()` — the exact composition `mul_inv` replaces — and `SEn3`
  overrode neither `lminus` nor `rminus`. The saving was offered to an external caller and
  declined for the first consumer's hot path, with nothing recorded about why.
- **the direction a tree walk takes was the one not added.** `rminus` is the default side
  (`0002`), its body is `(base.inverse() * *self).log()`, and a relative transform
  `T_w_a⁻¹ · T_w_b` is that shape. `mul_inv` was added because that is the name
  `tf_tree_math::Iso3` happens to spell.
- **`N = 1` was an arbitrary cut.** `0044` decision 3 restricts it to SE(3) and says the algebra
  "holds column by column" for general `N`. It does, in six characters of `array::from_fn`, and
  restricting it is what kept `rminus`/`lminus` — defined for every `N` — from using it.
- **its bound was measured against the wrong denominator.** `laws::e` divides by `max(‖want‖, 1)`.
  Both spellings carry an absolute error of order `u·max‖x‖`, so dividing their difference by the
  *difference* reports how much two nearby frames cancel, not how well the routine rounds. The
  recorded `7.052 u` is a fact about translations drawn in `[-1, 1)`: the same code reads
  `106.32 u` once they reach `10³`.

## Decision

1. **`Quat::dot` and `Quat::norm` are removed.** This supersedes `0044` decision 1. `0041`'s
   Wave 1 keeps `tf_tree_math`'s own bodies for them: they are not group operations, so D1 does not
   bite, and the adapter's `norm_squared`/`rotate`/`normalize` rows are unaffected. Nothing in
   either crate loses a caller, because neither had one.
2. **`SEn3::renormalize` stays, with no change in behaviour.** It closes a real hole — `q` is
   private, so the alternative is a `parts`/`from_parts` round trip around an encapsulation we
   chose — and `0044` decision 2 is reaffirmed as written. Its rustdoc drops the claim that an
   `SEn3` "could not be renormalized from outside at all", which `0044`'s own *Context* contradicts
   in the same words, names the [`Quat::renormalize`] its body actually calls rather than the
   `SO3::renormalize` layer it skips, and restates the domain as that method's: **none**, because
   the step is defined for every input and only its *accuracy* has a domain. That is how R6 is
   satisfied here, and `NUMERICS.md` §12's `renormalize` row already says it.
3. **`mul_inv` is defined at every `N`, gains its mirror `inv_mul`, and `SEn3`'s `rminus` and
   `lminus` *are* those two methods.** This supersedes `0044` decision 3's `N = 1` restriction.

   ```rust
   pub fn mul_inv(&self, other: &Self) -> Self   // self · other⁻¹
   pub fn inv_mul(&self, other: &Self) -> Self   // other⁻¹ · self
   ```

   `inv_mul` is `(q_b* q_a, R(q_b*)(x_{a,i} − x_{b,i}))`: one rotation per column, and the
   subtraction moves *ahead* of the rotation, where for two nearby frames it is exact, instead of
   between two separately rounded rotations. Each is a `NUMERICS.md` §14 row whose reference twin
   is the composition it is not bit-identical to.
4. **A twin bound names the scale it is relative to.** `laws::e_at(a, b, scale)` is `laws::e`
   against a stated denominator, and the two spellings are scored against `max(‖x_a‖, ‖x_b‖, 1)`
   per column, at three translation scales per draw. A bound is a property of the routine only if
   its denominator is one.

## Rationale

The surface shrinks by one item and every item left on it is load-bearing: `renormalize` reaches a
private field, and `mul_inv`/`inv_mul` are the bodies of two trait methods this crate ships. That
is the test `0044` should have applied and did not — not "is it trivial to write" but "does
anything here depend on it".

Decision 3 is also the only part of this with a measurable consequence, and the corpus is where it
was judged (D7, `0006`), not on an operation count:

| `se3_geodesic`, `f64` | before | after | oracle #1 |
|---|---|---|---|
| `geo:consecutive` | 2.556 u | **1.572 u** | 2.34 u |
| `geo:generic` | 3.239 u | **2.721 u** | 2.50 u |
| `geo:near-pi` | 3.466 u | **2.429 u** | 3.25 u |

`SEn3::geodesic` reads `rminus`, so the geodesic ids are where this shows. After the change all
three maxima equal `so3_geodesic`'s to every digit, which is the sharper statement: the translation
block is **no longer the scored maximum on any geodesic stratum**, and what is left is the shared
rotation path. Two of the three strata moved from losing to oracle #1 to dominating it, so the six
`geo:*` rows carry **three** domination failures where they carried five, and `just envelope` goes
from 117 to **115**. The one `se3_geodesic` row still losing, `geo:generic` at 2.721 u against
2.502, is now the same number `so3_geodesic` loses that stratum by — which is §1.2's case, not
this record's.

**Two of the seven geodesic identity legs moved the other way**, and the record says so rather than
quoting the favourable half: `symmetry` 13.352 → 16.111 u (`f32`, `N = 1`) and `left` 9.775 →
10.731 u (`f32`, `N = 2`), against `velocity` 8.353 → 7.768 and `t=1` 6.591 → 6.578, with `f64` and
`Dual` at `N = 1` improving on every leg. The legs compare *two calls with each other* —
`γ(x₀, x₁, t)` against `γ(x₁, x₀, 1−t)`, `γ(h x₀, h x₁, t)` against `h γ(x₀, x₁, t)` — so they
measure how well two roundings line up, and `inv_mul` subtracts inside a different rotation on each
side of the swap. The corpus measures error against a 110-digit reference. `0006` makes the corpus
the bar, on the max, per stratum; the legs are consistency checks whose bounds are re-recorded at
twice the new worst (`[3, 14, 33, 16, 22, 19, 3]`).

An alternative considered and rejected: keep `dot` and `norm` and give them a §12 row, an assert
and a corpus stratum. That is the full price of a public numeric routine (D7) paid for two
expressions with no caller, and `0044` item 5 had already deferred the stratum to `0.1.0`, which
means they would have sat unverified through a release. Removing them discharges that debt instead
of dating it.

## Consequences

- Public surface: **−2** (`Quat::dot`, `Quat::norm`), **+1** (`SEn3::inv_mul`), and `mul_inv`
  widens from `SE3` to every `SEn3`. `cargo-semver-checks` gates two fewer items from 1.0.
- `0044` item 5's `quat_*` corpus debt loses `norm` and `dot` and keeps `renormalize` and
  `SEn3::renormalize`, still due before `0.1.0`.
- `SEn3`'s `rminus`/`lminus` change their bits. Nothing is published, no baseline is blessed yet,
  and `conformance/results/` is not committed, so the only recorded figures to move are the
  `laws::Bounds` ones above and `se3_geodesic`'s three strata.
- `0041`'s Wave 1 list keeps `Quat::dot` and `Quat::norm` as `tf_tree_math`'s own; its Wave 2
  `Iso3::mul_inv` row now has a mirror to delegate `rminus`-shaped calls to.
- `laws::e_at` is available to any other twin whose answer can cancel. `laws::e` stays the default
  and `laws::e_at`'s rustdoc says when it is the wrong question.
- **`Bounds::plus_minus` was tight by luck, independently of this record.** Re-measuring it on a
  `laws::Rng` stream reads 9.666 u at this record's parent commit, where the committed bound was
  16 — twice the 7.561 u a 10^6-case *proptest* run had found. Two samplers of one distribution
  give different maxima, so a bound set from one of them is twice the luckier sample, not twice the
  worst, and `0006`'s convention ("the sample count is part of the bound") is what that violates.
  The bound is now 21, from 10.346 u, and `measure_plus_minus` makes the figure rerunnable the way
  `measure_geodesic` already made the legs'. This is reported separately from the 1.07x the
  decision itself cost so that neither hides behind the other.
- **`just lint` gained a `[symbols]` check**, because two of the citations this record corrects
  would have been caught by the rule the project already has. "Cite a symbol, never a line number"
  is only half a rule while a rename passes silently, and a stale citation is worse than a line
  number: `grep` returns nothing, so the reader cannot tell whether the claim moved or was deleted.
  A `<module>_tests::<name>` citation now resolves against that module's file. It found a **third**
  stale citation on its first run, in `laws::geodesic_legs`, which no review had named.

## Implementation plan

1. This record, `API.md` §3's two rows, `NUMERICS.md` §14's two rows, `0041`'s Wave 1 note —
   verified by `just lint`.
2. `Quat::{dot, norm}` and their test removed — verified by `just test` and `just doc`.
3. `mul_inv` at every `N`, `inv_mul`, the two `rminus`/`lminus` overrides, `laws::e_at`, the
   scale-relative twin proptest at `f64`/`f32`/`Dual` and `N ∈ {1, 2}` with its measured bound,
   and the quaternion/column split test — verified by `just test` and by
   `cargo xtask conformance --fn se3_geodesic` reading the table above.
4. `Bounds::geodesic` and `plus_minus` re-recorded from their `#[ignore]`d measurements — verified
   by `just test` at the 10^6-case protocol.

## Open questions

None.

## Further work

1. **`mul_inv`/`inv_mul` as corpus ids.** They are twinned and bounded but unscored; `0044` item 5
   is where that debt lives and this record adds two routines to it whose arithmetic is now on
   every `rminus`.
2. **`symmetry` and `left` at `f32`.** The two legs this cost are both swap/conjugation
   identities at binary32 and neither is the accuracy bar, but neither is explained to the digit
   either. `0038`'s paired instrument is the tool; it is not a blocker for the surface.
3. Whether `Product`'s `rminus`/`lminus` owe the same treatment. Its factors' overrides are
   already used, so an `SEn3` factor gets this for free; a `Product<SO3, Rn<3>>` under the SE(3)
   reading does not, which is `0045`'s GE.5(b) territory and not this record's.
