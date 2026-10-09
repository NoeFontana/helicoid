# 0060: the charts go first, and name their frame

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #118, #119, #120, #123

## Context

locus-tag's pose solver updates its iterate with `Pose::retract`: `(R Exp(ω), t + R v)`, which
is `Decoupled` (`PHASE5.md` §1.3). It is hand-written on nalgebra and has no `local` and no chart
Jacobian. Every one of its LM loops (`nielsen_lm`, the board solvers, `refine_aw_lm`) calls it,
and the Jacobians those loops form are checked only by finite differences through it. `0012`
names its target, `WithChart<SE3<f64>, Decoupled>`, and `0049` moves `PHASE5.md` §1 onto the
locus-tag path ahead of the rest of Phase 4. No spec gates §1, §2, §4 or §5 on Phase 4 or on
`PHASE3.md` §6 (SO(2), SE(2)): the only gate in `PHASE5.md` is §3's on `NUMERICS.md` §9. Still,
`PROJECT.md` §4 lists the phases in sequence, and Phase 5's appendix opens with the charts without
saying what may run beside them.

Every formula the charts need is derived and checked in `docs/maths/charts.md` (CH.1–CH.7, at
60 to 110 digits). Two things stand between those results and code.

**A sentence that is wrong.** `0012`'s Context says the three SE(3) retractions "agree to first
order and differ at second". `Screw` and `Decoupled` do: they have equal `D ret(0)`, and their
trial points differ by `R(J_l(φ) − I)ρ`, of norm at most `θ/2 ‖ρ‖` (CH.5(c), CH.7(c)).
`WorldTranslation` does not, unless `R = I`. Its `ρ` is a world-frame translation, so
`DΦ^{Dec→WT}(0) = diag(I, R)`, and for the same `δ` the retracted translations differ at first
order, by `2 sin(θ(R)/2) ‖ρ‖` (CH.5(c)). What is true is that it equals `Decoupled` exactly after
that linear map (CH.4(b)). `PHASE5.md` §1.3 carried the same sentence until #104 removed it;
`docs/maths/index.md` still lists it as open.

**Four readings `PHASE5.md` §1 leaves to the implementer.** A spike (rustc 1.87, the MSRV, and
stable; a 70-line crate reproducing `Tangent`, `Jac`, `LieGroup`, `ProductJac` and §1.1's traits
as declared, not committed) settles what the language allows.

1. *Where the base lives.* `Chart::at(base) -> Self` stores the base, so a chart type is generic
   in `S`. §1.2's `WithChart<SE3<f64>, Decoupled>` cannot be written as is.
2. *"Every `LieGroup` is a `Manifold`".* A blanket `impl<S, G: LieGroup<S>> Manifold<S> for G`
   beside `impl Manifold<S> for WithChart<M, C>` is **E0119**: a downstream crate may write
   `impl LieGroup<TheirReal> for WithChart<..>`, so coherence cannot rule out the overlap.
3. *`WithChart`'s chart.* If `C` implemented `Chart<S, WithChart<M, C>>` as well as `Chart<S, M>`,
   every concrete `c.retract(..)` would have two candidate impls and need a type annotation.
4. *§1.4's `ProductJac<Mat3<S>, Mat3<S>>`.* That type is a `Jac<S, (SO3Tangent, SO3Tangent)>`
   through `product.rs`'s blanket impl, not a `Jac<S, Twist<S>>`. Adding the second impl compiles,
   but then every concrete call whose signature does not name the tangent (`inverse`, `mul`,
   `neg`, `identity`, `write_dense`, `sandwich`) is **E0283** (type annotations needed). A consumer
   would write `Jac::<f64, Twist<f64>>::write_dense(&j, ..)`.

## Decision

1. **Order.** `PHASE5.md` §1, §2, §4 and §5 may proceed beside `PHASE4.md` §2 and `PHASE3.md` §6.
   §3 keeps its gate on `NUMERICS.md` §9. §1 comes first: it is the locus-tag dependency, and S²
   (§2) and `Gaussian` (§5) are written against `Chart` and `Jac`. `PROJECT.md` §4 gains one line
   saying so.
2. **The correction.** `0012`'s Context gets an amendment banner citing this record: `Screw` and
   `Decoupled` agree to first order; `WorldTranslation` agrees with them after `diag(I, R)` on the
   translation tangent, which is exact, not first order. `0012`'s Decision does not change: three
   charts, `Decoupled` for locus-tag. `docs/maths/index.md`'s open item is marked answered by
   `0060`, and CH.5(c)'s parenthetical points here.
3. **Chart types carry `S`.** In `crates/helicoid/src/chart.rs`:
   - `Chart<S, M>`, `Manifold<S>` and `WithChart<M, C>` as `PHASE5.md` §1.1, plus
     `WithChart::new(m)`, and one provided method on `Chart`,
     `local_with_jacobian(&self, other) -> (Tangent, Jac)`: a residual and its Jacobian, the pair
     every Gauss–Newton residual in a chart needs. The default calls the two methods. Every chart
     here overrides it so that the relative element and its `Log`, the most expensive part, are
     formed once, and the result is bit-identical to the two calls. `Clone`, `Copy` and `Debug` are written on `M` alone, so `C` need not be
     `Copy`.
   - `RightChart<G>(G)` and `LeftChart<G>(G)`, the frozen base. `RightChart`: `retract = X Exp(δ)`,
     `local = Y ⊖_R X`, `rj = J_r(δ)`, `lj = J_r⁻¹(Y ⊖_R X)`. `LeftChart`: `Exp(δ) X`, `Y ⊖_L X`,
     `J_l`, `J_l⁻¹` (CH.3). Their `Tangent` and `Jac` are the group's.

   In `crates/helicoid/src/se3_charts.rs`: `Screw<S> = RightChart<SE3<S>>`, and
   `Decoupled<S>(SE3<S>)` and `WorldTranslation<S>(SE3<S>)`, newtypes over the frozen base. Their
   `retract`, `local`, `rj` and `lj` are the table and CH.5's matrices, built from `SO3::{rplus,
   rminus, act, jr, jr_inv, to_matrix}`. `Decoupled::local` is `SEn3::inv_mul`'s `X⁻¹Y` read as
   `(Log R_rel, t_rel)`: one conjugate, and the routine's own reference twin. Like the CH.5 note, they use `exp_coeffs`, `log_ratio`,
   `jr_coeffs` and `jr_inv_coeff` through those methods and evaluate no coefficient themselves.
   §1.2's example becomes `WithChart<SE3<f64>, Decoupled<f64>>`.
4. **One `Manifold` impl per group.** `SO3<S>`, `SEn3<S, N>`, `Rn<S, N>` and
   `Product<A, B>` (for `A`, `B` groups) each get a five-line impl: `Tangent` and `DOF` the
   group's, `Chart = RightChart<Self>`. A comment at the impls cites decision 4 of this record for
   why there is no blanket. A group added later adds its impl in the same PR (`API.md` §6).
   With `LieGroup` and `Manifold` both in scope, `SE3::<f64>::DOF` is ambiguous (E0034).
   `Manifold`'s rustdoc says to write `<SE3<f64> as LieGroup<f64>>::DOF` there.
5. **`WithChart`'s chart is `Lifted<C>`**, a `#[repr(transparent)]` newtype implementing
   `Chart<S, WithChart<M, C>>` by delegating to `C: Chart<S, M>` and wrapping or unwrapping `.0`.
   It is the only `Chart` impl over `WithChart`, so calls infer. `Manifold for WithChart<M, C>`
   has `Tangent = C::Tangent`, `Chart = Lifted<C>` and `DOF = C::Tangent::DOF`, and `Blend` is
   `M`'s.
6. **`TwistBlockJac<S>`** is `#[repr(transparent)] pub struct TwistBlockJac<S>(ProductJac<Mat3<S>,
   Mat3<S>>)`, implementing `Jac<S, Twist<S>>` and nothing else. `mul`, `inverse`, `neg`,
   `identity`, `write_dense` and `sandwich` delegate to the pair impl: the same dense layout,
   rotation block first, so the same bits. `apply` and `apply_transpose` split the `Twist` into
   `phi` and `rho[0]`. `Decoupled` and `WorldTranslation` have `Tangent = Twist<S>` and
   `Jac = TwistBlockJac<S>`, so all three SE(3) charts take the same `δ`. `PHASE5.md` §1.4's
   `ProductJac<Mat3<S>, Mat3<S>>` becomes this newtype. Two accessors, `rotation_block()` and
   `translation_block()`, return the two `Mat3`s.
7. **The transition at δ = 0.** A sealed trait `Se3Chart<S>: Chart<S, SE3<S>, Tangent = Twist<S>>`
   with `const WORLD_TRANSLATION: bool` is implemented by `Screw<S>` (`false`), `Decoupled<S>`
   (`false`) and `WorldTranslation<S>` (`true`). `LeftChart<SE3<S>>` is not: its transition is
   `Ad_X` (CH.6). One inherent method:

   ```rust
   impl<S: Real> SE3<S> {
       /// `DΦ^{From→To}(0)` at this base (CH.4(b), CH.5(c)): a covariance in `From`'s
       /// coordinates moves to `To`'s as `j.sandwich(&cov)` (CH.6).
       pub fn chart_transition<From: Se3Chart<S>, To: Se3Chart<S>>(&self) -> TwistBlockJac<S>;
   }
   ```

   It returns `diag(I, R)` from a body-frame chart to `WorldTranslation`, `diag(I, Rᵀ)` the other
   way, and `diag(I, I)` otherwise, written as `Mat3::identity()`, so `Screw → Decoupled` is
   exactly `I`. `R` is `self.rotation().to_matrix()`.
8. **Corpus ids.** `se3_{screw,decoupled,world}_{retract,local}`, six ids, binary64 only, as every
   vector id is until a record extends `0016`.
   - `retract`: inputs are a base `X` and `δ`. The strata are `SEN3_STRATA` (the `theta:*` and
     `rho:*` families) on `δ`, and the base is drawn generic from the stratum's own stream: a
     unit quaternion uniform on `S³` and `t` a uniform direction at unit norm, both rounded to
     binary64. The reference is the table of `charts.md` §3 at 60 digits, with
     `Exp` and `J_l` by their series.
   - `local`: inputs are `X` and `Y = ret_X(δ)` rounded to binary64, over the same strata. The
     reference uses the geometric `Log` of the quaternion (`atan2`, `0045`), not `mp.logm`, which is
     complex from `θ ≈ 3.03` (`docs/maths/index.md`). `Y` is each chart's own retract, so a
     `local` record's answer is its stratum's `δ` up to the rounding of `Y`.
   - The metric is `NUMERICS.md` §11's forward error, as `se3_geodesic`'s: a `retract`'s
     quaternion sign-aligned with floor 1 and its translation with floor `‖x₀‖`. A `local` is held
     as `phi` and `rho`, with floors 1 and `‖x₀‖`. Its inputs are two absolute poses, so a nearby
     pair cancels and the answer's own size is the wrong denominator: scored relative, the
     `rho:1e-6/theta=1e-8` stratum reads `2e6 u` from `q₀* q₁`'s absolute `u`, which is the
     problem's conditioning and not the routine's. No oracle answers these ids, so their bars are
     no-regress and coverage.

   The ids go into `PHASE1.md` §4.3's definitions table and `xtask/src/envelope/coverage.rs`'s
   `required()` in the PR that adds them. `Screw`'s pair duplicates `sen3_exp_n1`/`sen3_log_n1`
   composed with a base, which is what a consumer calls, and it costs one generator function.

## Rationale

- **The newtype, not the second impl.** Decision 6 keeps §1.4's intent, a block-diagonal `Jac`
  whose blocks are `ProductJac`'s, and removes the E0283 cost from every consumer call. Two
  alternatives lost:
  - A product tangent per chart (`(SO3Tangent, RnTangent<S, 3>)` with `ProductJac<Mat3, RnJac>`)
    cannot express `Decoupled`'s translation block `Exp(−φ)`, since `RnJac` is `kI`. Making `Mat3`
    a `Jac` over `RnTangent<S, 3>` as well moves the ambiguity onto `Mat3` itself.
  - A different tangent type per chart would give `δ` a different type in each chart, and
    `chart_transition` and a consumer's switch between charts would both need conversions.

  One tangent for the three charts is what makes the chart a choice visible only in the type, which
  is `0012`'s point.
- **Per-group impls, not a blanket.** The blanket does not compile beside `WithChart`. A macro
  saves four five-line impls and hides the place a new group must be added. Dropping `WithChart`'s
  `Manifold` impl instead would undo `0012` decision 3.
- **`Lifted<C>`, not a second `Chart` impl on `C`.** The latter makes every concrete chart call
  need an annotation, for no gain.
- **A public transition.** locus-tag carries covariances between body and world frames today
  (`covariance_body_to_camera`) and would otherwise write `diag(I, R)` by hand. That is the frame
  mistake `0049`'s test asks the boundary to remove. A const on a sealed trait follows `Side`.
- **Ids for `Decoupled` and `WorldTranslation`.** `API.md` §6 item 2 accepts a stratum or a twin.
  These are compositions of corpus-checked SO(3) primitives with no cancelling coefficient and no
  fast path, so a twin would restate the definition. The ids cost one generator function each and
  are what locus-tag's migration gate (`PHASE5.md` §6) will compare against.

## Consequences

- `PHASE5.md` §1.1–§1.4 and §0.0, `PROJECT.md` §4, `API.md` §2 (`Chart`, `Manifold`, `WithChart`,
  `Lifted`, `RightChart`, `LeftChart`, `Screw`, `Decoupled`, `WorldTranslation`, `TwistBlockJac`,
  `Se3Chart`, `SE3::chart_transition`) and `PHASE1.md` §4.3 change as above.
- The SE(3) charts' types are part of the API; a fourth chart is a record (`0012` decision 4).
- Every chart is checked against `Dual` (`laws::chart_jacobians_match_dual`, step 2 below), like
  the group Jacobians.

## Implementation plan

1. This record and the documentation edits of *Consequences*. Verified by `just lint` and
   `just doc`.
2. `chart.rs`: decisions 3 (group charts), 4 and 5. `laws::chart_jacobians_match_dual<S, M, C, D>`
   checks `rj` as the `Dual` lanes of `C::at(&y).local(&c.retract(&δ.add(&η)))` and `lj` as those of
   `c.local(&C::at(&y).retract(&η))`, with `y = c.retract(&δ)` built from a constant `δ` so that
   only `η` carries lanes (CH.1). `laws::transition_dual` takes `DΦ(δ)` from
   `C2::at(x).local(&C1::at(x).retract(&(δ + η)))`. `chart_tests.rs` runs a `chart_laws_for!` macro
   over `SO3`, `SEn3<_, 1..=3>`, `Rn<_, 3>` and `Product<SO3, Rn<_, 3>>`, right and left:
   `retract_zero_is_base`, `local_retract_roundtrip` and `retract_local_roundtrip` (θ ≤ √3 < π,
   inside `U₀`), `chart_jacobians_match_dual`, and `right_to_left_change_of_chart` (`DΦ = Ad_X`,
   `Ψ = Ad_Y`, CH.6), plus `with_chart_is_bit_identical_to_its_chart`. The bounds are measured by
   `#[ignore]` tests and recorded as `laws::Bounds` does. Verified by `just test`, `just lint`,
   `just doc`, `just no-std` and `just msrv`.
3. `se3_charts.rs`: decisions 6 and 7, and `se3_charts_tests.rs`:
   - `chart_jacobians_match_dual_{screw,decoupled,world_translation}`;
   - `chart_transition_matches_dual_{from}_{to}` over the nine ordered pairs (CH.4(b));
   - `chart_transition_screw_decoupled_is_identity`, checked bit for bit;
   - `change_of_chart_{screw_decoupled,decoupled_world_translation,screw_world_translation}`, both
     identities of CH.6;
   - `twist_block_jac_is_the_pair_bit_for_bit`, and its `sandwich` against
     `reference::productjac_sandwich`.

   Verified as step 2.
4. The six ids of decision 8. Verified by `just corpus-check`, `just conformance` and
   `just envelope` (no `--bless` without the owner).
5. `PHASE5.md` §0.0 rows 1 and 2 set to Done, and the `Implementation:` lines here and on `0012`.

## Open questions

None.

## Further work

- CH.7(d)'s LM illustration is not a test. A trial-point comparison between `Screw` and
  `Decoupled` at a fixed `δ` would pin the `θ/2 ‖ρ‖` bound in code.
- `0012`'s Context reports a `JᵀWJ` conditioning improvement for `Decoupled` without saying
  against which chart, and CH.7(a) shows `Screw` and `Decoupled` have the same `J` at the base. If
  locus-tag's rustdoc carries the claim, its migration is where to measure it.
