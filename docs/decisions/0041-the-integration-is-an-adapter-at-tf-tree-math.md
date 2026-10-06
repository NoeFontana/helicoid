# 0041: the integration is an adapter at `tf_tree_math`, not a type substitution

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** steps 2 and 4's `helicoid` halves have landed (`SO3`, `SEn3`); the rest is open

> **Amended 2026-10-06 by [`0049`](./0049-the-boundary-is-what-removes-a-way-to-be-wrong.md), which is the signed edit the index's rules allow.** Decision 4's
> Wave 1 and Wave 2 lists lose two further items (`Quat::normalize`, `Iso3::normalized` — a
> *different* operation from `Quat::renormalize`, not a non-group one), decision 4's
> `quat_from_rot3` row is a recorded behaviour change rather than a port, and **decision 7's
> publish moves ahead of the remainder of plan step 4**: every Wave 1 and Wave 2 item maps onto
> shipped surface, so the publish is the only thing between here and a consumer compiling against
> this crate. The decision this record takes is unchanged; what changed is its order and three rows
> of its inventory.

## Context

`PROJECT.md` §4 names Phase 4 "geodesics and the `tf_tree` migration" and D1 says `helicoid`
replaces `tf_tree_math`'s SE(3), so the roadmap reads `tf_tree`-first. An evaluation (2026-10-02)
read three constraints as blockers and proposed locus-tag first instead; **no record carried
either the evaluation or the reordering.** Re-examined against the code, with `tf_tree` owned by
the same owner, two of the three are not blockers, and they were artefacts of one assumed plan —
substituting `helicoid`'s types for `tf_tree_math`'s.

- **The same-item assertion does not constrain the backend.** `tf_tree`'s
  `tests/math_reexports.rs` uses `same_item` to reject *a facade wrapper or a missing re-export*;
  its own doc comment says so ("What can fail is a missing re-export or a facade wrapper"). It
  constrains `tf_tree` against `tf_tree_math` and says nothing about what `tf_tree_math` is
  implemented on. Replacing `tf_tree_math`'s **bodies** leaves every assertion passing.
- **`Pod` is what D2 provides for, not what it forbids.** `Iso3`, `Quat`, `Vec3` and `Twist` are
  `repr(C) + Pod` and `PoseSlot` is `[AtomicU64; 7]` through `Iso3::to_bits`. D2's own words are
  "consumers own their storage formats (`tf_tree`'s arena records stay `tf_tree`'s)". The conflict
  exists only if `helicoid`'s types must *become* the stored ones.
- **The conventions already agree, exactly.** `tf_tree_math` is Hamilton `w`-first, its tangent is
  `[ω, v]` (rotation-first) and its perturbation is `T̂ · exp(ξ^)` (right) — `0002`'s three
  choices, with no conversion anywhere. Contrast locus-tag, which publishes a translation-first
  covariance needing `Π Σ Πᵗ` that no API item provides: on conventions `tf_tree` is the *easier*
  consumer, not the harder one.
- **The budget is already `helicoid`'s.** `tf_tree_math`'s dependencies are `libm` and `bytemuck`
  only — `0007` plus the storage crate. No `nalgebra`, no `alloc` to remove.
- **The types line up field for field.** `Iso3 { q: Quat { w, x, y, z }, t: Vec3 }` is seven
  contiguous `f64`; `PHASE3.md` §5 specifies `SEn3<S, N> { q: Quat<S>, x: [Vec3<S>; N] }` with
  `from_quat_translation`, `rotation()` and `translation()`. The conversion is a safe seven-field
  copy in the same order, no transmute and no `unsafe`.

What is genuinely constraining:

1. **Release sequencing.** `tf_tree` is published at `0.0.6`; `helicoid` is unpublished `0.0.x`
   where every release may break every other (`CLAUDE.md`, Status). A published crate cannot take
   that as a hard dependency.
2. **Latency, and one published constant derived from it.** `PHASE1.md` §9 gates a depth-3 hot
   lookup at p50 **under 300 ns with `ScLerp`** *and* **within 25 % of the committed baseline per
   percentile**; the committed median of `lookup/depth3/sclerp` is **192.7 ns**, so the second
   clause is the tighter one (≈241 ns, not 300). And `NS_PER_STEP_ESTIMATE = 64` in
   `tf_tree_py::tree` is *derived from that same benchmark* (192.7 / 3 steps) and published through
   `tf_tree`'s `API.md` §3.4, so a material change to it is a published-API change.
3. **~~`helicoid` has no group type.~~ Cleared.** When this record was drafted `PHASE3.md` §0.0
   read `SO3` "not started", and every item below landed on a type that did not exist; that was
   the keystone. `PHASE3.md` §0.0 now reads "**`SO3` done**" and "**`SEn3<S, N>` done**, with
   `SE3`/`SE23`", each scored over its corpus ids, so **steps 2 and 4's `helicoid` halves are
   complete** and the live front is step 3. What is left on this side is Phase 4 (§1–§5), which
   the implementation plan below now names.

## Decision

1. **The integration is an adapter at `tf_tree_math`.** It keeps every public type (`repr(C) +
   Pod`), every signature, `to_bits`/`from_bits` and its `bytemuck` dependency; its function
   bodies delegate to `helicoid`. The `tf_tree` facade is untouched and its `same_item` tests keep
   passing. **`helicoid` adds `Pod`, `Zeroable` or `serde` to no value type** (D2, immovable).
2. **Behind a non-default `helicoid` feature on `tf_tree_math`**, so the published default path is
   unchanged while `helicoid` is `0.0.x`. Both paths build in CI.
3. **That feature is the parity instrument.** `helicoid` already treats `tf_tree_math` as oracle #1
   (D15, `0010`) with a runner and `just oracle-tf-tree-math`, so parity is measured per function
   over the committed corpus, in units of `u`, by machinery that exists. This is what produces
   `PHASE4.md` §5's parity table rather than a separate exercise.
4. **Three waves, ordered by what `helicoid` ships**, each a PR on each side:
   - **Wave 1 — SO(3)** (`PHASE3.md` §4): `exp_so3`, `log_so3`, `quat_from_rot3`, `slerp`,
     `Quat::rotate`. Four of the nine items first listed here stay `tf_tree_math`'s own bodies, for
     two different reasons: `dot` and `norm` are not group operations and D1 does not reach them
     ([`0048`](./0048-the-relative-transform-pair-earns-the-surface-dot-and-norm-do-not.md)), and
     `normalize` with `norm_squared` is a *different operation* from the one `helicoid` ships —
     the exact projection `q/‖q‖` against `renormalize`'s single Newton step, which is accurate
     only inside the drift band ([`0049`](./0049-the-boundary-is-what-removes-a-way-to-be-wrong.md) decision 3). `quat_from_rot3` delegates and **changes its
     output**: it returns an un-normalized quaternion where `SO3::from_matrix` normalizes, so that
     row of the parity table is a recorded behaviour change ([`0049`](./0049-the-boundary-is-what-removes-a-way-to-be-wrong.md) decision 4).
   - **Wave 2 — SE(3)** (`PHASE3.md` §5): `exp_se3`, `log_se3`, `Iso3::{inverse, mul_inv}` and
     `Mul`, `Twist::{adjoint, adjoint_inv, to_spatial}`. `Iso3::normalized` is Wave 1's
     `normalize` one level up and stays with it.
   - **Wave 3 — the screw path** (`PHASE4.md` §1.2's dual-quaternion fast twin): `sclerp`,
     `screw_pow`, `screw_twist`, `screw_pow_with_twist`, `ScLerp`, `LerpSlerp`.
   **The whole latency risk is in Wave 3**: Waves 1 and 2 prove the mechanism on paths the 300 ns
   gate does not time. Wave 3 is *not* the only wave needing Phase 4 — the first draft of this
   record said so and its own Wave 1 contradicts it, because `slerp` is `SO3::geodesic`
   (`PHASE4.md` §5.1) and `PHASE4.md` §0.0 reads "Not started". Wave 1 therefore depends on
   `PHASE4.md` §1.1's provided `geodesic`, which is one expression over `rplus`/`rminus`/`scale`
   and needs no fast twin: SO(3) does not override it, because `docs/maths/geodesics.md` GE.14
   proves the provided body **is** shortest-arc slerp, and in the better-conditioned spelling
   (`atan2`, where `tf_tree_math::slerp` has an `acos` whose slope is infinite at π). The *screw*
   twin, which is a genuine Phase 4 implementation, stays Wave 3's alone.
5. **Consumer order: `tf_tree` first, locus-tag second.** locus-tag's blockers (a translation-first
   covariance with no API item, `0029`) are real but are *its* blockers; `tf_tree` needs no
   convention conversion. The covariance API is deferred by the owner and is not on this path.
6. A per-function fallback is kept: a function that regresses on accuracy or latency stays on the
   old path, and the feature is not all-or-nothing.
7. **`helicoid` and `helicoid-linalg` publish at `0.0.1` before Wave 1 lands on `tf_tree`'s
   `main`.** This was open question 1 and it is not optional: `tf_tree_math` is published at
   `0.0.6`, and `cargo publish` rejects a manifest whose dependency — optional or not — has a
   path or git source. So decision 2's feature flag removes the *semantic* dependency on
   `helicoid`'s stability but not the *packaging* one. `0.0.x` costs nothing semantically
   (`CLAUDE.md`'s Status: every release may break every other), and the alternative is a
   long-lived `tf_tree` branch carrying `[patch.crates-io]`, which defers this integration's one
   real risk — the latency gate — to the very end. Publishing `0.1.0` and making the feature
   default is a separate, later decision and is not taken here.

## Rationale

The alternative the first evaluation assumed — **type substitution**, `helicoid`'s `SE3` replacing
`Iso3` — is what produced all three "blockers": it breaks the facade's `same_item` tests, it needs
`Pod` on a `helicoid` type against D2, and it is a big-bang change to a published ABI. Nothing
required it. An adapter at the lowest layer keeps the ABI, the storage contract and the facade
fixed, and reduces the integration to "do these two implementations agree", which is a question
`helicoid`'s instrument is already built to answer.

**locus-tag first** stays viable but is strictly harder at the boundary: it needs a tangent-order
conversion, and a covariance crossing that boundary has no API item at all. Doing `tf_tree` first
also exercises `helicoid` against the oracle it was seeded from (`0010`), which is the comparison
with the most existing machinery behind it.

**Waiting for `helicoid` to publish** serialises every consumer behind a publish decision that
Phase 3 is not ready to make; the feature flag removes that dependency without pre-empting it.

## Consequences

- `tf_tree_math` carries two implementations of its math under a feature until the old path is
  deleted: one CI matrix entry and a maintenance cost, bounded by the wave it is in.
- `NS_PER_STEP_ESTIMATE` is re-derived if Wave 3 moves the median, and that is a `tf_tree`
  published-API change (its `API.md` §3.4 row 10), not an internal one.
- `tf_tree_math` keeps `bytemuck`; `helicoid` never sees it. D2 is satisfied by construction.
- The defensible claim from this migration is verified numerics against the conformance bars and
  D16 determinism — **not** speed. The gate is no-regress, and `sclerp` has ≈25 % of headroom.
- `PROJECT.md` §4's Phase 4 line is edited to separate the geodesic work `helicoid` owes regardless
  from the `tf_tree` waves, once this record is `ready`.

## Implementation plan

1. This record, and its twin in `tf_tree`'s own `docs/decisions/` — verified by `just lint` and by
   the row in `decisions/README.md`.
2. ~~`helicoid` `SO3`~~ **landed**, with `0029`'s `from_parts`/`parts`; every `so3_*` corpus id is
   scored. ~~`helicoid` `SEn3`/`SE3`~~ **landed** likewise (what was step 4's first half).
3. **The parity instrument before the parity claim** (D7): the two `PHASE4.md` §4 geodesic corpus
   ids, and the `tf_tree_math` runner extended from four ids to the rows §5.1 actually names —
   `so3_from_matrix`, `so3_act`, `sen3_ad_n1` and the two geodesics — verified by `just
   oracle-tf-tree-math` writing them and `just envelope` pairing them.
4. **`helicoid` Phase 4 §1 and §3**: the provided `geodesic`/`geodesic_velocity`,
   `reference::geodesic` and the invariance tests — verified by `just test` and the new corpus ids
   under `just conformance`. ~~`geodesic_jacobians`~~ was step 4's third item and is **step 8**:
   [`0049`](./0049-the-boundary-is-what-removes-a-way-to-be-wrong.md) decision 5 found that no `tf_tree_math` function takes a geodesic Jacobian, so it gates
   no wave and was serialising the integration behind work the consumer does not call.
5. `tf_tree_math`'s `helicoid` feature scaffold and **Wave 1** behind it, after decision 7's
   publish — verified by `just oracle-tf-tree-math` parity on the `so3_*` ids and `tf_tree`'s own
   reference twins passing under both features.
6. **Wave 2** — verified by the `sen3_*` ids and `tf_tree`'s bench gate showing no percentile
   regression.
7. `PHASE4.md` §1.2's dual-quaternion twin and **Wave 3** — verified by the 300 ns / 25 % gate, and
   a re-derived `NS_PER_STEP_ESTIMATE` if the median moves.
8. `0043`'s `geodesic_jacobians` and the first blessed envelope, which gate no wave of this
   integration and serve Phase 5's consumers — verified by `just test` and `just envelope --bless`.

Steps 3 and 4 are `helicoid`'s alone and gate every later step, and both have **landed** (#92, #93,
#94); §5.2's domination precondition is measured against oracle #1 as part of step 3, and what it
finds is settled before step 5, not after — with three geodesic strata still losing to oracle #1
and [`0049`](./0049-the-boundary-is-what-removes-a-way-to-be-wrong.md) *Further work* 1 holding what the split by stratum says about them.

## Open questions

None. Both of the draft's questions are answered above or moved below, per
[`0040`](./0040-a-draft-is-not-a-parking-space.md) item 1.

- *Does `helicoid` publish?* — decision 7. It has to, and at `0.0.1`, not `0.1.0`.
- *Do `quat_from_rot3` and `SO3::from_matrix` agree on their reading?* — this is not a question to
  settle by reading two sources. `so3_from_matrix` is a committed corpus id that the
  `tf_tree_math` runner does not yet answer although the pin exports the function, so the
  comparison is **one runner row away** and is measured rather than argued. Wave 1's scope does
  not depend on the answer: either the rows agree, or the disagreement is a parity row with a
  number on it, which is what §5.1 is for. Scheduled with the geodesic instrument, ahead of
  Wave 1.

## Further work

1. locus-tag's covariance-order gap: a translation-first `[t, ω]` 6×6 covariance needs `Π Σ Πᵗ`
   (`docs/maths/se3.md` SE.14(d)), which no API item provides (`docs/maths/gamma-gaussian.md`
   GG.13) — the one error in GG.13's NEES table with no counterpart in the type system. Deferred by
   the owner; an `API.md` §6 item and its own record when locus-tag comes up.
2. `mint` covers `Vector`, `Point` and `Matrix` to 4×4 (`PHASE2.md` §7), so locus-tag's 6×6 and 8×8
   blocks are not covered by it. Which conversion they take is unanswered.
3. Whether `0030`'s consumer comparison is this migration's gate or a separate trigger.
4. `tf_tree_math`'s `reference.rs` twins and `helicoid`'s §14 twins overlap once the backend is
   shared: whether both are kept (D6 says `helicoid`'s are) or `tf_tree_math`'s become the adapter
   test is unanswered.
