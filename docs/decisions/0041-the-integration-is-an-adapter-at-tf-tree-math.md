# 0041: the integration is an adapter at `tf_tree_math`, not a type substitution

**Status:** draft (awaiting a decision)
**Owner:** @NoeFontana
**Implementation:** none yet

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
3. **`helicoid` has no group type.** `PHASE3.md` §0.0: `Quat` is done, but `SO3` is "not started",
   and `SEn3`, `SE3`, `SE23` and "everything of §5 that needs the group" with them. Every item
   below lands on a type that does not exist yet. This is the keystone, and it is unchanged by
   anything in this record.

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
     `Quat::{normalize, rotate, dot, norm, norm_squared}`.
   - **Wave 2 — SE(3)** (`PHASE3.md` §5): `exp_se3`, `log_se3`, `Iso3::{inverse, mul_inv,
     normalized}` and `Mul`, `Twist::{adjoint, adjoint_inv, to_spatial}`.
   - **Wave 3 — the screw path** (`PHASE4.md` §1.2's dual-quaternion fast twin): `sclerp`,
     `screw_pow`, `screw_twist`, `screw_pow_with_twist`, `ScLerp`, `LerpSlerp`.
   **The whole latency risk is in Wave 3**, which is also the only wave needing Phase 4: Waves 1
   and 2 prove the mechanism on paths the 300 ns gate does not time.
5. **Consumer order: `tf_tree` first, locus-tag second.** locus-tag's blockers (a translation-first
   covariance with no API item, `0029`) are real but are *its* blockers; `tf_tree` needs no
   convention conversion. The covariance API is deferred by the owner and is not on this path.
6. A per-function fallback is kept: a function that regresses on accuracy or latency stays on the
   old path, and the feature is not all-or-nothing.

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
2. `helicoid` `SO3` (`PHASE3.md` §4), with `0029`'s `from_parts`/`parts` — verified by the `so3_*`
   corpus ids under `just conformance` and the envelope's first `helicoid` rows.
3. `tf_tree_math`'s `helicoid` feature scaffold and **Wave 1** behind it — verified by `just
   oracle-tf-tree-math` parity on the `so3_*` ids and `tf_tree`'s own reference twins passing under
   both features.
4. `helicoid` `SEn3`/`SE3` (`PHASE3.md` §5) and **Wave 2** — verified by the `sen3_*` ids and
   `tf_tree`'s bench gate showing no percentile regression.
5. `PHASE4.md` §1.2's dual-quaternion twin and **Wave 3** — verified by the 300 ns / 25 % gate, and
   a re-derived `NS_PER_STEP_ESTIMATE` if the median moves.

## Open questions

1. Does `helicoid` publish (`0.1.0`) at the end of Phase 3, which would let the feature become
   default and the old path be deleted? Waves 1 and 2 do not depend on the answer; decision 2's
   exit does.
2. Do `tf_tree_math::quat_from_rot3(&[f64; 9])` and `helicoid`'s `SO3::from_matrix` agree on their
   reading — the branch choice and `PHASE3.md` §4's `from_matrix_never_iterates`? Wave 1 includes
   that function, so its scope depends on this.

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
