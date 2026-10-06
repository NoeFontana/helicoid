# 0049: the boundary is what removes a way to be wrong, and the consumer comes before the rest of Phase 4

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** step 1 lands with this record; steps 2–7 are the integration

## Context

[`0009`](./0009-what-helicoid-does-not-own.md)'s ownership test asks whether an operation needs
knowledge of a sensor, a solver, or a consumer's state or storage. It is the right test for what
`helicoid` may **not** own, and it is silent on what `helicoid` should **decline** to own although
it may. `Quat::dot` passes `0009` cleanly. So does `Quat::normalize`. That silence is what
[`0044`](./0044-four-primitives-the-first-consumer-names-and-no-spec-does.md) walked into — it
shipped four items because [`0041`](./0041-the-integration-is-an-adapter-at-tf-tree-math.md)'s wave
lists name them — and what
[`0048`](./0048-the-relative-transform-pair-earns-the-surface-dot-and-norm-do-not.md) undid for two
of them one PR later.

Three instances now, read off the consumers rather than off a spec. The verdicts do not agree, and
that is the point:

1. **`Quat::{dot, norm}` — removed** (`0048`). No caller in either library crate, nothing private
   reached (`w, x, y, z` are `pub`, `norm_sq` is public), and `dot`'s stated purpose was a
   worse-conditioned spelling of a shipped routine.
2. **`Quat::normalize` and `Iso3::normalized` — not a delegation at all, and the `0048` reasoning
   does not reach them.** A projection onto the unit sphere is as group-ish an operation as exists;
   `0048`'s "not a group operation" argument says nothing here. What disqualifies them is that
   **they are a different operation from the one `helicoid` ships**. `tf_tree_math::Quat::normalize`
   is the exact projection `q/‖q‖`, with a `1e-300` guard returning a near-zero quaternion
   unchanged rather than producing infinities. `helicoid::Quat::renormalize` is **one Newton step**,
   a drift repair whose accuracy statement holds only while `|‖q‖² − 1|` is at most `2^-26.29`
   (`f64`) or `2^-11.79` (`f32`). The two agree on the band a composition chain drifts into and
   nowhere else, and the names hide it: delegating `normalize` to `renormalize` is a silent
   behaviour change for every input outside that band, and `Iso3::normalized` — whose whole body is
   `self.q.normalize()` — inherits it. Nothing in `tf_tree_math`'s signature restricts its argument
   to the band.
3. **locus-tag's `Pose::retract` — the opposite verdict, on an item that looks more local than
   either of the above.** It is `pub(crate)`, it has one call-site family (every pose LM plus
   `board.rs`), its delta is translation-first `[ρ, ω]` where `0002` is rotation-first, and it is
   four lines over public `nalgebra` fields. Every surface indicator says *leave it there*. It must
   still come from `helicoid`, because it is a **retraction**, and
   [`0012`](./0012-a-retraction-is-a-chart.md) is NORMATIVE that a retraction is a `Chart` here and
   never an inline formula in a consumer. Its own rustdoc is already drifting in the way `0012`
   predicts: it describes itself as holding "to first order (`V(ω) ≈ I`)", as though it were an
   approximation of the screw retraction, when it is `PHASE5.md` §1.3's `Decoupled` chart — an
   **exact** retraction in its own chart, which is why the LM iterates it produces are exact and
   not first-order.

So neither "could the consumer write it" nor "is it small" nor "is it a group operation" separates
the three. The project has decided this question three times by review and has no rule to cite.

## Decision

1. **The boundary test, NORMATIVE, two clauses.**

   `helicoid` **owns** an operation when a consumer writing it itself could get a **convention or a
   numerical choice** wrong in a way no test of that consumer would catch: a tangent order, a side,
   a chart, a switch point, a series, a branch cut, a quaternion sign, a Jacobian's frame. Those are
   D1, `0002`, `0004`, `0005` and `0012`'s subject matter, and a consumer's own test suite cannot
   see them, because it compares the consumer against itself.

   `helicoid` **declines** an operation whose only failure mode is arithmetic over its own public
   fields, and whose consumer-side version is already exact or already has a domain the consumer's
   storage defines: `dot`, `norm`, an exact `normalize`, `to_bits`/`from_bits`, a `Vec3` sum.

   **The question is not whether `helicoid` may own it (`0009`) and not whether the consumer could
   write it. It is whether `helicoid`'s version removes a way to be wrong.** A relocation that
   removes no failure mode is not free: it adds a conversion at the boundary, an edge in the
   consumer's dependency graph, a second spelling against `API.md` R3, and — as both of `0044`'s
   removed items showed — a `# Domain` section R6 has no stratum to enforce.

2. **`API.md` §6 gains question 8**: *"Does it remove a way to be wrong that the consumer cannot
   see, or does it only relocate code the consumer already has right?"* Question 6 stays as it is;
   it asks a different thing.

3. **`0041`'s Wave 1 and Wave 2 lists are corrected.** `Quat::normalize` and `Iso3::normalized`
   keep `tf_tree_math`'s own bodies, by Context item 2 — a different operation, not a
   non-group one. **`helicoid` adds no exact `normalize`**: `API.md` R3 is one spelling per
   operation, and `renormalize`'s entire value is that it has neither a `sqrt` nor a division, which
   is what `Iso3::normalized`'s own rustdoc says it wants kept off the hot path. A consumer that
   needs the exact projection writes it over `w, x, y, z` and `norm_sq`.

   `norm_sq` stays public, and it is the exception decision 1 predicts rather than one against it:
   it is how a consumer evaluates `renormalize`'s **stated** domain, `|‖q‖² − 1|` at most
   `2^-26.29`, and a consumer that cannot read that quantity invents a tolerance over `norm`
   instead. Clause two is about operations a consumer would otherwise get right, and a domain check
   against the wrong quantity is not one. `0048`'s argument for removing `dot` and `norm` rests on
   `norm_sq` being public, so the two readings are one reading.

4. **`quat_from_rot3` delegates, and that *changes its output*.** Its rustdoc says "the result is
   not normalized"; `SO3::from_matrix` normalizes, written out rather than taken from
   `from_wxyz_normalized` so a garbage matrix stays garbage instead of panicking. So this row is a
   **recorded behaviour change with a parity number**, from the `so3_from_matrix` corpus id the
   oracle now answers (`0041`'s scheduled reading, landed with #92) — not a port. The un-normalized
   output has no caller contract to break: `tf_tree_c::layout::read` rejects `|det R − 1| > 1e-6`
   upstream of it, so every reachable argument is already a rotation to that tolerance.

5. **`0041`'s implementation plan is re-ordered: its decision 7 — publish `0.0.1` — moves ahead of
   the remainder of plan step 4, whose `geodesic_jacobians` item becomes step 8.** Every Wave 1 and Wave 2 item maps onto surface that has already shipped:

   | `tf_tree_math` | `helicoid`, shipped today |
   |---|---|
   | `exp_so3(Vec3) -> Quat` | `SO3::exp(&τ).quat()` |
   | `log_so3(Quat) -> Vec3` | `SO3::from_quat_unchecked(q).log()` |
   | `quat_from_rot3(&[f64; 9])` | `SO3::from_matrix(&Mat3)` — decision 4 |
   | `slerp(qa, qb, s)` | `SO3::geodesic`, the provided body (GE.14, `0045`) |
   | `Quat::rotate(Vec3)` | `SO3::act` |
   | `Quat::{dot, norm, norm_squared}` | **stays** `tf_tree_math`'s (`0048`) |
   | `Quat::normalize` | **stays** `tf_tree_math`'s (decision 3) |
   | `exp_se3([f64; 6])` | `SE3::exp` — `[ω, v]` on both sides, no conversion (`0002`) |
   | `log_se3(Iso3)` | `SE3::log` |
   | `Iso3::inverse` | `SE3::inverse` |
   | `Iso3::mul_inv` | `SEn3::mul_inv`, every `N` (`0048`) |
   | `Iso3 * Iso3` | `SE3`'s `Mul` |
   | `Iso3::adjoint(&Twist)` | `SE3::adjoint()` then `Jac::apply` — no dense 6×6 (`0005`) |
   | `Iso3::adjoint_inv(&Twist)` | the same through `Jac::inverse` |
   | `Twist::to_spatial(&Iso3)` | a composition over the row above, no new maths |
   | `Iso3::normalized` | **stays** `tf_tree_math`'s (decision 3) |

   What is left of Phase 4 gates only Wave 3: §1.2's screw twin is the dual-quaternion path
   `sclerp`/`screw_pow` need, and
   [`0043`](./0043-the-geodesic-jacobian-ships-the-cancellation-free-form.md)'s `geodesic_jacobians`
   gates **no `tf_tree` wave at all** — no function in `tf_tree_math` takes a geodesic Jacobian. The
   order becomes **publish → Wave 1 → Wave 2 → screw twin → Wave 3**, with `geodesic_jacobians` and
   the first blessed envelope moved onto the Phase 5 / locus-tag path they actually serve.

   The reason is not only schedule. Wave 1 is the first time a crate `helicoid` does not own
   compiles against it. `0041` decision 3 already makes the feature the **parity** instrument; it is
   equally the **surface** instrument, and this record is the evidence: reading two consumers found
   three wrong rows in `0041`'s own wave lists and one behaviour change nobody had written down.
   A surface defect found after the screw twin is written is found at the wrong time, and D7 —
   instrument before algorithm — is the same argument one layer up.

6. **locus-tag's gate is named, so it stops being a vague second.** Its entire *production* Lie
   surface is five items: `Pose::retract` (Context item 3), `Pose::inverse`, `Pose::adjoint`,
   `covariance_{camera_to_body, body_to_camera}`, and `quat_from_so3` — whose only non-test callers
   are two in `locus-py`. Everything else matching `retract`/`Rotation3`/`adjoint` under
   `crates/locus-core/src/` is a test. Two of the five need `helicoid` items that **do not exist**:

   - **`Decoupled`** (`PHASE5.md` §1.3), which is `Pose::retract` exactly. Blocked behind
     `PHASE5.md` §1's `Chart`/`Manifold`/`WithChart`.
   - **The translation-first covariance permutation** `Π Σ Πᵗ` (`docs/maths/se3.md` SE.14(d),
     `docs/maths/gamma-gaussian.md` GG.13, `0041` *Further work* 1). locus-tag publishes a `[ρ, ω]`
     6×6 covariance; `helicoid` is rotation-first, and `mint` stops at 4×4 so it covers neither this
     6×6 nor locus-tag's 8×8 blocks.

   And a third constraint no document names: **`Pose { pub rotation: Matrix3<f64>, pub translation:
   Vector3<f64> }` is matrix-backed and public**, where `helicoid` is quaternion-backed. Either the
   retract converts matrix ↔ quaternion per LM iteration — up to twenty per pose, each paying
   `from_matrix`'s pivot selection and its rounding — or `Pose`'s representation changes, which is a
   locus-tag published-API change. `PHASE5.md` §6 says "nalgebra types convert through `mint`",
   which answers *how* and not *how many per iteration*. That is measured on locus-tag's paired
   bench instrument before the migration is designed, not argued. (`PHASE5.md` §6 names
   `scripts/bench_gate.py`; the instrument is `tools/bench/` — `compare/`, `metrics.py`,
   `strata.py` — driven through `tools/cli.py bench`. The cited path does not exist and §6 is
   corrected with this record.)

## Rationale

`0044`'s instinct — ship what the first consumer names — is the one this record replaces, and its
own *Rationale* had already conceded the problem ("for two of the four that is right") without a
rule that said which two. The rule had to separate three cases that every surface heuristic gets
wrong, and "does it remove a way to be wrong" is the only formulation that does: it keeps
`Pose::retract` (a convention a consumer's tests cannot see), releases `Quat::dot` (arithmetic over
public fields), and releases `Quat::normalize` for a reason neither of the others supplies (the
operation is not the one `helicoid` has).

The inverse rule — *the consumer writes anything it could write* — is what `0012` forbids for
retractions and D1 for group operations, and it is precisely how a tangent order drifts: locus-tag's
translation-first delta and its `V(ω) ≈ I` rustdoc are both live examples, in a repository whose
tests all pass.

**Finishing Phase 4 first**, which is the order the approved plan carries, loses on two counts. It
serialises the integration behind `geodesic_jacobians`, which no `tf_tree` wave consumes, and it
designs the screw twin and the Jacobian against a specification instead of against a caller. The
counter-argument is that publishing `0.0.1` with a surface still moving is a commitment — and it is
not, because `CLAUDE.md`'s Status already says every `0.0.x` may break every other, which is the
whole reason `0041` decision 7 chose `0.0.1` over `0.1.0`.

## Consequences

- `tf_tree_math` keeps `dot`, `norm`, `norm_squared`, `normalize`, `Iso3::normalized`,
  `to_bits`/`from_bits`, `Vec3`'s arithmetic and `LerpSlerp` (`0009`) **permanently**, not until a
  later wave. Its `helicoid` feature is a partial delegation by design even after Wave 3, which
  `0041` decision 6's per-function fallback already provides for.
- `tf_tree_math::reference.rs` stays, with a changed job — `0041` *Further work* 4, answered.
  `helicoid`'s `NUMERICS.md` §14 twins are the **numerics** twins (D6); `tf_tree_math`'s become the
  **adapter** twins, which answer "did delegation change the answer". Two different questions, both
  worth a file.
- `helicoid` gains no exact `normalize`, and a later request for one is answered by this record
  rather than re-litigated.
- The publish decision is now the single critical-path item for the whole integration, and nothing
  in Phase 4 blocks it.
- One row of the parity table is a recorded behaviour change rather than a parity measurement
  (decision 4), and `PHASE4.md` §5.1 has to say so, or the table claims agreement it is not
  measuring.

## Implementation plan

1. This record; `0041`'s amendment banner and its corrected wave lists; `API.md` §6 question 8;
   `PHASE4.md` §5.1's note on decision 4; `PHASE5.md` §6's corrected instrument path — verified by
   `just lint` and the row in `decisions/README.md`.
2. Release prep: `0.0.1`, the `readme` keys and per-crate READMEs both manifests lack, then
   `cargo package --list` and `cargo publish --dry-run` for both crates — verified by those two
   commands. **`cargo publish` itself waits for the owner's word.**
3. Publish, then Wave 1 behind `tf_tree_math`'s `helicoid` feature — verified by
   `just oracle-tf-tree-math` parity on the `so3_*` ids and `tf_tree`'s own reference twins passing
   under both feature arms.
4. Wave 2 — verified by the `sen3_*` ids and `tf_tree`'s bench gate showing no percentile
   regression.
5. `PHASE4.md` §1.2's screw twin, then Wave 3 — verified by the 300 ns / 25 % gate and a re-derived
   `NS_PER_STEP_ESTIMATE` if the median moves.
6. `0043`'s `geodesic_jacobians` and the first blessed envelope, off the `tf_tree` path — verified
   by `just test` and `just envelope --bless`.
7. locus-tag: `PHASE5.md` §1's charts, `Decoupled`, the permutation, then its migration — verified
   by its paired bench instrument with rotation p95 included, per `PHASE5.md` §6.

## Open questions

None.

## Further work

1. **The two `so3_geodesic` domination failures are not the screw twin's to fix, and the plan reads
   as though they are.** Measured over the committed corpus, `helicoid` **wins** `geo:consecutive`
   — 1.572 `u` against `tf_tree_math::slerp`'s 2.187 — and **loses** `geo:generic` (2.721 against
   1.834) and `geo:near-pi` (2.429 against 1.642). §1.2's twin addresses `se3_geodesic`'s
   translation path, which `0048` already took off the scored maximum, and leaves the shared
   rotation path untouched. Note what the measurement does to GE.14's argument: it predicts
   `helicoid` should win **near π**, where slerp's `acos` has infinite slope, and the corpus says
   the oracle wins there by 1.48×. So either the conditioning gain is real and the
   `Log` → `scale` → `Exp` round trip's own roundings exceed it, or the argument is incomplete — and
   the split by stratum says which, because the round trip is cheapest exactly where `helicoid`
   wins. The candidate is a direct-blend arm above a swept switch: slerp's weights in the `atan2`
   spelling, no round trip, the series/exact choice made at a generated switch (`0004`). It is also
   the one `helicoid` optimization that lands directly on `tf_tree`'s benched hot path. Measure
   first (D7); this is not an authorisation to write it.
2. The permutation's shape — a `Jac` method, a free function, or a `Gaussian` constructor — and
   whether locus-tag's 8×8 blocks take the same route (`0041` *Further work* 2).
3. Whether `Pose` becomes quaternion-backed, which decision 6 measures rather than decides.
4. `0028` chose option A, "narrow `SEn3Jac` to `pub(crate)`", and `SEn3Jac` is `pub` and re-exported
   from `lib.rs`. The adapter needs it nameable, so the code is right and the record is stale;
   whether `0028` is amended or superseded is not this record's call.
