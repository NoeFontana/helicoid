# 0011: Continuous time waits for a consumer

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

Three consumers may want Lie-group splines: locus-calib (camera–IMU calibration), locus_fusion
(continuous-time or rolling-shutter estimation), and `tf_tree` Phase 6 (cumulative B-spline
interpolation, kept by `tf_tree` 0009 §4). None has stated its problem size, knot spacing, order,
derivative orders or precision. The design choices that matter — uniform vs non-uniform knots,
SE(3) vs split, Kalibr parity — depend on exactly those answers. Kalibr's rotation is a Euclidean
spline on a rotation-vector parameterization, so a locus-calib comparison against Kalibr confounds
the parameterization with everything else unless the protocol accounts for it.

## Decision

1. `helicoid-spline` is **gated**, not scheduled. `PHASE7.md` is a requirements artifact; its §0.0
   lists the gate.
2. The gate opens on a consumer record answering `PHASE7.md` §0.0's rows.
3. The core owes splines only what it owes everyone: `Ad`, `ad`, `Jr`, `Jr⁻¹`, `Jl`, `Jl⁻¹` on
   every group and product (Phase 3), geodesics (Phase 4), `Dual` (Phase 2). **No spline-shaped hook
   is added to the core in advance.**

## Rationale

A spline designed without a consumer encodes guesses as API; `tf_tree` D21 applied the same
reasoning to its `tf2` shim. The primitives the Sommer et al. recursions need are general-purpose
and ship anyway, so gating costs nothing on the critical path.

## Consequences

- `PROJECT.md` §5.1 carries the row; `PHASE7.md` existing is not permission to build it.

## Implementation plan

None until the gate opens.

## Open questions

None for this record; `PHASE7.md` §3 holds the spline's own questions.
