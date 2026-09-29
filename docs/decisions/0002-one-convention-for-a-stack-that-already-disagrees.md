# 0002: One convention, for a stack that already disagrees

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

The stack already holds two tangent orders and two quaternion layouts:

- **`tf_tree_math`** (published on crates.io): Hamilton, `w` first, active,
  `a * b = T_a_x · T_x_b`, right perturbation `T = T̂·exp(ξ^)`, twist `ξ = [ω, v]`. `tf_tree`'s
  published surfaces expose that order: the C ABI's `tft_plan_at_with_derivatives` (unstable tier),
  Python's `Layout::QuatTwist`, and the `w`-first pose records of the `FORMAT_VERSION = 3` arena.
- **locus-tag** `Pose::retract`: nalgebra, delta `[v, ω]` (`delta[0..3]` translation), and a
  decoupled retraction `(R·Exp(ω), R·v + t)`.
- **locus_fusion**: Sophus, `[υ; ω]`, Eigen `xyzw` storage.
- The literature this project cites (Barfoot, Solà, Sophus, manif) is translation-first; GTSAM and
  OpenVINS are rotation-first.

An earlier design discussion proposed translation-first and `xyzw` storage. Reading `tf_tree_math`
and locus-tag's `pose.rs` before adoption reversed both: the only published public surface in the stack
is rotation-first and `w`-first.

## Decision

1. **Quaternions:** Hamilton, stored `w` first (`Quat { w, x, y, z }`), active rotations.
2. **Composition:** `a * b` is `T_a_x · T_x_b`; `X * p` maps a point from the child frame into the
   parent frame.
3. **Tangents are rotation-first:** `[φ; ρ₁; …; ρ_N]` for SE_N(3), `[ω; v]` for a twist, `[θ; ρ]`
   for SE(2), `[φ; ρ; σ]` for Sim(3). Named fields everywhere; the order is visible only through
   `write_dense`/`read_dense`.
4. **Right perturbation is the default; both sides are first-class and always spelled**
   (`rplus`, `lplus`, `rminus`, `lminus`, `Side`). No `Add`/`Sub` for ⊕/⊖.
5. **`Log`'s canonical branch** is $\theta \in [0, \pi]$ via the `w ≥ 0` flip; at `w = +0` it is a
   function of the quaternion's sign (`NUMERICS.md` §3.2).
6. **Converters for every other order**, named after it (API R3): `Quat::from_xyzw`, `to_xyzw`,
   `from_jpl` (Hamilton $(w, -x, -y, -z)$ from JPL $(x, y, z, w)$, same matrix — Sommer et al.
   2018), `Twist::from_translation_first`, `to_translation_first`.
7. `NUMERICS.md` states every formula already permuted into these conventions.

## Rationale

- **The published surface decides.** `tf_tree`'s ABI, Python layouts and arena format cannot change
  order without a breaking change on published artifacts (an arena change costs a `FORMAT_VERSION`);
  locus-tag's `[v, ω]` is internal to its LM. Matching the published one moves the permutation to
  the one place where it is free.
- **Order is cheap to convert and expensive to get wrong.** Named fields make the dense order a
  one-function concern; the failure mode (a covariance permuted silently across a library boundary)
  is prevented by the type, not by care.
- **Rotation-first costs nothing structurally.** SE_N(3) Jacobians are lower instead of upper
  block-triangular; the dual-matrix algebra is identical (0005).
- **Alternatives lost:** translation-first (Solà/Sophus/manif) would have matched the literature
  and two oracles, at the price of a permutation at `tf_tree`'s published boundary forever; `xyzw`
  would have matched Eigen, at the price of `tf_tree`'s arena and C ABI (whose C++ wrapper already
  transposes at the Eigen boundary).

## Consequences

- Every oracle runner converts inside the runner, in a named function with a hand-case unit test
  (`PHASE1.md` §7).
- Phase 1's seeded defect "translation-first `Exp`" must be caught by every `rho:*` stratum.
- locus-tag's migration converts at its LM boundary (`PHASE5.md` §6).

## Implementation plan

1. Conventions paragraph in `crates/helicoid/src/lib.rs` rustdoc, verbatim from this record —
   verified by `just doc`.
2. The converters of item 6 with hand-case tests — verified by `quat_converters_round_trip`,
   `twist_translation_first_round_trip`, `from_jpl_matches_matrix`.
3. Seeded defect — verified by `just conformance --self-test`.

## Open questions

None.
