# 0009: What `helicoid` does not own

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

A foundation library grows by accretion: every consumer has one adjacent thing it would like the
shared crate to hold. Each such addition widens the dependency, review and stability surface of
every other consumer.

## Decision

**The ownership test:** a feature belongs in `helicoid` iff it needs no knowledge of a sensor, a
solver, or a consumer's state or storage layout. By that test, `helicoid` does **not** own:

| Not owned | Owner |
|---|---|
| Camera, IMU, LiDAR models (projection, distortion, noise, intrinsics) | `locus-camera`, `locus-imu`, `locus-lidar` |
| Solvers: dense LM (omnisac's LM core), batch sparse NLLS, Schur, robust losses, covariance recovery | omnisac (tiny dense LM); the separate solver crate (batch sparse) |
| IMU preintegration, marginalization, FEJ, sliding windows | locus_fusion (Γ₁/Γ₂ are `helicoid`'s primitives) |
| Storage formats, `Pod` records, arenas, serialization | each consumer (D2) |
| `tf2`-compatibility behaviour (`LerpSlerp`'s fallback, `tf2` quirks) | `tf_tree` |
| Hand-eye, AX = XB, calibration pipelines | locus-calib |
| Point-cloud deskew | `tf_tree` cut it (its D8); `locus-lidar` if ever |
| RNG, sampling from `Gaussian` | consumers |

`helicoid` owns: groups, charts, sides, Jacobians, coefficients, geodesics, Γ functions,
`Gaussian` as a mean/covariance/side triple with its propagation, and the fixed-size linear algebra
those need (including `eig3`/`svd3`/`solve_cubic`, which every consumer above needs identically).

## Rationale

Every row fails the test by needing something only its owner has: sensor physics, a problem
structure, a state layout, `tf2`'s observed behaviour. The consequence of getting one wrong is the
same each time: a consumer's release cadence becomes `helicoid`'s.

## Consequences

- A request for a listed item is refused by citing this record; a request not listed goes through
  the test in a new record.
- The solver crate depends on `helicoid`, never the reverse.

## Implementation plan

None: this record constrains; it builds nothing.

## Open questions

None.
