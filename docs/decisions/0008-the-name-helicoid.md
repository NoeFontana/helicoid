# 0008: The name `helicoid`

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

The library needs a crate name, a family rule for sibling crates, and a line against the sensor
models that will sit beside it. Availability was checked on 2026-09-25 (crates.io and PyPI).

## Decision

1. **`helicoid`**: the surface swept by a screw motion — the geometry of SE(3)'s one-parameter
   subgroups. Free on crates.io and PyPI at the check (an unpublished Helix-editor WIP uses the
   word; no registry conflict).
2. **Crates:** `helicoid-linalg` (leaf), `helicoid` (groups). Reserved by convention, unpublished:
   `helicoid-spline` (0011), `helicoid-py` (ledger).
3. **The family rule:** `helicoid*` contains nothing that knows a sensor exists. Sensor models are
   `locus-camera`, `locus-imu`, `locus-lidar` (the locus family; `locus-lidar` deferred, no
   consumer). The solver crate is separate and not named here.
4. **Type names** follow the mathematics: `SO2`, `SO3`, `SE2`, `SE3`, `SE23`, `SEn3`, `Sim3`, `S2`;
   clippy's `upper_case_acronyms` is allowed workspace-wide.

## Rationale

Rejected at the check: `cartan` (an active Riemannian-geometry crate), `lieu` (taken), a
`sensor-*` prefix (unowned and generic; `sensor-scd30` is an unrelated driver), `cameras` (a capture
library). A name tied to a stack product would misstate the scope; `helicoid` names what the crate
computes.

## Consequences

- Publishing reserves `helicoid` and `helicoid-linalg` at `0.0.x` when Phase 3 lands.

## Implementation plan

1. Workspace and crate names — verified by `cargo publish --dry-run` in CI from Phase 3.

## Open questions

None.
