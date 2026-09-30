# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [Unreleased]

### Added

- Workspace skeleton: `helicoid-linalg` and `helicoid` (empty, `no_std`, `forbid(unsafe_code)`),
  `xtask` stub, lints, `justfile`, CI, `deny.toml`. No group code (`docs/PHASE1.md`).
- `docs/maths/`: non-normative derivations. Notation, labelling and verification convention, and the
  Lie-group foundations (`Ad`, `ad`, the Jacobians of `Exp`, every row of `NUMERICS.md` §2.3, a sign
  audit with exact gaps, the conditioning of `J` and `Ad`). Documentation only; no code, no formula
  change.
- `docs/maths/so3.md`: SO(3) derivations (quaternion action and double cover, `Exp`, `Log` and the
  conditioning of the `acos` form it avoids, `J` and `J⁻¹` solved in the algebra of `{I, W, W²}`,
  `Ad`/`ad`, the action Jacobians, Shepperd's `from_matrix` and the nearest rotation, the Newton
  renormalization and its bounds). Documentation only; no code, no formula change. Five open items
  for the normative documents are recorded in the maths index.
- `docs/maths/coefficients.md`: the coefficient catalogue of `NUMERICS.md` §4 (general terms and radii
  of convergence, the identities between `k, a, b, c, d, e, r` and the `cos(θ/2)` of §3.1, rounding
  error of the exact arms and of Horner, truncation bounds, the crossing that predicts switch-point
  magnitudes, derivatives through `Dual`, the safe argument and the underflow, `π` and `w → 0`
  cases, the reference precision budget, the cost of one shared switch per call-site group).
  Documentation only; no code, no formula change. Nine open items for the normative documents are
  recorded in the maths index.
