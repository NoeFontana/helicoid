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
- `docs/maths/se3.md`: SE_N(3) derivations (the matrix group and the `x_1 = v`, `x_2 = p` convention of
  SE₂(3), `Exp`/`Log` from the hat matrix's series, `Ad`/`ad`, Barfoot's `Q` block derived from
  `J_l = Σ adⁿ/(n+1)!` by the block-triangular powers of `ad` and `W³ = −θ²W` and identified with `b, d, e`,
  `J_r` and two block relations, the dual-matrix algebra with proofs of closure, product (27 + 54N
  multiplications), inverse and `apply_transpose`, `J⁻¹` without a second closed form, the SE(3) action
  Jacobians, the rotation-first/translation-first permutation of tangents, Jacobians and covariances,
  what the exact arms of `b, d, e` cost `Exp` and `Q` (and the error of `Q` at a switch), and the rounding of
  `Log`'s translation part and of the blocks of `J⁻¹`). Documentation only; no code, no formula
  change. One open item for the normative documents is extended and one added in the maths index.
- `docs/maths/so2-se2.md`: SO(2) and SE(2) derivations (unit-complex `SO2`, `atan2` sensitivity; `Exp` with
  `V = αI + βK = s(θ)R(θ/2)`, `Log` and the three forms of `V⁻¹`, the two preimages at `|θ| = π`; `Ad`,
  `ad` and `ad³ = −θ²ad`; the rotation-first `J_r`, `J_l` and their inverses derived as `I ± a·ad + b·ad²`
  and `I ∓ ½ad + c·ad²`, with the dense entries written out and marked **Proposed** for `NUMERICS.md` §6
  and §2.4; the reduction of `Q` on the plane and agreement with Solà et al. App. C through the
  translation-first permutation; the scalar functions with series and singularities; rounding, the
  conditioning of `J` and a sign audit). Documentation only; no code, no formula change. Four open items
  for the normative documents are added in the maths index.
