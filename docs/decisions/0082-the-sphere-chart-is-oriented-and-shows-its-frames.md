# 0082: The sphere chart is oriented and shows its frames

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

`S2Chart` freezes the Householder frame `b₁ = He_x`, `b₂ = He_y`, `ν = n + ςe_z`
(`NUMERICS.md` §8, `0063`, `0067`), and keeps it `pub(crate)`. Two needs of locus-tag's sensors
meet there (owner, 2026-10-12):

- **A sensor needs the frame.** A residual on a unit vector (a bearing, a gravity direction) is
  formed in chart coordinates, and its Jacobian chains an ambient `3`-vector derivative through
  the chart. Neither frame is reachable. And the two candidates differ: `δ` is a rotation vector
  `Bδ`, a quarter turn about `n` from the displacement it causes (CH.9(a)), so `∂retract/∂δ` at `0`
  is `[b₁ × n, b₂ × n]`, not `B`. `PHASE6.md`'s `AmbientChart` covers the quaternion and SE(3) only.
- **The frame flips handedness.** `H` is a reflection, so `b₁ × b₂ = ςn` (CH.8(c)), with `ς = -1`
  for `n_z < 0`: in a camera frame, exactly the bearings behind the lens. A consumer's signed
  `det[t₁ t₂ n] ≥ δ` mask then rejects the whole back hemisphere of a lens wider than 180°. It
  would have to take `|det|`, which also accepts a frame that is really mirrored.

S² is unreleased (CHANGELOG `[Unreleased]`), so its coordinates can change now at no consumer's cost.

## Decision

1. **The frame is oriented:** `b₂ = ς He_y`, so `b₁ × b₂ = n` on both hemispheres
   (`docs/maths/charts.md` CH.8(f)). `ς` keeps its role in `ν` (and its reading by comparison,
   `0063` decision 1), so the frame is still continuous on each open hemisphere and discontinuous
   at `n_z = 0`, as any global frame of S² is somewhere. The factor is exact. Every formula of
   CH.9–CH.11 holds with `ς = 1`: `retract` is `cos θ n + (sin θ/θ)(δ₂b₁ - δ₁b₂)`, and
   `retract_jacobian`'s normal part loses its `ς`. `local` and the frame change of CH.10 read the
   frames as computed, so their code is unchanged.
2. **`S2Chart::basis(&self) -> Matrix<S, 3, 2>`** returns `B = [b₁ b₂]`, the frame of `δ`'s
   rotation vector.
3. **`S2Chart::displacement_basis(&self) -> Matrix<S, 3, 2>`** returns
   `T = [b₁ × n, b₂ × n] = [-b₂, b₁]`. That is `∂retract/∂δ` at `δ = 0`, and `Tᵀ = ∂local/∂m` at
   `m = n`. It is orthonormal and oriented, `t₁ × t₂ = n`, and exact: `B`'s columns, one negated.
   Both names are kept, because they are two quantities (R3's "one spelling" is per operation).

## Rationale

- **Rejected: an oriented accessor over an unoriented chart.** A residual's coordinates would then
  live in a frame the chart's `local` and Jacobians do not use.
- **Rejected: telling consumers to use `|det|`.** It cannot tell a mirrored frame from an oriented
  one.
- **Rejected: `T` as a 3×2 `AmbientChart` row.** That trait is Phase 6 and not started; `T` at the
  base is all a residual's chain rule needs now.
- **Return type:** `Matrix<S, 3, 2>` and not a pair of `Vec3`, so that `B * δ` and `Tᵀ * v` are the
  products they are.

## Consequences

- For `n_z < 0`, `δ₂` changes sign; `δ₁` and every `n_z ≥ 0` output are unchanged. `retract(δ)` is, bit for bit, the old `retract(δ₁, -δ₂)` there (the sign moves exactly from
  `g` into `b₂`); 228 of each id's 384 records are regenerated. `just conformance` at both precisions:
  every `s2_*` stratum's max is unchanged or lower except `s2_retract/s2:near-antipode/f64`,
  5.965 → 6.029 `u`, the same arithmetic over resampled points. The
  generator's `s2.basis` multiplies `b₂` by `σ`, its CH.9(a) check drops `σ`, and `s2_retract` and
  `s2_local` are regenerated.
- API §6:
  1. Rules: R1 and R2 are untouched (values, no side); R3 holds, since two names are two
     quantities.
  2. Corpus: no new arithmetic. `basis` returns stored values, and `displacement_basis` is a
     negation of them. Both are covered by the `s2_*` strata and by
     `the_basis_is_orthonormal_and_oriented` and
     `displacement_basis_is_the_derivative_of_the_chart_at_its_base` (against `Dual`, within 4 `ε`).
  3. Domain: `S2Chart::at`'s.
  4. The return is a labelled `Matrix`, not a `bool` or an unlabelled array.
  5. No second spelling.
  6. Ownership: a chart's frame is the chart's.

## Implementation plan

1. One PR: the frame, the two accessors, `NUMERICS.md` §8, CH.8(f), the generator, the regenerated
   `s2_*` corpus and its manifest, the tests. Verified by `just test`, `just corpus-check` and
   `cargo xtask conformance --fn s2_retract|s2_local` at both precisions within `0067`'s bounds.

## Open questions

None.

## Further work

- `T` at a point other than the base, `∂retract/∂δ` at `δ ≠ 0`: `retract_jacobian` read in the
  ambient frame, if a sensor needs it.
