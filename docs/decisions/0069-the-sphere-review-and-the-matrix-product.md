# 0069: the S² review, and the matrix product

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #135

## Context

A review of #133 (S², [`0067`](./0067-the-sphere-takes-plane-coordinates-and-both-precisions.md))
found no correctness defect: the chart's formulas, the direction metric and `Vector::dot`'s
bit-identity were checked by hand. It raised seven points.

1. **No test reaches the regimes the Jacobian forms were chosen for.** The proptests draw `θ < √2`
   and uniform points; `local_jacobian`'s `Exp` form exists for `π - α ≤ 10⁻⁴` and
   `retract_jacobian`'s shared `step` for an `m` on the equator. A regression in either would pass.
2. **`retract_jacobian` runs `exp_coeffs` and `jr_coeffs` on the same `θ²`**, two branch trees,
   two roots, three sines; a fused group was proposed.
3. **`local_jacobian`'s `jr_coeffs(α²)`** recomputes what `log_parts` holds (`s = sin α`, `w`);
   closed forms from them were proposed.
4. **`Matrix × Matrix` sums through `sum((0..K).map(..))`**, the form `0067` rejected for `dot`.
5. **`selftest_s2`'s planted subject re-implements `shipped::Sphere::answer`.**
6. **`Vector::dot`'s `first().zip(..)` guard** says less plainly what `N == 0` says.
7. **`Mat2`'s `sandwich` helper is generic over `D`** and correct only because of the
   `D == 2` assertion beside it.

Points 2 and 3 rest on a premise the switches refute. `K`, `COS_HALF`, `A` and `B` all switch to
their series at `θ² = 9.647` (`f64`; `π²` at `f32`, `coeffs/generated.rs`): below `θ ≈ 3.1` both
groups are Horner polynomials, with no root and no sine, and at `θ = 1` — the row the review
priced — neither runs an exact arm. Point 4 was disassembled: in an `#[inline(never)]` wrapper in a
scratch crate, `Mat3 * Mat3` was 215 instructions at `f64` (205 at `f32`), **87 stores and 39
loads to the stack**, for 27 multiplications and 18 additions. The stack traffic is the nested
`array::from_fn` building the output, not `sum`: `Matrix × Vector`, the same `sum` over one
`from_fn`, has none.

## Decision

1. **Edge tests** (`s2_tests::jacobians_hold_near_the_antipode_and_on_the_equator`): both
   Jacobians against `Dual` and the inverse pair at `θ = π - 10⁻ᵏ`, `k = 1..6`, on three bases (one
   at `n_z = -0`), and at three equatorial bases with steps along the equator, where `m_z` is a
   rounded zero. `Dual` is scored over `κ = α/sin α`; the inverse pair over `κ²`, since `α/s`
   inherits `u/(π - α)` from the rounded `s` and `lj · rj` multiplies it by `‖lj‖ = κ`.
2. **No fused group for `retract_jacobian`.** Below `θ ≈ 3.1` there is nothing to share; above,
   fusing `exp_coeffs` with `jr_coeffs` changes which short prefix `k` and `cos(θ/2)` take between
   the two groups' second switches, so `retract_jacobian`'s `m` could part from `retract`'s by an
   ulp — the bit-identity `0067` decision 4 rests on for the frame of `m`.
3. **No closed forms for `local_jacobian`'s `a, b`.** `b = (α - s)/α³` and `a = (1 - w)/α²` cancel,
   and only `helicoid::coeffs` evaluates a cancelling coefficient (`0004`); below the switch they
   are not evaluated at all.
4. **`Matrix × Matrix` is a loop**: `cols` zeroed, each entry the first product then
   `acc + a[k] b[k]`, as `Vector::dot`. The same sums in the same order: bit-identical.
5. **`Sphere::answer_with(chart, record)`** answers both subjects; the sign-bit planting supplies
   its chart, the ratio planting doubles the clean `delta`, which is exact.
6. **`dot`** reads `if N == 0 { return zero } acc = a[0] b[0]; for i in 1..N`.
7. **`Mat2`'s `sandwich`** reads `cov` as a `Mat2` after the assertion and returns `(J Σ) Jᵀ`,
   `Mat3`'s order, the same products and sums.

## Measured

The edge cases (`measure_s2_edges`): `Dual` at most `2.027 κ`, the inverse pair `2.547 κ²`; the
near-antipode `Dual` cases at most `0.102 κ`. Bounds 5 and 6.

`Mat3 * Mat3` after decision 4: 81 instructions at `f64`, 79 at `f32`, no stack traffic.

The bench gate (`cargo xtask bench-gate --bench groups`, one pinned core, the baseline #134's
binary), over every row whose code multiplies two matrices or calls S²'s Jacobians, 70 rows,
concurrent A/A floor median 0.48%:

| rows | candidate / baseline |
|---|---|
| `{so3, se3, se23}/jac_mul` | 0.089–0.113 |
| `{se3, se23}/jr_inv`, `rminus_jacobians` | 0.40–0.60 |
| `f32` `s2/local_with_jacobian` | 0.79–0.90 |
| `f64` `s2/local_with_jacobian`, `se3/gaussian_propagate` | 0.96–0.99 |
| every `so3` `jr_inv`, `rminus_jacobians`, the other Gaussian rows | within the floor |
| `f32` `s2/retract_jacobian` near `π` | 1.023 (past its floor) |

`retract_jacobian` multiplies no two matrices, so its row is layout (`0059`'s finding).

## Consequences

- No output bit moves: `helicoid.csv` and `helicoid--f32.csv` match `main` row for row (the new
  ids aside), and `--self-test` passes.
- The S² edge cases are pinned; `0067`'s *Further work* item on `lj` near the antipode is answered
  against `Dual`, though still by no corpus id.

## Implementation plan

1. One PR (#135), on top of #134. Verified by `just lint`, `just test`, `just conformance` at both
   precisions against `main`, and the bench gate above.

## Open questions

None.
