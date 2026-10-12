# 0073: Sim(3)'s two sides share one coefficient evaluation

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #141

## Context

[`0070`](./0070-sim3-a-structured-jacobian-and-a-series-on-delta.md) measured `rminus_jacobians`
at 387–468 ns (`f64`), most of Sim(3)'s Jacobian budget. It runs the coefficient group three
times: `Log`'s `v` at `σ`, then `jr_inv` at `-σ` and `jl_inv` at `σ`, each with its own
`sinc_a_b_d`, `e^t` and `expm1 t`, and SO(3)'s `J_l`, `J_l⁻¹` twice. `rplus_jacobians` and
`lplus_jacobians` run it twice at one `σ`: `Exp`'s `V`, then `J`'s eight scalars. `0070`'s
*Further work* proposed `J_l⁻¹ = J_r⁻¹ Ad_{Exp(-τ)}`, which saves about 100 ns but makes the
`NUMERICS.md` §2.3 rows law a measured bound instead of exactly `0`.

The group at `-σ` is mostly the group at `σ`:

- `Δ = σ² + z` is even in `σ`, so are the arm, `form_b` and `σ_w`'s masks, and the safe argument.
- The closed arm's transcendentals are `sin θ/θ, a, b, d` at `z` and `e^t`, `expm1 t` at
  `t = -|σ|`: even, except at `σ = ±0`, where `t` is `σ` itself.
- Each quantity of the series' recurrence is homogeneous in `(σ, θ)`: `σʲ`, `ξ_j`, `ζ_j`, `ξ'_j`,
  `ζ'_j` of degrees `j` to `j - 4`. Negation commutes with rounding, so the recurrence at `-σ`
  holds each one at `σ` times `(-1)ⁿ` **exactly**, `n` its degree; only the sums differ.

## Decision

1. **`coeffs::sim3_coeffs_scaled_pair(σ, z)`** is `sim3_coeffs_scaled` at `σ` and `-σ`, each to the
   bit, from one branch: the series runs its recurrence once into sixteen sums (`series_at`,
   `-σ`'s sum subtracting a term of odd degree), and the closed arm forms its transcendentals once
   (`Transcendentals`) and assembles twice. At `σ = ±0` the mirrored `e^t`, `expm1 t` are evaluated
   again, a branch no other `σ` takes: there the two `t` differ in sign, and on `Dual` in
   derivative. The first version reused them, which is exact on `f64` and `f32` and gave the
   mirrored side the wrong `∂σ` on `Dual`; the test of decision 4 caught it.
2. **`⊖`'s Jacobians** (`sim3::log_inverses`): `Log`'s `v̂` is the pair's `σ` side, the same bits as
   `sim3_v_scaled`; `J_l⁻¹(τ)` takes the `σ` side and `J_r⁻¹(τ) = J_l⁻¹(-τ)` the `-σ` side; SO(3)'s
   `(J_l(φ), J_l⁻¹(φ))` is formed once and transposed for `-φ`, which is exact
   (`so3_tests::jl_is_jr_transposed_to_the_bit`). Three group evaluations become one pair, two
   SO(3) pairs one.
3. **`⊕`'s Jacobians** (`sim3::exp_and_jl`): `Exp`'s `V` is the full group's `v`, the same bits as
   `sim3_v`'s. Two group evaluations become one.
4. **Exact, so tested to the bit.** No formula, switch or output bit moves; the rows law stays `0`.
   `sim3_coeff_tests::the_pair_is_two_calls_to_the_bit` (both arms, `σ = ±0`, across `form_b`,
   `σ_w` and the scaled arm, every recorded ray point, `f64`, `f32` and `Dual<f64, 2>`) and
   `sim3_tests::fused_rows_are_the_composed_rows_to_the_bit` (the four rows against their
   compositions, `|σ|` to 8 and `σ = 0` with a nonzero derivative, `f64`, `f32`, `Dual<f64, 7>`).

`0070`'s `J_l⁻¹ = J_r⁻¹ Ad` is not taken: it trades exactness for a bound, and decision 2 removes
more work while keeping it.

## Measured

The bench gate (`cargo xtask bench-gate --bench groups --only /sim3/`, one pinned core, the
baseline #140's binary with the `rplus_jacobians` row added), candidate / baseline over 6 pairs,
every row within its concurrent floor or below:

| rows | `f64` | `f32` |
|---|---|---|
| `rminus_jacobians` | 0.65–0.72 (270–294 ns, was 387–468) | 0.78–0.84 (258–276 ns) |
| `rplus_jacobians` | 0.78–0.85 (118–145 ns) | 0.73–0.88 (102–125 ns) |
| `jr_inv` | 0.78–0.84 (97–132 ns) | 0.74–0.84 (92–110 ns) |
| `exp`, `log`, `jr`, `jac_mul`, `dense7_mul` | 0.95–1.01, within their floors | 0.95–1.01, within their floors |

`jr_inv` gains without a pair: `jl_blocks_with` is inlined into the inverse's path, which called
`jl_blocks` out of line before. The first version passed SO(3)'s `J_l` into `jl_blocks` as an
argument, and `jr` read 1.02–1.09 in every pair (floors 0.6–2.3%), the `3 x 3` carried through
memory into the out-of-line call; `jl_blocks` forms it again, and `jr` reads 0.95–1.00.

## Consequences

- `rminus_jacobians`, `lminus_jacobians`, `rplus_jacobians` and `lplus_jacobians` are faster and
  their outputs unchanged; no corpus row moves.
- The kernel's closed arm is `closed_from` over `Transcendentals`; the sweep's `closed` is the
  same bits.

## Implementation plan

1. One PR (#141), on top of #140. Verified by `just lint`, `just test`, `just conformance` at both
   precisions (no row moves), and the bench gate above.

## Open questions

None.

## Further work

None.
