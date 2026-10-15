# 0093: The inverse Jacobian's column is four words

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

`SEn3::jr_inv` evaluates `NUMERICS.md` §5.4 as written: `J_r⁻¹(φ) − ε J_r⁻¹(φ) Q(−ρᵢ, −φ) J_r⁻¹(φ)`.
- **The cost.** Per column: Barfoot's `Q` from seven `3 x 3` hat products, then two `3 x 3`
  products for the sandwich. `jl_inv` is the provided `J_r⁻¹(−τ)`, and `rminus_jacobians` fuses
  the two sides through one set of `Q`'s words.
- **The losses.** On the deep tier, after `0091` and `0092` (both drafts), 17 of the 35 resolved
  domination failures were `sen3_jl_inv_n*` and `sen3_jr_inv_n*`, all near `π` (`θ` from `π − 10⁻²`
  to `π − 10⁻¹²`). In those strata `helicoid` read 1.7–2.5 `u` and GTSAM 1.4–1.7 `u`. On the
  committed tier the worst stratum reads 5.0 `u`, and on the deep tier 6.8.
- **GTSAM's form.** `ExtendedPose3::LogmapDerivative` builds the column as `−J_l⁻¹ H J_r⁻¹`, `H`
  being the derivative of `J_l(φ)ρ`. That is the same sandwich, with fewer products in `H`.

`docs/maths/se3.md` SE.8 remark (ii) already reduces `Q` itself to five words. The inverse reduces
further.

## Decision

1. **`J⁻¹(τ)`'s column is evaluated in four words** (`NUMERICS.md` §5.4, `docs/maths/se3.md`
   SE.17):
   - `J_r⁻¹`: `½ρ^ + S(ρ)`; `J_l⁻¹`: `−½ρ^ + S(ρ)`;
   - `S(ρ) = c(φ^ρ^ + ρ^φ^) + 2c′(φ·ρ)φ^φ^`, with `c` from `NUMERICS.md` §4 and `c′ = dc/dθ²`.

   Each diagonal entry of `S` is `−(2c(φⱼρⱼ + φₖρₖ) + 2c′(φ·ρ)(φⱼ² + φₖ²))`. No `Q` and no product
   of `3 x 3` matrices is formed.
2. **`c′` is the `Dual` lane of `c`'s own evaluation**: `coeffs::jr_inv_coeff_dz`, `0062`'s
   construction for `sinc`. It has no switch, series or group of its own, and `c` is
   `jr_inv_coeff`'s to the bit.
3. **`jl_inv` is overridden** with the same code at the opposite sign of `½ρ^` and the transposed
   diagonal. It is the provided `J_r⁻¹(−τ)` to the bit, and
   `the_fused_inverses_are_the_separate_ones_to_the_bit` compares the two on every block.
4. **Unchanged:**
   - the reference twin, the dense Gauss–Jordan inverse of `jr` (`NUMERICS.md` §14);
   - `jr` and `jl`, which still evaluate `Q` as §5.3 writes it.

## Rationale

**The identity is exact, not fitted.** `−PQP`, with `P = I ∓ ½W + cW²`, maps to itself when every
word is reversed. By SE.7 the reversal-invariant words span `{Ξ, WΞ + ΞW, pW, pW²}`. With the
closed forms of `b, d, e, c`, the four coefficients are `−½`, `c`, `0` and `2c′`.
- **Checked symbolically:** each coefficient identity, with denominators cleared, reduces modulo
  `sin² + cos² − 1` to the zero polynomial (sympy).
- **Checked numerically:** at 110 digits against the sandwich and against `mp.inverse` of the dense
  series (SE.17's *Checked:* line).

**Accuracy.** Per-stratum maxima of the six ids, against `main` at `aacff0e`:

| tier | strata | better | worse | geometric mean of the ratio | worst stratum |
|---|---|---|---|---|---|
| committed | 300 | 300 | 0 | 0.59–0.62 | 4.2–5.0 → 1.36–1.37 `u` |
| deep | 300 | 300 | 0 | 0.57–0.60 | 5.2–6.8 → 1.40–1.46 `u` |

- **Against the oracles:** the six ids are at or below the best oracle on 619 of 624 strata over
  the two tiers. The other five are within 1.17×, and none is resolved.
- **The deep envelope** goes from 35 domination failures to 18. The 17 removed are this family's,
  and none is added.

**Latency.** Consumer loop, 4096 tangents per band, best of 40, one pinned core, `ρᵢ ~ U[−2, 2]`,
`f64`:

| call | `N` | `θ ~ 10⁻⁶..10⁻²` | `0.5..2` | `π − 0.1..π` |
|---|---|---|---|---|
| `jr_inv` | 1 | 52.7 → 11.8 ns | 70.4 → 24.0 | 86.5 → 30.6 |
| `jr_inv` | 3 | 132.2 → 22.6 | 149.1 → 39.6 | 170.4 → 51.2 |
| `jl_inv` | 3 | 148.0 → 23.5 | 165.8 → 40.3 | 184.5 → 54.5 |
| `rminus_jacobians` | 1 | 87.5 → 41.7 | 131.5 → 80.7 | 152.3 → 91.4 |
| `rminus_jacobians` | 3 | 209.7 → 61.2 | 256.5 → 101.2 | 282.2 → 120.5 |

At `f32`, `jr_inv` at `N = 3` goes from 148.9, 157.1 and 175.6 ns to 23.3, 39.7 and 43.3. The
remaining cost is `c`'s branch and the `Dual` evaluation of its series.

**Rejected:**
- **GTSAM's `−J_l⁻¹ H J_r⁻¹`.** It keeps both sandwich products, so it cannot be cheaper than
  dropping them. In the 17 strata GTSAM won, the four words read 1.16–1.39 `u` against its
  1.39–1.66.
- **`Q` in SE.8 remark (ii)'s five words, then the sandwich.** It keeps the two products and their
  rounding.
- **`c′` from the ODE `zc′ = ⅛ − (3/2)c + ½zc²` on every arm.** Below `θ ≈ 1` its right side is `O(z)`
  and cancels. The `Dual` lane takes the derivative of the arm `c` actually runs, a polynomial
  below `θ² = 3.65` (`f64`).

## Consequences

- **Bits.** `cargo xtask determinism` against `main` moves exactly `sen3_jr_inv_n1..3` and
  `sen3_jl_inv_n1..3` at `f64`. No other scored id calls them.
- **`NUMERICS.md` §5.4 now names the evaluation.** The sandwich stays the definition and SE.12's
  statement.
  - **SE.16(b)** bounds the sandwich, and says the program evaluates SE.17 instead.
  - **SE.12's consequence (i)** keeps its cost of `2N` products for the sandwich as written.
- **`SEn3::rminus_jacobians` and `lminus_jacobians`**, the hottest Jacobian rows of a solver on
  SE(3) and SE₂(3), share one `S(ρ)` per column between the two sides.
- **`sen3.rs` loses `Q`'s shared-words plumbing** (`Words`, `q_words`, `q_assemble`), which only
  the fused inverses used. `q_block` evaluates the same words for `jr` in the same order, so `jr`'s
  bits do not move.

## Implementation plan

1. `coeffs::jr_inv_coeff_dz`; `SEn3::{jr_inv, jl_inv}` and the fused `inverses` through
   `InvCols`; `NUMERICS.md` §4 and §5.4; `docs/maths/se3.md` SE.17. Verified by:
   - `just test`, including `sen3_jr_inv_matches_reference_*` and
     `the_fused_inverses_are_the_separate_ones_to_the_bit`, now over the `jl_inv` override too;
   - `just lint`, `just doc`;
   - `cargo xtask determinism` against `main`;
   - `cargo xtask envelope` with and without `--deep`.

## Open questions

None of substance; it is a `draft` until the owner accepts the new evaluation.

## Further work

- **`jr` and `jl`'s `Q` in SE.8 remark (ii)'s five words.** The forward Jacobians still form seven
  hat products per column. Near `π` they read up to 1.7× the best oracle on the committed tier
  (`sen3_jl_n1/theta:pi-1e-3`: 2.31 `u` against Sophus's 1.34), unresolved.
- **`Sim3`'s `log_inverses`.** It inverts the same kind of block and may reduce the same way.
