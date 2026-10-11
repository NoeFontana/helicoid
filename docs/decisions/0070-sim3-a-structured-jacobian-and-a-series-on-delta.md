# 0070: Sim(3) — a structured Jacobian, and a series on Δ

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #136 (this record and the spec edits); the plan below

## Context

`NUMERICS.md` §9 owes Sim(3)'s `Exp`, `Log`, `J_r`, their coefficients as functions of `(σ, θ)`
and the two-dimensional series at the joint limit `σ² + θ² → 0`, "as an edit to this section
plus a record", before `PHASE5.md` §3 starts. `docs/maths/sim3.md` derives every piece (SM.1–SM.13,
checked to 60–150 digits) and leaves four choices and four corrections to this record
(`docs/maths/index.md`, open items).

**The choices.**

1. **The Jacobian type.** `PHASE5.md` §3 and [`0005`](./0005-the-jacobian-is-a-dual-matrix.md)
   say dense `7 × 7`, because the scale column breaks `SEn3Jac`'s dual structure. SM.9 shows every
   Sim(3) Jacobian — `Ad`, `ad`, `J_l`, `J_r`, their inverses and the §2.3 rows built from them —
   has the pattern `[A 0 0; B C f; 0 0 e]` (`A`, `B`, `C` `3 × 3`, `f` a column, `e` a scalar),
   closed under product (`117` multiplications against `343` dense) and inverse.
2. **The joint limit.** The closed forms give `Exp`'s value and `Log` to `O(u)` everywhere
   (SM.7(b)), but the Jacobian blocks `Q` and `Γ₂(M)ρ` lose `≈ 4.5u/r_Δ`, `r_Δ = √Δ`,
   `Δ = σ² + θ²`: `2e12 u` at `r_Δ = 10⁻¹²` (SM.7(c)). Two arms remove it (SM.7(d)): the series
   of SM.6(b) in total degree, one recurrence for all eight scalars behind one mask on `Δ`
   (`≤ 5.5u` at `K = 16`), or the exact split of SM.6(d) into functions of `σ` (`≤ 6.1u`, more
   masks).
3. **Large `σ`.** The direct closed forms overflow silently long before `e^σ` does: `Log`'s `i₂`
   from `σ ≈ 242` (`f64`; `33` at `f32`), `v₁`, `v₂` from `703` (`84`). The scaled arrangement of
   SM.7(g) holds to `ln MAX`.
4. **What a sweep minimizes** (`0004` item 1): a coefficient's own error (no digit of `v₂` below
   `r_Δ ≈ 10⁻⁸` in closed form) or the consumer's, the blocks', where each scalar is multiplied by
   a power of `θ` that suppresses its error (SM.7(b)–(c), SM.11(e)).

**The corrections.** `NUMERICS.md` §12's row "`Sim3` | `σ` finite" promises no finite result
beyond `ln MAX` and has no row for `log`, `jr_inv`, `jl_inv`. `PHASE5.md` §3 takes `sim3_log`'s
reference from `mp.logm`, complex and non-principal from `θ ≈ 3.03` at `σ = 0` and earlier as
`|σ|` grows (`2.19` at `σ = 6`, SM.8). Its ids omit `jl`, `jl_inv` and the action. `Real` lacked
`exp` and `exp_m1`; [`0068`](./0068-real-owes-exp-and-exp-m1.md) adds them.

## Decision

1. **Types** (`helicoid::sim3`, re-exported):
   - `Sim3<S>`, stored `(q, t, σ)` with `s = e^σ` (§9), fields private; `from_parts(SO3<S>, Vec3<S>,
     S)`, `parts`, `rotation`, `translation`, `log_scale` (`σ`), `scale` (`e^σ`), `renormalize`
     (the quaternion's step), `mul_inv`, `inv_mul`, `Mul<Point3<S>>` (`sRp + t`), `act_many` and
     `act_jacobians::<Sd>` (SM.12), as `SEn3`'s.
   - `Sim3Tangent<S> { phi, rho, sigma }`, dense order `[φ; ρ; σ]`, `DOF = 7`.
   - **`Sim3Jac<S>`**, the pattern `[A 0 0; B C f; 0 0 e]`, implementing `Jac<S, Sim3Tangent<S>>`:
     `mul` is `(AA', BA' + CB', CC', Cf' + f e', e e')`, `inverse` is
     `(A⁻¹, -C⁻¹BA⁻¹, C⁻¹, -C⁻¹f/e, 1/e)` by `Mat3::inverse_adj`, `sandwich` skips the structural
     zeros. Dense twins `reference::sim3jac_{mul, inverse, apply, apply_transpose, sandwich}`.
     **`PHASE5.md` §3's "dense `Mat<7>`" and `0005`'s note on Sim(3) are amended by this record.**
   - `Sim3` is a `LieGroup` and a `Manifold` (`RightChart`).
2. **Formulas, adopted into `NUMERICS.md` §9** as `docs/maths/sim3.md` states them:
   - the group law, inverse and action (SM.2);
   - `Exp τ = (Exp φ, Vρ, σ)`, `V = v₀I + v₁W + v₂W²` (SM.4), the closed forms of SM.5(b);
   - `Log X = (φ, V⁻¹t, σ)`, `φ` through the quaternion's `atan2` (§3.2), `V⁻¹ = i₀I + i₁W + i₂W²`
     (SM.8);
   - `Ad`, `ad` (SM.9);
   - `J_l = [J_l(φ) 0 0; Q V −Γ₂(M)ρ; 0 0 1]`, `J_r(τ) = J_l(−τ)`,
     `Q = [y]_× J_l(φ) + ∂_φ y` with `y = Vρ` (SM.10(b)), `Γ₂(M) = w₀I + w₁W + w₂W²` (SM.11(b));
     the inverses by blocks (SM.10(c)), not a second closed form.
3. **The coefficients** `v₀, v₁, v₂, w₀, w₁, w₂, ∂_z v₁, ∂_z v₂` are one group in `helicoid::coeffs`
   (`sim3_coeffs(σ, z)`, with a value-only `sim3_v(σ, z)` for `Exp` and `Log`), three arms:
   - **`Δ = 0`**: the limits `1, ½, ⅙, ½, ⅙, 1/24, −1/24, −1/120` (SM.7(h)), selected.
   - **`Δ < r_×²`: the series of SM.6(b)** — one recurrence (SM.3, SM.11(c)) at **per-sequence
     depths**, `ξ_k` to `K`, `ζ_k` to `K + 1`, `ξ'_k` to `K + 2`, `ζ'_k` to `K + 3`, so every scalar
     is truncated at total degree `< K` and the tail is a function of `r_Δ` alone (SM.7(d)(ii)):
     one mask, on `Δ`, for every direction. `v₀` and `w₀` by their series in `σ` to the same degree.
     The constants `1/(k + ℓ)!` are generated, correctly rounded per precision.
   - **`Δ ≥ r_×²`: the closed forms** of SM.5(b), SM.11(b)–(c), with `a, b, d` from the kernel's
     arms; `v₁` by form b for `Δ ≥ Δ_b` and form a below; `w₀ = g₂(σ)` by its series in `σ` below
     `|σ| < σ_w` (its exact arm loses `2u/|σ|`); `v₀ = 1`, `w₀ = ½` at `σ = 0`. For `σ > 0` the
     **scaled arrangement** of SM.7(g): `e^σ → 1`, `expm1 σ → −expm1(−σ)`, each scalar formed as
     `e^{-σ}c` and multiplied by `e^σ` last; `Log` applies `V⁻¹ = e^{-σ}V̂⁻¹` to `t` as
     `e^{-σ/2}(e^{-σ/2}t)`.
4. **Generated, not typed** (`0004`): `K` and `r_×²`, `Δ_b`, `σ_w` and `g₂`'s series, per
   precision, written by `cargo xtask thresholds` from a sweep over rays `σ = r_Δ cos λ`,
   `θ = r_Δ sin λ` (`λ` over `[0, π]`, so both signs of `σ`) and, for `σ_w`, over `|σ|` at fixed `θ`.
   **The objective is the consumer's** (choice 4): each scalar's error weighted by the smallest
   power of `θ` that multiplies it in `V`, `Q` and `Γ₂(M)ρ` — `θ⁰` for `v₀`, `w₀` and `v₁` (`Q`'s
   `-v₁Ξ`, which nothing suppresses), `θ¹` for `v₂` (`ΞW`) and `w₁`, `θ²` for `w₂` and `∂_z v₁`, `θ³`
   for `∂_z v₂` — over `u` times the block's scale `max(1, v₀)`, the maximum over the rays: the
   quantity SM.7(c)–(d) measured. Candidates within `max(0.25, 2%)` of the best are equal and the
   cheapest wins: fewer terms, then the larger series region (`w₀`: the smaller, since its exact
   form reuses `expm1 σ`).
5. **Domains, `NUMERICS.md` §12**, SM.13's rows replacing "`σ` finite":
   - `exp`, `adjoint`, `ad`, `jr`, `jl`, `act`, `act_jacobians`: every `τ`, `X` with
     `|σ| < ln MAX` (`709.78`, `88.72`), by the scaled arrangement.
   - `log`: `θ(X) ≤ π`, `|σ| < ln MAX`; a function of `q`'s sign at `θ = π` (§3.2).
   - `jr_inv`, `jl_inv`: `θ < 2π`, `|σ| < ln MAX`. The derivative of `log` and `⊖`: `θ(X) < π`.
6. **Corpus** (`PHASE1.md` §4.3), both precisions with `@f32` twins: `sim3_exp`, `sim3_log`,
   `sim3_jr`, `sim3_jl`, `sim3_jr_inv`, `sim3_jl_inv`, `sim3_ad`, `sim3_act` — the four Jacobians,
   as every group's, where `PHASE5.md` §3 listed two.
   - Strata: a grid `σ ∈ {0, ±10⁻⁸, ±10⁻⁴, ±10⁻², ±0.3, ±1, ±3}` × `θ ∈ {10⁻⁸, 10⁻⁴, 10⁻², 0.3, 1,
     3, π − 10⁻⁶}` named `sigma:<σ>/theta:<θ>`, six records each, random axes and `ρ ~ N(0, I)`;
     the joint-limit rays `ray:<r_Δ>` for `r_Δ ∈ {10⁻¹², 10⁻⁸, 10⁻⁴, 10⁻², 0.3, 0.8, 1.6}` over
     eight angles; and large `σ`, `±10, ±30, ±100, ±700` (`f64`), `±10, ±30, ±80` (`f32`), at
     `θ ∈ {10⁻⁴, 1, 3}`.
   - References: `mp.expm` of the `4 × 4` for `exp`; **the geometric SO(3) logarithm and a `V`
     solve for `log`, never `mp.logm`**; the dense series of `ad` for the Jacobians (the top-right
     block of `mp.expm([ad I; 0 0])`, SM.7's) and `mp.inverse` for their inverses; at
     `60 + |σ|` digits where the series cancel.
   - Scoring: `sim3_exp` and `sim3_act` as `sen3_exp` and `so3_act`, `σ` exact; `sim3_log` as
     `sen3_log`; the Jacobians and `Ad` by Frobenius error over `max(‖J‖, 1)`.
7. **Plantings** (`PHASE1.md` §10): `V` built as SE(3)'s `J_l(φ)` (the scale ignored), which must
   fail every `σ ≠ 0` stratum of `sim3_exp`; `Q` without its `∂_z` terms, which must fail
   `sim3_jr` wherever `θ ≥ 10⁻⁴`.
8. **`Dual`** (`jacobians_match_dual_*`): the derivative of a Δ-series scalar loses one power of
   `Δ` through the closed arm (SM.7(f)); the tolerance is per stratum, recorded from a measurement.

## Measured (plan step 2, #137)

The sweep (`xtask/src/thresholds/sim3.rs`) over 3549 points per precision, 16 radii per decade of
`Δ ∈ [10⁻¹⁶, 10]` × 13 angles, references by the total-degree recurrence summed to convergence and
cross-checked against the closed forms at 250 digits to `10⁻⁹⁰` (`sim3_ref.jsonl`):

| precision | `K` | `r_×²` | `Δ_b` | `σ_w²`, terms | objective | series | closed | `w₀` |
|---|---|---|---|---|---|---|---|---|
| `f64` | 13 | 0.2371 | `= r_×²` | 0.5829, 15 | 4.72 u | 4.72 | 3.04 | 1.63 |
| `f32` | 8 | 0.4217 | `= r_×²` | 0.5829, 8 | 2.93 u | 2.83 | 2.93 | 1.63 |

`Δ_b = r_×²`: form b of `v₁` wins wherever the closed arm runs. Past the rays, the scaled
arrangement keeps every scalar finite to `σ = 700`; six are within `2.5u`, and `w₁ = v₂ - σw₂`,
`∂_z v₁ = e^σ a' - σ ∂_z v₂` grow as they subtract `σ` times their own size, to `1620u` and `183u`
at `σ = 700` (`2.3σ u`). A rearrangement would remove it at a generated `|σ|` switch, but the `Q`
block's own assembly already loses `0.6σ² u` there (SM.7(g)(iii)), so no consumer would see it;
*Further work*.

## Rationale

- **`Sim3Jac`.** Every Jacobian a solver forms is a product or an inverse of `Ad` and `J`; at `343`
  against `117` multiplications and with `Matrix × Matrix` just measured at a tenth of its old cost
  ([`0069`](./0069-the-sphere-review-and-the-matrix-product.md)), the structure is the hot path. The
  dense twins keep D6.
- **The `Δ`-series over the split.** One branch variable that also covers `Δ = 0` without a
  division, one recurrence, and a tail bound that depends on `r_Δ` alone; the split needs a select
  at `Δ = 0`, three one-variable rows each with its own `|σ|` branch, and loses `145u` at
  `σ = 100` (SM.7(h)). Its cost is a new sweep shape, which `0004` requires either way.
- **The consumer's objective.** The coefficients exist for the blocks; a `v₂` that keeps no digit at
  `r_Δ = 10⁻⁸` is harmless where `θ²` multiplies it, and minimizing its own error would push the
  switch where no consumer is better off.
- **The scaled arrangement** costs one `exp(−σ)` and a select on the sign of `σ`, and is the only
  arrangement measured finite to `ln MAX`.

## Consequences

- `NUMERICS.md` §9 is stated (no longer *owed*), §12's `Sim3` row is replaced by §5's rows, §14
  gains `Sim3Jac::{mul, inverse, apply, apply_transpose, sandwich}` and `Sim3::jr_inv` rows.
- `PHASE5.md` §3 is amended: `Sim3Jac`, the eight ids, the `log` reference.
- `xtask thresholds` gains a two-variable sweep (rays in `(σ, θ)`); `coeffs/generated.rs` gains
  the Sim(3) group's switches and constants.
- Sim(3) has no in-process oracle: `sophus_lie` 0.15 has none. The container runners (`PHASE1.md`
  §0.0, owed) would bring Sophus's `Sim3`; until then the ids carry the no-regress bar alone.

## Implementation plan

1. This record and the spec edits (#136). Verified by `just lint`, `just doc`.
2. The coefficient group: `sim3_coeffs`, `sim3_v`, the three arms and the scaled arrangement; the
   ray sweep in `xtask thresholds`; the generated switches and constants; tests against mpmath
   points on the rays at both precisions. Verified by `just thresholds-check`, `just test`.
3. `Sim3`, `Sim3Tangent`, `Sim3Jac`, the twins; `LieGroup`, `Manifold`; `laws_for!` and the
   `Dual` laws. Verified by `just test`, `just no-std`.
4. The corpus ids, the generator, the subject, coverage, the plantings. Verified by
   `just corpus-check`, `just conformance` at both precisions, `--self-test`, `just envelope`.
5. Benchmarks `sim3/*` per stratum and a status update of `PHASE5.md` §0.0.

## Open questions

None.

## Further work

- A Sophus `Sim3` runner, for a domination bar.
- `w₁` and `∂_z v₁` for `σ ≫ 1` from the other rows of SM.5(a), SM.11's systems
  (`w₁ = (v₁ - w₀ + z w₂)/σ`), behind a generated `|σ|` switch, if a consumer at large positive `σ`
  ever needs the coefficient and not the block.
- `sim3_geodesic`, and the Sim(3) charts, as `PHASE4.md` and `PHASE5.md` §1 extend.
