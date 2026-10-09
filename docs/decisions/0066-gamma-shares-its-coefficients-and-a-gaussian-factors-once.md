# 0066: Γ shares its coefficients, and a Gaussian factors once

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #132

## Context

A review of the Phase 5b stack (#126 → #131, `0064` and `0065`) found one false statement, one
undeclared domain, and five costs on the paths a fusion filter runs per IMU sample or per
association hypothesis:

1. **The switches are on `θ²`.** `0064` *Measured* and `NUMERICS.md` §7 say that on `θ ≤ π` "the
   generated switches put every call on the series arm". `generated.rs` puts `a` at `θ² < 8.35`, `b`
   and `d` at `θ² < 9.65` in binary64, and `b` at `θ² < 7.77` in binary32. So the exact arm runs
   from `θ = 2.89` (`M = 1`) or `3.11` (`M = 2`) in binary64, and from `2.79` in binary32. The
   `≤ π` column of `0064`'s table was measured with those draws included.
2. **`reference::gamma_apply(m, …)`** passes `m!` to `S::lit`, which requires an exact value and
   `debug_assert!`s it in binary32. `14!` is not a binary32, so `m = 14` panics a debug build, and
   in binary64 `m ≥ 19` rounds silently. Its `# Domain` stated `θ ≤ π` only.
3. **`gamma_apply_jacobian` ran the coefficient kernel on `Dual<S, 3>`.** Each Horner step carried
   three derivative lanes, and `v` was lifted to constants whose zero lanes the cross products then
   multiplied. The coefficients depend on `z = θ²` alone, so one lane, seeded at `z` (the pattern
   `0062`'s `sinc` uses), carries the same derivative.
4. **A preintegration step evaluated the coefficients twice.** `Δp = Δt² Γ₂ f` needs `Γ₂ f`,
   `∂/∂φ` and `∂/∂f = Γ₂`. The API made it call `gamma2` and `gamma_apply_jacobian`, two grouped
   branches on the same `θ²`.
5. **`J_r`, `J_r⁻¹` and `Γ₂` were assembled as three dense matrices**, `(I s0 + W s1) + W² s2`.
   LLVM may not fold the structural zeros (`0 · x` is not `0` for a NaN `x`, `x + 0` is not `x` for
   `−0`, and D16 forbids fast-math), and `W²`'s symmetric half was formed twice: 36
   multiplications and 27 additions.
6. **`Gaussian::mahalanobis_sq` factored `Σ` on every call.** A `χ²` gate of many candidates
   against one track paid `chol` per candidate. `0065` *Further work* named a factored form.
7. **`cov` was a `pub` field**, so the exact symmetry of `0065` decision 2 held only for values the
   methods returned, and every method copied the mirror of `cov` before reading it (`D²` moves per
   call).

The review also proposed a sandwich that computes and writes only the lower triangle of `J Σ Jᵀ`.

### Measured

**Accuracy.** `gamma_tests::measure_gamma_twin`: 20 000 draws per row, `θ` log-uniform on
`[1e-8, π]`, `0064`'s protocol. Worst Jacobian error against the §14 twin, in `u`:

| `M`, precision | `θ < 1e-2` | `< 1e-1` | `< 1` | `≤ π` |
|---|---|---|---|---|
| 1, `f64`: `Dual<S, 3>` → chain rule | 1.94 → 1.82 | 2.30 → 2.30 | 3.63 → 3.68 | 4.75 → 4.87 |
| 2, `f64` | 1.03 → 0.99 | 1.15 → 1.12 | 1.37 → 1.33 | 1.42 → 1.33 |
| 1, `f32` | 1.28 → 1.07 | 1.16 → 1.04 | 1.35 → 1.48 | 4.93 → 4.67 |
| 2, `f32` | 0.76 → 0.67 | 0.70 → 0.69 | 0.72 → 0.66 | 1.86 → 1.66 |

- The value column and the `Γ₁`, `Γ₂` matrix columns of `0064`'s table are unchanged to the bit.
- Three groupings of the chain rule's `φᵀ` terms were measured. They differ by less than 0.1 `u`
  in every row except `M = 1` near `π`, where they read 5.53, 4.96 and 4.87 at `f64`. The last is
  taken: `c φᵀ − 2σ₂ vφᵀ` kept apart, and added after `σ₂ φvᵀ`.

**Time.** `cargo xtask bench-gate --bench groups --against <baseline>`, pinned to one core, three
bracketed replicates. The baseline is this tree with the old bodies under the new signatures: the
`Dual<S, 3>` path plus a `gamma2` call; a `Whitener` that factors per point; and a gate that
calls `mahalanobis_sq` per point. Median candidate/baseline ratio:

| Rows (three strata where the row has them) | `f64` | `f32` |
|---|---|---|
| `so3/gamma_apply_jacobians_1` (against `Dual<S, 3>` plus `gamma1`) | 0.47–0.71 | 0.51–0.85 |
| `so3/gamma_apply_jacobians_2` (against `Dual<S, 3>` plus `gamma2`) | 0.58–0.85 | 0.49–0.79 |
| `so3/jr` | 0.71–0.87 | 0.77–0.97 |
| `so3/jr_inv` | 0.81–0.98 | 0.74–0.99 |
| `so3/rminus_jacobians` | 0.80–1.00 | 0.70–0.95 |
| `so3/gamma2` | 0.77–0.91 | 0.92–0.94, and 1.20 near identity (below) |
| `se3/jr`, `se23/jr` | 0.93–1.03 | 0.94–1.00 |
| `se3/jr_inv`, `se23/jr_inv`, `*/rminus_jacobians` | 0.98–1.00 | 0.96–1.01 |
| `se3/gaussian_{to_left,to_right,propagate}` | 0.79–0.86 | 0.85–0.90 |
| `se23/gaussian_{to_left,to_right,propagate}` | 0.95–0.97 | 0.97–0.99 |
| `{se3,se23}/gaussian_mahalanobis_sq` | 0.99–1.00 | 0.98–1.00 |
| `{se3,se23}/gaussian_gate_16` (16 points, one law) | 0.56–0.59 | 0.51–0.56 |

The A/A floor measured beside the comparison is 0.72% (median). Three rows are past theirs:
`f32` `se3/jr_inv/generic` at 1.005, `f32` `se3/rminus_jacobians/near-identity` at 1.008, and
`f32` `so3/gamma2/near-identity` at 1.20. The first two run `SO3::jr_inv` inside `SEn3`, which
reads 0.98–1.00 at `f64` and 0.74–0.99 on its own, so they are below 1% and are layout, not
arithmetic. The third is the following paragraph's.

**`f32` `so3/gamma2/near-identity` reads 1.19× to 1.21×** in every run. It does not measure
`Γ₂`. The routine is `gamma2_of`, an `#[inline(never)]`
wrapper that returns a 36-byte matrix through memory, about 7 ns in all. A scratch crate (not
committed, `__sweep`'s `gamma2_coeffs`, both assemblies verbatim, pinned core) times the two
three ways, median ns, old → new:

| | inlined, matrix returned | out-of-line call | inlined, `Γ₂ f` applied |
|---|---|---|---|
| `f32`, near-identity | 7.13 → 17.28 | 6.75 → 8.10 | 18.67 → 6.45 |
| `f64`, near-identity | 12.01 → 19.12 | 8.98 → 8.10 | 23.57 → 8.79 |
| `f32`, generic | 9.69 → 23.88 | 9.16 → 8.39 | 25.06 → 9.83 |

The first two columns move by how criterion's harness handles a returned matrix: an inlined
routine reads slower than the same routine behind a call. That is the codegen-context trap, not
the arithmetic. The third is a consumer's use, `Δp = Δt² Γ₂ f` inlined into its step, and the new
assembly is 0.35× to 0.39× there at both precisions.

**The lower-triangle sandwich.** It is not taken, by count. `SEn3Jac::sandwich` applies `J` by
`act`, which already skips `J`'s structural zeros: `15 + 33N` flops per application. Its second
stage is therefore 288 flops at `D = 6` and 729 at `D = 9`. A lower-triangle second stage cannot use
`act`, which returns whole rows. It needs `J` dense (`write_dense`, `D²` stores) and
`D(D+1)/2` dot products of length `D`: 231 flops at `D = 6` and 765 at `D = 9`. That is 20% less
arithmetic for SE(3) before the dense write, and more for SE₂(3). It is also a different rounding,
which would owe the twins a new measurement. The real waste in that path was the `D²` mirror copy
of item 7, which item 7's fix removes.

## Decision

1. **`NUMERICS.md` §7** states the exact-arm range above, and `0064` carries an amendment note.
   `gamma_apply_jacobians`' rustdoc states the same range.
2. **`reference::gamma_apply`**'s `# Domain` adds "`m!` exact at `S`": `m ≤ 18` in binary64 and
   `m ≤ 13` in binary32, `debug_assert!`ed. `gamma_tests::twin_domain` holds both ends.
3. **`HatSq`**, `pub(crate)` in `so3`. `HatSq::of(φ)` holds `W²`'s six distinct entries, and
   `HatSq::poly(s0, s1, s2) = s0 I + s1 W + s2 W²` assembles entry by entry: 15 multiplications
   and 12 additions in all.
   - `SO3::jr` is `poly(1, −a, b)`, `SO3::jr_inv` is `poly(1, ½, c)`, and `so3::gamma2` is
     `poly(½, b, d)`.
   - Each caller forms `HatSq::of(φ)` **before** its coefficient `branch`. With the products after
     the out-of-line call they were on the result's critical path, and `f32`
     `so3/gamma2/near-identity` read 1.19×. Before it, the core overlaps them with the call.
   - Each entry is the three-matrix sum's to the bit for finite operands. `W²` is `mul_hat`'s
     rounding: its diagonal `(−z)z + y(−y)` is exactly `−(z² + y²)`, and its off-diagonal pair is
     `xy` both ways. The dropped terms are exact zeros.
   - It differs where a dropped `±0` set the sign of a zero result, and where a non-finite `s1`
     reached the diagonal through `0 · s1`. `so3_tests::hat_poly_is_the_three_matrix_sum` holds
     both.
   - `jl` stays `jr(−φ)` and is `jr`'s transpose to the bit (`jl_is_jr_transposed_to_the_bit`).
     The closed `φφᵀ − θ²I` stays refused, for the reason `SO3::jr`'s comment gives.
4. **`so3::gamma_apply_jacobians::<M, S>(φ, v) -> (Γ_M v, ∂/∂φ, ∂/∂v = Γ_M)`** replaces
   `gamma_apply_jacobian`, which no release has shipped. The Jacobians come in argument order, as
   every `*_jacobians` returns them.
   - `(σ₁, σ₂)` and their `z`-derivatives come from one `jr_coeffs` (`M = 1`) or `gamma2_coeffs`
     (`M = 2`) evaluation on `Dual<S, 1>` seeded at `z = norm_sq(φ)`.
   - The value is `(s0 v + σ₁ φ×v) + σ₂ φ×(φ×v)`. It is the old value to the bit, because a `Dual`'s
     value lane is the plain evaluation (`0003`).
   - The Jacobian is
     `σ₂(φ·v) I − σ₁ v^ + ((σ₂φ) vᵀ + c φᵀ) − (2σ₂ v) φᵀ`, with `c = 2(σ₁′ φ×v + σ₂′ φ×(φ×v))`.
   - `Γ_M` is `HatSq::poly(s0, σ₁, σ₂)`: `gamma1` (`= SO3::jl`) or `gamma2` to the bit
     (`gamma_tests::the_third_output_is_gamma_to_the_bit`).
   - The twin proptests keep `0064`'s bounds and are renamed `gamma_apply_jacobians_matches_reference_*`.
5. **`Gaussian`'s `cov` is private**, read through `cov(&self) -> &Matrix<S, D, D>`.
   - `new` remains the only way in. It reads the lower triangle and mirrors it, so every
     `Gaussian`'s `cov` is exactly symmetric.
   - `to_left`, `to_right` and `propagate` sandwich `cov` as stored. Their outputs are unchanged to
     the bit for every value `new` built.
   - `mean` stays `pub`: it has no invariant.
   - `PHASE5.md` §5's normative struct changes with it.
6. **`Gaussian::whitener(&self) -> (Whitener<S, G, Sd, D>, S::Mask)`** factors `Σ` once, and
   **`Whitener::mahalanobis_sq(&self, x: &G) -> S`** is `⊖` and a forward substitution.
   - `Gaussian::mahalanobis_sq` is `0065`'s body unchanged: `⊖`, then `chol`, then the solve.
     `Whitener::mahalanobis_sq` is the same expression on the stored `L`, so the two agree to the
     bit, and the corpus ids `gaussian_mahalanobis_*` do not move.
   - `Gaussian::mahalanobis_sq` does not build a `Whitener`. Building one copies `L`, which
     measured 1.28× to 1.44× that call's time, and running `chol` before `⊖` measured 1.05× to
     1.16×.
   - The mask is `chol`'s, returned once. A clear mask voids every `d²` the `Whitener` returns.
   - `Whitener` is `Copy`; its fields are private.
7. **`PROJECT.md` §5.1** drops "Closed-form Γ directional Jacobians". The `Dual<S, 3>` path it
   gated is gone. A closed form for `σ′` would replace one `Dual<S, 1>` lane, which no measurement
   asks for.
8. **Benches.** `benches/groups.rs` gains `so3/gamma2`, `so3/gamma_apply_jacobians_{1,2}` per
   stratum, and `{se3,se23}/gaussian_{to_left,to_right,propagate,mahalanobis_sq,gate_16}`, each
   behind an `#[inline(never)]` call where the routine is small enough to be inlined.

## Rationale

- **The chain rule over `Dual<S, 1>`, not a closed `σ′`.** The lane is the derivative of the arm
  that runs, which `0062` measured to be the better reading for `sinc` (up to 2.3× against
  `(b − a)/2`). It needs no `τ₆`, which `𝒥₂`'s closed form does (`0064`). It also costs one lane
  in place of three.
- **One function with three outputs, not a fourth function.** A preintegration step needs all
  three. A caller that drops `Γ_M` pays a 27-flop assembly at most, and nothing once inlined. Two
  functions over the same evaluation would be two spellings (`API.md` §6.5).
- **A separate `Whitener`, not `L` cached in `Gaussian`.** Caching would double `Gaussian`'s size
  and make `to_left` and `propagate`, which never read `L`, pay a `chol`.
- **A private `cov`, not more copies.** The invariant is `0065`'s, and only a private field can
  hold it for every value. That also removes the reason for the copies.

## Consequences

- `gamma_apply_jacobian` is gone, before any release shipped it.
- `Gaussian { cov, .. }` and `g.cov = …` no longer compile outside the crate. `cov()` reads it,
  and `Gaussian::new` builds a new one.
- `SO3::jr`, `jr_inv`, `jl`, `jl_inv`, every `SEn3` Jacobian built on them, and `gamma2` change
  only the sign of some exact-zero entries. No corpus row moves (`just conformance`, then
  `just envelope --check`).

## Implementation plan

1. This record with every code and doc edit, one PR stacked on #131. Verified by `just test`,
   `just lint`, `just doc`, `just no-std`, `just msrv`, `just conformance` (both precisions),
   `just envelope --check`, and the bench gate above.

## Open questions

None.

## Further work

- A lower-triangle `Jac::sandwich` variant, if a consumer's bench at `D = 6` shows the ~20% of
  the second stage matters. It owes the twins a re-measurement.
