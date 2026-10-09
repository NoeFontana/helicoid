# 0065: A Gaussian names its side and is stored symmetric

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** (filled in as work lands)

## Context

`PHASE5.md` §5 specifies `Gaussian<S, G, Sd, D>` with `to_left`/`to_right`
(`Σ_L = Ad_μ Σ_R Ad_μᵀ` by `Jac::sandwich`), `propagate(&self, j: &G::Jac, mean: G)` and
`mahalanobis_sq(&self, x) -> (S, S::Mask)`. Its maths is GG.7–GG.14
(`docs/maths/gamma-gaussian.md`). The change of side is exact (GG.8), so its error is rounding
alone. At a large translation that rounding is not small on the translation block
(`≈ u‖t‖²σ_φ²/σ_ρ²`, GG.13(c)).

`docs/maths/index.md` leaves four questions open on §5:

1. Does `sandwich` symmetrise? `(Ad Σ) Adᵀ` is asymmetric by up to `1.2u` of its largest entry,
   and `chol` reads one triangle.
2. `propagate`'s `G::Jac` is one type for `J_r`, `J_l` and `Ad`, so the side of `j` goes
   unchecked.
3. The twin of `to_left`/`to_right` (§14) compares two roundings, and needs a per-block tolerance
   rather than one `κu`.
4. No corpus id exists for `Gaussian`.

### Measured

**The round trip.** `to_right(to_left(Σ))` is `B (A Σ Aᵀ) Bᵀ` with `A = Ad_μ` and `B = Ad_{μ⁻¹}`,
each formed from the group element. Each sandwich obeys Higham's componentwise product bound,
`γ_{2D}` (Thm 3.5 twice). `B A` is `I` only to the rounding of the two `Ad` formations: the
quaternion's matrix, `μ⁻¹`'s translation, and `[t]_× R`, 21 roundings in all (counted in
`gaussian_tests`' `ROUND_TRIP_AD`). The bound is therefore

    |to_right(to_left(Σ)) − Σ| ≤ γ_{4D+42} · (|B||A|) |Σ| (|B||A|)ᵀ,   componentwise.

It is order-independent, so it holds for the structured `SEn3Jac::sandwich`. `|B||A|` has entries
of size `‖t‖`, so it reproduces GG.13(c)'s `‖t‖²σ_φ²` on the translation block. The worst fraction
of the bound used, over 10⁵ draws per group (poses up to `‖t‖ = 10³`, `Σ = D₀ H D₀` with
`σ_φ ≪ σ_ρ` among them):

| | SE(3) | SE₂(3) | SO(3) | SO(3) × ℝ³ |
|---|---|---|---|---|
| `f64` | 0.212 | 0.175 | 0.269 | 0.212 |
| `f32` | 0.207 | 0.195 | 0.245 | 0.207 |

Without the `+ 42`, SO(3) reaches 1.21: the `Ad` formation is part of the bound, not slack.

**The twins.** `to_left`/`to_right` against `reference::sen3jac_sandwich` with `J = Ad`, bounded
componentwise by `2γ_{2D} |A||Σ||A|ᵀ`. That bound is per entry, so per block: a translation entry
of size `σ_ρ²` is not held to the scale of a `‖t‖²σ_φ²` one. The worst ratio is 0.345 (SE(3)) and
0.277 (SE₂(3)).

**`d²`.** Two tests:
- Against a dense shadow at `f64`, from `f32` inputs: within 0.15 of GG.12's `c_D u/(1 − ρ)`.
- Side invariance (GG.10(a)) at lever 1, `σ_φ = σ_ρ`: worst 32.8 `u`, held under 64.

**The corpus.** Ids `gaussian_mahalanobis_se3` and `_se23`, 32 records per stratum. The table
gives the worst relative error of `d²` in `u`. Beside it, the stratum's largest **componentwise
condition number** `κ = Σᵢ |∂d²/∂xᵢ · xᵢ| / d²`, over every input `xᵢ`: both quaternions, both
translations and `Σ`'s lower triangle. It is computed by central differences at 60 digits (a
script that is not committed).

| Stratum | `se3` `f64` (`κ`) | `se3` `f32` | `se23` `f64` (`κ`) | `se23` `f32` |
|---|---|---|---|---|
| `gauss:right` | 340 (3.45e3) | 271 | 447 (3.03e3) | 329 |
| `gauss:left` | 2 072 (3.06e4) | 3 369 | 2 570 (8.47e3) | 1 477 |
| `gauss:corr-1e-2` | — | 2 249 | — | 391 |
| `gauss:corr-1e-4` | 1.10e4 (4.86e4) | 4 991 | 1 536 (2.42e4) | 4 711 |
| `gauss:corr-1e-8` | 2.54e7 (1.54e8) | — | 3.98e6 (1.33e8) | — |
| `gauss:lever` | 3.35e4 (1.42e5) | — | 3.21e4 (2.62e5) | — |
| `gauss:near-cut` | 5.07 (11.4) | 5.10 | 5.34 (14.9) | 6.20 |
| `gauss:indefinite` | mask clear, as the reference | same | same | same |

Every `f64` maximum is at most **0.44 of the stratum's largest `κ`**. The large numbers are the
problem's, not the algorithm's:
- A residual of two absolute poses carries an absolute error of order `u` (rotation) and
  `u‖t‖` (translation), against a `δ` of the size of `σ`, which is down to `10⁻³` here.
- On `corr-*`, `λ_min(H)` enters as GG.12 says.
- At `near-cut`, where `κ ≈ 11`, the algorithm reads 5 `u`.

No oracle answers these ids, so domination has nothing to compare against and they are gated by
no-regress and coverage alone.

## Decision

1. **The type** is `PHASE5.md` §5's.
   - `new(mean, cov)` is the only constructor (`_side` is private). It asserts at monomorphization
     that `D == G::DOF == G::Tangent::DOF`, through the fully qualified paths, since `Manifold`
     also has a `DOF` (`0060`).
   - The type is `Copy`.
2. **Only the lower triangle of `cov` is read, and every value a method returns has its upper
   triangle copied from its lower one**, so a returned `cov` is exactly symmetric.
   - This is `D(D−1)/2` moves and no arithmetic.
   - It answers question 1 here, not in `Jac::sandwich`, whose bits other callers keep.
   - A caller that writes `cov` directly stays correct, because the methods read the lower
     triangle.
3. **`to_left`** is defined on `Gaussian<_, _, Right, _>` and returns `Gaussian<_, _, Left, _>`. It
   is `mean.adjoint().sandwich(cov)`.
4. **`to_right`** is defined on `Left` and returns `Right`. It is
   `mean.inverse().adjoint().sandwich(cov)`: the adjoint of the computed inverse, not
   `Jac::inverse`, which would round a second time.
5. **`propagate(&self, j, mean)`** is `j.sandwich(cov)` about the new mean.
   - Its documentation states the contract the type cannot check, answering question 2: `j` is the
     Jacobian on side `Sd`, `D^R F` for `Right`.
   - A Jacobian of the other side converts both ends silently (GG.9, V7 of GG.14).
   - A typed Jacobian side would be a change to `LieGroup` and is not made here.
6. **`mahalanobis_sq(&self, x) -> (S, S::Mask)`.**
   - `δ = Sd::minus(x, mean)`, `(L, mask) = chol(cov)`, `y = solve_lower(L, δ)`, and the result is
     `‖y‖²`. `Σ⁻¹` is never formed.
   - The mask is `chol`'s, unaltered. A clear mask means `d²` means nothing (GG.11(c)).
   - The domain is `θ(μ⁻¹x) < π` (left: `θ(xμ⁻¹) < π`), stated and not asserted (GG.10(d)).
   - The rustdoc recommends gating on the right side (GG.10(a), GG.13(b)).
7. **Tests** (`gaussian_tests.rs`):
   - the round trip within `γ_{4D+42} (|B||A|)|Σ|(|B||A|)ᵀ` on SE(3), SE₂(3), SO(3) and
     SO(3) × ℝ³ at both precisions;
   - the §14 twins within `2γ_{2D}|A||Σ||A|ᵀ` (question 3);
   - `d²` against GG.12, and side invariance;
   - mask passthrough, exact symmetry, and `propagate` by the identity changing nothing;
   - a `compile_fail` doctest for `D ≠ DOF`.
8. **Corpus ids** `gaussian_mahalanobis_se3` and `_se23` (question 4).
   - **Inputs:** the mean `(q0, x0)`, the point `(q1, x1)`, `side` (0 right, 1 left) and `Sigma`.
     The point is `μ Exp(ξ)` or `Exp(ξ) μ` with `ξ ~ N(0, Σ)` drawn at 120 digits, formed and
     rounded, so `d²` is a `χ²` draw.
   - **`Σ`:** `D₀ Q diag(s) Qᵀ D₀`, with `σ_φ`, `σ_ρ` log-uniform in `[10⁻³, 1]` and `Q` and the
     rounding `gen.chol`'s. `valid` is decided exactly, by an `LDLᵀ` in rationals.
   - **Reference:** `δ` by the geometric `Log` (never `mp.logm`, `0045`), then `δᵀΣ⁻¹δ` by
     `mp.lu_solve`. `d2 = 0` where `valid = 0`.
   - **Strata:**
     - `gauss:right` and `gauss:left`, with `s` in `10^U[−1, 0]`;
     - `gauss:corr-1e-{4, 8}` (binary64) and `gauss:corr-1e-{2, 4}` (binary32);
     - `gauss:lever` (binary64): the left form of a right `Σ` at `‖t‖ = 10³`, GG.13;
     - `gauss:near-cut`: residual `θ ∈ [2.5, 3]`;
     - `gauss:indefinite`.
   - **Scoring:** `Rule::Masked { mask: "valid", fields: [d2 relative] }`, as `chol`'s. Answered
     by `helicoid` at both precisions and by `helicoid:host-std` at binary64. Both ids are in
     `required()` with their files.
9. **`NUMERICS.md`.** §11 gains the `d²` rule, with the condition number above as its reading. The
   §14 row names the twin's bound.

## Rationale

- **Symmetry in `Gaussian`, not in `sandwich`.** `sandwich` is a `Jac` method with other callers
  and corpus-scored bits. A `Gaussian` is the value whose consumers factor and eigen-decompose
  `cov`, and the copy costs no rounding.
- **A componentwise bound, not `κu`.** GG.13(c) shows the round trip loses `‖t‖²σ_φ²/σ_ρ²` on the
  translation block. A normwise `κ(Ad)u` would either be too loose for the rotation block or
  falsely fail the translation one. Higham's componentwise bound tracks both, and is the bound a
  proof gives.
- **Relative `d²` with `κ` beside it, not a `κ`-scaled metric.** A `κ`-scaled metric would make the
  envelope read the algorithm alone. It would also need `κ` in every record (35 to 57
  derivatives per record) and a new metric rule, for ids that no oracle answers. The plain relative
  error makes no-regress meaningful, which is the only bar that can apply, and this record states
  the reading.
- **`propagate` as specified.** A side-typed Jacobian would change `LieGroup`'s associated types for
  every group. The documented contract and GG.14's numbers are the proportionate answer until a
  consumer's error says otherwise.

## Consequences

- `Gaussian` holds an exact-symmetry invariant, and every new method must keep it.
- `gaussian_mahalanobis_*` are no-regress-only. Their baseline rows will read up to `2.5·10⁷ u`
  at `corr-1e-8`, which is conditioning (0.16 of `κ`), not a defect.
- `mahalanobis_sq` factors `Σ` on every call. A consumer gating many candidates against one law
  pays `D³/6` per candidate (Further work).

## Implementation plan

1. This record, with the `PHASE5.md` §5, `NUMERICS.md` §11 and §14, `API.md` §2 and
   `docs/maths/index.md` edits. Verified by `just lint` and `just doc`.
2. `gaussian.rs` and `gaussian_tests.rs`. Verified by `just test`, `just no-std` and `just msrv`.
3. `gaussian_mahalanobis_{se3,se23}`: generator, cross-check (`d²_R = d²_L` through the exact `Ad`
   at 120 digits), corpus, subject, metric, coverage and `PHASE1.md` §4.3. Verified by
   `just corpus-check`, `just conformance` at both precisions, and `just envelope`.

## Open questions

None.

## Further work

- A factored form (`Gaussian::whitener() -> (L, mask)` or similar) for many-candidate gating, if a
  consumer's bench shows the per-call `chol` matters.
- A covariance converter between rotation-first and translation-first orders (`ΠΣΠᵀ`, SE.14(d)),
  the misreading GG.14 calls V2. It is the one with no counterpart in the type system.
