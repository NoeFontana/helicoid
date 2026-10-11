# 0067: the sphere takes plane coordinates, and both precisions

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #133

## Context

[`0063`](./0063-the-sphere-reads-its-sign-by-comparison-and-is-held-unit.md) settled S²'s four gaps
(the sign of zero, the unit norm, the near-antipode reference, the direction metric), and
`PHASE5.md` §2, `NUMERICS.md` §8 and `docs/maths/charts.md` CH.8–CH.11 give the formulas. Building
it left seven choices that none of them fixes, and its benchmarks found one defect outside it.

1. **Coordinates.** `PHASE5.md` §1.4 names S²'s Jacobian, `Mat2<S>`, and not its tangent. A
   `Tangent` impl on `Vector<S, 2>` would give the coordinates `Vector`'s `Add` and `Sub`, which
   `API.md` R1 keeps off every tangent (`RnTangent`'s `compile_fail` doctests).
2. **`renormalize`'s receiver.** `0063` decision 2 writes `S2::renormalize(&self)` "as
   `renormalize` is for the quaternion", whose step is `&mut self`. R2 makes `renormalize` the one
   `&mut self` numeric method and R3 keeps one spelling per operation.
3. **`retract`'s coefficients.** CH.9(a) reads `sin θ/θ = 2k cos(θ/2)` and `cos θ = 1 - 2(kθ)²`
   from `exp_coeffs`, "the reading assumed here": `PHASE5.md` §2 names no source.
4. **The Jacobian forms.** CH.10 gives `rj = B'ᵀ J_l(Bδ) B = Zᵀ A(δ)` and `lj = A(δ)⁻¹ Z`, and
   offers CH.10(c)'s `(v, w)` form of `Z` as "the form for `lj`". That form divides by
   `1 + n·m ≈ (π - α)²/2`: at `π - α = 10⁻⁴` it measured `1.5e9 u`, where `E = Exp(Bδ)` read
   `9.6 u` (CH.10, *Checked*). And `rj` reads the perturbed point in the frame frozen at
   `m = retract(δ)`: a basis computed from an `m` that differs from `retract`'s by an ulp can sit
   on the other side of `n_z = 0`, where the frames differ by the reflection `T` (CH.11(a)).
5. **Outside the domain.** `NUMERICS.md` §12 promises `local` an unspecified *finite* value at
   `m = -n`. There `s² = 0` and `w = -1`, and `log_ratio`'s exact arm returns `2π/0 = ∞`, times
   `n × m = 0`: NaN. On the far hemisphere `‖n × m‖²` also underflows before `m` reaches `-n`:
   at `n = e_z`, `m = (10⁻²⁰⁰, 0, -1)` is a unit double, not `-n`, and has `s² = 0`.
6. **The corpus.** `PHASE5.md` §2 names three strata and no counts and no precision. `0016` item 1
   gives an id `@f32` strata only when a record asks. `0063` decision 5 scores `s2_local` by the
   2-norm and states no floor.
7. **The plantings.** `PHASE1.md` §10 owes each corpus id a planted defect its bars detect.

8. **A defect outside S².** In a consumer loop over 1024 charts (one pinned core, `--release`),
   `S2Chart::local` took 46–61 ns at `f32` and 7.5–18 ns at `f64`, while `SO3::log`, which calls
   the same `log_ratio`, is faster at `f32` than at `f64` (3.6–10.0 ns against 4.1–14.8). The
   disassembly of `local::<f32>` shows each `Vector::dot` keeping its two `array::IntoIter`s on the
   stack: the iterator's alive range is stored, then every lane is reloaded through it, a
   store-forwarding chain per dot. At `f64` LLVM unrolled the same loop. `Vector::dot` was
   `sum(self.0.into_iter().zip(o.0).map(..))`; `Matrix`'s `Mul` sums the same way by index and
   does not show it.

## Decision

1. **`S2Chart<S>: Chart<S, S2<S>>` with `Tangent = RnTangent<S, 2>` and `Jac = Mat2<S>`**, through
   `impl Jac<S, RnTangent<S, 2>> for Mat2<S>` (`inverse` by `Mat2::inverse_adj`, `debug_assert!`ed
   non-singular as `Mat3`'s is; `sandwich` in `Mat3`'s order). `S2<S>: Manifold<S>` with
   `Chart = S2Chart<S>`. No new tangent type.
2. **`S2::renormalize(&mut self)`**, the quaternion's receiver; `0063`'s `&self` is read as that.
3. **`retract(δ) = cos θ n + ς (2k cos(θ/2)) (δ₂ b₁ - δ₁ b₂)`** with `cos θ = 1 - 2(k·k·θ²)`, one
   `exp_coeffs` at `θ² = ‖δ‖²`: CH.9(a)'s reading. The basis is the `H` form with `νᵀν` computed
   (CH.8(e)).
4. **The Jacobians.**
   - `local(m) = Bᵀ p`, `p = ½ log_ratio(s², w) · (n × m)`, CH.11(e)'s more accurate arrangement
     (1.4–1.9 against 2.3–3.8 `u/(π - α)`).
   - `retract_jacobian(δ) = B'ᵀ J_l(Bδ) B`, `B'` the basis at the `m` that `retract` returns, from
     the one private `step` both call, so bit-identical. `J_l bᵢ` is written out by
     `b₁ × b₂ = ς n` (CH.8(c)); `(a, b) = jr_coeffs(θ²)`.
   - `local_jacobian(m) = (α/s)(I - b δδᵀ) Z`, `Zᵢⱼ = (E bᵢ)·b'ⱼ`, with
     `E bᵢ = w bᵢ - (bᵢ·m) n + a δᵢ p`: Rodrigues for `E = Exp(p)`, reduced by `bᵢ ⊥ n`. It never
     divides by `1 + w`. `(a, b) = jr_coeffs(‖δ‖²)`.
   - `local_with_jacobian` forms `n × m`, `log_ratio` and `p` once.
5. **The safe argument.** `log_ratio` takes `s² = 1` where `s² = 0` and `w ≤ 0`, so `local(-n) = 0`
   and its Jacobian is finite. `NUMERICS.md` §12's row gains the far hemisphere's bound: where
   `n·m ≤ 0`, `‖n × m‖²` normal (`m` farther than about `1.5e-154` (`f64`), `1.1e-19` (`f32`)
   from `-n`); closer, an unspecified finite value. Nothing is asserted.
6. **The corpus.** 64 records per stratum. `s2:generic` and `s2:nz0` draw `θ` log-uniform in
   `[10⁻⁸, π - 0.1]` and a uniform direction; `s2:nz0`'s base is `(cos λ, sin λ, ±0)`, the sign of
   zero alternating; `s2:near-antipode` takes `θ = π - 10⁻ᵏ`, `k ∈ {1, 2, 3, 6, 9, 12}` at `f64` and
   `{1, …, 6}` at `f32` (binary32 cannot hold `π - 10⁻⁷` below `π`). Each stratum has its `@f32`
   twin, the same draws rounded once. `s2_local`'s `m` is the reference `retract` of the same
   record, rounded. `s2_retract` is scored by `0063` decision 5's direction metric; `s2_local` by
   the 2-norm with a floor of 1, as `se3_*_local`'s rotation is: two absolute unit vectors in, so a
   nearby pair cancels to a `δ` whose own size is the wrong denominator.
7. **Two plantings** (`xtask/src/conformance/selftest_s2.rs`): `ς` by the sign bit must fail both
   `s2:nz0` strata of both ids and leave every other S² stratum's maximum unchanged to the bit;
   `§8`'s ratio as first worded, `r` for `½ r`, must fail every `s2_local` stratum where `1024`
   times the clean maximum is below `1/u`.
8. **`Vector::dot` sums by index**, `sum((0..N).map(|i| self.0[i] * o.0[i]))`: the same terms in the
   same order, so every output of the workspace is unchanged to the bit.

## Measured

**Chart laws** (`s2_tests::measure_s2`, 10⁶ seeded cases, 10⁵ for `Dual`; bounds twice these):

| law | `f64` | `f32` |
|---|---|---|
| `retract(0) = n`, `local(n) = 0` | 0 | 0 |
| `local(retract(δ)) = δ` | 6.000 | 6.083 |
| `retract(local(m)) = m` (angle) | 16.521 | 17.051 |
| basis: `BᵀB = I`, `Bᵀn = 0`, `b₁ × b₂ = ςn` | 6.000 | 5.786 |
| `retract`'s drift, `‖m‖² - 1` | 10.877 | 8.766 |
| `local_jacobian · retract_jacobian = I` | 8.834 | 9.035 |
| twin, over `local`'s conditioning `α/sin α` | 5.870 | 6.034 |
| both Jacobians against `Dual` | 9.833 | — |

The twin is scored over the conditioning because it halves the angle through `1 + n·m`, the form
`0063` decision 3 rejects: unscaled it read 1684 `u` on the draw nearest `-n`.

**`renormalize`** from `| ‖n‖² - 1 | = 2^ℓ U(½, 1)`, `ℓ` up to the band's edge, 20 000 draws:
`5.005 u` (`f64`) and `4.447 u` (`f32`). `0063`'s table gives the `f64` worst as `5`, rounded; the
step's own truncation `¾η²` is `1.0 u` at the edge. The test bound is `5.5 u`.

**Conformance**, `max_u` per stratum, `helicoid` and its host-`std` twin identical:

| id | precision | `s2:nz0` | `s2:near-antipode` | `s2:generic` |
|---|---|---|---|---|
| `s2_retract` | `f64` | 4.159 | 5.965 | 2.740 |
| `s2_retract` | `f32` | 2.181 | 4.750 | 1.859 |
| `s2_local` | `f64` | 1.587 | 2.494e11 | 1.062 |
| `s2_local` | `f32` | 1.777 | 3.115e5 | 1.887 |

`s2_local`'s near-antipode maxima sit at `π - α = 10⁻¹²` and `10⁻⁶`, where the conditioning
`α/sin α` is `3.1e12` and `3.1e6`: `0.08` and `0.10` of it. The plantings read `1.8e16` and
`3.3e7 u` (sign bit, `s2:nz0`) and `9.0e15` and `1.7e7 u` (full ratio); `s2:near-antipode@f32` is at
its conditioning and not judged. The envelope is unchanged: the same 148 domination failures, none
on S², and coverage holds with both ids required.

**Time**, `f64`/`f32` ns per call, a consumer loop over 1024 charts on one pinned core (the
criterion rows `s2/*` time each method behind a call):

| | `θ = 7.5e-8` | `θ = 1` | `θ = 3.14159` |
|---|---|---|---|
| `retract` | 4.8 / 3.6 | 8.7 / 5.7 | 16.7 / 9.6 |
| `retract_jacobian` | 17.0 / 14.1 | 27.0 / 21.5 | 43.2 / 33.1 |
| `local` | 7.8 / 7.4 | 17.8 / 13.9 | 14.7 / 8.8 |
| `local_with_jacobian` | 24.9 / 22.5 | 49.3 / 40.0 | 62.8 / 40.7 |

At `f32` before decision 8, `local` read 46.4 / 61.0 / 47.9 ns and `local_with_jacobian` 132.8 /
159.5 / 159.2. `reference::s2_local` reads 1.3–1.6 times `local` at `f64` (criterion rows).

**`Vector::dot`**, the bench gate (`cargo xtask bench-gate --bench groups`, baseline the `main`
binary, one pinned core): GATE.

## Rationale

- **`RnTangent<S, 2>`.** It is already the public `ℝⁿ` coordinate type, with R1's guarantee and its
  `read_dense`/`write_dense`; a second two-component tangent would be a type with nothing to add.
  Its field is named `rho`, a translation's name, which is the one cost.
- **`Exp` for `lj`, not `(v, w)`.** Both are on the maths page; only one keeps `lj` at its
  conditioning near the antipode, and it costs one `jr_coeffs`, which `A⁻¹` needs anyway.
- **One `step` for `retract` and `retract_jacobian`.** The frame of `m` is discontinuous at
  `n_z = 0`; the Jacobian is only meaningful in the frame a caller will freeze at `retract(δ)`, so
  both must read the same bits, which a shared function guarantees and a duplicated formula does
  not.
- **`@f32` twins.** Bearings in a camera pipeline are `f32` as often as not, and the precision
  that differs is the one with the switch points and the conditioning: the twins are 384 records.
- **`dot` by index, not a hand-written sum per arity.** The index form is `Matrix`'s, already
  measured clean, and fixes every `N` and every caller at once.

## Consequences

- `PHASE5.md` §0.0 and §2, `NUMERICS.md` §12 and §14, `API.md` §2, `PHASE1.md` §4.3 and §10 change
  as above; `0063`'s *Implementation* line names this record.
- `s2_retract` and `s2_local` are required by coverage, out of `OWED`.
- `0016` item 1's `@f32` strata extend to both S² ids.
- The criterion rows `s2/{at, retract, retract_jacobian, local, local_with_jacobian, local_twin}`
  per `θ` stratum and precision join `benches/groups.rs`.

## Implementation plan

1. One PR (#133): `S2`, `S2Chart`, `Mat2`'s `Jac`, the twin, the tests; the generator's
   `gen/s2.py` and its checks; the corpus files; the metric, subject, coverage and plantings in
   `xtask`; `Vector::dot`; the docs above. Verified by `just lint`, `just test`, `just doc`,
   `just no-std`, `just msrv`, `just corpus-check`, `just conformance` at both precisions,
   `just conformance --self-test` and `just envelope` (only its pre-existing failures).

## Open questions

None.

## Further work

- `lj` near the antipode is checked against `Dual` on generic draws only; CH.10 leaves its `f64`
  accuracy there open, and no corpus id scores a chart Jacobian.
- locus-tag's bearings can take `S2` and the direction metric (`0063` *Further work*).
