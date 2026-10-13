# 0085: `from_matrix` is measured on exact rotations

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** step 1, this record's PR

## Context

locus-sensors reported that `SO3::from_matrix` in 0.0.3
([`0079`](./0079-from-matrix-normalizes-a-column-of-k.md)) returns, for an exact quarter turn
(90° about `z`, a signed permutation matrix), a quaternion whose two `√½` components are each one
ulp below the nearest double: `‖q‖² − 1 ≈ −1.6u` at `f64`, which every later `SO3::act` carries as
an error of about `2|η|‖v‖`.

The corpus could not see it, for three reasons:

1. `so3_from_matrix`'s inputs are `R(q)` rounded to binary64 over the quaternion strata. No exact
   rotation matrix is among them.
2. The id has no `@f32` strata.
3. The id is not scored. `NUMERICS.md` §11 reported matrix inputs "by backward error only", which
   was never built, so the harness checked shape and finiteness (`metric::Rule::BackwardOnly`).

`0079` measured its kernel outside the harness for that reason and named this instrument as further
work. [`0006`](./0006-the-instrument-comes-first.md) puts the instrument before the kernel fix.

The backward-only reading had one justification: on a rounded input the reference is undetermined
below `O(d_M)` (`maths/error-analysis.md` EA.9(d)). Conditioning was not it. `M ↦ ±q⋆` is
`1/(2√2)`-Lipschitz at every θ, π included (EA.9(a)). Only the sign of `q` is ill defined near π,
and a sign-free metric does not read it.

## Decision

1. **The metric.** `NUMERICS.md` §11's direction metric, in ℝ⁴ and sign-free:
   `max(atan2(‖q̂ ∧ q‖, |q̂·q|)/u, |‖q̂‖² − 1|/(2u))`, against the quaternion of the input's polar
   factor. The wedge is formed exactly as `‖q̂‖²‖q‖² − (q̂·q)²`, the integer identity that is
   `‖q̂ × q‖²` in ℝ³, so `s2_retract`'s scores do not move. In `xtask::conformance::metric`,
   `Rule::Axis("q")` replaces `Rule::BackwardOnly` for `so3_from_matrix`. Nothing else was
   backward-only, so `BackwardOnly`, `finiteness` and `Score::Unscored` go. An unscored candidate
   row then fails the envelope for every id.
2. **`matrix:exact`.** After the quaternion strata: the 24 signed permutation matrices with
   determinant +1, built from integers. `R(q)` of an mpf quaternion would leave residues of
   `1e-120` that rounding keeps. Each reference is exactly 0, ±½, ±√½ or ±1, written to 30 digits,
   and the polar-factor cross-check passes with `d_M = 0`.
3. **`@f32` twins.** After every binary64 stratum, the `@f32` twin of each, in order. Each record's
   matrix is its binary64 namesake's rounded to nearest-even binary32, and its reference is the
   function at that matrix ([`0016`](./0016-f32-exact-strata-for-the-scalar-coefficient-ids.md)
   items 1–2). `theta:subnormal@f32` is its own binary32-subnormal decade, as 0016 makes it.
   `q:nonunit@f32` scales the same quaternions by `1 ± 2⁻¹⁹`, since binary32 rounds `2⁻⁴⁵` away.
   That is the largest power of two Sophus's binary32 `SO3(Matrix)` accepts, so the stratum keeps
   an oracle (amended after the first run, which read a plain rotation there).
   `matrix:exact@f32` is `matrix:exact` under its own name. No existing record moves.
4. **Subjects.** `helicoid` answers `so3_from_matrix` at `f32`. The Sophus image answers it with
   `Sophus::SO3<float>`. `tf_tree_math` and GTSAM stay binary64-only.

## Rationale

Backward error `B_M` would need the subject's `R(q̂)` against the input, and every routine is
floored at `d_M/u` (EA.9(c)). That floor is a few `u` here, so a one-ulp defect at `√½` would sit
inside it. A forward score against the polar factor has no floor at `d_M = 0`. On rounded inputs
every routine pays the same `O(d_M)` gap, so domination still compares routines. The gap is
largest on `q:nonunit`, a rotation scaled by `1 ± 2⁻⁴⁵`. There every Shepperd extraction scores
89–163 `u`, against at most 1.6 `u` elsewhere.

The ℝ⁴ angle is half the rotation angle, so a rotation error `Δ` reads `Δ/(2u)`. That is the same
scale as `so3_exp`'s sign-aligned quaternion score (EA.3).

## Measured

Protocol: the committed corpus after this PR. `cargo xtask conformance --subject helicoid --fn
so3_from_matrix [--precision f32]`, and `--oracle tf_tree_math | sophus | gtsam` over the same id
(Sophus also at `f32`). Maxima in `u` per stratum; "best" is the smallest oracle maximum.

| Stratum | `helicoid` f64 | best oracle f64 | `helicoid` f32 | best oracle f32 (Sophus) |
|---|---|---|---|---|
| `matrix:exact` | **0.7985** | 0.7071 | 0.2871 | 0.2871 |
| `theta:1e-8` | 0.9766 | 0.4911 | 2.0e-8 | 2.0e-8 |
| `theta:1e-3` | 1.4082 | 0.4935 | 0.9822 | 0.7816 |
| `theta:1e0` | 1.5671 | 1.4377 | 1.5926 | 1.2825 |
| `theta:pi-1e-12` | 1.5571 | 0.8839 | 1.3702 | 1.4096 |
| `theta:dense` | 1.3928 | 1.0337 | 1.3701 | 0.8248 |
| `q:w0` | 1.2894 | 0.5928 | 1.1223 | 1.0644 |
| `q:nonunit` | 89.41 | 162.35 | 11.087 (at `1 ± 2⁻¹⁹`) | 20.528 |
| all strata | 89.41 (1.567 without `q:nonunit`) | | 1.8342 (`theta:pi-1e-11@f32`) | |

`matrix:exact` at `f64` is the reported defect. The two `√½` both come out as `ℓ`, the double
below the nearest, so the angle is 0 and the drift is `(1 − 2ℓ²)/(2u) = 0.7985`. The oracles' worst is one component one ulp
off, an angle of `0.7071`. The `0079` kernel is dominated on 18 of 31 binary64 strata and 14 of
31 binary32 strata. The binary64 ones are every `theta:1e-k` from `1e-8` to `1e0`, `theta:pi-1e-5`,
`theta:pi-1e-8` to `1e-12`, `theta:dense`, `q:w0` and `matrix:exact`. The oracles divide by
`2√c_j` and do not normalize. `0079` traded that accuracy for speed against exact Shepperd. Against
the polar factor the trade costs domination.

## Consequences

- `so3_from_matrix` is scored, so the envelope judges it. No-regress reads its rows once the
  baseline is blessed. Until the kernel is fixed, domination fails on the 32 strata above.
- `metric::Rule::BackwardOnly` and the "not scored" envelope count are gone. A future backward
  error (`Log` near π, §11) comes back as its own rule.
- The corpus grows by 2642 records (1.2 MB), to 65 355 261 bytes, under
  [`0071`](./0071-the-corpus-budget-is-64-mib.md)'s 64 MiB.

## Implementation plan

1. The instrument: the strata, the metric, the subjects, `NUMERICS.md` §11 and `PHASE1.md` §4.
   Verified by `just corpus-check`, `metric::tests::an_axis_is_a_direction_in_four_dimensions_that_forgives_the_sign`,
   `test_so3.MatrixStrataTest` and the measurement above.
2. The kernel (owed). `SO3::from_matrix` must dominate on every `so3_from_matrix` stratum at both
   precisions, and `matrix:exact` must come back on the unit sphere. Verified by this instrument:
   `just conformance` plus the three oracles, read by `just envelope`. The sphere, not correct
   rounding: a correctly rounded quarter turn reads `fl(‖q‖²) = 1 + 2u`, and the 24 are pinned by
   unit tests, which a stratum maximum cannot do (amended by `0086`).

## Open questions

None.
