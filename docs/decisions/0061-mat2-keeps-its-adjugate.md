# 0061: `Mat2` keeps its adjugate

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #118, #121

## Context

`PHASE2.md` §0 promises "`Mat2`/`Mat3` adjugates" and its §0.0 says "`Mat2` adjugate pending".
`API.md` §2 lists only `Mat3::inverse_adj`. `0015` (draft) PH.2 reads the gap and recommends
"drop [`Mat2` from §0]; add with a consumer". `0015` decided nothing ("Decision: None yet"), so
PH.2 is a recommendation and this record answers it.

Two consumers now name it:

- **S².** `PHASE5.md` §1.4 types its chart Jacobian as `Mat2<S>`. A `Jac` needs `inverse`, and
  CH.10's `lj = A(δ)⁻¹ Z` is its inverse, so `Mat2` needs one.
- **locus-tag.** It inverts 2×2 matrices in nine places: the structure-tensor adjugate in
  `finalize_corner_covariance`, three `Matrix2::try_inverse` calls of a corner covariance, the
  `fuse_corner` symmetric inverse, two Newton steps and four Cramer line intersections. They use
  five unrelated singularity thresholds: `1e-4`, `1e-6`, `1e-12`, `ε²` and `1e-18`. The threshold
  is the consumer's decision, but it should apply to a stated `det`, not to an inverse that has
  already divided by it.

A 2×2 adjugate is exact: `adj [[a, b], [c, d]] = [[d, −b], [−c, a]]` is a permutation with sign
changes. All the rounding is in `det = a d − b c`, two products and a subtraction. D16 forbids the
fused multiply-add that Kahan's algorithm would use, so `det` has an absolute error of at most
about `2u (|a d| + |b c|)`. That is a relative error of `u` times the cancellation ratio
`(|a d| + |b c|)/|a d − b c|`, which is the matrix's conditioning and not the algorithm's.

## Decision

1. **`Mat2::inverse_adj(&self) -> (Self, S)`** in `crates/helicoid-linalg/src/matrix.rs` returns
   `(adj(self)/det, det)` with `det = m00 m11 − m01 m10`, in that order and unfused. Each entry is
   one division by `det`, with no reciprocal, as `Mat3::inverse_adj`. Its rustdoc states that the
   adjugate is exact and that the only rounding is `det`'s and the four divisions.
2. **The domain is `Mat3`'s, at the 2×2 degree.** It is copied into `# Domain`: the caller decides
   what `det` means; at `det = 0` the entries are non-finite and nothing is asserted, since a
   singular matrix is legal input. `det` is quadratic and unscaled, so entries of magnitude `m`
   need about `1e-154 < m < 1e154` (`f64`) and `1e-19 < m < 1e19` (`f32`). Beyond that, `det`
   overflows or underflows as in `Mat3`'s text.
3. **Corpus id `mat2_inverse_adj`**, binary64 and binary32, 6 records per stratum:
   - `cond:1`, `cond:1e4` and `cond:1e8`. `M = s · R(α) diag(1, 1/κ) R(β)ᵀ`, with `α` and `β`
     uniform, `s = 10^U(−3, 3)` and `κ` the stratum's, rounded to binary64.
     `cond:1e8` is binary64 only: `κ = 1e8` exceeds `1/u` at binary32 (about `1.7e7`), so those
     matrices are singular in that precision and `det` rounds to `0` (4 of 6 records when
     generated), which is `cond:singular`'s question, not this stratum's. The reference is the
     inverse of the **rounded** entries at 60 digits, so the stratum names the draw's conditioning
     and the reference is exact for the input actually given.
   - `cond:singular`. `M = [[x, 2^j x], [y, 2^j y]]` with `x`, `y` drawn as above and
     `j ∈ {−3, …, 3}`. The two products are the same power of two times `x y`, so they round
     equally and `det` is exactly `0` in both precisions. The reference is `det = 0` with the
     inverse undefined.
   - The metric is `NUMERICS.md` §11's forward error, Frobenius, on `adj/det`, plus the relative
     error of `det`. On `cond:singular`, `det` is scored absolutely, as `|det̂| / (‖M‖_F² u)`, and
     the inverse is **expected** to be non-finite: there a non-finite inverse is the answer, not
     §11's failing non-finite count, and a finite one scores `1/u`. This is a new rule in
     `xtask/src/conformance/metric.rs`, and §11 gains a bullet for it.
   - The oracle is nalgebra's `Matrix2::try_inverse`, added to `runners/nalgebra` (`0056`'s
     runner), with `None` read as the singular answer. That is what locus-tag calls today.
4. **`PHASE2.md` §0.0** changes "`Mat2` adjugate pending" to the id and this record. `API.md` §2's
   row becomes `Mat2::inverse_adj`, `Mat3::inverse_adj`. `PHASE1.md` §4.3 gains the definition.
   `0015` stays a draft and is cited as such.

## Rationale

- **Keep `Mat2` rather than drop it.** PH.2's condition was a consumer, and there are two. Without
  it, S²'s `Jac` would need a private 2×2 inverse, which is a second spelling of one (`API.md` §6
  item 5), and locus-tag keeps nine hand-written ones.
- **The pair, not an `Option`.** `Mat3`'s contract returns `det` and lets the caller choose. An
  `Option` would put a threshold, the thing locus-tag has five of, into `helicoid`, where no
  single value is right. It would also be a boolean decided from a float (`API.md` R4).
- **A singular stratum built to be exactly singular.** A drawn matrix with `κ = ∞` rounds to a
  nonsingular one, so a "singular" draw would score the conditioning of rounding noise. Equal
  rounding of the two products makes `det = 0` a property of the input.
- **Strata by conditioning.** For this routine conditioning is the only axis of difficulty. The
  scale draw covers magnitude inside decision 2's domain.

## Consequences

- `Mat2` stays in `helicoid-linalg`'s surface, with an inverse whose contract is `Mat3`'s.
- `NUMERICS.md` §11 gains the expected-non-finite rule, which applies only to strata that declare
  it. Today that is `cond:singular` alone.
- A `Mat2` `Jac` impl for S²'s tangent is `PHASE5.md` §2's and not this record's.

## Implementation plan

1. This record, `PHASE2.md` §0.0, `API.md` §2 and `NUMERICS.md` §11. Verified by `just lint` and
   `just doc`.
2. `Mat2::inverse_adj` and its tests in `linalg_tests.rs`, mirroring `Mat3`'s: `adj2_exact`
   (integer entries, exact to the bit), `adj2_residual` (`M · inv ≈ I` scaled by the cancellation
   ratio), `adj2_identity`, the power-of-two range test, and a `Dual` lane check. Then the
   generator, the metric rule, the `Linalg` subject in `xtask/src/shipped.rs`, the nalgebra runner's
   answer, `coverage.rs`'s `required()` and `PHASE1.md` §4.3. Verified by `just test`,
   `just corpus-check`, `just conformance`, `just oracle-nalgebra`, `just envelope`, `just no-std`
   and `just msrv`.

## Open questions

None.

## Further work

- `Mat3::inverse_adj` has unit tests but no corpus id. The same strata, with a rank-2 product
  construction for `cond:singular`, would give it one.
