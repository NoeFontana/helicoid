# 0023: `eig3` departs from omnisac, and the limits it inherits

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** —

## Context

`0017` decision 4 ports `eig3` "the same way" as `solve_cubic`: omnisac's algorithm, generic over
`S: Real`. `PHASE2.md` §6 gives it "as is" with the signature `eig3(a) -> (Vec3<S>, Mat3<S>)`,
eigenvalues ascending, `det = +1`. omnisac's `eigendecomp_sym3_signed` returns
`Option<EigDecomp3>`, descending, with a third column of sign `±(v_max × v_mid)` (`det = -1` for
50.0% of 10^6 random symmetric matrices), and gates that return `None`: a non-finite entry,
`p <= 1e-12 max|a_ij|`, a vanishing cross product, a non-finite `r`. The specified signature has
no failure channel, so a literal port does not exist; the port that landed (step 3 of `0017`)
made the choices below. None is in a spec, so each is proposed here. Measurements are `f64`, from
scratch harnesses that are not committed (`eig3_tests` says which).

## Decision (proposed)

1. **Every input has an answer; the eigenvalues are not finite where they have none.** No
   tolerance is ported. A non-finite entry of the lower triangle, `p^3` overflowing, and `p^3`
   underflowing while the spread of `A - q I` exceeds the rounding of `q` give non-finite
   eigenvalues (`NaN`, or infinite on overflow), never a plausible number; below that spread
   `l = q` is right. Non-finite is the report, as for a `Dual` derivative at `r = ±1`. `eig3`
   therefore asserts nothing (`API.md` R6 lists it beside `chol` and `solve_cubic`, which report
   through a mask).
2. **The frame is built from cross products of unit vectors and is orthonormal to a few `u` for
   every input.** `v2` is the longest cross product of two rows of `A - l2 I`, else `e_z`;
   `w` is `v2 × e`, `e` the axis of its smallest component, and `x = v2 × w`; `v1` is the
   longest cross product for `l1` projected on `(w, x)`, else `w`; `v0 = v1 × v2`, else `x`.
   Every normalisation divides by the largest entry first. omnisac's Gram-Schmidt reaches an
   orthogonality error of `4e-7` over 10^6 matrices (`|V^T V - I|`; the port: `4.7 u`), and its
   `det` is `-1` half the time.
3. **`acos` and `pi` were `0022`'s reading**, shared as `pub(crate)` items of `cubic.rs` — a
   reading `0022` has since reversed: `Real` gains `acos` and `cos`, and `pi` is already a literal,
   so these call sites move with its step 1 and `eig3`'s bounds are re-measured once. `cos` is
   `sin_cos(x).1` and `x / 3` replaces omnisac's `x * (1/3)`.
4. **The limits are documented in `# Domain` and pinned, none fixed here.**
   - *A tie of the top pair.* Where `l2 - l1` is below about `2 sqrt(u) |A|`, equal eigenvalues
     included, `v2` is rounding noise or `e_z`, and `v1` and `v0` are built from it: the vector of
     the isolated `l0` has a residual of the order of `|A|` (`0.8 |A|` for `diag(5, 1, 5)`,
     `0.99 |A|` for `diag(100, 100, 1)`). Inherited: omnisac returns `None` for these axis-aligned
     ties (a cross product vanishes) but answers with the same wrong vectors otherwise (`0.63 |A|`
     at `g = 0` and `0.25 |A|` at `g = 1e-9` over 200 rotations of `(1, 5, 5 + g)`).
   - *The double root.* `acos` at `r = ±1` gives eigenvalues off by up to `sqrt(u) |A|`
     (`6e-9 |A|` on a rank-1 matrix), and the vectors of a close pair by about
     `0.3 u (|A| / gap)^2` radians; measured in `eig3_tests`.
   - *The magnitude range.* `p^3` must be normal: about `1e-100 < m < 1e100` (`f64`), `1e-12 < m <
     1e12` (`f32`); the vectors need `|A|^4` normal.

## Rationale

The alternatives are a mask beside the result (a signature change to `PHASE2.md` §6, as `0017`
made for `solve_cubic`) or the literal port (an `Option`, which §6 excludes). A non-finite eigenvalue keeps
the specified signature, and a caller tests it with `is_finite`.

## Consequences

- `PHASE2.md` §9 gains an `eig3` adapter (order, gates, `det`), written with this record.
- `NUMERICS.md` states no `eig3` formula (§13 only cites Kopp), so until this record is `ready`
  and the formulas move there, the rustdoc of `eig3` is their only statement.
- The `eig:*` corpus ids (`PHASE2.md` §6, `0006`) and a generator for the fixture and the
  double-double table owe a `PHASE1.md` step; until it lands the bounds in `eig3_tests` are a fit,
  not a verdict.

## Implementation plan

1. `eig3` as landed, documented and pinned (the `eig3` PR).
2. Whichever of questions 1 to 3 is answered: the remedy for the tie, then the `eig:*` ids.

## Open questions

1. **Accept decisions 1 and 2**, or take a mask beside the result?
2. **Anchoring.** The tie failure is a consequence of anchoring the frame on the largest
   eigenvalue, as omnisac does. Anchoring on the more isolated end of the spectrum (the end whose
   gap is larger, chosen by a mask on the two computed gaps) keeps the residual at the eigenvalue
   error: a `numpy` prototype without fallbacks measures `7e-9 |A|` at `g = 0` and
   `6.6e-9 |A|` at `g = 1e-9` on 200 rotations of `(1, 5, 5 + g)`, against `0.8 |A|` for the port.
   It is not Kopp's hybrid and needs no iteration, but it is a formula, so a `NUMERICS.md` edit and
   this record; the eigenvalue's `sqrt(u)` at a double root still needs the hybrid of `PHASE2.md`
   §6. Anchor by mask, or build the hybrid?
3. **The bars' `K`.** `PHASE2.md` §6 gives the bars as `K u |A|` and `K u |A| / gap` with no `K`,
   so "the bar holds" has no meaning; `eig3_tests` asserts a fit. Name a `K`, or let the bars be
   `0006`'s envelope against `nalgebra`'s `SymmetricEigen` (§8)?
4. **Scaling.** Dividing by the largest entry first would remove the magnitude cliff on both
   sides at the price of one rounding per entry (up to `u |A|` on the eigenvalues) and another
   departure from omnisac; or keep non-finite eigenvalues outside the range.
