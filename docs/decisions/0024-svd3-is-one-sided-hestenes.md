# 0024: `svd3` is one-sided Hestenes, not McAdams

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

`PHASE2.md` §6 specifies `svd3(a: &Mat3<S>) -> (Mat3<S>, Vec3<S>, Mat3<S>)` as "McAdams et al.,
fixed Jacobi sweep count (no convergence loop), signed so `det U = det V = +1` and only `sigma_3`
may be negative (`nearest_rotation` is `U V^T` directly)". `0017` decision 4 adds that `eig3` and
`svd3` are ported from omnisac; omnisac has no `svd3`, so there is nothing to port and §6's
reference is the whole specification.

McAdams computes `V` from the symmetric eigenproblem of `A^T A`, then `Sigma` and `U` from a Givens
QR of `B = A V`. Two of its ingredients do not survive this workspace: its approximate Givens
quaternion (the `gamma = 3 + 2 sqrt 2` region test with a fixed `cos/sin(pi/8)` fallback) trades
accuracy for branchlessness that `S::select` already gives, and its reciprocal square root is a
hardware instruction, which D16 forbids (bit-identical output on x86_64, aarch64 and wasm32).

What forces the decision now is `A^T A`. `0015` PH.2 already records that §6 is silent on
`|sigma_3| <= min(sigma_1, sigma_2)`, which `nearest_rotation = U V^T` needs (`so3.md` SO.13(c):
for `Sigma = diag(1, 1, -5)` and `U = V = I`, `U V^T` is at squared distance 36 where the nearest
rotation is at 20), and lists `svd3` as blocked on it.

**Measured**, with the scratch harness and the 82-row mpmath fixture of the implementation PR
(`f64`, 4 sweeps, worst over three matrices per stratum), on `sigma = (1, s, s(1 + 10^-4))` — two
small singular values in a close pair under one large one, the family where squaring compresses a
relative gap by `sigma_1 / sigma_2`. Both routes use the *same* exactly-diagonalizing eigensolver,
so the difference is the formation of `A^T A` and nothing else:

```text
              worst angle of a right singular vector      worst sigma error, in u sigma_k
  s        one-sided        two-sided (A^T A)          one-sided      two-sided (A^T A)
 1e-3       4.3e-11            5.3e-8                    1.5e2            6.6e4
 1e-6        1.8e-9            1.5e-3                    1.7e5            2.9e11
```

The cost is the factor `sigma_1 / sigma_2` the theory predicts (`10^3`, `10^6`), and no
eigensolver recovers it: the information is gone once `C = A^T A` is rounded. Consumers reach
`svd3` through exactly this shape — a near-rank-deficient geometry matrix, an essential matrix, a
corrupted rotation to project (`NUMERICS.md` §3.4) — so the loss is not hypothetical.

## Decision

1. **`svd3` is a cyclic one-sided Hestenes Jacobi**, not McAdams. `PHASE2.md` §6 is amended: its
   reference becomes Hestenes' one-sided form with Rutishauser's Givens coefficients, `NUMERICS.md`
   §13 keeps McAdams and Kopp as cited alternatives. The signature, the fixed sweep count, the
   sign convention and the `mp.svd_r` reference of §6 are unchanged. `A^T A` is never formed.

2. **The rotation.** For the column pair `(j, k)` of `B` (`B = A` initially, `V = I`), with
   `alpha = b_j . b_j`, `beta = b_k . b_k`, `gamma = b_j . b_k`:

   $$\zeta = \frac{\alpha - \beta}{2\gamma},\quad
     t = \frac{\operatorname{sgn}\zeta}{|\zeta| + \sqrt{1 + \zeta^2}},\quad
     c = \frac{1}{\sqrt{1 + t^2}},\quad s = t\,c$$

   `t` is the root of `t^2 + 2 zeta t - 1 = 0` of smaller magnitude, so `|theta| <= pi/4`, and it
   is **paired with the update**

   $$b_j' = c\,b_j + s\,b_k, \qquad b_k' = c\,b_k - s\,b_j$$

   applied to the columns of `B` and of `V` alike (post-multiplication by
   `G = [[c, -s], [s, c]]`). The pairing is normative because the other one does not converge:
   with `b_j' = c b_j - s b_k` and `b_k' = s b_j + c b_k` the same `t` rotates by `-theta` and
   leaves the off-diagonal at `2 cos(2 theta) gamma`, measured at `1.0` after 8 sweeps against
   `4.8 u` for the pairing above. `sgn(+0) = +1` (`copysign`), so `alpha = beta` gives the
   `45`-degree rotation the isotropic case needs.

   Only `gamma = +-0` is guarded, by the safe-argument pattern: there `zeta` would be `0/0` for
   `alpha = beta`, and `(c, s) = (1, 0)` is right because the pair is already orthogonal. **No
   negligibility threshold is typed**: a tiny `gamma` yields a correct small `t`, or the correct
   `45`-degree rotation when `alpha = beta`, so the classic threshold-skip would only change
   results, never rescue them. `0017` decision 3 therefore has nothing to apply to here.

3. **Four sweeps** of the cycle `(0,1), (0,2), (1,2)`, no convergence loop, fixed at every
   precision. Measured over 20507 matrices in 17 strata (`f64`), worst `max_{j<k} |b_j . b_k| /
   sigma_1^2`: 3 sweeps leaves 1342 matrices above `100 u` and reaches `4.0e10 u`; 4 sweeps reaches
   `6.7 u` with none above `100 u`; 5, 6 and 8 sweeps do not improve it and cost orthogonality
   (`|det V - 1|` grows from `26 u` to `40 u`). The count is an algorithmic parameter, not a switch
   point, so `0004` does not send it to `coeffs::generated`; it is a private constant and the
   measurement is a committed test, not a scratch harness (`eig3`'s is not reproducible).

4. **The gauge is constructed, never repaired.** No determinant is evaluated and no sign is
   tested against one.
   - Sort the columns of `B` and of `V` together, descending by `|b_i|`, with a three-step
     compare-exchange network; track the permutation's parity as a scalar `+-1`. `A = sum_i b_i
     v_i^T` is permutation-invariant, so the sort needs no compensation of its own. This is what
     makes `sigma_1 >= sigma_2 >= |sigma_3|` hold, which `0015` PH.2 asks §6 to state.
   - `u_0 = unit(b_0)`, `u_1 = unit(b_1 - u_0 (u_0 . b_1))`, and `u_2 = u_0 x u_1`. `det U = +1`
     by construction. `V` is a product of plane rotations, so `det V = +1` by construction.
   - `sigma_3 = parity * copysign(|b_2|, u_2 . b_2)` and the third column of `V` is multiplied by
     `parity`. The first factor puts the reflection of a `det A < 0` input in `sigma_3`, where §6
     wants it; the second keeps `det V = +1` after an odd permutation, and multiplying `sigma_3`
     by the same `parity` leaves `sigma_3 v_3^T` unchanged.
   - **Both frames are built the same way**: `v_0, v_1` normalized and orthogonalized as `u_0,
     u_1` are, and `v_2 = v_0 x v_1`. `V` is a product of plane rotations, so `v_0 x v_1` is its
     own third column times the sort's parity — the cross product therefore *is* the parity
     correction, and `sigma_3` carries the matching factor. Building `V` from `unit_or` rather
     than scaling its accumulated columns is what keeps it finite: a non-finite entry poisons
     every rotation, and a scaled NaN is still NaN.
   - Degeneracy is `eig3`'s `unit_or`/`any_orthogonal`, made `pub(crate)`: each normalisation
     divides by the largest entry first and falls back to a fixed frame where a column is zero or
     not finite. No threshold decides the rank. Where `b_2` is zero the sign of `sigma_3` is the
     sign of a zero and carries no meaning.

5. **Every input has an answer and nothing is asserted** (`API.md` R6, as `chol`, `solve_cubic`
   and `eig3`). A non-finite entry gives non-finite singular values and the fallback frame.

   **What holds for every input** is that `U` and `V` are finite, their columns normalized, and
   `sigma_1 >= sigma_2 >= |sigma_3|`: each column is a normalization or a cross product of two of
   them, so no input can break it. **Orthonormality and `det = +1` to the recorded bound hold
   inside the domain**, which is each column's squared norm a normal number (`|A|^2` for the
   singular values, `|A|^4` for the vectors, as `eig3`). Outside it a column norm overflows or
   goes subnormal, `zeta` saturates, the rotation that would have separated two nearly parallel
   columns is skipped, and the projection that builds the second column cancels: measured, the
   orthogonality of `U` reaches about `sqrt(u)`, and `svd3_tests` pins two such rows. An earlier
   draft of this record claimed a few `u` for *every* input; the proptest over arbitrary bit
   patterns disproved it, and the claim is now range-qualified.

   Scaling the input by its largest entry first would remove the cliff at both ends for one
   rounding per entry. That is `0023` (draft) question 4, asked there for `eig3` and identical
   here; this record keeps the documented range instead and does not answer it. Nothing in the
   *Decision* depends on the answer: scaling would be additive and would only narrow the
   documented domain's exclusions.

## Rationale

McAdams is the right algorithm for the problem it was written for: thousands of deformation
gradients per frame, `det > 0`, entries near 1, accuracy spent freely for branchlessness. This
crate's callers are Lie-group consumers — noisy geometry, gauge transformations, rank
degeneracies — where the `sigma_1 / sigma_2` loss measured above lands exactly on the inputs that
matter, and where branchlessness is already paid for by `S::select` rather than by an approximate
rotation. The one-sided form also needs no second algorithm for `Sigma` and `U`: the columns of
`B = A V` *are* the scaled left singular vectors, so reconstruction is structural.

Alternatives that lost: McAdams with exact rotations (keeps the `A^T A` loss, which is the whole
point); McAdams verbatim (that, plus a hardware `rsqrt` D16 forbids and a fixed-count residual
larger than the exact rotations'); a Jacobi eigendecomposition through the existing `eig3` (the
same squaring, plus `eig3`'s own unreliable frame at a tie of the top pair, `0023` (draft)).

Deciding the ordering rather than leaving it free (`0015` PH.2's second option, "have
`nearest_rotation` flip the axis of the smallest `|sigma|`") keeps the flip out of every consumer:
`nearest_rotation` stays `U V^T`, as §6 promises.

## Consequences

- `PHASE2.md` §6 is amended (the reference, the ordering, the pairing, the sweep count) and its
  §0.0 row records what is measured. `NUMERICS.md` §13 keeps its McAdams citation, now as an
  alternative; §3.4's "`svd3` first" route is unaffected.
- `0015` PH.2's `svd3` row is answered: the ordering is stated and `svd:*` strata are named below.
  Its other rows stay open, and `0015` stays `draft`.
- `eig3::unit_or` and `eig3::any_orthogonal` become `pub(crate)` and gain a second caller. They are
  not public API.
- The `svd:*` corpus ids (`PHASE1.md` §4.3, `mp.svd_r`) are owed a `PHASE1.md` step, as `eig:*` and
  the cubic's ids are. Until then the strata are the committed mpmath fixture and the planted
  proptest, the bars are a measured fit, and `svd3` ships in no release (`0006`).
- `svd3` is the `nearest_rotation` half of `NUMERICS.md` §3.4; the other half is `SO3::from_matrix`
  (Phase 3, `0015` NU.8), which this record does not touch.

## Implementation plan

1. This record, the `PHASE2.md` §6 amendment and the §0.0 row.
2. `svd3` and `svd3_tests` — verified by: the 82-row mpmath fixture (`mp.svd_r`, 120 digits,
   binary32-exact entries so both precisions decompose the same matrix) at `f64`, `f32` and `Dual`
   against the stated error model (`K = 24` for the singular values and angles, 32 for
   orthogonality, `det` and the reconstruction; worst `value / bound` 0.32 over the fixture and
   0.67 over 10^6 planted cases); a planted-spectrum proptest; the invariants over random,
   degenerate, out-of-range and non-finite inputs; the sweep count, asserting that four converge
   and three do not, **in the repository rather than in a scratch harness**; the relative accuracy
   of the column-graded family; the `Dual` value path bit-equal and its derivative against
   `d sigma_k = u_k^T dA v_k`; two lanes.
3. The `svd:*` corpus ids and the conformance subject, with `eig:*`, when `PHASE1.md`'s generator
   step lands. Strata: `svd:random`, `svd:graded-1e-k` (`sigma = 1, 10^-k, 10^-2k`),
   `svd:col-graded-1e-k` (`A = B diag(d)`, the relatively-accurate family),
   `svd:small-pair-1e-k` (the two-small-close family of *Context*), `svd:tie`, `svd:isotropic`,
   `svd:rank{2,1,0}`, `svd:rotation`, `svd:reflection`, `svd:scaled-2e{+-k}`, `svd:parallel`.

## Open questions

None.
