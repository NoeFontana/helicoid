# 0101: `svd3` is one-sided Jacobi on `eig3`'s rotation

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** the PRs of its *Implementation plan*

## Context

`PHASE2.md` §6 specifies `svd3(a: &Mat3<S>) -> (Mat3<S>, Vec3<S>, Mat3<S>)` as "McAdams et al.,
fixed Jacobi sweep count, signed so $\det U = \det V = +1$". It is the last routine of §6 without
code, and its id is the last row of `coverage::OWED` for Phase 2
([`0056`](./0056-the-routines-d7-does-not-reach.md) *Further work* 1). Two consumers wait on it:
`NUMERICS.md` §3.4's projection onto SO(3), and locus-calib's M3 metric, whose Wahba alignment
takes the nearest rotation of a 3×3 cross-covariance ([`0102`](./0102-the-nearest-rotation-takes-newton-where-it-is-well-posed.md)).

Three facts shape the decision:

- **McAdams is not backward stable.** It diagonalises $S = \mathrm{fl}(A^\top A)$ by Jacobi, sorts,
  QR-factors $B = A\hat V$ by Givens rotations and keeps $\mathrm{diag}(R)$. The columns of $B$ are
  orthogonal only to $u\lVert A\rVert^2$, so QR leaves $r_{23} \approx (b_2\cdot b_3)/\sigma_2$, and
  dropping it costs up to $\sqrt u\,\lVert A\rVert$ when $\sigma_2 \approx \sigma_3 \approx \sqrt u\lVert A\rVert$
  (`maths/svd3.md` SV.2). A scratch numpy harness (not committed; $2\times10^4$ matrices per
  family, a longdouble reference) reads its backward error at **$4.5\times10^7\,u$** on
  $\sigma = (1, 10^{-8}, \tfrac12 10^{-8})$ at `f64`, and $839\,u$ on $(1, 10^{-4}, \tfrac12 10^{-4})$
  at `f32`; one-sided Jacobi reads $4.5\,u$ and $4.4\,u$ there. Its vectors inherit $A^\top A$'s
  conditioning, $u\lVert A\rVert^2/\lvert\sigma_i^2 - \sigma_j^2\rvert$, which in the gap-weighted
  metric of §11 is $\lVert A\rVert/(\sigma_i + \sigma_j)$: $10^k$ for a bottom pair at
  $10^{-k}\lVert A\rVert$. Its speed comes from approximate Givens rotations and `rsqrt`, which D16
  rules out.
- **The ordering is unstated.** `0015` PH.2: §6 does not say $\lvert\sigma_3\rvert \le \sigma_2$,
  and `nearest_rotation` $= UV^\top$ needs it (`maths/so3.md` SO.13(c)): $\Sigma = \mathrm{diag}(1,
  1, -5)$ with $U = V = I$ gives $UV^\top = I$ at squared distance 36, against 20.
- **`eig3`'s rotation is already the one-sided rotation.** One-sided Jacobi rotates a pair of
  columns $(b_p, b_q)$ of $B = AV$ to annihilate $\gamma = b_p\cdot b_q$, which is `eig3`'s
  problem on the Gram block $\begin{bmatrix}\alpha & \gamma\\ \gamma & \beta\end{bmatrix}$,
  $\alpha = \lVert b_p\rVert^2$, $\beta = \lVert b_q\rVert^2$ (SV.1). `NUMERICS.md` §17's form
  never divides by the entry, has three divisions and two roots with one division fewer on its
  chain than the textbook $\zeta = (\beta - \alpha)/2\gamma$ form, and an exactly diagonal input
  differentiates through it ($\dot t = \dot w/\mathrm{den}$ at $w = 0$). The $\zeta$ form needs a
  small-argument arm at $\lvert 2\gamma\rvert < 2^{-26}\lvert\delta\rvert$ to keep $\zeta^2$
  finite: a typed switch point, which D5 forbids, for a cost this form does not pay.

**The exit pays, as it did for `eig3`** (`0084`). The same harness, one-sided Jacobi on that
rotation, the exit of decision 3 with its cap of twelve rotations:

| family (`f64` / `f32`) | mean rotations | max | max error, exit / 8 fixed sweeps (`u`) |
|---|---|---|---|
| Gaussian entries | 9.1 / 7.4 | 12 / 11 | 6.9 / 6.9 (f64) |
| near-rotation, $10^{-8}$ / $10^{-4}$ | 7.7 / 5.8 | 11 / 9 | 7.4 / 7.4 |
| $\sigma = (1, 10^{-4}, \tfrac12 10^{-4})$ | 6.5 / 6.0 | 8 / 8 | 4.8 / 4.8 |
| columns graded $(1, 10^{-5}, 10^{-10})$ | 3.8 / 3.0 | 6 / 4 | 3.6 / 3.6 |
| rank 1 | 6.5 / 6.5 | 8 / 12 | 4.4 / 4.4 |

Five fixed sweeps are fifteen rotations on every input. A test at $u$ rather than $4u$ stalls on
$\gamma$'s own rounding: 2.2% of the Gaussian family then runs to the cap. Preconditioning by a
Householder QR saved at most one rotation on generic input and read 11–15 `u` of backward error
against 5–8 `u`; it is not built.

## Decision

`NUMERICS.md` gains §18, which states decisions 2–6; the signature of §6 is unchanged.

1. **The contract.** $A = U\Sigma V^\top$ with $U, V \in \mathrm{SO}(3)$ to a few $u$,
   $\sigma_1 \ge \sigma_2 \ge \lvert\sigma_3\rvert$, and $\operatorname{sign}\sigma_3 =
   \operatorname{sign}\det A$ ($+0$ where $\det A = 0$). Orthogonality, $\det = +1$ and the ordering
   hold by construction for every finite input, whether or not the exit was reached:
   non-convergence shows only as backward error, which the corpus measures. With it,
   `PHASE2.md` §6 states PH.2's ordering, and $UV^\top$ is the nearest rotation where that is unique
   (SO.13(c), SV.4). Accuracy is **absolute**: a columnwise backward error of a few
   $u\lVert a_j\rVert$, so singular values within a few $u\lVert A\rVert$ and singular vectors
   within a few $u\lVert A\rVert/\mathrm{gap}_i$ radians (Wedin). The relative accuracy one-sided
   Jacobi has on $A = BD$, about $u\,\kappa(B)$ (Demmel–Veselić), is measured by the
   `svd:col-graded` strata and held by their committed rows, not promised (the owner,
   2026-10-18).
2. **The sweep.** $B = A$, $V = I$. Pairs in `eig3`'s order, $(0, 1), (1, 2), (2, 0), \dots$; each
   rotation computes $\alpha$, $\beta$, $\gamma$ from the current columns, $d = \beta - \alpha$,
   $w_2 = (2\gamma)^2$, and `eig3`'s rotation of §17 with $a_{pq} = \gamma$: $(c, s)$, the safe arm
   included ($c = 1$ where $d^2 + w_2$ is below the smallest normal). Then
   $(b_p, b_q) \leftarrow (c\,b_p - s\,b_q,\ s\,b_p + c\,b_q)$ and the same on $V$'s columns. The two
   routines share one private `angle`, so `eig3`'s bits do not move.
3. **The exit.** After each rotation, if every pair satisfies
   $\gamma_{ij}^2 \le 16u^2\max(\alpha_i, f)\max(\alpha_j, f)$, $f = u^2\max_k\alpha_k$ (the relative
   test at $4u$, floored at $u\lVert B\rVert$ so that rounding-level null columns, which no
   relative test can settle, count as orthogonal), the loop stops; at most twelve rotations,
   `eig3`'s bound. The mask decides through `Real::branch` and `Mask::decide`, as in `eig3`: a
   scalar leaves, lanes run on held. Then the **polish**: the rotations of the two pairs other
   than the last one rotated, both computed from the exit state, applied in turn; their values move
   by at most the entries, their `Dual` lanes annihilate the entries' lanes to first order (`0084`
   decision 4).
4. **The finish.**
   - Sort the columns by $\alpha$, descending, with the network $(0,1), (1,2), (0,1)$, an exchange
     where the second is strictly larger. Each exchange is $(b_p, b_q) \leftarrow (b_q, -b_p)$ and
     the same on $V$: a quarter turn, so $\det V$ stays $+1$.
   - $U$ by Gram–Schmidt on the sorted columns: $\sigma_1 = \lVert b_1\rVert$, $u_1 = b_1/\sigma_1$
     ($e_1$ where $\sigma_1 = 0$); $w = b_2 - (u_1\cdot b_2)u_1$, $\sigma_2 = \lVert w\rVert$,
     $u_2 = w/\sigma_2$ (where $\sigma_2 = 0$, the normalised cross product of $u_1$ with the axis
     of its smallest component); $u_3 = u_1\times u_2$; $\sigma_3 = u_3\cdot b_3$. The residuals
     $u_i\cdot b_j$ are the converged cosines, so the columnwise bound survives (SV.3).
   - $V$ re-orthonormalised the same way, without roots: $v_1$ and $w = v_2 - (v_1\cdot v_2)v_1$
     each by one Newton step on $1/\lVert v\rVert$ (`eig3`'s, every column within about $17u$ of
     unit length), $v_3 = v_1\times v_2$.
   - One masked pass fixes rounding-level inversions at ties: where $\sigma_1 < \sigma_2$,
     $(u_1, u_2) \leftarrow (u_2, -u_1)$, the same on $V$, the values exchanged; then where
     $\sigma_2 < \lvert\sigma_3\rvert$, with $s = \operatorname{copysign}(1, \sigma_3)$:
     $(u_2, u_3) \leftarrow (s\,u_3, -s\,u_2)$, $(v_2, v_3) \leftarrow (v_3, -v_2)$,
     $(\sigma_2, \sigma_3) \leftarrow (s\,\sigma_3, s\,\sigma_2)$; then $(0, 1)$ again. Each keeps
     $U\Sigma V^\top$ and both determinants (SV.5).
5. **Scaling.** `eig3`'s exact power-of-two ladder, on the largest $\lvert a_{ij}\rvert$, with the
   windows $[2^{-128}, 2^{128})$ at `f64` and $[2^{-16}, 2^{28})$ at `f32`. The rotation's radicand
   is of degree four in $A$; inside these windows it cannot overflow, and a radicand below the
   smallest normal belongs to a block below $2^{-255}\lVert A\rVert^2$ (`f64`), $2^{-31}\lVert A\rVert^2$
   (`f32`), below $u\lVert A\rVert^2$ (SV.6). The `f32` window reaches $2^{28}$ so that a Kabsch covariance in mm² takes no
   ladder. $\Sigma$ is scaled back.
6. **Non-finite input** iterates on the zero matrix and reports $\Sigma = (\mathrm{NaN}, \mathrm{NaN},
   \mathrm{NaN})$, $U = V = I$; never a plausible number, never `unit_or`'s fallback. The zero
   matrix gives $\Sigma = 0$, $U = V = I$. Nothing panics.
7. **`Dual`.** A result differentiates the rotations taken, the polish's included. Where
   $\sigma_i$ is simple, $\dot\sigma_i = u_i^\top\dot Av_i$ to the convergence of the exit; where
   $\sigma_i^2 \ne \sigma_j^2$ for every $j$, $\dot u_i$, $\dot v_i$ are first-order perturbation
   theory. At a repeated singular value the lanes are finite and in general no directional
   derivative (`eig3`'s statement, `0081`). `svd3_tests` asserts $\dot\sigma_i = u_i^\top\dot Av_i$,
   a free and exact oracle.
8. **The corpus id `svd3`** (`0056` *Further work* 1).
   - `in`: `A` (3×3, column-major). `out`: `U`, `sigma[3]`, `V`. A record part may hold more than
     one matrix: `shape` is then an object keyed by field (`{"U": [3, 3], "V": [3, 3]}`). A part
     with one matrix keeps today's `"shape": [r, c]`, so no committed record moves.
   - Reference: `mp.svd_r` of the stored `A` at 120 digits, rechecked at 150, put in decision 1's
     convention (where $\det V < 0$, $v_3$ and $\sigma_3$ negated; where $\det U < 0$, $u_3$ and
     $\sigma_3$), each column pair signed so the largest $\lvert u_i\rvert$ component (the first of
     equals) is positive. Singular values that are zero, decided by the rank of the stored `A` in
     rationals, are exactly $0$. Cross-check: $\lVert AV - U\Sigma\rVert < 10^{-100}\lVert A\rVert$,
     $\lVert U^\top U - I\rVert$ and $\lVert V^\top V - I\rVert < 10^{-100}$, both determinants $+1$,
     the ordering.
   - **Metric `Rule::Svd`** (`NUMERICS.md` §11), the largest of: $\sigma$'s forward error in
     $u\lVert\sigma\rVert_2$; the vectors, a joint sign per column,
     $\max_i\min_{s=\pm1}\max(\lVert\hat u_i - s\,u_i\rVert, \lVert\hat v_i - s\,v_i\rVert)\,
     \mathrm{gap}_i/(\lVert\sigma\rVert_2\,u)$, $\mathrm{gap}_i = \min_{j\ne i}\min(\lvert\lvert\sigma_i\rvert -
     \lvert\sigma_j\rvert\rvert, \lvert\sigma_i\rvert + \lvert\sigma_j\rvert)$ from the reference; the
     backward error $\lVert A - \hat U\hat\Sigma\hat V^\top\rVert_F/(\lVert A\rVert_F\,u)$, exact; the
     orthogonality of each of $\hat U$, $\hat V$, $\lVert X^\top X - I\rVert_F/u$, exact. An
     ordering violated, or a determinant below zero (exactly), scores $1/u$. Each is a ratio of
     exact integers rounded once, as `metric::forward_u`.
   - **Strata** (16 records each, 64 in the deep tier, unless stated; binary32 ranges in
     parentheses; every draw from its own stream):

     | Stratum | Draw |
     |---|---|
     | `svd:random` (32) | Gaussian entries |
     | `svd:gap-1e-k/bottom`, `/top`, k ∈ {0, 2, 4, 6, 8, 10, 12} ({0, 2, 4, 6}) | $\sigma$ with the two smallest (largest) at relative gap $10^{-k}$; $U$, $V$ Haar |
     | `svd:graded-1e-k`, k ∈ {2, 4, 8, 12} ({2, 4, 6}) | $\sigma = (1, 10^{-k}, \tfrac12 10^{-k})$: McAdams' failure |
     | `svd:col-graded` (32) | Gaussian $B$ times $\mathrm{diag}(1, 10^{-5}, 10^{-10})$ ($(1, 10^{-2}, 10^{-4})$) |
     | `svd:rank2`, `svd:rank1` | $\sigma = (1, \sigma_2, 0)$, $(1, 0, 0)$, planted |
     | `svd:zero` (4) | the zero matrix, signed zeros alternating |
     | `svd:rotation` | a Haar rotation, rounded |
     | `svd:near-rotation-1e-k`, k ∈ {4, 8, 12} ({2, 4, 6}) | a Haar rotation plus Gaussian entries of size $10^{-k}$: the $\sqrt u$ cluster |
     | `svd:near-reflection-1e-k`, the same k | the same times $\mathrm{diag}(1, 1, -1)$: Kabsch with $d = -1$ |
     | `svd:essential` | $\sigma = (s, s, 0)$ |
     | `svd:ordering` (32) | $U\,\mathrm{diag}(a, a, -b)V^\top$, $b > a$: PH.2's case |
     | `svd:scale-up`, `svd:scale-down` | `svd:random` times $10^{\pm70}$ ($10^{\pm12}$): both ladder arms |

     That is 580 binary64 and 468 binary32 records, about 1.2 MB of the 1.68 MB 0071 leaves; the
     deep tier holds 64 per stratum.
   - **Oracles**: nalgebra's `SVD` and Eigen's `JacobiSVD<Matrix3, ComputeFullU | ComputeFullV>` in
     the Sophus image (`0083`'s precedent), each sorted and put in decision 1's convention by the
     runner, at both precisions (`0077`).
   - **Bars**: domination of both on every stratum at both precisions; the committed rows of
     `committed_linalg.rs` reproduced bit for bit (`0056` decision 5).
   - **Self-test**: `U` and `V` rotated by `1e-9` about their first pair fires on `svd:random`;
     $\sigma_3$'s sign dropped fires on `svd:ordering`.

## Rationale

- **One-sided Jacobi, not McAdams** (*Context*): backward stable by construction, with vectors at
  the SVD's own conditioning rather than $A^\top A$'s. McAdams plus a 2×2 SVD of $R$'s trailing
  block would be backward stable too ($r_{12}, r_{13} \lesssim u\lVert A\rVert$), and stays a bench
  variant of the measurement in step 2, not the shipped program: it computes a symmetric
  eigendecomposition, a QR and a 2×2 SVD where this computes one iteration.
- **`eig3`'s rotation, not the $\zeta$ arms** (*Context*): no typed switch point, one division
  fewer on the chain, one tested formula for two routines. Its price is the degree-four radicand,
  which decision 5's windows absorb.
- **An exit and a polish, not a fixed count**: `0084`'s argument, measured again above. The fixed
  count's reason, that a value-masked skip would freeze `Dual` tangents, is what the polish
  answers.
- **The finish by Gram–Schmidt and a cross product**: $\det U = +1$ and the sign of $\sigma_3$ come
  from the construction, not from a determinant test, and $V$'s rotations already have
  $\det = +1$; the quarter-turn exchanges keep it.
- **Absolute accuracy promised, relative measured**: every consumer named (Kabsch, Wahba, the
  essential matrix, the projection onto SO(3)) needs the absolute bound; the relative one is held
  by the committed rows without becoming a rustdoc promise.

## Consequences

- `PHASE2.md` §6 states decision 1, its bars and strata, and §0.0's row; `NUMERICS.md` §11 gains
  the `svd3` metric, §12 its domain row, §13 Demmel–Veselić's one-sided reference (already listed),
  §18 the formulas; `PHASE1.md` §4.3 and §4.4 the id and strata; `API.md`'s surface row is
  unchanged.
- `maths/svd3.md` carries SV.1–SV.6, checked.
- `coverage::OWED` loses its Phase 2 row; `svd3` is a required id.
- `helicoid-linalg/README.md` drops "Not yet: `svd3`".
- The corpus reader and writer accept a keyed `shape`; every committed record is byte-identical.

## Implementation plan

1. This record, `maths/svd3.md` and [`0102`](./0102-the-nearest-rotation-takes-newton-where-it-is-well-posed.md), drafts.
2. `svd3` and its instrument: the keyed `shape`, `gen/svd.py`, the metric, the self-test, both
   oracles, the kernel, `svd3_tests`, the committed rows, the bench rows and this record's
   *Measured* — verified by `svd3_tests`, `the_unverified_routines_reproduce_their_committed_rows`,
   `just corpus-check`, `just conformance --self-test`, `just oracle-nalgebra`,
   `just oracle-sophus`, and domination at both precisions.

## Open questions

None. The owner answered both of the proposal's on 2026-10-18: relative accuracy is measured,
not promised (decision 1); `nearest_rotation` is in scope, in its own record (`0102`).

## Measured

Filled in by step 2.

## Further work

1. Carried column norms ($\alpha' = \alpha - t\gamma$), which would take one dot product off each
   rotation's chain: measured against the rows first, since a shrinking column cancels.
2. A `Dual`-free entry without the polish (`0084` *Further work* 2), if a consumer's budget asks.
