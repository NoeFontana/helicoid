# 0101: `svd3` is one-sided Jacobi on `eig3`'s rotation

**Status:** ready
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
   - $U$ by Gram–Schmidt on the sorted columns, $N$ one Newton step on $1/\lVert v\rVert$,
     $v(3 - \lVert v\rVert^2)/2$ (`eig3`'s): $\sigma_1 = \lVert b_1\rVert$,
     $u_1 = N(b_1/\sigma_1)$ ($e_1$ where $\sigma_1 = 0$); $w = b_2 - (u_1\cdot b_2)u_1$,
     $\sigma_2 = \lVert w\rVert$, $u_2 = N(w/\sigma_2)$ (where $\sigma_2 = 0$, Gram–Schmidt of the
     axis of $u_1$'s smallest component); $u_3 = N(u_1\times u_2)$; $\sigma_3 = u_3\cdot b_3$. The
     residuals $u_i\cdot b_j$ are the converged cosines, so the columnwise bound survives (SV.3).
     Without the three $N$ the frames read 8–10 `u` and set the metric's maximum; with them, 4.5.
   - $V$ re-orthonormalised the same way, without roots: $v_1 = N(v_1)$,
     $v_2 = N(v_2 - (v_1\cdot v_2)v_1)$, $v_3 = N(v_1\times v_2)$. Recomputing $B = AV$ instead
     would lose the graded columns' accuracy.
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
     | `svd:gap-1e-k/bottom`, `/top`, k ∈ {1, 2, 4, 6, 8, 10, 12} ({1, 2, 4, 6}) | $\sigma$ with the two smallest (largest) at relative gap $10^{-k}$; $U$, $V$ Haar |
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

**Corpus.** Per stratum, the maximum in `u` of `Rule::Svd`, for three programs:
- this;
- nalgebra 0.35's `SVD`, which is McAdams: `svd_ordered3`, from $A^\top A$;
- Eigen 3.4's `JacobiSVD<Matrix3<T>>`, in the Sophus image.

The maximum is 6.03 `u` at `f64` (`gap-1e-2/top`) and 5.08 at `f32` (`random`). Domination holds on
all 61 strata. On each one, the better oracle is at least 2.06× this program's maximum at `f64` and
2.34× at `f32`. `svd:zero` ties at 0.

Eigen's smallest nonzero maximum is 8.03 `u`. nalgebra reads $4.2\times10^7$ at `graded-1e-8` and
772 at `graded-1e-4@f32`, the backward error measured by *Context*'s harness.

The `gap-1e-k` strata start at $k = 1$. Two nonnegative singular values differ by less than
$\lVert\sigma\rVert_2$ unless the third is 0, so no spectrum has a relative gap of $10^0$.

`svd:col-graded` is $A = B\,\mathrm{diag}(1, 10^{-5}, 10^{-10})$, or $(1, 10^{-2}, 10^{-4})$ at
`f32`. It scores 3.41 at `f64` and 3.57 at `f32` in the absolute metric. Its committed rows hold
whatever relative accuracy that shows (decision 1).

| stratum | this | nalgebra | Eigen |
|---|---|---|---|
| `svd:random` | 4.22 | 14.8 | 20.9 |
| `svd:gap-1e-1/bottom` | 5 | 16.4 | 16.2 |
| `svd:gap-1e-2/bottom` | 5.36 | 14 | 14.6 |
| `svd:gap-1e-4/bottom` | 4.15 | 11.6 | 15.1 |
| `svd:gap-1e-6/bottom` | 3.2 | 12.5 | 20.4 |
| `svd:gap-1e-8/bottom` | 4.26 | 12.7 | 13.3 |
| `svd:gap-1e-10/bottom` | 4.07 | 13.1 | 18.1 |
| `svd:gap-1e-12/bottom` | 3.4 | 14.2 | 22.4 |
| `svd:gap-1e-1/top` | 4.54 | 12.3 | 13.1 |
| `svd:gap-1e-2/top` | 6.03 | 14.2 | 20.9 |
| `svd:gap-1e-4/top` | 3.85 | 14.7 | 16.1 |
| `svd:gap-1e-6/top` | 4.14 | 13.2 | 17.5 |
| `svd:gap-1e-8/top` | 3.84 | 10.4 | 17.8 |
| `svd:gap-1e-10/top` | 3.57 | 11.8 | 12.9 |
| `svd:gap-1e-12/top` | 3.64 | 13.3 | 20.8 |
| `svd:graded-1e-2` | 4.46 | 141 | 12.8 |
| `svd:graded-1e-4` | 3.92 | 9.5e+03 | 13.7 |
| `svd:graded-1e-8` | 3.56 | 4.22e7 | 14 |
| `svd:graded-1e-12` | 4.27 | 4.26e+03 | 12.2 |
| `svd:col-graded` | 3.41 | 14.2 | 13.7 |
| `svd:rank2` | 3.83 | 11 | 16.4 |
| `svd:rank1` | 2.83 | 12.9 | 9.6 |
| `svd:zero` | 0 | 0 | 0 |
| `svd:rotation` | 3.36 | 9.97 | 8.03 |
| `svd:near-rotation-1e-4` | 3.5 | 11 | 15.5 |
| `svd:near-rotation-1e-8` | 4.14 | 12.5 | 17.4 |
| `svd:near-rotation-1e-12` | 4.96 | 12.4 | 13.7 |
| `svd:near-reflection-1e-4` | 3.76 | 14.3 | 23 |
| `svd:near-reflection-1e-8` | 3.54 | 13.7 | 17.1 |
| `svd:near-reflection-1e-12` | 3.46 | 12.2 | 14.4 |
| `svd:essential` | 4 | 8.26 | 15.4 |
| `svd:ordering` | 3.52 | 13.7 | 23 |
| `svd:scale-up` | 3.98 | 11.5 | 23.8 |
| `svd:scale-down` | 3.58 | 11.6 | 14.9 |
| `svd:random@f32` | 5.08 | 14.4 | 17.7 |
| `svd:gap-1e-1/bottom@f32` | 3.95 | 14.5 | 15.5 |
| `svd:gap-1e-2/bottom@f32` | 3.88 | 12.6 | 18.1 |
| `svd:gap-1e-4/bottom@f32` | 4.04 | 12.5 | 18.5 |
| `svd:gap-1e-6/bottom@f32` | 3.36 | 9.99 | 11.8 |
| `svd:gap-1e-1/top@f32` | 4.59 | 11.4 | 17.7 |
| `svd:gap-1e-2/top@f32` | 4.05 | 12.2 | 11.8 |
| `svd:gap-1e-4/top@f32` | 3.97 | 14 | 10.4 |
| `svd:gap-1e-6/top@f32` | 3.92 | 11.9 | 13.7 |
| `svd:graded-1e-2@f32` | 4.09 | 112 | 17 |
| `svd:graded-1e-4@f32` | 3.58 | 772 | 11 |
| `svd:graded-1e-6@f32` | 3.21 | 13 | 9.08 |
| `svd:col-graded@f32` | 3.57 | 83 | 15.9 |
| `svd:rank2@f32` | 3.28 | 9.85 | 14.5 |
| `svd:rank1@f32` | 3.54 | 10.2 | 9.01 |
| `svd:zero@f32` | 0 | 0 | 0 |
| `svd:rotation@f32` | 3.5 | 14.3 | 8.19 |
| `svd:near-rotation-1e-2@f32` | 4.49 | 12.3 | 14.8 |
| `svd:near-rotation-1e-4@f32` | 4.28 | 12.8 | 15.4 |
| `svd:near-rotation-1e-6@f32` | 3.76 | 14.1 | 10.8 |
| `svd:near-reflection-1e-2@f32` | 4.14 | 13.2 | 17.3 |
| `svd:near-reflection-1e-4@f32` | 3.86 | 15.4 | 14.3 |
| `svd:near-reflection-1e-6@f32` | 3.89 | 13.6 | 9.82 |
| `svd:essential@f32` | 3.6 | 10.3 | 9.12 |
| `svd:ordering@f32` | 3.65 | 11.8 | 17 |
| `svd:scale-up@f32` | 4.02 | 14.3 | 13.7 |
| `svd:scale-down@f32` | 3.85 | 12 | 11.5 |

**Large sample.** `measure_svd3`: $5\times10^4$ matrices per family and precision, scored against a
double-double reference of the stored matrix (max / p99.9 / mean).
- The overall maximum is 8.14 `u` at `f64` (`top-1e-0`) and 7.20 at `f32` (`random`).
- Every mean is between 1.6 and 2.6.
- Each proptest, over 10⁶ cases with two seeds, reads for the planted bounds: `back` 8.7, `sv`
  10.3, `ang` 6.4, and 4.5 on the frames. These are under the bounds `svd3_tests` records,
  `K = 21` and `FRAME = 10`.

| family | `f64` max | p99.9 | mean | `f32` max | p99.9 | mean |
|---|---|---|---|---|---|---|
| `random` | 7.02 | 6.04 | 2.60 | 7.20 | 5.30 | 2.41 |
| `col-graded` | 5.93 | 4.84 | 2.06 | 4.06 | 3.42 | 1.58 |
| `rank1` | 5.89 | 4.25 | 1.94 | 5.21 | 4.01 | 1.86 |
| `rank2` | 7.53 | 5.54 | 2.37 | 6.10 | 5.00 | 2.19 |
| `essential` | 5.57 | 4.64 | 2.13 | 4.82 | 3.95 | 1.88 |
| `ordering` | 6.48 | 4.91 | 2.04 | 5.66 | 4.46 | 1.89 |
| `rotation` | 5.46 | 4.44 | 2.15 | 4.12 | 3.64 | 1.81 |
| `bottom-1e-0` | 7.03 | 5.49 | 2.49 | 6.65 | 4.89 | 2.28 |
| `top-1e-0` | 8.14 | 5.58 | 2.43 | 6.10 | 5.07 | 2.24 |
| `bottom-1e-2` | 6.90 | 5.52 | 2.40 | 5.84 | 4.97 | 2.23 |
| `top-1e-2` | 6.62 | 5.37 | 2.48 | 5.54 | 4.52 | 2.28 |
| `bottom-1e-4` | 7.06 | 5.46 | 2.38 | 4.92 | 4.05 | 1.96 |
| `top-1e-4` | 6.75 | 4.96 | 2.42 | 6.05 | 4.73 | 2.27 |
| `bottom-1e-8` | 5.38 | 4.67 | 2.11 | 4.94 | 3.96 | 1.91 |
| `top-1e-8` | 5.56 | 4.71 | 2.24 | 4.37 | 3.71 | 1.85 |
| `bottom-1e-12` | 5.74 | 4.58 | 2.09 | 5.01 | 4.02 | 1.91 |
| `top-1e-12` | 5.68 | 4.64 | 2.18 | 4.36 | 3.78 | 1.84 |
| `graded-1e-4` | 6.27 | 5.10 | 2.17 | 5.03 | 3.93 | 1.86 |
| `graded-1e-8` | 5.83 | 4.43 | 1.98 | 5.41 | 3.96 | 1.87 |
| `near-rotation-1e-4` | 6.35 | 5.30 | 2.58 | 6.18 | 5.24 | 2.46 |
| `near-rotation-1e-8` | 6.45 | 5.41 | 2.55 | 4.40 | 3.64 | 1.85 |
| `near-reflection-1e-4` | 6.37 | 5.36 | 2.58 | 6.78 | 5.24 | 2.46 |
| `near-reflection-1e-8` | 6.32 | 5.36 | 2.55 | 4.57 | 3.66 | 1.85 |

**Cost.** In ns. This program's timings are criterion's estimates from
`cargo bench -p helicoid --bench linalg -- svd3`, run with `taskset -c 2` while other sessions held
the load at 5–9. The same four matrices were timed through two others:
- nalgebra's `SVD::new(m, true, true)`, in a paired scratch harness, interleaved, min of 31 runs of
  $2\times10^4$;
- Eigen's `JacobiSVD` with full `U` and `V`, built with `g++ -O2 -ffp-contract=off` in the Sophus
  image, min of 31.

The data-dependent exit costs nothing on mixed input: the four-call row is the sum of the single
rows.

| row | `f64` | nalgebra | Eigen | `f32` | nalgebra | Eigen |
|---|---|---|---|---|---|---|
| random | 416 | 565 | 662 | 351 | 451 | 456 |
| `gap-1e-8` | 260 | 505 | 664 | 189 | 312 | 408 |
| `near-rotation-1e-8` | 378 | 567 | 786 | 162 | 310 | 120 |
| rank 2 | 467 | 564 | 661 | 401 | 452 | 509 |
| the four in turn | 1523 | | | 1099 | | |

The scratch harness reads this program at 0.48–0.84 of nalgebra's time on every row. It is also
faster than Eigen on every row but one, `near-rotation-1e-8` at `f32`: 162 against 120. That row
was not looked into further.

## Further work

1. Carried column norms ($\alpha' = \alpha - t\gamma$), which would take one dot product off each
   rotation's chain: measured against the rows first, since a shrinking column cancels.
2. A `Dual`-free entry without the polish (`0084` *Further work* 2), if a consumer's budget asks.
