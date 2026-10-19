# 0102: the nearest rotation takes Newton where it is well posed

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** the PRs of its *Implementation plan*

## Context

`NUMERICS.md` §3.4 projects an arbitrary 3×3 onto SO(3) as "`svd3` first ($U\,\mathrm{diag}(1, 1,
\det UV^\top)\,V^\top$), then this", and `SO3::from_matrix`'s rustdoc tells a caller to project
with its own SVD until `svd3` exists. `0015` PH.2 and `maths/so3.md` SO.13(c) name the routine
`nearest_rotation` $= UV^\top$, with `svd3`'s signed convention. Consumers: Kabsch and Wahba
alignment (locus-calib's M3 metric), a drifted rotation matrix, a homography's or an essential
matrix's rotation.

**`svd3` alone is the wrong program near SO(3)**, which is where most inputs are:

- Its worst case is a cluster of singular values, and a near-rotation is one: $\sigma_i = 1 +
  O(\varepsilon)$. `svd3` pays its most rotations exactly on the input that is easiest for a polar
  iteration.
- **Its tangents there are wrong.** At an exact rotation $\sigma = (1, 1, 1)$, so `svd3`'s first
  rotation's $d$ and $\gamma$ are both rounding errors: the angle is arbitrary and its `Dual` lane
  is $O(1/u)$. $\dot U$ and $\dot V$ cancel in $UV^\top$ only to $O(1)$: the tangent of $UV^\top$ is
  $O(1)$ wrong, while the polar factor is perfectly conditioned there ($\kappa = 2/(\sigma_2 +
  \sigma_3) = 1$).

**Newton's polar iteration has neither problem** (PO.1): $X \leftarrow \tfrac12(\mu X +
X^{-\top}/\mu)$ converges quadratically to the orthogonal polar factor of a nonsingular $X$, which
is the nearest rotation when $\det X > 0$, and it is rational, so `Dual` differentiates it
cleanly. $X^{-\top} = \mathrm{cof}(X)/\det X$ needs no root. A scratch numpy harness (not
committed; $2\times10^4$ matrices per family, a longdouble reference; error is §11's weighted
forward error, decision 7):

| family | Newton: max error, orthogonality (`u`), steps mean / max | `svd3`'s $UV^\top$ |
|---|---|---|
| near-rotation $10^{-2}$ | 2.07, 6.15; 4.0 / 5 | 5.05, 14.3 |
| near-rotation $10^{-8}$ | 1.95, 8.28; 2.0 / 3 | 4.66, 14.5 |
| near-rotation $10^{-16}$ | 1.78, 6.99; 1.0 / 2 | 5.04, 15.2 |
| $cR$, $c$ log-uniform in $e^{\pm10}$ | 1.62, 6.21; 2.0 / 2 | 5.53, 15.7 |
| $\det > 0$, $\sigma$ half-normal | 2.08, 6.30; 6.0 / 7 | 3.92, 14.3 |
| $\sigma = (1, \tfrac12, 10^{-8})$ | 1.02, 6.20; 7.0 / 7 | 1.89, 16.0 |
| `f32`, near-rotation $10^{-4}$ | 1.72, 6.11; 2.0 / 2 | 4.60, 15.9 |

Complex-step derivatives of the same iteration (a stand-in for `Dual`) against PO.3's closed form:
$7.9\times10^{-15}$ relative after one step at exact rotations, $1.5\times10^{-15}$ after two at
$10^{-8}$, and $4\times10^{-13}$ after six on the $\det > 0$ family, one step behind the values.

**Where Newton fails.** It computes $\det X$ to an absolute $u\lVert X\rVert^3$, so its relative
error is $u/\rho$, $\rho = \det M/\lVert M\rVert_F^3$. Measured on $\sigma = (1, a, b)$: at `f64` it
reads at most 1.5 `u` down to $\rho = 10^{-16}$ and $9\times10^{11}$ at $10^{-18}$; at `f32` 1.4 `u`
at $3.5\times10^{-7}$, $1.2\times10^{7}$ (and NaN on 13% of the family) at $3.5\times10^{-9}$.
`svd3`'s $UV^\top$ reads at most 3 `u` on every one of them.

## Decision

`NUMERICS.md` gains §19, which states decisions 2–5.

1. **The surface.** `helicoid-linalg` gains, in `polar.rs`:

   ```rust
   pub fn nearest_rotation<S: Real>(m: &Mat3<S>) -> Mat3<S>;
   ```

   It returns $R \in \mathrm{SO}(3)$, orthonormal to a few $u$ with $\det R > 0$, for every finite
   input: the Frobenius-nearest rotation $U\,\mathrm{diag}(1, 1, \det UV^\top)V^\top$ where it is
   unique ($\sigma_2 + d\sigma_3 > 0$, SO.13(c)), one of the nearest where it is not. It is
   matrix analysis, not a group: no quaternion, no `SO3`. `SO3::from_matrix` is unchanged; §3.4
   reads "`nearest_rotation` first, then this".
2. **Scaling and the arm test.** The input is scaled by `svd3`'s exact power-of-two ladder and
   windows (`0101` decision 5); the rotation needs no scaling back. With $n = \lVert M\rVert_F$ and
   $\det M$ from the cofactors, the arms are:
   - **Newton**, where $\det M > \tau n^3$;
   - **reflect**, where $\det M < -\tau n^3$;
   - **`svd3`**, between, the zero matrix included,

   with $\tau = 2^{-26}$ (`f64`), $2^{-12}$ (`f32`), about $\sqrt u$. Both sides of $\tau$ are
   accurate by a factor of about $1/\sqrt u$ (*Context*). It is a choice between two complete
   programs, not a switch point between a series and a closed form, and moving it moves no bound.
3. **Newton**, from $X_0 = M$:
   $C = \mathrm{cof}(X) = [x_2\times x_3,\ x_3\times x_1,\ x_1\times x_2]$,
   $\delta = x_1\cdot c_1$, $\mu = \sqrt{\lVert C\rVert_F/(\lvert\delta\rvert\,\lVert X\rVert_F)}$
   (Higham's Frobenius scaling), $X \leftarrow \tfrac12(\mu X + C/(\mu\delta))$, until
   $\lVert X_{k+1} - X_k\rVert_F^2 \le 48u^2$ ($4u\sqrt3$) or eight steps. The step that
   detects convergence is applied: it is the one the tangents need, one behind the values. The
   exit decides through `Real::branch` and `Mask::decide`, as `eig3`'s.
4. **Reflect.** Newton from $X_0 = M$ converges to the orthogonal polar factor $Q$, $\det Q = -1$.
   $H = \tfrac12(Q^\top M + M^\top Q)$ has the eigenvalues $(\sigma_1, \sigma_2, \lvert\sigma_3\rvert)$;
   $\lambda$, its smallest, from `eig3`, then one Newton step on $p(\lambda) = \det(H - \lambda I)$,
   $\lambda \leftarrow \lambda + p/\operatorname{tr}\mathrm{adj}(H - \lambda I)$;
   $P = \mathrm{adj}(H - \lambda I)/\operatorname{tr}\mathrm{adj}(H - \lambda I)$, the projector onto
   $\lambda$'s eigenvector; $R = Q(I - 2P)$ (PO.2), then the polish of decision 5. Every step
   after `eig3` is rational in $M$ and $\lambda$, and the refined $\lambda$ has the root's derivative whatever
   `eig3`'s lanes, so the tangent is exact wherever $\sigma_3$ is simple, which is where the nearest
   rotation of a $\det < 0$ matrix is unique, $\sigma_1 = \sigma_2$ included. The trace
   $t = \operatorname{tr}\mathrm{adj}(H - \lambda I)$ is the product of the gaps from $\lambda$ to
   the other two eigenvalues; where $t \le e(g + e)$, $e = 16u\lVert M\rVert_F$, $g$ the largest
   of `eig3`'s eigenvalues less the smallest, the gap $\sigma_2 - \lvert\sigma_3\rvert$ is at
   rounding and the nearest rotation is not unique: decision 5's. The review replaced
   $t \le 16u\lVert M\rVert_F^2$, which compares a product of two gaps with a squared norm: it
   sent $-\mathrm{diag}(1, 1, 1 - \delta)$, $\delta$ below $1.7\times10^{-3}$ at `f32`, to
   `svd3`, whose tangent at the $\sigma_1 = \sigma_2$ tie is no derivative. `eig3`'s eigenvalues
   come without its frame, bit for bit.
5. **`svd3`.** $R = UV^\top$, from `svd3`'s signed convention (`0101` decision 1), then the
   **polish**: one Newton step, which squares the error of the product of two frames (7.8–8.6 `u`
   without it), and a second where the first was longer than $\sqrt u$, $\lVert\Delta\rVert_F^2 > u$.
   A step leaves about the square of its length. The second step serves the reflect arm, whose $P$
   loses accuracy as $\sigma_2$ and $\lvert\sigma_3\rvert$ close: with
   $\sigma = (1.95, 0.2513, -0.2495)$ one step left 19.5 `u` on the frame at `f32`. The large
   sample's figures did not move.
6. **Non-finite input** gives a matrix of NaN, never a plausible rotation; the zero matrix gives
   $I$ (the `svd3` arm, $U = V = I$). Nothing panics.
7. **`Dual`.** The Newton and reflect arms differentiate a rational iteration and are exact to the
   convergence of the exit wherever the nearest rotation is unique, ties among the singular values
   included. The `svd3` arm, below $\tau$, has `svd3`'s tangents: exact where the singular values
   are simple, finite and in general no directional derivative at a tie. `polar_tests` asserts the
   tangents against PO.3's $\dot R = U\Omega V^\top$, $\Omega_{ij} = (C_{ij} - C_{ji})/(\tilde{\sigma}_i +
   \tilde{\sigma}_j)$, $C = U^\top\dot MV$, $\tilde{\sigma} = (\sigma_1, \sigma_2, d\sigma_3)$.
8. **The corpus id `nearest_rotation`.**
   - `in`: `M` (3×3). `out`: `R` (3×3).
   - Reference: $U\,\mathrm{diag}(1, 1, d)V^\top$ from `mp.svd_r` of the stored `M` at 120 digits,
     rechecked at 150, $d$ the exact sign of $\det M$ in rationals. Cross-check:
     $\lVert R^\top R - I\rVert < 10^{-100}$, $\det R = 1$ to 100 digits, and no better rotation
     among 50 random ones.
   - **Metric `Rule::Polar`** (`NUMERICS.md` §11), the largest of: the forward error weighted by
     the conditioning,
     $\lVert\hat R - R\rVert_F\,(\sigma_2 + d\sigma_3)/(2\lVert M\rVert_F\,u)$, weight $0$ where
     the reference is not unique (PO.3: $\lVert\dot R\rVert \le 2\lVert\dot M\rVert/(\sigma_2 +
     d\sigma_3)$); $\lVert\hat R^\top\hat R - I\rVert_F/u$, exact; and $1/u$ where $\det\hat R \le 0$.
     The weight is $\operatorname{tr}H - \lambda_{\max}(H)$, $H = \mathrm{sym}(R^\top M)$, in integers
     from the 30-digit $R$, rounded once: a binary64 eigenvalue carries $u\lVert M\rVert$ absolute,
     all of the weight on `non-unique` (this record's review).
   - **Strata** (16 each, 64 deep, unless stated):

     | Stratum | Draw |
     |---|---|
     | `polar:near-1e-k`, k ∈ {2, 4, 8, 12, 16} ({2, 4, 6}) | a Haar rotation plus Gaussian entries of size $10^{-k}$ |
     | `polar:exact` (24) | the signed permutation rotations of `matrix:exact` (`0085`) |
     | `polar:scaled` | $cR$, $c$ log-uniform in $[10^{-3}, 10^3]$ |
     | `polar:random` (32) | Gaussian entries: both signs of $\det$ |
     | `polar:reflection` | $U\,\mathrm{diag}(1, a, -b)V^\top$, $b < a/2$: the reflect arm |
     | `polar:cond-1e-k`, k ∈ {4, 8, 12, 16} ({2, 4, 6, 8}) | $\sigma = (1, \tfrac12, 10^{-k})$: both sides of $\tau$ |
     | `polar:rank2` | $\sigma = (1, \tfrac12, 0)$ |
     | `polar:non-unique` | $U\,\mathrm{diag}(1, a, -a)V^\top$: weight 0, orthogonality and $\det$ scored |
     | `polar:scale-up`, `polar:scale-down` | `polar:random` times $10^{\pm70}$ ($10^{\pm12}$) |

     That is 296 binary64 and 264 binary32 records, about 0.3 MB.
   - **Oracles**: nalgebra's `Rotation3::from_matrix_eps` (Müller et al.'s iteration, from the
     identity, `eps` the precision's `u`, at most $10^4$ iterations: at 100 it stops 0.80 short of
     the nearest rotation of $[[1,2,3],[4,5,6],[7,8,10]]$ and is within $5\times10^{-14}$ by 1000,
     and an oracle is scored at its best) and Sophus's `SO3::fitToSO3`
     (`JacobiSVD`, then $U\,\mathrm{diag}(1, 1, \det UV^\top)V^\top$), at both precisions.
   - **Bars**: domination of both on every stratum at both precisions; committed rows reproduced.
   - **Self-test**: `R` rotated by `1e-9` about $e_3$ fires on `polar:near-1e-8`.

## Rationale

- **Newton where it is well posed, not `svd3` everywhere**: about a third of `svd3`'s cost near
  SO(3) (two or three steps of one division and two roots, against eight rotations of three and
  two), half its error, and the only one of the two whose tangent is right there (*Context*).
- **A projector, not a reflection through `svd3`'s $v_3$, for $\det < 0$.** The first build took
  $X_0 = M - 2(Mv_3)v_3^\top$ with `svd3`'s $v_3$, then Newton: at $\Sigma = (2, 2, -\tfrac12)$ its
  tangent was $0.18$ off PO.3, because `svd3`'s rotation inside the tied pair has an arbitrary angle
  and an $O(1/u)$ lane, which the next rotation with $v_3$'s column carries into $v_3$'s lane at
  $O(1)$. The projector's tangent passes through no tie; `polar_tests` holds it there.
- **Frobenius scaling, not $\lvert\det X\rvert^{-1/3}$**: no cube root, and it is the scaling
  whose step count Higham bounds. The alternatives are rows of step 2's measurement.
- **No `nearest_rotation` without `svd3`**: the reflect and singular arms need it, and a Newton
  iteration on a singular input does not converge.
- **In `helicoid-linalg`**: the polar factor is a matrix factorization; `SO3` keeps its one
  constructor from a matrix and its contract.

`API.md` §6: (1) R4 holds, the arms are masks; (2) the id `nearest_rotation` and its strata;
(3) `# Domain`: every input, nothing asserted; (4) returns a `Mat3`, no `bool`, no dense
Jacobian; (5) it is the projection §3.4 names and nothing else spells; (6) no sensor, solver or
consumer state; (7) `no_std`, no allocation, D16; (8) it removes a way to be wrong: $UV^\top$ of
an SVD whose order or signs are not the signed convention is the nearest rotation only by luck
(SO.13(c)'s counterexample), and its tangent at a rotation is $O(1)$ wrong.

## Consequences

- `NUMERICS.md` §3.4 cites `nearest_rotation`; §11 gains `Rule::Polar`, §12 its domain row,
  §19 the formulas. `PHASE2.md` §6 and §0.0, `PHASE1.md` §4.3–§4.4 and `API.md` §2 list it.
- `SO3::from_matrix`'s rustdoc points at `nearest_rotation`.
- `maths/svd3.md` carries PO.1–PO.3.

## Implementation plan

1. With [`0101`](./0101-svd3-is-one-sided-jacobi-on-eig3-s-rotation.md) step 1: this record.
2. After `0101` step 2: `nearest_rotation`, its id, strata, metric, self-test, both oracles,
   `polar_tests`, committed rows, bench rows and *Measured* — verified by `polar_tests`,
   `the_unverified_routines_reproduce_their_committed_rows`, `just corpus-check`,
   `just conformance --self-test`, `just oracle-nalgebra`, `just oracle-sophus`, domination at
   both precisions.

## Open questions

None.

## Measured

**Large sample** (`polar_tests::measure_nearest_rotation`, $5\times10^4$ per family and
precision, the corpus metric's weighted error and frame against `svd3_tests`' double-double
reference; max, p99.9, mean in `u`):

| family | `f64` | `f32` |
|---|---|---|
| random | 6.15, 4.58, 1.61 | 5.90, 4.60, 1.61 |
| $cR$ | 5.08, 4.39, 1.60 | 5.64, 4.49, 1.53 |
| reflection $(1, 0.6, -0.2)$ | 6.44, 4.72, 1.62 | 6.10, 4.65, 1.63 |
| rank 2 | 5.93, 4.55, 1.56 | 5.80, 4.48, 1.56 |
| near-rotation $10^{-2}$ … $10^{-16}$ | 5.22–5.70, 4.41–4.61, 1.63–1.66 | 5.27–5.82, 4.43–4.51, 1.42–1.57 |
| $\sigma = (1, \tfrac12, 10^{-k})$, k = 4 … 16 | 5.32–5.70, 4.50–4.59, 1.56–1.58 | 5.33–5.81, 4.48–4.57, 1.55–1.56 |

**Corpus.** Per stratum, the maximum in `u` of `Rule::Polar`, for three programs:
- this;
- nalgebra 0.35's `Rotation3::from_matrix_eps` (Müller's iteration), its input brought to unit
  scale by a power of two, since its stopping test is absolute;
- Sophus's `SO3<T>::fitToSO3`, Eigen's `JacobiSVD`.

The maximum is 5.08 `u` at both precisions (`near-1e-4`, `scale-up@f32`), every one of them the
frame's term, and `polar:exact` reads 0. Domination holds on all 32 strata: the better oracle is at
least 4.24× this program's maximum at `f64` (`near-1e-4`) and 3.24× at `f32` (`near-1e-6@f32`).
Sophus is 9.2–33× this on every stratum. nalgebra reads $1.47\times10^{16}$ on `polar:exact` and
about $7\times10^3$ on `non-unique`, where the nearest rotation is not unique.

| stratum | this | nalgebra | Sophus |
|---|---|---|---|
| `polar:random` | 4 | 45.5 | 69.6 |
| `polar:scaled` | 3.44 | 23.6 | 31.8 |
| `polar:exact` | 0 | $1.47\times10^{16}$ | 22.6 |
| `polar:near-1e-2` | 3.46 | 20.5 | 116 |
| `polar:near-1e-4` | 5.08 | 21.5 | 160 |
| `polar:near-1e-8` | 4 | 19.7 | 67.6 |
| `polar:near-1e-12` | 4.15 | 22 | 39.3 |
| `polar:near-1e-16` | 3.62 | 23 | 44.5 |
| `polar:reflection` | 3.38 | 34.3 | 33.5 |
| `polar:cond-1e-4` | 3.65 | 31.6 | 80.2 |
| `polar:cond-1e-8` | 4.09 | 25 | 71 |
| `polar:cond-1e-12` | 4.16 | 31.5 | 80.2 |
| `polar:cond-1e-16` | 4.28 | 34.6 | 69.1 |
| `polar:rank2` | 3.36 | 24.1 | 61.3 |
| `polar:non-unique` | 4.92 | $7.13\times10^3$ | 48.3 |
| `polar:scale-up` | 4.14 | 45.8 | 48.2 |
| `polar:scale-down` | 3.64 | 56.8 | 38.1 |
| `polar:random@f32` | 4.41 | 35.3 | 66.7 |
| `polar:scaled@f32` | 3.8 | 15.2 | 38.1 |
| `polar:exact@f32` | 0 | $2.74\times10^{7}$ | 11.3 |
| `polar:near-1e-2@f32` | 4.37 | 18.3 | 50.4 |
| `polar:near-1e-4@f32` | 3.8 | 13.3 | 57.1 |
| `polar:near-1e-6@f32` | 4.36 | 14.1 | 62.9 |
| `polar:reflection@f32` | 3.26 | 17.7 | 89.5 |
| `polar:cond-1e-2@f32` | 4.45 | 20.8 | 48 |
| `polar:cond-1e-4@f32` | 4.04 | 18.5 | 65.8 |
| `polar:cond-1e-6@f32` | 2.99 | 30.1 | 57.5 |
| `polar:cond-1e-8@f32` | 4.32 | 20.1 | 55.7 |
| `polar:rank2@f32` | 3.46 | 21.1 | 64.5 |
| `polar:non-unique@f32` | 3.19 | $7.27\times10^3$ | 38.3 |
| `polar:scale-up@f32` | 5.08 | 26.6 | 98.9 |
| `polar:scale-down@f32` | 3.7 | 87.3 | 78.7 |

**Cost** (`benches/linalg.rs`, one pinned core, the median of three runs interleaved with the
commit before this record's review; `f64`, `f32`): a rotation $10^{-8}$ off SO(3) 99.5 and 52.6 ns;
`random`, Newton's arm, 286 and 227; the reflection $(1, 0.6, -0.2)$ 603 and 486; rank 2, `svd3`'s
arm, 525 and 468. The review put the ½ in Newton's scalars, reused the arm test's cofactors in the
first step, took `eig3`'s eigenvalues without its frame and put the rare arms out of line: Newton's
arm runs 0.71× (`f64`) and 0.68× (`f32`) as long, the reflection 0.80× and 0.75×; `svd3`'s arm and
`eig3` did not move.

## Further work

1. GTSAM's `Rot3::ClosestTo` as a third oracle, at the next rebuild of its image.
2. Byers–Xu's suboptimal scaling, which needs no norm after the first step.
