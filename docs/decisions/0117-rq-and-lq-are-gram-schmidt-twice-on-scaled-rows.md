# 0117: RQ and LQ are Gram–Schmidt twice on scaled rows

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** the PRs of its *Implementation plan*

## Context

A sensor-model consumer factors 3×3 matrices in two ways. One is $M = UR$, with $U$ upper triangular
and $R$ a rotation. The other is $M = LQ$, with $L$ lower triangular. It needs them for an IMU's
intrinsics, whose gauge fixes one axis of the sensor triad with a triangular matrix, and to convert
calibrations between tools that fix different axes: Basalt's accelerometer scale is lower
triangular, and the consumer's own model is upper triangular. It has four copies of the
factorization. Each is classical Gram–Schmidt once, normalising through a reciprocal, with no
evidence behind it. `helicoid-linalg` has `svd3` and `nearest_rotation` but no triangular
factorization.

What the consumer needs:

- $R$ (or $Q$) a proper rotation, a positive diagonal where $\det M > 0$, and a mask that reports
  it;
- a triangle with a positive diagonal returned with the identity rotation, bit for bit, because its
  round trip through Basalt's format is bit-exact;
- `Dual` derivatives;
- both precisions, `no_std`, no `alloc`.

**The candidates** (`maths/rq3.md` RQ.3–RQ.8):

- classical Gram–Schmidt (CGS) once, with the cross product for $r_1$, as the consumer has it;
- CGS twice (CGS2) with `eig3`'s Newton step on each unit row;
- Givens rotations;
- Householder reflections.

With three rows, modified Gram–Schmidt projects each row against one earlier row and is CGS. Each
candidate is measured.

**Protocol.** A scratch crate (not committed) depends on `helicoid-linalg` by path and builds the
candidates generically over `Real`. It was built with `--release`, without `target-cpu`, and with
no row scaling: the families are $O(1)$, and scaling is RQ.9's separate argument. It draws
$5\times10^4$ matrices per family and precision, each rounded entrywise. The scores are:

- the rowwise backward error $\max_i\lVert m_i - (\hat U\hat R)_i\rVert/(\lVert m_i\rVert u)$ and
  the orthogonality $\lVert\hat R^\top\hat R - I\rVert_F/u$, both accumulated in double-double with
  error-free products;
- a flag for $\det\hat R \le 0$, $\hat u_{22} < 0$ or $\hat u_{33} < 0$ (none was raised by CGS2 or
  Givens).

The families:

- `random`, Gaussian entries;
- `near-identity-1e-3`, $I + 10^{-3}G$, the calibration case;
- `near-triangular-1e-k`, with the lower entries scaled by $10^{-k}$;
- `parallel-1e-k`, with $m_2$ at angle $10^{-k}$ from $m_3$;
- `row-graded` and `row-graded-rev`, rows scaled by $(1, 10^{-5}, 10^{-10})$ or its reverse
  (`f32`: $(1, 10^{-2}, 10^{-4})$);
- `col-graded`, the columns scaled the same way;
- `near-singular-1e-k`, with $m_1$ within $10^{-k}$ of the span of $m_2, m_3$;
- `negative-det`;
- `rotation`, Haar distributed;
- `rank-bottom`, with $m_2 = m_3/2$.

**Accuracy.** Each cell is the maximum backward error over the maximum orthogonality, in `u`.

| family (`f64`) | CGS once | CGS2 + Kahan | Givens | Householder |
|---|---|---|---|---|
| `random` | 121 / 176 | **4.57 / 5.22** | 8.33 / 11.0 | 14.2 / 15.7 |
| `near-identity-1e-3` | 10.7 / 10.7 | **5.16 / 4.53** | 9.85 / 11.8 | 14.4 / 18.4 |
| `near-triangular-1e-8` | 18.2 / 30.1 | **4.62 / 4.15** | 8.00 / 10.4 | 12.0 / 15.6 |
| `parallel-1e-4` | 5.8e4 / 8.7e4 | **5.09 / 4.99** | 8.31 / 10.5 | 14.6 / 17.8 |
| `parallel-1e-12` | 6.6e12 / 9.5e12 | **4.82 / 5.73** | 8.10 / 9.70 | 13.5 / 16.5 |
| `parallel-1e-16` | 9.0e15 / 1.6e16 | **5.50 / 6.40** | 8.76 / 11.3 | 15.0 / 17.2 |
| `row-graded` | 248 / 423 | **4.75 / 5.13** | 8.18 / 9.87 | 14.6 / 17.1 |
| `col-graded` | 3.1e9 / 4.4e9 | **4.28 / 4.14** | 6.09 / 8.66 | 10.7 / 16.3 |
| `near-singular-1e-12` | 562 / 794 | **4.40 / 4.98** | 8.26 / 11.1 | 14.0 / 16.8 |
| `negative-det` | 910 / 2804 | **4.37 / 5.64** | 8.68 / 10.6 | 14.5 / 16.2 |
| `rotation` | 10.7 / 12.8 | **5.11 / 4.76** | 8.76 / 12.8 | 12.2 / 16.7 |
| `rank-bottom` | 9.0e15 / 1.6e16 | **4.47 / 5.44** | 8.88 / 8.43 | 14.6 / 15.1 |

At `f32` (k ∈ {2, 4, 6}, and `parallel-1e-7` to `-1e-9` at rounding level) the order is the same.
CGS2 reads at most 5.56 backward and 6.67 orthogonality. Givens reads 8.99 and 12.5, Householder 15.3
and 20.7. CGS once reaches $1.7\times10^7$.

**Kahan's test is needed, and its factor is measured.** Projecting twice is not enough on its own.
When $m_2$ lies in the span of $m_3$ to rounding, the first residual $w_1$ is rounding noise. That
noise lies almost along $r_3$, because the largest component's rounding dominates it, and the
second projection cancels again. `rank-bottom` then read 2352 / 3328 at `f64` and 9084 / 13000 at
`f32`. One draw left $\lVert w_2\rVert/\lVert w_1\rVert = 1.2\times10^{-4}$ and
$r_2\cdot r_3 = 1.9\times10^{-12}$.

Kahan's test (RQ.4) keeps $w_2$ only where $\lVert w_2\rVert^2 \ge c\lVert w_1\rVert^2$. Otherwise
the row is dependent at working precision. Measured maxima at the angles where the test decides,
orthogonality in `u`:

| $c$ | `parallel-1e-16` (`f64`) | `parallel-1e-17` | `parallel-1e-7` (`f32`) | `parallel-1e-8` |
|---|---|---|---|---|
| none | 262 | 4853 | 30.3 | 1165 |
| $\tfrac1{16}$ | 24.6 | 21.2 | 20.3 | 22.4 |
| $\tfrac14$ | 9.74 | 7.73 | 10.8 | 8.49 |
| $\tfrac12$ | **6.40** | **5.46** | **6.67** | **6.03** |

RQ.4's bound is $\sqrt{1/c}\,\gamma_4$. A larger $c$ tightens it, but it also declares more rows
dependent, and $c = \tfrac12$ (Parlett's $\kappa = \sqrt2$) still declares only those within a few
$u$ of the span, as the proof shows. $c = 1$ would declare a row dependent for any rounding loss.

**Exactness.** On $10^4$ draws, $\mathrm{diag}(d)$ with $d_i = 10^{U(-30, 30)}$, and the same
diagonal with Gaussian entries above it, counting $(U, R) = (M, I)$ bit for bit:

| | CGS once (reciprocal) | CGS once (divisions) | CGS2 | Givens | Householder |
|---|---|---|---|---|---|
| diagonal | 7470 | 10000 | 10000 | 10000 | 0 |
| upper triangle | 7470 | 10000 | 10000 | 10000 | 0 |

A reciprocal breaks exactness: $d\cdot\mathrm{fl}(1/d) \ne 1$ for a quarter of the draws. Only a
division keeps it (RQ.8).

**`Dual`.** The candidates were evaluated at `Dual<f64, 1>` along a Gaussian direction and checked
against RQ.6's tangent, $\dot R = \Omega R$, $\dot U = \dot MR^\top - U\Omega$. Each cell is the
maximum relative error over 2000 draws:

| point | CGS2 | Givens | Householder |
|---|---|---|---|
| random | 1.8e-15 | 7.4e-15 | 3.6e-15 |
| diagonal | 3.9e-16 | 4.1e-16 | 3.3e-15 |
| upper triangle | 4.8e-16 | 7.2e-16 | 3.0e-15 |
| $m_3 = (0, b, 0)$, $m_2 = (0, 0, c)$ | 5.6e-16 | **1.9** | 2.7e-15 |

Every fixed Givens order has a set where a rotation's pair is $(0, 0)$. RQ is smooth there, but the
rotation's angle is not, and the safe arm freezes the tangent. Gram–Schmidt is smooth wherever RQ.1
is unique.

**Cost.** `llvm-mca-18 -mcpu=znver3` on each candidate inlined into an `extern "C"` wrapper (`f64` /
`f32`):

| | instructions | divisions | roots | block reciprocal throughput (cycles) |
|---|---|---|---|---|
| CGS once | 275 / 278 | 6 | 3 | 57.0 / 45.3 |
| CGS2 | 296 / 298 | 6 | 3 | 57.0 / 48.5 |
| Givens | 229 / 219 | 6 | 3 | 57.0 / 35.0 |
| Householder | 248 / 254 | 2 | 2 | 40.5 / 41.3 |

The divisions bound all four at `f64`. Times on a quiet host are owed to *Measured*: the host stood
at load 26–36 on 8 cores throughout the study. In this table the completion arm (decision 2) is
evaluated branch-free. Decision 2 takes it behind `S::branch`.

## Decision

`NUMERICS.md` gains §22, which states decisions 1–6. `docs/maths/rq3.md` proves RQ.1–RQ.9.

1. **The contract.**
   - `rq3<S: Real>(m: &Mat3<S>) -> (UpperTri3<S>, Mat3<S>, S::Mask)` returns $(U, R, \mathrm{valid})$.
     For every finite input, $M = UR$ to a few $u$ of each row's norm, and
     $R \in \mathrm{SO}(3)$ to a few $u$ by construction.
   - $u_{22} \ge 0$ and $u_{33} \ge 0$ by construction, and $u_{11}$ takes the sign of $\det M$.
     Where $m_2$ and $m_3$ are independent the factorization is the unique one of RQ.1, whatever the
     sign of $\det M$.
   - `valid` is $u_{11} > 0 \wedge u_{22} > 0 \wedge u_{33} > 0$ on the computed (scaled) values: a
     statement about the computed values, as `chol`'s is (§15.2). In exact arithmetic it is
     $\det M > 0$.
   - Under a clear mask the factorization is kept, not replaced (the owner, 2026-11-15). A
     consumer with a negative determinant still holds an exact $M = UR$.
   - `lq3<S: Real>(m: &Mat3<S>) -> (LowerTri3<S>, Mat3<S>, S::Mask)` is
     $(\mathrm{flip}(U), \mathrm{flip}(R), \mathrm{valid})$ of `rq3(&flip(m))`. That is RQ.2's
     convention: $l_{11} \ge 0$ and $l_{22} \ge 0$, with $l_{33}$ taking the sign of $\det M$. It
     holds **bit for bit** by construction.
   - Accuracy is rowwise and backward (RQ.5). Forward errors are within a few $u$ times RQ.6's
     $\kappa_R$ and $\kappa_U$.
2. **The kernel** (on scaled rows, decision 3):
   1. **Row 3.** $a_3 = m_3\cdot m_3$, $u_{33} = \sqrt{a_3}$, and $r_3 = N(m_3/u_{33})$, three
      divisions. $N$ is `eig3`'s Newton step, $v(3/2 - \lVert v\rVert^2/2)$. Where $a_3$ is below
      the smallest normal, the row is zero after scaling: $u_{33} = 0$ and $r_3 = e_3$.
   2. **Row 2.** $w_1 = m_2 - (m_2\cdot r_3)r_3$, $w_2 = w_1 - (w_1\cdot r_3)r_3$,
      $a_2 = w_2\cdot w_2$. The row is dependent where $a_2 < \tfrac12\,w_1\cdot w_1$ (RQ.4, Kahan's
      test with Parlett's $\kappa = \sqrt2$) or $a_2$ is below the smallest normal.
      - Where it is independent: $u_{22} = \sqrt{a_2}$ and $r_2 = N(w_2/u_{22})$.
      - Where it is dependent: $u_{22} = 0$, and $r_2 = N(x/\lVert x\rVert)$ with
        $x = e_k - (r_3)_k\,r_3$, where $k$ is the axis of $r_3$'s smallest $\lvert$component$\rvert$,
        the **last** of equals. It is computed behind `S::branch`. The last of equals makes the
        zero matrix give $R = I$.
   3. **Row 1.** $r_1 = N(r_2\times r_3)$.
   4. **The rest of U.** $u_{11} = m_1\cdot r_1$, $u_{12} = m_1\cdot r_2$, $u_{13} = m_1\cdot r_3$,
      $u_{23} = m_2\cdot r_3$. Dot products sum left to right, as `Vector::dot` does.
   5. **No loops.** There is no data-dependent loop. The only decisions are the safe arms above and
      decision 3's ladder, through `S::select` and `S::branch`.
3. **Scaling** (RQ.9).
   - Each row is brought into `eig3`'s window by `eig3::into_window` on its three entries, with its
     own exact power of two. The windows are $[2^{-128}, 2^{128})$ (`f64`) and
     $[2^{-32}, 2^{32})$ (`f32`); the radicand is of degree two.
   - Row $i$ of $U$ is multiplied back by its own `back`, and $R$ is not scaled.
   - An entry beyond the largest finite number after scaling back is infinite. That is the only way
     a finite input yields a non-finite output.
4. **Non-finite input** gives $(I, I, \mathrm{false})$: the rows run on zeros (the safe argument),
   and a final blend replaces the outputs. It is never a plausible factorization, and every output
   is finite. The zero matrix gives $(0, I, \mathrm{false})$. Nothing panics or asserts (R6).
5. **Exactness** (RQ.7–RQ.9). An upper triangle with a positive diagonal, $+0$ below it and no
   $-0$ entry returns
   $(M, I)$ bit for bit at `f64`, `f32` and in `Dual` values, within the reach of RQ.9. That
   reach is: each row's entries within $2^{894}$ (`f64`), $2^{94}$ (`f32`) of the row's largest,
   and $m_{22}$ within $2^{383}$, $2^{31}$ of it. Every diagonal is inside it, subnormal entries
   included. Beyond it $u_{22}$ reads 0, which is still backward stable. `lq3` returns the same for a lower triangle.
   A $-0$ off the diagonal leaves the same values, with zeros of either sign.
6. **`Dual`.** The result differentiates the construction. Wherever RQ.1 is unique it is RQ.6's
   tangent, diagonal and triangular inputs included. The dependent arm and a zero row are constants,
   so the lanes there are finite and in general no directional derivative.
7. **Triangle types.** `UpperTri3<S>` and `LowerTri3<S>`, public in `helicoid-linalg` and
   re-exported at the root, hold their six entries row by row in a private `[S; 6]`. The layout is
   not a contract (D2). Their public items:
   - `identity()`;
   - `from_row_major([S; 6])` and `to_row_major() -> [S; 6]` (R3: the order is in the name);
   - `diagonal() -> Vec3<S>`;
   - `to_matrix() -> Mat3<S>`, with $+0$ in the zero triangle;
   - `transpose()`, from upper to lower and back;
   - `Blend<S>`, `Clone`, `Copy`, `Debug`.

   There is no arithmetic and no `mint`, which has no triangular types.
8. **No `SO3` wrapper in `helicoid`.** `SO3::from_matrix(&r)` is the one spelling (API.md §6 Q5).
   Its value at an exact identity is the identity quaternion bit for bit (trace 3, $w = \sqrt4/2$).
   A `helicoid` test pins that, because the consumer's bit-exact path composes the two.
9. **The corpus id `rq3`.** `lq3` is the §14 row `flip(rq3(flip m))`, bit for bit, with no id
   (`0104`'s precedent, `0108`).
   - **Fields.** `in`: `M` (3×3). `out`: `U` (3×3, exact zeros below), `R`, `valid`, and the
     reference-only `kappa_U` and `kappa_R` (RQ.6's norms, as `eig_sym`'s `kappa`).
   - **Reference.** RQ.1's construction at 120 digits, rechecked at 150.
     - `valid` and the sign of $u_{11}$ come from $\det M$ in rationals; the rank of $(m_2, m_3)$
       comes from `dense.rank` in rationals.
     - Where that rank is below 2, the reference completes $r_2$ (and $r_3$) by decision 2's axis
       rule, and both κ are stored as 0: not unique, the forward error unweighted.
     - The κ are the 2-norms of the linear map $\dot M \mapsto (\dot U, \dot R)$, a 9×9 matrix, from
       `mp.svd_r`.
     - Cross-check: $\lVert M - UR\rVert < 10^{-100}\lVert M\rVert$, $\lVert R^\top R - I\rVert <
       10^{-100}$, $\det R = 1$, and the convention.
   - **Metric `Rule::Rq`** (`NUMERICS.md` §11): the largest of
     - the rowwise backward error $\max_i\lVert m_i - (\hat U\hat R)_i\rVert/(\lVert m_i\rVert u)$,
       where a zero row reads 0 if reproduced exactly and $1/u$ otherwise;
     - $\lVert\hat R^\top\hat R - I\rVert_F/u$;
     - $\lVert\hat R - R\rVert_F/(\kappa_Ru)$ and $\lVert\hat U - U\rVert_F/(\kappa_U\lVert M\rVert_Fu)$,
       each weighing 0 where the stored κ is 0;
     - $1/u$ where $\det\hat R \le 0$, an entry below $\hat U$'s diagonal is nonzero,
       $\hat u_{22} < 0$ or $\hat u_{33} < 0$, or `valid` differs from $\hat U$'s diagonal test.

     Each term is a ratio of exact integers rounded once, as in `metric::forward_u`, except
     $\kappa$, which is read rounded from the record. A wrong sign of $u_{11}$ reads in the forward
     error: no mask tolerance is typed.
   - **Strata.** 16 records each and 64 in the deep tier, unless stated. Each draw comes from its
     own stream, and binary32 ranges are in parentheses.

     | Stratum | Draw |
     |---|---|
     | `rq:random` (32) | Gaussian entries |
     | `rq:diagonal` | $\mathrm{diag}(10^{U(-30,30)})$ ($10^{U(-6,6)}$): RQ.8: bit-exact, reading only the 30-digit reference's rounding |
     | `rq:triangular` | the same diagonal, Gaussian entries above: RQ.8, the same |
     | `rq:near-triangular-1e-k`, k ∈ {4, 8, 12} ({2, 4, 6}) | a positive triangle plus Gaussian entries of $10^{-k}$ below |
     | `rq:near-identity-1e-3` | $I + 10^{-3}G$: the calibration case |
     | `rq:parallel-1e-k`, k ∈ {2, 4, 8, 12, 16} ({2, 4, 6, 7}) | $m_2$ at angle $10^{-k}$ from $m_3$; the last value at rounding, where Kahan's test decides |
     | `rq:row-graded` | rows scaled by $(1, 10^{-5}, 10^{-10})$ ($(1, 10^{-2}, 10^{-4})$) |
     | `rq:row-spread` | rows scaled by $(10^{150}, 1, 10^{-150})$ ($(10^{18}, 1, 10^{-18})$): RQ.9 |
     | `rq:col-graded` | columns scaled by $(1, 10^{-5}, 10^{-10})$ ($(1, 10^{-2}, 10^{-4})$) |
     | `rq:near-singular-1e-k`, k ∈ {4, 8, 12} ({2, 4, 6}) | $m_1$ within $10^{-k}$ of the span of $m_2, m_3$ |
     | `rq:negative-det` | Gaussian with $\det < 0$ |
     | `rq:rank-bottom` | $m_2 = m_3/2$, exactly: forward error unscored |
     | `rq:zero` (4) | the zero matrix, signed zeros alternating |
     | `rq:rotation` | a Haar rotation, rounded |
     | `rq:scale-up`, `rq:scale-down` | `rq:random` times $10^{\pm70}$ ($10^{\pm12}$) |

     That is 388 binary64 and 372 binary32 records, 0.86 MB (862 589 bytes) of the 2.24 MiB `0110` leaves.
   - **Oracles.** nalgebra's `QR` and Eigen 3.4's `HouseholderQR<Matrix3<T>>` in the Sophus image, at
     both precisions.
     - Each factors $M^\top J = \tilde Q\tilde R$, then takes $U = J\tilde R^\top J$ and
       $R = J\tilde Q^\top$.
     - The runner puts the result in decision 1's convention by exact sign changes: rows 3 and 2 of
       $R$, with the matching columns of $U$, where $u_{33}$ or $u_{22}$ is negative; then row 1
       and column 1 where $\det R < 0$.
     - `valid` is the oracle's own diagonal test.
   - **Bars.** Domination of both oracles on every stratum at both precisions, and the committed
     rows of `committed_linalg.rs` reproduced bit for bit.
   - **Self-test.** $R$ turned by $10^{-9}$ about $z$, with $U$ kept, fires on `rq:random`
     (backward error). Rows 1 and 2 of $R$ and columns 1 and 2 of $U$ negated fires on
     `rq:random` (forward error and the diagonal test).

## Rationale

- **CGS2 over Givens.** It is at or below Givens on every family: 1.4–2.0× lower in backward error
  and 1.5–2.6× in orthogonality, except `near-triangular-1e-12`, where both read 1.81 / 0.
  It is smooth on the whole uniqueness domain, while every Givens order has a spurious singular set where `Dual` reads a wrong tangent (*Context*, relative error 1.9). Givens
  is 23% fewer instructions, and the divisions they share bound the time.
- **CGS2 over Householder.** It is 3× more accurate, exact on triangles where Householder never was,
  and needs no sign fix to reach $\det = +1$ with a positive diagonal.
- **Over CGS once.** The loss is not the general $\kappa^2u$ but $u/\sin\angle(m_2, m_3)$ (RQ.3).
  It is still unbounded: $10^{12}$ `u` at an angle of $10^{-12}$, and $3\times10^9$ on graded
  columns, where the bottom rows are nearly parallel. Its reciprocal breaks exactness on a quarter
  of diagonals. The second projection costs 21 instructions.
- **Kahan's test is a branch point derived from RQ.4's bound, and its factor is measured.** No
  constant is fitted to the data. $\tfrac12$ is the largest power of two below 1 and is exact at
  both precisions.
- **One kernel and a flip, not two kernels.** The two functions cannot disagree. A mirrored kernel
  would sum its dot products in the other order and differ in the last bit. One corpus id proves
  both (`0108`).
- **Rows scaled one by one, not by one global power of two** (RQ.9). Row scaling is free, since $R$
  is invariant. A global scale loses $\mathrm{diag}(10^{200}, 1, 10^{-200})$'s last row entirely.
- **Triangle types, fixed at 3.** The triangle is the consumer's parameter block and the type states
  which half is meaningful. A `Mat3` would leave three entries that every caller must zero or
  ignore. A general `UpperTri<S, N>` needs `[S; N(N+1)/2]`, which stable Rust cannot express
  (`generic_const_exprs`). That and arithmetic are *Further work*.
- **No `SO3` wrapper** (decision 8). It would be a second spelling of `from_matrix`, and
  `helicoid-linalg` cannot name `SO3`.

## Consequences

- **Specs.**
  - `NUMERICS.md` gains §22 (decisions 1–6), the §11 `rq3` metric, a §12 domain row and a §14 row
    for `lq3`.
  - `PHASE2.md` gains a §0.0 row, `rq3.rs` and `tri.rs` in §1, the triangle types in §4 and the
    signatures in §6.
  - `PHASE1.md` §4.3 and §4.4 gain the id and strata.
  - `API.md` §2 gains the items, and R6's list of functions that assert nothing gains `rq3` and
    `lq3`, reporting through a mask.
- **The instrument.**
  - `coverage::required()` gains `rq3`, and `committed_linalg::IDS` its rows.
  - The corpus grows by 0.86 MB to 78.58 MiB, leaving 1.42 MiB under `0110`'s cap.
- **Bits.** No existing output moves its bits: `eig3`'s helpers are reused, not changed.

## Implementation plan

1. This record and `maths/rq3.md`, as drafts.
2. `rq3`, `lq3`, `UpperTri3` and `LowerTri3` with their instrument. That means the kernel, the
   types, the unit tests and proptests, the generator, the metric, both oracles, the self-test, the
   committed rows, the bench rows and this record's *Measured* section. It is verified by:
   - `rq3_tests`, `tri_tests`;
   - `the_unverified_routines_reproduce_their_committed_rows`;
   - `just corpus-check`, `just conformance --self-test`;
   - `just oracle-nalgebra`, `just oracle-sophus`;
   - domination at both precisions;
   - `just determinism`.

## Open questions

None. The owner answered the failure convention, the wrapper and the storage, and accepted the
record, on 2026-11-15 (decisions 1, 7, 8).

## Measured

**Corpus.** The table gives the maximum of `Rule::Rq` in `u`, per stratum, for three programs:

- this program;
- nalgebra 0.35's `QR`, which is Householder;
- Eigen 3.4's `HouseholderQR<Matrix3<T>>`, in the Sophus image.

Both oracles are reversed into RQ and put in decision 1's convention by their runners.

- This program's maximum is 4.33 `u` at `f64` (`parallel-1e-4`) and 4.00 at `f32` (`parallel-1e-7`).
- Domination holds on all 47 strata. Wherever this program reads above rounding, the better oracle
  is at least 1.52× its maximum at `f64` (`rotation`) and 1.50× at `f32` (`parallel-1e-7`).
- `rq:diagonal` ties with Eigen, which is exact there too. `rq:zero` ties at 0.
- `rq:diagonal` and `rq:triangular` read only the 30-digit reference's rounding: below `3e-14`.
  The bits themselves are the unit tests' (decision 5). Eigen is not exact on a triangle (1.17), and
  nalgebra on neither (8).

| stratum | this | nalgebra | Eigen |
|---|---|---|---|
| `rq:random` | 3.65 | 13.4 | 7.53 |
| `rq:diagonal` | 2.6e-14 | 8 | 2.6e-14 |
| `rq:triangular` | 1.4e-16 | 8 | 1.17 |
| `rq:near-triangular-1e-4` | 3.11 | 13.5 | 9.75 |
| `rq:near-triangular-1e-8` | 2.09 | 11 | 8.25 |
| `rq:near-triangular-1e-12` | 1.06 | 7.07 | 4 |
| `rq:near-identity-1e-3` | 2.35 | 12.2 | 6.5 |
| `rq:parallel-1e-2` | 3.39 | 12.2 | 8.12 |
| `rq:parallel-1e-4` | 4.33 | 12.6 | 6.7 |
| `rq:parallel-1e-8` | 3.41 | 10.1 | 6.33 |
| `rq:parallel-1e-12` | 3.22 | 12.1 | 7.1 |
| `rq:parallel-1e-16` | 3.56 | 12.2 | 9.15 |
| `rq:row-graded` | 3.67 | 13.3 | 6 |
| `rq:row-spread` | 3.41 | 11.8 | 7.88 |
| `rq:col-graded` | 2.95 | 10.7 | 6.12 |
| `rq:near-singular-1e-4` | 3.16 | 11.2 | 8.32 |
| `rq:near-singular-1e-8` | 3.23 | 10.2 | 8.18 |
| `rq:near-singular-1e-12` | 3.22 | 11.9 | 9.61 |
| `rq:negative-det` | 3.3 | 11.1 | 11 |
| `rq:rank-bottom` | 3.38 | 12.1 | 7.61 |
| `rq:zero` | 0 | 0 | 0 |
| `rq:rotation` | 3.04 | 11.2 | 4.62 |
| `rq:scale-up` | 2.81 | 13 | 8.37 |
| `rq:scale-down` | 2.86 | 12 | 6.78 |
| `rq:random@f32` | 3.82 | 12.4 | 6.47 |
| `rq:diagonal@f32` | 8.4e-26 | 8 | 8.4e-26 |
| `rq:triangular@f32` | 1.1e-31 | 10.9 | 1.75 |
| `rq:near-triangular-1e-2@f32` | 2.89 | 13.3 | 8.61 |
| `rq:near-triangular-1e-4@f32` | 1.72 | 11.4 | 5.14 |
| `rq:near-triangular-1e-6@f32` | 0.787 | 9.06 | 4.26 |
| `rq:near-identity-1e-3@f32` | 2.76 | 12.2 | 10.4 |
| `rq:parallel-1e-2@f32` | 3.54 | 14.2 | 6.2 |
| `rq:parallel-1e-4@f32` | 3.04 | 13.6 | 6.7 |
| `rq:parallel-1e-6@f32` | 3.2 | 11.4 | 7.47 |
| `rq:parallel-1e-7@f32` | 4 | 12.5 | 5.98 |
| `rq:row-graded@f32` | 2.53 | 13.9 | 7.22 |
| `rq:row-spread@f32` | 2.75 | 12.8 | 7.32 |
| `rq:col-graded@f32` | 2.28 | 11 | 6.85 |
| `rq:near-singular-1e-2@f32` | 3.03 | 13.8 | 6.54 |
| `rq:near-singular-1e-4@f32` | 3.48 | 16.5 | 8.68 |
| `rq:near-singular-1e-6@f32` | 3.02 | 15.1 | 8.96 |
| `rq:negative-det@f32` | 3.39 | 14.7 | 11.3 |
| `rq:rank-bottom@f32` | 3.38 | 11.5 | 5.33 |
| `rq:zero@f32` | 0 | 0 | 0 |
| `rq:rotation@f32` | 3.01 | 9.33 | 7 |
| `rq:scale-up@f32` | 3.21 | 13.7 | 8.03 |
| `rq:scale-down@f32` | 3.37 | 12.6 | 5.58 |


**Large sample.** `measure_rq3` ran 10⁶ draws of each of *Context*'s families at each precision,
with two seeds, scored in `f64` from the output (backward error and frame, not the forward terms).
The worst readings were:

- `f64`: 8.0 backward, 7.9 frame;
- `f32`: 6.6 backward, 7.5 frame.

At `f64` they include the metric's own rounding. `rq3_tests` records `K = FRAME = 16`.

**Exactness.** `a_triangle_returns_itself_and_the_identity_bit_for_bit` passes on 2×10⁴ draws per
precision, covering both functions and the `Dual` value path:

- a diagonal over every binade, subnormal entries included;
- the same diagonal with entries above it, within `2^30` of their row's diagonal;
- a lower triangle for `lq3`.

`lq3_matches_reference` holds bit for bit on 4096 proptest cases per precision, with non-finite
entries included.

**`Dual`.** `derivative_is_first_order_perturbation_theory` agrees with RQ.6 to `1e-12` relative on
every entry's lane at 204 matrices: random ones, a diagonal, a triangle, a reflection and Givens'
singular point.

**Self-test.** Both defects fire on `rq:random` against a clean maximum of 3.65. `TurnsTheFactor`
reads `9.0e6` and `FlipsTheGauge` reads `9.5e15`.

## Further work

- Arithmetic on the triangle types (`U v`, `U⁻¹ v`), once a consumer evaluates a model through them.
- A general `UpperTri<S, N>`, if `chol`'s factor ever wants a type and stable Rust can say its size.
