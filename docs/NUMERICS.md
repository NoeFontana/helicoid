# helicoid — Numerics: formulas, series, domains

> **Companions:** [`PROJECT.md`](./PROJECT.md) D3–D9, D16;
> [`0002`](./decisions/0002-one-convention-for-a-stack-that-already-disagrees.md) (conventions),
> [`0004`](./decisions/0004-switch-points-are-generated-not-typed.md) (switch points),
> [`0005`](./decisions/0005-the-jacobian-is-a-dual-matrix.md) (Jacobian structure),
> [`0006`](./decisions/0006-the-instrument-comes-first.md) (the instrument).

Every formula a phase implements, **already stated in this repository's conventions**
(rotation-first, right perturbation, Hamilton `w`-first). Literature formulas (Barfoot, Solà,
Sophus) are translation-first; they have been permuted here. **Never permute them again in your
head: implement what is written.** The whole document is **NORMATIVE** except §13.

A formula here is a claim; the corpus is its test. The generator never uses these closed forms: it
computes from definitions (`PHASE1.md` §4.3), so a wrong line in this file fails the envelope.

## 0. Status

| Section | Status |
|---|---|
| §1–§8, §10–§12, §14, §15 | **Ready** |
| §9 Sim(3) | **Owed**: group law and `Ad` stated; `Exp`, `Log`, `Jr` block must be derived and generator-verified before `PHASE5.md` §3 starts |

## 1. Conventions

- **Quaternions.** $q = (w, x, y, z) = w + x\,i + y\,j + z\,k$, Hamilton ($i^2 = j^2 = k^2 = ijk = -1$),
  stored `w` first. Active: $v' = q\,v\,q^*$. With $u = (x, y, z)$:
  $R(q) = (w^2 - \|u\|^2)\,I + 2\,u u^\top + 2w\,[u]_\times$.
- **SE_N(3).** $X = (R, x_1, \dots, x_N)$, i.e. the $(3+N)\times(3+N)$ matrix
  $\begin{bmatrix} R & x_1 \cdots x_N \\ 0 & I_N \end{bmatrix}$. Composition
  $XY = (R_X R_Y,\ R_X y_i + x_i)$; inverse $(R^\top, -R^\top x_i)$. SE(3) is $N = 1$ with $x_1 = t$;
  `a * b` is $T_{a x}\,T_{x b}$ and `X * p` is $R p + t$.
- **Tangents are rotation-first.** $\tau = [\varphi;\ \rho_1;\ \dots;\ \rho_N] \in \mathbb{R}^{3+3N}$;
  a twist is $[\omega; v]$. Hat:
  $\tau^\wedge = \begin{bmatrix} \varphi^\wedge & \rho_1 \cdots \rho_N \\ 0 & 0 \end{bmatrix}$,
  $\varphi^\wedge = [\varphi]_\times$. SE(2): $[\theta; \rho]$. Sim(3): $[\varphi; \rho; \sigma]$.
- **Exp/Log.** $\mathrm{Exp}(\tau) = \exp(\tau^\wedge)$. $\mathrm{Log}$ is its inverse on the canonical
  branch $\theta = \|\varphi\| \in [0, \pi]$.
- **Perturbations.** Right (default): $X \oplus_R \tau = X\,\mathrm{Exp}(\tau)$,
  $Y \ominus_R X = \mathrm{Log}(X^{-1} Y)$. Left: $X \oplus_L \tau = \mathrm{Exp}(\tau)\,X$,
  $Y \ominus_L X = \mathrm{Log}(Y X^{-1})$.
- **Adjoints.** $X\,\mathrm{Exp}(\tau)\,X^{-1} = \mathrm{Exp}(\mathrm{Ad}_X\,\tau)$;
  $\mathrm{ad}_\tau\,\sigma = [\tau^\wedge, \sigma^\wedge]^\vee$.
- **Jacobians of Exp.** $\mathrm{Exp}(\tau + \delta) \approx \mathrm{Exp}(\tau)\,\mathrm{Exp}(J_r(\tau)\,\delta)
  \approx \mathrm{Exp}(J_l(\tau)\,\delta)\,\mathrm{Exp}(\tau)$, and
  $J_l(\tau) = J_r(-\tau) = \mathrm{Ad}_{\mathrm{Exp}(\tau)}\,J_r(\tau)$.

## 2. Notation and generic identities

### 2.1 Notation

$\theta^2 = \varphi \cdot \varphi$ is always computed as a dot product, never as the square of a
`sqrt`. $\theta = \sqrt{\theta^2}$ is formed only inside an exact arm, at the safe argument
(`PROJECT.md` §7). $W = \varphi^\wedge$. $u$ is the unit roundoff of the precision: $2^{-53}$ (`f64`,
`Dual<f64, N>`) or $2^{-24}$ (`f32`).

### 2.2 Dual-matrix notation for SE_N(3) Jacobians

A structured Jacobian `SEn3Jac { diag: A, col: [B_1, …, B_N] }` denotes the rotation-first dense
matrix

$$
\begin{bmatrix} A & 0 & \cdots & 0 \\ B_1 & A & & \\ \vdots & & \ddots & \\ B_N & & & A \end{bmatrix},
$$

written $A + \epsilon B$. Products and inverses (0005):
$(A + \epsilon B)(C + \epsilon D) = AC + \epsilon\,(B_i C + A D_i)$ and
$(A + \epsilon B)^{-1} = A^{-1} - \epsilon\,A^{-1} B_i A^{-1}$.

### 2.3 Sides and their Jacobians

Each side's Jacobians are expressed **in that side's perturbation convention**. Right: perturb as
$X \oplus_R \delta$. Left: perturb as $X \oplus_L \delta$.

| Operation | Right: $\partial/\partial$ first | Right: $\partial/\partial$ second | Left: $\partial/\partial$ first | Left: $\partial/\partial$ second |
|---|---|---|---|---|
| $X \oplus \tau$ | $\mathrm{Ad}_{\mathrm{Exp}(\tau)}^{-1}$ | $J_r(\tau)$ | $\mathrm{Ad}_{\mathrm{Exp}(\tau)}$ | $J_l(\tau)$ |
| $\tau = Y \ominus X$ ($\partial/\partial Y$, $\partial/\partial X$) | $J_r^{-1}(\tau)$ | $-J_l^{-1}(\tau)$ | $J_l^{-1}(\tau)$ | $-J_r^{-1}(\tau)$ |
| $XY$ | $\mathrm{Ad}_Y^{-1}$ | $I$ | $I$ | $\mathrm{Ad}_X$ |
| $X^{-1}$ | $-\mathrm{Ad}_X$ | — | $-\mathrm{Ad}_X^{-1}$ | — |
| $\mathrm{Exp}(\tau)$ | $J_r(\tau)$ | — | $J_l(\tau)$ | — |
| $\mathrm{Log}(X)$ | $J_r^{-1}(\mathrm{Log}\,X)$ | — | $J_l^{-1}(\mathrm{Log}\,X)$ | — |

### 2.4 Action Jacobians (SE(3), rotation-first)

$\partial(X p)/\partial X$: right $[\,-R[p]_\times,\ R\,]$; left $[\,-[Rp + t]_\times,\ I\,]$.
$\partial(Xp)/\partial p = R$. SO(3): right $-R[p]_\times$, left $-[Rp]_\times$.

## 3. SO(3)

### 3.1 Exp

$q = \left(\cos\tfrac{\theta}{2},\ k(\theta)\,\varphi\right)$, with $k$ from §4. One `sin_cos` of
$\theta/2$ in the exact arm; the series arm needs no transcendental for $k$ and computes
$\cos\tfrac{\theta}{2}$ from its own series (generated alongside $k$).

### 3.2 Log

1. If $w < 0$, $q \leftarrow -q$ (implemented as `copysign`, not a branch). At $w = +0$ nothing flips:
   **at exactly $\theta = \pi$, `Log` is a function of the quaternion, not of the rotation** — $q$
   and $-q$ return $\pm\pi\hat n$.
2. $n^2 = x^2 + y^2 + z^2$; $\theta = 2\,\mathrm{atan2}(n, w)$.
3. $\varphi = r(n^2, w)\,u$ with $r = 2\,\mathrm{atan2}(n, w)/n$ (§4).

`Log` is exactly scale-invariant in $q$ ($\mathrm{atan2}$ and $u/n$ are), so a non-unit input within
rounding returns the `Log` of its normalization. **Never $\mathrm{acos}((\mathrm{tr}R - 1)/2)$**
(D5): it resolves $\theta$ only to $\sim\sqrt{2u}$ near $0$ and $\pi$.

### 3.3 Action

$v' = v + 2w\,(u \times v) + 2\,u \times (u \times v)$. `act_many` forms $R(q)$ once and applies it
per point; its reference twin is the per-point `act`.

### 3.4 From a matrix

Shepperd's method: select the largest of $\{\mathrm{tr}R, R_{00}, R_{11}, R_{22}\}$ through nested
`branch`es and extract from that pivot (e.g. trace pivot: $w = \tfrac12\sqrt{1 + \mathrm{tr}R}$,
$x = (R_{21} - R_{12})/(4w)$, …), then normalize. **Closed form, never iterative**:
avoids iterative convergence failures (such as Müller iteration hangs on degenerate input).
Projecting an arbitrary $3\times3$ onto SO(3) is `svd3` first
($U\,\mathrm{diag}(1, 1, \det UV^\top)\,V^\top$), then this.

### 3.5 Jacobians and adjoints

$$
J_r = I - a\,W + b\,W^2,\quad J_l = I + a\,W + b\,W^2,\quad
J_r^{-1} = I + \tfrac12 W + c\,W^2,\quad J_l^{-1} = I - \tfrac12 W + c\,W^2,
$$

$\mathrm{Ad}_R = R$, $\mathrm{ad}_\varphi = W$.

### 3.6 Renormalization

Composition never normalizes. The explicit `renormalize` applies the first-order Newton step
$q \leftarrow q\,(3 - \|q\|^2)/2$ (exact to $O((\|q\|^2 - 1)^2)$, no `sqrt`); it repairs drift, and
is a normalization only for $|\|q\|^2 - 1| \le 2^{-26.29}$ (`f64`) or $2^{-11.79}$ (`f32`)
(`maths/so3.md` SO.15), returning $0$ at $\|q\|^2 = 3$ and reversing $q$ beyond.

Construction from external data **divides by the norm**: $q \leftarrow q/\|q\|$, one `sqrt` and one
division per component, whose domain is $\|q\|^2$ normal
([`0027`](./decisions/0027-a-normalizing-constructor-normalizes.md) — every `*_normalized`
constructor normalizes, as `API.md` R6 reads it, and the step belongs to `renormalize` alone).
`from_wxyz_unchecked` debug-asserts $|\|q\|^2 - 1| \le 2^{-40}$ (`f64`; `tf_tree_math`'s `1e-12`)
or $2^{-16}$ (`f32`); every $q$ it admits is inside the step's accuracy domain by 13.7 bits in
`f64` and 4.2 in `f32`.

## 4. The coefficient catalogue

**Only `helicoid::coeffs` evaluates these.** Each is a function of a branch variable ($\theta^2$,
or $n^2$ for $r$) with an exact arm and a series arm; the switch point and the number of series
terms are **generated per precision** by `cargo xtask thresholds`
([`0004`](./decisions/0004-switch-points-are-generated-not-typed.md)). The series coefficients are
generated from mpmath's Taylor expansion as exact rationals; **the leading four below are asserted
by the generator**, as a cross-check, not typed into code.

| Name | Definition | Exact arm computes | Series (leading four) | Naive cancellation | Prior |
|---|---|---|---|---|---|
| $k$ | $\sin(\theta/2)/\theta$ | definition (0/0 only) | $\tfrac12 - \tfrac{\theta^2}{48} + \tfrac{\theta^4}{3840} - \tfrac{\theta^6}{645120}$ | none | — |
| $a$ | $(1 - \cos\theta)/\theta^2$ | $2k^2$ | $\tfrac12 - \tfrac{\theta^2}{24} + \tfrac{\theta^4}{720} - \tfrac{\theta^6}{40320}$ | $\sim 2u/\theta^2$ naive; none after rewrite | `tf_tree` D12: $\theta < 0.1$, 4 terms |
| $b$ | $(\theta - \sin\theta)/\theta^3$ | definition — **no rewrite exists** | $\tfrac16 - \tfrac{\theta^2}{120} + \tfrac{\theta^4}{5040} - \tfrac{\theta^6}{362880}$ | $\sim 6u/\theta^2$ | `tf_tree` D12 |
| $c$ | $1/\theta^2 - (1 + \cos\theta)/(2\theta\sin\theta)$ | $1/\theta^2 - \cot(\theta/2)/(2\theta)$ (exact at $\pi$) | $\tfrac1{12} + \tfrac{\theta^2}{720} + \tfrac{\theta^4}{30240} + \tfrac{\theta^6}{1209600}$ | $\sim 12u/\theta^2$ | `tf_tree` D12 |
| $d$ | $(\theta^2 + 2\cos\theta - 2)/(2\theta^4)$ | $(\theta^2 - 4\sin^2(\theta/2))/(2\theta^4)$ | $\tfrac1{24} - \tfrac{\theta^2}{720} + \tfrac{\theta^4}{40320} - \tfrac{\theta^6}{3628800}$ | $\sim 24u/\theta^4$ naive, $\sim 24u/\theta^2$ rewritten | — |
| $e$ | $(2\theta - 3\sin\theta + \theta\cos\theta)/(2\theta^5)$ | definition — **no rewrite exists** | $\tfrac1{120} - \tfrac{\theta^2}{2520} + \tfrac{\theta^4}{120960} - \tfrac{\theta^6}{9979200}$ | $\sim 360u/\theta^4$ | — |
| $r$ | $2\,\mathrm{atan2}(n, w)/n$ | definition | $\tfrac{2}{w}\left(1 - \tfrac{s}{3} + \tfrac{s^2}{5} - \tfrac{s^3}{7}\right)$, $s = n^2/w^2$ | none (0/0) | — |

What the "naive cancellation" column means: evaluating the definition in the exact arm keeps
$-\log_{10}$ of that relative error in digits. $e$ is the reason one global threshold is wrong by
construction: at $\theta = 10^{-2}$ it keeps about five digits, at $10^{-3}$ about one.

**Call sites evaluate coefficients in groups, inside one `branch`:** `exp_coeffs` → $(k, \cos\tfrac\theta2)$;
`jr_coeffs` → $(a, b)$; `jr_inv_coeff` → $c$; `q_coeffs` → $(b, d, e)$; `gamma2_coeffs` → $(b, d)$;
`log_ratio` → $r$. A call site never evaluates one coefficient from the catalogue on its own.

**Continuity.** At every generated switch point, $|\text{series} - \text{exact}|$ is at most the
recorded max error of that coefficient; `branch_continuity_*` tests assert it.

## 5. SE_N(3)

### 5.1 Exp and Log

$\mathrm{Exp}(\tau) = (\mathrm{Exp}(\varphi),\ J_l(\varphi)\rho_1,\ \dots,\ J_l(\varphi)\rho_N)$.
$\mathrm{Log}(X) = (\varphi = \mathrm{Log}(R),\ \rho_i = J_l^{-1}(\varphi)\,x_i)$.
$J_l(\varphi)$ and $J_l^{-1}(\varphi)$ are §3.5's (the `V` and `V⁻¹` of `tf_tree_math`).

### 5.2 Adjoints

$\mathrm{Ad}_X = R + \epsilon\,[x_i]_\times R$, i.e. `diag = R`, `col[i] = [x_i]× R`.
$\mathrm{ad}_\tau = W + \epsilon\,[\rho_i]_\times$.

### 5.3 Jacobians

$$
J_l(\tau) = J_l(\varphi) + \epsilon\,Q(\rho_i, \varphi),\qquad
J_r(\tau) = J_r(\varphi) + \epsilon\,Q(-\rho_i, -\varphi),
$$

with Barfoot's block (Barfoot 2017; Barfoot & Furgale 2014), written with §4's coefficients:

$$
\begin{aligned}
Q(\rho, \varphi) ={}& \tfrac12\,\rho^\wedge
+ b\,\big(\varphi^\wedge\rho^\wedge + \rho^\wedge\varphi^\wedge + \varphi^\wedge\rho^\wedge\varphi^\wedge\big) \\
&+ d\,\big(\varphi^\wedge\varphi^\wedge\rho^\wedge + \rho^\wedge\varphi^\wedge\varphi^\wedge - 3\,\varphi^\wedge\rho^\wedge\varphi^\wedge\big)
+ e\,\big(\varphi^\wedge\rho^\wedge\varphi^\wedge\varphi^\wedge + \varphi^\wedge\varphi^\wedge\rho^\wedge\varphi^\wedge\big).
\end{aligned}
$$

In Barfoot's translation-first order this block sits upper-right; rotation-first puts it lower-left,
which is what `col[i]` holds. The block for each $\rho_i$ is independent of the others.

### 5.4 Inverses

$J_r^{-1}(\tau)$ and $J_l^{-1}(\tau)$ are **not separate formulas**: they are the dual-matrix inverse
(§2.2) of §5.3 with $A^{-1}$ from §3.5's closed form:
$J_r^{-1}(\tau) = J_r^{-1}(\varphi) - \epsilon\,J_r^{-1}(\varphi)\,Q(-\rho_i, -\varphi)\,J_r^{-1}(\varphi)$.
The dense Gauss–Jordan inverse of $J_r(\tau)$ is the reference twin (§14).

### 5.5 Action (N = 1)

$X p = Rp + t$, with §2.4's Jacobians.

## 6. SO(2) and SE(2)

SO(2) is a unit complex $(c, s)$: $\mathrm{Exp}\,\theta = (\cos\theta, \sin\theta)$,
$\mathrm{Log} = \mathrm{atan2}(s, c)$, $J_r = J_l = 1$, $\mathrm{Ad} = 1$.

SE(2), tangent $[\theta; \rho]$: $\mathrm{Exp} = (R(\theta),\ V(\theta)\rho)$ with
$V(\theta) = \begin{bmatrix} \alpha & -\beta \\ \beta & \alpha \end{bmatrix}$,
$\alpha = \sin\theta/\theta$, $\beta = (1 - \cos\theta)/\theta = 2\sin^2(\theta/2)/\theta$ (both 0/0
only: series from the generator, switch generated). $\mathrm{Log}$: $\theta = \mathrm{atan2}(s, c)$,
$\rho = V(\theta)^{-1} t$. $J_r$, $J_l$ and their inverses: Solà et al. 2018, Appendix (SE(2)),
**permuted to rotation-first** in the Phase 3 PR that implements them, and verified by `Dual` and
the corpus; that PR adds the permuted matrices to this section.

## 7. Integrated exponentials

$\Gamma_m(\varphi) = \sum_{n \ge 0} W^n/(n+m)!$. $\Gamma_0 = \mathrm{Exp}$; $\Gamma_1 = J_l$;
$\Gamma_2 = \tfrac12 I + b\,W + d\,W^2$. With piecewise-constant body-frame $\omega$ and $a$ over
$\Delta t$: $\Delta R = \mathrm{Exp}(\omega\Delta t)$,
$\Delta v = \Gamma_1(\omega\Delta t)\,a\,\Delta t$, $\Delta p = \Gamma_2(\omega\Delta t)\,a\,\Delta t^2$
(Barrau & Bonnabel 2020; Brossard et al. 2022). Preintegration itself is not `helicoid`'s
([`0009`](./decisions/0009-what-helicoid-does-not-own.md)).

**Directional Jacobians** $\partial(\Gamma_m(\varphi)\,v)/\partial\varphi$ are computed by
evaluating $\Gamma_m$ on `Dual<S, 3>` through `coeffs` — exact to rounding, Taylor branches
included. A closed form replaces this only if a bench shows the `Dual` path is a consumer
bottleneck (`PROJECT.md` §5.1).

## 8. S²

Storage: unit $n \in \mathbb{R}^3$. The chart at $n$ is frozen at construction:

- **Basis.** $\nu = n + \mathrm{sgn}(n_z)\,e_z$ with $\mathrm{sgn}(0) = +1$;
  $H = I - 2\,\nu\nu^\top/(\nu^\top\nu)$; $b_1 = H e_x$, $b_2 = H e_y$. $H e_z = -\mathrm{sgn}(n_z)\,n$,
  so $b_1, b_2 \perp n$; $\nu^\top\nu = 2(1 + |n_z|) \ge 2$ never cancels.
- **Retract.** $n \oplus \delta = \mathrm{Exp}(B\delta)\,n$, $B = [b_1\ b_2]$.
- **Local.** $\delta = B^\top\,(\alpha\,\hat m)$ with $\hat m = (n \times m)/\|n \times m\|$,
  $\alpha = \mathrm{atan2}(\|n \times m\|, n \cdot m)$; the ratio $\alpha/\|n \times m\|$ is §4's $r$
  kernel (with $w = n\cdot m$, $n^2 = \|n\times m\|^2$).
- **The basis is discontinuous at $n_z = 0$**, and by the hairy ball theorem every global basis is
  discontinuous somewhere. The chart is therefore computed once per linearization and reused; it is
  never recomputed inside $\oplus$ ([`0012`](./decisions/0012-a-retraction-is-a-chart.md)).

## 9. Sim(3) — owed

Storage $(q, t, \sigma)$ with scale $s = e^\sigma$ (so the inverse negates $\sigma$ exactly).
$X = \begin{bmatrix} sR & t \\ 0 & 1\end{bmatrix}$; $XY = (R_X R_Y,\ s_X R_X t_Y + t_X,\ \sigma_X + \sigma_Y)$.
Tangent $[\varphi; \rho; \sigma]$, hat $\begin{bmatrix} \varphi^\wedge + \sigma I & \rho \\ 0 & 0\end{bmatrix}$.
$\mathrm{Ad}_X[\varphi; \rho; \sigma] = [R\varphi;\ sR\rho + [t]_\times R\varphi - \sigma t;\ \sigma]$.

**Owed before `PHASE5.md` §3 starts**, as an edit to this section plus a record: $\mathrm{Exp}$'s
$W(\varphi, \sigma)$ block, $\mathrm{Log}$, $J_r$, their coefficients as functions of $(\sigma, \theta)$,
and the two-dimensional series at the joint limit $\sigma^2 + \theta^2 \to 0$ (Eade 2014; Strasdat
2012). Jacobians are dense $7 \times 7$ ([`0005`](./decisions/0005-the-jacobian-is-a-dual-matrix.md)).

## 10. Geodesics

$X(t) = X_0\,\mathrm{Exp}(t\,d)$ with $d = X_1 \ominus_R X_0$ and $\Delta = X_0^{-1}X_1$:

- **Definition = reference twin.** For SE(3) this is ScLERP; the unit-dual-quaternion power is the
  fast twin (`PHASE4.md` §1). For `Product<SO3, R3>` it is slerp + lerp, i.e. `tf2`'s semantics.
- **SO(3) ships GE.14's blend, and the denominator's spelling is NORMATIVE** ([`0050`](./decisions/0050-the-geodesic-s-denominator-is-the-whole-domination-gap.md)):
  $$\gamma(q_0, q_1, t) = \frac{\sin((1-t)\alpha)\,q_0 + \sin(t\alpha)\,q_1}{\sin\alpha},
  \qquad (w, \mathbf v) = q_0^{*}q_1 \text{ with } w \ge 0,\quad
  \alpha = \mathrm{atan2}(\lVert\mathbf v\rVert, w),$$
  with the denominator **$\sin\alpha$ recomputed from $\alpha$, never $\lVert\mathbf v\rVert$**,
  although for a unit quaternion the two are the same number. GE.14 proves this equal to
  $q_0\mathrm{Exp}(t\,\mathrm{Log}(q_0^{*}q_1))$, so it is not a different function and the
  definition above stands as the reference twin; what is normative is the spelling, because only the
  recomputed form makes both $t = 0$ and $t = 1$ exact — one number divided by itself — and that is
  the whole difference between dominating oracle #1 on all three `so3_geodesic` strata (1.644 /
  1.738 / 1.642 `u`) and losing two of them (1.572 / 2.721 / 2.429). $\lVert\mathbf v\rVert = 0$
  is the two rotations being equal and returns $q_0$; the domain below is unchanged, as the blend
  has none of its own.
- **Below `r`'s second switch the provided body is taken instead** ([`0051`](./decisions/0051-two-arms-on-the-switch-the-sweep-already-chose.md)), because it is
  both cheaper there — its coefficients are on their series arms, polynomials with no transcendental
  — and *more* accurate, 1.572 against 1.644 `u` at `geo:consecutive`. The boundary is
  `s = \tan^2\alpha <` `short_below`, `r`'s own generated number and not a new one, and **exactly**
  $t = 1$ takes the blend on both sides of it, since only the blend is exact there. Not $t \ge 1$:
  extrapolation below the switch belongs to the provided body, whose error does not grow with $t$
  where the blend's weights do — it forms $\approx -(t-1)q_0 + t q_1$, and the cancellation reads
  about $1.2\times10^6\,u$ at $t = 10^6$, $\theta(d) = 10^{-9}$. The two spellings are the same
  function (GE.14), so this is one routine with two arms and not two routines. At
  $\lVert v\rVert = 0$ the blend's own limit $(1-t)q_0 + t q_1$ is taken, which is exact at both
  endpoints even where $n^2$ underflows on a pair that still differs.
- **SE(3) ships the dual-quaternion power under SO(3)'s rotation** ([`0054`](./decisions/0054-the-screw-twin-is-the-oracles-translation-under-the-shipped-rotation.md)).
  The rotation is SO(3)'s routine above, to the bit, with its two arms; each translation column is
  GE.12's dual part, $x(t) = x_0 + R_0\,x_\Delta^t$ with
  $$q_{\mathrm d} = \tfrac12(0, x_\Delta)\otimes(w, \mathbf v),\quad
  \varkappa = q_{\mathrm d,w}/\lVert\mathbf v\rVert^2,\quad
  \bar m = q_{\mathrm d,\mathbf v} + \varkappa w\,\mathbf v,$$
  $$q_{\mathrm d}^t = \big(t\varpi_t q_{\mathrm d,w},\ \varpi_t\bar m - t\varkappa\cos(t\alpha)\,\mathbf v\big),\quad
  (0, x_\Delta^t) = 2\,q_{\mathrm d}^t\otimes(\cos t\alpha, \varpi_t\mathbf v)^*,$$
  both products Hamilton's. **Below `r`'s second switch** $\varpi_t = k\,t\,r$ and
  $\cos t\alpha = c$, §4's `log_ratio` and `exp_coeffs` on their short arms (GE.15(a)): no
  transcendental and no coefficient outside §4. **Above it**, and at exactly $t = 1$,
  $\varpi_t = \sin(t\alpha)/\lVert\mathbf v\rVert$ — **not** the blend's $\sin(t\alpha)/\sin\alpha$,
  which errs more in the translation on every regime measured — with $\cos t\alpha$ from the same
  `sin_cos`. At $\lVert\mathbf v\rVert^2 = 0$, $\varkappa = 0$ (the safe argument in the division):
  the terms it multiplies are $O(\alpha^2)$ (GE.15(b)), so no range constant of the scalar type
  is needed. Other $N$ keep the definition until a stratum verifies the power there (`0006`).
- **Invariance.** The SE(3) geodesic is left- and right-invariant:
  $X_0 H\,\mathrm{Exp}(t\,\mathrm{Ad}_{H^{-1}} d) = X_0\,\mathrm{Exp}(t\,d)\,H$. For
  `Product<SO3, R3>` **the answer depends on which law $a\cdot H$ uses**: under the product's own
  componentwise law it is left- *and* right-invariant, exactly; read as an SE(3) pose
  ($a\cdot H = (R_aR_H, R_a\mathbf t_H + \mathbf t_a)$) it is left-invariant only, and its
  right-invariance test is **supposed to fail** (`tf_tree` D5), by
  $\lvert\mu_s(\theta)\rvert\,\lVert\mathbf t_{H\perp}\rVert$ with
  $\lvert\mu_s\rvert = \tfrac12 s(1-s)\theta^2 + O(\theta^4)$ and
  $\mathbf t_{H\perp}\perp$ the axis (`docs/maths/geodesics.md` GE.5,
  [`0045`](./decisions/0045-two-phase-4-checks-cannot-be-taken-as-written.md)).
- **Jacobians.** $\partial X(t)/\partial t = d$ (right, body frame) on either side, and per side,
  with $d = X_1\ominus_R X_0$ and $d_L = X_1\ominus_L X_0 = \mathrm{Ad}_{X_0}d$:

  | | right | left |
  |---|---|---|
  | $\partial X(t)/\partial X_1$ | $t\,J_r(td)\,J_r^{-1}(d)$ | $t\,J_l(t\,d_L)\,J_l^{-1}(d_L)$ |
  | $\partial X(t)/\partial X_0$ | $(1-t)\,J_l\big((1-t)d\big)\,J_l^{-1}(d)$ | $(1-t)\,J_r\big((1-t)d_L\big)\,J_r^{-1}(d_L)$ |

  Each $\partial X(t)/\partial X_0$ is $\partial X(t)/\partial X_1$ of the swapped pair
  $(X_1, X_0)$ at $1-t$: swapping negates $d$, and $J_l(-x) = J_r(x)$ exchanges the pair.
  **This is the form that ships**
  ([`0043`](./decisions/0043-the-geodesic-jacobian-ships-the-cancellation-free-form.md)). The
  equal right form $\partial X(t)/\partial X_0 = \mathrm{Ad}_{\mathrm{Exp}(-td)} - t\,J_r(td)\,J_r^{-1}(d)\,\mathrm{Ad}_{\Delta^{-1}}$
  is the definition it is derived from and is **not** evaluated: it is a difference of two $O(1)$
  matrices vanishing like $1-t$, with relative error $O(u/(1-t))$ — measured $19\,u$ at
  $1-t = 10^{-1}$ and $1.3\times10^{12}\,u$ at $10^{-12}$, against $\le 0.72\,u$ for the form
  above — and it is not exactly $0$ at $t = 1$, which the boundary values require. The left
  $\partial X(t)/\partial X_0 = I - \partial X(t)/\partial X_1$ cancels the same way and is not
  evaluated either (`docs/maths/geodesics.md` GE.7).
- **Velocity.** Body velocity $d/\Delta t$, constant along the geodesic.
- **Domain.** $\theta(d) < \pi$ for a unique geodesic; at $\pi$ the problem is ill-posed and the
  result is §3.2's function of the quaternion sign.

## 11. Error metrics and conditioning

- **Forward error** $\|\hat y - y\| / \max(\|y\|, \|y\|_{\text{floor}})$ in units of $u$, norm-wise
  (vectors: 2-norm; matrices: Frobenius; quaternions: after aligning the sign of $\hat y$ with $y$).
  $\|y\|_{\text{floor}}$ is the function's scale (1 for rotations and Jacobians, $\|\rho\|$-scale for
  translations) so an exact zero is not divided by.
- **Backward error** for `Log`: $\|\mathrm{Exp}(\widehat{\mathrm{Log}}\,X) \ominus_R X\|$ in $u$. Primary
  for the near-π strata.
- **Quaternion inputs make `Log`'s forward error well-defined** everywhere except $w = +0$ exactly,
  where the metric is sign-invariant: $\min(\|\hat\varphi - \varphi\|, \|\hat\varphi + \varphi\|)$.
  Matrix inputs (`from_matrix`) are ill-conditioned near $\pi$ in the axis; their strata report
  backward error only.
- **Per stratum:** max and p99; never a mean. Non-finite outputs are counted separately and any
  non-zero count fails.
- **Bars** ([`0006`](./decisions/0006-the-instrument-comes-first.md)): domination over the best
  oracle's max; no-regress against the committed baseline, **exact** (D16).

## 12. Domains and singularities

| Routine | Domain | Outside (release) |
|---|---|---|
| `SO3::log` | all unit $q$; at $w = +0$ a function of $q$'s sign (§3.2) | — |
| `jr_inv`, `jl_inv` (SO(3), SE_N(3)) | $\theta < 2\pi$ | unspecified finite value |
| `from_wxyz_unchecked` | $\lvert\|q\|^2 - 1\rvert \le 2^{-40}$ (`f64`), $2^{-16}$ (`f32`) | garbage in, garbage out |
| `from_wxyz_normalized` | $\|q\|^2$ normal (components within $\approx 10^{\pm154}$ `f64`, $10^{\pm19}$ `f32`); `debug_assert!` on the result (§3.6, [`0027`](./decisions/0027-a-normalizing-constructor-normalizes.md)) | NaN at $q = 0$; zero above the overflow |
| `renormalize` | none; a normalization only for $\lvert\|q\|^2 - 1\rvert \le 2^{-26.29}$ (`f64`), $2^{-11.79}$ (`f32`) (§3.6) | defined everywhere; $0$ at $\|q\|^2 = 3$, $q$ reversed beyond |
| `S2Chart::local` | $m \ne -n$ | unspecified finite value |
| `geodesic` | $\theta(d) < \pi$ | §10 |
| `Sim3` | $\sigma$ finite | — |
| `Dual::sqrt` derivative | $v > 0$ | $\pm\infty$ ($d \ne 0$), NaN ($d = 0$); value unaffected (0020) |
| `Real::cbrt` | every $x$ (odd; signed zeros and infinities are their own cube roots) | — |
| `Dual::cbrt` derivative | $v \ne 0$ | $\pm\infty$ ($d \ne 0$), NaN ($d = 0$); value unaffected |
| `Dual::atan2` derivative | $x_v^2 + y_v^2$ normal (larger argument in $\approx 10^{\pm154}$ `f64`, $10^{\pm19}$ `f32`) | $\pm\infty$, NaN at the origin, or $0$; value unaffected |
| `Dual` quotient derivative | $q = a_v/b_v$ finite | NaN (value $\pm\infty$) |
| `chol` | positive definite | `(L, mask = false)`; `L` finite for every input, nothing asserted (§15.3) |
| `solve_lower`, `solve_upper`, `chol_solve` | every diagonal entry nonzero and not NaN (`debug_assert!`) | $\pm\infty$ or NaN (§15.3) |
| Strided `get`, `set`, `block` | in bounds | **panic** (the one documented class, D11) |

## 13. References

- T. D. Barfoot, *State Estimation for Robotics*, 2017; T. D. Barfoot, P. T. Furgale,
  "Associating Uncertainty With Three-Dimensional Poses for Use in Estimation Problems", T-RO 2014.
- J. Solà, J. Deray, D. Atchuthan, "A micro Lie theory for state estimation in robotics", 2018.
- H. Sommer et al., "Why and How to Avoid the Flipped Quaternion Multiplication", 2018.
- A. Barrau, S. Bonnabel, "A Mathematical Framework for IMU Error Propagation with Applications to
  Preintegration", ICRA 2020; M. Brossard, A. Barrau, P. Chauchat, S. Bonnabel, "Associating
  Uncertainty to Extended Poses for on Lie Group IMU Preintegration With Rotating Earth", T-RO 2022.
- E. Eade, "Lie Groups for Computer Vision", 2014; H. Strasdat, *Local Accuracy and Global
  Consistency for Efficient Visual SLAM*, PhD thesis, 2012.
- S. W. Shepperd, "Quaternion from rotation matrix", J. Guidance and Control, 1978.
- A. McAdams et al., "Computing the Singular Value Decomposition of 3×3 matrices with minimal
  branching and elementary floating point operations", 2011; J. Kopp, "Efficient numerical
  diagonalization of hermitian 3×3 matrices", 2008.
- H. Sommer, J. Forbes, R. Siegwart, P. Furgale, "Continuous-Time Estimation of attitude using
  B-splines on Lie groups", 2016; C. Sommer, V. Usenko, D. Schubert, N. Demmel, D. Cremers,
  "Efficient Derivative Computation for Cumulative B-Splines on Lie Groups", CVPR 2020 (Phase 7).

## 14. Reference twins

Primitives (`Exp`, `Log`, coefficients, `act`) are checked against the corpus. Composites keep a
twin in `helicoid::reference` and a proptest `<name>_matches_reference`; `cargo xtask lint` checks
this table against the code. One exception: `helicoid-linalg` has no `reference` module, so a
composite there whose twin is a composition of its public functions has no `reference` item; the
row spells the composition and the proptest writes it inline
([`0019`](./decisions/0019-a-cholesky-solve-without-the-transpose.md)).

| Fast | Reference twin | Phase |
|---|---|---|
| `SEn3::jr_inv` (dual-matrix inverse) | dense Gauss–Jordan inverse of `jr` | 3 |
| `LieGroup::jl` (`jr(-τ)`) | $\mathrm{Ad}_{\mathrm{Exp}(\tau)}\,J_r(\tau)$, dense | 3 |
| `*_jacobians` (§2.3 closed forms) | chains of primitive Jacobians, dense | 3 |
| `SO3::act_many` | per-point `act` | 3 |
| `SEn3Jac::mul`, `inverse` | dense product / Gauss–Jordan | 3 |
| `SEn3::mul_inv` ($a\,b^{-1}$, one rotation per column) | the composition `a * b.inverse()`, which it is **not** bit-identical to; it is `lminus`'s body ([`0048`](./decisions/0048-the-relative-transform-pair-earns-the-surface-dot-and-norm-do-not.md)) | 4 |
| `SEn3::inv_mul` ($b^{-1}a$, subtracting before the rotation) | the composition `b.inverse() * a`; it is `rminus`'s body ([`0048`](./decisions/0048-the-relative-transform-pair-earns-the-surface-dot-and-norm-do-not.md)) | 4 |
| `SEn3Jac::apply`, `apply_transpose` | dense $J x$ and $J^\top x$ | 3 |
| `SEn3Jac::sandwich` | dense $J\,\Sigma\,J^\top$; also the twin the `Gaussian` row below asks for | 3 |
| `ProductJac::sandwich` | dense $J\,\Sigma\,J^\top$ of the block-diagonal $J$ | 3 |
| `SO3::geodesic` (GE.14's blend, §10) | $q_0\,\mathrm{Exp}(t\,\mathrm{Log}(q_0^{*}q_1))$ — `reference::geodesic`, which is also `LieGroup::geodesic`'s **provided** body, so the `twin` leg of `laws::geodesic` is this row's proptest and reads 7.213 `u` at binary64 and 7.089 at binary32, where it read `gerr`'s 1.118 floor while the two were one expression ([`0050`](./decisions/0050-the-geodesic-s-denominator-is-the-whole-domination-gap.md)) | 4 |
| `SE3::geodesic` (dual-quaternion power, §10) | $X_0\,\mathrm{Exp}(t\,\mathrm{Log}(X_0^{-1}X_1))$ — `reference::geodesic`; `se3_geodesic_matches_reference` ($10^5$ pairs per precision, three regimes, extrapolation its own row) and the `twin` leg of `laws::geodesic`, 8.051 `u` at binary64 and 8.587 at binary32 ([`0054`](./decisions/0054-the-screw-twin-is-the-oracles-translation-under-the-shipped-rotation.md)) | 4 |
| `geodesic_jacobians` | `Dual` through the reference geodesic | 4 |
| `Gaussian::to_left` / `to_right` | dense $\mathrm{Ad}\,\Sigma\,\mathrm{Ad}^\top$: `reference::sen3jac_sandwich` with $J = \mathrm{Ad}$ | 5 |
| `gamma_apply_jacobian` | `Dual` through `reference` $\Gamma_m$ series (dense sum) | 5 |
| `S2Chart::local` | $\mathrm{Log}$ of the minimal rotation taking $n$ to $m$, projected on $B$ | 5 |
| `chol_solve` (`helicoid-linalg`) | composition of public fns, no `reference` item: `solve_upper(&l.transpose(), solve_lower(&l, b))` | 2 |

## 15. Cholesky and the triangular solves

Fixed size $N$, over any `Real`; only one triangle of the input is read, the caller owns symmetry.

### 15.1 Recurrence

For symmetric positive definite $A$ there is a unique lower-triangular $L$ with positive diagonal
and $A = L L^\top$. Since $(L L^\top)_{ij} = \sum_{k \le j} L_{ik} L_{jk}$ for $i \ge j$, solving
for the last term column by column (Cholesky–Crout) gives, for $j = 0, \dots, N-1$:

$$
d_j = a_{jj} - \sum_{k<j} L_{jk}^2, \qquad L_{jj} = \sqrt{d_j}, \qquad
L_{ij} = \frac{a_{ij} - \sum_{k<j} L_{ik} L_{jk}}{L_{jj}}\quad (i > j).
$$

Every sum runs in increasing $k$ from its first term, and the pivot is $a_{jj}$ minus the
*finished* sum, not a running subtraction: the operation sequence is fixed (D16) and pinned to the
bit by `the_factor_sums_left_to_right_and_subtracts_the_finished_sum`. The solves are
$x_i = (b_i - \sum_{k<i} l_{ik} x_k)/l_{ii}$ (forward, $i$ ascending) and
$x_i = (b_i - \sum_{k>i} u_{ik} x_k)/u_{ii}$ (backward, $i$ descending, $k$ ascending); $A x = b$
from the factor is §15.6.

### 15.2 The mask and the failure convention

With $\mathrm{fin}(v) \Leftrightarrow 0 \cdot v = 0$ (false for NaN and $\pm\infty$), a pivot passes
iff $\mathrm{ok}_j = (0 < d_j) \wedge \mathrm{fin}(d_j)$, and an entry $\hat x_{ij}$ (the quotient
above, as computed) passes iff $\mathrm{fin}(\hat x_{ij})$. The mask is
$\bigwedge_j \mathrm{ok}_j \wedge \bigwedge_{\mathrm{ok}_j,\, i>j} \mathrm{fin}(\hat x_{ij})$: a statement
about the *computed* values, not a certificate about $A$. A zero pivot ($\pm 0$) is not positive:
a singular matrix is not positive definite. Nothing branches; every lane evaluates every arm, so
the failed arms are made harmless:

- a failed pivot stores $L_{jj} = 1$ ($\sqrt{\cdot}$ is fed $\mathrm{select}(\mathrm{ok}_j, d_j, 1)$,
  so a `Dual` derivative stays finite) and $L_{ij} = +0$ for every $i > j$, so garbage does not
  feed later columns (it squares per column: $M + M^\top$ with entries below 20 overflows `f32`
  at $N = 6$);
- an entry that is not finite (overflow, or NaN or $\infty$ read from $A$) is stored as $+0$.

Hence $L$ is finite for every input, entries above the diagonal are $+0$, and $\operatorname{diag}
L > 0$. Overflow is reported through the mask, not returned as $\infty$; the magnitude range of
$A$ is therefore not a domain of `chol`.

### 15.3 Domains

`chol` accepts every input, so it asserts nothing (`API.md` R4, R6: the mask is the report); its
domain is "positive definite", on which $L$ is the Cholesky factor. `solve_lower`, `solve_upper`
and `chol_solve` need every diagonal entry nonzero and not NaN, `debug_assert!`ed; a release build
divides and returns $\pm\infty$ or NaN.

### 15.4 Error bounds (Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed.)

With $\gamma_k = k u/(1 - k u)$ and $u$ as in §2.1, when the mask is set (the factorization ran to
completion):

- **Thm 10.3.** $\hat L \hat L^\top = A + E$, $\lvert E \rvert \le \gamma_{N+1}
  \lvert \hat L \rvert \lvert \hat L \rvert^\top$ entrywise, however ill-conditioned $A$ is.
- **Thm 8.5.** A triangular solve returns $(T + E)\hat x = b$, $\lvert E \rvert \le \gamma_N
  \lvert T \rvert$.
- **Thm 10.4.** Solving $A x = b$ through $\hat L$ returns $(A + E)\hat x = b$,
  $\lvert E \rvert \le \gamma_{3N+1} \lvert \hat L \rvert \lvert \hat L \rvert^\top$.

The converse fails: a positive definite $A$ with $\mathrm{cond}(A)$ near $1/u$, or a rank-deficient
positive semidefinite one, has a pivot that is zero or tiny in exact arithmetic and is rounding
noise as computed, so the mask may go either way there. A clear mask means "a computed pivot was
not positive, or an entry overflowed". The tests check each bound with its own product's roundings
added (`chol_tests`, header).

### 15.5 Forward mode

Through `Dual` the recurrence differentiates itself. It agrees with
$\mathrm{d}L = L\,\Phi(L^{-1}\,\mathrm{d}A\,L^{-\top})$, $\Phi(X)$ the lower triangle of $X$ with
the diagonal halved (checked by hand at $A = \left[\begin{smallmatrix}4&2\\2&5\end{smallmatrix}\right]$
in `the_dual_derivative_matches_a_hand_derivation`). Derivative lanes of a failed arm are those of
the constants $1$ and $0$; outside the domain of `Dual::sqrt` and the quotient (§12) they follow
those rows.

### 15.6 Solving $A x = b$ from the factor

$A x = b$ is $L y = b$ (forward solve, §15.1) then $L^\top x = y$. With $u = L^\top$, $u_{ik} =
l_{ki}$, the second is the backward solve of §15.1 read from $L$ directly:

$$
x_i = \frac{y_i - \sum_{k>i} l_{ki}\, x_k}{l_{ii}}, \qquad i = N-1, \dots, 0,
$$

$k$ ascending from $i+1$, the sum accumulated from its first term and $+0$ when empty. $l_{ki}$ for
$k > i$ is column $i$ below its diagonal, so `chol_solve` reads $L$ in storage order and forms no
$L^\top$. It is defined as, and tested against (`chol_solve_matches_reference`), the composition
`solve_upper(&l.transpose(), solve_lower(&l, b))`, its reference twin (§14).

**Bit identity.** `transpose` moves entries and computes nothing, so $u_{ik}$ and $l_{ki}$ are the
same bits. Both routines then run the same sequence of `S` operations on the same operands: the
products $l_{ki} \cdot x_k$ in that operand order, the left fold of `sum`, the subtraction from
$y_i$, the division by $l_{ii}$, and the same `debug_assert!` on the same diagonal. The output is
therefore equal to the bit in every lane of every `S` (`f64`, `f32`, `Dual`), $\pm 0$ and
$\pm\infty$ included, and so is the release behaviour at a zero or NaN diagonal. The one exception
is the sign and payload of a NaN produced by arithmetic, which Rust leaves unspecified and a
release build may swap by commuting operands (`PHASE2.md` §3): a NaN is a NaN in both. The error
bound is that of the composition (Thm 10.4, §15.4).
