# Ambient Jacobians: the quaternion and SE(3) charts in storage coordinates

> Non-normative companion to [`PHASE6.md`](../PHASE6.md) §1 (`AmbientChart`), [`NUMERICS.md`](../NUMERICS.md) §1 (the quaternion,
> SE_N(3)), §3.1, §3.2, §3.6 (`Exp`, `Log`, renormalization), §12 (`from_wxyz_unchecked`) and
> [`0012`](../decisions/0012-a-retraction-is-a-chart.md). **`NUMERICS.md` wins on any conflict; a disagreement is an open item in the
> [maths index](./index.md)**, which holds the notation and the `Checked:` convention. The charts: [`charts.md`](./charts.md); the
> quaternion: [`so3.md`](./so3.md); $J_l$ and the sides: [`lie-groups.md`](./lie-groups.md).

$q = (w, v)$, $v = (x, y, z)$, is a stored quaternion (as in `PHASE6.md` §1; `so3.md` writes $u$ for the vector part), $\eta = \lVert q\rVert^2 - 1$,
$\hat q = q/\lVert q\rVert$, $\theta$ the rotation angle of $\hat q$ (so $\lVert v\rVert = \lVert q\rVert\sin\frac\theta2$; the norm of a tangent is always written $\lVert\delta\rVert$, $\lVert\varphi\rVert$) and $R = R(\hat q)$. $Q_L(q)$, $Q_R(q)$ are the
$4\times4$ matrices of $p \mapsto q\otimes p$ and $p \mapsto p\otimes q$, $E = \begin{bmatrix}0\\ I_3\end{bmatrix}$ ($4\times3$, the pure quaternions), $P^+$ the Moore–Penrose inverse.
$P$, $M$ are `PlusJacobian` and `MinusJacobian` of the **right** chart ($q\otimes\mathrm{Exp}\,\delta$); $P^L$, $M^L$ those of the left chart; $P_C$, $M_C$, $\delta_C$ Ceres'.
The ambient dimension is $m$ ($4$, or $7$ for $(q, t)$), the tangent dimension $d$ ($3$, or $6$); $n$ counts retractions (AJ.5(b)).

## Results

| Label | Result | `NUMERICS.md`; `PHASE6.md` (P6) | Implemented by |
|---|---|---|---|
| AJ.1–AJ.2 | the ambient data of a chart, $P_x = \mathrm D\,\mathrm{ret}_x(0)$ and $M_x = \mathrm D\,\mathrm{loc}_x(x)$; $MP = I$ whatever the extension of $\mathrm{loc}$; $PM$ projects onto $T_x\mathcal M$, orthogonally iff $M = P^+$; $\mathrm D(f\circ\mathrm{ret}) = \mathrm Df\,P$ | P6 §1 | `AmbientChart` |
| AJ.3 | $P = \frac12Q_L(q)E$ and $M = 2E^\top Q_L(q)^\top = 4P^\top$, derived from the Hamilton product, first-order $\mathrm{Exp}$ and the scale invariance of $\mathrm{Log}$; the retraction to second order | §1, §3.1, §3.2; P6 §1 | `write_plus_jacobian`, `write_minus_jacobian` (quaternion) |
| AJ.4 | $MP = \lVert q\rVert^2I_3$, $PM = \lVert q\rVert^2I_4 - qq^\top$; for unit $q$, $MP = I$ and $PM$ the orthogonal projector onto $T_qS^3$; singular values $\frac12$ and $2$ | P6 §1 | (checks) |
| AJ.5 | non-unit $q$: $P$ exact, $\mathrm D\,\mathrm{loc} = P^+ = M/\lVert q\rVert^2$, $MP - I = \eta I$; the drift of $\eta$; three readings at a stored non-unit $q$; the double cover | §3.6, §12 | `renormalize`, `Quat::from_wxyz_unchecked`; the corpus ids of P6 §1 |
| AJ.6 | SE(3) as $(q, t)$: $P_{SE3} = \mathrm{diag}(P, R)$, $M_{SE3} = \mathrm{diag}(M, R^\top)$ for `Screw` and `Decoupled`, proved; the dependence on $R$; `WorldTranslation` and the left chart | §5.1; P5 §1.3; P6 §1 | `write_plus_jacobian`, `write_minus_jacobian` (`Screw`, `Decoupled`) |
| AJ.7 | Ceres' `QuaternionManifold` is the left chart in half-angle coordinates, exactly: $P_C = 2PR^\top$, $M_C = \frac12RM$; what is the same and what differs (side, factor $2$, no sign flip) | P6 §1, §7 | the container runner of P6 §7 |
| AJ.8 | GTSAM retracts on the storage in the right chart of the tangent and has no ambient Jacobian | P6 §1 | (`Chart::retract_jacobian` serves it) |
| AJ.9 | conditioning and rounding: $P$, $M$ are exact, $\kappa_2 = 1$ ($2$ for SE(3)) at every $q$; the product $J_{\mathrm{amb}}P$ | §2.1, §3.6 | — |
| AJ.10 | sign audit: the usual wrong variants and their exact gaps | — | — |

## 1. Storage larger than the tangent

A solver that stores $\mathrm{SO}(3)$ as a quaternion, or $\mathrm{SE}(3)$ as $(q, t)$, differentiates its residuals with respect to the stored numbers and moves them by a chart.
Two matrices connect the coordinates.

**Definition AJ.1 (ambient data of a chart).** Let $\mathcal M \subset \mathbb R^m$ be an embedded $d$-manifold (the storage constraint, $\lVert q\rVert = 1$), $(\mathrm{ret}_x, \mathrm{loc}_x)$ a chart at
$x \in \mathcal M$ (CH.1) with $\mathrm{ret}_x$ valued in $\mathbb R^m$, and $\mathrm{loc}_x$ defined and differentiable at $x$ on a neighbourhood of $x$ in $\mathbb R^m$: the code evaluates it at stored values
that are not exactly on $\mathcal M$ (for the quaternion, $\mathrm{Log}$ is defined for every $q \ne 0$, `NUMERICS.md` §3.2). With $T_x = T_x\mathcal M$ and $N_x = T_x^\perp$,

$$
P_x = \mathrm D\,\mathrm{ret}_x(0) \in \mathbb R^{m\times d}\ (\text{`PlusJacobian`}),\qquad M_x = \mathrm D\,\mathrm{loc}_x(x) \in \mathbb R^{d\times m}\ (\text{`MinusJacobian`}),
$$

in ambient coordinates: Ceres' $D_2\mathrm{Plus}(x, 0)$ and $D_1\mathrm{Minus}(x, x)$. Both depend on $x$ alone, like the chart frozen at $x$ (CH.1), so `write_plus_jacobian(&self, ..)` of `AmbientChart` reads the base
it was frozen at.

**Proposition AJ.2.**

- (a) $M_xP_x = I_d$, for every extension of $\mathrm{loc}_x$ off $\mathcal M$.
- (b) $P_x$ is injective with range $T_x$. $P_xM_x$ is the projector onto $T_x$ along $\ker M_x$; it is symmetric (an orthogonal projector) iff $\ker M_x = N_x$ iff $M_x = P_x^+$, and
  $M_x$ is then the unique left inverse of $P_x$ that vanishes on $N_x$. Every left inverse is $P_x^+ + Z(I - P_xP_x^+)$: **$MP = I$ does not determine $M$**.
- (c) For $f$ differentiable at $x$, $\mathrm D(f\circ\mathrm{ret}_x)(0) = \mathrm Df(x)\,P_x$. Two extensions of $f\vert_{\mathcal M}$ have derivatives that differ by a matrix vanishing on $T_x$, so the product does not depend on
  the extension although $\mathrm Df(x)$ does.
- (d) A tangent covariance $\Sigma$ maps to the ambient $P_x\Sigma P_x^\top$ (rank $d$) and back by $M_x(\cdot)M_x^\top$.

*Proof.* (a) Differentiate $\mathrm{loc}_x(\mathrm{ret}_x(\delta)) = \delta$ (CH.1, (Inv)) at $0$: the curve $\mathrm{ret}_x(t\delta)$ lies in $\mathcal M$, so only $\mathrm D\,\mathrm{loc}_x$ along $T_x$ enters. (b) $\mathrm{ret}_x$ maps into
$\mathcal M$, so $P_x\delta \in T_x$; injective by (a); $\dim T_x = d$. $(PM)^2 = P(MP)M = PM$ with range $\mathrm{range}\,P$ and kernel $\ker M$. A symmetric idempotent is an orthogonal projector, of kernel $N_x$;
$P^+ = (P^\top P)^{-1}P^\top$ has kernel $\ker P^\top = N_x$, and a left inverse vanishing on $N_x$ is determined on $T_x \oplus N_x = \mathbb R^m$. (c) The chain rule. (d) By (a). $\square$

*Example (c).* $f(q) = R_\square(q)a$ and the vector form $a + 2w\,v\times a + 2\,v\times(v\times a)$, $a \in \mathbb R^3$, agree on the sphere and differ by $(\lVert q\rVert^2 - 1)a$ off it (SO.2(d)), so their $\mathrm Df$ differ by $2aq^\top$;
$q^\top P = 0$ (AJ.4) removes it, and both give $\mathrm Df\,P = -R[a]_\times$, the right row of `NUMERICS.md` §2.4 (SO.10). A solver may differentiate whichever extension: $P$ discards the radial derivative (in exact arithmetic; for the rounding, AJ.9(c)).

**Checked:** AJ.2 is exercised in the instances AJ.3–AJ.6. The example: mpmath 1.4.1, 100 digits, 30 unit $q$ and random $a$, central differences ($h = 10^{-32}$) of both forms: $\mathrm Df_\square - \mathrm Df_{\mathrm{vec}} - 2aq^\top$
$\le 2.5\times10^{-69}$, $(\mathrm Df_\square - \mathrm Df_{\mathrm{vec}})P \le 1.3\times10^{-69}$, and $\mathrm Df\,P + R[a]_\times \le 1.1\times10^{-69}$ for both. Script not committed. **Permanent:** none specified for AJ.2;
`act_jacobians` (`NUMERICS.md` §2.4) is the tangent-side check of the example.

## 2. The quaternion

Storage $q \in \mathbb R^4\setminus\{0\}$, $\mathrm{ret}_q(\delta) = q\otimes\mathrm{Exp}(\delta)$, $\mathrm{loc}_q(q') = \mathrm{Log}(q^{-1}\otimes q')$: `RightChart<SO3>` (CH.3), with $\mathrm{Exp}\,\delta = (\cos\frac{\lVert\delta\rVert}2, k\delta)$ (`NUMERICS.md` §3.1)
and the $\mathrm{Log}$ of §3.2. That $\mathrm{Log}$ is defined for every $q' \ne 0$ and invariant under $q' \mapsto sq'$, $s > 0$ ($\mathrm{atan2}$ and $u/n$ are); this is the extension of $\mathrm{loc}$ off the sphere, and it
fixes $M$ on $N_q = \mathrm{span}\{q\}$ (AJ.2(b)).

**Proposition AJ.3 ($P$, $M$).** For every $q \ne 0$, with $Q_L(q) = \begin{bmatrix}w & -v^\top\\ v & wI + [v]_\times\end{bmatrix}$ (SO.1):

- (a) $\mathrm D\,\mathrm{ret}_q(0) = P = \tfrac12Q_L(q)E = \tfrac12\begin{bmatrix}-v^\top\\ wI + [v]_\times\end{bmatrix} = \tfrac12\begin{bmatrix}-x & -y & -z\\ w & -z & y\\ z & w & -x\\ -y & x & w\end{bmatrix}$.
- (b) $\mathrm D\,\mathrm{loc}_q(q) = M/\lVert q\rVert^2$, with $M = 2E^\top Q_L(q)^\top = 2\big[-v \;\big\vert\; wI - [v]_\times\big] = 2\begin{bmatrix}-x & w & z & -y\\ -y & -z & w & x\\ -z & y & -x & w\end{bmatrix}$.
- (c) $M = 4P^\top$.
- (d) $\mathrm{ret}_q(\delta) = q + P\delta - \tfrac18\lVert\delta\rVert^2q + O(\lVert\delta\rVert^3)$.

These are the matrices of `PHASE6.md` §1.

*Proof.* (a), (d) $\cos\frac{\lVert\delta\rVert}2 = 1 - \frac{\lVert\delta\rVert^2}8 + O(\lVert\delta\rVert^4)$ and $k = \frac12 - \frac{\lVert\delta\rVert^2}{48} + O(\lVert\delta\rVert^4)$ (`NUMERICS.md` §4), so $\mathrm{Exp}\,\delta = (1, \frac12\delta) - \frac{\lVert\delta\rVert^2}8(1, 0) + O(\lVert\delta\rVert^3)$.
$p \mapsto q\otimes p = Q_L(q)p$ is linear (SO.1: the first row is $(w, -v^\top)$, the vector part is $p_wv + wp_v + [v]_\times p_v$), so $q\otimes\mathrm{Exp}\,\delta = q + \frac12Q_L(q)E\delta - \frac{\lVert\delta\rVert^2}8q + O(\lVert\delta\rVert^3)$.
(b) $q^*\otimes q'$ is a positive multiple ($\lVert q\rVert^2$) of $q^{-1}\otimes q'$, so $\mathrm{loc}_q(q') = \mathrm{Log}(q^*\otimes q')$. At $q' = q$, $q^*\otimes q = (\lVert q\rVert^2, 0)$, and
$q^*\otimes(q + \mathrm dq') = (\lVert q\rVert^2, 0) + Q_L(q)^\top\mathrm dq'$ ($Q_L(q^*) = Q_L(q)^\top$). For $(s, a)$ with $s > 0$, $\mathrm{Log} = 2\,\mathrm{atan2}(\lVert a\rVert, s)\,a/\lVert a\rVert = \frac{2a}s\,(1 + O(\lVert a\rVert^2/s^2))$
(the series of $r$, `NUMERICS.md` §4; at $a = 0$, where $a/\lVert a\rVert$ is undefined, $\mathrm{Log} = r(\lVert a\rVert^2, s)\,a$ with the smooth kernel $r$, and that is what is differentiated), so $\mathrm{loc}_q(q + \mathrm dq') = \frac2{\lVert q\rVert^2}E^\top Q_L(q)^\top\mathrm dq' + O(\lVert\mathrm dq'\rVert^2)$; the rows $2$–$4$ of $Q_L(q)^\top = \begin{bmatrix}w & v^\top\\ -v & wI - [v]_\times\end{bmatrix}$
are $[-v \mid wI - [v]_\times]$. (c) $P^\top = \frac12E^\top Q_L(q)^\top$. $\square$

So $P$ sends a rotation vector to the quaternion increment $\frac12q\otimes(0, \delta)$: the **factor $\frac12$** is the half angle of $q = (\cos\frac\theta2, \sin\frac\theta2\,\hat n)$, and the factor $2$ of $M$ is its inverse. For unit $q$, $q + P\delta$ has
$\lVert\cdot\rVert^2 = 1 + \frac14\lVert\delta\rVert^2$ (AJ.4) and the term $-\frac18\lVert\delta\rVert^2q$ of (d) brings it back; $P$ is first-order data and does not identify the retraction (AJ.10).

**Checked:** mpmath 1.4.1, 100 digits, 45 unit $q$ (40 random; $w = 0$; $q = 1$; $\lVert v\rVert = 2\times10^{-40}$; $w = 10^{-30}$; $-q_0$). $\mathrm{Exp}$ by the quaternion series $\sum p^m/m!$ (SO.4(c)), against
`mp.expm` of $W$: $5.7\times10^{-101}$. $P$ by central differences ($h = 10^{-32}$) of $\delta \mapsto q\otimes\mathrm{Exp}\,\delta$ against (a): $2.1\times10^{-66}$ (65 digits); $M$ by central differences of $q' \mapsto \mathrm{Log}(q^{-1}\otimes q')$,
the relative quaternion projected on the group and $\mathrm{Log}$ by `mp.logm` of its matrix (the definition; near the identity, where `mp.logm` is principal, index open items), against (b): $6.7\times10^{-65}$. `mp.diff` (the reference of `PHASE6.md` §1,
default step, 60 digits) on 6 samples: $P$ equal to the last digit (5 entries each), $M$ (4 entries each) $\le 1.6\times10^{-61}$. (c): exactly $0$. (d): $\lVert\mathrm{ret} - (q + P\delta - \frac18\lVert\delta\rVert^2q)\rVert/s^3 \le 1.2$ at
$\delta = s\delta_0$, $s = 10^{-12}$, 30 samples. Exact, sympy 1.14.0: (a) as the Jacobian of $q\otimes(1, \delta/2)$ and (b) as that of $2\,\mathrm{vec}(q^*\otimes q')/\mathrm{sc}(q^*\otimes q')$ at $q' = q$ (equal to $M/\lVert q\rVert^2$),
as polynomial identities in $w, x, y, z$. Scripts not committed. **Permanent:** planned, corpus `so3_plus_jacobian` (`PHASE6.md` §1; $P$ only); **none specified for $M$** ($M = 4P^\top$ is (c); `PHASE6.md` §1 names no minus id, and
$MP = I$ does not pin $M$, AJ.2(b); index, open items).

**Proposition AJ.4 ($MP$ and $PM$).** For every $q$:
$$
MP = \lVert q\rVert^2I_3,\quad PM = \lVert q\rVert^2I_4 - qq^\top,\quad Mq = 0,\quad q^\top P = 0,\quad P^\top P = \tfrac{\lVert q\rVert^2}4I_3,\quad MM^\top = 4\lVert q\rVert^2I_3,\quad (PM)^2 = \lVert q\rVert^2PM.
$$
Hence for unit $q$: $MP = I_3$; $PM = I_4 - qq^\top$ is the orthogonal projector onto $q^\perp = T_qS^3$ (symmetric, idempotent, rank $3$, trace $3$, kernel $\mathrm{span}\{q\}$); the singular values of $P$ are $\frac12$ (three times) and of $M$ are $2$;
$M = P^+ = 4P^\top$, the unique left inverse of $P$ with $Mq = 0$ (AJ.2(b)).

*Proof.* $Q_L(q)^\top Q_L(q) = Q_L(q)Q_L(q)^\top = \lVert q\rVert^2I_4$ ($\lVert q\otimes p\rVert = \lVert q\rVert\lVert p\rVert$, polarized) and $Q_L(q)e_w = q$, $E^\top E = I_3$, $E^\top e_w = 0$, $EE^\top = I_4 - e_we_w^\top$. Then
$MP = E^\top Q_L^\top Q_LE = \lVert q\rVert^2I_3$; $PM = Q_LEE^\top Q_L^\top = Q_L(I - e_we_w^\top)Q_L^\top = \lVert q\rVert^2I_4 - qq^\top$; $Mq = 2E^\top Q_L^\top Q_Le_w = 2\lVert q\rVert^2E^\top e_w = 0$ and $q^\top P = \frac12e_w^\top Q_L^\top Q_LE = 0$;
the rest by (c) of AJ.3 and $P^\top P = \frac14MP$. For unit $q$, $P^+ = (P^\top P)^{-1}P^\top = 4P^\top = M$. $\square$

**Checked:** as AJ.3 (45 unit $q$, 100 digits): $MP - I$ $1.4\times10^{-101}$; $PM - (I - qq^\top)$ $2.1\times10^{-101}$; $PM$ symmetric $0$ and $(PM)^2 - PM$ $2.1\times10^{-101}$; $Mq$ and $q^\top P$ exactly $0$; $P^\top P - \frac14I$ $3.6\times10^{-102}$; $MM^\top - 4I$
$5.7\times10^{-101}$; singular values (`mp.svd_r`) against $\frac12$ and $2$: $2.1\times10^{-101}$, $8.6\times10^{-101}$. Exact, sympy: every identity above, symbolic in $w, x, y, z$ and without $\lVert q\rVert = 1$. Script not committed. **Permanent:** none specified beyond
`so3_plus_jacobian` ($P$); the tests $MP = I$ and $PM = I - qq^\top$ do not exist yet (`AmbientChart` is not started).

**Proposition AJ.5 (off the sphere, and the double cover).** Let $q \ne 0$, $\eta = \lVert q\rVert^2 - 1$.

- (a) $P(q)$ is the exact $\mathrm D\,\mathrm{ret}_q(0)$ for every $q$, and $\mathrm D\,\mathrm{loc}_q(q) = M/\lVert q\rVert^2 = P^+$. The printed $M$ is $\mathrm D\,\mathrm{loc}_q(q)$ iff $\lVert q\rVert = 1$; in general
  $\lVert MP - I_3\rVert_2 = \lvert\eta\rvert$, and $PM$ has the eigenvalues $\lVert q\rVert^2$ (three times) and $0$: a projector iff $\eta = 0$. This is exact arithmetic on $q \ne 0$; in floating point the domain is that of `from_wxyz_unchecked`,
  $\lvert\eta\rvert \le 2^{-40}$ (`f64`; $2^{-16}$, `f32`; `NUMERICS.md` §12). For a tiny or huge $q$ ($\lVert q\rVert^2$ underflowing or overflowing) the readings of (d) that divide by $\lVert q\rVert^2$ or $\lVert q\rVert$ are not defined by these results.
- (b) $\mathrm{ret}_q(\delta)$ has exactly the norm of $q$; in floating point every retraction rounds and $\eta$ drifts. In a model (`f64`, start $q = 1$, no renormalization, $n$ successive retractions) there are two regimes, set by whether the rounding errors cancel.
  *Independent steps* (random axis, angle in $[0.01, 0.1]$): a random walk, $\lvert\eta\rvert_{\mathrm{rms}} \approx 1.1\sqrt n\,u$ (a fit: the ratio to $\sqrt n$ is $1.25$, $1.17$, $1.15$, $1.08$ at $n = 10^2$ to $10^5$), the maximum over $2000$ chains being $3.4$ times that at $n = 10^5$.
  *One repeated step* (a constant rate: fixed axis and angle): a systematic bias, $\eta$ **linear** in $n$, about $0.7u$ per retraction, of a size and a sign that depend on the angle ($+0.78$, $-0.05$, $+0.71$, $+0.65$, $+0.75$ at step angles $10^{-3}$, $10^{-2}$, $0.05$, $0.1$, $1$; mean over $200$ axes).
  At four of these five angles $\lvert\eta\rvert$ reaches the bound of `from_wxyz_unchecked`, $2^{-40} = 8192u$, after $\approx10^4$ retractions; the random walk would need $\sim5\times10^7$ (an extrapolation of the $\sqrt n$ law, whose runs stop at $10^5$). In `f32` (unit roundoff $u_{32}$) the repeated
  step biases by $0.3$ to $0.8u_{32}$ per retraction (angles $10^{-2}$ to $1$) and the bound $2^{-16} = 256u_{32}$ is met after $3\times10^2$ to $9\times10^2$ retractions.
- (c) Within that domain the pair is a left-inverse pair only to $\lvert\eta\rvert \le 2^{-40}$ ($8192u$), not to rounding. `renormalize` sends $\eta \mapsto -\frac34\eta^2 + \frac14\eta^3$ exactly (SO.14), below $u$ from $2^{-40}$.
- (d) *Three readings of "the Jacobian at a stored non-unit $q$".* A: the printed formulas, $P(q)$ and $M(q)$. B: the derivative of the definition at $q$, $P(q)$ and $M(q)/\lVert q\rVert^2$. C: the derivative at the projected point
  $\hat q$ (the reading of `so3_act`, `PHASE1.md` §4.3: $R(q/\lVert q\rVert)$), $P(q)/\lVert q\rVert$ and $M(q)/\lVert q\rVert$. A and B agree for $P$; $M$ differs between A and B by $\lvert\eta\rvert$, and every other pair by $\lvert\eta\rvert/2$ (to first order): at the `q:nonunit` stratum,
  $\lvert\eta\rvert = 2^{-45}$ (`PHASE1.md` §4.4), $256u$ and $128u$ of pure disagreement about the reference.
- (e) *Double cover.* $P(-q) = -P(q)$ and $M(-q) = -M(q)$, each the derivative at the stored representative ($\mathrm{Log}$ flips, so $\mathrm{loc}_q(-q') = \mathrm{loc}_q(q')$); $PM$ is even. For a cost even in $q$
  (any function of $R(q)$) $\mathrm Df$ is odd and $\mathrm Df\,P$ is unchanged. Mixing representatives gives $M(-q)P(q) = -I$.

*Proof.* (a) AJ.3, AJ.4; $\lVert MP - I\rVert_2 = \lvert\eta\rvert$ since $MP - I = \eta I$; $PM = \lVert q\rVert^2(I - \hat q\hat q^\top)$. (b) Measured, in a model: not a bound. (c) SO.14. (d) Scale $q \mapsto \hat q$: $P$ is linear and $M$ homogeneous of degree $1$ in $q$. (e) $Q_L(-q) = -Q_L(q)$. $\square$

**Checked:** mpmath 1.4.1, 100 digits, $q = \lambda\hat q_0$, $\lambda \in \{0.7,\ 1 \pm 2^{-45},\ 1.3,\ 3\}$, 8 unit $q_0$ each: $P$ by differentiating $q\otimes\mathrm{Exp}\,\delta$ against (a) $5.7\times10^{-66}$; $\mathrm D\,\mathrm{loc}$ (`mp.logm` of the projected relative
quaternion) against $M/\lVert q\rVert^2$ $1.8\times10^{-64}$ and against $P^+ = (P^\top P)^{-1}P^\top$ $1.8\times10^{-64}$; $MP - \lVert q\rVert^2I$ and $PM - (\lVert q\rVert^2I - qq^\top)$ $1.1\times10^{-100}$; $P^+ - M/\lVert q\rVert^2$ $2.9\times10^{-101}$; $P(\hat q) - P(q)/\lVert q\rVert$
$3.6\times10^{-102}$. Oddness (e): differentiating at $-q$, the 45 unit samples of AJ.3: $2.1\times10^{-66}$, $6.7\times10^{-65}$. (b): numpy 2.5.3, `default_rng(20260930)`, $q = 1$ at the start, $\eta$ of the final $q$ exact in rationals from the doubles. Random steps ($2000$ chains, axis normal then normalized, angle uniform), $n = 10^2$, $10^4$, $10^5$: rms $12.5u$, $115u$, $343u$ (max $47u$, $455u$, $1173u$);
in `f32`, $n \le 10^4$, rms $1.1$ to $1.25\sqrt n\,u_{32}$ likewise. One repeated step ($200$ axes per angle): mean $\eta/(nu)$ at $n = 10^2$, $10^3$, $10^4$ is $0.68$, $0.70$, $0.71$ at step angle $0.05$ (the other angles, $n = 10^4$: the values of (b); at angle $1$ the axes spread from $-7\times10^2u$ to $1.5\times10^4u$);
in `f32`, $n = 10^4$, $-0.58$, $-0.79$, $+0.30$, $-0.36$ per retraction at angles $10^{-2}$, $0.05$, $0.1$, $1$ ($10^{-3}$: erratic). The model is the Hamilton product with terms summed left to right and $(\cos\frac\alpha2, \sin\frac\alpha2\,\hat n)$ ($\alpha$ the step angle) from platform `sin`, `cos` at `f64` (rounded for `f32`), not the `helicoid` code:
the constants and the signs are those of this model. Scripts not committed. **Permanent:** none specified; (d) is a reading the corpus generator must fix
before a `so3_plus_jacobian` stratum at non-unit $q$ (index, open items).

## 3. SE(3): $(q, t)$

Storage $x = (q, t) \in \mathbb R^4\times\mathbb R^3$ with unit $q$, $R = R(q)$ ($X = (R, t)$, `NUMERICS.md` §1), tangent $\delta = [\varphi; \rho]$, ambient dimension $7$, $d = 6$ (`PHASE6.md` §1). With `PHASE5.md` §1.3 and CH.4:

| Chart | $\mathrm{ret}_x(\delta)$ | $\mathrm{loc}_x(y)$, $y = (q_Y, t_Y)$, $\varphi_Y = \mathrm{Log}(q^{-1}\otimes q_Y)$ |
|---|---|---|
| `Screw` | $(q\otimes\mathrm{Exp}\,\varphi,\ t + RJ_l(\varphi)\rho)$ | $(\varphi_Y,\ J_l^{-1}(\varphi_Y)R^\top(t_Y - t))$ |
| `Decoupled` | $(q\otimes\mathrm{Exp}\,\varphi,\ t + R\rho)$ | $(\varphi_Y,\ R^\top(t_Y - t))$ |
| `WorldTranslation` | $(q\otimes\mathrm{Exp}\,\varphi,\ t + \rho)$ | $(\varphi_Y,\ t_Y - t)$ |

**Proposition AJ.6 (SE(3)).**

- (a) For `Screw` and for `Decoupled`, $P_{SE3} = \begin{bmatrix}P & 0\\ 0 & R\end{bmatrix}$ ($7\times6$) and $M_{SE3} = \begin{bmatrix}M & 0\\ 0 & R^\top\end{bmatrix}$ ($6\times7$).
- (b) $M_{SE3}P_{SE3} = I_6$ and $P_{SE3}M_{SE3} = \mathrm{diag}(I_4 - qq^\top, I_3)$, the orthogonal projector onto $T_x(S^3\times\mathbb R^3)$; $M_{SE3} = P_{SE3}^+$. The singular values of $P_{SE3}$ are $\frac12$ (three) and $1$
  (three), those of $M_{SE3}$ are $2$ and $1$: $\kappa_2 = 2$ for both.
- (c) *They agree to first order in storage coordinates.* $\mathrm{ret}^{\mathrm{Scr}}_x(\delta) - \mathrm{ret}^{\mathrm{Dec}}_x(\delta) = (0,\ R(J_l(\varphi) - I)\rho)$, of norm $\ell(\lVert\varphi\rVert)\lVert\rho_\perp\rVert \le \frac{\lVert\varphi\rVert}2\lVert\rho\rVert \le \frac14\lVert\delta\rVert^2$
  (CH.7(c); $\lVert\varphi\rVert\lVert\rho\rVert \le \frac12\lVert\delta\rVert^2$), and to leading order $\frac12R(\varphi\times\rho)$; likewise $\mathrm{loc}^{\mathrm{Scr}} - \mathrm{loc}^{\mathrm{Dec}} = (0,\ (J_l^{-1}(\varphi_Y) - I)R^\top(t_Y - t)) = O(\lVert y - x\rVert^2)$.
  So $P^{\mathrm{Scr}} = P^{\mathrm{Dec}}$ and $M^{\mathrm{Scr}} = M^{\mathrm{Dec}}$.
- (d) *Dependence on $R$.* Only the translation block: $\rho$ is a body-frame increment and $t$ a world-frame coordinate, so $\rho \mapsto R\rho$ and $t_Y - t \mapsto R^\top(t_Y - t)$. Both matrices depend on $q$ alone (through $P(q)$ and $R(q)$),
  not on $t$, and are block diagonal: the coupling of `Screw`, $R\,\partial_\varphi(J_l(\varphi)\rho)$, vanishes at $\delta = 0$ because $\rho = 0$ there (away from it, it is $R\big(Q(\rho, \varphi) - [J_l(\varphi)\rho]_\times J_l(\varphi)\big)$, SE.9(c)).
- (e) `WorldTranslation`: $P_{SE3} = \mathrm{diag}(P, I)$ and $M_{SE3} = \mathrm{diag}(M, I)$, which equal `Decoupled`'s only if $R = I$ (CH.5(c)); `PHASE6.md` §1 specifies `Screw` and `Decoupled`. The left chart, $\mathrm{Exp}(\delta)X = (\mathrm{Exp}\,\varphi\,R,\ \mathrm{Exp}\,\varphi\,t + J_l(\varphi)\rho)$,
  has $P^L_{SE3} = \begin{bmatrix}PR^\top & 0\\ -[t]_\times & I\end{bmatrix}$, which depends on $t$.

*Proof.* (a) Differentiate the table at $\delta = 0$. Rotation: $\partial_\varphi(q\otimes\mathrm{Exp}\,\varphi) = P$ (AJ.3) and it does not depend on $\rho$. Translation: $\partial_\rho = RJ_l(0) = R$; $\partial_\varphi(t + RJ_l(\varphi)\rho)$ at $\rho = 0$ is $0$ because the map is $0$
identically on $\{\rho = 0\}$; for `Decoupled` replace $J_l$ by $I$. $\mathrm{loc}$: at $y = x$, $\varphi_Y = 0$, $\partial\varphi_Y/\partial q_Y = M$ (AJ.3(b), unit $q$), $\partial\varphi_Y/\partial t_Y = 0$; in the translation part
$J_l^{-1}(0) = I$ and $R^\top(t_Y - t) = 0$ at $y = x$, so its derivative is $R^\top$ in $t_Y$ and $0$ in $q_Y$. (b) AJ.4, $R^\top R = I$; $P_{SE3}^\top P_{SE3} = \mathrm{diag}(\frac14I, I)$, so $P_{SE3}^+ = \mathrm{diag}(4P^\top, R^\top) = M_{SE3}$.
(c) The differences are the ones of CH.4 (the table above), with $J_l(\varphi) - I = aW + bW^2 = O(\lVert\varphi\rVert)$ and $J_l^{-1}(\varphi_Y) - I = O(\lVert\varphi_Y\rVert)$ smooth at $0$; a second-order gap has no first derivative at $0$. Independently: both $\mathrm{loc}$ vanish on
$N_x = \mathrm{span}\{(q, 0)\}$ (the scale invariance of $\mathrm{Log}$ in $q_Y$), so by AJ.2(b) each $M$ is the unique left inverse of its $P$ vanishing there, and equal $P$ give equal $M$. (d), (e) as (a): $\partial_\varphi(\mathrm{Exp}\,\varphi\,t) = -[t]_\times$
and $\partial_\varphi(q\otimes\mathrm{Exp}\,\varphi)$ for the left chart is $\frac12Q_R(q)E = PR^\top$ (AJ.7). $\square$

**Checked:** mpmath 1.4.1, 100 digits, 13 cases: 10 random $(q, t)$ ($t \sim \mathcal N(0, 2^2)$), one with $t \sim \mathcal N(0, (10^6)^2)$, $q = 1$, and $q = (0, 0, 0, 1)$ ($\theta = \pi$). Charts coded from their definitions: the translation of `Screw` from `mp.expm` of the
$4\times4$ hat matrix, its rotation consistent with $R\,\mathrm{Exp}\,\varphi$ ($5.7\times10^{-101}$), $\mathrm{loc}$ by `mp.logm` of $X^{-1}Y$ ($4\times4$ for `Screw`, $3\times3$ for the rotation of the others), $Y$ with a non-unit $q_Y$ (projected). $P_{SE3}$ ($7\times6$)
and $M_{SE3}$ ($6\times7$) by central differences ($h = 10^{-32}$) of the ambient maps against (a): $6.3\times10^{-64}$ and $6.5\times10^{-64}$ (`Screw`, `Decoupled`), $6.7\times10^{-64}$ (`WorldTranslation`, against $\mathrm{diag}(\cdot, I)$);
`Screw` minus `Decoupled`: $P$ exactly $0$, $M$ $7.1\times10^{-102}$. (b): $2.9\times10^{-101}$ ($MP - I_6$, $PM - \mathrm{diag}(I - qq^\top, I)$), singular values $\le 7.1\times10^{-101}$ off. (c), $\delta = s\delta_0$, $s = 10^{-12}$, one random $\delta_0$ per case:
the rotations agree to $10^{-90}$; $\lVert\mathrm{gap} - \frac12R(\varphi\times\rho)s^2\rVert/s^3 \le 1.2$; $\lVert\mathrm{gap}\rVert/(\frac14\lVert\delta\rVert^2) \le 0.9905$. (e) left chart, 30 $(q, t)$, $\mathrm{Exp}(\delta)X$ by `mp.expm`: $9.7\times10^{-65}$. Scripts not committed.
**Permanent:** planned, corpus `se3_plus_jacobian` (`PHASE6.md` §1; $P_{SE3}$ only); none for $M_{SE3}$, `WorldTranslation` or the left chart.

## 4. Ceres and GTSAM

**Ceres** ([Ceres]). Read from `QuaternionPlusImpl`, `QuaternionPlusJacobianImpl`, `QuaternionMinusImpl` and `QuaternionMinusJacobianImpl` in `internal/ceres/manifold.cc` and the comment of `QuaternionManifold` in `include/ceres/manifold.h` (Ceres at
commit `e17a9b4`, 2026-09-26). Storage $(w, x, y, z)$, Hamilton product, unit norm; $\mathrm{Plus}(x, \delta_C) = \big(\cos\lVert\delta_C\rVert,\ \tfrac{\sin\lVert\delta_C\rVert}{\lVert\delta_C\rVert}\delta_C\big)\otimes x$;
$\mathrm{Minus}(y, x) = \mathrm{atan2}(\lVert u\rVert, w)\,u/\lVert u\rVert$ for $(w, u) = y\otimes x^*$; `PlusJacobian` $= D_2\mathrm{Plus}(x, 0)$ ($4\times3$), `MinusJacobian` $= D_1\mathrm{Minus}(x, x)$ ($3\times4$), the definitions of AJ.1.

**Proposition AJ.7 (Ceres' quaternion manifold).**

- (a) $\mathrm{Plus}(x, \delta_C) = \mathrm{Exp}(2\delta_C)\otimes x$ for every $\delta_C$: Ceres' chart is `LeftChart<SO3>` in the coordinates $\delta_C = \varphi/2$, a **half angle**. For $y\cdot x \ge 0$, $\mathrm{Minus}(y, x) = \frac12\,\mathrm{Log}(y\otimes x^*) = \frac12(y \ominus_L x)$;
  for $y\cdot x < 0$ it returns the half-angle vector of the *unflipped* quaternion, $\frac12(y \ominus_L x) + \pi\hat u$ with $\hat u = \mathrm{vec}(y\otimes x^*)/\lVert\cdot\rVert$ (the same rotation on the long branch, LG.2(c); `NUMERICS.md` §3.2 flips), so
  $\mathrm{Minus}(-y, x) \ne \mathrm{Minus}(y, x)$. At $y = x$ the two agree. This is the `Minus` of the cited commit, read from the source and not run in a released binary; the tags 2.1.0 (the first release with `QuaternionManifold`) and 2.2.0 have the same
  $\mathrm{atan2}(\lVert u\rVert, w)$ and no flip.
- (b) The formulas of `QuaternionPlusJacobianImpl` and `QuaternionMinusJacobianImpl` are, for every $x = (w, v) \ne 0$, $P_C = Q_R(x)E = \begin{bmatrix}-v^\top\\ wI - [v]_\times\end{bmatrix}$ (the exact $D_2\mathrm{Plus}(x, 0)$) and $M_C = P_C^\top = \big[-v \;\big\vert\; wI + [v]_\times\big]$
  (the exact $D_1\mathrm{Minus}(x, x)$ for unit $x$; $\lVert x\rVert^2$ times it in general, as AJ.5(a)). With $P^L = \frac12Q_R(x)E$, $M^L = 2E^\top Q_R(x)^\top$ the Jacobians of the left chart (`lplus`) and $R = R(\hat x)$,
  $$P_C = 2\,P^L = 2\,P\,R^\top,\qquad M_C = \tfrac12\,M^L = \tfrac12\,R\,M,\qquad M_CP_C = \lVert x\rVert^2I_3,\qquad P_CM_C = PM.$$
  A step converts as $\varphi = 2R^\top\delta_C$, a Jacobian as $J_{\mathrm{hel}} = J_C\cdot\frac12R$, a covariance as $\Sigma_C = \frac14R\,\Sigma_{\mathrm{hel}}R^\top$.
- (c) The Eigen storage $(x, y, z, w)$ of `EigenQuaternionManifold` is the same map with the rows of $P_C$ (columns of $M_C$) permuted: $\Pi P_C$ for the row permutation $\Pi$ of `Quat::from_xyzw`/`to_xyzw` (`API.md` R3).

*Same as `PHASE6.md` §1:* the storage order $(w, x, y, z)$ and the Hamilton convention; unit-norm storage with no renormalization, so $M_CP_C = \lVert x\rVert^2I$ as in AJ.4; the definitions of the two Jacobians; the projector $P_CM_C = PM = I - xx^\top$; the atan2 form of the
inverse, invariant under a positive scale. *Different:* **the side** (left, $\mathrm{Exp}\otimes x$, against right, $x\otimes\mathrm{Exp}$: the sign of $[v]_\times$ in the lower block, $wI - [v]_\times$ against $wI + [v]_\times$, equivalently the factor $R^\top$); **the factor 2**
(Ceres' tangent is the half-angle vector: $P_C$ has no $\frac12$, $M_C$ no $2$, singular values $1$ against $\frac12$); the sign flip of $\mathrm{Minus}$. The gaps are $\lVert P_C - P\rVert_2 = \frac12\sqrt{5 - 4\cos\theta} \in [\frac12, \frac32]$ and $\lVert M_C - M\rVert_2 = \sqrt{5 - 4\cos\theta} \in [1, 3]$; at $R = I$
the whole difference is the factor $2$, and the side is invisible there, as CH.3(c) says of the charts. Ceres attaches a `Manifold` to one parameter block: `Screw` and `Decoupled` couple $t$ to $q$ and need a single $7$-block; two blocks, `QuaternionManifold` and a Euclidean $t$, give $\mathrm{diag}(P_C, I)$,
$(\mathrm{Exp}(2\delta_q)\otimes q,\ t + \delta_t)$, a left rotation with a world translation, none of the three charts of `0012`.

*Proof.* (a) $\mathrm{Exp}(2\delta_C) = (\cos\lVert\delta_C\rVert, \frac{\sin\lVert\delta_C\rVert}{\lVert\delta_C\rVert}\delta_C)$ (`NUMERICS.md` §3.1 for the tangent $2\delta_C$). For $r = y\otimes x^*$ with $r_w < 0$, $\mathrm{atan2}(\lVert u\rVert, r_w) = \pi - \theta'$ with $\theta' = \mathrm{atan2}(\lVert u\rVert, -r_w)$, while the flip
gives $-\theta'\hat u$ (half of the canonical $\mathrm{Log}$). (b) $(1, \delta)\otimes x = Q_R(x)(1, \delta)$ gives $P_C$; $x\otimes(0, R^\top\delta) = (0, \delta)\otimes x$ gives $Q_R(x)E = Q_L(x)E\,R^\top = 2PR^\top$ (for $\hat x$; both sides scale by $\lVert x\rVert$),
and $Q_R(x)E = 2P^L$ by AJ.3(a) with the sides exchanged. $\mathrm{Minus}(y, x)$ is, near $y = x$, $\mathrm{vec}(r)/\mathrm{sc}(r)$ for $r = y\otimes x^* = Q_R(x^*)y$ (the half of AJ.3(b)'s series), whose derivative at $y = x$ is $E^\top Q_R(x^*)/\lVert x\rVert^2$, and
$E^\top Q_R(x^*) = E^\top Q_R(x)^\top = P_C^\top$; the products are AJ.4 with $Q_R$ for $Q_L$. The gaps: $\lVert PA\rVert_2 = \frac12\lVert A\rVert_2$ ($P^\top P = \frac14I$) with $A = 2R^\top - I$, whose singular values are $1$ (the axis) and
$\lvert2e^{\mp i\theta} - 1\rvert = \sqrt{5 - 4\cos\theta}$; $M_C - M = (\frac12R - I)M$ and $MM^\top = 4I$. (c) A permutation. $\square$

**Checked:** Ceres, mpmath 1.4.1, 100 digits, $\mathrm{Plus}$, $\mathrm{Minus}$, and the two Jacobian formulas transcribed from `manifold.cc` (Order $= (w, x, y, z)$), 40 random unit $x$. (a): $\mathrm{Plus}(x, \delta_C) - \mathrm{Exp}(2\delta_C)\otimes x$, $\delta_C \sim \mathcal N(0, 1.2^2)^3$ (finite, up to
rotations of angle $6$): $3.0\times10^{-101}$; $\mathrm{Minus}$ against $\frac12\mathrm{Log}(y\otimes x^*)$ (`mp.logm`, $\lvert y\cdot x\rvert \ge 0.3$ so that it is principal) for $y\cdot x \ge 0.3$: $1.4\times10^{-101}$, and for $y\cdot x \le -0.3$ against $\frac12\mathrm{Log}(y\otimes x^*) + \pi\hat u$
(22 samples): $2.9\times10^{-101}$. (b): the transcribed Jacobians against central differences of the transcribed maps $1.6\times10^{-65}$, $3.3\times10^{-65}$; $P_C - 2PR^\top$ and $M_C - \frac12RM$ $2.1\times10^{-101}$; $P_C - 2P^L$ with $P^L$ by differencing $\mathrm{Exp}(\varphi)\otimes x$
$4.1\times10^{-66}$; $M_CP_C - I$ $1.4\times10^{-101}$; $P_CM_C - PM$ exactly $0$; the two gap norms against the closed forms $5.7\times10^{-101}$, $1.1\times10^{-100}$. Exact, sympy 1.14.0: $\lVert x\rVert^2P_C = 2PR_\square^\top$, $\lVert x\rVert^2M_C = \frac12R_\square M$, $P_CM_C = PM$, $M_CP_C = \lVert x\rVert^2I$, as polynomial identities.
The real library, `pyceres` 2.6 (bundled Ceres version not stated by the wheel), `f64`, 30 random unit $x$: the Jacobian returned by `Problem::Evaluate` for the residual $r(x) = x$ (identity ambient Jacobian, so it is $P_C$) equals the formula of `manifold.cc` exactly (difference $0$), for `QuaternionManifold` and, with rows permuted, for
`EigenQuaternionManifold` ($0$); against $2PR^\top$ $2.2\times10^{-16}$; one Levenberg–Marquardt step (`max_num_iterations` $= 1$, trust-region radius $10^{16}$) of $r(x) = x - c$, $c$ at $0.3$ to $1.0$ rad, reaches $\mathrm{Plus}(x_0, -P_C^\top(x_0 - c))$ as read from the source to $2.2\times10^{-16}$, and lies $\ge 0.075$ from
the same step applied by the right chart with the full angle. `Minus` and `MinusJacobian` are not exposed by `pyceres`: they are read from the source and checked as above, and no released binary was run on them; `QuaternionMinusImpl` at the tags 2.1.0 and 2.2.0 differs from the cited commit's only in how $\lVert u\rVert$ and $y\otimes x^*$ are evaluated (`sqrt` or `hypot`, an inline product), not in the map. Scripts not committed. **Permanent:** none specified; the oracle comparison of `PHASE6.md` §7 (Ceres' `QuaternionManifold` through a container runner) needs the conversion of (b), or the left-chart twin (AJ.9(a)).

**GTSAM** ([GTSAM]). Read from `LieGroup::retract`, `Rot3::ChartAtOrigin`, `Pose3::ChartAtOrigin` and `Quaternion`'s `Logmap` (`gtsam/base/Lie.h`, `gtsam/geometry/Rot3.h`, `Pose3.cpp`, `Quaternion.h`; develop at `c786d78`, 2026-09-29) and run through the 4.3.0 wheel.

**Proposition AJ.8 (GTSAM).**

- (a) `retract(xi)` is `compose(ChartAtOrigin::Retract(xi))`: a **right** perturbation with the full rotation vector (no half angle). In the default build (`GTSAM_ROT3_EXPMAP` and `GTSAM_POSE3_EXPMAP` on) `Rot3` (`ROT3_DEFAULT_COORDINATES_MODE` $= \mathrm{EXPMAP}$) is `RightChart<SO3>` and `Pose3`
  retracts by $T\,\mathrm{Exp}(\xi)$, $\xi = [\omega; \rho]$ rotation first: `Screw`, in this repository's tangent order, so `gtsam_pose3_to_helicoid` (`PHASE1.md` §7) has no tangent order to convert (unlike `sophus_to_helicoid_tangent`); what is left is the value, `Pose3` as $(R, t)$
  and `Rot3` as a matrix or a quaternion. CMake couples the two options (either on forces the other on), so `Pose3`'s is off only with `Rot3`'s, and `Rot3` then retracts by the Cayley map: $(R\,\mathrm{Cayley}(\omega),\ t + R\rho)$, which is not `Decoupled` (that needs $\mathrm{Exp}$); it is `Decoupled` only under
  `GTSAM_USE_QUATERNIONS`, which keeps `Rot3` on $\mathrm{Exp}$.
- (b) There is no ambient Jacobian. The variable is stored as it is (`Rot3`: a $3\times3$ matrix, a quaternion under `GTSAM_USE_QUATERNIONS`; `Pose3`: $(R, t)$) and `retract` is applied to the storage; a factor returns $\partial/\partial\delta$ of its error at $\delta = 0$, which is
  $J_{\mathrm{amb}}P$ already. The Jacobians of `retract(xi, H1, H2)` are tangent to tangent: $H_1 = \mathrm{Ad}_{\mathrm{Exp}(\xi)}^{-1}$ and $H_2$ the derivative of the chart at the origin, `ExpmapDerivative`, the right Jacobian $J_r(\xi)$ (`RightChart::retract_jacobian`, CH.3; for
  `Pose3` in the rotation-first order). A consumer with this structure needs `Chart`, not `AmbientChart`.
- (c) `Quaternion`'s `Logmap` is invariant under a positive rescaling of $q$ (its comment says so) and negates $q$ when $\lnot(w > 0)$: it differs from `NUMERICS.md` §3.2 only at $w = 0$ (index, open items), so AJ.3(b) holds for it.

**Checked:** GTSAM, `f64`, 4.3.0 wheel (default build), 30 random $(R, \omega)$ and $(T, \xi)$: `Rot3.retract(omega)` against $R\exp(\hat\omega)$ $3.3\times10^{-16}$ (against $\exp(\hat\omega)R$: $\ge 0.11$);
`Pose3.retract(xi)` against $T\exp(\hat\xi)$ ($\xi = [\omega; \rho]$) $2.8\times10^{-16}$ (its translation against $t + R\rho$: $\ge 0.014$; with the halves of $\xi$ swapped: $\ge 0.063$); `localCoordinates(retract(xi))` $- \xi$ $2.2\times10^{-16}$; `Rot3.ExpmapDerivative` and `Pose3.ExpmapDerivative` against $J_r$ by differencing $\mathrm{Log}(\mathrm{Exp}(\xi)^{-1}\mathrm{Exp}(\xi + d))$ ($h = 10^{-7}$, 10 samples): $1.6\times10^{-8}$, $1.4\times10^{-8}$.
The options-off modes of (a) and the coupling of the two options are read from `Rot3.h`, `Pose3.cpp`, `Rot3M.cpp` and `cmake/HandleGeneralOptions.cmake` at `c786d78`, not run. Scripts not committed. **Permanent:** planned, the `gtsam` runner of `PHASE1.md` §7 (`Expmap`/`Logmap` and their derivatives) for the `ExpmapDerivative` of (b); none for (a) and (c).

## 5. Conditioning and rounding

**Proposition AJ.9.**

- (a) *Exactness.* The entries of $P$ and $M$ are $\pm\frac12q_i$ and $\pm2q_i$: a scaling by a power of two and a negation, exact in binary floating point except for the underflow of $\frac12q_i$ below the normal range. The stored $P$ and $M$ are those of the stored $q$, with no rounding error of their own; the same for $P_{SE3}$ and $M_{SE3}$ except the $R$ block, whose entries carry the rounding of $R(q)$. The entries of Ceres' $P_C$, $M_C$ are $\pm x_i$, also exact, so $2P^L = P_C$ bit for bit, whereas $P_C = 2PR^\top$ rounds in the product with $R$.
- (b) *Conditioning.* For unit $q$, $\kappa_2 = 1$ on the range for $P$ and for $M$, whatever the unit $q$ (AJ.4): $\lVert P\rVert_2 = \frac12$, $\lVert M\rVert_2 = 2$, $\lVert PM\rVert_2 = 1$. There is no singular $q$: these are the Jacobians of the chart at $\delta = 0$, not $\mathrm{rj}$ or $\mathrm{lj}$, whose conditioning grows with $\lVert\varphi\rVert$
  (LG.16; $J^{-1}$ is singular at $2\pi$, `NUMERICS.md` §12). For SE(3), $\kappa_2(P_{SE3}) = \kappa_2(M_{SE3}) = 2$ (AJ.6(b)): the half angle of the quaternion against metres.
- (c) *The product.* $J_{\mathrm{tan}} = J_{\mathrm{amb}}P$ is a length-$4$ inner product per entry, one factor exact but every product $J_{\mathrm{amb}}\cdot\frac12q_i$ rounding: $\lvert\mathrm{fl}(J_{\mathrm{amb}}P) - J_{\mathrm{amb}}P\rvert \le \gamma_4\lvert J_{\mathrm{amb}}\rvert\lvert P\rvert$ ([Higham], ch. 3), $\gamma_4 = 4u/(1 - 4u)$, and it halves an error
  in $J_{\mathrm{amb}}$ ($\lVert P\rVert_2 = \frac12$). It inherits the ambient Jacobian's own error (for a `Dual` Jacobian, the derivative of the shipped code). The bound is componentwise in $\lvert J_{\mathrm{amb}}\rvert$, so it depends on the extension: a radial part $aq^\top$ of $\mathrm Df$ (AJ.2(c), the example) cancels in exact
  arithmetic ($q^\top P = 0$) but not in floating point, and adds $\gamma_4\lvert a\rvert\lvert q\rvert^\top\lvert P\rvert$ to the error of $J_{\mathrm{tan}}$, whatever its size: an extension whose radial derivative is small keeps the error at the size of $J_{\mathrm{tan}}$.
- (d) *Drift.* For a non-unit stored $q$ the pair is exact but $MP = \lVert q\rVert^2I$ (AJ.5): in the model of AJ.5(b), $\lvert\eta\rvert$ is $\approx1.1\sqrt n\,u$ (rms; a fit over $n \le 10^5$) after $n$ independent retractions and grows linearly, $\approx0.7u$ per retraction, for a repeated step; a consumer that needs $MP = I$ to rounding renormalizes (`renormalize`, `NUMERICS.md` §3.6) at a fixed cadence
  set by the linear case (about $10^4$ retractions from the bound of `from_wxyz_unchecked` in `f64`, a few hundred in `f32`), not by the random walk.

*Proof.* (a), (b) AJ.3, AJ.4, AJ.6. (c) The inner-product bound of [Higham], ch. 3, and $\lVert P\rVert_2 = \frac12$. (d) AJ.5. $\square$

**Checked:** (a)–(b) are AJ.3, AJ.4, AJ.6 (singular values to $10^{-100}$); the bit-exactness of $P_C$ is the $0$ difference of AJ.7's `pyceres` run. (c): no run (the bound, and the radial term it gives, are [Higham]'s). (d) AJ.5. **Permanent:** none specified.

## 6. Sign audit

**Table AJ.10.** The correct entry, the usual wrong variant and their **exact gap** (unit $q$; $\theta$ the rotation angle of $q$). Every gap vanishes at $q = (\pm1, 0)$ except the factors and signs; a test needs $q$ away from the identity, and for the translation block of $M$ (whose gap $2\sin\theta$ also vanishes at the half turn, where $R = R^\top$) $\theta \in (0, \pi)$.

| Where | Correct | Wrong variant | Gap | Reason |
|---|---|---|---|---|
| side of $P$ | $\frac12\begin{bmatrix}-v^\top\\ wI + [v]_\times\end{bmatrix}$ (right) | the lower block $wI - [v]_\times$ (Ceres' left) | $\lVert P - P^L\rVert_2 = \sin\frac\theta2$ | AJ.3, AJ.7; $P^L = PR^\top$ |
| $P$ and $M$ of different sides | $MP = I$ | $MP^L = R^\top$, $M^LP = R$ | $\lVert MP^L - I\rVert_2 = 2\sin\frac\theta2$ | AJ.7(b): caught by $MP = I$ |
| scale of the tangent | $P$, $M$ for the rotation vector | Ceres' $P_C = 2PR^\top$, $M_C = \frac12RM$ | $\lVert P_C - P\rVert_2 = \frac12\sqrt{5 - 4\cos\theta}$, $\lVert M_C - M\rVert_2 = \sqrt{5 - 4\cos\theta}$ | AJ.7; at $R = I$ the factor $2$ alone |
| $M$ from $P$ | $M = 4P^\top$ | $M = P^\top$ (Ceres' $M_C = P_C^\top$) | $MP = \frac14I$: gap $\frac34$ | the rotation vector, not the half angle |
| normal part of $M$ | $Mq = 0$ | $M + zq^\top$ | $\lVert zq^\top\rVert$, arbitrary; **$MP = I$ still holds** | AJ.2(b): a test needs $Mq = 0$ or $PM$ |
| $M$ at a non-unit $q$ | $M/\lVert q\rVert^2 = \mathrm D\,\mathrm{loc}$, the derivative of the definition (AJ.1; B) | the printed $M$ (A); $M/\lVert q\rVert$ (C): differ from B, and the reading a reference takes is an open item | $\lvert\eta\rvert$; $\lvert\eta\rvert/2$ | AJ.5(d) |
| representative | $P$, $M$ at the same $q$ | $M(-q)P(q)$ | $-I$: gap $2$ | AJ.5(e) |
| SE(3) translation block of $P$ | $R$ | $I$ (`WorldTranslation`) | $\lVert R - I\rVert_2 = 2\sin\frac\theta2$ | AJ.6(e) |
| SE(3) translation block of $M$ | $R^\top$ | $R$ | the block of $MP$ is $R^2$: $\lVert R^2 - I\rVert_2 = 2\sin\theta$ | AJ.6(a) |
| the retraction with the same $P$ | $q\otimes\mathrm{Exp}\,\delta$ | $(q + P\delta)/\lVert q + P\delta\rVert$ | $\frac1{12}\lVert\delta\rVert^2P\delta + O(\lVert\delta\rVert^4)$ | AJ.3(d): both equal to second order; $P$, $M$ cannot tell them apart |

**Checked:** mpmath 1.4.1, 100 digits, 30 unit $q$: $\lVert P - P^L\rVert_2 - \sin\frac\theta2$ $2.9\times10^{-101}$; $MP^L - R^\top$, $M^LP - R$ $1.4\times10^{-101}$; $\lVert MP^L - I\rVert_2 - 2\sin\frac\theta2$ $5.7\times10^{-101}$; $P_C - P$, $M_C - M$ 2-norms (AJ.7); $P^\top P - \frac14I$ $3.6\times10^{-102}$
($M := P^\top$ gives $MP = P^\top P = \frac14I$); $\lVert R^2 - I\rVert_2 - 2\lvert\sin\theta\rvert$ $1.0\times10^{-100}$ and $\lVert R - I\rVert_2 - 2\sin\frac\theta2$ $7.1\times10^{-101}$; the last row, $\delta = s\delta_0$, $s = 10^{-12}$, 30 samples: $\mathrm{ret} - (q + P\delta)/\lVert q + P\delta\rVert$ equals
$\frac1{12}\lVert\delta\rVert^2P\delta - \frac1{48}\lVert\delta\rVert^4q$ to $11s^5$. Script not committed. **Permanent:** as the propositions cited.
