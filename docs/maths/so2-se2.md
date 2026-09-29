# SO(2) and SE(2): Exp, Log, adjoints and the rotation-first Jacobians

> Non-normative companion to [`NUMERICS.md`](../NUMERICS.md) §6, with the SE(2) tangent of §1 and the rows
> of §2.3 it instantiates. **`NUMERICS.md` wins on any conflict; a disagreement is an open item in the
> [maths index](./index.md)**, which holds the notation and the `Checked:` convention. What §6 leaves to
> "the Phase 3 PR" (the SE(2) Jacobians, rotation-first) is **derived** here and marked **Proposed**
> (§6 below); nothing on this page edits a formula of `NUMERICS.md`. The general algebra is
> [`lie-groups.md`](./lie-groups.md), the algebra reused for $J$ and $J^{-1}$ is that of
> [`so3.md`](./so3.md) SO.7–SO.8, the coefficients are [`coefficients.md`](./coefficients.md)'s.

On this page $\theta\in\mathbb R$ is a **signed** scalar (the SE(2) rotation angle, not a norm; $\theta(X)\in[-\pi,\pi]$ is signed
too, and a domain of `NUMERICS.md` written "$\theta<2\pi$", "$\theta(X)<\pi$", canonical branch $[0,\pi]$, reads $\lvert\theta\rvert$ here), $z = \theta^2$,
$\rho = (\rho_x,\rho_y)^\top$, $E = \mathrm{Exp}(\tau)$, and
$K = \begin{bmatrix}0&-1\\1&0\end{bmatrix}$ is the quarter turn ($K^2 = -I$, $K^\top = -K$; on
$\mathbb R^2\cong\mathbb C$ it is multiplication by $i$). $R(\theta) = \cos\theta\,I + \sin\theta\,K = e^{\theta K}$
is the counter-clockwise (active) rotation; the sign of $\beta$ in $V$ of `NUMERICS.md` §6 forces it (PL.4).
$s(\theta) = 2\sin(\theta/2)/\theta$ is LG.16's ($s(0) = 1$); $\sigma_m$ is CO.1's; $u$ is the unit roundoff;
$\tau_{\mathrm{tf}} = [\rho;\theta] = \Pi\tau$ with $\Pi = \begin{bmatrix}0 & I_2\\ 1 & 0\end{bmatrix}$ (the index's $P$ is the
$6\times6$ case).

## Results

| Label | Result | Status | `NUMERICS.md` | Implemented by |
|---|---|---|---|---|
| PL.1–PL.2 | $\alpha,\beta,a,b,c,\gamma$: definitions, parity, series, singularities; $\beta = \theta a$, $\alpha = 1 - zb$, $\gamma = 1 - zc$ | stated ($\alpha,\beta$: §6; $a,b,c$: §4); $\gamma$ derived | §4, §6 | `coeffs::{se2_coeffs, jr_coeffs, jr_inv_coeff}` |
| PL.3 | SO(2): unit complex, $\mathrm{Exp}$, $\mathrm{Log} = \mathrm{atan2}$, $\mathrm{Ad} = J = 1$, $\mathrm{ad} = 0$; the $\mathrm{atan2}$ sensitivity | stated | §6 | `SO2::{exp, log, adjoint, ad, jr, jl, jr_inv, jl_inv}` |
| PL.4 | $\mathrm{Exp}\,\tau = (R(\theta), V\rho)$, $V = \alpha I + \beta K = s(\theta)R(\theta/2)$, $\det V = s^2$ | stated | §6 | `SE2::exp` |
| PL.5 | $\mathrm{Log}$; three forms of $V^{-1}$; the two preimages at $\lvert\theta\rvert = \pi$ | stated | §6 | `SE2::log` |
| PL.6 | $\mathrm{Ad}_X$, $\mathrm{ad}_\tau$; $\mathrm{ad}^3 = -\theta^2\mathrm{ad}$; $\mathrm{Ad}_E = I + \alpha\,\mathrm{ad} + a\,\mathrm{ad}^2$ | **Proposed** | §6 (SE(2) not stated) | `SE2::{adjoint, ad}` |
| PL.7 | $J_{l,r} = I \pm a\,\mathrm{ad} + b\,\mathrm{ad}^2$ and its blocks | **Proposed** | §6 | `SE2::{jl, jr}` |
| PL.8 | $J_{l,r}^{-1} = I \mp \tfrac12\mathrm{ad} + c\,\mathrm{ad}^2$ and its blocks | **Proposed** | §6, §12, §14 | `SE2::{jl_inv, jr_inv}` |
| PL.9 | SE(2) inside SE(3): §5.3's $Q$ reduces to $a\Xi + b(W\Xi + \Xi W)$; agreement with Solà et al. App. C | derived | §5.3, §6 | — |
| PL.10 | action Jacobians of SO(2), SE(2) | **Proposed** | §2.4 | `SO2::act_jacobians`, `SE2::act_jacobians` |
| PL.11 | numerical consequences: exact arms, `Log`'s $V^{-1}$, $\kappa_2(J)$, domains | derived | §4, §11, §12 | `coeffs`, `SE2::{log, jr_inv}` |

## 1. The scalar functions

**Definition PL.1.** For $\theta\in\mathbb R$:

| | definition | parity | series in $z$ | singularities |
|---|---|---|---|---|
| $\alpha$ | $\sin\theta/\theta$ | even | $\sigma_1 = 1 - \frac z6 + \frac{z^2}{120} - \frac{z^3}{5040}$ | none (entire) |
| $\beta$ | $\frac{1-\cos\theta}\theta = \frac{2\sin^2(\theta/2)}\theta$ | odd | $\theta\,\sigma_2 = \theta\big(\frac12 - \frac z{24} + \frac{z^2}{720} - \frac{z^3}{40320}\big)$ | none |
| $a$ | $(1-\cos\theta)/\theta^2$ | even | $\sigma_2$, the series of $\beta/\theta$ | none |
| $b$ | $(\theta-\sin\theta)/\theta^3$ | even | $\sigma_3 = \frac16 - \frac z{120} + \frac{z^2}{5040} - \frac{z^3}{362880}$ | none |
| $\theta b$ | $(\theta-\sin\theta)/\theta^2 = (1-\alpha)/\theta$ | odd | $\theta\,\sigma_3$ | none |
| $c$ | $\theta^{-2} - \cot(\theta/2)/(2\theta)$ | even | $\frac1{12} + \frac z{720} + \frac{z^2}{30240} + \frac{z^3}{1209600}$ | simple poles at $\theta\in2\pi\mathbb Z_{\ne0}$ (radius $2\pi$) |
| $\gamma$ | $\frac\theta2\cot\frac\theta2$ | even | $1 - \frac z{12} - \frac{z^2}{720} - \frac{z^3}{30240} - \frac{z^4}{1209600}$ | the same poles, $\gamma\approx2\pi/(\theta-2\pi)$ there; $\gamma(\pm\pi) = 0$ |
| $\theta c$ | $(1-\gamma)/\theta$ | odd | $\theta\big(\frac1{12} + \frac z{720} + \frac{z^2}{30240} + \frac{z^3}{1209600}\big)$ | the same poles |

$a$, $b$, $c$ are `NUMERICS.md` §4's; $\alpha$, $\beta$ are §6's. $\theta b$ and $\theta c$ are the scalars that
multiply $\rho$ in the Jacobians (PL.7, PL.8); they are products of an even coefficient and the signed
$\theta$, not new functions.

**Proposition PL.2.** (a) Every entire entry has the series $\sum_j(-z)^j/(2j+m)!$ of CO.1 ($m = 1, 2, 3$
for $\alpha, a, b$), and $c = \sum_j\lvert B_{2j+2}\rvert z^j/(2j+2)!$, $\gamma = 1 - zc$. (b)
$\beta = \theta a = 2\theta k^2$ ($k$ of §4); $\alpha = 1 - zb$; $\gamma = 1 - zc$; $\alpha^2 + \beta^2 = s^2 =
2a$; and $\alpha + i\beta = s\,e^{i\theta/2}$.

*Proof.* (a) CO.2 and CO.3(a); $zc = 1 - \frac\theta2\cot\frac\theta2$ from the definition of $c$. (b) The first
three are CO.4(a),(c) and the definition of $c$. $\alpha + i\beta = \big(\sin\theta + i(1-\cos\theta)\big)/\theta =
\frac{2\sin(\theta/2)}\theta\big(\cos\frac\theta2 + i\sin\frac\theta2\big)$; take the modulus, and
$s^2 = 2(1-\cos\theta)/\theta^2 = 2a$. $\square$

**Checked:** sympy 1.14.0, exact rationals: the Taylor series of $\alpha$, $\beta$, $a$, $\theta b$, $c$, $\theta c$, $\gamma$ equal the general
terms through $\theta^{13}$ (Bernoulli numbers for $c$, $\theta c$, $\gamma$; the leading four agree with the
series column of `NUMERICS.md` §4 for $a$, $b$, $c$); the identities of (b) hold exactly in
$\mathbb Q(\theta,\tan\frac\theta2)$. Script not committed. **Permanent:** planned, the generator asserts
the leading four of $a$, $b$, $c$ (0004 item 2); `PHASE1.md` §6 sweeps $\alpha$, $\beta$ but §4.3 lists no
corpus id for them (only `coeff_k`…`coeff_r`, and `so2_*`, `se2_*` as "analogous").

## 2. SO(2)

**Proposition PL.3.** $\mathrm{SO}(2) = \{R(\theta)\}$ is identified with the unit complex numbers
$(c, s) = c + is$, $R(\theta)v \leftrightarrow e^{i\theta}v$.

- (a) $\mathrm{Exp}\,\theta = (\cos\theta, \sin\theta)$; composition is complex multiplication,
  $(c_1,s_1)(c_2,s_2) = (c_1c_2 - s_1s_2,\ c_1s_2 + s_1c_2)$; the inverse is the conjugate $(c,-s)$.
  $\mathrm{Exp}(\theta_1 + \theta_2) = \mathrm{Exp}\,\theta_1\,\mathrm{Exp}\,\theta_2$ in exact arithmetic (in `f64`
  $\cos$, $\sin$ of $\theta_1 + \theta_2$ and the complex product round differently).
- (b) $\mathrm{Ad} = 1$, $\mathrm{ad} = 0$, $J_r = J_l = J_r^{-1} = J_l^{-1} = 1$: every row of
  `NUMERICS.md` §2.3 has entries $\pm1$ ($Y\ominus X$ is the $\mathrm{atan2}$ of $\bar x\,y$).
- (c) $\mathrm{Log}(c,s) = \mathrm{atan2}(s, c)\in[-\pi,\pi]$ inverts $\mathrm{Exp}$ on $(-\pi,\pi)$ and is
  scale-invariant in exact arithmetic (in `f64` the scaled pair $(kc, ks)$ is exact only for $k$ a power of two). At $c<0$, $s = \pm0$ it returns $\pm\pi$ by the sign bit of the zero: there, as in
  `NUMERICS.md` §3.2, `Log` is a function of the representation, not of the rotation, and it jumps by $2\pi$.
- (d) For unit $(c, s)$ with relative perturbations $\xi_c$, $\xi_s$:
  $\delta\theta = \tfrac12\sin2\theta\,(\xi_s - \xi_c) + O(\xi^2)$, so $\lvert\delta\theta\rvert\le\lvert\theta\rvert\,\lvert\xi_s - \xi_c\rvert$:
  no amplification anywhere, and $\delta\theta = 0$ at $\pm\frac\pi2$, $\pm\pi$. $\arccos c$ instead loses the sign
  of $\theta$ and, by SO.6(a)–(b), $\approx u/(2\theta^2)$ relative near $0$ and $\approx u/(2\pi(\pi - \lvert\theta\rvert))$ near $\pm\pi$.
- (e) The Newton step $w\leftarrow w(3 - \lvert w\rvert^2)/2$ on $w = c + is$ gives, with
  $\eta = \lvert w\rvert^2 - 1$, $\lvert w'\rvert^2 - 1 = -\tfrac34\eta^2 + \tfrac14\eta^3$ exactly (SO.14's computation is
  scalar, so it is dimension-independent). `NUMERICS.md` specifies no SO(2) normalization; this is not one.

*Proof.* (a) $e^{\theta K} = \sum_j(-1)^j\theta^{2j}/(2j)!\,I + \sum_j(-1)^j\theta^{2j+1}/(2j+1)!\,K$ by $K^2 = -I$;
$e^{i\theta}$ is the same series. (b) The group is abelian, so $\mathrm{Ad}_X\tau = X\tau^\wedge X^{-1} = \tau$,
$[\tau^\wedge,\sigma^\wedge] = 0$, and $J = \sum(\mp\mathrm{ad})^n/(n+1)! = 1$. (c) $\mathrm{atan2}$ is homogeneous of degree
$0$ and ranges over $[-\pi,\pi]$. (d) $d\theta = (c\,ds - s\,dc)/(c^2 + s^2) = sc\,(\xi_s - \xi_c)$ with $ds = s\xi_s$,
$dc = c\xi_c$, and $\lvert\sin2\theta\rvert/2\le\lvert\theta\rvert$. (e) $\lvert w'\rvert^2 = (1+\eta)(2-\eta)^2/4$. $\square$

**Checked:** mpmath 1.3.0, 80 digits: (d) $\mathrm{atan2}$ of $(c(1+\xi_c), s(1+\xi_s))$ at $\xi = 10^{-30}$, 11 angles
($10^{-10}$ to $\pm3.1$, $\pm\frac\pi2$) $\times$ 20 random $(\xi_c,\xi_s)$: second-order residual $\le0.49\,\xi^2$,
$\max\lvert\delta\theta\rvert/(\lvert\theta\rvert\lvert\xi_s - \xi_c\rvert) = 1.0$; (e) the identity, $\le2\times10^{-81}$
($\eta = \pm10^{-3}, 0.2, -0.5$). (a)–(c) are definitions. Script not committed. **Permanent:** planned, corpus
`so2_*` (`PHASE1.md` §4.3, "analogous"); proptests `group_axioms_*`, `exp_log_roundtrip_*`.

## 3. SE(2): the group, Exp and Log

$G = \mathrm{SE}(2)$ is the set of $\begin{bmatrix}R & t\\ 0 & 1\end{bmatrix}$, $R\in\mathrm{SO}(2)$, $t\in\mathbb R^2$, with
$(R,t)(S,y) = (RS,\ Ry + t)$ and $(R,t)^{-1} = (R^\top, -R^\top t)$ (SE.2 with $\mathrm{SO}(2)$; `a * b` is
$T_{ax}T_{xb}$). Tangent $\tau = [\theta;\rho_x;\rho_y]$, $\tau^\wedge = \begin{bmatrix}\theta K & \rho\\ 0 & 0\end{bmatrix}$.
$\theta(X) = \mathrm{atan2}(R_{10}, R_{00})\in[-\pi,\pi]$.

**Proposition PL.4 (Exp).** For every $\tau$,
$\mathrm{Exp}\,\tau = (R(\theta),\ V\rho)$ with
$$
V = \sum_{k\ge0}\frac{(\theta K)^k}{(k+1)!} = \alpha I + \beta K = \begin{bmatrix}\alpha & -\beta\\ \beta & \alpha\end{bmatrix} = s(\theta)\,R(\theta/2),
\qquad \det V = s(\theta)^2 .
$$
$V$ is invertible iff $\theta\notin2\pi\mathbb Z_{\ne0}$. $\mathrm{Exp}$ is an analytic bijection of $\{\lvert\theta\rvert<\pi\}\times\mathbb R^2$
onto $\{\lvert\theta(X)\rvert<\pi\}$. It does not wrap: for $m\in\mathbb Z_{\ne0}$ and $\rho\ne0$ the tangents $[\theta;\rho]$ and $[\theta - 2\pi m;\rho]$ have the
same rotation, but $V(\theta)\rho\ne V(\theta - 2\pi m)\rho$ unless $\theta$ and $\theta - 2\pi m$ both lie in
$2\pi\mathbb Z_{\ne0}$ (there $V = 0$ for both).

*Proof.* $A = \tau^\wedge$ has $A^k = \begin{bmatrix}\theta^kK^k & \theta^{k-1}K^{k-1}\rho\\ 0 & 0\end{bmatrix}$ for $k\ge1$ (SE.3(a)), so
$\exp A = (e^{\theta K}, \sum_{k\ge0}(\theta K)^k\rho/(k+1)!)$. $K^2 = -I$ splits the sum by parity:
$V = \sigma_1(z)\,I + \theta\sigma_2(z)\,K$, i.e. PL.1. $s\,R(\theta/2)$ has $s\cos\frac\theta2 = \alpha$ and
$s\sin\frac\theta2 = \beta$ (PL.2(b)); a similarity of ratio $s$ has determinant $s^2$, zero iff
$\sin\frac\theta2 = 0$, $\theta\ne0$. Bijection: $R(\theta)$ is injective on $(-\pi,\pi)$ and $V$ invertible
there. No wrap: with $\theta' = \theta - 2\pi m$ and both $\theta,\theta'\ne0$, $R(\theta'/2) = (-1)^mR(\theta/2)$ and
$\sin\frac{\theta'}2 = (-1)^m\sin\frac\theta2$ give $V(\theta')= \frac{2\sin(\theta/2)}{\theta'}R(\tfrac\theta2)$, so
$V(\theta)\rho - V(\theta')\rho = 2\sin\tfrac\theta2\big(\tfrac1\theta - \tfrac1{\theta'}\big)R(\tfrac\theta2)\rho$, zero iff $\sin\frac\theta2 = 0$; if
one of $\theta,\theta'$ is $0$ the other is in $2\pi\mathbb Z_{\ne0}$, and $V(0)\rho = \rho\ne0$ against $V = 0$. $\square$

*Reading.* A constant twist traces a circular arc; the chord $V\rho$ has length $s\lVert\rho\rVert = 2\lvert\sin\frac\theta2\rvert\,\lVert\rho\rVert/\lvert\theta\rvert$
and is the initial displacement $\rho$ turned by $\theta/2$: the sign of $\beta$ in `NUMERICS.md` §6 is
the counter-clockwise convention of $R(\theta)$.

**Checked:** mpmath 1.3.0, 110 digits, 81 tangents: $\theta\in\{0\}\cup\pm\{10^{-30}, 10^{-10}, 10^{-3}, 0.5, 1, 2, 3, 5, 6.5, 20, \pi\pm10^{-3},
\pi\pm10^{-6}, \pi, 2\pi - 10^{-3}\}$, 40 uniform on $[-6,6]$ and 8 with $\rho\sim N(0,(10^6)^2)$ at $\theta\in\{\pm0.3, \pm2.9,
\pm10^{-10}, \pm3.14\}$, $\rho\sim N(0,2^2)$ otherwise; $\alpha,\beta$ from their definitions with guard digits.
`mp.expm` of the hat matrix against the blocks: $1.6\times10^{-111}$ (relative to $\max(1,\lVert\cdot\rVert_{\max})$);
$V$ by its series against $\alpha I + \beta K$: $1.9\times10^{-105}$ (at $\theta = 20$); against $sR(\theta/2)$:
$1.7\times10^{-111}$; $\det V - s^2$: $1.7\times10^{-111}$. No wrap (80 digits, $\rho = (1.3, -0.7)$): $V(\theta)\rho - V(\theta')\rho$
agrees with the closed form to $1.6\times10^{-81}$ at $(\theta, m) = (7, 2)$, and is $1.3 = \lVert\rho\rVert_{\max}$ at $(2\pi, 1)$, $(0, 1)$ and $(4\pi, 2)$
but $0$ at $(4\pi, 1)$, $(2\pi, 2)$. Script not committed. **Permanent:** planned, corpus
`se2_*` (`mp.expm` of the $3\times3$ hat matrix), `coeff_*`; proptest `exp_log_roundtrip_*`.

**Proposition PL.5 (Log).** For $X = (R, t)$ let $\theta = \theta(X)$. Then $\mathrm{Log}\,X = (\theta,\ V(\theta)^{-1}t)$ with
$$
V(\theta)^{-1} = \gamma I - \tfrac\theta2K = \frac1{s}\,R(-\tfrac\theta2) = \frac1{\alpha^2 + \beta^2}\begin{bmatrix}\alpha & \beta\\ -\beta & \alpha\end{bmatrix},
$$
$\mathrm{Exp}(\mathrm{Log}\,X) = X$ for every $X$, and $\mathrm{Log}(\mathrm{Exp}\,\tau) = \tau$ for $\lvert\theta\rvert<\pi$ (not for
$\lvert\theta\rvert>\pi$). At $\lvert\theta(X)\rvert = \pi$ there are exactly two preimages with $\lvert\theta\rvert<2\pi$, $(\pi, -\frac\pi2Kt)$ and
$(-\pi, +\frac\pi2Kt)$, so $\rho' - \rho = \pi Kt$ (LG.2(d), $W\to\theta K$; without the bound there are infinitely many,
$\theta\in\pm\pi + 2\pi\mathbb Z$, since $V$ is invertible there); which is returned follows the sign
bit of $R_{10}$, so there $\rho$ is a function of the representation, and $\mathrm{Log}$ jumps.

*Proof.* $\mathrm{Exp}(\theta(X), V^{-1}t) = (R, VV^{-1}t) = X$ since $R(\theta(X)) = R$. $V^{-1} = s^{-1}R(-\theta/2)$
from $V = sR(\theta/2)$; $s^{-1}\cos\frac\theta2 = \frac\theta2\cot\frac\theta2 = \gamma$ and $s^{-1}\sin\frac\theta2 = \frac\theta2$; and
$V^\top/(\alpha^2+\beta^2) = sR(-\frac\theta2)/s^2$. Inverse on $\lvert\theta\rvert<\pi$: PL.4. At $\theta = \pi$,
$\gamma = 0$ gives $V(\pm\pi)^{-1} = \mp\frac\pi2K$, and $R(\pi) = R(-\pi)$. $\square$

**Checked:** as PL.4 (110 digits). $V^{-1}$ by `mp.inverse` of the series $V$ against all three forms (81 tangents,
$\lvert\theta\rvert<6.3$, $\bigl\lvert\lvert\theta\rvert - 2\pi\bigr\rvert\ge10^{-2}$): $1.2\times10^{-109}$.
$\mathrm{Log}(\mathrm{Exp}\,\tau) = \tau$ for $\lvert\theta\rvert<\pi - 10^{-7}$: $\theta$ $9.7\times10^{-122}$, $\rho$ $1.5\times10^{-111}$.
$\mathrm{Exp}(\mathrm{Log}\,X) = X$: 60 $X$ with $\theta(X)\sim U[-\pi,\pi]$, $t\sim N(0,2^2)$ plus ten edge $X$
($\theta(X)\in\pm(\pi - 10^{-3}), \pm10^{-10}, \pi - 10^{-12}$; $t\sim N(0,2^2)$ and $N(0,(10^6)^2)$): $1.1\times10^{-111}$.
Preimages at $\pi$ (10 $\rho$): $\mathrm{Exp}(-\pi,\rho') = X$ $2.4\times10^{-111}$, $\rho'-\rho = \pi Kt$ $4.9\times10^{-111}$.
Beyond $2\pi$ (80 digits, $X = (-I, (1.3, -0.7))$): $\mathrm{Exp}(\theta, V(\theta)^{-1}t) = X$ for $\theta\in\{\pm\pi, \pm3\pi, 5\pi\}$, $\le9.6\times10^{-81}$.
`mp.logm` of the $3\times3$ matrix agrees with $(\theta, V^{-1}t)$ for $\lvert\theta\rvert\le3.02$ ($9\times10^{-41}$, 40
digits, 40 $X$) and returns a complex, non-principal result from an onset in $(3.02, 3.03]$ (8 translations $\times$
both signs, grid $0.005$ on $[2.9, 3.1]$; identical in mpmath 1.3.0 and 1.4.1): the fault of `so3_log` and
`sen3_log_n*` (index, open items) applies to any `se2_log` defined by `mp.logm`. Script not committed.
**Permanent:** planned, corpus `se2_*` (`Log`), proptest `exp_log_roundtrip_*`.

## 4. Ad and ad

**Proposition PL.6 (`adjoint`, `ad`; Proposed for `NUMERICS.md` §6).** For $X = (R, t)$, $\tau = [\theta;\rho]$, $\psi = [\psi_\theta;\eta]$ (no domain):

- (a) $\mathrm{Ad}_X\psi = [\psi_\theta;\ R\eta - \psi_\theta Kt]$ and $\mathrm{ad}_\tau\psi = [0;\ \theta K\eta - \psi_\theta K\rho]$, i.e.
$$
\mathrm{Ad}_X = \begin{bmatrix}1 & 0\\ -Kt & R\end{bmatrix},\qquad
\mathrm{ad}_\tau = \begin{bmatrix}0 & 0\\ -K\rho & \theta K\end{bmatrix},
$$
  $X\mathrm{Exp}(\psi)X^{-1} = \mathrm{Exp}(\mathrm{Ad}_X\psi)$, $\mathrm{Ad}_{XY} = \mathrm{Ad}_X\mathrm{Ad}_Y$,
  $\mathrm{Ad}_X^{-1} = \mathrm{Ad}_{X^{-1}}$ ($-K(-R^\top t) = R^\top Kt$ in the coupling block).
- (b) $\mathrm{ad}_\tau^2 = \begin{bmatrix}0&0\\ \theta\rho & -\theta^2I\end{bmatrix}$ and $\mathrm{ad}_\tau^3 = -\theta^2\,\mathrm{ad}_\tau$.
- (c) $\mathrm{Ad}_E = \exp(\mathrm{ad}_\tau) = I + \alpha\,\mathrm{ad}_\tau + a\,\mathrm{ad}_\tau^2$.

*Proof.* (a) $X\psi^\wedge X^{-1}$ has top-left $R(\psi_\theta K)R^\top = \psi_\theta K$ (planar rotations commute: the
rotation entry of $\mathrm{Ad}$ is $1$, not $R$) and top-right $R\eta - R(\psi_\theta K)R^\top t = R\eta - \psi_\theta Kt$; the commutator $[\tau^\wedge,\psi^\wedge]$ has top-left $0$ ($K$ commutes with itself) and
top-right $\theta K\eta - \psi_\theta K\rho$. The identities are LG.4(a),(b), which use only matrix products. (b) With
$M = \theta K$, $u = -K\rho$: $\mathrm{ad}^2 = \begin{bmatrix}0&0\\ Mu & M^2\end{bmatrix}$, $Mu = -\theta K^2\rho = \theta\rho$,
$M^2 = -\theta^2I$, and $\mathrm{ad}^3 = \begin{bmatrix}0&0\\ \theta M\rho & -\theta^2M\end{bmatrix} = -\theta^2\mathrm{ad}$ since
$\theta M\rho = \theta^2K\rho = -\theta^2u$. (c) LG.5, then (b): $\mathrm{ad}^{2j+1} = (-z)^j\mathrm{ad}$,
$\mathrm{ad}^{2j+2} = (-z)^j\mathrm{ad}^2$, and $\sum(-z)^j/(2j+1)! = \alpha$, $\sum(-z)^j/(2j+2)! = a$. $\square$

$-Kt = (t_y, -t_x) = t\times e_z$: rotating about the origin moves a body at $t$. The minimal polynomial
$\lambda^3 + \theta^2\lambda$ of $\mathrm{ad}_\tau$ is that of $W = [\varphi]_\times$ (SO.3): everything below is
SO.7–SO.8 with $W\to\mathrm{ad}_\tau$.

**Checked:** mpmath 1.3.0, 110 digits. $\mathrm{Ad}_X$ from Definition LG.3 (images of the basis under
$\psi\mapsto(X\psi^\wedge X^{-1})^\vee$, `mp.inverse`), 25 random $X$ ($\theta(X)\in[-3.1,3.1]$, $t\sim N(0,2^2)$)
plus 8 edge ($\theta(X) = \pm10^{-10}, \pm(\pi - 10^{-3})$; $t\sim N(0,2^2)$, $N(0,(10^6)^2)$): block form exact (bitwise:
mpmath's products and inverse round correctly at this precision); images stay in $\mathfrak g$ ($5\times10^{-115}$);
$\mathrm{Ad}_{X^{-1}} = \mathrm{Ad}_X^{-1}$ $1.2\times10^{-111}$; $\mathrm{Ad}_{XY} = \mathrm{Ad}_X\mathrm{Ad}_Y$ exact;
$X\mathrm{Exp}(\psi)X^{-1} = \mathrm{Exp}(\mathrm{Ad}_X\psi)$ $2.5\times10^{-111}$. $\mathrm{ad}_\tau$ by commutators (basis images),
81 $\tau$ of PL.4: exact; $\mathrm{Ad}_E$ against `mp.expm` of $\mathrm{ad}_\tau$: exact; (b) $\mathrm{ad}^3 + \theta^2\mathrm{ad}$
$1.4\times10^{-111}$ (71 $\tau$), the block of $\mathrm{ad}^2$ exact; (c) $2.0\times10^{-111}$. Script not committed.
**Permanent:** planned, proptest `adjoint_identity_*`; corpus `se2_*` ($\mathrm{Ad}$ by basis images, as `sen3_ad_n*`).

## 5. The Jacobians

The general statements of `lie-groups.md` (LG.3–LG.8 without the block forms of LG.4(c), LG.9(a),(b), LG.11–LG.15)
use only matrix exponentials and products, so they hold for $\mathrm{SE}(2)$ and $\mathrm{SO}(2)$; what is
$\mathrm{SE}_N(3)$-specific (LG.2, LG.4(c), LG.9(c), LG.10, LG.16) is replaced by PL.4–PL.8. In particular
$\mathrm{Exp}(\tau+\delta) = \mathrm{Exp}(\tau)\mathrm{Exp}(J_r\delta) + O = \mathrm{Exp}(J_l\delta)\mathrm{Exp}(\tau) + O$ for
$J_{r,l} = g(\mp\mathrm{ad}_\tau)$ (LG.7–LG.8), $J_l(\tau) = J_r(-\tau) = \mathrm{Ad}_EJ_r(\tau)$ (LG.9(a),(b)), and
the rows of §2.3 hold with the $J$ and $\mathrm{Ad}$ below (LG.13(d) needs $J_r$ invertible, PL.8, and
$\mathrm{Log}\circ\mathrm{Exp} = \mathrm{id}$, PL.5, i.e. $\lvert\theta\rvert<\pi$).

**Proposition PL.7 (`jr`, `jl`; Proposed for `NUMERICS.md` §6).** For every $\tau$,
$$
J_l = I + a\,\mathrm{ad}_\tau + b\,\mathrm{ad}_\tau^2 = \begin{bmatrix}1 & 0\\ \theta b\,\rho - a\,K\rho & \alpha I + \beta K\end{bmatrix},\qquad
J_r = I - a\,\mathrm{ad}_\tau + b\,\mathrm{ad}_\tau^2 = \begin{bmatrix}1 & 0\\ \theta b\,\rho + a\,K\rho & \alpha I - \beta K\end{bmatrix},
$$
with $\det J_l = \det J_r = s(\theta)^2 = 2a$, and $J_l - J_r = 2a\,\mathrm{ad}_\tau$.

*Proof.* LG.7: $J_l = \sum\mathrm{ad}^n/(n+1)!$. By PL.6(b) the odd powers give $\mathrm{ad}\sum(-z)^j/(2j+2)! = a\,\mathrm{ad}$
and the even powers $\ge2$ give $\mathrm{ad}^2\sum(-z)^j/(2j+3)! = b\,\mathrm{ad}^2$. $J_r$: $\mathrm{ad}\to-\mathrm{ad}$.
Insert $\mathrm{ad}$ and $\mathrm{ad}^2$: the coupling block is $-aK\rho + b\theta\rho$, the lower-right block
$(1 - b\theta^2)I + a\theta K = \alpha I + \beta K$ (PL.2(b)). $\det J = \det V = s^2 = 2a$ (block triangular; PL.4,
PL.2(b)). $\square$

**Proposition PL.8 (`jr_inv`, `jl_inv`; Proposed for `NUMERICS.md` §6).** For $\theta\notin2\pi\mathbb Z_{\ne0}$,
$$
J_l^{-1} = I - \tfrac12\mathrm{ad}_\tau + c\,\mathrm{ad}_\tau^2 = \begin{bmatrix}1 & 0\\ \theta c\,\rho + \frac12K\rho & \gamma I - \frac\theta2K\end{bmatrix},\qquad
J_r^{-1} = I + \tfrac12\mathrm{ad}_\tau + c\,\mathrm{ad}_\tau^2 = \begin{bmatrix}1 & 0\\ \theta c\,\rho - \frac12K\rho & \gamma I + \frac\theta2K\end{bmatrix}.
$$
Hence $J_r^{-1} - J_l^{-1} = \mathrm{ad}_\tau$ exactly; the lower-right block of $J_l^{-1}$ is $V(\theta)^{-1}$ and that of
$J_r^{-1}$ is $V(-\theta)^{-1}$ (PL.5).

*Proof.* Take $\hat J = I + \xi\,\mathrm{ad} + c'\,\mathrm{ad}^2$. By the product rule (2) of SO.3 (which uses only
$\mathrm{ad}^3 = -\theta^2\mathrm{ad}$, PL.6(b)), $J_r\hat J = I$ is $\alpha\xi + a z c' = a$, $-a\xi + \alpha c' = -b$, with
determinant $\alpha^2 + a^2z = s^2 = 2a\ne0$ (PL.2(b)). Its solution is SO.8's: $\xi = \frac12$ and
$c' = (a^2 - b\alpha)/2a = c$. A right inverse of a square matrix is the inverse. $J_l^{-1}$: $\mathrm{ad}\to-\mathrm{ad}$
($c$ is even). Blocks: coupling $\pm\frac12K\rho + c\theta\rho$, lower-right $(1 - cz)I\mp\frac\theta2K$, and
$1 - zc = \gamma$. $\square$

**Checked:** mpmath 1.3.0, 110 digits, 81 tangents of PL.4 (the $J$ rows exercise $\theta$ to $20$ and $\rho$ to
$10^6$; $J^{-1}$ rows only $\lvert\theta\rvert<6.3$, $\bigl\lvert\lvert\theta\rvert - 2\pi\bigr\rvert\ge10^{-2}$, plus $\theta = 2\pi - 10^{-3}$ separately). $a,b,c,\alpha,\beta,\gamma$
from their definitions with guard digits. The dense series $\sum(\mp\mathrm{ad}_\tau)^n/(n+1)!$ of the
`mp`-built $\mathrm{ad}_\tau$ against the closed forms: $8.3\times10^{-106}$ (relative to $\max(1,\lVert J\rVert_{\max})$);
$J = I \mp a\,\mathrm{ad} + b\,\mathrm{ad}^2$ $1.3\times10^{-111}$, $J_l - J_r - 2a\,\mathrm{ad}$ $1.5\times10^{-111}$;
$J_l(\tau) = J_r(-\tau)$ exact; $J_l = \mathrm{Ad}_EJ_r$ $2.8\times10^{-111}$; $\det J - s^2$ $1.7\times10^{-111}$.
`mp.inverse` of the series against PL.8: $2.0\times10^{-109}$ (at $2\pi - 10^{-3}$: $5.8\times10^{-107}$, relative to $\max(1,\lVert J^{-1}\rVert_{\max})$);
$J^{-1} = I\mp\frac12\mathrm{ad} + c\,\mathrm{ad}^2$ $2.1\times10^{-111}$; $J_r^{-1} - J_l^{-1} - \mathrm{ad}$ $6.4\times10^{-112}$;
exactly, in $\mathbb Q(\theta,\tan\frac\theta2,\rho)$ (sympy): $J_lJ_l^{-1} = J_rJ_r^{-1} = I$, $J_r^{-1} - J_l^{-1} = \mathrm{ad}_\tau$, $V\cdot V^{-1} = I$,
$\Gamma_2(\theta K) = aI + \theta bK$. The derivative of the operation itself (LG.11, second form; $h = 10^{-30}$ central
differences, $\mathrm{Log}$ by `mp.logm` of the near-identity matrix), 20 $\tau$ ($\theta\in\{0, \pm10^{-10}, \pm0.4, 1.7, -2.6, \pm3.0, 5.0, -6.2\}$, three with
$\rho\sim10^6$ at $\theta = 0.3, -2.9, 10^{-10}$, six uniform on $[-6,6]$): $\mathrm D^R\mathrm{Exp}$ against $J_r$ $1.3\times10^{-62}$, $\mathrm D^L\mathrm{Exp}$ against $J_l$
$1.5\times10^{-62}$ (the differencing); $\mathrm{Log}(\mathrm{Exp}\,\tau\,\mathrm{Exp}\,\delta) = \tau + J_r^{-1}\delta + O$ and the left form,
$\theta\in\{0.5, -1.2, 2.9, -3.0\}$, $\delta = h\,d$ with $d\sim N(0,1)^3$ fixed per $\tau$, $h = 10^{-20}, 10^{-21}$: error $\le0.16\,h^2$. Scripts not committed.
**Permanent:** planned, corpus `se2_*` ($J$, $J^{-1}$: `mp.inverse` of the dense series, `PHASE1.md` §4.3); proptests
`jl_is_ad_jr_*`, `jacobians_match_dual_*`, `*_matches_reference`.

## 6. Proposed for `NUMERICS.md` §6

**Proposed for `NUMERICS.md` §6: to be adopted, verified by `Dual` and the corpus, in the Phase 3 SE(2) PR
(`PHASE3.md` §6).** `NUMERICS.md` contains none of the matrices below; this page does not edit it.

Rotation-first, $\tau = [\theta;\rho_x;\rho_y]$, $X = (R, t)$, $E = \mathrm{Exp}\,\tau$; $\alpha,\beta$ of §6, $a,b,c$ of §4,
$\gamma = 1 - \theta^2c$ (PL.1). Entry $(0,0)$ is the rotation.
$$
\mathrm{Ad}_X = \begin{bmatrix}1&0&0\\ t_y & R_{00} & R_{01}\\ -t_x & R_{10} & R_{11}\end{bmatrix},\qquad
\mathrm{ad}_\tau = \begin{bmatrix}0&0&0\\ \rho_y & 0 & -\theta\\ -\rho_x & \theta & 0\end{bmatrix},
$$
$$
J_l = \begin{bmatrix}1&0&0\\ \theta b\rho_x + a\rho_y & \alpha & -\beta\\ \theta b\rho_y - a\rho_x & \beta & \alpha\end{bmatrix},\qquad
J_r = \begin{bmatrix}1&0&0\\ \theta b\rho_x - a\rho_y & \alpha & \beta\\ \theta b\rho_y + a\rho_x & -\beta & \alpha\end{bmatrix},
$$
$$
J_l^{-1} = \begin{bmatrix}1&0&0\\ \theta c\rho_x - \frac12\rho_y & \gamma & \frac\theta2\\ \theta c\rho_y + \frac12\rho_x & -\frac\theta2 & \gamma\end{bmatrix},\qquad
J_r^{-1} = \begin{bmatrix}1&0&0\\ \theta c\rho_x + \frac12\rho_y & \gamma & -\frac\theta2\\ \theta c\rho_y - \frac12\rho_x & \frac\theta2 & \gamma\end{bmatrix}.
$$
Equivalently $J_{l,r} = I \pm a\,\mathrm{ad}_\tau + b\,\mathrm{ad}_\tau^2$ and $J_{l,r}^{-1} = I \mp\frac12\mathrm{ad}_\tau + c\,\mathrm{ad}_\tau^2$
(PL.7–PL.8); $J_l(\tau) = J_r(-\tau)$; $J_r^{-1} - J_l^{-1} = \mathrm{ad}_\tau$. Only $a$, $b$, $c$, $\alpha$ are
evaluated: $\theta b$, $\theta c$ are $\theta$ times a catalogue value and $\gamma = 1 - \theta^2c$ (absolute error $\le5u$, no
relative accuracy near $\lvert\theta\rvert = \pi$ where $\gamma\to0$: harmless norm-wise, not componentwise; a $\gamma$ from its own cot arm would
need `jr_inv_coeff` to return it, an open item, PL.11(b)); **never** $(1-\alpha)/\theta$ or $(1-\gamma)/\theta$ (PL.11(a)).
`SE2::log` evaluates $V^{-1}t$ by the $(\alpha,\beta)$ form of PL.5: it needs only $\alpha$, $\beta$ and is componentwise accurate (PL.11(b)). Domains: $J_{l,r}$ every $\tau$; $J^{-1}$ $\lvert\theta\rvert<2\pi$; the
derivative of $\mathrm{Log}$ is $J^{-1}(\mathrm{Log}\,X)$ only for $\lvert\theta(X)\rvert<\pi$. Two further rows, absent from `NUMERICS.md`: §12 (`SE2::jr_inv`, `jl_inv`: $\lvert\theta\rvert<2\pi$, outside (release) "unspecified finite value" as for
SO(3); `SE2::log`: every $X$, a function of the sign of $R_{10}$ at $\lvert\theta(X)\rvert = \pi$, outside "—") and §14 (`SE2::jr_inv`, closed form: twin = dense inverse of `SE2::jr`).
Action Jacobians (PL.10) for §2.4: SO(2) $[KRp]$ both sides, $\partial/\partial p = R$; SE(2) right
$[\,RKp,\ R\,]$, left $[\,K(Rp+t),\ I\,]$, $\partial/\partial p = R$.

## 7. Cross-checks

**Proposition PL.9 (SE(2) in SE(3); the literature).** Let $\iota(\theta,\rho) = [0;0;\theta;\ \rho_x;\rho_y;0]$
and $\iota_G(R(\theta), t) = (R_z(\theta), (t;0))\in\mathrm{SE}(3)$.

- (a) $\iota_G$ is an injective homomorphism and $\iota_G(\mathrm{Exp}\,\tau) = \mathrm{Exp}(\iota\tau)$.
- (b) $L = \iota(\mathbb R^3)$ (entries $3,4,5$ of the rotation-first $[\varphi;\rho]$) is invariant under $\mathrm{Ad}_{\iota_GX}$,
  $\mathrm{ad}_{\iota\tau}$, $J_{r,l}(\iota\tau)$ (`NUMERICS.md` §5.3), and the restrictions are the matrices of PL.6–PL.7.
- (c) For $\rho\perp\varphi = \theta e_z$, $Q(\rho,\varphi) = a\,\Xi + b\,(W\Xi + \Xi W)$: the $e$ words and $W\Xi W$ of §5.3
  vanish, and the $d$ word collapses to $W^2\Xi + \Xi W^2 = -\theta^2\Xi$ (not $0$), which turns $\frac12\Xi$ into
  $(\frac12 - zd)\Xi = a\,\Xi$; $Q\,e_z = \theta b\,\rho - a\,K\rho$ is PL.7's coupling block.
- (d) Solà et al. (arXiv:1812.01537v9, App. C) print, translation-first ($\tau = [\rho;\theta]$, $[1]_\times = K$),
  $\mathrm{Ad}$ (159), $J_r$, $J_l$ (163)–(164) and the action Jacobian (166): each is $\Pi M\Pi^\top$ (resp. $M\Pi^\top$)
  of the matrix here. The appendix prints **no inverse** of $J_r$ or $J_l$ for SE(2); PL.8 has no source there.

*Proof.* (a) Embedding of $2\times2$ blocks in $3\times3$. (b) For $\psi\in L$, $\mathrm{ad}_{\iota\tau}\psi = [\varphi\times\psi_\varphi;\
\varphi\times\eta + \rho\times\psi_\varphi] = [0;\ \theta e_z\times\eta + \psi_\theta\,\rho\times e_z]$ (SE.4), with $e_z\times\eta = K\eta$
and $\rho\times e_z = -K\rho$: this is PL.6(a); $\mathrm{Ad}$ likewise ($t\times e_z = -Kt$). $J = \sum(\mp\mathrm{ad})^n/(n+1)!$
inherits both. (c) $\varphi\cdot\rho = 0$ in SE.8, remark (ii), leaves $a\,\Xi + b\,(W\Xi + \Xi W)$ (the correction $-(\varphi\cdot\rho)[\cdots]$ vanishes, and
$a = \frac12 - zd$ already holds the $d$ word: $[u]_\times[v]_\times[u]_\times = -(u\cdot v)[u]_\times$ gives $W\Xi W = 0$ and
$W^2\Xi + \Xi W^2 = -\theta^2\Xi$); then
$\Xi e_z = \rho\times e_z = -K\rho$, $W\Xi e_z = \theta\,e_z\times(\rho\times e_z) = \theta\rho$, $\Xi We_z = 0$. $\square$

**Checked:** mpmath 1.3.0, 110 digits, 33 tangents ($\theta\in\{\pm10^{-10}, 0.3, -0.7, 1.9, -2.9, \pm3.1, 4.5, -6, 20\}$, two with
$\rho\sim10^6$ at $\theta = 0.3, -2.9$, 20 uniform on $[-6,6]$). (b),(c) §5.3 built as written with $b,d,e$ from definitions (guard digits),
$J_{l,r}(\iota\tau)$ restricted to entries $(3,4,5)$ against PL.7: $2.5\times10^{-111}$, $1.7\times10^{-111}$; all entries coupling
$L$ to the entries $(1,2,6)$ exactly $0$; $Q = a\,\Xi + b\,(W\Xi + \Xi W)$ (full $3\times3$) $2.5\times10^{-111}$;
$\iota_G\mathrm{Exp}$ against the $4\times4$ `mp.expm`: exact. The words of §5.3 for $\rho\perp\varphi$ (the same 33 $\theta$, fresh $\rho$, 80 digits):
$W\Xi W$ and the $e$ words exactly $0$, $W^2\Xi + \Xi W^2 + \theta^2\Xi$ $1.1\times10^{-80}$ (relative to $\max(1,\lVert\Xi\rVert_{\max})$), and
$\lVert W^2\Xi + \Xi W^2\rVert_{\max}/\lVert\Xi\rVert_{\max} = 400.0 = \theta^2$ at $\theta = 20$. (d) the printed formulas transcribed: $J_r$, $J_l$ on the 62 tangents of the 71 of PL.6 with
$\lvert\theta\rvert>10^{-3}$ (their small-$\theta$ forms cancel), $6.0\times10^{-109}$; $\mathrm{Ad}$ and the action Jacobian
on 20 random $(X, p)$: exact. Scripts not committed. **Permanent:** planned, corpus `sen3_jl_n1`, `sen3_jr_n1` at planar
$\tau$ against `se2_*`; proptest `jacobians_match_dual_*`.

**Proposition PL.10 (action Jacobians; Proposed for `NUMERICS.md` §2.4).** With $\mathrm D^s$ of LG.11 and $p\in\mathbb R^2$:
$\mathrm{SO}(2)$: $\mathrm D^R_R(Rp) = \mathrm D^L_R(Rp) = KRp$, $\partial(Rp)/\partial p = R$. $\mathrm{SE}(2)$
($Xp = Rp + t$, tangent $[\delta\theta;\delta\rho]$): $\mathrm D^R_X = [\,RKp,\ R\,]$, $\mathrm D^L_X = [\,K(Rp + t),\ I\,]$,
$\partial(Xp)/\partial p = R$. They are `NUMERICS.md` §2.4 restricted to $L$ with $p_z = 0$: $-R[p]_\times e_z = RKp$,
$-[Rp + t]_\times e_z = K(Rp + t)$.

*Proof.* $\mathrm{Exp}(\delta)p = p + \delta_\theta Kp + \delta_\rho + O(\lVert\delta\rVert^2)$; right:
$R(p + \delta_\theta Kp + \delta_\rho) + t$; left: $Xp + \delta_\theta K(Xp) + \delta_\rho$. $RK = KR$. $\square$

**Checked:** mpmath 1.3.0, 100 digits, central differences $h = 10^{-30}$ of the operation itself with `mp.expm`,
20 random $(X, p)$: right and left $1.7\times10^{-61}$, $\partial/\partial p$ $2.3\times10^{-81}$; $\mathrm{SO}(2)$ (10):
$1.7\times10^{-61}$ (the differencing). Script not committed. **Permanent:** planned, proptests
`jacobians_match_dual_*` (`act_jacobians`).

## 8. Numerical consequences

**Proposition PL.11.** (a) *Exact arms.* CO.5–CO.7 apply unchanged. $\alpha$ has no cancellation (0/0 only at
$\theta = 0$). $\beta = \theta a$ with $a = 2k^2$ (PL.2(b), CO.7(a)): no cancellation; $(1-\cos\theta)/\theta$ has that of
naive $a$, relative $\approx u\,\theta^{-2}$. $\theta b$ cancels like $b$ ($6u\theta^{-2}$) and $\theta c$ like $c$
($p = 2$, CO.6): both come from `coeffs` (the generated series arm), never from $(1-\alpha)/\theta$ or $(1-\gamma)/\theta$. The reverse
subtraction, $\gamma = 1 - zc$, is harmless where $\gamma\approx1$ but loses relative accuracy where $\gamma\to0$ ((b)). The
odd entries multiply the even coefficient by the signed $\theta$ after the branch, so the sign is exact.

(b) *`Log`'s $V^{-1}$.* $\lVert V^{-1}\rVert_2 = 1/s\in[1,\pi/2]$ on $\lvert\theta\rvert\le\pi$. The $(\alpha,\beta)$ form and $s^{-1}R(-\frac\theta2)$ are
products, quotients and sums of positive terms ($\alpha^2 + \beta^2$) of accurate factors: no cancellation, $\le6.5u$ componentwise (measured). The
$\gamma$ form has the entries $\frac\theta2$ (exact to rounding) and $\gamma$: accurate to $2.8u$ relative when evaluated directly as
$h\,(\cos h/\sin h)$, $h = \frac\theta2$, but **not** as $\gamma = 1 - zc$ (the form of §6 for the blocks of $J^{-1}$). Its absolute
error is $z\,\delta c\le48\,c\,u\le4.9u$ (SE.16(a); $\le3.4u$ measured), and near $\pm\pi$, where $\gamma\to0$ and $zc\to1$, it has no relative
accuracy ($1.8\times10^6u$ within $10^{-2}$ of $\pi$, $5.7\times10^9u$ within $10^{-6}$). Norm-wise this is harmless ($\lVert V^{-1}\rVert\ge1$;
`NUMERICS.md` §11's metric is norm-wise), componentwise it is not: the diagonal of the $2\times2$ block of $J^{-1}$ and $\gamma t$ in $\rho$ would need the
cot form, i.e. `jr_inv_coeff` returning $\gamma$ (open item, index). §6 therefore proposes the $(\alpha,\beta)$ form for `SE2::log`. To first order
$\delta\rho = \partial_\theta\rho\,\delta\theta + {}$rounding, $\lVert\partial_\theta\rho\rVert = \kappa_x\lVert t\rVert$,
$$
\kappa_x(\theta) = \tfrac12\sqrt{\Big(\frac{x\cos x - \sin x}{\sin^2x}\Big)^2 + \Big(\frac x{\sin x}\Big)^2},\quad x = \tfrac{\lvert\theta\rvert}2,
$$
increasing (proof below) from $\frac12$ at $0$ to $0.931$ at $\pi$ (SE.16(a) has $\le1.29$ for SE(3)); $\delta\theta$ is (d) of PL.3.

(c) *Conditioning of $J$.* The singular values of $J_l$ and $J_r$ are $\sigma_+$, $\lvert s\rvert$, $\lvert s\rvert/\sigma_+$ with
$\sigma_+^2 = \big(T + \sqrt{T^2 - 4s^2}\big)/2\ge\max(1, s^2)$, $T = 1 + s^2 + \lVert\rho\rVert^2(a^2 + \theta^2b^2)$; so
$\lVert J^{-1}\rVert_2 = \sigma_+/\lvert s\rvert$ and $\kappa_2(J) = \sigma_+^2/\lvert s\rvert\ge1/\lvert s\rvert$ (LG.16(b)). At
$\theta = 0$, $\kappa_2 = \big((\nu + \sqrt{\nu^2 + 4})/2\big)^2\approx\nu^2$, $\nu = \lVert\rho\rVert/2$ (LG.16(c)). The
closed forms of PL.8 involve no solve, so the coupling block errs like $\lVert\rho\rVert\cdot$(coefficient error);
but the dense-inverse reference twin (§14 style) errs like $\kappa_2u$: at $\lVert\rho\rVert = 10^4$,
$\kappa_2\approx2.5\times10^7$, seven digits, and the `*_matches_reference` tolerance must scale with it.

(d) *Domains.* $\mathrm{Exp}$, $J_{l,r}$: every finite $\theta$ ($z = \theta^2$ overflows for
$\lvert\theta\rvert>1.3\times10^{154}$, `f64`, CO.15; it is subnormal below $1.5\times10^{-154}$ and $0$ below $1.6\times10^{-162}$, where
the series arm still returns $\alpha = 1$, $a = \frac12$, $b = \frac16$, $c = \frac1{12}$ to full relative accuracy, CO.15(b), so $\beta = \theta a$, $\theta b$,
$\theta c$ keep the sign of $\theta$, $-0$ included). NaN and $\pm\infty$ in $\theta$ or $\rho$ are outside every domain here. $J^{-1}$: $\lvert\theta\rvert<2\pi$, with
$\lVert J^{-1}\rVert_2 = \sigma_+/\lvert s\rvert\sim\sigma_+\,2\pi/(2\pi-\lvert\theta\rvert)$ toward $2\pi$. $\mathrm{Log}$: every $X$, jumping by $2\pi$ at
$\lvert\theta(X)\rvert = \pi$ (PL.5); its derivative is $J^{-1}$ only for $\lvert\theta(X)\rvert<\pi$.

*Proof.* (b) $\rho(\theta) = s^{-1}R(-\frac\theta2)t$; $R$ and $K$ commute, so
$\partial_\theta\rho = R(-\tfrac\theta2)\big(-\frac{s'}{s^2}I - \frac1{2s}K\big)t$, and $pI + qK$ is a scaled rotation of ratio $\sqrt{p^2 + q^2}$, so
$\kappa_x = \sqrt{(s'/s^2)^2 + 1/(4s^2)}$ with $s = \sin x/x$, $s' = \frac{x\cos x - \sin x}{2x^2}$; substitute. Increasing: with
$g = (\sin x - x\cos x)/\sin^2x > 0$ on $(0,\pi)$, $4\kappa_x^2 = g^2 + (x/\sin x)^2$ and $(x/\sin x)' = g$; $g' = N/\sin^3x$ with
$N = x(1 + \cos^2x) - \sin2x$, $N(0) = 0$, $N' = \sin x\,(3\sin x - 2x\cos x) > 0$ on $(0,\pi/2]$ since $\sin x > x\cos x$: so $g$,
$x/\sin x$ and $\kappa_x$ increase on $x = \lvert\theta\rvert/2\in(0,\pi/2]$. (c) $J = \begin{bmatrix}1&0\\C&V\end{bmatrix}$, $V = sR(\pm\frac\theta2)$ (PL.4),
$C = (\theta bI \mp aK)\rho$ with $\lVert C\rVert = \lVert\rho\rVert\sqrt{a^2 + \theta^2b^2}$ ($\theta bI - aK$ is a scaled rotation).
$\mathrm{diag}(1, R(\pm\frac\theta2))^\top J$, then $\mathrm{diag}(1, Q)^\top(\cdot)\,\mathrm{diag}(1, Q)$ for a rotation $Q$ of the plane, reduce $J$ to
$\begin{bmatrix}1&0&0\\ \lVert C\rVert & s & 0\\ 0&0&s\end{bmatrix}$: singular value $\lvert s\rvert$ and those of the $2\times2$ block, whose
product is $\lvert s\rvert$ and squares sum to $T$. $\sigma_+\ge1$ since $T\ge1 + s^2$. $\square$

**Checked:** (a),(b) `f64` (CPython floats, glibc `sin`, `cos`; the `libm` crate may differ by an ulp), reference
mpmath 1.3.0 at 150 digits from the same doubles, 3000 log-uniform $\theta$ per decade, maximum relative error in units of
$u$: on $[10^{-4},10^{-3}]$, $[10^{-2},10^{-1}]$, $[10^{-1},1]$: $\alpha$ $1.5, 1.4, 1.4$; $\beta = 2\sin^2(\theta/2)/\theta$ $2.9, 3.4, 3.3$;
$\beta = (1-\cos\theta)/\theta$ $9.1\times10^7, 9.1\times10^3, 98$ (all digits lost below $\sim10^{-8}$: $9.0\times10^{15}$ on $[10^{-8},10^{-7}]$);
$\theta b = (\theta-\sin\theta)/\theta^2$ $3.8\times10^8, 4.5\times10^4, 3.7\times10^2$; $\theta c = (1-\gamma)/\theta$
$1.8\times10^9, 2.2\times10^5, 1.8\times10^3$; $\gamma = \frac\theta2\cot\frac\theta2$ as $h\,(\cos h/\sin h)$, $h = \frac\theta2$, $2.2, 2.3, 2.2$. $\gamma$ on
six windows (3000 samples each; $[0.5,1]$, $[2,3]$, $[3,\pi - 10^{-2}]$, $[\pi - 10^{-2},\pi]$, $[\pi - 10^{-4},\pi]$, $[\pi - 10^{-6},\pi]$), maximum
relative error in $u$: as $h\,(\cos h/\sin h)$ $2.2, 2.8, 2.7, 2.8, 2.6, 2.4$, as $h/\tan h$ $\le1.9$, as $1 - \theta^2c$ ($c$ by
$\theta^{-2} - \cos h/(2\theta\sin h)$) $3.5, 12, 131, 1.8\times10^6, 4.1\times10^7, 5.7\times10^9$, with absolute error $\le3.4u$ throughout. $V^{-1}$ norm-wise
(Frobenius), 4000 samples each, both signs, against 200 digits: $\theta\in[10^{-8},10^{-2}]$: the $(\alpha,\beta)$ form $2.2$, the $\gamma$ form
(cot) $2.3$; $[10^{-2},3]$ and $[\pi - 10^{-2},\pi]$ (re-run with the third form): $(\alpha,\beta)$ $3.4, 2.6$, $\gamma$ (cot) $2.0, 0.01$,
$\gamma = 1 - \theta^2c$ $3.0, 1.0$ (exact arm of $c$ throughout; below its switch the series arm applies, not measured here). The
$(\alpha,\beta)$ form componentwise (4000 samples; $[10^{-8},10^{-2}]$ log-uniform, $[10^{-2},3]$, $[\pi - 10^{-2},\pi]$, $[\pi - 10^{-6},\pi]$): $5.5, 5.5, 6.5, 5.1$. $\kappa_x$: `mp.diff` of $V^{-1}(\theta)t$ against the
formula on 2000 grid points of $(0,\pi]$, two random $t$ each: $3.1\times10^{-61}$ relative, increasing, $\kappa_x(\pi) = 0.931048$; sympy: $g'\sin^3x = N$ and $N' = \sin x\,(3\sin x - 2x\cos x)$ exactly.
(c) 100 digits, `mp.svd_r` of $J_l$, $J_r$ at $\theta\in\{0, 10^{-8}, 0.5, 2, 3, \pi, 4.5, 6.2\}$, $\rho\sim N(0,c^2)$,
$c\in\{1, 100, 10^6\}$: the three singular values $2.3\times10^{-105}$, $\kappa_2$ likewise; at $\theta = 0$, $\rho\sim10^6$:
$\kappa_2 = 8.6402\times10^{11} = (\lVert\rho\rVert/2)^2$ to the digits printed. Scripts not committed. **Permanent:** none for
(a)–(c) beyond the strata that measure the quantities (`theta:*`, `rho:*` for `se2_*`); (d) is `NUMERICS.md` §12 (missing rows: index, open items).

## 9. Sign audit

Where a sign or a side is easy to lose (`lie-groups.md` §6's method): the correct form, the usual wrong one, and
the exact difference between them.

| Where | Correct | Wrong variant | The two differ by | Reason |
|---|---|---|---|---|
| $\mathrm{Ad}_X$ coupling | $-Kt$, block $(1,0)$ | $+Kt$; or upper-right (translation-first order, $\Pi\,\cdot\,\Pi^\top$) | $2Kt$ in block $(1,0)$, norm $2\lVert t\rVert$, zero iff $t = 0$ | PL.6: rotating about the origin moves a body at $t$ by $t\times e_z = -Kt$ |
| $J_l$ vs $J_r$ | $J_l = I + a\,\mathrm{ad} + b\,\mathrm{ad}^2$ | swapped | $2a\,\mathrm{ad}_\tau$, $2a = 1 - \theta^2/12 + \cdots$ | PL.7; LG.9(a) |
| $J^{-1}$, the $\frac12$ term | $J_l^{-1}$: $+\frac12K\rho$; $J_r^{-1}$: $-\frac12K\rho$ | swapped | the sign of the $\frac12$ alone: $K\rho$ in block $(1,0)$; the whole swap: $\mathrm{ad}_\tau$ exactly | PL.8 |
| $V^{-1}$ in `Log` | $\gamma I - \frac\theta2K$ | $+\frac\theta2K$, i.e. $V^{-\top}$ | $\theta K$, hence $\theta Kt$ in $\rho$ | PL.5: $V^{-1} = s^{-1}R(-\theta/2)$ |
| $V$ in `Exp` | $\alpha I + \beta K$ | $\alpha I - \beta K = V^\top$: chord turned by $-\theta/2$ | $2\beta K$ | PL.4 |
| $\ominus$ w.r.t. the subtrahend | R: $-J_l^{-1}$; L: $-J_r^{-1}$ | $-J_r^{-1}$ (R) | $\mathrm{ad}_\tau$ | LG.14: $\mathrm{Exp}(-\delta)$ flips sign and side |

**Checked:** mpmath 1.3.0, 110 digits, $\tau = (0.7, 0.3, 1.2)$, $X = (R(0.9), (1.0, -2.0))$: the $\mathrm{Ad}$ variants differ by
$[\,4.0,\ 2.0\,]$ in block $(1,0)$ $= 2Kt$, all else $0$; $J_l - J_r - 2a\,\mathrm{ad}_\tau$ $4.2\times10^{-112}$
($2a = 0.9598$); $J_l^{-1}$ with the sign of the $\frac12$ flipped differs by $[1.2, -0.3] = -K\rho$ in block $(1,0)$ ($\rho = (0.3, 1.2)$), all else $0$;
the $V^{-1}$ and $V$ variants differ by $\theta K$, $2\beta K$ exactly. Script not committed. **Permanent:** as PL.7, PL.8.
