# Geodesics: the group geodesic, its Jacobians, invariance and the dual-quaternion power

> Non-normative companion to [`NUMERICS.md`](../NUMERICS.md) §10 (geodesics), §12 (the `geodesic` domain),
> §14 (its two twins) and [`PHASE4.md`](../PHASE4.md) §1–§3. **`NUMERICS.md` wins on any conflict; a
> disagreement is an open item in the [maths index](./index.md)**, which holds the notation and the `Checked:`
> convention. $\mathrm{Exp}$, $\mathrm{Ad}$ and the sides: [`lie-groups.md`](./lie-groups.md); the quaternion:
> [`so3.md`](./so3.md); the dual-matrix algebra: [`se3.md`](./se3.md); the coefficients:
> [`coefficients.md`](./coefficients.md).

Throughout $G = \mathrm{SE}_N(3)$ ($N = 0$: SO(3)), $\Delta = X_0^{-1}X_1$, $d = \mathrm{Log}\,\Delta = [\varphi;\rho_i]$,
$E_t = \mathrm{Exp}(t\,d)$, $\theta = \theta(\Delta) = \lVert\varphi\rVert$, $\alpha = \theta/2$ (as in SO.5), and
$X(t) = \gamma(X_0, X_1, t)$. Jacobians are right-sided unless a side is named. Symbols specific to this page are
in the index. In §4 the interpolation parameter is $s$, because $t$ is a translation there; in §6–§7 the pose
translation is $\mathbf t$.

## Results

| Label | Result | `NUMERICS.md`; `PHASE4.md` (P4) | Implemented by |
|---|---|---|---|
| GE.1–GE.2 | $X(t) = X_0\,\mathrm{Exp}(t\,\mathrm{Log}(X_0^{-1}X_1))$ for $\theta(\Delta) < \pi$; endpoints, symmetry, constant body **and** spatial velocity, uniqueness among constant-velocity curves with $\theta < \pi$ | §10, §12; P4 §1.1 | `LieGroup::geodesic`, `reference::geodesic` |
| GE.3 | "geodesic" = autoparallel of $\nabla_AB = \tfrac12[A,B]$ (every $N$) = Levi-Civita curve of every *nondegenerate* $\mathrm{Ad}$-invariant form: definite on SO(3), indefinite $(3,3)$ on SE(3), none on $\mathrm{SE}_N(3)$, $N \ge 2$ (rank $\le 6 < 3 + 3N$); no Riemannian length is claimed for $N \ge 1$ | §10 | — |
| GE.4 | left-invariance (every group), right-invariance and inversion (through $\mathrm{Log}\,H^{-1}\Delta H = \mathrm{Ad}_{H^{-1}}d$) | §10; P4 §3 | the invariance tests |
| GE.5 | `Product<SO3, R3>`: (slerp, lerp); bi-invariant for its own law; read as SE(3) poses, left- but not right-invariant, with the exact gap $\lvert\mu_s(\theta)\rvert\lVert t_{H\perp}\rVert$ | §10; P4 §1.3, §3 | `Product` (provided `geodesic`) |
| GE.6–GE.7 | $J_1 = t\,J_r(td)J_r^{-1}(d)$, $J_0 = \mathrm{Ad}_{\mathrm{Exp}(-td)} - J_1\mathrm{Ad}_{\Delta^{-1}}$, $\partial_t = d$; the equal cancellation-free $J_0 = (1-t)J_l((1-t)d)J_l^{-1}(d)$; boundary values; left forms, $J_0^L$ cancellation-free too; domain | §10, §2.3; P4 §2 | `geodesic_jacobians` |
| GE.8 | velocity $d$ and its Jacobians | §10; P4 §2 | `geodesic_velocity` |
| GE.9–GE.10 | unit dual quaternions: $\hat q_T\hat q_{T'} = \pm\hat q_{TT'}$; screw parameters $(\theta, \hat n, h, m)$; $d = [\theta\hat n;\ h\hat n + \theta m]$ | P4 §1.2 | (the fast twin's algebra) |
| GE.11–GE.12 | $\hat q^{\,t} = \exp(t\Lambda)$ is the dual quaternion of $\mathrm{Exp}(td)$; the sign rule; the grouped one-`atan2`, one-`sin_cos` evaluation | P4 §1.2 | `SE3::geodesic` (fast twin) |
| GE.13 | rounding (`f64`, given $\hat q_\Delta$): no cancellation in the value at small angle or near $\pi$; the $0/0$ arm and its per-type guard; the derivative through `Dual` (the translation loses $\approx10^2\alpha^{-1}u$); the endpoints; the branch at $\pi$ | §4, §12; P4 §1.2 | `coeffs`, `S::branch` |
| GE.14 | SO(3): the same formula without the dual part is slerp | P4 §1.3, §5.1 | `SO3::geodesic` |

## 1. The curve

**Definition GE.1.** For $\theta(\Delta) < \pi$ and $t \in \mathbb R$,
$\gamma(X_0, X_1, t) = X_0\,\mathrm{Exp}(t\,\mathrm{Log}(X_0^{-1}X_1))$ (`NUMERICS.md` §10, the reference twin). By LG.2(c),
$d$ is the **unique** $\xi$ with $\theta(\xi) < \pi$ and $\mathrm{Exp}\,\xi = \Delta$. The excluded set
$\theta(\Delta) = \pi$ is a hypersurface: the domain is open and dense.

**Proposition GE.2.** For $\theta(\Delta) < \pi$:

- (a) $X(0) = X_0$ and $X(1) = X_1$.
- (b) $X^{-1}\dot X = d^\wedge$ and $\dot XX^{-1} = (\mathrm{Ad}_{X(t)}d)^\wedge$ with $\mathrm{Ad}_{X(t)}d = \mathrm{Ad}_{X_0}d$: the body
  **and** the spatial twist are constant. $X$ is the only curve from $X_0$ to $X_1$ with constant body velocity
  $\xi$, $\theta(\xi) < \pi$.
- (c) $\gamma(X_1, X_0, 1 - t) = \gamma(X_0, X_1, t)$.

*Proof.* (a) $\mathrm{Exp}\,0 = I$, $\mathrm{Exp}\,d = \Delta$. (b) $\frac{d}{dt}e^{tA} = e^{tA}A$ gives the body twist;
$\mathrm{Ad}_{E_t}d = e^{t\,\mathrm{ad}_d}d = d$ (LG.5, $\mathrm{ad}_dd = 0$), so
$\mathrm{Ad}_{X(t)}d = \mathrm{Ad}_{X_0}\mathrm{Ad}_{E_t}d = \mathrm{Ad}_{X_0}d$. A curve with constant body velocity $\xi$ is
$X_0\exp(t\xi^\wedge)$ (a linear ODE), and $X(1) = X_1$ means $\mathrm{Exp}\,\xi = \Delta$, so $\xi = d$. (c)
$\mathrm{Log}(\Delta^{-1}) = -d$ (LG.2(c)), and $X_1\mathrm{Exp}(-(1-t)d) = X_0\Delta\,\mathrm{Exp}(-d)\mathrm{Exp}(td) = X_0E_t$. $\square$

**Checked:** mpmath 1.3.0, 110 digits, $N = 0..3$; 14 random $(X_0, d)$ per $N$ ($\theta(X_0) \le 0.99\pi$,
translations $\mathcal N(0, 2^2)$; $\theta(d)$ cycling $10^{-10}, 0.3, 1.1, 2.5$, $U(0.1, 3)$, $\pi - 10^{-3}$, $\pi - 10^{-6}$;
$X_1 = X_0\mathrm{Exp}\,d$), $t \in \{0.3, -0.4, 1.7\}$. $\gamma$ from `mp.expm` and the geometric $\mathrm{SO}(3)$
$\mathrm{Log}$ (quaternion $\mathrm{atan2}$, SO.5) with a $\Gamma_1$ solve. Errors absolute, or relative to the largest
entry above $1$. $\mathrm{Log}\,\mathrm{Exp}\,d = d$: $7.0\times10^{-111}$; (a) $1.0\times10^{-110}$; (b)
$X^{-1}\dot X - d^\wedge$ by central differences of $\gamma$ itself (step $10^{-30}$) $1.6\times10^{-60}$, the differencing;
$\mathrm{Ad}_{X(t)}d - \mathrm{Ad}_{X_0}d$ $4.0\times10^{-111}$ (15 pairs, $t \in \{0.3, -0.7, 1.9\}$); (c)
$9.8\times10^{-111}$. For $\theta(d) = 1.2\pi$ (3 pairs per $N$) $\gamma(\tfrac12)$ differs from $X_0\mathrm{Exp}(d/2)$ by $0.32$
to $1.96$, though both reach $X_1$: the excluded curves are real constant-velocity curves, not the one $\mathrm{Log}$
returns. Script not committed. **Permanent:** planned, corpus `so3_geodesic`, `se3_geodesic` (`PHASE4.md` §4) for
the value; none specified for (b), (c).

## 2. What "geodesic" means

Let $\nabla$ be the connection with $\nabla_AB = \tfrac12[A,B]$ on left-invariant fields (the canonical symmetric
connection of Cartan–Schouten: torsion-free, since $\tfrac12[A,B] - \tfrac12[B,A] - [A,B] = 0$; on a matrix group it is
invariant under both translations and inversion: standard, and not used below, since GE.4 proves both invariances directly).

**Proposition GE.3.**

- (a) For a curve $g(t)$ with body velocity $\xi = g^{-1}\dot g$, $\nabla_{\dot g}\dot g \leftrightarrow \dot\xi$: the
  $\nabla$-geodesics are exactly $g_0\exp(t\xi^\wedge)$, the curves of GE.1, for every $N$.
- (b) If $B$ is a nondegenerate symmetric form on $\mathfrak g$ with $B([X,Y],Z) = B(X,[Y,Z])$ ($\mathrm{Ad}$-invariant),
  $\nabla$ is the Levi-Civita connection of the bi-invariant (pseudo-)metric $B$.
- (c) $\mathrm{SO}(3)$, $\mathbb R^n$ and their direct products have positive-definite invariant forms
  ($\varphi\cdot\varphi'$, $\rho\cdot\rho'$): the curves are Riemannian geodesics, and on $\mathrm{SO}(3)$ for
  $\theta < \pi$ the unique shortest ones, of length $\theta$ (standard [Milnor], not re-proved).
- (d) Every invariant symmetric form $B$ on $\mathfrak{se}_N(3)$ vanishes on every pair of translation blocks, so
  $B = \begin{bmatrix}C & \Gamma\\ \Gamma^\top & 0\end{bmatrix}$ ($C$ is $3\times3$, $\Gamma$ is $3\times3N$) has rank $\le 6$: none is
  positive definite, and for $N \ge 2$ none is nondegenerate ($6 < 3 + 3N$). **$\mathrm{SE}_N(3)$, $N \ge 1$, has no
  bi-invariant Riemannian metric, and $\mathrm{SE}_N(3)$, $N \ge 2$, no bi-invariant pseudo-metric** (so (b) has no instance
  there). The curves of GE.1 are geodesics of a left-invariant Riemannian metric $\xi^\top\mathbb I\xi'$ for every
  $(X_0, \xi)$ iff $\mathbb I$ is invariant: never. For $N = 1$ the invariant forms are the two-parameter family
  $c\,\omega\cdot\omega' + \beta(\omega\cdot v' + \omega'\cdot v)$; for $\beta \ne 0$ each has signature $(3,3)$ (the Klein
  form is $c = 0$, $\beta = 1$), for $\beta = 0$ it is degenerate.

**What is claimed for $\mathrm{SE}_N(3)$, $N \ge 1$:** $X(t)$ is the constant-screw motion, the $\nabla$-autoparallel,
invariant under both translations (GE.4). For $N = 1$ it is also the geodesic of every nondegenerate invariant form (all of
signature $(3,3)$); for $N \ge 2$ no bi-invariant pseudo-metric exists and nothing beyond the autoparallel and GE.4 is
claimed. It is **not** claimed to minimise a length for any Riemannian metric.

*Proof.* (a) In the global frame of left-invariant fields, $\dot g = \sum\xi^iE_i$ and
$\nabla_{\dot g}\dot g = \sum\dot\xi^iE_i + \sum\xi^i\xi^j\nabla_{E_i}E_j = \dot\xi + \tfrac12[\xi,\xi] = \dot\xi$.
(b) Koszul: $2B(\nabla_AB', C) = B([A,B'],C) - B([B',C],A) + B([C,A],B')$; invariance gives $B([B',C],A) = B(B',[C,A])$,
so the last two terms cancel. (d) Differentiating $B(\mathrm{Ad}_X\psi, \mathrm{Ad}_X\psi') = B(\psi,\psi')$ gives
$B(\mathrm{ad}_\tau\psi,\psi') + B(\psi,\mathrm{ad}_\tau\psi') = 0$. With $\tau = [0;\rho]$ (one block), $\psi = [a;0]$, $\psi' = [0;\eta]$
(any block), SE.4 gives $\mathrm{ad}_\tau\psi = [0;\rho\times a]$, $\mathrm{ad}_\tau\psi' = 0$, so
$B([0;\rho\times a],[0;\eta]) = 0$, and $\rho\times a$ ranges over $\mathbb R^3$: $B$ has the block form above, of rank at most
$\mathrm{rank}\,[C\ \ \Gamma] + \mathrm{rank}\,\Gamma^\top \le 3 + 3$. For $N = 1$, $\tau = [w;0]$ has $\mathrm{ad}_\tau = \mathrm{diag}(W, W)$ and invariance gives
$C[w]_\times = [w]_\times C$ and $\Gamma[w]_\times = [w]_\times\Gamma$ for every $w$, so $C = cI$, $\Gamma = \beta I$ ($\mathfrak{so}(3)$ acts
irreducibly on $\mathbb R^3$); each such form is invariant (SE.4, $\beta[(\rho\times a)\cdot a' + a\cdot(\rho\times a')] = 0$). For a left-invariant metric $\mathbb I$ the geodesic
equation in body velocity is $\mathbb I\dot\xi = \mathrm{ad}_\xi^\top\mathbb I\xi$ (Euler–Poincaré [Arnold], quoted): testing with
$\eta$, $\dot\xi = 0$ for all $\xi$ iff $T(\xi,\eta,\xi) = 0$, $T(\xi,\eta,\zeta) = \langle[\xi,\eta],\zeta\rangle_{\mathbb I}$; polarising and
using $T(\xi,\eta,\zeta) = -T(\eta,\xi,\zeta)$ gives $T(\xi,\eta,\zeta) = T(\eta,\zeta,\xi)$, which is invariance. $\square$

**Checked:** mpmath 1.3.0, 110 digits, $N = 1$, 20 random $(X, \tau, \psi, \psi')$. The Klein form is
$\mathrm{Ad}_X$-invariant (LG.3's $\mathrm{Ad}$) to $2.7\times10^{-110}$ and $\mathrm{ad}$-skew exactly; the Euclidean inner
product is not ($\lVert\mathrm{Ad}^\top\mathrm{Ad} - I\rVert_{\max}$ up to $60$); $\mathrm{ad}_{[0;\rho]}[a;0] - [0;\rho\times a]$
$8.3\times10^{-112}$, $\mathrm{ad}_{[0;\rho]}[0;\eta] = 0$; the Euler–Poincaré acceleration $\mathrm{ad}_\xi^\top\xi$ (unit
$\mathbb I$) has max entry $5.5$ on $\mathfrak{se}(3)$ and $0$ on $\mathfrak{so}(3)\oplus\mathbb R^3$. (d) The space of symmetric $B$ with
$\mathrm{ad}_{e_k}^\top B + B\,\mathrm{ad}_{e_k} = 0$ for every basis vector (sympy 1.14.0, exact, $\mathrm{ad}$ from the commutators of the hat matrices) has
dimension $N + 1$ and maximal rank $3$ of $3$ ($N = 0$), $6$ of $6$ ($N = 1$), $6$ of $9$ ($N = 2$), $6$ of $12$ ($N = 3$). (c) and the
Euler–Poincaré equation are standard and only outlined. Scripts not committed. **Permanent:** none specified.

## 3. Invariance

**Proposition GE.4.** For $\theta(\Delta) < \pi$:

- (a) *Left, every group:* $\gamma(HX_0, HX_1, t) = H\,\gamma(X_0, X_1, t)$.
- (b) *Right, $\mathrm{SE}_N(3)$:* $\gamma(X_0H, X_1H, t) = \gamma(X_0, X_1, t)\,H$, equivalently
  $X_0H\,\mathrm{Exp}(t\,\mathrm{Ad}_{H^{-1}}d) = X_0\mathrm{Exp}(t\,d)\,H$ (`NUMERICS.md` §10).
- (c) *Inversion:* $\gamma(X_0^{-1}, X_1^{-1}, t) = \gamma(X_0, X_1, t)^{-1}$.

*Proof.* (a) $(HX_0)^{-1}HX_1 = \Delta$. (b) The relative element becomes $H^{-1}\Delta H$; conjugation preserves the
rotation angle ($R_H^\top RR_H$), and so does $\mathrm{Ad}_{H^{-1}}$ (rotation block $R_H^\top\varphi$). LG.4(a) gives
$\mathrm{Exp}(\mathrm{Ad}_{H^{-1}}d) = H^{-1}\Delta H$, so $\mathrm{Log}(H^{-1}\Delta H) = \mathrm{Ad}_{H^{-1}}d$ by LG.2(c), and
$X_0H\,\mathrm{Exp}(t\,\mathrm{Ad}_{H^{-1}}d) = X_0H\,H^{-1}E_tH$. (c) $X_0X_1^{-1} = X_0\Delta^{-1}X_0^{-1}$, so
$\mathrm{Log} = \mathrm{Ad}_{X_0}(-d)$ as in (b), and $X_0^{-1}\mathrm{Exp}(-t\,\mathrm{Ad}_{X_0}d) = E_t^{-1}X_0^{-1}$. $\square$

The proof of (b) uses only that $\mathrm{Log}$ commutes with $\mathrm{Ad}$ on $\theta < \pi$: it holds for every group here whose
conjugation preserves $\theta$, **`Product<SO3, R3>` with its own law included** (GE.5). Right-invariance is not what
separates SE(3) from `tf2`'s interpolation; the composition law is.

**Checked:** as GE.2, plus random $G$ ($\theta \le 0.99\pi$, translations $\mathcal N(0, 3^2)$) and $H$ ($\theta \le \pi$,
$\mathcal N(0, 5^2)$): (a) $3.6\times10^{-111}$; (b) $2.5\times10^{-110}$, the tangent form $1.1\times10^{-110}$; (c)
$1.8\times10^{-110}$ ($N = 0..3$, 14 pairs, 3 values of $t$). Script not committed. **Permanent:** planned, the
invariance tests of `PHASE4.md` §3 (unnamed there) for (a), (b); none for (c).

## 4. `Product<SO3, R3>`

`Product` composes componentwise, $(R,t)(R',t') = (RR', t + t')$ (the `Rn`, `Product` row of `PHASE3.md` §0.0: "`Mul` and `Blend`
componentwise"; §7 calls `Product<SO3, Rn<3>>` "the tf2-semantics pose" and does not state the law of $a\cdot H$), so
$\mathrm{Exp}(\varphi,\rho) = (\mathrm{Exp}\,\varphi, \rho)$, $\mathrm{Log}(R,t) = (\mathrm{Log}\,R, t)$ and, for $a = (R_0,t_0)$,
$b = (R_1,t_1)$,
$\gamma(a,b,s) = \big(R_0\mathrm{Exp}(s\,\mathrm{Log}(R_0^\top R_1)),\ t_0 + s(t_1 - t_0)\big)$: **slerp** (GE.14) and **lerp**,
which is `tf2`'s semantics without `tf2`'s small-angle fallback (`PHASE4.md` §0, §1.3). That reading is the statement
of `NUMERICS.md` §10; it is not compared here with `tf2`'s code.

**Proposition GE.5.** Let $E = R_0^\top R_1$, $\hat n$ its axis, $E_s = \mathrm{Exp}(s\,\mathrm{Log}\,E)$.

- (a) With its own law $\gamma$ is left- and right-invariant (GE.4's proof: $\mathrm{Ad}_H = \mathrm{diag}(R_H, I)$ preserves
  $\theta$); it is the Riemannian geodesic of the bi-invariant metric $\lVert\varphi\rVert^2 + \lVert\rho\rVert^2$ (GE.3(c)).
- (b) Read $(R,t)$ as an $\mathrm{SE}(3)$ pose, $a\cdot H = (R_aR_H,\ R_at_H + t_a)$, $H\cdot a = (R_HR_a,\ R_Ht_a + t_H)$.
  Then $\gamma$ is left-invariant and its rotation part right-invariant, but
  $$
  \mathrm{transl}\,\gamma(aH, bH, s) - \mathrm{transl}\big(\gamma(a,b,s)H\big) = R_0\,M_s\,t_H,\qquad M_s = (1-s)I + sE - E_s,
  $$
  of norm $\lvert\mu_s(\theta)\rvert\,\lVert t_{H\perp}\rVert$, with $t_{H\perp}\perp\hat n$ the part of $t_H$ orthogonal to the axis,
  $\mu_s(\theta) = (1-s) + s\,e^{i\theta} - e^{is\theta}$ and $\lvert\mu_s\rvert = \tfrac12s(1-s)\theta^2 + O(\theta^4)$. It vanishes
  iff $s \in \{0,1\}$, $\theta = 0$ or $t_H\parallel\hat n$.
- (c) The body twist of (b) is $[\varphi;\ R(s)^\top(t_1 - t_0)]$: not constant, the translation rotates with the body.

*Proof.* (a) GE.4. (b) $\mathrm{transl}\,\gamma(aH,bH,s) = (1-s)(R_0t_H + t_0) + s(R_1t_H + t_1)$ and
$\mathrm{transl}(\gamma H) = R_0E_st_H + (1-s)t_0 + st_1$; subtract, with $R_1 = R_0E$. $M_s$ commutes with rotations about
$\hat n$, is $0$ on $\hat n$ ($(1-s) + s - 1$) and is multiplication by $\mu_s$ on the orthogonal plane, a similarity. Left:
$\gamma(Ha,Hb,s)$ has translation linear in $R_Ht_a + t_H$, and the rotation part is (a). (c) $X^{-1}\dot X = (R^\top\dot R, R^\top\dot t)$. $\square$

This is why the right-invariance test of `PHASE4.md` §3 is **supposed to fail** and must not be fixed: the lerp is a world-frame
straight line, whereas right-invariance forces the screw coupling $\rho = J_l^{-1}(\varphi)t$ (SE.3); fixing it gives `ScLerp`,
a different function. **But the test fails only under the $\mathrm{SE}(3)$ composition of (b): with `Product`'s own law
it holds exactly (a), and `PHASE4.md` §3 does not say which law $a\cdot H$ uses** (open item, index). The size of the
failure is $\tfrac12s(1-s)\theta^2\lVert t_{H\perp}\rVert$, so `max_err > 1e-6` at $s = \tfrac12$ needs $\theta^2\lVert t_{H\perp}\rVert \gtrsim 10^{-5}$.

**Checked:** mpmath 1.3.0, 110 digits, 18 random pairs ($\theta \in \{10^{-10}, 10^{-3}, 0.7, 2, 3, \pi - 10^{-6}\}$, three each;
$R_H$, $R_G$ of angle $U(0.2, 3)$, $t_H, t_G \sim \mathcal N(0, 3^2)$), $s \in \{0.25, 0.5, -0.4, 1.6\}$. The geodesic against the
sign-fixed slerp formula (GE.14, $s \in [0,1]$) $4.0\times10^{-111}$; (a) left $5.4\times10^{-111}$, right $7.3\times10^{-111}$
(rotation) and $1.7\times10^{-111}$ (translation); (b) left $1.3\times10^{-110}$, right rotation $7.3\times10^{-111}$, the translation
gap against $R_0M_st_H$ $2.7\times10^{-110}$ and its norm against $\lvert\mu_s\rvert\lVert t_{H\perp}\rVert$ $9.6\times10^{-111}$ (gap up to
$17.5$); $\lvert\mu_{1/2}\rvert = 1.25\times10^{-5}, 1.2497\times10^{-3}, 0.1224, 0.4597, 0.9293$ at $\theta = 0.01, 0.1, 1, 2, 3$.
Script not committed. **Permanent:** planned, the positive failure test of `PHASE4.md` §3, once the law of $a\cdot H$ is fixed.

## 5. Jacobians and velocity

**Proposition GE.6 (`NUMERICS.md` §10).** Let $\theta(d) < \pi$, $t \in \mathbb R$, $E = E_t$. As $(\delta_0,\delta_1,\delta_t) \to 0$,
$$
\gamma(X_0\mathrm{Exp}\,\delta_0,\ X_1\mathrm{Exp}\,\delta_1,\ t + \delta_t) = X(t)\,\mathrm{Exp}\big(J_0\delta_0 + J_1\delta_1 + d\,\delta_t\big) + O(\lVert\cdot\rVert^2),
$$
$$
J_1 = t\,J_r(td)\,J_r^{-1}(d),\qquad
J_0 = \mathrm{Ad}_{\mathrm{Exp}(-td)} - t\,J_r(td)J_r^{-1}(d)\,\mathrm{Ad}_{\Delta^{-1}} = \mathrm{Ad}_E^{-1} - J_1\mathrm{Ad}_{\Delta^{-1}}.
$$

*Proof.* $X_1$: $\Delta \to \Delta\mathrm{Exp}\,\delta_1$, and $\mathrm{Log}(\Delta\mathrm{Exp}\,\delta_1) = d + J_r^{-1}(d)\delta_1 + O$ (LG.13(d),
$\theta < \pi$); by LG.8, $\mathrm{Exp}(td + tJ_r^{-1}\delta_1) = E\,\mathrm{Exp}(J_r(td)\,tJ_r^{-1}(d)\delta_1 + O)$. $X_0$: $\Delta \to
\mathrm{Exp}(-\delta_0)\Delta$ and $\mathrm{Log} = d - J_l^{-1}(d)\delta_0 + O$ (LG.13(d), left form), so by LG.13(b) and LG.8
$$
X_0\mathrm{Exp}(\delta_0)\,\mathrm{Exp}\big(td - tJ_l^{-1}\delta_0\big) = X_0E\,\mathrm{Exp}(\mathrm{Ad}_E^{-1}\delta_0)\,\mathrm{Exp}\big(-tJ_r(td)J_l^{-1}(d)\delta_0\big) + O,
$$
and $J_l^{-1}(d) = J_r^{-1}(d)\mathrm{Ad}_{\Delta^{-1}}$ (LG.9(b), $\mathrm{Exp}\,d = \Delta$). $t$: $X(t + \delta_t) = X(t)\mathrm{Exp}(\delta_td)$
exactly (one subgroup). The three terms add because $\gamma$ is smooth on the domain. $\square$

**Corollary GE.7.** Under the hypotheses of GE.6:

- (a) $J_0 = \mathrm{Ad}_E^{-1} - t\,J_r(td)J_l^{-1}(d) = (1-t)\,J_l\big((1-t)d\big)\,J_l^{-1}(d)$, which is $J_1$ of the swapped pair
  $(X_1, X_0)$ at $1-t$. Equivalently $\mathrm{Ad}_E^{-1}J_l(d) - tJ_r(td) = (1-t)J_l((1-t)d)$.
- (b) $(J_0, J_1) = (I, 0)$ at $t = 0$ and $(0, I)$ at $t = 1$.
- (c) Every factor lies in $\mathcal D_N$ (LG.10): $J_0, J_1$ are `SEn3Jac` values, $J_r^{-1}(d)$ the dual inverse (SE.12), and
  $J_0$, $J_1$ share the evaluations at $d$ ($J_l^{-1} = J_r^{-1}\mathrm{Ad}_{\Delta^{-1}}$). For `Product` they are block diagonal,
  with the $\mathbb R^3$ blocks $(1-t)I$ and $tI$.
- (d) *Left side:* $J_i^L = \mathrm{Ad}_{X(t)}J_i\mathrm{Ad}_{X_i}^{-1}$; $J_1^L = t\,\mathrm{Ad}_{X_0}J_l(td)J_l^{-1}(d)\mathrm{Ad}_{X_0}^{-1}$,
  $J_0^L = I - J_1^L = (1-t)\,\mathrm{Ad}_{X_1}J_r\big((1-t)d\big)J_r^{-1}(d)\mathrm{Ad}_{X_1}^{-1}$ (the second form has no cancellation; it costs one
  evaluation of $J_r$ at $(1-t)d$ that $I - J_1^L$ does not), $\partial_t = \mathrm{Ad}_{X(t)}d = \mathrm{Ad}_{X_0}d$. **`NUMERICS.md` §10 lists the right
  side only**; `PHASE4.md` §2 takes a side, so these are the forms it needs (derived here, not adopted).

*Proof.* (a) By GE.2(c), $\partial X(t)/\partial X_0$ is the derivative with respect to the second argument of the pair $(X_1, X_0)$
at $1-t$, with $d' = -d$; GE.6 gives $(1-t)J_r(-(1-t)d)J_r^{-1}(-d)$ and $J_r(-x) = J_l(x)$ (LG.9(a)). A derivative matrix is unique, so
the two expressions agree; the last identity is that equality times $J_l(d)$. (b) Substitute. (c) LG.10, SE.12, LG.9(b). (d)
$\mathrm{Exp}(\delta)Y = Y\mathrm{Exp}(\mathrm{Ad}_Y^{-1}\delta)$ (LG.4(b)) at $X_i$ and at $X(t)$ turns left perturbations into right ones and back.
Then $\mathrm{Ad}_{X(t)}J_1\mathrm{Ad}_{X_1}^{-1} = \mathrm{Ad}_{X_0}\,\mathrm{Ad}_E\,tJ_r(td)\,J_r^{-1}(d)\mathrm{Ad}_\Delta^{-1}\mathrm{Ad}_{X_0}^{-1}$ with
$\mathrm{Ad}_EJ_r(td) = J_l(td)$ and $J_r^{-1}(d)\mathrm{Ad}_\Delta^{-1} = J_l^{-1}(d)$ (LG.9(b)). Moving both ends by the same left $\delta$ moves
$X(t)$ to $\mathrm{Exp}(\delta)X(t)$ (GE.4(a)), so $J_0^L + J_1^L = I$. As in (a), $J_0^L$ is $J_1^L$ of the swapped pair $(X_1, X_0)$ at $1-t$, with
$d' = -d$ and $J_l(-x) = J_r(x)$. $\square$

*Numerical consequences.* $J_1$ has a vanishing factor and no cancellation. **$J_0$ of `NUMERICS.md` §10 is a difference of two
$O(1)$ matrices that tends to $0$ like $1-t$:** its absolute error is $O(u\lVert\mathrm{Ad}\rVert)$ and its relative error grows
like $1/(1-t)$; the equal form of (a) has none. Measured with both operands rounded to `f64` (SE(3), $\theta(d) = 1$), the difference form
errs by $19u$, $820u$, $1.5\times10^5u$, $9.7\times10^8u$, $1.3\times10^{12}u$ at $1-t = 10^{-1}, 10^{-3}, 10^{-5}, 10^{-9}, 10^{-12}$ (of order
$u/(1-t)$), the form (a) by at most $0.72u$. Shipping the second arrangement is a `NUMERICS.md` §10 edit and a record; $J_0^L = I - J_1^L$ cancels the same
way, and (d) gives its cancellation-free form. *Domain:* $\mathrm{Log}$ is differentiable only for $\theta(\Delta) < \pi$ (LG.2(d)); $J_r^{-1}$ itself stays bounded there
($\lVert J_r^{-1}\rVert_2 = 1/s(\theta) \le \pi/2$ on the rotation block, LG.16) but the linearisation holds only until a perturbation of $\Delta$
crosses $\theta = \pi$ (GE.13(d)). For $t \notin [0,1]$ the formulas hold unchanged in exact arithmetic (extrapolation along the same screw); their rounding error grows with $\lvert t\rvert$ (GE.13(a)).

**Corollary GE.8 (velocity).** `geodesic_velocity`$(X_0, X_1) = d = X_1\ominus_RX_0$, the body twist of GE.2(b) per unit $t$ ($d/\Delta t$ over a
span $\Delta t$, `NUMERICS.md` §10), with $\partial d/\partial X_1 = J_r^{-1}(d)$ and $\partial d/\partial X_0 = -J_l^{-1}(d)$ (the $\ominus$ rows of
LG.14, $\theta(d) < \pi$). The spatial twist $\mathrm{Ad}_{X_0}d/\Delta t$ is constant too; the body twist of `Product` read as poses is not (GE.5(c)).

**Checked (GE.6–GE.8):** mpmath 1.3.0, 110 digits; $\mathrm{Exp}$ by `mp.expm`, $\mathrm{Log}$ as GE.2, $\mathrm{Ad}$ from LG.3, $J_r$, $J_l$ from the
series of $\mathrm{ad}$ (LG.7), inverses by `mp.inverse`. Each derivative is the central difference (step $10^{-30}$) of the geodesic map itself, the
argument perturbed on its side and the output differenced on its side (LG.11): columns for $X_0$, $X_1$ and $t$, both sides. $N = 0,1,2$; one random
$(X_0,d)$ for each $\theta(d) \in \{10^{-10}, 0.5, 2, 3, \pi - 10^{-3}, \pi - 10^{-6}\}$; $t \in \{0, 0.25, 0.5, 1 - 10^{-9}, 1, 1.5, -0.5\}$ (126 cases).
Maximum deviation: right $J_1$ $1.1\times10^{-61}$, right $J_0$ (§10 form) $1.1\times10^{-61}$, left $J_1$ and $J_0$ $1.0\times10^{-61}$, $\partial_t$ (right
against $d$, left against $\mathrm{Ad}_{X(t)}d$) $1.7\times10^{-81}$, all the differencing. The equal forms of GE.7(a), (d) against the §10 form and
each other $\le 6.1\times10^{-110}$; $\mathrm{Ad}_E^{-1}J_l(d) - tJ_r(td) - (1-t)J_l((1-t)d)$ $2.5\times10^{-110}$; the zero blocks of $J_0, J_1$ are
exactly $0$, the diagonal blocks equal to $3.3\times10^{-111}$ ($N = 1, 2$). Near the cut ($\theta(d) = \pi - 10^{-3}$, $t = \tfrac12$,
$X_1 \to X_1\mathrm{Exp}(c[\hat n;0])$): the first-order error is $O(c^2)$ (at most $1.6\times10^{-7}$ for $\lvert c\rvert \le 0.9\times10^{-3}$) and $2.2$ to $3.2$
for $c \ge 1.1\times10^{-3}$, where $\theta(\Delta)$ passes $\pi$. The `f64` figures: one sample ($N = 1$), the two operands rounded from the exact matrices.
The cancellation-free $J_0^L$ of GE.7(d) against $\mathrm{Ad}_{X(t)}J_0\mathrm{Ad}_{X_0}^{-1}$ (form (a)), against $I - J_1^L$ and against the conjugated
§10 form: $5.5\times10^{-80}$ (80 digits, $N = 0,1,2$, $\theta(d) \in \{10^{-6}, 0.5, 2, \pi - 10^{-3}\}$, $t \in \{0.3, 0.999, -0.5, 1.5, 1\}$, 60 cases; algebra
only, the left forms' derivatives being those checked above).
Scripts not committed. **Permanent:** planned, `geodesic_jacobians_match_reference`, `jacobians_match_dual_geodesic` (`PHASE4.md` §2); `PHASE4.md` §4 has
no corpus id for the Jacobians.

## 6. Unit dual quaternions

**Definition GE.9.** A dual quaternion is $\hat q = q_{\mathrm r} + \epsilon q_{\mathrm d}$, $\epsilon^2 = 0$, with
$\hat q\hat q' = q_{\mathrm r}q_{\mathrm r}' + \epsilon(q_{\mathrm r}q_{\mathrm d}' + q_{\mathrm d}q_{\mathrm r}')$ (`NUMERICS.md` §2.2's algebra with quaternions for blocks). A pose $T = (q, \mathbf t)$,
$q$ a unit quaternion, has $\hat q_T = q + \epsilon\,\tfrac12(0,\mathbf t)\,q$; the pose is recovered as $R(q_{\mathrm r})$ and $(0,\mathbf t) = 2q_{\mathrm d}q_{\mathrm r}^*$. The pairs
satisfy $\lVert q_{\mathrm r}\rVert = 1$, $q_{\mathrm r}\cdot q_{\mathrm d} = 0$; $\hat q_T$ and $-\hat q_T$ are one pose.

**Proposition GE.10.**

- (a) $\hat q_T\hat q_{T'} = \pm\hat q_{TT'}$: the product is the composition $T_{ax}T_{xb}$ (the repository's order; this is why $q_{\mathrm d} = \tfrac12\mathbf t\,q$, not $\tfrac12q\,\mathbf t$).
- (b) Let $w \ge 0$, $\theta = 2\,\mathrm{atan2}(\lVert q_{\mathrm v}\rVert, w) \in [0,\pi]$, $\hat n = q_{\mathrm v}/\lVert q_{\mathrm v}\rVert$ ($\theta > 0$), $h = \mathbf t\cdot\hat n$ (the
  displacement along the axis; the screw-theory pitch is $h/\theta$), $o$ a point of the screw axis, $m = o\times\hat n$ (the moment). Then $q_{\mathrm r} = (\cos\alpha, \sin\alpha\,\hat n)$ and
  $$
  q_{\mathrm d} = \big({-\tfrac h2}\sin\alpha,\ \tfrac h2\cos\alpha\,\hat n + \sin\alpha\,m\big),\qquad\text{i.e.}\qquad
  \hat q = \cos\tfrac\Theta2 + \sin\tfrac\Theta2\,\nu,\quad \Theta = \theta + \epsilon h,\ \nu = \hat n + \epsilon m,
  $$
  with $f(a + \epsilon b) = f(a) + \epsilon bf'(a)$.
- (c) $d = \mathrm{Log}\,T = [\theta\hat n;\ h\hat n + \theta m]$, i.e. $\varphi + \epsilon\rho = \Theta\nu$: the Plücker line $(\hat n, m)$ scaled by the dual angle
  (at $\theta = \pi$, $\hat n$ is the sign the quaternion gives, SO.5(d), and $d$ is the preimage `Log` returns).

*Proof.* (a) With $\hat q_T = q + \epsilon\tfrac12(0,\mathbf t)q$ and $\hat q_{T'} = q' + \epsilon\tfrac12(0,\mathbf t')q'$, the dual part of the product is
$\tfrac12q(0,\mathbf t')q' + \tfrac12(0,\mathbf t)qq' = \tfrac12(0,\ R\mathbf t' + \mathbf t)\,qq'$, using $q(0,v) = (0,Rv)\,q$. (b) $Ro + \mathbf t = o + h\hat n$ gives $\mathbf t = h\hat n + (I - R)o$. Now
$\tfrac12(0,(I-R)o)\,q = \tfrac12\big[(0,o)q - q(0,o)\big] = (0,\ o\times q_{\mathrm v}) = (0,\sin\alpha\,m)$ and
$\tfrac12(0,h\hat n)(\cos\alpha,\sin\alpha\,\hat n) = \tfrac h2(-\sin\alpha,\ \cos\alpha\,\hat n)$. Expanding $\cos\frac\Theta2$, $\sin\frac\Theta2$ to first order in
$\epsilon$ and multiplying by $\nu$ ($\hat n\cdot m = 0$) gives the same. (c) By SE.3(b) the translation of $\mathrm{Exp}\,\tau$ is $J_l(\varphi)\rho$; with $\varphi = \theta\hat n$,
$W\hat n = 0$ and $\theta m = -Wo$: $J_l(\varphi)(h\hat n + \theta m) = h\hat n - J_lWo = h\hat n - (e^W - I)o = \mathbf t$. So $h\hat n + \theta m = J_l^{-1}(\varphi)\mathbf t$, the
translation of $\mathrm{Log}\,T$ (SE.3(c), $\theta \le \pi$), with $\varphi = \theta\hat n = \mathrm{Log}\,R$ (SO.5(a)). $\square$

**Checked:** mpmath 1.3.0, 110 digits. (a) 40 random $\mathrm{SE}(3)$ pairs: $\hat q_T\hat q_{T'}$ against $\pm\hat q_{TT'}$ $6.7\times10^{-111}$; the constraints and the
recovered pose $2.0\times10^{-111}$. (b), (c) $\theta \in \{10^{-10}, 0.01, 0.3, 1, 2.5, 3, \pi - 10^{-6}\}$, five poses each ($\mathbf t \sim \mathcal N(0, 3^2)$): $h$ from $\mathbf t\cdot\hat n$ and
from $-2q_{\mathrm d,w}/\lVert q_{\mathrm v}\rVert$ agree to $5.0\times10^{-111}$; the moment from $q_{\mathrm d}$ against $o\times\hat n$, $o$ solved independently from
$(I - R + \hat n\hat n^\top)o = \mathbf t - h\hat n$, $2.2\times10^{-102}$; $Ro + \mathbf t = o + h\hat n$ $1.2\times10^{-101}$; the geometric $\mathrm{Log}\,T$ against
$[\theta\hat n;\ h\hat n + \theta m]$ $8.4\times10^{-111}$. Script not committed. **Permanent:** none specified (GE.11's `se3_geodesic_matches_reference` exercises the consequence).

## 7. The power

**Proposition GE.11.** Let $\hat q = \hat q_\Delta$ ($w \ge 0$), $d = [\varphi;\rho]$, $\Lambda = \tfrac12\big((0,\varphi) + \epsilon(0,\rho)\big) = \tfrac12\Theta\nu$ and
$\hat q^{\,t} = \exp(t\Lambda)$. Then:

- (a) $\hat q^{\,t} = \cos\tfrac{t\Theta}2 + \sin\tfrac{t\Theta}2\,\nu$: $q_{\mathrm r}^t = (\cos t\alpha,\ \sin t\alpha\,\hat n)$,
  $q_{\mathrm d}^t = \big({-\tfrac{th}2}\sin t\alpha,\ \tfrac{th}2\cos t\alpha\,\hat n + \sin t\alpha\,m\big)$ (angle $t\theta$, displacement $th$, the same axis).
- (b) $\hat q^{\,t}$ is the dual quaternion of $\mathrm{Exp}(t\,d)$, so the pose of $X_0\hat q_\Delta^{\,t}$ is $\gamma(X_0,X_1,t)$, for every real $t$.
- (c) *Sign rule.* The same formulas on $-\hat q_\Delta$ ($w < 0$) extract $\theta' = 2\pi - \theta$ about $-\hat n$ with $h' = -h$ and return $X_0\mathrm{Exp}(t\,d')$, the
  other constant-velocity curve of GE.2(b) ($\theta(d') > \pi$): the same $X_1$ at $t = 1$, the long way round in between. The canonicalisation $w \ge 0$ selects $\mathrm{Log}$'s branch.

*Proof.* (a) $\nu\cdot\nu = 1$ ($\lVert\hat n\rVert = 1$, $\hat n\cdot m = 0$), so $\Lambda^2 = -\Theta^2/4$ is a dual scalar and Euler's formula holds over the dual numbers; expand as in GE.10(b).
$\Lambda = \tfrac12\Theta\nu$ is GE.10(c). (b) The unit dual quaternions cover the poses twice, so $t \mapsto \hat q_{\mathrm{Exp}(td)}$ has a continuous lift $Q(t)$ with $Q(0) = 1$, and $Q(t) = \pm\hat q_{\mathrm{Exp}(td)}$. By GE.10(a)
and $\mathrm{Exp}(td)\mathrm{Exp}(sd) = \mathrm{Exp}((t+s)d)$, $Q(t)Q(s) = \varepsilon(t,s)Q(t+s)$ with $\varepsilon = \pm1$; $\varepsilon$ is continuous on the connected plane
$\mathbb R^2$ (a ratio of continuous nonzero dual quaternions) and $\varepsilon(t, 0) = 1$, so $\varepsilon \equiv 1$. $Q$ is then a continuous homomorphism from $\mathbb R$ to the unit group of a finite-dimensional associative algebra,
hence smooth, and $Q(t) = \exp(tQ'(0))$. $Q'(0) = \Lambda$: the rotation part is $(\cos t\alpha, \sin t\alpha\,\hat n)$ (continuous, $1$ at $0$), derivative
$\tfrac12(0,\varphi)$, and the translation $J_l(t\varphi)\,t\rho$ (SE.3(b)) has derivative $\rho$, so $q_{\mathrm d}' = \tfrac12(0,\rho)$. (c) $w' = -\cos\alpha$, so $\alpha' = \mathrm{atan2}(\sin\alpha, -\cos\alpha) = \pi - \alpha$,
$\hat n' = -\hat n$, $h' = \mathbf t\cdot\hat n' = -h$; $(-\hat q)^1$ is the same pose as $\hat q$, but $\alpha' > \pi/2$ means $\theta' > \pi$. $\square$

**Proposition GE.12 (one `atan2`, one `sin_cos`).** For $\hat q_\Delta$ ($w \ge 0$) put $\varkappa = q_{\mathrm d,w}/\lVert q_{\mathrm v}\rVert^2 = -h/(2\sin\alpha)$, $\bar m = q_{\mathrm d,v} + \varkappa\cos\alpha\,q_{\mathrm v} = \sin\alpha\,m$,
$\alpha = \mathrm{atan2}(\lVert q_{\mathrm v}\rVert, w)$ and the slerp weight $\varpi_t = \sin(t\alpha)/\lVert q_{\mathrm v}\rVert$. Then
$$
q_{\mathrm r}^t = \big(\cos t\alpha,\ \varpi_t\,q_{\mathrm v}\big),\qquad
q_{\mathrm d}^t = \big(t\,\varpi_t\,q_{\mathrm d,w},\ \ \varpi_t\,\bar m - t\varkappa\cos(t\alpha)\,q_{\mathrm v}\big),\qquad (0,\mathbf t') = 2\,q_{\mathrm d}^t\,(q_{\mathrm r}^t)^*,
$$
and $X(t) = X_0\,(q_{\mathrm r}^t, \mathbf t')$. The cost is one `sqrt`, one `atan2`, one `sin_cos`, two divisions ($\lVert q_{\mathrm v}\rVert^2$, $\lVert q_{\mathrm v}\rVert$) and products.

*Proof.* Substitute GE.10(b) in GE.11(a): $\sin t\alpha\,\hat n = \varpi_tq_{\mathrm v}$, $\sin t\alpha\,m = \varpi_t\bar m$, $\tfrac{th}2\cos t\alpha\,\hat n = -t\varkappa\cos(t\alpha)q_{\mathrm v}$ (since
$\varkappa q_{\mathrm v} = -\tfrac h2\hat n$) and $-\tfrac{th}2\sin t\alpha = t\varpi_tq_{\mathrm d,w}$. $\square$

**Checked (GE.11, GE.12):** mpmath 1.3.0, 110 digits. The dual-quaternion exponential by its series (multiplication of GE.9) against $\hat q_{\mathrm{Exp}\,\tau}$ for
$\theta \in \{0.1, 1, 3, 4, 6\}$ (beyond $\pi$ too) $5.0\times10^{-111}$. The power against $\gamma$ of the definition ($\mathrm{Exp}$ by `mp.expm`, $\mathrm{Log}$ as GE.2), three ways:
the screw parameters of GE.10 in GE.11(a), the series $\exp(t\Lambda)$, and the grouped formulas of GE.12 (the arithmetic of the exact arm of `tf_tree_math::dualquat::screw_pow`): max
$9.4\times10^{-111}$, $1.3\times10^{-110}$, $1.1\times10^{-110}$. Strata $\theta \in \{10^{-30}, 10^{-10}, 10^{-3}, 0.29, 0.31, 1, 2, 3, \pi - 10^{-3}, \pi - 10^{-6}\}$, each with $\lVert\mathbf t\rVert \sim 2$ and
"far" ($X_0$ and $\mathbf t$ with entries $\mathcal N(0, (10^4)^2)$), two poses each, $t \in \{0, 0.125, 0.7314, 1 - 10^{-9}, 1, 1.5, -0.5\}$. (c) The grouped formulas fed $-q$, without the flip, at
$\theta \in \{0.5, 1, 2, 3\}$, $t \in \{0.25, 0.5\}$: the pose differs from the canonical one by $1.06$ to $12$ (relative to the largest entry) and agrees at $t = 1$ ($7.6\times10^{-111}$).
Script not committed. **Permanent:** planned, `se3_geodesic_matches_reference`, `tf_tree_math_differential` ([`0010`](../decisions/0010-seeded-from-tf-tree-math.md)), corpus `se3_geodesic`
(`geo:consecutive`, `geo:generic`, `geo:near-pi`).

## 8. Rounding and the branch at $\pi$

**Observation GE.13** ($u$ the unit roundoff, `NUMERICS.md` §2.1).

- (a) *Small angle: no cancellation in the value.* As $\alpha \to 0$, $\varkappa \approx -h/(2\alpha)$ and $m \approx \mathbf t_\perp/(2\alpha)$ diverge, but enter only as
  $\varkappa q_{\mathrm v} \to -\tfrac h2\hat n$ and $\bar m \to \tfrac12\mathbf t_\perp$ ($\mathbf t_\perp$ the part of $\mathbf t$ orthogonal to the axis). For an axial translation the two terms of $\bar m$ are each
  $\approx\tfrac h2\cos\alpha$ and cancel to an absolute $O(u\lVert\mathbf t\rVert)$, the error scale of the output translation anyway. **In `f64` the exact arm needs no coefficient of `NUMERICS.md` §4 and no
  series for the value**; the series arm of `tf_tree_math`'s `screw_pow` is a speed choice (no transcendental), as its comment says. A port of that arm needs two series, neither in §4:
  the slerp weight $\varpi_t = t\big[1 + \tfrac{1-t^2}6\alpha^2 + \tfrac{(1-t^2)(7-3t^2)}{360}\alpha^4 + \cdots\big]$, a function of $\alpha^2$ and $t$, and $\alpha^2$ from $h = 1 - \cos\alpha$
  (`theta_sq_from_chord`: eight rational coefficients of the squared $\arcsin$ series, fed $h = \sin^2\alpha/(1+\cos\alpha)$). Each is a catalogue entry, a §4 edit and a record (`0004`). Whether the exact arm
  alone meets `PHASE4.md` §5.3's bench precondition against a `tf_tree_math` whose small-angle arm exists for speed (its `docs/design/fast-path.md`, lever 1b: `ScLerp` $51.6 \to 43.6$ ns) is a bench question.
  *Scope:* these are `f64` figures for the exact arm **given** $\hat q_\Delta$ and $s \in [0.125, 1 - 10^{-9}]$; `f32` is unmeasured. The error grows with $\lvert s\rvert$ outside $[0,1]$ (up to $12.6u$ at
  $s \in \{1.5, -0.5, 3\}$ on the same strata). The composed `geodesic(X0, X1, t)` also rounds $\Delta = X_0^{-1}X_1$ (translation: absolute $O(u\lVert\mathbf t_0\rVert)$; $q_{\mathrm v}$ of a product of
  unit quaternions: absolute $O(u)$, so relative $O(u/\theta)$ at small $\theta$) and the product $X_0\cdot\mathrm{pow}$; on `geo:consecutive` ($\lVert d\rVert \in [10^{-9}, 10^{-3}]$, $\lVert t_0\rVert$ up to $10^4$)
  those can dominate, and the figure here is not the end-to-end one: that is judged on `NUMERICS.md` §11's translation floor.
- (b) *$\theta = 0$.* $\varkappa$ is $0/0$ at $\lVert q_{\mathrm v}\rVert = 0$, and $\lVert q_{\mathrm v}\rVert^2$ underflows to $0$ a little above it (`f64`: $\lVert q_{\mathrm v}\rVert \lesssim 10^{-162}$), so the fast twin needs a second arm there: the reference twin (a pure translation is
  $\mathrm{Exp}$ of $[0;\rho]$), whose coefficients are the catalogue's, $k$, $a$, $b$ for $\mathrm{Exp}$ and $r$, $c$ for $\mathrm{Log}$ (`exp_coeffs`, `jr_coeffs`, `log_ratio`, `jr_inv_coeff`). The
  switch is a range guard on $\lVert q_{\mathrm v}\rVert^2$, not the minimiser of a coefficient's error: **it is not a `0004` quantity as `PHASE4.md` §1.2 is worded**, and it is a constant of the
  scalar type, not `tf_tree_math`'s `f64` literal $10^{-290}$, which `f32` cannot represent (there $\lVert q_{\mathrm v}\rVert^2$ leaves the normal range near $10^{-38}$). No document chooses the `f32` value, and
  `f32` is unmeasured. In `f64` the exact arm is still accurate below $10^{-290}$: translation error $\le 1.7u$ at $\lVert q_{\mathrm v}\rVert \in \{10^{-144}, 2\times10^{-145}, 10^{-146}, 10^{-150}, 10^{-153}, 10^{-155}, 10^{-160}\}$
  ($\lVert q_{\mathrm v}\rVert^2$ from $10^{-288}$ to $10^{-320}$, subnormal at the end), so $10^{-290}$ is a margin, not an accuracy switch. *Safe argument* (`0003`): where the guard fires the exact arm
  receives a safe value (for example $\lVert q_{\mathrm v}\rVert^2$ selected to $1$) and its result is discarded by `S::select`; otherwise a masked lane carries the $0/0$, and under `Dual`
  `sqrt`'s derivative is $\pm\infty$ at $0$ (`NUMERICS.md` §12). The guard's mask comes from `S::branch`'s comparison, never from an `if` on the scalar.
- (c) *Through `Dual`.* The loss is in the translation of the pose, not only in $\varpi_t$. The slerp weight alone has value error $\le 2.2u$ and derivative error $(3.4$ to $6.1)\,\alpha^{-2}u$ relative to its own derivative, which is $O(\alpha)$ (the law of CO.14,
  $p' = 2$); it costs the rotation nothing, since $\varpi_t$ multiplies $q_{\mathrm v}$. But the translation of $X_0\hat q^{\,t}$, differentiated in a tangent direction of $q$ through the exact arm, errs by $9\times10^3u$, $1.3\times10^6u$, $8.5\times10^7u$ at $\alpha = 10^{-2}, 10^{-4}, 10^{-6}$
  (about $(50$ to $150)\,\alpha^{-1}u$), and the rotation by $\le 3u$ at every $\alpha$. Exact derivatives for $\varpi_t$ and $\cos t\alpha$ do not change this ($7.3\times10^3u$, $1.5\times10^6u$, $8.5\times10^7u$), nor do exact ones
  for $\varkappa$, $\bar m$ and $\varpi_t$ together ($2.5\times10^3u$, $7.9\times10^5u$, $1.4\times10^8u$): each of the three carries $O(\lVert\mathbf t\rVert/\alpha)$ derivative terms that cancel in the assembled
  $q_{\mathrm d}^t$. **So no single-coefficient series, one for $\varpi_t$ included, removes it**, and `jacobians_match_dual_geodesic` at small angle needs a per-stratum tolerance (CO.14(ii)) of order
  $10^2\alpha^{-1}u$ on the translation block; no other remedy is derived here (index, open item).
- (d) *Near $\pi$ nothing cancels; the branch is discontinuous.* With $\cos\alpha \to 0^+$, $\sin\alpha \to 1$ every term is $O(1)$. What is lost is the choice of geodesic: at $\theta(\Delta) = \pi$ the two preimages
  are $d$ and $d'$ with $\varphi' = -\varphi$, $\rho_i' - \rho_i = \varphi\times x_i$ (LG.2(d)); $X_0\mathrm{Exp}(td)$ and $X_0\mathrm{Exp}(td')$ agree at every integer $t$ (both are $X_0\Delta^t$, since $\mathrm{Exp}\,d' = \mathrm{Exp}\,d = \Delta$) and differ in between (by $0.83$ to $1.66$ at
  $t = \tfrac12$ in the samples below), and $\gamma$ just below the cut is one and just above the other. (i) The fast twin's flip predicate ($w < 0$) must be `SO3::log`'s, or the twins disagree at
  $w = \pm0$ (index, open item on `copysign`). (ii) $\Delta$ is a product, so its $w$ has an absolute error of a few $u$: within that of $\pi$ the branch is noise; each result is the exact geodesic of a
  nearby pair but the forward error is $O(1)$, so `se3_geodesic_matches_reference` can only be asked to agree outside it. (iii) GE.6's linearisation holds while $\Delta$'s perturbation stays in
  $\theta < \pi$: $O(c^2)$ before, $O(1)$ after.
- (e) *The reference twin* evaluates `exp_coeffs`, `jr_coeffs` ($\mathrm{Exp}$ at $t\varphi$) and `log_ratio`, `jr_inv_coeff` ($\mathrm{Log}\,\Delta$), each with its generated switch; their errors reach the output
  multiplied by the size of the word they scale (SE.15, SE.16). The differential proptest measures the fast twin against it, not against the truth: the corpus is the arbiter.
- (f) *Endpoints.* The grouped formulas are exact at $t = 0$ ($\varpi_0 = 0$, $\cos 0 = 1$: the identity, so $X(0) = X_0$ bit for bit) but not at $t = 1$: $\varpi_1 = \sin\alpha/\lVert q_{\mathrm v}\rVert$ is not $1$ in
  general, and $\hat q_\Delta^{\,1}$ differs from $\hat q_\Delta$ by up to $2.0u$ (a quaternion component) and $4.9u$ (translation, relative to $\max(1,\lVert\mathbf t\rVert)$) on the strata of (a). `tf_tree_math`'s
  `ScLerp` returns its endpoints by an early return, with a proptest ("exact by construction"); no `helicoid` document says whether `geodesic` must be exact at $t \in \{0,1\}$, for either twin (index, open item).
  Likewise $J_0$ of GE.7(a) is exactly $0$ at $t = 1$ ($1 - t = 0$), where the difference form of §10 leaves a residue of the order of its error (GE.7).

**Checked (GE.13):** (a) `f64` (Python floats, the system `libm`; the arithmetic of the exact arm of `screw_pow`, no series arm) against the definition from the exact binary64 inputs at 60 digits, at quaternion
level (no $3\times3$ matrix, so a tiny angle keeps its digits), mpmath 1.3.0: $\theta \in \{3, 1, 0.31, 0.29, 0.1, 10^{-2}, 10^{-3}, 10^{-5}, 10^{-8}, 10^{-12}, 10^{-20}, 10^{-50}, 10^{-100}, 10^{-140}\}$ and
$\pi - 10^{-k}$, $k \in \{1, 2, 3, 5, 8, 12, 15\}$; translations of norm $1$ and $10^4$ (generic), $10^4$ along the axis, $10^4$ orthogonal, $10^{-3}$ along; 2 poses each, $s \in \{0.125, 0.5, 0.7314, 1 - 10^{-9}\}$
(840 evaluations per draw, three independent draws): maximum quaternion error $2.4u$ ($\lVert\Delta q\rVert$), translation error $5.4u$ ($4.3u$ to $5.4u$ per draw; relative to $\max(1,\lVert\mathbf t\rVert)$), with no trend in $\theta$.
The same at $s = 1.5$, $-0.5$, $3$ (one draw each): quaternion $3.3u$, $1.1u$, $4.4u$, translation $5.0u$, $2.1u$, $11.4u$; the three together, another draw: $4.4u$, $12.6u$. (b) The same procedure with $q = (1, \lVert q_{\mathrm v}\rVert\hat n)$ built directly,
$\lVert q_{\mathrm v}\rVert \in \{10^{-144}, 2\times10^{-145}, 10^{-146}, 10^{-150}, 10^{-153}, 10^{-154}, 10^{-155}, 10^{-160}\}$ (40 evaluations each): translation error $\le 1.7u$ (largest at $10^{-144}$), quaternion error (absolute) $\le 1.5\,u\lVert q_{\mathrm v}\rVert$.
(c) A forward-mode `Dual` over floats with the textbook
rules (`sqrt`, `atan2`, `sin`, quotient $q' = (a' - qb')/b$), $\varpi_t$ differentiated in one component of $q_{\mathrm v}$, 40 random axes and $t \sim U(0.05, 0.95)$ per $\alpha \in \{1, 10^{-1}, \dots, 10^{-7}\}$, against
`mp.diff` at 80 digits: derivative relative error $6.1u$, $4.8\times10^2u$, $4.5\times10^4u$, $4.1\times10^6u$, $5.2\times10^8u$, $4.1\times10^{10}u$, $4.2\times10^{12}u$, $3.4\times10^{14}u$ (times $\alpha^2$: $3.4$ to $6.1$).
The whole exact arm through the same `Dual`, the seed $q\,\mathrm{Exp}(\epsilon\delta)$ with $\mathbf t$ fixed and both rounded to floats (random axis, $\mathbf t \sim \mathcal N(0,2^2)$, $\delta \sim \mathcal N(0,I)$, $t \sim U(0.05, 0.95)$; 8 draws per
$\alpha$), against `mp.diff` at 60 digits of the same arm on the same rounded seeds, relative to the largest derivative entry: translation $2.2\times10^2u$, $7.6\times10^2u$, $8.9\times10^3u$, $5.6\times10^4u$, $1.3\times10^6u$, $7.6\times10^6u$,
$8.5\times10^7u$ at $\alpha = 1, 10^{-1}, \dots, 10^{-6}$; rotation (four components) $\le 2.9u$ at all of them. With exact derivatives substituted for $\varpi_t$ and $\cos t\alpha$: $7.3\times10^3u$, $1.45\times10^6u$, $8.5\times10^7u$
at $\alpha = 10^{-2}, 10^{-4}, 10^{-6}$; for $\varkappa$, $\bar m$, $\varpi_t$ and $\cos t\alpha$ together (5 draws): $2.5\times10^3u$, $7.9\times10^5u$, $1.4\times10^8u$.
The series of $\varpi_t$ through $\alpha^8$ equals the first four correction terms of `tf_tree_math`'s `slerp_weight` (sympy 1.14.0, exact). (d) The two preimages at $\theta = \pi$ ($N = 0,1,2$, 110 digits):
$\mathrm{Exp}\,d' = \mathrm{Exp}\,d$ to $3.9\times10^{-111}$; at integer $t$ ($N = 0,1,2$, three pairs each, 60 digits) $X_0\mathrm{Exp}(td) - X_0\mathrm{Exp}(td')$ is $\le 2.5\times10^{-60}$ for $t \in \{1, 2, -1, 3\}$ and $1.4$ to $3.5$ at $t = \tfrac12$ (largest entry);
midpoints of $X_0\mathrm{Exp}(sd)$, $X_0\mathrm{Exp}(sd')$ differ by $1.17, 1.66, 1.17$ ($N = 0$), $1.25, 0.83, 0.37$ ($N = 1$), $1.42, 1.35, 0.73$ ($N = 2$) at
$s = \tfrac14, \tfrac12, \tfrac34$, and $\gamma$ at $\theta = \pi \mp 10^{-9}$ differs at $s = \tfrac12$ by $1.659$, $0.829$, $1.354$, the same. The reference of `PHASE4.md` §4 by `mp.expm(t*mp.logm(Delta))` (SE(3), $X_0 = I$,
$t = \tfrac12$, 6 poses per $\theta$, 40 digits, mpmath 1.3.0 and 1.4.1 identical) is within $9.2\times10^{-41}$ of the geometric form for $\theta \le 3.0$ and returns a complex, non-principal logarithm (error $1.0$ to $1.9$) from
$\theta = 3.03$ (4 of 6 poses) and for all 6 poses at $\theta \ge 3.05$: `geo:near-pi` cannot use `mp.logm` (index, open item on `so3_log`). Scripts not committed. **Permanent:** planned, corpus `se3_geodesic` strata
(values), `jacobians_match_dual_geodesic` (derivatives); none for the rounding figures of (a), (c).
(f) The exact arm at $s = 1$ and $s = 0$ on the $\theta$ of (a), norms $1$ and $10^4$ in a generic direction, 20 draws per cell (840 poses): $\hat q^{\,1}$ against $\hat q$ at most $2.0u$ per quaternion component and $4.9u$ in the
translation (relative to $\max(1,\lVert\mathbf t\rVert)$); at $s = 0$ exactly the identity ($0$ deviation). Scripts not committed.

## 9. SO(3): slerp

**Proposition GE.14.** For unit quaternions $q_0, q_1$ with $q_0\cdot q_1 \ge 0$ (that is $w(q_0^*q_1) \ge 0$), $q_0^*q_1 = (w, v)$, $\alpha = \mathrm{atan2}(\lVert v\rVert, w)$ and $t \in \mathbb R$, the rotation part
of GE.12 (no dual part) is
$$
q_0\,(\cos t\alpha,\ \varpi_t\,v) = \frac{\sin((1-t)\alpha)\,q_0 + \sin(t\alpha)\,q_1}{\sin\alpha} = q_0\,\mathrm{Exp}\big(t\,\mathrm{Log}(q_0^*q_1)\big):
$$
`SO3`'s reference geodesic is the shortest-arc slerp with $\mathrm{Log}$'s sign rule ($q_1 \to -q_1$ when $q_0\cdot q_1 < 0$); on $[0,1]$ it is the short arc of the great circle through $q_0, q_1$; at
$q_0\cdot q_1 = 0$ ($\theta = \pi$) there are two arcs, the two preimages of GE.13(d). Its small-angle behaviour is that of the rotation part in GE.13: no cancellation in the value (GE.13(a)) and a derivative through `Dual` within $3u$ (GE.13(c)).

*Proof.* $q_1 = q_0(\cos\alpha, \sin\alpha\,\hat n)$ and $\sin((1-t)\alpha) + \sin(t\alpha)\cos\alpha = \sin\alpha\cos t\alpha$ give
$[\sin((1-t)\alpha) + \sin(t\alpha)(\cos\alpha + \sin\alpha\,\hat n)]/\sin\alpha = \cos t\alpha + \sin t\alpha\,\hat n$. $\square$

`tf_tree_math::slerp` computes the same function from the two quaternions with a chord-based $\theta^2$ and a normalised-LERP fallback below $10^{-6}$; `helicoid` does not reproduce the fallback
(`PHASE4.md` §0).

**Which of the three spellings ships, and why it is the right-hand one:** `0050` scored six of them over the committed `so3_geodesic` corpus and found the whole
domination gap in the **denominator** -- $\sin\alpha$ recomputed from $\alpha$ against $\lVert v\rVert$, which are the same number for a unit quaternion and not the same
floating-point value. Only the recomputed form makes both endpoints exact, since each weight is then one number divided by itself at its own endpoint. The
grouped middle expression is **bit-identical** to $q_0\mathrm{Exp}(t\mathrm{Log}(q_0^*q_1))$ where both are best ($k_t\,t\,r\,v$ and $\varpi_t v$ are one product
re-associated), so GE.12's rotation part is not an alternative to it. `NUMERICS.md` §10 carries the normative spelling. The sign of the output quaternion does not matter: the rotation is the same.

**Checked:** as GE.5 (slerp against $R_0\mathrm{Exp}(s\,\mathrm{Log}(R_0^\top R_1))$ by matrices, $4.0\times10^{-111}$); the identity $\sin((1-t)\alpha) + \sin(t\alpha)\cos\alpha = \sin\alpha\cos t\alpha$ exactly in sympy.
Script not committed. **Permanent:** planned, corpus `so3_geodesic`, with `tf_tree_math`'s `slerp` as oracle (`PHASE4.md` §4).
