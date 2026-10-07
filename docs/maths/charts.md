# Charts: retractions, the SE(3) charts and S²

> Non-normative companion to [`NUMERICS.md`](../NUMERICS.md) §1, §2.3, §8 (S²), §12 (`S2Chart::local`),
> §14 (its twin) and [`PHASE5.md`](../PHASE5.md) §1–§2 ([`0012`](../decisions/0012-a-retraction-is-a-chart.md)).
> **`NUMERICS.md` wins on any conflict; a disagreement is an open item in the [maths index](./index.md)**,
> which holds the notation and the `Checked:` convention. $\mathrm{Exp}$, $J_r$, $J_l$ and the sides table:
> [`lie-groups.md`](./lie-groups.md); SO(3): [`so3.md`](./so3.md); the block $Q$: [`se3.md`](./se3.md); the
> kernel $r$: [`coefficients.md`](./coefficients.md) CO.16.

Throughout $X = (R, t)$ and $Y = (R_Y, t_Y)$ are in $\mathrm{SE}(3)$ (the group charts hold for every
$\mathrm{SE}_N(3)$), $\delta = [\varphi;\rho]$, $\theta = \lVert\varphi\rVert$, $\varphi_Y = \mathrm{Log}(R^\top R_Y)$.
$\mathrm{ret}_X$, $\mathrm{loc}_X$ are the `retract` and `local` of the chart frozen at $X$, and $\mathrm{rj}_X$,
$\mathrm{lj}_X$ its `retract_jacobian` and `local_jacobian`. Group Jacobians are right-sided unless a side is
named. In §5, $n, m$ are unit vectors of $S^2$ as in `NUMERICS.md` §8 (not a dimension) and $\mathrm{Exp}(v)$,
$v \in \mathbb R^3$, is the SO(3) exponential.

## Results

| Label | Result | `NUMERICS.md`; `PHASE5.md` (P5) | Implemented by |
|---|---|---|---|
| CH.1–CH.2 | a chart is a retraction with a local inverse, frozen at its base; $\mathrm{rj}$, $\mathrm{lj}$ read the perturbed point in the chart frozen at that point; $\mathrm{lj}(\mathrm{ret}\,\delta) = \mathrm{rj}(\delta)^{-1}$ | P5 §1.1 | `Chart` |
| CH.3 | `RightChart`: $\mathrm{rj} = J_r(\delta)$, $\mathrm{lj} = J_r^{-1}(\mathrm{Log}\,X^{-1}Y)$; `LeftChart`: $J_l$, $J_l^{-1}$; left $=$ right after the linear map $\mathrm{Ad}_X$ | §2.3; P5 §1.1 | `RightChart<G>`, `LeftChart<G>` |
| CH.4 | the three SE(3) charts are two retractions in three frames: transition maps, on $U_0 = \{\theta < \pi\}$, $(\varphi, J_l(\varphi)\rho)$ and $(\varphi, R\rho)$ | P5 §1.3 | `Screw`, `Decoupled`, `WorldTranslation` |
| CH.5 | $\mathrm{rj}$, $\mathrm{lj}$ of `Decoupled` ($\mathrm{diag}(J_r(\varphi), \mathrm{Exp}(-\varphi))$) and `WorldTranslation`, derived; why they are not dual matrices; which pairs agree to first order | P5 §1.3–§1.4 | `ProductJac`, `Chart` |
| CH.6 | change of chart: $\mathrm{rj}'(\Phi\delta)\,\mathrm D\Phi = \Psi\,\mathrm{rj}(\delta)$; it contains SE.9(b), (c) | P5 §1.4 | (checks) |
| CH.7 | what changes in an LM iterate: same step for `Screw`, `Decoupled`; the trial point differs by $\ell(\theta)\lVert\rho_\perp\rVert \le \tfrac\theta2\lVert\rho\rVert$ | P5 §1.3, `0012` | (the choice of chart) |
| CH.8 | the Householder basis: $B^\top B = I$, $B^\top n = 0$, $b_1 \times b_2 = \varsigma n$, $\nu^\top\nu = 2(1 + \lvert n_z\rvert)$, closed form | §8 | `S2Chart::at` |
| CH.9 | $n \oplus \delta = \cos\theta\,n - \varsigma\frac{\sin\theta}\theta BK\delta$; $\mathrm{loc}_n(m) = \varsigma\frac\alpha s KB^\top m$; $\frac\alpha s = \frac12 r(s^2, w) = r(s^2, 1 + w)$; domains | §8, §4, §12, §14; P5 §2 | `S2Chart::{retract, local}` |
| CH.10 | $\mathrm{rj} = Z^\top A(\delta)$, $\mathrm{lj} = A(\delta)^{-1}Z$, $A^{-1} = \frac\alpha s(I - b\,\delta\delta^\top)$, $\lVert\mathrm{lj}\rVert_2 = \alpha/\sin\alpha$, $\det Z = \varsigma\varsigma_m$; $\mathrm{rj}$ singular at $\theta = \pi$ | P5 §1.4, §2 | `S2Chart::{retract_jacobian, local_jacobian}` |
| CH.11 | the basis jumps by a reflection at $n_z = 0$; no continuous basis exists; why the chart is frozen; $m \ne -n$; `f64` figures | §8, §12; P5 §2 | `S2Chart::at` |
| CH.12 | sign audit: the usual wrong variant of each entry and its exact gap | — | — |

## 1. What a chart is

**Definition CH.1.** Let $\mathcal M$ be a smooth $d$-manifold and $X \in \mathcal M$. A *chart at $X$* is a pair of smooth maps
$\mathrm{ret}_X : U \to \mathcal M$ ($0 \in U \subseteq \mathbb R^d$ open) and $\mathrm{loc}_X : V \to \mathbb R^d$
($V \subseteq \mathcal M$ open), which depend on $X$ alone (`Chart::at` freezes them), with

- (Ret) $\mathrm{ret}_X(0) = X$;
- (Inv) $\mathrm{loc}_X \circ \mathrm{ret}_X = \mathrm{id}$ on a neighbourhood $U_0$ of $0$, hence
  $\mathrm{ret}_X \circ \mathrm{loc}_X = \mathrm{id}$ on $\mathrm{ret}_X(U_0)$.

Two charts $C, C'$ at $X$ *agree to first order* if $\mathrm D\,\mathrm{ret}_X(0) = \mathrm D\,\mathrm{ret}'_X(0)$. Their
*transition* is $\Phi = \mathrm{loc}'_X \circ \mathrm{ret}_X$, a local diffeomorphism with $\Phi(0) = 0$; they agree to first
order iff $\mathrm D\Phi(0) = I$. With $Y = \mathrm{ret}_X(\delta)$, the Jacobians of `PHASE5.md` §1.1 are

$$
\mathrm{rj}_X(\delta) = \partial_\eta\big[\mathrm{loc}_Y(\mathrm{ret}_X(\delta + \eta))\big]_{\eta = 0},\qquad
\mathrm{lj}_X(Y) = \partial_\eta\big[\mathrm{loc}_X(\mathrm{ret}_Y(\eta))\big]_{\eta = 0}:
$$

the perturbed point is read in the chart frozen at *that point* ("in the frame of the retracted point"), never in the
chart of $X$. This is what carries a tangent perturbation, a gradient ($\mathrm{rj}^{-\top}g$) or a covariance
($\mathrm{rj}\,\Sigma\,\mathrm{rj}^\top$) from the coordinates at $X$ to those at $Y$ after a step.

**Proposition CH.2.** For $\delta \in U_0$ and $Y = \mathrm{ret}_X(\delta)$: (a) $\mathrm{rj}_X(0) = I$; (b) $\mathrm{rj}_X(\delta)$
is invertible and $\mathrm{lj}_X(Y) = \mathrm{rj}_X(\delta)^{-1}$.

*Proof.* (a) $\mathrm{rj}_X(0) = \mathrm D(\mathrm{loc}_X \circ \mathrm{ret}_X)(0) = I$ by (Inv). (b) $\mathrm{ret}_X$ is a diffeomorphism
near $\delta$ and $\mathrm{loc}_Y$ one near $Y$, so $\mathrm{rj}$ is invertible. By (Inv) at $Y$,
$\mathrm{ret}_X(\delta + \eta) = \mathrm{ret}_Y(\zeta)$ with $\zeta = \mathrm{rj}\,\eta + O(\lVert\eta\rVert^2)$, so
$\eta = \mathrm{rj}^{-1}\zeta + O(\lVert\zeta\rVert^2)$, and $\mathrm{loc}_X(\mathrm{ret}_Y(\zeta)) = \delta + \eta$. $\square$

**Checked:** with CH.3 (the figure for $\mathrm{lj}\,\mathrm{rj} = I$ is there). **Permanent:** none specified for CH.2.

## 2. Group charts

**Proposition CH.3 (`RightChart`, `LeftChart`; every $\mathrm{SE}_N(3)$, $N \ge 0$).** Let
$\mathrm{ret}^R_X(\delta) = X\,\mathrm{Exp}\,\delta$, $\mathrm{loc}^R_X(Y) = \mathrm{Log}(X^{-1}Y)$ and
$\mathrm{ret}^L_X(\delta) = \mathrm{Exp}(\delta)X$, $\mathrm{loc}^L_X(Y) = \mathrm{Log}(YX^{-1})$.

- (a) Both are charts, (Inv) holding on $U_0 = \{\theta < \pi\}$ (LG.2(c)); $\mathrm{ret}\circ\mathrm{loc} = \mathrm{id}$ for
  every $Y$ (at $\theta(X^{-1}Y) = \pi$ either logarithm maps back, LG.2(d)); $\mathrm{loc}$ is discontinuous at
  $\theta(X^{-1}Y) = \pi$.
- (b) $\mathrm{rj}^R_X(\delta) = J_r(\delta)$ and $\mathrm{lj}^R_X(Y) = J_r^{-1}(\tau)$, $\tau = \mathrm{Log}(X^{-1}Y)$,
  $\theta(\tau) < \pi$; $\mathrm{rj}^L_X(\delta) = J_l(\delta)$ and $\mathrm{lj}^L_X(Y) = J_l^{-1}(\tau)$,
  $\tau = \mathrm{Log}(YX^{-1})$.
- (c) $\mathrm{ret}^L_X(\delta) = \mathrm{ret}^R_X(\mathrm{Ad}_X^{-1}\delta)$: the transition is $\Phi^{R \to L} = \mathrm{Ad}_X$, exact and linear,
  and $\Psi_Y = \mathrm{Ad}_Y$. The two charts agree to first order only at $X = I$.

*Proof.* (a) LG.2(c), LG.14's domain. (b) Right: $\mathrm{loc}_Y(X\,\mathrm{Exp}(\delta + \eta)) = \mathrm{Log}(\mathrm{Exp}(\delta)^{-1}\mathrm{Exp}(\delta + \eta))
= \mathrm{Log}\,\mathrm{Exp}(J_r(\delta)\eta + O) = J_r(\delta)\eta + O$ by LG.8; and
$\mathrm{loc}_X(Y\,\mathrm{Exp}\,\eta) = \mathrm{Log}(\mathrm{Exp}(\tau)\mathrm{Exp}\,\eta) = \tau + J_r^{-1}(\tau)\eta + O$ by LG.13(d). These are the
rows "$X \oplus_R \tau$ with respect to $\tau$" and "$Y \ominus_R X$ with respect to $Y$" of `NUMERICS.md` §2.3 (LG.14). Left: the same
with $Y = \mathrm{Exp}(\delta)X$, $\mathrm{loc}_Y(\mathrm{Exp}(\delta + \eta)X) = \mathrm{Log}(\mathrm{Exp}(\delta + \eta)\mathrm{Exp}(-\delta))$ and
LG.8, LG.13(d) on the left. (c) LG.4(b): $\mathrm{Exp}(\delta)X = X\,\mathrm{Exp}(\mathrm{Ad}_X^{-1}\delta)$. $\square$

`Screw` is `RightChart<SE3>`; its $\mathrm{rj}$ is the dual matrix $J_r(\varphi) + \epsilon\,Q(-\rho, -\varphi)$ (SE.9(a)).

**Checked:** mpmath 1.3.0, 110 digits, SE(3). $\mathrm{Exp}$ by `mp.expm`; $\mathrm{Log}$ by the geometric SO(3) inverse (SO.5) and a
$\Gamma_1$ solve from its series; $J_r$, $J_l$ summed from Definition LG.7; the four charts of CH.3–CH.4 coded from their definitions.
13 cases per chart: $\theta(\varphi) \in \{0, 10^{-10}, 0.3, 1.1, 2.5, 3.0, \pi - 10^{-3}\}$, two random $(X, \delta)$ each
($\theta(R) \sim U[0.05, 3]$, $t, \rho \sim \mathcal N(0, 2^2)$; one at $\theta = 0$; one of the $\theta = 2.5$ pair with $\rho \sim \mathcal N(0, (10^6)^2)$).
Round trips, both ways: $\le 9.1\times10^{-109}$. $\mathrm{rj}$, $\mathrm{lj}$ by central differences ($h = 10^{-30}$) of the composite
operation itself, in the chart frozen at the point (Definition CH.1), against $J_r(\delta)$ and $J_r^{-1}(\tau)$
(`Screw`; $J_l$ for `Left`): $\le 1.2\times10^{-61}$ (the differencing; relative to the largest entry, absolute below $1$).
$\mathrm{lj}\,\mathrm{rj} = I$ (CH.2(b)): $\le 2.7\times10^{-61}$, and $2.9\times10^{-56}$ for the $\rho \sim 10^6$ case. $\Phi^{R \to L} = \mathrm{Ad}_X$
and $\mathrm D\Phi(0) = \mathrm{Ad}_X$ (12 cases as CH.4): $5.7\times10^{-111}$, $3.3\times10^{-81}$. LG.14 covers $N = 0..3$ for the rows used.
Script not committed. **Permanent:** planned, proptest `chart_jacobians_match_dual_*` (`PHASE5.md` §1.4) for $\mathrm{rj}$, $\mathrm{lj}$ of
each chart; none for CH.2(b) or the transition.

## 3. The three SE(3) charts

By SE.3, $\mathrm{Exp}\,\delta = (\mathrm{Exp}\,\varphi,\ J_l(\varphi)\rho)$ and $X^{-1}Y = (R^\top R_Y,\ R^\top(t_Y - t))$. So
(`PHASE5.md` §1.3)

| Chart | $\mathrm{ret}_X(\delta)$ | $\mathrm{loc}_X(Y)$ |
|---|---|---|
| `Screw` | $(R\,\mathrm{Exp}\,\varphi,\ t + R\,J_l(\varphi)\rho)$ | $(\varphi_Y,\ J_l^{-1}(\varphi_Y)R^\top(t_Y - t))$ |
| `Decoupled` | $(R\,\mathrm{Exp}\,\varphi,\ t + R\rho)$ | $(\varphi_Y,\ R^\top(t_Y - t))$ |
| `WorldTranslation` | $(R\,\mathrm{Exp}\,\varphi,\ t + \rho)$ | $(\varphi_Y,\ t_Y - t)$ |

**Proposition CH.4 (charts and transitions).**

- (a) Each is a chart with (Inv) on $U_0 = \{\theta < \pi\}$, $\mathrm{loc}$ defined on $V = \{\theta(R^\top R_Y) < \pi\}$ (the same set
  for all three, and for `LeftChart`), $\mathrm{ret}\circ\mathrm{loc} = \mathrm{id}$ on every $Y$.
- (b) For $\delta \in U_0$ the transitions at one base $X$ are
  $\Phi^{\mathrm{Scr} \to \mathrm{Dec}}(\varphi, \rho) = (\varphi,\ J_l(\varphi)\rho)$,
  $\Phi^{\mathrm{Dec} \to \mathrm{WT}}(\varphi, \rho) = (\varphi,\ R\rho)$ (linear), hence
  $\Phi^{\mathrm{Scr} \to \mathrm{WT}}(\varphi, \rho) = (\varphi,\ RJ_l(\varphi)\rho)$. For $\pi < \theta < 2\pi$ the rotation part of each is
  $\mathrm{Log}\,\mathrm{Exp}\,\varphi = (1 - 2\pi/\theta)\varphi$ (LG.2(c)), so these maps and the $\mathrm D\Phi$ of CH.6 are those of $U_0$ only.

So there are two retractions: the SE(3) one (`Screw`) and the decoupled one; `WorldTranslation` is `Decoupled` with the translation
tangent expressed in the world frame instead of the body frame, and nothing else.

*Proof.* (a) $\mathrm{Log}\,\mathrm{Exp}\,\varphi = \varphi$ for $\theta < \pi$ (LG.2(c)) and $R^\top R\rho = \rho$; conversely
$(R\,\mathrm{Exp}\,\varphi_Y,\ t + R\,R^\top(t_Y - t)) = (R_Y, t_Y)$. (b) Substitute one $\mathrm{ret}$ into the other $\mathrm{loc}$, whose rotation part is $\mathrm{Log}\,\mathrm{Exp}\,\varphi$. $\square$

**Checked:** as CH.3 (13 cases per chart, all four charts). (b): 12 $(X, \delta)$, $\theta \in \{0.3, 1.1, 2.5, 3.0\}$ three each,
$\mathrm{loc}'(\mathrm{ret}\,\delta)$ formed from the charts themselves against the stated maps: $\le 5.4\times10^{-111}$ ($\mathrm{Scr}\to\mathrm{Dec}$),
$5.1\times10^{-111}$ ($\mathrm{Dec}\to\mathrm{WT}$), $5.3\times10^{-111}$ ($\mathrm{Scr}\to\mathrm{WT}$). Beyond $U_0$ (mpmath 1.3.0, 60 digits, $\mathrm{Scr}\to\mathrm{Dec}$, $\theta \in \{3.3, 4, 5, 6\}$, one
$(X, \delta)$ each, $\mathrm{Exp}$ by `mp.expm`, $\mathrm{Log}$ geometric): the rotation part against $(1 - 2\pi/\theta)\varphi$ $\le 6.0\times10^{-60}$ (against $\varphi$: $2\pi$), the translation part against $J_l(\varphi)\rho$ $\le 2.8\times10^{-61}$. Script not committed. **Permanent:** none
specified for the transitions; the round trips of (a) are the `exp_log_roundtrip_*` family applied to the charts (planned).

**Proposition CH.5 (chart Jacobians of `Decoupled`, `WorldTranslation`).** With $\varphi_Y = \mathrm{Log}(R^\top R_Y)$
($\theta(\varphi_Y) < \pi$ for $\mathrm{lj}$; no domain for $\mathrm{rj}$):

$$
\mathrm{rj}^{\mathrm{Dec}}_X(\delta) = \begin{bmatrix}J_r(\varphi) & 0\\ 0 & \mathrm{Exp}(-\varphi)\end{bmatrix},\quad
\mathrm{lj}^{\mathrm{Dec}}_X(Y) = \begin{bmatrix}J_r^{-1}(\varphi_Y) & 0\\ 0 & \mathrm{Exp}(\varphi_Y)\end{bmatrix},\quad
\mathrm{rj}^{\mathrm{WT}}_X = \begin{bmatrix}J_r(\varphi) & 0\\ 0 & I\end{bmatrix},\quad
\mathrm{lj}^{\mathrm{WT}}_X = \begin{bmatrix}J_r^{-1}(\varphi_Y) & 0\\ 0 & I\end{bmatrix}.
$$

- (a) $\mathrm{rj}$ depends on $\delta$ only through $\varphi$, and $\mathrm{lj}$ on $X$ and $Y$ only through $\varphi_Y$; neither depends on $\rho$ or $t$.
- (b) They are not dual matrices: `SEn3Jac` (`NUMERICS.md` §2.2) has equal diagonal blocks, and $J_r(\varphi) \ne \mathrm{Exp}(-\varphi)$
  and $J_r(\varphi) \ne I$ for every $\theta > 0$. The block-diagonal
  matrices are closed under product and inverse, which is what `ProductJac<Mat3, Mat3>` needs.
- (c) $\mathrm D\,\mathrm{ret}_X(0)\delta = X\delta^\wedge$ for `Screw` and `Decoupled`, and $X\big(\mathrm{diag}(I, R^\top)\delta\big)^\wedge$ for
  `WorldTranslation`. So `Screw` and `Decoupled` agree to first order at every $X$, and `WorldTranslation` agrees with them only if $R = I$:
  $\mathrm D\Phi^{\mathrm{Dec} \to \mathrm{WT}}(0) = \mathrm{diag}(I, R)$ (`PHASE5.md` §1.3 says all three agree; index, open items).

*Proof.* (`Decoupled`) $Y = \mathrm{ret}_X(\delta)$ has $R_Y = R\,\mathrm{Exp}\,\varphi$ and $\mathrm{ret}_X(\delta + \eta) = (R\,\mathrm{Exp}(\varphi + \eta_\varphi),\ t + R(\rho + \eta_\rho))$. In the chart at $Y$ its rotation
part is $\mathrm{Log}(\mathrm{Exp}(-\varphi)\mathrm{Exp}(\varphi + \eta_\varphi)) = J_r(\varphi)\eta_\varphi + O(\lVert\eta\rVert^2)$ (LG.8, $N = 0$) and its translation part is
$R_Y^\top R\,\eta_\rho = \mathrm{Exp}(-\varphi)\eta_\rho$, exactly. (`WorldTranslation`) $t_Y = t + \rho$ and the translation part is $\eta_\rho$. For $\mathrm{lj}$: $\mathrm{ret}_Y(\eta) = (R_Y\mathrm{Exp}\,\eta_\varphi,\ t_Y + R_Y\eta_\rho)$
gives, in $\mathrm{loc}_X$, the rotation part $\mathrm{Log}(\mathrm{Exp}(\varphi_Y)\mathrm{Exp}\,\eta_\varphi) = \varphi_Y + J_r^{-1}(\varphi_Y)\eta_\varphi + O$ (LG.13(d)) and the translation part
$R^\top(t_Y - t) + R^\top R_Y\eta_\rho$, whose derivative is $\mathrm{Exp}(\varphi_Y)$ ($I$ for `WorldTranslation`); these are the inverses of the $\mathrm{rj}$ blocks (CH.2(b)). (b) $J_r = I - aW + bW^2$ and $\mathrm{Exp}(-\varphi) = I - \frac{\sin\theta}\theta W + aW^2$, $a = \frac{1 - \cos\theta}{\theta^2}$, $b = \frac{\theta - \sin\theta}{\theta^3}$; $W$ and $W^2$ are independent for $\theta > 0$ (SO.3), so equality needs both
$a = \frac{\sin\theta}\theta$ and $b = a$. The first holds at $\theta \in 2\pi\mathbb Z_{>0}$ and where $\tan\frac\theta2 = \theta$ ($2.3311$, the first). The second, $b - a = \frac{\theta\cos\theta - \sin\theta}{\theta^3} = 0$, only where $\tan\theta = \theta$ ($4.4934$, the first);
it fails on $(0, \pi]$, where $b - a < 0$ (the numerator vanishes at $0$ and has derivative $-\theta\sin\theta \le 0$). The two sets are disjoint: $b - a = \theta^{-2}$ at $2\pi k$, and $\tan\frac\theta2 = \theta$ with $\tan\theta = \frac{2\tan(\theta/2)}{1 - \tan^2(\theta/2)} = \theta$ would need $\theta^2 = -1$. At $\theta = 2.3311$, where the $W$ terms agree, $b - a = -0.18402$.
$J_r(\varphi) = I$ would need $a = b = 0$, and $a = 0$ only at $\theta \in 2\pi\mathbb Z_{>0}$, where $b \ne 0$.
(c) Differentiate the table of §3 at $0$. $\square$

$\mathrm{rj}$ of `Screw` needs $Q$ (`q_coeffs`: $b, d, e$) and $\mathrm{lj}$ its dual inverse; those of `Decoupled` and `WorldTranslation` need only $J_r(\varphi)$ (`jr_coeffs`),
$J_r^{-1}(\varphi)$ (`jr_inv_coeff`) and the rotation matrix already formed by $\mathrm{ret}$, and their $\mathrm{ret}$, $\mathrm{loc}$ need `exp_coeffs`, `log_ratio` alone (no $J_l(\varphi)$, no $J_l^{-1}$). $e$ is
the catalogue's worst cancellation (`NUMERICS.md` §4).

**Checked:** as CH.3: `Decoupled` and `WorldTranslation` $\mathrm{rj}$ $\le 1.9\times10^{-62}$, $\mathrm{lj}$ $\le 4.5\times10^{-62}$, $\mathrm{lj}\,\mathrm{rj} = I$ $\le 7.9\times10^{-62}$.
(c): $\mathrm D\Phi(0)$ by central differences at 12 bases, against $I$ (`Screw` to `Decoupled`) $2.5\times10^{-81}$ and against $\mathrm{diag}(I, R)$ (`Decoupled` and `Screw` to `WorldTranslation`) $3.3\times10^{-81}$;
for $R = \mathrm{Exp}(0, 0, 1.2)$, $t = 0$ and $\delta = [0; 10^{-8}, 0, 0]$ the translations of `WorldTranslation` and of the other two differ by $1.1292849\times10^{-8} = 2\sin(0.6)\lVert\rho\rVert$. (b) (60 digits): the roots by
`mp.findroot`, $2.33112237$ ($\tan\frac\theta2 = \theta$) and $4.49340946$ ($\tan\theta = \theta$); $b - a = -0.18402$ at the first and $\max(b - a) = -0.1014$ on $(0, \pi]$ (grid $10^{-3}$); $\lVert J_r - \mathrm{Exp}(-\varphi)\rVert_1$, random axis, $\theta = 0.01, 0.02, \dots, 6.99$: minimum $7.1\times10^{-3}$ (at $\theta = 0.01$). Script not committed. **Permanent:** planned, `chart_jacobians_match_dual_*` (`PHASE5.md` §1.4); none for (b), (c).

## 4. Change of chart, and what an LM iterate sees

**Proposition CH.6 (change of chart).** Let $C, C'$ be charts at $X$ with transition $\Phi$, $\delta \in U_0$ (so $\Phi = \mathrm{loc}'_X \circ \mathrm{ret}_X$ is the local diffeomorphism of CH.1; the rows below
are those of $U_0 = \{\theta < \pi\}$, CH.4(b)), $Y = \mathrm{ret}_X(\delta)$, and
$\Psi_Y = \mathrm D\Phi_Y(0)$ the first-order transition of the same two chart types at $Y$. Then

$$
\mathrm{rj}'_X(\Phi(\delta))\ \mathrm D\Phi(\delta) = \Psi_Y\ \mathrm{rj}_X(\delta),\qquad
\mathrm{lj}'_X(Y)\ \Psi_Y = \mathrm D\Phi(\delta)\ \mathrm{lj}_X(Y).
$$

| Pair | $\mathrm D\Phi(\delta)$ | $\Psi_Y$ | What the first identity says |
|---|---|---|---|
| `Screw` $\to$ `Decoupled` | $\begin{bmatrix}I & 0\\ \partial_\varphi(J_l(\varphi)\rho) & J_l(\varphi)\end{bmatrix}$ | $I$ | $\mathrm{Exp}(-\varphi)J_l(\varphi) = J_r(\varphi)$ and $\mathrm{Exp}(-\varphi)\,\partial_\varphi(J_l(\varphi)\rho) = Q(-\rho, -\varphi)$: SE.9(b) and (c) together |
| `Decoupled` $\to$ `WorldTranslation` | $\mathrm{diag}(I, R)$ | $\mathrm{diag}(I, R_Y)$ | $R_Y\mathrm{Exp}(-\varphi)R^\top = I$ |
| `RightChart` $\to$ `LeftChart` | $\mathrm{Ad}_X$ | $\mathrm{Ad}_Y$ | $J_l(\mathrm{Ad}_X\delta) = \mathrm{Ad}_XJ_l(\delta)\mathrm{Ad}_X^{-1}$ (LG.9(b)) |

A covariance at $X$ changes chart by $\mathrm D\Phi(0)\,\Sigma\,\mathrm D\Phi(0)^\top$: $\mathrm{diag}(I, R)$ turns the body-frame translation block into the world-frame
one, and $\mathrm{Ad}_X$ is `Gaussian::to_left` (`PHASE5.md` §5).

*Proof.* $\mathrm{ret}'_X(\Phi(\delta + \eta)) = \mathrm{ret}_X(\delta + \eta)$ ($\Phi$ inverts $\mathrm{ret}'_X$ after $\mathrm{ret}_X$). Apply $\mathrm{loc}'_Y$ and differentiate in $\eta$ at $0$. The left side gives
$\mathrm{rj}'_X(\Phi(\delta))\,\mathrm D\Phi(\delta)$, the point $Y$ being the same. By (Inv) at $Y$, $\mathrm{ret}_X(\delta + \eta) = \mathrm{ret}_Y(\mathrm{rj}_X(\delta)\eta + O)$, so the right side is
$\Phi_Y(\mathrm{rj}_X(\delta)\eta + O)$, of derivative $\Psi_Y\mathrm{rj}_X(\delta)$. The second identity is the first inverted (CH.2(b)). The rows: CH.3–CH.5, with $R_Y = R\,\mathrm{Exp}\,\varphi$ in the second and
$\mathrm{ad}_{\mathrm{Ad}_X\delta} = \mathrm{Ad}_X\mathrm{ad}_\delta\mathrm{Ad}_X^{-1}$ in the third. $\square$

**Checked:** as CH.4, 12 cases per pair; $\mathrm D\Phi(\delta)$ and $\Psi_Y$ by central differences of the transitions built from the charts, $\mathrm{rj}$, $\mathrm{lj}$ of both
charts by differencing (CH.3); both identities, largest deviation: $1.3\times10^{-61}$, $1.9\times10^{-61}$ ($\mathrm{Scr}\to\mathrm{Dec}$), $5.3\times10^{-81}$, $1.1\times10^{-80}$ ($\mathrm{Dec}\to\mathrm{WT}$,
linear), $1.2\times10^{-61}$, $1.8\times10^{-61}$ ($\mathrm{Scr}\to\mathrm{WT}$), $3.2\times10^{-62}$, $8.7\times10^{-62}$ (right to left). $\mathrm{Exp}(-\varphi)\,\partial_\varphi(J_l(\varphi)\rho)$
(central differences of the SO(3) $J_l$ series, $h = 10^{-30}$) against $Q(-\rho, -\varphi)$ of `NUMERICS.md` §5.3 ($b, d, e$ from their definitions) and against the lower-left block of the dense
$J_r(\tau)$ series: $8.3\times10^{-62}$. Script not committed. **Permanent:** planned, `jl_is_ad_jr_*` for the $\mathrm{Ad}$ identity; the rest none specified.

**Proposition CH.7 (an LM iterate).** Let $f(X) = \tfrac12\lVert r(X)\rVert^2$, $J_C = \mathrm D(r \circ \mathrm{ret}_X)(0)$ in chart $C$, and a step
$(J_C^\top J_C + \lambda\Lambda)\delta_C = -J_C^\top r$ accepted as $X^+ = \mathrm{ret}_X(\delta_C)$, with $\lambda \ge 0$ the damping and $\Lambda$ a symmetric positive-definite scaling
($I$, or $\mathrm{diag}(J_C^\top J_C)$).

- (a) $J_C$ transforms with the first-order transition: $J_{C'} = J_C\,\mathrm D\Phi(0)^{-1}$. So $J_{\mathrm{Scr}} = J_{\mathrm{Dec}}$ (same normal equations at the base), and $J_{\mathrm{WT}} = J_{\mathrm{Dec}}\,\mathrm{diag}(I, R^\top)$.
- (b) With $\lambda = 0$, $\delta_{C'} = \mathrm D\Phi(0)\delta_C$ for every pair. With $\Lambda = I$ this holds when $\mathrm D\Phi(0)$ is orthogonal: `Decoupled` and `WorldTranslation` take the same step up to the frame
  ($\mathrm{diag}(I, R)$), so reach the same $X^+$. Not with $\Lambda = \mathrm{diag}(J^\top J)$, nor for a per-coordinate bound or prior, nor between right and left charts unless $t = 0$.
- (c) `Screw` and `Decoupled` take the *same* $\delta$ (for any $\Lambda$ built from $J^\top J$) and land on the same rotation and on translations that differ by $R(J_l(\varphi) - I)\rho$, of norm
  $\ell(\theta)\lVert\rho_\perp\rVert$, $\rho_\perp$ the part of $\rho$ orthogonal to $\varphi$ and
  $\ell(\theta)^2 = 1 - 2\frac{\sin\theta}\theta + \frac{2(1 - \cos\theta)}{\theta^2}$, so $\ell(\theta) = \frac\theta2\big(1 - \frac{\theta^2}{36} + O(\theta^4)\big) \le \frac\theta2$; to first order the gap is $\frac12R(\varphi\times\rho)$.
- (d) So changing the chart changes neither the critical points nor the minimum, nor the step of `Screw` and `Decoupled`. It changes the *trial point*, by at most $\frac\theta2\lVert\rho\rVert$ in the translation ($4.9986$ cm for
  a step of $0.1$ rad and $1$ m, $0.5$ mm for $0.01$ rad and $0.1$ m; the rotation is identical), hence its cost, the gain ratio, the accept or reject and the update of $\lambda$, hence every later iterate. Between `Decoupled` and
  `WorldTranslation` there is no such gap (CH.4(b) is linear). By (a) the Gauss–Newton matrix $J^\top J$ at the base (weights absorbed in $r$) is the same for `Screw` and `Decoupled`. `0012`'s Context reports a `JᵀWJ` conditioning improvement for decoupled retractions, without saying against which
  chart; this page does not reproduce it.

*Proof.* (a) $\mathrm{ret}'_X = \mathrm{ret}_X \circ \Phi^{-1}$ near $0$. (b) Substitute (a) in the normal equations; an orthogonal $D$ leaves $\lambda I$ alone, since $D^{-\top}(J^\top J + \lambda I)D^{-1}$ needs $D^\top D = I$.
(c) $\mathrm{ret}^{\mathrm{Scr}}_X(\delta) - \mathrm{ret}^{\mathrm{Dec}}_X(\delta) = (0,\ R(J_l(\varphi) - I)\rho)$ (CH.4). $J_l - I = aW + bW^2$ is zero on $\varphi$ and, on the plane $\perp \varphi$, multiplication by
$g(i\theta) - 1$, $g(z) = (e^z - 1)/z$ (LG.7), a similarity of ratio $\ell(\theta) = \lvert g(i\theta) - 1\rvert = \lvert\frac{2\sin(\theta/2)}\theta e^{i\theta/2} - 1\rvert$, whose square is the stated form; also $\ell(\theta) \le \int_0^1\lvert e^{i\theta\lambda} - 1\rvert d\lambda \le \frac\theta2$. $\square$

**Checked:** mpmath 1.3.0, 110 digits (sympy 1.14.0 for the series of $\ell$: $\ell^2 = \frac{\theta^2}4 - \frac{\theta^4}{72} + \frac{\theta^6}{2880}$). (c): $\lVert(J_l - I)\rho\rVert = \ell(\theta)\lVert\rho_\perp\rVert$, 35 random $(\varphi, \rho)$ with
$\theta \in \{10^{-6}, 10^{-3}, 0.1, 0.5, 1.5, 3, 5\}$: relative $7\times10^{-88}$; $\max\lVert(J_l - I)\rho\rVert/(\frac\theta2\lVert\rho\rVert) = 0.9997$. LM illustration (110 digits, no `f64`): 6 points $p_i \sim \mathcal N(0, 2^2)$,
$r_i(X) = Xp_i - q_i$, $q_i = X^\star p_i$ plus noise $\mathcal N(0, 0.05^2)$, start $0.508$ rad and $1.12$ m from $X^\star$, undamped Gauss–Newton, $J_C$ by differencing $r \circ \mathrm{ret}_X$. Iteration 1: $J_{\mathrm{Scr}} - J_{\mathrm{Dec}} = 0$;
$\lVert\delta_{\mathrm{Scr}} - \delta_{\mathrm{Dec}}\rVert = 4\times10^{-112}$; $\lVert J_{\mathrm{WT}} - J_{\mathrm{Dec}}\mathrm{diag}(I, R^\top)\rVert_{\max} = 2.7\times10^{-81}$; $\lVert\delta_{\mathrm{WT}} - \mathrm{diag}(I, R)\delta_{\mathrm{Dec}}\rVert = 1.6\times10^{-81}$;
$\lVert t_{\mathrm{Scr}} - t_{\mathrm{Dec}}\rVert = 0.272$ against $\frac12\lVert\varphi\rVert\lVert\rho\rVert = 0.276$; $\lVert r\rVert^2$ at the trial point $0.762$ (`Screw`) and $0.0900$ (`Decoupled`, `WorldTranslation`). The three sequences reach the same minimum
($\lVert r\rVert^2 = 0.0275156852676$); $\lVert X_{\mathrm{Scr}} - X_{\mathrm{Dec}}\rVert_{\max}$ after iterations 1 to 7: $0.18$, $2.8\times10^{-3}$, $6.9\times10^{-7}$, $2.7\times10^{-11}$, $3.5\times10^{-15}$, $4.0\times10^{-18}$, $1.3\times10^{-21}$ (from iteration 2
the bases differ, so the bound of (c) no longer applies); `WorldTranslation` and `Decoupled` agree to $3\times10^{-83}$ at every iteration. One damped step from the start, $\lambda = 1$, residual components weighted $(1, 2, 3)$:
$(J^\top J + \lambda I)$: $\lVert X_{\mathrm{WT}} - X_{\mathrm{Dec}}\rVert_{\max} = 5\times10^{-82}$; $(J^\top J + \lambda\,\mathrm{diag}(J^\top J))$: $0.17$ (unweighted, the translation block of $J^\top J$ is $6I$ and both agree: a coincidence of that problem).
Scripts not committed. **Permanent:** none specified; `chart_jacobians_match_dual_*` checks the $J_C$ of (a) individually.

## 5. S²

`NUMERICS.md` §8: $\nu = n + \varsigma e_z$ with $\varsigma = \mathrm{sgn}(n_z)$, $\mathrm{sgn}(0) = +1$; $H = I - 2\nu\nu^\top/\nu^\top\nu$; $b_1 = He_x$, $b_2 = He_y$, $B = [b_1\ b_2]$; $K = \begin{bmatrix}0 & -1\\ 1 & 0\end{bmatrix}$.

**Proposition CH.8 (the basis).** For unit $n$:

- (a) $\nu^\top\nu = 2(1 + \lvert n_z\rvert) \ge 2$, and $\nu_z = n_z + \varsigma$ has modulus $1 + \lvert n_z\rvert \ge 1$: the two terms have the same sign, nothing cancels.
- (b) $H$ is symmetric, orthogonal, $\det H = -1$, and $He_z = -\varsigma n$.
- (c) $B^\top B = I_2$, $B^\top n = 0$, $BB^\top = I - nn^\top$, and $b_1 \times b_2 = \varsigma n$; equivalently $[n]_\times B = \varsigma BK$ and $B^\top[n]_\times = \varsigma KB^\top$.
- (d) $H = I - \nu\nu^\top/(1 + \lvert n_z\rvert)$, so $b_1 = \big(1 - \frac{n_x^2}{1 + \lvert n_z\rvert},\ -\frac{n_xn_y}{1 + \lvert n_z\rvert},\ -\varsigma n_x\big)$ and $b_2 = \big(-\frac{n_xn_y}{1 + \lvert n_z\rvert},\ 1 - \frac{n_y^2}{1 + \lvert n_z\rvert},\ -\varsigma n_y\big)$.
- (e) *Norm.* (a)–(d) hold for $\lVert n\rVert = 1$; `NUMERICS.md` §12 states no unit-norm tolerance for $S^2$ (index, open items). With $\lVert n\rVert = 1 + \varepsilon$ the $H$ form, with $\nu^\top\nu$ computed, keeps $B^\top B = I$ to
  rounding and has $B^\top n = O(\varepsilon)$; the closed form (d), which spends $\nu^\top\nu = 2(1 + \lvert n_z\rvert)$, is off by $2\varepsilon$ in both.

*Proof.* (a) $\nu^\top\nu = 1 + 2\varsigma n_z + 1$ and $\varsigma n_z = \lvert n_z\rvert$. (b) $H$ is a Householder reflection. $\nu^\top e_z = n_z + \varsigma = \varsigma(1 + \lvert n_z\rvert)$, so
$He_z = e_z - 2\nu\,\varsigma(1 + \lvert n_z\rvert)/(2(1 + \lvert n_z\rvert)) = e_z - \varsigma n - \varsigma^2e_z = -\varsigma n$. (c) The columns of an orthogonal matrix are orthonormal, so $b_1, b_2 \perp He_z = -\varsigma n$; for
$\det H = -1$, $b_1 \times b_2 = \det(H)\,He_z = \varsigma n$. Then $n \times b_1 = \varsigma b_2$, $n \times b_2 = -\varsigma b_1$. (d) $2/\nu^\top\nu = 1/(1 + \lvert n_z\rvert)$ and $\nu_z/(1 + \lvert n_z\rvert) = \varsigma$. $\square$

**Checked:** mpmath 1.3.0, 110 digits, 22 unit $n$: 6 random; $n_z \in \{0, \pm10^{-30}, \pm10^{-60}, \pm10^{-100}, \pm10^{-300}\}$; $n_z = 0$ with $\varsigma = +1$ and with $\varsigma = -1$ (both are valid, CH.11(d));
$\pm e_z$, $e_x$; two within $10^{-40}$ of $\pm e_z$. Every property of (a)–(d) (and $\nu^\top\nu \ge 2$, minimum $2.0$): $\le 3.3\times10^{-111}$. Exactly, with sympy 1.14.0: (a)–(d) and $H$-form $=$ closed form, for $\varsigma = \pm1$, modulo $n_x^2 + n_y^2 + n_z^2 = 1$.
`f64`, $2\times10^5$ random $n$ (a quarter with $n_z$ scaled by $10^{-1}$ to $10^{-11}$; numpy 2.5.3), the products of $B^\top B - I$ and $B^\top n$ evaluated in 80-bit `longdouble` (unit roundoff $2^{-64}$): $\max\lvert B^\top B - I\rvert = 4.3u$ ($H$ form), $4.5u$ (closed form);
$\max\lvert B^\top n\rvert = 4.7u$, $5.6u$ (measured maxima, not bounds). (e): $n$ scaled by $1 + \varepsilon$, $\varepsilon = 10^{-12}, 10^{-8}, 10^{-4}$, $2\times10^4$ of those $n$: $H$ form $\max\lvert B^\top B - I\rvert \le 6.6\times10^{-16}$ at all three,
$\max\lvert B^\top n\rvert = \varepsilon$; closed form both $2\varepsilon$. Scripts not committed.
**Permanent:** planned, corpus strata `s2:nz0`, `s2:generic` of `s2_retract`, `s2_local` (`PHASE5.md` §2), whose outputs depend on $B$; no basis-only check is specified.

**Proposition CH.9 (retract and local).** Let $\theta = \lVert\delta\rVert$, and for $m \ne -n$: $s = \lVert n \times m\rVert$, $w = n\cdot m$, $\alpha = \mathrm{atan2}(s, w) \in [0, \pi)$.

- (a) $n \oplus \delta = \mathrm{Exp}(B\delta)\,n = \cos\theta\,n + \frac{\sin\theta}\theta\,(B\delta \times n) = \cos\theta\,n - \varsigma\,\frac{\sin\theta}\theta\,BK\delta$. This is the exponential map of the round sphere applied to the tangent vector
  $B\delta \times n = -\varsigma BK\delta$, an isometry of $\mathbb R^2$ onto $T_nS^2$: the coordinates are the *rotation vector* $B\delta$, not the displacement.
  Its coefficients are $\frac{\sin\theta}\theta = 2k\cos\frac\theta2$ and $\cos\theta = 1 - 2(k\theta)^2$, functions of the two outputs $(k, \cos\frac\theta2)$ of `exp_coeffs` (`NUMERICS.md` §4), finite at $\delta = 0$
  ($k = \frac12$): no other coefficient is needed (`PHASE5.md` §2 names no source; this is the reading assumed here).
- (b) $\mathrm{loc}_n(m) = B^\top\mathrm{Log}(E_{nm}) = \frac\alpha s\,B^\top(n \times m) = \varsigma\,\frac\alpha s\,KB^\top m$, $E_{nm} = \mathrm{Exp}(\alpha\,\hat m)$ the minimal rotation taking $n$ to $m$,
  $\hat m = n \times m/s$ ($\frac\alpha s := 1$ at $m = n$, the limit).
- (c) $\frac\alpha s = \frac12\,r(s^2, w) = r(s^2, 1 + w)$, with $r(s^2, w) = 2\,\mathrm{atan2}(s, w)/s$ (`NUMERICS.md` §4 names its argument $n^2$, here $s^2$). `NUMERICS.md` §8's "the ratio $\alpha/\lVert n \times m\rVert$ is §4's $r$" is half of $r$ at $w = n\cdot m$.
  `local` takes it from `log_ratio`, whose series arm gives the limit $\frac\alpha s = 1$ at $m = n$. Use the first form, masked as CO.16(c), (e) require: its series variable is $s^2/w^2 = \tan^2\alpha$, valid for $w > 0$ only, and the mask is on that ratio.
  The second form has a non-negative second argument and the series variable $s^2/(1 + w)^2 = \tan^2\frac\alpha2$, valid for $\alpha < \frac\pi2$; its mask must be on that ratio (CO.16(e)), not on $s^2$: at $\alpha = \pi - 10^{-4}$, $s^2 = 10^{-8}$ but the ratio is $4\times10^8$,
  and the 8-term series returns $-4.4\times10^{67}$ for $31414.93$. $1 + w$ is exact for $w \le -\frac12$ (Sterbenz), so it does not cancel; the second form loses because $\mathrm{atan2}(s, 1 + w) = \frac\alpha2$ needs $s^2 + w^2 = 1$, which the computed $(s, w)$ break by about $u$,
  and near the antipode $\mathrm{atan2}$ divides that by $s$: 3 to 5 times the relative error of the first form (CH.11(e)).
- (d) $\mathrm{loc}_n(\mathrm{ret}_n\,\delta) = \delta$ iff $\theta < \pi$, and $(1 - 2\pi/\theta)\delta$ for $\pi < \theta < 2\pi$; $\mathrm{ret}_n(\mathrm{loc}_n\,m) = m$ for every $m \ne -n$. At $m = -n$ every $\delta$ with $\lVert\delta\rVert = \pi$ is a preimage and $\hat m$ is undefined:
  the domain of `S2Chart::local` is $m \ne -n$ (`NUMERICS.md` §12), and $\mathrm{loc}$ has no limit there (along $m = -n\cos\epsilon + \hat t\sin\epsilon$, $\hat t \perp n$, it tends to $\pi B^\top(n \times \hat t)$).

*Proof.* (a) Rodrigues (SO.4(a)) with $\omega = B\delta \perp n$, $\lVert\omega\rVert = \theta$, $\omega\cdot n = 0$; $B\delta \times n = -[n]_\times B\delta = -\varsigma BK\delta$ by CH.8(c). (b) The minimal rotation is the one about $n \times m$ by the angle $\alpha$ between $n$ and $m$; its
$\mathrm{Log}$ is $\alpha\hat m \perp n$, so $BB^\top(\alpha\hat m) = \alpha\hat m$ and $\mathrm{ret}_n(B^\top\alpha\hat m) = m$. Conversely $\mathrm{ret}_n(\delta) = \cos\theta\,n + \sin\theta\,\hat\omega \times n$ ($\hat\omega = B\delta/\theta$) has $n \times m = \sin\theta\,\hat\omega$
and $n \cdot m = \cos\theta$, so $\alpha = \theta$, $\hat m = \hat\omega$ and $B^\top\alpha\hat m = \delta$ for $\theta < \pi$: (Inv). The last form is $B^\top[n]_\times = \varsigma KB^\top$. (c) $r(s^2, w) = 2\alpha/s$ by definition. The quaternion
$(1 + w, n \times m) = 2\cos\frac\alpha2\,(\cos\frac\alpha2, \sin\frac\alpha2\,\hat m)$ has $\mathrm{atan2}(s, 1 + w) = \frac\alpha2$, so $r(s^2, 1 + w) = \alpha/s$ (SO.5(c)). (d) For $\pi < \theta < 2\pi$ the rotation by $\theta$ about $\hat\omega$ is that by $2\pi - \theta$ about
$-\hat\omega$ (LG.2(c)); the limit: $n \times m = \sin\epsilon\,(n \times \hat t)$ and $\alpha \to \pi$. $\square$

**Checked:** mpmath 1.3.0, 110 digits. (a): 40 $(n, \delta)$, $\lVert\delta\rVert$ from $3\times10^{-8}$ to $9$: `mp.expm` of the hat matrix against both closed forms $1.4\times10^{-111}$. (b), (c): 60 $(n, m)$ (34 with $w < 0$, ten within $10^{-1..-8}$ of $-n$, ten within $10^{-3..-12}$ of $n$),
$\mathrm{loc}$ by its definition (`mp.expm` of $\alpha\hat m$, geometric $\mathrm{Log}$, projected on $B$; the twin of `NUMERICS.md` §14): the two closed forms of (b) $2.6\times10^{-104}$ and $2.0\times10^{-103}$ (near the antipode; conditioning, CH.10); (c) $r(s^2, 1 + w)$ against $\alpha/s$ $1.0\times10^{-103}$, $r(s^2, w) = 2\alpha/s$
exactly; e.g. $n = e_z$, $m = (\sin1, 0, \cos1)$: $\alpha/s = 1.188395106$, $r(s^2, w) = 2.376790212$. At $\alpha = \pi - 10^{-4}$ (60 digits): $s^2 = 1.0\times10^{-8}$, $1 + w = 5.0\times10^{-9}$, $s^2/(1 + w)^2 = 4.0\times10^8$, the exact arm $31414.92659 = \alpha/s$, the 8-term series arm of $r(s^2, 1 + w)$ $-4.37\times10^{67}$ and of $r(s^2, w)$ $-2.0000$ (CO.16(c)); `f64` $1 + w$ was exact for all of $2\times10^5$ random $w \in [-1, -\frac12]$. $\mathrm{ret}(\mathrm{loc}\,m) = m$ $8.2\times10^{-111}$; $\mathrm{loc}(\mathrm{ret}\,\delta) = \delta$, 40 samples $\lVert\delta\rVert \in [10^{-9}, 3.1]$, $1.9\times10^{-110}$; the rule of (d) beyond $\pi$, 20 samples,
$2.3\times10^{-110}$. `f64`, direct evaluation with $H$ and the platform $\sin$, $\cos$, 20000 samples: $\mathrm{ret}$ within $3.9u$ componentwise of the formula at the same doubles, $\lvert\lVert m\rVert - 1\rvert \le 4u$. `mp.logm` of the minimal rotation, the reference `PHASE5.md` §2 names
for `s2_local`: exact ($\le 10^{-50}$) for $\alpha \le 3.02$ and complex or non-principal (deviation $3.1$) in 20 of 20 samples per angle for $\alpha \in \{3.03, 3.04, 3.05, 3.1, 3.14\}$, at 50 digits, identically in mpmath 1.3.0 and 1.4.1 (index, open items).
Scripts not committed. **Permanent:** planned, corpus `s2_retract`, `s2_local` (strata `s2:generic`, `s2:near-antipode`, `theta:*`), and the `*_matches_reference` proptest of `S2Chart::local` (`NUMERICS.md` §14; no name is fixed).

**Proposition CH.10 (Jacobians on $S^2$).** Let $m = n \oplus \delta$, $\omega = B\delta$, $B_m$ the Householder basis at $m$ with sign $\varsigma_m$, $E = \mathrm{Exp}(\omega)$, $Z = (EB)^\top B_m$, and
$A(\delta) = \frac{\sin\theta}\theta I + b(\theta)\,\delta\delta^\top = I + b\,(\delta\delta^\top - \theta^2I)$, $b$ of `NUMERICS.md` §4.

- (a) For every $\delta$, $\mathrm{rj}_n(\delta) = B_m^\top J_l(\omega)B = Z^\top A(\delta)$, with $Z \in \mathrm O(2)$, $\det Z = \varsigma\varsigma_m$, so $\det\mathrm{rj} = \varsigma\varsigma_m\frac{\sin\theta}\theta$; $A$ is symmetric with eigenvalues $1$ (along $\delta$) and $\frac{\sin\theta}\theta$:
  positive definite for $\theta < \pi$ (the domain $U_0$ of CH.2(b); there $\mathrm{sgn}\det\mathrm{rj} = \varsigma\varsigma_m$), singular at $\theta = \pi$ ($\mathrm{ret}_n(\delta) = -n$ for every direction: $\mathrm{rj}$ has rank one), indefinite for $\pi < \theta < 2\pi$ (there $\mathrm{sgn}\det\mathrm{rj} = -\varsigma\varsigma_m$).
- (b) $\mathrm{lj}_n(m) = A(\delta)^{-1}Z$ with $\delta = \mathrm{loc}_n(m)$, $\theta = \alpha$, $A(\delta)^{-1} = \frac\alpha s\,(I - b\,\delta\delta^\top)$, $b = b(\alpha)$ (no $\hat\delta = \delta/\theta$, which is $0/0$ at $m = n$): singular values $1$ and $\frac\alpha{\sin\alpha} = \frac\alpha s$, so $\lVert\mathrm{lj}\rVert_2 = \alpha/\sin\alpha \to \infty$ as $m \to -n$.
- (c) $E = I + [v]_\times + [v]_\times^2/(1 + w)$ with $v = n \times m$, $w = n \cdot m$: $Z$ needs no transcendental function. It is the form for $\mathrm{lj}$, where only $(n, m)$ are at hand. $\mathrm{rj}$ has $\delta$: take $E = \mathrm{Exp}(B\delta)$, the rotation `retract` already forms;
  the form above divides by $1 + w \approx \frac12(\pi - \theta)^2$ and loses accordingly (Checked).
- (d) In first-order form, $\mathrm{loc}_m(m + \mathrm dz) = B_m^\top(m \times \mathrm dz)$ for $\mathrm dz \perp m$.

*Proof.* (d) $\mathrm D\,\mathrm{ret}_m(0)\xi = B_m\xi \times m$ (CH.9(a)); for $\mathrm dz = B_m\xi \times m$, $m \times \mathrm dz = B_m\xi$. (a) By LG.8 (left form, $N = 0$), $\mathrm{ret}_n(\delta + \eta) = \mathrm{Exp}(J_l(\omega)B\eta)\,m + O = m + \zeta \times m + O$, $\zeta = J_l(\omega)B\eta$, so by (d)
$\mathrm{loc}_m$ gives $B_m^\top(m \times (\zeta \times m)) = B_m^\top\zeta$ ($B_m^\top m = 0$). Now $EB$ is an orthonormal basis of $m^\perp$ ($Eb_i \perp En = m$), so $B_m = EBZ$, and $E^\top J_l(\omega) = J_r(\omega)$ (LG.9(b)): $B_m^\top J_lB = Z^\top B^\top J_r(\omega)B$.
With $W = [\omega]_\times$, $B^\top WB = 0$ (its entries are $\omega \cdot (b_j \times b_i) \in \{0, \pm\varsigma\,\omega\cdot n\} = \{0\}$) and $B^\top W^2B = B^\top(\omega\omega^\top - \theta^2I)B = \delta\delta^\top - \theta^2I$, so $B^\top J_rB = I + b(\delta\delta^\top - \theta^2I)$; and $b\theta^2 = 1 - \frac{\sin\theta}\theta$.
$\det Z$: $Eb_1 \times Eb_2 = \varsigma m$ and $b_1' \times b_2' = \varsigma_mm$ (CH.8(c)). (b) CH.2(b), $Z^{-\top} = Z$; eigenvalues of $A^{-1}$; $\theta = \alpha$ by CH.9(b). With $\sigma = \frac{\sin\theta}\theta$ and $b\theta^2 = 1 - \sigma$, $A\,(I - b\,\delta\delta^\top) = \sigma I + b\,\delta\delta^\top(1 - \sigma - b\theta^2) = \sigma I$, and $1/\sigma = \alpha/s$. (c) $E = \mathrm{Exp}(\alpha\hat m) = I + \sin\alpha[\hat m]_\times + (1 - \cos\alpha)[\hat m]_\times^2$ with
$\sin\alpha = s$, $\cos\alpha = w$, $\frac{1 - w}{s^2} = \frac1{1 + w}$. $\square$

$A = I - b\,(\theta^2I - \delta\delta^\top)$ has no cancellation and needs the one coefficient $b$ (`jr_coeffs`); $A^{-1}$ needs $b$ and the ratio $\alpha/s$ of CH.9(c), and no other coefficient. The `f64` accuracy of $\mathrm{lj}$ near the antipode is not established here (CH.11(e) is that of $\mathrm{loc}$).

**Checked:** mpmath 1.3.0, 110 digits, 40 $(n, \delta)$, $\lVert\delta\rVert \in [10^{-6}, 3]$, eight with $\lvert n_z\rvert \in [10^{-4}, 10^{-1}]$ of both signs; 23 have $\varsigma \ne \varsigma_m$. $\mathrm{rj}$ by central differences ($h = 10^{-30}$) of $\eta \mapsto \mathrm{loc}_m(\mathrm{ret}_n(\delta + \eta))$ ($B_m$ frozen)
against $B_m^\top J_l(B\delta)B$ and against $Z^\top A$: $9.2\times10^{-62}$ (both); $\mathrm{lj}$ by differencing $\eta \mapsto \mathrm{loc}_n(\mathrm{ret}_m(\eta))$ against $A^{-1}Z$: $2.9\times10^{-59}$ (conditioning, up to $\alpha/\sin\alpha = 21$); $\mathrm{lj}\,\mathrm{rj} = I$ $3.1\times10^{-58}$;
singular values of $\mathrm{lj}$ against $\{1, \theta/\sin\theta\}$ $4.1\times10^{-58}$; $\det Z - \varsigma\varsigma_m$ $2.5\times10^{-111}$, $Z^\top Z - I$ $5.0\times10^{-111}$, the sign of $\det\mathrm{rj}$ equal to $\varsigma\varsigma_m$ in all 40; (c) $2.2\times10^{-109}$. Across $\theta = \pi$ (70 digits, $\mathrm{rj}$ by central differences, $h = 10^{-25}$, of $\mathrm{loc}_m\circ\mathrm{ret}_n$ with $\mathrm{loc}$ by its definition; three random $n$ and directions each):
$\det\mathrm{rj} = \varsigma\varsigma_m\frac{\sin\theta}\theta$ to the $8$ digits printed at $\theta \in \{1, 3, \pi - 0.01, 3.4, 4, 5\}$ (with $\frac{\sin\theta}\theta = -0.189$ at $4$, $-0.192$ at $5$); at $\theta = \pi$, $\lVert\mathrm{ret}_n(\delta) + n\rVert \le 1.6\times10^{-71}$ for three directions and the singular values of $\mathrm{rj}$ are $1$ and $1.2\times10^{-47}$ (zero to the differencing).
$A\,\frac\alpha s(I - b\,\delta\delta^\top) - I \le 8.6\times10^{-60}$ ($\theta = \alpha \in \{10^{-3}, 0.7, 2, 3.1\}$, 60 digits; the eigenvalues of $A$ at $3.1$: $1$ and $0.0134131 = \sin3.1/3.1$). (c) in `f64` (numpy 2.5.3, $H$ bases; $n$ a random unit double, $\delta$ of random direction with $\lVert\delta\rVert = \pi - \varepsilon$, 1500 samples per $\varepsilon$),
$Z$ formed from the `f64` $m = \mathrm{ret}_n(\delta)$ and $B_m$ by the $(v, w)$ form, against $Z = (EB)^\top B_m$ in mpmath (50 digits) from the same doubles $n$, $B$, $\delta$ ($E = \mathrm{Exp}(B\delta)$, $B_m$ the $H$ basis at the exact $m = En$): max entry error $1.7\times10^3u$, $1.6\times10^5u$, $1.6\times10^7u$, $1.5\times10^9u$ for $\varepsilon = 10^{-1..-4}$;
formed from $E = \mathrm{Exp}(B\delta)$ (Rodrigues, platform $\sin$, $\cos$): $9.5u$, $9.2u$, $8.5u$, $9.6u$. Scripts not committed. **Permanent:** planned,
`chart_jacobians_match_dual_*` (`PHASE5.md` §1.4, "S²: `Mat2<S>`"); none for (b)'s singular values or (c).

**Proposition CH.11 (the frame, and why the chart is frozen).**

- (a) For fixed $\varsigma$, $n \mapsto B(n)$ is smooth except at $n = -\varsigma e_z$ (where $\nu = 0$: the price of choosing $\varsigma = \mathrm{sgn}\,n_z$, which gives $\nu^\top\nu \ge 2$). The chosen $B$ has, at $n_z = 0$, the one-sided limits $B_\pm$ (frames of $\varsigma = \pm1$) with
  $B_+ - B_- = -2\,e_z\,(n_x, n_y)$, $\lVert B_+ - B_-\rVert_F = 2$, and $B_- = B_+T$, $T = \begin{bmatrix}n_y^2 - n_x^2 & -2n_xn_y\\ -2n_xn_y & n_x^2 - n_y^2\end{bmatrix}$: a reflection, $\det T = -1$, $\lVert I - T\rVert_2 = 2$.
- (b) No continuous $B$ exists on $S^2$ with $B^\top B = I$, $B^\top n = 0$ (hairy ball theorem [MilnorHB]: $b_1$ would be a continuous unit tangent field).
- (c) A step $\omega \perp n$ has coordinates $\delta_\pm = B_\pm^\top\omega$, $\delta_- = T\delta_+$; $\mathrm{Exp}(B_-\delta)n = \mathrm{Exp}(B_+T\delta)n \ne \mathrm{Exp}(B_+\delta)n$ for $0 < \lVert\delta\rVert < \pi$ unless $T\delta = \delta$ (a line of $\delta$). Hence a basis recomputed from the moving point inside $\oplus$, a line search or a Jacobian
  changes the meaning of the coordinates of one step by $T$ whenever the point crosses the equator. Frozen at the base of a linearization, $\mathrm{ret}_n$ and $\mathrm{loc}_n$ are smooth in $\delta$ and $m$; the discontinuity appears in exactly one place, the frame factor $Z$ of CH.10, whose $\det Z = \varsigma\varsigma_m$ is the
  honest transition between the chart frozen at $n$ and the one frozen at $m$, and is what moves a gradient or a covariance across a relinearization.
- (d) `NUMERICS.md` §8 says $\mathrm{sgn}(0) = +1$ and `PHASE5.md` §2 has a stratum "$n_z = \pm0$". For $n_z = -0.0$ the reading "$\mathrm{sgn}(0) = +1$" (`n_z >= 0`) gives $\varsigma = +1$ and `copysign(1, n_z)` gives $-1$; both give valid frames, related by $T$,
  so `retract` and `local` at $n_z = \pm0$ differ by $O(1)$ between the two readings.
- (e) *Accuracy near the antipode.* $\lVert\mathrm{lj}\rVert_2 = \alpha/\sin\alpha \approx \pi/(\pi - \alpha)$ (CH.10(b)) is the conditioning of $\mathrm{loc}_n$ in $m$: a relative error $u$ in $s$ or in $m$ is $u/(\pi - \alpha)$ in $\alpha/s$ and in $\mathrm{loc}$.
  In `f64` the evaluation of CH.9(b) at the given doubles errs by $C\,u/(\pi - \alpha)$ with $C$ of $1.4$ to $3.8$ for the first form of CH.9(c) and $5.0$ to $8.3$ for the second (Checked: expression and sampling); these are measured maxima of one harness, not bounds,
  and $C$ moves by a factor of about $2$ with which of $B^\top m$, $n \times m$ enters and with how $m$ is drawn. The factor $1/(\pi - \alpha)$ is the conditioning.

*Proof.* (a) At $n_z = 0$, $\nu^\top\nu = 2$ and the closed forms of CH.8(d) differ only in the sign of the third component, $\mp n_x$, $\mp n_y$; $T = B_+^\top B_-$ entrywise, $\det T = -(n_x^2 + n_y^2)^2 = -1$, trace $0$. (b) The theorem. (c) Substitute. (d), (e) Measured. $\square$

**Checked:** mpmath 1.3.0, 110 digits: at $n = (\cos\lambda, \sin\lambda, 0)$, $\lambda \in \{0.3, 1.9, 4.0\}$, $\lVert B_+ - B_-\rVert_F = 2$, $T$ and $B_- = B_+T$ to $8.3\times10^{-112}$, $\det T = -1$; $B(n)$ at $n_z = \pm10^{-50}$ within $1.0\times10^{-50}$ of $B_\pm$ (one-sided limits). sympy 1.14.0, exactly: $T$, $B_- = B_+T$, $\lVert B_+ - B_-\rVert_F^2 = 4$.
`f64`, numpy: $n = (0.6, 0.8, -0.0)$ gives $\varsigma = +1$ by `n_z >= 0` and $-1$ by `copysign`, $\lVert B_+ - B_-\rVert_F = 2.0$. (e), numpy 2.5.3: $n$ a random unit double, $t \perp n$ a random unit vector, $m$ = $-n\cos\epsilon + t\sin\epsilon$ normalised, $\pi - \alpha \approx \epsilon = 10^{-1..-8}$, 3000 samples per $\epsilon$; $s$ = `norm(cross(n, m))`, $w$ = `dot(n, m)`;
$\alpha/s$ by `arctan2(s, w)/s` (first form) or `2 arctan2(s, 1 + w)/s` (second); $\mathrm{loc}$ as $\varsigma\frac\alpha sKB^\top m$ or $B^\top(\frac\alpha s\,n \times m)$, with $B$ of CH.8(d) or the $H$ form (four combinations per basis); reference: mpmath (50 digits) at the same doubles $n$, $m$, $B$.
Max component error $\times(\pi - \alpha)/u$ over all combinations, bases and $\epsilon$: $[1.4, 3.8]$ (first form; $1.4$ to $1.9$ for $B^\top(\frac\alpha s\,n\times m)$, $2.3$ to $3.8$ for $\varsigma\frac\alpha sKB^\top m$) and $[5.0, 8.3]$ (second). The relative error of $\alpha/s$ alone, times $(\pi - \alpha)/u$: $[0.50, 0.86]$ against $[2.0, 2.9]$,
that is $3$ to $5$ times per $\epsilon$. Scripts not committed. **Permanent:** planned, strata `s2:nz0`, `s2:near-antipode` (`PHASE5.md` §2); the reading of $\mathrm{sgn}(-0)$ is not fixed by any of them (index, open items).

## 6. Sign audit

**Table CH.12.** The correct entry, the usual wrong variant and their **exact gap** ($\zeta$ the tangent the entry depends on; every gap tends to $0$ with $\zeta$, so a test needs $\zeta$ of order one, except three rows: the last, a discrete ambiguity; the S² $\mathrm{rj}$ row across the equator;
and the S² `local` factor, whose gap $\alpha/s$ tends to $1$ as $m \to n$: it is the easiest variant to catch, since any input exposes it).

| Where | Correct | Wrong variant | Gap | Reason |
|---|---|---|---|---|
| `Decoupled` $\mathrm{rj}$, rotation block | $J_r(\varphi)$ | $J_l(\varphi)$ | $2aW = W + O(\theta^2)$ | CH.5, SO.7 |
| `Decoupled` $\mathrm{rj}$, translation block | $\mathrm{Exp}(-\varphi)$ | $\mathrm{Exp}(\varphi)$; or $J_r(\varphi)$, the dual-matrix habit | $2\frac{\sin\theta}\theta W$, 2-norm $2\sin\theta$; $\mathrm{Exp}(-\varphi)(J_l(\varphi) - I)$, 2-norm $\ell(\theta) \approx \frac\theta2$ | the frame of $\rho$ is that of the *retracted* point, and $R_Y^\top R = \mathrm{Exp}(-\varphi)$ |
| `Decoupled` $\mathrm{lj}$, translation block | $\mathrm{Exp}(+\varphi_Y)$ | $\mathrm{Exp}(-\varphi_Y)$ | $2\frac{\sin\theta_Y}{\theta_Y}W_Y$ | the inverse of the $\mathrm{rj}$ block (CH.2(b)) |
| `WorldTranslation` translation block | $I$ | $\mathrm{Exp}(-\varphi)$ | $-W + O(\theta^2)$, 2-norm $2\sin\frac\theta2$ | $\rho$ is in the world frame; no rotation of the frame between $X$ and $Y$ |
| translation tangent of `WorldTranslation` against `Decoupled` | $\rho_{\mathrm{WT}} = R\rho_{\mathrm{Dec}}$ | the same $\rho$ | $(R - I)\rho$, $2\sin\frac{\theta(R)}2\lVert\rho\rVert$ for $\rho \perp$ axis; $0$ at $R = I$ | CH.5(c) |
| `Screw` against `Decoupled`, same $\delta$ | translation $t + RJ_l(\varphi)\rho$ | $t + R\rho$ | $R(J_l - I)\rho$, $\ell(\theta)\lVert\rho_\perp\rVert$ | CH.7(c); invisible at $\varphi = 0$ or $\rho \parallel \varphi$ |
| S² `retract` | $\mathrm{Exp}(B\delta)n$: displacement $B\delta \times n = -\varsigma BK\delta$ | $\exp_n(B\delta)$, displacement $B\delta$ | $\sqrt2\sin\theta$ in the point; a `local` that drops the quarter turn ($B^\top\mathrm{Log}_n m$) returns $\pm K\delta$, off by $\sqrt2\,\theta$ | CH.9(a): $\delta$ is a rotation vector, a quarter turn from the displacement |
| S² `local` factor | $\alpha/s = \frac12r(s^2, w) = r(s^2, 1 + w)$ | $r(s^2, w)$ (`NUMERICS.md` §8 as worded) | exactly a factor $2$ | CH.9(c) |
| S² `rj` | $Z^\top A(\delta)$ | $A(\delta)$, the frame of $m$ dropped | $(Z^\top - I)A$; $\det Z = -1$ when $n$, $m$ are on opposite sides of $z = 0$ | CH.10: the frame of $m$ is its own |
| S² `sgn` at $-0$ | one reading, stated | the other | $B_- = B_+T$, $\lVert I - T\rVert_2 = 2$ | CH.11(d) |

**Checked:** mpmath 1.3.0, 110 digits, $\varphi$ random axis, $\theta \in \{0.01, 0.1, 0.5, 1, 2\}$: $J_r - \mathrm{Exp}(-\varphi) = (\frac{\sin\theta}\theta - a)W + (b - a)W^2 = \mathrm{Exp}(-\varphi)(J_l - I)$ (form $\le 6.5\times10^{-111}$; 2-norm $\ell(\theta)$: $0.0499861$ at $\theta = 0.1$,
$0.893743$ at $2$); $\lVert\mathrm{Exp}(\varphi) - \mathrm{Exp}(-\varphi)\rVert_2 = 2\sin\theta$ (all five to $8$ digits), $J_l - J_r - 2aW$ $\le 2.5\times10^{-111}$, $\lVert\mathrm{Exp}(-\varphi) - I\rVert_2 = 2\sin\frac\theta2$. S² (20 samples, $\theta \in [0.01, 3]$): $\lVert\mathrm{ret} - \exp_n(B\delta)\rVert = \sqrt2\sin\theta$, $2.5\times10^{-111}$; 20 $(n, m)$: $B^\top\mathrm{Log}_n(m) = \pm K\delta$ to $4.2\times10^{-111}$ and $\lVert B^\top\mathrm{Log}_n(m) - \delta\rVert = \sqrt2\lVert\delta\rVert$ to $6.7\times10^{-111}$. The other rows are CH.5, CH.7, CH.9–CH.11 (`Checked:` there).
Scripts not committed. **Permanent:** as the propositions cited.
