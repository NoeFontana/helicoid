# Maths: the derivations behind `NUMERICS.md`

> **Non-normative.** [`NUMERICS.md`](../NUMERICS.md) is the contract and wins on every conflict.
> These pages derive why its formulas hold, where they are numerically delicate and where a sign
> is easy to lose. Code implements `NUMERICS.md`, never these pages. A disagreement is reported,
> not settled here: the corpus ([`PHASE1.md`](../PHASE1.md) §4.3), not a page, arbitrates, and a
> formula change is a `NUMERICS.md` edit together with a draft decision record
> (`docs/decisions/README.md`).

Every page uses the repository's conventions
([`0002`](../decisions/0002-one-convention-for-a-stack-that-already-disagrees.md)): Hamilton, `w`
first, active rotations; rotation-first tangents $[\varphi;\rho_1;\dots;\rho_N]$; right
perturbation by default, every perturbation naming its side; $a * b = T_{ax}T_{xb}$. Literature
(Barfoot, Solà) is translation-first: permute, never copy. For $N = 1$ the permutation is
$\tau_{\mathrm{tf}} = P\tau$, $P = \begin{bmatrix}0 & I_3\\ I_3 & 0\end{bmatrix}$, and a matrix $M$
acting on rotation-first tangents acts on translation-first ones as $PMP$; so
$\mathrm{Ad}_X = \begin{bmatrix}R & 0\\ [t]_\times R & R\end{bmatrix}$ here is
$\begin{bmatrix}R & [t]_\times R\\ 0 & R\end{bmatrix}$ in Solà et al.

## Notation

| Symbol | Meaning |
|---|---|
| $N$, $n$ | number of translation-like blocks; tangent dimension $n = 3 + 3N$ |
| $G$, $\mathfrak g$ | $\mathrm{SE}_N(3)$ ($N = 0$: $\mathrm{SO}(3)$) as $(3+N)$-square matrices, and its Lie algebra |
| $X$, $Y$, $R$, $x_i$ | group elements; $X = \begin{bmatrix}R & x_1\cdots x_N\\ 0 & I_N\end{bmatrix}$ (`NUMERICS.md` §1) |
| $\tau$, $\sigma$, $\delta$ | tangent vectors $[\varphi;\rho_1;\dots;\rho_N] \in \mathbb R^n$; $\delta$ is small |
| $\theta$, $W$ | $\theta = \lVert\varphi\rVert$ (rotation angle of a tangent); $W = \varphi^\wedge = [\varphi]_\times$, $[a]_\times b = a \times b$ |
| $\theta(X)$ | rotation angle of a group element, $\operatorname{tr}R = 1 + 2\cos\theta(X)$, $\theta(X) \in [0,\pi]$. "$\theta < 2\pi$" (a tangent) and "$\theta(X) < \pi$" (an element) are different domains |
| $\tau^\wedge$, $A^\vee$ | hat $\mathbb R^n \to \mathfrak g$ (`NUMERICS.md` §1) and its inverse |
| $\exp$, $\log$ | matrix exponential and principal matrix logarithm |
| $\mathrm{Exp}$, $\mathrm{Log}$ | $\mathrm{Exp}(\tau) = \exp(\tau^\wedge)$; $\mathrm{Log}$ its inverse on the canonical branch $\theta \in [0,\pi]$. $\mathrm{Exp}\,\delta$ abbreviates $\mathrm{Exp}(\delta)$ where unambiguous |
| $E$ | $\mathrm{Exp}(\tau)$ for the tangent $\tau$ of the row at hand (LG.14, LG.15) |
| $\mathrm{Ad}_X$, $\mathrm{ad}_\tau$ | adjoints of the group and of the algebra, as $n \times n$ matrices on tangent coordinates |
| $J_r$, $J_l$ | right and left Jacobians of $\mathrm{Exp}$; $J^{-1}$ is the matrix inverse |
| $\oplus_R$, $\ominus_R$, $\oplus_L$, $\ominus_L$ | $X \oplus_R \tau = X\,\mathrm{Exp}(\tau)$, $Y \ominus_R X = \mathrm{Log}(X^{-1}Y)$, $X \oplus_L \tau = \mathrm{Exp}(\tau)X$, $Y \ominus_L X = \mathrm{Log}(YX^{-1})$. A subscript or label $R$, $L$ names a side, never a matrix |
| $\mathrm D^s F$ | derivative of $F$ in the convention of side $s \in \{R, L\}$ (Definition LG.11); the $\partial/\partial$ of `NUMERICS.md` §2.3 |
| $A + \epsilon B$ | dual-block matrix (`NUMERICS.md` §2.2): diagonal blocks $A$, first block column $B_i$ |
| $\Gamma_m(\varphi)$ | $\sum_{k \ge 0} W^k/(k+m)!$ (`NUMERICS.md` §7); $\Gamma_1 = J_l(\varphi)$ |
| $u$ | unit roundoff (`NUMERICS.md` §2.1), except on `so3.md`, where $u$ is the quaternion's vector part (`NUMERICS.md` §1), $n = \lVert u\rVert$ (not the tangent dimension of the first row) and the roundoff is written $\mathsf u$; `coefficients.md` and `error-analysis.md` keep $u$ for the roundoff and write the vector part $q_{\mathrm v}$ |
| $\Xi$, $\Pi$, $T_M$ | on `se3.md` only: $[\rho]_\times$ for one translation block (so the $P$ above is $\Pi$ there), the block permutation $\Pi_N$ (SE.14; $\Pi_1 = P$), and CO.1's $\tau_M$ (renamed because $\tau$ is the tangent); it uses $\sigma_m$ as in CO.1, $n$ for a power (the tangent dimension is written $3 + 3N$), $\mu$ for $\varphi\times\rho$ (SE.7) and $s(\theta) = 2\sin(\theta/2)/\theta$ (SE.16, LG.16) |
| $K$, $\gamma$, $\Pi$, $\nu$ | on `so2-se2.md` only: the quarter turn $\begin{bmatrix}0&-1\\1&0\end{bmatrix}$ (not `so3.md`'s $K(M)$), $\gamma = \frac\theta2\cot\frac\theta2$, the $3\times3$ shift $\tau_{\mathrm{tf}} = \Pi\tau$ (not `se3.md`'s $\Pi_N$) and $\nu = \lVert\rho\rVert/2$; $\theta$ and $\theta(X) = \mathrm{atan2}(R_{10}, R_{00})\in[-\pi,\pi]$ are signed there: the $\theta(X)$ row above (branch $[0,\pi]$) and `NUMERICS.md` §1, §12 ("$\theta<2\pi$", "$\theta(X)<\pi$") are read as $\lvert\theta\rvert$, $\lvert\theta(X)\rvert$, $\alpha$, $\beta$ are `NUMERICS.md` §6's (not an angle or a relative rounding) and $a, b, c$ are §4's; label prefix `PL` |
| $\sigma_m$, $\tau_m$, $\rho_m$, $E_x$, $E_s$, $E_\times$, $s_j$ | on `coefficients.md` (and $\sigma_m$ on `se3.md`, as in the row above), not the tangents, translation tangent or $\mathrm{Exp}$ above: the series families (CO.1), a relative truncation size, the errors of the exact arm, of the series arm and at the crossing, and the coefficient of $z^j$; its symbol note lists the rest |
| $\mathcal F$, $\mathcal B$, $\mathcal B_M$, $\beta$, $e$, $\underline y$, $\nu_g$, $\zeta$, $\chi$, $\mathsf G$ | on `error-analysis.md` only: forward and backward error in units of $u$ (EA.3), $\mathcal B_M$ the Frobenius residual of `from_matrix`; the residual rotation vector $\mathrm{Log}(X^{-1}\mathrm{Exp}\,\hat\varphi)$ (not `so2-se2.md`'s $\beta$ nor LG.16's singular value) and the tangent error $\hat\varphi - \varphi$ (not `NUMERICS.md` §4's coefficient $e$, which EA.20 also names); the floor of the §11 metric; the accuracy factor of a library function ($\lvert\delta\rvert \le \nu_gu$); $\pi - \theta$; the measured constant of EA.10(c); splitmix64's increment. $\eta = 2\pi - \theta$ there (not `so3.md`'s $\eta$), $\varepsilon$ is a relative perturbation size (not the dual unit $\epsilon_i$), $\lambda$ the azimuth of EA.22 (and $\lambda_{\max}(K)$ of EA.9 is `so3.md`'s eigenvalue), $\hat n$ the axis (a sample of EA.22), $n$ the vector norm $\lVert q_{\mathrm v}\rVert$ or a sample count; label prefix `EA` |
| $\gamma$, $X(t)$, $\Delta$, $d$, $E_t$, $E$, $E_s$, $\hat q$, $h$, $m$, $o$, $\Theta$, $\nu$, $\Lambda$, $\varkappa$, $\bar m$, $\varpi_t$, $\mu_s$, $M_s$, $\mathbf t$ | on `geodesics.md`, and not the same symbols elsewhere: $\gamma$ and $\nu$ are not `so2-se2.md`'s ($\frac\theta2\cot\frac\theta2$, $\lVert\rho\rVert/2$), $\Lambda$ is not `se3.md`'s $\Lambda_i$ (SE.16), $\mu_s(\theta)$ is not its $\mu = \varphi\times\rho$ (SE.7), $\Delta$ is not `error-analysis.md`'s rotation error, and $E_s = \mathrm{Exp}(s\,\mathrm{Log}\,E)$ (a rotation) is not `coefficients.md`'s error $E_s$ ($E = R_0^\top R_1$ in GE.5, not the $E$ of the row above, which GE.6 uses as $E = E_t$). The geodesic map $\gamma(X_0, X_1, t) = X(t)$, $\Delta = X_0^{-1}X_1$, $d = \mathrm{Log}\,\Delta$, $E_t = \mathrm{Exp}(t\,d)$ and $\alpha = \theta/2$ (as in SO.5); the dual quaternion $\hat q = q_{\mathrm r} + \epsilon q_{\mathrm d}$ (the $\epsilon$ of `NUMERICS.md` §2.2 for one block); the screw's axial displacement $h$ (not a differencing step), its axis point $o$, moment $m = o\times\hat n$, dual angle $\Theta = \theta + \epsilon h$, dual axis $\nu = \hat n + \epsilon m$ and $\Lambda = \tfrac12\Theta\nu$; $\varkappa = q_{\mathrm d,w}/\lVert q_{\mathrm v}\rVert^2$ (the `k` of `screw_pow`, not `NUMERICS.md` §4's $k$), $\bar m = m\sin\alpha$, the slerp weight $\varpi_t = \sin(t\alpha)/\sin\alpha$; $\mu_s(\theta)$, $M_s$ of GE.5; $\mathbf t$ the pose translation in §6–§7 (there $t$ is the parameter), $s$ the parameter in §4 (there $t$ is a translation), and $q_{\mathrm r}, q_{\mathrm d}$ the real and dual parts of $\hat q$ ($q_{\mathrm v}$ the vector part, as on `coefficients.md`); label prefix `GE` |
| $\mathcal M$, $\mathrm{ret}_X$, $\mathrm{loc}_X$, $\mathrm{rj}_X$, $\mathrm{lj}_X$, $\Phi$, $\Psi_Y$, $\ell(\theta)$, $\Lambda$ | on `charts.md`: a manifold; the `retract` and `local` of the chart frozen at $X$ and their Jacobians `retract_jacobian`, `local_jacobian` (CH.1); the transition between two charts at one base and its first-order value at the retracted point (CH.6); $\ell(\theta) = \lVert(J_l(\varphi) - I)\rho\rVert/\lVert\rho_\perp\rVert$ (CH.7; the $g(z) = (e^z - 1)/z$ of LG.7 is used there too); the LM damping scaling $\Lambda$ (a symmetric positive-definite matrix, with the scalar $\lambda$; not `geodesics.md`'s) |
| $n$, $m$, $\varsigma$, $\nu$, $H$, $B$, $K$, $s$, $w$, $\alpha$, $E$, $Z$, $A(\delta)$, $T$ | on `charts.md` §5 only, as in `NUMERICS.md` §8: $n$, $m$ unit vectors of $S^2$ (not the tangent dimension); $\varsigma = \mathrm{sgn}(n_z)$; $\nu$ the Householder vector and $H$ the Householder matrix (not `geodesics.md`'s dual axis and group element); $B = [b_1\ b_2]$; $K$ the quarter turn of `so2-se2.md`; $s = \lVert n\times m\rVert$ (not $s(\theta)$), $w = n\cdot m$ (not the quaternion's), $\alpha$ the angle between them; $E = \mathrm{Exp}(B\delta)$ (or the minimal rotation $E_{nm}$); $Z$ the frame change between the bases at $n$ and at $m$; $A(\delta) = \frac{\sin\theta}\theta I + b\,\delta\delta^\top$ (not a dual-block matrix); $T$ the frame jump at $n_z = 0$; label prefix `CH` |
| $\omega$, $f$, $\Delta t$, $\Delta R$, $\Delta v$, $\Delta p$, $\Omega$, $\mathcal J_m$, $\mathcal N_R$, $\mathcal N_L$, $\Sigma_R$, $\Sigma_L$, $\xi$, $D$, $H$, $E_\Sigma$, $\Delta_H$, $\mathsf s_\pm$ | on `gamma-gaussian.md` only: the body-frame angular rate and specific force (accelerometer output) of `NUMERICS.md` §7 ($f$ is its $a$, renamed because $a$ is §4's coefficient), the interval and the increments in the frame at its start (not `geodesics.md`'s $\Delta = X_0^{-1}X_1$), the $5\times5$ generator of GG.4 (not a tangent), the directional Jacobian $\partial(\Gamma_m(\varphi)v)/\partial\varphi$; the right and left Gaussians of GG.7 and the covariances $\Sigma_R$, $\Sigma_L$ (the subscript names a side), a random tangent $\xi$; $D = \mathrm{diag}(\Sigma_{ii}^{1/2})$ and $H = D^{-1}\Sigma D^{-1}$ the correlation matrix (not `charts.md`'s Householder matrix), the backward error $E_\Sigma$ of the Cholesky factorization and the scaled total perturbation $\Delta_H = D^{-1}\Delta D^{-1}$ of $\Sigma$ (GG.11–GG.12; $F$, $G$ there are the triangular solve's, not a function; the size of the matrix, `NUMERICS.md` §15's $N$, is $n$ there), the singular values $\mathsf s_\pm$ of $\mathrm{Ad}_X$ (not the series $\sigma_m$ of CO.1). $M$ is both CO.1's coefficient index (GG.2(d), GG.5) and the number of series terms (GG.6(c)); $s$ is the integration variable of GG.2(c) and GG.4, the covariance scale of GG.9(b) ($\Sigma = s^2\Sigma_0$) and, as $\mathsf s_\pm$, a singular value; and $\sigma_\varphi$, $\sigma_\rho$ are standard deviations (GG.13–GG.14), not CO.1's family; label prefix `GG` |
| $P$, $M$, $P^+$, $Q_L$, $Q_R$, $E$, $m$, $d$, $n$, $\eta$, $\hat q$, $P^L$, $M^L$, $P_C$, $M_C$, $\delta_C$, $\Pi$ | on `ambient-jacobians.md`, and not the same symbols elsewhere: `PlusJacobian` and `MinusJacobian` of the right chart at a stored quaternion or $(q, t)$ (AJ.1; this $P$ is not the block permutation at the top of this page, this $M$ not `so3.md`'s $K(M)$ nor `geodesics.md`'s $M_s$), the pseudoinverse, the $4\times4$ matrices of $p \mapsto q\otimes p$ and $p \mapsto p\otimes q$, the embedding $E = [0; I_3]$ of the pure quaternions (not `charts.md`'s $E = \mathrm{Exp}(B\delta)$), the ambient dimension $m$ ($4$ or $7$) and the tangent dimension $d$ ($3$ or $6$; not `charts.md`'s unit vector, nor `geodesics.md`'s moment and logarithm), $n$ a count of retractions (AJ.5(b)), $\eta = \lVert q\rVert^2 - 1$ as in SO.14 but of a stored $q$, its normalization $\hat q$ (not `geodesics.md`'s dual quaternion), the Jacobians of the left chart, and Ceres': its two Jacobians and its tangent $\delta_C = \varphi/2$ (AJ.7); $\Pi$ is the row permutation of `Quat::from_xyzw` (not `se3.md`'s $\Pi_N$); $v$ is the vector part of $q$, as in `PHASE6.md` §1 (not a twist's translation). Three of `NUMERICS.md` §4's $(a, b, c, d, e)$ are in play with other meanings: $a$ is the action point of AJ.2's example ($\mathrm Df\,P = -R[a]_\times$) and a radial size in AJ.9(c), while $a$, $b$ are §4's coefficients in AJ.6(c)'s proof ($J_l(\varphi) - I = aW + bW^2$), and $d$ is the tangent dimension above, never §4's coefficient; label prefix `AJ` |
| $M$, $x$, $s$, $\mathsf V$, $v_j$, $w_j$, $\xi_k$, $\zeta_k$, $\Delta$, $r_\Delta$, $\Lambda$, $\chi$, $g_\ell$, $y$, $\mathsf Q$, $\kappa$, $\ell$, $\lambda$, $\Omega$, $\eta$, $\alpha$, $\psi$, $\varsigma$, $\gamma^\ell_j$, $i_j$, $p_k$, $U_k$, $T_k$ | on `sim3.md`, and not the same symbols elsewhere: $M = W + \sigma I$ (the $3\times3$ block of the Sim(3) hat matrix; not `so3.md`'s $K(M)$ nor `ambient-jacobians.md`'s $M$), $x = \sigma + i\theta$, the scale $s = e^\sigma$ (not LG.16's $s(\theta)$), the block $\mathsf V$ of $\mathrm{Exp}$ (`NUMERICS.md` §9's $W(\varphi,\sigma)$) with $\mathsf V = v_0I + v_1W + v_2W^2$ and $\Gamma_2(M) = w_0I + w_1W + w_2W^2$, $M^k = \sigma^kI + \xi_kW + \zeta_kW^2$ (not `coefficients.md`'s $\sigma_m$, $\tau_m$), $\Delta = \sigma^2 + \theta^2$ and $r_\Delta = \sqrt\Delta$, $\Lambda = \lvert g_1(x)\rvert^2$ (not `se3.md`'s $\Lambda_i$, `geodesics.md`'s $\Lambda$ nor `charts.md`'s), $\chi = \mathrm{Re}\,g_1(x)$, $g_\ell(x) = \sum_kx^k/(k+\ell)!$ ($g_1$ is LG.7's $g$; $\Gamma_\ell(\varphi) = g_\ell(W)$ is GG.2's $\Gamma_m$), $y = \mathsf V\rho$ (the translation of $\mathrm{Exp}$), $\mathsf Q$ the $(\rho, \varphi)$ block of $J$ (SE(3): §5.3's $Q$), $\kappa = \sigma/\theta$ a ray (the condition number stays $\kappa_2$), $\ell$ the index of $g_\ell$; $\mu = W\rho$, $\nu = W\mu$ ($\mu = \varphi\times\rho$ as in SE.7); the second tangent is $\psi = [\alpha;\eta;\varsigma]$ as in SE.4 (SM.9 only), with $\alpha$ a rotation block (not `charts.md`'s angle), $\eta$ a translation block (not `ambient-jacobians.md`'s $\lVert q\rVert^2 - 1$), $\varsigma$ its scale (not `charts.md`'s $\mathrm{sgn}(n_z)$) and $\Omega = \alpha^\wedge$ (not `gamma-gaussian.md`'s $5\times5$ generator); $\lambda$ the angle of a ray (SM.6(c)), $\gamma^\ell_j$ the coefficients of $\Gamma_\ell(M)$ ($v_j = \gamma^1_j$, $w_j = \gamma^2_j$), $i_0, i_1, i_2$ those of $\mathsf V^{-1}$ (SM.8), $p_k = \mathrm{Re}\,x^k$, $U_k$, $T_k$ Chebyshev polynomials. $K$ is the number of series terms (SM.6(b), SM.7(d)), not `so2-se2.md`'s quarter turn, and $n$ is a weighted total degree (SM.6(b)), not the tangent dimension of the first row; label prefix `SM` |
| $O(\lVert\delta\rVert^2)$ | remainder bounded by $C\lVert\delta\rVert^2$ as $\delta \to 0$, $C$ locally uniform in the base point |

## Labels and the `Checked:` line

- A statement is labelled `<page>.<n>` (`LG` for [`lie-groups.md`](./lie-groups.md), `SO` for
  [`so3.md`](./so3.md), `CO` for [`coefficients.md`](./coefficients.md), `SE` for
  [`se3.md`](./se3.md), `PL` (planar) for [`so2-se2.md`](./so2-se2.md), `EA` (error analysis) for
  [`error-analysis.md`](./error-analysis.md), `GE` (geodesics) for [`geodesics.md`](./geodesics.md), `CH` (charts) for [`charts.md`](./charts.md), `GG` (Γ and Gaussians) for [`gamma-gaussian.md`](./gamma-gaussian.md), `AJ` (ambient Jacobians) for [`ambient-jacobians.md`](./ambient-jacobians.md), `SM` (Sim(3)) for [`sim3.md`](./sim3.md)), numbered in order of appearance. Labels are stable:
  later results are appended, never renumbered, so other pages and PR descriptions can cite them.
- Each result states its domain and carries a proof, or says which step is only outlined.
- Each proposition ends with one **Checked:** line: what was run, at which precision, over which
  sample, and whether the script is committed; then **Permanent:** the corpus function id
  (`PHASE1.md` §4.3) or named proptest (`PHASE3.md` §9, or the record that names it) that keeps
  checking it. It reads `planned` until the code it tests exists, and `none` when no permanent
  check is specified.
- A `Checked:` run computes from definitions (`mp.expm`, the geometric $\mathrm{SO}(3)$ inverse of
  $\mathrm{Exp}$, the defining series, a high-precision inverse, numerical differentiation of the
  operation itself), never from the closed form under test (`PHASE1.md` §2, item 2).
- The scripts are **not committed**. A figure records one scratch run at the stated
  precision and sampling: it is not a gate and cannot be rerun from the repository. Only a
  `Permanent:` item is a check. Where a sampled figure is not a proved bound, the line says so.
- The mpmath version is that of the run, stated in the `Checked:` line: `lie-groups.md` ran at
  1.4.1, `so3.md` and `coefficients.md` at 1.3.0 (the latter with sympy 1.14.0 for the exact
  rational series), `se3.md` at 1.4.1 (sympy 1.14.0, run with mpmath 1.3.0, for its exact
  identities), `so2-se2.md` at 1.3.0 (sympy 1.14.0), `error-analysis.md` at 1.3.0 (numpy 2.5.3; the
  `libm` 0.2.16 crate; wasmtime 49.0.0 for its wasm32 run), `geodesics.md` at 1.3.0 (sympy 1.14.0), `charts.md` at 1.3.0 (sympy 1.14.0 for its exact identities; numpy 2.5.3 for its `f64` rows), `gamma-gaussian.md` at 1.3.0 (sympy 1.14.0 for its exact identities; numpy 2.5.3 for its `f32` rows), `ambient-jacobians.md` at 1.4.1 (sympy 1.14.0 for its exact identities; numpy 2.5.3 for its `f64` rows; `pyceres` 2.6 and `gtsam` 4.3.0 for the two libraries), `sim3.md` at 1.3.0 (sympy 1.14.0 for its exact identities; numpy 2.5.3 for its `f64` rows). It is not a pin, and no claim on these
  pages depends on it; the claims about mpmath's own behaviour (`mp.logm` of an $\mathrm{SO}(3)$, an
  $\mathrm{SE}_N(3)$ (an $\mathrm{SE}(3)$ one again in GE.13, an $\mathrm{SO}(3)$ minimal rotation in CH.9) and an $\mathrm{SE}(2)$ matrix, below) were checked on both 1.3.0 and 1.4.1.
- A wrong variant is characterised by its exact gap to the right formula (`lie-groups.md` §6),
  never by a sampled minimum: every wrong variant tends to the right one as the tangent tends to
  $0$, so a sampled minimum measures the sampler.

## Map

| Page | `NUMERICS.md` | API items (`API.md` §3, `PHASE3.md` §2, §5) |
|---|---|---|
| [`lie-groups.md`](./lie-groups.md) | §1, §2.2, §2.3, §5.1, §5.2, §5.4 (the algebra, given §5.3), §12 (`jr_inv` domain), §14 (`jl`, `*_jacobians`) | `LieGroup::{exp, log, adjoint, ad, jr, jr_inv, jl, jl_inv, rplus, lplus, rminus, lminus, *_jacobians, compose_jacobians, inverse_jacobian}`, `Side`, `SEn3Jac` |
| [`so3.md`](./so3.md) | §1 (quaternion), §2.4, §3.1–§3.6, §4 (definitions of $k, a, b, c, r$; the $c$ series), §11 (`Log` conditioning, matrix input), §12 (`jr_inv`, `from_wxyz_unchecked`), §14 (`act_many`) | `SO3::{exp, log, act, act_many, to_matrix, from_matrix, renormalize, adjoint, ad, jr, jr_inv, jl, jl_inv, act_jacobians}`, `Quat::{from_wxyz_unchecked, from_wxyz_normalized}` |
| [`coefficients.md`](./coefficients.md) | §4 (all seven coefficients and the $\cos\frac\theta2$ of §3.1: series, cancellation, switch-point magnitudes, derivatives), §2.1, §7 ($\Gamma_2$'s coefficients), §12 (the $\theta = 0$, $\pi$, $2\pi$ behaviour); `PHASE1.md` §2 item 3, §4.3, §6 | `coeffs::{exp_coeffs, jr_coeffs, jr_inv_coeff, q_coeffs, gamma2_coeffs, log_ratio, se2_coeffs}` (`pub(crate)`, `PHASE3.md` §3), `coeffs::generated`, `xtask thresholds` |
| [`se3.md`](./se3.md) | §1 (SE_N(3), rotation-first tangents), §2.2 (the dual-matrix algebra), §2.4 (SE(3) action Jacobians), §5.1–§5.5 ($\mathrm{Exp}$, $\mathrm{Log}$, $\mathrm{Ad}$, $\mathrm{ad}$, the $Q$ block, $J^{-1}$), §14 (`SEn3Jac`, `jr_inv`, `jl`); the order conversion of `0002` | `SEn3::{exp, log, adjoint, ad, jr, jl, jr_inv, jl_inv}`, `SEn3Jac::{mul, inverse, apply, apply_transpose, write_dense}`, `SE3::act_jacobians`, `Twist::{from_translation_first, to_translation_first}`, `coeffs::q_coeffs` |
| [`so2-se2.md`](./so2-se2.md) | §1 (SE(2) tangent), §2.3 (its rows, instantiated), §2.4 (SO(2), SE(2) action rows: proposed), §4 ($a, b, c$ and $\alpha, \beta$: series, singularities), §6, §12, §14 (missing SE(2) rows) | `SO2::{exp, log, adjoint, ad, jr, jl, jr_inv, jl_inv, act_jacobians}`, `SE2::{exp, log, adjoint, ad, jr, jl, jr_inv, jl_inv, act_jacobians}`, `coeffs::se2_coeffs` |
| [`error-analysis.md`](./error-analysis.md) | §2.1 ($u$), §3.2, §11 (the metric, its floors, the conditioning of `Log`, `from_matrix`, $J^{-1}$), §12 ($J^{-1}$ near $2\pi$, `Exp` at large $\theta$), D8, D16; `PHASE1.md` §2, §4.3, §4.4, §5, §8; `PHASE2.md` §3, §8 | `Dual<S, N>`, `Real::{select, branch}`, `xtask::conformance` metrics, `xtask envelope`, `conformance/generate` (splitmix64, stratum sampling, the precision budget) |
| [`geodesics.md`](./geodesics.md) | §10 (all), §12 (`geodesic`), §14 (the `SE3::geodesic` and `geodesic_jacobians` rows); `PHASE4.md` §1–§3 | `LieGroup::{geodesic, geodesic_jacobians, geodesic_velocity}`, `reference::geodesic`, `SE3::geodesic` (fast twin), `SO3::geodesic`, `Product` (provided methods), `Side` |
| [`charts.md`](./charts.md) | §1, §2.3 (the rows the group charts use), §4 ($b$, $r$), §5.3 ($Q$), §8 (all), §12 (`S2Chart::local`), §14 (its twin); `PHASE5.md` §1–§2; `0012` | `Chart::{at, base, retract, local, retract_jacobian, local_jacobian}`, `RightChart<G>`, `LeftChart<G>`, `Screw`, `Decoupled`, `WorldTranslation`, `ProductJac`, `S2`, `S2Chart` |
| [`gamma-gaussian.md`](./gamma-gaussian.md) | §7 (all), §2.3 (the rows the sides use), §12 and §15 (`chol`, the mask, the rounding of $d^2$), §14 (the `gamma_apply_jacobian` and `Gaussian::{to_left, to_right}` rows); `PHASE5.md` §4–§5 | `so3::{gamma1, gamma2, gamma_apply_jacobian}`, `coeffs::gamma2_coeffs`, `Gaussian::{to_left, to_right, propagate, mahalanobis_sq}`, `chol`, `solve_lower`, `Jac::sandwich` |
| [`ambient-jacobians.md`](./ambient-jacobians.md) | §1 (quaternion, SE_N(3)), §2.4 (the action Jacobian, times $P$), §3.1–§3.2 (`Exp` to first order, the scale invariance of `Log`), §3.6 and §12 (`renormalize`, `from_wxyz_unchecked`: the non-unit $q$), §5.1 (`Screw`'s `Exp`); `PHASE5.md` §1.3; `PHASE6.md` §1, §7; `0012` | `AmbientChart::{write_plus_jacobian, write_minus_jacobian}` (quaternion, `Screw`, `Decoupled`), `Chart` (what they extend), `Quat::{from_wxyz_unchecked, renormalize}`, `StridedMut` |
| [`sim3.md`](./sim3.md) | §9 (the group law and $\mathrm{Ad}$ checked; `Exp`, `Log`, the Jacobians, their coefficients and the joint limit derived and **proposed**), §1, §2.3 (every row holds for Sim(3): SM.2(d)), §2.4 (Sim(3) action rows: proposed), §4 ($a, b, d$), §12 (the `Sim3` row: proposed); `PHASE5.md` §3; `0005` | `Sim3::{exp, log, adjoint, ad, jr, jl, jr_inv, jl_inv, act, act_jacobians}`, a coefficient group for $v_1, v_2$ (and $w_j$, $\partial_zv_j$), `coeffs` (one-variable functions of $\sigma$), `xtask thresholds`, `Real::{exp, exp_m1}` (needed, absent from `PHASE2.md` §2) |

Not derived yet, so `NUMERICS.md` alone states them: every section outside the Map rows, notably
§9 (Sim(3), owed: `sim3.md` derives it and proposes the formulas, SM.1–SM.13, none adopted); §7 is placed by `gamma-gaussian.md` (GG.1–GG.6), §8 by `charts.md` (CH.8–CH.11), §10 by `geodesics.md` (GE.1–GE.15); §11 is placed by `error-analysis.md`
(EA.3–EA.13) and LG.16 (the conditioning of the exact $J$ and $\mathrm{Ad}$); `PHASE5.md` §5's `Gaussian` by GG.7–GG.14 and §15 by GG.11–GG.12 (the mask and the rounding of $d^2$; §15's own bounds are Higham's); `PHASE6.md` §1's ambient Jacobians by AJ.1–AJ.9 (§12 has no row for them: below). The switch points and
series lengths of §4 are generated (0004), not derived: `coefficients.md` derives the magnitudes to check
them against (CO.10).

## Open items for the normative documents

The gaps the corpus generator, `helicoid-linalg` and the lint found in the specs, with most of the items
below, are collected for decision, with a recommendation each, in the draft record
[`0015`](../decisions/0015-specification-gaps-found-while-building-the-instrument.md).

- `NUMERICS.md` §12 has no row for the Jacobians of `rminus`, `lminus` and `Log`
  (`rminus_jacobians`, `lminus_jacobians`), which exist only for $\theta(X) < \pi$ (LG.2(c),
  LG.14). `API.md` R6 asks every restricted-domain function to name its §12 entry. Adding the row
  is a `NUMERICS.md` edit and needs a record: it is gap NU.9 of the draft `0015`. The `local_jacobian` of a group chart is that Jacobian (CH.3), `S2Chart`'s is defined for $m \ne -n$ only (CH.10); `retract_jacobian` is a formula defined for every $\delta$, invertible on $U_0 = \{\theta < \pi\}$ (CH.2(b); `S2Chart`'s is singular at $\theta = \pi$, CH.10(a)). The same cut bounds `Gaussian::mahalanobis_sq` ($\theta(\mu^{-1}x) < \pi$, GG.10(d)); `NUMERICS.md` §12 has no row for it either.
- `NUMERICS.md` §3.2: the flip is "implemented as `copysign`", yet "at $w = +0$ nothing flips" and
  "$q$ and $-q$ return $\pm\pi\hat n$". A sign-bit `copysign` flips at $w = -0$ (`copysign(1, -0.0) = -1`;
  `w < 0` does not), so $q = (+0, u)$ and $-q = (-0, -u)$ both return $+\pi\hat n$; only the `w < 0`
  reading returns $\pm\pi\hat n$. Both are logarithms (SO.5(d)); §3.2 should pick one. GTSAM's `Logmap` for a quaternion negates when $\lnot(w > 0)$, so it flips at $w = +0$ too (AJ.8(c)).
- `NUMERICS.md` §3.3 and §1 ($R(q)$, used by `to_matrix`, `act_many`) differ by $\eta v$ on a non-unit
  $q$ (SO.2(d)), up to $2^{-40}\lVert v\rVert$ at the `from_wxyz_unchecked` bound: the `act_many` twin
  proptest (§14) needs normalized inputs.
- `PHASE1.md` §4.3 defines `so3_log` by `mp.logm` of $R(q/\lVert q\rVert)$ and `so3_from_matrix` as the
  "`mp.logm`-consistent quaternion". `mp.logm` returns a complex, non-principal result for a rotation of
  angle $\theta \ge 3.03$ (onset in $(3.02, 3.03]$ depending on the axis; 30 axes, 30 and 80 digits,
  identical in mpmath 1.3.0 and 1.4.1), so both ids are wrong at every `theta:pi-1e-k` stratum
  and at the top of `theta:1e0` ($[1, \pi - 0.1)$). An eigendecomposition-based principal logarithm
  agrees with the `atan2` form there (SO.5, `Checked:`). `sen3_log_n{1,2,3}` (`mp.logm` of the
  $(3+N)$-square matrix) fails the same way, with onset in $(3.02, 3.04]$ ($N = 1, 2, 3$, 8 axes each,
  40 digits, identical in mpmath 1.3.0 and 1.4.1; SE.3, `Checked:`), where the geometric $\mathrm{SO}(3)$
  logarithm with a $\Gamma_1$ solve agrees.
- `NUMERICS.md` §3.4 says "then normalize" without saying how, and gives the output no sign rule.
  (i) §3.6's Newton step reaches rounding level only for $\lvert\eta\rvert \lesssim 2^{-26}$ (`f64`,
  SO.15(a)): a matrix further than about $2\times10^{-9}$ (max entry; $\lvert\eta\rvert \le 6\varepsilon$,
  SO.12(d)) from $\mathrm{SO}(3)$ needs the `svd3` projection or a division by the norm, and from
  $\eta = 2$ the step returns $0$ (SO.14). (ii) The extraction fixes only the pivot component positive,
  so the result jumps between $q$ and $-q$ at pivot ties, for an exact rotation too (SO.12, *Sign*): a
  canonicalization (for instance $w \ge 0$, which moves the jump to $w = 0$) is a §3.4 decision.
- `PHASE2.md` §6 says `svd3` is signed with $\det U = \det V = +1$ and "only $\sigma_3$ may be negative",
  but not that $\lvert\sigma_3\rvert \le \min(\sigma_1, \sigma_2)$. `nearest_rotation` $= UV^\top$ needs it
  (SO.13(c)): with $\Sigma = \mathrm{diag}(1, 1, -5)$ it is at squared distance $36$ against $20$.
- `NUMERICS.md` §4 "Naive cancellation" ($\sim2u/\theta^2$, $6u/\theta^2$, $12u/\theta^2$,
  $24u/\theta^4$ and $24u/\theta^2$, $360u/\theta^4$): the exponents are right, the constants are
  not bounds, and the column is not built by one rule ($a$, $b$, $c$, $d$ count the dominant rounded
  operand once, $e$'s 360 counts all its operands, which would give 48 for $d$). Sampled maxima
  (6e4 samples, `f64`): $a$ 1.0, $b$ 6.0, $c$ 48 (spec 12), $d$ naive 36 and rewritten 45 (24), $e$
  367 (360) (CO.6); `so3.md` SO.8 quotes the $12/\theta^2$ of $c$ as an amplification, with a
  pointer here. Separately, the column is about values: through
  `Dual` every coefficient, $k$, $a$ and $r$ included, loses $\theta^{-(p+2)}$ (CO.14), which the
  spec does not state and which sets the switch of $k$, $a$, $r$ and the objective of 0004 item 1.
- `PHASE1.md` §6 spans the switch grid over "$\theta \in [10^{-8}, 1]$", but the predicted and
  measured optimum of $k$, $b$, $d$, $e$ (`f64`, $m = 8$; $e$: $\theta_s \approx 1.4$) and, in `f32`,
  of all six lies above $1$ (CO.10); a switch capped at $1$ costs $e$ a factor $6$ (value) and $7.5$
  (derivative) in `f64`, and up to $8\times10^3$ (derivative of $e$) in `f32`. `PHASE1.md` §6 also
  fixes $m \le 8$; with it, the derivative of $e$ cannot fall below $\sim2\times10^3u$ (`f64`). The
  committed sweep (`conformance/sweeps/thresholds-seeded.csv`, `f64`; $m = 8$, the cap, for every row but $e$)
  puts $k$ and $d$ at the top of the grid, $b$ one step below it (a unique optimum $1.3\%$ under the
  top, decided by two records), $a$ at $\theta = 0.85$ and $c$ at $0.75$; for $k$, $b$, $d$ its objective
  is the exact arm's error over the records with $\theta > 1$, which no switch on the grid changes.
  $e$ is one step below the top at $8.7\times10^3u$, the exact arm's derivative error at the one record
  at $\theta = 1$, which the strict $z < \text{switch}$ leaves on the exact arm; a grid top of
  $\operatorname{nextUp}(1)$ gives it $2.6\times10^3u$.
- `NUMERICS.md` §4 (Continuity) says $\lvert\text{series} - \text{exact}\rvert$ at a switch is at
  most the recorded max error; the two arms can each err by that much in opposite directions, so it
  is at most twice, and $1.4$–$1.7\times$ was measured at the crossing (CO.12):
  `branch_continuity_*` needs the sum.
- `PHASE1.md` §4.3 computes `coeff_*` "from the definition at 120 digits", and §4.4 has
  `theta:exact0` ($0/0$) and `theta:subnormal` ($\theta \sim10^{-310}$: $4L + 3 \approx 1243$ digits
  lost). §2 item 3's check (1% at 150 digits agreeing to 40 digits) fails for $e$'s definition from
  $\theta \approx 10^{-19}$ (39.2 digits at $10^{-20}$) and for its `mp.diff` derivative in $\theta$ from
  $\theta \approx 3\times10^{-18}$ ($37.0$ digits at $10^{-18}$) (CO.17): below, the reference must be
  the defining series.
- `PHASE1.md` §4.3 ("derivative by `mp.diff`" of `coeff_*`, input $\theta$), §6 (`Dual<S, 1>`) and
  `NUMERICS.md` §4 (branch variable $\theta^2$, $n^2$ for $r$) do not say in which variable the
  derivative is taken; $\mathrm d/\mathrm d\theta = 2\theta\,\mathrm d/\mathrm dz$ (CO.13), and the digits
  `mp.diff` keeps differ with it ($67$ in $\theta$, $79$ in $z$ at $10^{-12}$, CO.17).
- `NUMERICS.md` §4 names $n^2$ as the branch variable of $r$ and evaluates only the exact arm at the
  safe argument, but the series is in $s = n^2/w^2$ and valid for $w > 0$ only; its division by
  $w^2$ is non-finite at $w = 0$ (0003 item 3: "no arm is ever non-finite"), and §8's
  `S2Chart::local` has $w = n\cdot m$ of either sign, where the series returns $-2$ for $2\pi/n$
  ($n^2 = 10^{-6}$, $w = -1$) (CO.16). §8's "$\alpha/\lVert n \times m\rVert$ is §4's $r$" is $r/2$ (CH.9(c): it equals $r(\lVert n \times m\rVert^2, 1 + n\cdot m)$ exactly, a non-negative second argument).
- `NUMERICS.md` §12 has no row for the domain of `Exp` (and the other coefficient users), while
  `PHASE2.md` §2 promises a finite result for finite in-domain input: $\theta^2$ overflows for
  $\lVert\varphi\rVert > 1.3\times10^{154}$ (`f64`), $1.8\times10^{19}$ (`f32`) and the exact arm
  returns `NaN` (CO.15). Accuracy is a separate matter: $\theta u/2$ for `Exp` and $2\pi u/\eta$ for
  $J^{-1}$ of SO(3) near $2\pi$ are the conditioning of the input (EA.6, EA.10) and no stratum lies beyond
  $\pi$.
- `NUMERICS.md` §4 and `PHASE3.md` §3 make each call site "one `S::branch` over a tuple" while
  `generated.rs` holds "a `Switch` per coefficient": a group with one mask shares one $\theta_s$, which
  costs the worst member up to $\sim3\times$ ($m = 8$) or $\sim10\times$ ($m = 4$) over its own optimum
  (CO.18). Separately, §3.1 leaves the series of $\cos\frac\theta2$ to the generator: it is
  $\sigma_0(z/4)$ with $\rho_m$ $2m+1$ times $k$'s (CO.2, CO.8), and it, not $k$, limits the
  derivative switch of `exp_coeffs` ($1.1$ against $1.5$ at $m = 8$, CO.10).
- `PHASE1.md` §10's seeded defect "`b` by its definition" fits $\theta^{-p}$, $p \in [1.8, 2.2]$,
  over `theta:1e-8`…`theta:1e-2`, but $b$ is $0$ (relative error exactly $1$) below
  $\sqrt{6u} = 2.6\times10^{-8}$, inside `theta:1e-8`: simulated fits give $p = 1.91$–$1.94$ ($2.00$
  without that stratum), inside the window by $0.1$ (CO.6, *Checked*). The window should start above
  $2.6\times10^{-8}$, or the fit skip saturated strata. Measured by the self-test on the committed
  corpus: $p = 1.937$ (`value` field, seven strata), inside by $0.137$; the record maximum over `value`
  and `d_branch` fits $3.98$, so the field is a reading too (0014 (draft), question 7).
- `0004` item 1 and `PHASE1.md` §6 minimize each coefficient's own value and `Dual`-derivative error, but
  a consumer multiplies a coefficient by a word of size $\theta^p$: with the exact arms of $b, d, e$ alone,
  the value of $Q$ errs by $\lesssim(16.3/\theta + 20.8)u$ (measured $\le6.5u/\theta$, `f64`; $39u$ at a
  switch of $0.9$, measured $\le 9.4u$) and the translation column of `Exp` by $O(u)\lVert\rho\rVert$ at every
  $\theta \le 1$ (SE.15), against $6u\theta^{-2}$ for $b$ and $360u\theta^{-4}$ for $e$. The value strata of
  `sen3_exp_n*`, `sen3_jl_n*` are far less sensitive to the
  switch of $b$ and $e$ than `coeff_b`, `coeff_e`, so a switch that minimizes the coefficient is
  conservative for those values; whether the objective should be the consumer's is a `0004` question. The
  derivative through `Dual` was not examined.
- `NUMERICS.md` §7 gives $\Delta R$, $\Delta v$, $\Delta p$ without their frames: they hold in the frame at the start of the interval, with gravity, the initial velocity and $\frac12g\Delta t^2$ removed, for an inertial world frame (no Earth rate), no bias and $\omega$, $f$ constant (GG.4); it
  also uses $a$ for the acceleration and for §4's coefficient. It says the `Dual` directional Jacobians are "exact to rounding, Taylor branches included": the derivative through a Taylor branch is that of the truncated polynomial, worse than the value by $\approx 2M/\theta$ (GG.6(c)), and through an exact arm it loses $1.3$ to $1.8\,u/\theta$ ($m = 1$) or
  $6$ to $10\,u/\theta^2$ ($m = 2$): measured $28u$ at an illustrative switch of $1$ (8 terms), $620u$ at four terms and $0.12$, $10^{13}u$ at $10^{-8}$. The sweep of `0004` item 1 minimizes the coefficient's derivative error, not the Jacobian's (as for $Q$, SE.15).
- `PHASE5.md` §4: `gamma_apply_jacobian::<M>` does not say which `M` are supported; the catalogue gives $m \le 3$ as coefficients ($\sigma_5 = (d - 2e)/3$, GG.2), not as accuracy (on the exact arm $\Gamma_m$ loses $\theta^{-(m-1)}$: $\approx 2u/\theta^2$ for $\Gamma_3$, GG.3(c)), and a closed-form $\mathcal J_2$ needs $\tau_6$, outside §4 (GG.5), while $\mathcal J_1$ is two groups of §4 and a product — `q_coeffs` $(b, d, e)$ for $Q$ and `jr_coeffs` $(a, b)$ for $J_l$ (SE.9(c), GG.5(c)). The only corpus id is `so3_gamma2`; the Jacobian has a twin row (§14) but no id, and
  the twin's "dense sum" states no length or tolerance ($\lVert W^n\rVert/(n+m)!$ has to be summed to a $\theta$-dependent $n$). `NUMERICS.md` §12 has no row for $\Gamma_m$ (entire; the overflow of $z$ as for $\mathrm{Exp}$, CO.15(c)).
- `PHASE5.md` §5: `propagate(&self, j: &G::Jac, mean: G)` takes a `G::Jac` that is one type for $J_r$, $J_l$ and $\mathrm{Ad}$: the type checks the side of the `Gaussian`, not of `j`; a $J_l$ for a $J_r$ converts the output side for $\mathrm{Exp}$ (LG.9(b); V7 of GG.14), the input side for $\mathrm{Log}$ and both for a composition (GG.9), silently. `to_left`/`to_right` are exact (GG.8), so their twin
  (§14) compares two roundings, and its tolerance has to be per block, not one common $\kappa u$ (GG.13(c)). No API item converts a covariance between rotation-first and translation-first orders ($\Pi\Sigma\Pi^\top$, SE.14(d)); `Twist` converts tangents only. `Jac::sandwich` is specified as $J\Sigma J^\top$ (`PHASE3.md` §2) and `chol` reads one triangle, the caller owning symmetry (`NUMERICS.md` §15): whether `sandwich` symmetrises is not said, and $(\mathrm{Ad}\,\Sigma)\mathrm{Ad}^\top$ is asymmetric by up to $1.2u$ of its largest entry (GG.11, GG.13).
- `PHASE6.md` §1 and §7: the formulas of $P$, $M$, $P_{SE3}$, $M_{SE3}$ are right (AJ.3, AJ.6: 60 to 100 digits), but (i) nothing says what `write_minus_jacobian` returns for a non-unit stored $q$: the printed $M$ gives $MP = \lVert q\rVert^2I$ and is $\lVert q\rVert^2$ times the
  derivative of $\mathrm{loc}$ (AJ.5), off by up to $2^{-40}$ ($8192u$) inside the domain of `from_wxyz_unchecked`, and by an $\eta$ that a model puts at $\approx1.1\sqrt n\,u$ (rms) after $n$ independent unnormalized retractions but at $\approx0.7\,nu$ after $n$ identical ones, which meets that bound after $\approx10^4$ in `f64` (AJ.5(b)); `NUMERICS.md` §12 has no row for either method, which `API.md` R6 asks for,
  and §14 has no twin row for them (D6; whether they are primitives checked by the corpus, as the ids of `PHASE6.md` §1 suggest, is not said). (ii) Nothing says which of three readings a corpus
  reference takes at a non-unit $q$: the printed formulas, the derivative of the definition at $q$, or at $q/\lVert q\rVert$ (the reading of `so3_act`) differ by $256u$ or $128u$ at the `q:nonunit` stratum, $\lvert\eta\rvert = 2^{-45}$ (AJ.5(d)). **Which reading `so3_act` takes is now measured**: the shipped `SO3::act`, which is §3.3 as written, scores $510.3u$ at that stratum and `to_matrix()` $\times v$, the scaled rotation of §1, scores $257.8u$ — i.e. $2\cdot2^{-45}/u$ and $2^{-45}/u$ exactly, so the corpus's reference is the **normalized** reading, $R(q/\lVert q\rVert)v$, and neither shipped form is it. Adopting it would be a §3.3 edit and a record; until then the stratum is outside `act`'s unit-$q$ domain and that row is the open question, not a defect. (iii) `MinusJacobian` has no corpus id (`so3_plus_jacobian` and `se3_plus_jacobian`
  only, and neither is in the table of `PHASE1.md` §4.3), and $MP = I$ does not pin $M$: $M + zq^\top$ passes it; $Mq = 0$ (or $PM = I - qq^\top$) does (AJ.2(b)). (iv) Ceres' `PlusJacobian` is a different matrix, left and half angle ($P_C = 2PR^\top$, AJ.7), so "dominated where oracles exist (Ceres' `QuaternionManifold`)" needs a runner
  that converts with $R$ (the rounding of the product belongs to the converter) or compares with `LeftChart<SO3>` times $2$ (exact, AJ.9(a)); its `Minus` does not flip the sign at $y\cdot x < 0$ (read from the source at the cited commit and the tags 2.1.0 and 2.2.0, not run in a released binary: `pyceres` does not expose it, AJ.7). (v) The implementing types are not named: the quaternion chart, `Screw` and `Decoupled` have the matrices of AJ.3 and AJ.6; `WorldTranslation` and the left charts have their own (AJ.6(e)).
- `NUMERICS.md` §15.4 and the rustdoc of `chol` say a positive definite matrix "with $\mathrm{cond}(A)$ near $1/u$" may clear the mask. Measured (GG.11): the mask and the rounding of $d^2$ follow $\lambda_{\min}(H)$ of the correlation matrix $H = D^{-1}AD^{-1}$: no failure in $200$ samples at $\lambda_{\min}(H) \ge 10^{-15}$ with $\mathrm{cond}(A)$ up to $10^{38}$, failures from $10^{-16}$.
- `NUMERICS.md` §6 cites Solà et al. 2018, Appendix (SE(2)), for "$J_r$, $J_l$ and their inverses". Appendix C
  of arXiv:1812.01537v9 prints $\mathrm{Ad}$, $J_r$, $J_l$ and the action Jacobian, translation-first, and **no
  inverse**; the printed matrices equal PL.6, PL.7, PL.10 through $\Pi$ (PL.9(d)), while the rotation-first
  $J^{-1}$ of PL.8 is derived, not permuted from a source.
- `PHASE3.md` §3 has `se2_coeffs(θ²) -> (α, β)`, but $\beta = \theta a$ is odd in $\theta$ (PL.1): a function of
  $\theta^2$ cannot return it (the kernel needs $\theta$, or returns $a$). The SE(2) Jacobians also need
  $a, b$ (`jr_coeffs`) and $c$ (`jr_inv_coeff`) besides $\alpha$ (PL.7–PL.8), and no listed call-site group
  carries them; nor is there a $\gamma$, though $1 - \theta^2c$ has no relative accuracy near $\lvert\theta\rvert = \pi$ (PL.11(b): does `jr_inv_coeff` return
  $\gamma$ from the cot arm of $c$?). `PHASE1.md` §6 sweeps $\alpha$, $\beta$ but §4.3 has no coefficient id for them.
- `NUMERICS.md` §12 has no row for `SE2::log`, `SE2::jr_inv`, `SE2::jl_inv` (PL.11(d)); §14 has no twin row for a
  closed-form `SE2::jr_inv` (PL.8: the dense inverse of `SE2::jr`); §2.4 has no SO(2), SE(2) action rows, though
  `API.md` §3 lists `act_jacobians` for both (PL.10); §6 and `PHASE3.md` §6 specify no SO(2) normalization,
  as §3.6 does for quaternions (PL.3(e)).
- `PHASE1.md` §4.3 defines `so2_*`, `se2_*` as "analogous"; `se2_log` by `mp.logm` fails as `so3_log` does: the
  result is complex and non-principal from an onset in $(3.02, 3.03]$ (PL.5, `Checked:`).
- `NUMERICS.md` §11 gives no floor for a tangent output: `Log`'s $\varphi$ has norm $\theta$, and
  with the floor of "1 for rotations" the strata `theta:1e-k` measure absolute error (a total loss
  at $\theta = 10^{-12}$ reads $9\times10^{3}u$, not $9\times10^{15}u$; strata with
  $\theta < 10^{7}u$ cannot show the seeded `acos` defect of `PHASE1.md` §10; EA.4(b)). Nor for
  `theta:subnormal`, where a subnormal output has fewer than $p$ significant bits and $\mathcal F$
  measures its quantization unless the floor is at least the smallest normal number (EA.4(b)).
- `NUMERICS.md` §11 calls matrix inputs "ill-conditioned near $\pi$ in the axis". The map
  $M \mapsto \pm q_\star$ is uniformly conditioned, $\frac1{2\sqrt2}$ for every $\theta$ (EA.9(a));
  what is ill conditioned is the orientation of the axis (the `Log` cut, distance
  $\sqrt2(\pi - \theta)$) and the trace pivot. Reporting backward error only is justified by the
  reference being undetermined below $O(d_M)$ (`PHASE1.md` §4.3's "`mp.logm`-consistent quaternion"
  has no value for a non-orthogonal double matrix) and by the floor $\mathcal B_M \ge d_M/u$ (EA.9(c),
  (d)), not by conditioning. §11 does not define this backward error; $\mathcal B_M$ (Frobenius,
  absolute, EA.3) is a reading of "backward error only".
- `NUMERICS.md` §11 and `PHASE1.md` §5 do not say how $\mathcal B$ is evaluated (it needs
  $\mathrm{Exp}$ at the subject's output, which the corpus cannot hold; EA.5 gives a first-order
  form from the stored $\varphi$, valid only for $\lVert e\rVert \ll 1$ and after reducing the
  output to the logarithm nearest $\varphi$: an output on the other branch reads $2\pi$ otherwise
  and $\mathcal B$'s branch invariance, EA.8(b), is lost), nor how $\hat y - y$ is formed from a 30-digit reference: parsing
  it to `f64` first quantizes $\mathcal F$ to an ulp and blinds the exact no-regress bar (EA.19(b)).
- `PHASE1.md` §4.2, §4.4: an `f32` subject is not said to receive `f32`-exact inputs; if the
  binary64 corpus inputs are rounded, $\mathcal F$ includes $\kappa u_{32}$ of input rounding
  (EA.19(c)). `theta:1e-k` is read as $[10^{-k}, 10^{-k+1})$, and the realized $\pi - \theta$ of
  `theta:pi-1e-k` is nominal to $\approx\pi u$ (EA.19(a), EA.21).
- `PHASE1.md` §2 item 4: how a stratum's splitmix64 stream is derived is unspecified (counter seeds
  $a_0 + i\mathsf G$ overlap for every pair with $\lvert i - j\rvert < L$; one global stream shifts every later stratum when one
  changes; EA.23), and so is the p99 method (at $n = 64$ the nearest-rank p99 is the maximum;
  EA.11(d)).
- `PHASE2.md` §3 states no rule for $\times$ and $\div$; the two equal forms of the quotient rule
  differ in the last bit on $57\%$ of random inputs (EA.15), so the derivative column of the sweep
  (`0004` item 1) depends on the choice. §3 promises nested `Dual` "tested to second order on
  `sin_cos`, `atan2`" and §8 names no such test; through an exact arm the second derivative loses
  $\theta^{-6}$ (all digits at $\theta \approx 10^{-2}$ for $b$) and the sweep measures the first
  order only (EA.17).
- `NUMERICS.md` §10 and `PHASE4.md` §3 say `Product<SO3, R3>` "fails right-invariance, positively", but `Product` composes
  componentwise (the `Rn`, `Product` row of `PHASE3.md` §0.0; §7 calls `Product<SO3, Rn<3>>` "the tf2-semantics pose" and leaves the law of $a\cdot H$ open) and under that law the geodesic is bi-invariant exactly ($10^{-110}$, GE.4, GE.5(a)); it fails only
  when $(R, t)$ is composed as an SE(3) pose, $a\cdot H = (R_aR_H, R_at_H + t_a)$, by $\lvert\mu_s(\theta)\rvert\lVert t_{H\perp}\rVert$
  (GE.5(b)). §3 has to say which law $a\cdot H$ uses, or its `max_err > 1e-6` test cannot fail.
- `PHASE4.md` §1.2 has the fast twin's branches "through `S::branch` and the coefficients through `coeffs`", but its exact arm needs no
  catalogue coefficient in value (`f64`, given $\hat q_\Delta$: $\le 5.4u$ from $\theta = 10^{-140}$ to $\pi - 10^{-15}$, GE.13(a)). What it needs that §4 and `0004` do
  not provide is the guard for $0/0$ at $\lVert q_{\mathrm v}\rVert = 0$: a range constant of the scalar type and not a coefficient switch (`tf_tree_math`'s
  $10^{-290}$ is `f64`'s, and $10^{-290}$ is not representable in `f32`, whose value and rounding no document chooses or measures), applied with the safe argument
  (GE.13(b)). `tf_tree_math`'s transcendental-free small-angle arm is a speed lever and needs two series outside §4, the slerp weight $\varpi_t(\alpha^2; t)$ and
  $\alpha^2$ from $h = 1 - \cos\alpha$ (a §4 edit and a record each); whether the exact arm alone meets `PHASE4.md` §5.3's bench precondition is a bench question
  (GE.13(a)). No series repairs the derivative through `Dual`: the pose's translation derivative loses $\approx10^2\alpha^{-1}u$ in the assembly of $\varkappa$, $\bar m$ and
  $\varpi_t$, with or without exact derivatives for them (GE.13(c)), so `jacobians_match_dual_geodesic` needs a per-stratum tolerance (CO.14(ii)) at small angle. Nor does any
  document say whether `geodesic` is bit-exact at $t \in \{0,1\}$, which `tf_tree_math`'s `ScLerp` is by early return and the grouped formulas are not at $t = 1$ (GE.13(f)).
  *Answered by [`0054`](../decisions/0054-the-screw-twin-takes-two-arms-in-the-world-frame.md):* below $r$'s second switch the twin takes the
  definition's closed form on §4's $r$, $c$, $a$, $b$ (GE.15(b)), so there is no second series, no guard constant and no $0/0$; $t = 0$ is exact, $t = 1$ exact in
  the rotation and to rounding in the translation, and a consumer that needs both endpoints exact keeps its own early return.
- `NUMERICS.md` §10 lists the right-sided Jacobians only, while `PHASE4.md` §2 takes a side: the left forms are GE.7(d). Its $J_0$ is a
  difference of two $O(1)$ matrices that vanishes like $1 - t$ (relative error $\approx u/(1-t)$, measured $10^{12}u$ at $1 - t = 10^{-12}$);
  the equal form $(1-t)J_l((1-t)d)J_l^{-1}(d)$ has none (GE.7(a)), nor has $J_0^L$'s $(1-t)\mathrm{Ad}_{X_1}J_r((1-t)d)J_r^{-1}(d)\mathrm{Ad}_{X_1}^{-1}$ (GE.7(d)). Which to ship is a §10 edit.
- `PHASE4.md` §4 defines `so3_geodesic`, `se3_geodesic` with `mp.expm`/`mp.logm`; like `so3_log` (above) that is a complex, non-principal
  result from $\theta \approx 3.03$ (SE(3), $t = \tfrac12$: all 6 poses wrong by $1.0$ to $1.9$ at $\theta \ge 3.05$, mpmath 1.3.0 and 1.4.1,
  GE.13), so the `geo:near-pi` stratum needs another reference.
- `PHASE4.md` §1.2's `se3_geodesic_matches_reference` includes near-$\pi$ relative rotations, but within a few $u$ of $\pi$ the two twins may
  pick geodesics $O(1)$ apart (GE.13(d)): it needs a margin, or one shared flip predicate. Its "1e-14" cannot be an absolute bound at
  $\lVert t_0\rVert \sim 10^4$, where the spacing of `f64` is $1.8\times10^{-12}$; `NUMERICS.md` §11's translation floor is the scale to use
  (the exact arm errs by $\le 5.4u$ on that scale given $\hat q_\Delta$, GE.13(a); forming $\Delta$ is not in that figure).
  *Answered by `0054`:* the twins share the flip, since both read the sign of one $q_0^*q_1$, so no band is excluded; the bound is on the translation
  floor, per regime, with extrapolation a row of its own (`sen3_tests::group::se3_geodesic_matches_reference`).
- `PHASE5.md` §1.3 (and `0012`, Context) say the three SE(3) charts "agree to first order at $\delta = 0$ and differ at second order". `Screw` and `Decoupled` do (equal
  $\mathrm D\,\mathrm{ret}(0)$; the gap $R(J_l(\varphi) - I)\rho$, of norm $\ell(\theta)\lVert\rho_\perp\rVert \le \frac\theta2\lVert\rho\rVert$, CH.5(c), CH.7). `WorldTranslation` does not, unless $R = I$: its $\rho$ is the
  world-frame translation, $\mathrm D\Phi^{\mathrm{Dec}\to\mathrm{WT}}(0) = \mathrm{diag}(I, R)$, so for the same $\delta$ the retracted translations differ at first order ($1.13\lVert\delta\rVert$ at $\theta(R) = 1.2$,
  $\delta = [0; 10^{-8}, 0, 0]$, i.e. $2\sin\frac{\theta(R)}2\lVert\rho\rVert$). What is true: it agrees with `Decoupled` exactly after the linear map $\mathrm{diag}(I, R)$ on the translation tangent (CH.4(b)). An LM with
  identity damping then takes the same step in both; a diagonally scaled one does not (CH.7(b)). In storage coordinates the same holds: `WorldTranslation`'s ambient $P_{SE3}$ is $\mathrm{diag}(P, I)$ and `Decoupled`'s is $\mathrm{diag}(P, R)$ (AJ.6(e)), so
  `PHASE6.md` §1's "`Screw` and `Decoupled` alike" is right and does not extend to the third chart.
- `NUMERICS.md` §8 says $\mathrm{sgn}(0) = +1$, `PHASE5.md` §2 has a stratum "$n_z = \pm0$". For $n_z = -0.0$ the reading `n_z >= 0` gives $\varsigma = +1$ and `copysign(1, n_z)` gives $-1$; both bases are valid and differ by the
  reflection $T$ of CH.11(a) ($\lVert B_+ - B_-\rVert_F = 2$), so `s2_retract` and `s2_local` at that stratum differ by $O(1)$ between the readings. §8 should say which.
- `PHASE5.md` §2 gives the reference of `s2_local` as the minimal rotation "via `mp.expm`/`mp.logm`". Like `so3_log` above, `mp.logm` of it is complex and non-principal from $\alpha \approx 3.03$ ($20$ of $20$ samples per angle at
  $3.03, 3.04, 3.05, 3.1, 3.14$; exact at $3.02$; 50 digits, mpmath 1.3.0 and 1.4.1; CH.9), so the `s2:near-antipode` stratum needs a reference that does not use `mp.logm` (the geometric $\mathrm{SO}(3)$ inverse of CH.9's definition agrees with
  the closed form to $2.6\times10^{-104}$).
- `NUMERICS.md` §12 has no row for the unit-norm domain of `S2` (`n`, and `m` in `local`), while the quaternion has one (`from_wxyz_unchecked`, $2^{-40}$ in `f64`); §8 and `PHASE5.md` §2 do not say whether `S2` renormalises. With $\lVert n\rVert = 1 + \varepsilon$
  the basis of `Chart::at` keeps $B^\top B = I$ but has $B^\top n = O(\varepsilon)$ (the $H$ form), or loses both by $2\varepsilon$ (the closed form of CH.8(d), which spends $\nu^\top\nu = 2(1 + \lvert n_z\rvert)$) (CH.8(e)). A tolerance, or a renormalisation, is a `NUMERICS.md` §8, §12 edit and needs a record.
- `NUMERICS.md` §9, §12, §2.4 and `PHASE5.md` §3, on Sim(3) (`sim3.md`). (i) §12's row "`Sim3` | $\sigma$ finite | —" promises no finite result: beyond $\ln\mathrm{MAX}$ ($709.78$ `f64`, $88.72$ `f32`) $s^{\pm1}$ and $X^{-1}$ overflow, and the closed forms fail well before, silently: $v_0\Lambda$ overflows at $242.1$ ($33.1$ `f32`) and zeroes $i_2$, $v_1$ and $v_2$ overflow at $703.25$ ($84.3$); a scaled arrangement holds to $\ln\mathrm{MAX}$ (SM.7(g)); there is no row for `log`, `jr_inv`, `jl_inv` ($\theta < 2\pi$, $\lvert\sigma\rvert < \ln\mathrm{MAX}$) or the derivative of `log` ($\theta(X) < \pi$); the strata stop at $\lvert\sigma\rvert = 3$ (SM.13).
  (ii) `PHASE5.md` §3 takes the reference of `sim3_log` from `mp.logm` of the $4\times4$ matrix, which is complex and non-principal from an onset that falls with $\lvert\sigma\rvert$: $(3.02, 3.03]$ at $\sigma = 0$, $(3.00, 3.01]$ at $1$, $(2.96, 2.97]$ at $2$, $(2.85, 2.86]$ at $3$, $(2.88, 2.89]$ at $-3$, $(2.19, 2.20]$ at $6$ (SM.8: 8 axes, 40 digits; the rows tested at 80 digits agree), so most `theta:pi-1e-k` strata of the grid $\lvert\sigma\rvert \le 3$ would be wrong; the geometric $\mathrm{SO}(3)$ logarithm with a $\mathsf V$ solve agrees.
  (iii) The ids of §3 are `sim3_exp`, `sim3_log`, `sim3_jr`, `sim3_jr_inv`, `sim3_ad`: none for `jl`, `jl_inv`, the action Jacobians (`API.md` §3 lists `act_jacobians` for Sim(3); §2.4 has no Sim(3) rows, SM.12) or the scalar coefficients; a `sim3_jr` stratum at $\sigma = \theta = 10^{-12}$ cannot pass with closed forms alone ($\lesssim4.5u/r_\Delta$, $2\times10^{12}u$ measured, SM.7(c)). `0005` says "the scale column breaks the dual structure": the diagonal blocks $J_l(\varphi)$ and $\mathsf V$ already differ (SM.9), and the $(3,3,1)$ pattern is closed under product and inverse ($117$ against $343$ multiplications); §14 has no `jr_inv` row for `Sim3` (block algebra as the fast path, dense inverse as the twin: SM.10(c)).
  (iv) `Real` (`PHASE2.md` §2) has neither `exp` nor `exp_m1`, and every Sim(3) routine needs both ($s = e^\sigma$, $v_0 = \mathrm{expm1}\,\sigma/\sigma$, $g_\ell(\sigma)$, and $e^{-\sigma}$, $\mathrm{expm1}(-\sigma)$, $e^{-\sigma/2}$ of the scaled arrangement: SM.4, SM.7(g)): new public items, so a draft record with a `Dual` rule for each (`PHASE2.md` §3), a `libm` route (D16; the `f64` constants of `sim3.md` are glibc's), `real_*` strata (`PHASE2.md` §8) and the answers of `API.md` §6 —
  which is [`0022`](../decisions/0022-real-owes-acos-and-cos.md)'s pattern, already `ready` for `acos` and `cos`, so this is a second record of the same shape and not a new mechanism. §2's listed set (`sqrt`, `cbrt`,
  `sin_cos`, `atan2`, `abs`, `copysign`) has not been updated for `0022` either, so `Real` is short two pairs, not one; the single fix is a `PHASE2.md` §2 edit reconciling its method set with every owed transcendental.
- `NUMERICS.md` §9 lists "the two-dimensional series at the joint limit" as owed. `Exp`'s value and `Log` do not need it (the closed forms are $O(u)$ per consumer, SM.7(b)); the Jacobian blocks do ($u/r_\Delta$, SM.7(c)), by the series of SM.6(b) (one branch variable, $\Delta$, which also covers $\Delta = 0$, and one recurrence) or by the exact split of SM.6(d), which is not mask-free: its weights are $0/0$ at $\Delta = 0$ (a select beside the $\Delta$ mask) and it adds three one-variable rows ($g_\ell$, $v_1^\circ$, $v_2^\circ$, each with a $\lvert\sigma\rvert$ branch) to §4's (SM.7(h)). `0004` item 1 has to say which objective a sweep over these coefficients minimizes: the coefficient's own error (no digit of $v_2$ below $r_\Delta \approx 10^{-8}$ in closed form) or the consumer's ($u/r_\Delta$ in the Jacobian blocks only). A `Dual` derivative of $v_j$ loses one more power of $\Delta$ (SM.7(f), not measured). SM.7(d)(ii)'s block table was measured at a recurrence depth that leaves $\partial_zv_1$ and $\partial_zv_2$ two weighted degrees short ($28u$ and $401u$ against $2.1u$ and $2.5u$ at $K = 20$, $r_\Delta = 1.6$); the per-sequence depth is stated there and the table is owed at it.

## References

- **[Hall]** B. C. Hall, *Lie Groups, Lie Algebras, and Representations*, 2nd ed., Springer GTM 222,
  2015 (matrix groups; the differential of the exponential map).
- **[Wilcox]** R. M. Wilcox, "Exponential operators and parameter differentiation in quantum
  physics", J. Math. Phys. 8(4), 1967.
- **[AMH]** A. H. Al-Mohy, N. J. Higham, "Computing the Fréchet derivative of the matrix
  exponential", SIAM J. Matrix Anal. Appl. 30(4), 2009.
- **[Chirikjian]** G. S. Chirikjian, *Stochastic Models, Information Theory, and Lie Groups*, vol. 2,
  Birkhäuser, 2012.
- **[Barfoot]**, **[BF14]**, **[Solà]**, **[Shepperd]**: Barfoot 2017; Barfoot & Furgale 2014; Solà,
  Deray, Atchuthan 2018; Shepperd 1978 (full entries in `NUMERICS.md` §13).
- **[Horn]** B. K. P. Horn, "Closed-form solution of absolute orientation using unit quaternions",
  J. Opt. Soc. Am. A 4(4), 1987. **[Bar-Itzhack]** I. Y. Bar-Itzhack, "New method for extracting the
  quaternion from a rotation matrix", J. Guidance, Control, and Dynamics 23(6), 2000.
- **[Higham]** N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM, 2002
  (ch. 1: forward and backward error, conditioning; ch. 2: Sterbenz's lemma; ch. 3: the $\gamma_n$
  bound for inner products; §5.1: Horner's method; ch. 8 and 10: triangular solves and Cholesky, Thm 8.5, 10.3, the bounds of `NUMERICS.md` §15.4).
- **[DLMF]** *NIST Digital Library of Mathematical Functions*, https://dlmf.nist.gov, §4.19 and §4.22
  (the series and the partial fractions of $\cot$), §24.2 (Bernoulli numbers).
- **[SLF]** G. L. Steele Jr., D. Lea, R. Flood, "Fast splittable pseudorandom number generators",
  OOPSLA 2014. **[Vigna]** S. Vigna, `splitmix64.c`, public domain,
  https://prng.di.unimi.it/splitmix64.c (the reference outputs of EA.23).
- **[Milnor]** J. Milnor, "Curvatures of left invariant metrics on Lie groups", Adv. Math. 21, 1976 (bi-invariant metrics; the geodesics of $\mathrm{SO}(3)$). **[Arnold]** V. I. Arnold, "Sur la géométrie différentielle des groupes de Lie de dimension infinie et ses applications à l'hydrodynamique des fluides parfaits", Ann. Inst. Fourier 16, 1966 (the Euler–Poincaré equation of a left-invariant metric). **[Kavan]** L. Kavan, S. Collins, C. O'Sullivan, J. Žára, "Dual quaternions for rigid transformation blending", TCD-CS-2006-46, 2006 (the screw parameters of a dual quaternion; ScLERP).
- **[MilnorHB]** J. Milnor, "Analytic proofs of the 'hairy ball theorem' and the Brouwer fixed point theorem", Amer. Math. Monthly 85(7), 1978 (no continuous unit tangent field on $S^2$; CH.11(b)).
- **[Ceres]** S. Agarwal, K. Mierle et al., *Ceres Solver*, `include/ceres/manifold.h` and `internal/ceres/manifold.cc` (`QuaternionManifold`, `EigenQuaternionManifold`), read at commit `e17a9b41066a` (2026-09-26); the `pyceres` 2.6 wheel. **[GTSAM]** F. Dellaert et al., *GTSAM*, `gtsam/base/Lie.h`, `gtsam/geometry/Rot3.h`, `Pose3.cpp`, `Quaternion.h`, read at develop `c786d78c5a4e` (2026-09-29); the 4.3.0 wheel (the tools behind AJ.7, AJ.8).
- **[mpmath]** F. Johansson et al., *mpmath: a Python library for arbitrary-precision
  floating-point arithmetic* (the checking tool; version stated in each `Checked:` line).
