# Sim(3): Exp, Log, adjoints, Jacobians and the joint limit

> Non-normative companion to [`NUMERICS.md`](../NUMERICS.md) §9 (Sim(3), **owed**), with §1, §2.3, §2.4, §4, §12 and
> [`PHASE5.md`](../PHASE5.md) §3. **`NUMERICS.md` wins on any conflict; a disagreement is an open item in the
> [maths index](./index.md)**, which holds the notation and the `Checked:` convention. **Every result below that §9 does not
> already state is Proposed for `NUMERICS.md` §9: to be adopted by a record (`PHASE5.md` §3 does not start before that).** "Proposed" in a heading means exactly this.
> The composition of SM.2 and the $\mathrm{Ad}$ of SM.9 restate §9 and are only checked; the inverse and the action are proposed. General algebra:
> [`lie-groups.md`](./lie-groups.md); SO(3) and its coefficients: [`so3.md`](./so3.md), [`coefficients.md`](./coefficients.md);
> the same derivation for $\mathrm{SE}_N(3)$: [`se3.md`](./se3.md). Nothing here is copied from Eade 2014 or Strasdat 2012.

On this page $\tau = [\varphi;\rho;\sigma]$ (one translation block and the scale), $W = [\varphi]_\times$, $\Xi = [\rho]_\times$,
$\mu = W\rho$, $\nu = W\mu$, $z = \theta^2$, $\Delta = \sigma^2 + \theta^2$, $r_\Delta = \sqrt\Delta$, $s = e^\sigma$,
$M = W + \sigma I$, $x = \sigma + i\theta$; $\psi = [\alpha;\eta;\varsigma]$ is a second tangent (SM.9 only). $g_\ell(x) = \sum_{k\ge0}x^k/(k+\ell)!$ is entire
($g_1$ is LG.7's $g$, and $\Gamma_\ell(\varphi) = g_\ell(W)$ is GG.2's $\Gamma_m$); $\chi = \mathrm{Re}\,g_1(x) = v_0 - zv_2$ and $\Lambda = \lvert g_1(x)\rvert^2 = \chi^2 + zv_1^2$ (SM.5(c)). $a, b, d$ are `NUMERICS.md` §4's, with
$\mathrm{sinc}\,\theta = \sin\theta/\theta = 1 - zb$ and $a = \frac12 - zd$ (CO.4(a)). The block `NUMERICS.md` §9 calls $W(\varphi,\sigma)$ is written $\mathsf V$, because $W = \varphi^\wedge$
throughout, with the same arguments; $\mathsf Q(\rho, \varphi, \sigma)$ keeps §5.3's order for $Q(\rho, \varphi)$, $\sigma$ last; $s$ is the scale $e^\sigma$ (not LG.16's $s(\theta)$). In `Checked:` lines "max" is the largest absolute entry and
"norm-wise" the Frobenius error over $\max(1, \lVert\cdot\rVert_F)$ (`NUMERICS.md` §11); `f64` runs are CPython floats with glibc 2.39
(`exp`, `expm1`, `sin`, `cos`) and numpy 2.5.3 for the $7\times7$ assembly, referenced to mpmath 1.3.0 (sympy 1.14.0 for the exact
identities) from definitions: `mp.expm`, the defining series reduced by SM.3, the dense series of $\mathrm{ad}$, `mp.inverse`, central
differences of the operation itself. The switches inside those `f64` arms ($\theta = 1$ for $b, d$; $\lvert\sigma\rvert = 0.5$ for the
functions of $\sigma$; $r_\Delta \approx 1$ between closed form and series; $\Delta = 1$ between the two forms of $v_1$; the number $K$ of series terms)
are illustrations, not generated switches (0004). The `f64` constants are glibc's: the `libm` crate (D16) errs by more (`coefficients.md`) and its `exp`, `expm1` are not measured,
so they may move by an ulp or two; the `sim3_*` strata re-measure them.

## Results

| Label | Result | `NUMERICS.md` | Implemented by |
|---|---|---|---|
| SM.1–SM.2 | $\mathrm{Sim}(3)$ as $4\times4$ matrices; $XY$, $X^{-1} = (R^\top, -s^{-1}R^\top t, -\sigma)$, $Xp = sRp + t$ | §9 (composition stated, checked; inverse, action proposed) | `Sim3`, `Mul`, `inverse`, `act` |
| SM.3–SM.4 | $M^k = \sigma^kI + \xi_kW + \zeta_kW^2$; $\mathrm{Exp}\,\tau = (\mathrm{Exp}\,\varphi,\ \mathsf V\rho,\ \sigma)$, $\mathsf V = g_1(\sigma I + W) = \int_0^1e^{t\sigma}e^{tW}\mathrm dt = v_0I + v_1W + v_2W^2$ | §9 (proposed) | `Sim3::exp` |
| SM.5 | closed forms of $v_0, v_1, v_2$ (denominator $\Delta$); spectrum, $\det\mathsf V$, positivity; the limits $\sigma = 0$ ($J_l(\varphi)$) and $\theta = 0$ | §9, §4 (proposed) | `Sim3::exp`, a new coefficient group |
| SM.6 | the series in $(\sigma, z)$: general terms, total-degree form, recurrences, tail bound, rays; the exact split into one-variable functions | §9 (proposed) | series arm, `xtask thresholds` |
| SM.7 | rounding: $u/\Delta$ per coefficient, $O(u)$ for values, $u/r_\Delta$ for the Jacobian blocks; what removes it; crossing radii; large $\sigma$ (overflow onsets, the scaled arrangement); the evaluation in one table | §9, §4, §11 (proposed) | as SM.6 |
| SM.8 | $\mathrm{Log}$ and $\mathsf V^{-1}$ by algebra in $\{I, W, W^2\}$; the two preimages at $\theta = \pi$ | §9 (proposed) | `Sim3::log` |
| SM.9 | $\mathrm{Ad}$, $\mathrm{ad}$, $\mathrm{Ad}^{-1}$ rotation-first; why no dual structure | §9 ($\mathrm{Ad}$ stated), §5.2 analogue | `adjoint`, `ad` |
| SM.10 | $J_l$, $J_r$: the $(3,3,1)$ block structure; $\mathsf Q$ and the $\Gamma_2$ column from $y = \mathsf V\rho$; $\det$; inverse by blocks | §9, §12 (proposed) | `jr`, `jl`, `jr_inv`, `jl_inv` |
| SM.11 | the scalars of the blocks: $\Gamma_\ell(M)$, $w_0, w_1, w_2$, $\partial_zv_1, \partial_zv_2$; the joint-limit form of $e$ | §9, §4 (proposed) | a new coefficient group |
| SM.12 | action Jacobians, both sides | §2.4 (proposed rows) | `act_jacobians` |
| SM.13 | domains: $\theta$, $\sigma$, overflow, $\Delta = 0$ | §12 (proposed row) | — |

## 1. The group

**Definition SM.1.** $\mathrm{Sim}(3)$ is the group of matrices $X = \begin{bmatrix}sR & t\\ 0 & 1\end{bmatrix}$, $s = e^\sigma$,
$R \in \mathrm{SO}(3)$, $t \in \mathbb R^3$, stored as $(q, t, \sigma)$ (`NUMERICS.md` §9). It is closed in $\mathrm{GL}_4(\mathbb R)$; its
algebra is $\{\begin{bmatrix}\Omega + \sigma I & \rho\\ 0 & 0\end{bmatrix}\}$, and $\tau^\wedge = \begin{bmatrix}W + \sigma I & \rho\\ 0 & 0\end{bmatrix}$
is the rotation-first hat of §9 ($\supseteq$ by SM.4, $\subseteq$ by differentiating the block shape of a curve in the group at $0$, as LG.1).
$\mathrm{Exp}\,\tau = \exp(\tau^\wedge)$; the matrix logarithm, $\mathrm{Log}$ and the sides are those of `NUMERICS.md` §1.

**Proposition SM.2 (composition stated in §9; inverse and action proposed).** (a) $XY = (R_XR_Y,\ s_XR_Xt_Y + t_X,\ \sigma_X + \sigma_Y)$ and
$X^{-1} = (R^\top,\ -s^{-1}R^\top t,\ -\sigma)$; (b) $Xp = sRp + t$; (c) LG.4(a)–(b), LG.5, LG.8 and LG.9(a)–(b) hold verbatim, their proofs using only a matrix
group with a linear hat; LG.4(c), LG.9(c) and LG.10 do not (SM.9–SM.10); (d) LG.11–LG.15 hold, hence every row of §2.3 (both sides: $\oplus$, $\ominus$, $\mathrm{Log}$, the product, the inverse) with $\mathrm{Ad}$, $J_l$, $J_r$ of SM.9–SM.10, SM.10(c) in place of LG.9(c) (the invertibility of $J$ that LG.13(d) uses) and SM.8 in place of LG.2(c) ($\mathrm{Log}\circ\mathrm{Exp} = \mathrm{id}$ on $\theta < \pi$).

*Proof.* Block multiplication, $s_Xs_Y = e^{\sigma_X + \sigma_Y}$, and $\begin{bmatrix}sR & t\\ 0 & 1\end{bmatrix}\begin{bmatrix}s^{-1}R^\top & -s^{-1}R^\top t\\ 0 & 1\end{bmatrix} = I$. $\square$

The scale is additive, so $X^{-1}$ negates $\sigma$ exactly; only $Xp$, the composition and $X^{-1}$'s translation form $e^{\pm\sigma}$.

**Checked:** mpmath 1.3.0, 110 digits, 30 random pairs (rotation vector norm $\sim U[0,3]$, $\sigma \sim U[-3,3]$, $t \sim N(0,2^2)$): the $4\times4$
product and `mp.inverse` against the block formulas, max $1.9\times10^{-111}$ and $1.5\times10^{-111}$. (d) has no separate check: SM.9, SM.10 and SM.12 test the pieces its rows use (LG.4(b), LG.13(c)). **Permanent:** planned, proptest
`group_axioms_*` (`PHASE3.md` §9).

## 2. Exp

**Lemma SM.3 (powers of $M$; supports the proposed SM.4–SM.6).** $M^k = \sigma^kI + \xi_kW + \zeta_kW^2$ with $\xi_0 = \zeta_0 = 0$ and
$$
\xi_{k+1} = \sigma\xi_k + \sigma^k - z\zeta_k,\qquad \zeta_{k+1} = \sigma\zeta_k + \xi_k,
$$
i.e. $\xi_k = \sum_m\binom k{2m+1}(-z)^m\sigma^{k-2m-1} = \mathrm{Im}\,x^k/\theta$ and $\zeta_k = \sum_m\binom k{2m+2}(-z)^m\sigma^{k-2m-2} = (\sigma^k - \mathrm{Re}\,x^k)/z$:
polynomials in $(\sigma, z)$ of weighted degree $k-1$ and $k-2$ (weights $1$ and $2$). $\lvert\xi_k\rvert \le kr_\Delta^{k-1}$ and $\lvert\zeta_k\rvert \le \binom k2r_\Delta^{k-2}$, with equality at $\theta = 0$.

*Proof.* $\sigma I$ commutes with $W$; the binomial theorem with $W^{2m+1} = (-z)^mW$, $W^{2m+2} = (-z)^mW^2$ (SO.3) gives the sums, and
$(\sigma + i\theta)^k$ expanded gives the real and imaginary parts. The recurrence is $M^{k+1} = MM^k$ with $W^3 = -zW$. With $p_k = \mathrm{Re}\,x^k = \sigma^k - z\zeta_k$,
$\lvert p_k\rvert \le r_\Delta^k$ and $\xi_{k+1} = \sigma\xi_k + p_k$, so $\lvert\xi_{k+1}\rvert \le r_\Delta\lvert\xi_k\rvert + r_\Delta^k$; then
$\lvert\zeta_{k+1}\rvert \le r_\Delta\lvert\zeta_k\rvert + \lvert\xi_k\rvert$: induction from $\xi_1 = 1$, $\zeta_2 = 1$. $\square$

**Checked:** sympy 1.14.0: the recurrences equal the binomial sums for $k \le 14$; $M^k = \sigma^kI + \xi_kW + \zeta_kW^2$ with a symbolic $\varphi$ for $k \le 7$;
the two bounds on 4000 sampled $(\sigma, z)$, $k \le 11$: ratios $\le 1$, equal to $1$ on the $\sigma$ axis. **Permanent:** none specified (used by SM.6's series).

**Proposition SM.4 (Exp; proposed).** For every $\tau$,
$$
\mathrm{Exp}\,\tau = \begin{bmatrix}e^\sigma e^W & \mathsf V\rho\\ 0 & 1\end{bmatrix} = (\mathrm{Exp}\,\varphi,\ \mathsf V(\varphi,\sigma)\rho,\ \sigma),\qquad
\mathsf V = \sum_{k\ge0}\frac{M^k}{(k+1)!} = g_1(M) = \int_0^1e^{t\sigma}e^{tW}\mathrm dt = v_0I + v_1W + v_2W^2,
$$
$$
v_0 = \int_0^1e^{t\sigma}\mathrm dt = g_1(\sigma),\quad
v_1 = \int_0^1e^{t\sigma}\frac{\sin t\theta}\theta\,\mathrm dt,\quad
v_2 = \int_0^1e^{t\sigma}\frac{1 - \cos t\theta}{\theta^2}\,\mathrm dt,
$$
entire in $(\sigma, z)$. The rotation part is `NUMERICS.md` §3.1's quaternion, the scale is $e^\sigma$ exactly, and $\mathrm{Exp}$ costs one SO(3) $\mathrm{Exp}$, one
$\mathsf V$ and one matrix-vector product. $\mathsf V = J_l(\varphi)$ at $\sigma = 0$ and $\mathsf V = g_1(\sigma)I$ at $\varphi = 0$. It needs $e^\sigma$ and $\mathrm{expm1}\,\sigma$ (SM.5(b)), which `Real` (`PHASE2.md` §2) has neither of: a new public item (index, open items).

*Proof.* $(\tau^\wedge)^k = \begin{bmatrix}M^k & M^{k-1}\rho\\ 0 & 0\end{bmatrix}$ for $k \ge 1$ (induction), and summing $\sum(\tau^\wedge)^k/k!$ gives the top-right block
$\sum_{k\ge1}M^{k-1}\rho/k!$; $e^M = e^\sigma e^W$ since $\sigma I$ commutes with $W$. $g_1(M) = \int_0^1e^{tM}\mathrm dt$ termwise, $e^{tM} = e^{t\sigma}e^{tW}$ and
$e^{tW} = I + \frac{\sin t\theta}\theta W + \frac{1 - \cos t\theta}{\theta^2}W^2$ (SO.4). $\sigma = 0$: $g_1(W) = J_l(\varphi)$ (SO.7). $\square$

**Checked:** mpmath 1.3.0, 110 digits: `mp.expm` of the $4\times4$ hat matrix against the blocks, 111 cases (a $9\times9$ grid $\theta \in \{10^{-3}, 0.5, 1, 2, 3, \pi, 4, 6, 20\}$ ×
$\sigma \in \{0, 10^{-6}, -0.3, 1, -2, 3, -3, 10, -30\}$ and 30 random, random axes, $\rho \sim N(0,2^2)$): $\mathsf V\rho$ $5.1\times10^{-110}$ (relative to $\max(1, \lVert\cdot\rVert)$),
rotation block $1.7\times10^{-111}$; large $\lvert\sigma\rvert \in \{100, 300, 700\}$ (both signs) and $\pm10^{-30}$, $\theta \in \{0.5, 3, \pi\}$, 120 digits: $5.7\times10^{-119}$
relative to the column; the integrals by `mp.quad` (40 digits) $1.0\times10^{-35}$ (the quadrature). **Permanent:** planned, corpus `sim3_exp` (`PHASE5.md` §3), proptest `exp_log_roundtrip_*`.

## 3. The coefficients in closed form

**Proposition SM.5 (closed forms; proposed).**

- (a) In the basis $I, W, W^2$, $\mathsf VM = e^M - I$ reads $\sigma v_0 = e^\sigma - 1$, $\ v_0 + \sigma v_1 - zv_2 = e^\sigma\,\mathrm{sinc}\,\theta$, $\ v_1 + \sigma v_2 = e^\sigma a$.
- (b) Hence, with $\mathrm{expm1}\,\sigma = e^\sigma - 1$ (and $v_0 = 1$ at $\sigma = 0$), for $\Delta > 0$,
$$
v_0 = \frac{\mathrm{expm1}\,\sigma}\sigma,\quad
v_2 = \frac{v_0 + \sigma e^\sigma a - e^\sigma\,\mathrm{sinc}\,\theta}\Delta,\quad
v_1 = e^\sigma a - \sigma v_2 = \frac{\sigma e^\sigma\,\mathrm{sinc}\,\theta - \mathrm{expm1}\,\sigma + e^\sigma za}\Delta .
$$
  Call the two expressions of $v_1$ *form a* ($e^\sigma a - \sigma v_2$) and *form b* (the quotient).
- (c) $M$ has eigenvalues $\sigma$, $x$, $\bar x$; $\mathsf V$ is normal with eigenvalues $g_1(\sigma)$, $g_1(x)$, $\overline{g_1(x)}$, $v_1 = \mathrm{Im}\,g_1(x)/\theta$ and $\chi = v_0 - zv_2 = \mathrm{Re}\,g_1(x)$.
$\det\mathsf V = g_1(\sigma)\Lambda$, $\Lambda = \lvert g_1(x)\rvert^2 = \chi^2 + zv_1^2 = \big[\mathrm{expm1}(\sigma)^2 + 4e^\sigma\sin^2\tfrac\theta2\big]/\Delta$, a sum of two non-negative terms. $\mathsf V$ is singular iff $\sigma = 0$ and $\theta \in 2\pi\mathbb Z_{>0}$; its singular values are $g_1(\sigma), \sqrt\Lambda, \sqrt\Lambda$.
- (d) $v_0 > 0$, $v_2 > 0$ for all $(\sigma, \theta)$, and $v_1 > 0$ for $\theta \le \pi$.

*Proof.* (a) $\mathsf VM = e^M - I$ from $g_1(M)M = e^M - I$; multiply out with $W^3 = -zW$:
$(v_0 + v_1W + v_2W^2)(\sigma + W) = \sigma v_0 + (v_0 + \sigma v_1 - zv_2)W + (v_1 + \sigma v_2)W^2$, and $e^M - I = (e^\sigma - 1)I + e^\sigma\mathrm{sinc}\,\theta\,W + e^\sigma aW^2$ (SM.4).
$I, W, W^2$ are independent for $\varphi \ne 0$, and the equations extend to $\varphi = 0$ by continuity. (b) Insert the third equation into the second: $\Delta v_2 = v_0 + \sigma e^\sigma a - e^\sigma\mathrm{sinc}\,\theta$.
The second form of $v_1$ is $e^\sigma a - \sigma v_2$ with $\sigma v_0 = e^\sigma - 1$ and $\Delta - \sigma^2 = z$. (c) $\mathsf V$ is a real polynomial in the skew matrix $W$, so $\mathsf V^\top = v_0I - v_1W + v_2W^2$ commutes with it: $\mathsf V$ is normal. $W$ has spectrum $\{0, \pm i\theta\}$; on the plane $\perp\varphi$ it acts as $i\theta$, so
$v_0I + v_1W + v_2W^2$ acts there as $v_0 + i\theta v_1 - zv_2 = g_1(x)$, and on the axis as $v_0 = g_1(\sigma)$. $\lvert e^x - 1\rvert^2 = (e^\sigma - 1)^2 + 2e^\sigma(1 - \cos\theta)$, $\mathrm{Im}\,g_1(x) = \theta v_1$ gives $\Lambda = \chi^2 + zv_1^2$, and $g_1(x) = 0$ iff $e^x = 1$, $x \ne 0$.
(d) The integrands of SM.4 are $\ge 0$ ($e^{t\sigma} > 0$, $1 - \cos \ge 0$, $\sin t\theta \ge 0$ for $t\theta \le \pi$) and not identically $0$. $\square$

**Checked:** mpmath 1.3.0. (a), (b): SM.4's 111 cases (the blocks of `mp.expm` are the closed forms); (c) 70 digits, $\det\mathsf V$ against $g_1(\sigma)\Lambda$ over $84$ $(\sigma, \theta)$ ($8\times8$ grid with $\theta$ to $6$, $\sigma$ to $\pm30$, and 20 random), relative
$2.8\times10^{-59}$; 100 digits, $\Lambda$'s real form against $\chi^2 + zv_1^2$, $\chi = v_0 - zv_2$, at 55 points, $1.5\times10^{-100}$; 60 digits, normality and singular values (`mp.svd_r`) at 20 points, $10^{-59}$; (d) no violation at those 20 points.
**Permanent:** planned, corpus `sim3_exp`; none for (c), (d).

**Proposition SM.6 (series; proposed).**

- (a) *Double series.* $v_0 = \sum_i\dfrac{\sigma^i}{(i+1)!}$,
$$
v_1 = \sum_{m,i\ge0}\frac{(-z)^m\sigma^i}{(2m+1)!\,i!\,(2m+i+2)},\qquad
v_2 = \sum_{m,i\ge0}\frac{(-z)^m\sigma^i}{(2m+2)!\,i!\,(2m+i+3)},
$$
absolutely convergent for every $(\sigma, z)$; at $\sigma = 0$ they are $a = \sigma_2$ and $b = \sigma_3$ of CO.1.
- (b) *Total degree.* $v_1 = \sum_{n\ge0}\xi_{n+1}/(n+2)!$, $v_2 = \sum_{n\ge0}\zeta_{n+2}/(n+3)!$: homogeneous terms of weighted degree $n$ from the recurrences of SM.3 ($O(K)$ work, constants $1/(k+1)!$ only):
$$
\begin{aligned}
v_1 &= \tfrac12 + \tfrac\sigma3 + \tfrac{3\sigma^2 - z}{24} + \tfrac{\sigma(\sigma^2 - z)}{30} + \tfrac{5\sigma^4 - 10\sigma^2z + z^2}{720} + \cdots,\\
v_2 &= \tfrac16 + \tfrac\sigma8 + \tfrac{6\sigma^2 - z}{120} + \tfrac{\sigma(2\sigma^2 - z)}{144} + \tfrac{15\sigma^4 - 15\sigma^2z + z^2}{5040} + \cdots,
\end{aligned}
$$
and $v_0 = 1 + \frac\sigma2 + \frac{\sigma^2}6 + \frac{\sigma^3}{24} + \cdots$. Keeping degrees $n < K$, the truncation errors are at most $\sum_{k>K}kr_\Delta^{k-1}/(k+1)!$ ($v_1$) and $\sum_{k>K+1}\binom k2r_\Delta^{k-2}/(k+1)!$ ($v_2$): **a function of $r_\Delta = \sqrt\Delta$ alone, so one mask on $\Delta$ serves every direction**.
- (c) *Rays.* On $\sigma = r_\Delta\cos\lambda$, $\theta = r_\Delta\sin\lambda$ ($\lambda$ the angle of the ray): $\xi_k = r_\Delta^{k-1}U_{k-1}(\cos\lambda)$, $\zeta_k = r_\Delta^{k-2}(\cos^k\lambda - T_k(\cos\lambda))/\sin^2\lambda$ (Chebyshev). On $\sigma = \kappa\theta$,
$v_1 = \sum_n\mathrm{Im}(\kappa + i)^{n+1}\theta^n/(n+2)!$ and $v_2 = \sum_n[\kappa^{n+2} - \mathrm{Re}(\kappa + i)^{n+2}]\theta^n/(n+3)!$:

| $n$ | $0$ | $1$ | $2$ | $3$ | $4$ |
|---|---|---|---|---|---|
| $v_1$, coefficient of $\theta^n$ | $\frac12$ | $\frac\kappa3$ | $\frac{3\kappa^2 - 1}{24}$ | $\frac{\kappa(\kappa^2 - 1)}{30}$ | $\frac{5\kappa^4 - 10\kappa^2 + 1}{720}$ |
| $v_2$, coefficient of $\theta^n$ | $\frac16$ | $\frac\kappa8$ | $\frac{6\kappa^2 - 1}{120}$ | $\frac{\kappa(2\kappa^2 - 1)}{144}$ | $\frac{15\kappa^4 - 15\kappa^2 + 1}{5040}$ |

Every ray tends to $(\frac12, \frac16)$, uniformly: $v_1$, $v_2$ are entire in $(\sigma, z)$ (a) and $\lvert v_j - v_j(0)\rvert \le Cr_\Delta$ with $C$ independent of the ray (SM.3's bounds). What depends on the ray is the first correction, and the coefficients above are polynomials in $\kappa$: the expansion in $\theta$ at fixed $\kappa$ is not uniform in $\kappa$ ($\kappa \to \infty$ is the $\sigma$ axis), the uniform statement is (b), a function of $r_\Delta$ alone. The singularity of the closed forms at $\Delta = 0$ is removable.
- (d) *Exact split.* With $v_j^\circ(\sigma) = v_j(\sigma, 0)$,
$$
v_1 = \frac{\sigma^2}\Delta v_1^\circ + \frac z\Delta e^\sigma(a - \sigma b),\qquad
v_2 = \frac{\sigma^2}\Delta v_2^\circ + \frac z\Delta e^\sigma(b - \sigma d),
$$
$v_1^\circ = g_1 - g_2$, $v_2^\circ = \frac12g_1 - g_2 + g_3$ (functions of $\sigma$ alone, $g_\mu(\sigma) = [e^\sigma - \sum_{k<\mu}\sigma^k/k!]/\sigma^\mu$); the weights $\sigma^2/\Delta$, $z/\Delta$ lie in $[0,1]$ and sum to $1$, and are $0/0$ at $\Delta = 0$: the split needs the same select there as the closed forms (SM.7(h)).

*Proof.* (a) Expand $e^{t\sigma}$, $\sin t\theta/\theta = \sum(-z)^mt^{2m+1}/(2m+1)!$ and $(1 - \cos t\theta)/z = \sum(-z)^mt^{2m+2}/(2m+2)!$ in SM.4 and integrate $t^{2m+i+1}$, $t^{2m+i+2}$. (b) SM.3 in $\sum M^k/(k+1)!$; the tails by SM.3's bounds. (c) $x^k = r_\Delta^ke^{ik\lambda}$ and
$\mathrm{Im}(\kappa + i)^k = \sin k\lambda/\sin^k\lambda$ for $\kappa = \cot\lambda$. (d) With $\mathrm{sinc}\,\theta = 1 - zb$ and $a = \frac12 - zd$, SM.5(b) gives
$\Delta v_2 = [v_0 + \frac\sigma2e^\sigma - e^\sigma] + ze^\sigma(b - \sigma d)$ and $\Delta v_1 = [\sigma e^\sigma - \mathrm{expm1}\,\sigma] + ze^\sigma(a - \sigma b)$; the brackets are $\Delta v_j$ at $z = 0$, i.e. $\sigma^2v_j^\circ$.
$v_1^\circ = \int_0^1te^{t\sigma}\mathrm dt = g_1 - g_2$ because $\frac1{i!(i+2)} = \frac1{(i+1)!} - \frac1{(i+2)!}$, and $v_2^\circ = \int_0^1\frac{t^2}2e^{t\sigma}\mathrm dt$ with $\frac{(i+1)(i+2)}{2(i+3)!} = \frac1{2(i+1)!} - \frac1{(i+2)!} + \frac1{(i+3)!}$. $\square$

**Checked:** sympy 1.14.0, exact rationals: the closed-form numerators $\Delta v_1$, $\Delta v_2$ (SM.5(b)) expand to $\Delta\times$ the defining series through total degree $8$; the general terms of (a) and $\xi_k/(k+1)!$, $\zeta_k/(k+1)!$ agree for $k \le 8$ ($\zeta$: $9$); the table of (c) and the leading terms of (b) are those polynomials.
mpmath 1.3.0, 60 digits: the tails of (b) for $K \in \{8, 12, 16, 20\}$, $r_\Delta \in \{0.1, 0.4, 0.9, 1.5\}$, max over 400 directions, are within $3\%$ of the stated bounds (attained as $\theta \to 0$); the uniformity of (c): $\sup\lvert v_j - v_j(0)\rvert/r_\Delta$ over 401 directions is $0.333, 0.125$ at $r_\Delta = 10^{-3}$, $0.346, 0.130$ at $0.1$, $0.5, 0.192$ at $1$ (defining series); 100 digits: (d) at 60 random points ($\sigma \in [-6, 6]$, $\theta \in [0.01, 6.3]$) $5.2\times10^{-99}$
relative, and $v_1^\circ$, $v_2^\circ$ against $\int te^{t\sigma}$, $\int\frac{t^2}2e^{t\sigma}$ (`mp.quad`) and $\frac12g_1 - g_2 + g_3$ from the defining series, $\le 1.2\times10^{-98}$. **Permanent:** planned, corpus `sim3_exp` (the reference is the definition); no coefficient id for $v_1$, $v_2$ is specified in `PHASE5.md` §3.

## 4. Rounding and the joint limit

**Proposition SM.7 (`f64`, measured; proposed strategy).** Sampled maxima, not proved bounds, with glibc's `exp`, `expm1`, `sin`, `cos` (preamble).

- (a) *A coefficient by the closed forms of SM.5(b).* $v_0$ errs by at most $1.9u$ (`expm1` is accurate; only $\sigma = 0$ is $0/0$). $v_2$ sums terms of size $1$ to a value $\Delta/6$, and $v_1$ terms of size $\lvert\sigma\rvert, z/2$ to a value $\Delta/2$:
relative errors $\le\hat C_2u/\Delta$ and $\le\hat C_1u(\lvert\sigma\rvert + z)/\Delta$ with $\hat C_2 \le 26$, $\hat C_1 \le 8.7$ for $r_\Delta \le 1$ (56 and 16 on $[1, 3]$; $8$ and $3$ on the $\theta$ axis, $\sigma = 0$: $b$ and $a$), on every ray tested ($\kappa \in \{0, 0.01, 0.1, 1, 10, \infty\}$). They keep no digit of $v_2$ below $r_\Delta \approx 10^{-8}$. At $\Delta = 0$
(or once $\Delta$ underflows, $\lvert\sigma\rvert, \theta \lesssim 10^{-162}$) they are $0/0$, and at the safe argument $\Delta \to 1$ they return $0$ instead of $\frac16$: the limit values $(\frac12, \frac16)$ have to be selected there.
- (b) *A value.* $\lVert\hat{\mathsf V} - \mathsf V\rVert_F/\lVert\mathsf V\rVert_F \le 2.3u$ for $\Delta < 1$ and $\le 4.4u$ for $\Delta \ge 1$ from the closed forms **at every $(\sigma, \theta)$**, $r_\Delta$ down to $10^{-12}$: $W$ multiplies $v_1$ and $W^2$ multiplies $v_2$, so $\theta\lvert\delta v_1\rvert$ and $z\lvert\delta v_2\rvert$ are $O(u)$ for $r_\Delta \le 1$ (SE.15(a) again). The same holds for $\mathrm{Log}$'s $\mathsf V^{-1}t$ (SM.8: $\le 11.4u$, relative to $\lVert\rho\rVert$). For $\Delta \ge 1$ use form b of $v_1$ (a branch between the two forms; $\Delta = 1$ is an illustration, the switch is generated, 0004): form a reaches $660u$ for $\lvert\sigma\rvert$ up to $700$. Both forms overflow before $\ln\mathrm{MAX}$ (g).
- (c) *A Jacobian block.* $\mathsf Q$ contains $-v_1\Xi$, which nothing suppresses, and $\Gamma_2(M)\rho$ contains $w_1\mu$, suppressed by one power of $\theta$ only (SM.10–SM.11): with the closed forms alone both blocks err by $\lesssim4.5u/r_\Delta$ for $r_\Delta < 0.5$ ($\lVert\rho\rVert$ of order $1$): $2\times10^4u$ at $r_\Delta = 10^{-4}$, $2\times10^{12}u$ at $10^{-12}$ (measured), $\le28u$ on $[0.1, 1)$. $\mathrm{SE}(3)$'s $Q$ loses one power of $\theta$ where its $e$ loses four (SE.15); here it is one power of $r_\Delta$.
- (d) *What removes it.* Two arms, both measured on the same $160$ points ($\sigma \in \{0, \pm10^{-12}, \pm10^{-8}, \pm10^{-4}, \pm10^{-2}, \pm0.3, \pm1, \pm3, 10\}$ × $\theta \in \{0, 10^{-12}, 10^{-8}, 10^{-4}, 10^{-2}, 0.3, 1, 2, 3, \pi - 10^{-6}\}$, the origin included, three random $(\varphi, \rho)$ each, $\rho \sim N(0, I_3)$) and on $600$ random rays with $r_\Delta$ log-uniform in $[0.02, 2]$ ($711$ of the $1080$ samples have $\Delta < 1$):
  - (i) *Split*: $v_1, v_2$ by SM.6(d), $w_0$ by its own series, $w_1, w_2, \partial_zv_j$ by the closed forms of SM.11 from them: blocks $\le6.1u$ for $\Delta < 1$, including $\sigma = \theta = 10^{-12}$; at the origin the limits are selected and the blocks are exact. It needs only functions of one variable, $g_2$, $v_1^\circ$, $v_2^\circ$ (exact arms lose $\le2.0u/\lvert\sigma\rvert$, $5.0u/\lvert\sigma\rvert$, $18.5u/\sigma^2$: a series in $\sigma$, radius $\infty$, constants $1/(i+\ell)!$), and the catalogue's $a, b, d$.
  - (ii) *The series of SM.6(b)* for all eight scalars from one recurrence, $K$ terms, at a depth that is **per sequence**: $\xi_k$ is homogeneous of weighted degree $k - 1$ and $\zeta_k$ of $k - 2$ (SM.3), so each $\partial_z$ drops two,
    and a uniform total degree $< K$ — which is what makes SM.6(b)'s tail a function of $r_\Delta$ alone, hence one mask on $\Delta$ — needs $\xi_k$ up to $k = K$ ($v_1$, $w_1$), $\zeta_k$ up to $K + 1$ ($v_2$, $w_2$), and
    $\xi'_k$, $\zeta'_k$ up to $K + 2$ and $K + 3$ (SM.11(c)'s recurrence). Taking every sequence but $\xi_k$ to $K + 1$ instead leaves $\partial_zv_1$ at degree $K - 2$ and $\partial_zv_2$ at $K - 3$, outside that bound: at
    $K = 20$, $r_\Delta = 1.6$ the value scalars err $1.4u$ and $1.8u$ while $\partial_zv_1$ errs $28u$ and $\partial_zv_2$ $401u$, and the two extra steps — $O(1)$ work — bring them to $2.1u$ and $2.5u$. With $K = 20$ blocks
    $\le5.5u$ for $\Delta < 1$ (and $\le7.3u$ up to $r_\Delta = 1.6$), exact at the origin; the table below is at the depth as first stated, so what it measures for the two derivative scalars is the deficit, not the
    rule above. Its truncation error meets the closed form's $4.5u/r_\Delta$ where $(K+1)r_\Delta^K/(K+2)! \approx 4.5u/r_\Delta$: $r_\times \approx 0.08, 0.38, 0.91, 1.6$ for $K = 8, 12, 16, 20$ (error at the
    crossing $\approx54u, 12u, 4.9u, 2.8u$). Measured maxima of the blocks, in $u$:

| $r_\Delta$ | samples | closed | split | $K = 8$ | $12$ | $16$ | $20$ | $24$ |
|---|---|---|---|---|---|---|---|---|
| $[0, 0.1)$ | 354 | $\lesssim4.5/r_\Delta$ | 5.1 | 81 | 5.5 | 5.5 | 5.5 | 5.5 |
| $[0.1, 0.35)$ | 227 | 27 | 6.1 | $5\times10^6$ | 5.0 | 4.6 | 4.6 | 4.6 |
| $[0.35, 0.5)$ | 42 | 11 | 4.7 | $1\times10^8$ | 397 | 4.5 | 4.5 | 4.5 |
| $[0.5, 0.8)$ | 56 | 8.2 | 5.4 | $3\times10^9$ | $3\times10^4$ | 3.3 | 3.3 | 3.3 |
| $[0.8, 1)$ | 32 | 8.1 | 4.3 | $1\times10^{10}$ | $4\times10^5$ | 8.5 | 4.2 | 4.2 |
| $[1, 1.2)$ | 90 | 7.1 | 7.3 | $7\times10^{10}$ | $6\times10^6$ | 143 | 7.3 | 7.3 |
| $[1.2, 1.6)$ | 49 | 7.2 | 5.3 | $5\times10^{11}$ | $7\times10^7$ | $1\times10^4$ | 5.8 | 5.8 |
| $[1.6, 2)$ | 23 | 8.2 | 6.8 | $5\times10^{12}$ | $4\times10^9$ | $1\times10^6$ | 97 | 5.4 |

- (e) Above $r_\Delta \approx 1$ the closed forms give $\le8.2u$ for $r_\Delta < 3$ (table above) and $\le54u$ for $\sigma \le 10$ (the $168$ samples with $r_\Delta \ge 3$), growing with $\sigma > 0$ and not for $\sigma < 0$: the cancellation between $[y]_\times J_l(\varphi)$ and $\partial_\varphi y$ in SM.10; beyond $\sigma = 10$, see (g).
- (f) *What is not measured:* `f32` (except (g)'s $\mathrm{Log}$ at $\theta \ge 1$), and the derivative through `Dual`: through the split, or through any of these arms, the weights' derivatives $\pm\sigma^2/\Delta^2$ multiply a difference of two equal values of size $\frac16$, so a `Dual` derivative of $v_j$ loses one more power of $\Delta$ (CO.14): `jacobians_match_dual_*` needs a per-stratum tolerance.

- (g) *Large $\sigma$: overflow and the scaled arrangement.*
  - (i) *Where the closed forms fail*, `f64` (`f32` in brackets; $\theta \ge 1$), long before $e^\sigma$ overflows: $v_0\Lambda$, the denominator of $i_2$ (SM.8), overflows at $\sigma \approx 242.1$ ($33.1$) and $i_2$ becomes $0$ **silently**, so $\mathrm{Log}$'s translation is off by $10^{13}u$ to $10^{16}u$ ($10^7u_{32}$); $\Lambda$ overflows at $360.8$ ($48.3$); $\sigma e^\sigma\,\mathrm{sinc}\,\theta$ and $\sigma e^\sigma a$, hence $v_1$ and $v_2$ in both forms, overflow at $703.25$ ($84.3$), although $e^\sigma$ and $v_j \lesssim e^\sigma/\sigma$ are finite up to $\ln\mathrm{MAX}$.
  - (ii) *The scaled arrangement*, for $\sigma > 0$: factor $e^\sigma$ out, i.e. in SM.5(b) and SM.11(b)–(c) put $e^\sigma \to 1$, $\mathrm{expm1}\,\sigma \to -\mathrm{expm1}(-\sigma)$ and $w_0 \to e^{-\sigma}g_2(\sigma)$; every scalar is then $\hat c = e^{-\sigma}c$, and $c = e^\sigma\hat c$ is formed last (finite whenever $e^\sigma$ and the result are). For $\mathrm{Log}$, $\mathsf V^{-1} = e^{-\sigma}\hat{\mathsf V}^{-1}$ ($i_0, i_1, i_2$ of SM.8 on the $\hat v_j$, $\hat\Lambda = e^{-2\sigma}\Lambda$), applied to $t$ as $e^{-\sigma/2}(e^{-\sigma/2}t)$: $\hat{\mathsf V}^{-1}t$ alone overflows for $t \approx e^\sigma$, and $e^{-\sigma}$ is subnormal above $708.4$. For $\sigma < 0$ nothing overflows ($e^\sigma \le 1$).
  - (iii) *Measured*, one sample per point: $\mathsf V$ within $4.4u$ for $\sigma \in [300, 709.7]$ (unscaled: $\infty$ from $704$); $\mathrm{Log}$'s $\rho$ within $10.8u$ for $\sigma \in [10, 709.7]$ and $2.8u$ for $\sigma \in [-700, -1]$ ($\theta \in \{0.1, 1, 3, \pi - 10^{-6}\}$; `f32`, $\theta \ge 1$: $12u_{32}$ for $\sigma \in [1, 88.5]$). Blocks, $\lvert\sigma\rvert \ge 10$: $\mathsf V \le 7.2u$; for $\sigma > 0$ $\Gamma_2(M)\rho \lesssim 1.5\sigma u$ and $\mathsf Q \lesssim 0.6\sigma^2u$ ($56u$ at $10$, $4.2\times10^3u$ at $100$, $3.9\times10^4u$ at $300$, $1.1\times10^5u$ at $709.5$); $\le0.5u$ for $\sigma < 0$ down to $-709.5$. The assembly of SM.10(b) alone (exact coefficients rounded to `f64`) gives $\mathsf Q$ within $28u$, $454u$, $2.7\times10^3u$ at $\sigma = 10, 100, 700$, the cancellation of (e); the closed-form coefficients ($w_1 = v_2 - \sigma w_2$ among them) cost the rest. $i_2$ is noise below $\sigma \approx -45$, without consequence (SM.8). $e^{\pm\sigma}$, $\mathrm{expm1}(\pm\sigma)$ and $e^{-\sigma/2}$ are the two `Real` methods of SM.4.
- (h) *The evaluation in one table* (`f64`; every switch below is generated, 0004, and the values are illustrations):

| Region | Scalars: arm | Masks | Blocks, measured |
|---|---|---|---|
| $\Delta = 0$ | all eight: the limits $1, \frac12, \frac16, \frac12, \frac16, \frac1{24}, -\frac1{24}, -\frac1{120}$ ($v_0, v_1, v_2, w_0, w_1, w_2, \partial_zv_1, \partial_zv_2$) | select on $\Delta = 0$; safe argument $\Delta \to 1$ in the closed arm | exact |
| $r_\Delta < r_\times$ | all eight: **series** (d)(ii), one recurrence, $K = 16$ below $r_\times \approx 0.8$ (or $12$ below $0.35$; $20$ below $1.6$, $\le7.3u$) | $\Delta < r_\times^2$ only | $\le5.5u$ |
| same, alternative | **split** (d)(i): $v_1, v_2$ by SM.6(d); $g_\ell$, $w_0$ exact for $\lvert\sigma\rvert \ge 0.5$, series in $\sigma$ below; $a, b, d$ of §4; $w_1, w_2, \partial_zv_j$ by SM.11 | $\Delta < r_\times^2$, $\Delta = 0$ (the weights), $\lvert\sigma\rvert$, $\theta^2$ | $\le6.1u$ |
| $r_\Delta \ge r_\times$ | all eight: **closed forms** SM.5(b), SM.11(b)–(c), $v_1$ by form b for $\Delta \ge 1$ (form a below), $w_0 = g_2(\sigma)$, $v_0 = 1$ at $\sigma = 0$; $\sigma > 0$: scaled (g) | $\sigma > 0$ (arrangement), $\theta^2$ (for $a, b, d$), $\lvert\sigma\rvert$ (for $w_0$), $\sigma = 0$ | $\le12u$ for $0.35 \le r_\Delta < 3$; $\le54u$ for $\sigma \le 10$; $\mathsf Q \lesssim0.6\sigma^2u$ beyond (g) |

Both arms sit inside the closed arm's masks. Beyond them the series adds one branch variable ($\Delta$, which also covers $\Delta = 0$: no division) and one recurrence; the split adds two masks ($\Delta < r_\times^2$ and the select at $\Delta = 0$) and three one-variable rows ($g_\ell$, $v_1^\circ$, $v_2^\circ$, each with its own $\lvert\sigma\rvert$ branch) beside §4's. Which one `xtask thresholds` sweeps is `0004`'s question (index, open items).

**Checked:** `f64` (see the preamble); references from the defining series (60 digits plus $0.9\lvert\sigma\rvert$, since the series cancels by $e^{2\lvert\sigma\rvert}$ for $\sigma < 0$) and, for the blocks of (c)–(e), (g), $J_l = g_1(\mathrm{ad}_\tau)$ as the top-right block of `mp.expm` of $\begin{bmatrix}\mathrm{ad}_\tau & I\\ 0 & 0\end{bmatrix}$ ($90 + \lvert\sigma\rvert$ digits; it agrees with the dense series of $\mathrm{ad}$ to $3\times10^{-60}$ at one point), the error of each of $\mathsf V$, $\mathsf Q$, $\Gamma_2(M)\rho$ being the Frobenius error over $\max(1, \lVert\cdot\rVert_F)$, the largest of the three reported.
(a) $v_1$ by form a and unscaled; 400 samples per (ray, $r_\Delta$ range), $[10^{-6}, 10^{-3}]$, $[10^{-3}, 0.1]$, $[0.1, 1]$, $[1, 3]$; the constants are the largest of $\text{error}\cdot\Delta/u$ resp. $\text{error}\cdot\Delta/(u(\lvert\sigma\rvert + z))$.
(b) $41$ values of $\sigma$ ($0$, $\pm10^{-e}$ for $e = 1..12$, $\pm0.5, \pm1, \pm3, \pm10, \pm30, \pm100, \pm300, \pm700$) × $19$ of $\theta$ ($0$, $10^{-e}$ for $e = 1..12$, $0.3, 1, 2, 3, \pi - 10^{-9}, \pi$) without $(0, 0)$ (the select gives $\mathsf V = I$ there), and 300 random ($\sigma \in [-30, 30]$, $\theta \in [0, \pi]$); the split: $\le1.6u$ for $\Delta < 1$, $11u$ for $\Delta \ge 1$.
(c), (d) as stated (the table: `f64` runs of all arms on the same samples). (d)(ii)'s depth: mpmath 1.4.1, 120 digits, the recurrences of SM.3 and SM.11(c) against a $K = 110$ reference, relative error in $u$ over
$r_\Delta \in \{0.9, 1, 1.2, 1.6\}$ × four ray angles ($1^\circ, 30^\circ, 60^\circ, 89^\circ$); the quoted row is $r_\Delta = 1.6$, $1^\circ$ (the largest), where $v_1, v_2$ give $1.38u, 1.79u$, $\partial_zv_1, \partial_zv_2$
at $K + 1$ give $28.0u, 401u$, and at $K + 2$, $K + 3$ give $2.14u, 2.46u$; the recurrences were first checked to reproduce SM.6(b)'s and SM.11(d)'s printed numerators exactly (`xi_3`, `xi_5`, `zeta_6` and the five
leading-term rows, exact rationals). The crossing radii are the estimate $4.5u/r_\Delta = (K+1)r_\Delta^K/(K+2)!$, compared with the table. The one-variable functions: 300 samples per decade of $\lvert\sigma\rvert \in [10^{-6}, 1]$, both signs; the split at $\sigma = 100$, $\theta = \pi$ loses $145u$ in $v_1$ (cancellation $a - \sigma b$), which is why $\Delta \ge 1$ uses the direct forms. (g): $\sigma \in \{\pm10, \pm20, \pm30, \pm100, \pm300, \pm700, \pm709.5\}$ × $\theta \in \{0, 10^{-8}, 10^{-4}, 0.3, 1, 3, \pi - 10^{-6}\}$, two samples each, for the blocks (dps $400 + 2\lvert\sigma\rvert$ for the exact coefficients); $\mathsf V$: $\sigma \in \{300, 355, 400, 600, 700, 703, \dots, 709.7\}$ and $-100, -700, -709, -720$ against `mp.expm` of the $6\times6$ matrix $\begin{bmatrix}M & I\\ 0 & 0\end{bmatrix}$ ($300 + \lvert\sigma\rvert$ digits); $\mathrm{Log}$ as SM.8, $t = \mathrm{fl}(\mathsf V\rho_0)$; the overflow onsets by a scan in steps of $0.05$; `f32` with numpy `float32`. Scripts not committed; the `sim3_*` strata will replace them, through `libm`.
**Permanent:** none specified for the strategy; the strata of `sim3_exp`, `sim3_jr` (`PHASE5.md` §3: $\sigma \in \{0, \pm10^{-12}, \dots, \pm3\}$ × the `theta:*` decades, several rays) measure the quantities themselves, and (c) says they cannot pass with the closed forms alone.

## 5. Log

**Proposition SM.8 (Log; proposed).** Let $X = (q, t, \sigma)$, $R = R(q)$, $\varphi = \mathrm{Log}\,q$ (`NUMERICS.md` §3.2, $\theta \le \pi$). Then
$$
\mathrm{Log}\,X = (\varphi,\ \mathsf V(\varphi,\sigma)^{-1}t,\ \sigma),\qquad
\mathsf V^{-1} = i_0I + i_1W + i_2W^2,\quad i_0 = \frac1{v_0},\ \ i_1 = -\frac{v_1}\Lambda,\ \ i_2 = \frac{v_1^2 - v_2\chi}{v_0\Lambda},
$$
$\chi = v_0 - zv_2$, $\Lambda = \chi^2 + zv_1^2$ (SM.5(c)): no second closed form. It inverts $\mathrm{Exp}$ on $\{\theta < \pi\}$ for every $\sigma$ (the arithmetic is another matter: SM.7(g), SM.13); $\sigma$ is read off the storage. At $\sigma = 0$ it is $J_l^{-1}(\varphi) = I - \frac12W + cW^2$ (SO.8). $\mathsf V$ is normal, so $\lVert\mathsf V^{-1}\rVert_2 = 1/\min(g_1(\sigma), \sqrt\Lambda)$: about $\lvert\sigma\rvert$ as $\sigma \to -\infty$. At $\theta(X) = \pi$ there are two logarithms, $(\varphi, \rho)$ and $(-\varphi, \rho')$, $\rho' = \mathsf V(-\varphi, \sigma)^{-1}t$; for $\sigma = 0$, $\rho' - \rho = \varphi\times t$ (LG.2(d)); which is returned is the quaternion sign (§3.2).

*Proof.* Exp's rotation and scale blocks give $\varphi$ (SO.5) and $\sigma$; $\mathsf V\rho = t$ has the unique solution $\mathsf V^{-1}t$ (SM.5(c), $\theta \le \pi$). The product in $\mathrm{span}\{I, W, W^2\}$ is
$(x_0, x_1, x_2)(y_0, y_1, y_2) = (x_0y_0,\ x_0y_1 + x_1y_0 - z(x_1y_2 + x_2y_1),\ x_0y_2 + x_1y_1 + x_2y_0 - zx_2y_2)$, so $\mathsf V\mathsf V^{-1} = I$ is $v_0i_0 = 1$ and
$\begin{bmatrix}\chi & -zv_1\\ v_1 & \chi\end{bmatrix}\binom{i_1}{i_2} = -i_0\binom{v_1}{v_2}$, determinant $\Lambda$; Cramer with $\chi + zv_2 = v_0$. At $\sigma = 0$: $v = (1, a, b)$, $\Lambda = \mathrm{sinc}^2\theta + za^2 = 2a$ and $i_2 = (a^2 - b\,\mathrm{sinc}\,\theta)/2a = c$ (SO.8). $\square$

*Consequence.* Near the joint limit no cancellation beyond that of $v_0, v_1, v_2$: $\Lambda$ is a sum of squares and $v_1^2 - v_2\chi \to \frac1{12}$. Far from it there is one: for $\sigma \ll 0$ the leading terms of $v_1^2$ and $v_2\chi$ coincide, the numerator of $i_2$ cancels totally and $i_2$ is noise, of absolute error $\le2.3u\lvert\sigma\rvert/\Delta$ (measured, $\sigma \in [-300, -3]$, $\theta \in \{0.1, 1, 3\}$) against a true value $O(e^\sigma)$, which it exceeds below $\sigma \approx -45$ (`f64`). It multiplies $W^2$, of norm $z$, and $\lVert\mathsf V\rVert_2 = g_1(\sigma) \approx 1/\lvert\sigma\rvert$ there, so its contribution to $\rho$ is $\le2.3u\,(z/\Delta)\lVert\rho\rVert$: the translation is unaffected. $\mathrm{Log}$'s translation part is $\le11.4u$ (relative to $\lVert\rho\rVert$) from the closed forms of SM.5, and $\le2.9u$ for $\Delta < 1$ ($\le2.3u$ with the $K = 16$ series below $\Delta = 0.64$), given $\varphi$ (`f64`, $510$ cases: $17$ values of $\sigma$ in $\{0, \pm10^{-12}, 10^{-8}, \pm10^{-4}, 10^{-2}, \pm0.3, \pm1, \pm3, \pm10, \pm30\}$ × $\theta \in \{0, 10^{-12}, 10^{-8}, 10^{-4}, 10^{-2}, 0.3, 1, 2, 3, \pi - 10^{-6}\}$, three random $(\varphi, t)$ each, the origin included: exact, the limits are selected). For $\sigma > 0$ large the recipe needs the scaled arrangement (SM.7(g)).

**Checked:** mpmath 1.3.0, 100 digits, 55 $(\theta, \sigma)$ ($\theta \in \{10^{-3}, 0.5, 1.5, 3, 3.1\}$ × $\sigma \in \{0, 10^{-7}, -0.5, 2, -6, 25, -25\}$ and 20 random): $\mathrm{Log}(\mathrm{Exp}\,\tau)$ by the geometric $\mathrm{SO}(3)$ logarithm ($\mathrm{atan2}$ of the antisymmetric part and the trace) and `mp.lu_solve` against the series $\Gamma_1$: $\varphi$ $1.3\times10^{-100}$, $\rho$ $7.4\times10^{-92}$ (relative to $\max(1, \lVert\rho\rVert)$); the closed $(i_0, i_1, i_2)$ against `mp.inverse` $2.2\times10^{-101}$; 16 tangents at $\theta = \pi$
($\sigma \in \{0, 0.7, -2.5, 9\}$): $\mathrm{Exp}(-\varphi, \rho') = X$ and the $\sigma = 0$ relation $1.7\times10^{-100}$.
**`mp.logm` of the $4\times4$ matrix is not a reference** ($\mathrm{Exp}\,\tau$ built from $\theta \in \{0.1, \dots, 3.05\}$, $\sigma \in \{0, 0.1, 0.5, 1, 2, 3, -1, -3, 6\}$, 5 axes, 40 and 80 digits identically): it returns a complex non-principal result from an onset that falls with $\lvert\sigma\rvert$ (8 axes, step $0.01$, 40 digits):
$(3.02, 3.03]$ at $\sigma = 0$, $(3.00, 3.01]$ at $1$, $(2.96, 2.97]$ at $2$, $(2.85, 2.86]$ at $3$, $(2.88, 2.89]$ at $-3$, $(2.19, 2.20]$ at $6$ (as for $\mathrm{SO}(3)$, SO.5, and $\mathrm{SE}_N(3)$, SE.3). `f64` (the 510 cases above, and SM.7(g)): $\hat{\mathsf V}$ from the closed forms or the $K = 16$ series below $r_\Delta = 0.8$, error relative to $\lVert\rho\rVert$, reference `mp.lu_solve` of $\mathsf V$ from `mp.expm` of $\begin{bmatrix}M & I\\ 0 & 0\end{bmatrix}$; $i_2$ against the exact $i_2$ (200 + $\lvert\sigma\rvert$ digits) at 21 points, $\sigma \in \{-3, \dots, -300\}$.
**Permanent:** planned, corpus `sim3_log`, proptest `exp_log_roundtrip_*`; its reference must not be `mp.logm`.

## 6. Ad and ad

**Proposition SM.9 ($\mathrm{Ad}$ stated in §9; the rest proposed).** For $X = (R, t, \sigma)$, $\tau = [\varphi;\rho;\sigma]$, $\psi = [\alpha;\eta;\varsigma]$, $\Omega = \alpha^\wedge$ (no domain):
$$
\mathrm{Ad}_X\psi = \big[R\alpha;\ sR\eta + [t]_\times R\alpha - \varsigma t;\ \varsigma\big],\qquad
\mathrm{ad}_\tau\psi = \big[\varphi\times\alpha;\ \varphi\times\eta + \rho\times\alpha + \sigma\eta - \varsigma\rho;\ 0\big],
$$
$$
\mathrm{Ad}_X = \begin{bmatrix}R & 0 & 0\\ [t]_\times R & sR & -t\\ 0 & 0 & 1\end{bmatrix},\quad
\mathrm{ad}_\tau = \begin{bmatrix}W & 0 & 0\\ \Xi & W + \sigma I & -\rho\\ 0 & 0 & 0\end{bmatrix},\quad
\mathrm{Ad}_X^{-1} = \begin{bmatrix}R^\top & 0 & 0\\ -s^{-1}R^\top[t]_\times & s^{-1}R^\top & s^{-1}R^\top t\\ 0 & 0 & 1\end{bmatrix},
$$
$\det\mathrm{Ad}_X = e^{3\sigma}$, $\mathrm{Ad}_X^{-1} = \mathrm{Ad}_{X^{-1}}$, and $\kappa_2(\mathrm{Ad}_X) \ge e^{\lvert\sigma\rvert}$. **There is no dual structure** (SE.10): the diagonal blocks are $R, sR, 1$ for $\mathrm{Ad}$, $W, W + \sigma I, 0$ for $\mathrm{ad}$ and $J_l(\varphi), \mathsf V, 1$ for $J$ (SM.10),
unequal, so the shape of SE.10 is lost. The pattern $\begin{bmatrix}A & 0 & 0\\ B & C & \mathbf f\\ 0 & 0 & 1\end{bmatrix}$ ($\mathbf f$ a column) is closed under product, $(AA',\ BA' + CB',\ CC',\ C\mathbf f' + \mathbf f)$ ($4\cdot27 + 9 = 117$ multiplications against $343$ dense), and under inverse: an optimization, not a formula; `NUMERICS.md` §9, 0005 and `PHASE5.md` §3 say dense $7\times7$.
For a translation-first tangent $[\rho;\varphi;\sigma]$ (the scale kept last: an assumption to check per source), $\Pi = \mathrm{diag}(\Pi_1, 1)$ converts and $N \mapsto \Pi N\Pi^\top$ (SE.14(b)), e.g. $\mathrm{Ad}_{\mathrm{tf}} = \begin{bmatrix}sR & [t]_\times R & -t\\ 0 & R & 0\\ 0 & 0 & 1\end{bmatrix}$.

*Proof.* $X\psi^\wedge X^{-1} = \begin{bmatrix}R\Omega R^\top + \varsigma I & sR\eta - (R\alpha)^\wedge t - \varsigma t\\ 0 & 0\end{bmatrix}$ with $R\Omega R^\top = (R\alpha)^\wedge$ and $-(R\alpha)\times t = [t]_\times R\alpha$.
The commutator $[\tau^\wedge, \psi^\wedge]$ has top-left $[W, \Omega] = (\varphi\times\alpha)^\wedge$ (scalars commute) and top-right $(W + \sigma)\eta - (\Omega + \varsigma)\rho$. $\mathrm{Ad}_X^{-1}$ is $\mathrm{Ad}$ at $X^{-1}$, using $[-s^{-1}R^\top t]_\times R^\top = -s^{-1}R^\top[t]_\times$.
$\lVert\mathrm{Ad}_X\rVert_2 \ge \max(1, s)$ and $\lVert\mathrm{Ad}_X^{-1}\rVert_2 \ge \max(1, s^{-1})$ (a submatrix), whose product is $\ge e^{\lvert\sigma\rvert}$. $\square$

**Checked:** mpmath 1.3.0, 60 digits, 24 $(\tau, X)$: $\mathrm{ad}$ from the commutators of the $4\times4$ hat matrices, exact; $\mathrm{Ad}_X$ from $X\psi^\wedge X^{-1}$ $1.5\times10^{-61}$; $\mathrm{Ad}_{\mathrm{Exp}\,\tau} = \exp(\mathrm{ad}_\tau)$ (`mp.expm`) $1.1\times10^{-61}$. 110 digits, 30 random $(X, Y)$: $\mathrm{Ad}_X\mathrm{Ad}_Y = \mathrm{Ad}_{XY}$ $2.6\times10^{-111}$;
the block $\mathrm{Ad}^{-1}$, `mp.inverse` and $\mathrm{Ad}_{X^{-1}}$ agree to $1.7\times10^{-111}$. $\kappa_2$ by `mp.svd_r` (40 digits) at $\sigma \in \{0, 1, -1, 3, -3, 6\}$, one $X$ each: $7.1, 4.1, 14.8, 20.4, 117.8, 403.4$ against $e^{\lvert\sigma\rvert} = 1, 2.7, 2.7, 20.1, 20.1, 403.4$.
**Permanent:** planned, corpus `sim3_ad`, proptest `adjoint_identity_*`.

## 7. The Jacobians

**Proposition SM.10 (proposed).** For every $\tau$, $J_l = g_1(\mathrm{ad}_\tau) = \int_0^1\mathrm{Ad}_{\mathrm{Exp}(t\tau)}\mathrm dt$ and $J_r(\tau) = J_l(-\tau)$ (LG.7–LG.9, SM.2(c)), and
$$
J_l(\tau) = \begin{bmatrix}J_l(\varphi) & 0 & 0\\ \mathsf Q & \mathsf V & -\Gamma_2(M)\rho\\ 0 & 0 & 1\end{bmatrix},\qquad
J_r(\tau) = J_l(-\tau) = \begin{bmatrix}J_r(\varphi) & 0 & 0\\ \mathsf Q(-\rho, -\varphi, -\sigma) & \mathsf V(-\varphi, -\sigma) & \Gamma_2(-M)\rho\\ 0 & 0 & 1\end{bmatrix}.
$$
- (a) $\Gamma_2(M) = g_2(M) = \sum M^k/(k+2)! = w_0I + w_1W + w_2W^2$ (SM.11); $\mathsf Q = \mathsf Q(\rho, \varphi, \sigma)$ below.
- (b) With $y = \mathsf V\rho$, the translation part of $\mathrm{Exp}\,\tau$: $\partial_\sigma y = y - \Gamma_2(M)\rho$ and
$$
\mathsf Q = [y]_\times J_l(\varphi) + \partial_\varphi y,\qquad
\partial_\varphi y = 2\big(\partial_zv_1\,\mu + \partial_zv_2\,\nu\big)\varphi^\top - v_1\Xi + v_2(\Xi W - 2W\Xi),
$$
$\partial_z$ at fixed $\sigma$ ($v_0$ has no $z$ dependence). Equivalently $\mathsf Q = \sum_{i,j\ge0}M^i\,\Xi\,W^j/(i+j+2)! = \int_0^1\!\!\int_0^te^{u\sigma}e^{uW}\,\Xi\,e^{(t-u)W}\mathrm du\,\mathrm dt$; at $\sigma = 0$ it is $Q(\rho, \varphi)$ of §5.3 and the upper-left $6\times6$ block is $\mathrm{SE}(3)$'s $J_l$; at $\varphi = 0$ everything is a function of $\sigma$: $\mathsf Q = g_2(\sigma)\Xi$, $\mathsf V = g_1(\sigma)I$, $\Gamma_2(M) = g_2(\sigma)I$.
- (c) $\det J = \dfrac{2(1 - \cos\theta)}{\theta^2}\,g_1(\sigma)\,\Lambda$: invertible iff $\theta \notin 2\pi\mathbb Z_{>0}$ (the rotation block does not see $\sigma$), and
$$
J_l^{-1} = \begin{bmatrix}J_l^{-1}(\varphi) & 0 & 0\\ -\mathsf V^{-1}\mathsf QJ_l^{-1}(\varphi) & \mathsf V^{-1} & \mathsf V^{-1}\Gamma_2(M)\rho\\ 0 & 0 & 1\end{bmatrix},\qquad J_r^{-1} - J_l^{-1} = \mathrm{ad}_\tau,\quad J_l = \mathrm{Ad}_{\mathrm{Exp}\,\tau}J_r,
$$
with $J_l^{-1}(\varphi)$ from SO.8 and $\mathsf V^{-1}$ from SM.8: the inverse is block algebra, not a second closed form.

*Proof.* $\mathrm{ad}^{k+1} = \mathrm{ad}^k\mathrm{ad}$ with SM.9: the $\varphi$ block of $\mathrm{ad}^k$ is $W^k$ and the $\rho$ block $M^k$; the last column of $\mathrm{ad}^{k+1}$ is $[0; -M^k\rho; 0]$; block $(\rho, \varphi)$ obeys $B_{k+1} = B_kW + M^k\Xi$, $B_1 = \Xi$, so $B_k = \sum_{i+j=k-1}M^i\Xi W^j$. Summing with $1/(k+1)!$ gives the block form, (a) and the double sum; the double integral is block $(\rho,\varphi)$ of
$\int_0^1\mathrm{Ad}_{\mathrm{Exp}(t\tau)}\mathrm dt$ (SM.9 with the translation $y_t = \int_0^te^{uM}\rho\,\mathrm du$ of $\mathrm{Exp}(t\tau)$ and $R_t = e^{tW}$: $[y_t]_\times R_t = \int_0^te^{u\sigma}R_u\Xi R_u^\top R_t\,\mathrm du$).
(b) LG.8: $\mathrm{Exp}(\tau + \delta) = \mathrm{Exp}(J_l\delta)\mathrm{Exp}(\tau) + O(\lVert\delta\rVert^2)$, and the translation part of $\mathrm{Exp}(\eta)E$ is $(1 + \eta_\sigma)y + \eta_\varphi\times y + \eta_\rho + O$. With $\delta = [\delta\varphi; 0; 0]$: $\partial_\varphi y = \mathsf Q - [y]_\times J_l(\varphi)$; with $\delta = [0; 0; \delta\sigma]$: $\partial_\sigma y = y - \Gamma_2(M)\rho$.
The explicit $\partial_\varphi y$ differentiates $y = v_0\rho + v_1\mu + v_2\nu$ with $\partial_\varphi z = 2\varphi^\top$, $\partial_\varphi\mu = -\Xi$, $\partial_\varphi\nu\,\delta = -[\mu]_\times\delta - W\Xi\delta$ and $[\mu]_\times = [W, \Xi]$. At $\sigma = 0$, $M = W$ and the double sum is SE.5's.
(c) The determinant of a block triangular matrix is $\det J_l(\varphi)\cdot\det\mathsf V$ (SO.7, SM.5(c)); the inverse by block substitution; $J_r^{-1} - J_l^{-1} = \mathrm{ad}_\tau$ because $1/g_1(-\zeta) - 1/g_1(\zeta) = \zeta$ (LG §6), and $J_l = \mathrm{Ad}J_r$ is LG.9(b). $\square$

**Checked:** mpmath 1.3.0. Joint-limit rays, 150 digits: $\sigma = \kappa\theta$, $\kappa \in \{0, 0.1, 1, 10, 10^3\}$ and $r_\Delta \in \{10^{-1}, 10^{-3}, 10^{-6}, 10^{-9}\}$, both signs of $\sigma$ (40 points): the blocks of `mp.expm` against SM.5 $3.6\times10^{-143}$, $J_l$ and $J_r$ from the blocks against the dense series $1.8\times10^{-130}$; the total-degree series of SM.6(b) ($40$ levels) against the closed forms $9.2\times10^{-90}$ (its truncation at $r_\Delta = 0.1$); $\partial_zv_2$ against `mp.diff` of that series $5.4\times10^{-59}$; $\varphi = 0$ (seven $\sigma$ from $-40$ to $5$ and $10^{-12}$, 90 digits): the dense series against the form above $3.3\times10^{-77}$.
84 $(\theta, \sigma)$ ($8\times8$ grid $\theta \in \{10^{-3}, \dots, 6\}$, $\sigma \in \{0, 10^{-4}, -0.3, 1, -2, 3, -30, 20\}$, and 20 random), 70 digits: $J_l$ from the blocks against the dense series $\sum\mathrm{ad}^k/(k+1)!$ of the $7\times7$ $\mathrm{ad}$ built from commutators $8.8\times10^{-61}$ (relative to $\max(1, \lVert J\rVert)$), $J_r$ $5.5\times10^{-65}$, $J_l = \mathrm{Ad}_{\mathrm{Exp}\,\tau}J_r$ $8.8\times10^{-61}$, $\det$ $2.8\times10^{-59}$.
30 $(\theta, \sigma)$, 90 digits: central differences ($h = 10^{-25}$) of the operation itself, $(X - I)^\vee$ for $X = \mathrm{Exp}(\tau + he_j)\mathrm{Exp}(\tau)^{-1}$ and $\mathrm{Exp}(\tau)^{-1}\mathrm{Exp}(\tau + he_j)$: $J_l$ $3.0\times10^{-51}$, $J_r$ $1.7\times10^{-51}$ (the differencing); the block inverse against `mp.inverse` $2.3\times10^{-91}$; $J_r^{-1} - J_l^{-1} - \mathrm{ad}_\tau$ $5.4\times10^{-89}$; $\mathsf Q$ against the double sum truncated at total degree $70$ $9.1\times10^{-56}$; at $\sigma = 0$ (10 tangents, $\theta \in [10^{-3}, 20]$, 90 digits), $\mathsf Q$ against §5.3 with $b, d, e$ from definitions $9.1\times10^{-86}$, and the upper-left $6\times6$ block against the dense $\mathrm{SE}(3)$ series $2.5\times10^{-85}$.
**Permanent:** planned, corpus `sim3_jr`, `sim3_jr_inv` (`PHASE5.md` §3: no `sim3_jl`, `sim3_jl_inv`), proptests `jacobians_match_dual_*`, `jl_is_ad_jr_*`.

## 8. The scalars of the blocks

**Proposition SM.11 (proposed).**

- (a) *The family.* $\Gamma_\ell(M) = g_\ell(M) = \gamma^\ell_0I + \gamma^\ell_1W + \gamma^\ell_2W^2$, $\gamma^\ell_0 = g_\ell(\sigma)$, $\gamma^\ell_1 = \sum_k\xi_k/(k+\ell)!$, $\gamma^\ell_2 = \sum_k\zeta_k/(k+\ell)!$, i.e.
$\gamma^\ell_1 = \sum_{m,i}\frac{(-z)^m\sigma^i}{(2m+1)!\,i!\,(2m+i+2)_\ell}$, $\gamma^\ell_2 = \sum_{m,i}\frac{(-z)^m\sigma^i}{(2m+2)!\,i!\,(2m+i+3)_\ell}$ ($(y)_\ell$ the rising factorial). $v_j = \gamma^1_j$, $w_j = \gamma^2_j$;
$\Gamma_\ell = \frac{I}{\ell!} + M\Gamma_{\ell+1}$ and $\partial_\sigma\gamma^\ell_j = \gamma^\ell_j - \ell\gamma^{\ell+1}_j$. At $\sigma = 0$: $(\gamma^\ell_1, \gamma^\ell_2) = (\sigma_{\ell+1}, \sigma_{\ell+2})$ of CO.1 (GG.2(a)).
- (b) *$\Gamma_2$.* $w_0 = g_2(\sigma) = (\mathrm{expm1}\,\sigma - \sigma)/\sigma^2$, $\ w_2 = \dfrac{w_0 + \sigma v_2 - v_1}\Delta$, $\ w_1 = v_2 - \sigma w_2$.
- (c) *Derivatives.* With $a' = \mathrm da/\mathrm dz = -\frac12(b - 2d)$ and $\mathrm{sinc}' = -\frac12(a - b)$ (CO.4(b)):
$$
\partial_zv_2 = \frac{\sigma e^\sigma a' - e^\sigma\,\mathrm{sinc}' - v_2}\Delta,\qquad \partial_zv_1 = e^\sigma a' - \sigma\,\partial_zv_2,
$$
$\partial_zv_1 = -\sum\frac{(-1)^m(m+1)z^m\sigma^i}{(2m+3)!\,i!\,(2m+i+4)}$, $\partial_zv_2 = -\sum\frac{(-1)^m(m+1)z^m\sigma^i}{(2m+4)!\,i!\,(2m+i+5)}$, from $\xi'_{k+1} = \sigma\xi'_k - \zeta_k - z\zeta'_k$, $\zeta'_{k+1} = \sigma\zeta'_k + \xi'_k$. At $\sigma = 0$: $\partial_zv_1 = a'$ and $\partial_zv_2 = -e$: the joint-limit form of §4's $e$.
- (d) *Leading terms.*
$w_0 = \frac12 + \frac\sigma6 + \frac{\sigma^2}{24}$, $\ w_1 = \frac16 + \frac\sigma{12} + \frac{3\sigma^2 - z}{120}$, $\ w_2 = \frac1{24} + \frac\sigma{40} + \frac{6\sigma^2 - z}{720}$, $\ \partial_zv_1 = -\frac1{24} - \frac\sigma{30} + \frac{z - 5\sigma^2}{360}$, $\ \partial_zv_2 = -\frac1{120} - \frac\sigma{144} + \frac{2z - 15\sigma^2}{5040}$.
- (e) *Which cancels.* $w_0$: $2.0u/\lvert\sigma\rvert$ (SM.7(d)). $w_2$ and $\partial_zv_2$ divide by $\Delta$ a difference of $O(1)$ terms that is $O(\Delta)$: absolute error $\approx u/\Delta$ ($\lvert\sigma\rvert u/\Delta$ for $\partial_zv_1$), suppressed in the blocks by $z$ ($w_2\nu$), $\theta^3$ ($\partial_zv_2\,\nu\varphi^\top$) and $\theta^2$ ($\partial_zv_1\,\mu\varphi^\top$). $w_1 = v_2 - \sigma w_2$ multiplies $\mu$ (one power of $\theta$) and needs $v_2$ accurate to $u/\theta$, which the closed form is not: the split or the series (SM.7(c)).

*Proof.* (a) Insert SM.3 into $\sum M^k/(k+\ell)!$; $\binom k{2m+1}/(k+\ell)! = 1/[(2m+1)!\,i!\,(k+1)_\ell]$ for $k = 2m+1+i$, likewise for $\zeta$; the recurrence and $\partial_\sigma$ from $(k+1)/(k+1+\ell)! = 1/(k+\ell)! - \ell/(k+\ell+1)!$ and $\partial_\sigma M^k = kM^{k-1}$. (b) $\Gamma_2M = \Gamma_1 - I$ in coordinates: $\sigma w_0 = v_0 - 1$, $w_0 + \sigma w_1 - zw_2 = v_1$,
$w_1 + \sigma w_2 = v_2$; eliminate $w_1$. (c) $\mathrm d\sigma_m/\mathrm dz = -\tau_{m+2}$ (CO.1), so $a' = -\tau_4$, $\mathrm{sinc}' = -\tau_3$, and CO.4(b) gives the stated forms. Differentiate $\Delta v_2 = v_0 + \sigma e^\sigma a - e^\sigma\mathrm{sinc}\,\theta$ (SM.5(b)) in $z$, $\partial_z\Delta = 1$, and $v_1 = e^\sigma a - \sigma v_2$; the series by $\partial_z$ of SM.6(a); at $\sigma = 0$ $v_2 = b$ and $b' = -e$ (CO.4(b)). $\square$

**Checked:** sympy 1.14.0: the general terms of (a) for $\ell = 1, 2, 3$, $k \le 9$; the differentiated recurrences equal $\partial_z\xi_k$, $\partial_z\zeta_k$ ($k \le 13$); the derivative series ($k \le 10$, $11$); the terms of (d). mpmath 1.3.0, 100 digits, 60 random $(\sigma, \theta)$ ($\sigma \in [-6, 6]$, $\theta \in [0.01, 6.3]$):
$w_0, w_1, w_2$ against the defining series $\sum(\sigma^k, \xi_k, \zeta_k)/(k+2)!$ $4.7\times10^{-99}$; $\partial_zv_1$, $\partial_zv_2$ against `mp.diff` in $z$ of the defining series $1.5\times10^{-99}$, $6.0\times10^{-98}$; $g_2, g_3$ and $\partial_\sigma g_2 = g_2 - 2g_3$ $3.3\times10^{-98}$. The closed forms enter SM.10's checks.
**Permanent:** none specified (the corpus `sim3_jr` measures the blocks); the sampled bounds $\lvert\xi'_k\rvert \le \binom k3r_\Delta^{k-3}$, $\lvert\zeta'_k\rvert \le \binom k4r_\Delta^{k-4}$ (ratios $\le 1$ on 3000 points) are not proved.

## 9. Action and domains

**Proposition SM.12 (action Jacobians, both sides; proposed rows of §2.4).** For $F(X, p) = sRp + t$, $p' = Xp$, in the sense of LG.11:
$$
\mathrm D^{\mathrm{right}}_XF = \big[-sR[p]_\times,\ sR,\ sRp\big],\qquad
\mathrm D^{\mathrm{left}}_XF = \big[-[p']_\times,\ I,\ p'\big],\qquad
\partial F/\partial p = sR,\qquad \mathrm D^{\mathrm{right}} = \mathrm D^{\mathrm{left}}\mathrm{Ad}_X .
$$

*Proof.* $\mathrm{Exp}(\delta)[p; 1] = [p + \delta_\varphi\times p + \delta_\rho + \delta_\sigma p;\ 1] + O(\lVert\delta\rVert^2)$. Right: $X\,\mathrm{Exp}(\delta)[p;1]$ has translation $sR(p + \delta_\varphi\times p + \delta_\rho + \delta_\sigma p) + t$; left: $p' + \delta_\varphi\times p' + \delta_\rho + \delta_\sigma p'$. The relation is LG.4(b); explicitly
$[-[p']_\times, I, p']\mathrm{Ad}_X$ has first block $-[p']_\times R + [t]_\times R = -[sRp]_\times R = -sR[p]_\times$ and last $-t + p' = sRp$. $\square$

**Checked:** mpmath 1.3.0, 110 digits, 30 random $(X, p)$: central differences ($h = 10^{-30}$, $4\times4$ `mp.expm`) of $X\,\mathrm{Exp}(\delta)[p;1]$ and $\mathrm{Exp}(\delta)X[p;1]$: right and left $1.7\times10^{-61}$ (the differencing); the relation $\mathrm D^{\mathrm{right}} = \mathrm D^{\mathrm{left}}\mathrm{Ad}_X$ on the predicted forms $3.7\times10^{-111}$; $\partial F/\partial p$ is linear.
**Permanent:** planned, proptest `jacobians_match_dual_*` (`act_jacobians`).

**Proposition SM.13 (domains; proposed row of §12).**

| Routine | Domain | Note |
|---|---|---|
| `exp`, `adjoint`, `ad`, `jr`, `jl`, `act`, `act_jacobians` | every $\tau$, $X$ with $e^{\pm\sigma}$ finite: $\lvert\sigma\rvert < \ln\mathrm{MAX}$, $709.78$ (`f64`), $88.72$ (`f32`); `exp`, `jr`, `jl` with the scaled arrangement of SM.7(g) (measured to $709.5$, `f64`) | "$\sigma$ finite" (§12) is necessary ($s^{\pm1}$, $X^{-1}$'s translation) and not sufficient: the direct closed forms overflow at $703.25$ (`f64`; $84.3$ `f32`), $\mathsf V\rho \approx e^\sigma\rho/\sigma$ only when that value does; $e^\sigma$ is subnormal below $-708.40$ and $0$ below $-745.13$ (`f64`) |
| `log` | every $X$ ($\theta(X) \le \pi$) with $\lvert\sigma\rvert < \ln\mathrm{MAX}$; a function of $q$'s sign at $\theta(X) = \pi$ (SM.8) | inverse of $\mathrm{Exp}$ on $\theta < \pi$ for any $\sigma$; the direct recipe is silently wrong from $242.1$ (`f64`; $33.1$ `f32`), the scaled one measured to $709.7$; beyond $\ln\mathrm{MAX}$ not measured |
| `jr_inv`, `jl_inv` | $\theta < 2\pi$ (SM.10(c)), $\lvert\sigma\rvert < \ln\mathrm{MAX}$ | as $\mathrm{SO}(3)$, $\mathrm{SE}_N(3)$ (§12); the block products of SM.10(c) are not measured at large $\sigma$ |
| derivative of `log`, $\ominus$ | $\theta(X) < \pi$ (LG.2(c), LG.14) | as §12's missing row (index, open items) |
| coefficients $v_0, v_1, v_2, w_j, \partial_zv_j$ | **two** masks. $\Delta = 0$ is $0/0$ for $v_1, v_2, w_2, \partial_zv_2$ (and through them $w_1 = v_2 - \sigma w_2$, $\partial_zv_1$): select the limits (SM.7(a), SM.11(b)–(c)). $\sigma = 0$ is $0/0$ for $v_0 = \mathrm{expm1}(\sigma)/\sigma$ and $w_0 = g_2(\sigma)$ at **every** $\theta$: select $1$ and $\frac12$ (SM.5(b), SM.11(b)) | a mask on $\Delta$ is needed even where the consumer error is $O(u)$, and it does not cover $v_0, w_0$: at $\sigma = 0$, $\theta = 1$, $\Delta = 1$ is far from the joint limit and both are still $0/0$. SM.7(h) is the evaluation this row summarises |

**Checked:** SM.7(g) for the onsets and the errors; arithmetic on `f64`: $\ln\mathrm{MAX} = 709.7827$, $\ln(3.4028235\times10^{38}) = 88.7228$, `exp(-745.0)` $= 5\times10^{-324}$, `exp(-745.2)` $= 0$; the closed forms against `mp.expm` at $\lvert\sigma\rvert$ up to $700$ (SM.4), the scaled ones to $709.7$ (SM.7(g)). **Permanent:** none specified: the strata stop at $\lvert\sigma\rvert = 3$.
