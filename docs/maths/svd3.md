# The 3×3 SVD and the nearest rotation

> **Non-normative.** [`NUMERICS.md`](../NUMERICS.md) §18–§19 are the contract; this page derives
> why they hold. Records [`0101`](../decisions/0101-svd3-is-one-sided-jacobi-on-eig3-s-rotation.md)
> and [`0102`](../decisions/0102-the-nearest-rotation-takes-newton-where-it-is-well-posed.md)
> measure them. Labels: `SV` for the SVD, `PO` for the polar factor.

Notation: $A = U\Sigma V^\top$, $B = AV$ the working columns $b_j$, $\alpha_j = \lVert b_j\rVert^2$,
$\gamma_{pq} = b_p\cdot b_q$; $u$ the unit roundoff. **Signed convention**: $U, V \in \mathrm{SO}(3)$,
$\sigma_1 \ge \sigma_2 \ge \lvert\sigma_3\rvert$, $\operatorname{sign}\sigma_3 = \operatorname{sign}\det A$;
$\tilde{\sigma} = (\sigma_1, \sigma_2, \sigma_3)$ in it equals $(\sigma_1, \sigma_2, d\,\lvert\sigma_3\rvert)$,
$d = \operatorname{sign}\det A$, in SO.13(c)'s.

## 1. One-sided Jacobi

**SV.1 (the rotation).** For columns $b_p, b_q$ with $d = \beta - \alpha$, $\beta = \lVert b_q\rVert^2$,
$\alpha = \lVert b_p\rVert^2$, $\gamma = b_p\cdot b_q$, `NUMERICS.md` §17's rotation of the symmetric
block $\begin{bmatrix}\alpha & \gamma\\ \gamma & \beta\end{bmatrix}$, $w = 2\gamma\operatorname{copysign}(1, d)$,
$h = \sqrt{d^2 + w^2}$, $t = w/(\lvert d\rvert + h)$, $c^2 = (\lvert d\rvert + h)/2h$, $s = ct$, makes
$(c\,b_p - s\,b_q)\cdot(s\,b_p + c\,b_q) = 0$, with $\lvert t\rvert \le 1$. No step divides by
$\gamma$; at $\gamma = 0$, $t = 0$ and $\dot t = \dot w/(\lvert d\rvert + h)$.

*Proof.* The new inner product is $cs(\alpha - \beta) + (c^2 - s^2)\gamma = c^2(\gamma(1 - t^2) - td)$.
$t = \tan\theta$ with $\tan 2\theta = 2\gamma/d$ is the root of $\gamma t^2 + dt - \gamma = 0$ of
smaller magnitude, $t = 2\gamma/(d + \operatorname{sign}(d)\sqrt{d^2 + 4\gamma^2})$, which is the
stated $t$; $1 + t^2 = 2h/(\lvert d\rvert + h)$ gives $c^2$. $\square$

**SV.2 (McAdams is not backward stable).** Let $\hat V$ diagonalise $\mathrm{fl}(A^\top A)$ and
$B = A\hat V$, so $\lvert b_i\cdot b_j\rvert \lesssim u\lVert A\rVert^2$ for $i \ne j$, and let
$B = QR$ by Givens rotations with the columns sorted by norm. Then
$\lvert r_{23}\rvert \approx \lvert b_2\cdot b_3\rvert/\sigma_2$ and $\lvert r_{23}\rvert \le
\lVert b_3\rVert \le \lVert b_2\rVert \approx \sigma_2$, so keeping only $\mathrm{diag}(R)$ costs
$\min(u\lVert A\rVert^2/\sigma_2,\ \sigma_2) \le \sqrt u\,\lVert A\rVert$, attained at
$\sigma_2 \approx \sigma_3 \approx \sqrt u\lVert A\rVert$: $10^{-8}$ relative at `f64`,
$2.4\times10^{-4}$ at `f32`. Its vectors carry $A^\top A$'s conditioning, an angle up to
$u\lVert A\rVert^2/\lvert\sigma_i^2 - \sigma_j^2\rvert$: in §11's gap-weighted metric,
$\lVert A\rVert/(\sigma_i + \sigma_j)$.

*Proof.* $r_{12}, r_{13} \lesssim u\lVert A\rVert$ because $b_1$ is the largest column; $r_{23} =
q_2\cdot b_3$ with $q_2 = b_2/\lVert b_2\rVert + O(u)$; the bound is the smaller of the two
estimates. The angle is Davis–Kahan on $A^\top A$, whose gaps are $\lvert\sigma_i^2 - \sigma_j^2\rvert$. $\square$

**SV.3 (the finish keeps the columnwise bound).** If the sorted columns satisfy
$\lvert\gamma_{ij}\rvert \le \varepsilon\lVert b_i\rVert\lVert b_j\rVert$, then $u_1 = b_1/\lVert b_1\rVert$,
$u_2 = w/\lVert w\rVert$ with $w = b_2 - (u_1\cdot b_2)u_1$, $u_3 = u_1\times u_2$ and
$\sigma = (\lVert b_1\rVert, \lVert w\rVert, u_3\cdot b_3)$ satisfy
$\lVert b_j - \sigma_ju_j\rVert \le \sqrt2\,\varepsilon\lVert b_j\rVert + O(\varepsilon^2)$.

*Proof.* $b_1 - \sigma_1u_1 = 0$. $b_2 - \sigma_2u_2 = (u_1\cdot b_2)u_1$, of norm
$\lvert\gamma_{12}\rvert/\lVert b_1\rVert \le \varepsilon\lVert b_2\rVert$. $b_3 - \sigma_3u_3 =
(u_1\cdot b_3)u_1 + (u_2\cdot b_3)u_2$, each term at most $\varepsilon\lVert b_3\rVert + O(\varepsilon^2)$. $\square$

**SV.4 (the convention gives the nearest rotation).** In the signed convention $UV^\top$ is
SO.13(c)'s $R_\star = U'\mathrm{diag}(1, 1, d)V'^\top$, unique iff $\sigma_2 + \sigma_3 > 0$ (signed
$\sigma_3$). Without the ordering it is not: $\Sigma = \mathrm{diag}(1, 1, -5)$, $U = V = I$.

*Proof.* From an unsigned SVD $(U', \Sigma', V')$ the convention negates $v_3$ and $\sigma_3$ where
$\det V' < 0$, then $u_3$ and $\sigma_3$ where $\det U' < 0$, which leaves $U\Sigma V^\top$ and turns
$U'\mathrm{diag}(1, 1, d)V'^\top$ into $UV^\top$; SO.13(c) needs $\lvert\sigma_3\rvert$ to be the
smallest. $\square$

**SV.5 (the tie swaps).** With $s = \operatorname{copysign}(1, \sigma_3)$, the maps
$(u_2, u_3, v_2, v_3, \sigma_2, \sigma_3) \mapsto (s u_3, -s u_2, v_3, -v_2, s\sigma_3, s\sigma_2)$ and
$(u_1, u_2, v_1, v_2, \sigma_1, \sigma_2) \mapsto (u_2, -u_1, v_2, -v_1, \sigma_2, \sigma_1)$ keep
$U\Sigma V^\top$, $\det U$ and $\det V$, and the first makes $\sigma_2 \ge 0$ the larger magnitude.

*Proof.* $s\sigma_3\,(s u_3)v_3^\top + s\sigma_2(-s u_2)(-v_2)^\top = \sigma_3u_3v_3^\top +
\sigma_2u_2v_2^\top$, $s^2 = 1$. Each map exchanges two columns and negates one in each of $U$
and $V$ (the first negates $u_3$ for $s = 1$ or $u_2$ for $s = -1$): two sign changes per
determinant. $\square$

**SV.6 (the windows).** After a power-of-two scaling, $m = \max\lvert a_{ij}\rvert$ and
$\lVert B\rVert_F \le 3m$, so $\alpha_j \le 9m^2$ and the radicand $d^2 + w^2 \le 405m^4$. In
$[2^{-128}, 2^{128})$ (`f64`) it is below $2^{520.7}$; in $[2^{-16}, 2^{28})$ (`f32`) below
$2^{120.7}$: no overflow. At the window's floor, a radicand below the smallest normal belongs to
a block below $2^{-255}\lVert A\rVert^2$ (`f64`) and $2^{-31}\lVert A\rVert^2$ (`f32`), below
$u\lVert A\rVert^2$, so the safe arm ($c = 1$) skips only rotations that would move nothing scored.

*Proof.* Arithmetic: $(2^{-1022})^{1/2}/(2^{-128})^2$ and $(2^{-126})^{1/2}/(2^{-16})^2$. $\square$

## 2. The polar factor

**PO.1 (Newton).** For nonsingular $X_0$, $X_{k+1} = \tfrac12(\mu_kX_k + \mu_k^{-1}X_k^{-\top})$ with
$\mu_k = (\lVert X_k^{-1}\rVert_F/\lVert X_k\rVert_F)^{1/2}$ converges quadratically to the
orthogonal polar factor $Q$ of $X_0$, $\det Q = \operatorname{sign}\det X_0$; for $\det X_0 > 0$,
$Q$ is the nearest rotation. $X^{-\top} = \mathrm{cof}(X)/\det X$ needs no root, and
$\lVert X^{-1}\rVert_F = \lVert\mathrm{cof}(X)\rVert_F/\lvert\det X\rvert$. [Higham, *Functions of
Matrices*, 2008, §8.6 and Theorem 8.12.]

**PO.2 (reflect).** For $\det M < 0$ in the signed convention, $X_0 = M - 2(Mv_3)v_3^\top$ has
$\det X_0 = -\det M$ and polar factor $UV^\top$, the nearest rotation of $M$.

*Proof.* $M(I - 2v_3v_3^\top) = U\Sigma V^\top(I - 2v_3v_3^\top) = U\,\mathrm{diag}(\sigma_1, \sigma_2,
-\sigma_3)V^\top$, with $-\sigma_3 > 0$; then SV.4. $\square$

**PO.3 (the tangent).** Where $R_\star$ is unique, $\dot R_\star = U\Omega V^\top$,
$\Omega_{ij} = (C_{ij} - C_{ji})/(\sigma_i + \sigma_j)$, $C = U^\top\dot MV$, the
$\sigma_i$ signed; so $\lVert\dot R_\star\rVert_F \le 2\lVert\dot M\rVert_F/(\sigma_2 + \sigma_3)$.

*Proof.* Differentiate $R^\top M = H$ symmetric with $\dot R = R\,\Omega'$, $\Omega'$ skew:
$H\Omega' + \Omega'H = (R^\top\dot M) - (R^\top\dot M)^\top$; in the basis of $V$, $H =
V\tilde\Sigma V^\top$ and the Sylvester equation is diagonal. $\square$

**Checked (SV.1–SV.6, PO.2–PO.3):** mpmath 1.4.1, 100 digits unless stated. SV.1: 200 random column
pairs, the rotated inner product over $\sqrt{\alpha\beta}$ and $\lvert c^2 + s^2 - 1\rvert$ at most
$2.9\times10^{-101}$. SV.2: $\min(u/\sigma_2, \sigma_2)/\sqrt u$ is $1.0$ at $\sigma_2 = \sqrt u$
($u = 2^{-53}$) and $1.1\times10^{-4}$ at $10^{-4}$, $9.5\times10^{-5}$ at $10^{-12}$; a numpy
harness reads McAdams' backward error at $4.5\times10^7u$ on $\sigma = (1, 10^{-8}, \tfrac12 10^{-8})$
(`0101` *Context*). SV.3: 200 random $A$, the columns of $AV$ perturbed by $10^{-30}\lVert b_j\rVert$ per entry, $\varepsilon$
the largest $\lvert\gamma_{ij}\rvert/\lVert b_i\rVert\lVert b_j\rVert$ that leaves: worst residual
$1.41\,\varepsilon\lVert b_j\rVert$, the bound. SV.4: 60 random
$M$, no better rotation among 50 random ones each, $\lvert\det R - 1\rvert \le 1.1\times10^{-100}$.
SV.5: 100 random instances of both maps, $5.7\times10^{-101}$. SV.6: exact arithmetic. PO.2: 60
random $M$ with $\det M < 0$, $\lVert\mathrm{polar}(X_0) - R_\star\rVert$ and
$\lvert\det X_0 + \det M\rvert$ at most $9.1\times10^{-100}$. PO.3: 60 digits, 20 random $M$,
central difference at $h = 10^{-25}$: $1.3\times10^{-35}$ relative. Scripts not committed.
**Permanent:** corpus ids `svd3` and `nearest_rotation` (`0101`, `0102`); `svd3_tests` and
`polar_tests`.
