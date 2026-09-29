# SE_N(3): Exp, Log, adjoints, the Q block, the dual-matrix algebra

> Non-normative companion to [`NUMERICS.md`](../NUMERICS.md) §1 (SE_N(3), tangents), §2.2, §2.4, §5
> and §14 (`SEn3Jac`, `jr_inv`, `jl`). **`NUMERICS.md` wins on any conflict; a disagreement is an
> open item in the [maths index](./index.md)**, which holds the notation and the `Checked:`
> convention. The general algebra ($\mathrm{Ad}$, $\mathrm{ad}$, $J = g(\mathrm{ad})$, the sides
> table) is in [`lie-groups.md`](./lie-groups.md), the SO(3) closed forms in [`so3.md`](./so3.md),
> the series and rounding of $b, d, e$ in [`coefficients.md`](./coefficients.md). This page derives
> what those take as given: the block $Q$ of §5.3, the dual-matrix algebra of §2.2 and the order
> conversion of [`0002`](../decisions/0002-one-convention-for-a-stack-that-already-disagrees.md).

On this page $N \ge 1$, $\tau = [\varphi;\rho_1;\dots;\rho_N]$, $W = [\varphi]_\times$, $z = \theta^2$, and
$\Xi = [\rho]_\times = \rho^\wedge$ for one translation block $\rho$ (the index's $P$ is written $\Pi$
here). The tangent dimension $3 + 3N$ is written out: $n$ is a power (of $\mathrm{ad}$, LG.7), $\mu$ a
vector (SE.7), $u$ the unit roundoff and $s$ a summation index (SE.6) or $s(\theta) = 2\sin(\theta/2)/\theta$
(SE.16, LG.16). Block $(i,0)$ is the $i$-th $3\times3$ block of the first block
column. $\sigma_m$ is CO.1's series; CO.1's $\tau_m$ is written $T_m$, because $\tau$ is the tangent.
Norms are spectral ($\lVert W\rVert = \theta$, $\lVert\Xi\rVert = \lVert\rho\rVert$); "max" in a
`Checked:` line is the largest absolute entry.

## Results

| Label | Result | `NUMERICS.md` | Implemented by |
|---|---|---|---|
| SE.1–SE.2 | $\mathrm{SE}_N(3)$ as $(3+N)$-square matrices; $XY = (R_XR_Y,\ R_Xy_i + x_i)$, $X^{-1} = (R^\top, -R^\top x_i)$; $N=1$ SE(3), $N=2$ SE₂(3) with $x_1 = v$, $x_2 = p$ (a stated convention) | §1 | `SEn3`, `SE3`, `SE23`, `Mul`, `inverse` |
| SE.3 | $\mathrm{Exp}\,\tau = (\mathrm{Exp}\,\varphi,\ J_l(\varphi)\rho_i)$; $\mathrm{Log}\,X = (\mathrm{Log}\,R,\ J_l^{-1}(\varphi)x_i)$ for $\theta(X) \le \pi$, inverting $\mathrm{Exp}$ for $\theta < \pi$ | §5.1 | `exp`, `log` |
| SE.4 | $\mathrm{Ad}_X = R + \epsilon[x_i]_\times R$, $\mathrm{ad}_\tau = W + \epsilon[\rho_i]_\times$ | §5.2 | `adjoint`, `ad` |
| SE.5–SE.8 | $J_l = \sum\mathrm{ad}^n/(n+1)!$: dual powers of $\mathrm{ad}$, the words $W^j\Xi W^k$, $W^3 = -zW$, the reductions among them; $Q(\rho,\varphi)$ with $b, d, e$ | §5.3, §4 | `q_coeffs`, `jl`, `jr` |
| SE.9 | $J_r$ has $Q(-\rho,-\varphi)$ (expanded); $Q_l = [J_l\rho]_\times J_l + e^WQ_r$; $\partial_\varphi(J_l\rho) = Q - [J_l\rho]_\times J_l$; an equal five-term form | §5.3, §1 | `jr`, `jl` |
| SE.10–SE.11 | dual matrices $A + \epsilon B$: closure, product ($27 + 54N$ multiplications), inverse, `apply`, `apply_transpose`, determinant | §2.2, §14 | `SEn3Jac::{mul, inverse, apply, apply_transpose, write_dense}` |
| SE.12 | $J_r^{-1}$, $J_l^{-1}$ are the dual inverse of $J_r$, $J_l$: no second closed form | §5.4, §12 | `jr_inv`, `jl_inv` |
| SE.13 | SE(3) action Jacobians, both sides, and their relation through $\mathrm{Ad}$ | §2.4, §5.5 | `act_jacobians` |
| SE.14 | rotation-first vs translation-first: $\Pi$, $\Pi M\Pi^\top$, covariances, gradients | §1, `0002` | `Twist::{from_translation_first, to_translation_first}` |
| SE.15 | what the exact arms of $b, d, e$ cost `Exp` and $Q$: one power of $\theta$, not the coefficient's; the error at the switch | §4, §5.1, §5.3 | `exp_coeffs`, `q_coeffs` |
| SE.16 | rounding on the inverse side: `Log`'s translation part ($\lVert \Lambda_i\rVert \le 1.29\lVert x_i\rVert$), the blocks of $J^{-1}$ (the scale $\lVert Q\rVert/s^2$, $1/s$ toward $2\pi$) | §5.1, §5.4, §11 | `log`, `jr_inv`, `jl_inv` |

## 1. The group

**Definition SE.1.** $G = \mathrm{SE}_N(3)$ is the set of $(3+N)$-square matrices
$X = \begin{bmatrix}R & x_1\cdots x_N\\ 0 & I_N\end{bmatrix} = (R, x_1,\dots,x_N)$, $R \in \mathrm{SO}(3)$,
$x_i \in \mathbb R^3$ (`NUMERICS.md` §1). $N = 1$ is $\mathrm{SE}(3)$, $x_1 = t$, acting by
$Xp = Rp + t$, i.e. on $[p;1]$. $N = 2$ is $\mathrm{SE}_2(3)$, and **by convention $x_1 = v$ (velocity),
$x_2 = p$ (position)**,
$X = \begin{bmatrix}R & v & p\\ 0 & 1 & 0\\ 0 & 0 & 1\end{bmatrix}$ with tangent $[\varphi;\rho_v;\rho_p]$
(the invariant-EKF convention of `PHASE3.md` §5). The convention fixes labels and
the order of blocks in dense output only: the group law, $\mathrm{Exp}$ and every Jacobian below
treat the $x_i$ alike, and no result depends on which slot carries which quantity.

**Proposition SE.2.** $XY = (R_XR_Y,\ R_Xy_i + x_i)$ and $X^{-1} = (R^\top, -R^\top x_i)$. So $G$ is
closed under both, $G \cong \mathrm{SO}(3)\ltimes(\mathbb R^3)^N$ with $R$ acting on every slot, and
$X \mapsto (R, x_i)$ is a homomorphism onto $\mathrm{SE}(3)$ for each $i$: the columns are $N$ copies
of one rigid motion sharing the rotation.

*Proof.* Block multiplication,
$\begin{bmatrix}R & x\\ 0 & I\end{bmatrix}\begin{bmatrix}S & y\\ 0 & I\end{bmatrix} = \begin{bmatrix}RS & Ry + x\\ 0 & I\end{bmatrix}$,
and the stated inverse multiplies out to $I$. $\square$

**Checked:** mpmath 1.4.1, 110 digits, $N = 1, 2, 3$, 30 random pairs ($R, S$ from $N(0,1)$ rotation
vectors, $x_i, y_i \sim N(0, 2^2)$): the matrix product and `mp.inverse` against the block formulas,
max $1.3\times10^{-110}$; the last $N$ block rows stay $(0, I)$ exactly. Script not committed.
**Permanent:** planned, proptest `group_axioms_*` (`PHASE3.md` §9).

## 2. Exp and Log

**Proposition SE.3.** Let $\tau = [\varphi;\rho_i]$, $A = \tau^\wedge$.

- (a) $A^k = \begin{bmatrix}W^k & W^{k-1}\rho_1\cdots W^{k-1}\rho_N\\ 0 & 0\end{bmatrix}$ for $k \ge 1$.
- (b) $\mathrm{Exp}\,\tau = (e^W,\ \Gamma_1(\varphi)\rho_i)$, $\Gamma_1 = \sum_{k\ge0}W^k/(k+1)! = J_l(\varphi)$ (SO.7).
- (c) For $\theta(X) \le \pi$: $\mathrm{Log}(X) = (\varphi = \mathrm{Log}\,R,\ \rho_i = J_l^{-1}(\varphi)x_i)$ with
  $J_l^{-1} = I - \tfrac12W + cW^2$ (SO.8), and $\mathrm{Exp}(\mathrm{Log}\,X) = X$.
- (d) For $\theta(\tau) < \pi$: $\mathrm{Log}(\mathrm{Exp}\,\tau) = \tau$.

*Proof.* (a) $A^{k+1} = A\,A^k$ has top-left block $W\cdot W^k$, top-right $W\cdot W^{k-1}\rho$ (the
block $\rho\cdot0$ vanishes) and zero bottom rows. (b) Sum $I + \sum_{k\ge1}A^k/k!$: top-left $e^W$, top-right
$\sum_{k\ge1}W^{k-1}\rho_i/k! = \Gamma_1\rho_i$. (c) $\varphi = \mathrm{Log}\,R$ has $\theta \le \pi$ and
$e^W = R$ (by $\mathrm{atan2}$, SO.5); $\Gamma_1(\varphi)$ is invertible ($\det = 2(1-\cos\theta)/\theta^2$,
LG.2(b), and $\theta \le \pi < 2\pi$) with $\Gamma_1^{-1} = J_l^{-1}$ (SO.8), so (b) gives
$\mathrm{Exp}(\varphi, \rho_i) = X$. (d) For $\theta(\tau) < \pi$, $\varphi$ is the unique preimage of
$e^W$ in $\theta < \pi$ (LG.2(c)) and $\Gamma_1\rho_i = x_i$ has the unique solution $\rho_i$. $\square$

*Consequences.* Every translation column depends on its own $\rho_i$ through the one matrix
$J_l(\varphi)$: $\mathrm{Exp}$ costs one $\mathrm{SO}(3)$ exponential, one $J_l$ and $N$ matrix-vector
products. For $\mathrm{SE}_2(3)$ both columns use $\Gamma_1$; $\Gamma_2$ (`NUMERICS.md` §7) is not in
$\mathrm{Exp}$ but in the integrated quantities of a piecewise-constant $(\omega, a)$. At $\theta(X) = \pi$ both
signs of $\varphi$ are logarithms and (c) returns a valid preimage for either, but $\rho_i' - \rho_i = \varphi\times x_i$
(LG.2(d), with the jump of $\mathrm{Log}$): which is returned follows the sign of $\varphi$, the SO(3) `Log`'s tie at
$w = +0$ (`NUMERICS.md` §3.2; index, open items), so $\rho_i$ there is a function of the quaternion, not of $X$.
The rounding of the inverse side is SE.16.

**Checked:** mpmath 1.4.1, 110 digits, $N = 1, 2, 3$, 40 random $\tau$ ($\theta \sim U[0,6]$, $\rho_i \sim N(0, 2^2)$) plus edge
$\tau$, $\theta \in \{0, 10^{-10}, \pi\pm10^{-3}, 2\pi - 10^{-3}, 6.5, 20\}$, plus two with $\rho_i \sim N(0, (10^6)^2)$ at $\theta = 0.3, 2.9$;
errors relative to $\max(1, \lVert\cdot\rVert_{\max})$ of the exact quantity. (a) all powers $k \le 6$, $3.1\times10^{-111}$. (b) `mp.expm` of
the hat matrix against the blocks with $\Gamma_1$ summed from its series, $5.5\times10^{-105}$ (at $\theta = 20$; $4.2\times10^{-111}$ for
$\rho \sim 10^6$); $\det\Gamma_1$ against $2(1-\cos\theta)/\theta^2$ **evaluated at 250 digits** ($\theta = 0$: the limit $1$), relative
$2.5\times10^{-104}$ (at $\theta = 20$); at the working 110 digits the figure is LG.2's $1.2\times10^{-92}$, all of it the closed form's
cancellation at $\theta = 10^{-10}$. (c) $\mathrm{Log}(\mathrm{Exp}\,\tau) = \tau$ for $0 < \theta < 3.1$ (geometric $\mathrm{SO}(3)$ logarithm,
$\Gamma_1$ solved from its series), $8.3\times10^{-109}$ (the $\rho \sim 10^6$ pair, relative; $1.2\times10^{-109}$ otherwise); the closed
$J_l^{-1}$ ($c$ at 250 digits) against `mp.inverse`, $8.3\times10^{-111}$. `mp.logm` of the $(3+N)$-square matrix agrees with $\tau^\wedge$
for $\theta \le 3.02$ (60 digits, 84 $\tau$, $6.4\times10^{-61}$) and returns a complex, non-principal result from an onset in
$(3.02, 3.04]$ (grid $0.005$ on $[2.9, 3.1]$, 8 axes per $N$, 40 digits; the same onsets in mpmath 1.3.0 and 1.4.1), as for
$\mathrm{SO}(3)$ (SO.5): `sen3_log_n*` (`PHASE1.md` §4.3) shares the fault of `so3_log` (index, open items). Scripts not committed.
**Permanent:** planned, corpus `sen3_exp_n{1,2,3}`, `sen3_log_n{1,2,3}`; proptest `exp_log_roundtrip_*`.

## 3. Adjoints

**Proposition SE.4.** For all $X = (R, x_i)$, $\tau = [\varphi;\rho_i]$ and $\psi = [\alpha;\eta_i]$ (no domain):
$$
\mathrm{Ad}_X\psi = \big[R\alpha;\ R\eta_i + [x_i]_\times R\alpha\big],\qquad
\mathrm{ad}_\tau\psi = \big[\varphi\times\alpha;\ \varphi\times\eta_i + \rho_i\times\alpha\big],
$$
i.e. $\mathrm{Ad}_X = R + \epsilon[x_i]_\times R$ and $\mathrm{ad}_\tau = W + \epsilon[\rho_i]_\times$
(`NUMERICS.md` §5.2), and $X\,\mathrm{Exp}(\psi)X^{-1} = \mathrm{Exp}(\mathrm{Ad}_X\psi)$.

*Proof.* With $X^{-1} = (R^\top, -R^\top x_i)$, $X\psi^\wedge X^{-1}$ has top-left $R\alpha^\wedge R^\top =
(R\alpha)^\wedge$ and column $i$ of its top-right block $R\eta_i - R\alpha^\wedge R^\top x_i = R\eta_i -
(R\alpha)\times x_i = R\eta_i + [x_i]_\times R\alpha$; the bottom rows vanish. The commutator
$[\tau^\wedge,\psi^\wedge]$ has top-left $[W, \alpha^\wedge] = (\varphi\times\alpha)^\wedge$ and column $i$
$W\eta_i - \alpha^\wedge\rho_i = \varphi\times\eta_i + \rho_i\times\alpha$. The last claim is LG.4(a). $\square$

The sign of $[x_i]_\times R$ comes from $-(R\alpha)\times x_i = +[x_i]_\times R\alpha$; both blocks sit
**below** the diagonal because the rotation comes first (LG §6, "block position").

**Checked:** mpmath 1.4.1, 110 digits, $N = 1, 2, 3$: $\mathrm{Ad}_X$ and $\mathrm{ad}_\tau$ built from
their definitions (images of the basis under $\psi \mapsto (X\psi^\wedge X^{-1})^\vee$, $[\tau^\wedge,\cdot]^\vee$;
25 random $X$, rotation vector $\sim N(0, 0.6^2I)$, $x_i \sim N(0,2^2)$; the first 20 $\tau$ of SE.3):
images stay in $\mathfrak g$ (bottom rows exactly $0$), $\mathrm{Ad}$ against the blocks $6.7\times10^{-111}$,
$\mathrm{ad}$ exact, $X\,\mathrm{Exp}(\psi)X^{-1} = \mathrm{Exp}(\mathrm{Ad}_X\psi)$ $2.0\times10^{-110}$. Script not
committed. **Permanent:** planned, corpus `sen3_ad_n{1,2,3}`; proptest `adjoint_identity_*`.

## 4. The block Q

By LG.7–LG.8, $J_l(\tau) = \sum_{n\ge0}\mathrm{ad}_\tau^n/(n+1)!$. The diagonal blocks of $\mathrm{ad}^n$
are $W^n$, so those of $J_l$ are $J_l(\varphi)$ (LG.10); what is left is the first block column.

**Lemma SE.5 (dual powers of ad).** For every $\tau$ and $n \ge 1$, block $(i,0)$ of $\mathrm{ad}_\tau^n$ is
$\sum_{j+k=n-1}W^j\,\rho_i^\wedge\,W^k$, all other off-diagonal blocks vanish, and therefore
$J_l(\tau) = J_l(\varphi) + \epsilon\,Q(\rho_i,\varphi)$ with
$$
Q(\rho,\varphi) = \sum_{n\ge1}\frac{1}{(n+1)!}\sum_{j+k=n-1}W^j\,\Xi\,W^k
= \sum_{j,k\ge0}\frac{W^j\,\Xi\,W^k}{(j+k+2)!}. \tag{1}
$$
The series converges absolutely (each word has norm $\le\theta^{j+k}\lVert\rho\rVert$).

*Proof.* Block $(i,0)$ of a product of matrices of this shape is $M_{i0}M'_{00} + M_{ii}M'_{i0}$ and the
blocks $(i,j)$, $j\ne0$, $i\ne j$, stay $0$ (SE.11 states the general rule). So $\mathrm{ad}^{n+1} =
\mathrm{ad}^n\,\mathrm{ad}$ gives the column recursion $B_{n+1} = B_nW + W^n\Xi$, $B_1 = \Xi$, whose
solution is $B_n = \sum_{j+k=n-1}W^j\Xi W^k$ ($B_2 = \Xi W + W\Xi$, $B_3 = \Xi W^2 + W\Xi W + W^2\Xi$).
Sum with the weights $1/(n+1)!$. $\square$

**Lemma SE.6 (Cayley–Hamilton reduction).** With $W^3 = -zW$ (SO.3), $W^{2m+1} = (-z)^mW$ and
$W^{2m+2} = (-z)^mW^2$. Sorting the words of (1) by the class of $j$ and of $k$ (none: $0$;
odd; even $\ge2$) gives exactly
$$
Q = \tfrac12\Xi + b\,(W\Xi + \Xi W) + d\,(W^2\Xi + \Xi W^2) + T_4\,W\Xi W
+ e\,(W^2\Xi W + W\Xi W^2) + T_6\,W^2\Xi W^2, \tag{2}
$$
with $b = \sigma_3(z)$, $d = \sigma_4(z)$, $e = T_5(z)$ (CO.2) and $T_M = \sum_s(s+1)(-z)^s/(2s+M)!$.

*Proof.* Each class pair contributes one word with a scalar series:

| $j$, $k$ | word | $j+k+2$ | scalar |
|---|---|---|---|
| $0$, $0$ | $\Xi$ | $2$ | $\tfrac12$ |
| $2m+1$, $0$ | $W\Xi$ | $2m+3$ | $\sum_m(-z)^m/(2m+3)! = \sigma_3 = b$ |
| $2m+2$, $0$ | $W^2\Xi$ | $2m+4$ | $\sigma_4 = d$ |
| $2m+1$, $2l+1$ | $W\Xi W$ | $2(m+l)+4$ | $T_4$ |
| $2m+2$, $2l+1$ | $W^2\Xi W$ | $2(m+l)+5$ | $T_5 = e$ |
| $2m+2$, $2l+2$ | $W^2\Xi W^2$ | $2(m+l)+6$ | $T_6$ |

Here $(-z)^m(-z)^l$ has the exponent $s = m+l$ and there are $s+1$ pairs $(m,l)$ with that sum, which is
the factor $(s+1)$ of $T_M$. The pairs with $j$, $k$ exchanged give the mirrored words with the same
scalars. $\square$

**Lemma SE.7 (identities among the words).** $[\alpha]_\times[\beta]_\times = \beta\alpha^\top - (\alpha\cdot\beta)I$;
hence

- (a) $W\Xi W = -(\varphi\cdot\rho)\,W$;
- (b) $W^2\Xi W = W\Xi W^2 = -(\varphi\cdot\rho)\,W^2$ and $W^2\Xi W^2 = -z\,W\Xi W$;
- (c) $W^2\Xi + \Xi W^2 = -z\,\Xi + W\Xi W$.

*Proof.* The first identity is $\alpha\times(\beta\times v) = \beta(\alpha\cdot v) - v(\alpha\cdot\beta)$. (a)
$W\Xi W = (\rho\varphi^\top - (\varphi\cdot\rho)I)W = -(\varphi\cdot\rho)W$, since $\varphi^\top W = -(W\varphi)^\top = 0$.
(b) $W^2\Xi W = W(W\Xi W)$, $W\Xi W^2 = (W\Xi W)W$, and $W^2\Xi W^2 = W(W\Xi W)W = -(\varphi\cdot\rho)W^3 =
z(\varphi\cdot\rho)W = -zW\Xi W$. (c) $W^2 = \varphi\varphi^\top - zI$, so $W^2\Xi + \Xi W^2 =
\varphi\varphi^\top\Xi + \Xi\varphi\varphi^\top - 2z\Xi$, and with $\mu = \varphi\times\rho$,
$\varphi\varphi^\top\Xi = \varphi\mu^\top$, $\Xi\varphi\varphi^\top = -\mu\varphi^\top$, while
$\varphi\mu^\top - \mu\varphi^\top = -[\varphi\times\mu]_\times = -(\varphi\cdot\rho)W + z\Xi$; add (a). $\square$

**Proposition SE.8 (the block Q, `NUMERICS.md` §5.3).** For every $\rho$, $\varphi$ (any $\theta$), with $b, d, e$ of
`NUMERICS.md` §4,
$$
\begin{aligned}
Q(\rho,\varphi) ={}& \tfrac12\rho^\wedge
+ b\,\big(\varphi^\wedge\rho^\wedge + \rho^\wedge\varphi^\wedge + \varphi^\wedge\rho^\wedge\varphi^\wedge\big)
+ d\,\big(\varphi^\wedge\varphi^\wedge\rho^\wedge + \rho^\wedge\varphi^\wedge\varphi^\wedge - 3\,\varphi^\wedge\rho^\wedge\varphi^\wedge\big) \\
&+ e\,\big(\varphi^\wedge\rho^\wedge\varphi^\wedge\varphi^\wedge + \varphi^\wedge\varphi^\wedge\rho^\wedge\varphi^\wedge\big),
\end{aligned}
$$
and $J_l(\tau) = J_l(\varphi) + \epsilon\,Q(\rho_i,\varphi)$, in this repository's rotation-first order
(lower-left). Each block depends on its own $\rho_i$ only.

*Proof.* In (2) replace $W^2\Xi W^2$ by $-zW\Xi W$ (SE.7(b)): the coefficient of $W\Xi W$ becomes
$T_4 - zT_6$. By CO.4(b), $T_M = \tfrac12(\sigma_{M-1} - (M-2)\sigma_M)$, and by CO.4(a),
$z\sigma_M = \tfrac1{(M-2)!} - \sigma_{M-2}$. So $T_4 = \tfrac12(b - 2d)$ and
$zT_6 = \tfrac12(z\sigma_5 - 4z\sigma_6) = \tfrac12\big(\tfrac16 - b - \tfrac4{24} + 4d\big) = 2d - \tfrac b2$,
whence $T_4 - zT_6 = b - 3d$: the $b$-group carries $W\Xi W$ once and the $d$-group $-3$ times. The
remaining words and scalars of (2) are those of §5.3. $\square$

*Remarks.* (i) The nine words $W^j\Xi W^k$, $j,k\in\{0,1,2\}$, span only six dimensions (three relations,
SE.7; rank $6$ on 30 random $(\varphi,\rho)$), so the arrangement of §5.3 is one representative and (2)
the unreduced one; the corpus compares values, not arrangements. (ii) SE.7(a)–(c) turn §5.3 into the equal
$$
Q = a\,\Xi + b\,(W\Xi + \Xi W) - (\varphi\cdot\rho)\big[(b - 2d)\,W + 2e\,W^2\big],\qquad a = \tfrac12 - zd
$$
(CO.4(a)), five words in place of nine. It is *not* §5.3 and is not adopted: another arrangement is a §5.3
edit and a record.

**Checked:** mpmath 1.4.1, 110 digits; $b, d, e$ from their definitions with 100 guard digits. 200
random $(\varphi,\rho)$ ($\theta \sim U[0,6]$, $\rho_i \sim N(0,2^2)$) plus the 7 edge $\theta$ of SE.3 and its two
$\rho_i \sim 10^6$ cases, the words in $3\times3$: §5.3 against (1) summed by the recursion of SE.5,
$9.9\times10^{-105}$ (at $\theta = 20$; relative to $\max(1,\lVert Q\rVert)$; $7.8\times10^{-111}$ for
$\rho \sim 10^6$); (2) $1.0\times10^{-110}$; the double sum (1) truncated at $j+k\le90$ ($\theta\in[0.1,2]$,
60 random and one $\rho \sim 10^6$), $3.4\times10^{-109}$; the identities SE.7 $\le 2.6\times10^{-109}$ and exactly,
as polynomial identities in the six symbols $\varphi_i,\rho_i$ (sympy); $\sigma_3 = b$, $\sigma_4 = d$,
$T_5 = e$ $\le 1.1\times10^{-112}$; $T_4 - zT_6 = b - 3d$ $2.1\times10^{-112}$ and, in exact rationals,
through $z^{14}$ (sympy 1.14); the Taylor series of the closed forms of $b, d, e$ to $\theta^{13}$ equal those of
$\sigma_3$, $\sigma_4$, $T_5$; remark (ii) $8.6\times10^{-111}$. Against the **dense** series
$\sum(\pm\mathrm{ad}_\tau)^n/(n+1)!$ of the $3(N+1)$-square matrix (independent of SE.5), $N = 1, 2, 3$, 25
random $\tau$ ($\theta < 6$) plus the 9 edge: diagonal blocks equal $J_{l,r}(\varphi)$, blocks $(i,0)$ equal
$Q(\rho_i,\varphi)$ ($J_l$) and $Q(-\rho_i,-\varphi)$ ($J_r$) to $8.9\times10^{-105}$ (at $\theta = 20$;
$7.8\times10^{-111}$ for $\rho \sim 10^6$), every other block $0$.
Against the derivative of the operation itself, $N = 1, 2, 3$, 6 random $\tau$ ($\theta<2.9$) plus one at
$\theta = 10^{-10}$ and the two $\rho_i \sim 10^6$ cases: central differences $h = 10^{-30}$ of
$\mathrm{Log}(\mathrm{Exp}(\tau+h e_j)\mathrm{Exp}(\tau)^{-1})$
and $\mathrm{Log}(\mathrm{Exp}(\tau)^{-1}\mathrm{Exp}(\tau+h e_j))$ against the full $J_l$, $J_r$ with $Q$:
$3.4\times10^{-62}$ (relative to the largest entry, the differencing; $3.1\times10^{-62}$ for $\rho \sim 10^6$).
Scripts not committed. **Permanent:** planned,
corpus `sen3_jl_n{1,2,3}`, `sen3_jr_n{1,2,3}` (defined by exactly the dense series, `PHASE1.md` §4.3),
`coeff_b`, `coeff_d`, `coeff_e`; proptest `jacobians_match_dual_*`. The seeded defect "$Q$ with
$-\tfrac12\rho^\wedge$" (`PHASE1.md` §10) is the $n = 1$ term of (1) with the wrong sign.

**Corollary SE.9 ($J_r$, and two relations; every $\tau$).**

- (a) $J_r(\tau) = J_l(-\tau) = J_r(\varphi) + \epsilon\,Q(-\rho_i,-\varphi)$: a word with $m$ factors
  ($m-1$ copies of $W$ and one $\Xi$) changes sign by $(-1)^m$, so
  $Q(-\rho,-\varphi) = -\tfrac12\Xi + b(W\Xi + \Xi W - W\Xi W) - d(W^2\Xi + \Xi W^2 - 3W\Xi W)
  + e(W\Xi W^2 + W^2\Xi W)$. The **whole** tangent is negated.
- (b) $Q(\rho,\varphi) = [J_l(\varphi)\rho]_\times J_l(\varphi) + e^W\,Q(-\rho,-\varphi)$.
- (c) $\partial_\varphi\big(J_l(\varphi)\rho\big) = Q(\rho,\varphi) - [J_l(\varphi)\rho]_\times J_l(\varphi)$.

*Proof.* (a) $J_r(\tau) = J_l(-\tau)$ (LG.9(a)) and $\mathrm{ad}_{-\tau} = -\mathrm{ad}_\tau$. (b) LG.9(b),
$J_l = \mathrm{Ad}_EJ_r$ with $E = \mathrm{Exp}\,\tau$, $E = (e^W, x_i)$, $x_i = J_l(\varphi)\rho_i$ (SE.3), read in block
$(i,0)$ with the product rule of SE.11: $[x_i]_\times e^WJ_r(\varphi) + e^WQ_r$, and $e^WJ_r(\varphi) = J_l(\varphi)$.
(c) By LG.8, $\mathrm{Exp}(\tau+\delta) = \mathrm{Exp}(J_l\delta)\mathrm{Exp}(\tau) + O(\lVert\delta\rVert^2)$. Take
$\delta = [\delta\varphi;0]$: then $J_l\delta = [J_l(\varphi)\delta\varphi;\ \dots;\ Q\delta\varphi;\ \dots]$, and the translation
column $i$ of $\mathrm{Exp}(\eta)E$ is $x_i + \eta_\varphi\times x_i + \eta_{\rho_i} + O$. So
$\Gamma_1(\varphi+\delta\varphi)\rho_i = x_i - [x_i]_\times J_l\delta\varphi + Q\delta\varphi + O$. $\square$

(c) matters for tests: forward-mode differentiation of the *translation part* of `exp` gives the left side, not
the block of $J_l$; the difference is the $\mathrm{Ad}$-type term $[x_i]_\times J_l(\varphi)$.

**Checked:** as SE.8, with its two $\rho\sim10^6$ cases: (a) the expansion, $1.0\times10^{-110}$, and §5.3 at
$(-\rho,-\varphi)$ against the dense series, $8.9\times10^{-105}$ (SE.8); (b) $9.8\times10^{-105}$ (at $\theta = 20$;
$6.5\times10^{-111}$ for $\rho\sim10^6$); (c) 15 random $(\rho,\varphi)$, $\theta\in[0.05,5.5]$, and the two
$\rho\sim10^6$ cases, central differences of $\Gamma_1\rho$ (series), $h = 10^{-30}$, $9.4\times10^{-62}$
(relative to the largest entry, the differencing). Scripts not committed. **Permanent:**
proptest `jl_is_ad_jr_*` for (b); (c) none specified.

## 5. The dual-matrix algebra

**Definition SE.10.** $\mathcal E_N = \mathbb R[\epsilon_1,\dots,\epsilon_N]/(\epsilon_i\epsilon_j)$: the
commutative algebra with basis $1,\epsilon_1,\dots,\epsilon_N$ and *all* products $\epsilon_i\epsilon_j = 0$.
$\mathcal D_N = M_3(\mathcal E_N)$ is the set of $A + \sum_i\epsilon_iB_i$, $A, B_i \in \mathbb R^{3\times3}$, written
$A + \epsilon B$ (`NUMERICS.md` §2.2). Its **dense form** $\iota(A+\epsilon B)$ is the matrix of
$v\mapsto(A+\epsilon B)v$ on $\mathcal E_N^3 \cong \mathbb R^{3(N+1)}$ in the basis $(1,\epsilon_1,\dots,\epsilon_N)$:
with $v = v_0 + \sum_i\epsilon_iv_i$,
$$
(A+\epsilon B)v = Av_0 + \sum_i\epsilon_i\,(B_iv_0 + Av_i)
\quad\Longleftrightarrow\quad
\iota(A+\epsilon B) = \begin{bmatrix}A & & & \\ B_1 & A & & \\ \vdots & & \ddots & \\ B_N & & & A\end{bmatrix}.
$$

**Proposition SE.11.**

- (a) *Closure and product.* $\mathcal D_N$ is an associative unital algebra with
  $(A+\epsilon B)(C+\epsilon D) = AC + \epsilon\,(B_iC + AD_i)$, and $\iota$ is an injective algebra
  homomorphism onto the matrices of the displayed shape: **the dense product is the dual product**.
- (b) *Determinant and inverse.* $\det\iota(A+\epsilon B) = (\det A)^{N+1}$. If $\det A \ne 0$ then
  $(A+\epsilon B)^{-1} = A^{-1} - \epsilon\,A^{-1}B_iA^{-1}$, and it is the inverse of the dense matrix; if
  $\det A = 0$ there is no inverse.
- (c) *Action.* $\iota(A+\epsilon B)[\psi;\eta_i] = [A\psi;\ B_i\psi + A\eta_i]$ and
  $\iota(A+\epsilon B)^\top[\psi;\eta_i] = [A^\top\psi + \sum_iB_i^\top\eta_i;\ A^\top\eta_i]$. The transpose is block
  *upper* triangular, so it leaves $\mathcal D_N$: `apply_transpose` and `sandwich` are dense-shaped
  routines, not algebra operations.
- (d) *Cost.* The product takes $27 + 54N$ scalar multiplications ($AC$; $B_iC$ and $AD_i$ for each $i$) and
  $18 + 45N$ additions; the inverse $54N$ beyond $A^{-1}$; `apply` and `apply_transpose` $9 + 18N$.

| $N$ | dual $27+54N$ | generic block-lower-triangular $\tfrac{27}6(N+1)(N+2)(N+3)$ | full dense $27(N+1)^3$ |
|---|---|---|---|
| 1 | 81 | 108 | 216 |
| 2 | 135 | 270 | 729 |
| 3 | 189 | 540 | 1728 |

- (e) *Contents.* $\mathrm{Ad}_X$, $\mathrm{Ad}_X^{-1}$, $\mathrm{ad}_\tau$, $J_r$, $J_l$, $J_r^{-1}$, $J_l^{-1}$ lie in $\mathcal D_N$, with
  diagonal block the $\mathrm{SO}(3)$ counterpart (LG.10).

*Proof.* (a) Matrices over the commutative ring $\mathcal E_N$ multiply as
$(A + \sum\epsilon_iB_i)(C + \sum\epsilon_jD_j) = AC + \sum_i\epsilon_i(B_iC + AD_i) + \sum_{i,j}\epsilon_i\epsilon_jB_iD_j$,
and the last sum is $0$. $\iota$ is the action on the free module $\mathcal E_N^3$; it is a homomorphism
because that action is, injective because $v = e_c$ (a unit vector) returns column $c$ of $A$ and of every
$B_i$, and its image is read off the formula. (b) $\iota$ is block lower triangular with $N+1$ diagonal
blocks $A$. Multiplying out, $(A+\epsilon B)(A^{-1} - \epsilon A^{-1}B_iA^{-1}) = I + \epsilon(B_iA^{-1} - B_iA^{-1}) = I$
and likewise on the other side; $\iota$ is injective, so the dense inverse is unique and equal to it.
(c) The first is the formula above; the second is the transpose of the dense matrix. (d) Count the
$3\times3$ products (27 multiplications, 18 additions each, plus 9 for each block sum): $1 + 2N$ of them. A
generic block-lower-triangular product forms block $(i,j)$, $i \ge j$, from $i-j+1$ block products, so
$\sum_{d=0}^N(d+1)(N+1-d) = (N+1)(N+2)(N+3)/6$ of them ($4, 10, 20$); the dual product shares $AC$ across
the diagonal and has no other nonzero block. (e) LG.10. $\square$

The multiplication counts are `NUMERICS.md` §2.2, [`0005`](../decisions/0005-the-jacobian-is-a-dual-matrix.md)
and `PHASE3.md` §5 (81 against 216, 135 against 729).

**Checked:** mpmath 1.4.1, 110 digits, $N = 1, 2, 3$, 40 unstructured random block sets each: (a) dual
product against the dense product, exact ($0$), associativity exact; (b) `mp.inverse` of the dense matrix
against the dual inverse $3.0\times10^{-110}$ (relative), two-sided identity $1.9\times10^{-106}$,
$\det\iota - (\det A)^{N+1}$ $5.3\times10^{-111}$ (relative); (c) both formulas exact; (d) the multiplications
counted by instrumenting a naive $3\times3$ product equal $27 + 54N$ (the generic column by enumerating
the block products of a lower-triangular product, $N = 1, 2, 3$: $4, 10, 20$; the dense one is arithmetic); (e) $\mathrm{ad}^5 = W^5 + \epsilon\sum_{j=0}^4W^j\Xi_iW^{4-j}$, relative $2.1\times10^{-111}$; the dual
inverse of $\mathrm{Ad}_X$ equals $\mathrm{Ad}_{X^{-1}} = (R^\top, [-R^\top x_i]_\times R^\top)$, $1.3\times10^{-110}$. Rounding of the
code, mpmath 1.4.1 at 60 digits against the dense product and `mp.inverse`, `f64` and `f32`, $N = 1, 2, 3$, 200
seeded cases each, $A$ alternately Gaussian and $U\,\mathrm{diag}(1, s_2, s_3)V^\top$ with $s_2, s_3$ log-uniform on
$[10^{-4}, 1]$ ($\kappa_F$ up to $6\times10^9$), $B_i$ Gaussian of standard deviation 2: `SEn3Jac::mul` within
$1.6\,u$ and `sen3jac_mul` within $2.7\,u$ of the true product (Frobenius, over $\max(\lVert P\rVert_F, 1)$),
`SEn3Jac::inverse` within $0.14\,\kappa u$ and `sen3jac_inverse` within $0.32\,\kappa u$ of the true inverse
(over $\lVert M^{-1}\rVert_F$, $\kappa = \lVert M\rVert_F\lVert M^{-1}\rVert_F$); worst cases of one sample, not
bounds. Script not committed. **Permanent:** proptests `sen3jac_mul_matches_reference`,
`sen3jac_inverse_matches_reference` (`PHASE3.md` §0.0: the inverse on entries uniform on $[-1, 1)$ with
$\kappa u \le 10^{-3}$ only).

## 6. The inverses of the Jacobians

**Corollary SE.12.** For $\theta \notin 2\pi\mathbb Z_{>0}$ (`NUMERICS.md` §12: $\theta < 2\pi$),
$$
J_r^{-1}(\tau) = J_r^{-1}(\varphi) - \epsilon\,J_r^{-1}(\varphi)\,Q(-\rho_i,-\varphi)\,J_r^{-1}(\varphi),\qquad
J_l^{-1}(\tau) = J_l^{-1}(\varphi) - \epsilon\,J_l^{-1}(\varphi)\,Q(\rho_i,\varphi)\,J_l^{-1}(\varphi),
$$
with $J_{r}^{-1}(\varphi) = I + \tfrac12W + cW^2$ and $J_l^{-1}(\varphi) = I - \tfrac12W + cW^2$ (SO.8).
There is no separate closed form for the block: it is built from $Q$ and $c$.

*Proof.* SE.11(b) with $A = J_r(\varphi)$, invertible off $2\pi\mathbb Z_{>0}$ (LG.2(b)), and $B_i = Q(-\rho_i,-\varphi)$
(SE.9(a)); $J_l$ likewise. $\square$

*Consequences.* (i) The inverse block scales like $\lVert\rho_i\rVert$ and, near $2\pi$, like
$(2\pi/(2\pi-\theta))^2$ (LG.16(d)); it costs $2N$ products of $3\times3$ blocks and no dense factorization; its rounding is SE.16(b).
(ii) **The identity $J_r^{-1}J_r = I$ cannot detect a wrong $Q$**: SE.11(b) holds for *every* $B_i$, so a test
of `jr_inv` against `jr` checks $A^{-1}$ only. $Q$ is pinned down by the ids `sen3_jr_n*`,
`sen3_jl_n*` (defined by the series, not by $Q$) and by `jacobians_match_dual_*`; `jr_inv` by
`sen3_jr_inv_n*` (`mp.inverse` of the series). (iii) $J_r^{-1} - J_l^{-1} = \mathrm{ad}_\tau$ (LG §6) holds
blockwise and is a free consistency check of the pair.

**Checked:** mpmath 1.4.1, 110 digits, $N = 1, 2, 3$, 12 random $\tau$ ($\theta<6$) plus the 7 edge $\tau$ of
SE.3 and its two $\rho_i\sim10^6$ cases, excluding $\lvert\theta - 2\pi\rvert<0.02$ ($c$ at 250 digits): the two
formulas against `mp.inverse` of the dense series of $\mp\mathrm{ad}_\tau$, relative to $\max(1,\lVert\cdot\rVert_{\max})$,
$8.0\times10^{-104}$ ($J_r^{-1}$), $1.2\times10^{-103}$ ($J_l^{-1}$) (at $\theta = 20$; $7.0\times10^{-111}$ and
$6.4\times10^{-111}$ for $\rho\sim10^6$); $J_r^{-1} - J_l^{-1} - \mathrm{ad}_\tau$ $6.8\times10^{-109}$. Script not
committed. **Permanent:** planned, corpus `sen3_jr_inv_n{1,2,3}`,
`sen3_jl_inv_n{1,2,3}`; proptest `sen3jac_inverse_matches_reference` and the `jr_inv` twin (`NUMERICS.md` §14).

## 7. The action Jacobians

**Proposition SE.13 (`NUMERICS.md` §2.4, §5.5).** For $X = (R, t) \in \mathrm{SE}(3)$, $F(X, p) = Rp + t$ and
the derivatives $\mathrm D^{\mathrm{right}}$, $\mathrm D^{\mathrm{left}}$ of LG.11:

- (a) $\mathrm D^{\mathrm{right}}_XF = [\,-R[p]_\times,\ R\,] = R\,[\,-[p]_\times,\ I\,]$;
- (b) $\mathrm D^{\mathrm{left}}_XF = [\,-[Rp + t]_\times,\ I\,]$;
- (c) $\partial F/\partial p = R$;
- (d) $\mathrm D^{\mathrm{left}}_XF = \mathrm D^{\mathrm{right}}_XF\,\mathrm{Ad}_X^{-1}$, with
  $\mathrm{Ad}_X^{-1} = \begin{bmatrix}R^\top & 0\\ -R^\top[t]_\times & R^\top\end{bmatrix}$.

The $\mathrm{SO}(3)$ rows are (a), (b) with $t = 0$ and the first block only ($-R[p]_\times$, $-[Rp]_\times$).

*Proof.* $\mathrm{Exp}(\delta) = I + \delta^\wedge + O(\lVert\delta\rVert^2)$, so $\mathrm{Exp}(\delta)[q;1] =
[q + \delta_\varphi\times q + \delta_\rho;\ 1] + O$. Right: $X\,\mathrm{Exp}(\delta)[p;1]$ has translation part
$R(p + \delta_\varphi\times p + \delta_\rho) + t = Xp - R[p]_\times\delta_\varphi + R\delta_\rho$, using
$\delta_\varphi\times p = -[p]_\times\delta_\varphi$. Left: with $q = Rp + t$, $\mathrm{Exp}(\delta)X[p;1]$ gives
$q + \delta_\varphi\times q + \delta_\rho = q - [q]_\times\delta_\varphi + \delta_\rho$. (c) is linear. (d) By
$X\,\mathrm{Exp}(\delta) = \mathrm{Exp}(\mathrm{Ad}_X\delta)X$ (LG.4(b)), $\mathrm D^{\mathrm{right}} = \mathrm D^{\mathrm{left}}\mathrm{Ad}_X$.
Explicitly $[-R[p]_\times,\ R]\,\mathrm{Ad}_{X^{-1}} = [-R[p]_\times R^\top - RR^\top[t]_\times,\ I] = [-[Rp]_\times -
[t]_\times,\ I]$, using $R[p]_\times R^\top = [Rp]_\times$ and $[-R^\top t]_\times R^\top = -R^\top[t]_\times$ (SO.9;
the block of SE.4 at $X^{-1}$). $\square$

No $J$ appears: the perturbation is taken at $\delta = 0$, where $J_r(0) = I$, so the rows hold for every $X$
and $p$ (no domain).

**Checked:** mpmath 1.4.1, 110 digits, 30 random $(R, t, p)$ (rotation vector $\sim N(0,0.9^2I)$,
$t, p \sim N(0,2^2)$): central differences $[F(he_j) - F(-he_j)]/2h$, $h = 10^{-30}$, of the operation itself ($4\times4$ `mp.expm` for $X\,\mathrm{Exp}(\delta)$ and
$\mathrm{Exp}(\delta)X$): (a) $8.2\times10^{-61}$, (b) $1.3\times10^{-60}$, (c) $2.5\times10^{-81}$, $\mathrm{SO}(3)$
$8.2\times10^{-61}$ and $8.7\times10^{-61}$ (the differencing); (d) the predicted forms $6.7\times10^{-111}$, the
differenced ones $8.0\times10^{-61}$; $\mathrm{Ad}^{-1}$ block form $4.1\times10^{-111}$. Script not committed.
**Permanent:** planned, proptests `jacobians_match_dual_*` (`act_jacobians`).

## 8. Rotation-first and translation-first

Let $\Pi_N = \begin{bmatrix}0 & I_{3N}\\ I_3 & 0\end{bmatrix}$, so that $\Pi_N[\varphi;\rho_1;\dots;\rho_N] =
[\rho_1;\dots;\rho_N;\varphi]$; $\Pi_1 = \begin{bmatrix}0&I_3\\ I_3&0\end{bmatrix}$ is the $P$ of the index and
maps $[\omega;v]$ to the literature's $[v;\omega]$ (Barfoot $[\rho;\phi]$; Solà, Sophus, manif
$[\upsilon;\omega]$).

**Proposition SE.14.** Write $\tau_{\mathrm{tf}} = \Pi\tau$.

- (a) $\Pi$ is orthogonal ($\Pi^\top = \Pi^{-1}$; for $N = 1$ also symmetric). The hat matrix of $\tau_{\mathrm{tf}}$ (with
  the blocks read in that order) is the hat matrix of $\tau$: group elements, $\mathrm{Exp}$, $\mathrm{Log}$ and $\theta$
  do not depend on the order.
- (b) A linear map $\tau' = M\tau$ acts on translation-first coordinates as $M_{\mathrm{tf}} = \Pi M\Pi^\top$.
- (c) For $A + \epsilon B$, $\Pi\,\iota(A+\epsilon B)\,\Pi^\top$ is block *upper* triangular, $A$ on the diagonal and
  $B_i$ in the last block column, $\begin{bmatrix}A & B\\ 0 & A\end{bmatrix}$ for $N = 1$. So
  $\mathrm{Ad}_{\mathrm{tf}} = \begin{bmatrix}R & [t]_\times R\\ 0 & R\end{bmatrix}$,
  $\mathrm{ad}_{\mathrm{tf}} = \begin{bmatrix}W & \rho^\wedge\\ 0 & W\end{bmatrix}$ and
  $J_{l,\mathrm{tf}} = \begin{bmatrix}J_l(\varphi) & Q(\rho,\varphi)\\ 0 & J_l(\varphi)\end{bmatrix}$: the literature's
  $Q$ *is* §5.3's, moved from lower-left to upper-right.
- (d) A covariance $\Sigma = \mathbb E[\delta\delta^\top]$ becomes $\Pi\Sigma\Pi^\top$, an information matrix
  $\Pi\Sigma^{-1}\Pi^\top$, a gradient $g$ becomes $\Pi g$, and a Jacobian $J = \partial y/\partial x$ between
  manifolds of orders $\Pi_Y$, $\Pi_X$ becomes $\Pi_YJ\Pi_X^\top$. Propagation commutes:
  $(M\Sigma M^\top)_{\mathrm{tf}} = M_{\mathrm{tf}}\Sigma_{\mathrm{tf}}M_{\mathrm{tf}}^\top$.
- (e) Converting a tangent copies its $3+3N$ coordinates, with no arithmetic: it is exact.

*Proof.* (a) $\Pi$ is a permutation matrix; the hat matrix depends on the blocks, not on where they are stored.
(b) $\Pi\tau' = \Pi M\Pi^\top(\Pi\tau)$. (c) Permuting the block rows and columns of the $(N+1)\times(N+1)$ block
matrix moves index $0$ last: diagonal blocks stay diagonal, $(i,0)\mapsto(i,N)$. (d) $\mathbb E[\Pi\delta(\Pi\delta)^\top]$;
$(\Pi\Sigma\Pi^\top)^{-1} = \Pi\Sigma^{-1}\Pi^\top$; $g^\top\tau = (\Pi g)^\top\Pi\tau$ since $\Pi^\top\Pi = I$; the chain
rule. (e) is the definition. $\square$

*What it does not do.* $\Pi$ acts on tangent coordinates only; $\theta$, the coefficients and the group
element are unchanged. It converts an order, never a retraction: locus-tag's decoupled
$(R\,\mathrm{Exp}\,\omega,\ Rv + t)$ is a different map from $\mathrm{Exp}$, not $\Pi$ of it (`0012`). Using a
translation-first $M$ where a rotation-first one is expected is a large error, not a rounding:
$\mathrm{Ad}_{\mathrm{tf}} - \mathrm{Ad}$ has the block $[t]_\times R$ in both off-diagonal positions (zero iff
$t = 0$, LG §6), and a covariance $\mathrm{diag}(\sigma_\varphi^2 I, \sigma_\rho^2 I)$ (rad$^2$, m$^2$) read in the other order
attaches $\sigma_\rho^2$ to the rotation. For $N \ge 2$ a library may order the blocks differently
(`0002` states only the twist): a block order $\pi$ is the matrix $P_\pi\otimes I_3$, with $P_\pi$ the
permutation of the $N+1$ blocks; (a), (b), (d), (e) hold verbatim and in (c) the column of block $0$ moves
with it. API converters exist only for `Twist`, so the permutations of (c), (d) are the caller's, or the
oracle runner's (`PHASE1.md` §7, `PHASE5.md` §6).

**Checked:** mpmath 1.4.1, 110 digits, $N = 1, 2, 3$, 12 random $(X, \tau, \Sigma)$ each (rotation vector
of $R \sim N(0,0.7^2I)$, $x_i, \rho_i \sim N(0,2^2)$, $\theta(\tau)\in[0.2,2.8]$, $\Sigma = MM^\top + I$). The
translation-first objects were built **in that order from the definitions**, not by permuting: the hat and
vee reordered, $\mathrm{Ad}$ and $\mathrm{ad}$ as images of the translation-first basis, $J_{l,r}$ as the dense series of
$\mathrm{ad}_{\mathrm{tf}}$. They equal $\Pi(\cdot)\Pi^\top$ (rotation-first, built in SE.3–SE.8) exactly for
$\mathrm{Ad}$, $\mathrm{ad}$, $J$; $\Pi\Pi^\top = I$; $\tau_{\mathrm{tf}} = \Pi\tau$ and $\mathrm{Exp}_{\mathrm{tf}}(\Pi\tau) = \mathrm{Exp}(\tau)$ exactly;
$N = 1$: $\mathrm{Ad}_{\mathrm{tf}}$ against $[[R, [t]_\times R],[0,R]]$ $6.7\times10^{-111}$ and the
upper-right block of $J_{l,\mathrm{tf}}$ against §5.3 $8.3\times10^{-111}$; the zero pattern of $J_{\mathrm{tf}}$
(only the diagonal and the last block column) exact; $\mathrm{Ad}_{\mathrm{tf}}\Sigma_{\mathrm{tf}}\mathrm{Ad}_{\mathrm{tf}}^\top = \Pi(\mathrm{Ad}\,\Sigma\,\mathrm{Ad}^\top)\Pi^\top$ exact; the information
matrix $4.1\times10^{-112}$; $g^\top\tau$ exact; the wrong-order gap $\max\lvert\mathrm{Ad}_{\mathrm{tf}} - \mathrm{Ad}\rvert$ is $4.8$ to $5.4$
($x_i \sim N(0,2^2)$). Script not committed. **Permanent:** planned, `twist_translation_first_round_trip`
(`0002`, implementation plan 2); none for the Jacobian and covariance permutations.

## 9. What the exact arms cost

**Proposition SE.15.** Let $\delta b, \delta d, \delta e$ be the errors of the exact-arm values of $b, d, e$
(CO.6: relative $\hat C\,u\,\theta^{-p}$, $\hat C = 6, 45, 367$ sampled, $p = 2, 2, 4$) and $\delta a$ that of
$a = 2k^2$ ($\le5u$ relative). To first order, in the spectral norm, for $\theta \le 1$ (where CO.6's model was
sampled) and with $b, d, e$ at their values $\tfrac16,\tfrac1{24},\tfrac1{120}$ at $\theta = 0$:

- (a) the translation column of $\mathrm{Exp}$: $\lVert\delta(J_l\rho)\rVert \le (\lvert\delta a\rvert\theta + \lvert\delta b\rvert\theta^2)\lVert\rho\rVert
  = O(u)\lVert\rho\rVert$ **at every such $\theta$**: the cancellation of $b$ is cancelled by $W^2$;
- (b) the block: $\lVert\delta Q\rVert \le \lVert\rho\rVert\,[\lvert\delta b\rvert(2\theta + \theta^2) + 5\lvert\delta d\rvert\theta^2 + 2\lvert\delta e\rvert\theta^3]
  \approx \lVert\rho\rVert\,u\,(8.2/\theta + 10.4)$, and since $\lVert Q\rVert \to \lVert\rho\rVert/2$ as $\theta \to 0$,
  $\lVert\delta Q\rVert/\lVert Q\rVert \lesssim (16.3/\theta + 20.8)\,u$ (constants rounded up).

*Proof.* Perturb the coefficients of §5.3 and use $\lVert W\rVert = \theta$, $\lVert\Xi\rVert = \lVert\rho\rVert$:
$\lVert W\Xi + \Xi W + W\Xi W\rVert \le (2\theta + \theta^2)\lVert\rho\rVert$, the $d$-group has norm $\le 5\theta^2\lVert\rho\rVert$
and the $e$-group $\le2\theta^3\lVert\rho\rVert$. With $\lvert\delta b\rvert \approx u/\theta^2$,
$\lvert\delta d\rvert \approx 1.9u/\theta^2$, $\lvert\delta e\rvert \approx 3.1u/\theta^4$ the three terms are $(2/\theta + 1)u$,
$9.4u$, $6.1u/\theta$ (times $\lVert\rho\rVert$). (a) likewise with $\lVert W\rVert$, $\lVert W^2\rVert$. $\square$

*Where it applies.* $Q$ loses **one** power of $\theta$ where $e$ loses four, *if* an exact arm is evaluated there: at
$\theta = 10^{-2}$, where $e$'s exact arm keeps about five digits (`NUMERICS.md` §4), $Q$ would keep about thirteen. But
the exact arms serve only $\theta \ge \theta_s$ (generated, 0004; CO.10 predicts $0.1$ to $1.4$ for $b, d, e$, CO.11 a
best shared one of $0.9$ at $m = 8$), and below it the series arm errs by about $u$ times the word norms. The worst
error of $Q$ is at the switch: the bound there is $(16.3/\theta_s + 20.8)u$, $39u$ at $0.9$ ($170u$ at the $m = 4$ shared
switch $\approx0.11$, CO.11), measured $\le 9.4u$ on $[0.5,1)$; the $m = 8$ series arm errs by $\le 1.4u$ up to $0.9$ but by
$99u$ at $1.4$, past the value crossing of $b$ ($1.16$, CO.10). So `coeff_*` digits understate those of the block
(`sen3_jl_n*`), and 0004 item 1, which minimizes the coefficient's error, does not minimize the consumer's (index, open
items). The bound is relative to $\lVert\rho\rVert$; `NUMERICS.md` §11's floor for a Jacobian is $1$, so the corpus measures
$\lVert\delta Q\rVert/\max(\lVert Q\rVert, 1)$, the same only for $\lVert\rho\rVert \gtrsim 2$ and looser by about $\lVert\rho\rVert/2$ in the
`rho:1e-6`, `rho:1e-3`, `rho:1e0` strata; componentwise errors of small entries of $Q$ are not covered.

**Checked:** `f64` (CPython floats, glibc `sin`, `cos`; $z = \mathrm{fl}(\varphi\cdot\varphi)$, $\hat\theta = \sqrt z$), reference by
mpmath 1.4.1 (160 digits) from the same double inputs $(\varphi,\rho)$; $b$ and $e$ by their definitions,
$d = (\hat\theta^2 - 4\sin^2\frac{\hat\theta}2)/2\hat\theta^4$, $Q$ by §5.3. 400 samples log-uniform per decade of
$\theta\in[10^{-4},1)$, random axis, $\rho$ a random unit vector: $\lVert\hat Q - Q\rVert_F/\lVert Q\rVert_F$ at most $6.2, 5.9, 6.3, 6.5$
times $u/\theta$ (per decade, from $10^{-4}$; Frobenius as in `NUMERICS.md` §11, the bound is spectral: a factor $<\sqrt3$
apart), against the bound $16.3$; over $3000$ samples in each of the first three decades $b$, $d$, $e$ reached $6.0\,u\theta^{-2}$,
$41.6\,u\theta^{-2}$, $350\,u\theta^{-4}$. The translation column $J_l\rho$ ($a = 2k^2$, $b$ by its definition): $1.6, 1.7, 1.4,
1.6\,u\lVert\rho\rVert$. Beyond, 400 samples per bin, exact arms at most $49.9, 21.6, 9.4, 4.6, 7.9\,u$ on $[0.1,0.2)$,
$[0.2,0.5)$, $[0.5,1)$, $[1,1.5)$, $[1.5,3)$ (never above $0.42$ of the bound); series arm ($m = 8$ Horner in $z$ of CO.1's
rational coefficients of $b$, $d$, $e = T_5$ rounded to `f64`; an illustration, not a generated arm) $1.2, 1.2, 1.3, 1.4\,u$ on
$[10^{-6},10^{-3})$, $[10^{-3},0.1)$, $[0.1,0.5)$, $[0.5,0.9)$ and $99u$ on $[0.9,1.4)$ (at $1.39$). Sampled maxima, not proved
bounds; scripts not committed. **Permanent:** none for the bound; the strata of `sen3_exp_n*`, `sen3_jl_n*` (`rho:*` crossed
with $\theta$, `PHASE1.md` §4.4) measure the quantities themselves.

## 10. The inverse side

**Proposition SE.16 (rounding of `log`, `jr_inv`, `jl_inv`).** First order, spectral norms, $s = 2\sin(\theta/2)/\theta$
(LG.16), the exact arm of $c$ ($\hat C = 48$, $p = 2$, CO.6).

- (a) *The translation part of $\mathrm{Log}$, $\theta(X) < \pi$.* If $\hat\varphi = \varphi + \delta\varphi$ is the computed rotation part and
  $\hat\rho_i = \mathrm{fl}(\hat J_l^{-1}(\hat\varphi)x_i)$, then $\hat\rho_i - \rho_i = \Lambda_i\,\delta\varphi + \delta c\,W^2x_i + {}$rounding, with
  $\Lambda_i = J_l^{-1}\big([x_i]_\times J_l - Q(\rho_i,\varphi)\big)$, $\lVert \Lambda_i\rVert \le \kappa_x(\theta)\lVert x_i\rVert$,
  $\kappa_x = \tfrac12 + 2c\theta + 2c'\theta^3 \le 1.29$ ($c' = \mathrm dc/\mathrm dz$), and $\lvert\delta c\rvert\theta^2 \le 48\,c\,u \le 4.9u$: the
  $\theta^{-2}$ of the exact arm is cancelled by $W^2$, as in SE.15(a).
- (b) *The blocks of $J^{-1}$ (SE.12).* With $\varepsilon_J = \lVert\hat J^{-1}(\varphi) - J^{-1}(\varphi)\rVert$, $B_i = -J^{-1}QJ^{-1}$ and
  $\hat B_i$ its `f64` evaluation, $\lVert\hat B_i - B_i\rVert \le \dfrac{\lVert Q\rVert}{s^2}\Big(\dfrac{\lVert\delta Q\rVert}{\lVert Q\rVert} + 2s\,\varepsilon_J + c_pu\Big)$,
  $c_p \le 1$ (measured, the two $3\times3$ products) and $\varepsilon_J \le 3.1\,u/s^2$ (measured). $\lVert Q\rVert/s^2$ is the scale of LG.16(d).

*Proof.* (a) Differentiating $x_i = J_l(\varphi)\rho_i(\varphi)$ at fixed $x_i$ with SE.9(c) gives
$(Q - [x_i]_\times J_l)\,\delta\varphi + J_l\,\delta\rho = 0$. For the bound, $\rho = x - \tfrac12\varphi\times x + c(z)\,\varphi\times(\varphi\times x)$
differentiates to $-\tfrac12\delta\times x + 2c'(\varphi\cdot\delta)\,\varphi\times(\varphi\times x) + c\,(\delta\times(\varphi\times x) + \varphi\times(\delta\times x))$, of
norm $\le\kappa_x\lVert\delta\rVert\lVert x\rVert$. The series coefficients of $c$ are positive (SO.8), so $\kappa_x$ increases, to
$\kappa_x(\pi) = 1.285$ ($c(\pi) = \pi^{-2}$, $c'(\pi) = 2.4\times10^{-3}$); and $\lvert\delta c\rvert = \hat C\,c\,u\,\theta^{-2}$. (b) First order in
$\hat B_i - B_i = -J^{-1}\delta Q\,J^{-1} - E\,QJ^{-1} - J^{-1}Q\,E$, $E = \hat J^{-1} - J^{-1}$, $\lVert J^{-1}\rVert = 1/s$ (LG.16(a)), and the
rounding of the products. $\square$

*Consequences.* (a) The error of $\hat\rho_i$ is absolute, in units of $\lVert x_i\rVert$, and does not grow toward $\pi$ ($1/s \le \pi/2$
there). It inherits $\delta\varphi$: $\lesssim2.3u\theta$ from a quaternion (SO.6(c), no amplification; measured), the
ill-conditioned axis from a *matrix* near $\pi$ (`NUMERICS.md` §11: backward error only), times $\kappa_x \le 1.29$. (b) Against
$\lVert Q\rVert/s^2$ the sandwich adds only $2s\varepsilon_J$ to $Q$'s error; but $\lVert B_i\rVert \ge \lVert Q\rVert$ ($\lVert J\rVert = 1$), so
against the block itself the loss is at most $1/s^2$: $\le 2.5$ for $\theta \le \pi$, unbounded toward $2\pi$, where
$\varepsilon_J \approx 1.5\,u/s^2$ (measured): the conditioning of $c \approx 1/(2\pi(2\pi - \theta))$ (SO.8) at the rounded argument.

**Checked:** `f64`, references by mpmath 1.4.1 (160 digits) from the same double inputs, exact arms of $b, d, e, c$ at every $\theta$.
(a) `Log` of `NUMERICS.md` §3.2 ($\hat\varphi = 2\,\mathrm{atan2}(n,w)\,q_{\mathrm v}/n$, $q_{\mathrm v}$ the vector part) then $\hat J_l^{-1}(\hat\varphi)x_i$, from a double quaternion and
$x_i$ built from $\tau$ ($\rho\sim N(0,2^2)$); reference: $\Gamma_1$ (series) solved at the reference $\varphi$; 300 uniform $\theta$ in
each of $[0.05,0.5)$, $[0.5,1)$, $[1,2)$, $[2,3)$, $[3,\pi-10^{-3})$, $[\pi-10^{-3},\pi-10^{-6})$: $\lVert\delta\varphi\rVert \le 2.3\,u\theta$;
$\lVert\hat\rho - \rho\rVert/\lVert x\rVert \le 3.2, 3.1, 3.3, 4.1, 5.9, 5.1\,u$, never above $\kappa_x\lVert\delta\varphi\rVert + 4.9u\lVert x\rVert$ (largest ratio
$0.89$); $\lVert\Lambda\rVert/\lVert x\rVert \le 0.93$, under $\kappa_x$ without exception; $\Lambda$ against central differences ($h = 10^{-25}$, 60
digits) of $\varphi \mapsto \Gamma_1^{-1}x$, $1.5\times10^{-35}$ (the differencing). (b) $N = 1$, $\rho\sim N(0,2^2)$, 300 uniform
$\theta$ in each of $[0.05,0.5)$, $[0.5,1)$, $[1,2)$, $[2,\pi)$, $[\pi,5)$, $[5,6)$, $[6,2\pi-0.05)$, $[2\pi-0.05,2\pi-0.005)$; per bin, in $u$:
$\lVert\delta Q\rVert/\lVert Q\rVert \le 68.8, 9.6, 5.3, 6.4, 17.8, 22.5, 26.7, 24.4$; $2s\varepsilon_J \le 6.0, 6.4, 5.0, 3.8, 9.6, 44, 222, 1775$; the left
side over $\lVert Q\rVert/s^2$ $\le 68.6, 10.6, 7.0, 7.2, 13.6, 25.8, 124, 903$ ($c_p \le 0.91$: the largest per-sample excess over the
first two) and over $\lVert B_i\rVert$ $\le 68.6, 10.6, 7.0, 8.0, 15.9, 47.7, 325, 1765$; $\varepsilon_Js^2 \le 3.08u$. Sampled maxima, not
proved bounds; scripts not committed. **Permanent:** none specified; the strata of `sen3_jr_inv_n*`, `sen3_jl_inv_n*` measure the
quantities themselves.
