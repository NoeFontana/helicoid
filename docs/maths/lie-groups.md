# Lie groups, adjoints and the Jacobians of Exp

> Non-normative companion to [`NUMERICS.md`](../NUMERICS.md) §1, §2.2, §2.3, §5.1–§5.2 and the
> algebra of §5.4. Notation, the `Checked:` convention and the open items are in the
> [maths index](./index.md). The SE(2) case, where the same algebra gives the Jacobians, is
> [`so2-se2.md`](./so2-se2.md).

Throughout, $N \ge 0$, $n = 3 + 3N$, $G = \mathrm{SE}_N(3)$, rotation-first. A domain "$\theta < 2\pi$"
is on the tangent's $\theta = \lVert\varphi\rVert$; "$\theta(X) < \pi$" is on the rotation angle of $X$
($\operatorname{tr}R = 1 + 2\cos\theta$).

## Results

| Label | Result | `NUMERICS.md` | Implemented by |
|---|---|---|---|
| LG.2 | $\mathrm{Exp}(\tau) = (\mathrm{Exp}\,\varphi,\ \Gamma_1(\varphi)\rho_i)$; $\mathrm{Log}$ inverts it on $\theta < \pi$ and jumps at $\theta(X) = \pi$ | §1, §5.1 | `LieGroup::exp`, `log` |
| LG.4 | $X\,\mathrm{Exp}(\tau)X^{-1} = \mathrm{Exp}(\mathrm{Ad}_X\tau)$; blocks of $\mathrm{Ad}_X$, $\mathrm{ad}_\tau$ | §1, §5.2 | `adjoint`, `ad` |
| LG.5 | $\mathrm{Ad}_{\mathrm{Exp}\,\tau} = \exp(\mathrm{ad}_\tau)$ | §1, §14 | (used by LG.9) |
| LG.8 | $\mathrm{Exp}(\tau + \delta) = \mathrm{Exp}(\tau)\mathrm{Exp}(J_r\delta) + O = \mathrm{Exp}(J_l\delta)\mathrm{Exp}(\tau) + O$, with $J_r = \sum(-\mathrm{ad}_\tau)^k/(k+1)!$, $J_l = \sum \mathrm{ad}_\tau^k/(k+1)!$ | §1 | `jr`, `jl` |
| LG.9 | $J_l(\tau) = J_r(-\tau) = \mathrm{Ad}_{\mathrm{Exp}\,\tau}J_r(\tau)$; $\det J = (2(1-\cos\theta)/\theta^2)^{N+1}$; $J$ invertible iff $\theta \notin 2\pi\mathbb Z_{>0}$ | §1, §12, §14 | `jl`, `jl_inv`, `jr_inv` |
| LG.10 | dual-block matrices are closed under product and inverse and contain $\mathrm{Ad}$, $\mathrm{ad}$, $J$, $J^{-1}$ | §2.2, §5.2 | `SEn3Jac` |
| LG.13 | first-order identities for $\mathrm{Exp}$, $\mathrm{Log}$ and conjugation | §1 | — |
| LG.14 | every row of the §2.3 table, both sides | §2.3 | `rplus_jacobians`, `lplus_jacobians`, `rminus_jacobians`, `lminus_jacobians`, `compose_jacobians`, `inverse_jacobian` |
| LG.15 | the $\ominus$ rows again, by the chain rule | §2.3, §14 | reference twin of `*_jacobians` |
| LG.16 | conditioning of $J$ and $\mathrm{Ad}$: $\kappa_2 = \theta/(2\sin\tfrac\theta2)$ for SO(3); $\sim(\lVert\rho\rVert/2)^2$ at $\varphi = 0$ | §12 | — |

## 1. Setting

**Definition LG.1 (groups, hat, Exp, Log).** $G$ is the group of matrices
$X = \begin{bmatrix}R & x_1\cdots x_N\\ 0 & I_N\end{bmatrix}$, $R \in \mathrm{SO}(3)$. It is closed
in $\mathrm{GL}_{3+N}(\mathbb R)$, hence an embedded Lie group with Lie algebra
$\mathfrak g = \{A : \exp(tA) \in G\ \forall t\}$ [Hall, ch. 3]. For $\tau = [\varphi;\rho_1;\dots;\rho_N]$,

$$
\tau^\wedge = \begin{bmatrix}\varphi^\wedge & \rho_1\cdots\rho_N\\ 0 & 0\end{bmatrix}.
$$

$\tau \mapsto \tau^\wedge$ is a linear bijection $\mathbb R^n \to \mathfrak g$ with inverse $\vee$
($\supseteq$ by LG.2(a); $\subseteq$ by differentiating $R^\top R = I$ and the constant last block
rows of a curve in $G$ at $t = 0$). $\mathrm{Exp}(\tau) = \exp(\tau^\wedge)$, an everywhere convergent
series. $\mathrm{Log}$ is the inverse of $\mathrm{Exp}$ on $\{\theta < \pi\}$ (LG.2(c)); at
$\theta(X) = \pi$ it is one of two preimages (LG.2(d)), picked by the quaternion sign of
`NUMERICS.md` §3.2.

**Proposition LG.2 (Exp and Log).**

- (a) $\mathrm{Exp}(\tau) = \begin{bmatrix}e^W & \Gamma_1(\varphi)\rho_1\ \cdots\ \Gamma_1(\varphi)\rho_N\\ 0 & I_N\end{bmatrix}$,
  $\Gamma_1(\varphi) = \sum_{k\ge0} W^k/(k+1)!$.
- (b) $\det\Gamma_1(\varphi) = 2(1 - \cos\theta)/\theta^2 = \big(\sin(\theta/2)/(\theta/2)\big)^2$, which is $1$
  at $\theta = 0$ (removable), so $\Gamma_1$ is invertible iff $\theta \notin 2\pi\mathbb Z_{>0}$.
- (c) $\mathrm{Exp}$ maps $\{\theta < \pi\}$ bijectively and analytically onto $\{X : \theta(X) < \pi\}$,
  with inverse $\mathrm{Log}(X) = (\varphi = \mathrm{Log}\,R,\ \rho_i = \Gamma_1(\varphi)^{-1}x_i)$. The
  rotation part is not injective beyond: $\mathrm{Exp}(\varphi) = \mathrm{Exp}((1 - 2\pi/\theta)\varphi)$,
  which at $\theta = \pi$ is $\mathrm{Exp}(-\varphi)$.
- (d) In $\{\theta < 2\pi\}$, an $X$ with $\theta(X) = \pi$ has exactly two preimages,
  $(\varphi, \rho_i)$ and $(-\varphi, \rho_i')$ with $\rho_i' - \rho_i = \varphi \times x_i$: the translation
  blocks differ, not just the sign of $\varphi$. $\mathrm{Log}$ is discontinuous at every such $X$
  (its rotation block jumps by $2\pi$), so $\mathrm D\,\mathrm{Log}$ does not exist there.

*Proof.* (a) For $A = \tau^\wedge$ induction gives
$A^k = \begin{bmatrix}W^k & W^{k-1}\rho\\ 0 & 0\end{bmatrix}$, $k \ge 1$, $\rho = [\rho_1\cdots\rho_N]$;
sum $\sum_k A^k/k!$. (b) $W$ has spectrum $\{0, \pm i\theta\}$, and $\Gamma_1 = g(W)$ with the entire
function $g(z) = (e^z - 1)/z$, $g(0) = 1$. Its eigenvalues are $1$ and
$g(\pm i\theta) = e^{\pm i\theta/2}\,2\sin(\theta/2)/\theta$, so
$\det = \lvert g(i\theta)\rvert^2 = 2(1-\cos\theta)/\theta^2$ ($=1$ at $\theta = 0$, where $\Gamma_1 = I$), zero iff
$\theta \in 2\pi\mathbb Z_{>0}$.
(c) For $R$ this is the axis–angle parametrization of $\mathrm{SO}(3)$ [Hall, ch. 1]: $e^W$ is the
rotation by $\theta$ about $\varphi/\theta$, equal to the rotation by $2\pi - \theta$ about $-\varphi/\theta$,
and $\theta \mapsto \operatorname{tr}R$ is injective on $[0,\pi]$. The translation part is (a) with (b).
Analyticity: the inverse function theorem, $d\,\mathrm{Exp}$ being invertible (LG.9). (d) $\Gamma_1(\varphi)^{-1} = h(W)$ and
$\Gamma_1(-\varphi)^{-1} = h(-W)$ with $h(z) = z/(e^z - 1)$, holomorphic on $\lvert z\rvert < 2\pi$;
$h(-z) = ze^z/(e^z - 1)$, so $h(-z) - h(z) = z$ and $\rho_i' - \rho_i = Wx_i$. Both preimages map to $X$ by (a).
For the jump, $R^\pm_\epsilon = \mathrm{Exp}(\pm(\pi - \epsilon)\hat n)$ both tend to
$\mathrm{Exp}(\pi\hat n)$ while $\mathrm{Log}\,R^\pm_\epsilon = \pm(\pi - \epsilon)\hat n \to \pm\pi\hat n$;
the rotation block of $\mathrm{Log}\,X$ is $\mathrm{Log}\,R$. $\square$

**Checked:** mpmath 1.4.1. Sampled, 110 digits, each $N = 0..3$: 40 random $\tau$ (axis uniform on
$S^2$, $\theta \sim U[0, 6]$, $\rho_i \sim \mathcal N(0, 2^2)$) plus edge $\tau$: $\theta \in \{0, 10^{-10}, \pi \pm 10^{-3},
2\pi - 10^{-3}, 6.5, 20\}$ and two with $\rho_i \sim 10^6$. (a) the series of $\Gamma_1$ against `mp.expm`, max
$8.2\times10^{-105}$; (b) its determinant against the closed form, max $1.2\times10^{-92}$, all of it the
cancellation in evaluating $2(1 - \cos\theta)/\theta^2$ at $\theta = 10^{-10}$ (against a 200-digit closed form:
$1.4\times10^{-111}$; on the random samples: $3.2\times10^{-110}$). 100 digits, $N = 1..3$, 40 random $\tau$ with
$\theta < 3.1$: (c) $\Gamma_1^{-1}x_i$ against the $\rho_i$ used to build $X$, $2.2\times10^{-100}$;
$\mathrm{Log}(\mathrm{Exp}\,\tau) = \tau$ (geometric $\mathrm{SO}(3)$ inverse, $\Gamma_1$ solve), $1.2\times10^{-100}$;
the rotation identity for $\theta \in (\pi, 2\pi)$, $7.1\times10^{-101}$. (d) ten $\tau$ per $N$ at $\theta = \pi$:
$\mathrm{Exp}(-\varphi, \rho') = X$, $1.8\times10^{-100}$; $\rho'_i - \rho_i = \varphi \times x_i$, $4.0\times10^{-100}$.
Script not committed. **Permanent:** planned, corpus `sen3_exp_n{1,2,3}`, `sen3_log_n{1,2,3}`; proptest
`exp_log_roundtrip_*`.

## 2. Adjoints

**Definition LG.3.** For $X \in G$, $A \in \mathfrak g$: $\mathrm{Ad}_X A = XAX^{-1}$, which lies in
$\mathfrak g$ because $\exp(tXAX^{-1}) = X\exp(tA)X^{-1} \in G$. On coordinates,
$\mathrm{Ad}_X\tau = (X\tau^\wedge X^{-1})^\vee$ is an $n \times n$ matrix, and
$\mathrm{ad}_\tau\sigma = [\tau^\wedge, \sigma^\wedge]^\vee$.

**Proposition LG.4.** For all $X, Y \in G$ and $\tau, \delta \in \mathbb R^n$ (no domain):

- (a) $X\,\mathrm{Exp}(\tau)X^{-1} = \mathrm{Exp}(\mathrm{Ad}_X\tau)$.
- (b) $\mathrm{Ad}_{XY} = \mathrm{Ad}_X\mathrm{Ad}_Y$ and $\mathrm{Ad}_{X^{-1}} = \mathrm{Ad}_X^{-1}$; hence
  $X\,\mathrm{Exp}(\delta) = \mathrm{Exp}(\mathrm{Ad}_X\delta)\,X$ and
  $\mathrm{Exp}(\delta)\,X = X\,\mathrm{Exp}(\mathrm{Ad}_X^{-1}\delta)$.
- (c) $\mathrm{Ad}_X = R + \epsilon\,[x_i]_\times R$ and $\mathrm{ad}_\tau = W + \epsilon\,[\rho_i]_\times$
  (`NUMERICS.md` §5.2): $R$ (resp. $W$) on each of the $N+1$ diagonal blocks, $[x_i]_\times R$ (resp.
  $[\rho_i]_\times$) in block $(i, 0)$, zero elsewhere.

*Proof.* (a) $A \mapsto XAX^{-1}$ is linear and multiplicative, so $X A^k X^{-1} = (XAX^{-1})^k$;
sum the series. (b) $(XY)A(XY)^{-1} = X(YAY^{-1})X^{-1}$; the two consequences are (a) read
with $\tau \to \delta$ and with $X \to X^{-1}$. (c) With $X^{-1} = \begin{bmatrix}R^\top & -R^\top x\\ 0 & I\end{bmatrix}$,
$$
X\tau^\wedge X^{-1} = \begin{bmatrix}RWR^\top & R\rho - RWR^\top x\\ 0 & 0\end{bmatrix},\qquad
RWR^\top = (R\varphi)^\wedge,
$$
and column $i$ of the translation block is $R\rho_i - (R\varphi)\times x_i = R\rho_i + [x_i]_\times R\varphi$.
So $\mathrm{Ad}_X\tau = [R\varphi;\ R\rho_i + [x_i]_\times R\varphi]$. For $\sigma = [\psi;\eta_i]$ the
commutator has rotation block $[W, \psi^\wedge] = (\varphi\times\psi)^\wedge$ and translation columns
$W\eta_i - \psi^\wedge\rho_i = \varphi\times\eta_i + \rho_i\times\psi$, so
$\mathrm{ad}_\tau\sigma = [\varphi\times\psi;\ \varphi\times\eta_i + [\rho_i]_\times\psi]$. $\square$

**Checked:** mpmath 1.4.1, 110 digits, each $N = 0..3$: 40 random $(X, \tau)$ (rotation angle of $X$
uniform on $[0, 0.99\pi]$, $x_i \sim \mathcal N(0, 2^2)$, $\tau$ as in LG.2) and three edge $X$
($x_i \sim 10^6$; $\theta(X) = 10^{-10}$; $\theta(X) = \pi - 10^{-3}$); $\mathrm{Ad}_X$ and $\mathrm{ad}_\tau$ built
from Definition LG.3 (images of the basis). (a) max $5.8\times10^{-111}$ (relative to $\max(1, \lVert X\rVert^2)$);
(b) $\mathrm{Ad}_{XY}$ and $\mathrm{Ad}_{X^{-1}}$ max $3.1\times10^{-111}$; (c) block forms max
$1.9\times10^{-111}$ ($\mathrm{Ad}$), $0$ ($\mathrm{ad}$). Script not committed. **Permanent:** proptest `adjoint_identity_*` for (a); corpus `sen3_ad_n{1,2,3}` for (c)
(both planned).

**Proposition LG.5.** $\mathrm{Ad}_{\mathrm{Exp}\,\tau} = \exp(\mathrm{ad}_\tau)$.

*Proof.* On $(3+N)$-square matrices let $\mathrm{lm}_A M = AM$ and $\mathrm{rm}_A M = MA$. They commute, and
$e^{\mathrm{lm}_A} = \mathrm{lm}_{e^A}$, so $M \mapsto e^A M e^{-A}$ equals
$e^{\mathrm{lm}_A}e^{-\mathrm{rm}_A} = e^{\mathrm{lm}_A - \mathrm{rm}_A}$, where $\mathrm{lm}_A - \mathrm{rm}_A = [A, \cdot]$. Restrict to the invariant subspace $\mathfrak g$ and pass to coordinates:
$[\tau^\wedge, \cdot] \leftrightarrow \mathrm{ad}_\tau$, $e^{\tau^\wedge}(\cdot)e^{-\tau^\wedge} \leftrightarrow \mathrm{Ad}_{\mathrm{Exp}\,\tau}$. $\square$

**Checked:** mpmath 1.4.1, 110 digits, the $\tau$ of LG.2 (random and edge, $\theta$ up to $20$, $\rho_i$
up to $10^6$), each $N = 0..3$: `mp.expm` of the $\mathrm{ad}_\tau$ matrix against
$\mathrm{Ad}_{\mathrm{Exp}\,\tau}$ built from LG.3, max $1.9\times10^{-111}$ (relative to $\max(1, \lVert\mathrm{Ad}\rVert)$).
Script not committed. **Permanent:** none specified (LG.9 uses it; `jl_is_ad_jr_*` exercises the
consequence).

## 3. Jacobians of Exp

**Lemma LG.6 (Duhamel [Wilcox; AMH]).** For square matrices $A, B$ and a submultiplicative norm,
$e^{A+B} = e^A + \int_0^1 e^{sA}\,B\,e^{(1-s)A}\,ds + \mathcal R$, with
$\lVert\mathcal R\rVert \le \tfrac12\lVert B\rVert^2 e^{\lVert A\rVert + \lVert B\rVert}$.

*Proof.* Let $f(s) = e^{s(A+B)}e^{(1-s)A}$. Then $f'(s) = e^{s(A+B)}\,B\,e^{(1-s)A}$, so exactly
$$
e^{A+B} - e^A = f(1) - f(0) = \int_0^1 e^{s(A+B)}\,B\,e^{(1-s)A}\,ds. \tag{*}
$$
Apply $(*)$ to $(sA, sB)$: $\lVert e^{s(A+B)} - e^{sA}\rVert \le s\lVert B\rVert\int_0^1 e^{rs\lVert A+B\rVert + (1-r)s\lVert A\rVert}dr
\le s\lVert B\rVert e^{s(\lVert A\rVert + \lVert B\rVert)}$. Then
$\lVert\mathcal R\rVert = \lVert\int_0^1 (e^{s(A+B)} - e^{sA})Be^{(1-s)A}ds\rVert
\le \lVert B\rVert^2e^{\lVert A\rVert}\int_0^1 se^{s\lVert B\rVert}ds \le \tfrac12\lVert B\rVert^2e^{\lVert A\rVert+\lVert B\rVert}$. $\square$

**Definition LG.7.** With the entire function $g(z) = (e^z - 1)/z = \sum_{k\ge0}z^k/(k+1)!$:
$J_l(\tau) = g(\mathrm{ad}_\tau) = \sum_k \mathrm{ad}_\tau^k/(k+1)!$ and
$J_r(\tau) = g(-\mathrm{ad}_\tau) = \sum_k (-\mathrm{ad}_\tau)^k/(k+1)!$.

**Proposition LG.8 (Jacobians of Exp).** For every $\tau$ (any $\theta$), as $\delta \to 0$,
$$
\mathrm{Exp}(\tau + \delta) = \mathrm{Exp}(\tau)\,\mathrm{Exp}(J_r(\tau)\delta) + O(\lVert\delta\rVert^2)
= \mathrm{Exp}(J_l(\tau)\delta)\,\mathrm{Exp}(\tau) + O(\lVert\delta\rVert^2).
$$

*Proof.* Take $A = \tau^\wedge$, $B = \delta^\wedge$ in LG.6 and multiply by $e^{-A}$ on the left. With
$\lambda = 1 - s$, and LG.4(a), LG.5,
$$
e^{-A}\!\int_0^1 e^{sA}Be^{(1-s)A}ds = \int_0^1 e^{-\lambda A}Be^{\lambda A}d\lambda
= \Big(\int_0^1 \mathrm{Ad}_{\mathrm{Exp}(-\lambda\tau)}\delta\,d\lambda\Big)^{\wedge}
= \Big(\int_0^1 e^{-\lambda\,\mathrm{ad}_\tau}d\lambda\ \delta\Big)^{\wedge} = (J_r(\tau)\delta)^\wedge,
$$
using $\int_0^1 e^{-\lambda Z}d\lambda = \sum_k(-Z)^k/(k+1)!$. So $\mathrm{Exp}(\tau+\delta) = \mathrm{Exp}(\tau)(I + (J_r\delta)^\wedge) + O(\lVert\delta\rVert^2)$,
and $\mathrm{Exp}(J_r\delta) = I + (J_r\delta)^\wedge + O(\lVert\delta\rVert^2)$ (the $O$ of LG.6 is explicit). Multiplying by $e^{-A}$ on the
right instead gives $\int_0^1 e^{sA}Be^{-sA}ds = (J_l\delta)^\wedge$ and the second form. $\square$

**Checked:** mpmath 1.4.1, 110 digits, each $N = 0..3$: the numerical derivative of $\mathrm{Exp}$ itself
(rows "$\mathrm{Exp}(\tau)$" of LG.14, both sides, $\theta(\tau)$ up to $20$, beyond $\pi$ and $2\pi$)
against $J_r$, $J_l$ summed from Definition LG.7 until the term norm is below $10^{-120}$; max
$4.6\times10^{-62}$ (the differencing, see LG.14). Script not committed. **Permanent:** planned, corpus
`so3_jr`, `so3_jl`, `sen3_jr_n{1,2,3}`, `sen3_jl_n{1,2,3}` (defined by exactly these series in
`PHASE1.md` §4.3); proptest `jacobians_match_dual_*`.

**Proposition LG.9 (relations and invertibility).**

- (a) $J_l(\tau) = J_r(-\tau)$.
- (b) $J_l(\tau) = \mathrm{Ad}_{\mathrm{Exp}\,\tau}\,J_r(\tau)$.
- (c) $\det J_r(\tau) = \det J_l(\tau) = \big(2(1-\cos\theta)/\theta^2\big)^{N+1}$ (value $1$ at $\theta = 0$);
  $J_r$ and $J_l$ are invertible iff $\theta \notin 2\pi\mathbb Z_{>0}$. So $J^{-1}$ exists on the
  domain $\theta < 2\pi$ of `NUMERICS.md` §12, and also on the components beyond $2\pi$, which §12
  does not promise.

*Proof.* (a) $\mathrm{ad}_{-\tau} = -\mathrm{ad}_\tau$. (b) $g(z) = e^z g(-z)$, since
$e^z(1 - e^{-z})/z = (e^z - 1)/z$; put $z = \mathrm{ad}_\tau$ and use LG.5. (c) By LG.4(c), $\mathrm{ad}_\tau$
is block lower triangular with $N+1$ diagonal blocks $W$, so its spectrum is $\{0, \pm i\theta\}$, each
with multiplicity $N+1$. By spectral mapping the eigenvalues of $J_l$ are $g(0) = 1$ and
$g(\pm i\theta)$ ($J_r$: $g(0)$, $g(\mp i\theta)$), so the determinant is
$(\lvert g(i\theta)\rvert^2)^{N+1}$, as in LG.2(b) ($\lvert g(i\theta)\rvert^2 = (\sin(\theta/2)/(\theta/2))^2$,
$1$ at $\theta = 0$). $\square$

**Checked:** mpmath 1.4.1, 110 digits, the $\tau$ of LG.2 (random and edge), each $N = 0..3$: (a) max $0$;
(b) max $1.1\times10^{-104}$ (relative to $\max(1, \lVert J\rVert)$); (c) relative max $4.6\times10^{-92}$
(the closed form's cancellation at $\theta = 10^{-10}$, as in LG.2(b); $2.3\times10^{-109}$ on the random samples). Script not committed. **Permanent:** proptest
`jl_is_ad_jr_*` for (b) (planned); (c) none specified.

**Proposition LG.10 (dual-block structure).** Let $\mathcal D_N \subset \mathbb R^{n\times n}$ be the
matrices with $N+1$ equal diagonal $3\times3$ blocks $A$, blocks $B_1,\dots,B_N$ in the first block
column below the diagonal, and zero elsewhere ($A + \epsilon B$, `NUMERICS.md` §2.2). Then
$(A+\epsilon B)(C+\epsilon D) = AC + \epsilon(B_iC + AD_i)$, and if $A$ is invertible
$(A+\epsilon B)^{-1} = A^{-1} - \epsilon\,A^{-1}B_iA^{-1}$. Consequently $\mathrm{Ad}_X$, $\mathrm{Ad}_X^{-1}$,
$\mathrm{ad}_\tau$, $J_r$, $J_l$ and $J_r^{-1}$, $J_l^{-1}$ lie in $\mathcal D_N$, and the diagonal block of
$J_r(\tau)$ is the $\mathrm{SO}(3)$ Jacobian $J_r(\varphi)$ (likewise $J_l$).

*Proof.* Block $(i, j)$ of a product is $\sum_k M_{ik}N'_{kj}$. Column $j \ge 1$ of an element of
$\mathcal D_N$ has only its diagonal block, so blocks $(i, j)$, $j \ge 1$, $i \ne j$, vanish, and block $(0, j)$
vanishes; the diagonal blocks are $AC$; block $(i, 0)$ is $B_iC + AD_i$. The inverse formula follows by
multiplying out. $\mathcal D_N$ is thus a matrix algebra containing $I$; it is closed under limits
(finite-dimensional), so it contains every power series in $\mathrm{ad}_\tau$, hence $J_r$, $J_l$; it
contains $\mathrm{Ad}_X$, $\mathrm{ad}_\tau$ by LG.4(c), and the inverses by the formula. The diagonal-block
map $A + \epsilon B \mapsto A$ is an algebra homomorphism, and $\mathrm{ad}_\tau \mapsto W$. $\square$

**Checked:** mpmath 1.4.1, 110 digits, each $N = 0..3$: product and inverse formulas against dense
$3(N+1)$-square products and `mp.inverse` on six unstructured random block pairs, max $0$ and
$1.3\times10^{-111}$; zero pattern and equal diagonal blocks of dense $\mathrm{Ad}_X$, $\mathrm{ad}_\tau$,
$J_r$, $J_l$ and (for $\theta < 2\pi - 10^{-3}$) $J_r^{-1}$, $J_l^{-1}$ over the $\tau$ of LG.2, max
deviation $3.3\times10^{-108}$. Script not committed. **Permanent:** proptests `sen3jac_mul_matches_reference`,
`sen3jac_inverse_matches_reference` ([`0005`](../decisions/0005-the-jacobian-is-a-dual-matrix.md)).

## 4. First-order calculus with sides

**Definition LG.11 (side-relative derivative).** For $s \in \{R, L\}$ and a space $M$ that is either
$G$ or $\mathbb R^n$, define $x \oplus_s \delta$ as in the index on $G$, and $x + \delta$ on
$\mathbb R^n$ for both sides; $y \ominus_s x$ likewise, $y - x$ on $\mathbb R^n$. For a differentiable
$F$, $\mathrm D^s F(x)$ is the matrix with
$$
F(x \oplus_s \delta) = F(x) \oplus_s \big(\mathrm D^s F(x)\,\delta\big) + O(\lVert\delta\rVert^2),
\qquad\text{equivalently}\qquad
\mathrm D^s F(x) = \partial_\delta\big[F(x\oplus_s\delta)\ominus_s F(x)\big]_{\delta=0}.
$$
(The two agree because $\mathrm{Log}(\mathrm{Exp}(\eta)) = \eta$ for small $\eta$.) For several arguments
the partial derivative freezes the others. The input and the output carry the **same side**: this is
what "expressed in that side's perturbation convention" means in `NUMERICS.md` §2.3, and what
`Side::plus_jacobians` and `minus_jacobians` return.

**Lemma LG.12 (chain rule).** $\mathrm D^s(H\circ F)(x) = \mathrm D^s H(F(x))\,\mathrm D^s F(x)$.

*Proof.* $H(F(x\oplus_s\delta)) = H\big(F(x)\oplus_s(\mathrm D^sF\,\delta + O)\big) = H(F(x))\oplus_s\big(\mathrm D^sH\,\mathrm D^sF\,\delta\big) + O(\lVert\delta\rVert^2)$,
since $\oplus_s$ is smooth. $\square$

**Proposition LG.13 (first-order identities).** For $X \in G$ and $\tau, \delta \in \mathbb R^n$:

- (a) $\mathrm{Exp}(\delta)^{-1} = \mathrm{Exp}(-\delta)$ (exact).
- (b) $X\,\mathrm{Exp}(\delta) = \mathrm{Exp}(\mathrm{Ad}_X\delta)X$ and $\mathrm{Exp}(\delta)X = X\,\mathrm{Exp}(\mathrm{Ad}_X^{-1}\delta)$ (exact, LG.4).
- (c) $\mathrm{Exp}(\tau+\delta) = \mathrm{Exp}(\tau)\mathrm{Exp}(J_r\delta) + O = \mathrm{Exp}(J_l\delta)\mathrm{Exp}(\tau) + O$ (LG.8).
- (d) If $\theta < 2\pi$: $\mathrm{Exp}(\tau)\mathrm{Exp}(\delta) = \mathrm{Exp}(\tau + J_r^{-1}(\tau)\delta + O(\lVert\delta\rVert^2))$
  and $\mathrm{Exp}(\delta)\mathrm{Exp}(\tau) = \mathrm{Exp}(\tau + J_l^{-1}(\tau)\delta + O(\lVert\delta\rVert^2))$.
  If $\theta < \pi$: $\mathrm{Log}(\mathrm{Exp}(\tau)\mathrm{Exp}(\delta)) = \tau + J_r^{-1}(\tau)\delta + O(\lVert\delta\rVert^2)$
  and $\mathrm{Log}(\mathrm{Exp}(\delta)\mathrm{Exp}(\tau)) = \tau + J_l^{-1}(\tau)\delta + O(\lVert\delta\rVert^2)$.

*Proof.* (a) $\exp(-A)\exp(A) = I$. (d) $J_r$ is invertible for $\theta < 2\pi$ (LG.9), so $\mathrm{Exp}$ is a
local diffeomorphism at $\tau$ and $\mathrm{Exp}(\tau)\mathrm{Exp}(\delta) = \mathrm{Exp}(\tau + \eta(\delta))$ for a
smooth $\eta$ with $\eta(0) = 0$; by (c), $\mathrm{Exp}(\tau + \eta) = \mathrm{Exp}(\tau)\mathrm{Exp}(J_r\eta + O)$, so
$J_r\eta = \delta + O(\lVert\delta\rVert^2)$. For $\theta < \pi$ the point $\tau + \eta$ stays in
$\{\theta < \pi\}$, where $\mathrm{Log}\circ\mathrm{Exp} = \mathrm{id}$ (LG.2(c)). The left forms are identical
with $J_l$. $\square$

**Checked:** (a)–(c) are LG.4 and LG.8 above. (d): mpmath 1.4.1, 100 digits, each $N = 0..3$,
$\theta(\tau) \in \{0.5, 3.0, 3.5, 5, 6.2\}$, one random direction of $\delta$, $\lvert\delta\rvert = 10^{-6}, 10^{-7}$.
$\mathrm{Exp}$ forms (both sides): the error falls by exactly $100$ between the two $\lvert\delta\rvert$ at every
$\theta$, i.e. $O(\lvert\delta\rvert^2)$; at $10^{-6}$ it is $\le 5.6\times10^{-12}$ for $\theta \le 5$ and
$1.6\times10^{-8}$ at $\theta = 6.2$ (the constant grows with $\lVert J^{-1}\rVert^2$, LG.16). $\mathrm{Log}$ form
(right): $O(\lvert\delta\rvert^2)$, $\le 4.0\times10^{-13}$, for $\theta \le 3.0$; an $\lvert\delta\rvert$-independent error of $4.2$
to $6.3$ for $\theta \ge 3.5$. The $\mathrm{Log}$ and $\ominus$ rows of LG.14 exercise (d) on $\mathrm{Log}$ itself.
**Permanent:** as LG.14.

**The Log identities need $\theta < \pi$, not $\theta < 2\pi$.** For $\pi < \theta(\tau) < 2\pi$ the
$\mathrm{Exp}$ form of (d) still holds, but $\mathrm{Log}(\mathrm{Exp}(\tau + \eta))$ is the canonical
representative of $\mathrm{Exp}(\tau + \eta)$, not $\tau + \eta$ (LG.2(c)): the $\mathrm{Log}$ form fails by
$O(1)$. At $\theta(\tau) = \pi$ it fails for the reason of LG.2(d). In the rows of LG.14 the base
tangent is itself a $\mathrm{Log}$ output, hence canonical ($\theta \le \pi$): the window
$(\pi, 2\pi)$ never arises there, and the only excluded point is $\theta = \pi$. It arises when a caller
feeds an arbitrary $\tau$ with $\theta(\tau) > \pi$ into a $\oplus$-then-$\mathrm{Log}$ chain.

## 5. The sides table

**Proposition LG.14 (`NUMERICS.md` §2.3).** Let $E = \mathrm{Exp}(\tau)$, with $\tau$ the given tangent
for $\oplus$ and $\mathrm{Exp}$, $\tau = Y \ominus_s X$ for the $\ominus$ rows and $\tau = \mathrm{Log}\,X$
for the $\mathrm{Log}$ rows. Each row rewrites the perturbed operation as the unperturbed one with
$\delta$ moved to the outer end on the same side (Definition LG.11); the last column is $\mathrm D^s$.
Domain: $\theta(\tau) < \pi$ for the $\ominus$ and $\mathrm{Log}$ rows ($\tau$ is a $\mathrm{Log}$ output there, so
$\theta = \pi$ is the only excluded point, LG.13); every $\tau$ for the others.

| Operation, side | Perturbation | Rewriting | By | $\mathrm D^s$ |
|---|---|---|---|---|
| $X \oplus_R \tau$, w.r.t. $X$ | $X \to X\mathrm{Exp}\,\delta$ | $X\mathrm{Exp}(\delta)E = XE\,\mathrm{Exp}(\mathrm{Ad}_E^{-1}\delta)$ | LG.13(b) | $\mathrm{Ad}_E^{-1}$ |
| $X \oplus_R \tau$, w.r.t. $\tau$ | $\tau \to \tau + \delta$ | $X\,\mathrm{Exp}(\tau+\delta) = XE\,\mathrm{Exp}(J_r\delta)$ | LG.13(c) | $J_r(\tau)$ |
| $X \oplus_L \tau$, w.r.t. $X$ | $X \to \mathrm{Exp}\,\delta\,X$ | $E\,\mathrm{Exp}(\delta)X = \mathrm{Exp}(\mathrm{Ad}_E\delta)\,EX$ | LG.13(b) | $\mathrm{Ad}_E$ |
| $X \oplus_L \tau$, w.r.t. $\tau$ | $\tau \to \tau + \delta$ | $\mathrm{Exp}(\tau+\delta)X = \mathrm{Exp}(J_l\delta)\,EX$ | LG.13(c) | $J_l(\tau)$ |
| $Y \ominus_R X$, w.r.t. $Y$ | $Y \to Y\mathrm{Exp}\,\delta$ | $\mathrm{Log}(X^{-1}Y\mathrm{Exp}\,\delta) = \mathrm{Log}(E\,\mathrm{Exp}\,\delta) = \tau + J_r^{-1}\delta$ | LG.13(d) | $J_r^{-1}(\tau)$ |
| $Y \ominus_R X$, w.r.t. $X$ | $X \to X\mathrm{Exp}\,\delta$ | $\mathrm{Log}(\mathrm{Exp}(-\delta)X^{-1}Y) = \mathrm{Log}(\mathrm{Exp}(-\delta)E) = \tau - J_l^{-1}\delta$ | LG.13(a), (d) | $-J_l^{-1}(\tau)$ |
| $Y \ominus_L X$, w.r.t. $Y$ | $Y \to \mathrm{Exp}\,\delta\,Y$ | $\mathrm{Log}(\mathrm{Exp}(\delta)YX^{-1}) = \mathrm{Log}(\mathrm{Exp}(\delta)E) = \tau + J_l^{-1}\delta$ | LG.13(d) | $J_l^{-1}(\tau)$ |
| $Y \ominus_L X$, w.r.t. $X$ | $X \to \mathrm{Exp}\,\delta\,X$ | $\mathrm{Log}(YX^{-1}\mathrm{Exp}(-\delta)) = \mathrm{Log}(E\,\mathrm{Exp}(-\delta)) = \tau - J_r^{-1}\delta$ | LG.13(a), (d) | $-J_r^{-1}(\tau)$ |
| $XY$, R, w.r.t. $X$ | $X \to X\mathrm{Exp}\,\delta$ | $X\mathrm{Exp}(\delta)Y = XY\,\mathrm{Exp}(\mathrm{Ad}_Y^{-1}\delta)$ | LG.13(b) | $\mathrm{Ad}_Y^{-1}$ |
| $XY$, R, w.r.t. $Y$ | $Y \to Y\mathrm{Exp}\,\delta$ | $XY\,\mathrm{Exp}(\delta)$ | — | $I$ |
| $XY$, L, w.r.t. $X$ | $X \to \mathrm{Exp}\,\delta\,X$ | $\mathrm{Exp}(\delta)\,XY$ | — | $I$ |
| $XY$, L, w.r.t. $Y$ | $Y \to \mathrm{Exp}\,\delta\,Y$ | $X\mathrm{Exp}(\delta)Y = \mathrm{Exp}(\mathrm{Ad}_X\delta)\,XY$ | LG.13(b) | $\mathrm{Ad}_X$ |
| $X^{-1}$, R | $X \to X\mathrm{Exp}\,\delta$ | $\mathrm{Exp}(-\delta)X^{-1} = X^{-1}\mathrm{Exp}(-\mathrm{Ad}_X\delta)$ | LG.13(a), (b) at $X^{-1}$ | $-\mathrm{Ad}_X$ |
| $X^{-1}$, L | $X \to \mathrm{Exp}\,\delta\,X$ | $X^{-1}\mathrm{Exp}(-\delta) = \mathrm{Exp}(-\mathrm{Ad}_X^{-1}\delta)X^{-1}$ | LG.13(a), (b) at $X^{-1}$ | $-\mathrm{Ad}_X^{-1}$ |
| $\mathrm{Exp}(\tau)$, R | $\tau \to \tau + \delta$ | $\mathrm{Exp}(\tau+\delta) = E\,\mathrm{Exp}(J_r\delta)$ | LG.13(c) | $J_r(\tau)$ |
| $\mathrm{Exp}(\tau)$, L | $\tau \to \tau + \delta$ | $\mathrm{Exp}(\tau+\delta) = \mathrm{Exp}(J_l\delta)\,E$ | LG.13(c) | $J_l(\tau)$ |
| $\mathrm{Log}(X)$, R | $X \to X\mathrm{Exp}\,\delta$ | $\mathrm{Log}(\mathrm{Exp}(\tau)\mathrm{Exp}\,\delta) = \tau + J_r^{-1}\delta$ | LG.13(d) | $J_r^{-1}(\mathrm{Log}\,X)$ |
| $\mathrm{Log}(X)$, L | $X \to \mathrm{Exp}\,\delta\,X$ | $\mathrm{Log}(\mathrm{Exp}\,\delta\,\mathrm{Exp}(\tau)) = \tau + J_l^{-1}\delta$ | LG.13(d) | $J_l^{-1}(\mathrm{Log}\,X)$ |

Vector inputs ($\tau$ of $\oplus$ and $\mathrm{Exp}$) are perturbed additively on both sides, and vector
outputs ($\ominus$, $\mathrm{Log}$) are compared by subtraction (Definition LG.11); the side enters only
through the group-valued argument or result. These are exactly the entries of `NUMERICS.md` §2.3.

**Corollary LG.15 (chain-rule consistency).** Since $X \oplus_s \tau$ is a composition with
$\mathrm{Exp}(\tau)$ and $Y \ominus_s X$ is $\mathrm{Log}$ of a composition with an inverse, LG.12 derives
the $\oplus$ and $\ominus$ rows from the composition, inverse, $\mathrm{Exp}$ and $\mathrm{Log}$ rows. For
example, with $E = X^{-1}Y$ (right) and $E = YX^{-1}$ (left):
$$
\frac{\partial(Y\ominus_RX)}{\partial X} = J_r^{-1}(\tau)\,\mathrm{Ad}_Y^{-1}\,(-\mathrm{Ad}_X) = -J_r^{-1}\mathrm{Ad}_E^{-1} = -(\mathrm{Ad}_EJ_r)^{-1} = -J_l^{-1},
$$
$$
\frac{\partial(Y\ominus_LX)}{\partial X} = J_l^{-1}(\tau)\,\mathrm{Ad}_Y\,(-\mathrm{Ad}_X^{-1}) = -J_l^{-1}\mathrm{Ad}_E = -(\mathrm{Ad}_E^{-1}J_l)^{-1} = -J_r^{-1},
$$
by $\mathrm{Ad}_Y^{-1}\mathrm{Ad}_X = \mathrm{Ad}_{Y^{-1}X} = \mathrm{Ad}_E^{-1}$ (resp. $\mathrm{Ad}_Y\mathrm{Ad}_X^{-1} = \mathrm{Ad}_E$)
and LG.9(b). The remaining $\oplus$ and $\ominus$ rows follow the same way. This chain is the reference
twin of `*_jacobians` (`NUMERICS.md` §14). $\square$

**Checked (LG.14, LG.15):** mpmath 1.4.1, 110 digits. $\mathrm{Exp}$ by `mp.expm`; $\mathrm{Log}$ by the
geometric $\mathrm{SO}(3)$ inverse ($\mathrm{atan2}$ of the antisymmetric part and the trace) and a
$\Gamma_1(\varphi)$ solve; $\mathrm{Ad}$ and $J_{r,l}$ from Definitions LG.3 and LG.7; every inverse (group,
$\mathrm{Ad}$, $J$) by `mp.inverse`. $\mathrm D^s$ by the second form of Definition LG.11, the operation itself
perturbed column by column: central differences $[F(he_j) - F(-he_j)]/2h$, $h = 10^{-30}$, so the figures
are limited by the $O(h^2) = 10^{-60}$ truncation, not by the formulas. Each $N = 0..3$: 12 random
$(X, Y, \tau)$ ($X, Y$ as in LG.4, $\theta(\tau) \le 6$, rejected if $\theta(X^{-1}Y)$ or $\theta(X) > \pi - 0.02$); plus one
sample per stratum: $\theta(X^{-1}Y) = \pi - 10^{-3}$ (so $\theta(\tau) = \pi - 10^{-3}$ for the $\ominus$ rows),
$\theta(X) = \pi - 10^{-3}$ (the $\mathrm{Log}$ rows), $Y = X$ ($\tau = 0$), $\theta(\tau) = 10^{-10}$, $x_i, \rho_i \sim 10^6$;
and, for the domain-free rows only, $\theta(\tau) \in \{\pi + 10^{-3}, 2\pi - 10^{-3}, 6.5, 20\}$. All 18 rows, both
sides. Max discrepancy (absolute for entries $\le 1$, else relative to the largest predicted entry), over
all rows, sides and strata: $4.5\times10^{-62}$ (SO(3)), $1.3\times10^{-61}$ (SE(3)), $1.2\times10^{-61}$
(SE₂(3)), $1.1\times10^{-61}$ (SE₃(3)); the $J$-valued rows ($\oplus$ and $\mathrm{Exp}$ w.r.t. $\tau$, $\ominus$,
$\mathrm{Log}$) carry all of it, the $\mathrm{Ad}$-valued ones are below $2\times10^{-75}$. The chain of LG.15,
max $1.4\times10^{-108}$. Script not committed. **Permanent:** planned, proptests `jacobians_match_dual_*`
and `*_matches_reference`; corpus `sen3_jr_n*`, `sen3_jl_n*`, `sen3_jr_inv_n*`, `sen3_jl_inv_n*` for the
$J$ entries the rows use.

## 6. Sign audit

Where a sign or a side is easy to lose: the correct form, the usual wrong one, and the **exact gap**
between them ($\zeta$ is the tangent the entry depends on: $\tau$, $\mathrm{Log}\,X$ or $\mathrm{Log}\,Y$).
The gap of a swapped side or a swapped $J$ or $\mathrm{Ad}$ pair is a function of $\mathrm{ad}_\zeta$ that vanishes with
$\zeta$, so a test of a side needs $\lVert\zeta\rVert$ of order one; for $\mathrm{SO}(3)$ the $\mathrm{Ad}$ swaps also vanish at
$\theta(\zeta) = \pi$, where $R = R^\top$.

| Where | Correct | Wrong variant | Gap | Reason |
|---|---|---|---|---|
| $J_r$ vs $J_l$ | $J_r = \sum(-\mathrm{ad})^k/(k+1)!$ ($\mathrm{SO}(3)$: $I - aW + bW^2$, $a, b$ of `NUMERICS.md` §3.5, §4); $\mathrm{Exp}(\tau+\delta) = \mathrm{Exp}(\tau)\mathrm{Exp}(J_r\delta)$ | the two swapped, in $\mathrm{Exp}$, $\oplus$ or $\mathrm{Log}$ | $J_l - J_r = \mathrm{ad}_\tau\sum_k 2\,\mathrm{ad}_\tau^{2k}/(2k+2)! = \mathrm{ad}_\tau + O(\theta^3)$; $J_r^{-1} - J_l^{-1} = \mathrm{ad}_\tau$ exactly | LG.7, LG.8, LG.9(a) |
| $\ominus$ w.r.t. the subtrahend | R: $-J_l^{-1}(\tau)$; L: $-J_r^{-1}(\tau)$ | $-J_r^{-1}$ (R); or $+J_l^{-1}$ | $\mathrm{ad}_\tau$; $2J_l^{-1}$, of 2-norm $\ge 2$ (LG.16(b)) | $\mathrm{Exp}(-\delta)$ flips the sign **and** the side (LG.14) |
| $\ominus$ w.r.t. the minuend | R: $J_r^{-1}(\tau)$; L: $J_l^{-1}(\tau)$ | $J_l^{-1}$ (R) | $\mathrm{ad}_\tau$ | LG.13(d) |
| $\oplus$ w.r.t. the base | R: $\mathrm{Ad}_{\mathrm{Exp}\,\tau}^{-1}$; L: $\mathrm{Ad}_{\mathrm{Exp}\,\tau}$ | swapped; or $I$ | $2\sinh\mathrm{ad}_\tau = 2\,\mathrm{ad}_\tau + O(\theta^3)$; $\mp\mathrm{ad}_\tau + O(\theta^2)$ | LG.13(b) |
| $XY$ | R w.r.t. $X$: $\mathrm{Ad}_Y^{-1}$; L w.r.t. $Y$: $\mathrm{Ad}_X$ | $\mathrm{Ad}_Y$, $\mathrm{Ad}_X^{-1}$ | $2\sinh\mathrm{ad}_{\mathrm{Log}\,Y}$, $2\sinh\mathrm{ad}_{\mathrm{Log}\,X}$ | the moved factor is conjugated by the far factor (LG.13(b)) |
| $X^{-1}$ | R: $-\mathrm{Ad}_X$; L: $-\mathrm{Ad}_X^{-1}$ | swapped | $2\sinh\mathrm{ad}_{\mathrm{Log}\,X}$ | LG.14 |
| $\mathrm{Ad}^{-1}$ vs $\mathrm{Ad}^\top$ | $\mathrm{Ad}_X^{-1} = \mathrm{Ad}_{X^{-1}}$, exact | $\mathrm{Ad}_X^\top$ | $\pm R^\top[x_i]_\times$ in blocks $(i,0)$ and $(0,i)$: zero iff every $x_i = 0$ | block $(0,0)$ of $\mathrm{Ad}^\top\mathrm{Ad}$ is $I + \sum_i C_i^\top C_i$, $C_i = [x_i]_\times R$, so $\mathrm{Ad}_X$ is orthogonal iff every $x_i = 0$; $\mathrm{Ad}^{-1}$ is block lower-triangular ($-R^\top[x_i]_\times$ in $(i,0)$), $\mathrm{Ad}^\top$ block upper-triangular (same block in $(0,i)$) |
| $J_r(\tau)$ from $Q$ | $Q(-\rho, -\varphi)$: negate the **whole** tangent | negate only $\varphi$ or only $\rho$ | at $\tau = (0.7, -0.4, 1.1; 0.3, 1.2, -0.8)$: $0.778$ (only $\varphi$), $1.00$ (only $\rho$), $0.958$ (neither: the $J_l$ block) | LG.9(a); `NUMERICS.md` §5.3 (derived in [`se3.md`](./se3.md), SE.8–SE.9; measured here, see below) |
| Block position | lower-left: $\mathrm{Ad}_X$ column $[x_i]_\times R$, $\mathrm{ad}_\tau$ column $[\rho_i]_\times$ | upper-right (literature order); $-[x_i]_\times R$ | the block $C_i$ moves ($\ne 0$ iff $x_i \ne 0$); the sign variant differs by $2C_i$ | LG.4(c): $-(R\varphi)\times x_i = +[x_i]_\times R\varphi$ |
| Domain of $J^{-1}$ | exists for $\theta < 2\pi$ (LG.9) | assuming $\mathrm{Log}$ is differentiable at $\theta(X) = \pi$ | none: $\mathrm{Log}$ jumps there, by $2\pi$ in the rotation block | LG.2(d) |
| $\mathrm{Log}$ identity | $\mathrm{Log}(\mathrm{Exp}(\tau)\mathrm{Exp}(\delta)) = \tau + J_r^{-1}\delta$ needs $\theta(\tau) < \pi$ | using it for $\pi < \theta < 2\pi$ | $O(1)$: $4.2$ to $6.3$ measured at $\theta \ge 3.5$ (LG.13) | LG.13(d) |

**Checked:** mpmath 1.4.1. 110 digits, the $\tau$ and $X$ of LG.2, LG.4, each $N = 0..3$ ($N \ge 1$ for the
$\mathrm{Ad}^\top$ rows; inverses for $\theta < 2\pi - 10^{-3}$): $J_r^{-1} - J_l^{-1} = \mathrm{ad}_\tau$, max $9.2\times10^{-107}$;
the series for $J_l - J_r$, $1.1\times10^{-104}$; block $(0,0)$ of $\mathrm{Ad}^\top\mathrm{Ad}$, $1.5\times10^{-98}$; block
positions of $\mathrm{Ad}^{-1}$ and $\mathrm{Ad}^\top$, $3.5\times10^{-105}$. 100 digits, fixed vectors: the SE(3) $J_r$
block against $Q(-\rho,-\varphi)$, $1.4\times10^{-101}$, and the three wrong variants of the $Q$ row at
$\tau = (0.7, -0.4, 1.1; 0.3, 1.2, -0.8)$; for $R = \mathrm{Exp}(0.4, -0.3, 0.5)$, $x = (1.0, -2.0, 0.5)$: the lower-left
block of $\mathrm{Ad}_X$ minus $[x]_\times R$, $2.9\times10^{-101}$, plus it, $4.15$, the upper-right block exactly
$0$, and $\lVert\mathrm{Ad}^{-1} - \mathrm{Ad}^\top\rVert_{\max} = \lVert R^\top[x]_\times\rVert_{\max} = 2.08$; $R - R^\top$ at
$\theta = \pi$, $2.9\times10^{-102}$. Script not committed. **Permanent:** as LG.14.

## 7. Numerical consequences

**Proposition LG.16 (conditioning of $J$ and $\mathrm{Ad}$).** Let $0 < \theta < 2\pi$,
$s(\theta) = 2\sin(\theta/2)/\theta$, and $\kappa_2$ the 2-norm condition number.

- (a) $N = 0$: the singular values of $J_r(\varphi)$ and $J_l(\varphi)$ are $1, s, s$, so
  $\kappa_2 = \lVert J^{-1}\rVert_2 = 1/s(\theta) = \theta/(2\sin(\theta/2))$: $1 + \theta^2/24 + O(\theta^4)$ near $0$,
  $\pi/2$ at $\pi$, $\sim 2\pi/(2\pi - \theta)$ as $\theta \to 2\pi$.
- (b) Every $N$: $\lVert J^{-1}(\tau)\rVert_2 \ge 1/s(\theta)$ and $\kappa_2(J(\tau)) \ge 1/s(\theta)$.
- (c) $M = \begin{bmatrix}I & 0\\ B & I\end{bmatrix}$ ($B$ of size $3N \times 3$) has singular values
  $(\sqrt{\beta^2 + 4} \pm \beta)/2$ for each singular value $\beta$ of $B$, and $1$ otherwise; so
  $\kappa_2(M) = \lVert M\rVert_2^2$ with $\lVert M\rVert_2 = (\beta_{\max} + \sqrt{\beta_{\max}^2 + 4})/2$. At $\varphi = 0$,
  $J_{r,l}(\tau) = I \mp \tfrac12\mathrm{ad}_\tau$ is such an $M$ with $B = \mp\tfrac12\,[[\rho_1]_\times;\dots;[\rho_N]_\times]$,
  so $\kappa_2 \approx \beta_{\max}^2$, which is $(\lVert\rho\rVert/2)^2$ for $N = 1$, for large $\rho$. Likewise $\mathrm{Ad}_X = M\,\mathrm{diag}(R,\dots,R)$ with
  $B = [[x_1]_\times;\dots;[x_N]_\times]$, $\sigma = \lVert B\rVert_2 \in [\max_i\lVert x_i\rVert, (\sum_i\lVert x_i\rVert^2)^{1/2}]$:
  $\lVert\mathrm{Ad}_X\rVert_2 \le 1 + \sigma$ and $\kappa_2(\mathrm{Ad}_X) = \lVert\mathrm{Ad}_X\rVert_2^2$.
- (d) For $J = A + \epsilon B$, the $(i,0)$ block of $J^{-1}$ is $-A^{-1}B_iA^{-1}$ (LG.10), of norm at most
  $\lVert B_i\rVert/s(\theta)^2$, and
  $B_i$ is linear in $\rho_i$ (LG.4(c)): the dual algebra never forms or solves the dense $J$, but that block
  still scales like $\lVert\rho_i\rVert$, and like $(2\pi/(2\pi-\theta))^2$ near $2\pi$.

*Proof.* (a) $W$ is skew, hence normal, so $g(W)$ and $g(-W)$ are normal with eigenvalues $1$ and
$g(\pm i\theta)$ of modulus $s$ (proof of LG.2(b)); the singular values of a normal matrix are the moduli of
its eigenvalues. (b) By LG.10 the diagonal blocks of $J$ and $J^{-1}$ are $J(\varphi)$ and $J(\varphi)^{-1}$, and a
submatrix has no larger 2-norm; $\lVert J\rVert_2 \ge \lVert J(\varphi)\rVert_2 = 1$. (c) With $B = U\Sigma V^\top$,
$\mathrm{diag}(V^\top, U^\top)\,M\,\mathrm{diag}(V, U) = \begin{bmatrix}I & 0\\ \Sigma & I\end{bmatrix}$, a direct sum of
$\begin{bmatrix}1 & 0\\ \beta & 1\end{bmatrix}$, whose singular values have product $1$ and squares summing to
$2 + \beta^2$. $\mathrm{Ad}_X = M\,\mathrm{diag}(R,\dots,R)$ by LG.4(c), the orthogonal factor leaving singular values
unchanged. (d) Submultiplicativity, with $\lVert A^{-1}\rVert_2 = 1/s$ by (a). $\square$

**Checked:** mpmath 1.4.1, 100 digits, singular values by `mp.svd_r`. (a) $\theta \in \{10^{-3}, 0.5, 1, 2, 3, \pi, 4, 5, 6,
2\pi - 10^{-3}\}$: $\kappa_2$ and $\{1, s, s\}$, relative max $3.7\times10^{-97}$. (b), (d) ($(1,0)$ block) $N = 1, 2, 3$,
$\theta \in \{0.1, 1, \pi, 5, 2\pi - 10^{-2}\}$, $\rho_i \sim \mathcal N(0, 2^2)$ and $\mathcal N(0, 100^2)$: 30 cases, no
violation. (c) $N = 1, 2, 3$, $\rho_i, x_i \sim \mathcal N(0, c^2)$ with $c \in \{2, 100, 10^6\}$: $J_{r,l}(0,\rho)$ equals
$I \mp \mathrm{ad}_\tau/2$ exactly, and $\lVert\cdot\rVert_2$, $\kappa_2$ of $J$ and $\mathrm{Ad}_X$ match the formulas to
$1.9\times10^{-91}$ and $2.6\times10^{-92}$. Illustration, SE(3), $\tau = (0, 0, \theta; 100, 0, 0)$: $\kappa_2(J_r) = 2.5\times10^3$
($\theta = 10^{-8}$), $2.2\times10^3$ ($\theta = \pi$), $1.6\times10^6$ ($\theta = 2\pi - 10^{-3}$, where (b) gives
$\ge 6.3\times10^3$). Script not committed. **Permanent:** none specified; the rounding error of the closed forms
and of the dual algebra is measured by the corpus (`NUMERICS.md` §11), not here.

- The series of Definition LG.7 are **definitions and generator formulas** (`PHASE1.md` §4.3), not
  algorithms. Code evaluates the closed forms of `NUMERICS.md` §3.5, §4 and §5.3, whose coefficients
  cancel near $\theta = 0$; that analysis is not on this page.
- `jl` as `jr(-τ)` is a sign flip (LG.9(a)) and costs nothing; $\mathrm{Ad}_{\mathrm{Exp}\,\tau}J_r$ (LG.9(b))
  is the dense reference twin, not a cheaper route.
- $\mathrm{Ad}_X^{-1} = \mathrm{Ad}_{X^{-1}}$ (LG.4(b)): the inverse of an adjoint is another adjoint, so no
  matrix inverse or solve is needed. In the dual algebra (LG.10), $J^{-1}$ costs $A^{-1}$ (closed
  form) plus $2N$ products of $3\times3$ blocks, never a $3(N+1)$-square inversion.
- Two different domains: $J^{-1}$ exists for $\theta < 2\pi$ (`NUMERICS.md` §12), but the *derivative*
  of $\mathrm{Log}$ (hence of $\ominus$) is $J^{-1}(\mathrm{Log}\,X)$ only for $\theta(X) < \pi$. $J^{-1}$ stays
  bounded as $\theta \to \pi$ (LG.16(a): $\pi/2$) while $\mathrm{Log}$ jumps across $\theta(X) = \pi$ (LG.2(d)), so
  the derivative describes one side of the cut. `NUMERICS.md` §12 has no row for it (open item in the
  [maths index](./index.md)).
