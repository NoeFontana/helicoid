# SO(3): quaternion, Exp, Log, Jacobians, `from_matrix`

> Non-normative companion to [`NUMERICS.md`](../NUMERICS.md) §2.4 (SO(3) rows), §3 and the
> *definitions* of $k, a, b, c, r$ in §4 (their series, switch points and cancellation analysis
> beyond the amplification factors quoted here are not derived). **`NUMERICS.md` wins on any
> conflict; a disagreement is an open item in the [maths index](./index.md)**, which also holds the
> notation and the `Checked:` convention; the general algebra ($\mathrm{Ad}$,
> $\mathrm{ad}$, $J = g(\pm\mathrm{ad})$, the sides table) is in [`lie-groups.md`](./lie-groups.md).

On this page $q = (w, u)$, $w \in \mathbb R$, $u \in \mathbb R^3$ (the vector part, as in
`NUMERICS.md` §1, §3.3), $n = \lVert u\rVert$, $\hat n = u/n$ (so $u = n\hat n$),
$\eta = \lVert q\rVert^2 - 1$, $\theta = \lVert\varphi\rVert$, $W = [\varphi]_\times$,
$\mathcal A = \operatorname{span}\{I, W, W^2\}$. The unit roundoff is written $\mathsf u$ ($2^{-53}$
for `f64`, $2^{-24}$ for `f32`), because $u$ is the vector part; likewise $n$ is $\lVert u\rVert$ here,
not the tangent dimension of the index.

## Results

| Label | Result | `NUMERICS.md` | Implemented by |
|---|---|---|---|
| SO.2 | $qvq^* = \lVert q\rVert^2 v + 2w\,u\times v + 2\,u\times(u\times v)$; $R(q)$; $R(pq) = R(p)R(q)$; kernel $\{\pm1\}$; non-unit $q$ | §1, §3.3 | `SO3::act`, `act_many`, `to_matrix`, `Mul` |
| SO.4 | $\exp W = I + 2k\cos\tfrac\theta2\,W + 2k^2W^2 = R\big(\cos\tfrac\theta2, k\varphi\big)$; the corpus's power series | §3.1, §4 ($k$, $a$) | `SO3::exp` |
| SO.5 | $\mathrm{Log}$ by $\mathrm{atan2}$, the $w<0$ flip (off by exactly $2\pi\hat n$ without it), scale invariance, $w = 0$ | §3.2, §12 | `SO3::log` |
| SO.6 | why not $\arccos\frac{\operatorname{tr}R-1}{2}$: $\delta\theta\approx\delta\chi/\sin\theta$, floor $\sqrt{2\mathsf u}$; $\mathrm{atan2}$: $\delta\theta = \sin\theta\,(\xi_n - \xi_w)$; range and scaling | §3.2, §11 | (the design of `log`) |
| SO.7 | $J_r = I - aW + bW^2$, $J_l = I + aW + bW^2 = J_r^\top = RJ_r$ | §3.5, §4 ($a$, $b$) | `jr`, `jl` |
| SO.8 | $J_r^{-1} = I + \tfrac12W + cW^2$ solved in $\mathcal A$; $c = \theta^{-2} - \cot(\theta/2)/(2\theta)$; $\theta\notin2\pi\mathbb Z_{>0}$ | §3.5, §4 ($c$), §12 | `jr_inv`, `jl_inv` |
| SO.9 | $\mathrm{Ad}_R = R$, $\mathrm{ad}_\varphi = W$ | §3.5 | `adjoint`, `ad` |
| SO.10 | $\partial(Rp)/\partial R$ (both sides), $\partial/\partial p$; the SE(3) rows | §2.4 | `act_jacobians` |
| SO.11–SO.13 | $4qq^\top = K(R) + I$; Shepperd's four extractions, pivot rule, error bound $\tfrac74\varepsilon$, the output sign; non-orthogonal input; nearest rotation | §3.4 | `SO3::from_matrix` |
| SO.14–SO.15 | Newton step: $\lVert q'\rVert^2 - 1 = -\tfrac34\eta^2 + \tfrac14\eta^3$ exactly; the $2^{-40}$, $2^{-16}$ bounds | §3.6 | `renormalize`, `Quat::from_wxyz_unchecked` |

## 1. Quaternions and $R(q)$

**Definition SO.1.** On $\mathbb H = \mathbb R^4$,
$(w,u)(w',u') = (ww' - u\cdot u',\ wu' + w'u + u\times u')$, $q^* = (w,-u)$,
$qq^* = \lVert q\rVert^2$; $v \in \mathbb R^3$ is the pure quaternion $(0,v)$. This is Hamilton's
$i^2 = j^2 = k^2 = ijk = -1$ (`NUMERICS.md` §1). Unit quaternions $S^3$ form a group.

**Proposition SO.2.** For every $q$ and $v$, $qvq^*$ is pure and
$$
qvq^* = \lVert q\rVert^2 v + 2w\,u\times v + 2\,u\times(u\times v) = R_\square(q)\,v,\qquad
R_\square(q) = (w^2 - n^2)I + 2uu^\top + 2w[u]_\times. \tag{1}
$$
For $\lVert q\rVert = 1$ this is $v + 2w\,u\times v + 2\,u\times(u\times v)$ (`NUMERICS.md` §3.3),
$R(q) := R_\square(q)$ is the matrix of `NUMERICS.md` §1 and equals
$I + 2w[u]_\times + 2[u]_\times^2$, and:

- (a) $R(pq) = R(p)R(q)$, $R(q^*) = R(q)^\top$: `a * b` is $T_{ax}T_{xb}$.
- (b) $R(q)u = u$, $\cos\theta = w^2 - n^2 = 2w^2 - 1$, $\sin\theta = 2wn$ for the rotation angle
  $\theta$ about $\hat n$, and $\operatorname{tr}R = 4w^2 - 1$. So $w = \cos\frac\theta2$,
  $n = \sin\frac\theta2$.
- (c) $q \mapsto R(q)$ is a surjective homomorphism $S^3 \to \mathrm{SO}(3)$ with kernel $\{\pm1\}$.
- (d) For $\lVert q\rVert \ne 1$ the vector form $v + 2w\,u\times v + 2\,u\times(u\times v)$ equals
  $R_\square(q)v + (1 - \lVert q\rVert^2)v$, whereas
  $R_\square(q) = \lVert q\rVert^2 R(q/\lVert q\rVert)$. Against the true rotation $R(\hat q)$,
  $\hat q = q/\lVert q\rVert$, the vector form errs by $\eta\,(R(\hat q) - I)v$ (norm
  $\le 2\lvert\eta\rvert\sin\frac\theta2\lVert v\rVert$) and $R_\square$ by $\eta R(\hat q)v$ (norm
  $\lvert\eta\rvert\lVert v\rVert$): **the two forms of `NUMERICS.md` §3.3 differ by $\eta v$ on a
  non-unit $q$.**

*Proof.* $qv = (-u\cdot v,\ wv + u\times v)$, and multiplying by $q^* = (w, -u)$ gives scalar part
$-w\,u\cdot v + w\,u\cdot v + (u\times v)\cdot u = 0$ and vector part
$(u\cdot v)u + w^2v + 2w\,u\times v + u\times(u\times v)$, using
$-(u\times v)\times u = u\times(u\times v)$. With $u\times(u\times v) = (u\cdot v)u - n^2v$ this is
(1) in both forms. (a): $(pq)v(pq)^* = p(qvq^*)p^*$. (b): (1) with $v = u$ gives $R u = u$; for
$v\perp u$, $Rv = (w^2 - n^2)v + 2wn\,\hat n\times v$, and $\hat n\times v$ is $v$ turned by
$90^\circ$ in the plane $\perp\hat n$: a rotation with cosine $w^2 - n^2$ and sine $2wn$. (c):
$\lVert qvq^*\rVert = \lVert q\rVert^2\lVert v\rVert$ makes $R(q)$ orthogonal; $\det R$ is
continuous on the connected $S^3$ and $1$ at $q = 1$. Surjective: a rotation by $\theta$ about
$\hat n$ is $R(\cos\frac\theta2, \sin\frac\theta2\hat n)$ by (b). $R(q) = I$ gives $4w^2 - 1 = 3$,
so $w = \pm1$, $u = 0$. (d): degree-2 homogeneity of $R_\square$ and the first form. $\square$

**Checked:** mpmath 1.3.0, 80 digits. 200 random unit $q$ and $\lVert v\rVert \le 3\sqrt3$ (plus the
non-unit $q = \lambda\hat q$, $\lambda\in[0.5,1.5]$ of (d)): $qvq^*$ by the definition of the
product against both forms of (1), orthogonality, $\det$, trace, homomorphism, $R(-q) = R(q)$,
$R(q)u = u$, $\cos\theta$, and the error formula of (d): max $2.3\times10^{-80}$. Script not
committed. **Permanent:** planned, corpus `so3_act`; the twin of `act_many` (`NUMERICS.md` §14), on
*normalized* inputs (see (d)).

## 2. Exp

**Lemma SO.3.** $W^2 = \varphi\varphi^\top - \theta^2 I$ and $W^3 = -\theta^2W$; so
$W^{2j+1} = (-\theta^2)^jW$, $W^{2j+2} = (-\theta^2)^jW^2$, and $\mathcal A$ is a commutative
algebra:
$$
(x_0 + x_1W + x_2W^2)(y_0 + y_1W + y_2W^2) = z_0 + z_1W + z_2W^2,\quad
\begin{aligned}z_0 &= x_0y_0,\\ z_1 &= x_0y_1 + x_1y_0 - \theta^2(x_1y_2 + x_2y_1),\\
z_2 &= x_0y_2 + x_2y_0 + x_1y_1 - \theta^2x_2y_2.\end{aligned} \tag{2}
$$
*Proof.* $W^2v = \varphi\times(\varphi\times v) = (\varphi\cdot v)\varphi - \theta^2v$ and
$W\varphi = 0$. (2) uses $W\cdot W = W^2$, $W\cdot W^2 = W^3$, $W^2W^2 = W^4 = -\theta^2W^2$. For
$\theta > 0$, $I, W, W^2$ are independent ($W$ is skew and nonzero, $W^2$ symmetric of rank 2).
$\square$

**Proposition SO.4.** For every $\varphi$, with $k = \sin(\theta/2)/\theta$ and
$a = (1-\cos\theta)/\theta^2$ (both $\tfrac12$ at $0$):

- (a) $\exp W = I + \dfrac{\sin\theta}{\theta}W + aW^2$, with
  $\dfrac{\sin\theta}\theta = 2k\cos\dfrac\theta2$ and $a = 2k^2$.
- (b) $\exp W = R(q)$ for the unit $q = (\cos\frac\theta2,\ k\varphi)$.
- (c) The quaternion series $\sum_m p^m/m!$ at $p = (0, \varphi/2)$ equals that $q$ (the corpus
  definition of `so3_exp`, `PHASE1.md` §4.3).

*Proof.* (a) Split $\sum_m W^m/m!$ by parity with SO.3: the odd part is
$W\sum_j(-1)^j\theta^{2j}/(2j+1)! = W\sin\theta/\theta$ and the even part
$I + W^2\sum_j(-1)^j\theta^{2j}/(2j+2)! = I + aW^2$. The half-angle forms are
$\sin\theta = 2\sin\frac\theta2\cos\frac\theta2$ and $1 - \cos\theta = 2\sin^2\frac\theta2$. (b) In
(1) put $w = \cos\frac\theta2$, $u = \sin\frac\theta2\,\hat n = k\varphi$:
$R = \cos\theta\,I + (1-\cos\theta)\hat n\hat n^\top + \sin\theta[\hat n]_\times$, which is (a)
since $[\hat n]_\times^2 = \hat n\hat n^\top - I$. (c) $p^2 = -\theta^2/4$ is scalar, so
$e^p = \cos\frac\theta2 + \frac{\sin(\theta/2)}{\theta/2}\,p = \cos\frac\theta2 + k\varphi$.
$\square$

*Consequences.* $k$ is $0/0$ only at $\theta = 0$. The rewrite $a = 2k^2$ removes the cancellation
of $(1 - \cos\theta)/\theta^2$ (terms of size $\theta^{-2}$ for a value near $\tfrac12$:
amplification $2/\theta^2$). Exp is defined for all $\varphi$ (its accuracy degrades linearly in $\theta$: *Range and scaling*, after SO.6) and does **not** canonicalize: for
$\theta\in(\pi,2\pi)$, $w = \cos\frac\theta2 < 0$ (the same rotation as
$\mathrm{Exp}((1 - 2\pi/\theta)\varphi)$, `lie-groups.md` LG.2(c)); `Log` picks the branch (SO.5).

**Checked:** mpmath 1.3.0, 80 digits, 72 tangents
($\theta\in\{0, 10^{-10}, 10^{-30}, 1, \pi - 10^{-6}, \pi, \pi+10^{-3}, 5, 2\pi - 10^{-3}, 2\pi, 7, 20\}$
and 60 uniform on $[0,7]$; axes uniform). `mp.expm` of $W$ against (a) $7.4\times10^{-81}$, against
$R(q)$ of (b) $7.4\times10^{-81}$, against $R$ of the series (c) $3.8\times10^{-78}$;
$\lVert q\rVert = 1$, $a = 2k^2$, $\sin\theta/\theta = 2k\cos\frac\theta2$ ($k$, $a$ evaluated with
80 guard digits at tiny $\theta$): max $4.4\times10^{-82}$. `f64` (Python floats), $q$ of (b) from
$\theta = \sqrt{x^2+y^2+z^2}$, $k = \sin(\theta/2)/\theta$, against the exact $q$ of the same double
$\varphi$ (300 random axes at each $\theta$): $\max\lVert q_{\mathsf f} - q\rVert \le 0.58\,\theta\mathsf u$
for $\theta = 20, 10^3, 10^6, 10^9, 10^{12}$ ($0.83\,\mathsf u$ at $\theta = 1$). Script not committed.
**Permanent:** planned, corpus `so3_exp`, `coeff_k`, `coeff_a`.

## 3. Log

**Proposition SO.5.** Let $q$ be a quaternion with $n = \lVert u\rVert$; define
$L(q) = \dfrac{2\,\mathrm{atan2}(n, w)}{n}\,u$ ($0$ at $u = 0$).

- (a) If $\lVert q\rVert = 1$ and $w \ge 0$: $\theta(L) = 2\,\mathrm{atan2}(n,w) \in [0,\pi]$ and
  $\mathrm{Exp}\,L(q) = R(q)$; for $w > 0$, $L(q)$ is the canonical $\mathrm{Log}\,R(q)$ (the unique
  preimage in $\theta < \pi$, `lie-groups.md` LG.2(c)).
- (b) If $w < 0$ and $n > 0$, $L(q) = \mathrm{Log}(-q) + 2\pi\,\hat n$ ($\hat n = u/n$ of $q$, not of
  $-q$): the same rotation, $\theta > \pi$. The flip $q \leftarrow -q$ for $w < 0$ is what makes $L$
  canonical. (At $q = -1$, $n = 0$: $L = 0 = L(1)$, the identity, with no $\hat n$ to shift by.)
- (c) $L(\lambda q) = L(q)$ for $\lambda > 0$: for non-unit $q$, $L$ is exactly the $\mathrm{Log}$
  of $q/\lVert q\rVert$.
- (d) At $w = 0$, $\lVert q\rVert = 1$: $R(q) = 2\hat n\hat n^\top - I$ depends on $\hat n$ up to
  sign, so $\pm\pi\hat n$ are both logarithms. $\mathrm{Log}$ jumps by $2\pi$ across $w = 0$ (no
  continuous choice, `lie-groups.md` LG.2(d), $N = 0$); `NUMERICS.md` §3.2's "no flip at $w = +0$"
  returns $+\pi\hat n$ for the $u$ given, a function of $q$, not of $R$.

*Proof.* (a) $\alpha = \mathrm{atan2}(n,w) \in [0,\frac\pi2]$, and
$(w, n) = (\cos\alpha, \sin\alpha)$ for unit $q$, so $q = (\cos\alpha, \sin\alpha\,\hat n)$ and
SO.4(b) with $\theta = 2\alpha$, $\varphi = 2\alpha\hat n$ gives $R(q) = \mathrm{Exp}(L)$. (b)
$\mathrm{atan2}(n, -w) = \pi - \alpha'$ for $\alpha' = \mathrm{atan2}(n,w)$, so
$\mathrm{Log}(-q) = -\frac{2(\pi - \alpha')}{n}u$ and $L(q) - \mathrm{Log}(-q) = \frac{2\pi}{n}u$.
(c) $\mathrm{atan2}$ is homogeneous of degree $0$ and $u/n$ is invariant. (d) (1) at $w = 0$.
$\square$

**Proposition SO.6 (why not $\arccos$).** Let $\chi = \cos\theta = (\operatorname{tr}R - 1)/2$,
$0 < \theta < \pi$.

- (a) $\lvert d\theta/d\chi\rvert = 1/\sin\theta$, and
  $\arccos(1 - \delta) = \sqrt{2\delta}\,(1 + O(\delta))$,
  $\arccos(-1 + \delta) = \pi - \sqrt{2\delta}\,(1 + O(\delta))$. A perturbation $\delta\chi$
  therefore moves $\theta$ by $\approx\lvert\delta\chi\rvert/\sin\theta$ and, once
  $\theta^2 \lesssim \lvert\delta\chi\rvert$ (or $(\pi-\theta)^2 \lesssim \lvert\delta\chi\rvert$), by
  $\sqrt{2\lvert\delta\chi\rvert}$.
- (b) A double near $\pm1$ is uncertain by $\lvert\delta\chi\rvert \le \mathsf u/2$ (spacing $2^{-53}$
  on $[\frac12, 1)$) from rounding $\chi$ alone. So any $\theta$ read from $\chi$ has relative error
  $\approx \mathsf u/(2\theta^2)$ near $0$ and $\approx\mathsf u/(2\pi(\pi - \theta))$ near $\pi$: 8
  digits lost at $\theta = 10^{-4}$, none left below
  $\theta \approx \sqrt{\mathsf u} = 1.0\times10^{-8}$, and the resolution floor is
  $\arccos(1 - \mathsf u) = \sqrt{2\mathsf u}$ = $1.5\times10^{-8}$ (`f32`: $3.5\times10^{-4}$).
- (c) For unit $q$, relative perturbations $n \to n(1 + \xi_n)$, $w \to w(1 + \xi_w)$ move
  $\theta = 2\,\mathrm{atan2}(n,w)$ by
  $\delta\theta = \sin\theta\,(\xi_n - \xi_w) + O(\xi^2)$:
  $\lvert\delta\theta\rvert/\theta \le \lvert\xi_n - \xi_w\rvert$, no amplification
  anywhere for $n^2$ and $w$ normal (*Range and scaling* below), and $\delta\theta \to 0$ at $\pi$.
  ($2\arccos w$ fails at $0$ and $2\arcsin n$ at $\pi$; $\mathrm{atan2}$ uses both.)

*Proof.* (a) Differentiate; $1 - \cos\theta = \theta^2/2 - \theta^4/24 + \cdots$ inverts to
$\theta = \sqrt{2\delta}(1 + \delta/12 + \cdots)$. (b) (a) with
$\lvert\delta\chi\rvert = \mathsf u/2$; $\sin\theta \approx \theta$ near $0$, $\approx \pi - \theta$
near $\pi$; $\cos\theta$ rounds to $1$ iff $\theta^2/2 < \mathsf u/2$. (c)
$\partial_n\mathrm{atan2} = w/(n^2 + w^2)$, $\partial_w\mathrm{atan2} = -n/(n^2 + w^2)$, and
$2nw = \sin\theta$ for unit $q$. $\square$

The information about a small angle in a *matrix* is in $R - R^\top = 2\sin\theta\,[\hat n]_\times$
(relative precision $\mathsf u$), not in the trace (absolute precision $\mathsf u$); `from_matrix`
(SO.12) reads the former.

**Range and scaling.** SO.6(c) takes $n$ and $w$ with relative errors of order $\mathsf u$, which
needs $n^2 = x^2 + y^2 + z^2$ to be a normal number: below $\lvert u\rvert \approx 10^{-154}$ (`f64`)
or $10^{-19}$ (`f32`) the squares underflow, $n$ is inexact or $0$, and (c) says nothing. There
$L = 2u/w\,(1 + O(n^2/w^2))$ for $w > 0$ (from $\mathrm{atan2}(n,w)/n \to 1/w$), which does not use $n$
to first order: that is the series arm of $r$ (`NUMERICS.md` §4), whose accuracy this page does not
derive, and no `Log` check goes below $\theta = 10^{-15}$. `Exp` errs the other way: $\theta$ carries
an absolute error $\approx\mathsf u\theta$, so the rotation is off by $\approx\theta\mathsf u$ (SO.4,
`Checked:`), the effect of a relative $\mathsf u$ change of $\varphi$: conditioning linear in
$\theta$, not a defect, and negligible for the $\theta \le \pi$ that `Log` returns. $\theta^2$
overflows for $\theta \gtrsim 1.3\times10^{154}$ (`f64`).

**Checked (SO.5):** mpmath 1.3.0, 80 digits. (a): 86 unit $q$
($\theta\in\{10^{-10}, 10^{-8}, 10^{-3}, \pi - 10^{-3}, \pi - 10^{-6}, \pi - 10^{-12}\}$ and 80
uniform on $(0.001, \pi)$): $L$ against the principal matrix logarithm computed from its definition
by eigendecomposition ($V\,\mathrm{diag}(\log\lambda_i)V^{-1}$), $3.0\times10^{-69}$
(ill-conditioned eigenvectors at $10^{-10}$). 80 further $q$ ($\theta$ from $10^{-15}$ to
$\pi - 10^{-6}$) and their negatives: $\mathrm{Exp}\,L = R$ and $\lVert L\rVert \le \pi$,
$5.3\times10^{-81}$; (b) $8.4\times10^{-81}$; (c) with $\lambda\in[0.3,3]$, $4.2\times10^{-81}$.
`mp.logm` agrees for $\theta \le 3$ ($3\times10^{-81}$) but returns a complex, non-principal result
from an axis-dependent onset in $(3.02, 3.03]$ (30 random axes, grid $0.005$ on $[2.5, 3.1]$, 30 and 80
digits; $3.026$ for the axis $(1,2,3)$, $3.025$ for $(0,0,1)$), which persists up to $3.1$, identically
in mpmath 1.3.0 and 1.4.1 (see the maths index, open items).
**Checked (SO.6):** the sensitivity of (c), second-order residual $\le 0.93\,\xi^2$ at
$\xi = 10^{-30}$; (a), (b) in `f64`. Recipe: $\theta$ swept over 1500 equispaced values in
$[\theta_0, 2\theta_0]$ (in $[\pi - 2\delta, \pi - \delta]$ for the last two columns, $\delta = 10^{-4}$,
$10^{-6}$), one random axis each; $q$ the correctly rounded $(\cos\frac\theta2, \sin\frac\theta2\,\hat n)$,
so a double that is not exactly unit; the true $\theta$ is $2\,\mathrm{atan2}(n,w)$ of that double $q$
at 80 digits; the entry is the maximum over the sweep of $\lvert\theta_{\mathsf f} - \theta\rvert/\theta$
in units of $\mathsf u$. ($\chi$ depends on $\theta$ only, so an arccos row does not average over
axes: it is one rounding realization per $\theta$, and the sweep is what makes it a maximum.)

| $\theta$ from | $10^{-2}$ | $10^{-4}$ | $10^{-6}$ | $10^{-8}$ | $\pi - 10^{-4}$ | $\pi - 10^{-6}$ |
|---|---|---|---|---|---|---|
| $\arccos$ of $\chi$ correctly rounded | $4.8\times10^{3}$ | $4.8\times10^{7}$ | $4.8\times10^{11}$ | $9.0\times10^{15}$ (all lost) | $1.55\times10^{3}$ | $1.54\times10^{5}$ |
| model $1/(2\theta_0^2)$, $1/(2\pi\delta)$ | $5\times10^{3}$ | $5\times10^{7}$ | $5\times10^{11}$ | $5\times10^{15}$ | $1.6\times10^{3}$ | $1.6\times10^{5}$ |
| $\arccos$ of $(\operatorname{tr}\mathrm{fl}(R) - 1)/2$, entries of $R(q/\lVert q\rVert)$ rounded, `f64` trace | $1.7\times10^{4}$ | $1.7\times10^{8}$ | $1.7\times10^{12}$ | $9.0\times10^{15}$ | $3.3\times10^{3}$ | $3.6\times10^{5}$ |
| same, $R_{ii} = 1 - 2(\cdot^2 + \cdot^2)$ evaluated in `f64` from the double $q$ | $1.7\times10^{4}$ | $1.7\times10^{8}$ | $1.7\times10^{12}$ | $9.0\times10^{15}$ | $2.1\times10^{4}$ | $2.4\times10^{6}$ |
| $2\,\mathrm{atan2}(n, w)$, $n = \sqrt{x^2 + y^2 + z^2}$ in `f64` | 2.2 | 2.1 | 2.2 | 1.7 | 0.64 | 0.64 |

The model is the worst end of the sweep, and first order: at $10^{-8}$, $\chi$ rounds to $1$ and
$\theta$ is lost (relative error $1 = 9.0\times10^{15}\,\mathsf u$). Near $\pi$ the last row is the rounding of $\theta$ itself:
the sensitivity $\sin\theta \approx \pi - \theta$ has removed the input errors. Going through $R$ costs a
further factor $\approx3.5$ at small angles and, near $\pi$, $2.1$ to $2.3$ (rounded entries) or up
to $16$ (the diagonal formula). Scripts not committed. **Permanent:** planned, corpus `so3_log`
(strata `theta:1e-k`, `theta:pi-1e-k`, `q:w0`, `q:nonunit`), proptests `so3_log_w_flip`,
`so3_log_at_w0_is_a_function_of_the_sign`, `exp_log_roundtrip_*`; none for SO.6 beyond those strata.

## 4. Jacobians, their inverse, $\mathrm{Ad}$, $\mathrm{ad}$

**Proposition SO.7.** With $a = (1-\cos\theta)/\theta^2$, $b = (\theta - \sin\theta)/\theta^3$
($\tfrac12$, $\tfrac16$ at $0$):
$$
J_r(\varphi) = \sum_{m\ge0}\frac{(-W)^m}{(m+1)!} = I - aW + bW^2,\qquad J_l(\varphi) = J_r(-\varphi) = I + aW + bW^2 .
$$
Moreover $J_l = J_r^\top = R\,J_r$, $\det J = 2a$, and in axis form
$J_r = \frac{\sin\theta}\theta I + (1 - \frac{\sin\theta}\theta)\hat n\hat n^\top - \frac{1-\cos\theta}\theta[\hat n]_\times$.

*Proof.* By SO.3 the odd terms are $-W^{2j+1}/(2j+2)! = -(-1)^j\theta^{2j}W/(2j+2)!$, summing to
$-W(1 - \cos\theta)/\theta^2 = -aW$; the even terms
$W^{2j+2}/(2j+3)! = (-1)^j\theta^{2j}W^2/(2j+3)!$ sum to $W^2(\theta - \sin\theta)/\theta^3 = bW^2$.
$J_l$: $W \to -W$. $J_l = J_r^\top$ since $W^\top = -W$ and $(W^2)^\top = W^2$; $J_l = RJ_r$ is
LG.9(b) with SO.9; $\det J = 2a$ is LG.9(c) at $N = 0$. The axis form uses
$W^2 = \theta^2(\hat n\hat n^\top - I)$, $b\theta^2 = 1 - \sin\theta/\theta$. $\square$

**Proposition SO.8 ($J^{-1}$ by algebra in $\mathcal A$).** For $\theta\notin 2\pi\mathbb Z_{>0}$,
$$
J_r^{-1} = I + \tfrac12W + cW^2,\quad J_l^{-1} = I - \tfrac12W + cW^2,\quad
c = \frac{1}{\theta^2} - \frac{1+\cos\theta}{2\theta\sin\theta} = \frac{1}{\theta^2} - \frac{\cot(\theta/2)}{2\theta}
= \sum_{j\ge1}(-1)^{j-1}\frac{B_{2j}}{(2j)!}\theta^{2j-2},
$$
the series for $\theta < 2\pi$ ($B_{2j}$ Bernoulli:
$c = \frac1{12} + \frac{\theta^2}{720} + \frac{\theta^4}{30240} + \frac{\theta^6}{1209600} + \cdots$).
$c\to\frac1{12}$ at $0$, $c = \pi^{-2}$ at $\pi$, and $c \approx 1/(2\pi(2\pi - \theta))$ as
$\theta \to 2\pi$.

*Proof.* Take $\hat J = I + \alpha W + cW^2$ with unknowns $\alpha, c$ and
$\beta = \sin\theta/\theta = 1 - b\theta^2$. By (2) with $x = (1, -a, b)$, $y = (1, \alpha, c)$,
$J_r\hat J = I$ is
$$
\beta\alpha + a\theta^2 c = a,\qquad -a\alpha + \beta c = -b .
$$
Its determinant is
$\beta^2 + a^2\theta^2 = \big(\sin^2\theta + (1-\cos\theta)^2\big)/\theta^2 = 2a$, nonzero exactly
off $2\pi\mathbb Z_{>0}$ (no division by $\beta$, which vanishes at $\pi$). Cramer:
$\alpha = (a\beta + a\theta^2b)/2a = (\beta + b\theta^2)/2 = \tfrac12$ and $c = (a^2 - b\beta)/2a$,
with $a^2 - b\beta = (2 - 2\cos\theta - \theta\sin\theta)/\theta^4$, so
$c = \theta^{-2} - \sin\theta/(2\theta(1-\cos\theta))$; and
$\sin\theta/(1-\cos\theta) = (1+\cos\theta)/\sin\theta = \cot\frac\theta2$. Then $J_r\hat J = I$
gives $\hat J = J_r^{-1}$ (square). $J_l^{-1}$: $\varphi \to -\varphi$ ($c$ is even in $\theta$).
*Series.* $g(z) = (e^z - 1)/z$ has $1/g(z) = z/(e^z - 1) = \sum B_mz^m/m!$ ($B_1 = -\tfrac12$),
holomorphic for $\lvert z\rvert < 2\pi$; $W$ is normal with spectrum $\{0, \pm i\theta\}$, so for
$\theta < 2\pi$, $J_r^{-1} = (1/g)(-W) = I + \tfrac12W + \sum_{j\ge1}B_{2j}W^{2j}/(2j)!$ and
$W^{2j} = (-\theta^2)^{j-1}W^2$ gives the series, whose radius $2\pi$ is the domain of `NUMERICS.md`
§12. $\square$

*Consequences.* Terms of $c = \theta^{-2} - \cot(\theta/2)/2\theta$ are $\approx\theta^{-2}$ for a
value $\to\frac1{12}$: amplification $12/\theta^2$ (as $6/\theta^2$ for $b$ and $2/\theta^2$ for the
naive $a$, `NUMERICS.md` §4); the $\cot$ form is finite at $\pi$ where $(1+\cos\theta)/\sin\theta$
is $0/0$. On the plane $\perp\varphi$, $W$ acts as multiplication by $i\theta$ and $J_r^{-1}$ as
$\frac\theta2\cot\frac\theta2 + i\frac\theta2$, of modulus
$\theta/(2\sin\frac\theta2) = \lVert J^{-1}\rVert_2$ (LG.16(a)): the blow-up of $c$ at $2\pi$ is the
conditioning of $J$, not of the formula. $J_r^{-1} - J_l^{-1} = W$ exactly. In axis form
$J_r^{-1} = \frac\theta2\cot\frac\theta2\,I + (1 - \frac\theta2\cot\frac\theta2)\hat n\hat n^\top + \frac\theta2[\hat n]_\times$.

**Proposition SO.9.** $\mathrm{Ad}_R = R$ and $\mathrm{ad}_\varphi = W$ on tangent coordinates.
*Proof.* LG.4(c) at $N = 0$; directly, $R[\varphi]_\times R^\top = [R\varphi]_\times$ since
$R(a\times b) = Ra\times Rb$ for $\det R = 1$, so
$R\,\mathrm{Exp}(\varphi)R^\top = \mathrm{Exp}(R\varphi)$; and
$[W, \psi^\wedge] = (\varphi\times\psi)^\wedge$. $\square$

**Checked (SO.7–SO.9):** mpmath 1.3.0, 80 digits, 70 tangents
($\theta\in\{10^{-10}, 10^{-3}, 1, \pi-10^{-6}, \pi, 5, 2\pi - 10^{-3}, 7, 20, 30\}$ and 60 uniform
on $[0.01, 6.2]$). $J_{r,l}$: the series of the definition, summed to $10^{-90}$, against the closed
forms $3.2\times10^{-71}$ (relative; at $\theta = 30$ the terms reach $10^{13}$), the axis form
$2.7\times10^{-72}$, $J_l - J_r^\top$ $1.1\times10^{-81}$, $J_l - RJ_r$ $3.4\times10^{-81}$,
$\det J - 2a$ $3.2\times10^{-81}$. $J^{-1}$ (excluding $\lvert\theta - 2\pi\rvert < 10^{-2}$):
`mp.inverse` of $J_r$, $J_l$ against the formulas of SO.8 with $c$ from the closed form (80 guard
digits), $1.1\times10^{-79}$; the Cramer system ($\beta^2 + a^2\theta^2 = 2a$, $\alpha = \frac12$,
$c = (a^2 - b\beta)/2a$) $\le 1.2\times10^{-79}$ (also exact in sympy); the $\cot$ and
$(1+\cos)/\sin$ forms agree $3.8\times10^{-76}$; the Bernoulli series ($\theta < 2\pi - 1$,
relative) $1.6\times10^{-80}$; sympy: the four leading terms of $k, a, b, c, d, e, r$ of the
`NUMERICS.md` §4 table all match their Taylor expansions; the axis form of $J_r^{-1}$,
$5.3\times10^{-81}$; the two-term form of $c$ at $\theta = 10^{-10}$ shows the predicted
cancellation ($2\times10^{-61}$ at 80 digits). $\mathrm{Ad}_R$, $R\,\mathrm{Exp}\,R^\top$,
$\mathrm{ad}_\varphi$ against LG.3 definitions: $\le 3.4\times10^{-80}$. Scripts not committed.
**Permanent:** planned, corpus `so3_jr`, `so3_jl`, `so3_jr_inv`, `so3_jl_inv`, `coeff_a`, `coeff_b`,
`coeff_c`; proptests `jl_is_ad_jr_*`, `jacobians_match_dual_*`, `adjoint_identity_*`,
`branch_continuity_*`.

## 5. Action Jacobians

**Proposition SO.10 (`NUMERICS.md` §2.4).** For $F(X, p) = Xp = Rp$ with $X \in \mathrm{SO}(3)$ of
matrix $R$, and the derivatives $\mathrm D^{\mathrm{right}}_X$, $\mathrm D^{\mathrm{left}}_X$ of
`lie-groups.md` Definition LG.11 (vector output, compared by subtraction):
$$
\mathrm D^{\mathrm{right}}_X F = -R[p]_\times,\qquad \mathrm D^{\mathrm{left}}_X F = -[Rp]_\times,\qquad
\partial F/\partial p = R\ \text{(both sides)}.
$$
For $X = (R, t) \in \mathrm{SE}(3)$ and $F = Rp + t$ (tangent $[\delta\varphi;\delta\rho]$):
$\mathrm D^{\mathrm{right}}_X F = [\,-R[p]_\times,\ R\,]$ and
$\mathrm D^{\mathrm{left}}_X F = [\,-[Rp + t]_\times,\ I\,]$.

*Proof.* $\mathrm{Exp}(\delta) = I + [\delta]_\times + O(\lVert\delta\rVert^2)$. Right:
$R\,\mathrm{Exp}(\delta)p = Rp + R[\delta]_\times p + O = Rp - R[p]_\times\delta + O$ (as
$[\delta]_\times p = -[p]_\times\delta$). Left:
$\mathrm{Exp}(\delta)Rp = Rp + \delta\times Rp + O = Rp - [Rp]_\times\delta + O$. In $p$:
$R(p + \delta) = Rp + R\delta$. The two sides agree through
$X\,\mathrm{Exp}(\delta) = \mathrm{Exp}(R\delta)X$, i.e.
$\mathrm D^{\mathrm{left}} = \mathrm D^{\mathrm{right}}R^{-1}$, and
$R[p]_\times R^\top = [Rp]_\times$ (SO.9). SE(3): $\mathrm{Exp}(\delta)$ maps
$p \mapsto p + \delta\varphi\times p + \delta\rho + O$ (LG.2(a)), so the right perturbation gives
$R(p + \delta\varphi\times p + \delta\rho) + t$ and the left one
$Xp + \delta\varphi\times(Rp + t) + \delta\rho$. $\square$

**Checked:** mpmath 1.3.0, 100 digits, 30 random $(R, p)$ ($\theta \le 3$, $p\sim U[-3,3]^3$, $t$
uniform): central differences $[F(he_j) - F(-he_j)]/2h$, $h = 10^{-30}$, of the operation itself
with `mp.expm` for $\mathrm{Exp}$ (4x4 hat matrix for SE(3)); max $6.1\times10^{-61}$ (SO(3) both
sides, and SE(3) right), $1.0\times10^{-60}$ (SE(3) left), $1.3\times10^{-71}$
($\partial/\partial p$), all limited by the differencing. Script not committed. **Permanent:**
planned, proptests `jacobians_match_dual_*` (`act_jacobians`).

## 6. `from_matrix`

**Lemma SO.11 (the matrix $K$).** For $M\in\mathbb R^{3\times3}$ let $K(M)$ be the symmetric
$4\times4$ matrix, indices $(w,x,y,z)$,
$$
\begin{aligned}
K_{ww} &= M_{00} + M_{11} + M_{22}, & K_{wx} &= M_{21} - M_{12}, & K_{xy} &= M_{10} + M_{01},\\
K_{xx} &= M_{00} - M_{11} - M_{22}, & K_{wy} &= M_{02} - M_{20}, & K_{xz} &= M_{02} + M_{20},\\
K_{yy} &= -M_{00} + M_{11} - M_{22}, & K_{wz} &= M_{10} - M_{01}, & K_{yz} &= M_{21} + M_{12},\quad K_{zz} = -M_{00} - M_{11} + M_{22}.
\end{aligned}
$$
(a) $\operatorname{tr}(R_\square(q)^\top M) = q^\top K(M)\,q$ for all $q \in \mathbb R^4$. (b)
$K(R_\square(q)) = 4qq^\top - \lVert q\rVert^2I_4$; for unit $q$, $K(R(q)) + I_4 = 4qq^\top$. (c)
With $m = (m_w, m_x, m_y, m_z) = (\operatorname{tr}M, M_{00}, M_{11}, M_{22})$ and
$P_i = 1 + K_{ii}$: $P_i = 1 - m_w + 2m_i$, hence $P_i - P_j = 2(m_i - m_j)$ and $\sum_iP_i = 4$
**for every $M$**.

*Proof.* From (1), entrywise: $R_{00} = w^2 + x^2 - y^2 - z^2$ (cyclic for $R_{11}, R_{22}$),
$R_{21} = 2(yz + wx)$, $R_{12} = 2(yz - wx)$, $R_{02} = 2(xz + wy)$, $R_{20} = 2(xz - wy)$,
$R_{10} = 2(xy + wz)$, $R_{01} = 2(xy - wz)$. (a): $\sum R_{ij}M_{ij}$ collects to the quadratic
form with the stated coefficients ($w^2$: $M_{00}+M_{11}+M_{22}$; $wx$: $2(M_{21} - M_{12})$; $xy$:
$2(M_{10} + M_{01})$; …). (b): substitute the entries, e.g.
$K_{ww} = 3w^2 - n^2 = 4w^2 - \lVert q\rVert^2$, $K_{xx} = 4x^2 - \lVert q\rVert^2$, $K_{wx} = 4wx$,
$K_{xy} = 4xy$. (c): $K_{ii} = 2m_i - m_w$ for all four $i$ ($i = w$: $K_{ww} = m_w$); the diagonal
sums to $0$. $\square$

**Proposition SO.12 (Shepperd).** Let $R = R(\hat q)$ and $P_i = 4\hat q_i^2$ (SO.11(b)).

- (a) The four squares are $4w^2 = 1 + \operatorname{tr}R$, $4x^2 = 1 + R_{00} - R_{11} - R_{22}$,
  $4y^2 = 1 - R_{00} + R_{11} - R_{22}$, $4z^2 = 1 - R_{00} - R_{11} + R_{22}$; the products are
  $4wx = R_{21} - R_{12}$, $4wy = R_{02} - R_{20}$, $4wz = R_{10} - R_{01}$,
  $4xy = R_{10} + R_{01}$, $4xz = R_{02} + R_{20}$, $4yz = R_{21} + R_{12}$.
- (b) Column $j$ of $K + I$ is $4\hat q_j\hat q$. So for a pivot $j$ with $P_j > 0$, the extraction
  $q_j = +\frac12\sqrt{P_j}$, $q_i = (K + I)_{ij}/(4q_j)$ returns $+\hat q$ if $\hat q_j > 0$ and
  $-\hat q$ otherwise: **the pivot component is positive whatever the sign of $\hat q$**, so a pivot
  $x$, $y$ or $z$ can return $w < 0$ (*Sign*, below). Pivot $w$:
  $(\frac12\sqrt{P_w},\ \frac{R_{21} - R_{12}}{4w},\ \frac{R_{02} - R_{20}}{4w},\ \frac{R_{10} - R_{01}}{4w})$;
  pivot $x$: $w = \frac{R_{21} - R_{12}}{4x}$, $y = \frac{R_{10} + R_{01}}{4x}$,
  $z = \frac{R_{02} + R_{20}}{4x}$; pivot $y$: $w = \frac{R_{02} - R_{20}}{4y}$,
  $x = \frac{R_{10} + R_{01}}{4y}$, $z = \frac{R_{21} + R_{12}}{4y}$; pivot $z$:
  $w = \frac{R_{10} - R_{01}}{4z}$, $x = \frac{R_{02} + R_{20}}{4z}$,
  $y = \frac{R_{21} + R_{12}}{4z}$.
- (c) $\sum P_i = 4$ gives $\max_iP_i \ge 1$: the largest pivot has
  $\lvert\hat q_j\rvert \ge \frac12$. By SO.11(c) the largest $P$ is the largest of
  $m = (\operatorname{tr}R, R_{00}, R_{11}, R_{22})$: the nested `branch` of `NUMERICS.md` §3.4
  compares exactly these.
- (d) Let $M = R + E$, $\max_{ij}\lvert E_{ij}\rvert \le \varepsilon$, pivot $j$ the largest, and
  the sign of $\hat q$ chosen so that $\hat q_j > 0$. To first order in $\varepsilon$,
  $\lvert\delta\hat q_j\rvert \le \frac34\varepsilon$ and $\lvert\delta\hat q_i\rvert \le \frac74\varepsilon$,
  attained at $\hat q = (\frac12,\frac12,\frac12,\frac12)$; and the extracted $q = \hat q + \delta\hat q$ has
  $\lvert\eta\rvert = \lvert 2\hat q\cdot\delta\hat q\rvert \le 6\varepsilon$ (the tie attains $4.5\varepsilon$).
  With the trace pivot instead, $\delta\hat q_i \approx -\hat q_i\,\delta P_w/(8w^2)$, of size
  $\le\frac38\varepsilon/w^2$: unbounded as $\theta\to\pi$.

*Proof.* (a), (b): SO.11(b); (c): SO.11(c). (d) $P_j = 1 \pm M_{aa} \pm M_{bb} \pm M_{cc}$, so
$\lvert\delta P_j\rvert \le 3\varepsilon$ and $\delta\hat q_j = \delta P_j/(8\hat q_j)$,
$\lvert\delta\hat q_j\rvert \le 3\varepsilon/(8\cdot\frac12)$. For $i \ne j$,
$\hat q_i = N_i/(4\hat q_j)$ with $N_i$ a sum of two entries,
$\lvert\delta N_i\rvert \le 2\varepsilon$,
$\delta\hat q_i = \delta N_i/(4\hat q_j) - \hat q_i\,\delta\hat q_j/\hat q_j$ and
$\lvert\hat q_i\rvert \le \hat q_j$:
$\lvert\delta\hat q_i\rvert \le \varepsilon + \frac34\varepsilon$. Ties ($P_i = 1$ for all $i$) make
$\hat q_j = \frac12$ and every inequality an equality for suitable signs of $E$. For $\eta$:
$\lvert\eta\rvert \le 2\big(\frac34\hat q_j + \frac74\sum_{i\ne j}\lvert\hat q_i\rvert\big)\varepsilon$,
$\sum_{i\ne j}\lvert\hat q_i\rvert \le \sqrt{3(1 - \hat q_j^2)}$ (Cauchy–Schwarz), and
$\frac34\hat q_j + \frac74\sqrt{3(1 - \hat q_j^2)}$ decreases on $[\frac12, 1]$ from $3$. $\square$

*Consequences.* The radicand of the largest pivot is a sum of four terms of size $\le 1$ with value
$\ge 1$ (no cancellation beyond a factor $4$), and it is divided by a quantity $\ge\frac12$. A small
$P$ (the trace pivot at $\theta\to\pi$, $P_w = 4\cos^2\frac\theta2$) is a sum of $O(1)$ terms of
size $\to 0$: all digits lost, then divided by the small $w$. Rounding errors of the arithmetic
enter as perturbations of the same kind (numerators and radicand), so the conclusion of (d) carries
over with $\varepsilon$ a small multiple of $\mathsf u$; this is argued, not proved. Near $\pi$ the
axis is read from the symmetric part (the largest of $x, y, z$ has square $\ge\frac13(1 - w^2)$),
but its **sign** relative to $w$ is that of $R - R^\top = 2\sin\theta[\hat n]_\times$, of size
$\approx\pi - \theta$: for $\pi - \theta \lesssim \varepsilon$ it is noise, and `Log` (which flips
on the sign of $w$) returns $\pm\theta\hat n$ (SO.5(d)). This is why matrix-input strata report
backward error only (`NUMERICS.md` §11).

*Sign.* Pivot $w$ returns $w > 0$; pivots $x, y, z$ return their own component positive and $w$ of
either sign; `NUMERICS.md` §3.4 states no sign rule (no $w \ge 0$ canonicalization). So the sign of
the output is a **discontinuous function of $R$**: at a tie between two pivots one returns $\hat q$
and the other $-\hat q$ (a jump of $2$ in $\mathbb R^4$), for an exact rotation and for a
perturbation of any size. Example: $\hat q = (0.6, -0.6, 0.3, \sqrt{0.19})$ has
$\operatorname{tr}R = R_{00} = 0.44 > R_{22} = 0.1 > R_{11} = -0.1$; adding $+10^{-9}$ to $R_{11}$
selects $w$ and returns $(0.6, -0.6, 0.3, 0.4359)$, adding $-10^{-9}$ selects $x$ and returns exactly
its negative. `Log` is blind to the sign except at $w = 0$ (SO.5(b), (d)); a consumer that compares,
averages or interpolates successive outputs componentwise must align signs first. Canonicalizing to
$w \ge 0$ would move the jump to $w = 0$ ($\theta = \pi$), where `Log` already jumps: a §3.4
decision (open item, maths index).

**Proposition SO.13 (non-orthogonal input, nearest rotation).** Let $M \in \mathbb R^{3\times3}$ be
arbitrary and finite.

- (a) *Total, closed form.* By SO.11(c), $\max_iP_i \ge 1$ for every $M$: in exact arithmetic the
  selected arm has radicand $\ge1$, pivot $\ge\frac12$, no iteration and no division by zero
  ($M = 0$ gives $(\frac12, 0, 0, 0)$). Nothing loops on degenerate input
  (`from_matrix_never_iterates`). In floating point the result is finite while the sums of at most
  three entries of $M$ do not overflow ($M = 10^{308}I$: $\operatorname{tr}M = \infty$ in `f64`). The
  three unselected arms are not safe: a lane mask
  ([`0003`](../decisions/0003-the-scalar-that-cannot-say-less-than.md)) evaluates all four, and for a
  non-rotation $P_i < 0$ is possible ($M = -I$: $P_w = -2$, $P_x = P_y = P_z = 2$), so each arm's
  radicand and divisor need the safe argument of `0003` (item 3), whose form is the code's to choose.
- (b) *Not a projection.* It reads one column of $K(M) + I$ as if that matrix were
  $4\hat q\hat q^\top$ (rank one). On $\mathrm{SO}(3)$ the four columns agree up to sign; off it they
  do not, so the result depends on the pivot. **Modulo the global sign**, two pivots differ by
  $O(\lVert M - R\rVert)$; but each returns its own component positive (SO.12(b)), so at a pivot tie
  they can return $q$ and $-q$, a jump of $2$, for a perturbation of any size (SO.12, *Sign*).
  $\lVert\hat q\rVert^2 \ne 1$ in general. To first order the error is SO.12(d) for any $M$ within
  $\varepsilon$ of $R$.
- (c) *Nearest rotation.*
  $\min_{R\in\mathrm{SO}(3)}\lVert M - R\rVert_F^2 = \lVert M\rVert_F^2 + 3 - 2\lambda_{\max}(K(M))$,
  attained at $R(q^\star)$ for the unit top eigenvector $q^\star$ of $K(M)$ [Horn; Bar-Itzhack].
  With the SVD $M = U\Sigma V^\top$, $\sigma_1 \ge \sigma_2 \ge \sigma_3 \ge 0$,
  $d = \det(UV^\top)$: $R_\star = U\,\mathrm{diag}(1,1,d)\,V^\top$,
  $\lambda_{\max} = \sigma_1 + \sigma_2 + d\sigma_3$, and the minimizer is unique iff
  $\sigma_2 + d\sigma_3 > 0$. With `svd3`'s signed convention (`PHASE2.md` §6:
  $\det U = \det V = 1$, $\sigma_3$ signed) this is $R_\star = UV^\top$ **provided
  $\lvert\sigma_3\rvert \le \min(\sigma_1, \sigma_2)$**, i.e. the signed value is the smallest in
  magnitude, which `PHASE2.md` §6 does not state (open item, maths index). Without it,
  $\Sigma = \mathrm{diag}(1,1,-5)$, $U = V = I$ gives $UV^\top = I$ at squared distance $36$ from
  $M$, against $20$ for $R_\star$.

*Proof.* (a) SO.11(c). (b) SO.11(b), and SO.12(d). (c)
$\lVert M - R\rVert_F^2 = \lVert M\rVert_F^2 + 3 - 2\operatorname{tr}(R^\top M)$ and
$\operatorname{tr}(R(q)^\top M) = q^\top K(M)q$ (SO.11(a)) is maximized over unit $q$ by the top
eigenvector. For the SVD form:
$\operatorname{tr}(R^\top M) = \operatorname{tr}(V^\top R^\top U\,\Sigma) = \operatorname{tr}(Q'D)$
with $D = \mathrm{diag}(\sigma_1,\sigma_2,d\sigma_3)$ and
$Q' = V^\top R^\top U\,\mathrm{diag}(1,1,d) \in \mathrm{SO}(3)$ (the determinant is $d\cdot d = 1$),
and $Q'$ ranges over $\mathrm{SO}(3)$ with $R$. For diagonal $D = \mathrm{diag}(s_i)$, $K(D)$ is
*diagonal* ($M_{ij} \pm M_{ji} = 0$) with entries $s_1{+}s_2{+}s_3$, $s_1{-}s_2{-}s_3$,
$-s_1{+}s_2{-}s_3$, $-s_1{-}s_2{+}s_3$, and $\max_q q^\top K q$ is the largest of them. With
$s = (\sigma_1,\sigma_2,d\sigma_3)$ the first exceeds the others by $2(\sigma_2 + d\sigma_3)$,
$2(\sigma_1 + d\sigma_3)$, $2(\sigma_1 + \sigma_2)$, all $\ge 0$; it is attained at $q = (1,0,0,0)$,
i.e. $Q' = I$, i.e. $R^\top = V\,\mathrm{diag}(1,1,d)\,U^\top$; strictness gives uniqueness.
$\square$

If `from_matrix` normalizes with the step of SO.14 (`NUMERICS.md` §3.4 says "then normalize" without
saying how: open item, maths index), the step assumes $\lvert\eta\rvert \lesssim 2^{-26}$
(SO.15(a)). By SO.12(d) $\lvert\eta\rvert \le 6\varepsilon$ for $M$ within $\varepsilon$ of
$\mathrm{SO}(3)$ in max-entry norm, so that holds for $\varepsilon \lesssim 2\times10^{-9}$
(measured maximum $4.5\varepsilon$: $3\times10^{-9}$). Beyond that, project first (`NUMERICS.md`
§3.4: `svd3` then this): the step is not merely inaccurate but destructive once $\eta \ge 2$
(SO.14): $M = \frac{11}{3}I$ extracts $(\sqrt3, 0, 0, 0)$, $\eta = 2$, and the step returns $0$.

**Checked (SO.11–SO.13):** mpmath 1.3.0, 80 digits. SO.11: 50 random non-unit $q$ and random $M$:
(a) $3.4\times10^{-80}$; (b), (c) exactly $0$ (inputs short enough for exact arithmetic at 269
bits); $P_{\max} \ge 1$ on $2\times10^3$ random $M$ of scale $10^{-6}$ to $10^{6}$ (least
$\max_iP_i$: $1.000000124$). SO.12: 30 random rotations per dominant component, all four extractions
with $P_j > 10^{-30}$ against $\hat q$, signs aligned: $\le 2.7\times10^{-76}$; the pivot chosen by the
largest of $m$ is the dominant component in all 120. Axis $(1,2,3)/\sqrt{14}$, $\theta \in \{10^{-12},
\pi - 10^{-10}, \pi\}$: $\le 1.7\times10^{-81}$ for $P_j \ge 0.28$, and $\approx 10^{-81}/P_j$
below, to a factor $3$ (precision-limited, digits lost as $P_j \to 0$): $3.1\times10^{-62}$ at $P_w = 10^{-20}$
($\theta = \pi - 10^{-10}$), $8.3\times10^{-57}$ at $P_x = 7\times10^{-26}$ ($\theta = 10^{-12}$). (d): the
tie $(\frac12,\frac12,\frac12,\frac12)$ and 25 random unit quaternions, **all $2^9$ sign patterns**
$E_{ij} = \pm\varepsilon$, $\varepsilon = 10^{-30}$, largest-$m$ pivot, extracted $q$ sign-aligned with
$\hat q$: $\sup\lvert\delta\hat q_i\rvert/\varepsilon = 1.75$ at the tie (the bound; the pivot component
alone $0.75$), at most $1.44$ on the random ones (a sampled maximum: only $\frac74$ is a bound);
$\sup\lvert\eta\rvert/\varepsilon = 4.5$ at the tie (bound $6$), at most $2.6$ on the random ones.
Trace pivot forced, $\lvert w\rvert = 5\times10^{-9}, 5\times10^{-5}, 5\times10^{-3}$: for the axis
$(1,2,3)/\sqrt{14}$ $1.2\times10^{16}$, $1.2\times10^{8}$, $1.2\times10^4$ (the law
$3\max_i\lvert\hat q_i\rvert/(8w^2)$ to $2$ digits; the axis matters: $8.7\times10^{15}$,
$8.7\times10^{7}$, $8.8\times10^{3}$ for $(1,1,1)/\sqrt3$). SO.13: the counterexample of SO.12
(*Sign*), the two outputs sum to $6\times10^{-10}$; 6 random rotations, $E_{ij}$ uniform in $[-10^{-3}, 10^{-3}]$, the
four normalized extractions sign-aligned with $\hat q$ lie $2.3\times10^{-4}$ to $3.4\times10^{-3}$
from it for $P_j \ge 0.07$ and $7\times10^{-3}$ to $1.6\times10^{-2}$ for $P_j \le 0.012$ (one
instance has $P_y < 0$); unaligned, an extraction of opposite sign is $2.0$ away in 4 of the 6
instances. $M = -I$: $P = (-2, 2, 2, 2)$; $M = \frac{11}{3}I$: $\eta = 2$. (c) 60 random $M$
($\det<0$ half of them): `mp.svd_r` form against the top eigenvector
of $K$ (`mp.eigsy`) $1.8\times10^{-79}$, $\lambda_{\max}$ against
$\operatorname{tr}(R_\star^\top M)$ and against $\sigma_1 + \sigma_2 + d\sigma_3$
$1.0\times10^{-79}$, and no better rotation among 50 random ones per $M$; $M = \mathrm{diag}(1,1,-5)$:
$\lambda_{\max} = 5$, squared distance $20$, against $36$ for $I$. Scripts not committed.
**Permanent:** planned, corpus `so3_from_matrix` (backward error), proptest
`from_matrix_never_iterates`; the SVD route is `svd3`'s (PHASE2 §6), against `mp.svd_r`.

## 7. Renormalization

**Proposition SO.14 (the Newton step).** For $q \ne 0$ and $\eta = \lVert q\rVert^2 - 1$,
$q' = q\,(3 - \lVert q\rVert^2)/2 = q\,(1 - \eta/2)$ is the first Newton iterate, from
$\lambda_0 = 1$, of $g(\lambda) = \lambda^{-2} - \lVert q\rVert^2$ (root
$\lambda = \lVert q\rVert^{-1}$): $\lambda - g/g' = \lambda(3 - \lVert q\rVert^2\lambda^2)/2$. For all $\eta$
$$
\lVert q'\rVert^2 - 1 = -\tfrac34\eta^2 + \tfrac14\eta^3 \quad(\text{exactly}).
$$
For $\eta < 2$, $q'$ is the positive multiple $1 - \eta/2$ of $q$ (it keeps the direction) and
$\lVert q' - \hat q\rVert = \big\lvert\lVert q'\rVert - 1\big\rvert = \tfrac38\eta^2\,(1 + O(\eta))$.
At $\eta = 2$, $q' = 0$; for $\eta > 2$, $q'$ is a negative multiple of $q$ (reversed): the step is
destructive there, not merely inaccurate. For $\lvert\eta\rvert < 1$ the deviation decreases strictly
($\lvert\eta'\rvert \le \frac34\eta^2 + \frac14\lvert\eta\rvert^3 < \lvert\eta\rvert$); iterating
gives quadratic convergence with constant $-\frac34$. No `sqrt`, no division.

*Proof.* $\lVert q'\rVert^2 = (1 - \eta/2)^2(1 + \eta) = 1 - \frac34\eta^2 + \frac14\eta^3$.
$\lVert q'\rVert - 1 = (\lVert q'\rVert^2 - 1)/(\lVert q'\rVert + 1)$, and $q' - \hat q = (\lVert q'\rVert - 1)\hat q$
when $1 - \eta/2 > 0$. $\square$

**Corollary SO.15 (the bounds of `NUMERICS.md` §3.6).**

- (a) *Reach.* One step leaves a truncation term $\frac34\eta^2$; it is below rounding
  ($\le\mathsf u$) iff $\lvert\eta\rvert \le \sqrt{4\mathsf u/3} = 2^{-26.29}$ (`f64`) or
  $2^{-11.79}$ (`f32`). The `debug_assert!` bounds $2^{-40}$ and $2^{-16}$ are $2^{14}$ and $2^{4}$
  inside it:
a quaternion accepted by `from_wxyz_unchecked` is repaired to rounding level by one step (truncation
$\le\frac34\cdot2^{-80}$ in `f64`).
- (b) *Drift of composition* ("composition never normalizes"). Each component of $pq$ is a sum of 4
  signed products using each $p_i$, $q_j$ once, so
  $\lvert\mathrm{fl}(pq) - pq\rvert_k \le \gamma_4\lVert p\rVert\lVert q\rVert$ (Cauchy–Schwarz;
  $\gamma_4 = 4\mathsf u/(1 - 4\mathsf u)$ [Higham, ch. 3]),
$\lVert\mathrm{fl}(pq)\rVert \le \lVert p\rVert\lVert q\rVert(1 + 2\gamma_4)$, and $\eta$ grows by
at most $\approx16\mathsf u$ per product. The bound $2^{-40} = 2^{13}\mathsf u$ is therefore reached
after $\le 512$ compositions in the worst case (`f32`, $2^{-16} = 2^{8}\mathsf u$: 16); random data
drift like $\sqrt m\,\mathsf u$.
- (c) *Consumers.* Up to the bound, `Log` does not care (SO.5(c)); `act` errs by up to
  $2\lvert\eta\rvert\lVert v\rVert$ and `to_matrix` by $\lvert\eta\rvert\lVert v\rVert$
  ($\le 2^{14}\mathsf u\lVert v\rVert$ at the bound), and they differ from each other by $\eta v$
  (SO.2(d)). Applying $q \mapsto q(3 - \lVert q\rVert^2)/2$ before use removes it.

**Checked:** mpmath 1.3.0, 120 digits: the closed form of SO.14, 200 random $\eta\in[-0.9,1.9]$,
deviation $0$; iterates from $\eta_0 = 10^{-3}$:
$\eta_{k+1}/\eta_k^2 = -0.74975, -0.75000019, -0.75, -0.75$; $\lVert q'\rVert - 1$ at
$\eta = 10^{-6}$: $-3.7499988\times10^{-13}$ against $-\frac38\eta^2$. `f64` (Python floats), one
step from $\lvert\eta\rvert = 2^{-40}, 2^{-30}, 2^{-26}, 2^{-20}$ (2000 random $q$ each):
$\lvert\eta'\rvert$ at most $2.8, 2.6, 4.0, 6.1\times10^3$ $\mathsf u$ (the last is the truncation
term $\frac34\eta^2$); `f32` (numpy) at $2^{-16}, 2^{-12}, 2^{-8}$: $2.8, 3.6, 195$ $\mathsf u$.
(b): $2\times10^4$ random `f64` products, worst relative deviation of $\lVert pq\rVert^2$ from
$\lVert p\rVert^2\lVert q\rVert^2$ $4.6\mathsf u$ (bound $16\mathsf u$); 20 chains of 512 random
unit products: $\max\lvert\eta\rvert = 90\mathsf u$ ($2^{-40} = 8192\mathsf u$). Scripts not
committed. **Permanent:** none specified for the step; corpus `so3_log` stratum `q:nonunit`
($\lVert q\rVert^2 - 1 = \pm2^{-45}$) covers `Log`'s scale invariance.
