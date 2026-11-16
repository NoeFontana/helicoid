# The 3×3 RQ and LQ factorizations

> **Non-normative.** [`NUMERICS.md`](../NUMERICS.md) §22 is the contract; this page derives why it
> holds. Record [`0117`](../decisions/0117-rq-and-lq-are-gram-schmidt-twice-on-scaled-rows.md)
> measures it. Label: `RQ`.

Notation: $M$ with rows $m_1, m_2, m_3$; $M = UR$ with $U$ upper triangular and $R \in \mathrm{SO}(3)$
with rows $r_1, r_2, r_3$; $u$ the unit roundoff, $\gamma_k = ku/(1 - ku)$. $J$ is the reversal
permutation, $J_{ij} = [i + j = 4]$, and $\mathrm{flip}(X) = JXJ$. $\mathrm{low}(X)$ is the strictly
lower triangle of $X$.

## 1. The convention

**RQ.1 (uniqueness).** Where $m_2$ and $m_3$ are linearly independent there is exactly one $(U, R)$
with $M = UR$, $R \in \mathrm{SO}(3)$, $u_{22} > 0$, $u_{33} > 0$:
$r_3 = m_3/\lVert m_3\rVert$, $r_2 = w/\lVert w\rVert$ with $w = m_2 - (m_2\cdot r_3)r_3$,
$r_1 = r_2 \times r_3$, $u_{ij} = m_i\cdot r_j$ ($j \ge i$), and
$u_{11} = \det M/(u_{22}u_{33})$. So $u_{11} > 0 \Leftrightarrow \det M > 0$. $R$ depends on $m_2$,
$m_3$ only, and $(DU, R)$ is the factorization of $DM$ for every positive diagonal $D$.

*Proof.* Row 3 of $M = UR$ is $m_3 = u_{33}r_3$, which fixes $r_3$ by $u_{33} > 0$. Row 2 is
$m_2 = u_{22}r_2 + u_{23}r_3$ with $r_2 \perp r_3$, which fixes $u_{23} = m_2\cdot r_3$ and then
$r_2$ by $u_{22} > 0$. $r_1$ is orthogonal to both, and $\det R = +1$ chooses $r_2 \times r_3$.
Row 1 gives $u_{1j} = m_1\cdot r_j$. $\det M = \det U = u_{11}u_{22}u_{33}$. Left-multiplying by
$D$ scales the rows of $U$ and leaves the construction of $R$ unchanged. $\square$

Where $m_2 \parallel m_3$ (or $m_3 = 0$) every completion of $r_3$ (and $r_2$) to a rotation
factors $M$, with $u_{22} = 0$ (and $u_{33} = 0$). The routine picks one (§22).

**RQ.2 (the flip).** $\mathrm{flip}$ maps upper to lower triangular matrices and $\mathrm{SO}(3)$ to
itself, and $\mathrm{flip}(UR) = \mathrm{flip}(U)\,\mathrm{flip}(R)$. So the LQ factorization
$M = LQ$ in the mirrored convention ($l_{11} > 0$, $l_{22} > 0$, $Q \in \mathrm{SO}(3)$, $l_{33}$
carrying the sign of $\det M$) is $\mathrm{flip}$ of the RQ factorization of $\mathrm{flip}(M)$.
Being a permutation of entries, $\mathrm{flip}$ is exact in floating point.

*Proof.* $J^2 = I$, so $JURJ = (JUJ)(JRJ)$. $(JUJ)_{ij} = U_{4-i,4-j}$ is lower triangular with the
diagonal reversed; $JRJ$ is orthogonal with $\det(JRJ) = \det(J)^2\det R = +1$. Uniqueness is RQ.1
applied to $\mathrm{flip}(M)$, whose rows 2 and 3 are $J m_2$ and $J m_1$. $\square$

## 2. Rounding

**RQ.3 (one projection is not enough).** With $\hat r_3$ unit to $O(u)$ and
$\hat w = \mathrm{fl}(m_2 - \mathrm{fl}(m_2\cdot\hat r_3)\,\hat r_3)$,
$\lvert\hat r_3\cdot\hat w\rvert \le \gamma_4\lVert m_2\rVert$, so the computed $r_2$ makes an angle up
to $\gamma_4/\sin\angle(m_2, m_3)$ with the normal plane of $r_3$. The cross product transmits it to
$R$: $\lVert\hat R^\top\hat R - I\rVert$ grows like $u/\sin\angle(m_2, m_3)$. This is the classical
loss of Gram–Schmidt. With three rows and the cross product it is set by one angle, not by
$\kappa(M)^2$. Measured, it is $6.6\times10^{12}\,u$ at $\angle = 10^{-12}$ (`0117` *Context*).

*Proof.* The dot product carries an absolute error $\gamma_3\lVert m_2\rVert$, and the update one
more rounding per component. $\lVert\hat w\rVert = \lVert m_2\rVert\sin\angle(m_2, m_3)(1 + O(u))$. $\square$

**RQ.4 (twice is enough, with Kahan's test).** Let $w_1 = \hat w$ as in RQ.3 and $w_2$ the same
projection applied to $w_1$. If $\lVert w_2\rVert^2 \ge \tfrac12\lVert w_1\rVert^2$, then
$\lvert\hat r_3\cdot w_2\rvert \le \gamma_4\lVert w_1\rVert \le \sqrt2\,\gamma_4\lVert w_2\rVert$: $r_2$
is orthogonal to $r_3$ to a few $u$ whatever the angle. Otherwise $\lVert w_2\rVert \le
\gamma_4\lVert m_2\rVert$ to first order: $m_2$ lies in the span of $r_3$ to working precision, and
$u_{22} = 0$ with any unit $r_2 \perp r_3$ changes $m_2$ by at most that much (Parlett, *The
Symmetric Eigenvalue Problem*, §6-9, with $\kappa = \sqrt2$).

*Proof.* RQ.3 applied to $w_1$ gives the first bound. $w_2 = w_1 - (w_1\cdot\hat r_3)\hat r_3 + O(u\lVert w_1\rVert)$,
so $\lVert w_2\rVert^2 = \lVert w_1\rVert^2 - (w_1\cdot\hat r_3)^2 + O(u\lVert w_1\rVert^2)$. Failing the
test means $(w_1\cdot\hat r_3)^2 > \tfrac12\lVert w_1\rVert^2$, so $\lVert w_1\rVert < \sqrt2\,\lvert w_1\cdot\hat r_3\rvert
\le \sqrt2\,\gamma_4\lVert m_2\rVert$ by RQ.3, and $\lVert w_2\rVert < \lVert w_1\rVert/\sqrt2$. $\square$

The factor is measured, not only derived (`0117` *Context*). Rows at an angle near $u$ read 6.4 `u`
(`f64`) and 6.7 (`f32`) of orthogonality with $\tfrac12$, 9.7 and 10.8 with $\tfrac14$, 24.6 and 22.4
with $\tfrac1{16}$, and 4853 and 1835 with no test.

**RQ.5 (rowwise backward stability).** If $\lVert\hat R\hat R^\top - I\rVert \le \varepsilon$ and
$\hat u_{ij} = \mathrm{fl}(m_i\cdot\hat r_j)$ for $j > i$ and $j = i = 1$, with $\hat u_{22} =
\lVert w_2\rVert$, $\hat u_{33} = \lVert m_3\rVert$ and the strictly lower triangle zero, then
$\lVert m_i - (\hat U\hat R)_i\rVert \le (c\varepsilon + \gamma_{c'})\lVert m_i\rVert$ for each row,
with small constants $c$ and $c'$. A row is reproduced to a few $u$ of its own norm, so graded rows
lose nothing to the others.

*Proof.* $\sum_j (m_i\cdot\hat r_j)\hat r_j = m_i\hat R^\top\hat R = m_i(I + E)$ with
$\lVert E\rVert \le \varepsilon$. The dropped terms are $m_2\cdot\hat r_1$, $m_3\cdot\hat r_1$ and
$m_3\cdot\hat r_2$, each $O(u\lVert m_i\rVert)$ because $m_2$ and $m_3$ lie in the span of
$\hat r_2$ and $\hat r_3$ to that accuracy (RQ.4). $\lVert w_2\rVert - m_2\cdot\hat r_2 = O(u\lVert m_2\rVert)$. $\square$

## 3. Perturbation

**RQ.6 (the tangent).** Where RQ.1 is unique, $\dot R = \Omega R$ and
$\dot U = \dot M R^\top - U\Omega$, with $X = U^{-1}\dot M R^\top$ and
$\Omega = \mathrm{low}(X) - \mathrm{low}(X)^\top$. $\mathrm{low}(X)$ reads only rows 2 and 3 of
$U^{-1}$, which hold $u_{22}$, $u_{23}$, $u_{33}$. So $R$'s conditioning is that of the bottom
block, and $u_{11}$ (the determinant's sign) does not enter it.

*Proof.* Differentiating $M = UR$ gives $\dot MR^\top = \dot U + U\Omega$ with $\Omega = \dot RR^\top$
skew. Then $U^{-1}\dot MR^\top = U^{-1}\dot U + \Omega$. The first term is upper triangular, so
$\mathrm{low}(X) = \mathrm{low}(\Omega)$, and skewness gives $\Omega$. $\square$

The metric (`NUMERICS.md` §11) weighs forward errors by $\kappa_R = \lVert \dot M \mapsto \dot R\rVert_2
\lVert M\rVert_F$ and $\kappa_U = \lVert \dot M \mapsto \dot U\rVert_2$, the norms of this linear map
at the reference.

## 4. Exactness

**RQ.7 ($\sqrt{x^2}$).** In binary floating point with round-to-nearest, if $x^2$ is normal and finite,
then $\sqrt{\mathrm{fl}(x^2)} = \lvert x\rvert$ after the correctly rounded square root.

*Proof.* Let $\lvert x\rvert = \mu 2^e$ with $\mu \in [1, 2)$. If $\mu = 1$ then $x^2$ is exact. Otherwise
$\mathrm{fl}(x^2) = x^2(1 + \delta)$ with $\lvert\delta\rvert \le u/(1 + u)$. So
$\bigl\lvert\sqrt{\mathrm{fl}(x^2)} - \lvert x\rvert\bigr\rvert < \tfrac u2\lvert x\rvert$, while
both neighbours of $\lvert x\rvert$ lie $2u\,2^e$ away, so the midpoints are $u\,2^e = u\lvert x\rvert/\mu
> \tfrac u2\lvert x\rvert$ away. $\square$

**RQ.8 (a triangle factors exactly).** Let $M$ be upper triangular with a positive diagonal and
every row's largest entry inside the scaling window (RQ.9 extends it), with $m_{22}^2$ and $m_{33}^2$
normal. With entries $+0$ below the diagonal and no $-0$ entry, §22's routine returns $R = I$ and
$U = M$ bit for bit; a diagonal $M$ is the special case. Through $\mathrm{flip}$ (RQ.2), $\mathrm{lq3}$ returns $Q = I$ and $L = M$ for a
lower triangle.

*Proof.* Row 3 is $(0, 0, m_{33})$. Its squared norm is $\mathrm{fl}(m_{33}^2)$, and its root is
$m_{33}$ (RQ.7). The divisions give $0/m_{33} = 0$ and $m_{33}/m_{33} = 1$, so $r_3 = e_3$. The Newton
step leaves it there: $1.5 - 0.5\cdot1 = 1$. Row 2 is $(0, m_{22}, m_{23})$, and $m_2\cdot e_3 = 0 + 0 +
m_{23}$. Both projections give $(0 - 0, m_{22} - 0, m_{23} - m_{23}) = (0, m_{22}, 0)$. That passes
RQ.4's test and gives $u_{22} = m_{22}$ and $r_2 = e_2$. Then $e_2 \times e_3 = e_1$ with $+0$
entries. Each $m_i\cdot r_j$ has one nonzero product, $m_{ij}\cdot1$. $\square$

A $-0$ off the diagonal changes only signs of zeros: below it, $-0/m_{33} = -0$ in $R$; above it,
$(-0)\cdot1 + (+0) = +0$ in $U$. The values stay $I$ and $M$.

**RQ.9 (the row scaling is exact and free).** Scale each row by its own power of two $f_i$ into
`eig3`'s window, factor, and multiply row $i$ of $U$ by $1/f_i$. By RQ.1 the result is the
factorization of $M$. A scaling is exact where the scaled entries stay normal or zero, and RQ.8
needs $m_{22}^2$ normal after it. So RQ.8 holds for every finite triangle whose entries in each row
lie within $2^{894}$ (`f64`), $2^{94}$ (`f32`) of that row's largest, and whose $m_{22}$ lies within
$2^{383}$ (`f64`), $2^{31}$ (`f32`) of its row's largest. Every diagonal meets both, subnormal
entries included. Beyond the second, $u_{22}$ is read as 0. That stays backward stable, because the
row's residual is then below $2^{-383}$ ($2^{-31}$) of its norm, which is less than $u$. Inside the
window, $\lVert m_i\rVert^2$ cannot overflow, and the largest entry's
square is normal (`f64` $[2^{-128}, 2^{128})$, `f32` $[2^{-32}, 2^{32})$, squares within
$[2^{-256}, 3\cdot2^{256})$ and $[2^{-64}, 3\cdot2^{64})$). A row's own entries below its largest
lose digits only below $u$ of it. One global scale cannot do this. On
$\mathrm{diag}(10^{200}, 1, 10^{-200})$ it pushes the last row below the smallest subnormal.

*Proof.* RQ.1's $(DU, R)$ with $D = \mathrm{diag}(1/f_i)$. A product by a power of two is exact
whenever the result is normal or zero; scaling back yields $m_{ij}$, which is representable by
definition. The window's lower end is $2^{-128}$ ($2^{-32}$), the smallest normal $2^{-1022}$
($2^{-126}$). $\square$
