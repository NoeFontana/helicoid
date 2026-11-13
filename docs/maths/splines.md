# Uniform cumulative B-splines: why the recursions hold

> **Non-normative.** [`NUMERICS.md`](../NUMERICS.md) §21 is the contract. This page derives its
> formulas. Sommer, Usenko, Schubert, Demmel and Cremers, *Efficient derivative computation for
> cumulative B-splines on Lie groups* (CVPR 2020) is the source of the recursions; Qin, *General
> matrix representations for B-splines* (The Visual Computer, 2000) of the basis matrix.

Notation as [`index.md`](./index.md). The window is $X_0, \dots, X_{K-1}$, $d_j = \mathrm{Log}(X_{j-1}^{-1}X_j)$,
$A_j = \mathrm{Exp}(\lambda_j d_j)$, $R_j = \mathrm{Ad}_{A_j^{-1}}$ and $J_j = J_r(\lambda_j d_j)$.
Dots are time derivatives, $t = t_i + h u$.

**BS.1 (the cumulative form).** A B-spline in a vector space is $p = \sum_s b_s c_s$ with
$\sum_s b_s = 1$. Summing by parts, $p = c_0 + \sum_{j \ge 1}\lambda_j (c_j - c_{j-1})$ with
$\lambda_j = \sum_{s \ge j} b_s$. On a group, the difference becomes $d_j$ and the sum becomes the product
$X_0\prod_j\mathrm{Exp}(\lambda_j d_j)$. It is left-invariant, and it is $C^{K-2}$ because each
$\lambda_j$ is.

**BS.2 (the basis as integers).** Qin's matrix,
$m_{s,k} = C(K-1,k)\sum_{l=s}^{K-1}(-1)^{l-s}C(K,l-s)(K-1-l)^{K-1-k}/(K-1)!$, has integer numerators
over $(K-1)!$. For $K \le 6$, every numerator, every derivative numerator and every Gram numerator
and denominator is below $2^{24}$ ([`0115`](../decisions/0115-the-spline-tables-are-exact-integers.md)'s
table). At $K = 4$, $(K-1)!\,\tilde m$ has rows $(6,0,0,0)$, $(5,3,-3,1)$, $(1,3,3,-2)$ and $(0,0,0,1)$.

**BS.3 (velocity).** Let $Y_j = X_0A_1\cdots A_j$ and $\omega^{(j)} = (Y_j^{-1}\dot Y_j)^\vee$.
Then $Y_j = Y_{j-1}A_j$, and $A_j^{-1}\dot A_j = (\dot\lambda_j d_j)^\wedge$, because $A_j$ moves on a
one-parameter subgroup. So
$$\omega^{(j)} = \mathrm{Ad}_{A_j^{-1}}\omega^{(j-1)} + \dot\lambda_j d_j.$$

**BS.4 (acceleration).** Differentiate BS.3. Since $\mathrm{Ad}_{A_j^{-1}} = \exp(-\lambda_j\,\mathrm{ad}_{d_j})$,
its derivative is $-\dot\lambda_j\,\mathrm{ad}_{d_j}R_j$. So
$\dot\omega^{(j)} = R_j\dot\omega^{(j-1)} - \dot\lambda_j\,\mathrm{ad}_{d_j}R_j\omega^{(j-1)} + \ddot\lambda_j d_j$.
Substitute $R_j\omega^{(j-1)} = \omega^{(j)} - \dot\lambda_j d_j$, use $\mathrm{ad}_{d_j}d_j = 0$ and
$-\mathrm{ad}_{d}\,\omega = \mathrm{ad}_{\omega}\,d$:
$$\dot\omega^{(j)} = R_j\dot\omega^{(j-1)} + \ddot\lambda_j d_j + \dot\lambda_j\,\mathrm{ad}\big(\omega^{(j)}\big)d_j.$$

**BS.5 (the derivative of an adjoint action).** By the left Jacobian,
$\mathrm{Exp}(x + \delta) = \mathrm{Exp}(J_l(x)\delta)\,\mathrm{Exp}(x)$ to first order, so
$\partial(\mathrm{Ad}_{\mathrm{Exp}(x)}v)/\partial x = -\mathrm{ad}(\mathrm{Ad}_{\mathrm{Exp}(x)}v)\,J_l(x)$.
At $x = -\lambda d$, with $J_l(-y) = J_r(y)$:
$$\partial(R_j v)/\partial d_j = \lambda_j\,\mathrm{ad}(R_j v)\,J_j.$$
Applied to BS.3 and BS.4, this gives §21's step partials. The last two terms of
$\partial\dot\omega^{(j)}/\partial d_j$ come from $\mathrm{ad}(\omega)d = -\mathrm{ad}(d)\,\omega$, which
is linear in each argument.

**BS.6 (backward accumulation).** Step $j$ maps $(\omega^{(j-1)}, \dot\omega^{(j-1)})$ by
$S_j = \begin{bmatrix}R_j & 0\\ -\dot\lambda_j\,\mathrm{ad}(d_j)R_j & R_j\end{bmatrix}$. A product of such
block lower-triangular maps keeps equal diagonal blocks:
$S_{K-1}\cdots S_{j+1} = \begin{bmatrix}T_j & 0\\ C_j & T_j\end{bmatrix}$, and multiplying on the right by
$S_j$ gives §21's recurrence. $T_j = R_{K-1}\cdots R_{j+1} = \mathrm{Ad}_{P_j^{-1}}$, which is
also the value's carried map. Perturbing $d_j$ by $\delta$ moves
$X = Y_{j-1}A_jP_j$ to $X\,\mathrm{Exp}(\mathrm{Ad}_{P_j^{-1}}\lambda_j J_j\delta)$. One pass forward
and one backward give every Jacobian in $O(K)$ products. A forward pass per control would need
$O(K^2)$, and is the reference twin.

**BS.7 (increments and sides).** $X_j\mathrm{Exp}(\delta)$ changes $d_j$ to
$\mathrm{Log}(\mathrm{Exp}(d_j)\mathrm{Exp}(\delta))$, so $\partial d_j/\partial\delta_j = J_r^{-1}(d_j)$.
$X_{j-1}\mathrm{Exp}(\delta)$ changes it to $\mathrm{Log}(\mathrm{Exp}(-\delta)\mathrm{Exp}(d_j))$, so
$\partial d_j/\partial\delta_{j-1} = -J_l^{-1}(d_j)$. For a left perturbation,
$\mathrm{Exp}(\delta)X_m = X_m\mathrm{Exp}(\mathrm{Ad}_{X_m^{-1}}\delta)$ on the input, and
$\mathrm{Exp}(\epsilon)X = X\,\mathrm{Exp}(\mathrm{Ad}_{X^{-1}}\epsilon)$ on the value.

**BS.8 (Gram matrices).** On $\mathbb R^n$, $p^{(r)}(t) = h^{-r}\sum_s b_s^{(r)}(u)c_s$ and $dt = h\,du$,
so $\int\lVert p^{(r)}\rVert^2dt = h^{1-2r}\sum_{s,t}c_s^\top c_t\int_0^1 b_s^{(r)}b_t^{(r)}du$. The
integral of $u^{p-r}u^{q-r}$ is $1/(p+q-2r+1)$. A constant window, $c_s = c$, has zero
derivatives, so every row sums to zero for $r \ge 1$.

**Checked.** Every formula of §21 against central differences at 50 digits, on SO(3) and SE(3), at
$K = 6$, $u = 0.37$, $h = 1/4$, with random windows of relative angle up to about 1.5 and both
sides. The worst Jacobian entry differs by $2 \times 10^{-25}$. The worst velocity or acceleration
component differs by $8 \times 10^{-23}$, which is the finite differences' own truncation. The Gram
check values are exact rationals. The corpus ids of
[`0116`](../decisions/0116-the-spline-is-measured-before-it-ships.md) are the committed check.
