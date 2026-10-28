# The symmetric N×N eigendecomposition: why a relative threshold

> **Non-normative.** [`NUMERICS.md`](../NUMERICS.md) §20 is the contract. This page states the
> results §20 and §11's relative reading rest on, with their sources. It does not re-derive them.
> [`0111`](../decisions/0111-eig-sym-is-eig3-s-jacobi-at-general-n.md)'s convergence study and the
> corpus ids `eig_sym_n{4,9,16}` are what check them numerically.

$u$ is the unit roundoff. $A$ is symmetric $N \times N$ with eigenvalues
$\lambda_1 \le \dots \le \lambda_N$ and orthonormal eigenvectors $v_k$.

**ES.1 (absolute perturbation; Weyl, Davis–Kahan).** If $\tilde A = A + E$, then
$\lvert\tilde\lambda_k - \lambda_k\rvert \le \lVert E\rVert_2$. If $\lambda_k$ is simple with
$\mathrm{gap}_k = \min_{j\ne k}\lvert\lambda_k - \lambda_j\rvert$, then
$\sin\angle(\tilde v_k, v_k) \le \lVert E\rVert_2/\mathrm{gap}_k$ to first order. Within a cluster
only the invariant subspace is determined. Every backward stable method, Jacobi and tridiagonal QR
alike, has $\lVert E\rVert \le c(N)\,u\lVert A\rVert$. §11's absolute reading is this scale.

**ES.2 (relative perturbation; Demmel and Veselić 1992, Theorem 2.3).** Let $A = DHD$ be positive
definite, $D = \mathrm{diag}(\sqrt{a_{ii}})$, $H$ unit-diagonal. Let $\tilde A = D(H + \Delta H)D$
with $\lVert\Delta H\rVert_2 \le \eta < \lambda_{\min}(H)$. Then
$\lvert\tilde\lambda_k - \lambda_k\rvert \le \eta\,\kappa(H)\,\lambda_k$, where
$\kappa(H) \le \kappa(A)$, and for some diagonal $D$ by far: a graded $A$ has $\kappa(A) \sim 10^{16}$
and $\kappa(H) \sim 1$. For the vectors, $\sin\angle(\tilde v_k, v_k) \lesssim \eta\,\kappa(H)/\mathrm{relgap}_k$
with $\mathrm{relgap}_k = \min_{j\ne k}\lvert\lambda_k - \lambda_j\rvert/\sqrt{\lambda_k\lambda_j}$.

**ES.3 (two-sided Jacobi attains ES.2; ibid., §3).** Each rotation in floating point is a
perturbation $\Delta H$ of norm $O(u)$, provided:
- the rotation is computed as §17 computes it (no division by $a_{pq}$);
- the iteration stops on the relative test $\lvert a_{pq}\rvert \le u\sqrt{a_{pp}a_{qq}}$.

The computed eigenvalues are then within $c(N)\,u\,\kappa(H)\,\lambda_k$. An absolute test,
$\lvert a_{pq}\rvert \le u\lVert A\rVert$, may discard an entry that is large relative to
$\sqrt{a_{pp}a_{qq}}$. It then gives only ES.1. QR's tridiagonalisation mixes rows of very different
scale, and gives only ES.1 in general (ibid., §5).

**ES.4 (the test, without underflow).** §20 compares $\lvert a_{pq}\rvert$ with
$u\,r_p r_q$, $r_i = \sqrt{\lvert a_{ii}\rvert}$, not $a_{pq}^2$ with $u^2\lvert a_{pp}a_{qq}\rvert$.
The squared form underflows: at binary32, a graded diagonal of $2^{-91}$ squares to $2^{-182}$,
below the subnormals, and the test then never holds. Each $r_i$ stays normal wherever $a_{ii}$ is.
The floor $\max(\cdot, \mathrm{tiny})$ ends the iteration on two exactly zero diagonals. It discards
at most entries below the smallest normal, which is below $2^{-894}\lVert A\rVert$ inside §17's
binary64 window.

**ES.5 (the vector reading's cap).** A rounded unit vector carries an error of about $u$ in each
component, whatever ES.2 promises. Where $\mathrm{relgap}_k > 1$ the bound of ES.2 would be below
$u$, so §11 caps the relative gap at 1. The reading is then at most the angle in $u$, over $\kappa$.

**ES.6 (convergence).** Cyclic-by-row Jacobi converges (Forsythe and Henrici 1960), and
asymptotically quadratically per sweep (Schönhage 1964; Wilkinson 1962). Near-multiple eigenvalues
slow the early sweeps, not the quadratic phase. No closed bound on the sweep count is tight
enough to use. §20's cap is measured: `0111`'s study, 2e4 matrices of each of 18 families at
every $N$, then a hill climb.
