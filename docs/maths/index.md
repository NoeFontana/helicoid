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
| $u$ | unit roundoff (`NUMERICS.md` §2.1), except on `so3.md`, where $u$ is the quaternion's vector part (`NUMERICS.md` §1), $n = \lVert u\rVert$ (not the tangent dimension of the first row) and the roundoff is written $\mathsf u$; `coefficients.md` keeps $u$ for the roundoff and writes the vector part $q_{\mathrm v}$ |
| $\sigma_m$, $\tau_m$, $\rho_m$, $E_x$, $E_s$, $E_\times$, $s_j$ | on `coefficients.md` only, not the tangents, translation tangent or $\mathrm{Exp}$ above: the series families (CO.1), a relative truncation size, the errors of the exact arm, of the series arm and at the crossing, and the coefficient of $z^j$; its symbol note lists the rest |
| $O(\lVert\delta\rVert^2)$ | remainder bounded by $C\lVert\delta\rVert^2$ as $\delta \to 0$, $C$ locally uniform in the base point |

## Labels and the `Checked:` line

- A statement is labelled `<page>.<n>` (`LG` for [`lie-groups.md`](./lie-groups.md), `SO` for
  [`so3.md`](./so3.md), `CO` for [`coefficients.md`](./coefficients.md)), numbered in order of
  appearance. Labels are stable: later results are appended, never renumbered, so other pages and PR
  descriptions can cite them.
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
- The mpmath version is that of the run, stated in the `Checked:` line: `lie-groups.md` ran at 1.4.1,
  `so3.md` and `coefficients.md` at 1.3.0 (the latter with sympy 1.14.0 for the exact rational
  series). It is not a pin, and no claim on these pages depends on it; the one claim about mpmath's
  own behaviour (`mp.logm`, below) was checked on both.
- A wrong variant is characterised by its exact gap to the right formula (`lie-groups.md` §6),
  never by a sampled minimum: every wrong variant tends to the right one as the tangent tends to
  $0$, so a sampled minimum measures the sampler.

## Map

| Page | `NUMERICS.md` | API items (`API.md` §3, `PHASE3.md` §2, §5) |
|---|---|---|
| [`lie-groups.md`](./lie-groups.md) | §1, §2.2, §2.3, §5.1, §5.2, §5.4 (the algebra, given §5.3), §12 (`jr_inv` domain), §14 (`jl`, `*_jacobians`) | `LieGroup::{exp, log, adjoint, ad, jr, jr_inv, jl, jl_inv, rplus, lplus, rminus, lminus, *_jacobians, compose_jacobians, inverse_jacobian}`, `Side`, `SEn3Jac` |
| [`so3.md`](./so3.md) | §1 (quaternion), §2.4, §3.1–§3.6, §4 (definitions of $k, a, b, c, r$; the $c$ series), §11 (`Log` conditioning, matrix input), §12 (`jr_inv`, `from_wxyz_unchecked`), §14 (`act_many`) | `SO3::{exp, log, act, act_many, to_matrix, from_matrix, renormalize, adjoint, ad, jr, jr_inv, jl, jl_inv, act_jacobians}`, `Quat::{from_wxyz_unchecked, from_wxyz_normalized}` |
| [`coefficients.md`](./coefficients.md) | §4 (all seven coefficients and the $\cos\frac\theta2$ of §3.1: series, cancellation, switch-point magnitudes, derivatives), §2.1, §7 ($\Gamma_2$'s coefficients), §12 (the $\theta = 0$, $\pi$, $2\pi$ behaviour); `PHASE1.md` §2 item 3, §4.3, §6 | `coeffs::{exp_coeffs, jr_coeffs, jr_inv_coeff, q_coeffs, gamma2_coeffs, log_ratio, se2_coeffs}` (`pub(crate)`, `PHASE3.md` §3), `coeffs::generated`, `xtask thresholds` |

Not derived yet, so `NUMERICS.md` alone states them: every section outside the Map rows, notably §5.3
(the $Q$ block: `lie-groups.md` uses only $J = g(\mathrm{ad})$ and `coefficients.md` derives its
coefficients $b, d, e$, not the block) and §6–§11 (the SE(2) matrices, of which only the pair
$\alpha, \beta$ is placed in the coefficient family; the conditioning of the exact $J$ and $\mathrm{Ad}$
is LG.16; the metrics of §11 are used, not derived). The switch points and series lengths of §4 are
generated (0004), not derived: `coefficients.md` derives the magnitudes to check them against
(CO.10).

## Open items for the normative documents

- `NUMERICS.md` §12 has no row for the Jacobians of `rminus`, `lminus` and `Log`
  (`rminus_jacobians`, `lminus_jacobians`), which exist only for $\theta(X) < \pi$ (LG.2(c),
  LG.14). `API.md` R6 asks every restricted-domain function to name its §12 entry. Adding the row
  is a `NUMERICS.md` edit and needs a draft record; none is filed.
- `NUMERICS.md` §3.2: the flip is "implemented as `copysign`", yet "at $w = +0$ nothing flips" and
  "$q$ and $-q$ return $\pm\pi\hat n$". A sign-bit `copysign` flips at $w = -0$ (`copysign(1, -0.0) = -1`;
  `w < 0` does not), so $q = (+0, u)$ and $-q = (-0, -u)$ both return $+\pi\hat n$; only the `w < 0`
  reading returns $\pm\pi\hat n$. Both are logarithms (SO.5(d)); §3.2 should pick one.
- `NUMERICS.md` §3.3 and §1 ($R(q)$, used by `to_matrix`, `act_many`) differ by $\eta v$ on a non-unit
  $q$ (SO.2(d)), up to $2^{-40}\lVert v\rVert$ at the `from_wxyz_unchecked` bound: the `act_many` twin
  proptest (§14) needs normalized inputs.
- `PHASE1.md` §4.3 defines `so3_log` by `mp.logm` of $R(q/\lVert q\rVert)$ and `so3_from_matrix` as the
  "`mp.logm`-consistent quaternion". `mp.logm` returns a complex, non-principal result for a rotation of
  angle $\theta \ge 3.03$ (onset in $(3.02, 3.03]$ depending on the axis; 30 axes, 30 and 80 digits,
  identical in mpmath 1.3.0 and 1.4.1), so both ids are wrong at every `theta:pi-1e-k` stratum
  and at the top of `theta:1e0` ($[1, \pi - 0.1)$). An eigendecomposition-based principal logarithm
  agrees with the `atan2` form there (SO.5, `Checked:`).
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
  fixes $m \le 8$; with it, the derivative of $e$ cannot fall below $\sim2\times10^3u$ (`f64`).
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
  ($n^2 = 10^{-6}$, $w = -1$) (CO.16). §8's "$\alpha/\lVert n \times m\rVert$ is §4's $r$" is $r/2$.
- `NUMERICS.md` §12 has no row for the domain of `Exp` (and the other coefficient users), while
  `PHASE2.md` §2 promises a finite result for finite in-domain input: $\theta^2$ overflows for
  $\lVert\varphi\rVert > 1.3\times10^{154}$ (`f64`), $1.8\times10^{19}$ (`f32`) and the exact arm
  returns `NaN` (CO.15).
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
  $2.6\times10^{-8}$, or the fit skip saturated strata.

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
  (ch. 2: Sterbenz's lemma; ch. 3: the $\gamma_n$ bound for inner products; §5.1: Horner's method).
- **[DLMF]** *NIST Digital Library of Mathematical Functions*, https://dlmf.nist.gov, §4.19 and §4.22
  (the series and the partial fractions of $\cot$), §24.2 (Bernoulli numbers).
- **[mpmath]** F. Johansson et al., *mpmath: a Python library for arbitrary-precision
  floating-point arithmetic* (the checking tool; version stated in each `Checked:` line).
