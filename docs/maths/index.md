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
| $O(\lVert\delta\rVert^2)$ | remainder bounded by $C\lVert\delta\rVert^2$ as $\delta \to 0$, $C$ locally uniform in the base point |

## Labels and the `Checked:` line

- A statement is labelled `<page>.<n>` (`LG` for [`lie-groups.md`](./lie-groups.md), `SO` for
  [`so3.md`](./so3.md), `CO` for [`coefficients.md`](./coefficients.md), `SE` for
  [`se3.md`](./se3.md), `PL` (planar) for [`so2-se2.md`](./so2-se2.md), `EA` (error analysis) for
  [`error-analysis.md`](./error-analysis.md)), numbered in order of appearance. Labels are stable:
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
  `libm` 0.2.16 crate; wasmtime 49.0.0 for its wasm32 run). It is not a pin, and no claim on these
  pages depends on it; the claims about mpmath's own behaviour (`mp.logm` of an $\mathrm{SO}(3)$, an
  $\mathrm{SE}_N(3)$ and an $\mathrm{SE}(2)$ matrix, below) were checked on both 1.3.0 and 1.4.1.
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

Not derived yet, so `NUMERICS.md` alone states them: every section outside the Map rows, notably
§7–§10 (of these only $\Gamma_m$'s coefficients are placed, CO.4(d)); §11 is placed by `error-analysis.md`
(EA.3–EA.13) and LG.16 (the conditioning of the exact $J$ and $\mathrm{Ad}$). The switch points and
series lengths of §4 are generated (0004), not derived: `coefficients.md` derives the magnitudes to check
them against (CO.10).

## Open items for the normative documents

The gaps the corpus generator, `helicoid-linalg` and the lint found in the specs, with most of the items
below, are collected for decision, with a recommendation each, in the draft record
[`0015`](../decisions/0015-specification-gaps-found-while-building-the-instrument.md).

- `NUMERICS.md` §12 has no row for the Jacobians of `rminus`, `lminus` and `Log`
  (`rminus_jacobians`, `lminus_jacobians`), which exist only for $\theta(X) < \pi$ (LG.2(c),
  LG.14). `API.md` R6 asks every restricted-domain function to name its §12 entry. Adding the row
  is a `NUMERICS.md` edit and needs a record: it is gap NU.9 of the draft `0015`.
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
  $2.6\times10^{-8}$, or the fit skip saturated strata.
- `0004` item 1 and `PHASE1.md` §6 minimize each coefficient's own value and `Dual`-derivative error, but
  a consumer multiplies a coefficient by a word of size $\theta^p$: with the exact arms of $b, d, e$ alone,
  the value of $Q$ errs by $\lesssim(16.3/\theta + 20.8)u$ (measured $\le6.5u/\theta$, `f64`; $39u$ at a
  switch of $0.9$, measured $\le 9.4u$) and the translation column of `Exp` by $O(u)\lVert\rho\rVert$ at every
  $\theta \le 1$ (SE.15), against $6u\theta^{-2}$ for $b$ and $360u\theta^{-4}$ for $e$. The value strata of
  `sen3_exp_n*`, `sen3_jl_n*` are far less sensitive to the
  switch of $b$ and $e$ than `coeff_b`, `coeff_e`, so a switch that minimizes the coefficient is
  conservative for those values; whether the objective should be the consumer's is a `0004` question. The
  derivative through `Dual` was not examined.
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
  bound for inner products; §5.1: Horner's method).
- **[DLMF]** *NIST Digital Library of Mathematical Functions*, https://dlmf.nist.gov, §4.19 and §4.22
  (the series and the partial fractions of $\cot$), §24.2 (Bernoulli numbers).
- **[SLF]** G. L. Steele Jr., D. Lea, R. Flood, "Fast splittable pseudorandom number generators",
  OOPSLA 2014. **[Vigna]** S. Vigna, `splitmix64.c`, public domain,
  https://prng.di.unimi.it/splitmix64.c (the reference outputs of EA.23).
- **[mpmath]** F. Johansson et al., *mpmath: a Python library for arbitrary-precision
  floating-point arithmetic* (the checking tool; version stated in each `Checked:` line).
