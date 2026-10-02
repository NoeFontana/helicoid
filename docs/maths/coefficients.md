# The coefficient catalogue: series, cancellation, switch points

> Non-normative companion to [`NUMERICS.md`](../NUMERICS.md) §4, with §2.1 ($u$), §3.1–§3.2 and §3.5
> ($k, a, b, c, r$ at their call sites), §5.3 ($Q$), §7 ($\Gamma_2$), §12,
> [`0003`](../decisions/0003-the-scalar-that-cannot-say-less-than.md),
> [`0004`](../decisions/0004-switch-points-are-generated-not-typed.md) and
> [`PHASE1.md`](../PHASE1.md) §2 item 3, §4.3, §6. **`NUMERICS.md` wins on any conflict; a
> disagreement is an open item in the [maths index](./index.md)**, which holds the notation and the
> `Checked:` convention. Nothing here is a switch point, a series length or a coefficient to type:
> those are generated (0004); the numbers below are magnitudes to sanity-check them against.

On this page $z = \theta^2$ is the branch variable, $u$ the unit roundoff ($2^{-53}$ `f64`,
$2^{-24}$ `f32`), $L = \log_{10}(1/\theta)$, $s_j$ the coefficient of $z^j$ of the function at hand,
$\rho_m = \lvert s_m/s_0\rvert$; the *exact arm* and the *series arm* are the two arms of a branch.
For $r$, $n = \sqrt{n^2}$ and $w$ are the vector norm and the scalar part of a quaternion (as on
[`so3.md`](./so3.md)). Rounding statements are for `f64` unless a precision is named. Measured
constants are for the evaluation orders written in §3 and depend on `sin`/`cos` at the ulp level:
they were measured with glibc, whose `sin`/`cos` are essentially correctly rounded; the `libm` crate
errs by up to $1.09u$ on `sin`, `cos` and $1.96u$ on `atan2` (measured, glibc $\le1.00u$,
[`error-analysis.md`](./error-analysis.md) EA.2(d)), so the transcendental part of its constants can
be up to $\sim1.1\times$ (`sin`, `cos`) or $\sim2\times$ (`atan2`) these.

Symbols that the [index](./index.md) reserves are reused here with a local meaning: $\sigma_m$,
$\tau_m$ are the series families of CO.1 (not tangents); $\rho_m$ is a relative size (CO.8(e); not a
translation tangent); $E_x$, $E_s$, $E_\times$ are errors (not $\mathrm{Exp}$); $\epsilon_i$, $\eta_j$,
$\gamma_j$, $\delta_j$ are relative roundings, $\varepsilon$ an error target, and $\delta = \pi - \theta$
in CO.7 and CO.15; $u$ is always the roundoff, and the quaternion's vector part is written $q_{\mathrm v}$;
$s_j$ (the coefficient of $z^j$) is unrelated to $s = n^2/w^2$ (the variable of $r$; $s_\times$ its
crossing).

## Results

| Label | Result | `NUMERICS.md` | Implemented by |
|---|---|---|---|
| CO.1–CO.3 | $k, a, b, d, e$ are entire in $z$ with general terms $\frac{(-1)^j}{2^{2j+1}(2j+1)!}$, $\frac{(-1)^j}{(2j+2)!}$, $\frac{(-1)^j}{(2j+3)!}$, $\frac{(-1)^j}{(2j+4)!}$, $\frac{(-1)^j(j+1)}{(2j+5)!}$, and $\cos\frac\theta2 = \sigma_0(z/4)$; $c = \sum\lvert B_{2j+2}\rvert z^j/(2j+2)! = 2\sum_{m\ge1}(4\pi^2m^2 - z)^{-1}$, radius $4\pi^2$; $r$ has radius $1$ in $s = n^2/w^2$; the leading four of §4 | §4 | the series arms; the generator's rationals (0004 item 2) |
| CO.4 | $a = 2k^2$; $\sigma_{m-2} = \frac1{(m-2)!} - z\sigma_m$; $e = \frac12\big(d - \frac{1/2 - 3b}{z}\big)$; $\Gamma_m = \frac{I}{m!} + \sigma_{m+1}W + \sigma_{m+2}W^2$; $2ac = a^2 - b\sigma_1$; $rk = 1$; where each is used | §3.5, §4, §5.3, §7 | `exp_coeffs`, `jr_coeffs`, `jr_inv_coeff`, `q_coeffs`, `gamma2_coeffs`, `log_ratio` |
| CO.5–CO.6 | forward error of an exact arm is $\hat C\,u\,\theta^{-p}$: the exponents $p$, a rounding-count bound and the measured constants; the "naive cancellation" column | §4 | the exact arms |
| CO.7 | what the rewrites buy ($a$ everything, $c$ the behaviour at $\pi$ and nothing at small $\theta$, $d$ two powers); why $b, e$ have none (the spec's position) | §4, §12 | exact arms of `jr_coeffs`, `jr_inv_coeff`, `q_coeffs` |
| CO.8–CO.9 | truncation tails (alternating, monotone, Lagrange) with validity ranges; Horner rounding $\le 2.6u$ | §4 | the series arms |
| CO.10–CO.12 | the crossing $\theta_\times = (C_xu/\rho_m)^{1/(2m+p)}$ and predicted magnitudes; no single threshold (e); the jump at a switch | §4 (switch, continuity) | `coeffs::generated`, `xtask thresholds`, `branch_continuity_*` |
| CO.13–CO.14 | $\mathrm d/\mathrm dz$ of each coefficient; through the series arm; through the exact arm the exponent $p$ becomes $p + 2$ | §4, §7 | `Dual` through `coeffs` |
| CO.15 | safe argument, $z = 0$, underflow, overflow, $c$ at $\pi$ and $2\pi$ | §4, §12 | every arm |
| CO.16 | the $r$ kernel: domain, series, $w \to 0$, derivative | §3.2, §4, §8 | `log_ratio` |
| CO.17 | reference precision: digits lost by a definition, `dps = 120` | (`PHASE1.md` §2 item 3) | `conformance/generate` |
| CO.18 | what one shared mask per call-site group costs against a switch per coefficient | §4 (call sites), `PHASE3.md` §3 | `coeffs::generated`, the call-site kernels |

## 1. Definitions and series

**Definition CO.1.** For an integer $m \ge 0$ let
$$
\sigma_m(z) = \sum_{j\ge0}\frac{(-z)^j}{(2j+m)!},\qquad
\tau_m(z) = \sum_{j\ge0}\frac{(j+1)(-z)^j}{(2j+m)!},
$$
entire in $z$ (consecutive terms differ by a factor $O(j^{-2})$). $\sigma_0(\theta^2) = \cos\theta$,
$\sigma_1(\theta^2) = \sin\theta/\theta$.

**Proposition CO.2 (entire coefficients).** For $\theta > 0$ the definitions of `NUMERICS.md` §4 are

$$
k = \tfrac12\,\sigma_1(z/4) = \sum_{j\ge0}\frac{(-1)^jz^j}{2^{2j+1}(2j+1)!},\quad
a = \sigma_2,\quad b = \sigma_3,\quad d = \sigma_4,\quad e = \tau_5,\quad
\cos\tfrac\theta2 = \sigma_0(z/4) = \sum_{j\ge0}\frac{(-1)^jz^j}{4^j(2j)!}
$$

(the last is the series arm that `NUMERICS.md` §3.1 leaves to the generator, "alongside $k$").
Each is an entire function of $z$, hence even and analytic in $\theta$ with infinite radius, with
$(k, a, b, d, e)(0) = (\frac12, \frac12, \frac16, \frac1{24}, \frac1{120})$ and $\cos\frac\theta2 = 1$ at
$0$: **the singularity at $\theta = 0$ is removable** and the series are the extension. The general
terms at $j = 0, \dots, 3$ reproduce every entry of the series column of §4.

*Proof.* Multiply out the series of $\sin$ and $\cos$.
$\theta - \sin\theta = \sum_{j\ge1}(-1)^{j+1}\theta^{2j+1}/(2j+1)!$; divide by $\theta^3$ and put
$j \to j+1$: $\sigma_3$. Likewise $1 - \cos\theta$ gives $\sigma_2$, and
$\theta^2 + 2\cos\theta - 2 = 2\sum_{j\ge2}(-1)^j\theta^{2j}/(2j)!$ gives $\sigma_4$ after
$/2\theta^4$. For $e$: $\frac1{(2j)!} - \frac3{(2j+1)!} = \frac{2j-2}{(2j+1)!}$, so
$-3\sin\theta + \theta\cos\theta = \sum_j(-1)^j(2j-2)\theta^{2j+1}/(2j+1)!$; its $j = 0$ term
$-2\theta$ cancels $2\theta$, its $j = 1$ term vanishes, and $/2\theta^5$ with $i = j - 2$ leaves
$\sum_i(-1)^i(i+1)\theta^{2i}/(2i+5)!$. $k$: the series of $\sin(\theta/2)$ divided by $\theta$;
$\cos\frac\theta2$: the series of $\cos$ at $\theta/2$. The
leading four are the general terms at $j = 0..3$: $k$:
$\frac1{2}, \frac1{2^3\,3!}, \frac1{2^5\,5!}, \frac1{2^7\,7!}$; $a, b, d$: $1/n!$ for
$n = 2, 4, 6, 8$; $3, 5, 7, 9$; $4, 6, 8, 10$; $e$:
$\frac1{120}, \frac2{7!}, \frac3{9!}, \frac4{11!} = \frac1{120}, \frac1{2520}, \frac1{120960}, \frac1{9979200}$.
$\square$

**Proposition CO.3 ($c$ and $r$).**

- (a) $c = \sum_{j\ge0}\lvert B_{2j+2}\rvert\,z^j/(2j+2)!$, with $B_2 = \frac16$,
  $B_4 = -\frac1{30}$, $B_6 = \frac1{42}$, $B_8 = -\frac1{30}$:
  $\frac1{12} + \frac{z}{720} + \frac{z^2}{30240} + \frac{z^3}{1209600}$.
- (b) $c = 2\sum_{m\ge1}(4\pi^2m^2 - z)^{-1}$. So $c$ is analytic except for simple poles at
  $z = 4\pi^2m^2$: **radius $4\pi^2$ in $z$, $2\pi$ in $\theta$**; all coefficients are positive
  ($= 2\zeta(2j+2)/(2\pi)^{2j+2}$), $c$ is increasing and convex on $[0, 4\pi^2)$,
  $c \ge \frac1{12}$, $c(\pi^2) = 1/\pi^2$, and $c \approx 1/(2\pi(2\pi - \theta))$ as
  $\theta \to 2\pi$.
- (c) For $w > 0$, $s = n^2/w^2$:
  $r = \dfrac2w\,\dfrac{\arctan\sqrt s}{\sqrt s} = \dfrac2w\sum_{j\ge0}\dfrac{(-s)^j}{2j+1}$, radius
  exactly $1$ in $s$ (singularities at $s = -1$); the leading four are
  $\frac2w(1 - \frac s3 + \frac{s^2}5 - \frac{s^3}7)$.

*Proof.* (a) $x\cot x = \sum_{j\ge0}(-1)^jB_{2j}(2x)^{2j}/(2j)!$ [DLMF §4.19] at $x = \theta/2$
gives $1 - \frac\theta2\cot\frac\theta2 = -\sum_{j\ge1}(-1)^jB_{2j}\theta^{2j}/(2j)!$, and
$(-1)^{j+1}B_{2j} = \lvert B_{2j}\rvert$. (b) $x\cot x = 1 + 2\sum_m x^2/(x^2 - m^2\pi^2)$ [DLMF
§4.22] at $x = \theta/2$ gives $1 - \frac\theta2\cot\frac\theta2 = 2z\sum_m(4\pi^2m^2 - z)^{-1}$;
divide by $z$. Expanding each term in $z/4\pi^2m^2$ gives positive coefficients;
$\sum_m(4m^2 - 1)^{-1} = \frac12$ (telescoping) gives $c(\pi^2)$; the $m = 1$ term gives the pole.
(c) $\mathrm{atan2}(n, w) = \arctan(n/w)$ for $w > 0$. $\square$

**Checked (CO.1–CO.3):** sympy 1.14.0, exact rationals: the Taylor series of each definition (and of
the rewritten forms of $c$, $d$; and $\cos\frac\theta2$) to $z^{15}$ equals the general terms, odd powers vanish, the
leading four of all six match the §4 table; $r$: 8 terms. mpmath 1.3.0: (b) against the definition,
$10^{-15}$ (sum of $2\cdot10^4$ terms plus an integral tail, tail-limited),
$c(\pi) - 1/\pi^2 = 10^{-62}$; $c$ increasing on $(0, 2\pi - 0.01]$ with minimum $\frac1{12}$.
**Permanent:** planned, the generator asserts the leading four (0004 item 2); corpus
`coeff_k`…`coeff_r` for the values.

## 2. Identities between them, and where each is used

**Proposition CO.4.**

- (a) $\sigma_{m-2} = \frac1{(m-2)!} - z\,\sigma_m$ for $m \ge 2$. Hence $\cos\theta = 1 - za$,
  $\frac{\sin\theta}\theta = 1 - zb$, $a = \frac12 - zd$: the definitions of $a$, $b$ and $d$ are
  these three lines solved for the coefficient.
- (b) $\tau_m = \frac12\big(\sigma_{m-1} - (m-2)\sigma_m\big)$; so $e = \frac12(d - 3\sigma_5)$ and
  $e = \frac12\big(d - \frac{1/2 - 3b}{z}\big)$.
- (c) $a = 2k^2$ and $\frac{\sin\theta}\theta = 2k\cos\frac\theta2$.
- (d) $\Gamma_m(\varphi) = \sum_{n\ge0}W^n/(n+m)! = \frac{I}{m!} + \sigma_{m+1}W + \sigma_{m+2}W^2$.
- (e) $2ac = a^2 - b\,\sigma_1 = a^2 - b(1 - zb)$.
- (f) For a unit $q$ with $w \ge 0$: $r(n^2, w) = \theta/\sin\frac\theta2 = 1/k$.

*Proof.* (a) $z\sigma_m = \sum_{j}(-1)^jz^{j+1}/(2j+m)! = -\sum_{i\ge1}(-1)^iz^i/(2i+m-2)!$. (b)
$j + 1 = \frac12[(2j+m) - (m-2)]$ and $(2j+m)/(2j+m)! = 1/(2j+m-1)!$. (c)
$1 - \cos\theta = 2\sin^2\frac\theta2$, $\sin\theta = 2\sin\frac\theta2\cos\frac\theta2$. (d)
$W^{2j+1} = (-z)^jW$, $W^{2j+2} = (-z)^jW^2$ (SO.3). (e) the Cramer solution in SO.8. (f)
$n = \sin\alpha$, $w = \cos\alpha$, $\alpha = \theta/2$. $\square$

| Coefficient | Call site (`PHASE3.md` §3) | Where it enters |
|---|---|---|
| $k$, $\cos\frac\theta2$ | `exp_coeffs` | $q = (\cos\frac\theta2, k\varphi)$, §3.1 |
| $a$, $b$ | `jr_coeffs` | $J_{r,l} = I \mp aW + bW^2$, §3.5; $V = J_l$ in $\mathrm{Exp}$ of $\mathrm{SE}_N(3)$, §5.1 |
| $c$ | `jr_inv_coeff` | $J_{r,l}^{-1} = I \pm \frac12W + cW^2$, §3.5; $V^{-1}$ in $\mathrm{Log}$, §5.1 |
| $b, d, e$ | `q_coeffs` | the $Q$ block, §5.3 |
| $b, d$ | `gamma2_coeffs` | $\Gamma_2 = \frac12I + bW + dW^2$, §7 (CO.4(d) at $m = 2$) |
| $r$ | `log_ratio` | $\varphi = r\,q_{\mathrm v}$, §3.2 |
| $\sigma_1$, $\theta a$ | `se2_coeffs` | $\alpha = \sin\theta/\theta$, $\beta = (1 - \cos\theta)/\theta$, §6 |

**Checked (CO.4):** sympy, exact rational series to $z^{15}$ for (a), (b), (d), (e) and $\tau_m$ for
$m = 4..7$; (c) exact; (f) the 9-term series of $r$ in $s = \tan^2\frac\theta2$ against
$\theta/\sin\frac\theta2$ is $O(\theta^{17})$. **Permanent:** planned, corpus `coeff_*`,
`so3_gamma2` (`PHASE5.md` §4).

## 3. Rounding error of the exact arms

*Model.* Every $+, -, \times, \div, \sqrt{\ }$ is correctly rounded,
$\mathrm{fl}(x\circ y) = (x\circ y)(1 + \delta)$, $\lvert\delta\rvert \le u$; `sin`, `cos` likewise.
The argument is an exact input $\hat\theta = \mathrm{fl}(\sqrt z)$, and every operand of a row below
is a function of that one $\hat\theta$ ($\hat\theta^2$ is $\mathrm{fl}(\hat\theta\hat\theta)$): the
rounding of $\hat\theta$ then acts through the conditioning $\theta f'/f = 2zf_z/f = O(z)$ of the even
coefficients, negligible for $\theta \ll 1$. It does not act so if operands of a cancelling sum are
mixed between $z$ and $\hat\theta$ (end of CO.6).

**Lemma CO.5 (sums).** Let $N = \sum_it_i$ with $\hat t_i = t_i(1 + \epsilon_i)$,
$\lvert\epsilon_i\rvert \le \nu_iu$ ($\nu_i$ counts the roundings feeding $t_i$; an input or an
exact scaling has $\nu = 0$). To first order
$\lvert\hat N - N\rvert/\lvert N\rvert \le u\,\kappa_{\rm bd}$,
$\kappa_{\rm bd} = \sum_i\nu_i\lvert t_i\rvert/\lvert N\rvert$, plus one rounding per *inexact*
partial sum; a sum of two operands of opposite sign within a factor $2$ of each other is exact
(Sterbenz [Higham, ch. 2]). Multiplying or dividing by other factors adds $O(u)$, not amplified. A
sharper count uses that $\mathrm{fl}(x)$ errs by at most half the spacing of the binade of $x$,
$u\lvert x\rvert/\mu$ for a mantissa $\mu \in [1, 2)$.

*Proof.* Expand $\prod(1 + \epsilon_i)$ to first order and sum. $\square$

**Proposition CO.6 (the exact arms of the catalogue).** For the evaluation orders below, as
$\theta \to 0$ the relative error is at most $\kappa_{\rm bd}\,u\,\theta^{-p}(1 + O(\theta^2))$ to
first order, and its sampled maximum is $\hat C\,u\,\theta^{-p}$, with:

| | evaluated as | rounded operands (size $\times\,\nu$) | $p$ | $\kappa_{\rm bd}\,\theta^p$ | sampled max $\hat C$ | median | §4 column |
|---|---|---|---|---|---|---|---|
| $a$, naive | $(1 - \cos\theta)/\theta^2$ | $\cos\theta$ ($1 \times 1$) | 2 | 2 (sharp: 1) | 1.0 | 0.5 | 2 |
| $b$ | $(\theta - \sin\theta)/\theta^3$ | $\sin\theta$ ($\theta \times 1$) | 2 | 6 (sharp) | 6.0 | 2.1 | 6 |
| $c$ | $\theta^{-2} - \frac{1 + \cos\theta}{2\theta\sin\theta}$, or $\frac{\cos(\theta/2)}{2\theta\sin(\theta/2)}$ | $\theta^{-2}$ ($\times 2$), second term ($\times 4$) | 2 | 72 | 48 / 47 | 8.5 | 12 |
| $d$, naive | $(\theta^2 + 2\cos\theta - 2)/2\theta^4$ | $2\cos\theta$ ($2 \times 1$), $\mathrm{fl}(\theta^2 + 2\cos\theta)$ ($2 \times 1$) | 4 | 48 (sharp: 36) | 36.0 | 12 | 24 |
| $d$, rewritten | $(\theta^2 - 4\sin^2\frac\theta2)/2\theta^4$ | $\theta^2$ ($\theta^2 \times 1$), $4\sin^2$ ($\theta^2 \times 3$) | 2 | 48 | 45 | 9 | 24 |
| $e$ | $(2\theta - 3\sin\theta + \theta\cos\theta)/2\theta^5$ | $3\sin\theta$ ($3\theta \times 2$), $\theta\cos\theta$ ($\theta \times 2$) | 4 | 480 | 367 | 84 | 360 |
| $k$; $a = 2k^2$ | $\sin(\theta/2)/\theta$; $2k\cdot k$ | none: products and quotients of accurate factors | 0 | $\le 2$; $\le 5$ | 1.5; 3.3 | | none |

*Proof.* The exponent is the gap between the size $\theta^\alpha$ of the rounded operands and the
size $\theta^\beta$ of $N$ (the leading term of CO.2, CO.3): $p = \beta - \alpha$: $b$: $\theta$
against $\theta^3$; $c$: $\theta^{-2}$ against $1$; naive $d$: $1$ against $\theta^4$; $e$: $\theta$
against $\theta^5$, via
$2\theta - 3\sin\theta + \theta\cos\theta = 2(\theta - \sin\theta) - (\sin\theta - \theta\cos\theta)$,
two order-3 remainders ($\approx\theta^3/3$ each) that cancel again. Every subtraction is exact
(Sterbenz) except the partial sum $\theta^2 + 2\cos\theta$ of naive $d$. Counting roundings per
Lemma CO.5 gives the $\kappa_{\rm bd}$ column; e.g. $a$: $\cos\theta \in [\frac12, 1)$ errs by
$\le u/2$ and $N \approx \theta^2/2$, so $C = 1$; naive $d$: $2\cos\theta$ errs by $\le u$,
$\theta^2 + 2\cos\theta \in [2, 4)$ has spacing $4u$ so its rounding is $\le 2u$, total $3u$ over
$\theta^4/12$: $36$. $k$: $\sin$ and $\div$, $\le 2u$; $2k^2$: $4u + u$. The naive column of §4 is
not built by one rule. For $a$, $b$, $c$, $d$ it counts the dominant rounded operand once ($\cos\theta$,
$\sin\theta$, $\theta^{-2}$, $2\cos\theta$, over $N$: $2, 6, 12, 24$); $e$'s $360$ is not fitted but
counts all its operands, $(2\theta + 3\theta + \theta)/(\theta^5/60) = \kappa_{\rm sum}$ of CO.17
(equally $3\sin\theta$ counted twice, for $\sin$ and the $\times3$), and the same all-operands count
gives $48$ for $d$, not $24$. The measured maxima agree with the exponents throughout and are all
below $\kappa_{\rm bd}$. $\square$

*Numerical consequences.* Digits kept are $-\log_{10}(\hat Cu\theta^{-p})$: $b$ keeps 9.2 digits at
$\theta = 10^{-3}$ and 7.2 at $10^{-4}$ (`PROJECT.md` §2: "~7"); $e$ keeps 5.4 at $10^{-2}$, 1.4 at
$10^{-3}$ and none below $\theta = (\hat Cu)^{1/4} = 4.5\times10^{-4}$; $b$ evaluates to $0$ (all
digits lost) for $\theta \lesssim \sqrt{6u} = 2.6\times10^{-8}$. The column of §4 is a magnitude,
not a bound: the sampled maxima exceed it by a factor $4$ ($c$), $1.5$ ($d$ naive), $1.9$ ($d$
rewritten) and $1.02$ ($e$) and are half of it for $a$; a digit count moves by less than one digit
(open item in the index). Consistency of the argument matters: $z$ must not replace
$\mathrm{fl}(\hat\theta\hat\theta)$ in an operand of a cancelling sum (the numerator of $d$, the
$\theta^{-2}$ of $c$). The mismatch $z - \hat\theta^2 \sim uz$ (against $N \sim z^2/12$
in the numerator of $d$; as $u/z$ against $N \sim \frac1{12}$ in $c$) is an error of order $12u/z$, and raises $\hat C$ ($d$ rewritten $45 \to 58$, $6\times10^4$ samples; $c$
$43 \to 50$, $2\times10^4$); in a pure divisor ($2z^2$ against $2\hat\theta^4$) no difference was
measured.

**Proposition CO.7 (what each rewrite buys).**

- (a) $a = 2k^2$ removes the cancellation entirely: $1 - \cos\theta = 2\sin^2\frac\theta2$ is an
  identity between values of one function, and the small quantity becomes a product of accurate
  factors ($2u/\theta^2 \to$ at most $5u$).
- (b) $c$ by the $\cot$ form leaves the small-$\theta$ exponent unchanged ($\hat C$ 48 against 47)
  and fixes $\theta \to \pi$: with $\delta = \pi - \theta$ the definition form has
  $1 + \cos\theta = \delta^2/2$ with absolute error $\le u/2$ ($\cos\theta$ lies in the binade
  $[\frac12, 1)$ in absolute value), divided by $2\theta\sin\theta \approx 2\pi\delta$, against
  $c(\pi) = 1/\pi^2$: relative error $\le \pi u/(4\delta)$, attained (measured maximum
  $0.785\,u/\delta$, $\delta$ the actual distance to $\pi$, from $10^{-2}$ to $10^{-8}$: $7.9\times10^3u$
  at $\delta = 10^{-4}$, $7.9\times10^5u$ at $10^{-6}$, $7.8\times10^7u$, 8 digits, at $10^{-8}$;
  quoted at a nominal $\delta$ over a window reaching $\delta/2$ it reads $\pi/2$; below
  $\sim10^{-8}$ $1 + \cos\theta$ rounds to $0$ and the term, of size $\delta/4\pi$, is dropped). The
  $\cot$ form errs by $\le 2.3u$ for every $\delta$ measured.
- (c) $d$: $2\cos\theta - 2 = -4\sin^2\frac\theta2$ removes the constant pair and lowers $p$ from
  $4$ to $2$; what remains, $\theta^2 - 4\sin^2\frac\theta2 \approx \theta^4/12$, is a polynomial
  against a trigonometric value, the shape of $b$.
- (d) $b$ and $e$ have no rewrite in the spec. **This is the spec's position, not a theorem.**
  (a)–(c) are product and half-angle identities that turn a difference of trigonometric values (or a
  constant against one) into a product. The cancelling operands of $b$ ($\theta$ against
  $\sin\theta$) and of $e$ ($2\theta$, $\theta\cos\theta$, $3\sin\theta$) are polynomial against
  trigonometric: no closed-form identity turns $\theta - \sin\theta$ into a product, and `libm` has
  nothing for it (only `expm1`, `log1p`). A reconstruction by argument reduction does exist,
  $\sin3y = 3\sin y - 4\sin^3y$ giving $b = (3(y - \sin y) + 4\sin^3y)/\theta^3$ with $y = \theta/3$,
  whose terms do not cancel (max error $\approx7u$, 3000 samples in each of $[10^{-8}, 10^{-3}]$,
  $[10^{-3}, 10^{-1}]$, $[10^{-1}, 0.25]$; no growth in $1/\theta$), but it needs $y - \sin y$ from
  the series arm at $y$: a series-plus-recurrence design, which is why the spec keeps the series.
  CO.4(b) gives $e$ from $b$ and $d$, but the factor $1/z$ re-amplifies their errors by $1/z$: one
  cancellation is traded for another. The remedy is the series arm, which is why these two have the
  largest switch points (CO.10).

**Checked (CO.5–CO.7):** mpmath 1.3.0, 150 digits, references from the definitions; `f64` (CPython
floats, glibc `sin`, `cos`); 300 log-uniform $\theta$ per decade over $[10^{-9}, 1]$ for all forms:
the per-decade maximum of $\text{error}/(u\theta^{-p})$ is flat (within a factor $1.4$) from $10^{-7}$
to $1$ for $a$, $b$, $c$, $d$ rewritten, from $3\times10^{-4}$ for naive $d$ and from $10^{-8}$ for
$e$. The error saturates at exactly $1$ (the computed value is $0$) below $\sim10^{-8}$ for $a$, $b$,
$c$, $d$ rewritten and below $\sim3\times10^{-4}$ for naive $d$ (its numerator rounds to $0$); $e$ never
saturates, its relative error keeps growing like $\theta^{-4}$ ($3.7\times10^{10}$ at $10^{-6}$). Then
$6\times10^4$ over $[10^{-3}, 10^{-1}]$ ($[10^{-2}, 10^{-0.5}]$ for naive $d$) for the maxima (a
re-run gave $a$ 1.01, $b$ 5.96, $c$ 46.8 and 43.3, naive $d$ 36.3, rewritten $d$ 44.4, $e$ 362:
maxima are stable to about $10\%$); near $\pi$, 3000 samples per decade of $\delta$ over
$[10^{-10}, 10^{-1}]$ ($\delta$ the actual distance to $\pi$). Sampled maxima, not proved bounds,
except where "sharp" (attained). **Permanent:** planned, corpus `coeff_a`…`coeff_e` (strata
`theta:1e-k`); the seeded "`b` by its definition" (`PHASE1.md` §10: strata `theta:1e-8`…`theta:1e-2`,
fit $\theta^{-p}$, $p \in [1.8, 2.2]$) covers $b$ only, and its window contains the saturated stratum
`theta:1e-8` ($b = 0$ below $\sqrt{6u} = 2.6\times10^{-8}$, error exactly $1$ instead of $\sim6.7$):
in three simulated runs (64 samples per stratum, per-stratum maximum, least squares in
$\log\theta$) the fitted $p$ is $1.91$–$1.94$, inside the window but biased low by $0.07$ ($2.00$ without
that stratum), so a tighter window would fail (open item in the index). Implemented as
`seeded:b-no-series` (`just conformance --self-test`; test
`the_window_the_decades_and_the_measured_exponent_are_pinned`): on the committed corpus, the `value`
field's per-stratum maximum over the seven strata fits $p = 1.937$ (correct kernel: $-0.10$), as
predicted; the record maximum over `value` and `d_branch` fits $3.98$ (open question 7 of 0014 (draft)).

## 4. The series arm

**Proposition CO.8 (truncation).** Let $f \in \{k, a, b, d, e, \cos\frac\theta2\}$, $P_m = \sum_{j<m}s_jz^j$,
$R_m = f - P_m$.

- (a) *Alternating.* $R_m$ has the sign of $s_mz^m$ and $\lvert R_m\rvert \le \lvert s_m\rvert z^m$
  when $z < Z_m$: $a$: $(2m+3)(2m+4)$; $b$: $(2m+4)(2m+5)$; $d$: $(2m+5)(2m+6)$; $k$:
  $4(2m+2)(2m+3)$; $\cos\frac\theta2$: $4(2m+1)(2m+2)$; $e$: $(m+1)(2m+6)(2m+7)/(m+2)$. For $m \ge 1$,
  $Z_m \ge 30 > \pi^2$.
- (b) *Lagrange* (no monotonicity assumed). $R_m = f^{(m)}(\xi)z^m/m!$ for some $\xi \in (0, z)$,
  and for $f = \sigma_M$ ($a, b, d$: $M = 2, 3, 4$)
  $\lvert R_m\rvert \le \dfrac{z^m}{(2m+M)!}(1 - \varrho)^{-(m+1)}$,
  $\varrho = \dfrac{z}{(2m+M+1)(2m+M+2)} < 1$; for $k = \frac12\sigma_1(z/4)$ the same with $z/4$
  and a factor $\frac12$, and for $\cos\frac\theta2 = \sigma_0(z/4)$ with $z/4$ and $M = 0$.
- (c) *Monotone.* For $c$, with $t_m = \lvert B_{2m+2}\rvert z^m/(2m+2)!$ and $z < 4\pi^2$:
  $t_m \le R_m \le t_m/(1 - z/4\pi^2)$.
- (d) For $r$, relative to $2/w$: $\lvert R_m\rvert \le s^m/(2m+1)$ for $0 < s \le 1$.
- (e) Relative size $\rho_m$: $k$: $\frac1{4^m(2m+1)!}$; $\cos\frac\theta2$: $\frac1{4^m(2m)!}$, $2m+1$ times
  $k$'s ($17\times$ at $m = 8$); $a$: $\frac2{(2m+2)!}$; $b$:
  $\frac6{(2m+3)!}$; $d$: $\frac{24}{(2m+4)!}$; $e$: $\frac{120(m+1)}{(2m+5)!}$; $c$:
  $\frac{12\lvert B_{2m+2}\rvert}{(2m+2)!} \sim 24(2\pi)^{-2m-2}$; $r$: $\frac1{2m+1}$ (in $s$). The
  relative truncation error is $\rho_mz^m(1 + O(z))$.

*Proof.* (a) The terms $\lvert s_j\rvert z^j$ decrease for $j \ge m$ once
$\lvert s_{m+1}\rvert z \le \lvert s_m\rvert$ (the ratio $\lvert s_j/s_{j+1}\rvert$ grows with $j$),
and an alternating series with decreasing terms is bounded by its first term and has its sign. (b)
Taylor's theorem; $\sigma_M^{(m)}(\xi)/m! = (-1)^m\sum_i(-1)^i\binom{i+m}{m}\xi^i/(2i+2m+M)!$, with
$(2i+2m+M)! \ge (2m+M)!\,[(2m+M+1)(2m+M+2)]^i$ and
$\sum_i\binom{i+m}{m}\varrho^i = (1-\varrho)^{-m-1}$. (c)
$t_{j+1}/t_j = \frac{\zeta(2j+4)}{\zeta(2j+2)}\frac z{4\pi^2} \le \frac z{4\pi^2}$ ($\zeta$
decreases): a geometric tail over positive terms. (d) alternating, and $s(2j+1)/(2j+3) < 1$.
$\square$

*Consequences.* At $z \le 1$ every condition holds for every $m \ge 1$; the sign of the truncation
error is known (used in CO.12); the bound is tight as $z \to 0$. $r$ converges the slowest
($1/(2j+1)$, not factorially) and $c$ next ($(2\pi)^{-2j}$).

**Proposition CO.9 (Horner rounding).** Evaluate $\hat p_{m-1} = \hat s_{m-1}$,
$\hat p_j = \mathrm{fl}(\hat s_j + \mathrm{fl}(z\hat p_{j+1}))$ with $\hat s_j$ the correctly
rounded literal of $s_j$ ($\lvert\hat s_j - s_j\rvert \le u\lvert s_j\rvert$; measured $\le 0.78u$).
To first order
$$
\lvert\hat p_0 - P_m(z)\rvert \le 2u\sum_{i<m}(i+1)\lvert s_i\rvert z^i,\qquad\text{so}\qquad
\frac{\lvert\hat p_0 - P_m\rvert}{\lvert P_m\rvert} \le 2u\,\frac{1 + \sum_{i\ge1}(i+1)\rho_iz^i}{1 - \sum_{i\ge1}\rho_iz^i} \le 2.6\,u
\quad (z \le 1,\ k, a, b, c, d, e).
$$

*Proof.* $\hat p_j = (\hat s_j + z\hat p_{j+1}(1 + \eta_j))(1 + \delta_j)$ gives the first-order
recursion $E_j = s_j\gamma_j + zp_{j+1}\eta_j + p_j\delta_j + zE_{j+1}$
($\lvert\gamma_j\rvert, \lvert\eta_j\rvert, \lvert\delta_j\rvert \le u$); unroll, and use
$z^j\lvert p_j\rvert \le T_j = \sum_{i\ge j}\lvert s_i\rvert z^i$,
$\sum_jT_j = \sum_i(i+1)\lvert s_i\rvert z^i$. The bound at $z = 1$ is $2.57u$ for $a$ and below
$2.4u$ for the others. $\square$

The series arm has no cancellation of its own in $z \le 1$: the terms fall by a factor $\ge 12$
each, so the floor of the series arm is $\approx 2u$, not $\theta^{-p}u$. The $2.6u$ is proved for
$z \le 1$ only (and is $2.9u$ for $\cos\frac\theta2$ at $z = 1$, where its value has fallen to
$0.88$). Beyond, the alternating series does cancel and the bound is the first-order ratio
$2\sum_{i<m}(i+1)\lvert s_i\rvert z^i/\lvert P_m\rvert$ itself ($m = 8$): $2.1$–$3.2u$ at $\theta = 1.4$
and $2.3$–$6.1u$ at $2.8$, except $12u$ for $a$.

**Checked (CO.8–CO.9):** mpmath 200 digits, remainders from the *definitions* against $P_m$ from the
exact rationals: $k, a, b, c, d, e, \cos\frac\theta2$, $m = 1..10$, 60 log-uniform
$z \in [10^{-6}, 0.999Z_m]$ each: 0 violations of (a), (b), (c) or of the sign; $\max\lvert R_m\rvert/\lvert s_mz^m\rvert = 1.0000$ (as
$z \to 0$); $r$: 800 points, 0 violations. $Z_m$ is sufficient, not shown necessary (no
counterexample for $m \le 4$ up to $16Z_m$). `f64` Horner against the exact $P_m$, 3000 log-uniform
$z \in [10^{-8}, 1]$ per family and $m \in \{2, 4, 6, 8\}$: at most $1.33u$ (value) and $2.06u$ (the
derivative of dual Horner, CO.13(b)); the general ratio above evaluated at the stated $\theta$ from the
exact rationals. **Permanent:** the sweep CSV
(`conformance/sweeps/thresholds.csv`, `PHASE1.md` §6); nothing for the bounds themselves.

## 5. Where the arms meet

**Proposition CO.10 (crossing).** Let the exact side have error $E_x = C_xu\,\theta^{-p}$ (CO.6) and
the series side $E_s = \rho_m\theta^{2m} + \gamma u$, $\gamma \le 2.6$ (CO.8, CO.9: $\gamma \le 2.6$ is
proved for $z \le 1$ only). $E_x$ decreases
and $E_s$ increases in $\theta$, so the switch that minimizes the larger of the two errors over both
sides is their crossing
$$
\theta_\times = \Big(\frac{C_xu}{\rho_m}\Big)^{\frac1{2m+p}},\qquad
E_\times = C_xu\,\theta_\times^{-p} = (C_xu)^{\frac{2m}{2m+p}}\rho_m^{\frac p{2m+p}}
$$
(the floor $\gamma u$ ignored). A switch a factor $\lambda > 1$ *below* $\theta_\times$ costs
$\lambda^{+p}$ in error ($C_xu\theta^{-p}$ at $\theta_\times/\lambda$ is $E_\times\lambda^p$), a factor
$\lambda$ *above* costs $\lambda^{2m}$ ($218\times$ at $\lambda = 1.4$, $m = 8$). For the derivative
(CO.13, CO.14) the exact side is $C'_xu\,\theta^{-(p+2)}$ and the series side, relative to
$\lvert s_1\rvert$, is $m\rho'_m\theta^{2m-2}$ (differentiation lowers the power of $z$ by one),
$\rho'_m = \lvert s_m/s_1\rvert$, so
$$
\theta'_\times = \Big(\frac{C'_xu}{m\rho'_m}\Big)^{\frac1{2m+p}},\qquad
E'_\times = (C'_xu)^{\frac{2m-2}{2m+p}}(m\rho'_m)^{\frac{p+2}{2m+p}},
$$
the penalties being $\lambda^{p+2}$ below and $\lambda^{2m-2}$ above (substituting $C \to C'$,
$p \to p+2$, $\rho \to m\rho'$ into the boxed formulas would give the exponent $2m + p + 2$ and wrong
switch points). Raising $m$ raises $\theta_\times$ and lowers $E_\times$ (down to the floor $\gamma u$);
$\theta_\times \propto u^{1/(2m+p)}$, so `f32` against `f64` moves it by $2^{29/(2m+p)}$, about
$3\times$ at $m = 8$, $p = 2$.

Predicted against measured, `f64`, objective = the derivative (it is the larger, CO.14). $C_x$ from
CO.6 ($6, 48, 45, 370$), $C'_x$ from CO.14; "meas." is a scratch mini-sweep (cells of $1/30$ decade,
25 samples, $\theta \in [10^{-6}, \pi - 0.1]$, the evaluation orders of §3; it understates sampled
maxima by up to $1.5\times$). Neither column is a generated switch point.

| | $m = 4$: $\theta_\times$ pred / meas | $E_\times/u$ pred / meas | $m = 8$: $\theta_\times$ | $E_\times/u$ | meas., $m = 8$, $\theta_s \le 1$ |
|---|---|---|---|---|---|
| $k$ | 0.088 / 0.089 | $4.3\text{e}3$ / $3.7\text{e}3$ | 1.46 / 1.31 | 15 / 11 | 24 |
| $a$ | 0.059 / 0.056 | $1.0\text{e}4$ / $8.0\text{e}3$ | 0.88 / 0.89 | 46 / 32 | 32 |
| $b$ | 0.136 / 0.141 | $6.7\text{e}5$ / $6.7\text{e}5$ | 1.07 / 1.04 | 176 / 115 | 172 |
| $c$ | 0.155 / 0.153 | $7.6\text{e}6$ / $5.3\text{e}6$ | 0.77 / 0.76 | $1.2\text{e}4$ / $8.8\text{e}3$ | $8.8\text{e}3$ |
| $d$ | 0.184 / 0.178 | $2.1\text{e}6$ / $1.9\text{e}6$ | 1.30 / 1.21 | 840 / 630 | $1.2\text{e}3$ |
| $e$ | 0.284 / 0.282 | $3.8\text{e}7$ / $3.0\text{e}7$ | 1.38 / 1.41 | $2.9\text{e}3$ / $2.2\text{e}3$ | $1.7\text{e}4$ |

The *value* objective alone (only $b, c, d, e$ have one, $p > 0$) gives $m = 8$:
$\theta_\times = 1.16, 0.85, 1.41, 1.46$ with $E_\times = 4.5, 67, 23, 81\,u$ (measured
$1.04, 0.83, 1.41, 1.41$ and $3.4, 42, 13, 49\,u$). For $r$ (CO.16): $m = 8$ predicted
$s_\times = 0.012$, $E_\times = 494u$ (measured $0.010$, $332u$). For $\cos\frac\theta2$ ($p = p' = 0$:
no cancellation, a floor of $\le 0.6u$ for $\theta \le 1$ and up to $1.9u$ at $\theta = 2.5$, where the
rounding of $\hat\theta$ acts, for the value and for the derivative): $m = 8$ $\theta_\times \approx 1.3$
(value), $1.1$ (derivative; measured $1.04$, $1.8u$), $m = 4$: $0.08$ and $0.02$ (measured: the exact
arm is at the floor throughout, only $\theta = 0$ needs the series); its $\rho_m$
is $2m+1$ times that of $k$, so it, not $k$, limits `exp_coeffs`'s derivative. In `f32` ($m = 8$, predicted
$\theta_\times$, $E_\times/u$): $k$ 5.1, 1.3; $a$ 3.1, 3.8; $b$ 3.3, 2.0; $c$ 2.4, 142; $d$ 4.0,
9.7; $e$ 3.8, 7.0. Measured: the series arm serves all of $[10^{-3}, \pi - 0.1]$ for $b, d, e$
(errors $2$–$3u$), to $2.8$ for $k$, $a$, and to $2.1$ ($93u$) for $c$.

*Consequences.* (i) The predictions land within $15\%$ in $\theta_s$ and $1.8\times$ in error (the
sampled-maximum constants make them the larger), so a generated $\theta_s$ off by more than
$\sim30\%$, or an objective off by more than $\sim3\times$, beyond the effect of the grid, is a
reason to look at the sweep, not a proof of a defect. (ii) **For $k$, $b$, $d$, $e$ and, in `f32`,
for all six, the optimum lies above the $\theta = 1$ end of the switch grid of `PHASE1.md` §6
(spanning $\theta \in [10^{-8}, 1]$)**; the last column shows what a switch capped at $1$ costs
(`f64`, $m = 8$: $e$: $2.2\text{e}3 \to 1.7\text{e}4u$ on the derivative, $49 \to 290u$ on the
value; `f32`: $e$: $2\to1.6\text{e}4u$, $d$: $2.4 \to 1.8\text{e}3u$ on the derivative); it is not a
defect of the sweep but of the grid (open item in the index). The predictions above $\theta = 1$ use
$\gamma$ beyond its proof (CO.9): they are extrapolations, checked only by the mini-sweep.

**Corollary CO.11 (one threshold is wrong by construction; $e$).** $e$'s exact arm errs by
$367u/\theta^4$ (value) and $2.0\times10^4u/\theta^6$ (derivative). At the library folklore
$\theta_s = 10^{-8}$ it has no correct digit for $\theta < 4.5\times10^{-4}$; at `tf_tree` D12's
$\theta_s = 0.1$ the exact arm above $0.1$ keeps 9.7 digits of the value ($1.6\text{e}6u$ measured)
and 6 of the derivative ($9\text{e}9u$), against 15 of the value for $a$ and $k$; against $e$'s own
optimum at $m = 8$ ($49u$, $2.2\text{e}3u$) that is $3\times10^4\times$ and $4\times10^6\times$. Its
exact arm reaches $\varepsilon = 10u$ only for $\theta \ge (367/10)^{1/4} = 2.5$, beyond what an
$8$-term series serves. And the thresholds conflict: at $m = 8$, $e$ wants $\theta_s \approx 1.4$,
$a$'s value error is flat at its floor ($p = 0$, $3.3u$) up to $\theta \approx 1.0$, where the
truncation reaches it ($C = 3.3$ in CO.10: $\theta_\times = 1.01$), and its derivative wants $0.9$.
Moving $a$ to $1.4$ multiplies its value error by $\approx190$ ($\lambda^{2m}$, $\lambda = 1.4/1.01$;
$240$ measured) and its derivative error by $\approx670$ ($\lambda^{2m-2}$, $\lambda = 1.4/0.88$;
$690$ measured); moving $e$ down to $0.5$ multiplies its value error by $73$ ($\lambda^{+p}$,
$\lambda = 1.46/0.5$, $p = 4$; $77$ measured) and its derivative error by $500$ (measured). Over the
mini-sweep curves (re-run for this corollary) the best common $\theta_s$ at $m = 8$ is $0.9$ and costs
the worst coefficient, $e$, $7.6\times$ (value) or $15\times$ (derivative) over its own optimum (an
earlier run gave $5.9\times$ and $11\times$: the figure moves by $\sim30\%$ between runs); at D12's
four terms ($m = 4$) the best common $\theta_s$ is $0.11$–$0.12$ and costs $33\times$ (value: $a$
$33\times$, $e$ $29\times$) and $97\times$ (derivative: $e$ $97\times$, $a$ $79\times$). $\square$

**Remark CO.12 (the jump at a switch).** At $\theta_s$ the arms differ by at most $E_x + E_s$, since
each is within its own error of $f$: up to $2E_\times$ at the crossing, not $E_\times$. The
truncation part of $E_s$ has a known sign (CO.8), the rounding parts do not. Measured at the
predicted crossing ($m = 8$; $b, c, d, e$; 6000 samples within $\pm0.5\%$): the maximum jump is
$9.6, 127, 40, 135\,u$ against $\max(E_x, E_s) = 6.6, 73.5, 26.3, 97.7\,u$: $1.4$–$1.7\times$.
`NUMERICS.md` §4 says the jump is at most the recorded max error of the coefficient; a
`branch_continuity_*` tolerance has to be the sum of the two arms' recorded errors (open item).

**Checked (CO.10–CO.12):** the mini-sweep above (mpmath 1.3.0 200 digits, `mp.diff` for derivatives,
the Dual and Horner arms of CO.9, CO.14; `f32` emulated exactly by rounding every operation to
binary32, `sin`/`cos` correctly rounded); prediction and measurement compared as stated; the jump
with $6000$ samples per row. The explicit crossing formulas of CO.10 evaluated with the CO.6, CO.14
constants reproduce the predicted columns of the table. For CO.11 and CO.18 the `f64` mini-sweep was
re-run (glibc; 30 samples per cell of $1/30$ decade over $[10^{-6}, \pi - 0.1]$, $m = 4, 8$;
references from the definitions at 120 digits, derivatives as the symbolic $\mathrm d/\mathrm dz$ of
the definitions): its per-coefficient optima agree with the table within $20\%$ in $\theta_s$ and $1.5\times$
in error. Scripts not committed. **Permanent:** none for the magnitudes (a sanity
check, by design); the sweep CSV (`xtask thresholds`); planned `branch_continuity_*` (`PHASE3.md` §9).

## 6. Derivatives

The kernel takes $z$ (`exp_coeffs(θ²)`), so `Dual<S, 1>` seeded at $z$ differentiates in $z$;
$\mathrm d/\mathrm d\theta = 2\theta\,\mathrm d/\mathrm dz$ and
$\partial/\partial\varphi_i = 2\varphi_i\,\mathrm d/\mathrm dz$ (zero at $\varphi = 0$: correct, the
coefficients are even). `PHASE1.md` §4.3, §6 do not say in which variable `coeff_*` derivatives are
taken (open item).

**Proposition CO.13 (derivatives and the series arm).**

- (a) $\sigma_m' = -\tau_{m+2}$. In the catalogue: $b_z = -e$; $a_z = d - \frac b2$;
  $d_z = \frac{b - 4d}{2z}$; $k_z = \frac{\cos\frac\theta2 - 2k}{4z}$;
  $c_z = \sum_{j\ge1}j\lvert B_{2j+2}\rvert z^{j-1}/(2j+2)! > 0$;
  $r_z = \frac1z\big(\frac w{z + w^2} - \frac r2\big)$ ($z \to n^2$, $w$ fixed);
  $(\cos\frac\theta2)_z = -\frac k4$.
  $e_z = -\sum_i(-1)^i(i+1)(i+2)z^i/(2i+7)!$. Only $b_z = -e$, $a_z = d - b/2$ and
  $(\cos\frac\theta2)_z = -k/4$ are combinations of catalogue members with no division by $z$; $d_z$, $k_z$, $r_z$ divide by $z$ (a cancellation),
  $c_z$ and $e_z$ are new functions.
- (b) Dual Horner in $z$ is the recurrence $d_j = z\,d_{j+1} + p_{j+1}$: it computes *exactly the
  derivative of the polynomial* $P_m'(z) = \sum_{1\le j<m}js_jz^{j-1}$. Its error is the truncation
  of the derivative series, $\lvert R'_m\rvert \le m\lvert s_m\rvert z^{m-1}$ for
  $z < \frac m{m+1}Z_m$ (alternating; for $c$, $R'_m \ge$ the first omitted term), plus rounding
  $\le 2.1u$ measured. Relative to $\lvert s_1\rvert$ that is $m\rho'_mz^{m-1}$: **larger than the
  value's truncation error $\rho_mz^m$ by $m\lvert s_0/s_1\rvert/z$**, with
  $\lvert s_0/s_1\rvert = 24, 12, 20, 60, 30, 21$ for $k, a, b, c, d, e$.

*Proof.* (a) differentiate term by term and shift the index; the closed forms are CO.4 and the
quotient rule, e.g. $b_z = -e$ because $\tau_5 = e$; $k_z$ from $k = \sin(\theta/2)/\theta$,
$\theta_z = 1/2\theta$. (b) the dual product rule on the Horner step; the tail of the derivative
series is alternating with ratio $\frac{(j+1)}{j}\frac{s_{j+1}z}{s_j}$. $\square$

**Observation CO.14 (through the exact arm the exponent is $p + 2$).** `Dual` differentiates the
*formula*, and the relative error of the derivative it returns was measured to be
$\hat C'u\,\theta^{-(p+2)}$: $p' = 2, 2, 4, 4, 4, 6$ for $k, a, b, c, d, e$; and $p' = 0$ for
$\cos\frac\theta2$ (error $\le 1.9u$: $-\sin\frac\theta2/4\theta$ has no subtraction). The exponents are
explained by the outline below and the constants by counting roundings in one worked case; no bound is
proved.

*Outline.* $b, d, e, k$ have the form $N/D$, $D \sim z^{\ell}$.
$f_z = N_z/D - \ell N/(zD)$: the two terms are of size $\lvert f\rvert/z$ while
$f_z \sim \lvert f\rvert$, a cancellation by $1/z$; and $N_z$ carries the relative error of the sum
$N$ (its terms are the derivatives of the terms of $N$, all $\sim1/z$ times larger, and its value is
$\sim N/z$: the same $\kappa_N$). The factors multiply: $\kappa_{f_z} \sim \kappa_N/z$. For $k$, $a$
($\kappa_N = 1$) only the quotient-rule cancellation remains. For $c$, the difference of two terms
$\sim1/z$ has derivative a difference of two terms $\sim1/z^2$ of value $\frac1{720}$.

*Worked case $b$.* $b = N/D$, $N = \theta - \sin\theta$, $D = \theta^3$. The `Dual` rules are `sqrt`
and `sin_cos` of `PHASE2.md` §3 and the textbook product $a_vb_d + a_db_v$ and quotient
$q = a_v/b_v$, $q_d = (a_d - q\,b_d)/b_v$ (`PHASE2.md` §3 states neither), so
$b_z = (N_z - \hat bD_z)/D$ with $N_z = (1 - \cos\theta)/2\theta \approx \theta/4$ and
$D_z = 3\theta/2$: $\hat bD_z \approx \theta/4$ and the difference is $-\theta^3/120$, a cancellation
by $60/\theta^2$. The error of $\hat b$ ($\le 6u/\theta^2$, CO.6, attained) enters as
$6u\theta^{-2}\cdot\theta/4$, that is $180u/\theta^4$ of the result (measured with only $\hat b$
inexact: $179$); the roundings inside $N_z$, mainly $\cos\theta$, add $\approx85$ (measured alone):
a first-order sum of $265$ against $222$ measured (the two do not peak together). For $c$ the same
count is $\theta^{-4}$ against $\frac1{720}$ times the roundings, $720\nu$ with $\nu = 6$ as in its value
($\kappa_{\rm bd} = 72 = 6\cdot12$): $4300$ against $4450$ measured, a consistency check, not a
derivation.

Measured maxima of $\hat C' = \text{error}/u \cdot \theta^{p+2}$ over six decades:
$33, 36, 228, 4450, 2400, 19900$ (medians $7, 7, 63, 900, 470, 4600$), constant to $\pm20\%$ from
$10^{-6}$ to $1$ (a re-run for the fixes above: $33, 33, 222, 4480, 2600, 21000$). The rule of division
of `Dual` changes the constants, not the exponents.

*Numerical consequences.* (i) The derivative term dominates the objective of 0004 item 1: at the
optimum with $m = 8$ the derivative error is $11, 32, 115, 8.8\text{e}3, 630, 2.2\text{e}3\,u$
against value errors $1.2, 2.5, 3.4, 42, 13, 49\,u$ ($k$, $a$: the value has no cancellation at all,
and the derivative alone sets their switch). (ii) At the best switch with $m \le 8$,
$\mathrm d/\mathrm dz$ through the shipped code errs by $10^2$–$10^4u$ for $b, c, d, e$; a tolerance
in `jacobians_match_dual_*` (`PHASE3.md` §9) must be per stratum, not a common multiple of $u$.
(iii) A derivative in the series arm does not lose a digit to cancellation, only to truncation
(CO.13(b)).

**Checked (CO.13–CO.14):** sympy: (a) as exact series to $z^{15}$ and symbolically ($\partial_z$ of
the definitions); mpmath `mp.diff` at 200 digits agrees with $b_z = -e$, $a_z = d - b/2$, $d_z$ to
$10^{-149}$; (b): 1600 points ($k, a, b, d, e$; $m = 1..8$), 0 violations, plus 240 for $c$;
(CO.14): `f64` `Dual` class (rules as in the worked case) through the exact arms of §3, 1200
log-uniform $\theta$ per decade over $[10^{-6}, 1]$, against `mp.diff` of the definitions; the re-run
against the symbolic derivative of the definitions at 120 digits; the split of the worked case by
making one source of error exact at a time (20000 samples over $[10^{-3}, 0.03]$).
**Permanent:** the sweep CSV's derivative column (0004 item 1); planned `jacobians_match_dual_*`.

## 7. Safe argument, $\theta = 0$, underflow, $\pi$

**Proposition CO.15.**

- (a) *$z = 0$.* The exact arms are $0/0$ ($k, a, b, d, e$) or $\infty - \infty$ ($c$):
  $\mathtt{NaN}$. The values $s_0 \ne 0$ and the derivatives $s_1$ are finite (CO.2, CO.3), so the
  series arm is finite on $[0, z_s)$, derivative included. The exact arm is evaluated at the safe
  argument $z \mapsto \mathtt{select}(z < z_s, 1, z)$ (0003 item 3): $\theta = 1$ is in the domain
  of all six arms, and the derivative of the constant $1$ is $0$, so an unselected arm is finite in
  value and derivative.
- (b) *Underflow.* $z = \varphi\cdot\varphi$ is subnormal for
  $\lVert\varphi\rVert < 1.5\times10^{-154}$ and $0$ below $1.6\times10^{-162}$ (`f32`:
  $1.1\times10^{-19}$, $2.6\times10^{-23}$). Every coefficient is analytic at $0$ with $s_0 \ne 0$,
  so an absolute error of the size of the subnormal spacing changes it by
  $\lvert s_1\rvert\,\delta z \le 10^{-323}$: the series arm returns $s_0$ (and $s_1$) to full
  relative accuracy for any $z$ below the switch, subnormal or zero. Nothing in the catalogue
  divides by $\theta$, and the consumers form $k\varphi$, $aW$, $bW^2$, never $\varphi/\theta$. The
  exact arm would fail long before: $b = 0$ below $2.6\times10^{-8}$, and $\theta^3$ leaves the
  normal range below $2.8\times10^{-103}$.
- (c) *Overflow.* $z = \infty$ for $\lVert\varphi\rVert > 1.3\times10^{154}$ (`f64`; `f32`:
  $1.8\times10^{19}$), and the exact arm returns $\mathtt{NaN}$ ($\sin\infty$). `NUMERICS.md` §12
  has no row for the domain of `Exp`; `PHASE2.md` §2 promises a finite result for finite in-domain
  input (open item).
- (d) *$c$ at $\pi$.* $1 + \cos\theta$ and $\sin\theta$ both vanish: the definition form is $0/0$,
  removable (CO.7(b)). The $\cot$ form is regular: $c(\pi) = 1/\pi^2$ and
  $c'(\pi) = -2/\pi^3 + 1/(4\pi) = 0.0151$ in $\theta$: no special case, no branch. For
  $\mathrm{fl}(\pi) = \pi - 1.2\times10^{-16}$ the exact $c$ is within $0.2u$ of $1/\pi^2$ and both
  computed forms within $2u$. The series arm is not an option: at $z = \pi^2$ it converges with
  ratio $\frac14$ and eight terms leave $\approx10^{-5}$.
- (e) *$c$ at $2\pi$.* $c$ has a pole; the rounding of $\theta$ itself (relative $u$) is amplified
  by $\theta c'/c \approx 2\pi/(2\pi - \theta)$, i.e. $\log_{10}\frac{2\pi}{2\pi - \theta}$ digits:
  the conditioning of $J^{-1}$ (LG.16), not of the formula, and the reason for the domain
  $\theta < 2\pi$ of §12.

**Checked (CO.15):** `f64`: the series arm at $z = 0$ and $10^{-310}$ returns $s_0$ exactly for all
six; $b$ by its definition returns $0$ at $\theta \le 2\times10^{-8}$; $c$ near $\pi$: 400 samples
per $\delta$ as in CO.7; $c'(\pi)$ by `mp.diff` against the formula; magnitudes of $z$ from the
`f64`/`f32` limits. **Permanent:** planned, corpus strata `theta:exact0`, `theta:subnormal`,
`theta:pi-1e-k` and the `nonfinite` column (`PHASE1.md` §4.4, §5); the seeded defect "`sqrt` of
$\theta^2$ without the safe argument, under `Dual`" (§10), implemented as `seeded:k-sqrt-unsafe` and
the both-arms mask tests of EA.18.

## 8. The $r$ kernel

**Proposition CO.16.** $r(n^2, w) = 2\,\mathrm{atan2}(n, w)/n$, $n > 0$, any real $w$; `Log` is
$\varphi = r\,q_{\mathrm v}$ (`NUMERICS.md` §3.2), $q_{\mathrm v}$ the vector part.

- (a) $r(\lambda^2n^2, \lambda w) = r/\lambda$ ($\lambda > 0$): exactly the scale invariance of
  $\mathrm{Log}$. For a unit $q$ with $w \ge 0$, $r = \theta/\sin\frac\theta2 = 1/k \in [2, \pi]$
  (CO.4(f)).
- (b) *Series.* CO.3(c): $r = \frac2w\sum(-s)^j/(2j+1)$, $s = n^2/w^2$, for $w > 0$ and $s < 1$; for
  a unit $q$, $s = \tan^2\frac\theta2 < 1 \iff \theta < \frac\pi2$. Truncation CO.8(d): the slowest
  of the seven.
- (c) *Sign of $w$.* For $w < 0$, $\mathrm{atan2}(n, w) = \pi - \arctan(n/\lvert w\rvert)$ and
  $r \to 2\pi/n$, not $2/w$: at $n^2 = 10^{-6}$, $w = -1$ the exact value is $6281.19$ and the
  series arm gives $-2.0000$. `Log` never sees $w < 0$ (the flip, §3.2), but `S2Chart::local` (§8)
  has $w = n\cdot m$ of either sign, so the mask of `log_ratio` must include $w > 0$ (and §8's
  $\alpha/\lVert n\times m\rVert$ is $r/2$, not $r$).
- (d) *$w \to 0$* ($\theta \to \pi$). $s \to \infty$: the series arm is never selected. The exact
  arm has no cancellation, $\mathrm{atan2}(n, w) \to \frac\pi2$, $r \to \pi/n$, finite at $w = 0$;
  measured error $\le 0.9u$ for $\frac\pi2 - \frac\theta2 \in [10^{-15}, 10^{-1}]$. The hazard is
  the *series* arm: $s = n^2/w^2 = \infty$ at $w = 0$, and a lane mask evaluates both arms, so here
  it is the series arm that needs the safe argument ($w \to \mathtt{select}(\mathrm{small}, w, 1)$),
  which 0003 item 3 and §4 state only for the exact arm.
- (e) *Branch variable.* The series is in $s$, `NUMERICS.md` §4 names $n^2$; they agree only for
  $\lVert q\rVert \approx 1$ ($s = n^2/(1 - n^2)$). For $q = \lambda(\cos1, \sin1\,\hat n)$,
  $\lambda = 10^{-100}$ ($\theta = 2$): $n^2 = 7\times10^{-201}$ passes any $n^2$ threshold,
  $s = 2.43 > 1$, and the 8-term $\sum(-s)^j/(2j+1)$ returns $-21.7$ for $0.642$. §12 restricts
  `log` to unit $q$, so this bites outside the domain only, but it is the reason "exactly
  scale-invariant" holds for the exact arm alone; a mask on $s$ (`n2 < s_sw * w * w`, division-free)
  has no such gap.
- (f) *Value and derivative.* The exact arm's value error is $\le 1.9u$ (measured; a quotient of
  accurate factors, $p = 0$). Its derivative $r_z$ through `Dual` errs by $\hat C'u/s$,
  $\hat C' \approx 6$ (max $6.3$, median $1.2$), i.e. $p + 2 = 2$ in $s$. With the series derivative
  $\frac2{w^3}\sum_{j\ge1}(-1)^jj\,s^{j-1}/(2j+1)$ (first term $-\frac13$), truncated at $m$ terms,
  relative error $\frac{3m}{2m+1}s^{m-1}$, CO.10 gives
  $s_\times = \big(\hat C'u(2m+1)/(3m)\big)^{1/m}$, $E_\times = \hat C'u/s_\times$: $m = 8$:
  $s_\times \approx 0.012$ ($\theta \approx 0.22$), $E_\times \approx 500u$ (measured $0.010$,
  $330u$); $m = 4$: $1.5\times10^{-4}$, $4\times10^4u$ (measured $1.3\times10^{-4}$,
  $2.9\times10^4u$). The derivative of $r$ has the smallest switch point of the seven at $m = 8$
  ($\theta \approx 0.22$ against $0.77$–$1.46$ for $k, \dots, e$, CO.10), because its series
  converges only like $s^m$; its error there ($494u$) is in the middle of theirs ($15$ to
  $1.2\times10^4u$).

**Checked (CO.16):** sympy: 8 terms of $2\arctan(n/w)/n$; $rk = 1$ and the 9-term $r$ series in
$\tan^2\frac\theta2$ against $\theta/\sin\frac\theta2$ ($O(\theta^{17})$). mpmath 200 digits, `f64`
`Dual`: 1500 unit-$q$ samples per decade for $n \in [10^{-8}, 0.5]$ (value and $r_z$); 300 near
$w = 0$; the mini-sweep of CO.10 for $s_\times$; (c), (e) are single evaluations. **Permanent:**
planned, corpus `coeff_r` and `so3_log`; none for (c)–(e), which are questions to the spec, not to
the code.

## 9. The reference precision budget

**Proposition CO.17.** The corpus computes each coefficient from its definition at `mp.dps` $= p$
digits (`PHASE1.md` §4.3). Evaluating $N = \sum t_i$ at $p$ digits loses $\log_{10}\kappa_{\rm sum}$
digits, $\kappa_{\rm sum} = \sum\lvert t_i\rvert/\lvert N\rvert$ (each term rounded to $10^{-p}$):

| | $a$ | $b$ | $c$ | $d$ | $e$ |
|---|---|---|---|---|---|
| $\kappa_{\rm sum}$ | $4\theta^{-2}$ | $12\theta^{-2}$ | $24\theta^{-2}$ | $48\theta^{-4}$ | $360\theta^{-4}$ |
| digits lost | $2L + 0.6$ | $2L + 1.1$ | $2L + 1.4$ | $4L + 1.7$ | $4L + 2.6$ |

so the worst, $e$, loses $4L + 2.6 \le 4L + 3$ (`PHASE1.md` §2 item 3: 51 at $\theta = 10^{-12}$,
leaving 69). The definition can serve while $p - (4L+3) \ge 30$ (the 30 printed digits; 17 suffice
to judge a binary64 result), and the 1% recompute at $150$ digits agrees to 40 digits only while
$p - (4L+3) \ge 40$, i.e. $\theta \gtrsim 10^{-19}$ for the value (measured $43.5$ digits at $10^{-19}$,
$39.2$ at $10^{-20}$). Below that, and at `theta:exact0` ($0/0$) and
`theta:subnormal` ($4L + 3 \approx 1243$ at $10^{-310}$), a definition at $120$ digits is unusable:
the reference there must be the defining series (no cancellation, $s_0$-dominated). `mp.diff` raises
its working precision itself, but what the derivative of $e$ keeps depends on the variable: in $z$
(the variable of CO.13) it keeps *more* digits than the value ($78.8$ against $71.1$ at $10^{-12}$), in
$\theta$ (the variable `PHASE1.md` §4.3 names) *fewer* ($66.9$). The derivative in $\theta$ therefore
needs its own line: it keeps 40 digits only down to $\theta \approx 3\times10^{-18}$ ($42.6$ digits at
$10^{-17}$, $37.0$ at $10^{-18}$), against $10^{-19}$ for the value. The `dps = 120` conclusion stands
because the corpus strata stop at $10^{-12}$.

**Checked (CO.17):** mpmath 1.3.0, the six definitions at $\theta = 10^{-L}$ (rounded to binary64),
$L \in \{3, 6, 9, 12, 15, 18, 20, 22, 24\}$: correct digits at `dps = 120` against `dps = 700`,
value and `mp.diff` derivative (in $\theta$ and in $z$). For $e$: $71.1$ digits at $L = 12$ (loss $48.9$,
bound $50.6$), $39.2$ at $20$, $31.0$ at $22$, $23.6$ at $24$; the derivative at $L = 12, 15, 17, 18, 20$:
$66.9, 51.8, 42.6, 37.0, 28.0$ in $\theta$ and $78.8, 66.7, 59.0, 54.9, 47.0$ in $z$ (the value: $71.1,
58.6, 50.9, 47.9, 39.2$; $17.5$ gives $39.5$ in $\theta$); for $a, b, c, d$ at $L = 12$:
$99.4, 96.9, 96.4, 74.3$ (losses $20.6, 23.1, 23.6, 45.7$, bounds $24.6, 25.1, 25.4, 49.7$).
`dps = 120` against `dps = 150` at $10^{-12}$: agree to $71$ (value), $67$ ($\mathrm d/\mathrm d\theta$),
$79$ ($\mathrm d/\mathrm dz$) digits, above the 40 required. **Permanent:** planned, the 1% `dps = 150` recompute of `PHASE1.md` §2 item
3.

## 10. Call-site groups

**Remark CO.18 (one mask per group).** `NUMERICS.md` §4 says a call site evaluates its coefficients
"in groups, inside one `branch`", and `PHASE3.md` §3 that each group is *one* `S::branch` over a
tuple while `generated.rs` holds a `Switch` per coefficient. Neither says whether a group takes one
mask (one $\theta_s$ shared by its members) or one mask per member. CO.10 shows the optima differ
inside a group ($m = 8$, derivative: $a$ $0.9$, $b$ $1.1$; $b$, $d$, $e$ $1.1$, $1.3$, $1.4$; $k$ $1.5$
against $\cos\frac\theta2$ $1.1$). The cost of sharing is the largest, over the members, of the error at
the best common $\theta_s$ divided by the member's own optimum (`f64`, mini-sweep curves of CO.11):

| group (`NUMERICS.md` §4) | $m = 8$: best $\theta_s$; cost (value, derivative) | $m = 4$: best $\theta_s$; cost |
|---|---|---|
| $(a, b)$ `jr_coeffs` | $0.96$, $0.89$; $1.1\times$, $1.7\times$ | $0.076$; $2.7\times$, $7.8\times$ |
| $(b, d, e)$ `q_coeffs` | $1.1$; $2.6\times$, $3.0\times$ | $0.18$, $0.19$; $7.1\times$, $9.6\times$ |
| $(k, \cos\frac\theta2)$ `exp_coeffs` | $1.4$, $0.96$; $1.0\times$, $1.5\times$ | $0.065$, $0.028$; $1.0\times$, $10\times$ |
| $(b, d)$ `gamma2_coeffs` | $1.1$, $0.96$; $1.6\times$, $1.7\times$ | $0.15$; $1.4\times$, $2.4\times$ |

So at $m = 8$ a shared mask costs at most $\sim3\times$ over per-coefficient switches, and at $m = 4$
up to $\sim10\times$; a one-mask-per-group reading is affordable at the series lengths CO.10 favours
and not at D12's. Which reading is intended is an open item in the index. The figures are cell
maxima of a $30$-sample sweep and carry the $\sim30\%$ run-to-run spread of CO.11.

**Checked (CO.18):** the mini-sweep re-run of CO.10–CO.12; the best common $\theta_s$ minimizes, on the
cell grid ($1/30$ decade), the largest ratio over the members. Scripts not committed.
**Permanent:** none for the magnitudes; the sweep CSV (`PHASE1.md` §6) decides per group once the
reading is fixed.
