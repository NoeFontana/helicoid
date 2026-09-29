# Error analysis, conditioning, forward-mode AD and the reference corpus

> Non-normative companion to [`NUMERICS.md`](../NUMERICS.md) §2.1, §3.2, §11, §12, to
> [`PHASE1.md`](../PHASE1.md) §2, §4, §5, §8, [`PHASE2.md`](../PHASE2.md) §3, to
> [`0003`](../decisions/0003-the-scalar-that-cannot-say-less-than.md),
> [`0006`](../decisions/0006-the-instrument-comes-first.md) and to `PROJECT.md` D5–D8, D16.
> **`NUMERICS.md` wins on any conflict; a disagreement is an open item in the [maths
> index](./index.md)**, which holds the notation and the `Checked:` convention. This page derives
> what the instrument measures and why its bars are built the way they are; it fixes no threshold,
> no series length and no formula. Statements about the generator, the harness and `Dual` describe
> what the specs require of code that does not exist yet: nothing here is a check.

$u$ is the unit roundoff of the precision under discussion (as in
[`coefficients.md`](./coefficients.md); `so3.md` writes it $\mathsf u$), $\mathcal F$ and
$\mathcal B$ (in units of $u$) are forward and backward error (EA.3), $\kappa$ a condition number,
$s(\theta) = 2\sin(\theta/2)/\theta$ is LG.16's, $\theta = \lVert\varphi\rVert$,
$\zeta = \pi - \theta$ near $\pi$ and $\eta = 2\pi - \theta$ near $2\pi$. Rotations follow
[`so3.md`](./so3.md) except that the vector part is $q_{\mathrm v}$, as in `coefficients.md` (`so3.md`
writes it $u$): $q = (w, q_{\mathrm v})$, $n = \lVert q_{\mathrm v}\rVert$, $\lVert q\rVert = 1$
unless stated. $b(\theta) = (\theta - \sin\theta)/\theta^3$ is the coefficient of `NUMERICS.md` §4
as a function of $z = \theta^2$, $P_m$ its $m$-term series, $s_j$ the coefficient of $z^j$ and $E_x$
the exact arm's error (CO.1, CO.8, CO.10); errors quoted for these are relative. Measured figures are `f64`
unless a precision is named.

## Results

| Label | Result | `NUMERICS.md` / spec | Implemented by |
|---|---|---|---|
| EA.1–EA.2 | the rounding model; $u = 2^{-53}, 2^{-24}$, `Dual<f64, N>` inherits it; $\hat\theta^2 = \theta^2(1+\varepsilon)$, $\lvert\varepsilon\rvert \le \gamma_3$; `libm` is not correctly rounded (`atan2` up to $1.96u$ measured, glibc $1.00u$) | §2.1 | every exact arm |
| EA.3–EA.5 | forward, backward error, $\kappa$; the metric of §11: chord $= 2\sin(\Delta/4)$, when its floors bind, what a norm-wise metric cannot see, $\mathcal B$ from the stored reference (valid on the reference's branch) and why Jacobians are judged by $\mathcal F$ | §11 | `xtask::conformance` |
| EA.6 | `Exp`: forward error $O(u)$ for $\theta \le \pi$, $\approx\theta u/2$ beyond | §3.1, §12 | `SO3::exp` |
| EA.7–EA.8 | `Log`: $\kappa_{\mathrm{nw}} = 1/\sin\frac\theta2$, $\kappa_{\mathrm{cw}} \le 2\cos\frac\theta2 + \sin\frac\theta2$; $s\,\lVert e\rVert \le \lVert\beta\rVert \le \lVert e\rVert$; what $\mathcal B$ adds near $\pi$ (branch invariance) | §3.2, §11 | `SO3::log` |
| EA.9 | `from_matrix`: the nearest-rotation quaternion is $\frac1{2\sqrt2}$-Lipschitz for every $\theta$; $\mathcal B_M^2 = d_M^2 + \lVert R_\star - R(\hat q)\rVert_F^2$ | §3.4, §11 | `SO3::from_matrix` |
| EA.10 | $J^{-1}$ of SO(3) near $2\pi$: $\kappa_2 = 1/s$ and $\kappa_\varphi \sim 2\pi/\eta$; the closed form errs by $\approx\chi\kappa_\varphi u$; for SE$_N$(3) only $\kappa_2 \ge 1/s$ | §3.5, §12 | `jr_inv`, `jl_inv` |
| EA.11–EA.13 | max and p99 per stratum, never a mean; domination and exact no-regress; non-finite counts; what D16 buys (reproducibility) and what it does not (correct rounding) | §11, D8, D16 | `xtask envelope`, `libm` on every target |
| EA.14–EA.15 | dual numbers; the rules of `PHASE2.md` §3 and the two they leave out | PHASE2 §3 | `Dual<S, N>` |
| EA.16 | the value-path theorem: `Dual` evaluates the same rounding sequence | 0003 item 5 | `Dual<S, N>`, `Real::branch` |
| EA.17 | the derivative through branches and series arms; second order by nesting | §4, §7 | `Dual` through `coeffs` |
| EA.18 | $\sqrt{\ }$ at $0$ and the safe argument | 0003 item 3, §4 | `coeffs`, `Real::select` |
| EA.19–EA.20 | exact inputs, 30 stored digits, the precision budget and the `dps = 150` recheck | PHASE1 §2, §4 | `conformance/generate` |
| EA.21–EA.23 | log-uniform decades; $z$-uniform, $\lambda$-uniform axes are uniform on $S^2$; splitmix64 and per-stratum streams | PHASE1 §4.4 | `conformance/generate` |

## 1. The floating-point model

**Definition EA.1.** Precision $p \in \{53, 24\}$, $u = 2^{-p}$, round to nearest even. For
$\circ \in \{+,-,\times,\div\}$ and for $\sqrt{\ }$ (all correctly rounded, IEEE 754):
$\mathrm{fl}(x \circ y) = (x \circ y)(1+\delta)$, $\lvert\delta\rvert \le u$, while the exact result
is in the normal range (CO.15 treats the rest). A library function
$g \in \{\sin, \cos, \mathrm{atan2}\}$ satisfies $\mathrm{fl}(g) = g\,(1+\delta)$,
$\lvert\delta\rvert \le \nu_g u$, with $\nu_g \ge 1$ a property of the implementation ($1$ iff
correctly rounded). $\gamma_k = ku/(1-ku)$ [Higham, ch. 3].

**Proposition EA.2.**

- (a) For $x \in [2^e, 2^{e+1})$ one ulp is $2^{e+1-p}$, so
  $\mathrm{ulp}(x)/\lvert x\rvert \in (u, 2u]$: an error of $m$ ulps is a relative error between
  $mu$ and $2mu$; correct rounding is $\le\frac12$ ulp $\le u$.
- (b) $\theta^2 = \varphi\cdot\varphi$ evaluated left to right has
  $\hat\theta^2 = \theta^2(1+\varepsilon)$, $\lvert\varepsilon\rvert \le \gamma_3$, **whatever the
  data whose squares stay normal**, i.e. every component is $0$ or of magnitude
  $\gtrsim1.5\times10^{-154}$ (`f64`; `f32`: $1.1\times10^{-19}$), since all terms are non-negative
  (no cancellation); and $\hat\theta = \sqrt{\hat\theta^2}$ has relative error
  $\le \frac{\gamma_3}{2} + u = 2.5u + O(u^2)$. Below that the error of $\hat\theta^2$ is only
  bounded by the subnormal spacing, which the coefficients do not feel (CO.15(b)) but which the
  relative bound does not cover.
- (c) `Dual<S, N>` stores $N + 1$ values of type $S$ and every one is produced by $S$-operations:
  the model applies to each rounding with the same $u$, so `Dual<f64, N>` has $u = 2^{-53}$ and
  `Dual<f32, N>` $2^{-24}$. It bounds each *operation*; a derivative component is the value of an
  expression that may cancel (CO.14), so it is not accurate to $O(u)$ in general.
- (d) The `libm` crate is not correctly rounded, and its `atan2` is not under one ulp. `libm`
  0.2.16, default features off: $\nu_{\sin} = 1.07$, $\nu_{\cos} = 1.09$,
  $\nu_{\mathrm{atan2}} = 1.96$ (maxima; $\nu$ is in units of $u$, the scale of $\mathcal F$), i.e.
  $0.76$, $0.73$ and $1.38$ ulp, not correctly rounded on $0.53\%$, $0.49\%$ and $8.5\%$ of the
  inputs; glibc 2.39: $\nu \le 1.00$ ($0.52$ ulp), wrong on $\le 0.07\%$; `sqrt` exact. These are
  samples, not the crate's guarantees.

*Proof.* (a) Spacing of $[2^e, 2^{e+1})$ over $\lvert x\rvert \in [2^e, 2^{e+1})$. (b) The three
squares carry at most $3, 3, 2$ roundings (product, then the additions above them), so each term is
within a factor $1\pm\gamma_3$ of its exact value [Higham, Lemma 3.1] and so is their sum (all
non-negative); then $(1+\varepsilon)^{1/2}(1+\delta) \le 1 + \gamma_3/2 + u + O(u^2)$. (c) By
construction. $\square$

**Checked (EA.2):** (b) $3\times10^5$ random $\varphi$ (components log-uniform in magnitude over six
decades, mpmath 40 digits): max error of $\mathrm{fl}(\varphi\cdot\varphi)$ $2.81u$ (bound $3u$), of
$\hat\theta$ $2.21u$ (bound $2.5u$); (d) `libm` 0.2.16 (the crate of `Cargo.lock`, built with
`default-features = false` as the workspace sets it) for x86_64,
$2\times10^5$ inputs (splitmix64 seed $12345$; `sin`, `cos`, `sincos`, `sqrt` at $x$ log-uniform on
$[10^{-8}, \pi]$; `atan2(y, x)` with $y \in [0,1)$, $x \in [0,1)$ and, on every second draw,
$x \in [0, 10^{-6})$) against mpmath 1.3.0 at 40 digits, glibc 2.39 on the same inputs; the
not-correctly-rounded shares are over the first $6\times10^4$. Script not committed. **Permanent:** none for the model; the corpus measures the
composite (`PHASE1.md` §4.3).

## 2. Errors, condition numbers and the metric of §11

**Definition EA.3.** For $y = f(x)$ and a computed $\hat y$ at the *exact* input $x$:

- *forward error*
  $\mathcal F(\hat y) = \dfrac{\lVert\hat y - y\rVert}{\max(\lVert y\rVert, \underline y)\,u}$, the
  metric of `NUMERICS.md` §11, with floor $\underline y \ge 0$;
- *backward error* for `Log`:
  $\mathcal B(\hat\varphi) = \lVert\mathrm{Exp}(\hat\varphi)\ominus_R X\rVert/u$. With
  $\beta = \mathrm{Log}(X^{-1}\mathrm{Exp}\,\hat\varphi)$, $\hat\varphi$ is *an* exact logarithm of
  $X\,\mathrm{Exp}(\beta)$ (not necessarily the canonical one) and $\lVert\beta\rVert = \mathcal Bu$
  is the size of the right perturbation of the input that makes it exact: the residual angle;
- *backward error* for `from_matrix`, at a matrix $M$ that is not exactly a rotation:
  $\mathcal B_M(\hat q) = \lVert M - R(\hat q/\lVert\hat q\rVert)\rVert_F/u$, absolute (not divided
  by $\lVert M\rVert_F \approx \sqrt3$). It is this page's reading of `NUMERICS.md` §11's "report
  backward error only", not a definition of §11 (index, open items). For one rotation error $\Delta$
  the angle, the quaternion chord (EA.4(a)) and the Frobenius norm of the matrix difference are,
  to first order, $\Delta : \Delta/2 : \sqrt2\Delta$ ($2\sqrt2\sin(\Delta/2)$ exactly), so
  $\mathcal B$, the quaternion $\mathcal F$ and $\mathcal B_M$ of the same error read
  $1 : \frac12 : \sqrt2$ times $\Delta/u$;
- *condition number*
  $\kappa(x) = \lim_{\varepsilon\to0}\sup\big(\lVert f(x{+}\Delta x) - f(x)\rVert/\lVert f(x)\rVert\big)/\varepsilon$
  over a stated class of perturbations of relative size $\varepsilon$: **norm-wise**
  ($\lVert\Delta x\rVert \le \varepsilon\lVert x\rVert$: what an input built by arithmetic carries)
  or **componentwise** ($\lvert\Delta x_i\rvert \le \varepsilon\lvert x_i\rvert$: what rounding an
  input carries).

To first order, relative forward error $\le \kappa\cdot$ relative backward error, in the matching
class.

*Jacobians, $\mathrm{Ad}$ and $V$ are judged by $\mathcal F$ alone.* Their outputs are analytic in
the tangent for $\theta < 2\pi$ (no branch, no cut) and have floors $\ge 1$ (EA.4(b)), so $\mathcal F$
is a relative error and nothing needs the branch invariance that motivates $\mathcal B$ for
`Log` (EA.8(b)). A backward error would ask for a $\varphi'$ with $J(\varphi') = \hat J$ exactly, but
the image of the tangent under $\varphi \mapsto J(\varphi)$ has the dimension of the tangent (3, or
$3 + 3N$), far below that of the matrices, and a computed $\hat J$ is off it in general: its
"backward error" would be a projection, not a residual. What ties $\mathcal F$ to the input is the
conditioning: bounded for $J$, $\kappa_\varphi \sim 2\pi/\eta$ for $J^{-1}$ (EA.10).

**Proposition EA.4 (the metric).**

- (a) *Chord.* For unit $q, q'$ with $q\cdot q' \ge 0$ and $\Delta$ the rotation angle of
  $R(q)^\top R(q')$, $\lVert q - q'\rVert = 2\sin(\Delta/4)$. So
  $\lVert\hat q - q\rVert \approx \Delta/2$: $1\,u$ in §11's quaternion metric (sign aligned) is a
  rotation error of $2u$.
- (b) *Floors.* $\lVert q\rVert = 1$; for $\theta \le \pi$,
  $\lVert J\rVert_F \ge \sqrt{1 + 2s^2} \ge 1.34$, $\lVert J^{-1}\rVert_F \ge \sqrt3$,
  $\lVert\mathrm{Ad}\rVert_F \ge \sqrt3$, $\lVert V\rho\rVert \ge \frac2\pi\lVert\rho\rVert$: a
  floor of $1$ (resp. the $\rho$-scale) never binds on these outputs. It binds on `Log`'s tangent,
  whose norm is $\theta$: with $\underline y = 1$ a total loss ($\hat\varphi = 0$) at
  $\theta = 10^{-12}$ reads $\theta/(\underline y u) = 9\times10^{3}$, with
  $\underline y \ll \theta$ it reads $9\times10^{15}$. `NUMERICS.md` §11 gives the floor of a
  tangent output no value; the strata `theta:1e-k` see relative errors only if it is
  $\lesssim\theta$. On `theta:subnormal` the model of EA.1 does not apply (CO.15(b)) and an
  output below the smallest normal number has fewer than $p$ significant bits: $\mathcal F$
  measures its quantization, a half-spacing rounding reading $\le1$ only with a floor at least the
  smallest normal number ($2^{-1022}$; `f32` $2^{-126}$) and without bound below it.
- (c) *What a norm-wise metric cannot see.* $\mathcal F$ certifies absolute accuracy
  $u\lVert y\rVert$ per component. In $\mathrm{Exp}\,\varphi = (\cos\frac\theta2, k\varphi)$ the
  vector part has size $\theta/2$: a relative error $\varepsilon$ of $k$ (or of $aW$ in
  $J = I - aW + \cdots$) shows as $\varepsilon\theta/2 < u$, invisible for $\varepsilon < 2u/\theta$
  ($2\times10^{-8}$ at $\theta = 10^{-8}$). This is harmless to a consumer (the object is right to
  $u$) and is why the scalar `coeff_*` ids, which are relative, exist.

*Proof.* (a) $q^*q'$ has scalar part $\cos\psi = q\cdot q' \ge 0$ ($\psi \le \pi/2$) and rotation
angle $\Delta = 2\psi$; the chord is $2\sin(\psi/2)$. (b) Singular values $(1, s, s)$ of LG.16(a),
$s \ge 2/\pi$ on $[0, \pi]$, and $(1, 1/s, 1/s)$ for the inverse; $V = J_l$ for $N = 1$; a matrix
has at least the norm of its diagonal block. (c) arithmetic. $\square$

**Checked (EA.4):** mpmath 1.3.0, 60 digits: (a) 40 random pairs, $\theta \le 3.14$,
definition-based relative rotation,
$\lvert\text{chord} - 2\sin(\Delta/4)\rvert \le 2.3\times10^{-60}$; (b) arithmetic and LG.16(a).
Scripts not committed. **Permanent:** planned, `xtask::conformance` metrics (`PHASE1.md` §5).

**Proposition EA.5 ($\mathcal B$ from the stored reference).** Let $\varphi$ be the stored reference,
$\hat\varphi$ the subject's output and $\psi$ the nearer to $\varphi$ of $\hat\varphi$ and
$\hat\varphi\,(1 - 2\pi/\lVert\hat\varphi\rVert)$ (two logarithms of one rotation, LG.2(c)); let
$e = \psi - \varphi$, formed with the extra digits of EA.19(b). Then
$\mathcal B\,u = \lVert J_r(\varphi)e\rVert\,(1 + O(\lVert e\rVert))$ **for $\lVert e\rVert \ll 1$**:
no evaluation of $\mathrm{Exp}$ at the subject's output is needed, and $J_r$ is needed to three
digits (`f64` suffices). The reduction to $\psi$ is not optional. With $e = \hat\varphi - \varphi$
the form holds only for a subject on the branch of $\varphi$; a subject that returns the other
logarithm ($e \approx \mp2\pi\hat n$: $\hat\theta > \pi$, or the other sign at $w = 0$) reads
$\lVert J_re\rVert \approx 2\pi$ where $\mathcal Bu$ is $O(u)$, and the branch invariance of
EA.8(b) is lost. If $\lVert e\rVert$ is not small for both representatives the form is a
first-order estimate of a large residual, and the value is the exact
$\lVert\mathrm{Log}(\mathrm{Exp}(\varphi)^{-1}\mathrm{Exp}\,\hat\varphi)\rVert$, which needs
$\mathrm{Exp}$ at the output (index, open items).

*Proof.* $\mathrm{Exp}\,\psi = \mathrm{Exp}\,\hat\varphi$, and
$\mathrm{Exp}(\varphi+e) = \mathrm{Exp}(\varphi)\mathrm{Exp}(J_re + O(e^2))$ (LG.8), so
$\beta = \mathrm{Log}(\mathrm{Exp}(\varphi)^{-1}\mathrm{Exp}\,\hat\varphi) = J_re + O(e^2)$; for
$\lVert e\rVert \sim 10^{-16}$ the remainder is $\sim10^{-32}$. $\square$

**Checked (EA.5):** mpmath 60 digits, $\theta \in \{10^{-4}, 1, 3, \pi - 10^{-6}\}$,
$\lVert e\rVert = 10^{-16}\theta$, exact $\mathcal B$ from `mp.expm` and the definition-based
$\mathrm{Log}$ against $\lVert J_re\rVert$: relative deviation $\le 4\times10^{-18}$. Other branch: axis
$(1,2,3)/\sqrt{14}$, $\theta \in \{3.1415, \pi - 10^{-6}, \pi\}$, $\hat\varphi$ the other logarithm
$\varphi(1 - 2\pi/\theta)$ plus noise $\le10^{-16}$ per component: exact $\mathcal Bu = 5.6\times10^{-17},
6.5\times10^{-17}, 1.1\times10^{-16}$; unreduced $\lVert J_re\rVert = 6.28319$ ($= 2\pi$) in all three;
the reduced form agrees with the exact value to $2\times10^{-18}$ (relative). At $\theta = \pi$ with
$\hat\varphi = -\varphi$ plus noise: exact $1.2\times10^{-16}$, unreduced $2\pi$, reduced
$1.2\times10^{-16}$. `PHASE1.md` §5 does not say how the harness evaluates $\mathcal B$ (index,
open items). Script not committed. **Permanent:** planned, corpus `so3_log`.

## 3. Conditioning

**Proposition EA.6 (`Exp`).** (a) For every $\theta$, $\lVert D\mathrm{Exp}\rVert$ is at most
$\frac12$ in the quaternion metric: the rotation moves by $\lVert J_re\rVert \le \lVert e\rVert$
(singular values $1, s, s$) and the quaternion by half of it. (b) A relative perturbation
$\lvert\varepsilon\rvert \le u$ of $\theta$ moves the output by $\theta u/2$ (direction
$\hat\varphi$, singular value $1$). (c) So a routine that forms $\hat\theta$ in `S` (relative error
up to $2.5u$, EA.2(b)) has a $\theta$-proportional forward error of up to $2.5\,\theta u/2$
(measured $1.3$–$1.6\times\theta u/2$ for $\theta \ge 10$) on top of an $O(u)$ part independent of
$\theta$: $3.3u$ at $\theta = 3$, linear beyond. This is **backward stability**, not a defect: the
$\theta$-proportional part is a relative perturbation of $\varphi$'s magnitude, and only an
extended-precision $\theta$ (or a reduction mod $2\pi$) does better. §12 has no `Exp` domain row
(index, open items), and the corpus strata stop at $\pi$.

**Checked (EA.6):** `f64` exact arm (`math.sin`/`cos`, $\theta$ from the dot product) against mpmath
50 digits at the exact double $\varphi$, 2000 random axes per $\theta$: max $\mathcal F$
$= 0.33, 0.74, 3.3, 7.2, 81, 6.4\times10^{3}$ for $\theta = 0.1, 1, 3, 10, 100, 10^4$ (yardstick
$\theta/2$: $0.05, 0.5, 1.5, 5, 50, 5\times10^{3}$). Script not committed. **Permanent:** planned,
corpus `so3_exp` (strata up to $\pi - 10^{-12}$).

**Proposition EA.7 (`Log`: conditioning).** Let $\lVert q\rVert = 1$, $w > 0$,
$\theta = 2\,\mathrm{atan2}(n, w) < \pi$.

- (a) *Absolute.* $D\mathrm{Log} = J_r^{-1}(\varphi)$ on rotation perturbations,
  $\lVert\cdot\rVert_2 = 1/s(\theta)$: $1 + \theta^2/24$ at $0$, $\pi/2$ at $\pi$.
- (b) *Norm-wise relative.* $\kappa_{\mathrm{nw}}(\theta) = 1/\sin(\theta/2)$: $\approx 2/\theta$
  near $0$, $1$ at $\pi$.
- (c) *Componentwise relative.* For $\lvert\Delta q_i\rvert \le \varepsilon\lvert q_i\rvert$:
  $\lVert\Delta\varphi\rVert/\theta \le \varepsilon\,(2\cos\frac\theta2 + \sin\frac\theta2) + O(\varepsilon^2) \le \sqrt5\,\varepsilon$,
  and $\to 2\varepsilon$ as $\theta \to 0$.
- (d) *The cut.* $\mathrm{Log}$ jumps by $2\pi\hat n$ across $w = 0$ (SO.5(b), (d)), at quaternion
  distance $w = \sin(\zeta/2) \approx \zeta/2$ from $q$. Input noise $\ge w$ changes the output by
  up to $2\pi$; below it (c) holds.

*Proof.* (a) LG.8, LG.16(a). (b) $\mathrm{Log}$ is scale-invariant (SO.5(c)), so only
$\Delta q \perp q$ matters; then $\Delta q = \frac12 q\,(0, \xi)$,
$\lVert\Delta q\rVert = \lVert\xi\rVert/2$, $\Delta\varphi = J_r^{-1}\xi$, and
$\kappa = \lVert J_r^{-1}\rVert\cdot 2/\theta = 2/(s\theta) = 1/\sin\frac\theta2$. (c)
$\xi = 2\,\mathrm{Im}(q^*\Delta q) = 2(w\,\Delta q_{\mathrm v} - \Delta w\,q_{\mathrm v} - q_{\mathrm v}\times\Delta q_{\mathrm v})$,
$\lVert\xi\rVert \le 2\varepsilon n(2w + n)$, and
$\lVert\Delta\varphi\rVert/\theta \le \lVert\xi\rVert/(s\theta) = \lVert\xi\rVert/(2n)$. The maximum
of $2\cos a + \sin a$ is $\sqrt5$; at $\theta \to 0$ take $\Delta w = \varepsilon w$,
$\Delta q_{\mathrm v} = -\varepsilon q_{\mathrm v}$. (d) The distance from $q$ to the great sphere $w = 0$ is $\arcsin w$.
$\square$

*Consequence.* A quaternion carrying **rounding** noise (componentwise) keeps $\mathrm{Log}$ to
$O(\varepsilon)$ relative at every $\theta$ (components in the normal range, EA.2(b)), which is what
the $\mathrm{atan2}$ form delivers
(SO.6(c)); one carrying **absolute** noise (a matrix, a product) loses $2/\theta$: the same
phenomenon as the $\arccos$ form (SO.6).

**Checked (EA.7):** mpmath 1.3.0, 60 digits. $\mathrm{Log}$ is the solution of
$\mathrm{Exp}(\varphi) = R(q)$ by Newton on
$\mathrm{vee}\,\mathrm{skew}(R^\top\mathrm{Exp}\,\varphi)$ with the series $J_r$ (a definition; the
starting guess does not matter); $D\mathrm{Log}/\partial q$ by central differences, $h = 10^{-25}$;
(b) the largest singular value on $q^\perp$ over $\theta$, and (c) the maximum over the 16 sign
vertices of $\lVert(D\mathrm{Log})(\pm\lvert q_i\rvert)\rVert/\theta$ (a convex function on a box:
attained at a vertex). $\theta = 10^{-6}, 10^{-3}, 0.5, 1, 2, 3, \pi - 10^{-3}, \pi - 10^{-6}$, one
random axis each: $\kappa_{\mathrm{nw}}$ equals $1/\sin\frac\theta2$ to $10$ digits ($2\times10^6$,
$2000.0$, $4.04$, $2.09$, $1.19$, $1.0025$, $1.0000001$, $1$);
$\kappa_{\mathrm{cw}} = 2.0, 2.0, 1.92, 1.68, 1.12, 0.83, 0.97, 0.73$ against the bound
$2.0, 2.0, 2.19, 2.23, 1.92, 1.14, 1.00, 1.00$. $\lVert J_r^{-1}\rVert_2 = 1/s$ to $12$ digits at
$\theta = 10^{-6}, 0.5, 2, \pi - 10^{-6}, 2\pi - 10^{-2}$. Scripts not committed. **Permanent:**
planned, corpus `so3_log`.

**Proposition EA.8 (`Log`: backward against forward, and what $\mathcal B$ adds near $\pi$).**
(a) and (c) take $e = \hat\varphi - \varphi$ small, $\hat\varphi$ on the branch of $\varphi$; (b) is
about outputs that may be on the other branch.

- (a) $s(\theta)\lVert e\rVert \le \lVert\beta\rVert \le \lVert e\rVert$ to first order, both
  attained ($e \parallel \varphi$: $1$; $e \perp \varphi$: $s$). Hence
  $\lVert e\rVert \le \lVert\beta\rVert/s(\theta) \le \frac\pi2\lVert\beta\rVert$ on
  $\theta \le \pi$: the absolute forward error and $\mathcal B$ are equivalent within $\pi/2$, and
  *near $\pi$ (where $\lVert\varphi\rVert \approx \pi$, so absolute is relative) nothing is lost by
  reporting $\mathcal B$*. For an **exact** quaternion with $w > 0$ the $\mathrm{atan2}$ form is on
  the branch of $\varphi$ (SO.5(a)), so on such a stratum the choice between $\mathcal F$ and
  $\mathcal B$ is one of convention, not of information.
- (b) *Branch invariance.* $\mathcal B$ is unchanged by replacing $\hat\varphi$ with another
  logarithm of the same rotation (it is computed from $\mathrm{Exp}\,\hat\varphi$; from the stored
  reference, after the reduction of EA.5); $\lVert\hat\varphi - \varphi\rVert$ changes by up to
  $2\pi$. It matters where the output can be on either side of the cut: (i) a subject, or an
  oracle, that returns the other logarithm near $w = 0$ ($\lVert\hat\varphi - \varphi\rVert
  \approx 2\pi$, $\lVert\beta\rVert = O(u)$); (ii) input carrying noise $\varepsilon \gtrsim w$
  (EA.7(d): a composite, a product, a matrix), where $\lVert\hat\varphi - \varphi\rVert$ is
  $\approx0$ or $\approx2\pi$ by luck; (iii) $w = 0$ itself, where both $\pm\pi\hat n$ are
  logarithms (SO.5(d)) with $\beta = 0$ and $\lVert\hat\varphi - \varphi\rVert \in \{0, 2\pi\}$.
  Case (iii) is
  what §11's sign-invariant $\min(\lVert\hat\varphi - \varphi\rVert, \lVert\hat\varphi + \varphi\rVert)$
  removes, **at $w = +0$ exactly** (the only zero §11 names); for $0 < w \lesssim$ noise, (i) and
  (ii) remain, and there $\mathcal F$ conflates rounding with branch selection, which the stratum
  `q:w0` and the proptest `so3_log_w_flip` test on their own.
- (c) *Absolute against relative.* $\mathcal B$ is an angle: for $\theta \ll 1$,
  $\lVert\beta\rVert \approx \lVert e\rVert = \theta\cdot(\text{relative error})$: it is blind by a
  factor $\theta$ to exactly the relative errors that $\mathcal F$ (with a small floor, EA.4(b))
  sees.

The split of `NUMERICS.md` §11, $\mathcal F$ at small $\theta$ and $\mathcal B$ primary near $\pi$,
is (c) and, for the cases of (b), the branch invariance; on exact input on the branch it is (a).

*Proof.* (a) EA.5 and the singular values $(1, s, s)$. (b) SO.5(d), EA.5. (c) (a) with
$\theta \ll 1$. $\square$

**Checked (EA.8):** (a), (c): mpmath 60 digits, $\theta = 10^{-6}, 0.5, 2, \pi - 10^{-6}, 2\pi - 10^{-2}$,
$\lVert e\rVert = 10^{-25}$, 22 directions (parallel, perpendicular, 20 random) through the
definition-based $\mathrm{Log}$: $\min\lVert\beta\rVert/\lVert e\rVert$ $= s$ (to 12 digits:
$1, 0.9896, 0.8415, 0.6366, 1.6\times10^{-3}$), max $1$; (b): EA.5's other-branch runs. Script not committed. **Permanent:**
planned, corpus `so3_log`, proptests `so3_log_w_flip`, `so3_log_at_w0_is_a_function_of_the_sign`.

**Proposition EA.9 (`from_matrix`).** Let $R \in \mathrm{SO}(3)$, $M = R + E$, $q_\star(M)$ the
quaternion of the nearest rotation (SO.13(c)), signs aligned.

- (a) $\lVert\Delta q_\star\rVert \le \lVert E\rVert_F/(2\sqrt2) + O(\lVert E\rVert^2)$, attained
  ($E = R\,[\omega]_\times$), **for every $\theta$ including $\pi$**: the map
  $M \mapsto \pm q_\star$ is uniformly well conditioned.
- (b) What is ill conditioned near $\pi$ is not (a): it is (i) the *orientation* of the axis
  relative to the angle (the sign of $w$, read from $R - R^\top = 2\sin\theta[\hat n]_\times$, of
  size $\approx\zeta$): the distance from $R$ to the cut $\theta = \pi$ is $\zeta$ in rotation angle
  and $\sqrt2\,\zeta$ in $\lVert\cdot\rVert_F$, so noise $\gtrsim\sqrt2\zeta$ flips
  $\mathrm{Log}\circ\mathrm{from\_matrix}$ (EA.7(d)); (ii) the trace pivot,
  $\lvert\Delta\hat q_i\rvert \propto 1/w^2$ (SO.12(d)), an *algorithm* choice the largest pivot
  removes.
- (c) *Backward error* ($\mathcal B_M$ of EA.3, Frobenius, $R(\hat q)$ at the normalized $\hat q$).
  With
  $d_M = \mathrm{dist}_F(M, \mathrm{SO}(3)) = \sqrt{\lVert M\rVert_F^2 + 3 - 2\lambda_{\max}(K(M))}$
  (SO.13(c)) and $\hat q$ within $O(\lVert E\rVert)$ of $q_\star$,
  $\mathcal B_M^2u^2 = \lVert M - R(\hat q)\rVert_F^2 = d_M^2 + \lVert R_\star - R(\hat q)\rVert_F^2 + O(\cdot^3)$,
  and $\lVert R_\star - R(\hat q)\rVert_F = 2\sqrt2\,\lVert\hat q - q_\star\rVert$. So
  $\mathcal B_M \ge d_M/u$ **for every routine**; for a matrix rounded from a rotation
  $d_M \le \lVert E\rVert_F$ is a few $u$ (Frobenius units: $\sqrt2$ times the rotation angle it
  stands for, EA.3). The algorithm's share is $\sqrt{\mathcal B_M^2u^2 - d_M^2}$.
- (d) *Why backward only.* Under noise $\varepsilon$ the "correct" quaternion is not determined
  below $O(d_M)$: Shepperd's extraction and the polar factor differ by a few $\varepsilon$ (SO.12(d)
  bounds each component of the extraction's error by $\frac74\varepsilon$; (a) bounds the polar
  factor's by $1.06\,\varepsilon$, with $\lVert E\rVert_F \le 3\varepsilon$) by *definition*, not by
  rounding (SO.13(b)); (c) is the measure that needs no choice of reference, and the excess over
  $d_M$ isolates the routine's own error.

*Proof.* (a) $R^\top M = I + A$, $A = R^\top E$; the polar factor of $I + A$ is
$I + \mathrm{skew}(A) + O(A^2)$, so $\Delta R_\star = R\,S$, $S = \mathrm{skew}(A)$,
$\lVert S\rVert_F \le \lVert A\rVert_F = \lVert E\rVert_F$ (equality iff $A$ is skew),
$\lVert S\rVert_F = \sqrt2\lVert\omega\rVert$ and $\Delta q = \frac12 q\,(0, \omega)$. (b) SO.12,
EA.7(d); $\arccos\lvert q\cdot q'\rvert$ over the great sphere $q'_w = 0$ is minimal at
$q' = (0, q_{\mathrm v}/n)$, at angle $\zeta/2$; the rotation matrices differ by $2\sqrt2\sin(\zeta/2) \approx\sqrt2\zeta$. (c)
$M - R_\star = R_\star(P - I)$ with $P$ symmetric and
$R_\star - R(\hat q) = -R_\star[\omega']_\times + O(\omega'^2)$;
$\langle R_\star(P - I), R_\star[\omega']_\times\rangle_F = \mathrm{tr}((P-I)^\top[\omega']_\times) = 0$
(symmetric against skew). (d) SO.12(d) for the pivot extraction against $R$; (a) for the polar
factor, $\lVert E\rVert_F \le 3\varepsilon$. $\square$

**Checked (EA.9):** mpmath 1.3.0, 60 digits. (a) $q_\star$ as the top eigenvector of $K(M)$
(`mp.eigsy`), the $4\times9$ Jacobian by central differences ($h = 10^{-25}$): largest singular
value $0.353553390593 = 1/(2\sqrt2)$ at $\theta = 10^{-6}, 1, 3, \pi - 10^{-6}$ and exactly $\pi$.
(c) 18 random $M = R + E$, $\lVert E\rVert = 10^{-4}, 10^{-6}, 10^{-8}$ (entries),
$\hat q = q_\star$ plus a perturbation of the same size: relative residual of the Pythagoras
identity $5.7\times10^{-5}, 4.5\times10^{-7}, 4.6\times10^{-9}$ (first order,
$\propto\lVert E\rVert$); $d_M^2$ from $K$ against $\lVert M - R_\star\rVert_F^2$ from the SVD:
$\le10^{-44}$. (b) $\zeta = 10^{-2}, 10^{-6}, 10^{-12}$: quaternion distance to $w = 0$ equals
$\zeta/2$, no random point of the great sphere closer,
$\lVert R - R'\rVert_F/(\sqrt2\zeta) = 0.999996, 1, 1$. Scripts not committed. **Permanent:**
planned, corpus `so3_from_matrix` (backward error).

**Proposition EA.10 ($J^{-1}$ near $2\pi$).** Let $\theta = 2\pi - \eta$, $0 < \eta \ll 1$. (a)–(c)
are for SO(3), $N = 0$, where the singular values are $(1, s, s)$ (LG.16(a)); (d) says what
remains for SE$_N$(3).

- (a) $\kappa_2(J^{-1}) = \kappa_2(J) = 1/s = \theta/(2\sin\frac\eta2) \sim 2\pi/\eta$ (LG.16(a));
  the two singular directions are those perpendicular to $\hat n$, on which $\mathrm{Exp}$ is not
  injective at $2\pi$.
- (b) The relative sensitivity to a relative perturbation of $\varphi$,
  $\kappa_\varphi = \lVert\mathrm d J^{-1}/\mathrm d\ln\theta\rVert_2/\lVert J^{-1}\rVert_2$, is
  $2\pi/\eta\,(1 + O(\eta))$.
- (c) A closed form that forms $\hat\theta$ in `S` (relative error $\le2.5u$, EA.2(b)) therefore
  errs by $\approx\chi\,\kappa_\varphi u$, $\chi \lesssim 1.4$ measured: $6\times10^{3}u$ at
  $\eta = 10^{-3}$, and to keep $d$ digits $\eta \gtrsim 2\pi u\,10^{d}$.
- (d) $J$ itself is entire in $z = \theta^2$ and its $\mathrm d/\mathrm d\theta$ is bounded: no such
  loss. For $N \ge 1$ LG.16(b) gives only $\kappa_2(J) \ge 1/s$ and $\lVert J^{-1}\rVert_2 \ge 1/s$:
  $\kappa_2 = 1/s$ is *not* claimed (EA.20 measures $10^{8.3}$ to $10^{10.2}$ against $1/s = 74.6$ at
  $\theta = 6.2$, $\lVert\rho\rVert = 10^4$), and the $(i,0)$ block of $J^{-1}$ scales like
  $\lVert\rho_i\rVert/s^2$ (LG.16(d)).

*Proof.* (a) LG.16(a). (b) With $J^{-1} = I + \frac12W + cW^2$,
$c = \theta^{-2} - \cot(\theta/2)/(2\theta)$: $\cot\frac\theta2 = -\cot\frac\eta2 \approx -2/\eta$,
so $c \approx (\theta\eta)^{-1}$, $\mathrm dc/\mathrm d\ln\theta \approx \eta^{-2}$, and
$\lVert W^2\rVert_2 = \theta^2$: the derivative has size $\theta^2\eta^{-2}$ against
$\lVert J^{-1}\rVert \approx c\theta^2 = \theta/\eta$. (c) Rounding $\hat\theta$ is a relative
perturbation of $\varphi$'s magnitude; multiply by (b). (d) CO.2. $\square$

*Consequence.* §12 promises a finite value for $\theta < 2\pi$ and there is no corpus stratum beyond
$\pi$. A stratum at $\eta = 10^{-k}$ would measure $\kappa_\varphi u$, the conditioning of the
input, for every routine of this form; a domination bar there compares extended-precision $\theta$
handling, not the closed form.

**Checked (EA.10):** mpmath 60 digits, $J = \sum(-W)^n/(n+1)!$ and its `mp.inverse`,
$\eta = 10^{-1}, \dots, 10^{-5}$: $\kappa_\varphi$ measured $62.86, 628.32, 6283.19, 62831.9$
against $2\pi/\eta = 62.83, 628.32, 6283.19, 62831.9$; the `f64` closed form ($c$ from the $\cot$
arm of §4, $\hat\theta$ from the dot product), 60 random axes per $\eta$, at the exact double
$\varphi$: max $\mathcal F$ $= 83, 570, 8.5\times10^{3}, 4.5\times10^{4}, 5.8\times10^{5}$ against
$2\pi/\eta$ (ratios $1.33, 0.91, 1.36, 0.71, 0.93$). Script not committed. **Permanent:** none (no
stratum beyond $\pi$).

## 4. Statistics, domination, no-regress and determinism

**Proposition EA.11 (max and p99 per stratum, never a mean).** Let a stratum have errors
$\mathcal E_1, \dots, \mathcal E_n \ge 0$ (each an $\mathcal F$ or $\mathcal B$ of one record).

- (a) $\bar{\mathcal E} \le \max\mathcal E \le n\,\bar{\mathcal E}$: a bar $T$ on the mean admits a
  max of $nT$ ($64T$ on a 64-sample stratum), and a bar on the max implies one on the mean. Only the
  max bounds a single call.
- (b) A mean is an average over the *sampler*. For $\mathcal E(\theta) = C\theta^{-p}$ (CO.6) on a
  log-uniform decade, mean$/$sup $= (1 - 10^{-p})/(p\ln10)$: $0.215, 0.109, 0.072$ for
  $p = 2, 4, 6$; uniform-in-$\theta$ sampling changes it. The max of the corpus is a fixed,
  reproducible number (EA.13) that lower-bounds the stratum's supremum (EA.21).
- (c) A defect erring by $D$ on a fraction $f$ of the inputs reads $fD$ in the mean (a missing
  $w < 0$ flip: $f = \frac12$; a defect confined to `theta:exact0` or `q:w0`: $f \ll 1$) and $D$ in
  the max; a mean across strata dilutes it by the number of strata (twelve `theta:1e-k` decades).
  `PROJECT.md` §6 lists "one number aggregated over strata" as a smell.
- (d) *p99.* With $n = 64$ the nearest-rank p99 is rank $\lceil0.99\cdot64\rceil = 64$, **the max**;
  interpolation (position $0.99\cdot63 = 62.37$) gives $x_{(63)} + 0.37(x_{(64)} - x_{(63)})$. It
  differs from the max only from $n \ge 100$ (`theta:dense`, $800$ samples: the 9th largest). Its
  role is the distance $\max/\mathrm{p99}$, a flag for an isolated worst input; it is never a bar.
  The method is not specified (index, open items).
- (e) *Non-finite outputs.* A max over records is meaningful only if every $\mathcal E_i$ is a
  number: an IEEE `max` or a `>` fold silently drops a NaN, so a subject that returned NaN would pass
  a max bar. `NUMERICS.md` §11 counts non-finite outputs separately and any non-zero count fails
  (the `nonfinite` column of `PHASE1.md` §5): the bars below apply to strata whose count is $0$.

**Definition EA.12 (bars).** *Domination:* per (function, stratum, precision),
$\max_{\mathrm{helicoid}} \le \min_{o}\max_{o}$, ties passing. *No-regress:*
$\max_{\mathrm{helicoid}} \le \max_{\mathrm{baseline}}$, compared exactly.

**Proposition EA.13 (properties of the bars; what D16 buys).**

- (a) Domination compares maxima, not distributions: errors $(1, \dots, 1, 10)$ dominate
  $(0.1, \dots, 0.1, 12)$ though the second is $10\times$ better on $63$ of $64$ inputs. It is a
  worst-case guarantee, which is what a consumer's tail needs, and says nothing about the mean or
  the median. No-regress bounds the supremum, not each record: one record may worsen from $3$ to $9$
  while the max falls from $10$ to $9$.
- (b) With a tolerance $\tau$ a baseline sequence may drift, $B_k \le (1+\tau)^kB_0$ ($1\%$ over
  $100$ PRs is $2.7\times$); the exact bar is a ratchet, $B_k \le B_0$. Exactness is possible only
  if the max is a function of the source and the corpus alone: if outputs vary with the target by
  $d$ ulps the max varies by $\sim d$ and a tolerance is unavoidable.
- (c) **D16 buys reproducibility.** $+,-,\times,\div,\sqrt{\ }$ are correctly rounded and rustc
  neither fuses nor reassociates without an explicit intrinsic (`mul_add` is banned); a `libm`
  function is a program over these and integer operations; so every output is a function of the
  input bits alone (NaN payloads apart; a NaN is a failure, EA.11(e)), identical on x86_64, aarch64
  and wasm32, **with `libm`'s default features off** as `Cargo.toml` sets them (its optional `arch`
  fast paths would break bit identity). Measured, with them off: $2\times10^5$ inputs through
  `libm::{sin, cos, sincos, atan2, sqrt}` (six outputs each) built for x86_64 and for
  `wasm32-wasip1` (run under wasmtime): $0$ of $1.2\times10^6$ outputs differ, digests equal
  (aarch64 was not available: it is `just determinism`'s). The same inputs through `f64::sin`,
  `cos`, `atan2` on x86_64 (glibc 2.39) differ from `libm`'s bits in $0.53\%$, $0.48\%$ and
  $8.4\%$ of them; on `wasm32-wasip1` `std` and `libm` agree on all of them, so that difference is
  glibc's, not the target's.
- (d) **D16 does not buy accuracy.** EA.2(d), in the scale of $\mathcal F$ (units of $u$): `libm` has
  $\nu = 1.07, 1.09$ on `sin`, `cos` and $1.96$ on `atan2`, glibc $\le 1.00$. The constants measured
  on this repository are statements about the implementation they ran on: the transcendental part of
  a `libm` constant can be up to $\sim1.1\times$ (`sin`, `cos`) or $\sim2\times$ (`atan2`) glibc's.
  The envelope compares `helicoid` on `libm` with oracles on their platform's libm (containers
  pinned by `PHASE1.md` §7), so a `Log` stratum whose error is dominated by its final `atan2` may
  start up to $\approx1u$ behind an oracle with a nearly correctly rounded one: a domination failure
  there can be the library's and not the algorithm's, while no-regress, which compares `helicoid`
  with itself, is unaffected. Nor does D16 cover a change of the `libm` version (new bits: a new
  baseline) or `-C target-cpu`/fast-math (banned). The claim is *checked* by `just determinism`
  (Phase 6, planned), not by this page.

**Checked (EA.11–EA.13):** (b) numerical integration and $2\times10^6$ Monte Carlo draws; (d)
arithmetic on ranks. (c): `libm` 0.2.16 with `default-features = false`, a scratch program
(splitmix64 seed $12345$, `libm::pow` for the log-uniform draw) on x86_64-linux and `wasm32-wasip1`
under wasmtime 49.0.0 (Python `wasmtime` package): FNV-1a digests of $2\times10^5\times6$ outputs
equal ($\mathtt{b3394f65dc097b21}$; the same digest with default features on, on this host);
`std` against `libm` on the same inputs: on x86_64 (glibc 2.39) $1069$, $966$, $16879$ of
$2\times10^5$ differing bitwise (sin, cos, atan2), on `wasm32-wasip1` $0$, $0$, $0$. (An earlier version drew $x$ with `f64::powf` and differed in $104$ of
$2\times10^5$ inputs before any `libm` call: a `std` call in generator code is exactly what D16
forbids.) Scripts not committed. **Permanent:** planned, `xtask envelope` (bars), `just determinism`
(Phase 6).

## 5. Forward-mode AD

**Definition EA.14.** $\mathbb D_N(S) = S[\epsilon_1, \dots, \epsilon_N]/(\epsilon_i\epsilon_j)$:
elements $v + \sum_id_i\epsilon_i$ (`Dual { v, d }`); for $f$ differentiable at $v$,
$f(v + d\cdot\epsilon) = f(v) + f'(v)\,d\cdot\epsilon$ (every higher Taylor term contains some
$\epsilon_i\epsilon_j = 0$), and for several arguments $f(a + \epsilon b, \dots)$ adds
$\sum_k\partial_kf\,b_k$.

**Proposition EA.15 (the rules).** With $\mathrm{sgn}(\pm0) = +1$:

| operation | value | derivative part |
|---|---|---|
| $a + b$, $a - b$ | $a_v \pm b_v$ | $a_d \pm b_d$ |
| $a\,b$ | $a_vb_v$ | $a_vb_d + a_db_v$ |
| $a/b$ | $q = a_v/b_v$ | $(a_d - q\,b_d)/b_v$ |
| $\sqrt a$ | $\sqrt{a_v}$ | $a_d/(2\sqrt{a_v})$ |
| `sin_cos` | $\sin a_v,\ \cos a_v$ | $a_d\cos a_v,\ -a_d\sin a_v$ |
| $\mathrm{atan2}(y, x)$ | $\mathrm{atan2}(y_v, x_v)$ | $(x_vy_d - y_vx_d)/(x_v^2 + y_v^2)$ |
| $\lvert a\rvert$ | $\lvert a_v\rvert$ | $\mathrm{sgn}(a_v)\,a_d$ |
| $\mathrm{copysign}(x, s)$ | $\mathrm{copysign}(x_v, s_v)$ | $\mathrm{sgn}(x_v)\,\mathrm{sgn}(s_v)\,x_d$ |

The last five rows are `PHASE2.md` §3's. `PHASE2.md` §3 states no rule for $\times$ and $\div$; the
two forms of the quotient rule, $(a_d - q\,b_d)/b_v$ and $(a_db_v - a_vb_d)/b_v^2$, are equal in
exact arithmetic and differ in the last bit on $57\%$ of random inputs, so the choice moves every
derivative the sweep of `0004` item 1 measures (index, open items). At the kinks ($\lvert a\rvert$,
$\mathrm{copysign}$, $\mathrm{select}$) the rule returns *a* value, not a derivative:
$\lvert x\rvert$ has none at $0$, and $\mathrm{sgn}(\pm0) = +1$ is a choice.

*Proof.* Table rows 1–3: expand $(a_v + a_d\epsilon)(b_v + b_d\epsilon)$ and solve
$(q_v + q_d\epsilon)(b_v + b_d\epsilon) = a_v + a_d\epsilon$. Others: EA.14 with
$\frac{\mathrm d}{\mathrm dv}\sqrt v = \frac1{2\sqrt v}$, $(\sin)' = \cos$, $(\cos)' = -\sin$,
$\partial_y\mathrm{atan2} = x/(x^2 + y^2)$, $\partial_x\mathrm{atan2} = -y/(x^2 + y^2)$;
$\lvert\cdot\rvert$ and $\mathrm{copysign}(x, s) = \mathrm{sgn}(s)\lvert x\rvert$ are piecewise
linear in $x$ and constant in $s$. $\square$

**Checked (EA.15):** mpmath 1.3.0, a `Dual` over `mpf` (50 digits) implementing the table, against
`mp.diff`: $\sqrt{\ }$, $\sin$, $\cos$, $\mathrm{atan2}$, the exact arm of $b$ and its 4-term
polynomial: first derivatives to $\le 8\times10^{-49}$; the quotient forms on $2\times10^5$ random
`f64` draws. Script not committed. **Permanent:** planned, `dual_matches_mpmath_derivative` on
corpus ids `real_*` (`PHASE2.md` §8).

**Theorem EA.16 (value path).** Let $F$ be generic code over `S: Real` built from `Real` operations,
`select`, `branch` and masks. Run it with $S = T$ and with $S = \texttt{Dual<T, N>}$, and assume

- (H1) each `Dual` operation computes its value part with the *same* $T$-operation on the value
  parts (`sin_cos` calls $T$'s `sin_cos`, not `sin` and `cos` separately);
- (H2) masks are computed from value parts only;
- (H3) no derivative part is read into a value part (no `value_f64` in generic code outside tests
  and `debug_assert!`), and `select` moves value and derivative parts together.

Then both runs take the same branches, and every intermediate value part of the `Dual` run equals
the plain intermediate **bit for bit**.

*Proof.* Let $\mathrm{pr}(v + d\cdot\epsilon) = v$. (H1) says
$\mathrm{pr}(\mathrm{op}_D(a, b)) = \mathrm{op}_T(\mathrm{pr}\,a, \mathrm{pr}\,b)$ bitwise for every
operation, (H2) that $\mathrm{mask}_D(a, b) = \mathrm{mask}_T(\mathrm{pr}\,a, \mathrm{pr}\,b)$, (H3)
that no operation feeds $d$ into $v$. Induct on the evaluation trace: equal masks give the same
control path, and equal operands the same results. $\square$

*Consequences.* (i) Nested towers project to the innermost scalar ($\mathrm{pr}\circ\mathrm{pr}$),
and if each derivative rule is itself written in $S$-operations (EA.17) the first-order parts of
$\mathbb D_M(\mathbb D_N(S))$ equal those of $\mathbb D_N(S)$ bit for bit. (ii) A `Dual` run
differentiates the *arm the value run took*: "the derivative of the shipped code". (iii) **The
theorem is silent on derivative parts:** their rounding is not tied to the value's (CO.14), and they
can be non-finite while every value is finite (EA.18).

**Checked (EA.16):** a `Dual` (rules of EA.15, IEEE `f64` scalars, textbook $\times$ and $\div$)
through the exact arm, the 4-term series arm and the `branch` at $z < 0.25$ of $b$: $10^5$
log-uniform $z \in [10^{-6}, 5]$, value part against the plain evaluation: $0$ bit mismatches for
`Dual<f64, 3>` and for `Dual<Dual<f64, 3>, 3>`; $10^5$ draws of `atan2`, `sin`, `cos` on
`Dual<f64, 2>`: $0$. Separately, on the $2\times10^5$ inputs of EA.13, `libm::sincos` equals
`libm::sin`, `libm::cos` bitwise, so (H1)'s "same call" is not yet observable; nothing guarantees
it. Script not committed. **Permanent:** planned, `dual_value_is_plain_value` (`PHASE2.md` §8, §10).

**Proposition EA.17 (what the derivative is; second order).**

- (a) Let the masks partition the input space into cells on which $F = F_\sigma$, a composition of
  elementary functions. At an interior point the `Dual` run returns $\mathrm DF_\sigma$ exactly
  (chain rule, exact arithmetic): the derivative of the *computed* function.
- (b) It is the derivative of the mathematical function $f$ iff $F_\sigma' = f'$: true on an exact
  arm; on a series arm $F_\sigma = P_m$ and the result is $P_m'$, off by the truncation of the
  *derivative* series, larger than the value's (both relative) by $m\lvert s_0/s_1\rvert/z$
  (CO.13(b)). On a mask
  boundary the run returns the one-sided derivative of the arm the value mask selects; at a switch
  point the derivative is discontinuous by at most $\lvert P_m' - f'\rvert + \lvert E_x'\rvert$ (the
  derivative of the exact arm's error).
- (c) $\mathbb D_M(\mathbb D_N(S))$ carries $f$, $\partial_if$, $\partial_jf$ and
  $\partial_i\partial_jf$:
  $f(v + \epsilon_1g + \epsilon_2h + \epsilon_1\epsilon_2m) = f + f'(\epsilon_1g + \epsilon_2h) + \epsilon_1\epsilon_2(f'm + f''gh)$,
  so the seeds $g = h = 1$, $m = 0$ give $f''$ in the $\epsilon_1\epsilon_2$ part. This requires
  every rule of EA.15 to be written in operations of $S$ (its derivative part must itself be
  differentiable): `sqrt`'s $a_d/(2\sqrt{a_v})$ uses the inner tower's `sqrt` and division. The
  masked rules ($\lvert\cdot\rvert$, `copysign`) have second derivative $0$ almost everywhere.
- (d) In floating point each derivative order costs a factor $\theta^{-2}$ through an exact arm: for
  $b$ ($p = 2$, CO.6) value, first and second $z$-derivative err by $\sim6u\theta^{-2}$,
  $\sim10^2u\theta^{-4}$ (CO.14) and $\lesssim10^4u\theta^{-6}$ measured, all digits gone at
  $\theta \approx 10^{-2}$. The sweep of `0004` measures the first order only; a Hessian through
  `Dual<Dual<S, N>, N>` needs a larger switch.

*Proof.* (a) The chain rule on each cell. (b) Definition of the arms (CO.13). (c) Expand; the
$\epsilon_1\epsilon_2$ coefficient of $f(v + \epsilon_1g + \epsilon_2h + \epsilon_1\epsilon_2m)$
collects $f'm$ and $f''gh$. (d) measured. $\square$

**Checked (EA.17):** mpmath 50 digits, a `Dual<Dual<mpf>>` against `mp.diff` of order 2 (Hessian of
$\mathrm{atan2}$ in both arguments, all four entries): $\sqrt{\ }$, $\sin$, $\cos$,
$\mathrm{atan2}$, the exact arm of $b$ and its 4-term polynomial, at $z = 0.37$ or
$(y, x) = (0.7, -0.4)$: relative deviation $\le 10^{-46}$; the first-order parts of the nested tower
equal the first-order `Dual`'s exactly (asserted, `mpf`). (b): 4-term series of $b$ at $z = 0.2$:
relative value error $2.4\times10^{-10}$, relative derivative error $9.7\times10^{-8}$ (ratio
$400 = 4\cdot20/0.2$, CO.13(b)). (d): `f64` through `Dual<Dual<f64, 1>, 1>`, 300 samples per $\theta$: max error $/u$
(value, $\mathrm d/\mathrm dz$, $\mathrm d^2/\mathrm dz^2$) $= (3.7, 186, 1.0\times10^4)$ at
$\theta \approx 1$, $(346, 1.3\times10^6, 7.0\times10^9)$ at $0.1$,
$(4.2\times10^4, 1.3\times10^{10}, 6.6\times10^{15})$ at $0.01$. Scripts not committed.
**Permanent:** planned, `dual_matches_mpmath_derivative`; the nested case is promised by `PHASE2.md`
§3 and has no named test.

**Proposition EA.18 ($\sqrt{\ }$ at $0$ and the safe argument).**

- (a) $\sqrt v$ has an infinite derivative at $0$: the rule of EA.15 gives $a_d/0$, that is
  $\pm\infty$ for $a_d \ne 0$ and $\mathtt{NaN}$ for $a_d = 0$, which is the case for
  $z = \varphi\cdot\varphi$ at $\varphi = 0$ ($\mathrm dz = 2\varphi = 0$).
- (b) **Value finite does not imply derivative finite.** With a shared $\theta = \sqrt z$ at
  $z = 0$, $(\sin, \cos)(\theta/2)$ have finite values $(0, 1)$ and derivatives
  $(\infty, \mathtt{NaN})$. The exact arm of $b$ is $0/0$ in value and derivative. The composite
  $k(\theta^2)$ is analytic (CO.2); `Dual` differentiates its *factorization* through $\sqrt{\ }$,
  which is not.
- (c) *The safe argument* $\theta = \sqrt{\mathrm{select}(\mathrm{small}, 1, z)}$ evaluates the
  exact arm at a constant $1$ whenever the series arm is selected: value **and** derivative of the
  *exact* arm are finite at every input (the selected constant has derivative $0$). The *series*
  arm is a polynomial of degree $m - 1$ in $z$, finite where $z^{m-1}$ is, and a lane mask evaluates
  it at every input: for $z \gtrsim 10^{308/(m-1)}$ ($10^{44}$ at $m = 8$, i.e. $\theta \gtrsim 10^{22}$,
  inside the `Exp` domain of CO.15(c)) it overflows, and an arithmetic blend
  $0\cdot P_m + 1\cdot(\text{exact})$ is $0\cdot\infty = \mathtt{NaN}$. So "any blend is finite"
  holds for in-domain input ($\theta < 2\pi$), and beyond it only if the series arm also gets a
  safe argument (as CO.16(d) states for $r$). A scalar `bool` mask evaluates only the chosen arm, but
  not a subexpression hoisted out of it (the shared $\theta$ above); a lane mask evaluates both
  arms; reverse mode multiplies $0$ by the arm's derivative. The series arm at $z = 0$ returns
  $P_m'(0) = s_1$ (CO.15); in $\varphi$,
  $\partial_{\varphi_i} = 2\varphi_i\,\mathrm d/\mathrm dz = 0$ (CO.13).
- (d) `abs` at $0$ returns $+a_d$, $\mathrm{atan2}(0, 0)$ has derivative $0/0$: the same class,
  guarded the same way (CO.16(d) for $r$).

**Checked (EA.18):** `f64` (IEEE via numpy scalars), the rules of EA.15: exact arm of $b$ at $z = 0$
with $\mathrm dz = 1$: $(\mathtt{NaN}, \mathtt{NaN})$; shared $\theta$: values $(0.0, 1.0)$,
derivatives $(\infty, \mathtt{NaN})$; $z = \varphi\cdot\varphi$ at $\varphi = 0$:
$\mathrm d\sqrt z = (\mathtt{NaN}, \mathtt{NaN}, \mathtt{NaN})$; `branch` with the safe argument at
$z = 0$: $(1/6, -1/120)$ exactly; the same code with the series arm masked out:
$(\mathtt{NaN}, \mathtt{NaN})$; $\lvert x\rvert$ rule at $\pm0$: $+1$; $\mathrm{atan2}(0, 0)$:
$(\mathtt{NaN}, \mathtt{NaN})$; the 8-term series of $b$ at $z = 10^{45}$ is finite and at $10^{50}$
is $-\infty$, where $0\cdot P_8 + 1\cdot0.5$ is $\mathtt{NaN}$. Script not committed. **Permanent:** planned, the seeded defect
"`sqrt` of $\theta^2$ without the safe argument, under `Dual`" (`PHASE1.md` §10; `nonfinite > 0` in
`theta:exact0`).

## 6. The reference corpus

**Proposition EA.19 (exact inputs; why 30 digits suffice).**

- (a) An input is a binary64 written with `float.hex()` and the reference is computed at exactly
  that rational, so $\mathcal F$ measures the computation's error on the datum given: there is **no
  input-rounding term** and no dependence on $\kappa$ of the input. The value $10^{-8}$ is not a
  double; the corpus stores its nearest double and the stratum label is nominal. For `theta:pi-1e-k`
  rounding the components moves $\lVert\varphi\rVert$ by up to $\approx\pi u$ ($2.2\times10^{-16}$
  measured): the realized $\pi - \theta$ deviates from $10^{-k}$ by relative
  $2.6\times10^{-15}, 2.8\times10^{-13}, 2.4\times10^{-10}, 2.7\times10^{-7}, 2.2\times10^{-4}$ for
  $k = 1, 3, 6, 9, 12$ (measured maxima).
- (b) 30 significant digits round with relative error $\le\frac12\cdot10^{-29} = 5\times10^{-30}$
  componentwise, hence norm-wise, so a reference stored to 30 digits changes any $\mathcal F$ by at
  most $5\times10^{-30}/u = 4.5\times10^{-14}$ (`f64`) or $8.4\times10^{-23}$ (`f32`) in units of
  $u$: thirteen orders of magnitude below what a bar resolves. **This holds only if the harness
  forms $\hat y - y$ with the extra digits:** parsing $y$ to `f64` first quantizes $\mathcal F$ to
  an ulp ($\sim1u$) and an exact no-regress bar (EA.12) would be blind to anything below it. A
  binary64 needs $17$ digits to be identified and $\mathcal F$ is resolved to $10^{-3}$ of $u$ from
  $\sim20$ digits: $30$ leaves $10$ guard digits (double-double, $\sim32$, is enough).
- (c) *`f32`.* The corpus inputs are binary64. An `f32` subject that rounds them to `f32` evaluates
  at a different point, and $\mathcal F$ then includes $\kappa\,u_{32}$ of input rounding (EA.6(b),
  EA.7(c): up to $\theta/2$ for `Exp`, $\le\sqrt5\,\varepsilon$ for `Log` from a quaternion).
  `PHASE1.md` §4.4 does not say whether `f32` strata use `f32`-exact inputs (index, open items).

**Checked (EA.19):** (a) mpmath 40 digits, 400 axes per $k$, components rounded to `f64`; (b)
arithmetic. Scripts not committed. **Permanent:** planned, `just corpus-check` (`PHASE1.md` §2 item
1).

**Proposition EA.20 (precision budget; the `dps = 150` recheck).**

- (a) A definition evaluated at $p$ digits keeps $p - \ell$ digits,
  $\ell = \log_{10}\kappa_{\mathrm{alg}}$:
  $\kappa_{\mathrm{sum}} = \sum\lvert t_i\rvert/\lvert\sum t_i\rvert$ for a series (CO.17: the
  coefficient called $e$ in `NUMERICS.md` §4, not the tangent error of EA.5, loses $4L + 3$ at
  $\theta = 10^{-L}$), $\kappa_2(A)$ for a solve. For
  $J_{r,l}(\tau) = \sum(\mp\mathrm{ad}_\tau)^n/(n+1)!$ the terms grow only polynomially in $\rho$
  (linear in $\rho$, LG.16(d)). $\kappa_2(J) = (\lVert\rho\rVert/2)^2$ holds at $\varphi = 0$ only
  (LG.16(c)); for $\theta > 0$ the rotation block multiplies it, the direction of $\rho$ matters, and
  the general lower bound is $\kappa_2 \ge 1/s(\theta)$ (LG.16(b)). Measured ($N = 1$,
  $\lVert\rho\rVert = 10^4$, where $(\lVert\rho\rVert/2)^2$ is $\ell = 7.40$): $\ell = 7.35$ to $7.55$
  for $\theta \le 3.18$ and $8.30$ to $10.17$ at $\theta = 6.2$ ($1/s = 74.6$), the extremes at
  $\rho$ perpendicular and parallel to the axis. An unguarded solve at $120$ digits then keeps
  $\ge 109$; measured, mpmath's `inverse` (which carries guard digits) agrees between
  $\mathtt{dps} = 120$ and $150$ to $\ge117$ digits in all cases below.
- (b) *The recheck.* Recompute a sample at $p' = 150 > p = 120$ and demand agreement to $d = 40$
  digits. For one algorithm at two precisions,
  $\lvert\hat y_p - \hat y_{p'}\rvert \approx \lvert\hat y_p - y_{\mathrm{alg}}\rvert$ up to the
  smaller error, so agreement to $d$ digits certifies $\hat y_p$ to $\sim d$ digits as a value *of
  that algorithm*: $40 = 30 + 10$, and the printed $30$ digits are wrong only when the true value
  lies within $10^{-40}$ of a rounding boundary (probability $\sim10^{-10}$ per string, and then by
  one unit in the last digit).
- (c) **It cannot see an error that does not depend on $p$**: a truncated series, a wrong formula, a
  wrong branch. A $12$-term series of $b$ at $\theta = 3$ agrees between $120$ and $150$ digits to
  $120.8$ digits and is wrong by a relative $2.4\times10^{-16}$; `mp.logm` returns a complex, non-principal
  result for $\theta \ge 3.03$ at every precision (index, open items). The recheck guards against
  *precision loss* only; the *definition* is guarded by cross-checks between independent definitions
  (matrix exponential and quaternion series), the `Checked:` lines of these pages and the seeded
  defects (`PHASE1.md` §10).

**Checked (EA.20):** mpmath 1.3.0, $\tau = (\varphi; \rho)$, $N = 1$, axis $z$, $J$ by its series,
$J^{-1}$ by `mp.inverse`, $\mathrm{Exp}$ by `mp.expm`, agreement digits between
$\mathtt{dps} = 120$ and $150$ for $(\theta, \lVert\rho\rVert) = (10^{-8}, 1), (10^{-8}, 10^4),
(1, 10^4), (3.18, 10^4), (6.2, 10^4)$ with $\rho$ perpendicular / parallel to the axis: $J$
$121.3/121.4$, $120.8/121.5$, $120.5/120.5$, $120.5/120.7$, $120.5/119.1$; $J^{-1}$
$121.4/121.6$, $120.7/121.0$, $120.5/120.2$, $120.0/119.9$, $118.1/117.8$; $\mathrm{Exp}$ $\ge121$;
$\log_{10}\kappa_2(J) = 0.21$, $7.40$, $7.39/7.41$, $7.35/7.55$, $8.30/10.17$. Over axes $z$ and
$(1,1,1)/\sqrt3$ with $\rho$ perpendicular, parallel and four random directions
($\lVert\rho\rVert = 10^4$): $\ell \in [7.40, 7.40]$, $[7.39, 7.41]$, $[7.35, 7.55]$,
$[8.30, 10.17]$ for $\theta = 10^{-8}, 1, 3.18, 6.2$. The truncated series was summed at both
precisions. Scripts not committed. **Permanent:** planned, the 1% `dps = 150` recompute
(`PHASE1.md` §2 item 3).

**Proposition EA.21 (log-uniform decades).** Read `theta:1e-k` as $\theta \in [10^{-k}, 10^{-k+1})$
(`theta:1e0` is $[1, \pi - 0.1)$), $\log_{10}\theta$ uniform. The error curves of the catalogue are
$C\theta^{-p}$ (CO.6), straight lines in log-log, so uniform $\ln\theta$ gives equal resolution per
e-fold; uniform $\theta$ would put $5.6\%$ of the samples below $1.5\times$ the lower end of the
decade instead of $17.6\%$, where the maximum is. With $n$ samples the maximum error captures, in
expectation, $n\int_0^1(1-t)^{n-1}10^{-pt}\,\mathrm dt \approx n/(n + p\ln10)$ of the supremum over
the decade ($t$ the log-position of the smallest draw): $0.93, 0.88, 0.82$ for $p = 2, 4, 6$ at
$n = 64$; $0.98, 0.96, 0.94$ at $n = 200$. The corpus max is a lower bound on the stratum's supremum
by that margin; `theta:dense` ($200$ per decade over $[10^{-4}, 1]$) is the finer look.

*Proof.* $\Pr(t_{\min} > t) = (1 - t)^n$, so $t_{\min}$ has density $n(1-t)^{n-1}$ and the captured
fraction is $\mathbb E\,10^{-pt_{\min}}$. $\square$

**Checked (EA.21):** the integral against $2\times10^5$ Monte Carlo trials of $n = 64, 200$ draws:
$0.9338/0.9336$, $0.8757/0.8757$, $0.8243/0.8242$ ($n = 64$); the approximation $n/(n + p\ln10)$ is
within $0.002$. Script not committed. **Permanent:** none.

**Proposition EA.22 (axes uniform on $S^2$).** If $z \sim U[-1, 1]$ and $\lambda \sim U[0, 2\pi)$
are independent, $\hat n = (\sqrt{1 - z^2}\cos\lambda,\ \sqrt{1 - z^2}\sin\lambda,\ z)$ is uniform on
$S^2$. Hence $\hat n\cdot a \sim U[-1, 1]$ for every fixed unit $a$.

*Proof (Archimedes).* With $r = \sqrt{1 - z^2}$, $\hat n_z = (-z\cos\lambda/r, -z\sin\lambda/r, 1)$ and
$\hat n_\lambda = (-r\sin\lambda, r\cos\lambda, 0)$ (partial derivatives);
$\hat n_z\times\hat n_\lambda = -(r\cos\lambda, r\sin\lambda, z) = -\hat n$, of norm $1$. So the area element is
$\mathrm dA = \mathrm dz\,\mathrm d\lambda$, the product of the two uniform measures (total $4\pi$).
Rotating $a$ to $e_z$ preserves area, so $\hat n\cdot a$ has the law of $z$. $\square$

Choosing the polar angle $\vartheta$ uniformly is not uniform on the sphere
($\mathrm dA = \sin\vartheta\,\mathrm d\vartheta\,\mathrm d\lambda$): the two $10^\circ$ caps around
the poles hold $11.1\%$ of the draws instead of $1.5\%$. On the $53$-bit grid $z$ has spacing
$2^{-52}$, so an axis within $\approx2\times10^{-8}$ rad of a pole is drawn with probability
$\sim2^{-53}$; the axis is then rounded to binary64, which is harmless (EA.19(a)). The model
assumed for the generator is this discretized one, $z \in 2^{-52}\mathbb Z \cap [-1, 1)$ and
$\lambda \in 2\pi\,2^{-53}\mathbb Z$, and the independence of $z$ and $\lambda$ is an idealization:
successive splitmix64 outputs are consecutive values of a permutation (EA.23(a)), of the $2^{128}$
pairs only $2^{64}$ occur. It is checked below, not a property of the generator.

**Checked (EA.22):** $2\times10^6$ draws from splitmix64 streams (EA.23): Kolmogorov–Smirnov
statistic of $\hat n\cdot a$ against $U[-1,1]$ over $25$ random $a$: max $7.1\times10^{-4}$ (1% critical
value $1.15\times10^{-3}$); $\mathbb E[\hat n] \le 1.9\times10^{-4}$,
$\max\lvert\mathbb E[\hat n\hat n^\top] - I/3\rvert = 2.5\times10^{-4}$; the polar-angle sampler: KS $0.105$
for $a = e_z$, $11.1\%$ within $10^\circ$ of a pole (either). Script not committed. **Permanent:**
none (the generator's own tests are `PHASE1.md`'s to specify).

**Proposition EA.23 (splitmix64 and per-stratum streams).** splitmix64 is
$a \leftarrow a + \mathsf G$ ($\mathsf G = \mathtt{0x9E3779B97F4A7C15}$, odd), output
$\mu(a) = z_2 \oplus (z_2 \gg 31)$ with $z_2 = (z_1 \oplus (z_1 \gg 27))\,c_2$,
$z_1 = (a \oplus (a \gg 30))\,c_1$, $c_1 = \mathtt{0xBF58476D1CE4E5B9}$,
$c_2 = \mathtt{0x94D049BB133111EB}$, products mod $2^{64}$ [SLF; Vigna].

- (a) $\mu$ is a bijection of $\mathbb Z_{2^{64}}$ (a xor-shift and a multiplication by an odd
  constant are each invertible) and $\mathsf G$ is odd, so the state has period $2^{64}$ and the
  output takes every 64-bit value exactly once per period: perfect one-dimensional equidistribution
  over a period, at any prefix of bits (nothing is said of pairs of successive outputs, only
  $2^{64}$ of the $2^{128}$ of which occur). There is no bad seed (seed $0$ works; unlike xorshift).
- (b) *Random access.* $x_i = \mu(a_0 + (i+1)\mathsf G)$: record $i$ of a stream is a pure function
  of $(a_0, i)$.
- (c) *Streams from nearby seeds are shifted copies.* If $a' = a + k\mathsf G$ the stream of $a'$ is
  that of $a$ shifted by $k$: the windows of $L$ draws overlap iff
  $k = (a' - a)\mathsf G^{-1} \bmod 2^{64}$ (signed) lies in $(-L, L)$, by $L - \lvert k\rvert$
  draws. For seeds placed uniformly the probability is $\le(2L - 1)/2^{64}$ ($10^{-15}$ for
  $L = 10^4$). Seeds $a_0 + i\mathsf G$ ("counter" seeds, $k = i - j$) overlap for every pair with
  $\lvert i - j\rvert < L$: all pairs of $n \le L$ strata, a fraction of them beyond. Seeds that
  differ by a small *integer* $d$ do not: $k = d\,\mathsf G^{-1}$ is of order $10^{18}$
  ($\ge 7.7\times10^{15}$ for $d \le 1000$).
- (d) *Per-stratum streams.* `PHASE1.md` §2 item 4 fixes a seeded splitmix64 in the generator and
  byte-identical regeneration; how a stratum's stream is derived is not specified. What
  byte-identity with reviewable diffs needs: (R1) the stream is a function of (seed, stratum id)
  alone; (R2) adding, removing or resizing a stratum leaves every other stratum's bytes unchanged,
  which a single global stream cannot give (every later stratum shifts); (R3) a fixed draw order
  within a stratum. (c) excludes seeds spaced by a small multiple of $\mathsf G$ (a counter); a
  stratum's name has to become a number anyway, and for instance
  $\mu(\text{seed} \oplus h(\text{id}))$ with $h$ a fixed 64-bit string hash meets R1–R3 (an
  example, not a specification).

**Checked (EA.23):** Python integers, seed $0$:
$\mathtt{e220a8397b1dcdaf}, \mathtt{6e789e6aa1b965f4}, \mathtt{06c45d188009454f}$ (Vigna's reference
output); the stream of $a + 700\,\mathsf G$ equals the stream of $a$ shifted by $700$; $x_i$ by (b)
equals the sequential $x_i$; the vectorised (numpy `uint64`) and scalar versions agree; $\chi^2$ of
the top $12$ bits over $2\times10^6$ outputs $= 4144$ on $4095$ degrees of freedom (sd $90$); (c):
$k = d\,\mathsf G^{-1}$ for $d = 1, 2, 3, 1000$: $-1.02, -2.04, -3.05, -3.66\times10^{18}$, the
minimum of $\lvert k\rvert$ over $d \le 1000$ is $7.7\times10^{15}$; $100$ counter seeds with
$L = 192$: $4950$ of $4950$ pairs overlap, with $L = 64$: $4284$ ($86.5\%$). Script not committed. **Permanent:** planned, `just corpus-check` (byte identity); (c), (d) are
requirements on code that does not exist yet.
