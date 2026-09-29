# 0015: Specification gaps found while building the instrument

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** —

## Context

Three tracks have built against the specs without writing group code: the reference generator and
its committed corpus (`PHASE1.md` §4; 46 files, 31.8 of the 50 MB), `helicoid-linalg` (`PHASE2.md`
§2–§5, §7) with `cargo xtask lint` (`PHASE1.md` §3), and six derivation pages
([`maths/index.md`](../maths/index.md)) that recompute `NUMERICS.md` against mpmath. `CLAUDE.md` binds every
change to stop and ask where the documents leave a formula, a coefficient, a switch point, a tangent
order, a sign, a public item, a dependency or a feature open. Where a small reading could not wait,
a PR chose the smallest one and, in several PRs, edited the normative text to say so. Nothing below
is decided: those edits stand in the tree only until this record ratifies or reverses them (PR.1).

The gaps surfaced three ways. **The generator computes every reference from a definition** (0006 item
1), so it meets the places where a definition, a closed form and a sentence of the spec part company:
`mp.logm` near π, the precision budget at `theta:subnormal`, the tie at π. **The maths pages** check
each `NUMERICS.md` row and list what it leaves open (the "Open items" of the index, the source of most
rows). **Implementing `Dual`, `Vector`, `chol`, strided views, `mint` and the lint** forced signatures
and edge cases the specs omit. A label (`CO.n`, `SO.n`, `SE.n`, `PL.n`, `LG.n`, `EA.n`) names a
proposition, with its `Checked:` line, on those pages.

**Provenance.** The maths pages are non-normative: `NUMERICS.md` wins on any conflict. Their `Checked:`
figures, and those marked "re-run here" (each with its protocol where it appears), were run once, in
scratch scripts that are not committed; no permanent check exists until the corpus ids, the harness and
the proptests they name land. "Says" quotes the specs at `3da8149`, before the PRs that edited them
(`git show 3da8149:docs/PHASE1.md`); Found says where a passage has since changed.

**Notation.** $u = 2^{-53}$ (`f64`; $2^{-24}$ for `f32`, $u_{32}$). $\mathcal F$, $\mathcal B$ are the
forward and backward error of `NUMERICS.md` §11 in units of $u$ (EA.3; $\mathcal B_M$ for a matrix
input, $d_M$ its Frobenius distance to SO(3), EA.9). $\eta = \|q\|^2 - 1$. $K = \begin{bmatrix}0 &
-1\\ 1 & 0\end{bmatrix}$. For a coefficient, $\theta_s$ is its switch point, $E_x$, $E_s$ the errors of
its exact and series arms, $E_\times$ their common value at the crossing (CO.10). $\kappa$ is a
condition number.

No closed form of `NUMERICS.md` §3–§7, series coefficient or Jacobian was found wrong; one stated
identity is (§8's $\alpha/\|n\times m\|$ is $r/2$, NU.6(d)). False or unsatisfiable as written: §4's
continuity bound (NU.4) and "naive cancellation" constants (NU.7); `PHASE1.md` §2 item 3's 40-digit
recheck below $\theta \approx 10^{-19}$ (P1.8) and §10's "negated half of every `so3_log` stratum"
(P1.6(b)). Imprecise, not false: §11's matrix-input sentence, right about the axis and mis-motivated
(NU.1(e)); §3.2's "exactly scale-invariant", true of the exact arm only (NU.6(a)); §10's `b` rule,
which passes by a margin of 0.11 (P1.6(a)). Every other gap is a missing rule, an ambiguity, or a
requirement the instrument cannot measure as worded.

## Decision

None yet. Each gap carries a **Recommend.** line: the resolution its author would accept, and why. The
maintainer accepts, amends or rejects each one in *Open questions*; silence is not acceptance. A `ready`
version of this record then holds only the accepted resolutions, as the edits they imply (a formula
change is a `NUMERICS.md` edit and a record), and ratifies or reverts the provisional edits of PR.1 in
the same step. Gaps are numbered by document: **P1** `PHASE1.md` and the instrument, **NU**
`NUMERICS.md`, **PH** `PHASE2.md`/`PHASE3.md`, **PR** process. They come in two groups, each in the
order of the body. *Decisions* are ordered by the first upcoming work each blocks: the harness
(`PHASE1.md` §5), the sweep and `generated.rs` (§6), the seeded defects (§10), the rest of Phase 2 (§6),
Phase 3 (`PHASE3.md` §3–§9), Phase 5; "—" blocks nothing but is a real choice or amends a hard rule.
*Ratify* is the tree's own edits and readings, which block nothing and which the recommendation would
keep; each lists its options too.

**Decisions**

| Gap | About | Recommend | First blocked |
|---|---|---|---|
| NU.1 | §11: floors, backward error, p99, matrix input | relative error, floor = smallest normal for tangents and scalars | harness, sweep |
| P1.1 | `f32` has no corpus inputs | `f32`-exact companions for `coeff_*`, rounding elsewhere | harness, sweep |
| P1.2 | record schema, manifest identity | ratify; fix §4.1, §4.2 example, `shape` rule | harness |
| NU.2 | §3.3 against §1; `so3_act` at `q:nonunit` | unit `q` in domain; keep `q:nonunit`, record its 2\|η\| cost | harness, Phase 3 |
| NU.3 | §3.2 `copysign` against `w < 0` at `w = −0` | sign-bit flip; add `w = −0` records | harness, self-test, Phase 3 |
| P1.3 | §6 grid ends below the optimum | end at π − 0.1 (per coefficient); `m ≤ 8`; report `m ≤ 16` and a cap at 1 | sweep |
| P1.4 | derivative variable | the branch variable, as generated (`r`: n², w fixed) | sweep |
| P1.5 | no id for α, β, cos θ/2 | add three ids and series rows | sweep |
| NU.4 | §4 continuity bound | sum of the two arms' errors | sweep, Phase 3 |
| NU.5 | one mask per call-site group | a switch per coefficient | sweep, Phase 3 |
| NU.6 | the `r` kernel | mask on `s = n²/w²`, `w > 0`, safe argument | sweep, Phase 3, 5 |
| PH.1 | `Dual` conventions, `real_*` ids | ratify; add `real_*` ids | sweep |
| NU.7 | §4 naive-cancellation column | keep exponents, add `p + 2` column | sweep (predictions) |
| P1.6 | two seeded-defect rules: one thin, one with a hole | window from `theta:1e-7`; "but `q:w0`" | self-test |
| PH.2 | linalg readings: `vee`, `from_cols`, `Mat2`, `svd3`… | ratify; add the `svd3` ordering; drop `Mat2` | Phase 2, 3 |
| NU.8 | `from_matrix` normalization and sign | divide by the norm; `w ≥ 0` | Phase 3 |
| NU.9 | §12 rows: `Exp`, SE(2), `Log` Jacobians | add the rows | Phase 3 (R6) |
| NU.10 | SE(2) Jacobians, action rows, twin | adopt `so2-se2.md` §6, as derived | Phase 3 |
| NU.11 | SO(2), SE(2) `Log` at ±π | keep `atan2`'s sign-bit tie; sample both zeros | Phase 3, corpus |
| PH.3 | `se2_coeffs(θ²) -> (α, β)` | `(α, a)`; β formed by the caller | Phase 3 |
| P1.9 | §4.4 readings, sample counts, size | ratify, write the counts in | — |
| PH.5 | strided reads panic in release: D11, R6, §12 widened | accept, explicitly | — |

**Ratify**

| Gap | About | Recommend | First blocked |
|---|---|---|---|
| PH.4 | Cholesky (`NUMERICS.md` §15) | ratify | — |
| P1.7 | `mp.logm` as a reference | ratify the Newton reference | — |
| P1.8 | precision budget below 1e-12 | ratify guard digits | — |
| P1.10 | lint readings, 0003's item 5 and doctest name | ratify; correct 0003 while `ready` | — |
| PR.1 | edits to normative text by PRs | ratify or revert row by row | — |
| PR.2 | maths pages describe replaced text | resync after PR.1 | — |
| PR.3 | unmeasured: aarch64, macOS, CI time | keep visible | `PHASE1.md` §11 |

## Decisions, in blocking order

Each gap: the passage, what was found, the options, a recommendation and what it blocks.

**NU.1 The metric of §11 is underspecified.**

- **Says.** "$\|y\|_{\text{floor}}$ is the function's scale (1 for rotations and Jacobians,
  $\|\rho\|$-scale for translations) so an exact zero is not divided by." "**Backward error** for
  `Log`: $\|\mathrm{Exp}(\widehat{\mathrm{Log}}\,X) \ominus_R X\|$ in $u$." "Matrix inputs
  (`from_matrix`) are ill-conditioned near $\pi$ in the axis; their strata report backward error only."
  "**Per stratum:** max and p99; never a mean."
- **Found.** (a) *No floor for a tangent or a scalar.* `Log`'s φ has norm θ; with floor 1, `theta:1e-k`
  measure absolute error (a total loss at θ = 1e−12 reads 9e3u, not 9e15u) and, below θ ≈ 1e7u =
  1.1e−9, cannot show the seeded `acos` defect of §10 ("≥ $10^7\,u$"); on `theta:subnormal` the floor
  must reach the smallest normal number or F measures quantization (EA.4(b)). The `coeff_*` outputs are
  none of rotation, Jacobian or translation, and the sweep compares value and derivative errors: at
  branch value 0 the values run from 1/2 (k, a) to 1/120 (e) and the derivatives from 1/24 (a) to
  1/2520 (e) in magnitude, so with floor 1 a relative error of 2e3u on e's derivative reads 0.8u (this
  record's observation; EA.4(c) assumes them relative). (b) *B is not computable from the corpus:* it
  needs `Exp` at the subject's output; EA.5 gives ‖J_r(φ)e‖ from the stored φ, valid for ‖e‖ ≪ 1 and
  after reducing the output to the logarithm nearer φ (an output on the other branch reads 2π). (c) *F
  from a 30-digit string:* parsing the reference to `f64` first quantizes F to an ulp and blinds the
  exact no-regress bar (EA.19(b)); the difference needs ≥ 32 digits. (d) *p99* has no method; at n = 64
  the nearest-rank p99 is the max (EA.11(d)). (e) *Matrix input:* M ↦ nearest-rotation quaternion is
  1/(2√2)-Lipschitz at every θ, π included; the axis orientation is what is ill-conditioned, and
  backward-only scoring is right because the reference is undetermined below O(d_M) (EA.9), not because
  the map is ill-conditioned. §11 defines no backward error for a matrix (EA.3's B_M is one reading).
  (f) *Translations and exact zeros:* the ‖ρ‖-scale floor does not survive stacking ρ under φ in a
  tangent: one norm over [φ; ρ] hides an error in ρ where ‖ρ‖ ≪ ‖φ‖ (the `rho:1e-6` cells), and no rule
  says what a nonzero subject scores against an exactly zero reference.
- **Options.** (a) Relative error with floor = smallest normal (2⁻¹⁰²² `f64`, 2⁻¹²⁶ `f32`) for
  tangents and scalars, 1 kept for rotations and Jacobians, or floor 1 throughout. (b) EA.5's
  first-order form with the reduction, exact form through a double-double `Exp` for large residuals, or
  a stored second reference. (c) Double-double differences in the harness. (d) Nearest-rank p99
  reported beside n, never a bar, or interpolation. (e) Reword, or keep §11's sentence. (f) φ and ρ
  scored separately, or one norm over the stacked tangent.
- **Recommend.** (a) smallest normal; (b) first-order with reduction, exact fallback; (c); (d) nearest
  rank; (e) reword per EA.9, with B_M defined; (f) φ and ρ scored separately, ρ against its ‖ρ‖-scale
  floor, and an exactly zero reference scores 0 for an exactly zero subject and fails (F = +∞)
  otherwise.
- **Blocks.** The harness and its result CSV (§5); the `acos` defect (§10); the sweep objective (§6).

**P1.1 `f32` has no inputs in the corpus.**

- **Says.** §2 item 1: "Every corpus input is a binary64 value serialized with Python's `float.hex()`;
  the reference is computed at that exact value". §6: "Every coefficient in `NUMERICS.md` §4, plus
  SE(2)'s $\alpha, \beta$; each for `f64` and `f32`." §5: "`eval(fn_id, &Record, Precision) -> Output`".
- **Found.** An `f32` subject must round each input, so it evaluates at another point than the
  reference and its error includes κ·u₃₂ of input rounding (EA.19(c); up to θ/2 for `Exp`). Some strata
  do not survive: `theta:pi-1e-k` for k ≥ 8 round to one value, fl₃₂(π) = π + 8.7e−8, an angle above π
  (binary32 round-to-nearest of the binary64 π; the `f32` spacing at π is 2⁻²² = 2.4e−7);
  `theta:subnormal` becomes 0. No section says which inputs an `f32` subject gets.
- **Options.** (a) an `f32`-exact companion of each id, rounded by the generator with the reference at
  the rounded value: about 2 MB for the seven `coeff_*` ids at the committed record size, but doubling
  the group ids would take the 31.8 MB corpus past 50 MB. (b) Keep the binary64 corpus, let the
  subject round and report F with its input term. (c) (a) for the `coeff_*` ids the sweep needs, (b)
  elsewhere until Phase 3 has `f32` code to measure.
- **Recommend.** (c); file naming and schema for the companion are decided with it.
- **Blocks.** `f32` rows of the harness (§5), the `f32` half of the sweep (§6), `PHASE3.md` §9.

**P1.2 Record schema and manifest identity.**

- **Says.** §4.2, initial: "`MANIFEST.json` records generator git revision". "Matrices are
  column-major flat arrays with a sibling `"shape"`." The example lists keys `id, stratum, in, out`;
  §2 item 4 requires "sorted keys". §4.1: "`uv.lock` pinning the interpreter and `mpmath`".
- **Found.** A revision cannot appear in the commit it names, so it cannot be byte-stable. §4.2 now
  says "a SHA-256 over its sources": the sorted (path, length, bytes) of `pyproject.toml`, `uv.lock`,
  `.python-version` and `gen/**/*.py`, so any edit to them, a comment included, regenerates the corpus in
  the same PR. `uv.lock` pins mpmath but only `requires-python ==3.12.*`; byte identity was checked
  across CPython 3.12.3 and 3.12.14 and the manifest records "3.12". §4.2 does not say where `shape`
  lives or how two matrices coexist: it is a sibling key in the same `in`/`out` object, one matrix per
  object. Sorted keys give `id, in, out, stratum`, not the example's order.
- **Options.** (a) Ratify §4.2 as edited; correct §4.1, the example's key order, and state the `shape`
  rule. (b) Hash only the files that determine bytes (not separable from comments). (c) Drop the
  identity and rely on `corpus-check`.
- **Recommend.** (a).
- **Blocks.** The Rust reader of the corpus (§5) and the `kind` filter of the envelope (§8).

**NU.2 §3.3 and §1 disagree on a non-unit q, and the corpus applies `q:nonunit` to `act`.**

- **Says.** §3.3: "$v' = v + 2w\,(u \times v) + 2\,u \times (u \times v)$. `act_many` forms $R(q)$
  once and applies it per point; its reference twin is the per-point `act`." §1: "$R(q) = (w^2 -
  \|u\|^2)\,I + 2\,u u^\top + 2w\,[u]_\times$". `PHASE1.md` §4.3: `so3_act` is "$R(q/\|q\|)\,p$".
- **Found.** They differ by (1 − ‖q‖²)v = −ηv (SO.2(d)); the vector form errs against R(q̂)v by
  η(R(q̂) − I)v and the matrix form by ηR(q̂)v, up to 2⁻³⁹‖v‖ and 2⁻⁴⁰‖v‖ at the `from_wxyz_unchecked`
  bound (|η| ≤ 2⁻⁴⁰, §3.6; SO.15(c)), so the `act_many` twin (§14) needs normalized inputs. The corpus
  gives `q:nonunit` (|η| = 2⁻⁴⁵, inside that bound) to `so3_act`, which §4.4 at `3da8149` does not name
  it for (its purposes are `Log`'s scale invariance and `from_*` normalization); the SE_N(3) PR added
  `sen3_ad` to §4.4's list (in the tree; PR.1). Re-run here (the 64 committed `so3_act` `q:nonunit`
  records, stored q unnormalized, exact arithmetic at 60 digits, norm-wise error against the stored
  reference in units of 2⁻⁵³, max over the records; η = ±2⁻⁴⁵ to 0.4%): §3.3's form differs from the
  reference by up to 511u (2|η|) and §1's by 257u (|η|). That is the cost of not normalizing, not a
  defect of a subject: any implementation of either form as written misses the reference by that much.
- **Options.** (a) `act`, `act_many`, `to_matrix` take q with |η| ≤ 2⁻⁴⁰ (§3.6's bound,
  `debug_assert!`) and do not normalize; §3.3 says which form `act_many` forms; the twin runs on
  normalized q or carries the |η| term in its tolerance. `q:nonunit` stays on `so3_act`, and its 511u
  and 257u are recorded as that cost; the no-regress bar pins them. (b) `act` applies (3 − ‖q‖²)/2 first
  (one multiply, SO.15(c)): the error at the stratum returns to rounding level. (c) Define `act` by
  R(q) and change the `so3_act` reference (and NU.1's metric) to R(q)p, naming the form. (d) Drop
  `q:nonunit` from `so3_act` and `sen3_ad`: it removes committed strata because a spec form misses
  their reference, which D7 and 0006 (*Consequences*) forbid ("a stratum is never deleted or narrowed
  because it fails"); it needs an explicit exception or amendment from the maintainer.
- **Recommend.** (a): the stratum is in `act`'s domain, so it stays and the 2|η| cost stays visible;
  (b) if a consumer needs `act` exact on near-unit input. Not (d). `sen3_ad` follows the same reading;
  its inclusion in §4.4 is a PR edit for PR.1 to ratify or revert.
- **Blocks.** The `so3_act` rows of the harness; `SO3::act`, `act_many`, `to_matrix` and the twin
  (`PHASE3.md` §4; `NUMERICS.md` §14).

**NU.3 §3.2: `copysign` and `w < 0` differ at `w = −0`.**

- **Says.** "If $w < 0$, $q \leftarrow -q$ (implemented as `copysign`, not a branch). At $w = +0$
  nothing flips: **at exactly $\theta = \pi$, `Log` is a function of the quaternion, not of the
  rotation** — $q$ and $-q$ return $\pm\pi\hat n$."
- **Found** (SO.5(d)). `copysign(1, w)` is −1 at w = −0 while `w < 0` is false (IEEE 754). Under the
  sign-bit reading q = (+0, u) returns +πû and −q = (−0, −u) flips to (+0, u) and also returns +πû;
  only the `w < 0` reading returns ±πû as the sentence says. Both are logarithms of R(q). The generator
  implements "w < 0 negates, w = +0 stays" (its README); its `q:w0` has w = +0 only (64 records, none
  with −0), so the corpus cannot tell the readings apart. `so3_log_at_w0_is_a_function_of_the_sign`
  (`PHASE3.md` §9) tests the answer.
- **Options.** (a) The sign-bit flip as written: Log(q) = Log(−q) for every q, ±0 included; reword the
  sentence, add the `w = −0` (negated) records to `q:w0`. (b) A mask, `select(w.lt(S::zero()), −q, q)`: the
  sentence as written, but Log(q) = −Log(−q) at w = ±0 and the result then depends on the axis sign;
  drop "copysign".
- **Recommend.** (a): it is what "implemented as `copysign`" says, gives one property test for all q,
  and needs no comparison.
- **Blocks.** `SO3::log` (`PHASE3.md` §4) and its two tests; the corpus `q:w0`; §11's sign-invariant
  metric at w = +0; the flip defect of §10.

**P1.3 The switch grid ends below the optimum.**

- **Says.** §6: "Series terms $m \in \{1, \dots, 8\}$; switch point in the branch variable on a log
  grid of 64 points per decade, spanning $\theta \in [10^{-8}, 1]$." Objective: "the max over
  `theta:dense` and all `theta:*` strata of max(value error, derivative error)".
- **Found.** CO.10 (`f64`, m = 8): the optimum θ_s of k, b, d, e (predicted 1.07–1.46, measured
  1.04–1.41) lies above the grid's end, and so does that of all six in `f32` (predicted 2.4–5.1;
  measured, the series arm serves all of [1e−3, π − 0.1] for b, d, e). A cap at 1 costs e a factor 5.9
  on the value and 7.7 on the derivative in `f64` (2.2e3u to 1.7e4u) and a factor 8e3 on the derivative
  in `f32` (2u to 1.6e4u); with m ≤ 8 that `f64` derivative cannot go below ≈ 2e3u. `coeff_series`
  already holds 16 terms. The strata reach past 1: `theta:1e0` ends at π − 0.1 and `theta:pi-1e-k`
  (k = 1…12) run from there to π − 1e−12, and the objective takes "all `theta:*` strata". The switch
  is strict (§6's `below`, θ² < …), so θ = θ_s itself takes the exact arm. The grid is stated in θ, the
  switch in θ². The objective minimizes each coefficient's own error, but a consumer multiplies it by a
  word of size θᵖ: with the exact arms of b, d, e alone Q (the coupling block of SE(3)'s Jacobian) errs
  ≲ (16.3/θ + 20.8)u against 6u/θ² for b (SE.15); the derivative through Q was not examined.
- **Options.** Grid end: 1; π − 0.1, the end of the continuous strata; or each series' own radius
  (r: s < 1, θ < π/2 for a unit q, CO.16(b); c: θ < 2π; the others are entire). Terms: m ≤ 8, or ≤ 16,
  which lowers the floor of e's derivative (CO.10: a larger m raises θ_× and lowers E_× down to the
  floor) for twice the constants per coefficient and precision in `generated.rs` and a longer Horner
  chain, both unmeasured. Objective: the coefficient's (0004 item 1), or add the value strata of
  `sen3_exp`, `sen3_jl`.
- **Recommend.** The grid ends per coefficient in its own branch variable: π − 0.1 for all but r,
  whose grid ends at s = 1 (beyond it the series diverges and the extra candidates are noise). At
  π − 0.1 the `f32` optima of b, d, e sit on the boundary: the series serves every sampled θ ≤ π − 0.1
  and the `theta:pi-1e-k` strata take the exact arm; nothing above π is searched (P1.1: `f32` collapses
  those strata for k ≥ 8). Keep m ≤ 8 (§6's; 0004 item 1 fixes no bound) and the coefficient-level
  objective. Write into the sweep CSV, beside the chosen row, what a cap at 1 and m ≤ 16 would have
  cost, so the trade (a floor of ≈ 2e3u on e's derivative against constants and Horner length) is
  measured, not argued.
- **Blocks.** The sweep and `generated.rs` (0004 items 1, 3; `PHASE3.md` §3).

**P1.4 The derivative variable is unnamed.**

- **Says.** §4.3, initial: "the definition at 120 digits; derivative by `mp.diff`" (input θ). §6: "the
  derivative is `Dual<S, 1>` through the candidate and the reference is `mp.diff`". `NUMERICS.md` §4:
  "Each is a function of a branch variable ($\theta^2$, or $n^2$ for $r$)".
- **Found.** d/dθ = 2θ d/dz (CO.13), and the digits `mp.diff` keeps differ: at θ = 1e−12, 67 in θ
  against 79 in z (CO.17). The corpus takes z; §4.3 now says the outputs are "`value` and `d_branch`:
  the derivative with respect to the branch variable of `NUMERICS.md` §4, by `mp.diff` at the exact
  real branch value of the exact binary64 input (θ², not `fl(θ·θ)`)", and for `coeff_r` d/dn² at fixed
  w. For r that is not the variable NU.6 recommends for the switch (s = n²/w²): with r = (2/w)S(s),
  d/dn² at fixed w is w⁻² times d/ds at fixed w.
- **Options.** (a) The variable the kernel takes, as generated: kernels take θ², `Dual<S, 1>` is
  seeded on it; for r, on n² with w a constant. (b) θ: a √ before the kernel, `Dual` seeded on θ.
  (c) For r only, d/ds or d/dw.
- **Recommend.** (a), stated in §6 too, r included: the seed is n² with w constant, which is what
  `coeff_r`'s `d_branch` stores, so no corpus byte changes. NU.6's s is the variable of the switch mask,
  which reads value parts only and does not change what is differentiated. (c) would regenerate
  `coeff_r` (a w² factor).
- **Blocks.** The sweep: what `Dual` is seeded on, and so what "derivative error" measures.

**P1.5 Three coefficient functions have no corpus id.**

- **Says.** §6: "Every coefficient in `NUMERICS.md` §4, plus SE(2)'s $\alpha, \beta$". `NUMERICS.md`
  §3.1: the series arm "computes $\cos\tfrac{\theta}{2}$ from its own series (generated alongside
  $k$)"; §6: α, β "(both 0/0 only: series from the generator, switch generated)". §4.3 lists no id
  for them; `coeff_series` has seven rows. §0.0 lists them as "Missing".
- **Found.** α = sin θ/θ, β = θ·a (odd in θ, PL.1, PL.2) and cos θ/2 (the series σ₀(z/4), CO.2) have no
  reference and no series. cos θ/2, not k, limits the derivative switch of `exp_coeffs` (1.1 against
  1.5 at m = 8, CO.10).
- **Options.** (a) Ids `coeff_alpha`, `coeff_beta`, `coeff_cos_half` and three series rows (a corpus
  addition, so a record). (b) Sweep only what cancels; give α, β a fixed series switch. (c) Derive
  α from b: α = 1 − zb cancels near |θ| = π, as γ = 1 − zc does (PL.11(a)).
- **Recommend.** (a); PH.3 decides whether β is needed at all.
- **Blocks.** The sweep's scope (§6); `exp_coeffs` and `se2_coeffs` (`PHASE3.md` §3).

**NU.4 §4 continuity is stated as one arm's error.**

- **Says.** "At every generated switch point, $|\text{series} - \text{exact}|$ is at most the recorded
  max error of that coefficient; `branch_continuity_*` tests assert it."
- **Found.** Each arm is within its own error of f, so they differ by up to E_x + E_s, up to 2E_× at
  the crossing. Measured at the predicted crossing (m = 8; b, c, d, e; 6000 samples) the largest jump is
  9.6, 127, 40, 135u against max(E_x, E_s) = 6.6, 73.5, 26.3, 97.7u: 1.4–1.7× (CO.12). A test that
  asserts the sentence fails on a correct switch.
- **Options.** (a) The bound is the sum of the two arms' recorded errors at the switch (the sweep CSV
  carries both). (b) 2× the recorded maximum. (c) Read "recorded max error" as the maximum jump.
- **Recommend.** (a).
- **Blocks.** The columns of `thresholds.csv` (§6); `branch_continuity_*` (`PHASE3.md` §9).

**NU.5 One mask per call-site group, or one per coefficient.**

- **Says.** §4: "**Call sites evaluate coefficients in groups, inside one `branch`:** `exp_coeffs` →
  $(k, \cos\tfrac\theta2)$; `jr_coeffs` → $(a, b)$; …". `PHASE3.md` §3: "Each is **one** `S::branch`
  over a tuple", and "`generated.rs` holds a `Switch` per coefficient per precision".
- **Found.** One mask means one θ_s shared by the members, but their optima differ (m = 8, derivative:
  a 0.9, b 1.1; b, d, e 1.1, 1.3, 1.4; k 1.5, cos θ/2 1.1). Worst member against its own optimum, value
  and derivative (`f64`, m = 8): (a, b) 1.1×, 1.7×; (b, d, e) 2.6×, 3.0×; (k, cos) 1.0×, 1.5×; (b, d)
  1.6×, 1.7×; at m = 4 (D12's) up to ≈ 10× (CO.18). Neither text says which reading is meant.
- **Options.** (a) One mask, one `Switch` per group, minimized on the group's max. (b) A `Switch` per
  coefficient, each member selected by its own mask inside one call-site function. (c) The text as is.
- **Recommend.** (b): it is what 0004 measures and `generated.rs` holds. Its cost is more masks per
  call site (`q_coeffs` has three members): a `bool` mask runs one arm per mask, a lane mask both arms
  of each, and CO.18 finds a shared mask affordable at m = 8 (worst member ≈ 3×) but not at D12's four
  terms (≈ 10×). Shared exact-arm subexpressions (θ, sin, cos) are evaluated once and give the same
  bits (D16). The bench gate (§9) prices the extra masks in Phase 3; where the sweep CSV shows a
  group's optima within a grid step of each other, that group may share one mask. Reword "one
  `S::branch` over a tuple".
- **Blocks.** The shape of `generated.rs` (§6); the kernel (`PHASE3.md` §3).

**NU.6 The `r` kernel: variable, sign of w, safe argument, and §8.**

- **Says.** §4: "Each is a function of a branch variable ($\theta^2$, or $n^2$ for $r$)"; r's series is
  in "$s = n^2/w^2$". 0003 item 3: "the exact arm evaluates at `select(small, 1, θ²)`, so no arm is ever
  non-finite". §8: "the ratio $\alpha/\|n \times m\|$ is §4's $r$ kernel (with $w = n\cdot m$,
  $n^2 = \|n\times m\|^2$)".
- **Found** (CO.16). (a) The series holds for w > 0 and s < 1 only, so a mask on n² is not
  scale-invariant: for 1e−100(cos 1, sin 1 n̂), n² = 7e−201 passes any n² threshold while s = 2.43 and
  the 8-term series returns −21.7 for 0.642. §12 limits `Log` to unit q, so this bites outside the
  domain; but §3.2's "exactly scale-invariant" is a property of the exact arm alone (a near-unit q, the
  case §3.2 draws from it, is unaffected). (b) At w = 0, s = ∞: a lane mask evaluates the series arm,
  which needs the safe argument that 0003 and §4 state for the exact arm only. (c) `S2Chart::local` has
  w = n·m of either sign: at n² = 1e−6, w = −1 the exact value is 2·atan2(1e−3, −1)/1e−3 = 6281.19 and
  the 8-term series arm gives −2.0000 (evaluated here from the two formulas); `Log` never sees w < 0
  (the flip), so the mask must include w > 0. (d) With ν = ‖n×m‖ (the kernel's n) and w = n·m,
  α/‖n×m‖ = atan2(ν, w)/ν = r/2, not r. (e) The corpus stores d/dn² at fixed w (P1.4), which is not the
  variable of the series, s.
- **Options.** (a) Branch variable s = n²/w² (mask `n² < s_sw·w²`, division-free), `w > 0` in the mask,
  safe argument on both arms, §8 corrected to r/2. (b) Keep n² and state the domain (unit q, w > 0).
  (c) Split `log_ratio` (unit, w ≥ 0) from a chart kernel.
- **Recommend.** (a), with s the variable of the mask only: `coeff_r`'s `d_branch` (d/dn² at fixed w,
  P1.4) stays as committed. The corpus gains a `w < 0` stratum for `coeff_r` (CO.16(c)'s n² = 1e−6,
  w = −1 among them), so the `w > 0` clause of the mask has one (D7).
- **Blocks.** The sweep for r (the unit of its θ_s); `log_ratio` (`PHASE3.md` §3); `S2Chart::local`
  (`PHASE5.md` §2).

**PH.1 `Dual` conventions: edited by the PR, to be ratified.**

- **Says.** `PHASE2.md` §3, initial: "`copysign(x, s)` → $\mathrm{sgn}(x_v)\,\mathrm{sgn}(s_v)\,x_d$";
  "The value path of `Dual<S, N>` is bitwise identical to the plain `S` evaluation"; "`sqrt` at 0 has an
  infinite derivative by the rule above"; "tested to second order on `sin_cos`, `atan2`". No rule for
  × and ÷.
- **Found.** The `Dual` PR edited §3 and `NUMERICS.md` §12: sgn(s_v) = −1 iff the "**sign bit** of
  $s_v$ is set" (the derivative of the value returned at −0.0); the value path is bitwise identical
  "except the sign and payload of a NaN produced by arithmetic or a `libm` call" (unspecified in Rust;
  release commutes operands); `sqrt` at 0 is ±∞ where d ≠ 0 and NaN where d = 0; "`a / b` →
  $(a_d - q\,b_d)/b_v$ with $q = a_v/b_v$", the other form differing in the last bit on 57% of random
  inputs (EA.15), which moves every derivative the sweep measures (the product rule
  stays unstated); `atan2`'s derivative is accurate only while x_v² + y_v² is normal (about 1e±154
  `f64`, 1e±19 `f32`), `a / b` is NaN once q overflows, none asserted. Second order through an exact arm
  loses θ⁻⁶ (b: ≈ 1e4u·θ⁻⁶, no digit left near θ = 1e−2) while the sweep measures first order only
  (EA.17(d)). `PHASE2.md` §8 names `dual_matches_mpmath_derivative` on "corpus ids `real_*`", which
  `PHASE1.md` §4.3 does not list and the generator lacks: an inline 27-row fixture stands in.
- **Options.** (a) Ratify all; add `real_*` ids (`sqrt`, `sin_cos`, `atan2`, quotient) to §4.3.
  (b) Also scale `atan2`'s derivative against overflow: it changes the normative formula and the
  accuracy. (c) Revert the quotient rule to the textbook form.
- **Recommend.** (a); (b) when a consumer feeds `Dual<f32>` `atan2` outside 1e±19. While 0003 is
  `ready`, correct its item 5 and the README row for the NaN exception (P1.10).
- **Blocks.** The sweep's derivative column (the quotient rule fixes it); `real_*` conformance.

**NU.7 The "naive cancellation" column is not a bound, and is about values.**

- **Says.** The column: "$\sim 2u/\theta^2$ naive", "$\sim 6u/\theta^2$", "$\sim 12u/\theta^2$",
  "$\sim 24u/\theta^4$ naive, $\sim 24u/\theta^2$ rewritten", "$\sim 360u/\theta^4$", and "evaluating the
  definition in the exact arm keeps $-\log_{10}$ of that relative error in digits".
- **Found** (CO.6, CO.14; `f64`, 6e4 samples). The exponents are right; the constants are not bounds:
  sampled maxima a 1.0 (column 2), b 6.0 (6), c 48 (12), d naive 36 (24), d rewritten 45 (24), e 367
  (360). The column is not built by one rule (a, b, c, d count the dominant rounded operand once, e counts
  all its operands, which would give 48 for d). Through `Dual` every coefficient, k, a, r included,
  loses θ^−(p+2): measured p + 2 = 2, 2, 4, 4, 4, 6 for k, a, b, c, d, e, which the spec does not state
  and which sets the switch of k, a, r.
- **Options.** (a) Keep the exponents, replace the constants by a pointer to the measured sweep CSV, add the
  derivative column p + 2. (b) Leave the "∼". (c) Restate the constants as measured maxima.
- **Recommend.** (a).
- **Blocks.** Nothing hard; the sweep's predictions (CO.10) and readers' digit counts.

**P1.6 Two seeded-defect rules: one thin, one with a hole.**

- **Says.** §10: "`b` by its definition, no series | error curve in `theta:1e-8`…`1e-2` fits
  $\theta^{-p}$ with $p \in [1.8, 2.2]$"; "`Log` without the $w < 0$ flip | the negated-quaternion half
  of every `so3_log` stratum fails".
- **Found.** (a) *Thin.* b evaluates to 0, relative error exactly 1, below √(6u) = 2.6e−8, inside
  `theta:1e-8`: simulated fits give p = 1.91–1.94 (2.00 without that stratum). That is inside
  [1.8, 2.2], 0.11 above its lower edge, so the rule as written passes (CO.6, `Checked:`: biased low by
  0.07). (b) *A hole.* The negation of a `q:w0` sample has w = −0 and atan2(n, −0) = atan2(n, +0), so
  under the `w < 0` reading (the corpus's, NU.3) a `Log` without the flip returns what the flip would:
  nothing detects it. `q:w0` in the corpus has no negated half (64 records); every other `so3_log`
  stratum has one. Under the sign-bit reading the negation would detect the defect.
- **Options.** (a) Leave the window; or fit `theta:1e-7`…`1e-2`, or skip saturated records. (b) Restrict
  the rule to the strata other than `q:w0`, or add `w = −0` records (NU.3).
- **Recommend.** For (a), the window from `theta:1e-7`: a robustness choice, not a repair. D7 forbids
  deleting or narrowing a stratum; `theta:1e-8` stays in the corpus, the harness and the results, and
  only the self-test's fit window moves. If the maintainer reads D7 more widely, keep the window as
  written and accept the margin. For (b), "but `q:w0`" now (the stratum stays; only the rule's scope
  says where it can detect), revisit after NU.3.
- **Blocks.** `just conformance --self-test` (§10; appendix steps 3, 5).

**PH.2 Readings in the linear-algebra surface.**

| Passage | Reading in the tree | Options | Recommend |
|---|---|---|---|
| §4: "`vee(m) -> Vec3` (skew part, no symmetrization check)" | reads entries (2,1), (0,2), (1,0); on a non-skew m that is not ½(m − mᵀ); `vee(hat v)` is exact | three entries; or ½(m − mᵀ), six reads | three entries; reword |
| §4: `from_cols`, `from_rows`, `col`, `row`, no argument types | `Vector` in and out | `Vector`; or `[S; N]` | `Vector` |
| §0: "`Mat2`/`Mat3` adjugates"; §4, `API.md`: `Mat3::inverse_adj` only | `Mat3` only | drop "Mat2" from §0; or add `Mat2::inverse_adj` | drop; add with a consumer |
| §6: svd3 "signed so $\det U = \det V = +1$ and only $\sigma_3$ may be negative (`nearest_rotation` is $U V^\top$ directly)" | silent on \|σ₃\| ≤ min(σ₁, σ₂); for Σ = diag(1, 1, −5) UVᵀ is at squared distance 36 against 20 for the nearest rotation (SO.13(c)) | state the ordering and add a stratum; or leave the order free and have `nearest_rotation` flip the axis of the smallest \|σ\| | the ordering and a stratum |
| §7: `Point` at "`N` = 2, 3, 4" | `mint` has no `Point4`: N = 2, 3; §7, its status row and `API.md` §2 edited | N = 2, 3; or convert `Point<S, 4>` through `Vector4` | N = 2, 3 |
| §5: "Bounds are checked by safe slice indexing" | not enough: a write past a block's edge lands in its neighbour inside the same slice, so row and column extents are checked too (panics: PH.5); overlapping strides are memory-safe, documented, not rejected; strides are non-negative | as read; or reject overlapping strides; or signed strides | as read; signed strides not now |

- **Blocks.** `svd3` (Phase 2 §6) and every Phase 3 use of these types.

**NU.8 `from_matrix`: "then normalize" is only good near a rotation, and has no sign rule.**

- **Says.** §3.4: "select the largest of $\{\mathrm{tr}R, R_{00}, R_{11}, R_{22}\}$ … and extract from
  that pivot …, then normalize." §3.6: "Construction from external data and the explicit `renormalize`
  apply the first-order Newton step $q \leftarrow q\,(3 - \|q\|^2)/2$".
- **Found.** With η = ‖q‖² − 1 the step leaves ‖q′‖² − 1 = −¾η² + ¼η³ (exactly, SO.14), below rounding
  only for |η| ≤ 2^−26.3 (`f64`; 2^−11.8 `f32`; measured 4u at 2^−26, 6e3u at 2^−20). Shepperd's
  extraction from a matrix within ε of SO(3) has |η| ≤ 6ε (SO.12(d)), so `from_matrix` reaches rounding
  level only for ε ≲ 2e−9, and at η = 2 the step returns 0 (SO.14). `PHASE3.md` §4's fixtures, locus-tag's
  degenerate near-singular IPPE matrices, need not be near SO(3). The extraction fixes only the pivot
  component positive: the output jumps between q and −q at pivot ties, for an exact rotation too (SO.12,
  *Sign*). The corpus reference is w > 0, else the first nonzero of x, y, z positive; the harness scores
  backward error only.
- **Options.** (a) Divide by the norm in `from_matrix` (one √, one division; keep Newton for
  `renormalize` and `from_wxyz_*`). (b) Newton after an `svd3` projection, selected by a mask. (c) State
  the domain (ε ≲ 2e−9, `debug_assert!`). Sign: (i) none, (ii) canonical w ≥ 0 (moves the jump to w = 0,
  where `Log` already jumps), (iii) canonical as the corpus.
- **Recommend.** (a) and (ii).
- **Blocks.** `SO3::from_matrix` and `from_matrix_never_iterates` (`PHASE3.md` §4, §9).

**NU.9 §12 has no row for `Exp`, SE(2), or the `Log` and `⊖` Jacobians.**

- **Says.** §12 lists "`jr_inv`, `jl_inv` (SO(3), SE_N(3)) | $\theta < 2\pi$ | unspecified finite value".
  `PHASE2.md` §2: "no routine produces a non-finite value from finite in-domain input". `API.md` R6:
  "Every function with a restricted domain carries a `# Domain` rustdoc section naming it
  (`NUMERICS.md` §12)".
- **Found.** (a) `Exp`: θ² = φ·φ overflows for ‖φ‖ > 1.3e154 (`f64`), 1.8e19 (`f32`) and the exact arm
  returns NaN (CO.15(c)). The series arm overflows far below that, for ‖φ‖ ≳ 1e22 at m = 8 (b's
  8-term series is finite at z = 1e45 and −∞ at 1e50; 10^(308/(m−1)) = 1e44 is a lower bound;
  EA.18(c)), where a lane blend 0·∞ is NaN. Accuracy beyond π (θu/2) and of J⁻¹ near 2π (2πu/η) is
  the input's conditioning (EA.6, EA.10), not a domain. (b) SE(2): `SE2::log` (every X; jumps at |θ(X)| = π, PL.5) and
  `SE2::jr_inv`, `jl_inv` (|θ| < 2π; det = 4 sin²(θ/2)/θ², PL.11(d)). (c) The `Log` and `⊖` rows of
  §2.3: DLog(X) = J⁻¹(Log X) only for θ(X) < π (LG.13, LG.14), so `rminus_jacobians` and
  `lminus_jacobians` inherit it.
- **Options.** (a) Add the rows: Exp/`jr`/`jl` "φ·φ finite"; the SE(2) rows; `Log`, `rminus_jacobians`,
  `lminus_jacobians` "θ(X) < π". (b) State them in `# Domain` only, against R6.
- **Recommend.** (a). "φ·φ finite" is the scalar (`bool`-mask) domain, where only the chosen arm runs.
  A lane mask (0013, gated) evaluates the series arm too and is NaN from ‖φ‖ ≈ 1e22: its domain needs
  θ < 2π, or a safe argument on the series arm as NU.6(b) asks for r.
- **Blocks.** Every `# Domain` of Phase 3 (R6), the `debug_assert!`s, `jacobians_match_dual_*` domains.

**NU.10 SE(2) Jacobians, action rows and twin are not in the spec.**

- **Says.** §6: "$J_r$, $J_l$ and their inverses: Solà et al. 2018, Appendix (SE(2)), **permuted to
  rotation-first** in the Phase 3 PR that implements them". `PHASE3.md` §6: "the PR adds the
  rotation-first permuted SE(2) Jacobians to `NUMERICS.md` §6 before code".
- **Found.** Appendix C of arXiv:1812.01537v9 prints Ad, J_r, J_l and the action Jacobian,
  translation-first, and no inverse: the printed matrices equal the derived ones through the
  permutation (5.9e−109 for J, exact for Ad), but the rotation-first J⁻¹ is derived, not permuted
  (PL.9(d), PL.8). `so2-se2.md` §6 proposes Ad, ad, J_{l,r}, J⁻¹_{l,r} and the action Jacobians
  (PL.10). Also missing: §14's twin for a closed-form `SE2::jr_inv` (its dense inverse errs like
  κ₂(J) ≈ (‖ρ‖/2)², 2.5e7 at ‖ρ‖ = 1e4, so the tolerance must scale, PL.11(c)); §2.4's SO(2), SE(2)
  action rows, though `API.md` §3 lists `act_jacobians` for both; an SO(2) normalization (the Newton
  step is dimension-independent, PL.3(e)). The `se2_*` corpus could verify all of it today.
- **Options.** (a) Adopt `so2-se2.md` §6 into §6, §2.4, §12, §14 as the Phase 3 PR's first step, worded as derived
  and verified against the corpus (a formula change: edit and record). (b) Adopt Ad and J (permuted from
  Solà) and defer the inverses. (c) State an SO(2) `renormalize` with §3.6's bounds.
- **Recommend.** (a), with (c).
- **Blocks.** The SE(2) part of Phase 3 (`PHASE3.md` §6, §12); `se2_*` conformance.

**NU.11 SO(2), SE(2) `Log` at |θ| = π.**

- **Says.** §6: "$\mathrm{Log} = \mathrm{atan2}(s, c)$"; SE(2): "$\mathrm{Log}$: $\theta = \mathrm{atan2}(s,
  c)$, $\rho = V(\theta)^{-1} t$". `PHASE3.md` §6 says nothing of the branch.
- **Found.** At c = −1, `atan2` gives +π for s = +0 and −π for s = −0; SE(2) has exactly two preimages
  with |θ| < 2π at |θ(X)| = π, (π, −(π/2)Kt) and (−π, +(π/2)Kt), so ρ jumps by πKt with the sign bit
  (PL.3(c), PL.5). The corpus leaves z = (−1, ±0) unsampled (an mpmath value has no signed zero); §0.0
  lists it as undecided. The SO(3) rule of §3.2 does not extend to SO(2) by itself.
- **Options.** (a) Keep `atan2`'s sign-bit behaviour (a function of the representation, as §3.2), add
  both zeros as strata with the sign a generator input. (b) Canonicalize to (−π, π] with one extra
  select. (c) [−π, π).
- **Recommend.** (a): the tie is genuine and one rule for SO(2), SE(2) and SO(3) is cheaper to test.
- **Blocks.** `SO2::log`, `SE2::log` (`PHASE3.md` §6); the `so2_log`, `se2_log` strata.

**PH.3 The SE(2) kernel cannot return β from θ².**

- **Says.** `PHASE3.md` §3: "`se2_coeffs(θ²) -> (α, β)`", each call-site group "**one** `S::branch` over
  a tuple".
- **Found.** β = θ·a is odd in θ (PL.1, PL.2): a function of θ² cannot return it. The SE(2) Jacobians
  also need a, b (`jr_coeffs`), c (`jr_inv_coeff`) and α (PL.7, PL.8), and no listed group carries them.
  γ = 1 − θ²c has absolute accuracy only (≤ 5u) near |θ| = π where γ → 0 (PL.11(b)); a relative γ
  would need `jr_inv_coeff` to return it from its cot arm.
- **Options.** (a) `se2_coeffs(θ²) -> (α, a)`, β = θ·a formed by the caller (sign exact, PL.11(a)).
  (b) `se2_coeffs(θ)`, breaking the rule that kernels take θ². (c) Drop it; call `jr_coeffs`.
- **Recommend.** (a); Jacobian call sites reuse `jr_coeffs`, `jr_inv_coeff`; leave γ to its first
  componentwise consumer.
- **Blocks.** The Phase 3 kernel (`PHASE3.md` §3) with P1.5's ids.

**P1.9 The readings of §4.4, and the corpus size.**

- **Says.** §4.4: "`theta:1e-k`, $k = 12, \dots, 1$, and `theta:1e0` ($[1, \pi - 0.1)$): 64 samples
  each"; "`theta:dense`: 200 points per decade over $[10^{-4}, 1]$"; "`q:w0`: quaternions with $w = +0$
  exactly."; "`q:nonunit`: quaternions at $\|q\|^2 - 1 = \pm 2^{-45}$" (a purpose: "`Log`'s scale
  invariance, `from_*` normalization"); "Corpus v1 stays under 50 MB uncompressed; a family that would
  exceed it reduces its sample count by an edit to this list, never silently." The `rho:*` cells have
  no count (§4.4 has since gained "6 records" for SE_N(3) and the SO(2), SE(2) rows).
- **Found.** The tree's readings are in the table below. With 6 or 8 records a stratum's max
  under-samples its supremum, and the nearest-rank p99 is the max for every n < 100 (EA.11(d)): a defect
  touching a fraction f of the inputs is caught with probability 1 − (1 − f)ⁿ, 0.98 for f = ½ at n = 6.

| Left open by §4.4 | Reading in the tree |
|---|---|
| end of a decade | `theta:1e-k` = [10⁻ᵏ, 10⁻ᵏ⁺¹), as `theta:1e0` (EA.21 reads it so) |
| `theta:subnormal`, `theta:dense` | 64 log-uniform in [1e−310, 1e−309); 801 points 10^(j/200−4), j = 0…800, both ends |
| `theta:pi-1e-k` | scalar ids: the one value fl(π − 10⁻ᵏ); vector ids: 64 axes; π − θ is nominal to ≈ πu (2e−4 relative at k = 12, EA.19(a)) |
| `q:w0`, `q:nonunit` | (+0, u), u on S²; a Haar rotation scaled to ‖q‖² − 1 = ±2⁻⁴⁵ alternately. 64 records each for the `so3_*` ids (`so3_log`'s `q:nonunit` 128: each q and its negative; its `q:w0` has no negated half, P1.6), 6 for `sen3_ad` and `sen3_log` (12 for `sen3_log`'s `q:nonunit`). `coeff_r`'s `q:w0` is (n, +0), n = 1, 1e−3, 1e3: 3 records |
| who sees `q:*` | `so3_log`, `so3_act`, `so3_from_matrix`, `sen3_log`, `sen3_ad` (NU.2); act's p uniform on S² |
| `rho:*` cells | 6 records per stratum (`sen3_log` holds 12 in all but `q:w0`: each quaternion and its negative); 52 strata per file for the τ and matrix ids, 54 for `sen3_log` and `sen3_ad` (with `q:w0`, `q:nonunit`). The family is 18.4 MB; at 64 per stratum it would be 189 MB |
| `so2_*`, `se2_*` | both signs of each θ (3419 records); `se2_*` on the SE_N(3) strata without `theta:dense`, 8 records each, θ's sign alternating, ρ a direction on S¹; no non-unit z |
| streams (§2 item 4) | splitmix64 per stratum and purpose, seeded by SHA-256 of `<seed:016x>/<stratum>/<purpose>` (EA.23(d) R1–R3); seed `0x68656c69636f6964` |
| size, time | 31.8 of 50 MB; `just corpus` 17 CPU-minutes, `corpus-check` 3–5 minutes on 4 CPUs |

- **Options.** (a) Ratify and write the readings, counts included, into §4.4. (b) Also raise the
  vector-valued SE_N(3) ids (`sen3_exp`, `sen3_log`) to 64 records per stratum: +11.4 MB (exp 6.0, log
  5.5; per stratum, the committed bytes per record times (64 − count), summed over the six files),
  for 43.2 of the 50 MB and 6.8 MB left to every later id (Phase 2's `eig3`, `svd3`, `solve_cubic`;
  Sim(3), S², Γ; P1.1's `f32` companions, ≈ 2 MB).
- **Recommend.** (a). Counts in §4.4 make "reduces its sample count by an edit to this list" checkable.
  (b) is affordable but spends 11.4 of the 18.2 MB left before any later id has a size; the
  maintainer's call if the 0.98 above is judged too low.
- **Blocks.** Nothing; the counts bound what any stratum's max can prove.

**PH.5 Strided reads now panic in release: D11, R6 and §12 were widened.**

- **Says.** `PROJECT.md` D11 at `3da8149`: "The one documented panic class is out-of-bounds strided
  writes." `NUMERICS.md` §12: "Strided writes | in bounds | **panic** (the one documented class, D11)".
  `API.md` R6: "out-of-domain release behaviour is an unspecified value, never a panic, never UB (D11)".
  `PHASE2.md` §5: "Bounds are checked by safe slice indexing".
- **Found.** The strided-view PR changed D11 to "out-of-bounds strided access (`get`, `set`, `block`)",
  §12 to "Strided `get`, `set`, `block` | in bounds | **panic**" and R6 to "never a panic (bar
  out-of-bounds strided access, D11)": a read now panics in release too. That relaxes a hard rule of
  `CLAUDE.md` (no panic on valid input; a stated domain is a `debug_assert!`, never a release check),
  through an edit of a decision-log entry, not a record. The PR's reason (PH.2): slice indexing alone
  does not catch an access past a block's edge that lands in a neighbour inside the same slice, so the
  view checks row and column extents, on reads as on writes.
- **Options.** (a) Ratify: any out-of-bounds `get`, `set` or `block` panics in release, as edited.
  (b) D11 as written: only writes panic; an out-of-bounds read is a `debug_assert!` and returns the
  neighbouring element in release, silently. (c) Fallible reads (`get -> Option<&S>`), so only writes
  panic: a new public item.
- **Recommend.** (a), but as an explicit accept: it amends D11. A wrong read in release is a silent
  wrong answer, and the check is a comparison per access.
- **Blocks.** Nothing; `write_dense` (Phase 3) is the first consumer.

## Ratify the tree

The tree already holds these readings and edits and nothing waits on them. Each still needs an accept,
an amendment or a revert (*Open questions*).

**PH.4 Cholesky is specified by its PR.**

- **Says.** `PHASE2.md` §4: "`chol<S, const N: usize>(a: &Matrix<S, N, N>) -> (Matrix<S, N, N>,
  S::Mask)`: lower-triangular `L` and the positive-definiteness mask; `solve_lower`, `solve_upper`
  beside it." `API.md` R6 demanded a `debug_assert!` for every restricted domain.
- **Found.** The PR added `NUMERICS.md` §15 and edited §12 and R6: a pivot passes iff 0 < d_j and d_j is
  finite, an entry iff finite; a failed pivot stores L_jj = 1 and a zero column, a non-finite entry
  is stored as +0, so §12 has "`L` finite for every input"; `chol` asserts nothing (the R6 exception); the
  mask describes computed values, not a certificate about A; the solves read one triangle, take
  `b: Vector`, and `debug_assert!` a nonzero diagonal. No corpus id, stratum or twin exists for `chol`.
- **Options.** (a) Ratify §15 and the R6 exception. (b) Drop the finiteness promise: a failed pivot
  leaves `L` unspecified, and a consumer reads it only after the mask. (c) R6 as written: `chol`
  `debug_assert!`s positive definiteness and keeps the mask.
- **Recommend.** (a); decide with the Phase 2 ids whether `chol` needs strata (positive definite,
  near-singular, rank-deficient).
- **Blocks.** Nothing.

**P1.7 `mp.logm` is not a usable reference near π.**

- **Says.** §4.3, initial: "`mp.logm` of $R(q/\|q\|)$" for `so3_log`, "`mp.logm`-consistent
  quaternion" for `so3_from_matrix`, "`mp.logm`" for `sen3_log_n{1,2,3}`, "analogous" for `so2_*`,
  `se2_*`. §2 item 2, initial: "`Log` as a matrix logarithm".
- **Found.** `mp.logm` returns a non-real logarithm from an onset in (3.02, 3.03] for SO(3) and SE(2)
  and (3.02, 3.04] for SE_N(3), on mpmath 1.3.0 and 1.4.1 (SO.5, SE.3, PL.5; re-run here: `expm` of
  the hat matrix at 50 digits, two axes for SO(3) and two ρ for SE(2), both versions; max|Im| of `logm`
  is 0 at θ = 3.02 and 2.9 at 3.03 for SO(3), 0 and 3.1 for SE(2)), so it cannot serve any
  `theta:pi-1e-k` stratum nor the top of `theta:1e0`. For a non-orthogonal double matrix the
  "`mp.logm`-consistent quaternion" has no value (EA.9(d)). The generator uses Newton's method on the
  `Exp` definition, accepted by its residual and θ ≤ π, and `mp.logm` only as a test cross-check: the
  definition is kept, the tool replaced. The SO(3), SE_N(3) and SO(2)/SE(2) PRs edited §2 item 2 and
  §4.3 to say so, replaced "analogous" by nine rows, and read `so3_from_matrix` as the quaternion of
  the Frobenius-nearest rotation (the polar factor).
- **Options.** (a) Ratify the edited rows. (b) State only the definition (the preimage with θ ≤ π,
  by any method checked by its residual) and leave Newton to the generator's README. (c) Repair
  `mp.logm` by an eigendecomposition-based principal logarithm (agrees to 3e−69 at 80 digits, SO.5;
  untested at 120).
- **Recommend.** (a).
- **Blocks.** Nothing; it fixes what every log and `from_matrix` row of the corpus means.

**P1.8 The precision budget holds for θ ≥ 1e−12 only.**

- **Says.** §2 item 3, initial: "The worst catalogued cancellation (`e`, `NUMERICS.md` §4) loses about
  $4\log_{10}(1/\theta) + 3$ digits: 51 at $\theta = 10^{-12}$, leaving 69 against the 17 a binary64
  reference needs. A 1% sample is recomputed at `dps = 150` and must agree to 40 digits, or generation
  fails."
- **Found.** The bound holds (4L + 2.6; 71.1 digits kept at 1e−12) but §4.4 has strata below it. At
  `theta:subnormal` (≈ 1e−310) the loss is ≈ 1243 digits and a…e all evaluate to exactly 0 at 120
  digits; `theta:exact0` is 0/0; the 40-digit recheck of e fails from θ ≈ 1e−19 (39.2 digits at 1e−20),
  and of its `mp.diff` derivative in θ from 3e−18 (CO.17). The coefficient PR edited §2 item 3: 120 is
  the precision of the output; a cancelling definition is evaluated at guard digits, doubled until two
  evaluations agree, an exact 0 counting as disagreement (`coeff.stable`); at branch value 0 the
  extension is the definition at 2⁻²ᵖ. The recheck takes max(1, ⌈n/100⌉) records per stratum.
- **Options.** (a) Ratify. (b) Below a switch use the defining series as the reference (CO.17); a
  truncated series is the error class the recheck cannot see (EA.20(c)). (c) Narrow the strata to
  θ ≥ 1e−12, which D7 forbids.
- **Recommend.** (a).
- **Blocks.** Nothing.

**P1.10 Readings of `PHASE1.md` §3, 0007, 0003 and the decisions README.**

| Passage | Reading in the tree |
|---|---|
| §3: every `@generated` file "matches its generator" | `Cargo.lock` carries a foreign `@generated by Cargo` line: exempted by name |
| §3: the closure of each library crate "equals the allowed set" of 0007 | silent on an absent optional `mint` and on new members: `mint` is an optional direct dependency, never in `default`; a member with no budget entry is a violation |
| §3: "no crate other than `xtask` enables `__sweep`" | read to include forwarding (`__sweep = ["dep/__sweep"]`) and `default` reaching it; no crate has the feature yet, so only stubs test it |
| 0007 item 6: "`deny.toml` bans the list above for library crates" | cargo-deny cannot scope a ban to library crates; the closure lint is the per-crate guard; `just audit` runs bans without `xtask` and dev-dependencies |
| decisions README: no status table, banner or comment may cite a draft "as settled (`just lint` checks this)" | the lint reads status-table rows, Rust comments under `crates/`, `xtask/` and amendment banners; 0003 and 0007 name the draft 0013 in prose without saying so, which passes as worded |
| 0003 plan: "verified by the `compile_fail` doctest `real_has_no_partial_ord`" | rustdoc names doctests by item and line: two `compile_fail,E0369` doctests on `Real` |
| 0003 item 5: "its value path is bitwise identical to plain `S` evaluation"; README row: "`Dual` inherits the value path bitwise" | `PHASE2.md` §3 now excepts "the sign and payload of a NaN produced by arithmetic or a `libm` call" (PH.1); both texts lack the exception |

- **Options.** Per row, ratify the reading and amend the passage to say it, or hold the tool to the
  passage as written; rows 4 and 6 have no tool-side alternative (cargo-deny cannot scope a ban, rustdoc
  names a doctest by item and line).
- **Recommend.** Ratify the seven readings; correct 0003's plan, its item 5 and the README row while
  0003 is `ready` (afterwards a record is not edited to match code).
- **Blocks.** Nothing.

**PR.1 Edits to normative text made by implementing PRs.** Each edit below is a reading the PR could
not avoid and is in the tree now; it stands provisionally. This record's recommendations name the ones
that need a decision beyond ratification (P1.1–P1.9, NU.*, PH.*). Status-table rows (§0.0) are the
housekeeping every PR owes and are not listed.

| Document, section | Change | PR (commit subject) |
|---|---|---|
| `PHASE1.md` §2 item 2; §4.3 `so3_log`, `so3_from_matrix` | Newton reference; polar factor | "SO(3) corpus ids and rotation strata" |
| `PHASE1.md` §2 item 3; §4.2 `kind`; §4.3 `coeff_series` and the outputs paragraph | guard digits; series files | "corpus ids for the coefficient catalogue and its series" |
| `PHASE1.md` §4.2 identity; §4.3 `value`, `d_branch` | content hash; the derivative variable | "mpmath reference generator skeleton and coeff_k corpus" |
| `PHASE1.md` §4.3 `sen3_*`; §4.4 six records, `q:*` | references, inverse rule, counts; `sen3_log` and `sen3_ad` also see `q:w0`, `q:nonunit` (NU.2) | "SE_N(3) corpus ids (N = 1, 2, 3) and translation strata" |
| `PHASE1.md` §4.3 `so2_*`, `se2_*`; §4.4 SO(2), SE(2) strata | nine rows replace "analogous" | "SO(2) and SE(2) corpus ids" |
| `PHASE2.md` §3; `NUMERICS.md` §12 (`Dual` rows) | PH.1 | "forward-mode Dual<S, N> as a Real" |
| `NUMERICS.md` §12, §15; `API.md` R6, §2 | PH.4 | "fixed-size Cholesky with a positive-definiteness mask" |
| `PHASE2.md` §5; `NUMERICS.md` §12; `PROJECT.md` D11; `API.md` R6 | strided access, reads included, panics in release: D11 and R6 widened (PH.5); PH.2 | "strided read/write views over caller memory" |
| `PHASE2.md` §7; `API.md` §2 | `Point` at N = 2, 3 | "optional mint conversions" |

- **Options.** Ratify every row in one step; ratify or revert row by row (a revert redoes the PR's
  text, code or corpus).
- **Recommend.** Ratify or revert each row as its gap is answered; until then the corpus and the code
  are the readings, and no PR builds on a row marked "revert". The D11 row (PH.5) needs an explicit
  accept.
- **Blocks.** Nothing; it is the audit trail of what the tree assumed.

**PR.2 The maths pages describe text that has since changed.** They were written against the initial
specs: the open items on `mp.logm`, "`PHASE2.md` §3 states no rule for $\times$ and $\div$" and "the
nested case … has no named test" (EA.15, EA.17) describe passages that PR.1 has since edited.
**Options.** One resync PR after PR.1; a resync per page as its passages are decided; or a banner on each
page saying it was written against `3da8149`. **Recommend.** One PR resyncs the index and pages once PR.1
is decided, so a page is never ahead of the record; it also replaces the index's open items, which this
record now duplicates, by a pointer to the `ready` record. **Blocks.** Nothing.

**PR.3 Not measured, not decided.** `corpus-check` byte identity on macOS aarch64 and the aarch64
digests (`PHASE1.md` §11; EA.13) were not run (no such host), and the CI wall time of `corpus-check`
(3–5 minutes on 4 local CPUs; 60-minute timeout) is unmeasured on the runner. **Options.** Keep them as
`PHASE1.md` §11 items; or measure on the runner now. **Recommend.** Keep them as §11 items. **Blocks.**
`PHASE1.md` §11.

## Rationale

One record, not 29: the gaps share their evidence (the corpus and the pages), their owner and their
order, and answers depend on each other (P1.3, P1.4, P1.5, NU.5, NU.6 fix what the sweep is; NU.2, NU.3, NU.8
what `SO3` is; NU.10, NU.11, PH.3 what `SE2` is). A draft, not edits: formula, ordering, sign and
public-surface changes start as a record (`CONTRIBUTING.md`), and a draft cannot be cited as settled.
Recommendations, because several runs showed a cheaper option than the obvious one (NU.2's 511u,
P1.9's 11.4 MB). Lost: answering each gap in the PR that meets it, which is how PR.1 accumulated.

## Consequences

- No specification is edited by this record. Its PR adds the record's entries (`SUMMARY.md`, the
  decisions README, `CHANGELOG.md`) and, in the maths index, a pointer to it and the correction of one
  open item that said no record was filed. Until it is answered, the harness, sweep and Phase 3 PRs
  wait for the gaps that block them (the *Decisions* table's last column) or state the reading they take.
- A `ready` version changes the corpus contract: P1.1 (`f32` companions), P1.5 (three ids), NU.3
  (`w = −0` records), NU.6 (a `coeff_r` `w < 0` stratum), NU.11 (`z = (−1, ±0)`) and PH.1 (`real_*` ids)
  add ids or strata, and P1.9 (b) raises counts: `just corpus`, then `corpus-check`. P1.4 (a) keeps
  `coeff_r`'s `d_branch` as committed, and NU.2 (a) changes no byte; only NU.2 (d) removes strata, and
  needs a D7 exception. The bars do not change; NU.1 defines their unit.

## Implementation plan

Steps land in this order, each one PR with its verification, after the record is `ready`:

1. Ratification: apply or revert each PR.1 row and the accepted P1.2, P1.7–P1.10, PH.1, PH.2, PH.4 and
   PH.5 text (with 0003's item 5, plan and the README row, P1.10); update the §0.0 tables — verified by
   `just lint`, `just book`, and `just corpus-check` if a §4.4 count changed.
2. Instrument prerequisites: NU.1–NU.7, P1.1, P1.3–P1.6 and PH.1's `real_*` ids, as `NUMERICS.md` and
   `PHASE1.md` edits with the corpus changes they imply — verified by `just corpus-check` and the
   harness's own tests.
3. Phase 3 prerequisites: NU.8–NU.11 and PH.3 as `NUMERICS.md` and `PHASE3.md` edits, each with its
   stratum or twin row — verified by the corpus id that exercises it (`so3_log` `q:w0`, `se2_*`, `coeff_r`).
4. Resync the maths index and pages (PR.2) — verified by `just book`.

## Open questions

Each gap is open until the maintainer records accept, amend or reject beside its **Recommend.** line;
silence is not acceptance, the ratify gaps included. The decisions, in blocking order:

- [ ] NU.1: floors, `B` from the stored reference, double-double differences, the p99 method, the
  matrix sentence, φ and ρ scored apart (a)–(f)?
- [ ] P1.1: `f32` inputs (a), (b) or (c)? P1.2: ratify the schema and the manifest hash?
- [ ] NU.2: unit `q` in the domain of `act` with `q:nonunit` kept (a), `act` normalizes (b), or a D7
  exception to drop the stratum (d)?
- [ ] NU.3: sign-bit flip (a) or `w < 0` (b), and `w = −0` records in `q:w0`?
- [ ] P1.3: grid end per coefficient, `m ≤ 8`, coefficient objective? P1.4: the branch variable, `r`
  seeded on n²? P1.5: three new ids?
- [ ] NU.4: sum of the arms? NU.5: a switch per coefficient or per group? NU.6: `s = n²/w²` on the mask?
- [ ] PH.1: ratify `Dual`, add `real_*` ids? P1.6: window from `theta:1e-7`, "but `q:w0`"?
- [ ] PH.2: drop `Mat2`, add the `svd3` ordering, ratify the other four readings?
- [ ] NU.8: divide by the norm, canonical `w ≥ 0`? NU.9: add the §12 rows (scalar domain)?
- [ ] NU.10: adopt `so2-se2.md` §6 and an SO(2) `renormalize`? NU.11: keep `atan2`'s tie?
- [ ] PH.3: `se2_coeffs(θ²) -> (α, a)`? P1.9: 6 records per cell, or 64 for `sen3_exp`, `sen3_log`?
- [ ] PH.5: strided reads panic in release, amending D11 (a), or writes only (b)?

To accept, amend or revert (the tree already holds them):

- [ ] P1.7 the Newton reference; P1.8 guard digits; P1.10 the seven readings, with 0003 item 5 and the
  README row corrected; PH.4 Cholesky and the R6 exception.
- [ ] PR.1 row by row; PR.2 the resync after PR.1; PR.3 keep the unmeasured items in §11.
