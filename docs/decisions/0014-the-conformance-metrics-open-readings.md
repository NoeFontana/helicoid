# 0014: The conformance metric's open readings

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** none; a draft authorises nothing.

## Context

`NUMERICS.md` §11 defines the forward error and is silent on several things a harness must decide
before it can score a record. `xtask/src/conformance/metric.rs` (`TABLE` and its module docs) took
the smallest reading of each so that the harness core runs (`PHASE1.md` §5). The readings fix the
bar every envelope baseline is blessed against (`0006`); none is blessed yet (`PHASE1.md` §8 is not
started), so each is still free to change. `xtask/src/thresholds` (`PHASE1.md` §6) did the same for
the threshold sweep (questions 8 to 11) and for the file it writes, `xtask/src/seeded/generated.rs`
(12, 13); the self-test (`PHASE1.md` §10) did for the planted `c` (14, 15) and for the seeded SO(3)
subject and its `Log` defects (16 to 20), and for the seeded SE_N(3) subject and its defects (21
to 25); `f32` scoring and the self-test's `f32` half did for the curve, its kernel and its floor
(26 to 28).

## Decision

None yet. The record lists the readings so they are argued before the first baseline; each one
that survives is a `NUMERICS.md` §11 edit.

## Rationale

—

## Consequences

Until `ready`: the readings are provisional, `just envelope --bless` does not exist, and no baseline
is committed. A changed reading changes every `max_u` under it.

## Implementation plan

None until the open questions are resolved.

## Open questions

1. **Floor of a tangent or a coefficient.** The code uses `2^-1022`, so `theta:1e-k` reads a
   relative error (`EA.4(b)`, which argues it for `Log`'s tangent and for subnormals only). Is one
   floor right for the coefficients too?
2. **SE_N(3) tangent at `w = +0`.** The code gives it no sign forgiveness: `−τ` is not the logarithm
   the input's sign selects. `(−φ, J_l(−φ)⁻¹x)` is also a logarithm of the same `X`. Should §11's
   sign invariance cover it?
3. **One score per record.** The code takes the maximum over a record's output fields. Per-field
   rows would say which field scored.
4. **Whole-tangent norm.** §11 says vectors take the 2-norm, so a tangent's rotation and translation
   blocks share one norm: at `rho:1e4` a rotation error many `u` of `φ` reads under 1. Score the
   blocks apart?
5. **The ρ-scale.** The code takes the 2-norm of the input tangent's translation entries (all `N`
   blocks together) for `sen3_exp.x` and `se2_exp.t`. §11 says "‖ρ‖-scale".
6. **Overflow.** An `F` past binary64 (a gross error against a zero reference) counts as non-finite
   and fails the run. Is that the intended gate, or a failing score of its own?
7. **The field of the `b` defect's curve** (`PHASE1.md` §10, not §11). The self-test fits the `value`
   field's per-stratum `max_u` over `theta:1e-8`…`theta:1e-2` (`p` = 1.937). Through `Dual` the
   exact arm's derivative errs by θ⁻⁴, so a record's maximum over `value` and `d_branch`, which is
   what a results row holds, fits `p` = 3.98 and misses the window. Score the field, or widen the
   window and score the record?
8. **Per decade of what** (`PHASE1.md` §6). "A log grid of 64 points per decade" of the branch
   variable, "spanning θ ∈ [1e-8, 1]". The sweep takes decades of `z = θ²`: 1025 points,
   `10^((i−1024)/64)`, 128 per decade of θ. Decades of θ are 513 points, the even indices; they move
   `a` to θ = 0.835 (37.8 u against 25.8), `b` to the top (64.4 u against 63.5) and `e`'s switch to
   0.965 (the same objective), and leave `k`, `c`, `d`. Which is meant?
9. **The top of the grid.** It ends at `z = 1` and a candidate takes the series arm where
   `z < switch`, so the corpus record at θ = 1 exactly is always on the exact arm and no switch
   above θ = 1 is tried. That record is `e`'s objective (8.69e3 u, 7 terms); a top of `nextUp(1.0)`
   gives `e` 8 terms and 2.60e3 u, and leaves `k`, `a`, `b`, `c`, `d`. Is the range closed at 1, or
   does the grid extend to the first switch that puts that record on the series arm?
10. **The D12 prior for all six.** §6 describes it as one switch for `a`, `b`, `c` and evaluates it
    as a named candidate; the sweep scores it for all six coefficients. `seeded:correct` no longer
    runs it: it runs the generated switches, and D12 is scored only by the sweep (and, in tests,
    by a kernel run at it). Report it for `k`, `d`, `e`?
11. **Ties.** Exact `f64` equality of the objective (`tied` counts them). A tolerance would make
    near-equal candidates tie, and the cheaper one win. Exact, or a tolerance?
12. **The source revision of a generated file** (`PHASE1.md` §6: `Source sweep rev: <git rev>`).
    A git revision cannot name the commit that holds the file and differs on every commit, so no
    byte-for-byte drift check survives it. The emitter writes the SHA-256 of the sweep CSV and, on
    a line of its own, of `coeff_series.jsonl`, which supplies the literals. Amend §6 to say so, or
    is a revision meant?
13. **One objective line per constant.** §6 shows `Objective (max u): value 1.9, derivative 3.4.`
    on the revision's line, for one constant; the file holds six. The emitter writes one line
    above each constant, from the CSV's `value_max_u` and `deriv_max_u`. Which is meant?
14. **What `10^-8` is** (`PHASE1.md` §10, the planted `c`). §6 puts a switch in the branch
    variable, where `10^-8` is `z` and `θ = 1e-4`; `PROJECT.md` §2 names "`1e-8` in most
    libraries" without saying which variable. The planted `c` reads it as `θ`: `z = 1e-16`, the
    sweep grid's first point. Which is meant?
15. **"Ranks it dominated"** (§10). No margin is stated. The self-test says dominated when the
    candidate's objective exceeds the chosen candidate's by more than `10^6` (7.3e33 u against
    7.4e3 u, rank 7841 of 8200), and ranks the `(terms, switch)` the subject declares, which it
    checks against the objective the harness measures for the subject's `c`, to the bit. State a
    margin, a rank, or neither?
16. **The correct kernel's bar** (§10 names none for it). The self-test holds `seeded:correct`'s
    `Exp` and `Log` to 4 `u` on every stratum's max: 1.4 times the measured 2.91 `u`, a few
    roundings (`docs/maths/so3.md` SO.6, `docs/maths/coefficients.md` CO.16). `0006` blesses a
    baseline against a bar. State one, or leave it to the envelope?
17. **Which `acos` strata reach `10^7 u`** (§10: "`theta:1e-k` and `theta:pi-1e-k` max ≥ `10^7 u`").
    Not every stratum can: `θ` from `χ = cos θ` errs by `≈ u/(2θ²)` (`so3.md` SO.6(a)). The
    self-test gates every stratum of `theta:1e-k` from `k = 4` and of `theta:pi-1e-k` from `k = 7`,
    the boundaries as measured (8.2e7 `u` at `theta:1e-4`, 1.5e7 at `theta:pi-1e-7`), and prints
    the family maxima (9.0e15, 8.5e7). Every stratum from those `k`, the family max alone, or a
    stated `k`?
18. **`Log` without the flip: "every `so3_log` stratum"** (§10). Where every negated record has
    `n² = 0` (`theta:exact0`, `theta:subnormal`), `r`'s series arm `2/w` is analytic in the sign of
    `w`, so no defect can fail there. The self-test excuses a stratum only when that holds of its
    records and the defect answers it bit for bit as the correct kernel; every other negated half
    must fail (at least `10^7 u`; the missing flip is `2π` in `φ`), and one must. Amend §10 for the
    two strata, or plant a defect that survives the arm?
19. **`cos(θ/2)` in the seeded `Exp`.** §3.1 has the series arm compute it from its own series and
    one `sin_cos` in the exact arm; §4 evaluates `(k, cos θ/2)` in one `branch`. No series of
    `cos(θ/2)` is committed (`PHASE1.md` §0.0, *Missing*), so `Exp` takes `k` from `coefficient`
    alone and `cos(θ/2)` by the exact arm at every `z`, on a second `sin_cos`; the `sqrt(z)` under
    it has an infinite `Dual` derivative at `z = 0`, so `Exp` is a value subject. Acceptable until
    the series lands, or is the group rule owed first?
20. **The switch of `r`** (§4: generated; `PHASE1.md` §0.0: its branch variable and `w ≤ 0` domain
    are open, so it is not swept). The seeded `Log` takes the series arm where `n² = 0`, the one
    place the exact arm is `0/0`, and types no constant: elsewhere `atan2` and the division share
    the one `n`, so `Log` is scale-invariant (§3.2) to the bits a subnormal `n²` keeps, where a
    bound at the smallest normal number returns 1.0926 for the rotation of angle 1 scaled to
    `1.4·10^-154`. Is `0 < n²` the intended switch, or `s = n²/w²`?
21. **The SE_N(3) correct kernel's bar** (§10 names none). The self-test holds `seeded:correct`'s
    `Exp`, `J_r` and `J_l` to 12 `u` on every stratum's max of the nine
    `sen3_{exp,jr,jl}_n{1,2,3}` ids: 1.35 times the measured 8.90 (`J_l`, `rho:1e4/theta=pi-1e-6`;
    `Exp` at most 3.47), one bar for the nine as 16 has one for two. `Q` is several products deeper
    than SO(3)'s `Exp`, so it is not 16's 4. The maxima are of the per-coefficient kernel (24). State
    one bar, one per id, or leave it to the envelope?
22. **When an SE(3) defect "fails"** (§10: "every `rho:*` stratum fails", "`sen3_jr*` fails"). No
    magnitude is stated. The self-test takes the line of the SO(3) defects (17, 18): a stratum
    whose max reaches `10^7 u`, a million times 21's bar, or a non-finite output; a `NaN` or an
    unscored stratum does not fail, so it cannot fire a mechanism. The weakest stratum of the
    Jacobians is `rho:1e-6/theta=1e-8`, whose block is `5·10^-7` against §11's floor of 1, and it
    still reaches 5.2e9 `u`. Is `10^7 u` the line, a multiple of the bar, or any stratum over it?
23. **What "every stratum" counts** (§10). "Every `rho:*` stratum" is the 25 cells
    `rho:<scale>/theta=<θ>` (§4.4: 5 scales times 5 angles) of each of `sen3_exp_n1`, `_n2` and
    `_n3` (75, all three failing at least 8.9e15 `u`), translation-first being `[ρ₁; …; ρ_N; φ]` at
    every `N` (`docs/maths/se3.md` SE.14(c)); the `theta:*` strata are not read. "`sen3_jr*`
    fails" is every one of the 52 strata (27 `theta:*` without `theta:dense`, and the 25 cells) of
    every `sen3_jr_n<N>` and `sen3_jl_n<N>`, `Q` being in both (312 of 312); a file with fewer
    strata than that fires nothing. Must a mechanism have every stratum of a whole file, or a
    stated subset?
24. **Per-coefficient switches against §4's group rule.** §4 says a call site never evaluates one
    coefficient alone: `jr_coeffs` is (a, b) and `q_coeffs` is (b, d, e) in one `branch`. The
    sweep generates one switch per coefficient (`z < 0.7234` for `a`, `0.9647` for `b` and `e`, `1`
    for `d`), and the seeded `Exp` and `Q` run each at its own, so on `z` in `[0.9647, 1)` one `Q`
    has `b` and `e` on the exact arm and `d` on the series arm. The `helicoid` crate will not run
    that. Is a switch per call-site group owed to the sweep before 21's bar is measured, or is the
    per-coefficient kernel the subject until the `helicoid` one exists?
25. **The association of `Q`'s words** (§5.3 states the words, not an order of rounding). The
    subject forms 3×3 products, entry `(r, c)` as `(t₀ + t₁) + t₂`, each word associated left to
    right as written (`(φ^ρ^)φ^`, `(φ^φ^)ρ^`, `(φ^ρ^φ^)φ^`), and `W²` in `J_l(φ) = I + aW + bW²` as
    `W·W`, not `φφᵀ − θ²I`. Another association moves the max by a few roundings, and 21's bar
    reads it. State the association, or leave it to the `helicoid` subject?
26. **The `f32` curve** (`PHASE1.md` §10 fixes `theta:1e-8`…`1e-2` and `p ∈ [1.8, 2.2]`, at no
    precision). The `f32` half fits the value of `coeff_b` over `theta:1e-3@f32`…`theta:1e-1@f32`,
    the strata wholly above `√(6u) = 6·10^-4` at `u = 2^-24`: `p = 2.134`, against 1.937 at
    binary64. `theta:1e-4@f32` holds `6·10^-4`, where `b` by its definition is 0 and the error a
    plateau of `1/u`; with it the four fit `p = 1.644`, outside a window that does not depend on
    `u`, whereas binary64's range keeps its plateau stratum (one of seven). The top stratum is
    `[0.1, 1)`, above D12's switch, so the defect and the kernel run the same exact arm there
    (297 `u` both) and the two strata below it fit `p = 2.188`, 0.012 from the window's edge. The
    range was chosen after the fit. Which is meant: the strata above the plateau (as coded), the
    strata below D12's switch, or a range and window stated per precision?
27. **D12 as the `f32` "correct" kernel** (§10: "the correct seeded kernel"). Until the `f32`
    sweep (`0016` item 3), `seeded:correct` runs the `tf_tree` D12 prior at `f32`, four terms
    below `z = 0.01` for all six coefficients (`NUMERICS.md` §4 lists it for `a`, `b`, `c` only),
    its `subject_version` `d12`. Its errors on the `@f32` strata are large (up to 8.1e9 `u`,
    `coeff_e`'s `d_branch` at `theta:1e-1@f32`; 1.9e7 for `c`), printed and pinned by a test, not
    gated, and its `b` curve fits `p = -1.200` (`r² = 0.753`): silent because it rises with `θ`.
    Is a kernel that is not the correct one the subject the `f32` half is silent on, or does the
    half wait for the sweep's switches?
28. **The floor at `f32`** (question 1 at binary64). The code takes `2^-126`, the smallest normal
    binary32, as it takes `2^-1022` at binary64. No reference of an `@f32` stratum of the eight
    scalar ids is zero and the smallest is `1.6·10^-9` (`coeff_r`, `q:w0@f32`, `d_branch`), so the
    floor is not exercised today; it would read a subnormal by its quantization on the vector ids
    that a later record gives `@f32` strata. Is question 1's reading right at `f32`?
