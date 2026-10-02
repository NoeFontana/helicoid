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
(26 to 28); the `f32` sweep did for `r`'s branch variable, mask, domain and prior (29); the oracle runner did
for the file protocol and for what a runner owes (31, 32), and for
sophus-rs (33 to 35).

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
10. **The D12 prior for all eight.** §6 describes it as one switch for `a`, `b`, `c` and evaluates it
    as a named candidate; the sweep scores it for `k, a, b, c, d, e` and `cos θ/2` (`z < 0.01`,
    `θ < 0.1`) and, for `r`, as `(4 terms, 0.01)` in `s`, which is not D12's (29). `seeded:correct`
    no longer runs it: it runs the generated switches, and D12 is scored only by the sweep (and, in
    tests, by a kernel run at it). Report it for `k`, `d`, `e`, `cos θ/2`?
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
    are open, so it is swept only on the reading of 29). The seeded `Log` takes the series arm where `n² = 0`, the one
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
    `[0.1, 1)`, above D12's switch, so the defect and the D12 kernel the range was fitted against
    ran the same exact arm there (297 `u` both; the generated kernel that replaced D12 runs its
    series) and the two strata below it fit `p = 2.188`, 0.012 from the window's edge. The range
    was chosen after the fit. Which is meant: the strata above the plateau (as coded), the
    strata below D12's switch, or a range and window stated per precision?
27. **The `f32` "correct" kernel at the grid's edge** (§10: "the correct seeded kernel").
    `seeded:correct` ran the `tf_tree` D12 prior at `f32` until the `f32` sweep (`0016` item 3);
    it runs the sweep's constants now, as at binary64 (`subject_version` `generated`), and D12
    stays the sweep's prior. At `f32` each of `k`…`e` sits at the grid's top (`θ = 1`, its optimum
    above the grid: 9), so those switches are the grid's and not a measured optimum. The kernel's
    errors are printed and equal the CSV's, not gated (up to 4.4e3 `u`, `coeff_e`'s `d_branch` at
    `theta:1e0@f32`), and its `b` curve fits `p = -0.021` (`r² = 0.920`): silent. Is a kernel whose
    switches are the grid's edge the subject the `f32` half is silent on, or does the half wait for
    a grid that reaches them?
28. **The floor at `f32`** (question 1 at binary64). The code takes `2^-126`, the smallest normal
    binary32, as it takes `2^-1022` at binary64. No reference of an `@f32` stratum of the eight
    scalar ids is zero and the smallest is `1.6·10^-9` (`coeff_r`, `q:w0@f32`, `d_branch`), so the
    floor is not exercised today; it would read a subnormal by its quantization on the vector ids
    that a later record gives `@f32` strata. Is question 1's reading right at `f32`?
29. **The sweep of `r`** (`NUMERICS.md` §4 names `n²` as its branch variable and writes its series
    in `n²/w²`; `PHASE1.md` §0.0 leaves its branch variable, mask and `w ≤ 0` domain open; 0015
    (draft) NU.6 and P1.4 recommend a reading; no ready record specifies it, and `0016` item 1
    does not decide it). `xtask/src/thresholds` and `seeded/kernel.rs` sweep it as follows. The
    branch variable is `s = n²/w²`, `tan²(θ/2)`, formed by a division, the value the series arm's
    Horner takes, not NU.6's division-free `n² < s_sw·w²`. A candidate takes the series arm iff
    `w > 0` and `s < switch`; `w = 0` (`s` infinite) and `w < 0` take the exact arm. Only records
    with `w > 0` enter the objective (every `theta:*` record, no `q:w0` one); the corpus has no
    `w < 0` stratum, which NU.6 recommends, so that clause of the mask is not measured. The grid is
    the others' in `s`: 64 points per decade over `[1e-16, 1]` (`n/w` from `1e-8` to `1`), 1 to 8
    terms. `Dual<S, 1>` is seeded on `n²` at fixed `w`, the derivative the corpus stores (P1.4(a)),
    not on `s` (P1.4(c): a `w²` factor, regenerating `coeff_r`). The prior is `(4 terms, s < 0.01)`,
    the `z`-candidate of `k…e` taken in `s`: `θ = 0.2`, where D12's `θ < 0.1` is `s < 2.5e-3`, so
    its figures (1.1e10 `u` at `f64`, 229 at `f32`) are that candidate's and not D12's; `r` has no
    D12 value in `s`. Measured: 8 terms at `s = 1.07e-2` (`f64`, 260 `u`) and `0.143` (`f32`,
    26.7 `u`), both interior to the grid. Is this the reading: the division against the
    division-free mask, `s` against `n/w` or `θ` as the grid's variable, the seed, the prior, and
    the domain that leaves `w ≤ 0` unmeasured until a stratum exists?
30. **One maximum against per-stratum bars.** `PHASE1.md` §6 scores a sweep candidate by one
    maximum over all strata; §8 and D8/D15 compare per function, per stratum, on the max. They
    disagree where a candidate is best on the maximum and not on a stratum. Measured: `cos θ/2`'s
    objective is set by `theta:1e0`, where both arms are one expression, so the tie goes to two
    terms below `θ = 7.5e-8` (`z < 2.6e-6` at `f32`), and the exact arm above costs 1.4 to 1.9 `u`
    where the D12 series costs 0.85 to 1.7, on nine strata (`theta:1e-3`…`1e-8` at `f64`,
    `theta:1e-2`, `theta:1e-3`, `theta:dense` at `f32`). An oracle that behaves like that series
    would dominate `helicoid` there. Is the objective per stratum (a candidate must not lose to the
    prior on any stratum), or one maximum with the tie broken toward more terms, or is the
    measured result kept and the stratum bar of §8 the one to relax? Nothing is decided; the
    subject reports it and a test pins it.
31. **The oracle file protocol** (§7 says "reads corpus JSONL … writes `{"id":…,"out":{…}}` with
    hex-float outputs"). `xtask/src/conformance/oracle.rs` reads it as: `<runner> --version`
    prints the pin on one line, the rows' `subject_version`; `<runner> --out DIR FILE.jsonl…`
    writes `DIR/<fn>.jsonl` for each file whose function it answers; an answer's `id` is its line's,
    a scalar may be a bare string, `nan`, `inf` and `-inf` spell a non-finite value; a missing,
    extra or misnumbered record, a file of no corpus function and an owed function (`RUNNERS`) with
    no file are errors; a non-finite answer is a recorded row, not a failed run (§7: "oracles may
    be wrong"), where 6 fails an in-process subject on it. The three container runners will
    implement this in C++. Amend §7 with it, or take another shape?
32. **What `tf_tree_math` answers beyond §7's list** (exp, log, `a, b, c` via `V`/`V⁻¹`, `slerp`,
    ScLERP). At the pin the runner answers `so3_exp`, `so3_log`, `sen3_exp_n1` and `sen3_log_n1`.
    `coeff_*` and the Jacobians are not exported (`v_coeffs`, `vinv_c3` are private). Three exported
    functions have corpus ids and are not answered: `Iso3::adjoint` (`sen3_ad_n1`, `Ad(T)` on the
    six basis twists), `Quat::rotate` (`so3_act`, which assumes a unit `q`, so `q:nonunit` and
    `q:w0` would score that assumption) and `quat_from_rot3` (`so3_from_matrix`, checked for shape
    and finiteness only). For those ids `0010`'s "dominate on every paired stratum" pairs nothing.
    Answer them, with which strata, or leave the ids to sophus-rs and the containers? sophus-rs
    0.15.0 exports `LieGroup::adj` and `LieGroup::transform` too, and its runner leaves
    `sen3_ad_n1` and `so3_act` unanswered until this is answered.
33. **Which sophus-rs is oracle #2** (§7: "pinned version"; its example `sophus_to_helicoid_tangent`,
    `0002` and `PROJECT.md` §1 read Sophus as translation-first). `sophus_lie` 0.10.0 to 0.14.0
    are: tangent `[υ; ω]`, parameters `[p; q]`, no Jacobian. 0.15.0 (2025-08-09) is rotation-first,
    tangent `[ω; ν]`, parameters `[q; p]` with `w` first, and has `left_jacobian` and
    `inv_left_jacobian`. The runner pins 0.15.0, the latest, the only one with a Jacobian and the
    convention of `0002`, so its conversions change a layout, not an order, and it has no
    permutation and no permutation test; a bump fails its tests.
    Is 0.15.0 the oracle, is 0.14.0 a second runner (the permutation §7 has in mind), or does the
    version `locus_fusion` builds against decide? Say which `0002`'s "Sophus is translation-first"
    means, C++ Sophus and sophus-rs 0.14 or all of them.
34. **Right Jacobians of an oracle that has only the left one** (§7: "any exposed Jacobians").
    sophus-rs exposes `J_l` and its inverse; the runner answers `so3_jr`, `so3_jr_inv`,
    `sen3_jr_n1` and `sen3_jr_inv_n1` at the negated tangent, `J_r(τ) = J_l(−τ)` (`NUMERICS.md`
    §1), which repeats the left rows' errors under four more ids. Keep the four, or only the ids
    the oracle has?
35. **`sen3_log` at `w = +0` scores an oracle's branch.** sophus-rs 0.15.0 returns `−π û` for
    `q = (+0, u)` where `NUMERICS.md` §3.2 says `+π û`; `Exp(−π û)` is `−q`, the same rotation.
    `so3_log` aligns the sign at `w = +0`, so its row does not show it, while `sen3_log_n1`
    (`SignRule::Fixed`) scores 1.78e16 `u` in `q:w0`, since `ρ = J_l⁻¹(φ) x` follows the branch of
    `φ`. Align `sen3_log` there too (to which `ρ`?), or keep recording it as an error?
