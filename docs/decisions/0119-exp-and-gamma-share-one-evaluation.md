# 0119: Exp and Γ share one evaluation

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** the PR of its *Implementation plan*

## Context

A preintegration step over a constant rate and force needs three things: `ΔR = Exp φ`,
`Δv = Δt Γ₁(φ) f` and `Δp = Δt² Γ₂(φ) f` (GG.4). To linearise the step it also needs both
directional Jacobians. Today that takes three calls: `SO3::exp`, `gamma1(φ) * v` and
`gamma2(φ) * v`, or `exp` and `gamma_apply_jacobians::<1>` and `::<2>`. Each call forms its own `θ²`,
runs its own coefficient group (`exp_coeffs`, `jr_coeffs_sinc`, `gamma2_coeffs`) and builds its own
`W²`. The Jacobian pair also builds `φ×v`, `φ×(φ×v)` and `φ·v` twice. On the exact arm the three
groups take 3 roots, 1 `sin_cos` and 4 `sin` where 1, 1 and 1 suffice, and the Horner for `b`
runs twice.

**A single group would not be bit-identical.** `grouped` puts every member on its short prefix
below the *group's* second switch, which is the smallest of its members' (`0047` item 7). A member
moved into a larger group therefore takes its whole arm where its own group took the prefix. The
sweep admits a prefix only where it equals the whole arm "at every grid point and every corpus
record" (`thresholds::search::second`). The grid has 64 points per decade and the records are the
`coeff_*` ids', so that certifies a sample, not every float. Measured with one `grouped` over
`(k, cos θ/2, a, b, d)` against `exp_coeffs`, `jr_coeffs` and `gamma2_coeffs`, value and `d/dz` on
`Dual<_, 1>`:

| | `k` | `cos θ/2` | `a` | `b` | `d` |
|---|---|---|---|---|---|
| binary32, every float of `θ² ∈ [9e-4, 2e-3)` (9,905,269), `d/dz` bits | 28,558 | 99,468 | 0 | 27,912 | 22,395 |
| binary64, 5·10⁶ uniform `θ² ∈ [2.5e-3, 5.4e-3)`, value bits | 0 | 11 | 1 | 1 | 0 |
| binary64, the same draws, `d/dz` bits | 22,951 | 126,446 | 13,078 | 4,128 | 0 |

The `b` column compares against `jr_coeffs`. **The same split already exists between the shipped
groups.** `jr_coeffs`' `b` and `gamma2_coeffs`' `b` differ at the same 27,912 binary32 floats and
at 1 value and 4,128 `d/dz` of the binary64 draws, because their second switches are 2.738e-3 and
2.642e-3 (1.000e-3 and 1.382e-3 at binary32). `sinc_a_b_d`'s doc says its members are their own
groups' "to the bit", and its test checks a grid of ratio 1.0371, which these floats fall between.
The split is out of this record's scope: see *Further work*.

## Decision

1. **`so3::exp_gamma_apply(φ, v) -> (SO3<S>, Vec3<S>, Vec3<S>)`** returns `SO3::exp(φ)`,
   `gamma1(φ) * v` and `gamma2(φ) * v`, each to the bit. The `Γ_M v` here is the matrix applied by
   `Mul<Vector>`, as a caller writes it, not `gamma_apply_jacobians`' vector form.
2. **`so3::exp_gamma_apply_jacobians(φ, v) -> (SO3<S>, GammaJacobians<S>, GammaJacobians<S>)`**
   returns `SO3::exp(φ)`, `gamma_apply_jacobians::<1>(φ, v)` and `gamma_apply_jacobians::<2>(φ, v)`,
   each to the bit, `Γ₁`'s `0089` diagonal included. The two triples are separate tuple fields,
   not `[_; 2]`: an array indexed `m − 1` is an off-by-one a consumer cannot see (API.md §6 item 4).
3. **`so3::GammaJacobians<S> = (Vec3<S>, Mat3<S>, Mat3<S>)`** is a public alias for
   `(Γ_M v, ∂(Γ_M v)/∂φ, Γ_M)`, which `gamma_apply_jacobians` now names as its return type; the
   type itself is unchanged.
4. **`coeffs::exp_gamma_coeffs`**, with `exp_gamma_coeffs_dz` on `Dual<S, 1>`, evaluates `k` and
   `cos θ/2` as `exp_coeffs` does, `a`, `b₁`, `sin θ/θ` and `far` as `jr_coeffs_sinc` does, and
   `b₂` and `d` as `gamma2_coeffs` does. **Each member takes the arm its own group takes, by its
   own group's switches, in its own group's arithmetic**:
   - below the smallest of the three groups' second switches, every member takes its prefix and `b`
     is evaluated once;
   - between the smallest and the largest, each group's members take the prefix or the whole arm by
     that group's second switch. Here `b₁` and `b₂` can differ. This band lies below every first
     switch, so no exact arm runs in it (`exp_gamma_shorts_are_below_every_first_switch`);
   - above the largest, the evaluation is `grouped`'s: one Horner per member, and `b` once, while
     every member is on its series arm. Otherwise one exact arm runs at `select(all, 1, z)`, and each
     member chooses between it and its Horner by its own first switch. The divisions for `a` and
     `d` run only where one of `a`, `b`, `d` is on its exact arm, and `sin θ` only where `b` is,
     since `b` and `sin θ/θ` are its only readers. Where only `cos θ/2` is on its exact arm, this
     costs what `exp_coeffs` does. At binary64 `a` leaves its series at `θ² = 5.83` and `b` at
     `7.77`, so `θ ∈ [2.41, 2.79)` takes no `sin`.
5. `exp_gamma_coeffs_dz` inlines the body on `Dual<S, 1>` and keeps the value lanes of `k` and
   `cos θ/2`, which are the plain evaluation's bits (`0003` item 5). It never forms their
   derivative lanes.
6. Shared code, not repeated arithmetic: `SO3::exp`'s quaternion build is `so3::exp_from`, and
   `gamma_apply_jacobians`' body after its coefficients is `apply_jacobians_from` over
   `Operands`. Both new functions run those same bodies.

No formula, coefficient, switch or corpus record changes.

## Rationale

- **Per-group switches, not one group.** The table above is the reason. Bit identity holds by
  construction, so its tests are confirmation, not the proof: they take every switch and both of
  its neighbours, and a scan of the band. An exhaustive binary32 scan of the band (7,825,694 `θ²`,
  `measure_exp_gamma_band_every_f32`) found no difference in value or `d/dz`. With every slot
  forced onto the smallest switch, the kernel test, the binary32 test and both `Dual` subjects
  fail. The binary64 subject does not fail, since its value split is 1 in about 10⁶.
- **An explicit body, not a generic `grouped_as`.** A first version was `grouped_with` with a
  second switch per slot, evaluating a slot-indexed array inside the band. LLVM neither unrolled
  nor constant-folded it there: 650 instructions in the band against the three groups' 270
  (binary64). The explicit body returns scalars from every arm and captures by `move`
  (`.claude/rules/numerics.md`).
- **Forced inlining of the shared Jacobian body.** Left to LLVM, `apply_jacobians_from` was called
  out of line and returned its triple through memory, so the existing `gamma_apply_jacobians`
  took 105–155 instructions more per call (`exp + gamma_apply_jacobians::<1> + ::<2>`, 671 → 776
  at binary64, `θ = 0`). With `#[inline(always)]` it is `main`'s count to the instruction at
  binary64 and 3–6 fewer at binary32. The value-only separate calls are `main`'s to the instruction
  at both precisions.
- **`sin θ` behind `b`'s own mask.** Branching on it inside the exact arm saves 221 instructions
  per value call at binary64 `θ = 2.5` (89 more than without the branch; 131 more on the Jacobians). It
  costs up to 12 per call where every member is exact (`θ ≥ 2.79`), and that region is no wider.
  At binary32 `b` leaves its series first (`θ² = 8.98` against `a`'s `9.65`), so the branch is never
  taken there.
- **Why surface at all** (API.md §6 item 8). A consumer can call the three functions, so this
  removes no way to be wrong; it is a cost item. It returns the separate calls' bits, so moving
  from three calls to one changes no result.

**Measured** with callgrind instruction counts (`Ir`) per call. Each variant is
`#[inline(never)]` in a release binary at default flags, over 1000 calls at `φ = θ (0.48, 0.6,
−0.64)` and `v = (1, −2, 0.5)`. The separate calls are `main`'s (`f71a88d`), and `θ = 0.06` and
`0.04` fall in the binary64 and binary32 band:

| | θ | `exp`+`Γ₁v`+`Γ₂v` | fused | saved | `exp`+`gaj<1>`+`gaj<2>` | fused | saved |
|---|---|---|---|---|---|---|---|
| f64 | 0 | 241 | 187 | 54 (22%) | 671 | 559 | 112 (17%) |
| f64 | 0.06 | 273 | 247 | 26 (10%) | 866 | 683 | 183 (21%) |
| f64 | 0.3 | 289 | 235 | 54 (19%) | 882 | 691 | 191 (22%) |
| f64 | 2.35 | 410 | 371 | 39 (10%) | 1003 | 819 | 184 (18%) |
| f64 | 2.5 | 608 | 387 | 221 (36%) | 1229 | 799 | 430 (35%) |
| f64 | 2.9 | 760 | 441 | 319 (42%) | 1282 | 831 | 451 (35%) |
| f64 | 3.13 | 738 | 421 | 317 (43%) | 1260 | 811 | 449 (36%) |
| f32 | 0 | 262 | 206 | 56 (21%) | 573 | 545 | 28 (5%) |
| f32 | 0.04 | 298 | 234 | 64 (21%) | 653 | 612 | 41 (6%) |
| f32 | 0.3 | 310 | 244 | 66 (21%) | 665 | 624 | 41 (6%) |
| f32 | 2.9 | 373 | 318 | 55 (15%) | 728 | 694 | 34 (5%) |
| f32 | 3.13 | 511 | 323 | 188 (37%) | 957 | 713 | 244 (25%) |

The binary32 Jacobian rows save less because the three binary32 `Dual` groups cost little to begin
with (3-term prefixes): the fused group is 50 instructions at `θ = 0` against their 81. Wall-clock
rows (`so3/exp_gamma_apply*`, each beside `*_separate`) are in `benches/groups.rs`. They are owed
a quiet-host run of the bench gate, which is not recorded here.

## Consequences

- Neither function has a corpus id. Their promise is that no record of `so3_exp`, `so3_jl` or
  `so3_gamma2` moves. `xtask`'s `exp_gamma_apply_is_the_separate_calls_on_the_corpus` reads that
  promise on the committed `so3_gamma2` records (`so3_jl`'s strata and their `@f32` twins) and
  `so3_exp` records, at three `v` per record, plain and on `Dual`. `NUMERICS.md` §14 gains the row
  whose twin is the separate calls, to the bit.
- A change to any of the three groups' arms, switches or exact-arm arithmetic must be made in
  `exp_gamma_body` too. The arm lists are one definition each (`exp_arms`, `jr_arms`,
  `gamma2_arms`); the exact-arm arithmetic is pinned by the bit-identity tests, not shared.
- `apply_jacobians_from` is `#[inline(always)]`: the existing `gamma_apply_jacobians` depends on it
  for its instruction count.

## Implementation plan

1. Kernel, functions, tests, bench rows, docs. Verified by:
   - `coeffs::tests::exp_gamma_coeffs_are_their_groups_to_the_bit_{f64,f32}`;
   - `exp_gamma_shorts_are_below_every_first_switch`;
   - `exp_gamma_takes_one_root_and_two_sines`;
   - `exp_gamma_apply_matches_reference_f64` and its `f32`, `dual64` and `dual32` twins in `exp_gamma_tests`;
   - `xtask` `exp_gamma_apply_is_the_separate_calls_on_the_corpus`;
   - `measure_exp_gamma_band_every_f32` (ignored, about a minute under `ndebug`);
   - `just thresholds-check`, `just conformance` and `just determinism` unchanged.

## Open questions

None.

## Further work

- **The `b` split between `jr_coeffs` and `gamma2_coeffs`.** At about 3% of the binary32 `θ²` in
  `[1e-3, 1.382e-3)`, and at 1 value and 4,128 `d/dz` bits of 5·10⁶ binary64 draws in
  `[2.5e-3, 5.4e-3)`, `J_l`'s `b` and `Γ₂`'s `b` are different numbers. `sinc_a_b_d`'s "to the bit"
  claim fails at those floats. The fix belongs to the sweep: admit a prefix only where it is
  bit-identical on a scan dense enough to certify the band, or give each coefficient one second
  switch that every group takes, so a member's arm does not depend on its group. That needs its own
  record, and it moves corpus rows wherever the prefix now applies.
- `SEn3::exp` reads `Γ₁` as its translation column. An `SE_2(3)` step could reuse this evaluation
  if a consumer asks for one.
