# Reference generator

Computes the conformance corpus from **definitions** with mpmath at 120 digits, never from the
closed forms in `docs/NUMERICS.md` (`docs/PHASE1.md` §2, §4; `docs/decisions/0006`). Inputs are exact
binary64 values (exactly binary32 ones in the `@f32` strata below); the corpus in `../corpus/` is
committed and regenerates byte-identically.

```
just corpus         # regenerate ../corpus, deleting a stale *.jsonl   (python -m gen all --out DIR [--jobs N])
just corpus-check   # regenerate into a temp dir, compare every file byte for byte, run the tests
just corpus-test    # unit tests only (stdlib unittest)
uv run --frozen python -m gen list
```

## Files

| Path | Role |
|---|---|
| `pyproject.toml`, `uv.lock`, `.python-version` | CPython 3.12, `mpmath` pinned exactly; no numpy, no gmpy2 (`setup` refuses a non-Python backend) |
| `gen/precision.py` | the only place `mp.dps` is set (120; rechecks at 150); `to_f64`, `to_f32` |
| `gen/rng.py` | splitmix64, per-stratum streams, log-uniform, uniform on S² |
| `gen/strata.py` | the strata; `gen/registry.py` the function ids; `gen/coeff.py`, `gen/so3.py`, `gen/sen3.py`, `gen/so2.py`, `gen/se2.py` their definitions |
| `gen/series.py`, `gen/check.py`, `gen/check_sen3.py`, `gen/check_se2.py` | the exact Taylor series (`coeff_series`); the per-record cross-checks |
| `gen/fmt.py`, `gen/corpus.py`, `gen/manifest.py` | text formats and `Mat`, assembly (parallel by stratum) and recheck, `MANIFEST.json` |

## Record schema

One JSONL file per function id, one record per line, compact JSON with sorted keys, LF newlines:

```json
{"id":832,"in":{"theta":"0x0.0p+0"},"out":{"d_branch":"-2.08333333333333333333333333333e-2","value":"5.00000000000000000000000000000e-1"},"stratum":"theta:exact0"}
```

- `id`: 0-based line number in the file. Strata are emitted in the order of their catalogue, so a
  stratum appended to the end never renumbers the rest.
- `in`: exact binary64 as `float.hex()`. `out`: decimal strings, 30 significant digits, `d.ddd…e<n>`
  (no `+`, no zero padding), rounded half-to-even from the 120-digit value in integer arithmetic.
- Vectors are arrays. A matrix (`fmt.Mat`) is its column-major flat array plus a sibling
  `"shape":[rows,cols]` in the same `in` or `out` object, so an object holds at most one:
  `"out":{"J":[…9 entries…],"shape":[3,3]}`.
- `MANIFEST.json`: generator identity, Python (major.minor) and mpmath versions, `dps`, seed, recheck
  parameters, and per file the SHA-256, record count, `kind` and one check count: `rechecked`
  (records recomputed at 150 digits) for `kind: "corpus"`, the files of the schema above;
  `verified` (series equal to exact algebra) for `kind: "series"`, i.e. `coeff_series`, which has
  its own record schema below. A reader of function-id files skips every other `kind`.

## Strata (`SCALAR_THETA_STRATA`, `QUAT_STRATA`)

| Stratum | θ values (binary64) | Records |
|---|---|---|
| `theta:1e-k`, k = 12…1 | log-uniform in the decade [10⁻ᵏ, 10⁻ᵏ⁺¹) | 64 each |
| `theta:1e0` | log-uniform in [1, π − 0.1) | 64 |
| `theta:exact0` | 0 | 1 |
| `theta:subnormal` | log-uniform in [10⁻³¹⁰, 10⁻³⁰⁹), below the smallest normal | 64 |
| `theta:pi-1e-k`, k = 1…12 | the single value fl(π − 10⁻ᵏ) | 1 each |
| `theta:dense` | 10^(j/200 − 4), j = 0…800, rounded: a grid, no randomness | 801 |
| `<S>@f32`, S each of the above | S's values rounded to nearest-even binary32 (`theta:subnormal@f32`: below) | as S |

`@f32` strata (`docs/decisions/0016`) belong to the coefficient ids (`coeff_k`, `a`…`e`, `cos_half`,
`r`) and follow every stratum of the ids they belong to, in the order of the strata they twin: no
id, stratum name or record of the binary64 corpus moves. Each input is exactly a binary32, so a
binary32 subject receives it by a lossless cast, and the reference is the function at that rounded
input, at 120 digits (`to_f32` is integer arithmetic, and `corpus.require_binary32` asserts every
input round-trips through the platform's binary32 at generation). Readings:

- `theta:subnormal@f32` is log-uniform in [10⁻⁴⁰, 10⁻³⁹), rounded, from its own stream: rounding
  the binary64 decade [10⁻³¹⁰, 10⁻³⁰⁹) gives 64 zeros. Each decade lies below the smallest normal
  of its type, where θ² underflows and 1/θ overflows.
- `theta:pi-1e-k@f32`: fl(π − 10⁻ᵏ) rounded, for k = 8…12 the one binary32 `0x1.921fb6p+1`, which is
  above π by 8.7·10⁻⁸ (cos(θ/2) < 0 there). The records are kept, none deduplicated: five strata of
  the same input and reference in every θ id; `coeff_r` rounds n and w, which differ.
- `theta:exact0@f32` is `theta:exact0`'s record under its own name. A decade's rounded values may
  lie half a binary32 unit outside it.
- `coeff_r`: n and w of the binary64 quaternion each rounded (`q:w0@f32`: n ∈ {1, fl32(10⁻³), 10³},
  w = +0). `theta:subnormal@f32` is the quaternion of that stratum's unrounded draw, n and w each
  rounded, not of the rounded θ the θ ids score at the same index: n is half that θ in 29 of 64
  records.

A scalar id sees each fixed-θ stratum once. A vector id (`so3_*`) gives every θ an axis
(`Stratum.samples`), uniform on S² (z uniform, then φ uniform, from the stratum's `axis` stream):
one axis per random or dense θ, 64 for a fixed θ (`theta:pi-1e-k`), and (1, 0, 0) for
`theta:exact0`, where nothing reads it. `so3_exp` and the four Jacobians take φ = fl(θ · axis)
componentwise (‖φ‖ = θ(1 + O(u))): 2466 records. The quaternion ids (`so3_log`, `so3_act`,
`so3_from_matrix`) take the unit quaternion (cos θ/2, sin θ/2 · axis) rounded componentwise, and two
strata of their own, after the θ strata (`Stratum.quaternions`):

| Stratum | Quaternions | Samples |
|---|---|---|
| `q:w0` | (+0, u), u a unit vector, an angle of exactly π (`coeff_r`'s `q:w0` is its own, below) | 64 |
| `q:nonunit` | a Haar-random unit quaternion (Shoemake), scaled to ‖q‖² − 1 = +2⁻⁴⁵ (even samples) or −2⁻⁴⁵ (odd), rounded: the norm is that to ~1e-16 | 64 |

## SO(3)

| Id | In | Out | Reference, from the definition |
|---|---|---|---|
| `so3_exp` | `phi` | `q` | Σ pⁿ/n!, p = (0, φ/2), by Hamilton products |
| `so3_jr`, `so3_jl` | `phi` | `J` (3×3) | Σ (∓W)ⁿ/(n+1)!, W = [φ]× |
| `so3_jr_inv`, `so3_jl_inv` | `phi` | `J` | `mp.inverse` of the above |
| `so3_log` | `q` | `phi` | the φ with Exp(φ) = q/‖q‖ on θ ∈ [0, π], by Newton's method on the series |
| `so3_act` | `q`, `p` | `Rp` | R(q/‖q‖) p, R(q) of `docs/NUMERICS.md` §1; p uniform on S² |
| `so3_from_matrix` | `R` (3×3) | `q` | the quaternion of the rotation nearest to R in Frobenius norm |

- **Series** stop when two terms in a row are below 10^-(dps−10) of the partial sum, entry by entry
  (each part of the quaternion): an entry of 10⁻³¹¹ keeps its 30 digits, as at `theta:subnormal`.
  A zero component of φ (of vec q for `log`) stays exactly zero in the result: vec Exp(φ) ∥ φ.
- **`so3_log`** solves Exp(φ) = q̂ by φ += J_r(φ)⁻¹ · 2 vec(Exp(φ)* q̂) at precisions that follow
  the accuracy reached (float seed, 14 digits → 120), and stops when the residual, at the full
  precision, is below 10⁻¹¹⁰ of |vec q̂|. The series is the definition, so the float `atan2` that
  seeds it cannot bias the result, and the form of `docs/NUMERICS.md` §3.2 is never evaluated. q̂ is q/‖q‖ flipped
  to w ≥ 0 (§3.2: w < 0 negates, w = +0 stays): at w = +0, (+0, u) gives +πû and (+0, −u) gives
  −πû, the function of the quaternion, not of the rotation. Every stratum but `q:w0` holds each
  quaternion and then its negative (the half a `Log` without the flip gets wrong; `PHASE1.md` §10);
  the negative of (+0, u) has w = −0, not the +0 of `q:w0`, so that stratum has none.
- **`so3_from_matrix`**: R is R(q) of the stratum's quaternion rounded to binary64, so not exactly
  orthogonal, and a scaled rotation for `q:nonunit`. The reference is the polar factor Q of the
  exact input (Higham's iteration X ← (X + X⁻ᵀ)/2, the nearest rotation in Frobenius norm), then q
  with R(q) = Q by q ← q(1, d/2) with d = vee skew(R(q)ᵀQ) (cubic), sign w > 0, else the first
  nonzero of x, y, z positive. A rotation by exactly π (`q:w0` rounds to a symmetric matrix) has
  w = 0 and either sign; |w| below 10⁻¹¹⁰ is that zero. `docs/NUMERICS.md` §11 scores this id by
  backward error, so `q` is informational.

## SE_N(3)

Ids `sen3_<f>_n{1,2,3}`, one file per N. A tangent τ is the flat array `[φ; ρ₁; …; ρ_N]` (rotation
first, `docs/NUMERICS.md` §1); X is `q` (`[w, x, y, z]`) and `x = [x₁; …; x_N]`; both are hex floats.
A dense matrix is column-major in the tangent's order: column j is the image of the basis vector
e_j, m = 3 + 3N, `"shape":[m,m]`.

| Id | In | Out | Reference, from the definition |
|---|---|---|---|
| `sen3_exp` | `tau` | `q`, `x` | `mp.expm` of the (3+N)-square hat matrix [[φ^, ρ₁ … ρ_N], [0, 0]]: `x` is its translation block, `q` the quaternion series of `so3_exp` |
| `sen3_log` | `q`, `x` | `tau` | φ as `so3_log`; ρ_i solves J_l(φ) ρ_i = x_i, J_l the series of the same exponential's translation block. Each `q` is followed by −q, except w = +0 |
| `sen3_ad` | `q`, `x` | `Ad` | the images of the basis under σ ↦ (X σ^ X⁻¹)^∨, X = [[R(q/‖q‖), x], [0, I]], X⁻¹ from `mp.inverse` of R |
| `sen3_jr`, `sen3_jl` | `tau` | `J` | Σ (∓ad_τ)ⁿ/(n+1)! as dense m×m matrices; ad_τ is the images of the basis under σ ↦ [τ^, σ^]^∨ |
| `sen3_jr_inv`, `sen3_jl_inv` | `tau` | `J` | `mp.inverse` of the above, the blocks off the diagonal and the first block column exactly 0 |

- **Strata.** The θ strata at unit translation scale, minus `theta:dense` (a grid for the
  coefficients' switch points), 27; then 25 cells `rho:1e<e>/theta=<θ>`, e ∈ {−6, −3, 0, 3, 4} and
  θ ∈ {1e-8, 1e-4, 1e-1, 1, pi-1e-6} (fl(π − 10⁻⁶) as in `theta:pi-1e-k`), each with `rho_exp` e:
  52 strata. `sen3_log` and `sen3_ad`, which take a quaternion, also see `q:w0` (q = (+0, u): an
  angle of exactly π, the sign rule of `NUMERICS.md` §3.2, no −q) and `q:nonunit` (|q|² − 1 =
  ±2⁻⁴⁵) of `so3_log`, at unit translation scale: 54. `phi_of` and `quat_of` of `so3` turn
  (θ, axis) into φ and q; each ρ_i or x_i is a random direction at the stratum's scale, from the
  stream `rho<i>` of that stratum, so a file of N = 2 extends the records of N = 1. A random θ
  takes the first records of the stratum's stream; `theta:exact0` has φ = 0 and q = (1, 0, 0, 0),
  where a Log with `0/0` in it shows.
- **Samples.** `SEN3_SAMPLES` = 6 per stratum: an SE_N(3) file has 312 records; `ad` 324 and `log`
  642 (each `q` and its negative) with the two `q:*` strata. `PHASE1.md` §4.4's 64 would be 190
  MB, since each of `ad` and the four Jacobians is 5 dense matrices per record. 6 is a round
  budget trade-off, not the largest count that fits: the family is 18 of the 50 MB and about 12 of
  the 17 CPU-minutes of generation, the rest is left for the later phases' ids, and each further sample per stratum costs ~3 MB and ~2 CPU-minutes (one constant and a
  regeneration).
- **Dense inverses.** LU with pivoting leaves rounding residue in the blocks a dual matrix has no
  entry in, and the 150-digit recheck compares entrywise. `sen3.inverse` sets exactly those blocks
  to 0 (`sen3.structural_zero`, the pattern `check_sen3._dual` asserts) and keeps every other
  entry as the LU has it, down to `theta:subnormal`'s 1e-311: the 120- and 150-digit inverses
  agree on all of them.
- **Cross-checks** (`check_sen3.py`), all to 1e-100 of max(1, the largest entry): `exp`: R(q) is
  `mp.expm` of the 3×3 hat, x_i = J_l(φ) ρ_i, and Log(Exp τ) = τ. `log`: `so3_log`'s checks on φ,
  and `mp.expm` of the whole hat matrix returns X. `ad`: X Exp(t) X⁻¹ = Exp(Ad t) for one t, both
  sides `mp.expm`, X inverted whole. `jr`, `jl`: J_l = Ad_Exp(τ) J_r with Exp by `mp.expm`, against
  the other series; the inverses: J J⁻¹ = I against the series. Every dense output asserts the
  dual-matrix structure (`NUMERICS.md` §2.2): the blocks off the diagonal and the first block
  column are exactly 0, the diagonal blocks are one block. Not one check evaluates a block form of
  `NUMERICS.md` §5; `tests/test_sen3.py` compares the corpus functions to them.
- **`mp.logm`** of the (3+N)-square matrix is right to 1e-117 up to θ ≈ 2 and wrong near π (an
  imaginary part of order 1e4, 7e3 to 1.5e4 across the records of θ = π − 10⁻⁶, ρ = 10⁴), so it is
  a test cross-check only, as for `so3_log`.

## SO(2) and SE(2)

Ids `so2_exp`, `so2_log`, `se2_<f>`. A rotation is the unit complex `z = [c, s]`, an SO(2) tangent
the angle `theta`, an SE(2) tangent `tau = [theta; rho_x; rho_y]` (rotation first,
`docs/NUMERICS.md` §1), an element `X = (z, t)`. A dense matrix is 3×3, column-major in the
tangent's order, `"shape":[3,3]`.

| Id | In | Out | Reference, from the definition |
|---|---|---|---|
| `so2_exp` | `theta` | `z` | `mp.expm` of [[0, -θ], [θ, 0]] |
| `so2_log` | `z` | `theta` | the θ with \|θ\| ≤ π and Exp(θ) = z/\|z\|, by Newton's method on the complex series: θ += Im(conj(Exp θ) ẑ), so the error cubes |
| `se2_exp` | `tau` | `z`, `t` | `mp.expm` of the hat matrix [[0, -θ, ρ_x], [θ, 0, ρ_y], [0, 0, 0]] |
| `se2_log` | `z`, `t` | `tau` | θ as `so2_log`; ρ solves V(θ) ρ = t, V = Σ (θJ)ⁿ/(n+1)! the translation block of the same exponential (J the generator of rotations) |
| `se2_ad` | `z`, `t` | `Ad` | the images of the basis under σ ↦ (X σ^ X⁻¹)^∨, X = [[R(ẑ), t], [0, 1]], X⁻¹ from `mp.inverse` of R |
| `se2_jr`, `se2_jl` | `tau` | `J` | Σ (∓ad_τ)ⁿ/(n+1)!; ad_τ is the images of the basis under σ ↦ [τ^, σ^]^∨ |
| `se2_jr_inv`, `se2_jl_inv` | `tau` | `J` | `mp.inverse` of the above |

- **Strata.** `so2_*` see `SCALAR_THETA_STRATA`, each θ followed by its negative (0 has none): 3419
  records a file. `so2_log` takes z = (cos θ, sin θ) rounded componentwise, so |z| ≠ 1 by ~1e-16
  and the reference is that of z/|z|. `se2_*` see SE_N(3)'s 52 strata (`SE2_STRATA`, so
  `theta:dense` stays out for the same reason), `SE2_SAMPLES` = 8 records each, 416 a file: record
  i takes the i-th θ of the stratum's stream (a fixed θ repeats) with the sign (-1)ⁱ (`theta:exact0`
  stays +0), and a random direction on S¹ at the stratum's scale from the stream `rho` (`t` for
  `log` and `ad`). 8 is a budget trade-off (the nine ids are 2.3 MB), not the largest count that
  fits.
- **z = (-1, ±0), θ = π, is not sampled.** `NUMERICS.md` §6's `atan2(s, c)` gives π at s = +0 and
  -π at s = -0, and (-π, π] gives π at both, so only (-1, -0) is undecided. `so2.angle` refuses
  both: an `mpf` has no signed zero to tell them apart. No stratum reaches either (`theta:pi-1e-k`
  rounds cos θ to -1 but sin θ stays ≥ 1e-12); sampling (-1, +0) is listed as missing in
  `docs/PHASE1.md` §0.0.
- **No non-unit z.** No `q:w0` or `q:nonunit` analogue: the spec names none, and a stratum family
  is a record first. `so2_log`, `se2_log` and `se2_ad` renormalise z, but z is the rounded
  (cos θ, sin θ): `| |z|² - 1 | ≤ 1.4e-16` over their inputs, so a subject that skips the
  renormalisation errs by ~1 u and can pass.
- **Newton stops** at 1e-110 of |θ|, not of |Im z| as `so3.log_newton` does of |vec q|: near π the
  series of exp has an absolute error of 1e-120, which sin θ ≈ 1e-12 cannot bear.
- **No closed form** of `NUMERICS.md` §6 (V, α, β, `atan2`) is evaluated; it has none for the
  Jacobians yet. `tests/test_se2.py` compares Exp, Log and V⁻¹ to the closed forms and to
  `mp.atan2`, and J_r, J_l to what they mean, the derivative of `mp.expm`:
  Exp(τ + ε e_j)' = Exp(τ) (J_r e_j)^ = (J_l e_j)^ Exp(τ) (`NUMERICS.md` §1).
- **Structure**, asserted on every Jacobian and `Ad` (`check_se2._structure`): the first row is
  exactly (1, 0, 0) and the 2×2 block has the form p I + q J. `mp.inverse` leaves exactly 0 in the
  first row of all 832 series inverses in the corpus, so the SE_N(3) zeroing rule is not needed.
- **Cross-checks** (`check_se2.py`) to 1e-100 **of the size of the terms an entry sums**, not
  absolutely: at `theta:subnormal` angle-dependent entries are 1e-311 beside entries of 1, and an
  absolute 1e-100 would pass any of them. `so2_exp`: the complex series Σ (iθ)ⁿ/n! against
  `mp.expm`, |z| = 1. `so2_log`: |θ| ≤ π, `mp.expm` returns z/|z|. `se2_exp`: z as `so2_exp`,
  t = V ρ with V the series, Log(Exp τ) = τ. `se2_log`: `so2_log` on θ, `mp.expm` of the whole hat
  matrix returns t. `se2_ad`: the block is R(ẑ), X Exp(u) X⁻¹ = Exp(Ad u) for one u. `jr`, `jl`:
  J_l = Ad_Exp(τ) J_r, Exp by `mp.expm`, against the other series; the inverses: J J⁻¹ = I against
  the series.

## Coefficients

`coeff_k`, `coeff_a`…`coeff_e`, `coeff_cos_half` take θ; `coeff_r` takes `(n, w)`. Each is the
**raw** definition of `docs/NUMERICS.md` §4, never a rewrite (`a = (1 − cos θ)/θ²`,
`b = (θ − sin θ)/θ³`, `c = 1/θ² − (1 + cos θ)/(2θ sin θ)`, `d = (θ² + 2 cos θ − 2)/(2θ⁴)`,
`e = (2θ − 3 sin θ + θ cos θ)/(2θ⁵)`, `r = 2 atan2(n, w)/n`), extended analytically at its
removable singularity. `cos_half = cos(θ/2)`, the other half of `Exp`'s quaternion, is in no table
of `NUMERICS.md`; §3.1 has the series arm generate it alongside `k`, and `docs/decisions/0016`
item 1 names the id. It is the one of the three ids `docs/decisions/0015` (draft) P1.5 recommends
that is added (α and β are not). `out.value` is the coefficient at the exact input (binary64; binary32
in an `@f32` stratum). `out.d_branch` is its derivative with respect to the **branch variable** by
`mp.diff` at 120 digits, taken at the exact real θ² (not fl(θ·θ); the subject's own rounding of θ²
is part of what is measured), or n² at fixed w for `r`. The threshold sweep (`docs/PHASE1.md` §6)
compares `Dual<S, 1>` seeded on the branch variable against it.

Every coefficient id is 3420 records (`coeff_r` 3426): the binary64 strata, then their `@f32`
twins. `coeff_r` sees the θ strata as the unit quaternion of that angle,
`(n, w) = (sin θ/2, cos θ/2)` rounded to binary64, so n²/w² = tan²(θ/2) spans the decades θ² does (2.6e-25 … 4e24, and
2.5e-621 … 2.4e-619 in `theta:subnormal`): w > 0 throughout, w → 0 in `theta:pi-1e-k` (down to
5e-13; n is exactly 1 from k = 8, so those strata differ only in w), n = 0 in `theta:exact0`,
subnormal n in `theta:subnormal`. One stratum is its own: `q:w0`, `(n, +0)` for n = 1, 1e-3, 1e3
(PHASE1 §4.4's stratum name; θ = π exactly is no binary64 θ), where r = π/n and
d_branch = −π/(2n³): `Log` is scale-invariant (`NUMERICS.md` §3.2), so n ≠ 1 is a domain point, and
the norms exercise the 1/n and n⁻³ scalings. w < 0 is not sampled: `Log` flips it away, and
`r_inputs` and `r` refuse it.

## `coeff_series`

One record per coefficient (`k a b c d e r cos_half`, in that order), the committed Taylor series
the sweep and `helicoid::coeffs` are generated from (`docs/decisions/0004`):

```json
{"branch":"theta^2","coeff":"a","id":1,"prefactor":"1","series":["1/2","-1/24","1/720", …]}
```

Not a function-id file (`docs/PHASE1.md` §4.3): no `stratum`, `in` or `out`, and no hex-float
inputs; its manifest entry is `kind: "series"`. `value = prefactor · Σ series[j] · xʲ` in the branch
variable `x`; for `r`, `x = n²/w²`, prefactor `2/w`, and `d_branch = (2/w³) · Σ j · series[j] ·
xʲ⁻¹`. Sixteen terms, `"num/den"` strings, always with a denominator. Each series is mpmath's Taylor
expansion of the definition at 120 digits (`mp.taylor`), rationalized
(`Fraction.limit_denominator(10⁵⁰)`), reconstructed to 110 digits, and must equal an independent
exact derivation by power-series arithmetic on the Maclaurin series of sin, cos and atan
(`series.exact`), term by term; the leading four terms of `NUMERICS.md` §4 are asserted, never used
to produce data. `cos_half` (branch θ², prefactor 1) is in no table of `NUMERICS.md`: it is held to
the exact derivation alone. `verified` counts the series so verified; this file has no 150-digit
recheck.

## Cross-checks

Every `coeff_*` record is compared, at generation, to 100 digits, with (`check.py`):

- the committed rational series where the branch variable is ≤ `check.SMALL`, derived from
  `series.TERMS` (1e-8 for 16 terms) so that truncation stays 20 digits under the 100;
- a second formulation at 240 digits everywhere: `series.exact` summed (`k`, `a`…`e`, `cos_half`); a
  closed form and its calculus derivative (`r`, n > 0). An `@f32` record meets these at its rounded
  input, like any other.

Neither shares code with `mp.diff` on the definition. The 150-digit recheck still applies.

Every `so3_*` record is checked to 100 digits by a property that fixes its output without the
algorithm that produced it (`check.py`):

- `exp`, `log`: R(q) is `mp.expm` of the hat matrix, entries to 10⁻¹⁰⁰ and the skew part (which
  carries vec q) to 10⁻¹⁰⁰ of the off-diagonal size; R(q) is orthogonal with det +1; for `log`
  also the series residual (w > 0 side, to 10⁻¹⁰⁰ of |vec q|) and |φ| ≤ π.
- `act`: R(q/‖q‖) p is the sandwich q p q*. `from_matrix`: A = R(q) H with H symmetric positive
  definite, which pins the polar factor down.
- `jr`, `jl`: J_l = R J_r (`docs/NUMERICS.md` §1) against the other series; the inverses: J J⁻¹ = I
  against the series. The diagonal is held to 10⁻¹⁰⁰, the off-diagonal part to 10⁻¹⁰⁰ of the
  series' own size: at `theta:subnormal` J = I + 10⁻³¹⁰, and an absolute 10⁻¹⁰⁰ would pass any
  off-diagonal entry.
- `act` at `theta:subnormal` has no tiny part to check: Rp = p + O(10⁻³¹⁰) is p at the working
  precision (all 64 records), so 10⁻¹⁰⁰ is all its value carries.

## Decisions

- **Guard digits.** A definition that cancels is evaluated by `coeff.stable`: twice, at growing
  guard precision, until the two agree to the working precision (a result of exactly 0 counts as
  disagreement). 120 digits is the precision of the output (`PHASE1.md` §2): it absorbs the
  cancellation of θ ≥ 1e-12, while at `theta:subnormal` the definitions lose up to 4 · 310 digits
  and all five, `a`…`e`, evaluate to 0 at 120. `r` takes the same guard: the complex `atan` at the
  negative points of the `mp.diff` stencil is inexact at the working precision for small n.
  `mp.diff` and `mp.taylor` call the same evaluator at their own precision.
- **Analytic extension.** At branch value 0 the evaluator returns the definition at 2⁻²ᵖ (p the
  working precision in bits): the error is below the working precision. A stencil straddling 0 reaches
  negative x, where θ = √x is imaginary and the coefficients stay real; `r` continues `atan2(n, w)`
  as `atan(n/w)` there (w > 0).
- **Streams.** Every stratum draws from its own splitmix64 stream seeded by SHA-256 of
  `seed/stratum/purpose`; θ and axes are separate streams. Adding a function or stratum never
  reorders another's inputs. The seed is `"helicoid"` as a big-endian integer.
- **Rounding to binary64** is integer true division (correctly rounded, subnormals included), not
  `float(mpf)`.
- **Recheck.** At least 1% of every stratum, and at least one record of each, evenly spaced, is
  recomputed at 150 digits and must agree to 40 relative digits, or generation fails; the
  manifest's `rechecked` counts the records `recheck_stratum` checked.
- **Precision.** `Stratum.thetas` and `axes` draw at `precision.DPS` whatever `mp.dps` the caller
  set, so a stratum's inputs never depend on who drew it.
- **Identity.** `generator` in the manifest is a SHA-256 over the sorted (path, bytes) of
  `pyproject.toml`, `uv.lock`, `.python-version` and `gen/**/*.py`, not a git revision, which would
  change with every commit and cannot appear inside the commit it names. Any edit to those files
  changes the manifest, so the corpus is regenerated in the same PR.
- **Tests.** The unit tests regenerate a subset (`exact0`, `1e-6`, `q:w0` and their `@f32` twins)
  serially and with two workers and compare it with the committed records: the worker count cannot change a
  byte. The whole corpus is regenerated once, and compared to the committed one byte for byte,
  manifest included, by `just corpus-check`; a second full regeneration inside the tests would
  only double its 18 CPU-minutes.
- **Parallel.** One task per (id, stratum), joined in catalogue order, so the bytes are those of a
  serial run whatever `--jobs` (`build_all`; default every core). `just corpus` is 18 CPU-minutes
  (`so3_log` and the SE_N(3) Jacobians the most, N = 3 at 80 CPU-seconds a file); the wall time
  is that over the free cores, and `just corpus-check` adds the generator's unit tests.

## Adding a function id

1. Write `evaluate(inputs) -> outputs` in a math module: exact binary64 inputs in, mpf (or lists of
   mpf) out, computed from the definition at the current `mp.dps`, with no dependence on the closed
   forms under test.
2. Add one `FunctionSpec` to `FUNCTIONS` in `gen/registry.py` with its strata and input builder;
   a new stratum family goes in `gen/strata.py` and needs a decision record first (`q:w0` and
   `q:nonunit` are `PHASE1.md` §4.4's).
3. Add a test against an independent computation and, where a second formulation exists, a
   `check` for `FunctionSpec`; run `just corpus`, and commit the result.
