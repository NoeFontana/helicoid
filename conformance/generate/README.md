# Reference generator

Computes the conformance corpus from **definitions** with mpmath at 120 digits, never from the
closed forms in `docs/NUMERICS.md` (`docs/PHASE1.md` §2, §4; `docs/decisions/0006`). Inputs are exact
binary64 values; the corpus in `../corpus/` is committed and regenerates byte-identically.

```
just corpus         # regenerate ../corpus, deleting a stale *.jsonl   (python -m gen all --out DIR)
just corpus-check   # regenerate into a temp dir, compare every file byte for byte, run the tests
just corpus-test    # unit tests only (stdlib unittest)
uv run --frozen python -m gen list
```

## Files

| Path | Role |
|---|---|
| `pyproject.toml`, `uv.lock`, `.python-version` | CPython 3.12, `mpmath` pinned exactly; no numpy, no gmpy2 (`setup` refuses a non-Python backend) |
| `gen/precision.py` | the only place `mp.dps` is set (120; rechecks at 150); `to_f64` |
| `gen/rng.py` | splitmix64, per-stratum streams, log-uniform, uniform on S² |
| `gen/strata.py` | the strata; `gen/registry.py` the function ids; `gen/coeff.py` their definitions |
| `gen/series.py`, `gen/check.py` | the exact Taylor series (`coeff_series`); the per-record cross-checks |
| `gen/fmt.py`, `gen/corpus.py`, `gen/manifest.py` | text formats, assembly and recheck, `MANIFEST.json` |

## Record schema

One JSONL file per function id, one record per line, compact JSON with sorted keys, LF newlines:

```json
{"id":832,"in":{"theta":"0x0.0p+0"},"out":{"d_branch":"-2.08333333333333333333333333333e-2","value":"5.00000000000000000000000000000e-1"},"stratum":"theta:exact0"}
```

- `id`: 0-based line number in the file. Strata are emitted in the order of their catalogue, so a
  stratum appended to the end never renumbers the rest.
- `in`: exact binary64 as `float.hex()`. `out`: decimal strings, 30 significant digits, `d.ddd…e<n>`
  (no `+`, no zero padding), rounded half-to-even from the 120-digit value in integer arithmetic.
- Vectors are arrays. Matrices (column-major flat array plus a sibling `shape`) are not emitted yet;
  the first matrix id fixes the key name.
- `MANIFEST.json`: generator identity, Python (major.minor) and mpmath versions, `dps`, seed, recheck
  parameters, and per file the SHA-256, record count, `kind` and one check count: `rechecked`
  (records recomputed at 150 digits) for `kind: "corpus"`, the files of the schema above;
  `verified` (series equal to exact algebra) for `kind: "series"`, i.e. `coeff_series`, which has
  its own record schema below. A reader of function-id files skips every other `kind`.

## Strata for scalar-θ ids (`SCALAR_THETA_STRATA`)

| Stratum | θ values (binary64) | Records |
|---|---|---|
| `theta:1e-k`, k = 12…1 | log-uniform in the decade [10⁻ᵏ, 10⁻ᵏ⁺¹) | 64 each |
| `theta:1e0` | log-uniform in [1, π − 0.1) | 64 |
| `theta:exact0` | 0 | 1 |
| `theta:subnormal` | log-uniform in [10⁻³¹⁰, 10⁻³⁰⁹), below the smallest normal | 64 |
| `theta:pi-1e-k`, k = 1…12 | the single value fl(π − 10⁻ᵏ) | 1 each |
| `theta:dense` | 10^(j/200 − 4), j = 0…800, rounded: a grid, no randomness | 801 |

For vector-valued ids every stratum except `theta:exact0` and `theta:dense` (θ only) carries axes
(`Stratum.axes`: 64 uniform on S², z uniform then φ uniform, from the stratum's `axis` stream). A
scalar id sees each fixed-θ stratum once.

## Coefficients

`coeff_k`, `coeff_a`…`coeff_e` take θ; `coeff_r` takes `(n, w)`. Each is the **raw** definition of
`docs/NUMERICS.md` §4, never a rewrite (`a = (1 − cos θ)/θ²`, `b = (θ − sin θ)/θ³`,
`c = 1/θ² − (1 + cos θ)/(2θ sin θ)`, `d = (θ² + 2 cos θ − 2)/(2θ⁴)`,
`e = (2θ − 3 sin θ + θ cos θ)/(2θ⁵)`, `r = 2 atan2(n, w)/n`), extended analytically at its removable
singularity. `out.value` is the coefficient at the exact binary64 input. `out.d_branch` is its
derivative with respect to the **branch variable** by `mp.diff` at 120 digits, taken at the exact
real θ² (not fl(θ·θ); the subject's own rounding of θ² is part of what is measured), or n² at fixed
w for `r`. The threshold sweep (`docs/PHASE1.md` §6) compares `Dual<S, 1>` seeded on the branch
variable against it.

`coeff_r` sees the θ strata as the unit quaternion of that angle, `(n, w) = (sin θ/2, cos θ/2)`
rounded to binary64, so n²/w² = tan²(θ/2) spans the decades θ² does (2.6e-25 … 4e24, and
2.5e-621 … 2.4e-619 in `theta:subnormal`): w > 0 throughout, w → 0 in `theta:pi-1e-k` (down to
5e-13; n is exactly 1 from k = 8, so those strata differ only in w), n = 0 in `theta:exact0`,
subnormal n in `theta:subnormal`. One stratum is its own: `q:w0`, `(n, +0)` for n = 1, 1e-3, 1e3
(PHASE1 §4.4's stratum name; θ = π exactly is no binary64 θ), where r = π/n and
d_branch = −π/(2n³): `Log` is scale-invariant (`NUMERICS.md` §3.2), so n ≠ 1 is a domain point, and
the norms exercise the 1/n and n⁻³ scalings. w < 0 is not sampled: `Log` flips it away, and
`r_inputs` and `r` refuse it.

## `coeff_series`

One record per coefficient (`k a b c d e r`, in that order), the committed Taylor series the sweep
and `helicoid::coeffs` are generated from (`docs/decisions/0004`):

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
to produce data. `verified` counts the series so verified; this file has no 150-digit recheck.

## Cross-checks

Every `coeff_*` record is compared, at generation, to 100 digits, with (`check.py`):

- the committed rational series where the branch variable is ≤ `check.SMALL`, derived from
  `series.TERMS` (1e-8 for 16 terms) so that truncation stays 20 digits under the 100;
- a second formulation at 240 digits everywhere: `series.exact` summed (`k`, `a`…`e`); a closed form
  and its calculus derivative (`r`, n > 0).

Neither shares code with `mp.diff` on the definition. The 150-digit recheck still applies.

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
- **Serial.** A full run takes about ten seconds, most of it the cross-checks. Per-stratum streams make a process pool
  order-independent when a later id needs one.

## Adding a function id

1. Write `evaluate(inputs) -> outputs` in a math module: exact binary64 inputs in, mpf (or lists of
   mpf) out, computed from the definition at the current `mp.dps`, with no dependence on the closed
   forms under test.
2. Add one `FunctionSpec` to `FUNCTIONS` in `gen/registry.py` with its strata and input builder;
   a new stratum family goes in `gen/strata.py` and needs a decision record first.
3. Add a test against an independent computation and, where a second formulation exists, a
   `check` for `FunctionSpec`; run `just corpus`, and commit the result.
