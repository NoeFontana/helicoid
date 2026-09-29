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
  parameters, and per file the SHA-256, record count and number of rechecked records.

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

## `coeff_k`

`k = sin(θ/2)/θ`, extended analytically by 1/2 at 0. `out.value` is `k` at the exact θ.
`out.d_branch` is `dk/d(θ²)`, the derivative with respect to the **branch variable**, by `mp.diff` at
120 digits, taken at the exact real θ² (not fl(θ·θ); the subject's own rounding of θ² is part of
what is measured). The threshold sweep (`docs/PHASE1.md` §6) compares `Dual<S, 1>` seeded on the
branch variable against it. Every later `coeff_*` id uses the same two output names.

## Decisions

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
- **Serial.** A full run takes under a second. Per-stratum streams make a process pool
  order-independent when a later id needs one.

## Adding a function id

1. Write `evaluate(inputs) -> outputs` in a math module: exact binary64 inputs in, mpf (or lists of
   mpf) out, computed from the definition at the current `mp.dps`, with no dependence on the closed
   forms under test.
2. Add one `FunctionSpec` to `FUNCTIONS` in `gen/registry.py` with its strata and input builder;
   a new stratum family goes in `gen/strata.py` and needs a decision record first.
3. Add a test against an independent computation, run `just corpus`, and commit the result.
