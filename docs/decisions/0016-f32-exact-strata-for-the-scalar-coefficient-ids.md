# 0016: f32-exact strata for the scalar coefficient ids

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

`PHASE1.md` §6 sweeps every coefficient of `NUMERICS.md` §4 "for `f64` and `f32`", and `0004`
generates a `Switch` per coefficient **per precision**. The corpus has one kind of input: exact
binary64 values. An `f32` kernel cannot be evaluated on them without rounding the input to `f32`,
and rounding changes the true answer: the value moves by `f'·δ`, the **derivative** by `f''·δ`,
which is first order in `u₃₂` and as large as the error being measured. The sweep's objective
(`max(value error, derivative error)`) and the harness's `f32` scoring are therefore not defined by
the current corpus (`docs/maths/coefficients.md`, `docs/decisions/0015` gap on `f32`).

## Decision

1. For every scalar-θ coefficient id (`coeff_k`, `coeff_a` … `coeff_e`, `coeff_cos_half`,
   `coeff_r`) the generator emits, beside each existing stratum `S`, a stratum `S@f32`:
   its inputs are the corresponding binary64 inputs **rounded to nearest-even `f32`** (hence exactly
   representable in both precisions; `coeff_r`: `n` and `w` each) and the references (`value`,
   `d_branch`) are computed from their definitions **at those rounded inputs**, at 120 digits, like
   every other record. `coeff_cos_half` is cos(θ/2), the series arm of `Exp` that `NUMERICS.md`
   §3.1 generates alongside `k`, a corpus id with a `coeff_series` row (`PHASE1.md` §4.3). No
   ready record specifies the sweep of `coeff_r` (its branch variable, mask and `w ≤ 0` domain);
   `cargo xtask thresholds` sweeps it on a provisional reading, pending the maintainer (0014
   (draft) question 29), which this record does not decide.
2. The harness scores `Precision::F64` on strata without the suffix and `Precision::F32` on the
   `@f32` strata only; `u = 2^-24` for `f32`. An `f32` subject receives the inputs converted to
   `f32` exactly (a lossless cast by construction). No other id has an `@f32` stratum until a
   record asks for one.
3. `cargo xtask thresholds` sweeps each precision over its own strata with the same objective, grid
   and tie-breaks (`PHASE1.md` §6; `coeff_r`'s grid is over its own branch variable, item 1).
   `generated.rs` holds one `Switch` per coefficient per precision; the `f32` series literals are
   correctly rounded **from the exact rationals**, never from the `f64` literals.
4. The corpus size budget of `PHASE1.md` §4.4 applies to the sum; the `@f32` strata of the scalar
   ids add about the size of the existing scalar strata.

## Rationale

The only way to score an `f32` kernel exactly, value and derivative, is to give it inputs it can
represent and a reference computed at exactly those inputs. Correcting a binary64 reference to
first order fixes the value but not the derivative, and rounding the inputs silently makes the
sweep measure the input's conditioning. Not shipping `f32` until then would leave `S: Real`
generic code with a precision `Real` promises and the kernel cannot serve.

## Consequences

- The generator's `MANIFEST.json` gains files or records; `just corpus-check` stays byte-identical.
- `Precision::F32` becomes available in the harness for scalar ids; vector ids (`so3_*`, `sen3_*`)
  stay `f64`-only until a record extends `@f32` to them.
- Coefficient thresholds for `f32` move only when the corpus, the kernel or the objective moves
  (`0004` item 5 and `API.md` §5 apply).

## Implementation plan

1. Generator: the `@f32` strata for the scalar ids, its README table, byte-identical regeneration —
   verified by `just corpus-check` and an independent recompute of the `@f32` records.
2. Harness: `Precision::F32` scoring on `@f32` strata — verified by the rounded-reference sanity
   subject scoring `≤ 1` in units of `2^-24`.
3. Sweep and `generated.rs`: per-precision candidates and constants — verified by
   `just thresholds-check`.

## Open questions

None.
