# 0087: `Real` rounds a runtime `f64`

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

locus-sensors turns a runtime `f64` into its scalar `S: Real`: anchor-relative seconds, formed
exactly from integer nanoseconds. `Real`'s one constructor from `f64` is `lit`, whose contract is
**exactly representable constants**: the `f32` impl `debug_assert!`s it, and a release build
rounds silently. A runtime value at `f32` or `Dual<f32, N>` therefore has no supported path, so
locus converts time only when `S::PRECISION == F64` and says so under `# Domain`.

A consumer could write the conversion itself, `match S::PRECISION { F64 => S::lit(x), F32 =>
S::lit(f64::from(x as f32)) }`. It is correct only for the two precisions that exist today, it
leans on `lit` admitting a value that is a constant only by construction, and every consumer would
carry its own copy. A lane `Real` (`0013`) could not be served by it at all.

**`Real` is not sealed.** `pub trait Real` in `crates/helicoid-linalg/src/real.rs` has no private
supertrait, and `API.md` §2 lists its impls without closing them. `0068` added `exp` and `exp_m1` as
required methods and recorded the break for an out-of-tree impl, which `0.0.x` permits. In the
workspace, `Real` is implemented by `f64`, `f32` and `Dual`, the test lanes (`L2`, `Lane`,
`Counted`) and `xtask`'s `Host` and seeded `Lane`.

## Decision

1. **`Real::round_from_f64(x: f64) -> Self`**: `x` correctly rounded to nearest at `S`'s
   precision, ties to even. The name is the request's: `from_*` names a constructor in this
   workspace, and `round_` says it is not `lit`.
   - `f64`: the identity, bit for bit, NaN included.
   - `f32`: `x as f32`. The Rust Reference (*Numeric cast*) defines it as the nearest `f32`,
     roundTiesToEven, overflowing to the infinity of `x`'s sign. That is IEEE 754's
     binary64-to-binary32 conversion, one instruction on x86_64 (`cvtsd2ss`), aarch64 (`fcvt`) and
     wasm32 (`f32.demote_f64`), and not a `libm` call. So D16 holds, NaN's payload apart as D16
     already says.
   - `Dual<S, N>`: `Dual::constant(S::round_from_f64(x))`. A value read at run time carries no
     derivative.
2. **Provided, not required.** The trait's body is `match Self::PRECISION { F64 => Self::lit(x),
   F32 => Self::lit(f64::from(x as f32)) }`. It is correct for every impl that honours `lit`'s
   contract, because the argument is exactly representable at the precision. So an out-of-tree
   `impl Real` keeps compiling. `f64`, `f32` and `Dual` override it with the cast alone, without
   `lit`'s debug check, whose bit comparison a NaN widened and narrowed again need not pass on
   wasm32.
3. **Domain** (`NUMERICS.md` §12): every `f64`, and nothing is asserted (`API.md` R6 asks an
   assert only of a restricted domain). At `f32`:
   - everything from `(2 − 2⁻²⁴)·2¹²⁷` up is `±∞`. That point is a tie between `MAX`, whose
     significand is odd, and `2¹²⁸`;
   - the binary32 subnormals are reached below `2⁻¹²⁶`;
   - everything from `2⁻¹⁵⁰` down is `±0`, the sign kept. That point is a tie with `2⁻¹⁴⁹`;
   - NaN stays NaN.
4. **Corpus id `real_round_from_f64`** (`PHASE1.md` §4.3, §4.4), scored by `Forward` on `value`
   and on the `Dual` lane `d` (`REAL_UNARY`). 64 records per stratum, signs alternating:
   - binary64: `x:generic` (`|x|` log-uniform in `[1e-300, 1e300]`) and `x:subnormal`
     (`[1e-320, 1e-308]`);
   - `@f32`: `x:generic@f32` (`[1e-37, 3e38]`);
   - `x:tie@f32`: the midpoint `(2m + 1)·2^(e−24)` of two adjacent normal binary32s, below the
     edge;
   - `x:near-tie@f32`: a midpoint's binary64 neighbour, above and below alternately;
   - `x:overflow@f32`: `MAX + k·2⁷⁵`, from `MAX`'s lower neighbour to the binary64 below the edge.
     The edge itself rounds to `∞`, which the harness counts as non-finite, so the unit test holds
     it;
   - `x:underflow@f32`: log-uniform in `[2⁻¹⁵⁵, 2⁻¹²⁵]` interleaved with ties `k·2⁻¹⁵⁰`, `k` odd,
     with `2⁻¹⁵⁰ → 0` and `3·2⁻¹⁵⁰ → 2⁻¹⁴⁸` first.
5. **This id's `@f32` inputs are not binary32**, the one exception to `0016` item 1. It converts a
   binary64, so its input is a binary64 at every precision. Its reference is the input correctly
   rounded at the stratum's precision (`to_f32`, `gen/precision.py`), rather than the function at a
   rounded input. The generator's `FunctionSpec.converts` and the harness's `corpus::converts`
   name it, and they skip the binary32 input check for it alone. Outputs stay binary32, checked as
   before.
6. **The bar is the reference's own rounding.** A correct conversion reads only the 30-digit
   decimal's error: about `1e-14 u` at binary64 and `1e-23 u` at binary32. A binary32 one ulp off
   reads at least 1. A tie broken the wrong way reads at least 0.5, half a spacing over
   `|y|·2⁻²⁴`. The rows are committed in `xtask/src/conformance/committed_linalg.rs`, so
   `just test` reproduces them bit for bit on x86_64 and aarch64 before the envelope's first
   bless.

## Rationale

- **A trait method, not a free function.** The precision's rounding belongs to the scalar. A lane
  `Real`, or a future precision, overrides it once, where a consumer's `match` would silently fall
  through.
- **Provided** costs nothing in the shipped impls, which override it. It keeps the change additive
  for an out-of-tree impl. A required method would add one line to each of six in-tree impls whose
  body is the default anyway.
- **`as` over a hand-rolled rounding**: it is already IEEE's conversion on every D16 target. The
  corpus's reference rounds by integer arithmetic (`to_f32`), and the generator cross-checks it
  against the platform's C `(float)` conversion through `struct` (`check_real`), so neither side
  relies on the other.
- **Rejected: relaxing `lit`'s contract** to admit runtime values. `lit`'s debug check is what
  catches a constant that should have been generated (`coeffs::generated`), and that check would be
  lost.

## Measured

`cargo xtask conformance --fn real_round_from_f64`, both precisions, no non-finite output:

| stratum | `max_u` |
|---|---|
| `x:generic` | 3.04e-14 |
| `x:subnormal` | 1.68e-15 |
| `x:generic@f32` | 6.10e-23 |
| `x:tie@f32` | 5.07e-23 |
| `x:near-tie@f32` | 5.04e-23 |
| `x:overflow@f32` | 2.38e-23 |
| `x:underflow@f32` | 5.06e-23 |

Every other corpus file regenerates byte-identical. The manifest gains the file, 448 records, and
the generator's hash.

## Consequences

- `Real` gains a provided method. No public item is removed or changed, and an out-of-tree impl
  still compiles. A type that already has an inherent or other-trait method of the same name may
  see an ambiguity, which is the Rust API guidelines' "possibly breaking" minor change.
- `lit` keeps its contract and its check.
- locus-sensors can convert at every precision and drop its `F64`-only `# Domain` note.
- `real_round_from_f64` is required by coverage. Its rows are pinned in `committed_linalg`.
- The generator and the harness carry one named exception to `0016` item 1.

## Implementation plan

1. One PR, verified by `just test`, `just corpus-check` and `cargo xtask conformance --fn
   real_round_from_f64` at both precisions:
   - the method, its overrides and their tests:
     - `round_from_f64_rounds_to_nearest_even_at_binary32` (ties, both edges, signed zero, NaN);
     - `round_from_f64_is_the_identity_at_f64_and_a_constant_for_dual`;
     - `the_provided_round_from_f64_rounds_at_the_lane_s_precision`;
   - the corpus id, with its generator tests `test_the_round_strata_hold_what_they_are_named_for`
     and `test_round_cross_check_refuses_another_rounding`;
   - the conversion exemption;
   - the docs: `PHASE1.md`, `PHASE2.md`, `NUMERICS.md` §12, `API.md` §2.

## Open questions

None. Answered by the owner:

1. **Provided, not required.** Out-of-tree implementations keep compiling (Decision 2). A future
   lane `Real` (`0013`) overrides it as `f32` and `f64` do.
2. **The name is the request's,** `round_from_f64`. `value_f64` is its inverse at binary64 and
   stays tests-only (`API.md` R4), so the pair is asymmetric on purpose.
