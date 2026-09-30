# 0018: `libm`'s `arch` feature is bit-identical for exactly-rounded operations

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

The workspace pins `libm` with `default-features = false`, on the reading that the `arch` feature
"would break bit identity" (D16). In `libm` 0.2.16 `arch` is a *default* feature, so that line is what
switches it off. Every `Real::sqrt`, hence `Vector::norm`, `chol`, `Dual::sqrt` and each quaternion
normalisation, therefore runs `libm`'s generic `sqrt` (Goldschmidt iterations on integers, correctly
rounded), not the target's instruction.

**What `arch` changes in `libm` 0.2.16** (`src/math/arch/mod.rs`, and the `select_implementation!`
site of each function): it replaces the body of a function with the target's instruction, for
`sqrt`/`sqrtf`, `fma`/`fmaf` and `rint`/`rintf`, where the target has one. Read per target, on stable
Rust with `default-features = false, features = ["arch"]`:

| Target | Routed to hardware | Under |
|---|---|---|
| `x86_64-*`, `i686-*` | `sqrt`, `sqrtf` (`sqrtsd`, `sqrtss`); `fma`, `fmaf` (`vfmadd*`, or FMA4) | `target_feature = "sse2"` (baseline). The `fma` choice is made at run time by `cpuid`; a CPU with neither runs the generic `fma`, the code that runs with `arch` off. No `floor`, `ceil`, `trunc` or `rint`. |
| `aarch64-*` (Linux, Apple) | `sqrt`, `sqrtf` (`fsqrt`); `fma`, `fmaf` (`fmadd`); `rint`, `rintf` (`frintn`) | `target_feature = "neon"` (baseline). `aarch64-unknown-none-softfloat` has no `neon`: nothing is routed. No `floor`, `ceil`, `trunc`. |
| `wasm32-*` | **nothing** | The branch needs `intrinsics_enabled`, i.e. the nightly-only `unstable-intrinsics` feature. On stable `arch` is a no-op here. |
| `thumbv7em-none-eabihf` | **nothing** | `arch/mod.rs` has no Arm-M branch (`VSQRT` exists, `libm` does not use it). |
| `i586` (x86 without SSE2) | nothing new | `floor`, `ceil` (x87 `frndint`) and the x87 `exp` family are selected by `use_arch_required`, so they are on *with or without* `arch`. Outside the target set (D17). |

So `arch` reaches `sqrt` and `sqrtf` on x86_64 and aarch64 and, of what `helicoid` calls
(`sqrt`, `sincos`, `atan2` and their `f32` forms, `fabs`, `copysign`, `ldexp`), nothing else.

**Who consumes the routed primitives inside `libm`** (grep of `src/math`, the code that runs
whichever body is selected): `sqrt`/`sqrtf` are called by `hypot`, `hypotf`, `pow`, `powf`, `acos`,
`acosf`, `asin`, `asinf`, `acosh`, `acoshf`, `asinh`, `asinhf`, `j0`, `j1`, `jn` and the `f32` forms;
`fma` by `cbrt` alone (through `Float::fma`, which is the public `fma`); `floor` by `rem_pio2_large`
(behind `sin`, `cos`, `tan`, `sincos` for very large arguments), `lgamma_r`, `tgamma`; `rint` by
nothing internal (`rem_pio2` rounds with the `TO_INT` addition). `sin`, `cos`, `atan2`, `exp`,
`log`, `pow`, `cbrt` and the rest are software programs over `+ - * /`, integers, and those
primitives.

**Does an algorithm, not just a primitive, differ?** Read in full: no. `arch` only chooses which
function body computes `sqrt`, `fma` or `rint`; no caller forks on it. The `target_feature`
dependent forks in `libm` 0.2.16 are all x86-without-SSE2 (`rem_pio2` and `rem_pio2f` `force_eval`,
the x87 `exp` family, `floor`, `ceil`), selected by `cfg(x86_no_sse)` or `not(sse2)` and *not* by
`arch`. `libm` does not branch on `target_feature = "fma"` anywhere: x86_64 baseline has no FMA and
still reaches a hardware `vfmadd` by `cpuid`, aarch64 always has it, and both give the fused result.

**Measured, per entry point** (scratch crate, not committed). 95 digests over the entry points
(`f64` and `f32`: `sqrt cbrt sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh exp exp2
exp10 expm1 log log2 log10 log1p erf erfc tgamma lgamma j0 j1 y0 y1 floor ceil trunc round rint
roundeven fabs atan2 pow hypot fmod remainder copysign fmin fmax fdim sincos frexp modf fma fmaf`),
each over `4 x 10^5` inputs (splitmix64: a fifth raw bit patterns, a fifth log-uniform `2^-40..2^40`,
a fifth `[-1, 1)`, a fifth halves up to 1000, a fifth exponents over the whole range; `fma` also with
`z = -fl(x y)`, which exposes the fused low part). One FNV-1a digest per function, NaN mapped to one
value. Result: **every digest is identical** across all six builds,

- `x86_64-unknown-linux-gnu`, `arch` off and on (AMD EPYC Milan, FMA present, so the hardware `fma`
  ran),
- `wasm32-wasip1`, `arch` off and on (node 24 WASI host),
- `aarch64-unknown-linux-gnu`, `arch` off and on (`qemu-user` 11.0.3; the `arch` binary has `fsqrt`,
  `fmadd`, `frintn` in its disassembly, the other has none of the three).

**NaN bits are the one thing that differs**, and only for `sqrt` and `fma` (the exactly-rounded claim
below is about numbers). Bits of the result on the same inputs:

| | generic (`arch` off, every target) | `arch`, x86_64 | `arch`, aarch64 |
|---|---|---|---|
| `sqrt(-1)`, `sqrt(-inf)`, `sqrtf(-1)` | `7ff8000000000000`, `7fc00000` | `fff8000000000000`, `ffc00000` | `7ff8000000000000`, `7fc00000` |
| `sqrt` of a quiet NaN with payload `1` | `7ff8000000000000` (payload dropped) | `7ff8000000000001` | `7ff8000000000001` |
| `fma(inf, 0, 1)` | `7ff8000000000000` | `fff8000000000000` | `7ff8000000000000` |
| `inf - inf` (plain arithmetic, no `libm`) | `fff8000000000000` (x86_64) | | `7ff8000000000000` (aarch64) |

The sign of the default NaN is the instruction set's (x86 "real indefinite" is negative, Arm's is
positive), so plain `+ - * /` already produced target-dependent NaN bits before this record, and
D16 never covered them: `error-analysis.md` EA.11(c) states the claim "NaN payloads apart", and a NaN
output is a failure of the corpus (EA.11(e), the `nonfinite` column). `arch` adds `sqrt` and `fma` to
the operations whose NaN sign and payload are the target's.

The speed-up is not asserted here: the implementing PR measures it (*Consequences*).

## Decision

1. **Enable `libm`'s `arch` feature** in the workspace dependency, keeping `default-features = false`
   so that a default feature added by a later `libm` release does not enter unreviewed:
   `libm = { version = "0.2", default-features = false, features = ["arch"] }`. No other feature
   (`unstable-intrinsics`, `force-soft-floats`, `unstable-float` stay off), no `target-cpu`, no
   `target-feature`, no `mul_add` in `helicoid`.
2. **The claim, stated exactly.** With `arch`, every `libm` function that `helicoid` calls returns,
   for every input whose result is not NaN, the same bits on x86_64, aarch64 and wasm32. The sign
   and payload of a NaN are the target's, as those of plain arithmetic already were (Context), so a
   program is target-independent on every input inside its documented domains, where no NaN
   intermediate arises. A program that reads a NaN's sign bit sees the target's: `Real::copysign`
   and `Dual::copysign` do, so `copysign(1, sqrt(-1))` is `-1` in a release build on x86_64 with
   `arch`, and `+1` on aarch64, on wasm32 and before this record. That input is outside
   `Real::sqrt`'s domain (a `debug_assert!` in debug builds) and a NaN output fails the corpus
   (`error-analysis.md` EA.11(e)); no workspace code takes a NaN's sign as a value on a valid input.
3. **A permanent test is the tripwire**, target-independent: `Real::sqrt` for `f64` and `f32`
   against an integer square root written in the test (round to nearest even on `u128`), bit for
   bit, over a deterministic stream of at least `10^6` positive finite cases per type: subnormals,
   exact squares, powers of two, both sides of every `2^k`, values one ulp either side of perfect
   squares, the largest and smallest finite values. The special values are fixed-value tests beside
   the stream: `+-0` and `inf` exactly, NaN by `is_nan`, and a negative argument (release: NaN;
   debug: the documented `debug_assert!`). A **pinned digest** of the stream's results, equal on
   every target, pins the stream and the oracle against drift and is the cross-target evidence.
   The test passes with `arch` on and off, so it asserts the property and not the flag. A `libm`
   update that changes one bit of `sqrt` fails it; a routing change that leaves the bits alone is
   invisible to it, by design (Rationale).
4. **The `libm` version is reviewed against the dispatch table above** whenever `Cargo.lock` moves
   `libm`: a new `arch` route for a function that is not exactly rounded is the case that would
   break this record, and only the review, then `just determinism` (`PHASE6.md` §3, through
   `sincos` and `atan2` on the corpus), can see it.

## Rationale

D16 is a claim about outputs, and the outputs of an exactly-rounded operation are target-independent
by definition. Each routed primitive is one:

- `sqrt`: IEEE 754-2019 §5.4.1 requires it correctly rounded; `sqrtsd`/`sqrtss`, `fsqrt` are, in the
  default rounding mode Rust assumes; and the generic body is too (the test in Decision 3 checks it
  against an integer square root, with `arch` off as well as on).
- `fma`: §5.4.1 fused multiply-add, one rounding. `vfmadd`, `fmadd` and `libm`'s generic
  `fma_round` / `fma_wide_round` all return the correctly rounded value; the `4 x 10^5` cancelling
  and raw inputs per type above agree bit for bit between hardware and software, on this CPU.
  The only caller inside `libm` is `cbrt`, and it consumes the result, not the path.
- `rint`/`rintf` (`frintn`): round to nearest even, exact; no internal caller; the generic
  `rint_round(_, Nearest)` is the same function.
- `floor`, `ceil`, `trunc`: exact, and not routed on any target of D17 on stable.

A correctly-rounded function has one result per input, so the generic and the hardware paths cannot
differ on numbers; and every function above them (`sin`, `atan2`, `cbrt`, `pow`, `hypot`, ...)
runs the same program either way, so it cannot differ either. The digests above are the check of that
argument, not a substitute for it.

Alternatives lost: `core::arch` intrinsics (`unsafe`, `forbid(unsafe_code)`); `f64::sqrt` (`std`
only, and `std` maths is platform maths, D16); `libm`'s `unstable-intrinsics` (nightly, and a
different route on wasm32 only). `arch` is the only stable, safe, `no_std` route. The `asm!` lives
in `libm`; `forbid(unsafe_code)` in `0007` binds `helicoid`'s own roots.

**Q3, the budget.** `0007` item 1 counts crates (`cargo xtask lint` compares the normal-dependency
closure by package name; `deny.toml` bans by crate). A feature of a listed dependency adds no crate
and no edge, and `Cargo.lock` does not change; `0007`'s "nothing else" is unaffected. `PROJECT.md`
§6 lists "a feature to a library crate" (D10) as a stop-and-ask smell: that is a cargo feature of the
library itself (`mint`); a feature of a listed dependency changes the dependency set, and this
record is the review D10 asks for.

## Consequences

- `sqrt` is a few nanoseconds on x86_64 and aarch64, a larger gain than any algorithmic change in
  `helicoid-linalg`. The implementing PR records the measured effect on `Vector::norm`, `chol::<6>`
  and `Dual::sqrt`.
- **No change at all on wasm32 (stable) and thumbv7em-none-eabihf**, and none on
  `aarch64-unknown-none-softfloat`: nothing is routed there (Context). The cross-target digest
  stays trivially equal, and no speed-up should be expected on those targets.
- `Real::sqrt` of a negative argument returns a NaN whose sign is the target's in release builds
  (`Real::sqrt` already documents "NaN"). The permanent test compares NaN by `is_nan`, so its
  digest is target-independent.
- The comment in the workspace `Cargo.toml` states the exact-rounding argument and cites this
  record; `docs/maths/error-analysis.md` EA.11(c) and its `Checked` line stop saying that `arch`
  would break bit identity. Records `0007` and D16 need a citation, not a new rule.
- `i586` (x86 without SSE2) is unaffected by this record: its x87 paths are on with or without
  `arch`, and the x87 `exp` family is documented by `libm` as up to 1 ulp hardware-dependent. It is
  outside D17's target set and would need its own record.
- A `libm` release that routes a function that is not exactly rounded through `arch` would break
  this record (Decision 4).

## Implementation plan

1. **PR #32: ratify.** The record, the index rows. No code.
2. **Stacked PR: enable.** In this order, one commit, verified by the recipes named:
   - workspace `libm` dependency and its comment; the sentences of `error-analysis.md` that say
     `arch` breaks bit identity; a citation, with the NaN caveat, in `0007`, `PROJECT.md` D16,
     `PHASE2.md` §0.0 and `PHASE6.md` §3; `CHANGELOG.md` — verified by `just lint` and
     `just audit`;
   - the integer-`sqrt` test with the pinned digest in `helicoid-linalg` — verified by `just test`
     on x86_64, on aarch64 (CI's `ubuntu-24.04-arm` job), and under `wasm32-wasip1` on a WASI host;
   - the before/after measurement (protocol: scratch crate over `helicoid-linalg` built with and
     without `arch`, release, one core, 7 timed repeats per row, median of three interleaved
     runs; rows `Vector::norm` (3), `chol::<6>`, `Dual::sqrt`, bare `sqrt`, `f64` and `f32`,
     independent inputs (throughput) and a dependent chain (latency)) — recorded in *Consequences*;
   - `just no-std`, `just msrv`, `just wasm` — verified by the recipes (build only).

## Open questions

None. Resolved from `libm` 0.2.16's source and the measurements above:

1. **Does `arch` route anything not exactly rounded, on any target?** No. Per target, `sqrt`,
   `fma`, `rint` at most (Context table); nothing on wasm32 stable or thumbv7em. Nothing in the
   routed set is inexact, and the only algorithm forks by `target_feature` are x86-without-SSE2 and
   independent of `arch`.
2. **x86 without SSE2 (`i586`)?** Unaffected by `arch`; out of the target set (Consequences).
3. **Does 0007's budget count a feature flag?** No, it counts crates (Rationale).
