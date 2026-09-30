# 0018: `libm`'s `arch` feature is bit-identical for exactly-rounded operations

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** none; a draft authorises nothing.

## Context

0007 pins `libm` with `default-features = false`, on the reading that the `arch` feature "would
break bit identity" (D16). Every `Real::sqrt`, hence `Vector::norm`, `chol`, `Dual::sqrt` and each
quaternion normalisation, therefore runs `libm`'s bit-by-bit integer `sqrt`.

What `arch` changes in `libm` 0.2.16: it routes a call to the target's instruction only where one
exists, and only for `sqrt`/`sqrtf`, `fma`/`fmaf`, `floor`, `ceil`, `trunc`, `rint` (and their `f32`
forms). IEEE 754-2019 §5.4.1 requires `sqrt` and `fusedMultiplyAdd` to be correctly rounded, and
`floor`/`ceil`/`trunc`/`rint` are exact. A correctly-rounded function has one result per input, so
the generic and the hardware paths cannot differ. Every function that is *not* exactly rounded
(`sin`, `cos`, `atan2`, `exp`, `cbrt`, …) stays on `libm`'s software path either way.

Measurement (one machine, x86_64, release, `libm` 0.2.16, 50 M calls, protocol in the
*Implementation plan*):

| | `default-features = false` | `arch` |
|---|---|---|
| `f64::sqrt`, dependent chain (latency) | 73 ns | 37 ns |
| `f64::sqrt`, independent (throughput) | 88 ns | 9 ns |
| XOR-rotate digest of the results | `04f37f3d87af5116` | `04f37f3d87af5116` |

The digests agree, as the argument predicts. The absolute numbers are one machine's; the ratio is the
point.

## Decision

None yet. The record exists so 0007's reading can be replaced by a checked one.

Proposed: enable `libm`'s `arch` feature in the workspace dependency, and have `just determinism`
(PHASE6) compare `sqrt`, `Vector::norm` and `chol` digests across x86_64, aarch64 and wasm32 with the
feature on.

## Rationale

D16 is a claim about outputs, and outputs of exactly-rounded operations are target-independent by
definition. The alternative for a hardware `sqrt` without `arch` is `core::arch` intrinsics, which need
`unsafe` (`forbid(unsafe_code)`), or `f64::sqrt`, which is `std`-only. `arch` is the only stable,
safe, `no_std` route.

## Consequences

If accepted: `sqrt` moves from tens of nanoseconds to a few on every target with a hardware path, a
larger gain than any algorithmic change in `helicoid-linalg`. 0007's comment in the workspace
`Cargo.toml` is rewritten to state the exact-rounding argument. A new `libm` release that routes a
*non*-exactly-rounded function through `arch` would break the argument; the determinism digest is the
tripwire, and `libm` stays pinned to a minor version.

## Implementation plan

1. Bench protocol committed under `crates/helicoid-linalg/benches/`: `sqrt` latency and throughput,
   `Vector::norm`, `chol::<6>`, with and without `arch` — verified by `just bench`.
2. Enable `arch`; add the sqrt/norm/chol digests to `just determinism` — verified by the cross-target
   comparison being equal, on all three targets.
3. Rewrite the `Cargo.toml` comment and the 0007 text — verified by `just lint`.

## Open questions

1. **Does `libm` 0.2.16 `arch` route anything not exactly rounded on any of the three targets?**
   Read the dispatch table for each target, not only x86_64; the determinism digest is the check.
2. **x86 without SSE2 (`i586`)?** The dispatch table has an x87 path, which rounds in extended
   precision. Out of the target set (D17), but the record should say so.
3. **Does 0007's "nothing else" budget count a feature flag of an existing dependency?** It adds no
   crate; state that here.
