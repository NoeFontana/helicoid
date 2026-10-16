# 0098: `Exp`'s quaternion near `π` is a known limit

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

After [`0095`](./0095-exp-translation-takes-the-far-diagonal.md) and
[`0096`](./0096-sim3-exp-takes-chi-without-the-cancellation.md) (drafts), the deep envelope fails
one stratum: `sim3_exp/sigma:-3/theta:pi-1e-6@f32`, 2.49 `u` against Sophus's 2.13, resolved at
`p = 4.9 · 10⁻⁴`.

**It is the quaternion.** At `σ = −3`, `x`'s floor `‖[ρ; σ]‖ ≥ 3` makes its error small, so the
record's score is `q`'s. A binary32 emulation of `SO3::exp`, `(cos(θ/2), (sin(θ/2)/θ) φ)`, reproduces
the stratum's six largest scores to every digit (2.492, 2.479, 2.425, …).

**The error is radial.** Near `π`, `sin(θ/2)` rounds to `1`. What is left is a scale error on
`k = sin(θ/2)/θ`, from the roundings of `θ²`, its root, the quotient and the product, up to
about `2 u`.
- **Other programs:** sophus_rs and tf_tree_math use the same formula, and on `f64`
  `so3_exp/theta:pi-1e-9` they read the same 2.74 `u` as `helicoid` and its host-`std` twin. GTSAM
  (via its matrix), manif and Sophus read 1.90–2.01.
- **Means:** `helicoid`'s mean is at or below GTSAM's in 10 of the 12 `theta:pi-1e-k` strata. The
  loss is the maximum of 64 records.

**Measured alternatives**, by a binary64 emulation that reproduces `helicoid`'s maxima exactly, and
a consumer-loop timing of `SO3::exp`:

| form | worst, `π − 10⁻⁸` to `π − 10⁻¹¹` | cost |
|---|---|---|
| `k φ`, as shipped | 2.07–2.74 `u` | 10–13 ns |
| `(sin(θ/2) φᵢ)/θ`, one rounding where `sin(θ/2) = 1` | 1.84–2.27; worse at `π − 10⁻³`, `10⁻⁴`, `10⁻⁷` | +0.3 ns |
| `θ²` compensated (Dekker, no fused multiply-add, D16), then that division | 1.50–1.80, below every oracle | +6 ns (+50–60%) |

## Decision

1. **`SO3::exp` keeps `(cos(θ/2), k φ)`.** The one uniform improvement costs half again of the
   library's most-called operation, for under `1 u` on the maximum in a band of width `10⁻²` at
   `π`.
2. **`conformance/baseline/exceptions.toml` excepts
   `sim3_exp/sigma:-3/theta:pi-1e-6@f32`** from domination, citing this record (`0046`). It is
   watched by no-regress like every other stratum, and the row fails as stale the day an oracle
   stops winning there.

## Rationale

- **No cheap form wins:** each one moves the maximum from one near-`π` stratum to another.
- **Where it is accepted:** the owner accepted the limit on 2026-10-16, after the measurements
  above.

## Consequences

- **The deep envelope reads no domination failure** once the row is in; this record is the reason
  it does.
- **A future `Exp` that wants this back** pays the compensated `θ²`, or finds a cheaper exact
  root, and deletes the row.

## Implementation plan

1. This record and the row. Verified by `cargo xtask lint` (the row cites a `ready` record) and by
   `cargo xtask envelope --deep`, which lists the stratum as excepted and fails nothing for it.

## Open questions

None.
