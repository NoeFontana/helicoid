# 0089: The closed arms take θ² as given, and `J_r`'s diagonal does not cancel

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

Under [`0088`](./0088-the-domination-bar-must-tell-worse-from-equal.md)'s bar, ten of the
nineteen resolved domination failures are SO(3)'s Jacobian family near `π`:
- `so3_jr` and `so3_jl` at `theta:pi-1e-1`, `-5`, `-8` and `-12`;
- `so3_jr_inv` and `so3_jl_inv` at `pi-1e-6`;
- `so3_gamma2` at `pi-1e-3` and `pi-1e-8`.

Every oracle assembles these the way `helicoid` does, `s0 I + s1 W + s2 W²` with `W²` from the
product. The losses are in two places:

1. **The coefficients' closed arms re-square the rounded root.** `b = (θ − sin θ)/(θ̂·θ̂·θ̂)` with
   `θ̂ = fl(√z)`, `a = 2k²` with `k = sin(θ/2)/θ̂`, `d = (θ̂² − 4 sin²(θ/2))/(2θ̂⁴)`. `W²`'s entries are
   formed from `φ`, whose rounded square norm is `z`, so a coefficient over `θ̂²` carries the
   rounding of `√z` and of the square on top of `z`'s. Near `π`, `a`'s relative sensitivity to `θ`
   is `−2`. GTSAM divides every coefficient by `theta2 = ω·ω` (`SO3.cpp` at the pinned commit:
   `B = 2 sin²(θ/2)/θ²`, `C = (1 − sin θ/θ)/θ²`, `E = (1 − 2B)/(2θ²)`). Sophus forms
   `c = (1 − (θ/2) cot(θ/2))/(θ·θ)`.
2. **`J_r`'s diagonal cancels near `π`.** `1 + bW²ᵢᵢ = 1 − b(θ² − φᵢ²)`, and `bθ² → 1` at `π`, so the
   entry is `nᵢ²` out of two terms of opposite sign. The identity `1 − bθ² = sin θ/θ` gives the
   same entry as `sin θ/θ + bφᵢ²`, two terms of one sign.

## Decision

1. **The closed arms divide by the `z` they were given,** in GTSAM's arithmetic and order:
   - `a = 2 sin²(θ/2)/z`;
   - `b = (1 − sin θ/θ)/z`;
   - `d = (1 − 2a)/(2z)`, written `(½ − a)/z`: a factor of two is exact in both operands, so the
     same bits with two multiplies fewer.

   `c` becomes Sophus's `(1 − (θ/2) cot(θ/2))/(θ·θ)`. `e` is unchanged. Every group that holds a
   member computes it the same way, so each member is still its own group's to the bit. No new
   transcendental is introduced: each arm takes the root and the sines it took before.
2. **The switches are the sweep's,** regenerated against these arms by `cargo xtask thresholds`.
   No switch is typed.
3. **`SO3::jr` takes `sin θ/θ + bφᵢ²` on its diagonal where `b` is on its closed arm,** with the
   `sin θ/θ` that arm formed `b` from (`jr_coeffs_sinc`), and `1 + bW²ᵢᵢ` below it.
   - The choice is a `select` on `b`'s own mask (`HatSq::poly_far`), so there is no new switch.
   - `jr_coeffs` is `jr_coeffs_sinc`'s `(a, b)`, so every other reader of `(a, b)` takes the same
     bits by construction.
   - The off-diagonal is `W²`'s either way, so `jl = jr(−φ)` and `rminus_jacobians`'s transpose
     stay bit-exact.
   - Below the switch the identity's `1` dominates. There `1 + bW²ᵢᵢ` rounds once at the 1-scale,
     and the rewrite measured 1.0–1.3u against GTSAM's 0.43u on `theta:1e-8` to `1e-4`.

## Rationale

These are `cargo xtask conformance` rows at both precisions, judged by `cargo xtask envelope`
against the same oracle rows on `main` and on this record.

| stratum | `main` | this | best oracle |
|---|---|---|---|
| `so3_jr`/`jl` `pi-1e-1` | 3.111 | 2.397 (unresolved) | 1.867 (Sophus) |
| `so3_jr`/`jl` `pi-1e-5` | 2.805 | 1.713 | 1.933 (GTSAM) |
| `so3_jr`/`jl` `pi-1e-8` | 4.097 | 2.455 | 2.363 (sophus-rs) |
| `so3_jr`/`jl` `pi-1e-12` | 2.768 | 1.745 | 1.660 (Sophus) |
| `so3_jr_inv`/`jl_inv` `pi-1e-6` | 1.283 | 1.104 | 0.896 (Sophus) |
| `so3_gamma2` `pi-1e-3` | 1.459 | 1.091 | 1.091 (GTSAM) |
| `so3_gamma2` `pi-1e-8` | 1.416 | 1.001 | 1.070 (GTSAM) |

**The envelope:** 19 domination failures on `main` become **7**, all of them outside this record:
`sen3_log` five and `so3_log` two. No failure is added, resolved or gross.

**Reach.** The coefficients are shared, so 737 of 4424 rows move, in 35 ids, and every one of
those ids holds `a`, `b`, `c` or `d`; 46 of 116 determinism digests change. The functions'
worst strata fall almost everywhere:
- `so3_jr` 4.10 → 2.46;
- `sen3_jl_n1` 8.90 → 3.80 and `sen3_jl_inv_n1` 10.29 → 5.00;
- `sim3_exp` 8.17 → 4.69 and `sim3_log` 4.88 → 2.94 (`f64`);
- `coeff_c` 94.5 → 36.3 and `coeff_b` 5.67 → 3.57 (value and derivative).

The rises, none of them a resolved or gross loss:

- **`sim3_jr` and `sim3_jr_inv` at `sigma:-700/theta:3`**, 138.8 → 176.9 and 94.1 → 120.0.
  - **Attribution.** It is `d`'s formula alone. Restoring `d`'s old closed arm in `sinc_a_b_d`
    (Sim(3)'s group, read by nothing else) gives 137.4 and 93.2; restoring `a` or `b` changes
    nothing; and `b`'s move from series to closed arm at `θ² = 9` is not it.
  - **Where `d` enters.** Only in the derivative block, `da = −½(b − 2d)`. Every candidate's `da`
    is within 1.8 `u` of mpmath at `θ = 3`.
  - **Mechanism.** The assembly at `σ = −700` amplifies sub-ulp coefficient differences
    hundredsfold; `main` already reads 138.8 there. Forming `da` instead as `(2a − sin θ/θ)/θ²`,
    0.40 `u` and the most accurate, reads 402 on the deep tier, against 401 for `(1 − 2a)/(2θ²)`.
    So the reading is which way the roundings correlate with `v₂`'s, not any coefficient's
    accuracy.
  - **Deep tier, 64 records.** `main` reads 328.5 and this record 401.4. The swap test gives
    `p = 0.125`. Of the 122 deep Sim(3) strata, none is resolved worse, and 13 (`sim3_jr`) and 4
    (`sim3_jr_inv`) are resolved better. 49 of the 64 records read worse, but they share 3
    distinct `fl(θ²)`, so that is three draws of one rounding, not 64. The `σ = −700` strata at
    other `θ` do not move.
  - **Verdict.** A same-size reshuffle of an ill-conditioned assembly. The lever is Sim(3)'s
    closed form at large `|σ|`, outside this record.
- **`coeff_c`'s derivative lane.**
  - `f64`: 10 strata fall and 2 rise, `pi-1e-2` 2.17 → 13.62 and `pi-1e-4` 0.45 → 6.96; the
    function's maximum goes 94.5 → 36.3.
  - `f32`: 7 strata rise, to 5.84, against a maximum of 3.93; the value stays at most 1.25 `u`.
  - `c′` near `π` is itself a cancellation, `−2/θ³ + csc²(θ/2)/(4θ)`. Re-associating Sophus's form,
    `(θ/2)(cos/sin)`, reads worse at both precisions.
  - **The consumer does not regress.** Over 240 random rotations, near `π` and on `[2, 3.1]`,
    `∂J_r⁻¹/∂φ` through `Dual<S, 3>` against mpmath reads:

    | | `main` | this |
    |---|---|---|
    | `f64` near `π` | 1.69 | 1.75 |
    | `f64` on `[2, 3.1]` | 3.33 | 2.25 |
    | `f32` near `π` | 1.98 | 1.35 |
    | `f32` on `[2, 3.1]` | 2.36 | 1.83 |

    `SO3::jr_inv` and `SEn3`'s `jr_inv` and `Log` are the only readers of `c`'s derivative, and no
    corpus id differentiates them.
- **Per-stratum rises up to 0.9 `u`** in the SE_N(3) Jacobians.

**What lost:**
- **The cancellation-free diagonal at every `θ`** (`sin θ/θ = 2k cos(θ/2)`): it fails five
  small-angle strata of `so3_jr`.
- **Dividing by `z` without the diagonal rewrite:** `so3_jr` keeps `pi-1e-3`.
- **`c` over `z`:** 7 resolved failures in `so3_jr_inv`.
- **`b` as `(θ − sin θ)/(zθ)`:** `so3_gamma2` fails `pi-1e-7` and `pi-1e-9`.
- **Sophus's `(1 − cos θ)/z` for `a`:** `so3_jr` fails `pi-1e-7` and `pi-1e-9`.
- **`J_r⁻¹`'s diagonal as `(θ/2) cot(θ/2) + cφᵢ²`,** the same rewrite on `c`'s closed arm.
  - `so3_jr_inv` near `π` rises at every stratum: `pi-1e-6` 1.104 → 1.670 and `pi-1e-10`
    1.107 → 2.066, against Sophus's 0.896 and 0.822.
  - `c` divides by `θ̂·θ̂`, not by `z`, so `(θ/2) cot(θ/2)` is not `1 − cz`. That mismatch costs as
    much as the cancellation it removes, and the bar is normwise.

**Cost**, by x86-64 instruction count of `#[inline(never)]` wrappers, against `main`:

| function | instructions | divisions | multiplies | add/sub | jumps |
|---|---|---|---|---|---|
| `jr` | 181 → 189 | 2 → 2 | 39 → 36 | 36 → 34 | 6 → 6 |
| `jr_inv` | 139 → 135 | 1 → 2 | 28 → 28 | 26 → 26 | 4 → 4 |
| `gamma2` | 151 → 174 | 1 → 3 | 28 → 36 | 24 → 35 | 4 → 5 |

No `sqrt`, call or stack spill is added. The added divisions are on the exact arms, which now run
from `θ² ≥ 5.83` to `8.35` rather than `8.35` to `9.65`.

`jr`'s row is the operand-blend layout below. The shipped layout, which selects the diagonal's two
sums, reads 195 instructions, 39 multiplies, 38 add/sub and 9 jumps, and is the faster of the two.
`d`'s closed arm has since dropped two multiplies at the same bits (*Decision* 1).

**The latency gate,** `cargo xtask bench-gate --against` `main`'s binary, pinned to one core, each
row past its own concurrent A/A floor:
- `jr`: `f64` 1.06× near the identity, 0.93× generic, near `π` inside a 0.26 floor; `f32` 1.07×
  near the identity and generic, 0.97× near `π`;
- `jr_inv`: 1.03× at `f64` near `π`, else within noise;
- `gamma2`: 1.03–1.04× near `π`;
- `gamma_apply_jacobians::<2>`: `f32` 1.08× near `π`; `::<1>` within noise;
- `rminus_jacobians`: unchanged.

These microbenchmarks read codegen as much as arithmetic, so each layout that failed was
disassembled:
- **Blending the diagonal's operands** instead of its two sums: three jumps fewer, but `jr` reads
  1.09–1.13× near the identity. Selecting the sums stays.
- **A select-free diagonal,** `t + (vφᵢ² + uW²ᵢᵢ)` with the selects on scalars: the same bits, and
  1.17–1.29×.
- **A `branch` between `poly` and the rewritten matrix:** 2.27×; each arm's matrix goes through the
  stack.
- **`grouped` written as `grouped_with` with no extra value:** `gamma_apply_jacobians::<1>` at
  `f64` reads 1.30–1.58×. Out of line, `jr_coeffs_sinc::<Dual<f64, 1>>` loaded its argument with
  one 16-byte load across the caller's two 8-byte stores, which stalls store-to-load forwarding.
  Forcing it inline moved the stall to `f32` (1.58×). `grouped_with` stays a second copy of
  `grouped`'s branches, and `jr_coeffs_sinc_is_b_s_own` pins its `(a, b)` to `jr_coeffs`' bits.

## Consequences

- **Bits move** for SO(3)'s, SE_N(3)'s and Sim(3)'s Jacobians and logs near `π`, and for `Γ₂`.
- `generated.rs` moves:
  - binary64: `a` 8.35 → 5.83, `b` 9.65 → 7.77, `c` 3.92 → 3.65, `d` 9.65 → 8.35;
  - binary32: `b` 7.77 → 8.98.
- `NUMERICS.md` §3.5, §4 and §7 state the new arms and switches. `maths/coefficients.md` CO.6
  says the shipped arms are no longer its rows and why that is safe above the switches.
- **The seeded kernels take the same arms,** so `the_shipped_arms_are_the_seeded_arms_bit_for_bit`
  holds. Their own pins move with them:
  - the exact-arm bit pins at `z = 0.5` and `2` (`b`, `c`, `d`; `e` unchanged);
  - the D12 prior's per-coefficient maxima (`c` 2.14e3 → 1.11e3, `d` 2.83e3 → 2.73e3, `b`
    366 → 402);
  - `b-no-series`' fitted `p`, still detected: 1.937 → 1.919 at `f64`, 2.134 → 2.086 at `f32`;
  - the seeded SE(3) kernel's maxima (`sen3_jl_n1` 8.903 → 3.802);
  - `0038`'s `W²`-twin tally, 21/104 → 36/94 strata, every id still net worse;
  - the D12 prior now beats `helicoid` on one stratum, `c` at `f32` at `θ = π − 0.1`, in the
    derivative lane, where it beat it on two.
- **`gamma_apply_jacobians::<1>`** takes the same diagonal, so its third output stays `gamma1`
  to the bit.
- **The laws' `probe` inverts `J_r` at a tangent of sines,** `θ ≤ √3`, inside `Jac::inverse`'s
  domain, and `Ad_Y` beside it.
  - At `θ ~ 1e63`, `J_r` reads `nnᵗ`, rank 1.
  - The old diagonal's rounding residue had kept its determinant nonzero by accident.
- **`gamma_apply_jacobians::<1>`'s value keeps `v + a φ×v + b φ×(φ×v)`**, which cancels like the
  old diagonal near `π`. No corpus id scores that value, so it is not rewritten (`0006`).

## Implementation plan

1. This record, the kernel's arms, `jr_coeffs_sinc` and `HatSq::poly_far`, the regenerated
   switches and `NUMERICS.md`. Verified by:
   - `just thresholds-check`, `just test`, `just lint` and `cargo xtask conformance --self-test`;
   - `cargo xtask envelope` against `main` on the same oracle rows;
   - `cargo xtask bench-gate --only so3/jr` (*Cost*).

## Open questions

None. The latency gate and the two rises are measured above. Sim(3)'s conditioning at large
`|σ|` is the lever for `sigma:-700/theta:3` and is not this record's. Answered by the owner on
2026-10-14: the gate's 1.06–1.07× on `jr` near the identity is accepted.
