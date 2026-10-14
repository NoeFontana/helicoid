# 0089: The closed arms take θ² as given, and `J_r`'s diagonal does not cancel

**Status:** draft
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
   - `d = (1 − 2a)/(2z)`.

   `c` becomes Sophus's `(1 − (θ/2) cot(θ/2))/(θ·θ)`. `e` is unchanged. Every group that holds a
   member computes it the same way, so each member is still its own group's to the bit. No new
   transcendental is introduced: each arm takes the root and the sines it took before.
2. **The switches are the sweep's,** regenerated against these arms by `cargo xtask thresholds`.
   No switch is typed.
3. **`SO3::jr` takes `sin θ/θ + bφᵢ²` on its diagonal where `b` is on its closed arm,** with the
   `sin θ/θ` that arm formed `b` from (`jr_coeffs_sinc`), and `1 + bW²ᵢᵢ` below it.
   - The choice is a `select` on `b`'s own mask (`HatSq::poly_far`), so there is no new switch.
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

The rises:
- `sim3_jr` and `sim3_jr_inv` at `sigma:-700/theta:3`, 138.8 → 176.9 and 94.1 → 120.0. At
  `θ² = 9`, `b` moves from its series to its closed arm. That is not attributed further, and
  `σ = −700` is the extreme of Sim(3)'s scale range.
- `coeff_c` at `theta:pi-1e-4`, 0.45 → 6.96, and `coeff_c@f32` at `pi-1e-1`, 1.46 → 4.85. These
  are the derivative lane of a coefficient no `Dual` caller differentiates.
- Per-stratum rises up to 0.9 `u` in the SE_N(3) Jacobians.

None of them is a resolved or gross loss.

**What lost:**
- **The cancellation-free diagonal at every `θ`** (`sin θ/θ = 2k cos(θ/2)`): it fails five
  small-angle strata of `so3_jr`.
- **Dividing by `z` without the diagonal rewrite:** `so3_jr` keeps `pi-1e-3`.
- **`c` over `z`:** 7 resolved failures in `so3_jr_inv`.
- **`b` as `(θ − sin θ)/(zθ)`:** `so3_gamma2` fails `pi-1e-7` and `pi-1e-9`.
- **Sophus's `(1 − cos θ)/z` for `a`:** `so3_jr` fails `pi-1e-7` and `pi-1e-9`.

**Cost**, by x86-64 instruction count of `#[inline(never)]` wrappers, against `main`:

| function | instructions | divisions | jumps |
|---|---|---|---|
| `jr` | 181 → 192 | 2 → 2 | 6 → 9 |
| `jr_inv` | 139 → 135 | 1 → 2 | 4 → 4 |
| `gamma2` | 151 → 174 | 1 → 3 | 4 → 5 |

No `sqrt`, call or stack spill is added. The added divisions are on the exact arms, which now run
from `θ² ≥ 5.83` to `8.35` rather than `8.35` to `9.65`. `jr`'s three jumps are the diagonal's
`select`, lowered to a branch on a mask that is constant per call. The latency gate is owed (open
question 1).

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
- **The laws' `probe` inverts `Ad_Y`, not `J_r`.** `J_r`'s determinant `2(1 − cos θ)/θ²` reads 0 at
  `θ ~ 1e132`, which `Jac::inverse`'s domain excludes. The old diagonal's rounding residue had
  kept it nonzero by accident.

## Implementation plan

1. This record, the kernel's arms, `jr_coeffs_sinc` and `HatSq::poly_far`, the regenerated
   switches and `NUMERICS.md`. Verified by:
   - `just thresholds-check`, `just test`, `just lint` and `cargo xtask conformance --self-test`;
   - `cargo xtask envelope` against `main` on the same oracle rows;
   - `cargo xtask bench-gate --only so3/jr` (owed).

## Open questions

1. **Latency.** `jr`, `jr_inv` and `gamma2` against `main` on a quiet host.
2. **`sim3_jr` at `σ = −700`.** Whether the rise is `b`'s arm or its formula, and whether it is
   worth an exception or a fix.
