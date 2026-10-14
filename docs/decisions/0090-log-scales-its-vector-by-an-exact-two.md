# 0090: `Log` scales its vector by an exact two

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

Under [`0088`](./0088-the-domination-bar-must-tell-worse-from-equal.md)'s bar, `so3_log` lost two
strata at `f64` against sophus-rs, both resolved:

| stratum | `helicoid` | best oracle | p |
|---|---|---|---|
| `theta:1e-3` | 2.376 `u` | 1.631 (sophus-rs) | 6.1e-5 |
| `theta:dense` | 2.586 | 2.161 (sophus-rs) | 6.1e-5 |

Draft [`0032`](./0032-domination-charges-helicoid-for-d16.md) read `so3_log`'s failures as
`libm::atan2`'s cost (D16). That is not what loses these strata:
- **They are lost on `r`'s series arms.** Every one of `helicoid`'s worst records sits there
  (`s = n²/w²` from `1e-7` to `0.02`), where no `atan2` is called. On those records all three
  oracles agree with each other, at 0.75, 0.47 and 0.19 `u` where `helicoid` reads 2.38, 2.35 and
  2.00.
- **The cause is the arrangement.** `φ = r u` with `r = (2/w)·f(s)` rounds three times, at
  `2/w`, `·f` and `·u`, and `r = θ/sin(θ/2)` lies just above `2`. Just above a power of two a
  rounding costs a whole `u` relative; just below, half of one. sophus-rs computes
  `2·atan2(n, w)/n · u`, whose intermediates fall in arbitrary binades.
- **Where D16 does cost.** On the exact arm (`θ ≳ 0.58` at `f64`) `helicoid`'s arithmetic is
  sophus-rs's and `tf_tree_math`'s. There `helicoid` equals `tf_tree_math` on the same `libm`,
  and sophus-rs is better through glibc's `atan2`: `theta:dense`'s record at `θ = 1` reads 2.317
  `u` for both `libm` programs against 1.368 for sophus-rs. `0032`'s attribution holds for that
  arm alone.

## Decision

1. **`coeffs::log_ratio_scale(n², w, sign, u) -> sign · r u`**, with `log_ratio`'s arms, masks
   and `s = n²/w²`; `log_ratio` itself is unchanged.
   - **Series arms:** `t + (t/w)·(d + s·τ(s))`, with `t = 2 sign u` and `d = 1 − w`, where `τ` is
     the series less its constant term. `t` is exact, and so is `d` for `w ∈ [½, 2]` (Sterbenz).
     `d + s·τ` is `O(θ²)`, so the second term's roundings are `O(θ²)` of `u`, and each component
     rounds once at full size.
   - **Exact arm:** `r · (sign u)`, the old bits. There the quotient `atan2(n, w)/n` rounds
     whatever the arrangement, so a split would add a rounding without removing one: measured, it
     lost `theta:pi-1e-2` to `tf_tree_math`.
2. **`SO3::log` passes `w`'s sign flip in**, and each arm applies it where it is free: on `t`'s
   factor `2 sign`, and on `u` beside `atan2`.
3. **No new switch, no new coefficient.** The tail Horner runs over `R`'s committed terms from the
   second on, and `coeffs/generated.rs` is unchanged.

## Rationale

`cargo xtask conformance`, all ids at both precisions, against `main` (with `0089`). Maxima in `u`:

| stratum | `main` | this | best oracle |
|---|---|---|---|
| `so3_log/theta:1e-3` | 2.376 | **0.844** | 1.631 |
| `so3_log/theta:dense` | 2.586 | **2.317**, now its exact-arm record (unresolved) | 2.161 |
| `so3_log/theta:1e-2`, `1e-4` to `1e-8` | 1.80–2.63 | 0.69–0.91 | 1.59–2.09 |
| `sen3_log_n1`, `_n2`, `_n3` at `rho:1e-3/theta=1e-1` | 2.031 | **0.392** | |
| `sen3_log_n1`, `_n2`, `_n3` at `rho:1e-6/theta=1e-1` | 1.995 | **0.591** | |

- **The envelope's domination failures fall from 7 to 1.** It loses `so3_log`'s two and four of
  `sen3_log`'s five, which inherit `φ`, and adds none. The one left is
  `sen3_log_n2/theta:pi-1e-2`, which is `ρ`'s, not `φ`'s.
- **At `f32`** there is no `so3_log` stratum, so this was measured directly: 200 000 rotations per
  band against a binary64 `Log` of the same binary32 quaternion. Over `θ ∈ [10⁻³, 1)` the maximum
  falls from 2.40–3.00 to 0.98–1.01 `u` and the mean halves. `[1, 1.13)` reaches the exact arm,
  whose maximum, 3.09, is unchanged.
- **Every moved stratum.** 119 of 3 055 move: 97 fall and 22 rise. The rises:
  - the three SE(3) charts at `rho:1e-3/theta=1e-1` and `rho:1e-6/theta=1e-1`, 6-record strata,
    by 0.14 and 0.06 `u`; `sen3_log` at the same strata falls by 1.4–1.6 `u`;
  - `sim3_log` at `ray:0.3` and `theta:0.3`, by at most 0.15 `u`;
  - `gaussian_mahalanobis_se3` and `_se23`, ill-conditioned by design (maxima 3·10⁴ to 4·10⁶ `u`),
    by at most 3·10⁻³ relative;
  - ten more, each under 0.06 `u`.

  In each, `ρ`'s own rounding reshuffles which record is the maximum when `φ`'s last bit moves.
  None fails `0088`'s bar.
- **Latency.** A consumer loop over 4 096 rotations per band, best of 50, one pinned core of the
  bench host, in ns per `Log`:

  | band | `main` `f64` | this | `main` `f32` | this |
  |---|---|---|---|---|
  | `θ ∈ [10⁻⁸, 10⁻⁶)` | 3.96 | 4.34 (+9.5%) | 3.43 | 3.84 (+12%) |
  | `θ ∈ [0.1, 0.5)` | 10.13 | 10.61 (+4.7%) | 9.73 | 10.19 (+4.7%) |
  | `θ ∈ [0.5, 2)` | 13.18 | 12.53 (−4.9%) | 9.70 | 9.65 |
  | `θ ∈ [π − 0.1, π)` | 13.49 | 12.80 (−5.1%) | 9.73 | 9.11 (−6.4%) |

  The series cost is the fix itself: one multiply and one add per component that `r · u` does not
  have. The exact arm is faster because the flip is off `atan2`'s path. `bench-gate` agrees in
  direction (`so3_log` near the identity 1.12–1.13×); its single-input rows also move ±5% with
  code layout alone, so the loop is the figure.
- **Rejected layouts**, each measured against the shipped one:
  - **`1/w` before `s`, `s = n²·(1/w)²`** (one division fewer): 1.25× and 1.37× `main` near the
    identity. The division then heads the Horner's chain.
  - **A scalar `g = (d + s·τ)/w`, then `t + t·g`:** 1.20× near the identity at `f32`; one more
    operation after the Horner.
  - **The flip in the caller, as `main` has it:** 3–8% slower than this in every band.
  - **`(c, g)` returned, `t + t·g` in the caller** (this record's first form): 1.41× near the
    identity and 1.19× near `π`, where it added `t·0`.
  - **`2xᵢ/w · (1 + s·τ)`, three divisions:** 1.33–1.70 `u` against 0.91–0.98 in a bit-exact
    model of the `theta:1e-2` to `1e-4` strata.

## Consequences

- **Bits.** `cargo xtask determinism` against `main` moves 13 of 116 digests, each `SO3::log` or a
  caller of it:
  - `so3_log`;
  - `sen3_log_n1`, `_n2`, `_n3`;
  - `sim3_log` at both precisions;
  - `se3_decoupled_local`, `se3_screw_local`, `se3_world_local`;
  - `gaussian_mahalanobis_se3` and `_se23` at both precisions.

  `coeff_r`, the geodesics, S² and every other id keep their bits.
- **Scale invariance** holds within rounding on the series arms, where `1 − w` is not homogeneous
  in `q`; `log_is_scale_invariant_across_the_unchecked_domain` now covers both of them.
- **`theta:dense`'s remaining excess is D16's.** It is the exact arm's `atan2`, as above, and
  `0032`'s remedy, if any, is the one to apply there.

## Implementation plan

1. `log_ratio_scale`, `SO3::log`, `NUMERICS.md` §3.2, the tests `log_ratio_scale_is_r_u_{f64,f32}`
   and the series arms in `log_is_scale_invariant_across_the_unchecked_domain`. The new tests
   check:
   - the exact arm is `r u` to the bit;
   - a series arm is within `3u` of `r u` in value and derivative, so they also fail if `R`'s first
     term is not 1;
   - a `−1` sign negates the result to the bit.

   Verified by `just test`, `cargo xtask conformance` and `cargo xtask envelope`, `cargo xtask
   determinism` against `main`, and `just thresholds-check`.

## Open questions

1. **Latency.** `Log` near the identity costs +0.4 ns, about 10%, for its last `u`; generic and
   near-`π` rotations gain 5%. Accepting the trade is the owner's call.
