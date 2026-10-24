# 0106: `Log` rounds its vector once on the exact arm

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

[`0090`](./0090-log-scales-its-vector-by-an-exact-two.md) made `r`'s series arms round each
component of `φ` once and kept the exact arm as `r · (sign u)`, `r = 2 atan2(n, w)/n`: "there the
quotient rounds whatever the arrangement, so a split would add a rounding without removing one".
That holds for a split. A compensated product removes both roundings, and the exact arm is where
`so3_log`'s largest maxima now sit: `theta:dense` 2.317 `u` (its record at `θ = 1`, against
sophus-rs's 1.368), `theta:1e0` 1.741, `q:nonunit` 1.814, `theta:pi-1e-k` 1.24–1.79.

**Where the exact arm's error comes from.** A scratch attribution (not committed) drew about
$10^6$ exact-arm rotations per stratum and scored them with the conformance metric against a
double-double reference, replacing one stage at a time with its exact value. The median over
blocks of a stratum's size of the block maximum, the stand-in for a corpus maximum, in `u`:

| variant | dense | `1e0` | `pi-1e-k` | nonunit |
|---|---|---|---|---|
| shipped | 1.99 | 2.00 | 1.80 | 1.99 |
| `atan2` correctly rounded | 1.74 | 1.84 | 1.79 | 1.88 |
| `n` correctly rounded | 1.99 | 1.95 | 1.87 | 1.92 |
| `r` correctly rounded, then `r · uᵢ` | 1.42 | 1.37 | 1.14 | 1.32 |
| `r · u` compensated (F1) | 1.52 | 1.59 | 1.56 | 1.63 |
| F1 and `n` compensated (F2, this record) | 1.51 | 1.51 | 1.19 | 1.45 |
| F2 and `0105`'s `atan2` (draft) | 1.19 | 1.23 | 1.19 | 1.30 |
| `φ` correctly rounded (the floor) | 0.76 | 0.81 | 0.73 | 0.80 |

- **The arrangement is the largest share in every stratum.** `r ∈ [2.03, π]` sits just above a
  power of two, so `r` rounded costs up to a whole `u`, and `r · uᵢ` rounded another. Removing
  only the first moves nothing (`theta:dense` 1.99 → 1.83): both must go.
- **`atan2` is second**, 0.1–0.25 `u`, and nothing near `π`, where `libm` is already accurate.
  `0105` (draft, PR #200, not on `main`) measured that an
  `atan2` within 0.503 ulp alone moves `so3_log`'s maxima both ways.
- **`n` matters near `π`**, where `w → 0` and its error passes through `/n` whole; elsewhere it
  nearly cancels between `atan2(n, w)` and `/n`.

## Decision

1. **`coeffs::log_ratio_scale(w, sign, u)` forms `n² = (x² + y²) + z²` itself**, the one sum the
   exact arm recovers the error of, so the precondition is the code's and not a caller's.
   `SO3::log` passes `u` and the flip only. `log_ratio` and its arms are unchanged.
2. **The exact arm is `NUMERICS.md` §3.2's compensated form**, in `coeffs::exact_scale`:
   - `a = atan2(n, w)`, `rh = 2a·(1/n)`.
   - `n_ℓ = ((n² − fl(n)²) + δ(n²))/(2n)`, with `δ(n²)` the exact errors of the three squares and
     the two sums (Knuth's two-sum), and `fl(n)²`'s error exact (Dekker).
   - `r_ℓ = ((2a − rh·n) + n_ℓ(2w/(n² + w²) − rh))/n`, `2a − rh·n` exact (Dekker): the
     first-order change of `2 atan2(n, w)/n` with `n`. `2w/(n² + w²)` divides off `n_ℓ`'s path.
   - `φᵢ = rh uᵢ + (εᵢ + r_ℓ uᵢ)`, `εᵢ` the exact error of `rh uᵢ`: one rounding per component.
     A zero `uᵢ` returns `rh uᵢ`: the correction is `+0` there and would drop a `−0`.
   - The squares and the two sums are `log_ratio_scale`'s own, of `u` before the flip, which
     squares away; only the sums' errors are formed here.
3. **Error-free products without a fused multiply-add (D16):** Veltkamp halves with the splitter
   `2^⌈p/2⌉ + 1`, `2²⁷ + 1` at `f64` and `2¹² + 1` at `f32`, selected by `S::PRECISION`. The
   splitter is the format's, not a fitted value; nothing is typed that `0004` forbids.
4. **The compensation is a branch of its own**, taken for `n ≥ 2⁻⁴⁵⁰` (`2⁻³⁸` at `f32`). Every
   error term, down to `2⁻¹⁰⁶ n²` (`2⁻⁴⁸ n²`), stays normal from `2⁻⁴⁵⁸` (`2⁻³⁹`); the bound is
   that with a margin, the format's too, and no unit quaternion comes near it. Below it the arm is
   `exact_r`'s `(2a/n) u`, rounded as `main` rounds it. The branch is also what keeps LLVM's SLP vectorizer from packing
   the compensation's products with `rh`'s, which measured 14 ns longer in one block. The arm is
   a function of its own, not inlined, given scalars and the squares (*Consequences*, inlining).
5. **`helicoid:r-times-u`**, a planted twin: the shipped `SO3::log` with the exact arm `r · u`,
   dispatched on the masks `log_ratio_scale` branches on (`__sweep::log_ratio_takes_exact_arm`
   shares their one definition), so
   `cargo xtask compare` re-measures this choice (`0078`).
6. **No new switch, no new coefficient.** `generated.rs` is unchanged: the sweep reads `exact_r`,
   not `log_ratio_scale` (`just thresholds-check`).

## Rationale

**`so3_log`, `cargo xtask compare --a helicoid:r-times-u --b helicoid --fn so3_log`.** 921 of
5 124 records better, 71 worse (sign test $p = 2.6\times10^{-189}$); of 30 strata, 17 maxima fall
and none rises. The series strata keep their bits. Maxima in `u`:

| stratum | `r · u` | this |
|---|---|---|
| `theta:dense` | 2.317 | **1.578** |
| `theta:1e0` | 1.741 | 1.680 |
| `theta:1e-1` | 1.664 | 1.298 |
| `q:nonunit` | 1.814 | 1.405 |
| `q:w0` | 1.194 | 0.909 |
| `theta:pi-1e-1` to `pi-1e-12` | 1.24–1.79 | **0.68–0.88** |

Against the five oracles' rows (`cargo xtask conformance --oracle <name> --fn so3_log`), `helicoid`
now has the lowest maximum on every `f64` stratum. `r · u` lost `theta:dense` (2.317 against
sophus-rs's and manif's 2.161), sat 0.08–0.13 `u` above the best at `pi-1e-3`, `-6`, `-10` and `-12`
(Sophus or manif, unresolved by the swap test), and tied the best on the exact arm's other strata.

**Every caller, `just envelope` against `main`'s rows.** 180 of 4 565 maxima move: 150 fall and
30 rise. No domination failure is added; the bar reads none on this branch.
- **Falls**, by caller: `sim3_log` 49, `sen3_log_n1`–`_n3` 52, `so3_log` 17, the SE(3) chart
  `local`s 30, `gaussian_mahalanobis_se23` 2. The largest: `se3_screw_local` at `theta:pi-1e-8`
  4.835 → 2.663 and `pi-1e-2` 5.944 → 3.825, `sim3_log` at `sigma:-1/theta:3` 3.190 → 1.416,
  `sen3_log_n1`–`_n3` at `theta:pi-1e-1` 2.28–2.55 → 0.84–0.92.
- **Rises**, each in a stratum of 4–12 records where `ρ`'s or `σ`'s own rounding reshuffles which
  record is the maximum when `φ`'s last bit moves, as in `0090`: the three SE(3) charts at
  `rho:1e-6/theta=1` 1.805 → 2.389, `se3_screw_local` at `rho:1e3/theta=1` 2.145 → 2.686,
  `sim3_log` at `sigma:1/theta:pi-1e-6` 2.366 → 2.796, and 25 more under 0.46 `u`.

**At `f32`** there is no `so3_log` stratum (`0016`), so it was measured directly, as `0090` did:
200 000 rotations per band against a binary64 `Log` of the same binary32 quaternion. Maximum and
mean in `u`:

| band | `main` | this |
|---|---|---|
| `θ ∈ [0.5, 1)` (mostly a series arm) | 1.038 / 0.397 | 1.038 / 0.397 |
| `θ ∈ [1, 2)` | 3.133 / 0.698 | 2.358 / 0.551 |
| `θ ∈ [2, 3)` | 3.217 / 0.733 | 2.166 / 0.585 |
| `θ ∈ [3, π)` | 2.657 / 0.680 | 1.580 / 0.519 |

**Rejected arrangements**, each measured with `compare` against the twin:
- **`θ · (uᵢ/n)`, uncompensated (F0):** free, but near `π` it is worse than `r · u`:
  `theta:pi-1e-6` 1.255 → 2.266, `q:w0` 1.194 → 2.207; the id's maximum rises 1.13×.
- **`r · u` compensated, `n` not (F1):** every stratum the attribution named falls, but
  `theta:pi-1e-k` keeps up to 1.57 and one stratum rises 1.19×; F2's `n` takes `pi-1e-k` to
  0.68–0.88 and lets none rise.
- **`atan2` alone** (`0105`, draft): moves maxima both ways, as above.

**Cost.** A consumer loop over 4 096 rotations per band, independent calls (throughput) and
calls chained through `w` (latency), best of 50, one pinned core of the bench host with the load
under 2.5, one binary per variant, three interleaved rounds; ns per `Log`:

| band | `main` | F1 | F2 in one block | **F2, this** |
|---|---|---|---|---|
| `θ ∈ [10⁻⁸, 10⁻⁶)`, `f64` | 4.6 / 23.3 | 4.6 / 23.3 | 4.6 / 23.3 | 4.6 / 23.3 |
| `θ ∈ [0.1, 0.5)`, `f64` | 10.9 / 39.8 | 10.9 / 39.8 | 10.9 / 39.9 | 10.9 / 40.0 |
| `θ ∈ [0.5, 2)`, `f64` | 14.9 / 38.8 | 28.0 / 46.5 | 40.9 / 71.1 | 38.9 / 54.2 |
| `θ ∈ [π − 0.1, π)`, `f64` | 13.1 / 37.0 | 27.5 / 44.2 | 41.9 / 71.1 | 40.4 / 54.0 |
| `θ ∈ [0.5, 2)`, `f32` | 10.8 / 36.3 | 16.0 / 39.8 | 20.8 / 51.0 | 20.5 / 44.1 |
| `θ ∈ [π − 0.1, π)`, `f32` | 10.1 / 29.5 | 20.5 / 37.7 | 33.4 / 64.3 | 32.6 / 47.4 |

Throughput / latency. The series arms keep their time. On the exact arm F2 costs 2.6–3.1× `main`'s
throughput and 1.4–1.5× its latency; F1 1.9–2.1× and 1.2×. Against the first branched form
(`f64` 40.4 / 57.6 and 40.6 / 56.9), reusing `log_ratio_scale`'s squares and sums and taking the
`−0` fix cost nothing and saved about 2 ns of latency, and dividing `2w/(n² + w²)` off `n_ℓ`'s path
saved 2 more; both keep the bits on 2·10⁶ draws (`f64`, `f32` and `Dual`), and neither moves
throughput beyond ±2%. Without a fused multiply-add each exact
product is eight operations where `fma` would be one, and one `Log` then executes about twice
`main`'s instructions, `atan2` included, too many for consecutive calls to overlap. Measured and
not taken: the dispatch by `n² < T w²` instead of dividing `s` first (within ±2%), and the
square's error in five operations instead of Dekker's eight (within ±2%).

## Consequences

- **Bits.** `cargo xtask determinism` against `main` moves 13 of 123 digests, the set `0090`
  moved: `so3_log`; `sen3_log_n1`, `_n2`, `_n3`; `sim3_log` at both precisions;
  `se3_decoupled_local`, `se3_screw_local`, `se3_world_local`; `gaussian_mahalanobis_se3` and
  `_se23` at both precisions. Native, `ndebug` and `wasm32` agree on every digest.
- **`Dual` derivatives on the exact arm** carry the compensation's error terms, which vanish to
  first order, so their derivative parts are rounding noise: within 23.5 `u` of `r · u`'s
  derivative at `f64` and 27.7 at `f32` on `log_ratio_scale_is_r_u`'s grid, where the value is
  within 1.3. Whether that is closer to the true derivative is not measured.
- **Inlining.** With the exact arm inline, `SO3::log` grew past what LLVM inlines at a caller
  with several call sites, and an out-of-line `log` loads `(x, y)` as one 16-byte read across the
  two 16-byte stores of a quaternion just written (a non-inlined `compose`'s result, a by-value
  argument): the store does not forward, and consecutive calls stop overlapping. `log` of a
  composed rotation in a consumer loop (`f64`, ns per call, interleaved, one pinned core):

  | θ band | before `0106` | exact arm inline | exact arm out of line |
  |---|---|---|---|
  | `[10⁻⁸, 10⁻⁶)` | 10.5 | 25.3 | 10.4 |
  | `[0.1, 0.5)` | 19.0 | 42.0 | 15–18 |
  | `[0.5, 2)` | 21–29 | 64–67 | 48–49 |
  | `[π − 0.1, π)` | 18–26 | 63–64 | 47–48 |

  (Two rounds on a loaded host; the ranges are the rounds'.)

  Out of line, the arm costs at most 2–3 ns of isolated throughput and none of its latency; the
  bench gate (`so3/log`, and `so3/compose_log`, which times this pattern) finds no row above its
  A/A floor. The bits are unchanged (2·10⁶ draws at `f64`, `f32`, `Dual`; `just determinism`).
  A `log` that a caller does not inline still reads its argument that way, as it did before
  `0106`; making that load forward is not attempted.
- **Signed zeros.** A zero component keeps its sign on the exact arm, as `r · u` kept it, so
  `Log` of the conjugate is `−Log` to the bit there (`log_ratio_scale_keeps_signed_zeros`).
- **Scale invariance** holds within rounding, as on the series arms
  (`log_is_scale_invariant_across_the_unchecked_domain`).
- **The envelope baseline** moves with the maxima above; `--bless` is the owner's.

## Implementation plan

1. `log_ratio_scale`, `exact_scale` and its helpers, `SO3::log`, `NUMERICS.md` §3.2, the test
   `log_ratio_scale_is_r_u_{f64,f32}` (now within `tol` of `r u` rounded twice on every arm,
   value and derivative, the derivative's bound wider on the exact arm, either side of
   `compensated_from` and below it), `log_ratio_scale_keeps_signed_zeros`, and the twin
   `helicoid:r-times-u` with its control test (`the_exact_arm_is_taken_and_moves`). Verified by `just test`, `cargo xtask compare`,
   `just envelope`, `just determinism` (native, `ndebug` and `wasm32` identical) and
   `just thresholds-check`.

## Open questions

None. Answered by the owner on 2026-10-24 by taking the change: the exact arm's cost is accepted
for its accuracy, `n` is compensated (F2), and the callers' rises in 4–12-record strata are
accepted as `0090`'s were. `f32` near `π` is in *Further work*.

## Further work

- **`f32` near `π`** costs 32.6 ns against 20.5 at `θ ∈ [0.5, 2)`, where `main` is flat; a
  subnormal intermediate in `n`'s compensation is the suspect, not yet read.
- With `0105`'s `atan2` (draft) on top, the attribution reads about 1.2 `u` in every stratum.
- **Fewer instructions.** A bit-mask split (one `and`, one subtraction) would replace Veltkamp's
  four operations, but `Real` has no bit access; a lighter `n` compensation (the root's error
  only) is not measured.
- **`Dual` derivatives without the noise.** The compensation runs on the derivative parts too;
  compensating the value alone would leave `rh u̇ + ṙ u` there, but `Real` has no hook that tells a
  value from its derivative.
- **One error-free-transform kernel.** `halves`, Dekker's products and Knuth's sum also exist,
  same-shaped, in the test references (`s2_tests`, `eig3_tests`, `svd3_tests`, `sin_tests`);
  a shared private kernel in `helicoid-linalg`, beside D16's, would be a record of its own.
