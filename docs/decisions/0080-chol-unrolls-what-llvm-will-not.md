# 0080: `chol` unrolls what LLVM will not

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #155

## Context

[`0030`](./0030-a-consumer-comparison-is-a-trigger-0021-does-not-list.md) (draft) measured `chol`
2.6–3.4× slower than `nalgebra`'s `Cholesky` at `N = 6`, 1.5× at `N = 3` and 1.3–1.4× at `N = 8`,
in an ad-hoc harness whose control read up to +11%, and asked why the gap peaks at `N = 6`
(question 2). [`0021`](./0021-the-cholesky-finite-guard-stays.md) rule 1 sets the bar for rewriting
the kernel: at least 15% at `N = 6` on a pinned core, on the success path, no regression on failing
input or on `Dual`, and `(L, mask)` bit-identical on every input (NaN sign and payload aside). Its
rule 3 trigger has fired: `just bench` exists and `crates/helicoid/benches/linalg.rs` is the
`helicoid-linalg` bench binary (`0053`).

**Where the time goes.** Disassembled in an `#[inline(never)]` wrapper, `chol::<f64, 6>` is a
loop: one `sqrtsd`, a `memset` of the 288-byte column array, a `memcpy` of it into the result,
two `panic_bounds_check` paths, and the inner sums partially unrolled by four with remainder
branches. Every stored entry goes through the stack before the next column reads it. The loop
nest is not fully unrolled because its cost exceeds LLVM's full-unroll threshold: rebuilt with
`-C llvm-args=-unroll-threshold=T`, `chol` becomes straight-line (one `sqrtsd` per column) at
`f64` and `f32`, `N = 6`, from `T = 1000`, and at `N = 8` and `Dual<f64, 2>`, `N = 6`, from
`T = 1500`; at the default only `N <= 3` (`f64`, `f32`) is unrolled. That is the answer to `0030`
question 2: `N = 3` is straight-line code and `N = 6` is not. `nalgebra`'s `N = 6` factor is
straight-line, vectorized two lanes wide.

The mask design is not the cost (`0021` measured deleting the guard at −3% to −6%): fully
unrolled by the threshold flag, the shipped source ran `f64`, `N = 6` at 65 ns against 98 ns in a
scratch harness.

## Decision

1. **`chol` runs the recurrence fully unrolled for `4 <= N <= 9`.** `macro_rules!` in
   `chol.rs` expand one column per literal index: `rows!` writes `sum_{k<j} L_ik L_jk` for every row
   `i >= j` of the column, each from its first term, left to right (`+0` when empty), and `columns!`
   writes column `j` for the rows below it,
   each step guarded by `j < N`, `i < N`, constants after monomorphization. The operations, their
   operands and their order are those of the loop: the pivot `a_jj - dot`, `ok`, the fold into the
   mask, `sqrt(select(ok, d, 1))`, each entry's quotient, its `fin`, the fold `fin or not ok`, the
   store `select(ok and fin, x, 0)`. Output is written into the result matrix in place.
2. **The loop stays**, as `chol_loop` (`pub(crate)`), unchanged: it runs `N <= 3` and `N > 9`, and
   it is the reference the unrolled form is tested against.
   - **`N <= 3`** because LLVM already unrolls the loop there at `f64` and `f32`, and also
     unswitches each column on the loop-invariant pivot mask, so an entry's store is a single
     `and` with `fin`. The unrolled form has no loop to unswitch and broadcasts `ok` into the
     vector mask instead: `f32`, `N = 3`, positive definite measured **1.116×** (8.42 ns against
     7.55), failing in all six pairs against a floor of 0.0052. With the loop kept there, the
     `chol_of::<f32, 3>` machine code is identical to the baseline's.
   - **`N > 9`** because the unrolled code grows as `N^3`: `chol_of` is 1743 bytes at `f64`,
     `N = 6` (1157 looped), 3313 at `N = 8` (1782), 4628 at `N = 9` (2102), and in a scratch
     build 10051 at `N = 12` and 22035 at `N = 16`. `N = 9` covers every group tangent up to
     SE₂(3); a larger bound is a measurement on a consumer, not a microbenchmark.
3. **Bit identity is tested, not argued**: `chol_unrolled_is_the_loop_to_the_bit` (every special in
   every lane, `f64`, `f32`, `Dual<f64, 2>`, `N = 1..10`) and
   `chol_unrolled_is_the_loop_to_the_bit_on_plain_values` (`B Bᵀ + I` with one lower entry made
   special, so the mask is set and every rounding runs, then fails from a passing column). Both
   compare every lane of `L` and the mask. Two seeded defects — the quotient as
   `x · (1 / L_jj)`, and the entry's fold into the mask dropped — fail both tests.
4. **Bench rows** in `crates/helicoid/benches/linalg.rs`: `chol/{3,6,8,9,12}/{pd,failing}` at `f64`
   and `f32`, `chol/6/{pd,failing}` at `Dual<f64, 2>`. `pd` is `B Bᵀ + N I` from a fixed LCG;
   `failing` negates the pivot `N/2`. The factor and the mask are returned whole through an
   `#[inline(never)]` call, so every entry stays live (`0030`'s trap).

No formula changes: `NUMERICS.md` §15 is untouched, and no corpus row can move. `chol` has no §14
row and gains none; the bit-identity tests are its twin check, as `0069`'s are for `Matrix ×
Matrix`.

## Measurements

`cargo xtask bench-gate --bench linalg --only chol --against <baseline>` (`PHASE1.md` §9,
`0033`): baseline built at `origin/main` 41d28b5 with the bench rows of item 4, candidate this
change; `b₁ c b₂` per benchmark, three replicates a pass apart, pinned with `taskset -c 6` on the
AMD EPYC-Milan host (8 vCPU). Other sessions were benchmarking on other cores during the runs;
each row's own A/A floor is beside it. Times are the median ns per call of the least disturbed
replicate.

Two runs. The **first** was of a candidate that also unrolled `N <= 3`; it failed the gate on
`f32`, `N = 3`, `pd` (1.116×, every pair, floor 0.0052), which is Decision item 2's evidence, and
was otherwise quiet (median concurrent floor 6.75%). The **second** is this change: the gate
passed, "22 benchmarks, none above its own concurrent floor in all 6 pairs", on a noisier machine
(median floor 29.6%, worst 129%). The `chol_of` code sizes for `N >= 4` are equal in the two
candidates. *Baseline ns* is the mean of `b₁` and `b₂`, *ratio* candidate over that, both from the
replicate with the lowest total; *gate point (floor)* is the gate's own statistic and its
concurrent A/A floor for the second run.

| scalar | `N` | input | baseline ns | candidate ns | ratio | gate point (floor) | first run's ratio |
|---|---|---|---|---|---|---|---|
| `f64` | 3 | pd | 9.97 | 9.92 | **0.995** | 0.991 (0.035) | 1.016 (unrolled then) |
| `f64` | 3 | failing | 9.08 | 9.23 | **1.017** | 0.992 (0.309) | 1.022 (unrolled then) |
| `f64` | 6 | pd | 97.47 | 57.29 | **0.588** | 0.833 (0.214) | 0.588 |
| `f64` | 6 | failing | 85.30 | 53.72 | **0.630** | 0.638 (0.163) | 0.631 |
| `f64` | 8 | pd | 121.19 | 72.68 | **0.600** | 0.623 (0.331) | 0.602 |
| `f64` | 8 | failing | 109.42 | 70.91 | **0.648** | 0.649 (0.543) | 0.647 |
| `f64` | 9 | pd | 147.51 | 99.56 | **0.675** | 0.674 (0.272) | 0.680 |
| `f64` | 9 | failing | 142.31 | 93.42 | **0.656** | 0.486 (0.623) | 0.673 |
| `f64` | 12 | pd | 230.22 | 224.15 | **0.974** | 0.962 (1.290) | 0.981 |
| `f64` | 12 | failing | 230.52 | 222.28 | **0.964** | 0.961 (0.607) | 0.973 |
| `f32` | 3 | pd | 7.55 | 7.54 | **0.998** | 0.998 (0.022) | 1.116 (unrolled then) |
| `f32` | 3 | failing | 7.88 | 7.84 | **0.995** | 0.994 (0.296) | 1.000 (unrolled then) |
| `f32` | 6 | pd | 64.22 | 46.65 | **0.726** | 0.801 (0.268) | 0.709 |
| `f32` | 6 | failing | 59.60 | 49.29 | **0.827** | 0.582 (1.021) | 0.727 |
| `f32` | 8 | pd | 82.15 | 60.49 | **0.736** | 0.735 (0.293) | 0.741 |
| `f32` | 8 | failing | 78.34 | 60.43 | **0.771** | 0.750 (0.699) | 0.785 |
| `f32` | 9 | pd | 101.71 | 75.30 | **0.740** | 0.617 (0.629) | 0.733 |
| `f32` | 9 | failing | 101.45 | 72.58 | **0.715** | 0.615 (0.719) | 0.719 |
| `f32` | 12 | pd | 180.07 | 181.63 | **1.009** | 1.003 (0.253) | 1.007 |
| `f32` | 12 | failing | 179.91 | 181.43 | **1.008** | 1.022 (0.691) | 0.996 |
| `Dual<f64, 2>` | 6 | pd | 209.62 | 155.66 | **0.743** | 0.738 (0.060) | 0.738 |
| `Dual<f64, 2>` | 6 | failing | 209.83 | 172.63 | **0.823** | 0.669 (0.113) | 0.702 |

Against `0021`'s bar:

- **15% at `N = 6`, success path**: `f64` 0.588× in both runs (97.5 → 57.3 ns), `f32` 0.726×
  (0.709× in the first). Met.
- **No regression on failing input**: every unrolled failing row is faster, `f64` 0.630–0.656×,
  `f32` 0.715–0.827×; `N = 3` and `N = 12` run the unchanged loop and sit within their floors.
- **No regression on `Dual`**: `Dual<f64, 2>`, `N = 6`, 0.743× (`pd`) and 0.823× (failing; 0.702×
  in the quieter first run).
- **`(L, mask)` bit-identical on every input**: item 3 (cases below).

`N = 12` is the layout control: its `chol_of` is byte-for-byte the baseline's size and its ratio is
within its floor.

Bit identity, release, `PROPTEST_CASES=200000 PROPTEST_RNG_SEED=1`: both tests pass, 2 × 10⁵ cases
each × ten sizes × three scalars.

## Rationale

- **Why not raise LLVM's threshold.** `-C llvm-args` is a build flag, owed by every consumer and
  forbidden in spirit by D16's "no `target-cpu`, no fast-math" (a consumer's profile is not ours to
  set). Source that is straight-line by construction does not depend on a heuristic.
- **Why not a row-major or vectorized-column rewrite.** `0021`'s `rows` lost 37–52% on `Dual`; a
  whole-column left-looking form (every row of every column, a constant trip count, masked) was
  measured here at +50% to +98% in the scratch harness — it divides `N` entries per column where
  the triangle divides `N - 1 - j`.
- **Why not the reciprocal of `L_jj`.** It changes thousands of entries per batch (`0030`), so it
  breaks bit identity, and `nalgebra` divides too.
- **Why not `sqrt` before the `select`.** Bit-identical at `f64` and `f32` (`sqrt(1) = 1`) and a
  few percent at `f64`, but 5–8% slower at `f32` and `Dual` in the scratch harness, and it feeds a
  negative pivot to `sqrt`, against §15.2's stated arm.
- **What remains of the gap to `nalgebra`.** Not measured under the gate (`0030` question 3: it is
  not a subject). Advisory, in a scratch harness on `0030`'s protocol (`nalgebra` 0.35, 2048 dense
  SPD matrices, every entry and the mask checksummed, min of 41, pinned; the `nalgebra` row pays a
  conversion each way): `chol` over `nalgebra` went from 2.00× to 1.53× at `N = 6` and from 0.93×
  to 0.75× at `N = 8`, and stays 1.23× at `N = 3`, whose code is unchanged. The pinned order costs one addition per column on the critical path (`a - s`
  after the finished sum, where `nalgebra` subtracts as it goes), and the guard a multiply, a
  compare and two `and`s per entry, which `0021` keeps.

## Consequences

- `chol` is two code paths with one operation sequence; a change to the recurrence is made in
  `columns!`/`rows!` and `chol_loop` together, and the two bit-identity tests are what holds them
  equal. The unrolled column list is written once (`unrolled!`) and `UNROLLED` is its length, so
  the dispatch bound cannot outgrow the list.
- **Rows of a column accumulate together** (#155's review): `rows!` runs `k` outer and the rows
  inner, so each step reads column `k` of `L` contiguously and broadcasts `L_jk`, each row's terms
  still summed in `chol_loop`'s order. Against the per-row form, pinned on one core: `f32` `N = 8`
  0.93–0.97, `N = 9` 0.97–0.98, everything else at parity within its floor. LLVM's SLP pass already
  packed what it could. The remaining gap to nalgebra at `N = 6` is the column chain: column `j`
  waits on column `j - 1`'s `sqrt` and division, about 35 cycles a column. Shortening it with a
  reciprocal changes bits, so it is not taken here.
- A debug build (`opt-level = 0`) carries the unrolled body in every `chol::<S, N>` it
  instantiates, `N <= 3` and `N > 9` included: stable Rust has no conditional monomorphization on a
  const generic. Release builds drop the untaken path.
- `0021` stands: the guard stays, and its rule 1 is met by this change, not relaxed.
- `0030` (draft): its question 2 is answered here, and its proposed decision 1 ("do not touch
  `chol` yet") is overtaken by a measurement under the instrument; its questions 3 and 4 remain.

## Implementation plan

1. Bench rows, the unrolled `chol`, `chol_loop`, the two bit-identity tests, this record — verified
   by `just test`, `just lint`, `just no-std`, `cargo xtask conformance --fn chol_n{3,6}` and
   `chol_solve_n{3,6}` at both precisions byte-identical to the baseline's results, and the gate
   above.

## Further work

- Re-measure on aarch64 and wasm32 when `just bench` can time them (`0021` rule 3).
- Whether a reciprocal-per-column form (`1/l_jj` once, then multiplies), which moves bits and so
  needs the instrument's no-regress reading, closes the `N = 6` latency gap.
- Whether `N` from 10 to 16 pays in a consumer's loop, where its code size competes for the
  instruction cache.
