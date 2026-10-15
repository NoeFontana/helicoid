# 0094: A grid cell is a band

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** this record's PR

## Context

Two families of strata pin a value and vary the rest of the draw:
- **SE_N(3)'s `rho:<scale>/theta=<θ>` cells** (`PHASE1.md` §4.4): one `θ` and 64 axes, also read
  by `se2_*` and the SE(3) charts.
- **Sim(3)'s `sigma:<σ>/theta:<θ>` grid and large-`σ` cells** (`0070`, `0071`): one `σ` and one
  `θ`, with random axes and translations.

**What a pinned value does to the swap test.** The pinned value is one input every draw shares,
so a program's rounding of it is one error that every draw carries:
- `e^σ` is computed from the same binary64 `σ` in all 64 records of a `sim3` cell;
- the `w = cos(θ/2)` of a `rho:*/theta=1` quaternion is the same binary64 in all 64 records.

`0088`'s swap test assumes each draw's pair of errors is exchangeable on its own. It reads one
lucky or unlucky rounding as 64 independent wins. This is `0091`'s twin, across a whole stratum.

**Where the deep tier's losses are.** On `main`, the deep envelope fails 35 strata. 17 are
`sen3_j{r,l}_inv_n*` near `π`, which `0093` (draft) removes. Each of the other 18 pins a value:

| family | strata | example |
|---|---|---|
| `sen3_log_n*`, `rho:*/theta=1` | 7 | `rho:1e-6`: 2.75 `u` against GTSAM's 1.87 |
| `sim3_ad` grid cells | 4 | `sigma:1/theta:1e-4`: 1.28 against Sophus's 0.74 |
| `sim3_exp`, `sim3_log` grid cells | 6 | `sim3_log/sigma:1/theta:1e-8@f32`: 2.25 against 0.90 |
| `sen3_exp_n2/theta:pi-1e-9` | 1 | 2.74 against GTSAM's 1.93 |

None is in a log-uniform stratum.

**`sim3_ad` checked directly.** A binary64 emulation of both programs reproduces their corpus maxima
on `sim3_ad` exactly, and the whole gap is the `sR` block. Over 4000 random `σ`:
- **near the identity**, `helicoid`'s `R·s` is about twice as accurate as Sophus's
  homogeneous `(√s q)` form: at `θ = 10⁻⁸`, mean 0.46 against 0.95 `u` and maximum 1.30 against
  3.27;
- **at the pinned `σ = 1`**, Sophus's rounding of `√e` happens to beat `helicoid`'s rounding of
  `e`, and the corpus reads that one draw 64 times.

## Decision

1. **A grid cell's value is a band.** Each record draws the pinned value from the stratum's own
   stream: the offset from `0`, or from `π`, times `1 + 2⁻¹⁰(2u − 1)` (`gen.strata.band`).
   - **`rho:*/theta=<θ>`:** `θ`, 64 draws on the `theta` stream, of which the committed tier reads
     the first 6. `se2_*` and the charts read them too.
   - **Sim(3)'s grid and large-`σ` cells:** `σ` and `θ`, on a `band` stream.
   - **Unchanged:** `σ = 0` stays exactly `0`, as the formula gives. The rays already vary per
     record, the `theta:pi-1e-k` strata keep their point, and the names stay the same.
2. **`PHASE1.md` §4.4 and its Sim(3) row say so.**

## Rationale

**Why `2⁻¹⁰`.** The band must move a program's rounding of the value on every draw. It must also
keep the regime the cell names.
- At `2⁻¹⁰` the value moves by up to `2⁴²` ulps, or `2·10⁶` at `pi − 1e-6`, so every rounding is
  redrawn.
- **Regimes kept:** `pi − 1e-6` stays within `10⁻⁹` of its gap. `σ = 700` stays below `ln` of the
  largest binary64 (709.78), and `σ = 80` at `@f32` below `ln` of the largest binary32 (88.72).
- **Size:** the corpus grows by 89 bytes, to 65.4 MB.

**Measured.** `cargo xtask envelope`, every subject rescored on both tiers drawn under this record,
against `main` at `aacff0e`:

| | `main` | this |
|---|---|---|
| committed tier: domination failures | 0 | 0 |
| committed tier: unresolved excesses | 23 | 23 |
| committed tier: unresolvable excesses | 616 | 591 |
| deep tier: domination failures | 35 | 19 |
| deep tier: unresolved excesses | 99 | 91 |
| deep tier: unresolvable excesses | 193 | 174 |

- **16 of the 18 pinned losses go, and none is added.**
  - **Gone:** all 7 `sen3_log` losses at `theta=1`, the 4 of `sim3_ad`, and 5 of the 6 of
    `sim3_exp` and `sim3_log`.
  - **Left:** `sen3_exp_n2/theta:pi-1e-9`, not banded. Also
    `sim3_exp/sigma:-3/theta:pi-1e-6@f32`: 2.49 `u` against Sophus's 2.13, still resolved at
    `p = 4.9 · 10⁻⁴`. The 17 inverse-Jacobian losses stay too, in unbanded strata.
- **The program did not change**, so every move is the instrument's.
- **With `0093` as well**, the deep tier fails 2 strata: `sen3_exp_n2/theta:pi-1e-9` and
  `sim3_exp/sigma:-3/theta:pi-1e-6@f32`.
- **The twins read the same way.** `helicoid:homogeneous-r` against `helicoid`, by strata maxima:
  - `sim3_ad` at `f64`: 36 better and 51 worse on the pinned cells (`0076`'s tally), 60 and 38 on
    the bands;
  - `sim3_ad` at `f32`: 43 and 32, then 48 and 31.
  - `0092`'s unit diagonal wins where it lost.

**Rejected:**
- **Counting a pinned cell as one draw.** It is exact, but it leaves 1 draw per stratum, which no
  test can resolve. That throws away the 64 axes, which are independent.
- **Fixing the program to the pinned values.** Over a band of `σ`, `helicoid`'s `sR` is the more
  accurate form. A change that wins at `σ = 1` exactly would lose the regime.
- **Banding `theta:pi-1e-k` too.** Its point is the gap, which is shared by every SO(3), SE(3),
  SE(2) and coefficient id. At `k = 12`, `2⁻¹⁰` of the gap is 2 ulps, which redraws no rounding.
  Only one deep failure (`sen3_exp_n2/theta:pi-1e-9`) sits there, and `helicoid`'s host-`std`
  twin reproduces it, so it is the program's.

## Consequences

- **Corpus bytes:** 42 ids move: the 21 `sen3_*`, the 7 `se2_*`, the 6 SE(3) chart ids and the 8
  `sim3_*`. Every subject's rows over them are drawn again; no baseline is committed.
- **Four `xtask` tests pinned counts on the old draws.** They are re-pinned:
  - `the_task_reproduces_the_tally_of_0038` now asserts by strata maxima, as its doc states. By
    records, `sen3_exp_n1` was net worse by one (19 against 20), and now ties at 19.
  - `the_homogeneous_r_twin_reproduces_the_tally_of_0076`: the tallies above, now per width.
  - `the_committed_corpus_passes_and_the_report_says_how`: the seeded `Exp`'s maximum.
  - `the_swap_has_no_power_on_the_strata_near_pi`: `libm` and the host's `sin`/`cos` now
    disagree on 5 of `sen3_jl_n1`'s 312 records, all at `θ` near 1, where the pinned `θ = 1`
    showed none. Near `π` they still never disagree.
- **`MANIFEST-deep.json` is regenerated.**
- **A stratum's name is its cell, not its value.** A new pinned grid adopts the band, or says
  why not.

## Implementation plan

1. `gen.strata.band` and `CELL_THETAS`; `gen.sim3._cell`; `PHASE1.md` §4.4; the regenerated corpus
   and `MANIFEST-deep.json`. Verified by:
   - the generator's tests, including `test_thetas` and
     `test_a_grid_cell_draws_sigma_and_theta_in_its_band`;
   - `just corpus-check`, `just corpus-deep-check`;
   - `cargo xtask conformance --self-test`;
   - `cargo xtask envelope` with and without `--deep`, every subject rescored.

## Open questions

None of substance; it is a `draft` until the owner accepts the corpus change.

## Further work

- **`sen3_exp_n2/theta:pi-1e-9`**, the program's: the next fix.
