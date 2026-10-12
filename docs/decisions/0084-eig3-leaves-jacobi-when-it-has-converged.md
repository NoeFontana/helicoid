# 0084: `eig3` leaves Jacobi when it has converged

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** the PR that lands this record

## Context

[`0081`](./0081-eig3-is-four-jacobi-sweeps.md) made `eig3` four cyclic Jacobi sweeps, twelve
rotations with no convergence test, bought with latency. Pinned on one core of this host's AMD
EPYC-Milan VM, `#154` reads 323–387 ns (`f64`) and 281–327 ns (`f32`) on every fixture of
`benches/linalg.rs`, random spectrum, a pair `1e-8` apart and a triple eigenvalue alike (the gate
below). The coordinator's measurement of Eigen 3.4 on the same host (`SelfAdjointEigenSolver::compute`, the
image's flags, min of 31 runs): 304 / 159 / 22 ns (`f64`), 263 / 50 / 15 (`f32`); nalgebra
0.35: 411 / 258 / 160 and 360 / 129 / 121. Eigen and nalgebra are tridiagonal QR with a Wilkinson
shift and a deflation test: they stop when they have converged, `#154` does not.

Two costs, both measured below:

- **Every rotation is one long chain.** `0081` decision 2: `h = sqrt(d^2 + w^2)`, `den = |d| + h`,
  `g = sqrt(den^2 + w^2)`, `c = den / g`, `s = w / g`, then `c a_rp - s a_rq` is the next rotation's
  `a_pq`. Two square roots and a division sit between one rotation's `a_pq` and the next's.
- **Twelve rotations always run.** A random `f64` spectrum has every off-diagonal entry below
  `u |A|` after 9.4 of them on average (a scratch simulation of `2e5` uniform matrices); a
  multiple of the identity needs none.

## Decision

`NUMERICS.md` §17 states the formulas; the public signature and contract of `eig3` are unchanged.

1. **The rotation.** With `w2 = (2 a_pq)^2`: `h = sqrt(d^2 + w2)`, `den = |d| + h`, `t = w / den`,
   `c^2 = den / 2h`, `c = sqrt(c^2)`, `s = c t`; the diagonal moves by `ta = sign(d) (w2 / 2) / den`.
   The safe argument of `0081` decision 2 stays: below the smallest normal, `1` under the root and
   `c = 1`.
2. **The next entry reaches the next rotation squared.** Each rotation is followed by the pair that
   shares its second index, `(p, q) -> (q, r)`, so after the first rotation `a_rp` is always the entry
   the previous rotation annihilated, zero. The next rotation's entry is `c y` (`y = a_rq`). Its
   square is carried as `w2 = 4 y^2 c^2`, which needs `c^2` and not `c`. The first rotation, of
   `(0, 1)`, is general: `z = t a_02 + a_12`, `w2 = 4 z^2 c^2`, and the entries are `c z` and
   `c a_02 - s a_12`. `d` is recomputed from the stored diagonal.
3. **The exit.** After every rotation, if both entries left are at most `u |A|` (`|A|` the largest
   entry of the scaled matrix, `u = 2^-53`, `2^-24`), the loop stops. The mask decides through
   `Real::branch` and `Mask::decide`: a scalar (`f64`, `f32`, `Dual`) leaves the loop, a lane mask
   runs every iteration with each lane held by the branch. At most twelve rotations, `0081`'s
   bound.
4. **The polish.** The exit is followed by the rotations of the two remaining pairs, both computed
   from the state at the exit, in parallel. Their values move by at most the entries. Their
   `Dual` lanes annihilate the entries' lanes to first order, which the exit, a value test, cannot
   see: for a diagonal input with an off-diagonal lane, the values are already diagonal and only
   this step gives the vectors their derivative.
5. **The order.** The three eigenpairs leave in the frame of an exit after the first rotation,
   which is index order. A multiple of the identity and a diagonal matrix keep their axes. A later
   exit leaves another cyclic order, so `det` is unchanged; the sort makes the values ascending,
   and the order decides only between equal eigenvalues.
6. **The normalisation.** Each column is within about `17 u` of unit length (measured with the
   step removed). One Newton step on `1 / sqrt(n^2)`, `v (3 - n^2) / 2`, replaces `0081`'s root and
   three divisions per column.

## Rationale

**Measured against the alternatives.** Corpus: `cargo xtask conformance --fn eig3 [--precision
f32]`, the largest per-stratum `max_u`. Large sample: a scratch form of `eig3_tests::measure_eig3`
(the same double-double reference and metric, other seeds, `5e4` matrices per family), the same
seeds for every row; the `loguniform` family, the hardest, as p99.9 and mean. Speed: criterion on
`benches/linalg.rs`'s fixtures, `taskset -c 2`, one run each at a load below 1 unless marked; the
gate below is the verdict on the chosen row.

| candidate | corpus max f64 / f32 | `loguniform` p99.9, mean (f64) | random / gap / triple, f64 ns |
|---|---|---|---|
| `#154` | 4.76 / 3.24 | 4.28, 1.18 | 323 / 361 / 323 |
| 2–6, `0081`'s `c = den / g`, `s = w / g` (four divisions) | — | 3.95, 1.17 | 288 / 169 / 169 |
| 2–6, `c^2 = den / 2h`, `s = c t`, `ta = t a_pq` | — | 4.80, 1.19 | — |
| 2–6, `c = den / g`, `s = w / g`, `ta = t a_pq` | — | 5.22, 1.26 | — |
| 2–6, `c = den / g`, `s = c t`, `ta = t a_pq` | — | 5.39, 1.27 | — |
| **this** | **3.06 / 2.67** | **3.84, 1.14** | **228 / 139 / 85** |
| candidate 2: Eigen's QL ported to Rust | — | 11.37, 2.90 | 383 / 187 / 23 (scratch timer) |

- **Exit and polish, not a fixed count.** The exit is what makes ties and near ties cheap. The
  polish is what keeps `Dual` honest: without it, a diagonal input with an off-diagonal lane
  returns vectors with no derivative, because the loop sees a diagonal matrix and leaves.
- **`ta` from `w2`.** Reading the diagonal shift from the same `w2` the root saw, not from `t a_pq`:
  the second and fourth rows differ only there, p99.9 3.95 against 5.22 and mean 1.17 against 1.26.
  `s = c t` against `w / g` is worth less (the fourth and fifth rows).
- **`d` recomputed, not carried.** Carrying `d' = (a_rr - a_qq) - ta` saves an addition and loses
  the cancellation-free difference of a close pair: recomputing it took the `f64` maximum over
  `2e4` matrices from 4.22 to 3.64 (`graded`) and from 2.70 to 1.81 (`gap-1e-10/bottom`).
- **The first rotation's entry `c z`, not `s a_02 + c a_12`.** With `w2` from `z` but `a_pq` from the
  other form, the two disagree to `u |y| / |z|` relatively when `z` cancels, and the rotation that
  follows is not the one that annihilates its entry. That cost up to 30 `u` on a random matrix and
  580 `u` at `f32` (`gap-1e-10/bottom`) while the corpus read 3.05: 32 records per stratum did not
  see it.
- **The exit at `u |A|`.** Tightening it to `u |A| / 16` or `u |A| / 1024` changes no large-sample
  figure beyond its noise; loosening it to `2 u |A|` costs at ties (p99.9 at `f64`: rotated triple
  1.63 → 2.30, `near-triple-1e-16` 1.73 → 2.96) for 19 ns off the noisy-triple fixture.
- **A check before the first rotation, not taken.** Polishing all three pairs of an input already
  within the exit saves the first rotation on a diagonal or noisy-triple input. At a tie the three
  parallel rotations do not commute: `diag-tie+1e-17` went from 0 to 0.60 (p99.9), and every other
  fixture paid 9 ns for the branch.
- **Tracking the frame, not taken.** Carrying the cyclic phase through the loop, to return ties in
  index order after any exit, cost 30 ns at an `f32` triple (52 → 81) and 20 ns on a random
  spectrum. Decision 5 gives index order where it is observable, at a diagonal input, for free.
- **Candidate 2, Householder plus implicit QL with Wilkinson shift and deflation (Eigen's
  `computeFromTridiagonal_impl`, ported).** It wins the triple and `f32`'s close pair, by exiting
  before any rotation. It loses on accuracy in every family that is not exact for both: p99.9
  3.0–11.4 `u` (`f64`) where this reads at most 3.9, as Eigen's oracle rows in `0083` read 3.80–10.8
  `u`. It is also slower on a random spectrum. Not taken.
- **Candidate 3, closed-form seed and one sweep, not built.** The seed alone is the cost the closed
  form had: `0081` measured it at 146–154 ns (`f64`) before any sweep, the polish or the
  normalisation. At ties its vectors are not orthogonal, so it needs a fallback to the full
  iteration.
- **Greedy pivoting (largest entry next), not taken.** A scratch simulation: 9.11 rotations on
  average against 9.44, for a compare and a six-way select on every rotation's chain.

**The chain.** Per rotation, the next rotation's `w2` needs `den`, one division (`c^2`) and one
product. Its `d` needs `den`, one division (`ta`) and two additions. Neither waits for `c`'s root,
which feeds only the vectors and the entry two rotations later; `0081`'s chain ran through a second
root and a division. Divisions and roots per rotation: three and two, as `0081`'s, against four and
two for the second row. The divider is shared between division and root, and dropping the fourth
division was worth 288 → 228 ns at the random fixture, more than the shorter chain alone.

## Consequences

- `NUMERICS.md` §17 states decisions 1–6. `PHASE2.md` §6 and its §0.0 row name this record.
- `eig3`'s rustdoc: the cost and the `Dual` paragraph change, and the contract does not. At a
  repeated eigenvalue the lanes are finite and in general no directional derivative, as `0081` says.
  At a multiple of the identity a diagonal lane gives each eigenvalue one diagonal entry's
  derivative, in index order, and an off-diagonal lane gives the eigenvalues none and the vectors a
  rotation.
- **D16.** The exit is a comparison of values each target computes to the bit, with no
  transcendental, `libm` only through `Real::sqrt`, and no fused multiply-add. So the same rotations
  run on x86_64, aarch64 and wasm32, and the output is bit-identical. A lane mask runs the full
  loop with every lane held where its scalar would have left, so `two_lanes_agree_with_scalars`
  holds bit for bit.
- `committed_linalg.rs`'s 50 `eig3` rows are re-recorded (`measure_linalg`).
- `eig3_tests` gains `measure_eig3` (ignored), the large-sample table below, and
  `a_diagonal_input_differentiates_its_vectors` (decision 4). The recorded `ev`, `res`, `ang`
  bound goes from `32` to `36 u |A|`: the planted proptest's worst `ev` over 10^6 cases, two seeds,
  is 18.0 `u |A|` (`#154`, one seed: 14.7), most of it the rounding of the planted matrix, and the
  header's convention is twice the measured worst. The frame bound stays `16 u` (worst 7.0, was 7.7).

## Implementation plan

1. Decisions 1–6, `NUMERICS.md` §17, the rustdoc, the re-recorded rows and `measure_eig3`, verified
   by `eig3_tests`, `the_unverified_routines_reproduce_their_committed_rows`, `cargo xtask conformance
   --fn eig3` at both precisions against nalgebra's rows and `0083`'s Eigen maxima, and `cargo xtask
   bench-gate --bench linalg --only eig3 --against <#154's bench binary>`.

## Open questions

1. **The large-sample bar per family.** Of the 246 per-family figures of `measure_eig3` (p99.9,
   mean and frame, at each precision), 16 are above `#154`'s (table below). Each is small: frame in
   steps of `u / 2` at `f64`, means at ties up by 0.01 `u` on figures below 0.5 `u`. Both overall maxima
   and both frame maxima are lower. Is "every family at most `#154`'s" the bar, or overall maxima
   plus domination? This record is `draft` until that is answered.

## Measured

**Corpus.** Per stratum, `#154` → this, beside nalgebra 0.35. f64: 11 strata up, 19 down, one equal.
f32: 10 up, 7 down, two equal. Maximum 4.76 → 3.06 `u` (`f64`) and 3.24 → 2.67 `u` (`f32`).

The strata hold 32 records each, 64 for `random`, `rank1`, `scale-*` and `triple`. Their maxima move
by more than these two programs differ: on the large sample, which has `5e4` per family, the same
two programs differ by a few percent in p99.9 and mean. The largest moves are `gap-1e-2/bottom` 4.76
→ 2.10 down and `gap-1e-12/top` 1.08 → 1.65 up. Domination holds on all 50 against nalgebra. It
also holds against Eigen: every stratum is below the smallest Eigen maximum of its row in `0083`'s
table, which is 3.80 `u` at `f64` (`gap-1e-{8,…,12}`) and 3.39 at `f32` (`gap-1e-{3,…,6}`).
`eig:triple` ties at both precisions.

| stratum | `#154` | this | | nalgebra |
|---|---|---|---|---|
| `eig:gap-1e-0/bottom` | 2.3 | 2.12 | down | 8.02 |
| `eig:gap-1e-1/bottom` | 2.4 | 3.06 | up | 6.76 |
| `eig:gap-1e-2/bottom` | 4.76 | 2.1 | down | 5.88 |
| `eig:gap-1e-3/bottom` | 2.65 | 1.84 | down | 14.9 |
| `eig:gap-1e-4/bottom` | 2.41 | 2.28 | down | 11 |
| `eig:gap-1e-5/bottom` | 1.86 | 1.92 | up | 6.07 |
| `eig:gap-1e-6/bottom` | 1.89 | 1.67 | down | 5.81 |
| `eig:gap-1e-7/bottom` | 2.11 | 1.59 | down | 9.4 |
| `eig:gap-1e-8/bottom` | 1.83 | 2.03 | up | 7.97 |
| `eig:gap-1e-9/bottom` | 1.32 | 1.47 | up | 6.3 |
| `eig:gap-1e-10/bottom` | 1.61 | 1.79 | up | 7.85 |
| `eig:gap-1e-11/bottom` | 1.29 | 1.2 | down | 9.08 |
| `eig:gap-1e-12/bottom` | 1.36 | 1.59 | up | 9.19 |
| `eig:gap-1e-0/top` | 2.73 | 2.58 | down | 8.54 |
| `eig:gap-1e-1/top` | 2.8 | 2.45 | down | 8.78 |
| `eig:gap-1e-2/top` | 2.01 | 2.71 | up | 8.67 |
| `eig:gap-1e-3/top` | 2.84 | 2.52 | down | 9.96 |
| `eig:gap-1e-4/top` | 1.85 | 1.79 | down | 8.08 |
| `eig:gap-1e-5/top` | 1.74 | 1.77 | up | 4.9 |
| `eig:gap-1e-6/top` | 2.33 | 1.8 | down | 8.87 |
| `eig:gap-1e-7/top` | 2.13 | 1.71 | down | 5.33 |
| `eig:gap-1e-8/top` | 1.93 | 1.92 | down | 6.66 |
| `eig:gap-1e-9/top` | 1.53 | 1.69 | up | 6.28 |
| `eig:gap-1e-10/top` | 1.4 | 1.19 | down | 7.5 |
| `eig:gap-1e-11/top` | 1.92 | 1.35 | down | 12.7 |
| `eig:gap-1e-12/top` | 1.08 | 1.65 | up | 12.8 |
| `eig:random` | 3.53 | 2.76 | down | 8.19 |
| `eig:rank1` | 1.66 | 1.62 | down | 13.3 |
| `eig:scale-down` | 2.3 | 2.99 | up | 7.57 |
| `eig:scale-up` | 3.07 | 2.22 | down | 7.95 |
| `eig:triple` | 3.99e-14 | 3.99e-14 | = | 3.99e-14 |
| `eig:gap-1e-0/bottom@f32` | 2.11 | 2 | down | 6.32 |
| `eig:gap-1e-1/bottom@f32` | 1.82 | 2.45 | up | 6.79 |
| `eig:gap-1e-2/bottom@f32` | 1.96 | 2.24 | up | 6.8 |
| `eig:gap-1e-3/bottom@f32` | 1.66 | 1.99 | up | 6.29 |
| `eig:gap-1e-4/bottom@f32` | 1.56 | 1.27 | down | 7.28 |
| `eig:gap-1e-5/bottom@f32` | 1.39 | 2.39 | up | 8.56 |
| `eig:gap-1e-6/bottom@f32` | 1.72 | 1.72 | = | 10.3 |
| `eig:gap-1e-0/top@f32` | 1.84 | 2.35 | up | 6.27 |
| `eig:gap-1e-1/top@f32` | 1.79 | 2.02 | up | 6.33 |
| `eig:gap-1e-2/top@f32` | 1.58 | 1.81 | up | 6.65 |
| `eig:gap-1e-3/top@f32` | 1.54 | 1.44 | down | 7.17 |
| `eig:gap-1e-4/top@f32` | 1.57 | 2.15 | up | 7.66 |
| `eig:gap-1e-5/top@f32` | 1.81 | 1.3 | down | 10 |
| `eig:gap-1e-6/top@f32` | 1.82 | 1.91 | up | 10.4 |
| `eig:random@f32` | 2.36 | 2.12 | down | 6.15 |
| `eig:rank1@f32` | 1.99 | 1.56 | down | 12.8 |
| `eig:scale-down@f32` | 3.24 | 2.24 | down | 8.46 |
| `eig:scale-up@f32` | 1.99 | 2.67 | up | 8 |
| `eig:triple@f32` | 0 | 0 | = | 0 |

**Large sample** (`measure_eig3`, `5e4` per family and precision; p99.9 / mean / frame worst,
`#154` → this). Overall maximum 5.86 → 5.33 `u` (`f64`), 5.30 → 4.83 (`f32`). Worst frame 7.5 → 6.0
`u`, 7.72 → 5.90. The 16 figures above `#154`'s:

| family | figure | `#154` | this |
|---|---|---|---|
| `loguniform` | frame f32 | 5.152 | 5.645 |
| `bottom-1e-1` | frame f32 | 5.010 | 5.129 |
| `top-1e-1` | frame f32 | 4.937 | 5.902 |
| `near-triple-1e-1` | frame f64 | 5.500 | 6.000 |
| `near-triple-1e-2` | p99.9 f32 | 2.465 | 2.468 |
| `near-triple-1e-4` | frame f64 | 5.500 | 6.000 |
| `near-triple-1e-8` | p99.9 f32 | 0.137 | 0.149 |
| `near-triple-1e-8` | mean f32 | 0.100 | 0.103 |
| `near-triple-1e-10` | frame f64 | 5.500 | 6.000 |
| `bottom-1e-12` | frame f64 | 4.375 | 5.000 |
| `bottom-1e-14` | frame f64 | 4.000 | 5.000 |
| `near-triple-1e-16` | mean f64 | 0.393 | 0.405 |
| `triple` | p99.9 f64 | 1.883 | 2.058 |
| `triple` | mean f64 | 0.412 | 0.426 |
| `rank1` | frame f64 | 4.000 | 5.000 |
| `diag+1e-10` | p99.9 f64 | 0.010 | 0.029 |

**Cost.** `cargo xtask bench-gate --bench linalg --only eig3 --against <#154's bench binary>`, `taskset -c 2`
(six interleaved pairs per row over three replicates; load 7.5 at the start and 10.3 at the end, from
other sessions). The ratio is this over `#154`, the gate's point estimate; every effect is far
beyond its row's concurrent A/A floor of 0.5–2.7%. The ns are the range of the candidate's runs,
the baseline's in parentheses.

| row | `f64` ratio | `f64` ns | `f32` ratio | `f32` ns |
|---|---|---|---|---|
| random | 0.72 | 228–240 (325–334) | 0.52 | 163–168 (311–327) |
| `gap-1e-8` | 0.40 | 138–149 (358–387) | 0.25 | 71–75 (281–294) |
| triple | 0.28 | 84–95 (324–337) | 0.21 | 59–60 (282–293) |
| the three in turn | 0.47 | 482–492 (1004–1108) | 0.33 | 294–308 (871–916) |

Against Eigen's `compute` on this host, as the coordinator measured it (`-c 4`, min of 31): random
304 / 263, `gap-1e-8` 159 / 50, triple 22 / 15. This is faster on a random spectrum and at `f64`'s
close pair, and slower at `f32`'s close pair and at a triple. nalgebra 0.35: 411 / 360, 258 / 129,
160 / 121.

## Further work

1. Ties still cost rotations: entries of rounding size at a noisy triple shrink by a constant factor
   per rotation, not quadratically, before both fall below `u |A|`, and the polish runs after. Eigen leaves before any rotation, by a test
   that is looser by a factor of 2–3, and is 3.80–10.8 `u` there.
2. A `Dual`-free entry, one without the polish, would save one rotation's latency on every input.
   It is a new public item and would need its own record.
