# 0052: `Real` owes `sin`, and the corpus does not move

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** this record and `Real::sin` land together

> **Step 2 ran.** Against a baseline built in this tree, `bench-gate --bench groups` with
> `--only geodesic` and `--only jr`:
>
> | row | `f64` | `f32` |
> |---|---|---|
> | `so3/geodesic/generic-1` | 0.8387, **1.19x faster** | 0.7610, **1.31x faster** |
> | `so3/geodesic/near-pi` | 0.8680, **1.15x faster** | 0.7540, **1.33x faster** |
> | `so3/jr/near-pi` | 0.8770, **1.14x faster** | 0.9302, 1.08x faster |
> | `se3/jr/near-pi`, `se23/jr/near-pi` | 0.9751, 0.9792 | 0.9756, 0.9844 |
> | the rest of `jr`/`jr_inv` | 0.986 to 1.002 | 0.991 to 1.015 |
> | `so3/geodesic/near-identity` (the provided-body arm) | 1.0046 | 0.9961 |
> | every `se3`/`se23` geodesic row | 0.99 to 1.01 | 0.99 to 1.01 |
>
> The pattern is the one decision 4 predicts: where a cosine was discarded the routine is faster,
> and where none was the row does not move. The `se3`/`se23` geodesic rows are the control — they do
> not route through `SO3::geodesic` — and `so3/geodesic/near-identity` is the second, taking the
> provided body whose coefficients were not changed.
>
> **Two caveats, both about the run and not the change.** The `jr` run's A/A floor reached 0.6734 on
> one row (median 0.0125) because the machine was not quiet: builds and tests overlapped it. So the
> large ratios are resolvable and the small ones are not. And one row, `f32/se3_jr_inv/near-pi`,
> read **1.0147** and failed its 0.0067 floor — on a routine whose own coefficient (`exact_c`) uses
> both halves of its `sin_cos` and was not touched, whose eleven sibling `jr_inv` rows all came in
> at or below 1.00, and whose only changed code (`d_from`, through `q_coeffs`) *removes* a kernel.
> A 1.5% regression from removing work is what a code-layout shift looks like, which is real even
> when the arithmetic is not, so it was re-measured on a quiet machine rather than dismissed.
>
> **The quiet rerun passes the gate and does not clear the row.** All twelve `jr_inv` rows, nothing
> else running: `f32/se3_jr_inv/near-pi` reads **1.0125** with its least pair at 1.0022 against a
> floor of **0.0152**, so it is *inside* the floor and `bench-gate` reports no regression. That is
> not the same as no regression: the point estimate is ~1.01 in both runs, the row's own A/A floor
> moved from 0.0067 to 0.0152 between them, and a sub-2% effect on this row is simply **not
> resolvable** at its noise level on this host. Recorded as unresolved, against 1.08x to 1.33x on
> four rows that are resolvable several times over. `PHASE1.md` §9's quiet machine is what would
> settle it, and this host is not one.
>
> The rest of step 1's verification: `just lint`, `test`, `doc`, `no-std`, `msrv`,
> `thresholds-check` all **0**; `just envelope` **113**; `conformance/results/` byte-identical; and
> the eight exhaustive binary32 shards **8 passed in 138 s** of wall time. That last figure is the
> point of splitting them — the first two drafts of this module shipped a comment claiming
> concurrency that a single looping test did not have, and timed out at nextest's 180 s twice
> before the macro landed.

## Context

`Real` has `sin_cos` and no `sin`, so a caller that wants the sine alone writes `sin_cos().0` and
pays for a cosine nobody reads. Six **shipped** call sites do exactly that:

| site | what it computes | cosine |
|---|---|---|
| `coeffs::kernel::exact_a_b` (`SO3::jr`'s exact arm) | `sin(θ/2)` for `k`, and `sin θ` for `b` | discarded, twice |
| `coeffs::kernel::d_from` (`q_coeffs`'s exact arm) | `sin(θ/2)` for `d` | discarded |
| `SO3::geodesic`'s blend arm | `sin α`, `sin((1−t)α)`, `sin(tα)` | discarded, three times |

Three more are `cfg(test, __sweep)` (`exact_k`, `exact_b`, `exact_cos_half`) and are not on a hot
path; `exact_c`, `exact_b_d_e` and `e_from`'s caller use both halves and are left alone.

**Measured, this host, `libm` 0.2.16, 20 M calls over `(0, π/2]` in three alternating passes:**
`libm::sin` **3.84 ns**, `libm::sincos(..).0` **5.17 ns** — 1.33 ns, **26 %**. Reading the source
says why: both reduce with `rem_pio2`, then `sin` evaluates the one kernel its octant needs
(`k_sin, k_cos, -k_sin, -k_cos` on `n & 3`) where `sincos` evaluates **both** and permutes.

[`0022`](./0022-real-owes-acos-and-cos.md) is `ready` and decided this exact shape for the other
half — "`Real` gains two required methods, `acos` and `cos` … `cubic.rs` loses its private `acos`
and its three discarded sines", at a measured 4.4 ns from `cos`. Its step 1 is still unimplemented,
and `eig3` and `cubic` still each discard a sine. So this record is not a new idea; it is `0022`'s
argument applied to the half `0022` did not cover, for callers `0022` did not have.

## The identity this rests on, and why it is a test

Every figure recorded against a routine re-spelled from `sin_cos().0` to `sin` is valid only if the
two agree **to the bit**. They are not the same expression, and `libm` is a caret dependency. Three
places could part, read from 0.2.16:

1. **The general path agrees by construction.** `sincos`'s permutation
   `(s, c), (c, −s), (−s, −c), (−c, s)` has the same first component as `sin`'s four arms, from the
   same `rem_pio2` and the same kernels.
2. **The small-argument cuts differ.** `sin` returns `x` below `2^-26`, `sincos` below
   `2^-27 √2`; between them one returns `x` and the other a Horner. They agree anyway, but by a
   rounding argument: the correction is at most `x²/6 ≤ 2^-54.58` relative where half an ulp is at
   least `2^-53`, so the sum rounds back to `x`. That band is **reachable** — `tα` at
   `geo:consecutive` lands in it.
3. **At binary32 the octant arms are differently spelled.** Below `3π/4` and positive, `sinf`
   evaluates `k_cosf(x − S1_PIO2)` where `sincosf` evaluates `k_cosf(S1_PIO2 − x)` — exact
   negations. They agree only because `k_cosf` reads its argument solely through `x * x`. A kernel
   with a term linear in `x` would part them and nothing promises one will not appear.

So: **exhaustive over all `2^32` binary32 bit patterns, 0 disagreements** (3m47s), which settles
item 3 for good rather than sampling around it; and 0 disagreements over **520 018** binary64 arguments —
the `(0, π/2]` range, 200 000 points inside item 2's band, and 800 000 seeded draws across sixty
decades, both signs, past every reduction threshold.

**And the corpus does not move.** All 24 files under `conformance/results/` are byte-identical to
the pre-change run except for the git revision column, which gained `-dirty`; strip it and the
digests of `helicoid.csv` and `seeded-correct.csv` match exactly. That is `0047`'s standard and it
is the strongest statement available here: this change has no numerics, only cost.

## Decision

1. **`Real` gains one required method, `sin`**, routed through `libm::sin`/`libm::sinf` like every
   other transcendental (D16). `PHASE2.md` §2's trait block gains the line.

2. **It is NORMATIVE that `sin` is bit-identical to `sin_cos().0`.** This is a requirement on the
   impl, not an observation about today's `libm`, and `sin_tests` enforces it at both precisions —
   exhaustively at binary32. A `Real` impl that cannot meet it must not implement `sin` by a
   different route.

3. **`Dual::sin` goes through `sin_cos`.** A dual's derivative needs the cosine, so there is nothing
   to save, and routing it through `sin_cos` makes the value identical by construction rather than
   by a test. The saving is at `f64` and `f32`, which is where the hot paths are.

4. **Every shipped site that discards its cosine takes `sin`** — the six in the table above, and no
   others: a `cfg(test, __sweep)` site is not worth the churn, and a site that uses both halves is
   not a site. `coeffs::kernel::half_angle_sin` is the one spelling of the half-angle argument, so
   the `sin`-only form cannot drift from `half_angle`'s.

5. **The transcendental tally gains a `sin` column.** `coeffs::tests`'s counting scalar pins each
   group's `(sqrt, sin_cos, atan2)`; it becomes `(sqrt, sin_cos, sin, atan2)` and the rows are
   re-pinned: `jr` from `(1, 2, 0)` to **`(1, 0, 2, 0)`** — no `sin_cos` at all — and `q` from
   `(1, 2, 0)` to **`(1, 1, 1, 0)`**. Folding the two columns would hide exactly what this record
   changes, which is cost and not value, and `(1, 2, 0, 0)` for `jr` would say it regressed.

6. **`0022`'s `cos` and `acos` are still owed and are not taken here.** `eig3` and `cubic` each
   still discard a sine, at `0022`'s measured 4.4 ns. Its step 1 also deletes `cubic.rs`'s private
   `acos` and makes `bits_are_omnisacs` a statement about every arm, which is a larger change with
   its own verification; mixing it into this one would couple a trait addition to a bit-identity
   claim against another project.

## Rationale

The alternative is to leave the cosines being computed, which is what `0022` already declined for
the mirror case and on a smaller measured margin per call than this one.

The alternative *implementation* is a `sin` that is not `sin_cos().0` — for instance a cheaper
polynomial on a reduced range. That is refused by decision 2 and the reason is the project's own
verification strategy: the value of this change is that nothing needs re-measuring, and that value
is worth strictly more than any further nanosecond. A faster `sin` that moved a bit would owe a
re-measurement of every row it reaches, which is the whole corpus for `jr` and `q`.

**Why a trait method rather than a free function** over `S::sin_cos`: there is nothing to write. The
saving exists only inside `libm`'s choice of kernel, so it has to be reachable from the impl.

## Consequences

- `Real` is one method wider. Every implementor gains it: the two float impls, `Dual`, and the four
  test doubles (`linalg::tests`'s lane type, `xtask`'s `Host` and seeded-kernel scalars,
  `dualmat_tests`'s probe). `Host` uses `std`'s `sin`, which is the same call its `sin_cos` makes,
  so it is identical there for the same reason.
- `SO3::jr`'s exact arm computes **no cosine it does not use**, and `q_coeffs`'s computes one fewer.
  Those are `so3/jr`, `se3/jr`, `se23/jr` and the `q`-bearing `jr`/`jr_inv` rows at `generic` and
  `near-pi`, plus `so3/geodesic` above its switch.
- A `libm` patch release that re-shaped either function would now fail a test instead of silently
  moving a figure — which is a guard the crate did not have for `sin_cos` either, and gains only
  for `sin`.
- `0022` step 1 is unchanged and still owed; this record makes its case marginally stronger by
  demonstrating the shape end to end.

## Implementation plan

1. This record; `Real::sin` with its five impls; `sin_tests`; the six call sites;
   `coeffs::tests`'s `sin` column and re-pinned rows; `PHASE2.md` §2 — verified by `just lint`,
   `just test`, `just conformance` leaving `conformance/results/` byte-identical, and
   `just envelope` holding at **113**.
2. `cargo xtask bench-gate --bench groups --only geodesic` and `--only jr`, against a baseline built
   in this tree — the figure recorded here whichever way it goes, since decision 2 means accuracy
   cannot have moved and only cost can.

## Open questions

None.

## Further work

1. **`0022` step 1**: `Real::{acos, cos}`, `cubic.rs`'s private `acos` deleted, `eig3`'s and
   `cubic`'s discarded sines removed, `bits_are_omnisacs` extended to the trigonometric arm. Decided
   there, measured there, unimplemented.
2. **`sin_cos` has no bit-identity guard of its own.** `sin_tests` pins `sin` against `sin_cos`, so a
   `libm` change that moved *both* consistently would pass. What would catch that is a pinned digest
   of `sin_cos` over a seeded sweep, as `sqrt_tests` does for `sqrt` — the same shape, and owed for
   the same reason `0018` gave there.
3. Whether `atan2` has a cheaper single-output sibling worth the same treatment. `libm` exposes
   none, so this is a question about the kernel and not about `Real`.
