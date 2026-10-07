# 0056: the routines D7 does not reach

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** the four stacked PRs of its *Implementation plan*

## Context

D7 (`PROJECT.md` §5): "A public numeric routine without a corpus stratum is unverified and may not
ship in a release." Five public routines of the two library crates have none, and `0.0.1`
(`0041` decision 7, `0049` decision 5) is the first release:

| Routine | What verifies it today | The open item |
|---|---|---|
| `solve_cubic` | `cubic_tests`: a 46-row inline `mp.polyroots` fixture, planted-root proptests, 36 golden bit patterns | `0031` (draft) open question 1: "does `solve_cubic` get a corpus id and a conformance subject?" |
| `eig3` | `eig3_tests`: a 75-row inline `mp.eigsy` fixture whose generator is not committed (`0053` *Further work* 3), bounds that are fits | `0023` (draft) open question 3: `PHASE2.md` §6's bars name no constant |
| `chol`, `chol_solve` | `chol_tests`: the Higham bounds of `NUMERICS.md` §15.4 as proptests | `0015` (draft) PH.4: "decide with the Phase 2 ids whether `chol` needs strata" |
| `Quat`/`SO3`/`SEn3::renormalize` | `quat_tests` and `sen3_tests` proptests; `SO3::renormalize` has no test and no `# Domain` | `0044` item 5: the `quat_*` family, "owed before `0.1.0`, not before `0.0.1`" |
| `Dual<S, N>` | `dual_tests`: five inline mpmath tables | `PHASE2.md` §0.0: "corpus ids `real_*` pending generator"; `0015` PH.1 |

`coverage::OWED` already names `eig3`, `svd3`, `solve_cubic` and `real_*`. It names nothing for
`chol` or `renormalize`, so coverage reports nothing missing for either.

Three facts shape the instrument:

- **Nothing gates a new id today except finiteness.** The envelope's no-regress bar reads
  `conformance/baseline/helicoid.csv`, which has never been written: `--bless` writes nothing while
  a domination failure stands, and 112 stand. The envelope is binary64-only (`0016` *Consequences*).
- **No oracle answers any of these routines.** Each would be unpaired, so domination says nothing.
- **#104 removed the strata `PHASE2.md` §6 named for `solve_cubic` and the range of `eig:gap-1e-k`**,
  with no record. `0023`, `0031`, `0053` and the `eig3` rustdoc still cite the removed text. This
  record states those strata again, as its own.

## Decision

### 1. Function ids, records and references

Every input is exactly a binary64, or a binary32 in an `@f32` stratum (`0016`). Every reference is
the definition at the input as stored, at 120 digits, rechecked at 150 (`PHASE1.md` §4.2).
`gen.check` cross-checks each record by a property independent of the route that computed it.

| Id | `in` | `out` | Reference | Cross-check |
|---|---|---|---|---|
| `solve_cubic` | `a`, `b`, `c`, `d` | `re[3]`, `im[3]`: all three roots, ordered by `(re, im)` | The planted roots where the coefficients are exact (below). Otherwise `mp.polyroots(…, extraprec=300)` of the stored coefficients | `p(z) = 0` to 100 digits. A planted root is checked exactly in rationals |
| `eig3` | `A` (3×3, symmetric, column-major) | `lambda[3]` ascending, `V` (3×3, columns the eigenvectors) | `mp.eigsy` of the stored `A` | `‖AV − VΛ‖` and `‖VᵀV − I‖` below `10^-100 ‖A‖` |
| `chol_n3`, `chol_n6` | `A` (n×n, symmetric) | `valid[1]` (1 or 0); `L` (n×n, lower) | `valid` decided exactly, by an `LDLᵀ` of the stored `A` in rationals. `L` is `mp.cholesky(A, tol=0)` of the stored `A` (its default absolute tolerance refuses `chol:diag-scale`'s pivots), or zeros when `valid = 0` | `‖LLᵀ − A‖ < 10^-100 ‖A‖` |
| `chol_solve_n3`, `chol_solve_n6` | `A` (n×n), `b[n]` | `x[n]` | `mp.lu_solve` of the equilibrated system `D A D`, `D = diag(A)^-1/2` (`mp.lu_solve` of `A` itself calls `chol:diag-scale` singular) | `‖Ax − b‖ < 10^-100 ‖A‖‖x‖` |
| `quat_renormalize` | `q[4]` | `q[4]` | `q/‖q‖`, the exact projection (see below) | `‖q'‖ = 1` and `q' ∥ q` |
| `real_sqrt`, `real_cbrt`, `real_acos` | `x` | `value`, `d` | `value` is the function. `d` is its derivative, by calculus | `d` against `mp.diff` with a relative step |
| `real_sin_cos` | `x` | `sin`, `cos`, `d_sin`, `d_cos` | the same | the same |
| `real_atan2` | `y`, `x` | `value`, `d_y`, `d_x` | `atan2(y, x)` and its two partials | the same, per argument |
| `real_div` | `n`, `d` | `value`, `d_n`, `d_d` | `n/d` and its two partials | the same, per argument |

- **`quat_renormalize` is scored against the exact projection, not against the Newton step.** The
  claim being scored is `NUMERICS.md` §3.6's: "a normalization only for
  `|‖q‖² − 1| ≤ 2^-26.29` (`f64`) or `2^-11.79` (`f32`)". Its strata stay inside that band.
- **`real_*`'s derivative is the calculus derivative**: `1/(2√x)`, `1/(3 ∛x²)`,
  `−1/√(1 − x²)`, `(cos, −sin)`, `(x, −y)/(x² + y²)`, `(1/d, −n/d²)`. These are textbook
  derivatives, not closed forms of `NUMERICS.md`. They are evaluated at 120 digits, so §4.3's
  "never evaluates a closed form from `NUMERICS.md`" is untouched. `mp.diff` at a relative step
  `h = 10^-70 |x|`, its working precision raised until the difference's rounding is 60 digits below
  `|f'|`, is the generation-time cross-check, to 60 digits: the coefficients' absolute step
  underflows at `x = 1e-300`, and a step of `10^-40 |x|` truncates at `1e-37` (`sin` at `1e22`).
- **No id for `sin`, `cos`, or for `SO3`/`SEn3::renormalize`.** Each is bit-identical to an id
  that has strata:
  - `Real::sin` and `Real::cos` to `sin_cos` (`0052`, `0022`; `sin_tests`);
  - the two `renormalize` methods to `Quat::renormalize`, which they call.

  The two `renormalize` methods get bit-identity tests.
- **`svd3` stays owed** (`PHASE2.md` §6). It is not implemented.

### 2. Strata

**Naming and counts.** A stratum whose name ends `@f32` is binary32, and its name is otherwise the
binary64 stratum's wherever the binary64 one exists; a name never states a value that differs by
precision. The binary64 strata come
first, then the binary32 ones in the same order (`0016`'s layout). Where a range is given per
precision, the `@f32` strata draw their own values from their own streams rather than rounding
binary64 draws. Each stratum has 64 records unless the table says otherwise.

Every draw comes from `gen.rng.stream(seed, stratum, purpose)`, so no new stratum moves an existing
record.

| Id | Stratum | Draw, binary64 (binary32 where it differs) |
|---|---|---|
| `solve_cubic` | `cubic:distinct` | Three distinct planted roots of 6 significant bits sharing one scale `2^e`, `e ∈ [−4, 4]`; `a` of the same kind. The coefficients are then exact at both precisions (8 bits would put `d` at 32) |
| | `cubic:double` | The same, with two roots equal |
| | `cubic:triple` | The same, with all three roots equal |
| | `cubic:one-real` | `a (x − r)(x² + p x + q)` with `p² < 4q`, all of `r`, `p`, `q` planted |
| | `cubic:near-double-1e-k`, k ∈ {2, 4, 6, 8} | Roots `r`, `r(1 + 10^-k)`, `s` at 120 digits, coefficients rounded to the precision. The reference is the stored coefficients' roots: it may be a complex pair |
| | `cubic:one-real-p-small` | `x³ + p x − 1`, `p` log-uniform in `[1e-7, 1e-3]` (`[1e-4, 1e-1]`). This is the cancellation the rustdoc documents |
| | `cubic:coeff-scale-up`, `cubic:coeff-scale-down` | `cubic:distinct`'s cubics with every coefficient times `10^±100` (`10^±20`) |
| `eig3` | `eig:random` | `A = Q Λ Qᵀ`, `Q` from a Haar-random quaternion, `λ` uniform in `[−1, 1]` |
| | `eig:gap-1e-k/bottom`, `eig:gap-1e-k/top`, k = 0…12 (0…6) | The two smallest (or largest) eigenvalues at relative gap `10^-k` of `‖Λ‖`. 32 records each |
| | `eig:triple` | `A = c I` stored directly, `c` uniform in `[−1, 1]`: `Q (cI) Qᵀ` at 120 digits rounds to `cI` plus off-diagonal residue |
| | `eig:rank1` | `Λ = (0, 0, c)` |
| | `eig:scale-up`, `eig:scale-down` | `eig:random` times `10^±70` (`10^±8`): inside the rustdoc's range for the vectors (`‖A‖⁴` normal: `1e-75 < m < 1e75`, `1e-9 < m < 1e9`), the narrower of its two |
| `chol_n*`, `chol_solve_n*` | `chol:spd` | `A = Q D Qᵀ`, `Q` a product of `n` Householder reflections from uniform directions, `D` log-uniform in `[1, 10]`. 32 records |
| | `chol:cond-1e-k`, k ∈ {4, 8, 12} ({2, 4, 6}) | The same with `D` log-spaced from 1 to `10^-k`. 32 records |
| | `chol:diag-scale` | `S A S`, `A` from `chol:spd`, `S` diagonal and log-uniform in `[1e-100, 1e100]` (`[1e-15, 1e15]`). 32 records |
| | `chol:indefinite` (`chol_n*` only) | `chol:spd` with one eigenvalue set to `−10^-2 · max D`. 32 records |
| | (`chol_solve_n*` only) | `b` is uniform on the unit sphere |
| `quat_renormalize` | `renorm:eta-2^-k`, k ∈ {27, 30, 40, 52} ({12, 16, 20, 23}) | A Haar-random unit quaternion scaled to `‖q‖² − 1 = ±2^-k`, the sign alternating, rounded |
| | `renorm:eta-edge` | The same at `±2^-26.29` (`±2^-11.79`), §3.6's edge |
| `real_sqrt` | `x:1e{e}`, e ∈ {−300, −100, −10, −1, 0, 1, 10, 100, 300} ({−37, −10, −1, 0, 1, 10, 37}) | Log-uniform in the decade |
| | `x:subnormal` | Log-uniform in `[1e-310, 1e-309)` (`[1e-40, 1e-39)`) |
| `real_cbrt` | the strata of `real_sqrt` | The sign alternating |
| `real_sin_cos` | `x:tiny`, `x:small`, `x:moderate`, `x:large` | Log-uniform in `[1e-300, 1e-8)`, `[1e-8, 1)`, `[1, 1e3)`, `[1e3, 1e22)` (`[1e-30, 1e-4)`, `[1e-4, 1)`, `[1, 1e3)`, `[1e3, 1e9)`), the sign alternating |
| | `x:near-k-pi/2` | `fl(k·π/2)` for `k` log-uniform in `[1, 1e5]`: one of the outputs is tiny |
| `real_acos` | `x:interior` | Uniform in `(−0.99, 0.99)` |
| | `x:near+1`, `x:near-1` | `±(1 − 2^-j)`, `j` uniform in `[2, 52]` (`[2, 23]`), where the derivative diverges |
| | `x:tiny` | `±` log-uniform in `[1e-300, 1e-8)` (`[1e-30, 1e-4)`) |
| `real_atan2` | `yx:generic` | Angle uniform, radius log-uniform in `[1e-3, 1e3]` |
| | `yx:ratio-1e-k`, `yx:ratio-1e+k`, k ∈ {8, 100, 300} ({4, 15, 30}) | `‖y/x‖ = 10^∓k` in each quadrant |
| `real_div` | `nd:generic` | `n`, `d` log-uniform in `[1e-3, 1e3]`, signs alternating |
| | `nd:wide` | `n`, `d` log-uniform in `[1e-100, 1e100]` (`[1e-12, 1e12]`), so `n/d²` stays in range |

No stratum asks for a non-finite answer: there is no `acos(±1)`, no `atan2(0, 0)` and no division
by zero. A non-finite candidate output fails the conformance run (`PHASE1.md` §5), and those inputs
are the domain tests' job.

### 3. Metrics (`NUMERICS.md` §11)

Every score is a ratio of exact integers, rounded once (`metric::forward_u`'s construction), in
units of `u` of the stratum's precision.

1. **Forward error**, the existing `Rule::Forward`, applies to:
   - `lambda`: `Floor::Tiny`. Its norm is `‖λ‖₂ = ‖A‖_F`, so this is `PHASE2.md` §6's "eigenvalue
     error in `u‖A‖`", norm-wise;
   - `L`, `x`, `quat_renormalize`'s `q`, and every `real_*` field: `Floor::Tiny`, `SignRule::Fixed`.
     A `real_*` record scores the larger of its value and its derivatives.
2. **Eigenvector error, gap-weighted**, for `V`, with column `i`:

   `E_V = maxᵢ min(‖v̂ᵢ − vᵢ‖, ‖v̂ᵢ + vᵢ‖) · gapᵢ / (‖λ‖₂ u)`, where `gapᵢ = min_{j≠i} |λᵢ − λⱼ|`
   comes from the reference.

   - This is §6's "eigenvector angular error in `u‖A‖/gap`", the Davis–Kahan scale: a column whose
     error is at that scale scores O(1).
   - A column whose reference gap is exactly zero weighs 0. That is the right weight, since a
     repeated eigenvalue's vectors are not unique.
   - **This answers `0023` open question 3: the bar's constant is the measured, committed maximum
     (decision 5), not a typed `K`.**
3. **Root-set distance**, for `solve_cubic`. Let `Z` be the three reference roots and `R̂` the slots
   whose mask is set. Then:

   `E = max( max_{z ∈ Z real} min_{r̂ ∈ R̂} |r̂ − z|, max_{r̂ ∈ R̂} min_{z ∈ Z} |r̂ − z| ) / (‖Z‖₂ u)`

   with `|r̂ − z| = √((r̂ − Re z)² + (Im z)²)`.
   - A spurious real double root that stands for a complex pair costs that pair's `|Im z|`.
   - A real root no valid slot answers counts as a distance `‖Z‖₂`, so it scores exactly `1/u`:
     finite, so no-regress sees it and the run does not fail.
   - The slots need no ordering: the rustdoc says they are "not sorted".
4. **Mask agreement**, for `chol`'s `valid`. A reported mask that differs from the reference scores
   `1/u`. Where the reference `valid` is 0, `L` is not scored.

### 4. Oracle: `runners/nalgebra`

The runner has the shape of `runners/sophus_rs`:
- excluded from the workspace, with its own `Cargo.lock` and `deny.toml`;
- `Backend::HostStd`;
- pinned by `--version`;
- answering through `convert::TABLE`.

| Id | Its answer |
|---|---|
| `eig3` | `SymmetricEigen::new`, with the eigenpairs sorted ascending |
| `chol_n*` | `Cholesky::new`, where `None` gives `valid = 0` |
| `chol_solve_n*` | `Cholesky::solve` |
| `solve_cubic` | The companion matrix's real Schur form: each 1×1 block is a valid slot and each 2×2 block a complex pair, which is not. This is `numpy.roots`' route |
| `quat_renormalize` | `UnitQuaternion::renormalize_fast`: the same Newton step, so this is the one row that compares like with like |

It answers binary64 only, until the envelope's binary32 half exists. **Domination is scored as for
any other oracle.** `eig3` is a closed form and nalgebra iterates, so expect it to lose on the
narrow-gap strata. Those failures are the measurement `0023` open question 2, Kopp's hybrid, was
waiting for, and they are not entered in `exceptions.toml`, whose rows are ties (`0046`).

### 5. The gate: committed rows, reproduced bit for bit

`xtask/src/conformance/committed_linalg.rs` holds `(fn, stratum, precision, max_u)` for every
stratum of every id in decision 1. An `#[ignore]`d `measure_linalg` test writes the table, and
`the_unverified_routines_reproduce_their_committed_rows` requires the shipped subject to reproduce
each row **exactly**. The precedent is `measure_geodesic`'s `COMMITTED`.

- **No-regress, today.** The rows are per id, so they do not wait for the envelope's first bless.
- **An improvement is a re-recorded row**, as a regression is a failure: the table does not drift.
- **Cross-target identity.** CI runs the test on x86_64 and aarch64, so for these routines the
  committed rows are D16's bit identity across the two. That is the first slice of
  `just determinism` (`PHASE6.md` §3).

### 6. Coverage and the self-test

- `eig3`, `solve_cubic` and the `real_*` ids move from `coverage::OWED` to the required ids of
  `PHASE1.md` §4.3. `chol_*` and `quat_renormalize` are added to them. `svd3` stays owed.
- The self-test gains a half that plants one defect per new metric into the shipped program's
  output, and requires each to fire on its stratum while the clean program does not:

  | Defect | Fires on |
  |---|---|
  | a dropped valid root | `cubic:distinct` |
  | `V` rotated by `1e-9` about a column | `eig:random` |
  | `chol`'s mask flipped | `chol:indefinite` |
  | half a Newton step | `renorm:eta-2^-27` |
  | a doubled `sqrt` derivative | `x:1e0` |

### 7. Specification edits

- `PHASE1.md` §4.3 gets a row per id, and §4.4 the strata of decision 2.
- `PHASE2.md` §6 gets the strata and bars of decisions 2–3.
- `NUMERICS.md` §11 gets metrics 2–4.
- `NUMERICS.md` §12 gets rows for `solve_cubic`, `eig3` and `Dual::acos`, which `API.md` R6 points
  to and which §12 lacks.
- `SO3::renormalize` gets a `# Domain` section.
- `0016` decision 2 asks for a record before any id other than the coefficients gets an `@f32`
  stratum. This is that record, for the ids of decision 1.

## Rationale

**A record per routine** would repeat the instrument five times. The metrics, the f32 extension and
the gate are each one decision, so they are made once.

**A reference for `quat_renormalize` that is the Newton step at 120 digits** would score the
arithmetic of `q(3 − ‖q‖²)/2` and say nothing about the claim users rely on, which is that the
output is unit. The projection scores the claim. Inside §3.6's band, the step's own truncation error
`(3/8)η²` is under `u`.

**Sorted real roots for `solve_cubic`** fail exactly the cases the rustdoc documents as correct
behaviour: a complex pair reported as a real double root, and a near-double pair lost at `f32`. They
would need those strata removed, against D7's "do not narrow a stratum because it fails". The
set distance scores them by how far they actually are from the roots.

**A plain forward error on `V`** reads `10^12 u` on `eig:gap-1e-12` for a perfect algorithm. A
column-wise angle without the gap weight is the same number. Only the Davis–Kahan weight makes the
eigenvector score mean "how much worse than the conditioning".

**Blessing the envelope baseline first** is not possible while 112 domination failures stand
(`0006`), and the committed rows are per id, so they need no bless. They are also strictly stronger
than no-regress, which only checks one direction.

## Consequences

- `PHASE2.md` §0.0's corpus-id row and its `Dual` row, and `PHASE3.md` §0.0's `renormalize` row,
  become `Done`, `svd3` excepted.
- The corpus grows by an estimated 6–7 MB, to about 44 MB of `PHASE1.md` §4.4's 50 MB.
- Any change to the arithmetic of these routines now re-records rows. That is the intent: it was
  invisible before.
- The envelope gains paired strata for these routines, and with them new domination failures where
  nalgebra wins.
- `0023` open question 3 and `0031` open question 1 are answered here. `0015` PH.1 and PH.4 are
  answered for these routines. `0044` item 5's `quat_*` family is now one id, before `0.0.1`:
  `norm` and `dot` left with `0048`.

## Implementation plan

Stacked PRs, each verified by `just lint` and `just test`, plus:

1. This record and the specification edits of decision 7, except `PHASE1.md` §4.3–§4.4, which
   `coverage` parses and which therefore land with step 4 — verified by `just lint`.
2. The generator: ids, strata and cross-checks; `just corpus` — verified by `just corpus-check` and
   the generator's unit tests.
3. `runners/nalgebra` and `just oracle-nalgebra` — verified by that recipe.
4. xtask: the subject at both precisions with its `HostStd` mirror, metrics 2–4, coverage, the
   self-test half, and the committed rows of decision 5, with *Measured* below — verified by
   `cargo xtask conformance --self-test`, the sanity tests and
   `the_unverified_routines_reproduce_their_committed_rows` on both CI architectures.

## Open questions

None.

## Measured

**Protocol.** `cargo xtask conformance --subject helicoid --fn <id> --precision f64|f32` over the
committed corpus (`just corpus`, 2 min 37 s on 8 cores, +4.9 MB to 42.8 MB); `just oracle-nalgebra`
for the oracle's rows. Every figure is a per-stratum maximum in `u` of the stratum's precision, and
every one of them is a committed row of decision 5 (107 binary64, 91 binary32), reproduced bit for
bit in the dev and release profiles.

| Id | binary64 strata, max `u` | binary32 strata, max `u` | Dominates nalgebra on |
|---|---|---|---|
| `solve_cubic` | 0 … 9.0e15 | 42.5 … 1.7e7 | 8 of 11 |
| `eig3` | 2.26 … 1.4e16 | 2.38 … 2.1e7 | 4 of 31 |
| `chol_n3`, `chol_n6` | 0 … 2.6e6 | 0 … 1.0e3 | 9 of 12 |
| `chol_solve_n3`, `chol_solve_n6` | 1.88 … 3.1e11 | 1.93 … 5.1e5 | 6 of 10 |
| `quat_renormalize` | 1.23 … 2.05 | 0.98 … 1.75 | 2 of 5 |
| `real_sqrt`, `real_cbrt`, `real_acos`, `real_sin_cos`, `real_atan2`, `real_div` | ≤ 3.47 | ≤ 3.65 | (no oracle) |

What the instrument now sees, each a documented limit made a number:

- **`eig3` is the closed form `0023` describes.** At a pair gap of `10^-8` and below, at the top of
  the spectrum, and at `eig:rank1`, a column reads `1e15`–`1.4e16 u` on the gap-weighted scale: the
  rustdoc's "near a double eigenvalue no column is reliable", which a column weighted by an `O(1)`
  gap makes a total loss. At the bottom the same gaps read `5e7`–`1e8`, about `√u`. nalgebra's
  `SymmetricEigen` reads 5–13 on every one of these strata. This is the measurement `0023` open
  question 2 was waiting for.
- **`solve_cubic`'s `cubic:one-real-p-small` reads `1.8e10 u`** against nalgebra's 3.0: the
  rustdoc's "`x³ + p x − 1` off by up to `4e-6` at `p` near `1e-5`".
- **`solve_cubic`'s `cubic:coeff-scale-down` reads exactly `1/u`, no root at all**, against
  nalgebra's `9.3e4`. Every coefficient times `10^-100` is the same cubic, and the leading
  coefficient's floor `max(scale, 1) tol` is absolute below scale 1, so the whole stratum is "not a
  cubic". It is the rustdoc's domain, stated, and `0031`'s L2 scaling question, now measured.
- **`cubic:double` reads `1.5e8 u`, `cubic:near-double-1e-8` `9.1e7`**: `√u` of the scale, as
  documented; `cubic:triple`, planted exactly, reads 0.
- **`chol_solve` follows `κ(A)`**: `3.1e11 u` at `chol:cond-1e-12`, a few `u` at `chol:spd` and
  `chol:diag-scale`. Its losses to nalgebra, and `chol`'s, are within a factor of 2.6.
- **`quat_renormalize` loses three strata to `renormalize_fast` by at most 0.6 `u`** (1.56 against
  0.98): nalgebra sums `‖q‖²` as `(x² + z²) + (y² + w²)`.
- **`real_*` is at most 3.65 `u`** at either precision, value and derivative together.

The self-test half fires each planted defect at `1.1e7`–`1.1e16 u` against clean maxima of 0–37.
The envelope gains 69 paired strata, 40 of them domination failures, and 29 that `helicoid` wins.

## Further work

1. `svd3`, its id, its strata (`0015` PH.2) and its oracle (`SVD`).
2. Kopp's hybrid for `eig3` (`0023` open question 2), informed by decision 4's domination failures.
3. An independent derivative oracle for `real_*`, such as `num-dual`.
4. Bench rows for `eig3` and `solve_cubic` (`0053` *Further work* 2).
5. The oracle's binary32 rows, once the envelope's binary32 half exists.
6. `solve_cubic` at scaled coefficients: a relative floor would answer `cubic:coeff-scale-down`
   (`0031` L2), and `one-real-p-small`'s cancellation is `0031`'s other arm; both are now rows that
   a change must re-record.
