# 0103: The inverse action has its own Jacobians

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #195, #196, #197

## Context

$x = X^{-1}p$ maps a point from a parent frame into a body frame. It is the first step of every
reprojection residual and of every projection covariance with pose terms. Its consumers include
locus-sensors' `projection_covariance` (`PHASE-C1` §3.2), locus-calib's backend and the VIO.

helicoid 0.0.6 ships the forward action and pieces of the inverse, but not the inverse action
itself:
- `X * p` and `act_jacobians::<Sd>` (`NUMERICS.md` §2.4, §5.5), on SO(3) and SE(3);
- `inverse` and `inverse_jacobian::<Sd>` (§2.3).

So a consumer either composes `X.inverse() * p` with
`act_jacobians(X⁻¹, p) · inverse_jacobian(X)`, or derives the Jacobian by hand, which D1 forbids.

**The composition is exact in algebra, not in floating point.**
- **The value.** With $X^{-1} = (R^\top, -R^\top t)$, the composed value $R^\top p + (-R^\top t)$
  rounds two vectors of size $\lVert t\rVert$ and then cancels them down to $\lVert x\rVert$.
- **The right-side Jacobian.** Its rotation block is
  $(-R^\top[p]_\times)(-R) + R^\top(-[t]_\times R) = R^\top[p]_\times R - R^\top[t]_\times R$.
  That equals $[x]_\times$ only through the same cancellation.
- **The left-side Jacobian.** Its rotation block is
  $[x]_\times R^\top + [R^\top t]_\times R^\top = [R^\top p]_\times R^\top$, which cancels where
  $\lVert p\rVert \ll \lVert t\rVert$: a point near the parent's origin, seen from far away.
- **When it matters.** The value and the right side lose $u\lVert t\rVert/\lVert x\rVert$, the
  left side $u\lVert t\rVert/\lVert p\rVert$. The first is the common case: map or world frames
  far from the origin (UTM, a long trajectory), and a camera close to the points it sees.

A scratch Python harness measured it (binary64, sums left to right, no FMA; 2000 draws per row,
seed 20261021; $q$ uniform on $S^3$, $\lVert x\rVert$ uniform in $[1, 10]$, $t$ of norm
ratio·$\lVert x\rVert$ in a uniform direction, $p = Rx + t$ rounded as the input; mpmath at 60
digits from the rounded $(q, t, p)$; §11's forward error, Frobenius, floor 1, max over the draws):

| $\lVert t\rVert/\lVert x\rVert$ | $x$: direct / composed | right $J$: direct / composed | left $J$: direct / composed |
|---|---|---|---|
| 1 | 3.56 / 4.76 | 3.50 / 14.5 | 2.78 / 8.52 |
| $10^2$ | 3.30 / 325 | 2.99 / 287 | 2.86 / 9.98 |
| $10^4$ | 3.32 / $3.28\times10^4$ | 3.16 / $2.91\times10^4$ | 3.04 / 12.0 |
| $10^6$ | 3.49 / $3.14\times10^6$ | 3.42 / $2.93\times10^6$ | 3.00 / 10.3 |

With $p$ drawn as $x$ was above ($\lVert p\rVert$ uniform in $[1, 10]$), 1000 draws per row:

| $\lVert t\rVert/\lVert p\rVert$ | $x$: direct / composed | right $J$: direct / composed | left $J$: direct / composed |
|---|---|---|---|
| 1 | 3.74 / 8.18 | 3.65 / 11.1 | 2.89 / 9.51 |
| $10^2$ | 3.09 / 3.31 | 3.09 / 12.4 | 2.63 / 201 |
| $10^4$ | 3.67 / 3.51 | 3.67 / 10.8 | 3.06 / $2.05\times10^4$ |
| $10^6$ | 3.87 / 3.22 | 3.87 / 10.8 | 2.71 / $1.93\times10^6$ |

Each composed side has a case it loses; the direct forms hold in both, below 4 `u`. The harness is
not committed: decision 5's planted subject reruns the same comparison over the corpus, which is
what makes these tables citable (`0040` item 4).

**The formulas were checked against Exp.** On five random poses, `mp.diff` through §1's
$\mathrm{Exp}$ of $(X\,\mathrm{Exp}\,\tau)^{-1}p$ and $(\mathrm{Exp}(\tau)\,X)^{-1}p$, at 40
digits, matches decision 2's Jacobians to $9\times10^{-41}$.

## Decision

1. **The surface.** `SO3` and `SE3` (`SEn3<S, 1>`) gain, next to `act_jacobians`:

   ```rust
   impl<S: Real> SO3<S> {
       /// `R⁻¹ v = Rᵀ v`: `self.inverse().act(v)`, to the bit.
       pub fn inv_act(&self, v: Vec3<S>) -> Vec3<S>;
       /// `(∂(R⁻¹ v)/∂R, ∂(R⁻¹ v)/∂v)` in the convention of `Sd`.
       pub fn inv_act_jacobians<Sd: Side>(&self, v: Vec3<S>) -> (Mat3<S>, Mat3<S>);
   }
   impl<S: Real> SEn3<S, 1> {
       /// `X⁻¹ p = Rᵀ (p − t)`.
       pub fn inv_act(&self, p: Point3<S>) -> Point3<S>;
       /// `(∂(X⁻¹ p)/∂X, ∂(X⁻¹ p)/∂p)` in the convention of `Sd`.
       pub fn inv_act_jacobians<Sd: Side>(&self, p: Point3<S>) -> (Matrix<S, 3, 6>, Mat3<S>);
   }
   ```

   The pair mirrors `act` and `act_jacobians`. The 3×6 is dense under `0005` for
   `act_jacobians`' reason: no group's tangent is its row space.
2. **`NUMERICS.md` §2.4 gains**, with $x = X^{-1}p$:
   - $\partial(X^{-1}p)/\partial X$: right $[\,[x]_\times,\ -I\,]$; left $[\,R^\top[p]_\times,\ -R^\top\,]$;
   - $\partial(X^{-1}p)/\partial p = R^\top$;
   - SO(3), $x = R^\top v$: right $[x]_\times$, left $R^\top[v]_\times$, $\partial/\partial v = R^\top$.

   The derivations are one line each: right $(X\,\mathrm{Exp}\,\tau)^{-1}p = \mathrm{Exp}(-\tau)\,x$;
   left $(\mathrm{Exp}(\tau)\,X)^{-1}p = X^{-1}\,\mathrm{Exp}(-\tau)\,p$.
3. **`NUMERICS.md` §3.3 and §5.5 gain the value.** $R^\top v$ is §3.3's sandwich of the conjugate,
   which negates $u$ exactly. $X^{-1}p = R^\top(p - t)$ takes the subtraction first: it rounds
   once relative to its own result, $\lvert\mathrm{fl}(p - t) - (p - t)\rvert \le u\lvert p - t\rvert$
   componentwise (exactly, where $p$ and $t$ agree to within a factor of two), so the rotation sees
   $x$'s scale. It is never evaluated as $R^\top p - R^\top t$. It is `inv_mul`'s column
   arithmetic (`0048`): `SE3::from_rt(I, p).inv_mul(&X).translation()`, to the bit.
4. **Evaluation.**
   - **Right side.** The rotation block is `hat` of the value `inv_act` returns, the same
     sandwich, not $R^\top(p - t)$ through the matrix; $-I$ is written, not negated, so its
     off-diagonal zeros are `+0` (a negated identity carries `−0` there). So the block is
     `hat(X.inv_act(p))` to the bit, and a caller that asks for both evaluates the sandwich once:
     the two are one pure expression, which LLVM merges without fast-math (*Rationale*).
   - **Left side.** `mul_hat(&Rᵀ, p)` and $-R^\top$, with no $t$ to cancel.
   - **Both sides.** $R^\top$ is the transpose of one $R(q)$, a move; it is also $R(q^*)$ to the
     bit, since $R(q)$'s entries are symmetric in the sign of $u$. SE(3)'s rotation block is
     SO(3)'s `inv_act_jacobians` at $v = p - t$ (right) or $v = p$ (left).
   - **`Matrix::from_cols`, `scale` and `Neg` index** (`helicoid-linalg`): `array::from_fn` in
     place of the nested `map`s. The `map` over six columns was left as an out-of-line call and a
     `memcpy` in every `3 x 6` this record and `act_jacobians` build. The same moves and the same
     operations, so every output keeps its bits.
5. **Verification.**
   - **Corpus ids `so3_inv_act` and `se3_inv_act`.**
     - Inputs $q, x, p$ (`x` the translation, as every SE(3) id names it) and $q, p$. Outputs
       `Xinv_p` (`Rinv_p`), `Jr`, `Jl`, `Jp`, keyed shapes (`0101`).
     - `se3_inv_act` strata: `tx:1`, `tx:1e2`, `tx:1e4`, `tx:1e6` ($\lVert t\rVert/\lVert x\rVert$,
       $\lVert x\rVert$ uniform in $[1, 10]$, $p = Rx + t$ rounded), `tp:1` … `tp:1e6`
       ($\lVert t\rVert/\lVert p\rVert$, $\lVert p\rVert$ uniform in $[1, 10]$) and `generic`, at
       binary64 and binary32 (`@f32`, the same draws rounded once). `so3_inv_act` takes
       `so3_act`'s strata.
     - Reference: mpmath at 120 digits from the definitions, scored by §11's forward error with
       floor 1 on every field. Not `Floor::Scale` on the translation: dividing by $\lVert t\rVert$
       would excuse the composition the loss this record exists to measure.
     - Oracles at binary64 (`0077`): GTSAM's `Pose3::transformTo(p, H1, H2)` and
       `Rot3::unrotate(p, H1, H2)`, the direct form through a matrix; manif's
       `inverse().act(p, J, J)` chained with `inverse(J)`, the composition. Each oracle's left
       Jacobian is its right one times its own $\mathrm{Ad}_{X^{-1}}$, as `0074` and `0075` derive
       left Jacobians. manif also answers `se3_inv_act`'s `@f32` strata in `manif::SE3<float>`,
       the first id of its runner with any (`0077`'s `--precision f32`), so the binary32 rows
       have an oracle too; GTSAM is `double` throughout.
   - **The composition as a subject.** `helicoid:composed-inv-act`, planted as
     `helicoid:homogeneous-r` is, answers both ids by the composition; `cargo xtask compare` against
     `helicoid` reproduces *Context*'s tables over the corpus.
   - **`Dual` tests.** `inv_act_jacobians_differentiate_the_inverse_action` (in `so3_tests` and
     `sen3_tests`) differentiates `inv_act` through `rplus` and `lplus`, and in `p`.
   - **Reference twin (§14).** The composition $X^{-1} * p$ and
     $\mathrm{act\_jacobians}(X^{-1}, p)\cdot\mathrm{inverse\_jacobian}(X)$, by
     `so3_inv_act_matches_reference` and `se3_inv_act_matches_reference`. Its loss in *Context* is
     the twin's own; the proptest bounds the two's difference in $u\max(\lVert p\rVert,
     \lVert t\rVert, 1)$, `inv_mul`'s denominator (`0048` decision 4).

## Rationale

- **Not left to consumers.** Composition costs up to $u\lVert t\rVert/\lVert x\rVert$ on the value
  and the right side, and $u\lVert t\rVert/\lVert p\rVert$ on the left (*Context*). A hand
  derivation in each consumer is what D1 exists to prevent: three copies of $[\,[x]_\times, -I\,]$
  (locus-sensors, locus-calib, the VIO) with three sign conventions to get wrong.
- **Not a second spelling (`API.md` R3).** The inverse action has no operator, and the composition
  is not the same rounding: it is this pair's §14 twin, as `b.inverse() * a` is `inv_mul`'s.
- **SO(3) too** (the owner, 2026-10-21). Its value is the composition to the bit, so the pair wins
  no accuracy there on `inv_act`; it wins on the Jacobian, whose composed right block
  $R^\top[v]_\times R$ rounds two products that cancel on the diagonal (the twin reads up to
  23.6 `u` against the direct form), and it gives both groups that have an action one spelling of
  the side convention (`API.md` §6 question 8, `0049`).
- **Measured on the corpus.** Over the deep tier (64 records a stratum) the pair reads at most
  4.77 `u` at `f64` and 4.79 at `f32` on `se3_inv_act`, and its maximum is below GTSAM's direct
  form on every stratum (`tx:1e6` 3.66 against 4.47, `generic` 4.77 against 6.01); GTSAM's left
  Jacobian, its right one times `Ad_{X⁻¹}`, loses `tp:*` as *Context*'s composition does
  ($2.4\times10^6$ `u` at `tp:1e6`). The planted composition reads $4.0\times10^6$ on `tx:1e6` and
  $2.2\times10^6$ on `tp:1e6`, *Context*'s tables to their order. `so3_inv_act` reads at most
  5.03 `u` off `q:nonunit`, against GTSAM's 7.47 and manif's 11.6; on `q:nonunit` it carries
  `so3_act`'s drift (511 `u`), which neither oracle normalizes away (806 and 1469), so no row of
  `exceptions.toml` is owed. `just envelope --deep`: no domination failure on either id.
- **Not a chart (D12).** The inverse action is a group method, as the action is. It is not a
  retraction, and no chart reparametrizes it.
- **A pair, not a fused triple.** `act`/`act_jacobians` set the shape, and the measurement below
  says the fused form would buy nothing: the pair behind one call compiles to one sandwich (28
  multiplications at `f64`, where two sandwiches and `R(q)` would be about 51).
- **The right block takes the value's sandwich, not the matrix it already holds.** Alone,
  `inv_act_jacobians::<Right>` would be cheaper as `hat(Rᵀ(p − t))` through the `R(q)` it forms
  for $\partial/\partial p$ (the left side, which needs no value, is the measure: 9.2 against
  15.0 ns). But the caller this record is for asks for the value too, and would then pay a
  sandwich and a matrix product where it pays one sandwich, and the block would no longer be
  `hat(inv_act(p))` to the bit. The merge needs both calls inlined into the caller, which
  `#[inline]` on both gives; a caller that keeps them apart, or asks for the Jacobian alone, pays
  the second rotation, about 6 ns.
- **Cost** (`benches/groups.rs`, `generic-1`, one pinned core, `f64` / `f32`; *before* is the
  same code with `Matrix`'s `map`s, from `cargo xtask bench-gate --against`'s paired ratio, three
  replicates, every row inside its concurrent A/A floor but the ones below):

  | row | before | after |
  |---|---|---|
  | `se3/inv_act` | 12.4 / 13.1 ns | 12.4 / 13.1 ns |
  | `se3/inv_act_jacobians/right` | 31.6 / 31.6 ns | 15.0 / 16.3 ns |
  | `se3/inv_act_jacobians/left` | 27.4 / 24.6 ns | 9.2 / 9.3 ns |
  | `se3/inv_act_pair` (value and right Jacobians, one call) | 40.1 / 35.5 ns | 20.4 / 21.2 ns |
  | `se3/inv_act_composed` (the twin, `reference::se3_inv_act_jacobian`) | 128 / 124 ns | 90.5 / 90.7 ns |

  The pair takes 0.23× the reference composition's time at `f64`; a hand-fused dense product in
  place of the reference's took 56 ns. Without the same `map`s `Gaussian::propagate` takes 0.65×
  of its time, `to_left` 0.70×, `to_right` 0.75× and `mahalanobis_sq` 0.89×; no other row of
  `so3/*` and `se3/*` moved beyond its floor.

## Consequences

- `NUMERICS.md` §2.4, §3.3, §5.5, §12 (the carried row) and §14 gain decision 2's to 5's lines.
- `docs/API.md` lists the pair.
- `CHANGELOG.md` records it.
- **Two new corpus ids**, `so3_inv_act` and `se3_inv_act`, with `PHASE1.md` §4.3 and §4.4 rows,
  `coverage::required()` entries, envelope rows and the two oracles' answers. `so3_inv_act` on
  `q:nonunit` inherits `so3_act`'s carried drift (`0076`); no oracle wins there, so it owes no
  exception row.
- The manif runner answers in binary32 (`--precision f32`), for `se3_inv_act` alone.
- **locus-sensors' `PHASE-C1` §3.2** composes the pose Jacobian at the identity today. That
  composition is exact because it starts from $x$, not from $p$. It moves to
  `inv_act_jacobians::<Right>` when its pin reaches the release: same values, up to the sign of
  zeros.

## Implementation plan

1. `NUMERICS.md` §2.4, §3.3, §5.5, §12, §14; `inv_act` and `inv_act_jacobians` in `so3.rs` and
   `sen3.rs`; `Matrix::from_cols`; `docs/API.md`; the bench rows. Verified by both
   `inv_act_jacobians_differentiate_the_inverse_action`, the bit relations of decisions 3 and 4
   (`inv_act_is_the_inverse_s_action_and_the_right_block_its_hat`,
   `inv_act_is_inv_mul_s_column_and_its_jacobians_reuse_it`), and both `*_inv_act_matches_reference`.
2. Corpus ids `so3_inv_act` and `se3_inv_act` (generator, registry, cross-checks, harness answer,
   metric, coverage) and the planted `helicoid:composed-inv-act`. Verified by `just corpus-check`,
   `just conformance` and `cargo xtask compare`, whose rows reproduce *Context*'s tables.
3. GTSAM and manif answer both ids. Verified by `just oracle-gtsam`, `just oracle-manif` and
   `just envelope --deep`.

## Further work

- `inv_act_many(&mut [Point3])`, the per-slice form, as `act_many` is to `act`.
- Sim(3)'s inverse action, when a consumer needs it.
