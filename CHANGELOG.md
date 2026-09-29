# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [Unreleased]

### Added

- `docs/decisions/0020`: `Dual::sqrt` at zero keeps `d / (2 sqrt v)` (the NaN is the report and the
  `PHASE1.md` §10 row needs it); the guarded norm is a doctest, a zero-safe norm is a later record.
- `docs/decisions/0021`: the per-entry finite guard of `chol` stays; four bit-identical variants were
  measured, none clears the stated 15% bar, and `L` stays finite for every input. Documentation only;
  no code change.
- Tests pinning today's `Dual` derivative at a zero `sqrt` argument: `Vector<Dual>::norm` of the
  zero vector, and `chol` on a zero or a negative pivot.
- Workspace skeleton: `helicoid-linalg` and `helicoid` (empty, `no_std`, `forbid(unsafe_code)`),
  `xtask` stub, lints, `justfile`, CI, `deny.toml`. No group code (`docs/PHASE1.md`).
- `docs/maths/`: non-normative derivations. Notation, labelling and verification convention, and the
  Lie-group foundations (`Ad`, `ad`, the Jacobians of `Exp`, every row of `NUMERICS.md` §2.3, a sign
  audit with exact gaps, the conditioning of `J` and `Ad`). Documentation only; no code, no formula
  change.
- `docs/maths/so3.md`: SO(3) derivations (quaternion action and double cover, `Exp`, `Log` and the
  conditioning of the `acos` form it avoids, `J` and `J⁻¹` solved in the algebra of `{I, W, W²}`,
  `Ad`/`ad`, the action Jacobians, Shepperd's `from_matrix` and the nearest rotation, the Newton
  renormalization and its bounds). Documentation only; no code, no formula change. Five open items
  for the normative documents are recorded in the maths index.
- `docs/maths/coefficients.md`: the coefficient catalogue of `NUMERICS.md` §4 (general terms and radii
  of convergence, the identities between `k, a, b, c, d, e, r` and the `cos(θ/2)` of §3.1, rounding
  error of the exact arms and of Horner, truncation bounds, the crossing that predicts switch-point
  magnitudes, derivatives through `Dual`, the safe argument and the underflow, `π` and `w → 0`
  cases, the reference precision budget, the cost of one shared switch per call-site group).
  Documentation only; no code, no formula change. Nine open items for the normative documents are
  recorded in the maths index.
- `docs/maths/se3.md`: SE_N(3) derivations (the matrix group and the `x_1 = v`, `x_2 = p` convention of
  SE₂(3), `Exp`/`Log` from the hat matrix's series, `Ad`/`ad`, Barfoot's `Q` block derived from
  `J_l = Σ adⁿ/(n+1)!` by the block-triangular powers of `ad` and `W³ = −θ²W` and identified with `b, d, e`,
  `J_r` and two block relations, the dual-matrix algebra with proofs of closure, product (27 + 54N
  multiplications), inverse and `apply_transpose`, `J⁻¹` without a second closed form, the SE(3) action
  Jacobians, the rotation-first/translation-first permutation of tangents, Jacobians and covariances,
  what the exact arms of `b, d, e` cost `Exp` and `Q` (and the error of `Q` at a switch), and the rounding of
  `Log`'s translation part and of the blocks of `J⁻¹`). Documentation only; no code, no formula
  change. One open item for the normative documents is extended and one added in the maths index.
- `docs/maths/so2-se2.md`: SO(2) and SE(2) derivations (unit-complex `SO2`, `atan2` sensitivity; `Exp` with
  `V = αI + βK = s(θ)R(θ/2)`, `Log` and the three forms of `V⁻¹`, the two preimages at `|θ| = π`; `Ad`,
  `ad` and `ad³ = −θ²ad`; the rotation-first `J_r`, `J_l` and their inverses derived as `I ± a·ad + b·ad²`
  and `I ∓ ½ad + c·ad²`, with the dense entries written out and marked **Proposed** for `NUMERICS.md` §6
  and §2.4; the reduction of `Q` on the plane and agreement with Solà et al. App. C through the
  translation-first permutation; the scalar functions with series and singularities; rounding, the
  conditioning of `J` and a sign audit). Documentation only; no code, no formula change. Four open items
  for the normative documents are added in the maths index.
- `docs/maths/error-analysis.md`: the rounding model (`u = 2⁻⁵³`, `2⁻²⁴`; `libm` is not correctly rounded, measured),
  forward and backward error and the metric of `NUMERICS.md` §11 (chord, floors, what a norm-wise metric cannot
  see, `𝓑` from the stored reference and where it holds, why Jacobians are judged by forward error); the
  conditioning of `Log` (norm-wise `1/sin(θ/2)`, componentwise `≤ √5`, the cut), what backward error adds near π
  (branch invariance), `from_matrix` (uniform `1/(2√2)`, the backward-error floor `d_M`) and `J⁻¹` of SO(3) near
  `2π`; why bars are max and p99 and never a mean, non-finite counts, domination, exact no-regress and what D16
  buys (bit reproducibility across x86_64 and wasm32 with `libm`'s default features off, measured; not accuracy); forward-mode `Dual` (the
  rules, the value-path theorem, what the derivative of a branch or a series arm is, second order by nesting, the
  `sqrt`-at-0 hazard and the safe argument); the mpmath corpus (exact inputs, 30 digits, precision budget and the
  `dps = 150` recheck, log-uniform decades, uniform axes on S² by Archimedes, splitmix64 and per-stratum streams).
  Documentation only; no code, no formula change. Six open items for the normative documents are added in the maths index and one extended, and `coefficients.md`'s remark on the `libm` crate is corrected.
- `cargo xtask lint`, run by `just lint`: no line citations in markdown or comments, no draft
  decision record cited as settled (status-table rows, Rust comments, amendment banners; not
  prose elsewhere), and `@generated` files registered with their owning task and header. Owed:
  twin table (`docs/PHASE1.md`).
- `cargo xtask lint` reads `cargo metadata --locked` (all features): the normal-dependency
  closure of `helicoid-linalg` and `helicoid` must equal the `0007` set (`mint` only as an
  optional direct dependency, never enabled by `default`), reporting the path that brings in an
  offender; no crate but `xtask` may enable `__sweep`. `xtask` gains `serde` and `serde_json`;
  `deny.toml` allows `Unicode-3.0` for `unicode-ident` and bans the `0007` list, which
  `just audit` checks for every crate but `xtask`.
- `helicoid-linalg`: the scalar model (`docs/PHASE2.md` §2): `Mask`, `Real`, `Blend`, `Precision`;
  `f64`/`f32` as `Real` with every transcendental through `libm`; `Blend` for scalars, tuples up
  to arity 8 and arrays. `Real` has no `PartialOrd`/`PartialEq` (`compile_fail` doctests). New
  public API, nothing breaks.
- `helicoid-linalg`: forward-mode `Dual<S, N>` (`docs/PHASE2.md` §3), a `Real` and a `Blend<S>`
  with every derivative rule of §3; its value path is bitwise the plain `S` evaluation, NaN sign
  and payload of arithmetic outputs excepted. `copysign` takes `sgn(s)` from the sign bit. The
  derivative of `sqrt` at 0 is infinite or NaN, that of `atan2` is inexact once `x^2 + y^2` leaves
  the normal range, and the quotient's is NaN once `a / b` overflows (`NUMERICS.md` §12). New
  public API, nothing breaks.
- `helicoid-linalg`: `Vector<S, N>`, `Point<S, N>` and column-major `Matrix<S, R, C>` with the
  aliases `Vec2`/`Vec3`/`Mat2`/`Mat3`/`Point2`/`Point3` (`docs/PHASE2.md` §4): sums, `scale`,
  `dot`, `cross`, `norm`, generic matrix product, `transpose`, `identity`, `from_cols`,
  `from_rows`, `col`, `row`, `get`, `set`; `Point - Point` is a `Vector`, `Point + Vector` a
  `Point`. `hat`/`vee` and `Mat3::inverse_adj` (adjugate over determinant, plus the determinant).
  `Blend` for all three. No `PartialEq`. Every reduction sums left to right. New public API,
  nothing breaks. `proptest` is a dev-dependency only.
- `helicoid-linalg`: `chol` (lower-triangular `L` and a positive-definiteness `S::Mask`, set iff
  every computed pivot is strictly positive and every computed entry finite), `solve_lower` and
  `solve_upper` (`docs/PHASE2.md` §4; maths in `docs/NUMERICS.md` §15). Branch-free: `sqrt` gets a
  safe argument, a failed pivot gives `L_jj = 1` and a zero column, an overflowing entry is stored
  as `+0` and clears the mask, so `L` is finite for every input. The solves `debug_assert!` a
  nonzero diagonal. New public API, nothing breaks.
- `helicoid-linalg`: `Strided` and `StridedMut`, read-only and writable strided views over caller
  memory (`docs/PHASE2.md` §5): `col_major`, `row_major`, `with_strides` (faer and Ceres layouts),
  `block`, `get`, `set`. An out-of-bounds access panics in release, the index never wraps, and a
  write never leaves its view (D11's one panic class, widened from writes to `get`, `set` and
  `block` in `PROJECT.md` D11 and `NUMERICS.md` §12). New public API, nothing breaks.
- `helicoid-linalg`: optional feature `mint` (off by default; `mint` >= 0.5.7 joins the closure
  only with it): `From`/`Into` between `Vector<S, N>` (`N` = 2..=4), `Point<S, N>` (`N` = 2, 3; `mint`
  has no 4-D point) and `Matrix<S, N, N>` (`N` = 2..=4, `mint::ColumnMatrixN`), for `f32` and
  `f64` (`docs/PHASE2.md` §7). Scalars are moved, never computed on. `just lint` and `just test` cover the feature; `just msrv`, `just no-std`
  and `just wasm` build it. Bitwise round trips pin `-0`, infinities, subnormals and NaN payloads. New public API, nothing breaks.
- `docs/decisions/0018` (implemented): `libm`'s `arch` feature is bit-identical for exactly-rounded operations
  (dispatch table read for x86_64, aarch64, wasm32, thumbv7em; 95 digests over the `libm` entry points
  equal on x86_64, wasm32 (node WASI) and aarch64 (`qemu-user`); NaN sign and payload stay the target's,
  outside D16).
- Reference generator skeleton in `conformance/generate/` (uv, pinned CPython 3.12 and mpmath;
  splitmix64, per-stratum streams, `MANIFEST.json`, 150-digit recheck) and the first corpus file,
  `coeff_k`, all scalar-θ strata. New recipes `just corpus`, `corpus-check`, `corpus-test` and a
  `corpus-check` CI job. No library code changes; breaks nothing.
- Corpus files `coeff_a`…`coeff_e` and `coeff_r` (the θ strata, plus `q:w0` at three norms for `r`)
  and `coeff_series`, the exact 16-term rational Taylor series of every `NUMERICS.md` §4
  coefficient (manifest `kind: "series"`, `verified` in place of `rechecked`; `PHASE1.md` §4.3).
  Every record is cross-checked at generation against the series and a second formulation;
  cancelling definitions carry guard digits, so `theta:subnormal` is evaluated correctly. Manifest
  entries gain `kind`; `coeff_k` is unchanged. No library code changes; breaks nothing.
- Corpus files for SO(3): `so3_exp`, `so3_log`, `so3_act`, `so3_from_matrix`, `so3_jr`, `so3_jl`,
  `so3_jr_inv`, `so3_jl_inv` (22 642 records, 8.9 MB), computed from the definitions (quaternion and
  hat-matrix series, Newton's method on the exp series for `Log`, the polar factor for
  `from_matrix`) and each record checked to 100 digits by an independent property (`mp.expm`, the
  sandwich, polar uniqueness, Jacobian identities). Vector ids give each θ an axis; new strata `q:w0` and `q:nonunit`
  for the quaternion ids; every `so3_log` stratum but `q:w0` holds each quaternion and its
  negative. Matrices are column-major with a sibling `shape`. Generation is parallel by stratum
  (`--jobs`), byte-identical to a serial run. `docs/PHASE1.md` §2 and §4.3 state the `Log` and
  `from_matrix` references. No library code changes; breaks nothing.
- Corpus files for SE_N(3), N = 1, 2, 3: `sen3_exp`, `sen3_log`, `sen3_ad`, `sen3_jr`, `sen3_jl`,
  `sen3_jr_inv`, `sen3_jl_inv` (21 files, 7 578 records, 18 MB): `mp.expm` of the hat matrix,
  the inverse of it on θ ∈ [0, π], the conjugation definition of `Ad`, the defining series of
  `ad_τ` as dense rotation-first matrices, and `mp.inverse` of those (the dual matrix's zero
  blocks exactly 0, every other entry the LU's), each record checked to 100 digits (`mp.expm`,
  `J_l = Ad_Exp(τ) J_r`, `J J⁻¹ = I`) with the dual-matrix structure asserted on the data; no
  block form of `NUMERICS.md` §5 is evaluated. New strata: the 25 cells `rho:<scale>/theta=<θ>`
  and the `theta:*` strata at unit translation scale, 6 records each, and `q:w0` and `q:nonunit`
  for `sen3_log` and `sen3_ad`. `docs/PHASE1.md` §4.3 and §4.4 state the `sen3_log` reference
  (`mp.logm` is wrong near π), the inverse rule and the sample count. `just corpus` is 16
  CPU-minutes (3 minutes on eight cores); the `corpus-check` CI job gets a 60-minute timeout. No
  library code changes; breaks nothing.
- Corpus files for SO(2) and SE(2): `so2_exp`, `so2_log` (3 419 records each) and `se2_exp`,
  `se2_log`, `se2_ad`, `se2_jr`, `se2_jl`, `se2_jr_inv`, `se2_jl_inv` (416 each; 2.3 MB): `mp.expm`
  of the hat matrix, its inverse on θ ∈ (−π, π] by Newton's method on the series, the conjugation
  definition of `Ad`, the defining series of `ad_τ` as dense 3×3 rotation-first matrices and
  `mp.inverse` of those, each record checked to 100 digits of the size of its terms (`mp.expm`,
  the complex series, `J_l = Ad_Exp(τ) J_r`, `J J⁻¹ = I`) with the first row (1, 0, 0) and the
  rotation block asserted on the data. No closed form of `NUMERICS.md` §6 is evaluated. The `so2_*`
  ids see the θ strata with each θ in both signs, the `se2_*` ids SE_N(3)'s strata (no
  `theta:dense`) at 8 records each with the sign of θ alternating; z = (−1, ±0), θ = π, and a
  non-unit z are not sampled. `docs/PHASE1.md` §4.3 and §4.4 state the readings. `just corpus` is
  17 CPU-minutes. No library code changes; breaks nothing.
- `docs/decisions/0015-specification-gaps-found-while-building-the-instrument.md` (draft): the 29 gaps
  that the corpus generator, `helicoid-linalg`, `cargo xtask lint` and the maths pages found in
  `PHASE1.md`, `PHASE2.md`, `PHASE3.md` and `NUMERICS.md`, each with the passage, the evidence, the
  options, a recommended resolution and the work it blocks; the decisions in blocking order, then the
  readings and edits of the implementing PRs that await ratification (including the widening of D11 to
  strided reads). It decides nothing and edits no spec; `docs/maths/index.md` gains a pointer to it.
  Documentation only; breaks nothing.
- `helicoid`: the traits `Tangent`, `Jac` and `LieGroup`, `Side` with `Left` and `Right`
  (`docs/PHASE3.md` §2), the conventions of `0002` in the crate docs, and the trivial group
  `Rn<S, N>` with `RnTangent` and `RnJac` (§7): addition is `Mul`, there is no `Add` or `Sub`
  (`compile_fail` doctests), and every row of `NUMERICS.md` §2.3 is written for both sides.
  `Jac::sandwich` and the `DOF` tie are `const` assertions. Generic law checks (test-only
  `laws`) run for `Rn` and for a test-only non-abelian group (Heisenberg) under `f64`, `f32` and
  `Dual<f64, 3>` with recorded bounds. `docs/PHASE3.md` §7 and `docs/API.md` §3 name `RnTangent`
  and `RnJac`. New public API, nothing breaks. `proptest` becomes a dev-dependency of `helicoid`.
- `helicoid`: `Quat<S>` (`docs/PHASE3.md` §4): Hamilton product as `Mul`, `conjugate`, `norm_sq`,
  `to_matrix`, `from_wxyz_unchecked` (debug-asserts `|‖q‖² - 1| <= 2^-40` for `f64`, `2^-16` for
  `f32`), `from_wxyz_normalized` and `renormalize` (first-order Newton step, no domain
  asserted since `NUMERICS.md` states none; accuracy domain documented), and the converters
  `from_xyzw`, `to_xyzw`, `from_jpl`. New public API, nothing breaks; `SO3` is not here yet.
- `helicoid`: `SEn3Tangent<S, N>` with the `Twist` alias and `Twist::{omega, v,
  from_translation_first, to_translation_first}`, and `SEn3Jac<S, N>` with the dual-matrix algebra
  of `docs/PHASE3.md` §5 as a `Jac` (`mul` at `27 + 54N` multiplications, `inverse`, `neg`,
  `apply`, `apply_transpose`, `write_dense`, `sandwich`, `identity`). `helicoid::reference`
  (public) opens with the dense twins `sen3jac_mul` and `sen3jac_inverse` (Gauss-Jordan), writing
  into caller memory, and the proptests `sen3jac_mul_matches_reference` and
  `sen3jac_inverse_matches_reference` (`N = 1, 2, 3`; `f64`, `f32`, `Dual<f64, 3>`; the inverse's
  tolerance scales with the conditioning of the dense matrix). New public API, nothing breaks; the
  `SEn3` group is not here yet.
- `helicoid`: `Product<A, B>` and `ProductJac<JA, JB>` (`docs/PHASE3.md` §7): the product group
  with every `LieGroup` method the factors' method componentwise, tangent `(A::Tangent, B::Tangent)`
  (a `Tangent` impl for a pair, dense order `A` then `B`), `DOF = A::DOF + B::DOF`, and the
  block-diagonal Jacobian as a `Jac` (`sandwich` from the factors' `apply`). Nests. Tested by the
  generic laws on `Product<Rn<3>, Rn<2>>`, `Product<Rn<2>, Heis>`, `Product<Heis, Rn<2>>` and a
  nested product with two Heisenberg factors (test-only, the one non-abelian factor until SO(3)).
  New public API, nothing breaks. `Tangent::dot` documents that a composite tangent sums per factor.
- `just test` also runs the workspace tests in the release profile, where `debug_assert!` is
  compiled out: the tests of the documented out-of-domain behaviour are `cfg(not(debug_assertions))`
  and the dev profile skipped them. Nothing breaks.
- `cargo xtask conformance [--subject NAME] [--fn ID]` (`just conformance`), the harness core of
  `docs/PHASE1.md` §5: the corpus reader, the `Subject` trait with an empty in-process registry, the
  forward error of `NUMERICS.md` §11 in units of `u = 2⁻⁵³` formed exactly in integers and rounded
  once (the 30-digit reference is never parsed to `f64`), per `(fn, stratum, precision, subject)`
  `n`, `max_u`, `p99_u` (nearest rank), `argmax_id` and `nonfinite`, written to
  `conformance/results/<subject>.csv` and printed by `max_u` descending. A non-finite output, an
  error that overflows binary64 and a run that scored nothing (the registry is empty, so today's
  run) fail; `so3_from_matrix` is checked for shape and finiteness only. The floor of a tangent or a
  coefficient (`2⁻¹⁰²²`), the SE_N(3) tangent's sign, the maximum over a record's fields, one norm
  over a whole tangent and overflow as non-finite are readings recorded in `xtask/src/conformance/`
  and proposed by draft record 0014. Backward error, `f32` (refused), `--self-test` and the
  `helicoid` and seeded subjects are not implemented. `xtask` gains `num-bigint` and
  `helicoid-linalg`; no library code changes; breaks nothing.
- `cargo xtask conformance --self-test` (`just conformance --self-test`, a CI job) and the seeded
  subjects of `docs/PHASE1.md` §10 (`xtask/src/seeded`): the coefficient kernels `k, a, b, c, d, e`
  generic over `Real` (exact arm at the safe argument, one `S::branch`; series read at run time from
  `coeff_series.jsonl` and rounded once, in integers, at the precision they run at, `f64` or `f32`),
  evaluated over the `coeff_*` ids as value and `d/dz` of one `Dual<f64, 1>` at a
  `Candidate { terms, switch_z }`; `seeded:correct` is the `tf_tree` D12 candidate (for all six
  coefficients, though `NUMERICS.md` §4 lists it for `a`, `b`, `c`), and `b` by its definition and `k`
  with an unsafe `sqrt` of `θ²` under `Dual` are planted. The self-test fails unless the correct kernel
  fires no mechanism and each defect fires its own (`b`: the value's error curve over
  `theta:1e-8`…`theta:1e-2` fits `θ^-p`, `p` in [1.8, 2.2], 1.937 measured; `k`: `nonfinite > 0` in
  `theta:exact0`), and prints the correct kernel's errors and every mechanism's reading for every
  subject. Unit tests pin the correct kernel's `max_u` per coefficient and field, and check its safe
  argument on a mask that evaluates both arms. A plain `just conformance` now scores `seeded:correct`
  and skips planted subjects unless named. `xtask` gains `libm` (the fit's logarithm, the same bits on
  every host). No library code changes; breaks nothing.
- `cargo xtask thresholds [--check]` (`just thresholds`, `thresholds-check`, a CI job), the sweep of
  `docs/PHASE1.md` §6 over the seeded kernels at `f64` for `k, a, b, c, d, e`: 1 to 8 series terms
  times 1025 switch points (64 per decade of `z = θ²`, `θ` from `1e-8` to `1`, each the correctly
  rounded `10^(i/64)`), each arm's error formed once per record over every `theta:*` record (value
  and `d/dz` through `Dual<f64, 1>`, exact, in `u`), the maximum minimised, ties to fewer terms
  then the larger switch, the `tf_tree` D12 prior scored beside the choice. Writes
  `conformance/sweeps/thresholds.csv` (committed, byte-identical, compared with a fresh run by a
  test). The rows are of this corpus and grid, not general optima: `terms` is the cap of 8 for every row but `e`;
  `k` and `d` sit at the top of the grid (`θ = 1`), `b` one step below it (unique, 1.3% under the
  top), `a` and `c` inside it; `e`'s objective is one record at `θ = 1` exactly, which `z < switch`
  never puts on the series arm (a top of `nextUp(1)` gives it 8 terms and 3.3 times less). Columns
  name the record that attains each objective and the objective one grid step either side. The
  grid, top, prior and tie readings are `0014` (draft) questions 8 to 11. Not swept: `r`, `f32`,
  SE(2)'s `α`, `β`; `generated.rs` waits for Phase 3. No new dependency, no library code changes;
  breaks nothing.

### Changed

- `libm` is built with its `arch` feature (`0018`): `sqrt`, `sqrtf` (x86_64, aarch64), `fma`, `fmaf`
  and `rint`, `rintf` (aarch64; `fma` by `cpuid` on x86_64) run the target's instruction. All are exactly
  rounded, so every output that is not NaN is unchanged on every target; the sign and payload of a NaN
  from `sqrt` and `fma` are now the target's, as those of `+ - * /` always were (so `copysign` of the
  NaN of `sqrt(-1)` takes the target's sign; that input is outside `Real::sqrt`'s domain).
  `Real::sqrt` is about 4x faster on x86_64 (measured; `docs/decisions/0018`, Consequences). wasm32 and
  thumbv7em are unchanged: `libm` routes nothing there on stable. New test: `Real::sqrt` for `f64` and
  `f32` against an integer square root, bit for bit, over `2 x 10^6` cases with a pinned digest; checked
  by hand on x86_64, aarch64 (`qemu-user`) and wasm32 (node WASI), run by CI on x86_64 and aarch64.
- `helicoid-linalg`: `chol_solve(&l, b)`, `A x = b` from a Cholesky factor: `solve_lower`, then a
  back substitution against `l^T` read by column, with no transposed copy. Bit-identical to
  `solve_upper(&l.transpose(), solve_lower(&l, b))` (NaN sign and payload aside), same domain and
  release behaviour; `docs/decisions/0019`, `docs/NUMERICS.md` §15.6. New public API, nothing breaks.
