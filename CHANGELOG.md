# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [Unreleased]

### Changed

- `solve_cubic` internals, every one bit-identical (`bits_are_omnisacs` and the planted-root suite
  are unmoved): `pi` is a per-precision literal rather than `atan2(+0, -1)`, which was a `libm` call
  for a compile-time constant that no workspace LTO could fold — `eig3` paid it per call too; `tol`
  takes its exponent as a const parameter and evaluates `pow2` in a `const` block, so the thresholds
  are literals and `pow2`'s stated domain is a compile error rather than a silently wrong float;
  `pow2`, `tol`, `max`, `pi` and `acos` gain `#[inline]`; the three-root arm hoists `theta / 3`,
  `2 pi / 3` and `2 r` out of its per-phase closure; `term_q` loses an `abs` it cannot need and
  `a.abs()` is bound once.
- The single-and-double-root arm's unreachability is now an executable claim: a `debug_assert!` that
  `repeated` implies `triple`, exercised by the random-bit-pattern and planted-root proptests. It is
  reachable only if the leading-coefficient floor is lowered, which `0031` (draft) L3 recommends, so
  the assert fires before the arm silently starts answering.
- `pi_and_acos_are_within_their_ulps` keeps the `atan2` twin (D6): the literal is asserted bit-equal
  to `0.0.atan2(-1.0)` at both precisions, so the evidence that justified the swap does not vanish
  into the swap.
- `dual_tests`' `ulps`/`ulps32` and their assertion wrappers come from one `macro_rules!` per width,
  so a change to the total-order key cannot reach one precision and miss the other.
- `eig3`'s three private helpers gain `#[inline]`, the convention the rest of the crate keeps; it
  also picks up `pi` as a literal, which it was paying as a `libm::atan2` per call.

### Added

- `docs/decisions/0022` (**ready**, retitled *`Real` owes `acos` and `cos`*): the private
  `atan2(sqrt((1 - x)(1 + x)), x)` costs 2.61x `libm::acos` (11.42 ns against 4.37) and is
  marginally *less* accurate everywhere measured, the near-`±1` region it was chosen for included
  (1.28 u against 0.92 as `x -> 1`), and it forfeits bit-identity with omnisac in the trigonometric
  arm. `solve_cubic` and `eig3` are both callers — `eig3`'s own PR amended this record to say so —
  so "public surface for one caller" does not hold. `Real::cos` saves a further 26% per call over a
  `sin_cos` whose sine is discarded. `pi` as a literal has landed; the two methods are step 1.
- `docs/decisions/0031` (draft): the four numerical limits the cubic port inherits from omnisac,
  each measured, with a recommended order. The one-real-root cancellation first — `x^3 + p x - 1` at
  `p = 1e-5` returns exactly `1.0` for a root of `0.999996666666666679`, a relative error of 3.33e-6,
  which the `u v = -p/3` pairing takes to 3.42e-17 and which also removes a `-inf` `Dual` derivative
  at a simple root; it dominates on the max (4.67e-10 to 2.63e-12 over 16 010 random one-root
  cubics) but is *not* uniformly better per row. Then the underflow that reports three valid slots
  for a one-root cubic, the leading floor that rejects `(x - 100)(x - 101)(x + 99)` at `f32`, the
  discriminant band measured against the wrong quantity, and the `1/a` reciprocal. Nothing decided.
- `docs/decisions/0028` (draft): `SEn3Jac`'s `diag`/`col` are `pub` in `PHASE3.md` §5 and private
  under `0025` decision 5, and `ProductJac` took the other reading; three options, their costs, and a
  recommendation (narrow, the reversible direction). Nothing decided, no code change.
- `docs/decisions/0029` (draft): `Product` has private fields and no constructor, so the tf2 pose
  `Product<SO3, Rn<3>>` is reachable only through an `Exp`/`Log` round trip that is inexact, costs an
  `atan2` and a `sin_cos`, and is worst-conditioned at a half turn; `from_parts`/`parts` proposed, to
  land with the SO(3) PR. Nothing decided, no code change.
- `just lint` runs `cargo clippy --workspace --all-targets --release -- -D warnings` as well as the
  dev pass: `debug_assertions` is off in release, so the `cfg(not(debug_assertions))` test bodies that
  only `just test`'s release run executes were never checked against the denied `unwrap_used`,
  `expect_used`, `panic`, `todo`, `unimplemented` and `dbg_macro`. Clean today, so this closes a gap
  rather than fixing a violation.

- `docs/decisions/0025` (ready): the trait layer as built. Rⁿ's Jacobian is the structured `RnJac`,
  correcting `PHASE3.md` §7's `Jac = Mat<N>`, which cannot implement the exact `Jac::inverse` of
  `0005`; `Side` is sealed to `Left`/`Right` while the side selector stays with the SO(3) PR;
  `Tangent::dot_acc` becomes the required operation so that `Product` can meet the bound-`0`
  dense-order law; `read_dense` poisons a missing entry; and a group exposes an `Add`-carrying field
  only where it is abelian.
- `docs/decisions/0020`: `Dual::sqrt` at zero keeps `d / (2 sqrt v)` (the NaN is the report and the
  `PHASE1.md` §10 row needs it); the guarded norm is a doctest, a zero-safe norm is a later record.
- `docs/decisions/0021`: the per-entry finite guard of `chol` stays; four bit-identical variants were
  measured, none clears the stated 15% bar, and `L` stays finite for every input. Documentation only;
  no code change.
- Tests pinning today's `Dual` derivative at a zero `sqrt` argument: `Vector<Dual>::norm` of the
  zero vector, and `chol` on a zero or a negative pivot.
- The `helicoid` conformance subject (`docs/PHASE1.md` §5, `xtask/src/shipped.rs`): the shipped
  groups of `helicoid::coeffs` through `__sweep`, over the eight `coeff_*` ids at `f64` and `f32`,
  `value` and `d_branch` on `Dual<S, 1>`, registered as a plain subject, so `just conformance`
  prints and writes its per-`(fn, stratum, precision)` rows beside `seeded:correct`'s. It equals
  `seeded:correct` bit for bit on every record (the two share every candidate) and scores finite
  everywhere; a four-terms-below-`z = 0.01` prior (D12 for `a`, `b`, `c`; not D12's for `r`, `k`,
  `d`, `e`, `cos θ/2`: 0014 (draft) questions 10, 29) has a smaller maximum on nine `cos θ/2`
  strata, pinned by a test and raised as 0014 question 30. No baseline, no envelope. Tooling
  only; nothing breaks.
- The coefficient kernel (`docs/PHASE3.md` §3, `0004`): `helicoid::coeffs`, `pub(crate)`, with the
  groups `exp_coeffs`, `jr_coeffs`, `jr_inv_coeff`, `q_coeffs` and `log_ratio`, each one `S::branch`
  on "every member is on its series arm", off which the group's exact arms run once, at the safe
  argument, and each member selects by its own mask; the exact arms in the operand order of the
  seeded kernels, the series by Horner, the constants chosen by `S::PRECISION`. `cargo xtask
  thresholds helicoid` sweeps the shipped arms through the hidden feature `__sweep` (`exact_*`,
  `series_*`; only `xtask` enables it) and writes `crates/helicoid/src/coeffs/generated.rs` (a
  `Switch` per coefficient per precision, beside the eight swept series they take their terms from;
  registered with `cargo xtask lint`) and `conformance/sweeps/thresholds.csv`, a fixed point: a
  second run writes the same bytes, and `just thresholds-check` covers both this target and the
  seeded one (`cargo xtask thresholds [--check] [seeded|helicoid]`). At `f64`: `k` 8 terms below `θ²
  = 1`, `a` 8 below `0.72`, `b` 8 below `0.96`, `c` 8 below `0.56`, `d` 8 below `1`, `e` 7 below
  `0.96`, `cos θ/2` 2 below `5.6e-15`, `r` 8 below `n²/w² = 1.1e-2`; at `f32`: `k, a, b` 5 terms and
  `c, d, e` 4, all below `θ² = 1`, `cos θ/2` 2 below `2.6e-6`, `r` 8 below `0.143`. They are the
  seeded sweep's rows, and the shipped arms are the seeded arms bit for bit on every record where
  the series arm is defined (`w > 0` for `r`). Each objective is a maximum over the corpus's
  records, not over the domain, and eight of the sixteen switches (`k`, `d` at `f64`; `k`, `a`…`e`
  at `f32`) are `θ² = 1`, the top of the sweep's grid, not a measured optimum. The sweep CSVs gain
  four columns, each arm's error at the two records that bracket the switch, whose sum
  `branch_continuity_*` compares the jump between the arms with (`docs/maths/coefficients.md` CO.12;
  a sample at two records, not a bound over the interval; not `NUMERICS.md` §4's wording); the
  seeded CSV and the digest line of `xtask/src/seeded/generated.rs` change with them, its constants
  do not. Groups take their members' switches each (0015 (draft) NU.5): no exact arm runs below a
  group's smallest switch, which for `Exp` is `cos θ/2`'s (`θ² < 5.6e-15` at `f64`), so `Exp` is on
  its exact arm from `θ = 7.5e-8` at one `sqrt` and one `sin_cos` (`NUMERICS.md` §3.1); unmeasured.
  `r` is on its provisional reading (0014 (draft) question 29); `se2_coeffs` and `gamma2_coeffs` are
  not here. The kernel is built for the tests and `__sweep` only until SO(3) uses it, so `just
  no-std` and `just wasm` build it with `--features __sweep`. No public API; nothing breaks.
- Threshold sweep per precision (`0016` item 3): `cargo xtask thresholds` now sweeps `k, a, b, c, d, e`,
  `cos θ/2` and `r` at `f64` on the plain strata and at `f32` on the `@f32` strata, with one objective,
  grid (each point correctly rounded at the precision it is swept at) and tie-break. At `f32` each of
  `k`…`e` sits at the top of the grid (`θ = 1`), its optimum at or above it, as `docs/maths` CO.10
  predicts; `cos θ/2` sits at the exact arm's floor at both precisions, its switch the tie-break's.
  `r` is swept on a reading pending the maintainer, not a spec: branch variable `s = n²/w²`, only
  records with `w > 0`, series arm iff `w > 0` and `s < switch`, grid over `s` from `1e-16` to `1`,
  its prior `(4 terms, s < 0.01)` and not D12's (0014 (draft) question 29). The seeded sweep's
  outputs are renamed `conformance/sweeps/thresholds-seeded.csv` (16 rows) and
  `xtask/src/seeded/generated.rs` (a `Switch` per id and precision, the `f32` literals correctly
  rounded from the exact rationals, a test in integers), leaving `thresholds.csv` to
  `helicoid::coeffs`. `seeded:correct` runs the `f32` constants (where it ran the D12 prior;
  `subject_version` `generated`) and now also answers `coeff_cos_half` and `coeff_r`; the
  self-test's pinned `f32` figures follow (the correct kernel's `b` curve reads `p = -0.021`).
  `xtask` only; the earlier `f64` rows are unchanged; breaks nothing.
- Harness: `Precision::F32` scoring (`0016` item 2). `cargo xtask conformance --precision f32` scores
  `f32` subjects on the `@f32` strata alone, in units of `2^-24`, with inputs and outputs asserted
  exactly binary32. An id named with `--fn` that has no `@f32` strata is an error, not a skip; without
  `--fn` the run scores the eight ids that have them and names the others a selected subject supports.
  A subject with no `f32` kernel (the planted `c`, a sweep candidate) is an error at `f32`, never the
  D12 answer under its name. The seeded coefficient kernels run at `f32` through the same adapter (the
  D12 prior, not a correct kernel, until the `f32` sweep of `0016` item 3); the correctly rounded and
  neighbouring-ulp sanity subjects have `f32` variants; `--self-test` runs the error-curve and
  non-finite mechanisms at both precisions (`b` without its series fits p = 1.937 at binary64, 2.134
  at binary32; the `f32` range, kernel and floor are questions 26 to 28 of draft record 0014).
  `xtask` only; nothing breaks.
- Corpus: `@f32` strata for `coeff_k`, `coeff_a`…`coeff_e` and `coeff_r` (`0016` item 1: the binary64
  inputs rounded to nearest-even binary32, the reference recomputed at the rounded inputs), and the
  id `coeff_cos_half` (cos θ/2) with its series as the eighth row of `coeff_series`; `0016` item 1
  now names `coeff_cos_half` and `coeff_r` (`PHASE1.md` §4.3, §4.4 follow). Existing records are
  byte-identical; the harness skips `@f32` strata and the binary64 sweep excludes them until `0016`
  items 2 and 3.
- Decision records `0016` (`f32`-exact `@f32` strata for the scalar coefficient ids) and `0017`
  (`Real::cbrt`, mask-valued `solve_cubic` roots), both `ready`: docs only, no code change.
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
- `docs/maths/geodesics.md`: the geodesic of `NUMERICS.md` §10 (the curve, its constant body and spatial velocity,
  what "geodesic" claims on SE_N(3): the autoparallel of the canonical connection and, on SE(3) only, the geodesic of its indefinite invariant forms (no nondegenerate one for N >= 2, no Riemannian one);
  left-, right- and inversion-invariance; `Product<SO3, R3>` as (slerp, lerp), bi-invariant for its own law and
  failing right-invariance only when composed as SE(3) poses, with the exact gap; the Jacobians with respect to
  `X0`, `X1` and `t` derived from those of `Exp`, their left forms and cancellation-free forms of `J0` and its left form; unit
  dual quaternions, the screw parameters, the power, its sign rule and the grouped one-`atan2` evaluation of
  `tf_tree_math`'s `screw_pow`, with its equivalence to the geodesic; where the twin cancels or branches (no
  cancellation in `f64` value at small angle or near π, the per-type 0/0 guard, the `Dual` translation derivative that loses about 1/α whatever series is used, the endpoints, the cut at π); SO(3) slerp).
  Documentation only; no code, no formula change, no normative document changed. Five open items for the
  normative documents are added in the maths index.
- `docs/maths/charts.md`: charts of `PHASE5.md` §1–§2 (a chart as a frozen retraction with a local inverse and
  its two Jacobians read in the chart of the retracted point, `lj(ret δ) = rj(δ)⁻¹`; `RightChart` and `LeftChart` from
  the rows of `NUMERICS.md` §2.3; `Screw`, `Decoupled` and `WorldTranslation` as two retractions in three frames, with
  the block-diagonal `rj`, `lj` of the last two derived and why they are not dual matrices; the change-of-chart
  identity, which contains `se3.md` SE.9(b), (c); what a chart change alters in an LM iterate: not the step of `Screw`
  and `Decoupled`, but the trial point by at most `θ‖ρ‖/2`; S²: the Householder basis and its closed form, `retract`
  as the sphere's exponential of the rotation vector, `local` and the `r` kernel, `rj = ZᵀA(δ)`, `lj = A(δ)⁻¹Z`, the
  reflection at `n_z = 0`, the hairy-ball argument for freezing the chart, the antipode). Documentation only; no
  code, no formula change, no normative document changed. Four open items for the normative documents are added
  in the maths index (`PHASE5.md` §1.3's first-order claim is false for `WorldTranslation` unless `R = I`; `sgn(-0)`;
  `mp.logm` as the `s2_local` reference; the unit-norm domain of `S2`) and two extended.
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
- `helicoid`: the traits `Tangent`, `Jac` and `LieGroup`, the sealed `Side` with its only two
  implementations `Left` and `Right`
  (`docs/PHASE3.md` §2), the conventions of `0002` in the crate docs, and the trivial group
  `Rn<S, N>` with `RnTangent` and `RnJac` (§7): addition is `Mul`, there is no `Add` or `Sub`
  (`compile_fail` doctests), and every row of `NUMERICS.md` §2.3 is written for both sides.
  `Jac::sandwich` and the `DOF` tie are `const` assertions. Generic law checks (test-only
  `laws`) run for `Rn` and for a test-only non-abelian group (Heisenberg) under `f64`, `f32` and
  `Dual<f64, 3>` with recorded bounds. `docs/PHASE3.md` §7 and `docs/API.md` §3 name `RnTangent`
  and `RnJac`. New public API, nothing breaks. `proptest` becomes a dev-dependency of `helicoid`.
- `helicoid`: `Quat<S>` (`docs/PHASE3.md` §4): Hamilton product as `Mul`, `conjugate`, `norm_sq`,
  `to_matrix` (its operand grouping pinned by a bit-exact golden, D16), `from_wxyz_unchecked`
  (debug-asserts `|‖q‖² - 1| <= 2^-40` for `f64`, `2^-16` for `f32`), `from_wxyz_normalized`
  (divides by the norm, `0027`; domain `‖q‖²` normal, `debug_assert!` on the result) and
  `renormalize` (the first-order Newton step, drift repair, asserting nothing), no `Add`, `Sub`
  or `PartialEq` (`compile_fail` doctests), and the converters
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
- `just test` also runs the workspace tests in the release profile, where `debug_assert!` is
  compiled out: the tests of the documented out-of-domain behaviour are `cfg(not(debug_assertions))`
  and the dev profile skipped them. Nothing breaks.
- `helicoid`: `Product<A, B>` and `ProductJac<JA, JB>` (`docs/PHASE3.md` §7): the product group
  with every `LieGroup` method the factors' method componentwise, tangent `(A::Tangent, B::Tangent)`
  (a `Tangent` impl for a pair, dense order `A` then `B`, `dot_acc` threaded so that `dot` stays
  the one flat index-order sum, `0025`), `DOF = A::DOF + B::DOF`, and the block-diagonal Jacobian as
  a `Jac` (`sandwich` from the factors' `apply`). Nests. Neither type has a public field, as §7 and
  `0025` decision 5 state, so a `Product` is built and read through `exp`/`log` and
  `write_dense`/`read_dense`. `helicoid::reference::productjac_sandwich` is its dense twin, with
  `productjac_sandwich_matches_reference`. Tested by the generic laws on `Product<Rn<3>, Rn<2>>`,
  `Product<Rn<2>, Heis>`, `Product<Heis, Rn<2>>` and a nested product with two Heisenberg factors
  (test-only, the one non-abelian factor until SO(3)). New public API, nothing breaks.
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
  and proposed by draft record 0014. Backward error and the `helicoid` subject are not implemented
  (`f32`, `--self-test` and the seeded subjects follow in the entries below). `xtask` gains
  `num-bigint` and `helicoid-linalg`; no library code changes; breaks nothing.
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
- `cargo xtask thresholds` also writes `xtask/src/seeded/generated.rs` (`docs/PHASE1.md` §6, `0004`
  item 3), registered with `cargo xtask lint`: per coefficient a `Switch<f64, terms>` (the switch as a
  bit pattern with its decimal in a comment, the series as the shortest decimal that reads back as
  each exact rational rounded once at binary64, laid out as `rustfmt` leaves it, an objective line)
  under the `@generated` header and the SHA-256 of the sweep CSV and of `coeff_series.jsonl` as the
  source revisions. `just thresholds-check` (the CI job) now fails on any difference in the CSV or in
  that file and writes nothing; tests edit one character of each and drive the entry point under a
  scratch root, and the check names the file and the line (an edit that no longer compiles stops
  `xtask` building: `git restore` the file). `seeded:correct` runs the generated switches, where it
  ran the `tf_tree` D12 candidate (still the sweep's prior, its errors still pinned): its per-field
  maxima equal the CSV's, to the bit. `seeded:c-two-terms-1e-8`, `c` with switch `1e-8` and two
  terms, joins `--self-test`: detected when the sweep's objective for its candidate is more than
  `10^6` times the chosen one's (measured `10^30`, rank 7841 of 8200), while the correct kernel,
  whose `c` is the chosen candidate, stays silent; the self-test also fails when a subject's
  measured `c` objective is not the ranked one, to the bit. `xtask` gains `sha2` (the revisions;
  `xtask` only, `0007`); no library code changes; breaks nothing.
- Seeded SO(3) subject and its two `Log` defects (`docs/PHASE1.md` §10, `xtask/src/seeded/so3.rs`):
  `seeded:correct` now also runs `Exp` (`k` at the generated switch, `cos θ/2` by its exact arm: no
  series is committed for it) and `Log` (the `copysign` flip, `r = 2 atan2(n, w)/n` at the safe
  argument, its series arm `2/w` taken where `n² = 0`, no constant typed) over `so3_exp` and
  `so3_log`, at most 2.91 `u` on every stratum (`--self-test` gates 4), non-finite-free, its safe
  arguments checked on a mask that evaluates both arms. `seeded:log-acos` and `seeded:log-no-flip`
  join `--self-test`, which fails unless each fires its own mechanism and the correct kernel fires
  none: `acos` when every stratum of `theta:1e-k` from `k = 4` and of `theta:pi-1e-k` from `k = 7`
  reaches `10^7` (§10's per-`k` wording cannot hold for every `k`); the flip when every stratum with
  a negated half either fails or has nothing to detect (`n² = 0`, derived from the records: 27 of 29
  fail, `theta:exact0` and `theta:subnormal` cannot), and one fails. A stratum with no score, a
  non-finite output and a `NaN` fail the gate. The two defects run over `so3_log` only. The readings
  are 0014 (draft) questions 16 to 20. No new dependency, no library code changes; breaks nothing.
- Seeded SE_N(3) subject and its two defects (`docs/PHASE1.md` §10, `xtask/src/seeded/se3.rs`):
  `seeded:correct` now also runs `Exp` (`so3::exp` and `J_l(φ)ρ_i`) and the Jacobians `J_l` and
  `J_r = J_l(−τ)` with Barfoot's block `Q` of `NUMERICS.md` §5.3 over `sen3_{exp,jr,jl}_n{1,2,3}`
  (`k, a, b, d, e` each at its own generated switch, dense column-major output): at most 3.47 `u`
  (`Exp`) and 8.90 `u` (`J`) on every stratum (`--self-test` gates 12), no non-finite output; the
  errors are of this per-coefficient kernel, not of the grouped one of `NUMERICS.md` §4.
  `seeded:se3-exp-translation-first` and `seeded:q-minus-half` join `--self-test`, which fails
  unless each fires its own mechanism and the correct kernel fires none: every `rho:*` stratum of
  `sen3_exp_n{1,2,3}` reaches `10^7` `u`, and so does every stratum of every `sen3_jr_n*` and
  `sen3_jl_n*` (measured at least 8.9e15 and 5.2e9). §10's `Dual` comparison of the `Q` defect
  waits for the `helicoid` subject. The readings (the bar, "fails", the strata counted, the
  per-coefficient switches, the association of `Q`'s words) are 0014 (draft) questions 21 to 25. No
  new dependency, no library code changes; breaks nothing.
- `helicoid-linalg`: `Real::cbrt` (`docs/PHASE2.md` §2 and §3, `docs/API.md` §2, `docs/NUMERICS.md`
  §12; step 1 of `0017`, which `solve_cubic` needs): the real cube root, total and odd, through
  `libm::cbrt`/`cbrtf`. `Dual` differentiates it as `d / (3 c^2)`, `c = cbrt v`: infinite at 0
  like `sqrt`, with the sign of `d` at either zero, accurate elsewhere, `0` at `+-inf`.
  Checked against mpmath on chosen rows (both signs, subnormals; first derivative at `f64` and
  `f32`, second derivative and the Hessian of `cbrt(x y)` at `f64`) to 2 ulp; random inputs reach
  about 4 ulp on the first derivative and 10 on the second, the value 0 ulp. **Breaks** an
  out-of-tree `impl Real`, which gains a required method (none is known; the in-tree ones are
  updated); nothing else breaks.
- `helicoid-linalg`: `solve_cubic` (`docs/PHASE2.md` §6, step 2 of `0017`): the real roots of `a x³
  + b x² + c x + d` as `(Vec3<S>, [S::Mask; 3])`, slot `k` reported as a root iff its mask is set, a
  clear slot `+0`. omnisac's algorithm, generic over `S: Real` (Cardano, Viète's trigonometric form,
  a discriminant band), every branch a `branch`/`select` at a safe argument, so it runs on `f64`,
  `f32`, lanes and `Dual`. Differences from omnisac's `poly::solve_cubic`: the container (no
  `ArrayVec`, no compaction); the tolerances, `1e-14` to `2^-46` (128 `u`, from 90 `u`) and `1e-12`
  to `2^-40` (8192 `u`, from 9007 `u`) at `f64`, `2^-17` and `2^-11` at `f32`, the nearest
  power-of-two multiples of `u` (a polynomial whose leading coefficient or discriminant lies in the
  thin strip between old and new changes its root count); `acos` as `atan2(sqrt((1 - x)(1 + x)),
  x)`, within 2 ulp of `libm::acos`, and `pi` as `atan2(+0, -1)` — a reading `0022` has since
  reversed, and `pi` is already a literal: the roots of the trigonometric arm move by up to about
  11 `u` of the largest root, an observed maximum and not a bound, which `0022` step 1 removes. Everything else agrees to the bit with omnisac (on `libm`) on
  10^6 random and adversarial polynomials per seed (a scratch differential, not committed; `GOLDEN`
  in the tests is 36 of its rows). Tested against omnisac's cases, an mpmath fixture, planted roots
  against a measured error model, the slot order, the `Dual` derivative, two lanes and the mask.
  Limits it inherits, in the rustdoc `# Domain` and `0031` (draft), fixed by none of this: the
  one-real-root arm is not backward stable (`x^3 + p x - 1` is off by `4e-6` at `p` near `1e-5`,
  `f64`, and its `Dual` derivative is `-inf` below `p = 1.3e-5`); a repeated root can be dropped
  where `p` and `q` cancel; a cubic with one real root can read as three where `p^3` and `q^2`
  underflow. `docs/API.md` R6 names `solve_cubic` beside `chol` as the functions that assert
  nothing. New item: breaks nothing; omnisac's call sites need an adapter (`docs/PHASE2.md` §9).
- `helicoid-linalg`: `eig3` (`docs/PHASE2.md` §6, step 3 of `0017`): the eigenvalues, ascending, and
  an orthonormal basis of eigenvectors with `det = +1` of a symmetric `Mat3<S>`, as `(Vec3<S>,
  Mat3<S>)`; the lower triangle is read. omnisac's `eigendecomp_sym3_signed` (d3be7b7), generic over
  `S: Real` (Smith's eigenvalues, cross-product eigenvectors; `f64`, `f32`, lanes, `Dual`). What
  differs from omnisac is proposed in the draft record `0023`: no `Option` and no tolerance, so the
  eigenvalues are not finite for a non-finite entry and where `p^3` over- or underflows, and a frame
  built from cross products of unit vectors, orthonormal to a few `u` for every input. **Not done:
  Kopp's hybrid is owed** (§6). Where the two largest eigenvalues are within about `2 sqrt(u)|A|`
  (equal ones included) no column is reliable, the isolated smallest eigenvalue's included (residual
  of the order of `|A|`); omnisac's vectors are as wrong there, or `None` where a cross product
  vanishes, and anchoring the frame on the isolated end may fix it (`0023` (draft)). At a double
  root the eigenvalues lose half their digits (`0.6 sqrt(u)|A|` on a rank-1 matrix), which only the
  hybrid fixes. §6 names no constant for its bars, so none is claimed; the measurements (a scratch
  harness, not committed) are in `eig3_tests` and the §0.0 row. `docs/API.md` R6 names `eig3` beside
  `chol` and `solve_cubic`. New item: breaks nothing; omnisac's call sites need an adapter
  (`docs/PHASE2.md` §9).

### Changed

- `SEn3Jac::sandwich` builds both halves with `from_cols`/`from_rows` instead of filling and then
  overwriting two `D x D` scratch matrices: `M = J Σ` column by column, then `M Jᵀ` row by row, since
  row `r` of `M Jᵀ` is `J` applied to row `r` of `M`. It therefore forms neither `Aᵀ` nor the `N`
  transposes of `col`, reaches no entry through `Matrix::get`/`set`, and emits no `memset` where the
  block form emitted two of `D²` scalars that it then fully overwrote. The private `zeros`, `block`
  and `put` of `dualmat.rs` are deleted with it, which withdraws the `dualmat.rs` half of
  `PROJECT.md` §5.1's gate on `Matrix::block`/`set_block`.
  **Bit-identical**: every entry is the same products summed in the same order, because `Matrix`'s
  `Mul` and `Mul<Vector>` associate identically and `B_i w_0` is added to `A w_{i+1}` in `apply`'s
  order — checked over 150 000 random `(J, Σ)` at `f64`, `f32` and `Dual<f64, 3>` for `N = 1, 2, 3`
  with not one differing entry, so D16 holds and the recorded bounds are unmoved rather than
  re-measured. Release, x86_64, minimum of 7 runs of 200 000 calls: 545 → 48 ns (`f64`, `D = 6`),
  2406 → 448 ns (`f64`, `D = 12`), 2301 → 606 ns (`f32`, `D = 12`), 1747 → 680 ns and 8078 → 2328 ns
  (`Dual<f64, 3>`, `D = 6` and `12`). The result is fully packed: 162 `mulpd` and no `mulsd` for the
  324 scalar products at `N = 1`. A middle variant that built each block's `Vector` and gathered it
  was faster still on `f64`/`f32` (365 and 327 ns at `D = 12`) but 17% to 20% slower on `Dual`, and
  `0006`'s per-precision no-regress bar refuses it. Values unchanged, so nothing breaks.
- `SEn3Tangent::read_dense` and `write_dense` take one length test instead of one bounds check per
  component, on the path `ProductJac::sandwich` crosses once per column and once per row of its
  argument and again per nesting level. At `N = 3` in release on x86_64 the in-domain arm of
  `read_dense` is one `cmp`/`jb` and twelve loads, and `write_dense`'s is one `cmp`/`jb` and a
  `movups`-packed copy with no `memcpy` call — the copy vectorizes once the `Chain`/`FlatMap` of
  `comps` is gone, which
  `dot_acc` had already avoided for its own reason. Out-of-domain behaviour is unchanged and still
  pinned by `out_of_domain_does_not_panic_in_release`: a short `write_dense` view takes the prefix
  that fits, and `read_dense` poisons only the missing entries, never the present ones
  (`0025` decision 4).
- `dualmat_tests`' `kappa` takes `M⁻¹` from the closed form `A⁻¹ − ε A⁻¹ B_i A⁻¹` rather than from the
  `sen3jac_inverse` twin. The twin `debug_assert!`s a nonzero pivot, so the `κ u <= 1e-3` filter was
  computed from a call that panics on exactly the singular draws the filter exists to reject —
  `inverse_twin::<f64, 1, 6>(&[0.0; 18])` panicked instead of returning `None`, and because proptest
  shrinks toward `0`, that panic and not the bound violation is what a real regression would have
  reported. `Mat3::inverse_adj` asserts nothing, so a singular `A` now yields a non-finite `κ` that
  fails the filter like any other out-of-domain draw, and the case costs one 3x3 inverse instead of a
  second `D³` elimination. All twelve recorded bounds reproduce from the documented `10^6`-case run.
- `laws`' `Rng::vec` becomes `Rng::arr::<N>`, filling a caller-owned array with an explicit loop: the
  `10^6`-case measurements no longer allocate once per draw-set per iteration, and the draw order is
  a written loop rather than an unspecified visiting order, since that order is the stream every
  recorded figure is reproducible from.
- `laws`' dense reader is `reference::dense_oriented`, which `reference::dense` also calls: the
  NaN-poisoned scratch had two definitions differing only in orientation, and the poison semantics
  now have one place to be kept. `reference::sum`'s doc states what couples it to
  `helicoid_linalg`'s `vector::sum` and which tests assert the two associate identically.
- `zero3`, `put_view` and the `dualmat.rs` helpers that survive are `#[inline]`; the comment on
  `sandwich`'s scratch no longer claims a saving it did not make.

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
- `Tangent::dot` becomes a provided method over the new required `dot_acc`, and
  `Tangent::read_dense` reads a missing entry as NaN rather than `+0` (`0025`): `+0` is a valid
  component, so a wrongly sized buffer was indistinguishable from a zero tangent.
- The law harness divides every norm-wise error by a NaN-preserving denominator. `f64::max` returns
  the other operand for a NaN, so a NaN reference had been scaling away to a denominator of `1`,
  and a law could pass against it.
- The non-abelian test group's sample Jacobian is no longer a circulant. Circulants commute, so
  `jac_dense_order` — which pins `Jac::mul`'s operand order against the dense product — could not
  detect a `mul` that multiplied its operands the wrong way round. Row-scaling breaks the
  commutation; the recorded `f64` bounds move to `jac_order` 12 and `sandwich` 12, `f32` to 11 and 9.
- `Tangent::read_dense`/`write_dense` take the exact-length path through `first_chunk`, so the
  documented boundary costs one length test instead of one per component (measured at `DOF = 9` on
  the emitted release asm: 1 branch and vectorized moves, against 17 and 9 scalar ones).
