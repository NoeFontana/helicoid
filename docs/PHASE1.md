# helicoid — Phase 1 Implementation Specification: Skeleton and the Instrument

> Where this spec is silent, check `docs/PROJECT.md`'s decision log first. Sections marked
> **NORMATIVE** are requirements; code blocks are illustrative, but file layouts, schemas and
> command names are normative.

**Deliverable:** a workspace that builds on four targets; a committed mpmath corpus; a conformance
harness; a threshold sweep; oracle runners; the envelope and its bars; a bench gate — **all shown to
catch planted defects — and no group code.** The only numeric subjects in this phase are the seeded
defects and `tf_tree_math`, which becomes oracle #1 ([`0010`](./decisions/0010-seeded-from-tf-tree-math.md)).

## 0.0 Implementation status

| Area | Status |
|---|---|
| Workspace, lints, `justfile`, CI matrix (§3) | Partial: workspace, lints, `justfile` and CI for `build`/`test`/`lint`/`audit`/`msrv`/`no-std`/`wasm`/`doc`/`corpus-check`/`conformance --self-test`/`thresholds-check`/`envelope`. Determinism and container oracles pending. |
| Reference generator and committed corpus v1 (§4) | Partial: generator skeleton and committed corpus v1 with 150-digit mpmath reference generation across scalar and group strata. SE(2) coefficients pending. |
| Conformance harness and result schema (§5) | Partial: `cargo xtask conformance` runner, result CSV generation, and seeded defect checks. Backward error and f32 vector strata pending. |
| Threshold sweep and generated-file format (§6) | Partial: two-stage threshold sweep over full domain generating `coeffs` and series tables (`0039`, `0047`). |
| Oracle runners: `tf_tree_math`, sophus-rs (excluded crates); Sophus, manif, GTSAM (containers) (§7) | Partial: runners implemented for `tf_tree_math` and `sophus_rs` over file protocol. Containerized oracles pending. |
| Envelope and bars (§8) | Partial: domination and no-regression gating over baseline and oracles (`0006`). Baseline blessing awaiting Phase 3 candidate rows. |
| Bench harness and gate (§9) | Partial: criterion benchmarks and replicate-outermost gate (`0033` draft). Bench-check CI job pending. |
| Seeded-defect validation (§10) | Partial: `just conformance --self-test` runs defect injection for coefficient, SO(3), and SE(3) defects. |

## 0. Non-goals and guardrails — read first

**NORMATIVE.** Do not implement in Phase 1: any group, any Jacobian, any `Real` impl
(`helicoid-linalg` is Phase 2). The two library crates exist with their lints, `no_std` roots and
empty public surfaces. The instrument must exist, and must be shown to *detect*, before the first
thing it measures ([`0006`](./decisions/0006-the-instrument-comes-first.md)).

**Harness independence.** The generator shares no code with any Rust subject and never evaluates a
closed form from `NUMERICS.md` (§4.3). **Dependency isolation.** `xtask` and the runners may
depend on anything; no library crate depends on either, and the runners are excluded from the
workspace so their dependency graphs never enter `Cargo.lock`.

**If a design question is not answered by this document, stop and ask.**

## 1. Workspace layout

```
helicoid/
├── Cargo.toml              # workspace; [workspace.package], [workspace.lints]
├── justfile  deny.toml  rust-toolchain.toml  CONTRIBUTING.md  LICENSE-MIT  LICENSE-APACHE
├── crates/
│   ├── helicoid-linalg/    # Phase 2; empty lib in Phase 1
│   └── helicoid/           # Phase 3; empty lib in Phase 1
├── xtask/                  # corpus, conformance, thresholds, envelope, bench-gate, lint, determinism
│   └── src/seeded/         # planted defects (§10)
├── runners/                # [workspace] exclude — own Cargo.lock each
│   ├── tf_tree_math/       # git dep pinned to a tf_tree commit
│   └── sophus_rs/          # pinned crates.io version
├── conformance/
│   ├── generate/           # pyproject.toml + uv.lock (pinned Python, mpmath); gen/ package
│   ├── corpus/             # committed *.jsonl + MANIFEST.json
│   ├── sweeps/             # committed sweep CSVs (Phase 3 onward)
│   ├── baseline/           # committed per-stratum maxima: the no-regress bar
│   └── results/            # gitignored
├── docker/oracles/{sophus,manif,gtsam}/   # Dockerfile + runner, file protocol (§7)
├── baseline/               # bench baselines + HOST.md (§9)
└── docs/                   # this book; docs/evidence/ENVELOPE.md is generated (§8)
```

## 2. Load-bearing invariants

**NORMATIVE.** Encode each as a check in the named recipe.

1. **Exact inputs.** Every corpus input is a binary64 value serialized with Python's `float.hex()`;
   the reference is computed at that exact value (`corpus-check`).
2. **Definitions, not formulas** (§4.3). The generator computes `Exp` as a matrix exponential or a
   quaternion power series, `Log` as the inverse of that definition (Newton's method on it;
   `mp.logm` returns complex results near $\pi$, so it only cross-checks), Jacobians as their
   defining series, inverses by high-precision linear algebra, coefficients from their definitions.
3. **Precision budget.** `mp.dps = 120`. The worst catalogued cancellation (`e`, `NUMERICS.md` §4)
   loses about $4\log_{10}(1/\theta) + 3$ digits: 51 at $\theta = 10^{-12}$, leaving 69 against the
   17 a binary64 reference needs. Below $10^{-12}$ the budget is exceeded (`theta:subnormal`: every
   raw definition of $a, \dots, e$ evaluates to 0 at 120 digits), so 120 is the precision of the
   *output*: a cancelling definition is evaluated at guard digits, doubled until two evaluations
   agree (`coeff.stable`). A 1% sample is recomputed at `dps = 150` and must agree to 40 digits, or
   generation fails.
4. **Byte-identical regeneration.** Pinned Python and mpmath (`uv.lock`), a seeded splitmix64
   implemented in the generator (not `random`, whose algorithm is a CPython detail), sorted keys,
   fixed float formatting. `just corpus-check` regenerates into a temporary directory and `cmp`s.
5. **Per-stratum max and p99, never a mean.** Units of $u$ (`NUMERICS.md` §2.1, §11).

## 3. Toolchain, lints and CI

**NORMATIVE.**

- `[workspace.package]`: `version = "0.0.0"`, `edition = "2021"`, `rust-version = "1.87"`,
  `license = "MIT OR Apache-2.0"`, `repository`/`homepage` `https://github.com/NoeFontana/helicoid`.
- `[workspace.lints.clippy]`: `unwrap_used`, `expect_used`, `panic`, `todo`, `unimplemented`,
  `dbg_macro`, `float_cmp` = `deny`; `print_stdout`, `print_stderr`, `redundant_clone`,
  `needless_pass_by_value` = `warn`; `upper_case_acronyms` = `allow` (type names `SO3`, `SE3`,
  [`0008`](./decisions/0008-the-name-helicoid.md)). `[workspace.lints.rust]`: `missing_docs` =
  `warn`, `unreachable_pub` = `warn`, `unused_must_use` = `deny`. `xtask` and tests may `allow`
  `print_stdout`, `float_cmp` and `panic` locally.
- Library roots: `#![no_std]`, `#![forbid(unsafe_code)]`, `#![deny(missing_docs)]`. No `alloc`
  crate import.
- CI jobs are `just` recipes 1:1: `build`/`test` (x86_64-linux, aarch64-linux), `no-std`
  (`thumbv7em-none-eabihf`), `wasm` (`wasm32-wasip1` under wasmtime), `msrv`, `lint`, `doc`,
  `corpus-check`, `conformance --self-test`, `envelope` (Phase 3 onward), `determinism` (wired
  now, required from Phase 6), `oracles` (nightly, containers), `bench-check` (opt-in label).
- `cargo xtask lint` checks: no new `path.rs:NNN` citation in docs or comments; every
  `@generated` file matches its generator; the `NUMERICS.md` §14 twin table against
  `helicoid::reference` and proptest names (Phase 3 onward); **the normal-dependency closure of
  each library crate equals the allowed set** of [`0007`](./decisions/0007-the-budget-a-foundation-can-afford.md);
  no crate other than `xtask` enables `__sweep`.
- Owed files: `CONTRIBUTING.md`, `LICENSE-MIT`, `LICENSE-APACHE`, `deny.toml` (licenses; bans on
  `nalgebra`, `faer`, `num-traits`, `bytemuck`, `serde` for library crates).

## 4. The reference generator

**NORMATIVE.**

### 4.1 Layout and invocation

`conformance/generate/` is a `uv` project: `pyproject.toml`, `uv.lock` pinning the interpreter and
`mpmath` (no numpy). `just corpus` runs `uv run --frozen python -m gen all`; `just corpus-check`
runs it into a temporary directory and compares every file byte for byte with `conformance/corpus/`.

### 4.2 Record schema

One JSONL file per function id; one record per line (values below are illustrative, not a test vector):

```json
{"id":17,"stratum":"theta:1e-8","in":{"phi":["0x1.4f8b588e368f1p-27","-0x1.0c6f7a0b5ed8dp-28","0x1.b7cdfd9d7bdbbp-29"]},"out":{"q":["9.99999999999999982500000000000e-1","7.00000000000000009125814376540e-9","-2.49999999999999994652465432315e-9","1.62499999999999996822012103312e-9"]}}
```

Inputs are hex floats; outputs are decimal strings with 30 significant digits. Matrices are
column-major flat arrays with a sibling `"shape"`. `MANIFEST.json` records the generator identity (a
SHA-256 over its sources: a git revision changes with every commit and cannot appear in the commit
it names), Python and mpmath versions, `dps`, seed, and per file SHA-256, record count and `kind`:
`corpus` (the records above; `rechecked` counts those recomputed at 150 digits), `series`
(`coeff_series`, §4.3; `verified` counts its series equal to exact algebra) or `switch-ref`
(`coeff_switch_ref`, §6: the true value and `d/dz` of every coefficient at every grid point,
keyed by the point's binary64 bit pattern; `rechecked` counts those recomputed at 150 digits).
Neither of the last two is a function id, so neither is a row of §4.3 and neither is read by a
bar: `coeff_series` supplies the generated files' literals and `coeff_switch_ref` the right-hand
side of `docs/maths/coefficients.md` CO.12.

### 4.3 Definitions (v1 function ids)

| Function id | Input | Reference computed as |
|---|---|---|
| `coeff_k`, `coeff_a`, `coeff_b`, `coeff_c`, `coeff_d`, `coeff_e` | $\theta$ | the definition at 120 digits; derivative by `mp.diff` |
| `coeff_cos_half` | $\theta$ | $\cos(\theta/2)$, the other half of `Exp`'s quaternion (`NUMERICS.md` §3.1), at 120 digits; derivative by `mp.diff` |
| `coeff_r` | $(n, w)$ | $2\,\mathrm{atan2}(n, w)/n$ at 120 digits |
| `coeff_alpha` | $\theta$ | $\sin\theta/\theta$ at 120 digits; derivative by `mp.diff`; no series and no switch of its own ([`0062`](./decisions/0062-sin-theta-over-theta-is-public-and-differentiates-its-branch.md)) |
| `coeff_series` | — | the exact series of [`0004`](./decisions/0004-switch-points-are-generated-not-typed.md): `mp.taylor` of each definition at 120 digits, rationalized, equal term by term to an independent exact derivation |
| `so3_exp` | $\varphi$ | quaternion power series $\sum p^n/n!$, $p = (0, \varphi/2)$ |
| `so3_log` | $q$ | the $\varphi$ with $\|\varphi\| \le \pi$ and $\mathrm{Exp}(\varphi) = q/\|q\|$ ($w \ge 0$ after the flip), by Newton's method on the `so3_exp` series; `mp.logm` of $R(q/\|q\|)$ is a test cross-check only, since it returns complex results near $\pi$; stratum `q:w0` uses the sign rule of `NUMERICS.md` §3.2 |
| `so3_act` | $(q, p)$ | $R(q/\|q\|)\,p$ |
| `so3_from_matrix` | $R$ | the quaternion of the rotation nearest to $R$ in Frobenius norm (the polar factor): a rounded or scaled $R$ is not a rotation, so its `mp.logm` is not skew; backward error only (`NUMERICS.md` §11) |
| `so3_jr`, `so3_jl` | $\varphi$ | $\sum (\mp W)^n/(n+1)!$ |
| `so3_jr_inv`, `so3_jl_inv` | $\varphi$ | `mp.inverse` of the above |
| `so3_geodesic` | $(q_0, q_1, t)$ | $\hat q_0\,\mathrm{Exp}(t\,\mathrm{Log}\,\Delta)$ with $\hat q = q/\lVert q\rVert$ and $\Delta = \hat q_0^*\hat q_1$: the **geometric** $\mathrm{Log}$ (`so3_log`'s quaternion route) and the `so3_exp` series. **Not `mp.logm`**, which returns a complex, non-principal logarithm from $\theta = 3.03$ and so is wrong by $O(1)$ on `geo:near-pi` ([`0045`](./decisions/0045-two-phase-4-checks-cannot-be-taken-as-written.md)); it is the generation-time cross-check for $\theta \le 3$, where the two agree to $8.2\times10^{-113}$ over the committed corpus, and above it the group identities are the check — the endpoints, $\gamma(X_0,X_1,t) = \gamma(X_1,X_0,1-t)$, and $\mathrm{Log}(X(t)^{-1}X(t+h)) = h\,d$ at a fixed $h$, which is constant body velocity without a difference quotient. $X_1 = X_0\,\mathrm{Exp}(d)$ is **formed and rounded, never drawn**, so a stratum's relative motion is the $d$ it names |
| `se3_geodesic` | $(q_0, x_0, q_1, x_1, t)$ | the same over $\mathrm{SE}(3)$: $\Delta = (\hat q_0^*\hat q_1,\ R(\hat q_0)^\top(x_1 - x_0))$, $d$ its `sen3_log_n1` (the $\mathsf V$ solve), the answer $X_0\,\mathrm{Exp}(t\,d)$ with the translation block from `mp.expm` of the $4\times4$ hat matrix as `sen3_exp_n1`'s is. The same cross-checks, on the $4\times4$ |
| `se3_screw_retract`, `se3_decoupled_retract`, `se3_world_retract` | $(q_0, x_0, \tau)$, $\tau = [\varphi; \rho]$ | the chart's `retract` at $X = (\hat q_0, x_0)$ ([`0060`](./decisions/0060-the-charts-go-first-and-name-their-frame.md) decision 8, `docs/maths/charts.md` §3): $X\,\mathrm{Exp}(\tau)$ through `sen3_exp_n1`; $(\hat q_0\,\mathrm{Exp}\,\varphi,\ x_0 + R\rho)$; $(\hat q_0\,\mathrm{Exp}\,\varphi,\ x_0 + \rho)$. $\tau$ from the `SEN3_STRATA` twist, the base Haar-random with $\lVert x_0\rVert = 1$ from the stratum's own stream. Cross-checks: the chart's own `local` returns $\tau$; `screw`'s translation against `mp.expm` |
| `se3_screw_local`, `se3_decoupled_local`, `se3_world_local` | $(q_0, x_0, q_1, x_1)$ | the chart's `local` of $Y = (\hat q_1, x_1)$: $\mathrm{Log}(X^{-1}Y)$ through `sen3_log_n1`; $(\mathrm{Log}(\hat q_0^*\hat q_1),\ R^\top(x_1 - x_0))$; $(\mathrm{Log}(\hat q_0^*\hat q_1),\ x_1 - x_0)$, the geometric `Log`, never `mp.logm`. $Y$ is the chart's own `retract` of the same record's twist, rounded. Held as `phi` and `rho`, scored with floors 1 and $\lVert x_0\rVert$ (two absolute poses cancel). Cross-check: the chart's `retract` of the answer returns $Y$ |
| `sen3_exp_n{1,2,3}` | $\tau$ | `mp.expm` of the $(3+N)$-square hat matrix; the rotation as `so3_exp`'s quaternion series, asserted equal to the matrix's block |
| `sen3_log_n{1,2,3}` | $X$ | the $\tau$ with $\|\varphi\| \le \pi$ and $\mathrm{Exp}(\tau) = X$: `so3_log`'s $\varphi$, then $\rho_i$ solving $J_l(\varphi)\rho_i = x_i$ with $J_l$ the series of the same exponential's translation block; `mp.logm` of the $(3+N)$-square matrix is a test cross-check only, since it is wrong near $\pi$ |
| `sen3_ad_n{1,2,3}` | $X$ | images of the basis under $\sigma \mapsto (X\sigma^\wedge X^{-1})^\vee$ |
| `sen3_jr_n{1,2,3}`, `sen3_jl_n{1,2,3}` | $\tau$ | $\sum (\mp\,\mathrm{ad}_\tau)^n/(n+1)!$, as dense $(3+3N)$-square matrices |
| `sen3_jr_inv_n{1,2,3}`, `sen3_jl_inv_n{1,2,3}` | $\tau$ | `mp.inverse` of the above, the blocks off the diagonal and the first block column set to exactly 0 (the LU leaves rounding residue there; a dual matrix has no entry in them, `NUMERICS.md` §2.2), every other entry the LU's |
| `so2_exp` | $\theta$ | `mp.expm` of $\begin{bmatrix}0 & -\theta\\ \theta & 0\end{bmatrix}$, as the unit complex $(c, s)$ |
| `so2_log` | $z$ | the $\theta$ with $\lvert\theta\rvert \le \pi$ and $\mathrm{Exp}(\theta) = z/\lVert z\rVert$, by Newton's method on the complex exponential series ($z = (-1, \pm 0)$ is not sampled) |
| `se2_exp` | $\tau = [\theta; \rho]$ | `mp.expm` of the $3\times3$ hat matrix |
| `se2_log` | $X = (z, t)$ | `so2_log`'s $\theta$, then $\rho$ solving $V\rho = t$, $V$ the series of the same exponential's translation block |
| `se2_ad` | $X$ | images of the basis under $\sigma \mapsto (X\sigma^\wedge X^{-1})^\vee$ |
| `se2_jr`, `se2_jl` | $\tau$ | $\sum (\mp\,\mathrm{ad}_\tau)^n/(n+1)!$, as dense $3\times3$ matrices |
| `se2_jr_inv`, `se2_jl_inv` | $\tau$ | `mp.inverse` of the above |
| `solve_cubic` | $(a, b, c, d)$ | the three roots of the stored coefficients as `re`, `im`, ordered by $(\mathrm{Re}, \mathrm{Im})$: planted where the coefficients are exact, else `mp.polyroots` with `extraprec=300` ([`0056`](./decisions/0056-the-routines-d7-does-not-reach.md)) |
| `eig3` | $A$, symmetric | `mp.eigsy`: `lambda` ascending, `V` its eigenvectors as columns |
| `chol_n{3,6}` | $A$, symmetric | `valid` decided exactly by an $LDL^\top$ in rationals; `L` = `mp.cholesky`, zeros where `valid` is 0 |
| `chol_solve_n{3,6}` | $(A, b)$ | `mp.lu_solve` of the equilibrated system $DAD$, $D = \mathrm{diag}(A)^{-1/2}$ |
| `quat_renormalize` | $q$ | $q/\lVert q\rVert$, the projection `NUMERICS.md` §3.6 claims inside its band |
| `mat2_inverse_adj` | $A$, $2\times2$ | `det` and `inv` $= \mathrm{adj}(A)/\det$ of the stored entries in rationals, exactly; `inv` zeros where $\det = 0$ ([`0061`](./decisions/0061-mat2-keeps-its-adjugate.md)) |
| `real_sqrt`, `real_cbrt`, `real_acos`, `real_sin_cos`, `real_atan2`, `real_div` | $x$; $(y, x)$; $(n, d)$ | the value and its calculus derivative(s), what `Dual` carries; `mp.diff` at a relative step is the cross-check |

A `coeff_k`…`coeff_e`, `coeff_cos_half` or `coeff_r` record's outputs are `value` and `d_branch`:
the derivative with respect to the branch variable of `NUMERICS.md` §4 (θ² for `cos_half`), by
`mp.diff` at the exact real branch value of the exact input (θ², not `fl(θ·θ)`), binary64, or
binary32 in an `@f32` stratum (§4.4).

`coeff_series` is the one file that does not follow §4.2: one record per coefficient (the seven of
`NUMERICS.md` §4 and `cos_half`), `{"id","coeff","branch","prefactor","series"}` with sixteen
`"num/den"` terms in the branch variable, and no `stratum`, `in` or `out`. Its manifest `kind` is `series`; the conformance
harness, the envelope and the oracle runners read `kind: "corpus"` files only.

Series are summed until the term's norm is below $10^{-110}$ relative. Phase 2's `svd3`
(`mp.svd_r`) is owed; Phases 4–5 add their ids by the spec that needs them.

### 4.4 Strata (v1)

- `theta:1e-k`, $k = 12, \dots, 1$, and `theta:1e0` ($[1, \pi - 0.1)$): 64 samples each, $\theta$
  log-uniform in the decade, axis uniform on S².
- `theta:exact0` ($\varphi = 0$) and `theta:subnormal` (components $\sim 10^{-310}$).
- `theta:pi-1e-k`, $k = 1, \dots, 12$: 64 axes each. (`fl(π)` is below $\pi$ by
  $\approx 1.2\times10^{-16}$: there is no binary64 input with $\theta = \pi$ from $\varphi$.)
- `q:w0`: quaternions with $w = +0$ exactly.
- `theta:dense`: 200 points per decade over $[10^{-4}, 1]$ — switch-point continuity.
- `rho:1e-6`, `rho:1e-3`, `rho:1e0`, `rho:1e3`, `rho:1e4` (translation scale) crossed with
  $\theta \in \{10^{-8}, 10^{-4}, 10^{-1}, 1, \pi - 10^{-6}\}$ for SE_N(3), $N \in \{1, 2, 3\}$.
  Each cell is a stratum, `rho:1e4/theta=pi-1e-6`. The SE_N(3) ids also see the `theta:*` strata at
  unit translation scale, all but `theta:dense`; `sen3_log` and `sen3_ad`, which take a quaternion,
  also see `q:w0` and `q:nonunit`. Each SE_N(3) stratum is 6 records, not 64: five dense matrices
  per record at 64 would be 190 MB. 6 is a budget trade-off (the family is 18 of the 50 MB; the
  rest is left for the later phases' ids), not the largest count that fits.
- `so2_*` see the $\theta$ strata above, each $\theta$ in both signs. `se2_*` see the SE_N(3) strata
  (`theta:dense` excluded), 8 records each, the sign of $\theta$ alternating, $\rho$ a random
  direction on $S^1$ at the scale. They have no `q:w0` or `q:nonunit` analogue (a non-unit $z$, or
  $c = 0$): a new stratum family starts as a decision record, so `so2_log`, `se2_log` and
  `se2_ad` see $z$ only as the rounded $(\cos\theta, \sin\theta)$, $\lvert\lVert z\rVert^2 - 1\rvert \le
  1.4\cdot10^{-16}$, and their renormalisation $z/\lVert z\rVert$ is exercised at that level only.
- `geo:consecutive`, `geo:generic`, `geo:near-pi`
  ([`0045`](./decisions/0045-two-phase-4-checks-cannot-be-taken-as-written.md) item 3): the
  geodesic ids' strata, 6 pose pairs each crossed with
  $t \in \{0, 10^{-9}, \tfrac14, \tfrac12, 1 - 10^{-9}, 1\}$ and 4 uniform on $(0,1)$, so
  60 records per stratum. A stratum fixes the **relative** rotation $\theta(d)$, not a pose's:
  `geo:consecutive` log-uniform in $[10^{-9}, 10^{-3}]$ with the relative translation at
  $\theta$'s own scale, so $\lVert d\rVert$ is that band times $\sqrt2$ — `tf_tree`'s
  kilohertz edge, where the relative motion is **below the ulp of the poses it is between**: at
  $\lVert x_0\rVert = 10^4$ and $\theta = 10^{-9}$ the stored $x_1 - x_0$ keeps about three
  digits, which is the regime and not a defect, the records being the geodesic between the two
  poses they hold. `geo:generic` is log-uniform in $[10^{-3}, \pi - 0.1]$ at unit relative
  translation, `geo:near-pi` is $\theta = \pi - 10^{-k}$, $k \in \{1,2,3,6,9,12\}$, one per
  sample: it **stops strictly below $\pi$ and that $k$ is its margin**, because at $\pi$ the two
  preimages give geodesics $O(1)$ apart, so the stratum measures conditioning and not agreement.
  $\lVert x_0\rVert$ cycles $1, 10^2, 10^4$ within **every** `geo:*` stratum (§4 asks
  `geo:consecutive` for $10^4$; giving the other two the same span costs nothing and is what
  `PHASE4.md` §5.2's left-invariance bound is measured against). $t = 0$ and $t = 1$ are kept and
  are the cheapest rows in the corpus: they catch an endpoint error no interior $t$ sees.
- `q:nonunit`: quaternions at $\|q\|^2 - 1 = \pm 2^{-45}$ (`Log`'s scale invariance, `from_*` normalization).
- `<S>@f32` ([`0016`](./decisions/0016-f32-exact-strata-for-the-scalar-coefficient-ids.md)), for
  every stratum $S$ of `coeff_k`, `coeff_a`…`coeff_e`, `coeff_cos_half` and `coeff_r`, after all of
  them: $S$'s inputs rounded to nearest-even binary32 (`coeff_r`: $n$ and $w$ each), the reference
  recomputed at those inputs. `theta:subnormal@f32` is its own decade $[10^{-40}, 10^{-39})$
  (rounding $[10^{-310}, 10^{-309})$ gives 64 zeros); `theta:pi-1e-8`…`theta:pi-1e-12` round to one
  binary32, kept five times.
- The `cubic:*`, `eig:*`, `chol:*`, `renorm:*`, `x:*`, `yx:*` and `nd:*` strata of `0056`'s ids
  ([`0056`](./decisions/0056-the-routines-d7-does-not-reach.md) decision 2, which lists each):
  planted roots, spectra and conditioning, the drift band of `NUMERICS.md` §3.6, and `Dual`'s
  singular neighbourhoods, each with its `@f32` counterpart after every binary64 stratum, drawn
  from its own stream wherever binary32's range differs.

Corpus v1 stays under 50 MB uncompressed; a family that would exceed it reduces its sample count
by an edit to this list, never silently.

## 5. The conformance harness

**NORMATIVE.**

- **Subjects.** `xtask::conformance::Subject` has `name()`, `supports(fn_id)`, and
  `eval(fn_id, &Record, Precision) -> Output`. In-process subjects: `helicoid` (Phase 3 onward),
  `seeded:<defect>` (§10). Oracles are out of process (§7).
- **Metrics:** `NUMERICS.md` §11, per record; aggregated per `(fn, stratum, precision, subject)`.
- **Result CSV:** `fn,stratum,precision,subject,subject_version,n,max_u,p99_u,argmax_id,nonfinite,git_rev`.
- `just conformance` prints the table sorted by `max_u` descending and fails on any `nonfinite > 0`.

## 6. The threshold sweep

**NORMATIVE** ([`0004`](./decisions/0004-switch-points-are-generated-not-typed.md)).

- **Scope.** Every coefficient in `NUMERICS.md` §4, plus SE(2)'s $\alpha, \beta$; each for `f64`
  and `f32`.
- **Grid.** Series terms $m \in \{1, \dots, 16\}$ — the corpus's own series length (§4.3) — and the
  switch point in the branch variable on a log grid of 64 points per decade spanning
  $z \in [10^{-16}, 10]$, from the integer 64th root of a power of ten so that a point is exact at
  both precisions ([`0016`](./decisions/0016-f32-exact-strata-for-the-scalar-coefficient-ids.md)
  item 3). **The grid spans the domain and the selection rule stops below it**
  ([`0039`](./decisions/0039-the-sweeps-grid-stops-below-its-own-optimum.md) item 7): a switch at
  or above $\pi^2$, the top of `NUMERICS.md` §12's range, is refused. A grid is a geometric object
  and a domain is a semantic constraint; keeping the bound in the rule means a later domain change
  does not reshape the grid.
- **Objective (stage 1).** For each candidate, the max over `theta:dense` and all `theta:*` strata
  of max(value error, derivative error), where the derivative is `Dual<S, 1>` through the candidate
  and the reference is `mp.diff`. Minimize it; ties go to fewer terms, then to the larger switch
  point (the series arm is the cheaper one).
- **The second arm (stage 2).** A series arm's term count is set by the largest branch variable it
  serves and paid by the smallest, so **one arm makes every input pay the hardest one's price**
  ([`0039`](./decisions/0039-the-sweeps-grid-stops-below-its-own-optimum.md) item 5) — measured at
  1.87x the catalogue's corpus-weighted term count at binary64 and 1.99x at binary32 under the rule
  below, and at **1.14x to 4.25x** of its latency at near-identity `theta`
  ([`0047`](./decisions/0047-the-second-arm-is-admitted-by-agreement-not-by-the-objective.md)).
  Each `Switch`
  therefore carries a **second**, shorter arm: the first $m_0$ terms of the same series, taken below
  a second switch $z_0$.

  It is searched **after** the objective and never with it, because it cannot change the objective —
  a maximum the long arm already attains wherever it is selected, which a shorter arm underneath can
  only raise (item 6). A joint minimization returns a degenerate answer that reads as "a second arm
  buys nothing". So:

  1. stage 1 is the search above, unchanged;
  2. stage 2 takes stage 1's $(m, z)$ and finds the **cheapest prefix that agrees with the whole arm
     to the bit** — the largest grid point $z_0$ such that at every grid point and every corpus
     record below it the two arms give the same bits, in the value and in $\mathrm{d}/\mathrm{d}z$
     alike ([`0047`](./decisions/0047-the-second-arm-is-admitted-by-agreement-not-by-the-objective.md)
     item 1). "Holds the objective" is implied and is **not** the rule: the objective is one maximum
     and a bar is per stratum, and a prefix admitted by the objective alone was measured raising an
     individual record by up to **50224x** for 6-9% less cost.

  Minimize the corpus-weighted term count at **one term, one unit** over $m_0$, `O(m · n)` over a
  prefix maximum — there is no `grid²` (item 10). A timed cost model would make the generated file a
  function of the host and break D16. A second arm is a switch **and** a shorter prefix or neither;
  "neither" is a second switch of `0`, which no branch variable is below.

  The second arm is `helicoid::coeffs`' alone: `xtask/src/seeded`'s generated file keeps one arm,
  because a seeded twin of a prefix that is bit-identical by construction validates nothing
  (`0047` item 6).
- **Feasibility is reported with the objective, and a boundary optimum is not a result.** The two
  limits above are jointly binding — neither alone gets the worst swept row below 1846 $u$, both
  together reach 94.5 — so an objective quoted without the box it was found in reads as an optimum
  when it may be the edge ([`0039`](./decisions/0039-the-sweeps-grid-stops-below-its-own-optimum.md)
  items 2, 8, 9). Therefore: a choice on the grid's **first** point **fails the run**, because the
  search space is too small and nothing in the result would say so; a choice on the last point
  below $\pi^2$ is the **domain** binding and is printed; a choice at $m = 16$ is the **term cap**
  binding and is printed; and every generated switch carries its grid span, its term cap and which
  limit it sits against.
- **Prior.** `tf_tree` D12 (θ < 0.1, four terms, one switch for $a$, $b$, $c$) is evaluated as a
  named candidate and reported beside the generated choice.
- **Generated file.** `crates/helicoid/src/coeffs/generated.rs`:

  ```rust
  // @generated by `cargo xtask thresholds` from conformance/sweeps/thresholds.csv — do not edit.
  // Source sweep rev: sha256 <of the CSV>.
  // Source series: conformance/corpus/coeff_series.jsonl, sha256 <of the file>.

  // Objective (max u): value 1.9, derivative 3.4.
  // Feasibility: grid [1e-16, 1e1] of 64 per decade searched below π², terms <= 16; binding: none.
  // Second arm: 3 terms, 11735 term-units of corpus against 22087.
  pub(crate) const B_F64: Switch<f64, 5> = Switch::first(
      f64::from_bits(0x3f8d_...), // θ² < 1.44e-2
      f64::from_bits(0x3f67_...), // θ² < 2.84e-3
      3,
      &SWEPT_B_F64,
  );

  pub(crate) const SWEPT_B_F64: [f64; 16] = [1.6666666666666666e-1, /* correctly rounded rationals */ ];
  ```

  Each switch is preceded by its objective **and** its feasibility and second-arm lines, as above.
- **The reference at each switch.** `conformance/corpus/coeff_switch_ref.jsonl` holds the true
  value and $\mathrm{d}/\mathrm{d}z$ of every coefficient at **every grid point**, at 30 digits,
  recomputed at 150. `docs/maths/coefficients.md` CO.12 bounds the jump between a coefficient's two
  arms at its switch by the sum of the arms' errors there, and sampling that right-hand side at the
  two corpus records bracketing the switch understates it by up to 3.5x — the exact arm's error is
  a sawtooth, swinging 195x over 0.8% of $\theta$, so no density of records makes a point sample a
  bound. A reference at every grid point is a reference at every switch the sweep can return and is
  a function of the grid alone, never of the sweep, which keeps `just corpus-check` independent of
  `just thresholds`.

  Switch points are emitted as bit patterns with a decimal comment; series coefficients as
  correctly-rounded literals of the exact rationals, separately per precision (no `f64 → f32`
  double rounding). The eight swept terms sit beside the `Switch`es, not inside them: the sweep
  reads only those, a function of the corpus, so the file is a fixed point of the tool (a second
  run writes the same bytes) and a corpus whose series changed is found stale before sweeping. The
  revision is a SHA-256, not a git revision, which cannot name its own commit.
- **Drift gate.** `just thresholds-check` reruns the sweep and fails on any difference.
- **Phase 1 wiring.** The tool is built and validated against `xtask::seeded` kernels (a correct
  one and the §10 defects). Phase 3 points it at `helicoid::coeffs` through the private `__sweep`
  feature, so **the code measured is the code shipped**.

## 7. Oracle runners

**NORMATIVE.**

- **File protocol.** A runner reads corpus JSONL for the function ids it supports and writes
  `{"id":…,"out":{…}}` with hex-float outputs; the harness computes metrics. Convention conversion
  happens inside the runner, in functions named after the foreign convention
  (`sophus_to_helicoid_tangent`, `gtsam_pose3_to_helicoid`), each with a hand-computed unit test.
- **Runners.** `runners/tf_tree_math` (pinned tf_tree commit): SO(3)/SE(3) `exp`/`log`,
  the action (`Quat::rotate`), the adjoint (`Iso3::adjoint`, on the six basis twists), the
  quaternion of a matrix (`quat_from_rot3`), `slerp` and ScLERP (`so3_geodesic`,
  `se3_geodesic`), $a, b, c$ via its `V`/`V⁻¹`. A runner answers every id its oracle exports and the corpus holds, on every stratum
  of the file: selecting strata would make a paired comparison the runner's choice. `runners/sophus_rs` (pinned version):
  SO(3)/SE(3) `exp`/`log`, any exposed Jacobians. `docker/oracles/sophus` (Sophus C++),
  `manif` (header-only, small C++ runner), `gtsam` (`Rot3`/`Pose3` `Expmap`/`Logmap` and their
  derivatives).
- **Pinning.** Every runner pins its oracle; the envelope records oracle versions.
- **Oracles may be wrong.** Their errors are recorded; a stratum where every oracle is bad does not
  lower the bar — the baseline governs.

## 8. The envelope and the bars

**NORMATIVE.**

- **Domination.** For every `(fn, stratum, precision)` with at least one oracle row:
  `helicoid.max_u ≤ min(oracle.max_u)`. Ties pass.
- **No-regress.** `helicoid.max_u ≤ baseline.max_u`, compared **exactly** (D16 makes outputs
  deterministic). A PR that improves a max updates `conformance/baseline/` in the same PR with
  `just envelope --bless`; the diff is reviewed.
- **Coverage.** Every function id named in `NUMERICS.md` or `API.md` §3 has a corpus file once its
  phase lands; `xtask envelope` fails otherwise.
- **Evidence.** `--bless` regenerates `docs/evidence/ENVELOPE.md` (`@generated`): per function and
  stratum, helicoid's max and p99, the best oracle and its version.

## 9. Bench harness and gate

- `criterion` (pinned) in `crates/helicoid/benches` from Phase 3, fixtures by stratum
  (near-identity, generic, near-π).
- `cargo xtask bench-gate`: interleaved baseline/candidate runs, paired bootstrap 95% CI of the
  ratio; **fail when the whole CI lies above $1 + \delta$**, with $\delta$ the host's A/A noise
  floor measured by `bench-gate --aa` and recorded in `baseline/HOST.md`. A measured floor, not an
  asserted percentage.
- Oracle rows (`tf_tree_math`, sophus-rs) are **reported, never gated**.

## 10. Seeded-defect validation

**NORMATIVE.** `just conformance --self-test` runs every planted defect as a subject and **fails if
any is not detected** by its named mechanism.

| Planted defect (`xtask/src/seeded/`) | Must be detected by |
|---|---|
| `b` by its definition, no series | error curve in `theta:1e-8`…`1e-2` fits $\theta^{-p}$ with $p \in [1.8, 2.2]$ |
| `c` with switch $10^{-8}$ and two terms | sweep ranks it dominated; envelope fails against the correct seeded kernel |
| `Log` via `acos` of the trace | `theta:1e-k` and `theta:pi-1e-k` max ≥ $10^7\,u$ |
| SE(3) `Exp` reading the tangent translation-first | every `rho:*` stratum fails |
| `Q` with $-\tfrac12\rho^\wedge$ | `sen3_jr*` fails; `Dual` comparison fails |
| one `sqrt` of $\theta^2$ above the branch, no safe argument, the series arm at the rebuilt $\theta\cdot\theta$, under `Dual` ([`0020`](./decisions/0020-dual-sqrt-at-zero-keeps-its-derivative.md)) | `nonfinite > 0` in `theta:exact0` |
| `Log` without the $w < 0$ flip | the negated-quaternion half of every `so3_log` stratum fails |

## 11. Definition of done

- [ ] `just build test no-std wasm msrv lint doc` green on x86_64 and aarch64.
- [ ] Corpus v1 committed; `corpus-check` byte-identical on Linux x86_64 and macOS aarch64.
- [ ] `conformance --self-test` detects all §10 defects.
- [ ] `tf_tree_math` and sophus-rs runners produce envelope rows for SO(3)/SE(3) `exp`/`log`; the
      three containers build and produce rows (nightly).
- [ ] The sweep regenerates a seeded `generated.rs` byte-identically; `thresholds-check` fails on a
      hand edit (test).
- [x] `bench-gate --aa` noise floor recorded in `baseline/HOST.md`.
- [ ] `CONTRIBUTING.md`, licenses, `deny.toml`; the dependency-closure lint rejects a planted
      `nalgebra` dependency (test).

## Appendix: suggested implementation order

1. Workspace, lints, empty library crates, `justfile`, CI for `build`/`no-std`/`wasm`/`msrv`.
2. Generator: splitmix64, schema, the coefficient ids and `so3_*` ids, `MANIFEST.json`,
   `corpus-check`.
3. Harness: `Subject`, metrics, CSV, seeded defects for coefficients and `Log`.
4. Sweep against seeded kernels; generated-file format; `thresholds-check`.
5. `sen3_*`, `so2_*`, `se2_*` ids; the remaining seeded defects; `--self-test`.
6. Runners: `tf_tree_math`, sophus-rs; then the three containers.
7. Envelope, bars, `--bless`, `ENVELOPE.md`.
8. Bench gate and the A/A floor.
