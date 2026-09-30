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
| Workspace, lints, `justfile`, CI matrix (§3) | Partial: workspace, lints, `justfile` and CI for `build`/`test`/`lint`/`audit`/`msrv`/`no-std`/`wasm`/`doc`/`corpus-check`; `cargo xtask lint` Partial: line citations, draft-record citations (only `0.0`/`0.` status tables, Rust comments under `crates/` and `xtask/`, amendment banners; not prose, not `PROJECT.md` §5.1), `@generated` header and registry (registry empty, regeneration comparison arrives with each generator); normal-dependency closure of `helicoid-linalg` and `helicoid` against the `0007` set (from `cargo metadata`, all features; `mint` only as an optional direct dependency, never a default; by package name; any other workspace member must register a budget), `__sweep` (no member but `xtask` requests, forwards or defaults it), each with a planted-dependency test; `deny.toml` bans the `0007` list, checked by `just audit` for every crate but `xtask`; owed: twin table; `determinism`, `oracles`, `bench-check` jobs not started |
| Reference generator and committed corpus v1 (§4) | Partial: generator skeleton (per-stratum splitmix64 streams, schema, `MANIFEST.json`, 150-digit recheck, `just corpus`/`corpus-check`) and the scalar-θ strata; ids implemented: `coeff_k`, `coeff_a`…`coeff_e` (1710 records each) and `coeff_r` (1713: the θ strata as unit quaternions, plus `q:w0` at norms 1, 1e-3, 1e3), each record cross-checked at generation, and `coeff_series` (16-term exact rational Taylor series of the seven §4 coefficients; manifest `kind: "series"`, §4.3); readings in `conformance/generate/README.md`. Missing: `so3_*`, `sen3_*`, `so2_*`, `se2_*` (SE(2)'s α, β and the cos θ/2 series with them); the vector `q:w0`, `q:nonunit` and `rho:*` strata; matrix encoding; macOS aarch64 byte-identity unchecked (§11) |
| Conformance harness and result schema (§5) | Not started |
| Threshold sweep and generated-file format (§6) | Not started |
| Oracle runners: `tf_tree_math`, sophus-rs (excluded crates); Sophus, manif, GTSAM (containers) (§7) | Not started |
| Envelope and bars (§8) | Not started |
| Bench harness and gate (§9) | Not started |
| Seeded-defect validation (§10) | Not started |

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
   quaternion power series, `Log` as a matrix logarithm, Jacobians as their defining series,
   inverses by high-precision linear algebra, coefficients from their definitions.
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
`corpus` (the records above; `rechecked` counts those recomputed at 150 digits) or `series`
(`coeff_series`, §4.3; `verified` counts its series equal to exact algebra).

### 4.3 Definitions (v1 function ids)

| Function id | Input | Reference computed as |
|---|---|---|
| `coeff_k`, `coeff_a`, `coeff_b`, `coeff_c`, `coeff_d`, `coeff_e` | $\theta$ | the definition at 120 digits; derivative by `mp.diff` |
| `coeff_r` | $(n, w)$ | $2\,\mathrm{atan2}(n, w)/n$ at 120 digits |
| `coeff_series` | — | the exact series of [`0004`](./decisions/0004-switch-points-are-generated-not-typed.md): `mp.taylor` of each definition at 120 digits, rationalized, equal term by term to an independent exact derivation |
| `so3_exp` | $\varphi$ | quaternion power series $\sum p^n/n!$, $p = (0, \varphi/2)$ |
| `so3_log` | $q$ | `mp.logm` of $R(q/\|q\|)$; stratum `q:w0` uses the sign rule of `NUMERICS.md` §3.2 |
| `so3_act` | $(q, p)$ | $R(q/\|q\|)\,p$ |
| `so3_from_matrix` | $R$ | `mp.logm`-consistent quaternion; backward error only (`NUMERICS.md` §11) |
| `so3_jr`, `so3_jl` | $\varphi$ | $\sum (\mp W)^n/(n+1)!$ |
| `so3_jr_inv`, `so3_jl_inv` | $\varphi$ | `mp.inverse` of the above |
| `sen3_exp_n{1,2,3}` | $\tau$ | `mp.expm` of the $(3+N)$-square hat matrix |
| `sen3_log_n{1,2,3}` | $X$ | `mp.logm` |
| `sen3_ad_n{1,2,3}` | $X$ | images of the basis under $\sigma \mapsto (X\sigma^\wedge X^{-1})^\vee$ |
| `sen3_jr_n{1,2,3}`, `sen3_jl_n{1,2,3}` | $\tau$ | $\sum (\mp\,\mathrm{ad}_\tau)^n/(n+1)!$ |
| `sen3_jr_inv_n{1,2,3}`, `sen3_jl_inv_n{1,2,3}` | $\tau$ | `mp.inverse` of the above |
| `so2_*`, `se2_*` | analogous | analogous |

A `coeff_k`…`coeff_e` or `coeff_r` record's outputs are `value` and `d_branch`: the derivative with
respect to the branch variable of `NUMERICS.md` §4, by `mp.diff` at the exact real branch value of
the exact binary64 input (θ², not `fl(θ·θ)`).

`coeff_series` is the one file that does not follow §4.2: one record per coefficient,
`{"id","coeff","branch","prefactor","series"}` with sixteen `"num/den"` terms in the branch
variable, and no `stratum`, `in` or `out`. Its manifest `kind` is `series`; the conformance
harness, the envelope and the oracle runners read `kind: "corpus"` files only.

Series are summed until the term's norm is below $10^{-110}$ relative. Phase 2 adds `eig3`
(`mp.eigsy`), `svd3` (`mp.svd_r`), `solve_cubic` (`mp.polyroots`); Phases 4–5 add their ids by the
spec that needs them.

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
- `q:nonunit`: quaternions at $\|q\|^2 - 1 = \pm 2^{-45}$ (`Log`'s scale invariance, `from_*` normalization).

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
- **Grid.** Series terms $m \in \{1, \dots, 8\}$; switch point in the branch variable on a log grid
  of 64 points per decade, spanning $\theta \in [10^{-8}, 1]$.
- **Objective.** For each candidate, the max over `theta:dense` and all `theta:*` strata of
  max(value error, derivative error), where the derivative is `Dual<S, 1>` through the candidate and
  the reference is `mp.diff`. Minimize it; ties go to fewer terms, then to the larger switch point
  (the series arm is the cheaper one).
- **Prior.** `tf_tree` D12 (θ < 0.1, four terms, one switch for $a$, $b$, $c$) is evaluated as a
  named candidate and reported beside the generated choice.
- **Generated file.** `crates/helicoid/src/coeffs/generated.rs`:

  ```rust
  // @generated by `cargo xtask thresholds` from conformance/sweeps/thresholds.csv — do not edit.
  // Source sweep rev: <git rev>. Objective (max u): value 1.9, derivative 3.4.
  pub(crate) const B_F64: Switch<f64, 5> = Switch {
      below: f64::from_bits(0x3f8d_...), // θ² < 1.44e-2
      series: [1.0 / 6.0, /* correctly rounded rationals */ ],
  };
  ```

  Switch points are emitted as bit patterns with a decimal comment; series coefficients as
  correctly-rounded literals of the exact rationals, separately per precision (no `f64 → f32`
  double rounding).
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
  $a, b, c$ via its `V`/`V⁻¹`, `slerp`, ScLERP (Phase 4). `runners/sophus_rs` (pinned version):
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
- [ ] `bench-gate --aa` noise floor recorded in `baseline/HOST.md`.
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
