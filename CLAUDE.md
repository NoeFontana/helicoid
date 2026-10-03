# CLAUDE.md — agent guidance for helicoid

For AI agents. Humans: `CONTRIBUTING.md` (owed by [`PHASE1.md`](./docs/PHASE1.md) §3). `helicoid`
is the one Lie-group layer under `tf_tree`, omnisac, locus-tag, locus-calib, locus_fusion and
fuse-geometry, built along the seven-phase roadmap in [`docs/PROJECT.md`](./docs/PROJECT.md) §4.

**The documents in [`docs/`](./docs/) are the contract.** Read `docs/PROJECT.md`; then
`docs/NUMERICS.md` if you touch a formula; then the spec for the phase you touch — **before**
touching code. When they do not answer a question, **stop and ask**. Never invent a formula, a
series coefficient, a switch point, a tangent ordering or a sign.

## Status

Each spec's **§0.0 status table is the source of truth**, over its own prose, the README and this
file. Nothing is published. The workspace is on `0.0.x` (every release may break every other)
until [`PHASE6.md`](./docs/PHASE6.md) §6's 1.0 criteria hold. MSRV **1.87**, and never above the
lowest MSRV among consumers (D17).

## Read before changing X

| Doc | Role |
|---|---|
| [`PROJECT.md`](./docs/PROJECT.md) | Overview, roadmap, decision log D1–D18 (§5), gated-work ledger (§5.1), design smells (§6). |
| [`NUMERICS.md`](./docs/NUMERICS.md) | Every formula, stated rotation-first. Coefficient catalogue (§4), domains (§12), reference-twin table (§14). |
| [`API.md`](./docs/API.md) | Six rules (§1), crate surfaces (§2–§4), stability (§5), the new-surface checklist (§6). |
| [`PHASE1.md`](./docs/PHASE1.md) | Workspace, CI, reference generator, conformance harness, threshold sweep, oracle runners, bench gate. **No group code.** |
| [`PHASE2.md`](./docs/PHASE2.md) | `helicoid-linalg`: `Real`, `Mask`, `Dual`, fixed-size types, strided views, `eig3`/`svd3`/`solve_cubic`. |
| [`PHASE3.md`](./docs/PHASE3.md) | `helicoid`: coefficient kernel, SO(2), SO(3), SE(2), SE_N(3), Rⁿ, products, sides, Jacobians. |
| [`PHASE4.md`](./docs/PHASE4.md) | Geodesics and their Jacobians; what `helicoid` owes the `tf_tree` migration. |
| [`PHASE5.md`](./docs/PHASE5.md) | Charts, S², Sim(3), Γ functions, `Gaussian`; the omnisac / locus-tag / locus_fusion migrations. |
| [`PHASE6.md`](./docs/PHASE6.md) | Ambient Jacobians, `mint`, cross-target bit identity, locus-calib, the 1.0 criteria. |
| [`PHASE7.md`](./docs/PHASE7.md) | `helicoid-spline`. **GATED ([`0011`](./docs/decisions/0011-continuous-time-waits-for-a-consumer.md)), not scheduled**; §3 is questions, not answers. |
| [`decisions/`](./docs/decisions/) | A record's own `**Status:**` line says whether it is authoritative: `grep -m1 -H '^\*\*Status:' docs/decisions/0*.md`. [`README.md`](./docs/decisions/README.md) indexes them. [`0013`](./docs/decisions/0013-simd-lanes-owe-a-measurement.md) is a **draft** and authorises nothing. |

## Hard rules — do not relitigate

- **Conventions** ([`0002`](./docs/decisions/0002-one-convention-for-a-stack-that-already-disagrees.md)):
  Hamilton, **`w` first**, active; `a * b` is `T_a_x · T_x_b`; tangents are **rotation-first**
  `[φ; ρ₁; …; ρ_N]` (`[ω; v]` for a twist); right perturbation is the default and **every
  perturbation names its side** (`rplus`/`lplus`/`rminus`/`lminus`). No `Add`/`Sub` impl for ⊕/⊖.
- **Scalar model** ([`0003`](./docs/decisions/0003-the-scalar-that-cannot-say-less-than.md)):
  numeric code is generic over `S: Real`. `Real` has **no `PartialOrd`, on purpose** — generic
  code cannot write `if theta < eps`. Comparisons return `S::Mask`; branches go through
  `S::branch`/`S::select`, with the safe-argument pattern in the exact arm.
- **Measured numerics** ([`0004`](./docs/decisions/0004-switch-points-are-generated-not-typed.md)):
  `crates/helicoid/src/coeffs/generated.rs` is written by `cargo xtask thresholds` from a committed
  sweep; a hand edit fails `just thresholds-check`. Do not type a series coefficient or a switch
  point. **Only `helicoid::coeffs` evaluates a cancelling coefficient** (`(θ − sin θ)/θ³` and its
  relatives in `NUMERICS.md` §4).
- **`Log` goes through the quaternion `atan2`** (D5). No `acos` of a trace anywhere in the workspace.
- **Instrument before algorithm** ([`0006`](./docs/decisions/0006-the-instrument-comes-first.md)):
  a routine without a corpus stratum is unverified and does not ship. The bars are **domination**
  over the oracle envelope and **no-regress** against the committed baseline, per function, per
  stratum, per precision, on the **max** — never a mean.
- **Reference plus fast, forever** (D6): every row of `NUMERICS.md` §14 keeps its reference twin
  and its `*_matches_reference` proptest. Do not delete a twin because the fast version "is
  obviously right".
- **Jacobians are structured** ([`0005`](./docs/decisions/0005-the-jacobian-is-a-dual-matrix.md)):
  `G::Jac`; dense only through `write_dense`/`sandwich`. No group method returns `[[S; 6]; 6]`.
- **Budget** ([`0007`](./docs/decisions/0007-the-budget-a-foundation-can-afford.md)):
  `helicoid-linalg` = `libm`; `helicoid` = `helicoid-linalg` + `libm`; optional `mint`; nothing
  else. `no_std`, **no `alloc`**, `#![forbid(unsafe_code)]` on every library root. `xtask` and the
  oracle runners are unrestricted and are **never** a normal dependency of a library crate.
- **Values, not storage** (D2): no `Pod`, `Zeroable`, `serde`. Layout is not a semver contract;
  consumers own their storage formats (tf_tree's arena records stay tf_tree's).
- **Scope** ([`0009`](./docs/decisions/0009-what-helicoid-does-not-own.md)): no camera/IMU/LiDAR
  model (that is `locus-*`), no solver, no preintegration, no spline
  ([`0011`](./docs/decisions/0011-continuous-time-waits-for-a-consumer.md)), no SIMD `Real` impl
  ([`0013`](./docs/decisions/0013-simd-lanes-owe-a-measurement.md)).
- **A retraction is a chart** ([`0012`](./docs/decisions/0012-a-retraction-is-a-chart.md)): a
  consumer that needs a different retraction gets a `Chart` here, never an inline formula in its
  own crate.
- **Determinism** (D16): every transcendental through the `libm` crate; no `mul_add`, no
  `target-cpu`, no fast-math. Outputs are bit-identical on x86_64, aarch64 and wasm32.
- Lints: `unwrap_used`, `expect_used`, `panic`, `todo`, `unimplemented`, `dbg_macro` are denied.
  No panic on valid input; a stated domain is a `debug_assert!`, never a release check (D11).
- **Cite a symbol, never a line number**; `just lint` fails a new line citation. **A comment states
  the decision and its load-bearing evidence, then stops**; provenance belongs in the owning
  record, cited by number. Prefer deleting prose to adding it.

## Commands

Everything goes through `just`; CI invokes the recipes 1:1. Tests: **nextest, not `cargo test`**.
Single test: `cargo nextest run -p helicoid -- so3_exp_log_roundtrip`. The generator is Python
(`uv`), the container oracles are C++; everything else is Rust.

| Recipe | Purpose |
|---|---|
| `just build` / `test` | build; nextest (dev and release profiles) + doctests |
| `just lint` | fmt, clippy `-D warnings`, `cargo xtask lint` (line citations, draft citations, generated-file headers, dependency closure, `__sweep`; twin table owed), `cargo deny` |
| `just doc` | rustdoc, warnings denied |
| `just msrv` | build at `rust-version` |
| `just no-std` | build for `thumbv7em-none-eabihf` (no `std`, no `alloc`) |
| `just wasm` | build and run the conformance subject under `wasm32-wasip1` (wasmtime) |
| `just corpus` / `corpus-check` | regenerate the mpmath corpus; check it is byte-identical to the committed one |
| `just conformance` | run in-process subjects over the corpus → `conformance/results/`; `--self-test` runs the seeded defects |
| `just thresholds` / `thresholds-check` | run the sweep, regenerate `coeffs/generated.rs`; check for drift |
| `just oracle-tf-tree-math` | the `tf_tree_math` runner (workspace-excluded): its fmt, clippy, doc, deny, tests, then its rows over the corpus (`conformance --oracle`) |
| `just oracle-sophus-rs` | the same for the sophus-rs runner (`runners/sophus_rs`, audited with its own `deny.toml`) |
| `just oracles` | **owed** (`PHASE1.md` §0.0): container-only Sophus, manif, GTSAM runners plus the two excluded-crate runners. The two that exist are the `oracle-*` recipes above |
| `just envelope` | merge results, apply the domination and no-regress bars, then coverage; `--check` also fails on a baseline or `docs/evidence/ENVELOPE.md` that `--bless` would change; `--bless` (for `helicoid`) writes both, and nothing while a bar or coverage fails |
| `just determinism` | **owed** (`PHASE1.md` §0.0): compare output digests across x86_64, aarch64, wasm32 |
| `just bench` | criterion; `--against <bench-binary>` is the gate, `--aa` logs this host's noise (`PHASE1.md` §9). `bench-check` is **owed** |

## Decision workflow

A change the specs do not cover — a new public item, a new crate, a new dependency, a convention,
a new stratum family, a chart — starts as a **`draft` decision record** in
[`docs/decisions/`](./docs/decisions/), **not** as a PR (see its `README.md`). When a record is
`ready`, implement it as stated; its *Implementation plan* is the per-PR breakdown. A formula change
is a `NUMERICS.md` edit **and** a record.
