---
paths: ["runners/**", "docker/**", "conformance/**", "xtask/**"]
---

# Oracles, corpus and envelope

Oracle runners are workspace-excluded and never a normal dependency of a library crate. Containers,
noise bars and long debug sweeps stay local `just` recipes, never CI jobs.

| Recipe | Purpose |
|---|---|
| `just oracle-sophus-rs` | the `runners/sophus_rs` runner: its fmt, clippy, doc, deny (own `deny.toml`), tests, then its rows over the corpus (`conformance --oracle`) |
| `just oracle-nalgebra` | the same for `runners/nalgebra` (own `deny.toml`): `eig3`, `chol_*`, `solve_cubic`, `quat_renormalize` (`0056`), `svd3`, `nearest_rotation` (`0101`; `0102`) |
| `just oracle-sophus` | Sophus C++ container (`docker/oracles/sophus`, needs docker): builds the pinned image, runs its rows offline (`0072`); `eig3` and `svd3` are the image's Eigen (`0083`, `0101`), `nearest_rotation` Sophus's `fitToSO3` (`0102`) |
| `just oracle-gtsam` | GTSAM C++ (`docker/oracles/gtsam`), built in the image under the runner's flags (`0074`) |
| `just oracle-manif` | manif container: SO(2), SE(2), SO(3), SE(3), SE₂(3), screw and world charts (`0075`) |
| `just oracles` | every `oracle-*` above, in order; local only |
| `just corpus-deep` / `corpus-deep-check` / `conformance-deep` | `0088`'s 64-record deep tier, git-ignored, manifest committed; feeds `just envelope --deep` |
| `just envelope` | merge results, apply the domination and no-regress bars, then coverage; each domination failure says whether it is the program's, D16's or neither (`0037`, draft). `--check` also fails on a baseline or `docs/evidence/ENVELOPE.md` that `--bless` would change; `--bless` (for `helicoid`) writes both, and nothing while a bar or coverage fails. Never `--bless` unasked |
| `just compare` | `--a X --b Y [--fn ID] [--precision f64\|f32]`: paired comparison of two subjects, per stratum, per id, pooled (`0038`, `0078`). Printed only, **not a bar**. Choose between two programs with this, never by counting domination failures |
| `just bench` | criterion; `--against <bench-binary>` is the gate, `--aa` logs this host's noise (`PHASE1.md` §9). `bench-check` is **owed** |

- `conformance/corpus*`, `conformance/results*`, `baseline/bench` are data: never read or hand-edit
  them (`just corpus` regenerates; `corpus-check` verifies).
- Separate slow measurement from pure logic: replay recorded data, in the producer's own format,
  in the test suite; every expensive finding leaves a cheap regression test behind.
- A bench gate run is only trustworthy on a quiet host: check `ps` for other CPU load, and trust its
  concurrent A/A floor over single-run numbers.
