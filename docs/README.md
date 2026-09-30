# Docs

Normative specifications for whoever is changing `helicoid`, not tutorials.

| If you want to | Read |
|---|---|
| Call the API | [`API.md`](./API.md) §2–§4, after its six rules in §1 |
| Know which formula a routine implements, and its domain | [`NUMERICS.md`](./NUMERICS.md) §3–§12 |
| Know why a formula holds, or re-derive a sign or a side | [`maths/index.md`](./maths/index.md) (non-normative; `NUMERICS.md` wins), [`maths/lie-groups.md`](./maths/lie-groups.md) for $\mathrm{Ad}$, $\mathrm{ad}$, $J_r$, $J_l$ and every row of `NUMERICS.md` §2.3, [`maths/so3.md`](./maths/so3.md) for the SO(3) quaternion, `Exp`, `Log` (and why not $\arccos$), the closed forms of $J$ and $J^{-1}$, `from_matrix` and renormalization, [`maths/coefficients.md`](./maths/coefficients.md) for the series, rounding error, switch-point magnitudes and derivatives of the seven coefficients of `NUMERICS.md` §4 |
| Know whether a number is measured, and against what | [`NUMERICS.md`](./NUMERICS.md) §11, [`decisions/0006`](./decisions/0006-the-instrument-comes-first.md), `evidence/ENVELOPE.md` (generated from Phase 3) |
| Know whether a feature exists yet | The `§0.0` table heading each `PHASEn.md` |
| Know why it works the way it does | [`PROJECT.md`](./PROJECT.md) §5, the decision log |
| Know why something was *not* built | [`decisions/0009`](./decisions/0009-what-helicoid-does-not-own.md), [`PHASE7.md`](./PHASE7.md) §0.0, [`PROJECT.md`](./PROJECT.md) §5.1 |

## Reading order for changing the project

1. [`PROJECT.md`](./PROJECT.md): what, why, the seven-phase roadmap, **decision log D1–D18** (§5).
2. [`NUMERICS.md`](./NUMERICS.md): conventions (§1), the coefficient catalogue (§4), every
   formula, the error metrics (§11), domains (§12). Stated rotation-first; never permute it in
   your head.
3. [`PHASE1.md`](./PHASE1.md): the instrument. Nothing in Phases 2–7 is verified without it.
4. [`PHASE2.md`](./PHASE2.md): `helicoid-linalg`, the scalar model.
5. [`PHASE3.md`](./PHASE3.md): the groups.
6. [`PHASE4.md`](./PHASE4.md): geodesics; the `tf_tree` migration.
7. [`PHASE5.md`](./PHASE5.md): charts, S², Sim(3), Γ, `Gaussian`; the retraction migrations.
8. [`PHASE6.md`](./PHASE6.md): interop, determinism, locus-calib, 1.0.

Cross-cutting:

- [`API.md`](./API.md): the six rules (§1), the surfaces, the §6 checklist a new surface passes.
- [`PHASE7.md`](./PHASE7.md): `helicoid-spline`, **gated by [`0011`](./decisions/0011-continuous-time-waits-for-a-consumer.md) and not scheduled**.
- [`decisions/`](./decisions/): records; [`README.md`](./decisions/README.md) has the lifecycle.
