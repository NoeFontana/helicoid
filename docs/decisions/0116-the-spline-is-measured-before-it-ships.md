# 0116: The spline is measured before it ships

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #229, #230

## Context

[`0006`](./0006-the-instrument-comes-first.md) puts the instrument before the algorithm: no
`helicoid-spline` kernel ([`0114`](./0114-a-spline-evaluates-the-window-it-is-handed.md)) lands
before its corpus ids, references and seeded defects exist.

Three facts set this record's shape.

**The corpus has 4.00 MiB left.** The tracked files under `conformance/corpus/` total 79 694 370
bytes at `f6f6a44` (`git ls-tree -r -l HEAD conformance/corpus`), against
[`0110`](./0110-the-corpus-budget-is-80-mib-and-checked.md)'s 83 886 080.

At `0110`'s sizes (22 bytes per hex input, 38 per decimal output, about 100 per record), with
6 records per stratum and `K` averaged over 4 to 6, the ids below cost:

| Layout | binary64 | both precisions |
|---|---|---|
| every axis crossed: 5 increments × 4 `u` × 3 `K` × 3 `h` = 180 strata | 14.7 MiB | 29.3 MiB |
| one axis at a time, one default cell: 5 + 3 + 2 = 10 strata (Rⁿ: 8) | 0.8 MiB | 1.6 MiB |
| the same without SE(3) | 0.6 MiB | 1.2 MiB |

The crossed layout of the first sketch does not fit.

**No envelope is blessed.** `docs/evidence/ENVELOPE.md` holds no candidate rows, and `PHASE3.md`
§0.0 has the first blessing pending. The domination and no-regress bars apply to spline ids from the
workspace's first blessing (`PHASE3.md` §10), as they do to every other id.

**A strata cell under 10 records is judged on the deep tier**
([`0088`](./0088-the-domination-bar-must-tell-worse-from-equal.md)), which is git-ignored and has no
budget. SE_N(3) already runs at 6 records per committed stratum.

## Decision

1. **Function ids** (`PHASE1.md` §4.3 rows):

   | Id | Input | Output |
   |---|---|---|
   | `so3_spline_value`, `_vel`, `_acc` | `K` quaternions, `u`, `inv_h` | `X`; `ω h`; `ω̇ h²` |
   | `so3_spline_value_jac`, `_vel_jac`, `_acc_jac` | the same | `K` right Jacobians, dense, each output scaled as above |
   | `rn_spline_value`, `_vel`, `_acc` | `K` points of R³, `u`, `inv_h` | `p`; `ṗ h`; `p̈ h²` |
   | `rn_spline_value_jac`, `_vel_jac`, `_acc_jac` | the same | the `K` scalars `k` of each `RnJac = k I`, scaled as above |
   | `se3_spline_value`, `_vel`, `_acc` | `K` poses, `u`, `inv_h` | as SO(3), twists |
   | `rn_spline_gram` | `K`, `r`, `inv_h` | `G⁽ʳ⁾`, 0-ulp bar |

   - Scaling by `h` and `h²` makes every derivative output unit-scale.
   - Values use `NUMERICS.md` §11's forward error (quaternion sign-aligned). Tangents and Jacobians
     use the vector 2-norm and the Frobenius norm, floor 1.
   - Products are not ids. `0114`'s bitwise product-equals-components proptest covers them.
2. **Strata, one axis at a time.** The default cell is `inc:generic`, `u:random`, `h:1`; it is
   stratum `inc:generic` and is drawn once. Every other stratum moves one axis off it:
   - increments `‖d_j‖`:
     - `inc:tiny`: at most `10⁻⁸`, the series arms;
     - `inc:generic`: log-uniform over `10⁻³` to 1;
     - `inc:large`: 2.5 to `π − 10⁻³`;
     - `inc:near-pi`: within `10⁻⁶` of π;
     - `inc:mixed`: one tiny and one large increment per window;
   - `u:0`, `u:half`, `u:1-ulp` (the largest `u` below 1 at the stratum's precision: `1 − 2⁻⁵³`,
     and `1 − 2⁻²⁴` in the `@f32` twin, which is drawn there because rounding `1 − 2⁻⁵³` to
     binary32 gives 1);
   - `h:1e-3`, `h:1e3`.

   The `rn_spline_*` ids take `inc:generic`, `inc:large` and `inc:mixed` only: Rⁿ has no series arm
   and no π.

   Records cycle through the id's orders: `K = 4, 5, 6` for the groups, `2` to `6` for Rⁿ. Committed: 6 records per stratum, with an `@f32` twin
   ([`0113`](./0113-phase-7-opens-for-a-calibration-consumer.md)). Deep tier: 64.
   `rn_spline_gram` is one record per `(K, r, h)` with `h ∈ {2⁻¹⁰, 1, 2¹⁰}`, so the scaling is
   exact ([`0115`](./0115-the-spline-tables-are-exact-integers.md)).

   Budget: about 1.7 MiB at both precisions, with SE(3) and the Rⁿ Jacobian ids. The rest stays for later
   phases.
3. **References** (generator, independent of the recursions):
   - **Value**: the cumulative product with mpmath's exponential. `d_j` is the geometric logarithm,
     never `mp.logm` (as for S², Sim(3) and [`0045`](./0045-two-phase-4-checks-cannot-be-taken-as-written.md)).
   - **Derivatives**: `ω = (X⁻¹ dX/dt)^∨` and its derivative, by `mp.diff` of the reference curve at
     the generator's working precision.
   - **Jacobians**: `mp.diff` along `X_m Exp(ε e_k)` for every control `m` and axis `k`.
   - **Gram**: exact rationals, rounded once per precision.
4. **Oracle: basalt-headers** (Sommer and Usenko's implementation of these recursions). It is an
   image on [`0072`](./0072-a-container-oracle-is-an-image-on-the-same-protocol.md)'s protocol, and
   independent of `helicoid` ([`0107`](./0107-the-consumer-is-not-an-oracle.md)). The spike settles
   four things, each recorded per id before the runner merges:
   1. the orders it compiles;
   2. which Jacobians it exposes;
   3. whether its SE(3) spline is cumulative SE(3) or split. If it is split, it answers no
      `se3_spline_*` id;
   4. whether it accepts `u` directly or only integer nanoseconds. If only nanoseconds, it answers
      only records whose `(u, h)` map to integer nanoseconds exactly, so that time quantization is
      not scored as its error.

   Kalibr's splines use another parameterization and are not an oracle (`0113`).

   Where no oracle answers (binary32 and unexposed Jacobians), an id is held to its reference and
   its baseline, and `ENVELOPE.md` says so.
5. **Seeded defects** (`PHASE1.md` §10 rows). The clean max is the defect-free subject's on the same
   stratum.

   | Planted defect | Must be detected by |
   |---|---|
   | the non-cumulative basis in `λ` | every `so3_spline_value` and `rn_spline_value` stratum but `inc:tiny@f32` reads more than `1024×` the clean max |
   | `J_r⁻¹` for `J_l⁻¹` in `∂d_j/∂δ_{j−1}` | every `so3_spline_value_jac` stratum but `inc:tiny@f32` reads more than `1024×` the clean max |
   | the `ad` term of the `ω̇` recursion with its sign flipped | every `so3_spline_acc` stratum but `inc:tiny` and `inc:mixed` reads more than `1024×` the clean max |
   | `Ad_{A_j}` for `Ad_{A_j⁻¹}` | every `so3_spline_vel` stratum but `inc:tiny` and `inc:mixed` reads more than `1024×` the clean max |
   | `λ̈` without its `h⁻²` | `h:1e-3` and `h:1e3` of `so3_spline_acc` and `rn_spline_acc` read more than `1024×` the clean max; every other stratum (all at `h = 1`) reads the clean max to the bit |
   | a Gram matrix without its `h^{1−2r}` | `rn_spline_gram` fails its 0-ulp bar at `h ≠ 1` and passes at `h = 1` |

   The exceptions are where the defect's size, about `‖d‖` or `‖d‖²`, falls under the precision's
   `u`. The self-test confirms them against measured clean maxima. A row that does not hold as
   written amends this record, and is never quietly widened.
6. **Bench rows** (criterion, replicate-outermost like the group benches):
   - `spline/{so3, rn3, so3xr3, se3}/{4, 5, 6}/{value, kin1, kin2}/{plain, jac}`;
   - `spline/{so3, so3xr3}/{4, 5, 6}/segment_new` and `jac_segment_new`;
   - `spline/so3xr3/6/batch64/{kin2, kin2_jac}`, 64 samples on one segment.

   The gate is no-regress against the first committed baseline. No absolute number is set (`0113`).
7. **Bars.** Per `(id, stratum, precision)` on the max, as `PHASE1.md` §8 states. A real loss blocks
   the spline's envelope step, or gets a `ready` record and a row in
   [`0046`](./0046-explained-by-record-needs-a-record-to-point-at.md)'s exceptions table. If
   `inc:near-pi` or `inc:mixed` loses, that record says whether the loss is `Log`'s own or the
   spline's.

## Rationale

- **Crossing the axes** would need 29 MiB. The interactions it would catch are covered by the
  proptests, which sample all axes jointly at no corpus cost, and by the deep tier.
- **Raising the cap again** (as `0110` did) can come later. A spline that fits in half the headroom
  should not spend the rest of the phases' share.
- **One id per output** keeps each bar attributable. A single `so3_spline` id with three outputs
  would let a velocity regression hide behind a value maximum.

## Consequences

- `coverage::required()` lists every new id. `corpus-check` and `corpus-test` cover them as they
  cover the rest.
- `just determinism` digests them by construction ([`0109`](./0109-determinism-is-required-before-phase-6.md)).
- The basalt-headers runner is a local `just oracle-basalt` recipe and never a CI job.

## Implementation plan

1. Generator ids, strata and references. Verified by `just corpus-check` and `just corpus-test`, and
   by the budget lint.
2. The six seeded defects. Verified by `just conformance --self-test` failing on each, with
   the subject stubbed to `0114` step 1's reference twins until `0114` steps 2 and 3 land.
3. The basalt-headers image and runner, after the spike. Verified by `just oracle-basalt` answering
   the ids the spike lists.
4. Bench rows and the first baseline, after `0114`'s kernels. Verified by `just bench --against`.

## Open questions

None.
- If basalt-headers takes integer nanoseconds only, its rows stand on the integer-nanosecond subset
  (decision 4.4).
- The SE(3) Jacobians and every `::<Left>` Jacobian are not ids. They run the code that
  `so3_spline_*_jac` scores, through `SEn3Jac` and `Ad`, which are scored on their own. Each is
  held to its reference twin and to `Dual` by proptest. An SE(3) Jacobian record is about 8 KB,
  and three ids at both precisions would take 2.8 MB of the 4.00 MiB left.
