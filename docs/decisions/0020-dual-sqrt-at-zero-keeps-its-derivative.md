# 0020: `Dual::sqrt` at zero keeps its derivative

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** #35 (draft, authorises nothing). Two tests in `dual_tests` state today's
behaviour of `Vector<Dual>::norm` and `chol` and change none.

## Context

`Dual::sqrt` applies `d / (2 sqrt v)` in every lane. At `v = 0` a lane with `d != 0` is `+-inf` and
a lane with `d = 0` is NaN (`0 / 0`); the value is unaffected. `PHASE2.md` §3, `NUMERICS.md` §12,
the `# Domain` of `Dual::sqrt` and `Vector::norm`, and `EA.18` (`docs/maths/error-analysis.md`) all
read it that way. The question is whether that should change. This record recommends not, and
makes the recommendation conditional on Q1 and Q3.

**Who reaches `sqrt` at 0 on a `Dual`** (`crates/helicoid` is still empty, so the group rows are
the spec, not code):

| Site | Reaches `sqrt(0)`? |
|---|---|
| `Vector<Dual>::norm` of the zero vector | Yes: `d(n²) = 2 v·dv = 0` in every lane, so every lane is `0/0`. No consumer in the workspace; pinned by `dual_norm_of_the_zero_vector_is_nan_in_every_lane`. |
| Quaternion renormalise | No: the Newton step has no `sqrt` (`NUMERICS.md` §3.6). |
| `SO3::log` at θ = 0 | Not in the exact arm: `r(n², w)` is a `helicoid::coeffs` coefficient, whose exact arm sees `select(small, 1, n²)` (0003 item 3, `EA.18(c)`). |
| Shepperd `from_matrix` | The selected arm, no: its argument is the largest of four quantities that sum to 4, so at least 1 (`NUMERICS.md` §3.4). The other three arms are exactly 0 at the identity: a `bool` mask never evaluates them; a lane mask does, and only the value-part `select` keeps their NaN out. Not pinned. |
| `chol` pivot on `Dual` | No: `select(ok, d, 1).sqrt()`, so a zero, negative or NaN pivot reaches `sqrt` as the constant 1. Pinned by `a_zero_cholesky_pivot_leaves_every_dual_derivative_lane_finite` and `a_negative_cholesky_pivot_leaves_every_dual_derivative_lane_finite`. |

The `eig3`, `svd3` and `solve_cubic` ports (`PHASE2.md` §6) are not written; each must pass the safe
argument or say why not. So every routine written so far passes it, and the exposed surface is a
**consumer's own** `Vector<Dual>::norm` at an exactly-zero vector (a residual norm at a converged
point is the obvious case). Whether one exists is Q1.

**What the NaN is worth.** A `sqrt` hoisted out of a `branch` and shared with the selected arm, at
`φ = 0`, gives `θ.d = [NaN; 3]` and NaN in `sin(θ/2).d` and `cos(θ/2).d` (`EA.18(b)`). Under a
`d`-masking rule the same code gives all-zero derivatives: finite, the derivative of the analytic
factor `k(θ²)`, and for the non-differentiable `sin(θ/2)` alone a chosen subgradient. The
`PHASE1.md` §10 row ("`sqrt` of θ² without the safe argument, under `Dual`", `nonfinite > 0` in
`theta:exact0`) fires **only** in that hoisted shape. With the `sqrt` written inside the exact arm
and no safe argument, a `bool` mask evaluates that arm alone at `φ ≠ 0` and `Dual::select` blends
per lane, so the unselected NaN is dropped and every output is finite under any rule. The row does
not state the shape, and `xtask/src/seeded/` does not exist yet, so the row is not evidence for or
against a rule until Q3 is answered.

**A fact that shapes the options.** If the argument `z(x)` is `C¹`, non-negative in a neighbourhood
and `z(x0) = 0`, then `x0` is an interior minimum and `∇z(x0) = 0`. So a lane with `d != 0` at
`v = 0` exists only at a domain boundary, where `+inf` is the true one-sided slope (`sqrt(x)` at
`x = 0`); the `d = 0` lanes are the `|x|`-like interior case, where no derivative exists and `0`
is a member of the Clarke subdifferential (of the norm at 0: the unit ball, whose smallest element
is 0). Underflow breaks the premise (`z = x²` rounds to 0 with `x = 1e-200` and `d = 2e-200`, true
slope 1); that input is outside `norm`'s stated `# Domain` (`norm_sq` must be normal).

## Decision

**Keep the rule**, provided Q1 is answered "no consumer" and Q3 does not make NaN the only
detector. `Dual::sqrt` stays `d / (2 sqrt v)`; no `Real` method, no `Vector` method, no mask on
`v == 0` or `d == 0`. The one-line policy: *the derivative of the shipped code is reported as
computed, and a non-finite derivative at a non-differentiable point is the report.*

The follow-up (step 1 below) is documentation only: cite this record from the `PHASE2.md` §3
bullet once it is `ready`, and give `Vector::norm` a doctest of the guarded form
(`select(n2 <= 0, 1, n2).sqrt()`, then `select(n2 <= 0, 0, ·)`), the one existing path
(`API.md` §6 item 5).

## Rationale

Candidates, by what each does at `v = 0`:

| | rule | `d != 0` | `d = 0` | reads | verdict |
|---|---|---|---|---|---|
| (a) | as is | `+-inf` (true at a boundary) | NaN | value | **chosen** |
| (b) | per lane `d == 0 → 0` (ForwardDiff's opt-in NaN-safe mode) | `+-inf` | 0 | derivative lanes | rejected |
| (e) | `v == 0 →` all lanes 0 | 0 (wrong at a boundary) | 0 | value | rejected |
| (s) | `v == 0 →` lane `0 * d` (Stan `fvar`) | 0 (same) | 0; NaN or `inf` in `d` survives | value | rejected |
| (c) | leave `sqrt`; add a zero-safe norm | `sqrt` unchanged | norm derivative 0 | value | not now |

- **(e) and (s) are wrong at a real boundary inside the domain.** `sqrt(x)` at `x = 0` has slope
  `+inf` from the right; they return 0 there. (a) and (b) keep it.
- **(b) breaks a `Dual` invariant.** Masks and `select` read the value part only (`dual.rs`
  header; `0003` item 5). (b) makes a derivative lane steer a select. It keeps the value path
  bitwise (`r` is untouched) but the derivative lanes of `-0` become `+0`.
- **No candidate is chosen on the canary.** Hoisted, (b), (e) and (s) all blind the §10 row and (a)
  keeps it; unhoisted, no rule fires it. Which shape the seeded defect has is Q3. If the row must
  fire and the shape is the hoisted one, (a) is required; if the row does not need NaN, (b) is the
  candidate to re-argue, not (e) or (s).
- **(c) is the right shape for the real hazard and is not owed yet.** It is bit-equal to `norm`
  for `n² > 0` and cheap (Evidence), but `API.md` §6 item 5 asks whether it duplicates an existing
  path: `norm_sq` behind the safe argument is that path, and `0009` says the consumer owns its own
  conventions. Its form is fixed now so that a request has a place:
  `select(z, 0, select(z, 1, n²).sqrt())`, `z = n² <= 0`, built from existing `Real` operations.
  It shares the underflow hazard: for `x = 1e-200` as variable 0, `n²` rounds to 0, plain `norm`
  gives `d = [inf, NaN, NaN]` and the guarded form the silent `d = [0, 0, 0]` against a true
  `(1, 0, 0)`. That input is outside the stated domain, and a doctest must say so.
- **The strongest case for (b), and why it does not win here.** ForwardDiff's own argument: a lane
  whose perturbation is exactly 0 "should not propagate derivative information", yet `0 * inf`
  poisons it with NaN, which is what (a) does to every lane of a zero-vector norm. That is a real
  cost to a consumer. It loses here because every routine written so far routes its own `sqrt`
  through the safe argument, so the cost lands only on an external caller, and (b) is the one
  candidate that makes a derivative lane steer a select.
- **Other AD systems** (illustrative; Evidence, Sources). Ceres `Jet` `sqrt` is (a); its `hypot`
  says the function "is non-smooth at x=y=0, so the derivative is undefined there". Stan's `fvar`
  `sqrt` is (s). ForwardDiff's default is (a) and its opt-in `nansafe_mode` is (b), documented at
  "~5%-10%" cost. PyTorch pairs (a) for `sqrt` (`grad / (2 * result)`) with (c) for the 2-norm
  backward (`masked_fill_(norm == 0, 0)`). JAX's FAQ prescribes the double-`where`, the same as
  the safe argument.

## Consequences

- `PHASE2.md` §3, `NUMERICS.md` §12, `EA.18` and the `PHASE1.md` §10 row stand as written, the
  row with the Q3 caveat.
- A consumer's `Vector<Dual>::norm` at an exactly-zero vector returns NaN derivatives. Documented,
  pinned, not fixed. A Gauss-Newton residual norm at an exactly converged point would poison a
  whole Jacobian; that is why Q1 gates the decision.
- The value path claim ("bitwise identical to plain `S`") is untouched by every candidate.
- If Q1 is answered yes, a `Vector::norm_safe`-style surface needs its own record (or an amendment
  of this one) and the `API.md` §6 answers.

## Implementation plan

1. Docs only, once `ready`: cite 0020 from the `PHASE2.md` §3 `sqrt`-at-0 bullet and from
   `NUMERICS.md` §12's `Dual::sqrt` row; add the guarded-norm doctest to `Vector::norm` — verified
   by `just doc`, `just test`, `just lint`.

Reversing to (b), (e) or (s) would touch: the `Dual::sqrt` body, its `# Domain` and the header
table; `dual_sqrt_at_zero_has_an_infinite_derivative_and_the_plain_value` and
`dual_norm_of_the_zero_vector_is_nan_in_every_lane`; the `PHASE2.md` §3 bullet and §0.0 row; the
`NUMERICS.md` §12 row; `EA.18(a)`, `(b)` and its *Checked* block; `chol`'s `NUMERICS.md` §15.5 last
sentence; `Vector::norm`'s `# Domain`; and the `PHASE1.md` §10 row. (c) would add one method, one
corpus stratum or reference twin (`NUMERICS.md` §14), a `# Domain`, and the §6 checklist.

## Open questions

1. **Is there a consumer whose `Vector<Dual>::norm` meets an exactly-zero vector?** `locus-*` and
   `tf_tree` were not inspected. A yes moves (c) from "not now" to a proposal and puts (b) back on
   the table. Must be answered "no" before this record goes `ready`.
2. **Which subgradient, if (c) ships?** 0 (the smallest Clarke element; what PyTorch uses) versus a
   unit vector. `Dual::abs` picks `sgn(+-0) = +1`, not 0; for a 1-vector the two conventions would
   disagree. The reason `abs` chose `+1` is not recorded here.
3. **What shape does the §10 seeded defect have, and is NaN its only detector?** The row does not
   say (Context). A derivative-versus-mpmath comparison at `theta:exact0` would also catch it under
   (b), (e), (s). Must be answered before the canary can count as evidence for either side.

## Evidence

Illustrative, not reproducible from the repository: the harness is a scratch crate (path
dependency on `crates/helicoid-linalg`, `std`, release, `lto = "fat"`, `codegen-units = 1`,
`libm` 0.2 with `default-features = false`, and once with `arch`) that is not committed, per the
rule that benchmark harnesses stay out unless the docs prescribe one. What the repository does
reproduce is the two pin tests and `dual_sqrt_at_zero_has_an_infinite_derivative_and_the_plain_value`.

Protocol: one pinned core of a shared 8-core host (load average 5 to 13 from sibling jobs). 1024
`Dual<f64, N>` inputs, `v` uniform in `[0.1, 10.1]`, `d` uniform in `±0.5`; also with half of the
`d` lanes exactly 0 and with 10% of `v = 0`. Variants are generic and inlined; each sample runs all
variants once so contention is shared; p10 of 1500 samples, ns per element. Controls: the crate's
`Dual::sqrt` and a verbatim copy agree to 0.1 ns. The variants must use `Real::sqrt`; the inherent
`f64::sqrt` (hardware) gave a bogus 2x. Absolute numbers are this host's; only deltas carry over,
and a delta under about 3 ns has an unstable sign. Two full runs of the `arch`-off configuration
agreed on every delta within about 1 ns.

`libm` software `sqrt` (`arch` off, the shipped configuration), ns/op, current / (b) / (e) / (s),
first run:

| | N = 3 | N = 6 |
|---|---|---|
| sqrt throughput | 7.9 / 8.2 / 8.1 / 8.1 | 11.9 / 14.8 / 12.6 / 13.2 |
| sqrt latency (dependent chain) | 21.3 / 21.3 / 21.4 / 21.3 | 21.2 / 21.2 / 21.2 / 21.3 |
| `Vec3<Dual>` `norm_sq` + `sqrt` | 55.1 / 55.9 / 56.0 / 53.5 | 67.0 / 70.4 / 66.8 / 67.7 |

The second run read 10.9 / 12.2 / 11.5 / 12.0 and 12.5 / 15.5 / 13.1 / 13.7 (throughput), and
61.4 / 65.6 / 61.5 / 62.1 for the `N = 6` norm. With `arch` on (one run) the throughput rows read
4.3 / 4.7 / 4.3 / 4.3 and 5.4 / 7.6 / 5.5 / 5.4 (+41% for (b) at `N = 6`); norm 34.8 / 36.0 / 35.7 /
35.2 and 45.2 / 49.7 / 45.4 / 46.4. Latency at `N = 6` reads 13.1 / 8.5 / 13.0 / 8.0: codegen, not
the mask ((b) and (s) are faster there). Summary: (e) and (s) are within about 1.3 ns of the
current rule; (b) is +3 ns (+24%) on an isolated `N = 6` sqrt in both runs and +4 ns (+7%) on the
norm, just above the noise floor. Cost is small in every variant and does not decide the choice.

Bit checks: over 200000 inputs with `v > 0`, (b), (e), (s) equal the current rule in `v` and every
`d` lane bit for bit; Ceres' `d * (1 / 2r)` differs in 83336 of them and is not a candidate (D16).

Claims run in the same crate: `norm` of the zero vector with three variables gives `d = [NaN; 3]`;
the double-`where` norm gives `d = [0; 3]` and is bit-identical to `norm` at `(3, 4, 0)`; its cost
against `norm` (p10, ns) is `f64` 19.0 to 19.9, `Dual<f64, 3>` 56.1 to 58.2, `Dual<f64, 6>` 62.8
to 62.4 (software `sqrt`); hoisted θ at `φ = 0`: current rule all NaN, (b) all 0.

Sources, read on 2026-09-30 at these commits: `ceres-solver/ceres-solver@e17a9b4`
`include/ceres/jet.h` (`sqrt`, `hypot`); `stan-dev/math@8f96bdf`
`stan/math/fwd/fun/sqrt.hpp`; `pytorch/pytorch@5357b32`
`torch/csrc/autograd/FunctionsManual.cpp` (`norm_backward`) and `tools/autograd/derivatives.yaml`
(`sqrt`); `JuliaDiff/ForwardDiff.jl@a3c0f4f` `src/partials.jl` (`_mul_partial`, `_div_partial`)
and `docs/src/user/advanced.md` (NaN-safe mode); `jax-ml/jax@23e25ac` `docs/faq.rst`
("Gradients contain `NaN` where using `where`").
