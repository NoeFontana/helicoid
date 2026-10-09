# 0062: sin θ/θ is public, and differentiates its branch

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #118, #122

## Context

`α(θ²) = sin θ/θ` is the function a camera or bearing model needs near the optical axis. In an
equidistant or Kannala–Brandt unprojection, the bearing is `(α(θ²) θ_d û, cos θ)`, and S²'s retract
is `cos θ n + α(θ²) (Bδ × n)` (CH.9(a)). A consumer that writes it by hand needs a guard at
`θ = 0`. locus-tag's `camera.rs` already has one, `NEAR_AXIS_RADIUS = 1e-8`, chosen by hand, on
the Kannala–Brandt `θ_d/r`. A guard on sin θ/θ would need a second one, plus a series if the
derivative is to survive below it. The derivative is where it fails. `dα/dz = (θ cos θ − sin θ)/(2θ³)`
cancels as `θ → 0`, so a hand-written Jacobian loses digits exactly where a solver linearizes.
`0004` says switch points are generated, not typed, and that only `helicoid::coeffs` evaluates a
cancelling coefficient. A consumer that needs α today has to break one rule or the other in its
own crate.

`helicoid` has no α. `NUMERICS.md` §4's catalogue has `k`, `a`, `b`, `c`, `d`, `e` and `r`, and
its `a` is `(1 − cos θ)/θ²`, not sin θ/θ. §6 names sin θ/θ `α` for SE(2) ("both 0/0 only: series
from the generator, switch generated"), and `0015` (draft) P1.5 proposes an id `coeff_alpha`.
Neither is implemented. CH.9(a) already writes S²'s α from the two outputs of `exp_coeffs`:
`α = 2k cos(θ/2)`, `k = sin(θ/2)/θ`.

### Measured

**Protocol.** A scratch crate outside the workspace (not committed) depends on `helicoid` with
`__sweep` and calls `exp_coeffs` and `jr_coeffs` at 4 414 values of `z = θ²`: `θ` log-spaced over
`[1e-9, π]` (4 001 points), `π − 10^-k` for `k = 1..12`, and `[3.2, 6.2]` (401 points). It runs
at binary64, and at binary32 with `z` rounded. An mpmath 1.4.1 script at 60 digits gives the
reference `sin θ/θ` and `(θ cos θ − sin θ)/(2θ³)` at the exact `z` of each record. Errors are
relative, in `u` of each precision.

| Candidate | θ ≤ 1e-4 | 1e-4 to 1 | 1 to 2.5 | 2.5 to π |
|---|---|---|---|---|
| value `2k · cos(θ/2)` (`exp_coeffs`) | 1.00 | 1.39 | 2.02 | 1.24e12 |
| value `sin θ / θ`, `θ = √z` (naive) | 1.41 | 1.48 | 1.70 | 1.24e12 |
| derivative `(b − a)/2` (`jr_coeffs`) | 1.60 | 2.08 | 3.60 | 6.13 |
| derivative `Dual<S, 1>` through `2k · cos(θ/2)` | 1.98 | 2.01 | 1.88 | 2.61 |

The near-π column is the function's conditioning, not an algorithm's. The relative condition of
`α` in `z` is `|z α′/α|`, which grows like `1.6/(π − θ)`. Every form that forms `θ` loses it, the
naive one included. Divided by that condition (floored at 1), the worst errors over the whole
range, `θ` up to 6.2, are:

| | value | `Dual` derivative |
|---|---|---|
| binary64, `θ ≤ π` / `θ > π` | 1.87 / 2.29 `u` | 2.44 / 3.63 `u` |
| binary32, `θ ≤ π` / `θ > π` | 2.03 / 2.49 `u` | 2.31 / 4.07 `u` |

The identity `2k cos(θ/2) = sin θ/θ` agrees to `2.3e-61` at 60 digits, as it must:
`sin θ = 2 sin(θ/2) cos(θ/2)`.

So the derivative this record started from, the closed form `(b − a)/2` (exact:
`b − a = (θ cos θ − sin θ)/θ³`, CH.5(b)), loses to the `Dual` branch derivative by up to 2.3×
on `[1, π]`. The decision takes the measured winner.

## Decision

1. **`helicoid::sinc<S: Real>(theta_sq: S) -> (S, S)`** returns `(α, dα/dz)` at `z = θ²`, and
   **`helicoid::sinc_value<S: Real>(theta_sq: S) -> S`** returns `α` alone: decision 2's product on
   plain `S`, bit-identical to `sinc`'s value, for a per-point camera or bearing model that needs
   no derivative and should not pay for decision 3's `Dual` evaluation (two `Dual` Horner chains
   below the switch, a `Dual` `sqrt` and `sin_cos` above it). Both are re-exported at the crate
   root from `coeffs`. `coeffs` stays `pub(crate)`, so only these two functions cross the
   boundary.
2. **The value is `2k · cos(θ/2)`** from `exp_coeffs`, the group CH.9(a) and `Exp` already use. It
   needs no new switch, no new series and no new `generated.rs` row, so `cargo xtask thresholds`
   writes the same bytes.
3. **The derivative is the branch's:** `exp_coeffs` evaluated on `Dual<S, 1>::variable(z, 0)`,
   lane 0 of `2k · cos(θ/2)`. It is the derivative of the arm `exp_coeffs` takes at that `z`, what
   the corpus scores as `d_branch` for every `coeff_*` id, and the same `S: Real` code, so it holds
   for `S = Dual<_, _>` too. One `Dual` evaluation gives both outputs.
4. **`NUMERICS.md`.**
   - §4 gains a row, `α = sin θ/θ`, computed as `2k cos(θ/2)` on both arms of `exp_coeffs`, with
     no switch and no series of its own. It therefore adds nothing to `coeff_series` or the sweep.
     The row also says that α is evaluated only as `sinc`, inside `coeffs`. The rule "a call
     site never evaluates one coefficient on its own" holds because `sinc` is a group of one built
     from a group.
   - §6's α refers to it.
   - §12 gains the row `sinc`, `sinc_value`: defined for every `z ≥ 0`; the stated accuracy, the
     value and `dα/dz` within `5u` times the condition number, is for `θ ≤ π`, the range
     `coeff_alpha`'s strata score and every named consumer (a camera's angle off its axis, SE(2)'s
     `θ`) stays in. The scratch measurement above holds the same bound to `θ = 6.2`, where
     `dα/dz = 0` at `θ = 4.4934`, but no stratum verifies it there (`0006`), so the row does not
     promise it.
5. **Corpus id `coeff_alpha`** joins the `coeff_*` family: `COEFF_STRATA` (the scalar `theta:*`
   strata and their binary32 twins), the reference `sin θ/θ` and `d_branch = (θ cos θ − sin θ)/(2θ³)`
   at the corpus's precision, scored as the other coefficients. `xtask/src/shipped.rs`'s `Swept`
   subject answers it through `sinc`. `coverage.rs`'s `required()` gains it, and `PHASE1.md` §4.3's
   coefficient definitions gain α.

## Rationale

- **Public at all.** This removes a hand-typed switch point and a cancelling derivative that the
  consumer cannot see are wrong (`API.md` §6 item 8, `0049`). It does not relocate a formula
  locus-tag already has right: it has none, and its one guard is on a different function.
- **`2k cos(θ/2)` over a new generated α.** §6's sentence reads as a coefficient with its own
  switch. That would add a row to the sweep, which `0039` shows is still binding on its grid, for
  a value the table shows no better than the existing group's. If SE(2) (`PHASE3.md` §6) wants
  α's own switch, that is a measurement it can owe; this record does not foreclose it.
- **The `Dual` derivative over `(b − a)/2`.** It is measured better on `[1, π]` and equal below.
  It also differentiates the arm actually taken, so `sinc`'s two outputs are consistent with each
  other and with what a `Dual` caller of `S2Chart` would get.
- **One function returning both.** That is what was asked for, and the `Dual` evaluation computes
  both anyway. A value-only caller pays for one `Dual<S, 1>` evaluation, a cost this record has
  not measured. If a bench shows it matters, a value-only `sinc_value` is a later, additive record.

## Consequences

- `coeffs` has its first public item. `API.md` §3's sentence "`coeffs` is `pub(crate)`" gains "except
  `sinc` (`0062`)", and `API.md` §2 lists `sinc`.
- `exp_coeffs`'s switches now also bound a public function's accuracy. A future regeneration that
  moves them is scored on `coeff_alpha` as well, which is the point of the id.

## Implementation plan

1. This record and the `NUMERICS.md` §4, §6, §12 and `API.md` §2, §3 edits. Verified by
   `just lint` and `just doc`.
2. `coeffs::alpha` (the `Dual` evaluation of decision 3) and the re-export `sinc`, with tests in
   `coeffs/tests.rs`: `sinc_is_two_k_cos_half` (bit-identical to the group's product),
   `sinc_derivative_is_the_dual_lane`, and `sinc_at_zero_is_one_and_minus_a_sixth`. Then the
   `coeff_alpha` generator entry, the `Swept` subject arm, `required()` and `PHASE1.md` §4.3.
   Verified by `just test`, `just corpus-check`, `just conformance`, `just envelope`,
   `just thresholds-check` (unchanged bytes), `just no-std` and `just msrv`.

## Open questions

None.

## Further work

- locus-tag's L3 adopts `sinc` in its camera models; its migration gate compares unprojected
  bearings before and after.
- SE(2)'s `se2_coeffs` (`0015` (draft) PH.3) may reuse decision 2's value or sweep its own α. Either
  way, it is measured on `coeff_alpha`.
