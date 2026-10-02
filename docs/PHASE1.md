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
| Workspace, lints, `justfile`, CI matrix (§3) | Partial: workspace, lints, `justfile` and CI for `build`/`test`/`lint`/`audit`/`msrv`/`no-std`/`wasm`/`doc`/`corpus-check`/`conformance --self-test`/`thresholds-check`; `cargo xtask lint` Partial: line citations, draft-record citations (only `0.0`/`0.` status tables, Rust comments under `crates/` and `xtask/`, amendment banners; not prose, not `PROJECT.md` §5.1), `@generated` header and registry (`xtask/src/seeded/generated.rs` by `thresholds`; the regeneration comparison is each generator's own `--check`); normal-dependency closure of `helicoid-linalg` and `helicoid` against the `0007` set (from `cargo metadata`, all features; `mint` only as an optional direct dependency, never a default; by package name; any other workspace member must register a budget), `__sweep` (no member but `xtask` requests, forwards or defaults it), each with a planted-dependency test; `deny.toml` bans the `0007` list, checked by `just audit` for every crate but `xtask`; owed: twin table; the `wasm` job builds and does not run the conformance subject under wasmtime (owed to Phase 6's determinism job); `determinism`, `oracles`, `bench-check` jobs not started |
| Reference generator and committed corpus v1 (§4) | Partial: generator skeleton (per-stratum splitmix64 streams, schema, `MANIFEST.json`, 150-digit recheck, `just corpus`/`corpus-check`, parallel by stratum) and the scalar-θ strata; ids implemented: `coeff_k`, `coeff_a`…`coeff_e` and `coeff_cos_half` (cos θ/2 of `NUMERICS.md` §3.1; 3420 records each: the θ strata and their `@f32` twins) and `coeff_r` (3426: the θ strata as unit quaternions, plus `q:w0` at norms 1, 1e-3, 1e3, and their `@f32` twins), each record cross-checked at generation, and `coeff_series` (16-term exact rational Taylor series of the seven §4 coefficients and, as an eighth row, `cos_half`; manifest `kind: "series"`, §4.3); each `@f32` stratum (`0016` item 1) is the binary64 inputs rounded to nearest-even binary32 (`r`: n and w each) with the reference recomputed at the rounded input, every input asserted to round-trip through binary32, `theta:subnormal@f32` a binary32-subnormal decade and `theta:pi-1e-8`…`pi-1e-12` one input five times, none deduplicated; and `so3_exp`, `so3_jr`, `so3_jl`, `so3_jr_inv`, `so3_jl_inv` (2466 records each: the θ strata with an axis per sample), `so3_log` (5124: every quaternion and its negative), `so3_act` and `so3_from_matrix` (2594 each: plus the vector `q:w0` and `q:nonunit`), matrices column-major with a sibling `shape`, each record checked to 100 digits by `mp.expm`, a sandwich product, polar-factor uniqueness or a Jacobian identity; and `sen3_{exp,log,ad,jr,jl,jr_inv,jl_inv}_n{1,2,3}` (312 records per file; `ad` 324 and `log` 642, which also see `q:w0` and `q:nonunit`: the `theta:*` strata at unit translation scale except `theta:dense`, and the 25 `rho:*` × θ cells, 6 records each, not 64 (§4.4); dense matrices column-major, rotation-first), each record checked to 100 digits by `mp.expm`, the conjugation identity, J_l = Ad_Exp(τ) J_r or J J⁻¹ = I, and the dual-matrix structure asserted on the data; and `so2_{exp,log}` (3419 records each: the θ strata, each θ in both signs) and `se2_{exp,log,ad,jr,jl,jr_inv,jl_inv}` (416 records per file: SE_N(3)'s 52 strata, 8 records each, the sign of θ alternating; 3×3 matrices column-major, rotation-first), each record checked to 100 digits of the size of its terms by `mp.expm`, the complex series, the conjugation identity, J_l = Ad_Exp(τ) J_r or J J⁻¹ = I, and the first-row and rotation-block structure asserted on the data; corpus 35 MB of the 50 MB budget, `just corpus` 18 CPU-minutes; readings in `conformance/generate/README.md`. Missing: the coefficient ids and series of SE(2)'s α, β (0015 (draft) P1.5 recommends them with `cos_half`, the one added, on `NUMERICS.md` §3.1 and `0016`'s naming it); `so2_log` at z = (−1, +0), θ = π, unsampled, and at z = (−1, −0), where `atan2` and (−π, π] disagree, undecided; a non-unit or `c = 0` z for `so2_log`, `se2_log`, `se2_ad`, which need a stratum family and so a record (§4.4); macOS aarch64 byte-identity unchecked (§11) |
| Conformance harness and result schema (§5) | Partial: `cargo xtask conformance [--subject NAME] [--fn ID] [--precision f64\|f32]` (`just conformance`): the corpus reader (hex-float inputs, 30-digit decimal references as exact integers, `shape`; a wrong `id` or record count, an unlisted file, an unknown `kind` and a repeated key are errors), `Subject` (`name`, `supports`, `eval`) with an in-process registry (the seeded subjects below), the forward error of `NUMERICS.md` §11 in units of `u` formed exactly in integers and rounded once (norm-wise, sign-aligned quaternions, `so3_log` sign-invariant at `w = +0`, an explicit floor per output field, each row of the table pinned by a test), per `(fn, stratum, precision, subject)` `n`, `max_u`, `p99_u` (nearest rank), `argmax_id` and `nonfinite` (an error that overflows binary64 counts as non-finite), the result CSV under `conformance/results/`, the table by `max_u` descending, a non-zero exit on any non-finite output and when nothing was scored; `f32` (`--precision f32`, `0016` item 2) scores the `@f32` strata alone (a binary64 run skips them) in units of `u = 2^-24`, with the smallest-normal floor `2^-126` (`docs/maths/error-analysis.md` EA.4(b)), asserts every input and every output exactly binary32 (an input is a lossless cast, an output that is not a binary32 an error), makes a subject asked for `f32` on an id with no `@f32` stratum an error, and so a subject with no `f32` kernel (the planted `c`, below), and writes `<subject>[--<fn>]--f32.csv` (`precision` `f32`); without `--fn` it runs the ids that have `@f32` strata, the eight coefficient ids, and names on stderr the others that a selected subject supports (`so3_*`, `sen3_*` for `seeded:correct`), left out and not an error; `so3_from_matrix` is checked for shape and finiteness only (its rows have `max_u`, `p99_u` = NaN); checked over the committed corpus by a subject that returns the correctly rounded reference (at most 1 `u` on every scalar-output id) and one that returns the neighbouring ulp, and once against an independent computation (exact rationals and mpmath 1.4.1 written from §11 alone: 4540 perturbed outputs over all 44 forward ids, every score bit-identical to the correctly rounded exact value; script not committed); at `f32` the same two subjects, correctly rounded once to binary32 and one binary32 ulp above, score at most 1 `u` and in (0.5, 3] on every `@f32` stratum of the eight coefficient ids (225 rows), and 225 `f32` scores (one record per `@f32` stratum of each id; the seeded kernel's outputs, the neighbouring ulp for the ids it does not run) were recomputed once from the raw definitions (sympy limits and derivatives, mpmath at doubling precision, exact rationals; script not committed): every score equal to the harness's, the corpus's 30-digit references within 4.4e-30 of the definitions'. In-process subjects: `seeded:correct`, the kernels `k, a, b, c, d, e` of `NUMERICS.md` §4, `cos θ/2` (§3.1) and `r` (§3.2, on the sweep's reading below) generic over `Real`, run at the generated switches of `xtask/src/seeded/generated.rs` (per coefficient its series length, switch and terms, the sweep's choice over this corpus and grid, not a general optimum; §6), over the eight `coeff_*` ids the correct subject answers (a planted defect runs the six of `k`…`e`), value and `d/dz` of one `Dual<f64, 1>` (`r`: seeded on `n²` at fixed `w`); the kernels also run at the `tf_tree` D12 candidate, in tests and in the sweep, and at `f32` on the `@f32` strata, through the same adapter on `Dual<f32, 1>` (inputs cast exactly, `z = fl(θ·θ)` at `f32`), at the `_F32` constants of `generated.rs` (`0016` item 3; those rows' `subject_version` `generated`); the series rounded once in integers at the precision they run at, read from `coeff_series.jsonl`, `Exp` and `Log` of `NUMERICS.md` §3.1 and §3.2 over `so3_exp` and `so3_log` (`xtask/src/seeded/so3.rs`: `k` at the generated switch, `cos θ/2` by its exact arm at every `z` (its `coeff_series` row is committed and not read), `Log` with the flip as `copysign` and `r = 2 atan2(n, w)/n` at the safe argument, its series arm `2/w` taken where `n² = 0`, no constant typed, the generated `COS_HALF_F64` and `R_F64` not read by `Exp` and `Log`: 0014 (draft) questions 19 and 20), whose largest `max_u` is 2.743 (`so3_exp`, `theta:pi-1e-9`) and 2.901 (`so3_log`, `theta:dense`), no non-finite output, every stratum at most 4, `Exp`, `J_r` and `J_l` of `NUMERICS.md` §5.1 and §5.3 over `sen3_{exp,jr,jl}_n{1,2,3}` (`xtask/src/seeded/se3.rs`: `Exp` is `so3::exp` and `J_l(φ)ρ_i`, `J_r(τ)` is `J_l(−τ)`, `Q` is §5.3's nine words as 3×3 products associated left to right, each coefficient `k, a, b, d, e` at its own generated switch because the sweep has none per call-site group, the dense `J` column-major with exact zeros off the dual-matrix pattern: 0014 (draft) questions 24 and 25), whose largest `max_u` is 3.468 (`sen3_exp_n1`, `rho:1e3/theta=pi-1e-6`) and 8.903 (`sen3_jl_n1`, `rho:1e4/theta=pi-1e-6`), no non-finite output, every stratum at most 12; these are the errors of the per-coefficient kernel, not of the grouped one `NUMERICS.md` §4 prescribes and `helicoid` will run (question 24). Checked once, by a one-off script that is not committed (so not reproducible, and not a gate), against an independent computation (mpmath 60 digits, `ad_τ` from the commutator definition, `Σ(−ad)ⁿ/(n+1)!`; 50 `sen3_jr` outputs over `n = 1, 2, 3`: at most 2.89 `u`, every score above 0.01 `u` equal to the harness's to 1.7e-14 relative); `helicoid` (`xtask/src/shipped.rs`, a plain subject at both precisions) over the eight `coeff_*` ids: the shipped groups of `helicoid::coeffs` (`exp_coeffs`, `jr_coeffs`, `jr_inv_coeff`, `q_coeffs`, `log_ratio`) through `__sweep` on `Dual<S, 1>`, the seeded subject's arguments and fields (`value`, `d_branch`), `subject_version` the workspace version; it equals `seeded:correct` to the bit on every record, so on all 450 `(fn, stratum, precision)` rows (a test: the two share every candidate), no non-finite output, largest `max_u` over the strata, value and `d_branch`, at `f64` `k` 13.6, `a` 25.8, `b` 63.5, `c` 7.40e3, `d` 831, `e` 8.69e3, `cos θ/2` 1.95, `r` 260, at `f32` `k` 12.1, `a` 12.7, `b` 59.2, `c` 1.61e3, `d` 1.04e3, `e` 4.43e3, `cos θ/2` 1.87, `r` 26.7; the prior (four terms below `z = 0.01`: D12 for `a`, `b`, `c`; D12 applied to a coefficient it was not defined for for `k`, `d`, `e`, `cos θ/2`, 0014 (draft) question 10; `(4 terms, s < 0.01)`, `θ = 0.2` and not D12's, for `r`, question 29; true D12 for `r`, `s < 2.5e-3`, is no better than `helicoid` on any stratum either) has the smaller maximum on nine strata of `cos θ/2` and no other (`theta:1e-3`…`1e-8`; at `f32` `theta:1e-2`, `theta:1e-3`, `theta:dense`: 1.4 to 1.9 `u` against 0.85 to 1.7), the sweep's objective being one maximum over the strata, which `theta:1e0` sets and the exact arm from `θ = 7.5e-8` does not exceed (pinned by a test; whether the objective should be per stratum is open, 0014 (draft) question 30; no other stratum of any id has a prior below `helicoid`); the subject's tests are corpus-bound and do not cover the boundary (`s = switch`, `w = 0`) or the safe-argument pattern, which the `helicoid` crate's lane and boundary tests do; it is scored by `just conformance` and writes `conformance/results/helicoid[--f32].csv`, and has no envelope row, no baseline: no oracle yet (§7, §8); and the planted `seeded:b-no-series`, `seeded:k-sqrt-unsafe`, `seeded:c-two-terms-1e-8`, `seeded:log-acos`, `seeded:log-no-flip` (the last two over `so3_log` only), `seeded:se3-exp-translation-first` (over `sen3_exp_n*`) and `seeded:q-minus-half` (over `sen3_jr_n*` and `sen3_jl_n*`), which a run evaluates only when it names them (`--subject`) or under `--self-test`; a rule can score one output field (`COEFF_VALUE`, `COEFF_D_BRANCH`). Missing: backward error (`Log` near π, `from_matrix`), `f32` for the vector ids (they have no `@f32` stratum until a record extends `0016`), the `helicoid` subject on any id but the eight coefficient ones (it has no group code yet), a seeded SE_N(3) `Log`, `Ad`, `J_r⁻¹` and `J_l⁻¹`, the manifest's `sha256` (`just corpus-check` compares bytes); the readings §11 leaves open (floor of a tangent and of a coefficient, the SE_N(3) tangent at `w = +0`, one score per record, one norm over a whole tangent, overflow as non-finite) are in `xtask/src/conformance/`, proposed by 0014 (draft), not in `NUMERICS.md` |
| Threshold sweep and generated-file format (§6) | Partial: `cargo xtask thresholds [--check]` (`just thresholds`, `thresholds-check`, a CI job) over the seeded kernels, at `f64` and at `f32` (`0016` item 3), for `k, a, b, c, d, e`, `cos θ/2` and `r`: 1 to 8 terms times 1025 switch points (64 per decade of the branch variable `z = θ²`, `θ` from `1e-8` to `1`, ends included, each the correctly rounded `10^(i/64)` at the precision it is swept at, from the exact 64th root, so a binary32 point is never a binary64 one rounded again (`grid`'s `rounded` is tested where a second rounding differs, and `sweep_at` refuses a grid point its precision does not hold; on this corpus and grid none differs, so the outputs cannot show it); §6's "per decade" is read as of the branch variable, and the other reading, the even indices, moves `a`, `b`, `e`: 0014 (draft) question 8); the exact and every series arm's error is formed once per record (exact, in `u` of the precision, value and `d/dz` through `Dual<S, 1>`, all `theta:*` records of `coeff_<c>` including `theta:dense`, at `f32` the `@f32` strata), a candidate selects per record by `z < switch`, the objective is the larger maximum, ties go to fewer terms then the larger switch, and the `tf_tree` D12 prior is scored as a candidate beside it (equal, bit for bit, to the whole kernel run at that candidate: a test, at both precisions for all eight ids, at the prior, three grid switches and three at a record's own branch variable); `r` is swept on a reading pending the maintainer, not a spec (`xtask/src/thresholds`, `xtask/src/seeded/kernel.rs`): the branch variable is `s = n²/w²`, formed by a division (not the division-free mask 0015 (draft) NU.6 writes; 0014 (draft) question 29 holds the reading), the grid is over `s`, only a record with `w > 0` enters the objective (every `theta:*` record does, `q:w0` does not), a candidate takes the series arm iff `w > 0` and `s < switch` (`w = 0`, `s` infinite, and `w < 0`, which only S² charts produce, take the exact arm; the corpus has no `w < 0` stratum, so that clause is a unit test's and not measured), `Dual` is seeded on `n²` at fixed `w` (the derivative the corpus stores), and its prior is `(4 terms, s < 0.01)`, the candidate of the other ids taken in `s`, which is `θ = 0.2` where D12's `θ < 0.1` is `s < 2.5e-3`: not D12's, so `r`'s prior figures below are that candidate's; `conformance/sweeps/thresholds-seeded.csv` (16 rows, `f64` first; 25 columns in `xtask/src/thresholds`, the last four the two arms' errors at the two records that bracket the switch, which the jump between them is compared with (0015 (draft) NU.4; a sample at two records, not a bound over the interval)), byte-identical, and a test compares it with a fresh run; and `xtask/src/seeded/generated.rs`, written beside it in the format of §6 from that CSV and `coeff_series.jsonl` (`xtask/src/thresholds/emit.rs`): the `@generated` header, the source revision (the SHA-256 of the CSV's bytes and, on a line of its own, of `coeff_series.jsonl`'s, which supplies the literals; not a git revision, which cannot name its own commit: 0014 (draft) question 12), per coefficient and precision an objective line (question 13) and a `Switch<S, terms>` (`K_F64` … `R_F64`, then `K_F32` … `R_F32`) with the switch as a bit pattern (16 hex digits at `f64`, 8 at `f32`) and its decimal in a comment and the series as the shortest decimal that reads back as each exact rational rounded once at the constant's own precision (a test: every binary32 literal is the nearest binary32 to its rational, in integers; another, on a series term where a second rounding differs, that it is not a binary64 literal rounded again), laid out as `rustfmt` leaves it (a test runs `rustfmt` over the boundaries of its array width, its short-element packing and its line width); registered with `cargo xtask lint`; the seeded correct kernel runs both precisions' constants, and over the corpus scores each row's objective to the bit (a test, at both). The `helicoid` target (`cargo xtask thresholds helicoid`, both targets when none is named) sweeps the arms `helicoid::coeffs` ships, `exact_*` and `series_*` of the hidden feature `__sweep` (`docs/PHASE3.md` §3), and writes `conformance/sweeps/thresholds.csv` and `crates/helicoid/src/coeffs/generated.rs` (registered): a `Switch` per coefficient and precision taking its terms from the eight swept series that follow it in the file, a function of the corpus alone and all the sweep reads, so a second run writes the same bytes; a corpus whose series changed is found stale before sweeping, its series written beside placeholder switches by a run that then fails and is run again, and `--check` reports it. Its CSV is byte-identical to the seeded one (a test: the shipped arms are the seeded arms bit for bit, on every record where the series arm is defined: `w > 0` for `r`, so not the six `q:w0` records). `thresholds --check` (`just thresholds-check`, the CI job) fails on any difference in either file and writes nothing; a test edits one character of each and the check names the file and the line, and another drives the entry point, `run_at`, under a scratch root: it writes both files, passes them, and fails a hand edit without writing. The file is compiled into `xtask` (the seeded correct kernel runs it), so this holds for an edit that still compiles: one that does not stops the tool building, and `git restore` of the file is the recovery. Measured, of this corpus and grid, not general optima (each objective is a maximum over the corpus's records, not over the domain: between records an exact arm above a switch can exceed it), at `f64`: `terms` is 8, §6's cap, for every row but `e`, whose 7 ties with 8; `k` and `d` at the grid top (`θ = 1`); `b` one step below it, a unique optimum 1.3% under the top, decided by the two records between the top two grid points; `e` one step below it, tied with the top, its objective the exact arm's derivative error at the one record at `θ = 1` exactly, which the strict `z < switch` never puts on the series arm (a grid top of `nextUp(1)` gives `e` 8 terms and an objective 3.3 times lower: 0014 (draft) question 9); `a` at `θ = 0.85` and `c` at `0.75`, each sharp (one grid step costs `a` at least 46% and `c` at least 24%), `c`'s at a point where `theta:dense` meets the grid, so a record's arm there is the last bit of `θ·θ`; for `k, b, d` the objective is the exact arm's error over the records with `θ > 1`, which no switch on the grid changes; the D12 prior scores 8.7e3, 2.2e5, 1.2e6, 2.8e7, 1.2e7 and 9.2e9 `u` for `k`…`e`; `cos θ/2`: 2 terms at `θ = 7.5e-8`, objective 1.95 `u`, the exact arm's floor (attained at `theta:1e0`), 5426 candidates tied (D12 prior 2.6e4), so the switch is the tie-break's (only `z = 0`, where the exact arm's derivative is `0/0`, needs the series, and two terms reach the floor up to there); `r`: 8 terms at `s = 1.07e-2` (`n/w = 0.104`, grid index 898), objective 260 `u` (`theta:dense` at `s = 1.7e-2`, on the exact arm; a grid step costs at least 4%), against its prior's 1.1e10. At `f32`: `k, a, b` 5 terms and `c, d, e` 4, all six at the grid top (`θ = 1`, as CO.10 predicts, the optimum lying above the grid), the objective the exact arm's derivative error at a record with `z ≥ 1` (`c`'s at `z = 1` exactly), 12.1, 12.7, 59, 1.6e3, 1.0e3, 4.4e3 `u` (D12 prior 1.8e3, 1.8e3, 1.1e6, 1.9e7, 1.7e7, 8.1e9; one step below the top costs `e` 2.7 times); `cos θ/2`: 2 terms at `z = 2.6e-6`, objective 1.87 `u`, the floor again, tied with 6661 candidates and with the D12 prior; `r`: 8 terms at `s = 0.143` (`n/w = 0.379`, index 970), objective 26.7 `u`, its prior 229. No row is at the bottom of the grid. The CSV names the record that attains the objective (`argmax_*`) and the objective one grid step either side (`below_objective`, `above_objective`); the sweep's other readings, the prior for every id and ties, are 0014 (draft) questions 10 and 11. Missing: `r`'s branch variable, mask and `w ≤ 0` domain as a spec (the reading above is pending the maintainer; 0015 (draft) NU.6 recommends a `w < 0` stratum for `coeff_r`, so the `w > 0` clause has one), SE(2)'s `α`, `β` (no corpus id), and a switch shared by a call-site group (0015 (draft) NU.5) |
| Oracle runners: `tf_tree_math`, sophus-rs (excluded crates); Sophus, manif, GTSAM (containers) (§7) | Not started |
| Envelope and bars (§8) | Not started |
| Bench harness and gate (§9) | Not started |
| Seeded-defect validation (§10) | Partial: `just conformance --self-test` (a CI job) runs the correct seeded kernel, on which no mechanism may fire, and three planted coefficient defects (the two `Log` defects are below), each detected by its named mechanism or the run fails; the other mechanism's reading is printed, not gated (`b` by its definition is `0/0` at `z = 0` and also fires `nonfinite`). `b` by its definition: `coeff_b`'s value `max_u` over `theta:1e-8`…`theta:1e-2` (seven strata, the saturated first included) fits θ^-p, p ∈ [1.8, 2.2], log-log least squares; p = 1.937 on the corpus, −0.10 for the correct kernel. The field is a reading §10 leaves open (0014 (draft), question 7): the value alone, since a record's maximum over `value` and `d_branch` fits p = 3.98, the derivative through the exact arm erring by θ^-4. `k` with an unsafe `sqrt` of θ² under `Dual` (`nonfinite > 0` in `theta:exact0`): the shared root of `sqrt(z)` before the branch, the series arm at `θ·θ`. That is not §10's exact arm without the safe argument, which no `bool` subject can show, a `bool` mask never evaluating the unselected arm; a unit test covers it with a mask that evaluates both arms and counts non-finite operations, on the correct kernel (none, in value and derivative, at `z = 0`, subnormal and either side of the switch) and on the exact arm at `z` (seen). Also tested: each mechanism silent on the correct kernel and firing on its defect, the window, the seven strata and `p` pinned, the self-test failing when a defect is not planted or the window moves. The correct kernel is the generated one (§6), whose per-field maxima the self-test prints and a test equates with the CSV's; at `f32` (below) it runs the `_F32` constants; D12 stays the sweep's prior, its errors pinned by a test. `c` with switch `1e-8` (read as `θ`; `z = 1e-16`, the grid's first point: 0014 (draft) question 14) and two terms is the subject `seeded:c-two-terms-1e-8`, the generated kernels with `c` planted, detected when the sweep of `coeff_c` scores its candidate's objective more than `10^6` times the chosen candidate's (the margin is a reading §10 leaves open: `DOMINATED_BY`, 0014 (draft) question 15; measured 7.3e33 u against 7.4e3 u, rank 7841 of 8200), every subject being measured with its own `c`, so the correct kernel, whose `c` is the chosen candidate, is silent (rank 1); the ranked candidate is the one the subject declares, and the self-test fails when the objective the harness measures for the subject's `c` is not the ranked one, to the bit, so a subject that changes what it computes cannot keep its declaration (a test each way). §10's other half, the envelope failing against the correct seeded kernel, is not started. At `f32` (the `@f32` strata, units of `2^-24`) the two mechanisms that do not need a sweep run, over the same subjects, and the correct kernel there runs the `f32` constants of `generated.rs` (`0016` item 3): its errors reach 4.4e3 `u` (`coeff_e`'s `d_branch` at `theta:1e0@f32`; 1.6e3 for `c`), printed and equal to the CSV's `f32` rows (a test), not gated; each defect is detected by its own. `seeded:b-no-series`: the value of `coeff_b` over `theta:1e-3@f32` to `theta:1e-1@f32`, the strata wholly above `√(6u) = 6e-4`, fits θ^-p with p = 2.134 (r² = 1.000; the correct kernel −0.021, r² = 0.920), against 1.937 at binary64 (−0.10), the window [1.8, 2.2] at both. The `f32` range is a reading (0014 (draft) question 26), chosen after the fit and asymmetric with binary64's: the stratum `theta:1e-4@f32`, which holds `6e-4`, is a plateau of `2^24 u` and with the three above it would fit p = 1.644, outside the window, so the `f32` range omits it, where at binary64 the plateau stratum is one of seven and kept; the top stratum is above D12's switch, where the defect and the D12 kernel the reading was made against ran the same exact arm (297 `u`), and the two below it alone fit 2.188; the `f32` sweep leaves the range as it was. `seeded:k-sqrt-unsafe`: `nonfinite > 0` in `theta:exact0@f32` (`coeff_k`, `d_branch`, 1). `seeded:c-two-terms-1e-8` is judged by the sweep's ranking, which the self-test asks of the binary64 sweep only (ranking it against the `f32` sweep is not done), so it is not run at `f32` and a run of it there is an error, not the correct kernel's answer under its name, and the `Log` and SE(3) defects have no `@f32` stratum. The SO(3) half (`xtask/src/conformance/selftest_so3.rs`; readings 0014 (draft) questions 16 to 18): the correct `Exp` and `Log` must stay at most 4 `u` on every stratum of `so3_exp` and `so3_log` (a reading: a few roundings, 1.4 times the measured 2.91), with no non-finite output, `NaN` or unscored stratum, and fire no mechanism. `seeded:log-acos` is `θ = acos((tr R − 1)/2)` from the diagonal `1 − 2(a² + b²)` of the double `q`, the axis the vector part's, the flip kept; detected when every stratum of `theta:1e-k` from `k = 4` and of `theta:pi-1e-k` from `k = 7` reaches `10^7` (family maxima 9.0e15 and 8.5e7; from `k = 4` and `k = 7` as measured, 8.2e7 at `theta:1e-4` and 1.5e7 at `theta:pi-1e-7`). §10's "`theta:1e-k` and `theta:pi-1e-k` max ≥ `10^7 u`" cannot hold for every `k`: `θ` from `χ` errs by `≈ u/(2θ²)` (`docs/maths/so3.md` SO.6). `seeded:log-no-flip` is the correct `Log` without the flip; detected when every `so3_log` stratum with a negated half (`w < 0` in the input: 29 of 30, `q:w0` has none) either fails, its max at least `10^7` (27 strata, at least `1.8e16 u`: the error is `2π` in `φ`), or has nothing to detect, and one fails. §10's "every stratum" does not hold: in `theta:exact0` and `theta:subnormal` every negated record has `n² = 0` (derived from the records), where the series arm `2/w` is analytic in the sign of `w`, so the defect answers them bit for bit as the correct kernel does; a stratum answered as the correct kernel anywhere else is undecided and fails the gate. The corpus cannot pin `Exp`'s overall sign or `Log` at `w = ±0` (the metric aligns the sign): unit tests do. The SE(3) half (`xtask/src/conformance/selftest_se3.rs`; readings 0014 (draft) questions 21 to 25): the correct `Exp`, `J_r` and `J_l` must stay at most 12 `u` on every stratum of the nine ids (a reading: 1.35 times the measured 8.90 of the per-coefficient kernel, not of §4's grouped one; §10 states none) and fire no mechanism, and each planted defect must fire its own. "Fails" is a stratum at `10^7 u` or a non-finite output, the line of the SO(3) defects; a `NaN` does not fail, so it cannot fire a mechanism. `seeded:se3-exp-translation-first` reads the tangent as `[ρ₁; …; ρ_N; φ]`, the rotation block last at every `N` (`docs/maths/se3.md` SE.14); detected when every one of the 25 `rho:*` cells of each of `sen3_exp_n1`, `_n2` and `_n3` fails (75 of 75, 8.9e15 to 2.8e22 `u`; the 81 `theta:*` strata fail as well, at least 8.9e15). `seeded:q-minus-half` is `Q` with `−½` for `½` (so `J_r`'s `Q(−ρ, −φ)` too); detected when every stratum (52) of every `sen3_jr_n{1,2,3}` and `sen3_jl_n{1,2,3}` fails (312 of 312, 5.2e9 to 2.4e16 `u`), and a file with fewer strata fires nothing; the smallest is `rho:1e-6/theta=1e-8`, where the block is `5·10^-7` against §11's floor of 1 for a Jacobian. A mechanism reads only the ids the subject answers, so each defect is quiet on the other's. Missing: §10's `Dual` comparison of the `Q` defect (forward-mode differentiation of the shipped `Exp` is the `helicoid` subject, Phase 3), and `sen3_jr_inv*`, `sen3_jl_inv*`, which are built from `Q` (`docs/maths/se3.md` SE.12) and have no seeded subject |

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
`corpus` (the records above; `rechecked` counts those recomputed at 150 digits) or `series`
(`coeff_series`, §4.3; `verified` counts its series equal to exact algebra).

### 4.3 Definitions (v1 function ids)

| Function id | Input | Reference computed as |
|---|---|---|
| `coeff_k`, `coeff_a`, `coeff_b`, `coeff_c`, `coeff_d`, `coeff_e` | $\theta$ | the definition at 120 digits; derivative by `mp.diff` |
| `coeff_cos_half` | $\theta$ | $\cos(\theta/2)$, the other half of `Exp`'s quaternion (`NUMERICS.md` §3.1), at 120 digits; derivative by `mp.diff` |
| `coeff_r` | $(n, w)$ | $2\,\mathrm{atan2}(n, w)/n$ at 120 digits |
| `coeff_series` | — | the exact series of [`0004`](./decisions/0004-switch-points-are-generated-not-typed.md): `mp.taylor` of each definition at 120 digits, rationalized, equal term by term to an independent exact derivation |
| `so3_exp` | $\varphi$ | quaternion power series $\sum p^n/n!$, $p = (0, \varphi/2)$ |
| `so3_log` | $q$ | the $\varphi$ with $\|\varphi\| \le \pi$ and $\mathrm{Exp}(\varphi) = q/\|q\|$ ($w \ge 0$ after the flip), by Newton's method on the `so3_exp` series; `mp.logm` of $R(q/\|q\|)$ is a test cross-check only, since it returns complex results near $\pi$; stratum `q:w0` uses the sign rule of `NUMERICS.md` §3.2 |
| `so3_act` | $(q, p)$ | $R(q/\|q\|)\,p$ |
| `so3_from_matrix` | $R$ | the quaternion of the rotation nearest to $R$ in Frobenius norm (the polar factor): a rounded or scaled $R$ is not a rotation, so its `mp.logm` is not skew; backward error only (`NUMERICS.md` §11) |
| `so3_jr`, `so3_jl` | $\varphi$ | $\sum (\mp W)^n/(n+1)!$ |
| `so3_jr_inv`, `so3_jl_inv` | $\varphi$ | `mp.inverse` of the above |
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

A `coeff_k`…`coeff_e`, `coeff_cos_half` or `coeff_r` record's outputs are `value` and `d_branch`:
the derivative with respect to the branch variable of `NUMERICS.md` §4 (θ² for `cos_half`), by
`mp.diff` at the exact real branch value of the exact input (θ², not `fl(θ·θ)`), binary64, or
binary32 in an `@f32` stratum (§4.4).

`coeff_series` is the one file that does not follow §4.2: one record per coefficient (the seven of
`NUMERICS.md` §4 and `cos_half`), `{"id","coeff","branch","prefactor","series"}` with sixteen
`"num/den"` terms in the branch variable, and no `stratum`, `in` or `out`. Its manifest `kind` is `series`; the conformance
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
- `q:nonunit`: quaternions at $\|q\|^2 - 1 = \pm 2^{-45}$ (`Log`'s scale invariance, `from_*` normalization).
- `<S>@f32` ([`0016`](./decisions/0016-f32-exact-strata-for-the-scalar-coefficient-ids.md)), for
  every stratum $S$ of `coeff_k`, `coeff_a`…`coeff_e`, `coeff_cos_half` and `coeff_r`, after all of
  them: $S$'s inputs rounded to nearest-even binary32 (`coeff_r`: $n$ and $w$ each), the reference
  recomputed at those inputs. `theta:subnormal@f32` is its own decade $[10^{-40}, 10^{-39})$
  (rounding $[10^{-310}, 10^{-309})$ gives 64 zeros); `theta:pi-1e-8`…`theta:pi-1e-12` round to one
  binary32, kept five times.

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
  // Source sweep rev: sha256 <of the CSV>.
  // Source series: conformance/corpus/coeff_series.jsonl, sha256 <of the file>.

  // Objective (max u): value 1.9, derivative 3.4.
  pub(crate) const B_F64: Switch<f64, 5> = Switch::first(
      f64::from_bits(0x3f8d_...), // θ² < 1.44e-2
      &SWEPT_B_F64,
  );

  pub(crate) const SWEPT_B_F64: [f64; 8] = [1.6666666666666666e-1, /* correctly rounded rationals */ ];
  ```

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
