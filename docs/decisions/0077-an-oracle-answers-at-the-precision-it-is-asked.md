# 0077: an oracle answers at the precision it is asked

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** #148

## Context

The oracle file protocol (`PHASE1.md` §7, `xtask/src/conformance/oracle.rs`) has no precision.
Every runner answers every record in binary64, `@f32` records included, and `execute` registers
each with `no_f32`, so `cargo xtask conformance --oracle NAME --precision f32` is an error and
`FileSubject::eval` answers nothing at `f32`. `just envelope` reads the binary64 result files
alone (`xtask/src/envelope/mod.rs`, "Binary64 only"), although §8 states the bars per
`(fn, stratum, precision)`. So the `@f32` strata ([`0016`](./0016-f32-exact-strata-for-the-scalar-coefficient-ids.md)
and its successors) are scored for the candidate and judged against nothing: 1333 candidate rows,
`sim3_*`'s twins ([`0070`](./0070-sim3-a-structured-jacobian-and-a-series-on-delta.md)) among
them. [`0072`](./0072-a-container-oracle-is-an-image-on-the-same-protocol.md) left "a protocol
precision and `Sophus::Sim3<float>`" to *Further work*.

Of the runners' libraries, nalgebra is generic over `RealField` (implemented for `f32`), Sophus
and manif are templated on the scalar, GTSAM is `double` throughout, `tf_tree_math` is binary64
by construction (`tf_tree` D6), and `sophus_autodiff` 0.15.0 implements its scalar trait for `f64`
alone. Of the ids they answer, those with `@f32` strata are nalgebra's eight
(`eig3`, `chol_*`, `chol_solve_*`, `solve_cubic`, `quat_renormalize`, `mat2_inverse_adj`) and, of
Sophus's and GTSAM's, `sim3_{exp,log,ad,act}` and GTSAM's `so3_gamma2`. None of manif's 36 ids,
`tf_tree_math`'s or sophus-rs's has one.

Sophus's `Constants<float>::epsilon()` is `10⁻⁵` (`10⁻¹⁰` for `double`), so `RxSO3<float>`'s
constructors ensure a scale in `[10⁻⁵, 10⁵)`, `|σ| < 11.5`.

## Decision

1. **The protocol.** `<runner> [--precision f32] --out DIR FILE.jsonl…`. Without the flag it is
   today's binary64 protocol, unchanged, so a runner that answers binary64 only needs no change.
   With `--precision f32` the runner writes `DIR/<fn>.jsonl` for each file whose function it
   answers in binary32, one line per **`@f32` record** (the stratum's `@f32` suffix, `0016`), ids
   in order, and no line for any other record. It narrows each input to binary32 **losslessly**
   (an input that is not a binary32 is an error, never a rounding), computes with its library's
   binary32 instantiation, and writes the binary32 results widened to binary64, so any hex-float
   spelling the binary64 protocol accepts is one. `--version` is unchanged and precision-free:
   one pin, the `subject_version` of the rows at both precisions; the row's `precision` column
   names the instantiation.
2. **The declaration.** `Runner::answers_f32` lists the ids a runner must answer in binary32: every
   id of its `answers` whose corpus file has `@f32` strata and whose library has a binary32
   instantiation (`every_owed_function_is_a_corpus_function_of_a_runner_that_exists` holds the
   list to that rule). Empty means binary64 only: an `f32` run of the runner launches nothing,
   says "answers binary64 only; no f32 rows" and passes, so an `f32` pass over every runner does
   not fail on GTSAM; and the runner, asked anyway, refuses `--precision` as any unknown argument
   rather than answering in binary64.
3. **The harness.** At `f32`, `execute` passes the files that mention `@f32` and nothing else
   (`--fn` an id with none is an error), into `results/oracle/<name>[--<fn>]--f32`, and owes
   `answers_f32` where it owed `answers`. `FileSubject` holds one precision's answers by record
   `id`; ingestion at `f32` requires the file's ids to be the corpus's `@f32` ids in order (a
   skipped, extra or binary64 record is an error, not a score). It is then scored by `evaluate`
   exactly as an in-process `f32` subject: `u = 2⁻²⁴`, the inputs asserted binary32, every output
   asserted binary32 (an answer computed in binary64 fails the run), a non-finite answer recorded.
4. **The envelope.** Every subject's rows are its `<name>.csv` and its `<name>--f32.csv`, each read
   when it exists, a file holding its own precision's rows only; bars pair by
   `(fn, stratum, precision)` as before. An oracle with `answers_f32` owes its `f32` file; an
   oracle's `f32` rows with no candidate `f32` rows is an error, not unpaired strata. `just
   envelope` runs `conformance --precision f32`; `oracle-nalgebra` and `oracle-sophus` run both
   precisions. The twin has no `f32` kernel, so an `f32` failure no `libm`-crate oracle wins is
   reported unattributed.
5. **The runners.**
   - **nalgebra**: each id is one generic function at `f64` and `f32`; the iteration tolerance of
     `SymmetricEigen` and `Schur` is `T::default_epsilon()`, `f32::EPSILON` at binary32. simba
     with `std` computes `f32` with `f32`'s own `std` methods.
   - **Sophus**: `At<T>` evaluates every id over `Sophus::*<T>` and `convert.hpp` is templated on
     `T`; binary32 answers `sim3_{exp,log,ad,act}`. A scale outside `RxSO3<float>`'s domain is
     refused, `nan`, as `0072` decision 5; `Sim3<float>::exp` clamps, and is scored.
   - **GTSAM, manif, `tf_tree_math`, sophus-rs**: binary64 only, unchanged. manif could
     instantiate at `float`, but no id it answers has an `@f32` stratum, so `answers_f32` would be
     empty and the path never run.
   - Images keep `-ffp-contract=off -DEIGEN_DONT_VECTORIZE`; `float` arithmetic is binary32 on
     x86_64 and aarch64 (`FLT_EVAL_METHOD` 0), `std::exp(float)` and its kin are the `float`
     overloads, and the `sim3`, `rxso3`, `so3` and `sim_details` paths these ids take hold no
     unsuffixed `double` literal (the ones in `so3.hpp` and `se3.hpp` are in `Dx_exp_x`, not
     called).
6. **The hand tests.** nalgebra's `renormalize_in_binary32_rounds_in_binary32`: `q = (1 + 2⁻¹⁰)
   e_w` gives `1 − 3·2⁻²¹` in binary32 and `1 − 3·2⁻²¹ − 2⁻³¹` in binary64, so a promoted step
   fails; the narrowing's refusals and the quaternion order at `f32`. Sophus's `convert_test`: the
   narrowing (`0.1` refused, `0x1.99999ap-4` accepted and printed as such), the quaternion,
   tangent and matrix orders at `float`, `Sim3<float>`'s stored `√s·q`, `σ` read back, `Ad` and
   `exp` against the hand values of the `double` tests, every output a binary32, and `e^12`
   refused at `float` and not at `double`.

## Rationale

- **A flag, not a per-file mode.** The harness already runs one precision per invocation
  (`--precision`, `results_path`); a run that answered both precisions per file would need two
  answers per `@f32` record and a second answer grammar. The flag keeps each answer file one
  precision's, and the binary64 protocol byte-identical, so the four binary64-only runners and
  their images are untouched.
- **The `@f32` records alone at `f32`.** A binary64 record's inputs are not binary32, so answering
  it in binary32 would mean rounding them, which is a different function; answering it in binary64
  would be another kernel's answer under this precision. Which records are `@f32` is the corpus's
  definition, not the runner's choice, so §7's "every stratum of the file" is kept per precision.
- **The declaration in the harness, as `answers` is.** The registry already holds the owed ids
  "so that an id it stops answering fails the run"; a binary32 list beside it gets the same
  guarantee, and a test ties it to the corpus. A `--version` that listed precisions would be
  checked against the same list and say nothing more.
- **A binary64-only runner passes at `f32`.** Failing would make an `f32` sweep over the runners
  fail on GTSAM for a fact the registry already states; launching it would cost a container start
  for nothing.
- **Losslessly or not at all.** The metric asserts binary32 inputs and outputs for in-process
  subjects; a runner that rounded an input or promoted a result would be scored on another
  function, silently.
- **One pin.** `Sophus::Sim3<float>` and `<double>` are one source at one commit under one
  compiler; a second version string would only restate the `precision` column.

## Consequences

Measured on `8c53526` + this branch, with every oracle's rows at both precisions it answers,
`just conformance-twin`, `cargo xtask conformance` at both precisions, and `cargo xtask envelope`
(no baseline in the tree, so no no-regress verdict):

- **No binary64 row moves.** `nalgebra.csv` and `sophus.csv` are what they were, column for column
  but `git_rev`.
- **Pairing.** 489 binary32 strata newly have an oracle: `sim3_exp` 116, `sim3_log`, `sim3_ad`,
  `sim3_act` 104 each (Sophus), and nalgebra's 61 (`eig3` 19, `solve_cubic` 12, `chol_*` 6 each,
  `chol_solve_*` 5 each, `quat_renormalize` 5, `mat2_inverse_adj` 3). The other 12 of each of
  `sim3_{log,ad,act}` (`|σ| ∈ {30, 80}`, 144 records) are refused by `RxSO3<float>`'s scale, so
  unpaired; `sim3_exp` there is clamped and scored at `1.0·10⁷`–`1.4·10⁷ u`. 844 binary32 strata
  still have none, among them every `sim3_j*` and `so3_gamma2` stratum. Totals: 2105 → 2594
  paired, 917 → 1761 with none.
- **Domination failures**: 707 → 805, the 98 new ones all binary32 and all unattributed (no
  binary32 twin): `sim3_log` 25, `sim3_act` 19, `sim3_ad` 18, `sim3_exp` 12, `eig3` 16,
  `chol_solve_n6` 4, `chol_n6` 2, `chol_solve_n3` 1, `quat_renormalize` 1. The 3 exceptions are
  unchanged.
- **Sim(3) at binary32** fails 74 of its 428 paired strata, by 1.00×–2.53× on maxima of at most
  5.01 `u`: `sim3_exp` 12 of 116 (worst `sigma:3/theta:1e-4@f32`, 1.415 against 0.5904 `u`),
  `sim3_log` 25 of 104 (`sigma:10/theta:3@f32`, 1.121 against 0.4424), `sim3_ad` 18 of 104
  (`sigma:-3/theta:0.3@f32`, 0.9247 against 0.4216), `sim3_act` 19 of 104
  (`sigma:3/theta:0.3@f32`, 1.784 against 0.7059). The `R(q)` diagonal of
  [`0076`](./0076-the-oracles-first-verdict-r-of-q-stays-and-drift-is-excepted.md) decision 1 and
  `f32`'s transcendentals are each a candidate cause; neither is measured here. `PHASE5.md` §7's
  Sim(3) line stays open.
- **nalgebra at binary32** repeats its binary64 verdict: `eig3` fails 16 of 19 (to 481× at
  `eig:gap-1e-4/top@f32`, 3683 against 7.66 `u`, and `eig:triple@f32`, 2.379 against an exact
  0), as 27 of its binary64 strata already do; `chol_*` and `chol_solve_*` within 1.53×,
  `quat_renormalize` 1.19×.
- `just envelope` now requires `sophus--f32.csv` and `nalgebra--f32.csv`, which `oracle-sophus`
  and `oracle-nalgebra` write.

## Implementation plan

1. One PR: this record; `Runner::answers_f32`, the `f32` path of `execute`, `FileSubject` and
   `conformance --oracle`, and their tests (`a_binary32_run_asks_for_the_f32_files_and_scores_at_f32`,
   `a_binary32_file_holds_the_f32_records_and_no_other`, `a_binary64_runner_is_not_launched_at_f32`);
   the envelope's `f32` files (`f32_rows_are_judged_against_the_oracles_f32_rows`); the nalgebra
   runner and the Sophus image at `f32`; the recipes; `PHASE1.md` §7 and §0.0. Verified by
   `cargo nextest run -p xtask -- oracle envelope`, the Sophus image build (it runs
   `convert_test`), `just oracle-nalgebra`, `just oracle-sophus`, the other `oracle-*` recipes,
   `just envelope`, `just lint`, `just test`, `just doc`.

## Open questions

None.

## Further work

- A binary32 twin (`helicoid:host-std` at `f32`, or a `libm`-crate oracle that answers binary32)
  to say which of the 98 are D16's: the runners that answer binary32 are both on the host's `std`.
- An oracle for `sim3_j*` and `so3_gamma2` at binary32: GTSAM has `Γ₂` but no `float`, and no
  runner exports a Sim(3) Jacobian.
- manif in binary32, with the record that gives one of its ids an `@f32` stratum.
