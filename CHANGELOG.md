# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

**`0.0.x` is not ordinary semver.** Every release may break every other. What is implemented is
defined by the status tables in `docs/`; they win over this file.

## [Unreleased]

### Added

- **`LieGroup::geodesic` and `geodesic_velocity`** (`PHASE4.md` §1.1, `NUMERICS.md` §10), provided
  methods on every group, with `reference::geodesic` as the expression they call — so the provided
  body **is** the reference twin, as §1.1 asks, and a group's fast twin has something to be
  measured against that its own override cannot shadow.
  - `Product` delegates per factor, so **`Product<SO3, Rn<3>>` is slerp and lerp** by construction
    (§1.3): `tf2`'s semantics without `tf2`'s small-angle fallback, which `0009` leaves to
    `tf_tree`.
  - **`SO3::geodesic` is deliberately not overridden.** GE.14 proves the provided body is
    shortest-arc slerp with `Log`'s sign rule, in the better-conditioned spelling — `atan2` where
    slerp has `acos`, whose slope is infinite where the two quaternions are nearly equal — and the
    corpus now says so rather than the maths page alone: **1.572 `u` against
    `tf_tree_math::slerp`'s 2.187** on `geo:consecutive`, the kilohertz stratum. It also removes
    `Real::acos` (`0022`, unimplemented) from the path entirely.
  - **And where it is behind, with the numbers**, because that is the point of having the
    instrument first: `so3_geodesic` 2.721 / 2.429 `u` against 1.834 / 1.642 at `geo:generic` /
    `geo:near-pi`, and `se3_geodesic` 2.556 / 3.240 / 3.466 against `ScLerp`'s 2.336 / 2.502 /
    3.253. Five domination failures, and the motivation for §1.2's screw twin stated in numbers:
    the provided body is `Log` then `Exp` where `ScLerp` is one `atan2` and one `sin_cos`, so the
    reference pays for roundings the screw form does not have.
  - **`laws::geodesic`**, six legs in one law for every instantiated group (`Rn`, `SO3`, `SE3`,
    `SE_2(3)`, the test-only Heisenberg group and the four products) at `f64`, `f32` and
    `Dual<f64, D>`: both endpoints, GE.2(c)'s symmetry, the velocity, and GE.4's left **and**
    right invariance. Right invariance holds on *every* group, not only SE(3) —
    `Log(h⁻¹Δh) = Ad_{h⁻¹}d` and `h Exp(Ad_{h⁻¹}ξ) = Exp(ξ)h`, no commutativity used — so
    `Product<SO3, R3>`'s famous failure is against the **SE(3) reading** of `(R, t)` and not
    against `Product`'s own `Mul`, which `PHASE4.md` §3 asserts separately. `x₁` is `x₀ ⊕ d` and
    not a second draw, so `θ(d)` stays inside the domain: two draws of `[-1,1)^D` compose to up to
    `2√3 > π`, where the law would measure GE.13(d)'s conditioning instead of the identity.
  - Bounds are twice the worst of **10^6 draws** of a stream of its own, per group, printed per leg
    by an `#[ignore]`d `measure_geodesic`: `Rn` 10, `Product` 18, `Heis` 19, `SO3` 26, `SEn3` 71
    (its `symmetry` leg at `N = 2`, `f32`, is the binding one at 35.339). The `velocity` leg is
    exactly **0** on every group, and `t = 0` returns the left endpoint **bit for bit** on SO(3)
    and SE_2(3) — which the law's `gerr` cannot say, its floor for two bitwise equal quaternions
    being 1.118 `u`, so two dedicated tests say it instead.
  - `xtask/src/shipped.rs` gains a `Geodesic` family, so the two corpus ids are scored from the
    moment they exist.

- **The geodesic instrument: corpus ids `so3_geodesic` and `se3_geodesic`, 180 records each**
  (`PHASE4.md` §4, [`0045`](docs/decisions/0045-two-phase-4-checks-cannot-be-taken-as-written.md)
  plan steps 2 and 3). D7 puts the instrument before the routine, and Phase 4 had none: these are
  the first corpus ids of the geodesic, and the first in the repository whose reference is **not**
  `mp.expm`/`mp.logm` of the matrix.
  - **`mp.logm` is the cross-check, never the reference.** It returns a complex, non-principal
    logarithm of a rotation from θ = 3.03, which is exactly the `geo:near-pi` stratum — so a
    corpus generated from it would be wrong by `O(1)` on the one stratum that is hardest, and the
    hunt would have gone looking in `helicoid`. The reference is the geometric `Log` (the
    `so3_log` quaternion `atan2` route, with `sen3_log_n1`'s `V` solve for the translation block)
    and `sen3_exp_n1`'s matrix exponential of `t·d`. Measured over the whole corpus, the two agree
    to **8.2e-113** below θ = 3 (6.8e-121 at `geo:consecutive`), so the matrix route is held to the
    same 100 digits as every other cross-check here, with thirteen to spare.
  - **Every record is checked by identities as well**, each costing one further evaluation of the
    reference at a different argument: the two endpoints, the symmetry
    `γ(X₀,X₁,t) = γ(X₁,X₀,1−t)`, and `Log(X(t)⁻¹X(t+⅓)) = ⅓d` — constant body velocity without a
    difference quotient. `geo:near-pi` is above the `mp.logm` limit at every record, so there it is
    the whole of the check.
  - **Strata** (`PHASE1.md` §4.4): `geo:consecutive`, `geo:generic`, `geo:near-pi`, 6 pose pairs
    each crossed with 10 parameters (§4's six, including both endpoints, and four uniform), with
    `‖x₀‖` cycling 1, 10², 10⁴ in every one of them. `X₁ = X₀·Exp(d)` is **formed and rounded,
    never drawn**, so a stratum's relative motion is the `d` it names. `geo:near-pi` stops strictly
    below π at `π − 10⁻ᵏ`, `k ∈ {1,2,3,6,9,12}`: at π the two preimages give geodesics `O(1)`
    apart, so that stratum measures conditioning and not agreement.
  - **What the strata found before any routine exists.** At `‖x₀‖ = 10⁴` with `θ(d) = 10⁻⁹` the
    stored `x₁ − x₀` keeps about **three digits** — the relative motion is below the ulp of the
    poses it is between. That is `tf_tree`'s kilohertz edge exactly, and the reason `PHASE4.md` §4
    asks this stratum for that scale; the records are the geodesic between the two poses they
    hold, so it is the regime and not a defect, and `PHASE1.md` §4.4 now says so with the bound.
  - **Verified by 13 generator tests**: four hand-computed interpolations (a quarter turn halved,
    a pure translation lerped, a screw about `z` halved in both parts, a quarter turn at unit
    offset sweeping the arc), GE.4(a)'s left invariance and the matrix route as properties no
    cross-check uses, the stated margins and band — and **198 mutants** of the cross-checks, of
    which the only survivors are the 36 that negate the *whole* quaternion, which is the same
    rotation and what the metric's `SignRule::Align` forgives for this id too.
  - The `tf_tree_math` oracle answers both (`slerp` and `ScLerp`), hand-computed against the same
    anchors the generator uses, so `PHASE4.md` §5.1's parity rows are measurable the moment
    `LieGroup::geodesic` lands. The corpus is 37 MB of the 50 MB budget, from 34.

- **The `tf_tree_math` oracle answers three more corpus ids** — `so3_act` (`Quat::rotate`),
  `sen3_ad_n1` (`Iso3::adjoint` on the six basis twists) and `so3_from_matrix`
  (`quat_from_rot3`) — which the pin has always exported and `PHASE1.md` §7 did not list. It
  answers 278 rows where it answered 164, which is 114 strata that were paired with nothing. This
  settles [`0014`](docs/decisions/0014-the-conformance-metrics-open-readings.md) (draft) question
  32 for that runner and `0041`'s scheduled `quat_from_rot3` reading, and it is the measurement
  `PHASE4.md` §5.1's parity table is built from.
  - `so3_act` is **bit-identical to `helicoid` on all 30 strata** — the second clean control after
    `so3_log` (`0032`, draft). Both programs spell the action `p + 2w(u × p) + 2u × (u × p)`, so
    both read `q:nonunit` at 510.3 `u` and neither is the reference's `R(q/‖q‖)p`.
  - `sen3_ad_n1` differs on 42 of 54 strata: `helicoid` is smaller on 20, among them `q:nonunit`
    (256.1 `u` against 418.8 — `Quat::rotate` at a non-unit `q` is neither `R` nor `‖q‖²R`, where
    `to_matrix` is `‖q‖²R` exactly), and **larger on 22, at worst 3.30x** (`theta:1e-2`, 1.569 `u`
    against 0.475). `helicoid`'s rows did not move; the pairing is new, so these are 22 newly
    *visible* domination failures, not a regression — which is the whole of D7. `SEn3::adjoint`
    forms `R` once with `to_matrix` and multiplies; `Iso3::adjoint` applies the sandwich per basis
    twist. The obvious guess is already refuted: `docs/maths/so3.md` SO.6 measures the diagonal
    `1 - 2(·² + ·²)` as *worse* near π, by 16x through the trace. Left to the `0038` paired
    instrument with the numbers above.
  - `so3_from_matrix` is scored by backward error, which `NUMERICS.md` §11 owes, so its 30 rows
    record shape and finiteness (NaN `max_u`, no non-finite answer) and pair no number. The
    corpus's `3 x 3` is **column-major** and `quat_from_rot3`'s is **row-major**: that conversion
    is the one in this runner that is not a layout no-op, and its hand-computed test asserts that
    the wrong reading returns the conjugate — a unit quaternion, silently the inverse rotation.
  - Five stale `0014` question numbers in `PHASE1.md` §0.0 row 20 corrected (the record's
    numbering shifted by five when questions were inserted); sophus-rs's own `adj` and `transform`
    are now owed rather than blocked.

### Changed

- **The coefficient kernel got its speed back: every switch now has a second, shorter arm**
  ([`0039`](docs/decisions/0039-the-sweeps-grid-stops-below-its-own-optimum.md) items 5, 6 and 10,
  with [`0047`](docs/decisions/0047-the-second-arm-is-admitted-by-agreement-not-by-the-objective.md)
  deciding what admits one). A series arm's term count is set by the largest branch variable it
  serves and paid by the smallest, so lifting the search space to reach config E's accuracy left
  every near-identity input paying the hardest one's price — 11 to 16 terms where 3 to 6 give the
  same bits. `coeffs::Switch` now carries `short_below` and `short_terms`, the first `m0` terms of
  the same series below a second generated grid point, and the sweep is two stages: stage 1 is the
  objective, unchanged, and stage 2 is a **cost** search afterwards, because a shorter arm under a
  maximum the long arm already attains cannot lower it.
  - Corpus-weighted term count, one term one unit: **172427 -> 92264 at `f64` (1.87x)** and
    **119547 -> 60082 at `f32` (1.99x)**. Every one of the sixteen rows gets a second arm: `m0 = 5`
    for seven of eight at `f64` and 3 for six of eight at `f32`.
  - **Measured, at the stratum `PHASE3.md` §11 names first.** `cargo xtask bench-gate --bench coeffs
    --only near-identity-7.5e-8 --against` a baseline built in this tree from `main`, point ratio
    (`f32` / `f64`): `jr_inv_coeff` **0.2353 / 0.2605** (4.3x and 3.8x faster), `jr_coeffs`
    **0.5282** / 0.7476, `exp_coeffs` **0.5525** / 0.8760, `log_ratio` 0.5975 / 0.6365, `q_coeffs`
    0.7826 / 0.6585. **Ten of ten faster**, floors 0.0038 to 0.0537, every effect outside its own.
  - **And what it costs, which the near-identity table does not say.** At `near-pi`, where every
    group is on its exact arm and the extra comparison is all that is left: `jr_inv_coeff`
    1.0014 / 0.9997, `q_coeffs` 1.0190 / 1.0081, `jr_coeffs` 1.0204 / 1.0350, `exp_coeffs`
    1.0275 / 1.0269 — 1% to 3.5% for one comparison against a `sqrt` and a `sin_cos` — and
    `log_ratio` 1.1288 / **1.3036**. That last is 1.37 ns, about four cycles, on a mask chain that
    already waits on the division forming `s` and in front of a body that is one `sqrt` and one
    `atan2`. `#[inline]` on all five entry points was measured against it and changed nothing, so
    it is not shipped. Swapping `log_ratio`'s two switch tests would move the 1.37 ns onto its
    series arm instead; `0047` *Further work* holds it with the numbers that decide it.
  - **A call-site group reads one second switch, its members' smallest** (`0047` item 7) — and that
    is the measurement that changed the design. A mask and a branch *per member* costs one compare
    and one branch on every call, including every call above the second switch, and it measured
    **1.17x, 1.15x and 1.15x slower** than one arm on the three multi-member groups at `f32`
    near-identity, and 1.20x to 1.28x slower above the second switch. One switch for the group puts
    the short arms one compare from the entry point and leaves the rest of `grouped` as it was. It
    costs no accuracy: between the group's switch and a member's own, the member takes the *whole*
    arm, which is the arm the prefix was admitted to agree with there.
  - **Not one bit of any measured row moves.** A prefix is admitted only where it agrees with the
    whole arm **to the bit** — at every grid point and every corpus record below the second switch,
    in the value and in `d/dz` — which is why all nineteen files under `conformance/results/` are
    byte-identical to the committed ones apart from the git revision each run stamps, and why no
    stage-1 cell of `thresholds.csv` changes.
  - `0039` item 6's literal rule, "the cheapest arm that *holds* stage 1's objective", was measured
    and **refused**: 5.8% cheaper at `f64`, 9.3% at `f32`, and it raises individual corpus records
    by up to **50224x** (`e` at `f32`), because the objective is one maximum and `0006`'s bars are
    per stratum.
  - Between the points stage 2 checks, the two arms part by at most **1.875 u** at `f64` and
    **3.000 u** at `f32` — the last place or two of the shorter arm, against a catalogue whose
    gentlest objective is 1.747 u. Carried as a no-regress bound.
- **`log_ratio` takes its second arm at one division fewer.** Its short arm's Horner argument is
  `s = n^2/w^2` itself, made safe by a `select`: on a selected lane `small` implies `w > 0`, so the
  mask's own `select` took `w * w` and `s` is `n2 / (w * w)` to the bit, where the whole arm
  re-divides at a re-selected `w`. A division is 1.11 ns against the 1.41 ns of a five-term arm.
  The whole arm's path is left exactly as it was, deliberately: measured with the prefix selected
  *inside* the series arm — which also removed that division — the function was **1.16x slower** at
  `theta = 0.5`, where `s = 0.065` takes the whole arm. One comparison on the way in buys the
  prefix; a comparison under the series arm buys it and pays for two Horner bodies in one block.
- **The kernel's module doc was stale and is corrected.** It still described the switches the old
  grid produced — `cos θ/2` at `θ² < 5.6e-15`, so "`Exp` is on its exact arm from `θ = 7.5e-8`", and
  "eight of the sixteen switches are `θ² = 1`, the top of the sweep's grid". Under the committed
  sweep `exp_coeffs` is on its series arms to `θ = 2.29`, and six switches sit on the **domain**
  bound with the file printing which limit each sits against.
- **The threshold sweep searched a box too small, and lifting it moved every coefficient**
  ([`0039`](docs/decisions/0039-the-sweeps-grid-stops-below-its-own-optimum.md) plan step 0 and
  decisions 1-4, 7-9). `PHASE1.md` §6's grid stopped at `z = 1` while `NUMERICS.md` §12's domain
  runs to `pi^2`, and the cap stopped at 8 terms while the corpus commits 16. The two were
  **jointly** binding: neither alone got the worst swept row below 1846 u, both together reach
  **94.48**. Worst row per precision, against the D12 prior it replaced: `f64` `c` 2.8e7 -> 94.48 u,
  `e` 9.2e9 -> 19.87, `d` 1.2e7 -> 7.90; `f32` `e` 8.1e9 -> **15.65**, `d` 1.7e7 -> 12.96. No row is
  above 100 u at either precision now. The **D12 prior beats `helicoid` on two strata instead of
  nine**, and the two are near `pi` where both candidates are on their exact arms, so the margin is
  the assembly's and not a switch's.
  - The grid now **spans** the domain and the selection **rule** stops below it: a grid is a
    geometric object, a domain is a semantic constraint, and keeping the bound in the rule leaves
    `the_grid_has_64_points_per_decade_and_both_ends` a statement about the grid alone.
  - **A boundary optimum is no longer reported as a result.** A choice on the grid's first point
    **fails the run** — it is the edge of the box and nothing in the output would say so, which is
    how eight switches sat on a wall for the generated file's whole existence. A choice at the last
    point below `pi^2` is the *domain* binding and a choice at 16 terms is the *cap's*; both are
    printed, and every generated switch now carries its grid span, its cap and which limit it sits
    against. Under this configuration three `f64` and three `f32` rows sit at the domain bound and
    `c` at `f64` takes every term the corpus holds — exactly what `0039` predicted.
- **`coeff_switch_ref.jsonl`: the reference at every switch the sweep can return** (`0039` plan
  step 0; 17 424 records, 2.7 MB, `manifest kind: "switch-ref"`). `docs/maths/coefficients.md`
  CO.12 bounds the jump between a coefficient's two arms at its switch by the sum of the arms'
  errors *there*, and sampling that right-hand side at the two corpus records bracketing the switch
  understates it by up to **3.5x** — the exact arm's error is a sawtooth, swinging 195x over 0.8% of
  `theta`, so no density of records makes a point sample a bound. `coeffs::tests::branch_continuity`
  now checks CO.12 as stated at **both** precisions against the true value and `d/dz`, with the
  reference's own half-ulp charged to the bound once per arm.
- **Honest accounting of what this cost.** Paired with `tf_tree_math` at `f64` over 164 strata,
  `helicoid` is worse on **17 against 13** before: four new marginal regressions, all at
  `theta ~ 1`, where the switch has moved and `theta = 1` now takes a 13-to-16-term series instead
  of the exact arm — `sen3_log_n1` at `rho:1e3/theta=1` (1.045x) and `rho:1e4/theta=1` (1.042x),
  `sen3_exp_n1` at `theta:1e0` (1.043x), `so3_exp` at `theta:dense` (1.034x). A per-call-site switch
  (CO.18, `0015` draft NU.5) is the lead.
  And it **retires a diagnosis this project has been carrying**: `sen3_log_n1`'s 2.2x loss at
  `theta=1e-1` is unmoved to four digits although `coeff_c` there is now **1.54 u**, so it was never
  `c`'s generated switch against D12's typed prior. It is a two-program question about the assembly
  of `J_r^-1 = I - W/2 + cW^2`, which is what `0038`'s paired comparison exists to answer.

- **Phase 4's records: five decided, two specs corrected, and the first blessed envelope unblocked.**
  No code; the contract the geodesic work is implemented against.
  - [`0041`](docs/decisions/0041-the-integration-is-an-adapter-at-tf-tree-math.md) is **ready** and
    on `main` for the first time -- `PHASE3.md` §0.0 and this file have been citing a record that
    existed only on a branch. Two stale claims corrected on promotion: its keystone blocker
    ("`helicoid` has no group type") is cleared, so its plan steps 2 and 4 are done; and Wave 1
    contains `slerp`, so Wave 3 is not "the only wave needing Phase 4" -- though SO(3) does **not**
    override the provided geodesic, which `docs/maths/geodesics.md` GE.14 proves already *is*
    shortest-arc slerp, in the `atan2` spelling rather than `acos`. Its open question 1 is decided:
    `helicoid` and `helicoid-linalg` publish at `0.0.1` before Wave 1 lands, because `cargo publish`
    rejects a path or git dependency, optional or not.
  - [`0043`](docs/decisions/0043-the-geodesic-jacobian-ships-the-cancellation-free-form.md) edits
    `NUMERICS.md` §10. Its `J_0` is a difference of two `O(1)` matrices vanishing like `1-t`, with
    relative error `O(u/(1-t))`: **19 u at `1-t = 1e-1` and 1.3e12 u at 1e-12**, against **<= 0.72 u**
    for GE.7(a)'s equal form, which is also exactly `0` at `t = 1`. §4 samples `t = 1 - 1e-9`, so
    this is a stratum the spec requires. §10 now states both arrangements, which one ships, and the
    **left** pair it never had. `Jac` gains exactly one method, `scale` -- no `add`, no `sub`, no
    fused form, because `Ad - J Ad` is what an FMA contracts and D16 forbids the contraction.
  - [`0044`](docs/decisions/0044-four-primitives-the-first-consumer-names-and-no-spec-does.md): the
    four names `0041` delegates to that no spec lists. Two are holes, not conveniences -- **`SEn3`
    cannot be renormalized from outside at all** (its `q` is private while `Mul`'s rustdoc makes
    renormalizing the caller's step), and **`mul_inv` does one rotation of a vector where
    `a * b.inverse()` does two**, so the two are equal exactly and not bit-identical, which makes it
    a routine with its own `NUMERICS.md` §14 twin. `API.md` §3's `Quat` row was stale by four
    shipped methods and is corrected.
  - [`0045`](docs/decisions/0045-two-phase-4-checks-cannot-be-taken-as-written.md): two NORMATIVE
    Phase 4 checks were unexecutable. `PHASE4.md` §4's `mp.logm` reference returns a **complex,
    non-principal** logarithm from `theta = 3.03` -- error 1.0 to 1.9, mpmath 1.3.0 and 1.4.1 alike
    -- so the stratum §4 requires (`geo:near-pi`) is exactly the one its stated reference gets wrong,
    by `O(1)`; the reference becomes the geometric `Log`, with `mp.logm` kept as the cross-check only
    where it is sound. And §3's "`Product<SO3, R3>` fails right-invariance, `max_err > 1e-6`" is
    **vacuous** under the product's own law, where GE.5(a) proves exact bi-invariance at 1e-111: §3
    now names the SE(3) reading and carries the fixture inequality `theta^2 ||t_perp|| >~ 1e-5`.
  - [`0046`](docs/decisions/0046-explained-by-record-needs-a-record-to-point-at.md): `PHASE3.md` §10
    and `PHASE4.md` §5.2 both let a lost stratum be "explained by record", and **neither hatch
    existed** -- so `--bless` writes nothing, `docs/evidence/` is empty, and **no baseline has ever
    been blessed**, leaving D8's no-regress half switched off for 1529 rows. A committed
    `conformance/baseline/exceptions.toml` excepts **domination only**, needs a `ready` record, and
    fails the run when the stratum is dominated after all. Three strata are excepted on the measured
    ground that `helicoid` and `tf_tree_math` agree there to **ten significant digits** (5e-11
    relative, both under 1/2 u); the other ten paired losses are not.
  - [`0038`](docs/decisions/0038-a-program-comparison-is-not-a-bar.md) is **ready**: the max stays
    the gate, a program comparison becomes `cargo xtask compare` -- a reporting task, no result CSV,
    no baseline, read by no bar -- pooled per id then combined, reporting a magnitude and not only a
    direction. Ranking programs by counting domination failures is retired.

- `cargo xtask envelope` reports **which kind of oracle** beat the candidate. Each runner declares
  where its transcendentals come from (`conformance::Backend`): `tf_tree_math` is on the `libm` crate
  (`libm.workspace = true`), sophus-rs is Rust `std` and so the host's library. Every domination
  failure now says whether a `libm`-crate oracle also beats it -- in which case the gap is the
  program's and not D16's -- and the summary counts both classes, which on the current corpus reads
  30 and 45 of 75 (`0036`, draft, plan step 1). The test is whether **any** same-backend oracle beats
  the candidate, not whether the best one does; those questions differ here by 30 against 26. An
  oracle with no declared backend is left unclassified rather than assumed.
  Also corrects the envelope's module doc, which said the oracle runners "are measured on the host's
  `libm`": true of sophus-rs, false of `tf_tree_math`, and that difference is what makes it `0032`'s
  control.
- **`bench-gate --against` runs its replicate loop outermost**, a full pass of all 60 benchmarks per
  replicate. Measured: with it innermost, a benchmark's six candidate pairs were all taken inside
  about 30 s, so they were not independent, and **2 of 60** identical-code benchmarks failed with a
  quiet concurrent control -- the all-pairs rule survives one disturbed window but not a disturbance
  outlasting the whole span, which leaves the control quiet because it covers both control windows
  equally. A middle-window bias was the first explanation and the data refused it: over 60 benchmarks
  the median candidate-to-baseline ratio is 0.9994 and 0.9989 under two `criterion` versions. The
  windows within a triplet stay contiguous, the same total windows are measured, and a benchmark's
  three triplets are now about 20 min apart (`0033`, draft).
- `cargo xtask lint` checks that `libm::` appears only in `helicoid-linalg`'s `float.rs`, the one
  private kernel D16 routes every transcendental through. It holds today by construction -- six calls
  there, none in `Dual`, which inherits through `Real` -- and `no_std` enforced only half of it, since
  removing `std` makes `f64::sin` as an inherent method vanish but stops nobody adding a second call
  site on purpose. Tests and doc comments are exempt: a test pinning `Real::atan2` against
  `libm::atan2` is the twin comparison D6 asks for (`0034`, draft, plan step 1).
- `criterion` moves from `=0.5.1` to `=0.8.2`, `sha2` to 0.11 and `num-bigint` to 0.5 (dependabot).
  The pin stays exact: §9 pins it so the harness's behaviour is a deliberate choice, and 0.8.2 was
  checked against what `bench-gate` actually depends on -- `--save-baseline` slots persist,
  `--list`/`--exact` are unchanged, and `sample.json` still carries `iters`/`times`/`sampling_mode`
  with 100 Linear samples. `criterion::black_box` is deprecated there, so the benches take
  `core::hint::black_box`. The corpus digests are byte-identical under `sha2` 0.11. The new `alloca`
  -> `cc` chain is criterion's, so a C toolchain is needed to build the benches and nothing else:
  dev-dependencies do not reach a consumer or `just build`.
- `cargo xtask lint` fails on an **untracked, non-ignored file**. Every other check reads
  `git ls-files --cached`, so a file written and not yet added is invisible to all of them -- and
  every new decision record is untracked at the moment it is written. Measured here: a citation of
  `0033` (draft) in `docs/PHASE1.md` §0.0 was accepted by `cargo xtask lint` *and* `just lint`, and
  became a `drafts` violation the moment `git add` made the record visible, after the commit that
  claimed the lint was clean. `.gitignore` still excludes what it should.

- `solve_cubic` internals, every one bit-identical (`bits_are_omnisacs` and the planted-root suite
  are unmoved): `pi` is a per-precision literal rather than `atan2(+0, -1)`, which was a `libm` call
  for a compile-time constant that no workspace LTO could fold — `eig3` paid it per call too; `tol`
  takes its exponent as a const parameter and evaluates `pow2` in a `const` block, so the thresholds
  are literals and `pow2`'s stated domain is a compile error rather than a silently wrong float;
  `pow2`, `tol`, `max`, `pi` and `acos` gain `#[inline]`; the three-root arm hoists `theta / 3`,
  `2 pi / 3` and `2 r` out of its per-phase closure; `term_q` loses an `abs` it cannot need and
  `a.abs()` is bound once.
- The single-and-double-root arm's unreachability is now an executable claim: a `debug_assert!` that
  `repeated` implies `triple`, exercised by the random-bit-pattern and planted-root proptests. It is
  reachable only if the leading-coefficient floor is lowered, which `0031` (draft) L3 recommends, so
  the assert fires before the arm silently starts answering.
- `pi_and_acos_are_within_their_ulps` keeps the `atan2` twin (D6): the literal is asserted bit-equal
  to `0.0.atan2(-1.0)` at both precisions, so the evidence that justified the swap does not vanish
  into the swap.
- `dual_tests`' `ulps`/`ulps32` and their assertion wrappers come from one `macro_rules!` per width,
  so a change to the total-order key cannot reach one precision and miss the other.
- `eig3`'s three private helpers gain `#[inline]`, the convention the rest of the crate keeps; it
  also picks up `pi` as a literal, which it was paying as a `libm::atan2` per call.

### Added

- **`PHASE3.md` §8's second check: every §2.3 row against `Dual<S, DOF>` differentiation of the
  operation it is the Jacobian of** (`laws::jacobians_match_dual`), instantiated by `laws_for!`'s
  per-group arm on SO(3), SE(3), SE₂(3), `Rn`, the test-only Heisenberg group and the four product
  groups. The first check, `laws::jacobian_rows`, compares each closed form with a chain of the
  crate's *own* primitives, so a sign carried consistently through `jr`, `jl` and `Ad` satisfies it;
  this differentiates `X ⊕ τ`, `Y ⊖ X`, `X Y`, `X⁻¹`, `Exp` and `Log` in the side's own convention
  -- a group-valued row through `Log` of the correction, a tangent-valued one directly -- and so
  answers "is this matrix the derivative" from outside the closed forms entirely.

  Bounds twice the worst of 10 000--20 000 draws: `Rn` and `Rn x Rn` **exactly `0`** (every row is
  `±I` or `k I` and the operations are additions), Heisenberg 0.408, the products 0.316, SO(3)
  5.485, SE(3) 6.946, SE₂(3) 5.274. Not `0` for a group with a non-trivial `Exp`: the comparison
  carries `Log`'s conditioning, which the closed forms do not.

- `helicoid`: **`SO3::act_jacobians::<Sd>`** (`NUMERICS.md` §2.4, right `(−R [p]_×, R)`, left
  `(−[R p]_×, R)`). `API.md` §4 lists `act_jacobians` for this group and §2.4 states its rows, so it
  was the one specified item SO(3) was missing. Both groups' action rows are now checked against
  `Dual` differentiation of the action itself (`act_jacobians_differentiate_the_action`), bound 11
  from a measured 5.066 `u` over 20 000 draws -- which is not a differentiation error alone, since
  `act` is §3.3's quaternion sandwich where the closed form goes through `to_matrix`, the same
  difference §14's `act_many` row records.

### Changed

- **Every hat-structured product in the crate skips `hat`'s three structural zeros**, through
  `so3::mul_hat` and `hat_mul`: `Q`'s seven products (`NUMERICS.md` §5.3), the `W²` of `SO3::jr`
  and `jr_inv` (§3.5), and `SEn3::adjoint`'s `[x_i]_× R` (§5.2). 18 multiplications and 9 additions
  against 27 and 18. `Matrix::mul` cannot skip them itself: `0 · x` is not `0` for a non-finite
  `x`, so LLVM may not fold it and D16 forbids the fast-math that would let it.

  **Measured** against this tree at `966d715` (`f64`, `taskset`-pinned, baseline built in-tree),
  near-identity / generic / near-π: `so3/jr` **0.297 / 0.420 / 0.436**, `so3/jr_inv` 0.153 / 0.286 /
  0.326, `se3/jr` **0.165 / 0.217 / 0.225** (449 -> 125 ns), `se23/jr` 0.142 / 0.175 / 0.181
  (779 -> 167 ns), `se3/jr_inv` 0.273 / 0.301 / 0.310, `se3/adjoint` 0.370 / 0.371 / 0.371,
  `se3/rminus_jacobians` 0.459 / 0.501 / 0.508, `se23/rminus_jacobians` 0.439 / 0.463 / 0.466.
  Floors 0.003--0.13, every effect far outside; `se3/exp`, which reaches no matrix product, held at
  0.999--1.005 inside its own. With the fusion of this release, `se3/rminus_jacobians` is
  **1124 -> 340 ns**.

  **Not a rounding change.** The corpus is byte-identical -- 1559 rows, every maximum unmoved --
  which the finite-entry bit identity predicts and two tests pin: 20 000 random pairs agree exactly
  per helper, and the degenerate cases are counted rather than assumed (28 and 44 keep a `∓0` the
  dropped term would have normalized to `+0`; 216 and 225 carry a non-finite entry that `0 · x`
  would have spread across the column), with the invariant asserted that a finite difference is a
  zero's sign and nothing else. Those counts are stable across the dev and release profiles only
  under the crate's standing NaN exception -- which NaN a product returns is the implementation's
  choice of input payload, and the two instruction schedules choose differently, so on `to_bits`
  alone the non-finite counts read 216 and 234. `laws_for!`'s `dual_value_is_plain_value` arm
  already excepts NaN payloads for the same reason; these tests now do too. Every recorded law bound and `jl_is_jr_transposed_to_the_bit` hold
  unchanged.

  **The ratios exceed the arithmetic.** `SEn3::jr` saves `216 -> 144` multiplications, about 1.5x,
  and the rows moved 4--7x: a chained generic `Matrix::mul` runs at roughly 1 flop/ns where the
  written-out entries reach 3.2. `PROJECT.md` §5.1 carries that as a row -- whether the *generic*
  `3 x 3` product is worth specializing is a `helicoid-linalg` bench nobody has run, and this
  change does not touch it.

### Added

- **`PHASE3.md` §11's group benches** (`crates/helicoid/benches/groups.rs`), in a criterion target
  of their own so a change to a group does not pay for the coefficient kernel's 60 windows and the
  other way round. `bench-gate` gains `--bench <target>` and a repeatable `--only <substring>`: a
  subset is the same per-benchmark rule over fewer rows, so it decides one change, with the A/A
  control measured over that same subset -- what it cannot do is report the suite's worst, so a
  gate run claims nothing from it.
  The A/A floor over the 27 rows two pending optimisations would move was **worst 1.25%, median
  0.32%** on a `taskset`-pinned core, against the 0.6864 that one session of `0033` (draft)
  recorded. `--bless` and the report path cover a whole target, so `--only` is refused there rather
  than silently ignored, the committed samples are keyed on the target (blessing one would have
  emptied the other's directory), and the default path clears its criterion slot, which only
  mattered once a second target shared the tree. The rows then reordered both: `se3/jr` is **6.8x** `so3/jr` (445 ns against 65), and
  its near-identity row -- every coefficient on its series arm, so not one `sin_cos` -- is only
  30 ns cheaper, so the transcendentals are **~7% of that row** and Barfoot's `Q` is the cost.

- The library's own host-`std` twin, **`helicoid:host-std`**: the shipped subject's program at
  `seeded::Host`, whose transcendentals are Rust `std`'s. The envelope now reads
  `<candidate>:host-std` where it exists and falls back to the seeded stand-in's, which is what
  `bars::judge_with`'s contract already said a twin is -- "the candidate's own program with one
  variable changed". **Every domination failure is attributed: 13 / 5 / 81 / 0**, where the
  stand-in's twin left 27 with no row, because `PHASE1.md` §10's seeded kernel answers `exp`, `jr`
  and `jl` and nothing else. It also corrected rows in both directions: `so3_log`/`theta:dense` was
  called D16's cost by the stand-in and is not, four `sen3_log_n1` rows are -- the stand-in's `Log`
  is not the library's, so it was never the right program to attribute with.

### Changed

- **`SEn3`'s `rminus_jacobians`/`lminus_jacobians` fuse their two inversions**, 0.645--0.661 of the
  time at `N = 1` and 0.621--0.634 at `N = 2` -- on what `SO3::rminus_jacobians` calls the crate's
  hottest Jacobian path, one per residual per solver iteration. Both sides share `θ²`, `W`, `W²`,
  the `q_coeffs` branch, the `SO3::jr_inv` (the left diagonal block is the right one **transposed**,
  bit for bit, which `so3_tests::jl_is_jr_transposed_to_the_bit` pins over 4000 draws and
  `SO3::rminus_jacobians` has always taken) and **all seven of `Q`'s matrix products**, which is
  most of what the call costs: `Q(−ρ, −φ)`'s words are `Q(ρ, φ)`'s with two sign changes and no new
  product, so `q_words` forms them once and `q_assemble` reads the signs off a `const` side.
  **Bit-identical, and pinned as such:** `the_fused_inverses_are_the_separate_ones_to_the_bit`
  compares `to_bits` over 4002 blocks -- a random sweep plus the pure translation where `φ = 0`,
  every word of `Q` is a zero matrix and only the signs of those zeros are left -- and finds **0**
  differing. `laws::jacobian_rows` bounds the same pair at `0 u`, but that is an *error* bound and
  cannot see a signed zero, which is why the bit test exists: each side builds its own `X` with
  `hat`, because negating the other's gave the `½ρ^` diagonal `−0` where the unfused path has the
  `+0` `hat` writes. The whole corpus re-scores byte-identically -- 1559 rows, every failing row's
  maximum unmoved.
  Measured with `bench-gate --against`, the baseline binary built **in this tree from this bench
  file** with only the two method bodies reverted, so nothing but the implementation differs:
  `se3` **0.5953--0.5981**, `se23` **0.5945--0.5964**, floors 0.006--0.009. Every control --
  `so3/*`, which does not reach `Q`, and `se3/jr`/`jr_inv`, which the change does not touch -- sat
  at 0.993--1.024 inside its own floor.
  A first pass measured the baseline in a **separate git worktree** and read 0.97 on `se3/jr`,
  code that is identical in both trees: a build directory is enough to move a row 2--3%, which is
  the layout effect this module's docs already warn about. A baseline for a latency claim is built
  in the tree it is compared against.

### Added

- `helicoid`: **`SEn3<S, N>`**, with `SE3` and `SE23` (`docs/PHASE3.md` §5) — the group the whole
  `tf_tree` migration lands on. `exp`/`log` of `NUMERICS.md` §5.1, `adjoint`/`ad` of §5.2, `jr` of
  §5.3 with Barfoot's `Q`, `jr_inv` of §5.4, the six §2.3 rows, `Mul`, `Blend`, the accessors, and
  the action as `Mul<Point3<S>>` with `act_many` and `act_jacobians::<Sd>`. `jl`/`jl_inv` are not
  overridden: the provided `jr(−τ)` is §5.3's left form exactly, because negating the tangent
  negates both arguments of `Q` twice over.
  **`exp` and `log` apply `J_l(φ)` and `J_l⁻¹(φ)` by two cross products per column**, not by
  assembling `I + aW + bW²` and multiplying. `0036` (draft) measured that choice on the real corpus
  -- 16 domination failures against 36 for the assembled form on `sen3_exp_n1` -- and the shipped
  code reproduces that 16 exactly; it is also cheaper at every `N` this type has (`12N`
  multiplications against `27 + 9N`) and it is what `tf_tree_math::exp_se3` does, so the parity rows
  compare one program's rounding with its own (`0010`, `0041` draft).
  **`jr_inv` does not call `Jac::inverse`**: it composes `−A⁻¹ Q A⁻¹` with `A⁻¹` from `SO3::jr_inv`'s
  closed form, which §5.4 states and which has no division -- so `NUMERICS.md` §12's promise of a
  finite `jr_inv` holds by construction and `det A = 0` never arises on that path. That settles what
  `PHASE3.md` §0.0 parked "for the `jr_inv` PR". §14's twin backs it: `sen3_jr_inv_matches_reference`
  against the dense Gauss-Jordan inverse of `jr`, at 7 `u` -- twice the worst of 60 000 draws
  (3.437 at `N = 2` `f64`, 3.243/3.132/3.166 at `N = 1`, 3.258 at `N = 3`).
  Under every generic law of `laws.rs` at `f64`, `f32` and `Dual` for `N = 1, 2`, with bounds twice
  the worst of 60 000 draws per law per scalar: `adjoint` 8.279 `u`, `ad` 7.855, `plus_minus` 7.561,
  `axioms` 6.085 down to `sides` 1.118, and `tangent_order` and `rows` **exactly `0`**.
  `rows = 0` is an exactness claim -- every §2.3 row reproduced bit for bit, held over 200 000 cases
  at `Dual<f64, 9>` -- and it is what taking `Ad_Exp(τ)⁻¹` as **`Ad_Exp(−τ)`** buys, rather than as
  `Ad` of the group inverse: `Exp(−τ)`'s quaternion is `Exp(τ)`'s conjugate to the bit, since the
  coefficients are even in `θ` and a negation is exact, so the inverse's `N` sandwiches are saved
  *and* the row becomes bit-identical to `laws::jacobian_rows`'s own reference. Through the group
  inverse the same law read 2.182 `u` (`inverse` reaches its columns as `−Rᵗ(J_l ρ)` where
  `Ad_Exp(−τ)` reaches them as `−J_r ρ`). Review finding.
  Two more from the same review: `SE3::from_quat_translation` stated the §3.6 unit-norm domain and
  **enforced nothing** -- `SO3::from_quat_unchecked` is a move, and `Quat::from_wxyz_unchecked` is
  the item that carries that `debug_assert!`, so the constructor now goes through it; and `q_block`
  formed `W²` once per translation column, `(N − 1) · 27` multiplications for the same bits, so it
  is formed once per call and passed in.
  Scored over **all 21 `sen3_*` corpus ids** at `N = 1, 2, 3`, no non-finite output, which takes the
  envelope from 170 paired strata to **484**.
  `from_parts`/`parts` land with it (`0042`, ready): `PHASE3.md` §5 names the accessors per `N`, and
  the corpus's `sen3_log_n3` and `sen3_ad_n3` records *are* an element's parts, so without them the
  shipped code could not answer its own ids -- and making the conformance subject a privileged
  reader would retire `0006`'s external-consumer guarantee for every id at once.

### Fixed

- Two silent-mislabel hazards in the conformance subjects, both the shape `shipped.rs`'s own `So3`
  note was written about. `Sen3::Ad` built its output under the key `J` and renamed it, so a key
  change would have answered `{"Ad": []}` -- non-empty, therefore scored, therefore a committed
  conformance row holding no numbers; the field name is now a parameter of `sen3_jac_out`. And the
  Jacobian dispatch's `_ =>` arm answered any unlisted variant with `J_l⁻¹` under its own name;
  every arm is spelled, so a new variant is a build error. In the seeded twin, `so3_jac_at` read its
  side from `fn_id == "so3_jr"`, which answered every other id with `J_l`: both ids are now spelled
  and anything else answers with nothing. Review findings.

- `SEn3Jac::inverse`'s `debug_assert!` **panicked on a NaN determinant**. It read "every lane has
  `0 < |det| < ∞`"; it now reads "no lane has `det = 0` and none has `|det| = ∞`". The two spellings
  differ only on NaN, and a NaN determinant is a NaN matrix whose inverse is NaN -- the answer a
  value function owes, not a panic. `Mat3`'s `Jac::inverse` and `coeffs`'s `nonnegative` assert
  already read it that way and say why; this was the one place that did not. Found by instantiating
  `laws.rs` on `SEn3`: `dual_value_is_plain_value` draws from `f64::ANY`, so `laws::probe` inverts a
  `jr` built from garbage, which is exactly the caller the two spellings disagree about. The
  infinite case still fires, because `adj/det` is then a finite, wrong zero with a finite input to
  blame.

- The seeded subject and its host-`std` twin answer **`so3_jr` and `so3_jl`**, the two ids SO(3)
  brought that no subject but `helicoid` scored. The program is the rotation block of the SE_N(3)
  Jacobian (`se3::so3_jacobian`) and not a second reading of §3.5: it is, bit for bit and under both
  `W²` forms and both sides, the block `se3::jacobian` writes
  (`the_so3_jacobian_is_the_rotation_block_of_the_sen3_one`), so one `W²` twin and one backend twin
  serve both id families. With those rows, **every domination failure of `helicoid` is attributed**:
  5 the program's, 2 D16's, **16 neither** and **0 unattributed**, where the 16 `so3_j{r,l}` strata
  near π previously had no twin row to read (`0032`, draft, open question 1; `0037`, draft).
  They are not D16's for a measured reason, not an assumed one: of the 2466 records of each id's 28
  strata, 33 reach a `sin`/`cos` argument where the `libm` crate and glibc disagree, and **not one of
  them is on a `theta:pi-1e*` stratum** -- which is where all 16 failures are
  (`the_swap_has_no_power_near_pi_on_the_so3_jacobians`, the reading
  `the_swap_has_no_power_on_the_strata_near_pi` already records for `sen3_jl_n1`). The twin's rows on
  those strata are `helicoid`'s to the bit, so what is left is 1.3--1.7x against sophus-rs on 64
  records per stratum: a resolved comparison, and the assembly's rather than the backend's.
  `power` takes the record's input key, `phi` for these ids where the `sen3_*` ids hold `tau`.

- `helicoid`: **`SO3<S>` and `SO3Tangent<S>`** (`docs/PHASE3.md` §4), the first group with a
  non-trivial `Exp`, and the coefficient kernel's first consumer — so `mod coeffs` loses the
  `cfg(any(test, feature = "__sweep"))` it carried and a default build reaches `exp_coeffs`,
  `jr_coeffs`, `jr_inv_coeff` and `log_ratio` through it. `exp`/`log` are `NUMERICS.md` §3.1–§3.2
  (the grouped `(k, cos θ/2)` branch, the `copysign` flip, the swept `r` switch); `act`, `act_many`,
  `to_matrix`, `from_matrix` (Shepperd, §3.4, four nested `branch`es each at a safe `sqrt`
  argument), `renormalize`, `Mul`, `inverse`, `adjoint = R`, `ad = W`, `jr`, `jr_inv`, and every
  `NUMERICS.md` §2.3 row. `Jac = Mat3<S>`: SO(3)'s adjoint *is* a rotation matrix, so there is no
  structure below a dense `3 x 3`. `jl`/`jl_inv` are not overridden — the provided `jr(−τ)` is
  `I + aW + bW²` to the bit — so there is one code path and §14's twin stays the provided body.
- `helicoid`: `Side::IS_RIGHT`, the side selector `docs/PHASE3.md` §2 left to the SO(3) PR, taken
  as an associated **const** over a `TypeId` comparison: `match Sd::IS_RIGHT` resolves at
  monomorphization, so neither arm survives into the emitted code and the `'static` bound earns
  nothing. `SO3` is the first group whose two sides differ at all.
- `helicoid`: `Product::{from_parts, parts}` (`0029` option A, which specified them to land with
  this PR because SO(3) is what makes `Product<SO3, Rn<3>>` — the tf2 pose — constructible).
  Fields stay private (`0025` decision 5); before this the only route was an `Exp`/`Log` round trip
  through both factors, which for SO(3) is neither exact nor cheap and is worst-conditioned at a
  half turn.
- `helicoid`: `act_many` gains the `*_matches_reference` proptest D6 requires of its
  `NUMERICS.md` §14 twin — and it found **§14's `3 u` for that row is too tight**: the per-point
  `act` and `act_many` diverge by up to 7.587 `u` near `π`, 4.873 at `θ ~ 1` and 3.379 at
  `θ ~ 1e-6`. §3.3 prescribes forming `R(q)` once, so the nine entries round before any point is
  touched while `act` rounds a sandwich per point; the gap is the two algorithms, not a defect.
  Raising §14's figure is a normative edit and a record, recorded as owed in `PHASE3.md` §0.0.
- `cargo xtask conformance`: the `helicoid` subject answers **every `so3_*` corpus id** — `so3_exp`,
  `so3_log`, `so3_act`, `so3_from_matrix`, `so3_jr`, `so3_jl`, `so3_jr_inv`, `so3_jl_inv`, six of
  which no subject scored before. `f64`, no non-finite output, worst `max_u`: `so3_exp` 2.743
  (`theta:pi-1e-9`, equal to the seeded subject's), `so3_log` 2.627 (`theta:1e-4`), `so3_jr` and
  `so3_jl` 4.097 (`theta:pi-1e-8`), `so3_jr_inv` and `so3_jl_inv` 2.112 (`theta:dense`), `so3_act`
  4.899 (`theta:pi-1e-5`); `so3_from_matrix` is shape and finiteness only (its `max_u` is NaN by
  design). On `theta:dense` the shipped `log_ratio` scores **2.465 against the seeded reading's
  2.901**: the swept `r` switch on `s = n²/w²` beats the seeded `n² = 0` one.

### Changed

- `helicoid::coeffs`: `jr_coeffs`'s exact arm shares one `θ = sqrt z` between `a` and `b`
  (`exact_a_b`), which is what `grouped`'s contract asks and what `exact_k_cos_half` and
  `exact_b_d_e` already did. It was taking **two** square roots per call — the counting scalar of
  `a_group_runs_each_exact_arm_once` recorded `(2, 2, 0)` for the group and now records
  `(1, 2, 0)`, so every group takes one root. A common subexpression on the same `z`, so not one
  bit of any score moves; `libm::sqrt` is a single instruction only where the target maps it to
  hardware, so `just no-std` and `just wasm` were paying two software roots per `SO3::jr`.
- `helicoid`: `SO3::{rminus_jacobians, lminus_jacobians}` answer both §2.3 rows from **one**
  `jr_inv` and a transpose. `J_l(φ) = J_r(φ)ᵗ` and `J_l⁻¹(φ) = J_r⁻¹(φ)ᵗ` hold *bit for bit* —
  `hat(−φ)` is exactly `hat(φ)ᵗ`, `a` and `b` are functions of `θ²`, and `W²` as `Matrix::mul`
  forms it is exactly symmetric — pinned by `jl_is_jr_transposed_to_the_bit` over 4 000 samples.
  Each call had been paying a second `norm_sq`, a second grouped coefficient `branch` (two `sqrt`
  and two `sin_cos` on the exact arm), a second `hat` and a second 27-multiply 3×3 product for a
  matrix a transpose already held, on the crate's hottest Jacobian path.
- `helicoid`: `SO3::from_matrix` no longer panics out of domain, in debug or release, and reads
  each matrix entry once through `Matrix::get` rather than building and discarding a `Vec3` per
  scalar (the diagonal was read four times over). A matrix whose pivot candidate overflows, or an
  all-NaN one, gives a non-finite quaternion as `Quat::from_wxyz_unchecked` does; the
  normalization is written out so `Quat::from_wxyz_normalized`'s `debug_assert!` cannot turn such
  input into a panic. That matters because the method exists for locus-tag's degenerate
  near-singular matrices, and `from_matrix_does_not_panic_out_of_domain` pins it.
- `docs/maths/index.md`: the reading `so3_act` takes off the unit sphere is **measured**, not just
  listed. `SO3::act`, which is §3.3 as written, scores 510.3 `u` at the `q:nonunit` stratum and
  `to_matrix()` × v, the scaled rotation of §1, scores 257.8 — exactly `2·2^-45/u` and `2^-45/u` at
  that stratum's `η`, so the corpus reference is the **normalized** reading and neither shipped
  form is it. Adopting it would be a §3.3 edit and a record; until then the stratum is outside
  `act`'s unit-`q` domain and that row is the open question, not a defect.
- `helicoid::coeffs`: `nonnegative` now reads "no lane is negative" rather than "every lane is
  `>= 0`". The two differ only on NaN, which the assert is not there to reject — a *negative* `θ²`
  is a sign error worth a panic, a NaN one is a NaN input, and every arm returns NaN for it, which
  is the answer a value function owes. Reached by `laws::dual_value_is_plain_value`, which drives
  `SO3::exp` over `f64::ANY`. `SO3::exp` likewise builds its quaternion field by field instead of
  through `Quat::from_wxyz_unchecked`, whose `debug_assert!` is for a caller *claiming* unit and
  would turn a NaN tangent into a panic; `Exp`'s output is unit by construction.

- `docs/decisions/0036` (draft): partitions the envelope's 75 domination failures by *which* oracle
  won, which turns 67 unexplained into 30 worth an experiment and 45 worth a confirmation.
  `tf_tree_math` routes through the `libm` crate as we do, so the 30 it beats us on cannot be D16's
  cost; sophus-rs is glibc, so the 45 it alone beats us on are `0032`'s class. 29 of the 30 are
  `sen3_exp_n1`, and the cause is that `NUMERICS.md` §5.1 does not say how `J_l(φ) ρ` is *applied*:
  forming `I + aW + bW·W` as a matrix costs 36 failures, two cross products cost 16, and the
  `W² = φφᵀ − θ²I` identity costs 19, with every other id untouched in all three runs. No form
  dominates -- cross products win away from π (worst 1.5087 u against 2.1517) and lose near it
  (4.6856 against 3.2855) -- so the record proposes the §5.1 addition and leaves the choice between
  accepting the trade and sweeping the regime. Measured on the seeded stand-in, whose readings differ
  from the shipped kernel's in documented ways; `PHASE3.md` §10 owes the real subject.
- `cargo xtask bench-gate --against <binary> --record <dir>` and `--replay <dir>`, which split the
  gate's twenty-minute measurement from the arithmetic it feeds. The decision rule -- median, paired
  bootstrap CI, swap-invariant `δ`, every-pair-above-the-floor -- is a pure, seeded function of
  sample vectors, so validating it never needed a quiet machine; it needed data. `--record` copies
  every window criterion writes, `--replay` runs the identical rule and report over a recording and
  measures nothing, and the contamination properties are now ordinary tests: one disturbed window, a
  disturbance spanning one triplet, uniform drift, and a real slowdown. The 30-second independence
  failure that cost a twenty-minute run to find is four lines (`0035`, draft).
  The stored format is criterion's own `sample.json`, copied byte for byte: lossless by construction,
  no second parser, diffable against `target/criterion`, and a replay goes through
  `samples::read_file` so the `times`/`iters` division is replayed rather than assumed.
- `docs/decisions/0034` (draft): D16 states an implementation and claims a property, and the property
  holds only per *resolved* `libm` version -- the workspace declares `libm = "0.2"`, correctly for a
  library, so a patch release may legally move bits, and nothing records which version produced a
  committed baseline. The seam D16 names already exists (six `libm::` sites in `helicoid-linalg`'s
  `float.rs`, none in `Dual`, `sincos` so there is one range reduction), so what is proposed is a
  three-clause reword, a containment lint, a recorded version, and a trigger for a second backend:
  end-to-end speedup is `1/(1 - f + f/k)`, which is 1.08x at `f = 0.15`, `k = 2` -- below the bar
  tf_tree's `0016` used to reject a 12% gain -- and the only profile in the stack shows no `libm` row.
  Build the guard, not the knob. A bit change is a minor-version bump, by the owner's call, which is
  what makes a faster portable kernel available as a step at all.
- The bench harness and `cargo xtask bench-gate` (`docs/PHASE1.md` §9, `docs/PHASE3.md` §11's
  owed-first row), which did not exist: `criterion` pinned at `=0.8.2` in `crates/helicoid/benches`
  over the shipped coefficient kernel through `__sweep` (so the code timed is the code that ships),
  60 benchmarks by stratum -- §9's near-identity, generic and near-π plus the two switch
  neighbourhoods, since a group's cost is which arm ran. `just bench` runs it: `--against` gates,
  `--aa` measures this host's A/A noise, `--bless` writes what each mode measured, and the bare
  form reports against the committed samples.
- `cargo xtask bench-gate --against <bench-binary>`: §9's gate. Each benchmark runs **baseline,
  candidate, baseline**, three times, so the A/A floor is a control measured *beside* the
  comparison, on the same benchmark, in the same invocation -- and a benchmark fails only when all
  six candidate pairs put their whole CI above `1 + floor`. Validated by both controls: the
  identical binary as its own baseline gives **0 of 60** failures, and slowing
  `coeffs::kernel::exact_a` fails **5 of the 5** benchmarks that evaluate it at 1.38x to 1.57x.
- `baseline/HOST.md`: a **log** of what this host's A/A noise has shown, written by `--aa --bless`.
  It is not the gate's allowance, and `--against` does not read it. The samples under
  `baseline/bench/` are *not* committed (`.gitignore`): an accuracy maximum is a function of the
  source and the corpus alone (D16) and travels, a timing is one core under one load and does not.
- `docs/decisions/0033` (draft): the estimator choices §9 leaves open, each settled by a
  measurement that refuted the first guess, and the §9 deviation the gate now carries -- §9 records
  `δ` in `HOST.md`, and this host says a stored `δ` cannot be an allowance.

### Changed

- **A stored noise floor is not an allowance on this host, and the gate no longer uses one.** The
  same `--aa` protocol, the same binary against itself, measured a worst `δ` of 0.0161 in one
  session and **0.6864** a few hours later -- 23 of 60 benchmarks above 2% where the first run had
  none. A floor from the quiet session fails identical code (8 of 60); one from the noisy session
  would pass a 50% regression. Nothing in the protocol changed between them, so what moved is the
  machine. Per-benchmark floors, and an accumulated maximum over runs, were both implemented and
  both refuted by measurement before the concurrent control replaced them (`0033`, draft).
- **The cross-run comparison against committed samples is reported, never gated.** A stored
  baseline is not adjacent in time by construction: identical code against samples blessed earlier
  the same day read 0.49x to 2.04x, 27 of 60 past their floor. `--against` is what gates.
- **A gate failure says the two binaries differ at that benchmark, not that the change caused it.**
  The negative control also failed four `exp_coeffs` benchmarks at 1.016x to 1.041x, which have no
  path to `exact_a` -- it has one caller. Between the two binaries 2 of 3192 text symbols changed
  size, `.text` grew 64 bytes and 2997 symbols kept an identical size at a different address:
  byte-identical code, relocated. Source locality does not survive to machine level.
- `docs/decisions/0032` (draft): `so3_log`'s 8 domination failures are not an algorithm. `Log` is the
  same program in `helicoid` and both oracles (D5, the quaternion `atan2`); the gap is the `libm`
  crate's `atan2` against the host's glibc, which D16 buys and `error-analysis.md` EA.13(d)
  predicted ("a domination failure there can be the library's and not the algorithm's", `libm`'s
  nu 1.96 on `atan2` against glibc's <= 1.00). Measured over the corpus with one variable changed:
  the two `atan2`s disagree on 692 of 4994 records, the `libm` column reproduces the harness's
  `seeded:correct` on 8 of 8 losing strata and the glibc column reproduces sophus-rs's on 7 of 8,
  and `tf_tree_math` -- also on the `libm` crate -- is bit-identical to `helicoid` on all 30 strata.
  D16's price is at most 0.74 u (ratios 1.08 to 1.45). Proposes a third `libm-bound` verdict so the
  bar stops charging the candidate for choosing reproducibility; notes the other 67 failures are
  *not* explained by this, since `Exp` goes through `sin_cos` where the penalty is ~1.1x not ~2x.
  Four open questions, no code change.
- `docs/decisions/0030` (draft): `chol` factorises 2.6 to 3.4x slower than `nalgebra`'s at `N = 6`
  and 1.3 to 1.5x at `N = 3, 8`, measured advisorily beside the two solves `locus-tag` runs today;
  `chol_solve` is at parity and the residual favours `chol`. The attribution did not settle: a
  verbatim-copy control of the measured body itself read −0.4% to +11.3%, so the harness cannot
  adjudicate `0021`'s 15% bar. A bit-identical indexed reduction is nonetheless a live candidate
  (whole ten-run range below zero, midpoint near −20%) and is deliberately **not** taken, because
  the sweep tested one of the bar's three clauses — the `f64` success path — and not failing input
  or `Dual`, where `0021`'s own row-major variants lost 38% to 41%. Proposes an external subject as
  a `0021` rule-3 trigger, `PHASE1.md` §9's `criterion`/`bench-gate` harness extended to
  `helicoid-linalg`, and comparison rows in a workspace-excluded `runners/` crate. Five open
  questions, no code change.
- `docs/decisions/0022` (**ready**, retitled *`Real` owes `acos` and `cos`*): the private
  `atan2(sqrt((1 - x)(1 + x)), x)` costs 2.61x `libm::acos` (11.42 ns against 4.37) and is
  marginally *less* accurate everywhere measured, the near-`±1` region it was chosen for included
  (1.28 u against 0.92 as `x -> 1`), and it forfeits bit-identity with omnisac in the trigonometric
  arm. `solve_cubic` and `eig3` are both callers — `eig3`'s own PR amended this record to say so —
  so "public surface for one caller" does not hold. `Real::cos` saves a further 26% per call over a
  `sin_cos` whose sine is discarded. `pi` as a literal has landed; the two methods are step 1.
- `docs/decisions/0031` (draft): the four numerical limits the cubic port inherits from omnisac,
  each measured, with a recommended order. The one-real-root cancellation first — `x^3 + p x - 1` at
  `p = 1e-5` returns exactly `1.0` for a root of `0.999996666666666679`, a relative error of 3.33e-6,
  which the `u v = -p/3` pairing takes to 3.42e-17 and which also removes a `-inf` `Dual` derivative
  at a simple root; it dominates on the max (4.67e-10 to 2.63e-12 over 16 010 random one-root
  cubics) but is *not* uniformly better per row. Then the underflow that reports three valid slots
  for a one-root cubic, the leading floor that rejects `(x - 100)(x - 101)(x + 99)` at `f32`, the
  discriminant band measured against the wrong quantity, and the `1/a` reciprocal. Nothing decided.
- `docs/decisions/0028` (draft): `SEn3Jac`'s `diag`/`col` are `pub` in `PHASE3.md` §5 and private
  under `0025` decision 5, and `ProductJac` took the other reading; three options, their costs, and a
  recommendation (narrow, the reversible direction). Nothing decided, no code change.
- `docs/decisions/0029` (draft): `Product` has private fields and no constructor, so the tf2 pose
  `Product<SO3, Rn<3>>` is reachable only through an `Exp`/`Log` round trip that is inexact, costs an
  `atan2` and a `sin_cos`, and is worst-conditioned at a half turn; `from_parts`/`parts` proposed, to
  land with the SO(3) PR. Nothing decided, no code change.
- `just lint` runs `cargo clippy --workspace --all-targets --release -- -D warnings` as well as the
  dev pass: `debug_assertions` is off in release, so the `cfg(not(debug_assertions))` test bodies that
  only `just test`'s release run executes were never checked against the denied `unwrap_used`,
  `expect_used`, `panic`, `todo`, `unimplemented` and `dbg_macro`. Clean today, so this closes a gap
  rather than fixing a violation.

||||||| parent of 02e5b6e (docs(decisions): 0026 (draft) write_dense pays for zeros, not for checks)
- `docs/decisions/0026` (draft): `Jac::write_dense` pays for its structural zeros, not for its
  per-entry bounds checks. The 81 branches and 81 panic calls at `DOF = 9` are never-taken, cold
  branches and the loop already runs at about one store per cycle; two check-removing designs
  measured 53% to 326% slower, while bulk-zeroing the view and writing only the non-zeros wins
  36% to 63% at `DOF` 6 and 9 on a column-major destination and loses at `DOF` 3 and on
  row-major. No API change: the path is instrumented first, and adoption is deferred to the
  Jacobians a solver actually calls. Documentation only, plus one rustdoc line.
- `docs/decisions/0025` (ready): the trait layer as built. Rⁿ's Jacobian is the structured `RnJac`,
  correcting `PHASE3.md` §7's `Jac = Mat<N>`, which cannot implement the exact `Jac::inverse` of
  `0005`; `Side` is sealed to `Left`/`Right` while the side selector stays with the SO(3) PR;
  `Tangent::dot_acc` becomes the required operation so that `Product` can meet the bound-`0`
  dense-order law; `read_dense` poisons a missing entry; and a group exposes an `Add`-carrying field
  only where it is abelian.
- `docs/decisions/0020`: `Dual::sqrt` at zero keeps `d / (2 sqrt v)` (the NaN is the report and the
  `PHASE1.md` §10 row needs it); the guarded norm is a doctest, a zero-safe norm is a later record.
- `docs/decisions/0021`: the per-entry finite guard of `chol` stays; four bit-identical variants were
  measured, none clears the stated 15% bar, and `L` stays finite for every input. Documentation only;
  no code change.
- Tests pinning today's `Dual` derivative at a zero `sqrt` argument: `Vector<Dual>::norm` of the
  zero vector, and `chol` on a zero or a negative pivot.
- The `helicoid` conformance subject (`docs/PHASE1.md` §5, `xtask/src/shipped.rs`): the shipped
  groups of `helicoid::coeffs` through `__sweep`, over the eight `coeff_*` ids at `f64` and `f32`,
  `value` and `d_branch` on `Dual<S, 1>`, registered as a plain subject, so `just conformance`
  prints and writes its per-`(fn, stratum, precision)` rows beside `seeded:correct`'s. It equals
  `seeded:correct` bit for bit on every record (the two share every candidate) and scores finite
  everywhere; a four-terms-below-`z = 0.01` prior (D12 for `a`, `b`, `c`; not D12's for `r`, `k`,
  `d`, `e`, `cos θ/2`: 0014 (draft) questions 10, 29) has a smaller maximum on nine `cos θ/2`
  strata, pinned by a test and raised as 0014 question 30. No baseline, no envelope. Tooling
  only; nothing breaks.
- The coefficient kernel (`docs/PHASE3.md` §3, `0004`): `helicoid::coeffs`, `pub(crate)`, with the
  groups `exp_coeffs`, `jr_coeffs`, `jr_inv_coeff`, `q_coeffs` and `log_ratio`, each one `S::branch`
  on "every member is on its series arm", off which the group's exact arms run once, at the safe
  argument, and each member selects by its own mask; the exact arms in the operand order of the
  seeded kernels, the series by Horner, the constants chosen by `S::PRECISION`. `cargo xtask
  thresholds helicoid` sweeps the shipped arms through the hidden feature `__sweep` (`exact_*`,
  `series_*`; only `xtask` enables it) and writes `crates/helicoid/src/coeffs/generated.rs` (a
  `Switch` per coefficient per precision, beside the eight swept series they take their terms from;
  registered with `cargo xtask lint`) and `conformance/sweeps/thresholds.csv`, a fixed point: a
  second run writes the same bytes, and `just thresholds-check` covers both this target and the
  seeded one (`cargo xtask thresholds [--check] [seeded|helicoid]`). At `f64`: `k` 8 terms below `θ²
  = 1`, `a` 8 below `0.72`, `b` 8 below `0.96`, `c` 8 below `0.56`, `d` 8 below `1`, `e` 7 below
  `0.96`, `cos θ/2` 2 below `5.6e-15`, `r` 8 below `n²/w² = 1.1e-2`; at `f32`: `k, a, b` 5 terms and
  `c, d, e` 4, all below `θ² = 1`, `cos θ/2` 2 below `2.6e-6`, `r` 8 below `0.143`. They are the
  seeded sweep's rows, and the shipped arms are the seeded arms bit for bit on every record where
  the series arm is defined (`w > 0` for `r`). Each objective is a maximum over the corpus's
  records, not over the domain, and eight of the sixteen switches (`k`, `d` at `f64`; `k`, `a`…`e`
  at `f32`) are `θ² = 1`, the top of the sweep's grid, not a measured optimum. The sweep CSVs gain
  four columns, each arm's error at the two records that bracket the switch, whose sum
  `branch_continuity_*` compares the jump between the arms with (`docs/maths/coefficients.md` CO.12;
  a sample at two records, not a bound over the interval; not `NUMERICS.md` §4's wording); the
  seeded CSV and the digest line of `xtask/src/seeded/generated.rs` change with them, its constants
  do not. Groups take their members' switches each (0015 (draft) NU.5): no exact arm runs below a
  group's smallest switch, which for `Exp` is `cos θ/2`'s (`θ² < 5.6e-15` at `f64`), so `Exp` is on
  its exact arm from `θ = 7.5e-8` at one `sqrt` and one `sin_cos` (`NUMERICS.md` §3.1); unmeasured.
  `r` is on its provisional reading (0014 (draft) question 29); `se2_coeffs` and `gamma2_coeffs` are
  not here. The kernel is built for the tests and `__sweep` only until SO(3) uses it, so `just
  no-std` and `just wasm` build it with `--features __sweep`. No public API; nothing breaks.
- Threshold sweep per precision (`0016` item 3): `cargo xtask thresholds` now sweeps `k, a, b, c, d, e`,
  `cos θ/2` and `r` at `f64` on the plain strata and at `f32` on the `@f32` strata, with one objective,
  grid (each point correctly rounded at the precision it is swept at) and tie-break. At `f32` each of
  `k`…`e` sits at the top of the grid (`θ = 1`), its optimum at or above it, as `docs/maths` CO.10
  predicts; `cos θ/2` sits at the exact arm's floor at both precisions, its switch the tie-break's.
  `r` is swept on a reading pending the maintainer, not a spec: branch variable `s = n²/w²`, only
  records with `w > 0`, series arm iff `w > 0` and `s < switch`, grid over `s` from `1e-16` to `1`,
  its prior `(4 terms, s < 0.01)` and not D12's (0014 (draft) question 29). The seeded sweep's
  outputs are renamed `conformance/sweeps/thresholds-seeded.csv` (16 rows) and
  `xtask/src/seeded/generated.rs` (a `Switch` per id and precision, the `f32` literals correctly
  rounded from the exact rationals, a test in integers), leaving `thresholds.csv` to
  `helicoid::coeffs`. `seeded:correct` runs the `f32` constants (where it ran the D12 prior;
  `subject_version` `generated`) and now also answers `coeff_cos_half` and `coeff_r`; the
  self-test's pinned `f32` figures follow (the correct kernel's `b` curve reads `p = -0.021`).
  `xtask` only; the earlier `f64` rows are unchanged; breaks nothing.
- Harness: `Precision::F32` scoring (`0016` item 2). `cargo xtask conformance --precision f32` scores
  `f32` subjects on the `@f32` strata alone, in units of `2^-24`, with inputs and outputs asserted
  exactly binary32. An id named with `--fn` that has no `@f32` strata is an error, not a skip; without
  `--fn` the run scores the eight ids that have them and names the others a selected subject supports.
  A subject with no `f32` kernel (the planted `c`, a sweep candidate) is an error at `f32`, never the
  D12 answer under its name. The seeded coefficient kernels run at `f32` through the same adapter (the
  D12 prior, not a correct kernel, until the `f32` sweep of `0016` item 3); the correctly rounded and
  neighbouring-ulp sanity subjects have `f32` variants; `--self-test` runs the error-curve and
  non-finite mechanisms at both precisions (`b` without its series fits p = 1.937 at binary64, 2.134
  at binary32; the `f32` range, kernel and floor are questions 26 to 28 of draft record 0014).
  `xtask` only; nothing breaks.
- Corpus: `@f32` strata for `coeff_k`, `coeff_a`…`coeff_e` and `coeff_r` (`0016` item 1: the binary64
  inputs rounded to nearest-even binary32, the reference recomputed at the rounded inputs), and the
  id `coeff_cos_half` (cos θ/2) with its series as the eighth row of `coeff_series`; `0016` item 1
  now names `coeff_cos_half` and `coeff_r` (`PHASE1.md` §4.3, §4.4 follow). Existing records are
  byte-identical; the harness skips `@f32` strata and the binary64 sweep excludes them until `0016`
  items 2 and 3.
- Decision records `0016` (`f32`-exact `@f32` strata for the scalar coefficient ids) and `0017`
  (`Real::cbrt`, mask-valued `solve_cubic` roots), both `ready`: docs only, no code change.
- Workspace skeleton: `helicoid-linalg` and `helicoid` (empty, `no_std`, `forbid(unsafe_code)`),
  `xtask` stub, lints, `justfile`, CI, `deny.toml`. No group code (`docs/PHASE1.md`).
- `docs/maths/`: non-normative derivations. Notation, labelling and verification convention, and the
  Lie-group foundations (`Ad`, `ad`, the Jacobians of `Exp`, every row of `NUMERICS.md` §2.3, a sign
  audit with exact gaps, the conditioning of `J` and `Ad`). Documentation only; no code, no formula
  change.
- `docs/maths/so3.md`: SO(3) derivations (quaternion action and double cover, `Exp`, `Log` and the
  conditioning of the `acos` form it avoids, `J` and `J⁻¹` solved in the algebra of `{I, W, W²}`,
  `Ad`/`ad`, the action Jacobians, Shepperd's `from_matrix` and the nearest rotation, the Newton
  renormalization and its bounds). Documentation only; no code, no formula change. Five open items
  for the normative documents are recorded in the maths index.
- `docs/maths/coefficients.md`: the coefficient catalogue of `NUMERICS.md` §4 (general terms and radii
  of convergence, the identities between `k, a, b, c, d, e, r` and the `cos(θ/2)` of §3.1, rounding
  error of the exact arms and of Horner, truncation bounds, the crossing that predicts switch-point
  magnitudes, derivatives through `Dual`, the safe argument and the underflow, `π` and `w → 0`
  cases, the reference precision budget, the cost of one shared switch per call-site group).
  Documentation only; no code, no formula change. Nine open items for the normative documents are
  recorded in the maths index.
- `docs/maths/se3.md`: SE_N(3) derivations (the matrix group and the `x_1 = v`, `x_2 = p` convention of
  SE₂(3), `Exp`/`Log` from the hat matrix's series, `Ad`/`ad`, Barfoot's `Q` block derived from
  `J_l = Σ adⁿ/(n+1)!` by the block-triangular powers of `ad` and `W³ = −θ²W` and identified with `b, d, e`,
  `J_r` and two block relations, the dual-matrix algebra with proofs of closure, product (27 + 54N
  multiplications), inverse and `apply_transpose`, `J⁻¹` without a second closed form, the SE(3) action
  Jacobians, the rotation-first/translation-first permutation of tangents, Jacobians and covariances,
  what the exact arms of `b, d, e` cost `Exp` and `Q` (and the error of `Q` at a switch), and the rounding of
  `Log`'s translation part and of the blocks of `J⁻¹`). Documentation only; no code, no formula
  change. One open item for the normative documents is extended and one added in the maths index.
- `docs/maths/so2-se2.md`: SO(2) and SE(2) derivations (unit-complex `SO2`, `atan2` sensitivity; `Exp` with
  `V = αI + βK = s(θ)R(θ/2)`, `Log` and the three forms of `V⁻¹`, the two preimages at `|θ| = π`; `Ad`,
  `ad` and `ad³ = −θ²ad`; the rotation-first `J_r`, `J_l` and their inverses derived as `I ± a·ad + b·ad²`
  and `I ∓ ½ad + c·ad²`, with the dense entries written out and marked **Proposed** for `NUMERICS.md` §6
  and §2.4; the reduction of `Q` on the plane and agreement with Solà et al. App. C through the
  translation-first permutation; the scalar functions with series and singularities; rounding, the
  conditioning of `J` and a sign audit). Documentation only; no code, no formula change. Four open items
  for the normative documents are added in the maths index.
- `docs/maths/error-analysis.md`: the rounding model (`u = 2⁻⁵³`, `2⁻²⁴`; `libm` is not correctly rounded, measured),
  forward and backward error and the metric of `NUMERICS.md` §11 (chord, floors, what a norm-wise metric cannot
  see, `𝓑` from the stored reference and where it holds, why Jacobians are judged by forward error); the
  conditioning of `Log` (norm-wise `1/sin(θ/2)`, componentwise `≤ √5`, the cut), what backward error adds near π
  (branch invariance), `from_matrix` (uniform `1/(2√2)`, the backward-error floor `d_M`) and `J⁻¹` of SO(3) near
  `2π`; why bars are max and p99 and never a mean, non-finite counts, domination, exact no-regress and what D16
  buys (bit reproducibility across x86_64 and wasm32 with `libm`'s default features off, measured; not accuracy); forward-mode `Dual` (the
  rules, the value-path theorem, what the derivative of a branch or a series arm is, second order by nesting, the
  `sqrt`-at-0 hazard and the safe argument); the mpmath corpus (exact inputs, 30 digits, precision budget and the
  `dps = 150` recheck, log-uniform decades, uniform axes on S² by Archimedes, splitmix64 and per-stratum streams).
  Documentation only; no code, no formula change. Six open items for the normative documents are added in the maths index and one extended, and `coefficients.md`'s remark on the `libm` crate is corrected.
- `docs/maths/geodesics.md`: the geodesic of `NUMERICS.md` §10 (the curve, its constant body and spatial velocity,
  what "geodesic" claims on SE_N(3): the autoparallel of the canonical connection and, on SE(3) only, the geodesic of its indefinite invariant forms (no nondegenerate one for N >= 2, no Riemannian one);
  left-, right- and inversion-invariance; `Product<SO3, R3>` as (slerp, lerp), bi-invariant for its own law and
  failing right-invariance only when composed as SE(3) poses, with the exact gap; the Jacobians with respect to
  `X0`, `X1` and `t` derived from those of `Exp`, their left forms and cancellation-free forms of `J0` and its left form; unit
  dual quaternions, the screw parameters, the power, its sign rule and the grouped one-`atan2` evaluation of
  `tf_tree_math`'s `screw_pow`, with its equivalence to the geodesic; where the twin cancels or branches (no
  cancellation in `f64` value at small angle or near π, the per-type 0/0 guard, the `Dual` translation derivative that loses about 1/α whatever series is used, the endpoints, the cut at π); SO(3) slerp).
  Documentation only; no code, no formula change, no normative document changed. Five open items for the
  normative documents are added in the maths index.
- `docs/maths/charts.md`: charts of `PHASE5.md` §1–§2 (a chart as a frozen retraction with a local inverse and
  its two Jacobians read in the chart of the retracted point, `lj(ret δ) = rj(δ)⁻¹`; `RightChart` and `LeftChart` from
  the rows of `NUMERICS.md` §2.3; `Screw`, `Decoupled` and `WorldTranslation` as two retractions in three frames, with
  the block-diagonal `rj`, `lj` of the last two derived and why they are not dual matrices; the change-of-chart
  identity, which contains `se3.md` SE.9(b), (c); what a chart change alters in an LM iterate: not the step of `Screw`
  and `Decoupled`, but the trial point by at most `θ‖ρ‖/2`; S²: the Householder basis and its closed form, `retract`
  as the sphere's exponential of the rotation vector, `local` and the `r` kernel, `rj = ZᵀA(δ)`, `lj = A(δ)⁻¹Z`, the
  reflection at `n_z = 0`, the hairy-ball argument for freezing the chart, the antipode). Documentation only; no
  code, no formula change, no normative document changed. Four open items for the normative documents are added
  in the maths index (`PHASE5.md` §1.3's first-order claim is false for `WorldTranslation` unless `R = I`; `sgn(-0)`;
  `mp.logm` as the `s2_local` reference; the unit-norm domain of `S2`) and two extended.
- `docs/maths/gamma-gaussian.md`: `NUMERICS.md` §7 and `PHASE5.md` §4–§5 (the integrated exponentials `Γ_m` with their closed forms, recursion and integral form, derived, and what
  the exact arms cost `Γ_1` and `Γ_2`; the piecewise-constant IMU increments `ΔR`, `Δv`, `Δp` derived from the ODE under stated assumptions, with the `5×5` exponential and the composition
  rule; the directional Jacobians of `Γ_m v` in closed form (`m = 1` is `Q − [J_l ρ]× J_l`; `m = 2` needs a coefficient outside §4) and what the `Dual<S, 3>` path returns, exactly and
  measured, Taylor branches included; `Gaussian` as a mean, covariance and side, the change of side `Σ_L = Ad Σ_R Adᵀ` proved **exact**, the order of first-order propagation and
  composition, the Mahalanobis distance; what `chol`'s mask certifies and the rounding of `d²`, both set by the correlation matrix's smallest eigenvalue; `Ad` of SE(3) at a large
  translation (its conditioning by the unit-free correlation matrix, the round-trip loss), and seven misread covariances quantified). Documentation only; no code, no formula change, no normative document changed. Four open items for the normative
  documents are added in the maths index and one extended.
- `docs/maths/ambient-jacobians.md`: `PHASE6.md` §1 (over-parameterized storage and the ambient Jacobians of a chart, `PlusJacobian` = `D ret(0)` and `MinusJacobian` = `D loc(x)`; the quaternion's `P` (4×3) and
  `M` (3×4) derived from the Hamilton product, the first-order `Exp` and the scale invariance of `Log`; `MP = ‖q‖²I`, `PM = ‖q‖²I − qqᵀ`, the orthogonal projector onto the tangent space of the unit sphere for unit `q`,
  `M = 4Pᵀ = P⁺`; what a non-unit stored `q` does (`P` exact, `M` off by `η`, three readings of a reference, the drift of `η`: a random walk for independent steps, linear for a repeated one, the double cover); SE(3) as `(q, t)`: `diag(P, R)` and `diag(M, Rᵀ)` for `Screw` and
  `Decoupled`, proved equal to first order, the dependence on `R`, `WorldTranslation` and the left chart; Ceres' `QuaternionManifold` shown to be the left chart in half-angle coordinates, exactly (`P_C = 2PRᵀ`,
  `M_C = ½RM`; its `Minus` does not flip the sign), checked against the real library through `pyceres`; GTSAM's retraction on storage (default and option-off modes) and why it needs no ambient Jacobian; conditioning and rounding; a sign audit with
  exact gaps). Documentation only; no code, no formula change, no normative document changed, and no phase status table advances. One open item for the normative documents is added in the maths index (`PHASE6.md`
  §1, §7: the non-unit domain, the corpus reading, the missing `MinusJacobian` id and §14 twin row, the Ceres convention, the implementing types) and two are extended.
- `docs/maths/sim3.md`: Sim(3) (`NUMERICS.md` §9, owed). `Exp`, `Log`, the adjoints, the dense `7×7` `J_l`, `J_r` and the action Jacobians, derived; the coefficients of `𝖵 = v₀I + v₁W + v₂W²` in closed form
  (denominator `Δ = σ² + θ²`), as a two-variable series and as an exact split into functions of `σ` alone; the rounding at the joint limit, the overflow onsets and a scaled arrangement that holds to `ln MAX`.
  The page's own Results table carries the numbers. Every formula §9 does not state is marked proposed, to be adopted by a record; `PHASE5.md` §3 does not start before that. Documentation only; no code, no
  formula change, no normative document changed, and no phase status table advances. The open items it raises for the normative documents are in the maths index.
- `cargo xtask lint`, run by `just lint`: no line citations in markdown or comments, no draft
  decision record cited as settled (status-table rows, Rust comments, amendment banners; not
  prose elsewhere), and `@generated` files registered with their owning task and header. Owed:
  twin table (`docs/PHASE1.md`).
- `cargo xtask lint` reads `cargo metadata --locked` (all features): the normal-dependency
  closure of `helicoid-linalg` and `helicoid` must equal the `0007` set (`mint` only as an
  optional direct dependency, never enabled by `default`), reporting the path that brings in an
  offender; no crate but `xtask` may enable `__sweep`. `xtask` gains `serde` and `serde_json`;
  `deny.toml` allows `Unicode-3.0` for `unicode-ident` and bans the `0007` list, which
  `just audit` checks for every crate but `xtask`.
- `helicoid-linalg`: the scalar model (`docs/PHASE2.md` §2): `Mask`, `Real`, `Blend`, `Precision`;
  `f64`/`f32` as `Real` with every transcendental through `libm`; `Blend` for scalars, tuples up
  to arity 8 and arrays. `Real` has no `PartialOrd`/`PartialEq` (`compile_fail` doctests). New
  public API, nothing breaks.
- `helicoid-linalg`: forward-mode `Dual<S, N>` (`docs/PHASE2.md` §3), a `Real` and a `Blend<S>`
  with every derivative rule of §3; its value path is bitwise the plain `S` evaluation, NaN sign
  and payload of arithmetic outputs excepted. `copysign` takes `sgn(s)` from the sign bit. The
  derivative of `sqrt` at 0 is infinite or NaN, that of `atan2` is inexact once `x^2 + y^2` leaves
  the normal range, and the quotient's is NaN once `a / b` overflows (`NUMERICS.md` §12). New
  public API, nothing breaks.
- `helicoid-linalg`: `Vector<S, N>`, `Point<S, N>` and column-major `Matrix<S, R, C>` with the
  aliases `Vec2`/`Vec3`/`Mat2`/`Mat3`/`Point2`/`Point3` (`docs/PHASE2.md` §4): sums, `scale`,
  `dot`, `cross`, `norm`, generic matrix product, `transpose`, `identity`, `from_cols`,
  `from_rows`, `col`, `row`, `get`, `set`; `Point - Point` is a `Vector`, `Point + Vector` a
  `Point`. `hat`/`vee` and `Mat3::inverse_adj` (adjugate over determinant, plus the determinant).
  `Blend` for all three. No `PartialEq`. Every reduction sums left to right. New public API,
  nothing breaks. `proptest` is a dev-dependency only.
- `helicoid-linalg`: `chol` (lower-triangular `L` and a positive-definiteness `S::Mask`, set iff
  every computed pivot is strictly positive and every computed entry finite), `solve_lower` and
  `solve_upper` (`docs/PHASE2.md` §4; maths in `docs/NUMERICS.md` §15). Branch-free: `sqrt` gets a
  safe argument, a failed pivot gives `L_jj = 1` and a zero column, an overflowing entry is stored
  as `+0` and clears the mask, so `L` is finite for every input. The solves `debug_assert!` a
  nonzero diagonal. New public API, nothing breaks.
- `helicoid-linalg`: `Strided` and `StridedMut`, read-only and writable strided views over caller
  memory (`docs/PHASE2.md` §5): `col_major`, `row_major`, `with_strides` (faer and Ceres layouts),
  `block`, `get`, `set`. An out-of-bounds access panics in release, the index never wraps, and a
  write never leaves its view (D11's one panic class, widened from writes to `get`, `set` and
  `block` in `PROJECT.md` D11 and `NUMERICS.md` §12). New public API, nothing breaks.
- `helicoid-linalg`: optional feature `mint` (off by default; `mint` >= 0.5.7 joins the closure
  only with it): `From`/`Into` between `Vector<S, N>` (`N` = 2..=4), `Point<S, N>` (`N` = 2, 3; `mint`
  has no 4-D point) and `Matrix<S, N, N>` (`N` = 2..=4, `mint::ColumnMatrixN`), for `f32` and
  `f64` (`docs/PHASE2.md` §7). Scalars are moved, never computed on. `just lint` and `just test` cover the feature; `just msrv`, `just no-std`
  and `just wasm` build it. Bitwise round trips pin `-0`, infinities, subnormals and NaN payloads. New public API, nothing breaks.
- `docs/decisions/0018` (implemented): `libm`'s `arch` feature is bit-identical for exactly-rounded operations
  (dispatch table read for x86_64, aarch64, wasm32, thumbv7em; 95 digests over the `libm` entry points
  equal on x86_64, wasm32 (node WASI) and aarch64 (`qemu-user`); NaN sign and payload stay the target's,
  outside D16).
- Reference generator skeleton in `conformance/generate/` (uv, pinned CPython 3.12 and mpmath;
  splitmix64, per-stratum streams, `MANIFEST.json`, 150-digit recheck) and the first corpus file,
  `coeff_k`, all scalar-θ strata. New recipes `just corpus`, `corpus-check`, `corpus-test` and a
  `corpus-check` CI job. No library code changes; breaks nothing.
- Corpus files `coeff_a`…`coeff_e` and `coeff_r` (the θ strata, plus `q:w0` at three norms for `r`)
  and `coeff_series`, the exact 16-term rational Taylor series of every `NUMERICS.md` §4
  coefficient (manifest `kind: "series"`, `verified` in place of `rechecked`; `PHASE1.md` §4.3).
  Every record is cross-checked at generation against the series and a second formulation;
  cancelling definitions carry guard digits, so `theta:subnormal` is evaluated correctly. Manifest
  entries gain `kind`; `coeff_k` is unchanged. No library code changes; breaks nothing.
- Corpus files for SO(3): `so3_exp`, `so3_log`, `so3_act`, `so3_from_matrix`, `so3_jr`, `so3_jl`,
  `so3_jr_inv`, `so3_jl_inv` (22 642 records, 8.9 MB), computed from the definitions (quaternion and
  hat-matrix series, Newton's method on the exp series for `Log`, the polar factor for
  `from_matrix`) and each record checked to 100 digits by an independent property (`mp.expm`, the
  sandwich, polar uniqueness, Jacobian identities). Vector ids give each θ an axis; new strata `q:w0` and `q:nonunit`
  for the quaternion ids; every `so3_log` stratum but `q:w0` holds each quaternion and its
  negative. Matrices are column-major with a sibling `shape`. Generation is parallel by stratum
  (`--jobs`), byte-identical to a serial run. `docs/PHASE1.md` §2 and §4.3 state the `Log` and
  `from_matrix` references. No library code changes; breaks nothing.
- Corpus files for SE_N(3), N = 1, 2, 3: `sen3_exp`, `sen3_log`, `sen3_ad`, `sen3_jr`, `sen3_jl`,
  `sen3_jr_inv`, `sen3_jl_inv` (21 files, 7 578 records, 18 MB): `mp.expm` of the hat matrix,
  the inverse of it on θ ∈ [0, π], the conjugation definition of `Ad`, the defining series of
  `ad_τ` as dense rotation-first matrices, and `mp.inverse` of those (the dual matrix's zero
  blocks exactly 0, every other entry the LU's), each record checked to 100 digits (`mp.expm`,
  `J_l = Ad_Exp(τ) J_r`, `J J⁻¹ = I`) with the dual-matrix structure asserted on the data; no
  block form of `NUMERICS.md` §5 is evaluated. New strata: the 25 cells `rho:<scale>/theta=<θ>`
  and the `theta:*` strata at unit translation scale, 6 records each, and `q:w0` and `q:nonunit`
  for `sen3_log` and `sen3_ad`. `docs/PHASE1.md` §4.3 and §4.4 state the `sen3_log` reference
  (`mp.logm` is wrong near π), the inverse rule and the sample count. `just corpus` is 16
  CPU-minutes (3 minutes on eight cores); the `corpus-check` CI job gets a 60-minute timeout. No
  library code changes; breaks nothing.
- Corpus files for SO(2) and SE(2): `so2_exp`, `so2_log` (3 419 records each) and `se2_exp`,
  `se2_log`, `se2_ad`, `se2_jr`, `se2_jl`, `se2_jr_inv`, `se2_jl_inv` (416 each; 2.3 MB): `mp.expm`
  of the hat matrix, its inverse on θ ∈ (−π, π] by Newton's method on the series, the conjugation
  definition of `Ad`, the defining series of `ad_τ` as dense 3×3 rotation-first matrices and
  `mp.inverse` of those, each record checked to 100 digits of the size of its terms (`mp.expm`,
  the complex series, `J_l = Ad_Exp(τ) J_r`, `J J⁻¹ = I`) with the first row (1, 0, 0) and the
  rotation block asserted on the data. No closed form of `NUMERICS.md` §6 is evaluated. The `so2_*`
  ids see the θ strata with each θ in both signs, the `se2_*` ids SE_N(3)'s strata (no
  `theta:dense`) at 8 records each with the sign of θ alternating; z = (−1, ±0), θ = π, and a
  non-unit z are not sampled. `docs/PHASE1.md` §4.3 and §4.4 state the readings. `just corpus` is
  17 CPU-minutes. No library code changes; breaks nothing.
- `docs/decisions/0015-specification-gaps-found-while-building-the-instrument.md` (draft): the 29 gaps
  that the corpus generator, `helicoid-linalg`, `cargo xtask lint` and the maths pages found in
  `PHASE1.md`, `PHASE2.md`, `PHASE3.md` and `NUMERICS.md`, each with the passage, the evidence, the
  options, a recommended resolution and the work it blocks; the decisions in blocking order, then the
  readings and edits of the implementing PRs that await ratification (including the widening of D11 to
  strided reads). It decides nothing and edits no spec; `docs/maths/index.md` gains a pointer to it.
  Documentation only; breaks nothing.
- `helicoid`: the traits `Tangent`, `Jac` and `LieGroup`, the sealed `Side` with its only two
  implementations `Left` and `Right`
  (`docs/PHASE3.md` §2), the conventions of `0002` in the crate docs, and the trivial group
  `Rn<S, N>` with `RnTangent` and `RnJac` (§7): addition is `Mul`, there is no `Add` or `Sub`
  (`compile_fail` doctests), and every row of `NUMERICS.md` §2.3 is written for both sides.
  `Jac::sandwich` and the `DOF` tie are `const` assertions. Generic law checks (test-only
  `laws`) run for `Rn` and for a test-only non-abelian group (Heisenberg) under `f64`, `f32` and
  `Dual<f64, 3>` with recorded bounds. `docs/PHASE3.md` §7 and `docs/API.md` §3 name `RnTangent`
  and `RnJac`. New public API, nothing breaks. `proptest` becomes a dev-dependency of `helicoid`.
- `helicoid`: `Quat<S>` (`docs/PHASE3.md` §4): Hamilton product as `Mul`, `conjugate`, `norm_sq`,
  `to_matrix` (its operand grouping pinned by a bit-exact golden, D16), `from_wxyz_unchecked`
  (debug-asserts `|‖q‖² - 1| <= 2^-40` for `f64`, `2^-16` for `f32`), `from_wxyz_normalized`
  (divides by the norm, `0027`; domain `‖q‖²` normal, `debug_assert!` on the result) and
  `renormalize` (the first-order Newton step, drift repair, asserting nothing), no `Add`, `Sub`
  or `PartialEq` (`compile_fail` doctests), and the converters
  `from_xyzw`, `to_xyzw`, `from_jpl`. New public API, nothing breaks; `SO3` is not here yet.
- `helicoid`: `SEn3Tangent<S, N>` with the `Twist` alias and `Twist::{omega, v,
  from_translation_first, to_translation_first}`, and `SEn3Jac<S, N>` with the dual-matrix algebra
  of `docs/PHASE3.md` §5 as a `Jac` (`mul` at `27 + 54N` multiplications, `inverse`, `neg`,
  `apply`, `apply_transpose`, `write_dense`, `sandwich`, `identity`). `helicoid::reference`
  (public) opens with the dense twins `sen3jac_mul` and `sen3jac_inverse` (Gauss-Jordan), writing
  into caller memory, and the proptests `sen3jac_mul_matches_reference` and
  `sen3jac_inverse_matches_reference` (`N = 1, 2, 3`; `f64`, `f32`, `Dual<f64, 3>`; the inverse's
  tolerance scales with the conditioning of the dense matrix). New public API, nothing breaks; the
  `SEn3` group is not here yet.
- `just test` also runs the workspace tests in the release profile, where `debug_assert!` is
  compiled out: the tests of the documented out-of-domain behaviour are `cfg(not(debug_assertions))`
  and the dev profile skipped them. Nothing breaks.
- `helicoid`: `Product<A, B>` and `ProductJac<JA, JB>` (`docs/PHASE3.md` §7): the product group
  with every `LieGroup` method the factors' method componentwise, tangent `(A::Tangent, B::Tangent)`
  (a `Tangent` impl for a pair, dense order `A` then `B`, `dot_acc` threaded so that `dot` stays
  the one flat index-order sum, `0025`), `DOF = A::DOF + B::DOF`, and the block-diagonal Jacobian as
  a `Jac` (`sandwich` from the factors' `apply`). Nests. Neither type has a public field, as §7 and
  `0025` decision 5 state, so a `Product` is built and read through `exp`/`log` and
  `write_dense`/`read_dense`. `helicoid::reference::productjac_sandwich` is its dense twin, with
  `productjac_sandwich_matches_reference`. Tested by the generic laws on `Product<Rn<3>, Rn<2>>`,
  `Product<Rn<2>, Heis>`, `Product<Heis, Rn<2>>` and a nested product with two Heisenberg factors
  (test-only, the one non-abelian factor until SO(3)). New public API, nothing breaks.
- `cargo xtask conformance [--subject NAME] [--fn ID]` (`just conformance`), the harness core of
  `docs/PHASE1.md` §5: the corpus reader, the `Subject` trait with an empty in-process registry, the
  forward error of `NUMERICS.md` §11 in units of `u = 2⁻⁵³` formed exactly in integers and rounded
  once (the 30-digit reference is never parsed to `f64`), per `(fn, stratum, precision, subject)`
  `n`, `max_u`, `p99_u` (nearest rank), `argmax_id` and `nonfinite`, written to
  `conformance/results/<subject>.csv` and printed by `max_u` descending. A non-finite output, an
  error that overflows binary64 and a run that scored nothing (the registry is empty, so today's
  run) fail; `so3_from_matrix` is checked for shape and finiteness only. The floor of a tangent or a
  coefficient (`2⁻¹⁰²²`), the SE_N(3) tangent's sign, the maximum over a record's fields, one norm
  over a whole tangent and overflow as non-finite are readings recorded in `xtask/src/conformance/`
  and proposed by draft record 0014. Backward error and the `helicoid` subject are not implemented
  (`f32`, `--self-test` and the seeded subjects follow in the entries below). `xtask` gains
  `num-bigint` and `helicoid-linalg`; no library code changes; breaks nothing.
- `cargo xtask conformance --self-test` (`just conformance --self-test`, a CI job) and the seeded
  subjects of `docs/PHASE1.md` §10 (`xtask/src/seeded`): the coefficient kernels `k, a, b, c, d, e`
  generic over `Real` (exact arm at the safe argument, one `S::branch`; series read at run time from
  `coeff_series.jsonl` and rounded once, in integers, at the precision they run at, `f64` or `f32`),
  evaluated over the `coeff_*` ids as value and `d/dz` of one `Dual<f64, 1>` at a
  `Candidate { terms, switch_z }`; `seeded:correct` is the `tf_tree` D12 candidate (for all six
  coefficients, though `NUMERICS.md` §4 lists it for `a`, `b`, `c`), and `b` by its definition and `k`
  with an unsafe `sqrt` of `θ²` under `Dual` are planted. The self-test fails unless the correct kernel
  fires no mechanism and each defect fires its own (`b`: the value's error curve over
  `theta:1e-8`…`theta:1e-2` fits `θ^-p`, `p` in [1.8, 2.2], 1.937 measured; `k`: `nonfinite > 0` in
  `theta:exact0`), and prints the correct kernel's errors and every mechanism's reading for every
  subject. Unit tests pin the correct kernel's `max_u` per coefficient and field, and check its safe
  argument on a mask that evaluates both arms. A plain `just conformance` now scores `seeded:correct`
  and skips planted subjects unless named. `xtask` gains `libm` (the fit's logarithm, the same bits on
  every host). No library code changes; breaks nothing.
- `cargo xtask thresholds [--check]` (`just thresholds`, `thresholds-check`, a CI job), the sweep of
  `docs/PHASE1.md` §6 over the seeded kernels at `f64` for `k, a, b, c, d, e`: 1 to 8 series terms
  times 1025 switch points (64 per decade of `z = θ²`, `θ` from `1e-8` to `1`, each the correctly
  rounded `10^(i/64)`), each arm's error formed once per record over every `theta:*` record (value
  and `d/dz` through `Dual<f64, 1>`, exact, in `u`), the maximum minimised, ties to fewer terms
  then the larger switch, the `tf_tree` D12 prior scored beside the choice. Writes
  `conformance/sweeps/thresholds.csv` (committed, byte-identical, compared with a fresh run by a
  test). The rows are of this corpus and grid, not general optima: `terms` is the cap of 8 for every row but `e`;
  `k` and `d` sit at the top of the grid (`θ = 1`), `b` one step below it (unique, 1.3% under the
  top), `a` and `c` inside it; `e`'s objective is one record at `θ = 1` exactly, which `z < switch`
  never puts on the series arm (a top of `nextUp(1)` gives it 8 terms and 3.3 times less). Columns
  name the record that attains each objective and the objective one grid step either side. The
  grid, top, prior and tie readings are `0014` (draft) questions 8 to 11. Not swept: `r`, `f32`,
  SE(2)'s `α`, `β`; `generated.rs` waits for Phase 3. No new dependency, no library code changes;
  breaks nothing.
- `cargo xtask thresholds` also writes `xtask/src/seeded/generated.rs` (`docs/PHASE1.md` §6, `0004`
  item 3), registered with `cargo xtask lint`: per coefficient a `Switch<f64, terms>` (the switch as a
  bit pattern with its decimal in a comment, the series as the shortest decimal that reads back as
  each exact rational rounded once at binary64, laid out as `rustfmt` leaves it, an objective line)
  under the `@generated` header and the SHA-256 of the sweep CSV and of `coeff_series.jsonl` as the
  source revisions. `just thresholds-check` (the CI job) now fails on any difference in the CSV or in
  that file and writes nothing; tests edit one character of each and drive the entry point under a
  scratch root, and the check names the file and the line (an edit that no longer compiles stops
  `xtask` building: `git restore` the file). `seeded:correct` runs the generated switches, where it
  ran the `tf_tree` D12 candidate (still the sweep's prior, its errors still pinned): its per-field
  maxima equal the CSV's, to the bit. `seeded:c-two-terms-1e-8`, `c` with switch `1e-8` and two
  terms, joins `--self-test`: detected when the sweep's objective for its candidate is more than
  `10^6` times the chosen one's (measured `10^30`, rank 7841 of 8200), while the correct kernel,
  whose `c` is the chosen candidate, stays silent; the self-test also fails when a subject's
  measured `c` objective is not the ranked one, to the bit. `xtask` gains `sha2` (the revisions;
  `xtask` only, `0007`); no library code changes; breaks nothing.
- Seeded SO(3) subject and its two `Log` defects (`docs/PHASE1.md` §10, `xtask/src/seeded/so3.rs`):
  `seeded:correct` now also runs `Exp` (`k` at the generated switch, `cos θ/2` by its exact arm: no
  series is committed for it) and `Log` (the `copysign` flip, `r = 2 atan2(n, w)/n` at the safe
  argument, its series arm `2/w` taken where `n² = 0`, no constant typed) over `so3_exp` and
  `so3_log`, at most 2.91 `u` on every stratum (`--self-test` gates 4), non-finite-free, its safe
  arguments checked on a mask that evaluates both arms. `seeded:log-acos` and `seeded:log-no-flip`
  join `--self-test`, which fails unless each fires its own mechanism and the correct kernel fires
  none: `acos` when every stratum of `theta:1e-k` from `k = 4` and of `theta:pi-1e-k` from `k = 7`
  reaches `10^7` (§10's per-`k` wording cannot hold for every `k`); the flip when every stratum with
  a negated half either fails or has nothing to detect (`n² = 0`, derived from the records: 27 of 29
  fail, `theta:exact0` and `theta:subnormal` cannot), and one fails. A stratum with no score, a
  non-finite output and a `NaN` fail the gate. The two defects run over `so3_log` only. The readings
  are 0014 (draft) questions 16 to 20. No new dependency, no library code changes; breaks nothing.
- Seeded SE_N(3) subject and its two defects (`docs/PHASE1.md` §10, `xtask/src/seeded/se3.rs`):
  `seeded:correct` now also runs `Exp` (`so3::exp` and `J_l(φ)ρ_i`) and the Jacobians `J_l` and
  `J_r = J_l(−τ)` with Barfoot's block `Q` of `NUMERICS.md` §5.3 over `sen3_{exp,jr,jl}_n{1,2,3}`
  (`k, a, b, d, e` each at its own generated switch, dense column-major output): at most 3.47 `u`
  (`Exp`) and 8.90 `u` (`J`) on every stratum (`--self-test` gates 12), no non-finite output; the
  errors are of this per-coefficient kernel, not of the grouped one of `NUMERICS.md` §4.
  `seeded:se3-exp-translation-first` and `seeded:q-minus-half` join `--self-test`, which fails
  unless each fires its own mechanism and the correct kernel fires none: every `rho:*` stratum of
  `sen3_exp_n{1,2,3}` reaches `10^7` `u`, and so does every stratum of every `sen3_jr_n*` and
  `sen3_jl_n*` (measured at least 8.9e15 and 5.2e9). §10's `Dual` comparison of the `Q` defect
  waits for the `helicoid` subject. The readings (the bar, "fails", the strata counted, the
  per-coefficient switches, the association of `Q`'s words) are 0014 (draft) questions 21 to 25. No
  new dependency, no library code changes; breaks nothing.
- `runners/tf_tree_math`, oracle #1 (`docs/PHASE1.md` §7, `0010`): a workspace-excluded runner with
  its own `Cargo.lock` and `tf_tree` commit `20bc5a05` as a pinned git dependency. It reads corpus
  JSONL and writes hex-float answers for `so3_exp`, `so3_log`, `sen3_exp_n1` and `sen3_log_n1`,
  through `helicoid_to_tf_tree_*` / `tf_tree_to_helicoid_*` conversions with hand-computed tests on
  values whose components all differ. `cargo xtask conformance --oracle NAME` (`FileSubject`) builds
  and runs a runner and scores its answers with the harness's exact metric, `subject_version` being
  the pin; a skipped or extra record, a stray file and an owed id with no file are errors, a
  non-finite answer a recorded row. `just oracle-tf-tree-math` runs it, with the runner's `doc` and
  `cargo deny` (`deny.toml` allows its git source; Dependabot covers its registry dependencies).
  Measured: at most 4.69 `u` on every stratum but `so3_log` `theta:subnormal` (3.97e14, an
  underflow to 0). The pin exports no `V`/`V⁻¹` and no Jacobian, so there is no `coeff_*` row; the
  file protocol and the ids not answered are 0014 (draft) questions 26 and 27. No library code
  changes; breaks nothing.
- `runners/sophus_rs`, oracle #2 (`docs/PHASE1.md` §7): a workspace-excluded runner with its own
  `Cargo.lock` and `sophus_lie` 0.15.0 (with `sophus_autodiff` 0.15.0) as exact registry pins. It
  answers `so3_{exp,log,jl,jr,jl_inv,jr_inv}` and `sen3_{exp,log,jl,jr,jl_inv,jr_inv}_n1` through
  `helicoid_to_sophus_rs_*` / `sophus_rs_to_helicoid_*` conversions with hand-computed tests: 0.15.0
  is rotation-first and `w`-first like `0002` (0.14.0 and earlier are translation-first), so they
  change a layout, not an order. The right Jacobians are the left ones at `−τ`. `just
  oracle-sophus-rs` runs it, audited with the runner's own `deny.toml` (one `paste` advisory).
  Measured, 484 rows, no non-finite output: `so3_exp` and `so3_log` at most 2.74 `u`, and every
  Jacobian up to 7.3e12 `u` for `θ ≤ 1e-3`, where `left_jacobian` and `inv_left_jacobian` take a
  small-angle arm of the wrong sign; `sen3_exp_n1` at most 5.8e8 `u`; `sen3_log_n1` 1.78e16 `u` in
  `q:w0`, the other branch of `Log` at `θ = π`. The per-stratum rows are only in the untracked
  `conformance/results/sophus_rs.csv`. The pin, the right Jacobians and that branch are 0014 (draft)
  questions 28 to 30. No library code changes; breaks nothing.
- `cargo xtask envelope` (`just envelope`, a CI job; `docs/PHASE1.md` §8, `0006`): merges
  `conformance/results/<subject>.csv` of the candidate (`helicoid`) and the oracle runners and judges
  per `(fn, stratum, precision)` on the max: domination over the smallest oracle max (ties pass; a
  non-finite or unscored oracle row is no bar), no-regress against `conformance/baseline/`, compared
  exactly, and any non-finite or unexpectedly unscored output; then coverage of the corpus (the 45 ids
  of §4.3 required, each with a file that exists; the 16 that later phases name owed until the §0.0
  row that owns them is `Done`). `--bless` writes the baseline and the `@generated`
  `docs/evidence/ENVELOPE.md` (registered with `cargo xtask lint`) and writes nothing while a
  failure stands that a new baseline would not move (domination, coverage, a dropped or narrowed
  stratum); `--check` fails on a hand edit or an improvement not blessed, the page compared without
  its oracle columns (the runners' `libm` is not D16's); `--dry-run` writes nothing. `helicoid` is
  not an in-process subject before Phase 3, so today the run checks coverage and the page only, and
  no baseline is committed (blessing is `PHASE3.md` §10's, and 0014 (draft) says none before it is
  `ready`); once it is registered, no rows fails the run: the CI job becomes the gate with the
  subject. `--self-test` now also runs §10's envelope half: the planted `c` loses to the correct
  seeded kernel on 9 strata of `coeff_c`. Measured on the stand-in `seeded:correct`, not a bar: 214
  strata paired with an oracle, an oracle beats it on 75. The readings are 0014 (draft) questions 31
  to 35. No library code changes; breaks nothing.
- `helicoid-linalg`: `Real::cbrt` (`docs/PHASE2.md` §2 and §3, `docs/API.md` §2, `docs/NUMERICS.md`
  §12; step 1 of `0017`, which `solve_cubic` needs): the real cube root, total and odd, through
  `libm::cbrt`/`cbrtf`. `Dual` differentiates it as `d / (3 c^2)`, `c = cbrt v`: infinite at 0
  like `sqrt`, with the sign of `d` at either zero, accurate elsewhere, `0` at `+-inf`.
  Checked against mpmath on chosen rows (both signs, subnormals; first derivative at `f64` and
  `f32`, second derivative and the Hessian of `cbrt(x y)` at `f64`) to 2 ulp; random inputs reach
  about 4 ulp on the first derivative and 10 on the second, the value 0 ulp. **Breaks** an
  out-of-tree `impl Real`, which gains a required method (none is known; the in-tree ones are
  updated); nothing else breaks.
- `helicoid-linalg`: `solve_cubic` (`docs/PHASE2.md` §6, step 2 of `0017`): the real roots of `a x³
  + b x² + c x + d` as `(Vec3<S>, [S::Mask; 3])`, slot `k` reported as a root iff its mask is set, a
  clear slot `+0`. omnisac's algorithm, generic over `S: Real` (Cardano, Viète's trigonometric form,
  a discriminant band), every branch a `branch`/`select` at a safe argument, so it runs on `f64`,
  `f32`, lanes and `Dual`. Differences from omnisac's `poly::solve_cubic`: the container (no
  `ArrayVec`, no compaction); the tolerances, `1e-14` to `2^-46` (128 `u`, from 90 `u`) and `1e-12`
  to `2^-40` (8192 `u`, from 9007 `u`) at `f64`, `2^-17` and `2^-11` at `f32`, the nearest
  power-of-two multiples of `u` (a polynomial whose leading coefficient or discriminant lies in the
  thin strip between old and new changes its root count); `acos` as `atan2(sqrt((1 - x)(1 + x)),
  x)`, within 2 ulp of `libm::acos`, and `pi` as `atan2(+0, -1)` — a reading `0022` has since
  reversed, and `pi` is already a literal: the roots of the trigonometric arm move by up to about
  11 `u` of the largest root, an observed maximum and not a bound, which `0022` step 1 removes. Everything else agrees to the bit with omnisac (on `libm`) on
  10^6 random and adversarial polynomials per seed (a scratch differential, not committed; `GOLDEN`
  in the tests is 36 of its rows). Tested against omnisac's cases, an mpmath fixture, planted roots
  against a measured error model, the slot order, the `Dual` derivative, two lanes and the mask.
  Limits it inherits, in the rustdoc `# Domain` and `0031` (draft), fixed by none of this: the
  one-real-root arm is not backward stable (`x^3 + p x - 1` is off by `4e-6` at `p` near `1e-5`,
  `f64`, and its `Dual` derivative is `-inf` below `p = 1.3e-5`); a repeated root can be dropped
  where `p` and `q` cancel; a cubic with one real root can read as three where `p^3` and `q^2`
  underflow. `docs/API.md` R6 names `solve_cubic` beside `chol` as the functions that assert
  nothing. New item: breaks nothing; omnisac's call sites need an adapter (`docs/PHASE2.md` §9).
- `helicoid-linalg`: `eig3` (`docs/PHASE2.md` §6, step 3 of `0017`): the eigenvalues, ascending, and
  an orthonormal basis of eigenvectors with `det = +1` of a symmetric `Mat3<S>`, as `(Vec3<S>,
  Mat3<S>)`; the lower triangle is read. omnisac's `eigendecomp_sym3_signed` (d3be7b7), generic over
  `S: Real` (Smith's eigenvalues, cross-product eigenvectors; `f64`, `f32`, lanes, `Dual`). What
  differs from omnisac is proposed in the draft record `0023`: no `Option` and no tolerance, so the
  eigenvalues are not finite for a non-finite entry and where `p^3` over- or underflows, and a frame
  built from cross products of unit vectors, orthonormal to a few `u` for every input. **Not done:
  Kopp's hybrid is owed** (§6). Where the two largest eigenvalues are within about `2 sqrt(u)|A|`
  (equal ones included) no column is reliable, the isolated smallest eigenvalue's included (residual
  of the order of `|A|`); omnisac's vectors are as wrong there, or `None` where a cross product
  vanishes, and anchoring the frame on the isolated end may fix it (`0023` (draft)). At a double
  root the eigenvalues lose half their digits (`0.6 sqrt(u)|A|` on a rank-1 matrix), which only the
  hybrid fixes. §6 names no constant for its bars, so none is claimed; the measurements (a scratch
  harness, not committed) are in `eig3_tests` and the §0.0 row. `docs/API.md` R6 names `eig3` beside
  `chol` and `solve_cubic`. New item: breaks nothing; omnisac's call sites need an adapter
  (`docs/PHASE2.md` §9).

### Changed

- `SEn3Jac::sandwich` builds both halves with `from_cols`/`from_rows` instead of filling and then
  overwriting two `D x D` scratch matrices: `M = J Σ` column by column, then `M Jᵀ` row by row, since
  row `r` of `M Jᵀ` is `J` applied to row `r` of `M`. It therefore forms neither `Aᵀ` nor the `N`
  transposes of `col`, reaches no entry through `Matrix::get`/`set`, and emits no `memset` where the
  block form emitted two of `D²` scalars that it then fully overwrote. The private `zeros`, `block`
  and `put` of `dualmat.rs` are deleted with it, which withdraws the `dualmat.rs` half of
  `PROJECT.md` §5.1's gate on `Matrix::block`/`set_block`.
  **Bit-identical**: every entry is the same products summed in the same order, because `Matrix`'s
  `Mul` and `Mul<Vector>` associate identically and `B_i w_0` is added to `A w_{i+1}` in `apply`'s
  order — checked over 150 000 random `(J, Σ)` at `f64`, `f32` and `Dual<f64, 3>` for `N = 1, 2, 3`
  with not one differing entry, so D16 holds and the recorded bounds are unmoved rather than
  re-measured. Release, x86_64, minimum of 7 runs of 200 000 calls: 545 → 48 ns (`f64`, `D = 6`),
  2406 → 448 ns (`f64`, `D = 12`), 2301 → 606 ns (`f32`, `D = 12`), 1747 → 680 ns and 8078 → 2328 ns
  (`Dual<f64, 3>`, `D = 6` and `12`). The result is fully packed: 162 `mulpd` and no `mulsd` for the
  324 scalar products at `N = 1`. A middle variant that built each block's `Vector` and gathered it
  was faster still on `f64`/`f32` (365 and 327 ns at `D = 12`) but 17% to 20% slower on `Dual`, and
  `0006`'s per-precision no-regress bar refuses it. Values unchanged, so nothing breaks.
- `SEn3Tangent::read_dense` and `write_dense` take one length test instead of one bounds check per
  component, on the path `ProductJac::sandwich` crosses once per column and once per row of its
  argument and again per nesting level. At `N = 3` in release on x86_64 the in-domain arm of
  `read_dense` is one `cmp`/`jb` and twelve loads, and `write_dense`'s is one `cmp`/`jb` and a
  `movups`-packed copy with no `memcpy` call — the copy vectorizes once the `Chain`/`FlatMap` of
  `comps` is gone, which
  `dot_acc` had already avoided for its own reason. Out-of-domain behaviour is unchanged and still
  pinned by `out_of_domain_does_not_panic_in_release`: a short `write_dense` view takes the prefix
  that fits, and `read_dense` poisons only the missing entries, never the present ones
  (`0025` decision 4).
- `dualmat_tests`' `kappa` takes `M⁻¹` from the closed form `A⁻¹ − ε A⁻¹ B_i A⁻¹` rather than from the
  `sen3jac_inverse` twin. The twin `debug_assert!`s a nonzero pivot, so the `κ u <= 1e-3` filter was
  computed from a call that panics on exactly the singular draws the filter exists to reject —
  `inverse_twin::<f64, 1, 6>(&[0.0; 18])` panicked instead of returning `None`, and because proptest
  shrinks toward `0`, that panic and not the bound violation is what a real regression would have
  reported. `Mat3::inverse_adj` asserts nothing, so a singular `A` now yields a non-finite `κ` that
  fails the filter like any other out-of-domain draw, and the case costs one 3x3 inverse instead of a
  second `D³` elimination. All twelve recorded bounds reproduce from the documented `10^6`-case run.
- `laws`' `Rng::vec` becomes `Rng::arr::<N>`, filling a caller-owned array with an explicit loop: the
  `10^6`-case measurements no longer allocate once per draw-set per iteration, and the draw order is
  a written loop rather than an unspecified visiting order, since that order is the stream every
  recorded figure is reproducible from.
- `laws`' dense reader is `reference::dense_oriented`, which `reference::dense` also calls: the
  NaN-poisoned scratch had two definitions differing only in orientation, and the poison semantics
  now have one place to be kept. `reference::sum`'s doc states what couples it to
  `helicoid_linalg`'s `vector::sum` and which tests assert the two associate identically.
- `zero3`, `put_view` and the `dualmat.rs` helpers that survive are `#[inline]`; the comment on
  `sandwich`'s scratch no longer claims a saving it did not make.

- `libm` is built with its `arch` feature (`0018`): `sqrt`, `sqrtf` (x86_64, aarch64), `fma`, `fmaf`
  and `rint`, `rintf` (aarch64; `fma` by `cpuid` on x86_64) run the target's instruction. All are exactly
  rounded, so every output that is not NaN is unchanged on every target; the sign and payload of a NaN
  from `sqrt` and `fma` are now the target's, as those of `+ - * /` always were (so `copysign` of the
  NaN of `sqrt(-1)` takes the target's sign; that input is outside `Real::sqrt`'s domain).
  `Real::sqrt` is about 4x faster on x86_64 (measured; `docs/decisions/0018`, Consequences). wasm32 and
  thumbv7em are unchanged: `libm` routes nothing there on stable. New test: `Real::sqrt` for `f64` and
  `f32` against an integer square root, bit for bit, over `2 x 10^6` cases with a pinned digest; checked
  by hand on x86_64, aarch64 (`qemu-user`) and wasm32 (node WASI), run by CI on x86_64 and aarch64.
- `helicoid-linalg`: `chol_solve(&l, b)`, `A x = b` from a Cholesky factor: `solve_lower`, then a
  back substitution against `l^T` read by column, with no transposed copy. Bit-identical to
  `solve_upper(&l.transpose(), solve_lower(&l, b))` (NaN sign and payload aside), same domain and
  release behaviour; `docs/decisions/0019`, `docs/NUMERICS.md` §15.6. New public API, nothing breaks.
- `Tangent::dot` becomes a provided method over the new required `dot_acc`, and
  `Tangent::read_dense` reads a missing entry as NaN rather than `+0` (`0025`): `+0` is a valid
  component, so a wrongly sized buffer was indistinguishable from a zero tangent.
- The law harness divides every norm-wise error by a NaN-preserving denominator. `f64::max` returns
  the other operand for a NaN, so a NaN reference had been scaling away to a denominator of `1`,
  and a law could pass against it.
- The non-abelian test group's sample Jacobian is no longer a circulant. Circulants commute, so
  `jac_dense_order` — which pins `Jac::mul`'s operand order against the dense product — could not
  detect a `mul` that multiplied its operands the wrong way round. Row-scaling breaks the
  commutation; the recorded `f64` bounds move to `jac_order` 12 and `sandwich` 12, `f32` to 11 and 9.
- `Tangent::read_dense`/`write_dense` take the exact-length path through `first_chunk`, so the
  documented boundary costs one length test instead of one per component (measured at `DOF = 9` on
  the emitted release asm: 1 branch and vectorized moves, against 17 and 9 scalar ones).
