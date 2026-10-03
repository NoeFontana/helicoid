# 0034: A second transcendental backend owes a measurement

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** proposes a D16 reword and two guards; no backend change, no public surface.

## Context

D16 reads "every transcendental through the `libm` crate; no `mul_add`, no `target-cpu`, no
fast-math. Outputs are bit-identical on x86_64, aarch64 and wasm32." That states an
**implementation** and claims a **property**, and the two are not the same thing:

- The property holds **per resolved `libm` version**. `libm` 0.2.16's architecture-specific paths
  cover only IEEE-exact operations — `sqrt`, `fma`, `rint`, `floor`/`ceil`/`trunc`, `fabs` — and
  `sin`, `atan2`, `exp`, `expm1`, `log1p` contain no `fma`, so the cross-target claim is sound for a
  given version. The workspace declares `libm = { version = "0.2", default-features = false,
  features = ["arch"] }`, which is right for a library — an `=` pin would force every consumer's
  graph onto one patch release, and the lockfile belongs to applications — and it means **a patch
  release may legally move bits**.
- Nothing records which `libm` produced a committed baseline. When one moves, the envelope's
  no-regress bar and `baseline/HOST.md` report unexplained drift rather than "`libm` moved".

The mechanism D16 names is already a single chokepoint, which is what makes the reword cheap:
`helicoid-linalg`'s `float.rs` is one macro `impl Real for f32/f64` holding every `libm::` call —
`sqrt`, `cbrt`, `sincos`, `atan2`, `fabs`, `copysign`, six in total, `sincos` so there is one range
reduction — and `dual.rs` has none, so `Dual` inherits through `Real` ([`0003`](./0003-the-scalar-that-cannot-say-less-than.md)).
`no_std` removes `std`'s inherent `f64::sin`, so the compiler enforces half of the containment.
**Nothing enforces the other half**, and nothing records the version.

## Measurement

Measured here, on this repository:

- `libm::atan2` against the host's glibc, over the corpus's 30-digit references:
  `libm`'s $\nu$ is 1.96 against glibc's $\le 1.00$, and the whole of `so3_log`'s 8 domination
  failures is that difference, at most **0.74 u**
  ([`0032`](./0032-domination-charges-helicoid-for-d16.md), draft). So D16's accuracy cost is real,
  bounded, and already attributed. Its **latency** cost is unmeasured anywhere in this stack.
- The seam's shape, above: 6 `libm::` call sites in one private module, 0 in `Dual`.

Taken from the consumers' own evidence, not re-measured here:

- tf_tree's lookup profile is led by the bracket path — `sample.rs` at 27.1% of instructions and
  76.8% of mispredicts, then `quat.rs` and `buffer.rs`. **No `libm` row appears.** Building with
  `x86-64-v3` made `at_many` *slower*. tf_tree's `0016` rejected `pulp` at about +12% for a
  dependency count going 2 → 11.
- tf_tree already pays the `libm` cost today (`libm::sincos`, plus `acos` and `sin` in `LerpSlerp`).
- locus-tag is the one consumer that *changes* backend: it goes through nalgebra with `std`, so
  glibc on Linux, and `libm` appears only in its test utilities. Its migration is the first place
  D16 could cost measurable time.

The arithmetic that makes the ladder's order inevitable: end-to-end speedup is
$1/(1 - f + f/k)$ for a transcendental share $f$ and a backend speedup $k$. At $f = 0.15$, $k = 2$
that is **1.08×** — below the bar tf_tree's `0016` used to reject a 12% gain. A second backend has
to clear that bar *end to end*, not per call.

## Decision

Proposed, not settled. **Build the seam's guard and the measurement. Do not build the knob.**

1. **Reword D16 from implementation to property**, so that a later kernel swap is not a breaking
   change:
   - bit-identical on x86_64, aarch64 and wasm32 for a given `helicoid` version **and a given
     resolved `libm` version**;
   - per-function error bounded by the envelope;
   - implemented through one private kernel, currently the `libm` crate; no `mul_add`, no
     `target-cpu`, no fast-math.

   The three prohibitions stay literal. This is a `PROJECT.md` §5 edit and is **not** made by this
   record.
2. **Guard what already holds.** A `cargo xtask lint` check that `libm::` appears only in the
   kernel module, so the containment cannot erode silently — the compiler covers only the half
   `no_std` removes. And every generated evidence header records the **resolved `libm` version**,
   so a bump is named rather than observed as drift.
3. **The trigger for a second backend**, as `0013`'s sibling: a committed consumer benchmark where
   measured $f$ and $k$ clear tf_tree `0016`'s bar *end to end*, **and** a portable kernel cannot
   close the gap. Until both hold, this record closes with no change.
4. **A bit change is a version bump.** A kernel change that moves bits is a **minor-version
   release** of `helicoid` — never a patch — and a consumer that needs something stronger pins an
   exact version, which is the consumer's call and not a surface `helicoid` offers. Decided by the
   owner. This is what makes ladder step (iii) available at all: without it, a faster portable
   kernel has no release shape and the ladder ends at (ii). `API.md` §5 owes the sentence.
5. **The ladder, in order**, when the trigger fires:
   1. close with no change;
   2. algorithmic savings — two-stage geodesics, reusing `sincos`, half-angle forms. These cost
      nothing in determinism;
   3. faster **portable** kernels: pure Rust, basic operations only, a fixed operation sequence,
      restricted to the Lie argument domains, verified by the corpus that already exists, and
      upstreamed to `libm` where possible so `0007` stays at "`libm` only". Bits change once, as a
      minor-version event, and no knob is needed;
   4. only if (iii) cannot close the gap, a build-time `--cfg helicoid_math` owned by the binary,
      `libm` the default, and a `CONTRACT` const that consumers `const`-assert, so a conflict is a
      compile error. **Never a Cargo feature**, because features unify across a graph and a
      bit-identity requirement cannot survive that.

## Rationale

Ownership decides the shape. The application knows its platform and whether it replays; tf_tree and
fuse-geometry know their bit-identity promises. `helicoid` knows neither, so it owns the mechanism
and the strongest default and never the choice.

Reversibility decides the order. A private chokepoint is a two-way door — swap the implementation,
regenerate the baseline. A public knob is one-way: once a consumer builds on a `system` backend, it
is supported forever, and D16 is load-bearing in four places that each become a test matrix — the
envelope's no-regress bar, fuse-geometry's native/wasm rebuild proof, tf_tree's `at_many ≡ at` and
cross-process arena, and omnisac's `RunRecord` gates.

The measurement decides whether to start at all, and the only profile in the stack does not show
`libm` at the top. That is why (iii) precedes (iv) and why both sit behind a trigger.

A rejected alternative: `exp_with::<M>` per operation, a backend type parameter on group methods.
`helicoid` is generic over `S: Real`, so the seam already lives at the scalar; a second axis on
every group method fails `API.md` §6's checklist and has no consumer.

## Consequences

- D16's reword makes a kernel swap a minor-version event rather than a breach, which is the policy
  the owner has set. `API.md` §5 then owes the sentence: bits may change in a minor release, never in
  a patch, and the envelope never regresses. A consumer whose store or replay cannot absorb that
  pins an exact `helicoid` version; fuse-geometry is the likeliest to, and that is fuse's decision.
- **A requirement on omnisac's migration PR, not on this record.** The Phase 2 "pure move,
  `RunRecord`s bit-identical" gate holds only if omnisac's `eig3`/`svd3`/`solve_cubic` already call
  `libm`. omnisac is private and not readable from here, so that PR must check it and, if they call
  `std`, split into two — switch omnisac to `libm` in place (`RunRecord`s change once, explained),
  then move (bit-identical). One gate takes one variable.
- The containment lint makes `float.rs` the only place a transcendental may enter. A new
  transcendental is a `Real` method plus one line there, which is the shape `0022` already assumes.
- Recording the resolved `libm` version means a `libm` bump invalidates baselines deliberately, in
  its own PR, with the envelope saying so.
- Cross-target digests **stay in Phase 6**, by the owner's call. `PHASE6.md` already owns them, and
  the kernel question is not reopened before then, so the ordering costs nothing. Noted because
  `CLAUDE.md`'s recipe table claimed `just determinism`, `just oracles` and `just bench-check`, none
  of which exist in the `justfile`; `PHASE1.md` §0.0 is the accurate one ("`determinism`, `oracles`,
  `bench-check` jobs not started"), and the table is corrected alongside this record.

## Implementation plan

1. The `libm::`-containment lint, with a planted-violation test — verified by `cargo xtask lint`.
   **Owed.**
2. The resolved `libm` version in `baseline/HOST.md` and `docs/evidence/ENVELOPE.md` headers, and
   an envelope message that names a `libm` move — verified by the generators' own `--check`.
   **Owed.**
3. Transcendentals as conformance subjects in their own right: max and p99 ULP per stratum against
   mpmath for `libm` and for bench-only candidates (`std`/glibc, `core-math`). Runner dependencies
   are unrestricted under `0007`, so the budget does not move. `0032` has the `atan2` row already.
   **Owed.**
4. `PROJECT.md` §5's D16 reword and `API.md` §5's stability sentence ("bits may change in a minor
   release, never in a patch; the envelope never regresses") — only when this record is `ready`,
   which now waits on the reword's wording alone. **Owed.**
5. A two-stage geodesic form in Phase 4 — `Geodesic::new(a, b)` storing $\Omega$, `.at(τ)` costing
   one `exp` — if `PHASE4.md` only specifies `geodesic(a, b, τ)`. Caching the log is the one
   transcendental saving that never touches determinism; adopting it at bracket level is tf_tree's
   call. **Owed, and a `PHASE4.md` question before it is a code question.**

## Open questions

1. **Does the owner accept the three-clause D16 reword as worded in the Decision?** That is the only
   question left between this record and `ready`. D16 is a hard rule, so the wording is the owner's
   and is not changed by this record.

Answered by the owner, and kept here because the reasoning is load-bearing:

- *Is a bit change a minor-version event or a pin?* **A version bump**, with a pin available to a
  consumer that needs one. Folded into the Decision as item 4; it is what makes step 5(iii)
  available.
- *Does omnisac call `libm` or `std`?* **Unknown and not readable from here** — omnisac is private.
  It therefore constrains omnisac's own migration PR, recorded under Consequences, and does not block
  this record.
- *Should cross-target digests come forward from Phase 6?* **No, later.** They stay in `PHASE6.md`.
- locus_fusion is C++, so bits never cross that boundary; its contract is the envelope via the
  corpus in its own suite. Only a correctly-rounded backend would extend bit identity across
  languages. Recorded so it is not rediscovered as a gap.
