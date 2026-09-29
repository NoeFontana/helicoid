# 0004: Switch points are generated, not typed

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

`tf_tree` D12 measured one switch (θ < 0.1, four terms) for three coefficients ($a$, $b$, $c$) at
`f64`, against the `1e-8` most libraries use. `helicoid` has seven coefficients with different
cancellation classes (`NUMERICS.md` §4): $e$ loses $\sim 360u/\theta^4$ and keeps about one digit at
$\theta = 10^{-3}$ from its definition, while $k$ has no cancellation at all. `f32` moves every
switch. A single typed threshold is wrong by construction for some coefficient at some precision,
and a typed series coefficient is a transcription risk with no test that would notice a small
error.

## Decision

1. For every coefficient and precision, `cargo xtask thresholds` chooses the number of series terms
   and the switch point that minimize the max over the corpus strata of max(value error,
   derivative error through `Dual`) (`PHASE1.md` §6).
2. Series coefficients are generated as exact rationals from mpmath's Taylor expansion; the leading
   four in `NUMERICS.md` §4 are asserted by the generator.
3. `coeffs/generated.rs` is `@generated`, with switch points as bit patterns and literals
   correctly rounded per precision. `just thresholds-check` fails on drift or hand edits.
4. The sweep evaluates `helicoid::coeffs` itself (private `__sweep` feature): the measured code is
   the shipped code.
5. Every generated switch point has a `branch_continuity_*` test.
6. `tf_tree` D12's choice is kept in the sweep as a named prior and reported beside the result.

## Rationale

`tf_tree` D12 already forbids asserting thresholds; this generalizes its method from one measured
constant to a generator, which is the only form that survives seven coefficients, two precisions
and future kernels. Alternatives lost: a hand-maintained table (the thing D12 forbids); a single
conservative threshold (costs digits on $d$, $e$ in exactly the regime VIO lives in).

## Consequences

- A change to a switch point changes outputs without changing the API; it ships in a minor release
  with a changelog line (API §5).
- The sweep CSV is committed evidence; reviewing a threshold change means reviewing that CSV diff.

## Implementation plan

1. Sweep tool against seeded kernels (`PHASE1.md` §6) — verified by byte-identical regeneration.
2. `helicoid::coeffs` + `__sweep`; first `generated.rs` (`PHASE3.md` §3) — verified by
   `thresholds-check` and `branch_continuity_*`.

## Open questions

None.
