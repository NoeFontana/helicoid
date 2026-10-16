---
paths: ["crates/**/src/**"]
---

# Writing numeric code

- Read the formula's page in `docs/maths/` first: it is checked to 110 digits and already holds the
  cancellation-free forms and the refuted alternatives.
- `S::branch` arms that only *choose* among precomputed values return scalars (a tuple, not an
  array) and capture by `move`; a borrowing closure makes LLVM spill every candidate and load
  through a `cmov`-selected stack address (measured 0.88x vs 0.54x on `SO3::from_matrix`, same bits).
- Pass scalars, not `[S; N]`, across an `#[inline(never)]` boundary; return small matrices only
  when the consumer applies them. Keep rarely-taken ladders `#[cold] #[inline(never)]`.
- Index-form sums (`(0..N).map(|i| a[i]*b[i])`) not `into_iter().zip`: `array::IntoIter` spills
  at `f32`.
- Microbenches of ~50 ns routines measure inlining, not arithmetic. Adjudicate a latency claim in
  the consumer's loop (interleaved A/B, `taskset`-pinned, one binary per variant), and `objdump`
  the wrapper before theorizing. Fewer ops can be slower; single-row ±3–5% moves from spill slots
  alone are not a finding.
- Any scripted edit (python/sed) asserts the old text exists; `cargo fmt` reflows silently.
