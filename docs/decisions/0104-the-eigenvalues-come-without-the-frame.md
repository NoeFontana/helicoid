# 0104: The eigenvalues come without the frame

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** (none yet)

## Context

`helicoid-linalg` has a crate-private `eig3_values(a: &Mat3<S>) -> Vec3<S>`: `eig3`'s Jacobi
without the accumulated rotation, then `eig3`'s sort and scaling on the values alone. `0102`'s
review added it for `nearest_rotation`, whose arm test reads only the eigenvalues. Its rustdoc
promises `eig3`'s eigenvalues "bit for bit"; no test checks that.

A consumer summarises a 6 × 6 pose covariance by $\sqrt{\lambda_{\max}}$ of each 3 × 3 diagonal
block, through `eig3(&block).0`: two calls per summary, each discarding the frame. The public
surface gives no cheaper path.

**Cost.** `crates/helicoid/benches/linalg.rs`'s `eig3` fixtures, with `eig3_values/*` rows added
and the function made `pub` in a scratch tree at `c65da5e`; one pinned core, three runs of the
binary, criterion's estimate in ns per call. The host was shared (load 4–5); the three runs
agree within 2 ns except `f32` `gap-1e-8` (39–44).

| Fixture | `f64` `eig3` | `f64` values | ratio | `f32` `eig3` | `f32` values | ratio |
|---|---|---|---|---|---|---|
| `random` | 228 | 188 | 0.82 | 161 | 125 | 0.78 |
| `gap-1e-8` | 140 | 99 | 0.71 | 69 | 40 | 0.58 |
| `triple` | 85 | 46 | 0.54 | 56 | 26 | 0.46 |
| `mixed-three-calls` | 454 | 333 | 0.73 | 288 | 193 | 0.67 |

**Bit identity.** A scratch proptest compared `eig3_values(a)` with `eig3(a).0`, every lane's bits,
a NaN counted as one value, on 200 000 draws of `eig3_tests.rs`'s `entry()` (ties, signed zeros
and non-finite entries included) at `f64`, `f32` and `Dual<f64, 2>` (derivative lanes from a
second draw): no difference.

## Decision

1. **`helicoid_linalg::eig3_values<S: Real>(a: &Mat3<S>) -> Vec3<S>`** is public, re-exported at the
   crate root beside `eig3`. Its contract is `eig3(a).0`, bit for bit at every `S: Real`, a NaN
   being one value: the same ascending order, the same lower-triangle read, the same `NaN`
   eigenvalues for a non-finite entry. It has no domain to assert (R6).
2. **It is a §14 row, not a corpus id.** `NUMERICS.md` §14 gains the row: fast `eig3_values`
   (`helicoid-linalg`), twin the composition `eig3(a).0`, bit-identical, in `0019`'s form for a
   composite of public functions. The proptest `eig3_values_matches_reference` asserts the bits at
   `f64`, `f32` and `Dual`, over `entry()`.
3. **Surface text.** `API.md` §2's `eig3` row names `eig3_values`; R6's list of functions that
   assert nothing adds it. `PHASE2.md` §6 gains the signature under `eig3`'s bullet. The rustdoc
   states the contract against `eig3` and, in an example, reads the largest eigenvalue.
4. **Bench rows** `eig3_values/{random,gap-1e-8,triple,mixed-three-calls}` join
   `benches/linalg.rs` at `f64` and `f32`.

## Rationale

`API.md` §6, line by line:

1. **Rules.** It returns values in `eig3`'s stated order and no mask or `bool` (R3, R4), and has
   no domain (R6).
2. **Verification.** A composite with a bit-identical twin; the corpus id `eig3` covers both.
3. **Domain.** None, as `eig3`.
4. **Return.** A `Vec3` in the order `PHASE2.md` §6 states; no `bool`, no Jacobian.
5. **A second spelling of `eig3(a).0`.** Yes, and the table is why it earns one: the frame costs
   18–54 % of `eig3`'s time, and the consumer discards it at every call. `inv_act` (`0103`) and
   `mul_inv` (`0048`) are named paths beside a composition that stays their twin; this one is
   stricter, since the twin is the same bits.
6. **Ownership.** `eig3` is `helicoid-linalg`'s already (`PHASE2.md` §6).
7. **Budget.** The same kernel as `eig3`: `no_std`, no allocation, `libm` only, bit-identical
   across targets wherever `eig3` is.
8. **A way to be wrong.** The consumer's cheap alternative is a characteristic-cubic eigenvalue,
   the closed form `0081` replaced because its `acos` lost $\sqrt u\,\lVert A\rVert$ at a close
   pair. A values-only path with `eig3`'s accuracy removes the reason to write one.

**Alternatives.**
- *Document `eig3(a).0` and add nothing* (§6 question 5's default): the consumer pays 1.2–2.2×
  per call for a frame it drops.
- *`eig3_max`, the largest eigenvalue alone*: Jacobi converges all three together, so it saves
  only the sort's two exchanges, and it would be a third spelling.
- *Expose the kernel* (`jacobi::<S, false>`): it returns the unsorted, scaled diagonal and a
  finiteness mask, and every caller would rewrite the sort and the scaling back.

## Consequences

- `eig3` and `eig3_values` stay bit-identical: a change to `eig3`'s eigenvalue path (the exit,
  the scaling, the sort) is a change to both, and `eig3_values_matches_reference` fails if one
  moves alone.
- No bar changes: the corpus id `eig3` and its strata score the values both functions return.
- A public item at `0.0.x`; removing it later is a breaking change.

## Implementation plan

1. `pub` and the re-export, the rustdoc, `eig3_values_matches_reference`, the §14 row, the
   `API.md` §2 and R6 text, the `PHASE2.md` §6 signature, the bench rows and a changelog entry —
   verified by `cargo nextest run -p helicoid-linalg -- eig3_values_matches_reference` (`dev` and
   `ndebug`), `just lint` (the twin table), `just doc`, `just test`.

## Open questions

None.

## Further work

None.
