# 0115: The spline tables are exact integers

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** —

## Context

[`0114`](./0114-a-spline-evaluates-the-window-it-is-handed.md)'s `Weights::at` needs the cumulative basis
`M̃` of order `K`. An Rⁿ spline's derivative energy needs the Gram matrices `G⁽ʳ⁾`. `PHASE7.md` §2
asks for the basis "computed at compile time". [`0004`](./0004-switch-points-are-generated-not-typed.md)
has measured numbers written by a generator; these are not measured.

Both tables are rational with small denominators. In exact arithmetic (Python `fractions`, every
order from 2 to 6):

| `K` | basis denominator `(K−1)!` | largest basis numerator | `G⁽¹⁾` denominator, largest numerator | `G⁽²⁾` denominator, largest numerator |
|---|---|---|---|---|
| 2 | 1 | 1 | 1, 1 | — |
| 3 | 2 | 2 | 6, 2 | 1, 4 |
| 4 | 6 | 6 | 120, 34 | 6, 6 |
| 5 | 24 | 24 | 5 040, 950 | 120, 126 |
| 6 | 120 | 120 | 362 880, 47 248 | 5 040, 1 832 |

Every number is below `2²⁴`, so each is exact in binary32 and binary64. Check values at `h = 1`:
`G⁽¹⁾₀₀ = 1/20` at `K = 4` and `G⁽²⁾₀₀ = 1/252` at `K = 6`, and every row of `G⁽ʳ⁾`, `r ≥ 1`,
sums to zero. The K = 4 rows of `M̃` are `(6, 0, 0, 0)`, `(5, 3, −3, 1)`, `(1, 3, 3, −2)` and
`(0, 0, 0, 1)`, over 6.

## Decision

1. **`const fn`s in `helicoid-spline` compute the integer numerators** at compile time, in `i64`:
   - the basis `m_{s,n} = C(K−1, n) Σ_{l=s}^{K−1} (−1)^{l−s} C(K, l−s) (K−1−l)^{K−1−n}` (an
     integer before the division by `(K−1)!`);
   - its cumulative sums `m̃_{j,n} = Σ_{s≥j} m_{s,n}`;
   - the Gram numerators over their least common denominator.

   No generator, no generated file and no check recipe: a `const fn` cannot drift from itself.
2. **One rounding per entry.** A Gram entry is `S::lit(num) / S::lit(den)`: two exact conversions
   and one correctly rounded division, so binary32 never rounds through binary64. A weight is Horner
   on the `S::lit` numerators, then one division by `S::lit((K−1)!)`.
3. **Each `h` scaling is a separate multiply.**
   - `λ̇` multiplies by `inv_h` once.
   - `λ̈` multiplies by `inv_h * inv_h`, formed once.
   - `G⁽ʳ⁾` multiplies by `inv_h^{2r−1}`.
   - A power-of-two `inv_h` makes each of these exact, and the Gram corpus id is drawn on such `h`
     ([`0116`](./0116-the-spline-is-measured-before-it-ships.md)).
4. **The orders are `ValidOrder`'s**: the set [`0113`](./0113-phase-7-opens-for-a-calibration-consumer.md)
   settles, `K = 2` to `6`, and no other. A `const` assertion in each table fails the build
   if any numerator reaches `2²⁴`, so a later order cannot silently lose exactness at binary32.
5. **The public Gram function** is
   `pub fn derivative_gram<S: Real, const K: usize, const R: usize>(inv_h: S) -> Matrix<S, K, K>`,
   for `R` = 1 or 2, asserted at monomorphization. It is the matrix of the energy
   `∫ ‖p⁽ʳ⁾(t)‖² dt` over one segment, in the segment's `K` controls:
   `G⁽ʳ⁾ = h^{1−2r} M D_r H D_rᵀ Mᵀ`, `H_{ab} = 1/(a+b+1)`, with `b(u) = M [1, u, …]ᵀ` (rows are
   controls). A weighted energy `Q⁻¹` is the caller's Kronecker product.

## Rationale

- **A generator plus a check recipe** (the first sketch) costs a Python path, a generated file, a
  `just` recipe and a CI step to guard numbers that a `const fn` computes exactly with nothing to
  guard. `0004`'s generator exists because switch points are measured; these are combinatorics.
- **Literal tables** typed into the source would be the hand-typed coefficients CLAUDE.md forbids.
- **Rationals in `f64` then cast to `f32`** would round twice.
- **Independent check**: the mpmath corpus id `rn_spline_gram` (`0116`) at a 0-ulp bar, and a unit
  test against exact rationals for the basis, so the `const fn` is not its own oracle.

## Consequences

- `PHASE7.md` §2's "computed at compile time" holds as written.
- `NUMERICS.md`'s spline section states the integer forms and the one-division rule.

## Implementation plan

1. The `const fn`s and `derivative_gram`, in `0114`'s step 1. Verified by `weights_match_exact_rationals`
   and `gram_matches_exact_rationals` (exact `i128` rational arithmetic in the test, every order,
   both `r`), and `rn_spline_gram` at 0 ulp once `0116`'s corpus exists.

## Open questions

None.
