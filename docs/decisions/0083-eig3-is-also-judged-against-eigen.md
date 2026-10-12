# 0083: `eig3` is also judged against Eigen

**Status:** ready
**Owner:** @NoeFontana
**Implementation:** branch `oracle/eigen-eig3`

## Context

`eig3` has one oracle, nalgebra's `SymmetricEigen` ([`0056`](./0056-the-routines-d7-does-not-reach.md)
decision 4, [`0077`](./0077-an-oracle-answers-at-the-precision-it-is-asked.md)), and
[`0081`](./0081-eig3-is-four-jacobi-sweeps.md) made it four Jacobi sweeps that dominate nalgebra on
all 50 strata. A C++ consumer's symmetric 3 × 3 solver is Eigen's
`SelfAdjointEigenSolver<Matrix3<T>>`: Sophus, manif and GTSAM all sit on Eigen. The Sophus image
([`0072`](./0072-a-container-oracle-is-an-image-on-the-same-protocol.md) decision 2) already pins
Eigen 3.4.0 from its archive snapshot, compiled with `-ffp-contract=off -DEIGEN_DONT_VECTORIZE`,
and its runner already answers in binary32 (`0077`).

Eigen offers two routes: `compute` (a Householder tridiagonalisation, then implicit symmetric QR
steps until the off-diagonal is negligible) and `computeDirect` (the trigonometric closed form,
vectors by cross products: the kind of solver `0081` retired).

## Decision

1. **The Sophus runner answers `eig3`** with `Eigen::SelfAdjointEigenSolver<Matrix3<T>>(A)`,
   eigenvectors computed: `T = double` on every record, `T = float` on the `@f32` records.
   `Runner::answers` and `answers_f32` gain `eig3`. The subject stays `sophus`; its pin already names
   `eigen@3.4.0`, so `subject_version` does not change. The image, its flags and its other ids are
   unchanged.
2. **The conversion** (`convert.hpp`, `helicoid_to_eigen_matrix3`, `eigen_to_helicoid_eig3`): the
   corpus's column-major `A` is Eigen's storage order; `lambda` ascending (Eigen sorts, and a stable
   sort states it, as nalgebra's runner does); `V` the matching columns, column-major, each
   column's sign Eigen's (the `eigen` metric aligns each column). An `info()` other than `Success`
   is a refusal, `nan` (`0072` decision 5); the corpus has none.
3. **The hand test** (`convert_test.cpp`, `eig3`): `A`'s storage order; `A = [5 0 0; 0 2 1; 0 1 2]`,
   whose largest eigenvalue comes first on the diagonal, gives `λ = (1, 3, 5)` and
   `V`'s columns `(0, 1, −1)/√2`, `(0, 1, 1)/√2`, `e₀` up to sign, so an unsorted or a transposed
   `V` fails; the same at `float`, every output a binary32.
4. **`computeDirect` is not a subject.** A runner answers an id once, and a second answer would
   need a second runner name for one id of one image. Its rows are measured once, below, by a
   scratch build of this runner calling `computeDirect`.

## Rationale

- **The Sophus image, not a fourth one**: the same base, snapshot, compiler, flags and hex-float
  codec would be copied for one id. The price is a subject called `sophus` answering with Eigen,
  which this record and the pin state.
- **`compute`, not `computeDirect`**: it is what `SelfAdjointEigenSolver(A)` does, and the closed
  form loses to `helicoid` on every stratum (below), so as an oracle it would say nothing.
- **The image's flags** (`0072`'s rule): no contraction and no Eigen SIMD, one program on every
  host. Speed under a consumer's flags is measured below as context, not as a bar.

## Consequences

Measured on `0f5628e` (`fix/eig3-accuracy`) + this branch: `just oracle-sophus`,
`cargo xtask conformance --fn eig3 [--precision f32]` and `--oracle nalgebra --fn eig3`; Eigen's
`computeDirect` by the scratch runner of decision 4, scored by the same harness. `max_u` per
stratum:

| stratum | `helicoid` | Eigen `compute` | nalgebra | Eigen `computeDirect` |
|---|---|---|---|---|
| `eig:random` | 3.53 | 7.72 | 8.19 | 13.0 |
| `eig:gap-1e-{0,1}/{bottom,top}` | 2.30–2.80 | 5.94–7.67 | 6.76–8.78 | 3.21–5.37 |
| `eig:gap-1e-2/{bottom,top}` | 2.01–4.76 | 8.10–9.53 | 5.88–8.67 | 26.6–41.4 |
| `eig:gap-1e-{3,…,7}/{bottom,top}` | 1.74–2.84 | 5.20–10.1 | 4.90–14.9 | 193–3.67e6 |
| `eig:gap-1e-{8,…,12}/{bottom,top}` | 1.08–1.93 | 3.80–6.04 | 6.28–12.8 | 2.04e7–6.75e7 |
| `eig:triple` | 3.99e-14 | 3.99e-14 | 3.99e-14 | 3.99e-14 |
| `eig:rank1` | 1.66 | 5.65 | 13.3 | 5.93e7 |
| `eig:scale-up`, `-down` | 3.07, 2.30 | 10.0, 10.8 | 7.95, 7.57 | 26.8, 12.6 |
| `eig:random@f32` | 2.36 | 6.73 | 6.15 | 21.6 |
| `eig:gap-1e-{0,1}/{bottom,top}@f32` | 1.79–2.11 | 7.15–9.13 | 6.27–6.79 | 3.67–6.45 |
| `eig:gap-1e-2/{bottom,top}@f32` | 1.58–1.96 | 6.15–6.50 | 6.65–6.80 | 29.2–40.6 |
| `eig:gap-1e-{3,…,6}/{bottom,top}@f32` | 1.39–1.82 | 3.39–6.68 | 6.29–10.4 | 269–3.05e3 |
| `eig:triple@f32` | 0 | 0 | 0 | 0 |
| `eig:rank1@f32` | 1.99 | 4.23 | 12.8 | 3.08e3 |
| `eig:scale-up@f32`, `-down@f32` | 1.99, 3.24 | 7.42, 7.11 | 8.00, 8.46 | 7.83, 14.7 |

- **Pairing**: the 31 binary64 and 19 binary32 `eig3` strata gain a second oracle. No record is
  refused; no row is non-finite.
- **Domination holds against Eigen on all 50**, `eig:triple` a tie at both precisions; the ratio of
  maxima is 0.19–0.59 at binary64 and 0.21–0.51 at binary32. No domination failure is added.
  Eigen is the better oracle of the two on 20 of the 30 non-triple binary64 strata and 13 of 18 at
  binary32, so it sets most of `eig3`'s bars.
- **`computeDirect`** loses to `helicoid` on every non-triple stratum, as the retired closed form
  did: `1/gap` growth to `6.75e7 u` at binary64 and `3.08e3 u` at binary32.
- **Speed, as context.** Scratch harnesses (not committed): `benches/linalg.rs`'s three fixtures,
  bit for bit, and 1024 random symmetric matrices in turn; each call behind a non-inlined wrapper,
  every output escaped; minimum ns per call over 31 runs of about 20 ms; one pinned core of an AMD
  EPYC-Milan VM (load average 3); Eigen 3.4.0 under g++ 13.3 `-O3 -DNDEBUG`; `helicoid` and
  nalgebra 0.35.0 at the release profile. `criterion` on `benches/linalg.rs` reads `helicoid` at
  325, 364, 326 / 307, 279, 278 ns.

  | `f64` / `f32`, ns | random | `gap-1e-8` | triple | 1024 random |
  |---|---|---|---|---|
  | `helicoid` | 320 / 303 | 361 / 273 | 321 / 273 | 322 / 310 |
  | nalgebra | 411 / 360 | 258 / 129 | 160 / 121 | 388 / 304 |
  | Eigen `compute`, image flags | 304 / 263 | 159 / 50 | 22 / 15 | 281 / 200 |
  | Eigen `compute`, `-O3` | 355 / 300 | 181 / 89 | 56 / 44 | 324 / 236 |
  | Eigen `compute`, `-O3 -march=x86-64-v3` | 351 / 295 | 179 / 83 | 54 / 39 | 329 / 233 |
  | Eigen `computeDirect`, image flags | 90 / 77 | 114 / 58 | 93 / 76 | 92 / 87 |
  | Eigen `computeDirect`, `-O3` | 116 / 84 | 96 / 76 | 97 / 80 | 98 / 84 |

  On a random spectrum `compute` takes 0.95–1.11× `helicoid`'s time at `f64` and 0.87–0.99× at
  `f32` over the three flag sets (clang 18 reads 0.91–1.28× g++'s times); where its iteration stops
  early (a near pair, a multiple of the identity) it is 2–18× faster, as nalgebra is (`0081`
  *Cost*). `computeDirect` is 2.8–4.7× faster than `helicoid` and loses every stratum.

## Implementation plan

1. One PR: this record, the runner, `convert.hpp`, `convert_test.cpp`, the registry, the recipe's
   comment, `PHASE1.md` §0.0 and §7 — verified by the image build (it runs `convert_test`),
   `cargo nextest run -p xtask -- oracle`, `just oracle-sophus`, `just lint`, `just test`,
   `just doc`.

## Open questions

None.

## Further work

- An early exit for `eig3` where the off-diagonal is already negligible: Eigen's `compute` reads
  22 ns on `eig:triple`'s fixture against `helicoid`'s 321. A convergence test is a branch on
  data and moves bits, so it is a record of its own against these rows.
