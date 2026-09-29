# 0022: `solve_cubic`'s `acos` and the limits it inherits

**Status:** draft
**Owner:** @NoeFontana
**Implementation:** —

## Context

`0017` step 2 ports omnisac's `solve_cubic` "as is". Two things in it do not survive a literal
port, and testing the port against planted roots found three limits that omnisac has too.

1. omnisac's three-root arm calls `f64::acos`, `f64::cos` and `core::f64::consts::PI`. `Real` has
   `sqrt`, `cbrt`, `sin_cos`, `atan2`, `abs` and `copysign` only, and no constant but `lit`.
2. Its tolerances are `0017` decision 3 (settled). Its limits are not written down anywhere.

## Decision (proposed)

1. `acos x = atan2(sqrt((1 - x)(1 + x)), x)` and `pi = atan2(+0, -1)`, private to `cubic.rs`. The
   factored product keeps the sine accurate where `x` is near `+-1` (a double root), where
   `sqrt(1 - x²)` loses everything. Evidence: within 2 ulp of `libm::acos` and `acosf` on 4e6
   points, the ends included; `pi` bit-equal to `f64::consts::PI` and `f32::consts::PI`; the roots
   of the trigonometric arm move against omnisac by 5 `u` of the largest root on the 4e4 random rows
   the golden rows are drawn from and 11 `u` at the most observed over 5e6 random polynomials (a
   sample maximum, not a bound); every other arm is bit-equal.
2. The limits below are documented in the rustdoc `# Domain` and pinned by tests, not fixed here.

## Rationale

`Real::acos` is the alternative: one more required method on every `Real` (as `cbrt` was), a
`Dual` rule `-d / sqrt(1 - x²)` singular at `+-1` like the form above, and bit identity with omnisac
in the trigonometric arm. It costs public surface for one caller, so the private form is taken
until the owner decides.

## Consequences

`PHASE2.md` §9's gate falls under its second clause: the trigonometric arm's values move by ulps
against omnisac, and a polynomial between the old and the new tolerance can change its root count.

## Implementation plan

1. `solve_cubic`, as landed by `0017` step 2, with the reading of decision 1.
2. If a limit below is to be fixed: a `NUMERICS.md` edit and a record for the formula, then code.

## Open questions

1. Accept the `atan2` form of `acos` and `pi`, or add `Real::acos`?
2. Which of these limits, all measured on the port and present in omnisac, get a fix?
   - **The one-real-root arm cancels `h - s`** when `|p|³ ≪ q²`: `x³ + p x - 1` is off by 4e-6 at
     `p` near 1e-5 (`f64`) and 4e-3 at `p` near 1e-2 (`f32`). A `Dual` root has derivative `-inf`
     for every `p` below 1.3e-5 (`f64`) although the root is simple, because `h - s` is exactly 0.
     One option: pair the cube roots by `u v = -p/3` instead of subtracting them.
   - **The band is relative to the summands of `disc`, not to `B`, `C`, `D`.** Where `p` and `q`
     cancel (roots close together, far from the origin) the rounding error of `disc` exceeds it:
     `(x - 1.09375)² (x - 0.921875)` at `f32` reads `disc > 0` and drops its double root. The
     planted-root test (`cubic_tests.rs`, `regular`) states where the error model holds.
   - **`p³` and `q²` underflow** (depressed roots below about 4e-8 at `f32`, 2e-54 at `f64`): `disc`
     and its band read 0, and `p < 0` sends the cubic to the trigonometric arm whatever its roots,
     so one real root is reported as three valid slots that are not roots
     (`x³ - 1e-16 x + 1e-24` at `f32`). Options: scale `p` and `q` by a power of two before forming
     the summands, or mask an arm whose summands underflowed.
   - **The leading floor bounds `|B|`, `|C|`, `|D|` by `2^17` at `f32`** (`2^46` at `f64`): the
     fixed multiple of `u` of `0017` decision 3. `x³ - 1e6 x` has no valid slot at `f32`.
