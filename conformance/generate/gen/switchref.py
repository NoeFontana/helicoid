"""The true coefficient at every switch the sweep may choose (docs/decisions/0039 plan step 0).

docs/maths/coefficients.md CO.12 bounds the jump between a coefficient's two arms at its switch by
the sum of the arms' errors *there*. Sampling that right-hand side at the two corpus records
bracketing the switch underestimates it: the exact arm's error is a sawtooth, swinging 195x over
0.8% of theta, so no density of corpus records makes a point sample a bound (0039 refuted both
cheap repairs -- a denser `theta:dense` and a per-stratum maximum, the second vacuous for the
series arm at 6.3e6 u). The repair is a reference *at* the switch, and this file is it.

A switch is a point of the grid docs/PHASE1.md section 6 fixes, so a reference at every grid point
is a reference at every switch the sweep can return -- and it stays a function of the grid alone,
never of the sweep, which is what keeps `just corpus-check` independent of `just thresholds`.

The grid is rebuilt here by the same integer 64th-root construction as
xtask/src/thresholds/grid.rs, at both precisions, and a record is keyed by the point's **bit
pattern**. A divergence between the two constructions is then a missing key and a loud failure in
the test that reads it, never a silent mismatch of two nearly equal numbers.
"""

import struct

from mpmath import mp, mpf

from . import coeff, fmt, precision

# docs/PHASE1.md section 6, as xtask/src/thresholds/grid.rs builds it.
PER_DECADE = 64
DECADES_BELOW = 16
DECADES_ABOVE = 1
ROOT_PRECISION = 128  # `P` in grid.rs: the integer root is taken in units of 2^-P

#: The coefficients with a generated switch, in the order `xtask` sweeps them.
NAMES = ("k", "a", "b", "c", "d", "e", "cos_half", "r")

#: `r`'s branch variable is `n^2/w^2`, and its value scales as `1/w`: the arms this reference is
#: compared against are evaluated at `w = 1` (`coeffs::tests::arms`), so the reference is too.
R_W = 1


class SwitchRefError(Exception):
    pass


def _iroot64(n: int) -> int:
    """The integer 64th root of `n`, by bisection on the bit length: exact, and no float anywhere."""
    if n < 1:
        raise SwitchRefError(f"no 64th root of {n}")
    lo, hi = 1, 1 << (n.bit_length() // 64 + 1)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid**64 <= n:
            lo = mid
        else:
            hi = mid - 1
    return lo


def grid_point(e: int, bits: int) -> float:
    """`10^(e/64)` rounded to nearest at `bits` significand bits, held as a binary64.

    `m = floor(2^P * 10^(e/64))` is the integer 64th root of `2^(64P) * 10^e` (the floor of a root
    is the root of the floor); `(2m + 1) / 2^(P + 1)` lies strictly between `m` and `m + 1` in
    units of `2^-P`, so the one rounding is never a tie.
    """
    unit = 1 << (64 * ROOT_PRECISION)
    scaled = unit * 10**e if e >= 0 else unit // 10 ** (-e)
    m = _iroot64(scaled)
    with mp.workprec(64 * ROOT_PRECISION):
        exact = mpf(2 * m + 1) / mpf(1 << (ROOT_PRECISION + 1))
    return precision.to_f64(exact) if bits == 53 else precision.to_f32(exact)


def grid(bits: int) -> list[float]:
    """The switch points at `bits` significand bits, increasing, as `grid.rs` returns them."""
    lo, hi = DECADES_BELOW * PER_DECADE, DECADES_ABOVE * PER_DECADE
    return [grid_point(e, bits) for e in range(-lo, hi + 1)]


def _bits(x: float) -> int:
    """The binary64 bit pattern of `x`, as Rust's `f64::to_bits` gives it."""
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def _at(name: str, z: float) -> dict:
    """The true value and `d/dz` of `name` at the branch variable `z`, at 150 digits."""
    if name == "r":
        got = coeff._derive(lambda y: coeff.r_branch(y, R_W), mpf(z))
    else:
        got = coeff._derive(coeff.BRANCH[name], mpf(z))
    for field, x in got.items():
        if not mp.isfinite(x):
            raise SwitchRefError(f"{name} at z = {z!r}: {field} is {x}")
    return got


def build() -> tuple[bytes, int, int]:
    """(file bytes, record count, checked count): one record per coefficient, precision and point.

    Every record is checked by being recomputed at a higher precision and agreeing to
    `precision.RECHECK_DIGITS`, which is what makes a 30-digit reference a reference and not an
    evaluation: the cancelling definitions of `coeff.DEFINITIONS` lose digits at small `z`, and
    `coeff.stable` is the only thing standing between that and a wrong committed number.
    """
    lines, checked = [], 0
    for bits, name_of in ((53, "f64"), (24, "f32")):
        points = grid(bits)
        for name in NAMES:
            for z in points:
                precision.setup()
                got = _at(name, z)
                mp.dps = precision.RECHECK_DPS
                again = _at(name, z)
                precision.setup()
                for field in ("value", "d_branch"):
                    a, b = got[field], again[field]
                    if a != 0 and abs(a - b) / abs(b) > mpf(10) ** -precision.RECHECK_DIGITS:
                        raise SwitchRefError(f"{name} {name_of} at {z!r}: {field} not reproducible")
                checked += 1
                lines.append(
                    fmt.dumps(
                        {
                            # The key a reader matches on: the bit pattern of the binary64 that
                            # holds the point. An f32 point is the binary64 it equals (0016 item
                            # 3), so one key serves both precisions and a grid that drifts apart
                            # from grid.rs is a missing key, not a near miss.
                            "bits": f"{_bits(z):016x}",
                            "coeff": name,
                            "d_branch": fmt.dec30(got["d_branch"]),
                            "precision": name_of,
                            "value": fmt.dec30(got["value"]),
                            "z": fmt.hex_float(z),
                        }
                    )
                    + "\n"
                )
    return "".join(lines).encode(), len(lines), checked
