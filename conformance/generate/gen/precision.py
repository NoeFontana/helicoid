"""The one place the working precision is set (docs/PHASE1.md section 2, invariant 3)."""

import math
import struct

import mpmath
from mpmath import mp, mpf

DPS = 120
RECHECK_DPS = 150
RECHECK_DIGITS = 40


def setup() -> None:
    # gmpy2 would change mpmath's arithmetic backend; the corpus is pinned to the pure-Python one.
    if mpmath.libmp.BACKEND != "python":
        raise RuntimeError(f"mpmath backend is {mpmath.libmp.BACKEND!r}, expected 'python'")
    mp.dps = DPS


def signed_man_exp(x) -> tuple[int, int]:
    """x = m * 2**e exactly, m signed (`mpf.man_exp` drops the sign)."""
    m, e = x.man_exp
    return (-m if x < 0 else m), e


def to_f64(x) -> float:
    """Round an mpf to the nearest binary64 (ties to even, subnormals included).

    Integer true division is correctly rounded in CPython; `float(mpf)` is not relied on.
    """
    m, e = signed_man_exp(x)
    return float(m << e) if e >= 0 else m / (1 << -e)


F32_BITS = 24  # significand bits of a binary32
F32_QUANTUM = -149  # exponent of the spacing of the binary32 subnormals, 2^-149
F32_TOP = 128  # 2^128 is the first value beyond the largest binary32


def to_f32(x) -> float:
    """Round a float or an mpf to the nearest binary32 (ties to even, subnormals included),
    returned as the binary64 that holds it exactly (0016 item 1).

    Integer arithmetic on the exact value, like `to_f64`: the platform's `float32` is not relied
    on (a test compares them). A value beyond the binary32 range is an error, never an infinity.
    """
    m, e = signed_man_exp(mpf(x))
    n = abs(m)
    unit = max(n.bit_length() - 1 + e - (F32_BITS - 1), F32_QUANTUM)  # exponent of the last place
    if e >= unit:
        q = n << (e - unit)
    else:
        q, rest = divmod(n, 1 << (unit - e))
        half = 1 << (unit - e - 1)
        if rest > half or (rest == half and q & 1):
            q += 1
    if q.bit_length() + unit > F32_TOP:
        raise OverflowError(f"{x!r} is beyond the binary32 range")
    return math.ldexp(-q if m < 0 else q, unit)


def is_binary32(x: float) -> bool:
    """`x` is exactly a binary32: it survives the platform's round trip through one. The check
    that the corpus's `@f32` inputs need, made by a route independent of `to_f32`."""
    try:
        return struct.unpack("<f", struct.pack("<f", x))[0] == x
    except OverflowError:
        return False
