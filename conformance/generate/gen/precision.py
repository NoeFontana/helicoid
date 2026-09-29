"""The one place the working precision is set (docs/PHASE1.md section 2, invariant 3)."""

import mpmath
from mpmath import mp

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
