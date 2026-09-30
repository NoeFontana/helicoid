"""Generation-time cross-checks of coefficient records (docs/PHASE1.md section 4.3).

Every record's value and d_branch must agree to 100 digits with two computations that share
nothing with `mp.diff` on the definition: the committed rational Taylor series where the branch
variable is small, and a second formulation at 240 digits everywhere: the exact series of
`series.exact` (`k`, `a`...`e`), a closed form and its calculus derivative (`r`).
"""

from collections.abc import Callable
from functools import cache
from math import ceil

from mpmath import mp, mpf

from . import coeff, series
from .precision import DPS

DIGITS = 100
# The committed series truncate after `series.TERMS` terms; their relative truncation error at
# x <= SMALL is ~ TERMS * SMALL^(TERMS - 1) (`r`, the slowest), kept 20 digits under DIGITS.
SMALL = mpf(10) ** -ceil((DIGITS + 20) / (series.TERMS - 1))


class CrossCheckError(Exception):
    pass


@cache
def _mp(name: str, exact: bool) -> tuple:
    """The series as mpf at 240 digits: `series.exact`'s 232 terms, or the committed 16."""
    with mp.workdps(2 * DPS):
        rows = series.exact()[name] if exact else series.taylor(name)
        return tuple(mpf(q.numerator) / q.denominator for q in rows)


def _sum(rows: tuple, x) -> tuple:
    """(sum s_j x^j, sum j s_j x^(j-1)) by Horner."""
    value = deriv = mpf(0)
    for s in reversed(rows):
        deriv = deriv * x + value
        value = value * x + s
    return value, deriv


def _compare(name: str, inputs: dict, out: dict, want: tuple, how: str) -> None:
    for key, w in zip(("value", "d_branch"), want, strict=True):
        if abs(out[key] - w) > abs(w) * mpf(10) ** -DIGITS:
            raise CrossCheckError(f"{name} {key} vs {how} {inputs}: {out[key]} vs {w}")


def _r_closed(n, w) -> tuple:
    """r and dr/d(n^2) at fixed w, n > 0, by calculus on 2 atan2(n, w) / n."""

    def angle():
        return mp.atan(n / w) if w > 0 else mp.pi / 2

    def derivative(_):  # cancels ~ 1 / n^2: `stable` gives it the guard digits
        return w / ((w * w + n * n) * n * n) - angle() / n**3

    return 2 * angle() / n, coeff.stable(derivative, n * n)


def _check_r(inputs: dict, out: dict) -> None:
    n, w = mpf(inputs["n"]), mpf(inputs["w"])
    if w > 0 and (n / w) ** 2 <= SMALL:  # value (2/w) S(s), d/d(n^2) (2/w^3) S'(s), s = n^2/w^2
        v, d = _sum(_mp("r", False), (n / w) ** 2)
        _compare("r", inputs, out, (2 / w * v, 2 / w**3 * d), "Taylor")
    if n > 0:
        _compare("r", inputs, out, _r_closed(n, w), "closed form")


def coefficient(name: str) -> Callable[[dict, dict], None]:
    """The check of coefficient `name` for `FunctionSpec.check`: (inputs, out) -> None."""

    def check(inputs: dict, out: dict) -> None:
        with mp.workdps(2 * DPS):
            if name == "r":
                return _check_r(inputs, out)
            x = mpf(inputs["theta"]) ** 2
            if x <= SMALL:
                _compare(name, inputs, out, _sum(_mp(name, False), x), "Taylor")
            rows = _mp(name, True)[: 96 if x <= 1 else None]  # 96 terms: below 1e-150 at x <= 1
            _compare(name, inputs, out, _sum(rows, x), "exact series")

    return check
