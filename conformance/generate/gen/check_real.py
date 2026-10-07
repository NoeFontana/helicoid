"""Generation-time cross-checks of the `real_*` ids (docs/decisions/0056).

- Every derivative against `mp.diff` at a relative step h = 10^-70 |x|, per argument, to 60
  digits. 0056 states 10^-40 |x|: the central difference's truncation, h^2 f'''/6, is then 1e-37
  of `sin`'s derivative at x = 1e22 and 1e-49 of `acos`'s at 1 - 2^-52, where the singularity is
  2.2e-16 away; at 10^-70 both are below 1e-89. The precision is raised until the difference's
  rounding, 10^-p |f| / h, is 60 digits under |f'| (`cos` at 1e-300 needs p ~ 760).
- Every value and derivative to 100 digits by an identity that ties it to the input or to the
  value: v^2 = x and 2 v d = 1 (`sqrt`), v^3 = x and 3 v^2 d = 1 (`cbrt`), sin^2 + cos^2 = 1,
  (d_sin, d_cos) = (cos, -sin) and the Maclaurin series for |x| <= 1 (`sin_cos`), the half-angle
  `atan2` form and d^2 (1 - x^2) = 1, d < 0 (`acos`), `atan` of y/x by quadrant,
  y d_y + x d_x = 0 and x d_y - y d_x = 1 (`atan2`), v d = n, d d_n = 1, d d_d = -v (`div`).
"""

from collections.abc import Callable

from mpmath import mp, mpf

from .precision import DPS

DIGITS = 100
DIFF_DIGITS = 60
DIFF_STEP = 70  # h = 10^-DIFF_STEP |x|


class CrossCheckError(Exception):
    pass


def _within(what: str, err, scale) -> None:
    if not abs(err) <= mpf(10) ** -DIGITS * abs(scale):
        raise CrossCheckError(f"{what}: {mp.nstr(err, 10)} vs {mp.nstr(scale, 10)}")


def against_diff(what: str, f: Callable, x, want) -> None:
    """`want` is f'(x) to 60 digits by `mp.diff` at h = 10^-70 |x|."""
    x = mpf(x)
    h = mpf(10) ** -DIFF_STEP * abs(x)
    with mp.workdps(30):
        ratio = abs(f(x)) / (h * abs(want))
        lost = max(0, int(mp.ceil(mp.log10(ratio)))) if ratio > 0 else 0
    with mp.workdps(max(DPS, lost + DIFF_DIGITS + 30)):
        got = mp.diff(f, x, h=h)
    if not abs(got - want) <= mpf(10) ** -DIFF_DIGITS * abs(want):
        raise CrossCheckError(f"{what}: {mp.nstr(want, 40)} vs mp.diff {mp.nstr(got, 40)}")


def _cbrt(x):
    return mp.cbrt(abs(x)) * (-1 if x < 0 else 1)


def real_sqrt(inputs: dict, out: dict) -> None:
    x, v, d = mpf(inputs["x"]), out["value"], out["d"]
    with mp.workdps(2 * DPS):
        _within("real_sqrt: v^2 - x", v * v - x, x)
        _within("real_sqrt: 2 v d - 1", 2 * v * d - 1, 1)
    against_diff("real_sqrt: d", mp.sqrt, x, d)


def real_cbrt(inputs: dict, out: dict) -> None:
    x, v, d = mpf(inputs["x"]), out["value"], out["d"]
    with mp.workdps(2 * DPS):
        _within("real_cbrt: v^3 - x", v**3 - x, x)
        _within("real_cbrt: 3 v^2 d - 1", 3 * v * v * d - 1, 1)
    against_diff("real_cbrt: d", _cbrt, x, d)


def _maclaurin(x) -> tuple:
    """(sin x, cos x) by their series, for |x| <= 1."""
    s, c, term, n = mpf(0), mpf(0), mpf(1), 0
    while True:
        c += term
        term *= x / (n + 1)
        s += term
        term *= -x / (n + 2)
        n += 2
        if abs(term) <= mp.eps * min(abs(s), abs(c)):
            return s, c


def real_sin_cos(inputs: dict, out: dict) -> None:
    x, s, c = mpf(inputs["x"]), out["sin"], out["cos"]
    with mp.workdps(2 * DPS):
        _within("real_sin_cos: sin^2 + cos^2 - 1", s * s + c * c - 1, 1)
        _within("real_sin_cos: d_sin - cos", out["d_sin"] - c, c)
        _within("real_sin_cos: d_cos + sin", out["d_cos"] + s, s)
        if abs(x) <= 1:
            ms, mc = _maclaurin(x)
            _within("real_sin_cos: sin vs series", s - ms, ms)
            _within("real_sin_cos: cos vs series", c - mc, mc)
    against_diff("real_sin_cos: d_sin", mp.sin, x, out["d_sin"])
    against_diff("real_sin_cos: d_cos", mp.cos, x, out["d_cos"])


def real_acos(inputs: dict, out: dict) -> None:
    x, v, d = mpf(inputs["x"]), out["value"], out["d"]
    with mp.workdps(2 * DPS):
        _within("real_acos: half-angle", v - 2 * mp.atan2(mp.sqrt(1 - x), mp.sqrt(1 + x)), v)
        _within("real_acos: d^2 (1 - x^2) - 1", d * d * (1 - x) * (1 + x) - 1, 1)
    if not d < 0:
        raise CrossCheckError(f"real_acos: d = {d} is not negative")
    against_diff("real_acos: d", mp.acos, x, d)


def real_atan2(inputs: dict, out: dict) -> None:
    y, x, v = mpf(inputs["y"]), mpf(inputs["x"]), out["value"]
    dy, dx = out["d_y"], out["d_x"]
    with mp.workdps(2 * DPS):
        t = mp.atan(y / x)
        want = t if x > 0 else t + (mp.pi if y >= 0 else -mp.pi)
        _within("real_atan2: atan by quadrant", v - want, want)
        _within("real_atan2: y d_y + x d_x", y * dy + x * dx, abs(y * dy) + abs(x * dx))
        _within("real_atan2: x d_y - y d_x - 1", x * dy - y * dx - 1, 1)
    against_diff("real_atan2: d_y", lambda t: mp.atan2(t, x), y, dy)
    against_diff("real_atan2: d_x", lambda t: mp.atan2(y, t), x, dx)


def real_div(inputs: dict, out: dict) -> None:
    n, d, v = mpf(inputs["n"]), mpf(inputs["d"]), out["value"]
    with mp.workdps(2 * DPS):
        _within("real_div: v d - n", v * d - n, n)
        _within("real_div: d d_n - 1", d * out["d_n"] - 1, 1)
        _within("real_div: d d_d + v", d * out["d_d"] + v, v)
    against_diff("real_div: d_n", lambda t: t / d, n, out["d_n"])
    against_diff("real_div: d_d", lambda t: n / t, d, out["d_d"])
