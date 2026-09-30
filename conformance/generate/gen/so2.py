"""SO(2) ids (docs/PHASE1.md section 4.3, `so2_*`), computed from their definitions.

`so2_exp` is `mp.expm` of the 2x2 hat matrix [[0, -theta], [theta, 0]], `so2_log` the inverse of
the exponential on theta in (-pi, pi] by Newton's method on its power series. Neither `atan2` nor
`cos`/`sin` of docs/NUMERICS.md section 6 is evaluated; a test compares them to the corpus. A
unit complex is `z = [c, s]`, a rotation angle is `theta`, both exact binary64 in the inputs.
"""

import math

from mpmath import mp, mpc, mpf

from . import so3
from .precision import DPS, to_f64
from .strata import Stratum


def hat(theta):
    o = mpf(0)
    return [[o, -theta], [theta, o]]


def expm(theta):
    E = mp.expm(mp.matrix(hat(theta)))
    return [[E[i, j] for j in range(2)] for i in range(2)]


def exp_series(theta) -> mpc:
    """sum (i theta)^n / n!: the complex form of the exponential, a second computation of `exp`."""
    p, eps = mpc(0, theta), so3._eps()
    term = total = mpc(1, 0)
    n = quiet = 0
    while quiet < 2:  # even terms are real, odd terms imaginary: each part against its own sum
        n += 1
        if n > so3.MAX_ITER:
            raise ArithmeticError("exp series does not converge")
        term = term * p / n
        total += term
        part = (lambda w: w.real) if n % 2 == 0 else (lambda w: w.imag)
        quiet = quiet + 1 if abs(part(term)) <= eps * abs(part(total)) else 0
    return total


def angle(z) -> mpf:
    """theta with |theta| <= pi and exp(i theta) = z / |z|, by Newton's method on the series.

    Each step is theta += Im(conj(exp(i theta)) z), the sine of the remaining angle: the error
    cubes. The precision follows the accuracy reached, the float seed only starts the iteration,
    and it stops when the residual, at the full precision, is below `so3._eps()` of |theta|: the
    series of exp near pi has an absolute error of 1e-120, which sin(theta) ~ 1e-12 cannot bear.
    z = (-1, 0) is refused, never sampled: `atan2` gives pi at s = +0 and -pi at s = -0, while
    (-pi, pi] gives pi at both, and an mpf has no signed zero to tell the two apart.
    """
    z = so3.normalize(z)
    if z[1] == 0:
        if z[0] < 0:
            raise ValueError("z = (-1, 0): an mpf cannot carry the sign of the zero")
        return mpf(0)
    theta = mpf(math.atan2(float(z[1]), float(z[0])))
    digits, full, eps = 14, mp.dps, so3._eps()
    for _ in range(20):
        wp = min(full, 3 * digits + 15)
        with mp.workdps(wp):
            E = exp_series(theta)
            r = E.real * z[1] - E.imag * z[0]  # Im(conj(E) z)
            if wp == full and abs(r) <= eps * abs(theta):
                return theta
            theta += r
        digits = min(3 * digits - 2, wp - 5)
    raise ArithmeticError("angle does not converge")


def exp(inputs: dict) -> dict:
    E = expm(mpf(inputs["theta"]))
    return {"z": [E[0][0], E[1][0]]}


def log(inputs: dict) -> dict:
    return {"theta": angle([mpf(c) for c in inputs["z"]])}


def signed(thetas: list[float]) -> list[float]:
    """Every theta followed by its negative; 0 has none."""
    return [s for t in thetas for s in ((t, -t) if t else (t,))]


def z_of(theta: float) -> list[float]:
    """The unit complex of an angle, rounded componentwise to binary64."""
    with mp.workdps(DPS):
        return [to_f64(mp.cos(mpf(theta))), to_f64(mp.sin(mpf(theta)))]


def theta_inputs(stratum: Stratum) -> list[dict]:
    return [{"theta": t} for t in signed(stratum.thetas())]


def z_inputs(stratum: Stratum) -> list[dict]:
    return [{"z": z_of(t)} for t in signed(stratum.thetas())]
