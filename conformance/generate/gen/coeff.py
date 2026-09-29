"""Coefficient ids (docs/NUMERICS.md section 4), computed from their definitions.

`value` is the coefficient at the exact binary64 input. `d_branch` is its derivative with respect to
the branch variable (theta^2; n^2 at fixed w for `r`), taken by `mp.diff` at the exact real branch
value (not fl(theta * theta)). The definitions are the raw ones of the catalogue, never its
cancellation-free rewrites; a removable singularity takes its analytic extension.
"""

from mpmath import mp, mpf

from . import precision
from .strata import Q_W0_NORMS, Stratum


def theta_inputs(stratum: Stratum) -> list[dict]:
    return [{"theta": t} for t in stratum.thetas()]


def r_inputs(stratum: Stratum) -> list[dict]:
    """The unit quaternion of angle theta: (n, w) = (sin(theta/2), cos(theta/2)) in binary64.

    `q:w0` is a quaternion angle of exactly pi, which no binary64 theta is: (n, +0) at the norms
    of `Q_W0_NORMS`, since `Log` is scale-invariant (docs/NUMERICS.md section 3.2).
    """
    if stratum.name == "q:w0":
        return [{"n": n, "w": 0.0} for n in Q_W0_NORMS]
    return [
        {"n": precision.to_f64(mp.sin(mpf(t) / 2)), "w": precision.to_f64(mp.cos(mpf(t) / 2))}
        for t in stratum.thetas()
    ]


def _derive(f, x) -> dict:
    d = mp.diff(f, x)
    # A stencil straddling x = 0 reaches negative x, where `f` is complex-valued but real-analytic.
    if abs(mp.im(d)) > mpf(10) ** (10 - mp.dps):
        raise ArithmeticError(f"derivative has imaginary part {mp.im(d)}")
    return {"value": mp.re(f(x)), "d_branch": mp.re(d)}


def _with_branch_derivative(f, theta: float) -> dict:
    return _derive(f, mpf(theta) ** 2)


def _k(x):
    """k = sin(theta / 2) / theta as a function of x = theta^2; analytic extension 1/2 at 0."""
    if x == 0:
        return mpf(1) / 2
    s = mp.sqrt(x)
    return mp.sin(s / 2) / s


def k(inputs: dict) -> dict:
    return _with_branch_derivative(_k, inputs["theta"])


def stable(g, x):
    """`g(x)` to the caller's precision when the definition cancels.

    Two evaluations at growing guard precision until they agree to that precision. A result of
    exactly 0 counts as disagreement: no coefficient here vanishes, and a starved evaluation of a
    cancelling difference returns exactly 0 (at theta = 1e-310 every definition of `a`...`e` does
    at 120 digits; `e` loses 4 * 310 of them).
    """
    bits = mp.prec
    guard = 2 * max(0, -mp.mag(abs(x))) + 64
    while True:
        with mp.workprec(bits + guard):
            a = g(x)
        with mp.workprec(bits + 2 * guard):
            b = g(x)
        if a != 0 and abs(a - b) <= abs(b) * mpf(2) ** -(bits + 4):
            return +b
        guard *= 2


def _branch(g):
    """The coefficient with definition `g(theta)` as a function of x = theta^2."""

    def f(x):
        if x == 0:  # removable: the analytic extension, evaluated below the working precision
            x = mp.ldexp(mpf(1), -2 * mp.prec)
        return stable(lambda x: g(mp.sqrt(x)), x)

    return f


# The raw definitions of docs/NUMERICS.md section 4, one per coefficient (`k` is above).
DEFINITIONS = {
    "a": lambda t: (1 - mp.cos(t)) / t**2,
    "b": lambda t: (t - mp.sin(t)) / t**3,
    "c": lambda t: 1 / t**2 - (1 + mp.cos(t)) / (2 * t * mp.sin(t)),
    "d": lambda t: (t**2 + 2 * mp.cos(t) - 2) / (2 * t**4),
    "e": lambda t: (2 * t - 3 * mp.sin(t) + t * mp.cos(t)) / (2 * t**5),
}
BRANCH = {"k": _k, **{name: _branch(g) for name, g in DEFINITIONS.items()}}


def evaluator(name: str):
    return lambda inputs: _with_branch_derivative(BRANCH[name], inputs["theta"])


def r_branch(y, w):
    """r = 2 atan2(n, w) / n as a function of y = n^2 at fixed w >= 0; 2 / w at n = 0."""
    if y == 0:
        y = mp.ldexp(mpf(1), -2 * mp.prec)

    def g(y):
        n = mp.sqrt(y)
        # `atan2` is real-only. The stencil of `mp.diff` reaches y < 0 only for n ~ 0, w > 0, where
        # atan2(n, w) = atan(n / w) continues analytically. That complex `atan` is inexact at the
        # working precision (2e-26 relative at y = -1e-200), which left `d_branch` at ~73 digits
        # for n in 1e-150..1e-80: `stable` restores them.
        return 2 * (mp.atan2(n, w) if y > 0 else mp.atan(n / w)) / n

    return stable(g, y)


def r(inputs: dict) -> dict:
    n, w = mpf(inputs["n"]), mpf(inputs["w"])
    if n < 0 or w < 0 or n + w == 0:
        raise ValueError(f"r is defined for n >= 0, w >= 0, not both 0: {inputs}")
    return _derive(lambda y: r_branch(y, w), n * n)
