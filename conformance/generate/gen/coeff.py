"""Coefficient ids (docs/NUMERICS.md section 4), computed from their definitions.

`value` is the coefficient at the exact binary64 theta. `d_branch` is its derivative with respect to
the branch variable theta^2, taken by `mp.diff` at the exact real theta^2 (not fl(theta * theta)).
"""

from mpmath import mp, mpf

from .strata import Stratum


def theta_inputs(stratum: Stratum) -> list[dict]:
    return [{"theta": t} for t in stratum.thetas()]


def _with_branch_derivative(f, theta: float) -> dict:
    x = mpf(theta) ** 2
    d = mp.diff(f, x)
    # A stencil straddling x = 0 reaches negative x, where `f` is complex-valued but real-analytic.
    if abs(mp.im(d)) > mpf(10) ** (10 - mp.dps):
        raise ArithmeticError(f"derivative has imaginary part {mp.im(d)}")
    return {"value": f(x), "d_branch": mp.re(d)}


def _k(x):
    """k = sin(theta / 2) / theta as a function of x = theta^2; analytic extension 1/2 at 0."""
    if x == 0:
        return mpf(1) / 2
    s = mp.sqrt(x)
    return mp.sin(s / 2) / s


def k(inputs: dict) -> dict:
    return _with_branch_derivative(_k, inputs["theta"])
