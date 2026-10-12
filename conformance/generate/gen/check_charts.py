"""Generation-time cross-checks of the SE(3) chart records (0060 decision 8) and the S2 ones (0063),
to 100 digits.

Each holds a property that fixes the output without the evaluator that produced it:

- `retract`: the chart's own `local` returns the twist (docs/maths/charts.md CH.4(a), theta < pi,
  so the `theta:pi-1e-k` strata, which stop strictly below pi, are inside it), and for `screw` the
  pose is X times `mp.expm` of the 4x4 hat matrix where `mp.logm`-free `mp.expm` is sound;
- `local`: the chart's own `retract` of the answer returns Y (CH.4(a) holds for every Y), the
  rotation compared after sign alignment and the translation relative to its own size.
"""

from collections.abc import Callable

from mpmath import mp, mpf

from . import charts, check, s2, sen3, so3
from .check_geodesic import _agree_pose


def _rel(what: str, got, want) -> None:
    got, want = list(got), list(want)
    check._agree(what, got, want, max(1, so3.maxabs(want), so3.maxabs(got)))


def retract(chart: str) -> Callable[[dict, dict], None]:
    def check_(inputs: dict, out: dict) -> None:
        back = charts.local(chart)({**inputs, "q1": out["q"], "x1": out["x"]})
        back = back["phi"] + back["rho"]
        _rel(f"se3_{chart}_retract: local round trip", back, [mpf(c) for c in inputs["tau"]])
        if chart == "screw":
            n0, r0 = charts._frame(inputs["q0"])
            E = sen3.expm([mpf(c) for c in inputs["tau"]], 1)
            x = [
                mpf(inputs["x0"][r]) + sum((r0[r][c] * E[c][3] for c in range(3)), mpf(0))
                for r in range(3)
            ]
            _rel("se3_screw_retract: expm translation", out["x"], x)

    return check_


def local(chart: str) -> Callable[[dict, dict], None]:
    def check_(inputs: dict, out: dict) -> None:
        tau = list(out["phi"]) + list(out["rho"])
        y = charts.retract(chart)({"q0": inputs["q0"], "x0": inputs["x0"], "tau": tau})
        want = {
            "q": so3.normalize([mpf(c) for c in inputs["q1"]]),
            "x": [mpf(c) for c in inputs["x1"]],
        }
        _agree_pose(f"se3_{chart}_local: retract round trip", y, want)

    return check_


# S2 (0063, docs/maths/charts.md CH.9): each check fixes the output by a route the evaluator does
# not take.


def s2_retract(inputs: dict, out: dict) -> None:
    """CH.9(a)'s closed form cos(theta) n - sin(theta)/theta B K delta (the oriented frame, 0082), the unit norm, and
    `local` of the answer returning delta (CH.9(d): every stratum has theta < pi)."""
    nh, _, b1, b2 = s2.basis(inputs["n"])
    d1, d2 = (mpf(c) for c in inputs["delta"])
    theta = mp.sqrt(d1 * d1 + d2 * d2)
    sinc = mp.sin(theta) / theta if theta != 0 else mpf(1)
    bk = [-d2 * a + d1 * b for a, b in zip(b1, b2, strict=True)]  # B K delta
    closed = [mp.cos(theta) * x - sinc * y for x, y in zip(nh, bk, strict=True)]
    m = list(out["m"])
    check._agree("s2_retract: CH.9(a) closed form", m, closed)
    check._agree("s2_retract: |m|", [s2.norm(m)], [mpf(1)])
    back = s2.local({"n": inputs["n"], "m": m})["delta"]
    check._agree("s2_retract: local round trip", back, [d1, d2])


def s2_local(inputs: dict, out: dict) -> None:
    """`retract` of the answer returning m (CH.9(d), every m != -n), and 0063 decision 4's route:
    the quaternion (cos(alpha/2), sin(alpha/2) v/s) through `so3.log` gives the rotation vector,
    and its rotation maps n to m."""
    delta = list(out["delta"])
    p, nh, mh = s2.rotation_vector(inputs["n"], inputs["m"])
    check._agree("s2_local: retract round trip", s2._retract(inputs["n"], delta), mh)
    alpha = s2.norm(p)
    if alpha == 0:
        return
    axis = [c / alpha for c in p]
    q = [mp.cos(alpha / 2), *(mp.sin(alpha / 2) * c for c in axis)]
    check._agree("s2_local: so3.log of the minimal rotation", so3.log({"q": q})["phi"], p)
    R = so3.rot(q)
    check._agree(
        "s2_local: the minimal rotation maps n to m",
        [sum(R[i][k] * nh[k] for k in range(3)) for i in range(3)],
        mh,
    )
