"""Generation-time cross-checks of the SE(3) chart records (0060 decision 8), to 100 digits.

Each holds a property that fixes the output without the evaluator that produced it:

- `retract`: the chart's own `local` returns the twist (docs/maths/charts.md CH.4(a), theta < pi,
  so the `theta:pi-1e-k` strata, which stop strictly below pi, are inside it), and for `screw` the
  pose is X times `mp.expm` of the 4x4 hat matrix where `mp.logm`-free `mp.expm` is sound;
- `local`: the chart's own `retract` of the answer returns Y (CH.4(a) holds for every Y), the
  rotation compared after sign alignment and the translation relative to its own size.
"""

from collections.abc import Callable

from mpmath import mpf

from . import charts, check, sen3, so3
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
