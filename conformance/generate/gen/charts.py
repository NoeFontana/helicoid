"""The SE(3) chart ids (docs/PHASE5.md section 1.3, 0060 decision 8), computed from their definitions.

`se3_<chart>_retract` maps a base X = (q0, x0) and a twist delta = [phi; rho] to a pose, and
`se3_<chart>_local` maps (X, Y) back to a twist, held as `phi` and `rho` so that each part is
scored on its own scale (0060 decision 8; docs/maths/charts.md section 3):

    screw       X Exp(delta)               Log(X^-1 Y)
    decoupled   (R Exp phi, t + R rho)     (Log(R^T R_Y), R^T (t_Y - t))
    world       (R Exp phi, t + rho)       (Log(R^T R_Y), t_Y - t)

`mp.logm` appears nowhere: the rotation's `Log` is `gen.so3`'s geometric one and `screw`'s
translation is `gen.sen3`'s J_l solve, as for `sen3_log` and `se3_geodesic` (0045). A record's
quaternions are rounded unit quaternions and denote the rotations they are nearest to, so both are
normalized first, as `gen.geodesic` does.

The twist is the stratum's (`gen.sen3.tau_inputs(1)`, the `theta:*` and `rho:*` families). The base
is drawn generic from the stratum's own stream: a Haar-random rotation and a translation of unit
norm. `local`'s Y is the chart's own retract of the same twist, rounded, so its answer is the
stratum's twist up to the rounding of Y.
"""

from collections.abc import Callable

from mpmath import mp, mpf

from . import sen3, so3
from .geodesic import delta
from .precision import DPS, to_f64
from .rng import SEED, stream, unit_quaternion_s3, unit_vector_s2
from .strata import Stratum

CHARTS = ("screw", "decoupled", "world")


def _frame(q) -> tuple[list, list]:
    n = so3.normalize([mpf(c) for c in q])
    return n, so3.rot(n)


def _apply(R, v) -> list:
    return [sum((R[r][c] * v[c] for c in range(3)), mpf(0)) for r in range(3)]


def _retract(chart: str, inputs: dict) -> dict:
    n0, r0 = _frame(inputs["q0"])
    tau = [mpf(c) for c in inputs["tau"]]
    x0 = [mpf(c) for c in inputs["x0"]]
    if chart == "screw":
        e = sen3.exp(1)({"tau": tau})
        q, dx = so3.qmul(n0, e["q"]), _apply(r0, e["x"])
    else:
        q = so3.qmul(n0, so3.exp({"phi": tau[:3]})["q"])
        dx = _apply(r0, tau[3:]) if chart == "decoupled" else tau[3:]
    return {"q": q, "x": [a + b for a, b in zip(x0, dx, strict=True)]}


def retract(chart: str) -> Callable[[dict], dict]:
    def evaluate(inputs: dict) -> dict:
        return _retract(chart, inputs)

    return evaluate


def local(chart: str) -> Callable[[dict], dict]:
    def evaluate(inputs: dict) -> dict:
        d = delta(inputs)  # (q0* q1, R0^T (x1 - x0)), both quaternions normalized
        if list(inputs["q0"]) == list(inputs["q1"]):
            # The rotation step rounded away (`theta:exact0`, `theta:subnormal`): q0* q1 is the
            # identity exactly, where the product at finite precision leaves a 1e-122 vector part
            # that the two-precision recheck rightly rejects.
            d["q"] = [mpf(1), mpf(0), mpf(0), mpf(0)]
        if chart == "screw":
            tau = sen3.log(1)(d)["tau"]
            return {"phi": list(tau[:3]), "rho": list(tau[3:])}
        phi = so3.log({"q": d["q"]})["phi"]
        if chart == "decoupled":
            rho = d["x"]
        else:
            rho = [mpf(b) - mpf(a) for a, b in zip(inputs["x0"], inputs["x1"], strict=True)]
        return {"phi": list(phi), "rho": list(rho)}

    return evaluate


# Inputs: exact binary64, one dict per record, in the order of the stratum's twists.


def _bases(stratum: Stratum, count: int) -> list[tuple[list[float], list[float]]]:
    rng = stream(SEED, stratum.name, "chart-base")
    with mp.workdps(DPS):
        out = []
        for _ in range(count):
            q = [to_f64(c) for c in unit_quaternion_s3(rng)]
            x = [to_f64(mpf(c)) for c in unit_vector_s2(rng)]
            out.append((q, x))
        return out


def retract_inputs(stratum: Stratum) -> list[dict]:
    taus = sen3.tau_inputs(1)(stratum)
    return [
        {"q0": q0, "x0": x0, "tau": t["tau"]}
        for (q0, x0), t in zip(_bases(stratum, len(taus)), taus, strict=True)
    ]


def local_inputs(chart: str) -> Callable[[Stratum], list[dict]]:
    def inputs(stratum: Stratum) -> list[dict]:
        out = []
        with mp.workdps(DPS):
            for rec in retract_inputs(stratum):
                y = _retract(chart, rec)
                out.append(
                    {
                        "q0": rec["q0"],
                        "x0": rec["x0"],
                        "q1": [to_f64(c) for c in y["q"]],
                        "x1": [to_f64(c) for c in y["x"]],
                    }
                )
        return out

    return inputs
