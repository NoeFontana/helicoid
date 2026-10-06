"""The geodesic ids (docs/PHASE1.md section 4.3), computed from their definition.

`gamma(X_0, X_1, t) = X_0 Exp(t Log(X_0^-1 X_1))` (docs/NUMERICS.md section 10). `mp.logm` appears
**nowhere** in the forward path: it returns a complex, non-principal logarithm of a rotation from
theta = 3.03, which is exactly the `geo:near-pi` stratum, so the reference is the geometric `Log`
of `gen.so3` (the quaternion `atan2` route by Newton on the exp series) with `gen.sen3`'s J_l solve
for the translation block, then `gen.sen3`'s matrix exponential of t*d. `gen/check_geodesic.py`
holds `mp.logm` where 0045 measured it sound, below theta = 3.

A record's q_0 and q_1 are rounded unit quaternions and denote the rotations they are nearest to,
so both are normalized first, as `so3.act` and `check_sen3.log` do. At t = 0 that makes the answer
X_0 with the quaternion normalized, which is where a program carrying the stored q_0 reads about
half an ulp; it is the geodesic between the two *rotations* the records denote, and the alternative
-- the stored, scaled quaternions -- is not a geodesic on the group at all.

`X_1 = X_0 Exp(d)` is **formed** here and rounded, never drawn: the relative motion of a record is
then the `d` its stratum names, up to two roundings, which is what makes `geo:consecutive` the
regime it claims to be.
"""

from mpmath import mp, mpf

from . import sen3, so3
from .precision import DPS, to_f64
from .rng import SEED, stream, unit_quaternion_s3, unit_vector_s2
from .strata import GEO_SAMPLES, Stratum

# The interpolation parameters every `geo:*` stratum is crossed with: the six of PHASE4.md section
# 4, which include both endpoints (0045 decision 3 keeps them: the cheapest rows in the corpus and
# the only ones that catch an endpoint error), then GEO_T_UNIFORM drawn uniform on (0, 1) from the
# stratum's own stream. A record's `t` is exact binary64, like every other input.
GEO_T_FIXED = (
    lambda: mpf(0),
    lambda: mpf(10) ** -9,
    lambda: mpf(1) / 4,
    lambda: mpf(1) / 2,
    lambda: 1 - mpf(10) ** -9,
    lambda: mpf(1),
)
GEO_T_UNIFORM = 4

# `geo:consecutive` is tf_tree's kilohertz edge: the relative rotation and the relative translation
# are the same size, so the whole of `d` is small. The other two strata put the relative
# translation at unit scale, where the rotation alone carries the conditioning.
CONSECUTIVE = "geo:consecutive"

# ||x_0|| of the six samples of every `geo:*` stratum: PHASE4.md section 4 asks `geo:consecutive`
# for ||t_0|| up to 1e4, and giving the other two the same span costs nothing and is what
# PHASE4.md section 5.2's left-invariance bound is measured against.
GEO_POSE_EXP = (0, 0, 2, 2, 4, 4)


def _samples(stratum: Stratum) -> list:
    """(theta of `d`, axis of `d`) per sample, as `gen.sen3._samples` pairs them."""
    thetas = stratum.thetas()
    axes = stratum.axes(GEO_SAMPLES)
    return [(thetas[i % len(thetas)], axes[i]) for i in range(GEO_SAMPLES)]


def _ts(stratum: Stratum) -> list[float]:
    rng = stream(SEED, stratum.name, "t")
    with mp.workdps(DPS):
        fixed = [to_f64(value()) for value in GEO_T_FIXED]
        return fixed + [to_f64(mpf(rng.uniform())) for _ in range(GEO_T_UNIFORM)]


def _pairs(stratum: Stratum) -> list[tuple]:
    """One `(q_0, x_0, q_1, x_1)` per sample, every value exact binary64.

    One stream per role, so `so3_geodesic` and `se3_geodesic` see the same poses: the rotation of
    `Exp[phi; rho]` does not depend on `rho`, so `q_1` is the same whether or not a record carries
    a translation, and the two ids are one problem measured twice.
    """
    streams = {k: stream(SEED, stratum.name, k) for k in ("q0", "x0", "rho")}
    out = []
    with mp.workdps(DPS):
        for i, (theta, axis) in enumerate(_samples(stratum)):
            q0 = [to_f64(c) for c in unit_quaternion_s3(streams["q0"])]
            scale = mpf(10) ** GEO_POSE_EXP[i]
            x0 = [to_f64(scale * mpf(c)) for c in unit_vector_s2(streams["x0"])]
            relative = mpf(theta) if stratum.name == CONSECUTIVE else mpf(1)
            phi = [mpf(theta) * mpf(a) for a in axis]
            rho = [relative * mpf(c) for c in unit_vector_s2(streams["rho"])]
            n0, r0 = _frame([mpf(c) for c in q0])
            e = sen3.exp(1)({"tau": phi + rho})
            q1 = [to_f64(c) for c in so3.qmul(n0, e["q"])]
            x1 = [to_f64(mpf(x0[r]) + _dot(r0[r], e["x"])) for r in range(3)]
            out.append((q0, x0, q1, x1))
    return out


def _dot(row, v) -> mpf:
    return sum((a * b for a, b in zip(row, v, strict=True)), mpf(0))


def _frame(q) -> tuple[list, list]:
    """(q/||q||, R(q/||q||)): the rotation a record's quaternion denotes, and its matrix."""
    n = so3.normalize(q)
    return n, so3.rot(n)


def delta(inputs: dict) -> dict:
    """`X_0^-1 X_1 = (q_0^* q_1, R_0^T (x_1 - x_0))`, the relative transform a record holds.

    Shared with `check_geodesic`, which needs theta(Delta) to decide whether `mp.logm` applies.
    `x_1`/`x_0` absent (`so3_geodesic`) gives the rotation alone.
    """
    n0, r0 = _frame([mpf(c) for c in inputs["q0"]])
    n1 = so3.normalize([mpf(c) for c in inputs["q1"]])
    out = {"q": so3.qmul(so3.qconj(n0), n1)}
    if "x0" in inputs:
        dx = [mpf(b) - mpf(a) for a, b in zip(inputs["x0"], inputs["x1"], strict=True)]
        out["x"] = [_dot([r0[k][i] for k in range(3)], dx) for i in range(3)]
    return out


def so3_evaluate(inputs: dict) -> dict:
    phi = so3.log(delta(inputs))["phi"]
    t = mpf(inputs["t"])
    n0 = so3.normalize([mpf(c) for c in inputs["q0"]])
    return {"q": so3.qmul(n0, so3.exp_series([t * c for c in phi]))}


def se3_evaluate(inputs: dict) -> dict:
    tau = sen3.log(1)(delta(inputs))["tau"]
    t = mpf(inputs["t"])
    n0, r0 = _frame([mpf(c) for c in inputs["q0"]])
    e = sen3.exp(1)({"tau": [t * c for c in tau]})
    return {
        "q": so3.qmul(n0, e["q"]),
        "x": [mpf(inputs["x0"][r]) + _dot(r0[r], e["x"]) for r in range(3)],
    }


# Inputs: exact binary64, one dict per record, each pose pair crossed with every `t` in order.


def so3_inputs(stratum: Stratum) -> list[dict]:
    return [
        {"q0": q0, "q1": q1, "t": t}
        for q0, _, q1, _ in _pairs(stratum)
        for t in _ts(stratum)
    ]


def se3_inputs(stratum: Stratum) -> list[dict]:
    return [
        {"q0": q0, "x0": x0, "q1": q1, "x1": x1, "t": t}
        for q0, x0, q1, x1 in _pairs(stratum)
        for t in _ts(stratum)
    ]
