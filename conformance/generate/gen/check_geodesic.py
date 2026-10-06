"""Generation-time cross-checks of the geodesic records (docs/PHASE1.md section 4.3), 0045 item 2.

`mp.logm` is the **cross-check and not the reference**: it agrees with the geometric form to
9.2e-41 for theta(Delta) <= 3.0 and returns a complex, non-principal logarithm of a rotation from
theta = 3.03 (0045, mpmath 1.3.0 and 1.4.1 identical), which is exactly the `geo:near-pi` stratum.
So the matrix route runs below MP_LOGM_LIMIT and the group identities run at **every** record: no
record is checked by the reference's own algorithm alone.

The identities are GE.2's: the two endpoints, the symmetry gamma(X_0, X_1, t) =
gamma(X_1, X_0, 1 - t), and constant body velocity in the finite form
`Log(X(t)^-1 X(t + h)) = h d`, which needs no differencing. Each needs one further evaluation of
the reference at a different argument, so a defect that is a function of `t` alone cannot pass all
three.
"""

from collections.abc import Callable

from mpmath import mp, mpf

from . import check, geodesic, sen3, so3
from .check import CrossCheckError

# 0045: the limit below which `mp.logm` of a rotation is the principal logarithm. The onset
# measured for SO(3) alone is axis-dependent in (3.02, 3.03] (`docs/maths/so3.md` SO.5).
MP_LOGM_LIMIT = mpf(3)
# The matrix route is held to `check.DIGITS` like every other cross-check here: measured over the
# whole committed corpus it agrees to **8.2e-113** (worst, `geo:generic`; 6.8e-121 at
# `geo:consecutive`), so `mp.logm` and `mp.expm` -- inverse scaling and squaring, not the
# 120-digit series the reference is -- cost thirteen digits of the hundred and leave thirteen
# spare. `geo:near-pi` is above MP_LOGM_LIMIT at every record (pi - 0.1 > 3), so that stratum is
# checked by the identities alone, which is what 0045 item 2 states.
# `Log(X(t)^-1 X(t + h)) = h d`: one fixed step, whose `h * theta` stays below pi for every theta a
# `geo:*` stratum holds, so the identity's own `Log` takes no branch.
BODY_STEP = mpf(1) / 3
# `geo:near-pi` stops strictly below pi (0045 item 3): its smallest nominal margin is 1e-12 and the
# two roundings of X_1 move theta(Delta) by about 1e-16, so this leaves a decade of room. At pi the
# two preimages give geodesics O(1) apart, so a record there would measure the sign convention.
NEAR_PI_MARGIN = mpf(10) ** -13


def _align(got, want) -> list:
    """`want` or its negative, whichever `got` is nearer: a quaternion and its negative are one
    rotation, and which one the geodesic returns follows `Log`'s branch at the far endpoint."""
    near = sum((g - w) ** 2 for g, w in zip(got, want, strict=True))
    far = sum((g + w) ** 2 for g, w in zip(got, want, strict=True))
    return list(want) if near <= far else [-w for w in want]


def _agree_pose(what: str, got: dict, want: dict) -> None:
    """Both poses to `check.DIGITS`, the quaternion after sign alignment and the translation
    relative to its own size (`||x_0||` reaches 1e4, where an absolute 1e-100 is not a
    reading)."""
    tol = mpf(10) ** -check.DIGITS
    for key in got:
        w = _align(got["q"], want["q"]) if key == "q" else want[key]
        scale = mpf(1) if key == "q" else max(mpf(1), so3.maxabs(w), so3.maxabs(got[key]))
        for g, x in zip(got[key], w, strict=True):
            if abs(g - x) > tol * scale:
                raise CrossCheckError(
                    f"{what}: {key} {mp.nstr(g, 40)} vs {mp.nstr(x, 40)}"
                )


def _endpoints(inputs: dict, out: dict, t) -> None:
    """gamma(0) = X_0 and gamma(1) = X_1, the rotations the records denote (GE.2(a))."""
    at = {mpf(0): "q0", mpf(1): "q1"}.get(t)
    if at is None:
        return
    want = {"q": so3.normalize([mpf(c) for c in inputs[at]])}
    if "x0" in inputs:
        want["x"] = [mpf(c) for c in inputs["x0" if at == "q0" else "x1"]]
    _agree_pose(f"geodesic: endpoint t = {at[1]}", out, want)


def _symmetry(evaluate: Callable, inputs: dict, out: dict, t) -> None:
    """gamma(X_0, X_1, t) = gamma(X_1, X_0, 1 - t) (GE.2(c)): the pair swapped, which negates `d`
    and so takes the other side of `Log`'s branch at the same record."""
    swapped = {"q0": inputs["q1"], "q1": inputs["q0"], "t": 1 - t}
    if "x0" in inputs:
        swapped["x0"], swapped["x1"] = inputs["x1"], inputs["x0"]
    _agree_pose("geodesic: symmetry", out, evaluate(swapped))


def _body_velocity(evaluate: Callable, inputs: dict, out: dict, t, d: dict) -> None:
    """`Log(X(t)^-1 X(t + h)) = h d` (GE.2(b)): the body displacement over a fixed parameter
    interval is the same at every `t`, which is what "constant body velocity" is without a
    difference quotient."""
    later = evaluate({**inputs, "t": t + BODY_STEP})
    between = {"q0": out["q"], "q1": later["q"]}
    if "x" in out:
        between["x0"], between["x1"] = out["x"], later["x"]
    rel = geodesic.delta(between)
    if "x" in out:
        got, want = sen3.log(1)(rel)["tau"], d["tau"]
    else:
        got, want = so3.log(rel)["phi"], d["phi"]
    check._agree(
        "geodesic: body velocity",
        got,
        [BODY_STEP * c for c in want],
        max(mpf(1), so3.maxabs(want)),
    )


def _homogeneous(q, x) -> list:
    """`[[R, x], [0, 1]]` of the rotation a record's quaternion denotes; the `3 x 3` when `x` is
    `None`."""
    r = so3.rot(so3.normalize([mpf(c) for c in q]))
    if x is None:
        return r
    m = sen3.eye(4)
    for i in range(3):
        m[i][:3] = r[i]
        m[i][3] = mpf(x[i])
    return m


def _inverse(m: list) -> list:
    """`[[R^T, -R^T x], [0, 1]]`, written out: `mp.inverse` of a `4 x 4` would be an LU where the
    structure gives the answer exactly."""
    n = len(m)
    out = [[m[j][i] for j in range(3)] + [mpf(0)] * (n - 3) for i in range(3)]
    if n == 4:
        for i in range(3):
            out[i][3] = -sum((m[k][i] * m[k][3] for k in range(3)), mpf(0))
        out.append([mpf(0)] * 3 + [mpf(1)])
    return out


def _matrix_route(inputs: dict, out: dict, t, theta) -> None:
    """X_0 expm(t logm(X_0^-1 X_1)) against the record, where `mp.logm` is sound."""
    if theta > MP_LOGM_LIMIT:
        return
    x0 = inputs.get("x0")
    m0, m1 = _homogeneous(inputs["q0"], x0), _homogeneous(inputs["q1"], inputs.get("x1"))
    log = _real(mp.logm(mp.matrix(sen3.mul(_inverse(m0), m1))), "geodesic")
    e = mp.expm(mp.matrix([[t * c for c in row] for row in log]))
    want = sen3.mul(m0, [list(row) for row in e.tolist()])
    got = {"q": so3.entries(so3.rot(out["q"]))}
    reference = {"q": so3.entries([row[:3] for row in want[:3]])}
    if x0 is not None:
        got["x"], reference["x"] = out["x"], [want[r][3] for r in range(3)]
    # `R(q)`, not `q`: the matrix route has no quaternion, so there is no sign to align.
    _agree_pose("geodesic: mp.logm", got, reference)


def _real(m, what: str) -> list:
    """`m` with every entry real, or the complex logarithm 0045 names, reported as one."""
    tol = mpf(10) ** -check.DIGITS
    out = []
    for i in range(m.rows):
        row = []
        for j in range(m.cols):
            z = m[i, j]
            if abs(mp.im(z)) > tol:
                raise CrossCheckError(
                    f"{what}: mp.logm is complex at ({i}, {j}): {mp.nstr(z, 20)}"
                )
            row.append(mp.re(z))
        out.append(row)
    return out


def _check(evaluate: Callable, inputs: dict, out: dict) -> None:
    d = geodesic.delta(inputs)
    phi = so3.log(d)["phi"]
    theta = so3.norm(phi)
    if theta > mp.pi - NEAR_PI_MARGIN:
        raise CrossCheckError(f"geodesic: theta(Delta) = {mp.nstr(theta, 40)} is at pi")
    tangent = {"phi": phi} if "x0" not in inputs else {"tau": sen3.log(1)(d)["tau"]}
    t = mpf(inputs["t"])
    _endpoints(inputs, out, t)
    _symmetry(evaluate, inputs, out, t)
    _body_velocity(evaluate, inputs, out, t, tangent)
    _matrix_route(inputs, out, t, theta)


def so3_geodesic(inputs: dict, out: dict) -> None:
    _check(geodesic.so3_evaluate, inputs, out)


def se3_geodesic(inputs: dict, out: dict) -> None:
    _check(geodesic.se3_evaluate, inputs, out)
