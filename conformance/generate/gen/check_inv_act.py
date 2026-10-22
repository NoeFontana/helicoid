"""Generation-time cross-checks of the inverse-action records (0103), to 100 digits.

The value is pinned by R x + t = p, with the rotation as the sandwich q x q^* of the normalized
quaternion rather than `so3.rot`'s matrix; R^T is the sandwich of q^* on the basis. The Jacobians,
which `gen.inv_act` forms from the generators, are checked against docs/NUMERICS.md section 2.4's
rows: right [[x]x, -I], left [R^T [p]x, -R^T] (SO(3): their first blocks). Tolerances are
1e-100 of the largest of 1 and the entries compared.
"""

from mpmath import mpf

from . import so3
from .check import _agree


def _sandwich(q, v) -> list:
    return so3.qmul(so3.qmul(q, [mpf(0)] + list(v)), so3.qconj(q))[1:]


def _cols(m) -> list:
    rows, cols = m.shape
    return [m.data[c * rows : (c + 1) * rows] for c in range(cols)]


def _close(what: str, got, want) -> None:
    got, want = list(got), list(want)
    _agree(what, got, want, max(1, so3.maxabs(got), so3.maxabs(want)))


def _check(name: str, inputs: dict, out: dict, value: str, dof: int) -> None:
    q = so3.normalize([mpf(c) for c in inputs["q"]])
    t = [mpf(c) for c in inputs.get("x", [0, 0, 0])]
    p = [mpf(c) for c in inputs["p"]]
    x = out[value]
    _close(f"{name}: R x + t = p", [a + b for a, b in zip(_sandwich(q, x), t, strict=True)], p)
    qc = so3.qconj(q)
    rt_cols = [_sandwich(qc, e) for e in so3.identity()]  # R^T's columns
    rt_rows = [[rt_cols[c][r] for c in range(3)] for r in range(3)]
    _close(f"{name}: Jp = R^T", [c for col in _cols(out["Jp"]) for c in col], sum(rt_cols, []))
    hx, hp = so3.hat(x), so3.hat(p)
    right = [[hx[r][c] for r in range(3)] for c in range(3)]
    left = [[sum(rt_rows[r][k] * hp[k][c] for k in range(3)) for r in range(3)] for c in range(3)]
    if dof == 6:
        right += [[-mpf(int(r == c)) for r in range(3)] for c in range(3)]
        left += [[-rt_cols[c][r] for r in range(3)] for c in range(3)]
    for key, want in (("Jr", right), ("Jl", left)):
        _close(f"{name}: {key}", [c for col in _cols(out[key]) for c in col], sum(want, []))


def so3_inv_act(inputs: dict, out: dict) -> None:
    _check("so3_inv_act", inputs, out, "Rinv_p", 3)


def se3_inv_act(inputs: dict, out: dict) -> None:
    _check("se3_inv_act", inputs, out, "Xinv_p", 6)
