"""The inverse-action ids `so3_inv_act` and `se3_inv_act` (0103), computed from their definitions.

The value is R^T (p - t) of the normalized quaternion. The Jacobians are the derivatives at 0 of
the definitions, (X Exp(tau))^-1 p = Exp(-tau) X^-1 p (right) and (Exp(tau) X)^-1 p =
X^-1 Exp(-tau) p (left): d/ds Exp(s e_j) at 0 is the generator e_j^, so column j is
-(e_j^ [x; 1]) on the right and X^-1 (-(e_j^ [p; 1])) on the left, in homogeneous coordinates. No
row of docs/NUMERICS.md section 2.4 appears; `check_inv_act` checks those rows against these.

`se3_inv_act`'s strata draw the two cases where the composition X.inverse() * p cancels: `tx:<r>`
puts the body close to the point and far from the origin, |t| = r |x| with |x| uniform in
[1, 10] and p = R x + t rounded; `tp:<r>` puts the point near the parent's origin, |t| = r |p|
with |p| uniform in [1, 10]. `generic` draws t and p uniform in [-10, 10]^3. q is uniform on S3,
rounded. Each has an `@f32` twin, the same draws rounded once to binary32; rounding p there moves
x by up to about 2 sqrt(3) u r |x|, so `tx:1e6@f32` holds its ratio to about 4% (measured 3.3%),
the stratum being the regime and not the exact ratio. `so3_inv_act` reads `so3_act`'s records.
"""

from mpmath import mp, mpf

from . import sen3, so3
from .fmt import Mat
from .rng import unit_quaternion_s3
from .strata import SEN3_SAMPLES, Draw, drawn_strata

RATIOS = ("1", "1e2", "1e4", "1e6")


def _point(v) -> list:
    return [mpf(c) for c in v]


def _cols(cols: list, rows: int) -> Mat:
    """Columns of mpf as a column-major `Mat`."""
    return Mat((rows, len(cols)), [c for col in cols for c in col])


def _rotate(R, v) -> list:
    return [sum((R[i][k] * v[k] for k in range(3)), mpf(0)) for i in range(3)]


def _evaluate(q, t, p, dof: int) -> tuple:
    """(R^T (p - t), the right and left Jacobians' columns, R^T), from the generators of the
    `dof`-dimensional group: 3 for SO(3) (no translation), 6 for SE(3)."""
    R = so3.rot(so3.normalize(q))
    Rt = so3.mt(R)
    x = _rotate(Rt, [a - b for a, b in zip(p, t, strict=True)])
    right, left = [], []
    for j in range(dof):
        G = sen3.hat(sen3.basis(6, j), 1)  # e_j^ in se(3), rotation first
        right.append([-sum((G[i][k] * c for k, c in enumerate(x + [1])), mpf(0)) for i in range(3)])
        gp = [-sum((G[i][k] * c for k, c in enumerate(p + [1])), mpf(0)) for i in range(3)]
        left.append(_rotate(Rt, gp))
    return x, right, left, Rt


def so3_evaluate(inputs: dict) -> dict:
    zero = [mpf(0)] * 3
    x, right, left, Rt = _evaluate(_point(inputs["q"]), zero, _point(inputs["p"]), 3)
    return {
        "Rinv_p": x,
        "Jr": _cols(right, 3),
        "Jl": _cols(left, 3),
        "Jp": so3.to_mat(Rt),
    }


def se3_evaluate(inputs: dict) -> dict:
    q, t, p = _point(inputs["q"]), _point(inputs["x"]), _point(inputs["p"])
    x, right, left, Rt = _evaluate(q, t, p, 6)
    return {
        "Xinv_p": x,
        "Jr": _cols(right, 3),
        "Jl": _cols(left, 3),
        "Jp": so3.to_mat(Rt),
    }


def _direction(rng) -> list:
    """Uniform on S2, unrounded."""
    z = 2 * mpf(rng.uniform()) - 1
    a = 2 * mp.pi * mpf(rng.uniform())
    r = mp.sqrt(1 - z * z)
    return [r * mp.cos(a), r * mp.sin(a), z]


def _make(kind: str, ratio: str | None):
    def make(draw: Draw, count: int) -> list[dict]:
        rq, rn, rd, rt = (draw.stream(s) for s in ("q", "norm", "dir", "t"))
        out = []
        for _ in range(count):
            q = [draw.round(c) for c in unit_quaternion_s3(rq)]
            if kind == "generic":
                t = [draw.round(20 * mpf(rt.uniform()) - 10) for _ in range(3)]
                p = [draw.round(20 * mpf(rd.uniform()) - 10) for _ in range(3)]
            else:
                n, scale = 1 + 9 * mpf(rn.uniform()), mpf(ratio)
                t = [draw.round(scale * n * c) for c in _direction(rt)]
                v = [n * c for c in _direction(rd)]
                if kind == "tx":  # v is x, the body-frame point
                    R = so3.rot(so3.normalize(_point(q)))
                    p = [draw.round(a + b) for a, b in zip(_rotate(R, v), t, strict=True)]
                else:  # v is p itself
                    p = [draw.round(c) for c in v]
            out.append({"q": q, "x": t, "p": p})
        return out

    return make


def _catalogue() -> list:
    out = [(f"{kind}:{r}", _make(kind, r), SEN3_SAMPLES) for kind in ("tx", "tp") for r in RATIOS]
    return out + [("generic", _make("generic", None), SEN3_SAMPLES)]


SE3_STRATA = drawn_strata(_catalogue(), [(n, m, c, True) for n, m, c in _catalogue()])
