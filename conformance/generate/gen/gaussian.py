"""`gaussian_mahalanobis_se3` and `gaussian_mahalanobis_se23` (docs/PHASE5.md section 5, 0065).

A record is a mean X = (q0, x0), a point Y = (q1, x1), a `side` (0 right, 1 left) and a symmetric
`Sigma` of the tangent's size. The reference is

    delta = Log(X^-1 Y)  (right)   or   Log(Y X^-1)  (left),
    d2    = delta^T Sigma^-1 delta,   `valid` = Sigma positive definite,

decided as `gen.chol` decides it: `valid` exactly, by an LDL^T of the stored Sigma in rationals,
and `d2` by `mp.lu_solve` of the stored Sigma, 0 where `valid` is 0 (the field is then not scored,
docs/maths/gamma-gaussian.md GG.11(c)). `Log` is the geometric one of `gen.sen3.log` (the SO(3)
quaternion route and a J_l solve), never `mp.logm` (0045); both quaternions are normalized first,
as `gen.charts` does.

Y is drawn from the law the record's Gaussian states and rounded: Y = X Exp(xi) or Exp(xi) X with
xi = L z, Sigma = L L^T the stored matrix's factor and z standard normal, so d2 is a chi-square
draw (GG.10(b)). `gauss:near-cut` rescales the rotation of xi to theta in [2.5, 3], where Log is
still on its branch (GG.10(d)); `gauss:indefinite` draws xi = D0 z, Sigma having no factor.

Sigma = D0 (Q diag(s) Q^T) D0 with D0 = diag(sigma_phi I3, sigma_rho I_3n), each sigma
log-uniform in [1e-3, 1] per record, Q a product of Householder reflections (`gen.chol`'s), and s:

    gauss:right, gauss:left   10^U[-1, 0]
    gauss:corr-1e-k           log-spaced from 1 to 10^-k (`chol:cond-1e-k`'s)
    gauss:lever               right's, with sigma_phi = 1e-2 and sigma_rho = 1e-1 fixed, and the
                              *left* form Ad_X Sigma Ad_X^T of it stored at |t| = 1e3 (GG.13)
    gauss:near-cut            right's
    gauss:indefinite          right's with one entry -max(s)/100 (`chol:indefinite`'s)

The base's translation columns are unit vectors times 10^U[-1, 1] (1e3 at `gauss:lever`).
"""

from collections.abc import Callable

from mpmath import mp, mpf

from . import sen3, so3
from .chol import INDEFINITE, _householder, _log10_uniform
from .dense import exact, from_mat, gaussian, ldlt_pivots, round_symmetric, unit_vector
from .rng import unit_quaternion_s3
from .strata import Draw, drawn_strata

N_GAUSS = 32
COND_K = {False: (4, 8), True: (2, 4)}
SIGMA = (-3, 0)  # log10 range of sigma_phi and sigma_rho
LEVER = 1000
NEAR_CUT = (mpf("2.5"), mpf(3))

LEFT = {"gauss:right": False, "gauss:left": True, "gauss:lever": True}


def _dof(n: int) -> int:
    return 3 + 3 * n


def _spectrum(name: str, rd, m: int) -> list:
    if name.startswith("gauss:corr-1e-"):
        k = int(name.rsplit("-", 1)[1])
        return [mpf(10) ** (-mpf(k) * i / (m - 1)) for i in range(m)]
    s = [_log10_uniform(rd, -1, 0) for _ in range(m)]
    if name == "gauss:indefinite":
        s[rd.next_u64() % m] = -max(s) / INDEFINITE
    return s


def _pose(draw: Draw, rp, n: int, lever) -> tuple[list, list]:
    q = [draw.round(c) for c in unit_quaternion_s3(rp)]
    x = []
    for _ in range(n):
        scale = lever if lever is not None else _log10_uniform(rp, -1, 1)
        x.extend(draw.round(scale * c) for c in unit_vector(rp, 3))
    return q, x


def _compose(a: dict, b: dict, n: int) -> dict:
    """(q_a q_b, R_a x_b + x_a), column by column."""
    Ra = so3.rot(so3.normalize(a["q"]))
    xs = [
        sum((Ra[r][c] * xb[c] for c in range(3)), mpf(0)) + xa[r]
        for xa, xb in zip(sen3.chunks(a["x"], n), sen3.chunks(b["x"], n), strict=True)
        for r in range(3)
    ]
    return {"q": so3.qmul(a["q"], b["q"]), "x": xs}


def _records(name: str, n: int):
    m = _dof(n)

    def make(draw: Draw, count: int) -> list[dict]:
        rq, rd, rp, rx = (draw.stream(p) for p in ("householder", "spectrum", "pose", "xi"))
        left = LEFT.get(name, False)
        out = []
        for _ in range(count):
            lever = LEVER if name == "gauss:lever" else None
            q0, x0 = _pose(draw, rp, n, lever)
            if name == "gauss:lever":
                sphi, srho = mpf("1e-2"), mpf("1e-1")
            else:
                sphi, srho = (_log10_uniform(rd, *SIGMA) for _ in range(2))
            d0 = [sphi] * 3 + [srho] * (3 * n)
            Q = _householder(rq, m)
            A = mp.diag(d0) * Q * mp.diag(_spectrum(name, rd, m)) * Q.T * mp.diag(d0)
            if name == "gauss:lever":
                R = so3.rot(so3.normalize([mpf(c) for c in q0]))
                Ad = mp.matrix(sen3.adjoint(R, [mpf(c) for c in x0], n))
                A = Ad * A * Ad.T
            sigma = round_symmetric(A, draw.round)
            z = mp.matrix([gaussian(rx) for _ in range(m)])
            if ldlt_pivots(exact(sigma)) is None:
                xi = mp.diag(d0) * z
            else:
                xi = mp.cholesky(from_mat(sigma), tol=0) * z
            xi = [xi[i] for i in range(m)]
            if name == "gauss:near-cut":
                theta = NEAR_CUT[0] + (NEAR_CUT[1] - NEAR_CUT[0]) * mpf(rx.uniform())
                r = so3.norm(xi[:3])
                xi[:3] = [c * theta / r for c in xi[:3]]
            e = sen3.exp(n)({"tau": xi})
            x = {"q": [mpf(c) for c in q0], "x": [mpf(c) for c in x0]}
            y = _compose(e, x, n) if left else _compose(x, e, n)
            out.append(
                {
                    "q0": q0,
                    "x0": x0,
                    "q1": [draw.round(c) for c in y["q"]],
                    "x1": [draw.round(c) for c in y["x"]],
                    "side": float(left),
                    "Sigma": sigma,
                }
            )
        return out

    return make


def strata(n: int):
    base = ["gauss:right", "gauss:left"]
    tail = ["gauss:near-cut", "gauss:indefinite"]
    wide = [*base, *(f"gauss:corr-1e-{k}" for k in COND_K[False]), "gauss:lever", *tail]
    narrow = [*base, *(f"gauss:corr-1e-{k}" for k in COND_K[True]), *tail]
    return drawn_strata(
        [(s, _records(s, n), N_GAUSS) for s in wide],
        [(s, _records(s, n), N_GAUSS, s in wide) for s in narrow],
    )


# --- the reference ---------------------------------------------------------------------------


def residual(inputs: dict, n: int, left: bool) -> list:
    """Y minus X on `left`'s side, by the geometric Log."""
    n0 = so3.normalize([mpf(c) for c in inputs["q0"]])
    n1 = so3.normalize([mpf(c) for c in inputs["q1"]])
    R0, R1 = so3.rot(n0), so3.rot(n1)
    x0 = sen3.chunks([mpf(c) for c in inputs["x0"]], n)
    x1 = sen3.chunks([mpf(c) for c in inputs["x1"]], n)
    pairs = list(zip(x0, x1, strict=True))
    if left:  # Y X^-1 = (R1 R0^T, x1 - R1 R0^T x0)
        q = so3.qmul(n1, so3.qconj(n0))
        M = so3.mm(R1, so3.mt(R0))
        x = [
            b[r] - sum((M[r][c] * a[c] for c in range(3)), mpf(0)) for a, b in pairs for r in range(3)
        ]
    else:  # X^-1 Y = (R0^T R1, R0^T (x1 - x0))
        q = so3.qmul(so3.qconj(n0), n1)
        x = [
            sum((R0[c][r] * (b[c] - a[c]) for c in range(3)), mpf(0))
            for a, b in pairs
            for r in range(3)
        ]
    if list(inputs["q0"]) == list(inputs["q1"]):
        # The rotation step rounded away: the product leaves a 1e-122 vector part that the
        # two-precision recheck rightly rejects, as in `gen.charts.local`.
        q = [mpf(1), mpf(0), mpf(0), mpf(0)]
    return list(sen3.log(n)({"q": q, "x": x})["tau"])


def mahalanobis(n: int) -> Callable[[dict], dict]:
    def evaluate(inputs: dict) -> dict:
        side = inputs["side"]
        valid = ldlt_pivots(exact(inputs["Sigma"])) is not None
        if not valid:
            return {"valid": [mpf(0)], "d2": [mpf(0)]}
        delta = residual(inputs, n, bool(side))
        x = mp.lu_solve(from_mat(inputs["Sigma"]), mp.matrix(delta))
        return {"valid": [mpf(1)], "d2": [sum((d * x[i] for i, d in enumerate(delta)), mpf(0))]}

    return evaluate
