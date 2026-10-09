"""SO(3) ids (docs/PHASE1.md section 4.3), computed from their definitions at the working precision.

`so3_exp` is the quaternion power series, `so3_jr`/`so3_jl` the defining series in the hat matrix,
their inverses `mp.inverse`, `so3_log` the inverse of the exp series by Newton's method on it,
`so3_act` the matrix of the normalized quaternion, `so3_from_matrix` the quaternion of the
Frobenius-nearest rotation. No closed form of docs/NUMERICS.md appears; a float seed only starts an
iteration whose result is verified by the definition (residual, `check.py`).

Matrices are rows of mpf here, column-major `fmt.Mat` in records. A series stops when two terms in a
row are below `_eps()` relative to the partial sum, entry by entry: an entry of 1e-311 keeps 30
digits, as it does in the corpus.
"""

import math
from functools import lru_cache

from mpmath import mp, mpf

from .fmt import Mat
from .precision import to_f32, to_f64
from .strata import Stratum

MAX_ITER = 1000  # a series or an iteration that has not stopped by now is a bug, not slow


def _eps():
    return mpf(10) ** -(mp.dps - 10)


def norm(v):
    return mp.sqrt(sum(c * c for c in v))


def maxabs(xs):
    return max(abs(x) for x in xs)


def hat(v):
    x, y, z = v
    o = mpf(0)
    return [[o, -z, y], [z, o, -x], [-y, x, o]]


def vee_skew(M):
    """vee of the skew part (M - M^T) / 2."""
    return [(M[2][1] - M[1][2]) / 2, (M[0][2] - M[2][0]) / 2, (M[1][0] - M[0][1]) / 2]


def identity():
    return [[mpf(int(i == j)) for j in range(3)] for i in range(3)]


def mt(A):
    return [[A[j][i] for j in range(3)] for i in range(3)]


def mm(A, B):
    return [
        [A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j] for j in range(3)]
        for i in range(3)
    ]


def entries(A):
    return [x for row in A for x in row]


def qmul(a, b):
    """Hamilton product, w first."""
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return [
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    ]


def qconj(q):
    return [q[0], -q[1], -q[2], -q[3]]


def rot(q):
    """R(q) = (w^2 - |u|^2) I + 2 u u^T + 2 w [u]x (docs/NUMERICS.md section 1); for a unit q."""
    w, u = q[0], q[1:]
    d, s = w * w - sum(c * c for c in u), hat([2 * w * c for c in u])
    return [[(d if i == j else 0) + 2 * u[i] * u[j] + s[i][j] for j in range(3)] for i in range(3)]


def to_mat(A) -> Mat:
    return Mat((3, 3), [A[i][j] for j in range(3) for i in range(3)])


def from_mat(m: Mat):
    return [[mpf(m.data[3 * j + i]) for j in range(3)] for i in range(3)]


def exp_series(phi):
    """sum p^n / n!, p = (0, phi/2), by Hamilton products."""
    p, eps = [mpf(0)] + [c / 2 for c in phi], _eps()
    term = total = [mpf(1), mpf(0), mpf(0), mpf(0)]
    n = quiet = 0
    while quiet < 2:  # even terms are scalar, odd terms pure: each part against its own sum
        n += 1
        if n > MAX_ITER:
            raise ArithmeticError("exp series does not converge")
        term = [t / n for t in qmul(term, p)]
        total = [s + t for s, t in zip(total, term, strict=True)]
        part = (0,) if n % 2 == 0 else (1, 2, 3)
        small = maxabs(term[i] for i in part) <= eps * maxabs(total[i] for i in part)
        quiet = quiet + 1 if small else 0
    return total


def matrix_series(A, m: int = 1):
    """sum A^n / (n+m)!"""
    eps = _eps()
    term = total = [[x / mp.factorial(m) for x in row] for row in identity()]
    n = quiet = 0
    while quiet < 2:
        n += 1
        if n > MAX_ITER:
            raise ArithmeticError("matrix series does not converge")
        term = [[x / (n + m) for x in row] for row in mm(term, A)]
        total = [
            [s + t for s, t in zip(rs, rt, strict=True)] for rs, rt in zip(total, term, strict=True)
        ]
        small = all(
            abs(t) <= eps * abs(s) for t, s in zip(entries(term), entries(total), strict=True)
        )
        quiet = quiet + 1 if small else 0
    return total


def jr(phi):
    return matrix_series([[-x for x in row] for row in hat(phi)])


def jl(phi):
    return matrix_series(hat(phi))


def inverse(A):
    B = mp.inverse(mp.matrix(A))
    return [[B[i, j] for j in range(3)] for i in range(3)]


def normalize(q):
    n = norm(q)
    return [c / n for c in q]


def canonical(q):
    """The quaternion of NUMERICS section 3.2 after its flip: w < 0 negates; w = +0 stays."""
    q = normalize(q)
    return [-c for c in q] if q[0] < 0 else q


def log_newton(q):
    """phi with |phi| <= pi and exp_series(phi) = q (unit, w >= 0), by Newton's method on it.

    Each step is phi += J_r(phi)^-1 * 2 vec(exp(phi)^* q), as Exp(phi + d) = Exp(phi) Exp(J_r d +
    O(d^2)): the residual squares. The precision follows the accuracy reached (the step's own
    precision bounds it), and J_r^-1 needs half of it. The float seed only starts the iteration: it
    stops when the exp residual of phi, at the full precision, is below `_eps()` of |vec q|. At
    w = +0 that selects +pi u (Exp = q, not -q); a zero of vec q stays a zero of phi (`parallel`).
    """
    vmax = maxabs(q[1:])
    if vmax == 0:
        return [mpf(0)] * 3
    fv = [float(c) for c in q[1:]]
    nv = math.hypot(*fv)
    scale = 2 * math.atan2(nv, float(q[0])) / nv  # applied last: atan2 * c underflows at 1e-310
    phi = [mpf(scale * c) for c in fv]
    digits, full, eps = 14, mp.dps, _eps()  # digits: of the float seed, relative to |vec q|
    for _ in range(20):
        wp = min(full, 2 * digits + 15)
        with mp.workdps(wp):
            r = qmul(qconj(exp_series(phi)), q)
            eps_vec = [2 * c for c in r[1:]]
            if wp == full and maxabs(eps_vec) <= eps * vmax:
                return parallel(phi, q[1:])
            with mp.workdps(min(wp, digits + 15)):
                jinv = inverse(jr(phi))
            phi = [p + sum(jinv[i][k] * eps_vec[k] for k in range(3)) for i, p in enumerate(phi)]
        digits = min(2 * digits - 2, wp - 5)
    raise ArithmeticError("log does not converge")


def polar(A):
    """The rotation nearest to A in Frobenius norm: Higham's iteration X <- (X + X^-T) / 2."""
    X, eps = A, _eps()
    for _ in range(50):
        Y = [
            [(x + y) / 2 for x, y in zip(rx, ry, strict=True)]
            for rx, ry in zip(X, mt(inverse(X)), strict=True)
        ]
        change = maxabs(y - x for x, y in zip(entries(X), entries(Y), strict=True))
        X = Y
        if change <= eps * maxabs(entries(Y)):
            return X
    raise ArithmeticError("polar iteration does not converge")


def _seed(Q):
    """A float quaternion of the rotation Q (the pivot of Shepperd's method: exact for I)."""
    R = [[float(x) for x in row] for row in Q]
    t = R[0][0] + R[1][1] + R[2][2]
    d = [t, R[0][0], R[1][1], R[2][2]]
    k = d.index(max(d))
    if k == 0:
        w = math.sqrt(1 + t) / 2
        return [w, *(c / (2 * w) for c in vee_skew(R))]
    i, j, l = k - 1, k % 3, (k + 1) % 3  # (i, j, l) cyclic; q_i = s / 2 is the pivot
    s = math.sqrt(1 + 2 * R[i][i] - t)
    q = [0.0] * 4
    q[0], q[1 + i] = (R[l][j] - R[j][l]) / (2 * s), s / 2
    q[1 + j], q[1 + l] = (R[i][j] + R[j][i]) / (2 * s), (R[i][l] + R[l][i]) / (2 * s)
    return q


def quaternion_of(Q):
    """q with R(q) = Q by q <- q (1, d/2) / |.|, d = vee skew(R(q)^T Q): the error cubes (BCH).

    Canonical sign: w > 0, else the first nonzero of (x, y, z) positive. A rotation by exactly pi
    (a symmetric input, e.g. R(+0, u) rounded) has w = 0 and either sign; w below `_eps()` is that
    zero, not a value (a binary64 matrix cannot be nearer pi without being on it).
    """
    q, eps = normalize([mpf(c) for c in _seed(Q)]), _eps()
    for _ in range(20):
        d = vee_skew(mm(mt(rot(q)), Q))
        if maxabs(d) <= eps * maxabs(q[1:]):
            if abs(q[0]) <= eps:
                q = normalize([mpf(0)] + q[1:])
            first = next((c for c in q if c != 0), 1)
            return [-c for c in q] if first < 0 else q
        q = normalize(qmul(q, [mpf(1)] + [c / 2 for c in d]))
    raise ArithmeticError("quaternion iteration does not converge")


def parallel(x, v):
    """`x` with a zero wherever `v` has one: vec(Exp phi) and Log q are parallel to phi and vec q,
    and a rounding residue of 1e-122 is not the zero it should be."""
    return [c if k != 0 else mpf(0) for c, k in zip(x, v, strict=True)]


def exp(inputs: dict) -> dict:
    phi = [mpf(c) for c in inputs["phi"]]
    q = exp_series(phi)
    return {"q": q[:1] + parallel(q[1:], phi)}


def jacobian(name: str):
    """`so3_jr`, `so3_jl` and their inverses (`mp.inverse` of the series)."""
    base = jr if name.startswith("jr") else jl

    def evaluate(inputs: dict) -> dict:
        J = base([mpf(c) for c in inputs["phi"]])
        return {"J": to_mat(inverse(J) if name.endswith("inv") else J)}

    return evaluate


def gamma2(inputs: dict) -> dict:
    """`so3_gamma2` (docs/decisions/0064): Gamma_2(phi) = sum W^n / (n+2)! (NUMERICS section 7)."""
    return {"G": to_mat(matrix_series(hat([mpf(c) for c in inputs["phi"]]), 2))}


@lru_cache(maxsize=2)
def _log(q: tuple, dps: int) -> tuple:
    return tuple(log_newton(list(q)))


def log(inputs: dict) -> dict:
    """q and -q are one problem after the flip: the second of a pair is a cache hit."""
    return {"phi": _log(tuple(canonical([mpf(c) for c in inputs["q"]])), mp.dps)}


def act(inputs: dict) -> dict:
    R = rot(normalize([mpf(c) for c in inputs["q"]]))
    p = [mpf(c) for c in inputs["p"]]
    return {"Rp": [sum(R[i][k] * p[k] for k in range(3)) for i in range(3)]}


def from_matrix(inputs: dict) -> dict:
    return {"q": quaternion_of(polar(from_mat(inputs["R"])))}


# Inputs: exact binary64, one dict per record, in the order of the stratum's samples.


def phi_of(theta: float, axis, f32: bool = False) -> list[float]:
    """theta times the axis, each component rounded once: to binary32 in an `@f32` stratum."""
    rnd = to_f32 if f32 else to_f64
    return [rnd(mpf(theta) * mpf(a)) for a in axis]


def quat_of(theta: float, axis) -> list[float]:
    h = mpf(theta) / 2
    return [to_f64(mp.cos(h))] + [to_f64(mp.sin(h) * mpf(a)) for a in axis]


def quaternions(stratum: Stratum) -> list[list[float]]:
    if stratum.draw_quats is not None:
        return [list(q) for q in stratum.quaternions()]
    return [quat_of(t, a) for t, a in stratum.samples()]


def phi_inputs(stratum: Stratum) -> list[dict]:
    return [{"phi": phi_of(t, a, stratum.f32)} for t, a in stratum.samples()]


def log_inputs(stratum: Stratum) -> list[dict]:
    """Every quaternion and its negative (the half `Log` must flip); w = +0 has no negative."""
    out = []
    for q in quaternions(stratum):
        out.append({"q": q})
        if q[0] != 0:
            out.append({"q": [-c for c in q]})
    return out


def act_inputs(stratum: Stratum) -> list[dict]:
    qs = quaternions(stratum)
    return [{"q": q, "p": list(p)} for q, p in zip(qs, stratum.points(len(qs)), strict=True)]


def matrix_inputs(stratum: Stratum) -> list[dict]:
    """R(q) rounded to binary64: not exactly a rotation, and a scaled one for a non-unit q."""
    return [
        {"R": Mat((3, 3), [to_f64(x) for x in to_mat(rot([mpf(c) for c in q])).data])}
        for q in quaternions(stratum)
    ]
