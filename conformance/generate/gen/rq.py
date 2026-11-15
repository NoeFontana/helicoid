"""`rq3` (docs/decisions/0117 decision 9): M = U R of the stored 3x3 `M`, U upper triangular,
R in SO(3), u22 >= 0, u33 >= 0, u11 of the sign of det M (maths/rq3.md RQ.1).

A planted M is drawn at the working precision and rounded once entrywise to the stratum's
precision, as `svd3`'s are: the reference is the stored matrix's factorization. `lq3` has no id:
it is `rq3` of the reversed matrix, reversed, bit for bit (docs/NUMERICS.md §14).
"""

from mpmath import mp, mpf

from .dense import det3, exact, from_mat, gaussian, rank, round_matrix, to_mat
from .fmt import Mat
from .strata import Draw, deep, drawn_strata
from .svd import _ZERO_SIGNS, gaussian_matrix, rotation, uniform

N_RQ = deep(16)
N_WIDE = deep(32)  # `rq:random`
N_ZERO = 4  # the four sign patterns of `svd:zero`, at every tier
NEAR_K = {False: (4, 8, 12), True: (2, 4, 6)}  # near-triangular, near-singular
PARALLEL_K = {False: (2, 4, 8, 12, 16), True: (2, 4, 6, 7)}  # the last at rounding: Kahan's test
DIAG_EXP = {False: 30, True: 6}  # rq:diagonal, rq:triangular: 10^U(-e, e)
GRADES = {False: (0, 5, 10), True: (0, 2, 4)}  # rows or columns times 10^-g
SPREAD = {False: 150, True: 18}  # rq:row-spread: rows times 10^(s, 0, -s)
SCALE = {False: 70, True: 12}  # rq:scale-up/down: 10^(+-SCALE)


def _random(sign: int = 0):
    """`rq:random`, Gaussian entries, from its stream, times 10^(sign s), s the precision's
    `SCALE`."""

    def make(draw: Draw, n: int) -> list[dict]:
        rng = draw.stream("entries", of="rq:random")
        factor = mpf(10) ** (sign * SCALE[draw.f32])
        return [{"M": round_matrix(factor * gaussian_matrix(rng), draw.round)} for _ in range(n)]

    return make


def _diagonal(upper: bool):
    """A positive diagonal, log-uniform in 10^(-e, e); `upper`: Gaussian entries above it."""

    def make(draw: Draw, n: int) -> list[dict]:
        rd, re = draw.stream("diagonal"), draw.stream("entries")
        e = DIAG_EXP[draw.f32]
        out = []
        for _ in range(n):
            d = [mpf(10) ** uniform(rd, -e, e) for _ in range(3)]
            M = mp.diag(d)
            if upper:
                for i in range(3):
                    for j in range(i + 1, 3):
                        M[i, j] = gaussian(re)
            out.append({"M": round_matrix(M, draw.round)})
        return out

    return make


def _near_triangular(k: int):
    """A triangle with diagonal |g| + 1/2 and Gaussian entries above, plus Gaussian entries of
    size 10^-k below."""

    def make(draw: Draw, n: int) -> list[dict]:
        rng = draw.stream("entries")
        out = []
        for _ in range(n):
            G = gaussian_matrix(rng)
            for i in range(3):
                G[i, i] = abs(G[i, i]) + mpf(1) / 2
                for j in range(i):
                    G[i, j] *= mpf(10) ** -k
            out.append({"M": round_matrix(G, draw.round)})
        return out

    return make


def _near_identity(draw: Draw, n: int) -> list[dict]:
    """I + 10^-3 G: a calibration's scale and misalignment."""
    rng = draw.stream("entries")
    return [
        {"M": round_matrix(mp.eye(3) + mpf(10) ** -3 * gaussian_matrix(rng), draw.round)}
        for _ in range(n)
    ]


def _parallel(k: int):
    """m2 at angle 10^-k from m3, in a random plane, of length in [1/2, 3/2)."""

    def make(draw: Draw, n: int) -> list[dict]:
        rng = draw.stream("entries")
        out = []
        for _ in range(n):
            G = gaussian_matrix(rng)
            a = [G[2, j] for j in range(3)]
            t = [gaussian(rng) for _ in range(3)]
            aa = sum(x * x for x in a)
            ta = sum(x * y for x, y in zip(t, a, strict=True)) / aa
            p = [x - ta * y for x, y in zip(t, a, strict=True)]
            na, np_ = mp.sqrt(aa), mp.sqrt(sum(x * x for x in p))
            e, s = mpf(10) ** -k, uniform(rng, mpf(1) / 2, mpf(3) / 2)
            for j in range(3):
                G[1, j] = s * (mp.cos(e) * a[j] / na + mp.sin(e) * p[j] / np_)
            out.append({"M": round_matrix(G, draw.round)})
        return out

    return make


def _graded(kind: str):
    """Gaussian entries, rows (`row-graded`) or columns (`col-graded`) times 10^-g, or rows times
    10^(s, 0, -s) (`row-spread`)."""

    def make(draw: Draw, n: int) -> list[dict]:
        rng = draw.stream("entries")
        if kind == "row-spread":
            s = SPREAD[draw.f32]
            scale = [mpf(10) ** s, mpf(1), mpf(10) ** -s]
        else:
            scale = [mpf(10) ** -g for g in GRADES[draw.f32]]
        D = mp.diag(scale)
        out = []
        for _ in range(n):
            G = gaussian_matrix(rng)
            out.append({"M": round_matrix(G * D if kind == "col-graded" else D * G, draw.round)})
        return out

    return make


def _near_singular(k: int):
    """m1 = a m2 + b m3 + 10^-k t, Gaussian a, b, t: u11 small against |m1|."""

    def make(draw: Draw, n: int) -> list[dict]:
        rng = draw.stream("entries")
        out = []
        for _ in range(n):
            G = gaussian_matrix(rng)
            a, b = gaussian(rng), gaussian(rng)
            for j in range(3):
                G[0, j] = a * G[1, j] + b * G[2, j] + mpf(10) ** -k * gaussian(rng)
            out.append({"M": round_matrix(G, draw.round)})
        return out

    return make


def _negative_det(draw: Draw, n: int) -> list[dict]:
    """Gaussian entries, m1 negated where the stored det is positive: det M < 0 exactly."""
    rng = draw.stream("entries")
    out = []
    for _ in range(n):
        m = round_matrix(gaussian_matrix(rng), draw.round)
        if det3(exact(m)) > 0:
            m = Mat((3, 3), [-x if i % 3 == 0 else x for i, x in enumerate(m.data)])
        out.append({"M": m})
    return out


def _rank_bottom(draw: Draw, n: int) -> list[dict]:
    """m2 = m3 / 2 exactly: rows 2 and 3 dependent, the factorization not unique."""
    rng = draw.stream("entries")
    out = []
    for _ in range(n):
        m = round_matrix(gaussian_matrix(rng), draw.round)
        data = list(m.data)
        for j in range(3):
            data[3 * j + 1] = data[3 * j + 2] / 2
        out.append({"M": Mat((3, 3), data)})
    return out


def _zero(draw: Draw, n: int) -> list[dict]:
    pattern = [_ZERO_SIGNS[i % len(_ZERO_SIGNS)] for i in range(n)]
    return [{"M": Mat((3, 3), [0.0 if s > 0 else -0.0 for s in signs])} for signs in pattern]


def _haar_rotation(draw: Draw, n: int) -> list[dict]:
    rng = draw.stream("quat")
    return [{"M": round_matrix(rotation(rng), draw.round)} for _ in range(n)]


def _named(f32: bool) -> list[tuple[str, object, int, bool]]:
    """(name, make, count, shared) in the catalogue's order; `shared`: the binary64 namesake's
    range, so the `@f32` stratum reads its stream."""
    near = [
        (f"rq:near-triangular-1e-{k}", _near_triangular(k), N_RQ, k in NEAR_K[False])
        for k in NEAR_K[f32]
    ]
    parallel = [
        (f"rq:parallel-1e-{k}", _parallel(k), N_RQ, k in PARALLEL_K[False]) for k in PARALLEL_K[f32]
    ]
    singular = [
        (f"rq:near-singular-1e-{k}", _near_singular(k), N_RQ, k in NEAR_K[False])
        for k in NEAR_K[f32]
    ]
    return [
        ("rq:random", _random(), N_WIDE, True),
        ("rq:diagonal", _diagonal(False), N_RQ, False),
        ("rq:triangular", _diagonal(True), N_RQ, False),
        *near,
        ("rq:near-identity-1e-3", _near_identity, N_RQ, True),
        *parallel,
        ("rq:row-graded", _graded("row-graded"), N_RQ, False),
        ("rq:row-spread", _graded("row-spread"), N_RQ, False),
        ("rq:col-graded", _graded("col-graded"), N_RQ, False),
        *singular,
        ("rq:negative-det", _negative_det, N_RQ, True),
        ("rq:rank-bottom", _rank_bottom, N_RQ, True),
        ("rq:zero", _zero, N_ZERO, True),
        ("rq:rotation", _haar_rotation, N_RQ, True),
        ("rq:scale-up", _random(1), N_RQ, False),
        ("rq:scale-down", _random(-1), N_RQ, False),
    ]


def _strata():
    return drawn_strata(
        [(name, make, count) for name, make, count, _ in _named(False)],
        _named(True),
    )


STRATA = _strata()


# --- the reference ---------------------------------------------------------------------------


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b, strict=True))


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _unit(x):
    r = mp.sqrt(_dot(x, x))
    return [c / r for c in x]


def _axis(k: int) -> list:
    return [mpf(int(i == k)) for i in range(3)]


def completion(r3) -> tuple[list, list]:
    """docs/NUMERICS.md §22's rule: Gram-Schmidt against the unit r3 of the axis e_k of its
    smallest |component|, the last of equals, x = e_k - r3_k r3, r2 = x / |x|; and r1 = r2 x r3,
    formed as (e_k x r3) / |x|, which it equals exactly: its k-th component is 0 at every
    precision, not the rounding of a cancellation."""
    k = min(range(3), key=lambda i: (abs(r3[i]), -i))
    e = _axis(k)
    x = [p - r3[k] * q for p, q in zip(e, r3, strict=True)]
    n = mp.sqrt(_dot(x, x))
    return [c / n for c in x], [c / n for c in _cross(e, r3)]


def tangent(U, R, E):
    """maths/rq3.md RQ.6 for the direction E: (dU, dR), from rows 2 and 3 of U^-1 only, so u11 may
    be 0. X = U^-1 E R^T, Omega = low(X) - low(X)^T, dR = Omega R, dU = E R^T - U Omega."""
    ERt = E * R.T
    inv2 = [mpf(0), 1 / U[1, 1], -U[1, 2] / (U[1, 1] * U[2, 2])]
    inv3 = [mpf(0), mpf(0), 1 / U[2, 2]]
    x = lambda row, j: sum(row[k] * ERt[k, j] for k in range(3))
    x10, x20, x21 = x(inv2, 0), x(inv3, 0), x(inv3, 1)
    W = mp.matrix([[0, -x10, -x20], [x10, 0, -x21], [x20, x21, 0]])
    return ERt - U * W, W * R


def kappas(M, U, R):
    """(kappa_U, kappa_R): the 2-norms of E -> dU and E -> dR over the nine unit directions,
    kappa_R times |M|_F (docs/NUMERICS.md §11)."""
    JU, JR = mp.matrix(9, 9), mp.matrix(9, 9)
    for c in range(9):
        E = mp.matrix(3, 3)
        E[c // 3, c % 3] = 1
        dU, dR = tangent(U, R, E)
        for r in range(9):
            JU[r, c], JR[r, c] = dU[r // 3, r % 3], dR[r // 3, r % 3]
    top = lambda J: max(mp.svd_r(J, compute_uv=False))
    norm_m = mp.sqrt(sum(M[i, j] ** 2 for i in range(3) for j in range(3)))
    return top(JU), top(JR) * norm_m


def factor(M, independent: bool):
    """RQ.1's construction of the exact M: rows of R from the last row, `completion` where a row
    is dependent, as docs/NUMERICS.md §22 completes it."""
    m = [[M[i, j] for j in range(3)] for i in range(3)]
    zero3 = all(x == 0 for x in m[2])
    r3 = _axis(2) if zero3 else _unit(m[2])
    w = [x - _dot(m[1], r3) * y for x, y in zip(m[1], r3, strict=True)]
    if independent:
        r2 = _unit(w)
        r1 = _cross(r2, r3)
    else:
        r2, r1 = completion(r3)
    U = mp.matrix(3, 3)
    U[0, 0], U[0, 1], U[0, 2] = _dot(m[0], r1), _dot(m[0], r2), _dot(m[0], r3)
    U[1, 1] = mp.sqrt(_dot(w, w)) if independent else mpf(0)
    U[1, 2] = _dot(m[1], r3)
    U[2, 2] = mpf(0) if zero3 else mp.sqrt(_dot(m[2], m[2]))
    return U, mp.matrix([r1, r2, r3])


def rq3(inputs: dict) -> dict:
    """RQ.1's construction of the stored M at the working precision. `valid` is det M > 0 in
    rationals. Where m2 and the reference's r3 are dependent (decided in rationals: rank of
    (m2, m3) below 2, or m3 = 0 and m2 a multiple of e3 or zero) the factorization is not unique,
    `completion` completes it, and both kappas are 0, weight 0 in the metric."""
    E = exact(inputs["M"])
    zero3 = all(x == 0 for x in E[2])
    independent = rank([E[1], E[2] if not zero3 else [0, 0, 1]]) == 2
    M = from_mat(inputs["M"])
    U, R = factor(M, independent)
    unique = independent and not zero3
    ku, kr = kappas(M, U, R) if unique else (mpf(0), mpf(0))
    return {
        "U": to_mat(U),
        "R": to_mat(R),
        "valid": [mpf(int(det3(E) > 0))],
        "kappa_U": [ku],
        "kappa_R": [kr],
    }
