"""`svd3` (docs/decisions/0101 decision 8): A = U diag(sigma) V^T of the stored 3x3 `A` in the
signed convention: U, V in SO(3), sigma1 >= sigma2 >= |sigma3|, sign sigma3 = sign det A.

A planted A is U diag(sigma) V^T at the working precision, U and V the rotations of Haar-random
quaternions, rounded once entrywise to the stratum's precision: the reference is the stored
matrix's decomposition, not the drawn sigma.
"""

from mpmath import mp, mpf

from . import so3
from .dense import exact, from_mat, gaussian, rank, round_matrix, to_mat
from .fmt import Mat
from .rng import unit_quaternion_s3
from .strata import Draw, deep, drawn_strata

N_SVD = deep(16)
N_WIDE = deep(32)  # `svd:random`, `svd:col-graded`, `svd:ordering`
N_ZERO = 4  # the four sign patterns of `_ZERO_SIGNS`, at every tier
# Two nonnegative singular values differ by less than |sigma|_2 unless the other two are 0, so a
# relative gap of 10^0 has no positive spectrum: the gaps start at 10^-1.
GAP_K = {False: (1, 2, 4, 6, 8, 10, 12), True: (1, 2, 4, 6)}
GRADED_K = {False: (2, 4, 8, 12), True: (2, 4, 6)}
COL_GRADES = {False: (0, 5, 10), True: (0, 2, 4)}  # B diag(10^-g) for `svd:col-graded`
NEAR_K = {False: (4, 8, 12), True: (2, 4, 6)}
SCALE = {False: 70, True: 12}  # svd:scale-up/down: 10^(+-SCALE)


def uniform(rng, lo, hi):
    return lo + (hi - lo) * mpf(rng.uniform())


def rotation(rng):
    return mp.matrix(so3.rot(list(unit_quaternion_s3(rng))))


def gaussian_matrix(rng):
    return mp.matrix([[gaussian(rng) for _ in range(3)] for _ in range(3)])


def with_sigma(spectrum):
    """A planted stratum: `spectrum(rng)` draws sigma from the stream `sigma`."""

    def make(draw: Draw, n: int) -> list[dict]:
        rq_u, rq_v, rs = draw.stream("U"), draw.stream("V"), draw.stream("sigma")
        out = []
        for _ in range(n):
            U, V = rotation(rq_u), rotation(rq_v)
            out.append({"A": round_matrix(U * mp.diag(spectrum(rs)) * V.T, draw.round)})
        return out

    return make


def _random(sign: int = 0):
    """`svd:random`, Gaussian entries, from its stream, times 10^(sign s), s the precision's
    `SCALE`."""

    def make(draw: Draw, n: int) -> list[dict]:
        rng = draw.stream("entries", of="svd:random")
        factor = mpf(10) ** (sign * SCALE[draw.f32])
        return [{"A": round_matrix(factor * gaussian_matrix(rng), draw.round)} for _ in range(n)]

    return make


def gap_spectrum(c, x, t, top: bool) -> list:
    """Positive singular values whose close pair (c - g/2, c + g/2) is g = t |sigma|_2 apart
    exactly (a quadratic in g). `bottom`: the third is c + g/2 + x beyond the pair; `top`: it is
    x (c - g/2) below it."""
    if top:
        a, b = [x * c, c, c], [-x / 2, -mpf(1) / 2, mpf(1) / 2]
    else:
        a, b = [c, c, c + x], [-mpf(1) / 2, mpf(1) / 2, mpf(1) / 2]
    dot = lambda u, v: sum(p * q for p, q in zip(u, v, strict=True))
    aa, ab, bb = dot(a, a), dot(a, b), dot(b, b)
    k = 1 - t * t * bb  # g^2 = t^2 (aa + 2 ab g + bb g^2)
    g = (t * t * ab + mp.sqrt(t**4 * ab * ab + k * t * t * aa)) / k
    return [p + q * g for p, q in zip(a, b, strict=True)]


def _gap(k: int, top: bool):
    def spectrum(rng):
        c = uniform(rng, mpf(1) / 2, 1)
        x = uniform(rng, mpf(1) / 4, mpf(3) / 4) if top else uniform(rng, mpf(1) / 4, 1)
        return gap_spectrum(c, x, mpf(10) ** -k, top)

    return with_sigma(spectrum)


def _graded(k: int):
    return with_sigma(lambda rng: [mpf(1), mpf(10) ** -k, mpf(10) ** -k / 2])


def _col_graded(draw: Draw, n: int) -> list[dict]:
    rng = draw.stream("entries")
    D = mp.diag([mpf(10) ** -g for g in COL_GRADES[draw.f32]])
    return [{"A": round_matrix(gaussian_matrix(rng) * D, draw.round)} for _ in range(n)]


# `svd:zero`'s sign patterns, column-major: all +0, all -0, then alternating from +0 and from -0.
_ZERO_SIGNS = ((1,) * 9, (-1,) * 9, (1, -1) * 4 + (1,), (-1, 1) * 4 + (-1,))


def _zero(draw: Draw, n: int) -> list[dict]:
    pattern = [_ZERO_SIGNS[i % len(_ZERO_SIGNS)] for i in range(n)]
    return [{"A": Mat((3, 3), [0.0 if s > 0 else -0.0 for s in signs])} for signs in pattern]


def _haar_rotation(draw: Draw, n: int) -> list[dict]:
    rng = draw.stream("quat")
    return [{"A": round_matrix(rotation(rng), draw.round)} for _ in range(n)]


def near_rotation(k: int, reflect: bool):
    """A Haar rotation plus Gaussian entries of size 10^-k; `reflect` times diag(1, 1, -1)."""

    def make(draw: Draw, n: int) -> list[dict]:
        rq, re = draw.stream("quat"), draw.stream("entries")
        flip = mp.diag([1, 1, -1 if reflect else 1])
        out = []
        for _ in range(n):
            R, E = rotation(rq), gaussian_matrix(re)
            out.append({"A": round_matrix((R + mpf(10) ** -k * E) * flip, draw.round)})
        return out

    return make


def _rank2(rng):
    return [mpf(1), uniform(rng, mpf(1) / 4, 1), mpf(0)]


def _essential(rng):
    s = uniform(rng, mpf(1) / 4, 1)
    return [s, s, mpf(0)]


def _ordering(rng):
    """diag(a, a, -b), b = a (1 + h) > a: PH.2's case, sigma = (b, a, -a)."""
    a = uniform(rng, mpf(1) / 4, 1)
    return [a, a, -a * (1 + uniform(rng, mpf(1) / 4, 4))]


def _named(f32: bool) -> list[tuple[str, object, int, bool]]:
    """(name, make, count, shared) in the catalogue's order; `shared`: the binary64 namesake's
    range, so the `@f32` stratum reads its stream."""
    gaps = [
        (f"svd:gap-1e-{k}/{side}", _gap(k, side == "top"), N_SVD, True)
        for side in ("bottom", "top")
        for k in GAP_K[f32]
    ]
    graded = [
        (f"svd:graded-1e-{k}", _graded(k), N_SVD, k in GRADED_K[False]) for k in GRADED_K[f32]
    ]
    near64 = NEAR_K[False]
    near = [
        (f"svd:near-{kind}-1e-{k}", near_rotation(k, kind == "reflection"), N_SVD, k in near64)
        for kind in ("rotation", "reflection")
        for k in NEAR_K[f32]
    ]
    return [
        ("svd:random", _random(), N_WIDE, True),
        *gaps,
        *graded,
        ("svd:col-graded", _col_graded, N_WIDE, False),
        ("svd:rank2", with_sigma(_rank2), N_SVD, True),
        ("svd:rank1", with_sigma(lambda rng: [mpf(1), mpf(0), mpf(0)]), N_SVD, True),
        ("svd:zero", _zero, N_ZERO, True),
        ("svd:rotation", _haar_rotation, N_SVD, True),
        *near,
        ("svd:essential", with_sigma(_essential), N_SVD, True),
        ("svd:ordering", with_sigma(_ordering), N_WIDE, True),
        ("svd:scale-up", _random(1), N_SVD, False),
        ("svd:scale-down", _random(-1), N_SVD, False),
    ]


def _strata():
    return drawn_strata(
        [(name, make, count) for name, make, count, _ in _named(False)],
        _named(True),
    )


STRATA = _strata()


# --- the reference ---------------------------------------------------------------------------


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _unit(x):
    r = mp.sqrt(sum(c * c for c in x))
    return [c / r for c in x]


def _complement(x):
    """A unit vector orthogonal to the unit `x`: x times the axis of its smallest component (the
    first of equals), normalized. A null space of dimension two has no vector svd_r would return
    the same at every precision."""
    axis = min(range(3), key=lambda r: (abs(x[r]), r))
    return _unit(_cross(x, [mpf(int(r == axis)) for r in range(3)]))


def _signed(u, v):
    """(u, v) jointly signed so that u's largest-magnitude component (the first of equals) is
    positive."""
    big = max(range(3), key=lambda r: (abs(u[r]), -r))
    return (u, v) if u[big] > 0 else ([-c for c in u], [-c for c in v])


def _columns(*cols) -> Mat:
    return Mat((3, 3), [c for col in cols for c in col])


def svd3(inputs: dict) -> dict:
    """`mp.svd_r` of the stored A, sorted descending. (u1, v1) and (u2, v2) signed by `_signed`,
    u3 = u1 x u2, v3 = v1 x v2, sigma3 = u3^T A v3: U, V in SO(3) and sign sigma3 = sign det A by
    construction. The singular values that are 0 by the rank of the stored A in rationals are 0
    (svd_r leaves 1e-121); at rank 1 (u2, v2) are `_complement`s, at rank 0 U = V = I."""
    r = rank(exact(inputs["A"]))
    if r == 0:
        eye = to_mat(mp.eye(3))
        return {"U": eye, "sigma": [mpf(0)] * 3, "V": eye}
    A = from_mat(inputs["A"])
    W, S, Vt = mp.svd_r(A)
    order = sorted(range(3), key=lambda i: -S[i])
    left = [[W[k, i] for k in range(3)] for i in order]
    right = [[Vt[i, k] for k in range(3)] for i in order]
    u1, v1 = _signed(left[0], right[0])
    u2, v2 = _signed(*((left[1], right[1]) if r >= 2 else (_complement(u1), _complement(v1))))
    u3, v3 = _cross(u1, u2), _cross(v1, v2)
    Av3 = A * mp.matrix(v3)
    sigma = [
        S[order[0]],
        S[order[1]] if r >= 2 else mpf(0),
        sum(u3[k] * Av3[k] for k in range(3)) if r == 3 else mpf(0),
    ]
    return {"U": _columns(u1, u2, u3), "sigma": sigma, "V": _columns(v1, v2, v3)}
