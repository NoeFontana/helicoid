"""`nearest_rotation` (docs/decisions/0102 decision 8): the Frobenius-nearest rotation R of
the stored 3x3 `M`, U diag(1, 1, d) V^T, d = sign det M.

A planted M is drawn at the working precision and rounded once entrywise to the stratum's
precision, as `svd3`'s are: the reference is the stored matrix's nearest rotation, not the drawn
one. `polar:rank2` is built exact instead, and `polar:exact` is exact already.
"""

from fractions import Fraction

from mpmath import mp, mpf

from . import so3
from .cubic import discriminant, multiple_roots
from .dense import det3, exact, from_mat, round_matrix, to_mat
from .fmt import Mat
from .strata import Draw, deep, drawn_strata
from .svd import gaussian_matrix, near_rotation, rotation, uniform, with_sigma

N_POLAR = deep(16)
N_RANDOM = deep(32)  # `polar:random`
N_EXACT = so3.N_EXACT  # the 24 signed permutation rotations, at every tier
NEAR_K = {False: (2, 4, 8, 12, 16), True: (2, 4, 6)}
# sigma = (1, 1/2, +-10^-k): det M / |M|_F^3 is 0.36 10^-k, so k = 4 is above 0102's tau = 2^-26
# and k = 8 below it at binary64, k = 2 above 2^-12 and k = 4 below at binary32.
COND_K = {False: (4, 8, 12, 16), True: (2, 4, 6, 8)}
SCALE = {False: 70, True: 12}  # polar:scale-up/down: 10^(+-SCALE)
# `polar:rank2`'s grid: X Y^T with entries of X and Y in 2^-p Z sums two products below 2^(2p+1)
# 2^-2p, exact in the precision's 53 or 24 bits.
RANK2_BITS = {False: 26, True: 11}


def _as_m(make):
    """An `svd3` stratum's draw, its `A` stored as `M`."""

    def renamed(draw: Draw, n: int) -> list[dict]:
        return [{"M": record["A"]} for record in make(draw, n)]

    return renamed


def _exact(draw: Draw, n: int) -> list[dict]:
    return [{"M": record["R"]} for record in so3.exact_rotations(draw, n)]


def _scaled(draw: Draw, n: int) -> list[dict]:
    """c R, R Haar, c log-uniform in [10^-3, 10^3] at the working precision."""
    rq, rc = draw.stream("quat"), draw.stream("scale")
    out = []
    for _ in range(n):
        R, c = rotation(rq), mpf(10) ** (6 * mpf(rc.uniform()) - 3)
        out.append({"M": round_matrix(c * R, draw.round)})
    return out


def _random(sign: int = 0):
    """`polar:random`, Gaussian entries, from its stream, times 10^(sign s), s the precision's
    `SCALE`."""

    def make(draw: Draw, n: int) -> list[dict]:
        rng = draw.stream("entries", of="polar:random")
        factor = mpf(10) ** (sign * SCALE[draw.f32])
        return [{"M": round_matrix(factor * gaussian_matrix(rng), draw.round)} for _ in range(n)]

    return make


def _reflection(rng):
    """diag(1, a, -b), a in [1/2, 1], b = h a, h in [1/8, 1/2): det M / |M|_F^3 > 0.009, past
    0102's tau at both precisions, and sigma2 - |sigma3| > a/2."""
    a = uniform(rng, mpf(1) / 2, 1)
    return [mpf(1), a, -a * uniform(rng, mpf(1) / 8, mpf(1) / 2)]


def _cond(k: int):
    """sigma = (1, 1/2, s 10^-k), s = +1 on even records and -1 on odd: both signs of det."""

    def make(draw: Draw, n: int) -> list[dict]:
        rq_u, rq_v = draw.stream("U"), draw.stream("V")
        out = []
        for i in range(n):
            U, V = rotation(rq_u), rotation(rq_v)
            sigma = [mpf(1), mpf(1) / 2, (-1) ** i * mpf(10) ** -k]
            out.append({"M": round_matrix(U * mp.diag(sigma) * V.T, draw.round)})
        return out

    return make


def _rank2(draw: Draw, n: int) -> list[dict]:
    """sigma = (1, 1/2, 0) exactly in rank: X Y^T, X = [u1, u2] and Y = [v1, v2/2] of Haar U, V,
    each entry rounded to the grid 2^-p, the product exact (`RANK2_BITS`). Its singular values are
    (1, 1/2) to about 2^-p."""
    p = RANK2_BITS[draw.f32]
    rq_u, rq_v = draw.stream("U"), draw.stream("V")
    grid = lambda x: int(mp.nint(x * 2**p))
    out = []
    for _ in range(n):
        U, V = rotation(rq_u), rotation(rq_v)
        X = [[grid(U[i, 0]), grid(U[i, 1])] for i in range(3)]
        Y = [[grid(V[j, 0]), grid(V[j, 1] / 2)] for j in range(3)]
        data = []
        for j in range(3):
            for i in range(3):
                num = X[i][0] * Y[j][0] + X[i][1] * Y[j][1]
                x = draw.round(mp.ldexp(mpf(num), -2 * p))
                if Fraction(x) != Fraction(num, 1 << 2 * p):
                    raise ArithmeticError(f"polar:rank2: {num} 2^-{2 * p} is not exact")
                data.append(x)
        out.append({"M": Mat((3, 3), data)})
    return out


def _non_unique(rng):
    """diag(1, a, -a), a in [1/4, 1]: sigma2 + d sigma3 = 0, a circle of nearest rotations. The
    stored M's rounding breaks the tie, by 1e-18 to 1e-16 (3e-9 to 4e-8 at binary32), so its
    nearest rotation is unique, a point of the circle the rounding picks; `unique` holds it."""
    a = uniform(rng, mpf(1) / 4, 1)
    return [mpf(1), a, -a]


def _named(f32: bool) -> list[tuple[str, object, int, bool]]:
    """(name, make, count, shared) in the catalogue's order; `shared`: the binary64 namesake's
    range, so the `@f32` stratum reads its stream."""
    near = [
        (f"polar:near-1e-{k}", _as_m(near_rotation(k, False)), N_POLAR, k in NEAR_K[False])
        for k in NEAR_K[f32]
    ]
    cond = [(f"polar:cond-1e-{k}", _cond(k), N_POLAR, k in COND_K[False]) for k in COND_K[f32]]
    return [
        *near,
        ("polar:exact", _exact, N_EXACT, True),
        ("polar:scaled", _scaled, N_POLAR, True),
        ("polar:random", _random(), N_RANDOM, True),
        ("polar:reflection", _as_m(with_sigma(_reflection)), N_POLAR, True),
        *cond,
        ("polar:rank2", _rank2, N_POLAR, True),
        ("polar:non-unique", _as_m(with_sigma(_non_unique)), N_POLAR, True),
        ("polar:scale-up", _random(1), N_POLAR, False),
        ("polar:scale-down", _random(-1), N_POLAR, False),
    ]


def _strata():
    return drawn_strata(
        [(name, make, count) for name, make, count, _ in _named(False)],
        _named(True),
    )


STRATA = _strata()


# --- the reference ---------------------------------------------------------------------------


def _sign(x) -> int:
    return (x > 0) - (x < 0)


def _is_rotation(M: list[list[Fraction]]) -> bool:
    gram = lambda i, j: sum(M[k][i] * M[k][j] for k in range(3))
    return all(gram(i, j) == int(i == j) for i in range(3) for j in range(3)) and det3(M) == 1


def unique(M: list[list[Fraction]]) -> bool:
    """Whether the nearest rotation of the exact M is unique: unless det M <= 0 and the two
    smallest eigenvalues of M^T M tie (sigma2 + d sigma3 = 0, SO.13(c)), decided in rationals by
    the discriminant of its characteristic polynomial and, where that is 0, its rational roots."""
    S = [[sum(M[k][i] * M[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    minors = sum(S[i][i] * S[j][j] - S[i][j] * S[j][i] for i in range(3) for j in range(i + 1, 3))
    coeffs = [Fraction(1), -(S[0][0] + S[1][1] + S[2][2]), minors, -det3(S)]
    if det3(M) > 0 or discriminant(*coeffs) != 0:
        return True
    lo, mid, _ = sorted(multiple_roots(coeffs))
    return lo != mid


def nearest_rotation(inputs: dict) -> dict:
    """U diag(1, 1, d) V^T from `mp.svd_r` of the stored M, sorted descending, d the sign of det M
    in rationals; where det M is 0 (`polar:rank2`), d = det U det V, which makes R the one
    nearest rotation whatever signs svd_r gives the null vectors. An M that is exactly a rotation
    is R itself: svd_r of a signed permutation leaves 1e-121 where R has its exact zeros. An M
    whose nearest rotation is not unique is refused: svd_r's choice on the circle is not the same
    at every precision, and no stratum stores one (`polar:non-unique`'s rounding breaks the tie)."""
    E = exact(inputs["M"])
    if _is_rotation(E):
        return {"R": Mat((3, 3), [mpf(x) for x in inputs["M"].data])}
    if not unique(E):
        raise ValueError(f"nearest_rotation: not unique for {inputs['M'].data}")
    W, S, Vt = mp.svd_r(from_mat(inputs["M"]))
    order = sorted(range(3), key=lambda i: -S[i])
    U = mp.matrix([[W[k, i] for i in order] for k in range(3)])
    V = mp.matrix([[Vt[i, k] for i in order] for k in range(3)])
    d = _sign(det3(E)) or _sign(mp.det(U) * mp.det(V))
    return {"R": to_mat(U * mp.diag([1, 1, d]) * V.T)}
