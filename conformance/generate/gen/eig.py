"""`eig3` (docs/decisions/0056): the eigen-decomposition of the stored symmetric 3x3 `A`, by
`mp.eigsy`, eigenvalues ascending, `V`'s columns the eigenvectors, each signed so that its
largest-magnitude component (the first of equals) is positive.

A = Q diag(lambda) Q^T at the working precision, Q the rotation of a Haar-random quaternion, then
rounded once to the stratum's precision (the upper triangle, mirrored): the stored matrix is
exactly symmetric, and its spectrum is the reference's, not the drawn lambda. `eig:triple` is c I
itself: Q c I Q^T at 120 digits would round to c I plus off-diagonal entries of 1e-121.
"""

from fractions import Fraction

from mpmath import mp, mpf

from . import so3
from .dense import exact, from_mat, round_symmetric, to_mat
from .fmt import Mat
from .rng import unit_quaternion_s3
from .strata import Draw, drawn_strata

N_RANDOM = 64
N_GAP = 32
GAP_K = {False: range(13), True: range(7)}
SCALE = {False: 70, True: 8}  # eig:scale-up/down: 10^(+-SCALE), the rustdoc's narrower range


def _uniform(rng, lo, hi):
    return lo + (hi - lo) * mpf(rng.uniform())


def _rotation(rng):
    return mp.matrix(so3.rot(list(unit_quaternion_s3(rng))))


def _sandwich(Q, lam):
    return Q * mp.diag(lam) * Q.T


def _random(sign: int = 0):
    """`eig:random`, from its stream, times 10^(sign s), s the precision's `SCALE`."""

    def make(draw: Draw, n: int) -> list[dict]:
        rq, rl = draw.stream("quat", of="eig:random"), draw.stream("lambda", of="eig:random")
        factor = mpf(10) ** (sign * SCALE[draw.f32])
        out = []
        for _ in range(n):
            Q, lam = _rotation(rq), [_uniform(rl, -1, 1) for _ in range(3)]
            out.append({"A": round_symmetric(factor * _sandwich(Q, lam), draw.round)})
        return out

    return make


def gap_spectrum(c, h, t, top: bool) -> list:
    """(c - g/2, c + g/2, c + g/2 + h) with g = t |lambda|_2 exactly (a quadratic in g); `top`
    negates it, so that the close pair is the largest two."""
    s = c + h
    k = 1 - 3 * t * t / 4  # g^2 = t^2 (2 c^2 + s^2 + s g + 3 g^2 / 4)
    g = (t * t * s + mp.sqrt(t**4 * s * s + 4 * k * t * t * (2 * c * c + s * s))) / (2 * k)
    lam = [c - g / 2, c + g / 2, c + g / 2 + h]
    return sorted(-x for x in lam) if top else lam


def _gap(k: int, top: bool):
    def make(draw: Draw, n: int) -> list[dict]:
        rq, rl = draw.stream("quat"), draw.stream("lambda")
        out = []
        for _ in range(n):
            Q = _rotation(rq)
            c, h = _uniform(rl, -1, 1), _uniform(rl, mpf(1) / 4, 1)
            lam = gap_spectrum(c, h, mpf(10) ** -k, top)
            out.append({"A": round_symmetric(_sandwich(Q, lam), draw.round)})
        return out

    return make


def _triple(draw: Draw, n: int) -> list[dict]:
    rl = draw.stream("lambda")
    out = []
    for _ in range(n):
        c = draw.round(_uniform(rl, -1, 1))
        out.append({"A": Mat((3, 3), [c if i % 4 == 0 else 0.0 for i in range(9)])})
    return out


def _rank1(draw: Draw, n: int) -> list[dict]:
    rq, rl = draw.stream("quat"), draw.stream("lambda")
    out = []
    for _ in range(n):
        Q, c = _rotation(rq), _uniform(rl, -1, 1)
        out.append({"A": round_symmetric(_sandwich(Q, [0, 0, c]), draw.round)})
    return out


def _gap_names(k_range) -> list[tuple[str, object]]:
    return [
        (f"eig:gap-1e-{k}/{side}", _gap(k, side == "top"))
        for side in ("bottom", "top")
        for k in k_range
    ]


def _strata():
    head = [("eig:random", _random(), N_RANDOM, True)]
    tail = [
        ("eig:triple", _triple, N_RANDOM, True),
        ("eig:rank1", _rank1, N_RANDOM, True),
        ("eig:scale-up", _random(1), N_RANDOM, False),
        ("eig:scale-down", _random(-1), N_RANDOM, False),
    ]
    gaps64 = [(name, make, N_GAP, True) for name, make in _gap_names(GAP_K[False])]
    gaps32 = [(name, make, N_GAP, True) for name, make in _gap_names(GAP_K[True])]
    return drawn_strata(
        [(name, make, count) for name, make, count, _ in (*head, *gaps64, *tail)],
        [*head, *gaps32, *tail],
    )


STRATA = _strata()


# --- the reference ---------------------------------------------------------------------------


def _rank(A: list[list[Fraction]]) -> int:
    rows = [row[:] for row in A]
    rank = 0
    for col in range(3):
        pivot = next((r for r in range(rank, 3) if rows[r][col] != 0), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for r in range(3):
            if r != rank and rows[r][col] != 0:
                f = rows[r][col] / rows[rank][col]
                rows[r] = [x - f * y for x, y in zip(rows[r], rows[rank], strict=True)]
        rank += 1
    return rank


def eig3(inputs: dict) -> dict:
    """`mp.eigsy` of the stored A. An eigenvalue that is exactly 0 (A exactly singular, decided in
    rationals) is set to 0: `eigsy` leaves 1e-121 there, which no 30-digit value agrees with."""
    A = from_mat(inputs["A"])
    E, Q = mp.eigsy(A)
    order = sorted(range(3), key=lambda i: E[i])
    lam = [E[i] for i in order]
    zeros = 3 - _rank(exact(inputs["A"]))
    for i in sorted(range(3), key=lambda i: abs(lam[i]))[:zeros]:
        lam[i] = mpf(0)
    V = mp.matrix(3, 3)
    for j, i in enumerate(order):
        col = [Q[r, i] for r in range(3)]
        big = max(range(3), key=lambda r: (abs(col[r]), -r))
        sign = -1 if col[big] < 0 else 1
        for r in range(3):
            V[r, j] = sign * col[r]
    return {"lambda": lam, "V": to_mat(V)}
