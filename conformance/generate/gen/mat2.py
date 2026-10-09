"""`mat2_inverse_adj` (docs/decisions/0061): `(adj(A)/det, det)` of the stored 2x2 `A`.

`cond:1`, `cond:1e4`, `cond:1e8`: A = s R(alpha) diag(1, 1/kappa) R(beta)^T at the working precision,
alpha and beta uniform on [0, 2 pi), s = 10^U(-3, 3), each entry rounded once to the stratum's
precision. The reference is the inverse of the rounded entries, exactly: the stratum names the
draw's conditioning and the reference is the function at the input given.

`cond:singular`: A = [[x, 2^j x], [y, 2^j y]], x and y rounded first, j uniform in {-3, ..., 3}.
The two products of `det` are the same power of two times `x y`, so they round equally and `det`
is exactly 0 in both precisions. The reference is `det = 0` and `inv` zero, an inverse that does
not exist (NUMERICS.md section 11 expects the subject's to be non-finite).
"""

from fractions import Fraction

from mpmath import mp, mpf

from .dense import exact, to_mat
from .fmt import Mat
from .strata import Draw, drawn_strata

N_MAT2 = 6
# cond:1e{k}, cond:1 at k = 0. binary32 stops at 1e4: at kappa = 1e8 > 1/u the two products of `det`
# round equal on a nonsingular input (4 of 6 records), and a zero `det` there is a non-finite
# inverse that section 11 counts as a failure; chol's `COND_K` narrows its range at f32 the same way.
KAPPA = {False: (0, 4, 8), True: (0, 4)}


def _uniform(rng, lo, hi):
    return lo + (hi - lo) * mpf(rng.uniform())


def _rotation(angle):
    c, s = mp.cos(angle), mp.sin(angle)
    return mp.matrix([[c, -s], [s, c]])


def _scale(rng):
    return mpf(10) ** _uniform(rng, -3, 3)


def _conditioned(k: int):
    def make(draw: Draw, count: int) -> list[dict]:
        ra, rs = draw.stream("angles"), draw.stream("scale")
        out = []
        for _ in range(count):
            alpha, beta = _uniform(ra, 0, 2 * mp.pi), _uniform(ra, 0, 2 * mp.pi)
            A = _scale(rs) * _rotation(alpha) * mp.diag([1, mpf(10) ** -k]) * _rotation(beta).T
            out.append({"A": Mat((2, 2), [draw.round(A[i, j]) for j in range(2) for i in range(2)])})
        return out

    return make


def _singular(draw: Draw, count: int) -> list[dict]:
    rx, rs, rj = draw.stream("entries"), draw.stream("scale"), draw.stream("power")
    out = []
    for _ in range(count):
        x, y = (draw.round(_scale(rs) * _uniform(rx, -1, 1)) for _ in range(2))
        p = 2.0 ** (int(rj.next_u64() % 7) - 3)
        out.append({"A": Mat((2, 2), [x, y, p * x, p * y])})
    return out


def _named(f32: bool):
    named = [(f"cond:{'1' if k == 0 else f'1e{k}'}", _conditioned(k)) for k in KAPPA[f32]]
    return [*named, ("cond:singular", _singular)]


def _strata():
    return drawn_strata(
        [(name, make, N_MAT2) for name, make in _named(False)],
        [(name, make, N_MAT2, True) for name, make in _named(True)],
    )


STRATA = _strata()


def determinant(A: list[list[Fraction]]) -> Fraction:
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


def inverse_adj(inputs: dict) -> dict:
    """The stored A's determinant and inverse in rationals, exactly; zeros where `det` is 0."""
    A = exact(inputs["A"])
    det = determinant(A)
    adj = [[A[1][1], -A[0][1]], [-A[1][0], A[0][0]]]
    inv = mp.matrix(2, 2)
    if det != 0:
        for i in range(2):
            for j in range(2):
                q = adj[i][j] / det
                inv[i, j] = mpf(q.numerator) / q.denominator
    return {"det": [mpf(det.numerator) / det.denominator], "inv": to_mat(inv)}
