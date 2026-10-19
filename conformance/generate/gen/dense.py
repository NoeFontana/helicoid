"""n x n matrices of the 0056 ids: `mp.matrix` in the math, column-major `fmt.Mat` in records,
and the exact rational arithmetic their references decide things in."""

from collections.abc import Callable
from fractions import Fraction

from mpmath import mp, mpf

from .fmt import Mat
from .rng import SplitMix64


def to_mat(A) -> Mat:
    """An `mp.matrix` as a record's column-major `Mat`."""
    return Mat((A.rows, A.cols), [A[i, j] for j in range(A.cols) for i in range(A.rows)])


def from_mat(m: Mat):
    rows, cols = m.shape
    return mp.matrix([[mpf(m.data[rows * j + i]) for j in range(cols)] for i in range(rows)])


def exact(m: Mat) -> list[list[Fraction]]:
    """The stored matrix's rows, exactly."""
    rows, cols = m.shape
    return [[Fraction(m.data[rows * j + i]) for j in range(cols)] for i in range(rows)]


def round_symmetric(A, rnd: Callable) -> Mat:
    """A symmetric matrix rounded by `rnd`: the upper triangle once, mirrored, so the stored
    matrix is exactly symmetric."""
    n = A.rows
    upper = {(i, j): rnd(A[i, j]) for i in range(n) for j in range(i, n)}
    return Mat((n, n), [upper[min(i, j), max(i, j)] for j in range(n) for i in range(n)])


def round_matrix(A, rnd: Callable) -> Mat:
    """A general matrix rounded entrywise by `rnd`, once."""
    return Mat((A.rows, A.cols), [rnd(A[i, j]) for j in range(A.cols) for i in range(A.rows)])


def rank(A: list[list[Fraction]]) -> int:
    """The rank of an exact matrix, by Gauss-Jordan elimination in rationals."""
    rows = [row[:] for row in A]
    n, out = len(rows[0]) if rows else 0, 0
    for col in range(n):
        pivot = next((r for r in range(out, len(rows)) if rows[r][col] != 0), None)
        if pivot is None:
            continue
        rows[out], rows[pivot] = rows[pivot], rows[out]
        for r in range(len(rows)):
            if r != out and rows[r][col] != 0:
                f = rows[r][col] / rows[out][col]
                rows[r] = [x - f * y for x, y in zip(rows[r], rows[out], strict=True)]
        out += 1
    return out


def frobenius(A):
    return mp.sqrt(sum(abs(A[i, j]) ** 2 for i in range(A.rows) for j in range(A.cols)))


def gaussian(rng: SplitMix64):
    """A standard normal at the working precision (Box-Muller, the cosine half)."""
    u1, u2 = 1 - mpf(rng.uniform()), mpf(rng.uniform())  # u1 in (0, 1]
    return mp.sqrt(-2 * mp.log(u1)) * mp.cos(2 * mp.pi * u2)


def unit_vector(rng: SplitMix64, n: int) -> list:
    """Uniform on S^(n-1), unrounded: normalized Gaussians."""
    while True:
        v = [gaussian(rng) for _ in range(n)]
        r = mp.sqrt(sum(c * c for c in v))
        if r > mpf(10) ** -10:
            return [c / r for c in v]


def ldlt_pivots(A: list[list[Fraction]]) -> list[Fraction] | None:
    """The pivots of an exact LDL^T of a symmetric A, or None at the first one that is not
    positive: A is positive definite exactly when this is not None."""
    n = len(A)
    L = [[Fraction(0)] * n for _ in range(n)]
    d: list[Fraction] = []
    for k in range(n):
        dk = A[k][k] - sum(L[k][j] ** 2 * d[j] for j in range(k))
        if dk <= 0:
            return None
        d.append(dk)
        for i in range(k + 1, n):
            L[i][k] = (A[i][k] - sum(L[i][j] * L[k][j] * d[j] for j in range(k))) / dk
    return d
