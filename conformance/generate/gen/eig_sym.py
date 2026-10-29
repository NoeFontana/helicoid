"""`eig_sym_n{4,9,16}` (docs/decisions/0111): the eigendecomposition of the stored symmetric n x n
`A` by `mp.eigsy`, eigenvalues ascending, `V`'s columns the eigenvectors, each signed so that its
largest-magnitude component (the first of equals) is positive, as `eig3`'s; and `kappa`, the
2-norm condition number of D^-1 A D^-1, D = diag(sqrt(a_ii)), where the stored A is positive
definite (decided exactly), else 0. `kappa` is what the relative reading of docs/NUMERICS.md
section 11 divides by.

A = Q diag(lambda) Q^T at the working precision, Q a product of n Householder reflections from
directions uniform on S^(n-1) (`chol`'s), then rounded once to the stratum's precision (the upper
triangle, mirrored): the stored matrix is exactly symmetric, and its spectrum is the reference's,
not the drawn lambda. Every stream's purpose names n, so the three ids draw independently.
"""

from fractions import Fraction

from mpmath import mp, mpf

from .chol import _householder
from .dense import exact, from_mat, ldlt_pivots, rank, round_symmetric, to_mat
from .fmt import Mat
from .strata import Draw, deep, drawn_strata

NS = (4, 9, 16)
# Records per stratum (docs/decisions/0110 allots the family 8.2 MB): the code path is the same at
# every n, so the small n carry the counts and n = 16 its size.
COUNT = {4: deep(32), 9: deep(16), 16: deep(4)}
GAP_K = {False: (1, 4, 8, 12), True: (1, 3, 5)}
# The steepest reach below sqrt(tiny), where `eig_sym` scales a pair before `angle` (0111): a
# diagonal spanning 10^-200 at binary64, 10^-30 at binary32.
GRADED_K = {False: (4, 8, 16, 200), True: (4, 8, 30)}
UNIT_DIAG_K = {False: (3, 6, 9, 12), True: (2, 4)}
SCALE = {False: 70, True: 8}  # `eig3`'s
RANK_DEF = (1, 3)
INTEGER = 8  # rank-def: B's entries are integers in [-INTEGER, INTEGER]


def _uniform(rng, lo, hi):
    return lo + (hi - lo) * mpf(rng.uniform())


def _sandwich(Q, lam):
    return Q * mp.diag(lam) * Q.T


def _spectrum(rng, n: int) -> list:
    return sorted(_uniform(rng, -1, 1) for _ in range(n))


def _planted(n: int, spectrum, of: str | None = None):
    """Q diag(spectrum(rng)) Q^T, rounded; `of` reads another stratum's streams."""

    def make(draw: Draw, count: int) -> list[dict]:
        rq, rl = draw.stream(f"householder-n{n}", of), draw.stream(f"lambda-n{n}", of)
        return [
            {"A": round_symmetric(_sandwich(_householder(rq, n), spectrum(rl)), draw.round)}
            for _ in range(count)
        ]

    return make


def _gap(n: int, k: int, top: bool):
    """A uniform spectrum with its two smallest (`top`: largest) 10^-k apart."""

    def spectrum(rng):
        lam = _spectrum(rng, n)
        g = mpf(10) ** -k
        if top:
            lam[-2] = lam[-1] - g
        else:
            lam[1] = lam[0] + g
        return lam

    return _planted(n, spectrum)


def _cluster(n: int, m: int):
    """A uniform spectrum whose smallest m are one value: an m-fold eigenvalue, exact before the
    rounding of A."""

    def spectrum(rng):
        lam = _spectrum(rng, n)
        return [lam[0]] * m + lam[m:]

    return _planted(n, spectrum)


def _rank_def(n: int, r: int):
    """A = B B^T, B an integer n x (n - r): exactly representable at both precisions (entries at
    most n INTEGER^2), so the nullity is exactly r, which `rank` checks."""

    def make(draw: Draw, count: int) -> list[dict]:
        rb = draw.stream(f"integers-n{n}")
        out = []
        for _ in range(count):
            B = [[int(rb.next_u64() % (2 * INTEGER + 1)) - INTEGER for _ in range(n - r)] for _ in range(n)]
            A = [[sum(B[i][k] * B[j][k] for k in range(n - r)) for j in range(n)] for i in range(n)]
            m = Mat((n, n), [float(A[i][j]) for j in range(n) for i in range(n)])
            if rank(exact(m)) != n - r:
                raise ArithmeticError(f"eig_sym_n{n} {draw.name}: B is not of full column rank")
            out.append({"A": m})
        return out

    return make


def _equilibrated(A):
    """D^-1 A D^-1, D = diag(sqrt(a_ii)): the unit-diagonal matrix of A's relative perturbation
    theory (Demmel and Veselic)."""
    n = A.rows
    d = [1 / mp.sqrt(A[i, i]) for i in range(n)]
    return mp.matrix([[d[i] * A[i, j] * d[j] for j in range(n)] for i in range(n)])


def _graded(n: int, k: int):
    """D A D: A positive definite with a unit diagonal and eigenvalues in [1/4, 7/4] before the
    equilibration, D geometric from 1 to 10^-(k/2), so A's diagonal spans 10^-k."""

    def make(draw: Draw, count: int) -> list[dict]:
        rq, rl = draw.stream(f"householder-n{n}"), draw.stream(f"lambda-n{n}")
        d = [mpf(10) ** (-mpf(k) * i / (2 * (n - 1))) for i in range(n)]
        out = []
        for _ in range(count):
            H = _equilibrated(_sandwich(_householder(rq, n), [_uniform(rl, 0.25, 1.75) for _ in range(n)]))
            G = mp.matrix([[d[i] * H[i, j] * d[j] for j in range(n)] for i in range(n)])
            out.append({"A": round_symmetric(G, draw.round)})
        return out

    return make


def correlation(Q, lam):
    """The unit-diagonal matrix with spectrum `lam` (trace n), from Q diag(lam) Q^T by n - 1 plane
    rotations, each setting one diagonal entry to 1 (Bendel and Mickey; Davies and Higham 2000,
    section 3): the rotation of (i, j), a_ii < 1 < a_jj, by t, the root of
    t^2 (a_jj - 1) - 2 t a_ij + (a_ii - 1) = 0 taken without cancellation."""
    A = _sandwich(Q, lam)
    n = A.rows
    done = mpf(10) ** (10 - mp.dps)  # an entry already 1 to the working precision
    for _ in range(n - 1):
        below = [i for i in range(n) if A[i, i] < 1 - done]
        above = [i for i in range(n) if A[i, i] > 1 + done]
        if not below or not above:
            break
        i, j = below[0], above[0]
        aii, ajj, aij = A[i, i] - 1, A[j, j] - 1, A[i, j]
        root = mp.sqrt(aij * aij - aii * ajj)
        t = aii / (aij + root) if aij >= 0 else aii / (aij - root)
        c = 1 / mp.sqrt(1 + t * t)
        s = c * t
        R = mp.eye(n)
        R[i, i], R[i, j], R[j, i], R[j, j] = c, -s, s, c
        A = R * A * R.T
    return A


def _unit_diag(n: int, k: int):
    """A correlation matrix of condition number 10^k: spectrum log-spaced from 1 to 10^-k, scaled
    to trace n, made unit-diagonal by `correlation`. Its diagonal rounds to 1 exactly."""

    def make(draw: Draw, count: int) -> list[dict]:
        rq = draw.stream(f"householder-n{n}")
        lam = [mpf(10) ** (-mpf(k) * i / (n - 1)) for i in range(n)]
        lam = [x * n / sum(lam) for x in lam]
        return [
            {"A": round_symmetric(correlation(_householder(rq, n), lam), draw.round)} for _ in range(count)
        ]

    return make


def _diagonal(n: int):
    """A uniform diagonal with two exact ties, at positions 0 and 2, and 1 and n - 1."""

    def make(draw: Draw, count: int) -> list[dict]:
        rl = draw.stream(f"lambda-n{n}")
        out = []
        for _ in range(count):
            d = [draw.round(_uniform(rl, -1, 1)) for _ in range(n)]
            d[2], d[n - 1] = d[0], d[1]
            out.append({"A": Mat((n, n), [d[i] if i == j else 0.0 for j in range(n) for i in range(n)])})
        return out

    return make


def _identity_multiple(n: int):
    def make(draw: Draw, count: int) -> list[dict]:
        rl = draw.stream(f"lambda-n{n}")
        out = []
        for _ in range(count):
            c = draw.round(_uniform(rl, -1, 1))
            out.append({"A": Mat((n, n), [c if i == j else 0.0 for j in range(n) for i in range(n)])})
        return out

    return make


def _scaled(n: int, sign: int):
    """`eig:random`'s draws times 10^(sign s), s the precision's `SCALE`."""

    def make(draw: Draw, count: int) -> list[dict]:
        factor = mpf(10) ** (sign * SCALE[draw.f32])
        return _planted(n, lambda rng: [factor * x for x in _spectrum(rng, n)], of="eig:random")(draw, count)

    return make


def _strata_of(n: int, f32: bool) -> list[tuple[str, object, bool]]:
    """(name, make, shared) in the catalogue's order."""
    return [
        ("eig:random", _planted(n, lambda rng: _spectrum(rng, n)), True),
        *(
            (f"eig:gap-1e-{k}/{side}", _gap(n, k, side == "top"), k in GAP_K[False])
            for side in ("bottom", "top")
            for k in GAP_K[f32]
        ),
        *((f"eig:cluster-{m}", _cluster(n, m), True) for m in (2, 3)),
        *((f"eig:rank-def-{r}", _rank_def(n, r), True) for r in RANK_DEF),
        *((f"eig:graded-1e-{k}", _graded(n, k), k in GRADED_K[False]) for k in GRADED_K[f32]),
        *(
            (f"eig:unit-diag-cond-1e-{k}", _unit_diag(n, k), k in UNIT_DIAG_K[False])
            for k in UNIT_DIAG_K[f32]
        ),
        ("eig:diagonal", _diagonal(n), True),
        ("eig:identity-multiple", _identity_multiple(n), True),
        ("eig:scale-up", _scaled(n, 1), False),
        ("eig:scale-down", _scaled(n, -1), False),
    ]


def strata(n: int):
    return drawn_strata(
        [(name, make, COUNT[n]) for name, make, _ in _strata_of(n, False)],
        [(name, make, COUNT[n], shared) for name, make, shared in _strata_of(n, True)],
    )


# --- the reference ---------------------------------------------------------------------------


def kappa(m: Mat):
    """The 2-norm condition number of the stored A's equilibration where A is exactly positive
    definite, else 0."""
    if ldlt_pivots(exact(m)) is None:
        return mpf(0)
    E = mp.eigsy(_equilibrated(from_mat(m)), eigvals_only=True)
    lo, hi = min(E), max(E)
    return hi / lo


def nullspace(A: list[list[Fraction]]) -> list[list[Fraction]]:
    """A basis of the exact null space of A, one vector per free column of its reduced row echelon
    form, by Gauss-Jordan elimination in rationals: the same basis at every working precision."""
    rows = [row[:] for row in A]
    n = len(rows[0]) if rows else 0
    pivots: list[int] = []
    for col in range(n):
        r = len(pivots)
        pivot = next((i for i in range(r, len(rows)) if rows[i][col] != 0), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        rows[r] = [x / rows[r][col] for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][col] != 0:
                f = rows[i][col]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[r], strict=True)]
        pivots.append(col)
    basis = []
    for free in (c for c in range(n) if c not in pivots):
        v = [Fraction(0)] * n
        v[free] = Fraction(1)
        for r, p in enumerate(pivots):
            v[p] = -rows[r][free]
        basis.append(v)
    return basis


def _orthonormal(vectors: list[list[Fraction]]) -> list:
    """Gram-Schmidt, twice, at the working precision."""
    out: list = []
    for v in vectors:
        w = [mpf(x) for x in v]
        for _ in range(2):
            for u in out:
                d = sum(a * b for a, b in zip(w, u, strict=True))
                w = [a - d * b for a, b in zip(w, u, strict=True)]
        r = mp.sqrt(sum(a * a for a in w))
        out.append([a / r for a in w])
    return out


def _span(m: Mat) -> int:
    """The decimal digits between the largest and smallest nonzero |a_ii|: `eigsy` is accurate to
    the working precision of |A|, and a positive definite A's smallest eigenvalue is at least its
    smallest diagonal entry over the equilibrated condition number, which the working precision
    already covers."""
    d = [abs(mpf(m.data[i * (m.shape[0] + 1)])) for i in range(m.shape[0])]
    d = [x for x in d if x != 0]
    return int(mp.ceil(mp.log10(max(d) / min(d)))) if d else 0


def eig_sym(inputs: dict) -> dict:
    """`mp.eigsy` of the stored A, at the working precision plus the digits its diagonal spans, so
    a graded A's smallest eigenvalues are as accurate as its largest. Where A is exactly singular (decided in rationals; only when the
    smallest |eigenvalue| is below 10^-(DPS/2) of the largest, so a regular A pays no rational
    elimination), the eigenvalues nearest 0 are set to 0: `eigsy` leaves 1e-121 there, which no
    30-digit value agrees with. A null space of two or more dimensions has no unique basis, and
    `eigsy` picks another at every precision, so its columns are `nullspace`'s, orthonormalized,
    in the zero eigenvalues' order."""
    m = inputs["A"]
    with mp.workdps(mp.dps + _span(m)):
        return _eig_sym(m)


def _eig_sym(m: Mat) -> dict:
    n = m.shape[0]
    E, Q = mp.eigsy(from_mat(m))
    order = sorted(range(n), key=lambda i: E[i])
    lam = [E[i] for i in order]
    cols = [[Q[r, i] for r in range(n)] for i in order]
    big = max(abs(x) for x in lam)
    if big == 0 or min(abs(x) for x in lam) < mpf(10) ** (-mp.dps // 2) * big:
        A = exact(m)
        zeros = sorted(sorted(range(n), key=lambda i: abs(lam[i]))[: n - rank(A)])
        for i in zeros:
            lam[i] = mpf(0)
        if len(zeros) > 1:
            for i, v in zip(zeros, _orthonormal(nullspace(A)), strict=True):
                cols[i] = v
    V = mp.matrix(n, n)
    for j, col in enumerate(cols):
        top = max(range(n), key=lambda r: (abs(col[r]), -r))
        sign = -1 if col[top] < 0 else 1
        for r in range(n):
            V[r, j] = sign * col[r]
    return {"lambda": lam, "V": to_mat(V), "kappa": [kappa(m)]}
