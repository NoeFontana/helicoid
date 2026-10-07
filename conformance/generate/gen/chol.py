"""`chol_n{3,6}` and `chol_solve_n{3,6}` (docs/decisions/0056).

`chol`: `valid` is decided exactly, by an LDL^T of the stored A in rationals; `L` is `mp.cholesky`
of the stored A, or zeros where `valid` is 0. `chol_solve`: `mp.lu_solve` of the stored A and b.

A = Q D Q^T at the working precision, Q a product of n Householder reflections from directions
uniform on S^(n-1), rounded once to the stratum's precision (the upper triangle, mirrored).
"""

from mpmath import mp, mpf

from .dense import exact, from_mat, ldlt_pivots, round_symmetric, to_mat, unit_vector
from .strata import Draw, drawn_strata

N_CHOL = 32
COND_K = {False: (4, 8, 12), True: (2, 4, 6)}
DIAG_SCALE = {False: 100, True: 15}  # S log-uniform in [10^-s, 10^s]
INDEFINITE = 100  # the negative eigenvalue is -max D / INDEFINITE


def _householder(rng, n: int):
    Q = mp.eye(n)
    for _ in range(n):
        v = mp.matrix(unit_vector(rng, n))
        Q = Q * (mp.eye(n) - 2 * v * v.T)
    return Q


def _sandwich(Q, D):
    return Q * mp.diag(D) * Q.T


def _log10_uniform(rng, lo: int, hi: int):
    return mpf(10) ** (lo + (hi - lo) * mpf(rng.uniform()))


def _spd_draws(draw: Draw, n: int, count: int, of: str | None = None) -> list:
    rq, rd = draw.stream("householder", of), draw.stream("spectrum", of)
    out = []
    for _ in range(count):
        Q = _householder(rq, n)
        out.append(_sandwich(Q, [_log10_uniform(rd, 0, 1) for _ in range(n)]))
    return out


def _spd(n: int):
    def make(draw: Draw, count: int) -> list:
        return [round_symmetric(A, draw.round) for A in _spd_draws(draw, n, count)]

    return make


def _cond(n: int, k: int):
    """D log-spaced from 1 to 10^-k, Q random."""

    def make(draw: Draw, count: int) -> list:
        D = [mpf(10) ** (-mpf(k) * i / (n - 1)) for i in range(n)]
        rq = draw.stream("householder")
        return [round_symmetric(_sandwich(_householder(rq, n), D), draw.round) for _ in range(count)]

    return make


def _diag_scale(n: int):
    """S A S, A the `chol:spd` record of the same index (from its stream), S from this stratum's."""

    def make(draw: Draw, count: int) -> list:
        rs, s = draw.stream("scale"), DIAG_SCALE[draw.f32]
        out = []
        for A in _spd_draws(draw, n, count, of="chol:spd"):
            S = mp.diag([_log10_uniform(rs, -s, s) for _ in range(n)])
            out.append(round_symmetric(S * A * S, draw.round))
        return out

    return make


def _indefinite(n: int):
    def make(draw: Draw, count: int) -> list:
        rq, rd = draw.stream("householder"), draw.stream("spectrum")
        out = []
        for _ in range(count):
            Q = _householder(rq, n)
            D = [_log10_uniform(rd, 0, 1) for _ in range(n)]
            D[rd.next_u64() % n] = -max(D) / INDEFINITE
            out.append(round_symmetric(_sandwich(Q, D), draw.round))
        return out

    return make


def _with_b(n: int, make):
    """`chol_solve`'s records: the matrix stratum's A, exactly positive definite, and b uniform on
    the unit sphere from this stratum's stream."""

    def records(draw: Draw, count: int) -> list[dict]:
        rb = draw.stream("b")
        out = []
        for A in make(draw, count):
            if ldlt_pivots(exact(A)) is None:
                raise ArithmeticError(f"chol_solve_n{n} {draw.name}: A is not positive definite")
            out.append({"A": A, "b": [draw.round(c) for c in unit_vector(rb, n)]})
        return out

    return records


def _matrix_strata(n: int, f32: bool, solve: bool) -> list[tuple[str, object, bool]]:
    """(name, make, shared) in the catalogue's order."""
    strata = [
        ("chol:spd", _spd(n), True),
        *((f"chol:cond-1e-{k}", _cond(n, k), k in COND_K[False]) for k in COND_K[f32]),
        ("chol:diag-scale", _diag_scale(n), False),
    ]
    if solve:
        return [(name, _with_b(n, make), shared) for name, make, shared in strata]
    strata.append(("chol:indefinite", _indefinite(n), True))
    return [(name, lambda d, c, m=make: [{"A": A} for A in m(d, c)], s) for name, make, s in strata]


def strata(n: int, solve: bool):
    return drawn_strata(
        [(name, make, N_CHOL) for name, make, _ in _matrix_strata(n, False, solve)],
        [(name, make, N_CHOL, shared) for name, make, shared in _matrix_strata(n, True, solve)],
    )


# --- the references --------------------------------------------------------------------------


def chol(inputs: dict) -> dict:
    """`tol=0`: `mp.cholesky` refuses a pivot below an absolute eps, which every pivot of
    `chol:diag-scale` (down to 1e-200) is; positive definiteness is already decided exactly."""
    n = inputs["A"].shape[0]
    valid = ldlt_pivots(exact(inputs["A"])) is not None
    L = mp.cholesky(from_mat(inputs["A"]), tol=0) if valid else mp.zeros(n, n)
    return {"valid": [mpf(int(valid))], "L": to_mat(L)}


def equilibrator(A):
    """diag(|A_ii|^-1/2) (1 where A_ii = 0): D A D has a unit diagonal and A's inertia."""
    return mp.diag([1 / mp.sqrt(abs(A[i, i])) if A[i, i] else mpf(1) for i in range(A.rows)])


def chol_solve(inputs: dict) -> dict:
    """`mp.lu_solve` of the stored A and b, equilibrated: x = D (D A D)^-1 D b. mpmath's LU pivots
    against an absolute ||A||_1 eps and loses cond(A) eps, 1e400 eps at `chol:diag-scale`; D A D
    is within n of the best diagonal scaling (van der Sluis), so there within n of `chol:spd`."""
    A = from_mat(inputs["A"])
    D = equilibrator(A)
    y = mp.lu_solve(D * A * D, D * mp.matrix([mpf(c) for c in inputs["b"]]))
    x = D * y
    return {"x": [x[i] for i in range(x.rows)]}
