"""SE_N(3) ids (docs/PHASE1.md section 4.3, N = 1, 2, 3), computed from their definitions.

`sen3_exp` is `mp.expm` of the (3+N)-square hat matrix (its rotation part as the quaternion power
series of `gen.so3`), `sen3_log` the inverse of that on theta in [0, pi] (Newton on the quaternion
series for phi, then rho_i from the translation block of the same exponential), `sen3_ad` the
images of the basis under sigma -> (X sigma^ X^-1)^vee, `sen3_jr` and `sen3_jl` the defining series
of ad_tau, dense, and their inverses `mp.inverse` of those (the blocks a dual matrix has no entry
in exactly 0). No block formula of docs/NUMERICS.md section 5 appears: the generator never
evaluates the Q block or a dual-matrix inverse.

Tangents are flat and rotation-first, [phi; rho_1; ...; rho_N]; X is (q, x), x = [x_1; ...; x_N].
Dense matrices are rows of mpf here and column-major `fmt.Mat` in records, in the tangent's order.
"""

from mpmath import mp, mpf

from . import so3
from .fmt import Mat
from .precision import DPS, to_f64
from .rng import SEED, stream, unit_vector_s2
from .strata import SEN3_SAMPLES, Stratum


def hat(tau, n: int):
    """The (3+N)-square hat matrix [[phi^, rho_1 ... rho_N], [0, 0]]."""
    m, o = 3 + n, mpf(0)
    H = [[o] * m for _ in range(m)]
    for i, row in enumerate(so3.hat(tau[:3])):
        H[i][:3] = row
    for k in range(n):
        for r in range(3):
            H[r][3 + k] = tau[3 + 3 * k + r]
    return H


def vee(M, n: int) -> list:
    """The tangent of a (3+N)-square matrix of the Lie algebra: skew part of phi^, columns rho_i."""
    top = so3.vee_skew([row[:3] for row in M[:3]])
    return top + [M[r][3 + k] for k in range(n) for r in range(3)]


def eye(m: int):
    return [[mpf(int(i == j)) for j in range(m)] for i in range(m)]


def mul(A, B):
    """Dense product, skipping exact zeros of A: every entry is the same sum as without."""
    m = len(B[0])
    out = []
    for row in A:
        acc = [mpf(0)] * m
        for k, a in enumerate(row):
            if a != 0:
                for j, b in enumerate(B[k]):
                    if b != 0:
                        acc[j] += a * b
        out.append(acc)
    return out


def sub(A, B):
    return [[a - b for a, b in zip(ra, rb, strict=True)] for ra, rb in zip(A, B, strict=True)]


def expm(tau, n: int):
    return [list(row) for row in mp.expm(mp.matrix(hat(tau, n))).tolist()]


def basis(m: int, j: int) -> list:
    return [mpf(int(i == j)) for i in range(m)]


def ad(tau, n: int):
    """ad_tau: the images of the basis under sigma -> [tau^, sigma^]^vee, as columns."""
    m, T = 3 + 3 * n, hat(tau, n)
    cols = []
    for j in range(m):
        S = hat(basis(m, j), n)
        cols.append(vee(sub(mul(T, S), mul(S, T)), n))
    return [[cols[j][i] for j in range(m)] for i in range(m)]


def adjoint(R, xs: list, n: int):
    """Ad_X: the images of the basis under sigma -> (X sigma^ X^-1)^vee, X = [[R, x], [0, I]]."""
    m, Rinv, o = 3 + 3 * n, so3.inverse(R), mpf(0)
    X, Xinv = eye(3 + n), eye(3 + n)
    for i in range(3):
        for j in range(3):
            X[i][j], Xinv[i][j] = R[i][j], Rinv[i][j]
    for k in range(n):
        for r in range(3):
            X[r][3 + k] = xs[3 * k + r]
            Xinv[r][3 + k] = -sum((Rinv[r][c] * xs[3 * k + c] for c in range(3)), o)
    cols = [vee(mul(mul(X, hat(basis(m, j), n)), Xinv), n) for j in range(m)]
    return [[cols[j][i] for j in range(m)] for i in range(m)]


def matrix_series(A):
    """sum A^k / (k+1)!, stopping like `so3.matrix_series`: two terms in a row below `so3._eps()`
    of the partial sum, entry by entry."""
    eps = so3._eps()
    term = total = eye(len(A))
    quiet = 0
    for k in range(1, so3.MAX_ITER):
        term = [[x / (k + 1) for x in row] for row in mul(term, A)]
        total = [[s + t for s, t in zip(rs, rt, strict=True)] for rs, rt in zip(total, term)]
        small = all(abs(t) <= eps * abs(s) for rt, rs in zip(term, total) for t, s in zip(rt, rs))
        quiet = quiet + 1 if small else 0
        if quiet == 2:
            return total
    raise ArithmeticError("Jacobian series does not converge")


def series(tau, n: int, sign: int):
    """sum (sign * ad_tau)^k / (k+1)!. sign -1 is J_r, +1 J_l."""
    return matrix_series([[sign * a for a in row] for row in ad(tau, n)])


def structural_zero(bi: int, bj: int) -> bool:
    """Block (bi, bj) of a dual matrix has no entry: off the diagonal and the first block column."""
    return bi != bj and not (bj == 0 and bi > 0)


def inverse(J):
    """`mp.inverse`, with the blocks of `structural_zero` exactly 0: the LU leaves rounding residue
    there, which the 150-digit recheck cannot compare. Every other entry is the LU's, however small
    (the recheck agrees on all of them, down to 1e-311)."""
    B = mp.inverse(mp.matrix(J))
    return [
        [mpf(0) if structural_zero(i // 3, j // 3) else B[i, j] for j in range(B.cols)]
        for i in range(B.rows)
    ]


def to_mat(A) -> Mat:
    m = len(A)
    return Mat((m, m), [A[i][j] for j in range(m) for i in range(m)])


def chunks(v, n: int) -> list:
    return [v[3 * k : 3 * k + 3] for k in range(n)]


def exp(n: int):
    def evaluate(inputs: dict) -> dict:
        tau = [mpf(c) for c in inputs["tau"]]
        E = expm(tau, n)
        q = so3.exp_series(tau[:3])
        x = [E[r][3 + k] for k in range(n) for r in range(3)]
        return {"q": q[:1] + so3.parallel(q[1:], tau[:3]), "x": x}

    return evaluate


def log(n: int):
    def evaluate(inputs: dict) -> dict:
        phi = so3.log({"q": inputs["q"]})["phi"]
        Jinv = so3.inverse(so3.jl(phi))
        rho = [
            sum((Jinv[r][c] * mpf(x[c]) for c in range(3)), mpf(0))
            for x in chunks(inputs["x"], n)
            for r in range(3)
        ]
        return {"tau": list(phi) + rho}

    return evaluate


def ad_of(n: int):
    def evaluate(inputs: dict) -> dict:
        R = so3.rot(so3.normalize([mpf(c) for c in inputs["q"]]))
        return {"Ad": to_mat(adjoint(R, [mpf(c) for c in inputs["x"]], n))}

    return evaluate


def jacobian(n: int, name: str):
    """`sen3_jr`, `sen3_jl` and their inverses."""

    def evaluate(inputs: dict) -> dict:
        J = series([mpf(c) for c in inputs["tau"]], n, -1 if name.startswith("jr") else 1)
        return {"J": to_mat(inverse(J) if name.endswith("inv") else J)}

    return evaluate


# Inputs: exact binary64, one dict per record. A stratum's theta (fixed, or random for
# `theta:1e-k`) and axis give phi or q as in `gen.so3`, a quaternion stratum gives q; each x_i or
# rho_i is a random direction at the stratum's translation scale (1 for a `theta:*` or `q:*`
# stratum), its own stream per index i so that the N = 2 records extend the N = 1 ones.


def _vectors(stratum: Stratum, n: int, count: int) -> list[list[float]]:
    streams = [stream(SEED, stratum.name, f"rho{k + 1}") for k in range(n)]
    with mp.workdps(DPS):
        scale = mpf(10) ** (stratum.rho_exp or 0)
        return [
            [to_f64(scale * mpf(c)) for r in streams for c in unit_vector_s2(r)]
            for _ in range(count)
        ]


def _samples(stratum: Stratum, count: int = SEN3_SAMPLES) -> list:
    thetas = stratum.thetas()
    axes = stratum.axes(count) or [(1.0, 0.0, 0.0)] * count
    return [(thetas[i % len(thetas)], axes[i]) for i in range(count)]


def _quaternions(stratum: Stratum) -> list[list[float]]:
    """The q of a stratum: a quaternion stratum's own, or the unit quaternion of (theta, axis)."""
    if stratum.draw_quats is not None:
        return [list(q) for q in stratum.quaternions()[:SEN3_SAMPLES]]
    return [so3.quat_of(t, a) for t, a in _samples(stratum)]


def tau_inputs(n: int):
    def inputs(stratum: Stratum) -> list[dict]:
        rho = _vectors(stratum, n, SEN3_SAMPLES)
        return [
            {"tau": so3.phi_of(t, a) + r} for (t, a), r in zip(_samples(stratum), rho, strict=True)
        ]

    return inputs


def x_inputs(n: int, both_signs: bool):
    """(q, x): q the unit quaternion of (theta, axis), or a quaternion stratum's. `log` holds every
    q and its negative, the half a `Log` without the flip gets wrong; w = +0 has no negative."""

    def inputs(stratum: Stratum) -> list[dict]:
        rho = _vectors(stratum, n, SEN3_SAMPLES)
        out = []
        for q, x in zip(_quaternions(stratum), rho, strict=True):
            out.append({"q": q, "x": x})
            if both_signs and q[0] != 0:
                out.append({"q": [-c for c in q], "x": x})
        return out

    return inputs
