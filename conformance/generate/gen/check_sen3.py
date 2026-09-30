"""Generation-time cross-checks of the SE_N(3) records (docs/PHASE1.md section 4.3), to 100 digits.

Each holds a property that fixes the output without the algorithm that produced it: `mp.expm` of the
hat matrix (a second implementation of the exponential), the conjugation identity, J_l =
Ad_Exp(tau) J_r, J J^-1 = I against the series. Every dense output also asserts the dual-matrix
structure of docs/NUMERICS.md section 2.2 on the generated data: the blocks other than the
diagonal and the first block column are exactly 0, and the diagonal blocks are one block.
Tolerances are absolute, 1e-100 of max(1, the largest entry): the metric of docs/NUMERICS.md
section 11 floors a Jacobian's norm at 1, and an entry of 1e-310 next to an entry of 1e4 carries no
information in it.
"""

from collections.abc import Callable

from mpmath import mp, mpf

from . import check, sen3, so3
from .check import CrossCheckError


def _entries(A) -> list:
    return [x for row in A for x in row]


def _agree(what: str, got, want) -> None:
    got, want = list(got), list(want)
    check._agree(what, got, want, max(1, so3.maxabs(want), so3.maxabs(got)))


def _rows(out_mat) -> list:
    m = out_mat.shape[0]
    return [[out_mat.data[j * m + i] for j in range(m)] for i in range(m)]


def _dual(what: str, J, n: int) -> None:
    """`J` is [[A, 0], [B_i, A]] blockwise: exactly 0 off the diagonal and the first block column
    (a dual matrix has no other entry), and the diagonal blocks all equal."""
    for bi in range(n + 1):
        for bj in range(n + 1):
            block = [J[3 * bi + r][3 * bj + c] for r in range(3) for c in range(3)]
            if sen3.structural_zero(bi, bj) and any(x != 0 for x in block):
                raise CrossCheckError(f"{what}: block ({bi}, {bj}) is not 0")
            if bi == bj:
                _agree(f"{what}: diagonal block {bi}", block, _entries([r[:3] for r in J[:3]]))


def _translation(E, n: int) -> list:
    return [E[r][3 + k] for k in range(n) for r in range(3)]


def exp(n: int) -> Callable[[dict, dict], None]:
    """R(q) is `mp.expm` of the 3x3 hat; x_i = J_l(phi) rho_i, the series of the same block; and
    Log(Exp(tau)) = tau."""

    def check_(inputs: dict, out: dict) -> None:
        tau = [mpf(c) for c in inputs["tau"]]
        check.so3_exp({"phi": inputs["tau"][:3]}, {"q": out["q"]})
        Jl = so3.jl(tau[:3])
        want = [
            sum((Jl[r][c] * t[c] for c in range(3)), mpf(0))
            for t in sen3.chunks(tau[3:], n)
            for r in range(3)
        ]
        _agree("sen3_exp: x", out["x"], want)
        _agree("sen3_exp: log round trip", sen3.log(n)(out)["tau"], tau)

    return check_


def log(n: int) -> Callable[[dict, dict], None]:
    """|phi| <= pi with Exp(phi) = q (`so3_log`'s check), and `mp.expm` of the whole hat matrix
    returns X: the rotation of q/|q| and x."""

    def check_(inputs: dict, out: dict) -> None:
        tau, x = out["tau"], [mpf(c) for c in inputs["x"]]
        check.so3_log({"q": inputs["q"]}, {"phi": tau[:3]})
        E = sen3.expm(tau, n)
        R = so3.rot(so3.canonical([mpf(c) for c in inputs["q"]]))
        _agree("sen3_log: Exp(Log X) rotation", _entries([r[:3] for r in E[:3]]), _entries(R))
        _agree("sen3_log: Exp(Log X) translation", _translation(E, n), x)

    return check_


def ad(n: int) -> Callable[[dict, dict], None]:
    """X Exp(t) X^-1 = Exp(Ad_X t) for one t: both sides `mp.expm`, X inverted whole."""

    def check_(inputs: dict, out: dict) -> None:
        A, m, x = _rows(out["Ad"]), 3 + 3 * n, [mpf(c) for c in inputs["x"]]
        _dual("sen3_ad", A, n)
        R = so3.rot(so3.normalize([mpf(c) for c in inputs["q"]]))
        X = sen3.eye(3 + n)
        for i in range(3):
            X[i][:3] = R[i]
            for k in range(n):
                X[i][3 + k] = x[3 * k + i]
        Xinv = mp.inverse(mp.matrix(X))
        t = [mpf((-1) ** j * (j + 1)) / (m + 1) for j in range(m)]
        lhs = sen3.mul(sen3.mul(X, sen3.expm(t, n)), Xinv.tolist())
        rhs = sen3.expm([sum(a * c for a, c in zip(row, t, strict=True)) for row in A], n)
        _agree("sen3_ad: X Exp(t) X^-1 = Exp(Ad t)", _entries(lhs), _entries(rhs))

    return check_


def jacobian(n: int, name: str) -> Callable[[dict, dict], None]:
    """J_l = Ad_Exp(tau) J_r, Exp by `mp.expm`, against the other series; the inverses: J J^-1 = I
    against the series."""
    right, m = name.startswith("jr"), 3 + 3 * n

    def check_(inputs: dict, out: dict) -> None:
        tau, J = [mpf(c) for c in inputs["tau"]], _rows(out["J"])
        _dual(f"sen3_{name}", J, n)
        if name.endswith("inv"):
            S = sen3.series(tau, n, -1 if right else 1)
            _agree(f"sen3_{name}", _entries(sen3.mul(J, S)), _entries(sen3.eye(m)))
            return
        E = sen3.expm(tau, n)
        Ad = sen3.adjoint([r[:3] for r in E[:3]], _translation(E, n), n)
        other = sen3.series(tau, n, 1 if right else -1)
        got, want = (sen3.mul(Ad, J), other) if right else (sen3.mul(Ad, other), J)
        _agree(f"sen3_{name}", _entries(got), _entries(want))

    return check_
