"""Generation-time cross-checks of the SO(2) and SE(2) records (docs/PHASE1.md section 4.3).

Each holds a property that fixes the output without the algorithm that produced it: the complex
series against `mp.expm`, the conjugation identity, J_l = Ad_Exp(tau) J_r, J J^-1 = I against the
series. A tolerance is 1e-100 of the size of the terms an entry sums (`_close`), not an absolute
one: at `theta:subnormal` the angle-dependent entries are 1e-311 next to entries of 1, and an
absolute 1e-100 would pass any of them. Every dense output also asserts the structure of a 3x3
Jacobian on the generated data: the first row is (1, 0, 0) and the rotation-translation block has
the form p I + q J, J the generator of rotations.
"""

from collections.abc import Callable

from mpmath import mp, mpf

from . import check, check_sen3, se2, sen3, so2, so3
from .check import CrossCheckError

_flat, _rows = check_sen3._entries, check_sen3._rows


def _close(what: str, got, want, size) -> None:
    """Every entry of `got` within 1e-100 of `size[i]` of `want[i]`: exactly equal where 0."""
    for g, w, s in zip(got, want, size, strict=True):
        if abs(g - w) > mpf(10) ** -check.DIGITS * s:
            raise CrossCheckError(f"{what}: {mp.nstr(g, 40)} vs {mp.nstr(w, 40)}")


def _abs(A):
    return [[abs(x) for x in row] for row in A]


def _product(what: str, A, B, want) -> None:
    """A B = want, entry by entry against the size of the terms of the sum."""
    size = [s + abs(w) for s, w in zip(_flat(sen3.mul(_abs(A), _abs(B))), _flat(want), strict=True)]
    _close(what, _flat(sen3.mul(A, B)), _flat(want), size)


def _structure(what: str, J) -> None:
    _close(f"{what}: first row", J[0], [1, 0, 0], [1, 0, 0])
    _close(f"{what}: diagonal of the block", [J[1][1]], [J[2][2]], [max(1, abs(J[1][1]))])
    _close(f"{what}: off-diagonal of the block", [J[1][2]], [-J[2][1]], [abs(J[2][1])])


def so2_exp(inputs: dict, out: dict) -> None:
    """z is the unit complex of the complex series (c to 1e-100, s to 1e-100 of itself)."""
    E = so2.exp_series(mpf(inputs["theta"]))
    _close("so2_exp", out["z"], [E.real, E.imag], [1, abs(E.imag)])
    _close("so2_exp: |z|", [sum(c * c for c in out["z"])], [1], [1])


def so2_log(inputs: dict, out: dict) -> None:
    """|theta| <= pi, and `mp.expm` of the hat matrix is z / |z|."""
    z, theta = so3.normalize([mpf(c) for c in inputs["z"]]), out["theta"]
    if abs(theta) > mp.pi + mpf(10) ** -check.DIGITS:
        raise CrossCheckError(f"so2_log: |theta| > pi: {mp.nstr(theta, 40)}")
    E = so2.expm(theta)
    _close("so2_log: Exp(theta) = z", [E[0][0], E[1][0]], z, [1, abs(z[1])])


def se2_exp(inputs: dict, out: dict) -> None:
    """z as `so2_exp`; t = V(theta) rho with V the series of the same block; Log(Exp tau) = tau."""
    tau = [mpf(c) for c in inputs["tau"]]
    so2_exp({"theta": inputs["tau"][0]}, {"z": out["z"]})
    V, rho = se2.translation_block(tau[0]), [[r] for r in tau[1:]]
    _close("se2_exp: t", out["t"], _flat(sen3.mul(V, rho)), _flat(sen3.mul(_abs(V), _abs(rho))))
    _close("se2_exp: log round trip", se2.log(out)["tau"], tau, [abs(c) for c in tau])


def se2_log(inputs: dict, out: dict) -> None:
    """`so2_log` on theta, and `mp.expm` of the whole hat matrix returns X."""
    tau, t = out["tau"], [mpf(c) for c in inputs["t"]]
    so2_log(inputs, {"theta": tau[0]})
    E = se2.expm(tau)
    _close(
        "se2_log: Exp(Log X) translation",
        [E[0][2], E[1][2]],
        t,
        [mp.sqrt(sum(c * c for c in t))] * 2,
    )


def se2_ad(inputs: dict, out: dict) -> None:
    """The block is R(z); X Exp(u) X^-1 = Exp(Ad u) for one u: both sides `mp.expm`, X inverted."""
    A, t = _rows(out["Ad"]), [mpf(c) for c in inputs["t"]]
    c, s = so3.normalize([mpf(x) for x in inputs["z"]])
    _structure("se2_ad", A)
    _close(
        "se2_ad: block is R", _flat([r[1:] for r in A[1:]]), [c, -s, s, c], [1, abs(s), abs(s), 1]
    )
    X = se2._affine([[c, -s], [s, c]], t)
    u = [mpf((-1) ** j * (j + 1)) / 4 for j in range(3)]
    lhs = sen3.mul(sen3.mul(X, se2.expm(u)), se2.inv(X))
    rhs = se2.expm([sum(a * x for a, x in zip(row, u, strict=True)) for row in A])
    check_sen3._agree("se2_ad: X Exp(u) X^-1 = Exp(Ad u)", _flat(lhs), _flat(rhs))


def jacobian(name: str) -> Callable[[dict, dict], None]:
    """J_l = Ad_Exp(tau) J_r, Exp by `mp.expm`, against the other series; the inverses: J J^-1 = I
    against the series."""
    right = name.startswith("jr")

    def check_(inputs: dict, out: dict) -> None:
        tau, J = [mpf(c) for c in inputs["tau"]], _rows(out["J"])
        _structure(f"se2_{name}", J)
        if name.endswith("inv"):
            return _product(f"se2_{name}", J, se2.series(tau, -1 if right else 1), sen3.eye(3))
        E = se2.expm(tau)
        Ad = se2.adjoint([r[:2] for r in E[:2]], [E[0][2], E[1][2]])
        other = se2.series(tau, 1 if right else -1)
        A, B, want = (Ad, J, other) if right else (Ad, other, J)
        _product(f"se2_{name}", A, B, want)

    return check_
