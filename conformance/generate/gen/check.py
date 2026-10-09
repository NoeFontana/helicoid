"""Generation-time cross-checks of records (docs/PHASE1.md section 4.3).

Every coefficient record's value and d_branch must agree to 100 digits with two computations that
share nothing with `mp.diff` on the definition: the committed rational Taylor series where the
branch variable is small, and a second formulation at 240 digits everywhere: the exact series of
`series.exact` (`k`, `a`...`e`), a closed form and its calculus derivative (`r`), and for `alpha`
(0062), which has no committed series, its own alternating series.

Every SO(3) record is checked to 100 digits by a property that pins its output down without the
algorithm that produced it: `mp.expm` of the hat matrix, uniqueness of the polar decomposition,
the sandwich product, the Jacobian identities of docs/NUMERICS.md section 1.
"""

from collections.abc import Callable
from functools import cache, lru_cache
from math import ceil

from mpmath import mp, mpf

from . import coeff, series, so3
from .precision import DPS

DIGITS = 100
# The committed series truncate after `series.TERMS` terms; their relative truncation error at
# x <= SMALL is ~ TERMS * SMALL^(TERMS - 1) (`r`, the slowest), kept 20 digits under DIGITS.
SMALL = mpf(10) ** -ceil((DIGITS + 20) / (series.TERMS - 1))


class CrossCheckError(Exception):
    pass


@cache
def _mp(name: str, exact: bool) -> tuple:
    """The series as mpf at 240 digits: `series.exact`'s 232 terms, or the committed 16."""
    with mp.workdps(2 * DPS):
        rows = series.exact()[name] if exact else series.taylor(name)
        return tuple(mpf(q.numerator) / q.denominator for q in rows)


def _sum(rows: tuple, x) -> tuple:
    """(sum s_j x^j, sum j s_j x^(j-1)) by Horner."""
    value = deriv = mpf(0)
    for s in reversed(rows):
        deriv = deriv * x + value
        value = value * x + s
    return value, deriv


def _compare(name: str, inputs: dict, out: dict, want: tuple, how: str) -> None:
    for key, w in zip(("value", "d_branch"), want, strict=True):
        if abs(out[key] - w) > abs(w) * mpf(10) ** -DIGITS:
            raise CrossCheckError(f"{name} {key} vs {how} {inputs}: {out[key]} vs {w}")


def _r_closed(n, w) -> tuple:
    """r and dr/d(n^2) at fixed w, n > 0, by calculus on 2 atan2(n, w) / n."""

    def angle():
        return mp.atan(n / w) if w > 0 else mp.pi / 2

    def derivative(_):  # cancels ~ 1 / n^2: `stable` gives it the guard digits
        return w / ((w * w + n * n) * n * n) - angle() / n**3

    return 2 * angle() / n, coeff.stable(derivative, n * n)


@cache
def _alpha_series() -> tuple:
    """sin(theta) / theta = sum_j (-x)^j / (2j + 1)!, x = theta^2, to 232 terms at 240 digits
    (0062): alpha has no committed series, so its second formulation is this one everywhere."""
    with mp.workdps(2 * DPS):
        return tuple(mpf(-1) ** j / mp.factorial(2 * j + 1) for j in range(232))


def _check_alpha(inputs: dict, out: dict) -> None:
    _compare("alpha", inputs, out, _sum(_alpha_series(), mpf(inputs["theta"]) ** 2), "series")


def _check_r(inputs: dict, out: dict) -> None:
    n, w = mpf(inputs["n"]), mpf(inputs["w"])
    if w > 0 and (n / w) ** 2 <= SMALL:  # value (2/w) S(s), d/d(n^2) (2/w^3) S'(s), s = n^2/w^2
        v, d = _sum(_mp("r", False), (n / w) ** 2)
        _compare("r", inputs, out, (2 / w * v, 2 / w**3 * d), "Taylor")
    if n > 0:
        _compare("r", inputs, out, _r_closed(n, w), "closed form")


def coefficient(name: str) -> Callable[[dict, dict], None]:
    """The check of coefficient `name` for `FunctionSpec.check`: (inputs, out) -> None."""

    def check(inputs: dict, out: dict) -> None:
        with mp.workdps(2 * DPS):
            if name == "r":
                return _check_r(inputs, out)
            if name == "alpha":
                return _check_alpha(inputs, out)
            x = mpf(inputs["theta"]) ** 2
            if x <= SMALL:
                _compare(name, inputs, out, _sum(_mp(name, False), x), "Taylor")
            rows = _mp(name, True)[: 96 if x <= 1 else None]  # 96 terms: below 1e-150 at x <= 1
            _compare(name, inputs, out, _sum(rows, x), "exact series")

    return check


def _agree(what: str, got, want, scale=1) -> None:
    """Every entry within 1e-100 of `want`, times `scale`: 1 for O(1) entries, or the size of the
    tiny part being checked (1e-311 at `theta:subnormal`)."""
    tol = mpf(10) ** -DIGITS * scale
    for g, w in zip(got, want, strict=True):
        if abs(g - w) > tol:
            raise CrossCheckError(f"{what}: {mp.nstr(g, 40)} vs {mp.nstr(w, 40)}")


def _offdiag(A) -> mpf:
    """The size of a rotation matrix's tiny part: the largest off-diagonal entry."""
    return so3.maxabs(A[i][j] for i in range(3) for j in range(3) if i != j)


def _agree_mat(what: str, got, want, tiny) -> None:
    """3x3 matrices: the diagonal to 1e-100, the off-diagonal entries to 1e-100 of `tiny`, the size
    of the off-diagonal part of the series (1e-310 at `theta:subnormal`, where the diagonal is 1
    to 1e-620 and the information sits off it)."""
    diag = [i for i in range(9) if i % 4 == 0]
    g, w = so3.entries(got), so3.entries(want)
    _agree(what, [g[i] for i in diag], [w[i] for i in diag])
    off = [i for i in range(9) if i % 4]
    _agree(f"{what}: off-diagonal", [g[i] for i in off], [w[i] for i in off], tiny)


def _det3(A):
    return (
        A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
        - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
        + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])
    )


def _rotation(what: str, q: list) -> list:
    """R(q) of a unit q: orthogonal with determinant +1."""
    R = so3.rot(q)
    _agree(f"{what}: |q|", [so3.norm(q)], [1])
    _agree(f"{what}: R^T R", so3.entries(so3.mm(so3.mt(R), R)), so3.entries(so3.identity()))
    _agree(f"{what}: det R", [_det3(R)], [1])
    return R


@lru_cache(maxsize=2)  # a record and its negative share one phi
def _expm_of(phi: tuple, dps: int):
    E = mp.expm(mp.matrix(so3.hat(phi)))
    return [[E[i, j] for j in range(3)] for i in range(3)]


def _expm(phi: list):
    return _expm_of(tuple(phi), mp.dps)


def _matches_expm(what: str, phi: list, q: list) -> None:
    """R(q) is `mp.expm` of the hat matrix: entries to 1e-100, and the skew part, which carries
    the vector part of q, to 1e-100 of the off-diagonal size (1e-311 at `theta:subnormal`)."""
    R, E = _rotation(what, q), _expm(phi)
    _agree(what, so3.entries(R), so3.entries(E))
    _agree(f"{what}: skew part", so3.vee_skew(R), so3.vee_skew(E), _offdiag(E))


def so3_exp(inputs: dict, out: dict) -> None:
    _matches_expm("so3_exp", [mpf(c) for c in inputs["phi"]], out["q"])


def so3_log(inputs: dict, out: dict) -> None:
    """Exp(phi) is q (the residual, w > 0 side, to 1e-100 of |vec q|), and |phi| <= pi, and
    `mp.expm` agrees."""
    q, phi = so3.canonical([mpf(c) for c in inputs["q"]]), out["phi"]
    r = so3.qmul(so3.qconj(so3.exp_series(phi)), q)
    _agree("so3_log: exp residual", [2 * c for c in r[1:]], [0] * 3, so3.maxabs(q[1:]))
    _agree("so3_log: exp residual w", [r[0]], [1])
    if so3.norm(phi) > mp.pi + mpf(10) ** -DIGITS:
        raise CrossCheckError(f"so3_log: |phi| > pi: {mp.nstr(so3.norm(phi), 40)}")
    _matches_expm("so3_log", phi, q)


def so3_act(inputs: dict, out: dict) -> None:
    """R(q/|q|) p is the sandwich q p q^* of the normalized quaternion (Hamilton, active)."""
    q = so3.normalize([mpf(c) for c in inputs["q"]])
    _rotation("so3_act", q)
    p = [mpf(0)] + [mpf(c) for c in inputs["p"]]
    _agree("so3_act", out["Rp"], so3.qmul(so3.qmul(q, p), so3.qconj(q))[1:])


def so3_from_matrix(inputs: dict, out: dict) -> None:
    """A = R(q) H with H symmetric positive definite: R(q) is the polar factor of A, so the
    nearest rotation in Frobenius norm."""
    A, q = so3.from_mat(inputs["R"]), out["q"]
    H = so3.mm(so3.mt(_rotation("so3_from_matrix", q)), A)
    _agree("so3_from_matrix: H skew", so3.vee_skew(H), [0] * 3, _offdiag(A))
    if not (H[0][0] > 0 and H[0][0] * H[1][1] - H[0][1] * H[1][0] > 0 and _det3(H) > 0):
        raise CrossCheckError("so3_from_matrix: H is not positive definite")


def so3_gamma2(inputs: dict, out: dict) -> None:
    """GG.2(b) I + Gamma_2 W = Gamma_1 against `so3_jl`'s series, which pins Gamma_2 off the kernel
    of W, and GG.2(d) Gamma_2 phi = phi / 2, which pins it on that kernel."""
    phi, G = [mpf(c) for c in inputs["phi"]], so3.from_mat(out["G"])
    shifted = [
        [(1 if i == j else 0) + w for j, w in enumerate(row)]
        for i, row in enumerate(so3.mm(G, so3.hat(phi)))
    ]
    jl = so3.jl(phi)
    _agree_mat("so3_gamma2: I + Gamma_2 W", shifted, jl, _offdiag(jl))
    scale = so3.maxabs(phi) or mpf(1)
    G_phi = [sum(G[i][k] * phi[k] for k in range(3)) for i in range(3)]
    _agree("so3_gamma2: Gamma_2 phi", G_phi, [c / 2 for c in phi], scale)


def so3_jacobian(name: str) -> Callable[[dict, dict], None]:
    """J_l = R J_r (NUMERICS section 1) for `jr`/`jl`; J^-1 J = I, against the series, for the
    inverses. The off-diagonal part is held to 1e-100 of the series' own: it is all there is to
    check at `theta:subnormal`."""

    def check(inputs: dict, out: dict) -> None:
        phi, J = [mpf(c) for c in inputs["phi"]], so3.from_mat(out["J"])
        right = name.startswith("jr")
        if name.endswith("inv"):
            series = (so3.jr if right else so3.jl)(phi)
            got, want, tiny = so3.mm(J, series), so3.identity(), _offdiag(series)
        else:
            R = so3.rot(so3.exp_series(phi))
            got, want = so3.mm(R if right else so3.mt(R), J), (so3.jl if right else so3.jr)(phi)
            tiny = _offdiag(want)
        _agree_mat(f"so3_{name}", got, want, tiny)

    return check
