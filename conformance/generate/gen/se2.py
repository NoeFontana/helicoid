"""SE(2) ids (docs/PHASE1.md section 4.3, `se2_*`), computed from their definitions.

`se2_exp` is `mp.expm` of the 3x3 hat matrix [[0, -theta, rho_x], [theta, 0, rho_y], [0, 0, 0]],
`se2_log` the inverse of it on theta in (-pi, pi] (`gen.so2`'s angle, then rho solving V(theta)
rho = t, V the series of the same exponential's translation block), `se2_ad` the images of the
basis under sigma -> (X sigma^ X^-1)^vee, `se2_jr` and `se2_jl` the defining series of ad_tau,
dense, and their inverses `mp.inverse` of those. No closed form appears (docs/NUMERICS.md section
6 has none for the Jacobians, and its V, alpha and beta are not evaluated).

Tangents are [theta; rho_x; rho_y], rotation first; X is (z, t), z = [c, s]. Dense matrices are
rows of mpf here and column-major `fmt.Mat` in records, in the tangent's order.
"""

from mpmath import mp, mpf

from . import sen3, so2, so3
from .precision import DPS, to_f64
from .rng import SEED, stream, unit_vector_s1
from .strata import SE2_SAMPLES, Stratum

BASIS = [sen3.basis(3, j) for j in range(3)]


def hat(tau):
    theta, x, y = tau
    o = mpf(0)
    return [[o, -theta, x], [theta, o, y], [o, o, o]]


def vee(M) -> list:
    """The tangent of a 3x3 matrix of the Lie algebra: skew part of the rotation block, column 2."""
    return [(M[1][0] - M[0][1]) / 2, M[0][2], M[1][2]]


def inv(A):
    B = mp.inverse(mp.matrix(A))
    return [[B[i, j] for j in range(B.cols)] for i in range(B.rows)]


def expm(tau):
    return [list(row) for row in mp.expm(mp.matrix(hat(tau))).tolist()]


def translation_block(theta):
    """V(theta) = sum (theta J)^n / (n+1)!, J the generator of rotations: the translation block of
    Exp, t = V rho, as a series of the same exponential."""
    return sen3.matrix_series(so2.hat(theta))


def ad(tau):
    """ad_tau: the images of the basis under sigma -> [tau^, sigma^]^vee, as columns."""
    T = hat(tau)
    cols = [vee(sen3.sub(sen3.mul(T, S), sen3.mul(S, T))) for S in map(hat, BASIS)]
    return [[cols[j][i] for j in range(3)] for i in range(3)]


def _affine(R, t):
    return [[R[0][0], R[0][1], t[0]], [R[1][0], R[1][1], t[1]], [mpf(0), mpf(0), mpf(1)]]


def adjoint(R, t):
    """Ad_X: the images of the basis under sigma -> (X sigma^ X^-1)^vee, X = [[R, t], [0, 1]]."""
    Rinv = inv(R)
    X = _affine(R, t)
    Xinv = _affine(Rinv, [-(Rinv[i][0] * t[0] + Rinv[i][1] * t[1]) for i in range(2)])
    cols = [vee(sen3.mul(sen3.mul(X, hat(b)), Xinv)) for b in BASIS]
    return [[cols[j][i] for j in range(3)] for i in range(3)]


def series(tau, sign: int):
    """sum (sign * ad_tau)^k / (k+1)!. sign -1 is J_r, +1 J_l."""
    return sen3.matrix_series([[sign * a for a in row] for row in ad(tau)])


def exp(inputs: dict) -> dict:
    E = expm([mpf(c) for c in inputs["tau"]])
    return {"z": [E[0][0], E[1][0]], "t": [E[0][2], E[1][2]]}


def log(inputs: dict) -> dict:
    theta = so2.angle([mpf(c) for c in inputs["z"]])
    Vinv = inv(translation_block(theta))
    t = [mpf(c) for c in inputs["t"]]
    return {"tau": [theta] + [sum((Vinv[i][k] * t[k] for k in range(2)), mpf(0)) for i in range(2)]}


def ad_of(inputs: dict) -> dict:
    R = so3.normalize([mpf(c) for c in inputs["z"]])
    rot = [[R[0], -R[1]], [R[1], R[0]]]
    return {"Ad": sen3.to_mat(adjoint(rot, [mpf(c) for c in inputs["t"]]))}


def jacobian(name: str):
    """`se2_jr`, `se2_jl` and their inverses."""

    def evaluate(inputs: dict) -> dict:
        J = series([mpf(c) for c in inputs["tau"]], -1 if name.startswith("jr") else 1)
        return {"J": sen3.to_mat(inv(J) if name.endswith("inv") else J)}

    return evaluate


# Inputs: exact binary64, one dict per record. Sample i of a stratum takes the i-th theta of its
# stream (a fixed theta repeats) with the sign (-1)^i, and a translation, a random direction at
# the stratum's scale (1 for a `theta:*` stratum) from its own stream; se2_log and se2_ad take
# (z, t) with z the unit complex of theta rounded to binary64, t the translation.


def _samples(stratum: Stratum) -> list[tuple[float, list[float]]]:
    thetas, rng = stratum.thetas(), stream(SEED, stratum.name, "rho")
    with mp.workdps(DPS):
        scale = mpf(10) ** (stratum.rho_exp or 0)
        rho = [[to_f64(scale * mpf(c)) for c in unit_vector_s1(rng)] for _ in range(SE2_SAMPLES)]
    out = []
    for i, r in enumerate(rho):
        t = thetas[i % len(thetas)]
        out.append((-t if i % 2 and t else t, r))
    return out


def tau_inputs(stratum: Stratum) -> list[dict]:
    return [{"tau": [t, *r]} for t, r in _samples(stratum)]


def x_inputs(stratum: Stratum) -> list[dict]:
    return [{"z": so2.z_of(t), "t": r} for t, r in _samples(stratum)]
