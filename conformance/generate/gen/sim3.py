"""The Sim(3) ids (docs/decisions/0070, docs/PHASE5.md section 3, docs/maths/sim3.md), from their
definitions.

A tangent is [phi; rho; sigma], rotation first; an element X is (q, x, sigma), x its translation and
s = e^sigma its scale. The references:

    sim3_exp              (Exp phi as gen.so3's quaternion series, V rho, sigma), V = sum M^k/(k+1)!
                          with M = phi^ + sigma I, its scalars from gen.sim3coeff (the total-degree
                          series to convergence for |sigma| <= 3, the closed forms of SM.5(b) beyond)
    sim3_log              the geometric SO(3) logarithm of the normalized q (gen.so3.log, never
                          mp.logm, SM.8), and rho = V^-1 x by mp.lu_solve
    sim3_jr, sim3_jl      the dense series sum (-+ad)^k/(k+1)! of the 7 x 7 ad (SM.9) for |sigma| <= 3;
                          beyond, SM.10's blocks at |sigma| more digits, the series being e^|sigma| long
    sim3_jr_inv, _jl_inv  mp.inverse of those, the structural zeros of [A 0 0; B C f; 0 0 e] exact
    sim3_ad               the images of the basis under psi -> (X psi^ X^-1)^vee
    sim3_act              s R p + t

Strata (N_SIM3 records each, eight per ray; the budget's trade, docs/decisions/0071):

    sigma:<s>/theta:<t>   s in GRID_SIGMA, t in GRID_THETA: random axis, rho ~ N(0, I)
    ray:<r>               the joint limit, sigma = r cos(lambda), theta = r sin(lambda), 8 angles
    sigma:<s>/theta:<t>   large |sigma| (LARGE_SIGMA), theta in LARGE_THETA

each with its `@f32` twin, the same draws rounded once, where binary32 holds the range (LARGE_SIGMA
stops at 80 there). The X of `log`, `ad` and `act` is the reference Exp of the record's tangent,
rounded.
"""

from collections.abc import Callable

from mpmath import mp, mpf

from . import sim3coeff, so3
from .fmt import Mat
from .strata import Draw, deep, drawn_strata

N_SIM3 = deep(4)  # 0071: six would be 17 MB, past the budget; 0088: 64 deep
GRID_SIGMA = ("0", "1e-8", "-1e-8", "1e-4", "-1e-4", "1e-2", "-1e-2", "0.3", "-0.3", "1", "-1", "3", "-3")
GRID_THETA = ("1e-8", "1e-4", "1e-2", "0.3", "1", "3", "pi-1e-6")
RAYS = ("1e-12", "1e-8", "1e-4", "1e-2", "0.3", "0.8", "1.6")
N_RAY = 8
LARGE_SIGMA = {False: ("10", "-10", "30", "-30", "100", "-100", "700", "-700"),
               True: ("10", "-10", "30", "-30", "80", "-80")}
LARGE_THETA = ("1e-4", "1", "3")
SERIES_SIGMA = 3  # the dense series of ad, and gen.sim3coeff's series, up to |sigma| = 3
EXPM_SIGMA = 30  # the blocks are checked against mp.expm([ad I; 0 0]) up to |sigma| = 30


def value(s: str):
    return mp.pi - mpf(10) ** -6 if s == "pi-1e-6" else mpf(s)


# --- the algebra, rows of mpf ---------------------------------------------------------------------


def eye(n: int):
    return [[mpf(int(i == j)) for j in range(n)] for i in range(n)]


def mul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), mpf(0)) for j in range(len(B[0]))]
            for i in range(len(A))]


def cross(a, b) -> list:
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def hat7(tau):
    """The 4 x 4 hat of SM.1: [[phi^ + sigma I, rho], [0, 0]]."""
    W = so3.hat(tau[:3])
    H = [[W[i][j] + (tau[6] if i == j else 0) for j in range(3)] + [tau[3 + i]] for i in range(3)]
    return H + [[mpf(0)] * 4]


def vee7(H) -> list:
    """The tangent of a 4 x 4 algebra element: phi from the skew part, sigma from the trace."""
    sigma = (H[0][0] + H[1][1] + H[2][2]) / 3
    return [H[2][1], H[0][2], H[1][0], H[0][3], H[1][3], H[2][3], sigma]


def ad(tau):
    """ad_tau (SM.9): the images of the basis under psi -> [tau^, psi^]^vee, as columns."""
    T = hat7(tau)
    cols = []
    for j in range(7):
        S = hat7([mpf(int(i == j)) for i in range(7)])
        D = [[a - b for a, b in zip(ra, rb, strict=True)] for ra, rb in zip(mul(T, S), mul(S, T), strict=True)]
        cols.append(vee7(D))
    return [[cols[j][i] for j in range(7)] for i in range(7)]


def series(tau, sign: int):
    """sum (sign ad_tau)^k/(k+1)!, entrywise until two terms in a row are below so3's eps."""
    A = [[sign * a for a in row] for row in ad(tau)]
    eps = so3._eps()
    term = total = eye(7)
    quiet = 0
    for k in range(1, 100_000):
        term = [[x / (k + 1) for x in row] for row in mul(term, A)]
        total = [[s + t for s, t in zip(rs, rt, strict=True)] for rs, rt in zip(total, term, strict=True)]
        small = all(abs(t) <= eps * max(abs(s), eps) for rt, rs in zip(term, total) for t, s in zip(rt, rs))
        quiet = quiet + 1 if small else 0
        if quiet == 2:
            return total
    raise ArithmeticError("Sim(3) Jacobian series does not converge")


def coeffs(sigma, z) -> dict:
    """The eight scalars at (sigma, z): the series to convergence near the joint limit, the closed
    forms beyond |sigma| = SERIES_SIGMA (their cancellation is O(1) there)."""
    if abs(sigma) <= SERIES_SIGMA or z == 0:
        return sim3coeff.scalars(sigma, z)
    return sim3coeff.closed(sigma, z)


def V(phi, sigma, c=None):
    c = c or coeffs(sigma, sum(p * p for p in phi))
    W = so3.hat(phi)
    WW = mul(W, W)
    return [[(c["v0"] if i == j else 0) + c["v1"] * W[i][j] + c["v2"] * WW[i][j] for j in range(3)]
            for i in range(3)]


def blocks(tau, sign: int):
    """J_l(sign tau) by SM.10's blocks, dense: [[J_l(phi), 0, 0], [Q, V, -Gamma_2 rho], [0, 0, 1]]."""
    t = [sign * x for x in tau]
    phi, rho, sigma = t[:3], t[3:6], t[6]
    c = coeffs(sigma, sum(p * p for p in phi))
    W, X = so3.hat(phi), so3.hat(rho)
    WW = mul(W, W)
    Jl = so3.jl(phi)
    Vm = V(phi, sigma, c)
    mu = cross(phi, rho)
    nu = cross(phi, mu)
    y = [c["v0"] * rho[i] + c["v1"] * mu[i] + c["v2"] * nu[i] for i in range(3)]
    gam = [c["w0"] * rho[i] + c["w1"] * mu[i] + c["w2"] * nu[i] for i in range(3)]
    Y = so3.hat(y)
    YJ, XW, WX = mul(Y, Jl), mul(X, W), mul(W, X)
    dy = [2 * (c["dv1"] * mu[i] + c["dv2"] * nu[i]) for i in range(3)]
    Q = [[YJ[i][j] + dy[i] * phi[j] - c["v1"] * X[i][j] + c["v2"] * (XW[i][j] - 2 * WX[i][j])
          for j in range(3)] for i in range(3)]
    z = mpf(0)
    rows = [Jl[i] + [z, z, z, z] for i in range(3)]
    rows += [Q[i] + Vm[i] + [-gam[i]] for i in range(3)]
    rows += [[z] * 6 + [mpf(1)]]
    return rows


def structural(J):
    """J with the zeros of [A 0 0; B C f; 0 0 e] exact."""
    zero = lambda i, j: (i < 3 and j >= 3) or (i == 6 and j < 6)  # noqa: E731
    return [[mpf(0) if zero(i, j) else J[i][j] for j in range(7)] for i in range(7)]


def to_mat(A) -> Mat:
    return Mat((len(A), len(A[0])), [A[i][j] for j in range(len(A[0])) for i in range(len(A))])


# --- the ids ---------------------------------------------------------------------------------------


def _tau(inputs) -> list:
    return [mpf(c) for c in inputs["tau"]]


def _extra(sigma) -> int:
    return int(abs(sigma)) + 10


def exp_ref(tau) -> dict:
    phi, rho, sigma = tau[:3], tau[3:6], tau[6]
    q = so3.exp_series(phi)
    Vm = V(phi, sigma)
    x = [sum((Vm[i][j] * rho[j] for j in range(3)), mpf(0)) for i in range(3)]
    return {"q": q[:1] + so3.parallel(q[1:], phi), "x": x, "sigma": sigma}


def exp(inputs: dict) -> dict:
    tau = _tau(inputs)
    with mp.workdps(mp.dps + _extra(tau[6])):
        return exp_ref(tau)


def _x(inputs):
    q = so3.normalize([mpf(c) for c in inputs["q"]])
    return q, [mpf(c) for c in inputs["x"]], mpf(inputs["sigma"])


def log(inputs: dict) -> dict:
    q, x, sigma = _x(inputs)
    with mp.workdps(mp.dps + _extra(sigma)):
        phi = so3.log({"q": q})["phi"]
        rho = mp.lu_solve(mp.matrix(V(phi, sigma)), mp.matrix(x))
        return {"tau": list(phi) + [rho[i] for i in range(3)] + [sigma]}


def jacobian(name: str) -> Callable[[dict], dict]:
    sign = -1 if name.startswith("jr") else 1

    def evaluate(inputs: dict) -> dict:
        tau = _tau(inputs)
        with mp.workdps(mp.dps + _extra(tau[6])):
            J = series(tau, sign) if abs(tau[6]) <= SERIES_SIGMA else blocks(tau, sign)
            if name.endswith("inv"):
                B = mp.inverse(mp.matrix(J))
                J = [[B[i, j] for j in range(7)] for i in range(7)]
            return {"J": to_mat(structural(J))}

    return evaluate


def adjoint_of(q, x, sigma):
    """Ad_X as the images of the basis under psi -> (X psi^ X^-1)^vee."""
    R = so3.rot(q)
    s = mp.exp(sigma)
    X = [[s * R[i][j] for j in range(3)] + [x[i]] for i in range(3)] + [[mpf(0)] * 3 + [mpf(1)]]
    Rt = [[R[j][i] for j in range(3)] for i in range(3)]
    tinv = [-sum((Rt[i][k] * x[k] for k in range(3)), mpf(0)) / s for i in range(3)]
    Xi = [[Rt[i][j] / s for j in range(3)] + [tinv[i]] for i in range(3)] + [[mpf(0)] * 3 + [mpf(1)]]
    cols = [vee7(mul(mul(X, hat7([mpf(int(i == j)) for i in range(7)])), Xi)) for j in range(7)]
    return [[cols[j][i] for j in range(7)] for i in range(7)]


def ad_of(inputs: dict) -> dict:
    q, x, sigma = _x(inputs)
    with mp.workdps(mp.dps + _extra(sigma)):
        return {"Ad": to_mat(structural(adjoint_of(q, x, sigma)))}


def act(inputs: dict) -> dict:
    q, x, sigma = _x(inputs)
    with mp.workdps(mp.dps + _extra(sigma)):
        R, s = so3.rot(q), mp.exp(sigma)
        p = [mpf(c) for c in inputs["p"]]
        return {"Xp": [s * sum((R[i][j] * p[j] for j in range(3)), mpf(0)) + x[i] for i in range(3)]}


# --- inputs ------------------------------------------------------------------------------------------


def _axis(rng):
    z = 2 * mpf(rng.uniform()) - 1
    a = 2 * mp.pi * mpf(rng.uniform())
    r = mp.sqrt(1 - z * z)
    return [r * mp.cos(a), r * mp.sin(a), z]


def _normal(rng):
    """Box-Muller: two uniforms, one standard normal."""
    u1, u2 = mpf(rng.uniform()), mpf(rng.uniform())
    return mp.sqrt(-2 * mp.log(1 - u1)) * mp.cos(2 * mp.pi * u2)


def _taus(draw: Draw, count: int, cells: Callable[[int], tuple]) -> list[list[float]]:
    ra, rr = draw.stream("axis"), draw.stream("rho")
    out = []
    for i in range(count):
        sigma, theta = cells(i)
        axis = _axis(ra)
        phi = [draw.round(theta * a) for a in axis]
        rho = [draw.round(_normal(rr)) for _ in range(3)]
        out.append(phi + rho + [draw.round(sigma)])
    return out


def _cell(s: str, t: str):
    return lambda i: (value(s), value(t))


def _ray(r: str):
    return lambda i: (value(r) * mp.cos(mp.pi * (i + mpf(1) / 2) / N_RAY),
                      value(r) * mp.sin(mp.pi * (i + mpf(1) / 2) / N_RAY))


def _catalogue(f32: bool) -> list:
    out = []
    for s in GRID_SIGMA:
        for t in GRID_THETA:
            out.append((f"sigma:{s}/theta:{t}", _cell(s, t), N_SIM3, True))
    for r in RAYS:
        out.append((f"ray:{r}", _ray(r), N_RAY, True))
    for s in LARGE_SIGMA[f32]:
        for t in LARGE_THETA:
            out.append((f"sigma:{s}/theta:{t}", _cell(s, t), N_SIM3, s in LARGE_SIGMA[False] and s in LARGE_SIGMA[True]))
    return out


def _strata(make: Callable):
    return drawn_strata(
        [(name, make(cells), n) for name, cells, n, _ in _catalogue(False)],
        [(name, make(cells), n, shared) for name, cells, n, shared in _catalogue(True)],
    )


def _tau_make(cells):
    def make(draw: Draw, count: int) -> list[dict]:
        return [{"tau": t} for t in _taus(draw, count, cells)]

    return make


def _x_make(cells, with_point: bool):
    def make(draw: Draw, count: int) -> list[dict]:
        rp = draw.stream("point")
        out = []
        for t in _taus(draw, count, cells):
            tau = [mpf(c) for c in t]
            with mp.workdps(mp.dps + _extra(tau[6])):
                e = exp_ref(tau)
            rec = {
                "q": [draw.round(c) for c in e["q"]],
                "x": [draw.round(c) for c in e["x"]],
                "sigma": t[6],
            }
            if with_point:
                rec["p"] = [draw.round(c) for c in _axis(rp)]
            out.append(rec)
        return out

    return make


TAU_STRATA = _strata(_tau_make)
X_STRATA = _strata(lambda cells: _x_make(cells, False))
ACT_STRATA = _strata(lambda cells: _x_make(cells, True))
