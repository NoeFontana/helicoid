"""Generation-time cross-checks of the Sim(3) records (docs/decisions/0070), each by a route the
reference does not take:

- `exp`: `mp.expm` of the 4 x 4 hat for |sigma| <= EXPM_SIGMA (its blocks are s R and the
  translation); beyond, SM.5(a)'s three identities in V's coordinates, which tie v0, v1, v2 to
  e^sigma, sinc and a.
- `log`: Exp of the answer returns X (the rotation up to the quaternion's sign, the translation).
- `jr`, `jl`: the dense series against SM.10's blocks for |sigma| <= 3; the blocks against the
  top-right block of `mp.expm([[ad, I], [0, 0]])` (SM.7's reference) to EXPM_SIGMA; beyond,
  J_l = Ad_Exp(tau) J_r (SM.10(c)).
- the inverses: J J^-1 = I with J from the other route.
- `ad`: SM.9's closed form [R 0 0; [t]x R  sR  -t; 0 0 1].
- `act`: the quaternion sandwich q p q* (gen.so3.qmul), not the rotation matrix.
"""

from collections.abc import Callable

from mpmath import mp, mpf

from . import sim3, sim3coeff, so3
from .check import CrossCheckError

TOL = 90


def _close(what: str, got, want, scale=None) -> None:
    got, want = [mpf(g) for g in got], [mpf(w) for w in want]
    big = scale if scale is not None else max([mpf(1)] + [abs(w) for w in want])
    err = max(abs(g - w) for g, w in zip(got, want, strict=True))
    if not err <= mpf(10) ** -TOL * big:
        raise CrossCheckError(f"{what}: {mp.nstr(err, 5)} relative to {mp.nstr(big, 5)}")


def _flat(rows) -> list:
    return [x for r in rows for x in r]


def _rows(m) -> list:
    """A `fmt.Mat` (column-major) back to rows."""
    r, c = m.shape
    return [[m.data[j * r + i] for j in range(c)] for i in range(r)]


def exp(inputs: dict, out: dict) -> None:
    tau = [mpf(c) for c in inputs["tau"]]
    phi, sigma = tau[:3], tau[6]
    with mp.workdps(mp.dps + sim3._extra(sigma)):
        if abs(sigma) <= sim3.EXPM_SIGMA:
            E = mp.expm(mp.matrix(sim3.hat7(tau)))
            s = mp.exp(sigma)
            R = so3.rot(out["q"])
            _close("sim3_exp: s R", [E[i, j] for i in range(3) for j in range(3)],
                   [s * R[i][j] for i in range(3) for j in range(3)])
            _close("sim3_exp: translation", [E[i, 3] for i in range(3)], out["x"])
        else:
            z = sum(p * p for p in phi)
            th = mp.sqrt(z)
            c = sim3.coeffs(sigma, z)
            e, sinc, a = mp.exp(sigma), mp.sin(th) / th, (1 - mp.cos(th)) / z
            _close("sim3_exp: SM.5(a)",
                   [sigma * c["v0"], c["v0"] + sigma * c["v1"] - z * c["v2"], c["v1"] + sigma * c["v2"]],
                   [e - 1, e * sinc, e * a], abs(e))
    if out["sigma"] != tau[6]:
        raise CrossCheckError("sim3_exp: sigma is not the input's")


def log(inputs: dict, out: dict) -> None:
    tau = out["tau"]
    with mp.workdps(mp.dps + sim3._extra(tau[6])):
        e = sim3.exp_ref(tau)
        q = so3.normalize([mpf(c) for c in inputs["q"]])
        sign = 1 if sum(a * b for a, b in zip(q, e["q"], strict=True)) >= 0 else -1
        _close("sim3_log: rotation", e["q"], [sign * c for c in q])
        _close("sim3_log: translation", e["x"], [mpf(c) for c in inputs["x"]])


def jacobian(name: str) -> Callable[[dict, dict], None]:
    sign = -1 if name.startswith("jr") else 1

    def check(inputs: dict, out: dict) -> None:
        tau = [mpf(c) for c in inputs["tau"]]
        with mp.workdps(mp.dps + sim3._extra(tau[6])):
            got = _rows(out["J"])
            blocks = sim3.blocks(tau, sign)
            if abs(tau[6]) <= sim3.SERIES_SIGMA:
                other = sim3.series(tau, sign)
            elif abs(tau[6]) <= sim3.EXPM_SIGMA:
                A = sim3.ad([sign * t for t in tau])
                big = [[A[i][j] if j < 7 else mpf(int(j - 7 == i)) for j in range(14)] for i in range(7)]
                big += [[mpf(0)] * 14 for _ in range(7)]
                E = mp.expm(mp.matrix(big))
                other = [[E[i, 7 + j] for j in range(7)] for i in range(7)]
            else:
                # J_l(tau) = Ad_Exp(tau) J_r(tau) (SM.10(c)): the blocks at both signs and Ad.
                other = blocks
                e = sim3.exp_ref(tau)
                Ad = sim3.adjoint_of(so3.normalize(e["q"]), e["x"], tau[6])
                lhs = sim3.blocks(tau, 1)
                rhs = sim3.mul(Ad, sim3.blocks(tau, -1))
                _close(f"sim3_{name}: J_l = Ad J_r", _flat(lhs), _flat(rhs))
            _close(f"sim3_{name}: two routes", _flat(blocks), _flat(other))
            if name.endswith("inv"):
                eye = sim3.mul(got, blocks)
                _close(f"sim3_{name}: J J^-1", _flat(eye), _flat(sim3.eye(7)))
            else:
                _close(f"sim3_{name}: the record", _flat(got), _flat(other))

    return check


def ad(inputs: dict, out: dict) -> None:
    q = so3.normalize([mpf(c) for c in inputs["q"]])
    x, sigma = [mpf(c) for c in inputs["x"]], mpf(inputs["sigma"])
    with mp.workdps(mp.dps + sim3._extra(sigma)):
        R, s = so3.rot(q), mp.exp(sigma)
        T = so3.hat(x)
        TR = sim3.mul(T, R)
        z = mpf(0)
        want = [R[i] + [z] * 4 for i in range(3)]
        want += [TR[i] + [s * r for r in R[i]] + [-x[i]] for i in range(3)]
        want += [[z] * 6 + [mpf(1)]]
        _close("sim3_ad: SM.9", _flat(_rows(out["Ad"])), _flat(want))


def act(inputs: dict, out: dict) -> None:
    q = so3.normalize([mpf(c) for c in inputs["q"]])
    x, sigma = [mpf(c) for c in inputs["x"]], mpf(inputs["sigma"])
    with mp.workdps(mp.dps + sim3._extra(sigma)):
        p = [mpf(0)] + [mpf(c) for c in inputs["p"]]
        rp = so3.qmul(so3.qmul(q, p), so3.qconj(q))[1:]
        s = mp.exp(sigma)
        _close("sim3_act: quaternion sandwich", out["Xp"], [s * rp[i] + x[i] for i in range(3)])


_ = sim3coeff
