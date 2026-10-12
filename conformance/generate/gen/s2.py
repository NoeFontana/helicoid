"""`s2_retract` and `s2_local` (docs/PHASE5.md section 2, 0063), from their definitions.

The chart at a point n of S2 is frozen at n (NUMERICS.md section 8, docs/maths/charts.md CH.8):
sigma = -1 iff the stored n_z < 0, a comparison, so +0 and -0 both give +1 (0063 decision 1);
nu = n + sigma e_z, H = I - 2 nu nu^T / nu^T nu, b1 = H e_x, b2 = sigma H e_y, B = [b1 b2]: the
oriented frame, b1 x b2 = n on both hemispheres (0082).

    retract   m = Exp(B delta) n                      Rodrigues, CH.9(a)
    local     delta = B^T (alpha v / s)               v = n x m, s = |v|, alpha = atan2(s, n . m)

`local` is the geometric logarithm of the minimal rotation, never `mp.logm`, which is complex and
non-principal from alpha ~ 3.03 (CH.9, 0063 decision 4). A record's n and m are rounded unit
vectors and denote the points they are nearest to, so both are normalized first, as `gen.charts`
normalizes its quaternions; sigma is read from the stored n_z, which normalizing does not change.

Each stratum is N_S2 records of (n, delta), delta = theta (cos beta, sin beta) with beta uniform:

    s2:nz0            n = (cos lambda, sin lambda, z0), z0 = +0 on even records, -0 on odd ones,
                      the two readings of sgn(0) that 0063 decision 1 makes one; theta generic
    s2:near-antipode  n uniform, theta = pi - 10^-k, k cycling through NEAR_ANTIPODE_K: retract
                      lands near -n, and local's m is there, where |lj| = alpha / sin alpha (CH.10)
    s2:generic        n uniform, theta log-uniform in [1e-8, pi - 0.1]

Every component is drawn at the working precision and rounded once, to binary32 in an `@f32`
stratum (0067), which reads its binary64 namesake's streams. binary32 holds no pi - 1e-7 below pi,
so its antipode margins stop at 1e-6. `local`'s m is the reference retract of the same (n, delta),
rounded once, as `gen.charts.local_inputs` forms its Y.
"""

from collections.abc import Callable

from mpmath import mp, mpf

from .strata import Draw, drawn_strata

N_S2 = 64
NEAR_ANTIPODE_K = {False: (1, 2, 3, 6, 9, 12), True: (1, 2, 3, 4, 5, 6)}
GENERIC_THETA = (lambda: mpf(10) ** -8, lambda: mp.pi - mpf(1) / 10)


def norm(v):
    return mp.sqrt(sum(c * c for c in v))


def dot(a, b):
    return sum((x * y for x, y in zip(a, b, strict=True)), mpf(0))


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def unit(v) -> list:
    n = norm(v)
    return [mpf(c) / n for c in v]


def sigma(n_z: float) -> int:
    """0063 decision 1: -1 iff the stored n_z < 0; +0, -0 (and NaN) give +1."""
    return -1 if n_z < 0 else 1


def basis(n: list[float]) -> tuple[list, int, list, list]:
    """(n normalized, sigma, b1, b2): CH.8's Householder frame at the stored n, b2 times sigma so
    that b1 x b2 = n (0082)."""
    nh = unit([mpf(c) for c in n])
    s = sigma(n[2])
    nu = [nh[0], nh[1], nh[2] + s]
    k = 2 / dot(nu, nu)
    b1 = [(1 if i == 0 else 0) - k * nu[i] * nu[0] for i in range(3)]
    b2 = [s * ((1 if i == 1 else 0) - k * nu[i] * nu[1]) for i in range(3)]
    return nh, s, b1, b2


def rodrigues(omega, x) -> list:
    """Exp(omega) x."""
    theta = norm(omega)
    if theta == 0:
        return list(x)
    k = [c / theta for c in omega]
    c, s = mp.cos(theta), mp.sin(theta)
    kx, kdx = cross(k, x), dot(k, x)
    return [c * xi + s * a + (1 - c) * kdx * ki for xi, a, ki in zip(x, kx, k, strict=True)]


def _retract(n: list[float], delta: list[float]) -> list:
    nh, _, b1, b2 = basis(n)
    d1, d2 = (mpf(c) for c in delta)
    return rodrigues([d1 * a + d2 * b for a, b in zip(b1, b2, strict=True)], nh)


def retract(inputs: dict) -> dict:
    return {"m": _retract(inputs["n"], inputs["delta"])}


def rotation_vector(n: list[float], m: list[float]) -> tuple[list, list, list]:
    """(p, n normalized, m normalized): p = alpha v / s, the minimal rotation's Log; 0 at s = 0."""
    nh, mh = unit([mpf(c) for c in n]), unit([mpf(c) for c in m])
    v = cross(nh, mh)
    s = norm(v)
    if s == 0:
        return [mpf(0)] * 3, nh, mh
    alpha = mp.atan2(s, dot(nh, mh))
    return [alpha * c / s for c in v], nh, mh


def local(inputs: dict) -> dict:
    _, _, b1, b2 = basis(inputs["n"])
    p, _, _ = rotation_vector(inputs["n"], inputs["m"])
    return {"delta": [dot(b1, p), dot(b2, p)]}


# Inputs: exact binary64 (or binary32 in an `@f32` stratum), one dict per record.


def _uniform(rng, lo, hi):
    return lo + (hi - lo) * mpf(rng.uniform())


def _uniform_point(rng) -> list:
    z = 2 * mpf(rng.uniform()) - 1
    phi = 2 * mp.pi * mpf(rng.uniform())
    r = mp.sqrt(1 - z * z)
    return [r * mp.cos(phi), r * mp.sin(phi), z]


def _draws(name: str) -> Callable[[Draw, int], list[dict]]:
    def make(draw: Draw, count: int) -> list[dict]:
        rp, ra, rt = draw.stream("point"), draw.stream("angle"), draw.stream("theta")
        out = []
        for i in range(count):
            if name == "s2:nz0":
                lam = 2 * mp.pi * mpf(rp.uniform())
                n = [draw.round(mp.cos(lam)), draw.round(mp.sin(lam)), -0.0 if i % 2 else 0.0]
            else:
                n = [draw.round(c) for c in _uniform_point(rp)]
            beta = 2 * mp.pi * mpf(ra.uniform())
            if name == "s2:near-antipode":
                ks = NEAR_ANTIPODE_K[draw.f32]
                theta = mp.pi - mpf(10) ** -ks[i % len(ks)]
            else:
                lo, hi = (b() for b in GENERIC_THETA)
                theta = lo * (hi / lo) ** mpf(rt.uniform())
            delta = [draw.round(theta * mp.cos(beta)), draw.round(theta * mp.sin(beta))]
            if not norm([mpf(c) for c in delta]) < mp.pi:
                raise ArithmeticError(f"{name}: |delta| rounded to pi or past it: {delta}")
            out.append({"n": n, "delta": delta})
        return out

    return make


def _local_draws(name: str) -> Callable[[Draw, int], list[dict]]:
    draws = _draws(name)

    def make(draw: Draw, count: int) -> list[dict]:
        return [
            {"n": r["n"], "m": [draw.round(c) for c in _retract(r["n"], r["delta"])]}
            for r in draws(draw, count)
        ]

    return make


NAMES = ("s2:nz0", "s2:near-antipode", "s2:generic")


def _strata(make: Callable[[str], Callable]):
    return drawn_strata(
        [(name, make(name), N_S2) for name in NAMES],
        [(name, make(name), N_S2, True) for name in NAMES],
    )


RETRACT_STRATA = _strata(_draws)
LOCAL_STRATA = _strata(_local_draws)
