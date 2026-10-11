"""The Sim(3) coefficient group's two data files (docs/decisions/0070, docs/maths/sim3.md SM.6, SM.11).

`sim3_series.jsonl` holds 1/n! for n = 0..N_FACT, correctly rounded at binary64 and binary32: the
constants of the total-degree series, which `cargo xtask thresholds` copies into
`coeffs/generated_sim3.rs` (0004: a series coefficient is generated, never typed).

`sim3_ref.jsonl` holds the eight scalars v0, v1, v2, w0, w1, w2, dz v1, dz v2 at points on rays
sigma = r cos(lambda), theta = r sin(lambda): RADII_PER_DECADE radii per decade of Delta = r^2 over
the switch grid's span [1e-16, 10], and ANGLES angles over [0, pi], so both signs of sigma and the
theta axis. A point's (sigma, z = theta^2) is rounded once to the precision, and its reference is
the scalars at those exact inputs: the sweep scores the shipped arms against it.

The reference is the total-degree series of SM.6(b) and SM.11(c), the recurrences of SM.3 summed
until a level is below 10^-(DPS + 10) of the first, at DPS digits plus what the series cancels
(e^(2|sigma|), SM.7's reference rule). The cross-check is independent: the closed forms of SM.5(b),
SM.11(b)-(c) at CHECK_DPS digits, where Delta > CHECK_DELTA, to 10^-90 of the scalar.
"""

from math import factorial

from mpmath import mp, mpf

from . import fmt, precision
from .precision import to_f32, to_f64

N_FACT = 32
RADII_PER_DECADE = 16
DECADES = (-16, 1)  # Delta from 1e-16 to 1e1, the switch grid's span
ANGLES = 13  # lambda = j pi / 12, j = 0..12
CHECK_DPS = 250
CHECK_DELTA = mpf(10) ** -20
NAMES = ("v0", "v1", "v2", "w0", "w1", "w2", "dv1", "dv2")


class Sim3RefError(Exception):
    pass


def build_series() -> tuple[bytes, int, int]:
    """1/n!, one record per n; `verified` counts the records whose two roundings were re-derived."""
    lines = []
    for n in range(N_FACT + 1):
        with mp.workdps(precision.DPS):
            x = mpf(1) / factorial(n)
            f64, f32 = to_f64(x), to_f32(x)
        if f64 != 1 / factorial(n) and n <= 22:
            # CPython's int / int is correctly rounded: an independent rounding of the same value.
            raise Sim3RefError(f"1/{n}! rounds to {f64!r}, CPython says {1 / factorial(n)!r}")
        lines.append(fmt.dumps({"f32": fmt.hex_float(f32), "f64": fmt.hex_float(f64), "n": n}))
    data = ("\n".join(lines) + "\n").encode()
    return data, len(lines), len(lines)


def scalars(sigma, z) -> dict:
    """The eight scalars at (sigma, z) by the total-degree recurrences, at the current mp.dps."""
    eps = mpf(10) ** -(mp.dps + 10)
    fact = [mpf(1)]
    out = dict.fromkeys(NAMES, mpf(0))
    xi, zeta, dxi, dzeta, p = mpf(0), mpf(0), mpf(0), mpf(0), mpf(1)  # xi_0 .. ; p = sigma^0
    k, quiet = 0, 0
    while quiet < 3:
        while len(fact) < k + 4:
            fact.append(fact[-1] * len(fact))
        terms = (
            p / fact[k + 1],
            xi / fact[k + 1],
            zeta / fact[k + 1],
            p / fact[k + 2],
            xi / fact[k + 2],
            zeta / fact[k + 2],
            dxi / fact[k + 1],
            dzeta / fact[k + 1],
        )
        for name, t in zip(NAMES, terms, strict=True):
            out[name] += t
        big = max(abs(t) for t in terms)
        quiet = quiet + 1 if k > 4 and big <= eps else 0
        # SM.3 and SM.11(c): xi_{k+1} = s xi + s^k - z zeta, zeta_{k+1} = s zeta + xi,
        # xi'_{k+1} = s xi' - zeta - z zeta', zeta'_{k+1} = s zeta' + xi'.
        xi, zeta, dxi, dzeta = (
            sigma * xi + p - z * zeta,
            sigma * zeta + xi,
            sigma * dxi - zeta - z * dzeta,
            sigma * dzeta + dxi,
        )
        p *= sigma
        k += 1
        if k > 4000:
            raise Sim3RefError(f"the series did not converge at ({sigma}, {z})")
    return out


def closed(sigma, z) -> dict:
    """SM.5(b), SM.11(b)-(c), unscaled, for Delta > 0 and sigma != 0, theta != 0."""
    th = mp.sqrt(z)
    e, m1 = mp.exp(sigma), mp.expm1(sigma)
    delta = sigma * sigma + z
    sinc = mp.sin(th) / th
    a = (1 - mp.cos(th)) / z
    b = (th - mp.sin(th)) / (z * th)
    d = (z - 4 * mp.sin(th / 2) ** 2) / (2 * z * z)
    v0 = m1 / sigma
    v2 = (v0 + sigma * e * a - e * sinc) / delta
    v1 = e * a - sigma * v2
    w0 = (m1 - sigma) / (sigma * sigma)
    w2 = (w0 + sigma * v2 - v1) / delta
    w1 = v2 - sigma * w2
    da, dsinc = -(b - 2 * d) / 2, -(a - b) / 2
    dv2 = (sigma * e * da - e * dsinc - v2) / delta
    dv1 = e * da - sigma * dv2
    return dict(zip(NAMES, (v0, v1, v2, w0, w1, w2, dv1, dv2), strict=True))


def points(f32: bool) -> list[tuple[float, float]]:
    """(sigma, z) on the rays, rounded once to the precision; duplicates (theta = 0) removed."""
    rnd = to_f32 if f32 else to_f64
    seen, out = set(), []
    with mp.workdps(precision.DPS):
        lo, hi = DECADES
        for i in range((hi - lo) * RADII_PER_DECADE + 1):
            r = mp.sqrt(mpf(10) ** (lo + mpf(i) / RADII_PER_DECADE))
            for j in range(ANGLES):
                lam = mp.pi * j / (ANGLES - 1)
                s, th = r * mp.cos(lam), r * mp.sin(lam)
                if j in (0, ANGLES - 1):
                    th = mpf(0)
                if j == (ANGLES - 1) // 2:
                    s = mpf(0)
                key = (rnd(s), rnd(th * th))
                if key not in seen:
                    seen.add(key)
                    out.append(key)
    return out


def _check(sigma, z, ref: dict) -> None:
    if z == 0 or sigma == 0 or sigma * sigma + z <= CHECK_DELTA:
        return
    with mp.workdps(CHECK_DPS):
        want = closed(mpf(sigma), mpf(z))
    for name in NAMES:
        if not abs(ref[name] - want[name]) <= mpf(10) ** -90 * max(abs(want[name]), mpf(10) ** -30):
            raise Sim3RefError(f"{name} at ({sigma!r}, {z!r}): {ref[name]} vs closed {want[name]}")


def build_ref() -> tuple[bytes, int, int]:
    lines, checked = [], 0
    for f32 in (False, True):
        for sigma, z in points(f32):
            with mp.workdps(precision.DPS + int(2 * abs(sigma)) + 20):
                ref = scalars(mpf(sigma), mpf(z))
                _check(sigma, z, ref)
                checked += z != 0 and sigma != 0 and sigma * sigma + z > CHECK_DELTA
            record = {"precision": "f32" if f32 else "f64"}
            record["sigma"], record["z"] = fmt.hex_float(sigma), fmt.hex_float(z)
            for name in NAMES:
                record[name] = fmt.dec30(ref[name])
            lines.append(fmt.dumps(record))
    data = ("\n".join(lines) + "\n").encode()
    return data, len(lines), checked
