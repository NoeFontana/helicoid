"""The strata of docs/PHASE1.md section 4.4, each on its own splitmix64 stream.

A stratum draws at the working precision (`precision.DPS`) whatever `mp.dps` the caller has set, and
bounds are evaluated when it is drawn, never at import.
"""

from collections.abc import Callable
from dataclasses import dataclass

from mpmath import mp, mpf

from .precision import DPS, to_f64
from .rng import SEED, SplitMix64, log_uniform, stream, unit_quaternion_s3, unit_vector_s2

N_RANDOM = 64  # samples per random stratum
N_AXES = 64  # axes per axis-carrying stratum, for vector-valued ids
DENSE_PER_DECADE = 200
DENSE_DECADES = 4  # theta in [1e-4, 1]


@dataclass(frozen=True)
class Stratum:
    name: str
    count: int  # distinct theta values a scalar function sees
    draw: Callable[[SplitMix64], list[float]]  # binary64 thetas, in stream order
    carries_axes: bool = True  # PHASE1 section 4.4 lists an axis; `exact0` has none
    draw_quats: Callable[[SplitMix64], list[tuple]] | None = None  # a quaternion stratum's samples

    def thetas(self, seed: int = SEED) -> list[float]:
        with mp.workdps(DPS):
            return self.draw(stream(seed, self.name, "theta"))

    def axes(self, n: int = N_AXES, seed: int = SEED) -> list[tuple[float, float, float]]:
        if not self.carries_axes:
            return []
        rng = stream(seed, self.name, "axis")
        with mp.workdps(DPS):
            return [unit_vector_s2(rng) for _ in range(n)]

    def samples(self, seed: int = SEED) -> list[tuple[float, tuple[float, float, float]]]:
        """(theta, axis) of a theta stratum for a vector id: a random or dense theta has its own
        axis, a fixed one takes `N_AXES`; `exact0` has the axis (1, 0, 0), which nothing reads."""
        thetas = self.thetas(seed)
        if not self.carries_axes:
            return [(t, (1.0, 0.0, 0.0)) for t in thetas]
        n = max(self.count, N_AXES)
        return list(zip(thetas * (n // len(thetas)), self.axes(n, seed), strict=True))

    def quaternions(self, seed: int = SEED) -> list[tuple[float, float, float, float]]:
        """Binary64 `[w, x, y, z]` of a quaternion stratum (`q:w0`, `q:nonunit`)."""
        with mp.workdps(DPS):
            return self.draw_quats(stream(seed, self.name, "quat"))

    def points(self, n: int, seed: int = SEED) -> list[tuple[float, float, float]]:
        """`n` unit vectors for `so3_act`, uniform on S2."""
        rng = stream(seed, self.name, "point")
        with mp.workdps(DPS):
            return [unit_vector_s2(rng) for _ in range(n)]


def _log_uniform(name: str, bounds: Callable[[], tuple]) -> Stratum:
    def draw(rng: SplitMix64) -> list[float]:
        lo, hi = bounds()
        return [log_uniform(rng, lo, hi) for _ in range(N_RANDOM)]

    return Stratum(name, N_RANDOM, draw)


def _decade(exp: int) -> Callable[[], tuple]:
    return lambda: (mpf(10) ** exp, mpf(10) ** (exp + 1))


def _fixed(name: str, value: Callable[[], object], carries_axes: bool = True) -> Stratum:
    return Stratum(name, 1, lambda rng: [to_f64(value())], carries_axes)


def _dense(rng: SplitMix64) -> list[float]:
    n = DENSE_DECADES * DENSE_PER_DECADE
    return [to_f64(mpf(10) ** (mpf(j - n) / DENSE_PER_DECADE)) for j in range(n + 1)]


# Canonical order: it fixes record ids in every scalar-theta file.
SCALAR_THETA_STRATA = (
    *(_log_uniform(f"theta:1e{e}", _decade(e)) for e in range(-12, 0)),
    _log_uniform("theta:1e0", lambda: (mpf(1), mp.pi - mpf(1) / 10)),
    _fixed("theta:exact0", lambda: mpf(0), carries_axes=False),
    _log_uniform("theta:subnormal", _decade(-310)),
    *(_fixed(f"theta:pi-1e-{k}", lambda k=k: mp.pi - mpf(10) ** -k) for k in range(1, 13)),
    Stratum("theta:dense", DENSE_DECADES * DENSE_PER_DECADE + 1, _dense),
)

# `r` also sees a quaternion with w = +0 exactly (PHASE1 section 4.4): an angle of exactly pi, which
# no binary64 theta is, so `coeff.r_inputs` builds its records itself, one per norm n = |v|. At
# w = 0, r = pi / n and d_branch = -pi / (2 n^3): n != 1 exercises the 1 / n and n^-3 scalings.
Q_W0_NORMS = (1.0, 1e-3, 1e3)
Q_W0 = Stratum("q:w0", len(Q_W0_NORMS), lambda rng: [], carries_axes=False)
R_STRATA = (*SCALAR_THETA_STRATA, Q_W0)


# The quaternion ids (`so3_log`, `so3_act`, `so3_from_matrix`) see the theta strata as the unit
# quaternion of (theta, axis), then two strata of their own, in the order that keeps every id.
# `q:w0` is (+0, u) for a unit u, an angle of exactly pi. `q:nonunit` is a Haar-random unit
# quaternion scaled to |q|^2 - 1 = +-2^-45 (alternately), each component rounded to binary64.
NONUNIT_EXP = -45


def _q_w0(rng: SplitMix64) -> list[tuple]:
    return [(0.0, *unit_vector_s2(rng)) for _ in range(N_AXES)]


def _q_nonunit(rng: SplitMix64) -> list[tuple]:
    out = []
    for i in range(N_AXES):
        scale = mp.sqrt(1 + (-1) ** i * mpf(2) ** NONUNIT_EXP)
        out.append(tuple(to_f64(scale * c) for c in unit_quaternion_s3(rng)))
    return out


QUAT_STRATA = (
    *SCALAR_THETA_STRATA,
    Stratum("q:w0", N_AXES, lambda rng: [], carries_axes=False, draw_quats=_q_w0),
    Stratum("q:nonunit", N_AXES, lambda rng: [], carries_axes=False, draw_quats=_q_nonunit),
)
