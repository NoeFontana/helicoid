"""The strata of docs/PHASE1.md section 4.4, each on its own splitmix64 stream.

A stratum draws at the working precision (`precision.DPS`) whatever `mp.dps` the caller has set, and
bounds are evaluated when it is drawn, never at import.
"""

from collections.abc import Callable
from dataclasses import dataclass, replace

from mpmath import mp, mpf

from .precision import DPS, to_f32, to_f64
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
    rho_exp: int | None = None  # a `rho:*` SE_N(3) stratum's translation scale is 10^rho_exp
    f32: bool = False  # an `@f32` stratum: its binary64 draw rounded to nearest-even binary32
    source: "Stratum | None" = None  # the binary64 stratum an `@f32` one rounds; None: own draw

    def binary64(self) -> "Stratum":
        """The stratum whose binary64 draw an `@f32` stratum rounds (0016 item 1)."""
        return self.source or replace(self, f32=False)

    def thetas(self, seed: int = SEED) -> list[float]:
        """Binary64 thetas; for an `@f32` stratum also exactly binary32 ones."""
        if self.source is not None:
            thetas = self.source.thetas(seed)
        else:
            with mp.workdps(DPS):
                thetas = self.draw(stream(seed, self.name, "theta"))
        return [to_f32(t) for t in thetas] if self.f32 else thetas

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


def _log_uniform(name: str, bounds: Callable[[], tuple], count: int = N_RANDOM) -> Stratum:
    def draw(rng: SplitMix64) -> list[float]:
        lo, hi = bounds()
        return [log_uniform(rng, lo, hi) for _ in range(count)]

    return Stratum(name, count, draw)


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

# 0016 item 1: beside each stratum S of a scalar-theta id, `S@f32`, the same draw rounded to
# nearest-even binary32 (exactly representable in both precisions). `theta:subnormal` is the one
# whose rounding is not the reading: its decade [1e-310, 1e-309) is 0 in binary32, so its analogue
# is a binary32-subnormal decade, drawn from its own stream. In both, theta^2 underflows and 1 /
# theta overflows (binary32's smallest normal is 2^-126 ~ 1.2e-38, its largest 2^128 ~ 3.4e38).
F32_SUFFIX = "@f32"
F32_SUBNORMAL_DECADE = -40


def f32_twin(s: Stratum) -> Stratum:
    if s.name == "theta:subnormal":
        name = s.name + F32_SUFFIX
        return replace(_log_uniform(name, _decade(F32_SUBNORMAL_DECADE)), f32=True)
    return replace(s, name=s.name + F32_SUFFIX, f32=True, source=s)


SCALAR_F32_STRATA = tuple(f32_twin(s) for s in SCALAR_THETA_STRATA)
# The strata of `coeff_k`, `coeff_a`...`coeff_e` and `coeff_cos_half`: every `S`, then every
# `S@f32` in the same order, so no existing record moves.
COEFF_STRATA = (*SCALAR_THETA_STRATA, *SCALAR_F32_STRATA)

# `r` also sees a quaternion with w = +0 exactly (PHASE1 section 4.4): an angle of exactly pi, which
# no binary64 theta is, so `coeff.r_inputs` builds its records itself, one per norm n = |v|. At
# w = 0, r = pi / n and d_branch = -pi / (2 n^3): n != 1 exercises the 1 / n and n^-3 scalings.
Q_W0_NORMS = (1.0, 1e-3, 1e3)
Q_W0 = Stratum("q:w0", len(Q_W0_NORMS), lambda rng: [], carries_axes=False)
R_STRATA = (*SCALAR_THETA_STRATA, Q_W0)
COEFF_R_STRATA = (*R_STRATA, *SCALAR_F32_STRATA, f32_twin(Q_W0))


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


Q_STRATA = (
    Stratum("q:w0", N_AXES, lambda rng: [], carries_axes=False, draw_quats=_q_w0),
    Stratum("q:nonunit", N_AXES, lambda rng: [], carries_axes=False, draw_quats=_q_nonunit),
)
QUAT_STRATA = (*SCALAR_THETA_STRATA, *Q_STRATA)


# SE_N(3) ids (`sen3_*`) see the theta strata at unit translation scale, without `theta:dense` (a
# grid for the coefficients' switch points), then the `rho:*` strata: five translation scales
# crossed with five fixed thetas, each named `rho:<scale>/theta=<value>`. The two ids that take a
# quaternion, `sen3_log` and `sen3_ad`, also see `Q_STRATA`. A stratum is SEN3_SAMPLES records: 64
# would be 190 MB of dense matrices against a 50 MB corpus (docs/PHASE1.md section 4.4).
SEN3_SAMPLES = 6
RHO_EXPONENTS = (-6, -3, 0, 3, 4)
CELL_THETAS = (
    ("1e-8", lambda: mpf(10) ** -8),
    ("1e-4", lambda: mpf(10) ** -4),
    ("1e-1", lambda: mpf(10) ** -1),
    ("1", lambda: mpf(1)),
    ("pi-1e-6", lambda: mp.pi - mpf(10) ** -6),
)
SEN3_STRATA = (
    *(s for s in SCALAR_THETA_STRATA if s.name != "theta:dense"),
    *(
        Stratum(f"rho:1e{e}/theta={name}", 1, lambda rng, v=value: [to_f64(v())], rho_exp=e)
        for e in RHO_EXPONENTS
        for name, value in CELL_THETAS
    ),
)
SEN3_QUAT_STRATA = (*SEN3_STRATA, *Q_STRATA)


# SE(2) ids (`se2_*`) see SE_N(3)'s strata, `theta:dense` excluded for the same reason, with the
# angle's sign alternating along the samples: even records +theta, odd -theta. A stratum is
# SE2_SAMPLES records. The `so2_*` ids see `SCALAR_THETA_STRATA`, each theta in both signs.
SE2_SAMPLES = 8
SE2_STRATA = SEN3_STRATA


# The geodesic ids (`so3_geodesic`, `se3_geodesic`) see three strata of their own, PHASE4.md
# section 4's, adopted by 0045 item 3. A stratum's `theta` is the **relative** rotation theta(d),
# not a pose's; GEO_SAMPLES pose pairs, each crossed with every `t` of `gen.geodesic`. The pose
# translation scale cycles over `gen.geodesic.GEO_POSE_EXP` within every stratum, so each reaches
# the 1e4 that section 4 asks `geo:consecutive` for.
#
# `geo:consecutive` is tf_tree's kilohertz edge: the relative translation is at theta's own scale,
# so the whole of `d` is small (section 4 states the band on ||d||; ||d|| is then theta*sqrt(2),
# which is this band up to that factor, and the band on theta is what the draw holds exactly).
# `geo:generic` is the rotation range the other two leave, below pi - 0.1 as `theta:1e0` is.
# `geo:near-pi` stops strictly below pi and its margin is GEO_NEAR_PI_K: at pi the two preimages
# give geodesics O(1) apart, so the stratum measures conditioning and not agreement (0045 item 3).
GEO_SAMPLES = 6
GEO_NEAR_PI_K = (1, 2, 3, 6, 9, 12)
GEO_STRATA = (
    _log_uniform(
        "geo:consecutive", lambda: (mpf(10) ** -9, mpf(10) ** -3), count=GEO_SAMPLES
    ),
    _log_uniform(
        "geo:generic", lambda: (mpf(10) ** -3, mp.pi - mpf(1) / 10), count=GEO_SAMPLES
    ),
    Stratum(
        "geo:near-pi",
        len(GEO_NEAR_PI_K),
        lambda rng: [to_f64(mp.pi - mpf(10) ** -k) for k in GEO_NEAR_PI_K],
    ),
)


# The ids of docs/decisions/0056 (`solve_cubic`, `eig3`, `chol_*`, `quat_renormalize`, `real_*`)
# draw whole records, not thetas: a `Drawn` stratum's `make(draw, count)` builds its inputs from
# `draw.stream(purpose)` and rounds by `draw.round`: to binary64, or once to binary32 in an `@f32`
# stratum. An `@f32` stratum whose range is its binary64 namesake's reads that stratum's stream
# (`stream_of`), so its records are the same draws rounded once to binary32; one whose range 0056
# gives per precision, or that has no binary64 namesake, reads its own.
@dataclass(frozen=True)
class Draw:
    name: str  # the stratum whose streams are read
    f32: bool
    seed: int = SEED

    def stream(self, purpose: str, of: str | None = None) -> SplitMix64:
        """The stream of `purpose`; `of` names another stratum whose draws this one reuses."""
        return stream(self.seed, of or self.name, purpose)

    def round(self, x) -> float:
        return to_f32(x) if self.f32 else to_f64(x)


@dataclass(frozen=True)
class Drawn:
    name: str
    count: int
    make: Callable[[Draw, int], list[dict]]
    f32: bool = False
    stream_of: str | None = None  # a range shared with the binary64 namesake: its stream

    def records(self, seed: int = SEED) -> list[dict]:
        with mp.workdps(DPS):
            return self.make(Draw(self.stream_of or self.name, self.f32, seed), self.count)


def drawn_inputs(stratum: Drawn) -> list[dict]:
    return stratum.records()


def drawn_strata(
    binary64: list[tuple[str, Callable, int]], binary32: list[tuple[str, Callable, int, bool]]
) -> tuple[Drawn, ...]:
    """Every binary64 stratum `(name, make, count)`, then every binary32 one
    `(name, make, count, shared)`, named `name@f32` and reading `name`'s stream where `shared`."""
    return (
        *(Drawn(name, count, make) for name, make, count in binary64),
        *(
            Drawn(name + F32_SUFFIX, count, make, True, name if shared else None)
            for name, make, count, shared in binary32
        ),
    )


def log_uniform_at(rng: SplitMix64, lo, hi, f32: bool) -> float:
    """`log_uniform`'s draw, rounded once to binary32 when `f32`. A binary64 value must stay in
    [lo, hi); a binary32 one may lie half a unit outside, as in `theta:subnormal@f32`."""
    if not f32:
        return log_uniform(rng, lo, hi)
    return to_f32(lo * (hi / lo) ** rng.uniform())
