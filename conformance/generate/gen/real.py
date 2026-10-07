"""The `real_*` ids (docs/decisions/0056): a `Dual<S, N>` method's value and its derivative, at
the working precision. The derivative is the calculus one (1/(2 sqrt x), 1/(3 cbrt(x)^2),
-1/sqrt(1 - x^2), (cos, -sin), (x, -y)/(x^2 + y^2), (1/d, -n/d^2)): textbook, not a closed form
of docs/NUMERICS.md. `check_real` holds it against `mp.diff`.
"""

from mpmath import mp, mpf

from .rng import SplitMix64
from .strata import Draw, drawn_strata, log_uniform_at

N_RANDOM = 64


def _sign(i: int) -> int:
    return -1 if i % 2 else 1


def _decades(lo: str, hi: str, alternate: bool):
    """Log-uniform in [lo, hi), the sign alternating along the records if `alternate`."""

    def make(draw: Draw, n: int) -> list[dict]:
        rng = draw.stream("x")
        a, b = mpf(lo), mpf(hi)
        return [
            {"x": (_sign(i) if alternate else 1) * log_uniform_at(rng, a, b, draw.f32)}
            for i in range(n)
        ]

    return make


# --- real_sqrt, real_cbrt --------------------------------------------------------------------

DECADES = {False: (-300, -100, -10, -1, 0, 1, 10, 100, 300), True: (-37, -10, -1, 0, 1, 10, 37)}
SUBNORMAL = {False: ("1e-310", "1e-309"), True: ("1e-40", "1e-39")}


def _root_strata(alternate: bool):
    def catalogue(f32: bool) -> list:
        return [
            *(
                (f"x:1e{e}", _decades(f"1e{e}", f"1e{e + 1}", alternate), e in DECADES[False])
                for e in DECADES[f32]
            ),
            ("x:subnormal", _decades(*SUBNORMAL[f32], alternate), False),
        ]

    return drawn_strata(
        [(name, make, N_RANDOM) for name, make, _ in catalogue(False)],
        [(name, make, N_RANDOM, shared) for name, make, shared in catalogue(True)],
    )


SQRT_STRATA = _root_strata(False)
CBRT_STRATA = _root_strata(True)


def sqrt(inputs: dict) -> dict:
    v = mp.sqrt(mpf(inputs["x"]))
    return {"value": v, "d": 1 / (2 * v)}


def cbrt(inputs: dict) -> dict:
    x = mpf(inputs["x"])
    v = mp.cbrt(abs(x)) * (-1 if x < 0 else 1)
    return {"value": v, "d": 1 / (3 * v * v)}


# --- real_sin_cos ----------------------------------------------------------------------------

SIN_COS = {
    "x:tiny": {False: ("1e-300", "1e-8"), True: ("1e-30", "1e-4")},
    "x:small": {False: ("1e-8", "1"), True: ("1e-4", "1")},
    "x:moderate": {False: ("1", "1e3"), True: ("1", "1e3")},
    "x:large": {False: ("1e3", "1e22"), True: ("1e3", "1e9")},
}
NEAR_K_MAX = 5  # k log-uniform in [1, 10^5]


def _near_k_pi_2(draw: Draw, n: int) -> list[dict]:
    """fl(k pi/2), k an integer log-uniform in [1, 10^5]; rounded once to the precision."""
    rng = draw.stream("k")
    return [
        {"x": draw.round(mp.floor(mpf(10) ** (NEAR_K_MAX * mpf(rng.uniform()))) * mp.pi / 2)}
        for _ in range(n)
    ]


def _sin_cos_strata():
    def catalogue(f32: bool) -> list:
        return [
            *(
                (name, _decades(*r[f32], alternate=True), r[True] == r[False])
                for name, r in SIN_COS.items()
            ),
            ("x:near-k-pi/2", _near_k_pi_2, True),
        ]

    return drawn_strata(
        [(name, make, N_RANDOM) for name, make, _ in catalogue(False)],
        [(name, make, N_RANDOM, shared) for name, make, shared in catalogue(True)],
    )


SIN_COS_STRATA = _sin_cos_strata()


def sin_cos(inputs: dict) -> dict:
    x = mpf(inputs["x"])
    s, c = mp.sin(x), mp.cos(x)
    return {"sin": s, "cos": c, "d_sin": c, "d_cos": -s}


# --- real_acos -------------------------------------------------------------------------------

NEAR_ONE_J = {False: (2, 52), True: (2, 23)}
ACOS_TINY = {False: ("1e-300", "1e-8"), True: ("1e-30", "1e-4")}


def _interior(draw: Draw, n: int) -> list[dict]:
    rng = draw.stream("x")
    return [{"x": draw.round(mpf("-0.99") + mpf("1.98") * mpf(rng.uniform()))} for _ in range(n)]


def _near_one(sign: int):
    """sign (1 - 2^-j), j uniform in the precision's range: exact at it."""

    def make(draw: Draw, n: int) -> list[dict]:
        rng = draw.stream("j")
        lo, hi = NEAR_ONE_J[draw.f32]
        return [
            {"x": sign * (1 - 2.0 ** -(lo + rng.next_u64() % (hi - lo + 1)))} for _ in range(n)
        ]

    return make


def _acos_strata():
    def catalogue(f32: bool) -> list:
        return [
            ("x:interior", _interior, True),
            ("x:near+1", _near_one(1), False),
            ("x:near-1", _near_one(-1), False),
            ("x:tiny", _decades(*ACOS_TINY[f32], alternate=True), False),
        ]

    return drawn_strata(
        [(name, make, N_RANDOM) for name, make, _ in catalogue(False)],
        [(name, make, N_RANDOM, shared) for name, make, shared in catalogue(True)],
    )


ACOS_STRATA = _acos_strata()


def acos(inputs: dict) -> dict:
    x = mpf(inputs["x"])
    return {"value": mp.acos(x), "d": -1 / mp.sqrt(1 - x * x)}


# --- real_atan2 ------------------------------------------------------------------------------

RATIO_K = {False: (8, 100, 300), True: (4, 15, 30)}
QUADRANTS = ((1, 1), (-1, 1), (-1, -1), (1, -1))  # (sign x, sign y), cycling along the records


def _radius(rng: SplitMix64):
    return mpf(10) ** (-3 + 6 * mpf(rng.uniform()))


def _generic_yx(draw: Draw, n: int) -> list[dict]:
    """The angle uniform in [-pi, pi), the radius log-uniform in [1e-3, 1e3]."""
    ra, rr = draw.stream("angle"), draw.stream("radius")
    out = []
    for _ in range(n):
        phi, r = mp.pi * (2 * mpf(ra.uniform()) - 1), _radius(rr)
        out.append({"y": draw.round(r * mp.sin(phi)), "x": draw.round(r * mp.cos(phi))})
    return out


def _ratio(exp: int):
    """|y/x| = 10^exp, the larger of the two log-uniform in [1e-3, 1e3], in each quadrant."""

    def make(draw: Draw, n: int) -> list[dict]:
        rr = draw.stream("radius")
        out = []
        for i in range(n):
            big = _radius(rr)
            small = big * mpf(10) ** -abs(exp)
            y, x = (small, big) if exp < 0 else (big, small)
            sx, sy = QUADRANTS[i % 4]
            out.append({"y": draw.round(sy * y), "x": draw.round(sx * x)})
        return out

    return make


def _atan2_strata():
    def catalogue(f32: bool) -> list:
        return [
            ("yx:generic", _generic_yx, True),
            *((f"yx:ratio-1e-{k}", _ratio(-k), False) for k in RATIO_K[f32]),
            *((f"yx:ratio-1e+{k}", _ratio(k), False) for k in RATIO_K[f32]),
        ]

    return drawn_strata(
        [(name, make, N_RANDOM) for name, make, _ in catalogue(False)],
        [(name, make, N_RANDOM, shared) for name, make, shared in catalogue(True)],
    )


ATAN2_STRATA = _atan2_strata()


def atan2(inputs: dict) -> dict:
    y, x = mpf(inputs["y"]), mpf(inputs["x"])
    r2 = x * x + y * y
    return {"value": mp.atan2(y, x), "d_y": x / r2, "d_x": -y / r2}


# --- real_div --------------------------------------------------------------------------------

DIV_WIDE = {False: ("1e-100", "1e100"), True: ("1e-12", "1e12")}


def _nd(lo: str, hi: str):
    """n, d log-uniform in [lo, hi); the signs of (n, d) cycle (+, +), (-, +), (+, -), (-, -)."""

    def make(draw: Draw, n: int) -> list[dict]:
        rn, rd = draw.stream("n"), draw.stream("d")
        a, b = mpf(lo), mpf(hi)
        return [
            {
                "n": _sign(i) * log_uniform_at(rn, a, b, draw.f32),
                "d": _sign(i // 2) * log_uniform_at(rd, a, b, draw.f32),
            }
            for i in range(n)
        ]

    return make


def _div_strata():
    def catalogue(f32: bool) -> list:
        return [("nd:generic", _nd("1e-3", "1e3"), True), ("nd:wide", _nd(*DIV_WIDE[f32]), False)]

    return drawn_strata(
        [(name, make, N_RANDOM) for name, make, _ in catalogue(False)],
        [(name, make, N_RANDOM, shared) for name, make, shared in catalogue(True)],
    )


DIV_STRATA = _div_strata()


def div(inputs: dict) -> dict:
    n, d = mpf(inputs["n"]), mpf(inputs["d"])
    return {"value": n / d, "d_n": 1 / d, "d_d": -n / (d * d)}
