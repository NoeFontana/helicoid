"""The `real_*` ids (docs/decisions/0056): a `Dual<S, N>` method's value and its derivative, at
the working precision. The derivative is the calculus one (1/(2 sqrt x), 1/(3 cbrt(x)^2),
-1/sqrt(1 - x^2), (cos, -sin), (x, -y)/(x^2 + y^2), (1/d, -n/d^2), e^x for both `exp` and `expm1`):
textbook, not a closed form of docs/NUMERICS.md. `check_real` holds it against `mp.diff`.
`round_from_f64` is a constant (0087): its derivative is 0 and its value the precision's rounding.
"""

import math
import struct
from collections.abc import Callable

from mpmath import mp, mpf

from .precision import to_f32
from .rng import SplitMix64, log_uniform
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


# 0105 (draft): where a portable kernel's formulas meet, and the ends of the range. A record is
# (n, d), 0 <= n <= d, placed by `_place`: the quadrant cycles along the records, and in every
# other run of four |y| > |x| (the swap). `Dual::atan2`'s derivative is promised only where
# x^2 + y^2 is normal (docs/NUMERICS.md section 12), d in [2^-511, 2^511) ([2^-63, 2^63)), and
# every field of a record is scored: nothing beyond is drawn, so not the binary64 kernel's upper
# window edge 2^990, nor a d at 1e+-300, at FLT_MAX or subnormal.
#
# yx:centre-mid   binary64 only: n/d = (2j + 1)/16 + k ulp(n), j = 0..7 in turn every 8 records,
#                 k uniform in [-3, 3], d with its low 4 bits clear so that d (2j + 1)/16 is exact:
#                 the midpoints between the binary64 kernel's centres j/8, hit exactly at k = 0
# yx:fold         n = d/2 + k ulp(n), k uniform in [-4, 4]: the binary32 kernel's fold (one of the
#                 binary64 kernel's centres)
# yx:diagonal     n = d - k ulp(d), k uniform in [0, 4]: k pi/4, and the swap's tie at k = 0
# yx:axis         (y, x) = (+0, d), (+0, -d), (+-d, +0), (+-d, -0); never y = -0, which mpmath
#                 cannot tell from +0 (atan2(-0, -1) = -pi)
# yx:window-lo    binary64 only: d = 2^-460 + k ulp, k uniform in [-4, 4], n/d uniform in
#                 [1/16, 1]: the lower edge of the binary64 kernel's window
# yx:domain-low   d log-uniform in [2^-511, 2^-460) ([2^-63, 2^-50)), n/d uniform in [1/16, 1]:
#                 the binary64 kernel's rescale, as far down as the derivative is promised
# yx:domain-high  d log-uniform in [2^460, 2^511) ([2^50, 2^63)), n/d as above: the top of it
# yx:subnormal    n log-uniform in [1e-310, 1e-309) ([1e-40, 1e-39)), d in yx:domain-low's range
#
# d is `_radius`'s in the first four. Each stratum draws at its precision on its own stream, but
# yx:axis@f32, which rounds yx:axis's.

ATAN2_D = {False: (-511, -460, 460, 511), True: (-63, -50, 50, 63)}  # log2 of d's bounds
WINDOW_LO = -460


def _step(v: float, k: int, f32: bool) -> float:
    """`v > 0` moved `k` numbers of the precision, up for `k > 0`."""
    fmt, ints = ("<f", "<I") if f32 else ("<d", "<Q")
    bits = struct.unpack(ints, struct.pack(fmt, v))[0] + k
    return struct.unpack(fmt, struct.pack(ints, bits))[0]


def _place(i: int, n: float, d: float) -> dict:
    """Record `i` of `(n, d)`: in the quadrant `QUADRANTS[i % 4]`, with |y| > |x| for odd `i // 4`."""
    sx, sy = QUADRANTS[i % 4]
    y, x = (d, n) if (i // 4) % 2 else (n, d)
    return {"y": sy * y, "x": sx * x}


def _k(rng: SplitMix64, lo: int, hi: int) -> int:
    """Uniform in [lo, hi]."""
    return lo + rng.next_u64() % (hi - lo + 1)


def _near(draw: Draw, n: int, at: Callable[[int, float], float], k: tuple[int, int]) -> list[dict]:
    """d = `_radius`, and n `k` numbers of the precision from `at(i, d)`."""
    rr, rk = draw.stream("radius"), draw.stream("k")
    out = []
    for i in range(n):
        d = draw.round(_radius(rr))
        out.append(_place(i, _step(at(i, d), _k(rk, *k), draw.f32), d))
    return out


def _clear4(d: float) -> float:
    """`d` with the low 4 bits of its significand clear, so that `d (2j + 1)/16` is exact."""
    bits = struct.unpack("<Q", struct.pack("<d", d))[0]
    return struct.unpack("<d", struct.pack("<Q", bits >> 4 << 4))[0]


def _centre_mid(draw: Draw, n: int) -> list[dict]:
    rr, rk = draw.stream("radius"), draw.stream("k")
    out = []
    for i in range(n):
        d = _clear4(draw.round(_radius(rr)))
        mid = d * (2 * ((i // 8) % 8) + 1) / 16  # exact
        out.append(_place(i, _step(mid, _k(rk, -3, 3), False), d))
    return out


def _fold(draw: Draw, n: int) -> list[dict]:
    return _near(draw, n, lambda i, d: d / 2, (-4, 4))


def _diagonal(draw: Draw, n: int) -> list[dict]:
    return _near(draw, n, lambda i, d: d, (-4, 0))


def _axis(draw: Draw, n: int) -> list[dict]:
    rr = draw.stream("radius")
    out = []
    for i in range(n):
        d, s = draw.round(_radius(rr)), _sign(i // 4)
        y, x = ((0.0, d), (0.0, -d), (s * d, 0.0), (s * d, -0.0))[i % 4]
        out.append({"y": y, "x": x})
    return out


def _n_of(draw: Draw, rng: SplitMix64, d: float) -> float:
    """n = t d, t uniform in [1/16, 1], rounded once."""
    return draw.round(d * (1 + 15 * mpf(rng.uniform())) / 16)


def _window_lo(draw: Draw, n: int) -> list[dict]:
    rk, rt = draw.stream("k"), draw.stream("ratio")
    out = []
    for i in range(n):
        d = _step(2.0**WINDOW_LO, _k(rk, -4, 4), False)
        out.append(_place(i, _n_of(draw, rt, d), d))
    return out


def _d_at(draw: Draw, rng: SplitMix64, high: bool) -> float:
    lo, hi = ATAN2_D[draw.f32][2:] if high else ATAN2_D[draw.f32][:2]
    return log_uniform_at(rng, mpf(2) ** lo, mpf(2) ** hi, draw.f32)


def _domain(high: bool):
    def make(draw: Draw, n: int) -> list[dict]:
        rd, rt = draw.stream("d"), draw.stream("ratio")
        out = []
        for i in range(n):
            d = _d_at(draw, rd, high)
            out.append(_place(i, _n_of(draw, rt, d), d))
        return out

    return make


def _subnormal_n(draw: Draw, n: int) -> list[dict]:
    rn, rd = draw.stream("n"), draw.stream("d")
    lo, hi = (mpf(b) for b in SUBNORMAL[draw.f32])
    return [
        _place(i, log_uniform_at(rn, lo, hi, draw.f32), _d_at(draw, rd, False)) for i in range(n)
    ]


def _atan2_strata():
    def catalogue(f32: bool) -> list:
        return [
            ("yx:generic", _generic_yx, True),
            *((f"yx:ratio-1e-{k}", _ratio(-k), False) for k in RATIO_K[f32]),
            *((f"yx:ratio-1e+{k}", _ratio(k), False) for k in RATIO_K[f32]),
            # 0105's, after the earlier ones of their precision (0056 decision 2)
            *([] if f32 else [("yx:centre-mid", _centre_mid, False)]),
            ("yx:fold", _fold, False),
            ("yx:diagonal", _diagonal, False),
            ("yx:axis", _axis, True),
            *([] if f32 else [("yx:window-lo", _window_lo, False)]),
            ("yx:domain-low", _domain(False), False),
            ("yx:domain-high", _domain(True), False),
            ("yx:subnormal", _subnormal_n, False),
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


# --- real_exp, real_expm1 (0068) ---------------------------------------------------------------
#
# x:tiny      |x| log-uniform in [1e-300, 1e-8], both signs: expm1 must not cancel there
# x:unit      x uniform in [-1, 1]
# x:moderate  |x| uniform in [1, 30], both signs (binary32: [1, 15])
# x:large     x uniform in [30, ln MAX), the overflow onset (binary32: [15, ln MAX))
# x:negative  x uniform in (ln of the smallest normal, -30]: a normal e^x far below 1 (`exp` only)
# x:subnormal x uniform in (ln of half the smallest subnormal, ln of the smallest normal): e^x
#             subnormal, where `exp` is judged by its quantization (`exp` only)
# x:saturated x uniform in [-700, -40] (binary32: [-80, -20]): expm1 is -1 exactly and its
#             derivative e^x a normal number (`expm1` only)

EXP_RANGES = {
    "x:moderate": {False: ("1", "30"), True: ("1", "15")},
    "x:large": {False: ("30", "709.78"), True: ("15", "88.72")},
    "x:negative": {False: ("-708.39", "-30"), True: ("-87.33", "-15")},
    "x:subnormal": {False: ("-745.13", "-708.40"), True: ("-103.97", "-87.34")},
    "x:saturated": {False: ("-700", "-40"), True: ("-80", "-20")},
}


def _uniform_in(name: str, alternate: bool):
    def make(draw: Draw, n: int) -> list[dict]:
        rng = draw.stream("x")
        lo, hi = (mpf(b) for b in EXP_RANGES[name][draw.f32])
        out = []
        for i in range(n):
            x = draw.round(lo + (hi - lo) * mpf(rng.uniform()))
            out.append({"x": (_sign(i) if alternate else 1) * x})
        return out

    return make


def _unit(draw: Draw, n: int) -> list[dict]:
    rng = draw.stream("x")
    return [{"x": draw.round(2 * mpf(rng.uniform()) - 1)} for _ in range(n)]


def _exp_strata(names: tuple[str, ...]):
    def catalogue(f32: bool) -> list:
        out = [
            ("x:tiny", _decades(*ACOS_TINY[f32], alternate=True), False),
            ("x:unit", _unit, True),
        ]
        for name in names:
            out.append((name, _uniform_in(name, name == "x:moderate"), False))
        return out

    return drawn_strata(
        [(name, make, N_RANDOM) for name, make, _ in catalogue(False)],
        [(name, make, N_RANDOM, shared) for name, make, shared in catalogue(True)],
    )


EXP_STRATA = _exp_strata(("x:moderate", "x:large", "x:negative", "x:subnormal"))
EXPM1_STRATA = _exp_strata(("x:moderate", "x:large", "x:saturated"))


def exp(inputs: dict) -> dict:
    v = mp.exp(mpf(inputs["x"]))
    return {"value": v, "d": v}


def expm1(inputs: dict) -> dict:
    x = mpf(inputs["x"])
    return {"value": mp.expm1(x), "d": mp.exp(x)}


# --- real_round_from_f64 (0087) ----------------------------------------------------------------
#
# Inputs are binary64 at both precisions; an `@f32` stratum's are not binary32, by design. The
# reference is the input correctly rounded to the stratum's precision, ties to even: the input
# itself at binary64, `to_f32` of it at binary32.
#
# x:generic    |x| log-uniform in [1e-300, 1e300], both signs (binary32: [1e-37, 3e38])
# x:subnormal  |x| log-uniform in the binary64 subnormals [1e-320, 1e-308], both signs
# x:tie        binary32 only: the midpoint (2m + 1) 2^(e - 24) of two adjacent normal binary32,
#              m uniform in [2^23, 2^24), e uniform in [-126, 127], below the overflow edge
# x:near-tie   binary32 only: such a midpoint's binary64 neighbour, above and below alternately
# x:overflow   binary32 only: MAX + k 2^75 for k uniform in [-2^29, 2^28), 2^75 being the binary64
#              spacing there: from MAX's lower neighbour up to the edge (2 - 2^-24) 2^127, a tie
#              that rounds to inf and so is not scored; the first record is the binary64 below it
# x:underflow  binary32 only: |x| log-uniform in [2^-155, 2^-125] and, every other record, a tie
#              k 2^-150 for an odd k in [1, 2^24); the first two are 2^-150 (to 0) and 3 2^-150

ROUND_GENERIC = {False: ("1e-300", "1e300"), True: ("1e-37", "3e38")}
ROUND_SUBNORMAL = ("1e-320", "1e-308")
F32_MAX = (2**24 - 1) * 2**104
UNDERFLOW = (-155, -125)


def _signed(one: Callable[[SplitMix64, int], float]):
    """`one(rng, i)` for record `i`, the sign alternating along the records."""

    def make(draw: Draw, n: int) -> list[dict]:
        rng = draw.stream("x")
        return [{"x": _sign(i) * one(rng, i)} for i in range(n)]

    return make


def _log_uniform(lo, hi):
    a, b = mpf(lo), mpf(hi)
    return lambda rng, i: log_uniform(rng, a, b)


def _tie(rng: SplitMix64, i: int) -> float:
    """(2m + 1) 2^(e - 24): exact in binary64, which has 29 more bits than it needs."""
    e = -126 + rng.next_u64() % 254
    top = 2**24 - 1 if e == 127 else 2**24
    m = 2**23 + rng.next_u64() % (top - 2**23)
    return math.ldexp(2 * m + 1, e - 24)


def _near_tie(rng: SplitMix64, i: int) -> float:
    return math.nextafter(_tie(rng, i), math.inf if (i // 2) % 2 == 0 else 0.0)


def _overflow(rng: SplitMix64, i: int) -> float:
    k = 2**28 - 1 if i == 0 else -(2**29) + rng.next_u64() % (2**29 + 2**28)
    return float(F32_MAX + k * 2**75)


def _underflow(rng: SplitMix64, i: int) -> float:
    if i < 2:
        return math.ldexp(2 * i + 1, -150)
    if i % 2:
        return math.ldexp(2 * (rng.next_u64() % 2**23) + 1, -150)
    return log_uniform(rng, mpf(2) ** UNDERFLOW[0], mpf(2) ** UNDERFLOW[1])


ROUND_STRATA = drawn_strata(
    [
        ("x:generic", _signed(_log_uniform(*ROUND_GENERIC[False])), N_RANDOM),
        ("x:subnormal", _signed(_log_uniform(*ROUND_SUBNORMAL)), N_RANDOM),
    ],
    [
        ("x:generic", _signed(_log_uniform(*ROUND_GENERIC[True])), N_RANDOM, False),
        ("x:tie", _signed(_tie), N_RANDOM, False),
        ("x:near-tie", _signed(_near_tie), N_RANDOM, False),
        ("x:overflow", _signed(_overflow), N_RANDOM, False),
        ("x:underflow", _signed(_underflow), N_RANDOM, False),
    ],
)


def round_from_f64(inputs: dict, f32: bool) -> dict:
    x = mpf(inputs["x"])
    return {"value": mpf(to_f32(x)) if f32 else x, "d": mpf(0)}
