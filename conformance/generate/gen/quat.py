"""`quat_renormalize` (docs/decisions/0056): q / |q| at the working precision, the exact
projection, not the Newton step it is scored against.

A Haar-random unit quaternion scaled to |q|^2 - 1 = +-2^-k (the sign alternating), rounded once to
the stratum's precision. Every stored q is inside docs/NUMERICS.md section 3.6's band,
| |q|^2 - 1 | <= 2^-26.29 (binary64) or 2^-11.79 (binary32), decided in rationals: a draw that
rounding takes outside it (only `renorm:eta-edge` can) is drawn again.
"""

from fractions import Fraction

from mpmath import mp, mpf

from .rng import unit_quaternion_s3
from .strata import Draw, drawn_strata

N_RANDOM = 64
ETA_K = {False: (27, 30, 40, 52), True: (12, 16, 20, 23)}
EDGE = {False: "26.29", True: "11.79"}  # the band is |eta| <= 2^-EDGE


def band(f32: bool):
    return mpf(2) ** -mpf(EDGE[f32])


def eta(q: list[float]) -> Fraction:
    """|q|^2 - 1 of the stored q, exactly."""
    return sum((Fraction(c) ** 2 for c in q), Fraction(0)) - 1


def _make(exponent=None):
    """`exponent` None is the band's edge."""

    def make(draw: Draw, n: int) -> list[dict]:
        rng, edge = draw.stream("quat"), band(draw.f32)
        size = edge if exponent is None else mpf(2) ** -exponent
        out = []
        for i in range(n):
            scale = mp.sqrt(1 + (-1) ** i * size)
            while True:
                q = [draw.round(scale * c) for c in unit_quaternion_s3(rng)]
                e = eta(q)
                if abs(mpf(e.numerator) / e.denominator) <= edge:
                    break
            out.append({"q": q})
        return out

    return make


def _strata(f32: bool) -> list:
    return [
        *((f"renorm:eta-2^-{k}", _make(k)) for k in ETA_K[f32]),
        ("renorm:eta-edge", _make()),
    ]


STRATA = drawn_strata(
    [(name, make, N_RANDOM) for name, make in _strata(False)],
    [(name, make, N_RANDOM, False) for name, make in _strata(True)],
)


def renormalize(inputs: dict) -> dict:
    q = [mpf(c) for c in inputs["q"]]
    r = mp.sqrt(sum(c * c for c in q))
    return {"q": [c / r for c in q]}
