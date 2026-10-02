"""Taylor series of the coefficient catalogue as exact rationals (docs/NUMERICS.md section 4).

The committed series come from mpmath's Taylor expansion of the raw definitions
(docs/decisions/0004). Each coefficient is rationalized, reconstructed to 110 digits, and must
equal an independent exact derivation from the Maclaurin series of sin, cos and atan. The leading
four terms of the catalogue are asserted here, never used to produce data. `cos_half`, cos(theta/2),
is the eighth row: docs/NUMERICS.md section 3.1 asks for its series "generated alongside k", and
the catalogue lists none, so it is held to the exact derivation alone.
"""

from fractions import Fraction
from functools import cache
from math import factorial

from mpmath import mp, mpf

from . import coeff, fmt, precision

TERMS = 16
MAX_DEN = 10**50  # denominators reach ~2e43 at 16 terms; no wrong rational matches to 110 digits
RECONSTRUCT_DIGITS = 110
EXACT_TERMS = 232  # 1e-130 at x <= pi^2, even for `c` (radius of convergence 4 pi^2)

NAMES = ("k", "a", "b", "c", "d", "e", "r", "cos_half")
BRANCH = {"r": "n^2/w^2"}  # every other coefficient: theta^2
PREFACTOR = {"r": "2/w"}  # every other coefficient: 1

# docs/NUMERICS.md section 4, "Series (leading four)"; `cos_half` is not in that table.
DOC_LEADING = {
    "k": ("1/2", "-1/48", "1/3840", "-1/645120"),
    "a": ("1/2", "-1/24", "1/720", "-1/40320"),
    "b": ("1/6", "-1/120", "1/5040", "-1/362880"),
    "c": ("1/12", "1/720", "1/30240", "1/1209600"),
    "d": ("1/24", "-1/720", "1/40320", "-1/3628800"),
    "e": ("1/120", "-1/2520", "1/120960", "-1/9979200"),
    "r": ("1", "-1/3", "1/5", "-1/7"),
}


class SeriesError(Exception):
    pass


def to_fraction(x) -> Fraction:
    m, e = precision.signed_man_exp(x)
    return Fraction(m << e) if e >= 0 else Fraction(m, 1 << -e)


def rationalize(x) -> Fraction:
    value = to_fraction(x)
    q = value.limit_denominator(MAX_DEN)
    if abs(value - q) > abs(value) / 10**RECONSTRUCT_DIGITS:
        raise SeriesError(f"{mp.nstr(x, 30)} is not a rational with denominator <= 10^50")
    return q


def _quotient(p: list[Fraction], q: list[Fraction], n: int) -> list[Fraction]:
    out: list[Fraction] = []
    for j in range(n):
        out.append((p[j] - sum((out[i] * q[j - i] for i in range(j)), Fraction(0))) / q[0])
    return out


@cache
def exact() -> dict[str, tuple[Fraction, ...]]:
    """Every series by exact power-series arithmetic in x = theta^2 (s = n^2/w^2 for `r`).

    Each numerator is the definition rewritten in cos(t) and sinc(t) = sin(t) / t, whose Maclaurin
    series are the definitions of sin and cos, over a power of x; `r` is the series of atan(t) / t.
    """
    n = EXACT_TERMS
    m = n + 2
    sinc = [Fraction((-1) ** j, factorial(2 * j + 1)) for j in range(m)]
    cos = [Fraction((-1) ** j, factorial(2 * j)) for j in range(m)]
    one = [Fraction(int(j == 0)) for j in range(m)]
    x = [Fraction(int(j == 1)) for j in range(m)]
    # Numerators in x, each with its cancelling leading terms, which must be exactly zero.
    num = {
        "a": [one[j] - cos[j] for j in range(m)],  # (1 - cos) / x
        "b": [one[j] - sinc[j] for j in range(m)],  # (1 - sinc) / x
        "d": [x[j] + 2 * cos[j] - 2 * one[j] for j in range(m)],  # (x + 2 cos - 2) / 2 x^2
        "e": [2 * one[j] - 3 * sinc[j] + cos[j] for j in range(m)],  # (2 - 3 sinc + cos) / 2 x^2
        "c": [2 * sinc[j] - one[j] - cos[j] for j in range(m)],  # over x * 2 sinc
    }
    shift = {"a": 1, "b": 1, "d": 2, "e": 2, "c": 1}
    for name, series in num.items():
        if any(series[: shift[name]]):
            raise SeriesError(f"{name}: numerator does not vanish to order {shift[name]}")
    out = {
        "k": tuple(Fraction(1, 2) * sinc[j] / 4**j for j in range(n)),  # sin(t/2) / t
        "a": tuple(num["a"][1 : n + 1]),
        "b": tuple(num["b"][1 : n + 1]),
        "d": tuple(v / 2 for v in num["d"][2 : n + 2]),
        "e": tuple(v / 2 for v in num["e"][2 : n + 2]),
        "c": tuple(_quotient(num["c"][1:], [2 * s for s in sinc], n)),
        "r": tuple(Fraction((-1) ** j, 2 * j + 1) for j in range(n)),  # atan(t) / t
        "cos_half": tuple(cos[j] / 4**j for j in range(n)),  # cos(t / 2)
    }
    return out


@cache
def taylor(name: str) -> tuple[Fraction, ...]:
    """The first `TERMS` coefficients of `name` from mpmath, verified; the committed series."""
    precision.setup()
    f = coeff.BRANCH[name] if name != "r" else lambda y: coeff.r_branch(y, mpf(1)) / 2
    with mp.workdps(precision.DPS):
        raw = mp.taylor(f, 0, TERMS - 1)
        if any(abs(mp.im(c)) > mpf(10) ** (10 - mp.dps) for c in raw):
            raise SeriesError(f"{name}: Taylor coefficient with an imaginary part")
        series = tuple(rationalize(mp.re(c)) for c in raw)
    for j, (got, want) in enumerate(zip(series, exact()[name], strict=False)):
        if got != want:
            raise SeriesError(f"{name}: term {j} is {got} from mpmath, {want} by exact algebra")
    if name in DOC_LEADING and series[:4] != tuple(Fraction(t) for t in DOC_LEADING[name]):
        raise SeriesError(f"{name}: leading terms differ from docs/NUMERICS.md section 4")
    return series


def _str(q: Fraction) -> str:
    return f"{q.numerator}/{q.denominator}"


def build() -> tuple[bytes, int, int]:
    """(file bytes, record count, verified count): one record per coefficient, in `NAMES` order."""
    lines = [
        fmt.dumps(
            {
                "branch": BRANCH.get(name, "theta^2"),
                "coeff": name,
                "id": i,
                "prefactor": PREFACTOR.get(name, "1"),
                "series": [_str(q) for q in taylor(name)],
            }
        )
        + "\n"
        for i, name in enumerate(NAMES)
    ]
    return "".join(lines).encode(), len(lines), len(lines)
