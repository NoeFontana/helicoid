"""The polynomial and the constants of the portable `atan2` candidate (docs/decisions/0105, draft).

The kernel folds its quotient to |u| <= umax and evaluates atan(u) = u + u s P(s), s = u^2. P is
the weighted minimax polynomial of h(s) = (g(s) - 1)/s, g(s) = atan(sqrt s)/sqrt s, on [0, umax^2],
the weight w(s) = s/g(s) making the error the one atan(u) carries: relative, u s (P - h)/atan(u).
No coefficient is typed anywhere; this module is their only source.

pi/4 is split as hi + lo with hi rounded to p - 2 bits, so k hi is exact for k in {0, ..., 4}.
"""

import json
import struct
from fractions import Fraction

from mpmath import mp, mpf

from gen.precision import setup, signed_man_exp, to_f32, to_f64

PRECISIONS = {"f64": (53, to_f64), "f32": (24, to_f32)}

# Points of the Chebyshev-distributed grid the exchange and the reported error read.
GRID = 2000


def g(s):
    """atan(sqrt s)/sqrt s, the series where the quotient cancels."""
    return mpf(1) + s * h(s)


def h(s):
    """(g(s) - 1)/s = sum_j (-1)^(j+1) s^j/(2j + 3), summed directly below 2^-20."""
    if s < mpf(2) ** -20:
        total, term, j = mpf(0), mpf(-1), 0
        while abs(term) > mpf(2) ** (-mp.prec - 8):
            total += term / (2 * j + 3)
            term *= -s
            j += 1
        return total
    r = mp.sqrt(s)
    return (mp.atan(r) / r - 1) / s


def horner(c, s):
    acc = mpf(0)
    for x in reversed(c):
        acc = acc * s + x
    return acc


def grid(top):
    """`GRID` Chebyshev-distributed points of (0, top], 0 excluded: the weight vanishes there."""
    return [top / 2 * (1 - mp.cos(mp.pi * i / (GRID - 1))) for i in range(1, GRID)]


def _alternating(xs, e, count):
    """`count` points of `xs` where `e` alternates in sign, each the largest |e| of its run."""
    peaks = []
    for i in range(len(xs)):
        left = e[i - 1] if i > 0 else None
        right = e[i + 1] if i + 1 < len(xs) else None
        if all(n is None or abs(e[i]) >= abs(n) for n in (left, right)):
            peaks.append(i)
    runs = []
    for i in peaks:
        if runs and (e[i] > 0) == (e[runs[-1]] > 0):
            if abs(e[i]) > abs(e[runs[-1]]):
                runs[-1] = i
        else:
            runs.append(i)
    while len(runs) > count:
        runs.pop(0 if abs(e[runs[0]]) < abs(e[runs[-1]]) else -1)
    if len(runs) < count:
        raise ArithmeticError(f"{len(runs)} alternation points, {count} needed")
    return [xs[i] for i in runs]


def remez(terms, top, iterations=40):
    """`terms` coefficients of P, lowest first, and the levelled weighted error."""
    xs = grid(top)
    hv = [h(x) for x in xs]
    wv = [x / g(x) for x in xs]
    ref = [top / 2 * (1 - mp.cos(mp.pi * (i + 1) / (terms + 1))) for i in range(terms + 1)]
    for _ in range(iterations):
        a = mp.matrix(terms + 1, terms + 1)
        b = mp.matrix(terms + 1, 1)
        for i, x in enumerate(ref):
            for j in range(terms):
                a[i, j] = x**j
            a[i, terms] = (-1) ** i * g(x) / x
            b[i] = h(x)
        sol = mp.lu_solve(a, b)
        c = [sol[j] for j in range(terms)]
        level = abs(sol[terms])
        e = [w * (horner(c, x) - hx) for x, hx, w in zip(xs, hv, wv)]
        worst = max(abs(v) for v in e)
        if worst <= level * (1 + mpf(2) ** -20):
            return c, worst
        ref = _alternating(xs, e, terms + 1)
    raise ArithmeticError(f"no levelling after {iterations} exchanges")


def weighted_error(c, top):
    """max over the grid of |u s (P(s) - h(s))/atan(u)|, P's coefficients exact as given."""
    return max(abs(x / g(x) * (horner(c, x) - h(x))) for x in grid(top))


def round_bits(x, bits):
    """`x` rounded to `bits` significant bits, ties to even."""
    m, e = signed_man_exp(mpf(x))
    n, sign = abs(m), -1 if m < 0 else 1
    drop = n.bit_length() - bits
    if drop <= 0:
        return mpf(x)
    q, rest = divmod(n, 1 << drop)
    half = 1 << (drop - 1)
    if rest > half or (rest == half and q & 1):
        q += 1
    return mpf(sign * q) * mpf(2) ** (e + drop)


def bits_hex(x: float) -> str:
    """The binary64 bits of `x`, which holds a binary32 exactly where the precision is `f32`."""
    return "0x" + struct.pack(">d", x).hex()


def fit(precision: str, terms: int, umax: Fraction) -> dict[str, object]:
    p, rnd = PRECISIONS[precision]
    top = mpf(umax.numerator) ** 2 / mpf(umax.denominator) ** 2
    exact, level = remez(terms, top)
    rounded = [rnd(x) for x in exact]
    pi4 = mp.pi / 4
    hi = rnd(round_bits(pi4, p - 2))
    lo = rnd(pi4 - mpf(hi))
    log2 = lambda v: round(float(mp.log(v, 2)), 2)
    return {
        "umax": f"{umax.numerator}/{umax.denominator}",
        "terms": terms,
        "coeffs": [bits_hex(x) for x in rounded],
        "minimax_log2": log2(level),
        "rounded_log2": log2(weighted_error([mpf(x) for x in rounded], top)),
        "pi4_hi": bits_hex(hi),
        "pi4_lo": bits_hex(lo),
    }


def write(specs: list[tuple[str, int, Fraction]]) -> str:
    """The committed file's text, one fit per precision."""
    setup()
    out: dict[str, object] = {
        "source": "conformance/generate: python -m fit atan2 (docs/decisions/0105, draft)",
    }
    for precision, terms, umax in specs:
        out[precision] = fit(precision, terms, umax)
    return json.dumps(out, indent=2) + "\n"
