"""The polynomials and tables of the portable `atan2` candidate (docs/decisions/0105, draft).

Both kernels reduce to atan(v) = v + v s P(s), s = v^2, P the weighted minimax polynomial of
h(s) = (g(s) - 1)/s, g(s) = atan(sqrt s)/sqrt s, on [0, vmax^2], the weight w(s) = s/g(s) making
the error the one atan(v) carries: relative, v s (P - h)/atan(v). P's coefficients are rounded one
at a time, lowest first, each rounding followed by a refit of the rest: naive rounding of the
binary64 fit on |v| <= 1/2 loses five of its 62 bits.

- binary64: t = n/d in [0, 1] is reduced about the nearest centre c = j/C, j in 0..C, by
  atan(t) = atan(c) + atan(v), v = (n - c d)/(d + c n), so |v| <= 1/(2C). atan(c) is a pair hi + lo.
- binary32: evaluated in binary64, folded at t = 1/2 by (n - d)/(n + d), so |v| <= 1/2.

pi/4 is split as hi + lo with hi rounded to 51 bits, so k hi is exact for k in {0, ..., 4}.
No coefficient is typed anywhere; this module is their only source.
"""

import json
import struct

from mpmath import mp, mpf

from gen.precision import setup, signed_man_exp, to_f64

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


def remez(terms, top, fixed=(), iterations=40):
    """`terms` coefficients of P, lowest first, and the levelled weighted error, the first
    `len(fixed)` of them held at `fixed`."""
    xs = grid(top)
    m = len(fixed)
    free = terms - m
    hv = [h(x) - horner(list(fixed), x) for x in xs]
    wv = [x / g(x) for x in xs]
    if free == 0:
        return list(fixed), max(abs(w * hx) for w, hx in zip(wv, hv))
    ref = [top / 2 * (1 - mp.cos(mp.pi * (i + 1) / (free + 1))) for i in range(free + 1)]
    best = None
    for _ in range(iterations):
        a = mp.matrix(free + 1, free + 1)
        b = mp.matrix(free + 1, 1)
        for i, x in enumerate(ref):
            for j in range(free):
                a[i, j] = x ** (m + j)
            a[i, free] = (-1) ** i * g(x) / x
            b[i] = h(x) - horner(list(fixed), x)
        sol = mp.lu_solve(a, b)
        c = list(fixed) + [sol[j] for j in range(free)]
        level = abs(sol[free])
        e = [w * (horner(c[m:], x) * x**m - hx) for x, hx, w in zip(xs, hv, wv)]
        worst = max(abs(v) for v in e)
        if worst <= level * (1 + mpf(2) ** -20):
            return c, worst
        if best is None or worst < best[1]:
            best = (c, worst)
        try:
            ref = _alternating(xs, e, free + 1)
        except ArithmeticError:
            if not fixed:
                raise
            break
    # With rounded coefficients held, the rest may level only to the grid's resolution; the best
    # exchange still fixes the next coefficient to far below its rounding.
    if fixed and best is not None:
        return best
    raise ArithmeticError(f"no levelling after {iterations} exchanges")


def rounded(terms, top, rnd):
    """P's coefficients rounded by `rnd` one at a time, lowest first, each rounding followed by a
    refit of the rest; and the weighted error the rounded coefficients carry."""
    fixed = []
    for j in range(terms):
        c, _ = remez(terms, top, fixed)
        fixed.append(mpf(rnd(c[j])))
    return [float(x) for x in fixed], weighted_error(fixed, top)


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
    """The binary64 bits of `x`. Every value is a binary64, the binary32 kernel's coefficients too:
    that kernel evaluates in binary64."""
    return "0x" + struct.pack(">d", x).hex()


# The fold and the centres' comparisons round, so |v| may pass its bound by an ulp or so.
MARGIN = 1 + mpf(2) ** -40


def log2(v) -> float:
    return round(float(mp.log(v, 2)), 2)


def polynomial(terms: int, vmax) -> dict[str, object]:
    """P at binary64 on |v| <= vmax."""
    top = (vmax * MARGIN) ** 2
    _, level = remez(terms, top)
    coeffs, error = rounded(terms, top, to_f64)
    return {
        "terms": terms,
        "coeffs": [bits_hex(x) for x in coeffs],
        "minimax_log2": log2(level),
        "rounded_log2": log2(error),
    }


def binary64(centres: int, terms: int) -> dict[str, object]:
    """The binary64 kernel: `centres` a power of two, |v| <= 1/(2 centres)."""
    if centres < 1 or centres & (centres - 1):
        raise ValueError(f"{centres} centres: a power of two keeps n - c d exact")
    pi4 = mp.pi / 4
    hi = to_f64(round_bits(pi4, 51))
    lo = to_f64(pi4 - mpf(hi))
    atan_hi, atan_lo = [], []
    for j in range(centres + 1):
        a = mp.atan(mpf(j) / centres)
        a_hi = to_f64(a)
        atan_hi.append(bits_hex(a_hi))
        atan_lo.append(bits_hex(to_f64(a - mpf(a_hi))))
    return {
        "centres": centres,
        **polynomial(terms, mpf(1) / (2 * centres)),
        "atan_hi": atan_hi,
        "atan_lo": atan_lo,
        "pi4_hi": bits_hex(hi),
        "pi4_lo": bits_hex(lo),
    }


def binary32(terms: int) -> dict[str, object]:
    """The binary32 kernel, evaluated in binary64 and folded at 1/2."""
    return {"evaluated_at": "f64", "fold": "1/2", **polynomial(terms, mpf(1) / 2)}


def write(centres: int, terms64: int, terms32: int) -> str:
    """The committed file's text."""
    setup()
    out = {
        "source": "conformance/generate: python -m fit atan2 (docs/decisions/0105, draft)",
        "f64": binary64(centres, terms64),
        "f32": binary32(terms32),
    }
    return json.dumps(out, indent=2) + "\n"
