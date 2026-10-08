"""`solve_cubic` (docs/decisions/0056): the three roots of a x^3 + b x^2 + c x + d at the stored
coefficients, ordered by (re, im), a real root's `im` exactly 0.

Which roots are real is decided exactly, by the sign of the discriminant in rationals. A real root
that is rational (every planted root) is found exactly: `mp.polyroots`' value is snapped to the
nearest fraction of denominator at most 2^64 and kept only if it is a root in rationals. A zero
discriminant has rational roots in closed form. Otherwise the roots are `mp.polyroots`' at the
working precision plus 300 bits, its error estimate held below the working precision.
"""

from fractions import Fraction

from mpmath import mp, mpf

from .precision import is_binary32, signed_man_exp
from .rng import SplitMix64
from .strata import Draw, drawn_strata, log_uniform_at

N_RANDOM = 64
# A planted root or leading coefficient is m 2^(e-6), 1 <= m <= 63, e in [-4, 4]: six significant
# bits, so d = -a r1 r2 r3 has at most 24 and every coefficient is exact at binary32. The roots of
# one cubic share e.
PLANTED_MAX = 63
SCALE_EXPS = range(-4, 5)
NEAR_DOUBLE_K = (2, 4, 6, 8)
P_SMALL = {False: ("1e-7", "1e-3"), True: ("1e-4", "1e-1")}
COEFF_SCALE = {False: 100, True: 20}
EXTRAPREC = 300  # bits
MAXSTEPS = 500


def _bits(rng: SplitMix64, n: int) -> int:
    return rng.next_u64() % n


def _planted(rng: SplitMix64, e: int) -> Fraction:
    m = 1 + _bits(rng, PLANTED_MAX)
    return Fraction(-m if _bits(rng, 2) else m, 2 ** (6 - e))


def _exponent(rng: SplitMix64) -> int:
    return SCALE_EXPS[_bits(rng, len(SCALE_EXPS))]


def expand(a, roots) -> list:
    """a (x - r1)(x - r2)(x - r3), descending."""
    r1, r2, r3 = roots
    return [a, -a * (r1 + r2 + r3), a * (r1 * r2 + r1 * r3 + r2 * r3), -a * r1 * r2 * r3]


def _exact_record(coeffs: list[Fraction]) -> dict:
    """Coefficients that are exact at both precisions, as the record's floats."""
    out = {}
    for key, c in zip("abcd", coeffs, strict=True):
        x = float(c)  # integer true division: correctly rounded
        if Fraction(x) != c or not is_binary32(x):
            raise ArithmeticError(f"planted coefficient {key} = {c} is not exact at binary32")
        out[key] = x
    return out


def planted_distinct(rng: SplitMix64) -> tuple[Fraction, list[Fraction]]:
    a, e = _planted(rng, _exponent(rng)), _exponent(rng)
    roots: list[Fraction] = []
    while len(roots) < 3:
        r = _planted(rng, e)
        if r not in roots:
            roots.append(r)
    return a, roots


def planted_double(rng: SplitMix64) -> tuple[Fraction, list[Fraction]]:
    a, e = _planted(rng, _exponent(rng)), _exponent(rng)
    r = _planted(rng, e)
    while (s := _planted(rng, e)) == r:
        pass
    return a, [r, r, s]


def planted_triple(rng: SplitMix64) -> tuple[Fraction, list[Fraction]]:
    a, e = _planted(rng, _exponent(rng)), _exponent(rng)
    r = _planted(rng, e)
    return a, [r, r, r]


def planted_one_real(rng: SplitMix64) -> tuple[Fraction, list[Fraction]]:
    """a (x - r)(x^2 + p x + q), p^2 < 4q: r, p as a root, q = m_q 2^(2e-8), so the complex pair
    has modulus below 2^(e-1); returns (a, [r, p, q])."""
    a, e = _planted(rng, _exponent(rng)), _exponent(rng)
    r = _planted(rng, e)
    m_q = 1 + _bits(rng, PLANTED_MAX)
    q = m_q * Fraction(2) ** (2 * e - 8)
    m_p = _bits(rng, 64)
    while m_p * m_p >= 64 * m_q:  # p^2 < 4q is m_p^2 < 64 m_q
        m_p = _bits(rng, 64)
    p = Fraction(-m_p if _bits(rng, 2) else m_p, 2 ** (6 - e))
    return a, [r, p, q]


def one_real_coefficients(a, r, p, q) -> list:
    return [a, a * (p - r), a * (q - r * p), -a * r * q]


def _planted_make(planted, one_real: bool = False):
    def make(draw: Draw, n: int) -> list[dict]:
        rng = draw.stream("cubic")
        out = []
        for _ in range(n):
            a, roots = planted(rng)
            coeffs = one_real_coefficients(a, *roots) if one_real else expand(a, roots)
            out.append(_exact_record(coeffs))
        return out

    return make


def _signed_log_uniform(rng: SplitMix64):
    """A root of scale 2^e, e in [-4, 4]: |r| log-uniform in [2^-4, 2^4], the sign a bit."""
    r = mpf(2) ** (8 * mpf(rng.uniform()) - 4)
    return -r if _bits(rng, 2) else r


def _near_double(k: int):
    """Roots r, r(1 + 10^-k), s at the working precision, |s - r| >= max(|r|, |s|) / 4, and a
    planted `a`; each coefficient rounded once to the stratum's precision."""

    def make(draw: Draw, n: int) -> list[dict]:
        rng = draw.stream("cubic")
        out = []
        for _ in range(n):
            a = mpf(_planted(rng, _exponent(rng)))
            r = _signed_log_uniform(rng)
            while abs((s := _signed_log_uniform(rng)) - r) < max(abs(r), abs(s)) / 4:
                pass
            coeffs = expand(a, (r, r * (1 + mpf(10) ** -k), s))
            out.append({key: draw.round(c) for key, c in zip("abcd", coeffs, strict=True)})
        return out

    return make


def _p_small(draw: Draw, n: int) -> list[dict]:
    """x^3 + p x - 1, p log-uniform in the precision's range."""
    rng = draw.stream("p")
    lo, hi = (mpf(b) for b in P_SMALL[draw.f32])
    return [
        {"a": 1.0, "b": 0.0, "c": log_uniform_at(rng, lo, hi, draw.f32), "d": -1.0}
        for _ in range(n)
    ]


def _coeff_scale(sign: int):
    """`cubic:distinct`'s cubics, from its stream, every coefficient times 10^(sign s), s the
    precision's, rounded once."""

    def make(draw: Draw, n: int) -> list[dict]:
        rng = draw.stream("cubic", of="cubic:distinct")
        factor = mpf(10) ** (sign * COEFF_SCALE[draw.f32])
        out = []
        for _ in range(n):
            a, roots = planted_distinct(rng)
            coeffs = expand(a, roots)
            out.append(
                {key: draw.round(mpf(c) * factor) for key, c in zip("abcd", coeffs, strict=True)}
            )
        return out

    return make


def _strata():
    planted = [
        ("cubic:distinct", _planted_make(planted_distinct)),
        ("cubic:double", _planted_make(planted_double)),
        ("cubic:triple", _planted_make(planted_triple)),
        ("cubic:one-real", _planted_make(planted_one_real, one_real=True)),
    ]
    near = [(f"cubic:near-double-1e-{k}", _near_double(k)) for k in NEAR_DOUBLE_K]
    rest = [
        ("cubic:one-real-p-small", _p_small, False),
        ("cubic:coeff-scale-up", _coeff_scale(1), False),
        ("cubic:coeff-scale-down", _coeff_scale(-1), False),
    ]
    shared = [(name, make, True) for name, make in planted + near]
    every = [*shared, *rest]
    return drawn_strata(
        [(name, make, N_RANDOM) for name, make, _ in every],
        [(name, make, N_RANDOM, same) for name, make, same in every],
    )


STRATA = _strata()


# --- the reference ---------------------------------------------------------------------------


def fraction(x) -> Fraction:
    """An mpf or a float, exactly."""
    m, e = signed_man_exp(mpf(x))
    return Fraction(m << e) if e >= 0 else Fraction(m, 1 << -e)


def horner(coeffs, z):
    value = 0 * z
    for c in coeffs:
        value = value * z + c
    return value


def discriminant(a, b, c, d):
    return 18 * a * b * c * d - 4 * b**3 * d + b * b * c * c - 4 * a * c**3 - 27 * a * a * d * d


def is_root(coeffs: list[Fraction], x: Fraction) -> bool:
    return horner(coeffs, x) == 0


def _snap(coeffs: list[Fraction], x) -> Fraction | None:
    q = fraction(x).limit_denominator(2**64)
    return q if is_root(coeffs, q) else None


def _multiple(coeffs: list[Fraction]) -> list[Fraction]:
    """The roots of a cubic with a zero discriminant, all rational."""
    a, b, c, d = coeffs
    s = b * b - 3 * a * c
    if s == 0:
        roots = [-b / (3 * a)] * 3
    else:
        double = (9 * a * d - b * c) / (2 * s)
        roots = [double, double, (4 * a * b * c - 9 * a * a * d - b**3) / (a * s)]
    if not all(is_root(coeffs, r) for r in roots):
        raise ArithmeticError(f"zero discriminant, roots {roots} are not roots of {coeffs}")
    return roots


def _polyroots(coeffs: list[Fraction]) -> list:
    approx, err = mp.polyroots(
        [mpf(c) for c in reversed(coeffs)],
        maxsteps=MAXSTEPS,
        extraprec=EXTRAPREC,
        error=True,
        asc=True,
    )
    scale = max(abs(z) for z in approx)
    if err > 16 * mp.eps * max(1, scale):  # converged: `err` is the working precision's eps
        raise ArithmeticError(f"polyroots error {mp.nstr(err, 5)} at {coeffs}")
    return list(approx)


def roots(coeffs: list[Fraction]) -> list[tuple]:
    """The three roots as (re, im), sorted."""
    a, b, c, d = coeffs
    disc = discriminant(a, b, c, d)
    if disc == 0:
        out = [(mpf(r), mpf(0)) for r in _multiple(coeffs)]
    elif disc > 0:
        out = []
        for z in _polyroots(coeffs):
            exact = _snap(coeffs, mp.re(z))
            out.append((mp.re(z) if exact is None else mpf(exact), mpf(0)))
    else:
        approx = _polyroots(coeffs)
        real = min(approx, key=lambda z: abs(mp.im(z)))
        pair = [z for z in approx if z is not real]
        r = _snap(coeffs, mp.re(real))
        if r is None:  # the pair made exactly conjugate
            re = (mp.re(pair[0]) + mp.re(pair[1])) / 2
            im = (abs(mp.im(pair[0])) + abs(mp.im(pair[1]))) / 2
            r_mp = mp.re(real)
        else:  # the pair from Vieta, exactly but for one square root
            re_q = (-b / a - r) / 2
            modulus2 = -d / (a * r) if r != 0 else c / a
            re, im, r_mp = mpf(re_q), mp.sqrt(mpf(modulus2 - re_q * re_q)), mpf(r)
        out = [(r_mp, mpf(0)), (re, -im), (re, im)]
    return sorted(out)


def solve_cubic(inputs: dict) -> dict:
    coeffs = [fraction(inputs[k]) for k in "abcd"]
    zs = roots(coeffs)
    return {"re": [z[0] for z in zs], "im": [z[1] for z in zs]}
