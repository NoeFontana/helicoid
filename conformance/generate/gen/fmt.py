"""Corpus text: hex-float inputs, 30-significant-digit decimal outputs, compact sorted JSON."""

import json

from .precision import signed_man_exp

SIG = 30


def hex_float(x: float) -> str:
    return float(x).hex()


def dec30(x) -> str:
    """An mpf as d.ddd...e<exp> with 30 significant digits, rounded half-to-even in integers."""
    m, e = signed_man_exp(x)
    if m == 0:
        return "0." + "0" * (SIG - 1) + "e0"
    num, den = (abs(m) << e, 1) if e >= 0 else (abs(m), 1 << -e)
    e10 = len(str(num)) - len(str(den))
    if num * 10 ** max(-e10, 0) < den * 10 ** max(e10, 0):
        e10 -= 1
    k = SIG - 1 - e10
    num, den = (num * 10**k, den) if k >= 0 else (num, den * 10**-k)
    q, r = divmod(num, den)
    if 2 * r > den or (2 * r == den and q & 1):
        q += 1
    if q == 10**SIG:
        q, e10 = q // 10, e10 + 1
    s = str(q)
    return f"{'-' if m < 0 else ''}{s[0]}.{s[1:]}e{e10}"


def dumps(record: dict) -> str:
    return json.dumps(record, separators=(",", ":"), sort_keys=True)
