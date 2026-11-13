"""The spline ids of docs/decisions/0116 (docs/PHASE1.md section 4.3), from their definition.

A record holds a window of `K` controls `X_0 .. X_{K-1}`, the segment parameter `u` and `inv_h`.
The value is `X_0 Exp(lambda_1 d_1) ... Exp(lambda_{K-1} d_{K-1})` (docs/NUMERICS.md section 21),
and nothing else of section 21 runs here: `gen/check_spline.py` evaluates its recursions.

- `lambda_j = sum_{s >= j} b_s`, `b_s` the B-spline basis by **Cox-de Boor's recursion** on the
  integer knots, as polynomials in `u` over exact rationals (`basis`), never section 21's integer
  table, which the check holds to this one.
- `d_j` is the **geometric** logarithm of `X_{j-1}^-1 X_j`: `gen.so3.log` (Newton on the
  quaternion series) and `gen.sen3.log`, as `gen.geodesic` takes it; `mp.logm` appears nowhere.
- `Exp` is mpmath's: the cosine and sine of the half angle for SO(3), `gen.sen3.exp` for SE(3).
- Every derivative output is the curve's derivative in `u`, which is `h` (`h^2`) times its time
  derivative: the ids' scaled outputs, so the reference does not read `inv_h`. `omega h` is
  `2 vec(q^* q')`, `alpha h^2` is `2 vec(q^* q'')` (`q^* q'` is pure, so `q'^* q'` is real), and
  SE(3)'s body twist is `[omega; R^T x']`, its derivative `[alpha; R^T x'' - omega x R^T x']`,
  with `q'`, `q''`, `x'`, `x''` by `mp.diff` at the working precision (it raises the precision
  itself).
- The SO(3) Jacobians are `mp.diff` along `X_m Exp(eps e_k)`, the value's right-trivialized,
  `2 vec(X^* X(eps))`, which is `Log(X^-1 X(eps))` to first order. Column `3m + k` of `jac`.
- R^n is exact: the window, `u` and the basis are rationals, so every output is a rational rounded
  once to the working precision, and the Jacobian `k` is `b_m` and its derivatives.
- `rn_spline_gram` is `h^{1-2r} int_0^1 b_s^(r) b_t^(r) du` in exact rationals.

A factor whose weight vanishes at `u` with every derivative up to the output's order is left out:
it contributes nothing up to that order (at `u = 0`, `lambda_{K-1} = u^{K-1}/(K-1)!`), and left in,
a control the output does not reach would read difference noise instead of an exact zero, which
the 150-digit recheck cannot compare.
"""

from fractions import Fraction
from functools import lru_cache

from mpmath import mp, mpf

from . import sen3, so3
from .fmt import Mat
from .rng import log_uniform, unit_quaternion_s3, unit_vector_s2
from .strata import BAND, F32_SUFFIX, SEN3_SAMPLES, Draw, Drawn

# The orders a record cycles through: 0113's, `ValidOrder`'s.
GROUP_ORDERS = (4, 5, 6)
RN_ORDERS = (2, 3, 4, 5, 6)


# The basis, exactly.


def _pmul(a: list, b: list) -> list:
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def _padd(a: list, b: list) -> list:
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def deriv(p: list) -> list:
    """`d/du` of a polynomial, coefficients lowest first."""
    return [k * c for k, c in enumerate(p)][1:] or [Fraction(0)]


def deriv_n(p: list, r: int) -> list:
    for _ in range(r):
        p = deriv(p)
    return p


@lru_cache(maxsize=None)
def basis(k: int) -> tuple:
    """`b_s(u)`, s = 0 .. K-1: the degree-(K-1) B-splines on the knots 0, 1, 2, ... active on
    [K-1, K), at `x = K - 1 + u`, by Cox-de Boor's recursion. Coefficients lowest first."""
    p = k - 1
    n = {i: [Fraction(int(i == p))] for i in range(2 * p + 1)}  # degree 0: 1 on [p, p + 1)
    for d in range(1, p + 1):
        n = {
            i: _padd(
                _pmul([Fraction(p - i, d), Fraction(1, d)], n[i]),  # (x - i) / d
                _pmul([Fraction(i + d + 1 - p, d), Fraction(-1, d)], n[i + 1]),  # (i+d+1-x)/d
            )
            for i in range(2 * p + 1 - d)
        }
    return tuple(tuple(n[s] + [Fraction(0)] * (k - len(n[s]))) for s in range(k))


@lru_cache(maxsize=None)
def cumulative(k: int) -> tuple:
    """`lambda_j(u) = sum_{s >= j} b_s(u)`, j = 0 .. K-1 (`lambda_0 = 1`)."""
    b = basis(k)
    return tuple(tuple(sum(b[s][c] for s in range(j, k)) for c in range(k)) for j in range(k))


def at(p, u):
    """A rational polynomial at `u`: exactly for a `Fraction`, at the precision for an `mpf`."""
    acc = Fraction(0) if isinstance(u, Fraction) else mpf(0)
    for c in reversed(p):
        acc = acc * u + (c if isinstance(u, Fraction) else mpf(c.numerator) / c.denominator)
    return acc


def _active(k: int, u: float, order: int) -> list[int]:
    """The increments j whose weight, or one of its derivatives up to `order`, is not 0 at `u`."""
    exact = Fraction(u)
    lam = cumulative(k)
    return [
        j
        for j in range(1, k)
        if any(at(deriv_n(list(lam[j]), r), exact) != 0 for r in range(order + 1))
    ]


# Derivatives of a vector curve in `u`.


def derivatives(curve, u, orders) -> dict:
    """`{r: [component d^r/du^r of curve at u]}` by `mp.diff` per component. The curve is
    evaluated once per point and precision, which `mp.diff` sets itself."""
    memo: dict = {}

    def point(t):
        key = (t, mp.prec)
        if key not in memo:
            memo[key] = curve(t)
        return memo[key]

    width = len(point(mpf(u)))
    out = {}
    for r in orders:
        out[r] = [
            point(mpf(u))[i] if r == 0 else mp.diff(lambda t, i=i: point(t)[i], mpf(u), r)
            for i in range(width)
        ]
    return out


# SO(3).


def qexp(phi) -> list:
    """mpmath's exponential of the pure quaternion `(0, phi/2)`: `(cos(theta/2), sin(theta/2) n)`,
    `n` the unit axis."""
    theta = so3.norm(phi)
    if theta == 0:
        return [mpf(1), mpf(0), mpf(0), mpf(0)]
    s = mp.sin(theta / 2) / theta
    return [mp.cos(theta / 2)] + [s * c for c in phi]


def _vec2(q) -> list:
    return [2 * c for c in q[1:]]


def quats(flat) -> list:
    """The unit quaternions a record's `q` denotes, one per control."""
    return [so3.normalize([mpf(c) for c in flat[4 * i : 4 * i + 4]]) for i in range(len(flat) // 4)]


def so3_increment(a, b) -> list:
    """`d = Log(a^-1 b)`, the geometric logarithm."""
    return list(so3.log({"q": so3.qmul(so3.qconj(a), b)})["phi"])


class So3Window:
    def __init__(self, controls: list, active: list[int]):
        self.n = controls
        self.active = active
        self.d = {j: so3_increment(controls[j - 1], controls[j]) for j in active}

    def perturbed(self, m: int, k: int, eps) -> "So3Window":
        """`X_m Exp(eps e_k)`: only `d_m` and `d_{m+1}` (and `X_0`) change."""
        e = [mpf(0)] * 3
        e[k] = eps
        out = So3Window.__new__(So3Window)
        out.n = list(self.n)
        out.n[m] = so3.qmul(self.n[m], qexp(e))
        out.active = self.active
        out.d = dict(self.d)
        for j in (m, m + 1):
            if j in self.active:
                out.d[j] = so3_increment(out.n[j - 1], out.n[j])
        return out

    def curve(self, k: int):
        lam = cumulative(k)

        def q(u):
            acc = self.n[0]
            for j in self.active:
                acc = so3.qmul(acc, qexp([at(lam[j], u) * c for c in self.d[j]]))
            return acc

        return q

    def kinematics(self, u, order: int) -> dict:
        """`{0: q, 1: omega h, 2: alpha h^2}` up to `order`."""
        ds = derivatives(self.curve(len(self.n)), u, range(order + 1))
        conj = so3.qconj(ds[0])
        out = {0: ds[0]}
        for r in range(1, order + 1):
            out[r] = _vec2(so3.qmul(conj, ds[r]))
        return out


def _so3_window(inputs: dict, order: int) -> So3Window:
    controls = quats(inputs["q"])
    return So3Window(controls, _active(len(controls), inputs["u"], order))


def so3_output(order: int):
    key = ("q", "omega", "alpha")[order]

    def evaluate(inputs: dict) -> dict:
        return {key: _so3_window(inputs, order).kinematics(inputs["u"], order)[order]}

    return evaluate


def so3_jacobian(order: int):
    """`[J_0 ... J_{K-1}]`, 3 x 3K: column `3m + k` is `d/deps` of the output at
    `X_m Exp(eps e_k)`, the value right-trivialized at the unperturbed `X`."""

    def evaluate(inputs: dict) -> dict:
        w = _so3_window(inputs, order)
        u, k = inputs["u"], len(w.n)
        base = w.kinematics(u, 0)[0] if order == 0 else None

        def output(window: So3Window) -> list:
            if order == 0:
                return _vec2(so3.qmul(so3.qconj(base), window.curve(k)(mpf(u))))
            return window.kinematics(u, order)[order]

        cols = []
        for m in range(k):
            for axis in range(3):
                memo: dict = {}

                def at_eps(eps, m=m, axis=axis, memo=memo):
                    key = (eps, mp.prec)
                    if key not in memo:
                        memo[key] = output(w.perturbed(m, axis, eps))
                    return memo[key]

                if m not in w.active and m + 1 not in w.active and m != 0:
                    cols.append([mpf(0)] * 3)  # no factor reads X_m: an exact zero (above)
                    continue
                cols.append([mp.diff(lambda e, i=i: at_eps(e)[i], 0) for i in range(3)])
        return {"jac": Mat((3, 3 * k), [c for col in cols for c in col])}

    return evaluate


# SE(3): a pose is (q, x), q unit.


def _compose(a, b):
    R = so3.rot(a[0])
    x = [a[1][r] + sum(R[r][c] * b[1][c] for c in range(3)) for r in range(3)]
    return (so3.qmul(a[0], b[0]), x)


def _rt(R, v) -> list:
    return [sum(R[c][r] * v[c] for c in range(3)) for r in range(3)]


def se3_increment(a, b) -> list:
    """`d = Log(a^-1 b)` through `gen.sen3.log`, the geometric logarithm and the J_l solve."""
    R = so3.rot(a[0])
    rel = {"q": so3.qmul(so3.qconj(a[0]), b[0]), "x": _rt(R, [y - x for x, y in zip(a[1], b[1])])}
    return list(sen3.log(1)(rel)["tau"])


def _sexp(tau):
    e = sen3.exp(1)({"tau": tau})
    return (e["q"], e["x"])


def poses(inputs: dict) -> list:
    qs = quats(inputs["q"])
    xs = [[mpf(c) for c in inputs["x"][3 * i : 3 * i + 3]] for i in range(len(qs))]
    return list(zip(qs, xs, strict=True))


def se3_output(order: int):
    def evaluate(inputs: dict) -> dict:
        ps = poses(inputs)
        k, u = len(ps), inputs["u"]
        active = _active(k, u, order)
        lam = cumulative(k)
        d = {j: se3_increment(ps[j - 1], ps[j]) for j in active}

        def curve(t):
            acc = ps[0]
            for j in active:
                acc = _compose(acc, _sexp([at(lam[j], t) * c for c in d[j]]))
            return acc[0] + acc[1]

        ds = derivatives(curve, u, range(order + 1))
        q, x = ds[0][:4], ds[0][4:]
        if order == 0:
            return {"q": q, "x": x}
        conj, R = so3.qconj(q), so3.rot(q)
        omega = _vec2(so3.qmul(conj, ds[1][:4]))
        v = _rt(R, ds[1][4:])
        if order == 1:
            return {"twist": omega + v}
        alpha = _vec2(so3.qmul(conj, ds[2][:4]))
        cross = [
            omega[(i + 1) % 3] * v[(i + 2) % 3] - omega[(i + 2) % 3] * v[(i + 1) % 3]
            for i in range(3)
        ]
        return {"twist": alpha + [a - c for a, c in zip(_rt(R, ds[2][4:]), cross)]}

    return evaluate


# R^n (n = 3), exactly.


def _points(inputs: dict) -> list:
    p = inputs["p"]
    return [[Fraction(c) for c in p[3 * i : 3 * i + 3]] for i in range(len(p) // 3)]


def _mp(x: Fraction):
    return mpf(x.numerator) / x.denominator


def rn_output(order: int):
    """`c_0 + sum_j lambda_j d_j` and its derivatives in `u`, the cumulative form in rationals."""
    key = ("p", "v", "a")[order]

    def evaluate(inputs: dict) -> dict:
        c, u = _points(inputs), Fraction(inputs["u"])
        lam = cumulative(len(c))
        acc = [x if order == 0 else Fraction(0) for x in c[0]]
        for j in range(1, len(c)):
            w = at(deriv_n(list(lam[j]), order), u)
            acc = [a + w * (y - x) for a, x, y in zip(acc, c[j - 1], c[j], strict=True)]
        return {key: [_mp(a) for a in acc]}

    return evaluate


def rn_jacobian(order: int):
    """The scalar `k` of each control's `k I`: `b_m^(order)(u)`, Cox-de Boor's basis."""

    def evaluate(inputs: dict) -> dict:
        k, u = len(inputs["p"]) // 3, Fraction(inputs["u"])
        return {"k": [_mp(at(deriv_n(list(b), order), u)) for b in basis(k)]}

    return evaluate


def gram_exact(k: int, r: int, inv_h: Fraction) -> list:
    """`G^(r)`, rows and columns the controls, exactly."""
    db = [deriv_n(list(b), r) for b in basis(k)]

    def integral(p):
        return sum((c / (i + 1) for i, c in enumerate(p)), Fraction(0))

    scale = inv_h ** (2 * r - 1)
    return [[integral(_pmul(db[s], db[t])) * scale for t in range(k)] for s in range(k)]


def gram(inputs: dict) -> dict:
    k, r = int(inputs["k"]), int(inputs["r"])
    g = gram_exact(k, r, Fraction(inputs["inv_h"]))
    return {"g": Mat((k, k), [_mp(g[s][t]) for t in range(k) for s in range(k)])}


# Strata (0116 decision 2): one axis at a time off the default cell `inc:generic`.


def _tiny(rng) -> mpf:
    return mpf(log_uniform(rng, mpf(10) ** -10, mpf(10) ** -8))


def _generic(rng) -> mpf:
    return mpf(log_uniform(rng, mpf(10) ** -3, mpf(1)))


def _large(rng) -> mpf:
    return 2.5 + (mp.pi - mpf(10) ** -3 - mpf(2.5)) * mpf(rng.uniform())


def _near_pi(rng) -> mpf:
    """Within 1e-6 of pi, never at it: `pi - 1e-6 (1 - 2^-10 v)`, the `pi-1e-6` cell's band."""
    return mp.pi - mpf(10) ** -6 * (1 - BAND * mpf(rng.uniform()))


INCREMENTS = {"generic": _generic, "tiny": _tiny, "large": _large, "near-pi": _near_pi}

# The u and h of a stratum that does not move them: the default cell's.
U_ONE_ULP = {False: 1 - mpf(2) ** -53, True: 1 - mpf(2) ** -24}


def _kinds(inc: str, i: int, k: int) -> list[str]:
    """The kind of each increment j = 1 .. K-1 of record `i`. `inc:mixed` has one tiny and one
    large increment, on slots that move with the record, the others generic; at K = 2 its one
    increment is tiny on even records and large on odd ones."""
    if inc != "mixed":
        return [inc] * (k - 1)
    if k == 2:
        return ["tiny" if i % 2 == 0 else "large"]
    kinds = ["generic"] * (k - 1)
    kinds[i % (k - 1)] = "tiny"
    kinds[(i + 1) % (k - 1)] = "large"
    return kinds


# The eight quaternions with one non-zero component, +-1: a binary32 anchor a tiny increment
# survives. A Haar-random binary32 control is O(1) in every component, and rounding its neighbour
# `X Exp(d)` to binary32 moves `d` by about u = 6e-8, more than the 1e-8 the stratum names; a
# product with one of these permutes and signs the components of `Exp(d)`, which round alone.
def _exact_unit(rng) -> list:
    i = int(8 * rng.uniform())
    q = [mpf(0)] * 4
    q[i % 4] = mpf(1 if i < 4 else -1)
    return q


def _window(group: str, inc: str, draw: Draw, streams: dict, i: int, k: int) -> dict:
    """The controls of one record, formed outward from an anchor and rounded at each step, so the
    increments the record holds are the stratum's (up to the roundings): `X_j = X_{j-1} Exp(d_j)`
    forward and `X_{j-1} = X_j Exp(-d_j)` backward. The anchor is the left end of the tiny
    increment where there is one, else `X_0`."""
    kinds = _kinds(inc, i, k)
    incs = []
    for kind in kinds:
        size = INCREMENTS[kind](streams["norm"])
        axis = [mpf(c) for c in unit_vector_s2(streams["axis"])]
        if group == "rn":
            incs.append([size * c for c in axis])
            continue
        phi = [size * c for c in axis]
        if group == "se3":
            scale = size if kind == "tiny" else mpf(1)
            phi += [scale * mpf(c) for c in unit_vector_s2(streams["rho"])]
        incs.append(phi)
    anchor = kinds.index("tiny") if "tiny" in kinds else 0
    # Every record draws its anchor whether it keeps it or not, so a twin's streams stay aligned.
    exact = draw.f32 and "tiny" in kinds
    if group == "rn":
        a = [mpf(c) for c in unit_vector_s2(streams["anchor"])]
        a = [mpf(0)] * 3 if exact else a
        ctl = {anchor: [draw.round(c) for c in a]}
        for j in range(anchor + 1, k):
            ctl[j] = [draw.round(mpf(x) + c) for x, c in zip(ctl[j - 1], incs[j - 1])]
        for j in range(anchor, 0, -1):
            ctl[j - 1] = [draw.round(mpf(x) - c) for x, c in zip(ctl[j], incs[j - 1])]
        return {"p": [c for j in range(k) for c in ctl[j]]}
    q = list(unit_quaternion_s3(streams["anchor"]))
    x = [mpf(c) for c in unit_vector_s2(streams["anchor-x"])]
    if exact:
        q, x = _exact_unit(streams["anchor-exact"]), [mpf(0)] * 3
    pose = {anchor: ([draw.round(c) for c in q], [draw.round(c) for c in x])}

    def step(stored, tau):
        n = so3.normalize([mpf(c) for c in stored[0]])
        if group == "so3":
            moved = (so3.qmul(n, so3.exp({"phi": tau})["q"]), [mpf(0)] * 3)
        else:
            moved = _compose((n, [mpf(c) for c in stored[1]]), _sexp(tau))
        return ([draw.round(c) for c in moved[0]], [draw.round(c) for c in moved[1]])

    for j in range(anchor + 1, k):
        pose[j] = step(pose[j - 1], incs[j - 1])
    for j in range(anchor, 0, -1):
        pose[j - 1] = step(pose[j], [-c for c in incs[j - 1]])
    out = {"q": [c for j in range(k) for c in pose[j][0]]}
    if group == "se3":
        out["x"] = [c for j in range(k) for c in pose[j][1]]
    return out


def _make(group: str, inc: str, u: str, h: str):
    orders = RN_ORDERS if group == "rn" else GROUP_ORDERS

    def make(draw: Draw, count: int) -> list[dict]:
        names = ("norm", "axis", "rho", "anchor", "anchor-x", "anchor-exact", "u")
        streams = {s: draw.stream(s) for s in names}
        out = []
        for i in range(count):
            k = orders[i % len(orders)]
            rec = _window(group, inc, draw, streams, i, k)
            uniform = mpf(streams["u"].uniform())
            us = {"random": uniform, "0": mpf(0), "half": mpf(1) / 2, "1-ulp": U_ONE_ULP[draw.f32]}
            rec["u"] = draw.round(us[u])
            rec["inv_h"] = draw.round({"1": mpf(1), "1e-3": mpf(1000), "1e3": mpf(1) / 1000}[h])
            out.append(rec)
        return out

    return make


def _catalogue(group: str) -> list:
    incs = ("generic", "large", "mixed") if group == "rn" else (
        "generic", "tiny", "large", "near-pi", "mixed"
    )
    cells = [(f"inc:{i}", i, "random", "1") for i in incs]
    cells += [(f"u:{u}", "generic", u, "1") for u in ("0", "half", "1-ulp")]
    cells += [(f"h:{h}", "generic", "random", h) for h in ("1e-3", "1e3")]
    return [(name, _make(group, i, u, h), SEN3_SAMPLES) for name, i, u, h in cells]


def strata(group: str) -> tuple[Drawn, ...]:
    """Every binary64 stratum, then its `@f32` twin on the same streams (0116 decision 2)."""
    cat = _catalogue(group)
    return (
        *(Drawn(name, count, make) for name, make, count in cat),
        *(Drawn(name + F32_SUFFIX, count, make, True, name) for name, make, count in cat),
    )


SO3_STRATA = strata("so3")
SE3_STRATA = strata("se3")
RN_STRATA = strata("rn")

# rn_spline_gram: one record per (K, r), r <= K - 1, per stratum `h:<h>`, h a power of two so the
# scaling is exact (0115).
GRAM_H = (("2^-10", 10), ("1", 0), ("2^10", -10))  # (name, log2 inv_h)
GRAM_CASES = tuple((k, r) for k in RN_ORDERS for r in (1, 2) if r <= k - 1)


def _gram_make(log2_inv_h: int):
    def make(draw: Draw, count: int) -> list[dict]:
        return [
            {"k": float(k), "r": float(r), "inv_h": draw.round(mpf(2) ** log2_inv_h)}
            for k, r in GRAM_CASES[:count]
        ]

    return make


GRAM_STRATA = (
    *(Drawn(f"h:{name}", len(GRAM_CASES), _gram_make(e)) for name, e in GRAM_H),
    *(Drawn(f"h:{name}{F32_SUFFIX}", len(GRAM_CASES), _gram_make(e), True) for name, e in GRAM_H),
)
