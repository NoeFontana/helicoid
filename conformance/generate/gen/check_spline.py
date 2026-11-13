"""Generation-time cross-checks of the spline records (docs/decisions/0116 decision 3).

The reference (`gen.spline`) is a product of exponentials differentiated by `mp.diff`, on Cox-de
Boor's basis. Every record is held to what shares neither:

- **docs/NUMERICS.md section 21 itself**, evaluated here at the working precision: its integer
  basis `n_{s,k}`, the velocity and acceleration recursions, and for the SO(3) Jacobians the
  forward and backward passes with `J_r^-1` and `J_l^-1` as `mp.inverse` of the series. `Exp` is
  the quaternion series (`gen.so3.exp_series`, `gen.sen3.exp`) where the reference takes mpmath's
  sine and cosine, and `d_j` the closed `atan2` form where the reference takes Newton's method.
- **Continuity across a knot** (section 21's domains): the window shifted by one control, at
  `u = 0`, against the record's window at `u = 1`, to order `min(K - 2, 2)`, both by the reference.
- **R^n**: the non-cumulative sum `sum_s b_s c_s` on section 21's integers, exactly; each `k` as
  the exact difference of the value at `c_m = e_1` and `c_m = 0`, which is `d p / d c_m` since the
  value is linear in the controls.
- **Gram**: section 21's sum over `n_{s,p} n_{t,q}`, exactly; symmetric rows summing to zero; and
  the 30-digit reference rounding to the exact value's binary64 and binary32.
"""

from fractions import Fraction
from math import comb, factorial

from mpmath import mp, mpf

from . import sen3, so3, spline
from .check import DIGITS, CrossCheckError
from .fmt import dec30
from .precision import to_f32, to_f64


# Section 21's integers.


def numerators(k: int) -> list:
    """`n_{s,c} = C(K-1, c) sum_{l=s}^{K-1} (-1)^{l-s} C(K, l-s) (K-1-l)^{K-1-c}`."""
    return [
        [
            comb(k - 1, c)
            * sum(
                (-1) ** (l - s) * comb(k, l - s) * (k - 1 - l) ** (k - 1 - c) for l in range(s, k)
            )
            for c in range(k)
        ]
        for s in range(k)
    ]


def weights(k: int, u, order: int) -> list:
    """`lambda_j^(order)` in `u`, j = 0 .. K-1, by Horner on the cumulative numerators over
    `(K-1)!`; `u` a `Fraction` gives them exactly."""
    n = numerators(k)
    exact = isinstance(u, Fraction)
    out = []
    for j in range(k):
        cum = [sum(n[s][c] for s in range(j, k)) for c in range(k)]
        poly = [cum[c] * factorial(c) // factorial(c - order) for c in range(order, k)]
        acc = Fraction(0) if exact else mpf(0)
        for c in reversed(poly):
            acc = acc * u + c
        out.append(acc / factorial(k - 1))
    return out


def _agree(what: str, got, want, scale=None) -> None:
    """Every entry within 1e-DIGITS of the larger of 1 and the output's own size."""
    tol = mpf(10) ** -DIGITS * max(mpf(1), so3.maxabs(want), scale or 0)
    for g, w in zip(got, want, strict=True):
        if abs(mpf(g) - mpf(w)) > tol:
            raise CrossCheckError(f"{what}: {mp.nstr(mpf(g), 40)} vs {mp.nstr(mpf(w), 40)}")


# Matrices of mpf, rows.


def _mv(A, v) -> list:
    return [sum((A[i][c] * v[c] for c in range(len(v))), mpf(0)) for i in range(len(A))]


def _add(A, B):
    return [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def _scale(A, s):
    return [[s * a for a in row] for row in A]


def _zero(m: int):
    return [[mpf(0)] * m for _ in range(m)]


def log_atan2(q) -> list:
    """`2 atan2(|v|, w) v / |v|` of `q / |q|` after the flip: the closed form of section 3.2."""
    q = so3.canonical(q)
    n = so3.norm(q[1:])
    if n == 0:
        return [mpf(0)] * 3
    return [2 * mp.atan2(n, q[0]) / n * c for c in q[1:]]


class Recursion:
    """Section 21's forward pass on a group given by its `exp`, `compose`, `Ad`, `ad`."""

    def __init__(self, group: str, inputs: dict, order: int):
        self.group = group
        if group == "so3":
            self.x = spline.quats(inputs["q"])
            self.d = [None] + [
                log_atan2(so3.qmul(so3.qconj(a), b)) for a, b in zip(self.x, self.x[1:])
            ]
        else:
            self.x = spline.poses(inputs)
            self.d = [None] + [self._se3_log(a, b) for a, b in zip(self.x, self.x[1:])]
        k, u = len(self.x), mpf(inputs["u"])
        self.k = k
        self.lam = [weights(k, u, r) for r in range(3)]
        m = 3 if group == "so3" else 6
        self.m = m
        self.R, self.om, self.al = [None] * k, [[mpf(0)] * m], [[mpf(0)] * m]
        value = self.x[0]
        for j in range(1, k):
            a = self._exp([self.lam[0][j] * c for c in self.d[j]])
            value = self._compose(value, a)
            self.R[j] = self._adjoint(self._inverse(a))
            om = _add([_mv(self.R[j], self.om[-1])], [[self.lam[1][j] * c for c in self.d[j]]])[0]
            ad_d = _mv(self._ad(om), self.d[j])
            al = [
                x + self.lam[2][j] * dc + self.lam[1][j] * a_
                for x, dc, a_ in zip(_mv(self.R[j], self.al[-1]), self.d[j], ad_d)
            ]
            self.om.append(om)
            self.al.append(al)
        self.value = value

    @staticmethod
    def _se3_log(a, b) -> list:
        R = so3.rot(a[0])
        rel_q = so3.qmul(so3.qconj(a[0]), b[0])
        rel_x = [sum(R[c][r] * (b[1][c] - a[1][c]) for c in range(3)) for r in range(3)]
        phi = log_atan2(rel_q)
        rho = mp.lu_solve(mp.matrix(so3.jl(phi)), mp.matrix(rel_x))
        return phi + [rho[i] for i in range(3)]

    def _exp(self, tau):
        if self.group == "so3":
            return so3.exp_series(tau)
        e = sen3.exp(1)({"tau": tau})
        return (e["q"], e["x"])

    def _compose(self, a, b):
        return so3.qmul(a, b) if self.group == "so3" else spline._compose(a, b)

    def _inverse(self, a):
        if self.group == "so3":
            return so3.qconj(a)
        Rt = so3.mt(so3.rot(a[0]))
        return (so3.qconj(a[0]), [-c for c in _mv(Rt, a[1])])

    def _adjoint(self, a):
        if self.group == "so3":
            return so3.rot(a)
        return sen3.adjoint(so3.rot(a[0]), a[1], 1)

    def _ad(self, tau):
        return so3.hat(tau) if self.group == "so3" else sen3.ad(tau, 1)

    def jacobians(self, order: int) -> list:
        """`[J_0 ... J_{K-1}]` of the output of `order`, rows of a 3 x 3K matrix (SO(3))."""
        k, I = self.k, so3.identity()
        lam, dlam, ddlam = self.lam
        J = [None] + [so3.jr([lam[j] * c for c in self.d[j]]) for j in range(1, k)]
        jr_inv = [None] + [so3.inverse(so3.jr(self.d[j])) for j in range(1, k)]
        jl_inv = [None] + [so3.inverse(so3.jl(self.d[j])) for j in range(1, k)]
        p_om, p_al = [None] * k, [None] * k
        for j in range(1, k):
            hd = so3.hat(self.d[j])
            p_om[j] = _add(
                _scale(so3.mm(so3.hat(_mv(self.R[j], self.om[j - 1])), J[j]), lam[j]),
                _scale(I, dlam[j]),
            )
            p_al[j] = _add(
                _add(
                    _scale(so3.mm(so3.hat(_mv(self.R[j], self.al[j - 1])), J[j]), lam[j]),
                    _scale(I, ddlam[j]),
                ),
                _add(
                    _scale(so3.hat(self.om[j]), dlam[j]),
                    _scale(so3.mm(hd, p_om[j]), -dlam[j]),
                ),
            )
        T, C = I, _zero(3)
        by_d = [None] * k
        for j in range(k - 1, 0, -1):
            by_d[j] = (
                _scale(so3.mm(T, J[j]), lam[j]),
                so3.mm(T, p_om[j]),
                _add(so3.mm(C, p_om[j]), so3.mm(T, p_al[j])),
            )[order]
            C = _add(
                so3.mm(C, self.R[j]),
                _scale(so3.mm(so3.mm(T, so3.hat(self.d[j])), self.R[j]), -dlam[j]),
            )
            T = so3.mm(T, self.R[j])
        blocks = []
        for m in range(k):
            acc = T if (m == 0 and order == 0) else _zero(3)
            if m >= 1:
                acc = _add(acc, so3.mm(by_d[m], jr_inv[m]))
            if m + 1 <= k - 1:
                acc = _add(acc, _scale(so3.mm(by_d[m + 1], jl_inv[m + 1]), -1))
            blocks.append(acc)
        return blocks


def _knot(evaluate, inputs: dict, order: int, k: int, fields: tuple) -> None:
    """The window shifted by one control (the last repeated: at `u = 0` its weight vanishes to
    order K - 2) at `u = 0` against the record's at `u = 1`, where order <= K - 2."""
    if order > k - 2:
        return
    shifted = {}
    for key, width in fields:
        v = list(inputs[key])
        shifted[key] = v[width:] + v[-width:]
    one = evaluate({**inputs, "u": 1.0})
    zero = evaluate({**inputs, **shifted, "u": 0.0})
    for key in one:
        a, b = one[key], zero[key]
        if key == "q" and sum(x * y for x, y in zip(a, b)) < 0:
            b = [-c for c in b]
        _agree(f"spline: continuity across a knot, {key}", a, b)


def group_output(group: str, order: int):
    evaluate = (spline.so3_output if group == "so3" else spline.se3_output)(order)
    fields = (("q", 4),) if group == "so3" else (("q", 4), ("x", 3))

    def check(inputs: dict, out: dict) -> None:
        rec = Recursion(group, inputs, order)
        if order == 0:
            if group == "so3":
                _agree("spline: value against section 21", out["q"], rec.value)
            else:
                _agree("spline: value against section 21", out["q"], rec.value[0])
                _agree("spline: translation against section 21", out["x"], rec.value[1])
        else:
            want = (rec.om if order == 1 else rec.al)[-1]
            key = {"so3": ("omega", "alpha")[order - 1], "se3": "twist"}[group]
            _agree(f"spline: order {order} against section 21", out[key], want)
        _knot(evaluate, inputs, order, len(inputs["q"]) // 4, fields)

    return check


def so3_jacobian(order: int):
    def check(inputs: dict, out: dict) -> None:
        blocks = Recursion("so3", inputs, order).jacobians(order)
        want = [blocks[c // 3][r][c % 3] for c in range(3 * len(blocks)) for r in range(3)]
        _agree(f"spline: Jacobian of order {order} against section 21", out["jac"].data, want)

    return check


def rn_output(order: int):
    key = ("p", "v", "a")[order]
    evaluate = spline.rn_output(order)

    def check(inputs: dict, out: dict) -> None:
        p = inputs["p"]
        c = [[Fraction(x) for x in p[3 * i : 3 * i + 3]] for i in range(len(p) // 3)]
        k, u = len(c), Fraction(inputs["u"])
        b = [w - v for w, v in zip(weights(k, u, order), weights(k, u, order)[1:] + [0])]
        want = [sum((b[s] * c[s][i] for s in range(k)), Fraction(0)) for i in range(3)]
        _agree("spline: R^n against sum b_s c_s", out[key], [spline._mp(w) for w in want])
        _knot(evaluate, inputs, order, k, (("p", 3),))

    return check


def rn_jacobian(order: int):
    value = spline.rn_output(order)
    key = ("p", "v", "a")[order]

    def check(inputs: dict, out: dict) -> None:
        k, u = len(inputs["p"]) // 3, Fraction(inputs["u"])
        lam = weights(k, u, order) + [0]
        want = [spline._mp(lam[m] - lam[m + 1]) for m in range(k)]
        _agree("spline: k against section 21", out["k"], want)
        for m in range(k):
            at = [list(inputs["p"]) for _ in range(2)]
            at[0][3 * m : 3 * m + 3] = [1.0, 0.0, 0.0]
            at[1][3 * m : 3 * m + 3] = [0.0, 0.0, 0.0]
            hi, lo = (value({**inputs, "p": p})[key][0] for p in at)
            _agree(f"spline: k_{m} against d p / d c_{m}", [out["k"][m]], [hi - lo])

    return check


def gram(inputs: dict, out: dict) -> None:
    k, r, inv_h = int(inputs["k"]), int(inputs["r"]), Fraction(inputs["inv_h"])
    n, f = numerators(k), factorial(k - 1)
    exact = spline.gram_exact(k, r, inv_h)
    for s in range(k):
        for t in range(k):
            want = inv_h ** (2 * r - 1) * sum(
                (
                    Fraction(
                        n[s][p] * n[t][q] * factorial(p) * factorial(q),
                        factorial(p - r) * factorial(q - r),
                    )
                    / (p + q - 2 * r + 1)
                    for p in range(r, k)
                    for q in range(r, k)
                ),
                Fraction(0),
            ) / (f * f)
            if want != exact[s][t] or exact[s][t] != exact[t][s]:
                what = f"spline: Gram ({s}, {t}) {exact[s][t]} vs section 21 {want}"
                raise CrossCheckError(what)
        if sum(exact[s]) != 0:
            raise CrossCheckError(f"spline: Gram row {s} sums to {sum(exact[s])}")
    for i, g in enumerate(out["g"].data):
        e = exact[i % k][i // k]
        _agree(f"spline: Gram entry {i} against section 21", [g], [spline._mp(e)])
        stored = mpf(dec30(g))
        if to_f64(stored) != e.numerator / e.denominator or to_f32(stored) != to_f32(spline._mp(e)):
            raise CrossCheckError(f"spline: Gram entry {i}: 30 digits round away from {e}")
