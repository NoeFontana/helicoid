"""Generation-time cross-checks of the 0056 linear-algebra ids, each by a property independent of
the route that computed the record, at twice the working precision, to 100 digits:

- `solve_cubic`: each root's residual to 1e-100 of the size of the terms it sums, and Vieta's
  three relations to 1e-100 of theirs; the real roots are the discriminant's count, a complex pair
  is conjugate, the order is (re, im); a zero discriminant's roots (planted double and triple) are
  roots exactly, in rationals.
- `eig3`: |A V - V Lambda| <= 1e-100 |A| and |V^T V - I| <= 1e-100, Frobenius; lambda ascending.
- `chol`: `valid` is the sign of the smallest eigenvalue (`mp.eigsy` of the equilibrated A, whose
  eigenvalues have A's signs); L is lower triangular with a
  positive diagonal and |L L^T - A| <= 1e-100 |A|, or zero where `valid` is 0.
- `chol_solve`: |A x - b| <= 1e-100 |A| |x|.
- `quat_renormalize`: |q'| = 1, q' parallel to q (every 2x2 minor of [q, q'] to 1e-100 |q|),
  and q'.q > 0.
"""

from mpmath import mp, mpf

from .chol import equilibrator
from .cubic import discriminant, fraction, horner, is_root
from .dense import frobenius, from_mat
from .precision import DPS

DIGITS = 100


class CrossCheckError(Exception):
    pass


def _within(what: str, err, scale) -> None:
    if not err <= mpf(10) ** -DIGITS * scale:
        raise CrossCheckError(f"{what}: {mp.nstr(err, 10)} vs {mp.nstr(scale, 10)}")


def solve_cubic(inputs: dict, out: dict) -> None:
    coeffs = [fraction(inputs[k]) for k in "abcd"]
    a, b, c, d = coeffs
    disc = discriminant(a, b, c, d)
    if list(zip(out["re"], out["im"])) != sorted(zip(out["re"], out["im"])):
        raise CrossCheckError(f"solve_cubic: roots not ordered by (re, im): {out}")
    real = [i for i in range(3) if out["im"][i] == 0]
    if len(real) != (1 if disc < 0 else 3):
        raise CrossCheckError(f"solve_cubic: {len(real)} real roots, discriminant {disc}")
    if disc < 0:
        p, q = (i for i in range(3) if i not in real)
        if out["re"][p] != out["re"][q] or out["im"][p] != -out["im"][q]:
            raise CrossCheckError(f"solve_cubic: the pair is not conjugate: {out}")
    if disc == 0 and not all(is_root(coeffs, fraction(out["re"][i])) for i in range(3)):
        raise CrossCheckError(f"solve_cubic: a multiple root that is not exact: {out}")
    with mp.workdps(2 * DPS):
        cs = [mpf(x) for x in coeffs]
        zs = [mp.mpc(re, im) for re, im in zip(out["re"], out["im"], strict=True)]
        for z in zs:
            terms = sum(abs(cs[k]) * abs(z) ** (3 - k) for k in range(4))
            _within(f"solve_cubic: residual at {mp.nstr(z, 20)}", abs(horner(cs, z)), terms)
        z1, z2, z3 = zs
        n1, n2, n3 = (abs(z) for z in zs)
        vieta = [
            (z1 + z2 + z3, -cs[1] / cs[0], n1 + n2 + n3),
            (z1 * z2 + z1 * z3 + z2 * z3, cs[2] / cs[0], n1 * n2 + n1 * n3 + n2 * n3),
            (z1 * z2 * z3, -cs[3] / cs[0], n1 * n2 * n3),
        ]
        for i, (got, want, scale) in enumerate(vieta):
            _within(f"solve_cubic: Vieta {i + 1}", abs(got - want), scale)


def eig3(inputs: dict, out: dict) -> None:
    lam = out["lambda"]
    if lam != sorted(lam):
        raise CrossCheckError(f"eig3: eigenvalues not ascending: {lam}")
    with mp.workdps(2 * DPS):
        A, V = from_mat(inputs["A"]), from_mat(out["V"])
        _within("eig3: A V - V Lambda", frobenius(A * V - V * mp.diag(lam)), frobenius(A))
        _within("eig3: V^T V - I", frobenius(V.T * V - mp.eye(3)), 1)


def chol(inputs: dict, out: dict) -> None:
    (valid,) = out["valid"]
    with mp.workdps(2 * DPS):
        A, L = from_mat(inputs["A"]), from_mat(out["L"])
        n = A.rows
        D = equilibrator(A)  # a congruence: the signs of the eigenvalues stay, the scale goes
        smallest = min(mp.eigsy(D * A * D, eigvals_only=True))
        if valid != int(smallest > 0):
            raise CrossCheckError(f"chol: valid {valid}, smallest eigenvalue {smallest}")
        if valid == 0:
            if any(L[i, j] != 0 for i in range(n) for j in range(n)):
                raise CrossCheckError("chol: L is not zero where valid is 0")
            return
        if any(L[i, j] != 0 for i in range(n) for j in range(i + 1, n)):
            raise CrossCheckError("chol: L is not lower triangular")
        if not all(L[i, i] > 0 for i in range(n)):
            raise CrossCheckError("chol: L has a diagonal entry that is not positive")
        _within("chol: L L^T - A", frobenius(L * L.T - A), frobenius(A))


def chol_solve(inputs: dict, out: dict) -> None:
    with mp.workdps(2 * DPS):
        A = from_mat(inputs["A"])
        b = mp.matrix([mpf(c) for c in inputs["b"]])
        x = mp.matrix(out["x"])
        _within("chol_solve: A x - b", frobenius(A * x - b), frobenius(A) * frobenius(x))


def quat_renormalize(inputs: dict, out: dict) -> None:
    with mp.workdps(2 * DPS):
        q, p = [mpf(c) for c in inputs["q"]], out["q"]
        norm_q = mp.sqrt(sum(c * c for c in q))
        _within("quat_renormalize: |q'| - 1", abs(mp.sqrt(sum(c * c for c in p)) - 1), 1)
        for i in range(4):
            for j in range(i + 1, 4):
                _within(f"quat_renormalize: minor {i}{j}", abs(p[i] * q[j] - p[j] * q[i]), norm_q)
        if not sum(a * b for a, b in zip(p, q, strict=True)) > 0:
            raise CrossCheckError("quat_renormalize: q' is not on q's side")
