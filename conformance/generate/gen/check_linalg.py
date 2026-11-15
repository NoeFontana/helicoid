"""Generation-time cross-checks of the 0056 linear-algebra ids, each by a property independent of
the route that computed the record, at twice the working precision, to 100 digits:

- `solve_cubic`: each root's residual to 1e-100 of the size of the terms it sums, and Vieta's
  three relations to 1e-100 of theirs; the real roots are the discriminant's count, a complex pair
  is conjugate, the order is (re, im); a zero discriminant's roots (planted double and triple) are
  roots exactly, in rationals.
- `eig3`: |A V - V Lambda| <= 1e-100 |A| and |V^T V - I| <= 1e-100, Frobenius; lambda ascending.
- `eig_sym_n*` (0111): `eig3`'s at n; `kappa` is 0 exactly where the stored A is not positive
  definite (an exact LDL^T), else at least 1 and within 1e-100 of lambda_max / lambda_min of the
  equilibrated A, recomputed by `mp.eig` (a route independent of `eigsy`).
- `svd3` (0101): |A V - U Sigma| <= 1e-100 |A|, |U^T U - I| and |V^T V - I| <= 1e-100, Frobenius;
  det U and det V within 1e-100 of 1; sigma1 >= sigma2 >= |sigma3|; sign sigma3 = sign det A, the
  determinant of the stored A in rationals.
- `nearest_rotation` (0102): |R^T R - I| <= 1e-100, Frobenius, and det R within 1e-100 of 1;
  no rotation of `_rivals` nearer M, Frobenius.
- `rq3` (0117): |M - U R| <= 1e-100 |M|, |R^T R - I| <= 1e-100, det R within 1e-100 of 1; U is
  zero below its diagonal and u22, u33 >= 0; `valid` is det M > 0 and u11 has the sign of det M
  where the factorization is unique; and where it is, `kappa_U` and `kappa_R` bound the change of
  (U, R) along a random direction, recomputed by a central difference of the construction at
  twice the working precision, and RQ.6's closed form agrees with that difference to 1e-30
  relative.
- `chol`: `valid` is the sign of the smallest eigenvalue (`mp.eigsy` of the equilibrated A, whose
  eigenvalues have A's signs); L is lower triangular with a
  positive diagonal and |L L^T - A| <= 1e-100 |A|, or zero where `valid` is 0.
- `chol_solve`: |A x - b| <= 1e-100 |A| |x|.
- `quat_renormalize`: |q'| = 1, q' parallel to q (every 2x2 minor of [q, q'] to 1e-100 |q|),
  and q'.q > 0.
- `mat2_inverse_adj` (0061): |det - mp.det(A)| <= 1e-100 |A|^2; where det is 0, inv is zero; else
  |A inv - I| <= 1e-100 |A| |inv|.
"""

from mpmath import mp, mpf

from . import so3
from .chol import equilibrator
from .cubic import discriminant, fraction, horner, is_root
from .dense import det3, exact, frobenius, from_mat, ldlt_pivots, unit_vector
from .precision import DPS
from .rng import SEED, stream, unit_quaternion_s3

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


def eig_sym(inputs: dict, out: dict) -> None:
    lam, n = out["lambda"], inputs["A"].shape[0]
    if lam != sorted(lam):
        raise CrossCheckError(f"eig_sym_n{n}: eigenvalues not ascending: {lam}")
    with mp.workdps(2 * DPS):
        A, V = from_mat(inputs["A"]), from_mat(out["V"])
        _within(f"eig_sym_n{n}: A V - V Lambda", frobenius(A * V - V * mp.diag(lam)), frobenius(A))
        _within(f"eig_sym_n{n}: V^T V - I", frobenius(V.T * V - mp.eye(n)), 1)
        (k,) = out["kappa"]
        if ldlt_pivots(exact(inputs["A"])) is None:
            if k != 0:
                raise CrossCheckError(f"eig_sym_n{n}: kappa {k} for an A that is not positive definite")
            return
        D = equilibrator(A)
        E = [mp.re(x) for x in mp.eig(D * A * D, left=False, right=False)]
        _within(f"eig_sym_n{n}: kappa", abs(k - max(E) / min(E)), k)


def _sign(x) -> int:
    return (x > 0) - (x < 0)


def svd3(inputs: dict, out: dict) -> None:
    s1, s2, s3 = out["sigma"]
    if not s1 >= s2 >= abs(s3):
        raise CrossCheckError(f"svd3: singular values out of order: {out['sigma']}")
    det = det3(exact(inputs["A"]))
    if _sign(s3) != _sign(det):
        raise CrossCheckError(f"svd3: sigma3 = {mp.nstr(s3, 10)} against det A = {float(det)}")
    with mp.workdps(2 * DPS):
        A, U, V = from_mat(inputs["A"]), from_mat(out["U"]), from_mat(out["V"])
        _within("svd3: A V - U Sigma", frobenius(A * V - U * mp.diag(out["sigma"])), frobenius(A))
        for name, X in (("U", U), ("V", V)):
            _within(f"svd3: {name}^T {name} - I", frobenius(X.T * X - mp.eye(3)), 1)
            _within(f"svd3: det {name} - 1", abs(mp.det(X) - 1), 1)


def rq3(inputs: dict, out: dict) -> None:
    from . import rq

    E = exact(inputs["M"])
    det = det3(E)
    if out["valid"][0] != int(det > 0):
        raise CrossCheckError(f"rq3: valid {out['valid'][0]} against det M = {float(det)}")
    with mp.workdps(2 * DPS):
        M, U, R = from_mat(inputs["M"]), from_mat(out["U"]), from_mat(out["R"])
        if any(U[i, j] != 0 for i in range(3) for j in range(i)):
            raise CrossCheckError("rq3: U is not upper triangular")
        if U[1, 1] < 0 or U[2, 2] < 0:
            raise CrossCheckError("rq3: a negative u22 or u33")
        _within("rq3: M - U R", frobenius(M - U * R), frobenius(M))
        _within("rq3: R^T R - I", frobenius(R.T * R - mp.eye(3)), 1)
        _within("rq3: det R - 1", abs(mp.det(R) - 1), 1)
        ku, kr = out["kappa_U"][0], out["kappa_R"][0]
        if ku == 0:
            return
        if _sign(U[0, 0]) != _sign(det):
            raise CrossCheckError(f"rq3: u11 = {mp.nstr(U[0, 0], 10)} against det M = {float(det)}")
        rng = stream(SEED, "check:rq3", repr(inputs["M"].data))
        # Each row of the direction at its row's scale: a step absolute in M would swamp the small
        # rows of `rq:row-spread`.
        norms = [mp.sqrt(sum(M[i, j] ** 2 for j in range(3))) or 1 for i in range(3)]
        D = mp.matrix([[norms[i] * unit_vector(rng, 3)[j] for j in range(3)] for i in range(3)])
        h = mpf(10) ** -(3 * DPS // 4)
        (Up, Rp), (Um, Rm) = (rq.factor(M + s * h * D, True) for s in (1, -1))
        dU, dR = (Up - Um) / (2 * h), (Rp - Rm) / (2 * h)
        norm_m = frobenius(M)
        bound = (1 + mpf(10) ** -30) * frobenius(D)
        if frobenius(dU) > ku * bound or frobenius(dR) * norm_m > kr * bound:
            raise CrossCheckError("rq3: a direction exceeds kappa")
        tU, tR = rq.tangent(U, R, D)
        # Relative: a step of 10^-90 at 240 digits leaves a difference of 10^-60 relative on the
        # worst-conditioned strata (`rq:diagonal`'s 10^60 spread); a wrong formula is O(1) off.
        for what, t, d in (("dU", tU, dU), ("dR", tR, dR)):
            if frobenius(t - d) > mpf(10) ** -30 * frobenius(d):
                raise CrossCheckError(f"rq3: {what} against RQ.6: {mp.nstr(frobenius(t - d), 5)}")


N_RIVALS = 25  # of each kind in `_rivals`


def _rivals(R):
    """50 rotations to hold `nearest_rotation`'s R against: 25 Haar, and R turned by 10^-k,
    k = 1..25, about random axes; one stream, the same for every record."""
    rng = stream(SEED, "nearest_rotation", "rivals")
    haar = [mp.matrix(so3.rot(list(unit_quaternion_s3(rng)))) for _ in range(N_RIVALS)]
    near = []
    for k in range(1, N_RIVALS + 1):
        half, axis = mpf(10) ** -k / 2, unit_vector(rng, 3)
        q = [mp.cos(half), *(mp.sin(half) * c for c in axis)]
        near.append(R * mp.matrix(so3.rot(q)))
    return haar + near


def nearest_rotation(inputs: dict, out: dict) -> None:
    with mp.workdps(2 * DPS):
        M, R = from_mat(inputs["M"]), from_mat(out["R"])
        _within("nearest_rotation: R^T R - I", frobenius(R.T * R - mp.eye(3)), 1)
        _within("nearest_rotation: det R - 1", abs(mp.det(R) - 1), 1)
        best = frobenius(M - R)
        for i, Q in enumerate(_rivals(R)):
            if frobenius(M - Q) < best:
                raise CrossCheckError(f"nearest_rotation: rival {i} is nearer M than R")


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


def mat2_inverse_adj(inputs: dict, out: dict) -> None:
    (det,) = out["det"]
    with mp.workdps(2 * DPS):
        A, inv = from_mat(inputs["A"]), from_mat(out["inv"])
        size = frobenius(A) ** 2
        _within("mat2_inverse_adj: det - mp.det(A)", abs(det - mp.det(A)), size)
        if det == 0:
            if any(inv[i, j] != 0 for i in range(2) for j in range(2)):
                raise CrossCheckError("mat2_inverse_adj: inv is not zero where det is 0")
            return
        _within(
            "mat2_inverse_adj: A inv - I",
            frobenius(A * inv - mp.eye(2)),
            frobenius(A) * frobenius(inv),
        )
