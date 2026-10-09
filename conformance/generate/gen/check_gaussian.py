"""Generation-time cross-check of the `gaussian_mahalanobis_*` records (0065), to 100 digits.

GG.10(a): d2 is side-invariant. The other side's residual, formed by its own `Log`, and the other
side's covariance, Ad_X Sigma Ad_X^T (right to left) or its inverse image (left to right), give
the same d2 in exact arithmetic, so a record's d2 is recomputed that way without `gaussian.
mahalanobis`'s own residual. `valid` is recomputed from the eigenvalues of the equilibrated Sigma,
as `check_linalg.chol` does.
"""

from collections.abc import Callable

from mpmath import mp, mpf

from . import check, gaussian, sen3, so3
from .chol import equilibrator
from .dense import from_mat
from .precision import DPS


def mahalanobis(n: int) -> Callable[[dict, dict], None]:
    def check_(inputs: dict, out: dict) -> None:
        (valid,) = out["valid"]
        with mp.workdps(2 * DPS):
            A = from_mat(inputs["Sigma"])
            D = equilibrator(A)
            smallest = min(mp.eigsy(D * A * D, eigvals_only=True))
            if valid != int(smallest > 0):
                raise check.CrossCheckError(
                    f"gaussian_mahalanobis: valid {valid}, smallest eigenvalue {smallest}"
                )
            if valid == 0:
                return
            left = bool(inputs["side"])
            R = so3.rot(so3.normalize([mpf(c) for c in inputs["q0"]]))
            Ad = mp.matrix(sen3.adjoint(R, [mpf(c) for c in inputs["x0"]], n))
            other = Ad**-1 * A * (Ad**-1).T if left else Ad * A * Ad.T
            delta = mp.matrix(gaussian.residual(inputs, n, not left))
            d2 = (delta.T * mp.lu_solve(other, delta))[0]
            (got,) = out["d2"]
            check._agree("gaussian_mahalanobis: the other side's d2", [got], [d2], max(1, abs(d2)))

    return check_
