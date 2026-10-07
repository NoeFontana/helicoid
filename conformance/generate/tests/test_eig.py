"""`eig3` against what its strata draw and the characteristic polynomial in rationals."""

import unittest
from fractions import Fraction

from mpmath import mp, mpf

from gen import check_linalg, eig, precision
from gen.dense import exact
from gen.fmt import Mat
from gen.precision import is_binary32

precision.setup()
BY_NAME = {s.name: s for s in eig.STRATA}


def charpoly(A: list[list[Fraction]]) -> list[Fraction]:
    """det(lambda I - A) = lambda^3 + c2 lambda^2 + c1 lambda + c0, exactly."""
    tr = A[0][0] + A[1][1] + A[2][2]
    minors = sum(A[i][i] * A[j][j] - A[i][j] * A[j][i] for i in range(3) for j in range(i + 1, 3))
    det = (
        A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
        - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
        + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])
    )
    return [Fraction(1), -tr, minors, -det]


class StrataTest(unittest.TestCase):
    def test_stored_matrices_are_exactly_symmetric_and_binary32_at_f32(self):
        for s in eig.STRATA:
            for record in s.records()[:4]:
                A = exact(record["A"])
                self.assertTrue(all(A[i][j] == A[j][i] for i in range(3) for j in range(3)))
                if s.f32:
                    self.assertTrue(all(is_binary32(x) for x in record["A"].data), s.name)

    def test_the_gap_strata_hold_their_relative_gap(self):
        """The close pair is 10^-k |lambda|_2 apart, to the precision's rounding of A."""
        for name, k, tol in (("eig:gap-1e-12/bottom", 12, 1e-2), ("eig:gap-1e-6/top@f32", 6, 0.2)):
            for record in BY_NAME[name].records()[:8]:
                lam = eig.eig3(record)["lambda"]
                close = lam[1] - lam[0] if name.split("/")[1].startswith("bottom") else lam[2] - lam[1]
                rel = close / mp.sqrt(sum(x * x for x in lam)) / mpf(10) ** -k
                self.assertLess(abs(rel - 1), tol, (name, rel))

    def test_triple_is_c_times_the_identity(self):
        for record in BY_NAME["eig:triple@f32"].records()[:4]:
            c = record["A"].data[0]
            self.assertEqual(record["A"], Mat((3, 3), [c, 0.0, 0.0, 0.0, c, 0.0, 0.0, 0.0, c]))
            out = eig.eig3(record)
            self.assertEqual(out["lambda"], [c] * 3)

    def test_scale_strata_are_eig_random_scaled(self):
        for name, f32, s in (("eig:scale-up", False, 70), ("eig:scale-down@f32", True, -8)):
            base = BY_NAME["eig:random@f32" if f32 else "eig:random"].records()[0]["A"].data
            scaled = BY_NAME[name].records()[0]["A"].data
            for a, b in zip(base, scaled, strict=True):
                self.assertLess(abs(mpf(b) / mpf(10) ** s / mpf(a) - 1), 2e-7 if f32 else 3e-16)


class ReferenceTest(unittest.TestCase):
    def test_eigenvalues_are_roots_of_the_exact_characteristic_polynomial(self):
        """To 1e-110 |A|, by one Newton step on the exact polynomial (`eigsy`'s accuracy is
        norm-wise: `eig:rank1`'s 1e-19 eigenvalues carry 1e-124 absolute)."""
        for name in ("eig:random", "eig:gap-1e-9/top", "eig:rank1", "eig:scale-down"):
            for record in BY_NAME[name].records()[:4]:
                cs = [mpf(c) for c in charpoly(exact(record["A"]))]
                norm = mp.sqrt(sum(mpf(x) ** 2 for x in record["A"].data))
                for lam in eig.eig3(record)["lambda"]:
                    value = ((lam + cs[1]) * lam + cs[2]) * lam + cs[3]
                    slope = (3 * lam + 2 * cs[1]) * lam + cs[2]
                    self.assertLess(abs(value / slope), mpf(10) ** -110 * norm, (name, lam))

    def test_eigenvectors_are_signed_by_their_largest_component(self):
        for record in BY_NAME["eig:random"].records()[:8]:
            V = eig.eig3(record)["V"].data
            for j in range(3):
                col = V[3 * j : 3 * j + 3]
                self.assertGreater(max(col, key=abs), 0)


class CrossCheckTest(unittest.TestCase):
    def test_a_perturbation_beyond_100_digits_fails(self):
        record = BY_NAME["eig:random"].records()[0]
        good = eig.eig3(record)
        check_linalg.eig3(record, good)
        for rel, fails in ((mpf(10) ** -95, True), (mpf(10) ** -105, False)):
            for i in range(3):
                lam = list(good["lambda"])
                lam[i] *= 1 + rel
                V = list(good["V"].data)
                V[4 * i] *= 1 + rel
                for bad in ({**good, "lambda": lam}, {**good, "V": Mat((3, 3), V)}):
                    if fails:
                        with self.assertRaises(check_linalg.CrossCheckError, msg=(i, bad)):
                            check_linalg.eig3(record, bad)
                    else:
                        check_linalg.eig3(record, bad)


if __name__ == "__main__":
    unittest.main()
