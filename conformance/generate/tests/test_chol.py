"""`chol_n*` and `chol_solve_n*` against exact rational arithmetic."""

import unittest
from fractions import Fraction

from mpmath import mpf

from gen import check_linalg, chol, precision
from gen.dense import exact, ldlt_pivots
from gen.fmt import Mat
from gen.precision import is_binary32

precision.setup()
STRATA = {(n, solve): {s.name: s for s in chol.strata(n, solve)} for n in (3, 6) for solve in (0, 1)}


def solve_exact(A: list[list[Fraction]], b: list[Fraction]) -> list[Fraction]:
    """Gauss-Jordan in rationals: the exact A^-1 b."""
    n = len(A)
    M = [[*row, bi] for row, bi in zip(A, b, strict=True)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        M[c] = [x / M[c][c] for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                M[r] = [x - M[r][c] * y for x, y in zip(M[r], M[c], strict=True)]
    return [M[r][n] for r in range(n)]


class StrataTest(unittest.TestCase):
    def test_catalogue(self):
        names = list(STRATA[3, 0])
        self.assertEqual(
            names,
            [
                "chol:spd",
                "chol:cond-1e-4",
                "chol:cond-1e-8",
                "chol:cond-1e-12",
                "chol:diag-scale",
                "chol:indefinite",
                "chol:spd@f32",
                "chol:cond-1e-2@f32",
                "chol:cond-1e-4@f32",
                "chol:cond-1e-6@f32",
                "chol:diag-scale@f32",
                "chol:indefinite@f32",
            ],
        )
        self.assertEqual(list(STRATA[6, 1]), [n for n in names if "indefinite" not in n])

    def test_valid_is_decided_exactly(self):
        self.assertIsNone(ldlt_pivots([[Fraction(1), Fraction(2)], [Fraction(2), Fraction(4)]]))
        self.assertEqual(ldlt_pivots([[Fraction(4), Fraction(2)], [Fraction(2), Fraction(3)]]), [4, 2])
        for n in (3, 6):
            for name, s in STRATA[n, 0].items():
                for record in s.records()[:8]:
                    want = 0 if name.startswith("chol:indefinite") else 1
                    self.assertEqual(chol.chol(record)["valid"], [want], (n, name))

    def test_diag_scale_is_spd_scaled_and_binary32_at_f32(self):
        spd = STRATA[3, 0]["chol:spd"].records()[0]["A"]
        scaled = STRATA[3, 0]["chol:diag-scale"].records()[0]["A"]
        d = [mpf(scaled.data[4 * i]) / mpf(spd.data[4 * i]) for i in range(3)]
        for i in range(3):
            for j in range(3):  # S A S with S^2 = d: entries d_i^1/2 d_j^1/2 A_ij to rounding
                want = mpf(spd.data[3 * j + i]) * (d[i] * d[j]) ** 0.5
                self.assertLess(abs(mpf(scaled.data[3 * j + i]) / want - 1), 1e-15)
        for record in STRATA[6, 1]["chol:diag-scale@f32"].records()[:4]:
            self.assertTrue(all(is_binary32(x) for x in [*record["A"].data, *record["b"]]))


class ReferenceTest(unittest.TestCase):
    def test_chol_solve_is_the_exact_solution(self):
        """Equilibrated LU at 120 digits against Gauss-Jordan in rationals, at every scale."""
        for n in (3, 6):
            for name in ("chol:cond-1e-12", "chol:diag-scale", "chol:diag-scale@f32"):
                for record in STRATA[n, 1][name].records()[:2]:
                    x = chol.chol_solve(record)["x"]
                    want = solve_exact(exact(record["A"]), [Fraction(c) for c in record["b"]])
                    for got, w in zip(x, want, strict=True):
                        w = mpf(w.numerator) / w.denominator
                        self.assertLess(abs(got - w), mpf(10) ** -105 * abs(w), (n, name))


class CrossCheckTest(unittest.TestCase):
    def test_a_perturbation_beyond_100_digits_fails(self):
        for n in (3, 6):
            record = STRATA[n, 0]["chol:spd"].records()[0]
            good = chol.chol(record)
            check_linalg.chol(record, good)
            solve = STRATA[n, 1]["chol:spd"].records()[0]
            x = chol.chol_solve(solve)
            check_linalg.chol_solve(solve, x)
            for rel, fails in ((mpf(10) ** -95, True), (mpf(10) ** -105, False)):
                L = list(good["L"].data)
                L[0] *= 1 + rel
                xs = list(x["x"])
                xs[-1] *= 1 + rel
                for check, rec, bad in (
                    (check_linalg.chol, record, {**good, "L": Mat((n, n), L)}),
                    (check_linalg.chol_solve, solve, {"x": xs}),
                ):
                    if fails:
                        with self.assertRaises(check_linalg.CrossCheckError, msg=n):
                            check(rec, bad)
                    else:
                        check(rec, bad)

    def test_a_flipped_mask_fails(self):
        good = STRATA[3, 0]["chol:indefinite"].records()[0]
        with self.assertRaises(check_linalg.CrossCheckError):
            check_linalg.chol(good, {**chol.chol(good), "valid": [mpf(1)]})
        spd = STRATA[3, 0]["chol:spd"].records()[0]
        with self.assertRaises(check_linalg.CrossCheckError):
            check_linalg.chol(spd, {**chol.chol(spd), "valid": [mpf(0)]})


if __name__ == "__main__":
    unittest.main()
