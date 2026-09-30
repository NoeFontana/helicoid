import unittest
from fractions import Fraction
from math import factorial

from mpmath import mp, mpf

from gen import coeff, precision

precision.setup()


def series(theta: float, terms: int = 80) -> tuple[Fraction, Fraction]:
    """k and dk/d(theta^2), exactly, from k = 1/2 sum (-x/4)^n / (2n+1)!."""
    x = Fraction(theta) ** 2
    value = sum(Fraction(1, 2) * (-x / 4) ** n / factorial(2 * n + 1) for n in range(terms))
    deriv = sum(
        Fraction(1, 2) * n * Fraction(-1, 4) ** n * x ** (n - 1) / factorial(2 * n + 1)
        for n in range(1, terms)
    )
    return value, deriv


def close(a, b: Fraction, rel: Fraction) -> bool:
    m, e = precision.signed_man_exp(a)
    return abs(m * Fraction(2) ** e - b) <= rel * abs(b)


class CoeffKTest(unittest.TestCase):
    def test_removable_singularity(self):
        out = coeff.k({"theta": 0.0})
        self.assertEqual(out["value"], mpf(1) / 2)
        self.assertLess(abs(out["d_branch"] + mpf(1) / 48), mpf(10) ** -110)

    def test_matches_exact_series_to_100_digits(self):
        for theta in (1e-310, 1e-12, 1e-4, 0.3, 1.0, 2.5, 3.141592653589793 - 1e-3):
            out = coeff.k({"theta": theta})
            value, deriv = series(theta)
            self.assertTrue(close(out["value"], value, Fraction(1, 10**100)), theta)
            self.assertTrue(close(out["d_branch"], deriv, Fraction(1, 10**100)), theta)

    def test_a_complex_derivative_is_an_error_not_a_value(self):
        with self.assertRaises(ArithmeticError):
            coeff._with_branch_derivative(lambda x: mp.mpc(0, 1) * x, 0.5)
        out = coeff._with_branch_derivative(lambda x: x * x, 0.5)  # x = 1/4: f' = 2x
        self.assertEqual((out["value"], out["d_branch"]), (mpf(1) / 16, mpf(1) / 2))


if __name__ == "__main__":
    unittest.main()
