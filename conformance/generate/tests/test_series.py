import json
import unittest
from fractions import Fraction
from math import factorial
from unittest import mock

from mpmath import mp, mpf

from gen import precision, series

precision.setup()


class ExactSeriesTest(unittest.TestCase):
    def test_closed_forms_of_the_coefficients(self):
        """Textbook expansions (Taylor of cos, sin; Bernoulli numbers for c), not from the code."""
        exact = series.exact()
        for j in range(30):
            sign = (-1) ** j
            self.assertEqual(exact["a"][j], Fraction(sign, factorial(2 * j + 2)))
            self.assertEqual(exact["b"][j], Fraction(sign, factorial(2 * j + 3)))
            self.assertEqual(exact["d"][j], Fraction(sign, factorial(2 * j + 4)))
            self.assertEqual(exact["e"][j], Fraction(sign * (j + 1), factorial(2 * j + 5)))
            self.assertEqual(exact["k"][j], Fraction(sign, 2 * 4**j * factorial(2 * j + 1)))
            self.assertEqual(exact["r"][j], Fraction(sign, 2 * j + 1))
            self.assertEqual(exact["cos_half"][j], Fraction(sign, 4**j * factorial(2 * j)))
            bernoulli = series.rationalize(mp.bernoulli(2 * j + 2))
            self.assertEqual(exact["c"][j], sign * bernoulli / factorial(2 * j + 2))


class TaylorTest(unittest.TestCase):
    def test_mpmath_series_equal_the_exact_algebra_and_the_document(self):
        for name in series.NAMES:
            got = series.taylor(name)
            self.assertEqual(len(got), series.TERMS)
            self.assertEqual(got, series.exact()[name][: series.TERMS], name)

    def test_a_wrong_document_value_fails_generation(self):
        wrong = {**series.DOC_LEADING, "e": ("1/120", "-1/2520", "1/120960", "-1/9979201")}
        series.taylor.cache_clear()
        self.addCleanup(series.taylor.cache_clear)
        with mock.patch.dict(series.DOC_LEADING, wrong):
            with self.assertRaises(series.SeriesError):
                series.taylor("e")
            self.assertEqual(series.taylor("a")[0], Fraction(1, 2))

    def test_rationalize_reconstructs_or_refuses(self):
        self.assertEqual(series.rationalize(mpf(1) / 3), Fraction(1, 3))
        self.assertEqual(series.rationalize(-mpf(7) / 645120), Fraction(-7, 645120))
        with self.assertRaises(series.SeriesError):
            series.rationalize(mp.pi)
        with self.assertRaises(series.SeriesError):  # a rational off by 1e-90 is not identified
            series.rationalize(mpf(1) / 3 + mpf(10) ** -90)

    def test_record_format(self):
        data, n, checked = series.build()
        rows = [json.loads(line) for line in data.decode().splitlines()]
        self.assertEqual((n, checked, len(rows)), (8, 8, 8))
        self.assertEqual([r["coeff"] for r in rows], [*"kabcder", "cos_half"])
        self.assertEqual([r["id"] for r in rows], list(range(8)))
        for r in rows:
            self.assertEqual(len(r["series"]), series.TERMS)
            self.assertTrue(all("/" in s for s in r["series"]))
            self.assertEqual([Fraction(s) for s in r["series"]], list(series.taylor(r["coeff"])))
        self.assertEqual(
            [(r["branch"], r["prefactor"]) for r in rows[-2:]],
            [("n^2/w^2", "2/w"), ("theta^2", "1")],
        )

    def test_cos_half_is_held_to_the_exact_algebra_and_to_nothing_from_the_document(self):
        self.assertNotIn("cos_half", series.DOC_LEADING)
        self.assertEqual(series.taylor("cos_half")[:3], (1, Fraction(-1, 8), Fraction(1, 384)))
        series.taylor.cache_clear()
        self.addCleanup(series.taylor.cache_clear)
        wrong = {**series.exact(), "cos_half": (Fraction(1), Fraction(-1, 9))}
        with mock.patch.object(series, "exact", lambda: wrong):
            with self.assertRaises(series.SeriesError):
                series.taylor("cos_half")

if __name__ == "__main__":
    unittest.main()
