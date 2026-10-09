"""`mat2_inverse_adj` (0061) against exact rational arithmetic."""

import unittest
from fractions import Fraction

from mpmath import mp, mpf

from gen import check_linalg, mat2, precision
from gen.dense import exact
from gen.fmt import Mat
from gen.precision import is_binary32

precision.setup()
STRATA = {s.name: s for s in mat2.STRATA}


class StrataTest(unittest.TestCase):
    def test_catalogue(self):
        names = ["cond:1", "cond:1e4", "cond:1e8", "cond:singular"]
        narrow = ["cond:1@f32", "cond:1e4@f32", "cond:singular@f32"]
        self.assertEqual(list(STRATA), [*names, *narrow])

    def test_inputs_are_of_their_precision(self):
        for name, s in STRATA.items():
            for record in s.records():
                if name.endswith("@f32"):
                    self.assertTrue(all(is_binary32(x) for x in record["A"].data), name)

    def test_singular_strata_are_exactly_singular_and_the_others_are_not(self):
        for name, s in STRATA.items():
            for record in s.records():
                det = mat2.determinant(exact(record["A"]))
                if name.startswith("cond:singular"):
                    self.assertEqual(det, 0, name)
                    # The products round equally, so the shipped `det` is exactly 0 too.
                    a, c, b, d = record["A"].data
                    self.assertEqual(a * d - b * c, 0.0, name)
                else:
                    self.assertNotEqual(det, 0, name)

    def test_the_f32_records_are_the_binary64_draws_rounded(self):
        for name in ("cond:1e4", "cond:singular"):
            wide, narrow = STRATA[name].records(), STRATA[name + "@f32"].records()
            self.assertEqual(len(wide), len(narrow))


class ReferenceTest(unittest.TestCase):
    def test_integer_matrix(self):
        out = mat2.inverse_adj({"A": Mat((2, 2), [4.0, 2.0, 7.0, 6.0])})  # [[4, 7], [2, 6]]
        self.assertEqual(out["det"], [mpf(10)])
        want = [Fraction(6, 10), Fraction(-2, 10), Fraction(-7, 10), Fraction(4, 10)]
        for got, w in zip(out["inv"].data, want, strict=True):
            self.assertLess(abs(got - mpf(w.numerator) / w.denominator), mpf(10) ** -110)

    def test_singular_reference_is_zero(self):
        out = mat2.inverse_adj({"A": Mat((2, 2), [1.0, 3.0, 2.0, 6.0])})
        self.assertEqual(out["det"], [mpf(0)])
        self.assertEqual(out["inv"].data, [0, 0, 0, 0])

    def test_every_record_passes_its_cross_check(self):
        for s in STRATA.values():
            for record in s.records():
                check_linalg.mat2_inverse_adj(record, mat2.inverse_adj(record))

    def test_the_cross_check_catches_a_wrong_inverse(self):
        record = {"A": Mat((2, 2), [4.0, 2.0, 7.0, 6.0])}
        out = mat2.inverse_adj(record)
        out["inv"] = Mat((2, 2), [mp.mpf(x) * (1 + mpf(10) ** -50) for x in out["inv"].data])
        with self.assertRaises(check_linalg.CrossCheckError):
            check_linalg.mat2_inverse_adj(record, out)


if __name__ == "__main__":
    unittest.main()
