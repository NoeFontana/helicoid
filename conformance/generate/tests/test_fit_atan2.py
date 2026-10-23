import struct
import unittest
from fractions import Fraction

from mpmath import mp, mpf

from fit import atan2 as atanfit
from gen.precision import setup


def value(bits_hex: str) -> float:
    return struct.unpack(">d", bytes.fromhex(bits_hex[2:]))[0]


class AtanFitTest(unittest.TestCase):
    def setUp(self):
        setup()

    def test_the_series_below_the_switch_is_the_quotient(self):
        s = mpf(2) ** -20 * (1 - mpf(2) ** -40)
        series = atanfit.h(s)
        with mp.workdps(4 * mp.dps):
            r = mp.sqrt(s)
            direct = (mp.atan(r) / r - 1) / s
        self.assertLess(abs(series - direct), mpf(2) ** -300)

    def test_the_fit_levels_and_shortens_with_terms(self):
        top = mpf(1) / 4
        levels = []
        for terms in (4, 5, 6):
            c, level = atanfit.remez(terms, top)
            self.assertAlmostEqual(float(c[0]), -1 / 3, places=4)
            self.assertLessEqual(atanfit.weighted_error(c, top), level * (1 + mpf(2) ** -20))
            levels.append(level)
        self.assertTrue(levels[0] > levels[1] > levels[2])

    def test_pi_over_4_splits_with_two_spare_bits(self):
        # A binary32 sits in the top 24 bits of the binary64 significand, so its last two bits are
        # binary64 bits 29 and 30.
        for precision, p, spare in (("f64", 53, 0b11), ("f32", 24, 0b11 << 29)):
            fit = atanfit.fit(precision, 4, Fraction(1, 2))
            self.assertEqual(int(str(fit["pi4_hi"]), 16) & spare, 0)
            total = mpf(value(str(fit["pi4_hi"]))) + mpf(value(str(fit["pi4_lo"])))
            self.assertLess(abs(total - mp.pi / 4), mpf(2) ** (-2 * p + 2))


if __name__ == "__main__":
    unittest.main()
