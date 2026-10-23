import struct
import unittest

from mpmath import mp, mpf

from fit import atan2 as atanfit
from gen.precision import setup, to_f64


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

    def test_rounding_one_at_a_time_keeps_what_naive_rounding_loses(self):
        top = mpf(1) / 4
        exact, level = atanfit.remez(7, top)
        naive = atanfit.weighted_error([mpf(to_f64(x)) for x in exact], top)
        _, refit = atanfit.rounded(7, top, to_f64)
        self.assertLessEqual(refit, naive * (1 + mpf(2) ** -20))
        self.assertLess(refit, level * 2)

    def test_the_centres_tables_and_pi_over_4_split(self):
        fit = atanfit.binary64(4, 4)
        self.assertEqual(int(str(fit["pi4_hi"]), 16) & 0b11, 0)
        total = mpf(value(str(fit["pi4_hi"]))) + mpf(value(str(fit["pi4_lo"])))
        self.assertLess(abs(total - mp.pi / 4), mpf(2) ** -104)
        for j, (h, l) in enumerate(zip(fit["atan_hi"], fit["atan_lo"])):
            exact = mp.atan(mpf(j) / 4)
            self.assertLess(abs(mpf(value(h)) + mpf(value(l)) - exact), mpf(2) ** -104)
        with self.assertRaises(ValueError):
            atanfit.binary64(6, 4)


if __name__ == "__main__":
    unittest.main()
