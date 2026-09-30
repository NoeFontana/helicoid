import unittest
from fractions import Fraction

from mpmath import mp, mpf

from gen import fmt, precision

precision.setup()


class HexTest(unittest.TestCase):
    def test_round_trip_including_subnormal_and_negative(self):
        for x in (
            0.0,
            -0.0,
            1.0,
            -0.1,
            5e-324,
            2.2250738585072009e-308,
            1.7976931348623157e308,
            3.141592653589793,
        ):
            self.assertEqual(float.fromhex(fmt.hex_float(x)).hex(), x.hex())

    def test_to_f64_rounds_to_nearest_even_and_keeps_sign(self):
        for x in (0.1, -0.1, 1e-310, -1e-310, 5e-324, 1e300, 0.0):
            self.assertEqual(precision.to_f64(mpf(x)), x)
        half_ulp_above_one = mpf(1) + mpf(2) ** -53  # tie: rounds to even, i.e. 1.0
        self.assertEqual(precision.to_f64(half_ulp_above_one), 1.0)
        self.assertEqual(precision.to_f64(half_ulp_above_one + mpf(2) ** -100), 1.0 + 2.0**-52)
        self.assertEqual(
            precision.to_f64(mpf(2) ** -1075), 0.0
        )  # half the smallest subnormal: ties to even
        self.assertEqual(precision.to_f64(mpf(2) ** -1075 * 1.5), 5e-324)


class Dec30Test(unittest.TestCase):
    def test_shape_and_known_values(self):
        self.assertEqual(fmt.dec30(mpf(0)), "0.00000000000000000000000000000e0")
        self.assertEqual(fmt.dec30(mpf(1) / 2), "5.00000000000000000000000000000e-1")
        self.assertEqual(fmt.dec30(-mpf(1) / 3), "-3.33333333333333333333333333333e-1")
        self.assertEqual(fmt.dec30(mpf(2) ** -1074), "4.94065645841246544176568792868e-324")
        self.assertEqual(fmt.dec30(mpf(1000)), "1.00000000000000000000000000000e3")

    def test_ties_round_half_even_and_carry(self):
        self.assertEqual(
            fmt.dec30(mpf(1000000000000000000000000000005)), "1.00000000000000000000000000000e30"
        )
        self.assertEqual(
            fmt.dec30(mpf(1000000000000000000000000000015)), "1.00000000000000000000000000002e30"
        )
        self.assertEqual(
            fmt.dec30(mpf(10**30) - mpf(1) / 2), "1.00000000000000000000000000000e30"
        )  # carry

    def test_round_trip_within_half_unit_of_30th_digit(self):
        for x in (
            mp.pi,
            -mp.e,
            mp.pi * mpf(10) ** -200,
            mpf(2) ** 300 / 7,
            mpf(1) - mpf(2) ** -100,
        ):
            s = fmt.dec30(x)
            mant, exp = s.split("e")
            self.assertEqual(len(mant.lstrip("-")), 31)
            back = Fraction(mant) * Fraction(10) ** int(exp)
            m, e = precision.signed_man_exp(x)
            exact = m * Fraction(2) ** e
            self.assertLessEqual(abs(back - exact), abs(exact) * Fraction(1, 10**29) / 2)


class DumpsTest(unittest.TestCase):
    def test_compact_sorted(self):
        self.assertEqual(
            fmt.dumps({"b": [1, "x"], "a": {"z": 0, "y": 1}}), '{"a":{"y":1,"z":0},"b":[1,"x"]}'
        )


if __name__ == "__main__":
    unittest.main()
