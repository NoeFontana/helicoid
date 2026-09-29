import random
import struct
import unittest
from unittest import mock

import mpmath
from mpmath import mp, mpf

from gen import precision


def platform_f32(x: float) -> float:
    """The C cast to binary32 and back: nearest-even on every IEEE platform, and independent of
    `precision.to_f32`, which does it in integers."""
    return struct.unpack("<f", struct.pack("<f", x))[0]


class SetupTest(unittest.TestCase):
    def test_sets_the_working_precision(self):
        with mp.workdps(15):
            precision.setup()
            self.assertEqual(mp.dps, precision.DPS)

    def test_refuses_a_non_python_backend(self):
        with (
            mock.patch.object(mpmath.libmp, "BACKEND", "gmpy"),
            self.assertRaisesRegex(RuntimeError, "gmpy"),
        ):
            precision.setup()


class ToF32Test(unittest.TestCase):
    def test_known_values(self):
        ulp = 2.0**-23
        cases = {
            0.1: float.fromhex("0x1.99999ap-4"),
            1 + ulp / 2: 1.0,  # a tie: 1.0 has the even significand
            1 + 3 * ulp / 2: 1 + 2 * ulp,  # a tie: 1 + ulp is odd
            1 + ulp / 2 + 2.0**-40: 1 + ulp,  # just above the tie
            2.0**-149: 2.0**-149,  # the smallest subnormal
            2.0**-150: 0.0,  # a tie between 0 and 2^-149: 0 is even
            3 * 2.0**-150: 2.0**-148,  # a tie between 2^-149 and 2^-148: the second is even
            1e-310: 0.0,  # a binary64 subnormal, beneath a binary32
            1e-40: 71362 * 2.0**-149,  # 1e-40 / 2^-149 = 71362.4
            3.4028234663852886e38: 3.4028234663852886e38,  # the largest binary32
        }
        for x, want in cases.items():
            self.assertEqual(precision.to_f32(x), want, x.hex())
            self.assertEqual(precision.to_f32(-x), -want, x.hex())
        with mp.workdps(120):  # an mpf is rounded once, from its exact value
            self.assertEqual(precision.to_f32(mpf(1) + mpf(2) ** -24 + mpf(2) ** -100), 1 + ulp)

    def test_beyond_the_range_is_an_error_not_an_infinity(self):
        for x in (2.0**128, (2**24 - 0.5) * 2.0**104, 1e300, -1e39):
            with self.assertRaises(OverflowError, msg=x.hex()):
                precision.to_f32(x)

    def test_agrees_with_the_platform_cast_on_random_and_tie_inputs(self):
        rng = random.Random(0x32)
        xs = [
            rng.choice((-1, 1)) * rng.uniform(1, 2) * 2.0 ** rng.randint(-160, 127)
            for _ in range(20000)  # every binade of a binary32, its subnormals, and below
        ]
        # Exact ties: an odd multiple of half a unit in the last place.
        for _ in range(2000):
            xs.append((2 * rng.randrange(1, 2**23) + 1) * 2.0 ** (rng.randint(-149, 100) - 1))
        for x in xs:
            self.assertEqual(precision.to_f32(x), platform_f32(x), x.hex())

    def test_is_binary32_refuses_what_a_binary32_cannot_hold(self):
        exact = (0.0, 1.0, 2.0**-149, 3.4028234663852886e38)
        self.assertTrue(all(map(precision.is_binary32, exact)))
        self.assertFalse(any(map(precision.is_binary32, (0.1, 1 + 2.0**-24, 1e-310, 1e39))))


if __name__ == "__main__":
    unittest.main()
