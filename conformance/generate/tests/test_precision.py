import unittest
from unittest import mock

import mpmath
from mpmath import mp

from gen import precision


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


if __name__ == "__main__":
    unittest.main()
