"""The Sim(3) coefficient files (docs/decisions/0070): the constants, the reference's two routes,
and the rays."""

import unittest
from fractions import Fraction
from math import factorial

from mpmath import mp, mpf

from gen import precision, sim3coeff
from gen.precision import to_f32


class Sim3CoeffTest(unittest.TestCase):
    def setUp(self):
        precision.setup()

    def test_the_constants_are_one_over_n_factorial_rounded_once(self):
        data, n, _ = sim3coeff.build_series()
        self.assertEqual(n, sim3coeff.N_FACT + 1)
        for i, line in enumerate(data.decode().splitlines()):
            f64 = float.fromhex(line.split('"f64":"')[1].split('"')[0])
            f32 = float.fromhex(line.split('"f32":"')[1].split('"')[0])
            self.assertEqual(f64, Fraction(1, factorial(i)).__float__())
            self.assertEqual(f32, to_f32(mpf(1) / factorial(i)))

    def test_the_series_is_the_limits_at_the_origin(self):
        c = sim3coeff.scalars(mpf(0), mpf(0))
        want = (1, mpf(1) / 2, mpf(1) / 6, mpf(1) / 2, mpf(1) / 6, mpf(1) / 24, -mpf(1) / 24, -mpf(1) / 120)
        for name, w in zip(sim3coeff.NAMES, want, strict=True):
            self.assertLess(abs(c[name] - w), mpf(10) ** -100, name)

    def test_the_two_routes_agree_off_the_axes(self):
        for s, z in ((0.3, 0.7), (-2.5, 1e-3), (1.7, 4.0), (-1e-3, 2.0)):
            with mp.workdps(200):
                a, b = sim3coeff.scalars(mpf(s), mpf(z)), sim3coeff.closed(mpf(s), mpf(z))
                for name in sim3coeff.NAMES:
                    self.assertLess(abs(a[name] - b[name]), mpf(10) ** -150 * abs(b[name]), name)

    def test_the_rays_cover_both_signs_and_the_axes(self):
        pts = sim3coeff.points(False)
        self.assertEqual(len(pts), len(set(pts)))
        self.assertTrue(any(s < 0 and z == 0 for s, z in pts))
        self.assertTrue(any(s > 0 and z == 0 for s, z in pts))
        self.assertTrue(any(s == 0 and z > 0 for s, z in pts))
        self.assertTrue(all(z >= 0 for _, z in pts))

    def test_a_planted_error_fails_the_check(self):
        ref = sim3coeff.scalars(mpf(0.5), mpf(0.25))
        ref["dv2"] *= 1 + mpf(10) ** -50
        with self.assertRaises(sim3coeff.Sim3RefError):
            sim3coeff._check(0.5, 0.25, ref)


if __name__ == "__main__":
    unittest.main()
