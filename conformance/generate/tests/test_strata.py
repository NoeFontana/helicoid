import unittest

from mpmath import mp, mpf

from gen import precision
from gen.strata import SCALAR_THETA_STRATA as STRATA

precision.setup()
BY_NAME = {s.name: s for s in STRATA}
NO_AXES = {"theta:exact0", "theta:dense"}


class StrataTest(unittest.TestCase):
    def test_catalogue(self):
        names = [s.name for s in STRATA]
        self.assertEqual(len(set(names)), len(names))
        self.assertEqual(names[:3], ["theta:1e-12", "theta:1e-11", "theta:1e-10"])
        self.assertEqual(
            names[11:15], ["theta:1e-1", "theta:1e0", "theta:exact0", "theta:subnormal"]
        )
        self.assertEqual(names[15], "theta:pi-1e-1")
        self.assertEqual(names[-2:], ["theta:pi-1e-12", "theta:dense"])
        self.assertEqual(sum(s.count for s in STRATA), 1710)
        for s in STRATA:
            self.assertEqual(len(s.thetas()), s.count, s.name)

    def test_decades(self):
        for e in range(-12, 0):
            xs = BY_NAME[f"theta:1e{e}"].thetas()
            self.assertTrue(all(mpf(10) ** e <= mpf(x) < mpf(10) ** (e + 1) for x in xs), e)
        xs = BY_NAME["theta:1e0"].thetas()
        self.assertTrue(all(1 <= x < float(mp.pi - mpf(1) / 10) for x in xs))

    def test_fixed_and_special(self):
        self.assertEqual(BY_NAME["theta:exact0"].thetas(), [0.0])
        sub = BY_NAME["theta:subnormal"].thetas()
        self.assertTrue(all(0 < x < 2.2250738585072014e-308 for x in sub))
        for k in range(1, 13):
            (t,) = BY_NAME[f"theta:pi-1e-{k}"].thetas()
            self.assertLess(mpf(t), mp.pi)
            self.assertLess(abs(mpf(t) - (mp.pi - mpf(10) ** -k)), mpf(2) ** -51)

    def test_dense_grid(self):
        xs = BY_NAME["theta:dense"].thetas()
        self.assertEqual((len(xs), xs[0], xs[-1]), (801, 1e-4, 1.0))
        self.assertEqual([xs[i] for i in (200, 400, 600)], [1e-3, 1e-2, 1e-1])
        self.assertTrue(all(a < b for a, b in zip(xs, xs[1:], strict=False)))

    def test_streams_are_per_stratum_and_reproducible(self):
        a = BY_NAME["theta:1e-8"]
        first = a.thetas()
        BY_NAME["theta:1e-7"].thetas()
        self.assertEqual(a.thetas(), first)
        self.assertNotEqual(a.thetas(seed=1), first)
        self.assertEqual(a.axes(4), a.axes(4))
        self.assertNotEqual(a.axes(4), BY_NAME["theta:1e-7"].axes(4))

    def test_axes_only_where_the_spec_lists_one(self):
        for s in STRATA:
            self.assertEqual(len(s.axes()), 0 if s.name in NO_AXES else 64, s.name)
        self.assertEqual(sum(s.carries_axes for s in STRATA), len(STRATA) - len(NO_AXES))

    def test_axes_known_answers(self):
        # Independent implementation (own splitmix64 and Fraction rounding), 120 digits.
        got = [tuple(c.hex() for c in v) for v in BY_NAME["theta:pi-1e-3"].axes(3)]
        self.assertEqual(
            got,
            [
                ("0x1.4843d0ad0eea0p-1", "0x1.884a99e66c7f0p-1", "-0x1.6379844f98520p-5"),
                ("-0x1.a197690d61960p-1", "-0x1.14556c753bc4fp-1", "-0x1.ab1a72ba599a0p-3"),
                ("0x1.32d0845ed5e0dp-3", "-0x1.d72ea40874e01p-1", "-0x1.721fadb80e448p-2"),
            ],
        )

    def test_draws_do_not_depend_on_the_ambient_precision(self):
        for name in ("theta:1e-8", "theta:1e0", "theta:pi-1e-5", "theta:dense"):
            s = BY_NAME[name]
            expected = (s.thetas(), s.axes(4))
            for dps in (15, 53, 400):
                with mp.workdps(dps):
                    self.assertEqual((s.thetas(), s.axes(4)), expected, (name, dps))
            self.assertEqual(mp.dps, precision.DPS)


if __name__ == "__main__":
    unittest.main()
