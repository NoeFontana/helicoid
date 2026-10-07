"""`quat_renormalize`: its strata stay inside NUMERICS section 3.6's band; the reference is unit."""

import unittest

from mpmath import mp, mpf

from gen import check_linalg, precision, quat
from gen.precision import is_binary32

precision.setup()
BY_NAME = {s.name: s for s in quat.STRATA}


class StrataTest(unittest.TestCase):
    def test_catalogue(self):
        self.assertEqual(
            list(BY_NAME),
            [
                *(f"renorm:eta-2^-{k}" for k in (27, 30, 40, 52)),
                "renorm:eta-edge",
                *(f"renorm:eta-2^-{k}@f32" for k in (12, 16, 20, 23)),
                "renorm:eta-edge@f32",
            ],
        )

    def test_every_stored_q_is_inside_the_band(self):
        for s in quat.STRATA:
            band = quat.band(s.f32)
            for record in s.records():
                eta = quat.eta(record["q"])
                self.assertLessEqual(abs(mpf(eta.numerator) / eta.denominator), band, s.name)
                if s.f32:
                    self.assertTrue(all(is_binary32(c) for c in record["q"]))

    def test_eta_is_the_stratum_s_with_the_sign_alternating(self):
        """Where rounding (2^-52 of |q|^2, or 2^-23) is far below 2^-k, eta is +-2^-k."""
        for name, k, tol in (("renorm:eta-2^-40", 40, 1e-3), ("renorm:eta-2^-16@f32", 16, 0.05)):
            for i, record in enumerate(BY_NAME[name].records()):
                eta = float(quat.eta(record["q"]))
                self.assertLess(abs(eta / ((-1) ** i * 2.0**-k) - 1), tol, (name, i))
        edge = [float(quat.eta(r["q"])) for r in BY_NAME["renorm:eta-edge"].records()]
        self.assertTrue(all(abs(abs(e) / 2**-26.29 - 1) < 1e-7 for e in edge))


class ReferenceTest(unittest.TestCase):
    def test_the_projection_is_unit_and_a_perturbation_beyond_100_digits_fails(self):
        record = BY_NAME["renorm:eta-edge"].records()[0]
        good = quat.renormalize(record)
        with mp.workdps(200):
            self.assertLess(abs(mp.sqrt(sum(c * c for c in good["q"])) - 1), mpf(10) ** -115)
        check_linalg.quat_renormalize(record, good)
        for i in range(4):
            for rel, fails in ((mpf(10) ** -95, True), (mpf(10) ** -105, False)):
                q = list(good["q"])
                q[i] *= 1 + rel
                if fails:
                    with self.assertRaises(check_linalg.CrossCheckError, msg=i):
                        check_linalg.quat_renormalize(record, {"q": q})
                else:
                    check_linalg.quat_renormalize(record, {"q": q})
        with self.assertRaises(check_linalg.CrossCheckError):  # -q' is unit and parallel
            check_linalg.quat_renormalize(record, {"q": [-c for c in good["q"]]})


if __name__ == "__main__":
    unittest.main()
