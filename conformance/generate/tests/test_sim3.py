"""The Sim(3) ids (docs/decisions/0070): the references against the hat algebra, the checks against
planted errors, and the strata's draws."""

import unittest

from mpmath import mp, mpf

from gen import check, check_sim3, precision, sim3, so3
from gen.precision import is_binary32


class Sim3Test(unittest.TestCase):
    def setUp(self):
        precision.setup()

    def test_at_zero_sigma_and_rho_exp_is_the_rotation(self):
        out = sim3.exp({"tau": [0.3, -0.2, 0.5, 0.0, 0.0, 0.0, 0.0]})
        self.assertEqual(out["x"], [0, 0, 0])
        q = so3.exp_series([mpf(0.3), mpf(-0.2), mpf(0.5)])
        self.assertLess(max(abs(a - b) for a, b in zip(out["q"], q)), mpf(10) ** -100)

    def test_log_inverts_exp_and_the_two_jacobian_routes_agree(self):
        tau = [0.4, 0.1, -0.7, 1.5, -0.3, 0.8, -0.6]
        e = sim3.exp({"tau": tau})
        back = sim3.log({"q": e["q"], "x": e["x"], "sigma": tau[6]})["tau"]
        self.assertLess(max(abs(a - mpf(b)) for a, b in zip(back, tau)), mpf(10) ** -100)
        t = [mpf(c) for c in tau]
        for sign in (1, -1):
            a, b = sim3.series(t, sign), sim3.blocks(t, sign)
            self.assertLess(max(abs(x - y) for ra, rb in zip(a, b) for x, y in zip(ra, rb)), mpf(10) ** -100)

    def test_large_sigma_is_checked_by_the_identities(self):
        rec = {"tau": [0.6, 0.0, 0.8, 0.3, -1.2, 0.5, 700.0]}
        out = sim3.exp(rec)
        check_sim3.exp(rec, out)
        jr = sim3.jacobian("jr")(rec)
        check_sim3.jacobian("jr")(rec, jr)

    def test_planted_errors_fail_their_checks(self):
        rec = {"tau": [0.4, 0.1, -0.7, 1.5, -0.3, 0.8, -0.6]}
        out = sim3.exp(rec)
        out["x"][0] *= 1 + mpf(10) ** -50
        with self.assertRaises(check.CrossCheckError):
            check_sim3.exp(rec, out)
        jr = sim3.jacobian("jr")(rec)
        jr["J"].data[6 * 7 + 3] = -jr["J"].data[6 * 7 + 3]
        with self.assertRaises(check.CrossCheckError):
            check_sim3.jacobian("jr")(rec, jr)

    def test_the_f32_strata_hold_binary32_inputs_below_their_range(self):
        for s in sim3.TAU_STRATA:
            if not s.f32:
                continue
            for r in s.records()[:2]:
                self.assertTrue(all(is_binary32(x) for x in r["tau"]), s.name)
                self.assertLess(abs(r["tau"][6]), 88.72, s.name)


if __name__ == "__main__":
    unittest.main()
