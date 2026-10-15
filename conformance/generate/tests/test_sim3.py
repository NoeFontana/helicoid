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

    def test_a_grid_cell_draws_sigma_and_theta_in_its_band(self):
        """0094 (draft): no two records of a cell share sigma or theta; each is within 2^-10 of
        the cell's offset; the `@f32` cell rounds the binary64 draw."""
        by = {s.name: s for s in sim3.TAU_STRATA}
        for name, s, t in (("sigma:1/theta:1e-4", 1, 1e-4), ("sigma:-700/theta:3", -700, 3)):
            taus = [r["tau"] for r in by[name].records()]
            sigmas = [r[6] for r in taus]
            thetas = [so3.norm([mpf(c) for c in r[:3]]) for r in taus]
            self.assertEqual(len(set(sigmas)), len(taus), name)
            self.assertTrue(all(abs(x - s) <= abs(s) * 2**-10 for x in sigmas), name)
            self.assertTrue(all(abs(x - t) <= t * (2**-10 + 1e-15) for x in thetas), name)
        near_pi = [so3.norm([mpf(c) for c in r["tau"][:3]]) for r in by["sigma:0/theta:pi-1e-6"].records()]
        self.assertTrue(all(abs(mp.pi - x - mpf(10) ** -6) <= 10**-6 * (2**-10 + 1e-9) for x in near_pi))
        self.assertTrue(all(r["tau"][6] == 0 for r in by["sigma:0/theta:1"].records()))
        f64 = [r["tau"][6] for r in by["sigma:1/theta:1e-4"].records()]
        f32 = [r["tau"][6] for r in by["sigma:1/theta:1e-4@f32"].records()]
        self.assertEqual(f32, [precision.to_f32(x) for x in f64])

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
