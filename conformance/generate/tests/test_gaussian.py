"""The `gaussian_mahalanobis_*` ids (0065) against hand-computed cases, planted errors for the
side-invariance cross-check of `gen.check_gaussian`, and the input draw."""

import unittest

from mpmath import mp, mpf

from gen import check_gaussian, gaussian, precision, so3
from gen.check import CrossCheckError
from gen.fmt import Mat
from gen.registry import FUNCTIONS

precision.setup()
TOL = mpf(10) ** -100
IDENTITY = [1.0, 0.0, 0.0, 0.0]


def quarter_z():
    return [float(mp.cos(mp.pi / 4)), 0.0, 0.0, float(mp.sin(mp.pi / 4))]


def diag(values):
    n = len(values)
    return Mat((n, n), [values[i] if i == j else 0.0 for j in range(n) for i in range(n)])


class HandComputedTest(unittest.TestCase):
    def test_a_pure_translation_is_whitened_by_its_own_variance(self):
        # Y = X Exp([0; rho]) at X = identity: delta = rho, d2 = sum rho_i^2 / s_i.
        rec = {
            "q0": IDENTITY,
            "x0": [0.0, 0.0, 0.0],
            "q1": IDENTITY,
            "x1": [1.0, 2.0, 0.0],
            "side": 0.0,
            "Sigma": diag([1.0, 1.0, 1.0, 4.0, 0.25, 1.0]),
        }
        (d2,) = gaussian.mahalanobis(1)(rec)["d2"]
        self.assertLess(abs(d2 - (mpf(1) / 4 + mpf(4) / mpf("0.25"))), TOL)

    def test_the_left_residual_is_in_the_world_frame(self):
        # X = (Rz(90), 0), Y = (Rz(90), (1, 0, 0)): the right residual is R^T t = (0, -1, 0), the
        # left one t = (1, 0, 0).
        rec = {
            "q0": quarter_z(),
            "x0": [0.0, 0.0, 0.0],
            "q1": quarter_z(),
            "x1": [1.0, 0.0, 0.0],
            "Sigma": diag([1.0, 1.0, 1.0, 1.0, 4.0, 1.0]),
        }
        right = gaussian.mahalanobis(1)({**rec, "side": 0.0})["d2"][0]
        left = gaussian.mahalanobis(1)({**rec, "side": 1.0})["d2"][0]
        self.assertLess(abs(right - mpf(1) / 4), 1e-30)  # rho_y = -1 against a variance of 4
        self.assertLess(abs(left - 1), 1e-30)  # rho_x = 1 against a variance of 1

    def test_an_indefinite_sigma_is_not_valid_and_d2_is_zero(self):
        rec = {
            "q0": IDENTITY,
            "x0": [0.0, 0.0, 0.0],
            "q1": IDENTITY,
            "x1": [1.0, 0.0, 0.0],
            "side": 0.0,
            "Sigma": diag([1.0, 1.0, 1.0, 1.0, 1.0, -1.0]),
        }
        out = gaussian.mahalanobis(1)(rec)
        self.assertEqual((out["valid"], out["d2"]), ([0], [0]))
        check_gaussian.mahalanobis(1)(rec, out)


class CrossCheckTest(unittest.TestCase):
    def test_a_wrong_answer_fails_the_other_side(self):
        for name in ("se3", "se23"):
            spec = FUNCTIONS[f"gaussian_mahalanobis_{name}"]
            for stratum in spec.strata[:2]:
                rec = spec.inputs(stratum)[0]
                out = spec.evaluate(rec)
                spec.check(rec, out)
                with self.assertRaises(CrossCheckError):
                    spec.check(rec, {**out, "d2": [out["d2"][0] * (1 + mpf(10) ** -50)]})
                with self.assertRaises(CrossCheckError):
                    spec.check(rec, {**out, "valid": [mpf(0)]})


class InputTest(unittest.TestCase):
    def test_every_stratum_holds_its_side_and_shapes(self):
        for name, n in (("se3", 1), ("se23", 2)):
            spec = FUNCTIONS[f"gaussian_mahalanobis_{name}"]
            names = [s.name for s in spec.strata]
            self.assertIn("gauss:lever", names)
            self.assertNotIn("gauss:lever@f32", names)
            for stratum in spec.strata:
                rec = spec.inputs(stratum)[0]
                self.assertEqual(len(rec["x0"]), 3 * n)
                self.assertEqual(rec["Sigma"].shape, (3 + 3 * n, 3 + 3 * n))
                left = stratum.name.split("@")[0] in ("gauss:left", "gauss:lever")
                self.assertEqual(rec["side"], float(left), stratum.name)

    def test_the_lever_stratum_has_its_translation_at_a_kilometre(self):
        spec = FUNCTIONS["gaussian_mahalanobis_se3"]
        stratum = next(s for s in spec.strata if s.name == "gauss:lever")
        for rec in spec.inputs(stratum)[:4]:
            self.assertAlmostEqual(float(so3.norm([mpf(c) for c in rec["x0"]])), 1000.0, places=9)

    def test_the_near_cut_residual_turns_between_two_and_a_half_and_three(self):
        spec = FUNCTIONS["gaussian_mahalanobis_se3"]
        stratum = next(s for s in spec.strata if s.name == "gauss:near-cut")
        for rec in spec.inputs(stratum)[:4]:
            phi = gaussian.residual(rec, 1, False)[:3]
            self.assertTrue(mpf("2.49") < so3.norm(phi) < mpf("3.01"))


if __name__ == "__main__":
    unittest.main()
