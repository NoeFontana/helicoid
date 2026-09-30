import unittest

from mpmath import mp, mpf

from gen import check, coeff, precision, series
from gen.registry import FUNCTIONS

precision.setup()


def record(name: str, inputs: dict) -> dict:
    return coeff.r(inputs) if name == "r" else coeff.evaluator(name)(inputs)


class CrossCheckTest(unittest.TestCase):
    def test_correct_records_pass(self):
        thetas = (0.0, 1e-310, 1e-9, 1e-5, 0.7, 2.0, 3.141592653589793 - 1e-9)
        for name in "kabcde":
            for theta in thetas:
                check.coefficient(name)({"theta": theta}, record(name, {"theta": theta}))
        for inputs in (
            {"n": 0.0, "w": 1.0},
            {"n": 1e-9, "w": 1.0},
            {"n": 0.6, "w": 0.8},
            {"n": 1.0, "w": 0.0},
            {"n": 1.0, "w": 5e-13},
        ):
            check.coefficient("r")(inputs, record("r", inputs))

    def test_a_disagreement_beyond_100_digits_fails_in_either_output(self):
        cases = [
            ("e", {"theta": 0.7}),
            ("e", {"theta": 1e-5}),
            ("c", {"theta": 3.1}),
            ("r", {"n": 0.6, "w": 0.8}),
            ("r", {"n": 1e-9, "w": 1.0}),
        ]
        for name, inputs in cases:
            good = record(name, inputs)
            for key in ("value", "d_branch"):
                bad = {**good, key: good[key] * (1 + mpf(10) ** -95)}
                with self.assertRaises(check.CrossCheckError, msg=(name, inputs, key)):
                    check.coefficient(name)(inputs, bad)
                check.coefficient(name)(inputs, {**good, key: good[key] * (1 + mpf(10) ** -105)})

    def test_the_independent_check_covers_theta_the_taylor_check_cannot(self):
        inputs = {"theta": 2.0}  # theta^2 = 4 > SMALL: only the exact series applies
        good = record("e", inputs)
        with self.assertRaises(check.CrossCheckError):
            check.coefficient("e")(inputs, {**good, "d_branch": good["d_branch"] * 1.001})

    def test_the_taylor_branch_reaches_100_digits_at_its_upper_limit(self):
        """`SMALL` is derived from `series.TERMS`: the truncated series must still match."""
        with mp.workdps(2 * precision.DPS):
            for name in series.NAMES:
                short, long = check._mp(name, False), check._mp(name, True)
                got_want = zip(
                    check._sum(short, check.SMALL), check._sum(long, check.SMALL), strict=True
                )
                for got, want in got_want:
                    self.assertLess(abs(got - want), abs(want) * mpf(10) ** -check.DIGITS, name)


class RegistryTest(unittest.TestCase):
    def test_every_coefficient_id_carries_its_own_cross_check(self):
        for name, spec in FUNCTIONS.items():
            self.assertTrue(name.startswith("coeff_"), name)
            inputs = {"n": 0.6, "w": 0.8} if name == "coeff_r" else {"theta": 0.7}
            good = spec.evaluate(inputs)
            self.assertIsNotNone(spec.check, name)
            spec.check(inputs, good)  # not another coefficient's check
            for key in ("value", "d_branch"):
                with self.assertRaises(check.CrossCheckError, msg=(name, key)):
                    spec.check(inputs, {**good, key: good[key] * (1 + mpf(10) ** -95)})


if __name__ == "__main__":
    unittest.main()
