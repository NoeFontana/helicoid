import unittest
from fractions import Fraction
from math import factorial

from mpmath import mp, mpf

from gen import coeff, precision
from gen import series as gen_series
from gen.strata import COEFF_R_STRATA, R_STRATA

precision.setup()


def series(theta: float, terms: int = 80) -> tuple[Fraction, Fraction]:
    """k and dk/d(theta^2), exactly, from k = 1/2 sum (-x/4)^n / (2n+1)!."""
    x = Fraction(theta) ** 2
    value = sum(Fraction(1, 2) * (-x / 4) ** n / factorial(2 * n + 1) for n in range(terms))
    deriv = sum(
        Fraction(1, 2) * n * Fraction(-1, 4) ** n * x ** (n - 1) / factorial(2 * n + 1)
        for n in range(1, terms)
    )
    return value, deriv


def close(a, b: Fraction, rel: Fraction) -> bool:
    m, e = precision.signed_man_exp(a)
    return abs(m * Fraction(2) ** e - b) <= rel * abs(b)


class CoeffKTest(unittest.TestCase):
    def test_removable_singularity(self):
        out = coeff.k({"theta": 0.0})
        self.assertEqual(out["value"], mpf(1) / 2)
        self.assertLess(abs(out["d_branch"] + mpf(1) / 48), mpf(10) ** -110)

    def test_matches_exact_series_to_100_digits(self):
        for theta in (1e-310, 1e-12, 1e-4, 0.3, 1.0, 2.5, 3.141592653589793 - 1e-3):
            out = coeff.k({"theta": theta})
            value, deriv = series(theta)
            self.assertTrue(close(out["value"], value, Fraction(1, 10**100)), theta)
            self.assertTrue(close(out["d_branch"], deriv, Fraction(1, 10**100)), theta)

    def test_a_complex_derivative_is_an_error_not_a_value(self):
        with self.assertRaises(ArithmeticError):
            coeff._with_branch_derivative(lambda x: mp.mpc(0, 1) * x, 0.5)
        out = coeff._with_branch_derivative(lambda x: x * x, 0.5)  # x = 1/4: f' = 2x
        self.assertEqual((out["value"], out["d_branch"]), (mpf(1) / 16, mpf(1) / 2))


# Local copies of the definitions, so the tests do not share them with `coeff.DEFINITIONS`.
RAW = {
    "a": lambda t: (1 - mp.cos(t)) / t**2,
    "b": lambda t: (t - mp.sin(t)) / t**3,
    "c": lambda t: 1 / t**2 - (1 + mp.cos(t)) / (2 * t * mp.sin(t)),
    "d": lambda t: (t**2 + 2 * mp.cos(t) - 2) / (2 * t**4),
    "e": lambda t: (2 * t - 3 * mp.sin(t) + t * mp.cos(t)) / (2 * t**5),
}
REL = mpf(10) ** -100


class CoefficientTest(unittest.TestCase):
    def test_matches_the_raw_definition_at_high_precision(self):
        # 1200 digits absorb the cancellation of every definition down to theta = 1e-12.
        for name, g in RAW.items():
            for theta in (1e-12, 1e-4, 0.3, 1.0, 2.5, 3.141592653589793 - 1e-12):
                out = coeff.evaluator(name)({"theta": theta})
                with mp.workdps(1200):
                    x = mpf(theta) ** 2
                    value, deriv = (
                        g(mp.sqrt(x)),
                        mp.diff(lambda x, g=g: g(mp.sqrt(x)), x),
                    )
                self.assertLess(abs(out["value"] - value), REL * abs(value), (name, theta))
                self.assertLess(abs(out["d_branch"] - deriv), REL * abs(deriv), (name, theta))

    def test_zero_and_subnormal_take_the_analytic_extension(self):
        exact = gen_series.exact()
        for name in "kabcde":
            for theta in (0.0, 1e-310):
                out = coeff.evaluator(name)({"theta": theta})
                for key, want in zip(("value", "d_branch"), exact[name][:2], strict=True):
                    self.assertTrue(close(out[key], want, Fraction(1, 10**100)), (name, theta, key))

    def test_stable_grows_its_guard_until_the_evaluations_agree(self):
        x = mpf(2) ** -1000
        self.assertEqual((1 + x) - 1, 0)  # the starved evaluation is exactly 0 ...
        self.assertEqual(coeff.stable(lambda x: (1 + x) - 1, x), x)  # ... and is not accepted

    def test_r_definition(self):
        for n, w in ((1.0, 0.0), (0.0, 2.0), (0.3, 0.8), (1e-5, 1.0), (2.0, 1e-3)):
            out = coeff.r({"n": n, "w": w})
            n_, w_ = mpf(n), mpf(w)
            want = mp.pi / n_ if w == 0 else (2 / w_ if n == 0 else 2 * mp.atan2(n_, w_) / n_)
            self.assertLess(abs(out["value"] - want), REL * want, (n, w))

    def test_r_derivative_is_taken_at_fixed_w_with_respect_to_n_squared(self):
        for n, w in ((1.0, 0.0), (0.3, 0.8), (1e-5, 1.0), (2.0, 1e-3)):
            y, w_ = mpf(n) ** 2, mpf(w)
            angle = mp.atan2(mpf(n), w_)
            want = w_ / ((w_ * w_ + y) * y) - angle / mpf(n) ** 3
            got = coeff.r({"n": n, "w": w})["d_branch"]
            self.assertLess(abs(got - want), mpf(10) ** -90 * abs(want), (n, w))
        self.assertLess(abs(coeff.r({"n": 0.0, "w": 2.0})["d_branch"] + mpf(1) / 12), REL)

    def test_r_derivative_keeps_every_digit_at_tiny_n(self):
        """Regression: the stencil's y < 0 points cost ~47 digits for n in 1e-150..1e-80."""
        exact = gen_series.exact()["r"]
        for n in (1e-310, 1e-150, 1e-100, 1e-80, 1e-60):
            y = mpf(n) ** 2  # w = 1: d_branch = 2 S'(y), S = sum exact[j] y^j
            with mp.workdps(300):
                want = 2 * sum(
                    j * mpf(q.numerator) / q.denominator * y ** (j - 1)
                    for j, q in enumerate(exact[:8])
                    if j
                )
            got = coeff.r({"n": n, "w": 1.0})["d_branch"]
            self.assertLess(abs(got - want), REL * abs(want), n)

    def test_r_is_defined_for_w_nonnegative_and_not_both_zero(self):
        for bad in ({"n": 0.0, "w": 0.0}, {"n": 1.0, "w": -1.0}, {"n": -1.0, "w": 1.0}):
            with self.assertRaises(ValueError):
                coeff.r(bad)

    def test_r_inputs_are_the_unit_quaternion_of_each_angle(self):
        by_name = {s.name: s for s in R_STRATA}
        self.assertEqual(coeff.r_inputs(by_name["theta:exact0"]), [{"n": 0.0, "w": 1.0}])
        self.assertEqual(
            coeff.r_inputs(by_name["q:w0"]), [{"n": n, "w": 0.0} for n in (1.0, 1e-3, 1e3)]
        )
        (near_pi,) = coeff.r_inputs(by_name["theta:pi-1e-3"])
        self.assertLess(abs(near_pi["w"] - 5e-4), 1e-10)
        for rec in coeff.r_inputs(by_name["theta:subnormal"]):
            self.assertTrue(0 < rec["n"] < 2.2250738585072014e-308 and rec["w"] == 1.0)


class CosHalfTest(unittest.TestCase):
    def test_value_and_derivative_are_those_of_the_cosine(self):
        """d/dx cos(sqrt(x) / 2) = -sin(sqrt(x) / 2) / (4 sqrt(x)), -1/8 at 0, by calculus."""
        for theta in (0.0, 1e-310, 1e-12, 1e-4, 0.3, 1.0, 3.141592653589793, 6.0):
            out = coeff.evaluator("cos_half")({"theta": theta})
            with mp.workdps(400):
                t = mpf(theta)
                value = mp.cos(t / 2)
                deriv = -mpf(1) / 8 if theta == 0 else -mp.sin(t / 2) / (4 * t)
            self.assertLess(abs(out["value"] - value), REL * abs(value), theta)
            self.assertLess(abs(out["d_branch"] - deriv), REL * abs(deriv), theta)


class F32InputsTest(unittest.TestCase):
    """An `@f32` record's inputs are binary32 values, and its reference is at those inputs."""

    def test_the_reference_is_the_definition_at_the_rounded_input(self):
        theta32 = precision.to_f32(0.1)
        at32 = coeff.evaluator("b")({"theta": theta32})
        with mp.workdps(1200):
            want = (mpf(theta32) - mp.sin(mpf(theta32))) / mpf(theta32) ** 3
        self.assertLess(abs(at32["value"] - want), REL * abs(want))
        at64 = coeff.evaluator("b")({"theta": 0.1})
        self.assertGreater(abs(at32["value"] - at64["value"]), 1e-12)  # not the binary64 one

    def test_r_rounds_n_and_w_each_and_q_w0_keeps_w_zero(self):
        by_name = {s.name: s for s in COEFF_R_STRATA}
        for name in ("theta:1e-3", "theta:1e0", "theta:pi-1e-9", "theta:exact0", "q:w0"):
            wide, narrow = coeff.r_inputs(by_name[name]), coeff.r_inputs(by_name[name + "@f32"])
            self.assertEqual(
                narrow, [{k: precision.to_f32(v) for k, v in a.items()} for a in wide], name
            )
        second = coeff.r_inputs(by_name["q:w0@f32"])[1]
        self.assertEqual(second, {"n": precision.to_f32(1e-3), "w": 0.0})
        for rec in coeff.r_inputs(by_name["theta:subnormal@f32"]):  # from the binary32 decade
            self.assertTrue(0 < rec["n"] < 2.0**-126 and rec["w"] == 1.0, rec)


if __name__ == "__main__":
    unittest.main()
