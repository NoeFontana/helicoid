"""The `real_*` ids: their strata's ranges and signs, and the cross-checks' reach."""

import unittest

from mpmath import mp, mpf

from gen import check_real, precision, real
from gen.precision import is_binary32
from gen.registry import FUNCTIONS

precision.setup()
REAL = [name for name in FUNCTIONS if name.startswith("real_")]


def by_name(strata) -> dict:
    return {s.name: s for s in strata}


class StrataTest(unittest.TestCase):
    def test_decades_and_their_f32_layout(self):
        names = list(by_name(real.SQRT_STRATA))
        self.assertEqual(names, list(by_name(real.CBRT_STRATA)))
        self.assertEqual(
            names,
            [
                *(f"x:1e{e}" for e in (-300, -100, -10, -1, 0, 1, 10, 100, 300)),
                "x:subnormal",
                *(f"x:1e{e}@f32" for e in (-37, -10, -1, 0, 1, 10, 37)),
                "x:subnormal@f32",
            ],
        )
        strata = by_name(real.SQRT_STRATA)
        for e in (-300, 0, 300):
            xs = [r["x"] for r in strata[f"x:1e{e}"].records()]
            self.assertTrue(all(mpf(10) ** e <= x < mpf(10) ** (e + 1) for x in xs), e)
        # a shared range: the binary64 draw rounded once; an own range: its own stream
        x64 = [r["x"] for r in strata["x:1e0"].records()]
        x32 = [r["x"] for r in strata["x:1e0@f32"].records()]
        self.assertTrue(all(abs(a / b - 1) < 2**-24 for a, b in zip(x64, x32, strict=True)))
        sub = [r["x"] for r in strata["x:subnormal@f32"].records()]
        self.assertTrue(all(1e-40 * (1 - 2**-10) <= x < 1e-39 * (1 + 2**-10) for x in sub))
        self.assertTrue(all(0 < r["x"] < 2.2250738585072014e-308 for r in strata["x:subnormal"].records()))

    def test_signs_alternate(self):
        cbrt = [r["x"] for r in by_name(real.CBRT_STRATA)["x:1e-1"].records()]
        self.assertEqual([x < 0 for x in cbrt[:4]], [False, True, False, True])
        nd = by_name(real.DIV_STRATA)["nd:wide@f32"].records()
        self.assertEqual(
            [(r["n"] < 0, r["d"] < 0) for r in nd[:4]],
            [(False, False), (True, False), (False, True), (True, True)],
        )

    def test_acos_never_reaches_plus_or_minus_one(self):
        strata = by_name(real.ACOS_STRATA)
        for name, j_max in (("x:near+1", 52), ("x:near-1@f32", 23)):
            for r in strata[name].records():
                j = -mp.log(1 - abs(mpf(r["x"])), 2)
                self.assertTrue(j == int(j) and 2 <= j <= j_max, (name, r))

    def test_atan2_ratios_and_quadrants(self):
        strata = by_name(real.ATAN2_STRATA)
        for name, exp in (("yx:ratio-1e-300", -300), ("yx:ratio-1e+30@f32", 30)):
            records = strata[name].records()
            for r in records:
                ratio = mp.log10(abs(mpf(r["y"]) / mpf(r["x"])))
                self.assertLess(abs(ratio - exp), 1e-6, name)
            signs = [(r["x"] > 0, r["y"] > 0) for r in records[:4]]
            self.assertEqual(signs, [(True, True), (False, True), (False, False), (True, False)])

    def test_f32_inputs_are_binary32_and_every_output_is_a_normal_number(self):
        for name in REAL:
            spec = FUNCTIONS[name]
            for s in spec.strata:
                top, low = (mpf(2) ** 128, 2.0**-126) if s.f32 else (mpf(2) ** 1024, 2.0**-1022)
                for r in s.records()[:8]:
                    if s.f32:
                        self.assertTrue(all(is_binary32(x) for x in r.values()), (name, s.name))
                    for v in spec.evaluate(r).values():
                        self.assertTrue(low <= abs(v) < top, (name, s.name, r, v))


EXAMPLES = {
    "real_sqrt": {"x": 0.7},
    "real_cbrt": {"x": -3e-5},
    "real_sin_cos": {"x": 0.7},
    "real_acos": {"x": 1 - 2.0**-40},
    "real_atan2": {"y": -1e-30, "x": -2.5},
    "real_div": {"n": -3.0, "d": 7e-3},
}


class CrossCheckTest(unittest.TestCase):
    def test_a_perturbation_beyond_100_digits_fails_in_every_output(self):
        for name, inputs in EXAMPLES.items():
            spec = FUNCTIONS[name]
            good = spec.evaluate(inputs)
            spec.check(inputs, good)
            for key in good:
                with self.assertRaises(check_real.CrossCheckError, msg=(name, key)):
                    spec.check(inputs, {**good, key: good[key] * (1 + mpf(10) ** -95)})
                spec.check(inputs, {**good, key: good[key] * (1 + mpf(10) ** -105)})

    def test_mp_diff_holds_the_derivative_to_60_digits(self):
        d = real.sin_cos({"x": 1e22})["d_sin"]
        check_real.against_diff("cos", mp.sin, 1e22, d)
        check_real.against_diff("cos", mp.sin, 1e22, d * (1 + mpf(10) ** -65))
        with self.assertRaises(check_real.CrossCheckError):
            check_real.against_diff("cos", mp.sin, 1e22, d * (1 + mpf(10) ** -55))

    def test_0056_s_step_would_not_reach_60_digits(self):
        """h = 1e-40 |x| leaves h^2 sin'''/6 ~ 1e-37 of cos at x = 1e22, and 1e-49 of acos' at
        1 - 2^-52: why `check_real` steps at 1e-70 |x|."""
        for f, x, want in (
            (mp.sin, mpf(1e22), mp.cos(mpf(1e22))),
            (mp.acos, 1 - mpf(2) ** -52, -1 / mp.sqrt(1 - (1 - mpf(2) ** -52) ** 2)),
        ):
            with mp.workdps(300):
                got = mp.diff(f, x, h=mpf(10) ** -40 * abs(x))
            self.assertGreater(abs(got / want - 1), mpf(10) ** -60)


if __name__ == "__main__":
    unittest.main()
