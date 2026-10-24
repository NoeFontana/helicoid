"""The `real_*` ids: their strata's ranges and signs, and the cross-checks' reach."""

import math
import unittest

from mpmath import mp, mpf

from gen import check_real, precision, real
from gen.precision import is_binary32, to_f32, to_f64
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
            if spec.converts:
                continue  # binary64 inputs and edge outputs by design (0087), held below
            for s in spec.strata:
                top, low = (mpf(2) ** 128, 2.0**-126) if s.f32 else (mpf(2) ** 1024, 2.0**-1022)
                if name == "real_exp" and s.name.startswith("x:subnormal"):
                    continue  # its outputs are subnormal by design (0068), held below
                # `yx:axis` (0105, draft) has exact zeros by design, which the metric does not
                # divide by; held below
                zero_ok = s.name.startswith("yx:axis")
                for r in s.records()[:8]:
                    if s.f32:
                        self.assertTrue(all(is_binary32(x) for x in r.values()), (name, s.name))
                    for v in spec.evaluate(r).values():
                        if zero_ok and v == 0:
                            continue
                        self.assertTrue(low <= abs(v) < top, (name, s.name, r, v))

    def test_the_0105_atan2_strata_hold_what_they_are_named_for(self):
        strata = by_name(real.ATAN2_STRATA)
        added = ["centre-mid", "fold", "diagonal", "axis", "window-lo", "domain-low"]
        added += ["domain-high", "subnormal"]
        names = list(strata)
        # after the earlier strata of their precision
        self.assertEqual(names[7:15], [f"yx:{n}" for n in added])
        self.assertEqual(
            names[22:], [f"yx:{n}@f32" for n in added if n not in ("centre-mid", "window-lo")]
        )
        for name, s in strata.items():
            f32 = s.f32
            lo, hi = (mpf(2) ** -63, mpf(2) ** 63) if f32 else (mpf(2) ** -511, mpf(2) ** 511)
            records = s.records()
            self.assertEqual(len(records), 64, name)
            for r in records:
                y, x = mpf(r["y"]), mpf(r["x"])
                # the derivative's domain (docs/NUMERICS.md section 12): the larger argument
                self.assertTrue(lo <= max(abs(x), abs(y)) < hi, (name, r))
                self.assertFalse(r["y"] == 0 and math.copysign(1, r["y"]) < 0, (name, r))
                if f32:
                    self.assertTrue(is_binary32(r["y"]) and is_binary32(r["x"]), (name, r))
            nd = [(min(abs(r["y"]), abs(r["x"])), max(abs(r["y"]), abs(r["x"]))) for r in records]
            ulp = 2.0**-23 if f32 else 2.0**-52
            base = name.split("@")[0]
            if base == "yx:centre-mid":
                mids = [d * (2 * ((i // 8) % 8) + 1) / 16 for i, (_, d) in enumerate(nd)]
                for (n, d), mid in zip(nd, mids, strict=True):
                    self.assertLessEqual(abs(n - mid), 3 * ulp * mid, (n, d))
                self.assertTrue(any(n == mid for (n, _), mid in zip(nd, mids, strict=True)))
            if base == "yx:fold":
                self.assertTrue(all(abs(n - d / 2) <= 4 * ulp * d / 2 for n, d in nd), name)
                self.assertTrue({n < d / 2 for n, d in nd} == {True, False}, name)
            if base == "yx:diagonal":
                self.assertTrue(all(d - 4 * ulp * d <= n <= d for n, d in nd), name)
                self.assertTrue(any(n == d for n, d in nd), name)
            if base[3:] not in added:
                continue
            if base != "yx:axis":
                for i, r in enumerate(records):
                    y, x = abs(r["y"]), abs(r["x"])
                    self.assertTrue(y >= x if (i // 4) % 2 else y <= x, (name, i, r))
                signs = [(r["x"] > 0, r["y"] > 0) for r in records[:4]]
                self.assertEqual(signs, [(True, True), (False, True), (False, False), (True, False)])
            else:
                first = [(r["y"], math.copysign(1, r["x"])) for r in records[:4]]
                self.assertEqual([y for y, _ in first[:2]], [0.0, 0.0])
                self.assertEqual([s for _, s in first], [1.0, -1.0, 1.0, -1.0])
                self.assertEqual([r["x"] for r in records[2:4]], [0.0, 0.0])
            if base == "yx:window-lo":
                self.assertTrue(all(abs(d / 2.0**-460 - 1) <= 4 * ulp for n, d in nd), name)
                self.assertTrue({d < 2.0**-460 for n, d in nd} == {True, False}, name)
            if base == "yx:subnormal":
                tiny = 2.0**-126 if f32 else 2.0**-1022
                self.assertTrue(all(0 < n < tiny for n, d in nd), name)


    def test_the_exp_strata_hold_what_they_are_named_for(self):
        """0068: `x:subnormal` gives a subnormal e^x, `x:saturated` an expm1 of -1 at the precision
        with a normal derivative, and `x:large` stays below ln MAX."""
        for f32 in (False, True):
            tiny, top = (2.0**-126, mpf(2) ** 128) if f32 else (2.0**-1022, mpf(2) ** 1024)
            suffix = "@f32" if f32 else ""
            exp = by_name(real.EXP_STRATA)
            for r in exp["x:subnormal" + suffix].records():
                v = real.exp(r)["value"]
                self.assertTrue(0 < v < tiny, r)
            for r in exp["x:large" + suffix].records():
                self.assertLess(real.exp(r)["value"], top)
            for r in by_name(real.EXPM1_STRATA)["x:saturated" + suffix].records():
                out = real.expm1(r)
                rounded = to_f32(out["value"]) if f32 else to_f64(out["value"])
                self.assertEqual(rounded, -1.0, r)
                self.assertTrue(tiny <= out["d"], r)

    def test_the_round_strata_hold_what_they_are_named_for(self):
        """0087: at binary32 the inputs are binary64 that are not binary32, the reference is
        `to_f32` of them and finite, and each stratum reaches its edge."""
        strata = by_name(real.ROUND_STRATA)
        self.assertEqual(
            list(strata),
            [
                "x:generic",
                "x:subnormal",
                *(
                    f"x:{n}@f32"
                    for n in ("generic", "tie", "near-tie", "overflow", "underflow")
                ),
            ],
        )
        for r in strata["x:subnormal"].records():
            self.assertTrue(0 < abs(r["x"]) < 2.0**-1022, r)
        edge = mpf(2) ** 128 - mpf(2) ** 103
        values = {}
        for name, s in strata.items():
            records = s.records()
            self.assertEqual([r["x"] < 0 for r in records[:2]], [False, True], name)
            values[name] = [real.round_from_f64(r, s.f32)["value"] for r in records]
            if not s.f32:
                self.assertEqual(values[name], [mpf(r["x"]) for r in records])
                continue
            for r, v in zip(records, values[name], strict=True):
                self.assertEqual(v, mpf(to_f32(r["x"])))
                self.assertLess(abs(mpf(r["x"])), edge, (name, r))
            self.assertFalse(any(is_binary32(r["x"]) for r in records), name)
        tie = [mpf(r["x"]) for r in strata["x:tie@f32"].records()]
        for x, v in zip(tie, values["x:tie@f32"], strict=True):
            # halfway: the result is one of the two neighbours, at half a spacing either way
            self.assertEqual(abs(v - x), abs(mpf(to_f32(2 * x - v)) - x), x)
        over = [abs(v) for v in values["x:overflow@f32"]]
        self.assertEqual(max(over), (2**24 - 1) * mpf(2) ** 104)
        self.assertEqual(min(over), (2**24 - 2) * mpf(2) ** 104)
        under = [abs(v) for v in values["x:underflow@f32"]]
        self.assertEqual(under[:2], [0, mpf(2) ** -148])
        self.assertTrue(any(0 < v < 2.0**-126 for v in under))

    def test_round_cross_check_refuses_another_rounding(self):
        x = 1 + 2.0**-24 + 2.0**-40  # nearest binary32 is 1 + 2^-23
        check_real.real_round_from_f64({"x": x}, real.round_from_f64({"x": x}, True))
        check_real.real_round_from_f64({"x": x}, real.round_from_f64({"x": x}, False))
        for bad in ({"value": mpf(1), "d": mpf(0)}, {"value": mpf(x), "d": mpf(1)}):
            with self.assertRaises(check_real.CrossCheckError):
                check_real.real_round_from_f64({"x": x}, bad)


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

    def test_the_atan2_cross_check_holds_on_the_axes(self):
        """0105 (draft): a zero argument and a zero derivative, at each of `yx:axis`'s four."""
        spec = FUNCTIONS["real_atan2"]
        for inputs in ({"y": 0.0, "x": 3.0}, {"y": 0.0, "x": -3.0}, {"y": -3.0, "x": -0.0}):
            good = spec.evaluate(inputs)
            spec.check(inputs, good)
            for key in good:
                bad = good[key] * (1 + mpf(10) ** -95) or mpf(10) ** -95
                with self.assertRaises(check_real.CrossCheckError, msg=(inputs, key)):
                    spec.check(inputs, {**good, key: bad})

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
