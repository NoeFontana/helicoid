import math
import unittest

from mpmath import mp, mpf

from gen import precision
from gen.rng import stream, unit_quaternion_s3
from gen.strata import (
    COEFF_R_STRATA,
    COEFF_STRATA,
    QUAT_STRATA,
    R_STRATA,
    SCALAR_F32_STRATA,
    f32_twin,
)
from gen.strata import SCALAR_THETA_STRATA as STRATA

precision.setup()
BY_NAME = {s.name: s for s in STRATA}
NO_AXES = {"theta:exact0"}


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

    def test_r_sees_the_theta_strata_and_then_q_w0(self):
        self.assertEqual([s.name for s in R_STRATA], [*(s.name for s in STRATA), "q:w0"])
        self.assertEqual(sum(s.count for s in R_STRATA), 1713)

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

    def test_samples_pair_a_theta_with_an_axis_for_the_vector_ids(self):
        for s in STRATA:
            samples = s.samples()
            self.assertEqual(len(samples), {"theta:exact0": 1, "theta:dense": 801}.get(s.name, 64))
            if s.name == "theta:exact0":
                self.assertEqual(samples, [(0.0, (1.0, 0.0, 0.0))])
                continue
            self.assertEqual([a for _, a in samples], s.axes(len(samples)), s.name)
            self.assertEqual(len({a for _, a in samples}), len(samples), s.name)
            if s.count == 1:  # a fixed theta takes 64 axes
                self.assertEqual({t for t, _ in samples}, set(s.thetas()))
            else:  # a random or dense theta has its own
                self.assertEqual([t for t, _ in samples], s.thetas())

    def test_axes_are_uniform_on_s2(self):
        axes = BY_NAME["theta:dense"].samples()
        for k in range(3):
            xs = [a[k] for _, a in axes]
            self.assertLess(abs(sum(xs) / len(xs)), 0.09)  # sigma 0.020
            self.assertLess(abs(sum(x * x for x in xs) / len(xs) - 1 / 3), 0.053)  # sigma 0.0105
        self.assertTrue(all(abs(sum(c * c for c in a) - 1) < 4e-16 for _, a in axes))

    def test_theta_is_log_uniform_in_its_decade(self):
        positions = sorted(
            (mpf(t) / mpf(10) ** e).__float__()
            for e in range(-12, 0)
            for t in BY_NAME[f"theta:1e{e}"].thetas()
        )
        logs = [math.log10(x) for x in positions]  # position in the decade, uniform on [0, 1)
        n = len(logs)
        ks = max(max(x - i / n, (i + 1) / n - x) for i, x in enumerate(logs))
        self.assertLess(ks, 1.63 / math.sqrt(n))  # Kolmogorov-Smirnov, 1% level

    def test_the_quaternion_strata(self):
        self.assertEqual(
            [s.name for s in QUAT_STRATA], [*(s.name for s in STRATA), "q:w0", "q:nonunit"]
        )
        w0, nonunit = QUAT_STRATA[-2:]
        self.assertEqual((len(w0.quaternions()), len(nonunit.quaternions())), (64, 64))
        self.assertEqual(w0.quaternions(), w0.quaternions())
        self.assertNotEqual(w0.quaternions(), nonunit.quaternions())
        self.assertEqual(len(set(w0.points(64) + w0.points(64))), 64)

    def test_haar_quaternions_are_uniform_on_s3(self):
        rng = stream(1, "s3")
        qs = [unit_quaternion_s3(rng) for _ in range(2000)]
        self.assertTrue(all(abs(sum(c * c for c in q) - 1) < mpf(10) ** -100 for q in qs))
        for k in range(4):
            xs = [float(q[k]) for q in qs]
            self.assertLess(abs(sum(xs) / len(xs)), 0.03)  # sigma 0.011
            self.assertLess(abs(sum(x * x for x in xs) / len(xs) - 0.25), 0.02)  # sigma 0.0056

    def test_draws_do_not_depend_on_the_ambient_precision(self):
        for name in ("theta:1e-8", "theta:1e0", "theta:pi-1e-5", "theta:dense"):
            s = BY_NAME[name]
            expected = (s.thetas(), s.axes(4))
            for dps in (15, 53, 400):
                with mp.workdps(dps):
                    self.assertEqual((s.thetas(), s.axes(4)), expected, (name, dps))
            self.assertEqual(mp.dps, precision.DPS)


F32_BY_NAME = {s.name: s for s in SCALAR_F32_STRATA}


class F32StrataTest(unittest.TestCase):
    """0016 item 1: beside every stratum S of a scalar id, `S@f32`, appended after all of them."""

    def test_catalogue_appends_every_twin_in_the_order_of_its_source(self):
        twins = [f"{s.name}@f32" for s in STRATA]
        self.assertEqual([s.name for s in COEFF_STRATA], [*(s.name for s in STRATA), *twins])
        self.assertEqual(
            [s.name for s in COEFF_R_STRATA], [*(s.name for s in R_STRATA), *twins, "q:w0@f32"]
        )
        self.assertEqual(sum(s.count for s in COEFF_STRATA), 2 * 1710)
        self.assertEqual(sum(s.count for s in COEFF_R_STRATA), 2 * 1713)
        self.assertFalse(any(s.f32 for s in (*STRATA, *R_STRATA, *QUAT_STRATA)))

    def test_a_twin_is_its_source_rounded_to_nearest_even_binary32(self):
        for s in SCALAR_F32_STRATA:
            thetas = s.thetas()
            self.assertEqual(len(thetas), s.count, s.name)
            self.assertTrue(all(precision.is_binary32(t) for t in thetas), s.name)
            if s.name != "theta:subnormal@f32":
                base = BY_NAME[s.name[: -len("@f32")]]
                self.assertEqual(thetas, [precision.to_f32(t) for t in base.thetas()], s.name)
        self.assertEqual(F32_BY_NAME["theta:exact0@f32"].thetas(), [0.0])
        xs = F32_BY_NAME["theta:dense@f32"].thetas()
        self.assertTrue(len(xs) == 801 and all(a < b for a, b in zip(xs, xs[1:], strict=False)))

    def test_subnormal_is_a_binary32_subnormal_decade_not_the_rounding_of_the_binary64_one(self):
        literal = {precision.to_f32(t) for t in BY_NAME["theta:subnormal"].thetas()}
        self.assertEqual(literal, {0.0})  # rounding [1e-310, 1e-309) is 64 zeros
        xs = F32_BY_NAME["theta:subnormal@f32"].thetas()
        self.assertTrue(all(1e-40 * (1 - 1e-4) < x < 1e-39 * (1 + 1e-4) for x in xs))
        self.assertTrue(all(x < 2.0**-126 for x in xs))  # below the smallest normal binary32
        self.assertGreater(len(set(xs)), 60)
        self.assertTrue(all(precision.to_f32(x * x) == 0.0 for x in xs))  # theta^2 underflows
        self.assertTrue(all(1 / x > 3.4028234663852886e38 for x in xs))  # 1 / theta overflows

    def test_several_pi_strata_round_to_one_binary32_and_every_record_is_kept(self):
        values = {k: F32_BY_NAME[f"theta:pi-1e-{k}@f32"].thetas() for k in range(1, 13)}
        self.assertEqual(len({v[0] for v in values.values()}), 8)  # k = 1..7, and 8..12 as one
        self.assertEqual({v[0] for k, v in values.items() if k >= 8}, {precision.to_f32(math.pi)})
        self.assertGreater(values[8][0], math.pi)  # the binary32 nearest pi is above pi

    def test_a_twin_is_named_by_the_suffix_and_draws_at_any_ambient_precision(self):
        twin = f32_twin(BY_NAME["theta:1e-8"])
        self.assertEqual((twin.name, twin.source), ("theta:1e-8@f32", BY_NAME["theta:1e-8"]))
        for s in (twin, F32_BY_NAME["theta:subnormal@f32"]):
            with mp.workdps(400):
                self.assertEqual(s.thetas(), F32_BY_NAME[s.name].thetas(), s.name)


if __name__ == "__main__":
    unittest.main()
