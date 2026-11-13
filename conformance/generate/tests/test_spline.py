"""The spline ids (0116) against what shares nothing with `gen.spline`: the basis's published rows
and check values, linear precision (a window on one geodesic interpolates it at constant speed),
the identity window's Jacobians, the strata's own properties, planted errors for the
cross-checks, and the committed files."""

import json
import unittest
from fractions import Fraction

from mpmath import mp, mpf

from gen import check_spline, manifest, precision, so3, spline
from gen.check import CrossCheckError
from gen.precision import is_binary32
from gen.registry import FUNCTIONS

precision.setup()
COMMITTED = manifest.ROOT.parent / "corpus"
TOL = mpf(10) ** -100


def near(got, want):
    for g, w in zip(got, want, strict=True):
        if abs(mpf(g) - mpf(w)) > TOL:
            raise AssertionError(f"{mp.nstr(mpf(g), 40)} vs {mp.nstr(mpf(w), 40)}")


def about_z(angle) -> list[float]:
    """The rotation by `angle` about z, rounded to binary64."""
    return [float(mp.cos(mpf(angle) / 2)), 0.0, 0.0, float(mp.sin(mpf(angle) / 2))]


class BasisTest(unittest.TestCase):
    def test_cox_de_boor_is_section_21s_integers_at_every_order(self):
        for k in spline.RN_ORDERS:
            n = check_spline.numerators(k)
            f = Fraction(1, 1)
            for i in range(2, k):
                f *= i
            for s in range(k):
                self.assertEqual(list(spline.basis(k)[s]), [Fraction(c) / f for c in n[s]], k)

    def test_the_cumulative_rows_at_k_4_are_qins(self):
        rows = [[6, 0, 0, 0], [5, 3, -3, 1], [1, 3, 3, -2], [0, 0, 0, 1]]
        got = [[c * 6 for c in row] for row in spline.cumulative(4)]
        self.assertEqual(got, rows)

    def test_the_basis_is_a_partition_of_unity(self):
        for k in spline.RN_ORDERS:
            total = [sum(b[c] for b in spline.basis(k)) for c in range(k)]
            self.assertEqual(total, [1] + [0] * (k - 1), k)

    def test_the_gram_check_values(self):
        self.assertEqual(spline.gram_exact(4, 1, Fraction(1))[0][0], Fraction(1, 20))
        self.assertEqual(spline.gram_exact(6, 2, Fraction(1))[0][0], Fraction(1, 252))
        # h^{1 - 2r}: r = 2 at h = 2 scales by 2^-3.
        self.assertEqual(spline.gram_exact(6, 2, Fraction(1, 2))[0][0], Fraction(1, 252 * 8))


class HandComputedTest(unittest.TestCase):
    """Linear precision, `sum_s s b_s(u) = u + (K - 2)/2`: a window `X_s = Exp(s delta e_z)` is
    the rotation by `(u + (K - 2)/2) delta` about z, at body rate `delta e_z`, unaccelerated."""

    def test_a_window_on_one_axis_turns_at_constant_rate(self):
        delta, u = mpf(3) / 10, 0.375
        for k in spline.GROUP_ORDERS:
            q = [c for s in range(k) for c in about_z(s * delta)]
            inputs = {"q": q, "u": u, "inv_h": 1.0}
            out = [spline.so3_output(r)(inputs) for r in range(3)]
            # The rounded controls are evenly spaced to 1e-16: the rate is delta to that.
            for g, w in zip(out[1]["omega"] + out[2]["alpha"], [0, 0, delta, 0, 0, 0]):
                self.assertLess(abs(g - w), mpf(10) ** -15)
            # Rotations about one axis commute: the angle is sum lambda_j theta_j, exactly, with
            # the increments the rounded controls hold.
            stored = spline.quats(q)
            theta = [spline.so3_increment(a, b)[2] for a, b in zip(stored, stored[1:])]
            lam = [
                [spline.at(spline.deriv_n(list(c), r), Fraction(u)) for c in spline.cumulative(k)]
                for r in range(3)
            ]
            angle = mp.atan2(stored[0][3], stored[0][0]) * 2
            angle += sum(spline._mp(lam[0][j]) * theta[j - 1] for j in range(1, k))
            near(out[0]["q"], [mp.cos(angle / 2), 0, 0, mp.sin(angle / 2)])
            for r, key in ((1, "omega"), (2, "alpha")):
                rate = sum(spline._mp(lam[r][j]) * theta[j - 1] for j in range(1, k))
                near(out[r][key], [0, 0, rate])

    def test_points_on_a_line_move_at_constant_velocity(self):
        v = [0.5, -0.25, 2.0]
        for k in spline.RN_ORDERS:
            inputs = {"p": [s * c for s in range(k) for c in v], "u": 0.375, "inv_h": 1.0}
            at = Fraction(3, 8) + Fraction(k - 2, 2)
            near(spline.rn_output(0)(inputs)["p"], [at * Fraction(c) for c in v])
            near(spline.rn_output(1)(inputs)["v"], v)
            near(spline.rn_output(2)(inputs)["a"], [0, 0, 0])

    def test_at_the_identity_each_control_moves_the_value_by_its_basis_function(self):
        k, u = 4, 0.25
        inputs = {"q": [1.0, 0.0, 0.0, 0.0] * k, "u": u, "inv_h": 1.0}
        jac = spline.so3_jacobian(0)(inputs)["jac"]
        self.assertEqual(jac.shape, (3, 3 * k))
        for m, b in enumerate(spline.basis(k)):
            weight = spline.at(b, Fraction(u))
            block = jac.data[9 * m : 9 * m + 9]
            near(block, [weight if i % 4 == 0 else 0 for i in range(9)])

    def test_a_control_the_output_does_not_reach_reads_an_exact_zero(self):
        """At `u = 0` the last control's weight vanishes to order K - 2: its block is 0, not
        noise."""
        record = FUNCTIONS["so3_spline_vel_jac"].strata[5].records()[0]  # u:0
        self.assertEqual(record["u"], 0.0)
        k = len(record["q"]) // 4
        jac = spline.so3_jacobian(1)(record)["jac"].data
        self.assertTrue(all(x == 0 for x in jac[9 * (k - 1) :]))
        self.assertTrue(any(x != 0 for x in jac[: 9 * (k - 1)]))


def increments(record: dict, group: str) -> list:
    if group == "rn":
        p = [[mpf(c) for c in record["p"][3 * i : 3 * i + 3]] for i in range(len(record["p"]) // 3)]
        return [so3.norm([b - a for a, b in zip(x, y)]) for x, y in zip(p, p[1:])]
    q = spline.quats(record["q"])
    return [
        so3.norm(check_spline.log_atan2(so3.qmul(so3.qconj(a), b))) for a, b in zip(q, q[1:])
    ]


class StrataTest(unittest.TestCase):
    """Each stratum holds the increments, `u` and `inv_h` it names, at both precisions: the
    controls are formed from the increments and rounded, never drawn."""

    def test_the_strata_hold_what_they_name(self):
        groups = (("so3", spline.SO3_STRATA), ("se3", spline.SE3_STRATA), ("rn", spline.RN_STRATA))
        for group, strata in groups:
            orders = spline.RN_ORDERS if group == "rn" else spline.GROUP_ORDERS
            for s in strata:
                name = s.name.removesuffix("@f32")
                records = s.records()
                self.assertEqual(len(records), 6, s.name)
                for i, r in enumerate(records):
                    width = 3 if group == "rn" else 4
                    key = "p" if group == "rn" else "q"
                    self.assertEqual(len(r[key]) // width, orders[i % len(orders)], s.name)
                    if s.f32:
                        values = [x for v in r.values() for x in (v if type(v) is list else [v])]
                        self.assertTrue(all(is_binary32(x) for x in values), s.name)
                    d = increments(r, group)
                    if name == "inc:tiny":
                        self.assertTrue(all(x <= 1.01e-8 for x in d), (s.name, d))
                    elif name == "inc:near-pi":
                        self.assertTrue(all(0 < mp.pi - x < 1.05e-6 for x in d), (s.name, d))
                    elif name == "inc:large":
                        top = mp.pi - 1e-3 + 1e-6
                        self.assertTrue(all(2.5 - 1e-6 < x < top for x in d), s.name)
                    elif name == "inc:mixed" and len(d) == 1:  # K = 2: tiny, then large
                        self.assertTrue(d[0] <= 1.01e-8 if i % 2 == 0 else d[0] > 2.5 - 1e-6)
                    elif name == "inc:mixed":
                        self.assertEqual(sum(x <= 1.01e-8 for x in d), 1, (s.name, d))
                        self.assertEqual(sum(x >= 2.5 - 1e-6 for x in d), 1, (s.name, d))
                    else:
                        self.assertTrue(all(1e-3 - 1e-6 < x < 1 + 1e-6 for x in d), (s.name, d))
                    u = {"u:0": 0.0, "u:half": 0.5, "u:1-ulp": 1 - 2.0 ** (-24 if s.f32 else -53)}
                    if name in u:
                        self.assertEqual(r["u"], u[name], s.name)
                    else:
                        self.assertTrue(0 <= r["u"] <= 1, s.name)  # binary32 may round to 1
                    want_h = {"h:1e-3": 1000.0, "h:1e3": 1e-3}.get(name, 1.0)
                    self.assertAlmostEqual(r["inv_h"], want_h, delta=want_h * 1e-7)

    def test_the_f32_twin_draws_the_same_increments(self):
        for s, t in zip(spline.SO3_STRATA[:10], spline.SO3_STRATA[10:], strict=True):
            self.assertEqual(t.name, s.name + "@f32")
            for a, b in zip(s.records(), t.records(), strict=True):
                for x, y in zip(increments(a, "so3"), increments(b, "so3"), strict=True):
                    # A tiny increment survives binary32 relatively (its anchor is exact); the
                    # others move by binary32's rounding of O(1) components.
                    self.assertLess(abs(x - y), 1e-6 * x if x < 1e-6 else 1e-6, s.name)

    def test_the_gram_strata_hold_every_order_and_derivative_at_their_h(self):
        for s in spline.GRAM_STRATA:
            records = s.records()
            want = [(float(k), float(r)) for k, r in spline.GRAM_CASES]
            self.assertEqual([(r["k"], r["r"]) for r in records], want)
            h = {"h:2^-10": 2.0**10, "h:1": 1.0, "h:2^10": 2.0**-10}[s.name.removesuffix("@f32")]
            self.assertTrue(all(r["inv_h"] == h for r in records), s.name)


class CrossCheckTest(unittest.TestCase):
    """A planted error in each output fails its cross-check."""

    @classmethod
    def setUpClass(cls):
        cls.so3 = spline.SO3_STRATA[0].records()[0]
        cls.se3 = spline.SE3_STRATA[0].records()[0]
        cls.rn = spline.RN_STRATA[0].records()[2]

    def planted(self, name: str, inputs: dict, key: str, k: int = 0):
        spec = FUNCTIONS[name]
        out = spec.evaluate(inputs)
        spec.check(inputs, out)
        bad = dict(out)
        value = out[key]
        if hasattr(value, "data"):
            data = list(value.data)
            data[k] += mpf(10) ** -90
            bad[key] = type(value)(value.shape, data)
        else:
            bad[key] = [c + (mpf(10) ** -90 if i == k else 0) for i, c in enumerate(value)]
        with self.assertRaises(CrossCheckError, msg=name):
            spec.check(inputs, bad)

    def test_so3(self):
        for name, key in (("value", "q"), ("vel", "omega"), ("acc", "alpha")):
            self.planted(f"so3_spline_{name}", self.so3, key, 1)

    def test_so3_jacobians(self):
        for order, name in enumerate(("value", "vel", "acc")):
            self.planted(f"so3_spline_{name}_jac", self.so3, "jac", 7 + order)

    def test_se3(self):
        self.planted("se3_spline_value", self.se3, "x", 2)
        for name in ("vel", "acc"):
            self.planted(f"se3_spline_{name}", self.se3, "twist", 4)

    def test_rn(self):
        for name, key in (("value", "p"), ("vel", "v"), ("acc", "a")):
            self.planted(f"rn_spline_{name}", self.rn, key, 0)
            self.planted(f"rn_spline_{name}_jac", self.rn, "k", 1)

    def test_gram(self):
        self.planted("rn_spline_gram", {"k": 5.0, "r": 2.0, "inv_h": 1024.0}, "g", 3)


class CommittedTest(unittest.TestCase):
    def test_every_spline_file_holds_its_strata_in_order(self):
        for name, spec in FUNCTIONS.items():
            if "spline" not in name:
                continue
            lines = (COMMITTED / f"{name}.jsonl").read_text().splitlines()
            strata = [json.loads(line)["stratum"] for line in lines]
            want = [s.name for s in spec.strata for _ in range(s.count)]
            self.assertEqual(strata, want, name)


if __name__ == "__main__":
    unittest.main()
