"""The SE(3) chart ids (0060 decision 8) against hand-computed poses, the relations between the
three charts that `gen.charts` never uses (docs/maths/charts.md CH.4(b)), planted errors for the
cross-checks, and the input draw."""

import unittest

from mpmath import mp, mpf

from gen import charts, check_charts, precision, so3
from gen.check import CrossCheckError
from gen.registry import FUNCTIONS
from gen.strata import SEN3_SAMPLES, SEN3_STRATA

precision.setup()
TOL = mpf(10) ** -100
IDENTITY = [1.0, 0.0, 0.0, 0.0]


def quarter_z():
    return [float(mp.cos(mp.pi / 4)), 0.0, 0.0, float(mp.sin(mp.pi / 4))]


def near(got, want):
    for g, w in zip(got, want, strict=True):
        if abs(mpf(g) - mpf(w)) > TOL * max(1, abs(mpf(w))):
            raise AssertionError(f"{mp.nstr(mpf(g), 40)} vs {mp.nstr(mpf(w), 40)}")


class HandComputedTest(unittest.TestCase):
    def test_the_translation_moves_in_the_frame_each_chart_names(self):
        rec = {"q0": quarter_z(), "x0": [1.0, 2.0, 3.0], "tau": [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]}
        n = so3.normalize([mpf(c) for c in quarter_z()])
        R = so3.rot(n)
        body = [mpf(c) + R[r][0] for r, c in enumerate(rec["x0"])]
        near(charts.retract("decoupled")(rec)["x"], body)
        near(charts.retract("screw")(rec)["x"], body)  # phi = 0: J_l = I
        near(charts.retract("world")(rec)["x"], [2, 2, 3])

    def test_screw_and_decoupled_differ_by_the_jl_term(self):
        # CH.4(b): the translations differ by R (J_l(phi) - I) rho, which is nonzero off the axis.
        rec = {"q0": IDENTITY, "x0": [0.0, 0.0, 0.0], "tau": [0.0, 0.0, 1.0, 1.0, 0.0, 0.0]}
        s = charts.retract("screw")(rec)["x"]
        d = charts.retract("decoupled")(rec)["x"]
        jl = so3.jl([mpf(0), mpf(0), mpf(1)])
        near(s, [jl[r][0] for r in range(3)])
        near(d, [1, 0, 0])


class TransitionTest(unittest.TestCase):
    def test_local_of_another_chart_is_the_ch4b_transition(self):
        # Phi^{Dec -> WT}(phi, rho) = (phi, R rho), exactly: the world local of a decoupled retract.
        rec = {"q0": quarter_z(), "x0": [0.5, -1.0, 2.0], "tau": [0.1, -0.2, 0.3, 0.4, 0.5, -0.6]}
        y = charts.retract("decoupled")(rec)
        got = charts.local("world")({"q0": rec["q0"], "x0": rec["x0"], "q1": y["q"], "x1": y["x"]})
        R = so3.rot(so3.normalize([mpf(c) for c in rec["q0"]]))
        rho = [mpf(c) for c in rec["tau"][3:]]
        want = rec["tau"][:3] + [sum(R[r][c] * rho[c] for c in range(3)) for r in range(3)]
        near(got["phi"] + got["rho"], want)


class CrossCheckTest(unittest.TestCase):
    def test_a_wrong_answer_fails_each_check(self):
        rec = {"q0": quarter_z(), "x0": [0.5, -1.0, 2.0], "tau": [0.1, -0.2, 0.3, 0.4, 0.5, -0.6]}
        for chart in charts.CHARTS:
            out = charts.retract(chart)(rec)
            check_charts.retract(chart)(rec, out)
            bad = {"q": out["q"], "x": [out["x"][0] + mpf(10) ** -50, *out["x"][1:]]}
            with self.assertRaises(CrossCheckError):
                check_charts.retract(chart)(rec, bad)
            loc = {"q0": rec["q0"], "x0": rec["x0"], "q1": out["q"], "x1": out["x"]}
            ans = charts.local(chart)(loc)
            check_charts.local(chart)(loc, ans)
            with self.assertRaises(CrossCheckError):
                bad = {"phi": ans["phi"], "rho": [*ans["rho"][:2], ans["rho"][2] + 1e-50]}
                check_charts.local(chart)(loc, bad)


class InputTest(unittest.TestCase):
    def test_every_id_sees_every_se3_stratum_with_binary64_inputs(self):
        for chart in charts.CHARTS:
            for op in ("retract", "local"):
                spec = FUNCTIONS[f"se3_{chart}_{op}"]
                self.assertEqual([s.name for s in spec.strata], [s.name for s in SEN3_STRATA])
        stratum = SEN3_STRATA[0]
        recs = charts.local_inputs("decoupled")(stratum)
        self.assertEqual(len(recs), SEN3_SAMPLES)
        for rec in recs:
            for key in ("q0", "x0", "q1", "x1"):
                self.assertTrue(all(isinstance(c, float) for c in rec[key]))

    def test_the_base_draw_is_the_stratum_s_own_and_repeatable(self):
        a = charts.retract_inputs(SEN3_STRATA[3])
        b = charts.retract_inputs(SEN3_STRATA[3])
        c = charts.retract_inputs(SEN3_STRATA[4])
        self.assertEqual(a, b)
        self.assertNotEqual([r["q0"] for r in a], [r["q0"] for r in c])


if __name__ == "__main__":
    unittest.main()
