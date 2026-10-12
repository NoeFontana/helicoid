"""The S2 ids (0063) against hand-computed points, the frame of CH.8 at n_z = +-0, the input draw,
and planted errors for the cross-checks."""

import unittest

from mpmath import mp, mpf

from gen import check_charts, corpus, fmt, precision, s2
from gen.check import CrossCheckError
from gen.precision import is_binary32
from gen.registry import FUNCTIONS

precision.setup()
TOL = mpf(10) ** -100


def near(got, want):
    for g, w in zip(got, want, strict=True):
        if abs(mpf(g) - mpf(w)) > TOL * max(1, abs(mpf(w))):
            raise AssertionError(f"{mp.nstr(mpf(g), 40)} vs {mp.nstr(mpf(w), 40)}")


def stratum(name: str, which: str = "s2_retract"):
    return next(s for s in FUNCTIONS[which].strata if s.name == name)


class BasisTest(unittest.TestCase):
    POINTS = (
        [0.6, 0.8, 0.0],
        [0.6, 0.8, -0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.48, -0.6, 0.64],
        [-0.48, 0.6, -0.64],
    )

    def test_the_frame_is_orthonormal_tangent_and_oriented(self):
        # CH.8(c), oriented (0082): B^T B = I, B^T n = 0, b1 x b2 = n on both hemispheres.
        for n in self.POINTS:
            nh, _, b1, b2 = s2.basis(n)
            near([s2.dot(b1, b1), s2.dot(b2, b2), s2.dot(b1, b2)], [1, 1, 0])
            near([s2.dot(b1, nh), s2.dot(b2, nh)], [0, 0])
            near(s2.cross(b1, b2), nh)

    def test_the_sign_of_zero_picks_no_frame(self):
        # 0063 decision 1: +0 and -0 are one point and one frame, sigma = +1.
        plus, minus = s2.basis([0.6, 0.8, 0.0]), s2.basis([0.6, 0.8, -0.0])
        self.assertEqual((plus[1], minus[1]), (1, 1))
        self.assertEqual(plus, minus)
        self.assertEqual(s2.basis([0.48, -0.6, -0.64])[1], -1)


class HandComputedTest(unittest.TestCase):
    def test_retract_at_the_pole_rotates_about_b1(self):
        # n = e_z: H = diag(1, 1, -1), B = [e_x e_y]; Exp(theta e_x) e_z = (0, -sin, cos).
        m = s2.retract({"n": [0.0, 0.0, 1.0], "delta": [1.0, 0.0]})["m"]
        near(m, [0, -mp.sin(1), mp.cos(1)])

    def test_local_is_half_of_r(self):
        # CH.9(c): n = e_z, m = (sin 1, 0, cos 1): n x m = (0, sin 1, 0), alpha = 1, so delta = (0, 1)
        # and alpha / s = 1.188395106, half of r(s^2, w) = 2.376790212.
        m = [float(mp.sin(1)), 0.0, float(mp.cos(1))]
        p, _, mh = s2.rotation_vector([0.0, 0.0, 1.0], m)
        self.assertAlmostEqual(float(p[1] / s2.norm(s2.cross([0, 0, 1], mh))), 1.188395106, 9)
        near(s2.local({"n": [0.0, 0.0, 1.0], "m": m})["delta"], [0, s2.norm(p)])

    def test_local_inverts_retract(self):
        for n in BasisTest.POINTS:
            for delta in ([1e-8, -2e-8], [0.3, -1.2], [-2.0, 2.2]):
                m = s2.retract({"n": n, "delta": delta})["m"]
                near(s2.local({"n": n, "m": m})["delta"], delta)


class DrawTest(unittest.TestCase):
    def test_every_delta_is_inside_pi_and_every_f32_input_is_binary32(self):
        for which in ("s2_retract", "s2_local"):
            for st in FUNCTIONS[which].strata:
                recs = st.records()
                self.assertEqual(len(recs), s2.N_S2, st.name)
                for r in recs:
                    if st.f32:
                        self.assertTrue(all(is_binary32(x) for v in r.values() for x in v))
                    if "delta" in r:
                        self.assertLess(s2.norm([mpf(c) for c in r["delta"]]), mp.pi)

    def test_near_antipode_reaches_its_margins_at_both_precisions(self):
        for name, k in (("s2:near-antipode", 12), ("s2:near-antipode@f32", 6)):
            thetas = [s2.norm([mpf(c) for c in r["delta"]]) for r in stratum(name).records()]
            self.assertLess(min(mp.pi - t for t in thetas), mpf(10) ** -(k - 1), name)

    def test_nz0_alternates_the_sign_of_zero_and_the_hex_keeps_it(self):
        for name in ("s2:nz0", "s2:nz0@f32"):
            recs = stratum(name).records()
            zs = [r["n"][2] for r in recs[:4]]
            self.assertEqual([mp.sign(z) for z in zs], [0] * 4)
            self.assertEqual([str(z) for z in zs], ["0.0", "-0.0", "0.0", "-0.0"])
            text = corpus._encode(recs[1], fmt.hex_float)["n"][2]
            self.assertEqual(text, "-0x0.0p+0")
            self.assertEqual(str(float.fromhex(text)), "-0.0")

    def test_local_reads_the_rounded_retract_of_the_same_draw(self):
        for name in ("s2:generic", "s2:near-antipode@f32"):
            ret, loc = stratum(name).records(), stratum(name, "s2_local").records()
            for a, b in zip(ret, loc, strict=True):
                self.assertEqual(a["n"], b["n"])
                m = s2.retract(a)["m"]
                rnd = precision.to_f32 if name.endswith("@f32") else precision.to_f64
                self.assertEqual(b["m"], [rnd(c) for c in m])


class PlantedTest(unittest.TestCase):
    REC = {"n": [0.6, 0.8, -0.0], "delta": [0.4, -1.1]}

    def test_retract_check_catches_the_sign_bit_reading(self):
        # `copysign(1, n_z)` gives sigma = -1 at -0: the frame B_- = B_+ T (CH.11(d)).
        nh, _, _, _ = s2.basis(self.REC["n"])
        nu = [nh[0], nh[1], nh[2] - 1]
        k = 2 / s2.dot(nu, nu)
        b1 = [(1 if i == 0 else 0) - k * nu[i] * nu[0] for i in range(3)]
        b2 = [(1 if i == 1 else 0) - k * nu[i] * nu[1] for i in range(3)]
        d1, d2 = (mpf(c) for c in self.REC["delta"])
        wrong = s2.rodrigues([d1 * a + d2 * b for a, b in zip(b1, b2, strict=True)], nh)
        check_charts.s2_retract(self.REC, s2.retract(self.REC))
        with self.assertRaises(CrossCheckError):
            check_charts.s2_retract(self.REC, {"m": wrong})

    def test_local_check_catches_the_factor_of_two(self):
        # NUMERICS section 8 as first worded: r(s^2, w) where alpha / s is half of it.
        rec = {"n": self.REC["n"], "m": s2.retract(self.REC)["m"]}
        out = s2.local(rec)
        check_charts.s2_local(rec, out)
        with self.assertRaises(CrossCheckError):
            check_charts.s2_local(rec, {"delta": [2 * c for c in out["delta"]]})


if __name__ == "__main__":
    unittest.main()
