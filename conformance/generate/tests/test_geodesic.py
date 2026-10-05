"""The geodesic ids against computations that share nothing with `gen.geodesic`: hand-computed
interpolations, the matrix route, GE.4's left invariance (which no generation-time check uses),
the strata's own properties, planted errors for the cross-checks, and the committed files."""

import json
import unittest

from mpmath import mp, mpf

from gen import check_geodesic, geodesic, manifest, precision, sen3, so3
from gen.check import CrossCheckError
from gen.registry import FUNCTIONS
from gen.strata import GEO_NEAR_PI_K, GEO_SAMPLES, GEO_STRATA

precision.setup()
COMMITTED = manifest.ROOT.parent / "corpus"
IDENTITY = [1.0, 0.0, 0.0, 0.0]
TOL = mpf(10) ** -100


def quarter_z():
    """The quarter turn about `z`, `(cos pi/4, 0, 0, sin pi/4)`, rounded to binary64."""
    return [float(mp.cos(mp.pi / 4)), 0.0, 0.0, float(mp.sin(mp.pi / 4))]


def near(got, want, scale=1):
    for g, w in zip(got, want, strict=True):
        if abs(mpf(g) - mpf(w)) > TOL * scale:
            raise AssertionError(f"{mp.nstr(mpf(g), 40)} vs {mp.nstr(mpf(w), 40)}")


def near30(got, want):
    """A committed reference is 30 significant decimal digits (`fmt.dec30`), so it is compared
    relative to its own size and not to the 1e-100 the cross-checks use."""
    for g, w in zip(got, want, strict=True):
        if abs(mpf(g) - mpf(w)) > mpf(10) ** -29 * max(abs(mpf(w)), mpf(10) ** -300):
            raise AssertionError(f"{mp.nstr(mpf(g), 40)} vs {mp.nstr(mpf(w), 40)}")


def strata(name):
    return {s.name: s for s in FUNCTIONS[name].strata}


class HandComputedTest(unittest.TestCase):
    def test_half_of_a_quarter_turn_is_an_eighth_turn(self):
        out = geodesic.so3_evaluate({"q0": IDENTITY, "q1": quarter_z(), "t": 0.5})
        near(out["q"], [mp.cos(mp.pi / 8), 0, 0, mp.sin(mp.pi / 8)])

    def test_a_pure_translation_interpolates_linearly(self):
        out = geodesic.se3_evaluate(
            {
                "q0": IDENTITY,
                "x0": [0.0, 0.0, 0.0],
                "q1": IDENTITY,
                "x1": [1.0, 2.0, 3.0],
                "t": 0.25,
            }
        )
        near(out["q"], IDENTITY)
        near(out["x"], [mpf(1) / 4, mpf(1) / 2, mpf(3) / 4])

    def test_a_screw_about_z_through_the_origin_halves_both_parts(self):
        """`Exp[(0, 0, pi/2); (0, 0, 1)]` is the quarter turn with `x = (0, 0, 1)`: `V` is the
        identity along the axis. At `t = 1/2` the rotation halves and so does the rise."""
        out = geodesic.se3_evaluate(
            {
                "q0": IDENTITY,
                "x0": [0.0, 0.0, 0.0],
                "q1": quarter_z(),
                "x1": [0.0, 0.0, 1.0],
                "t": 0.5,
            }
        )
        near(out["q"], [mp.cos(mp.pi / 8), 0, 0, mp.sin(mp.pi / 8)])
        near(out["x"], [0, 0, mpf(1) / 2])

    def test_a_quarter_turn_about_z_at_unit_offset_sweeps_the_circle(self):
        """`X_0 = (I, (1, 0, 0))`, `X_1` the quarter turn about the origin's `z` carrying it to
        `(0, 1, 0)`: the geodesic is the arc, so at `t = 1/2` the point is
        `(cos pi/4, sin pi/4, 0)` and the rotation is the eighth turn."""
        q1 = quarter_z()
        out = geodesic.se3_evaluate(
            {
                "q0": IDENTITY,
                "x0": [1.0, 0.0, 0.0],
                "q1": q1,
                "x1": [0.0, 1.0, 0.0],
                "t": 0.5,
            }
        )
        near(out["q"], [mp.cos(mp.pi / 8), 0, 0, mp.sin(mp.pi / 8)])
        near(out["x"], [mp.cos(mp.pi / 4), mp.sin(mp.pi / 4), 0])


class PropertyTest(unittest.TestCase):
    """Properties the generation-time checks do **not** use, so a defect common to the reference
    and its checks still shows: GE.4's two invariances, and `mp.expm` of the whole matrix."""

    def pairs(self):
        for name in ("so3_geodesic", "se3_geodesic"):
            spec = FUNCTIONS[name]
            for stratum in spec.strata:
                ins = spec.inputs(stratum)
                for i in (0, 3, 17, 42):
                    yield name, spec, ins[i]

    def test_left_invariance(self):
        """`gamma(G X_0, G X_1, t) = G gamma(X_0, X_1, t)` (GE.4(a)), `G` a fixed pose off the
        origin: `0045` item 5 measures the bound per `||t_G||`, and this is the 110-digit half."""
        g = so3.normalize([mpf(c) for c in (0.3, -0.5, 0.7, 0.4)])
        tg = [mpf(2), mpf(-3), mpf(5)]
        rg = so3.rot(g)
        for name, spec, inputs in self.pairs():
            with self.subTest(name=name, t=inputs["t"]):
                out = spec.evaluate(inputs)
                moved = dict(inputs)
                for q, x in (("q0", "x0"), ("q1", "x1")):
                    moved[q] = so3.qmul(g, so3.normalize([mpf(c) for c in inputs[q]]))
                    if x in inputs:
                        moved[x] = [
                            tg[r] + sum(rg[r][k] * mpf(inputs[x][k]) for k in range(3))
                            for r in range(3)
                        ]
                got = spec.evaluate(moved)
                want = {"q": so3.qmul(g, out["q"])}
                if "x" in out:
                    want["x"] = [
                        tg[r] + sum(rg[r][k] * out["x"][k] for k in range(3)) for r in range(3)
                    ]
                near(got["q"], check_geodesic._align(got["q"], want["q"]))
                if "x" in want:
                    near(got["x"], want["x"], scale=max(mpf(1), so3.maxabs(want["x"])))

    def test_the_matrix_route_agrees_where_mp_logm_is_the_principal_logarithm(self):
        """`X_0 expm(t logm(X_0^-1 X_1))`, the reference `PHASE4.md` §4 used to state. It is the
        cross-check and not the reference because it is wrong above `MP_LOGM_LIMIT`; below it the
        two agree, and `0045` item 1 is that distinction."""
        checked = 0
        for name, spec, inputs in self.pairs():
            theta = so3.norm(so3.log(geodesic.delta(inputs))["phi"])
            if theta > check_geodesic.MP_LOGM_LIMIT:
                continue
            out = spec.evaluate(inputs)
            check_geodesic._matrix_route(inputs, out, mpf(inputs["t"]), theta)
            checked += 1
        self.assertGreater(checked, 0)

    def test_mp_logm_is_not_the_principal_logarithm_on_geo_near_pi(self):
        """Why the reference cannot be `mp.logm` (`0045` item 1): on this stratum it returns a
        complex logarithm, and `_real` reports it as one rather than silently taking the real
        part."""
        spec = FUNCTIONS["se3_geodesic"]
        ins = spec.inputs(strata("se3_geodesic")["geo:near-pi"])
        complex_rows = 0
        for inputs in ins[::10]:
            m0 = check_geodesic._homogeneous(inputs["q0"], inputs["x0"])
            m1 = check_geodesic._homogeneous(inputs["q1"], inputs["x1"])
            rel = mp.matrix(sen3.mul(check_geodesic._inverse(m0), m1))
            try:
                check_geodesic._real(mp.logm(rel), "test")
            except CrossCheckError:
                complex_rows += 1
        self.assertGreater(complex_rows, 0)


class CrossCheckTest(unittest.TestCase):
    def test_every_perturbation_but_the_quaternion_s_sign_is_caught(self):
        """One component moved by `1e-95`, a swap of two, and `w` negated, at every stratum: each
        raises. Negating the **whole** quaternion does not, and must not -- `q` and `-q` are one
        rotation, which is the metric's `SignRule::Align` for this id as well."""
        tiny = mpf(10) ** -95
        for name in ("so3_geodesic", "se3_geodesic"):
            spec = FUNCTIONS[name]
            for stratum in spec.strata:
                ins = spec.inputs(stratum)
                for i in (0, 3, 5, 25):
                    inputs = ins[i]
                    out = spec.evaluate(inputs)
                    spec.check(inputs, out)
                    q = out["q"]
                    bad = {
                        "q moved": [q[0] + tiny] + list(q[1:]),
                        "q scaled": [c * (1 + tiny) for c in q],
                        "q w negated": [-q[0]] + list(q[1:]),
                        "q swapped": [q[0], q[2], q[1], q[3]],
                    }
                    for label, value in bad.items():
                        with self.subTest(name=name, stratum=stratum.name, i=i, bad=label):
                            with self.assertRaises(CrossCheckError):
                                spec.check(inputs, {**out, "q": value})
                    if "x" in out:
                        x = out["x"]
                        scale = max(mpf(1), so3.maxabs(x))
                        for label, value in {
                            "x moved": [x[0] + tiny * scale] + list(x[1:]),
                            "x scaled": [c * (1 + tiny) for c in x],
                            "x swapped": [x[1], x[0], x[2]],
                        }.items():
                            with self.subTest(name=name, stratum=stratum.name, i=i, bad=label):
                                with self.assertRaises(CrossCheckError):
                                    spec.check(inputs, {**out, "x": value})
                    spec.check(inputs, {**out, "q": [-c for c in q]})


class StrataTest(unittest.TestCase):
    def test_the_three_strata_are_phase4_section_4_s(self):
        self.assertEqual(
            [s.name for s in GEO_STRATA], ["geo:consecutive", "geo:generic", "geo:near-pi"]
        )
        for s in GEO_STRATA:
            self.assertEqual(s.count, GEO_SAMPLES)

    def test_geo_near_pi_stops_strictly_below_pi_at_its_stated_margins(self):
        """`0045` item 3: the stratum records its margin. The relative rotation of a record is
        `pi - 10^-k` up to the two roundings of `X_1`, which is 1e-16 against a smallest nominal
        margin of 1e-12."""
        spec = FUNCTIONS["so3_geodesic"]
        pairs = geodesic._pairs(strata("so3_geodesic")["geo:near-pi"])
        for k, (q0, _, q1, _) in zip(GEO_NEAR_PI_K, pairs, strict=True):
            theta = so3.norm(so3.log(geodesic.delta({"q0": q0, "q1": q1}))["phi"])
            margin = mp.pi - theta
            self.assertGreater(margin, 0)
            self.assertLess(abs(margin - mpf(10) ** -k), mpf(10) ** -15)
        del spec

    def test_geo_consecutive_holds_its_band_and_reaches_the_kilohertz_scales(self):
        """The relative rotation is in `[1e-9, 1e-3]` and the relative translation is at its own
        scale, so `||d||` is that band times `sqrt(2)`; `||x_0||` reaches 1e4.

        The relative translation is only that to the **quantization of the absolute poses**, and
        that is the regime, not a defect: a record stores `x_0` and `x_1`, so an increment of
        1e-9 at `||x_0|| = 1e4` is below the ulp of `x_0` and survives to about three digits. The
        bound is `||x_0|| u / theta` per unit, and the measured worst over the stratum is 1.02 of
        it (at `||x_0|| = 1`, where the quantization is the quaternion's), so the factor 4 here is
        room and not a fitted constant. This is exactly `tf_tree`'s kilohertz edge, which is why
        `PHASE4.md` §4 asks this stratum for `||t_0||` up to 1e4.
        """
        unit = mpf(2) ** -53
        pairs = geodesic._pairs(strata("se3_geodesic")["geo:consecutive"])
        norms = []
        for q0, x0, q1, x1 in pairs:
            tau = sen3.log(1)(geodesic.delta({"q0": q0, "x0": x0, "q1": q1, "x1": x1}))["tau"]
            theta = so3.norm(tau[:3])
            self.assertTrue(mpf(10) ** -9 <= theta < mpf(10) ** -3, mp.nstr(theta, 10))
            scale = max(mpf(1), so3.norm([mpf(c) for c in x0]))
            self.assertLess(abs(so3.norm(tau[3:]) / theta - 1), 4 * scale * unit / theta)
            norms.append(scale)
        self.assertAlmostEqual(float(max(norms)), 1e4, delta=1.0)
        self.assertAlmostEqual(float(min(norms)), 1.0, delta=1e-9)

    def test_the_t_values_are_phase4_section_4_s_six_plus_the_uniform_ones(self):
        ts = geodesic._ts(GEO_STRATA[0])
        self.assertEqual(len(ts), len(geodesic.GEO_T_FIXED) + geodesic.GEO_T_UNIFORM)
        self.assertEqual(ts[0], 0.0)
        self.assertEqual(ts[-1 - geodesic.GEO_T_UNIFORM], 1.0)
        self.assertTrue(all(0 < t < 1 for t in ts[-geodesic.GEO_T_UNIFORM :]))
        # Every sample of a stratum is crossed with the same `t` values, in the same order.
        self.assertEqual(ts, geodesic._ts(GEO_STRATA[0]))


class CommittedTest(unittest.TestCase):
    def test_the_committed_files_hold_what_the_generator_builds(self):
        """The first record of each stratum, re-evaluated: a file that drifted from the generator
        fails here as well as in `just corpus-check`."""
        for name in ("so3_geodesic", "se3_geodesic"):
            path = COMMITTED / f"{name}.jsonl"
            if not path.exists():  # before the first `just corpus`
                continue
            spec = FUNCTIONS[name]
            lines = path.read_text().splitlines()
            seen = set()
            for line in lines:
                record = json.loads(line)
                if record["stratum"] in seen:
                    continue
                seen.add(record["stratum"])
                inputs = {
                    k: (float.fromhex(v) if isinstance(v, str) else [float.fromhex(c) for c in v])
                    for k, v in record["in"].items()
                }
                out = spec.evaluate(inputs)
                for key, want in record["out"].items():
                    near30(out[key], [mpf(w) for w in want])
            self.assertEqual(seen, {s.name for s in spec.strata})


if __name__ == "__main__":
    unittest.main()
