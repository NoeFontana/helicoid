"""The SO(2) and SE(2) ids against computations that share nothing with `gen.so2` and `gen.se2`:
the closed forms of docs/NUMERICS.md section 6 at high precision (the generator never evaluates
them), `mp.atan2`, the derivative of the matrix exponential (what J_r and J_l mean), hand-computed
layouts, the committed files, and planted errors for the generation-time checks."""

import json
import unittest

from mpmath import mp, mpf

from gen import check, manifest, precision, se2, so2, so3
from gen.fmt import Mat
from gen.registry import FUNCTIONS
from gen.strata import SCALAR_THETA_STRATA, SE2_SAMPLES, SE2_STRATA, SEN3_STRATA

precision.setup()
COMMITTED = manifest.ROOT.parent / "corpus"
NEAR_PI = 3.141592653589793 - 1e-9
TAUS = (
    [0.3, -0.5, 0.7],
    [0.01, 30.0, -7.0],
    [-2.0, 1.5, -1e3],
    [NEAR_PI, 2e3, 5.0],
    [-1e-8, 1.0, 2.0],
)
IDS = ("exp", "log", "ad", "jr", "jl", "jr_inv", "jl_inv")


def close(got, want, tol=None) -> bool:
    tol = mpf(10) ** -100 if tol is None else tol
    return all(abs(g - w) <= tol * max(1, abs(w)) for g, w in zip(got, want, strict=True))


def col_major(mat: Mat) -> list:
    return [[mat.data[3 * j + i] for j in range(3)] for i in range(3)]


def hat_of_column(J, j):
    """The hat matrix of the vector J e_j."""
    return mp.matrix(se2.hat([J[i][j] for i in range(3)]))


class ClosedFormTest(unittest.TestCase):
    """The definitions give NUMERICS section 6's Exp and Log, so a wrong line there fails this."""

    def test_exp_is_R_and_V_rho(self):
        for tau in TAUS:
            with mp.workdps(500):
                th, x, y = (mpf(c) for c in tau)
                a, b = mp.sin(th) / th, (1 - mp.cos(th)) / th
                want = [mp.cos(th), mp.sin(th), a * x - b * y, b * x + a * y]
            got = se2.exp({"tau": tau})
            self.assertTrue(close([*got["z"], *got["t"]], want), tau)

    def test_log_is_atan2_and_V_inverse_t(self):
        for tau in TAUS:
            x = se2.exp({"tau": tau})
            z, t = [precision.to_f64(c) for c in x["z"]], [precision.to_f64(c) for c in x["t"]]
            with mp.workdps(500):
                c, s = so3.normalize([mpf(v) for v in z])
                th = mp.atan2(s, c)
                a, b = mp.sin(th) / th, (1 - mp.cos(th)) / th
                want = [
                    th,
                    *mp.inverse(mp.matrix([[a, -b], [b, a]])) * mp.matrix([mpf(v) for v in t]),
                ]
            self.assertTrue(close(se2.log({"z": z, "t": t})["tau"], want, mpf(10) ** -95), tau)

    def test_angle_at_the_extremes(self):
        for z in (
            [1.0, 0.0],
            [1.0, 1e-310],
            [-1.0, 1e-12],
            [-1.0, -1e-12],
            [0.0, 1.0],
            [3.0, -4.0],
        ):
            with mp.workdps(500):
                want = mp.atan2(mpf(z[1]), mpf(z[0]))
            got = so2.angle([mpf(c) for c in z])
            self.assertLessEqual(abs(got - want), mpf(10) ** -100 * abs(want), z)
        with self.assertRaises(ValueError):  # atan2 splits +0 from -0 here; an mpf cannot
            so2.angle([mpf(-1), mpf(0)])

    def test_subnormal_theta_keeps_its_digits(self):
        """exp(i theta) = 1 + i theta to 1e-100 of each part, so the sine, 1e-310 beside a cosine of
        1, is kept to its own digits and not to an absolute 1e-100."""
        theta, tol = float.fromhex("0x0.031af1add4eb6p-1022"), mpf(10) ** -100
        for z in (so2.exp({"theta": theta})["z"], se2.exp({"tau": [theta, 1.0, 1.0]})["z"]):
            self.assertLess(abs(z[0] - 1), tol)
            self.assertLess(abs(z[1] - theta), tol * theta)
        self.assertLess(abs(so2.log({"z": [1.0, theta]})["theta"] - theta), tol * theta)


class JacobianTest(unittest.TestCase):
    def test_jr_and_jl_are_the_derivatives_of_exp(self):
        """d/de Exp(tau + e e_j) = Exp(tau) (J_r e_j)^ = (J_l e_j)^ Exp(tau): NUMERICS section 1."""
        h = mpf(10) ** -25
        with mp.workdps(90):
            for tau in TAUS[:4]:
                t, E = [mpf(c) for c in tau], mp.matrix(se2.expm([mpf(c) for c in tau]))
                Jr, Jl = (col_major(se2.jacobian(n)({"tau": tau})["J"]) for n in ("jr", "jl"))
                for j in range(3):
                    up, dn = ([t[k] + s * h * (k == j) for k in range(3)] for s in (1, -1))
                    D = (mp.matrix(se2.expm(up)) - mp.matrix(se2.expm(dn))) / (2 * h)
                    right, left = E * hat_of_column(Jr, j), hat_of_column(Jl, j) * E
                    for name, want in (("jr", right), ("jl", left)):
                        self.assertLess(
                            mp.mnorm(D - want, "inf"),
                            mpf(10) ** -35 * (1 + abs(t[1]) + abs(t[2])),
                            (name, tau, j),
                        )

    def test_hand_computed_layouts_at_zero_rotation(self):
        # tau = [0; (3, 5)]: ad has only its first column [0, 5, -3], ad^2 = 0, so J_l = I + ad/2
        # and J_l^-1 = I - ad/2; Ad of a pure translation (1, 2) has first column [1, 2, -1].
        first = {
            "jl": [1, 2.5, -1.5],
            "jr": [1, -2.5, 1.5],
            "jl_inv": [1, -2.5, 1.5],
            "jr_inv": [1, 2.5, -1.5],
        }
        for name, col in first.items():
            J = se2.jacobian(name)({"tau": [0.0, 3.0, 5.0]})["J"]
            self.assertEqual((J.shape, list(J.data)), ((3, 3), [*col, 0, 1, 0, 0, 0, 1]), name)
        Ad = se2.ad_of({"z": [1.0, 0.0], "t": [1.0, 2.0]})["Ad"]
        self.assertEqual(list(Ad.data), [1, 2, -1, 0, 1, 0, 0, 0, 1])

    def test_a_committed_record_has_that_layout(self):
        for name, scale in (("jl", 0.5), ("jr", -0.5), ("ad", None)):
            records = map(json.loads, (COMMITTED / f"se2_{name}.jsonl").read_text().splitlines())
            record = next(r for r in records if r["stratum"] == "theta:exact0")
            v = [float.fromhex(h) for h in record["in"].get("tau", record["in"].get("t"))[-2:]]
            want = [1, v[1] * (scale or 1), -v[0] * (scale or 1)] if scale else [1, v[1], -v[0]]
            got = [mpf(c) for c in record["out"].get("J", record["out"].get("Ad"))[:3]]
            self.assertLess(
                max(abs(g - w) for g, w in zip(got, want, strict=True)), mpf(10) ** -29, name
            )


class InputsTest(unittest.TestCase):
    def test_strata_and_sample_counts(self):
        self.assertIs(SE2_STRATA, SEN3_STRATA)
        self.assertEqual(len(SE2_STRATA), 52)
        for s in SE2_STRATA:
            self.assertEqual(len(se2.tau_inputs(s)), SE2_SAMPLES, s.name)
            self.assertEqual(len(se2.x_inputs(s)), SE2_SAMPLES, s.name)
        self.assertEqual(sum(len(so2.theta_inputs(s)) for s in SCALAR_THETA_STRATA), 3419)

    def test_every_theta_has_both_signs_in_so2_and_the_sign_alternates_in_se2(self):
        s = next(s for s in SCALAR_THETA_STRATA if s.name == "theta:1e-5")
        thetas = [r["theta"] for r in so2.theta_inputs(s)]
        self.assertEqual(thetas[1::2], [-t for t in thetas[::2]])
        self.assertEqual(
            [r["theta"] for r in so2.theta_inputs(SCALAR_THETA_STRATA[13])], [0.0]
        )  # exact0
        for stratum in (s, next(s for s in SE2_STRATA if s.name == "rho:1e3/theta=1")):
            signs = [r["tau"][0] > 0 for r in se2.tau_inputs(stratum)]
            self.assertEqual(signs, [i % 2 == 0 for i in range(SE2_SAMPLES)], stratum.name)
        exact0 = next(s for s in SE2_STRATA if s.name == "theta:exact0")
        self.assertEqual({repr(r["tau"][0]) for r in se2.tau_inputs(exact0)}, {"0.0"})  # never -0.0

    def test_records_share_theta_and_translation_and_the_translation_has_the_scale(self):
        for s in SE2_STRATA:
            for tau, x in zip(se2.tau_inputs(s), se2.x_inputs(s), strict=True):
                self.assertEqual(
                    (so2.z_of(tau["tau"][0]), tau["tau"][1:]), (x["z"], x["t"]), s.name
                )
                scale = mpf(10) ** (s.rho_exp or 0)
                self.assertLess(
                    abs(so3.norm([mpf(c) for c in x["t"]]) / scale - 1), mpf(10) ** -15, s.name
                )
        self.assertEqual(len({tuple(r["t"]) for r in se2.x_inputs(SE2_STRATA[0])}), SE2_SAMPLES)


def _entries(value):
    """The entries of a scalar, a vector or a matrix output, and how to put them back."""
    if isinstance(value, Mat):
        return list(value.data), lambda v: Mat(value.shape, v)
    return (list(value), list) if isinstance(value, list) else ([value], lambda v: v[0])


class CheckTest(unittest.TestCase):
    """A cross-check passes a correct record, fails an error of 1e-95 of an entry and passes 1e-105;
    an entry that must be exactly 0 (or is, at theta = 0) fails any error."""

    def cases(self):
        for tau in TAUS[:2]:
            x = se2.exp({"tau": tau})
            z, t = [precision.to_f64(c) for c in x["z"]], [precision.to_f64(c) for c in x["t"]]
            yield "so2_exp", {"theta": tau[0]}
            yield "so2_log", {"z": z}
            for name in IDS:
                yield f"se2_{name}", ({"z": z, "t": t} if name in ("log", "ad") else {"tau": tau})

    def test_every_id_has_its_own_check_and_it_detects_a_planted_error(self):
        for name, inputs in self.cases():
            spec = FUNCTIONS[name]
            good = spec.evaluate(inputs)
            spec.check(inputs, good)
            for key, value in good.items():
                entries, rebuild = _entries(value)
                for i, v in enumerate(entries):
                    for rel, ok in ((mpf(10) ** -95, False), (mpf(10) ** -105, True)):
                        bad = list(entries)
                        bad[i] = v * (1 + rel) if v else rel  # an entry that is 0 is exactly 0
                        if ok and v:
                            spec.check(inputs, {**good, key: rebuild(bad)})
                        else:
                            with self.assertRaises(check.CrossCheckError, msg=(name, key, i)):
                                spec.check(inputs, {**good, key: rebuild(bad)})

    def test_a_theta_a_full_turn_off_is_refused(self):
        """Exp is 2 pi periodic, so Exp(theta) = z cannot see the branch: only the range can.
        se2_log gets the rho that keeps Exp(tau) = X, so nothing else in its check objects."""
        for tau in TAUS[:2]:
            x = se2.exp({"tau": tau})
            z, t = [precision.to_f64(c) for c in x["z"]], [precision.to_f64(c) for c in x["t"]]
            theta = so2.log({"z": z})["theta"]
            for turn in (2, -2):
                bad = theta + turn * mp.pi
                Vinv = se2.inv(se2.translation_block(bad))
                rho = [sum((Vinv[i][k] * mpf(t[k]) for k in range(2)), mpf(0)) for i in range(2)]
                for name, inputs, out in (
                    ("so2_log", {"z": z}, {"theta": bad}),
                    ("se2_log", {"z": z, "t": t}, {"tau": [bad, *rho]}),
                ):
                    with self.assertRaisesRegex(check.CrossCheckError, "> pi", msg=(name, turn)):
                        FUNCTIONS[name].check(inputs, out)

    def test_translation_first_exp_is_refused(self):
        tau = TAUS[0]
        with self.assertRaises(check.CrossCheckError):
            FUNCTIONS["se2_exp"].check({"tau": tau}, se2.exp({"tau": tau[1:] + tau[:1]}))

    def test_the_wrong_sign_is_refused(self):
        for name, other in (("jr", "jl"), ("jl", "jr")):
            with self.assertRaises(check.CrossCheckError, msg=name):
                FUNCTIONS[f"se2_{name}"].check(
                    {"tau": TAUS[0]}, se2.jacobian(other)({"tau": TAUS[0]})
                )

    def test_a_tiny_entry_is_kept_to_100_digits(self):
        f = 1e-310  # ad has (1, 2) = -f and (2, 1) = f, and no product of ad's entries lands there
        J = se2.jacobian("jl")({"tau": [f, 1.0, 0.0]})["J"]
        self.assertLess(abs(J.data[2 * 3 + 1] + mpf(f) / 2), mpf(f) * mpf(10) ** -100)
        self.assertLess(abs(J.data[1 * 3 + 2] - mpf(f) / 2), mpf(f) * mpf(10) ** -100)


if __name__ == "__main__":
    unittest.main()
