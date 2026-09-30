"""The SO(3) ids against computations that share nothing with `gen.so3`: Rodrigues and the closed
forms of the Jacobians at high precision, `mp.logm`, a sandwich product, exact polar factors."""

import math
import unittest
from unittest import mock

from mpmath import mp, mpf

from gen import check, precision, so3
from gen.fmt import Mat
from gen.registry import FUNCTIONS
from gen.strata import QUAT_STRATA, SCALAR_THETA_STRATA

precision.setup()
REL = mpf(10) ** -100
NEAR_PI = 3.141592653589793 - 1e-9
PHIS = (
    [0.0, 0.0, 0.0],
    [1e-310, 2e-310, -1e-310],
    [3e-9, -4e-9, 1e-9],
    [0.3, -0.5, 0.7],
    [2.0, 1.5, -1.9],
    [0.6 * NEAR_PI, 0.0, 0.8 * NEAR_PI],  # a zero component stays exactly zero
)
BY_NAME = {s.name: s for s in QUAT_STRATA}


def norm(v):
    return mp.sqrt(sum(c * c for c in v))


def digits_for(phi) -> int:
    """Guard digits for a closed form that cancels like theta^-4 (`e`: theta^-5 over theta)."""
    with mp.workdps(50):
        t = norm([mpf(c) for c in phi])
    return 400 + (4 * max(0, -int(mp.floor(mp.log10(t)))) if t else 0)


def hat(v):
    x, y, z = v
    return mp.matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])


def close(got, want, tol=REL) -> bool:
    return abs(got - want) <= tol * abs(want)


class ExpAndJacobianTest(unittest.TestCase):
    def test_exp_is_rodrigues(self):
        for phi in PHIS:
            got = so3.exp({"phi": phi})["q"]
            with mp.workdps(digits_for(phi)):
                p = [mpf(c) for c in phi]
                t = norm(p)
                want = [mp.cos(t / 2)] + [(mp.sin(t / 2) / t if t else 0) * c for c in p]
                self.assertTrue(all(close(g, w) for g, w in zip(got, want, strict=True)), phi)

    def test_jacobians_are_their_closed_forms(self):
        names = ("jr", "jl", "jr_inv", "jl_inv")
        for phi in PHIS:
            got = {n: so3.jacobian(n)({"phi": phi})["J"] for n in names}
            with mp.workdps(digits_for(phi)):
                t, W, eye = norm([mpf(c) for c in phi]), hat([mpf(c) for c in phi]), mp.eye(3)
                a, b, c = (
                    (0, 0, 0)
                    if t == 0
                    else (
                        (1 - mp.cos(t)) / t**2,
                        (t - mp.sin(t)) / t**3,
                        1 / t**2 - (1 + mp.cos(t)) / (2 * t * mp.sin(t)),
                    )
                )
                want = {
                    "jr": eye - a * W + b * W * W,
                    "jl": eye + a * W + b * W * W,
                    "jr_inv": eye + W / 2 + c * W * W,
                    "jl_inv": eye - W / 2 + c * W * W,
                }
                for name in names:
                    self.assertEqual(got[name].shape, (3, 3))
                    for j in range(3):  # column-major
                        for i in range(3):
                            self.assertTrue(
                                close(got[name].data[3 * j + i], want[name][i, j]),
                                (name, phi, i, j),
                            )

    def test_a_series_keeps_a_tiny_entry_to_100_digits(self):
        J = so3.jacobian("jl")({"phi": [1e-310, 0.0, 0.0]})["J"]  # I + W / 2: (2, 1) = +phi_x / 2
        self.assertTrue(close(J.data[3 * 1 + 2], mpf(1e-310) / 2))


class LogTest(unittest.TestCase):
    def test_log_is_the_inverse_of_exp_on_the_canonical_branch(self):
        for phi in PHIS[:5]:
            q = [precision.to_f64(c) for c in so3.exp({"phi": phi})["q"]]
            with mp.workdps(500):
                v = [mpf(c) for c in q[1:]]
                t = 2 * mp.atan2(norm(v), mpf(q[0])) / norm(v) if norm(v) else 0
                want = [t * c for c in v]
            variants = (q, [-c for c in q], [2 * c for c in q])  # the flip and the scale invariance
            for signed in variants:
                got = so3.log({"q": signed})["phi"]
                self.assertTrue(
                    all(close(g, w) for g, w in zip(got, want, strict=True)), (phi, signed)
                )

    def test_log_matches_logm(self):
        for phi in (PHIS[2], PHIS[3], [1.2, -0.7, 1.4]):  # `mp.logm` is complex near pi
            q = [precision.to_f64(c) for c in so3.exp({"phi": phi})["q"]]
            got = so3.log({"q": q})["phi"]
            L = mp.logm(mp.matrix(so3.rot(so3.normalize([mpf(c) for c in q]))))
            for g, w in zip(got, (L[2, 1], L[0, 2], L[1, 0]), strict=True):
                self.assertLess(abs(g - w), REL * norm(got), phi)

    def test_w_zero_is_the_function_of_the_quaternion(self):
        for signs in (1, -1):
            q = [0.0, 0.6 * signs, 0.0, 0.8 * signs]
            got = so3.log({"q": q})["phi"]
            self.assertLess(abs(norm(got) - mp.pi), REL)
            self.assertEqual(got[1], 0)
            v = [mpf(c) for c in q[1:]]  # phi = pi u for the unit u of vec q, sign and all
            for g, c in zip(got, v, strict=True):
                self.assertLess(abs(g - mp.pi * c / norm(v)), REL)

    def test_the_identity_and_its_negative_are_zero_exactly(self):
        for q in ([1.0, 0.0, 0.0, 0.0], [-1.0, -0.0, -0.0, -0.0]):
            self.assertEqual(so3.log({"q": q})["phi"], (0, 0, 0))


class ActAndMatrixTest(unittest.TestCase):
    def test_act_is_the_sandwich_of_the_normalized_quaternion(self):
        q, p = [0.3, -1.1, 0.4, 2.0], [0.2, 0.5, -0.7]
        with mp.workdps(300):
            n = norm([mpf(c) for c in q])
            a, b = [mpf(c) / n for c in q], [mpf(0)] + [mpf(c) for c in p]
            want = so3.qmul(so3.qmul(a, b), so3.qconj(a))[1:]
        for g, w in zip(so3.act({"q": q, "p": p})["Rp"], want, strict=True):
            self.assertLess(abs(g - w), REL)

    def test_to_mat_is_column_major(self):
        A = [[mpf(3 * i + j) for j in range(3)] for i in range(3)]
        self.assertEqual(so3.to_mat(A), Mat((3, 3), [0, 3, 6, 1, 4, 7, 2, 5, 8]))
        self.assertEqual(so3.from_mat(so3.to_mat(A)), A)

    def test_from_matrix_returns_the_polar_factor_of_a_perturbed_rotation(self):
        # R has entries in {0, +-1} and H is dyadic: A = R H is exact, and R is its polar factor.
        R = so3.rot([mpf(0.5), mpf(0.5), mpf(0.5), mpf(-0.5)])
        H = [[1.25, 0.125, 0.0], [0.125, 0.75, 0.0625], [0.0, 0.0625, 1.0]]
        A = so3.mm(R, [[mpf(x) for x in row] for row in H])
        self.assertTrue(all(x == float(x) for x in so3.entries(A)))
        got = so3.from_matrix({"R": so3.to_mat(A)})["q"]
        for g, w in zip(got, (0.5, 0.5, 0.5, -0.5), strict=True):
            self.assertLess(abs(g - w), REL)
        scaled = so3.to_mat([[2 * x for x in row] for row in R])  # a scaled rotation: same factor
        self.assertLess(abs(so3.from_matrix({"R": scaled})["q"][0] - mpf(0.5)), REL)

    def test_from_matrix_at_pi_has_w_zero_and_a_canonical_sign(self):
        for signs in ((1, -1, -1), (-1, 1, -1), (-1, -1, 1)):
            A = Mat(
                (3, 3), [float(s) if i == j else 0.0 for j, s in enumerate(signs) for i in range(3)]
            )
            q = so3.from_matrix({"R": A})["q"]
            axis = signs.index(1)
            self.assertEqual([c == 0 for c in q], [True] + [i != axis for i in range(3)])
            self.assertLess(abs(q[1 + axis] - 1), REL)

    def test_from_matrix_of_identity_plus_a_subnormal_skew_part(self):
        s = [1e-310, -2e-310, 3e-310]
        A = so3.to_mat([[1, -s[2], s[1]], [s[2], 1, -s[0]], [-s[1], s[0], 1]])
        A = Mat((3, 3), [float(x) for x in A.data])
        q = so3.from_matrix({"R": A})["q"]
        self.assertEqual(q[0], 1)
        for g, c in zip(q[1:], s, strict=True):
            self.assertTrue(close(g, mpf(c) / 2))


class CheckTest(unittest.TestCase):
    """Each SO(3) cross-check passes a correct record and fails one off by 1e-95, not by 1e-105."""

    def cases(self):
        phi = [0.3, -0.5, 0.7]
        q = [precision.to_f64(c) for c in so3.exp({"phi": phi})["q"]]
        A = so3.to_mat(so3.rot(so3.normalize([mpf(c) for c in q])))
        A = Mat((3, 3), [precision.to_f64(x) for x in A.data])
        yield "so3_exp", {"phi": phi}
        yield "so3_log", {"q": q}
        yield "so3_act", {"q": q, "p": [0.2, 0.5, -0.7]}
        yield "so3_from_matrix", {"R": A}
        for name in ("jr", "jl", "jr_inv", "jl_inv"):
            yield f"so3_{name}", {"phi": phi}

    def test_every_id_carries_its_own_check_and_it_detects_a_planted_error(self):
        cases = list(self.cases())
        self.assertEqual({n for n, _ in cases}, {n for n in FUNCTIONS if n.startswith("so3_")})
        for name, inputs in cases:
            spec = FUNCTIONS[name]
            good = spec.evaluate(inputs)
            spec.check(inputs, good)
            key = next(iter(good))
            values = list(good[key].data if isinstance(good[key], Mat) else good[key])
            for i in range(len(values)):
                for size, ok in ((mpf(10) ** -95, False), (mpf(10) ** -105, True)):
                    bad = list(values)
                    bad[i] = values[i] + size
                    out = {key: Mat((3, 3), bad) if isinstance(good[key], Mat) else bad}
                    if ok:
                        spec.check(inputs, out)
                    else:
                        with self.assertRaises(check.CrossCheckError, msg=(name, i)):
                            spec.check(inputs, out)

    def test_a_tiny_vector_part_is_checked_to_its_own_size(self):
        inputs = {"phi": [1e-310, 2e-310, -1e-310]}
        good = so3.exp(inputs)["q"]
        check.so3_exp(inputs, {"q": good})
        for i in (1, 2, 3):
            bad = list(good)
            bad[i] *= 1 + mpf(10) ** -95
            with self.assertRaises(check.CrossCheckError, msg=i):
                check.so3_exp(inputs, {"q": bad})

    def test_a_tiny_off_diagonal_part_is_checked_to_its_own_size(self):
        """At theta:subnormal a Jacobian is I plus 1e-310 and `Log` returns phi of that size: the
        wrong side, or an error of 1e-95 of an entry (not 1e-105), fails, off the diagonal."""
        inputs = {"phi": [1e-310, 2e-310, -1e-310]}
        other = {"jr": "jl", "jl": "jr", "jr_inv": "jl_inv", "jl_inv": "jr_inv"}
        good = {n: so3.jacobian(n)(inputs)["J"] for n in other}
        for name, spec in ((n, FUNCTIONS[f"so3_{n}"]) for n in other):
            spec.check(inputs, {"J": good[name]})
            with self.assertRaises(check.CrossCheckError, msg=name):
                spec.check(inputs, {"J": good[other[name]]})
            for i in (1, 2, 3, 5, 6, 7):  # column-major, off the diagonal
                for size, ok in ((mpf(10) ** -95, False), (mpf(10) ** -105, True)):
                    data = list(good[name].data)
                    data[i] *= 1 + size
                    if ok:
                        spec.check(inputs, {"J": Mat((3, 3), data)})
                    else:
                        with self.assertRaises(check.CrossCheckError, msg=(name, i)):
                            spec.check(inputs, {"J": Mat((3, 3), data)})
        q = [precision.to_f64(c) for c in so3.exp(inputs)["q"]]
        phi = so3.log({"q": q})["phi"]
        check.so3_log({"q": q}, {"phi": phi})
        for i in range(3):
            bad = list(phi)
            bad[i] *= 1 + mpf(10) ** -95
            with self.assertRaisesRegex(check.CrossCheckError, "exp residual:", msg=i):
                check.so3_log({"q": q}, {"phi": bad})

    def test_log_off_the_canonical_branch_or_sign_is_refused(self):
        inputs = {"q": [0.5, 0.5, 0.5, 0.5]}
        good = so3.log(inputs)["phi"]
        check.so3_log(inputs, {"phi": good})
        t = norm(good)
        wrapped = [c * (t - 2 * mp.pi) / t for c in good]  # Exp is -q: the rotation, not the q
        with self.assertRaisesRegex(check.CrossCheckError, "exp residual:"):
            check.so3_log(inputs, {"phi": [-c for c in good]})
        with self.assertRaisesRegex(check.CrossCheckError, "exp residual w:"):
            check.so3_log(inputs, {"phi": wrapped})

    def test_log_beyond_pi_is_refused_when_only_the_bound_fails(self):
        inputs = {"q": [0.5, 0.5, 0.5, 0.5]}
        good = so3.log(inputs)["phi"]
        t = norm(good)
        long = [c * (t + 4 * mp.pi) / t for c in good]  # two more turns: same q, same R, |phi| > pi
        with self.assertRaisesRegex(check.CrossCheckError, r"\|phi\| > pi"):
            check.so3_log(inputs, {"phi": long})

    def test_from_matrix_with_a_half_turned_polar_factor_is_refused(self):
        """A = R H, H = diag(2, 3, 5); R(q') = R Rot(pi, x) has H' = diag(2, -3, -5): symmetric,
        det > 0, indefinite. Only the positive-definiteness test refuses it."""
        q = [0.5, 0.5, 0.5, -0.5]
        R = so3.rot([mpf(c) for c in q])
        A = so3.to_mat(so3.mm(R, [[mpf(2), 0, 0], [0, mpf(3), 0], [0, 0, mpf(5)]]))
        inputs = {"R": Mat((3, 3), [float(x) for x in A.data])}
        check.so3_from_matrix(inputs, {"q": so3.from_matrix(inputs)["q"]})
        q_half_turned = so3.qmul([mpf(c) for c in q], [mpf(0), mpf(1), mpf(0), mpf(0)])
        with self.assertRaisesRegex(check.CrossCheckError, "not positive definite"):
            check.so3_from_matrix(inputs, {"q": q_half_turned})

    def test_a_reflection_is_not_a_rotation(self):
        turn = so3.rot
        with mock.patch.object(so3, "rot", lambda q: [[-x for x in row] for row in turn(q)]):
            with self.assertRaisesRegex(check.CrossCheckError, "det R"):
                check.so3_act({"q": [1.0, 0.0, 0.0, 0.0], "p": [0.0, 0.0, 1.0]}, {"Rp": [0, 0, 1]})


class InputsTest(unittest.TestCase):
    def test_phi_inputs_carry_an_axis_per_sample(self):
        by = {s.name: s for s in SCALAR_THETA_STRATA}
        self.assertEqual(so3.phi_inputs(by["theta:exact0"]), [{"phi": [0.0, 0.0, 0.0]}])
        recs = so3.phi_inputs(by["theta:pi-1e-3"])
        self.assertEqual(len(recs), 64)
        self.assertEqual(len({tuple(r["phi"]) for r in recs}), 64)
        for r in recs:
            self.assertLess(abs(norm([mpf(c) for c in r["phi"]]) - (mp.pi - mpf(10) ** -3)), 3e-15)
        self.assertEqual(len(so3.phi_inputs(by["theta:dense"])), 801)

    def test_log_inputs_pair_every_quaternion_with_its_negative_except_at_w_zero(self):
        recs = so3.log_inputs(BY_NAME["theta:1e-3"])
        self.assertEqual(len(recs), 128)
        for q, n in zip(recs[::2], recs[1::2], strict=True):
            self.assertEqual(n["q"], [-c for c in q["q"]])
        self.assertEqual(len(so3.log_inputs(BY_NAME["q:w0"])), 64)
        self.assertEqual(len(so3.log_inputs(BY_NAME["theta:exact0"])), 2)

    def test_quaternion_strata(self):
        for q in so3.quaternions(BY_NAME["q:w0"]):
            self.assertEqual(q[0], 0.0)
            self.assertEqual(math.copysign(1, q[0]), 1)  # +0, never -0
            self.assertLess(abs(norm([mpf(c) for c in q]) - 1), 1e-15)
        for i, q in enumerate(so3.quaternions(BY_NAME["q:nonunit"])):
            n2 = sum(mpf(c) ** 2 for c in q)
            self.assertLess(abs(n2 - 1 - (-1) ** i * mpf(2) ** -45), 1e-15, i)  # +-2^-45, rounded

    def test_act_and_matrix_inputs(self):
        stratum = BY_NAME["theta:pi-1e-3"]
        for r in so3.act_inputs(stratum):
            self.assertLess(abs(norm([mpf(c) for c in r["p"]]) - 1), 1e-15)
        for r, q in zip(so3.matrix_inputs(stratum), so3.quaternions(stratum), strict=True):
            self.assertEqual(r["R"].shape, (3, 3))
            with mp.workdps(30):  # a rounded rotation: orthogonal to 1e-15
                R = so3.from_mat(r["R"])
                err = so3.entries(so3.mm(so3.mt(R), R))
                self.assertLess(max(abs(e - int(i % 4 == 0)) for i, e in enumerate(err)), 1e-15)


if __name__ == "__main__":
    unittest.main()
