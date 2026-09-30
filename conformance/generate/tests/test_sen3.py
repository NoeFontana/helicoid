"""The SE_N(3) ids against computations that share nothing with `gen.sen3`: the block forms of
docs/NUMERICS.md section 5 at high precision (the generator never evaluates them), `mp.logm`,
hand-computed layouts, the committed files, and planted errors for the generation-time checks."""

import json
import math
import unittest

from mpmath import mp, mpf

from gen import check, check_sen3, manifest, precision, sen3, so3
from gen.fmt import Mat
from gen.registry import FUNCTIONS
from gen.strata import SEN3_QUAT_STRATA, SEN3_SAMPLES, SEN3_STRATA

precision.setup()
COMMITTED = manifest.ROOT.parent / "corpus"
NEAR_PI = 3.141592653589793 - 1e-9
PHI = ([0.3, -0.5, 0.7], [0.01, 0.02, -0.01], [2.0, 1.5, -1.9], [0.6 * NEAR_PI, 0.0, 0.8 * NEAR_PI])
RHO = ([1.0, -2.0, 0.5], [30.0, 0.25, -7.0], [-1e3, 2e3, 5.0])  # one per rho_i, up to N = 3
IDS = ("exp", "log", "ad", "jr", "jl", "jr_inv", "jl_inv")


def taus(n):
    return [p + [c for r in RHO[:n] for c in r] for p in PHI]


def hat(v):
    return mp.matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])


def blocks(A, B, n):
    """The dense [[A, 0], [B_i, A]] of a dual matrix, as a nested list."""
    m = 3 + 3 * n
    M = mp.zeros(m, m)
    for i in range(n + 1):
        M[3 * i : 3 * i + 3, 3 * i : 3 * i + 3] = A
        if i:
            M[3 * i : 3 * i + 3, 0:3] = B[i - 1]
    return M


def closed_forms(phi, rho):
    """A_r, A_l, A_r^-1, A_l^-1 and Q(rho_i, phi), Q(-rho_i, -phi): NUMERICS sections 3.5, 5.3."""
    t, W, I3 = mp.sqrt(sum(c * c for c in phi)), hat(phi), mp.eye(3)
    a, b = (1 - mp.cos(t)) / t**2, (t - mp.sin(t)) / t**3
    c = 1 / t**2 - (1 + mp.cos(t)) / (2 * t * mp.sin(t))
    d = (t * t + 2 * mp.cos(t) - 2) / (2 * t**4)
    e = (2 * t - 3 * mp.sin(t) + t * mp.cos(t)) / (2 * t**5)

    def Q(r, W):
        P = hat(r)
        return (
            P / 2
            + b * (W * P + P * W + W * P * W)
            + d * (W * W * P + P * W * W - 3 * W * P * W)
            + e * (W * P * W * W + W * W * P * W)
        )

    return {
        "Ar": I3 - a * W + b * W * W,
        "Al": I3 + a * W + b * W * W,
        "Ari": I3 + W / 2 + c * W * W,
        "Ali": I3 - W / 2 + c * W * W,
        "Ql": [Q(r, W) for r in rho],
        "Qr": [Q([-x for x in r], -W) for r in rho],
    }


def dense(mat) -> mp.matrix:
    m = mat.shape[0]
    return mp.matrix([[mat.data[j * m + i] for j in range(m)] for i in range(m)])


def close(got, want, tol=None) -> bool:
    tol = mpf(10) ** -100 if tol is None else tol
    return mp.mnorm(got - want, "inf") <= tol * max(1, mp.mnorm(want, "inf"))


class BlockFormsTest(unittest.TestCase):
    """The definitions give NUMERICS section 5's blocks, so a wrong line there fails this test."""

    def test_jacobians_and_their_dual_inverses(self):
        for n in (1, 2, 3):
            for tau in taus(n):
                with mp.workdps(500):
                    phi, rho = [mpf(c) for c in tau[:3]], sen3.chunks([mpf(c) for c in tau[3:]], n)
                    f = closed_forms(phi, rho)
                    want = {
                        "jl": blocks(f["Al"], f["Ql"], n),
                        "jr": blocks(f["Ar"], f["Qr"], n),
                        # section 5.4: A^-1 - eps A^-1 B_i A^-1
                        "jr_inv": blocks(f["Ari"], [-f["Ari"] * q * f["Ari"] for q in f["Qr"]], n),
                        "jl_inv": blocks(f["Ali"], [-f["Ali"] * q * f["Ali"] for q in f["Ql"]], n),
                    }
                    for name, w in want.items():
                        got = dense(sen3.jacobian(n, name)({"tau": tau})["J"])
                        self.assertTrue(close(got, w), (n, name, tau))

    def test_exp_and_the_adjoint(self):
        for n in (1, 2, 3):
            for tau in taus(n):
                with mp.workdps(500):
                    phi, rho = [mpf(c) for c in tau[:3]], sen3.chunks([mpf(c) for c in tau[3:]], n)
                    f, t = closed_forms(phi, rho), mp.sqrt(sum(c * c for c in phi))
                    got = sen3.exp(n)({"tau": tau})
                    q = [mp.cos(t / 2)] + [mp.sin(t / 2) / t * c for c in phi]
                    x = [f["Al"] * mp.matrix(r) for r in rho]  # section 5.1: x_i = J_l(phi) rho_i
                    want = [*q, *(v for xi in x for v in xi)]
                    self.assertTrue(close(mp.matrix([*got["q"], *got["x"]]), mp.matrix(want)))
                    # section 5.2, Ad_X = R + eps [x_i]x R, of the binary64 X the record holds
                    inputs = {
                        "q": [precision.to_f64(c) for c in q],
                        "x": [precision.to_f64(v) for xi in x for v in xi],
                    }
                    R = mp.matrix(so3.rot(so3.normalize([mpf(c) for c in inputs["q"]])))
                    xs = sen3.chunks([mpf(c) for c in inputs["x"]], n)
                    Ad = blocks(R, [hat(xi) * R for xi in xs], n)
                    self.assertTrue(close(dense(sen3.ad_of(n)(inputs)["Ad"]), Ad), (n, tau))


class LayoutTest(unittest.TestCase):
    """Column-major, rotation-first: column j is the image of e_j, rows [phi; rho_1; ...]."""

    def test_hand_computed_jacobians_at_zero_rotation(self):
        # tau = [0; (1, 2, 3)]: ad_tau^2 = 0, J = I -+ ad_tau / 2; the lower-left block of ad_tau
        # is rho^ = [[0,-3,2],[3,0,-1],[-2,1,0]], and the layout is column-major.
        rho_hat = [[0, -3, 2], [3, 0, -1], [-2, 1, 0]]
        eye = [1 if i == j else 0 for j in range(6) for i in range(6)]
        for name, sign in (("jr", -1), ("jl", 1)):
            want = list(eye)
            for i in range(3):
                for j in range(3):
                    want[j * 6 + 3 + i] = mpf(sign * rho_hat[i][j]) / 2
            got = sen3.jacobian(1, name)({"tau": [0.0] * 3 + [1.0, 2.0, 3.0]})["J"]
            self.assertEqual((got.shape, list(got.data)), ((6, 6), want), name)

    def test_hand_computed_adjoint_of_a_translation(self):
        # X = (I, x), x = (1, 2, 3): Ad_X = [[I, 0], [x^, I]], x^ = [[0,-3,2],[3,0,-1],[-2,1,0]].
        got = sen3.ad_of(1)({"q": [1.0, 0.0, 0.0, 0.0], "x": [1.0, 2.0, 3.0]})["Ad"]
        col = lambda j: [got.data[6 * j + i] for i in range(6)]
        self.assertEqual(col(0), [1, 0, 0, 0, 3, -2])
        self.assertEqual(col(1), [0, 1, 0, -3, 0, 1])
        self.assertEqual(col(2), [0, 0, 1, 2, -1, 0])
        self.assertEqual(
            [col(j) for j in (3, 4, 5)], [[0] * i + [1] + [0] * (5 - i) for i in (3, 4, 5)]
        )

    def test_a_committed_record_has_that_layout(self):
        """The lower-left block of the theta:exact0 records of `sen3_jr_n1`, `sen3_jl_n1` and the
        adjoint, read as hex floats and strings from the committed files: -+rho^/2 and x^."""
        for name, scale, key in (("jr", -0.5, "J"), ("jl", 0.5, "J"), ("ad", 1.0, "Ad")):
            lines = (COMMITTED / f"sen3_{name}_n1.jsonl").read_text().splitlines()
            record = next(r for r in map(json.loads, lines) if r["stratum"] == "theta:exact0")
            v = [float.fromhex(h) for h in record["in"].get("tau", record["in"].get("x", []))[-3:]]
            hat_v = [[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]
            self.assertEqual(record["out"]["shape"], [6, 6])
            for i in range(3):
                for j in range(3):
                    got = mpf(record["out"][key][j * 6 + 3 + i])
                    self.assertLess(
                        abs(got - scale * mpf(hat_v[i][j])), mpf(10) ** -29, (name, i, j)
                    )
                    self.assertEqual(mpf(record["out"][key][j * 6 + i]), int(i == j))


class TinyTest(unittest.TestCase):
    def test_a_series_keeps_a_tiny_entry_to_100_digits(self):
        """phi = (f, 0, 0), rho = (0, 1, 0), f = 1e-310: rho^ has no (0, 1) entry, so that entry of
        the Q block of J_l is (phi^ rho^ + rho^ phi^) / 6 = f / 6 alone."""
        f, tau = 1e-310, [1e-310, 0.0, 0.0, 0.0, 1.0, 0.0]
        J = sen3.jacobian(1, "jl")({"tau": tau})["J"]
        self.assertLess(abs(J.data[1 * 6 + 3] - mpf(f) / 6), mpf(f) * mpf(10) ** -100)
        self.assertLess(abs(J.data[1 * 6 + 2] - mpf(f) / 2), mpf(f) * mpf(10) ** -100)  # A: (2, 1)

    def test_a_dense_inverse_zeroes_the_structural_blocks_and_keeps_every_other_entry(self):
        J = sen3.series([mpf(c) for c in taus(3)[2]], 3, -1)
        raw, Jinv = mp.inverse(mp.matrix(J)), sen3.inverse(J)
        cells = [(i, j) for i in range(12) for j in range(12)]
        zero = {(i, j) for i, j in cells if sen3.structural_zero(i // 3, j // 3)}
        self.assertTrue(any(raw[i, j] != 0 for i, j in zero))  # the LU leaves residue here
        for i, j in cells:
            self.assertEqual(Jinv[i][j], 0 if (i, j) in zero else raw[i, j], (i, j))
        self.assertEqual(sen3.inverse(sen3.eye(3))[0][1], 0)
        # no magnitude floor: at phi = (f, 0, 0), f = 1e-310, J_r^-1 = I + phi^ / 2 + O(f^2), so its
        # (2, 1) entry is f / 2 and the (1, 2) entry -f / 2, next to the rho block's entries of 1
        f = 1e-310
        tiny = sen3.inverse(sen3.series([mpf(c) for c in [f, 0, 0, 0, 1, 0]], 1, -1))
        for got, want in ((tiny[2][1], mpf(f) / 2), (tiny[1][2], -mpf(f) / 2)):
            self.assertLess(abs(got - want), mpf(f) * mpf(10) ** -100)


class LogTest(unittest.TestCase):
    def test_log_matches_logm(self):
        for n in (1, 2):  # `mp.logm` is complex near pi, and wrong for the (3+N)-square matrix
            for tau in [*taus(n)[:2], [1.2, -0.7, 1.4] + taus(n)[0][3:]]:
                E = sen3.expm([mpf(c) for c in tau], n)
                q = [precision.to_f64(c) for c in so3.exp_series([mpf(c) for c in tau[:3]])]
                x = [precision.to_f64(E[r][3 + k]) for k in range(n) for r in range(3)]
                got = sen3.log(n)({"q": q, "x": x})["tau"]
                R = so3.rot(so3.normalize([mpf(c) for c in q]))
                X = sen3.eye(3 + n)
                for r in range(3):
                    X[r][:3] = R[r]
                    for k in range(n):
                        X[r][3 + k] = mpf(x[3 * k + r])
                L = mp.logm(mp.matrix(X))
                want = sen3.vee([[L[i, j] for j in range(3 + n)] for i in range(3 + n)], n)
                self.assertTrue(close(mp.matrix(got), mp.matrix(want), mpf(10) ** -90), (n, tau))

    def test_a_record_and_its_negative_have_one_log(self):
        q = [precision.to_f64(c) for c in so3.exp_series([mpf(0.3), mpf(-0.5), mpf(0.7)])]
        a = sen3.log(2)({"q": q, "x": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})["tau"]
        b = sen3.log(2)({"q": [-c for c in q], "x": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})["tau"]
        self.assertEqual(a, b)


class CheckTest(unittest.TestCase):
    """A cross-check passes a correct record and fails an error of 1e-95, not of 1e-105."""

    def cases(self, n):
        tau = taus(n)[0]
        q = [precision.to_f64(c) for c in so3.exp_series([mpf(c) for c in tau[:3]])]
        for name in IDS:
            yield name, ({"q": q, "x": tau[3:]} if name in ("log", "ad") else {"tau": tau})

    def test_every_id_has_its_own_check_and_it_detects_a_planted_error(self):
        for n in (1, 2):
            m = 3 + 3 * n
            for name, inputs in self.cases(n):
                spec, key = (
                    FUNCTIONS[f"sen3_{name}_n{n}"],
                    {"exp": "x", "log": "tau", "ad": "Ad"}.get(name, "J"),
                )
                good = spec.evaluate(inputs)
                spec.check(inputs, good)
                mat = isinstance(good[key], Mat)
                values = list(good[key].data if mat else good[key])
                for i in range(0, len(values), 11):
                    (bi, bj), zero = divmod(i, m)[::-1], False
                    if mat:  # a structural zero is exactly zero: any error fails
                        zero = sen3.structural_zero(bi // 3, bj // 3)
                    for size, ok in ((mpf(10) ** -95, False), (mpf(10) ** -105, True)):
                        bad = list(values)
                        bad[i] = values[i] + size
                        out = {**good, key: Mat(good[key].shape, bad) if mat else bad}
                        if ok and not zero:
                            spec.check(inputs, out)
                        else:
                            with self.assertRaises(check.CrossCheckError, msg=(name, n, i)):
                                spec.check(inputs, out)

    def test_translation_first_exp_is_refused(self):
        tau = taus(1)[0]
        out = sen3.exp(1)({"tau": tau[3:] + tau[:3]})
        with self.assertRaises(check.CrossCheckError):
            FUNCTIONS["sen3_exp_n1"].check({"tau": tau}, out)

    def test_the_wrong_sign_of_the_q_block_is_refused(self):
        tau = taus(2)[1]
        for name, other in (("jr", "jl"), ("jl", "jr")):
            wrong = sen3.jacobian(2, other)({"tau": tau})
            with self.assertRaises(check.CrossCheckError, msg=name):
                FUNCTIONS[f"sen3_{name}_n2"].check({"tau": tau}, wrong)

    def test_the_dual_structure_is_asserted_on_the_data(self):
        n, tau = 2, taus(2)[0]
        J = sen3.series([mpf(c) for c in tau], n, 1)
        check_sen3._dual("ok", J, n)
        for bi, bj, size in (
            (0, 1, 1e-200),
            (1, 2, 1e-200),
        ):  # a block above, and off the first column
            bad = [list(r) for r in J]
            bad[3 * bi][3 * bj] = mpf(size)
            with self.assertRaises(check.CrossCheckError):
                check_sen3._dual("bad", bad, n)
        bad = [list(r) for r in J]
        bad[3][3] += mpf(10) ** -95  # diagonal blocks are one block
        with self.assertRaises(check.CrossCheckError):
            check_sen3._dual("bad", bad, n)


class InputsTest(unittest.TestCase):
    def test_strata(self):
        names = [s.name for s in SEN3_STRATA]
        self.assertEqual((len(names), len(set(names))), (52, 52))
        self.assertNotIn("theta:dense", names)
        self.assertEqual(names[0], "theta:1e-12")
        self.assertEqual(names[26], "theta:pi-1e-12")
        self.assertEqual(
            names[27:32], [f"rho:1e-6/theta={t}" for t in ("1e-8", "1e-4", "1e-1", "1", "pi-1e-6")]
        )
        self.assertEqual(names[-1], "rho:1e4/theta=pi-1e-6")

    def test_every_stratum_has_the_same_number_of_records(self):
        for s in SEN3_STRATA:
            self.assertEqual(len(sen3.tau_inputs(2)(s)), SEN3_SAMPLES, s.name)
            self.assertEqual(len(sen3.x_inputs(2, False)(s)), SEN3_SAMPLES, s.name)
            self.assertEqual(len(sen3.x_inputs(2, True)(s)), 2 * SEN3_SAMPLES, s.name)

    def test_the_quaternion_ids_also_see_w0_and_nonunit(self):
        self.assertEqual(SEN3_QUAT_STRATA[:-2], SEN3_STRATA)
        w0, nonunit = SEN3_QUAT_STRATA[-2:]
        self.assertEqual((w0.name, nonunit.name), ("q:w0", "q:nonunit"))
        for name, spec in FUNCTIONS.items():
            if name.startswith("sen3_"):
                takes_q = name.startswith(("sen3_log_", "sen3_ad_"))
                self.assertIs(spec.strata is SEN3_QUAT_STRATA, takes_q, name)
        for n in (1, 2, 3):
            recs = sen3.x_inputs(n, True)(w0)
            self.assertEqual(len(recs), SEN3_SAMPLES)  # w = +0 has no negative
            self.assertTrue(all(r["q"][0] == 0 and math.copysign(1, r["q"][0]) > 0 for r in recs))
            self.assertEqual(len(sen3.x_inputs(n, True)(nonunit)), 2 * SEN3_SAMPLES)
            self.assertEqual(len(sen3.x_inputs(n, False)(w0)), SEN3_SAMPLES)
        for rec in sen3.x_inputs(2, True)(w0):  # Log at an angle of exactly pi
            phi = sen3.log(2)(rec)["tau"][:3]
            self.assertLess(abs(so3.norm(phi) - mp.pi), mpf(10) ** -100)
        dev = [so3.norm([mpf(c) for c in r["q"]]) ** 2 - 1 for r in sen3.x_inputs(1, False)(nonunit)]
        self.assertTrue(all(abs(abs(d) - mpf(2) ** -45) < 1e-15 for d in dev))
        self.assertEqual([d > 0 for d in dev], [i % 2 == 0 for i in range(SEN3_SAMPLES)])

    def test_the_translations_of_a_record_are_distinct(self):
        """A block i mixed up with a block j goes unseen if rho_i = rho_j."""
        for s in SEN3_QUAT_STRATA:
            vectors = [r["x"] for r in sen3.x_inputs(3, False)(s)]
            if s in SEN3_STRATA:
                vectors += [r["tau"][3:] for r in sen3.tau_inputs(3)(s)]
            for v in vectors:
                self.assertEqual(len({tuple(r) for r in sen3.chunks(v, 3)}), 3, s.name)

    def test_translations_are_random_directions_at_the_scale_and_nest_in_n(self):
        by = {s.name: s for s in SEN3_STRATA}
        for e in (-6, -3, 0, 3, 4):
            s = by[f"rho:1e{e}/theta=1"]
            two = sen3.tau_inputs(2)(s)
            for rec in two:
                for r in sen3.chunks([mpf(c) for c in rec["tau"][3:]], 2):
                    self.assertLess(abs(so3.norm(r) / mpf(10) ** e - 1), mpf(10) ** -15, e)
            self.assertEqual([r["tau"][:6] for r in two], [r["tau"] for r in sen3.tau_inputs(1)(s)])
        unit = sen3.tau_inputs(1)(by["theta:1e-5"])
        self.assertTrue(
            all(abs(so3.norm([mpf(c) for c in r["tau"][3:]]) - 1) < 1e-15 for r in unit)
        )

    def test_thetas(self):
        by = {s.name: s for s in SEN3_STRATA}
        want = {
            "1e-8": 1e-8,
            "1e-4": 1e-4,
            "1e-1": 0.1,
            "1": 1.0,
            "pi-1e-6": 3.141592653589793 - 1e-6,
        }
        for name, theta in want.items():
            (rec, *_) = sen3.tau_inputs(1)(by[f"rho:1e0/theta={name}"])
            self.assertLess(
                abs(so3.norm([mpf(c) for c in rec["tau"][:3]]) - mpf(theta)), theta * 1e-15
            )
        (rec, *_) = sen3.tau_inputs(1)(by["theta:exact0"])
        self.assertEqual(rec["tau"][:3], [0.0, 0.0, 0.0])
        self.assertEqual(
            sen3.x_inputs(1, True)(by["theta:exact0"])[1]["q"], [-1.0, -0.0, -0.0, -0.0]
        )


if __name__ == "__main__":
    unittest.main()
