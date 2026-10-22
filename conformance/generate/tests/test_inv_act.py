"""The inverse-action ids (docs/decisions/0103): the references against differentiation through
Exp, the checks against planted errors, and the strata's draws."""

import unittest

from mpmath import mp, mpf

from gen import check, check_inv_act, inv_act, precision, sen3, so3
from gen.fmt import Mat
from gen.precision import is_binary32
from gen.strata import drawn_inputs

Q = [0.5, -0.1, 0.7, 0.4]  # not unit: the reference normalizes
T, P = [3.0, -1.5, 0.25], [-0.5, 2.0, 1.25]


def _homogeneous(q, t):
    X = sen3.eye(4)
    R = so3.rot(so3.normalize([mpf(c) for c in q]))
    for i in range(3):
        X[i][:3] = R[i]
        X[i][3] = mpf(t[i])
    return mp.matrix(X)


def _by_diff(q, t, p, right: bool, dof: int) -> list:
    """Column j of d/ds (X Exp(s e_j))^-1 p (right) or (Exp(s e_j) X)^-1 p (left) at 0, by
    `mp.diff` through `mp.expm` of the hat: the definition, with no generator taken by hand."""
    X, pt = _homogeneous(q, t), mp.matrix([mpf(c) for c in p] + [1])

    def moved(j: int, i: int):
        def f(s):
            E = mp.expm(mp.matrix(sen3.hat([s if k == j else 0 for k in range(6)], 1)))
            return (mp.inverse(X * E if right else E * X) * pt)[i]

        return f

    return [[mp.diff(moved(j, i), 0) for i in range(3)] for j in range(dof)]


class InvActTest(unittest.TestCase):
    def setUp(self):
        precision.setup()

    def test_the_jacobians_are_the_derivatives_through_exp(self):
        with mp.workdps(60):
            out = inv_act.se3_evaluate({"q": Q, "x": T, "p": P})
            for key, right in (("Jr", True), ("Jl", False)):
                want = _by_diff(Q, T, P, right, 6)
                got = out[key].data
                err = max(abs(g - w) for g, w in zip(got, (c for col in want for c in col)))
                self.assertLess(err, mpf(10) ** -40, key)
            so = inv_act.so3_evaluate({"q": Q, "p": P})
            for key, right in (("Jr", True), ("Jl", False)):
                want = _by_diff(Q, [0, 0, 0], P, right, 3)
                err = max(abs(g - w) for g, w in zip(so[key].data, (c for col in want for c in col)))
                self.assertLess(err, mpf(10) ** -40, key)

    def test_every_output_is_checked_and_a_planted_error_is_caught(self):
        cases = (
            ("so3_inv_act", {"q": Q, "p": P}, inv_act.so3_evaluate, check_inv_act.so3_inv_act),
            ("se3_inv_act", {"q": Q, "x": T, "p": P}, inv_act.se3_evaluate, check_inv_act.se3_inv_act),
        )
        for name, rec, evaluate, check_ in cases:
            out = evaluate(rec)
            check_(rec, out)
            for key, val in out.items():
                for i in range(len(val.data) if isinstance(val, Mat) else len(val)):
                    bad = dict(out)
                    if isinstance(val, Mat):
                        data = list(val.data)
                        data[i] += mpf(10) ** -90
                        bad[key] = Mat(val.shape, data)
                    else:
                        bad[key] = [v + (mpf(10) ** -90 if k == i else 0) for k, v in enumerate(val)]
                    with self.assertRaises(check.CrossCheckError, msg=(name, key, i)):
                        check_(rec, bad)

    def test_the_strata_draw_their_ratio_and_the_f32_twins_round_the_same_draws(self):
        by_name = {s.name: s for s in inv_act.SE3_STRATA}
        norm = lambda v: mp.sqrt(sum(mpf(c) ** 2 for c in v))
        for r in inv_act.RATIOS:
            # Rounding p and t to binary64 moves x by up to 2 sqrt(3) u r |x|: the ratio `tx:<r>`
            # holds to that, and `tp:<r>`'s, which rounds nothing it divides, to a few u.
            for rec in drawn_inputs(by_name[f"tx:{r}"]):
                x = inv_act.se3_evaluate(rec)["Xinv_p"]
                ratio = norm(rec["x"]) / norm(x) / mpf(r)
                self.assertLess(abs(ratio - 1), 8 * mpf(2) ** -53 * (mpf(r) + 1), (r, rec))
            for rec in drawn_inputs(by_name[f"tp:{r}"]):
                ratio = norm(rec["x"]) / norm(rec["p"]) / mpf(r)
                self.assertLess(abs(ratio - 1), 8 * mpf(2) ** -53, (r, rec))
            for a, b in zip(drawn_inputs(by_name[f"tp:{r}"]), drawn_inputs(by_name[f"tp:{r}@f32"])):
                self.assertTrue(all(is_binary32(c) for v in b.values() for c in v))
                self.assertLess(abs(norm(a["p"]) / norm(b["p"]) - 1), mpf(2) ** -22)


if __name__ == "__main__":
    unittest.main()
