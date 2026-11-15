"""`rq3` against what its strata draw, RQ.1's construction and its cross-check
(docs/decisions/0117 decision 9)."""

import unittest

from mpmath import mp, mpf

from gen import check_linalg, precision, rq
from gen.dense import det3, exact, from_mat, rank
from gen.fmt import Mat
from gen.precision import is_binary32

precision.setup()
BY_NAME = {s.name: s for s in rq.STRATA}


def mat(rows) -> Mat:
    """A record's `M` from its rows."""
    return Mat((3, 3), [float(rows[i][j]) for j in range(3) for i in range(3)])


def det(m: Mat):
    return det3(exact(m))


class StrataTest(unittest.TestCase):
    def test_the_catalogue_is_0117s_table(self):
        names = [s.name for s in rq.STRATA]
        half = names.index("rq:random@f32")
        self.assertTrue(all(n.endswith("@f32") == (i >= half) for i, n in enumerate(names)))
        self.assertEqual(sum(s.count for s in rq.STRATA[:half]), 388)
        self.assertEqual(sum(s.count for s in rq.STRATA[half:]), 372)
        self.assertEqual(names[:3], ["rq:random", "rq:diagonal", "rq:triangular"])
        self.assertEqual(names[half - 2 : half], ["rq:scale-up", "rq:scale-down"])
        self.assertEqual(BY_NAME["rq:zero"].count, 4)
        self.assertEqual(BY_NAME["rq:random@f32"].stream_of, "rq:random")
        self.assertEqual(BY_NAME["rq:parallel-1e-4@f32"].stream_of, "rq:parallel-1e-4")
        self.assertIsNone(BY_NAME["rq:parallel-1e-7@f32"].stream_of)  # no binary64 namesake
        self.assertIsNone(BY_NAME["rq:diagonal@f32"].stream_of)  # its own range
        self.assertIsNone(BY_NAME["rq:row-spread@f32"].stream_of)

    def test_inputs_are_binary32_at_f32(self):
        for s in rq.STRATA:
            if s.f32:
                for record in s.records()[:4]:
                    self.assertTrue(all(is_binary32(x) for x in record["M"].data), s.name)

    def test_the_exact_strata_are_positive_triangles(self):
        for name in ("rq:diagonal", "rq:triangular", "rq:diagonal@f32", "rq:triangular@f32"):
            for record in BY_NAME[name].records():
                E = exact(record["M"])
                self.assertTrue(all(E[i][j] == 0 for i in range(3) for j in range(i)), name)
                self.assertTrue(all(E[i][i] > 0 for i in range(3)), name)
                if "diagonal" in name:
                    self.assertTrue(all(E[i][j] == 0 for i in range(3) for j in range(3) if i != j))

    def test_parallel_rows_are_at_their_angle(self):
        for name, k in (("rq:parallel-1e-4", 4), ("rq:parallel-1e-12", 12), ("rq:parallel-1e-6@f32", 6)):
            for record in BY_NAME[name].records()[:4]:
                M = from_mat(record["M"])
                a, b = [M[2, j] for j in range(3)], [M[1, j] for j in range(3)]
                c = rq._cross(a, b)
                sin = mp.sqrt(rq._dot(c, c) / (rq._dot(a, a) * rq._dot(b, b)))
                self.assertLess(abs(sin / mpf(10) ** -k - 1), 1e-2 if "f32" not in name else 0.5, name)

    def test_the_signs_of_det_are_planted(self):
        for record in BY_NAME["rq:negative-det"].records():
            self.assertLess(det(record["M"]), 0)
        for name in ("rq:rank-bottom", "rq:rank-bottom@f32"):
            for record in BY_NAME[name].records():
                self.assertEqual(rank(exact(record["M"])[1:]), 1, name)

    def test_scale_strata_are_rq_random_scaled(self):
        for name, f32, s in (("rq:scale-up", False, 70), ("rq:scale-down@f32", True, -12)):
            base = BY_NAME["rq:random@f32" if f32 else "rq:random"].records()[0]["M"].data
            scaled = BY_NAME[name].records()[0]["M"].data
            for a, b in zip(base, scaled, strict=True):
                self.assertLess(abs(mpf(b) / mpf(10) ** s / mpf(a) - 1), 2e-7 if f32 else 3e-16)


class ReferenceTest(unittest.TestCase):
    def test_a_triangle_is_its_own_factor(self):
        out = rq.rq3({"M": mat([[2, 1, -1], [0, 3, 0.5], [0, 0, 0.25]])})
        self.assertEqual(out["U"].data, [2, 0, 0, 1, 3, 0, -1, 0.5, 0.25])
        self.assertEqual(out["R"].data, [1, 0, 0, 0, 1, 0, 0, 0, 1])
        self.assertEqual(out["valid"], [1])

    def test_a_reflection_keeps_the_factorization(self):
        """diag(-1, 2, 3): R in SO(3) with u22, u33 > 0 forces R = I and u11 = det M / 6 = -1, the
        factorization kept under a clear mask (0117 decision 1)."""
        out = rq.rq3({"M": mat([[-1, 0, 0], [0, 2, 0], [0, 0, 3]])})
        self.assertEqual(out["U"].data, [-1, 0, 0, 0, 2, 0, 0, 0, 3])
        self.assertEqual(out["valid"], [0])

    def test_the_tangent_is_the_derivative_of_the_construction(self):
        record = BY_NAME["rq:random"].records()[0]
        check_linalg.rq3(record, rq.rq3(record))

    def test_dependent_rows_are_completed_and_unscored(self):
        for rows in ([[1, 2, 3], [1, 1, 2], [2, 2, 4]], [[0, 0, 0]] * 3, [[1, 0, 0], [0, 1, 0], [0, 0, 0]]):
            out = rq.rq3({"M": mat(rows)})
            self.assertEqual(out["kappa_U"], [0])
            check_linalg.rq3({"M": mat(rows)}, out)
        zero = rq.rq3({"M": mat([[0, 0, 0]] * 3)})
        self.assertEqual(zero["R"].data, [1, 0, 0, 0, 1, 0, 0, 0, 1])

    def test_every_stratum_passes_the_cross_check(self):
        for s in rq.STRATA:
            for record in s.records()[:2]:
                check_linalg.rq3(record, rq.rq3(record))


class CrossCheckTest(unittest.TestCase):
    def test_a_perturbation_beyond_100_digits_fails(self):
        record = BY_NAME["rq:random"].records()[0]
        good = rq.rq3(record)
        check_linalg.rq3(record, good)
        bad = dict(good)
        bad["R"] = Mat((3, 3), [x * (1 + mpf(10) ** -90) for x in good["R"].data])
        with self.assertRaises(check_linalg.CrossCheckError):
            check_linalg.rq3(record, bad)

    def test_a_wrong_mask_or_kappa_fails(self):
        record = BY_NAME["rq:random"].records()[0]
        good = rq.rq3(record)
        with self.assertRaises(check_linalg.CrossCheckError):
            check_linalg.rq3(record, {**good, "valid": [1 - good["valid"][0]]})
        with self.assertRaises(check_linalg.CrossCheckError):
            check_linalg.rq3(record, {**good, "kappa_U": [good["kappa_U"][0] / 4]})


if __name__ == "__main__":
    unittest.main()
