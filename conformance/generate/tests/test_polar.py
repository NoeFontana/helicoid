"""`nearest_rotation` against what its strata draw, Newton's polar iteration, the symmetric
H = R^T M of SO.13(c), and its cross-check (docs/decisions/0102 decision 8)."""

import json
import unittest

from mpmath import mp, mpf

from gen import check_linalg, corpus, polar, precision
from gen.dense import det3, exact, from_mat, rank
from gen.fmt import Mat
from gen.precision import is_binary32
from gen.registry import FUNCTIONS

precision.setup()
BY_NAME = {s.name: s for s in polar.STRATA}
TAU = {False: mpf(2) ** -26, True: mpf(2) ** -12}  # 0102 decision 2


def mat(rows) -> Mat:
    """A record's `M` from its rows."""
    return Mat((3, 3), [float(rows[i][j]) for j in range(3) for i in range(3)])


def det(m: Mat):
    return det3(exact(m))


def sigma(m: Mat) -> list:
    return sorted(mp.svd_r(from_mat(m), compute_uv=False), reverse=True)


def rho(m: Mat):
    """det M / |M|_F^3, the quantity 0102's arm test holds against tau."""
    return mpf(det(m)) / mp.sqrt(sum(mpf(x) ** 2 for x in m.data)) ** 3


class StrataTest(unittest.TestCase):
    def test_the_catalogue_is_0102s_table(self):
        names = [s.name for s in polar.STRATA]
        half = names.index("polar:near-1e-2@f32")
        self.assertTrue(all(n.endswith("@f32") == (i >= half) for i, n in enumerate(names)))
        self.assertEqual(sum(s.count for s in polar.STRATA[:half]), 296)
        self.assertEqual(sum(s.count for s in polar.STRATA[half:]), 264)
        self.assertEqual(names[5:8], ["polar:exact", "polar:scaled", "polar:random"])
        self.assertEqual(names[half - 2 : half], ["polar:scale-up", "polar:scale-down"])
        self.assertEqual(BY_NAME["polar:exact"].count, 24)
        self.assertEqual(BY_NAME["polar:random@f32"].count, 32)
        self.assertEqual(BY_NAME["polar:cond-1e-8@f32"].stream_of, "polar:cond-1e-8")
        self.assertIsNone(BY_NAME["polar:cond-1e-2@f32"].stream_of)  # no binary64 namesake
        self.assertIsNone(BY_NAME["polar:near-1e-6@f32"].stream_of)
        self.assertIsNone(BY_NAME["polar:scale-up@f32"].stream_of)  # its own scale

    def test_inputs_are_binary32_at_f32(self):
        for s in polar.STRATA:
            if s.f32:
                for record in s.records()[:4]:
                    self.assertTrue(all(is_binary32(x) for x in record["M"].data), s.name)

    def test_near_strata_are_within_their_size_of_so3(self):
        cases = (("polar:near-1e-2", 2), ("polar:near-1e-12", 12), ("polar:near-1e-6@f32", 6))
        for name, k in cases:
            for record in BY_NAME[name].records()[:4]:
                s = sigma(record["M"])
                self.assertLess(max(abs(x - 1) for x in s), 20 * mpf(10) ** -k, name)
                self.assertGreater(det(record["M"]), 0, name)

    def test_exact_is_the_24_signed_permutation_rotations(self):
        records = BY_NAME["polar:exact@f32"].records()
        self.assertEqual(len({tuple(r["M"].data) for r in records}), 24)
        for r in records:
            self.assertEqual(det(r["M"]), 1)
            self.assertEqual(polar.nearest_rotation(r)["R"].data, r["M"].data)

    def test_scaled_is_a_scaled_rotation(self):
        for record in BY_NAME["polar:scaled"].records():
            s = sigma(record["M"])
            self.assertLess((s[0] - s[2]) / s[0], 1e-15)
            self.assertTrue(mpf(10) ** -3 <= s[1] <= mpf(10) ** 3, s)

    def test_random_has_both_signs_of_det(self):
        for name in ("polar:random", "polar:random@f32"):
            signs = {det(r["M"]) > 0 for r in BY_NAME[name].records()}
            self.assertEqual(signs, {False, True}, name)

    def test_reflection_is_the_reflect_arm(self):
        """det M <= -tau |M|^3 at the precision and |sigma3| < sigma2 / 2."""
        for name, f32 in (("polar:reflection", False), ("polar:reflection@f32", True)):
            for record in BY_NAME[name].records():
                s = sigma(record["M"])
                self.assertLess(rho(record["M"]), -TAU[f32], name)
                self.assertLess(s[2], s[1] / 2, name)

    def test_cond_holds_both_signs_and_both_sides_of_tau(self):
        for name, k, f32, newton in (
            ("polar:cond-1e-4", 4, False, True),
            ("polar:cond-1e-8", 8, False, False),
            ("polar:cond-1e-2@f32", 2, True, True),
            ("polar:cond-1e-4@f32", 4, True, False),
        ):
            records = BY_NAME[name].records()
            for i, record in enumerate(records):
                self.assertEqual(det(record["M"]) > 0, i % 2 == 0, (name, i))
                self.assertEqual(abs(rho(record["M"])) > TAU[f32], newton, (name, i))
                self.assertLess(abs(sigma(record["M"])[2] / mpf(10) ** -k - 1), 1e-3, name)

    def test_rank2_is_exactly_rank_2(self):
        for name, tol in (("polar:rank2", 1e-7), ("polar:rank2@f32", 1e-2)):
            for record in BY_NAME[name].records():
                self.assertEqual(rank(exact(record["M"])), 2, name)
                s = sigma(record["M"])
                self.assertLess(abs(s[0] - 1) + abs(s[1] - 0.5), tol, (name, s))

    def test_non_unique_is_tied_to_rounding_and_not_exactly(self):
        """sigma2 - |sigma3| is the rounding of M, but never 0: the nearest rotation is unique."""
        for name, tol in (("polar:non-unique", 1e-14), ("polar:non-unique@f32", 1e-5)):
            for record in BY_NAME[name].records():
                s = sigma(record["M"])
                self.assertLess(det(record["M"]), 0, name)
                self.assertLess(s[1] - s[2], tol, (name, s))
                self.assertTrue(polar.unique(exact(record["M"])), name)

    def test_scale_strata_are_polar_random_scaled(self):
        for name, f32, s in (("polar:scale-up", False, 70), ("polar:scale-down@f32", True, -12)):
            base = BY_NAME["polar:random@f32" if f32 else "polar:random"].records()[0]["M"].data
            scaled = BY_NAME[name].records()[0]["M"].data
            for a, b in zip(base, scaled, strict=True):
                self.assertLess(abs(mpf(b) / mpf(10) ** s / mpf(a) - 1), 2e-7 if f32 else 3e-16)


def newton(m: Mat, steps: int = 60):
    """Newton's polar iteration X <- (X + X^-T) / 2 (PO.1, unscaled), an independent route to the
    nearest rotation of a det > 0 matrix."""
    X = from_mat(m)
    for _ in range(steps):
        X = (X + mp.inverse(X).T) / 2
    return X


class ReferenceTest(unittest.TestCase):
    def test_newton_agrees_where_det_is_positive(self):
        for name in ("polar:near-1e-2", "polar:random", "polar:scaled", "polar:cond-1e-4"):
            for record in BY_NAME[name].records()[:8]:
                if det(record["M"]) <= 0:
                    continue
                R = from_mat(polar.nearest_rotation(record)["R"])
                diff = newton(record["M"]) - R
                self.assertLess(max(abs(x) for x in diff), mpf(10) ** -100, name)

    def test_r_t_m_is_symmetric_with_a_positive_smallest_pair(self):
        """SO.13(c): R is the unique nearest rotation when H = R^T M is symmetric and its two
        smallest eigenvalues sum to more than 0; `polar:non-unique`'s sum is the rounding's."""
        for s in polar.STRATA:
            for record in s.records()[:2]:
                M, R = from_mat(record["M"]), from_mat(polar.nearest_rotation(record)["R"])
                H = R.T * M
                norm = max(abs(x) for x in M)
                self.assertLess(max(abs(x) for x in H - H.T), mpf(10) ** -100 * norm, s.name)
                lam = sorted(mp.eigsy((H + H.T) / 2, eigvals_only=True))
                self.assertGreater(lam[0] + lam[1], 0, s.name)

    def test_a_reflection_is_turned_by_its_smallest_singular_vector(self):
        """diag(1, 2, -5): sigma = (5, 2, 1), d = -1, so R = diag(-1, 1, -1), not I."""
        out = polar.nearest_rotation({"M": mat([[1, 0, 0], [0, 2, 0], [0, 0, -5]])})
        want = [-1, 0, 0, 0, 1, 0, 0, 0, -1]
        self.assertLess(max(abs(a - b) for a, b in zip(out["R"].data, want, strict=True)), 1e-110)

    def test_an_exactly_singular_m_has_one_nearest_rotation(self):
        record = {"M": mat([[1, 2, 3], [4, 5, 6], [7, 8, 9]])}
        out = polar.nearest_rotation(record)
        check_linalg.nearest_rotation(record, out)
        with mp.workdps(precision.RECHECK_DPS):
            hi = polar.nearest_rotation(record)
        for a, b in zip(out["R"].data, hi["R"].data, strict=True):
            self.assertLess(abs(a - b), mpf(10) ** -110)

    def test_an_exact_tie_is_refused(self):
        """diag(1, 1, -5) (sigma = (5, 1, -1)), a rank-1 M and 0: a circle of nearest rotations."""
        for rows in (
            [[1, 0, 0], [0, 1, 0], [0, 0, -5]],
            [[2, 4, -6], [1, 2, -3], [-3, -6, 9]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        ):
            self.assertFalse(polar.unique(exact(mat(rows))), rows)
            with self.assertRaises(ValueError):
                polar.nearest_rotation({"M": mat(rows)})
        for rows in ([[1, 0, 0], [0, 2, 0], [0, 0, -5]], [[2, 0, 0], [0, 2, 0], [0, 0, -1]]):
            self.assertTrue(polar.unique(exact(mat(rows))), rows)  # sigma3 simple, or det > 0

    def test_every_stratum_passes_the_cross_check(self):
        for s in polar.STRATA:
            for record in s.records()[:2]:
                check_linalg.nearest_rotation(record, polar.nearest_rotation(record))


class CrossCheckTest(unittest.TestCase):
    def test_a_perturbation_beyond_100_digits_fails(self):
        record = BY_NAME["polar:random"].records()[0]
        good = polar.nearest_rotation(record)
        check_linalg.nearest_rotation(record, good)
        for rel, fails in ((mpf(10) ** -95, True), (mpf(10) ** -105, False)):
            for i in (0, 4, 8):
                data = list(good["R"].data)
                data[i] *= 1 + rel
                bad = {"R": Mat((3, 3), data)}
                if fails:
                    with self.assertRaises(check_linalg.CrossCheckError, msg=i):
                        check_linalg.nearest_rotation(record, bad)
                else:
                    check_linalg.nearest_rotation(record, bad)

    def test_a_rotation_that_is_not_the_nearest_fails(self):
        """R turned by 10^-3 about e3, and R diag(1, -1, -1): rotations, but not M's nearest."""
        for name in ("polar:random", "polar:near-1e-8", "polar:reflection"):
            record = BY_NAME[name].records()[3]
            R = from_mat(polar.nearest_rotation(record)["R"])
            c, s = mp.cos(mpf(10) ** -3), mp.sin(mpf(10) ** -3)
            for G in (mp.matrix([[c, -s, 0], [s, c, 0], [0, 0, 1]]), mp.diag([1, -1, -1])):
                Q = R * G
                bad = {"R": Mat((3, 3), [Q[i, j] for j in range(3) for i in range(3)])}
                with self.assertRaises(check_linalg.CrossCheckError, msg=name):
                    check_linalg.nearest_rotation(record, bad)


class EncodingTest(unittest.TestCase):
    def test_one_matrix_per_part_keeps_its_shape(self):
        data, n, _ = corpus.build(FUNCTIONS["nearest_rotation"], only={"polar:exact@f32"})
        records = [json.loads(line) for line in data.decode().splitlines()]
        self.assertEqual(n, 24)
        for r in records:
            self.assertEqual(sorted(r["in"]), ["M", "shape"])
            self.assertEqual(sorted(r["out"]), ["R", "shape"])
            self.assertEqual(r["in"]["shape"], [3, 3])
            self.assertEqual(r["out"]["shape"], [3, 3])
        self.assertEqual(records[0]["out"]["R"][1], "0.00000000000000000000000000000e0")


if __name__ == "__main__":
    unittest.main()
