"""`svd3` against what its strata draw, the convention of docs/decisions/0101 decision 1, the
eigenvalues of A^T A, and its keyed `shape`."""

import json
import math
import unittest

from mpmath import mp, mpf

from gen import check_linalg, corpus, precision, svd
from gen.dense import exact, from_mat
from gen.fmt import Mat
from gen.precision import is_binary32
from gen.registry import FUNCTIONS

precision.setup()
BY_NAME = {s.name: s for s in svd.STRATA}


def mat(rows) -> Mat:
    """A record's `A` from its rows."""
    return Mat((3, 3), [float(rows[i][j]) for j in range(3) for i in range(3)])


def det(m: Mat):
    A = exact(m)
    return (
        A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
        - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
        + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])
    )


class StrataTest(unittest.TestCase):
    def test_the_catalogue_is_0101s_table(self):
        names = [s.name for s in svd.STRATA]
        half = names.index("svd:random@f32")
        self.assertTrue(all(n.endswith("@f32") == (i >= half) for i, n in enumerate(names)))
        self.assertEqual(sum(s.count for s in svd.STRATA[:half]), 580)
        self.assertEqual(sum(s.count for s in svd.STRATA[half:]), 468)
        self.assertEqual(names[:2], ["svd:random", "svd:gap-1e-1/bottom"])
        self.assertEqual(names[half - 2 : half], ["svd:scale-up", "svd:scale-down"])
        self.assertEqual(BY_NAME["svd:zero"].count, 4)
        self.assertIn("svd:graded-1e-6@f32", BY_NAME)
        self.assertIsNone(BY_NAME["svd:graded-1e-6@f32"].stream_of)  # no binary64 namesake
        self.assertEqual(BY_NAME["svd:graded-1e-4@f32"].stream_of, "svd:graded-1e-4")
        self.assertIsNone(BY_NAME["svd:col-graded@f32"].stream_of)  # its own grades

    def test_inputs_are_binary32_at_f32(self):
        for s in svd.STRATA:
            if s.f32:
                for record in s.records()[:4]:
                    self.assertTrue(all(is_binary32(x) for x in record["A"].data), s.name)

    def test_the_gap_strata_hold_their_relative_gap(self):
        """The close pair is 10^-k |sigma|_2 apart, to the precision's rounding of A."""
        for name, k, tol in (
            ("svd:gap-1e-12/bottom", 12, 1e-2),
            ("svd:gap-1e-1/top", 1, 1e-12),
            ("svd:gap-1e-6/top@f32", 6, 0.3),
        ):
            for record in BY_NAME[name].records()[:8]:
                s1, s2, s3 = svd.svd3(record)["sigma"]
                close = s2 - abs(s3) if name.split("/")[1].startswith("bottom") else s1 - s2
                rel = close / mp.sqrt(s1 * s1 + s2 * s2 + s3 * s3) / mpf(10) ** -k
                self.assertLess(abs(rel - 1), tol, (name, rel))

    def test_gap_spectrum_is_positive_and_exact(self):
        for top in (False, True):
            for c, x in ((mpf(1) / 2, mpf(1) / 4), (mpf(1), mpf(1))):
                sigma = svd.gap_spectrum(c, x, mpf(10) ** -1, top)
                self.assertTrue(all(s > 0 for s in sigma), sigma)
                pair = sigma[1:] if top else sigma[:2]
                norm = mp.sqrt(sum(s * s for s in sigma))
                self.assertLess(abs((pair[1] - pair[0]) / norm - mpf(10) ** -1), mpf(10) ** -110)

    def test_planted_spectra(self):
        for name, want, tol in (
            ("svd:graded-1e-8", [1, 1e-8, 0.5e-8], 1e-15),
            ("svd:rank1", [1, 0, 0], 1e-15),
            ("svd:essential@f32", None, 1e-6),
        ):
            for record in BY_NAME[name].records()[:4]:
                sigma = svd.svd3(record)["sigma"]
                expect = want or [sigma[0], sigma[0], 0]
                for got, w in zip(sigma, expect, strict=True):
                    self.assertLess(abs(abs(got) - w), tol, (name, sigma))

    def test_ordering_is_ph2s_case(self):
        """U diag(a, a, -b) V^T, b > a: sigma = (b, a, -a), det A < 0."""
        for record in BY_NAME["svd:ordering"].records()[:8]:
            s1, s2, s3 = svd.svd3(record)["sigma"]
            self.assertLess(det(record["A"]), 0)
            self.assertGreater(s1, s2 * 1.2)
            self.assertLess(abs(s2 + s3), 1e-15)

    def test_near_reflection_is_a_reflection(self):
        for name in ("svd:near-reflection-1e-4", "svd:near-reflection-1e-2@f32"):
            for record in BY_NAME[name].records()[:4]:
                self.assertLess(det(record["A"]), 0, name)
                self.assertLess(abs(svd.svd3(record)["sigma"][2] + 1), 1e-1, name)

    def test_col_graded_scales_the_columns(self):
        for name, grades in (("svd:col-graded", (0, 5, 10)), ("svd:col-graded@f32", (0, 2, 4))):
            for record in BY_NAME[name].records()[:4]:
                for j, g in enumerate(grades):
                    norm = math.hypot(*record["A"].data[3 * j : 3 * j + 3])
                    self.assertLess(norm, 10.0 ** (1 - g), (name, j))
                    self.assertGreater(norm, 10.0 ** (-2 - g), (name, j))

    def test_zero_alternates_its_signed_zeros(self):
        records = BY_NAME["svd:zero@f32"].records()
        signs = [[math.copysign(1, x) for x in r["A"].data] for r in records]
        self.assertTrue(all(x == 0 for r in records for x in r["A"].data))
        self.assertEqual(signs[0], [1] * 9)
        self.assertEqual(signs[1], [-1] * 9)
        self.assertEqual(signs[2], [(-1) ** i for i in range(9)])
        self.assertEqual(signs[3], [-((-1) ** i) for i in range(9)])

    def test_scale_strata_are_svd_random_scaled(self):
        for name, f32, s in (("svd:scale-up", False, 70), ("svd:scale-down@f32", True, -12)):
            base = BY_NAME["svd:random@f32" if f32 else "svd:random"].records()[0]["A"].data
            scaled = BY_NAME[name].records()[0]["A"].data
            for a, b in zip(base, scaled, strict=True):
                self.assertLess(abs(mpf(b) / mpf(10) ** s / mpf(a) - 1), 2e-7 if f32 else 3e-16)


class ReferenceTest(unittest.TestCase):
    def assert_convention(self, record, out):
        """U, V in SO(3), A = U Sigma V^T, the ordering and sign sigma3 = sign det A."""
        check_linalg.svd3(record, out)

    def test_diag_1_1_minus_5_is_ordered_by_magnitude(self):
        """PH.2's example: sigma = (5, 1, -1), so U V^T is the nearest rotation, not I."""
        record = {"A": mat([[1, 0, 0], [0, 1, 0], [0, 0, -5]])}
        out = svd.svd3(record)
        self.assert_convention(record, out)
        for got, want in zip(out["sigma"], (5, 1, -1), strict=True):
            self.assertLess(abs(got - want), mpf(10) ** -110)
        U, V = from_mat(out["U"]), from_mat(out["V"])
        self.assertLess(abs(abs(U[2, 0]) - 1), mpf(10) ** -110)  # u1 is e3's axis
        self.assertGreater(U[2, 0], 0)  # its largest component positive

    def test_exactly_rank_deficient_matrices_have_exact_zeros(self):
        for rows, zeros in (
            ([[1, 2, 3], [4, 5, 6], [7, 8, 9]], 1),
            ([[2, 4, -6], [1, 2, -3], [-3, -6, 9]], 2),
            ([[0, 0, 0], [0, 0, 0], [0, 0, 0]], 3),
        ):
            record = {"A": mat(rows)}
            out = svd.svd3(record)
            self.assert_convention(record, out)
            self.assertEqual(out["sigma"][3 - zeros :], [0] * zeros, rows)
            self.assertTrue(all(s > 0 for s in out["sigma"][: 3 - zeros]), rows)
        self.assertEqual(out["U"], Mat((3, 3), [mpf(int(i % 4 == 0)) for i in range(9)]))
        self.assertEqual(out["V"], out["U"])

    def test_the_rank1_complement_is_the_same_at_every_precision(self):
        record = {"A": mat([[2, 4, -6], [1, 2, -3], [-3, -6, 9]])}
        lo = svd.svd3(record)
        with mp.workdps(precision.RECHECK_DPS):
            hi = svd.svd3(record)
        for key in ("U", "V"):
            for a, b in zip(lo[key].data, hi[key].data, strict=True):
                self.assertLess(abs(a - b), mpf(10) ** -110, key)

    def test_singular_values_square_to_the_eigenvalues_of_a_t_a(self):
        """An independent route: mp.eigsy of the stored A^T A, to 1e-100 |A|^2."""
        for name in ("svd:random", "svd:gap-1e-8/top", "svd:ordering", "svd:scale-down"):
            for record in BY_NAME[name].records()[:4]:
                A = from_mat(record["A"])
                eig = sorted(mp.eigsy(A.T * A, eigvals_only=True), reverse=True)
                norm2 = sum(mpf(x) ** 2 for x in record["A"].data)
                for s, lam in zip(svd.svd3(record)["sigma"], eig, strict=True):
                    self.assertLess(abs(s * s - lam), mpf(10) ** -100 * norm2, name)

    def test_columns_one_and_two_are_signed_by_their_largest_component(self):
        for record in BY_NAME["svd:random"].records()[:8]:
            U = svd.svd3(record)["U"].data
            for j in range(2):
                self.assertGreater(max(U[3 * j : 3 * j + 3], key=abs), 0)

    def test_every_stratum_holds_the_convention(self):
        for s in svd.STRATA:
            for record in s.records()[:2]:
                self.assert_convention(record, svd.svd3(record))


class CrossCheckTest(unittest.TestCase):
    def test_a_perturbation_beyond_100_digits_fails(self):
        record = BY_NAME["svd:random"].records()[0]
        good = svd.svd3(record)
        check_linalg.svd3(record, good)
        for rel, fails in ((mpf(10) ** -95, True), (mpf(10) ** -105, False)):
            for i in range(3):
                sigma = list(good["sigma"])
                sigma[i] *= 1 + rel
                bads = [{**good, "sigma": sigma}]
                for key in ("U", "V"):
                    data = list(good[key].data)
                    data[4 * i] *= 1 + rel
                    bads.append({**good, key: Mat((3, 3), data)})
                for bad in bads:
                    if fails:
                        with self.assertRaises(check_linalg.CrossCheckError, msg=(i, bad)):
                            check_linalg.svd3(record, bad)
                    else:
                        check_linalg.svd3(record, bad)

    def test_the_convention_is_enforced(self):
        record = BY_NAME["svd:random"].records()[1]
        good = svd.svd3(record)
        s1, s2, s3 = good["sigma"]
        U = good["U"].data
        reflected = Mat((3, 3), [*U[:6], *(-c for c in U[6:])])  # A V = U Sigma, det U = -1
        for bad in (
            {**good, "U": reflected, "sigma": [s1, s2, -s3]},
            {**good, "sigma": [s2, s1, s3]},
        ):
            with self.assertRaises(check_linalg.CrossCheckError, msg=bad):
                check_linalg.svd3(record, bad)


class EncodingTest(unittest.TestCase):
    def test_two_matrices_in_one_part_key_their_shape(self):
        data, n, _ = corpus.build(FUNCTIONS["svd3"], only={"svd:zero"})
        records = [json.loads(line) for line in data.decode().splitlines()]
        self.assertEqual(n, 4)
        for r in records:
            self.assertEqual(r["in"]["shape"], [3, 3])
            self.assertEqual(r["out"]["shape"], {"U": [3, 3], "V": [3, 3]})
            self.assertEqual(sorted(r["out"]), ["U", "V", "shape", "sigma"])
        self.assertEqual(records[1]["in"]["A"], ["-0x0.0p+0"] * 9)


if __name__ == "__main__":
    unittest.main()
