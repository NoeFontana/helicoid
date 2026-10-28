"""`eig_sym_n{4,9,16}` (0111) against what its strata draw: symmetric and exact at binary32, the
planted structure of each stratum, `kappa` where it applies, and the reference's convention."""

import unittest

from mpmath import mp, mpf

from gen import check_linalg, eig_sym, precision
from gen.dense import exact, from_mat, ldlt_pivots, rank
from gen.fmt import Mat
from gen.precision import DPS, is_binary32
from gen.strata import Drawn

precision.setup()


def records(n: int, name: str, count: int = 2) -> list[dict]:
    s = next(s for s in eig_sym.strata(n) if s.name == name)
    return Drawn(s.name, count, s.make, s.f32, s.stream_of).records()


class StrataTest(unittest.TestCase):
    def test_stored_matrices_are_exactly_symmetric_and_binary32_at_f32(self):
        for n in (4, 9):
            for s in eig_sym.strata(n):
                for record in Drawn(s.name, 1, s.make, s.f32, s.stream_of).records():
                    A = exact(record["A"])
                    self.assertTrue(all(A[i][j] == A[j][i] for i in range(n) for j in range(n)))
                    if s.f32:
                        self.assertTrue(all(is_binary32(x) for x in record["A"].data), s.name)

    def test_the_catalogue_is_the_records(self):
        names = [s.name for s in eig_sym.strata(4)]
        self.assertEqual(len([s for s in names if not s.endswith("@f32")]), 24)
        self.assertEqual(len([s for s in names if s.endswith("@f32")]), 19)

    def test_the_rank_deficient_strata_are_exactly_singular(self):
        for n, r in ((4, 1), (9, 3), (16, 3)):
            for record in records(n, f"eig:rank-def-{r}"):
                self.assertEqual(rank(exact(record["A"])), n - r)

    def test_the_unit_diagonal_strata_have_a_unit_diagonal_and_their_condition_number(self):
        for name, k in (("eig:unit-diag-cond-1e-6", 6), ("eig:unit-diag-cond-1e-4@f32", 4)):
            for record in records(9, name):
                A = record["A"]
                self.assertTrue(all(A.data[10 * i] == 1.0 for i in range(9)), name)
                with mp.workdps(DPS):
                    kappa = eig_sym.kappa(A)
                self.assertLess(abs(mp.log10(kappa) - k), 0.05, name)

    def test_the_graded_strata_are_positive_definite_and_graded(self):
        for record in records(16, "eig:graded-1e-16"):
            A = record["A"]
            self.assertIsNotNone(ldlt_pivots(exact(A)))
            self.assertLess(abs(mp.log10(A.data[-1] / A.data[0]) + 16), 0.5)

    def test_the_diagonal_stratum_has_its_two_ties(self):
        for record in records(9, "eig:diagonal"):
            d = [record["A"].data[10 * i] for i in range(9)]
            self.assertEqual((d[2], d[8]), (d[0], d[1]))


class ReferenceTest(unittest.TestCase):
    def test_a_reference_passes_its_cross_check_and_signs_its_columns(self):
        for name in ("eig:random", "eig:rank-def-3", "eig:graded-1e-8", "eig:diagonal@f32"):
            (record,) = records(9, name, 1)
            with mp.workdps(DPS):
                out = eig_sym.eig_sym(record)
                check_linalg.eig_sym(record, out)
            V = out["V"]
            for j in range(9):
                col = V.data[9 * j : 9 * j + 9]
                top = max(range(9), key=lambda r: (abs(col[r]), -r))
                self.assertGreater(col[top], 0, name)

    def test_kappa_is_zero_off_the_positive_definite_cone(self):
        (record,) = records(4, "eig:rank-def-1", 1)
        with mp.workdps(DPS):
            self.assertEqual(eig_sym.eig_sym(record)["kappa"], [0])

    def test_a_wrong_kappa_fails_the_cross_check(self):
        (record,) = records(4, "eig:graded-1e-4", 1)
        with mp.workdps(DPS):
            out = eig_sym.eig_sym(record)
            out["kappa"] = [out["kappa"][0] * (1 + mp.mpf(10) ** -50)]
            with self.assertRaises(check_linalg.CrossCheckError):
                check_linalg.eig_sym(record, out)

    def test_a_null_space_takes_one_basis_at_every_precision(self):
        (record,) = records(4, "eig:rank-def-3", 1)
        with mp.workdps(DPS):
            low = eig_sym.eig_sym(record)["V"]
        with mp.workdps(DPS + 30):
            high = eig_sym.eig_sym(record)["V"]
        self.assertTrue(all(abs(a - b) < mpf(10) ** -100 for a, b in zip(low.data, high.data)))

    def test_an_exactly_singular_matrix_has_exact_zeros(self):
        A = Mat((4, 4), [1.0, 1.0, 0.0, 0.0] * 2 + [0.0] * 8)
        with mp.workdps(DPS):
            lam = eig_sym.eig_sym({"A": A})["lambda"]
        self.assertEqual(lam[:3], [0, 0, 0])
        self.assertEqual(from_mat(A).rows, 4)


if __name__ == "__main__":
    unittest.main()
