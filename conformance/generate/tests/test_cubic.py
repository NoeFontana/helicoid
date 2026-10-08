"""`solve_cubic` against what its strata plant, Cardano's formula and the discriminant."""

import unittest
from fractions import Fraction

from mpmath import mp, mpf

from gen import check_linalg, cubic, precision
from gen.precision import is_binary32
from gen.strata import Draw

precision.setup()
BY_NAME = {s.name: s for s in cubic.STRATA}
PLANTED = {
    "cubic:distinct": cubic.planted_distinct,
    "cubic:double": cubic.planted_double,
    "cubic:triple": cubic.planted_triple,
}


def planted(name: str):
    """The stratum's planted (a, roots) again, from its own stream."""
    with mp.workdps(precision.DPS):
        rng = Draw(name, False).stream("cubic")
        draw = cubic.planted_one_real if name == "cubic:one-real" else PLANTED[name]
        return [draw(rng) for _ in range(cubic.N_RANDOM)]


class PlantedTest(unittest.TestCase):
    def test_planted_coefficients_are_exact_at_both_precisions(self):
        for name in (*PLANTED, "cubic:one-real"):
            records = BY_NAME[name].records()
            self.assertEqual(records, BY_NAME[name + "@f32"].records(), name)
            for record in records:
                self.assertTrue(all(is_binary32(x) for x in record.values()), (name, record))

    def test_the_reference_is_the_planted_roots_exactly(self):
        for name in PLANTED:
            for inputs, (a, roots) in zip(BY_NAME[name].records(), planted(name), strict=True):
                self.assertEqual([Fraction(inputs[k]) for k in "abcd"], cubic.expand(a, roots))
                out = cubic.solve_cubic(inputs)
                self.assertEqual([cubic.fraction(r) for r in out["re"]], sorted(roots), inputs)
                self.assertEqual(out["im"], [0, 0, 0])

    def test_one_real_is_r_and_the_planted_pair(self):
        for inputs, (a, (r, p, q)) in zip(
            BY_NAME["cubic:one-real"].records(), planted("cubic:one-real"), strict=True
        ):
            self.assertLess(p * p, 4 * q)
            out = cubic.solve_cubic(inputs)
            real = [i for i in range(3) if out["im"][i] == 0]
            self.assertEqual([cubic.fraction(out["re"][i]) for i in real], [r])
            pair = [i for i in range(3) if i not in real]
            self.assertEqual([cubic.fraction(out["re"][i]) for i in pair], [-p / 2] * 2)
            want = mp.sqrt(mpf(q - p * p / 4))
            self.assertEqual(out["im"][pair[0]], -out["im"][pair[1]])
            self.assertLess(abs(out["im"][pair[1]] - want), mpf(10) ** -110 * want)


class ReferenceTest(unittest.TestCase):
    def test_p_small_matches_cardano(self):
        """x^3 + p x - 1 has one real root u - p/(3u), u^3 = 1/2 + sqrt(1/4 + p^3/27)."""
        for s in (BY_NAME["cubic:one-real-p-small"], BY_NAME["cubic:one-real-p-small@f32"]):
            for inputs in s.records()[:8]:
                out = cubic.solve_cubic(inputs)
                with mp.workdps(2 * precision.DPS):
                    p = mpf(inputs["c"])
                    u = mp.cbrt(mpf(1) / 2 + mp.sqrt(mpf(1) / 4 + p**3 / 27))
                    want = u - p / (3 * u)
                (i,) = [i for i in range(3) if out["im"][i] == 0]
                self.assertLess(abs(out["re"][i] - want), mpf(10) ** -115, inputs)

    def test_real_roots_follow_the_discriminant(self):
        """`near-double-1e-8` holds both a pair of close real roots and a complex pair."""
        counts = set()
        for inputs in BY_NAME["cubic:near-double-1e-8"].records():
            out = cubic.solve_cubic(inputs)
            disc = cubic.discriminant(*(Fraction(inputs[k]) for k in "abcd"))
            n_real = sum(1 for x in out["im"] if x == 0)
            self.assertEqual(n_real, 1 if disc < 0 else 3)
            counts.add(n_real)
        self.assertEqual(counts, {1, 3})

    def test_roots_are_ordered_by_re_then_im(self):
        out = cubic.solve_cubic({"a": 1.0, "b": 1.0, "c": 1.0, "d": 1.0})  # (x + 1)(x^2 + 1)
        self.assertEqual([mp.nint(x) for x in out["re"]], [-1, 0, 0])
        self.assertEqual([mp.nint(x) for x in out["im"]], [0, -1, 1])

    def test_a_zero_discriminant_has_its_roots_in_closed_form(self):
        out = cubic.solve_cubic({"a": 2.0, "b": -3.0, "c": 0.0, "d": 1.0})  # 2 (x - 1)^2 (x + 1/2)
        self.assertEqual([cubic.fraction(x) for x in out["re"]], [Fraction(-1, 2), 1, 1])


class CrossCheckTest(unittest.TestCase):
    CASES = (
        {"a": 1.0, "b": -6.0, "c": 11.0, "d": -6.0},
        {"a": 1.0, "b": 1.0, "c": 1.0, "d": 1.0},
        {"a": 3.0, "b": 0.0, "c": 1e-5, "d": -1.0},
    )

    def test_a_perturbation_beyond_100_digits_fails(self):
        """A real root alone, or a pair's real or imaginary parts together (a pair that is not
        conjugate fails at any size)."""
        for inputs in self.CASES:
            good = cubic.solve_cubic(inputs)
            check_linalg.solve_cubic(inputs, good)
            real = [i for i in range(3) if good["im"][i] == 0]
            pair = [i for i in range(3) if i not in real]
            moves = [("re", [i]) for i in real] + ([("re", pair), ("im", pair)] if pair else [])
            for key, which in moves:
                if good[key][which[0]] == 0:  # the pair of (x + 1)(x^2 + 1) has re = 0
                    continue
                for rel, fails in ((mpf(10) ** -95, True), (mpf(10) ** -105, False)):
                    bad = {**good, key: list(good[key])}
                    for i in which:
                        bad[key][i] *= 1 + rel
                    if fails:
                        with self.assertRaises(check_linalg.CrossCheckError, msg=(inputs, key)):
                            check_linalg.solve_cubic(inputs, bad)
                    else:
                        check_linalg.solve_cubic(inputs, bad)

    def test_a_dropped_or_spurious_root_fails(self):
        inputs = {"a": 1.0, "b": 1.0, "c": 1.0, "d": 1.0}
        good = cubic.solve_cubic(inputs)
        with self.assertRaises(check_linalg.CrossCheckError):
            check_linalg.solve_cubic(inputs, {**good, "im": [mpf(0)] * 3})


if __name__ == "__main__":
    unittest.main()
