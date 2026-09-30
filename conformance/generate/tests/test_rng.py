import math
import unittest

from mpmath import mpf

from gen import precision, rng

precision.setup()

M = 0xFFFFFFFFFFFFFFFF


def vigna_reference(seed: int, n: int) -> list[int]:
    """A line-by-line port of splitmix64.c (Vigna, public domain), uint64_t wrap via masking."""
    x = seed
    out = []
    for _ in range(n):
        x = (x + 0x9E3779B97F4A7C15) & M
        z = x
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M
        out.append(z ^ (z >> 31))
    return out


class SplitMix64Test(unittest.TestCase):
    def test_known_answers_seed_zero(self):
        g = rng.SplitMix64(0)
        got = [g.next_u64() for _ in range(3)]
        self.assertEqual(got, [0xE220A8397B1DCDAF, 0x6E789E6AA1B965F4, 0x06C45D188009454F])

    def test_matches_reference_port(self):
        for seed in (0, 1, rng.SEED, M, 0x9E3779B97F4A7C15):
            g = rng.SplitMix64(seed)
            self.assertEqual([g.next_u64() for _ in range(1000)], vigna_reference(seed, 1000))

    def test_uniform_is_the_top_53_bits(self):
        # Known answers from seed 0: the first outputs above, shifted right by 11, over 2**53.
        g = rng.SplitMix64(0)
        got = [g.uniform().hex() for _ in range(3)]
        self.assertEqual(
            got, ["0x1.c4415072f63b9p-1", "0x1.b9e279aa86e58p-2", "0x1.b117462002500p-6"]
        )

    def test_uniform_is_53_bit_and_in_unit_interval(self):
        g = rng.SplitMix64(7)
        for _ in range(1000):
            u = g.uniform()
            self.assertTrue(0.0 <= u < 1.0)
            self.assertEqual((u * 2**53) % 1, 0)

    def test_streams_are_reproducible_and_independent(self):
        a = [rng.stream(rng.SEED, "theta:1e-8", "theta").next_u64() for _ in range(2)]
        self.assertEqual(a[0], a[1])
        others = {
            rng.stream(rng.SEED, s, p).next_u64()
            for s in ("theta:1e-8", "theta:1e-7")
            for p in ("theta", "axis")
        }
        self.assertEqual(len(others), 4)


def counts(values, edges) -> list[int]:
    return [sum(lo <= v < hi for v in values) for lo, hi in zip(edges, edges[1:], strict=False)]


class SamplingTest(unittest.TestCase):
    def test_log_uniform_stays_in_its_decade(self):
        g = rng.SplitMix64(3)
        lo, hi = mpf(10) ** -5, mpf(10) ** -4
        xs = [rng.log_uniform(g, lo, hi) for _ in range(500)]
        self.assertTrue(all(lo <= mpf(x) < hi for x in xs))
        self.assertGreater(
            sum(x < 10**-4.5 for x in xs), 200
        )  # half the mass below the geometric mean

    def test_unit_vector_is_unit_to_rounding(self):
        g = rng.SplitMix64(11)
        for _ in range(200):
            x, y, z = rng.unit_vector_s2(g)
            self.assertLess(abs(x * x + y * y + z * z - 1), 4 * 2.0**-53)

    def test_unit_vector_is_uniform_on_the_sphere(self):
        """z uniform on [-1, 1) and longitude uniform on [-pi, pi): 4000 fixed draws, so the
        bounds below (five binomial sigmas of 1000 per bin, about 27 each) are not flaky."""
        g = rng.SplitMix64(5)
        points = [rng.unit_vector_s2(g) for _ in range(4000)]
        z_bins = counts([z for _, _, z in points], [-1.0, -0.5, 0.0, 0.5, 1.0000001])
        phi = [math.atan2(y, x) for x, y, _ in points]
        phi_bins = counts(phi, [-math.pi, -math.pi / 2, 0.0, math.pi / 2, math.pi + 1e-9])
        for bins in (z_bins, phi_bins):
            self.assertEqual(sum(bins), 4000)
            self.assertTrue(all(abs(c - 1000) < 135 for c in bins), bins)

    def test_unit_vector_known_answers(self):
        # Two draws from stream (SEED, "theta:1e-8", "axis"), checked against an independent
        # implementation (its own splitmix64 and Fraction rounding): draw order z, then phi.
        g = rng.stream(rng.SEED, "theta:1e-8", "axis")
        got = [tuple(c.hex() for c in rng.unit_vector_s2(g)) for _ in range(2)]
        self.assertEqual(
            got,
            [
                ("-0x1.19a4252741229p-3", "0x1.ce4b596136a4fp-1", "-0x1.a0fda7a43d2d8p-2"),
                ("-0x1.b44adad1ef172p-2", "-0x1.551294eee90fdp-1", "-0x1.3968eb2c783a0p-1"),
            ],
        )


if __name__ == "__main__":
    unittest.main()
