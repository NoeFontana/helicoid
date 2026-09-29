"""Deterministic sampling. splitmix64 lives here: `random`'s algorithm is a CPython detail."""

import hashlib

from mpmath import mp, mpf

from .precision import to_f64

SEED = 0x68656C69636F6964  # "helicoid"
_MASK = (1 << 64) - 1


class SplitMix64:
    def __init__(self, seed: int) -> None:
        self.state = seed & _MASK

    def next_u64(self) -> int:
        self.state = (self.state + 0x9E3779B97F4A7C15) & _MASK
        z = self.state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & _MASK
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & _MASK
        return z ^ (z >> 31)

    def uniform(self) -> float:
        """Uniform in [0, 1) from the top 53 bits; exactly representable."""
        return (self.next_u64() >> 11) * 2.0**-53


def stream(seed: int, *labels: str) -> SplitMix64:
    """An independent stream per (seed, labels): a new stratum never reorders another's inputs."""
    digest = hashlib.sha256("/".join((f"{seed:016x}", *labels)).encode()).digest()
    return SplitMix64(int.from_bytes(digest[:8], "big"))


def log_uniform(rng: SplitMix64, lo, hi) -> float:
    """Log-uniform in [lo, hi) (mpf bounds), rounded to binary64; a rounding escape is an error."""
    theta = to_f64(lo * (hi / lo) ** rng.uniform())
    if not lo <= mpf(theta) < hi:
        raise ArithmeticError(f"{theta!r} rounded outside [{lo}, {hi})")
    return theta


def unit_vector_s2(rng: SplitMix64) -> tuple[float, float, float]:
    """Uniform on S2 (z uniform, phi uniform) at working precision, then rounded to binary64."""
    z = 2 * mpf(rng.uniform()) - 1
    phi = 2 * mp.pi * rng.uniform()
    r = mp.sqrt(1 - z * z)
    return to_f64(r * mp.cos(phi)), to_f64(r * mp.sin(phi)), to_f64(z)


def unit_quaternion_s3(rng: SplitMix64) -> tuple:
    """Uniform on S3 (Shoemake 1992): a Haar-random rotation, at working precision, unrounded."""
    u1, u2, u3 = (mpf(rng.uniform()) for _ in range(3))
    a, b, s, t = mp.sqrt(1 - u1), mp.sqrt(u1), 2 * mp.pi * u2, 2 * mp.pi * u3
    return a * mp.sin(s), a * mp.cos(s), b * mp.sin(t), b * mp.cos(t)
