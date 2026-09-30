//! `Real::sqrt` for `f64` and `f32` against a correctly rounded integer square root, bit for bit,
//! and a pinned digest of the results (`0018`).
//!
//! `libm`'s `arch` feature swaps the body of `sqrt` for the target's instruction (x86_64, aarch64).
//! Both are exactly rounded, so the bits cannot move; this test is what would notice if they did.
//! The oracle uses integers only (no float square root, no `libm`): the value `m 2^e` is put in the
//! form `m' 2^(2e')`, `m' 2^72` gets an integer root `r` on `u128`, and `r` with its remainder as a
//! sticky bit is rounded to `p` bits, ties to even. The stream is a fixed splitmix64 sequence, so the
//! digest is the same on every target; it does not depend on the `arch` feature either.
//!
//! NaN is compared by `is_nan` only: its sign and payload are the target's (`0018`, Context).
//!
//! The oracle's sticky path (a root whose bits below the guard bit are zero, with a nonzero
//! remainder) is reached by the `f64` stream only; the `f32` stream does not aim at it.

use crate::Real;

/// The format of a binary float, by its significand precision and exponent width.
#[derive(Clone, Copy)]
struct Fmt {
    /// Significand precision `p`, the leading bit included.
    p: u32,
    /// Exponent field width.
    ebits: u32,
}

const F64: Fmt = Fmt { p: 53, ebits: 11 };
const F32: Fmt = Fmt { p: 24, ebits: 8 };

impl Fmt {
    fn frac_bits(self) -> u32 {
        self.p - 1
    }
    fn bias(self) -> i32 {
        (1 << (self.ebits - 1)) - 1
    }
    /// The all-ones exponent field: infinity and NaN.
    fn emax_field(self) -> u64 {
        (1 << self.ebits) - 1
    }
    fn sign(self) -> u64 {
        1 << (self.frac_bits() + self.ebits)
    }
    fn pack(self, exp_field: u64, frac: u64) -> u64 {
        (exp_field << self.frac_bits()) | frac
    }
    /// `(m, e)` with `x = m 2^e`, `m > 0`, for a positive finite `x`.
    fn decode(self, bits: u64) -> (u64, i32) {
        let frac = bits & ((1 << self.frac_bits()) - 1);
        let field = bits >> self.frac_bits();
        let smallest = 1 - self.bias() - self.frac_bits() as i32;
        if field == 0 {
            (frac, smallest)
        } else {
            (frac | (1 << self.frac_bits()), field as i32 + smallest - 1)
        }
    }
}

/// `floor(sqrt(n))`, digit by digit.
fn isqrt(n: u128) -> u128 {
    let (mut rem, mut root, mut bit) = (n, 0_u128, 1_u128 << 126);
    while bit > rem {
        bit >>= 2;
    }
    while bit != 0 {
        if rem >= root + bit {
            rem -= root + bit;
            root = (root >> 1) + bit;
        } else {
            root >>= 1;
        }
        bit >>= 2;
    }
    root
}

fn bit_len(v: u128) -> u32 {
    128 - v.leading_zeros()
}

/// The bits of `(v + s) 2^q` rounded to `fmt`, ties to even, for `0 <= s < 1` (`sticky` is `s > 0`).
/// The result must be normal.
fn round_pack(v: u128, sticky: bool, q: i32, fmt: Fmt) -> u64 {
    let len = bit_len(v);
    let (mut sig, mut q) = if len <= fmt.p {
        assert!(!sticky, "an inexact value narrower than the format");
        ((v as u64) << (fmt.p - len), q - (fmt.p - len) as i32)
    } else {
        let sh = len - fmt.p;
        let mut sig = (v >> sh) as u64;
        let guard = (v >> (sh - 1)) & 1 == 1;
        let low = v & ((1 << (sh - 1)) - 1) != 0 || sticky;
        if guard && (low || sig & 1 == 1) {
            sig += 1;
        }
        (sig, q + sh as i32)
    };
    if sig == 1 << fmt.p {
        sig >>= 1;
        q += 1;
    }
    let field = q + fmt.p as i32 - 1 + fmt.bias();
    assert!(
        (1..fmt.emax_field() as i32).contains(&field),
        "result not normal"
    );
    fmt.pack(field as u64, sig - (1 << fmt.frac_bits()))
}

/// The correctly rounded square root of a positive finite `x`, as bits.
fn oracle(bits: u64, fmt: Fmt) -> u64 {
    let (m, e) = fmt.decode(bits);
    let mut s = fmt.p + 1 - (64 - m.leading_zeros());
    if (e - s as i32) & 1 != 0 {
        s += 1;
    }
    let scaled = u128::from(m) << (s + 72);
    let r = isqrt(scaled);
    round_pack(r, scaled != r * r, (e - s as i32) / 2 - 36, fmt)
}

struct SplitMix64(u64);

impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform in `lo..=hi` (modulo bias below `2^-40`, which a fixed stream does not care about).
    fn range(&mut self, lo: i64, hi: i64) -> i64 {
        lo + (self.next() % (hi - lo + 1) as u64) as i64
    }
}

/// One positive finite bit pattern of `fmt`, of the class `i % 12`.
fn case(i: u64, rng: &mut SplitMix64, fmt: Fmt) -> u64 {
    let fb = fmt.frac_bits();
    let frac_mask = (1_u64 << fb) - 1;
    let emax = fmt.emax_field() as i64;
    // Half the exponent range, less a margin, in which a `k^2 2^(2j)` stays normal.
    let jmax = i64::from((fmt.bias() - 1) / 2) - 32;
    let kmax = isqrt((1 << fmt.p) - 1) as i64;
    match i % 12 {
        // Every normal exponent equally likely; and raw random bits with the sign cleared.
        0 => fmt.pack(rng.range(1, emax - 1) as u64, rng.next() & frac_mask),
        1 => {
            let b = rng.next() & (fmt.sign() - 1);
            if b >> fb == fmt.emax_field() {
                b ^ (1 << (fb + fmt.ebits - 1))
            } else {
                b
            }
        }
        // Subnormals: a random fraction, and a few set bits.
        2 => fmt.pack(0, rng.next() & frac_mask | 1),
        3 => fmt.pack(0, (1 << rng.range(0, i64::from(fb) - 1)) | (rng.next() & 3)),
        // Exact squares `k^2 2^(2j)`, and one ulp either side.
        4 | 5 => {
            let k = rng.range(1, kmax) as u128;
            let j = rng.range(-jmax, jmax) as i32;
            let exact = round_pack(k * k, false, 2 * j, fmt);
            match (i % 12, rng.next() % 3) {
                (4, _) | (_, 0) => exact,
                (_, 1) => exact - 1,
                _ => exact + 1,
            }
        }
        // Powers of two and their neighbours: both sides of every `2^k`.
        6 => {
            let p2 = fmt.pack(rng.range(1, emax - 1) as u64, 0);
            match rng.next() % 3 {
                0 => p2,
                1 => p2 - 1,
                _ => p2 + 1,
            }
        }
        // A midpoint squared, `(n + 1/2)^2`, rounded to `fmt`: a root near a rounding boundary.
        7 => {
            let n = (1_u64 << fb) | (rng.next() & frac_mask);
            let sq = u128::from(2 * n + 1) * u128::from(2 * n + 1);
            let want = rng.range(1 - i64::from(fmt.bias()) + 8, i64::from(fmt.bias()) - 8) as i32;
            let q = (want - (bit_len(sq) as i32 - 1)).div_euclid(2) * 2 - 2;
            round_pack(sq, false, q, fmt)
        }
        // The largest and the smallest exponents.
        8 => fmt.pack((emax - 1 - rng.range(0, 7)) as u64, rng.next() & frac_mask),
        9 => fmt.pack(rng.range(1, 8) as u64, rng.next() & frac_mask),
        // A fraction of all ones below `t` bits: just under the next power of two.
        10 => {
            let t = rng.range(1, i64::from(fb)) as u32;
            fmt.pack(rng.range(1, emax - 1) as u64, (1 << t) - 1)
        }
        // Small integers and their exact squares, one ulp above.
        _ => {
            let k = rng.range(1, kmax.min(1 << 16)) as u128;
            let exact = round_pack(k * k, false, 0, fmt);
            if rng.next() & 1 == 0 {
                exact
            } else {
                exact + 1
            }
        }
    }
}

/// Fixed edge cases, in front of the stream: every power of two (with both neighbours), the extreme
/// finite values, and the squares of the first integers (with the neighbour above).
fn edges(fmt: Fmt) -> impl Iterator<Item = u64> {
    let fb = fmt.frac_bits();
    let powers = (1..fmt.emax_field()).flat_map(move |f| {
        let p2 = fmt.pack(f, 0);
        [p2 - 1, p2, p2 + 1]
    });
    let subnormal_powers = (0..fb).map(move |t| 1_u64 << t);
    let extremes = [
        1,
        (1 << fb) - 1,
        1 << fb,
        (1 << fb) + 1,
        fmt.pack(fmt.emax_field() - 1, (1 << fb) - 1),
        fmt.pack(fmt.emax_field() - 1, (1 << fb) - 2),
        fmt.pack(fmt.bias() as u64, 0),
        fmt.pack(fmt.bias() as u64, 1),
        fmt.pack(fmt.bias() as u64 - 1, (1 << fb) - 1),
    ];
    let squares = (1_u128..=4096).flat_map(move |k| {
        let sq = round_pack(k * k, false, 0, fmt);
        [sq, sq + 1, sq - 1]
    });
    powers
        .chain(subnormal_powers)
        .chain(extremes)
        .chain(squares)
}

const CASES: u64 = 1_000_000;

/// The `edges`, then `CASES` of the stream; each result is compared with the oracle and folded
/// into an FNV-1a digest of the result bits.
fn run(fmt: Fmt, seed: u64, sqrt: fn(u64) -> u64) -> (u64, u64) {
    let mut rng = SplitMix64(seed);
    let mut digest = 0xCBF2_9CE4_8422_2325_u64;
    let mut count = 0_u64;
    let mut visit = |bits: u64| {
        let got = sqrt(bits);
        let want = oracle(bits, fmt);
        assert_eq!(
            got, want,
            "sqrt({bits:#x}): got {got:#x}, correctly rounded {want:#x}"
        );
        digest = (digest ^ got).wrapping_mul(0x0000_0100_0000_01B3);
        count += 1;
    };
    edges(fmt).for_each(&mut visit);
    for i in 0..CASES {
        visit(case(i, &mut rng, fmt));
    }
    (digest, count)
}

fn sqrt64(bits: u64) -> u64 {
    Real::sqrt(f64::from_bits(bits)).to_bits()
}

fn sqrt32(bits: u64) -> u64 {
    u64::from(Real::sqrt(f32::from_bits(bits as u32)).to_bits())
}

/// The digest of the results over the stream; the same on every target, with `arch` on or off.
const DIGEST_F64: u64 = 0xaa04_43e9_153d_eb6d;
const DIGEST_F32: u64 = 0x14c6_bc71_6b7e_3361;

#[test]
fn sqrt_f64_is_correctly_rounded_over_the_stream() {
    let (digest, count) = run(F64, 0x5EED_0018_F640_0001, sqrt64);
    assert!(count >= CASES);
    assert_eq!(
        digest, DIGEST_F64,
        "digest of the f64 results: {digest:#018x}"
    );
}

#[test]
fn sqrt_f32_is_correctly_rounded_over_the_stream() {
    let (digest, count) = run(F32, 0x5EED_0018_F320_0001, sqrt32);
    assert!(count >= CASES);
    assert_eq!(
        digest, DIGEST_F32,
        "digest of the f32 results: {digest:#018x}"
    );
}

/// The oracle on values whose root is known, by hand: `sqrt(4)`, `sqrt(1)`, `sqrt(2)` (the
/// correctly rounded constants `0x3FF6A09E667F3BCD`, `0x3FB504F3`) and `sqrt` of the smallest
/// subnormal: `2^-537` exactly for `f64` (`2^-1074`), `sqrt(2) 2^-75` for `f32` (`2^-149`).
#[test]
fn the_oracle_agrees_with_known_roots() {
    let (one, two, four) = (
        0x3FF0_0000_0000_0000,
        0x4000_0000_0000_0000,
        0x4010_0000_0000_0000,
    );
    assert_eq!(oracle(four, F64), two);
    assert_eq!(oracle(one, F64), one);
    assert_eq!(oracle(two, F64), 0x3FF6_A09E_667F_3BCD);
    assert_eq!(oracle(0x4000_0000, F32), 0x3FB5_04F3);
    assert_eq!(oracle(1, F64), 486 << 52);
    assert_eq!(oracle(1, F32), (52 << 23) | 0x35_04F3);
}

/// Zeros, infinity and NaN. `sqrt(-0) = -0` (IEEE 754-2019 §5.4.1); a NaN stays a NaN.
#[test]
fn sqrt_of_zeros_infinity_and_nan() {
    for x in [0.0_f64, -0.0, f64::INFINITY] {
        assert_eq!(Real::sqrt(x).to_bits(), x.to_bits());
    }
    for x in [0.0_f32, -0.0, f32::INFINITY] {
        assert_eq!(Real::sqrt(x).to_bits(), x.to_bits());
    }
    for bits in [
        0x7FF8_0000_0000_0000_u64,
        0xFFF8_0000_0000_0001,
        0x7FF0_0000_0000_0001,
    ] {
        assert!(Real::sqrt(f64::from_bits(bits)).is_nan());
    }
    for bits in [0x7FC0_0000_u32, 0xFFC0_0001, 0x7F80_0001] {
        assert!(Real::sqrt(f32::from_bits(bits)).is_nan());
    }
}

/// A negative argument is a domain error: NaN in a release build (its sign is the target's).
#[cfg(not(debug_assertions))]
#[test]
fn sqrt_of_a_negative_is_nan_in_release() {
    for x in [-1.0_f64, -f64::MIN_POSITIVE, -5e-324, f64::NEG_INFINITY] {
        assert!(Real::sqrt(x).is_nan());
    }
    for x in [-1.0_f32, -f32::MIN_POSITIVE, -1e-45, f32::NEG_INFINITY] {
        assert!(Real::sqrt(x).is_nan());
    }
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "negative argument")]
fn sqrt_of_a_negative_f64_asserts_in_debug() {
    let _ = Real::sqrt(-1.0_f64);
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "negative argument")]
fn sqrt_of_a_negative_f32_asserts_in_debug() {
    let _ = Real::sqrt(f32::NEG_INFINITY);
}

/// A sign-reading operation on a NaN from `sqrt` gets the target's sign (`0018`, Decision 2): the
/// result is `+-|x|`, never anything else, and which sign is not pinned.
#[cfg(not(debug_assertions))]
#[test]
fn copysign_of_a_sqrt_nan_is_plus_or_minus_x() {
    let c = Real::copysign(3.0_f64, Real::sqrt(-1.0_f64));
    assert_eq!(c.abs().to_bits(), 3.0_f64.to_bits());
    let c = Real::copysign(3.0_f32, Real::sqrt(-1.0_f32));
    assert_eq!(c.abs().to_bits(), 3.0_f32.to_bits());
}
