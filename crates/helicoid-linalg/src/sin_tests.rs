//! `Real::sin` is bit-identical to `Real::sin_cos().0`, which is the whole of [`0052`]'s licence to
//! substitute one for the other.
//!
//! Every accuracy figure recorded against a routine that was re-spelled from `sin_cos().0` to `sin`
//! rests on this, so it is a test and not a remark. `libm` is a caret dependency: a patch release
//! may re-shape either function, and the two are **not** identical by construction. Reading
//! `libm` 0.2.16, three places could part:
//!
//! 1. The general path agrees by construction. Both call `rem_pio2` and then select from
//!    `k_sin`/`k_cos` on `n & 3`, and `sincos`'s permutation `(s, c), (c, -s), (-s, -c), (-c, s)`
//!    has the same first component as `sin`'s `k_sin, k_cos, -k_sin, -k_cos`.
//! 2. **The small-argument cuts differ.** `sin` returns `x` below `2^-26`, `sincos` below
//!    `2^-27 sqrt 2`, so in between one returns `x` and the other a Horner. They agree anyway: the
//!    correction is at most `|x|^3/6`, i.e. `x^2/6 <= 2^-54.58` relative, and half an ulp is at
//!    least `2^-53` relative there, so the sum rounds back to `x`. That is an argument about
//!    rounding, not a shared expression, and [`the_band_where_the_cuts_differ`] is the check.
//! 3. **At binary32 the octant arms are differently spelled.** Below `3 pi/4` and positive, `sinf`
//!    evaluates `k_cosf(x - S1_PIO2)` where `sincosf` evaluates `k_cosf(S1_PIO2 - x)` — exact
//!    negations of each other, so they agree only because `k_cosf` reads its argument solely
//!    through `x * x`. A kernel that gained a term linear in `x` would part them, and nothing in
//!    `libm` promises it will not.
//!
//! So binary32 is checked **exhaustively**, over all `2^32` bit patterns, which settles item 3 for
//! good rather than sampling around it. Binary64 is checked over the band of item 2 and a wide
//! seeded sweep; exhaustive is not available there, and item 1 makes the rest structural.
//!
//! [`0052`]: ../../../docs/decisions/0052-real-owes-sin-and-the-corpus-does-not-move.md

use crate::Real;

/// The seeded stream, as `sqrt_tests` uses: fixed, so the sweep is the same on every target.
fn splitmix(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// `sin` against `sin_cos().0`, by bits, with NaN compared by `is_nan` — its payload is the
/// target's (`0018`).
fn agree<S: Real>(x: S) -> bool {
    let (a, b) = (x.sin(), x.sin_cos().0);
    let (av, bv) = (a.value_f64(), b.value_f64());
    if av.is_nan() || bv.is_nan() {
        return av.is_nan() && bv.is_nan();
    }
    av.to_bits() == bv.to_bits()
}

/// Item 2: the band where `sin` returns `x` and `sin_cos` does not, `[2^-27 sqrt 2, 2^-26)`.
///
/// Walked at every `f64` representable in it would be `2^52` values, so it is walked by bits from
/// each end and across the middle instead — the rounding argument in this module's header is
/// uniform over the band, so a defect would not hide in one corner of it.
#[test]
fn the_band_where_the_cuts_differ() {
    let (lo, hi) = (0x3e46_a09e_0000_0000u64, 0x3e50_0000_0000_0000u64);
    let step = (hi - lo) / 20_000;
    let mut bad = 0usize;
    for i in 0..=20_000u64 {
        for b in [lo + i * step, lo + i, hi - 1 - i] {
            let x = f64::from_bits(b);
            bad += usize::from(!agree(x) || !agree(-x));
        }
    }
    assert_eq!(
        bad, 0,
        "`sin` and `sin_cos().0` part inside the differing cuts"
    );
}

/// Binary64 over sixty decades, both signs, past every reduction threshold.
#[test]
fn a_wide_sweep_agrees_at_binary64() {
    let mut s = 0x5369_6e5f_7377_6565u64;
    for _ in 0..200_000 {
        let r = splitmix(&mut s);
        // `m 2^e` with `m` in `[1, 2)`: the exponent sweeps the reduction's regimes, including
        // arguments far above `2^53` where `rem_pio2` takes its slow path.
        let e = ((r >> 52) % 120) as i32 - 60;
        let m = 1.0 + ((r & 0x000f_ffff_ffff_ffff) as f64) / 4_503_599_627_370_496.0;
        let x = m * libm::pow(2.0, f64::from(e));
        assert!(agree(x) && agree(-x), "parted at {x:e}");
    }
    // The named points: zero, the cuts, pi/4, pi/2, pi, and the non-finite arguments.
    for x in [
        0.0,
        -0.0,
        f64::from_bits(0x3e46_a09e_0000_0000),
        f64::from_bits(0x3e50_0000_0000_0000),
        core::f64::consts::FRAC_PI_4,
        core::f64::consts::FRAC_PI_2,
        core::f64::consts::PI,
        f64::MIN_POSITIVE,
        f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ] {
        assert!(agree(x), "parted at the named point {x:e}");
    }
}

/// One eighth of the binary32 domain, by bit pattern, for [`binary32_agrees_exhaustively`].
///
/// `libm::sinf`/`sincosf` directly rather than through `Real`: the `f32` impl of each *is* that one
/// call, so this is the same claim one layer down, and going through `Real` costs a generic
/// dispatch and two `value_f64` round trips per argument — about 4x, on a loop that is already
/// minutes long. [`binary32_agrees_on_a_dense_sample`] covers the `Real` layer.
fn shard(k: u32) -> (u64, Option<u32>) {
    let (lo, hi) = (k << 29, u64::from(k << 29) + (1u64 << 29));
    let (mut bad, mut first) = (0u64, None);
    for b in u64::from(lo)..hi {
        let x = f32::from_bits(b as u32);
        let (a, c) = (libm::sinf(x), libm::sincosf(x).0);
        // NaN by `is_nan`: its payload is the target's (`0018`).
        let parts = if a.is_nan() || c.is_nan() {
            a.is_nan() != c.is_nan()
        } else {
            a.to_bits() != c.to_bits()
        };
        if parts {
            bad += 1;
            first = first.or(Some(b as u32));
        }
    }
    (bad, first)
}

/// Binary32, **exhaustively**: all `2^32` bit patterns across the eight tests below, which is the
/// only way to settle the differently-spelled octant arms of item 3.
///
/// Eight **separate tests** of `2^29` and not one loop over eight shards: the whole domain is about
/// four minutes in release — most patterns are `|x| > 2^24`, where `rem_pio2f` takes its slow path
/// — and a single test of that length exceeds nextest's 180 s timeout, which is how the first two
/// drafts of this module failed. Separate tests run concurrently, so the wall time is one shard's;
/// a loop inside one test does not, which is the trap the second draft fell into after the comment
/// had already been written to claim otherwise.
///
/// All `#[ignore]`d; [`binary32_agrees_on_a_dense_sample`] is what `just test` runs. Rerun these
/// whenever `libm` moves:
/// `cargo nextest run --release -p helicoid-linalg -E 'test(exhaustively)' --run-ignored all`.
macro_rules! exhaustive_shard {
    ($($name:ident = $k:expr;)+) => {$(
        #[test]
        #[ignore = "exhaustive over 2^29 of 2^32; the dense sample is what `just test` runs"]
        fn $name() {
            let (bad, first) = shard($k);
            assert_eq!(bad, 0, "{bad} arguments part, first at {first:?}");
        }
    )+};
}

exhaustive_shard! {
    binary32_agrees_exhaustively_0 = 0;
    binary32_agrees_exhaustively_1 = 1;
    binary32_agrees_exhaustively_2 = 2;
    binary32_agrees_exhaustively_3 = 3;
    binary32_agrees_exhaustively_4 = 4;
    binary32_agrees_exhaustively_5 = 5;
    binary32_agrees_exhaustively_6 = 6;
    binary32_agrees_exhaustively_7 = 7;
}

/// Binary32 on a dense sample: every `2^11`th bit pattern, plus the whole of the two octant bands
/// item 3 names, where the spellings differ.
#[test]
fn binary32_agrees_on_a_dense_sample() {
    let mut bad = 0u64;
    let mut b = 0u32;
    loop {
        bad += u64::from(!agree(f32::from_bits(b)));
        match b.checked_add(1 << 11) {
            Some(n) => b = n,
            None => break,
        }
    }
    // `(pi/4, 3pi/4]` and `(3pi/4, 5pi/4]`, the two bands whose arms are spelled as negations, at
    // every representable `f32`: `0x3f490fda` to `0x407b53d1` is under 1.5e7 values.
    for b in 0x3f49_0fdau32..=0x407b_53d1 {
        let x = f32::from_bits(b);
        bad += u64::from(!agree(x) || !agree(-x));
    }
    assert_eq!(bad, 0, "`sinf` and `sincosf().0` part");
}

/// `Dual`'s `sin` is its `sin_cos().0` by construction, and the derivative still rides the cosine.
#[test]
fn dual_keeps_its_derivative() {
    use crate::Dual;
    let x = Dual::<f64, 2>::variable(0.7, 0);
    let (s, c) = x.sin_cos();
    let got = x.sin();
    assert_eq!(
        got.v.to_bits(),
        s.v.to_bits(),
        "the value is `sin_cos().0`'s"
    );
    assert_eq!(
        got.d[0].to_bits(),
        s.d[0].to_bits(),
        "and so is the derivative"
    );
    // The seeded lane is `cos(0.7)`, which is `sin_cos`'s other half: the chain rule is intact.
    assert_eq!(got.d[0].to_bits(), c.v.to_bits(), "d sin / dx is cos x");
    assert_eq!(
        got.d[1].to_bits(),
        0.0_f64.to_bits(),
        "the unseeded lane is zero"
    );
}
