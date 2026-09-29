//! Tests of `solve_cubic`: omnisac's own cases, the strata of `docs/PHASE2.md` §6 against mpmath and
//! against planted roots, the slot order and the bits of omnisac's algorithm, the value path of
//! `Dual`, its derivative, two lanes, and the mask.
//!
//! The mpmath rows (`MP`): each polynomial is the exact expansion of dyadic roots (or a listed exact
//! cubic) with every coefficient a binary64 number, so no row rounds its input; a row with a
//! coefficient that is not an `f32` is skipped at `f32`, and both precisions solve the polynomial
//! mpmath solved. mpmath (120 digits) finds the roots with `mp.polyroots` (Durand-Kerner), an
//! independent method, and checks each simple root by Horner evaluation of the exact coefficients
//! (residual below `1e-80` of its terms); a polynomial with an exactly zero discriminant (a repeated
//! root, where Durand-Kerner converges only linearly) takes sympy's exact algebraic roots to 130
//! digits, not `mp.polyroots`. A row is `(stratum, [a, b, c, d], real parts, imaginary parts)`,
//! rounded to binary64; the generator is not committed. The planted rows (`planted`) have exact
//! dyadic roots by construction.
//!
//! A valid slot must lie within `min(K u ((W + Wd) / P + R), cap) + loss + pair` of a root `z`, and
//! every real root within it of a slot:
//!
//! - `(W + Wd) / P` is the first-order effect of rounding at two stages, `P = |p'(z)| = prod |z -
//!   z_k|` (monic): `W = |z|^3 + |c_1| |z|^2 + |c_2| |z| + |c_3|` the normalisation `b / a`, `c / a`,
//!   `d / a`; `Wd = |t| pw + qw` the formation of `p = C - B^2/3` and `q = 2 B^3/27 - B C/3 + D` of
//!   the depressed cubic `t^3 + p t + q` (`t = z + B/3`, `pw = |C| + B^2/3`, `qw = 2 |B|^3/27 + |B C|/3
//!   + |D|`), which `W` misses where the roots lie close together far from the origin;
//! - `R`, the largest root, is the absolute error of the shift `B/3` and of the trigonometric arm's
//!   `2 r cos`;
//! - `cap` is the cluster rounding spreads a multiple root into, `4 sqrt(u) R` for a double and
//!   `4 u^(1/3) R` for a triple root or a real root with a pair near it;
//! - `loss`, in the one-real-root arm: the root is `cbrt(h + s) + cbrt(h - s)`, and the rounding
//!   error `4 u (h + s)` of `h - s` reaches it as `4 u (h + s) / (3 v^2)`, `v = p / (3 cbrt(h +
//!   s))`, saturating at `|v|` plus the noise `cbrt(4 u (h + s))` (`cardano_loss`; `h = |q|/2`,
//!   `s = sqrt(h^2 + p^3/27)`, `p`, `q` of the depressed cubic);
//! - `pair`: `|Im z|` of a complex pair inside the band's radius `2 sqrt(band) R`, which is
//!   reported as a real double root.
//!
//! `K` is 16. Over the fixture and 3e5 planted rows per scalar the worst `err / (u ((W + Wd) / P +
//! R))` is 4.24 (three real roots at `f64`) where no other term applies, and the worst `err /
//! bound` is 0.88 (a double root at its `cap`, `f32`). The band is relative to the summands of
//! `disc`, not to `B`, `C`, `D`, so where the rounding error of `disc` reaches the band or `disc`
//! itself (`regular`) the model claims nothing: a real pair is answered as one root and a complex
//! pair as a real double root, which moves the real root too (`a_repeated_root_can_be_dropped`).
//!
//! `GOLDEN` is 36 polynomials with random binary64 coefficients and the bits of their roots. They
//! were computed by omnisac d3be7b7 `poly::solve_cubic` on `libm` 0.2.16 with this port's tolerances
//! (`2^-46`, `2^-40`, `2^-46`), drawn from 4e4 random rows: the bits agree outside the trigonometric
//! arm, and inside it the `acos` form differs by up to 5 `u` of the largest root on those rows, and
//! by 11 `u` at the most observed over 5e6 random polynomials (a sample maximum, not a bound).
//! The differential of omnisac's own `poly::solve_cubic` against this port (1e6 polynomials per seed,
//! 3 seeds, 20 families of random, planted, strip, subnormal, huge and non-finite coefficients, on
//! `libm` and std math) is not committed: 95.2% bit-identical, the rest the leading-strip and
//! band-strip tolerance changes and the `acos` ulps.
//!
//! The proptests run 512 cases (`cfg`); 10^6 from a fixed seed is
//!
//! ```text
//! PROPTEST_CASES=1000000 PROPTEST_RNG_SEED=1 cargo nextest run --release -p helicoid-linalg cubic_tests
//! ```

use crate::cubic::{acos, pi};
use crate::linalg_tests::{cfg, unit, Lane, D};
use crate::tests::L2;
use crate::{solve_cubic, Dual, Mask, Precision, Real};
use core::f64::consts::{PI, SQRT_2};
use proptest::prelude::*;
use std::vec::Vec;

const K: f64 = 16.0;

/// Every slot of `solve_cubic`'s answer for `c` as `f64`, and its mask. Checks what holds for every
/// input: a slot that is not valid is `+0`, a valid one is finite and never `-0`.
fn solve<S: Lane>(c: [f64; 4]) -> ([f64; 3], [bool; 3]) {
    let [a, b, cc, d] = c.map(|x| S::make(x, [0.0; 2]));
    let (r, m) = solve_cubic(a, b, cc, d);
    let (r, m) = (r.0.map(|x| x.lane(0)), m.map(Mask::all));
    for (&x, &v) in r.iter().zip(&m) {
        let ok = if v {
            x.is_finite() && (x.abs() > 0.0 || x.to_bits() == 0)
        } else {
            x.to_bits() == 0
        };
        assert!(ok, "{c:?}: {r:?} {m:?}");
    }
    (r, m)
}

/// The valid roots, ascending.
fn valid<S: Lane>(c: [f64; 4]) -> Vec<f64> {
    let (r, m) = solve::<S>(c);
    let mut v: Vec<f64> = r.iter().zip(&m).filter(|p| *p.1).map(|p| *p.0).collect();
    v.sort_by(f64::total_cmp);
    v
}

/// The leading-coefficient floor and the discriminant band of the rustdoc table, written here
/// independently of the code: the tests pin the table, not the implementation.
fn table<S: Real>() -> (f64, f64) {
    match S::PRECISION {
        Precision::F64 => (2f64.powi(-46), 2f64.powi(-40)),
        Precision::F32 => (2f64.powi(-17), 2f64.powi(-11)),
    }
}

fn is_cubic<S: Lane>(c: [f64; 4]) -> bool {
    let c = c.map(|x| S::make(x, [0.0; 2]).lane(0));
    c[0].abs() > table::<S>().0 * c.iter().fold(1.0_f64, |m, x| m.max(x.abs()))
}

/// omnisac's `poly` and `fundamental` cases: three real roots, one, none, an exact triple root and
/// `x^3` in every sign of zero (the root is `+0`), non-finite input, a leading coefficient below its
/// floor; at `f64` also its repeated roots, which must not be dropped, and a pair `1e-7` apart.
fn omnisac_cases<S: Lane>() {
    let same = |c: [f64; 4], want: &[f64]| {
        let got = valid::<S>(c);
        let close = |(g, w): (&f64, &f64)| (g - w).abs() <= 64.0 * unit::<S>() * w.abs().max(1.0);
        assert!(
            got.len() == want.len() && got.iter().zip(want).all(close),
            "{c:?}: {got:?}"
        );
    };
    same([1.0, -6.0, 11.0, -6.0], &[1.0, 2.0, 3.0]);
    same([1.0, -2.0, -5.0, 6.0], &[-2.0, 1.0, 3.0]);
    same([1.0, 1.0, 1.0, 1.0], &[-1.0]);
    same([1.0, 0.0, 0.0, -1.0], &[1.0]);
    same([1.0, -6.0, 12.0, -8.0], &[2.0]);
    for c in [
        [1e-20, 1.0, 2.0, 3.0],
        [0.0, 1.0, -3.0, 2.0],
        [0.0; 4],
        [1e-310; 4],
        [1e-40; 4],
    ] {
        assert!(valid::<S>(c).is_empty(), "{c:?} is not a cubic");
    }
    for (i, bad) in
        (0..4).flat_map(|i| [f64::NAN, f64::INFINITY, f64::NEG_INFINITY].map(|b| (i, b)))
    {
        let mut c = [1.0, 2.0, 3.0, 4.0];
        c[i] = bad;
        assert!(valid::<S>(c).is_empty(), "{c:?}");
    }
    for bits in 0..16u32 {
        let z = |k: u32| if bits >> k & 1 == 0 { 0.0 } else { -0.0 };
        let c = [if bits & 8 == 0 { 1.0 } else { -1.0 }, z(0), z(1), z(2)];
        let (r, m) = solve::<S>(c);
        assert!(
            m == [true, false, false] && r[0].to_bits() == 0,
            "{c:?}: {r:?} {m:?}"
        );
    }
    if S::PRECISION == Precision::F64 {
        // `2 p r` underflows: the `acos` argument is NaN, and the mask says so instead of a NaN slot.
        assert!(solve::<S>([1.0, 0.0, -1e-250, 0.0]).1 == [false; 3]);
        let eval = |c: [f64; 4], x: f64| ((c[0] * x + c[1]) * x + c[2]) * x + c[3];
        for (r, s) in [
            (0.1, 2.7),
            (0.3, 5.9),
            (0.3, 3.1),
            (0.1, 11.7),
            (0.3, 1.618),
        ] {
            let c = [1.0, -(2.0 * r + s), r * r + 2.0 * r * s, -(r * r * s)];
            let got = valid::<S>(c);
            let has = |t: f64| got.iter().any(|x| (x - t).abs() < 1e-6);
            assert!(got.len() >= 2 && has(r) && has(s), "{c:?}: {got:?}");
            assert!(got
                .iter()
                .all(|&x| eval(c, x).abs() < 1e-6 * (1.0 + s).powi(3)));
        }
        // `p^3` and `q^2` underflow: the triple-root arm, `cbrt(-q)`, not the `3q/p` of a double root.
        let (r, m) = solve::<S>([1.0, 0.0, 1e-110, 1e-200]);
        assert!(
            m == [true, false, false] && (r[0] + 2.1544e-67).abs() < 1e-70,
            "{r:?}"
        );
        let (r, e, s) = (1.0, 1e-7, 2.0);
        let c = [
            1.0,
            -(2.0 * r + e + s),
            r * (r + e) + (2.0 * r + e) * s,
            -(r * (r + e) * s),
        ];
        let got = valid::<S>(c);
        let has = |t: f64, tol: f64| got.iter().any(|x| (x - t).abs() < tol);
        assert!(got.len() >= 2 && has(r, 1e-3) && has(s, 1e-9), "{got:?}");
    }
}

#[test]
fn omnisac_cases_hold_at_every_scalar() {
    omnisac_cases::<f64>();
    omnisac_cases::<f32>();
    omnisac_cases::<D>();
}

/// The leading coefficient must exceed `max(scale, 1) 2^-46` (`2^-17` at `f32`; omnisac's `1e-14`
/// was 0.7 of it at `f64`), `scale` the largest coefficient, whichever it is. The band of `disc` is `2^-40` (`2^-11`): `x^3 - 3x + q` with
/// `q = 2 sqrt(1 + delta)` has `disc = delta`, and inside the band the trigonometric arm reports
/// three roots, above it Cardano one; omnisac's `1e-12` is 1.1 of the `f64` band, so
/// `delta = 0.98e-12` moved from three roots to one.
fn thresholds<S: Lane>() {
    let (floor, band) = table::<S>();
    let n = |c: [f64; 4]| valid::<S>(c).len();
    assert!(n([floor * 1.001, 1.0, 1.0, -1.0]) >= 1 && n([floor * 0.999, 1.0, 1.0, -1.0]) == 0);
    let q = |delta: f64| [1.0, 0.0, -3.0, 2.0 * (1.0 + delta).sqrt()];
    assert_eq!((n(q(0.9 * band)), n(q(1.1 * band))), (3, 1));
    if S::PRECISION == Precision::F64 {
        assert_eq!((n([1.2e-14, 1.0, 1.0, -1.0]), n(q(0.98e-12))), (0, 1));
    }
    // A monic coefficient beyond `1 / floor` is rejected: `x^3 - c x` has the roots `0` and `+-sqrt(c)`.
    let (ok, over) = (
        n([1.0, 0.0, -0.9 / floor, 0.0]),
        n([1.0, 0.0, -1.1 / floor, 0.0]),
    );
    assert!(ok >= 1 && over == 0, "{ok} {over}");
    // Each of `b`, `c`, `d` in turn the largest coefficient (8): the floor is `8 floor`, and a
    // leading coefficient at it is rejected (`<=`) and just above it accepted.
    for i in 1..4 {
        let mut c = [0.0, 1.0, 1.0, 1.0];
        c[i] = -8.0;
        let with = |a: f64| n([a, c[1], c[2], c[3]]);
        let at = 8.0 * floor;
        assert_eq!((with(at), with(at * 1.001) >= 1), (0, true), "{c:?}");
    }
}

#[test]
fn thresholds_sit_at_their_table_values() {
    thresholds::<f64>();
    thresholds::<f32>();
}

/// `(stratum, coefficients, real parts, imaginary parts)` of the roots.
type Row = (&'static str, [f64; 4], [f64; 3], [f64; 3]);

#[rustfmt::skip]
const MP: &[Row] = &[
    ("distinct", [1.0, -6.0, 11.0, -6.0], [1.0, 2.0, 3.0], [0.0, 0.0, 0.0]),
    ("distinct", [1.0, -2.0, -5.0, 6.0], [-2.0, 1.0, 3.0], [0.0, 0.0, 0.0]),
    ("distinct", [1.0, -0.75, 0.171875, -0.01171875], [0.125, 0.25, 0.375], [0.0, 0.0, 0.0]),
    ("distinct", [1.0, -3.75, -81.625, 20.625], [-7.5, 0.25, 11.0], [0.0, 0.0, 0.0]),
    ("distinct", [1.0, 0.0009765625, -9.5367431640625e-06, 7.450580596923828e-09], [-0.00390625, 0.0009765625, 0.001953125], [0.0, 0.0, 0.0]),
    ("distinct", [1.0, -2560.0, -5767168.0, 7516192768.0], [-2048.0, 1024.0, 3584.0], [0.0, 0.0, 0.0]),
    ("distinct", [1.0, 0.0, -1.0, 0.0], [-1.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
    ("distinct", [1.0, -102.0, -9799.0, 999900.0], [-99.0, 100.0, 101.0], [0.0, 0.0, 0.0]),
    ("distinct", [1.0, -0.3125, -0.109375, 0.0078125], [-0.25, 0.0625, 0.5], [0.0, 0.0, 0.0]),
    ("distinct", [1.0, -4.0078125, 5.0234375, -2.015625], [1.0, 1.0078125, 2.0], [0.0, 0.0, 0.0]),
    ("near-double", [1.0, -4.125, 5.375, -2.25], [1.0, 1.125, 2.0], [0.0, 0.0, 0.0]),
    ("near-double", [1.0, -4.03125, 5.09375, -2.0625], [1.0, 1.03125, 2.0], [0.0, 0.0, 0.0]),
    ("near-double", [1.0, -4.00390625, 5.01171875, -2.0078125], [1.0, 1.00390625, 2.0], [0.0, 0.0, 0.0]),
    ("near-double", [1.0, -4.000244140625, 5.000732421875, -2.00048828125], [1.0, 1.000244140625, 2.0], [0.0, 0.0, 0.0]),
    ("near-double", [1.0, -4.0000152587890625, 5.0000457763671875, -2.000030517578125], [1.0, 1.0000152587890625, 2.0], [0.0, 0.0, 0.0]),
    ("near-double", [1.0, -4.000000953674316, 5.000002861022949, -2.000001907348633], [1.0, 1.0000009536743164, 2.0], [0.0, 0.0, 0.0]),
    ("near-double", [1.0, -4.000000007450581, 5.000000022351742, -2.000000014901161], [1.0, 1.0000000074505806, 2.0], [0.0, 0.0, 0.0]),
    ("near-double", [1.0, -4.000000000116415, 5.000000000349246, -2.0000000002328306], [1.0, 1.0000000001164153, 2.0], [0.0, 0.0, 0.0]),
    ("near-double", [1.0, -4.0000000000009095, 5.0000000000027285, -2.000000000001819], [1.0, 1.0000000000009095, 2.0], [0.0, 0.0, 0.0]),
    ("near-double", [1.0, 5.500002861022949, 6.000007152557373, -4.500004291534424], [-3.000002861022949, -3.0, 0.5], [0.0, 0.0, 0.0]),
    ("double", [1.0, -3.0, 0.703125, -0.04296875], [0.125, 0.125, 2.75], [0.0, 0.0, 0.0]),
    ("double", [1.0, -6.625, 4.546875, -0.826171875], [0.375, 0.375, 5.875], [0.0, 0.0, 0.0]),
    ("double", [1.0, -5.0, 8.0, -4.0], [1.0, 2.0, 2.0], [0.0, 0.0, 0.0]),
    ("double", [1.0, -1.0, -9.75, -9.0], [-1.5, -1.5, 4.0], [0.0, 0.0, 0.0]),
    ("triple", [1.0, -6.0, 12.0, -8.0], [2.0, 2.0, 2.0], [0.0, 0.0, 0.0]),
    ("triple", [1.0, -0.375, 0.046875, -0.001953125], [0.125, 0.125, 0.125], [0.0, 0.0, 0.0]),
    ("triple", [1.0, 11.25, 42.1875, 52.734375], [-3.75, -3.75, -3.75], [0.0, 0.0, 0.0]),
    ("triple", [1.0, -3072.0, 3145728.0, -1073741824.0], [1024.0, 1024.0, 1024.0], [0.0, 0.0, 0.0]),
    ("triple", [1.0, -0.0029296875, 2.86102294921875e-06, -9.313225746154785e-10], [0.0009765625, 0.0009765625, 0.0009765625], [0.0, 0.0, 0.0]),
    ("triple", [1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
    ("one-real", [1.0, 1.0, 1.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 1.0]),
    ("one-real", [1.0, -3.0, 2.2500000009313226, -0.5000000018626451], [2.0, 0.5, 0.5], [0.0, -3.0517578125e-05, 3.0517578125e-05]),
    ("one-real", [1.0, -3.0, 2.2500009536743164, -0.5000019073486328], [2.0, 0.5, 0.5], [0.0, -0.0009765625, 0.0009765625]),
    ("one-real", [1.0, -3.0, 2.2509765625, -0.501953125], [2.0, 0.5, 0.5], [0.0, -0.03125, 0.03125]),
    ("one-real", [1.0, -3.0, 2.28125, -0.5625], [2.0, 0.5, 0.5], [0.0, -0.1767766952966369, 0.1767766952966369]),
    ("one-real", [1.0, 7.5, 26.0, 42.0], [-3.5, -2.0, -2.0], [0.0, -2.8284271247461903, 2.8284271247461903]),
    ("one-real", [1.0, 99.5, 6546.0, -3298.0], [0.5, -50.0, -50.0], [0.0, -64.0, 64.0]),
    ("one-real", [1.0, -1003.0, 3004.25, -4250.0], [1000.0, 1.5, 1.5], [0.0, -SQRT_2, SQRT_2]),
    ("one-real", [1.0, 1.0, 0.3282470703125, 0.035186767578125], [-0.25, -0.375, -0.375], [0.0, -0.011048543456039806, 0.011048543456039806]),
    ("one-real", [1.0, 2.0, 2.0, 0.0], [0.0, -1.0, -1.0], [0.0, -1.0, 1.0]),
    ("one-real", [1.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, -1.0, 1.0]),
    ("one-real", [1.0, 0.0, 0.125, -1.0], [0.9583584482351548, -0.4791792241175774, -0.4791792241175774], [0.0, -0.9021298057806152, 0.9021298057806152]),
    ("one-real", [1.0, 0.0, 0.0078125, -1.0], [0.9973958392355421, -0.49869791961777105, -0.49869791961777105], [0.0, -0.8682806833583703, 0.8682806833583703]),
    ("one-real", [1.0, 0.0, 0.0009765625, -1.0], [0.9996744791781682, -0.4998372395890841, -0.4998372395890841], [0.0, -0.8663073131055206, 0.8663073131055206]),
    ("one-real", [1.0, 0.0, 6.103515625e-05, -1.0], [0.9999796549479195, -0.49998982747395976, -0.49998982747395976], [0.0, -0.8660430231163866, 0.8660430231163866]),
    ("one-real", [1.0, 0.0, 9.5367431640625e-07, -1.0], [0.9999996821085612, -0.4999998410542806, -0.4999998410542806], [0.0, -0.8660256790865003, 0.8660256790865003]),
];

/// `(p, q, pw, qw)` of the depressed cubic `t^3 + p t + q`, `t = x + B/3`, of the monic `c`:
/// `p = C - B^2/3` and `q = 2 B^3/27 - B C/3 + D` are sums of terms of the sizes `pw = |C| + B^2/3`
/// and `qw = 2 |B|^3/27 + |B C|/3 + |D|`.
fn depressed(c: [f64; 4]) -> (f64, f64, f64, f64) {
    let (b, cc, d) = (c[1] / c[0], c[2] / c[0], c[3] / c[0]);
    let (cube, bc) = (2.0 * b.powi(3) / 27.0, b * cc / 3.0);
    let (p, q) = (cc - b * b / 3.0, cube - bc + d);
    (
        p,
        q,
        cc.abs() + b * b / 3.0,
        cube.abs() + bc.abs() + d.abs(),
    )
}

/// `(W + Wd) / P + R` for root `j` of the monic polynomial `c` with roots `z` (real and imaginary
/// part), the first-order effect of rounding at two stages: `W = |z|^3 + |c_1| |z|^2 + |c_2| |z| +
/// |c_3|` for the normalisation `b / a`, `c / a`, `d / a`, and `Wd = |t| pw + qw`, `t = z + B/3`,
/// for the formation of `p` and `q`, which `W` misses where the roots lie close together far from
/// the origin; `P = |p'(z)|`, `R` the largest root.
fn sensitivity(c: [f64; 4], z: [(f64, f64); 3], j: usize) -> f64 {
    let modulus = z.map(|z| z.0.hypot(z.1));
    let w = ((modulus[j] + c[1].abs()) * modulus[j] + c[2].abs()) * modulus[j] + c[3].abs();
    let (_, _, pw, qw) = depressed(c);
    let t = (z[j].0 + c[1] / 3.0).hypot(z[j].1);
    let p: f64 = (0..3)
        .filter(|&k| k != j)
        .map(|k| (z[j].0 - z[k].0).hypot(z[j].1 - z[k].1))
        .product();
    (w + t * pw + qw) / p + modulus.into_iter().fold(0.0, f64::max)
}

/// The rounding error of `h - s` of the one-real-root arm as it reaches the root, `0` if the
/// discriminant is not positive; it saturates at the size of the root `v` and of the noise
/// `cbrt(4 u (h + s))`, which a `h - s` that cancels to nothing turns into a cube root.
fn cardano_loss(c: [f64; 4], u: f64) -> f64 {
    let (p, q, ..) = depressed(c);
    let h = q.abs() / 2.0;
    let disc = h * h + p * p * p / 27.0;
    if disc <= 0.0 {
        return 0.0;
    }
    let s = disc.sqrt();
    let (v, noise) = (p.abs() / (3.0 * (h + s).cbrt()), (4.0 * u * (h + s)).cbrt());
    (4.0 * u * (h + s) / (3.0 * v * v)).min(v + noise)
}

/// The rounding error of `disc = q^2/4 + p^3/27` for the monic `c` at unit roundoff `u`, relative to
/// its larger summand `big`, and `big`. `p` and `q` are formed to `u` of their terms' sizes (`pw`,
/// `qw` of `depressed`), so they cancel when the roots lie close together far from the origin; the
/// second order in `dq` matters where `q` is smaller than its own error.
fn disc_noise(c: [f64; 4], u: f64) -> (f64, f64) {
    let (p, q, pw, qw) = depressed(c);
    let (dp, dq) = (u * pw, u * qw);
    let (tq, tp) = (q * q / 4.0, (p.powi(3) / 27.0).abs());
    let big = tq.max(tp).max(f64::MIN_POSITIVE);
    let abs =
        (q.abs() / 2.0 + dq / 4.0) * dq + (p * p / 9.0 + p.abs() * dp / 9.0 + dp * dp / 27.0) * dp;
    ((abs + u * (tq + tp)) / big, big)
}

/// `disc` of the monic cubic with roots `z`: `-Delta / 108`, `Delta` the product of the squared
/// differences.
fn disc_exact(z: [(f64, f64); 3]) -> f64 {
    let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let mut delta = (1.0, 0.0);
    for (i, j) in [(0, 1), (0, 2), (1, 2)] {
        let d = (z[i].0 - z[j].0, z[i].1 - z[j].1);
        delta = mul(delta, mul(d, d));
    }
    -delta.0 / 108.0
}

/// Whether the model of the header holds for a row at `S`: the band is relative to `disc`'s summands,
/// not to `B`, `C`, `D`, so the arm is only certain where `disc` and the band exceed the rounding
/// error of `disc`. Inside it a pair of real roots is answered as one root, and a complex pair as a
/// real double root, which moves the real root too (`a_repeated_root_can_be_dropped`).
fn regular<S: Real>(&(stratum, c, re, im): &Row) -> bool {
    let (noise, big) = disc_noise(c, unit::<S>());
    let z = [0, 1, 2].map(|j| (re[j], im[j]));
    let band = table::<S>().1;
    match stratum {
        "triple" => true,
        "one-real" => disc_exact(z) / big > band + noise,
        _ => noise <= band,
    }
}

/// One row at scalar `S` against the bound of the header, and the number of valid slots its stratum
/// allows. A row with a coefficient that is not an `f32` is skipped at `f32`.
fn check<S: Lane>(&(stratum, c, re, im): &Row) {
    let (u, band) = (unit::<S>(), table::<S>().1);
    if c.iter()
        .any(|&v| S::make(v, [0.0; 2]).lane(0).to_bits() != v.to_bits())
    {
        return;
    }
    let (x, m) = solve::<S>(c);
    if !is_cubic::<S>(c) {
        assert!(m == [false; 3], "{c:?} is below the leading floor");
        return;
    }
    let z = [0, 1, 2].map(|j| (re[j], im[j]));
    let scale = z
        .iter()
        .map(|z| z.0.hypot(z.1))
        .fold(f64::MIN_POSITIVE, f64::max);
    // The cluster rounding spreads a multiple root into; a real root and a pair near it are a triple.
    let cluster = if matches!(stratum, "triple" | "one-real") {
        u.cbrt()
    } else {
        u.sqrt()
    };
    let cap = 4.0 * cluster * scale;
    let loss = if stratum == "one-real" {
        cardano_loss(c, u)
    } else {
        0.0
    };
    let pair = |j: usize| {
        if im[j].abs() <= 2.0 * band.sqrt() * scale {
            im[j].abs()
        } else {
            0.0
        }
    };
    let bound = [0, 1, 2].map(|j| (K * u * sensitivity(c, z, j)).min(cap) + loss + pair(j));
    let dist = |v: f64, j: usize| (v - re[j]).hypot(im[j]);
    let slots: Vec<f64> = x.iter().zip(&m).filter(|p| *p.1).map(|p| *p.0).collect();
    for j in (0..3).filter(|&j| im[j].abs() <= 0.0) {
        let err = slots
            .iter()
            .map(|&v| dist(v, j))
            .fold(f64::INFINITY, f64::min);
        assert!(
            err <= bound[j],
            "{stratum} {c:?}: {x:?}, root {j} missed by {err:e}"
        );
    }
    for &v in &slots {
        assert!(
            (0..3).any(|j| dist(v, j) <= bound[j]),
            "{stratum} {c:?}: {x:?} {m:?}"
        );
    }
    let far = im.iter().all(|&i| i.abs() <= 0.0 || i.abs() >= 0.1 * scale);
    match stratum {
        "distinct" => assert_eq!(slots.len(), 3, "{c:?}"),
        "one-real" if far => assert_eq!(slots.len(), 1, "{c:?}"),
        "near-double" | "double" => assert!(slots.len() >= 2, "{c:?}"),
        _ => assert!(!slots.is_empty(), "{c:?}"),
    }
}

#[test]
fn strata_match_mpmath() {
    for row in MP {
        check::<f64>(row);
        check::<f32>(row);
        check::<D>(row);
    }
}

/// A planted row: the roots are `n / 64` for integers `|n| <= 256`, times `2^e`, so the coefficients
/// are exact binary64 numbers (exact `f32` ones when they fit 24 bits, else `check` skips `f32`) and
/// the stratum and every root are known without a solver. The families are three roots, two a `1/64`
/// or `2/64` apart, a double root, a triple root, and one real root with a pair `s +- i b`,
/// `b >= 1/64`.
fn planted() -> impl Strategy<Value = Row> {
    let n = || -256..=256i64;
    (0..5u8, n(), n(), n(), 1..=256i64, -10..=10i32).prop_map(|(family, r, s, t, b, e)| {
        let (stratum, re, im, [e1, e2, e3]) = if family == 4 {
            // (x - r)(x^2 - 2 s x + s^2 + b^2).
            let (w, q) = (2 * s, s * s + b * b);
            ("one-real", [r, s, s], [0, -b, b], [r + w, q + r * w, r * q])
        } else {
            let roots = match family {
                0 => [r, s, t],
                1 => [r, r + 1 + b % 2, s],
                2 => [r, r, s],
                _ => [r, r, r],
            };
            let [x, y, z] = roots;
            let stratum = match (x == y, x == z, y == z) {
                (true, true, _) => "triple",
                _ if family == 1 => "near-double",
                (false, false, false) => "distinct",
                _ => "double",
            };
            (
                stratum,
                roots,
                [0; 3],
                [x + y + z, x * y + x * z + y * z, x * y * z],
            )
        };
        let k = 2f64.powi(e);
        let real = |n: i64| n as f64 / 64.0 * k;
        let c = [
            1.0,
            -real(e1),
            e2 as f64 / 4096.0 * k * k,
            -(e3 as f64) / 262_144.0 * k * k * k,
        ];
        (stratum, c, re.map(real), im.map(real))
    })
}

/// The fixture, the special cases, the strips of the two thresholds, and 200 seeded rows over a wide
/// range of scales.
fn pool() -> Vec<[f64; 4]> {
    let mut rows: Vec<[f64; 4]> = MP.iter().map(|r| r.1).collect();
    let (inf, nan) = (f64::INFINITY, f64::NAN);
    let strip = |delta: f64| [1.0, 0.0, -3.0, 2.0 * (1.0 + delta).sqrt()];
    rows.extend([
        [-1.0, -0.0, -0.0, -0.0],
        [0.0, 1.0, -3.0, 2.0],
        [1.2e-14, 1.0, 1.0, -1.0],
        [1e-310; 4],
        [nan, 0.0, 0.0, 0.0],
        [1.0, inf, 0.0, 0.0],
        [1.0, 0.0, -inf, 1.0],
        [1.0, 0.0, 0.0, nan],
        strip(0.5e-12),
        strip(0.98e-12),
        strip(2e-12),
        strip(-1e-12),
        [1.0, 0.0, -1e-250, 0.0],
        [1.0, 1e-300, 1e-300, 1e-300],
        [1e300, 1e300, 1e300, 1e300],
        [1.0, -1e13, 1e26, -1e39],
    ]);
    let mut state = 1_u64;
    let mut next = move || {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let z = (state ^ (state >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        let z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        (z ^ (z >> 31)) as f64 / u64::MAX as f64 * 2.0 - 1.0
    };
    rows.extend((0..200).map(|_| [(); 4].map(|()| next() * 2f64.powi((next() * 20.0) as i32))));
    rows
}

/// The value path of `Dual` is the plain value path: roots and masks to the bit.
fn value_path<S: Real<Mask = bool>, T: Real<Mask = bool>>(
    conv: fn(f64) -> S,
    lift: fn(S, usize) -> T,
) {
    for c in pool() {
        let x = c.map(conv);
        let (plain, pm) = solve_cubic(x[0], x[1], x[2], x[3]);
        let y = [0, 1, 2, 3].map(|i| lift(x[i], i));
        let (dual, dm) = solve_cubic(y[0], y[1], y[2], y[3]);
        assert_eq!(pm, dm, "{c:?}");
        for (p, d) in plain.0.iter().zip(dual.0) {
            assert_eq!(p.value_f64().to_bits(), d.value_f64().to_bits(), "{c:?}");
        }
    }
}

#[test]
fn dual_value_is_plain_value() {
    value_path::<f64, Dual<f64, 4>>(|x| x, Dual::variable);
    value_path::<f32, Dual<f32, 4>>(|x| x as f32, Dual::variable);
    value_path::<f64, Dual<Dual<f64, 4>, 1>>(|x| x, |x, i| Dual::variable(Dual::variable(x, i), 0));
}

/// Where the root is simple and the arm regular, `d root / d coefficient` is `-x^k / p'(root)` (the
/// implicit function theorem) to `1e-10` of the largest component (measured worst: `2e-13` in the
/// trigonometric arm, `3e-11` in the one-real-root arm); an invalid slot is a constant. The
/// one-real-root arm differentiates `cbrt(h - s)`, which loses more than the value as `p -> 0`
/// (`4e-7` at `p = 2^-10`, `1e-3` at `2^-14`, for `x^3 + p x - 1`; `h - s` is 0 and the derivative
/// `-inf` for every `p` below `1.3e-5`), so rows whose `cardano_loss` exceeds `1e-12` are left out;
/// and `x^3 - 1` has `p = 0`, where `cbrt` is evaluated at 0 and the derivative is not finite though
/// the root is simple.
#[test]
fn derivative_is_the_implicit_function_derivative() {
    type D4 = Dual<f64, 4>;
    for &(stratum, c, _, im) in MP {
        let x = [0, 1, 2, 3].map(|i| D4::variable(c[i], i));
        let (r, m) = solve_cubic(x[0], x[1], x[2], x[3]);
        let real = im.iter().filter(|i| i.abs() <= 0.0).count();
        let regular = cardano_loss(c, unit::<f64>()) <= 1e-12;
        if !matches!(stratum, "distinct" | "one-real")
            || !regular
            || m.iter().filter(|&&v| v).count() != real
        {
            continue;
        }
        for (k, &valid) in m.iter().enumerate() {
            let (v, dv) = (r.0[k].v, r.0[k].d);
            let fp = (3.0 * c[0] * v + 2.0 * c[1]) * v + c[2];
            let want = [v * v * v, v * v, v, 1.0].map(|t| if valid { -t / fp } else { 0.0 });
            let top = want.iter().fold(f64::MIN_POSITIVE, |m, w| m.max(w.abs()));
            let err = dv
                .iter()
                .zip(want)
                .map(|(d, w)| (d - w).abs())
                .fold(0.0, f64::max)
                / top;
            assert!(err < 1e-10, "{stratum} {c:?} slot {k}: {dv:?} vs {want:?}");
        }
    }
    let x = [0, 1, 2, 3].map(|i| D4::variable([1.0, 0.0, 0.0, -1.0][i], i));
    let (r, m) = solve_cubic(x[0], x[1], x[2], x[3]);
    assert!(m == [true, false, false] && r.0[0].v.to_bits() == 1.0_f64.to_bits());
    assert!(r.0[0].d.iter().all(|d| !d.is_finite()));
}

/// Both arms of every branch run in every lane, so a lane that does not take an arm must still
/// hand it a safe argument: `L2`'s `sqrt` asserts a non-negative argument in debug builds.
#[test]
fn two_lanes_agree_with_scalars() {
    let rows = pool();
    let n = rows.len();
    let pairs = (0..n).flat_map(|i| [(i, (i * 37 + 11) % n), (i, i), (i, n - 1 - i)]);
    for (i, j) in pairs.chain((0..48).flat_map(|i| (0..48).map(move |j| (i, j)))) {
        let (ci, cj) = (rows[i], rows[j]);
        let x = [0, 1, 2, 3].map(|k| L2([ci[k], cj[k]]));
        let (r, m) = solve_cubic(x[0], x[1], x[2], x[3]);
        for (lane, c) in [(0, ci), (1, cj)] {
            let (want, want_mask) = solve::<f64>(c);
            for k in 0..3 {
                assert_eq!(
                    r.0[k].0[lane].to_bits(),
                    want[k].to_bits(),
                    "{c:?} slot {k}"
                );
                assert_eq!(m[k].0[lane], want_mask[k], "{c:?} slot {k}");
            }
        }
    }
}

/// Every bit pattern (NaN, infinities, subnormals, the largest numbers), wide-range finite values
/// and the special values.
fn any_coef() -> impl Strategy<Value = f64> {
    let special = std::vec![
        0.0,
        -0.0,
        1.0,
        -1.0,
        f64::INFINITY,
        f64::NAN,
        f64::MIN_POSITIVE,
        5e-324
    ];
    prop_oneof![
        any::<u64>().prop_map(f64::from_bits),
        (-1.0..1.0f64, -40..40i32).prop_map(|(m, e)| m * 2f64.powi(e)),
        prop::sample::select(special),
    ]
}

proptest! {
    #![proptest_config(cfg())]

    /// A slot that is not valid is `+0`, a valid one finite and never `-0` (checked by `solve`),
    /// and nothing is valid when a coefficient is not finite or the leading one is below half its
    /// floor.
    #[test]
    fn the_mask_is_honest_for_any_bits(c in prop::array::uniform4(any_coef())) {
        let finite = c.iter().all(|x| x.is_finite());
        let masks = [(table::<f64>(), solve::<f64>(c)), (table::<f32>(), solve::<f32>(c)), (table::<f64>(), solve::<D>(c))];
        for ((floor, _), (_, m)) in masks {
            if !finite || c[0].abs() < floor / 2.0 {
                prop_assert!(m == [false; 3], "{c:?}: {m:?}");
            }
        }
    }

    /// Planted roots against the error model of the header, at `f64`, `f32` and `Dual`.
    #[test]
    fn planted_roots_match_the_error_model(row in planted()) {
        if regular::<f64>(&row) {
            check::<f64>(&row);
            check::<D>(&row);
        }
        if regular::<f32>(&row) {
            check::<f32>(&row);
        }
    }
}

/// omnisac d3be7b7's algorithm on random binary64 coefficients, to the bit: `(coefficients, roots,
/// mask)` with the bits of each `f64`, the mask bit `k` being slot `k` (see the header). Pins the
/// operation order, the slot order and `libm`.
#[rustfmt::skip]
const GOLDEN: &[([u64; 4], [u64; 3], u8)] = &[
    // three real roots (mask bit `k` is slot `k`)
    ([0x3FE2_69C3_B502_F4E2, 0x400A_0594_A6FA_BF9D, 0x4005_D238_7D2F_EAFB, 0xBFEE_3B4D_379C_1B1A], [0x3FD0_B7E1_0C66_3914, 0xBFF6_3931_A76A_8750, 0xC012_19BC_FD0C_0596], 7),
    ([0x3FF0_0000_0000_0000, 0x3FF9_4DA8_4713_A0C1, 0xC00C_3E58_C8A8_3FE2, 0xC016_A6F5_B317_F577], [0x3FFE_2937_8CA7_4282, 0xBFFB_A6F5_0895_C2A0, 0xBFFB_CFEA_CB25_20A6], 7),
    ([0x3FE5_922F_1FCF_62B4, 0x3FD0_A338_ACFA_1C24, 0xC000_EFB7_4066_D698, 0x3FED_7C8F_E6C2_404A], [0x3FF4_0501_C97A_EDE1, 0x3FE0_4A7C_D6B2_1AAA, 0xC001_2AED_E8F2_76DD], 7),
    ([0x3E55_11EC_4935_A76D, 0xBE74_4499_F41B_C315, 0xBE79_A805_1B37_F37B, 0x3E5D_B12E_CBCC_AF92], [0x4013_3452_984B_0B34, 0x3FCF_5A5F_0950_9F98, 0xBFF3_2C3F_F187_BD16], 7),
    ([0x3FF0_0000_0000_0000, 0x4015_3B5C_B1F3_F09D, 0x4020_58E5_981E_10AC, 0x4006_EF1F_EE00_A2C7], [0xBFDF_AE96_4D47_B500, 0xC003_3635_DC7B_9F05, 0xC003_4AB0_BDC3_4B97], 7),
    ([0x3FF0_0000_0000_0000, 0xC011_07EB_1FA7_A014, 0x4017_CBA9_2158_49F6, 0xC005_B71F_9325_6E8F], [0x3FF9_89CE_BE71_43E3, 0x3FF9_89C4_4208_6DE3, 0x3FF1_0C19_7E24_CE8A], 7),
    ([0x3FF0_0000_0000_0000, 0xC007_94B2_DBD0_B4F8, 0x3FF6_4777_1106_439E, 0xBFC7_22F6_9B05_1EC4], [0x4003_3011_DCC2_30F7, 0x3FD1_E46F_8158_C906, 0x3FD1_4098_771B_5704], 7),
    ([0x3FE5_E935_4F23_618A, 0x3FED_02DF_13FC_D82B, 0xC00F_EBEE_3CAB_BD36, 0x4002_CB75_BB76_9834], [0x3FF2_6182_5343_289C, 0x3FEC_6ED0_453A_5364, 0xC00A_E41C_AF46_860C], 7),
    ([0x3FE6_3990_9158_C3FE, 0xBFCF_486B_A9EA_BA76, 0xC022_6FB1_E2F9_5022, 0x4025_C19F_F0BE_4D5B], [0x4008_37A4_73B1_150D, 0x3FF4_D199_8662_55FF, 0xC00F_CFC5_1F38_6EE2], 7),
    ([0x3FF2_2653_A220_71F9, 0x4021_5A0D_31A5_8510, 0x4030_2E79_A5C9_868F, 0x4000_FF05_0DC9_983F], [0xBFC2_2955_EA90_6410, 0xC006_7E29_F9ED_4550, 0xC012_C777_4F4C_BBF5], 7),
    ([0x3FF0_0000_0000_0000, 0x401A_16D9_AAC8_D016, 0x402A_2126_72C7_2AF0, 0x401D_9769_57F9_0C41], [0xBFEE_8BD5_5CF1_EFCA, 0xC006_455E_FA23_5532, 0xC006_455F_0431_CF08], 7),
    ([0x3D93_87E2_C901_3D56, 0x3DA2_0558_9F0C_7E79, 0xBD96_66FF_9D99_FB69, 0xBD83_ABE8_E349_2AD6], [0x3FE7_0915_B6ED_8D7B, 0xBFD3_DAEB_88CD_0AD9, 0xC002_0A3E_E2B3_3BFB], 7),
    ([0x3FF0_0000_0000_0000, 0xBFD6_377B_F745_27A4, 0xBFF0_EC2D_F11B_1C2E, 0xBFD4_4DF6_B21E_FA28], [0x3FF5_358F_C773_A329, 0xBFDF_4F61_86C5_DE41, 0xBFDF_4F61_9FC3_86BD], 7),
    ([0x3FF0_0000_0000_0000, 0x0000_0000_0000_0000, 0xC004_6769_7B0B_252E, 0xBFEB_1AD6_182D_6390], [0x3FFB_E197_5A04_5A0E, 0xBFD6_519B_9409_1B32, 0xBFF6_4D30_7502_1341], 7),
    ([0xBED4_FE70_8CF4_5E8C, 0xBE9A_9A67_F56E_4B88, 0x3EE1_1DC7_6B0F_D0D6, 0xBEC2_CF0D_95A0_780E], [0x3FF0_F6C6_3813_F176, 0x3FD2_DB0B_86D2_7F14, 0xBFF6_F1EF_4C30_29A4], 7),
    ([0xC060_C0FE_8EFA_015F, 0x408B_8B6D_1833_0C6F, 0xC083_7850_9895_9C10, 0xC08E_C1C0_DEBF_C496], [0x4015_F046_3541_2FD4, 0x3FFD_3405_E9D8_F81C, 0xBFE7_7935_3A33_F68C], 7),
    ([0x3FF0_0000_0000_0000, 0x0000_0000_0000_0000, 0xC003_E9D1_FC49_42B8, 0x3F9B_BF71_C4FD_A380], [0x3FF9_27E3_F310_7A41, 0x3F86_4BC0_82DC_A4B4, 0xBFF9_547B_7416_338B], 7),
    ([0x3FED_8D3B_EAA6_31B4, 0x4000_0672_86C3_A4F5, 0xC023_296B_9DD3_3B6E, 0xC034_263F_EA2A_47A6], [0x4009_9BDF_DAEE_BA53, 0xC000_75E5_0990_1B7E, 0xC00A_804E_B623_602E], 7),
    ([0x3FE5_6572_0C41_39E0, 0x3FE7_7774_B694_0CD8, 0xBFED_35B2_0482_684E, 0xBFD1_136B_4161_2560], [0x3FEC_E881_D45D_2376, 0xBFD0_2DBB_8CE7_2830, 0xBFFB_F528_B03B_D0A0], 7),
    ([0x3FF0_0000_0000_0000, 0x4000_0508_D822_B31C, 0xBFFF_956B_2D17_1FD0, 0xC010_EC6C_3B8E_7C22], [0x3FF6_EF98_CD57_A558, 0xBFFB_7CA0_D109_0A9D, 0xBFFB_7D09_AC94_00F7], 7),
    // one real root (mask bit `k` is slot `k`)
    ([0x3FE6_8122_FE4D_E24C, 0x4007_2D71_262E_16A2, 0x4017_3E0B_CF22_46BF, 0x4015_D48E_3C5E_9A52], [0xBFFE_BEA7_2E67_BFC7, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    ([0xBFD8_4C06_E4E2_8C48, 0x3FDD_CE81_1E0C_A3A8, 0x3FE5_1325_9D58_2268, 0x3FEA_5AA1_7D4E_2AF6], [0x4002_D679_1C5B_D408, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    ([0x3FC6_B517_C7C1_1C50, 0xBFE6_7619_FC92_666C, 0x3FE9_5D03_F501_6C4C, 0x3FEC_B1AB_B57E_A49C], [0xBFE5_64FC_F1E6_82AE, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    ([0x40F0_3356_85FD_2D73, 0x4100_3C8F_D40E_0913, 0x4128_2385_60D4_6FA7, 0xC0DD_8E46_2919_D716], [0x3FA3_76BB_2DD0_9360, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    ([0x3FF0_0000_0000_0000, 0x0000_0000_0000_0000, 0xC003_DF7D_8951_6414, 0xC003_303E_6DB7_B958], [0x3FFE_E319_A57D_A823, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    ([0x3FF0_0000_0000_0000, 0x0000_0000_0000_0000, 0xC001_2BF8_8F3E_3639, 0x4005_02C1_EC18_31AB], [0xBFFE_1CB8_B08F_23A4, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    ([0x3FD1_27F1_54DF_834C, 0xBFCD_0D71_CC27_D078, 0xBFD6_5DA7_787B_CE4C, 0xBFCC_D1AD_9D22_51A0], [0x3FFD_165B_4528_EA3C, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    ([0x3FF7_B278_DD23_8C1E, 0xBFEC_2CE2_EEEF_BCF7, 0xBFF1_8ECB_45C1_18C6, 0x3FEC_080A_8E76_5C00], [0xBFED_57BC_34FC_C64F, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    ([0x4007_D9D0_C4DC_777C, 0xBFD5_3065_E8A6_E94A, 0x3FE5_1B49_ADAD_9DE4, 0xC000_EEB7_6440_0552], [0x3FEB_0635_4BAD_1511, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    ([0x3FE5_CF01_C71C_D7A6, 0xBFD8_73A1_9701_0328, 0x3FCF_DF92_9B67_93D8, 0xBFE2_D340_9BAC_8002], [0x3FF0_67BE_39EA_85A4, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    ([0x3FF2_2335_1EAC_F265, 0x4009_3DE0_B629_1683, 0x3FF9_713F_2799_1C4F, 0x3FF0_A02C_6B14_1656], [0xC002_D26D_65FF_AE78, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    ([0x3FF0_0000_0000_0000, 0x0000_0000_0000_0000, 0x3FFC_E02F_3056_86B4, 0x4000_6E4D_1DB2_5D44], [0xBFEA_6DB7_C011_F3B9, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    // a triple root (mask bit `k` is slot `k`)
    ([0x3FF0_0000_0000_0000, 0xBCC5_DFCB_9B20_E65E, 0x3983_EFCB_2891_055F, 0xB628_3A4C_DCD6_7DF2], [0x3CAD_2A64_CED6_887D, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    ([0x3FF0_0000_0000_0000, 0xBFF4_951F_141B_3D46, 0x3FE1_A6D2_58E4_A1AF, 0xBFB4_2F25_FF71_8A3C], [0x3FDB_717E_C579_A708, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    ([0x3FF0_0000_0000_0000, 0xBB4F_2C22_D3B7_26FE, 0x3694_3E80_3C14_A70E, 0xB1C1_8787_49A9_A514], [0x3B34_C817_37CF_6F54, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
    ([0x3FF0_0000_0000_0000, 0xBF01_256B_8C89_A433, 0x3DF8_7FD0_13D7_32B4, 0xBCD7_564F_6749_ED12], [0x3EE6_DC8F_660C_DAEF, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000], 1),
];

#[test]
fn bits_are_omnisacs() {
    for &(c, r, m) in GOLDEN {
        let c = c.map(f64::from_bits);
        let (got, mask) = solve::<f64>(c);
        assert!(
            got.map(f64::to_bits) == r && mask == [m & 1 != 0, m & 2 != 0, m & 4 != 0],
            "{c:?}: {got:?} {mask:?}"
        );
    }
}

/// Slot `k` of the trigonometric arm is the `k`-th phase `2 r cos(acos(..)/3 - 2 pi k / 3)`, so three
/// roots come largest first, up to rounding; omnisac's adapter pushes them in slot order
/// (`docs/PHASE2.md` §9) and its callers read them so. The other tests sort.
fn slot_order<S: Lane>() {
    for (c, want) in [
        ([1.0, -6.0, 11.0, -6.0], [3.0, 2.0, 1.0]),
        ([1.0, -2.0, -5.0, 6.0], [3.0, 1.0, -2.0]),
        ([-2.0, 12.0, -22.0, 12.0], [3.0, 2.0, 1.0]),
    ] {
        let (r, m) = solve::<S>(c);
        let close = |(x, w): (&f64, &f64)| (x - w).abs() <= 64.0 * unit::<S>() * w.abs().max(1.0);
        assert!(
            m == [true; 3] && r.iter().zip(&want).all(close),
            "{c:?}: {r:?}"
        );
    }
}

#[test]
fn slots_are_in_the_order_of_the_phases() {
    slot_order::<f64>();
    slot_order::<f32>();
    slot_order::<D>();
}

/// Where `p^3` and `q^2` of the depressed cubic underflow (roots below about `4e-8` at `f32`,
/// `2e-54` at `f64`) `disc` reads 0. `2 p r` underflowing first leaves no root (the `acos` argument
/// is NaN); with `p >= 0` the triple-root arm answers `cbrt(-q)`; with `p < 0` the trigonometric arm
/// reports three valid slots for a cubic with one real root, none of them a root (`# Domain`).
fn underflow<S: Lane>() {
    let rows = match S::PRECISION {
        Precision::F64 => [
            [1.0, 0.0, -1e-250, 0.0],
            [1.0, 0.0, 1e-110, 1e-200],
            [1.0, 0.0, -1e-110, 1e-165],
        ],
        Precision::F32 => [
            [1.0, 0.0, -1e-32, 0.0],
            [1.0, 0.0, 1e-16, 1e-30],
            [1.0, 0.0, -1e-16, 1e-24],
        ],
    };
    assert!(solve::<S>(rows[0]).1 == [false; 3]);
    let (r, m) = solve::<S>(rows[1]);
    assert!(m == [true, false, false] && (r[0] + rows[1][3].cbrt()).abs() <= 1e-3 * r[0].abs());
    let c = rows[2].map(|x| S::make(x, [0.0; 2]).lane(0));
    let (r, m) = solve::<S>(c);
    let f = |x: f64| (x * x + c[2]) * x + c[3];
    assert!(
        m == [true; 3] && r.iter().all(|&x| f(x).abs() > 0.5 * c[3]),
        "{r:?}"
    );
}

#[test]
fn underflow_is_pinned_at_both_precisions() {
    underflow::<f64>();
    underflow::<f32>();
    underflow::<D>();
}

/// Where `disc` cancels the model of the header claims nothing (`regular`): `(x - 1.09375)^2 (x -
/// 0.921875)` has `disc = 0`, formed with an error that may exceed the band at either precision. At
/// `f32` it does: `disc` reads as positive, the arm has one root and the double root is dropped; at
/// `f64` all three are found, this row's error being below its bound.
#[test]
fn a_repeated_root_can_be_dropped() {
    let c = [1.0, -3.109375, 3.212890625, -1.1028289794921875];
    let row = ("double", c, [1.09375, 1.09375, 0.921875], [0.0; 3]);
    assert!(!regular::<f32>(&row));
    assert_eq!(valid::<f64>(c).len(), 3);
    let got = valid::<f32>(c);
    assert!(
        got.len() == 1 && (got[0] - 0.921875).abs() < 1e-4,
        "{got:?}"
    );
}

/// `pi` is `atan2(+0, -1)`, correctly rounded at each precision and a constant for `Dual`; `acos`
/// by `atan2` is within 4 `u` of `libm::acos` and `acosf` (2 ulp, on 4 million random points too)
/// on a grid, at the ends and near them (`1 - 2^-k`, where `sqrt(1 - x^2)` would lose everything).
#[test]
fn pi_and_acos_are_within_their_ulps() {
    assert_eq!(pi::<f64>().to_bits(), PI.to_bits());
    assert_eq!(pi::<f32>().to_bits(), core::f32::consts::PI.to_bits());
    assert!(pi::<Dual<f64, 2>>().d.iter().all(|d| d.abs() <= 0.0));
    let grid = (0..=4096).map(|i| -1.0 + f64::from(i) / 2048.0);
    let ends = (1..53).flat_map(|k| [1.0 - 2f64.powi(-k), 2f64.powi(-k - 1) - 1.0, 2f64.powi(-k)]);
    for x in grid.chain(ends) {
        let (a, w) = (acos(x), libm::acos(x));
        assert!(
            (a - w).abs() <= 4.0 * unit::<f64>() * w.abs(),
            "f64 at {x:e}"
        );
        let (a, w) = (f64::from(acos(x as f32)), f64::from(libm::acosf(x as f32)));
        assert!(
            (a - w).abs() <= 4.0 * unit::<f32>() * w.abs(),
            "f32 at {x:e}"
        );
    }
    assert!(acos(1.0_f64).to_bits() == 0 && acos(-1.0_f64).to_bits() == PI.to_bits());
    assert!(acos(f64::NAN).is_nan());
}
