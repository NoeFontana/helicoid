//! Tests of the kernel (`docs/PHASE3.md` §3, §9). The corpus is `xtask`'s to read: `cargo xtask
//! thresholds` scores these arms and groups over it. Here the groups are checked against the arms
//! and the arms against each other, at both precisions, and the generated switches against
//! `conformance/sweeps/thresholds.csv`, the sweep that chose them.

use core::cell::Cell;
use core::marker::PhantomData;
use core::ops::{Add, Div, Mul, Neg, Sub};
use std::string::{String, ToString};
use std::vec::Vec;
use std::{format, thread_local, vec};

use helicoid_linalg::{Dual, Mask, Precision, Real};

use super::generated::{
    A_F32, A_F64, B_F32, B_F64, COS_HALF_F32, COS_HALF_F64, C_F32, C_F64, D_F32, D_F64, E_F32,
    E_F64, K_F32, K_F64, R_F32, R_F64,
};
use super::kernel::{
    exact_a, exact_b, exact_c, exact_cos_half, exact_d, exact_e, exact_k, exact_r, series_a,
    series_b, series_c, series_cos_half, series_d, series_e, series_k, series_r,
};
use super::{exp_coeffs, jr_coeffs, jr_inv_coeff, log_ratio, q_coeffs, Switch};

/// The order of the coefficients everywhere below, and of `coeff_series.jsonl`.
const NAMES: [&str; 8] = ["k", "a", "b", "c", "d", "e", "cos_half", "r"];

type D<S> = Dual<S, 1>;

/// `(exact, series)` of `NAMES[i]` at `z` and `terms` series terms, `r` at `w = 1`.
fn arms<S: Real>(i: usize, z: S, terms: usize) -> (S, S) {
    let one = S::one();
    match i {
        0 => (exact_k(z), series_k(z, terms)),
        1 => (exact_a(z), series_a(z, terms)),
        2 => (exact_b(z), series_b(z, terms)),
        3 => (exact_c(z), series_c(z, terms)),
        4 => (exact_d(z), series_d(z, terms)),
        5 => (exact_e(z), series_e(z, terms)),
        6 => (exact_cos_half(z), series_cos_half(z, terms)),
        _ => (exact_r(z, one), series_r(z, one, terms)),
    }
}

/// Every group at `(z, w)`: `NAMES`, then `b` again, from `q_coeffs`.
fn groups<S: Real>(z: S, w: S) -> [S; 9] {
    let (k, cos_half) = exp_coeffs(z);
    let (a, b) = jr_coeffs(z);
    let (b_q, d, e) = q_coeffs(z);
    [
        k,
        a,
        b,
        jr_inv_coeff(z),
        d,
        e,
        cos_half,
        log_ratio(z, w),
        b_q,
    ]
}

/// A generated constant's `(switch, bits, terms)`.
fn f64_of<const M: usize>(s: &Switch<f64, M>) -> (f64, u64, usize) {
    (s.below, s.below.to_bits(), M)
}

fn f32_of<const M: usize>(s: &Switch<f32, M>) -> (f64, u64, usize) {
    (f64::from(s.below), u64::from(s.below.to_bits()), M)
}

/// The generated `(switch, bits, terms)` of every coefficient at `S`'s precision.
fn chosen<S: Real>() -> [(f64, u64, usize); 8] {
    match S::PRECISION {
        Precision::F64 => [
            f64_of(&K_F64),
            f64_of(&A_F64),
            f64_of(&B_F64),
            f64_of(&C_F64),
            f64_of(&D_F64),
            f64_of(&E_F64),
            f64_of(&COS_HALF_F64),
            f64_of(&R_F64),
        ],
        Precision::F32 => [
            f32_of(&K_F32),
            f32_of(&A_F32),
            f32_of(&B_F32),
            f32_of(&C_F32),
            f32_of(&D_F32),
            f32_of(&E_F32),
            f32_of(&COS_HALF_F32),
            f32_of(&R_F32),
        ],
    }
}

/// `x` rounded to a value `S` holds.
fn at<S: Real>(x: f64) -> S {
    match S::PRECISION {
        Precision::F64 => S::lit(x),
        Precision::F32 => S::lit(f64::from(x as f32)),
    }
}

/// The values of `z` that matter, each one `S` holds: zero, subnormal, tiny, log-uniform over
/// `1e-12..40`, and every switch with its two neighbours.
fn samples<S: Real>() -> Vec<S> {
    let mut state = 7u64;
    let mut next = || {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let z = (state ^ (state >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        let z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut zs = vec![0.0, 5e-324, 1e-310, 1e-40, 1e-30, 1.0, 6.0, 39.0];
    zs.extend((0..600).map(|_| 10f64.powf(-12.0 + 13.6 * next())));
    for (below, ..) in chosen::<S>() {
        let near = match S::PRECISION {
            Precision::F64 => [below.next_down(), below.next_up()],
            Precision::F32 => [(below as f32).next_down(), (below as f32).next_up()].map(f64::from),
        };
        zs.extend([below, near[0], near[1]]);
    }
    zs.into_iter().map(at).collect()
}

fn bits<S: Real>(x: S) -> u64 {
    x.value_f64().to_bits()
}

/// The generated switches: what the sweep recorded for `NAMES[i]` at `precision` in
/// `conformance/sweeps/thresholds.csv`.
struct Recorded {
    terms: usize,
    bits: u64,
    /// The exact arm's and the series arm's `(value, derivative)` maxima, in units of `u`.
    exact: (f64, f64),
    series: (f64, f64),
}

fn recorded(name: &str, precision: &str) -> Result<Recorded, String> {
    const SWEEP: &str = include_str!("../../../../conformance/sweeps/thresholds.csv");
    let mut lines = SWEEP.lines();
    let header: Vec<&str> = lines.next().ok_or("empty sweep")?.split(',').collect();
    let row: Vec<&str> = SWEEP
        .lines()
        .map(|l| l.split(',').collect::<Vec<_>>())
        .find(|c| c[0] == name && c[1] == precision)
        .ok_or_else(|| format!("no {name} {precision} row"))?;
    let cell = |column: &str| -> Result<&str, String> {
        let i = header.iter().position(|h| *h == column);
        i.and_then(|i| row.get(i).copied())
            .ok_or_else(|| format!("no {column}"))
    };
    let number = |column: &str| -> Result<f64, String> {
        cell(column)?
            .parse::<f64>()
            .map_err(|e| format!("{column}: {e}"))
    };
    let bits = cell("switch_bits")?.trim_start_matches("0x");
    Ok(Recorded {
        terms: cell("terms")?.parse().map_err(|e| format!("terms: {e}"))?,
        bits: u64::from_str_radix(bits, 16).map_err(|e| e.to_string())?,
        exact: (
            number("at_switch_exact_value_u")?,
            number("at_switch_exact_deriv_u")?,
        ),
        series: (
            number("at_switch_series_value_u")?,
            number("at_switch_series_deriv_u")?,
        ),
    })
}

#[test]
fn a_switch_takes_the_first_terms_of_the_swept_series() {
    let s = Switch::<f64, 2>::first(7.0, &[1.0, 2.0, 3.0]);
    assert_eq!(
        (s.below.to_bits(), s.series.map(f64::to_bits)),
        (7f64.to_bits(), [1f64, 2.0].map(f64::to_bits))
    );
}

#[test]
fn the_generated_switches_are_the_sweeps_choice() -> Result<(), String> {
    for (precision, want) in [("f64", chosen::<f64>()), ("f32", chosen::<f32>())] {
        for (name, (_, bits, terms)) in NAMES.iter().zip(want) {
            let r = recorded(name, precision)?;
            assert_eq!((r.bits, r.terms), (bits, terms), "{name} {precision}");
        }
    }
    Ok(())
}

/// Which coefficient of `NAMES` each output of `groups` is.
const OF: [usize; 9] = [0, 1, 2, 3, 4, 5, 6, 7, 2];

/// A group member is its series arm below its own switch and its exact arm elsewhere, bit for bit,
/// value and derivative.
fn groups_are_their_members_own_arm<S: Real + Into<f64>>() {
    let table = chosen::<S>();
    let piece = |i: usize, z: D<S>| {
        let (below, _, terms) = table[i];
        let (exact, series) = arms(i, z, terms);
        D::select(z.lt(D::lit(below)), series, exact)
    };
    for z in samples::<S>() {
        let z = D::variable(z, 0);
        for (got, i) in groups(z, D::one()).iter().zip(OF) {
            let want = piece(i, z);
            let bits = |x: D<S>| (bits(x.v), bits(x.d[0]));
            assert_eq!(
                bits(*got),
                bits(want),
                "{} at {}",
                NAMES[i],
                z.v.value_f64()
            );
        }
    }
}

#[test]
fn groups_are_their_members_own_arm_f64() {
    groups_are_their_members_own_arm::<f64>();
}

#[test]
fn groups_are_their_members_own_arm_f32() {
    groups_are_their_members_own_arm::<f32>();
}

/// The value path of every group under `Dual` is the plain value, bit for bit (`0003` item 5).
fn dual_value_is_plain_value<S: Real + Into<f64>>() {
    for z in samples::<S>() {
        let plain = groups(z, S::one());
        let dual = groups(D::variable(z, 0), D::one());
        for (j, (p, d)) in plain.into_iter().zip(dual).enumerate() {
            assert_eq!(bits(p), bits(d.v), "{} at {}", NAMES[OF[j]], z.value_f64());
        }
        // `r` at other `w`, of either sign: the mask on `w` is on the value part alone.
        for w in [10.0, 0.5, 1e-3, 0.0, -1.0].map(at::<S>) {
            let (r, d) = (
                log_ratio(z, w),
                log_ratio(D::variable(z, 0), D::constant(w)),
            );
            assert_eq!(
                bits(r),
                bits(d.v),
                "r at ({}, {})",
                z.value_f64(),
                w.value_f64()
            );
        }
    }
}

#[test]
fn dual_value_is_plain_value_f64() {
    dual_value_is_plain_value::<f64>();
}

#[test]
fn dual_value_is_plain_value_f32() {
    dual_value_is_plain_value::<f32>();
}

/// `r` takes its series iff `w > 0` and `n²/w² < switch`, else the exact arm, at every sign of `w`.
fn log_ratio_takes_the_series_iff_w_is_positive_and_s_is_below_the_switch<
    S: Real<Mask = bool> + Into<f64>,
>() {
    let (below, _, terms) = chosen::<S>()[7];
    for &w in &[10.0, 1.0, 0.5, 0.05, 1e-3, 0.0, -1e-3, -1.0] {
        for &n2 in &[0.0, 1e-30, 1e-8, 1e-4, 0.01, 0.5, 4.0] {
            // `q = 0` is not a rotation, and `0/0` has no one NaN to compare.
            if n2 <= 0.0 && w <= 0.0 {
                continue;
            }
            let (n2, w) = (D::variable(at::<S>(n2), 0), D::constant(at::<S>(w)));
            let small = S::zero()
                .lt(w.v)
                .and((n2.v / (w.v * w.v)).lt(S::lit(below)));
            let want = match small {
                true => series_r(n2, w, terms),
                false => exact_r(n2, w),
            };
            let got = log_ratio(n2, w);
            let at = (n2.v.value_f64(), w.v.value_f64());
            assert_eq!(
                (bits(got.v), bits(got.d[0])),
                (bits(want.v), bits(want.d[0])),
                "{at:?}"
            );
        }
    }
    // `w` so small that `w²` underflows (`1e-200` at `f64`, `1e-30` at `f32`): `s` is infinite, the
    // mask false, and the exact arm answers, finite.
    let tiny = match S::PRECISION {
        Precision::F64 => 1e-200,
        Precision::F32 => 1e-30,
    };
    for n2 in [1.0, 1e-8] {
        let (n2, w) = (D::variable(at::<S>(n2), 0), D::constant(at::<S>(tiny)));
        let (got, want) = (log_ratio(n2, w), exact_r(n2, w));
        assert_eq!(
            (bits(got.v), bits(got.d[0])),
            (bits(want.v), bits(want.d[0])),
            "n2 = {}, w = {tiny:e}",
            n2.v.value_f64()
        );
        assert!(got.v.value_f64().is_finite());
    }
}

#[test]
fn log_ratio_takes_the_series_iff_w_is_positive_and_s_is_below_the_switch_f64() {
    log_ratio_takes_the_series_iff_w_is_positive_and_s_is_below_the_switch::<f64>();
}

#[test]
fn log_ratio_takes_the_series_iff_w_is_positive_and_s_is_below_the_switch_f32() {
    log_ratio_takes_the_series_iff_w_is_positive_and_s_is_below_the_switch::<f32>();
}

/// The two arms are one function: they agree to far below the accuracy of the exact arm where both
/// are good, in value and derivative (a wrong formula or a wrong series does not). `tol` is ten
/// times the largest disagreement measured, the exact arm's own error: 4.9e-11 at `f64`, 1.2e-4 at
/// `f32`. `r`'s series converges for `n²/w² < 1` only, and the others' for every `z` the eight
/// terms reach. A dropped top term is seen at `f64` alone: at `f32` the eighth term changes no bit
/// of any series arm for `z <= 1`.
fn the_arms_agree_where_both_are_good<S: Real + Into<f64>>(zs: [f64; 2], tol: f64) {
    for (i, name) in NAMES.iter().enumerate() {
        for z in if i == 7 { [0.005, 0.02] } else { zs } {
            let (exact, series) = arms(i, D::variable(at::<S>(z), 0), 8);
            for (e, s) in [(exact.v, series.v), (exact.d[0], series.d[0])] {
                let (e, s) = (e.value_f64(), s.value_f64());
                assert!(((e - s) / s).abs() <= tol, "{name} at {z}: {e} vs {s}");
            }
        }
    }
}

#[test]
fn the_arms_agree_where_both_are_good_f64() {
    the_arms_agree_where_both_are_good::<f64>([0.3, 1.0], 5e-10);
}

#[test]
fn the_arms_agree_where_both_are_good_f32() {
    the_arms_agree_where_both_are_good::<f32>([1.0, 2.0], 1.2e-3);
}

/// The bits of `x` at its own precision.
fn pattern<S: Real>(x: S) -> u64 {
    match S::PRECISION {
        Precision::F64 => x.value_f64().to_bits(),
        Precision::F32 => u64::from((x.value_f64() as f32).to_bits()),
    }
}

/// An exact arm and its name.
type Named<S> = (&'static str, fn(S) -> S);

/// The operand order of the exact arms (`docs/maths/coefficients.md` CO.6) is pinned here, where
/// `cargo nextest run -p helicoid` reaches it: `b`, `c`, `d`, `e` at `z = 0.5` and `2`, whose `θ̂²`
/// differs from `z` in the last bit, are the seeded kernels' bits (checked against them once, and at
/// every corpus record by `xtask`'s `the_shipped_arms_are_the_seeded_arms_bit_for_bit`).
fn the_exact_arms_keep_their_operand_order<S: Real>(want: [[u64; 2]; 4]) {
    let arms: [Named<S>; 4] = [
        ("b", exact_b),
        ("c", exact_c),
        ("d", exact_d),
        ("e", exact_e),
    ];
    for ((name, arm), want) in arms.into_iter().zip(want) {
        for (z, want) in [0.5, 2.0].into_iter().zip(want) {
            assert_eq!(pattern(arm(at::<S>(z))), want, "{name} at {z}");
        }
    }
}

#[test]
fn the_exact_arms_keep_their_operand_order_f64() {
    the_exact_arms_keep_their_operand_order::<f64>([
        [0x3fc4_ce6a_066c_004e, 0x3fc3_4c7f_360a_49ab],
        [0x3fb5_8364_a25d_38c0, 0x3fb6_1480_c86b_8714],
        [0x3fa4_fb1e_8b34_487d, 0x3fa3_f5f6_8684_70ba],
        [0x3f80_aa1e_4aa8_9f22, 0x3f7f_03c7_2b08_aaec],
    ]);
}

#[test]
fn the_exact_arms_keep_their_operand_order_f32() {
    the_exact_arms_keep_their_operand_order::<f32>([
        [0x3e26_7352, 0x3e1a_63fa],
        [0x3dac_1b40, 0x3db0_a408],
        [0x3d27_d901, 0x3d1f_afb5],
        [0x3c05_51ed, 0x3bf8_1e37],
    ]);
}

/// At `z = 0` the exact arms are `0/0` and the series arm answers: `k`, `a`, `b`, `c`, `d`, `e`,
/// `cos(θ/2)` and `r` are `1/2, 1/2, 1/6, 1/12, 1/24, 1/120, 1, 2` there (`docs/maths/coefficients.md`
/// CO.2, CO.3), each the correctly rounded quotient.
fn the_series_answers_at_zero<S: Real>() {
    let inv = |n: f64| S::one() / S::lit(n);
    let want = [
        inv(2.0),
        inv(2.0),
        inv(6.0),
        inv(12.0),
        inv(24.0),
        inv(120.0),
        S::one(),
        S::lit(2.0),
    ];
    let got = groups(S::zero(), S::one());
    for (i, w) in want.into_iter().chain([inv(6.0)]).enumerate() {
        assert_eq!(bits(got[i]), bits(w), "{}", NAMES[OF[i]]);
    }
}

#[test]
fn the_series_answers_at_zero_f64() {
    the_series_answers_at_zero::<f64>();
}

#[test]
fn the_series_answers_at_zero_f32() {
    the_series_answers_at_zero::<f32>();
}

/// `docs/maths/coefficients.md` CO.12: at a switch the arms differ by at most the sum of their
/// errors there, in the value and in the derivative through `Dual`. The errors are those
/// `conformance/sweeps/thresholds.csv` records for the chosen switch: each arm's larger error at the
/// two corpus records that bracket it, a sample and not a bound over the interval between them, and
/// not `NUMERICS.md` §4's wording, that the jump is at most the recorded error of the coefficient,
/// an arm's alone (0015 (draft) NU.4). The jump is taken relative to the exact arm's value and the
/// errors to the true one, so the plain sum can be exceeded, by `E_x u` of itself, and the bound is
/// `(E_x + E_s)/(1 - E_x u)`, `E` in units of `u`, and a few roundings of this test's own. It is
/// tight, not slack: the jump is the sum to three decimals for `a` at `f64` in the derivative (47.234
/// against 47.238 `u`) and exceeds the plain sum for `d` at `f32` (638.380 against 638.360 `u`).
fn branch_continuity<S: Real + Into<f64>>(precision: &str) -> Result<(), String> {
    let u = match S::PRECISION {
        Precision::F64 => 2f64.powi(-53),
        Precision::F32 => 2f64.powi(-24),
    };
    let bound = |ex: f64, es: f64| (ex + es) / (1.0 - ex * u) * (1.0 + 8.0 * f64::EPSILON);
    for (i, (below, _, terms)) in chosen::<S>().into_iter().enumerate() {
        let r = recorded(NAMES[i], precision)?;
        let (exact, series) = arms(i, D::variable(at::<S>(below), 0), terms);
        let jump = |e: S, s: S| {
            let (e, s) = (e.value_f64(), s.value_f64());
            ((e - s) / e).abs() / u
        };
        let (value, deriv) = (jump(exact.v, series.v), jump(exact.d[0], series.d[0]));
        let (max_value, max_deriv) = (bound(r.exact.0, r.series.0), bound(r.exact.1, r.series.1));
        let name = NAMES[i];
        assert!(
            value <= max_value,
            "{name} {precision}: value jump {value} u > {max_value} u"
        );
        assert!(
            deriv <= max_deriv,
            "{name} {precision}: derivative jump {deriv} u > {max_deriv} u"
        );
    }
    Ok(())
}

#[test]
fn branch_continuity_f64() -> Result<(), String> {
    branch_continuity::<f64>("f64")
}

#[test]
fn branch_continuity_f32() -> Result<(), String> {
    branch_continuity::<f32>("f32")
}

/// What a scalar did since the last [`counted`]: the operations with a non-finite result and the
/// calls of each transcendental.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Counts {
    nonfinite: usize,
    sqrt: usize,
    sin_cos: usize,
    atan2: usize,
}

thread_local!(static COUNTS: Cell<Counts> = const {
    Cell::new(Counts { nonfinite: 0, sqrt: 0, sin_cos: 0, atan2: 0 })
});

fn tally(f: impl FnOnce(&mut Counts)) {
    let mut c = COUNTS.get();
    f(&mut c);
    COUNTS.set(c);
}

/// `f`'s result and what it did.
fn counted<T>(f: impl FnOnce() -> T) -> (T, Counts) {
    COUNTS.set(Counts::default());
    let out = f();
    (out, COUNTS.get())
}

/// An `f64` that counts, at the precision `F32` names (its tables, not its arithmetic), under the
/// mask `M`.
#[derive(Clone, Copy)]
struct Lane<const F32: bool, M>(f64, PhantomData<M>);

/// A lane mask: both arms are evaluated and blended (`docs/PHASE2.md` §2).
type Wide<const F32: bool> = Lane<F32, Both>;
/// A scalar mask: one arm, as `f64` and `f32` do.
type Narrow<const F32: bool> = Lane<F32, bool>;

impl<const F32: bool, M> Lane<F32, M> {
    fn new(x: f64) -> Self {
        Lane(x, PhantomData)
    }

    /// `x`, counted when it is not finite.
    fn note(x: f64) -> Self {
        if !x.is_finite() {
            tally(|c| c.nonfinite += 1);
        }
        Self::new(x)
    }
}

#[derive(Clone, Copy)]
struct Both(bool);

impl From<bool> for Both {
    fn from(b: bool) -> Self {
        Both(b)
    }
}

impl Mask for Both {
    fn and(self, o: Self) -> Self {
        Both(self.0 & o.0)
    }
    fn or(self, o: Self) -> Self {
        Both(self.0 | o.0)
    }
    fn not(self) -> Self {
        Both(!self.0)
    }
    fn all(self) -> bool {
        self.0
    }
    fn any(self) -> bool {
        self.0
    }
    fn decide<T>(
        self,
        t: impl FnOnce() -> T,
        f: impl FnOnce() -> T,
        blend: impl FnOnce(Self, T, T) -> T,
    ) -> T {
        let (on, off) = (t(), f());
        blend(self, on, off)
    }
}

macro_rules! lane_op {
    ($($tr:ident $f:ident $op:tt;)*) => {$(
        impl<const F32: bool, M: Mask> $tr for Lane<F32, M> {
            type Output = Self;
            fn $f(self, o: Self) -> Self {
                Self::note(self.0 $op o.0)
            }
        }
    )*};
}

lane_op! { Add add +; Sub sub -; Mul mul *; Div div /; }

impl<const F32: bool, M: Mask> Neg for Lane<F32, M> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.0)
    }
}

impl<const F32: bool, M: Mask + From<bool> + 'static> Real for Lane<F32, M> {
    type Mask = M;
    const PRECISION: Precision = if F32 { Precision::F32 } else { Precision::F64 };
    fn lit(x: f64) -> Self {
        Self::new(x)
    }
    fn zero() -> Self {
        Self::new(0.0)
    }
    fn one() -> Self {
        Self::new(1.0)
    }
    fn lt(self, rhs: Self) -> M {
        M::from(self.0 < rhs.0)
    }
    fn le(self, rhs: Self) -> M {
        M::from(self.0 <= rhs.0)
    }
    fn select(m: M, t: Self, f: Self) -> Self {
        if m.all() {
            t
        } else {
            f
        }
    }
    fn sqrt(self) -> Self {
        tally(|c| c.sqrt += 1);
        Self::note(libm::sqrt(self.0))
    }
    fn sin_cos(self) -> (Self, Self) {
        tally(|c| c.sin_cos += 1);
        let (s, c) = libm::sincos(self.0);
        (Self::note(s), Self::note(c))
    }
    fn atan2(self, x: Self) -> Self {
        tally(|c| c.atan2 += 1);
        Self::note(libm::atan2(self.0, x.0))
    }
    fn abs(self) -> Self {
        Self::new(libm::fabs(self.0))
    }
    fn copysign(self, s: Self) -> Self {
        Self::new(libm::copysign(self.0, s.0))
    }
    fn value_f64(self) -> f64 {
        self.0
    }
}

/// The generated switches of `F32`'s precision.
fn table<const F32: bool>() -> [(f64, u64, usize); 8] {
    if F32 {
        chosen::<f32>()
    } else {
        chosen::<f64>()
    }
}

/// A lane evaluates every arm of every group: none may see a non-finite operation, in value or
/// derivative, at `θ² = 0`, subnormal, tiny, at each switch or above (`0003` item 3).
fn a_lane_that_evaluates_both_arms_sees_no_non_finite_operation<const F32: bool>() {
    type L<const F32: bool> = D<Wide<F32>>;
    let one = L::<F32>::one();
    let mut zs: Vec<f64> = vec![
        0.0, 5e-324, 1e-310, 1e-40, 1e-30, 1e-16, 1e-8, 0.3, 1.0, 6.0, 39.0,
    ];
    zs.extend(
        table::<F32>()
            .iter()
            .flat_map(|&(below, ..)| [below * 0.999, below, below * 1.001]),
    );
    for z in zs {
        let z = L::<F32>::variable(Lane::new(z), 0);
        let (out, c) = counted(|| groups(z, one));
        assert_eq!(c.nonfinite, 0, "a non-finite operation at z = {}", z.v.0);
        assert!(
            out.iter()
                .all(|c| c.v.0.is_finite() && c.d[0].0.is_finite()),
            "z = {}",
            z.v.0
        );
    }
    // `r` at every sign of `w`, `n²` from zero to subnormal, `w = 0` included (`CO.16(d)`), and `w`
    // small enough that `w²` is a small normal number (`1e-100`): the divisor of `s` is `w²` there.
    for (n2, w) in [
        (0.0, 1.0),
        (5e-324, 1.0),
        (1.0, 0.0),
        (1e-4, 0.5),
        (1e-6, -1.0),
        (1e-30, 1e-3),
        (1.0, 1e-100),
    ] {
        let (r, c) = counted(|| {
            log_ratio(
                L::<F32>::variable(Lane::new(n2), 0),
                L::<F32>::constant(Lane::new(w)),
            )
        });
        assert_eq!(c.nonfinite, 0, "r at ({n2}, {w})");
        assert!(
            r.v.0.is_finite() && r.d[0].0.is_finite(),
            "r at ({n2}, {w})"
        );
    }
    // The detector sees what it should: every exact arm is not finite at `z = 0`.
    for (i, name) in NAMES.iter().enumerate() {
        let (_, c) = counted(|| arms(i, L::<F32>::variable(Lane::new(0.0), 0), 1));
        assert!(c.nonfinite > 0, "{name}: the exact arm at 0 is finite");
    }
}

#[test]
fn a_lane_that_evaluates_both_arms_sees_no_non_finite_operation_f64() {
    a_lane_that_evaluates_both_arms_sees_no_non_finite_operation::<false>();
}

#[test]
fn a_lane_that_evaluates_both_arms_sees_no_non_finite_operation_f32() {
    a_lane_that_evaluates_both_arms_sees_no_non_finite_operation::<true>();
}

/// Each group's `(sqrt, sin_cos, atan2)` calls, in the order `exp`, `jr`, `jr_inv`, `q`, `log`,
/// one call each under a scalar mask, at `z` (`n²` for `log`, at `w = 1`).
fn calls<const F32: bool>(z: f64) -> [(usize, usize, usize); 5] {
    let (z, w) = (Narrow::<F32>::new(z), Narrow::<F32>::one());
    [
        counted(|| exp_coeffs(z)).1,
        counted(|| jr_coeffs(z)).1,
        counted(|| jr_inv_coeff(z)).1,
        counted(|| q_coeffs(z)).1,
        counted(|| log_ratio(z, w)).1,
    ]
    .map(|c| (c.sqrt, c.sin_cos, c.atan2))
}

/// The coefficients of each group of `calls`, as indices of `NAMES`.
const MEMBERS: [&[usize]; 5] = [&[0, 6], &[1, 2], &[3], &[2, 4, 5], &[7]];

/// While every member is on its series arm no exact arm runs: no `sqrt`, `sin_cos` or `atan2`,
/// at `θ² = 0` and below the smallest switch of the group.
fn a_group_below_its_smallest_switch_runs_no_exact_arm<const F32: bool>() {
    let table = table::<F32>();
    for (members, g) in MEMBERS.iter().zip(0..) {
        let least = members
            .iter()
            .map(|&i| table[i].0)
            .fold(f64::INFINITY, f64::min);
        for z in [0.0, 1e-30, 0.5 * least, least.next_down()] {
            assert_eq!(calls::<F32>(z)[g], (0, 0, 0), "group {g} at {z:e}");
        }
    }
}

#[test]
fn a_group_below_its_smallest_switch_runs_no_exact_arm_f64() {
    a_group_below_its_smallest_switch_runs_no_exact_arm::<false>();
}

#[test]
fn a_group_below_its_smallest_switch_runs_no_exact_arm_f32() {
    a_group_below_its_smallest_switch_runs_no_exact_arm::<true>();
}

/// Above the switches the exact arms of a group run once and share what they can: `Exp` one `sqrt` and
/// one `sin_cos` (`NUMERICS.md` §3.1) from `cos θ/2`'s switch up, `b` and `e` one `sin_cos` of `θ`.
fn a_group_runs_each_exact_arm_once<const F32: bool>() {
    // (sqrt, sin_cos, atan2) of `exp`, `jr`, `jr_inv`, `q`, `log`.
    let above = [(1, 1, 0), (2, 2, 0), (1, 1, 0), (1, 2, 0), (1, 0, 1)];
    for z in [2.0, 6.0, 9.8] {
        assert_eq!(calls::<F32>(z), above, "z = {z}");
    }
    // `cos θ/2` alone is exact between its switch and `k`'s: at both precisions, at `z = 1e-3`.
    assert_eq!(calls::<F32>(1e-3)[0], (1, 1, 0));
    if !F32 {
        // The one `sin_cos` for nothing (`super`): `b` on its series arm, `a` not, and `d` likewise.
        assert_eq!(calls::<false>(0.8)[1], (2, 2, 0));
        assert_eq!(calls::<false>(0.97)[3], (1, 2, 0));
    }
}

#[test]
fn a_group_runs_each_exact_arm_once_f64() {
    a_group_runs_each_exact_arm_once::<false>();
}

#[test]
fn a_group_runs_each_exact_arm_once_f32() {
    a_group_runs_each_exact_arm_once::<true>();
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "the branch variable is not >= 0")]
fn a_negative_branch_variable_is_refused_in_debug() {
    let _ = exp_coeffs(-1.0_f64);
}
