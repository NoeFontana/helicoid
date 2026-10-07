//! Behavioural tests of the scalar model. Floats are compared by bit pattern, so signed zeros
//! and NaN payloads are distinguished. Signed-zero and argument-order cases go through generic
//! helpers: on a concrete float `x.atan2(y)` resolves to the inherent method, not to `Real`.

use crate::{Blend, Mask, Precision, Real};
use core::cell::Cell;
use core::ops::{Add, Div, Mul, Neg, Sub};

fn bits(x: f64) -> u64 {
    x.to_bits()
}

const F64_GRID: [f64; 14] = [
    0.0,
    -0.0,
    1.0,
    -1.0,
    0.5,
    2.0,
    3.0,
    1e-300,
    1e300,
    core::f64::consts::PI,
    -core::f64::consts::FRAC_PI_2,
    0.1,
    123.456,
    -7.0e-5,
];

#[test]
fn bool_mask_truth_tables() {
    for a in [false, true] {
        assert_eq!(Mask::not(a), !a);
        assert_eq!(a.all(), a);
        assert_eq!(a.any(), a);
        for b in [false, true] {
            assert_eq!(Mask::and(a, b), a && b);
            assert_eq!(Mask::or(a, b), a || b);
        }
    }
}

#[test]
fn bool_decide_evaluates_exactly_one_arm() {
    for m in [false, true] {
        let (t_calls, f_calls, blend_calls) = (Cell::new(0), Cell::new(0), Cell::new(0));
        let out = m.decide(
            || {
                t_calls.set(t_calls.get() + 1);
                1
            },
            || {
                f_calls.set(f_calls.get() + 1);
                2
            },
            |_, a, _| {
                blend_calls.set(blend_calls.get() + 1);
                a
            },
        );
        assert_eq!(out, if m { 1 } else { 2 });
        assert_eq!(
            (t_calls.get(), f_calls.get()),
            (u32::from(m), u32::from(!m))
        );
        assert_eq!(blend_calls.get(), 0);
    }
}

#[test]
fn comparisons_and_nan() {
    assert!(1.0_f64.lt(2.0));
    assert!(!2.0_f64.lt(2.0));
    assert!(2.0_f64.le(2.0));
    assert!(!3.0_f64.le(2.0));
    assert!(0.0_f64.le(-0.0) && !0.0_f64.lt(-0.0));
    let nan = f64::NAN;
    assert!(!nan.lt(1.0) && !nan.le(1.0) && !1.0_f64.lt(nan) && !1.0_f64.le(nan));
    assert!(1.0_f32.lt(2.0) && !nan.lt(nan));
}

#[test]
fn select_and_branch() {
    assert_eq!(bits(f64::select(true, 1.0, 2.0)), bits(1.0));
    assert_eq!(bits(f64::select(false, 1.0, 2.0)), bits(2.0));
    // A NaN in the unselected operand does not leak.
    assert_eq!(bits(f64::select(true, 1.0, f64::NAN)), bits(1.0));
    assert_eq!(bits(f64::select(false, f64::NAN, -0.0)), bits(-0.0));
    assert_eq!(bits(f64::branch(true, || 3.0, || 4.0)), bits(3.0));
    assert_eq!(bits(f64::branch(false, || 3.0, || 4.0)), bits(4.0));
}

#[test]
fn branch_runs_one_arm_for_scalars() {
    let calls = Cell::new(0);
    let out = f64::branch(
        true,
        || 1.0,
        || {
            calls.set(calls.get() + 1);
            2.0
        },
    );
    assert_eq!((bits(out), calls.get()), (bits(1.0), 0));
}

#[test]
fn blend_tuples_and_arrays() {
    type T = (f64, f64, [f64; 2], (f64,));
    let t: T = (1.0, 2.0, [3.0, 4.0], (5.0,));
    let f: T = (-1.0, -2.0, [-3.0, -4.0], (-5.0,));
    let flat = |x: T| [x.0, x.1, x.2[0], x.2[1], x.3 .0].map(f64::to_bits);
    assert_eq!(
        flat(<T as Blend<f64>>::blend(true, t, f)),
        [1.0, 2.0, 3.0, 4.0, 5.0].map(f64::to_bits)
    );
    assert_eq!(
        flat(<T as Blend<f64>>::blend(false, t, f)),
        [-1.0, -2.0, -3.0, -4.0, -5.0].map(f64::to_bits)
    );
}

/// Compiles only if `Blend<S>` is implemented for the tuple; one instantiation per arity.
fn assert_blend<S: Real, T: Blend<S>>() {}

#[test]
fn blend_exists_for_arities_one_to_eight() {
    type F = f64;
    assert_blend::<F, (F,)>();
    assert_blend::<F, (F, F)>();
    assert_blend::<F, (F, F, F)>();
    assert_blend::<F, (F, F, F, F)>();
    assert_blend::<F, (F, F, F, F, F)>();
    assert_blend::<F, (F, F, F, F, F, F)>();
    assert_blend::<F, (F, F, F, F, F, F, F)>();
    assert_blend::<F, (F, F, F, F, F, F, F, F)>();
}

#[test]
fn branch_returns_tuple_of_eight() {
    let t = [1.0_f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let f = [-1.0_f32, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0];
    let tup = |v: [f32; 8]| (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
    let arr = |x: (f32, f32, f32, f32, f32, f32, f32, f32)| {
        [x.0, x.1, x.2, x.3, x.4, x.5, x.6, x.7].map(f32::to_bits)
    };
    assert_eq!(
        arr(f32::branch(true, || tup(t), || tup(f))),
        t.map(f32::to_bits)
    );
    assert_eq!(
        arr(f32::branch(false, || tup(t), || tup(f))),
        f.map(f32::to_bits)
    );
}

#[test]
fn constants_and_precision() {
    assert_eq!(f64::PRECISION, Precision::F64);
    assert_eq!(f32::PRECISION, Precision::F32);
    assert_eq!(bits(f64::zero()), bits(0.0));
    assert_eq!(bits(f64::one()), bits(1.0));
    assert_eq!(f32::zero().to_bits(), 0.0_f32.to_bits());
    assert_eq!(f32::one().to_bits(), 1.0_f32.to_bits());
    assert_eq!(bits(f32::lit(0.5).value_f64()), bits(0.5));
    assert_eq!(bits(f32::lit(-3.0).value_f64()), bits(-3.0));
    assert_eq!(bits(f64::lit(0.1)), bits(0.1));
    assert_eq!(bits(f32::lit(-0.0).value_f64()), bits(-0.0));
}

// `f32::lit(0.1)` rounds: a `debug_assert!` fires, since generated constants own that job.
#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "exactly representable")]
fn lit_rejects_inexact_f32_constant() {
    let _ = f32::lit(0.1);
}

#[test]
fn value_f64_widens_exactly() {
    assert_eq!(bits(0.1_f32.value_f64()), bits(f64::from(0.1_f32)));
    assert_eq!(bits(0.1_f64.value_f64()), bits(0.1));
}

// Generic helpers force dispatch through `Real`.
fn sin_cos<S: Real>(x: S) -> (S, S) {
    x.sin_cos()
}
fn atan2<S: Real>(y: S, x: S) -> S {
    y.atan2(x)
}
fn abs<S: Real>(x: S) -> S {
    x.abs()
}
fn sqrt<S: Real>(x: S) -> S {
    x.sqrt()
}
fn cbrt<S: Real>(x: S) -> S {
    x.cbrt()
}
fn copysign<S: Real>(x: S, sign: S) -> S {
    x.copysign(sign)
}

/// IEEE signed-zero and argument-order results, at any precision. `pi` and `half_pi` are the
/// correctly rounded constants of `S`, widened to `f64`.
fn ieee_corners<S: Real>(pi: f64, half_pi: f64) {
    let z = |x: f64| S::lit(x);
    let b = |x: S| bits(x.value_f64());
    let (s, c) = sin_cos(z(-0.0));
    assert_eq!((b(s), b(c)), (bits(-0.0), bits(1.0)));
    let (s, c) = sin_cos(z(0.0));
    assert_eq!((b(s), b(c)), (bits(0.0), bits(1.0)));
    // atan2(y, x) at the zero corners.
    assert_eq!(b(atan2(z(0.0), z(0.0))), bits(0.0));
    assert_eq!(b(atan2(z(-0.0), z(0.0))), bits(-0.0));
    assert_eq!(b(atan2(z(0.0), z(-0.0))), bits(pi));
    assert_eq!(b(atan2(z(-0.0), z(-0.0))), bits(-pi));
    assert_eq!(b(atan2(z(-0.0), z(1.0))), bits(-0.0));
    assert_eq!(b(atan2(z(0.0), z(-1.0))), bits(pi));
    assert_eq!(b(atan2(z(-0.0), z(-1.0))), bits(-pi));
    // The angle of the point (x, y) is atan2(y, x): the receiver is y.
    assert_eq!(b(atan2(z(1.0), z(0.0))), bits(half_pi));
    assert_eq!(b(atan2(z(-1.0), z(0.0))), bits(-half_pi));
    assert_eq!(b(atan2(z(0.0), z(1.0))), bits(0.0));
    assert_eq!(b(atan2(z(0.0), z(-1.0))), bits(pi));
    assert_eq!(b(abs(z(-0.0))), bits(0.0));
    assert_eq!(b(abs(z(-2.0))), bits(2.0));
    assert_eq!(b(copysign(z(1.0), z(-0.0))), bits(-1.0));
    assert_eq!(b(copysign(z(-1.0), z(0.0))), bits(1.0));
    assert_eq!(b(copysign(z(0.0), z(-1.0))), bits(-0.0));
    assert_eq!(b(copysign(z(-0.0), z(1.0))), bits(0.0));
    assert_eq!(b(sqrt(z(-0.0))), bits(-0.0));
    assert_eq!(b(sqrt(z(4.0))), bits(2.0));
    // `cbrt` is total: a zero and an infinity keep their sign, a cube is recovered exactly.
    let inf = S::one() / S::zero();
    assert_eq!(b(cbrt(z(-0.0))), bits(-0.0));
    assert_eq!(b(cbrt(z(0.0))), bits(0.0));
    assert_eq!(b(cbrt(z(-27.0))), bits(-3.0));
    assert_eq!(b(cbrt(z(0.125))), bits(0.5));
    assert_eq!(b(cbrt(inf)), bits(f64::INFINITY));
    assert_eq!(b(cbrt(-inf)), bits(f64::NEG_INFINITY));
    assert!(cbrt(S::zero() / S::zero()).value_f64().is_nan());
}

#[test]
fn ieee_corners_f64() {
    ieee_corners::<f64>(core::f64::consts::PI, core::f64::consts::FRAC_PI_2);
}

#[test]
fn ieee_corners_f32() {
    ieee_corners::<f32>(
        f64::from(core::f32::consts::PI),
        f64::from(core::f32::consts::FRAC_PI_2),
    );
}

/// A deterministic argument set: the grid, subnormals and extremes, arguments straddling
/// multiples of pi/2, and 4096 xorshift points over +-3000. The libm/std spread on `sin_cos` is
/// about 3% of random arguments, so a std or architecture substitution fails this sweep.
fn for_each_point(mut visit: impl FnMut(f64)) {
    let half_pi = core::f64::consts::FRAC_PI_2;
    for x in F64_GRID {
        visit(x);
    }
    for x in [
        5e-324,
        f64::MIN_POSITIVE,
        f64::MAX,
        f64::EPSILON,
        1e22,
        1.4e-45,
        1.2e-38,
        3.4e38,
    ] {
        visit(x);
        visit(-x);
    }
    let mut k = 1.0_f64;
    while k < 2e6 {
        let x = k * half_pi;
        visit(x);
        visit(f64::from_bits(x.to_bits() + 1));
        visit(f64::from_bits(x.to_bits() - 1));
        k = k * 7.0 + 1.0;
    }
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    for _ in 0..4096 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        visit((state >> 11) as f64 / 9_007_199_254_740_992.0 * 6000.0 - 3000.0);
    }
}

#[test]
fn f64_routes_through_libm_bit_for_bit() {
    let mut prev = 0.5_f64;
    for_each_point(|x| {
        let (s, c) = libm::sincos(x);
        let (rs, rc) = sin_cos(x);
        assert_eq!((bits(rs), bits(rc)), (bits(s), bits(c)), "sin_cos({x:e})");
        if x >= 0.0 {
            assert_eq!(bits(sqrt(x)), bits(libm::sqrt(x)), "sqrt({x:e})");
        }
        assert_eq!(bits(cbrt(x)), bits(libm::cbrt(x)), "cbrt({x:e})");
        assert_eq!(bits(abs(x)), bits(libm::fabs(x)));
        for y in [prev, 1.0, -1.0, 0.0] {
            assert_eq!(
                bits(atan2(x, y)),
                bits(libm::atan2(x, y)),
                "atan2({x:e}, {y:e})"
            );
            assert_eq!(bits(copysign(x, y)), bits(libm::copysign(x, y)));
        }
        prev = x;
    });
}

#[test]
fn f32_routes_through_libm_bit_for_bit() {
    let mut prev = 0.5_f32;
    for_each_point(|x64| {
        let x = x64 as f32;
        let (s, c) = libm::sincosf(x);
        let (rs, rc) = sin_cos(x);
        assert_eq!(
            (rs.to_bits(), rc.to_bits()),
            (s.to_bits(), c.to_bits()),
            "sin_cos({x:e})"
        );
        if x >= 0.0 {
            assert_eq!(sqrt(x).to_bits(), libm::sqrtf(x).to_bits(), "sqrt({x:e})");
        }
        assert_eq!(cbrt(x).to_bits(), libm::cbrtf(x).to_bits(), "cbrt({x:e})");
        assert_eq!(abs(x).to_bits(), libm::fabsf(x).to_bits());
        for y in [prev, 1.0, -1.0, 0.0] {
            assert_eq!(
                atan2(x, y).to_bits(),
                libm::atan2f(x, y).to_bits(),
                "atan2({x:e}, {y:e})"
            );
            assert_eq!(copysign(x, y).to_bits(), libm::copysignf(x, y).to_bits());
        }
        prev = x;
    });
}

/// `cbrt` is odd to the bit at every point of the sweep, subnormals and extremes included, and
/// `libm` rounds `-x` to the negation of `x`'s result, so `Dual::cbrt`'s derivative is even.
fn cbrt_is_odd<S: Real>(conv: fn(f64) -> S) {
    for_each_point(|x| {
        let (p, n) = (cbrt(conv(x)).value_f64(), cbrt(conv(-x)).value_f64());
        assert_eq!(n.to_bits(), (-p).to_bits(), "cbrt(-x) at {x:e}");
    });
}

#[test]
fn cbrt_is_odd_f64() {
    cbrt_is_odd::<f64>(|x| x);
}

#[test]
fn cbrt_is_odd_f32() {
    cbrt_is_odd::<f32>(|x| x as f32);
}

#[test]
fn sin_cos_pythagoras() {
    for &x in &F64_GRID {
        let (s, c) = Real::sin_cos(x);
        assert!((s * s + c * c - 1.0).abs() < 1e-15);
    }
}

/// The safe-argument pattern, written once for every `Real`: finite in the selected arm even
/// where the exact arm would divide by zero.
fn sinc<S: Real>(x: S) -> S {
    let small = x.abs().lt(S::lit(0.5));
    let safe = S::select(small, S::one(), x);
    S::branch(
        small,
        || S::one() - x * x / S::lit(6.0),
        || safe.sin_cos().0 / safe,
    )
}

#[test]
fn generic_code_runs_at_both_precisions() {
    assert_eq!(bits(sinc(0.0_f64)), bits(1.0));
    assert!((sinc(2.0_f64) - libm::sin(2.0) / 2.0).abs() < 1e-15);
    assert_eq!(sinc(0.0_f32).to_bits(), 1.0_f32.to_bits());
    assert!((sinc(2.0_f32).value_f64() - libm::sin(2.0) / 2.0).abs() < 1e-6);
}

/// Two lanes: a test-only `Mask` and `Real` that evaluate both arms and blend, so the lane path
/// of `Real::branch` (never taken by `bool`) is exercised.
#[derive(Clone, Copy)]
pub(crate) struct M2(pub(crate) [bool; 2]);

impl Mask for M2 {
    fn and(self, o: Self) -> Self {
        Self([self.0[0] & o.0[0], self.0[1] & o.0[1]])
    }
    fn or(self, o: Self) -> Self {
        Self([self.0[0] | o.0[0], self.0[1] | o.0[1]])
    }
    fn not(self) -> Self {
        Self([!self.0[0], !self.0[1]])
    }
    fn all(self) -> bool {
        self.0[0] && self.0[1]
    }
    fn any(self) -> bool {
        self.0[0] || self.0[1]
    }
    fn decide<T>(
        self,
        t: impl FnOnce() -> T,
        f: impl FnOnce() -> T,
        blend: impl FnOnce(Self, T, T) -> T,
    ) -> T {
        blend(self, t(), f())
    }
}

#[derive(Clone, Copy)]
pub(crate) struct L2(pub(crate) [f64; 2]);

impl L2 {
    fn map(self, g: impl Fn(f64) -> f64) -> Self {
        Self([g(self.0[0]), g(self.0[1])])
    }
    fn zip(self, o: Self, g: impl Fn(f64, f64) -> f64) -> Self {
        Self([g(self.0[0], o.0[0]), g(self.0[1], o.0[1])])
    }
}

impl Add for L2 {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        self.zip(o, |a, b| a + b)
    }
}
impl Sub for L2 {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        self.zip(o, |a, b| a - b)
    }
}
impl Mul for L2 {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        self.zip(o, |a, b| a * b)
    }
}
impl Div for L2 {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        self.zip(o, |a, b| a / b)
    }
}
impl Neg for L2 {
    type Output = Self;
    fn neg(self) -> Self {
        self.map(|a| -a)
    }
}

impl Real for L2 {
    type Mask = M2;
    const PRECISION: Precision = Precision::F64;
    fn lit(x: f64) -> Self {
        Self([x, x])
    }
    fn zero() -> Self {
        Self::lit(0.0)
    }
    fn one() -> Self {
        Self::lit(1.0)
    }
    fn lt(self, r: Self) -> M2 {
        M2([self.0[0] < r.0[0], self.0[1] < r.0[1]])
    }
    fn le(self, r: Self) -> M2 {
        M2([self.0[0] <= r.0[0], self.0[1] <= r.0[1]])
    }
    fn select(m: M2, t: Self, f: Self) -> Self {
        Self([
            f64::select(m.0[0], t.0[0], f.0[0]),
            f64::select(m.0[1], t.0[1], f.0[1]),
        ])
    }
    fn sqrt(self) -> Self {
        self.map(sqrt)
    }
    fn cbrt(self) -> Self {
        self.map(cbrt)
    }
    fn sin_cos(self) -> (Self, Self) {
        (self.map(|a| sin_cos(a).0), self.map(|a| sin_cos(a).1))
    }
    fn sin(self) -> Self {
        self.map(Real::sin)
    }
    fn cos(self) -> Self {
        self.map(Real::cos)
    }
    fn acos(self) -> Self {
        self.map(Real::acos)
    }
    fn atan2(self, x: Self) -> Self {
        self.zip(x, atan2)
    }
    fn abs(self) -> Self {
        self.map(abs)
    }
    fn copysign(self, s: Self) -> Self {
        self.zip(s, copysign)
    }
    fn value_f64(self) -> f64 {
        self.0[0]
    }
}

#[test]
fn lane_branch_evaluates_both_arms_and_blends_per_lane() {
    let m = M2([true, false]);
    let (t_calls, f_calls) = (Cell::new(0), Cell::new(0));
    let out = L2::branch(
        m,
        || {
            t_calls.set(t_calls.get() + 1);
            (L2([1.0, 2.0]), [L2([3.0, 4.0])])
        },
        || {
            f_calls.set(f_calls.get() + 1);
            (L2([-1.0, -2.0]), [L2([-3.0, -4.0])])
        },
    );
    assert_eq!((t_calls.get(), f_calls.get()), (1, 1));
    assert_eq!(out.0 .0.map(f64::to_bits), [1.0, -2.0].map(f64::to_bits));
    assert_eq!(out.1[0].0.map(f64::to_bits), [3.0, -4.0].map(f64::to_bits));
}

#[test]
fn lane_sinc_takes_the_right_arm_per_lane() {
    let out = sinc(L2([0.0, 2.0]));
    assert_eq!(bits(out.0[0]), bits(1.0));
    assert_eq!(bits(out.0[1]), bits(libm::sin(2.0) / 2.0));
}
