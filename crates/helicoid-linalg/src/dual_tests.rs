//! Tests of `Dual`. The value path is compared with plain evaluation by bit pattern; each
//! derivative rule against an independent derivative (exact arithmetic, mpmath, a five-point
//! stencil of the plain function); nesting supplies second derivatives.

use crate::tests::{L2, M2};
use crate::{chol, Dual, Mask, Matrix, Precision, Real, Vector};
use core::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, PI};

type D1 = Dual<f64, 1>;
type D2 = Dual<f64, 2>;
type Hess = Dual<Dual<f64, 2>, 2>;

#[rustfmt::skip]
const POINTS: &[f64] = &[
    0.0, -0.0, 1.0, -1.0, 0.5, -2.0, 3.0, 0.1, -123.456, 7e-5, PI, -FRAC_PI_2, f64::EPSILON,
    1.0 + f64::EPSILON, 1e-300, -1e300, 5e-324, -5e-324, f64::MIN_POSITIVE, f64::MAX, f64::MIN,
    1e22, 1.4e-45, 3.4e38, f64::INFINITY, f64::NEG_INFINITY, f64::NAN,
];

/// The output names of `probe`, and whether a NaN there is compared by bit pattern. The sign
/// and payload of a NaN from arithmetic or `libm` are unspecified in Rust, and a release build may
/// commute operands (`NaN * -NaN`), so those outputs treat every NaN as one value. Sign-bit
/// operations and selection are deterministic and compared exactly.
#[rustfmt::skip]
/// How many expressions [`probe`] returns and [`OUTPUTS`] names.
///
/// One constant for both, because the two were `[S; 29]` against `[(&str, bool); 26]` for one
/// commit and `zip` truncates in silence: three expressions went uncompared and every name from
/// `sin` on was attached to the wrong one, so `x.sin()` was reported as `lt` and compared *exactly*
/// where its flag says otherwise. A shared length makes that a compile error (`0053`).
const PROBES: usize = 29;

const OUTPUTS: [(&str, bool); PROBES] = [
    ("add", false),
    ("sub", false),
    ("mul", false),
    ("div", false),
    ("neg", true),
    ("lit", true),
    ("lit", true),
    ("zero", true),
    ("one", true),
    ("sin", false),
    ("cos", false),
    ("sqrt", false),
    ("cbrt", false),
    ("atan2", false),
    ("abs", true),
    ("copysign", true),
    ("sin", false),
    ("cos", false),
    ("acos", false),
    ("lt", true),
    ("le", true),
    ("not", true),
    ("and", true),
    ("or", true),
    ("select", true),
    ("branch", false),
    ("branch", false),
    ("chain", false),
    ("chain", false),
];

const _: () = assert!(
    OUTPUTS.len() == PROBES,
    "every `probe` expression needs its name and exactness flag"
);

/// One output per `Real` method and operator, so a mismatch names the method. `sqrt` takes the
/// `abs` of its argument: its domain is non-negative.
#[rustfmt::skip]
fn probe<S: Real>(x: S, y: S, z: S) -> [S; PROBES] {
    let (s, c) = x.sin_cos();
    let (lt, le) = (x.lt(y), y.le(z));
    let pick = |m: S::Mask| S::select(m, S::one(), S::zero());
    let (bt, bf) = S::branch(lt, || (x + z, y * z), || (x - z, y / z));
    [
        x + y, x - y, x * y, x / y, -x, S::lit(0.5), S::lit(-3.0), S::zero(), S::one(),
        s, c, x.abs().sqrt(), x.cbrt(), x.atan2(y), x.abs(), x.copysign(y),
        // `sin`, `cos` and `acos`: this list is written by hand, so a `Real` method added without
        // extending it is one `dual_value_is_plain_value` silently stops covering -- which is what
        // happened to `sin` in `0052` and is corrected here along with `0022`'s two. Most `POINTS`
        // are outside `acos`'s domain, so most of its rows are NaN, which is the comparison this
        // test is built to make ("up to NaN sign and payload", `PHASE2.md` §3).
        x.sin(), x.cos(), x.acos(),
        pick(lt), pick(le), pick(lt.not()), pick(lt.and(le)), pick(lt.or(le)),
        S::select(le, x, z), bt, bf,
        (x * y - z).atan2(x + z), (x.abs() + y.abs()).sqrt() * z.copysign(x),
    ]
}

fn value_path<S: Real, D: Real>(label: &str, conv: fn(f64) -> S, lift: fn(S, usize) -> D) {
    for &a in POINTS {
        for &b in POINTS {
            for &c in POINTS {
                let (x, y, z) = (conv(a), conv(b), conv(c));
                let plain = probe(x, y, z);
                let dual = probe(lift(x, 0), lift(y, 1), lift(z, 2));
                for ((p, d), (name, exact)) in plain.iter().zip(&dual).zip(OUTPUTS) {
                    let (p, d) = (p.value_f64(), d.value_f64());
                    let same = p.to_bits() == d.to_bits() || (!exact && p.is_nan() && d.is_nan());
                    assert!(
                        same,
                        "{label}: {name} at ({a:e}, {b:e}, {c:e}): {p:e} vs {d:e}"
                    );
                }
            }
        }
    }
}

#[test]
fn dual_value_is_plain_value() {
    let poison = |v: f64, _: usize| Dual::<f64, 3> {
        v,
        d: [f64::NAN, f64::INFINITY, -1e300],
    };
    value_path("f64", |x| x, Dual::<f64, 3>::variable);
    value_path("f64, poisoned derivatives", |x| x, poison);
    value_path("f32", |x| x as f32, Dual::<f32, 3>::variable);
    value_path(
        "nested",
        |x| x,
        |v, i| Dual::<Dual<f64, 2>, 3>::variable(Dual::variable(v, i % 2), i),
    );
}

#[test]
fn constants_variables_and_precision() {
    assert_eq!(Dual::<f32, 2>::PRECISION, Precision::F32);
    assert_eq!(Dual::<Dual<f64, 1>, 2>::PRECISION, Precision::F64);
    let v = Dual::<f64, 3>::variable(2.0, 1);
    assert_eq!(v.v.to_bits(), 2.0_f64.to_bits());
    assert_eq!(v.d.map(f64::to_bits), [0, 1.0_f64.to_bits(), 0]);
    for (c, value) in [
        (D2::constant(-0.5), -0.5),
        (D2::lit(-0.5), -0.5),
        (D2::zero(), 0.0),
        (D2::one(), 1.0),
    ] {
        assert_eq!(
            (c.v.to_bits(), c.d.map(f64::to_bits)),
            (f64::to_bits(value), [0, 0])
        );
    }
}

// The release build of `variable` with `i >= N` returns a constant instead.
#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "index out of range")]
fn variable_index_is_checked_in_debug() {
    let _ = D2::variable(1.0, 2);
}

fn flat(r: D2) -> [u64; 3] {
    [r.v, r.d[0], r.d[1]].map(f64::to_bits)
}

#[test]
fn arithmetic_rules_are_exact_on_dyadic_inputs() {
    let (x, y) = (D2::variable(3.0, 0), D2::variable(0.5, 1));
    let want = |v: [f64; 3]| v.map(f64::to_bits);
    assert_eq!(flat(x + y), want([3.5, 1.0, 1.0]));
    assert_eq!(flat(x - y), want([2.5, 1.0, -1.0]));
    assert_eq!(flat(x * y), want([1.5, 0.5, 3.0]));
    assert_eq!(flat(x / y), want([6.0, 2.0, -12.0]));
    assert_eq!(flat(-x), want([-3.0, -1.0, -0.0]));
    let (x, y) = (D2::variable(3.0, 0), D2::variable(4.0, 1));
    assert_eq!(flat((x * x + y * y).sqrt()), want([5.0, 0.6, 0.8]));
}

/// The distance in representable values, and the assertion over it, once per width.
///
/// Both precisions share one total-order key (negatives inverted, non-negatives sign-bit set) so
/// that a change to how it treats signed zeros cannot reach one and miss the other.
macro_rules! ulps_for {
    ($ulps:ident, $near:ident, $f:ty, $u:ty, $sign:expr) => {
        /// The distance in representable `
        #[doc = stringify!($f)]
        /// `s.
        fn $ulps(a: $f, b: $f) -> u64 {
            let key = |x: $f| {
                let bits = x.to_bits();
                if bits >> $sign == 1 {
                    !bits
                } else {
                    bits | (1 << $sign)
                }
            };
            u64::from(key(a).abs_diff(key(b)))
        }

        fn $near(got: $f, want: $f, max_ulps: u64, what: core::fmt::Arguments<'_>) {
            assert!(
                $ulps(got, want) <= max_ulps,
                "{what}: got {got:e}, want {want:e}"
            );
        }
    };
}

ulps_for!(ulps, near, f64, u64, 63);
ulps_for!(ulps32, near32, f32, u32, 31);

// The reference values below are `mp.diff` (central, step 2^-50 times the input scale) of
// `mp.sqrt`, `mp.sin`/`mp.cos`, `mp.atan2` and the quotient and product at 80 digits from exact
// binary64 inputs, rounded to nearest; each was cross-checked against the closed form to 20
// digits. Rows: `(x, sqrt' x)`, `(x, sin' x, cos' x)`, `(y, x, d/dy, d/dx, d2/dy2, d2/dydx,
// d2/dx2)` of `atan2(y, x)`, `(y, x, d/dy, d/dx)` of `atan2(y, x)` at the edge of its derivative
// domain, `(a, b, d/da, d/db)` of `a / b`, and `(t, k1, k2, d/dt)` of `(k1 t) (k2 t)`.
#[rustfmt::skip]
const SQRT: [(f64, f64); 10] = [
    (5e-324, 2.2494568972715982e+161), (1e-300, 5e+149), (1e-08, 5000.0), (0.25, 1.0),
    (0.5, FRAC_1_SQRT_2), (1.0, 0.5), (2.0, 0.3535533905932738),
    (3.0, 0.28867513459481287), (100000.0, 0.0015811388300841897), (1e+300, 5e-151),
];
#[rustfmt::skip]
const SIN_COS: [(f64, f64, f64); 8] = [
    (-2.5, -0.8011436155469337, 0.5984721441039565), (1e-08, 1.0, -1e-08),
    (0.7, 0.7648421872844885, -0.644217687237691),
    (FRAC_PI_2, 6.123233995736766e-17, -1.0),
    (PI, -1.0, -1.2246467991473532e-16),
    (100.0, 0.8623188722876839, 0.5063656411097588),
    (100000.0, -0.9993608074382124, -0.03574879797201651),
    (-7.5, 0.3466353178350258, 0.9379999767747389),
];
#[rustfmt::skip]
const ATAN2: [[f64; 7]; 9] = [
    [1.0, 2.0, 0.4, -0.2, -0.16, -0.12, 0.16],
    [-3.0, 0.5, 0.05405405405405406, 0.32432432432432434, 0.03506208911614317, 0.10226442658875091, -0.03506208911614317],
    [2.0, -7.0, -0.1320754716981132, -0.03773584905660377, 0.009967960128159488, -0.016019935920256318, -0.009967960128159488],
    [-0.25, -0.125, -1.6, 3.2, -10.24, 7.68, 10.24],
    [1e-05, 3.0, 0.3333333333296296, -1.1111111110987655e-06, -7.407407407242799e-07, -0.1111111111074074, 7.407407407242799e-07],
    [3.0, 1e-05, 1.1111111110987655e-06, -0.3333333333296296, -7.407407407242799e-07, 0.1111111111074074, 7.407407407242799e-07],
    [1e-100, 2e-100, 4e+99, -2e+99, -1.6e+199, -1.1999999999999999e+199, 1.6e+199],
    [-1e+150, 3e+150, 2.9999999999999998e-151, 1e-151, 6e-302, -8e-302, -6e-302],
    [0.6, 0.8, 0.8, -0.6, -0.96, -0.2800000000000001, 0.96],
];

#[rustfmt::skip]
const ATAN2_EDGE: [[f64; 4]; 4] = [
    [1e+153, 2e+153, 4e-154, -2e-154], [-3e+153, 4e+153, 1.6e-154, 1.2000000000000001e-154],
    [3e-154, 4e-154, 1.5999999999999999e+153, -1.2000000000000001e+153],
    [-4e-154, 3e-154, 1.2000000000000001e+153, 1.5999999999999999e+153],
];
#[rustfmt::skip]
const DIV: [[f64; 4]; 9] = [
    [1.0, 3.0, 0.3333333333333333, -0.1111111111111111],
    [5.52, -7.397, -0.135189941868325, -0.1008852885106332],
    [-1.35, -2.392, -0.4180602006688963, 0.23594534736747916],
    [1.56, 3.7, 0.27027027027027023, -0.1139517896274653],
    [0.766, -3.77, -0.26525198938992045, -0.05389470129248781],
    [-5.94, -0.65, -1.5384615384615383, 14.059171597633137],
    [1e-05, 3.0, 0.3333333333333333, -1.1111111111111112e-06],
    [1e+100, -7e-100, -1.4285714285714286e+99, -2.0408163265306123e+298],
    [1e-100, 3e-100, 3.333333333333333e+99, -1.111111111111111e+99],
];
#[rustfmt::skip]
const MUL: [[f64; 4]; 6] = [
    [1.1, 3.0, 0.7, 4.62], [0.3, 0.1, 7.0, 0.42], [-2.7, 1.3, 0.9, -6.3180000000000005],
    [123.456, 0.02, 5.1, 25.185024], [1e-3, 3e4, -1.7, -102.0], [5.5, -0.6, -0.35, 2.3099999999999996],
];

// The bound is 1 ulp everywhere; the measured maximum against these rows is 1 ulp for `sqrt`,
// `sin_cos`, `atan2`, the quotient and the product, and 0 on the `atan2` edge rows. Results are
// deterministic (D16). A quotient by a rounded reciprocal, for one, reaches 2 ulp on the `DIV` rows.
#[test]
fn dual_matches_mpmath_derivative() {
    for (x, want) in SQRT {
        let r = D1::variable(x, 0).sqrt();
        assert_eq!(r.v.to_bits(), libm::sqrt(x).to_bits());
        near(r.d[0], want, 1, format_args!("sqrt' at {x:e}"));
    }
    for (x, sin_d, cos_d) in SIN_COS {
        let (s, c) = D1::variable(x, 0).sin_cos();
        near(s.d[0], sin_d, 1, format_args!("sin' at {x:e}"));
        near(c.d[0], cos_d, 1, format_args!("cos' at {x:e}"));
    }
    for [y, x, dy, dx, ..] in ATAN2 {
        let r = D2::variable(y, 0).atan2(D2::variable(x, 1));
        assert_eq!(r.v.to_bits(), libm::atan2(y, x).to_bits());
        near(r.d[0], dy, 1, format_args!("atan2 d/dy at ({y:e}, {x:e})"));
        near(r.d[1], dx, 1, format_args!("atan2 d/dx at ({y:e}, {x:e})"));
    }
    for [a, b, da, db] in DIV {
        let r = D2::variable(a, 0) / D2::variable(b, 1);
        assert_eq!(r.v.to_bits(), (a / b).to_bits());
        near(
            r.d[0],
            da,
            1,
            format_args!("(a / b)' wrt a at ({a:e}, {b:e})"),
        );
        near(
            r.d[1],
            db,
            1,
            format_args!("(a / b)' wrt b at ({a:e}, {b:e})"),
        );
    }
    // Both factors depend on `t`, so each product rounds.
    for [t, k1, k2, want] in MUL {
        let x = D1::variable(t, 0);
        let r = (x * D1::constant(k1)) * (x * D1::constant(k2));
        assert_eq!(r.v.to_bits(), ((k1 * t) * (k2 * t)).to_bits());
        near(r.d[0], want, 1, format_args!("((k1 t)(k2 t))' at {t:e}"));
    }
}

// The derivative of `atan2` divides by `x^2 + y^2` (`Dual::atan2`, `# Domain`): accurate while the
// larger argument lies in about `1e-154 .. 1e154`, and silently wrong outside. The value is the
// plain one throughout.
#[test]
fn dual_atan2_derivative_is_accurate_to_the_edge_of_its_domain_and_only_there() {
    let at = |y: f64, x: f64| D2::variable(y, 0).atan2(D2::variable(x, 1));
    for [y, x, dy, dx] in ATAN2_EDGE {
        let r = at(y, x);
        assert_eq!(r.v.to_bits(), libm::atan2(y, x).to_bits());
        near(r.d[0], dy, 1, format_args!("atan2 d/dy at ({y:e}, {x:e})"));
        near(r.d[1], dx, 1, format_args!("atan2 d/dx at ({y:e}, {x:e})"));
    }
    let (over, under, origin) = (at(1e200, 2e200), at(1e-200, 2e-200), at(0.0, 0.0));
    // The true derivatives are `4e-201`, `4e199` (up to sign) and undefined.
    assert_eq!(over.d.map(|d| d.abs().to_bits()), [0, 0]);
    assert!(under.d.iter().all(|d| d.is_infinite()));
    assert!(origin.d.iter().all(|d| d.is_nan()));
    for (r, (y, x)) in [
        (over, (1e200, 2e200)),
        (under, (1e-200, 2e-200)),
        (origin, (0.0, 0.0)),
    ] {
        assert_eq!(r.v.to_bits(), libm::atan2(y, x).to_bits());
    }
}

// The quotient rule `(a' - q b') / b` overflows with `q = a / b` (`Dual::div` has no `b^2`): the
// value is infinite and the derivative NaN, though the true derivative `1e10` is finite.
#[test]
fn dual_quotient_derivative_is_nan_once_the_quotient_overflows() {
    let r = D1::variable(1e300, 0) / D1::constant(1e-10);
    assert!(r.v.is_infinite() && r.d[0].is_nan());
}

#[test]
fn nesting_gives_second_derivatives_of_sin_cos() {
    type D11 = Dual<Dual<f64, 1>, 1>;
    for (x, sin_d, cos_d) in SIN_COS {
        let (s, c) = D11::variable(Dual::variable(x, 0), 0).sin_cos();
        // s' = cos, s'' = -sin = cos'; c' = -sin, c'' = -cos = -sin'.
        near(s.v.d[0], sin_d, 1, format_args!("sin' at {x:e}"));
        near(s.d[0].v, sin_d, 1, format_args!("sin' (outer) at {x:e}"));
        near(s.d[0].d[0], cos_d, 1, format_args!("sin'' at {x:e}"));
        near(c.v.d[0], cos_d, 1, format_args!("cos' at {x:e}"));
        near(c.d[0].d[0], -sin_d, 1, format_args!("cos'' at {x:e}"));
    }
}

#[test]
fn nesting_gives_the_hessian_of_atan2() {
    let var = |v: f64, i: usize| Hess::variable(Dual::variable(v, i), i);
    for [y, x, dy, dx, hyy, hyx, hxx] in ATAN2 {
        let r = var(y, 0).atan2(var(x, 1));
        let close = |got: f64, want: f64, what: &str| {
            let ok = (got - want).abs() <= 1e-15 * want.abs();
            assert!(ok, "{what} at ({y:e}, {x:e}): got {got:e}, want {want:e}");
        };
        assert_eq!(r.v.v.to_bits(), libm::atan2(y, x).to_bits());
        close(r.v.d[0], dy, "d/dy");
        close(r.d[1].v, dx, "d/dx");
        close(r.d[0].d[0], hyy, "d2/dy2");
        close(r.d[0].d[1], hyx, "d2/dydx");
        close(r.d[1].d[0], hyx, "d2/dxdy");
        close(r.d[1].d[1], hxx, "d2/dx2");
    }
}

// Rows: `CBRT` is `(x, cbrt x, cbrt' x)`, `CBRT_D2` is `(x, cbrt' x, cbrt'' x)`, `CBRT_HESS` is
// `(x, y, f_x, f_y, f_xx, f_xy, f_yy)` of `f = cbrt(x y)` and `CBRT32` is `CBRT` at `f32`, rounded
// to nearest. From mpmath at 120 digits and exact inputs: the real `y` with `y^3 = x` by Newton's
// method on that equation (checked by cubing to 70 digits) and its derivatives by central
// differences with a step of `2^-40` times the argument. The closed forms `1 / (3 c^2)`,
// `-2 / (9 c^5)` and the partials of `cbrt(x) cbrt(y)` only cross-check them, to at least 18
// digits. Both signs; `f64` and `f32` subnormals, the smallest normal and the largest finite
// number. The generator is not committed.
#[rustfmt::skip]
const CBRT: [(f64, f64, f64); 23] = [
    (5e-324, 1.7031839360032603e-108, 1.1490942214795806e+215),
    (1e-310, 4.641588833612774e-104, 1.5471962778709295e+206),
    (2.2250738585072014e-308, 2.812644285236262e-103, 4.2135594353157945e+204),
    (1e-300, 1e-100, 3.3333333333333334e+199),
    (1e-08, 0.002154434690031884, 71814.48966772946),
    (0.001, 0.1, 33.333333333333336),
    (0.125, 0.5, 1.3333333333333333),
    (0.5, 0.7937005259840998, 0.5291336839893999),
    (1.0, 1.0, 0.3333333333333333),
    (2.0, 1.2599210498948732, 0.20998684164914552),
    (3.0, 1.4422495703074083, 0.1602499522563787),
    (27.0, 3.0, 0.037037037037037035),
    (100000.0, 46.415888336127786, 0.00015471962778709262),
    (1e+22, 21544346.90031884, 7.181448966772946e-16),
    (1e+300, 1e+100, 3.3333333333333335e-201),
    (1.7976931348623157e+308, 5.643803094122362e+102, 1.046489893941144e-206),
    (-5e-324, -1.7031839360032603e-108, 1.1490942214795806e+215),
    (-1e-310, -4.641588833612774e-104, 1.5471962778709295e+206),
    (-1e-300, -1e-100, 3.3333333333333334e+199),
    (-0.3, -0.6694329500821695, 0.7438143889801884),
    (-1.0, -1.0, 0.3333333333333333),
    (-27.0, -3.0, 0.037037037037037035),
    (-1e+300, -1e+100, 3.3333333333333335e-201),
];
#[rustfmt::skip]
const CBRT_D2: [(f64, f64, f64); 10] = [
    (1e-08, 71814.48966772946, -4787632644515.297),
    (0.125, 1.3333333333333333, -7.111111111111111),
    (0.5, 0.5291336839893999, -0.7055115786525331),
    (1.0, 0.3333333333333333, -0.2222222222222222),
    (2.0, 0.20998684164914552, -0.06999561388304851),
    (27.0, 0.037037037037037035, -0.0009144947416552355),
    (100000.0, 0.00015471962778709262, -1.0314641852472842e-09),
    (-0.125, 1.3333333333333333, 7.111111111111111),
    (-1.0, 0.3333333333333333, 0.2222222222222222),
    (-27.0, 0.037037037037037035, 0.0009144947416552355),
];
#[rustfmt::skip]
const CBRT_HESS: [[f64; 7]; 5] = [
    [2.0, 3.0, 0.3028534321386899, 0.20190228809245997, -0.10095114404622998, 0.033650381348743326, -0.04486717513165777],
    [0.5, 8.0, 1.0582673679787997, 0.06614171049867498, -1.4110231573050662, 0.04409447366578332, -0.005511809208222915],
    [-1.5, 4.0, 0.40380457618491994, -0.15142671606934496, 0.17946870052663108, 0.033650381348743326, 0.025237786011557496],
    [-0.25, -2.0, -1.0582673679787997, -0.13228342099734997, -2.8220463146101324, 0.17637789466313328, -0.04409447366578332],
    [0.001, 7.0, 63.764372759079635, 0.009109196108439948, -42509.581839386425, 3.036398702813316, -0.0008675424865180903],
];
#[rustfmt::skip]
const CBRT32: [(f32, f32, f32); 13] = [
    (1e-45_f32, 1.1190347e-15_f32, 2.6618994e+29_f32),
    (1e-40_f32, 4.6415806e-14_f32, 1.5472018e+26_f32),
    (1.1754944e-38_f32, 2.2737368e-13_f32, 6.4476046e+24_f32),
    (1e-30_f32, 1e-10_f32, 3.3333333e+19_f32),
    (0.125_f32, 0.5_f32, 1.3333334_f32),
    (1.0_f32, 1.0_f32, 0.33333334_f32),
    (3.0_f32, 1.4422495_f32, 0.16024995_f32),
    (27.0_f32, 3.0_f32, 0.037037037_f32),
    (1e+10_f32, 2154.4346_f32, 7.181449e-08_f32),
    (3.4028235e+38_f32, 6.9814636e+12_f32, 6.8388925e-27_f32),
    (-1e-30_f32, -1e-10_f32, 3.3333333e+19_f32),
    (-27.0_f32, -3.0_f32, 0.037037037_f32),
    (-0.3_f32, -0.66943294_f32, 0.74381435_f32),
];

// The bounds, 1 ulp on the value and 2 on every derivative, are those of these rows, not of the
// rule. Over these rows the value is 0 (`libm::cbrt` is correctly rounded, `cbrtf` agrees here),
// the derivative `d / (3 c^2)` at most 2, the second derivative and the Hessian entries at most 2.
// Random inputs against mpmath, `f64` and `f32` alike, reach about 4 ulp on the first derivative
// and about 10 on the second (the roundings of `c`, `c c`, `3 *` and `/`, and the nesting); the
// value stays at 0. First order runs at both precisions, second order at `f64`.
// `(x, acos x, d/dx acos x)` from mpmath at 120 digits, at exact inputs of each precision, each
// rounded to what its type holds. `0053`: `Dual::acos`'s rule shipped without one of these, which
// is the standard this module's own header sets ("each derivative rule against an independent
// derivative"), and `eig3`'s ambient Jacobian reaches it at a double eigenvalue.
#[rustfmt::skip]
const ACOS: [(f64, f64, f64); 5] = [
    // `acos` at `0` and `+-0.5` is `pi/2`, `pi/3` and `2 pi/3` exactly, so those rows would be
    // named constants spelled badly -- clippy's `approx_constant` says so -- and they test the
    // value, which `pi_and_acos_are_within_their_ulps` already pins against `libm` to the bit.
    // What is left is what the *derivative* needs: three generic points and the near-singular one.
    (0.25,  1.318_116_071_652_818,   -1.032_795_558_988_644_4),
    (-0.9,  2.690_565_841_793_530_8, -2.294_157_338_705_618),
    (0.9,   0.451_026_811_796_262_4, -2.294_157_338_705_618),
    (0.75,  0.722_734_247_813_415_6, -1.511_857_892_036_909_3),
    // `1 - 2^-20`: the derivative is -724, and the factored denominator earns its spelling here --
    // `sqrt(1 - x^2)` would lose half its digits.
    (0.999_999_046_325_683_6, 0.001_381_068_041_762_417_1, -724.077_516_568_578),
];

#[rustfmt::skip]
const ACOS32: [(f32, f32, f32); 5] = [
    // `acos` of the **binary32-rounded** argument, not of the binary64 one rounded afterwards:
    // `0.9_f32` is `0.899_999_976...`, and `acosf` of that differs from `acos(0.9)` in the last
    // place. A first draft of this table got it the other way round and failed against
    // `libm::acosf` by one ulp. Each literal is the shortest decimal that round-trips at binary32,
    // which is what clippy's `excessive_precision` asks for.
    (0.25,  1.318_116_1,  -1.032_795_5),
    (-0.9,  2.690_565_8,  -2.294_157),
    (0.9,   0.451_026_86, -2.294_157),
    (0.75,  0.722_734_3,  -1.511_857_9),
    (0.999_999_05, 0.001_381_068, -724.077_5),
];

/// `Dual::acos`'s value and derivative against mpmath, and the two degenerate cases its `# Domain`
/// names: `-+inf` where `d` is nonzero at `v = +-1` and NaN (`0 / 0`) where it is zero.
#[test]
fn dual_acos_matches_mpmath_derivative() {
    for (x, a, want) in ACOS {
        let r = D1::variable(x, 0).acos();
        assert_eq!(r.v.to_bits(), libm::acos(x).to_bits());
        near(r.v, a, 1, format_args!("acos at {x:e}"));
        near(r.d[0], want, 2, format_args!("acos' at {x:e}"));
    }
    for (x, a, want) in ACOS32 {
        let r = Dual::<f32, 1>::variable(x, 0).acos();
        assert_eq!(r.v.to_bits(), libm::acosf(x).to_bits());
        near32(r.v, a, 1, format_args!("f32 acos at {x:e}"));
        near32(r.d[0], want, 2, format_args!("f32 acos' at {x:e}"));
    }
    // `v = 1`: the seeded lane is `-inf` and the untouched ones NaN, which is the pair of cases
    // `Dual::cbrt` documents and `acos`'s domain did not until `0053`.
    let at_one = Dual::<f64, 2>::variable(1.0, 0).acos();
    assert_eq!(at_one.v.to_bits(), 0.0_f64.to_bits());
    assert!(at_one.d[0] == f64::NEG_INFINITY, "{:?}", at_one.d[0]);
    assert!(at_one.d[1].is_nan(), "{:?}", at_one.d[1]);
    // A constant at `v = 1`: every lane is `0 / 0`.
    let c = Dual::<f64, 2>::constant(1.0).acos();
    assert!(c.d.iter().all(|d| d.is_nan()), "{:?}", c.d);
    // Out of domain: the value is NaN and the derivative follows the same two cases, no panic even
    // in debug -- the `select` that feeds `sqrt` a safe argument (`0053` decision 7).
    let out = Dual::<f64, 2>::variable(2.0, 0).acos();
    assert!(out.v.is_nan() && out.d[0] == f64::NEG_INFINITY && out.d[1].is_nan());
}

/// Nesting gives the second derivative of `acos`: `-x (1 - x^2)^(-3/2)`.
#[test]
fn nesting_gives_second_derivatives_of_acos() {
    type D11 = Dual<Dual<f64, 1>, 1>;
    for (x, _, _) in ACOS {
        let inner = Dual::<f64, 1>::variable(x, 0);
        let r = D11::variable(inner, 0).acos();
        let s = 1.0 - x * x;
        let want = -x / (s * libm::sqrt(s));
        near(r.d[0].d[0], want, 8, format_args!("acos'' at {x:e}"));
    }
}

#[test]
fn dual_cbrt_matches_mpmath_derivative() {
    for (x, c, want) in CBRT {
        let r = D1::variable(x, 0).cbrt();
        assert_eq!(r.v.to_bits(), libm::cbrt(x).to_bits());
        near(r.v, c, 1, format_args!("cbrt at {x:e}"));
        near(r.d[0], want, 2, format_args!("cbrt' at {x:e}"));
    }
    for (x, c, want) in CBRT32 {
        let r = Dual::<f32, 1>::variable(x, 0).cbrt();
        assert_eq!(r.v.to_bits(), libm::cbrtf(x).to_bits());
        near32(r.v, c, 1, format_args!("f32 cbrt at {x:e}"));
        near32(r.d[0], want, 2, format_args!("f32 cbrt' at {x:e}"));
    }
}

#[test]
fn nesting_gives_second_derivatives_of_cbrt() {
    type D11 = Dual<Dual<f64, 1>, 1>;
    for (x, first, second) in CBRT_D2 {
        let r = D11::variable(Dual::variable(x, 0), 0).cbrt();
        near(r.v.d[0], first, 2, format_args!("cbrt' at {x:e}"));
        near(r.d[0].v, first, 2, format_args!("cbrt' (outer) at {x:e}"));
        near(r.d[0].d[0], second, 2, format_args!("cbrt'' at {x:e}"));
    }
}

// `cbrt(x y)` puts the chain rule through the product into the second order.
#[test]
fn nesting_gives_the_hessian_of_cbrt_of_a_product() {
    let var = |v: f64, i: usize| Hess::variable(Dual::variable(v, i), i);
    for [x, y, fx, fy, fxx, fxy, fyy] in CBRT_HESS {
        let r = (var(x, 0) * var(y, 1)).cbrt();
        let at = format_args!("at ({x:e}, {y:e})");
        near(r.v.d[0], fx, 2, format_args!("d/dx {at}"));
        near(r.d[1].v, fy, 2, format_args!("d/dy {at}"));
        near(r.d[0].d[0], fxx, 2, format_args!("d2/dx2 {at}"));
        near(r.d[0].d[1], fxy, 2, format_args!("d2/dxdy {at}"));
        near(r.d[1].d[0], fxy, 2, format_args!("d2/dydx {at}"));
        near(r.d[1].d[1], fyy, 2, format_args!("d2/dy2 {at}"));
    }
}

// `cbrt` is odd and `libm` rounds symmetrically, so the derivative `d / (3 c^2)` is even to the
// bit, zeros and infinities included; only NaN is left out, whose bits are unspecified.
#[test]
fn dual_cbrt_value_is_odd_and_derivative_is_even() {
    for &x in POINTS.iter().filter(|x| !x.is_nan()) {
        let (p, n) = (D1::variable(x, 0).cbrt(), D1::variable(-x, 0).cbrt());
        assert_eq!(n.v.to_bits(), (-p.v).to_bits(), "value at {x:e}");
        assert_eq!(n.d[0].to_bits(), p.d[0].to_bits(), "derivative at {x:e}");
    }
}

// Unlike `sqrt`, whose derivative at a zero has that zero's sign, `c^2` is `+0` for both zeros, so
// the infinity has the sign of `d` alone.
#[test]
fn dual_cbrt_at_zero_has_an_infinite_derivative_and_the_plain_value() {
    for zero in [0.0_f64, -0.0] {
        let up = D2::variable(zero, 0).cbrt();
        assert_eq!(up.v.to_bits(), zero.to_bits());
        assert!(up.d[0].is_infinite() && up.d[0] > 0.0);
        assert!(up.d[1].is_nan(), "0 / 0");
        let down = (-D1::variable(zero, 0)).cbrt();
        assert!(down.d[0].is_infinite() && down.d[0] < 0.0);
    }
    // The safe argument keeps the selected arm finite.
    let t = D1::variable(0.0, 0);
    let safe = D1::select(t.lt(D1::lit(0.25)), D1::one(), t);
    assert!(safe.cbrt().d[0].is_finite());
    // Past the finite range the derivative is the limit `0`, and NaN stays NaN.
    let inf = D1::variable(f64::INFINITY, 0).cbrt();
    assert_eq!(
        (inf.v.to_bits(), inf.d[0].to_bits()),
        (f64::INFINITY.to_bits(), 0)
    );
    let nan = D1::variable(f64::NAN, 0).cbrt();
    assert!(nan.v.is_nan() && nan.d[0].is_nan());
}

/// Every rule in one function, so the chain rule between them is exercised too.
fn composite<S: Real>(x: S, y: S) -> S {
    let (s, c) = x.sin_cos();
    let r = (x * x + y * y).sqrt();
    y.atan2(x) * r + s / (c + S::lit(2.0)) - x.copysign(y) * y.abs()
        + (x - y) / (r + S::one())
        + (x * y).cbrt()
}

/// Five-point stencil: error `h^4 f^(5) / 30`, about `1e-13` here.
fn stencil(f: impl Fn(f64) -> f64, t: f64) -> f64 {
    let h = 1e-3;
    (f(t - 2.0 * h) - 8.0 * f(t - h) + 8.0 * f(t + h) - f(t + 2.0 * h)) / (12.0 * h)
}

#[test]
fn dual_gradient_matches_a_finite_difference_of_the_plain_function() {
    let mut state = 0x2545_F491_4F6C_DD1D_u64;
    let mut coord = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let mag = 0.3 + (state >> 12) as f64 / 4_503_599_627_370_496.0 * 2.7;
        if state & 1 == 0 {
            mag
        } else {
            -mag
        }
    };
    for _ in 0..400 {
        let (x, y) = (coord(), coord());
        let g = composite(D2::variable(x, 0), D2::variable(y, 1));
        assert_eq!(g.v.to_bits(), composite(x, y).to_bits());
        let want = [
            stencil(|t| composite(t, y), x),
            stencil(|t| composite(x, t), y),
        ];
        for (got, want) in g.d.into_iter().zip(want) {
            let ok = (got - want).abs() <= 1e-9 * (1.0 + want.abs());
            assert!(ok, "at ({x}, {y}): got {got:e}, want {want:e}");
        }
    }
}

#[test]
fn dual_f32_gradient_tracks_f64() {
    for (x, y) in [(1.5, -0.75), (-2.0, 0.5), (0.25, 2.5)] {
        let g32 = composite(
            Dual::<f32, 2>::variable(x as f32, 0),
            Dual::variable(y as f32, 1),
        );
        let g64 = composite(D2::variable(x, 0), D2::variable(y, 1));
        for (got, want) in g32.d.into_iter().zip(g64.d) {
            let err = (f64::from(got) - want).abs();
            assert!(err <= 1e-5 * (1.0 + want.abs()), "at ({x}, {y})");
        }
    }
}

#[test]
fn abs_and_copysign_derivatives_follow_the_sign_conventions() {
    // `abs` and the copied-onto value of `copysign` use sgn(+-0) = +1.
    for (x, want) in [(2.0_f64, 1.0_f64), (-2.0, -1.0), (0.0, 1.0), (-0.0, 1.0)] {
        let r = D1::variable(x, 0).abs();
        assert_eq!(
            (r.v.to_bits(), r.d[0].to_bits()),
            (libm::fabs(x).to_bits(), want.to_bits())
        );
    }
    // `s` contributes its sign bit, so -0 flips.
    let cases: [(f64, f64, f64); 10] = [
        (2.0, 1.0, 1.0),
        (2.0, -1.0, -1.0),
        (-2.0, 1.0, -1.0),
        (-2.0, -1.0, 1.0),
        (0.0, 1.0, 1.0),
        (-0.0, 1.0, 1.0),
        (0.0, -1.0, -1.0),
        (2.0, -0.0, -1.0),
        (-2.0, -0.0, 1.0),
        (2.0, 0.0, 1.0),
    ];
    for (x, s, want) in cases {
        let r = D1::variable(x, 0).copysign(D1::constant(s));
        let got = (r.v.to_bits(), r.d[0].to_bits());
        assert_eq!(
            got,
            (libm::copysign(x, s).to_bits(), want.to_bits()),
            "({x}, {s})"
        );
    }
    // No dependence on the sign source.
    let r = D1::constant(2.0).copysign(D1::variable(-1.0, 0));
    assert_eq!(r.d[0].abs().to_bits(), 0);
}

#[test]
fn dual_sqrt_at_zero_has_an_infinite_derivative_and_the_plain_value() {
    let r = D2::variable(0.0, 0).sqrt();
    assert_eq!(r.v.to_bits(), 0.0_f64.to_bits());
    assert!(r.d[0].is_infinite() && r.d[0] > 0.0);
    assert!(r.d[1].is_nan(), "0 / 0");
    let n = D1::variable(-0.0, 0).sqrt();
    assert_eq!(n.v.to_bits(), (-0.0_f64).to_bits());
    assert!(n.d[0].is_infinite() && n.d[0] < 0.0);
    // The safe argument keeps the selected arm finite.
    let t2 = D1::variable(0.0, 0) * D1::variable(0.0, 0);
    let safe = D1::select(t2.lt(D1::lit(0.25)), D1::one(), t2);
    assert!(safe.sqrt().d[0].is_finite());
}

/// `d(n^2) = 2 v . dv = 0` in every lane at the zero vector, so the `sqrt` rule is `0 / 0`, and
/// a `sqrt` shared with a selected arm poisons every derivative built on it.
#[test]
fn dual_norm_of_the_zero_vector_is_nan_in_every_lane() {
    type D3 = Dual<f64, 3>;
    let z = Vector([
        D3::variable(0.0, 0),
        D3::variable(0.0, 1),
        D3::variable(0.0, 2),
    ]);
    let n = z.norm();
    assert_eq!(n.v.to_bits(), 0.0_f64.to_bits());
    assert!(n.d.iter().all(|d| d.is_nan()));
    let (s, c) = (n * D3::lit(0.5)).sin_cos();
    assert!(s.d.iter().chain(&c.d).all(|d| d.is_nan()));
}

/// `chol` of `diag(p, 1, 1)` and `diag(1, 1, p)` with every diagonal entry a variable: a failed
/// pivot reaches `sqrt` as the constant 1, so no value or derivative lane is non-finite.
fn chol_of_a_variable_diagonal_is_finite(p: usize, value: f64) {
    type D3 = Dual<f64, 3>;
    let o = D3::constant(0.0);
    let x: [D3; 3] = core::array::from_fn(|i| D3::variable(if i == p { value } else { 1.0 }, i));
    let (l, pd) = chol(&Matrix::from_cols([
        Vector([x[0], o, o]),
        Vector([o, x[1], o]),
        Vector([o, o, x[2]]),
    ]));
    assert!(!pd);
    for r in 0..3 {
        for k in 0..3 {
            let e = l.get(r, k);
            assert!(
                e.v.is_finite() && e.d.iter().all(|d| d.is_finite()),
                "L_{r}{k}"
            );
        }
    }
}

#[test]
fn a_zero_cholesky_pivot_leaves_every_dual_derivative_lane_finite() {
    chol_of_a_variable_diagonal_is_finite(0, 0.0);
}

#[test]
fn a_negative_cholesky_pivot_leaves_every_dual_derivative_lane_finite() {
    chol_of_a_variable_diagonal_is_finite(2, -1.0);
}

#[test]
fn select_and_branch_carry_the_derivative_of_the_arm_taken() {
    let (a, b) = (D2::variable(3.0, 0), D2::variable(0.5, 1));
    let (t, f) = (a * b, a / b);
    for m in [true, false] {
        let want = flat(if m { t } else { f });
        let (one_arm,) = D2::branch(m, || (t,), || (f,));
        for r in [D2::select(m, t, f), one_arm, f64::branch(m, || t, || f)] {
            assert_eq!(flat(r), want);
        }
    }
}

#[test]
fn dual_over_lanes_selects_per_lane() {
    type DL = Dual<L2, 1>;
    let x = DL::variable(L2([-2.0, 3.0]), 0);
    let lanes = |l: L2| l.0.map(f64::to_bits);
    let a = x.abs();
    assert_eq!(lanes(a.v), [2.0, 3.0].map(f64::to_bits));
    assert_eq!(lanes(a.d[0]), [-1.0, 1.0].map(f64::to_bits));
    let r = DL::branch(M2([true, false]), || x * x, || -x);
    assert_eq!(lanes(r.v), [4.0, -3.0].map(f64::to_bits));
    assert_eq!(lanes(r.d[0]), [-4.0, -1.0].map(f64::to_bits));
    let c = DL::variable(L2([-8.0, 27.0]), 0).cbrt();
    assert_eq!(lanes(c.v), [-2.0, 3.0].map(f64::to_bits));
    assert_eq!(lanes(c.d[0]), [1.0 / 12.0, 1.0 / 27.0].map(f64::to_bits));
}
