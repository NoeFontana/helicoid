//! S² and its chart (`docs/PHASE5.md` §2, `0063`, `0067`): the chart laws of `chart_tests` with a
//! direction's distance, since S² is not a group and has no `Log` to compare points through.
//!
//! Points are drawn as `laws::sample::<3>()` normalized, tangents as `laws::sample::<2>()`, so
//! `θ < √2 < π`: every case is inside `U₀`. The bounds are twice the worst of `measure_s2`.

use crate::laws::{self, e, unit, worst, Rng, Sample};
use crate::{reference, Chart, Jac, Manifold, RnTangent, S2Chart, Tangent, S2};
use core::array;
use helicoid_linalg::{Dual, Mat2, Matrix, Real, StridedMut, Vec3, Vector};
use proptest::prelude::*;
use std::vec::Vec;

type D2 = Dual<f64, 2>;

fn point<S: Sample>(v: &[f64; 3]) -> S2<S> {
    S2::from_vec_normalized(Vector(array::from_fn(|i| S::constant(v[i]))))
}

fn tangent<S: Sample>(v: &[f64; 2]) -> RnTangent<S, 2> {
    RnTangent {
        rho: Vector(array::from_fn(|i| S::constant(v[i]))),
    }
}

fn f(v: Vec3<impl Real>) -> [f64; 3] {
    v.0.map(Real::value_f64)
}

fn ft<S: Real>(t: &RnTangent<S, 2>) -> [f64; 2] {
    t.rho.0.map(S::value_f64)
}

fn fm<S: Real>(m: &Mat2<S>) -> [f64; 4] {
    let c = [m.col(0).0, m.col(1).0];
    [c[0][0], c[0][1], c[1][0], c[1][1]].map(S::value_f64)
}

fn bits<const N: usize>(a: [f64; N]) -> [u64; N] {
    a.map(f64::to_bits)
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The angle between two directions in `u`, `atan2(‖a × b‖, a·b)` (`NUMERICS.md` §11).
fn angle<S: Real>(a: &S2<S>, b: &S2<S>) -> f64 {
    let (a, b) = (f(a.vec()), f(b.vec()));
    let c = cross(a, b);
    dot(c, c).sqrt().atan2(dot(a, b)) / unit::<S>()
}

/// `‖x‖² - 1` in double-double: Veltkamp splits make every partial product exact without a fused
/// multiply-add (D16), and Knuth's two-sum carries the rest, so `5u` at binary64 is resolved.
fn excess(x: [f64; 3]) -> f64 {
    let split = |a: f64| {
        let c = 134_217_729.0 * a; // 2^27 + 1
        let hi = c - (c - a);
        (hi, a - hi)
    };
    let two_sum = |a: f64, b: f64| {
        let s = a + b;
        let bb = s - a;
        (s, (a - (s - bb)) + (b - bb))
    };
    let (mut hi, mut lo) = (-1.0_f64, 0.0_f64);
    for a in x {
        let (h, l) = split(a);
        for t in [h * h, 2.0 * h * l, l * l] {
            let (s, err) = two_sum(hi, t);
            hi = s;
            lo += err;
        }
    }
    hi + lo
}

/// CH.2's round trips at `n`, legs `retract(0) = n`, `local(retract(δ)) = δ`, `retract(local(m)) =
/// m`, `local(n) = 0`; points by [`angle`], tangents by `laws::e`.
fn legs<S: Sample>(a: &[f64; 3], b: &[f64; 2], c: &[f64; 3]) -> [f64; 4] {
    let (n, delta, m) = (point::<S>(a), tangent::<S>(b), point::<S>(c));
    let chart = S2Chart::at(&n);
    let zero = RnTangent::<S, 2>::zero();
    [
        angle(&chart.retract(&zero), &n),
        e::<S>(&ft(&chart.local(&chart.retract(&delta))), &ft(&delta)),
        angle(&chart.retract(&chart.local(&m)), &m),
        e::<S>(&ft(&chart.local(&n)), &[0.0; 2]),
    ]
}

/// `max(|BᵀB - I|, |Bᵀn|, ‖b₁ × b₂ - ςn‖)` in `u` (CH.8(c)), from the values of `S`.
fn basis_err<S: Sample>(a: &[f64; 3]) -> f64 {
    let n = point::<S>(a);
    let c = S2Chart::at(&n);
    let (b1, b2, nn) = (f(c.b1), f(c.b2), f(n.vec()));
    let sigma = if nn[2] < 0.0 { -1.0 } else { 1.0 };
    let k = cross(b1, b2);
    let errs = [
        dot(b1, b1) - 1.0,
        dot(b2, b2) - 1.0,
        dot(b1, b2),
        dot(b1, nn),
        dot(b2, nn),
        k[0] - sigma * nn[0],
        k[1] - sigma * nn[1],
        k[2] - sigma * nn[2],
    ];
    errs.iter().fold(0.0, |w, x| worst(w, x.abs())) / unit::<S>()
}

/// `|‖retract(δ)‖² - 1|` in `u`: the drift `retract` adds to a base on the sphere.
fn drift<S: Sample>(a: &[f64; 3], b: &[f64; 2]) -> f64 {
    let chart = S2Chart::at(&point::<S>(a));
    let m = chart.retract(&tangent::<S>(b)).vec().0.map(S::value_f64);
    excess(m).abs() / unit::<S>()
}

/// `S2Chart`'s two Jacobians against `Dual`, CH.1 read literally (`laws::chart_jacobians_match_dual`
/// with S²'s points): the frozen points are constants and only `η` carries lanes.
fn jacobians_match_dual(a: &[f64; 3], b: &[f64; 2], c: &[f64; 3]) -> f64 {
    let (n, delta, m) = (point::<D2>(a), tangent::<D2>(b), point::<D2>(c));
    let eta = RnTangent {
        rho: Vector([Dual::variable(0.0, 0), Dual::variable(0.0, 1)]),
    };
    let chart = S2Chart::at(&n);
    let lanes = |t: &RnTangent<D2, 2>| {
        let r = t.rho.0;
        [r[0].d[0], r[1].d[0], r[0].d[1], r[1].d[1]]
    };
    let y = chart.retract(&delta);
    let rj = lanes(&S2Chart::at(&y).local(&chart.retract(&delta.add(&eta))));
    let lj = lanes(&chart.local(&S2Chart::at(&m).retract(&eta)));
    worst(
        e::<f64>(&fm(&chart.retract_jacobian(&delta)), &rj),
        e::<f64>(&fm(&chart.local_jacobian(&m)), &lj),
    )
}

/// `local_jacobian(retract(δ)) · retract_jacobian(δ) = I` (CH.2(b)).
fn inverse_pair<S: Sample>(a: &[f64; 3], b: &[f64; 2]) -> f64 {
    let chart = S2Chart::at(&point::<S>(a));
    let delta = tangent::<S>(b);
    let p = chart
        .local_jacobian(&chart.retract(&delta))
        .mul(&chart.retract_jacobian(&delta));
    e::<S>(&fm(&p), &[1.0, 0.0, 0.0, 1.0])
}

/// `local` against `reference::s2_local`, relative to `max(‖δ‖, 1)` and divided by `local`'s
/// conditioning `α/sin α` (CH.11(e)): the twin halves the angle through `1 + n·m`, which near the
/// antipode errs by `u/(π - α)` times a constant 3 to 5 times the fast path's (CH.9(c)), so an
/// unscaled bound would be set by how close one draw came to `-n`.
fn twin<S: Sample>(a: &[f64; 3], c: &[f64; 3]) -> f64 {
    let (chart, m) = (S2Chart::at(&point::<S>(a)), point::<S>(c));
    let (n, mm) = (f(chart.base().vec()), f(m.vec()));
    let k = cross(n, mm);
    let (s, alpha) = (dot(k, k).sqrt(), dot(k, k).sqrt().atan2(dot(n, mm)));
    let cond = if s > 0.0 { alpha / s } else { 1.0 };
    e::<S>(&ft(&chart.local(&m)), &ft(&reference::s2_local(&chart, &m))) / cond
}

// `measure_s2`, 10^6 cases (10^5 for `Dual`): legs `f64` [0, 6.000, 16.521, 0], `f32`
// [0, 6.083, 17.051, 0]; basis `f64` 6.000, `f32` 5.786; drift 10.877, 8.766; inverse pair 8.834,
// 9.035; twin 5.870, 6.034; `Dual` 9.833.
const LEGS_F64: [f64; 4] = [0.0, 12.0, 34.0, 0.0];
const LEGS_F32: [f64; 4] = [0.0, 13.0, 35.0, 0.0];
const BASIS: [f64; 2] = [12.0, 12.0];
const DRIFT: [f64; 2] = [22.0, 18.0];
const INVERSE: [f64; 2] = [18.0, 19.0];
const TWIN: [f64; 2] = [12.0, 13.0];
const DUAL: f64 = 20.0;

proptest! {
    #[test]
    fn chart_legs(a in laws::sample::<3>(), b in laws::sample::<2>(), c in laws::sample::<3>()) {
        for (bounds, l) in [(LEGS_F64, legs::<f64>(&a, &b, &c)), (LEGS_F32, legs::<f32>(&a, &b, &c))] {
            for ((name, v), bound) in laws::CHART_LEGS.iter().zip(l).zip(bounds) {
                laws::within_leg(name, v, bound)?;
            }
        }
    }

    #[test]
    fn the_basis_is_orthonormal_and_right_handed_up_to_sign(a in laws::sample::<3>()) {
        laws::within(basis_err::<f64>(&a), BASIS[0])?;
        laws::within(basis_err::<f32>(&a), BASIS[1])?;
    }

    #[test]
    fn retract_stays_on_the_sphere(a in laws::sample::<3>(), b in laws::sample::<2>()) {
        laws::within(drift::<f64>(&a, &b), DRIFT[0])?;
        laws::within(drift::<f32>(&a, &b), DRIFT[1])?;
    }

    #[test]
    fn chart_jacobians_match_dual(
        a in laws::sample::<3>(), b in laws::sample::<2>(), c in laws::sample::<3>(),
    ) {
        laws::within(jacobians_match_dual(&a, &b, &c), DUAL)?;
    }

    #[test]
    fn local_jacobian_inverts_retract_jacobian(a in laws::sample::<3>(), b in laws::sample::<2>()) {
        laws::within(inverse_pair::<f64>(&a, &b), INVERSE[0])?;
        laws::within(inverse_pair::<f32>(&a, &b), INVERSE[1])?;
    }

    #[test]
    fn s2_local_matches_reference(a in laws::sample::<3>(), c in laws::sample::<3>()) {
        laws::within(twin::<f64>(&a, &c), TWIN[0])?;
        laws::within(twin::<f32>(&a, &c), TWIN[1])?;
    }

    /// `local_with_jacobian` is the two calls to the bit (`0060` decision 3).
    #[test]
    fn local_with_jacobian_is_the_two_calls(a in laws::sample::<3>(), c in laws::sample::<3>()) {
        let (chart, m) = (S2Chart::at(&point::<f64>(&a)), point::<f64>(&c));
        let (t, j) = chart.local_with_jacobian(&m);
        prop_assert_eq!(ft(&t).map(f64::to_bits), ft(&chart.local(&m)).map(f64::to_bits));
        prop_assert_eq!(fm(&j).map(f64::to_bits), fm(&chart.local_jacobian(&m)).map(f64::to_bits));
    }

    /// Every chart method under `Dual<f64, 2>` has the value of the plain `f64` one, bit for bit:
    /// the lanes ride along and change no branch.
    #[test]
    fn dual_value_is_plain_value(
        a in laws::sample::<3>(), b in laws::sample::<2>(), c in laws::sample::<3>(),
    ) {
        let plain = |n: S2<f64>, d: RnTangent<f64, 2>, m: S2<f64>| {
            let ch = S2Chart::at(&n);
            let mut v: Vec<f64> = f(ch.retract(&d).vec()).to_vec();
            v.extend(ft(&ch.local(&m)));
            v.extend(fm(&ch.retract_jacobian(&d)));
            v.extend(fm(&ch.local_jacobian(&m)));
            v
        };
        let lifted = |n: S2<D2>, d: RnTangent<D2, 2>, m: S2<D2>| {
            let ch = S2Chart::at(&n);
            let mut v: Vec<f64> = f(ch.retract(&d).vec()).to_vec();
            v.extend(ft(&ch.local(&m)));
            v.extend(fm(&ch.retract_jacobian(&d)));
            v.extend(fm(&ch.local_jacobian(&m)));
            v
        };
        let (n, m) = (point::<f64>(&a), point::<f64>(&c));
        let var = |p: S2<f64>| S2::from_vec_unchecked(Vector(array::from_fn(|i| {
            Dual::variable(p.vec().0[i], i % 2)
        })));
        let dv = RnTangent { rho: Vector([Dual::variable(b[0], 0), Dual::variable(b[1], 1)]) };
        let want = plain(n, tangent::<f64>(&b), m);
        let got = lifted(var(n), dv, var(m));
        prop_assert_eq!(want.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                        got.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
    }
}

/// The worst of each law over seeded cases, the figures the bounds above are recorded from.
#[test]
#[ignore = "measurement: prints the figures the bounds are recorded from"]
#[allow(clippy::print_stdout)]
fn measure_s2() {
    let mut rng = Rng(0x7332_5f63_6861_7274);
    let (mut l64, mut l32) = ([0.0_f64; 4], [0.0_f64; 4]);
    let mut w = [0.0_f64; 8];
    for _ in 0..1_000_000 {
        let (a, b, c) = (rng.shaped::<3>(), rng.shaped::<2>(), rng.shaped::<3>());
        for (acc, l) in [
            (&mut l64, legs::<f64>(&a, &b, &c)),
            (&mut l32, legs::<f32>(&a, &b, &c)),
        ] {
            for (x, y) in acc.iter_mut().zip(l) {
                *x = worst(*x, y);
            }
        }
        let v = [
            basis_err::<f64>(&a),
            basis_err::<f32>(&a),
            drift::<f64>(&a, &b),
            drift::<f32>(&a, &b),
            inverse_pair::<f64>(&a, &b),
            inverse_pair::<f32>(&a, &b),
            twin::<f64>(&a, &c),
            twin::<f32>(&a, &c),
        ];
        for (x, y) in w.iter_mut().zip(v) {
            *x = worst(*x, y);
        }
    }
    let mut dual = 0.0_f64;
    for _ in 0..100_000 {
        let (a, b, c) = (rng.shaped::<3>(), rng.shaped::<2>(), rng.shaped::<3>());
        dual = worst(dual, jacobians_match_dual(&a, &b, &c));
    }
    std::println!(
        "legs f64 {l64:.3?} f32 {l32:.3?}; basis {:.3} {:.3}; drift {:.3} {:.3}; inverse pair {:.3} {:.3}; twin {:.3} {:.3}; dual {dual:.3}",
        w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7]
    );
}

/// `ς` is read by comparison (`0063` decision 1): at `n_z = +0` and `-0` the basis, and so every
/// chart output, is the same bits.
#[test]
fn basis_ignores_the_sign_of_zero() {
    fn case<S: Sample>(lambda: f64) {
        let at = |z: S| {
            let n = S2::from_vec_unchecked(Vector([
                S::constant(lambda.cos()),
                S::constant(lambda.sin()),
                z,
            ]));
            let c = S2Chart::at(&n);
            let d = tangent::<S>(&[0.3, -0.7]);
            let m = point::<S>(&[0.2, -0.5, 0.4]);
            let mut v = f(c.b1).to_vec();
            v.extend(f(c.b2));
            v.extend(f(c.retract(&d).vec()));
            v.extend(ft(&c.local(&m)));
            v.extend(fm(&c.retract_jacobian(&d)));
            v.extend(fm(&c.local_jacobian(&m)));
            v.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        };
        assert_eq!(at(S::zero()), at(-S::zero()), "lambda {lambda}");
    }
    for lambda in [0.0, 0.3, 1.9, 4.0, -2.5] {
        case::<f64>(lambda);
        case::<f32>(lambda);
    }
}

/// `0063`'s measurement, as a test: one step from `| ‖n‖² - 1 | = 2^ℓ U(½, 1)`, `ℓ` up to the band's
/// edge, lands within about `5u` of the sphere at both precisions. Measured `5.005 u` (`f64`) and
/// `4.447 u` (`f32`) over these draws: `0063`'s table gives the `f64` worst as `5`, rounded, and
/// the step's own truncation `¾η²` is already `1.0 u` at the edge (SO.14). The bound is the
/// measurement rounded up to the half `u`.
#[test]
fn renormalize_lands_within_five_u_in_the_band() {
    fn worst_step<S: Sample>(edge: f64, lowest: i32) -> f64 {
        let mut rng = Rng(0x7332_5f72_656e_6f72);
        let mut worst_u = 0.0_f64;
        for i in 0..20_000 {
            let l = f64::from(lowest) + (edge - f64::from(lowest)) * (0.5 + 0.5 * rng.unif());
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            let eta = (0.75 + 0.25 * rng.unif()) * 2_f64.powf(l) * sign;
            let dir = rng.arr::<3>();
            let r = (1.0 + eta).sqrt() / dot(dir, dir).sqrt();
            let v = Vector(dir.map(|x| S::constant(x * r)));
            if excess(v.0.map(S::value_f64)).abs() > 2_f64.powf(edge) {
                continue;
            }
            // A drifted point, as repeated retracts leave one: past what `from_vec_unchecked`
            // vouches.
            let mut n = S2(v);
            n.renormalize();
            worst_u = worst_u.max(excess(n.vec().0.map(S::value_f64)).abs() / unit::<S>());
        }
        worst_u
    }
    let (w64, w32) = (
        worst_step::<f64>(-26.29, -40),
        worst_step::<f32>(-11.79, -16),
    );
    assert!(w64 <= 5.5 && w32 <= 5.5, "f64 {w64} u, f32 {w32} u");
}

/// `from_vec_unchecked` vouches to `2^-40` (`f64`): a vector beyond it fails the assert.
#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "S2::from_vec_unchecked")]
fn from_vec_unchecked_asserts_the_unit_norm() {
    let s = (1.0 + 2_f64.powi(-38)).sqrt();
    let _ = S2::from_vec_unchecked(Vector([s, 0.0, 0.0]));
}

/// Positive control: inside the bound, and NaN fails it (`unit_tolerance`'s reading).
#[test]
fn from_vec_unchecked_accepts_its_bound() {
    let s = (1.0 + 2_f64.powi(-42)).sqrt();
    let v = S2::from_vec_unchecked(Vector([s, 0.0, 0.0])).vec().0;
    assert_eq!(bits(v), bits([s, 0.0, 0.0]));
    let nan = std::panic::catch_unwind(|| S2::from_vec_unchecked(Vector([f64::NAN, 0.0, 0.0])));
    assert_eq!(nan.is_err(), cfg!(debug_assertions));
}

#[test]
fn from_vec_normalized_divides_by_the_norm() {
    let n = S2::from_vec_normalized(Vector([3.0_f64, 0.0, 4.0]));
    assert_eq!(bits(n.vec().0), bits([0.6, 0.0, 0.8]));
    let n = S2::from_vec_normalized(Vector([0.0_f32, -2.0, 0.0]));
    assert_eq!(
        n.vec().0.map(f32::to_bits),
        [0.0_f32, -1.0, 0.0].map(f32::to_bits)
    );
}

/// At `m = -n`, outside the domain, `local` is `0` and its Jacobian finite (`NUMERICS.md` §12);
/// closer to `-n` than `‖n × m‖²` can hold, finite too.
#[test]
fn local_is_finite_at_and_beyond_the_antipode() {
    for (n, m) in [
        ([0.0, 0.0, 1.0], [0.0, 0.0, -1.0]),
        ([0.6, 0.0, -0.8], [-0.6, 0.0, 0.8]),
        ([0.0, 0.0, 1.0], [1e-200, 0.0, -1.0]),
    ] {
        let chart = S2Chart::at(&S2::from_vec_unchecked(Vector(n)));
        let (t, j) = chart.local_with_jacobian(&S2::from_vec_unchecked(Vector(m)));
        assert!(
            ft(&t).iter().chain(&fm(&j)).all(|x| x.is_finite()),
            "{n:?} {m:?}"
        );
    }
    let chart = S2Chart::at(&S2::from_vec_unchecked(Vector([0.0_f64, 0.0, 1.0])));
    let t = chart.local(&S2::from_vec_unchecked(Vector([0.0, 0.0, -1.0])));
    assert_eq!(bits(ft(&t).map(f64::abs)), [0; 2]);
}

/// At `n = e_z` the basis is `(e_x, e_y)` exactly (`ν = 2e_z`), so `retract` and `local` are the
/// textbook ones: `δ = (θ, 0)` rotates about `e_x`, taking `e_z` towards `-e_y`.
#[test]
fn the_north_pole_is_the_textbook_chart() {
    let chart = S2Chart::at(&S2::from_vec_unchecked(Vector([0.0_f64, 0.0, 1.0])));
    // `-(kx * y)` at `x = y = 0` is `-0`: the off-diagonal zeros are negative.
    assert_eq!(bits(f(chart.b1)), bits([1.0, -0.0, -0.0]));
    assert_eq!(bits(f(chart.b2)), bits([-0.0, 1.0, -0.0]));
    let t = 0.25_f64;
    let m = chart
        .retract(&RnTangent {
            rho: Vector([t, 0.0]),
        })
        .vec()
        .0;
    assert!(
        (m[0]).abs() == 0.0 && (m[1] + t.sin()).abs() < 1e-16 && (m[2] - t.cos()).abs() < 1e-16
    );
    let back = chart.local(&S2::from_vec_unchecked(Vector(m)));
    assert!((ft(&back)[0] - t).abs() < 1e-16 && ft(&back)[1].abs() < 1e-16);
}

/// `Mat2` as a `Jac` on `RnTangent<S, 2>`: integer entries make every product exact, so each
/// method is checked bit for bit against its dense definition.
#[test]
fn mat2_jac_is_the_dense_matrix() {
    let a: Mat2<f64> = Matrix::from_rows([Vector([2.0, -3.0]), Vector([5.0, 7.0])]);
    let b: Mat2<f64> = Matrix::from_rows([Vector([-1.0, 4.0]), Vector([6.0, 0.5])]);
    let t = RnTangent {
        rho: Vector([3.0, -2.0]),
    };
    assert_eq!(bits(a.apply(&t).rho.0), bits([12.0, 1.0]));
    assert_eq!(bits(a.apply_transpose(&t).rho.0), bits([-4.0, -23.0]));
    assert_eq!(bits(fm(&a.mul(&b))), bits([-20.0, 37.0, 6.5, 23.5]));
    let neg = Jac::<f64, RnTangent<f64, 2>>::neg(&a);
    assert_eq!(bits(fm(&neg)), bits([-2.0, -5.0, 3.0, -7.0]));
    let id = <Mat2<f64> as Jac<f64, RnTangent<f64, 2>>>::identity();
    assert_eq!(bits(fm(&id)), bits([1.0, 0.0, 0.0, 1.0]));
    // det = 1: the adjugate is the inverse, exactly.
    let u: Mat2<f64> = Matrix::from_rows([Vector([2.0, 3.0]), Vector([1.0, 2.0])]);
    let inv = Jac::<f64, RnTangent<f64, 2>>::inverse(&u);
    assert_eq!(bits(fm(&inv)), bits([2.0, -1.0, -3.0, 2.0]));
    let mut out = [0.0; 4];
    Jac::<f64, RnTangent<f64, 2>>::write_dense(&a, &mut StridedMut::col_major(&mut out, 2, 2));
    assert_eq!(bits(out), bits([2.0, 5.0, -3.0, 7.0]));
    // `J Σ Jᵀ` with an asymmetric `Σ`, so a transposed sandwich would not pass.
    let cov: Matrix<f64, 2, 2> = Matrix::from_rows([Vector([1.0, 2.0]), Vector([-1.0, 3.0])]);
    let s = a.sandwich::<2>(&cov);
    let want = a * cov * a.transpose();
    assert_eq!(bits(fm(&s)), bits(fm(&want)));
}

#[test]
fn s2_is_a_manifold_of_dimension_two() {
    assert_eq!(<S2<f64> as Manifold<f64>>::DOF, 2);
    assert_eq!(<RnTangent<f64, 2> as Tangent<f64>>::DOF, 2);
}
