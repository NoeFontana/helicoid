//! Sim(3) under the group laws, its structured Jacobian against the dense twins, its reduction to
//! SE(3) at `σ = 0`, and `Exp`/`Log` to `ln MAX` (`NUMERICS.md` §9, `0070`).
//!
//! Samples are `laws::sample`'s, so `|σ| < 1` and `θ < √3`; the bounds are twice the worst of
//! `measure_sim3`, rounded up.

use crate::laws::{self, tangent as law_tangent, Bounds, Rng, Sample};
use crate::{reference, Jac, LieGroup, SEn3Tangent, Sim3, Sim3Jac, Sim3Tangent, Tangent, SE3};
use helicoid_linalg::{Dual, Matrix, StridedMut, Vector};
use proptest::prelude::*;

// `measure_sim3` (2·10⁵ draws, 2·10⁴ for `Dual`) and `measure_geodesic` (10⁶), worst of each law:
// `f64` axioms 6.021, exp_log 4.556, adjoint 8.546, jl_ad_jr 14.695, plus_minus 9.879, rows 0
// (exact), ad 10.375, sides 1.089, tangent_order 0, jac_order 3.316, sandwich 4.431; `Dual` within
// those; `dual_rows` 93.818; geodesic legs 1.057, 6.708, 33.426, 15.376, 10.720, 10.791, 1.118.
// `dual_rows` is the `Dual` derivative through the closed arm near `r_×`, where it divides by `Δ`
// once more than the value (SM.7(f)); it read 6720 before `v₀` took its series below `σ_w`.
const F64: Bounds = Bounds {
    axioms: 13.0,
    exp_log: 10.0,
    adjoint: 18.0,
    jl_ad_jr: 30.0,
    plus_minus: 20.0,
    rows: 0.0,
    ad: 21.0,
    sides: 3.0,
    tangent_order: 0.0,
    jac_order: 7.0,
    sandwich: 9.0,
    dual_rows: 188.0,
    geodesic: [3.0, 14.0, 67.0, 31.0, 22.0, 22.0, 3.0],
};
// `f32`: axioms 6.191, exp_log 3.634, adjoint 8.424, jl_ad_jr 7.735, plus_minus 11.277, rows 0,
// ad 7.368, sides 1.079, tangent_order 3.284, jac_order 3.466, sandwich 3.342; geodesic legs
// 1.048, 6.230, 33.755, 9.337, 9.379, 10.916, 1.031.
const F32: Bounds = Bounds {
    axioms: 13.0,
    exp_log: 8.0,
    adjoint: 17.0,
    jl_ad_jr: 16.0,
    plus_minus: 23.0,
    ad: 15.0,
    tangent_order: 7.0,
    sandwich: 7.0,
    geodesic: [3.0, 13.0, 68.0, 19.0, 19.0, 22.0, 3.0],
    ..F64
};

/// `J_r(τ)` at the sample: every block nonzero, so `apply_transpose` and `f` are exercised.
fn jac<S: Sample>(v: &[f64; 7]) -> Sim3Jac<S> {
    <Sim3<S> as LieGroup<S>>::jr(&law_tangent::<S, Sim3<S>, 7>(v))
}

mod group {
    use super::{jac, F32, F64};
    use crate::laws::laws_for;
    use crate::Sim3;
    laws_for!(sim3, Sim3, jac, 7, F64, F32);
}

fn t<S: Sample>(v: &[f64; 7]) -> Sim3Tangent<S> {
    law_tangent::<S, Sim3<S>, 7>(v)
}

/// The worst of every law over seeded draws, the figures `F64`/`F32` are recorded from.
#[test]
#[ignore = "measurement: prints the figures the bounds are recorded from"]
#[allow(clippy::print_stdout)]
fn measure_sim3() {
    fn run<S: Sample>(n: usize) -> [f64; 11] {
        let mut rng = Rng(0x7369_6d33_6c61_7773);
        let mut w = [0.0_f64; 11];
        for _ in 0..n {
            let (a, b, c) = (rng.shaped::<7>(), rng.shaped::<7>(), rng.shaped::<7>());
            let g = |v: &[f64; 7]| Sim3::<S>::exp(&t::<S>(v));
            let k = S::sample(c[0], 0);
            let v = [
                laws::group_axioms::<S, Sim3<S>, 7>(g(&a), g(&b), g(&c)),
                laws::exp_log_roundtrip::<S, Sim3<S>, 7>(&t(&a), &g(&b)),
                laws::adjoint_identity::<S, Sim3<S>, 7>(&g(&a), &t(&b)),
                laws::jl_is_ad_jr::<S, Sim3<S>, 7>(&t(&a)),
                laws::plus_minus::<S, Sim3<S>, 7>(&g(&a), &g(&b), &t(&c)),
                laws::jacobian_rows::<S, Sim3<S>, 7>(&g(&a), &g(&b), &t(&c)),
                laws::ad_consistency::<S, Sim3<S>, 7>(&t(&a), &t(&b), &t(&c)),
                laws::side_delegation::<S, Sim3<S>, 7>(&g(&a), &g(&b), &t(&c)),
                laws::tangent_dense_order::<S, Sim3<S>, 7>(&t(&a), &t(&b), k),
                laws::jac_dense_order::<S, Sim3<S>, 7>(&jac(&a), &jac(&b), &t(&c)),
                laws::sandwich_matches_dense::<S, Sim3<S>, 7>(
                    &jac(&a),
                    &laws::cov::<S, 7>(&a, &b, &c),
                ),
            ];
            for (x, y) in w.iter_mut().zip(v) {
                *x = laws::worst(*x, y);
            }
        }
        w
    }
    let names = [
        "axioms",
        "exp_log",
        "adjoint",
        "jl_ad_jr",
        "plus_minus",
        "rows",
        "ad",
        "sides",
        "tangent_order",
        "jac_order",
        "sandwich",
    ];
    for (label, w) in [
        ("f64", run::<f64>(200_000)),
        ("f32", run::<f32>(200_000)),
        ("dual", run::<Dual<f64, 7>>(20_000)),
    ] {
        std::print!("{label}:");
        for (n, v) in names.iter().zip(w) {
            std::print!(" {n} {v:.3}");
        }
        std::println!();
    }
    let mut rng = Rng(0x7369_6d33_6475_616c);
    let mut d = 0.0_f64;
    for _ in 0..20_000 {
        let (a, b, c) = (rng.shaped::<7>(), rng.shaped::<7>(), rng.shaped::<7>());
        d = laws::worst(
            d,
            laws::jacobians_match_dual::<f64, Sim3<Dual<f64, 7>>, 7>(&a, &b, &c),
        );
    }
    std::println!("dual_rows {d:.3}");
}

/// The dense matrix of `j`, row-major.
fn dense(j: &Sim3Jac<f64>) -> [[f64; 7]; 7] {
    let mut m = [[0.0; 7]; 7];
    j.write_dense(&mut StridedMut::with_strides(
        m.as_flattened_mut(),
        7,
        7,
        7,
        1,
    ));
    m
}

/// `‖a - b‖_F / max(‖b‖_F, 1)` in `u`.
fn err(a: &[f64], b: &[f64]) -> f64 {
    laws::e::<f64>(a, b)
}

// `measure_twins`, 10⁵ draws: mul 1.281, inverse 3.165, apply 2.983, apply_transpose 2.004,
// sandwich 4.150, jr_inv 2.683.
const TWINS: [f64; 6] = [3.0, 7.0, 6.0, 5.0, 9.0, 6.0];

/// Each `Sim3Jac` operation and `jr_inv`'s block algebra against its dense twin
/// (`NUMERICS.md` §14), `[mul, inverse, apply, apply_transpose, sandwich, jr_inv]` in `u`.
fn twins(a: &[f64; 7], b: &[f64; 7], c: &[f64; 7]) -> [f64; 6] {
    let (ja, jb) = (jac::<f64>(a), <Sim3<f64>>::exp(&t(b)).adjoint());
    let x = t::<f64>(c);
    let mut out = [0.0; 49];
    reference::sim3jac_mul(
        &ja,
        &jb,
        &mut StridedMut::with_strides(&mut out, 7, 7, 7, 1),
    );
    let mul = err(dense(&ja.mul(&jb)).as_flattened(), &out);
    reference::sim3jac_inverse(&ja, &mut StridedMut::with_strides(&mut out, 7, 7, 7, 1));
    let inv = err(dense(&ja.inverse()).as_flattened(), &out);
    let mut xs = [0.0; 7];
    x.write_dense(&mut xs);
    let (mut y, mut d) = ([0.0; 7], [0.0; 7]);
    reference::sim3jac_apply(&ja, &xs, &mut y);
    ja.apply(&x).write_dense(&mut d);
    let apply = err(&d, &y);
    reference::sim3jac_apply_transpose(&ja, &xs, &mut y);
    ja.apply_transpose(&x).write_dense(&mut d);
    let apply_t = err(&d, &y);
    let cov: Matrix<f64, 7, 7> = laws::cov::<f64, 7>(a, b, c);
    reference::sim3jac_sandwich(
        &ja,
        &cov,
        &mut StridedMut::with_strides(&mut out, 7, 7, 7, 1),
    );
    let s = ja.sandwich::<7>(&cov);
    let sw: [f64; 49] = core::array::from_fn(|i| s.get(i / 7, i % 7));
    let sandwich = err(&sw, &out);
    let jr = <Sim3<f64>>::jr(&t(a));
    reference::sim3jac_inverse(&jr, &mut StridedMut::with_strides(&mut out, 7, 7, 7, 1));
    let jr_inv = err(dense(&<Sim3<f64>>::jr_inv(&t(a))).as_flattened(), &out);
    [mul, inv, apply, apply_t, sandwich, jr_inv]
}

proptest! {
    #[test]
    fn sim3jac_matches_reference(
        a in laws::sample::<7>(), b in laws::sample::<7>(), c in laws::sample::<7>(),
    ) {
        for (v, bound) in twins(&a, &b, &c).into_iter().zip(TWINS) {
            laws::within(v, bound)?;
        }
    }

    /// At `σ = 0` Sim(3) is SE(3) with a scale of one: `Exp`'s rotation and translation and the
    /// upper-left `6 x 6` of `J_r` are SE(3)'s within a few `u` (SM.10(b)).
    #[test]
    fn sigma_zero_is_se3(a in laws::sample::<7>()) {
        let mut v = a;
        v[6] = 0.0;
        let x = <Sim3<f64>>::exp(&t(&v));
        let se = SEn3Tangent { phi: Vector([v[0], v[1], v[2]]), rho: [Vector([v[3], v[4], v[5]])] };
        let y = SE3::exp(&se);
        let (_, tx, _) = x.parts();
        laws::within(err(&tx.0, &y.translation().0), SIGMA_ZERO)?;
        let js = dense(&<Sim3<f64>>::jr(&t(&v)));
        let je = reference::dense::<f64, _, _, 6>(&SE3::jr(&se));
        let a6: [f64; 36] = core::array::from_fn(|i| js[i / 6][i % 6]);
        laws::within(err(&a6, je.as_flattened()), SIGMA_ZERO)?;
    }
}

// 3.408 over 10⁵ draws.
const SIGMA_ZERO: f64 = 7.0;

#[test]
#[ignore = "measurement: prints the figures `TWINS` is recorded from"]
#[allow(clippy::print_stdout)]
fn measure_twins() {
    let mut rng = Rng(0x7369_6d33_7477_696e);
    let mut w = [0.0_f64; 6];
    for _ in 0..100_000 {
        let (a, b, c) = (rng.shaped::<7>(), rng.shaped::<7>(), rng.shaped::<7>());
        for (x, y) in w.iter_mut().zip(twins(&a, &b, &c)) {
            *x = laws::worst(*x, y);
        }
    }
    std::println!("twins {w:.3?}");
}

/// `Log(Exp τ) = τ` to `ln MAX`, where the unscaled closed forms would overflow from `σ ≈ 242`
/// (`Log`'s `Λ`) and `703` (`v₁`, `v₂`) and fail silently (SM.7(g)): every value finite and the
/// round trip relative to `max(‖τ‖, 1)` within `LARGE` `u`.
#[test]
fn exp_log_hold_to_ln_max() {
    for sigma in [10.0, -10.0, 100.0, -100.0, 300.0, -300.0, 700.0, -700.0] {
        for theta in [1e-4, 1.0, 3.0] {
            let tau = Sim3Tangent {
                phi: Vector([theta * 0.6, 0.0, theta * 0.8]),
                rho: Vector([0.3, -1.2, 0.5]),
                sigma,
            };
            let back = <Sim3<f64>>::exp(&tau).log();
            let (mut a, mut b) = ([0.0; 7], [0.0; 7]);
            back.write_dense(&mut a);
            tau.write_dense(&mut b);
            let e = err(&a, &b);
            assert!(
                e.is_finite() && e <= LARGE,
                "σ = {sigma}, θ = {theta}: {e} u"
            );
        }
    }
}

// At most 0.580 (`σ = 10`, `θ = 3`).
const LARGE: f64 = 2.0;
