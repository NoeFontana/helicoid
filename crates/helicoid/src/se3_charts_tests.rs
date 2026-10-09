//! The SE(3) charts under the chart laws (`docs/PHASE5.md` §1.3–§1.4, `0060`), their transitions
//! (CH.4(b)) and changes of chart (CH.6), and `TwistBlockJac` against the pair it wraps.
//!
//! The bounds are recorded as described at [`ChartBounds`], from `measure_charts`.

use crate::chart_tests::chart_laws_for;
use crate::laws::{self, ChartBounds, CHART_LEGS};
use crate::{
    Chart, Decoupled, Jac, LieGroup, ProductJac, SO3Tangent, Screw, Se3Chart, Tangent, Twist,
    TwistBlockJac, WorldTranslation, SE3,
};
use helicoid_linalg::{Mat3, Matrix, Real, Vector};

type G<S> = SE3<S>;
type Dl = helicoid_linalg::Dual<f64, 6>;

/// `SE3::chart_transition::<F, T>` as a function of the base, the closed form the `Dual`
/// transition is checked against.
fn tr<S: Real, F: Se3Chart<S>, T: Se3Chart<S>>(x: &SE3<S>) -> TwistBlockJac<S> {
    x.chart_transition::<F, T>()
}

// `decoupled`, worst of `measure_charts`: legs `f64` [1.116, 6.06, 8.015, 1.118], `f32` [1.068, 6.103, 8.17, 1.118]; `Dual` 5.908, transition 1.658,
// change of chart 7.561.
const DECOUPLED_F64: ChartBounds = ChartBounds {
    legs: [3.0, 13.0, 17.0, 3.0],
    dual: 12.0,
    transition: 4.0,
    change: 16.0,
};
const DECOUPLED_F32: ChartBounds = ChartBounds {
    legs: [3.0, 13.0, 17.0, 3.0],
    ..DECOUPLED_F64
};

// `screw`, worst of `measure_charts`: legs `f64` [1.116, 6.212, 9.106, 1.118], `f32` [1.068, 6.183, 7.984, 1.118]; `Dual` 4.849, transition 2.219,
// change of chart 9.906.
const SCREW_F64: ChartBounds = ChartBounds {
    legs: [3.0, 13.0, 19.0, 3.0],
    dual: 10.0,
    transition: 5.0,
    change: 20.0,
};
const SCREW_F32: ChartBounds = ChartBounds {
    legs: [3.0, 13.0, 16.0, 3.0],
    ..SCREW_F64
};

// `world_translation`, worst of `measure_charts`: legs `f64` [1.116, 6.06, 8.013, 1.118], `f32` [1.068, 6.102, 7.794, 1.118]; `Dual` 3.489, transition 1.658,
// change of chart 7.543.
const WORLD_TRANSLATION_F64: ChartBounds = ChartBounds {
    legs: [3.0, 13.0, 17.0, 3.0],
    dual: 7.0,
    transition: 4.0,
    change: 16.0,
};
const WORLD_TRANSLATION_F32: ChartBounds = ChartBounds {
    legs: [3.0, 13.0, 16.0, 3.0],
    ..WORLD_TRANSLATION_F64
};

chart_laws_for!(
    decoupled,
    G,
    Decoupled,
    WorldTranslation,
    tr::<Dl, Decoupled<Dl>, WorldTranslation<Dl>>,
    6,
    DECOUPLED_F64,
    DECOUPLED_F32
);
chart_laws_for!(
    world_translation,
    G,
    WorldTranslation,
    Decoupled,
    tr::<Dl, WorldTranslation<Dl>, Decoupled<Dl>>,
    6,
    WORLD_TRANSLATION_F64,
    WORLD_TRANSLATION_F32
);
chart_laws_for!(
    screw,
    G,
    Screw,
    Decoupled,
    tr::<Dl, Screw<Dl>, Decoupled<Dl>>,
    6,
    SCREW_F64,
    SCREW_F32
);

/// `DΦ(0)` of every ordered pair of the three charts against `chart_transition` (CH.4(b)):
/// `I` within the body frame and between a chart and itself, `diag(I, R)` into the world frame,
/// `diag(I, Rᵀ)` out of it.
macro_rules! every_pair {
    ($($name:ident: $F:ident -> $T:ident;)*) => {
        proptest::proptest! {$(
            #[test]
            fn $name(a in laws::sample::<6>()) {
                let v = laws::transition_matches::<f64, G<Dl>, $F<Dl>, $T<Dl>, _, 6>(
                    &a, tr::<Dl, $F<Dl>, $T<Dl>>);
                laws::within(v, TRANSITION)?;
            }
        )*}
    };
}

// Twice the worst of `measure_pairs` (10^5 cases, all nine pairs), rounded up: 2.273.
const TRANSITION: f64 = 5.0;

every_pair! {
    chart_transition_matches_dual_screw_screw: Screw -> Screw;
    chart_transition_matches_dual_screw_decoupled: Screw -> Decoupled;
    chart_transition_matches_dual_screw_world: Screw -> WorldTranslation;
    chart_transition_matches_dual_decoupled_screw: Decoupled -> Screw;
    chart_transition_matches_dual_decoupled_decoupled: Decoupled -> Decoupled;
    chart_transition_matches_dual_decoupled_world: Decoupled -> WorldTranslation;
    chart_transition_matches_dual_world_screw: WorldTranslation -> Screw;
    chart_transition_matches_dual_world_decoupled: WorldTranslation -> Decoupled;
    chart_transition_matches_dual_world_world: WorldTranslation -> WorldTranslation;
}

// The third pair CH.6 names; the other two are `decoupled` and `screw` above. Twice the worst of
// `measure_pairs`, rounded up: 4.819.
const CHANGE_SCREW_WORLD: f64 = 10.0;

proptest::proptest! {
    #[test]
    fn change_of_chart_screw_world_translation(a in laws::sample::<6>(), b in laws::sample::<6>()) {
        let v = laws::change_of_chart::<f64, G<Dl>, Screw<Dl>, WorldTranslation<Dl>, 6>(&a, &b);
        laws::within(v, CHANGE_SCREW_WORLD)?;
    }
}

/// The worst `DΦ(0)` error per pair, and the `Screw → WorldTranslation` change of chart, over
/// 10^5 seeded cases: the figures `TRANSITION` and `CHANGE_SCREW_WORLD` are recorded from.
#[test]
#[ignore = "measurement: prints the figures the bounds are recorded from"]
#[allow(clippy::print_stdout)]
fn measure_pairs() {
    fn pair<F: Se3Chart<Dl>, T: Se3Chart<Dl>>(a: &[f64; 6]) -> f64 {
        laws::transition_matches::<f64, G<Dl>, F, T, _, 6>(a, tr::<Dl, F, T>)
    }
    let mut rng = laws::Rng(0x6368_6172_7473_0003);
    let (mut worst, mut change) = (0.0_f64, 0.0_f64);
    for _ in 0..100_000 {
        let (a, b) = (rng.shaped::<6>(), rng.shaped::<6>());
        for v in [
            pair::<Screw<Dl>, Screw<Dl>>(&a),
            pair::<Screw<Dl>, Decoupled<Dl>>(&a),
            pair::<Screw<Dl>, WorldTranslation<Dl>>(&a),
            pair::<Decoupled<Dl>, Screw<Dl>>(&a),
            pair::<Decoupled<Dl>, Decoupled<Dl>>(&a),
            pair::<Decoupled<Dl>, WorldTranslation<Dl>>(&a),
            pair::<WorldTranslation<Dl>, Screw<Dl>>(&a),
            pair::<WorldTranslation<Dl>, Decoupled<Dl>>(&a),
            pair::<WorldTranslation<Dl>, WorldTranslation<Dl>>(&a),
        ] {
            worst = laws::worst(worst, v);
        }
        change = laws::worst(
            change,
            laws::change_of_chart::<f64, G<Dl>, Screw<Dl>, WorldTranslation<Dl>, 6>(&a, &b),
        );
    }
    std::println!("transition {worst:.3} change screw->world {change:.3}");
}

/// `Screw → Decoupled` is the identity to the bit, at any base (CH.4(b): `D(J_l(φ)ρ)/Dρ = I` at
/// `φ = 0`), and so is a chart to itself.
#[test]
fn chart_transition_screw_decoupled_is_identity() {
    let id = laws::dense_bits::<Twist<f64>, _, 6>(
        &<TwistBlockJac<f64> as Jac<f64, Twist<f64>>>::identity(),
    );
    let mut rng = laws::Rng(0x6368_6172_7473_0004);
    for _ in 0..10_000 {
        let x = G::<f64>::exp(&laws::tangent::<f64, G<f64>, 6>(&rng.shaped::<6>()));
        for j in [
            x.chart_transition::<Screw<f64>, Decoupled<f64>>(),
            x.chart_transition::<Decoupled<f64>, Screw<f64>>(),
            x.chart_transition::<WorldTranslation<f64>, WorldTranslation<f64>>(),
        ] {
            assert_eq!(laws::dense_bits::<Twist<f64>, _, 6>(&j), id);
        }
    }
}

/// `TwistBlockJac` is the pair it wraps, to the bit, for every operation (`0060` decision 6), so
/// the pair's `productjac_sandwich_matches_reference` twin covers its `sandwich`.
#[test]
fn twist_block_jac_is_the_pair_bit_for_bit() {
    type P = ProductJac<Mat3<f64>, Mat3<f64>>;
    type T2 = (SO3Tangent<f64>, SO3Tangent<f64>);
    let mut rng = laws::Rng(0x6368_6172_7473_0005);
    let mat = |rng: &mut laws::Rng| {
        Matrix::from_cols([
            Vector(rng.arr::<3>()),
            Vector(rng.arr::<3>()),
            Vector(rng.arr::<3>()),
        ])
    };
    let pair_bits = |p: &P| laws::dense_bits::<T2, P, 6>(p);
    let new_bits = |j: &TwistBlockJac<f64>| laws::dense_bits::<Twist<f64>, _, 6>(j);
    for _ in 0..10_000 {
        let (a, b, c, d) = (mat(&mut rng), mat(&mut rng), mat(&mut rng), mat(&mut rng));
        let (j, k) = (
            TwistBlockJac::from_blocks(a, b),
            TwistBlockJac::from_blocks(c, d),
        );
        let (p, q): (P, P) = (ProductJac(a, b), ProductJac(c, d));
        assert_eq!(new_bits(&j), pair_bits(&p));
        assert_eq!(
            new_bits(&j.mul(&k)),
            pair_bits(&Jac::<f64, T2>::mul(&p, &q))
        );
        assert_eq!(
            new_bits(&j.inverse()),
            pair_bits(&Jac::<f64, T2>::inverse(&p))
        );
        assert_eq!(new_bits(&j.neg()), pair_bits(&Jac::<f64, T2>::neg(&p)));
        let v = [rng.arr::<3>(), rng.arr::<3>()];
        let tw = Twist {
            phi: Vector(v[0]),
            rho: [Vector(v[1])],
        };
        let pt = (
            SO3Tangent { phi: Vector(v[0]) },
            SO3Tangent { phi: Vector(v[1]) },
        );
        let dense = |t: &Twist<f64>| {
            let mut o = [0.0; 6];
            t.write_dense(&mut o);
            o.map(f64::to_bits)
        };
        let pdense = |t: &T2| {
            let mut o = [0.0; 6];
            t.write_dense(&mut o);
            o.map(f64::to_bits)
        };
        assert_eq!(
            dense(&j.apply(&tw)),
            pdense(&Jac::<f64, T2>::apply(&p, &pt))
        );
        assert_eq!(
            dense(&j.apply_transpose(&tw)),
            pdense(&Jac::<f64, T2>::apply_transpose(&p, &pt))
        );
        let cov = Matrix::<f64, 6, 6>::from_cols(core::array::from_fn(|_| Vector(rng.arr::<6>())));
        let (s1, s2) = (
            j.sandwich::<6>(&cov),
            Jac::<f64, T2>::sandwich::<6>(&p, &cov),
        );
        for i in 0..36 {
            assert_eq!(
                s1.get(i % 6, i / 6).to_bits(),
                s2.get(i % 6, i / 6).to_bits()
            );
        }
        assert_eq!(
            j.rotation_block().get(1, 2).to_bits(),
            a.get(1, 2).to_bits()
        );
        assert_eq!(
            j.translation_block().get(2, 0).to_bits(),
            b.get(2, 0).to_bits()
        );
    }
}

/// locus-tag's `Pose::retract` is `Decoupled` (`0012`): `t + R ρ` in the body frame, against a
/// hand-built `(R Exp φ, t + R ρ)` at one point, and `WorldTranslation` is `t + ρ`.
#[test]
fn decoupled_moves_the_translation_in_the_body_frame() {
    let r = G::<f64>::exp(&Twist {
        phi: Vector([0.0, 0.0, 1.2]),
        rho: [Vector([0.0; 3])],
    });
    let x = SE3::from_rt(r.rotation(), Vector([1.0, 2.0, 3.0]));
    let d = Twist {
        phi: Vector([0.0; 3]),
        rho: [Vector([1e-8, 0.0, 0.0])],
    };
    let body = Decoupled::at(&x).retract(&d).translation();
    let world = WorldTranslation::at(&x).retract(&d).translation();
    // `R e_x 1e-8 = (cos 1.2, sin 1.2, 0) 1e-8`.
    let want = [1.0 + 1.2_f64.cos() * 1e-8, 2.0 + 1.2_f64.sin() * 1e-8, 3.0];
    for (got, want) in body.0.iter().zip(want) {
        assert!((got - want).abs() <= 4.0 * f64::EPSILON * 3.0);
    }
    assert_eq!(
        world.0.map(f64::to_bits),
        [1.0 + 1e-8, 2.0, 3.0].map(f64::to_bits)
    );
    // CH.5(c)'s figure: the two differ at first order by `2 sin(0.6) ‖ρ‖`.
    let gap = laws::norm(&core::array::from_fn::<f64, 3, _>(|i| {
        body.0[i] - world.0[i]
    }));
    assert!((gap / 1e-8 - 2.0 * 0.6_f64.sin()).abs() < 1e-7);
}
