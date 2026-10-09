//! The group charts (`RightChart`, `LeftChart`) under the chart laws of [`crate::laws`]
//! (`docs/PHASE5.md` §1.1, `0060`), on every group.
//!
//! Samples are `laws::sample`'s, so a rotation part has `θ < √3 < π`: every case is inside `U₀`,
//! where `local(retract(δ)) = δ` holds (CH.2). The bounds are recorded as described at
//! [`ChartBounds`], from `measure_charts`.

use crate::laws::{self, ChartBounds, ChartLegs, CHART_LEGS};
use crate::{
    Chart, LeftChart, LieGroup, Lifted, Manifold, Product, RightChart, Rn, SEn3, Tangent,
    WithChart, SO3,
};

/// The chart laws for the group `$G<S>` in the chart `$C<S>`, at `f64` and `f32`; the `Dual` laws
/// at `f64` (they differentiate at one scalar and read the closed forms at the same one). The
/// transition and change-of-chart laws take `$C` to `$C2`, whose first-order transition at `x` is
/// `$phi0(&x)`.
macro_rules! chart_laws_for {
    ($p:ident, $G:ident, $C:ident, $C2:ident, $phi0:expr, $D:literal, $b64:ident, $b32:ident) => {
        mod $p {
            use super::*;
            chart_laws_for!(@case as_f64, f64, $b64.legs, $G, $C, $D);
            chart_laws_for!(@case as_f32, f32, $b32, $G, $C, $D);

            type Dl = helicoid_linalg::Dual<f64, $D>;

            proptest::proptest! {
                #[test]
                fn chart_jacobians_match_dual(
                    a in laws::sample::<$D>(), b in laws::sample::<$D>(), c in laws::sample::<$D>(),
                ) {
                    let v = laws::chart_jacobians_match_dual::<f64, $G<Dl>, $C<Dl>, $D>(&a, &b, &c);
                    laws::within(v, $b64.dual)?;
                }
                #[test]
                fn transition_at_zero(a in laws::sample::<$D>()) {
                    let v = laws::transition_matches::<f64, $G<Dl>, $C<Dl>, $C2<Dl>, _, $D>(&a, $phi0);
                    laws::within(v, $b64.transition)?;
                }
                #[test]
                fn change_of_chart(a in laws::sample::<$D>(), b in laws::sample::<$D>()) {
                    let v = laws::change_of_chart::<f64, $G<Dl>, $C<Dl>, $C2<Dl>, $D>(&a, &b);
                    laws::within(v, $b64.change)?;
                }
            }

            /// The worst error of each chart law over 10^6 seeded cases (10^5 for the `Dual` ones),
            /// the figures the bounds are recorded from. Its own stream (`laws::Rng`).
            #[test]
            #[ignore = "measurement: prints the figures the bounds are recorded from"]
            #[allow(clippy::print_stdout)]
            fn measure_charts() {
                let mut rng = laws::Rng(0x6368_6172_7473_0001);
                let mut legs = [[0.0_f64; CHART_LEGS.len()]; 2];
                for _ in 0..1_000_000 {
                    let (a, b, c) = (rng.shaped::<$D>(), rng.shaped::<$D>(), rng.shaped::<$D>());
                    let l64 = as_f64::legs(&a, &b, &c);
                    let l32 = as_f32::legs(&a, &b, &c);
                    for i in 0..CHART_LEGS.len() {
                        legs[0][i] = laws::worst(legs[0][i], l64[i]);
                        legs[1][i] = laws::worst(legs[1][i], l32[i]);
                    }
                }
                let mut dual = [0.0_f64; 3];
                for _ in 0..100_000 {
                    let (a, b, c) = (rng.shaped::<$D>(), rng.shaped::<$D>(), rng.shaped::<$D>());
                    let v = [
                        laws::chart_jacobians_match_dual::<f64, $G<Dl>, $C<Dl>, $D>(&a, &b, &c),
                        laws::transition_matches::<f64, $G<Dl>, $C<Dl>, $C2<Dl>, _, $D>(&a, $phi0),
                        laws::change_of_chart::<f64, $G<Dl>, $C<Dl>, $C2<Dl>, $D>(&a, &b),
                    ];
                    for (w, x) in dual.iter_mut().zip(v) {
                        *w = laws::worst(*w, x);
                    }
                }
                std::println!(
                    "{} legs f64 {:.3?} f32 {:.3?} dual {:.3} transition {:.3} change {:.3}",
                    module_path!(), legs[0], legs[1], dual[0], dual[1], dual[2],
                );
            }
        }
    };
    (@case $m:ident, $S:ty, $B:expr, $G:ident, $C:ident, $D:literal) => {
        mod $m {
            use super::*;

            type S = $S;

            pub(super) fn legs(a: &[f64; $D], b: &[f64; $D], c: &[f64; $D]) -> [f64; CHART_LEGS.len()] {
                let t = |v: &[f64; $D]| laws::tangent::<S, $G<S>, $D>(v);
                let (x, y) = (<$G<S>>::exp(&t(a)), <$G<S>>::exp(&t(c)));
                laws::chart_legs::<S, $G<S>, $C<S>, $D>(&x, &t(b), &y)
            }

            proptest::proptest! {
                #[test]
                fn chart_legs(
                    a in laws::sample::<$D>(), b in laws::sample::<$D>(), c in laws::sample::<$D>(),
                ) {
                    for ((name, v), bound) in CHART_LEGS.iter().zip(legs(&a, &b, &c)).zip($B) {
                        laws::within_leg(name, v, bound)?;
                    }
                }
                #[test]
                fn local_with_jacobian_is_the_two_calls(a in laws::sample::<$D>(), c in laws::sample::<$D>()) {
                    let t = |v: &[f64; $D]| laws::tangent::<S, $G<S>, $D>(v);
                    let (x, y) = (<$G<S>>::exp(&t(&a)), <$G<S>>::exp(&t(&c)));
                    proptest::prop_assert!(laws::local_with_jacobian_is_the_two_calls::<S, $G<S>, $C<S>, $D>(&x, &y));
                }
            }
        }
    };
}

pub(crate) use chart_laws_for;

/// The right chart's first-order transition to the left one, `Ad_X` (CH.6).
fn ad<S: helicoid_linalg::Real, G: LieGroup<S>>(x: &G) -> G::Jac {
    x.adjoint()
}

/// The left chart's to the right one, `Ad_{X⁻¹}`.
fn ad_inv<S: helicoid_linalg::Real, G: LieGroup<S>>(x: &G) -> G::Jac {
    x.inverse().adjoint()
}

// `product_left`, worst of `measure_charts`: legs `f64` [1.116, 6.247, 8.072, 1.118], `f32` [1.068, 6.134, 8.094, 1.118]; `Dual` 3.422, transition 3.304,
// change of chart 5.174.
const PRODUCT_LEFT_F64: ChartBounds = ChartBounds {
    legs: [3.0, 13.0, 17.0, 3.0],
    dual: 7.0,
    transition: 7.0,
    change: 11.0,
};
const PRODUCT_LEFT_F32: ChartLegs = [3.0, 13.0, 17.0, 3.0];

// `product_right`, worst of `measure_charts`: legs `f64` [1.116, 6.06, 8.013, 1.118], `f32` [1.068, 6.102, 7.802, 1.118]; `Dual` 3.489, transition 2.739,
// change of chart 4.779.
const PRODUCT_RIGHT_F64: ChartBounds = ChartBounds {
    legs: [3.0, 13.0, 17.0, 3.0],
    dual: 7.0,
    transition: 6.0,
    change: 10.0,
};
const PRODUCT_RIGHT_F32: ChartLegs = [3.0, 13.0, 16.0, 3.0];

// `rn_left`, worst of `measure_charts`: legs `f64` [0.0, 1.5, 1.5, 0.0], `f32` [0.0, 1.732, 1.5, 0.0]; `Dual` 0.0, transition 0.0,
// change of chart 0.0.
const RN_LEFT_F64: ChartBounds = ChartBounds {
    legs: [0.0, 3.0, 3.0, 0.0],
    dual: 0.0,
    transition: 0.0,
    change: 0.0,
};
const RN_LEFT_F32: ChartLegs = [0.0, 4.0, 3.0, 0.0];

// `rn_right`, worst of `measure_charts`: legs `f64` [0.0, 1.5, 1.5, 0.0], `f32` [0.0, 1.732, 1.5, 0.0]; `Dual` 0.0, transition 0.0,
// change of chart 0.0.
const RN_RIGHT_F64: ChartBounds = ChartBounds {
    legs: [0.0, 3.0, 3.0, 0.0],
    dual: 0.0,
    transition: 0.0,
    change: 0.0,
};
const RN_RIGHT_F32: ChartLegs = [0.0, 4.0, 3.0, 0.0];

// `se23_left`, worst of `measure_charts`: legs `f64` [1.092, 7.412, 12.629, 1.92], `f32` [1.052, 8.154, 10.577, 1.804]; `Dual` 7.239, transition 2.561,
// change of chart 6.793.
const SE23_LEFT_F64: ChartBounds = ChartBounds {
    legs: [3.0, 15.0, 26.0, 4.0],
    dual: 15.0,
    transition: 6.0,
    change: 14.0,
};
const SE23_LEFT_F32: ChartLegs = [3.0, 17.0, 22.0, 4.0];

// `se23_right`, worst of `measure_charts`: legs `f64` [1.092, 6.188, 8.262, 1.118], `f32` [1.052, 6.16, 8.009, 1.118]; `Dual` 5.896, transition 2.425,
// change of chart 7.559.
const SE23_RIGHT_F64: ChartBounds = ChartBounds {
    legs: [3.0, 13.0, 17.0, 3.0],
    dual: 12.0,
    transition: 5.0,
    change: 16.0,
};
const SE23_RIGHT_F32: ChartLegs = [3.0, 13.0, 17.0, 3.0];

// `se3_left`, worst of `measure_charts`: legs `f64` [1.116, 7.282, 10.193, 1.871], `f32` [1.068, 8.762, 9.949, 1.734]; `Dual` 8.484, transition 3.592,
// change of chart 6.845.
const SE3_LEFT_F64: ChartBounds = ChartBounds {
    legs: [3.0, 15.0, 21.0, 4.0],
    dual: 17.0,
    transition: 8.0,
    change: 14.0,
};
const SE3_LEFT_F32: ChartLegs = [3.0, 18.0, 20.0, 4.0];

// `se3_right`, worst of `measure_charts`: legs `f64` [1.116, 6.212, 9.106, 1.118], `f32` [1.068, 6.183, 7.984, 1.118]; `Dual` 4.849, transition 3.082,
// change of chart 7.477.
const SE3_RIGHT_F64: ChartBounds = ChartBounds {
    legs: [3.0, 13.0, 19.0, 3.0],
    dual: 10.0,
    transition: 7.0,
    change: 15.0,
};
const SE3_RIGHT_F32: ChartLegs = [3.0, 13.0, 16.0, 3.0];

// `sen3_n3_left`, worst of `measure_charts`: legs `f64` [1.033, 7.911, 12.455, 1.961], `f32` [1.048, 8.756, 10.912, 1.939]; `Dual` 7.914, transition 2.606,
// change of chart 7.699.
const SEN3_N3_LEFT_F64: ChartBounds = ChartBounds {
    legs: [3.0, 16.0, 25.0, 4.0],
    dual: 16.0,
    transition: 6.0,
    change: 16.0,
};
const SEN3_N3_LEFT_F32: ChartLegs = [3.0, 18.0, 22.0, 4.0];

// `sen3_n3_right`, worst of `measure_charts`: legs `f64` [1.033, 6.039, 9.24, 1.118], `f32` [1.048, 6.294, 7.939, 1.118]; `Dual` 5.826, transition 2.346,
// change of chart 8.021.
const SEN3_N3_RIGHT_F64: ChartBounds = ChartBounds {
    legs: [3.0, 13.0, 19.0, 3.0],
    dual: 12.0,
    transition: 5.0,
    change: 17.0,
};
const SEN3_N3_RIGHT_F32: ChartLegs = [3.0, 13.0, 16.0, 3.0];

// `so3_left`, worst of `measure_charts`: legs `f64` [1.08, 6.032, 8.162, 1.118], `f32` [1.079, 6.088, 8.644, 1.118]; `Dual` 5.257, transition 3.391,
// change of chart 7.293.
const SO3_LEFT_F64: ChartBounds = ChartBounds {
    legs: [3.0, 13.0, 17.0, 3.0],
    dual: 11.0,
    transition: 7.0,
    change: 15.0,
};
const SO3_LEFT_F32: ChartLegs = [3.0, 13.0, 18.0, 3.0];

// `so3_right`, worst of `measure_charts`: legs `f64` [1.08, 6.63, 8.27, 1.118], `f32` [1.079, 7.022, 9.074, 1.118]; `Dual` 4.954, transition 3.464,
// change of chart 6.565.
const SO3_RIGHT_F64: ChartBounds = ChartBounds {
    legs: [3.0, 14.0, 17.0, 3.0],
    dual: 10.0,
    transition: 7.0,
    change: 14.0,
};
const SO3_RIGHT_F32: ChartLegs = [3.0, 15.0, 19.0, 3.0];

type So3<S> = SO3<S>;
type So3R<S> = RightChart<SO3<S>>;
type So3L<S> = LeftChart<SO3<S>>;
chart_laws_for!(
    so3_right,
    So3,
    So3R,
    So3L,
    ad,
    3,
    SO3_RIGHT_F64,
    SO3_RIGHT_F32
);
chart_laws_for!(
    so3_left,
    So3,
    So3L,
    So3R,
    ad_inv,
    3,
    SO3_LEFT_F64,
    SO3_LEFT_F32
);

type Se3<S> = SEn3<S, 1>;
type Se3R<S> = RightChart<Se3<S>>;
type Se3L<S> = LeftChart<Se3<S>>;
chart_laws_for!(
    se3_right,
    Se3,
    Se3R,
    Se3L,
    ad,
    6,
    SE3_RIGHT_F64,
    SE3_RIGHT_F32
);
chart_laws_for!(
    se3_left,
    Se3,
    Se3L,
    Se3R,
    ad_inv,
    6,
    SE3_LEFT_F64,
    SE3_LEFT_F32
);

type Se23<S> = SEn3<S, 2>;
type Se23R<S> = RightChart<Se23<S>>;
type Se23L<S> = LeftChart<Se23<S>>;
chart_laws_for!(
    se23_right,
    Se23,
    Se23R,
    Se23L,
    ad,
    9,
    SE23_RIGHT_F64,
    SE23_RIGHT_F32
);
chart_laws_for!(
    se23_left,
    Se23,
    Se23L,
    Se23R,
    ad_inv,
    9,
    SE23_LEFT_F64,
    SE23_LEFT_F32
);

type Sen33<S> = SEn3<S, 3>;
type Sen33R<S> = RightChart<Sen33<S>>;
type Sen33L<S> = LeftChart<Sen33<S>>;
chart_laws_for!(
    sen3_n3_right,
    Sen33,
    Sen33R,
    Sen33L,
    ad,
    12,
    SEN3_N3_RIGHT_F64,
    SEN3_N3_RIGHT_F32
);
chart_laws_for!(
    sen3_n3_left,
    Sen33,
    Sen33L,
    Sen33R,
    ad_inv,
    12,
    SEN3_N3_LEFT_F64,
    SEN3_N3_LEFT_F32
);

type R3<S> = Rn<S, 3>;
type R3R<S> = RightChart<R3<S>>;
type R3L<S> = LeftChart<R3<S>>;
chart_laws_for!(rn_right, R3, R3R, R3L, ad, 3, RN_RIGHT_F64, RN_RIGHT_F32);
chart_laws_for!(rn_left, R3, R3L, R3R, ad_inv, 3, RN_LEFT_F64, RN_LEFT_F32);

type Pr<S> = Product<SO3<S>, Rn<S, 3>>;
type PrR<S> = RightChart<Pr<S>>;
type PrL<S> = LeftChart<Pr<S>>;
chart_laws_for!(
    product_right,
    Pr,
    PrR,
    PrL,
    ad,
    6,
    PRODUCT_RIGHT_F64,
    PRODUCT_RIGHT_F32
);
chart_laws_for!(
    product_left,
    Pr,
    PrL,
    PrR,
    ad_inv,
    6,
    PRODUCT_LEFT_F64,
    PRODUCT_LEFT_F32
);

/// A group's default chart is its right chart (`docs/PHASE5.md` §1.1): this compiles only if
/// `Manifold::Chart` is `RightChart<Self>` for each group.
#[test]
fn every_group_defaults_to_its_right_chart() {
    fn right<S: helicoid_linalg::Real, G: LieGroup<S> + Manifold<S, Chart = RightChart<G>>>() {}
    right::<f64, So3<f64>>();
    right::<f64, Se3<f64>>();
    right::<f64, Se23<f64>>();
    right::<f64, Sen33<f64>>();
    right::<f64, R3<f64>>();
    right::<f32, Pr<f32>>();
    assert_eq!(
        <Se3<f64> as Manifold<f64>>::DOF,
        <Se3<f64> as LieGroup<f64>>::DOF
    );
}

/// `WithChart<M, C>` linearizes exactly as `C` does: the same bits for `retract`, `local` and both
/// Jacobians (`0060` decision 5).
#[test]
fn with_chart_is_bit_identical_to_its_chart() {
    type W = WithChart<Se3<f64>, Se3L<f64>>;
    let mut rng = laws::Rng(0x6368_6172_7473_0002);
    let bits = |t: &<Se3<f64> as LieGroup<f64>>::Tangent| {
        let mut buf = [0.0; 6];
        t.write_dense(&mut buf);
        buf.map(f64::to_bits)
    };
    for _ in 0..10_000 {
        let (a, b, c) = (rng.shaped::<6>(), rng.shaped::<6>(), rng.shaped::<6>());
        let t = |v: &[f64; 6]| laws::tangent::<f64, Se3<f64>, 6>(v);
        let (x, d, y) = (Se3::<f64>::exp(&t(&a)), t(&b), Se3::<f64>::exp(&t(&c)));
        let lifted: Lifted<Se3L<f64>> = <W as Manifold<f64>>::Chart::at(&W::new(x));
        let plain = Se3L::<f64>::at(&x);
        assert_eq!(
            bits(&lifted.retract(&d).0.log()),
            bits(&plain.retract(&d).log())
        );
        assert_eq!(bits(&lifted.local(&W::new(y))), bits(&plain.local(&y)));
        assert_eq!(
            laws::dense_bits::<_, _, 6>(&lifted.retract_jacobian(&d)),
            laws::dense_bits::<_, _, 6>(&plain.retract_jacobian(&d)),
        );
        assert_eq!(
            laws::dense_bits::<_, _, 6>(&lifted.local_jacobian(&W::new(y))),
            laws::dense_bits::<_, _, 6>(&plain.local_jacobian(&y)),
        );
        assert_eq!(bits(&lifted.base().0.log()), bits(&x.log()));
    }
}
