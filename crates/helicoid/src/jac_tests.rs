//! `Jac::add` and `Jac::scale` on every impl, against the entrywise dense operation (`0112`).
//!
//! Each impl is built from samples uniform on `[-1, 1)`, one per stored entry, so every stored
//! entry is nonzero and a `+0` in the dense image is a structural zero
//! ([`laws::jac_add_and_scale_match_dense`]). The factor `k` includes `-0`, `-1`, `±∞` and NaN.

use crate::dualmat_tests::Lanes;
use crate::heis_tests::HJac;
use crate::laws;
use crate::{
    Jac, ProductJac, RnJac, RnTangent, SEn3Jac, SEn3Tangent, SO3Tangent, Sim3Jac, Sim3Tangent,
    Tangent, Twist, TwistBlockJac,
};
use core::array;
use helicoid_linalg::{Mat2, Mat3, Matrix, Vector};
use proptest::prelude::*;

/// A `Jac` with `LEN` stored entries, built from as many samples; `Dual` sample `i` is variable
/// `i mod 3`, so every lane is exercised.
trait Stored<S: Lanes>: Sized {
    const LEN: usize;
    fn of(v: &[f64]) -> Self;
}

fn s<S: Lanes>(v: &[f64], i: usize) -> S {
    S::sample(v[i], i % 3)
}

/// `3 x 3`, column-major from sample `off`.
fn mat3<S: Lanes>(v: &[f64], off: usize) -> Mat3<S> {
    Matrix::from_cols(array::from_fn(|c| {
        Vector(array::from_fn(|r| s(v, off + 3 * c + r)))
    }))
}

impl<S: Lanes> Stored<S> for Mat3<S> {
    const LEN: usize = 9;
    fn of(v: &[f64]) -> Self {
        mat3(v, 0)
    }
}

impl<S: Lanes> Stored<S> for Mat2<S> {
    const LEN: usize = 4;
    fn of(v: &[f64]) -> Self {
        Matrix::from_cols(array::from_fn(|c| {
            Vector(array::from_fn(|r| s(v, 2 * c + r)))
        }))
    }
}

impl<S: Lanes, const N: usize> Stored<S> for RnJac<S, N> {
    const LEN: usize = 1;
    fn of(v: &[f64]) -> Self {
        RnJac::scalar(s(v, 0))
    }
}

impl<S: Lanes, const N: usize> Stored<S> for SEn3Jac<S, N> {
    const LEN: usize = 9 * (N + 1);
    fn of(v: &[f64]) -> Self {
        SEn3Jac {
            diag: mat3(v, 0),
            col: array::from_fn(|i| mat3(v, 9 * (i + 1))),
        }
    }
}

impl<S: Lanes> Stored<S> for Sim3Jac<S> {
    const LEN: usize = 31;
    fn of(v: &[f64]) -> Self {
        Sim3Jac {
            a: mat3(v, 0),
            b: mat3(v, 9),
            c: mat3(v, 18),
            f: Vector(array::from_fn(|i| s(v, 27 + i))),
            e: s(v, 30),
        }
    }
}

impl<S: Lanes, JA: Stored<S>, JB: Stored<S>> Stored<S> for ProductJac<JA, JB> {
    const LEN: usize = JA::LEN + JB::LEN;
    fn of(v: &[f64]) -> Self {
        ProductJac(JA::of(v), JB::of(&v[JA::LEN..]))
    }
}

impl<S: Lanes> Stored<S> for TwistBlockJac<S> {
    const LEN: usize = 18;
    fn of(v: &[f64]) -> Self {
        TwistBlockJac::from_blocks(mat3(v, 0), mat3(v, 9))
    }
}

impl<S: Lanes> Stored<S> for HJac<S> {
    const LEN: usize = 9;
    fn of(v: &[f64]) -> Self {
        HJac(mat3(v, 0))
    }
}

/// `a` from the first `LEN` samples, `b` from the next `LEN`.
fn check<S, T, J, const D: usize>(v: &[f64], k: f64) -> bool
where
    S: Lanes,
    T: Tangent<S>,
    J: Jac<S, T> + Stored<S>,
{
    let (a, b) = (J::of(v), J::of(&v[J::LEN..]));
    laws::jac_add_and_scale_match_dense::<S, T, J, D>(&a, &b, S::sample(k, 0))
}

macro_rules! impls {
    ($($m:ident: $J:ty, $T:ty, $D:literal;)*) => {$(
        mod $m {
            use super::*;
            impls!(@case as_f64, f64, $J, $T, $D);
            impls!(@case as_f32, f32, $J, $T, $D);
            impls!(@case as_dual, helicoid_linalg::Dual<f64, 3>, $J, $T, $D);
        }
    )*};
    (@case $c:ident, $S:ty, $J:ty, $T:ty, $D:literal) => {
        mod $c {
            use super::*;
            type S = $S;
            proptest! {
                #[test]
                fn jac_add_and_scale_match_dense(
                    v in proptest::collection::vec(-1.0_f64..1.0, 2 * <$J as Stored<S>>::LEN),
                    k in laws::factor(),
                ) {
                    prop_assert!(check::<S, $T, $J, $D>(&v, k));
                }
            }
        }
    };
}

impls! {
    mat3: Mat3<S>, SO3Tangent<S>, 3;
    mat2: Mat2<S>, RnTangent<S, 2>, 2;
    rnjac: RnJac<S, 3>, RnTangent<S, 3>, 3;
    sen3jac_1: SEn3Jac<S, 1>, SEn3Tangent<S, 1>, 6;
    sen3jac_2: SEn3Jac<S, 2>, SEn3Tangent<S, 2>, 9;
    sim3jac: Sim3Jac<S>, Sim3Tangent<S>, 7;
    productjac: ProductJac<SEn3Jac<S, 1>, RnJac<S, 2>>, (SEn3Tangent<S, 1>, RnTangent<S, 2>), 8;
    twistblockjac: TwistBlockJac<S>, Twist<S>, 6;
    hjac: HJac<S>, RnTangent<S, 3>, 3;
}
