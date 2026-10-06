//! `Rn` under the generic laws of [`crate::laws`], and its hand cases.
//!
//! The bounds are recorded as described at [`Bounds`]. The worst errors, in `u`: `f64` and
//! `Dual<f64, 3>` (same samples, same values) 3.42 for `group_axioms`, 1.5 for `adjoint_identity`
//! and `plus_minus`, 1 for `jac_dense_order` (`k (1/k)`) and 0 for the rest; `f32` 2.72, 1.47, 1.5,
//! 1, 1.61 for `tangent_dense_order`, 1.66 for `sandwich_matches_dense` and 0 for the rest. The
//! zeros pin the operation order of `f64`: the `f32` figures are one rounding against an `f64`
//! reference.

use crate::laws::{laws_for, Bounds, Sample};
use crate::{Jac, Left, LieGroup, Right, Rn, RnJac, RnTangent, Side, Tangent};
use core::array;
use helicoid_linalg::{Blend, StridedMut, Vector};
use std::vec::Vec;

type G<S> = Rn<S, 3>;

const F64: Bounds = Bounds {
    axioms: 7.0,
    exp_log: 0.0,
    adjoint: 3.0,
    jl_ad_jr: 0.0,
    plus_minus: 3.0,
    rows: 0.0,
    ad: 0.0,
    sides: 0.0,
    tangent_order: 0.0,
    jac_order: 2.0,
    sandwich: 0.0,
    // `PHASE3.md` §8's second check, exactly `0`: every row is `±I` or `k I` and the operations are additions, so the dual derivative is the closed form.
    dual_rows: 0.0,
    // `PHASE4.md` §1 and §3's seven legs, in `GEODESIC_LEGS`'s order, each twice the worst of 10^6 draws of `laws::Rng::shaped` -- `laws::sample`'s own distribution, which is what the proptest draws -- rounded up: symmetry 3.288 / 3.269, velocity 3.538 / 3.288, left and right 4.031 / 3.947 (`f64` and `Dual` / `f32`). `t=0`, `t=1` and `twin` are **exactly 0**: `Rn`'s inverse is exact, so both endpoints come back unmoved and the provided body is `reference::geodesic` to the bit.
    geodesic: [0.0, 0.0, 7.0, 8.0, 9.0, 9.0, 0.0],
};
const F32: Bounds = Bounds {
    axioms: 6.0,
    tangent_order: 4.0,
    sandwich: 4.0,
    dual_rows: 0.0,
    ..F64
};

// `k I` with `0.5 <= |k| < 1.5`: invertible, of either sign.
fn rn_jac<S: Sample>(v: &[f64; 3]) -> RnJac<S, 3> {
    RnJac::scalar(S::sample(v[0] + 0.5_f64.copysign(v[0]), 0))
}

laws_for!(rn, G, rn_jac, 3, F64, F32);

fn bits(x: &[f64]) -> Vec<u64> {
    x.iter().map(|v| v.to_bits()).collect()
}

fn dense(j: &RnJac<f64, 3>) -> [[u64; 3]; 3] {
    let mut buf = [[f64::NAN; 3]; 3];
    j.write_dense(&mut StridedMut::col_major(buf.as_flattened_mut(), 3, 3));
    buf.map(|c| c.map(f64::to_bits))
}

/// `s I`, with `+0` off the diagonal.
fn scaled_identity(s: f64) -> [[u64; 3]; 3] {
    array::from_fn(|c| array::from_fn(|r| if r == c { s } else { 0.0 }.to_bits()))
}

#[test]
fn the_group_law_is_addition() {
    let (a, b) = (Rn(Vector([1.0, -2.0, 0.5])), Rn(Vector([4.0, 2.0, -0.5])));
    assert_eq!(bits(&(a * b).0 .0), bits(&[5.0, 0.0, 0.0]));
    assert_eq!(bits(&G::<f64>::identity().0 .0), bits(&[0.0; 3]));
    assert_eq!(bits(&a.inverse().0 .0), bits(&[-1.0, 2.0, -0.5]));
    assert_eq!(bits(&a.rplus(&a.log()).0 .0), bits(&[2.0, -4.0, 1.0]));
    assert_eq!(bits(&a.rminus(&b).rho.0), bits(&[-3.0, -4.0, 1.0]));
    // Not sign-symmetric, so the left forms are pinned by value too.
    assert_eq!(bits(&a.lplus(&b.log()).0 .0), bits(&[5.0, 0.0, 0.0]));
    assert_eq!(bits(&a.lminus(&b).rho.0), bits(&[-3.0, -4.0, 1.0]));
    assert_eq!(
        (G::<f64>::DOF, <RnTangent<f64, 3> as Tangent<f64>>::DOF),
        (3, 3)
    );
}

#[test]
fn dense_order_is_index_order() {
    let t = RnTangent {
        rho: Vector([1.0, 2.0, 3.0]),
    };
    let mut out = [0.0; 3];
    t.write_dense(&mut out);
    assert_eq!(bits(&out), bits(&[1.0, 2.0, 3.0]));
    assert_eq!(
        bits(&RnTangent::<f64, 3>::read_dense(&out).rho.0),
        bits(&out)
    );
    // Row-major and padded views: the block is written and its padding is not.
    let j = RnJac::<f64, 3>::scalar(2.0);
    let mut padded = [f64::NAN; 15];
    j.write_dense(&mut StridedMut::with_strides(&mut padded, 3, 3, 1, 5));
    for (i, v) in padded.iter().enumerate() {
        let (r, c) = (i % 5, i / 5);
        if r >= 3 {
            assert!(v.is_nan(), "padding {i} was written");
        } else {
            assert_eq!(
                v.to_bits(),
                if r == c { 2.0_f64 } else { 0.0 }.to_bits(),
                "entry {i}"
            );
        }
    }
}

#[test]
fn jacobian_rows_are_the_abelian_case() {
    let (x, y, tau) = (
        Rn(Vector([1.0, 2.0, 3.0])),
        Rn(Vector([-4.0, 5.0, 6.0])),
        RnTangent {
            rho: Vector([0.5; 3]),
        },
    );
    let (i, m) = (scaled_identity(1.0), scaled_identity(-1.0));
    let rows = |p: (RnJac<f64, 3>, RnJac<f64, 3>)| (dense(&p.0), dense(&p.1));
    // NUMERICS.md §2.3 with `Ad = J = J⁻¹ = I`.
    assert_eq!(rows(x.rplus_jacobians(&tau)), (i, i));
    assert_eq!(rows(x.lplus_jacobians(&tau)), (i, i));
    assert_eq!(rows(x.rminus_jacobians(&y)), (i, m));
    assert_eq!(rows(x.lminus_jacobians(&y)), (i, m));
    assert_eq!(rows(x.compose_jacobians::<Right>(&y)), (i, i));
    assert_eq!(rows(x.compose_jacobians::<Left>(&y)), (i, i));
    assert_eq!(dense(&x.inverse_jacobian::<Right>()), m);
    assert_eq!(dense(&x.inverse_jacobian::<Left>()), m);
    assert_eq!(dense(&x.adjoint()), i);
    assert_eq!(dense(&G::jr(&tau)), i);
    assert_eq!(dense(&G::jl_inv(&tau)), i);
    assert_eq!(dense(&G::ad(&tau)), scaled_identity(0.0));
    let (jr, jl) = (
        Right::plus_jacobians(&x, &tau),
        Left::plus_jacobians(&x, &tau),
    );
    assert_eq!((dense(&jr.1), dense(&jl.1)), (i, i));
}

#[test]
fn blend_selects_each_value_type() {
    let (a, b) = (Rn(Vector([1.0; 3])), Rn(Vector([2.0; 3])));
    let (ta, tb) = (a.log(), b.log());
    let (ja, jb) = (RnJac::<f64, 3>::scalar(1.0), RnJac::scalar(-1.0));
    for (m, x, t, j) in [(true, a, ta, ja), (false, b, tb, jb)] {
        assert_eq!(bits(&Rn::blend(m, a, b).0 .0), bits(&x.0 .0));
        assert_eq!(bits(&RnTangent::blend(m, ta, tb).rho.0), bits(&t.rho.0));
        assert_eq!(dense(&RnJac::blend(m, ja, jb)), dense(&j));
    }
}

/// A domain is a `debug_assert!` (D11): the debug profile that `just test` runs panics ...
#[cfg(debug_assertions)]
mod out_of_domain {
    use super::*;

    #[test]
    #[should_panic(expected = "Tangent::write_dense")]
    fn write_dense_of_the_wrong_length() {
        RnTangent::<f64, 3>::zero().write_dense(&mut [0.0; 4]);
    }

    #[test]
    #[should_panic(expected = "Tangent::read_dense")]
    fn read_dense_of_the_wrong_length() {
        let _ = RnTangent::<f64, 3>::read_dense(&[0.0; 2]);
    }

    // A view larger than `DOF x DOF` is no strided-access panic, so only the assertion sees it.
    #[test]
    #[should_panic(expected = "Jac::write_dense")]
    fn write_dense_into_a_larger_view() {
        let mut buf = [0.0; 16];
        RnJac::<f64, 3>::identity().write_dense(&mut StridedMut::col_major(&mut buf, 4, 4));
    }

    #[test]
    #[should_panic(expected = "singular")]
    fn ad_has_no_inverse() {
        let _ = G::<f64>::ad(&RnTangent::zero()).inverse();
    }
}

/// ... and a release build (`cargo nextest run --release`) never panics: `write_dense` writes
/// `min(len, DOF)` entries, `read_dense` reads a missing entry as NaN, a larger `Jac` view keeps
/// its outside and `inverse` of `ad` divides by zero.
#[cfg(not(debug_assertions))]
#[test]
fn out_of_domain_does_not_panic_in_release() {
    let t = RnTangent {
        rho: Vector([1.0, 2.0, 3.0]),
    };
    let (mut short, mut long) = ([9.0; 2], [9.0; 5]);
    t.write_dense(&mut short);
    t.write_dense(&mut long);
    assert_eq!(bits(&short), bits(&[1.0, 2.0]));
    assert_eq!(bits(&long), bits(&[1.0, 2.0, 3.0, 9.0, 9.0]));
    // The missing component is poisoned, not zeroed, so a wrongly sized buffer cannot pass for a
    // tangent whose last component happens to be zero.
    let from_short = RnTangent::<f64, 3>::read_dense(&short).rho.0;
    assert_eq!(bits(&from_short[..2]), bits(&[1.0, 2.0]));
    assert!(from_short[2].is_nan(), "a missing entry must be NaN");
    assert_eq!(
        bits(&RnTangent::<f64, 3>::read_dense(&[1.0, 2.0, 3.0, 4.0]).rho.0),
        bits(&[1.0, 2.0, 3.0])
    );
    let mut view = [f64::NAN; 16];
    RnJac::<f64, 3>::identity().write_dense(&mut StridedMut::col_major(&mut view, 4, 4));
    for (i, v) in view.iter().enumerate() {
        assert_eq!(v.is_nan(), i % 4 == 3 || i >= 12, "entry {i}");
    }
    let inv = G::<f64>::ad(&RnTangent::zero()).inverse();
    let mut d = [[0.0; 3]; 3];
    inv.write_dense(&mut StridedMut::col_major(d.as_flattened_mut(), 3, 3));
    assert_eq!(bits(d.as_flattened())[0], f64::INFINITY.to_bits());
}
