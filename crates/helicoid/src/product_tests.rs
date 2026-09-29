//! `Product` under the generic laws of [`crate::laws`], and its hand cases.
//!
//! Four groups, each with `Rn` and the test-only Heisenberg group of `heis_tests.rs` (the one
//! non-abelian factor until SO(3) lands, so the non-abelian block structure is exercised here
//! only on a 2-step nilpotent group): `Product<Rn<3>, Rn<2>>`, `Product<Rn<2>, Heis>` (the
//! non-abelian block at offset 2), `Product<Heis, Rn<2>>` (the non-abelian block first, so the
//! operand order of the first block of `Jac::mul` is seen) and the nested
//! `Product<Product<Heis, Rn<2>>, Heis>` (`D = 8`, blocks at offsets 0, 3 and 5).
//!
//! The bounds are recorded as described at [`Bounds`], one pair for the four groups: the
//! larger of their worst errors. The worst errors, in `u`: `f64` and `Dual<f64, D>` (same
//! samples, same values) 5.07 for `group_axioms`, 3.79 for `adjoint_identity`, 3.93 for
//! `plus_minus`, 4.38 for `jac_dense_order`, 5.69 for `sandwich_matches_dense`, 3.99 for
//! `tangent_dense_order` (`dot` is two partial sums) and 0 for the rest; `f32` 4.34, 3.91, 3.93,
//! 4.24, 3.74, 2.94 and 0 for the rest.

use crate::heis_tests::{heis_jac, HJac, Heis};
use crate::laws::{laws_for, Bounds, Sample};
use crate::{Jac, Left, LieGroup, Product, ProductJac, Right, Rn, RnJac, Tangent};
use core::array;
use helicoid_linalg::{Blend, Matrix, StridedMut, Vector};
use std::vec::Vec;

type P1<S> = Product<Rn<S, 3>, Rn<S, 2>>;
type P2<S> = Product<Rn<S, 2>, Heis<S>>;
type P3<S> = Product<Product<Heis<S>, Rn<S, 2>>, Heis<S>>;
type H2 = Product<Heis<f64>, Rn<f64, 2>>;
type H2S<S> = Product<Heis<S>, Rn<S, 2>>;
type H2Jac = ProductJac<HJac<f64>, RnJac<f64, 2>>;

// `k I` with `0.5 <= |k| < 1.5`: invertible, of either sign.
fn rn_jac<S: Sample, const N: usize>(x: f64, i: usize) -> RnJac<S, N> {
    RnJac::scalar(S::sample(x + 0.5_f64.copysign(x), i))
}

fn heis_at<S: Sample>(v: &[f64], o: usize) -> HJac<S> {
    heis_jac::<S>(&[v[o], v[o + 1], v[o + 2]])
}

fn p1_jac<S: Sample>(v: &[f64; 5]) -> ProductJac<RnJac<S, 3>, RnJac<S, 2>> {
    ProductJac(rn_jac(v[0], 0), rn_jac(v[3], 3))
}

fn p2_jac<S: Sample>(v: &[f64; 5]) -> ProductJac<RnJac<S, 2>, HJac<S>> {
    ProductJac(rn_jac(v[0], 0), heis_at(v, 2))
}

fn h2_jac<S: Sample>(v: &[f64; 5]) -> ProductJac<HJac<S>, RnJac<S, 2>> {
    ProductJac(heis_at(v, 0), rn_jac(v[3], 3))
}

type P3Jac<S> = ProductJac<ProductJac<HJac<S>, RnJac<S, 2>>, HJac<S>>;

fn p3_jac<S: Sample>(v: &[f64; 8]) -> P3Jac<S> {
    ProductJac(ProductJac(heis_at(v, 0), rn_jac(v[3], 3)), heis_at(v, 5))
}

const F64: Bounds = Bounds {
    axioms: 11.0,
    exp_log: 0.0,
    adjoint: 8.0,
    jl_ad_jr: 0.0,
    plus_minus: 8.0,
    rows: 0.0,
    ad: 0.0,
    sides: 0.0,
    tangent_order: 8.0,
    jac_order: 9.0,
    sandwich: 12.0,
};
const F32: Bounds = Bounds {
    axioms: 9.0,
    sandwich: 8.0,
    tangent_order: 6.0,
    ..F64
};

laws_for!(rn3_rn2, P1, p1_jac, F64, F32, 5);
laws_for!(rn2_heis, P2, p2_jac, F64, F32, 5);
laws_for!(heis_rn2, H2S, h2_jac, F64, F32, 5);
laws_for!(nested, P3, p3_jac, F64, F32, 8);

fn bits(x: &[f64]) -> Vec<u64> {
    x.iter().map(|v| v.to_bits()).collect()
}

/// The dense matrix of `j` in column-major order, as bits, through a NaN-poisoned view.
fn dense<const D: usize, T: Tangent<f64>, J: Jac<f64, T>>(j: &J) -> [[u64; D]; D] {
    let mut buf = [[f64::NAN; D]; D];
    j.write_dense(&mut StridedMut::col_major(buf.as_flattened_mut(), D, D));
    buf.map(|c| c.map(f64::to_bits))
}

/// `diag(a, b)` for column-major blocks, as bits: `+0` off the diagonal blocks.
fn block_diag<const A: usize, const B: usize, const D: usize>(
    a: [[u64; A]; A],
    b: [[u64; B]; B],
) -> [[u64; D]; D] {
    array::from_fn(|c| {
        array::from_fn(|r| match (r < A, c < A) {
            (true, true) => a[c][r],
            (false, false) => b[c - A][r - A],
            _ => 0.0_f64.to_bits(),
        })
    })
}

fn rn2(x: f64, y: f64) -> Rn<f64, 2> {
    Rn(Vector([x, y]))
}

fn heis(x: f64, y: f64, z: f64) -> Heis<f64> {
    Heis(Vector([x, y, z]))
}

#[test]
fn the_group_law_is_componentwise() {
    let a = Product(heis(1.0, 2.0, 3.0), rn2(1.0, -2.0));
    let b = Product(heis(4.0, 5.0, 6.0), rn2(0.5, 2.0));
    let ab = a * b;
    // The Heisenberg cross term stays in the first factor.
    assert_eq!(bits(&ab.0 .0 .0), bits(&[5.0, 7.0, 7.5]));
    assert_eq!(bits(&ab.1 .0 .0), bits(&[1.5, 0.0]));
    let i = H2::identity();
    assert_eq!(bits(&i.0 .0 .0), bits(&[0.0; 3]));
    assert_eq!(bits(&i.1 .0 .0), bits(&[0.0; 2]));
    let inv = a.inverse();
    assert_eq!(bits(&inv.0 .0 .0), bits(&[-1.0, -2.0, -3.0]));
    assert_eq!(bits(&inv.1 .0 .0), bits(&[-1.0, 2.0]));
    // `Exp` and `Log` on the pair, and the two sides of `⊕`, factor by factor.
    let tau = a.log();
    assert_eq!(bits(&tau.0.rho.0), bits(&[1.0, 2.0, 3.0]));
    assert_eq!(bits(&tau.1.rho.0), bits(&[1.0, -2.0]));
    assert_eq!(bits(&H2::exp(&tau).1 .0 .0), bits(&[1.0, -2.0]));
    let tb = b.log();
    let (r, l) = (a.rplus(&tb), a.lplus(&tb));
    assert_eq!(bits(&r.0 .0 .0), bits(&[5.0, 7.0, 7.5]));
    assert_eq!(bits(&l.0 .0 .0), bits(&[5.0, 7.0, 10.5]));
    assert_eq!(bits(&r.1 .0 .0), bits(&[1.5, 0.0]));
    let (mr, ml) = (b.rminus(&a), b.lminus(&a));
    assert_eq!(bits(&mr.0.rho.0), bits(&[3.0, 3.0, 4.5]));
    assert_eq!(bits(&ml.0.rho.0), bits(&[3.0, 3.0, 1.5]));
    assert_eq!(bits(&mr.1.rho.0), bits(&[-0.5, 4.0]));
}

#[test]
fn dof_is_the_sum() {
    assert_eq!(P1::<f64>::DOF, 5);
    assert_eq!(
        <<P1<f64> as LieGroup<f64>>::Tangent as Tangent<f64>>::DOF,
        5
    );
    assert_eq!(P2::<f64>::DOF, 5);
    assert_eq!(P3::<f64>::DOF, 8);
    assert_eq!(
        <<P3<f64> as LieGroup<f64>>::Tangent as Tangent<f64>>::DOF,
        8
    );
}

#[test]
fn the_tangent_is_a_then_b() {
    type T = <P3<f64> as LieGroup<f64>>::Tangent;
    let src: [f64; 8] = array::from_fn(|i| (i + 1) as f64);
    let t = T::read_dense(&src);
    assert_eq!(bits(&t.0 .0.rho.0), bits(&[1.0, 2.0, 3.0]));
    assert_eq!(bits(&t.0 .1.rho.0), bits(&[4.0, 5.0]));
    assert_eq!(bits(&t.1.rho.0), bits(&[6.0, 7.0, 8.0]));
    let mut out = [f64::NAN; 8];
    t.write_dense(&mut out);
    assert_eq!(bits(&out), bits(&src));
}

#[test]
fn dot_is_nested_partial_sums() {
    type T = <P3<f64> as LieGroup<f64>>::Tangent;
    // Cancellation makes the grouping visible: `((a0 + a1) + a2) + (a3 + a4)` and `(a5 + a6) + a7`
    // are the factors' sums, and the pair adds them; one flat sum would give exactly 4.
    let a = [1e16, 1.0, 1.0, -1e16, 1.0, 1.0, 1.0, 1.0];
    let (t, ones) = (T::read_dense(&a), T::read_dense(&[1.0; 8]));
    let nested = (((a[0] + a[1]) + a[2]) + (a[3] + a[4])) + ((a[5] + a[6]) + a[7]);
    assert_eq!(t.dot(&ones).to_bits(), nested.to_bits());
    assert_ne!(nested.to_bits(), 4.0_f64.to_bits());
}

#[test]
fn write_dense_is_block_diagonal() {
    let x = heis(1.0, 2.0, 3.0);
    // `Ad_x = I + ad_x` is not symmetric: rows `[1 0 0]`, `[0 1 0]`, `[-2 1 1]`.
    let ad = dense::<3, _, _>(&x.adjoint());
    assert_eq!(
        ad,
        [[1.0, 0.0, -2.0], [0.0, 1.0, 1.0], [0.0, 0.0, 1.0]].map(|c| c.map(f64::to_bits))
    );
    let j = ProductJac(RnJac::<f64, 2>::scalar(2.0), x.adjoint());
    let want: [[f64; 5]; 5] = [
        [2.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 2.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, -2.0],
        [0.0, 0.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 0.0, 0.0, 1.0],
    ];
    assert_eq!(dense::<5, _, _>(&j), want.map(|c| c.map(f64::to_bits)));
    // A padded row-major view: the block is written and its padding is not.
    let mut padded = [f64::NAN; 35];
    j.write_dense(&mut StridedMut::with_strides(&mut padded, 5, 5, 7, 1));
    for (i, v) in padded.iter().enumerate() {
        let (r, c) = (i / 7, i % 7);
        if c >= 5 {
            assert!(v.is_nan(), "padding {i} was written");
        } else {
            assert_eq!(v.to_bits(), want[c][r].to_bits(), "entry {i}");
        }
    }
    // Nested: `diag(diag(Ad_x, 3 I_2), Ad_x)`.
    let n = ProductJac(
        ProductJac(x.adjoint(), RnJac::<f64, 2>::scalar(3.0)),
        x.adjoint(),
    );
    let (a, k, e) = (
        dense::<3, _, _>(&x.adjoint()),
        dense::<2, _, _>(&RnJac::<f64, 2>::scalar(3.0)),
        dense::<8, _, _>(&n),
    );
    assert_eq!(e, block_diag::<5, 3, 8>(block_diag::<3, 2, 5>(a, k), a));
}

#[test]
fn every_row_is_the_factors_rows() {
    let (x, y) = (
        Product(heis(1.0, 2.0, 3.0), rn2(1.0, -2.0)),
        Product(heis(4.0, 5.0, 6.0), rn2(0.5, 2.0)),
    );
    let tau = y.log();
    let bd = |p: &(H2Jac, H2Jac)| (dense::<5, _, _>(&p.0), dense::<5, _, _>(&p.1));
    let factors = |a: (HJac<f64>, HJac<f64>), b: (RnJac<f64, 2>, RnJac<f64, 2>)| {
        (
            block_diag::<3, 2, 5>(dense::<3, _, _>(&a.0), dense::<2, _, _>(&b.0)),
            block_diag::<3, 2, 5>(dense::<3, _, _>(&a.1), dense::<2, _, _>(&b.1)),
        )
    };
    let (hx, hy, ht) = (x.0, y.0, tau.0);
    let (rx, ry, rt) = (x.1, y.1, tau.1);
    assert_eq!(
        bd(&x.rplus_jacobians(&tau)),
        factors(hx.rplus_jacobians(&ht), rx.rplus_jacobians(&rt))
    );
    assert_eq!(
        bd(&x.lplus_jacobians(&tau)),
        factors(hx.lplus_jacobians(&ht), rx.lplus_jacobians(&rt))
    );
    assert_eq!(
        bd(&y.rminus_jacobians(&x)),
        factors(hy.rminus_jacobians(&hx), ry.rminus_jacobians(&rx))
    );
    assert_eq!(
        bd(&y.lminus_jacobians(&x)),
        factors(hy.lminus_jacobians(&hx), ry.lminus_jacobians(&rx))
    );
    assert_eq!(
        bd(&x.compose_jacobians::<Right>(&y)),
        factors(
            hx.compose_jacobians::<Right>(&hy),
            rx.compose_jacobians::<Right>(&ry)
        )
    );
    assert_eq!(
        bd(&x.compose_jacobians::<Left>(&y)),
        factors(
            hx.compose_jacobians::<Left>(&hy),
            rx.compose_jacobians::<Left>(&ry)
        )
    );
    // The sides differ on the non-abelian factor, so the two rows above are not the same test.
    assert_ne!(
        bd(&x.compose_jacobians::<Right>(&y)).0,
        bd(&x.compose_jacobians::<Left>(&y)).0
    );
    let single = |a: HJac<f64>, b: RnJac<f64, 2>| {
        block_diag::<3, 2, 5>(dense::<3, _, _>(&a), dense::<2, _, _>(&b))
    };
    assert_eq!(
        dense::<5, _, _>(&x.inverse_jacobian::<Right>()),
        single(
            hx.inverse_jacobian::<Right>(),
            rx.inverse_jacobian::<Right>()
        )
    );
    assert_eq!(
        dense::<5, _, _>(&x.inverse_jacobian::<Left>()),
        single(hx.inverse_jacobian::<Left>(), rx.inverse_jacobian::<Left>())
    );
    // `Ad`, `ad`, `J_r`, `J_l` and their inverses.
    assert_eq!(
        dense::<5, _, _>(&x.adjoint()),
        single(hx.adjoint(), rx.adjoint())
    );
    let ad = H2::ad(&tau);
    assert_eq!(dense::<5, _, _>(&ad), single(Heis::ad(&ht), Rn::ad(&rt)));
    let jr = H2::jr(&tau);
    assert_eq!(dense::<5, _, _>(&jr), single(Heis::jr(&ht), Rn::jr(&rt)));
    let jl = H2::jl(&tau);
    assert_eq!(dense::<5, _, _>(&jl), single(Heis::jl(&ht), Rn::jl(&rt)));
    let jri = H2::jr_inv(&tau);
    assert_eq!(
        dense::<5, _, _>(&jri),
        single(Heis::jr_inv(&ht), Rn::jr_inv(&rt))
    );
    let jli = H2::jl_inv(&tau);
    assert_eq!(
        dense::<5, _, _>(&jli),
        single(Heis::jl_inv(&ht), Rn::jl_inv(&rt))
    );
}

/// The dense product `a b` of column-major `n x n` matrices, in the order `a` then `b`.
fn matmul<const N: usize>(a: [[u64; N]; N], b: [[u64; N]; N]) -> [[u64; N]; N] {
    let f = |x: u64| f64::from_bits(x);
    array::from_fn(|c| {
        array::from_fn(|r| {
            (0..N)
                .map(|k| f(a[k][r]) * f(b[c][k]))
                .sum::<f64>()
                .to_bits()
        })
    })
}

#[test]
fn mul_is_the_dense_product_in_operand_order() {
    // Dyadic entries, so every sum is exact and the bits are the dense product's in any order.
    let (a, b) = ([0.5, -1.0, 0.25], [1.0, 0.5, -0.75]);
    let ja = ProductJac(heis_jac::<f64>(&a), RnJac::<f64, 2>::scalar(2.0));
    let jb = ProductJac(heis_jac::<f64>(&b), RnJac::<f64, 2>::scalar(3.0));
    let (da, db) = (dense::<5, _, _>(&ja), dense::<5, _, _>(&jb));
    let (ab, ba) = (matmul(da, db), matmul(db, da));
    assert_ne!(ab, ba, "the fixture commutes");
    assert_eq!(dense::<5, _, _>(&ja.mul(&jb)), ab);
    assert_eq!(dense::<5, _, _>(&jb.mul(&ja)), ba);
    // The same on the second block: `Product<Rn<2>, Heis>`.
    let (ka, kb) = (
        ProductJac(RnJac::<f64, 2>::scalar(2.0), heis_jac::<f64>(&a)),
        ProductJac(RnJac::<f64, 2>::scalar(3.0), heis_jac::<f64>(&b)),
    );
    let (ea, eb) = (dense::<5, _, _>(&ka), dense::<5, _, _>(&kb));
    assert_eq!(dense::<5, _, _>(&ka.mul(&kb)), matmul(ea, eb));
    assert_eq!(dense::<5, _, _>(&kb.mul(&ka)), matmul(eb, ea));
}

#[test]
fn sandwich_scales_each_block_pair() {
    // `diag(2 I_2, 3 I_3)` on all ones: `2 * 2`, `2 * 3` across and `3 * 3`.
    let j = ProductJac(RnJac::<f64, 2>::scalar(2.0), RnJac::<f64, 3>::scalar(3.0));
    let ones = Matrix::from_cols([Vector([1.0; 5]); 5]);
    let s = j.sandwich::<5>(&ones);
    for c in 0..5 {
        for r in 0..5 {
            let k = |i: usize| if i < 2 { 2.0_f64 } else { 3.0 };
            assert_eq!(s.get(r, c).to_bits(), (k(r) * k(c)).to_bits(), "({r}, {c})");
        }
    }
}

#[test]
fn blend_selects_each_value_type() {
    let a = Product(rn2(1.0, 1.0), heis(1.0, 2.0, 3.0));
    let b = Product(rn2(2.0, 2.0), heis(4.0, 5.0, 6.0));
    let (ta, tb) = (a.log(), b.log());
    let (ja, jb) = (a.adjoint(), ProductJac(RnJac::scalar(-1.0), b.1.adjoint()));
    for (m, x, t, j) in [(true, a, ta, ja), (false, b, tb, jb)] {
        let g = Product::blend(m, a, b);
        assert_eq!(bits(&g.0 .0 .0), bits(&x.0 .0 .0));
        assert_eq!(bits(&g.1 .0 .0), bits(&x.1 .0 .0));
        let u = <(_, _)>::blend(m, ta, tb);
        assert_eq!(bits(&u.0.rho.0), bits(&t.0.rho.0));
        assert_eq!(bits(&u.1.rho.0), bits(&t.1.rho.0));
        let k = ProductJac::blend(m, ja, jb);
        assert_eq!(dense::<5, _, _>(&k), dense::<5, _, _>(&j));
    }
}

/// A domain is a `debug_assert!` (D11): the debug profile that `just test` runs panics ...
#[cfg(debug_assertions)]
mod out_of_domain {
    use super::*;

    #[test]
    #[should_panic(expected = "Tangent::write_dense")]
    fn write_dense_of_the_wrong_length() {
        <P1<f64> as LieGroup<f64>>::Tangent::zero().write_dense(&mut [0.0; 6]);
    }

    #[test]
    #[should_panic(expected = "Tangent::read_dense")]
    fn read_dense_of_the_wrong_length() {
        let _ = <P1<f64> as LieGroup<f64>>::Tangent::read_dense(&[0.0; 4]);
    }

    #[test]
    #[should_panic(expected = "Jac::write_dense")]
    fn write_dense_into_a_larger_view() {
        let mut buf = [0.0; 36];
        let j: <P1<f64> as LieGroup<f64>>::Jac = Jac::identity();
        j.write_dense(&mut StridedMut::col_major(&mut buf, 6, 6));
    }
}

/// ... and a release build (`cargo nextest run --release`) never panics: `write_dense` writes
/// `min(len, DOF)` entries and `read_dense` reads a missing entry as `+0`.
#[cfg(not(debug_assertions))]
#[test]
fn out_of_domain_does_not_panic_in_release() {
    type T = <P1<f64> as LieGroup<f64>>::Tangent;
    let t = T::read_dense(&[1.0, 2.0, 3.0, 4.0, 5.0]);
    let (mut short, mut long) = ([9.0; 4], [9.0; 7]);
    t.write_dense(&mut short);
    t.write_dense(&mut long);
    assert_eq!(bits(&short), bits(&[1.0, 2.0, 3.0, 4.0]));
    assert_eq!(bits(&long), bits(&[1.0, 2.0, 3.0, 4.0, 5.0, 9.0, 9.0]));
    let s = T::read_dense(&short);
    assert_eq!(bits(&s.0.rho.0), bits(&[1.0, 2.0, 3.0]));
    assert_eq!(bits(&s.1.rho.0), bits(&[4.0, 0.0]));
    let l = T::read_dense(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    assert_eq!(bits(&l.1.rho.0), bits(&[4.0, 5.0]));
    // A slice shorter than the first factor: `min(len, DOF)` entries written, the rest read as +0.
    let mut tiny = [9.0; 2];
    t.write_dense(&mut tiny);
    assert_eq!(bits(&tiny), bits(&[1.0, 2.0]));
    let z = T::read_dense(&[1.0]);
    assert_eq!(bits(&z.0.rho.0), bits(&[1.0, 0.0, 0.0]));
    assert_eq!(bits(&z.1.rho.0), bits(&[0.0, 0.0]));
}
