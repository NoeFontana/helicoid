//! `Product` under the generic laws of [`crate::laws`], its hand cases and the `sandwich` twin.
//!
//! Four groups, each with `Rn` and the test-only Heisenberg group of `heis_tests.rs` (the one
//! non-abelian factor until SO(3) lands, so the non-abelian block structure is exercised here
//! only on a 2-step nilpotent group): `Product<Rn<3>, Rn<2>>`, `Product<Rn<2>, Heis>` (the
//! non-abelian block at offset 2), `Product<Heis, Rn<2>>` (the non-abelian block first, so the
//! operand order of the first block of `Jac::mul` is seen) and the nested
//! `Product<Product<Heis, Rn<2>>, Heis>` (`D = 8`, blocks at offsets 0, 3 and 5).
//!
//! The bounds are recorded as described at [`Bounds`], one pair for the four groups: the larger of
//! their worst errors. The worst errors, in `u`: `f64` and `Dual<f64, D>` (same samples, same
//! values) 5.07 for `group_axioms`, 3.79 for `adjoint_identity`, 3.93 for `plus_minus`, 4.38 for
//! `jac_dense_order`, 5.60 for `sandwich_matches_dense` and for
//! `productjac_sandwich_matches_reference`, and 0 for the rest; `f32` 4.34, 3.91, 3.93, 4.38, 4.21,
//! 5.28, and 3.38 for `tangent_dense_order` and 0 for the rest.
//!
//! `tangent_dense_order` is `0` for `f64` and `Dual`, as it is for every other group: `dot_acc`
//! threads one flat sum over the concatenated dense order rather than adding the factors' partial
//! sums, which is why `0025` made it the required operation.

use crate::heis_tests::{heis_jac_at, HJac, Heis};
use crate::laws::{cov, dense_bits, e, laws_for, sample, within, worst, Bounds, Sample};
use crate::reference::productjac_sandwich;
use crate::{Jac, Left, LieGroup, Product, ProductJac, Right, Rn, RnJac, Tangent};
use core::array;
use helicoid_linalg::{Blend, Matrix, Real, StridedMut, Vector};
use std::vec::Vec;

type P1<S> = Product<Rn<S, 3>, Rn<S, 2>>;
type P2<S> = Product<Rn<S, 2>, Heis<S>>;
type H2<S> = Product<Heis<S>, Rn<S, 2>>;
type P3<S> = Product<Product<Heis<S>, Rn<S, 2>>, Heis<S>>;
type H2Jac<S> = ProductJac<HJac<S>, RnJac<S, 2>>;
/// An `H2` Jacobian beside the two factors' Jacobians it must be the block diagonal of.
type Blocked = (H2Jac<f64>, (HJac<f64>, RnJac<f64, 2>));
type P3Jac<S> = ProductJac<ProductJac<HJac<S>, RnJac<S, 2>>, HJac<S>>;

/// `k I` with `0.5 <= |k| < 1.5`: invertible, of either sign.
///
/// `lane` is the block's offset in the dense order, which is also the `Dual` variable the sample
/// becomes: two blocks of one fixture then never claim a lane while holding different values.
fn rn_jac<S: Sample, const N: usize>(x: f64, lane: usize) -> RnJac<S, N> {
    RnJac::scalar(S::sample(x + 0.5_f64.copysign(x), lane))
}

/// The Heisenberg block whose dense order starts at `o`: its samples and its `Dual` lanes are `o`,
/// `o + 1` and `o + 2`.
fn heis_at<S: Sample>(v: &[f64], o: usize) -> HJac<S> {
    heis_jac_at::<S>(&[v[o], v[o + 1], v[o + 2]], o)
}

fn p1_jac<S: Sample>(v: &[f64; 5]) -> ProductJac<RnJac<S, 3>, RnJac<S, 2>> {
    ProductJac(rn_jac(v[0], 0), rn_jac(v[3], 3))
}

fn p2_jac<S: Sample>(v: &[f64; 5]) -> ProductJac<RnJac<S, 2>, HJac<S>> {
    ProductJac(rn_jac(v[0], 0), heis_at(v, 2))
}

fn h2_jac<S: Sample>(v: &[f64; 5]) -> H2Jac<S> {
    ProductJac(heis_at(v, 0), rn_jac(v[3], 3))
}

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
    tangent_order: 0.0,
    jac_order: 9.0,
    sandwich: 12.0,
    // `PHASE3.md` §8's second check, measured 0.316 over the four product groups, 10 000 draws each; `Rn x Rn` is exactly `0`.
    dual_rows: 1.0,
    // `PHASE4.md` §1 and §3's six legs, over 10^6 draws of `laws::Rng`'s own geodesic stream (`measure_geodesic`), twice the worst of the three scalars, rounded up: the worst of the four groups: symmetry 5.038 / 6.345, left 5.993 / 7.280, right 7.430 / 8.821. The `velocity` leg is exactly 0 on every group -- it is `rminus` against itself until a group overrides -- and `t=0` is `gerr`'s floor for two **bitwise equal** elements, not a geodesic error: 0 for `Rn x Rn`, 0 for the Heisenberg factors too.
    geodesic: 18.0,
};
const F32: Bounds = Bounds {
    axioms: 9.0,
    tangent_order: 7.0,
    sandwich: 9.0,
    dual_rows: 1.0,
    ..F64
};

laws_for!(rn3_rn2, P1, p1_jac, 5, F64, F32);
laws_for!(rn2_heis, P2, p2_jac, 5, F64, F32);
laws_for!(heis_rn2, H2, h2_jac, 5, F64, F32);
laws_for!(nested, P3, p3_jac, 8, F64, F32);

/// The bound of `productjac_sandwich_matches_reference`, in `u`, recorded as [`Bounds`] describes;
/// `TWIN32` is the `f32` figure, `TWIN64` the `f64` and `Dual` one.
const TWIN64: f64 = 12.0;
const TWIN32: f64 = 11.0;

/// `sandwich` against [`productjac_sandwich`], the dense `D⁴` product written into caller memory,
/// as a norm-wise error in units of `u` on the value part.
///
/// The `Dual` derivative lanes of `sandwich` are not compared here: `dual_value_is_plain_value`
/// pins the value path bit for bit, and the lanes are what `jacobians_match_dual_*` owes
/// (`PHASE3.md` §8), which is not started.
fn twin_err<S, TA, TB, JA, JB, const D: usize>(j: &ProductJac<JA, JB>, c: &Matrix<S, D, D>) -> f64
where
    S: Real,
    TA: Tangent<S>,
    TB: Tangent<S>,
    JA: Jac<S, TA>,
    JB: Jac<S, TB>,
{
    // NaN-poisoned, as the harness's own readers are: an entry the twin leaves unwritten fails
    // the comparison instead of passing as a structural zero.
    let mut twin = [[S::zero() / S::zero(); D]; D];
    productjac_sandwich::<S, TA, TB, JA, JB, D>(
        j,
        c,
        &mut StridedMut::row_major(twin.as_flattened_mut(), D, D),
    );
    let got = j.sandwich(c);
    let fast: [[f64; D]; D] = array::from_fn(|r| array::from_fn(|k| got.get(r, k).value_f64()));
    let want = twin.map(|row| row.map(S::value_f64));
    e::<S>(fast.as_flattened(), want.as_flattened())
}

/// `productjac_sandwich_matches_reference` for the groups of one dimension, over `f64`, `f32` and
/// `Dual<f64, $D>`, on the covariance of [`cov`] and the `Jac` fixtures above.
macro_rules! twins_for {
    ($m:ident, $D:literal, $($jac:ident),+) => {
        mod $m {
            use super::*;
            twins_for!(@case as_f64, f64, TWIN64, $D, $($jac),+);
            twins_for!(@case as_f32, f32, TWIN32, $D, $($jac),+);
            twins_for!(@case as_dual, helicoid_linalg::Dual<f64, $D>, TWIN64, $D, $($jac),+);
        }
    };
    (@case $n:ident, $S:ty, $B:ident, $D:literal, $($jac:ident),+) => {
        mod $n {
            use super::*;
            proptest::proptest! {
                #[test]
                fn productjac_sandwich_matches_reference(
                    a in sample::<$D>(), b in sample::<$D>(), c in sample::<$D>(),
                ) {
                    let sigma = cov::<$S, $D>(&a, &b, &c);
                    let mut w = 0.0;
                    $( w = worst(w, twin_err::<$S, _, _, _, _, $D>(&$jac::<$S>(&a), &sigma)); )+
                    within(w, $B)?;
                }
            }
        }
    };
}

twins_for!(twin_d5, 5, p1_jac, p2_jac, h2_jac);
twins_for!(twin_d8, 8, p3_jac);

fn bits(x: &[f64]) -> Vec<u64> {
    x.iter().map(|v| v.to_bits()).collect()
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

/// A `Product` of two `Heis`-shaped tangents, for a hand case that needs one by value.
fn h2(a: Heis<f64>, b: Rn<f64, 2>) -> H2<f64> {
    Product(a, b)
}

#[test]
fn the_group_law_is_componentwise() {
    let a = h2(heis(1.0, 2.0, 3.0), rn2(1.0, -2.0));
    let b = h2(heis(4.0, 5.0, 6.0), rn2(0.5, 2.0));
    let ab = a * b;
    // The Heisenberg cross term stays in the first factor.
    assert_eq!(bits(&ab.0 .0 .0), bits(&[5.0, 7.0, 7.5]));
    assert_eq!(bits(&ab.1 .0 .0), bits(&[1.5, 0.0]));
    let i = H2::<f64>::identity();
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
fn dot_is_one_flat_sum() {
    type T = <P3<f64> as LieGroup<f64>>::Tangent;
    // Cancellation makes the grouping visible. Flat, the four `1.0`s before `-1e16` are lost and
    // the three after it survive, for `4`; the factors' partial sums — `((a0 + a1) + a2) +
    // (a3 + a4)`, then `(a5 + a6) + a7` — keep `a4` instead, for `3`. `0025` rejects the second.
    let a = [1e16, 1.0, 1.0, -1e16, 1.0, 1.0, 1.0, 1.0];
    let (t, ones) = (T::read_dense(&a), T::read_dense(&[1.0; 8]));
    let flat = a.iter().fold(0.0_f64, |s, x| s + x);
    assert_eq!(t.dot(&ones).to_bits(), flat.to_bits());
    assert_eq!(flat.to_bits(), 4.0_f64.to_bits());
    assert_ne!(flat.to_bits(), 3.0_f64.to_bits());
    // The accumulator is threaded into `A` first, so a seed joins the same left-to-right sum.
    let seeded = a.iter().fold(-1e16_f64, |s, x| s + x);
    assert_eq!(t.dot_acc(&ones, -1e16).to_bits(), seeded.to_bits());
}

#[test]
fn write_dense_is_block_diagonal() {
    let x = heis(1.0, 2.0, 3.0);
    // `Ad_x = I + ad_x` is not symmetric: rows `[1 0 0]`, `[0 1 0]`, `[-2 1 1]`.
    let ad = dense_bits::<_, _, 3>(&x.adjoint());
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
    assert_eq!(dense_bits::<_, _, 5>(&j), want.map(|c| c.map(f64::to_bits)));
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
        dense_bits::<_, _, 3>(&x.adjoint()),
        dense_bits::<_, _, 2>(&RnJac::<f64, 2>::scalar(3.0)),
        dense_bits::<_, _, 8>(&n),
    );
    assert_eq!(e, block_diag::<5, 3, 8>(block_diag::<3, 2, 5>(a, k), a));
}

#[test]
fn every_row_is_the_factors_rows() {
    let (x, y) = (
        h2(heis(1.0, 2.0, 3.0), rn2(1.0, -2.0)),
        h2(heis(4.0, 5.0, 6.0), rn2(0.5, 2.0)),
    );
    let tau = y.log();
    let bd =
        |p: &(H2Jac<f64>, H2Jac<f64>)| (dense_bits::<_, _, 5>(&p.0), dense_bits::<_, _, 5>(&p.1));
    let factors = |a: (HJac<f64>, HJac<f64>), b: (RnJac<f64, 2>, RnJac<f64, 2>)| {
        (
            block_diag::<3, 2, 5>(dense_bits::<_, _, 3>(&a.0), dense_bits::<_, _, 2>(&b.0)),
            block_diag::<3, 2, 5>(dense_bits::<_, _, 3>(&a.1), dense_bits::<_, _, 2>(&b.1)),
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
        block_diag::<3, 2, 5>(dense_bits::<_, _, 3>(&a), dense_bits::<_, _, 2>(&b))
    };
    assert_eq!(
        dense_bits::<_, _, 5>(&x.inverse_jacobian::<Right>()),
        single(
            hx.inverse_jacobian::<Right>(),
            rx.inverse_jacobian::<Right>()
        )
    );
    assert_eq!(
        dense_bits::<_, _, 5>(&x.inverse_jacobian::<Left>()),
        single(hx.inverse_jacobian::<Left>(), rx.inverse_jacobian::<Left>())
    );
    // `Ad`, `ad`, `J_r`, `J_l` and their inverses.
    assert_eq!(
        dense_bits::<_, _, 5>(&x.adjoint()),
        single(hx.adjoint(), rx.adjoint())
    );
    let rows: [Blocked; 5] = [
        (H2::ad(&tau), (Heis::ad(&ht), Rn::ad(&rt))),
        (H2::jr(&tau), (Heis::jr(&ht), Rn::jr(&rt))),
        (H2::jl(&tau), (Heis::jl(&ht), Rn::jl(&rt))),
        (H2::jr_inv(&tau), (Heis::jr_inv(&ht), Rn::jr_inv(&rt))),
        (H2::jl_inv(&tau), (Heis::jl_inv(&ht), Rn::jl_inv(&rt))),
    ];
    for (whole, (a, b)) in rows {
        assert_eq!(dense_bits::<_, _, 5>(&whole), single(a, b));
    }
}

/// The dense product `a b` of column-major `n x n` bit matrices, through `Matrix`'s own
/// left-to-right sum, so the test's reference is not a second multiplication routine.
fn matmul<const N: usize>(a: [[u64; N]; N], b: [[u64; N]; N]) -> [[u64; N]; N] {
    let m = |d: [[u64; N]; N]| Matrix::from_cols(d.map(|c| Vector(c.map(f64::from_bits))));
    let p = m(a) * m(b);
    array::from_fn(|c| array::from_fn(|r| p.get(r, c).to_bits()))
}

#[test]
fn mul_is_the_dense_product_in_operand_order() {
    // Dyadic entries, so every sum is exact and the bits are the dense product's in any order.
    let (a, b) = ([0.5, -1.0, 0.25], [1.0, 0.5, -0.75]);
    let ja = ProductJac(heis_jac_at::<f64>(&a, 0), RnJac::<f64, 2>::scalar(2.0));
    let jb = ProductJac(heis_jac_at::<f64>(&b, 0), RnJac::<f64, 2>::scalar(3.0));
    let (da, db) = (dense_bits::<_, _, 5>(&ja), dense_bits::<_, _, 5>(&jb));
    let (ab, ba) = (matmul(da, db), matmul(db, da));
    assert_ne!(ab, ba, "the fixture commutes");
    assert_eq!(dense_bits::<_, _, 5>(&ja.mul(&jb)), ab);
    assert_eq!(dense_bits::<_, _, 5>(&jb.mul(&ja)), ba);
    // The same on the second block: `Product<Rn<2>, Heis>`.
    let (ka, kb) = (
        ProductJac(RnJac::<f64, 2>::scalar(2.0), heis_jac_at::<f64>(&a, 0)),
        ProductJac(RnJac::<f64, 2>::scalar(3.0), heis_jac_at::<f64>(&b, 0)),
    );
    let (ea, eb) = (dense_bits::<_, _, 5>(&ka), dense_bits::<_, _, 5>(&kb));
    assert_eq!(dense_bits::<_, _, 5>(&ka.mul(&kb)), matmul(ea, eb));
    assert_eq!(dense_bits::<_, _, 5>(&kb.mul(&ka)), matmul(eb, ea));
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
        assert_eq!(dense_bits::<_, _, 5>(&k), dense_bits::<_, _, 5>(&j));
    }
}

/// A domain is a `debug_assert!` (D11): the debug profile that `just test` runs panics, including
/// on a `Jac` view that is not `DOF x DOF` in either direction ...
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

    #[test]
    #[should_panic(expected = "Jac::write_dense")]
    fn write_dense_into_a_smaller_view() {
        let mut buf = [0.0; 16];
        let j: <P1<f64> as LieGroup<f64>>::Jac = Jac::identity();
        j.write_dense(&mut StridedMut::col_major(&mut buf, 4, 4));
    }
}

/// ... and a release build (`cargo nextest run --release`) never panics: `write_dense` writes
/// `min(len, DOF)` entries, `read_dense` reads a missing entry as NaN in either factor, and a
/// larger `Jac` view keeps its outside.
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
    // The missing component is poisoned, not zeroed, so a wrongly sized buffer cannot pass for a
    // tangent whose last component happens to be zero (`0025`).
    let s = T::read_dense(&short);
    assert_eq!(bits(&s.0.rho.0), bits(&[1.0, 2.0, 3.0]));
    assert_eq!(bits(&s.1.rho.0[..1]), bits(&[4.0]));
    assert!(s.1.rho.0[1].is_nan(), "a missing entry must be NaN");
    // A slice shorter than the first factor: `A` keeps what fits, `B` is wholly poisoned.
    let mut tiny = [9.0; 2];
    t.write_dense(&mut tiny);
    assert_eq!(bits(&tiny), bits(&[1.0, 2.0]));
    let z = T::read_dense(&[1.0]);
    assert_eq!(bits(&z.0.rho.0[..1]), bits(&[1.0]));
    assert!(z.0.rho.0[1..].iter().all(|x| x.is_nan()), "A is poisoned");
    assert!(z.1.rho.0.iter().all(|x| x.is_nan()), "B is poisoned");
    // A longer slice drops what is past `DOF`.
    let l = T::read_dense(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    assert_eq!(bits(&l.1.rho.0), bits(&[4.0, 5.0]));
    // A larger `Jac` view: the top-left `DOF x DOF` is written and the rest keeps its NaN.
    let mut view = [f64::NAN; 36];
    let j: <P1<f64> as LieGroup<f64>>::Jac = Jac::identity();
    j.write_dense(&mut StridedMut::col_major(&mut view, 6, 6));
    for (i, v) in view.iter().enumerate() {
        let (r, c) = (i % 6, i / 6);
        assert_eq!(v.is_nan(), r >= 5 || c >= 5, "entry {i}");
    }
}

/// A smaller `Jac` view is the one documented panic class (D11): without the `debug_assert!` the
/// strided access itself refuses, so the write neither wraps nor lands outside the view.
#[cfg(not(debug_assertions))]
#[test]
#[should_panic(expected = "StridedMut::set: index out of range")]
fn a_smaller_jac_view_panics_in_the_strided_access() {
    let mut buf = [0.0; 16];
    let j: <P1<f64> as LieGroup<f64>>::Jac = Jac::identity();
    j.write_dense(&mut StridedMut::col_major(&mut buf, 4, 4));
}
