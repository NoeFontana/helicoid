//! Law checks written once over `G: LieGroup<S>`; every group's tests call them
//! (`docs/PHASE3.md` §9). `D` is `G::DOF`, spelled as a parameter because stable Rust cannot size
//! an array from an associated const (`0005`); a law asserts `D == G::DOF == Tangent::DOF` at
//! build time, so a group whose tests call one has its `DOF` tie checked.
//!
//! A law returns its worst error in units of `u` (`NUMERICS.md` §11): a norm-wise difference over
//! `max(‖reference‖, 1)`, NaN propagating, compared by the caller with a bound recorded from a
//! measured run. Group elements are compared through `Log(a ⊖_R b)`: a group has no equality (R4),
//! so `Log` is trusted wherever a group's tests do not check it first.

use crate::{Jac, Left, LieGroup, Right, Side, Tangent};
use core::array;
use helicoid_linalg::{Dual, Matrix, Precision, Real, StridedMut};
use proptest::prelude::*;
use std::{format, vec::Vec};

/// A scalar a test can build from a sample: `Dual` gets one variable per component.
pub(crate) trait Sample: Real {
    fn sample(v: f64, i: usize) -> Self;
}

impl Sample for f64 {
    fn sample(v: f64, _: usize) -> Self {
        v
    }
}

impl Sample for f32 {
    fn sample(v: f64, _: usize) -> Self {
        v as f32
    }
}

impl<const N: usize> Sample for Dual<f64, N> {
    fn sample(v: f64, i: usize) -> Self {
        Dual::variable(v, i)
    }
}

/// The tangent whose dense components are the samples `v`.
pub(crate) fn tangent<S: Sample, G: LieGroup<S>, const D: usize>(v: &[f64; D]) -> G::Tangent {
    let s: [S; D] = array::from_fn(|i| S::sample(v[i], i));
    <G::Tangent as Tangent<S>>::read_dense(&s)
}

fn unit<S: Real>() -> f64 {
    match S::PRECISION {
        Precision::F64 => f64::EPSILON / 2.0,
        Precision::F32 => f64::from(f32::EPSILON) / 2.0,
    }
}

/// The larger error; NaN wins.
fn worst(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// `‖a - b‖ / max(‖b‖, 1)` in units of `u`.
fn e<S: Real>(a: &[f64], b: &[f64]) -> f64 {
    let diff = a.iter().zip(b).map(|(x, y)| (x - y) * (x - y));
    diff.sum::<f64>().sqrt() / norm(b).max(1.0) / unit::<S>()
}

fn nan<S: Real>() -> S {
    S::zero() / S::zero()
}

/// The dense components of `t`, from a NaN-poisoned buffer so that an unwritten entry fails a law.
fn dt<S: Real, G: LieGroup<S>, const D: usize>(t: &G::Tangent) -> [f64; D] {
    const { assert!(D == G::DOF && D == <G::Tangent as Tangent<S>>::DOF) };
    let mut buf = [nan::<S>(); D];
    t.write_dense(&mut buf);
    buf.map(S::value_f64)
}

/// The dense matrix of `j`, column-major (`m[c][r]`), through a NaN-poisoned `StridedMut`.
fn dj<S: Real, G: LieGroup<S>, const D: usize>(j: &G::Jac) -> [[f64; D]; D] {
    const { assert!(D == G::DOF && D == <G::Tangent as Tangent<S>>::DOF) };
    let mut buf = [[nan::<S>(); D]; D];
    j.write_dense(&mut StridedMut::col_major(buf.as_flattened_mut(), D, D));
    buf.map(|col| col.map(S::value_f64))
}

fn terr<S: Real, G: LieGroup<S>, const D: usize>(a: &G::Tangent, b: &G::Tangent) -> f64 {
    e::<S>(&dt::<S, G, D>(a), &dt::<S, G, D>(b))
}

fn jerr<S: Real, G: LieGroup<S>, const D: usize>(a: &G::Jac, b: &G::Jac) -> f64 {
    e::<S>(
        dj::<S, G, D>(a).as_flattened(),
        dj::<S, G, D>(b).as_flattened(),
    )
}

/// `‖Log(a ⊖_R b)‖` over `max(‖Log b‖, 1)`.
fn gerr<S: Real, G: LieGroup<S>, const D: usize>(a: &G, b: &G) -> f64 {
    let (d, r) = (dt::<S, G, D>(&a.rminus(b)), dt::<S, G, D>(&b.log()));
    norm(&d) / norm(&r).max(1.0) / unit::<S>()
}

/// `(x y) z = x (y z)`, `1 x = x 1 = x`, `x x⁻¹ = x⁻¹ x = 1`.
pub(crate) fn group_axioms<S: Real, G: LieGroup<S>, const D: usize>(x: G, y: G, z: G) -> f64 {
    let id = G::identity();
    let pairs = [
        ((x * y) * z, x * (y * z)),
        (id * x, x),
        (x * id, x),
        (x * x.inverse(), id),
        (x.inverse() * x, id),
    ];
    pairs
        .iter()
        .fold(0.0, |w, (a, b)| worst(w, gerr::<S, G, D>(a, b)))
}

/// `Log(Exp τ) = τ` for `θ < π`, and the backward error `Exp(Log x) ⊖ x` of `Log`.
pub(crate) fn exp_log_roundtrip<S: Real, G: LieGroup<S>, const D: usize>(
    tau: &G::Tangent,
    x: &G,
) -> f64 {
    worst(
        terr::<S, G, D>(&G::exp(tau).log(), tau),
        gerr::<S, G, D>(&G::exp(&x.log()), x),
    )
}

/// `x Exp(τ) x⁻¹ = Exp(Ad_x τ)`.
pub(crate) fn adjoint_identity<S: Real, G: LieGroup<S>, const D: usize>(
    x: &G,
    tau: &G::Tangent,
) -> f64 {
    gerr::<S, G, D>(
        &(*x * G::exp(tau) * x.inverse()),
        &G::exp(&x.adjoint().apply(tau)),
    )
}

/// `J_l = Ad_Exp(τ) J_r`, `J_r J_r⁻¹ = J_l J_l⁻¹ = I`.
pub(crate) fn jl_is_ad_jr<S: Real, G: LieGroup<S>, const D: usize>(tau: &G::Tangent) -> f64 {
    let (jr, jl) = (G::jr(tau), G::jl(tau));
    let id = <G::Jac as Jac<S, G::Tangent>>::identity();
    [
        (jl, G::exp(tau).adjoint().mul(&jr)),
        (jr.mul(&G::jr_inv(tau)), id),
        (jl.mul(&G::jl_inv(tau)), id),
    ]
    .iter()
    .fold(0.0, |w, (a, b)| worst(w, jerr::<S, G, D>(a, b)))
}

/// The worst error of two Jacobians of a pair.
fn perr<S: Real, G: LieGroup<S>, const D: usize>(
    a: &(G::Jac, G::Jac),
    b: &(G::Jac, G::Jac),
) -> f64 {
    worst(jerr::<S, G, D>(&a.0, &b.0), jerr::<S, G, D>(&a.1, &b.1))
}

/// `Side::*` is the `r*`/`l*` method of its side, for `plus`, `minus` and both `*_jacobians`.
pub(crate) fn side_delegation<S: Real, G: LieGroup<S>, const D: usize>(
    x: &G,
    y: &G,
    tau: &G::Tangent,
) -> f64 {
    [
        gerr::<S, G, D>(&Right::plus(x, tau), &x.rplus(tau)),
        gerr::<S, G, D>(&Left::plus(x, tau), &x.lplus(tau)),
        terr::<S, G, D>(&Right::minus(y, x), &y.rminus(x)),
        terr::<S, G, D>(&Left::minus(y, x), &y.lminus(x)),
        perr::<S, G, D>(&Right::plus_jacobians(x, tau), &x.rplus_jacobians(tau)),
        perr::<S, G, D>(&Left::plus_jacobians(x, tau), &x.lplus_jacobians(tau)),
        perr::<S, G, D>(&Right::minus_jacobians(y, x), &y.rminus_jacobians(x)),
        perr::<S, G, D>(&Left::minus_jacobians(y, x), &y.lminus_jacobians(x)),
    ]
    .into_iter()
    .fold(0.0, worst)
}

/// `(x ⊕_R τ) ⊖_R x = τ`, `(x ⊕_L τ) ⊖_L x = τ`, `y ⊖_L x = Ad_y (y ⊖_R x)` and
/// `x ⊕_L τ = x ⊕_R Ad_{x⁻¹} τ`: the values of the four perturbations, which `side_delegation`
/// only compares with each other.
///
/// # Domain
///
/// `θ < π` for the tangents involved.
pub(crate) fn plus_minus<S: Real, G: LieGroup<S>, const D: usize>(
    x: &G,
    y: &G,
    tau: &G::Tangent,
) -> f64 {
    [
        terr::<S, G, D>(&x.rplus(tau).rminus(x), tau),
        terr::<S, G, D>(&x.lplus(tau).lminus(x), tau),
        terr::<S, G, D>(&y.lminus(x), &y.adjoint().apply(&y.rminus(x))),
        gerr::<S, G, D>(&x.lplus(tau), &x.rplus(&x.inverse().adjoint().apply(tau))),
    ]
    .into_iter()
    .fold(0.0, worst)
}

/// The rows of `NUMERICS.md` §2.3, each against `Ad`, `J_r`, `J_l` and their inverses: `⊕`, `⊖`,
/// `X Y` and `X⁻¹`, both sides. The rows are not compared with `Dual` differentiation here
/// (`docs/PHASE3.md` §8).
pub(crate) fn jacobian_rows<S: Real, G: LieGroup<S>, const D: usize>(
    x: &G,
    y: &G,
    tau: &G::Tangent,
) -> f64 {
    let id = <G::Jac as Jac<S, G::Tangent>>::identity();
    let ad_x = x.adjoint();
    let (ad_exp, ad_exp_inv) = (G::exp(tau).adjoint(), G::exp(&tau.neg()).adjoint());
    let (tr, tl) = (y.rminus(x), y.lminus(x));
    let pairs = [
        (x.rplus_jacobians(tau), (ad_exp_inv, G::jr(tau))),
        (x.lplus_jacobians(tau), (ad_exp, G::jl(tau))),
        (
            y.rminus_jacobians(x),
            (G::jr_inv(&tr), G::jl_inv(&tr).neg()),
        ),
        (
            y.lminus_jacobians(x),
            (G::jl_inv(&tl), G::jr_inv(&tl).neg()),
        ),
        (x.compose_jacobians::<Right>(y), (y.inverse().adjoint(), id)),
        (x.compose_jacobians::<Left>(y), (id, ad_x)),
    ];
    let singles = [
        (x.inverse_jacobian::<Right>(), ad_x.neg()),
        (x.inverse_jacobian::<Left>(), x.inverse().adjoint().neg()),
    ];
    let w = pairs
        .iter()
        .fold(0.0, |w, (a, b)| worst(w, perr::<S, G, D>(a, b)));
    singles
        .iter()
        .fold(w, |w, (a, b)| worst(w, jerr::<S, G, D>(a, b)))
}

/// The terms of the series in [`ad_consistency`]: `‖ad‖^40 / 40!` is below `1e-20` for `‖ad‖ <= 4`.
const SERIES: usize = 40;

/// `ad` is a Lie bracket, `ad_a a = 0`, `ad_a b = -ad_b a` and `ad_{ad_a b} c = ad_a ad_b c - ad_b
/// ad_a c`, and it defines `Ad` and the Jacobians of `Exp` through their series:
/// `Ad_{Exp a} c = Σ ad_a^k c / k!`, `J_r(a) c = Σ (-ad_a)^k c / (k + 1)!` and
/// `J_l(a) c = Σ ad_a^k c / (k + 1)!`, cut at [`SERIES`] terms.
///
/// # Domain
///
/// `‖ad_a‖ <= 4`.
pub(crate) fn ad_consistency<S: Real, G: LieGroup<S>, const D: usize>(
    a: &G::Tangent,
    b: &G::Tangent,
    c: &G::Tangent,
) -> f64 {
    type T<S, G> = <G as LieGroup<S>>::Tangent;
    let (ada, adb) = (G::ad(a), G::ad(b));
    let ab = ada.apply(b);
    let zero = <T<S, G> as Tangent<S>>::zero();
    // `p = ad_a^k c` and `f = 1 / k!`.
    let (mut p, mut f) = (*c, S::one());
    let (mut ex, mut jr, mut jl) = (zero, zero, zero);
    for k in 0..SERIES {
        let f1 = f / S::lit((k + 1) as f64);
        let (term, next) = (p.scale(f), p.scale(f1));
        ex = ex.add(&term);
        jl = jl.add(&next);
        jr = if k % 2 == 0 {
            jr.add(&next)
        } else {
            jr.sub(&next)
        };
        (p, f) = (ada.apply(&p), f1);
    }
    [
        (ada.apply(a), zero),
        (ab.add(&adb.apply(a)), zero),
        (
            G::ad(&ab).apply(c),
            ada.apply(&adb.apply(c)).sub(&adb.apply(&ada.apply(c))),
        ),
        (G::exp(a).adjoint().apply(c), ex),
        (G::jr(a).apply(c), jr),
        (G::jl(a).apply(c), jl),
    ]
    .iter()
    .fold(0.0, |w, (x, y)| worst(w, terr::<S, G, D>(x, y)))
}

/// The tangent operations against the dense arrays `write_dense` produces.
pub(crate) fn tangent_dense_order<S: Real, G: LieGroup<S>, const D: usize>(
    a: &G::Tangent,
    b: &G::Tangent,
    k: S,
) -> f64 {
    type T<S, G> = <G as LieGroup<S>>::Tangent;
    let (da, db) = (dt::<S, G, D>(a), dt::<S, G, D>(b));
    let dot = (0..D).map(|i| da[i] * db[i]).sum::<f64>();
    let mut buf = [S::zero(); D];
    a.write_dense(&mut buf);
    [
        (dt::<S, G, D>(&a.add(b)), array::from_fn(|i| da[i] + db[i])),
        (dt::<S, G, D>(&a.sub(b)), array::from_fn(|i| da[i] - db[i])),
        (dt::<S, G, D>(&a.neg()), da.map(|x| -x)),
        (dt::<S, G, D>(&a.scale(k)), da.map(|x| x * k.value_f64())),
        (dt::<S, G, D>(&<T<S, G> as Tangent<S>>::zero()), [0.0; D]),
        (
            dt::<S, G, D>(&<T<S, G> as Tangent<S>>::read_dense(&buf)),
            da,
        ),
    ]
    .iter()
    .fold(e::<S>(&[a.dot(b).value_f64()], &[dot]), |w, (x, y)| {
        worst(w, e::<S>(x, y))
    })
}

/// The `Jac` operations against the dense matrices `write_dense` produces: `apply`,
/// `apply_transpose`, `mul`, `neg`, `identity` and `a.mul(a.inverse())`, so the dense order of `Jac`
/// agrees with that of `Tangent`.
///
/// # Domain
///
/// `a` is invertible.
pub(crate) fn jac_dense_order<S: Real, G: LieGroup<S>, const D: usize>(
    a: &G::Jac,
    b: &G::Jac,
    t: &G::Tangent,
) -> f64 {
    let (da, db, dv) = (dj::<S, G, D>(a), dj::<S, G, D>(b), dt::<S, G, D>(t));
    let apply: [f64; D] = array::from_fn(|r| (0..D).map(|c| da[c][r] * dv[c]).sum());
    let transpose: [f64; D] = array::from_fn(|r| (0..D).map(|c| da[r][c] * dv[c]).sum());
    let prod: [[f64; D]; D] =
        array::from_fn(|c| array::from_fn(|r| (0..D).map(|k| da[k][r] * db[c][k]).sum()));
    let id = <G::Jac as Jac<S, G::Tangent>>::identity();
    let eye = array::from_fn(|c| array::from_fn(|r| if r == c { 1.0 } else { 0.0 }));
    let vec = [
        (dt::<S, G, D>(&a.apply(t)), apply),
        (dt::<S, G, D>(&a.apply_transpose(t)), transpose),
    ];
    let mat = [
        (dj::<S, G, D>(&a.mul(b)), prod),
        (dj::<S, G, D>(&a.neg()), da.map(|c| c.map(|x| -x))),
        (dj::<S, G, D>(&id), eye),
        (dj::<S, G, D>(&a.mul(&a.inverse())), eye),
    ];
    let w = vec.iter().fold(0.0, |w, (x, y)| worst(w, e::<S>(x, y)));
    mat.iter().fold(w, |w, (x, y)| {
        worst(w, e::<S>(x.as_flattened(), y.as_flattened()))
    })
}

/// `j.sandwich(cov)` against `J Σ Jᵀ` formed densely in `f64` from `write_dense`.
pub(crate) fn sandwich_matches_dense<S: Real, G: LieGroup<S>, const D: usize>(
    j: &G::Jac,
    cov: &Matrix<S, D, D>,
) -> f64 {
    let c: [[f64; D]; D] = array::from_fn(|i| array::from_fn(|k| cov.get(i, k).value_f64()));
    let (d, got) = (dj::<S, G, D>(j), j.sandwich(cov));
    let want: [[f64; D]; D] = array::from_fn(|i| {
        array::from_fn(|k| {
            (0..D)
                .map(|p| (0..D).map(|q| d[p][i] * c[p][q] * d[q][k]).sum::<f64>())
                .sum()
        })
    });
    let got: [[f64; D]; D] = array::from_fn(|i| array::from_fn(|k| got.get(i, k).value_f64()));
    e::<S>(got.as_flattened(), want.as_flattened())
}

/// Every method of `G` at the tangents `a`, `b`, as value parts: run over `f64` and over
/// `Dual<f64, D>`, the two must agree to the bit (`dual_value_is_plain_value`).
pub(crate) fn probe<S: Real, G: LieGroup<S>, const D: usize>(
    a: &[S; D],
    b: &[S; D],
    k: S,
) -> Vec<f64> {
    type T<S, G> = <G as LieGroup<S>>::Tangent;
    let (ta, tb) = (
        <T<S, G> as Tangent<S>>::read_dense(a),
        <T<S, G> as Tangent<S>>::read_dense(b),
    );
    let (x, y) = (G::exp(&ta), G::exp(&tb));
    let (ja, jb) = (G::jr(&ta), G::adjoint(&y));
    let cov = Matrix::from_cols([helicoid_linalg::Vector(*a); D]);
    let mut out = Vec::new();
    let g = [
        G::identity(),
        x.inverse(),
        x * y,
        x.rplus(&tb),
        x.lplus(&tb),
    ];
    let t = [
        x.log(),
        x.rminus(&y),
        x.lminus(&y),
        ta.add(&tb),
        ta.sub(&tb),
        ta.neg(),
        ta.scale(k),
        <T<S, G> as Tangent<S>>::zero(),
        ja.apply(&tb),
        jb.apply_transpose(&tb),
    ];
    let mut j = Vec::from([x.adjoint(), G::ad(&ta), ja, G::jr_inv(&ta), G::jl(&ta)]);
    j.extend([G::jl_inv(&ta), ja.mul(&jb), ja.inverse(), jb.neg()]);
    j.push(<G::Jac as Jac<S, G::Tangent>>::identity());
    for (p, q) in [
        x.rplus_jacobians(&tb),
        x.lplus_jacobians(&tb),
        x.rminus_jacobians(&y),
        x.lminus_jacobians(&y),
        x.compose_jacobians::<Right>(&y),
        x.compose_jacobians::<Left>(&y),
    ] {
        j.extend([p, q]);
    }
    j.extend([x.inverse_jacobian::<Right>(), x.inverse_jacobian::<Left>()]);
    for v in g {
        out.extend(dt::<S, G, D>(&v.log()));
    }
    for v in &t {
        out.extend(dt::<S, G, D>(v));
    }
    out.push(ta.dot(&tb).value_f64());
    for v in &j {
        out.extend(dj::<S, G, D>(v).as_flattened());
    }
    let s = ja.sandwich::<D>(&cov);
    out.extend((0..D * D).map(|i| s.get(i % D, i / D).value_f64()));
    out
}

/// The bound of each law, in `u` (`NUMERICS.md` §11): twice the worst error of one run of
///
/// ```text
/// PROPTEST_CASES=1000000 PROPTEST_RNG_SEED=1 cargo nextest run --release -p helicoid <group>_tests
/// ```
///
/// rounded up to an integer, and `0` where that error is `0`. The worst errors sit in the header
/// of each group's tests.
pub(crate) struct Bounds {
    pub(crate) axioms: f64,
    pub(crate) exp_log: f64,
    pub(crate) adjoint: f64,
    pub(crate) jl_ad_jr: f64,
    pub(crate) plus_minus: f64,
    pub(crate) rows: f64,
    pub(crate) ad: f64,
    pub(crate) sides: f64,
    pub(crate) tangent_order: f64,
    pub(crate) jac_order: f64,
    pub(crate) sandwich: f64,
}

/// Tangent samples with entries `m 2^e`, `|m| < 1`, `-6 <= e <= 0`.
pub(crate) fn sample<const D: usize>() -> impl Strategy<Value = [f64; D]> {
    let entry = (-1.0_f64..1.0, -6_i32..=0).prop_map(|(m, e)| m * 2_f64.powi(e));
    proptest::array::uniform(entry)
}

/// `Ok` when the worst error `v` of a law is within `bound`; a NaN is not.
pub(crate) fn within(v: f64, bound: f64) -> Result<(), TestCaseError> {
    if v <= bound {
        Ok(())
    } else {
        Err(TestCaseError::fail(format!(
            "worst error {v} u exceeds {bound} u"
        )))
    }
}

/// Every law above as a proptest in a module `$p`, for `f64`, `f32` and `Dual<f64, D>` (which takes
/// the bounds of `f64`, its values being the same), for the `D`-dimensional group `$G<S>` (`D = 3`
/// when omitted) whose tests define `$jac::<S>(&[f64; D])`, an invertible `Jac` built from a sample.
macro_rules! laws_for {
    ($p:ident, $G:ident, $jac:ident, $b64:ident, $b32:ident) => {
        $crate::laws::laws_for!($p, $G, $jac, $b64, $b32, 3);
    };
    ($p:ident, $G:ident, $jac:ident, $b64:ident, $b32:ident, $D:literal) => {
        mod $p {
            use super::*;
            $crate::laws::laws_for!(@case as_f64, f64, $b64, $G, $jac, $D);
            $crate::laws::laws_for!(@case as_f32, f32, $b32, $G, $jac, $D);
            $crate::laws::laws_for!(@case as_dual, helicoid_linalg::Dual<f64, $D>, $b64, $G, $jac, $D);
            $crate::laws::laws_for!(@plain $G, $D);
        }
    };
    (@case $m:ident, $S:ty, $B:ident, $G:ident, $jac:ident, $D:literal) => {
        mod $m {
            use super::*;
            use $crate::laws::{self, sample, tangent, within, Sample};
            use $crate::LieGroup;
            use core::array;
            use helicoid_linalg::{Matrix, Vector};
            use proptest::prelude::*;

            type S = $S;
            type Gp = $G<S>;

            fn t(v: &[f64; $D]) -> <Gp as LieGroup<S>>::Tangent {
                tangent::<S, Gp, $D>(v)
            }
            fn g(v: &[f64; $D]) -> Gp {
                Gp::exp(&t(v))
            }
            fn j(v: &[f64; $D]) -> <Gp as LieGroup<S>>::Jac {
                $jac::<S>(v)
            }
            fn cov(a: &[f64; $D], b: &[f64; $D], c: &[f64; $D]) -> Matrix<S, $D, $D> {
                // Column `k` is sample `k % 3` rotated by `k / 3`: distinct columns, so no
                // column index of the sandwich aliases another for any `D`.
                let cols = [a, b, c];
                Matrix::from_cols(array::from_fn(|k| {
                    Vector(array::from_fn(|r| S::sample(cols[k % 3][(r + k / 3) % $D], r)))
                }))
            }

            proptest! {
                #[test]
                fn group_axioms(a in sample::<$D>(), b in sample::<$D>(), c in sample::<$D>()) {
                    within(laws::group_axioms::<S, Gp, $D>(g(&a), g(&b), g(&c)), $B.axioms)?;
                }
                #[test]
                fn exp_log_roundtrip(a in sample::<$D>(), b in sample::<$D>()) {
                    within(laws::exp_log_roundtrip::<S, Gp, $D>(&t(&a), &g(&b)), $B.exp_log)?;
                }
                #[test]
                fn adjoint_identity(a in sample::<$D>(), b in sample::<$D>()) {
                    within(laws::adjoint_identity::<S, Gp, $D>(&g(&a), &t(&b)), $B.adjoint)?;
                }
                #[test]
                fn jl_is_ad_jr(a in sample::<$D>()) {
                    within(laws::jl_is_ad_jr::<S, Gp, $D>(&t(&a)), $B.jl_ad_jr)?;
                }
                #[test]
                fn plus_minus(a in sample::<$D>(), b in sample::<$D>(), c in sample::<$D>()) {
                    within(laws::plus_minus::<S, Gp, $D>(&g(&a), &g(&b), &t(&c)), $B.plus_minus)?;
                }
                #[test]
                fn jacobian_rows(a in sample::<$D>(), b in sample::<$D>(), c in sample::<$D>()) {
                    within(laws::jacobian_rows::<S, Gp, $D>(&g(&a), &g(&b), &t(&c)), $B.rows)?;
                }
                #[test]
                fn ad_consistency(a in sample::<$D>(), b in sample::<$D>(), c in sample::<$D>()) {
                    within(laws::ad_consistency::<S, Gp, $D>(&t(&a), &t(&b), &t(&c)), $B.ad)?;
                }
                #[test]
                fn side_delegation(a in sample::<$D>(), b in sample::<$D>(), c in sample::<$D>()) {
                    within(laws::side_delegation::<S, Gp, $D>(&g(&a), &g(&b), &t(&c)), $B.sides)?;
                }
                #[test]
                fn tangent_dense_order(a in sample::<$D>(), b in sample::<$D>(), c in sample::<$D>()) {
                    let k = S::sample(c[0], 0);
                    within(laws::tangent_dense_order::<S, Gp, $D>(&t(&a), &t(&b), k), $B.tangent_order)?;
                }
                #[test]
                fn jac_dense_order(a in sample::<$D>(), b in sample::<$D>(), c in sample::<$D>()) {
                    within(laws::jac_dense_order::<S, Gp, $D>(&j(&a), &j(&b), &t(&c)), $B.jac_order)?;
                }
                #[test]
                fn sandwich_matches_dense(a in sample::<$D>(), b in sample::<$D>(), c in sample::<$D>()) {
                    within(laws::sandwich_matches_dense::<S, Gp, $D>(&j(&a), &cov(&a, &b, &c)), $B.sandwich)?;
                }
            }
        }
    };
    // `probe` reads every method of the group; over any `f64` and with NaN, infinite and huge
    // derivative lanes the value parts must be the plain result, bit for bit (NaN sign and payload
    // of arithmetic excepted, as for `Dual`).
    (@plain $G:ident, $D:literal) => {
        proptest::proptest! {
            #[test]
            fn dual_value_is_plain_value(
                a in proptest::array::uniform::<_, $D>(proptest::num::f64::ANY),
                b in proptest::array::uniform::<_, $D>(proptest::num::f64::ANY),
                k in proptest::num::f64::ANY,
            ) {
                type D3 = helicoid_linalg::Dual<f64, $D>;
                let plain = $crate::laws::probe::<f64, $G<f64>, $D>(&a, &b, k);
                for poison in [false, true] {
                    let lift = |v: &[f64; $D]| -> [D3; $D] {
                        core::array::from_fn(|i| if poison {
                            helicoid_linalg::Dual {
                                v: v[i],
                                d: core::array::from_fn(|l| [f64::NAN, f64::INFINITY, -1e300][l % 3]),
                            }
                        } else {
                            helicoid_linalg::Dual::variable(v[i], i)
                        })
                    };
                    let dual = $crate::laws::probe::<D3, $G<D3>, $D>(
                        &lift(&a), &lift(&b), D3::constant(k));
                    proptest::prop_assert_eq!(plain.len(), dual.len());
                    for (p, d) in plain.iter().zip(&dual) {
                        proptest::prop_assert!(p.to_bits() == d.to_bits() || (p.is_nan() && d.is_nan()));
                    }
                }
            }
        }
    };
}
pub(crate) use laws_for;
