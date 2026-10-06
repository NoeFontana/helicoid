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
use helicoid_linalg::{Dual, Matrix, Precision, Real};
use proptest::prelude::*;
use std::{format, vec::Vec};

/// A scalar a test can build from a sample: `Dual` gets one variable per component.
pub(crate) trait Sample: Real {
    fn sample(v: f64, i: usize) -> Self;
    /// The same value with no derivative, for a quantity a law holds fixed. `Real::lit` will not
    /// do: it debug-asserts that its argument is exactly representable, which a sample is not.
    fn constant(v: f64) -> Self;
}

impl Sample for f64 {
    fn sample(v: f64, _: usize) -> Self {
        v
    }
    fn constant(v: f64) -> Self {
        v
    }
}

impl Sample for f32 {
    fn sample(v: f64, _: usize) -> Self {
        v as f32
    }
    fn constant(v: f64) -> Self {
        v as f32
    }
}

impl<const N: usize> Sample for Dual<f64, N> {
    fn sample(v: f64, i: usize) -> Self {
        Dual::variable(v, i)
    }
    fn constant(v: f64) -> Self {
        Dual::constant(v)
    }
}

/// The tangent whose dense components are the samples `v`.
pub(crate) fn tangent<S: Sample, G: LieGroup<S>, const D: usize>(v: &[f64; D]) -> G::Tangent {
    let s: [S; D] = array::from_fn(|i| S::sample(v[i], i));
    <G::Tangent as Tangent<S>>::read_dense(&s)
}

pub(crate) fn unit<S: Real>() -> f64 {
    match S::PRECISION {
        Precision::F64 => f64::EPSILON / 2.0,
        Precision::F32 => f64::from(f32::EPSILON) / 2.0,
    }
}

/// The larger error; NaN wins.
pub(crate) fn worst(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

pub(crate) fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// `max(‖v‖, 1)`, the denominator of a norm-wise error. `f64::max` returns the *other* operand for
/// a NaN, so it would turn a NaN reference into a denominator of `1`; here NaN wins, as in
/// [`worst`]. Every error goes through this, so no law can pass against a NaN reference.
fn den(v: &[f64]) -> f64 {
    let n = norm(v);
    if n.is_nan() {
        f64::NAN
    } else {
        n.max(1.0)
    }
}

/// `‖a - b‖ / max(‖b‖, 1)` in units of `u`.
pub(crate) fn e<S: Real>(a: &[f64], b: &[f64]) -> f64 {
    let diff = a.iter().zip(b).map(|(x, y)| (x - y) * (x - y));
    diff.sum::<f64>().sqrt() / den(b) / unit::<S>()
}

/// `‖a - b‖ / max(scale, 1)` in units of `u`: [`e`] against a *stated* scale instead of the
/// reference's own norm.
///
/// `e` divides by `max(‖b‖, 1)`, which is the right question for a routine whose answer is the
/// size of its inputs and the wrong one for a routine whose answer can cancel to nothing. Both
/// spellings of a relative transform carry an absolute error of order `u max‖x‖`, so scoring their
/// difference against the *difference* reports the cancellation and not the routine: on one set of
/// draws `6.80 u` at `‖x‖ = 1` becomes `106.32 u` at `10^3` through [`e`] and stays `7.94 u`
/// through this. A bound is a property of the routine only if its denominator is one
/// (`0048` decision 4).
pub(crate) fn e_at<S: Real>(a: &[f64], b: &[f64], scale: f64) -> f64 {
    let diff = a.iter().zip(b).map(|(x, y)| (x - y) * (x - y));
    let den = if scale.is_nan() {
        f64::NAN
    } else {
        scale.max(1.0)
    };
    diff.sum::<f64>().sqrt() / den / unit::<S>()
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

/// The dense matrix of `j`, column-major (`m[c][r]`), through `reference::dense_oriented` so that
/// the NaN poison a hand case relies on has one definition.
///
/// Generic over the `Jac` and not over a group, so a hand case can read a block of a composite
/// Jacobian — a factor's — which belongs to no group of the composite's dimension.
fn dense<S: Real, T: Tangent<S>, J: Jac<S, T>, const D: usize>(j: &J) -> [[S; D]; D] {
    crate::reference::dense_oriented(j, 1, D)
}

/// [`dense`] as bits, so a hand case tells `+0` from `-0` and sees a NaN.
pub(crate) fn dense_bits<T: Tangent<f64>, J: Jac<f64, T>, const D: usize>(j: &J) -> [[u64; D]; D] {
    dense::<f64, T, J, D>(j).map(|c| c.map(f64::to_bits))
}

/// [`dense`] of a group's Jacobian as `f64`, with the dimension tie asserted.
fn dj<S: Real, G: LieGroup<S>, const D: usize>(j: &G::Jac) -> [[f64; D]; D] {
    const { assert!(D == G::DOF && D == <G::Tangent as Tangent<S>>::DOF) };
    dense::<S, G::Tangent, G::Jac, D>(j).map(|col| col.map(S::value_f64))
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
    // Unlike `e`, the numerator and the denominator are independent quantities: a NaN in `Log b`
    // does not reach `Log(b⁻¹ a)`, so only `den` keeps it from being scaled away.
    let (d, r) = (dt::<S, G, D>(&a.rminus(b)), dt::<S, G, D>(&b.log()));
    norm(&d) / den(&r) / unit::<S>()
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

/// `PHASE4.md` §1 and §3's names for [`geodesic_legs`]'s seven, in the order it returns them.
///
/// One bound per leg, not one over their maximum: `0006` asks for the bar on the max *per
/// quantity*, and folding seven identities into one number lets the tight ones rot -- at `SEn3`
/// the `symmetry` leg sets 71 `u` while `velocity` is exactly `0`, so a broken endpoint could
/// grow sixty-fold under a single folded bound and pass.
pub(crate) const GEODESIC_LEGS: [&str; 7] = [
    "t=0", "t=1", "symmetry", "velocity", "left", "right", "twin",
];

/// The reciprocal of the parameter step the `velocity` leg takes: the step is `1 / 3`, formed as a
/// division because `Real::lit` takes **exactly representable** constants and `1.0 / 3.0` is not
/// one at either precision (it fired on `f32` in the dev profile, which is what that assert is
/// for). `step * theta(d) < pi` for every drawn tangent, so the identity's own `Log` takes no
/// branch, and a non-dyadic step keeps a rounding in `scale` that `0.25` would remove.
const VELOCITY_STEP_RECIP: f64 = 3.0;

/// `PHASE4.md` §1 and §3's identities on one fixture, each leg measured and bounded apart: the
/// two endpoints, GE.2(c)'s symmetry, the velocity, GE.4's left and right invariance, and D6's tie
/// to `reference::geodesic`.
///
/// [`GEODESIC_LEGS`] names them and the return type is its length, so a leg added without a name
/// fails the build instead of being dropped silently by a `zip`.
///
/// `x1` is `x0 ⊕_R d` and not a second draw, so `θ(d)` is the drawn tangent's and the curve stays
/// inside [`LieGroup::geodesic`]'s domain: two independent draws of `sample`'s distribution
/// compose to a relative rotation of up to `2√3 > π`, where the answer is the quaternion sign's
/// and the law would be measuring GE.13(d)'s conditioning instead of the identity. It is also how
/// `geo:consecutive` builds its records, for the same reason.
///
/// Right invariance holds for **every** group, not only SE(3): `γ(x0 h, x1 h, t) = γ(x0, x1, t) h`
/// follows from `Log(h⁻¹ Δ h) = Ad_{h⁻¹} d` and `h Exp(Ad_{h⁻¹} ξ) = Exp(ξ) h`, with no
/// commutativity used. `Product<SO3, R3>`'s famous failure is against the **SE(3) reading** of
/// `(R, t)`, a different composition from `Product`'s own `Mul`, which
/// `product_tests::so3_r3::the_se3_reading_fails_right_invariance_by_ge5bs_closed_form` asserts
/// on the group §1.3 names, against GE.5(b)'s closed form rather than a threshold (`0045` item 4).
///
/// `h` is drawn at the fixture's own scale here; `0045` item 5 asks for the left-invariance bound
/// at `‖t_G‖ ∈ {0, 1, 1e4}`, which `product_tests::so3_r3::left_invariance_at_three_translation_scales`
/// measures, since the scales are an `SEn3` quantity and this law is generic.
///
/// # Domain
///
/// `θ(d) < π`, `θ(Ad_{h⁻¹} d) < π` (the same angle on every group that ships), and
/// `θ((t + 1/3) d) < π` for the velocity leg.
pub(crate) fn geodesic_legs<S: Real, G: LieGroup<S>, const D: usize>(
    x0: &G,
    d: &G::Tangent,
    h: &G,
    t: S,
) -> [f64; GEODESIC_LEGS.len()] {
    let x1 = x0.rplus(d);
    let at = |a: &G, b: &G, s: S| G::geodesic(a, b, s);
    let one = S::one();
    // One evaluation of the curve at `t` for the four legs that compare against it: recomputing
    // it was two thirds of the measurement's work.
    let here = at(x0, &x1, t);
    let step = S::one() / S::lit(VELOCITY_STEP_RECIP);
    [
        gerr::<S, G, D>(&at(x0, &x1, S::zero()), x0),
        gerr::<S, G, D>(&at(x0, &x1, one), &x1),
        gerr::<S, G, D>(&here, &at(&x1, x0, one - t)),
        // GE.2(b) in finite form: the body displacement over a fixed parameter interval is
        // `step * velocity` wherever the interval starts. Comparing `geodesic_velocity` against
        // `rminus` instead would be a tautology -- `rminus` *is* its body, and `Product`'s
        // override is `rminus` per factor -- so it would read exactly `0` however wrong the
        // returned tangent was.
        terr::<S, G, D>(
            &at(x0, &x1, t + step).rminus(&here),
            &G::geodesic_velocity(x0, &x1).scale(step),
        ),
        gerr::<S, G, D>(&at(&(*h * *x0), &(*h * x1), t), &(*h * here)),
        gerr::<S, G, D>(&at(&(*x0 * *h), &(x1 * *h), t), &(here * *h)),
        // D6's tie: a group that overrides `geodesic` is measured against the expression the
        // provided body is. Exactly `0` for a group that does not override, the two being the
        // same call -- which is what makes it free to carry everywhere.
        gerr::<S, G, D>(&here, &crate::reference::geodesic(x0, &x1, t)),
    ]
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

/// Every row of `NUMERICS.md` §2.3 against `Dual<S, D>` differentiation of the operation it is the
/// Jacobian *of* — the second of the two checks `PHASE3.md` §8 asks for, where
/// [`jacobian_rows`] is the first.
///
/// The difference matters: `jacobian_rows` compares each closed form with a chain of the crate's
/// own primitives, so a sign carried consistently through `jr`, `jl` and `Ad` would satisfy it.
/// This differentiates the *operation* — `X ⊕ τ`, `Y ⊖ X`, `X Y`, `X⁻¹` — and so answers "is this
/// matrix the derivative" from outside the closed forms entirely.
///
/// Everything runs at one scalar, `Dual<S, D>`: the elements are lifted as constants, the
/// perturbation `δ` is seeded with [`Dual::variable`], and the closed forms are evaluated at the
/// same scalar, where their *value* lanes are the answer `S` would have given. A row whose output
/// is a group element is differentiated through `Log` of the correction, which is what the
/// perturbation *means* (§2.3 states each row in its own side's convention):
///
/// ```text
/// right: J[i] = d/dδ_i Log( f(X)⁻¹ · f(X ⊕_R δ) )   left: d/dδ_i Log( f(X ⊕_L δ) · f(X)⁻¹ )
/// ```
///
/// A row whose output is already a tangent (`⊖`) is differentiated directly.
///
/// # Domain
///
/// `D == G::DOF`. The error is `laws::e`'s, so a bound is in `u`; it is **not** `0` for any group,
/// because the two sides compute different expressions — the closed form is one formula and this is
/// a difference of two `Log`s divided by nothing, so the comparison carries `Log`'s conditioning
/// near a half turn.
pub(crate) fn jacobians_match_dual<S: Real, G: LieGroup<Dual<S, D>>, const D: usize>(
    a: &[f64; D],
    b: &[f64; D],
    c: &[f64; D],
) -> f64 {
    type T<S, G> = <G as LieGroup<S>>::Tangent;
    let lift = |v: &[f64; D]| -> T<Dual<S, D>, G> {
        <T<Dual<S, D>, G> as Tangent<Dual<S, D>>>::read_dense(&v.map(|x| Dual::constant(S::lit(x))))
    };
    // `δ`, seeded: lane `i` is `∂/∂δ_i`, and the value is `0`, so every element below is its own
    // unperturbed self in the value lane.
    let seeds: [Dual<S, D>; D] = array::from_fn(|i| Dual::variable(S::zero(), i));
    let delta = <T<Dual<S, D>, G> as Tangent<Dual<S, D>>>::read_dense(&seeds);
    let (x, y, tau) = (G::exp(&lift(a)), G::exp(&lift(b)), lift(c));
    // The `D x D` in a tangent's derivative lanes, column-major as `dj` reads a `Jac`.
    let lanes = |t: &T<Dual<S, D>, G>| -> [[f64; D]; D] {
        let mut buf = [Dual::constant(S::zero()); D];
        t.write_dense(&mut buf);
        array::from_fn(|col| array::from_fn(|row| buf[row].d[col].value_f64()))
    };
    // `Log` of the correction, in the side's own convention.
    let right_of = |moved: &G, base: &G| lanes(&moved.rminus(base));
    let left_of = |moved: &G, base: &G| lanes(&moved.lminus(base));
    let mut worst = 0.0_f64;
    let mut check = |closed: &G::Jac, got: [[f64; D]; D]| {
        let want = dj::<Dual<S, D>, G, D>(closed);
        worst = crate::laws::worst(
            worst,
            e::<Dual<S, D>>(want.as_flattened(), got.as_flattened()),
        );
    };

    // `X ⊕ τ`: the first row is `∂/∂X`, the second `∂/∂τ`.
    let (rp, lp) = (x.rplus(&tau), x.lplus(&tau));
    let (jr_x, jr_t) = x.rplus_jacobians(&tau);
    check(&jr_x, right_of(&x.rplus(&delta).rplus(&tau), &rp));
    check(&jr_t, right_of(&x.rplus(&tau.add(&delta)), &rp));
    let (jl_x, jl_t) = x.lplus_jacobians(&tau);
    check(&jl_x, left_of(&x.lplus(&delta).lplus(&tau), &lp));
    check(&jl_t, left_of(&x.lplus(&tau.add(&delta)), &lp));

    // `Y ⊖ X`, whose output is a tangent already.
    let (rm_y, rm_x) = y.rminus_jacobians(&x);
    check(&rm_y, lanes(&y.rplus(&delta).rminus(&x)));
    check(&rm_x, lanes(&y.rminus(&x.rplus(&delta))));
    let (lm_y, lm_x) = y.lminus_jacobians(&x);
    check(&lm_y, lanes(&y.lplus(&delta).lminus(&x)));
    check(&lm_x, lanes(&y.lminus(&x.lplus(&delta))));

    // `X Y` and `X⁻¹`, both sides.
    let xy = x * y;
    let (cr_x, cr_y) = x.compose_jacobians::<Right>(&y);
    check(&cr_x, right_of(&(x.rplus(&delta) * y), &xy));
    check(&cr_y, right_of(&(x * y.rplus(&delta)), &xy));
    let (cl_x, cl_y) = x.compose_jacobians::<Left>(&y);
    check(&cl_x, left_of(&(x.lplus(&delta) * y), &xy));
    check(&cl_y, left_of(&(x * y.lplus(&delta)), &xy));
    let inv = x.inverse();
    check(
        &x.inverse_jacobian::<Right>(),
        right_of(&x.rplus(&delta).inverse(), &inv),
    );
    check(
        &x.inverse_jacobian::<Left>(),
        left_of(&x.lplus(&delta).inverse(), &inv),
    );

    // `Exp(τ)` and `Log(X)`: `J_r(τ)` and `J_r⁻¹(Log X)` are those rows (§2.3), and the left pair
    // with them.
    check(
        &G::jr(&tau),
        right_of(&G::exp(&tau.add(&delta)), &G::exp(&tau)),
    );
    check(
        &G::jl(&tau),
        left_of(&G::exp(&tau.add(&delta)), &G::exp(&tau)),
    );
    check(&G::jr_inv(&x.log()), lanes(&x.rplus(&delta).log()));
    check(&G::jl_inv(&x.log()), lanes(&x.lplus(&delta).log()));
    worst
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
    /// [`jacobians_match_dual`]'s bound, `PHASE3.md` §8's second check.
    pub(crate) dual_rows: f64,
    /// One bound per leg of [`geodesic_legs`], in [`GEODESIC_LEGS`]'s order (`PHASE4.md` §1
    /// and §3).
    ///
    /// Two legs read a fixed number on every group, and the comments that record a group's own
    /// figures do not repeat why. **`twin`** is exactly `0` unless the group overrides
    /// [`LieGroup::geodesic`]: the provided body and `reference::geodesic` are then the same call.
    /// **`t=0`** is [`gerr`]'s floor for two *bitwise equal* elements -- about `1.118` on a
    /// quaternion group, where `q* q` leaves one `u` in the vector part, and `0` where the inverse
    /// is exact -- so it is not a geodesic error. The statement it cannot make, that the answer
    /// *is* the left endpoint bit for bit, is made by
    /// `geodesic_at_zero_is_the_left_endpoint_bit_for_bit` in the groups whose representation a
    /// test can read.
    pub(crate) geodesic: [f64; GEODESIC_LEGS.len()],
}

/// splitmix64, seeded; `unif` is uniform on `[-1, 1)`.
///
/// The stream every measurement draws from, so no method may change: the figures the bounds are
/// recorded from would stop being reproducible. A measurement that gains a law gives it a stream
/// of its own rather than drawing from one the others already use.
pub(crate) struct Rng(pub(crate) u64);

impl Rng {
    pub(crate) fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    pub(crate) fn unif(&mut self) -> f64 {
        (self.next() >> 11) as f64 / 2_f64.powi(52) - 1.0
    }
    /// `N` draws of [`sample`]'s distribution -- `m 2^e`, `|m| < 1`, `-6 <= e <= 0` -- which is
    /// **not** [`arr`](Self::arr)'s uniform `[-1, 1)`.
    ///
    /// A bound is a measurement of what the bar samples (`0006`), so a measurement feeding a
    /// proptest's bound has to draw the proptest's distribution: `sample`'s `2^e` factor reaches
    /// `theta = 0.02` in one draw of seven per component, where uniform `[-1, 1)^3` reaches it
    /// about once in `10^6`. Two draws per entry, the mantissa then the exponent, in that order.
    pub(crate) fn shaped<const N: usize>(&mut self) -> [f64; N] {
        let mut out = [0.0; N];
        for x in &mut out {
            let m = self.unif();
            *x = m * 2_f64.powi(-((self.next() % 7) as i32));
        }
        out
    }

    /// `N` draws into a caller-owned array, the stream [`unif`](Self::unif) gives in order:
    /// uniform `[-1, 1)`, which is **not** [`shaped`](Self::shaped)'s `m 2^e`.
    ///
    /// The loop is explicit rather than `array::from_fn` because the order of the draws *is* the
    /// stream: a helper whose visiting order is unspecified would put the reproducibility of every
    /// recorded figure at the mercy of its implementation. No allocation, so a `10^6`-case
    /// measurement does not pay one per draw-set per iteration.
    pub(crate) fn arr<const N: usize>(&mut self) -> [f64; N] {
        let mut out = [0.0; N];
        for x in &mut out {
            *x = self.unif();
        }
        out
    }
}

/// The covariance `sandwich_matches_dense` and a `sandwich` twin are measured on: column `k` is
/// sample `k % 3` rotated down by `k / 3`.
///
/// `Σ` is not symmetric, so `sandwich` cannot pass with `J Σ Jᵀ` transposed. Cycling the three
/// samples alone repeats a column every three: at `D = 5` columns 0 and 3 would be equal, and
/// `sandwich` could read one where the other was meant. At `D = 3` the rotation is the identity,
/// so the recorded bounds of the three-dimensional groups stand.
pub(crate) fn cov<S: Sample, const D: usize>(
    a: &[f64; D],
    b: &[f64; D],
    c: &[f64; D],
) -> Matrix<S, D, D> {
    let cols = [a, b, c];
    Matrix::from_cols(array::from_fn(|k| {
        helicoid_linalg::Vector(array::from_fn(|r| {
            S::sample(cols[k % 3][(r + k / 3) % D], r)
        }))
    }))
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

/// [`within`], naming the leg it came from: a folded bound says only that *something* moved.
pub(crate) fn within_leg(name: &str, v: f64, bound: f64) -> Result<(), TestCaseError> {
    within(v, bound).map_err(|e| TestCaseError::fail(format!("leg `{name}`: {e}")))
}

/// Every law above as a proptest in a module `$p`, for `f64`, `f32` and `Dual<f64, $D>` (which takes
/// the bounds of `f64`, its values being the same), for the `$D`-dimensional group `$G<S>` whose
/// tests define `$jac::<S>(&[f64; $D])`, an invertible `Jac` built from a sample.
///
/// `$D` is the group's `DOF`. The laws themselves are generic over it (`D` is a parameter of each);
/// only the proptest bodies need it spelled, since a `Strategy` has to name its array length.
macro_rules! laws_for {
    ($p:ident, $G:ident, $jac:ident, $D:literal, $b64:ident, $b32:ident) => {
        mod $p {
            use super::*;
            $crate::laws::laws_for!(@case as_f64, f64, $b64, $G, $jac, $D);
            $crate::laws::laws_for!(@case as_f32, f32, $b32, $G, $jac, $D);
            $crate::laws::laws_for!(
                @case as_dual, helicoid_linalg::Dual<f64, $D>, $b64, $G, $jac, $D);
            $crate::laws::laws_for!(@plain $G, $D, $b64);
        }
    };
    (@case $m:ident, $S:ty, $B:ident, $G:ident, $jac:ident, $D:literal) => {
        mod $m {
            use super::*;
            use $crate::laws::{self, cov, sample, tangent, within, Sample};
            use $crate::LieGroup;
            use proptest::prelude::*;

            type S = $S;
            type Gp = $G<S>;
            const D: usize = $D;

            fn t(v: &[f64; D]) -> <Gp as LieGroup<S>>::Tangent {
                tangent::<S, Gp, D>(v)
            }
            fn g(v: &[f64; D]) -> Gp {
                Gp::exp(&t(v))
            }
            fn j(v: &[f64; D]) -> <Gp as LieGroup<S>>::Jac {
                $jac::<S>(v)
            }

            proptest! {
                #[test]
                fn group_axioms(a in sample::<D>(), b in sample::<D>(), c in sample::<D>()) {
                    within(laws::group_axioms::<S, Gp, D>(g(&a), g(&b), g(&c)), $B.axioms)?;
                }
                #[test]
                fn exp_log_roundtrip(a in sample::<D>(), b in sample::<D>()) {
                    within(laws::exp_log_roundtrip::<S, Gp, D>(&t(&a), &g(&b)), $B.exp_log)?;
                }
                #[test]
                fn adjoint_identity(a in sample::<D>(), b in sample::<D>()) {
                    within(laws::adjoint_identity::<S, Gp, D>(&g(&a), &t(&b)), $B.adjoint)?;
                }
                #[test]
                fn jl_is_ad_jr(a in sample::<D>()) {
                    within(laws::jl_is_ad_jr::<S, Gp, D>(&t(&a)), $B.jl_ad_jr)?;
                }
                #[test]
                fn plus_minus(a in sample::<D>(), b in sample::<D>(), c in sample::<D>()) {
                    within(laws::plus_minus::<S, Gp, D>(&g(&a), &g(&b), &t(&c)), $B.plus_minus)?;
                }
                #[test]
                fn jacobian_rows(a in sample::<D>(), b in sample::<D>(), c in sample::<D>()) {
                    within(laws::jacobian_rows::<S, Gp, D>(&g(&a), &g(&b), &t(&c)), $B.rows)?;
                }
                #[test]
                fn ad_consistency(a in sample::<D>(), b in sample::<D>(), c in sample::<D>()) {
                    within(laws::ad_consistency::<S, Gp, D>(&t(&a), &t(&b), &t(&c)), $B.ad)?;
                }
                #[test]
                fn side_delegation(a in sample::<D>(), b in sample::<D>(), c in sample::<D>()) {
                    within(laws::side_delegation::<S, Gp, D>(&g(&a), &g(&b), &t(&c)), $B.sides)?;
                }
                #[test]
                fn tangent_dense_order(a in sample::<D>(), b in sample::<D>(), c in sample::<D>()) {
                    let k = S::sample(c[0], 0);
                    within(laws::tangent_dense_order::<S, Gp, D>(&t(&a), &t(&b), k), $B.tangent_order)?;
                }
                #[test]
                fn jac_dense_order(a in sample::<D>(), b in sample::<D>(), c in sample::<D>()) {
                    within(laws::jac_dense_order::<S, Gp, D>(&j(&a), &j(&b), &t(&c)), $B.jac_order)?;
                }
                #[test]
                fn sandwich_matches_dense(a in sample::<D>(), b in sample::<D>(), c in sample::<D>()) {
                    let s = laws::sandwich_matches_dense::<S, Gp, D>(&j(&a), &cov::<S, D>(&a, &b, &c));
                    within(s, $B.sandwich)?;
                }
                #[test]
                fn geodesic(a in sample::<D>(), b in sample::<D>(), c in sample::<D>()) {
                    let at = S::sample(c[0], 0);
                    let legs = laws::geodesic_legs::<S, Gp, D>(&g(&a), &t(&b), &g(&c), at);
                    for ((name, v), bound) in
                        laws::GEODESIC_LEGS.iter().zip(legs).zip($B.geodesic)
                    {
                        laws::within_leg(name, v, bound)?;
                    }
                }
            }

            /// The worst error of the geodesic laws over 10^6 seeded cases, the figures
            /// `Bounds::geodesic` is recorded from. Its own stream, so adding it left every other
            /// law's recorded figure reproducible (`laws::Rng`).
            #[test]
            #[ignore = "measurement: prints the figures the bounds are recorded from"]
            #[allow(clippy::print_stdout)]
            fn measure_geodesic() {
                let mut rng = $crate::laws::Rng(0x6765_6F64_6573_6963);
                let mut w = [0.0_f64; laws::GEODESIC_LEGS.len()];
                for _ in 0..1_000_000 {
                    // `shaped`, not `arr`, and `t` from the third draw as the proptest takes it:
                    // the measurement has to sample what the bound bounds.
                    let (a, b, c) = (rng.shaped::<D>(), rng.shaped::<D>(), rng.shaped::<D>());
                    let s = S::sample(c[0], 0);
                    let legs = laws::geodesic_legs::<S, Gp, D>(&g(&a), &t(&b), &g(&c), s);
                    for (acc, v) in w.iter_mut().zip(legs) {
                        *acc = laws::worst(*acc, v);
                    }
                }
                std::print!("{} geodesic", module_path!());
                for (n, v) in laws::GEODESIC_LEGS.iter().zip(w) {
                    std::print!("  {n} {v:.3}");
                }
                std::println!();
            }
        }
    };
    // `probe` reads every method of the group; over any `f64` and with NaN, infinite and huge
    // derivative lanes the value parts must be the plain result, bit for bit (NaN sign and payload
    // of arithmetic excepted, as for `Dual`).
    //
    // `jacobians_match_dual` sits here and not in `@case` because it names one scalar of its own:
    // it differentiates at `Dual<f64, D>` and reads the closed forms' value lanes at the same
    // scalar, so running it per precision would measure nothing extra (`PHASE3.md` §8).
    (@plain $G:ident, $D:literal, $B:ident) => {
        proptest::proptest! {
            #[test]
            fn jacobians_match_dual(
                a in $crate::laws::sample::<$D>(),
                b in $crate::laws::sample::<$D>(),
                c in $crate::laws::sample::<$D>(),
            ) {
                $crate::laws::within(
                    $crate::laws::jacobians_match_dual::<
                        f64, $G<helicoid_linalg::Dual<f64, $D>>, $D>(&a, &b, &c),
                    $B.dual_rows,
                )?;
            }
        }
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
                        // The three poison lanes cycle over `DOF`; at `DOF = 3` this is the array.
                        let p = [f64::NAN, f64::INFINITY, -1e300];
                        core::array::from_fn(|i| if poison {
                            helicoid_linalg::Dual {
                                v: v[i],
                                d: core::array::from_fn(|l| p[l % 3]),
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

/// The harness's own rules, which no group's laws would show.
mod self_test {
    use super::{den, e, within, worst};

    #[test]
    fn a_nan_reference_does_not_become_a_denominator_of_one() {
        // The trap: `f64::max` returns the *other* operand for a NaN.
        assert_eq!(f64::NAN.max(1.0).to_bits(), 1.0_f64.to_bits());
        assert!(den(&[f64::NAN]).is_nan());
        assert_eq!(den(&[0.5]).to_bits(), 1.0_f64.to_bits());
        assert_eq!(den(&[3.0, 4.0]).to_bits(), 5.0_f64.to_bits());
        assert!(e::<f64>(&[1.0], &[f64::NAN]).is_nan());
        assert!(worst(f64::NAN, 0.0).is_nan() && worst(0.0, f64::NAN).is_nan());
    }

    #[test]
    fn a_nan_error_is_within_no_bound() {
        assert!(within(f64::NAN, f64::INFINITY).is_err());
        assert!(within(0.0, 0.0).is_ok());
    }
}
