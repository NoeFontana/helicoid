//! Tests of the fixed-size types. Layout and small cases are exact; the algebra runs under proptest
//! for `f64`, `f32` and `Dual<f64, 2>` on entries `m 2^e`, `|m| < 1`, `|e| <= 6`.
//!
//! An approximate identity is checked entrywise, per lane (the value and each derivative), against
//! `bound(lane) * scale`, where `scale` is the same computation on the absolute values of every
//! lane (the "shadow": a sum of non-negative terms, so it carries no cancellation) and `bound` is
//! the standard forward-error model, `gamma_k = k u / (1 - k u)` (Higham, *Accuracy and Stability
//! of Numerical Algorithms*, Thm 3.5) with one more rounding per multiplication level in a
//! derivative lane. A test's doc records the worst `|diff| / (bound * scale)` for `f64`, `f32` and `D`, all
//! lanes, over 10^6 cases from a fixed seed, so a bound that stops being tight shows up in review:
//!
//! ```text
//! PROPTEST_CASES=1000000 PROPTEST_RNG_SEED=1 cargo nextest run --release -p helicoid-linalg linalg_tests
//! ```
//!
//! Each figure is that single run's maximum, not a bound: the derivative lanes never exceeded the
//! value lane's ratio, so `D` reads as `f64`.

use crate::{hat, vee, Blend, Dual, Mat2, Mat3, Matrix, Point, Precision, Real, Vec3, Vector};
use core::array;
use proptest::prelude::*;
use std::vec::Vec;

pub(crate) type D = Dual<f64, 2>;

/// A scalar whose lanes tests can read and build: `f64` and `f32` have one, `D` has three.
pub(crate) trait Lane: Real {
    const LANES: usize;
    fn make(v: f64, d: [f64; 2]) -> Self;
    fn lane(self, k: usize) -> f64;
}

impl Lane for f64 {
    const LANES: usize = 1;
    fn make(v: f64, _: [f64; 2]) -> Self {
        v
    }
    fn lane(self, k: usize) -> f64 {
        if k == 0 {
            self
        } else {
            0.0
        }
    }
}

impl Lane for f32 {
    const LANES: usize = 1;
    fn make(v: f64, _: [f64; 2]) -> Self {
        v as f32
    }
    fn lane(self, k: usize) -> f64 {
        if k == 0 {
            f64::from(self)
        } else {
            0.0
        }
    }
}

impl Lane for D {
    const LANES: usize = 3;
    fn make(v: f64, d: [f64; 2]) -> Self {
        Dual { v, d }
    }
    fn lane(self, k: usize) -> f64 {
        match k {
            0 => self.v,
            k => self.d.get(k - 1).copied().unwrap_or(0.0),
        }
    }
}

pub(crate) fn unit<S: Real>() -> f64 {
    match S::PRECISION {
        Precision::F64 => f64::EPSILON / 2.0,
        Precision::F32 => f64::from(f32::EPSILON) / 2.0,
    }
}

pub(crate) fn gamma<S: Real>(k: usize) -> f64 {
    let ku = k as f64 * unit::<S>();
    ku / (1.0 - ku)
}

/// Extra roundings a derivative lane carries per multiplication level: `a.d b.v + a.v b.d` is two
/// products and a sum where the value lane has one product.
pub(crate) fn extra(lane: usize, levels: usize) -> usize {
    if lane == 0 {
        0
    } else {
        levels
    }
}

pub(crate) fn shadow<S: Lane>(x: S) -> D {
    Dual {
        v: x.lane(0).abs(),
        d: [x.lane(1).abs(), x.lane(2).abs()],
    }
}

pub(crate) fn shadow_m<S: Lane, const R: usize, const C: usize>(
    m: &Matrix<S, R, C>,
) -> Matrix<D, R, C> {
    Matrix::from_cols(array::from_fn(|c| {
        Vector(array::from_fn(|r| shadow(m.get(r, c))))
    }))
}

pub(crate) fn shadow_v<S: Lane, const N: usize>(v: Vector<S, N>) -> Vector<D, N> {
    Vector(v.0.map(shadow))
}

pub(crate) fn flat<S: Real, const R: usize, const C: usize>(m: &Matrix<S, R, C>) -> Vec<S> {
    (0..R * C).map(|i| m.get(i % R, i / R)).collect()
}

fn col_of<S: Real, const N: usize>(v: Vector<S, N>) -> Matrix<S, N, 1> {
    Matrix::from_cols([v])
}

/// The worst `|x - y| / (bound(lane) scale)` over entries and lanes; NaN propagates.
pub(crate) fn ratio<S: Lane>(x: &[S], y: &[S], scale: &[D], bound: &dyn Fn(usize) -> f64) -> f64 {
    let mut worst = 0.0_f64;
    for ((&a, &b), s) in x.iter().zip(y).zip(scale) {
        for k in 0..S::LANES {
            let diff = (a.lane(k) - b.lane(k)).abs();
            let q = if diff.to_bits() == 0 {
                0.0
            } else {
                diff / (bound(k) * s.lane(k))
            };
            worst = if q.is_nan() || worst.is_nan() {
                f64::NAN
            } else {
                worst.max(q)
            };
        }
    }
    worst
}

/// Bitwise equality of every lane; with `zeros`, `+0` and `-0` agree.
fn same<S: Lane>(x: &[S], y: &[S], zeros: bool) -> bool {
    x.iter().zip(y).all(|(&a, &b)| {
        (0..S::LANES).all(|k| {
            let (a, b) = (a.lane(k), b.lane(k));
            a.to_bits() == b.to_bits()
                || (zeros && a.abs().to_bits() == 0 && b.abs().to_bits() == 0)
        })
    })
}

pub(crate) fn within(name: &str, worst: f64) -> Result<(), TestCaseError> {
    prop_assert!(worst <= 1.0, "{name}: {worst} of the bound");
    Ok(())
}

/// The entries the algebra tests draw from: `(m 2^e, d0, d1)`.
pub(crate) struct Pool<'a>(pub(crate) core::slice::Iter<'a, [f64; 3]>);

impl Pool<'_> {
    pub(crate) fn next<S: Lane>(&mut self) -> S {
        let entry = self.0.next().copied();
        debug_assert!(entry.is_some(), "the pool of {POOL} entries is exhausted");
        let [v, a, b] = entry.unwrap_or([0.5, 0.0, 0.0]);
        S::make(v, [a, b])
    }
    pub(crate) fn vec<S: Lane, const N: usize>(&mut self) -> Vector<S, N> {
        Vector(array::from_fn(|_| self.next()))
    }
    pub(crate) fn mat<S: Lane, const R: usize, const C: usize>(&mut self) -> Matrix<S, R, C> {
        Matrix::from_cols(array::from_fn(|_| self.vec()))
    }
}

pub(crate) const POOL: usize = 64;

pub(crate) fn pool() -> impl Strategy<Value = Vec<[f64; 3]>> {
    let entry = (-1.0_f64..1.0, -6_i32..=6, -1.0_f64..1.0, -1.0_f64..1.0)
        .prop_map(|(m, e, a, b)| [libm::ldexp(m, e), a, b]);
    prop::collection::vec(entry, POOL)
}

/// `max_global_rejects` is raised for the `prop_assume!` in `adj_exact` (singular integer
/// matrices, about 10%) and `adj_residual` (nearly singular ones, about 5%): proptest's default of
/// 1024 aborts a run of a few 10^4 cases.
pub(crate) fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 512,
        max_global_rejects: 1 << 24,
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

/// `(AB)C` against `A(BC)`. Each side is within `gamma_K + gamma_L + gamma_K gamma_L` of the exact
/// product times `|A||B||C|`, so they differ by twice that. Measured worst: 0.59 (f64, D) and 0.63
/// (f32) for 3x3x3x3; 0.54 (f64, D) and 0.55 (f32) for 2x3x4x2.
fn assoc<S: Lane, const R: usize, const K: usize, const L: usize, const C: usize>(
    p: &mut Pool,
) -> Result<(), TestCaseError> {
    let (a, b, c) = (p.mat::<S, R, K>(), p.mat::<S, K, L>(), p.mat::<S, L, C>());
    let (lhs, rhs) = ((a * b) * c, a * (b * c));
    let scale = (shadow_m(&a) * shadow_m(&b)) * shadow_m(&c);
    let bound = |k| {
        let (gk, gl) = (gamma::<S>(K + extra(k, 1)), gamma::<S>(L + extra(k, 1)));
        2.0 * (gk + gl + gk * gl)
    };
    within(
        "assoc",
        ratio(&flat(&lhs), &flat(&rhs), &flat(&scale), &bound),
    )
}

/// `(Ax) . y` against `x . (A^T y)`: each is `R C` triple products with `R + C` roundings on the
/// deepest path, and two multiplication levels. Measured worst: 0.56 (f64, D), 0.54 (f32).
fn adjoint<S: Lane, const R: usize, const C: usize>(p: &mut Pool) -> Result<(), TestCaseError> {
    let (a, x, y) = (p.mat::<S, R, C>(), p.vec::<S, C>(), p.vec::<S, R>());
    let (lhs, rhs) = ((a * x).dot(y), x.dot(a.transpose() * y));
    let scale = (shadow_m(&a) * shadow_v(x)).dot(shadow_v(y));
    let bound = |k| 2.0 * gamma::<S>(R + C + extra(k, 2));
    within("adjoint", ratio(&[lhs], &[rhs], &[scale], &bound))
}

/// Rounding is symmetric in the operands, so these hold to the bit (`+0`/`-0` aside where noted).
fn exact_symmetries<S: Lane>(p: &mut Pool) -> Result<(), TestCaseError> {
    let (a, b, c) = (p.mat::<S, 2, 3>(), p.mat::<S, 3, 4>(), p.vec::<S, 3>());
    prop_assert!(same(
        &flat(&(a * b).transpose()),
        &flat(&(b.transpose() * a.transpose())),
        false
    ));
    prop_assert!(same(&flat(&a.transpose().transpose()), &flat(&a), false));
    prop_assert!(same(&flat(&col_of(a * c)), &flat(&(a * col_of(c))), false));
    let (u, v, w) = (p.vec::<S, 4>(), p.vec::<S, 4>(), p.vec::<S, 3>());
    let z = Vector([S::zero(); 4]);
    prop_assert!(same(&(u + v).0, &(v + u).0, false) && same(&(u - v).0, &(u + -v).0, false));
    prop_assert!(same(&((u - v) + (v - u)).0, &z.0, true) && same(&[u.dot(v)], &[v.dot(u)], false));
    prop_assert!(
        same(&[u.norm_sq()], &[u.dot(u)], false) && same(&[u.norm()], &[u.norm_sq().sqrt()], false)
    );
    let (x, y) = (p.vec::<S, 3>(), p.vec::<S, 3>());
    let o = Vector([S::zero(); 3]);
    prop_assert!(same(&x.cross(y).0, &(-y.cross(x)).0, true) && same(&x.cross(x).0, &o.0, true));
    // hat(a) b is a x b; the zero terms of hat only ever add `+0`.
    prop_assert!(same(&(hat(x) * y).0, &x.cross(y).0, true));
    prop_assert!(same(&(hat(x) * x).0, &o.0, true));
    prop_assert!(same(&flat(&hat(x).transpose()), &flat(&-hat(x)), true));
    prop_assert!(same(&vee(hat(w)).0, &w.0, false));
    Ok(())
}

/// `a . (a x b) = 0`. `c = a x b` has componentwise error `gamma_2 (|p| + |q|)` and the dot adds
/// `gamma_3`, so `gamma_5` of `sum_i |a_i| (|p_i| + |q_i|)`. Measured worst: 0.52 (f64, D), 0.49
/// (f32).
fn cross_orthogonal<S: Lane>(p: &mut Pool) -> Result<(), TestCaseError> {
    let (a, b) = (p.vec::<S, 3>(), p.vec::<S, 3>());
    let ([a0, a1, a2], [b0, b1, b2]) = (shadow_v(a).0, shadow_v(b).0);
    let abs_cross = Vector([a1 * b2 + a2 * b1, a2 * b0 + a0 * b2, a0 * b1 + a1 * b0]);
    let scale = shadow_v(a).dot(abs_cross);
    let bound = |k| gamma::<S>(5 + extra(k, 2));
    within(
        "cross",
        ratio(&[a.dot(a.cross(b))], &[S::zero()], &[scale], &bound),
    )
}

/// `hat(a x b) = [hat a, hat b]`, the Jacobi identity of `so(3)`, holds to the bit: an off-diagonal
/// entry of either product is one product of components, the diagonal ones are `-(a_i b_i + a_j
/// b_j)` in both orders, and the remaining terms are `+0`.
fn hat_commutator<S: Lane>(p: &mut Pool) -> Result<(), TestCaseError> {
    let (a, b) = (p.vec::<S, 3>(), p.vec::<S, 3>());
    let (ha, hb) = (hat(a), hat(b));
    prop_assert!(same(
        &flat(&hat(a.cross(b))),
        &flat(&(ha * hb - hb * ha)),
        true
    ));
    Ok(())
}

fn small_ints() -> impl Strategy<Value = [i32; 9]> {
    prop::array::uniform9(-4_i32..=4)
}

/// `Mat3::inverse_adj` on integer matrices: `det` and `adj` are exact integers (independent
/// cofactor expansion in `i64`), so the result is the correctly rounded quotient, to the bit.
fn adj_exact<S: Lane>(m: [i32; 9]) -> Result<(), TestCaseError> {
    let e = |r: usize, c: usize| i64::from(m[3 * c + r]);
    let minor = |r: usize, c: usize| {
        let (r0, r1) = ((r + 1) % 3, (r + 2) % 3);
        let (c0, c1) = ((c + 1) % 3, (c + 2) % 3);
        e(r0, c0) * e(r1, c1) - e(r0, c1) * e(r1, c0)
    };
    // Cyclic minors carry the cofactor sign already.
    let det = e(0, 0) * minor(0, 0) + e(0, 1) * minor(0, 1) + e(0, 2) * minor(0, 2);
    let a = Mat3::from_cols(array::from_fn(|c| {
        Vector(array::from_fn(|r| S::make(e(r, c) as f64, [0.0; 2])))
    }));
    let (inv, d) = a.inverse_adj();
    prop_assert!(same(&[d], &[S::make(det as f64, [0.0; 2])], true));
    prop_assume!(det != 0);
    for r in 0..3 {
        for c in 0..3 {
            // adj(A)[r][c] is the cofactor of entry (c, r).
            let want = S::make(minor(c, r) as f64, [0.0; 2]) / S::make(det as f64, [0.0; 2]);
            prop_assert!(same(&[inv.get(r, c)], &[want], true), "({r}, {c})");
        }
    }
    Ok(())
}

/// `(s, t)` with `s = |A| adjabs`, `adjabs` the adjugate of `|A|` with every cofactor's two products
/// added in absolute value, so `s` bounds the rounding errors of `A adj`, and `t = tr s`.
fn adj_shadow<S: Lane>(a: &Mat3<S>) -> (Mat3<D>, D) {
    let sa = shadow_m(a);
    let (c0, c1, c2) = (sa.col(0), sa.col(1), sa.col(2));
    let abs_cross = |x: Vec3<D>, y: Vec3<D>| {
        let ([x0, x1, x2], [y0, y1, y2]) = (x.0, y.0);
        Vector([x1 * y2 + x2 * y1, x2 * y0 + x0 * y2, x0 * y1 + x1 * y0])
    };
    let adjabs = Mat3::from_rows([abs_cross(c1, c2), abs_cross(c2, c0), abs_cross(c0, c1)]);
    let s = sa * adjabs;
    (s, s.get(0, 0) + s.get(1, 1) + s.get(2, 2))
}

/// `A * inv - I` on a matrix with `|det| >= tr(|A| |adj|) / 16` (not nearly singular). With
/// `inv = adj/det`, `A inv - I = (A E - e_det I) / det` for the errors `E` of `adj` (`gamma_2`
/// of `|A|`-cross terms, so `|A| adjabs`) and `e_det` of `det` (`gamma_3` of `tr(|A| adjabs)`),
/// plus the division (`gamma_1`) and the product (`gamma_3`): `gamma_9` of
/// `(|A| adjabs + tr I) / |det|`, the quotient taken lane-wise as the shadow of `x / y`.
/// Measured worst: 0.33 (f64, D), 0.33 (f32).
fn adj_residual<S: Lane>(p: &mut Pool) -> Result<(), TestCaseError> {
    let a = p.mat::<S, 3, 3>();
    let (inv, det) = a.inverse_adj();
    let (s, t) = adj_shadow(&a);
    let sd = shadow(det);
    prop_assume!(sd.v * 16.0 >= t.v);
    let quot = |x: D| D {
        v: x.v / sd.v,
        d: array::from_fn(|i| x.d[i] / sd.v + x.v * sd.d[i] / (sd.v * sd.v)),
    };
    let scale: Mat3<D> = Matrix::from_cols(array::from_fn(|c| {
        Vector(array::from_fn(|r| {
            quot(s.get(r, c) + if r == c { t } else { D::constant(0.0) })
        }))
    }));
    let bound = |k| gamma::<S>(9 + extra(k, 7));
    within(
        "adjugate",
        ratio(
            &flat(&(a * inv)),
            &flat(&Mat3::<S>::identity()),
            &flat(&scale),
            &bound,
        ),
    )
}

/// `A adj(A) = det I` without a division, on singular and nearly singular `A`: the third column is
/// `alpha c0 + beta c1 + eps n` with `eps` of `sqrt(u)`, `u` and `0`, so `det` is nearly or exactly
/// zero. `adj` is read back as `inv det` (two more roundings), and an exactly zero `det` is
/// skipped, for `inv` is then not finite. The errors are those of `adj_residual` without the
/// division: `gamma_8` of `s + tr(s) I` (`adj_shadow`). Value lanes only: a derivative lane of
/// `inv det` carries `1 / det^2` factors that no bound of this form covers. Measured worst: 0.44
/// (f64), 0.44 (f32).
fn adj_identity<S: Lane>(p: &mut Pool) -> Result<(), TestCaseError> {
    let (c0, c1) = (p.vec::<S, 3>(), p.vec::<S, 3>());
    let (alpha, beta, noise) = (p.next::<S>(), p.next::<S>(), p.vec::<S, 3>());
    for eps in [unit::<S>().sqrt(), unit::<S>(), 0.0] {
        let c2 = c0.scale(alpha) + c1.scale(beta) + noise.scale(S::make(eps, [0.0; 2]));
        let a = Mat3::from_cols([c0, c1, c2]);
        let (inv, det) = a.inverse_adj();
        if det.lane(0).abs().to_bits() == 0 {
            continue;
        }
        let residual = a * inv.scale(det) - Mat3::identity().scale(det);
        let (s, t) = adj_shadow(&a);
        let scale: Mat3<D> = Matrix::from_cols(array::from_fn(|c| {
            Vector(array::from_fn(|r| {
                s.get(r, c) + if r == c { t } else { D::constant(0.0) }
            }))
        }));
        let bound = |_| gamma::<S>(8);
        within(
            "adjugate identity",
            ratio(&flat(&residual), &[S::zero(); 9], &flat(&scale), &bound),
        )?;
    }
    Ok(())
}

fn small_ints2() -> impl Strategy<Value = [i32; 4]> {
    prop::array::uniform4(-9_i32..=9)
}

/// `Mat2::inverse_adj` on integer matrices: `det` is an exact integer and the adjugate is a
/// permutation with signs, so the result is the correctly rounded quotient, to the bit (`0061`).
fn adj2_exact<S: Lane>(m: [i32; 4]) -> Result<(), TestCaseError> {
    let e = |r: usize, c: usize| i64::from(m[2 * c + r]);
    let det = e(0, 0) * e(1, 1) - e(0, 1) * e(1, 0);
    let a = Mat2::from_cols(array::from_fn(|c| {
        Vector(array::from_fn(|r| S::make(e(r, c) as f64, [0.0; 2])))
    }));
    let (inv, d) = a.inverse_adj();
    prop_assert!(same(&[d], &[S::make(det as f64, [0.0; 2])], true));
    prop_assume!(det != 0);
    let adj = [[e(1, 1), -e(0, 1)], [-e(1, 0), e(0, 0)]];
    for (r, row) in adj.iter().enumerate() {
        for (c, &x) in row.iter().enumerate() {
            let want = S::make(x as f64, [0.0; 2]) / S::make(det as f64, [0.0; 2]);
            prop_assert!(same(&[inv.get(r, c)], &[want], true), "({r}, {c})");
        }
    }
    Ok(())
}

/// `(s, t)` with `s = |A| |adj A|` and `t = |m00 m11| + |m01 m10|`, which bounds the rounding
/// error of `det` (`gamma_2` of it).
fn adj2_shadow<S: Lane>(a: &Mat2<S>) -> (Mat2<D>, D) {
    let sa = shadow_m(a);
    let (p, q, r, w) = (sa.get(0, 0), sa.get(0, 1), sa.get(1, 0), sa.get(1, 1));
    let adjabs = Mat2::from_rows([Vector([w, q]), Vector([r, p])]);
    (sa * adjabs, p * w + q * r)
}

/// `A * inv - I` on a matrix with `|det| >= t / 16` (not nearly singular). The adjugate is exact,
/// so the errors are `det`'s (`gamma_2` of `t`), the division (`gamma_1`) and the two-term product
/// (`gamma_2`): `gamma_5` of `(|A| |adj| + t I) / |det|`, the quotient taken lane-wise as the
/// shadow of `x / y`. Measured worst: 0.40 (f64, D), 0.40 (f32).
fn adj2_residual<S: Lane>(p: &mut Pool) -> Result<(), TestCaseError> {
    let a = p.mat::<S, 2, 2>();
    let (inv, det) = a.inverse_adj();
    let (s, t) = adj2_shadow(&a);
    let sd = shadow(det);
    prop_assume!(sd.v * 16.0 >= t.v);
    let quot = |x: D| D {
        v: x.v / sd.v,
        d: array::from_fn(|i| x.d[i] / sd.v + x.v * sd.d[i] / (sd.v * sd.v)),
    };
    let scale: Mat2<D> = Matrix::from_cols(array::from_fn(|c| {
        Vector(array::from_fn(|r| {
            quot(s.get(r, c) + if r == c { t } else { D::constant(0.0) })
        }))
    }));
    let bound = |k| gamma::<S>(5 + extra(k, 5));
    within(
        "adjugate 2x2",
        ratio(
            &flat(&(a * inv)),
            &flat(&Mat2::<S>::identity()),
            &flat(&scale),
            &bound,
        ),
    )
}

/// `A adj(A) = det I` without a division, on singular and nearly singular `A`: the second column
/// is `alpha c0 + eps n` with `eps` of `sqrt(u)`, `u` and `0`. `adj` is read back as `inv det`
/// (two more roundings) and an exactly zero `det` is skipped. `gamma_6` of `s + t I`. Value lanes
/// only, as `adj_identity`. Measured worst: 0.42 (f64), 0.40 (f32).
fn adj2_identity<S: Lane>(p: &mut Pool) -> Result<(), TestCaseError> {
    let c0 = p.vec::<S, 2>();
    let (alpha, noise) = (p.next::<S>(), p.vec::<S, 2>());
    for eps in [unit::<S>().sqrt(), unit::<S>(), 0.0] {
        let c1 = c0.scale(alpha) + noise.scale(S::make(eps, [0.0; 2]));
        let a = Mat2::from_cols([c0, c1]);
        let (inv, det) = a.inverse_adj();
        if det.lane(0).abs().to_bits() == 0 {
            continue;
        }
        let residual = a * inv.scale(det) - Mat2::identity().scale(det);
        let (s, t) = adj2_shadow(&a);
        let scale: Mat2<D> = Matrix::from_cols(array::from_fn(|c| {
            Vector(array::from_fn(|r| {
                s.get(r, c) + if r == c { t } else { D::constant(0.0) }
            }))
        }));
        let bound = |_| gamma::<S>(6);
        within(
            "adjugate identity 2x2",
            ratio(&flat(&residual), &[S::zero(); 4], &flat(&scale), &bound),
        )?;
    }
    Ok(())
}

/// `(s A) x` against `s (A x)` and `(s u) . v` against `s (u . v)`: each side has one more product
/// than the plain reduction, so `gamma_{n+1}` of the shadow on each, and two multiplication levels
/// in a derivative lane. Measured worst: 0.60 (f64, D) and 0.62 (f32) for the matrix
/// form; 0.60 (f64, D) and 0.73 (f32) for the vector form.
fn scale_is_bilinear<S: Lane>(p: &mut Pool) -> Result<(), TestCaseError> {
    let (a, x, s) = (p.mat::<S, 3, 4>(), p.vec::<S, 4>(), p.next::<S>());
    let (u, v) = (p.vec::<S, 4>(), p.vec::<S, 4>());
    let sh = shadow(s);
    let bound = |k| 2.0 * gamma::<S>(5 + extra(k, 2));
    let scale = (shadow_m(&a) * shadow_v(x)).scale(sh);
    within(
        "scale, matrix",
        ratio(&(a.scale(s) * x).0, &((a * x).scale(s)).0, &scale.0, &bound),
    )?;
    let scale = shadow_v(u).dot(shadow_v(v)) * sh;
    within(
        "scale, vector",
        ratio(&[u.scale(s).dot(v)], &[u.dot(v) * s], &[scale], &bound),
    )
}

/// A blend returns the selected operand whole, lanes included, and a NaN in the other one does not
/// leak.
fn blends<S: Lane>(p: &mut Pool) -> Result<(), TestCaseError> {
    let nan = S::make(f64::NAN, [f64::NAN; 2]);
    let (m1, m2) = (p.mat::<S, 3, 3>(), p.mat::<S, 3, 3>());
    let (v1, v2) = (p.vec::<S, 3>(), p.vec::<S, 3>());
    let (q1, q2) = (Point(v1.0), Point(v2.0));
    let (mn, vn, qn) = (
        Mat3::from_cols([Vector([nan; 3]); 3]),
        Vector([nan; 3]),
        Point([nan; 3]),
    );
    let (yes, no) = (S::zero().lt(S::one()), S::one().lt(S::zero()));
    for (m, first) in [(yes, true), (no, false)] {
        let (fm, fv, fq) = S::branch(m, || (m1, v1, q1), || (m2, v2, q2));
        let (em, ev, eq) = if first { (m1, v1, q1) } else { (m2, v2, q2) };
        prop_assert!(same(&flat(&fm), &flat(&em), false) && same(&fv.0, &ev.0, false));
        prop_assert!(same(&fq.0, &eq.0, false));
        let (t, f) = if first { (m1, mn) } else { (mn, m2) };
        prop_assert!(same(&flat(&Blend::<S>::blend(m, t, f)), &flat(&em), false));
        let (t, f) = if first { (v1, vn) } else { (vn, v2) };
        prop_assert!(same(&Blend::<S>::blend(m, t, f).0, &ev.0, false));
        let (t, f) = if first { (q1, qn) } else { (qn, q2) };
        prop_assert!(same(&Blend::<S>::blend(m, t, f).0, &eq.0, false));
    }
    Ok(())
}

fn bits(x: &[f64]) -> Vec<u64> {
    x.iter().map(|v| v.to_bits()).collect()
}

fn v3(x: f64, y: f64, z: f64) -> Vec3<f64> {
    Vector([x, y, z])
}

#[test]
fn matrix_is_column_major() {
    let m = Matrix::from_rows([Vector([1.0, 2.0, 3.0]), Vector([4.0, 5.0, 6.0])]);
    assert_eq!(bits(&flat(&m)), bits(&[1.0, 4.0, 2.0, 5.0, 3.0, 6.0]));
    let n = Matrix::from_cols([Vector([1.0, 4.0]), Vector([2.0, 5.0]), Vector([3.0, 6.0])]);
    assert_eq!(bits(&flat(&n)), bits(&flat(&m)));
    assert_eq!(bits(&m.col(1).0), bits(&[2.0, 5.0]));
    assert_eq!(bits(&m.row(1).0), bits(&[4.0, 5.0, 6.0]));
    assert_eq!(m.get(1, 2).to_bits(), 6.0_f64.to_bits());
    let mut t = m.transpose();
    assert_eq!(t.get(2, 1).to_bits(), 6.0_f64.to_bits());
    t.set(0, 1, -9.0);
    assert_eq!(bits(&flat(&t)), bits(&[1.0, 2.0, 3.0, -9.0, 5.0, 6.0]));
}

/// Out of range is a `debug_assert!` (D11), so the debug profile that `just test` runs panics ...
#[cfg(debug_assertions)]
mod out_of_range {
    use super::*;
    use core::hint::black_box;

    #[test]
    #[should_panic(expected = "Matrix::get")]
    fn get() {
        black_box(Matrix::<f64, 2, 2>::identity().get(0, 2));
    }

    #[test]
    #[should_panic(expected = "Matrix::set")]
    fn set() {
        Matrix::<f64, 2, 2>::identity().set(2, 0, 7.0);
    }

    #[test]
    #[should_panic(expected = "Matrix::col")]
    fn col() {
        black_box(Matrix::<f64, 2, 2>::identity().col(2));
    }

    #[test]
    #[should_panic(expected = "Matrix::row")]
    fn row() {
        black_box(Matrix::<f64, 2, 2>::identity().row(2));
    }
}

/// ... and a release build (`cargo nextest run --release`) never panics and leaves the in-range
/// entries alone; the value it returns outside is unspecified (`docs/API.md` R6).
#[cfg(not(debug_assertions))]
#[test]
fn out_of_range_access_does_not_panic_in_release() {
    let mut m = Matrix::<f64, 2, 2>::identity();
    m.set(2, 0, 7.0);
    core::hint::black_box((m.get(0, 2), m.col(2), m.row(2)));
    assert_eq!(bits(&flat(&m)), bits(&[1.0, 0.0, 0.0, 1.0]));
}

/// Every reduction adds left to right from the first term, and the bits show it where the order
/// is visible: `(1 + 1e16) - 1e16` is `0` and `(1e16 - 1e16) + 1` is `1`. A lone `-0` stays `-0`
/// because the sum starts at the first term, not at `+0` (D16).
#[test]
fn reductions_add_left_to_right_from_the_first_term() {
    let big = 1e16;
    let (up, down, ones) = (v3(1.0, big, -big), v3(big, -big, 1.0), v3(1.0, 1.0, 1.0));
    assert_eq!(bits(&[up.dot(ones), down.dot(ones)]), bits(&[0.0, 1.0]));
    let m = Matrix::from_rows([up, down]);
    assert_eq!(bits(&(m * ones).0), bits(&[0.0, 1.0]));
    assert_eq!(
        bits(&flat(&(m * Matrix::from_cols([ones, ones])))),
        bits(&[0.0, 1.0, 0.0, 1.0])
    );
    // `c1 x c2 = (1, -1, 1)`, so the terms of `det` are `1, 1e16, -1e16` in that order.
    let a = Mat3::from_cols([v3(1.0, -big, -big), v3(1.0, 1.0, 0.0), v3(0.0, 1.0, 1.0)]);
    assert_eq!(a.inverse_adj().1.to_bits(), 0);
    let z = -0.0_f64;
    let (lone, pair) = (Vector([z]), Vector([z, z]));
    assert_eq!(lone.dot(Vector([1.0])).to_bits(), z.to_bits());
    assert_eq!(pair.dot(Vector([1.0, 1.0])).to_bits(), z.to_bits());
    assert_eq!(
        (Matrix::from_cols([lone]) * Vector([1.0])).0[0].to_bits(),
        z.to_bits()
    );
    let id = Mat3::<f64>::identity();
    assert_eq!(
        bits(&flat(&id)),
        bits(&[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])
    );
}

/// The magnitude ranges documented on `norm` and `inverse_adj` are met exactly on powers of two:
/// the last exponents at which `norm_sq`, `det` and `adj` stay normal.
fn documented_range<S: Lane>(norm_exp: i32, inv_exp: i32) -> Result<(), TestCaseError> {
    let pow = |e: i32| S::make(libm::ldexp(1.0, e), [0.0; 2]);
    for e in [norm_exp, -norm_exp] {
        let n = Vector([pow(e), S::zero(), S::zero()]).norm();
        prop_assert!(same(&[n], &[pow(e)], false), "norm 2^{e}");
    }
    for e in [inv_exp, -inv_exp] {
        let (inv, det) = Mat3::<S>::identity().scale(pow(e)).inverse_adj();
        prop_assert!(same(&[det], &[pow(3 * e)], false), "det 2^{e}");
        prop_assert!(same(&[inv.get(1, 1)], &[pow(-e)], false), "inverse 2^{e}");
        prop_assert!(same(&[inv.get(0, 1)], &[S::zero()], true), "inverse 2^{e}");
    }
    Ok(())
}

/// `Mat2::inverse_adj`'s range (`det` quadratic): the last exponents at which `det` stays normal.
fn documented_range2<S: Lane>(inv_exp: i32) -> Result<(), TestCaseError> {
    let pow = |e: i32| S::make(libm::ldexp(1.0, e), [0.0; 2]);
    for e in [inv_exp, -inv_exp] {
        let (inv, det) = Mat2::<S>::identity().scale(pow(e)).inverse_adj();
        prop_assert!(same(&[det], &[pow(2 * e)], false), "det 2^{e}");
        prop_assert!(same(&[inv.get(1, 1)], &[pow(-e)], false), "inverse 2^{e}");
        prop_assert!(same(&[inv.get(0, 1)], &[S::zero()], true), "inverse 2^{e}");
    }
    Ok(())
}

#[test]
fn documented_ranges_hold() -> Result<(), TestCaseError> {
    documented_range::<f64>(500, 330)?;
    documented_range::<f32>(60, 40)?;
    documented_range2::<f64>(510)?;
    documented_range2::<f32>(63)
}

#[test]
fn small_products_are_exact() {
    let a = Matrix::from_rows([Vector([1.0, 2.0]), Vector([3.0, 4.0])]);
    let b = Matrix::from_rows([Vector([5.0, 6.0]), Vector([7.0, 8.0])]);
    assert_eq!(bits(&flat(&(a * b))), bits(&[19.0, 43.0, 22.0, 50.0]));
    assert_eq!(bits(&(a * Vector([5.0, 6.0])).0), bits(&[17.0, 39.0]));
    assert_eq!(bits(&flat(&(a * Matrix::identity()))), bits(&flat(&a)));
    assert_eq!(bits(&flat(&(a + b - a))), bits(&flat(&b)));
    assert_eq!(
        bits(&flat(&(-a).scale(2.0))),
        bits(&[-2.0, -6.0, -4.0, -8.0])
    );
    assert_eq!(
        bits(&v3(1.0, 2.0, 3.0).scale(2.0).0),
        bits(&[2.0, 4.0, 6.0])
    );
    let (e1, e2) = (v3(1.0, 0.0, 0.0), v3(0.0, 1.0, 0.0));
    assert_eq!(bits(&e1.cross(e2).0), bits(&[0.0, 0.0, 1.0]));
    assert_eq!(
        bits(&hat(v3(1.0, 2.0, 3.0)).row(0).0),
        bits(&[0.0, -3.0, 2.0])
    );
    assert_eq!(
        bits(&[v3(3.0, 4.0, 0.0).norm(), v3(1.0, 2.0, 3.0).norm_sq()]),
        bits(&[5.0, 14.0])
    );
}

#[test]
fn points_subtract_to_vectors_and_translate() {
    let (p, q) = (Point([1.0, 2.0, 3.0]), Point([4.0, 6.0, 8.0]));
    let d = q - p;
    assert_eq!(bits(&d.0), bits(&[3.0, 4.0, 5.0]));
    assert_eq!(bits(&(p + d).0), bits(&q.0));
}

/// The `Dual` rules reach through the fixed-size types: `|(3 x, 4 y)|` has gradient `(3, 4) / 5`.
#[test]
fn norm_differentiates_through_dual() {
    let v = Vector([D::variable(3.0, 0), D::variable(4.0, 1)]);
    let n = v.norm();
    assert_eq!(bits(&[n.v, n.d[0], n.d[1]]), bits(&[5.0, 0.6, 0.8]));
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn matmul_associates_to_the_rounding_bound(pl in pool()) {
        assoc::<f64, 3, 3, 3, 3>(&mut Pool(pl.iter()))?;
        assoc::<f32, 3, 3, 3, 3>(&mut Pool(pl.iter()))?;
        assoc::<D, 3, 3, 3, 3>(&mut Pool(pl.iter()))?;
        assoc::<f64, 2, 3, 4, 2>(&mut Pool(pl.iter()))?;
        assoc::<f32, 2, 3, 4, 2>(&mut Pool(pl.iter()))?;
        assoc::<D, 2, 3, 4, 2>(&mut Pool(pl.iter()))?;
    }

    #[test]
    fn matvec_is_adjoint_to_transpose(pl in pool()) {
        adjoint::<f64, 3, 4>(&mut Pool(pl.iter()))?;
        adjoint::<f32, 3, 4>(&mut Pool(pl.iter()))?;
        adjoint::<D, 3, 4>(&mut Pool(pl.iter()))?;
    }

    #[test]
    fn products_and_cross_have_exact_symmetries(pl in pool()) {
        exact_symmetries::<f64>(&mut Pool(pl.iter()))?;
        exact_symmetries::<f32>(&mut Pool(pl.iter()))?;
        exact_symmetries::<D>(&mut Pool(pl.iter()))?;
    }

    #[test]
    fn cross_is_orthogonal_to_its_factors(pl in pool()) {
        cross_orthogonal::<f64>(&mut Pool(pl.iter()))?;
        cross_orthogonal::<f32>(&mut Pool(pl.iter()))?;
        cross_orthogonal::<D>(&mut Pool(pl.iter()))?;
    }

    #[test]
    fn hat_of_cross_is_the_commutator(pl in pool()) {
        hat_commutator::<f64>(&mut Pool(pl.iter()))?;
        hat_commutator::<f32>(&mut Pool(pl.iter()))?;
        hat_commutator::<D>(&mut Pool(pl.iter()))?;
    }

    #[test]
    fn adjugate_identity_holds_near_singular(pl in pool()) {
        adj_identity::<f64>(&mut Pool(pl.iter()))?;
        adj_identity::<f32>(&mut Pool(pl.iter()))?;
    }

    #[test]
    fn scale_is_bilinear_to_the_rounding_bound(pl in pool()) {
        scale_is_bilinear::<f64>(&mut Pool(pl.iter()))?;
        scale_is_bilinear::<f32>(&mut Pool(pl.iter()))?;
        scale_is_bilinear::<D>(&mut Pool(pl.iter()))?;
    }

    #[test]
    fn blend_selects_whole_values(pl in pool()) {
        blends::<f64>(&mut Pool(pl.iter()))?;
        blends::<f32>(&mut Pool(pl.iter()))?;
        blends::<D>(&mut Pool(pl.iter()))?;
    }

    #[test]
    fn inverse_adj_is_exact_on_integers(m in small_ints()) {
        adj_exact::<f64>(m)?;
        adj_exact::<f32>(m)?;
        adj_exact::<D>(m)?;
    }

    #[test]
    fn inverse_adj2_is_exact_on_integers(m in small_ints2()) {
        adj2_exact::<f64>(m)?;
        adj2_exact::<f32>(m)?;
        adj2_exact::<D>(m)?;
    }

    #[test]
    fn matrix2_times_inverse_adj_is_identity(pl in pool()) {
        adj2_residual::<f64>(&mut Pool(pl.iter()))?;
        adj2_residual::<f32>(&mut Pool(pl.iter()))?;
        adj2_residual::<D>(&mut Pool(pl.iter()))?;
    }

    #[test]
    fn adjugate2_identity_holds_near_singular(pl in pool()) {
        adj2_identity::<f64>(&mut Pool(pl.iter()))?;
        adj2_identity::<f32>(&mut Pool(pl.iter()))?;
    }

    #[test]
    fn matrix_times_inverse_adj_is_identity(pl in pool()) {
        adj_residual::<f64>(&mut Pool(pl.iter()))?;
        adj_residual::<f32>(&mut Pool(pl.iter()))?;
        adj_residual::<D>(&mut Pool(pl.iter()))?;
    }
}
