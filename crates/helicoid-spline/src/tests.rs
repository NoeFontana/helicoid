//! The crate's laws: exact tables, twins, `Dual`, continuity, products and Rⁿ (`0114` plan).

extern crate std;

use std::vec::Vec;

use helicoid::{LieGroup, Product, Rn, Tangent, SE3, SO3};
use helicoid_linalg::{Dual, Real, Vector};
use proptest::prelude::*;

use crate::order::Tables;
use crate::reference::{spline_kinematics, spline_value};
use crate::{derivative_gram, Kinematics, Order, Segment, ValidOrder, Weights};

// ---- exact rationals ---------------------------------------------------------------------------

/// A reduced fraction in `i128`: enough for every table at `K <= 6` and `u` on a `2^-4` grid.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Q(i128, i128);

fn gcd(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

impl Q {
    fn new(n: i128, d: i128) -> Self {
        let g = gcd(n, d).max(1) * d.signum();
        Q(n / g, d / g)
    }
    fn int(n: i128) -> Self {
        Q(n, 1)
    }
    fn add(self, o: Self) -> Self {
        Q::new(self.0 * o.1 + o.0 * self.1, self.1 * o.1)
    }
    fn mul(self, o: Self) -> Self {
        Q::new(self.0 * o.0, self.1 * o.1)
    }
    fn f64(self) -> f64 {
        self.0 as f64 / self.1 as f64
    }
}

/// The cardinal B-spline of order `k` on `[0, k)` by the Cox–de Boor recursion: a route that
/// shares nothing with Qin's matrix in `order.rs`.
fn cardinal(k: usize, x: Q) -> Q {
    if k == 1 {
        return Q::int(i128::from(x.0 >= 0 && x.0 < x.1));
    }
    let km1 = Q::int(k as i128 - 1);
    let left = x.mul(cardinal(k - 1, x));
    let shifted = x.add(Q::int(-1));
    let right = Q::int(k as i128)
        .add(Q(-x.0, x.1))
        .mul(cardinal(k - 1, shifted));
    left.add(right).mul(Q(km1.1, km1.0))
}

fn table_basis<const K: usize>(s: usize, u: Q) -> Q
where
    Order<K>: ValidOrder,
{
    // `b_s = λ_s − λ_{s+1}` from the cumulative table, `λ_K = 0`.
    let row = |j: usize| -> Q {
        if j >= K {
            return Q::int(0);
        }
        let mut acc = Q::int(0);
        let mut p = Q::int(1);
        for k in 0..K {
            acc = acc.add(p.mul(Q::int(i128::from(Tables::<K>::CUMULATIVE[j][k]))));
            p = p.mul(u);
        }
        acc.mul(Q(1, i128::from(Tables::<K>::FACTORIAL)))
    };
    row(s).add(row(s + 1).mul(Q::int(-1)))
}

fn basis_is_cox_de_boor<const K: usize>()
where
    Order<K>: ValidOrder,
{
    for m in 0..16 {
        let u = Q::new(m, 16);
        for s in 0..K {
            // On segment `[0, 1)`, control `s` carries the cardinal spline at `u + K - 1 - s`.
            let want = cardinal(K, u.add(Q::int((K - 1 - s) as i128)));
            assert_eq!(table_basis::<K>(s, u), want, "K = {K}, s = {s}, u = {m}/16");
        }
    }
}

#[test]
fn the_basis_is_cox_de_boor_at_every_order() {
    basis_is_cox_de_boor::<2>();
    basis_is_cox_de_boor::<3>();
    basis_is_cox_de_boor::<4>();
    basis_is_cox_de_boor::<5>();
    basis_is_cox_de_boor::<6>();
}

/// `|got − want| / max(|want|, 1)` in units of `eps`.
fn ulps(got: f64, want: f64, eps: f64) -> f64 {
    (got - want).abs() / want.abs().max(1.0) / eps
}

fn weights_match<const K: usize>()
where
    Order<K>: ValidOrder,
{
    let f = Q::int(i128::from(Tables::<K>::FACTORIAL));
    let inv_f = Q(f.1, f.0);
    // `u` on a `2^-4` grid over `[-1, 2]`: inside the segment and extrapolated (R6).
    for m in -16..=32 {
        let u = Q::new(m, 16);
        let w64 = Weights::<f64, K, 2>::at(u.f64(), 4.0);
        let w32 = Weights::<f32, K, 2>::at(u.f64() as f32, 4.0);
        for j in 1..K {
            let (mut l, mut ld, mut ldd) = (Q::int(0), Q::int(0), Q::int(0));
            for k in 0..K {
                let c = Q::int(i128::from(Tables::<K>::CUMULATIVE[j][k]));
                let pow = |e: usize| (0..e).fold(Q::int(1), |p, _| p.mul(u));
                l = l.add(c.mul(pow(k)));
                if k >= 1 {
                    ld = ld.add(c.mul(Q::int(k as i128)).mul(pow(k - 1)));
                }
                if k >= 2 {
                    ldd = ldd.add(c.mul(Q::int((k * (k - 1)) as i128)).mul(pow(k - 2)));
                }
            }
            let want = [
                l.mul(inv_f),
                ld.mul(inv_f).mul(Q::int(4)),
                ldd.mul(inv_f).mul(Q::int(16)),
            ];
            let got64 = [w64.lambda[j], w64.dot[j], w64.ddot[j]];
            let got32 = [w32.lambda[j], w32.dot[j], w32.ddot[j]];
            for r in 0..3 {
                // Horner over at most 5 steps on a `|u| <= 2` grid: a few ulps of the row's scale.
                let e64 = ulps(got64[r], want[r].f64(), f64::EPSILON);
                let e32 = ulps(f64::from(got32[r]), want[r].f64(), f64::from(f32::EPSILON));
                assert!(
                    e64 <= 64.0 && e32 <= 64.0,
                    "K={K} j={j} r={r} u={m}/16: {e64} {e32}"
                );
            }
        }
    }
}

#[test]
fn weights_match_exact_rationals() {
    weights_match::<2>();
    weights_match::<3>();
    weights_match::<4>();
    weights_match::<5>();
    weights_match::<6>();
}

#[test]
fn weights_above_d_are_zero() {
    let w = Weights::<f64, 5, 0>::at(0.3, 7.0);
    assert!(w.dot.iter().chain(w.ddot.iter()).all(|x| x.to_bits() == 0));
    let w = Weights::<f64, 5, 1>::at(0.3, 7.0);
    assert!(w.ddot.iter().all(|x| x.to_bits() == 0));
    assert_eq!(w.lambda[0].to_bits(), 1.0_f64.to_bits());
}

/// `inv_h`'s domain is a `debug_assert!` (D11): `+∞` passes a bare `0 < inv_h`.
#[cfg(debug_assertions)]
mod out_of_domain {
    use super::*;

    #[test]
    #[should_panic(expected = "finite and positive")]
    fn weights_at_an_infinite_inv_h() {
        let _ = Weights::<f64, 4, 1>::at(0.5, f64::INFINITY);
    }

    #[test]
    #[should_panic(expected = "finite and positive")]
    fn gram_at_an_infinite_inv_h() {
        let _ = derivative_gram::<f64, 4, 1>(f64::INFINITY);
    }
}

fn gram_laws<const K: usize, const R: usize>()
where
    Order<K>: ValidOrder,
{
    let g = derivative_gram::<f64, K, R>(1.0);
    for s in 0..K {
        let row: f64 = (0..K).map(|t| g.get(s, t)).sum();
        assert!(
            row.abs() <= 64.0 * f64::EPSILON,
            "K={K} R={R} row {s}: {row}"
        );
        for t in 0..K {
            assert_eq!(g.get(s, t).to_bits(), g.get(t, s).to_bits(), "symmetric");
        }
    }
    // A power-of-two spacing scales every entry exactly (`0115` decision 3).
    let scaled = derivative_gram::<f64, K, R>(1024.0);
    let f = 1024.0_f64.powi((2 * R - 1) as i32);
    for s in 0..K {
        for t in 0..K {
            assert_eq!(scaled.get(s, t).to_bits(), (g.get(s, t) * f).to_bits());
        }
    }
}

#[test]
fn gram_rows_sum_to_zero_and_scale_exactly() {
    gram_laws::<2, 1>();
    gram_laws::<3, 1>();
    gram_laws::<3, 2>();
    gram_laws::<4, 1>();
    gram_laws::<4, 2>();
    gram_laws::<5, 1>();
    gram_laws::<5, 2>();
    gram_laws::<6, 1>();
    gram_laws::<6, 2>();
    let g = derivative_gram::<f64, 4, 1>(1.0);
    assert_eq!(g.get(0, 0).to_bits(), (1.0_f64 / 20.0).to_bits());
    let g = derivative_gram::<f32, 6, 2>(1.0);
    assert_eq!(g.get(0, 0).to_bits(), (1.0_f32 / 252.0).to_bits());
}

// ---- groups ------------------------------------------------------------------------------------

type R3<S> = Rn<S, 3>;
type Split<S> = Product<SO3<S>, R3<S>>;

/// Dense `[S]` of a tangent.
fn dense<S: Real, T: Tangent<S>>(t: &T) -> Vec<S> {
    let mut v = std::vec![S::zero(); T::DOF];
    t.write_dense(&mut v);
    v
}

fn norm<T: Tangent<f64>>(t: &T) -> f64 {
    dense(t).iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// `‖a − b‖ / max(‖b‖, 1)` in `f64::EPSILON`.
fn tangent_err<T: Tangent<f64>>(a: &T, b: &T) -> f64 {
    norm(&a.sub(b)) / norm(b).max(1.0) / f64::EPSILON
}

fn group_err<G: LieGroup<f64>>(a: &G, b: &G) -> f64 {
    norm(&a.rminus(b)) / f64::EPSILON
}

/// A window from dense tangents: `X_0 = Exp(τ_0)`, `X_j = X_{j-1} Exp(τ_j)`.
fn window<S: Real, G: LieGroup<S>, const K: usize>(taus: &[Vec<S>]) -> [G; K] {
    let mut out = [G::identity(); K];
    let mut x = G::identity();
    for j in 0..K {
        x = x * G::exp(&G::Tangent::read_dense(&taus[j]));
        out[j] = x;
    }
    out
}

fn taus(raw: &[f64], dof: usize, k: usize, step: f64) -> Vec<Vec<f64>> {
    (0..k)
        .map(|j| {
            let s = if j == 0 { 2.0 } else { step };
            (0..dof)
                .map(|i| s * raw[(j * dof + i) % raw.len()])
                .collect()
        })
        .collect()
}

fn twins<G: LieGroup<f64>, const K: usize>(raw: &[f64], u: f64, inv_h: f64, step: f64)
where
    Order<K>: ValidOrder,
{
    let w = Weights::<f64, K, 2>::at(u, inv_h);
    let ctrl = window::<f64, G, K>(&taus(raw, G::Tangent::DOF, K, step));
    let seg = Segment::new(&ctrl);
    let fast = seg.kinematics(&w);
    let slow = spline_kinematics(&ctrl, &w);
    // `slow.value` is `spline_value`: one comparison covers both twins.
    let ev = group_err(&fast.value, &spline_value(&ctrl, &w));
    let evel = tangent_err(&fast.vel, &slow.vel) / inv_h.max(1.0);
    let eacc = tangent_err(&fast.acc, &slow.acc) / (inv_h * inv_h).max(1.0);
    let bits = |g: &G| {
        dense(&g.log())
            .into_iter()
            .map(f64::to_bits)
            .collect::<Vec<_>>()
    };
    assert_eq!(
        bits(&seg.value(&w)),
        bits(&fast.value),
        "value and kinematics share one product"
    );
    // The twin associates the product and the sums differently: bounded, not bitwise.
    assert!(
        ev <= 1e3 && evel <= 1e4 && eacc <= 1e5,
        "K={K} u={u} inv_h={inv_h}: value {ev} vel {evel} acc {eacc}"
    );
}

proptest! {
    #[test]
    fn spline_kinematics_matches_reference(
        raw in prop::collection::vec(-1.0_f64..1.0, 36),
        u in -0.25_f64..1.25,
        inv_h in prop::sample::select(std::vec![1e-3_f64, 1.0, 1e3]),
        step in prop::sample::select(std::vec![1e-9_f64, 1e-3, 0.5, 1.5]),
    ) {
        twins::<SO3<f64>, 4>(&raw, u, inv_h, step);
        twins::<SO3<f64>, 6>(&raw, u, inv_h, step);
        twins::<SE3<f64>, 5>(&raw, u, inv_h, step);
        twins::<R3<f64>, 2>(&raw, u, inv_h, step);
        twins::<Split<f64>, 6>(&raw, u, inv_h, step);
    }
}

// ---- `Dual`: velocity is d/dt of the value, acceleration d/dt of the velocity -----------------

type D1 = Dual<f64, 1>;

fn dual_window<G: LieGroup<D1>, const K: usize>(t: &[Vec<f64>]) -> [G; K] {
    let lifted: Vec<Vec<D1>> = t
        .iter()
        .map(|v| v.iter().map(|&x| D1::constant(x)).collect())
        .collect();
    window::<D1, G, K>(&lifted)
}

fn lanes<T: Tangent<D1>>(t: &T) -> Vec<(f64, f64)> {
    dense(t).iter().map(|x| (x.v, x.d[0])).collect()
}

fn derivatives_by_dual<G: LieGroup<D1>, const K: usize>(raw: &[f64], u: f64, inv_h: f64)
where
    Order<K>: ValidOrder,
{
    let t = taus(raw, G::Tangent::DOF, K, 0.7);
    let ctrl = dual_window::<G, K>(&t);
    let seg = Segment::<D1, G, K>::new(&ctrl);
    let (uv, ih) = (D1::variable(u, 0), D1::constant(inv_h));
    let k1: Kinematics<D1, G> = seg.kinematics(&Weights::<D1, K, 1>::at(uv, ih));
    let k2: Kinematics<D1, G> = seg.kinematics(&Weights::<D1, K, 2>::at(uv, ih));
    // `ω = h⁻¹ d/du Log(X(u₀)⁻¹ X(u))` at `u₀`, where the logarithm's argument is the identity.
    let x0 = {
        let w = Weights::<D1, K, 0>::at(D1::constant(u), D1::constant(inv_h));
        seg.value(&w)
    };
    let local = (x0.inverse() * k1.value).log();
    for ((_, du), (v, _)) in lanes(&local).into_iter().zip(lanes(&k1.vel)) {
        let want = du * inv_h;
        assert!(
            (want - v).abs() <= 1e-12 * v.abs().max(1.0) * inv_h.max(1.0),
            "vel {want} vs {v}"
        );
    }
    // `ω̇ = h⁻¹ dω/du`.
    for ((_, dvel), (a, _)) in lanes(&k1.vel).into_iter().zip(lanes(&k2.acc)) {
        let want = dvel * inv_h;
        let tol = 1e-12 * a.abs().max(1.0) * (inv_h * inv_h).max(1.0);
        assert!((want - a).abs() <= tol, "acc {want} vs {a}");
    }
}

proptest! {
    #[test]
    fn velocity_and_acceleration_are_time_derivatives(
        raw in prop::collection::vec(-1.0_f64..1.0, 36),
        u in 0.0_f64..1.0,
        inv_h in prop::sample::select(std::vec![0.5_f64, 1.0, 4.0]),
    ) {
        derivatives_by_dual::<SO3<D1>, 4>(&raw, u, inv_h);
        derivatives_by_dual::<SO3<D1>, 6>(&raw, u, inv_h);
        derivatives_by_dual::<SE3<D1>, 5>(&raw, u, inv_h);
        derivatives_by_dual::<Product<SO3<D1>, Rn<D1, 3>>, 4>(&raw, u, inv_h);
    }
}

// ---- continuity, products, Rⁿ ------------------------------------------------------------------

fn continuity<G: LieGroup<f64>, const K: usize>(raw: &[f64])
where
    Order<K>: ValidOrder,
{
    // `K + 1` controls: segment A is controls 0..K, segment B is controls 1..K+1.
    let t = taus(raw, G::Tangent::DOF, K + 1, 0.6);
    let mut all = Vec::new();
    let mut x = G::identity();
    for tau in &t {
        x = x * G::exp(&G::Tangent::read_dense(tau));
        all.push(x);
    }
    let a: [G; K] = core::array::from_fn(|i| all[i]);
    let b: [G; K] = core::array::from_fn(|i| all[i + 1]);
    let ka = Segment::new(&a).kinematics(&Weights::<f64, K, 2>::at(1.0, 2.0));
    let kb = Segment::new(&b).kinematics(&Weights::<f64, K, 2>::at(0.0, 2.0));
    assert!(group_err(&ka.value, &kb.value) <= 1e3, "C0 at K={K}");
    if K >= 3 {
        assert!(tangent_err(&ka.vel, &kb.vel) <= 1e4, "C1 at K={K}");
    }
    if K >= 4 {
        assert!(tangent_err(&ka.acc, &kb.acc) <= 1e5, "C2 at K={K}");
    }
}

proptest! {
    #[test]
    fn segments_meet_to_order_k_minus_2(raw in prop::collection::vec(-1.0_f64..1.0, 42)) {
        continuity::<SO3<f64>, 2>(&raw);
        continuity::<SO3<f64>, 3>(&raw);
        continuity::<SO3<f64>, 4>(&raw);
        continuity::<SE3<f64>, 6>(&raw);
        continuity::<R3<f64>, 5>(&raw);
    }

    #[test]
    fn a_product_equals_its_components_bitwise(
        raw in prop::collection::vec(-1.0_f64..1.0, 36),
        u in -0.25_f64..1.25,
    ) {
        const K: usize = 5;
        let t = taus(&raw, 6, K, 0.8);
        let rot: Vec<Vec<f64>> = t.iter().map(|v| v[..3].to_vec()).collect();
        let pos: Vec<Vec<f64>> = t.iter().map(|v| v[3..].to_vec()).collect();
        let r = window::<f64, SO3<f64>, K>(&rot);
        let p = window::<f64, R3<f64>, K>(&pos);
        let both: [Split<f64>; K] = core::array::from_fn(|i| Product::from_parts(r[i], p[i]));
        let w = Weights::<f64, K, 2>::at(u, 3.0);
        let kp = Segment::new(&both).kinematics(&w);
        let kr = Segment::new(&r).kinematics(&w);
        let kt = Segment::new(&p).kinematics(&w);
        let bits = |v: Vec<f64>| v.into_iter().map(f64::to_bits).collect::<Vec<_>>();
        let mut want = dense(&kr.vel);
        want.extend(dense(&kt.vel));
        prop_assert_eq!(bits(dense(&kp.vel)), bits(want));
        let mut want = dense(&kr.acc);
        want.extend(dense(&kt.acc));
        prop_assert_eq!(bits(dense(&kp.acc)), bits(want));
        let (vr, vp) = kp.value.parts();
        prop_assert_eq!(bits(dense(&vr.log())), bits(dense(&kr.value.log())));
        prop_assert_eq!(bits(vp.0 .0.to_vec()), bits(kt.value.0 .0.to_vec()));
    }

    #[test]
    fn rn_is_the_non_cumulative_sum(
        raw in prop::collection::vec(-10.0_f64..10.0, 18),
        u in -0.25_f64..1.25,
    ) {
        const K: usize = 6;
        let c: [R3<f64>; K] = core::array::from_fn(|i| Rn(Vector([raw[3 * i], raw[3 * i + 1], raw[3 * i + 2]])));
        let w = Weights::<f64, K, 2>::at(u, 2.0);
        let k = Segment::new(&c).kinematics(&w);
        let b = |l: &[f64; K], m: usize| l[m] - if m + 1 < K { l[m + 1] } else { 0.0 };
        for axis in 0..3 {
            let sum = |l: &[f64; K]| (0..K).map(|m| b(l, m) * c[m].0 .0[axis]).sum::<f64>();
            let scale = c.iter().map(|x| x.0 .0[axis].abs()).fold(1.0_f64, f64::max);
            let tol = 64.0 * f64::EPSILON * scale;
            prop_assert!((k.value.0 .0[axis] - sum(&w.lambda)).abs() <= tol);
            prop_assert!((dense(&k.vel)[axis] - sum(&w.dot)).abs() <= 4.0 * tol);
            prop_assert!((dense(&k.acc)[axis] - sum(&w.ddot)).abs() <= 16.0 * tol);
        }
    }
}

// ---- control Jacobians ---------------------------------------------------------------------------

use crate::reference::spline_jacobians;
use crate::{ControlJacobians, JacSegment};
use helicoid::{Jac, Left, Right, Side};
use helicoid_linalg::StridedMut;

fn dense_jac<T: Tangent<f64>, J: Jac<f64, T>>(j: &J) -> Vec<f64> {
    let n = T::DOF;
    let mut buf = std::vec![0.0; n * n];
    j.write_dense(&mut StridedMut::col_major(&mut buf, n, n));
    buf
}

/// Frobenius `‖a − b‖ / max(‖b‖, 1)` in `f64::EPSILON`.
fn jac_err(a: &[f64], b: &[f64]) -> f64 {
    let diff: f64 = a
        .iter()
        .zip(b)
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt();
    let scale: f64 = b.iter().map(|y| y * y).sum::<f64>().sqrt();
    diff / scale.max(1.0) / f64::EPSILON
}

fn jac_twins<G: LieGroup<f64>, Sd: Side, const K: usize>(raw: &[f64], u: f64, inv_h: f64, step: f64)
where
    Order<K>: ValidOrder,
{
    let w = Weights::<f64, K, 2>::at(u, inv_h);
    let ctrl = window::<f64, G, K>(&taus(raw, G::Tangent::DOF, K, step));
    let seg = JacSegment::new(&ctrl);
    let (kin, fast) = seg.kinematics_jacobians::<Sd, 2>(&w);
    let slow: ControlJacobians<f64, G, K> = spline_jacobians::<f64, G, Sd, K, 2>(&ctrl, &w);
    let plain = seg.segment().kinematics(&w);
    let bits = |t: &G::Tangent| dense(t).into_iter().map(f64::to_bits).collect::<Vec<_>>();
    assert_eq!(
        bits(&kin.vel),
        bits(&plain.vel),
        "kinematics_jacobians repeats kinematics"
    );
    assert_eq!(bits(&kin.acc), bits(&plain.acc));
    for m in 0..K {
        let ev = jac_err(&dense_jac(&fast.value[m]), &dense_jac(&slow.value[m]));
        let ew = jac_err(&dense_jac(&fast.vel[m]), &dense_jac(&slow.vel[m]));
        let ea = jac_err(&dense_jac(&fast.acc[m]), &dense_jac(&slow.acc[m]));
        assert!(
            ev <= 1e3 && ew <= 1e4 && ea <= 1e5,
            "K={K} m={m} right={}: value {ev} vel {ew} acc {ea}",
            Sd::IS_RIGHT
        );
    }
}

proptest! {
    #[test]
    fn spline_jacobians_matches_reference(
        raw in prop::collection::vec(-1.0_f64..1.0, 36),
        u in -0.25_f64..1.25,
        inv_h in prop::sample::select(std::vec![1e-3_f64, 1.0, 1e3]),
        step in prop::sample::select(std::vec![1e-9_f64, 1e-3, 0.5, 1.5]),
    ) {
        jac_twins::<SO3<f64>, Right, 4>(&raw, u, inv_h, step);
        jac_twins::<SO3<f64>, Left, 6>(&raw, u, inv_h, step);
        jac_twins::<SE3<f64>, Right, 5>(&raw, u, inv_h, step);
        jac_twins::<SE3<f64>, Left, 4>(&raw, u, inv_h, step);
        jac_twins::<R3<f64>, Right, 3>(&raw, u, inv_h, step);
        jac_twins::<Split<f64>, Right, 6>(&raw, u, inv_h, step);
        jac_twins::<SO3<f64>, Left, 2>(&raw, u, inv_h, step);
        jac_twins::<Split<f64>, Left, 3>(&raw, u, inv_h, step);
    }
}

/// Below `D = 2`, the formed rows are the `D = 2` rows to the bit and the rest are zero.
fn jacobian_rows_follow_d<G: LieGroup<f64>, Sd: Side, const K: usize>(raw: &[f64], u: f64)
where
    Order<K>: ValidOrder,
{
    let ctrl = window::<f64, G, K>(&taus(raw, G::Tangent::DOF, K, 0.7));
    let seg = JacSegment::new(&ctrl);
    let full = seg.kinematics_jacobians::<Sd, 2>(&Weights::at(u, 2.0)).1;
    let d1 = seg.kinematics_jacobians::<Sd, 1>(&Weights::at(u, 2.0)).1;
    let d0 = seg.kinematics_jacobians::<Sd, 0>(&Weights::at(u, 2.0)).1;
    let bits = |j: &G::Jac| {
        dense_jac(j)
            .into_iter()
            .map(f64::to_bits)
            .collect::<Vec<_>>()
    };
    let is_zero = |j: &G::Jac| dense_jac(j).iter().all(|&x| x == 0.0);
    for m in 0..K {
        assert_eq!(bits(&d1.value[m]), bits(&full.value[m]));
        assert_eq!(bits(&d1.vel[m]), bits(&full.vel[m]));
        assert!(is_zero(&d1.acc[m]));
        assert_eq!(bits(&d0.value[m]), bits(&full.value[m]));
        assert!(is_zero(&d0.vel[m]) && is_zero(&d0.acc[m]));
    }
}

proptest! {
    #[test]
    fn jacobian_rows_above_d_are_zero(
        raw in prop::collection::vec(-1.0_f64..1.0, 36),
        u in 0.0_f64..1.0,
    ) {
        jacobian_rows_follow_d::<SO3<f64>, Right, 4>(&raw, u);
        jacobian_rows_follow_d::<SE3<f64>, Left, 5>(&raw, u);
    }
}

/// `Dual` along `X_m Exp(ε)` (right) or `Exp(ε) X_m` (left), `N` = the group's DOF lanes.
fn jacobians_by_dual<G, GD, Sd: Side, const K: usize, const N: usize>(
    raw: &[f64],
    u: f64,
    inv_h: f64,
) where
    G: LieGroup<f64>,
    GD: LieGroup<Dual<f64, N>>,
    Order<K>: ValidOrder,
{
    type DN<const N: usize> = Dual<f64, N>;
    let t = taus(raw, N, K, 0.7);
    let ctrl = window::<f64, G, K>(&t);
    let (kin, jac) =
        JacSegment::new(&ctrl).kinematics_jacobians::<Sd, 2>(&Weights::<f64, K, 2>::at(u, inv_h));
    let lifted: Vec<Vec<DN<N>>> = t
        .iter()
        .map(|v| v.iter().map(|&x| DN::<N>::constant(x)).collect())
        .collect();
    let base = window::<DN<N>, GD, K>(&lifted);
    let eps = GD::Tangent::read_dense(&core::array::from_fn::<DN<N>, N, _>(|i| {
        DN::<N>::variable(0.0, i)
    }));
    let x0 = {
        let v: Vec<DN<N>> = dense(&kin.value.log())
            .into_iter()
            .map(DN::<N>::constant)
            .collect();
        GD::exp(&GD::Tangent::read_dense(&v))
    };
    let wd = Weights::<DN<N>, K, 2>::at(DN::<N>::constant(u), DN::<N>::constant(inv_h));
    for m in 0..K {
        let mut moved = base;
        moved[m] = if Sd::IS_RIGHT {
            base[m] * GD::exp(&eps)
        } else {
            GD::exp(&eps) * base[m]
        };
        let k = Segment::new(&moved).kinematics(&wd);
        let local = if Sd::IS_RIGHT {
            (x0.inverse() * k.value).log()
        } else {
            (k.value * x0.inverse()).log()
        };
        // Column `c` of a Jacobian is lane `c` of each output component; dense is column-major.
        let cols = |t: &GD::Tangent| {
            let rows = dense(t);
            let mut out = std::vec![0.0; N * N];
            for (r, x) in rows.iter().enumerate() {
                for c in 0..N {
                    out[c * N + r] = x.d[c];
                }
            }
            out
        };
        let ev = jac_err(&dense_jac(&jac.value[m]), &cols(&local));
        let ew = jac_err(&dense_jac(&jac.vel[m]), &cols(&k.vel));
        let ea = jac_err(&dense_jac(&jac.acc[m]), &cols(&k.acc));
        assert!(
            ev <= 1e4 && ew <= 1e4 && ea <= 1e5,
            "m={m} right={}: {ev} {ew} {ea}",
            Sd::IS_RIGHT
        );
    }
}

proptest! {
    #[test]
    fn jacobians_match_dual(
        raw in prop::collection::vec(-1.0_f64..1.0, 36),
        u in 0.0_f64..1.0,
        inv_h in prop::sample::select(std::vec![0.5_f64, 1.0, 4.0]),
    ) {
        jacobians_by_dual::<SO3<f64>, SO3<Dual<f64, 3>>, Right, 4, 3>(&raw, u, inv_h);
        jacobians_by_dual::<SO3<f64>, SO3<Dual<f64, 3>>, Left, 6, 3>(&raw, u, inv_h);
        jacobians_by_dual::<SE3<f64>, SE3<Dual<f64, 6>>, Right, 5, 6>(&raw, u, inv_h);
        jacobians_by_dual::<SE3<f64>, SE3<Dual<f64, 6>>, Left, 4, 6>(&raw, u, inv_h);
        jacobians_by_dual::<R3<f64>, Rn<Dual<f64, 3>, 3>, Right, 5, 3>(&raw, u, inv_h);
    }
}
