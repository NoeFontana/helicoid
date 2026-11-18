//! `so3::{exp_gamma_apply, exp_gamma_apply_jacobians}` against their definitions, the separate
//! calls, **to the bit** (`0119`, `NUMERICS.md` §14): `SO3::exp`, `gamma1(φ) * v`, `gamma2(φ) * v`
//! and `gamma_apply_jacobians::<1 | 2>`, at `f64`, `f32`, `Dual<f64, 1>` and `Dual<f32, 1>`, every
//! lane compared.
//!
//! The draws are `gamma_tests`' (`θ` log-uniform on `[1e-8, π]`), `θ²` at every switch of the
//! members and its two neighbours along two axes, and `θ` across the band where a member's arm
//! depends on the group it reproduces, the band a single group would part from them in.

use std::vec::Vec;

use helicoid_linalg::{Dual, Mat3, Precision, Real, Vec3, Vector};

use crate::coeffs::exp_gamma_switches;
use crate::so3::{
    exp_gamma_apply, exp_gamma_apply_jacobians, gamma1, gamma2, gamma_apply_jacobians,
};
use crate::{LieGroup, SO3Tangent, SO3};

/// Every lane's bits: one for a float, the value's and then the derivative's for `Dual`.
trait Lanes: Real {
    fn lanes(self, out: &mut Vec<u64>);
}

impl Lanes for f64 {
    fn lanes(self, out: &mut Vec<u64>) {
        out.push(self.to_bits());
    }
}

impl Lanes for f32 {
    fn lanes(self, out: &mut Vec<u64>) {
        out.push(u64::from(self.to_bits()));
    }
}

impl<S: Lanes> Lanes for Dual<S, 1> {
    fn lanes(self, out: &mut Vec<u64>) {
        self.v.lanes(out);
        self.d[0].lanes(out);
    }
}

fn of_vec<S: Lanes>(v: Vec3<S>, out: &mut Vec<u64>) {
    v.0.iter().for_each(|x| x.lanes(out));
}

fn of_mat<S: Lanes>(m: &Mat3<S>, out: &mut Vec<u64>) {
    (0..9).for_each(|k| m.get(k / 3, k % 3).lanes(out));
}

fn of_rot<S: Lanes>(r: &SO3<S>, out: &mut Vec<u64>) {
    let q = r.quat();
    [q.w, q.x, q.y, q.z].iter().for_each(|x| x.lanes(out));
}

fn of_triple<S: Lanes>((gv, j, g): &(Vec3<S>, Mat3<S>, Mat3<S>), out: &mut Vec<u64>) {
    of_vec(*gv, out);
    of_mat(j, out);
    of_mat(g, out);
}

/// Both fused functions against the separate calls at one `(φ, v)`.
fn fused_is_separate<S: Lanes>(phi: &SO3Tangent<S>, v: Vec3<S>) {
    let (mut got, mut want) = (Vec::new(), Vec::new());
    let (r, g1v, g2v) = exp_gamma_apply(phi, v);
    of_rot(&r, &mut got);
    of_vec(g1v, &mut got);
    of_vec(g2v, &mut got);
    let r = SO3::exp(phi);
    of_rot(&r, &mut want);
    of_vec(gamma1(phi) * v, &mut want);
    of_vec(gamma2(phi) * v, &mut want);
    let at = phi.phi.0.map(|x| x.value_f64());
    assert_eq!(got, want, "exp_gamma_apply at φ = {at:?}");

    let (mut got, mut want) = (Vec::new(), Vec::new());
    let (r, t1, t2) = exp_gamma_apply_jacobians(phi, v);
    of_rot(&r, &mut got);
    of_triple(&t1, &mut got);
    of_triple(&t2, &mut got);
    of_rot(&SO3::exp(phi), &mut want);
    of_triple(&gamma_apply_jacobians::<1, S>(phi, v), &mut want);
    of_triple(&gamma_apply_jacobians::<2, S>(phi, v), &mut want);
    assert_eq!(got, want, "exp_gamma_apply_jacobians at φ = {at:?}");
}

/// A scalar of `S`'s precision: `x` rounded to it, and for `Dual` the derivative lane `dx`.
trait Seed: Lanes {
    fn seed(x: f64, dx: f64) -> Self;
}

fn rounded<S: Real>(x: f64) -> S {
    match S::PRECISION {
        Precision::F64 => S::lit(x),
        Precision::F32 => S::lit(f64::from(x as f32)),
    }
}

impl Seed for f64 {
    fn seed(x: f64, _: f64) -> Self {
        x
    }
}

impl Seed for f32 {
    fn seed(x: f64, _: f64) -> Self {
        x as f32
    }
}

impl<S: Seed> Seed for Dual<S, 1> {
    fn seed(x: f64, dx: f64) -> Self {
        Dual {
            v: S::seed(x, 0.0),
            d: [rounded::<S>(dx)],
        }
    }
}

/// `x` the next float of `S`'s precision up (`dir > 0`) or down.
fn step<S: Real>(x: f64, dir: i32) -> f64 {
    match (S::PRECISION, dir > 0) {
        (Precision::F64, true) => x.next_up(),
        (Precision::F64, false) => x.next_down(),
        (Precision::F32, true) => f64::from((x as f32).next_up()),
        (Precision::F32, false) => f64::from((x as f32).next_down()),
    }
}

/// `(φ, v)` from `f64` draws: `φ`'s derivative lane along a fixed direction, `v` constant.
fn point<S: Seed>(phi: [f64; 3], v: [f64; 3]) -> (SO3Tangent<S>, Vec3<S>) {
    const DIR: [f64; 3] = [0.6, -0.48, 0.64];
    let phi = SO3Tangent {
        phi: Vector(core::array::from_fn(|i| S::seed(phi[i], DIR[i]))),
    };
    (phi, Vector(v.map(|x| S::seed(x, 0.0))))
}

const VS: [[f64; 3]; 2] = [[1.0, -2.0, 0.5], [-0.3, 0.7, 0.9]];
const AXES: [[f64; 3]; 2] = [[1.0, 0.0, 0.0], [0.48, 0.6, -0.64]];

fn draws(count: usize) -> Vec<[f64; 3]> {
    let mut st = 0x2545_F491_4F6C_DD1D_u64;
    let mut next = || {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        (st >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut out = std::vec![[0.0; 3], [5e-324, -5e-324, 1e-310]];
    for _ in 0..count {
        let theta = 10f64.powf(-8.0 + (8.0 + core::f64::consts::PI.log10()) * next());
        let (u, w) = (2.0 * next() - 1.0, 2.0 * core::f64::consts::PI * next());
        let r = (1.0 - u * u).sqrt();
        out.push([r * w.cos(), r * w.sin(), u].map(|a| a * theta));
    }
    for theta in [
        3.0,
        3.1,
        3.13,
        core::f64::consts::PI - 1e-6,
        core::f64::consts::PI,
    ] {
        out.extend(AXES.map(|a| a.map(|x| x * theta)));
    }
    out
}

fn fused_is_separate_at<S: Seed>(stride: usize) {
    let check = |phi: [f64; 3]| {
        for v in VS {
            let (p, w) = point::<S>(phi, v);
            fused_is_separate(&p, w);
        }
    };
    draws(4_000).into_iter().for_each(check);
    let (switches, (lo, hi)) = exp_gamma_switches::<S>();
    for theta in switches.iter().map(|z| z.sqrt()) {
        for axis in AXES {
            let (mut up, mut down) = (theta, theta);
            for _ in 0..8 {
                check(axis.map(|a| a * up));
                check(axis.map(|a| a * down));
                up = step::<S>(up, 1);
                down = step::<S>(down, -1);
            }
        }
    }
    // The band along the first axis, `θ² = θ̂²` there up to one rounding: every `stride`-th
    // binary32, or 4096 binary64 in geometric steps.
    let (mut theta, top) = (lo.sqrt(), hi.sqrt());
    let r = (top / theta).powf(1.0 / 4096.0);
    while theta < top {
        check([theta, 0.0, 0.0]);
        theta = match S::PRECISION {
            Precision::F64 => theta * r,
            Precision::F32 => (0..stride).fold(theta, |t, _| step::<S>(t, 1)),
        };
    }
}

#[test]
fn exp_gamma_apply_matches_reference_f64() {
    fused_is_separate_at::<f64>(1);
}

#[test]
fn exp_gamma_apply_matches_reference_f32() {
    fused_is_separate_at::<f32>(97);
}

#[test]
fn exp_gamma_apply_matches_reference_dual64() {
    fused_is_separate_at::<Dual<f64, 1>>(1);
}

#[test]
fn exp_gamma_apply_matches_reference_dual32() {
    fused_is_separate_at::<Dual<f32, 1>>(389);
}
