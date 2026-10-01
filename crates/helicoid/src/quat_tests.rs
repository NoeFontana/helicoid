//! `Quat`: Hamilton identities, the converters, `R(q)`, renormalization and the asserted domains.
//!
//! A bound is recorded: twice the worst error, in units of `u`, of `measure_worst_errors` (a
//! seeded splitmix64 stream, 10^6 cases per scalar; `cargo nextest run -p helicoid --run-ignored
//! only measure_worst_errors --no-capture`), rounded up to an integer. `Dual<f64, 4>` reuses the
//! `f64` bounds: its value path is bit-identical (`dual_value_is_plain_value`). The worst
//! errors, `f64` then `f32`: `conj` 6.00 and 3.00, `assoc` 5.00 and 5.00,
//! `orth` 13.00 and 6.90, `det` 18.00 and 6.94, `matrix` 7.00 and 3.35, `sandwich` 6.00 and 5.48,
//! `jpl` 7.00 and 3.37, `step` (one Newton step from `|η| = 2^-40` or `2^-16`, the edge of what
//! `from_wxyz_unchecked` accepts) 5.00 and 3.73. `orth` and `det` include the rounding of the
//! `f64` check itself.

// `u` is `laws::unit`, the unit roundoff of `NUMERICS.md` §2.1; it is defined there, not again here.
use crate::laws::{unit as u, Sample};
use crate::Quat;
use core::array;
use helicoid_linalg::{Blend, Dual, Precision, Real};
use proptest::prelude::*;
use std::vec::Vec;

struct Bounds {
    conj: f64,
    assoc: f64,
    orth: f64,
    det: f64,
    matrix: f64,
    sandwich: f64,
    jpl: f64,
    step: f64,
}
const F64: Bounds = Bounds {
    conj: 12.0,
    assoc: 10.0,
    orth: 26.0,
    det: 36.0,
    matrix: 14.0,
    sandwich: 12.0,
    jpl: 14.0,
    step: 10.0,
};
const F32: Bounds = Bounds {
    conj: 6.0,
    assoc: 10.0,
    orth: 14.0,
    det: 14.0,
    matrix: 7.0,
    sandwich: 11.0,
    jpl: 7.0,
    step: 8.0,
};

/// A struct literal: no domain, for the quaternions that are not unit.
fn raw<S: Real>(w: S, x: S, y: S, z: S) -> Quat<S> {
    Quat { w, x, y, z }
}

/// The quaternion whose components are the samples `v` (`Dual`: one variable per component).
fn quat<S: Sample>(v: &[f64; 4]) -> Quat<S> {
    let s: [S; 4] = array::from_fn(|i| S::sample(v[i], i));
    Quat::from_wxyz_unchecked(s[0], s[1], s[2], s[3])
}

fn vals<S: Real>(q: &Quat<S>) -> [f64; 4] {
    [q.w, q.x, q.y, q.z].map(S::value_f64)
}

/// The largest componentwise difference, in `u`.
fn diff<S: Real>(a: &[f64], b: &[f64]) -> f64 {
    let d = a.iter().zip(b).map(|(x, y)| (x - y).abs());
    d.fold(0.0, |m, d| if d > m || d.is_nan() { d } else { m }) / u::<S>()
}

fn normalized(v: [f64; 4]) -> [f64; 4] {
    let n = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    v.map(|x| x / n)
}

/// `m 2^e`, `|m| < 1`, `-6 <= e <= 0`, four of them, with `‖.‖ >= 1/16`, then normalized.
fn unit_sample() -> impl Strategy<Value = [f64; 4]> {
    let entry = (-1.0_f64..1.0, -6_i32..=0).prop_map(|(m, e)| m * 2_f64.powi(e));
    proptest::array::uniform4(entry)
        .prop_filter("away from zero", |v| {
            v.iter().map(|x| x * x).sum::<f64>() > 0.004
        })
        .prop_map(normalized)
}

fn vec3() -> impl Strategy<Value = [f64; 3]> {
    proptest::array::uniform3(-1.0_f64..1.0)
}

/// `R` as rows, from the entries of `to_matrix`.
fn rows<S: Real>(q: &Quat<S>) -> [[f64; 3]; 3] {
    let m = q.to_matrix();
    array::from_fn(|r| array::from_fn(|c| m.get(r, c).value_f64()))
}

fn matmul(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    array::from_fn(|r| array::from_fn(|c| (0..3).map(|k| a[r][k] * b[k][c]).sum()))
}

fn transpose(a: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    array::from_fn(|r| array::from_fn(|c| a[c][r]))
}

const EYE: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// `‖q q* - 1‖` and `‖q* q - 1‖`.
fn conj_is_inverse<S: Sample>(v: &[f64; 4]) -> f64 {
    let q = quat::<S>(v);
    let one = [1.0, 0.0, 0.0, 0.0];
    diff::<S>(&vals(&(q * q.conjugate())), &one).max(diff::<S>(&vals(&(q.conjugate() * q)), &one))
}

/// `(p q) r - p (q r)`.
fn associativity<S: Sample>(a: &[f64; 4], b: &[f64; 4], c: &[f64; 4]) -> f64 {
    let (p, q, r) = (quat::<S>(a), quat::<S>(b), quat::<S>(c));
    diff::<S>(&vals(&((p * q) * r)), &vals(&(p * (q * r))))
}

/// `‖RᵀR - I‖` and `|det R - 1|`, computed in `f64` from the entries.
fn orthogonal<S: Sample>(v: &[f64; 4]) -> (f64, f64) {
    let r = rows(&quat::<S>(v));
    let g = matmul(&transpose(&r), &r);
    let det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
    (
        diff::<S>(g.as_flattened(), EYE.as_flattened()),
        (det - 1.0).abs() / u::<S>(),
    )
}

/// `I + 2w [u]× + 2 [u]×²` (unit `q`, `docs/maths/so3.md` SO.2) from the entries of `q`.
fn rodrigues(q: [f64; 4]) -> [[f64; 3]; 3] {
    let [w, x, y, z] = q;
    let k = [[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]];
    let k2 = matmul(&k, &k);
    array::from_fn(|r| array::from_fn(|c| EYE[r][c] + 2.0 * w * k[r][c] + 2.0 * k2[r][c]))
}

/// `R(q)` against `rodrigues`.
fn matrix_form<S: Sample>(v: &[f64; 4]) -> f64 {
    let q = quat::<S>(v);
    diff::<S>(rows(&q).as_flattened(), rodrigues(vals(&q)).as_flattened())
}

/// `R(q) p` against the definition `q (0, p) q*`, in `S`.
fn sandwich<S: Sample>(v: &[f64; 4], p: &[f64; 3]) -> f64 {
    let q = quat::<S>(v);
    let z = S::zero();
    // The point is a constant, not a sample: `quat` has already taken `Dual`'s four variable lanes
    // for the components of `q`, so sampling it would alias three of them. The law compares values
    // only, so it needs no lane of its own — but the collision would be silent.
    let pq = raw(z, S::constant(p[0]), S::constant(p[1]), S::constant(p[2]));
    let (r, p) = (rows(&q), vals(&pq));
    let out = vals(&(q * pq * q.conjugate()));
    let rp: [f64; 3] = array::from_fn(|i| (0..3).map(|k| r[i][k] * p[k + 1]).sum());
    diff::<S>(&rp, &out[1..]).max(out[0].abs() / u::<S>())
}

/// `to_matrix(from_jpl([x, y, z, w]))` against `(2w² - 1) I - 2w [u]× + 2 u uᵀ`, the matrix `C(q)`
/// of a JPL quaternion (Trawny and Roumeliotis 2005; Sommer et al. 2018). `docs/maths` does not
/// state it; in mpmath it equals `expm(-θ [k]×)` for `θ = 2 atan2(‖u‖, w)`, `k = u/‖u‖` (100 random
/// unit quaternions, to 1e-45 u), and the matrix of `from_jpl` equals `q (0, e_c) q*` of `(w, -u)`
/// (200 random quaternions, few u).
fn jpl_matrix<S: Sample>(v: &[f64; 4]) -> f64 {
    let jpl: [S; 4] = array::from_fn(|i| S::sample(v[i], i));
    let q = Quat::from_jpl(jpl);
    let [x, y, z, w] = jpl.map(S::value_f64);
    let k = [[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]];
    let c: [[f64; 3]; 3] = array::from_fn(|r| {
        array::from_fn(|c| {
            let uu = [x, y, z][r] * [x, y, z][c];
            let d = if r == c { 2.0 * w * w - 1.0 } else { 0.0 };
            d - 2.0 * w * k[r][c] + 2.0 * uu
        })
    });
    diff::<S>(rows(&q).as_flattened(), c.as_flattened())
}

/// `|η'|` after one step from a unit `q` scaled to `η = ±2^-40` (`f64`) or `±2^-16` (`f32`), the
/// edge of what `from_wxyz_unchecked` accepts.
fn newton_step<S: Sample>(v: &[f64; 4], negative: bool) -> f64 {
    let edge = match S::PRECISION {
        Precision::F64 => 2_f64.powi(-40),
        Precision::F32 => 2_f64.powi(-16),
    };
    let lambda = (1.0 + if negative { -edge } else { edge }).sqrt();
    let s: [S; 4] = array::from_fn(|i| S::sample(v[i] * lambda, i));
    let mut q = raw(s[0], s[1], s[2], s[3]);
    q.renormalize();
    eta(&q).abs() / u::<S>()
}

/// `‖q‖² - 1`, summed in `f64` from the entries.
fn eta<S: Real>(q: &Quat<S>) -> f64 {
    vals(q).iter().map(|x| x * x).sum::<f64>() - 1.0
}

macro_rules! props {
    ($m:ident, $S:ty, $B:ident) => {
        mod $m {
            use super::*;
            type S = $S;
            fn within(v: f64, bound: f64) -> Result<(), TestCaseError> {
                prop_assert!(v <= bound, "worst error {} u exceeds {} u", v, bound);
                Ok(())
            }
            proptest! {
                #[test]
                fn conjugate_is_the_inverse_of_a_unit(a in unit_sample()) {
                    within(conj_is_inverse::<S>(&a), $B.conj)?;
                }
                #[test]
                fn hamilton_associativity(a in unit_sample(), b in unit_sample(), c in unit_sample()) {
                    within(associativity::<S>(&a, &b, &c), $B.assoc)?;
                }
                #[test]
                fn matrix_is_orthogonal_with_det_one(a in unit_sample()) {
                    let (o, d) = orthogonal::<S>(&a);
                    within(o, $B.orth)?;
                    within(d, $B.det)?;
                }
                #[test]
                fn matrix_matches_the_rodrigues_form(a in unit_sample()) {
                    within(matrix_form::<S>(&a), $B.matrix)?;
                }
                #[test]
                fn matrix_matches_the_sandwich(a in unit_sample(), p in vec3()) {
                    within(sandwich::<S>(&a, &p), $B.sandwich)?;
                }
                #[test]
                fn from_jpl_matches_matrix(a in unit_sample()) {
                    within(jpl_matrix::<S>(&a), $B.jpl)?;
                }
                #[test]
                fn a_newton_step_repairs_the_accepted_range(a in unit_sample(), neg in any::<bool>()) {
                    within(newton_step::<S>(&a, neg), $B.step)?;
                }
            }
        }
    };
}

props!(as_f64, f64, F64);
props!(as_f32, f32, F32);
props!(as_dual, Dual<f64, 4>, F64);

/// splitmix64, seeded; `unif` is uniform on `[-1, 1)`.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn unif(&mut self) -> f64 {
        (self.next() >> 11) as f64 / 2_f64.powi(52) - 1.0
    }
    /// The distribution of `unit_sample`.
    ///
    /// The draws are an explicit loop, not `array::from_fn`: two per component come off the stream
    /// and `from_fn` does not promise the order it calls its closure in, so the recorded bounds of
    /// the module header would be reproducible only by accident.
    fn unit(&mut self) -> [f64; 4] {
        loop {
            let mut v = [0.0; 4];
            for e in &mut v {
                let m = self.unif();
                *e = m * 2_f64.powi(-((self.next() % 7) as i32));
            }
            if v.iter().map(|x| x * x).sum::<f64>() > 0.004 {
                return normalized(v);
            }
        }
    }
}

/// The worst error of every law over 10^6 seeded cases; prints them in the order of the header.
#[test]
#[ignore = "measurement: prints the figures the bounds are recorded from"]
#[allow(clippy::print_stdout)]
fn measure_worst_errors() {
    fn run<S: Sample>() -> [f64; 8] {
        let mut rng = Rng(0x0123_4567_89AB_CDEF);
        let mut w = [0.0_f64; 8];
        let mut up = |i: usize, v: f64| w[i] = if v.is_nan() || v > w[i] { v } else { w[i] };
        for _ in 0..1_000_000 {
            let (a, b, c) = (rng.unit(), rng.unit(), rng.unit());
            let p = [rng.unif(), rng.unif(), rng.unif()];
            up(0, conj_is_inverse::<S>(&a));
            up(1, associativity::<S>(&a, &b, &c));
            let (o, d) = orthogonal::<S>(&a);
            up(2, o);
            up(3, d);
            up(4, matrix_form::<S>(&a));
            up(5, sandwich::<S>(&a, &p));
            up(6, jpl_matrix::<S>(&a));
            up(7, newton_step::<S>(&a, rng.next() & 1 == 1));
        }
        w
    }
    std::println!("f64  {:.2?}", run::<f64>());
    std::println!("f32  {:.2?}", run::<f32>());
}

fn fb<const N: usize>(a: [f64; N]) -> [u64; N] {
    a.map(f64::to_bits)
}

fn bits<S: Real>(q: &Quat<S>) -> [u64; 4] {
    fb(vals(q))
}

fn mbits(r: [[f64; 3]; 3]) -> [[u64; 3]; 3] {
    r.map(fb)
}

#[test]
fn hamilton_identities() {
    let (one, i, j, k) = (
        raw(1.0_f64, 0.0, 0.0, 0.0),
        raw(0.0, 1.0, 0.0, 0.0),
        raw(0.0, 0.0, 1.0, 0.0),
        raw(0.0, 0.0, 0.0, 1.0),
    );
    let minus = |q: Quat<f64>| raw(-q.w, -q.x, -q.y, -q.z);
    // `-0 + 0 = +0`: the identities are about values, and `-k` has zeros of either sign.
    let canon = |q: Quat<f64>| vals(&q).map(|x| (x + 0.0).to_bits());
    let same = |a: Quat<f64>, b: Quat<f64>| assert_eq!(canon(a), canon(b));
    // Hamilton: `i² = j² = k² = ijk = -1`, `ij = k`, `ji = -k` and their cycles.
    for sq in [i * i, j * j, k * k, i * j * k] {
        same(sq, minus(one));
    }
    same(i * j, k);
    same(j * k, i);
    same(k * i, j);
    same(j * i, minus(k));
    same(k * j, minus(i));
    same(i * k, minus(j));
    same(one * k, k);
    // By hand: (1 + 2i + 3j + 4k)(5 + 6i + 7j + 8k), and the other order.
    let (p, q) = (raw(1.0_f64, 2.0, 3.0, 4.0), raw(5.0_f64, 6.0, 7.0, 8.0));
    assert_eq!(bits(&(p * q)), fb([-60.0, 12.0, 30.0, 24.0]));
    assert_eq!(bits(&(q * p)), fb([-60.0, 20.0, 14.0, 32.0]));
    assert_eq!(p.norm_sq().to_bits(), 30.0_f64.to_bits());
    same(p.conjugate(), raw(1.0, -2.0, -3.0, -4.0));
}

#[test]
fn matrix_hand_cases() {
    // `(1/2)(1, 1, 1, 1)` is the cyclic permutation `e_x -> e_y -> e_z -> e_x`, exactly.
    let cyc = Quat::from_wxyz_unchecked(0.5_f64, 0.5, 0.5, 0.5);
    let r = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    assert_eq!(mbits(rows(&cyc)), mbits(r));
    // Active: `R e_x` is `q i q*`, the first column of the rows above.
    let i = raw(0.0, 1.0, 0.0, 0.0);
    assert_eq!(bits(&(cyc * i * cyc.conjugate())), fb([0.0, 0.0, 1.0, 0.0]));
    let diag = |a: f64, b: f64, c: f64| [[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]];
    let half = |x, y, z| rows(&Quat::from_wxyz_unchecked(0.0_f64, x, y, z));
    assert_eq!(mbits(half(1.0, 0.0, 0.0)), mbits(diag(1.0, -1.0, -1.0)));
    assert_eq!(mbits(half(0.0, 1.0, 0.0)), mbits(diag(-1.0, 1.0, -1.0)));
    assert_eq!(mbits(half(0.0, 0.0, 1.0)), mbits(diag(-1.0, -1.0, 1.0)));
    // Not unit: `R(q) = ‖q‖² R(q / ‖q‖)`.
    assert_eq!(
        mbits(rows(&raw(2.0_f64, 0.0, 0.0, 0.0))),
        mbits(diag(4.0, 4.0, 4.0))
    );
    assert_eq!(
        mbits(rows(&raw(0.0_f64, 2.0, 0.0, 0.0))),
        mbits(diag(4.0, -4.0, -4.0))
    );
}

#[test]
fn jpl_hand_case() {
    // The JPL matrix of `[1/2, 1/2, 1/2, 1/2]` is the transpose of the Hamilton one of the same
    // numbers: `C(q) = R(q)ᵀ`, `e_y -> e_x`.
    let a = [0.5_f64; 4];
    let (h, j) = (Quat::from_xyzw(a), Quat::from_jpl(a));
    assert_eq!(mbits(rows(&j)), mbits(transpose(&rows(&h))));
    assert_eq!(
        mbits(rows(&j)),
        mbits([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
    );
}

proptest! {
    #[test]
    fn quat_converters_round_trip(bits in proptest::array::uniform4(any::<u64>())) {
        let a = bits.map(f64::from_bits);
        let same = |x: [f64; 4], y: [f64; 4]| x.map(f64::to_bits) == y.map(f64::to_bits);
        let q = Quat::from_xyzw(a);
        // Named fields, in the order `[x, y, z, w]`.
        prop_assert!(same(vals(&q), [a[3], a[0], a[1], a[2]]));
        prop_assert!(same(q.to_xyzw(), a));
        // `from_jpl` conjugates: `(w, -x, -y, -z)`, and is undone by conjugating again.
        let j = Quat::from_jpl(a);
        prop_assert!(same(j.to_xyzw(), [-a[0], -a[1], -a[2], a[3]]));
        prop_assert!(same(vals(&j), vals(&q.conjugate())));
        prop_assert!(same(Quat::from_jpl(j.to_xyzw()).to_xyzw(), a));
    }
}

/// `η_0, η_1, ...` of repeated steps from `η_0` (the scale of a unit quaternion).
fn iterates<S: Sample>(eta0: f64, steps: usize) -> Vec<f64> {
    let s = 0.5 * (1.0 + eta0).sqrt();
    let mut q = raw(
        S::sample(s, 0),
        S::sample(s, 1),
        S::sample(s, 2),
        S::sample(s, 3),
    );
    let mut out = Vec::new();
    for _ in 0..=steps {
        out.push(eta(&q));
        q.renormalize();
    }
    out
}

/// `η' = -3/4 η² + 1/4 η³` exactly (`docs/maths/so3.md` SO.14): each step squares the deviation
/// with constant `-3/4`, checked while the truncation term is above rounding.
fn quadratic_convergence<S: Sample>(eta0: f64, live: usize, tol: f64) {
    let e = iterates::<S>(eta0, live + 1);
    for k in 0..live {
        let predicted = -0.75 * e[k] * e[k] + 0.25 * e[k].powi(3);
        let off = (e[k + 1] - predicted).abs() / u::<S>();
        assert!(
            off <= tol,
            "step {k}: {off} u off the closed form (eta = {e:?})"
        );
        assert!(
            (e[k + 1] / (e[k] * e[k]) + 0.75).abs() < 0.01,
            "step {k}: not quadratic ({e:?})"
        );
    }
    assert!(
        e[live + 1].abs() <= tol * u::<S>(),
        "not at rounding level: {e:?}"
    );
}

#[test]
fn renormalize_converges_quadratically() {
    quadratic_convergence::<f64>(1e-3, 2, 4.0);
    quadratic_convergence::<Dual<f64, 4>>(1e-3, 2, 4.0);
    quadratic_convergence::<f32>(2_f64.powi(-8), 1, 4.0);
}

#[test]
fn newton_step_hand_case() {
    // `η = 0.5625`, `k = (3 - 1.5625)/2 = 0.71875`: exact in binary.
    let mut q = raw(1.25_f64, 0.0, 0.0, 0.0);
    q.renormalize();
    assert_eq!(bits(&q), [0.898_437_5, 0.0, 0.0, 0.0].map(f64::to_bits));
    // `η = -1/4`, `k = 9/8`, and `η' = -3/4 η² + 1/4 η³ = -0.05078125`; all exact.
    let mut q = raw(0.5_f64, 0.5, 0.5, 0.0);
    q.renormalize();
    assert_eq!(bits(&q), [0.5625, 0.5625, 0.5625, 0.0].map(f64::to_bits));
    assert_eq!(eta(&q).to_bits(), (-0.050_781_25_f64).to_bits());
    let n = Quat::from_wxyz_normalized(0.5_f64, 0.5, 0.5, 0.0);
    assert_eq!(bits(&n), bits(&q));
}

#[test]
fn blend_selects_both_ways() {
    let (a, b) = (raw(1.0_f64, 2.0, 3.0, 4.0), raw(5.0, 6.0, 7.0, 8.0));
    assert_eq!(bits(&Quat::blend(true, a, b)), bits(&a));
    assert_eq!(bits(&Quat::blend(false, a, b)), bits(&b));
}

proptest! {
    /// Every method under `Dual<f64, 4>` has the value of the plain `f64` one, bit for bit.
    #[test]
    fn dual_value_is_plain_value(a in unit_sample(), b in unit_sample()) {
        let (pa, pb) = (quat::<f64>(&a), quat::<f64>(&b));
        let (da, db) = (quat::<Dual<f64, 4>>(&a), quat::<Dual<f64, 4>>(&b));
        let plain = |q: &Quat<Dual<f64, 4>>| vals(q).map(f64::to_bits);
        prop_assert_eq!(plain(&(da * db)), bits(&(pa * pb)));
        prop_assert_eq!(plain(&da.conjugate()), bits(&pa.conjugate()));
        prop_assert_eq!(da.norm_sq().value_f64().to_bits(), pa.norm_sq().to_bits());
        prop_assert_eq!(mbits(rows(&da)), mbits(rows(&pa)));
        let (mut ra, mut rp) = (da, pa);
        ra.renormalize();
        rp.renormalize();
        prop_assert_eq!(plain(&ra), bits(&rp));
        let w = da.to_xyzw().map(|x| x.value_f64().to_bits());
        prop_assert_eq!(w, pa.to_xyzw().map(f64::to_bits));
        prop_assert_eq!(plain(&Quat::from_jpl(da.to_xyzw())), bits(&Quat::from_jpl(pa.to_xyzw())));
    }

    /// The derivative lanes: `∂(p q)/∂q` is the left-multiplication matrix of `p`, `∂‖q‖²/∂q` is
    /// `2 q`.
    #[test]
    fn dual_derivatives(a in unit_sample(), b in unit_sample()) {
        let (p, q) = (quat::<f64>(&a), quat::<Dual<f64, 4>>(&b));
        let pq = raw(Dual::constant(p.w), Dual::constant(p.x), Dual::constant(p.y), Dual::constant(p.z)) * q;
        let left = [
            [p.w, -p.x, -p.y, -p.z],
            [p.x, p.w, -p.z, p.y],
            [p.y, p.z, p.w, -p.x],
            [p.z, -p.y, p.x, p.w],
        ];
        for (c, comp) in [pq.w, pq.x, pq.y, pq.z].iter().enumerate() {
            prop_assert!(diff::<f64>(&comp.d, &left[c]) <= 8.0, "row {c}");
        }
        prop_assert!(diff::<f64>(&q.norm_sq().d, &b.map(|x| 2.0 * x)) <= 2.0);
    }
}

/// A domain is a `debug_assert!` (D11): the debug profile that `just test` runs panics ...
#[cfg(debug_assertions)]
mod out_of_domain {
    use super::*;

    // `w = 1 + 2^-41` has `‖q‖² = 1 + 2^-40` after rounding: the last value accepted.
    #[test]
    fn the_tolerance_is_inclusive_in_f64() {
        let _ = Quat::from_wxyz_unchecked(1.0 + 2_f64.powi(-41), 0.0, 0.0, 0.0);
        let _ = Quat::from_wxyz_unchecked(
            Dual::<f64, 4>::constant(1.0 + 2_f64.powi(-41)),
            Dual::constant(0.0),
            Dual::constant(0.0),
            Dual::constant(0.0),
        );
    }

    #[test]
    fn the_tolerance_is_inclusive_in_f32() {
        let _ = Quat::from_wxyz_unchecked(1.0 + 2_f32.powi(-17), 0.0, 0.0, 0.0);
    }

    #[test]
    #[should_panic(expected = "Quat::from_wxyz_unchecked")]
    fn unchecked_rejects_just_outside_in_f64() {
        let _ = Quat::from_wxyz_unchecked(1.0 + 2_f64.powi(-40), 0.0, 0.0, 0.0);
    }

    #[test]
    #[should_panic(expected = "Quat::from_wxyz_unchecked")]
    fn unchecked_rejects_just_outside_in_dual() {
        let w = Dual::<f64, 4>::constant(1.0 + 2_f64.powi(-40));
        let z = Dual::constant(0.0);
        let _ = Quat::from_wxyz_unchecked(w, z, z, z);
    }

    #[test]
    #[should_panic(expected = "Quat::from_wxyz_unchecked")]
    fn unchecked_rejects_just_outside_in_f32() {
        let _ = Quat::from_wxyz_unchecked(1.0 + 2_f32.powi(-16), 0.0, 0.0, 0.0);
    }

    #[test]
    #[should_panic(expected = "Quat::from_wxyz_unchecked")]
    fn unchecked_rejects_nan() {
        let _ = Quat::from_wxyz_unchecked(f64::NAN, 0.0, 0.0, 0.0);
    }

    // The lower half of the domain: `‖q‖² < 1` and the zero quaternion.
    #[test]
    #[should_panic(expected = "Quat::from_wxyz_unchecked")]
    fn unchecked_rejects_just_below_in_f64() {
        let _ = Quat::from_wxyz_unchecked(1.0 - 2_f64.powi(-40), 0.0, 0.0, 0.0);
    }

    #[test]
    #[should_panic(expected = "Quat::from_wxyz_unchecked")]
    fn unchecked_rejects_just_below_in_dual() {
        let w = Dual::<f64, 4>::constant(1.0 - 2_f64.powi(-40));
        let z = Dual::constant(0.0);
        let _ = Quat::from_wxyz_unchecked(w, z, z, z);
    }

    #[test]
    #[should_panic(expected = "Quat::from_wxyz_unchecked")]
    fn unchecked_rejects_just_below_in_f32() {
        let _ = Quat::from_wxyz_unchecked(1.0 - 2_f32.powi(-16), 0.0, 0.0, 0.0);
    }

    #[test]
    #[should_panic(expected = "Quat::from_wxyz_unchecked")]
    fn unchecked_rejects_the_zero_quaternion() {
        let _ = Quat::from_wxyz_unchecked(0.0_f64, 0.0, 0.0, 0.0);
    }

    // `w = 1 - 2^-41` squares to `1 - 2^-40 + 2^-82`, which rounds to `1 - 2^-40`.
    #[test]
    fn the_tolerance_is_inclusive_below() {
        let _ = Quat::from_wxyz_unchecked(1.0 - 2_f64.powi(-41), 0.0, 0.0, 0.0);
        let _ = Quat::from_wxyz_unchecked(1.0 - 2_f32.powi(-17), 0.0, 0.0, 0.0);
    }
}

/// ... and a release build never panics: `from_wxyz_unchecked` returns its four values.
#[cfg(not(debug_assertions))]
#[test]
fn out_of_domain_does_not_panic_in_release() {
    let q = Quat::from_wxyz_unchecked(3.0_f64, 0.0, 0.0, f64::NAN);
    assert_eq!(q.w.to_bits(), 3.0_f64.to_bits());
    assert!(q.z.is_nan());
}

/// The step asserts nothing (`NUMERICS.md` §3.6 states no domain): it is `q (3 - ‖q‖²)/2` for
/// every input in every profile, zero at `η = 2` and reversing `q` beyond it (SO.14).
#[test]
fn the_step_is_defined_for_every_input() {
    let step = |w, x, y, z| bits(&Quat::from_wxyz_normalized(w, x, y, z));
    assert_eq!(step(0.0, 1.0, 1.0, 1.0), fb([0.0; 4]));
    assert_eq!(step(2.0, 0.0, 0.0, 0.0), fb([-1.0, -0.0, -0.0, -0.0]));
    assert_eq!(step(0.0, 0.0, 0.0, 0.0), fb([0.0; 4]));
    // In `f32`: `‖q‖² = 2.25`, `k = 0.375`, `w' = 0.5625`.
    let q = Quat::from_wxyz_normalized(1.5_f32, 0.0, 0.0, 0.0);
    assert_eq!(q.w.to_bits(), 0.5625_f32.to_bits());
}

/// The operation order is part of the contract (D16): a bit-exact golden of random bit
/// patterns (computed in IEEE doubles, left to right), which a reassociation changes.
#[test]
fn product_and_norm_sum_left_to_right() {
    const CASES: [([u64; 4], [u64; 4], [u64; 4]); 3] = [
        (
            [
                0xbfb8_61b5_0897_adb0,
                0x3fbe_9a7c_76d6_d7e0,
                0x3feb_2644_263f_ef40,
                0xbfb1_9650_6315_2df0,
            ],
            [
                0x3f90_0f15_dae4_4540,
                0x3fc6_5eda_2613_ad68,
                0xbfe4_2e86_6191_9cb8,
                0x3f98_638e_79cb_9e80,
            ],
            [
                0x3fe0_7567_578b_6be9,
                0xbfa3_6671_4826_087c,
                0x3fad_f611_adc6_a7f6,
                0xbfcd_0ef0_dd4b_e212,
            ],
        ),
        (
            [
                0x3fd0_9fff_39d2_c67c,
                0x3fe2_c022_117d_48a0,
                0xbfe9_f9e1_9c6d_6092,
                0xbfd9_2a25_bdbb_3508,
            ],
            [
                0xbfea_3274_33b0_c176,
                0x3fe3_d137_4f1f_7fce,
                0x3fd8_c297_9a26_b7f4,
                0xbfed_51d5_233c_8bca,
            ],
            [
                0xbfe3_e55d_08fb_b516,
                0x3fe2_7737_65fe_c812,
                0x3ff0_ef46_d045_680d,
                0x3fea_0676_ecd5_5048,
            ],
        ),
        (
            [
                0x3fee_dc41_cb2b_55ea,
                0x3fed_be97_6b9f_bbc4,
                0x3fd3_b3bb_cbd2_9974,
                0x3fcd_9584_74cd_61a0,
            ],
            [
                0xbfe5_eb9d_e303_d53c,
                0xbfef_0a3a_5958_acd2,
                0x3fad_0ffa_8815_c9c0,
                0xbfec_3050_8fb6_fbae,
            ],
            [
                0x3fdb_5632_d8be_5b88,
                0xbffd_b426_bcb0_cc0a,
                0x3fdc_1010_398d_5c9a,
                0xbfe5_01e0_2673_3b06,
            ],
        ),
    ];
    let of = |a: [u64; 4]| {
        raw(
            f64::from_bits(a[0]),
            f64::from_bits(a[1]),
            f64::from_bits(a[2]),
            f64::from_bits(a[3]),
        )
    };
    for (p, q, r) in CASES {
        assert_eq!(bits(&(of(p) * of(q))), bits(&of(r)));
    }
    // `w² + x² + y² + z²` left to right; `(w² + y²) + (x² + z²)` ends one ulp lower.
    let p = of(CASES[0].0);
    assert_eq!(p.norm_sq().to_bits(), 0x3fe7_eecb_d52a_fdd5);
}

/// `R(q)`'s grouping is part of the same contract, and `matrix_hand_cases` cannot pin it: those
/// inputs (`0`, `±0.5`, `±1`, `±2`) make every intermediate exact, so they hold under any
/// regrouping — including the `1 - 2(y² + z²)` diagonal that assumes a unit `q`.
///
/// The quaternion below is unit to `2^-52` (inside the domain of `from_wxyz_unchecked`) and chosen
/// so that all three other groupings of `d = w² - ((x² + y²) + z²)` land one ulp away.
#[test]
fn matrix_sums_group_as_written() {
    let q = Quat::from_wxyz_unchecked(
        f64::from_bits(0x3fdb_04b1_d1e9_1bff),
        f64::from_bits(0xbfe9_74eb_499d_dec2),
        f64::from_bits(0x3fc4_c888_1ff2_cebf),
        f64::from_bits(0xbfd9_cdd1_8fcb_04d6),
    );
    assert_eq!(
        mbits(rows(&q)),
        [
            [
                0x3fe3_e8ac_0881_c992,
                0x3fb5_035f_615d_dd40,
                0x3fe8_ea23_6537_6268
            ],
            [
                0xbfe3_2906_030c_e70c,
                0xbfe2_e817_44c2_c5d3,
                0x3fe1_4dcc_d231_0ee0
            ],
            [
                0x3fe0_2403_b077_ebf8,
                0xbfe9_aef8_6040_c380,
                0xbfd4_6150_026a_1566
            ],
        ]
    );
}
