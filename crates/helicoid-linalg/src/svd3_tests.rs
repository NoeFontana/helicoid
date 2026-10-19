//! Tests of `svd3`: planted spectra, any matrix, the value path of `Dual` and its derivative, two
//! lanes, the signed convention and the magnitude range.
//!
//! **Metrics**, in `f64` from the `S` output, `u` the unit roundoff: `frame = max(|U^T U - I|,
//! |det U - 1|, |V^T V - I|, |det V - 1|)`; the backward error `back = |A - U diag(sigma) V^T|_F /
//! |A|_F`; the singular value error `sv = |sigma - sigma*| / |sigma*|`; the angle of a singular pair
//! `ang_i = min_s max(|u_i - s u_i*|, |v_i - s v_i*|)`, weighted by `gap_i / |sigma*|`, `gap_i` the
//! distance of `|sigma_i*|` to the nearest other `|sigma_j*|` or `-|sigma_j*|` (the corpus metric
//! of `0101` decision 8).
//!
//! **The recorded bounds** are `frame <= FRAME u` and `back`, `sv`, `ang * gap / |sigma|` at most
//! `K u`, one constant at every gap, ties included, about twice the measured worst over 10^6 cases
//! of each proptest (the config's `cases` raised; `SVD3_REPORT=1` prints each new largest), two
//! runs with different seeds: `back` 8.7, `sv` 10.3, `ang` 6.4, and 4.5 `u` on the frames
//! (`measure_svd3`'s families). A planted truth is that of the exact
//! `U diag(sigma) V^T`, so `sv` and `ang` include the rounding of the planted matrix to `S`;
//! `measure_svd3` scores against a double-double reference of the stored matrix instead, per family
//! (`0101`'s large-sample table).

use crate::linalg_tests::{cfg, unit, Lane};
use crate::tests::L2;
use crate::{svd3, Dual, Mat3, Real, Vector};
use core::array;
use proptest::prelude::*;
use std::vec::Vec;

/// A matrix by rows, as `f64`.
type M = [[f64; 3]; 3];
type V3 = [f64; 3];

/// The recorded bound on `back`, `sv` and the weighted angle, in `u`.
const K: f64 = 21.0;
/// The frames' bound, in `u`.
const FRAME: f64 = 10.0;

/// A decomposition: `U` and `V` by rows, `sigma`.
#[derive(Clone, Copy, Debug)]
struct Svd {
    u: M,
    s: V3,
    v: M,
}

fn dot(x: V3, y: V3) -> f64 {
    x[0] * y[0] + x[1] * y[1] + x[2] * y[2]
}

fn cross(x: V3, y: V3) -> V3 {
    [
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0],
    ]
}

fn norm(x: V3) -> f64 {
    dot(x, x).sqrt()
}

fn col(m: &M, c: usize) -> V3 {
    array::from_fn(|r| m[r][c])
}

fn fro(m: &M) -> f64 {
    m.iter().flatten().map(|x| x * x).sum::<f64>().sqrt()
}

fn det(m: &M) -> f64 {
    dot(col(m, 0), cross(col(m, 1), col(m, 2)))
}

/// `a` with every entry rounded to `S`.
fn round<S: Lane>(a: &M) -> M {
    a.map(|r| r.map(|x| S::make(x, [0.0; 2]).lane(0)))
}

fn mat<S: Lane>(a: &M) -> Mat3<S> {
    Mat3::from_rows(a.map(|r| Vector(r.map(|x| S::make(x, [0.0; 2])))))
}

/// `svd3` at `S` of `a`, read back as `f64`.
fn run<S: Lane>(a: &M) -> Svd {
    let (u, s, v) = svd3(&mat::<S>(a));
    let rows = |m: &Mat3<S>| array::from_fn(|r| array::from_fn(|c| m.get(r, c).lane(0)));
    Svd {
        u: rows(&u),
        s: s.0.map(|x| x.lane(0)),
        v: rows(&v),
    }
}

/// `max(|X^T X - I|, |det X - 1|)`.
fn frame_of(x: &M) -> f64 {
    // In double-double: an `f64` metric rounds at the ulp of 1, `2u`, which is the size of what it
    // measures.
    let c = [0, 1, 2].map(|k| col(x, k).map(Dd::of));
    let one = |i: usize| Dd::of(f64::from(u8::from(i / 3 == i % 3)));
    let cr = |a: [Dd; 3], b: [Dd; 3]| -> [Dd; 3] {
        [
            a[1].mul(b[2]).sub(a[2].mul(b[1])),
            a[2].mul(b[0]).sub(a[0].mul(b[2])),
            a[0].mul(b[1]).sub(a[1].mul(b[0])),
        ]
    };
    let det = dd_dot(&c[0], &cr(c[1], c[2])).sub(Dd::of(1.0)).f().abs();
    (0..9)
        .map(|i| dd_dot(&c[i / 3], &c[i % 3]).sub(one(i)).f().abs())
        .fold(det, f64::max)
}

/// `U diag(sigma) V^T`, each entry a dot product of three terms.
fn product(d: &Svd) -> M {
    array::from_fn(|i| array::from_fn(|j| (0..3).map(|k| d.u[i][k] * d.s[k] * d.v[j][k]).sum()))
}

/// The gap of `|s_i|`: its distance to the nearest other `|s_j|` and to `-|s_j|`.
fn gap(s: V3, i: usize) -> f64 {
    (0..3)
        .filter(|&j| j != i)
        .map(|j| (s[i].abs() - s[j].abs()).abs().min(s[i].abs() + s[j].abs()))
        .fold(f64::INFINITY, f64::min)
}

/// `(sv, ang)` of `d` against the truth `t`, both over `|sigma*|`, in units of `u`.
fn forward(d: &Svd, t: &Svd, u: f64) -> (f64, f64) {
    let sn = norm(t.s);
    if sn == 0.0 {
        let sv = if d.s.iter().all(|x| *x == 0.0) {
            0.0
        } else {
            f64::INFINITY
        };
        return (sv, 0.0);
    }
    let sv = norm(array::from_fn(|k| d.s[k] - t.s[k])) / (sn * u);
    let ang = (0..3).fold(0.0_f64, |w, i| {
        let (du, dv) = (col(&d.u, i), col(&d.v, i));
        let (tu, tv) = (col(&t.u, i), col(&t.v, i));
        let e = [1.0, -1.0]
            .map(|s| {
                let a = norm(array::from_fn(|r| du[r] - s * tu[r]));
                let b = norm(array::from_fn(|r| dv[r] - s * tv[r]));
                a.max(b)
            })
            .into_iter()
            .fold(f64::INFINITY, f64::min);
        w.max(e * gap(t.s, i) / (sn * u))
    });
    (sv, ang)
}

/// The convention, truth-free: both frames rotations, the ordering, `sigma_3` of the sign of
/// `det a` (where `det a` is clearly not zero), and the backward error; the worst `back` in `u`.
fn convention<S: Lane>(what: &str, a: &M, d: &Svd) -> Result<f64, TestCaseError> {
    let u = unit::<S>();
    let f = frame_of(&d.u).max(frame_of(&d.v));
    prop_assert!(f <= FRAME * u, "{what}: frame {:.3} u, {d:?}", f / u);
    prop_assert!(
        d.s[0] >= d.s[1] && d.s[1] >= d.s[2].abs(),
        "{what}: order {:?}",
        d.s
    );
    let da = det(a);
    let an = fro(a);
    if da.abs() > 1e3 * u * an * an * an {
        prop_assert!(
            (d.s[2] < 0.0) == (da < 0.0),
            "{what}: sign {:?} det {da:e}",
            d.s
        );
    }
    let p = product(d);
    let r = fro(&array::from_fn(|i| array::from_fn(|j| a[i][j] - p[i][j])));
    let back = if an == 0.0 { 0.0 } else { r / (an * u) };
    prop_assert!(back <= K, "{what}: back {back:.3} u");
    Ok(back)
}

/// The worst of each recorded figure (`back`, `sv`, `ang`), printed under `SVD3_REPORT=1`.
#[allow(clippy::print_stdout)]
fn note(what: &str, w: [f64; 3]) {
    use std::sync::Mutex;
    static WORST: Mutex<[f64; 3]> = Mutex::new([0.0; 3]);
    if std::env::var_os("SVD3_REPORT").is_none() {
        return;
    }
    let mut g = WORST
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    for k in 0..3 {
        if w[k] > g[k] {
            g[k] = w[k];
            std::println!("{what}: new worst {k} {:.3}", w[k]);
        }
    }
}

/// Checks `svd3` at `S` on the planted `a` against its truth.
fn check<S: Lane>(what: &str, a: &M, t: &Svd) -> Result<(), TestCaseError> {
    let a = round::<S>(a);
    let d = run::<S>(&a);
    let back = convention::<S>(what, &a, &d)?;
    let (sv, ang) = forward(&d, t, unit::<S>());
    prop_assert!(
        sv <= K && ang <= K,
        "{what}: sv {sv:.3} ang {ang:.3}, {d:?} {t:?}"
    );
    note(what, [back, sv, ang]);
    Ok(())
}

fn rotation([w, x, y, z]: [f64; 4]) -> M {
    let k = 2.0 / (w * w + x * x + y * y + z * z);
    let (xx, yy, zz) = (k * x * x, k * y * y, k * z * z);
    let (xy, xz, yz) = (k * x * y, k * x * z, k * y * z);
    let (wx, wy, wz) = (k * w * x, k * w * y, k * w * z);
    [
        [1.0 - yy - zz, xy - wz, xz + wy],
        [xy + wz, 1.0 - xx - zz, yz - wx],
        [xz - wy, yz + wx, 1.0 - xx - yy],
    ]
}

/// `U diag(s) V^T` for the rotations of `p`, `q`, and its decomposition in the signed convention
/// (`s` given as `s_1 >= s_2 >= |s_3|`).
fn planted(p: [f64; 4], q: [f64; 4], s: V3) -> (M, Svd) {
    let (u, v) = (rotation(p), rotation(q));
    let a = array::from_fn(|i| array::from_fn(|j| (0..3).map(|k| u[i][k] * s[k] * v[j][k]).sum()));
    (a, Svd { u, s, v })
}

/// A planted matrix: `(what, matrix, truth)`.
fn planted_case() -> impl Strategy<Value = (std::string::String, M, Svd)> {
    let parts = (
        prop::array::uniform4(-1.0..1.0f64),
        prop::array::uniform4(-1.0..1.0f64),
        0..9usize,
        0..=12i32,
        -6..=6i32,
        prop::array::uniform3(0.0..1.0f64),
        any::<bool>(),
    );
    parts
        .prop_filter("two rotations", |p| {
            p.0.iter().map(|x| x * x).sum::<f64>() > 0.05
                && p.1.iter().map(|x| x * x).sum::<f64>() > 0.05
        })
        .prop_map(|(p, q, kind, k, e, r, neg)| {
            let g = 10f64.powi(-k);
            let scale = 2f64.powi(e);
            let mut sorted = r;
            sorted.sort_by(|x, y| y.total_cmp(x));
            let (name, mut s) = match kind {
                0 => ("random", sorted),
                1 => ("bottom", [1.0, 0.3 * (1.0 + g), 0.3]),
                2 => ("top", [1.0 + g, 1.0, 0.3]),
                3 => ("graded", [1.0, g, 0.5 * g]),
                4 => ("rank1", [1.0, 0.0, 0.0]),
                5 => ("rank2", [1.0, 0.5, 0.0]),
                6 => ("rotation", [1.0, 1.0, 1.0]),
                7 => ("near-rotation", [1.0 + 2.0 * g, 1.0 + g, 1.0]),
                _ => ("essential", [1.0, 1.0, 0.0]),
            };
            if neg {
                s[2] = -s[2];
            }
            let (a, t) = planted(p, q, s.map(|x| x * scale));
            (std::format!("{name} k {k} e {e} neg {neg}"), a, t)
        })
}

proptest! {
    #![proptest_config(cfg())]

    /// The recorded bounds on planted spectra (random, a narrow gap at either end, McAdams'
    /// graded failure, rank 1 and 2, a rotation and its neighbourhood, an essential matrix), both
    /// signs of the determinant, in random rotations at a scale of `2^-6 ..= 2^6`, at `f64` and
    /// `f32`.
    #[test]
    fn planted_spectra_meet_the_recorded_bounds((what, a, t) in planted_case()) {
        check::<f64>(&what, &a, &t)?;
        check::<f32>(&what, &a, &t)?;
    }
}

/// Entries `m 2^e` with `|m| < 1`, `|e| <= 6`, and the values that make ties, zeros and non-finite
/// input.
fn entry() -> impl Strategy<Value = f64> {
    prop_oneof![
        12 => (-1.0..1.0f64, -6..=6i32).prop_map(|(m, e)| m * 2f64.powi(e)),
        4 => prop::sample::select(std::vec![0.0, -0.0, 1.0, -1.0, 0.5, 2.0]),
        1 => prop::sample::select(std::vec![f64::NAN, f64::INFINITY, f64::NEG_INFINITY]),
    ]
}

fn identity_frames(d: &Svd) -> bool {
    let i: M = array::from_fn(|r| array::from_fn(|c| f64::from(u8::from(r == c))));
    d.u == i && d.v == i
}

fn any_at<S: Lane>(a: &M) -> Result<(), TestCaseError> {
    let a = round::<S>(a);
    let d = run::<S>(&a);
    if a.iter().flatten().any(|x| !x.is_finite()) {
        prop_assert!(
            d.s.iter().all(|x| x.is_nan()) && identity_frames(&d),
            "{a:?}: {d:?}"
        );
        return Ok(());
    }
    prop_assert!(d.s.iter().all(|x| x.is_finite()), "{a:?}: {d:?}");
    convention::<S>("any", &a, &d).map(|_| ())
}

proptest! {
    #![proptest_config(cfg())]

    /// Any matrix, degenerate and non-finite ones included: rotations for frames, the ordering,
    /// the sign of `sigma_3`, the backward error; a non-finite entry gives NaN singular values and
    /// the identity frames.
    #[test]
    fn any_matrix_meets_the_convention(e in prop::array::uniform9(entry())) {
        let a = array::from_fn(|r| array::from_fn(|c| e[3 * r + c]));
        any_at::<f64>(&a)?;
        any_at::<f32>(&a)?;
    }
}

/// Equal to the bit.
fn bits(x: V3, y: V3) -> bool {
    (0..3).all(|k| x[k].to_bits() == y[k].to_bits())
}

fn diag(d: V3) -> M {
    array::from_fn(|i| array::from_fn(|j| if i == j { d[i] } else { 0.0 }))
}

/// `0015` (draft) PH.2's case: `diag(1, 1, -5)` gives `sigma = (5, 1, -1)`, and `U V^T` is the nearest
/// rotation, at squared distance 20 from the input, not 36.
#[test]
fn the_ordering_puts_the_sign_on_the_smallest() {
    fn at<S: Lane>() {
        let a = diag([1.0, 1.0, -5.0]);
        let d = run::<S>(&a);
        assert!(bits(d.s, [5.0, 1.0, -1.0]), "{d:?}");
        let r: M =
            array::from_fn(|i| array::from_fn(|j| (0..3).map(|k| d.u[i][k] * d.v[j][k]).sum()));
        let dist: f64 = (0..3)
            .flat_map(|i| (0..3).map(move |j| (i, j)))
            .map(|(i, j)| (a[i][j] - r[i][j]).powi(2))
            .sum();
        assert!((dist - 20.0).abs() < 1e-5, "{dist}");
    }
    at::<f64>();
    at::<f32>();
}

/// The zero matrix: `sigma = 0`, the identity frames; a non-finite entry: NaN, the identity frames.
#[test]
fn the_zero_matrix_and_non_finite_input() {
    fn at<S: Lane>() {
        let d = run::<S>(&diag([0.0; 3]));
        assert!(bits(d.s, [0.0; 3]) && identity_frames(&d), "{d:?}");
        for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut a = diag([1.0, 2.0, 3.0]);
            a[1][2] = x;
            let d = run::<S>(&a);
            assert!(
                d.s.iter().all(|s| s.is_nan()) && identity_frames(&d),
                "{d:?}"
            );
        }
    }
    at::<f64>();
    at::<f32>();
}

/// A diagonal matrix with distinct entries keeps its axes (up to the sort's quarter turns) and its
/// values exactly; a rotation gives `sigma = (1, 1, 1)` and `U V^T` itself.
#[test]
fn exact_cases() {
    fn at<S: Lane>() {
        let d = run::<S>(&diag([2.0, -3.0, 0.5]));
        assert!(bits(d.s, [3.0, 2.0, -0.5]), "{d:?}");
        let r = rotation([0.9, 0.3, -0.2, 0.25]);
        let r = round::<S>(&r);
        let d = run::<S>(&r);
        let u = unit::<S>();
        assert!(d.s.iter().all(|s| (s - 1.0).abs() <= 8.0 * u), "{d:?}");
        let p: M =
            array::from_fn(|i| array::from_fn(|j| (0..3).map(|k| d.u[i][k] * d.v[j][k]).sum()));
        let e = fro(&array::from_fn(|i| array::from_fn(|j| p[i][j] - r[i][j])));
        assert!(e <= 16.0 * u, "{e:e}");
    }
    at::<f64>();
    at::<f32>();
}

/// Every finite magnitude: a fixed matrix scaled by `2^e` over the whole normal range keeps its
/// decomposition, the singular values scaled exactly by `2^e` where they stay normal.
#[test]
fn every_finite_magnitude_is_in_the_domain() {
    fn at<S: Lane>(lo: i32, hi: i32) {
        let a = [[0.3, -0.7, 0.2], [0.5, 0.1, -0.9], [-0.4, 0.6, 0.8]];
        let a = round::<S>(&a);
        let base = run::<S>(&a);
        for e in lo..=hi {
            let f = 2f64.powi(e);
            let d = run::<S>(&a.map(|r| r.map(|x| x * f)));
            assert!(
                frame_of(&d.u).max(frame_of(&d.v)) <= FRAME * unit::<S>(),
                "2^{e}: {d:?}"
            );
            for k in 0..3 {
                let w = base.s[k] * f;
                assert!(
                    (d.s[k] - w).abs() <= 4.0 * unit::<S>() * w.abs().max(f * 1e-300),
                    "2^{e}: {:?} {:?}",
                    d.s,
                    base.s
                );
            }
        }
    }
    at::<f64>(-1000, 1000);
    at::<f32>(-120, 120);
}

/// A seeded stream of numbers in `[-1, 1]` (splitmix64).
fn lcg(seed: u64) -> impl FnMut() -> f64 {
    let mut state = seed;
    move || {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let z = (state ^ (state >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        let z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        (z ^ (z >> 31)) as f64 / u64::MAX as f64 * 2.0 - 1.0
    }
}

/// Special cases and 100 seeded random matrices.
fn pool() -> Vec<M> {
    let (inf, nan) = (f64::INFINITY, f64::NAN);
    let mut rows: Vec<M> = std::vec![
        diag([0.0; 3]),
        diag([2.0; 3]),
        diag([1.0, 1.0, -5.0]),
        diag([0.0, 0.0, 9.0]),
        diag([1e-310; 3]),
        diag([1e300, -1e300, 1e300]),
        rotation([0.9, 0.3, -0.2, 0.25]),
        [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [-1.0, 0.5, nan]],
        [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [3.0, 6.0, 9.0]],
        [[inf, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    ];
    let mut next = lcg(1);
    rows.extend((0..100).map(|_| {
        array::from_fn(|_| array::from_fn(|_| next() * 2f64.powi((next() * 8.0) as i32)))
    }));
    rows
}

/// Equal to the bit, except that a `NaN` matches any `NaN`.
fn same(x: f64, y: f64) -> bool {
    x.to_bits() == y.to_bits() || x.is_nan() && y.is_nan()
}

fn value_path<S: Real<Mask = bool>, T: Real<Mask = bool>>(
    conv: fn(f64) -> S,
    lift: fn(S, usize) -> T,
) {
    for a in pool() {
        let x = Mat3::from_rows(a.map(|r| Vector(r.map(conv))));
        let y = Mat3::from_rows(
            [0, 1, 2].map(|i| Vector([0, 1, 2].map(|j| lift(x.get(i, j), 3 * i + j)))),
        );
        let ((pu, ps, pv), (du, ds, dv)) = (svd3(&x), svd3(&y));
        for k in 0..3 {
            assert!(same(ps.0[k].value_f64(), ds.0[k].value_f64()), "{a:?}");
        }
        for i in 0..9 {
            let (r, c) = (i % 3, i / 3);
            assert!(
                same(pu.get(r, c).value_f64(), du.get(r, c).value_f64()),
                "{a:?}"
            );
            assert!(
                same(pv.get(r, c).value_f64(), dv.get(r, c).value_f64()),
                "{a:?}"
            );
        }
    }
}

/// The value path of `Dual` is the plain value path, to the bit.
#[test]
fn dual_value_is_plain_value() {
    value_path::<f64, Dual<f64, 9>>(|x| x, Dual::variable);
    value_path::<f32, Dual<f32, 9>>(|x| x as f32, Dual::variable);
    value_path::<f64, Dual<Dual<f64, 9>, 1>>(|x| x, |x, i| Dual::variable(Dual::variable(x, i), 0));
}

/// Both arms of every branch run in every lane: two lanes agree with the scalars, a degenerate lane
/// beside a regular one included.
#[test]
fn two_lanes_agree_with_scalars() {
    let rows = pool();
    let n = rows.len();
    for (i, j) in (0..n).flat_map(|i| [(i, (i * 37 + 11) % n), (i, i), (i, n - 1 - i)]) {
        let x = Mat3::from_rows(
            [0, 1, 2].map(|r| Vector([0, 1, 2].map(|c| L2([rows[i][r][c], rows[j][r][c]])))),
        );
        let (u, s, v) = svd3(&x);
        for (lane, a) in [(0, &rows[i]), (1, &rows[j])] {
            let w = run::<f64>(a);
            assert!((0..3).all(|k| same(s.0[k].0[lane], w.s[k])), "{a:?}");
            assert!(
                (0..9).all(|k| same(u.get(k % 3, k / 3).0[lane], w.u[k % 3][k / 3])
                    && same(v.get(k % 3, k / 3).0[lane], w.v[k % 3][k / 3])),
                "{a:?}"
            );
        }
    }
}

type D9 = Dual<f64, 9>;

/// `svd3` of `a` with one derivative lane per entry.
fn lanes(a: &M) -> (Mat3<D9>, [D9; 3], Mat3<D9>) {
    let x = Mat3::from_rows(
        [0, 1, 2].map(|i| Vector([0, 1, 2].map(|j| D9::variable(a[i][j], 3 * i + j)))),
    );
    let (u, s, v) = svd3(&x);
    (u, s.0, v)
}

/// Where the singular values are simple and no two `|sigma|` are equal, the derivative is the
/// unique first-order one: `d sigma_i = u_i^T E v_i` (the free and exact oracle of `0101` decision
/// 7), `U^T dU` and `V^T dV` skew, and `d(U Sigma V^T) = E`; to `1e-11` on every entry's lane.
#[test]
#[allow(clippy::needless_range_loop)]
fn derivative_is_first_order_perturbation_theory() {
    let mut next = lcg(7);
    let mut cases: Vec<M> = std::vec![
        diag([3.0, 2.0, -0.5]),
        [[0.3, -0.7, 0.2], [0.5, 0.1, -0.9], [-0.4, 0.6, 0.8]]
    ];
    cases.extend((0..50).map(|_| array::from_fn(|_| array::from_fn(|_| next()))));
    for a in cases {
        let d = run::<f64>(&a);
        let min_gap = (0..3).map(|i| gap(d.s, i)).fold(f64::INFINITY, f64::min);
        if min_gap < 0.05 {
            continue;
        }
        let (u, s, v) = lanes(&a);
        for l in 0..9 {
            let (i0, j0) = (l / 3, l % 3);
            let du: M = array::from_fn(|r| array::from_fn(|c| u.get(r, c).d[l]));
            let dv: M = array::from_fn(|r| array::from_fn(|c| v.get(r, c).d[l]));
            for k in 0..3 {
                let want = d.u[i0][k] * d.v[j0][k];
                assert!(
                    (s[k].d[l] - want).abs() <= 1e-11,
                    "{a:?} lane {l} sigma {k}"
                );
            }
            for (x, dx) in [(&d.u, &du), (&d.v, &dv)] {
                for p in 0..3 {
                    for q in 0..3 {
                        let xp = dot(col(x, p), col(dx, q));
                        let xq = dot(col(x, q), col(dx, p));
                        assert!((xp + xq).abs() <= 1e-11, "{a:?} lane {l}: not skew");
                    }
                }
            }
            for r in 0..3 {
                for c in 0..3 {
                    let got: f64 = (0..3)
                        .map(|k| {
                            du[r][k] * d.s[k] * d.v[c][k]
                                + d.u[r][k] * s[k].d[l] * d.v[c][k]
                                + d.u[r][k] * d.s[k] * dv[c][k]
                        })
                        .sum();
                    let want = f64::from(u8::from(r == i0 && c == j0));
                    assert!(
                        (got - want).abs() <= 1e-11,
                        "{a:?} lane {l} entry {r}{c}: {got}"
                    );
                }
            }
        }
    }
}

/// A diagonal input with off-diagonal lanes: the values are already a decomposition, so only the
/// polish gives the vectors their derivative (`0084` decision 4); it is the one of
/// [`derivative_is_first_order_perturbation_theory`].
#[test]
#[allow(clippy::needless_range_loop)]
fn a_diagonal_input_differentiates_its_vectors() {
    let a = diag([3.0, 2.0, 0.5]);
    let d = run::<f64>(&a);
    let (u, s, v) = lanes(&a);
    let mut moved = 0.0_f64;
    for l in [1, 2, 3, 5, 6, 7] {
        let du: M = array::from_fn(|r| array::from_fn(|c| u.get(r, c).d[l]));
        let dv: M = array::from_fn(|r| array::from_fn(|c| v.get(r, c).d[l]));
        moved = moved.max(fro(&du)).max(fro(&dv));
        let (i0, j0) = (l / 3, l % 3);
        for r in 0..3 {
            for c in 0..3 {
                let got: f64 = (0..3)
                    .map(|k| {
                        du[r][k] * d.s[k] * d.v[c][k]
                            + d.u[r][k] * s[k].d[l] * d.v[c][k]
                            + d.u[r][k] * d.s[k] * dv[c][k]
                    })
                    .sum();
                let want = f64::from(u8::from(r == i0 && c == j0));
                assert!((got - want).abs() <= 1e-12, "lane {l} entry {r}{c}: {got}");
            }
        }
    }
    assert!(moved > 0.1, "{moved}");
}

/// At a multiple of a rotation every lane stays finite.
#[test]
fn the_derivative_at_a_tie_is_finite() {
    for a in [
        diag([2.0; 3]),
        rotation([0.9, 0.3, -0.2, 0.25]),
        diag([0.0; 3]),
    ] {
        let (u, s, v) = lanes(&a);
        for k in 0..9 {
            let (r, c) = (k % 3, k / 3);
            assert!(
                u.get(r, c)
                    .d
                    .iter()
                    .chain(&v.get(r, c).d)
                    .all(|x| x.is_finite()),
                "{a:?}"
            );
        }
        assert!(s.iter().flat_map(|x| x.d).all(|x| x.is_finite()), "{a:?}");
    }
}

/// A double-double: `hi + lo`, from the error-free transformations (no fused multiply-add, D16).
#[derive(Clone, Copy, Debug)]
struct Dd(f64, f64);

impl Dd {
    fn of(x: f64) -> Self {
        Self(x, 0.0)
    }
    fn fast(a: f64, b: f64) -> Self {
        let s = a + b;
        Self(s, b - (s - a))
    }
    fn sum(a: f64, b: f64) -> Self {
        let s = a + b;
        let bb = s - a;
        Self(s, (a - (s - bb)) + (b - bb))
    }
    fn prod(a: f64, b: f64) -> Self {
        let split = |x: f64| {
            let c = 134_217_729.0 * x;
            let hi = c - (c - x);
            (hi, x - hi)
        };
        let (p, (ah, al), (bh, bl)) = (a * b, split(a), split(b));
        Self(p, ((ah * bh - p) + ah * bl + al * bh) + al * bl)
    }
    fn add(self, o: Self) -> Self {
        let (s, t) = (Self::sum(self.0, o.0), Self::sum(self.1, o.1));
        let s = Self::fast(s.0, s.1 + t.0);
        Self::fast(s.0, s.1 + t.1)
    }
    fn neg(self) -> Self {
        Self(-self.0, -self.1)
    }
    fn sub(self, o: Self) -> Self {
        self.add(o.neg())
    }
    fn mul(self, o: Self) -> Self {
        let p = Self::prod(self.0, o.0);
        Self::fast(p.0, p.1 + (self.0 * o.1 + self.1 * o.0))
    }
    fn div(self, o: Self) -> Self {
        let q1 = self.0 / o.0;
        let r = self.sub(o.mul(Self::of(q1)));
        let q2 = r.0 / o.0;
        let r = r.sub(o.mul(Self::of(q2)));
        Self::fast(q1, q2).add(Self::of(r.0 / o.0))
    }
    fn sqrt(self) -> Self {
        if self.0 <= 0.0 {
            return Self::of(0.0);
        }
        let x = self.0.sqrt();
        let r = self.sub(Self::prod(x, x));
        Self::fast(x, r.0 / (2.0 * x))
    }
    fn f(self) -> f64 {
        self.0 + self.1
    }
}

fn dd_dot(x: &[Dd; 3], y: &[Dd; 3]) -> Dd {
    (0..3).fold(Dd::of(0.0), |s, k| s.add(x[k].mul(y[k])))
}

/// The decomposition of `a` by one-sided Jacobi in double-double, to about `1e-30 |A|`: the
/// reference of [`measure_svd3`], independent of `svd3`'s rotation (the textbook `zeta`, `t`, `c`),
/// put in the signed convention.
fn reference(a: &M) -> Svd {
    let mut b: [[Dd; 3]; 3] = array::from_fn(|j| array::from_fn(|i| Dd::of(a[i][j])));
    let mut v: [[Dd; 3]; 3] =
        array::from_fn(|j| array::from_fn(|i| Dd::of(f64::from(u8::from(i == j)))));
    let one = Dd::of(1.0);
    for _ in 0..40 {
        let mut off = 0.0_f64;
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            let (al, be, ga) = (
                dd_dot(&b[p], &b[p]),
                dd_dot(&b[q], &b[q]),
                dd_dot(&b[p], &b[q]),
            );
            let scale = (al.0 * be.0).sqrt();
            if scale == 0.0 || ga.0.abs() <= 1e-33 * scale {
                continue;
            }
            off = off.max(ga.0.abs() / scale);
            let zeta = be.sub(al).div(ga.mul(Dd::of(2.0)));
            let az = if zeta.0 < 0.0 { zeta.neg() } else { zeta };
            let t = one.div(az.add(zeta.mul(zeta).add(one).sqrt()));
            let t = if zeta.0 < 0.0 { t.neg() } else { t };
            let c = one.div(t.mul(t).add(one).sqrt());
            let s = t.mul(c);
            for m in [&mut b, &mut v] {
                let (x, y) = (m[p], m[q]);
                m[p] = array::from_fn(|k| c.mul(x[k]).sub(s.mul(y[k])));
                m[q] = array::from_fn(|k| s.mul(x[k]).add(c.mul(y[k])));
            }
        }
        if off <= 1e-32 {
            break;
        }
    }
    let n = [0, 1, 2].map(|j| dd_dot(&b[j], &b[j]).sqrt());
    let mut order = [0, 1, 2];
    order.sort_by(|&i, &j| n[j].0.total_cmp(&n[i].0).then(n[j].1.total_cmp(&n[i].1)));
    let mut s = order.map(|j| n[j].f());
    let mut u: M = array::from_fn(|r| {
        array::from_fn(|c| {
            let j = order[c];
            if n[j].0 == 0.0 {
                f64::from(u8::from(r == c))
            } else {
                b[j][r].div(n[j]).f()
            }
        })
    });
    let mut vv: M = array::from_fn(|r| array::from_fn(|c| v[order[c]][r].f()));
    if det(&vv) < 0.0 {
        (0..3).for_each(|r| vv[r][2] = -vv[r][2]);
        s[2] = -s[2];
    }
    if det(&u) < 0.0 {
        (0..3).for_each(|r| u[r][2] = -u[r][2]);
        s[2] = -s[2];
    }
    Svd { u, s, v: vv }
}

/// The corpus metric (`0101` decision 8) of `d` against `t`, in `u`.
fn corpus_metric(a: &M, d: &Svd, t: &Svd, u: f64) -> f64 {
    let (sv, ang) = forward(d, t, u);
    let an = fro(a);
    let p = product(d);
    let back = if an == 0.0 {
        0.0
    } else {
        fro(&array::from_fn(|i| array::from_fn(|j| a[i][j] - p[i][j]))) / (an * u)
    };
    sv.max(ang)
        .max(back)
        .max(frame_of(&d.u).max(frame_of(&d.v)) / u)
}

/// The families of [`measure_svd3`].
#[allow(clippy::type_complexity)]
fn families() -> Vec<(
    std::string::String,
    std::boxed::Box<dyn Fn(&mut dyn FnMut() -> f64) -> M>,
)> {
    use std::boxed::Box;
    use std::format;
    let q = |r: &mut dyn FnMut() -> f64| [r(), r(), r(), r()];
    let spec = |s: V3| -> Box<dyn Fn(&mut dyn FnMut() -> f64) -> M> {
        Box::new(move |r: &mut dyn FnMut() -> f64| planted(q(r), q(r), s).0)
    };
    let mut f: Vec<(
        std::string::String,
        Box<dyn Fn(&mut dyn FnMut() -> f64) -> M>,
    )> = std::vec![
        (
            "random".into(),
            Box::new(|r: &mut dyn FnMut() -> f64| array::from_fn(|_| array::from_fn(|_| r())))
        ),
        (
            "col-graded".into(),
            Box::new(|r: &mut dyn FnMut() -> f64| {
                let d = [1.0, 1e-5, 1e-10];
                array::from_fn(|_| array::from_fn(|c| r() * d[c]))
            }),
        ),
        ("rank1".into(), spec([1.0, 0.0, 0.0])),
        ("rank2".into(), spec([1.0, 0.5, 0.0])),
        ("essential".into(), spec([1.0, 1.0, 0.0])),
        ("ordering".into(), spec([5.0, 1.0, -1.0])),
        ("rotation".into(), spec([1.0, 1.0, 1.0])),
    ];
    for k in [0, 2, 4, 8, 12] {
        let g = 10f64.powi(-k);
        f.push((format!("bottom-1e-{k}"), spec([1.0, 0.3 * (1.0 + g), 0.3])));
        f.push((format!("top-1e-{k}"), spec([1.0 + g, 1.0, 0.3])));
    }
    for k in [4, 8] {
        let g = 10f64.powi(-k);
        f.push((format!("graded-1e-{k}"), spec([1.0, g, 0.5 * g])));
        let e = 10f64.powi(-k);
        f.push((
            format!("near-rotation-1e-{k}"),
            Box::new(move |r: &mut dyn FnMut() -> f64| {
                let m = rotation(q(r));
                m.map(|row| row.map(|x| x + e * r()))
            }),
        ));
        f.push((
            format!("near-reflection-1e-{k}"),
            Box::new(move |r: &mut dyn FnMut() -> f64| {
                let m = rotation(q(r));
                m.map(|row| [row[0], row[1], -row[2]].map(|x| x + e * r()))
            }),
        ));
    }
    f
}

fn summary(mut s: Vec<f64>) -> [f64; 3] {
    s.sort_by(f64::total_cmp);
    let n = s.len();
    [
        s[n - 1],
        s[(n * 999) / 1000],
        s.iter().sum::<f64>() / n as f64,
    ]
}

#[test]
#[ignore = "measurement: prints the large-sample table of `0101` and the bounds are recorded from"]
#[allow(clippy::print_stdout)]
fn measure_svd3() {
    fn at<S: Lane>(a: &M) -> f64 {
        let a = round::<S>(a);
        let d = run::<S>(&a);
        corpus_metric(&a, &d, &reference(&a), unit::<S>())
    }
    std::println!("family                 f64 max  p99.9   mean | f32 max  p99.9   mean");
    for (i, (name, make)) in families().into_iter().enumerate() {
        let mut next = lcg(0x7376_6433 + i as u64);
        let (mut s64, mut s32) = (Vec::new(), Vec::new());
        for _ in 0..50_000 {
            let a = make(&mut next);
            s64.push(at::<f64>(&a));
            s32.push(at::<f32>(&a));
        }
        let ([a0, a1, a2], [b0, b1, b2]) = (summary(s64), summary(s32));
        std::println!("{name:22} {a0:7.3} {a1:6.3} {a2:6.3} | {b0:7.3} {b1:6.3} {b2:6.3}");
    }
}
