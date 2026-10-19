//! Tests of `nearest_rotation`: every arm against a double-double reference, any matrix, the value
//! path of `Dual` and its derivative against `maths/svd3.md` PO.3, two lanes.
//!
//! **Metrics**, in `f64` from the `S` output, `u` the unit roundoff: the frame `max(|R^T R - I|,
//! |det R - 1|)`, and the forward error weighted by the conditioning, `|R - R*|_F (sigma_2 +
//! sigma_3) / (2 |sigma|_2)`, `sigma` the signed singular values of the stored input, `R*` its
//! nearest rotation (the corpus metric of `0102` decision 8). The reference is `svd3_tests`'
//! one-sided Jacobi in double-double, so a planted input is scored as stored.
//!
//! **The recorded bounds** are `frame <= FRAME u` and the weighted error at most `K u`, about
//! twice the measured worst over 10^6 cases of each proptest (the config's `cases` raised;
//! `POLAR_REPORT=1` prints each new largest), three runs with different seeds: 3.2 on the weighted
//! error (the `non-unique` family), 7.0 on the frame.

use crate::linalg_tests::{cfg, unit, Lane};
use crate::svd3_tests::{frame_of, fro, lcg, planted, reference, rotation, round, Svd, M, V3};
use crate::tests::L2;
use crate::{nearest_rotation, Dual, Mat3, Real, Vector};
use core::array;
use proptest::prelude::*;
use std::vec::Vec;

/// The recorded bound on the weighted error, in `u`.
const K: f64 = 8.0;
/// The frame's bound, in `u`.
const FRAME: f64 = 14.0;

fn mat<S: Lane>(a: &M) -> Mat3<S> {
    Mat3::from_rows(a.map(|r| Vector(r.map(|x| S::make(x, [0.0; 2])))))
}

fn run<S: Lane>(a: &M) -> M {
    let r = nearest_rotation(&mat::<S>(a));
    array::from_fn(|i| array::from_fn(|j| r.get(i, j).lane(0)))
}

fn uvt(t: &Svd) -> M {
    array::from_fn(|i| array::from_fn(|j| (0..3).map(|k| t.u[i][k] * t.v[j][k]).sum()))
}

/// The weighted forward error of `r` against the reference of `a`, in `u`.
fn weighted(a: &M, r: &M, u: f64) -> f64 {
    let t = reference(a);
    let sn = t.s.iter().map(|x| x * x).sum::<f64>().sqrt();
    let w = (t.s[1] + t.s[2]).max(0.0);
    if sn == 0.0 || w == 0.0 {
        return 0.0;
    }
    let star = uvt(&t);
    fro(&array::from_fn(|i| {
        array::from_fn(|j| r[i][j] - star[i][j])
    })) * w
        / (2.0 * sn * u)
}

#[allow(clippy::print_stdout)]
fn note(what: &str, w: [f64; 2]) {
    use std::sync::Mutex;
    static WORST: Mutex<[f64; 2]> = Mutex::new([0.0; 2]);
    if std::env::var_os("POLAR_REPORT").is_none() {
        return;
    }
    let mut g = WORST
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    for k in 0..2 {
        if w[k] > g[k] {
            g[k] = w[k];
            std::println!("{what}: new worst {k} {:.3}", w[k]);
        }
    }
}

fn check<S: Lane>(what: &str, a: &M) -> Result<(), TestCaseError> {
    let (a, u) = (round::<S>(a), unit::<S>());
    let r = run::<S>(&a);
    let f = frame_of(&r) / u;
    let e = weighted(&a, &r, u);
    prop_assert!(
        f <= FRAME && e <= K,
        "{what}: frame {f:.3} error {e:.3}, {a:?} -> {r:?}"
    );
    note(what, [e, f]);
    Ok(())
}

/// The families of the header: a rotation and its neighbourhood, a scaled rotation, random
/// entries (both signs of `det`), a reflection with an isolated `sigma_3`, `sigma_3` small enough to
/// cross `tau`, rank 2, and the non-unique `(1, a, -a)`.
fn case() -> impl Strategy<Value = (std::string::String, M)> {
    (
        prop::array::uniform4(-1.0..1.0f64),
        prop::array::uniform4(-1.0..1.0f64),
        prop::array::uniform9(-1.0..1.0f64),
        0..8usize,
        0..=16i32,
        -6..=6i32,
        0.05..1.0f64,
    )
        .prop_filter("two rotations", |p| {
            p.0.iter().map(|x| x * x).sum::<f64>() > 0.05
                && p.1.iter().map(|x| x * x).sum::<f64>() > 0.05
        })
        .prop_map(|(p, q, e, kind, k, ex, a)| {
            let g = 10f64.powi(-k);
            let scale = 2f64.powi(ex);
            let near =
                |m: M| -> M { array::from_fn(|i| array::from_fn(|j| m[i][j] + g * e[3 * i + j])) };
            let (name, m) = match kind {
                0 => ("near-rotation", near(rotation(p))),
                1 => ("scaled", rotation(p).map(|r| r.map(|x| x * scale))),
                2 => (
                    "random",
                    array::from_fn(|i| array::from_fn(|j| e[3 * i + j])),
                ),
                3 => ("reflection", planted(p, q, [1.0, a, -0.4 * a]).0),
                4 => ("cond+", planted(p, q, [1.0, 0.5, g]).0),
                5 => ("cond-", planted(p, q, [1.0, 0.5, -g]).0),
                6 => ("rank2", planted(p, q, [1.0, a, 0.0]).0),
                _ => ("non-unique", planted(p, q, [1.0, a, -a]).0),
            };
            (std::format!("{name} k {k} e {ex}"), m)
        })
}

proptest! {
    #![proptest_config(cfg())]

    /// The recorded bounds on every arm, at `f64` and `f32`.
    #[test]
    fn every_arm_meets_the_recorded_bounds((what, a) in case()) {
        check::<f64>(&what, &a)?;
        check::<f32>(&what, &a)?;
    }
}

fn entry() -> impl Strategy<Value = f64> {
    prop_oneof![
        12 => (-1.0..1.0f64, -6..=6i32).prop_map(|(m, e)| m * 2f64.powi(e)),
        4 => prop::sample::select(std::vec![0.0, -0.0, 1.0, -1.0, 0.5, 2.0]),
        1 => prop::sample::select(std::vec![f64::NAN, f64::INFINITY, f64::NEG_INFINITY]),
    ]
}

fn any_at<S: Lane>(a: &M) -> Result<(), TestCaseError> {
    let a = round::<S>(a);
    let r = run::<S>(&a);
    if a.iter().flatten().any(|x| !x.is_finite()) {
        prop_assert!(r.iter().flatten().all(|x| x.is_nan()), "{a:?}: {r:?}");
        return Ok(());
    }
    let f = frame_of(&r) / unit::<S>();
    note("any", [0.0, f]);
    prop_assert!(f <= FRAME, "{a:?}: frame {f:.3}, {r:?}");
    Ok(())
}

proptest! {
    #![proptest_config(cfg())]

    /// Any matrix, degenerate and non-finite ones included: a rotation, or NaN for a non-finite
    /// entry.
    #[test]
    fn any_matrix_gives_a_rotation(e in prop::array::uniform9(entry())) {
        let a = array::from_fn(|r| array::from_fn(|c| e[3 * r + c]));
        any_at::<f64>(&a)?;
        any_at::<f32>(&a)?;
    }
}

fn identity() -> M {
    array::from_fn(|i| array::from_fn(|j| f64::from(u8::from(i == j))))
}

/// The zero matrix gives `I`; a rotation is returned to a few `u`; `0015` PH.2's `diag(1, 1, -5)`
/// gives the nearest rotation, `diag(-1, 1, -1)` or another at squared distance 20.
#[test]
fn exact_cases() {
    fn at<S: Lane>() {
        let u = unit::<S>();
        assert!(run::<S>(&[[0.0; 3]; 3]) == identity());
        let r = round::<S>(&rotation([0.9, 0.3, -0.2, 0.25]));
        let got = run::<S>(&r);
        let e = fro(&array::from_fn(|i| array::from_fn(|j| got[i][j] - r[i][j])));
        assert!(e <= 8.0 * u, "{e:e}");
        let a = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -5.0]];
        let got = run::<S>(&a);
        let d: f64 = (0..9)
            .map(|k| (a[k / 3][k % 3] - got[k / 3][k % 3]).powi(2))
            .sum();
        assert!((d - 20.0).abs() <= 1e-4, "{d} {got:?}");
        // The reflect arm with its two smallest singular values 0.0019 apart, `sigma = (1.95,
        // 0.2513, 0.2495)`: one polishing step left 19.5 `u` on the frame at `f32`.
        let a = round::<S>(&[
            [0.5, 1.0, -1.588_624_211_521_113],
            [
                0.012_303_650_824_949_865,
                0.288_420_627_821_753_8,
                0.012_199_184_420_463_211,
            ],
            [-0.245_913_889_249_208_93, 0.0, -0.049_008_597_635_288_24],
        ]);
        let f = frame_of(&run::<S>(&a)) / u;
        assert!(f <= FRAME, "{f}");
    }
    at::<f64>();
    at::<f32>();
}

/// Special cases, one per arm and its edges, and 100 seeded random matrices.
fn pool() -> Vec<M> {
    let (inf, nan) = (f64::INFINITY, f64::NAN);
    let mut rows: Vec<M> = std::vec![
        [[0.0; 3]; 3],
        identity(),
        rotation([0.9, 0.3, -0.2, 0.25]),
        [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, -0.5]],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -5.0]],
        [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [3.0, 6.0, 9.0]],
        [[1.0, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1e-12]],
        [[1e-300, 0.0, 0.0], [0.0, 1e-300, 0.0], [0.0, 0.0, 1e-300]],
        [[1e300, 0.0, 0.0], [0.0, -1e300, 0.0], [0.0, 0.0, -1e300]],
        [[nan, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[1.0, inf, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    ];
    let mut next = lcg(3);
    rows.extend((0..100).map(|_| {
        array::from_fn(|_| array::from_fn(|_| next() * 2f64.powi((next() * 8.0) as i32)))
    }));
    rows
}

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
        let (p, d) = (nearest_rotation(&x), nearest_rotation(&y));
        assert!(
            (0..9).all(|k| same(
                p.get(k % 3, k / 3).value_f64(),
                d.get(k % 3, k / 3).value_f64()
            )),
            "{a:?}"
        );
    }
}

/// The value path of `Dual` is the plain value path, to the bit.
#[test]
fn dual_value_is_plain_value() {
    value_path::<f64, Dual<f64, 9>>(|x| x, Dual::variable);
    value_path::<f32, Dual<f32, 9>>(|x| x as f32, Dual::variable);
    value_path::<f64, Dual<Dual<f64, 9>, 1>>(|x| x, |x, i| Dual::variable(Dual::variable(x, i), 0));
}

/// Both arms of every branch run in every lane: two lanes agree with the scalars, lanes in
/// different arms included.
#[test]
fn two_lanes_agree_with_scalars() {
    let rows = pool();
    let n = rows.len();
    for (i, j) in (0..n).flat_map(|i| [(i, (i * 37 + 11) % n), (i, i), (i, n - 1 - i)]) {
        let x = Mat3::from_rows(
            [0, 1, 2].map(|r| Vector([0, 1, 2].map(|c| L2([rows[i][r][c], rows[j][r][c]])))),
        );
        let r = nearest_rotation(&x);
        for (lane, a) in [(0, &rows[i]), (1, &rows[j])] {
            let w = run::<f64>(a);
            assert!(
                (0..9).all(|k| same(r.get(k % 3, k / 3).0[lane], w[k % 3][k / 3])),
                "{a:?}"
            );
        }
    }
}

/// `maths/svd3.md` PO.3: `dR = U Omega V^T`, `Omega_ij = (C_ij - C_ji) / (sigma_i + sigma_j)`, `C
/// = U^T E V`, signed `sigma`, for the perturbation `E` of entry `l`.
fn closed_form(t: &Svd, l: usize) -> M {
    let (i0, j0) = (l / 3, l % 3);
    let c: M = array::from_fn(|i| array::from_fn(|j| t.u[i0][i] * t.v[j0][j]));
    let om: M = array::from_fn(|i| {
        array::from_fn(|j| {
            if i == j {
                0.0
            } else {
                (c[i][j] - c[j][i]) / (t.s[i] + t.s[j])
            }
        })
    });
    array::from_fn(|r| {
        array::from_fn(|q| {
            (0..3)
                .flat_map(|i| (0..3).map(move |j| (i, j)))
                .map(|(i, j)| t.u[r][i] * om[i][j] * t.v[q][j])
                .sum()
        })
    })
}

/// The tangent where the nearest rotation is unique, against PO.3, to `1e-10` relative: an exact
/// rotation (Newton's arm, where `svd3`'s `U V^T` is `O(1)` wrong), its neighbourhood, a scaled
/// rotation, generic input of both signs of `det`, a reflection whose two largest singular values
/// tie (the reflection arm), and `sigma_3 = 1e-12` (the `svd3` arm, simple singular values).
#[test]
fn the_tangent_is_the_derivative_of_the_projection() {
    type D9 = Dual<f64, 9>;
    let mut next = lcg(11);
    let mut cases: Vec<(&str, M)> = std::vec![
        ("rotation", rotation([0.9, 0.3, -0.2, 0.25])),
        (
            "near-rotation",
            rotation([0.2, -0.7, 0.4, 0.5]).map(|r| r.map(|x| x + 1e-8))
        ),
        (
            "scaled",
            rotation([0.3, 0.3, -0.8, 0.1]).map(|r| r.map(|x| 3.0 * x))
        ),
        (
            "reflection-tie",
            planted(
                [0.4, 0.1, 0.6, -0.3],
                [0.2, 0.9, 0.1, 0.3],
                [2.0, 2.0, -0.5]
            )
            .0
        ),
        (
            "svd-arm",
            planted(
                [0.4, 0.1, 0.6, -0.3],
                [0.2, 0.9, 0.1, 0.3],
                [1.0, 0.5, 1e-12]
            )
            .0
        ),
    ];
    for _ in 0..20 {
        cases.push(("random", array::from_fn(|_| array::from_fn(|_| next()))));
    }
    for (what, a) in cases {
        let t = reference(&a);
        if t.s[1] + t.s[2] < 0.05 {
            continue;
        }
        let x = Mat3::from_rows(
            [0, 1, 2].map(|i| Vector([0, 1, 2].map(|j| D9::variable(a[i][j], 3 * i + j)))),
        );
        let r = nearest_rotation(&x);
        for l in 0..9 {
            let want = closed_form(&t, l);
            let got: M = array::from_fn(|i| array::from_fn(|j| r.get(i, j).d[l]));
            let e = fro(&array::from_fn(|i| {
                array::from_fn(|j| got[i][j] - want[i][j])
            }));
            assert!(e <= 1e-10 * fro(&want).max(1.0), "{what} lane {l}: {e:e}");
        }
    }
}

/// The families of [`measure_nearest_rotation`].
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
            "scaled".into(),
            Box::new(move |r: &mut dyn FnMut() -> f64| {
                let c = 10f64.powf(3.0 * r());
                rotation(q(r)).map(|row| row.map(|x| c * x))
            })
        ),
        ("reflection".into(), spec([1.0, 0.6, -0.2])),
        ("rank2".into(), spec([1.0, 0.5, 0.0])),
    ];
    for k in [2, 4, 8, 12, 16] {
        let e = 10f64.powi(-k);
        f.push((
            format!("near-1e-{k}"),
            Box::new(move |r: &mut dyn FnMut() -> f64| {
                let m = rotation(q(r));
                m.map(|row| row.map(|x| x + e * r()))
            }),
        ));
    }
    for k in [4, 8, 12, 16] {
        f.push((format!("cond-1e-{k}"), spec([1.0, 0.5, 10f64.powi(-k)])));
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
#[ignore = "measurement: prints the large-sample table of `0102`"]
#[allow(clippy::print_stdout)]
fn measure_nearest_rotation() {
    fn at<S: Lane>(a: &M) -> f64 {
        let (a, u) = (round::<S>(a), unit::<S>());
        let r = run::<S>(&a);
        weighted(&a, &r, u).max(frame_of(&r) / u)
    }
    std::println!("family                 f64 max  p99.9   mean | f32 max  p99.9   mean");
    for (i, (name, make)) in families().into_iter().enumerate() {
        let mut next = lcg(0x706f_6c61 + i as u64);
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
