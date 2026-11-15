//! Tests of `rq3` and `lq3`: the contract on any matrix, planted families, the exact triangles,
//! the flip, the value path of `Dual` and its derivative, two lanes, and the magnitude range.
//!
//! **Metrics**, in `f64` from the `S` output, `u` the unit roundoff: the rowwise backward error
//! `back = max_i |m_i - (U R)_i| / |m_i|` and the frame `frame = max(|R^T R - I|_F, |det R - 1|)`
//! (`0117` decision 9's metric without the forward terms, which the corpus scores against its
//! reference).
//!
//! **The recorded bounds** are `back <= K u` and `frame <= FRAME u`, about twice the worst measured
//! by `measure_rq3` over 10^6 draws of each family at each precision, two seeds: `back` 8.0 and
//! `frame` 7.9 (`f64`), 6.6 and 7.5 (`f32`). The metric itself rounds in `f64`, which costs a few
//! `u` at `f64` and nothing at `f32`.

use crate::linalg_tests::{unit, Lane};
use crate::tests::L2;
use crate::{lq3, rq3, Dual, LowerTri3, Mat3, Real, UpperTri3, Vector};
use core::array;
use proptest::prelude::*;
use std::vec::Vec;

type M = [[f64; 3]; 3];
type V3 = [f64; 3];

/// The recorded bound on `back`, in `u`.
const K: f64 = 16.0;
/// The recorded bound on `frame`, in `u`.
const FRAME: f64 = 16.0;

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

fn mm(a: &M, b: &M) -> M {
    array::from_fn(|i| array::from_fn(|j| (0..3).map(|k| a[i][k] * b[k][j]).sum()))
}

fn tr(a: &M) -> M {
    array::from_fn(|i| array::from_fn(|j| a[j][i]))
}

fn fro(m: &M) -> f64 {
    m.iter().flatten().map(|x| x * x).sum::<f64>().sqrt()
}

fn diag(d: V3) -> M {
    array::from_fn(|i| array::from_fn(|j| if i == j { d[i] } else { 0.0 }))
}

const I3: M = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// `a` with every entry rounded to `S`.
fn round<S: Lane>(a: &M) -> M {
    a.map(|r| r.map(|x| S::make(x, [0.0; 2]).lane(0)))
}

fn mat<S: Lane>(a: &M) -> Mat3<S> {
    Mat3::from_rows(a.map(|r| Vector(r.map(|x| S::make(x, [0.0; 2])))))
}

fn get<S: Real>(m: &Mat3<S>) -> M {
    array::from_fn(|i| array::from_fn(|j| m.get(i, j).value_f64()))
}

/// A factorization: the triangle and the rotation by rows, and the mask.
#[derive(Clone, Copy, Debug)]
struct Fac {
    t: M,
    r: M,
    valid: bool,
}

fn run<S: Lane + Real<Mask = bool>>(a: &M) -> Fac {
    let (u, r, valid) = rq3(&mat::<S>(a));
    Fac {
        t: get(&u.to_matrix()),
        r: get(&r),
        valid,
    }
}

fn run_lq<S: Lane + Real<Mask = bool>>(a: &M) -> Fac {
    let (l, q, valid) = lq3(&mat::<S>(a));
    Fac {
        t: get(&l.to_matrix()),
        r: get(&q),
        valid,
    }
}

fn det(m: &M) -> f64 {
    dot(m[0], cross(m[1], m[2]))
}

/// `(back, frame)` in `u`.
fn metric(a: &M, f: &Fac, u: f64) -> (f64, f64) {
    let p = mm(&f.t, &f.r);
    let back = (0..3)
        .map(|i| {
            let e = norm(array::from_fn(|c| a[i][c] - p[i][c]));
            let n = norm(a[i]);
            if n == 0.0 {
                if e == 0.0 {
                    0.0
                } else {
                    f64::INFINITY
                }
            } else {
                e / n
            }
        })
        .fold(0.0, f64::max);
    let g = mm(&tr(&f.r), &f.r);
    let off = fro(&array::from_fn(|i| array::from_fn(|j| g[i][j] - I3[i][j])));
    (back / u, off.max((det(&f.r) - 1.0).abs()) / u)
}

/// The contract on a finite input, beyond the bounds: a rotation, the structure, the signs and the
/// mask as `u11 > 0 && u22 > 0 && u33 > 0`.
fn contract<S: Lane + Real<Mask = bool>>(
    what: &str,
    a: &M,
    f: &Fac,
    upper: bool,
) -> Result<(f64, f64), TestCaseError> {
    let (back, frame) = metric(a, f, unit::<S>());
    prop_assert!(back <= K, "{what} {a:?}: back {back}");
    prop_assert!(frame <= FRAME, "{what} {a:?}: frame {frame}");
    prop_assert!(det(&f.r) > 0.0, "{what} {a:?}: det R");
    let zeros = if upper {
        [f.t[1][0], f.t[2][0], f.t[2][1]]
    } else {
        [f.t[0][1], f.t[0][2], f.t[1][2]]
    };
    prop_assert!(zeros.iter().all(|z| z.to_bits() == 0), "{what} {a:?}");
    let d = [f.t[0][0], f.t[1][1], f.t[2][2]];
    let (free, signed) = if upper { ([1, 2], 0) } else { ([0, 1], 2) };
    prop_assert!(free.iter().all(|&k| d[k] >= 0.0), "{what} {a:?}: {d:?}");
    prop_assert_eq!(f.valid, d.iter().all(|&x| x > 0.0), "{} {:?}", what, a);
    let _ = signed;
    Ok((back, frame))
}

// ---------- draws ----------

fn lcg(seed: u64) -> impl FnMut() -> f64 {
    let mut state = seed;
    move || {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let z = (state ^ (state >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        let z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        (z ^ (z >> 31)) as f64 / u64::MAX as f64 * 2.0 - 1.0
    }
}

fn gauss(next: &mut impl FnMut() -> f64) -> f64 {
    let a = 0.5 * (next() + 1.0) + f64::MIN_POSITIVE;
    let b = next();
    (-2.0 * a.ln()).sqrt() * (core::f64::consts::PI * b).cos()
}

fn gaussian(next: &mut impl FnMut() -> f64) -> M {
    array::from_fn(|_| array::from_fn(|_| gauss(next)))
}

/// A family: its name and a draw. `f32` takes the binary32 parameters.
type Draw = std::boxed::Box<dyn Fn(&mut dyn FnMut() -> f64) -> M>;
type Family = (std::string::String, Draw);

fn families(f32: bool) -> Vec<Family> {
    let mut v: Vec<Family> = Vec::new();
    let mut add = |name: std::string::String, f: Draw| v.push((name, f));
    add(
        "random".into(),
        std::boxed::Box::new(|n| gaussian(&mut || n())),
    );
    add(
        "near-identity-1e-3".into(),
        std::boxed::Box::new(|n| {
            let g = gaussian(&mut || n());
            array::from_fn(|i| array::from_fn(|j| I3[i][j] + 1e-3 * g[i][j]))
        }),
    );
    let ks: &[i32] = if f32 { &[2, 4, 6] } else { &[4, 8, 12] };
    for &k in ks {
        add(
            std::format!("near-triangular-1e-{k}"),
            std::boxed::Box::new(move |n| {
                let g = gaussian(&mut || n());
                array::from_fn(|i| {
                    array::from_fn(|j| match i.cmp(&j) {
                        core::cmp::Ordering::Equal => g[i][j].abs() + 0.5,
                        core::cmp::Ordering::Less => g[i][j],
                        core::cmp::Ordering::Greater => g[i][j] * 10f64.powi(-k),
                    })
                })
            }),
        );
    }
    let pk: &[i32] = if f32 {
        &[2, 4, 6, 7, 8, 9]
    } else {
        &[2, 4, 8, 12, 15, 16, 17]
    };
    for &k in pk {
        add(
            std::format!("parallel-1e-{k}"),
            std::boxed::Box::new(move |n| {
                let mut m = gaussian(&mut || n());
                let a = m[2];
                let t: V3 = array::from_fn(|_| gauss(&mut || n()));
                let p: V3 = array::from_fn(|i| t[i] - dot(t, a) / dot(a, a) * a[i]);
                let (na, np, e) = (norm(a), norm(p), 10f64.powi(-k));
                let s = 0.5 * (n() + 1.0) + 0.5;
                m[1] = array::from_fn(|i| s * (e.cos() * a[i] / na + e.sin() * p[i] / np));
                m
            }),
        );
    }
    let grade = if f32 {
        [1.0, 1e-2, 1e-4]
    } else {
        [1.0, 1e-5, 1e-10]
    };
    add(
        "row-graded".into(),
        std::boxed::Box::new(move |n| {
            gaussian(&mut || n()).map(|r| r).iter().enumerate().fold(
                [[0.0; 3]; 3],
                |mut m, (i, r)| {
                    m[i] = r.map(|x| x * grade[i]);
                    m
                },
            )
        }),
    );
    add(
        "col-graded".into(),
        std::boxed::Box::new(move |n| {
            gaussian(&mut || n()).map(|r| array::from_fn(|j| r[j] * grade[j]))
        }),
    );
    for &k in ks {
        add(
            std::format!("near-singular-1e-{k}"),
            std::boxed::Box::new(move |n| {
                let mut m = gaussian(&mut || n());
                let (a, b) = (gauss(&mut || n()), gauss(&mut || n()));
                let t: V3 = array::from_fn(|_| gauss(&mut || n()));
                m[0] = array::from_fn(|i| a * m[1][i] + b * m[2][i] + 10f64.powi(-k) * t[i]);
                m
            }),
        );
    }
    add(
        "negative-det".into(),
        std::boxed::Box::new(|n| {
            let mut m = gaussian(&mut || n());
            if det(&m) > 0.0 {
                m[0] = m[0].map(|x| -x);
            }
            m
        }),
    );
    add(
        "rank-bottom".into(),
        std::boxed::Box::new(|n| {
            let mut m = gaussian(&mut || n());
            m[1] = m[2].map(|x| x * 0.5);
            m
        }),
    );
    v
}

// ---------- the contract ----------

/// Entries `m 2^e` with `|m| < 1`, `|e| <= 6`, and the values that make ties, zeros and non-finite
/// input.
fn entry() -> impl Strategy<Value = f64> {
    prop_oneof![
        12 => (-1.0..1.0f64, -6..=6i32).prop_map(|(m, e)| m * 2f64.powi(e)),
        4 => prop::sample::select(std::vec![0.0, -0.0, 1.0, -1.0, 0.5, 2.0]),
        1 => prop::sample::select(std::vec![f64::NAN, f64::INFINITY, f64::NEG_INFINITY]),
    ]
}

fn any_at<S: Lane + Real<Mask = bool>>(a: &M) -> Result<(), TestCaseError> {
    let a = round::<S>(a);
    for (f, upper) in [(run::<S>(&a), true), (run_lq::<S>(&a), false)] {
        if a.iter().flatten().all(|x| x.is_finite()) {
            prop_assert!(f
                .t
                .iter()
                .chain(f.r.iter())
                .flatten()
                .all(|x| x.is_finite()));
            contract::<S>("any", &a, &f, upper)?;
        } else {
            prop_assert!(!f.valid);
            prop_assert!(f.t == I3 && f.r == I3, "{a:?}: {f:?}");
        }
    }
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 4096, failure_persistence: None, ..ProptestConfig::default() })]

    /// Every input meets the contract at both precisions, `lq3` with the rows reversed.
    #[test]
    fn any_matrix_meets_the_contract(e in prop::array::uniform9(entry())) {
        let a: M = array::from_fn(|i| array::from_fn(|j| e[3 * i + j]));
        any_at::<f64>(&a)?;
        any_at::<f32>(&a)?;
    }

    /// The families of `0117`'s study, both functions, both precisions.
    #[test]
    fn every_family_meets_the_recorded_bounds(seed in any::<u64>()) {
        for f32 in [false, true] {
            for (name, draw) in families(f32) {
                let mut next = lcg(seed);
                let a = draw(&mut next);
                if f32 {
                    let a = round::<f32>(&a);
                    contract::<f32>(&name, &a, &run::<f32>(&a), true)?;
                    contract::<f32>(&name, &a, &run_lq::<f32>(&a), false)?;
                } else {
                    contract::<f64>(&name, &a, &run::<f64>(&a), true)?;
                    contract::<f64>(&name, &a, &run_lq::<f64>(&a), false)?;
                }
            }
        }
    }

    /// `lq3` is `rq3` of the reversed matrix, reversed, bit for bit (`NUMERICS.md` §14).
    #[test]
    fn lq3_matches_reference(e in prop::array::uniform9(entry())) {
        fn at<S: Lane + Real<Mask = bool>>(a: &M) -> Result<(), TestCaseError> {
            let m = mat::<S>(a);
            let flip = |x: &Mat3<S>| Mat3::from_rows([2, 1, 0].map(|i| Vector([x.get(i, 2), x.get(i, 1), x.get(i, 0)])));
            let (l, q, valid) = lq3(&m);
            let (u, r, v2) = rq3(&flip(&m));
            let (lt, ut, qt, rt) = (l.to_matrix(), flip(&u.to_matrix()), q, flip(&r));
            prop_assert_eq!(valid, v2);
            for k in 0..9 {
                let (i, j) = (k / 3, k % 3);
                prop_assert_eq!(lt.get(i, j).value_f64().to_bits(), ut.get(i, j).value_f64().to_bits());
                prop_assert_eq!(qt.get(i, j).value_f64().to_bits(), rt.get(i, j).value_f64().to_bits());
            }
            Ok(())
        }
        let a: M = array::from_fn(|i| array::from_fn(|j| e[3 * i + j]));
        at::<f64>(&a)?;
        at::<f32>(&a)?;
    }
}

// ---------- exactness ----------

fn bits(x: &M, y: &M) -> bool {
    x.iter()
        .flatten()
        .zip(y.iter().flatten())
        .all(|(a, b)| a.to_bits() == b.to_bits())
}

/// An upper triangle with a positive diagonal and `+0` below is its own factor and `R = I`, bit for
/// bit, at every magnitude (`maths/rq3.md` RQ.8, RQ.9); `lq3` the same for a lower triangle; `Dual`
/// values the same.
#[test]
fn a_triangle_returns_itself_and_the_identity_bit_for_bit() {
    fn at<S: Lane + Real<Mask = bool>>(lo: i32, hi: i32, seed: u64) {
        let mut next = lcg(seed);
        for _ in 0..20_000 {
            // A diagonal over the whole range, subnormal included, and entries off it within
            // `2^30` of their row's diagonal, inside RQ.9's reach at both precisions.
            let d: V3 = array::from_fn(|_| {
                let e = lo + ((next() + 1.0) * 0.5 * f64::from(hi - lo)) as i32;
                (1.0 + 0.999 * (next() + 1.0) * 0.5) * 2f64.powi(e)
            });
            let (mut t, mut lower) = (diag(d), diag(d));
            for i in 0..3 {
                for j in 0..3 {
                    let x = next() * d[i] * 2f64.powi((30.0 * next()) as i32);
                    if j > i {
                        t[i][j] = x;
                    } else if j < i {
                        lower[i][j] = x;
                    }
                }
            }
            for (a, l) in [(diag(d), diag(d)), (t, lower)] {
                // `+ 0.0` turns an underflowed `-0` into `+0`: the exactness is stated without `-0`.
                let plus = |m: &M| round::<S>(m).map(|r| r.map(|x| x + 0.0));
                let (a, l) = (plus(&a), plus(&l));
                let bad = |a: &M| {
                    a.iter().flatten().any(|x| !x.is_finite()) || (0..3).any(|i| a[i][i] <= 0.0)
                };
                if bad(&a) || bad(&l) {
                    continue;
                }
                let f = run::<S>(&a);
                assert!(bits(&f.t, &a) && bits(&f.r, &I3) && f.valid, "{a:?}: {f:?}");
                let f = run_lq::<S>(&l);
                assert!(bits(&f.t, &l) && bits(&f.r, &I3) && f.valid, "{l:?}: {f:?}");
                let dual = rq3(&Mat3::from_rows(a.map(|r| {
                    Vector(r.map(|x| Dual::<S, 2>::variable(S::make(x, [0.0; 2]), 0)))
                })));
                assert!(bits(&get(&dual.0.to_matrix()), &a) && bits(&get(&dual.1), &I3));
            }
        }
    }
    at::<f64>(-1074, 1023, 1);
    at::<f32>(-149, 127, 2);
}

/// Equal as values: the bits, or two zeros of any sign.
fn value_eq(x: f64, y: f64) -> bool {
    x.to_bits() == y.to_bits() || (x.abs() + y.abs()).to_bits() == 0
}

/// A `-0` off the diagonal keeps the values; only the signs of zeros may change.
#[test]
fn a_negative_zero_off_the_diagonal_keeps_the_values() {
    let a = [[2.0, -0.0, -1.0], [-0.0, 3.0, -0.0], [-0.0, -0.0, 0.25]];
    let f = run::<f64>(&a);
    assert!(f.valid);
    assert!(
        (0..9).all(|k| value_eq(f.r[k / 3][k % 3], I3[k / 3][k % 3])),
        "{f:?}"
    );
    assert!(
        (0..9).all(|k| value_eq(f.t[k / 3][k % 3], a[k / 3][k % 3])),
        "{f:?}"
    );
}

// ---------- degenerate input ----------

/// The zero matrix factors as `U = 0`, `R = I`; a non-finite entry gives `(I, I)`, a clear mask.
#[test]
fn the_zero_matrix_and_non_finite_input() {
    fn at<S: Lane + Real<Mask = bool>>() {
        let f = run::<S>(&[[0.0; 3]; 3]);
        assert!(
            bits(&f.t, &[[0.0; 3]; 3]) && bits(&f.r, &I3) && !f.valid,
            "{f:?}"
        );
        let f = run_lq::<S>(&[[0.0; 3]; 3]);
        assert!(
            bits(&f.t, &[[0.0; 3]; 3]) && bits(&f.r, &I3) && !f.valid,
            "{f:?}"
        );
        for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for k in 0..9 {
                let mut a = diag([1.0, 2.0, 3.0]);
                a[k / 3][k % 3] = x;
                for f in [run::<S>(&a), run_lq::<S>(&a)] {
                    assert!(
                        bits(&f.t, &I3) && bits(&f.r, &I3) && !f.valid,
                        "{a:?}: {f:?}"
                    );
                }
            }
        }
    }
    at::<f64>();
    at::<f32>();
}

/// Dependent rows: `u22 = 0` (or `u33 = 0`) and a completed rotation; the factorization holds.
#[test]
fn dependent_rows_complete_the_rotation() {
    fn at<S: Lane + Real<Mask = bool>>() {
        let cases: [M; 5] = [
            [[1.0, 2.0, 3.0], [0.5, -1.0, 2.0], [1.0, -2.0, 4.0]],
            [[1.0, 2.0, 3.0], [0.0, 0.0, 0.0], [1.0, -2.0, 4.0]],
            [[1.0, 2.0, 3.0], [0.5, -1.0, 2.0], [0.0, 0.0, 0.0]],
            [[1.0, 2.0, 3.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 7.0]],
        ];
        for a in cases {
            let f = run::<S>(&a);
            assert!(!f.valid, "{a:?}");
            assert!(f.t[1][1] == 0.0 || f.t[2][2] == 0.0, "{a:?}: {f:?}");
            let (back, frame) = metric(&a, &f, unit::<S>());
            assert!(back <= K && frame <= FRAME, "{a:?}: {back} {frame}");
        }
        // `m3 = e3` and nothing else: `R = I` by the completion's last of equals.
        let f = run::<S>(&[[0.0; 3], [0.0; 3], [0.0, 0.0, 7.0]]);
        assert!(bits(&f.r, &I3), "{f:?}");
    }
    at::<f64>();
    at::<f32>();
}

/// A negative determinant keeps the factorization: `u11 < 0`, `u22, u33 > 0`, a clear mask.
#[test]
fn a_negative_determinant_keeps_the_factorization() {
    fn at<S: Lane + Real<Mask = bool>>() {
        let mut next = lcg(11);
        for _ in 0..1000 {
            let mut a = round::<S>(&gaussian(&mut next));
            if det(&a) > 0.0 {
                a[0] = a[0].map(|x| -x);
            }
            let f = run::<S>(&a);
            assert!(!f.valid && f.t[0][0] < 0.0 && f.t[1][1] > 0.0 && f.t[2][2] > 0.0);
            let (back, frame) = metric(&a, &f, unit::<S>());
            assert!(back <= K && frame <= FRAME, "{a:?}: {back} {frame}");
        }
    }
    at::<f64>();
    at::<f32>();
}

/// Rows at every magnitude, a row's own scale independent of the others' (`maths/rq3.md` RQ.9):
/// `diag(1e200, 1, 1e-200)` and graded rows to the ends of the range lose nothing.
#[test]
fn every_finite_magnitude_is_in_the_domain() {
    fn at<S: Lane + Real<Mask = bool>>(lo: i32, hi: i32) {
        let mut next = lcg(5);
        for _ in 0..2000 {
            let g = gaussian(&mut next);
            let e: [i32; 3] =
                array::from_fn(|_| lo + ((next() + 1.0) * 0.5 * f64::from(hi - lo)) as i32);
            let a = round::<S>(&array::from_fn(|i| g[i].map(|x| x * 2f64.powi(e[i]))));
            if a.iter().any(|r| r.iter().all(|x| *x == 0.0)) {
                continue;
            }
            let f = run::<S>(&a);
            assert!(
                f.t.iter()
                    .chain(f.r.iter())
                    .flatten()
                    .all(|x| x.is_finite()),
                "{a:?}"
            );
            let (back, frame) = metric(&a, &f, unit::<S>());
            // Subnormal rows keep only their own digits.
            let sub = a
                .iter()
                .flatten()
                .any(|x| *x != 0.0 && x.abs() < 2f64.powi(lo + 60));
            assert!(sub || back <= K, "{a:?}: back {back}");
            assert!(frame <= FRAME, "{a:?}: frame {frame}");
        }
    }
    at::<f64>(-1060, 1020);
    at::<f32>(-140, 124);
}

// ---------- Dual ----------

fn pool() -> Vec<M> {
    let (inf, nan) = (f64::INFINITY, f64::NAN);
    let mut rows: Vec<M> = std::vec![
        [[0.0; 3]; 3],
        diag([2.0; 3]),
        diag([1.0, 1.0, -5.0]),
        diag([1e-310; 3]),
        diag([1e300, 1e-300, 1e300]),
        [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [-1.0, 0.5, nan]],
        [[1.0, 2.0, 3.0], [0.5, -1.0, 2.0], [1.0, -2.0, 4.0]],
        [[-1.0, 0.3, 0.2], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
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
        let ((pu, pr, pv), (du, dr, dv)) = (rq3(&x), rq3(&y));
        assert_eq!(pv, dv, "{a:?}");
        let (pu, du) = (pu.to_row_major(), du.to_row_major());
        assert!(
            (0..6).all(|k| same(pu[k].value_f64(), du[k].value_f64())),
            "{a:?}"
        );
        assert!(
            (0..9).all(|k| same(
                pr.get(k / 3, k % 3).value_f64(),
                dr.get(k / 3, k % 3).value_f64()
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

/// Both arms of every branch run in every lane: two lanes agree with the scalars.
#[test]
fn two_lanes_agree_with_scalars() {
    let rows = pool();
    let n = rows.len();
    for (i, j) in (0..n).flat_map(|i| [(i, (i * 37 + 11) % n), (i, i), (i, n - 1 - i)]) {
        let x = Mat3::from_rows(
            [0, 1, 2].map(|r| Vector([0, 1, 2].map(|c| L2([rows[i][r][c], rows[j][r][c]])))),
        );
        let (u, r, valid) = rq3(&x);
        let (u, r) = (u.to_matrix(), r);
        for (lane, a) in [(0, &rows[i]), (1, &rows[j])] {
            let w = run::<f64>(a);
            assert_eq!(valid.0[lane], w.valid, "{a:?}");
            assert!(
                (0..9).all(|k| same(u.get(k / 3, k % 3).0[lane], w.t[k / 3][k % 3])
                    && same(r.get(k / 3, k % 3).0[lane], w.r[k / 3][k % 3])),
                "{a:?}"
            );
        }
    }
}

type D9 = Dual<f64, 9>;

/// The first-order tangent of `maths/rq3.md` RQ.6 for the direction `dm`: `(dU, dR)`.
fn tangent(f: &Fac, dm: &M) -> (M, M) {
    let u = f.t;
    let mut ui = [[0.0; 3]; 3];
    ui[0][0] = 1.0 / u[0][0];
    ui[1][1] = 1.0 / u[1][1];
    ui[2][2] = 1.0 / u[2][2];
    ui[0][1] = -u[0][1] / (u[0][0] * u[1][1]);
    ui[1][2] = -u[1][2] / (u[1][1] * u[2][2]);
    ui[0][2] = (u[0][1] * u[1][2] - u[0][2] * u[1][1]) / (u[0][0] * u[1][1] * u[2][2]);
    let x = mm(&mm(&ui, dm), &tr(&f.r));
    let w: M = array::from_fn(|i| {
        array::from_fn(|j| match i.cmp(&j) {
            core::cmp::Ordering::Greater => x[i][j],
            core::cmp::Ordering::Less => -x[j][i],
            core::cmp::Ordering::Equal => 0.0,
        })
    });
    let (a, b) = (mm(dm, &tr(&f.r)), mm(&u, &w));
    (
        array::from_fn(|i| array::from_fn(|j| a[i][j] - b[i][j])),
        mm(&w, &f.r),
    )
}

/// Where the factorization is unique the lanes are RQ.6's tangent, to `1e-12` relative on every
/// entry's lane: at random matrices, a diagonal, a triangle, a negative determinant and the set
/// where a Givens order would freeze its tangent, `m3 = (0, b, 0)`.
#[test]
fn derivative_is_first_order_perturbation_theory() {
    let mut next = lcg(7);
    let mut cases: Vec<M> = std::vec![
        diag([3.0, 2.0, 0.5]),
        [[3.0, -1.0, 0.2], [0.0, 2.0, 0.7], [0.0, 0.0, 0.5]],
        [[-1.0, 0.3, 0.2], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
    ];
    cases.extend((0..200).map(|_| gaussian(&mut next)));
    for a in cases {
        let f = run::<f64>(&a);
        let x = Mat3::from_rows(
            [0, 1, 2].map(|i| Vector([0, 1, 2].map(|j| D9::variable(a[i][j], 3 * i + j)))),
        );
        let (u, r, _) = rq3(&x);
        let u = u.to_matrix();
        for l in 0..9 {
            let dm: M = array::from_fn(|i| array::from_fn(|j| f64::from(u8::from(3 * i + j == l))));
            let (du, dr) = tangent(&f, &dm);
            let scale = fro(&du).max(fro(&dr)).max(1.0);
            for k in 0..9 {
                let (i, j) = (k / 3, k % 3);
                let (eu, er) = (
                    (u.get(i, j).d[l] - du[i][j]).abs(),
                    (r.get(i, j).d[l] - dr[i][j]).abs(),
                );
                assert!(
                    eu.max(er) <= 1e-12 * scale,
                    "{a:?} lane {l} entry {i}{j}: {eu} {er}"
                );
            }
        }
    }
}

/// The derivative at dependent rows is finite.
#[test]
fn the_derivative_at_dependent_rows_is_finite() {
    for a in [
        [[1.0, 2.0, 3.0], [0.5, -1.0, 2.0], [1.0, -2.0, 4.0]],
        [[0.0; 3]; 3],
    ] {
        let x = Mat3::from_rows(
            [0, 1, 2].map(|i| Vector([0, 1, 2].map(|j| D9::variable(a[i][j], 3 * i + j)))),
        );
        let (u, r, _) = rq3(&x);
        let (u, r) = (u.to_matrix(), r);
        assert!((0..9).all(|k| u
            .get(k / 3, k % 3)
            .d
            .iter()
            .chain(r.get(k / 3, k % 3).d.iter())
            .all(|x| x.is_finite())));
    }
}

// ---------- the gauge conversion ----------

/// Converting a factorization between the two gauges: from `M = L Q`, `rq3(L Q)` and back through
/// `lq3(U R)` returns `L` and `Q` to the recorded bound times the conditioning of the near-identity
/// matrices of a calibration, within `4 K u`.
#[test]
fn an_lq_product_round_trips_through_rq3() {
    fn at<S: Lane + Real<Mask = bool>>() {
        let mut next = lcg(13);
        for _ in 0..2000 {
            let g = gaussian(&mut next);
            let l: M = array::from_fn(|i| {
                array::from_fn(|j| {
                    if j > i {
                        0.0
                    } else {
                        I3[i][j] + 1e-2 * g[i][j]
                    }
                })
            });
            let q = mm(&diag([1.0; 3]), &{
                let (w, x, y, z) = (1.0, 0.1 * next(), 0.1 * next(), 0.1 * next());
                let n = (w * w + x * x + y * y + z * z).sqrt();
                let (w, x, y, z) = (w / n, x / n, y / n, z / n);
                [
                    [
                        1. - 2. * (y * y + z * z),
                        2. * (x * y - w * z),
                        2. * (x * z + w * y),
                    ],
                    [
                        2. * (x * y + w * z),
                        1. - 2. * (x * x + z * z),
                        2. * (y * z - w * x),
                    ],
                    [
                        2. * (x * z - w * y),
                        2. * (y * z + w * x),
                        1. - 2. * (x * x + y * y),
                    ],
                ]
            });
            let m = round::<S>(&mm(&l, &q));
            let f = run::<S>(&m);
            let back = run_lq::<S>(&round::<S>(&mm(&f.t, &f.r)));
            let lq = run_lq::<S>(&m);
            let e = fro(&array::from_fn(|i| {
                array::from_fn(|j| back.t[i][j] - lq.t[i][j])
            }))
            .max(fro(&array::from_fn(|i| {
                array::from_fn(|j| back.r[i][j] - lq.r[i][j])
            })));
            assert!(f.valid && back.valid);
            assert!(e <= 4.0 * K * unit::<S>(), "{m:?}: {}", e / unit::<S>());
        }
    }
    at::<f64>();
    at::<f32>();
}

// ---------- measurement ----------

fn summary(mut v: Vec<f64>) -> [f64; 3] {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    [
        v[n - 1],
        v[n * 999 / 1000],
        v.iter().sum::<f64>() / n as f64,
    ]
}

#[test]
#[ignore = "measurement: prints the worst readings the bounds are recorded from; run with --release"]
#[allow(clippy::print_stdout)]
fn measure_rq3() {
    let n: usize = std::env::var("RQ3_DRAWS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1_000_000);
    let seed: u64 = std::env::var("RQ3_SEED")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    for f32 in [false, true] {
        std::println!(
            "{}: family             back    frame",
            if f32 { "f32" } else { "f64" }
        );
        for (i, (name, draw)) in families(f32).into_iter().enumerate() {
            let mut next = lcg(0x7271_3300 + i as u64 + (seed << 32));
            let (mut b, mut f) = (Vec::new(), Vec::new());
            for _ in 0..n {
                let (x, y) = if f32 {
                    let a = round::<f32>(&draw(&mut next));
                    metric(&a, &run::<f32>(&a), unit::<f32>())
                } else {
                    let a = draw(&mut next);
                    metric(&a, &run::<f64>(&a), unit::<f64>())
                };
                b.push(x);
                f.push(y);
            }
            std::println!("     {name:22} {:7.3} {:7.3}", summary(b)[0], summary(f)[0]);
        }
    }
}

/// The triangle types: the row-major order, the exact zeros, the transpose and the blend.
#[test]
fn triangle_types_round_trip() {
    let a = [1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0];
    let (u, l) = (UpperTri3::from_row_major(a), LowerTri3::from_row_major(a));
    let b6 = |x: [f64; 6]| x.map(f64::to_bits);
    assert_eq!(b6(u.to_row_major()), b6(a));
    assert_eq!(
        get(&u.to_matrix()),
        [[1.0, 2.0, 3.0], [0.0, 4.0, 5.0], [0.0, 0.0, 6.0]]
    );
    assert_eq!(
        get(&l.to_matrix()),
        [[1.0, 0.0, 0.0], [2.0, 3.0, 0.0], [4.0, 5.0, 6.0]]
    );
    assert_eq!(get(&u.transpose().to_matrix()), tr(&get(&u.to_matrix())));
    assert_eq!(get(&l.transpose().to_matrix()), tr(&get(&l.to_matrix())));
    assert_eq!(b6(u.transpose().transpose().to_row_major()), b6(a));
    assert_eq!(
        u.diagonal().0.map(f64::to_bits),
        [1.0, 4.0, 6.0].map(f64::to_bits)
    );
    assert_eq!(
        l.diagonal().0.map(f64::to_bits),
        [1.0, 3.0, 6.0].map(f64::to_bits)
    );
    assert_eq!(get(&UpperTri3::<f64>::identity().to_matrix()), I3);
    assert_eq!(get(&LowerTri3::<f64>::identity().to_matrix()), I3);
    use crate::Blend;
    let id = UpperTri3::<f64>::identity();
    assert_eq!(b6(UpperTri3::blend(true, u, id).to_row_major()), b6(a));
    assert_eq!(
        b6(UpperTri3::blend(false, u, id).to_row_major()),
        b6(id.to_row_major())
    );
}

/// Beyond RQ.9's reach, `u22` reads 0 and the factorization stays backward stable: the residual is
/// below `u` of its row.
#[test]
fn beyond_the_reach_u22_reads_zero_and_stays_backward_stable() {
    fn at<S: Lane + Real<Mask = bool>>(e: i32) {
        let a = [
            [1.0, 0.5, 0.25],
            [0.0, 2f64.powi(e - 1), 2f64.powi(-30)],
            [0.0, 0.0, 1.0],
        ];
        let f = run::<S>(&a);
        assert!(f.t[1][1] == 0.0 && !f.valid, "{f:?}");
        let (back, frame) = metric(&a, &f, unit::<S>());
        assert!(back <= K && frame <= FRAME, "{back} {frame}");
    }
    at::<f64>(-30 - 482);
    at::<f32>(-30 - 34);
}
