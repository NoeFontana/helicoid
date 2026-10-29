//! Tests of `eig_sym`: exact cases, planted spectra against a double-double reference at every
//! `N`, relative accuracy on graded matrices, the value path of `Dual` and its derivative, two
//! lanes, and the convergence study that sets the sweep cap (`0111`).
//!
//! **Metrics**, in `f64` from the `S` output against [`reference`] of the matrix rounded to `S`:
//! - `abs`, the corpus's (`NUMERICS.md` §11, `Rule::Eigen`): the larger of `|l - l*|_2 / |l*|_2`
//!   and, over the columns, `min |v_k -+ v_k*| gap_k / |l*|_2`, in `u`;
//! - `rel`, for a positive definite matrix (§11's relative reading): the larger of
//!   `|l_k - l_k*| / (kappa l_k*)` and `min |v_k -+ v_k*| min(relgap_k, 1) / kappa`, in `u`,
//!   `kappa` the condition number of `D^-1 A D^-1`, `D = diag(sqrt(a_ii))`, `relgap_k =
//!   min_j |l_k - l_j| / sqrt(l_k l_j)`. A relative gap above 1 claims an angle below the `u` a
//!   unit vector's components are rounded to, so it is capped there;
//! - the frame, `max(|V^T V - I|, |det V - 1|)` in `u`.
//!
//! **The recorded bounds** are twice the worst `measure_eig_sym` reads over its families, at every
//! `N` in `4..=16` and both precisions (`0111` *Measured*; `2e4` matrices of each of 20 families,
//! 300 of them scored): `abs` 12.63 `u` (`N = 16`, `f64`), `rel` 4.25 (`N = 8`, `graded-1e-8`),
//! the frame 20.66 (`N = 15`), so `abs <= ABS`, `rel <= REL`, the frame `<= FRAME`. The relative
//! reading is taken where `kappa u <= 2^-10`, the corpus's (`0111` decision 5), inside Demmel and
//! Veselić's hypothesis.

// Index loops read as the formulas do, as in `svd3_tests`.
#![allow(clippy::needless_range_loop)]

use crate::eig3_tests::Dd;
use crate::eig_sym::jacobi;
use crate::linalg_tests::{unit, Lane};
use crate::tests::L2;
use crate::{eig_sym, Dual, Matrix, Precision, Real, Vector};
use core::array;
use std::vec::Vec;

/// The recorded bounds of the header.
const ABS: f64 = 25.3;
const REL: f64 = 9.0;
const FRAME: f64 = 42.0;

/// A matrix by rows, as `f64`.
type M<const N: usize> = [[f64; N]; N];

/// A seeded stream (splitmix64): uniform in `[-1, 1]`, and standard normal.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let z = (self.0 ^ (self.0 >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        let z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        (z ^ (z >> 31)) as f64 / u64::MAX as f64 * 2.0 - 1.0
    }
    fn gauss(&mut self) -> f64 {
        let u = 0.5 * (self.next() + 1.0);
        let v = self.next();
        (-2.0 * (1.0 - u).max(f64::MIN_POSITIVE).ln()).sqrt() * (core::f64::consts::PI * v).cos()
    }
}

/// A Haar-random orthogonal matrix: Gram–Schmidt, twice, on Gaussian columns (by rows).
fn haar<const N: usize>(r: &mut Rng) -> M<N> {
    let mut q: [[f64; N]; N] = array::from_fn(|_| array::from_fn(|_| r.gauss()));
    for c in 0..N {
        for _ in 0..2 {
            for j in 0..c {
                let d: f64 = (0..N).map(|i| q[c][i] * q[j][i]).sum();
                for i in 0..N {
                    q[c][i] -= d * q[j][i];
                }
            }
        }
        let n = q[c].iter().map(|x| x * x).sum::<f64>().sqrt();
        q[c] = q[c].map(|x| x / n);
    }
    array::from_fn(|i| array::from_fn(|j| q[j][i]))
}

/// `Q diag(l) Q^T`, exactly symmetric.
fn planted<const N: usize>(q: &M<N>, l: &[f64; N]) -> M<N> {
    let mut a = [[0.0; N]; N];
    for i in 0..N {
        for j in 0..=i {
            let x = (0..N).map(|k| q[i][k] * l[k] * q[j][k]).sum();
            a[i][j] = x;
            a[j][i] = x;
        }
    }
    a
}

/// `a` scaled to a unit diagonal, `D^-1 a D^-1`, one root at a time: `a_ii a_jj` underflows on
/// `graded-1e-200`.
fn equilibrated<const N: usize>(a: &M<N>) -> M<N> {
    array::from_fn(|i| array::from_fn(|j| a[i][j] / a[i][i].sqrt() / a[j][j].sqrt()))
}

/// The eigenpairs, `lam` ascending and `vec[k]` a unit eigenvector of `lam[k]`.
#[derive(Debug)]
struct Truth<const N: usize> {
    lam: [f64; N],
    vec: [[f64; N]; N],
}

/// The eigenpairs of `a` by cyclic Jacobi in double-double, independent of `eig_sym`'s rotation
/// (the textbook `theta`, `t`, `c` and the full two-sided update), with a relative threshold, so
/// a graded matrix's small eigenvalues are as accurate as its large ones: about `1e-32 kappa` of
/// each (Demmel and Veselić's bound at double-double's `u`).
fn reference<const N: usize>(a: &M<N>) -> Truth<N> {
    let mut m = a.map(|r| r.map(Dd::of));
    let mut v: [[Dd; N]; N] =
        array::from_fn(|i| array::from_fn(|j| Dd::of(f64::from(u8::from(i == j)))));
    let an = a.iter().flatten().fold(0.0_f64, |s, x| s.max(x.abs()));
    let (one, two) = (Dd::of(1.0), Dd::of(2.0));
    for _ in 0..100 {
        let mut rotated = false;
        for p in 0..N {
            for q in p + 1..N {
                let scale = (m[p][p].0.abs() * m[q][q].0.abs()).sqrt();
                if m[p][q].0.abs() <= (1e-33 * scale).max(1e-300 * an) {
                    continue;
                }
                rotated = true;
                let theta = m[q][q].sub(m[p][p]).div(m[p][q].mul(two));
                let at = if theta.0 < 0.0 { theta.neg() } else { theta };
                // Past Dekker's range `theta^2` overflows; `t = 1 / (2 theta)` there.
                let t = if at.0 > 1e140 {
                    one.div(at.mul(two))
                } else {
                    one.div(at.add(theta.mul(theta).add(one).sqrt()))
                };
                let t = if theta.0 < 0.0 { t.neg() } else { t };
                let c = one.div(t.mul(t).add(one).sqrt());
                let s = t.mul(c);
                let rot = |x: Dd, y: Dd| (c.mul(x).sub(s.mul(y)), s.mul(x).add(c.mul(y)));
                for row in &mut m {
                    (row[p], row[q]) = rot(row[p], row[q]);
                }
                let (mp, mq) = (m[p], m[q]);
                (m[p], m[q]) = (
                    array::from_fn(|k| rot(mp[k], mq[k]).0),
                    array::from_fn(|k| rot(mp[k], mq[k]).1),
                );
                for row in &mut v {
                    (row[p], row[q]) = rot(row[p], row[q]);
                }
            }
        }
        if !rotated {
            break;
        }
    }
    let mut order: [usize; N] = array::from_fn(|k| k);
    order.sort_by(|&i, &j| {
        m[i][i]
            .0
            .total_cmp(&m[j][j].0)
            .then(m[i][i].1.total_cmp(&m[j][j].1))
    });
    Truth {
        lam: order.map(|k| m[k][k].0 + m[k][k].1),
        vec: order.map(|k| {
            let n = (0..N)
                .fold(Dd::of(0.0), |s, r| s.add(v[r][k].mul(v[r][k])))
                .sqrt();
            array::from_fn(|r| v[r][k].div(n).0)
        }),
    }
}

/// `a` with every entry rounded to `S`.
fn round<S: Lane, const N: usize>(a: &M<N>) -> M<N> {
    a.map(|r| r.map(|x| S::make(x, [0.0; 2]).lane(0)))
}

fn matrix<S: Lane, const N: usize>(a: &M<N>) -> Matrix<S, N, N> {
    Matrix::from_rows(a.map(|r| Vector(r.map(|x| S::make(x, [0.0; 2])))))
}

/// `eig_sym` at `S` of `a`: the eigenvalues and the eigenvectors as columns, `v[k]`.
fn run<S: Lane, const N: usize>(a: &M<N>) -> ([f64; N], M<N>) {
    let (l, v) = eig_sym(&matrix::<S, N>(a));
    (
        l.0.map(|x| x.lane(0)),
        array::from_fn(|k| array::from_fn(|r| v.get(r, k).lane(0))),
    )
}

fn dot<const N: usize>(x: &[f64; N], y: &[f64; N]) -> f64 {
    (0..N).map(|i| x[i] * y[i]).sum()
}

fn norm<const N: usize>(x: &[f64; N]) -> f64 {
    dot(x, x).sqrt()
}

/// `min |x -+ y|`.
fn aligned<const N: usize>(x: &[f64; N], y: &[f64; N]) -> f64 {
    norm::<N>(&array::from_fn(|r| x[r] - y[r])).min(norm::<N>(&array::from_fn(|r| x[r] + y[r])))
}

/// The determinant of `v`'s columns, by Gaussian elimination with partial pivoting.
fn det<const N: usize>(v: &M<N>) -> f64 {
    let mut a = *v;
    let mut d = 1.0;
    for c in 0..N {
        let p = (c..N)
            .max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))
            .unwrap_or(c);
        if p != c {
            a.swap(p, c);
            d = -d;
        }
        d *= a[c][c];
        for i in c + 1..N {
            let f = a[i][c] / a[c][c];
            for j in c..N {
                a[i][j] -= f * a[c][j];
            }
        }
    }
    d
}

/// The frame, `max(|V^T V - I|, |det V - 1|)` in `u`.
fn frame<S: Real, const N: usize>(v: &M<N>) -> f64 {
    let mut w = (det(v) - 1.0).abs();
    for i in 0..N {
        for j in 0..N {
            w = w.max((dot(&v[i], &v[j]) - f64::from(u8::from(i == j))).abs());
        }
    }
    w / unit::<S>()
}

/// The corpus metric, in `u`.
fn abs_metric<const N: usize>(t: &Truth<N>, (l, v): &([f64; N], M<N>), u: f64) -> f64 {
    let ln = norm(&t.lam);
    if ln == 0.0 {
        return if l.iter().all(|x| *x == 0.0) {
            0.0
        } else {
            f64::INFINITY
        };
    }
    let ev = norm::<N>(&array::from_fn(|k| l[k] - t.lam[k])) / (ln * u);
    (0..N).fold(ev, |w, k| {
        let gap = (0..N)
            .filter(|&j| j != k)
            .map(|j| (t.lam[k] - t.lam[j]).abs())
            .fold(f64::INFINITY, f64::min);
        w.max(aligned(&v[k], &t.vec[k]) * gap / (ln * u))
    })
}

/// The relative reading of a positive definite `a`, in `u`; `kappa` that of the equilibrated `a`.
fn rel_metric<const N: usize>(t: &Truth<N>, kappa: f64, (l, v): &([f64; N], M<N>), u: f64) -> f64 {
    (0..N).fold(0.0_f64, |w, k| {
        let relgap = (0..N)
            .filter(|&j| j != k)
            .map(|j| (t.lam[k] - t.lam[j]).abs() / (t.lam[k] * t.lam[j]).sqrt())
            .fold(1.0, f64::min);
        let ev = (l[k] - t.lam[k]).abs() / (kappa * t.lam[k] * u);
        w.max(ev)
            .max(aligned(&v[k], &t.vec[k]) * relgap / (kappa * u))
    })
}

/// The condition number of the equilibrated `a`, from its reference eigenvalues.
fn kappa<const N: usize>(a: &M<N>) -> f64 {
    let t = reference(&equilibrated(a));
    t.lam[N - 1] / t.lam[0]
}

/// A family of the convergence study: a name, whether its matrices are positive definite, and
/// its generator.
type Family<const N: usize> = (&'static str, bool, fn(&mut Rng) -> M<N>);

fn spectrum<const N: usize>(r: &mut Rng) -> [f64; N] {
    let mut l: [f64; N] = array::from_fn(|_| r.next());
    l.sort_by(f64::total_cmp);
    l
}

fn gap<const N: usize>(r: &mut Rng, g: f64, top: bool) -> M<N> {
    let mut l = spectrum::<N>(r);
    if top {
        l[N - 1] = l[N - 2] + g;
    } else {
        l[1] = l[0] + g;
    }
    planted(&haar(r), &l)
}

fn cluster<const N: usize>(r: &mut Rng, m: usize) -> M<N> {
    let mut l = spectrum::<N>(r);
    let c = l[0];
    for x in l.iter_mut().take(m) {
        *x = c;
    }
    planted(&haar(r), &l)
}

fn rank_def<const N: usize>(r: &mut Rng, m: usize) -> M<N> {
    let mut l = spectrum::<N>(r);
    for x in l.iter_mut().take(m) {
        *x = 0.0;
    }
    planted(&haar(r), &l)
}

/// `D A D`, `A` positive definite with a unit diagonal and eigenvalues in `[1/4, 7/4]`, `D` graded
/// over `10^-k` geometrically.
fn graded<const N: usize>(r: &mut Rng, k: f64) -> M<N> {
    let l: [f64; N] = array::from_fn(|_| 1.0 + 0.75 * r.next());
    let a = equilibrated(&planted(&haar(r), &l));
    let d: [f64; N] = array::from_fn(|i| 10f64.powf(-k * i as f64 / (2.0 * (N - 1) as f64)));
    array::from_fn(|i| array::from_fn(|j| d[i] * a[i][j] * d[j]))
}

/// Positive definite with a unit diagonal, from a log-spaced spectrum over `10^-k`.
fn unit_diag<const N: usize>(r: &mut Rng, k: f64) -> M<N> {
    let l: [f64; N] = array::from_fn(|i| 10f64.powf(-k * i as f64 / (N - 1) as f64));
    equilibrated(&planted(&haar(r), &l))
}

fn families<const N: usize>() -> Vec<Family<N>> {
    std::vec![
        ("random", false, |r| {
            let mut a = [[0.0; N]; N];
            for i in 0..N {
                for j in 0..=i {
                    a[i][j] = r.next();
                    a[j][i] = a[i][j];
                }
            }
            a
        }),
        ("loguniform", false, |r| {
            let mut a = [[0.0; N]; N];
            for i in 0..N {
                for j in 0..=i {
                    a[i][j] = r.next() * 10f64.powf(-8.0 * r.next().abs());
                    a[j][i] = a[i][j];
                }
            }
            a
        }),
        ("spectrum", false, |r| {
            let l = spectrum::<N>(r);
            planted(&haar(r), &l)
        }),
        ("gap-1e-4/bottom", false, |r| gap(r, 1e-4, false)),
        ("gap-1e-8/top", false, |r| gap(r, 1e-8, true)),
        ("gap-1e-12/bottom", false, |r| gap(r, 1e-12, false)),
        ("gap-1e-16/top", false, |r| gap(r, 1e-16, true)),
        ("cluster-2", false, |r| cluster(r, 2)),
        ("cluster-3", false, |r| cluster(r, 3)),
        ("cluster-half", false, |r| cluster(r, N / 2)),
        ("rank-def-1", false, |r| rank_def(r, 1)),
        ("rank-def-3", false, |r| rank_def(r, 3)),
        ("graded-1e-8", true, |r| graded(r, 8.0)),
        ("graded-1e-16", true, |r| graded(r, 16.0)),
        ("graded-1e-30", true, |r| graded(r, 30.0)),
        ("graded-1e-200", true, |r| graded(r, 200.0)),
        ("unit-diag-1e-6", true, |r| unit_diag(r, 6.0)),
        ("unit-diag-1e-12", true, |r| unit_diag(r, 12.0)),
        ("diag+1e-10", false, |r| {
            let mut a = [[0.0; N]; N];
            for i in 0..N {
                for j in 0..i {
                    a[i][j] = 1e-10 * r.next();
                    a[j][i] = a[i][j];
                }
                a[i][i] = r.next();
            }
            a
        }),
        ("diag-tie+1e-17", false, |r| {
            let mut a = [[0.0; N]; N];
            let x = r.next();
            for i in 0..N {
                for j in 0..i {
                    a[i][j] = 1e-17 * r.next();
                    a[j][i] = a[i][j];
                }
                a[i][i] = if i % 2 == 0 { x } else { r.next() };
            }
            a
        }),
    ]
}

/// The sweeps `eig_sym` runs at `S` on `a` before the polish.
fn sweeps<S: Lane, const N: usize>(a: &M<N>) -> usize {
    jacobi::<S, N, true>(&matrix::<S, N>(a)).3
}

/// `(abs, rel, frame)` of `eig_sym` at `S` on `a` (rounded to `S` first); `rel` 0 unless `pd`.
fn score<S: Lane, const N: usize>(a: &M<N>, pd: bool) -> [f64; 3] {
    let a = round::<S, N>(a);
    let out = run::<S, N>(&a);
    let (t, u) = (reference(&a), unit::<S>());
    // Rounded to `S`, a family may leave the hypothesis: `unit-diag-1e-12` is no longer positive
    // definite at `f32`, and `graded-1e-200` keeps no normal diagonal there.
    let normal = match S::PRECISION {
        Precision::F64 => f64::MIN_POSITIVE,
        Precision::F32 => f64::from(f32::MIN_POSITIVE),
    };
    let pd = pd && (0..N).all(|i| a[i][i] >= normal);
    let k = if pd { kappa(&a) } else { f64::INFINITY };
    let rel = if 0.0 < k && k * u <= 2f64.powi(-10) {
        rel_metric(&t, k, &out, u)
    } else {
        0.0
    };
    [abs_metric(&t, &out, u), rel, frame::<S, N>(&out.1)]
}

/// The worst `(sweeps, abs, rel, frame)` at `S` over `n` matrices of each family, seeded by `seed`;
/// the accuracy on the first `scored` of them.
fn study<S: Lane, const N: usize>(
    n: usize,
    scored: usize,
    seed: u64,
) -> Vec<(&'static str, [f64; 4])> {
    families::<N>()
        .into_iter()
        .enumerate()
        .map(|(i, (name, pd, make))| {
            let mut r = Rng(seed ^ (i as u64) << 32 ^ N as u64);
            let mut w = [0.0_f64; 4];
            for j in 0..n {
                let a = make(&mut r);
                w[0] = w[0].max(sweeps::<S, N>(&round::<S, N>(&a)) as f64);
                if j < scored {
                    let s = score::<S, N>(&a, pd);
                    for (x, y) in w[1..].iter_mut().zip(s) {
                        *x = x.max(y);
                    }
                }
            }
            (name, w)
        })
        .collect()
}

/// Every family at `N`, at `f64` and `f32`, within the recorded bounds and the cap.
fn bounded<const N: usize>(n: usize) {
    for (prec, rows) in [
        ("f64", study::<f64, N>(n, n, 0x6569_6753)),
        ("f32", study::<f32, N>(n, n, 0x6569_6753)),
    ] {
        for (name, [s, a, r, f]) in rows {
            let what = std::format!("N = {N} {prec} {name}");
            assert!(
                s < crate::eig_sym::CAP[N - 4] as f64,
                "{what}: {s} sweeps hit the cap"
            );
            assert!(
                a <= ABS && r <= REL && f <= FRAME,
                "{what}: abs {a} rel {r} frame {f}"
            );
        }
    }
}

#[test]
fn every_family_is_within_its_bounds_at_n_4_9_16() {
    bounded::<4>(40);
    bounded::<9>(20);
    bounded::<16>(6);
}

#[test]
fn every_family_is_within_its_bounds_at_the_other_n() {
    bounded::<5>(6);
    bounded::<6>(6);
    bounded::<7>(6);
    bounded::<8>(6);
    bounded::<10>(4);
    bounded::<11>(4);
    bounded::<12>(4);
    bounded::<13>(3);
    bounded::<14>(3);
    bounded::<15>(3);
}

/// Prints `study` at every `N`: the figures `0111` records and the bounds are set from. Run with
/// `cargo nextest run --release -p helicoid-linalg --run-ignored only measure_eig_sym
/// --no-capture`; `EIG_SYM_N`, `EIG_SYM_SCORED` override the sample counts.
#[test]
#[ignore = "measurement: prints the figures `0111` and the bounds are recorded from"]
#[allow(clippy::print_stdout)]
fn measure_eig_sym() {
    let get = |k: &str, d: usize| {
        std::env::var(k)
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(d)
    };
    let (n, scored) = (get("EIG_SYM_N", 20_000), get("EIG_SYM_SCORED", 500));
    fn at<const N: usize>(n: usize, scored: usize) {
        let (a, b) = (study::<f64, N>(n, scored, 1), study::<f32, N>(n, scored, 1));
        for ((name, x), (_, y)) in a.into_iter().zip(b) {
            std::println!(
                "N={N:2} {name:16} f64 sweeps {:2} abs {:6.3} rel {:6.3} frame {:6.3} | f32 sweeps {:2} abs {:6.3} rel {:6.3} frame {:6.3}",
                x[0], x[1], x[2], x[3], y[0], y[1], y[2], y[3]
            );
        }
    }
    at::<4>(n, scored);
    at::<5>(n, scored);
    at::<6>(n, scored);
    at::<7>(n, scored);
    at::<8>(n, scored);
    at::<9>(n, scored);
    at::<10>(n, scored);
    at::<11>(n, scored);
    at::<12>(n, scored);
    at::<13>(n, scored);
    at::<14>(n, scored);
    at::<15>(n, scored);
    at::<16>(n, scored);
}

/// Equal to the bit, except that a `NaN` matches any `NaN`.
fn same(x: f64, y: f64) -> bool {
    x.to_bits() == y.to_bits() || x.is_nan() && y.is_nan()
}

#[test]
fn a_diagonal_input_returns_its_sorted_diagonal_bit_for_bit() {
    let d = [3.0, -1.0, 2.0, -1.0, 0.0, 5.0, 2.0];
    let a: M<7> = array::from_fn(|i| array::from_fn(|j| if i == j { d[i] } else { 0.0 }));
    let (l, v) = run::<f64, 7>(&a);
    let mut sorted = d;
    sorted.sort_by(f64::total_cmp);
    assert_eq!(l.map(f64::to_bits), sorted.map(f64::to_bits));
    // Stable: the two `-1`s and the two `2`s keep their axes in index order.
    let axes = [1, 3, 4, 2, 6, 0, 5];
    for (k, &e) in axes.iter().enumerate() {
        assert_eq!(
            v[k][e].abs().to_bits(),
            1f64.to_bits(),
            "column {k}: {:?}",
            v[k]
        );
    }
    assert_eq!(det(&v).to_bits(), 1f64.to_bits());
}

#[test]
fn a_multiple_of_the_identity_is_its_own_decomposition() {
    for c in [0.0, 1.0, -2.5, 1e-300, 3e300] {
        let a: M<9> = array::from_fn(|i| array::from_fn(|j| if i == j { c } else { 0.0 }));
        let (l, v) = run::<f64, 9>(&a);
        assert!(l.iter().all(|x| x.to_bits() == c.to_bits()), "{c}: {l:?}");
        for k in 0..9 {
            for r in 0..9 {
                assert_eq!(
                    v[k][r].to_bits(),
                    f64::from(u8::from(k == r)).to_bits(),
                    "{c}"
                );
            }
        }
    }
}

#[test]
fn an_exact_power_of_two_scales_the_values_and_keeps_the_vectors() {
    let mut r = Rng(5);
    for _ in 0..20 {
        let a = families::<8>()[0].2(&mut r);
        let (l, v) = run::<f64, 8>(&a);
        for k in [-600, -40, 3, 500] {
            let s = 2f64.powi(k);
            let (ls, vs) = run::<f64, 8>(&a.map(|row| row.map(|x| x * s)));
            assert!((0..8).all(|i| same(ls[i], l[i] * s)), "{k}");
            assert!(
                (0..8).all(|i| (0..8).all(|j| same(vs[i][j], v[i][j]))),
                "{k}"
            );
        }
    }
}

#[test]
fn a_non_finite_entry_has_no_eigenvalues_and_the_identity_frame() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut a: M<5> = [[0.5; 5]; 5];
        a[3][1] = bad;
        let (l, v) = run::<f64, 5>(&a);
        assert!(l.iter().all(|x| x.is_nan()), "{bad}: {l:?}");
        for k in 0..5 {
            for r in 0..5 {
                assert_eq!(v[k][r].to_bits(), f64::from(u8::from(k == r)).to_bits());
            }
        }
    }
}

#[test]
fn the_upper_triangle_is_ignored() {
    let mut r = Rng(9);
    let a = families::<6>()[0].2(&mut r);
    let mut b = a;
    for i in 0..6 {
        for j in i + 1..6 {
            b[i][j] = 7.0 * r.next();
        }
    }
    let (x, y) = (run::<f64, 6>(&a), run::<f64, 6>(&b));
    assert!((0..6).all(|i| same(x.0[i], y.0[i])));
    assert!((0..6).all(|i| (0..6).all(|j| same(x.1[i][j], y.1[i][j]))));
}

#[test]
fn a_permutation_moves_no_eigenvalue_beyond_its_bound() {
    let mut r = Rng(11);
    let perm = [4, 0, 7, 2, 9, 1, 3, 8, 5, 6];
    for _ in 0..50 {
        let a = families::<10>()[2].2(&mut r);
        let p: M<10> = array::from_fn(|i| array::from_fn(|j| a[perm[i]][perm[j]]));
        let (x, y) = (run::<f64, 10>(&a).0, run::<f64, 10>(&p).0);
        let ln = norm(&x);
        let d = norm::<10>(&array::from_fn(|k| x[k] - y[k])) / (ln * unit::<f64>());
        assert!(d <= 2.0 * ABS, "{d}");
    }
}

/// The relative threshold is what keeps a graded matrix's small eigenvalues: the absolute
/// threshold `|a_pq| <= u |A|` leaves them to `u |A|`, a relative error of `1e16 u` on
/// `graded-1e-16`, while the relative one reads within `REL` (`0111`).
#[test]
fn an_absolute_threshold_loses_the_graded_eigenvalues() {
    let mut r = Rng(13);
    let (mut relative, mut absolute) = (0.0_f64, 0.0_f64);
    for _ in 0..20 {
        let a = graded::<9>(&mut r, 16.0);
        let (t, k, u) = (reference(&a), kappa(&a), unit::<f64>());
        relative = relative.max(rel_metric(&t, k, &run::<f64, 9>(&a), u));
        let (l, v, back, _) = jacobi::<f64, 9, false>(&matrix::<f64, 9>(&a));
        let mut order: [usize; 9] = array::from_fn(|i| i);
        order.sort_by(|&i, &j| l[i].total_cmp(&l[j]));
        let out = (order.map(|i| l[i] * back), order.map(|i| v[i]));
        absolute = absolute.max(rel_metric(&t, k, &out, u));
    }
    assert!(relative <= REL, "{relative}");
    assert!(absolute >= 1e6, "{absolute}");
}

/// The matrices the scalar tests share with the lane and `Dual` tests: every family once, and the
/// special inputs.
fn pool<const N: usize>() -> Vec<M<N>> {
    let mut r = Rng(17);
    let mut rows: Vec<M<N>> = families::<N>().iter().map(|f| f.2(&mut r)).collect();
    let mut nan = [[0.25; N]; N];
    nan[N - 1][0] = f64::NAN;
    rows.extend([[[0.0; N]; N], [[1e-310; N]; N], [[1e300; N]; N], nan]);
    rows
}

/// The value path of `Dual` is the plain value path, bit for bit.
#[test]
fn dual_value_is_plain_value() {
    fn at<const N: usize>() {
        for a in pool::<N>() {
            let x = matrix::<f64, N>(&a);
            let y: Matrix<Dual<f64, 2>, N, N> = Matrix::from_rows(array::from_fn(|i| {
                Vector(array::from_fn(|j| Dual {
                    v: a[i][j],
                    d: [0.5 * a[j][i], f64::from(u8::from(i == j))],
                }))
            }));
            let ((pl, pv), (dl, dv)) = (eig_sym(&x), eig_sym(&y));
            assert!((0..N).all(|k| same(pl.0[k], dl.0[k].v)), "{a:?}");
            assert!(
                (0..N).all(|r| (0..N).all(|c| same(pv.get(r, c), dv.get(r, c).v))),
                "{a:?}"
            );
        }
    }
    at::<4>();
    at::<11>();
}

/// `x` and `y` as the two lanes of one call agree with each run as a scalar, bit for bit.
fn lanes_agree<const N: usize>(x: &M<N>, y: &M<N>) {
    let m: Matrix<L2, N, N> = Matrix::from_rows(array::from_fn(|r| {
        Vector(array::from_fn(|c| L2([x[r][c], y[r][c]])))
    }));
    let (l, v) = eig_sym(&m);
    for (lane, a) in [(0, x), (1, y)] {
        let (wl, wv) = run::<f64, N>(a);
        assert!((0..N).all(|k| same(l.0[k].0[lane], wl[k])), "{a:?}");
        assert!(
            (0..N).all(|k| (0..N).all(|r| same(v.get(r, k).0[lane], wv[k][r]))),
            "{a:?}"
        );
    }
}

/// A lane mask runs every rotation with its lanes held: two lanes agree with the scalars, a
/// degenerate lane beside a regular one included.
#[test]
fn two_lanes_agree_with_scalars() {
    fn at<const N: usize>() {
        let rows = pool::<N>();
        let n = rows.len();
        for (i, j) in (0..n).flat_map(|i| [(i, (i * 7 + 3) % n), (i, i), (i, n - 1 - i)]) {
            lanes_agree(&rows[i], &rows[j]);
        }
    }
    at::<4>();
    at::<9>();
}

/// A held lane keeps the sign of a zero: rotated by `(1, 0)`, the `-0` beside a negative entry
/// would become `+0` and part from the scalar run.
#[test]
fn a_held_lane_keeps_the_sign_of_zero() {
    let x = [
        [1.0, 0.0, -0.5, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [-0.5, 0.0, 0.5, 0.0],
        [0.0, 1.0, 0.0, -0.0],
    ];
    let y = [
        [0.0, -0.0, -0.0, -0.0],
        [-0.0, -0.0, -1.0, 0.0],
        [-0.0, -1.0, -0.0, -0.0],
        [-0.0, 0.0, -0.0, -1.0],
    ];
    lanes_agree(&x, &y);
    lanes_agree(&y, &x);
}

/// `d l_k = v_k^T E v_k` and `d v_k = sum_{j != k} v_j (v_j^T E v_k) / (l_k - l_j)` (first-order
/// perturbation theory) on spectra at least `0.2` apart, to `1e-9`; and on a diagonal input with an
/// off-diagonal lane, where only the polish turns a rotation (`0084` decision 4 at `N`).
#[test]
fn derivative_is_first_order_perturbation_theory() {
    const N: usize = 6;
    let mut r = Rng(19);
    let e: M<N> = {
        let mut e = [[0.0; N]; N];
        for i in 0..N {
            for j in 0..=i {
                e[i][j] = r.next();
                e[j][i] = e[i][j];
            }
        }
        e
    };
    let spectra: [[f64; N]; 2] = [
        [-1.0, -0.6, 0.0, 0.3, 1.1, 2.0],
        [0.5, 1.0, 1.5, 2.0, 2.5, 4.0],
    ];
    let diag = array::from_fn(|i| array::from_fn(|j| if i == j { 1.0 + i as f64 } else { 0.0 }));
    for a in [
        planted(&haar(&mut r), &spectra[0]),
        planted(&haar(&mut r), &spectra[1]),
        diag,
    ] {
        let x: Matrix<Dual<f64, 1>, N, N> = Matrix::from_rows(array::from_fn(|i| {
            Vector(array::from_fn(|j| Dual {
                v: a[i][j],
                d: [e[i][j]],
            }))
        }));
        let (l, v) = eig_sym(&x);
        let vs: [[f64; N]; N] = array::from_fn(|k| array::from_fn(|i| v.get(i, k).v));
        for k in 0..N {
            let ev: [f64; N] = array::from_fn(|i| dot(&e[i], &vs[k]));
            assert!((l.0[k].d[0] - dot(&vs[k], &ev)).abs() <= 1e-9, "{k} {a:?}");
            let dv = (0..N).filter(|&j| j != k).fold([0.0; N], |acc, j| {
                let c = dot(&vs[j], &ev) / (l.0[k].v - l.0[j].v);
                array::from_fn(|i| acc[i] + c * vs[j][i])
            });
            assert!(
                (0..N).all(|i| (v.get(i, k).d[0] - dv[i]).abs() <= 1e-9),
                "{k} {a:?}"
            );
        }
    }
}

/// A hill climb on the entries for the most sweeps, at `N = 4, 9, 16` and `f64`: from each
/// family's matrices, perturbations of a shrinking size, kept where the sweeps do not fall. Prints
/// the most found, beside `study`'s, for `0111`'s cap. `EIG_SYM_RESTARTS` overrides the restarts
/// per family.
#[test]
#[ignore = "measurement: the adversarial half of the convergence study (`0111`)"]
#[allow(clippy::print_stdout)]
fn hill_climb_eig_sym() {
    let restarts = std::env::var("EIG_SYM_RESTARTS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(8);
    fn at<const N: usize>(restarts: usize) {
        let mut best = (0, "");
        for (i, (name, _, make)) in families::<N>().into_iter().enumerate() {
            let mut r = Rng(0x6869_6c6c ^ (i as u64) << 20 ^ N as u64);
            for _ in 0..restarts {
                let mut a = make(&mut r);
                let mut s = sweeps::<f64, N>(&a);
                for step in 0..300 {
                    let size = 10f64.powf(-1.0 - 15.0 * step as f64 / 300.0);
                    let mut b = a;
                    for p in 0..N {
                        for q in 0..=p {
                            let x = b[p][q] + size * r.next() * (b[p][q].abs() + b[p][p].abs());
                            b[p][q] = x;
                            b[q][p] = x;
                        }
                    }
                    let t = sweeps::<f64, N>(&b);
                    if t >= s {
                        (a, s) = (b, t);
                    }
                }
                if s > best.0 {
                    best = (s, name);
                }
            }
        }
        std::println!("N={N:2} hill climb: {} sweeps (from {})", best.0, best.1);
    }
    at::<4>(restarts);
    at::<9>(restarts);
    at::<16>(restarts);
}

/// A positive definite block far below `|A|`, where `angle` unscaled takes its safe arm: `2^-70`
/// at `f32`, `2^-600` at `f64`, keeps its digits.
#[test]
fn a_block_below_the_safe_arm_keeps_its_digits() {
    fn at<S: Lane>(e: f64) {
        let mut a: M<4> = [[0.0; 4]; 4];
        (a[0][0], a[1][1]) = (1.0, 1.0);
        (a[2][2], a[3][3]) = (e, e);
        (a[3][2], a[2][3]) = (0.5 * e, 0.5 * e);
        let a = round::<S, 4>(&a);
        let u = unit::<S>();
        let w = rel_metric(&reference(&a), kappa(&a), &run::<S, 4>(&a), u);
        assert!(w <= REL, "{e}: {w}");
    }
    at::<f32>(2f64.powi(-70));
    at::<f64>(2f64.powi(-600));
}
