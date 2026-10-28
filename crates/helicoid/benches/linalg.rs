//! `helicoid-linalg`'s solvers by stratum (`docs/PHASE1.md` §9): `solve_cubic` (a closed form) and
//! `eig3` (Jacobi with an exit, `0084`) with its values alone, `eig3_values` (`0104`), at `f64` and
//! `f32` (`0053` *Further work* 2), and
//! `chol` at the three clauses of `0021`'s bar (`0080`): positive definite and failing input at
//! `f64` and `f32`, and `Dual`; `eig_sym` (`0111`) at `N = 4, 9, 16` on a random spectrum, a
//! graded positive definite matrix and a triple eigenvalue, at `f64` and `f32`.
//!
//! **A third binary**, for the reason `benches/groups.rs` gives: `bench-gate --against <binary>`
//! takes one bench binary, so a change to a solver pays for these rows alone.
//!
//! Each routine is timed behind an `#[inline(never)]` call, which pins one compilation of it for
//! baseline and candidate alike (`benches/groups.rs`, `geodesic_of`). Fixtures are one input per
//! stratum of the corpus ids (`0056`), so a row's arm is the stratum's: `solve_cubic` at three
//! distinct roots (the trigonometric arm), one real root, a small root beside a complex pair, and a
//! triple root; `eig3` and `eig3_values` at a random spectrum, a pair `1e-8` apart, and a triple
//! eigenvalue, and the three in turn; `svd3` (one-sided Jacobi with an exit, `0101`) at generic
//! singular values, a pair `1e-8` apart, a rotation `1e-8` off SO(3) and a rank-2 matrix, and the
//! four in turn; `nearest_rotation` (`0102`) at one input per arm: Newton near SO(3) and at
//! generic singular values, the reflection, and `svd3`'s at a rank-2 matrix.
//!
//! A baseline binary for `--against` is built in this tree with only the bodies under test
//! reverted, as `benches/groups.rs` says.

// `criterion_main!` generates an undocumented `main`; the workspace warns on `missing_docs` and
// clippy runs with `-D warnings`.
#![allow(
    missing_docs,
    reason = "criterion_main! generates an undocumented `main`"
)]

use core::hint::black_box;
use core::time::Duration;

use criterion::measurement::WallTime;
use criterion::{criterion_group, criterion_main, BenchmarkGroup, Criterion};
use helicoid_linalg::{
    chol, eig3, eig3_values, eig_sym, nearest_rotation, solve_cubic, svd3, Dual, Mat3, Matrix,
    Real, Vec3, Vector,
};

/// Four scalars, not `[S; 4]`: an array by value reaches the callee through the stack, and a wide
/// load of two of its narrow stores is a store-to-load forwarding stall that read 10-20 ns of a
/// 14-45 ns routine and moved with the routine's own load order. A caller of `solve_cubic` holds
/// its coefficients in registers.
#[inline(never)]
fn cubic_of<S: Real>(a: S, b: S, c: S, d: S) -> (Vec3<S>, [S::Mask; 3]) {
    solve_cubic(a, b, c, d)
}

#[inline(never)]
fn eig3_of<S: Real>(a: &Mat3<S>) -> (Vec3<S>, Mat3<S>) {
    eig3(a)
}

#[inline(never)]
fn eig3_values_of<S: Real>(a: &Mat3<S>) -> Vec3<S> {
    eig3_values(a)
}

#[inline(never)]
fn svd3_of<S: Real>(a: &Mat3<S>) -> (Mat3<S>, Vec3<S>, Mat3<S>) {
    svd3(a)
}

#[inline(never)]
fn nearest_rotation_of<S: Real>(a: &Mat3<S>) -> Mat3<S> {
    nearest_rotation(a)
}

#[inline(never)]
fn eig_sym_of<S: Real, const N: usize>(a: &Matrix<S, N, N>) -> (Vector<S, N>, Matrix<S, N, N>) {
    eig_sym(a)
}

/// The whole factor and the mask come back, and criterion's `black_box` keeps every entry live:
/// consuming one entry lets the optimizer delete most of the factorisation (`0030`).
#[inline(never)]
fn chol_of<S: Real, const N: usize>(a: &Matrix<S, N, N>) -> (Matrix<S, N, N>, S::Mask) {
    chol(a)
}

/// A fixture at the precision: the `f64` value, rounded (`benches/groups.rs`, `Fixture`).
trait Fixture: Real {
    fn of(x: f64) -> Self;
}

impl Fixture for f64 {
    fn of(x: f64) -> Self {
        x
    }
}

impl Fixture for f32 {
    fn of(x: f64) -> Self {
        x as f32
    }
}

/// `(name, a, b, c, d)`: each arm of `solve_cubic`, the one-root arm at both signs of `p`, and a
/// cubic outside the homogenisation window (`(x - 1e40)(x - 2e40)(x - 3e40)` at `f64`, `1e6` at
/// `f32`), whose cold ladder the others never reach.
fn cubics<S: Real>() -> [(&'static str, [f64; 4]); 6] {
    let far = match S::PRECISION {
        helicoid_linalg::Precision::F64 => 1e40,
        helicoid_linalg::Precision::F32 => 1e6,
    };
    [
        ("distinct", [1.0, -6.0, 11.0, -6.0]),
        ("one-real", [1.0, 1.0, 1.0, 1.0]),
        ("one-real-p-negative", [1.0, 0.0, -3.0, 3.0]),
        ("one-real-small-root", [1.0, 0.0, 1.0, 1e-6]),
        ("triple", [1.0, -6.0, 12.0, -8.0]),
        (
            "out-of-window",
            [1.0, -6.0 * far, 11.0 * far * far, -6.0 * far * far * far],
        ),
    ]
}

/// `Q diag(λ) Qᵀ`, `Q` the rotation by `1` rad about `(1, 2, 2)/3`, rounded once at `f64`.
fn sym<S: Fixture>(lambda: [f64; 3]) -> Mat3<S> {
    let (s, c) = (1.0f64.sin(), 1.0f64.cos());
    let n = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
    let q = |i: usize, j: usize| {
        let cross = match (i, j) {
            (0, 1) => -n[2],
            (1, 0) => n[2],
            (0, 2) => n[1],
            (2, 0) => -n[1],
            (1, 2) => -n[0],
            (2, 1) => n[0],
            _ => 0.0,
        };
        let delta = if i == j { 1.0 } else { 0.0 };
        c * delta + (1.0 - c) * n[i] * n[j] + s * cross
    };
    let a = |i: usize, j: usize| (0..3).map(|k| q(i, k) * lambda[k] * q(j, k)).sum::<f64>();
    Matrix::from_cols(core::array::from_fn(|j| {
        Vector(core::array::from_fn(|i| S::of(a(i.max(j), i.min(j)))))
    }))
}

/// The rotation by `t` rad about the unit `n`, in `f64`.
fn rotation(t: f64, n: [f64; 3]) -> [[f64; 3]; 3] {
    let (s, c) = (t.sin(), t.cos());
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let cross = match (i, j) {
                (0, 1) => -n[2],
                (1, 0) => n[2],
                (0, 2) => n[1],
                (2, 0) => -n[1],
                (1, 2) => -n[0],
                (2, 1) => n[0],
                _ => 0.0,
            };
            let delta = if i == j { 1.0 } else { 0.0 };
            c * delta + (1.0 - c) * n[i] * n[j] + s * cross
        })
    })
}

/// `U diag(σ) Vᵀ + ε E`, `U`, `V` fixed rotations and `E` a fixed matrix of entries in `[-1, 1]`,
/// rounded once at `f64`.
fn general<S: Fixture>(sigma: [f64; 3], eps: f64) -> Mat3<S> {
    let u = rotation(1.0, [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]);
    let v = rotation(2.0, [0.6, 0.0, -0.8]);
    let e = [[0.3, -0.9, 0.1], [0.7, 0.2, -0.4], [-0.5, 0.8, 0.6]];
    let a = |i: usize, j: usize| {
        (0..3).map(|k| u[i][k] * sigma[k] * v[j][k]).sum::<f64>() + eps * e[i][j]
    };
    Matrix::from_cols(core::array::from_fn(|j| {
        Vector(core::array::from_fn(|i| S::of(a(i, j))))
    }))
}

const SINGULAR: &[(&str, [f64; 3], f64)] = &[
    ("random", [0.9, 0.5, -0.2], 0.0),
    ("gap-1e-8", [0.9, 0.4 + 1e-8, 0.4], 0.0),
    ("near-rotation-1e-8", [1.0, 1.0, 1.0], 1e-8),
    ("rank2", [0.9, 0.5, 0.0], 0.0),
];

/// One input per arm of `nearest_rotation`: `det` above `tau |M|^3` (Newton), below `-tau |M|^3`
/// (the reflection), and between (`svd3`'s `U V^T`).
const NEAREST: &[(&str, [f64; 3], f64)] = &[
    ("near-rotation-1e-8", [1.0, 1.0, 1.0], 1e-8),
    ("random", [0.9, 0.5, 0.2], 0.0),
    ("reflection", [1.0, 0.6, -0.2], 0.0),
    ("rank2", [0.9, 0.5, 0.0], 0.0),
];

impl Fixture for Dual<f64, 2> {
    fn of(x: f64) -> Self {
        Dual { v: x, d: [x, 0.5] }
    }
}

/// `B^T B + N I`, `B` entries in `[-1, 1)` from a fixed LCG, so the mask is set; `failing` negates
/// the middle pivot, so column `N / 2` runs the failed arm (unit pivot, zeros below) and the later
/// columns factor the matrix with that row and column removed.
fn spd<S: Fixture, const N: usize>(failing: bool) -> Matrix<S, N, N> {
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
    };
    let b: [[f64; N]; N] = core::array::from_fn(|_| core::array::from_fn(|_| next()));
    let a = |i: usize, j: usize| {
        let g = (0..N).map(|k| b[k][i] * b[k][j]).sum::<f64>();
        let shift = if i == j { N as f64 } else { 0.0 };
        let flip = if failing && i == N / 2 && j == N / 2 {
            -1.0
        } else {
            1.0
        };
        flip * (g + shift)
    };
    Matrix::from_cols(core::array::from_fn(|j| {
        Vector(core::array::from_fn(|i| S::of(a(i, j))))
    }))
}

fn chol_row<S: Fixture, const N: usize>(g: &mut BenchmarkGroup<'_, WallTime>) {
    for (name, failing) in [("pd", false), ("failing", true)] {
        let a = spd::<S, N>(failing);
        g.bench_function(format!("chol/{N}/{name}"), |b| {
            b.iter(|| chol_of(black_box(&a)));
        });
    }
}

/// `eig_sym`'s fixtures at `N`, `Q diag(l) Q^T` with `Q` two Householder reflections from a fixed
/// LCG: `random` (`l` uniform in `[-1, 1]`), `graded` (`D Q diag(l) Q^T D`, `l` in `[1/2, 3/2]`,
/// `D` geometric from 1 to `1e-4`, so the diagonal spans `1e-8` and `f32` keeps it normal), and
/// `cluster-3` (three equal eigenvalues).
fn sym_n<S: Fixture, const N: usize>(kind: &str) -> Matrix<S, N, N> {
    let mut state = 0x6569_6753_7973_6d00_u64;
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
    };
    let mut q: [[f64; N]; N] =
        core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
    for _ in 0..2 {
        let v: [f64; N] = core::array::from_fn(|_| next());
        let n2 = v.iter().map(|x| x * x).sum::<f64>();
        let qv: [f64; N] = core::array::from_fn(|i| (0..N).map(|k| q[i][k] * v[k]).sum());
        for i in 0..N {
            for j in 0..N {
                q[i][j] -= 2.0 * qv[i] * v[j] / n2;
            }
        }
    }
    let mut l: [f64; N] = core::array::from_fn(|_| next());
    let mut d = [1.0; N];
    match kind {
        "graded" => {
            l = l.map(|x| 1.0 + 0.5 * x);
            d = core::array::from_fn(|i| 10f64.powf(-4.0 * i as f64 / (N - 1) as f64));
        }
        "cluster-3" => {
            let c = l[0];
            l[1] = c;
            l[2] = c;
        }
        _ => {}
    }
    let a =
        |i: usize, j: usize| d[i] * (0..N).map(|k| q[i][k] * l[k] * q[j][k]).sum::<f64>() * d[j];
    Matrix::from_cols(core::array::from_fn(|j| {
        Vector(core::array::from_fn(|i| S::of(a(i.max(j), i.min(j)))))
    }))
}

fn eig_sym_row<S: Fixture, const N: usize>(g: &mut BenchmarkGroup<'_, WallTime>) {
    for kind in ["random", "graded", "cluster-3"] {
        let a = sym_n::<S, N>(kind);
        g.bench_function(format!("eig_sym/{N}/{kind}"), |b| {
            b.iter(|| eig_sym_of(black_box(&a)));
        });
    }
}

const SPECTRA: &[(&str, [f64; 3])] = &[
    ("random", [-0.7, 0.2, 0.9]),
    ("gap-1e-8", [-0.7, 0.4, 0.4 + 1e-8]),
    ("triple", [0.3, 0.3, 0.3]),
];

fn rows<S: Fixture>(g: &mut BenchmarkGroup<'_, WallTime>) {
    let all = cubics::<S>().map(|(name, c)| (name, c.map(S::of)));
    for (name, [a, b, c, d]) in all {
        g.bench_function(format!("solve_cubic/{name}"), |bench| {
            bench.iter(|| cubic_of(black_box(a), black_box(b), black_box(c), black_box(d)));
        });
    }
    // Every fixture in turn, so a branch on the arm or on the sign of `p` is taken both ways: the
    // cost a caller with mixed cubics pays, six calls a row, that the rows above predict perfectly.
    g.bench_function("solve_cubic/mixed-six-calls", |bench| {
        bench.iter(|| {
            for (_, [a, b, c, d]) in black_box(&all) {
                black_box(cubic_of(*a, *b, *c, *d));
            }
        });
    });
    for &(name, lambda) in SPECTRA {
        let a = sym::<S>(lambda);
        g.bench_function(format!("eig3/{name}"), |b| {
            b.iter(|| eig3_of(black_box(&a)));
        });
        g.bench_function(format!("eig3_values/{name}"), |b| {
            b.iter(|| eig3_values_of(black_box(&a)));
        });
    }
    for &(name, sigma, eps) in SINGULAR {
        let a = general::<S>(sigma, eps);
        g.bench_function(format!("svd3/{name}"), |b| {
            b.iter(|| svd3_of(black_box(&a)));
        });
    }
    for &(name, sigma, eps) in NEAREST {
        let a = general::<S>(sigma, eps);
        g.bench_function(format!("nearest_rotation/{name}"), |b| {
            b.iter(|| nearest_rotation_of(black_box(&a)));
        });
    }
    chol_row::<S, 3>(g);
    chol_row::<S, 6>(g);
    chol_row::<S, 8>(g);
    chol_row::<S, 9>(g);
    chol_row::<S, 12>(g);
    eig_sym_row::<S, 4>(g);
    eig_sym_row::<S, 9>(g);
    eig_sym_row::<S, 16>(g);
    // The three fixtures in turn, as `solve_cubic/mixed-six-calls`: a data-dependent exit is
    // predicted perfectly by the rows above and not here.
    let mixed: Vec<Mat3<S>> = SPECTRA.iter().map(|&(_, l)| sym::<S>(l)).collect();
    g.bench_function("eig3/mixed-three-calls", |b| {
        b.iter(|| {
            for a in black_box(&mixed) {
                black_box(eig3_of(a));
            }
        });
    });
    g.bench_function("eig3_values/mixed-three-calls", |b| {
        b.iter(|| {
            for a in black_box(&mixed) {
                black_box(eig3_values_of(a));
            }
        });
    });
    let mixed: Vec<Mat3<S>> = SINGULAR
        .iter()
        .map(|&(_, s, e)| general::<S>(s, e))
        .collect();
    g.bench_function("svd3/mixed-four-calls", |b| {
        b.iter(|| {
            for a in black_box(&mixed) {
                black_box(svd3_of(a));
            }
        });
    });
}

fn linalg(c: &mut Criterion) {
    let mut g = c.benchmark_group("linalg/f64");
    rows::<f64>(&mut g);
    g.finish();
    let mut g = c.benchmark_group("linalg/f32");
    rows::<f32>(&mut g);
    g.finish();
    let mut g = c.benchmark_group("linalg/dual2");
    chol_row::<Dual<f64, 2>, 6>(&mut g);
    g.finish();
}

// The time budget of `benches/coeffs.rs`, for its reasons: 100 samples, about 2 s per benchmark.
criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_millis(1500))
        .sample_size(100);
    targets = linalg
}
criterion_main!(benches);
