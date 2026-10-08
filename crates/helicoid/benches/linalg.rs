//! `helicoid-linalg`'s closed forms by stratum (`docs/PHASE1.md` §9): `solve_cubic` and `eig3`, at
//! `f64` and `f32` (`0053` *Further work* 2).
//!
//! **A third binary**, for the reason `benches/groups.rs` gives: `bench-gate --against <binary>`
//! takes one bench binary, so a change to a solver pays for these rows alone.
//!
//! Each routine is timed behind an `#[inline(never)]` call, which pins one compilation of it for
//! baseline and candidate alike (`benches/groups.rs`, `geodesic_of`). Fixtures are one input per
//! stratum of the corpus ids (`0056`), so a row's arm is the stratum's: `solve_cubic` at three
//! distinct roots (the trigonometric arm), one real root, a small root beside a complex pair, and a
//! triple root; `eig3` at a random spectrum, a pair `1e-8` apart, and a triple eigenvalue.
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
use helicoid_linalg::{eig3, solve_cubic, Mat3, Matrix, Real, Vec3, Vector};

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

/// `(name, a, b, c, d)`.
const CUBICS: &[(&str, [f64; 4])] = &[
    ("distinct", [1.0, -6.0, 11.0, -6.0]),
    ("one-real", [1.0, 1.0, 1.0, 1.0]),
    ("one-real-small-root", [1.0, 0.0, 1.0, 1e-6]),
    ("triple", [1.0, -6.0, 12.0, -8.0]),
];

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

const SPECTRA: &[(&str, [f64; 3])] = &[
    ("random", [-0.7, 0.2, 0.9]),
    ("gap-1e-8", [-0.7, 0.4, 0.4 + 1e-8]),
    ("triple", [0.3, 0.3, 0.3]),
];

fn rows<S: Fixture>(g: &mut BenchmarkGroup<'_, WallTime>) {
    for &(name, c) in CUBICS {
        let [a, b, c, d] = c.map(S::of);
        g.bench_function(format!("solve_cubic/{name}"), |bench| {
            bench.iter(|| cubic_of(black_box(a), black_box(b), black_box(c), black_box(d)));
        });
    }
    for &(name, lambda) in SPECTRA {
        let a = sym::<S>(lambda);
        g.bench_function(format!("eig3/{name}"), |b| {
            b.iter(|| eig3_of(black_box(&a)));
        });
    }
}

fn linalg(c: &mut Criterion) {
    let mut g = c.benchmark_group("linalg/f64");
    rows::<f64>(&mut g);
    g.finish();
    let mut g = c.benchmark_group("linalg/f32");
    rows::<f32>(&mut g);
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
