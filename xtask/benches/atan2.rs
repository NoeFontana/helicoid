//! The portable `atan2` candidate of `0105` (draft) against `libm` and the host's `atan2`, and
//! `so3_log` end to end at the shipped scalar and at the twin's: the time per call that draft
//! record cites. `cargo bench -p xtask --bench atan2`.
//!
//! Rows are batches of 256 calls: `latency` feeds each result into the next call's `y` through
//! `r·0`, `throughput` sums independent calls. The regimes are random octants (`mixed`), one octant
//! (`octant`), and `y = 0` (`axis`, where `libm` returns early).

// `criterion_main!` generates an undocumented `main`; the workspace warns on `missing_docs` and
// clippy runs with `-D warnings`.
#![allow(
    missing_docs,
    reason = "criterion_main! generates an undocumented `main`"
)]

use core::hint::black_box;
use core::time::Duration;

use criterion::measurement::WallTime;
use criterion::{criterion_group, criterion_main, BenchmarkGroup, Criterion, Throughput};
use helicoid::{LieGroup, SO3Tangent, SO3};
use helicoid_linalg::{Real, Vector};

#[allow(
    dead_code,
    reason = "the twin's module, of which the bench reads the kernels"
)]
#[path = "../src/shipped/portable_atan2/kernel.rs"]
mod kernel;

use kernel::{atan2_f32, atan2_f64, fits, Portable};

const BATCH: usize = 256;

/// splitmix64, so a batch is the same on every host.
fn draws(seed: u64) -> impl FnMut() -> f64 {
    let mut state = seed;
    move || {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        // `±m·2^e`, `m` in `[1, 2)`, `e` in `[-2, 2]`.
        let m = 1.0 + (z >> 11) as f64 / (1u64 << 53) as f64;
        let x = libm::ldexp(m, (z % 5) as i32 - 2);
        if z & 1 == 0 {
            x
        } else {
            -x
        }
    }
}

fn batch(regime: &str) -> (Vec<f64>, Vec<f64>) {
    let mut draw = draws(9);
    (0..BATCH)
        .map(|_| {
            let (y, x) = (draw(), draw());
            match regime {
                "octant" => (y.abs().min(x.abs()), y.abs().max(x.abs())),
                "axis" => (0.0, x.abs()),
                _ => (y, x),
            }
        })
        .unzip()
}

/// One subject's latency and throughput rows in one regime.
fn rows<T, F>(g: &mut BenchmarkGroup<'_, WallTime>, label: &str, ys: &[T], xs: &[T], f: F)
where
    T: Real,
    F: Fn(T, T) -> T,
{
    g.bench_function(format!("{label}/latency"), |b| {
        b.iter(|| {
            let mut r = T::zero();
            for (&y, &x) in ys.iter().zip(xs) {
                r = f(y + r * T::zero(), x);
            }
            black_box(r)
        });
    });
    g.bench_function(format!("{label}/throughput"), |b| {
        b.iter(|| {
            let mut acc = T::zero();
            for (&y, &x) in black_box(ys).iter().zip(black_box(xs)) {
                acc = acc + f(y, x);
            }
            black_box(acc)
        });
    });
}

fn group<'a>(c: &'a mut Criterion, name: &str) -> BenchmarkGroup<'a, WallTime> {
    let mut g = c.benchmark_group(name);
    g.throughput(Throughput::Elements(BATCH as u64));
    g.warm_up_time(Duration::from_millis(500));
    g.measurement_time(Duration::from_secs(2));
    g
}

fn atan2(c: &mut Criterion) {
    let Ok(fits) = fits() else {
        return;
    };
    let (f64_fit, f32_fit) = (&fits.f64, &fits.f32);
    let pi4 = (f64_fit.pi4_hi, f64_fit.pi4_lo);
    let mut g = group(c, "atan2/f64");
    for regime in ["mixed", "octant", "axis"] {
        let (ys, xs) = batch(regime);
        rows(&mut g, &format!("libm/{regime}"), &ys, &xs, libm::atan2);
        rows(&mut g, &format!("candidate/{regime}"), &ys, &xs, |y, x| {
            atan2_f64(y, x, f64_fit)
        });
        rows(&mut g, &format!("host/{regime}"), &ys, &xs, f64::atan2);
    }
    g.finish();
    let mut g = group(c, "atan2/f32");
    for regime in ["mixed", "octant", "axis"] {
        let (ys, xs) = batch(regime);
        let ys: Vec<f32> = ys.iter().map(|&v| v as f32).collect();
        let xs: Vec<f32> = xs.iter().map(|&v| v as f32).collect();
        rows(&mut g, &format!("libm/{regime}"), &ys, &xs, libm::atan2f);
        rows(&mut g, &format!("candidate/{regime}"), &ys, &xs, |y, x| {
            atan2_f32(y, x, f32_fit, pi4)
        });
        rows(&mut g, &format!("host/{regime}"), &ys, &xs, f32::atan2);
    }
    g.finish();
}

/// `BATCH` rotations at `θ` in `[lo, hi)` about random axes.
fn rotations<S: Real>(lo: f64, hi: f64) -> Vec<SO3<S>> {
    let mut draw = draws(3);
    (0..BATCH)
        .map(|_| {
            let axis = [draw(), draw(), draw()];
            let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
            let theta = lo + (hi - lo) * (draw().abs() - 0.25) / 3.75;
            SO3::exp(&SO3Tangent {
                phi: Vector(axis.map(|a| S::round_from_f64(a / norm * theta))),
            })
        })
        .collect()
}

fn so3_log(c: &mut Criterion) {
    let mut g = group(c, "so3_log");
    // `r`'s exact arm calls `atan2` (`θ ≳ 0.58`, `0090`); below it a series arm does not.
    for (arm, lo, hi) in [("exact", 0.6, 3.1), ("series", 0.01, 0.5)] {
        let (a, b) = (rotations::<f64>(lo, hi), rotations::<Portable<f64>>(lo, hi));
        g.bench_function(format!("f64/libm/{arm}"), |bench| {
            bench.iter(|| black_box(&a).iter().map(|x| x.log().phi.0[0]).sum::<f64>());
        });
        g.bench_function(format!("f64/candidate/{arm}"), |bench| {
            bench.iter(|| {
                black_box(&b)
                    .iter()
                    .map(|x| x.log().phi.0[0].0)
                    .sum::<f64>()
            });
        });
    }
    g.finish();
}

criterion_group!(benches, atan2, so3_log);
criterion_main!(benches);
