//! The coefficient kernel's latency, by stratum (`docs/PHASE1.md` §9, `docs/PHASE3.md` §11).
//!
//! §11 owes `exp` at near-identity `θ` from `7.5e-8` to `1` first, where the kernel runs one `sqrt`
//! and one `sin_cos`; the whole shipped kernel is benched here because `Exp` does not exist yet and
//! the kernel is what it will call. The subject is the code that ships, reached through the hidden
//! `__sweep` feature, as the conformance subject reaches it (`0004` item 4).
//!
//! **Fixtures are by stratum**, §9's three plus the two switches the kernel branches on, because a
//! group's cost is which arm ran: `a` takes its series below `θ² = 0.723`, `b` below `0.965`
//! (`coeffs/generated.rs`), so `θ = 0.9` is the one point where `jr_coeffs` runs one exact and one
//! series member. The branch variable is `z = θ²` for every group but `log_ratio`, whose is
//! `s = n²/w²`.
//!
//! Timing is a function of the host, not of the output bits, so nothing here is a D16 claim and
//! nothing is compared against a reference. The gate is `cargo xtask bench-gate`, whose `--aa`
//! floor says how much of a ratio is this machine rather than the code.

// `criterion_main!` generates an undocumented `main`; the workspace warns on `missing_docs` and
// `just lint` denies warnings. The attribute is the macro's, not this file's code.
#![allow(
    missing_docs,
    reason = "criterion_main! generates an undocumented `main`"
)]

use core::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use helicoid::__sweep as k;

/// `(name, θ)`: §9's three strata, then the two switch neighbourhoods.
const THETA: &[(&str, f64)] = &[
    // §11's near-identity end: below every switch, so every group is on its series arm and `Exp`
    // runs one `sqrt` and one `sin_cos`.
    ("near-identity-7.5e-8", 7.5e-8),
    ("small-1e-4", 1e-4),
    ("series-0.5", 0.5),
    // At binary64 `a` is exact here and `b` is not, the one fixture where `jr_coeffs`
    // straddles its group. At binary32 it is not a straddle: `A_F32` and `B_F32` both switch
    // at `θ² < 1`, so `z = 0.81` is on both series arms and no exact arm runs. Measured by
    // negative control — slowing `exact_a` moves the binary64 row 1.38x and the binary32 row
    // not at all.
    ("straddle-0.9", 0.9),
    ("generic-1", 1.0),
    ("near-pi", core::f64::consts::PI - 1e-6),
];

/// `z = fl(θ·θ)`, the branch variable the kernel takes, formed at the bench's precision.
fn z_f64(theta: f64) -> f64 {
    theta * theta
}

fn z_f32(theta: f64) -> f32 {
    let t = theta as f32;
    t * t
}

/// `(n², w)` of the unit quaternion at `θ`: `log_ratio`'s arguments.
fn nw_f64(theta: f64) -> (f64, f64) {
    let (s, c) = (theta * 0.5).sin_cos();
    (s * s, c)
}

fn nw_f32(theta: f64) -> (f32, f32) {
    let (n2, w) = nw_f64(theta);
    (n2 as f32, w as f32)
}

fn groups(c: &mut Criterion) {
    let mut g = c.benchmark_group("coeffs/f64");
    for &(name, theta) in THETA {
        let z = z_f64(theta);
        let (n2, w) = nw_f64(theta);
        g.bench_function(format!("exp_coeffs/{name}"), |b| {
            b.iter(|| k::exp_coeffs(black_box(z)));
        });
        g.bench_function(format!("jr_coeffs/{name}"), |b| {
            b.iter(|| k::jr_coeffs(black_box(z)));
        });
        g.bench_function(format!("jr_inv_coeff/{name}"), |b| {
            b.iter(|| k::jr_inv_coeff(black_box(z)));
        });
        g.bench_function(format!("q_coeffs/{name}"), |b| {
            b.iter(|| k::q_coeffs(black_box(z)));
        });
        g.bench_function(format!("log_ratio/{name}"), |b| {
            b.iter(|| k::log_ratio(black_box(n2), black_box(w)));
        });
    }
    g.finish();

    let mut g = c.benchmark_group("coeffs/f32");
    for &(name, theta) in THETA {
        let z = z_f32(theta);
        let (n2, w) = nw_f32(theta);
        g.bench_function(format!("exp_coeffs/{name}"), |b| {
            b.iter(|| k::exp_coeffs(black_box(z)));
        });
        g.bench_function(format!("jr_coeffs/{name}"), |b| {
            b.iter(|| k::jr_coeffs(black_box(z)));
        });
        g.bench_function(format!("jr_inv_coeff/{name}"), |b| {
            b.iter(|| k::jr_inv_coeff(black_box(z)));
        });
        g.bench_function(format!("q_coeffs/{name}"), |b| {
            b.iter(|| k::q_coeffs(black_box(z)));
        });
        g.bench_function(format!("log_ratio/{name}"), |b| {
            b.iter(|| k::log_ratio(black_box(n2), black_box(w)));
        });
    }
    g.finish();
}

// A gate has to be runnable: criterion's defaults (3 s warm-up, 5 s measurement) are 16 minutes
// over two passes of these 60 benchmarks. 0.5 s and 1.5 s keep 100 samples each and put a pass at
// about two minutes; whether that is enough is not asserted here, it is *measured* — the A/A floor
// `cargo xtask bench-gate --aa` reports is exactly the cost of this choice.
criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_millis(1500))
        .sample_size(100);
    targets = groups
}
criterion_main!(benches);
