//! The groups' latency, by stratum (`docs/PHASE1.md` §9, `docs/PHASE3.md` §11).
//!
//! §11's list — `exp`, `log`, `compose`, `act`, `act_many`, `jr`, `jr_inv`, `adjoint`, `Jac::mul`
//! for SO(3), SE(3) and SE₂(3), per stratum, at `f64` and `f32` — plus `rminus_jacobians`, which
//! §11 does not list and which `SO3::rminus_jacobians`'s own rustdoc calls the crate's hottest
//! Jacobian path: one per residual per solver iteration. Two deferred optimisations are decided by
//! these rows (the grouped-coefficient fusion of `jr`/`jr_inv`, and fusing the two inversions of
//! `rminus_jacobians`), so the rows they would move are benched before either is attempted.
//!
//! **A second binary, not more rows in `coeffs`.** `bench-gate --against <binary>` takes one bench
//! binary, so a target of its own is what keeps a change to the groups from paying for the
//! coefficient kernel's 60 benchmarks and the other way round.
//!
//! Fixtures are the three strata of §9 — near-identity, generic, near-π — with `|φ| = θ` along
//! `(1, 2, 2)/3` and the translation columns fixed, so a row differs from its neighbour in `θ`
//! alone. A `⊖` row is the one that takes care: `x.rminus_jacobians(&base)` evaluates at
//! `τ = Log(base⁻¹ x)`, so two elements of the labelled `θ` about *different* axes would put it on
//! another stratum entirely — at `θ = π − 1e-6` about axes 63.6° apart, `τ` comes out near
//! 2.22 rad and the near-π arm is never reached. Those rows take `base = x · Exp(−τ)`, which makes
//! `Log(base⁻¹ x)` the labelled tangent itself. The `f32` fixture is the `f64` one rounded, which
//! is what a caller at `f32` holds;
//! `Real::lit` is not used, because `7.5e-8` and `π − 1e-6` are no binary32 and `lit` states its
//! argument is exactly representable.
//!
//! A baseline binary for `--against` is built **in this tree**, with only the bodies under test
//! reverted. Built in a separate git worktree instead, identical `se3/jr` code read 0.97 — a build
//! directory alone moves a row 2–3%, which is the layout effect `xtask::bench`'s docs warn about.
//!
//! Timing is a function of the host, not of the output bits, so nothing here is a D16 claim. The
//! gate is `cargo xtask bench-gate --against`, whose bracketed A/A control says how much of a ratio
//! is this machine rather than the code (`0033`, draft).

// `criterion_main!` generates an undocumented `main`; the workspace warns on `missing_docs` and
// `just lint` denies warnings. The attribute is the macro's, not this file's code.
#![allow(
    missing_docs,
    reason = "criterion_main! generates an undocumented `main`"
)]

use core::hint::black_box;
use core::time::Duration;

use criterion::measurement::WallTime;
use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkGroup, Criterion};
use helicoid::{Jac, LieGroup, SEn3Tangent, SO3Tangent, Tangent, SE3, SO3};
use helicoid_linalg::{Point, Point3, Real, Vec3, Vector};

/// `(name, θ)`: §9's three strata.
const THETA: &[(&str, f64)] = &[
    ("near-identity-7.5e-8", 7.5e-8),
    ("generic-1", 1.0),
    ("near-pi", core::f64::consts::PI - 1e-6),
];

/// `G::geodesic` behind a call LLVM may not inline, which is how every geodesic row times it.
///
/// Inlined into criterion's closure, a routine of a few hundred instructions is timed with
/// whatever LLVM made of it *there*: an unchanged `screw_pow` read 44.4 and 64.6 ns in two builds
/// differing only in other bench functions, and the provided body 79.5 and 59.5 (`0054`). A call
/// pins one compilation of the routine for baseline and candidate alike, so a ratio compares
/// routines and not two inlining decisions.
#[inline(never)]
fn geodesic_of<S: Real, G: LieGroup<S>>(x0: &G, x1: &G, t: S) -> G {
    G::geodesic(x0, x1, t)
}

/// How many points `act_many` moves per call: one patch of a frame's worth of landmarks.
const POINTS: usize = 64;

/// A bench fixture at the precision: the `f64` value, rounded. Not [`Real::lit`], whose argument is
/// stated to be exactly representable and whose `f32` impl checks it.
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

/// `φ` with `‖φ‖ = θ` along `(1, 2, 2)/3`, and a second direction for the two-element rows.
fn phi<S: Fixture>(theta: f64, other: bool) -> Vec3<S> {
    let t = theta / 3.0;
    let v = match other {
        false => [t, 2.0 * t, 2.0 * t],
        true => [-2.0 * t, t, 2.0 * t],
    };
    Vector(v.map(S::of))
}

/// `[ρ₁; …; ρ_N]`, fixed: a row differs from its neighbour in `θ` alone.
fn rho<S: Fixture, const N: usize>(other: bool) -> [Vec3<S>; N] {
    let base: [[f64; 3]; 3] = match other {
        false => [[0.5, -2.0, 0.25], [-1.0, 4.0, 0.125], [3.0, 0.75, -0.5]],
        true => [[-1.5, 0.25, 2.0], [0.5, -3.0, 1.25], [-0.25, 1.0, 4.0]],
    };
    core::array::from_fn(|i| Vector(base[i].map(S::of)))
}

/// The points `act_many` moves, spread over three decades so no row is a special case.
fn points<S: Fixture>() -> [Point3<S>; POINTS] {
    core::array::from_fn(|i| {
        let k = i as f64;
        Point([1.0 + k, 0.5 - 0.25 * k, -2.0 + 0.125 * k].map(S::of))
    })
}

/// SO(3)'s rows at one stratum.
fn so3<S: Fixture>(g: &mut BenchmarkGroup<'_, WallTime>, name: &str, theta: f64) {
    let (ta, tb) = (
        SO3Tangent {
            phi: phi::<S>(theta, false),
        },
        SO3Tangent {
            phi: phi::<S>(theta, true),
        },
    );
    let (x, y) = (SO3::<S>::exp(&ta), SO3::<S>::exp(&tb));
    // `x ⊖ base = Log(base⁻¹ x)`, so this `base` makes the `⊖` row's tangent `ta` and its stratum
    // the label's; `y` is a second element about another axis, which is what the other rows want.
    let base = x * SO3::<S>::exp(&SO3Tangent { phi: -ta.phi });
    let (ja, jb) = (SO3::<S>::jr(&ta), y.adjoint());
    let v = phi::<S>(1.0, true);
    let t = S::of(1.0 / 3.0);
    let pts = points::<S>();
    g.bench_function(format!("so3/exp/{name}"), |b| {
        b.iter(|| SO3::<S>::exp(black_box(&ta)));
    });
    g.bench_function(format!("so3/log/{name}"), |b| {
        b.iter(|| black_box(&x).log());
    });
    g.bench_function(format!("so3/compose/{name}"), |b| {
        b.iter(|| *black_box(&x) * *black_box(&y));
    });
    g.bench_function(format!("so3/act/{name}"), |b| {
        b.iter(|| black_box(&x).act(black_box(v)));
    });
    g.bench_function(format!("so3/act_many/{name}"), |b| {
        b.iter_batched_ref(
            || pts.map(|p| Vector(p.0)),
            |p| black_box(&x).act_many(p),
            BatchSize::SmallInput,
        );
    });
    g.bench_function(format!("so3/jr/{name}"), |b| {
        b.iter(|| SO3::<S>::jr(black_box(&ta)));
    });
    g.bench_function(format!("so3/jr_inv/{name}"), |b| {
        b.iter(|| SO3::<S>::jr_inv(black_box(&ta)));
    });
    g.bench_function(format!("so3/adjoint/{name}"), |b| {
        b.iter(|| black_box(&x).adjoint());
    });
    g.bench_function(format!("so3/jac_mul/{name}"), |b| {
        b.iter(|| Jac::<S, SO3Tangent<S>>::mul(black_box(&ja), black_box(&jb)));
    });
    g.bench_function(format!("so3/rminus_jacobians/{name}"), |b| {
        b.iter(|| black_box(&x).rminus_jacobians(black_box(&base)));
    });
    // `base` to `x`, so the relative motion is `ta` and `θ(d)` is the label's, not the angle
    // between two same-`θ` elements about different axes. `t` is non-dyadic for `laws`'s reason:
    // `0.25` would remove a rounding the shipped path has.
    g.bench_function(format!("so3/geodesic/{name}"), |b| {
        b.iter(|| geodesic_of::<S, SO3<S>>(black_box(&base), black_box(&x), black_box(t)));
    });
}

/// SE_N(3)'s rows at one stratum; `act` and `act_many` are `N = 1`'s alone (`PHASE3.md` §5).
fn sen3<S: Fixture, const N: usize>(
    g: &mut BenchmarkGroup<'_, WallTime>,
    tag: &str,
    name: &str,
    theta: f64,
) {
    let (ta, tb) = (
        SEn3Tangent {
            phi: phi::<S>(theta, false),
            rho: rho::<S, N>(false),
        },
        SEn3Tangent {
            phi: phi::<S>(theta, true),
            rho: rho::<S, N>(true),
        },
    );
    let (x, y) = (
        <helicoid::SEn3<S, N> as LieGroup<S>>::exp(&ta),
        <helicoid::SEn3<S, N> as LieGroup<S>>::exp(&tb),
    );
    // As in `so3`: the `⊖` row's tangent is `ta` itself, not `Log(y⁻¹ x)` of two same-`θ`
    // elements about different axes, which lands on another stratum.
    let base = x * <helicoid::SEn3<S, N> as LieGroup<S>>::exp(&ta.neg());
    let (ja, jb) = (<helicoid::SEn3<S, N> as LieGroup<S>>::jr(&ta), x.adjoint());
    let t = S::of(1.0 / 3.0);
    g.bench_function(format!("{tag}/exp/{name}"), |b| {
        b.iter(|| <helicoid::SEn3<S, N> as LieGroup<S>>::exp(black_box(&ta)));
    });
    g.bench_function(format!("{tag}/log/{name}"), |b| {
        b.iter(|| black_box(&x).log());
    });
    g.bench_function(format!("{tag}/compose/{name}"), |b| {
        b.iter(|| *black_box(&x) * *black_box(&y));
    });
    g.bench_function(format!("{tag}/jr/{name}"), |b| {
        b.iter(|| <helicoid::SEn3<S, N> as LieGroup<S>>::jr(black_box(&ta)));
    });
    g.bench_function(format!("{tag}/jr_inv/{name}"), |b| {
        b.iter(|| <helicoid::SEn3<S, N> as LieGroup<S>>::jr_inv(black_box(&ta)));
    });
    g.bench_function(format!("{tag}/adjoint/{name}"), |b| {
        b.iter(|| black_box(&x).adjoint());
    });
    g.bench_function(format!("{tag}/jac_mul/{name}"), |b| {
        b.iter(|| Jac::<S, SEn3Tangent<S, N>>::mul(black_box(&ja), black_box(&jb)));
    });
    g.bench_function(format!("{tag}/rminus_jacobians/{name}"), |b| {
        b.iter(|| black_box(&x).rminus_jacobians(black_box(&base)));
    });
    // `base` to `x`, as in `so3`: the relative motion is `ta`, so the row's `θ(d)` is its label's.
    // SE(3)'s row is `0054`'s screw twin; SE₂(3)'s is still the provided body.
    g.bench_function(format!("{tag}/geodesic/{name}"), |b| {
        b.iter(|| {
            geodesic_of::<S, helicoid::SEn3<S, N>>(black_box(&base), black_box(&x), black_box(t))
        });
    });
}

/// SE(3)'s action, which SE₂(3) does not have.
fn se3_action<S: Fixture>(g: &mut BenchmarkGroup<'_, WallTime>, name: &str, theta: f64) {
    let tau = SEn3Tangent {
        phi: phi::<S>(theta, false),
        rho: rho::<S, 1>(false),
    };
    let x = SE3::<S>::exp(&tau);
    let p = Point(phi::<S>(1.0, true).0);
    let pts = points::<S>();
    g.bench_function(format!("se3/act/{name}"), |b| {
        b.iter(|| *black_box(&x) * black_box(p));
    });
    g.bench_function(format!("se3/act_many/{name}"), |b| {
        b.iter_batched_ref(|| pts, |q| black_box(&x).act_many(q), BatchSize::SmallInput);
    });
}

/// SE(3)'s geodesic at `θ = 1e-3`, the relative motion of `tf_tree`'s `lookup/depth3/sclerp`
/// fixture (`0054`): both of the twin's coefficients on their short arms, which no §9 stratum
/// reaches -- `near-identity` is four decades lower and `generic` on the exact arms.
fn se3_consecutive<S: Fixture>(g: &mut BenchmarkGroup<'_, WallTime>) {
    let ta = SEn3Tangent {
        phi: phi::<S>(1e-3, false),
        rho: rho::<S, 1>(false),
    };
    let x = SE3::<S>::exp(&SEn3Tangent {
        phi: phi::<S>(1.0, true),
        rho: rho::<S, 1>(true),
    });
    let base = x * SE3::<S>::exp(&ta.neg());
    let t = S::of(1.0 / 3.0);
    g.bench_function("se3/geodesic/consecutive-1e-3", |b| {
        b.iter(|| geodesic_of::<S, SE3<S>>(black_box(&base), black_box(&x), black_box(t)));
    });
}

fn groups(c: &mut Criterion) {
    let mut g = c.benchmark_group("groups/f64");
    for &(name, theta) in THETA {
        so3::<f64>(&mut g, name, theta);
        sen3::<f64, 1>(&mut g, "se3", name, theta);
        se3_action::<f64>(&mut g, name, theta);
        sen3::<f64, 2>(&mut g, "se23", name, theta);
    }
    se3_consecutive::<f64>(&mut g);
    g.finish();

    let mut g = c.benchmark_group("groups/f32");
    for &(name, theta) in THETA {
        so3::<f32>(&mut g, name, theta);
        sen3::<f32, 1>(&mut g, "se3", name, theta);
        se3_action::<f32>(&mut g, name, theta);
        sen3::<f32, 2>(&mut g, "se23", name, theta);
    }
    se3_consecutive::<f32>(&mut g);
    g.finish();
}

// The time budget of `benches/coeffs.rs`, for its reasons: 100 samples, about 2 s per benchmark.
criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_millis(1500))
        .sample_size(100);
    targets = groups
}
criterion_main!(benches);
