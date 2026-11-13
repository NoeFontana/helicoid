//! The spline's latency (`0116` decision 6), the gate's no-regress rows; no absolute number is
//! set (`0113`).
//!
//! One sample is what a caller pays per evaluation: `Weights::at` plus the evaluation, on a
//! segment formed outside the timed loop. `segment_new` rows time the per-segment work alone, and
//! `batch64` times 64 samples on one segment, the shape of an IMU stream. Every routine sits
//! behind an `#[inline(never)]` call, for `helicoid`'s `benches/groups.rs` reason: a ratio then
//! compares routines, not two inlining decisions. binary64 only: the consumer's precision.
//!
//! Fixtures are deterministic: generic increments of about 0.6 rad about fixed axes, `u = 0.37`,
//! `h = 0.01 s` (a 100 Hz knot rate).

#![allow(
    missing_docs,
    reason = "criterion_main! generates an undocumented `main`"
)]

use core::hint::black_box;
use core::time::Duration;

use criterion::measurement::WallTime;
use criterion::{criterion_group, criterion_main, BenchmarkGroup, Criterion};
use helicoid::{LieGroup, Product, Right, Rn, Tangent, SE3, SO3};
use helicoid_spline::{JacSegment, Order, Segment, ValidOrder, Weights};

const U: f64 = 0.37;
const INV_H: f64 = 100.0;

fn window<G: LieGroup<f64>, const K: usize>() -> [G; K] {
    let mut x = G::identity();
    core::array::from_fn(|j| {
        let dense: [f64; 6] =
            core::array::from_fn(|i| 0.6 * ((j * 7 + i * 3) % 11) as f64 / 11.0 - 0.25);
        x = x * G::exp(&G::Tangent::read_dense(&dense[..G::Tangent::DOF]));
        x
    })
}

#[inline(never)]
fn eval<G: LieGroup<f64>, const K: usize, const D: usize>(
    seg: &Segment<f64, G, K>,
    u: f64,
) -> G::Tangent
where
    Order<K>: ValidOrder,
{
    let w = Weights::<f64, K, D>::at(u, INV_H);
    if D == 0 {
        seg.value(&w).log()
    } else {
        let k = seg.kinematics(&w);
        k.acc.add(&k.vel)
    }
}

#[inline(never)]
fn eval_jac<G: LieGroup<f64>, const K: usize, const D: usize>(
    seg: &JacSegment<f64, G, K>,
    u: f64,
) -> (G::Tangent, G::Jac)
where
    Order<K>: ValidOrder,
{
    let w = Weights::<f64, K, D>::at(u, INV_H);
    let (k, j) = seg.kinematics_jacobians::<Right, D>(&w);
    (k.vel, j.value[K / 2])
}

#[inline(never)]
fn segment_new<G: LieGroup<f64>, const K: usize>(ctrl: &[G; K]) -> Segment<f64, G, K>
where
    Order<K>: ValidOrder,
{
    Segment::new(ctrl)
}

#[inline(never)]
fn jac_segment_new<G: LieGroup<f64>, const K: usize>(ctrl: &[G; K]) -> JacSegment<f64, G, K>
where
    Order<K>: ValidOrder,
{
    JacSegment::new(ctrl)
}

fn rows<G: LieGroup<f64>, const K: usize>(
    g: &mut BenchmarkGroup<'_, WallTime>,
    name: &str,
    with_new: bool,
) where
    Order<K>: ValidOrder,
{
    let ctrl = window::<G, K>();
    let seg = Segment::new(&ctrl);
    let jseg = JacSegment::new(&ctrl);
    let p = format!("{name}/{K}");
    g.bench_function(format!("{p}/value/plain"), |b| {
        b.iter(|| eval::<G, K, 0>(black_box(&seg), black_box(U)))
    });
    g.bench_function(format!("{p}/kin1/plain"), |b| {
        b.iter(|| eval::<G, K, 1>(black_box(&seg), black_box(U)))
    });
    g.bench_function(format!("{p}/kin2/plain"), |b| {
        b.iter(|| eval::<G, K, 2>(black_box(&seg), black_box(U)))
    });
    g.bench_function(format!("{p}/value/jac"), |b| {
        b.iter(|| eval_jac::<G, K, 0>(black_box(&jseg), black_box(U)))
    });
    g.bench_function(format!("{p}/kin1/jac"), |b| {
        b.iter(|| eval_jac::<G, K, 1>(black_box(&jseg), black_box(U)))
    });
    g.bench_function(format!("{p}/kin2/jac"), |b| {
        b.iter(|| eval_jac::<G, K, 2>(black_box(&jseg), black_box(U)))
    });
    if with_new {
        g.bench_function(format!("{p}/segment_new"), |b| {
            b.iter(|| segment_new::<G, K>(black_box(&ctrl)))
        });
        g.bench_function(format!("{p}/jac_segment_new"), |b| {
            b.iter(|| jac_segment_new::<G, K>(black_box(&ctrl)))
        });
    }
}

type Split = Product<SO3<f64>, Rn<f64, 3>>;

fn spline(c: &mut Criterion) {
    let mut g = c.benchmark_group("spline/f64");
    rows::<SO3<f64>, 4>(&mut g, "so3", true);
    rows::<SO3<f64>, 5>(&mut g, "so3", true);
    rows::<SO3<f64>, 6>(&mut g, "so3", true);
    rows::<Rn<f64, 3>, 4>(&mut g, "rn3", false);
    rows::<Rn<f64, 3>, 5>(&mut g, "rn3", false);
    rows::<Rn<f64, 3>, 6>(&mut g, "rn3", false);
    rows::<Split, 4>(&mut g, "so3xr3", true);
    rows::<Split, 5>(&mut g, "so3xr3", true);
    rows::<Split, 6>(&mut g, "so3xr3", true);
    rows::<SE3<f64>, 4>(&mut g, "se3", false);
    rows::<SE3<f64>, 5>(&mut g, "se3", false);
    rows::<SE3<f64>, 6>(&mut g, "se3", false);

    // 64 samples on one segment: a 6.4 kHz IMU against a 100 Hz knot rate.
    let ctrl = window::<Split, 6>();
    let seg = Segment::new(&ctrl);
    let jseg = JacSegment::new(&ctrl);
    let us: [f64; 64] = core::array::from_fn(|i| (i as f64 + 0.5) / 64.0);
    g.bench_function("so3xr3/6/batch64/kin2", |b| {
        b.iter(|| {
            for &u in &us {
                black_box(eval::<Split, 6, 2>(black_box(&seg), u));
            }
        })
    });
    g.bench_function("so3xr3/6/batch64/kin2_jac", |b| {
        b.iter(|| {
            for &u in &us {
                black_box(eval_jac::<Split, 6, 2>(black_box(&jseg), u));
            }
        })
    });
    g.finish();
}

// `helicoid`'s time budget: 100 samples, about 2 s per benchmark.
criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_millis(1500))
        .sample_size(100);
    targets = spline
}
criterion_main!(benches);
