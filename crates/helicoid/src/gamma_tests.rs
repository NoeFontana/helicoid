//! `so3::{gamma1, gamma2, gamma_apply_jacobian}` (`0064`, `docs/PHASE5.md` §4) against their
//! definitions: `Γ₁` is `J_l` to the bit, `Γ₂` at zero and its identities GG.2(b), (d), and the
//! directional Jacobian against `reference::gamma_apply` on `Dual<_, 3>`, the §14 twin.
//!
//! Draws are deterministic: `θ` log-uniform on `[1e-8, π]` in four bands, a uniform axis, `v`
//! uniform in `[-1, 1)³`. The twin's error is `laws::e`, the norm of the difference over
//! `max(‖reference‖, 1)`, in `u` of the subject's precision. An `f32` subject is scored against the
//! twin at `f64` from the same `f32` inputs; an `f64` subject against the `f64` twin, which is the
//! twin comparison D6 asks for and not an error against the exact value.

use std::vec;
use std::vec::Vec;

use helicoid_linalg::{Dual, Mat3, Precision, Real, Vec3, Vector};

use crate::laws::{e, worst};
use crate::reference::gamma_apply;
use crate::so3::{gamma1, gamma2, gamma_apply_jacobian};
use crate::{LieGroup, SO3Tangent, SO3};

/// The upper ends of the `θ` bands the bounds are stated on.
const BANDS: [f64; 4] = [1e-2, 1e-1, 1.0, core::f64::consts::PI];

/// `(φ, v)` draws, `θ` log-uniform over `[1e-8, π]`, plus `φ = 0` and a subnormal `φ`.
fn draws(count: usize) -> Vec<([f64; 3], [f64; 3])> {
    let mut st = 0x9E37_79B9_7F4A_7C15_u64;
    let mut next = || {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        (st >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut out = vec![
        ([0.0; 3], [1.0, -2.0, 0.5]),
        ([5e-324, -5e-324, 1e-310], [0.3, 0.7, -1.0]),
    ];
    for _ in 0..count {
        let theta = 10f64.powf(-8.0 + (8.0 + core::f64::consts::PI.log10()) * next());
        let (u, w) = (2.0 * next() - 1.0, 2.0 * core::f64::consts::PI * next());
        let r = (1.0 - u * u).sqrt();
        let axis = [r * w.cos(), r * w.sin(), u];
        let v = [2.0 * next() - 1.0, 2.0 * next() - 1.0, 2.0 * next() - 1.0];
        out.push((axis.map(|a| a * theta), v));
    }
    out
}

fn band(phi: &[f64; 3]) -> usize {
    let theta = phi.iter().map(|x| x * x).sum::<f64>().sqrt();
    BANDS
        .iter()
        .position(|&b| theta < b)
        .unwrap_or(BANDS.len() - 1)
}

/// `x` rounded to a value `S` holds: `Real::lit` takes only exact constants.
fn at<S: Real>(x: f64) -> S {
    match S::PRECISION {
        Precision::F64 => S::lit(x),
        Precision::F32 => S::lit(f64::from(x as f32)),
    }
}

fn tan<S: Real>(phi: [f64; 3]) -> SO3Tangent<S> {
    SO3Tangent {
        phi: Vector(phi.map(at)),
    }
}

fn f64s<S: Real>(v: Vec3<S>) -> [f64; 3] {
    v.0.map(|x| x.value_f64())
}

/// Column-major, as `laws::e` reads a matrix.
fn cols<S: Real>(m: &Mat3<S>) -> [f64; 9] {
    core::array::from_fn(|k| m.get(k % 3, k / 3).value_f64())
}

/// The twin at `f64` from the inputs as `S` holds them: value and `∂/∂φ`, column-major.
fn twin<S: Real>(m: usize, phi: &SO3Tangent<S>, v: Vec3<S>) -> ([f64; 3], [f64; 9]) {
    type D3 = Dual<f64, 3>;
    let p = SO3Tangent {
        phi: Vector(core::array::from_fn(|i| {
            D3::variable(phi.phi.0[i].value_f64(), i)
        })),
    };
    let g = gamma_apply(m, &p, Vector(v.0.map(|x| D3::constant(x.value_f64()))));
    (
        g.0.map(|x| x.v),
        core::array::from_fn(|k| g.0[k % 3].d[k / 3]),
    )
}

/// The worst `(value, Jacobian)` error per band of `gamma_apply_jacobian::<M>` against the twin.
fn twin_errors<const M: usize, S: Real>() -> [(f64, f64); 4] {
    let mut out = [(0.0, 0.0); 4];
    for (phi, v) in draws(20_000) {
        let (p, w) = (tan::<S>(phi), Vector(v.map(at::<S>)));
        let (gv, j) = gamma_apply_jacobian::<M, S>(&p, w);
        let (want_v, want_j) = twin(M, &p, w);
        let b = band(&phi);
        out[b].0 = worst(out[b].0, e::<S>(&f64s(gv), &want_v));
        out[b].1 = worst(out[b].1, e::<S>(&cols(&j), &want_j));
    }
    out
}

/// Twice the worst of the measurement below, per band, `(value, Jacobian)`, rounded up.
/// The `f64` rows are a comparison of two `f64` evaluations, so they carry the twin's own rounding:
/// its value at `M = 1` sums terms up to `θ²/2` against a result near 1 and is the larger side.
const BOUNDS_1_F64: [(f64, f64); 4] = [(5.0, 4.0), (7.0, 5.0), (11.0, 8.0), (16.0, 10.0)];
const BOUNDS_2_F64: [(f64, f64); 4] = [(3.0, 3.0), (5.0, 3.0), (7.0, 3.0), (9.0, 3.0)];
const BOUNDS_1_F32: [(f64, f64); 4] = [(3.0, 3.0), (3.0, 3.0), (4.0, 3.0), (8.0, 10.0)];
const BOUNDS_2_F32: [(f64, f64); 4] = [(2.0, 2.0), (2.0, 2.0), (2.0, 2.0), (4.0, 4.0)];

fn assert_within(got: [(f64, f64); 4], bounds: [(f64, f64); 4], what: &str) {
    for (b, (g, w)) in got.iter().zip(bounds).enumerate() {
        assert!(
            g.0 <= w.0 && g.1 <= w.1,
            "{what}, θ < {}: (value, Jacobian) = {g:?} u against {w:?}",
            BANDS[b]
        );
    }
}

#[test]
fn gamma_apply_jacobian_matches_reference_m1_f64() {
    assert_within(twin_errors::<1, f64>(), BOUNDS_1_F64, "M = 1, f64");
}

#[test]
fn gamma_apply_jacobian_matches_reference_m2_f64() {
    assert_within(twin_errors::<2, f64>(), BOUNDS_2_F64, "M = 2, f64");
}

#[test]
fn gamma_apply_jacobian_matches_reference_m1_f32() {
    assert_within(twin_errors::<1, f32>(), BOUNDS_1_F32, "M = 1, f32");
}

#[test]
fn gamma_apply_jacobian_matches_reference_m2_f32() {
    assert_within(twin_errors::<2, f32>(), BOUNDS_2_F32, "M = 2, f32");
}

/// `Γ₂`'s matrix against the twin's columns `Γ₂ e_k`, and `Γ_M v` from the matrices against the
/// vector form `gamma_apply_jacobian` evaluates: the worst per band, `(matrix, value)`.
fn matrix_errors<const M: usize, S: Real>() -> [(f64, f64); 4] {
    let mut out = [(0.0, 0.0); 4];
    for (phi, v) in draws(20_000) {
        let (p, w) = (tan::<S>(phi), Vector(v.map(at::<S>)));
        let g = if M == 1 { gamma1(&p) } else { gamma2(&p) };
        let pf = SO3Tangent {
            phi: Vector(p.phi.0.map(|x| x.value_f64())),
        };
        let want: [f64; 9] = core::array::from_fn(|k| {
            let mut unit = [0.0; 3];
            unit[k / 3] = 1.0;
            gamma_apply(M, &pf, Vector(unit)).0[k % 3]
        });
        let b = band(&phi);
        out[b].0 = worst(out[b].0, e::<S>(&cols(&g), &want));
        let (gv, _) = gamma_apply_jacobian::<M, S>(&p, w);
        out[b].1 = worst(out[b].1, e::<S>(&f64s(gv), &f64s(g * w)));
    }
    out
}

const MATRIX_1_F64: [(f64, f64); 4] = [(2.0, 5.0), (3.0, 5.0), (5.0, 7.0), (10.0, 9.0)];
const MATRIX_2_F64: [(f64, f64); 4] = [(2.0, 3.0), (3.0, 3.0), (4.0, 3.0), (6.0, 5.0)];
const MATRIX_1_F32: [(f64, f64); 4] = [(1.0, 5.0), (1.0, 6.0), (2.0, 6.0), (7.0, 8.0)];
const MATRIX_2_F32: [(f64, f64); 4] = [(1.0, 3.0), (1.0, 3.0), (2.0, 3.0), (4.0, 4.0)];

#[test]
fn gamma_matrices_match_reference_and_the_vector_form() {
    assert_within(matrix_errors::<1, f64>(), MATRIX_1_F64, "Γ₁, f64");
    assert_within(matrix_errors::<2, f64>(), MATRIX_2_F64, "Γ₂, f64");
    assert_within(matrix_errors::<1, f32>(), MATRIX_1_F32, "Γ₁, f32");
    assert_within(matrix_errors::<2, f32>(), MATRIX_2_F32, "Γ₂, f32");
}

/// `Γ₁` is `SO3::jl`, bit for bit (GG.2(a)): one evaluation for `Δv` and `Exp`'s translation.
#[test]
fn gamma1_is_jl_to_the_bit() {
    for (phi, _) in draws(2_000) {
        let p = tan::<f64>(phi);
        assert_eq!(
            cols(&gamma1(&p)).map(f64::to_bits),
            cols(&SO3::jl(&p)).map(f64::to_bits)
        );
    }
}

/// `Γ₂(0) = ½I` exactly (GG.2(d)), `½I + bW` with `W` subnormal at a subnormal `φ`, and the
/// directional Jacobian at both is `−v^/(M+1)!` (GG.5(b)) to the rounding of `1/(M+1)!`, finite
/// (GG.6(b)).
#[test]
fn gamma_at_zero_is_the_leading_term() {
    for phi in [[0.0; 3], [5e-324, -5e-324, 1e-310]] {
        let p = tan::<f64>(phi);
        let half: [f64; 9] = core::array::from_fn(|k| if k % 4 == 0 { 0.5 } else { 0.0 });
        for (g, h) in cols(&gamma2(&p)).iter().zip(half) {
            assert!((g - h).abs() <= 1e-300, "Γ₂ at {phi:?}: {g} against {h}");
        }
        let v = [0.3, -0.7, 1.1];
        let hat = |s: f64| {
            [
                0.0,
                v[2] * s,
                -v[1] * s,
                -v[2] * s,
                0.0,
                v[0] * s,
                v[1] * s,
                -v[0] * s,
                0.0,
            ]
        };
        let (_, j1) = gamma_apply_jacobian::<1, f64>(&p, Vector(v));
        let (_, j2) = gamma_apply_jacobian::<2, f64>(&p, Vector(v));
        // `−v^` column-major: column `c` is `−e_c × v`... written out as `v × e_c`.
        for (got, want) in [(cols(&j1), hat(-0.5)), (cols(&j2), hat(-1.0 / 6.0))] {
            for (g, w) in got.iter().zip(want) {
                assert!(
                    g.is_finite() && (g - w).abs() <= 1e-16,
                    "{got:?} against {want:?}"
                );
            }
        }
    }
}

/// GG.2(b) `Γ₁ = I + Γ₂W` within twice the measured 4.32 `u` (`laws::e`), and GG.2(d)
/// `Γ₂(−φ) = Γ₂(φ)ᵀ` exactly: `b`, `d` read `θ²` alone, `hat(−φ) = −hat(φ)`, and the structured
/// `W²` is symmetric entry for entry (`x y` against `y x`).
#[test]
fn gamma_identities_hold() {
    let (mut shift, mut transpose) = (0.0_f64, 0.0_f64);
    for (phi, _) in draws(5_000) {
        let p = tan::<f64>(phi);
        let w = helicoid_linalg::hat(p.phi);
        let lhs = gamma1(&p);
        let rhs = Mat3::identity() + gamma2(&p) * w;
        shift = worst(shift, e::<f64>(&cols(&lhs), &cols(&rhs)));
        let minus = tan::<f64>(phi.map(|x| -x));
        transpose = worst(
            transpose,
            e::<f64>(&cols(&gamma2(&minus)), &cols(&gamma2(&p).transpose())),
        );
    }
    assert!(
        shift <= 9.0 && transpose == 0.0,
        "GG.2(b) {shift} u, GG.2(d) {transpose} u"
    );
}

/// The measurement the bounds above are twice of: `cargo nextest run -p helicoid --run-ignored
/// only -- measure_gamma_twin --nocapture`.
#[test]
#[ignore = "measurement"]
#[allow(clippy::print_stdout)]
fn measure_gamma_twin() {
    std::println!("M=1 f64 {:?}", twin_errors::<1, f64>());
    std::println!("M=2 f64 {:?}", twin_errors::<2, f64>());
    std::println!("M=1 f32 {:?}", twin_errors::<1, f32>());
    std::println!("M=2 f32 {:?}", twin_errors::<2, f32>());
    std::println!("Γ1 f64 {:?}", matrix_errors::<1, f64>());
    std::println!("Γ2 f64 {:?}", matrix_errors::<2, f64>());
    std::println!("Γ1 f32 {:?}", matrix_errors::<1, f32>());
    std::println!("Γ2 f32 {:?}", matrix_errors::<2, f32>());
}
