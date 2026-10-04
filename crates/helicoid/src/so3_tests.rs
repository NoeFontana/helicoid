//! `SO3` under the generic laws of [`crate::laws`], and its hand cases.
//!
//! The bounds are recorded as described at [`Bounds`], which asks for twice the worst error,
//! rounded up. The protocol: **60 000** tangents from an xorshift64 stream, each component
//! `m · 2^-e` with `m` uniform in `(-1, 1)` and `e = state mod 7` (the shape of `laws::sample`,
//! wider), every law on each sliding window of three, the maximum kept per law and per scalar.
//!
//! The worst errors, in `u` — `f64` and `Dual<f64, 3>` agreeing **to every digit**:
//!
//! | law | `f64`, `Dual` | `f32` |
//! |---|---|---|
//! | `ad_consistency` | 7.693 | 5.977 |
//! | `plus_minus` | 7.434 | 7.478 |
//! | `group_axioms` | 6.187 | 6.401 |
//! | `adjoint_identity` | 6.083 | 5.514 |
//! | `sandwich_matches_dense` | 3.849 | 2.783 |
//! | `jl_is_ad_jr` | 3.838 | 3.847 |
//! | `exp_log_roundtrip` | 3.551 | 4.072 |
//! | `jac_dense_order` | 2.951 | 3.464 |
//! | `side_delegation` | 1.118 | 1.118 |
//! | `tangent_dense_order` | **0** | 1.333 |
//! | `jacobian_rows` | **0** | **0** |
//!
//! The two zeros are the load-bearing ones: `jacobian_rows` at `0` says every row of
//! `NUMERICS.md` §2.3 is reproduced bit for bit by the closed forms, and the `f64`/`Dual`
//! agreement says the derivative lanes cost the value path nothing (D16).
//!
//! A first pass set these from 600 samples and got two of them wrong in opposite directions:
//! `f32` `group_axioms` was 6.0 against a measured 4.29 and `proptest` found 6.32, and `f64`
//! `ad_consistency` and `plus_minus` then sat at 8.0 against maxima of 7.69 and 7.43 — 4% and 7%
//! of headroom on laws `proptest` explores every run. A sample count is part of a bound.

use crate::laws::{laws_for, tangent, Bounds, Sample};
use crate::{LieGroup, Quat, Tangent, SO3};
use helicoid_linalg::{hat, Mat3, Matrix, Vector};

type G<S> = SO3<S>;

// Twice the larger of the two measured maxima, rounded up, as `Bounds` asks. One set serves both
// precisions because the two agree within 25% on every law; only `tangent_dense_order` differs in
// kind, being exactly `0` at `f64` and one rounding at `f32`. The two zeros keep no headroom at
// all, because they are exactness claims and not a margin.
const F64: Bounds = Bounds {
    axioms: 13.0,
    exp_log: 9.0,
    adjoint: 13.0,
    jl_ad_jr: 8.0,
    plus_minus: 15.0,
    rows: 0.0,
    ad: 16.0,
    sides: 3.0,
    tangent_order: 0.0,
    jac_order: 7.0,
    sandwich: 8.0,
};
const F32: Bounds = Bounds {
    tangent_order: 3.0,
    ..F64
};

/// `J_r(φ)` at the sample: invertible for every `|φ| < 2π`, which `laws::sample` is
/// (`|φ| < √3`), and not symmetric, so `apply_transpose` is a real second case.
fn so3_jac<S: Sample>(v: &[f64; 3]) -> Mat3<S> {
    <SO3<S> as LieGroup<S>>::jr(&tangent::<S, SO3<S>, 3>(v))
}

laws_for!(so3, G, so3_jac, 3, F64, F32);

fn tan(phi: [f64; 3]) -> <SO3<f64> as LieGroup<f64>>::Tangent {
    <SO3<f64> as LieGroup<f64>>::Tangent::read_dense(&phi)
}

/// Bit equality: the crate states an exactness claim this way, never with `==` on floats.
fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

fn close(a: &[f64], b: &[f64], tol: f64) {
    let err = a
        .iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max);
    assert!(err <= tol, "{a:?} vs {b:?} (tol {tol})");
}

/// A unit quaternion of angle `theta` about `(1, 2, 2)/3`, from the standard library.
fn unit(theta: f64) -> Quat<f64> {
    let (s, c) = (theta / 2.0).sin_cos();
    Quat::from_wxyz_unchecked(c, s / 3.0, 2.0 * s / 3.0, 2.0 * s / 3.0)
}

#[test]
fn exp_is_the_half_angle_quaternion() {
    let q = SO3::<f64>::exp(&tan([0.0; 3])).quat();
    // Exact at the origin: the series arms are constants there.
    assert_eq!(
        [q.w, q.x, q.y, q.z].map(f64::to_bits),
        [1.0, 0.0, 0.0, 0.0].map(f64::to_bits)
    );
    for theta in [1e-9, 1e-3, 0.3, 1.2, 3.0, core::f64::consts::PI] {
        let phi = [theta / 3.0, 2.0 * theta / 3.0, 2.0 * theta / 3.0];
        let got = SO3::<f64>::exp(&tan(phi)).quat();
        let want = unit(theta);
        close(
            &[got.w, got.x, got.y, got.z],
            &[want.w, want.x, want.y, want.z],
            4e-16,
        );
    }
}

#[test]
fn log_is_a_function_of_the_quaternion_at_w_plus_zero_and_of_the_rotation_elsewhere() {
    // §3.2: `q` and `-q` are one rotation and `Log` flips, except at `w = +0`.
    for theta in [1e-9, 0.3, 2.0, 3.0] {
        let q = unit(theta);
        let neg = Quat::from_wxyz_unchecked(-q.w, -q.x, -q.y, -q.z);
        let bits = |r: SO3<f64>| {
            let mut out = [0.0; 3];
            r.log().write_dense(&mut out);
            out.map(f64::to_bits)
        };
        assert_eq!(
            bits(SO3::from_quat_unchecked(q)),
            bits(SO3::from_quat_unchecked(neg)),
            "theta {theta}"
        );
    }
    // At exactly `theta = pi`, `w = +0` does not flip and `w = -0` does: `±pi n̂`.
    let (u, pi) = ([0.6, 0.0, 0.8], core::f64::consts::PI);
    let log_of = |w: f64, s: f64| {
        let mut out = [0.0; 3];
        SO3::from_quat_unchecked(Quat::from_wxyz_unchecked(w, s * u[0], s * u[1], s * u[2]))
            .log()
            .write_dense(&mut out);
        out
    };
    close(&log_of(0.0, 1.0), &u.map(|c| pi * c), 4e-16);
    close(&log_of(0.0, -1.0), &u.map(|c| -pi * c), 4e-16);
    close(&log_of(-0.0, -1.0), &u.map(|c| pi * c), 4e-16);
}

#[test]
fn log_is_scale_invariant_across_the_unchecked_domain() {
    // §3.2: `atan2` and `u/n` are exactly scale-invariant, so a `q` that is unit only within
    // rounding returns the `Log` of its normalization. The scales tested are the ones
    // `from_wxyz_unchecked` admits — `|‖q‖² − 1| <= 2^-40` (§3.6) — because a quaternion outside
    // that cannot be built at all in a debug build, which is the domain working as specified.
    let (s, c) = 0.5_f64.sin_cos();
    for lambda in [1.0, 1.0 + 2f64.powi(-42), 1.0 - 2f64.powi(-42)] {
        let mut out = [0.0; 3];
        SO3::from_quat_unchecked(Quat::from_wxyz_unchecked(c * lambda, s * lambda, 0.0, 0.0))
            .log()
            .write_dense(&mut out);
        close(&out, &[1.0, 0.0, 0.0], 1e-15);
    }
    // Through the normalizing constructor, every scale is legal and must give the same tangent.
    // The quaternion is built from its fields, not by `Quat::from_wxyz_normalized`, which would
    // divide `lambda` out before `SO3::from_quat_normalized` -- the method under test -- ever saw
    // a non-unit input, and would handle the subnormal `‖q‖²` case itself.
    for lambda in [0.5, 2.0, 1e3, 1.4e-154] {
        let mut out = [0.0; 3];
        let q = Quat {
            w: c * lambda,
            x: s * lambda,
            y: 0.0,
            z: 0.0,
        };
        SO3::from_quat_normalized(q).log().write_dense(&mut out);
        close(&out, &[1.0, 0.0, 0.0], 1e-15);
    }
}

#[test]
fn from_matrix_round_trips_and_reaches_every_pivot() {
    // The four pivots: the trace one near the identity, and one per axis near a half turn about
    // it, where `1 + tr R` is the smallest of the four candidates.
    let cases = [
        [0.0, 0.0, 0.0],
        [0.2, -0.1, 0.05],
        [core::f64::consts::PI - 1e-6, 0.0, 0.0],
        [0.0, core::f64::consts::PI - 1e-6, 0.0],
        [0.0, 0.0, core::f64::consts::PI - 1e-6],
        [1.2, -2.0, 0.7],
    ];
    for phi in cases {
        let r = SO3::<f64>::exp(&tan(phi));
        let back = SO3::<f64>::from_matrix(&r.to_matrix());
        // Up to the double cover: compare the rotation, not the quaternion.
        let (a, b) = (r.to_matrix(), back.to_matrix());
        for c in 0..3 {
            close(&a.col(c).0, &b.col(c).0, 8e-16);
        }
    }
}

#[test]
fn from_matrix_never_iterates_on_a_degenerate_matrix() {
    // `PHASE3.md` §4: the closed form completes on input where nalgebra's Müller iteration does
    // not. Every one of these is handed to `from_matrix` and must return a finite quaternion.
    let degenerate = [
        Matrix::<f64, 3, 3>::identity().scale(0.0),
        Matrix::identity().scale(-1.0),
        Matrix::from_cols([
            Vector([1.0, 0.0, 0.0]),
            Vector([0.0, -1.0, 0.0]),
            Vector([0.0, 0.0, -1.0]),
        ]),
        Matrix::from_cols([Vector([1e-300, 0.0, 0.0]); 3]),
        Matrix::from_cols([Vector([1.0, 1.0, 1.0]); 3]),
    ];
    for m in degenerate {
        let q = SO3::<f64>::from_matrix(&m).quat();
        for c in [q.w, q.x, q.y, q.z] {
            assert!(c.is_finite(), "non-finite {c} from {m:?}");
        }
        // Either a unit quaternion, or the all-zero one a zero matrix can only give.
        let n = q.norm_sq();
        assert!((n - 1.0).abs() < 1e-12 || n < 1e-300, "norm_sq {n}");
    }
}

#[test]
fn act_agrees_with_the_matrix_and_act_many() {
    let r = SO3::<f64>::exp(&tan([0.3, -1.1, 0.7]));
    let m = r.to_matrix();
    let pts = [
        Vector([1.0, 0.0, 0.0]),
        Vector([0.0, 2.0, -1.0]),
        Vector([-3.0, 0.5, 1.5]),
    ];
    let mut many = pts;
    r.act_many(&mut many);
    for (i, p) in pts.into_iter().enumerate() {
        close(&r.act(p).0, &(m * p).0, 8e-16);
        close(&many[i].0, &(m * p).0, 0.0);
    }
}

#[test]
fn the_jacobians_are_the_closed_forms_of_3_5() {
    for phi in [
        [0.4, -0.2, 0.9],
        [1e-9, 0.0, 0.0],
        [0.0; 3],
        [2.9, 0.4, -0.2],
    ] {
        let t = tan(phi);
        let (jr, jl) = (SO3::<f64>::jr(&t), SO3::<f64>::jl(&t));
        let w = hat(Vector(phi));
        // `J_r = I − aW + bW²`, `J_l = I + aW + bW²` (§3.5): the sum is symmetric and the
        // difference is `−2aW`, antisymmetric and a multiple of `W`. Both to the bit, since `a`
        // and `b` are functions of `θ²` and `(−W)² = W²` in the same products.
        for c in 0..3 {
            for r in 0..3 {
                let (s, d) = (
                    jr.row(r).0[c] + jl.row(r).0[c],
                    jr.row(r).0[c] - jl.row(r).0[c],
                );
                assert!(
                    same(s, jr.row(c).0[r] + jl.row(c).0[r]),
                    "sum not symmetric at ({r}, {c})"
                );
                // Antisymmetry as a sum, not a negation: `d` and its transpose are equal and
                // opposite, so they cancel to `+0` exactly, whereas `same(d, -d_t)` would compare
                // `+0` against `-0` on the diagonal and fail for the sign alone.
                assert!(
                    same((d + (jr.row(c).0[r] - jl.row(c).0[r])).abs(), 0.0),
                    "difference not antisymmetric at ({r}, {c})"
                );
                // Where `W` is zero the difference `−2aW` must be too, to the bit.
                if same(w.row(r).0[c].abs(), 0.0) {
                    assert!(same(d.abs(), 0.0), "difference is not a multiple of W");
                }
            }
        }
        // `J_r J_r⁻¹ = I` and `J_l J_l⁻¹ = I` (`NUMERICS.md` §3.5, domain `θ < 2π`).
        for (j, inv) in [(jr, SO3::<f64>::jr_inv(&t)), (jl, SO3::<f64>::jl_inv(&t))] {
            let p = j * inv;
            for c in 0..3 {
                let want: [f64; 3] = core::array::from_fn(|r| f64::from(u8::from(r == c)));
                close(&p.col(c).0, &want, 8e-16);
            }
        }
    }
}

/// `J_l = J_rᵗ` and `J_l⁻¹ = J_r⁻¹ᵗ` to the bit, which is what lets `rminus_jacobians` and
/// `lminus_jacobians` answer both rows from one `jr_inv` and a transpose.
#[test]
fn jl_is_jr_transposed_to_the_bit() {
    let mut st = 0x1234_5678_9abc_def0_u64;
    let mut next = || {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        (st >> 11) as f64 / (1u64 << 53) as f64 * 4.0 - 2.0
    };
    for _ in 0..4000 {
        let t = tan([next(), next(), next()]);
        let (jr, jl) = (SO3::<f64>::jr(&t), SO3::<f64>::jl(&t));
        let (jri, jli) = (SO3::<f64>::jr_inv(&t), SO3::<f64>::jl_inv(&t));
        for r in 0..3 {
            for c in 0..3 {
                assert!(same(jl.row(r).0[c], jr.row(c).0[r]), "J_l != J_rᵗ");
                assert!(same(jli.row(r).0[c], jri.row(c).0[r]), "J_l⁻¹ != J_r⁻¹ᵗ");
            }
        }
    }
    // And the two `⊖` rows are exactly the independent forms, to the bit.
    let (x, y) = (
        SO3::<f64>::exp(&tan([0.3, -1.1, 0.7])),
        SO3::<f64>::exp(&tan([-0.4, 0.2, 0.9])),
    );
    let (rt, rb) = x.rminus_jacobians(&y);
    let (lt, lb) = x.lminus_jacobians(&y);
    let (tr, tl) = (x.rminus(&y), x.lminus(&y));
    for (got, want) in [
        (rt, SO3::<f64>::jr_inv(&tr)),
        (rb, -SO3::<f64>::jl_inv(&tr)),
        (lt, SO3::<f64>::jl_inv(&tl)),
        (lb, -SO3::<f64>::jr_inv(&tl)),
    ] {
        for c in 0..3 {
            close(&got.col(c).0, &want.col(c).0, 0.0);
        }
    }
}

/// `from_matrix` is closed-form and returns for every input, in a **debug** build too: a matrix
/// whose pivot candidate overflows, or an all-NaN one, is out of domain and gives a non-finite
/// quaternion rather than a panic (`Quat::from_wxyz_normalized`'s `debug_assert!` would).
#[test]
fn from_matrix_does_not_panic_out_of_domain() {
    let huge = Matrix::from_cols([
        Vector([1e308, 0.0, 0.0]),
        Vector([0.0, 1e308, 0.0]),
        Vector([0.0, 0.0, -1e308]),
    ]);
    let nan = Matrix::from_cols([Vector([f64::NAN; 3]); 3]);
    for m in [
        huge,
        nan,
        Matrix::from_cols([Vector([f64::INFINITY; 3]); 3]),
    ] {
        let q = SO3::<f64>::from_matrix(&m).quat();
        // Out of domain, so the value means nothing; what is asserted is that it returns at all.
        let _ = [q.w, q.x, q.y, q.z];
    }
}

/// D6: `act_many`'s `NUMERICS.md` §14 twin is the per-point `act`. Nothing in `laws.rs` touches
/// either, so without this the row had a twin and no `*_matches_reference` proptest, and the
/// twin-table lint that would have caught that is itself owed (`PHASE3.md` §0.0).
///
/// **§14 states a tolerance of `3 u` for this row and the two diverge by more.** Measured over
/// 200 000 samples per regime: `7.587 u` for `θ` near `π`, `4.873 u` for `θ` of order 1, and
/// `3.379 u` for `θ ~ 1e-6`, with points of order 1. It is not an implementation choice to fix —
/// §3.3 *prescribes* that `act_many` form `R(q)` once, so its nine entries round before any point
/// is touched while `act` rounds a sandwich per point, and the gap is the two algorithms. The
/// bound below is the measured figure with headroom; §14's `3` was never exercised, because
/// nothing implemented `act_many` until now, and correcting it is a normative edit and a record
/// (`CLAUDE.md`), not this PR's to make.
mod act_twin {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn act_many_matches_reference(
            a in prop::array::uniform3(-2.0_f64..2.0),
            p in prop::array::uniform9(-4.0_f64..4.0),
        ) {
            let r = SO3::<f64>::exp(&tan(a));
            let pts = [
                Vector([p[0], p[1], p[2]]),
                Vector([p[3], p[4], p[5]]),
                Vector([p[6], p[7], p[8]]),
            ];
            let mut many = pts;
            r.act_many(&mut many);
            for (i, q) in pts.into_iter().enumerate() {
                let err = crate::laws::e::<f64>(&many[i].0, &r.act(q).0);
                prop_assert!(err <= 10.0, "act_many vs act: {err} u at point {i}");
            }
        }
    }
}

#[test]
fn adjoint_is_the_rotation_and_ad_is_the_hat() {
    let r = SO3::<f64>::exp(&tan([0.3, -1.1, 0.7]));
    let (ad, m) = (r.adjoint(), r.to_matrix());
    for c in 0..3 {
        close(&ad.col(c).0, &m.col(c).0, 0.0);
    }
    let phi = [0.4, -0.2, 0.9];
    let w = helicoid_linalg::hat(Vector(phi));
    let got = SO3::<f64>::ad(&tan(phi));
    for c in 0..3 {
        close(&got.col(c).0, &w.col(c).0, 0.0);
    }
}
