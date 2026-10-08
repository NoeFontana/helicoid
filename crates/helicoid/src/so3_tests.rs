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
use crate::{Left, LieGroup, Quat, Right, SO3Tangent, Tangent, SO3};
use helicoid_linalg::{hat, Dual, Mat3, Matrix, Vector};

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
    // `PHASE3.md` §8's second check, measured 5.485 over 20 000 draws.
    dual_rows: 11.0,
    // `PHASE4.md` §1 and §3's seven legs, in `GEODESIC_LEGS`'s order, each twice the worst of 10^6 draws of `laws::Rng::shaped` -- `laws::sample`'s own distribution, which is what the proptest draws -- rounded up, `f64` / `f32`: symmetry 5.000 / 6.013, velocity 7.024 / 7.107, left 9.276 / 8.140, right 10.979 / 9.289, twin 7.213 / 7.233. **`t=0` and `t=1` are both 1.118**, `gerr`'s floor for two bitwise equal quaternions (`Bounds::geodesic` says why). `t=1` reading the floor is `0050`'s property and `0051` keeps it on both arms, which is why that record routes exactly `t = 1` to the blend even below the switch: the provided body is exact at `t=0` alone and read 7.160 here. `twin` is 1.118 no longer -- `SO3::geodesic` is two arms of which one is not `reference::geodesic`, so the leg compares genuinely different expressions where it once compared a call with itself. **Every figure is `0050`'s to the digit**, so `0051`'s second arm costs these laws nothing: its draws that fall below the switch are too few to move a maximum. An earlier `0051` read symmetry 5.919 / 7.670 and that was the extrapolation defect its own review found, not the dispatch -- `t >= 1` sent large-`t` draws below the switch to the blend, whose weights grow like `t`, and `t = 1` does not.
    geodesic: [3.0, 3.0, 13.0, 15.0, 19.0, 22.0, 15.0],
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
    // `from_wxyz_unchecked` admits, `|‖q‖² − 1| <= 2^-40` (§3.6). A carried quaternion beyond
    // that bound, inside the drift band, is `carried::each_operation_errs_by_its_first_order_term_across_the_drift_band`'s
    // case (`0058`).
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
/// Bit equality with the crate's standing NaN exception: which NaN an arithmetic operation
/// returns is the implementation's choice of input payload, so two NaNs match.
///
/// `laws_for!`'s `dual_value_is_plain_value` arm reads the same way. Without it, the structured
/// products' degenerate counts differ between the dev and release profiles -- 216 against 234 --
/// because the two instruction schedules pick a different one of two NaN operands.
fn same_or_nan(a: f64, b: f64) -> bool {
    same(a, b) || (a.is_nan() && b.is_nan())
}

/// §2.4's action rows against `Dual<f64, 3>` differentiation of the action itself, both sides and
/// `∂/∂p` — `PHASE3.md` §8's second check for the action, as `jacobians_match_dual` is for §2.3.
///
/// `11` is twice the worst of 20 000 draws, 5.066 `u` over both sides and both arguments. It is
/// not a differentiation error alone: `act` is §3.3's quaternion sandwich where the closed form
/// goes through `to_matrix`, so the comparison carries that difference too — the same one §14's
/// `act_many` row records.
#[test]
fn act_jacobians_differentiate_the_action() {
    type D3 = Dual<f64, 3>;
    let mut st = 0x6A09_E667_F3BC_C908_u64;
    let mut next = || {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        (st >> 11) as f64 / (1u64 << 53) as f64 * 4.0 - 2.0
    };
    let mut worst = 0.0_f64;
    for _ in 0..20_000 {
        let (phi, p) = ([next(), next(), next()], [next(), next(), next()]);
        let x = SO3::<f64>::exp(&tan(phi));
        let xd = SO3::<D3>::exp(&SO3Tangent {
            phi: Vector(phi.map(D3::constant)),
        });
        let pd = Vector(p.map(D3::constant));
        // `δ`, seeded: lane `i` is `∂/∂δ_i`, and the value is `0`.
        let delta = SO3Tangent {
            phi: Vector(core::array::from_fn(|i| D3::variable(0.0, i))),
        };
        for right in [true, false] {
            let moved = match right {
                true => xd.rplus(&delta),
                false => xd.lplus(&delta),
            };
            let got = moved.act(pd);
            let (j, dp) = match right {
                true => x.act_jacobians::<Right>(Vector(p)),
                false => x.act_jacobians::<Left>(Vector(p)),
            };
            // Column-major, as `Jac::write_dense` and `laws::e` read a matrix.
            let lanes: [[f64; 3]; 3] =
                core::array::from_fn(|c| core::array::from_fn(|r| got.0[r].d[c]));
            let want: [[f64; 3]; 3] =
                core::array::from_fn(|c| core::array::from_fn(|r| j.get(r, c)));
            worst = crate::laws::worst(
                worst,
                crate::laws::e::<f64>(lanes.as_flattened(), want.as_flattened()),
            );
            // `∂(R p)/∂p = R`, by differentiating in `p` instead.
            let pv = Vector(core::array::from_fn(|i| D3::variable(p[i], i)));
            let moved = xd.act(pv);
            let lanes: [[f64; 3]; 3] =
                core::array::from_fn(|c| core::array::from_fn(|r| moved.0[r].d[c]));
            let want: [[f64; 3]; 3] =
                core::array::from_fn(|c| core::array::from_fn(|r| dp.get(r, c)));
            worst = crate::laws::worst(
                worst,
                crate::laws::e::<f64>(lanes.as_flattened(), want.as_flattened()),
            );
        }
    }
    assert!(worst <= 11.0, "{worst} u");
}

/// `mul_hat` against the generic `Mat3 * hat(v)` it replaces, entry by entry on the bits.
///
/// Finite entries must agree exactly — that is the whole licence for skipping `hat`'s zeros — and
/// the two stated exceptions are counted, not assumed: a partial sum of `∓0`, which the dropped
/// `0 · x` term would have normalized to `+0`, and a non-finite entry of `a`, which `0 · x` would
/// have spread across the column as NaN.
#[test]
fn the_structured_product_is_the_generic_one_to_the_bit() {
    use crate::so3::mul_hat;
    let mut st = 0x9E37_79B9_7F4A_7C15_u64;
    let mut next = || {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        (st >> 11) as f64 / (1u64 << 53) as f64 * 4.0 - 2.0
    };
    let agree = |x: &Mat3<f64>, y: &Mat3<f64>| {
        (0..9).all(|i| same_or_nan(x.get(i / 3, i % 3), y.get(i / 3, i % 3)))
    };
    // Finite, generic: exact agreement, and no headroom.
    for _ in 0..20_000 {
        let a = Matrix::from_rows(core::array::from_fn(|_| Vector([next(), next(), next()])));
        let v = Vector([next(), next(), next()]);
        assert!(agree(&mul_hat(&a, v), &(a * hat(v))), "finite");
    }
    // Degenerate and non-finite, where the dropped term did something: counted.
    let odd = [
        0.0,
        -0.0,
        1.0,
        -1.0,
        f64::INFINITY,
        -f64::INFINITY,
        f64::NAN,
    ];
    let (mut zeros, mut nonfinite, mut total) = (0usize, 0usize, 0usize);
    for &e in &odd {
        for &f in &odd {
            let a = Matrix::from_rows([
                Vector([e, f, 0.0]),
                Vector([-0.0, e, f]),
                Vector([f, 0.0, e]),
            ]);
            for &g in &odd {
                let v = Vector([g, -0.0, 0.0]);
                total += 1;
                let (s, m) = (mul_hat(&a, v), a * hat(v));
                if !agree(&s, &m) {
                    match [e, f, g].iter().any(|x| !x.is_finite()) {
                        true => nonfinite += 1,
                        false => {
                            zeros += 1;
                            // The invariant, not just the count: with every entry finite, what
                            // differs is a zero's sign and nothing else.
                            for i in 0..9 {
                                let (a, b) = (s.get(i / 3, i % 3), m.get(i / 3, i % 3));
                                assert!(
                                    same(a, b) || (a == 0.0 && b == 0.0),
                                    "{e} {f} {g}: {a} vs {b}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }
    // Pinned, not bounded: 28 of the finite degenerate cases keep a `∓0` the dropped term would
    // have normalized, and 200 carry a non-finite entry that `0 · x` would have spread as NaN.
    // Every one of the 20 000 generic finite cases above agrees, which is the licence. The counts
    // are the same in both profiles *because* of `same_or_nan`: on `to_bits` alone they read 216
    // and 234, since the two instruction schedules pick a different one of two NaN operands.
    assert_eq!(
        (zeros, nonfinite, total),
        (28, 200, 343),
        "signed-zero / non-finite / cases"
    );
}

/// `hat_mul` against the generic `hat(v) * Mat3`, the mirror of
/// `the_structured_product_is_the_generic_one_to_the_bit` and held to the same licence: finite
/// entries exactly, and the degenerate cases counted.
#[test]
fn the_structured_left_product_is_the_generic_one_to_the_bit() {
    use crate::so3::hat_mul;
    let mut st = 0x2545_F491_4F6C_DD1D_u64;
    let mut next = || {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        (st >> 11) as f64 / (1u64 << 53) as f64 * 4.0 - 2.0
    };
    let agree = |x: &Mat3<f64>, y: &Mat3<f64>| {
        (0..9).all(|i| same_or_nan(x.get(i / 3, i % 3), y.get(i / 3, i % 3)))
    };
    for _ in 0..20_000 {
        let b = Matrix::from_rows(core::array::from_fn(|_| Vector([next(), next(), next()])));
        let v = Vector([next(), next(), next()]);
        assert!(agree(&hat_mul(v, &b), &(hat(v) * b)), "finite");
    }
    let odd = [
        0.0,
        -0.0,
        1.0,
        -1.0,
        f64::INFINITY,
        -f64::INFINITY,
        f64::NAN,
    ];
    let (mut zeros, mut nonfinite, mut total) = (0usize, 0usize, 0usize);
    for &e in &odd {
        for &f in &odd {
            let b = Matrix::from_rows([
                Vector([e, f, 0.0]),
                Vector([-0.0, e, f]),
                Vector([f, 0.0, e]),
            ]);
            for &g in &odd {
                let v = Vector([g, -0.0, 0.0]);
                total += 1;
                let (h, m) = (hat_mul(v, &b), hat(v) * b);
                if !agree(&h, &m) {
                    match [e, f, g].iter().any(|x| !x.is_finite()) {
                        true => nonfinite += 1,
                        false => {
                            zeros += 1;
                            for i in 0..9 {
                                let (a, c) = (h.get(i / 3, i % 3), m.get(i / 3, i % 3));
                                assert!(
                                    same(a, c) || (a == 0.0 && c == 0.0),
                                    "{e} {f} {g}: {a} vs {c}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }
    // Pinned as its mirror is, and the same in both profiles for its reason.
    assert_eq!(
        (zeros, nonfinite, total),
        (44, 212, 343),
        "signed-zero / non-finite / cases"
    );
}

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

/// `γ(x₀, x₁, 1)` is `x₁` **bit for bit**, up to `Log`'s sign — the property `0050`'s override
/// adds and the whole reason it dominates oracle #1.
///
/// It holds by the arithmetic, not by luck: the right weight is `sin(1·α)/sin α`, and `1·α` is `α`
/// exactly, so it is one number divided by **itself** and therefore exactly `1`; the left weight is
/// `sin(0·α)/sin α`, an exact `+0`. Dividing by `‖v‖` instead -- the same number for a unit
/// quaternion, and the spelling `0050` measured and rejected -- gives `1 + ε` here, which is the
/// entire 2.721 -> 1.738 `u` at `geo:generic`.
///
/// The sign is `Log`'s: at `q₀·q₁ < 0` the answer is `−q₁`, the same rotation on the short arc,
/// which is what `tf_tree_math::slerp` documents for its own `s = 1`.
///
/// **Up to the sign of a zero**, which is the `t = 0` test's exception and reaches `t = 1` for the
/// mirror reason: the left weight is an exact `+0`, so `a·q₀[i] + b·q₁[i]` adds `+0.0` to a `−0.0`
/// component of `q₁` and `−0.0 + 0.0` is `+0.0`. The flip makes it reachable — `Exp`'s unused
/// vector components are `+0.0`, so negating `q₁` makes them `−0.0` — which is why the case below
/// at `‖φ‖ = 3.5` is the one that found this, and why the five cases that were here before (all
/// with `q₀·q₁ > 0`) did not.
#[test]
fn geodesic_at_one_is_the_right_endpoint_bit_for_bit() {
    let q =
        |v: [f64; 4]| SO3::from_quat_unchecked(Quat::from_wxyz_normalized(v[0], v[1], v[2], v[3]));
    let cases = [
        ([1.0, 0.0, 0.0, 0.0], [0.3, -0.5, 0.7]),
        ([0.2, 0.3, -0.5, 0.78], [1.0, -2.0, 0.5]),
        ([0.0, 0.6, 0.0, 0.8], [1e-9, 0.0, -1e-9]),
        // `θ(d) = 0`: the arc is a point, which takes the `‖v‖ = 0` branch and not the blend.
        ([-0.4, 0.1, 0.9, -0.2], [0.0, 0.0, 0.0]),
        // Past `π`, so `Log`'s flip fires and the expected endpoint is `−q₁`. `3.5` and not
        // `2.6`: `‖φ‖ = 2.6` gives `w = cos(1.3) = +0.267`, which does **not** flip, so the case
        // that was there exercised none of the sign rule and neither did the other four.
        ([1.0, 0.0, 0.0, 0.0], [3.5, 0.0, 0.0]),
    ];
    let mut flipped = 0;
    for (w, d) in cases {
        let x0 = q(w);
        let x1 = x0.rplus(&SO3Tangent { phi: Vector(d) });
        let got = SO3::geodesic(&x0, &x1, 1.0).quat();
        let (a, b) = (x0.quat(), x1.quat());
        // The sign `Log` would choose, read off the relative quaternion as `geodesic` reads it.
        let dot = ((a.w * b.w + a.x * b.x) + a.y * b.y) + a.z * b.z;
        let f = 1.0_f64.copysign(dot);
        flipped += usize::from(f < 0.0);
        let want = [f * b.w, f * b.x, f * b.y, f * b.z];
        // Bits, with `±0.0` equal: the rustdoc's one exception. A `float_cmp` on a zero is what
        // `partial_cmp` is for, as the `t = 0` test's own `w` check uses.
        for (i, (g, v)) in [got.w, got.x, got.y, got.z].iter().zip(want).enumerate() {
            let zeros = g.partial_cmp(&0.0) == Some(core::cmp::Ordering::Equal)
                && v.partial_cmp(&0.0) == Some(core::cmp::Ordering::Equal);
            assert!(
                zeros || g.to_bits() == v.to_bits(),
                "t = 1 moved {w:?} along {d:?}: component {i} is {g:?}, not {v:?}"
            );
        }
    }
    assert!(
        flipped > 0,
        "no case has `q₀·q₁ < 0`, so the `−q₁` half of this property is untested and the claim \
         `NUMERICS.md` §10 and `PHASE4.md` §0.0 both make is uncovered"
    );

    // `n2` underflows while the quaternions still differ: `1e-170` squares to zero, so the blend's
    // `‖v‖ = 0` arm fires on a pair that is *not* a point arc. It must still be exact here, which
    // is why that arm is the blend's limit `(1 − t) q₀ + t q₁` and not the constant `q₀`.
    let x0 = SO3::from_quat_unchecked(crate::Quat {
        w: 1.0_f64,
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let x1 = SO3::from_quat_unchecked(crate::Quat {
        w: 1.0_f64,
        x: 1e-170,
        y: 0.0,
        z: 0.0,
    });
    let got = SO3::geodesic(&x0, &x1, 1.0).quat();
    assert_eq!(
        [got.w, got.x].map(f64::to_bits),
        [1.0_f64, 1e-170].map(f64::to_bits),
        "an underflowed `n2` returned the left endpoint at t = 1"
    );
}

/// Extrapolation below the switch belongs to the **provided body**, whose error does not grow with
/// `t`, and not to the blend, whose weights do.
///
/// `LieGroup::geodesic` states that `t` outside `[0, 1]` extrapolates along the same curve, so this
/// is a reachable argument. The blend forms `≈ -(t-1) q₀ + t q₁` there and the two `O(t)` terms
/// cancel: at `θ(d) = 1e-9` it reads a relative `2.6e-10` at `t = 1e6`, about `1.2e6 u`, where the
/// provided body is exact to its own roundings. So `0051`'s clause is `t = 1` and not `t ≥ 1`, and
/// this test is what says so — it compares against `Exp(t·Log Δ)` computed from the same inputs,
/// which is the curve both arms claim to be on.
#[test]
fn extrapolation_below_the_switch_does_not_grow_with_t() {
    let d = SO3Tangent {
        phi: Vector([1e-9, 0.0, 0.0]),
    };
    let x0 = SO3::<f64>::identity();
    let x1 = x0.rplus(&d);
    for t in [1.0, 2.0, 1e3, 1e6] {
        let got = SO3::geodesic(&x0, &x1, t).quat();
        // The curve itself: `x₀ ⊕ t·d`, which is `reference::geodesic`'s body at this `t`.
        let want = x0.rplus(&d.scale(t)).quat();
        let err = [
            got.w - want.w,
            got.x - want.x,
            got.y - want.y,
            got.z - want.z,
        ]
        .iter()
        .map(|e| e * e)
        .sum::<f64>()
        .sqrt()
            / f64::EPSILON;
        assert!(
            err <= 4.0,
            "t = {t}: {err} u from the curve; the blend's cancellation reads 1.2e6 u at t = 1e6"
        );
    }
}

/// `γ(x₀, x₁, 0)` is `x₀` **bit for bit**, which `laws::geodesic`'s `t=0` leg cannot say: `gerr`
/// compares two elements through `Log(x₀⁻¹ x₀)`, and `q* q` leaves about one `u` in the vector
/// part whatever the curve did (1.118 `u`, the floor every quaternion group's law reads there).
///
/// It holds by the arithmetic and not by luck: `d.scale(0)` is `±0` per component, `Exp` of that
/// is `(1, ±0, ±0, ±0)` — `θ² = 0`, so `k` is its series value and `cos(θ/2)` is exactly `1` —
/// and multiplying a quaternion by it adds signed zeros to each component, which is exact
/// (`docs/maths/geodesics.md` GE.7(b), `0045` item 3).
#[test]
fn geodesic_at_zero_is_the_left_endpoint_bit_for_bit() {
    let q =
        |v: [f64; 4]| SO3::from_quat_unchecked(Quat::from_wxyz_normalized(v[0], v[1], v[2], v[3]));
    let cases = [
        ([1.0, 0.0, 0.0, 0.0], [0.3, -0.5, 0.7]),
        ([0.2, 0.3, -0.5, 0.78], [1.0, -2.0, 0.5]),
        ([0.0, 0.6, 0.0, 0.8], [1e-9, 0.0, -1e-9]),
        ([-0.4, 0.1, 0.9, -0.2], [0.0, 0.0, 0.0]),
    ];
    for (w, d) in cases {
        let x0 = q(w);
        let x1 = x0.rplus(&SO3Tangent { phi: Vector(d) });
        let got = SO3::geodesic(&x0, &x1, 0.0).quat();
        let want = x0.quat();
        assert_eq!(
            [got.w, got.x, got.y, got.z].map(f64::to_bits),
            [want.w, want.x, want.y, want.z].map(f64::to_bits),
            "t = 0 moved {w:?} along {d:?}"
        );
    }
    // The one exception, which `LieGroup::geodesic`'s rustdoc states: a `−0.0` in the
    // representation can come back `+0.0`, because `w₀ − x₀·(−0) − …` sums signed zeros and a sum
    // of zeros is negative only when every term is. The value is unchanged; for `w` the bit is
    // load-bearing, since `NUMERICS.md` §3.2 keeps `w = +0` and `Log`'s flip reads it.
    let x0 = q([-0.0, 0.6, 0.0, 0.8]);
    let x1 = x0.rplus(&SO3Tangent {
        phi: Vector([-0.3, -0.5, -0.7]),
    });
    let got = SO3::geodesic(&x0, &x1, 0.0).quat();
    assert_eq!(x0.quat().w.to_bits(), (-0.0_f64).to_bits());
    // A value comparison, not a float `==`: `+0.0` and `-0.0` compare equal and clippy's
    // `float_cmp` is denied.
    assert_eq!(
        got.w.partial_cmp(&x0.quat().w),
        Some(core::cmp::Ordering::Equal),
        "the value is unchanged"
    );
    assert_eq!(
        got.w.to_bits(),
        0.0_f64.to_bits(),
        "a `−0.0` `w` is expected back as `+0.0`; if this ever keeps the sign, the rustdoc caveat \
         on `LieGroup::geodesic` can be dropped"
    );
}

/// `SO3::renormalize` is `Quat::renormalize` on the stored quaternion, bit for bit, at both
/// precisions and through `Dual`: `0056` scores the step as `quat_renormalize` alone, so the
/// group's spelling must not be a second one.
#[test]
fn renormalize_is_the_quaternions_step_bit_for_bit() {
    fn both<S: helicoid_linalg::Real>(q: Quat<S>) -> [Quat<S>; 2] {
        let (mut group, mut quat) = (SO3::from_quat_unchecked(q), q);
        group.renormalize();
        quat.renormalize();
        [group.quat(), quat]
    }
    let parts = |q: Quat<f64>| [q.w, q.x, q.y, q.z].map(f64::to_bits);
    // `‖q‖² − 1` near 2^-27 and 2^-12, inside each precision's band (`NUMERICS.md` §3.6).
    let drift = |e: f64| [0.5 * (1.0 + e), 0.5, -0.5, 0.5 * (1.0 - e / 3.0)];
    let [w, x, y, z] = drift(2f64.powi(-27));
    let [a, b] = both(Quat { w, x, y, z });
    assert_eq!(parts(a), parts(b));
    let [w, x, y, z] = drift(2f64.powi(-12)).map(|c| c as f32);
    let [a, b] = both(Quat { w, x, y, z });
    let parts32 = |q: Quat<f32>| [q.w, q.x, q.y, q.z].map(f32::to_bits);
    assert_eq!(parts32(a), parts32(b));
    let [w, x, y, z] = drift(2f64.powi(-27)).map(|c| Dual::<f64, 1>::variable(c, 0));
    let [a, b] = both(Quat { w, x, y, z });
    let dual =
        |q: Quat<Dual<f64, 1>>| [q.w, q.x, q.y, q.z].map(|c| [c.v.to_bits(), c.d[0].to_bits()]);
    assert_eq!(dual(a), dual(b));
}

/// `0058`'s carried path: `λq`, `λ² = 1 + η`, built by struct literal and moved into
/// [`SO3::from_quat_unchecked`], so no assertion sees it. `η` is read back as `‖λq‖² − 1`.
mod carried {
    use super::*;
    use crate::laws::{unit as ulp, Rng};

    pub(super) fn drifted(q: Quat<f64>, eta: f64) -> (SO3<f64>, f64) {
        let l = (1.0 + eta).sqrt();
        let d = Quat {
            w: q.w * l,
            x: q.x * l,
            y: q.y * l,
            z: q.z * l,
        };
        (SO3::from_quat_unchecked(d), d.norm_sq() - 1.0)
    }

    fn norm(v: [f64; 3]) -> f64 {
        ((v[0] * v[0] + v[1] * v[1]) + v[2] * v[2]).sqrt()
    }

    fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }

    /// The rotation angle between two quaternions, whatever their norms.
    pub(super) fn angle(a: Quat<f64>, b: Quat<f64>) -> f64 {
        let d = a.conjugate() * b;
        2.0 * norm([d.x, d.y, d.z]).atan2(d.w.abs())
    }

    /// `|η|` at the five points of the drift band the record measures, `2^-40` to `2^-26.29`.
    pub(super) const BAND: [f64; 5] = [-40.0, -36.0, -32.0, -28.0, -26.29];

    /// `0058` decision 2's rows for `log`, `act`, `to_matrix`, `Mul` and `inverse`, each against
    /// its first-order statement, across the band. The protocol: 4 000 draws per point of the band
    /// from `Rng(0x0058)`, `η` uniform in `±2^e`, rotations `Exp` of a uniform `[-3, 3)^3`, vectors
    /// uniform `[-1, 1)^3`. The worst errors, in `u`: `log` (relative) 3.959, `act` 12.961,
    /// `to_matrix` 10.885 (both against the unit `act`, so two roundings of the sandwich), `Mul`'s
    /// `η` 6.000. `inverse` keeps `η` to the bit. The bounds are twice those, rounded up.
    #[test]
    fn each_operation_errs_by_its_first_order_term_across_the_drift_band() {
        let u = ulp::<f64>();
        let mut rng = Rng(0x0058);
        let mut worst = [0.0_f64; 4];
        for e in BAND {
            for _ in 0..4000 {
                let tau = |rng: &mut Rng| tan(rng.arr::<3>().map(|c| 3.0 * c));
                let (qa, qb) = (SO3::exp(&tau(&mut rng)), SO3::exp(&tau(&mut rng)));
                let (ea, eb) = (rng.unif() * e.exp2(), rng.unif() * e.exp2());
                let ((da, ea), (db, eb)) = (drifted(qa.quat(), ea), drifted(qb.quat(), eb));

                // `Log` is scale-invariant (§3.2): the tangent of `q/‖q‖`.
                let (l, l0) = (da.log().phi.0, qa.log().phi.0);
                worst[0] = worst[0].max(norm(sub(l, l0)) / norm(l0) / u);

                // §3.3's sandwich at `λq` is `Rv + η(Rv − v)`, §1's matrix `(1 + η)R`.
                let v = rng.arr::<3>();
                let r = qa.act(Vector(v)).0;
                let sandwich = [0, 1, 2].map(|i| r[i] + ea * (r[i] - v[i]));
                worst[1] = worst[1].max(norm(sub(da.act(Vector(v)).0, sandwich)) / norm(v) / u);
                let scaled = r.map(|c| (1.0 + ea) * c);
                let m = (da.to_matrix() * Vector(v)).0;
                worst[2] = worst[2].max(norm(sub(m, scaled)) / norm(v) / u);

                // Norms multiply: `η_ab = η_a + η_b + η_a η_b`.
                let eab = (da * db).quat().norm_sq() - 1.0;
                worst[3] = worst[3].max((eab - (ea + eb + ea * eb)).abs() / u);

                let ei = da.inverse().quat().norm_sq() - 1.0;
                assert!(same(ei, ea), "inverse moved η: {ei} from {ea}");
            }
        }
        let bound = [8.0, 26.0, 22.0, 12.0];
        for (w, b) in worst.iter().zip(bound) {
            assert!(*w <= b, "worst {worst:?} against {bound:?}");
        }
    }

    /// `0058` decision 2's `geodesic` row. Below `r`'s second switch the provided body is
    /// scale-invariant and carries no tilt; the blend tilts toward the longer endpoint by
    /// `sin((1-t)α) sin(tα)/sin α · |η₁ − η₀|` in rotation angle, `α` the half-angle of
    /// `q₀*q₁`. The endpoints stay bit-exact. The protocol: 4 000 draws per point of the band from
    /// `Rng(0x0158)`, `q₀` `Exp` of `[-3, 3)^3`, the relative tangent `[-1.7, 1.7)^3` (the blend)
    /// or `[-1e-4, 1e-4)^3` (the provided body), `t` in `[0, 1)`. Worst: the blend's error is
    /// `1.0027` times its first-order term; the provided body's is `5.612 u` against the geodesic
    /// of the unit endpoints, at every `η`.
    #[test]
    fn the_geodesic_tilts_by_the_norm_difference_only_on_the_blend() {
        let u = ulp::<f64>();
        let mut rng = Rng(0x0158);
        let (mut ratio, mut short) = (0.0_f64, 0.0_f64);
        for e in BAND {
            for k in 0..4000 {
                let scale = if k % 2 == 0 { 1.7 } else { 1e-4 };
                let q0 = SO3::exp(&tan(rng.arr::<3>().map(|c| 3.0 * c)));
                let phi = rng.arr::<3>().map(|c| scale * c);
                let q1 = q0 * SO3::exp(&tan(phi));
                let (e0, e1) = (rng.unif() * e.exp2(), rng.unif() * e.exp2());
                let ((d0, e0), (d1, e1)) = (drifted(q0.quat(), e0), drifted(q1.quat(), e1));
                let t = 0.5 * (rng.unif() + 1.0);
                let got = SO3::geodesic(&d0, &d1, t).quat();
                let err = angle(got, SO3::geodesic(&q0, &q1, t).quat());
                if scale > 1.0 {
                    let a = 0.5 * norm(phi);
                    let term = ((1.0 - t) * a).sin() * (t * a).sin() / a.sin() * (e1 - e0).abs();
                    assert!(err <= 1.01 * term + 32.0 * u, "{err} over the term {term}");
                    if term > 1e3 * u {
                        ratio = ratio.max(err / term);
                    }
                } else {
                    short = short.max(err / u);
                }
                let parts = |q: Quat<f64>| [q.w, q.x, q.y, q.z].map(f64::to_bits);
                assert_eq!(parts(SO3::geodesic(&d0, &d1, 0.0).quat()), parts(d0.quat()));
                // `‖φ‖ < π` here, so `q₀*q₁` has `w > 0` and nothing flips.
                assert_eq!(parts(SO3::geodesic(&d0, &d1, 1.0).quat()), parts(d1.quat()));
            }
        }
        assert!(
            ratio <= 1.01 && short <= 12.0,
            "ratio {ratio}, short {short} u"
        );
    }

    /// The carried path propagates NaN and never panics, in debug or release: no assertion is on
    /// it, which is what lets a consumer whose functions are total in NaN take it.
    #[test]
    fn the_carried_path_propagates_nan() {
        let nan = SO3::from_quat_unchecked(Quat {
            w: f64::NAN,
            x: 0.0,
            y: 0.0,
            z: 0.0,
        });
        let one = SO3::exp(&tan([0.3, -0.2, 0.1]));
        assert!(nan.log().phi.0.iter().all(|c| c.is_nan()));
        assert!(nan
            .act(Vector([1.0, 2.0, 3.0]))
            .0
            .iter()
            .any(|c| c.is_nan()));
        assert!((nan * one).quat().w.is_nan());
        assert!(SO3::geodesic(&one, &nan, 0.5).quat().w.is_nan());
        assert!(SO3::geodesic(&nan, &one, 0.5).quat().w.is_nan());
    }

    /// `0058` decision 3's budget: a chain of one repeated step drifts linearly, and one
    /// `renormalize` takes it back to rounding. The protocol: `10^4` products of each step below;
    /// the drift per product, in `u`, is `0.608`, `0.730` and `0.000`. A repeated step is the
    /// worst case: fresh random steps measured at most `0.024 u` per product over `10^4` to
    /// `10^6` products.
    #[test]
    fn a_repeated_step_drifts_linearly_and_one_step_repairs_it() {
        let u = ulp::<f64>();
        let n = 10_000;
        let mut worst = 0.0_f64;
        for phi in [[0.3, -0.2, 0.5], [1e-3, 2e-3, -1e-3], [0.01, -0.02, 0.015]] {
            let step = SO3::exp(&tan(phi));
            let mut x = step;
            for _ in 1..n {
                x = x * step;
            }
            let eta = x.quat().norm_sq() - 1.0;
            worst = worst.max(eta.abs() / u / f64::from(n));
            x.renormalize();
            let fixed = (x.quat().norm_sq() - 1.0).abs();
            assert!(fixed <= 4.0 * u, "renormalize left η = {fixed} from {eta}");
        }
        assert!(worst <= 1.5, "{worst} u per product");
    }
}
