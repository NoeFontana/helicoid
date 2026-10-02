//! `SEn3Jac`: the dense twins, the layout of `write_dense`, the algebra, `sandwich` and the
//! multiplication counts.
//!
//! Entries are uniform on `[-1, 1)`. A bound is recorded: twice the worst error of
//! `measure_worst_errors` (a seeded splitmix64 stream, 10^6 cases per scalar and `N`; `cargo
//! nextest run --release -p helicoid -E 'test(/dualmat_tests.*measure_worst_errors/)'
//! --run-ignored only --no-capture`), rounded up to an integer, with `N = 1, 2, 3` sharing one
//! bound.
//!
//! An error is `‖fast - twin‖_F` in units of `u`, over `max(‖twin‖_F, 1)` for a product, for
//! `sandwich` and for `apply` (`‖·‖` of a vector there), and the worst over the value and the
//! three derivative lanes of `Dual<f64, 3>` (a plain scalar has zero lanes). The inverse is
//! scaled by its conditioning instead: `κ = ‖M‖_F ‖M⁻¹‖_F` of the dense matrix `M`, and then
//! `‖M⁻¹‖` for the value or `‖M⁻¹‖² ‖dM‖` for a derivative lane (the size of `M⁻¹ dM M⁻¹`). The
//! twin loses `κ u` on `M`, the structured inverse loses it on `A` and again in `A⁻¹ B A⁻¹`, so
//! the two agree to a few `κ u`, not to a few `u`; the recorded figure is that ratio, on the
//! inputs with `κ u <= 1e-3`, where first-order error analysis holds (`κ` up to about `1.7e4` for
//! `f32`). The others are rejected.
//!
//! The worst figures, for `N = 1, 2, 3`: `f64`, product 2.82, 2.37, 2.04, inverse 1.25, 0.92, 0.74,
//! `sandwich` 5.74, 3.87, 2.94, `apply` 3.61, 4.89, 7.20; `f32`, product 2.67, 2.41, 2.14, inverse
//! 1.06, 0.80, 0.83, `sandwich` 6.48, 3.58, 2.93, `apply` 3.62, 4.72, 6.63; `Dual<f64, 3>`, product
//! 3.34, 2.56, 2.27, inverse 1.25, 0.92, 0.74, `sandwich` 8.62, 4.06, 3.47, `apply` 4.05, 4.89,
//! 7.20. `apply` grows with `N` because the dense sum it is compared against grows with `D`, not
//! because the structured form loses more.

// `Rng` shadows `proptest::prelude`'s re-export of the `rand` trait of that name.
use crate::laws::{norm, unit as u, within, worst, Rng, Sample};
use crate::reference::{
    dense, sen3jac_apply, sen3jac_apply_transpose, sen3jac_inverse, sen3jac_mul, sen3jac_sandwich,
};
use crate::{Jac, SEn3Jac, SEn3Tangent, Tangent};
use core::array;
use core::cell::Cell;
use core::ops::{Add, Div, Mul, Neg, Sub};
use helicoid_linalg::{Blend, Dual, Mat3, Matrix, Precision, Real, StridedMut, Vector};
use proptest::prelude::*;
use std::vec::Vec;

struct Bounds {
    mul: f64,
    inverse: f64,
    sandwich: f64,
    apply: f64,
}
const F64: Bounds = Bounds {
    mul: 6.0,
    inverse: 3.0,
    sandwich: 12.0,
    apply: 15.0,
};
const F32: Bounds = Bounds {
    mul: 6.0,
    inverse: 3.0,
    sandwich: 13.0,
    apply: 14.0,
};
const DUAL: Bounds = Bounds {
    mul: 7.0,
    inverse: 3.0,
    sandwich: 18.0,
    apply: 15.0,
};

/// The value and the three derivative lanes of a scalar as `f64`; a plain scalar has zero lanes.
pub(crate) trait Lanes: Sample {
    fn lanes(self) -> [f64; 4];
}
impl Lanes for f64 {
    fn lanes(self) -> [f64; 4] {
        [self, 0.0, 0.0, 0.0]
    }
}
impl Lanes for f32 {
    fn lanes(self) -> [f64; 4] {
        [f64::from(self), 0.0, 0.0, 0.0]
    }
}
impl Lanes for Dual<f64, 3> {
    fn lanes(self) -> [f64; 4] {
        [self.v, self.d[0], self.d[1], self.d[2]]
    }
}

/// The `k`-th of `9 (N + 1)` samples is entry `k % 9` (column-major) of block `k / 9`: `A`, then
/// `B₁, …, B_N`. `mk` builds the scalar from a sample and its index.
fn jac_of<S: Real, const N: usize>(v: &[f64], mk: impl Fn(f64, usize) -> S) -> SEn3Jac<S, N> {
    let mat = |b: usize| {
        Matrix::from_cols(array::from_fn(|c| {
            Vector(array::from_fn(|r| {
                mk(v[9 * b + 3 * c + r], 9 * b + 3 * c + r)
            }))
        }))
    };
    SEn3Jac {
        diag: mat(0),
        col: array::from_fn(|i| mat(i + 1)),
    }
}

/// `Dual`: sample `k` is variable `k mod 3`, so every lane is exercised.
fn jac<S: Lanes, const N: usize>(v: &[f64]) -> SEn3Jac<S, N> {
    jac_of(v, |x, k| S::sample(x, k % 3))
}

fn lane<S: Lanes, const D: usize>(m: &[[S; D]; D], l: usize) -> [[f64; D]; D] {
    array::from_fn(|r| array::from_fn(|c| m[r][c].lanes()[l]))
}

/// `‖m‖_F` of each of the four lanes of a dense matrix.
fn norms<S: Lanes, const D: usize>(m: &[[S; D]; D]) -> [f64; 4] {
    array::from_fn(|l| norm(lane(m, l).as_flattened()))
}

/// `‖a - b‖_F` of each of the four lanes, in units of `u`.
fn diffs<S: Lanes, const D: usize>(a: &[[S; D]; D], b: &[[S; D]; D]) -> [f64; 4] {
    array::from_fn(|l| {
        let (x, y) = (lane(a, l), lane(b, l));
        let d: [f64; D] =
            array::from_fn(|r| norm(&array::from_fn::<f64, D, _>(|c| x[r][c] - y[r][c])));
        norm(&d) / u::<S>()
    })
}

/// The worst over the four lanes of `‖fast - twin‖_F / scale`, where a lane that agrees exactly is
/// zero whatever its scale (a plain scalar has zero lanes) and a NaN difference stays NaN.
fn worst_lane(d: [f64; 4], scale: [f64; 4]) -> f64 {
    (0..4).fold(0.0, |w, l| {
        worst(w, if d[l] == 0.0 { 0.0 } else { d[l] / scale[l] })
    })
}

/// `fast.mul(b)` against the dense product.
fn mul_twin<S: Lanes, const N: usize, const D: usize>(v: &[f64]) -> f64 {
    let (a, b) = (jac::<S, N>(&v[..9 * (N + 1)]), jac(&v[9 * (N + 1)..]));
    let mut twin = [[S::zero(); D]; D];
    sen3jac_mul::<S, N, D>(
        &a,
        &b,
        &mut StridedMut::row_major(twin.as_flattened_mut(), D, D),
    );
    let fast = dense::<S, _, _, D>(&a.mul(&b));
    worst_lane(diffs(&fast, &twin), norms(&twin).map(|n| n.max(1.0)))
}

fn to_f64<S: Lanes, const N: usize>(j: &SEn3Jac<S, N>) -> SEn3Jac<f64, N> {
    let m = |m: &Mat3<S>| {
        Matrix::from_cols(array::from_fn(|c| {
            Vector(array::from_fn(|r| m.get(r, c).lanes()[0]))
        }))
    };
    SEn3Jac {
        diag: m(&j.diag),
        col: j.col.map(|b| m(&b)),
    }
}

/// `κ_F` of the dense matrix of `j`, with `M⁻¹` from the twin in `f64`.
fn kappa<S: Lanes, const N: usize, const D: usize>(j: &SEn3Jac<S, N>) -> f64 {
    let j64 = to_f64(j);
    let mut inv = [[0.0; D]; D];
    sen3jac_inverse::<f64, N, D>(
        &j64,
        &mut StridedMut::row_major(inv.as_flattened_mut(), D, D),
    );
    norm(dense::<f64, _, _, D>(&j64).as_flattened()) * norm(inv.as_flattened())
}

/// The error of `inverse` against the twin over `κ u`, or `None` where `κ u > 1e-3`. The value
/// is scaled by `‖M⁻¹‖ κ` and a derivative lane by `‖M⁻¹‖² ‖dM‖ κ`, the size of `M⁻¹ dM M⁻¹`
/// that the same rounding acts on.
fn inverse_twin<S: Lanes, const N: usize, const D: usize>(v: &[f64]) -> Option<f64> {
    let j = jac::<S, N>(v);
    let k = kappa::<S, N, D>(&j);
    let fits = k * u::<S>() <= 1e-3;
    if !fits {
        return None;
    }
    let mut twin = [[S::zero(); D]; D];
    sen3jac_inverse::<S, N, D>(
        &j,
        &mut StridedMut::row_major(twin.as_flattened_mut(), D, D),
    );
    let (fast, m) = (dense::<S, _, _, D>(&j.inverse()), dense::<S, _, _, D>(&j));
    let (inv, dm) = (norms(&twin)[0], norms(&m));
    let scale = array::from_fn(|l| k * inv * if l == 0 { 1.0 } else { inv * dm[l] });
    Some(worst_lane(diffs(&fast, &twin), scale))
}

/// `fast.sandwich(Σ)` against the twin, on `A + εB` and a dense `Σ` (not symmetric: the formula is
/// linear in `Σ`).
fn sandwich_twin<S: Lanes, const N: usize, const D: usize>(v: &[f64]) -> f64 {
    let (j, s) = (jac::<S, N>(&v[..9 * (N + 1)]), &v[9 * (N + 1)..]);
    let cov = Matrix::<S, D, D>::from_cols(array::from_fn(|c| {
        Vector(array::from_fn(|r| S::sample(s[D * c + r], (D * c + r) % 3)))
    }));
    let mut twin = [[S::zero(); D]; D];
    sen3jac_sandwich::<S, N, D>(
        &j,
        &cov,
        &mut StridedMut::row_major(twin.as_flattened_mut(), D, D),
    );
    let got = j.sandwich::<D>(&cov);
    let got: [[S; D]; D] = array::from_fn(|r| array::from_fn(|c| got.get(r, c)));
    worst_lane(diffs(&got, &twin), norms(&twin).map(|n| n.max(1.0)))
}

/// `‖a - b‖` of each of the four lanes of a dense vector, in units of `u`.
fn vdiffs<S: Lanes, const D: usize>(a: &[S; D], b: &[S; D]) -> [f64; 4] {
    array::from_fn(|l| {
        let d: [f64; D] = array::from_fn(|i| a[i].lanes()[l] - b[i].lanes()[l]);
        norm(&d) / u::<S>()
    })
}

/// `‖a‖` of each of the four lanes of a dense vector.
fn vnorms<S: Lanes, const D: usize>(a: &[S; D]) -> [f64; 4] {
    array::from_fn(|l| norm(&array::from_fn::<f64, D, _>(|i| a[i].lanes()[l])))
}

/// `apply` and `apply_transpose` against their twins, the dense `J x` and `Jᵀ x`.
///
/// The entries must not be exactly representable for this to say anything: the structured forms do
/// not associate the sum the way a dense row does, and
/// `apply_and_apply_transpose_are_the_dense_matrix_and_its_transpose` compares them on integers,
/// which hides it.
fn apply_twin<S: Lanes, const N: usize, const D: usize>(v: &[f64]) -> f64 {
    let w = 9 * (N + 1);
    let j = jac::<S, N>(&v[..w]);
    let t =
        SEn3Tangent::<S, N>::read_dense(&array::from_fn::<S, D, _>(|i| S::sample(v[w + i], i % 3)));
    let mut x = [S::zero(); D];
    t.write_dense(&mut x);
    let (mut a, mut b) = ([S::zero(); D], [S::zero(); D]);
    sen3jac_apply::<S, N, D>(&j, &x, &mut a);
    sen3jac_apply_transpose::<S, N, D>(&j, &x, &mut b);
    let run = |got: SEn3Tangent<S, N>, twin: [S; D]| {
        let mut y = [S::zero(); D];
        got.write_dense(&mut y);
        worst_lane(vdiffs(&y, &twin), vnorms(&twin).map(|n| n.max(1.0)))
    };
    worst(run(j.apply(&t), a), run(j.apply_transpose(&t), b))
}

/// The four lanes of every entry of a dense matrix as bits, signed zeros canonicalized: a
/// structural zero is `+0` from `write_dense` and `-0` from a negated dense matrix.
fn lane_bits<S: Lanes, const D: usize>(m: &[[S; D]; D]) -> [[[u64; 4]; D]; D] {
    m.map(|row| row.map(|x| x.lanes().map(canon)))
}

/// `neg` and `Blend` carry every derivative lane.
///
/// Both are per-entry, so the dense comparison is exact whatever the input and no bound is
/// recorded. They are tested over `Dual` and not only `f64` because the mask is a `bool` for every
/// `Real` of this crate: at `f64` `Blend` cannot do worse than pick one of two whole structs,
/// while at `Dual` it runs `S::select` per value *and* per lane, and `Real::branch` — what every
/// future switch-point arm of `SEn3` goes through — is built on it.
fn neg_and_blend_carry_lanes<S: Lanes, const N: usize, const D: usize>(v: &[f64]) -> bool {
    let w = 9 * (N + 1);
    let (a, b) = (jac::<S, N>(&v[..w]), jac::<S, N>(&v[w..]));
    let (yes, no) = (S::zero().lt(S::one()), S::one().lt(S::zero()));
    let minus = dense::<S, _, _, D>(&a).map(|row| row.map(|x| -x));
    lane_bits(&dense::<S, _, _, D>(&SEn3Jac::blend(yes, a, b)))
        == lane_bits(&dense::<S, _, _, D>(&a))
        && lane_bits(&dense::<S, _, _, D>(&SEn3Jac::blend(no, a, b)))
            == lane_bits(&dense::<S, _, _, D>(&b))
        && lane_bits(&dense::<S, _, _, D>(&a.neg())) == lane_bits(&minus)
}

fn entries(n: usize) -> impl Strategy<Value = Vec<f64>> {
    proptest::collection::vec(-1.0_f64..1.0, n)
}

macro_rules! props {
    ($m:ident, $S:ty, $B:ident) => {
        mod $m {
            use super::*;
            type S = $S;
            props!(@n n1, 1, 6, $B);
            props!(@n n2, 2, 9, $B);
            props!(@n n3, 3, 12, $B);
        }
    };
    (@n $m:ident, $N:literal, $D:literal, $B:ident) => {
        mod $m {
            use super::*;
            proptest! {
                #[test]
                fn sen3jac_mul_matches_reference(v in entries(18 * ($N + 1))) {
                    within(mul_twin::<S, $N, $D>(&v), $B.mul)?;
                }
                #[test]
                fn sen3jac_inverse_matches_reference(v in entries(9 * ($N + 1))) {
                    let Some(r) = inverse_twin::<S, $N, $D>(&v) else {
                        return Err(TestCaseError::reject("kappa u > 1e-3"));
                    };
                    within(r, $B.inverse)?;
                }
                #[test]
                fn sen3jac_sandwich_matches_reference(v in entries(9 * ($N + 1) + $D * $D)) {
                    within(sandwich_twin::<S, $N, $D>(&v), $B.sandwich)?;
                }
                #[test]
                fn sen3jac_apply_matches_reference(v in entries(9 * ($N + 1) + $D)) {
                    within(apply_twin::<S, $N, $D>(&v), $B.apply)?;
                }
                #[test]
                fn neg_and_blend_carry_every_lane(v in entries(18 * ($N + 1))) {
                    prop_assert!(neg_and_blend_carry_lanes::<S, $N, $D>(&v));
                }
            }

            /// The worst error of every law over 10^6 seeded cases: mul, inverse, sandwich, apply.
            #[test]
            #[ignore = "measurement: prints the figures the bounds are recorded from"]
            #[allow(clippy::print_stdout)]
            fn measure_worst_errors() {
                let mut rng = Rng(0x0123_4567_89AB_CDEF);
                // A separate stream for `apply`, so that adding it left mul, inverse and sandwich
                // drawing exactly what they drew when their figures were recorded.
                let mut ap = Rng(0xDEAD_BEEF_0BAD_F00D);
                let mut w = [0.0_f64; 4];
                let mut up = |i: usize, v: f64| w[i] = worst(w[i], v);
                for _ in 0..1_000_000 {
                    up(0, mul_twin::<S, $N, $D>(&rng.vec(18 * ($N + 1))));
                    if let Some(r) = inverse_twin::<S, $N, $D>(&rng.vec(9 * ($N + 1))) {
                        up(1, r);
                    }
                    up(2, sandwich_twin::<S, $N, $D>(&rng.vec(9 * ($N + 1) + $D * $D)));
                    up(3, apply_twin::<S, $N, $D>(&ap.vec(9 * ($N + 1) + $D)));
                }
                std::println!("{} {} {:.2?}", std::any::type_name::<S>(), stringify!($m), w);
            }
        }
    };
}

props!(as_f64, f64, F64);
props!(as_f32, f32, F32);
props!(as_dual, Dual<f64, 3>, DUAL);

fn canon(x: f64) -> u64 {
    // `-0 + 0 = +0`: the exact tests are about values, and a sum of products has zeros of either sign.
    (x + 0.0).to_bits()
}

fn canon_dense<S: Lanes, const D: usize>(m: &[[S; D]; D]) -> [[u64; D]; D] {
    lane(m, 0).map(|row| row.map(canon))
}

fn rows3<S: Real>(r: [[f64; 3]; 3]) -> Mat3<S> {
    Matrix::from_rows(r.map(|row| Vector(row.map(S::lit))))
}

/// `A = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]` and `B_i = 10 i + A`: distinct integers, so every dense
/// entry says where it came from and every product and sum below is exact. `A` is singular.
fn ints<S: Real, const N: usize>() -> SEn3Jac<S, N> {
    let m = |o: f64| {
        Matrix::from_rows(array::from_fn(|r| {
            Vector(array::from_fn(|c| S::lit(o + (1 + 3 * r + c) as f64)))
        }))
    };
    SEn3Jac {
        diag: m(0.0),
        col: array::from_fn(|i| m(10.0 * (i + 1) as f64)),
    }
}

/// `ints` with the unimodular `A = [[1, 2, 3], [0, 1, 4], [0, 0, 1]]`: its inverse is integral.
fn unimodular<S: Real, const N: usize>() -> SEn3Jac<S, N> {
    let a = rows3([[1.0, 2.0, 3.0], [0.0, 1.0, 4.0], [0.0, 0.0, 1.0]]);
    SEn3Jac { diag: a, ..ints() }
}

fn eye<const D: usize>() -> [[u64; D]; D] {
    array::from_fn(|r| array::from_fn(|c| if r == c { 1.0_f64 } else { 0.0 }.to_bits()))
}

#[test]
fn write_dense_is_the_matrix_of_numerics_2_2() {
    #[rustfmt::skip]
    let n1: [[f64; 6]; 6] = [
        [ 1.0,  2.0,  3.0, 0.0, 0.0, 0.0],
        [ 4.0,  5.0,  6.0, 0.0, 0.0, 0.0],
        [ 7.0,  8.0,  9.0, 0.0, 0.0, 0.0],
        [11.0, 12.0, 13.0, 1.0, 2.0, 3.0],
        [14.0, 15.0, 16.0, 4.0, 5.0, 6.0],
        [17.0, 18.0, 19.0, 7.0, 8.0, 9.0],
    ];
    #[rustfmt::skip]
    let n2: [[f64; 9]; 9] = [
        [ 1.0,  2.0,  3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [ 4.0,  5.0,  6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [ 7.0,  8.0,  9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [11.0, 12.0, 13.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0],
        [14.0, 15.0, 16.0, 4.0, 5.0, 6.0, 0.0, 0.0, 0.0],
        [17.0, 18.0, 19.0, 7.0, 8.0, 9.0, 0.0, 0.0, 0.0],
        [21.0, 22.0, 23.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0],
        [24.0, 25.0, 26.0, 0.0, 0.0, 0.0, 4.0, 5.0, 6.0],
        [27.0, 28.0, 29.0, 0.0, 0.0, 0.0, 7.0, 8.0, 9.0],
    ];
    let bits = |m: &[[f64; 6]; 6]| m.map(|r| r.map(f64::to_bits));
    assert_eq!(
        canon_dense(&dense::<f64, _, _, 6>(&ints::<f64, 1>())),
        bits(&n1)
    );
    let bits = |m: &[[f64; 9]; 9]| m.map(|r| r.map(f64::to_bits));
    assert_eq!(
        canon_dense(&dense::<f64, _, _, 9>(&ints::<f64, 2>())),
        bits(&n2)
    );
    // `N = 3`: `B_3` (entries 31 to 39) alone fills block row 3 of block column 0.
    let d = lane(&dense::<f64, _, _, 12>(&ints::<f64, 3>()), 0);
    for (r, c) in [(9, 0), (11, 2), (10, 1)] {
        assert_eq!(
            d[r][c].to_bits(),
            (31.0 + (3 * (r - 9) + c) as f64).to_bits()
        );
    }
    for (r, c) in [(9, 3), (9, 8), (10, 6), (3, 9), (6, 11), (0, 9)] {
        assert_eq!(canon(d[r][c]), 0.0_f64.to_bits(), "({r}, {c})");
    }
    for (k, row) in d.iter().enumerate() {
        assert_eq!(
            row[k].to_bits(),
            f64::from(1 + 4 * (k as u32 % 3)).to_bits()
        );
    }
}

#[test]
fn write_dense_fills_the_view_and_only_the_view() {
    let j = ints::<f64, 1>();
    let mut padded = [f64::NAN; 6 * 8];
    j.write_dense(&mut StridedMut::with_strides(&mut padded, 6, 6, 1, 8));
    let mut rowmajor = [f64::NAN; 36];
    j.write_dense(&mut StridedMut::row_major(&mut rowmajor, 6, 6));
    let want = dense::<f64, _, _, 6>(&j);
    for (i, v) in padded.iter().enumerate() {
        let (r, c) = (i % 8, i / 8);
        assert!(if r < 6 {
            v.to_bits() == want[r][c].to_bits()
        } else {
            v.is_nan()
        });
    }
    for (i, v) in rowmajor.iter().enumerate() {
        assert_eq!(v.to_bits(), want[i / 6][i % 6].to_bits());
    }
}

/// The tangent `[1, 2, ..., 3 + 3N]`.
fn ramp<S: Real, const N: usize, const D: usize>() -> SEn3Tangent<S, N> {
    SEn3Tangent::read_dense(&array::from_fn::<S, D, _>(|i| S::lit((i + 1) as f64)))
}

/// On the integer `ints`, where every product and partial sum is exact whatever the association,
/// `apply` and `apply_transpose` are their twins entry for entry. The association *is* visible on
/// other inputs; `sen3jac_apply_matches_reference` is what bounds it.
#[test]
fn apply_and_apply_transpose_are_the_dense_matrix_and_its_transpose() {
    fn run<const N: usize, const D: usize>() {
        let (j, t) = (ints::<f64, N>(), ramp::<f64, N, D>());
        let mut x = [0.0; D];
        t.write_dense(&mut x);
        let (mut a, mut b) = ([0.0; D], [0.0; D]);
        sen3jac_apply::<f64, N, D>(&j, &x, &mut a);
        sen3jac_apply_transpose::<f64, N, D>(&j, &x, &mut b);
        let (mut y, mut z) = ([0.0; D], [0.0; D]);
        j.apply(&t).write_dense(&mut y);
        j.apply_transpose(&t).write_dense(&mut z);
        assert_eq!(y.map(canon), a.map(canon));
        assert_eq!(z.map(canon), b.map(canon));
    }
    run::<1, 6>();
    run::<2, 9>();
    run::<3, 12>();
}

#[test]
fn integer_products_and_inverses_are_exact() {
    fn run<const N: usize, const D: usize>() {
        let (a, b) = (unimodular::<f64, N>(), ints::<f64, N>());
        let c = SEn3Jac {
            diag: rows3([[1.0, 0.0, 0.0], [5.0, 1.0, 0.0], [-2.0, 3.0, 1.0]]),
            ..b
        };
        let d = |j: &SEn3Jac<f64, N>| canon_dense(&dense::<f64, _, _, D>(j));
        // Closure and the product formula: the dense product has the structure and the values.
        for (x, y) in [(a, b), (b, a), (a, c), (b, b)] {
            let mut twin = [[0.0; D]; D];
            let mut view = StridedMut::row_major(twin.as_flattened_mut(), D, D);
            sen3jac_mul::<f64, N, D>(&x, &y, &mut view);
            assert_eq!(canon_dense(&twin), d(&x.mul(&y)));
        }
        assert_eq!(d(&a.mul(&b).mul(&c)), d(&a.mul(&b.mul(&c))));
        let id = SEn3Jac::<f64, N>::identity();
        assert_eq!((d(&id.mul(&b)), d(&b.mul(&id))), (d(&b), d(&b)));
        assert_eq!(d(&id), eye());
        // `A` is unimodular, so `A⁻¹` and `A⁻¹ B A⁻¹` are integral and `J J⁻¹ = J⁻¹ J = I`.
        assert_eq!(
            (d(&a.mul(&a.inverse())), d(&a.inverse().mul(&a))),
            (eye(), eye())
        );
        assert_eq!(d(&a.inverse().inverse()), d(&a));
        let minus = dense::<f64, _, _, D>(&b).map(|row| row.map(|x| -x));
        assert_eq!(d(&b.neg()), canon_dense(&minus));
    }
    run::<1, 6>();
    run::<2, 9>();
    run::<3, 12>();
}

/// A permutation `A` (`det = -1`, `A⁻¹ = A`) has a zero in the pivot position `(0, 0)`, so the
/// twin must exchange rows, and the inverse is `[[A, 0], [-A B A, A]]`, worked by hand for
/// `B = [[2, 3, 4], [5, 6, 7], [8, 9, 10]]`. Both inverses are exact here.
#[test]
fn the_inverses_of_a_permutation_are_the_hand_inverse() {
    let a = rows3([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
    let b = rows3([[2.0, 3.0, 4.0], [5.0, 6.0, 7.0], [8.0, 9.0, 10.0]]);
    let j = SEn3Jac { diag: a, col: [b] };
    #[rustfmt::skip]
    let want: [[f64; 6]; 6] = [
        [ 0.0,  1.0,   0.0, 0.0, 0.0, 0.0],
        [ 1.0,  0.0,   0.0, 0.0, 0.0, 0.0],
        [ 0.0,  0.0,   1.0, 0.0, 0.0, 0.0],
        [-6.0, -5.0,  -7.0, 0.0, 1.0, 0.0],
        [-3.0, -2.0,  -4.0, 1.0, 0.0, 0.0],
        [-9.0, -8.0, -10.0, 0.0, 0.0, 1.0],
    ];
    let mut twin = [[0.0; 6]; 6];
    sen3jac_inverse::<f64, 1, 6>(
        &j,
        &mut StridedMut::row_major(twin.as_flattened_mut(), 6, 6),
    );
    // Partial pivoting takes the pivot of column 0 from `B` (entries up to 8), so the twin divides
    // by non-powers of two and is not exact: its largest entry error is 4 u (bound: twice that);
    // the structured inverse is exact.
    let gap = (0..36).map(|i| (twin[i / 6][i % 6] - want[i / 6][i % 6]).abs());
    assert!(gap.fold(0.0, f64::max) <= 8.0 * f64::EPSILON / 2.0);
    let want = want.map(|r| r.map(f64::to_bits));
    assert_eq!(canon_dense(&dense::<f64, _, _, 6>(&j.inverse())), want);
}

/// The same for `sandwich`: on integers it is its twin entry for entry.
#[test]
fn sandwich_of_integers_is_exact() {
    fn run<const N: usize, const D: usize>() {
        let j = ints::<f64, N>();
        // A non-symmetric integer `Σ`, entry `(i + 2k) mod 7 - 3`.
        let cov = Matrix::<f64, D, D>::from_cols(array::from_fn(|c| {
            Vector(array::from_fn(|r| ((r + 2 * c) % 7) as f64 - 3.0))
        }));
        let mut twin = [[0.0; D]; D];
        sen3jac_sandwich::<f64, N, D>(
            &j,
            &cov,
            &mut StridedMut::row_major(twin.as_flattened_mut(), D, D),
        );
        let got = j.sandwich::<D>(&cov);
        for (i, row) in twin.iter().enumerate() {
            for (k, want) in row.iter().enumerate() {
                assert_eq!(canon(got.get(i, k)), canon(*want), "({i}, {k})");
            }
        }
    }
    run::<1, 6>();
    run::<2, 9>();
    run::<3, 12>();
}

std::thread_local! {
    static MULS: Cell<u32> = const { Cell::new(0) };
}

fn tick() {
    MULS.with(|c| c.set(c.get() + 1));
}

/// An `f64` that counts its multiplications, so the cost of `0005` is measured, not asserted.
#[derive(Clone, Copy, Debug)]
struct Counted(f64);
impl Add for Counted {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self(self.0 + o.0)
    }
}
impl Sub for Counted {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self(self.0 - o.0)
    }
}
impl Div for Counted {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        Self(self.0 / o.0)
    }
}
impl Neg for Counted {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0)
    }
}
impl Mul for Counted {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        tick();
        Self(self.0 * o.0)
    }
}
impl Real for Counted {
    type Mask = bool;
    const PRECISION: Precision = Precision::F64;
    fn lit(x: f64) -> Self {
        Self(x)
    }
    fn zero() -> Self {
        Self(0.0)
    }
    fn one() -> Self {
        Self(1.0)
    }
    fn lt(self, r: Self) -> bool {
        self.0 < r.0
    }
    fn le(self, r: Self) -> bool {
        self.0 <= r.0
    }
    fn select(m: bool, t: Self, f: Self) -> Self {
        if m {
            t
        } else {
            f
        }
    }
    fn sqrt(self) -> Self {
        Self(self.0.sqrt())
    }
    fn sin_cos(self) -> (Self, Self) {
        let (s, c) = self.0.sin_cos();
        (Self(s), Self(c))
    }
    fn atan2(self, x: Self) -> Self {
        Self(self.0.atan2(x.0))
    }
    fn abs(self) -> Self {
        Self(self.0.abs())
    }
    fn copysign(self, s: Self) -> Self {
        Self(self.0.copysign(s.0))
    }
    fn value_f64(self) -> f64 {
        self.0
    }
}

/// The multiplications `f` performs.
fn count<R>(f: impl FnOnce() -> R) -> u32 {
    MULS.with(|c| c.set(0));
    let _ = f();
    MULS.with(Cell::get)
}

#[test]
fn the_multiplication_counts_of_0005() {
    fn run<const N: usize, const D: usize>() {
        let (a, b) = (unimodular::<Counted, N>(), ints::<Counted, N>());
        let t = ramp::<Counted, N, D>();
        let n = N as u32;
        let adj = count(|| a.diag.inverse_adj());
        let mut out = [Counted(0.0); 200];
        let twin = count(|| {
            sen3jac_mul::<Counted, N, D>(&a, &b, &mut StridedMut::row_major(&mut out, D, D));
        });
        assert_eq!(count(|| a.mul(&b)), 27 + 54 * n);
        assert_eq!(count(|| a.inverse()) - adj, 54 * n);
        assert_eq!(count(|| a.apply(&t)), 9 + 18 * n);
        assert_eq!(count(|| a.apply_transpose(&t)), 9 + 18 * n);
        assert_eq!(twin, 27 * (n + 1).pow(3));
    }
    // 81 against 216, 135 against 729, 189 against 1728: the table of `docs/maths/se3.md` SE.11.
    run::<1, 6>();
    run::<2, 9>();
    run::<3, 12>();
}

/// The value parts of `mul`, `inverse`, `neg`, `apply`, `apply_transpose` and `sandwich`, as bits.
fn probe<S: Lanes, const N: usize, const D: usize>(
    v: &[f64],
    mk: &dyn Fn(f64, usize) -> S,
) -> Vec<u64> {
    let w = 9 * (N + 1);
    let (a, b) = (
        jac_of::<S, N>(&v[..w], mk),
        jac_of::<S, N>(&v[w..2 * w], mk),
    );
    let t = SEn3Tangent::<S, N>::read_dense(&array::from_fn::<S, D, _>(|i| mk(v[2 * w + i], i)));
    let db = dense::<S, _, _, D>(&b);
    let cov =
        Matrix::<S, D, D>::from_cols(array::from_fn(|c| Vector(array::from_fn(|r| db[r][c]))));
    let s = j_sandwich(&a, &cov);
    let (mut x, mut y) = ([S::zero(); D], [S::zero(); D]);
    a.apply(&t).write_dense(&mut x);
    a.apply_transpose(&t).write_dense(&mut y);
    let mut out: Vec<u64> = [a.mul(&b), a.inverse(), a.neg()]
        .iter()
        .flat_map(|j| {
            dense::<S, _, _, D>(j)
                .as_flattened()
                .iter()
                .map(|e| e.lanes()[0].to_bits())
                .collect::<Vec<_>>()
        })
        .collect();
    out.extend(x.iter().chain(&y).map(|e| e.lanes()[0].to_bits()));
    out.extend(s.iter().map(|e| e.lanes()[0].to_bits()));
    out
}

fn j_sandwich<S: Real, const N: usize, const D: usize>(
    j: &SEn3Jac<S, N>,
    cov: &Matrix<S, D, D>,
) -> Vec<S> {
    let s = j.sandwich::<D>(cov);
    (0..D * D).map(|i| s.get(i / D, i % D)).collect()
}

proptest! {
    /// Over `Dual<f64, 3>`, with variable or poisoned (NaN, infinite, huge) derivative lanes, the
    /// value parts are the plain `f64` results, bit for bit (`PHASE3.md` §9).
    #[test]
    fn sen3jac_dual_value_is_plain_value(v in entries(2 * 27 + 12)) {
        type D3 = Dual<f64, 3>;
        let plain = probe::<f64, 2, 9>(&v, &|x, _| x);
        let vars = probe::<D3, 2, 9>(&v, &|x, k| Dual::variable(x, k % 3));
        let poisoned = probe::<D3, 2, 9>(&v, &|x, _| Dual { v: x, d: [f64::NAN, f64::INFINITY, -1e300] });
        prop_assert_eq!(&plain, &vars);
        prop_assert_eq!(&plain, &poisoned);
    }
}

/// `A` with a zero row: singular, and the third pivot of the twin is exactly `0`.
fn singular() -> SEn3Jac<f64, 1> {
    let a = rows3([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [0.0, 0.0, 0.0]]);
    SEn3Jac { diag: a, col: [a] }
}

/// `A = 1e120 I`: `det A = 1e360` overflows to `+∞` while `adj A = 1e240` does not, so every entry
/// of `A⁻¹` is a finite `0`.
fn overflowing() -> SEn3Jac<f64, 1> {
    let a = rows3([[1e120, 0.0, 0.0], [0.0, 1e120, 0.0], [0.0, 0.0, 1e120]]);
    SEn3Jac { diag: a, col: [a] }
}

/// A domain is a `debug_assert!` (D11): the debug profile that `just test` runs panics ...
#[cfg(debug_assertions)]
mod out_of_domain {
    use super::*;

    #[test]
    #[should_panic(expected = "det A")]
    fn inverse_of_a_singular_a() {
        let _ = singular().inverse();
    }

    #[test]
    #[should_panic(expected = "infinite")]
    fn inverse_where_det_a_overflows() {
        let _ = overflowing().inverse();
    }

    #[test]
    #[should_panic(expected = "pivot")]
    fn the_twin_of_a_singular_a() {
        let mut out = [0.0; 36];
        sen3jac_inverse::<f64, 1, 6>(&singular(), &mut StridedMut::row_major(&mut out, 6, 6));
    }

    // A view larger than `DOF x DOF` is no strided-access panic, so only the assertion sees it.
    #[test]
    #[should_panic(expected = "Jac::write_dense")]
    fn write_dense_into_a_larger_view() {
        let mut buf = [0.0; 49];
        ints::<f64, 1>().write_dense(&mut StridedMut::col_major(&mut buf, 7, 7));
    }

    // The same for each twin's view assertion, which a smaller view never reaches either.
    #[test]
    #[should_panic(expected = "sen3jac_mul")]
    fn the_twin_of_mul_into_a_larger_view() {
        let (j, mut buf) = (ints::<f64, 1>(), [0.0; 49]);
        sen3jac_mul::<f64, 1, 6>(&j, &j, &mut StridedMut::col_major(&mut buf, 7, 7));
    }

    #[test]
    #[should_panic(expected = "sen3jac_inverse")]
    fn the_twin_of_inverse_into_a_larger_view() {
        let mut buf = [0.0; 49];
        sen3jac_inverse::<f64, 1, 6>(&unimodular(), &mut StridedMut::col_major(&mut buf, 7, 7));
    }
}

/// ... and a release build never panics: both inverses divide by zero, and where `det A`
/// overflows `inverse` returns a finite, wrong zero (`Mat3::inverse_adj`).
#[cfg(not(debug_assertions))]
#[test]
fn out_of_domain_does_not_panic_in_release() {
    let wrong = dense::<f64, _, _, 6>(&overflowing().inverse());
    assert!(wrong.as_flattened().iter().all(|&x| x == 0.0));
    let j = singular();
    let fast = dense::<f64, _, _, 6>(&j.inverse());
    let mut twin = [[0.0; 6]; 6];
    sen3jac_inverse::<f64, 1, 6>(
        &j,
        &mut StridedMut::row_major(twin.as_flattened_mut(), 6, 6),
    );
    assert!(fast.as_flattened().iter().any(|x| !x.is_finite()));
    assert!(twin.as_flattened().iter().any(|x| !x.is_finite()));
}
