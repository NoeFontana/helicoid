//! `SEn3Tangent`: the dense order, its operations against the dense arrays, the twist converters
//! and the asserted domains.
//!
//! Every comparison is on bits, over the value and the three derivative lanes: each operation is
//! the same `S` operation as its dense counterpart, so nothing is measured, nothing is bounded.

use crate::dualmat_tests::Lanes;
use crate::{SEn3Tangent, Tangent, Twist};
use core::array;
use helicoid_linalg::{Blend, Dual, Real, Vector};
use proptest::prelude::*;

fn tangent<S: Lanes, const N: usize, const D: usize>(v: &[f64]) -> SEn3Tangent<S, N> {
    SEn3Tangent::read_dense(&array::from_fn::<S, D, _>(|i| S::sample(v[i], i % 3)))
}

fn dense_of<S: Real, const N: usize, const D: usize>(t: &SEn3Tangent<S, N>) -> [S; D] {
    let mut out = [S::zero(); D];
    t.write_dense(&mut out);
    out
}

fn bits<S: Lanes, const D: usize>(a: &[S; D]) -> [[u64; 4]; D] {
    a.map(|x| x.lanes().map(f64::to_bits))
}

/// The operations of `Tangent` against the same operations on the dense arrays, and `dot` against
/// one left-to-right sum of all `3 + 3N` products seeded with `+0`, the seed `0025` makes
/// normative for the provided `dot`.
fn ops_match_dense<S: Lanes, const N: usize, const D: usize>(v: &[f64]) -> bool {
    let (a, b) = (
        tangent::<S, N, D>(&v[..D]),
        tangent::<S, N, D>(&v[D..2 * D]),
    );
    let k = S::sample(v[2 * D], 0);
    let (da, db) = (dense_of::<S, N, D>(&a), dense_of::<S, N, D>(&b));
    let dot = (0..D).fold(S::zero(), |acc, i| acc + da[i] * db[i]);
    // Seeded with `k`, not `+0`: nothing here composes tangents, so `dot` alone would leave
    // `dot_acc`'s accumulator — the operation `0025` makes the required one — unexecuted.
    let dot_k = (0..D).fold(k, |acc, i| acc + da[i] * db[i]);
    let round = dense_of::<S, N, D>(&SEn3Tangent::<S, N>::read_dense(&da));
    // `Blend` selects one whole tangent. It is here, not in a test of its own at `f64`, because
    // the mask is a `bool` for every `Real` of this crate: a lost derivative lane of `Dual` is the
    // only way the impl can differ from a copy, and `Real::branch` is built on it.
    let (yes, no) = (S::zero().lt(S::one()), S::one().lt(S::zero()));
    bits(&dense_of::<S, N, D>(&SEn3Tangent::blend(yes, a, b))) == bits(&da)
        && bits(&dense_of::<S, N, D>(&SEn3Tangent::blend(no, a, b))) == bits(&db)
        && bits(&dense_of::<S, N, D>(&a.add(&b))) == bits(&array::from_fn(|i| da[i] + db[i]))
        && bits(&dense_of::<S, N, D>(&a.sub(&b))) == bits(&array::from_fn(|i| da[i] - db[i]))
        && bits(&dense_of::<S, N, D>(&a.neg())) == bits(&da.map(|x| -x))
        && bits(&dense_of::<S, N, D>(&a.scale(k))) == bits(&da.map(|x| x * k))
        && bits(&dense_of::<S, N, D>(&SEn3Tangent::<S, N>::zero())) == bits(&[S::zero(); D])
        && bits(&round) == bits(&da)
        && bits(&[a.dot(&b)]) == bits(&[dot])
        && bits(&[a.dot_acc(&b, k)]) == bits(&[dot_k])
}

/// The converters are a permutation of six scalars: what goes in comes out, bit for bit, and the
/// rotation-first dense order is `[ω; v]`.
fn twist_round_trips<S: Lanes>(a: [f64; 6]) -> bool {
    let s: [S; 6] = array::from_fn(|i| S::sample(a[i], i % 3));
    let (v, w) = ([s[0], s[1], s[2]], [s[3], s[4], s[5]]);
    let t = Twist::from_translation_first(s);
    let r = Twist {
        phi: Vector(v),
        rho: [Vector(w)],
    };
    let (fw, dense) = (r.to_translation_first(), dense_of::<S, 1, 6>(&t));
    bits(&t.to_translation_first()) == bits(&s)
        && bits(&t.v().0) == bits(&v)
        && bits(&t.omega().0) == bits(&w)
        && bits(&dense) == bits(&[s[3], s[4], s[5], s[0], s[1], s[2]])
        && bits(&fw) == bits(&[w[0], w[1], w[2], v[0], v[1], v[2]])
        && bits(&Twist::from_translation_first(fw).to_translation_first()) == bits(&fw)
}

macro_rules! props {
    ($m:ident, $S:ty) => {
        mod $m {
            use super::*;
            type S = $S;
            proptest! {
                #[test]
                fn tangent_ops_match_the_dense_arrays(v in proptest::collection::vec(-1.0_f64..1.0, 25)) {
                    prop_assert!(ops_match_dense::<S, 1, 6>(&v));
                    prop_assert!(ops_match_dense::<S, 2, 9>(&v));
                    prop_assert!(ops_match_dense::<S, 3, 12>(&v));
                }
                #[test]
                fn twist_translation_first_round_trip(a in proptest::array::uniform6(proptest::num::f64::ANY)) {
                    prop_assert!(twist_round_trips::<S>(a));
                }
            }
        }
    };
}

props!(as_f64, f64);
props!(as_f32, f32);
props!(as_dual, Dual<f64, 3>);

fn fb<const D: usize>(a: [f64; D]) -> [u64; D] {
    a.map(f64::to_bits)
}

#[test]
fn the_dense_order_is_phi_then_rho() {
    let v = |a: f64| Vector([a, a + 1.0, a + 2.0]);
    let t = SEn3Tangent::<f64, 3> {
        phi: v(1.0),
        rho: [v(4.0), v(7.0), v(10.0)],
    };
    let mut out = [0.0; 12];
    t.write_dense(&mut out);
    let ramp: [f64; 12] = array::from_fn(|i| 1.0 + i as f64);
    assert_eq!(fb(out), fb(ramp));
    let back = SEn3Tangent::<f64, 3>::read_dense(&ramp);
    assert_eq!(fb(back.phi.0), fb(t.phi.0));
    assert_eq!(back.rho.map(|r| fb(r.0)), t.rho.map(|r| fb(r.0)));
    assert_eq!(
        (SEn3Tangent::<f64, 1>::DOF, SEn3Tangent::<f64, 3>::DOF),
        (6, 12)
    );
}

/// Barfoot's `[ρ; φ]`, Solà's, Sophus's and manif's `[υ; ω]` and locus-tag's `[v, ω]` against the
/// rotation-first `[ω; v]`, by hand.
#[test]
fn twist_translation_first_by_hand() {
    let t = Twist::from_translation_first([1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0]);
    assert_eq!(
        (fb(t.v().0), fb(t.omega().0)),
        (fb([1.0, 2.0, 3.0]), fb([4.0, 5.0, 6.0]))
    );
    let mut out = [0.0; 6];
    t.write_dense(&mut out);
    assert_eq!(fb(out), fb([4.0, 5.0, 6.0, 1.0, 2.0, 3.0]));
    assert_eq!(
        fb(t.to_translation_first()),
        fb([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    );
    // A pure rotation about `z` is `[0, 0, 0, 0, 0, w]` translation-first and `[0, 0, w, 0, 0, 0]`
    // rotation-first.
    let z = Twist::from_translation_first([0.0_f64, 0.0, 0.0, 0.0, 0.0, 7.0]);
    assert_eq!(
        (fb(z.omega().0), fb(z.v().0)),
        (fb([0.0, 0.0, 7.0]), fb([0.0; 3]))
    );
}

/// A domain is a `debug_assert!` (D11): the debug profile that `just test` runs panics ...
#[cfg(debug_assertions)]
mod out_of_domain {
    use super::*;

    #[test]
    #[should_panic(expected = "Tangent::write_dense")]
    fn write_dense_of_the_wrong_length() {
        SEn3Tangent::<f64, 2>::zero().write_dense(&mut [0.0; 10]);
    }

    #[test]
    #[should_panic(expected = "Tangent::read_dense")]
    fn read_dense_of_the_wrong_length() {
        let _ = SEn3Tangent::<f64, 2>::read_dense(&[0.0; 8]);
    }
}

/// ... and a release build never panics: `write_dense` writes `min(len, DOF)` entries and
/// `read_dense` reads a missing entry as NaN, never as `+0`, which would be a valid component and
/// so a plausible tangent (`Tangent::read_dense`).
#[cfg(not(debug_assertions))]
#[test]
fn out_of_domain_does_not_panic_in_release() {
    let t = SEn3Tangent::<f64, 1>::read_dense(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    let (mut short, mut long) = ([9.0; 4], [9.0; 8]);
    t.write_dense(&mut short);
    t.write_dense(&mut long);
    assert_eq!(fb(short), fb([1.0, 2.0, 3.0, 4.0]));
    assert_eq!(fb(long), fb([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 9.0, 9.0]));
    let back = SEn3Tangent::<f64, 1>::read_dense(&short);
    assert_eq!(fb(back.phi.0), fb([1.0, 2.0, 3.0]));
    let rho = back.rho[0].0;
    assert_eq!(fb([rho[0]]), fb([4.0]));
    assert!(rho[1].is_nan() && rho[2].is_nan());
}

/// The group `SEn3`: the generic laws at `N = 1` and `N = 2`, §14's `jr_inv` twin at `N = 1, 2, 3`,
/// and the hand cases of §5 the laws do not reach.
///
/// The bounds are twice the worst error of **60 000 draws** of `laws::sample`'s distribution per
/// law per scalar, rounded up, as [`Bounds`] asks; the worst of the six runs (`N = 1, 2` at `f64`,
/// `f32` and `Dual`) is the one recorded, and the figures are in the comment on [`F64`]. One
/// measurement at 20 000 draws first read `adjoint` 7.01 where 60 000 reads 8.28 and `ad` 7.06
/// where it reads 7.86, so the sample count is part of the claim, as it is for SO(3).
mod group {
    // Not `use super::*`: the file's glob of `proptest::prelude` would shadow the one
    // `laws_for!` imports into each case module, and an unused import is a denied warning.
    use crate::laws::{laws_for, tangent as law_tangent, Bounds, Sample};
    use crate::{reference, Jac, Left, LieGroup, Right, SEn3, SEn3Jac, SE23, SE3};
    use crate::{SEn3Tangent, Tangent, Twist};
    use core::array;
    use helicoid_linalg::{Dual, Vector};
    use helicoid_linalg::{Mat3, Point, StridedMut};

    // Measured maxima over the six runs, in `u`: axioms 6.085, exp_log 4.171, adjoint 8.279,
    // jl_ad_jr 4.490, plus_minus 7.561, rows 2.182, ad 7.855, sides 1.118, tangent_order 0 at
    // `f64` and `Dual` / 2.694 at `f32`, jac_order 3.882, sandwich 3.644.
    //
    // `rows` is **not** `0` here, where SO(3) records exactly `0`: that group's `Ad_Y⁻¹` is a
    // transpose, which reproduces the dense chain bit for bit, while this one's is `Ad` of the
    // group inverse — a conjugate and `N` sandwiches, whose rounding is its own. An exactness
    // claim would be false, so a measured bound stands in its place.
    const F64: Bounds = Bounds {
        axioms: 13.0,
        exp_log: 9.0,
        adjoint: 17.0,
        jl_ad_jr: 9.0,
        plus_minus: 16.0,
        rows: 5.0,
        ad: 16.0,
        sides: 3.0,
        tangent_order: 0.0,
        jac_order: 8.0,
        sandwich: 8.0,
    };
    // `tangent_order` is one rounding at `f32` where it is exact at `f64`, as for SO(3); every
    // other law agrees within 25% across the precisions, so one set serves them.
    const F32: Bounds = Bounds {
        tangent_order: 6.0,
        ..F64
    };

    type G1<S> = SE3<S>;
    type G2<S> = SE23<S>;

    /// `J_r(τ)` at the sample: invertible for `|φ| < 2π`, which `laws::sample` is, and neither
    /// symmetric nor block-diagonal, so `apply_transpose` and the `col` blocks are both exercised.
    fn jac1<S: Sample>(v: &[f64; 6]) -> SEn3Jac<S, 1> {
        <G1<S> as LieGroup<S>>::jr(&law_tangent::<S, G1<S>, 6>(v))
    }

    fn jac2<S: Sample>(v: &[f64; 9]) -> SEn3Jac<S, 2> {
        <G2<S> as LieGroup<S>>::jr(&law_tangent::<S, G2<S>, 9>(v))
    }

    laws_for!(se3, G1, jac1, 6, F64, F32);
    laws_for!(se23, G2, jac2, 9, F64, F32);

    /// `jr_inv` against §14's twin, the dense Gauss-Jordan inverse of `jr` (`reference::
    /// sen3jac_inverse`), which is what makes §5.4's closed-form `A⁻¹` an answer and not a claim.
    ///
    /// `7` is twice the worst of 60 000 draws per scalar (15 000 at `N = 3`): 3.437 at `N = 2`
    /// `f64`, 3.243/3.132/3.166 at `N = 1` for `f64`/`f32`/`Dual`, 3.258 at `N = 3`.
    fn jr_inv_vs_reference<S: Sample, const N: usize, const D: usize>(v: &[f64; D]) -> f64 {
        let t = law_tangent::<S, SEn3<S, N>, D>(v);
        let (j, fast) = (
            <SEn3<S, N> as LieGroup<S>>::jr(&t),
            <SEn3<S, N> as LieGroup<S>>::jr_inv(&t),
        );
        let (mut got, mut want) = ([[S::zero(); D]; D], [[S::zero(); D]; D]);
        Jac::<S, SEn3Tangent<S, N>>::write_dense(
            &fast,
            &mut StridedMut::col_major(got.as_flattened_mut(), D, D),
        );
        reference::sen3jac_inverse::<S, N, D>(
            &j,
            &mut StridedMut::col_major(want.as_flattened_mut(), D, D),
        );
        let f = |m: &[[S; D]; D]| -> [[f64; D]; D] { m.map(|r| r.map(|x| x.value_f64())) };
        let (gv, wv) = (f(&got), f(&want));
        crate::laws::e::<S>(gv.as_flattened(), wv.as_flattened())
    }

    proptest::proptest! {
        #[test]
        fn sen3_jr_inv_matches_reference_n1(v in crate::laws::sample::<6>()) {
            crate::laws::within(jr_inv_vs_reference::<f64, 1, 6>(&v), 7.0)?;
            crate::laws::within(jr_inv_vs_reference::<f32, 1, 6>(&v), 7.0)?;
            crate::laws::within(jr_inv_vs_reference::<Dual<f64, 6>, 1, 6>(&v), 7.0)?;
        }
        #[test]
        fn sen3_jr_inv_matches_reference_n2(v in crate::laws::sample::<9>()) {
            crate::laws::within(jr_inv_vs_reference::<f64, 2, 9>(&v), 7.0)?;
            crate::laws::within(jr_inv_vs_reference::<f32, 2, 9>(&v), 7.0)?;
        }
        #[test]
        fn sen3_jr_inv_matches_reference_n3(v in crate::laws::sample::<12>()) {
            crate::laws::within(jr_inv_vs_reference::<f64, 3, 12>(&v), 7.0)?;
        }
    }

    fn twist(v: [f64; 6]) -> Twist<f64> {
        Twist::read_dense(&v)
    }

    fn close(a: &[f64], b: &[f64], tol: f64) {
        let err = a
            .iter()
            .zip(b)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max);
        assert!(
            a.len() == b.len() && err <= tol,
            "{a:?} vs {b:?} (tol {tol})"
        );
    }

    /// `Exp` of a pure translation is the translation itself, to the bit: `W = 0`, so `J_l(0)` is
    /// `I` and each column passes through two exactly-zero cross products. `Log` inverts it.
    #[test]
    fn exp_of_a_pure_translation_is_exact() {
        let rho = [0.25, -1.5, 3.0];
        let x = SE3::<f64>::exp(&twist([0.0, 0.0, 0.0, rho[0], rho[1], rho[2]]));
        let (r, t) = (x.rotation().quat(), x.translation());
        assert_eq!(
            [r.w, r.x, r.y, r.z].map(f64::to_bits),
            [1.0, 0.0, 0.0, 0.0].map(f64::to_bits)
        );
        assert_eq!(t.0.map(f64::to_bits), rho.map(f64::to_bits));
        let mut back = [0.0; 6];
        x.log().write_dense(&mut back);
        assert_eq!(
            back.map(f64::to_bits),
            [0.0, 0.0, 0.0, 0.25, -1.5, 3.0].map(f64::to_bits)
        );
    }

    /// The action is `R p + t`, and composition acts as the composition of the actions.
    #[test]
    fn the_action_composes() {
        let (a, b) = (
            SE3::<f64>::exp(&twist([0.3, -0.7, 1.1, 0.5, -2.0, 0.25])),
            SE3::<f64>::exp(&twist([-1.2, 0.4, 0.9, 3.0, 1.0, -0.5])),
        );
        let p = Point([1.5, -0.5, 2.0]);
        let direct = (a * b) * p;
        let staged = a * (b * p);
        close(&direct.0, &staged.0, 1e-14);
        // `R p + t` written out, against the operator. Not bit equality: the operator is the
        // quaternion sandwich of §3.3 and this is the matrix form, which rounds `R`'s nine entries
        // first -- the difference `act_many_is_the_matrix_action` pins from the other side.
        let m = a.rotation().to_matrix();
        let want = m * Vector(p.0) + a.translation();
        close(&(a * p).0, &want.0, 1e-14);
    }

    /// `act_many` is the matrix action per point, where `act` is the quaternion sandwich: the two
    /// differ in the bits, which is why §14 keeps SO(3)'s `act_many` row.
    #[test]
    fn act_many_is_the_matrix_action() {
        let x = SE3::<f64>::exp(&twist([0.3, -0.7, 1.1, 0.5, -2.0, 0.25]));
        let pts = [Point([1.5, -0.5, 2.0]), Point([0.0, 1.0, -3.0])];
        let mut many = pts;
        x.act_many(&mut many);
        let m = x.rotation().to_matrix();
        for (got, p) in many.iter().zip(&pts) {
            let want = m * Vector(p.0) + x.translation();
            assert_eq!(got.0.map(f64::to_bits), want.0.map(f64::to_bits));
        }
    }

    /// §2.4's action Jacobians against `Dual` differentiation of the action itself, both sides.
    ///
    /// The only check these rows have: `PHASE3.md` §8's laws are owed, and no corpus id covers
    /// them. The tolerance is a hand case's, not a recorded bound.
    #[test]
    fn act_jacobians_differentiate_the_action() {
        type D6 = Dual<f64, 6>;
        let base = [0.3, -0.7, 1.1, 0.5, -2.0, 0.25];
        let p = [1.5, -0.5, 2.0];
        let x64 = SE3::<f64>::exp(&twist(base));
        let xd = SE3::<D6>::exp(&SEn3Tangent::read_dense(&base.map(D6::constant)));
        let delta =
            SEn3Tangent::<D6, 1>::read_dense(&array::from_fn::<D6, 6, _>(|i| D6::variable(0.0, i)));
        for right in [true, false] {
            let moved = match right {
                true => xd.rplus(&delta),
                false => xd.lplus(&delta),
            };
            let got = moved * Point(p.map(D6::constant));
            let (j, dp): (_, Mat3<f64>) = match right {
                true => x64.act_jacobians::<Right>(Point(p)),
                false => x64.act_jacobians::<Left>(Point(p)),
            };
            for r in 0..3 {
                for c in 0..6 {
                    let (lane, want) = (got.0[r].d[c], j.get(r, c));
                    assert!(
                        (lane - want).abs() <= 1e-13,
                        "{right} ({r}, {c}): {lane} vs {want}"
                    );
                }
            }
            // `∂(X p)/∂p = R` on both sides.
            let m = x64.rotation().to_matrix();
            for r in 0..3 {
                close(&dp.row(r).0, &m.row(r).0, 0.0);
            }
        }
    }

    /// `from_parts` and `parts` are each other's inverse, exactly, and the `N = 2` accessors name
    /// the columns `0012`'s order gives them.
    #[test]
    fn parts_round_trip_and_the_se23_accessors_name_their_columns() {
        let x = SE23::<f64>::exp(&SEn3Tangent::read_dense(&[
            0.3, -0.7, 1.1, 0.5, -2.0, 0.25, -1.0, 4.0, 0.125,
        ]));
        let (r, cols) = x.parts();
        let again = SE23::from_parts(r, cols);
        let (r2, cols2) = again.parts();
        assert_eq!(
            [r.quat().w, r.quat().x].map(f64::to_bits),
            [r2.quat().w, r2.quat().x].map(f64::to_bits)
        );
        for i in 0..2 {
            assert_eq!(cols[i].0.map(f64::to_bits), cols2[i].0.map(f64::to_bits));
        }
        assert_eq!(
            x.velocity().0.map(f64::to_bits),
            cols[0].0.map(f64::to_bits)
        );
        assert_eq!(
            x.position().0.map(f64::to_bits),
            cols[1].0.map(f64::to_bits)
        );
    }

    /// `from_quat_translation` is `from_rt` of the unchecked constructor, and `rotation` reads it
    /// back: the boundary where a caller's claim of unitness is taken on trust (`0027`).
    #[test]
    fn from_quat_translation_keeps_the_quaternion_it_was_given() {
        let q = crate::Quat::from_wxyz_unchecked(1.0, 0.0, 0.0, 0.0);
        let x = SE3::from_quat_translation(q, Vector([1.0, 2.0, 3.0]));
        let got = x.rotation().quat();
        assert_eq!(
            [got.w, got.x, got.y, got.z].map(f64::to_bits),
            [q.w, q.x, q.y, q.z].map(f64::to_bits)
        );
        assert_eq!(
            x.translation().0.map(f64::to_bits),
            [1.0, 2.0, 3.0].map(f64::to_bits)
        );
    }
}
