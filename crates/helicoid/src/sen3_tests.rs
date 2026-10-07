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

/// Barfoot's `[ρ; φ]` and Solà's, Sophus's, and manif's `[υ; ω]` against
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
    // jl_ad_jr 4.490, rows 0 (exact), ad 7.855, sides 1.118, tangent_order 0 at `f64` and `Dual`
    // / 2.694 at `f32`, jac_order 3.882, sandwich 3.644.
    //
    // `plus_minus` is 21, from `measure_plus_minus`: 7.110 / 7.898 / 7.110 / 10.346 / 7.419 over
    // `se3` at `f64`/`f32`/`Dual` and `se23` at `f64`/`f32`. The bound was 16, recorded from
    // 7.561 on the *proptest* stream, and the convention is twice the worst of a 10^6-draw sample
    // of the distribution the bar draws -- so 16 was never twice a worst, it was twice the
    // luckier of two samplers: at `0048`'s parent commit this same measurement reads 9.666, which
    // already wants 20. `0048` then moved the figure by 1.07x (9.666 -> 10.346, `se23` `f64`,
    // against `se3` improving at every scalar), so the two causes are separated and neither is
    // hidden behind the other.
    //
    // `rows` is `0`: an exactness claim, and it keeps no headroom. Every §2.3 row is reproduced
    // bit for bit, which is what taking `Ad_Exp(τ)⁻¹` as `Ad_Exp(−τ)` buys — through `Ad` of the
    // group inverse the same law read 2.182 `u`, because `inverse` reaches its columns as
    // `−Rᵗ(J_l ρ)` where `Ad_Exp(−τ)` reaches them as `−J_r ρ`.
    const F64: Bounds = Bounds {
        axioms: 13.0,
        exp_log: 9.0,
        adjoint: 17.0,
        jl_ad_jr: 9.0,
        plus_minus: 21.0,
        rows: 0.0,
        ad: 16.0,
        sides: 3.0,
        tangent_order: 0.0,
        jac_order: 8.0,
        sandwich: 8.0,
        // `PHASE3.md` §8's second check, measured 6.946 at `N = 1` and 5.274 at `N = 2` over 10 000 draws each.
        dual_rows: 14.0,
        // `PHASE4.md` §1 and §3's seven legs, in `GEODESIC_LEGS`'s order, each twice the worst of 10^6 draws of `laws::Rng::shaped` -- `laws::sample`'s own distribution, which is what the proptest draws -- rounded up: the worst of `N = 1` and `N = 2`: `t=0` 1.118, `t=1` 6.578, symmetry 14.694, velocity 7.509, left 10.731, right 8.316, `twin` 8.579.
        //
        // Re-recorded for `0054`, whose screw twin makes `N = 1`'s `twin` leg compare two genuinely different expressions (1.118 -> 8.579 at `f32`, 8.051 at `f64`, against SO(3)'s 7.213 for the same reason). The world-frame translation tightened `N = 1`'s own legs -- symmetry 16.111 -> 11.415, right 9.360 -> 8.211 -- so `N = 2`'s figures now bind every leg but `twin` and `velocity`. Accuracy against the 110-digit corpus is `0054`'s table, not these.
        geodesic: [3.0, 14.0, 30.0, 16.0, 22.0, 17.0, 18.0],
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
    ///
    /// # Validity
    ///
    /// The twin's own forward error is `O(κ u)` in the conditioning of the *dense* matrix, which
    /// is why `sen3jac_inverse_matches_reference` rejects draws with `κ u > 1e-3`. A flat bound is
    /// sound here only because `laws::sample` caps every entry at `|m| 2^0`, so `θ <= √3` and
    /// `J_r` stays far from the `θ = 2π` singularity that `jr_inv`'s *Domain* names. Widening the
    /// sample toward `2π` would need the `κ` filter, and would be measuring the twin.
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
            crate::laws::within(jr_inv_vs_reference::<Dual<f64, 9>, 2, 9>(&v), 7.0)?;
        }
        // `N = 3` at every scalar: it is the only width where `col` holds more than two blocks, so
        // the only one where a per-block index swap or a per-lane `Dual` error cannot hide behind
        // symmetry (`0042` says so itself), and the generic laws are not instantiated there.
        #[test]
        fn sen3_jr_inv_matches_reference_n3(v in crate::laws::sample::<12>()) {
            crate::laws::within(jr_inv_vs_reference::<f64, 3, 12>(&v), 7.0)?;
            crate::laws::within(jr_inv_vs_reference::<f32, 3, 12>(&v), 7.0)?;
            crate::laws::within(jr_inv_vs_reference::<Dual<f64, 12>, 3, 12>(&v), 7.0)?;
        }
    }

    /// Both relative-transform spellings against the compositions they are *not* bit-identical
    /// to, `a * b.inverse()` and `b.inverse() * a` (`NUMERICS.md` §14, `0048` decision 4): one
    /// rotation per column against two.
    ///
    /// Scored on the columns against the **input** scale, `max(‖x_a‖, ‖x_b‖, 1)`, and at three
    /// translation scales per draw — `1`, `10³`, `10⁶`, each applied to both elements. That is
    /// the only denominator in which the figure is a property of the routine: both forms carry an
    /// absolute error of order `u max‖x‖`, so `laws::e`'s `max(‖want‖, 1)` divides by the
    /// *difference* and reports the cancellation of two nearby frames instead of the rounding.
    /// On one set of draws it reads 6.80 u at scale `1` and 106.32 u at `10³`, where against the
    /// input scale the same draws stay 6.80 u and 7.94 u — flat, which is the claim.
    ///
    /// The quaternion half is not scored because it is bit-identical by construction — both
    /// spellings reach it through the one Hamilton product of a conjugate — and
    /// `the_relative_spellings_share_the_quaternion_and_part_on_the_columns` asserts exactly that,
    /// so a figure folding the two together would be reporting a zero.
    fn rel_vs_reference<S: Sample, const N: usize, const D: usize>(
        v: &[f64],
        inv_mul: bool,
    ) -> f64 {
        let half = |o: usize| -> SEn3<S, N> {
            SEn3::exp(&law_tangent::<S, SEn3<S, N>, D>(&array::from_fn(|i| {
                v[o + i]
            })))
        };
        let (a, b) = (half(0), half(D));
        let cols = |g: &SEn3<S, N>| g.parts().1.map(|c| c.0.map(|s| s.value_f64()));
        let mut out = 0.0;
        for e in [1.0, 1e3, 1e6] {
            let k = S::lit(e);
            let lift = |g: &SEn3<S, N>| {
                let (r, c) = g.parts();
                SEn3::from_parts(r, c.map(|v| v.scale(k)))
            };
            let (x, y) = (lift(&a), lift(&b));
            let (got, want) = if inv_mul {
                (x.inv_mul(&y), y.inverse() * x)
            } else {
                (x.mul_inv(&y), x * y.inverse())
            };
            let (gc, wc, xc, yc) = (cols(&got), cols(&want), cols(&x), cols(&y));
            // Per column: each is an independent translation, so each carries its own scale, and
            // folding `N` of them into one vector would score the smallest against the largest.
            for i in 0..N {
                let scale =
                    crate::laws::worst(crate::laws::norm(&xc[i]), crate::laws::norm(&yc[i]));
                out = crate::laws::worst(out, crate::laws::e_at::<S>(&gc[i], &wc[i], scale));
            }
        }
        out
    }

    proptest::proptest! {
        /// `18` is twice the worst of 10^6 `laws::shaped` draws at three scales each, rounded
        /// up (`rel_twin_measure`): `mul_inv` 8.119/8.119/8.415 u at `N = 1` for
        /// `f64`/`Dual`/`f32` and 8.893/8.994 at `N = 2`; `inv_mul` 7.802/7.802/7.512 and
        /// 8.088/7.951. One bound serves both rows and both widths: the spread is 18% and a
        /// tighter one per row would be recording the sampler, not the routines.
        #[test]
        fn sen3_mul_inv_matches_reference_n1(v in crate::laws::sample::<12>()) {
            crate::laws::within(rel_vs_reference::<f64, 1, 6>(&v, false), 18.0)?;
            crate::laws::within(rel_vs_reference::<f32, 1, 6>(&v, false), 18.0)?;
            crate::laws::within(rel_vs_reference::<Dual<f64, 6>, 1, 6>(&v, false), 18.0)?;
        }
        #[test]
        fn sen3_inv_mul_matches_reference_n1(v in crate::laws::sample::<12>()) {
            crate::laws::within(rel_vs_reference::<f64, 1, 6>(&v, true), 18.0)?;
            crate::laws::within(rel_vs_reference::<f32, 1, 6>(&v, true), 18.0)?;
            crate::laws::within(rel_vs_reference::<Dual<f64, 6>, 1, 6>(&v, true), 18.0)?;
        }
        /// `N = 2`: the column loop, which `N = 1` cannot distinguish from a single `act`.
        #[test]
        fn sen3_mul_inv_matches_reference_n2(v in crate::laws::sample::<18>()) {
            crate::laws::within(rel_vs_reference::<f64, 2, 9>(&v, false), 18.0)?;
            crate::laws::within(rel_vs_reference::<f32, 2, 9>(&v, false), 18.0)?;
        }
        #[test]
        fn sen3_inv_mul_matches_reference_n2(v in crate::laws::sample::<18>()) {
            crate::laws::within(rel_vs_reference::<f64, 2, 9>(&v, true), 18.0)?;
            crate::laws::within(rel_vs_reference::<f32, 2, 9>(&v, true), 18.0)?;
        }
    }

    /// `plus_minus`'s figure, which `0048` moved by making `rminus`/`lminus` the relative
    /// spellings. Its own `Rng` stream, so every other recorded figure stays reproducible.
    /// `cargo nextest run -p helicoid -- measure_plus_minus --ignored`, release.
    #[test]
    #[allow(clippy::print_stdout)]
    #[ignore = "measurement: prints the figure `Bounds::plus_minus` is recorded from"]
    fn measure_plus_minus() {
        fn run<S: Sample, const N: usize, const D: usize>(name: &str) {
            let mut rng = crate::laws::Rng(0x706C_7573_5F6D_696E);
            let mut w = 0.0;
            for _ in 0..1_000_000 {
                let (a, b, c) = (rng.shaped::<D>(), rng.shaped::<D>(), rng.shaped::<D>());
                // `(g(a), g(b), t(c))`, the order the generated proptest passes, so the figure
                // is of the same law on the same arguments and only the stream differs.
                let at = |v: &[f64; D]| law_tangent::<S, SEn3<S, N>, D>(v);
                let (x, y) = (SEn3::<S, N>::exp(&at(&a)), SEn3::<S, N>::exp(&at(&b)));
                w = crate::laws::worst(
                    w,
                    crate::laws::plus_minus::<S, SEn3<S, N>, D>(&x, &y, &at(&c)),
                );
            }
            std::println!("{name} plus_minus {w:.3}");
        }
        run::<f64, 1, 6>("se3 f64");
        run::<f32, 1, 6>("se3 f32");
        run::<Dual<f64, 6>, 1, 6>("se3 dual");
        run::<f64, 2, 9>("se23 f64");
        run::<f32, 2, 9>("se23 f32");
    }

    /// The figures the twin bound and `mul_inv`'s rustdoc quote, and the two denominators side by
    /// side (`0048` decision 4). `cargo nextest run -p helicoid -- rel_twin_measure --ignored
    /// --no-capture`, release.
    #[test]
    #[allow(clippy::print_stdout)]
    #[ignore = "a measurement, not a law; its result is the bound above"]
    fn rel_twin_measure() {
        const DRAWS: usize = 1_000_000;
        fn run<S: Sample, const N: usize, const D: usize>(name: &str, inv_mul: bool) {
            let mut rng = crate::laws::Rng(0x7265_6C5F_7477_696E);
            let mut worst = 0.0;
            for _ in 0..DRAWS {
                let v: [f64; 18] = rng.shaped();
                worst = crate::laws::worst(worst, rel_vs_reference::<S, N, D>(&v, inv_mul));
            }
            std::println!("{name}: {worst} u of the input scale");
        }
        run::<f64, 1, 6>("mul_inv n1 f64", false);
        run::<f32, 1, 6>("mul_inv n1 f32", false);
        run::<Dual<f64, 6>, 1, 6>("mul_inv n1 dual", false);
        run::<f64, 2, 9>("mul_inv n2 f64", false);
        run::<f32, 2, 9>("mul_inv n2 f32", false);
        run::<f64, 1, 6>("inv_mul n1 f64", true);
        run::<f32, 1, 6>("inv_mul n1 f32", true);
        run::<Dual<f64, 6>, 1, 6>("inv_mul n1 dual", true);
        run::<f64, 2, 9>("inv_mul n2 f64", true);
        run::<f32, 2, 9>("inv_mul n2 f32", true);
        // The two denominators on the same draws: `laws::e` divides by the result, `e_at` by the
        // inputs. The first is what a naive twin bound would have recorded.
        let mut rng = crate::laws::Rng(0x7265_6C5F_7477_696E);
        for e in [1.0_f64, 1e3, 1e6] {
            let (mut res, mut inp) = (0.0, 0.0);
            for _ in 0..DRAWS {
                let v: [f64; 12] = rng.shaped();
                let half = |o: usize| -> SE3<f64> {
                    let g = SE3::exp(&law_tangent::<f64, SE3<f64>, 6>(&array::from_fn(|i| {
                        v[o + i]
                    })));
                    let (r, c) = g.parts();
                    SE3::from_parts(r, c.map(|u| u.scale(e)))
                };
                let (x, y) = (half(0), half(6));
                let (got, want) = (x.mul_inv(&y), x * y.inverse());
                let c = |g: &SE3<f64>| g.parts().1[0].0;
                let scale =
                    crate::laws::worst(crate::laws::norm(&c(&x)), crate::laws::norm(&c(&y)));
                res = crate::laws::worst(res, crate::laws::e::<f64>(&c(&got), &c(&want)));
                inp = crate::laws::worst(inp, crate::laws::e_at::<f64>(&c(&got), &c(&want), scale));
            }
            std::println!("scale {e:e}: result-relative {res} u, input-relative {inp} u");
        }
    }

    /// What the two spellings share with the compositions they replace and where they part: the
    /// quaternion to the bit, the columns not (`0048` decision 4).
    ///
    /// Both halves are load-bearing and neither is a restatement of a body. The quaternion claim
    /// compares two *different* expressions — `q_a q_b*` against the `q` that `inverse` then `Mul`
    /// produce — and it is what lets the twin proptest score the columns alone. The column claim
    /// is `0048`'s case for the routines existing: if they agreed to the bit they would be
    /// spellings, so a sweep requires a disagreement and prints how often, since the *rate* is a
    /// distribution statistic and not an invariant — one pair is the assertion.
    #[test]
    fn the_relative_spellings_share_the_quaternion_and_part_on_the_columns() {
        let mut rng = crate::laws::Rng(0x6D75_6C5F_696E_7600);
        let (mut mul, mut inv) = (0usize, 0usize);
        const DRAWS: usize = 20_000;
        for _ in 0..DRAWS {
            let (a, b) = (rng.shaped::<9>(), rng.shaped::<9>());
            let x = SE23::<f64>::exp(&SEn3Tangent::read_dense(&a));
            let y = SE23::<f64>::exp(&SEn3Tangent::read_dense(&b));
            let quat = |g: &SE23<f64>| {
                let q = g.parts().0.quat();
                [q.w, q.x, q.y, q.z].map(f64::to_bits)
            };
            let cols = |g: &SE23<f64>| g.parts().1.map(|c| c.0.map(f64::to_bits));
            for (got, want, n) in [
                (x.mul_inv(&y), x * y.inverse(), &mut mul),
                (x.inv_mul(&y), y.inverse() * x, &mut inv),
            ] {
                assert_eq!(quat(&got), quat(&want), "the quaternion halves parted");
                if cols(&got) != cols(&want) {
                    *n += 1;
                }
            }
        }
        assert!(
            mul > 0 && inv > 0,
            "no column differed in {DRAWS} draws ({mul} / {inv}): these would be spellings, \
             not routines"
        );
    }

    /// `renormalize` repairs the drift a chain of compositions accumulates, on the rotation only
    /// (`0044` item 2).
    ///
    /// The drift is *composed*, not scaled on: `Mul` never renormalizes (`0027`), so a chain is
    /// where `η = ‖q‖² − 1` actually comes from and the only case that shows the method is for
    /// the domain it claims. The step itself is `Quat::renormalize`'s and
    /// `quat_tests::renormalize_converges_quadratically` is where its `−¾η² + ¼η³` is checked
    /// across scalars; here the claim is the delegation, pinned on an exactly representable `η`
    /// rather than against a second call of the same body, which could not fail.
    #[test]
    fn renormalize_repairs_a_composed_drift_and_leaves_the_columns_alone() {
        let seed: [f64; 9] = [0.3, -0.7, 1.1, 0.5, -2.0, 0.25, -1.0, 4.0, 0.125];
        let step = SE23::<f64>::exp(&SEn3Tangent::read_dense(&seed));
        let mut x = step;
        for _ in 0..200 {
            x = x * step;
        }
        let eta = |g: &SE23<f64>| g.parts().0.quat().norm_sq() - 1.0;
        let drifted = eta(&x);
        // The composed drift is real and inside the one-step domain, `|η| <= 2^-26.29` (§3.6).
        assert!(
            drifted != 0.0 && drifted.abs() <= f64::powf(2.0, -26.29),
            "200 compositions drifted by {drifted}, outside the one-step domain"
        );
        let (_, before) = x.parts();
        let mut repaired = x;
        repaired.renormalize();
        let (_, after) = repaired.parts();
        // `η' = −¾η² + ¼η³` to rounding, so one step lands at `u` from a drift this size.
        let predicted = -0.75 * drifted * drifted + 0.25 * f64::powi(drifted, 3);
        let off = (eta(&repaired) - predicted).abs() / crate::laws::unit::<f64>();
        assert!(off <= 4.0, "{off} u off the closed form");
        assert!(eta(&repaired).abs() <= 4.0 * crate::laws::unit::<f64>());
        for (a, b) in after.iter().zip(&before) {
            assert_eq!(a.0.map(f64::to_bits), b.0.map(f64::to_bits));
        }
        // The delegation, on the golden of `quat_tests::newton_step_hand_case`: `η = 0.5625`,
        // `k = 0.71875`, both exact in binary, so `1.25` must come back `0.8984375`.
        // The struct literal, not `from_wxyz_unchecked`: `η = 0.5625` is far outside that
        // constructor's `2^-40` debug assert, which is the whole point of the case.
        let far = crate::Quat {
            w: 1.25_f64,
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let mut hand = SE3::from_rt(
            crate::SO3::from_quat_unchecked(far),
            Vector([1.0, 2.0, 3.0]),
        );
        hand.renormalize();
        let q = hand.parts().0.quat();
        assert_eq!(
            [q.w, q.x, q.y, q.z].map(f64::to_bits),
            [0.898_437_5_f64, 0.0, 0.0, 0.0].map(f64::to_bits)
        );
    }

    /// The fused `rminus_jacobians` against the two inversions it replaces, **on the bits**.
    ///
    /// `laws::jacobian_rows` bounds this row at `0 u`, which is an *error* bound and so cannot see
    /// a signed zero: `laws::e` scores `+0` against `−0` as zero. This compares `to_bits`, over a
    /// random sweep and over the degenerate case the fusion is most likely to part company on — a
    /// pure translation, `φ = 0`, where every word of `Q` is a zero matrix and only its signs are
    /// left to disagree about.
    #[test]
    fn the_fused_inverses_are_the_separate_ones_to_the_bit() {
        let mut st = 0x2545_F491_4F6C_DD1D_u64;
        let mut next = || {
            st ^= st << 13;
            st ^= st >> 7;
            st ^= st << 17;
            (st >> 11) as f64 / (1u64 << 53) as f64 * 4.0 - 2.0
        };
        let dense = |j: &SEn3Jac<f64, 1>| -> [[f64; 6]; 6] {
            let mut out = [[0.0; 6]; 6];
            Jac::<f64, Twist<f64>>::write_dense(
                j,
                &mut StridedMut::col_major(out.as_flattened_mut(), 6, 6),
            );
            out
        };
        // Case 0 is a pair with the *same* rotation, whose `rminus` is a pure translation:
        // `y⁻¹x = (I, Rᵗ(c_x − c_y))`, so `φ = 0`, every word of `Q` is a zero matrix, and only
        // the signs of those zeros are left to disagree about.
        let r = SE3::<f64>::exp(&twist([0.3, -0.7, 1.1, 0.0, 0.0, 0.0])).rotation();
        let mut differing = 0usize;
        for case in 0..2001 {
            let (x, y) = match case {
                0 => (
                    SE3::from_rt(r, helicoid_linalg::Vector([0.5, -2.0, 0.25])),
                    SE3::from_rt(r, helicoid_linalg::Vector([-1.0, 4.0, 0.125])),
                ),
                _ => {
                    let a = twist(array::from_fn(|_| next()));
                    let b = twist(array::from_fn(|_| next()));
                    (SE3::exp(&a), SE3::exp(&b))
                }
            };
            let tau = x.rminus(&y);
            let want = (
                SE3::<f64>::jr_inv(&tau),
                Jac::<f64, Twist<f64>>::neg(&SE3::<f64>::jl_inv(&tau)),
            );
            let got = x.rminus_jacobians(&y);
            for (g, w) in [(&got.0, &want.0), (&got.1, &want.1)] {
                let (dg, dw) = (dense(g), dense(w));
                let same = dg
                    .as_flattened()
                    .iter()
                    .zip(dw.as_flattened())
                    .all(|(a, b)| a.to_bits() == b.to_bits());
                if !same {
                    differing += 1;
                    for (a, b) in dg.as_flattened().iter().zip(dw.as_flattened()) {
                        assert!(a.to_bits() == b.to_bits(), "case {case}: {a} vs {b}");
                    }
                }
            }
        }
        // Recorded, not asserted at zero: the random sweep agrees to the bit, and the pure
        // translation is where a zero's sign is all that is left of `Q`.
        // An exactness claim, and it keeps no headroom: every one of the 4002 blocks matches.
        // It is `0` *because* each side builds its own `X` with `hat`: negating the other's gave
        // the `½ρ^` diagonal `−0` where the unfused path has the `+0` `hat` writes, which this
        // test sees and `laws::jacobian_rows`'s error bound cannot.
        assert_eq!(differing, 0, "of 4002 blocks");
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

    /// `act_many` against its two references: the per-point action, which is what §14's row for
    /// `SO3::act_many` names and the only one that can catch a swapped order, a skipped point or a
    /// translation applied before the rotation; and `SO3::act_many` plus the translation, which is
    /// the same matrix path through a different function and so holds to the bit.
    ///
    /// `NUMERICS.md` §14 has no row for `SEn3::act_many` and the corpus no id, so D6 is owed one;
    /// re-deriving this function's own loop body and asserting bit equality would have proved
    /// nothing, which is what a first version of this test did.
    #[test]
    fn act_many_is_the_matrix_action_and_agrees_with_the_per_point_one() {
        let x = SE3::<f64>::exp(&twist([0.3, -0.7, 1.1, 0.5, -2.0, 0.25]));
        let pts = [
            Point([1.5, -0.5, 2.0]),
            Point([0.0, 1.0, -3.0]),
            Point([-4.25, 0.125, 0.0]),
        ];
        let mut many = pts;
        x.act_many(&mut many);
        // The per-point action rounds a sandwich per point where this rounds `R(q)` once, so the
        // two differ in the bits: `3 u` is §14's tolerance for the SO(3) row this follows.
        let u = 3.0 * f64::EPSILON;
        for (got, p) in many.iter().zip(&pts) {
            let want = x * *p;
            let scale = want.0.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
            close(&got.0, &want.0, u * scale);
        }
        // And bit for bit against the rotation's own batch path plus the translation.
        let mut rotated = pts.map(|p| Vector(p.0));
        x.rotation().act_many(&mut rotated);
        for (got, r) in many.iter().zip(&rotated) {
            let want = *r + x.translation();
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

    /// `γ(x₀, x₁, 0)` returns `x₀` **bit for bit** at both widths, which `laws::geodesic`'s `t=0`
    /// leg cannot say: `gerr` compares two elements through `Log(x₀⁻¹ x₀)` and the quaternion's
    /// `q* q` leaves about one `u` in the vector part whatever the curve did.
    ///
    /// It holds by the arithmetic: `d.scale(0)` is `±0` per component, `Exp` of that is the
    /// identity quaternion with signed zeros in `x`, and composing adds those signed zeros to each
    /// component of `x₀` — exact, at `‖x₀‖ = 1e4` as at `0`.
    ///
    /// **With one exception, which the last fixture holds down:** a sum of signed zeros is `−0.0`
    /// only when every term is, so a `−0.0` *component of `x₀`* can come back `+0.0`. The value is
    /// unchanged and `gerr` reads its floor either way, but `to_bits` does not — and for the
    /// quaternion's `w` that sign is load-bearing, since `NUMERICS.md` §3.2 keeps `w = +0` and
    /// `Log`'s flip reads it. So the claim is bit identity *given no negative zero in the
    /// representation*, and `LieGroup::geodesic`'s rustdoc says it that way
    /// (`docs/maths/geodesics.md` GE.7(b), `0045` item 3).
    #[test]
    fn geodesic_at_zero_is_the_left_endpoint_bit_for_bit() {
        // `N = 1` is the width `tf_tree` migrates onto, so it is tested first and not only `N = 2`.
        let v1: [[f64; 6]; 4] = [
            [0.3, -0.7, 1.1, 0.5, -2.0, 0.25],
            [0.0, 0.0, 0.0, 1e4, -1e4, 1e4],
            [3.0, 0.1, -0.2, 1e-9, 1e-9, 1e-9],
            [1e-9, 0.0, -1e-9, 0.0, 0.0, 0.0],
        ];
        for v in v1 {
            // The step is a third of the fixture, negated: `d` is neither `0` nor `v`, which the
            // all-zero fixture of an earlier draft could not say.
            let d = SEn3Tangent::<f64, 1>::read_dense(&v.map(|c| -c / 3.0));
            let x0 = SE3::<f64>::exp(&SEn3Tangent::read_dense(&v));
            same_pose::<1>(&SE3::geodesic(&x0, &x0.rplus(&d), 0.0), &x0, &v);
        }
        let v2: [[f64; 9]; 3] = [
            [0.3, -0.7, 1.1, 0.5, -2.0, 0.25, -1.0, 4.0, 0.125],
            [0.0, 0.0, 0.0, 1e4, -1e4, 1e4, 1e-9, 0.0, -1e-9],
            [3.0, 0.1, -0.2, 1e-9, 1e-9, 1e-9, 1.0, -1.0, 0.5],
        ];
        for v in v2 {
            let d = SEn3Tangent::<f64, 2>::read_dense(&v.map(|c| -c / 3.0));
            let x0 = SE23::<f64>::exp(&SEn3Tangent::read_dense(&v));
            same_pose::<2>(&SE23::geodesic(&x0, &x0.rplus(&d), 0.0), &x0, &v);
        }

        // The exception: `w = −0.0` with a negative step, where `aw − ax·(−0) − …` sums to `+0.0`.
        let q = crate::Quat::<f64>::from_wxyz_normalized(-0.0, 0.6, 0.0, 0.8);
        let x0 = SE3::from_parts(
            crate::SO3::from_quat_unchecked(q),
            [Vector([1.0, 2.0, 3.0])],
        );
        let d = SEn3Tangent::<f64, 1>::read_dense(&[-0.3, -0.5, -0.7, -0.1, -0.2, -0.3]);
        let got = SE3::geodesic(&x0, &x0.rplus(&d), 0.0);
        let (a, b) = (got.parts().0.quat(), q);
        // A value comparison, not a float `==`: the two zeros compare equal.
        assert_eq!(
            a.w.partial_cmp(&b.w),
            Some(core::cmp::Ordering::Equal),
            "the value is unchanged"
        );
        assert_eq!(b.w.to_bits(), (-0.0_f64).to_bits());
        assert_eq!(
            a.w.to_bits(),
            0.0_f64.to_bits(),
            "a `−0.0` `w` is expected to come back `+0.0`; if this ever holds the sign instead, \
             `LieGroup::geodesic`'s rustdoc caveat can be dropped"
        );
    }

    /// Two poses equal component for component, on the bits.
    fn same_pose<const N: usize>(got: &SEn3<f64, N>, want: &SEn3<f64, N>, at: &[f64]) {
        let (rg, cg) = got.parts();
        let (rw, cw) = want.parts();
        let q =
            |r: crate::SO3<f64>| [r.quat().w, r.quat().x, r.quat().y, r.quat().z].map(f64::to_bits);
        assert_eq!(q(rg), q(rw), "t = 0 moved the rotation of {at:?}");
        for i in 0..N {
            assert_eq!(
                cg[i].0.map(f64::to_bits),
                cw[i].0.map(f64::to_bits),
                "t = 0 moved column {i} of {at:?}"
            );
        }
    }

    /// The screw twin's draws (`0054`): `(X₀, X₁ = X₀ Exp(d), t)` with the relative angle from one
    /// of three regimes -- `0`: log-uniform `1e-9 ..= 1e-3`, the `geo:consecutive` range; `1`:
    /// uniform `1e-3 ..= 3`; `2`: `π − 10^-k`, `k` uniform in `1 ..= 6` -- `‖x₀‖` log-uniform to
    /// `1e4`, and `t` cycling through the corpus's endpoints and extrapolations before a uniform
    /// draw. Returns the pose pair and `t`; the caller reads errors against the pair's scale.
    fn screw_draw<S: Sample>(
        rng: &mut crate::laws::Rng,
        regime: u64,
        i: usize,
    ) -> (SE3<S>, SE3<S>, S) {
        let mut u01 = || 0.5 * (rng.unif() + 1.0);
        let dir = |u: &mut dyn FnMut() -> f64| {
            let a = [u() - 0.5, u() - 0.5, u() - 0.5];
            let n = crate::laws::norm(&a).max(1e-3);
            a.map(|c| c / n)
        };
        let (a0, a, r0, r) = (dir(&mut u01), dir(&mut u01), dir(&mut u01), dir(&mut u01));
        let theta = match regime {
            0 => 10_f64.powf(-9.0 + 6.0 * u01()),
            1 => 1e-3 + (3.0 - 1e-3) * u01(),
            _ => core::f64::consts::PI - 10_f64.powf(-1.0 - 5.0 * u01()),
        };
        let (big, step, phi0) = (
            10_f64.powf(4.0 * u01()),
            10_f64.powf(2.0 * u01() - 1.0),
            3.0 * u01(),
        );
        const TS: [f64; 9] = [0.0, 1e-9, 0.25, 0.5, 1.0 - 1e-9, 1.0, -0.5, 1.5, 3.0];
        let t = TS.get(i % 12).copied().unwrap_or_else(u01);
        let c = |v: f64| S::constant(v);
        let x0 = SE3::<S>::exp(&Twist {
            phi: Vector(a0.map(|v| c(phi0 * v))),
            rho: [Vector(r0.map(|v| c(big * v)))],
        });
        let d = Twist {
            phi: Vector(a.map(|v| c(theta * v))),
            rho: [Vector(r.map(|v| c(step * v)))],
        };
        (x0, x0.rplus(&d), c(t))
    }

    /// The worst of `SE3::geodesic` against `reference::geodesic` over `n` draws per regime, in
    /// `u` of `S`: the quaternion sign-aligned, the translation against `max(‖x₀‖, ‖x₁‖, 1)`
    /// (`NUMERICS.md` §11's translation floor), since a fixed absolute bound cannot hold at
    /// `‖x₀‖ = 1e4` (`docs/maths/index.md`).
    fn screw_vs_reference<S: Sample>(seed: u64, n: usize) -> [[f64; 2]; 3] {
        let mut rng = crate::laws::Rng(seed);
        let f = |v: S| v.value_f64();
        array::from_fn(|regime| {
            let mut worst = [0.0; 2];
            for i in 0..n {
                let (x0, x1, t) = screw_draw::<S>(&mut rng, regime as u64, i);
                let (got, want) = (SE3::geodesic(&x0, &x1, t), reference::geodesic(&x0, &x1, t));
                let q = |g: &SE3<S>| {
                    let q = g.rotation().quat();
                    [f(q.w), f(q.x), f(q.y), f(q.z)]
                };
                let (qg, mut qw) = (q(&got), q(&want));
                if qg.iter().zip(&qw).map(|(a, b)| a * b).sum::<f64>() < 0.0 {
                    qw = qw.map(|v| -v);
                }
                let x = |g: &SE3<S>| g.translation().0.map(f);
                let scale = crate::laws::norm(&x(&x0)).max(crate::laws::norm(&x(&x1)));
                let e = crate::laws::worst(
                    crate::laws::e_at::<S>(&qg, &qw, 1.0),
                    crate::laws::e_at::<S>(&x(&got), &x(&want), scale),
                );
                let k = usize::from(!(0.0..=1.0).contains(&t.value_f64()));
                worst[k] = crate::laws::worst(worst[k], e);
            }
            worst
        })
    }

    /// `PHASE4.md` §1.2's `se3_geodesic_matches_reference`, `10⁵` pairs per precision (`0054`).
    ///
    /// The two twins share `Log`'s flip -- both read the sign of the one `q₀* q₁` -- so no band
    /// near `π` is excluded (`docs/maths/geodesics.md` GE.13(d)(i)). The bounds are twice the
    /// worst of `10⁶` draws per regime, rounded up (`measure_se3_geodesic_vs_reference`).
    #[test]
    fn se3_geodesic_matches_reference() {
        let seed = 0x0073_6372_6577;
        let (f64s, f32s) = (
            screw_vs_reference::<f64>(seed, 33_334),
            screw_vs_reference::<f32>(seed, 33_334),
        );
        for (name, got, bound) in [("f64", f64s, SCREW_TWIN_F64), ("f32", f32s, SCREW_TWIN_F32)] {
            for (regime, (g, b)) in got.iter().zip(bound).enumerate() {
                for (part, (g, b)) in ["t in [0, 1]", "t outside"].iter().zip(g.iter().zip(b)) {
                    assert!(
                        *g <= b,
                        "{name}, regime {regime}, {part}: {g} u exceeds {b} u"
                    );
                }
            }
        }
    }

    /// Twice the worst of `10⁶` draws per regime, rounded up -- consecutive, generic, near `π`,
    /// each as `[t in [0, 1], t outside]`: extrapolation is its own row because both twins' errors
    /// grow with `|t|`. Measured `[[23.516, 54.860], [24.147, 54.220], [28.194, 54.311]]` at `f64`
    /// and `[[35.925, 56.855], [21.285, 62.507], [24.221, 86.546]]` at `f32`.
    ///
    /// The difference is mostly the **reference's**: it rotates the translation out of `x₀`'s frame
    /// and back, which the world-frame twin does not. At `t = 1`, where the truth is `x₁` itself,
    /// the draws that differ by more than `12 u` read `2.6 u` mean (11.1 max) for the twin and
    /// `12.6 u` (35.9) for the reference (`0054`).
    const SCREW_TWIN_F64: [[f64; 2]; 3] = [[48.0, 110.0], [49.0, 109.0], [57.0, 109.0]];
    const SCREW_TWIN_F32: [[f64; 2]; 3] = [[72.0, 114.0], [43.0, 126.0], [49.0, 174.0]];

    /// `cargo nextest run -p helicoid --release --run-ignored only -- measure_se3_geodesic_vs_reference`.
    #[test]
    #[allow(clippy::print_stdout)]
    #[ignore = "measurement: prints the figures `SCREW_TWIN_*` are recorded from"]
    fn measure_se3_geodesic_vs_reference() {
        std::println!(
            "f64 {:?}",
            screw_vs_reference::<f64>(0x6d65_6173_7572, 1_000_000)
        );
        std::println!(
            "f32 {:?}",
            screw_vs_reference::<f32>(0x6d65_6173_7572, 1_000_000)
        );
    }

    /// The twin where `‖v‖²` is exact `0` (equal rotations), `0` by **underflow** while the
    /// rotations still differ (`n = 1e-170` at `f64`, `1e-25` at `f32`), and subnormal just above
    /// it: the short arm, so the value is the reference's to a few `u` and finite, at `f64`, `f32`
    /// and in `Dual`'s value lanes (`0054`). The derivative lanes are
    /// `the_screw_twin_differentiates_like_the_reference`'s.
    #[test]
    fn the_screw_twin_at_zero_and_underflowing_angles() {
        fn check<S: Sample>(n: f64, bound: f64) {
            let c = |v: f64| S::constant(v);
            let x0 =
                SE3::<S>::from_parts(crate::SO3::identity(), [Vector([c(2.0), c(-3.0), c(0.5)])]);
            // `q₁ = (1, n, 0, 0)` exactly, so `q₀* q₁` carries `n` without a rounding.
            let q1 = crate::Quat {
                w: c(1.0),
                x: c(n),
                y: c(0.0),
                z: c(0.0),
            };
            let x1 = SE3::<S>::from_parts(
                crate::SO3::from_quat_unchecked(q1),
                [Vector([
                    S::sample(1.5, 0),
                    S::sample(-2.0, 1),
                    S::sample(4.0, 2),
                ])],
            );
            for t in [0.25, 0.5, 1.0 - 1e-6, 1.0, 1.5, -0.5] {
                let t = c(t);
                let (got, want) = (SE3::geodesic(&x0, &x1, t), reference::geodesic(&x0, &x1, t));
                let (xg, xw) = (got.translation().0, want.translation().0);
                for (a, b) in xg.iter().zip(&xw) {
                    let (a, b) = (a.value_f64(), b.value_f64());
                    assert!(a.is_finite(), "n = {n:e}, t = {}: {a}", t.value_f64());
                    assert!(
                        (a - b).abs() <= bound * crate::laws::unit::<S>() * 5.0,
                        "n = {n:e}, t = {}: {a} against {b}",
                        t.value_f64()
                    );
                }
                let (qg, qw) = (got.rotation().quat(), want.rotation().quat());
                for (a, b) in [(qg.w, qw.w), (qg.x, qw.x)] {
                    assert!(
                        (a.value_f64() - b.value_f64()).abs() <= bound * crate::laws::unit::<S>()
                    );
                }
            }
        }
        for n in [0.0, 1e-150, 1e-160, 1e-170] {
            check::<f64>(n, 4.0);
            check::<Dual<f64, 3>>(n, 4.0);
        }
        for n in [0.0, 1e-15, 1e-20, 1e-25] {
            check::<f32>(n, 4.0);
        }
    }

    /// At `‖v‖ = 0` exactly the twin's derivative in `x₁`'s translation is the reference's.
    #[test]
    fn the_screw_twin_at_zero_angle_differentiates_like_the_reference() {
        type D = Dual<f64, 3>;
        let x0 = SE3::<D>::from_parts(
            crate::SO3::identity(),
            [Vector([
                D::constant(2.0),
                D::constant(-3.0),
                D::constant(0.5),
            ])],
        );
        let x1 = SE3::<D>::from_parts(
            crate::SO3::identity(),
            [Vector([
                D::variable(1.5, 0),
                D::variable(-2.0, 1),
                D::variable(4.0, 2),
            ])],
        );
        for t in [0.25, 0.5, 1.5] {
            let t = D::constant(t);
            let (got, want) = (SE3::geodesic(&x0, &x1, t), reference::geodesic(&x0, &x1, t));
            for (a, b) in got.translation().0.iter().zip(&want.translation().0) {
                for (da, db) in a.d.iter().zip(&b.d) {
                    assert!(
                        (da - db).abs() <= 4.0 * f64::EPSILON,
                        "{:?} against {:?}",
                        a.d,
                        b.d
                    );
                }
            }
        }
    }

    /// The twin's derivative against the reference's, through `Dual` with lanes on the relative
    /// tangent, per relative angle: the worst lane error relative to the largest reference lane, in
    /// `u`. `cargo nextest run -p helicoid --run-ignored only --no-capture -- measure_screw_twin_derivative`.
    fn screw_derivative_error(theta: f64, t: f64) -> f64 {
        type D = Dual<f64, 6>;
        let x0 = SE3::<D>::exp(&Twist {
            phi: Vector([0.3, -0.5, 0.4].map(D::constant)),
            rho: [Vector([2.0, -3.0, 0.5].map(D::constant))],
        });
        let axis = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
        let rho = [0.45, -0.6, 1.2];
        let d = Twist {
            phi: Vector(array::from_fn(|i| D::variable(theta * axis[i], i))),
            rho: [Vector(array::from_fn(|i| D::variable(rho[i], 3 + i)))],
        };
        let x1 = x0.rplus(&d);
        let t = D::constant(t);
        let (got, want) = (SE3::geodesic(&x0, &x1, t), reference::geodesic(&x0, &x1, t));
        let lanes = |g: &SE3<D>| {
            let q = g.rotation().quat();
            let mut out = std::vec::Vec::new();
            for c in [q.w, q.x, q.y, q.z].iter().chain(g.translation().0.iter()) {
                out.extend_from_slice(&c.d);
            }
            out
        };
        let (lg, lw) = (lanes(&got), lanes(&want));
        let scale = lw.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        lg.iter()
            .zip(&lw)
            .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()))
            / scale
            / f64::EPSILON
            * 2.0
    }

    /// The twin differentiates like the reference through `Dual`, at every angle (`0054`).
    ///
    /// Below `r`'s second switch the translation is the definition's closed form, so the lanes are
    /// the provided body's: measured at most 3.41 `u` from `θ = 1e-12` to `1e-2`. Above it, GE.12's
    /// `ϰ` grouping loses `≈ 10² α⁻¹ u` (GE.13(c)), which the switch bounds: 33.53 at `θ = 0.1`,
    /// 4.69 at `1`. Bounds twice those, rounded up. Before `0054`'s review the short side was GE.12
    /// too, and read `1.7e12 u` at `θ = 1e-12`.
    #[test]
    fn the_screw_twin_differentiates_like_the_reference() {
        for (thetas, bound) in [
            (&[1e-12, 1e-9, 1e-6, 1e-3, 1e-2][..], 7.0),
            (&[0.1, 1.0, 3.0][..], 68.0),
        ] {
            for &theta in thetas {
                for t in [0.25, 1.0 / 3.0, 0.7314, 1.0] {
                    let e = screw_derivative_error(theta, t);
                    assert!(
                        e <= bound,
                        "θ = {theta:e}, t = {t}: {e} u exceeds {bound} u"
                    );
                }
            }
        }
    }

    /// A pure screw about `z` through the origin: angle `θ`, displacement `h` along the axis. Its
    /// geodesic turns by `tθ` and advances `t h` -- the pitch is kept, which is the property the
    /// translation-first `LerpSlerp` lacks and the reason `ScLerp` is `tf_tree`'s default (D5).
    #[test]
    fn the_screw_twin_keeps_an_axial_pitch() {
        for theta in [1e-6, 0.1, 1.0, 3.0] {
            let x1 = SE3::<f64>::from_parts(
                crate::SO3::exp(&crate::SO3Tangent {
                    phi: Vector([0.0, 0.0, theta]),
                }),
                [Vector([0.0, 0.0, 2.5])],
            );
            for t in [0.25, 0.5, 0.75, 1.5] {
                let got = SE3::geodesic(&SE3::identity(), &x1, t);
                let want_q = crate::SO3::exp(&crate::SO3Tangent {
                    phi: Vector([0.0, 0.0, t * theta]),
                })
                .quat();
                let q = got.rotation().quat();
                let [x, y, z] = got.translation().0;
                assert!(
                    x.abs() <= 1e-15 && y.abs() <= 1e-15,
                    "θ = {theta}, t = {t}: ({x}, {y}) off the axis"
                );
                assert!(
                    (z - 2.5 * t).abs() <= 8.0 * f64::EPSILON,
                    "θ = {theta}, t = {t}: advanced {z}"
                );
                assert!(
                    (q.w - want_q.w).abs() <= 2.0 * f64::EPSILON
                        && (q.z - want_q.z).abs() <= 2.0 * f64::EPSILON
                );
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
        for (a, b) in cols.iter().zip(&cols2) {
            assert_eq!(a.0.map(f64::to_bits), b.0.map(f64::to_bits));
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
