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
/// one left-to-right sum of all `3 + 3N` products.
fn ops_match_dense<S: Lanes, const N: usize, const D: usize>(v: &[f64]) -> bool {
    let (a, b) = (
        tangent::<S, N, D>(&v[..D]),
        tangent::<S, N, D>(&v[D..2 * D]),
    );
    let k = S::sample(v[2 * D], 0);
    let (da, db) = (dense_of::<S, N, D>(&a), dense_of::<S, N, D>(&b));
    let dot = (1..D).fold(da[0] * db[0], |acc, i| acc + da[i] * db[i]);
    let round = dense_of::<S, N, D>(&SEn3Tangent::<S, N>::read_dense(&da));
    bits(&dense_of::<S, N, D>(&a.add(&b))) == bits(&array::from_fn(|i| da[i] + db[i]))
        && bits(&dense_of::<S, N, D>(&a.sub(&b))) == bits(&array::from_fn(|i| da[i] - db[i]))
        && bits(&dense_of::<S, N, D>(&a.neg())) == bits(&da.map(|x| -x))
        && bits(&dense_of::<S, N, D>(&a.scale(k))) == bits(&da.map(|x| x * k))
        && bits(&dense_of::<S, N, D>(&SEn3Tangent::<S, N>::zero())) == bits(&[S::zero(); D])
        && bits(&round) == bits(&da)
        && bits(&[a.dot(&b)]) == bits(&[dot])
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

#[test]
fn blend_selects_the_whole_tangent() {
    let a = SEn3Tangent::<f64, 2>::read_dense(&array::from_fn::<f64, 9, _>(|i| i as f64));
    let b = a.neg();
    for (m, want) in [(true, a), (false, b)] {
        let got = SEn3Tangent::blend(m, a, b);
        assert_eq!(
            fb(dense_of::<f64, 2, 9>(&got)),
            fb(dense_of::<f64, 2, 9>(&want))
        );
    }
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
/// `read_dense` reads a missing entry as `+0`.
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
    assert_eq!(
        (fb(back.phi.0), fb(back.rho[0].0)),
        (fb([1.0, 2.0, 3.0]), fb([4.0, 0.0, 0.0]))
    );
}
