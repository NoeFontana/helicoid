//! `Gaussian` (`docs/PHASE5.md` §5, `0065`): the round trip `to_right ∘ to_left` within the
//! sandwich bound, the `to_left`/`to_right` twins of `NUMERICS.md` §14, `mahalanobis_sq` against
//! a dense shadow and across sides (GG.10), the mask, and the stored symmetry.
//!
//! A case is a pose `Exp(τ)` with its translation scaled by a lever up to `10³` (GG.13's regime),
//! and `Σ = D₀ C D₀`, `D₀ = diag(σ_φ I, σ_ρ I)` over three decades each,
//! `C = (1 − ρ) I + ρ v vᵀ` with `‖v‖ = 1`, `ρ < 0.9`: `λ_min(C) = 1 − ρ` and
//! `λ_min(H) ≥ λ_min(C)` for the correlation matrix `H` (`diag C ≤ 1`), which is what GG.12's
//! bound for `d²` needs. The measured figures are `measure_gaussian`'s.

use crate::laws::{self, Sample};
use crate::{reference, Gaussian, Jac, Left, LieGroup, Product, Right, Rn, SEn3, Side, SO3};
use core::array;
use helicoid_linalg::{Mask, Matrix, Real, StridedMut, Vector};
use proptest::prelude::*;
use std::format;

type Se3<S> = SEn3<S, 1>;
type Se23<S> = SEn3<S, 2>;
type So3R3<S> = Product<SO3<S>, Rn<S, 3>>;

/// What one case is built from (module docs).
#[derive(Clone, Copy, Debug)]
struct Case<const D: usize> {
    pose: [f64; D],
    lever: f64,
    sphi: f64,
    srho: f64,
    v: [f64; D],
    rho: f64,
}

fn case<const D: usize>() -> impl Strategy<Value = Case<D>> {
    (
        laws::sample::<D>(),
        0_i32..=3,
        -3_i32..=0,
        -3_i32..=0,
        laws::sample::<D>(),
        0.0_f64..0.9,
    )
        .prop_map(|(pose, l, p, r, v, rho)| Case {
            pose,
            lever: 10_f64.powi(l),
            sphi: 10_f64.powi(p),
            srho: 10_f64.powi(r),
            v,
            rho,
        })
}

impl<const D: usize> Case<D> {
    fn from_rng(rng: &mut laws::Rng) -> Self {
        let pose = rng.shaped::<D>();
        let v = rng.shaped::<D>();
        let mut e = || (rng.next() % 4) as i32;
        let (l, p, r) = (e(), -e(), -e());
        Case {
            pose,
            lever: 10_f64.powi(l),
            sphi: 10_f64.powi(p),
            srho: 10_f64.powi(r),
            v,
            rho: 0.45 * (rng.unif() + 1.0),
        }
    }

    fn mean<S: Sample, G: LieGroup<S>>(&self) -> G {
        let tau: [f64; D] = array::from_fn(|i| {
            let s = if i < 3 { 1.0 } else { self.lever };
            self.pose[i] * s
        });
        G::exp(&laws::tangent::<S, G, D>(&tau))
    }

    /// `D₀ C D₀` in `f64`, then rounded to `S`.
    fn cov<S: Sample>(&self) -> Matrix<S, D, D> {
        let n = laws::norm(&self.v).max(f64::MIN_POSITIVE);
        let v: [f64; D] = array::from_fn(|i| self.v[i] / n);
        let d0 = |i: usize| if i < 3 { self.sphi } else { self.srho };
        Matrix::from_cols(array::from_fn(|c| {
            Vector(array::from_fn(|r| {
                let (i, j) = if r >= c { (r, c) } else { (c, r) };
                let cij = if i == j { 1.0 - self.rho } else { 0.0 } + self.rho * v[i] * v[j];
                S::constant(d0(i) * cij * d0(j))
            }))
        }))
    }

    fn gaussian<S: Sample, G: LieGroup<S>>(&self) -> Gaussian<S, G, Right, D> {
        Gaussian::new(self.mean::<S, G>(), self.cov::<S>())
    }
}

fn f64_of<S: Real, const D: usize>(m: &Matrix<S, D, D>) -> [[f64; D]; D] {
    array::from_fn(|r| array::from_fn(|c| m.get(r, c).value_f64()))
}

/// The dense `Ad_x` of `x`, row-major, in `f64`.
fn ad<S: Real, G: LieGroup<S>, const D: usize>(x: &G) -> [[f64; D]; D] {
    reference::dense::<S, G::Tangent, G::Jac, D>(&x.adjoint()).map(|r| r.map(S::value_f64))
}

fn abs_mul<const D: usize>(a: &[[f64; D]; D], b: &[[f64; D]; D]) -> [[f64; D]; D] {
    array::from_fn(|i| array::from_fn(|j| (0..D).map(|k| a[i][k].abs() * b[k][j].abs()).sum()))
}

/// `|m| |s| |m|ᵀ`.
fn abs_sandwich<const D: usize>(m: &[[f64; D]; D], s: &[[f64; D]; D]) -> [[f64; D]; D] {
    let ms = abs_mul(m, s);
    array::from_fn(|i| array::from_fn(|j| (0..D).map(|k| ms[i][k] * m[j][k].abs()).sum()))
}

/// `γ_k = k u/(1 − k u)` of `S`.
fn gamma<S: Real>(k: usize) -> f64 {
    let ku = k as f64 * laws::unit::<S>();
    ku / (1.0 - ku)
}

/// `max |got − want| / bound` over the lower triangle; `0/0` is `0`, NaN wins.
fn ratio<const D: usize>(got: &[[f64; D]; D], want: &[[f64; D]; D], bound: &[[f64; D]; D]) -> f64 {
    let mut worst = 0.0_f64;
    for i in 0..D {
        for j in 0..=i {
            let diff = (got[i][j] - want[i][j]).abs();
            let r = if diff == 0.0 { 0.0 } else { diff / bound[i][j] };
            worst = laws::worst(worst, r);
        }
    }
    worst
}

/// The `k` of the round-trip bound, `4D + ROUND_TRIP_AD`.
///
/// Componentwise (Higham Thm 3.5 for each product): `fl(fl(AΣ)Aᵀ) = AΣAᵀ + E`,
/// `|E| ≤ γ_{2D} |A||Σ||A|ᵀ`, for every order of the sums and with structural zeros only making it
/// smaller, so `SEn3Jac::sandwich` and `ProductJac::sandwich` are inside it. Two of those, with
/// `B = Ad(μ̂⁻¹)` after `A = Ad(μ̂)`, give `γ_{4D} M|Σ|Mᵀ`, `M = |B||A|`, about `BA Σ (BA)ᵀ`.
/// `BA = I + F` is not `I`: both are formed from rounded entries — `R` from the quaternion, the
/// translation's `[t]× R`, and the inverse's `−Rᵀ t` — and `|F| ≤ c u M` with `c` the roundings
/// of one entry of `A` (`R`: 3, `[t]× R`: 5), of `B` (`R`: 3, `−Rᵀ t`: 3, then `[·]× Rᵀ`: 5) and
/// the unit-norm defect of `q` (2): `c = 8 + 11 + 2 = 21`. `F` enters as `FΣ + ΣFᵀ` and
/// `M|Σ| ≤ M|Σ|Mᵀ` (`diag M ≥ 1`), so `k = 4D + 2c`.
const ROUND_TRIP_AD: usize = 42;

/// The round trip `to_right(to_left(Σ))` against `Σ`, as a fraction of its bound.
fn round_trip<S: Sample, G: LieGroup<S>, const D: usize>(c: &Case<D>) -> f64 {
    let g = c.gaussian::<S, G>();
    let back = g.to_left().to_right();
    let m = abs_mul(&ad::<S, G, D>(&g.mean.inverse()), &ad::<S, G, D>(&g.mean));
    let s = f64_of(g.cov());
    let k = gamma::<S>(4 * D + ROUND_TRIP_AD);
    let bound = abs_sandwich(&m, &s).map(|r| r.map(|x| k * x));
    ratio(&f64_of(back.cov()), &s, &bound)
}

/// `to_left` and `to_right` against `reference::sen3jac_sandwich` with `J = Ad`, as a fraction of
/// `2γ_{2D} |A||Σ||A|ᵀ`: both are within `γ_{2D}` of the exact product (see [`ROUND_TRIP_AD`]),
/// so their difference is within twice that; one bound for every block, and the per-block loss
/// of GG.13(c) is what it reduces to at a large translation.
fn twins<S: Sample, const N: usize, const D: usize>(c: &Case<D>) -> f64 {
    let g = c.gaussian::<S, SEn3<S, N>>();
    let twin = |x: &SEn3<S, N>, cov: &Matrix<S, D, D>| {
        let mut out = [[S::zero(); D]; D];
        reference::sen3jac_sandwich(
            &x.adjoint(),
            cov,
            &mut StridedMut::with_strides(out.as_flattened_mut(), D, D, D, 1),
        );
        let want = out.map(|r| r.map(S::value_f64));
        let bound = abs_sandwich(&ad::<S, SEn3<S, N>, D>(x), &f64_of(cov))
            .map(|r| r.map(|v| 2.0 * gamma::<S>(2 * D) * v));
        (want, bound)
    };
    let left = g.to_left();
    let (want, bound) = twin(&g.mean, g.cov());
    let l = ratio(&f64_of(left.cov()), &want, &bound);
    let (want, bound) = twin(&left.mean.inverse(), left.cov());
    let r = ratio(&f64_of(left.to_right().cov()), &want, &bound);
    laws::worst(l, r)
}

/// `δᵀ Σ⁻¹ δ` in `f64` by Cholesky, from `Σ`'s lower triangle.
fn shadow_d2<const D: usize>(s: &[[f64; D]; D], delta: &[f64; D]) -> f64 {
    let mut l = [[0.0_f64; D]; D];
    for j in 0..D {
        let d = s[j][j] - (0..j).map(|k| l[j][k] * l[j][k]).sum::<f64>();
        l[j][j] = d.sqrt();
        for i in j + 1..D {
            l[i][j] = (s[i][j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>()) / l[j][j];
        }
    }
    let mut y = [0.0_f64; D];
    for i in 0..D {
        y[i] = (delta[i] - (0..i).map(|k| l[i][k] * y[k]).sum::<f64>()) / l[i][i];
    }
    y.iter().map(|x| x * x).sum()
}

/// GG.12's `c_D`: `|d̂² − d²|/d² ≤ c_D u/λ_min(H)`, `c_D u = Dγ_{D+1} + 2Dγ_D + γ_D`.
fn c_d(d: usize) -> f64 {
    (3 * d * d + 2 * d) as f64
}

/// The point at the residual `ξ`, `‖ξ‖ = 1/2` along the case's `v` turned by one index. Not
/// smaller: `⊖` carries an absolute `O(u)` from the composition it takes the `Log` of, so a short
/// `ξ` measures that and not `mahalanobis_sq`.
fn point<S: Sample, G: LieGroup<S>, const D: usize>(c: &Case<D>, mean: &G) -> G {
    let n = 2.0 * laws::norm(&c.v).max(f64::MIN_POSITIVE);
    let xi: [f64; D] = array::from_fn(|i| c.v[(i + 1) % D] / n);
    mean.rplus(&laws::tangent::<S, G, D>(&xi))
}

/// `mahalanobis_sq` against [`shadow_d2`] of the same `δ` and the same stored `Σ`, as a fraction
/// of GG.12's bound (doubled at `f64`, where the shadow's own rounding is the subject's).
fn shadow<S: Sample, G: LieGroup<S>, Sd: Side, const D: usize>(
    g: &Gaussian<S, G, Sd, D>,
    x: &G,
    lambda: f64,
) -> f64 {
    let (d2, ok) = g.mahalanobis_sq(x);
    assert!(ok.all(), "mask clear on a positive definite Σ");
    let mut delta = [S::zero(); D];
    crate::Tangent::write_dense(&Sd::minus(x, &g.mean), &mut delta);
    let want = shadow_d2(&f64_of(g.cov()), &delta.map(S::value_f64));
    let twice = if laws::unit::<S>() < 1e-10 { 2.0 } else { 1.0 };
    let bound = twice * c_d(D) * laws::unit::<S>() / lambda;
    (d2.value_f64() - want).abs() / want / bound
}

/// `d²` of a point on the right and on the left (GG.10(a)), in `u`, relative. Well conditioned on
/// both sides: lever `1` and `σ_φ = σ_ρ`, so `λ_min(H_L)` is not GG.13(b)'s `σ_ρ²/(2σ_φ²‖t‖²)`.
fn sides<S: Sample, G: LieGroup<S>, const D: usize>(c: &Case<D>) -> f64 {
    let c = Case {
        lever: 1.0,
        sphi: c.srho,
        ..*c
    };
    let g = c.gaussian::<S, G>();
    let x = point::<S, G, D>(&c, &g.mean);
    let (r, l) = (g.mahalanobis_sq(&x).0, g.to_left().mahalanobis_sq(&x).0);
    (r.value_f64() - l.value_f64()).abs() / r.value_f64() / laws::unit::<S>()
}

fn symmetric<S: Real, const D: usize>(m: &Matrix<S, D, D>) -> bool {
    (0..D).all(|i| {
        (0..i).all(|j| m.get(i, j).value_f64().to_bits() == m.get(j, i).value_f64().to_bits())
    })
}

// Worst of `measure_gaussian` (10⁵ cases each), as a fraction of each bound — round trip `f64` /
// `f32`: SE3 0.212 / 0.207, SE23 0.175 / 0.195, SO3 0.269 / 0.245, SO3×R3 0.212 / 0.207 (without
// `ROUND_TRIP_AD` SO3 would read 1.21); twins SE3 0.345, SE23 0.277; shadow `f32` at most 0.151
// (`f64` 0: the shadow is the same recurrence in the same order). Sides, in `u`: SE3 32.8, SE23
// 25.1, SO3 16.9, SO3×R3 12.5.
const SIDES: f64 = 64.0;

macro_rules! gaussian_for {
    ($p:ident, $G:ident, $D:literal) => {
        mod $p {
            use super::*;
            proptest! {
                #[test]
                fn round_trip_within_sandwich_bound(c in case::<$D>()) {
                    laws::within(round_trip::<f64, $G<f64>, $D>(&c), 1.0)?;
                    laws::within(round_trip::<f32, $G<f32>, $D>(&c), 1.0)?;
                }
                #[test]
                fn mahalanobis_matches_dense_shadow(c in case::<$D>()) {
                    let lambda = 1.0 - c.rho;
                    let g = c.gaussian::<f64, $G<f64>>();
                    let x = point::<f64, $G<f64>, $D>(&c, &g.mean);
                    laws::within(shadow(&g, &x, lambda), 1.0)?;
                    let g = c.gaussian::<f32, $G<f32>>();
                    let x = point::<f32, $G<f32>, $D>(&c, &g.mean);
                    laws::within(shadow(&g, &x, lambda), 1.0)?;
                }
                #[test]
                fn mahalanobis_is_side_invariant(c in case::<$D>()) {
                    laws::within(sides::<f64, $G<f64>, $D>(&c), SIDES)?;
                }
                #[test]
                fn every_method_stores_a_symmetric_cov(c in case::<$D>()) {
                    let mut raw = c.cov::<f64>();
                    raw.set(0, $D - 1, 7.0); // an upper triangle that disagrees is not read
                    let g = Gaussian::<f64, $G<f64>, Right, $D>::new(c.mean::<f64, $G<f64>>(), raw);
                    prop_assert!(symmetric(g.cov()));
                    prop_assert!(g.cov().get(0, $D - 1).to_bits() == c.cov::<f64>().get($D - 1, 0).to_bits());
                    prop_assert!(symmetric(g.to_left().cov()));
                    prop_assert!(symmetric(g.to_left().to_right().cov()));
                    let j = g.mean.adjoint();
                    prop_assert!(symmetric(g.propagate(&j, g.mean).cov()));
                }
                /// `Whitener::mahalanobis_sq` is `mahalanobis_sq`'s value to the bit, and the
                /// mask is `chol`'s, on both sides and at both precisions (`0066`).
                #[test]
                fn the_whitener_is_mahalanobis_sq_to_the_bit(c in case::<$D>()) {
                    let g = c.gaussian::<f64, $G<f64>>();
                    let x = point::<f64, $G<f64>, $D>(&c, &g.mean);
                    let (w, ok) = g.whitener();
                    let (d2, mask) = g.mahalanobis_sq(&x);
                    prop_assert!(w.mahalanobis_sq(&x).to_bits() == d2.to_bits() && ok == mask);
                    let l = g.to_left();
                    let (w, _) = l.whitener();
                    prop_assert!(w.mahalanobis_sq(&x).to_bits() == l.mahalanobis_sq(&x).0.to_bits());
                    let g = c.gaussian::<f32, $G<f32>>();
                    let x = point::<f32, $G<f32>, $D>(&c, &g.mean);
                    prop_assert!(g.whitener().0.mahalanobis_sq(&x).to_bits() == g.mahalanobis_sq(&x).0.to_bits());
                }
                #[test]
                fn propagate_by_the_identity_keeps_cov(c in case::<$D>()) {
                    let g = c.gaussian::<f64, $G<f64>>();
                    let p = g.propagate(&<$G<f64> as LieGroup<f64>>::Jac::identity(), g.mean);
                    for i in 0..$D {
                        for k in 0..$D {
                            // `(a − b)` is `+0` exactly when `a == b`, signed zeros included.
                            prop_assert!((p.cov().get(i, k) - g.cov().get(i, k)).abs().to_bits() == 0);
                        }
                    }
                }
            }

            #[test]
            fn an_indefinite_cov_clears_the_mask() {
                let mut cov = Matrix::<f64, $D, $D>::identity();
                cov.set($D - 1, $D - 1, -1.0);
                let g = Gaussian::<f64, $G<f64>, Left, $D>::new(<$G<f64>>::identity(), cov);
                assert!(!g.mahalanobis_sq(&<$G<f64>>::identity()).1);
            }
        }
    };
}

gaussian_for!(se3, Se3, 6);
gaussian_for!(se23, Se23, 9);
gaussian_for!(so3, SO3, 3);
gaussian_for!(so3_r3, So3R3, 6);

proptest! {
    #[test]
    fn to_left_and_to_right_match_reference_se3(c in case::<6>()) {
        laws::within(twins::<f64, 1, 6>(&c), 1.0)?;
        laws::within(twins::<f32, 1, 6>(&c), 1.0)?;
    }
    #[test]
    fn to_left_and_to_right_match_reference_se23(c in case::<9>()) {
        laws::within(twins::<f64, 2, 9>(&c), 1.0)?;
        laws::within(twins::<f32, 2, 9>(&c), 1.0)?;
    }
}

/// The worst fraction of each bound over `10⁵` seeded cases; the figures the header records.
#[test]
#[ignore = "measurement: prints the figures the bounds are recorded from"]
#[allow(clippy::print_stdout)]
fn measure_gaussian() {
    fn group<G64: LieGroup<f64>, G32: LieGroup<f32>, const D: usize>(
        name: &str,
    ) -> std::string::String {
        let mut rng = laws::Rng(0x6761_7573_7300_0001);
        let mut w = [0.0_f64; 5];
        for _ in 0..100_000 {
            let c = Case::<D>::from_rng(&mut rng);
            w[0] = laws::worst(w[0], round_trip::<f64, G64, D>(&c));
            w[1] = laws::worst(w[1], round_trip::<f32, G32, D>(&c));
            let lambda = 1.0 - c.rho;
            let g = c.gaussian::<f64, G64>();
            w[2] = laws::worst(w[2], shadow(&g, &point::<f64, G64, D>(&c, &g.mean), lambda));
            let g = c.gaussian::<f32, G32>();
            w[3] = laws::worst(w[3], shadow(&g, &point::<f32, G32, D>(&c, &g.mean), lambda));
            w[4] = laws::worst(w[4], sides::<f64, G64, D>(&c));
        }
        format!(
            "{name}: round trip f64 {:.4} f32 {:.4}; shadow f64 {:.4} f32 {:.4}; sides {:.3} u",
            w[0], w[1], w[2], w[3], w[4]
        )
    }
    std::println!("{}", group::<Se3<f64>, Se3<f32>, 6>("se3"));
    std::println!("{}", group::<Se23<f64>, Se23<f32>, 9>("se23"));
    std::println!("{}", group::<SO3<f64>, SO3<f32>, 3>("so3"));
    std::println!("{}", group::<So3R3<f64>, So3R3<f32>, 6>("so3_r3"));
    let mut rng = laws::Rng(0x6761_7573_7300_0002);
    let mut t = [0.0_f64; 2];
    for _ in 0..100_000 {
        let c6 = Case::<6>::from_rng(&mut rng);
        let c9 = Case::<9>::from_rng(&mut rng);
        t[0] = laws::worst(
            t[0],
            laws::worst(twins::<f64, 1, 6>(&c6), twins::<f32, 1, 6>(&c6)),
        );
        t[1] = laws::worst(
            t[1],
            laws::worst(twins::<f64, 2, 9>(&c9), twins::<f32, 2, 9>(&c9)),
        );
    }
    std::println!("twins: se3 {:.4}, se23 {:.4}", t[0], t[1]);
}
