//! The seeded SE_N(3) subject: `Exp` (`docs/NUMERICS.md` §5.1) and the Jacobians `J_l`, `J_r` with
//! Barfoot's block `Q` (§5.3) over `sen3_{exp,jr,jl}_n{1,2,3}`, and the two planted defects of
//! `docs/PHASE1.md` §10 (`Order::TranslationFirst`, the `half` of [`jacobian`]). Generic over
//! `Real`; the harness runs it at `f64`. Tangents are rotation-first `[φ; ρ₁; …; ρ_N]`, the dense
//! `J` is column-major in that order (`conformance/generate/README.md`).
//!
//! Readings where §5.1, §5.3 and §4 are silent or differ, each the smallest:
//!
//! - **Coefficients** `a, b, d, e` (and `k` for the rotation) are each evaluated at their own
//!   generated switch, where §4 groups them per call site inside one `branch` (`jr_coeffs`,
//!   `q_coeffs`): the sweep generates one switch per coefficient, not per group (`docs/PHASE1.md`
//!   §0.0, *Missing*). Every error measured here is of this kernel, not of the grouped one:
//!   0014 (draft) question 24.
//! - **The rotation** of `Exp` is `so3::exp` itself, `cos θ/2` reading included (0014 (draft)
//!   question 19).
//! - **Products** are 3×3 matrix products, entry `(r, c)` summed `(t₀ + t₁) + t₂`, each word of §5.3
//!   associated left to right as it is written, and `J_l(φ) = I + aW + bW²` is formed as a matrix
//!   (`W²` is `W·W`, not `φφᵀ − θ²I`) before it multiplies `ρ_i` (question 25).
//! - **`J_r(τ)` is `J_l(−τ)`** (`docs/maths/se3.md` SE.9(a)): §5.3's `Q(−ρ_i, −φ)` and
//!   `J_r(φ) = J_l(−φ)`, so the two sides share every operation and differ by the sign of the
//!   tangent, which is exact.
//! - **Translation-first** is `[ρ₁; …; ρ_N; φ]`, the rotation block moved last for every `N`
//!   (`docs/maths/se3.md` SE.14(c)), so that defect is planted at `N = 1, 2, 3`.

use helicoid_linalg::Real;

use super::kernel::{coefficient, Candidate, Coeff};
use super::so3::{self, norm_sq};

/// What a `sen3_*` id computes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Op {
    Exp,
    /// `J_r`, block `Q(−ρ_i, −φ)`.
    Jr,
    /// `J_l`, block `Q(ρ_i, φ)`.
    Jl,
}

/// `(op, N)` of `sen3_<op>_n<N>` for the ids this subject computes.
pub(crate) fn parse(fn_id: &str) -> Option<(Op, usize)> {
    let (op, n) = fn_id.strip_prefix("sen3_")?.rsplit_once("_n")?;
    let op = match op {
        "exp" => Op::Exp,
        "jr" => Op::Jr,
        "jl" => Op::Jl,
        _ => return None,
    };
    let n = match n {
        "1" => 1,
        "2" => 2,
        "3" => 3,
        _ => return None,
    };
    Some((op, n))
}

/// Where the rotation block sits in the tangent a routine reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Order {
    /// `[φ; ρ₁; …; ρ_N]` (`docs/NUMERICS.md` §1).
    RotationFirst,
    /// `[ρ₁; …; ρ_N; φ]` (`docs/maths/se3.md` SE.14): the planted defect reads it so.
    TranslationFirst,
}

/// The generated kernel of each coefficient, by [`Coeff::index`].
pub(crate) struct Kernels<'a, S> {
    pub(crate) arms: [(Candidate<S>, &'a [S]); 6],
}

impl<S: Real> Kernels<'_, S> {
    fn at(&self, c: Coeff, z: S) -> S {
        let (candidate, series) = self.arms[c.index()];
        coefficient(c, z, candidate, series)
    }
}

/// A 3×3 matrix by rows.
type Mat<S> = [[S; 3]; 3];

fn hat<S: Real>([x, y, z]: [S; 3]) -> Mat<S> {
    let o = S::zero();
    [[o, -z, y], [z, o, -x], [-y, x, o]]
}

fn identity<S: Real>() -> Mat<S> {
    let (o, i) = (S::zero(), S::one());
    [[i, o, o], [o, i, o], [o, o, i]]
}

fn mul<S: Real>(a: &Mat<S>, b: &Mat<S>) -> Mat<S> {
    let entry = |r: usize, c: usize| (a[r][0] * b[0][c] + a[r][1] * b[1][c]) + a[r][2] * b[2][c];
    std::array::from_fn(|r| std::array::from_fn(|c| entry(r, c)))
}

fn apply<S: Real>(a: &Mat<S>, v: [S; 3]) -> [S; 3] {
    a.map(|row| (row[0] * v[0] + row[1] * v[1]) + row[2] * v[2])
}

fn zip<S: Real>(a: &Mat<S>, b: &Mat<S>, f: impl Fn(S, S) -> S) -> Mat<S> {
    std::array::from_fn(|r| std::array::from_fn(|c| f(a[r][c], b[r][c])))
}

fn add<S: Real>(a: &Mat<S>, b: &Mat<S>) -> Mat<S> {
    zip(a, b, |x, y| x + y)
}

fn scale<S: Real>(s: S, a: &Mat<S>) -> Mat<S> {
    a.map(|row| row.map(|x| s * x))
}

/// `J_l(φ) = I + a W + b W²` (`NUMERICS.md` §3.5), `W = φ^`.
fn jl_so3<S: Real>(w: &Mat<S>, a: S, b: S) -> Mat<S> {
    add(&add(&identity(), &scale(a, w)), &scale(b, &mul(w, w)))
}

/// `Q(ρ, φ)` of `NUMERICS.md` §5.3 with `half` for its `½` (`±½`: the planted defect).
fn q_block<S: Real>(w: &Mat<S>, rho: [S; 3], [b, d, e]: [S; 3], half: S) -> Mat<S> {
    let x = hat(rho);
    let (wx, xw, ww) = (mul(w, &x), mul(&x, w), mul(w, w));
    let (wxw, wwx, xww) = (mul(&wx, w), mul(&ww, &x), mul(&xw, w));
    let b_words = add(&add(&wx, &xw), &wxw);
    let d_words = zip(&add(&wwx, &xww), &wxw, |s, t| s - S::lit(3.0) * t);
    let e_words = add(&mul(&wxw, w), &mul(&wwx, w));
    let head = add(&scale(half, &x), &scale(b, &b_words));
    add(&add(&head, &scale(d, &d_words)), &scale(e, &e_words))
}

/// `φ` and the `ρ_i` of `tau`; `None` unless `tau` holds `3 + 3N` entries, `N ≥ 1`.
fn split<S: Real>(tau: &[S], order: Order) -> Option<([S; 3], Vec<[S; 3]>)> {
    let chunks = tau.chunks_exact(3);
    if tau.len() < 6 || !chunks.remainder().is_empty() {
        return None;
    }
    let blocks: Vec<[S; 3]> = chunks.filter_map(|c| <[S; 3]>::try_from(c).ok()).collect();
    let (phi, rho) = match order {
        Order::RotationFirst => blocks.split_first()?,
        Order::TranslationFirst => blocks.split_last()?,
    };
    Some((*phi, rho.to_vec()))
}

/// `Exp(τ) = (Exp(φ), J_l(φ) ρ₁, …, J_l(φ) ρ_N)` (§5.1): `q` and the flat `[x₁; …; x_N]`.
pub(crate) fn exp<S: Real>(
    tau: &[S],
    kernels: &Kernels<S>,
    order: Order,
) -> Option<([S; 4], Vec<S>)> {
    let (phi, rho) = split(tau, order)?;
    let (candidate, series) = kernels.arms[Coeff::K.index()];
    let q = so3::exp(phi, candidate, series);
    let z = norm_sq(phi);
    let j = jl_so3(&hat(phi), kernels.at(Coeff::A, z), kernels.at(Coeff::B, z));
    Some((q, rho.iter().flat_map(|&r| apply(&j, r)).collect()))
}

/// `J_l(τ)` or, `right`, `J_r(τ) = J_l(−τ)` (§5.3) as a dense column-major `(3+3N)`-square matrix:
/// `J(φ)` on every diagonal block, `Q` in block `(i, 0)`, exact zeros elsewhere.
pub(crate) fn jacobian<S: Real>(
    tau: &[S],
    right: bool,
    kernels: &Kernels<S>,
    half: S,
) -> Option<Vec<S>> {
    let signed: Vec<S> = tau.iter().map(|&t| if right { -t } else { t }).collect();
    let (phi, rho) = split(&signed, Order::RotationFirst)?;
    let z = norm_sq(phi);
    let [a, b, d, e] = [Coeff::A, Coeff::B, Coeff::D, Coeff::E].map(|c| kernels.at(c, z));
    let w = hat(phi);
    let diagonal = jl_so3(&w, a, b);
    let m = 3 + 3 * rho.len();
    let mut dense = vec![S::zero(); m * m];
    let mut put = |block: (usize, usize), values: &Mat<S>| {
        for (r, row) in values.iter().enumerate() {
            for (c, &x) in row.iter().enumerate() {
                dense[(3 * block.1 + c) * m + 3 * block.0 + r] = x;
            }
        }
    };
    put((0, 0), &diagonal);
    for (i, &r) in rho.iter().enumerate() {
        put((i + 1, i + 1), &diagonal);
        put((i + 1, 0), &q_block(&w, r, [b, d, e], half));
    }
    Some(dense)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seeded::{Seeded, D1};

    const BAD_LENGTH: &str = "a tangent of the wrong length";

    fn constants(v: &[f64]) -> Vec<D1> {
        v.iter().map(|&x| D1::constant(x)).collect()
    }

    fn values(v: &[D1]) -> Vec<f64> {
        v.iter().map(|c| c.v).collect()
    }

    fn exp64(tau: &[f64], order: Order) -> Result<([f64; 4], Vec<f64>), String> {
        let kernels = Seeded::generated();
        let (q, x) = exp(&constants(tau), &kernels.kernels(), order).ok_or(BAD_LENGTH)?;
        Ok((q.map(|c| c.v), values(&x)))
    }

    fn jac64(tau: &[f64], right: bool, half: f64) -> Result<Vec<f64>, String> {
        let generated = Seeded::generated();
        let half = D1::constant(half);
        let j = jacobian(&constants(tau), right, &generated.kernels(), half);
        Ok(values(&j.ok_or(BAD_LENGTH)?))
    }

    /// `[φ; ρ₁; ρ₂]` with `|φ| = theta`.
    fn tangent(theta: f64) -> Vec<f64> {
        let phi = [theta / 3.0, 2.0 * theta / 3.0, 2.0 * theta / 3.0];
        [phi, [0.3, -1.2, 0.5], [2.0, 0.1, -0.7]].concat()
    }

    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    fn close(a: &[f64], b: &[f64], tol: f64) {
        let err = a.iter().zip(b).map(|(x, y)| (x - y).abs());
        assert!(
            a.len() == b.len() && err.fold(0.0, f64::max) <= tol,
            "{a:?} vs {b:?}"
        );
    }

    #[test]
    fn exp_is_the_so3_exp_and_j_l_of_phi_applied_to_every_rho() -> Result<(), String> {
        for theta in [0.7, 2.0] {
            let tau = tangent(theta);
            let (q, x) = exp64(&tau, Order::RotationFirst)?;
            let (s, c) = (theta / 2.0).sin_cos();
            close(&q, &[c, s / 3.0, 2.0 * s / 3.0, 2.0 * s / 3.0], 4e-16);
            // `J_l ρ = ρ + a φ×ρ + b φ×(φ×ρ)`, `a` and `b` by the standard library at these θ.
            let a = (1.0 - theta.cos()) / theta.powi(2);
            let b = (theta - theta.sin()) / theta.powi(3);
            let phi = [tau[0], tau[1], tau[2]];
            for (i, rho) in tau[3..].chunks_exact(3).enumerate() {
                let rho = [rho[0], rho[1], rho[2]];
                let (c1, c2) = (cross(phi, rho), cross(phi, cross(phi, rho)));
                let want: Vec<f64> = (0..3).map(|k| rho[k] + a * c1[k] + b * c2[k]).collect();
                close(&x[3 * i..3 * i + 3], &want, 1e-14);
            }
        }
        let at_zero = exp64(&[0.0, 0.0, 0.0, 1.5, -2.0, 0.25], Order::RotationFirst)?;
        assert_eq!(at_zero, ([1.0, 0.0, 0.0, 0.0], vec![1.5, -2.0, 0.25]));
        for wrong in [&[0.0; 5][..], &[0.0; 3], &[0.0; 7]] {
            assert!(exp64(wrong, Order::RotationFirst).is_err());
            assert!(jac64(wrong, false, 0.5).is_err());
        }
        Ok(())
    }

    /// Column-major `m`-square product.
    fn matmul(a: &[f64], b: &[f64], m: usize) -> Vec<f64> {
        let entry = |k: usize| (0..m).map(|j| a[j * m + k % m] * b[(k / m) * m + j]).sum();
        (0..m * m).map(entry).collect()
    }

    /// `Σ_{n<40} (±ad_τ)ⁿ/(n+1)!` of the dense rotation-first `ad_τ = W + ε[ρ_i]×` (`NUMERICS.md`
    /// §5.2) in plain `f64`, column-major: no `Q`, no block bookkeeping.
    fn series(tau: &[f64], sign: f64) -> Vec<f64> {
        let m = tau.len();
        let mut ad = vec![0.0; m * m];
        let block = |i: usize| [tau[3 * i], tau[3 * i + 1], tau[3 * i + 2]];
        let mut put = |v: [f64; 3], (br, bc): (usize, usize)| {
            for (r, row) in hat(v).iter().enumerate() {
                for (c, &x) in row.iter().enumerate() {
                    ad[(3 * bc + c) * m + 3 * br + r] = sign * x;
                }
            }
        };
        for i in 0..m / 3 {
            put(block(0), (i, i));
            if i > 0 {
                put(block(i), (i, 0));
            }
        }
        let mut power = vec![0.0; m * m];
        (0..m).for_each(|i| power[i * m + i] = 1.0);
        let (mut sum, mut factorial) = (vec![0.0; m * m], 1.0);
        for n in 0..40 {
            sum.iter_mut()
                .zip(&power)
                .for_each(|(s, p)| *s += p / factorial);
            power = matmul(&power, &ad, m);
            factorial *= f64::from(n + 2);
        }
        sum
    }

    #[test]
    fn the_jacobians_are_the_defining_series_of_the_dense_ad() -> Result<(), String> {
        // Both arms of `b, d, e` (their switches are at θ = 0.98) and both sides, N = 2.
        for theta in [0.5, 0.95, 1.4, 2.6] {
            let tau = tangent(theta);
            close(&jac64(&tau, false, 0.5)?, &series(&tau, 1.0), 1e-13);
            close(&jac64(&tau, true, 0.5)?, &series(&tau, -1.0), 1e-13);
        }
        Ok(())
    }

    #[test]
    fn the_defects_change_one_thing_each() -> Result<(), String> {
        let tau = tangent(1.3)[..6].to_vec();
        // Translation-first is the correct `Exp` of the tangent with its two blocks swapped.
        let swapped = [&tau[3..], &tau[..3]].concat();
        let bad = exp64(&tau, Order::TranslationFirst)?;
        let good = exp64(&tau, Order::RotationFirst)?;
        let bits = |e: ([f64; 4], Vec<f64>)| {
            let all = [e.0.to_vec(), e.1].concat();
            all.into_iter().map(f64::to_bits).collect::<Vec<_>>()
        };
        assert_eq!(
            bits(bad.clone()),
            bits(exp64(&swapped, Order::RotationFirst)?)
        );
        assert!(bad.1.iter().zip(&good.1).any(|(a, b)| (a - b).abs() > 0.1));
        // `-½` for `+½` moves block `(1, 0)` by `∓ρ^` (`J_r` has `Q(-ρ, -φ)`) and nothing else.
        let rho = [tau[3], tau[4], tau[5]];
        let hat = hat(rho);
        for (right, sign) in [(false, -1.0), (true, 1.0)] {
            let (bad, good) = (jac64(&tau, right, -0.5)?, jac64(&tau, right, 0.5)?);
            for (k, (b, g)) in bad.iter().zip(&good).enumerate() {
                let (r, c) = (k % 6, k / 6);
                let want = if r >= 3 && c < 3 {
                    sign * hat[r - 3][c]
                } else {
                    0.0
                };
                assert!((b - g - want).abs() < 1e-15, "{right} ({r}, {c})");
            }
        }
        Ok(())
    }
}
