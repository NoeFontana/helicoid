//! The seeded spline (`0116` decision 5): `NUMERICS.md` §21 written out on `helicoid`'s groups,
//! the clean subject the planted spline defects are measured against, and those defects.
//!
//! - **Value, velocity, acceleration** run one generic recursion over `LieGroup`'s `exp`, `log`
//!   (through `rminus`), `adjoint` and `ad`, on SO(3), SE(3) and `Rn<_, 3>`.
//! - **The SO(3) Jacobians** are §21's forward and backward passes in dense `Mat3` arithmetic, which
//!   SO(3)'s `Jac = Mat3` allows, so they share nothing with `helicoid-spline`'s generic
//!   accumulation through `Jac::add` (`0112`).
//! - **R^n's Jacobian** is `k = λ_m − λ_{m+1}`, the recursion's per-control sum at `Ad = J = I`.
//! - **The Gram matrix** is §21's sum reduced to one fraction per entry, `S::lit(num) / S::lit(den)`,
//!   then `h⁻¹` raised to `2r − 1` by repeated products.
//!
//! A derivative output is scaled to unit size as its id states: `ω / h⁻¹`, `ω̇ / (h⁻¹ h⁻¹)`, the
//! product formed as `λ̈`'s is. Every input reaches the group exactly (`shipped::scalars`).

use helicoid::{Jac, LieGroup, Tangent, SO3};
use helicoid_linalg::{Mat3, Precision, Real, StridedMut};

use crate::conformance::corpus::Record;
use crate::conformance::subject::{Output, Registered, Subject};
use crate::shipped::scalars;
use crate::shipped::spline::{dense, points, poses, quat_of, rotations, scaled, segment};

/// A planted defect of `0116` decision 5, each one change to the clean evaluator.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Planted {
    /// `λ_j` from the non-cumulative basis `b_j` (the rows `n_{j,·}`, not their sums).
    NonCumulative,
    /// `∂d_j/∂δ_{j−1} = −J_r⁻¹(d_j)` where §21 has `−J_l⁻¹(d_j)`.
    JrInvForJlInv,
    /// The `ad` term of the `ω̇` recursion subtracted.
    AdSignFlipped,
    /// `R_j = Ad_{A_j}` where §21 has `Ad_{A_j⁻¹}`.
    AdNotInverted,
    /// `λ̈` without its `h⁻²`.
    NoInvHSquared,
    /// The Gram matrix without its `h^{1−2r}`.
    GramUnscaled,
}

impl Planted {
    pub(crate) const ALL: [Planted; 6] = [
        Planted::NonCumulative,
        Planted::JrInvForJlInv,
        Planted::AdSignFlipped,
        Planted::AdNotInverted,
        Planted::NoInvHSquared,
        Planted::GramUnscaled,
    ];

    pub(crate) fn name(self) -> &'static str {
        match self {
            Planted::NonCumulative => "seeded:spline-non-cumulative",
            Planted::JrInvForJlInv => "seeded:spline-jr-inv-for-jl-inv",
            Planted::AdSignFlipped => "seeded:spline-ad-sign-flipped",
            Planted::AdNotInverted => "seeded:spline-ad-not-inverted",
            Planted::NoInvHSquared => "seeded:spline-no-inv-h-squared",
            Planted::GramUnscaled => "seeded:spline-gram-unscaled",
        }
    }

    /// The ids the defect is run over, decision 5's.
    pub(crate) fn ids(self) -> &'static [&'static str] {
        match self {
            Planted::NonCumulative => &["so3_spline_value", "rn_spline_value"],
            Planted::JrInvForJlInv => &["so3_spline_value_jac"],
            Planted::AdSignFlipped => &["so3_spline_acc"],
            Planted::AdNotInverted => &["so3_spline_vel"],
            Planted::NoInvHSquared => &["so3_spline_acc", "rn_spline_acc"],
            Planted::GramUnscaled => &["rn_spline_gram"],
        }
    }
}

/// Every id the clean evaluator answers.
pub(crate) const IDS: [&str; 16] = [
    "so3_spline_value",
    "so3_spline_vel",
    "so3_spline_acc",
    "so3_spline_value_jac",
    "so3_spline_vel_jac",
    "so3_spline_acc_jac",
    "rn_spline_value",
    "rn_spline_vel",
    "rn_spline_acc",
    "rn_spline_value_jac",
    "rn_spline_vel_jac",
    "rn_spline_acc_jac",
    "se3_spline_value",
    "se3_spline_vel",
    "se3_spline_acc",
    "rn_spline_gram",
];

/// `C(n, k)`.
fn choose(n: i64, k: i64) -> i64 {
    (0..k).fold(1, |acc, i| acc * (n - i) / (i + 1))
}

fn factorial(n: i64) -> i64 {
    (1..=n).product()
}

/// §21's `n_{s,c} = C(K−1, c) Σ_{l=s}^{K−1} (−1)^{l−s} C(K, l−s) (K−1−l)^{K−1−c}`, an integer.
fn numerators(k: usize) -> Vec<Vec<i64>> {
    let k = k as i64;
    let pow = |b: i64, e: i64| (0..e).fold(1, |acc, _| acc * b);
    (0..k)
        .map(|s| {
            (0..k)
                .map(|c| {
                    let sum: i64 = (s..k)
                        .map(|l| {
                            let sign = if (l - s) % 2 == 0 { 1 } else { -1 };
                            sign * choose(k, l - s) * pow(k - 1 - l, k - 1 - c)
                        })
                        .sum();
                    choose(k - 1, c) * sum
                })
                .collect()
        })
        .collect()
}

/// `λ_j`, `λ̇_j`, `λ̈_j`, `j = 0 … K−1`, by §21's Horner on `S::lit` numerators.
pub(crate) struct Weights<S> {
    pub(crate) lam: Vec<S>,
    pub(crate) dot: Vec<S>,
    pub(crate) ddot: Vec<S>,
}

pub(crate) fn weights<S: Real>(k: usize, u: S, inv_h: S, planted: Option<Planted>) -> Weights<S> {
    let n = numerators(k);
    let fact = S::lit(factorial(k as i64 - 1) as f64);
    let (mut lam, mut dot, mut ddot) = (vec![S::one()], vec![S::zero()], vec![S::zero()]);
    for j in 1..k {
        // The cumulative numerators, or the defect's single row.
        let row: Vec<i64> = match planted {
            Some(Planted::NonCumulative) => n[j].clone(),
            _ => (0..k).map(|c| (j..k).map(|s| n[s][c]).sum()).collect(),
        };
        let horner = |coeff: &dyn Fn(usize) -> i64, lowest: usize| {
            let mut acc = S::zero();
            for c in (lowest..k).rev() {
                acc = acc * u + S::lit(coeff(c) as f64);
            }
            acc
        };
        lam.push(horner(&|c| row[c], 0) / fact);
        dot.push(inv_h * (horner(&|c| c as i64 * row[c], 1) / fact));
        let second = horner(&|c| (c as i64) * (c as i64 - 1) * row[c], 2) / fact;
        ddot.push(match planted {
            Some(Planted::NoInvHSquared) => second,
            _ => (inv_h * inv_h) * second,
        });
    }
    Weights { lam, dot, ddot }
}

/// `(X, ω, ω̇)` of the window, §21's recursion.
fn kinematics<S: Real, G: LieGroup<S>>(
    x: &[G],
    w: &Weights<S>,
    planted: Option<Planted>,
) -> (G, G::Tangent, G::Tangent) {
    let (mut value, mut om, mut al) = (x[0], G::Tangent::zero(), G::Tangent::zero());
    for j in 1..x.len() {
        let d = x[j].rminus(&x[j - 1]);
        let a = G::exp(&d.scale(w.lam[j]));
        value = value * a;
        let r = match planted {
            Some(Planted::AdNotInverted) => a.adjoint(),
            _ => a.inverse().adjoint(),
        };
        om = r.apply(&om).add(&d.scale(w.dot[j]));
        let ad = G::ad(&om).apply(&d).scale(w.dot[j]);
        let rest = r.apply(&al).add(&d.scale(w.ddot[j]));
        al = match planted {
            Some(Planted::AdSignFlipped) => rest.sub(&ad),
            _ => rest.add(&ad),
        };
    }
    (value, om, al)
}

fn zero3<S: Real>() -> Mat3<S> {
    Mat3::identity().scale(S::zero())
}

fn hat<S: Real>(t: &<SO3<S> as LieGroup<S>>::Tangent) -> Mat3<S> {
    SO3::<S>::ad(t)
}

/// `[J_0 … J_{K−1}]` of the output of `order` (0 value, 1 `ω`, 2 `ω̇`), §21's right Jacobians.
fn so3_jacobians<S: Real>(
    x: &[SO3<S>],
    w: &Weights<S>,
    order: usize,
    planted: Option<Planted>,
) -> Vec<Mat3<S>> {
    type T<S> = <SO3<S> as LieGroup<S>>::Tangent;
    let k = x.len();
    let one = Mat3::<S>::identity();
    // Forward: R_j, ω^(j), ω̇^(j), J_j and the per-increment partials.
    let mut d: Vec<T<S>> = vec![T::<S>::zero()];
    let mut r = vec![one];
    let (mut om, mut al) = (vec![T::<S>::zero()], vec![T::<S>::zero()]);
    let (mut p_om, mut p_al, mut jj) = (vec![zero3()], vec![zero3()], vec![one]);
    for j in 1..k {
        let dj = x[j].rminus(&x[j - 1]);
        let ld = dj.scale(w.lam[j]);
        let a = SO3::exp(&ld);
        let rj = a.inverse().adjoint();
        let jr = SO3::jr(&ld);
        let rw = rj.apply(&om[j - 1]);
        let ra = rj.apply(&al[j - 1]);
        let omj = rw.add(&dj.scale(w.dot[j]));
        let alj = ra
            .add(&dj.scale(w.ddot[j]))
            .add(&hat(&omj).apply(&dj).scale(w.dot[j]));
        let po = (hat(&rw) * jr).scale(w.lam[j]) + one.scale(w.dot[j]);
        let pa = (hat(&ra) * jr).scale(w.lam[j]) + one.scale(w.ddot[j]) + hat(&omj).scale(w.dot[j])
            - (hat(&dj) * po).scale(w.dot[j]);
        d.push(dj);
        r.push(rj);
        om.push(omj);
        al.push(alj);
        p_om.push(po);
        p_al.push(pa);
        jj.push(jr);
    }
    // Backward: T_j, C_j and the partials in each increment.
    let (mut t, mut c) = (one, zero3::<S>());
    let mut by_d = vec![zero3::<S>(); k];
    for j in (1..k).rev() {
        by_d[j] = match order {
            0 => (t * jj[j]).scale(w.lam[j]),
            1 => t * p_om[j],
            _ => c * p_om[j] + t * p_al[j],
        };
        c = c * r[j] - ((t * hat(&d[j])) * r[j]).scale(w.dot[j]);
        t = t * r[j];
    }
    // Per control: the term at j = m, then the one at j = m + 1 added to it.
    (0..k)
        .map(|m| {
            let mut acc = match (m, order) {
                (0, 0) => Some(t),
                _ => None,
            };
            if m >= 1 {
                let term = by_d[m] * SO3::jr_inv(&d[m]);
                acc = Some(acc.map_or(term, |a| a + term));
            }
            if m + 1 < k {
                let inv = match planted {
                    Some(Planted::JrInvForJlInv) => SO3::jr_inv(&d[m + 1]),
                    _ => SO3::jl_inv(&d[m + 1]),
                };
                let term = -(by_d[m + 1] * inv);
                acc = Some(acc.map_or(term, |a| a + term));
            }
            acc.unwrap_or_else(zero3)
        })
        .collect()
}

fn order_of(fn_id: &str) -> usize {
    match fn_id {
        f if f.contains("_vel") => 1,
        f if f.contains("_acc") => 2,
        _ => 0,
    }
}

/// §21's Gram entry `(s, t)` as one reduced fraction `num / den`.
fn gram_fraction(n: &[Vec<i64>], k: usize, r: usize, s: usize, t: usize) -> (i64, i64) {
    let falling = |p: usize| factorial(p as i64) / factorial((p - r) as i64);
    let lcm = (1..2 * k as i64).fold(1i64, |a, b| a / gcd(a, b) * b);
    let mut num = 0i64;
    for p in r..k {
        for q in r..k {
            let d = (p + q - 2 * r + 1) as i64;
            num += n[s][p] * n[t][q] * falling(p) * falling(q) * (lcm / d);
        }
    }
    let den = lcm * factorial(k as i64 - 1).pow(2);
    let g = gcd(num.abs(), den).max(1);
    (num / g, den / g)
}

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn gram<S: Real + Into<f64>>(record: &Record, planted: Option<Planted>) -> Option<Vec<f64>> {
    let k = scalars::<S>(record, "k")?.first()?.value_f64() as usize;
    let r = scalars::<S>(record, "r")?.first()?.value_f64() as usize;
    let inv_h = *scalars::<S>(record, "inv_h")?.first()?;
    let n = numerators(k);
    let mut power = inv_h;
    for _ in 1..(2 * r - 1) {
        power = power * inv_h;
    }
    let mut out = Vec::with_capacity(k * k);
    for t in 0..k {
        for s in 0..k {
            let (num, den) = gram_fraction(&n, k, r, s, t);
            let g = S::lit(num as f64) / S::lit(den as f64);
            out.push(match planted {
                Some(Planted::GramUnscaled) => g,
                _ => g * power,
            });
        }
    }
    Some(out.into_iter().map(Into::into).collect())
}

/// The answer to `fn_id` at `S`, nothing for an id it does not answer or a record it cannot read.
pub(crate) fn answer<S: Real + Into<f64>>(
    fn_id: &str,
    record: &Record,
    planted: Option<Planted>,
) -> Output {
    let Some(out) = evaluate::<S>(fn_id, record, planted) else {
        return Output::new();
    };
    out.into_iter().map(|(k, v)| (k.to_string(), v)).collect()
}

fn evaluate<S: Real + Into<f64>>(
    fn_id: &str,
    record: &Record,
    planted: Option<Planted>,
) -> Option<Vec<(&'static str, Vec<f64>)>> {
    if fn_id == "rn_spline_gram" {
        return Some(vec![("g", gram::<S>(record, planted)?)]);
    }
    let (u, inv_h) = segment::<S>(record)?;
    let order = order_of(fn_id);
    let tangent = ["omega", "alpha"];
    match fn_id {
        "so3_spline_value" | "so3_spline_vel" | "so3_spline_acc" => {
            let x = rotations::<S>(record)?;
            let w = weights(x.len(), u, inv_h, planted);
            let (v, om, al) = kinematics(&x, &w, planted);
            Some(match order {
                0 => vec![("q", quat_of(&v))],
                1 => vec![(tangent[0], scaled(dense(&om, 3), 1, inv_h))],
                _ => vec![(tangent[1], scaled(dense(&al, 3), 2, inv_h))],
            })
        }
        "so3_spline_value_jac" | "so3_spline_vel_jac" | "so3_spline_acc_jac" => {
            let x = rotations::<S>(record)?;
            let w = weights(x.len(), u, inv_h, planted);
            let mut flat = Vec::with_capacity(9 * x.len());
            for block in so3_jacobians(&x, &w, order, planted) {
                let mut buf = [S::zero(); 9];
                Jac::<S, <SO3<S> as LieGroup<S>>::Tangent>::write_dense(
                    &block,
                    &mut StridedMut::col_major(&mut buf, 3, 3),
                );
                flat.extend(buf);
            }
            Some(vec![("jac", scaled(flat, order, inv_h))])
        }
        "rn_spline_value" | "rn_spline_vel" | "rn_spline_acc" => {
            let x = points::<S>(record)?;
            let w = weights(x.len(), u, inv_h, planted);
            let (v, om, al) = kinematics(&x, &w, planted);
            Some(match order {
                0 => vec![("p", v.0 .0.map(Into::into).to_vec())],
                1 => vec![("v", scaled(dense(&om, 3), 1, inv_h))],
                _ => vec![("a", scaled(dense(&al, 3), 2, inv_h))],
            })
        }
        "rn_spline_value_jac" | "rn_spline_vel_jac" | "rn_spline_acc_jac" => {
            let k = points::<S>(record)?.len();
            let w = weights(k, u, inv_h, planted);
            let row = [&w.lam, &w.dot, &w.ddot][order];
            let next = |m: usize| if m + 1 < k { row[m + 1] } else { S::zero() };
            let ks = (0..k).map(|m| row[m] - next(m)).collect();
            Some(vec![("k", scaled(ks, order, inv_h))])
        }
        "se3_spline_value" | "se3_spline_vel" | "se3_spline_acc" => {
            let x = poses::<S>(record)?;
            let w = weights(x.len(), u, inv_h, planted);
            let (v, om, al) = kinematics(&x, &w, planted);
            Some(match order {
                0 => {
                    let (r, t) = v.parts();
                    vec![("q", quat_of(&r)), ("x", t[0].0.map(Into::into).to_vec())]
                }
                1 => vec![("twist", scaled(dense(&om, 6), 1, inv_h))],
                _ => vec![("twist", scaled(dense(&al, 6), 2, inv_h))],
            })
        }
        _ => None,
    }
}

/// The clean evaluator, or one planted defect over its ids.
pub(crate) struct Spline(pub(crate) Option<Planted>);

impl Subject for Spline {
    fn name(&self) -> &str {
        self.0.map_or("seeded:spline", Planted::name)
    }

    fn supports(&self, fn_id: &str) -> bool {
        match self.0 {
            None => IDS.contains(&fn_id),
            Some(p) => p.ids().contains(&fn_id),
        }
    }

    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        if !self.supports(fn_id) {
            return Output::new();
        }
        match precision {
            Precision::F64 => answer::<f64>(fn_id, record, self.0),
            Precision::F32 => answer::<f32>(fn_id, record, self.0),
        }
    }
}

/// The clean evaluator as a subject a run names (`--subject seeded:spline`), never a plain run's.
pub(crate) fn registered() -> Registered {
    Registered {
        version: "numerics-21".to_string(),
        version_f32: "numerics-21".to_string(),
        no_f32: None,
        subject: Box::new(Spline(None)),
        planted: true,
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)] // the weights at u = 0 and 1 are exact
mod tests {
    use super::*;

    #[test]
    fn the_integer_tables_are_section_21s_and_below_two_to_the_24() {
        // `0115`'s K = 4 cumulative rows, over 3! = 6.
        let n = numerators(4);
        let cum: Vec<Vec<i64>> = (0..4)
            .map(|j| (0..4).map(|c| (j..4).map(|s| n[s][c]).sum()).collect())
            .collect();
        assert_eq!(
            cum,
            [[6, 0, 0, 0], [5, 3, -3, 1], [1, 3, 3, -2], [0, 0, 0, 1]]
        );
        for k in 2..=6 {
            let n = numerators(k);
            for r in 1..k.min(3) {
                for s in 0..k {
                    let row: i64 = (0..k).map(|t| gram_fraction(&n, k, r, s, t).0).sum();
                    let den = gram_fraction(&n, k, r, s, s).1;
                    assert!(den < 1 << 24, "{k} {r}");
                    // Rows sum to zero: every reduced numerator over its own denominator.
                    let sum: f64 = (0..k)
                        .map(|t| {
                            let (a, b) = gram_fraction(&n, k, r, s, t);
                            assert!(a.abs() < 1 << 24 && b < 1 << 24, "{k} {r} {s} {t}");
                            a as f64 / b as f64
                        })
                        .sum();
                    assert!(sum.abs() < 1e-12, "{k} {r} {s}: {sum} {row}");
                }
            }
        }
        // `NUMERICS.md` §21's check values at h = 1.
        assert_eq!(gram_fraction(&numerators(4), 4, 1, 0, 0), (1, 20));
        assert_eq!(gram_fraction(&numerators(6), 6, 2, 0, 0), (1, 252));
    }

    #[test]
    fn the_weights_are_a_partition_and_vanish_at_the_segment_start() {
        for k in 2..=6 {
            let w = weights::<f64>(k, 0.0, 1.0, None);
            assert_eq!(w.lam[0], 1.0);
            // λ_{K−1}(0) = 0 with its first K − 2 derivatives.
            assert_eq!(w.lam[k - 1], 0.0);
            if k >= 3 {
                assert_eq!(w.dot[k - 1], 0.0);
            }
            // At u = 1, λ_1 = 1: the first control has left the segment.
            assert_eq!(weights::<f64>(k, 1.0, 1.0, None).lam[1], 1.0);
        }
    }
}
