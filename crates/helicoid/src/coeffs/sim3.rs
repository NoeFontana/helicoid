//! The Sim(3) coefficient group (`NUMERICS.md` §9, `0070`, `docs/maths/sim3.md` SM.5–SM.11):
//! `v₀, v₁, v₂` of `V`, `w₀, w₁, w₂` of `Γ₂(M)` and `∂_z v₁, ∂_z v₂` of `Q`, as functions of
//! `(σ, z = θ²)`, one call per site.
//!
//! Two arms, chosen by one mask on `Δ = σ² + z`. Below `r_×²` the **total-degree series** of SM.6(b):
//! the recurrences of SM.3 and SM.11(c) at per-sequence depths, so every scalar is truncated at
//! total degree `< K` and the tail is a function of `r_Δ` alone; it covers `Δ = 0`, where it is the
//! limits, without a division. Above, the **closed forms** of SM.5(b), SM.11(b)–(c) on `a, b, d`
//! and `sin θ/θ` from the kernel's swept arms, `v₁` by form b from `Δ_b`, `w₀ = g₂(σ)` and
//! `v₀ = g₁(σ)` by their series in `σ` below `σ_w`, and for `σ > 0` the **scaled arrangement** of SM.7(g): every scalar is formed
//! as `e^{-σ} c` from `expm1(-σ)` and `e^{-σ}`, then multiplied by `e^σ`, so nothing overflows
//! before the answer does. `K`, `r_×²`, `Δ_b`, `σ_w²` and `g₂`'s terms are generated
//! (`generated_sim3`, `0004`).

use helicoid_linalg::{Blend, Mask, Precision, Real};

use super::generated_sim3::{RECIP_FACT_F32, RECIP_FACT_F64, SIM3_F32, SIM3_F64};
use super::kernel::sinc_a_b_d;

/// The generated parameters of the group at one precision (`0070` decision 4).
#[derive(Clone, Copy)]
pub(crate) struct Sim3Switch<T> {
    /// `r_×²`: the series arm below it, in `Δ = σ² + z`.
    pub(crate) series_below: T,
    /// `K`: the series' total degree is `< K`.
    pub(crate) series_terms: usize,
    /// `Δ_b`: `v₁` by form a below it, form b at and above.
    pub(crate) form_b_below: T,
    /// `σ_w²`: `w₀ = g₂(σ)` by its series below it, in `σ²`.
    pub(crate) w0_below: T,
    /// The terms of `g₂`'s series, `Σ σⁱ/(i + 2)!`.
    pub(crate) w0_terms: usize,
}

/// The eight scalars: `v = [v₀, v₁, v₂]`, `w = [w₀, w₁, w₂]`, `dv = [∂_z v₁, ∂_z v₂]`; and
/// `chi_delta = χΔ`, `χ = v₀ − zv₂` (SM.5(c)), `Δ = σ² + z`, which the closed arm forms without
/// that cancellation. Its reader divides by `Δ` where it reads it (`0096`, draft).
#[derive(Clone, Copy, Debug)]
pub(crate) struct Sim3Coeffs<S> {
    pub(crate) v: [S; 3],
    pub(crate) w: [S; 3],
    pub(crate) dv: [S; 2],
    pub(crate) chi_delta: S,
}

impl<S: Real> Blend<S> for Sim3Coeffs<S> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self {
            v: <[S; 3]>::blend(m, t.v, f.v),
            w: <[S; 3]>::blend(m, t.w, f.w),
            dv: <[S; 2]>::blend(m, t.dv, f.dv),
            chi_delta: S::blend(m, t.chi_delta, f.chi_delta),
        }
    }
}

#[inline]
fn switch<S: Real>() -> Sim3Switch<f64> {
    match S::PRECISION {
        Precision::F64 => SIM3_F64,
        Precision::F32 => {
            let s = SIM3_F32;
            Sim3Switch {
                series_below: f64::from(s.series_below),
                series_terms: s.series_terms,
                form_b_below: f64::from(s.form_b_below),
                w0_below: f64::from(s.w0_below),
                w0_terms: s.w0_terms,
            }
        }
    }
}

/// `1/n!`, correctly rounded at `S`'s precision (`sim3_series.jsonl`); `NaN` past the table.
#[inline]
fn recip_fact<S: Real>(n: usize) -> S {
    let x = match S::PRECISION {
        Precision::F64 => RECIP_FACT_F64.get(n).copied(),
        Precision::F32 => RECIP_FACT_F32.get(n).map(|&x| f64::from(x)),
    };
    x.map_or_else(|| S::zero() / S::zero(), S::lit)
}

/// The total-degree series at depth `k` (SM.6(b), SM.11(c)). `FULL` adds `w` and `dv`; without it
/// the recurrence stops at what `v` needs, and `v` is the same bits either way.
///
/// Level `j` holds `p = σʲ`, `ξ_j`, `ζ_j`, `ξ'_j`, `ζ'_j`, of total degrees `j`, `j - 1`, `j - 2`,
/// `j - 3`, `j - 4`, so `v₀, w₀` sum `j < k`, `v₁, w₁` `j ≤ k`, `v₂, w₂` `j ≤ k + 1`, `∂_z v₁`
/// `j ≤ k + 2` and `∂_z v₂` `j ≤ k + 3`: every scalar to total degree `k - 1` (SM.7(d)(ii)).
#[inline]
pub(crate) fn series<S: Real, const FULL: bool>(sigma: S, z: S, k: usize) -> Sim3Coeffs<S> {
    let [c, _] = series_at::<S, FULL, false>(sigma, z, k);
    c
}

/// [`series`] at `σ` and, with `PAIR`, at `-σ`, each to the bit, from one recurrence.
///
/// A level of total degree `n` is odd or even in `σ` with `n`, and rounding commutes with
/// negation, so the recurrence at `-σ` holds the level at `σ` times `(-1)ⁿ` exactly: only the sums
/// differ, and `-σ`'s adds a term where `σ`'s does or subtracts it (`0073`).
#[inline]
fn series_at<S: Real, const FULL: bool, const PAIR: bool>(
    sigma: S,
    z: S,
    k: usize,
) -> [Sim3Coeffs<S>; 2] {
    let zero = S::zero();
    let (mut p, mut xi, mut zeta, mut dxi, mut dzeta) = (S::one(), zero, zero, zero, zero);
    // `v₀, v₁, v₂, w₀, w₁, w₂, ∂_z v₁, ∂_z v₂`, each at `σ` and `-σ`.
    let mut acc = [[zero; 2]; 8];
    let last = if FULL { k + 3 } else { k + 1 };
    for j in 0..=last {
        let (r1, r2) = (recip_fact::<S>(j + 1), recip_fact::<S>(j + 2));
        // `p`, `ζ`, `ζ'` have the parity of `j`; `ξ`, `ξ'` the other.
        let odd = j % 2 == 1;
        if j < k {
            gain::<S, PAIR>(&mut acc[0], p * r1, odd);
            if FULL {
                gain::<S, PAIR>(&mut acc[3], p * r2, odd);
            }
        }
        if j <= k {
            gain::<S, PAIR>(&mut acc[1], xi * r1, !odd);
            if FULL {
                gain::<S, PAIR>(&mut acc[4], xi * r2, !odd);
            }
        }
        if j <= k + 1 {
            gain::<S, PAIR>(&mut acc[2], zeta * r1, odd);
            if FULL {
                gain::<S, PAIR>(&mut acc[5], zeta * r2, odd);
            }
        }
        if FULL {
            if j <= k + 2 {
                gain::<S, PAIR>(&mut acc[6], dxi * r1, !odd);
            }
            gain::<S, PAIR>(&mut acc[7], dzeta * r1, odd);
        }
        // SM.3, SM.11(c): ξ_{j+1} = σξ + σʲ - zζ, ζ_{j+1} = σζ + ξ,
        // ξ'_{j+1} = σξ' - ζ - zζ', ζ'_{j+1} = σζ' + ξ'.
        let (nxi, nzeta) = ((sigma * xi + p) - z * zeta, sigma * zeta + xi);
        if FULL {
            (dxi, dzeta) = ((sigma * dxi - zeta) - z * dzeta, sigma * dzeta + dxi);
        }
        (xi, zeta, p) = (nxi, nzeta, p * sigma);
    }
    // Below `r_×²` nothing cancels in `v₀ − zv₂`: `v₀ ≈ 1` and `zv₂ < r_×²/6`. No reader of
    // `χΔ` takes this arm (`b`'s exact arm starts above `r_×²`); it is the definition.
    let at = |i: usize| Sim3Coeffs {
        v: [acc[0][i], acc[1][i], acc[2][i]],
        w: [acc[3][i], acc[4][i], acc[5][i]],
        dv: [acc[6][i], acc[7][i]],
        chi_delta: (acc[0][i] - z * acc[2][i]) * (sigma * sigma + z),
    };
    [at(0), at(1)]
}

/// `σ`'s sum gains `t`; with `PAIR`, `-σ`'s gains `-t` where `flip` and `t` elsewhere.
#[inline]
fn gain<S: Real, const PAIR: bool>(acc: &mut [S; 2], t: S, flip: bool) {
    acc[0] = acc[0] + t;
    if PAIR {
        acc[1] = if flip { acc[1] - t } else { acc[1] + t };
    }
}

/// `g_ℓ(σ) = Σ σⁱ/(i + ℓ)!` by Horner from the last of `terms` terms.
#[inline]
pub(crate) fn g_series<S: Real>(sigma: S, l: usize, terms: usize) -> S {
    (0..terms)
        .rev()
        .fold(S::zero(), |acc, i| recip_fact::<S>(i + l) + sigma * acc)
}

/// The closed forms at a safe argument (`Δ > 0`) and **scaled** for `σ > 0`: every scalar is
/// `e^{-σ} c` there, and [`unscale`] forms `c`; `form_b` is `Δ ≥ Δ_b`, `w0_small` is `σ² < σ_w²`.
#[inline]
pub(crate) fn closed<S: Real, const FULL: bool>(
    sigma: S,
    z: S,
    form_b: S::Mask,
    w0_small: S::Mask,
    w0_terms: usize,
) -> Sim3Coeffs<S> {
    closed_from::<S, FULL>(
        sigma,
        z,
        Transcendentals::at(sigma, z),
        form_b,
        w0_small,
        w0_terms,
    )
}

/// What [`closed`] evaluates through `libm`: `sin θ/θ, a, b, d` at `z`, and `e^t`, `expm1 t` at
/// `t = -|σ|`, which is even in `σ` but at `σ = ±0`, where `t` is `σ` itself.
#[derive(Clone, Copy)]
struct Transcendentals<S> {
    sinc_a_b_d: (S, S, S, S),
    et: S,
    mt: S,
}

impl<S: Real> Transcendentals<S> {
    #[inline]
    fn at(sigma: S, z: S) -> Self {
        // `t = -|σ|`: `e^t` and `expm1 t` are `E` and `m₁` for `σ ≤ 0`; for `σ > 0` the scaled
        // arrangement has `E = 1`, `m₁ = -expm1(-σ)` and `e^σ` applied last (SM.7(g)).
        let t = S::select(S::zero().lt(sigma), -sigma, sigma);
        Self {
            sinc_a_b_d: sinc_a_b_d(z),
            et: t.exp(),
            mt: t.exp_m1(),
        }
    }

    /// [`at`](Self::at)`(-σ, z)` to the bit, given `at(σ, z)`. At `σ = ±0` the two `t` are `σ` and
    /// `-σ`, which differ in sign and, on `Dual`, in derivative, so they are evaluated again: a
    /// branch no other `σ` takes.
    #[inline]
    fn mirrored(self, sigma: S) -> Self {
        let zero = S::zero();
        let at_zero = sigma.le(zero).and(zero.le(sigma));
        let (et, mt) = S::branch(
            at_zero,
            || ((-sigma).exp(), (-sigma).exp_m1()),
            || (self.et, self.mt),
        );
        Self { et, mt, ..self }
    }
}

/// [`closed`] from its transcendentals.
#[inline]
fn closed_from<S: Real, const FULL: bool>(
    sigma: S,
    z: S,
    tr: Transcendentals<S>,
    form_b: S::Mask,
    w0_small: S::Mask,
    w0_terms: usize,
) -> Sim3Coeffs<S> {
    let (zero, one) = (S::zero(), S::one());
    let (sinc, a, b, d) = tr.sinc_a_b_d;
    let delta = sigma * sigma + z;
    let pos = zero.lt(sigma);
    let at_zero = sigma.le(zero).and(zero.le(sigma));
    let (et, mt) = (tr.et, tr.mt);
    let e = S::select(pos, one, et);
    let m1 = S::select(pos, -mt, mt);
    let f = S::select(pos, et, one);
    let ss = S::select(at_zero, one, sigma);
    // Below `σ_w`, `v₀ = g₁(σ)` by its series too (SM.6(a)), with one term more than `w₀`'s, so
    // its tail is at most `|σ|` times the one the sweep measured for `w₀`. The value `expm1 σ/σ`
    // is accurate there, but its `Dual` derivative `(σe^σ − expm1 σ)/σ²` cancels to `u/σ`
    // (SM.7(f)): `jacobians_match_dual` read 6720 `u` at `σ = 9e-7` before (`0070`).
    let v0 = S::select(
        w0_small,
        f * g_series(sigma, 1, w0_terms + 1),
        S::select(at_zero, one, m1 / ss),
    );
    let v2 = ((v0 + sigma * e * a) - e * sinc) / delta;
    let v1a = e * a - sigma * v2;
    let v1b = ((sigma * e * sinc - m1) + e * z * a) / delta;
    let v1 = S::select(form_b, v1b, v1a);
    // `χΔ = Δ Re g₁(σ + iθ) = σ(e^σ cos θ − 1) + θe^σ sin θ` (SM.5(c)), with `e^σ cos θ − 1 =
    // m₁ − eza` and `σm₁ = σ²v₀`, scaled as the rest: the two non-negative terms first, then the
    // one of `−σ`'s sign. `v₀ − zv₂` cancels to `sin θ/θ` at `σ = 0` near `π` (`0096`, draft).
    let chi_delta = (sigma * (sigma * v0) + e * (z * sinc)) - sigma * (e * (z * a));
    let mut c = Sim3Coeffs {
        v: [v0, v1, v2],
        w: [zero; 3],
        dv: [zero; 2],
        chi_delta,
    };
    if FULL {
        let w0 = S::select(
            w0_small,
            f * g_series(sigma, 2, w0_terms),
            (m1 - sigma * f) / (ss * ss),
        );
        let w2 = ((w0 + sigma * v2) - v1) / delta;
        let w1 = v2 - sigma * w2;
        let half = S::lit(0.5);
        let (da, dsinc) = (-(half * (b - S::lit(2.0) * d)), -(half * (a - b)));
        let dv2 = ((sigma * e * da - e * dsinc) - v2) / delta;
        let dv1 = e * da - sigma * dv2;
        c.w = [w0, w1, w2];
        c.dv = [dv1, dv2];
    }
    c
}

/// `c = e^σ ĉ` where `scaled`, `ĉ` itself elsewhere: the last step of the scaled arrangement.
#[inline]
pub(crate) fn unscale<S: Real>(c: Sim3Coeffs<S>, sigma: S, scaled: S::Mask) -> Sim3Coeffs<S> {
    let scale = S::branch(scaled, || sigma.exp(), S::one);
    Sim3Coeffs {
        v: c.v.map(|x| scale * x),
        w: c.w.map(|x| scale * x),
        dv: c.dv.map(|x| scale * x),
        chi_delta: scale * c.chi_delta,
    }
}

/// The group at the generated parameters: the series below `r_×²`, else the closed forms at the
/// safe argument `(1, 1)` for a lane that is not on them (`0003` item 3).
#[inline]
fn group<S: Real, const FULL: bool>(sigma: S, z: S) -> (Sim3Coeffs<S>, S::Mask) {
    let sw = switch::<S>();
    let delta = sigma * sigma + z;
    let small = delta.lt(S::lit(sw.series_below));
    let scaled = small.not().and(S::zero().lt(sigma));
    let c = S::branch(
        small,
        || series::<S, FULL>(sigma, z, sw.series_terms),
        || {
            let one = S::one();
            let (s, zz) = (S::select(small, one, sigma), S::select(small, one, z));
            let d = s * s + zz;
            closed::<S, FULL>(
                s,
                zz,
                S::lit(sw.form_b_below).le(d),
                (s * s).lt(S::lit(sw.w0_below)),
                sw.w0_terms,
            )
        },
    );
    (c, scaled)
}

/// [`group`] at `σ` and at `-σ`, each to the bit: `Δ`, so the branch, the safe argument's
/// masks and the closed arm's [`Transcendentals`] are even in `σ` and formed once, and the series'
/// recurrence runs once ([`series_at`]).
#[inline]
fn group_pair<S: Real>(sigma: S, z: S) -> [(Sim3Coeffs<S>, S::Mask); 2] {
    let sw = switch::<S>();
    let delta = sigma * sigma + z;
    let small = delta.lt(S::lit(sw.series_below));
    let zero = S::zero();
    let c = S::branch(
        small,
        || series_at::<S, true, true>(sigma, z, sw.series_terms),
        || {
            let one = S::one();
            let (s, zz) = (S::select(small, one, sigma), S::select(small, one, z));
            let d = s * s + zz;
            let (form_b, w0_small) = (
                S::lit(sw.form_b_below).le(d),
                (s * s).lt(S::lit(sw.w0_below)),
            );
            let tr = Transcendentals::at(s, zz);
            [
                closed_from::<S, true>(s, zz, tr, form_b, w0_small, sw.w0_terms),
                closed_from::<S, true>(-s, zz, tr.mirrored(s), form_b, w0_small, sw.w0_terms),
            ]
        },
    );
    let big = small.not();
    [
        (c[0], big.and(zero.lt(sigma))),
        (c[1], big.and(zero.lt(-sigma))),
    ]
}

/// All eight scalars at `(σ, z = θ²)`.
///
/// # Domain
///
/// `z ≥ 0`; `|σ| < ln MAX` for a finite answer (`NUMERICS.md` §12).
#[inline]
pub(crate) fn sim3_coeffs<S: Real>(sigma: S, z: S) -> Sim3Coeffs<S> {
    let (c, scaled) = group::<S, true>(sigma, z);
    unscale(c, sigma, scaled)
}

/// `(ĉ, scaled)`: the eight scalars as the arm forms them, `ĉ = e^{-σ} c` where `scaled` (the closed
/// arm at `σ > 0`) and `c` elsewhere. For a caller that would square or invert `c` and overflow
/// first: `Log`'s `V⁻¹ = e^{-σ} V̂⁻¹` (SM.7(g)).
#[inline]
pub(crate) fn sim3_coeffs_scaled<S: Real>(sigma: S, z: S) -> (Sim3Coeffs<S>, S::Mask) {
    group::<S, true>(sigma, z)
}

/// [`sim3_coeffs_scaled`] at `σ` and at `-σ`, each to the bit, for less than two calls cost
/// ([`group_pair`]): `⊖`'s Jacobians want `J_l⁻¹` and `J_r⁻¹` at one `τ`.
#[inline]
pub(crate) fn sim3_coeffs_scaled_pair<S: Real>(sigma: S, z: S) -> [(Sim3Coeffs<S>, S::Mask); 2] {
    group_pair(sigma, z)
}

/// [`sim3_coeffs_scaled`]'s `v` alone: the same bits.
#[inline]
pub(crate) fn sim3_v_scaled<S: Real>(sigma: S, z: S) -> ([S; 3], S::Mask) {
    let (c, scaled) = group::<S, false>(sigma, z);
    (c.v, scaled)
}

/// `([v₀, v₁, v₂], χΔ)` alone, for `Exp`: [`sim3_coeffs`]'s `v` and `chi_delta` to the bit.
///
/// # Domain
///
/// [`sim3_coeffs`]'s.
#[inline]
pub(crate) fn sim3_v<S: Real>(sigma: S, z: S) -> ([S; 3], S) {
    let (c, scaled) = group::<S, false>(sigma, z);
    let c = unscale(c, sigma, scaled);
    (c.v, c.chi_delta)
}
