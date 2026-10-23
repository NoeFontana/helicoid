//! The hidden `__sweep` surface (`0004` item 4, `docs/PHASE3.md` §3), for `cargo xtask thresholds`
//! and nothing else: the arms it measures, the groups it scores against its own choice, and the
//! swept series it checks against the corpus. Not API; only `xtask` enables the feature.
//!
//! An exact arm has no safe argument: the sweep measures it at the record's own `z`, and it is
//! not finite at `0`, as `docs/maths/coefficients.md` CO.15 says.

use helicoid_linalg::Real;

use super::generated::{
    SWEPT_A_F32, SWEPT_A_F64, SWEPT_B_F32, SWEPT_B_F64, SWEPT_COS_HALF_F32, SWEPT_COS_HALF_F64,
    SWEPT_C_F32, SWEPT_C_F64, SWEPT_D_F32, SWEPT_D_F64, SWEPT_E_F32, SWEPT_E_F64, SWEPT_K_F32,
    SWEPT_K_F64, SWEPT_R_F32, SWEPT_R_F64,
};
use super::kernel;

macro_rules! arms {
    ($($exact:ident, $series:ident;)*) => {$(
        /// The exact arm at `z = θ² > 0`.
        pub fn $exact<S: Real>(z: S) -> S {
            kernel::$exact(z)
        }

        /// The first `terms` (1 to 8) series terms at `z = θ²`.
        pub fn $series<S: Real>(z: S, terms: usize) -> S {
            kernel::$series(z, terms)
        }
    )*};
}

arms! {
    exact_k, series_k;
    exact_a, series_a;
    exact_b, series_b;
    exact_c, series_c;
    exact_d, series_d;
    exact_e, series_e;
    exact_cos_half, series_cos_half;
}

/// The exact arm of `r` at `n² = n2 > 0`.
pub fn exact_r<S: Real>(n2: S, w: S) -> S {
    kernel::exact_r(n2, w)
}

/// The first `terms` (1 to 8) series terms of `r` at `w > 0`: `2/w` times the series in `n²/w²`.
pub fn series_r<S: Real>(n2: S, w: S, terms: usize) -> S {
    kernel::series_r(n2, w, terms)
}

/// `(k, cos(θ/2))` as the kernel evaluates them at `θ² = z`.
pub fn exp_coeffs<S: Real>(z: S) -> (S, S) {
    super::exp_coeffs(z)
}

/// `(a, b)` as the kernel evaluates them at `θ² = z`.
pub fn jr_coeffs<S: Real>(z: S) -> (S, S) {
    super::jr_coeffs(z)
}

/// `c` as the kernel evaluates it at `θ² = z`.
pub fn jr_inv_coeff<S: Real>(z: S) -> S {
    super::jr_inv_coeff(z)
}

/// `(b, d, e)` as the kernel evaluates them at `θ² = z`.
pub fn q_coeffs<S: Real>(z: S) -> (S, S, S) {
    super::q_coeffs(z)
}

/// `(b, d)` as the kernel evaluates them at `θ² = z`.
pub fn gamma2_coeffs<S: Real>(z: S) -> (S, S) {
    super::gamma2_coeffs(z)
}

/// `r` as the kernel evaluates it at `n² = n2` and `w`.
pub fn log_ratio<S: Real>(n2: S, w: S) -> S {
    super::log_ratio(n2, w)
}

/// Whether `log_ratio_scale` takes its exact arm at `n² = n2` and `w`, its own dispatch: for
/// `helicoid:r-times-u`, the arrangement `0106` (draft) replaces.
pub fn log_ratio_takes_exact_arm<S: Real>(n2: S, w: S) -> S::Mask {
    kernel::log_ratio_takes_exact_arm(n2, w)
}

/// The swept series of `coeff` (`k`, `a`, `b`, `c`, `d`, `e`, `cos_half` or `r`) at binary64, lowest
/// power first: what the sweep measures, for the tool to check against the corpus.
pub fn swept_f64(coeff: &str) -> Option<&'static [f64]> {
    Some(match coeff {
        "k" => &SWEPT_K_F64,
        "a" => &SWEPT_A_F64,
        "b" => &SWEPT_B_F64,
        "c" => &SWEPT_C_F64,
        "d" => &SWEPT_D_F64,
        "e" => &SWEPT_E_F64,
        "cos_half" => &SWEPT_COS_HALF_F64,
        "r" => &SWEPT_R_F64,
        _ => return None,
    })
}

/// [`swept_f64`] at binary32.
pub fn swept_f32(coeff: &str) -> Option<&'static [f32]> {
    Some(match coeff {
        "k" => &SWEPT_K_F32,
        "a" => &SWEPT_A_F32,
        "b" => &SWEPT_B_F32,
        "c" => &SWEPT_C_F32,
        "d" => &SWEPT_D_F32,
        "e" => &SWEPT_E_F32,
        "cos_half" => &SWEPT_COS_HALF_F32,
        "r" => &SWEPT_R_F32,
        _ => return None,
    })
}

/// The Sim(3) group's arms (`0070`), as `[v₀, v₁, v₂, w₀, w₁, w₂, ∂_z v₁, ∂_z v₂]`: the total-degree
/// series at depth `k`.
pub fn sim3_series<S: Real>(sigma: S, z: S, k: usize) -> [S; 8] {
    flat(super::sim3::series::<S, true>(sigma, z, k))
}

/// The closed arm with its three choices handed in: `v₁` by form b, `w₀` by `w0_terms` terms of its
/// series (`0` for the exact form). `Δ > 0`.
pub fn sim3_closed<S: Real<Mask = bool>>(sigma: S, z: S, form_b: bool, w0_terms: usize) -> [S; 8] {
    let c = super::sim3::closed::<S, true>(sigma, z, form_b, w0_terms > 0, w0_terms);
    flat(super::sim3::unscale(c, sigma, S::zero().lt(sigma)))
}

/// The group at its generated parameters.
pub fn sim3_coeffs<S: Real>(sigma: S, z: S) -> [S; 8] {
    flat(super::sim3::sim3_coeffs(sigma, z))
}

/// `1/n!` as the kernel holds it, at binary64 and binary32, for the tool to check against
/// `sim3_series.jsonl`.
pub fn sim3_recip_fact() -> (&'static [f64], &'static [f32]) {
    (
        &super::generated_sim3::RECIP_FACT_F64,
        &super::generated_sim3::RECIP_FACT_F32,
    )
}

fn flat<S: Real>(c: super::sim3::Sim3Coeffs<S>) -> [S; 8] {
    [
        c.v[0], c.v[1], c.v[2], c.w[0], c.w[1], c.w[2], c.dv[0], c.dv[1],
    ]
}
