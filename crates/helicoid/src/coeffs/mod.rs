//! The coefficient kernel (`docs/PHASE3.md` §3, `docs/NUMERICS.md` §4, `0004`): the only module that
//! evaluates a cancelling coefficient.
//!
//! Each call site takes its coefficients as one group: `exp_coeffs`, `jr_coeffs`, `jr_inv_coeff`,
//! `q_coeffs`, `log_ratio`. A group is one `S::branch` on "every member is on its series arm", the
//! series by Horner in the branch variable; off it the members' exact arms run once, at the safe
//! argument, and each member selects its own. Not here: `se2_coeffs` (`docs/decisions/0015` (draft)
//! PH.3) and `gamma2_coeffs` (Phase 5).
//!
//! **Constants** are `generated`, one [`Switch`] per coefficient per precision, chosen by
//! `S::PRECISION` (a `const` match) and converted with `S::lit`, exact because the table is that
//! precision's own (`0003` item 6). `cargo xtask thresholds` writes the file and
//! `conformance/sweeps/thresholds.csv` from `conformance/corpus/coeff_series.jsonl` and a sweep of the
//! arms below, through the hidden feature `__sweep` (`sweep`): the code measured is the code shipped
//! (`0004` item 4). The file has two parts. The swept series, the first eight terms of each row
//! rounded once at each precision, is a function of the corpus alone. The `Switch` constants are a
//! function of the sweep. The sweep reads only the first part, never a `Switch`, so a second run
//! writes the same bytes. When the corpus's series change, the first run finds the compiled-in
//! series stale, writes them beside placeholder switches and stops; the second sweeps.
//!
//! **Readings the documents leave open**, provisional:
//!
//! - A group's members each have their own switch, as `generated.rs` holds one per coefficient and
//!   the sweep measures each alone (0015 (draft) NU.5, option (b)). No exact arm runs while every
//!   member is on its series arm, i.e. below the group's smallest switch: `θ² <` `5.6e-15` (`f64`) or
//!   `2.6e-6` (`f32`) for `exp_coeffs`, where `cos θ/2` switches, so `Exp` is on its exact arm from
//!   `θ = 7.5e-8` (`f64`) or `1.6e-3` (`f32`), and `0.72` and `0.96` (`f64`; `1` at `f32`) for
//!   `jr_coeffs` and `q_coeffs`.
//!   Above it the group's exact arms run once, at `select(all, 1, θ²)`, sharing `θ` and each
//!   `sin_cos` they can, and each member selects by its own mask; a scalar mask runs the Horner of
//!   the members on their series arm alone. `sqrt` and `sin_cos` calls: `exp_coeffs` one and one
//!   (`NUMERICS.md` §3.1), `jr_coeffs` two and two, `jr_inv_coeff` one and one, `q_coeffs` one and two
//!   (`b` and `e` share `sin θ`, `cos θ`). Evaluated for nothing: one `sin_cos` of `jr_coeffs` for
//!   `θ² ∈ [0.72, 0.96)` and of `q_coeffs` for `[0.96, 1)` (`f64` only), a division of `exp_coeffs`.
//!   Unmeasured: no bench yet.
//! - `r` (`log_ratio`) takes its series iff `w > 0` and `s = n²/w² < switch`, `s` formed by a division
//!   (0014 (draft) question 29, 0015 (draft) NU.6): the reading the sweep measures. Where `s`
//!   overflows (a unit `q` with `w` below about `1e-154` at `f64`, `5e-20` at `f32`) the mask is false
//!   and the value right, but a lane that evaluates both arms sees a non-finite operation; NU.6's
//!   division-free mask has no such case.
//!
//! **What a switch row is.** Its `objective` is the maximum over the corpus's records, not over the
//! domain: between records an exact arm just above a switch can exceed it. Eight of the sixteen
//! switches (`k` and `d` at `f64`; `k`, `a`…`e` at `f32`) are `θ² = 1`, the top of the sweep's grid
//! (`grid_index` 1024), not a measured optimum (0014 (draft) question 9).
//!
//! SO(3) is the first consumer (`lib.rs`): `exp_coeffs`, `jr_coeffs`, `jr_inv_coeff` and
//! `log_ratio` are on a default build's path from `SO3`.

mod generated;
mod kernel;
pub(crate) use kernel::{exp_coeffs, jr_coeffs, jr_inv_coeff, log_ratio};
// `q_coeffs` is SE_N(3)'s alone, so only the tests and the sweep name it until §5 lands; the
// function itself stays compiled (`kernel.rs` says why), so its exact arms do not split.
#[cfg(any(test, feature = "__sweep"))]
pub(crate) use kernel::q_coeffs;
#[cfg(feature = "__sweep")]
pub mod sweep;
#[cfg(test)]
mod tests;

/// `M` series terms, in the branch variable, below the switch `below`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Switch<S, const M: usize> {
    /// The series arm is taken where the branch variable is `< below`: `θ²`, and `n²/w²` for `r`.
    pub(crate) below: S,
    /// The series in the branch variable, lowest power first.
    pub(crate) series: [S; M],
}

impl<S: Copy, const M: usize> Switch<S, M> {
    /// The first `M` of the `N` swept terms (a compile error when `M > N`).
    pub(crate) const fn first<const N: usize>(below: S, swept: &[S; N]) -> Self {
        let mut series = [below; M];
        let mut i = 0;
        while i < M {
            series[i] = swept[i];
            i += 1;
        }
        Self { below, series }
    }
}
