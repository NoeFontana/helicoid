//! The coefficient kernel (`docs/PHASE3.md` §3, `docs/NUMERICS.md` §4, `0004`): the only module that
//! evaluates a cancelling coefficient.
//!
//! Each call site takes its coefficients as one group: `exp_coeffs`, `jr_coeffs`, `jr_inv_coeff`,
//! `q_coeffs`, `log_ratio`. A group is one `S::branch` on "every member is on its series arm", the
//! series by Horner in the branch variable; off it the members' exact arms run once, at the safe
//! argument, and each member selects its own. Not here: `se2_coeffs` (`docs/decisions/0015` (draft)
//! PH.3) and `gamma2_coeffs` (Phase 5).
//!
//! **Grouping is a cost and never a value.** A member takes the arm its own switches select,
//! whatever it is grouped with: the group's one branch decides only whether the shared exact
//! closure runs, and in the `!all` arm a member below its own switch still takes its Horner. So
//! `tests::groups_are_their_members_own_arm` compares each output against that member's arms alone,
//! and a regrouping is measurable without risking a bit (`0047`).
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
//!   member is on its series arm, i.e. below the group's smallest switch: at binary64 `θ² < 5.23`
//!   for `exp_coeffs`, where `cos θ/2` switches, so `Exp` is on its series arms to `θ = 2.29`;
//!   `8.35` for `jr_coeffs` (`a`), `8.98` for `q_coeffs` (`e`), `3.92` for `jr_inv_coeff`.
//!   Above it the group's exact arms run once, at `select(all, 1, θ²)`, sharing `θ` and each
//!   transcendental they can, and each member selects by its own mask; a scalar mask runs the
//!   Horner of the members on their series arm alone. `(sqrt, sin_cos, sin)` calls, as
//!   `coeffs::tests`'s counting scalar pins them: `exp_coeffs` `(1, 1, 0)` (`NUMERICS.md` §3.1),
//!   `jr_coeffs` **`(1, 0, 2)`** — no `sin_cos` at all, since both of its sines discard their
//!   cosine (`0052`) — `jr_inv_coeff` `(1, 1, 0)`, `q_coeffs` **`(1, 1, 1)`** (`b` and `e` share
//!   `sin θ` and `cos θ`; `d`'s half-angle sine is the `sin`). Evaluated for nothing: one `sin` of
//!   `jr_coeffs` for `θ² ∈ [8.35, 9.65)` and `b`'s and `d`'s Horners beside `q_coeffs`' exact arm
//!   for `[8.98, 9.65)` — and, across groups, `b` twice in `SEn3::jr` (`0047` *Further work*).
//! - `r` (`log_ratio`) takes its series iff `w > 0` and `s = n²/w² < switch`, `s` formed by a division
//!   (0014 (draft) question 29, 0015 (draft) NU.6): the reading the sweep measures. Where `s`
//!   overflows (a unit `q` with `w` below about `1e-154` at `f64`, `5e-20` at `f32`) the mask is false
//!   and the value right, but a lane that evaluates both arms sees a non-finite operation; NU.6's
//!   division-free mask has no such case.
//!
//! **What a switch row is.** Its `objective` is the maximum over the corpus's records, not over the
//! domain: between records an exact arm just above a switch can exceed it. Six of the sixteen
//! switches sit on the **domain** bound — the largest grid point below `π²` — which `generated.rs`
//! prints as `binding: domain` beside each one, and `c` at binary64 sits on the term cap
//! (`0039` items 8 and 9). A boundary optimum is reported, never silent: eight switches sat on the
//! grid's old ceiling for the file's whole existence with nothing in the file saying so.
//!
//! **A switch row has two arms.** `below` is the series/exact choice stage 1 makes; `short_below`
//! and `short_terms` are the prefix stage 2 adds under it, admitted only where it agrees with the
//! whole arm to the bit at every grid point and every corpus record (`PHASE1.md` §6, `0047`). The
//! catalogue's corpus-weighted term count is 1.87x lower at binary64 and 1.99x at binary32 for it,
//! and every conformance row is unchanged byte for byte.
//!
//! Every grouped entry point is on a default build's path: `exp_coeffs`, `jr_coeffs`,
//! `jr_inv_coeff` and `log_ratio` from `SO3`, and `q_coeffs` — Barfoot's `Q` — from `SEn3`, which
//! is `PHASE3.md` §5 and the last of the five to arrive.

use helicoid_linalg::Real;

mod generated;
mod kernel;
pub(crate) use kernel::{
    exp_coeffs, jr_coeffs, jr_inv_coeff, log_ratio, log_ratio_short, log_ratio_takes_short_arm,
    q_coeffs,
};
#[cfg(feature = "__sweep")]
pub mod sweep;
#[cfg(test)]
mod tests;

/// `(sin θ/θ, d(sin θ/θ)/d(θ²))` at `theta_sq = θ²`: with [`sinc_value`], the only items of
/// `coeffs` outside the crate (`docs/decisions/0062`).
///
/// The value is `2k cos(θ/2)` of the group `Exp` uses (`docs/NUMERICS.md` §4), so it adds no
/// switch point and no series; the derivative is the `Dual` lane of that same evaluation, the
/// derivative of the arm taken at `theta_sq`. Taking `θ²` spares a camera or bearing model the
/// `sqrt` and the hand-chosen guard at `θ = 0`: at `theta_sq = 0` it returns `(1, -1/6)`.
///
/// ```
/// let (alpha, d) = helicoid::sinc(0.25_f64); // θ = 0.5
/// assert!((alpha - 0.5_f64.sin() / 0.5).abs() < 1e-16);
/// assert!((d - (0.5 * 0.5_f64.cos() - 0.5_f64.sin()) / (2.0 * 0.125)).abs() < 1e-16);
/// ```
///
/// # Domain
///
/// Defined for every `theta_sq >= 0`. For `θ <= π`, the range the corpus id `coeff_alpha` scores,
/// the value and the derivative are each within `5u` times their own condition number in `θ²`
/// (measured: 1.9u and 2.4u at `f64`, 2.0u and 2.3u at `f32`); beyond `π` nothing is promised.
/// NaN returns NaN; a negative `theta_sq` is a sign error upstream and fails a `debug_assert!`.
pub fn sinc<S: Real>(theta_sq: S) -> (S, S) {
    kernel::alpha(theta_sq)
}

/// `sin θ/θ` at `theta_sq = θ²`, the value of [`sinc`] to the bit, without its derivative.
///
/// For a caller that needs no Jacobian, such as a camera model unprojecting every point of a
/// frame: `sinc` evaluates its group once on `Dual<S, 1>` to get the derivative, and this runs it
/// on plain `S`.
///
/// ```
/// assert_eq!(helicoid::sinc_value(0.25_f64).to_bits(), helicoid::sinc(0.25_f64).0.to_bits());
/// ```
///
/// # Domain
///
/// [`sinc`]'s.
pub fn sinc_value<S: Real>(theta_sq: S) -> S {
    kernel::alpha_value(theta_sq)
}

/// `M` series terms, in the branch variable, below the switch `below`, and the first
/// `short_terms` of them alone below `short_below`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Switch<S, const M: usize> {
    /// The series arm is taken where the branch variable is `< below`: `θ²`, and `n²/w²` for `r`.
    pub(crate) below: S,
    /// The series in the branch variable, lowest power first.
    pub(crate) series: [S; M],
    /// The **second** switch (`0039` items 5, 6 and 10): below it the first `short_terms` terms
    /// alone. A term count is set by the largest branch variable an arm serves and paid by the
    /// smallest, so one arm makes every input pay the hardest one's price; the sweep chooses the
    /// cheapest prefix that is **bit-identical** to all `M` terms everywhere below this point, so
    /// which arm runs is a cost and never a value. `0` where no prefix is cheaper, and no branch
    /// variable is below that (`z >= 0`).
    pub(crate) short_below: S,
    /// `M` where [`Self::short_below`] is `0`: then the short arm is the long one and unreachable.
    pub(crate) short_terms: usize,
}

impl<S: Copy, const M: usize> Switch<S, M> {
    /// The first `M` of the `N` swept terms (a compile error when `M > N`), and the `short_terms`
    /// of them the second arm takes below `short_below`.
    pub(crate) const fn first<const N: usize>(
        below: S,
        short_below: S,
        short_terms: usize,
        swept: &[S; N],
    ) -> Self {
        let mut series = [below; M];
        let mut i = 0;
        while i < M {
            series[i] = swept[i];
            i += 1;
        }
        Self {
            below,
            series,
            short_below,
            short_terms,
        }
    }
}
