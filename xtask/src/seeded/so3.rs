//! The seeded SO(3) subject: `Exp` and `Log` of `docs/NUMERICS.md` §3.1 and §3.2 over the
//! coefficient kernels, and the two planted `Log` defects of `docs/PHASE1.md` §10. Generic over
//! `Real`, as the coefficient kernels are; the harness runs them at `f64`. Tangents are `φ`,
//! quaternions `[w, x, y, z]`.
//!
//! Readings where §3.1, §3.2 and §4 are silent, each the smallest:
//!
//! - **`cos(θ/2)`** has a committed series (`coeff_series.jsonl`) this subject does not read, and
//!   no cancellation (`docs/maths/coefficients.md` CO.10), so `Exp` evaluates it at every `z` by
//!   the exact arm, on its own `sin_cos` and outside `k`'s branch: two transcendentals in the
//!   exact arm where §3.1 says one, and `k` alone where §4's group rule says `(k, cos θ/2)` in
//!   one `branch`. That is 0014 (draft) question 19. `θ = sqrt(z)` is finite at `z = 0` and its `Dual` derivative there
//!   is not, so `Exp` is a value subject: the harness scores no derivative of it.
//! - **`r`** has no swept switch (`docs/PHASE1.md` §0.0: its branch variable and `w ≤ 0` domain are
//!   open), and no constant is typed for it (0014 (draft) question 20). Its series arm is taken
//!   where `n² = 0`, the one place the exact arm is `0/0`, and is the limit `2/w` of the
//!   definition, its leading term. Elsewhere the exact arm answers, even where `n²` is subnormal:
//!   `atan2` and the division share the one `n`, so `Log` stays scale-invariant (§3.2) to the bits
//!   a subnormal `n²` keeps.
//! - **The flip** is a multiplication by `copysign(1, w)`, exact, so `w = -0` flips and
//!   `w = +0` does not, as §3.2 says of `copysign`.
//!
//! The corpus cannot pin `Exp`'s overall sign or the flip at `w = ±0`: the metric aligns the
//! quaternion's sign (§11). Unit tests do.

use helicoid_linalg::Real;

use super::kernel::{coefficient, Candidate, Coeff};

/// `x·x + y·y + z·z`, left to right: `θ²` and `n²` are dot products (`NUMERICS.md` §2.1).
pub(super) fn norm_sq<S: Real>([x, y, z]: [S; 3]) -> S {
    (x * x + y * y) + z * z
}

/// Where `Log` takes the series arm of `r`: `n² = 0`, so `2 atan2(n, w)/n` is not `0/0`.
fn series_arm<S: Real>(n2: S) -> S::Mask {
    n2.le(S::zero())
}

/// Whether `Log` answers `q` by the series arm of `r`. There `2/w` is analytic in the sign of `w`,
/// so a missing flip changes nothing (`conformance::selftest_so3`).
pub(crate) fn takes_series_arm(q: [f64; 4]) -> bool {
    series_arm(norm_sq([q[1], q[2], q[3]]))
}

/// `Exp(φ) = (cos(θ/2), k(θ) φ)` (§3.1), `k` from the generated switch `cand`/`k_series`.
pub(crate) fn exp<S: Real>(phi: [S; 3], cand: Candidate<S>, k_series: &[S]) -> [S; 4] {
    let z = norm_sq(phi);
    let k = coefficient(Coeff::K, z, cand, k_series);
    let cos_half = (S::lit(0.5) * z.sqrt()).sin_cos().1;
    let [x, y, z] = phi;
    [cos_half, k * x, k * y, k * z]
}

/// `Log(q)` (§3.2): the flip, then `φ = r(n², w) u`, `r = 2 atan2(n, w)/n`.
pub(crate) fn log<S: Real>(q: [S; 4]) -> [S; 3] {
    let flip = S::one().copysign(q[0]);
    log_ratio(q.map(|c| flip * c))
}

/// Defect: `Log` without the `w < 0` flip (`docs/PHASE1.md` §10).
pub(crate) fn log_no_flip<S: Real>(q: [S; 4]) -> [S; 3] {
    log_ratio(q)
}

/// Steps 2 and 3 of §3.2 at the safe arguments (0003 item 3): `n²` in the exact arm, `w` in the
/// series arm, so no unselected arm is non-finite in any lane.
fn log_ratio<S: Real>([w, x, y, z]: [S; 4]) -> [S; 3] {
    let n2 = norm_sq([x, y, z]);
    let small = series_arm(n2);
    let two = S::lit(2.0);
    let r = S::branch(
        small,
        || two / S::select(small, w, S::one()),
        || {
            let n = S::select(small, S::one(), n2).sqrt();
            two * n.atan2(w) / n
        },
    );
    [r * x, r * y, r * z]
}

/// Defect: `Log` through `acos` of the trace (`docs/PHASE1.md` §10), `f64` only. `θ` is
/// `acos((tr R − 1)/2)` with `R`'s diagonal `1 − 2(a² + b²)` evaluated from the double `q`; the
/// axis is the vector part's, so the loss is the `acos`'s and not also the matrix's skew part
/// (`docs/maths/so3.md` SO.6, SO.12 *Sign*). The flip keeps `q` and `-q` the same rotation.
pub(crate) fn log_acos(q: [f64; 4]) -> [f64; 3] {
    let [w, x, y, z] = q;
    let flip = libm::copysign(1.0, w);
    let diagonal = |a: f64, b: f64| 1.0 - 2.0 * (a * a + b * b);
    let trace = (diagonal(y, z) + diagonal(x, z)) + diagonal(x, y);
    let theta = libm::acos(((trace - 1.0) / 2.0).clamp(-1.0, 1.0));
    let n = libm::sqrt(norm_sq([x, y, z]));
    let scale = if n > 0.0 { flip * theta / n } else { 0.0 };
    [scale * x, scale * y, scale * z]
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::conformance::testkit::splitmix;
    use crate::seeded::generated::K_F64;
    use crate::seeded::kernel::tests::{Lane, NONFINITE};

    fn k() -> (Candidate<f64>, &'static [f64]) {
        let (below, series) = K_F64.parts();
        let cand = Candidate {
            terms: series.len(),
            switch_z: below,
        };
        (cand, series)
    }

    fn exp64(phi: [f64; 3]) -> [f64; 4] {
        let (cand, series) = k();
        exp(phi, cand, series)
    }

    /// A unit quaternion of angle `theta` about the axis `(1, 2, 2)/3`, from the standard library.
    fn unit(theta: f64) -> [f64; 4] {
        let (s, c) = (theta / 2.0).sin_cos();
        [c, s / 3.0, 2.0 * s / 3.0, 2.0 * s / 3.0]
    }

    fn close(a: &[f64], b: &[f64], tol: f64) {
        let err = a
            .iter()
            .zip(b)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max);
        assert!(err <= tol, "{a:?} vs {b:?}");
    }

    #[test]
    fn exp_is_the_half_angle_quaternion_and_log_inverts_it() {
        assert_eq!(exp64([0.0; 3]), [1.0, 0.0, 0.0, 0.0]);
        for theta in [1e-9, 0.3, 1.2, 3.0] {
            let phi = [theta / 3.0, 2.0 * theta / 3.0, 2.0 * theta / 3.0];
            let q = exp64(phi);
            close(&q, &unit(theta), 4e-16);
            close(&log(q), &phi, 8e-16 * theta);
        }
    }

    #[test]
    fn log_is_a_function_of_the_rotation_except_at_w_plus_zero() {
        let mut state = 7;
        for _ in 0..64 {
            let mut q = [0.0; 4];
            for c in &mut q {
                *c = (splitmix(&mut state) >> 11) as f64 / (1u64 << 53) as f64 - 0.5;
            }
            let neg = q.map(|c| -c);
            let bits = |p: [f64; 3]| p.map(f64::to_bits);
            assert_eq!(bits(log(q)), bits(log(neg)));
            assert_ne!(bits(log_no_flip(q)), bits(log_no_flip(neg)));
        }
        // §3.2: at w = +0 nothing flips, so q and -q return +-pi n; at w = -0 `copysign` flips.
        let (u, pi) = ([0.6, 0.0, 0.8], std::f64::consts::PI);
        close(&log([0.0, u[0], u[1], u[2]]), &u.map(|c| pi * c), 4e-16);
        close(&log([0.0, -u[0], -u[1], -u[2]]), &u.map(|c| -pi * c), 4e-16);
        close(&log([-0.0, -u[0], -u[1], -u[2]]), &u.map(|c| pi * c), 4e-16);
    }

    #[test]
    fn the_series_arm_is_taken_where_n_squared_is_zero_and_nowhere_else() {
        let tiny = 3e-310;
        // n² underflows to 0: `2/w`, the analytic continuation, so the missing flip does not show.
        assert!(takes_series_arm([1.0, tiny, 0.0, -tiny]));
        assert!(takes_series_arm([-1.0, 0.0, -0.0, 0.0]));
        assert_eq!(log([1.0, tiny, 0.0, -tiny]), [2.0 * tiny, 0.0, -2.0 * tiny]);
        let neg = [-1.0, -tiny, -0.0, tiny];
        assert_eq!(
            log(neg).map(f64::to_bits),
            log([1.0, tiny, 0.0, -tiny]).map(f64::to_bits)
        );
        assert_eq!(
            log_no_flip(neg).map(f64::to_bits),
            log(neg).map(f64::to_bits)
        );
        assert_eq!(log([1.0, 0.0, 0.0, 0.0]), [0.0; 3]);
        // The arm is `2/w`, and not `2` or `2w`, which agree with it for a unit `q` (`w = ±1` there):
        // a non-unit `q` is outside `NUMERICS.md` §12's domain and tells them apart, `Log` being scale-invariant.
        assert_eq!(log([0.5, tiny, 0.0, -tiny]), [4.0 * tiny, 0.0, -4.0 * tiny]);
        assert_eq!(
            log([-0.5, tiny, 0.0, -tiny]),
            [-4.0 * tiny, 0.0, 4.0 * tiny]
        );
        // The last normal `n²` and a subnormal one, positive: the exact arm, at its full accuracy,
        // `atan2` and the division sharing the one `n`.
        for n in [1.6e-154, 1e-160] {
            assert!(!takes_series_arm([1.0, n, 0.0, 0.0]), "{n}");
            close(&log([1.0, n, 0.0, 0.0]), &[2.0 * n, 0.0, 0.0], 1e-15 * n);
        }
        // Scale invariance across the smallest normal `n²`, where no bound is typed: the rotation of
        // angle 1 about x, scaled to `1.4e-154`, is still `φ = (1, 0, 0)`.
        let (s, c) = 0.5_f64.sin_cos();
        let lambda = 1.4e-154;
        close(
            &log([c * lambda, s * lambda, 0.0, 0.0]),
            &[1.0, 0.0, 0.0],
            1e-15,
        );
        let f = log::<f32>([1.0, 0.0, 0.0, 0.0]);
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn the_defects_are_the_correct_log_except_where_they_are_planted() {
        // Without the flip: off by exactly 2 pi n (`docs/maths/so3.md` SO.5(b)) for w < 0.
        let q = unit(2.0).map(|c| -c);
        let (good, bad) = (log(q), log_no_flip(q));
        let gap = (0..3)
            .map(|i| (bad[i] - good[i]).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!((gap - 2.0 * std::f64::consts::PI).abs() < 1e-14, "{gap}");
        let q = unit(2.0);
        assert_eq!(log(q).map(f64::to_bits), log_no_flip(q).map(f64::to_bits));
        // Through acos: fine at theta = 1, and at 1e-9 the trace rounds to 3 and theta to 0.
        close(&log_acos(unit(1.0)), &log(unit(1.0)), 1e-13);
        assert_eq!(log_acos(unit(1e-9)), [0.0; 3]);
        assert_eq!(log_acos([1.0, 0.0, 0.0, 0.0]), [0.0; 3]);
        assert_eq!(
            log_acos(unit(2.0)).map(f64::to_bits),
            log_acos(unit(2.0).map(|c| -c)).map(f64::to_bits)
        );
    }

    #[test]
    fn a_lane_that_evaluates_both_arms_never_sees_a_non_finite_operation() {
        // The safe arguments: `n²` in the exact arm, `w` in the series arm, at `w = 0`, `n = 0`,
        // a subnormal and either side of the mask, with the blend giving what the scalar path does.
        let (cand, series) = k();
        let lane_cand = Candidate {
            terms: cand.terms,
            switch_z: Lane(cand.switch_z),
        };
        let lane_series: Vec<Lane> = series.iter().map(|&s| Lane(s)).collect();
        for q in [
            [1.0, 0.0, 0.0, 0.0],
            [-1.0, -0.0, -0.0, -0.0],
            [0.0, 0.6, 0.0, 0.8],
            [1.0, 3e-310, 0.0, -3e-310],
            [0.5, 3e-310, 0.0, -3e-310],
            [1.0, 1e-160, 0.0, 0.0],
            [1.0, 1.6e-154, 0.0, 0.0],
            [0.8, 0.6, 0.0, 0.0],
            [-0.6, 0.0, 0.0, 0.8],
        ] {
            NONFINITE.set(0);
            let lanes = log(q.map(Lane));
            assert_eq!(NONFINITE.get(), 0, "log at {q:?}");
            assert_eq!(lanes.map(|c| c.0.to_bits()), log(q).map(f64::to_bits));
            let phi = [q[1], q[2], q[3]];
            let lanes = exp(phi.map(Lane), lane_cand, &lane_series);
            assert_eq!(NONFINITE.get(), 0, "exp at {phi:?}");
            assert_eq!(lanes.map(|c| c.0.to_bits()), exp64(phi).map(f64::to_bits));
        }
    }
}
