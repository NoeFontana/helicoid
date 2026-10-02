//! The coefficient kernels `k, a, b, c, d, e` of `docs/NUMERICS.md` §4, `cos θ/2` (§3.1) and `r`
//! (§3.2) as test subjects, generic over `S: Real` in the branch variable `z = θ²`, so
//! `Dual<f64, 1>` seeded at `z` gives `d/dz` (`docs/maths/coefficients.md` CO.13). They are not the
//! shipped kernel: nothing is typed that `helicoid::coeffs` will own. The series constants are read
//! from `coeff_series.jsonl` and rounded at `S` (`series`); the switch is the caller's
//! ([`Candidate`]), rounded at `S` too.
//!
//! The exact arm evaluates the "Exact arm computes" column at the safe argument, one `S::branch`
//! per call. Its operand order is the one CO.6 measured, so `θ²` inside a cancelling sum is
//! `fl(θ̂·θ̂)`, not `z` (CO.6: `z` there raises the constants); `NUMERICS.md` §4 fixes no order.
//!
//! **`r` is a reading, pending the maintainer, not a spec** (`NUMERICS.md` §4 names `n²` as its
//! branch variable and its series is in `n²/w²`; 0014 (draft) question 29, and 0015 (draft) NU.6,
//! which writes the mask division-free): the switch compares `s = n²/w²`, formed by a division,
//! the series arm is taken iff `w > 0` and `s < switch`, and `w = 0` (`s` infinite)
//! or `w < 0` (S² charts only) takes the exact arm. `Dual` is seeded on `n²` at fixed `w`, the
//! derivative the corpus stores. Each arm is at its safe argument: `n²` in the exact arm, `w` in
//! the series arm (CO.16(d)).

use helicoid_linalg::{Dual, Mask, Precision, Real};

use super::series::Series;

/// The dual number the harness evaluates at: value and `d/dz`.
pub(crate) type D1 = Dual<f64, 1>;

/// [`D1`] at binary32.
pub(crate) type D1F32 = Dual<f32, 1>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Coeff {
    K,
    A,
    B,
    C,
    D,
    E,
}

impl Coeff {
    pub(crate) const ALL: [Coeff; 6] = [Coeff::K, Coeff::A, Coeff::B, Coeff::C, Coeff::D, Coeff::E];

    /// The name in `coeff_series.jsonl` (`"k"`) and after `coeff_` in a corpus id.
    pub(crate) fn name(self) -> &'static str {
        match self {
            Coeff::K => "k",
            Coeff::A => "a",
            Coeff::B => "b",
            Coeff::C => "c",
            Coeff::D => "d",
            Coeff::E => "e",
        }
    }

    /// The position in [`Coeff::ALL`].
    pub(crate) fn index(self) -> usize {
        self as usize
    }

    pub(crate) fn of_name(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|c| c.name() == name)
    }

    pub(crate) fn of_fn(fn_id: &str) -> Option<Self> {
        Self::of_name(fn_id.strip_prefix("coeff_")?)
    }
}

/// What the sweep generates a switch for: [`Coeff`], `cos θ/2` (`Exp`'s quaternion, `NUMERICS.md`
/// §3.1) and `r` (`Log`'s ratio, §3.2), the eight rows of `coeff_series.jsonl`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Swept {
    Coeff(Coeff),
    CosHalf,
    R,
}

impl Swept {
    pub(crate) const ALL: [Swept; 8] = [
        Swept::Coeff(Coeff::K),
        Swept::Coeff(Coeff::A),
        Swept::Coeff(Coeff::B),
        Swept::Coeff(Coeff::C),
        Swept::Coeff(Coeff::D),
        Swept::Coeff(Coeff::E),
        Swept::CosHalf,
        Swept::R,
    ];

    /// The name in `coeff_series.jsonl` and after `coeff_` in a corpus id.
    pub(crate) fn name(self) -> &'static str {
        match self {
            Swept::Coeff(c) => c.name(),
            Swept::CosHalf => "cos_half",
            Swept::R => "r",
        }
    }

    /// The position in [`Swept::ALL`], [`Coeff::index`] for the six.
    pub(crate) fn index(self) -> usize {
        match self {
            Swept::Coeff(c) => c.index(),
            Swept::CosHalf => 6,
            Swept::R => 7,
        }
    }

    pub(crate) fn of_fn(fn_id: &str) -> Option<Self> {
        let name = fn_id.strip_prefix("coeff_")?;
        Self::ALL.into_iter().find(|s| s.name() == name)
    }
}

/// A record's arguments at `S`: `z = θ²` (`n²` for `r`) and, for `r` alone, `w`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Input<S> {
    pub(crate) z: S,
    pub(crate) w: S,
}

impl<S: Real> Input<S> {
    /// The same at `Dual<S, 1>` seeded on `z` (`w` fixed).
    pub(crate) fn seed(self) -> Input<Dual<S, 1>> {
        Input {
            z: Dual::variable(self.z, 0),
            w: Dual::constant(self.w),
        }
    }
}

/// The variable a switch is compared with: `z`, and `s = n²/w²` for `r`.
pub(crate) fn branch_variable<S: Real>(id: Swept, x: Input<S>) -> S {
    match id {
        Swept::R => x.z / (x.w * x.w),
        _ => x.z,
    }
}

/// The kernel of `id` at `x`.
pub(crate) fn evaluate<S: Real>(id: Swept, x: Input<S>, cand: Candidate<S>, series: &[S]) -> S {
    match id {
        Swept::Coeff(c) => coefficient(c, x.z, cand, series),
        Swept::CosHalf => cos_half(x.z, cand, series),
        Swept::R => log_ratio(x, cand, series),
    }
}

/// `terms` series terms below `switch_z` in the branch variable `θ²`, both at the scalar `S`. Build
/// it with [`Candidate::new`].
#[derive(Clone, Copy, Debug)]
pub(crate) struct Candidate<S> {
    pub(crate) terms: usize,
    pub(crate) switch_z: S,
}

impl<S: Real> Candidate<S> {
    /// `terms` of `series` below the decimal `switch`, rounded once at `S` by the standard
    /// library's parser: `Real::lit` takes only what `S` holds, and a decimal such as `0.01` is
    /// exact at neither precision (0003 item 6, no `f64 -> f32` double rounding).
    pub(crate) fn new(series: &Series<S>, terms: usize, switch: &str) -> Result<Self, String> {
        if !(1..=series.terms()).contains(&terms) {
            return Err(format!("{terms} terms of a {}-term series", series.terms()));
        }
        let bad = |e: std::num::ParseFloatError| format!("switch `{switch}`: {e}");
        let z = match S::PRECISION {
            Precision::F64 => switch.parse::<f64>().map_err(bad)?,
            Precision::F32 => f64::from(switch.parse::<f32>().map_err(bad)?),
        };
        if !(z.is_finite() && z >= 0.0) {
            return Err(format!("switch `{switch}` is not a finite z >= 0"));
        }
        Ok(Self {
            terms,
            switch_z: S::lit(z),
        })
    }
}

/// `tf_tree` D12 (`θ < 0.1`, four terms), the named prior of `docs/PHASE1.md` §6, which
/// `NUMERICS.md` §4 lists for `a`, `b`, `c` only. Here the evaluation candidate of all six: for
/// `k`, `d`, `e` it is D12 applied to a coefficient it was not defined for. The sweep also takes
/// it for `cos θ/2` and, in `s`, for `r`, where `(4, 0.01)` is `θ = 0.2` and not D12's
/// (0014 (draft) question 29). Not a generated switch (0004) and not a claim of optimality.
pub(crate) fn d12<S: Real>(series: &Series<S>) -> Result<Candidate<S>, String> {
    Candidate::new(series, 4, "0.01")
}

/// `Σ series[j] z^j` over the first `terms`, `p_j = s_j + z·p_{j+1}` (CO.9): no `mul_add`, and
/// its `Dual` is exactly the derivative of the polynomial.
pub(crate) fn horner<S: Real>(series: &[S], terms: usize, z: S) -> S {
    let step = |p: S, &s: &S| s + z * p;
    series.iter().take(terms).rev().fold(S::zero(), step)
}

/// The exact arm at `z > 0`.
fn exact<S: Real>(c: Coeff, z: S) -> S {
    let (half, two, three, four) = (S::lit(0.5), S::lit(2.0), S::lit(3.0), S::lit(4.0));
    let th = z.sqrt();
    let t2 = th * th;
    let sin_half = || (half * th).sin_cos().0;
    match c {
        Coeff::K => sin_half() / th,
        Coeff::A => {
            let k = sin_half() / th;
            two * k * k
        }
        Coeff::B => (th - th.sin_cos().0) / (t2 * th),
        Coeff::C => {
            let (s, co) = (half * th).sin_cos();
            S::one() / t2 - co / (two * th * s)
        }
        Coeff::D => {
            let s = sin_half();
            (t2 - four * s * s) / (two * t2 * t2)
        }
        Coeff::E => {
            let (s, co) = th.sin_cos();
            ((two * th - three * s) + th * co) / (two * th * (t2 * t2))
        }
    }
}

/// The coefficient: the series below `cand.switch_z`, else the exact arm at the safe argument
/// (0003 item 3), in one `S::branch`.
pub(crate) fn coefficient<S: Real>(c: Coeff, z: S, cand: Candidate<S>, series: &[S]) -> S {
    let small = z.lt(cand.switch_z);
    S::branch(
        small,
        || horner(series, cand.terms, z),
        || exact(c, S::select(small, S::one(), z)),
    )
}

/// `cos(θ/2)`: the series below `cand.switch_z`, else the definition at the safe argument.
fn cos_half<S: Real>(z: S, cand: Candidate<S>, series: &[S]) -> S {
    let small = z.lt(cand.switch_z);
    S::branch(
        small,
        || horner(series, cand.terms, z),
        || {
            (S::lit(0.5) * S::select(small, S::one(), z).sqrt())
                .sin_cos()
                .1
        },
    )
}

/// `r = 2 atan2(n, w)/n`, `n = sqrt(n²)`: the series in `s = n²/w²` iff `w > 0` and
/// `s < cand.switch_z`, else the definition.
fn log_ratio<S: Real>(x: Input<S>, cand: Candidate<S>, series: &[S]) -> S {
    let Input { z: n2, w } = x;
    let small = S::zero()
        .lt(w)
        .and(branch_variable(Swept::R, x).lt(cand.switch_z));
    let two = S::lit(2.0);
    S::branch(
        small,
        || {
            let w = S::select(small, w, S::one());
            two / w * horner(series, cand.terms, n2 / (w * w))
        },
        || {
            let n = S::select(small, S::one(), n2).sqrt();
            two * n.atan2(w) / n
        },
    )
}

/// Defect: `b` by its definition everywhere, no series and no branch (`docs/PHASE1.md` §10).
pub(crate) fn b_no_series<S: Real>(z: S) -> S {
    exact(Coeff::B, z)
}

/// Defect: `k` with one `sqrt(z)` shared by both arms and no safe argument, its series arm
/// evaluated at `θ·θ` rebuilt from that root. Under `Dual` the root's derivative is infinite at
/// `z = 0`, so the selected series arm returns a finite value and a NaN derivative
/// (`docs/maths/error-analysis.md` EA.18(b)). It does not model an exact arm at `z` without the
/// safe argument: a `bool` mask never evaluates the unselected arm, so that omission is invisible
/// to any `bool` subject. `a_lane_that_evaluates_both_arms_sees_the_missing_safe_argument`
/// covers it.
pub(crate) fn k_sqrt_unsafe<S: Real>(z: S, cand: Candidate<S>, series: &[S]) -> S {
    let th = z.sqrt();
    S::branch(
        z.lt(cand.switch_z),
        || horner(series, cand.terms, th * th),
        || (S::lit(0.5) * th).sin_cos().0 / th,
    )
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
pub(super) mod tests {
    use std::cell::Cell;

    use helicoid_linalg::Mask;

    use super::*;
    use crate::conformance::corpus_dir;
    use crate::conformance::testkit::splitmix;

    thread_local!(pub(in crate::seeded) static NONFINITE: Cell<usize> = const { Cell::new(0) });

    /// An `f64` that counts every operation with a non-finite result.
    #[derive(Clone, Copy, Debug)]
    pub(in crate::seeded) struct Lane(pub(in crate::seeded) f64);

    fn note(x: f64) -> Lane {
        if !x.is_finite() {
            NONFINITE.set(NONFINITE.get() + 1);
        }
        Lane(x)
    }

    /// A mask that evaluates both arms and blends, as a SIMD lane does (`docs/PHASE2.md` §2).
    #[derive(Clone, Copy)]
    pub(in crate::seeded) struct Both(bool);

    impl Mask for Both {
        fn and(self, o: Self) -> Self {
            Both(self.0 & o.0)
        }
        fn or(self, o: Self) -> Self {
            Both(self.0 | o.0)
        }
        fn not(self) -> Self {
            Both(!self.0)
        }
        fn all(self) -> bool {
            self.0
        }
        fn any(self) -> bool {
            self.0
        }
        fn decide<T>(
            self,
            t: impl FnOnce() -> T,
            f: impl FnOnce() -> T,
            blend: impl FnOnce(Self, T, T) -> T,
        ) -> T {
            let (on, off) = (t(), f());
            blend(self, on, off)
        }
    }

    impl std::ops::Add for Lane {
        type Output = Self;
        fn add(self, o: Self) -> Self {
            note(self.0 + o.0)
        }
    }
    impl std::ops::Sub for Lane {
        type Output = Self;
        fn sub(self, o: Self) -> Self {
            note(self.0 - o.0)
        }
    }
    impl std::ops::Mul for Lane {
        type Output = Self;
        fn mul(self, o: Self) -> Self {
            note(self.0 * o.0)
        }
    }
    impl std::ops::Div for Lane {
        type Output = Self;
        fn div(self, o: Self) -> Self {
            note(self.0 / o.0)
        }
    }
    impl std::ops::Neg for Lane {
        type Output = Self;
        fn neg(self) -> Self {
            Lane(-self.0)
        }
    }

    impl Real for Lane {
        type Mask = Both;
        const PRECISION: Precision = Precision::F64;
        fn lit(x: f64) -> Self {
            Lane(x)
        }
        fn zero() -> Self {
            Lane(0.0)
        }
        fn one() -> Self {
            Lane(1.0)
        }
        fn lt(self, rhs: Self) -> Both {
            Both(self.0 < rhs.0)
        }
        fn le(self, rhs: Self) -> Both {
            Both(self.0 <= rhs.0)
        }
        fn select(m: Both, t: Self, f: Self) -> Self {
            if m.0 {
                t
            } else {
                f
            }
        }
        fn sqrt(self) -> Self {
            note(libm::sqrt(self.0))
        }
        fn cbrt(self) -> Self {
            note(libm::cbrt(self.0))
        }
        fn sin_cos(self) -> (Self, Self) {
            let (s, c) = libm::sincos(self.0);
            (note(s), note(c))
        }
        fn atan2(self, x: Self) -> Self {
            note(libm::atan2(self.0, x.0))
        }
        fn abs(self) -> Self {
            Lane(libm::fabs(self.0))
        }
        fn copysign(self, s: Self) -> Self {
            Lane(libm::copysign(self.0, s.0))
        }
        fn value_f64(self) -> f64 {
            self.0
        }
    }

    type DLane = Dual<Lane, 1>;

    /// The series and the D12 candidate at `S`.
    fn kernel<S: Real>() -> Result<(Series<S>, Candidate<S>), String> {
        let series = Series::<S>::load(&corpus_dir()?)?;
        let cand = d12(&series)?;
        Ok((series, cand))
    }

    fn at(z: f64) -> D1 {
        D1::variable(z, 0)
    }

    fn near(got: f64, want: f64, tol: f64) -> bool {
        (got - want).abs() <= tol * want.abs()
    }

    #[test]
    fn the_exact_arms_agree_with_the_committed_series_in_value_and_derivative() -> Result<(), String>
    {
        let (series, _) = kernel::<D1>()?;
        for c in Coeff::ALL {
            for z in [0.3, 1.0, 3.0, 6.0] {
                let (e, s) = (exact(c, at(z)), horner(series.of(c), 16, at(z)));
                assert!(
                    near(e.v, s.v, 1e-9),
                    "{c:?} value at {z}: {} vs {}",
                    e.v,
                    s.v
                );
                assert!(
                    near(e.d[0], s.d[0], 1e-9),
                    "{c:?} d/dz at {z}: {} vs {}",
                    e.d[0],
                    s.d[0]
                );
            }
        }
        Ok(())
    }

    #[test]
    fn the_exact_arm_is_non_finite_at_zero_and_the_series_arm_answers_there() -> Result<(), String>
    {
        let (series, cand) = kernel::<D1>()?;
        for c in Coeff::ALL {
            let unsafe_arm = exact(c, at(0.0));
            assert!(
                !unsafe_arm.v.is_finite() && !unsafe_arm.d[0].is_finite(),
                "{c:?}"
            );
            // The series arm answers at 0 with its first two terms, the value and the derivative.
            let k = coefficient(c, at(0.0), cand, series.of(c));
            let s = series.of(c);
            assert_eq!(
                (k.v.to_bits(), k.d[0].to_bits()),
                (s[0].v.to_bits(), s[1].v.to_bits())
            );
        }
        Ok(())
    }

    #[test]
    fn the_exact_arms_keep_the_digits_the_naive_forms_lose() -> Result<(), String> {
        // The exact arm errs by `c u theta^-p` for theta in [1e-3, 1e-1] (CO.6: `c` = 1, 6, 47, 44,
        // 362 for `a, b, c, d, e`, p = 0, 2, 2, 2, 4), here against the 16-term series in plain
        // `f64`, good to a few `u`; measured 3.0, 5.5, 38, 31, 324, and 1.0 for `k`. `a` as
        // `(1 - cos)/theta^2` errs by `2/theta^2`, `d` without the `sin^2` rewrite likewise.
        let (series, _) = kernel::<f64>()?;
        let u = 2f64.powi(-53);
        for (c, p, bound) in [
            (Coeff::K, 0, 2.0),
            (Coeff::A, 0, 5.0),
            (Coeff::B, 2, 8.0),
            (Coeff::C, 2, 60.0),
            (Coeff::D, 2, 50.0),
            (Coeff::E, 4, 450.0),
        ] {
            let mut worst = 0.0f64;
            for i in 0..=60 {
                let theta = 10f64.powf(f64::from(i) / 30.0 - 3.0);
                let z = theta * theta;
                let want = horner(series.of(c), 16, z);
                let scaled = ((exact(c, z) - want) / want).abs() / u * theta.powi(p);
                worst = worst.max(scaled);
            }
            assert!(worst <= bound, "{c:?}: {worst} u theta^-{p}");
        }
        Ok(())
    }

    /// `coefficient` with the exact arm at `z` itself: the safe argument left out.
    fn without_the_safe_argument<S: Real>(c: Coeff, z: S, cand: Candidate<S>, series: &[S]) -> S {
        S::branch(
            z.lt(cand.switch_z),
            || horner(series, cand.terms, z),
            || exact(c, z),
        )
    }

    /// The `z` that matter: zero, subnormal, tiny, either side of the switch, and above.
    const ZS: [f64; 10] = [
        0.0, 5e-324, 1e-310, 1e-16, 0.005, 0.0099, 0.01, 0.5, 6.0, 39.0,
    ];

    #[test]
    fn a_lane_that_evaluates_both_arms_never_sees_a_non_finite_operation() -> Result<(), String> {
        let (lanes, lane_cand) = kernel::<DLane>()?;
        let (plain, plain_cand) = kernel::<D1>()?;
        for c in Coeff::ALL {
            for z in ZS {
                NONFINITE.set(0);
                let lane = coefficient(c, DLane::variable(Lane(z), 0), lane_cand, lanes.of(c));
                let bits = (lane.v.0.to_bits(), lane.d[0].0.to_bits());
                assert_eq!(NONFINITE.get(), 0, "{c:?} at {z}: a non-finite operation");
                // The blend selects what the scalar `bool` path computes.
                let want = coefficient(c, at(z), plain_cand, plain.of(c));
                assert_eq!(
                    bits,
                    (want.v.to_bits(), want.d[0].to_bits()),
                    "{c:?} at {z}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn a_lane_that_evaluates_both_arms_sees_the_missing_safe_argument() -> Result<(), String> {
        let (lanes, lane_cand) = kernel::<DLane>()?;
        let (plain, plain_cand) = kernel::<D1>()?;
        for c in Coeff::ALL {
            NONFINITE.set(0);
            let unsafe_lane =
                without_the_safe_argument(c, DLane::variable(Lane(0.0), 0), lane_cand, lanes.of(c));
            assert!(NONFINITE.get() > 0, "{c:?}: the exact arm at 0 is not seen");
            // The blend still selects the series arm, and a `bool` never runs the exact arm: the
            // answer is right and the omission invisible to every `bool` subject.
            let bool_mask = without_the_safe_argument(c, at(0.0), plain_cand, plain.of(c));
            assert_eq!(unsafe_lane.v.0.to_bits(), bool_mask.v.to_bits(), "{c:?}");
            assert!(
                bool_mask.v.is_finite() && bool_mask.d[0].is_finite(),
                "{c:?}"
            );
        }
        Ok(())
    }

    #[test]
    fn the_series_arm_is_the_horner_of_exactly_terms_terms_below_the_switch() -> Result<(), String>
    {
        let (series, cand) = kernel::<f64>()?;
        for c in Coeff::ALL {
            let s = series.of(c);
            for z in [0.0, 1e-9, 0.0099, 0.01f64.next_down()] {
                let four = s[0] + z * (s[1] + z * (s[2] + z * s[3]));
                assert_eq!(
                    coefficient(c, z, cand, s).to_bits(),
                    four.to_bits(),
                    "{c:?}"
                );
            }
            // At the switch the exact arm answers, at the safe argument: `z < switch`, not `<=`.
            let z = 0.01;
            assert_eq!(
                coefficient(c, z, cand, s).to_bits(),
                exact(c, z).to_bits(),
                "{c:?}"
            );
        }
        Ok(())
    }

    #[test]
    fn a_candidate_needs_terms_the_series_has_and_a_switch_it_can_round() -> Result<(), String> {
        let (series, cand) = kernel::<f64>()?;
        assert_eq!((cand.terms, cand.switch_z), (4, 0.01));
        for (terms, switch) in [
            (0, "0.01"),
            (series.terms() + 1, "0.01"),
            (4, "-1e-3"),
            (4, "inf"),
            (4, "nan"),
            (4, ""),
            (4, "1/100"),
        ] {
            assert!(
                Candidate::new(&series, terms, switch).is_err(),
                "{terms} {switch}"
            );
        }
        assert!(Candidate::new(&series, series.terms(), "0").is_ok());
        Ok(())
    }

    #[test]
    fn the_kernels_run_at_binary32_on_constants_rounded_there() -> Result<(), String> {
        let (series, cand) = kernel::<f32>()?;
        let (wide, wide_cand) = kernel::<f64>()?;
        // 0.01 is exact at neither precision: the candidate holds binary32's nearest.
        assert_eq!(cand.switch_z.to_bits(), 0.01f32.to_bits());
        assert_ne!(f64::from(cand.switch_z), 0.01);
        for c in Coeff::ALL {
            for z in [0.0, 1e-30, 1e-3, 0.0099, 0.5, 2.0, 6.0] {
                let (narrow, exact) = (
                    coefficient(c, z as f32, cand, series.of(c)),
                    coefficient(c, z, wide_cand, wide.of(c)),
                );
                // The exact arms at binary32 err by `c u theta^-p`: 6e-8 * 362 * 4 at worst here.
                assert!(near(f64::from(narrow), exact, 1e-3), "{c:?} at {z}");
            }
        }
        Ok(())
    }

    #[test]
    fn the_dual_value_path_is_the_plain_value() -> Result<(), String> {
        let (plain, plain_cand) = kernel::<f64>()?;
        let (dual, dual_cand) = kernel::<D1>()?;
        let mut state = 9;
        for i in 0..4000 {
            // Log-uniform z over 1e-12 .. 10, the series and the exact arm alike.
            let u = (splitmix(&mut state) >> 11) as f64 / (1u64 << 53) as f64;
            let z = if i == 0 {
                0.0
            } else {
                10f64.powf(-12.0 + 13.0 * u)
            };
            for c in Coeff::ALL {
                let p = coefficient(c, z, plain_cand, plain.of(c));
                let d = coefficient(c, at(z), dual_cand, dual.of(c));
                assert_eq!(p.to_bits(), d.v.to_bits(), "{c:?} at {z}");
            }
        }
        Ok(())
    }

    #[test]
    fn cos_half_is_the_series_below_its_switch_and_the_definition_above() -> Result<(), String> {
        let series = Series::<D1>::load(&corpus_dir()?)?;
        let (s, cand) = (
            series.swept(Swept::CosHalf),
            Candidate::new(&series, 8, "0.5")?,
        );
        let x = |z| {
            evaluate(
                Swept::CosHalf,
                Input {
                    z: at(z),
                    w: D1::one(),
                },
                cand,
                s,
            )
        };
        // At `z = 0` the definition is `0/0` in the derivative and the series answers: `1`, `-1/8`.
        let zero = x(0.0);
        assert_eq!(
            (zero.v.to_bits(), zero.d[0].to_bits()),
            (1f64.to_bits(), (-0.125f64).to_bits())
        );
        // Above the switch: `cos(θ/2)` and `d/dz = -sin(θ/2)/(4θ)`, `θ = sqrt(z)`.
        for z in [0.5, 2.0, 6.0] {
            let (theta, got) = (f64::sqrt(z), x(z));
            let (sin, cos) = (theta / 2.0).sin_cos();
            assert!(
                near(got.v, cos, 1e-15) && near(got.d[0], -sin / (4.0 * theta), 1e-14),
                "{z}"
            );
        }
        Ok(())
    }

    #[test]
    fn r_takes_the_series_iff_w_is_positive_and_s_is_below_the_switch() -> Result<(), String> {
        let series = Series::<D1>::load(&corpus_dir()?)?;
        let (s, cand) = (series.swept(Swept::R), Candidate::new(&series, 8, "0.01")?);
        let r = |n2: f64, w: f64| {
            evaluate(
                Swept::R,
                Input {
                    z: at(n2),
                    w: D1::constant(w),
                },
                cand,
                s,
            )
        };
        // `2 atan2(n, w)/n` and its `d/dn²` at fixed `w`, `(w/(n² + w²) - r/2)/n²` (CO.13(a)).
        let definition = |n2: f64, w: f64| {
            let (n, v) = (n2.sqrt(), 2.0 * n2.sqrt().atan2(w) / n2.sqrt());
            (v, (w / (n2 + w * w) - v / 2.0) / (n * n))
        };
        // The series arm (`s = 2.5e-5`): the series in `s` times `2/w`, no `atan2`.
        let (got, (v, d)) = (r(1e-4, 2.0), definition(1e-4, 2.0));
        assert!(near(got.v, v, 1e-15) && near(got.d[0], d, 1e-8));
        let (two, w) = (D1::constant(2.0), D1::constant(2.0));
        let plain = two / w * horner(s, 8, at(1e-4) / (w * w));
        assert_eq!(got.v.to_bits(), plain.v.to_bits());
        assert_eq!(got.d[0].to_bits(), plain.d[0].to_bits());
        // `n² = 0` is `s = 0`: the series arm, `2/w` and `-2/(3w³)`, finite where `0/0` is not.
        let zero = r(0.0, 2.0);
        assert_eq!((zero.v, zero.d[0]), (1.0, -1.0 / 12.0));
        // Above the switch, at `w = 0` (`s` infinite) and at `w < 0` with `s` below it (CO.16(c):
        // the series would give -2.0000 for 6281.19) the exact arm answers: `π/n`, `-π/(2n³)`.
        let (pi, exact) = (std::f64::consts::PI, r(1e-2, 0.0));
        assert!(near(exact.v, pi / 0.1, 1e-15) && near(exact.d[0], -pi / 2e-3, 1e-14));
        let negative = r(1e-6, -1.0);
        assert!(
            near(negative.v, 6281.19, 1e-6) && near(negative.v, definition(1e-6, -1.0).0, 1e-15)
        );
        // The mask is on `s = n²/w²` and not on `n²`: at `w = 10`, `n² = 0.5` is far above the
        // switch and `s = 5e-3` below it, so the series arm answers; at `w = 0.05`, `n² = 1e-4` is
        // far below it and `s = 0.04` above, so the exact arm does. Each is the arm's own
        // operations bit for bit, and the two arms differ there.
        let series_arm = |n2: f64, w: f64| {
            let w = D1::constant(w);
            D1::constant(2.0) / w * horner(s, 8, at(n2) / (w * w))
        };
        let exact_arm = |n2: f64, w: f64| {
            let n = at(n2).sqrt();
            D1::constant(2.0) * n.atan2(D1::constant(w)) / n
        };
        for (n2, w, series) in [(0.5, 10.0, true), (1e-4, 0.05, false)] {
            let (got, want, other) = match series {
                true => (r(n2, w), series_arm(n2, w), exact_arm(n2, w)),
                false => (r(n2, w), exact_arm(n2, w), series_arm(n2, w)),
            };
            let bits = |x: D1| [x.v.to_bits(), x.d[0].to_bits()];
            assert_eq!(bits(got), bits(want), "n² = {n2}, w = {w}");
            assert_ne!(bits(got), bits(other), "n² = {n2}, w = {w}");
        }
        Ok(())
    }

    #[test]
    fn the_planted_k_defect_is_a_dual_defect_at_zero_only() -> Result<(), String> {
        let (plain, plain_cand) = kernel::<f64>()?;
        let (dual, dual_cand) = kernel::<D1>()?;
        let (s, sd) = (plain.of(Coeff::K), dual.of(Coeff::K));
        let bare = k_sqrt_unsafe(0.0f64, plain_cand, s);
        assert_eq!(bare.to_bits(), s[0].to_bits());
        let d = k_sqrt_unsafe(at(0.0), dual_cand, sd);
        assert_eq!(d.v.to_bits(), s[0].to_bits());
        assert!(d.d[0].is_nan());
        for z in [1e-300, 1e-8, 1e-3, 0.5] {
            let (bad, good) = (
                k_sqrt_unsafe(at(z), dual_cand, sd),
                coefficient(Coeff::K, at(z), dual_cand, sd),
            );
            assert!(
                near(bad.v, good.v, 1e-12) && near(bad.d[0], good.d[0], 1e-9),
                "{z}"
            );
        }
        Ok(())
    }

    #[test]
    fn b_by_its_definition_is_the_exact_arm_everywhere() {
        // Below sqrt(6u) = 2.6e-8 in theta the difference rounds to 0: all digits lost.
        assert_eq!(b_no_series(1e-16_f64).to_bits(), 0);
        assert!(near(b_no_series(1.0_f64), 1.0 - 1.0_f64.sin(), 1e-15));
    }
}
