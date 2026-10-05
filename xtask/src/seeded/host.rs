//! The host-`std` twin: the correct seeded kernels at a scalar whose transcendentals are Rust
//! `std`'s — the host's libm, so glibc on Linux — where every other subject in this crate routes
//! them through the `libm` crate (D16).
//!
//! It exists to **attribute** a domination failure, never to score one. `0032` (draft) measured
//! `so3_log`'s 8 failures by swapping `atan2` in a script that was not kept; `0036` (draft) plan
//! step 2 owes the same swap for `sen3_jl_n1` and `sen3_jr_n1`, and `0032` (draft) open question 1 asks
//! whether D16's price is a per-function constant or a per-stratum measurement. One subject
//! answers all three, for every id the correct kernel runs, through the harness's own metric.
//!
//! **One variable.** [`Host`] differs from `f64` in `sin_cos`, `atan2` and `cbrt` only: `sqrt`,
//! `abs` and `copysign` are IEEE-exact operations, the same function in either library
//! (`0018`). The ids this twin answers reach `sin_cos` (every coefficient's exact arm,
//! `kernel::exact`) and `atan2` (`r`, so `so3_log`); none reaches `cbrt`. The coefficients, the
//! switches, the series, the associations and the application form are `seeded:correct`'s,
//! shared by construction — the two subjects run the same generic code over the same
//! `generated.rs`.
//!
//! **It is planted** (`Registered::planted`), so a plain `just conformance` skips it and no bar
//! ever reads its rows: a twin is not a candidate, and the envelope comparing `helicoid` with a
//! differently-rounded copy of itself would be a measurement of nothing.
//!
//! The twin has no `f32` kernel. `std`'s `f32` transcendentals are not the host's `f32` libm on
//! every target (LLVM may widen), so an `f32` twin would change two things, and `0016`'s `@f32`
//! strata are the coefficient ids' alone.

use helicoid_linalg::{Dual, Mask, Precision, Real};

use crate::conformance::corpus::Record;
use crate::conformance::subject::{Output, Registered, Subject};

use super::generated::{A_F64, B_F64, COS_HALF_F64, C_F64, D_F64, E_F64, K_F64, R_F64};
use super::kernel::{Candidate, Coeff, Swept};
use super::{answer, exp_at, kernels_of, log_at, se3, sen3_at, so3_jac_at, Arm};

/// The subject's name, and the value of `--subject`.
pub(super) const NAME: &str = "seeded:host-std";

/// The suffix a candidate's own host-`std` twin carries: `<candidate>:host-std`.
///
/// `seeded:correct`'s twin is [`NAME`] for historical reasons — it was the only candidate when
/// this module was written. The envelope looks for the suffixed name first, so a candidate that
/// has a twin of its own program is attributed by *that* and not by a stand-in's
/// (`judge_with`'s own words: "the candidate's own program with one variable changed").
pub(crate) const SUFFIX: &str = ":host-std";

/// An `f64` whose transcendentals come from Rust `std`, i.e. from the host's libm.
///
/// `std` has no `sincos`, so [`Real::sin_cos`] is two calls where the `libm` crate's is one.
/// That is the right shape for this comparison: sophus-rs, the oracle whose wins this twin is
/// here to explain, is Rust `std` and calls `sin()` and `cos()` separately too
/// (`conformance::Backend::HostStd`).
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub(crate) struct Host(f64);

impl From<f64> for Host {
    fn from(x: f64) -> Self {
        Host(x)
    }
}

impl From<Host> for f64 {
    fn from(x: Host) -> f64 {
        x.0
    }
}

macro_rules! binary {
    ($trait:ident, $method:ident, $op:tt) => {
        impl core::ops::$trait for Host {
            type Output = Self;
            fn $method(self, o: Self) -> Self {
                Host(self.0 $op o.0)
            }
        }
    };
}

binary!(Add, add, +);
binary!(Sub, sub, -);
binary!(Mul, mul, *);
binary!(Div, div, /);

impl core::ops::Neg for Host {
    type Output = Self;
    fn neg(self) -> Self {
        Host(-self.0)
    }
}

impl Real for Host {
    type Mask = bool;
    const PRECISION: Precision = Precision::F64;
    fn lit(x: f64) -> Self {
        Host(x)
    }
    fn zero() -> Self {
        Host(0.0)
    }
    fn one() -> Self {
        Host(1.0)
    }
    fn lt(self, rhs: Self) -> bool {
        self.0 < rhs.0
    }
    fn le(self, rhs: Self) -> bool {
        self.0 <= rhs.0
    }
    fn select(m: bool, t: Self, f: Self) -> Self {
        if m.all() {
            t
        } else {
            f
        }
    }
    // The six below are the whole of this subject. Three are IEEE-exact and agree with the `libm`
    // crate's bit for bit; `sin_cos`, `atan2` and `cbrt` are the measurement.
    fn sqrt(self) -> Self {
        Host(self.0.sqrt())
    }
    fn cbrt(self) -> Self {
        Host(self.0.cbrt())
    }
    fn sin_cos(self) -> (Self, Self) {
        (Host(self.0.sin()), Host(self.0.cos()))
    }
    fn atan2(self, x: Self) -> Self {
        Host(self.0.atan2(x.0))
    }
    fn abs(self) -> Self {
        Host(self.0.abs())
    }
    fn copysign(self, sign: Self) -> Self {
        Host(self.0.copysign(sign.0))
    }
    fn value_f64(self) -> f64 {
        self.0
    }
}

/// The correct kernels at [`Host`]: the generated switches, series and candidates of
/// `generated.rs`, converted exactly (a newtype of `f64`).
pub(super) struct Twin {
    /// By `Swept::index`, as `Seeded`'s are.
    arms: [Arm<Dual<Host, 1>>; 8],
}

/// One generated switch at [`Host`]. `Arm::generated` needs a `&[S]`, and the committed series is
/// `&[f64]`, so the terms are converted here instead.
fn arm((below, series): (f64, &'static [f64])) -> Arm<Dual<Host, 1>> {
    Arm {
        candidate: Candidate {
            terms: series.len(),
            switch_z: Dual::constant(Host(below)),
        },
        series: series.iter().map(|&x| Dual::constant(Host(x))).collect(),
    }
}

impl Twin {
    pub(super) fn registered() -> Registered {
        let arms = [
            K_F64.parts(),
            A_F64.parts(),
            B_F64.parts(),
            C_F64.parts(),
            D_F64.parts(),
            E_F64.parts(),
            COS_HALF_F64.parts(),
            R_F64.parts(),
        ]
        .map(arm);
        Registered {
            version: "generated@host-std".to_string(),
            version_f32: "generated@host-std".to_string(),
            no_f32: Some(
                "the twin measures the host's binary64 transcendentals; `0016`'s `@f32` strata \
                 are the coefficient ids' own"
                    .to_string(),
            ),
            planted: true,
            subject: Box::new(Twin { arms }),
        }
    }
}

impl Subject for Twin {
    fn name(&self) -> &str {
        NAME
    }

    /// Every id `seeded:correct` answers, so that the price of D16 is measured wherever the
    /// correct kernel is scored and not only where a bar happens to fail today.
    fn supports(&self, fn_id: &str) -> bool {
        Swept::of_fn(fn_id).is_some()
            || se3::parse(fn_id).is_some()
            || matches!(fn_id, "so3_exp" | "so3_log" | "so3_jr" | "so3_jl")
    }

    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        if !self.supports(fn_id) || precision != Precision::F64 {
            return Output::new();
        }
        match Swept::of_fn(fn_id) {
            Some(id) => answer(None, id, &self.arms[id.index()], record),
            None if fn_id == "so3_exp" => exp_at(record, &self.arms[Coeff::K.index()]),
            None if fn_id == "so3_log" => log_at::<Host>(record),
            None if matches!(fn_id, "so3_jr" | "so3_jl") => {
                so3_jac_at(fn_id, record, &kernels_of(&self.arms), se3::W2::Product)
            }
            None => sen3_at(
                fn_id,
                record,
                &kernels_of(&self.arms),
                se3::Order::RotationFirst,
                0.5,
                se3::W2::Product,
            ),
        }
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;
    use crate::conformance::corpus;
    use crate::conformance::{corpus_dir, evaluate};
    use crate::seeded::so3::norm_sq;

    /// The measurement table of `0032` (draft), its host-glibc column: the `so3_log` strata it
    /// reports, to the four decimals it prints. `theta:1e-7` is the value its prose gives
    /// (1.8492), not sophus-rs's 1.9396, which that stratum does not reproduce.
    const PUBLISHED: [(&str, f64); 8] = [
        ("theta:dense", 2.1605),
        ("theta:1e-1", 1.6643),
        ("theta:1e-2", 1.9876),
        ("theta:1e-3", 1.6312),
        ("theta:1e-4", 2.0877),
        ("theta:1e-5", 1.6585),
        ("theta:1e-6", 1.9923),
        ("theta:1e-7", 1.8492),
    ];

    /// The positive control. `0032` (draft) measured its glibc column with a script that was
    /// not kept; this subject must reproduce it, or it is measuring something else.
    #[test]
    fn the_twin_reproduces_the_glibc_column_of_0032() -> Result<(), String> {
        let dir = corpus_dir()?;
        let mut entries = corpus::manifest(&dir)?;
        entries.retain(|e| e.fn_id == "so3_log");
        let subjects = [Twin::registered()];
        let rows = evaluate(&dir, &entries, &subjects, Precision::F64)?;
        let got: BTreeMap<&str, f64> = rows
            .first()
            .ok_or("the twin answered no rows")?
            .iter()
            .map(|r| (r.stratum.as_str(), r.max_u))
            .collect();
        for (stratum, want) in PUBLISHED {
            let max = *got.get(stratum).ok_or(format!("no row for {stratum}"))?;
            assert!(
                (max - want).abs() < 5e-5,
                "{stratum}: {max} u, `0032` reports {want}"
            );
        }
        Ok(())
    }

    /// Whether the two libraries' `sin` and `cos` agree in every bit at `x`.
    fn agree(x: f64) -> bool {
        let (s, c) = libm::sincos(x);
        s == x.sin() && c == x.cos()
    }

    /// Per stratum of `fn_id`: how many records reach a `sin`/`cos` argument where the two
    /// libraries disagree, out of how many. The arguments are `kernel::exact`'s: `θ = sqrt(z)`
    /// and `θ/2`, `z = |φ|²` formed as `so3::norm_sq` forms it, read from the record's `key` —
    /// `tau` for the `sen3_*` ids, `phi` for `so3_j{r,l}`.
    fn power(fn_id: &str, key: &str) -> Result<BTreeMap<String, (usize, usize)>, String> {
        let path = corpus_dir()?.join(format!("{fn_id}.jsonl"));
        let text =
            std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        let mut out: BTreeMap<String, (usize, usize)> = BTreeMap::new();
        for line in text.lines() {
            let r = corpus::parse_line(line)?;
            let Some(&[x, y, z]) = r.input(key).and_then(|t| t.first_chunk::<3>()) else {
                continue;
            };
            let th = norm_sq([x, y, z]).sqrt();
            let cell = out.entry(r.stratum.clone()).or_default();
            cell.1 += 1;
            if !(agree(th) && agree(0.5 * th)) {
                cell.0 += 1;
            }
        }
        Ok(out)
    }

    /// The twin's power, which is what makes "no change" a reading and not an absence.
    ///
    /// `sen3_jl_n1` is 52 strata of 6 records (`docs/PHASE1.md` §4.4, a stated budget trade-off),
    /// and at one of its 312 records do the two libraries disagree on `sin` or `cos` at all. Near
    /// `π` they never do — so on the `theta:pi-*` and `theta=pi-1e-6` strata, which is where 24 of
    /// this id's 25 domination failures are, the swap **cannot** change an answer, and a failure
    /// there is not D16's whatever else it is (`0037`, draft).
    ///
    /// Pinned, not bounded: a `libm` or glibc version that moves these counts is a change this
    /// reading depends on, and `0034` (draft) asks to hear about it.
    #[test]
    fn the_swap_has_no_power_on_the_strata_near_pi() -> Result<(), String> {
        let per_stratum = power("sen3_jl_n1", "tau")?;
        assert_eq!(per_stratum.len(), 52);
        let near_pi = |name: &str| name.contains("pi-1e");
        for (stratum, (diverged, records)) in &per_stratum {
            assert_eq!(*records, 6, "{stratum}");
            if near_pi(stratum) {
                assert_eq!(*diverged, 0, "{stratum}");
            }
        }
        let total: usize = per_stratum.values().map(|&(d, _)| d).sum();
        assert_eq!(total, 1, "divergent records of sen3_jl_n1's 312");
        assert_eq!(per_stratum.get("theta:1e0").map(|&(d, _)| d), Some(1));
        Ok(())
    }

    /// The same reading for `so3_j{r,l}`, the two ids this twin gained with SO(3).
    ///
    /// Of the 2466 records of each id's 28 strata, 33 reach a `sin`/`cos` argument where the two
    /// libraries disagree — and **not one is on a `theta:pi-1e*` stratum**, which is where all 16
    /// of the pair's domination failures are (2026-10-04). So the swap cannot change an answer
    /// where they fail, and none of the 16 is D16's cost, whatever else it is (`0037`, draft).
    ///
    /// Pinned, not bounded, for the reason `the_swap_has_no_power_on_the_strata_near_pi` is.
    #[test]
    fn the_swap_has_no_power_near_pi_on_the_so3_jacobians() -> Result<(), String> {
        for fn_id in ["so3_jr", "so3_jl"] {
            let per_stratum = power(fn_id, "phi")?;
            assert_eq!(per_stratum.len(), 28, "{fn_id}");
            let mut near_pi = 0usize;
            for (stratum, (diverged, _)) in &per_stratum {
                if stratum.contains("pi-1e") {
                    near_pi += 1;
                    assert_eq!(*diverged, 0, "{fn_id} {stratum}");
                }
            }
            let total: usize = per_stratum.values().map(|&(d, _)| d).sum();
            let records: usize = per_stratum.values().map(|&(_, n)| n).sum();
            assert_eq!((near_pi, total, records), (12, 33, 2466), "{fn_id}");
        }
        Ok(())
    }

    /// The module's central claim: three of the six are the same function in either library, and
    /// three are the measurement. A scalar that differed in more than the transcendentals would
    /// make every reading above two variables instead of one.
    #[test]
    fn the_scalar_differs_from_the_libm_crate_in_the_transcendentals_only() {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let mut next = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            3.2 * ((state >> 11) as f64 / (1u64 << 53) as f64)
        };
        let (mut sin_cos, mut atan2, mut cbrt) = (0usize, 0usize, 0usize);
        for _ in 0..100_000 {
            let (x, y) = (next(), next());
            let h = Host(x);
            // Exact operations: the same bits, every time.
            assert_eq!(h.sqrt().0, libm::sqrt(x));
            assert_eq!(h.abs().0, libm::fabs(x));
            assert_eq!(h.copysign(Host(-y)).0, libm::copysign(x, -y));
            // The three that can differ.
            let (s, c) = h.sin_cos();
            let (ls, lc) = libm::sincos(x);
            sin_cos += usize::from(s.0 != ls || c.0 != lc);
            atan2 += usize::from(h.atan2(Host(y)).0 != libm::atan2(x, y));
            let wide = (x - 1.6) * libm::exp2(libm::round(40.0 * y - 64.0));
            cbrt += usize::from(Host(wide).cbrt().0 != libm::cbrt(wide));
        }
        // `sin_cos` and `atan2` disagree, so the swap is a change and not a relabelling.
        assert!(sin_cos > 0 && atan2 > 0, "sin_cos {sin_cos}, atan2 {atan2}");
        // `cbrt` does not, over signed arguments spanning `2^-64` to `2^64`: the two libraries
        // are the same function there, so D16 costs `0017`'s `cbrt` nothing on this host. No id
        // this twin answers reaches it; the sample is what makes that a reading.
        assert_eq!(cbrt, 0);
    }
}
