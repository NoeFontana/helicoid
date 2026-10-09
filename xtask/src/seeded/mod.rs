//! Seeded subjects (`docs/PHASE1.md` §10): the correct coefficient kernels and the planted
//! defects, run over the `coeff_*` corpus ids as in-process subjects. The correct subject runs all
//! eight ids of [`Swept`] (`k, a, b, c, d, e`, `cos θ/2`, `r`); the planted defects run the six of
//! [`Coeff`] only. `Log` (`so3`) still evaluates `r` inline and does not read `R_F64`.
//!
//! The correct subject also runs `Exp` and `Log` of SO(3) over `so3_exp` and `so3_log` (`so3`), on
//! the generated `k`; the planted `Log` defects run over `so3_log` only. SO(3)'s `J_r` and `J_l`
//! over `so3_jr` and `so3_jl` are the rotation block of the SE_N(3) ones, on the generated `a, b`
//! (`se3::so3_jacobian`), and no planted defect reaches them. It runs SE_N(3)'s `Exp`,
//! `J_r` and `J_l` over `sen3_{exp,jr,jl}_n{1,2,3}` (`se3`) on the generated `k, a, b, d, e`; the
//! planted SE(3) defects run over `sen3_exp_n{1,2,3}` and over the two Jacobians.
//!
//! The correct kernel runs the **generated** switches of each coefficient: its series length,
//! switch and series terms are `generated.rs`'s, which `cargo xtask thresholds` writes
//! (`crate::thresholds`), and nothing here is typed. The sweep measures the same kernel through
//! [`evaluate`] with the candidate it scores, so the constants the subject runs are the
//! candidate the sweep chose (`the_generated_kernel_scores_the_objective_the_sweep_chose`).
//!
//! At `f32` the same kernels run through the same adapter on `Dual<f32, 1>` over the `@f32` strata
//! (`docs/decisions/0016` item 2): the inputs cast exactly, `z = fl(θ·θ)` (`n²` for `r`) formed at
//! `f32`, the series and switches the `_F32` constants of `generated.rs`, which the `f32` sweep
//! writes (item 3) from the exact rationals rounded once at `f32`.
//!
//! The `f32` kernels model the correct kernels and the defects `b` and `k`. A subject whose
//! defect they do not model (the planted `c`, a candidate of the binary64 sweep; `uniform`) has no
//! `f32` kernel (`Registered::no_f32`): a run at `f32` is an error, not the correct answer under its
//! name.
//!
//! The planted defects are the correct kernels with one coefficient changed. `c` with a switch of
//! `1e-8` and two terms ([`C_PLANTED`]) is named by the sweep's ranking (`conformance::selftest`);
//! its other mechanism, the envelope, is not started. The subject receives `θ`, forms the branch
//! variable `z = fl(θ·θ)` itself (`conformance/generate/README.md`) and reports the value and
//! `d/dz` of one `Dual<f64, 1>` evaluation, so the value path is the plain one
//! (`the_dual_value_path_is_the_plain_value`).

mod generated;
mod host;
mod kernel;
pub(crate) mod se3;
mod series;
mod so3;
mod switch;

use helicoid_linalg::{Dual, Precision, Real};

use crate::conformance::corpus::{exact_f32, Record};
use crate::conformance::subject::{Output, Registered, Subject};
use generated::{
    A_F32, A_F64, B_F32, B_F64, COS_HALF_F32, COS_HALF_F64, C_F32, C_F64, D_F32, D_F64, E_F32,
    E_F64, K_F32, K_F64, R_F32, R_F64,
};
use kernel::D1F32;
use kernel::{b_no_series, k_sqrt_unsafe};
pub(crate) use kernel::{branch_variable, d12, evaluate, Candidate, Coeff, Input, Swept, D1};
pub(crate) use series::{Series, FILE as SERIES_FILE};
pub(crate) use so3::takes_series_arm;

/// The one-variable twin's subject name (`0037`, draft): `seeded:correct`'s program on the host's
/// transcendentals, which the envelope reads to attribute a domination failure.
pub(crate) const TWIN: &str = host::NAME;

/// The scalar those twins run on, and the suffix that names one: an `f64` whose transcendentals
/// are the host's. Public to the crate so a twin of *any* candidate generic over `Real` — the
/// shipped library included — is that candidate's own program at this scalar.
pub(crate) use host::{Host, SUFFIX as TWIN_SUFFIX};

/// The assembly-form twin's subject name (`0037`, draft): `seeded:correct` with `J_l`'s `W²` by
/// `φφᵀ − θ²I` instead of `W·W`, and nothing else changed.
pub(crate) const TWIN_W2: &str = "seeded:w2-identity";

/// The planted `c` (`docs/PHASE1.md` §10): two terms below `z = 1e-16`, §10's `1e-8` read as `θ`
/// (0014 (draft) question 14), which is the first point of the sweep's grid.
pub(crate) const C_PLANTED: (usize, f64) = (2, 1e-16);

/// A planted defect (`docs/PHASE1.md` §10), applied to one coefficient; the others stay correct.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Defect {
    /// `b` by its definition, no series.
    BNoSeries,
    /// `k` with an unsafe `sqrt` of `θ²`, under `Dual`.
    KSqrtUnsafe,
    /// `c` with switch `1e-8` and two terms ([`C_PLANTED`]).
    CTwoTermsEarly,
    /// `Log` through `acos` of the trace (`so3::log_acos`).
    LogAcos,
    /// `Log` without the `w < 0` flip (`so3::log_no_flip`).
    LogNoFlip,
    /// SE(3) `Exp` reading the tangent translation-first (`se3::Order::TranslationFirst`).
    Se3ExpTranslationFirst,
    /// `Q` with `-½ρ^` where `NUMERICS.md` §5.3 has `+½ρ^` (`se3::jacobian`'s `half`).
    QMinusHalf,
}

impl Defect {
    /// The defects of one coefficient, run over the `coeff_*` ids.
    pub(crate) const COEFFICIENT: [Defect; 3] = [
        Defect::BNoSeries,
        Defect::KSqrtUnsafe,
        Defect::CTwoTermsEarly,
    ];
    /// The defects of `Log`, run over `so3_log`.
    pub(crate) const LOG: [Defect; 2] = [Defect::LogAcos, Defect::LogNoFlip];
    /// The defects of SE(3), run over `sen3_exp_n*` and `sen3_{jr,jl}_n*`.
    pub(crate) const SE3: [Defect; 2] = [Defect::Se3ExpTranslationFirst, Defect::QMinusHalf];

    pub(crate) fn name(self) -> &'static str {
        match self {
            Defect::BNoSeries => "b-no-series",
            Defect::KSqrtUnsafe => "k-sqrt-unsafe",
            Defect::CTwoTermsEarly => "c-two-terms-1e-8",
            Defect::LogAcos => "log-acos",
            Defect::LogNoFlip => "log-no-flip",
            Defect::Se3ExpTranslationFirst => "se3-exp-translation-first",
            Defect::QMinusHalf => "q-minus-half",
        }
    }
}

/// One coefficient's kernel: how many terms below which switch, and the series they are the first
/// terms of.
struct Arm<S> {
    candidate: Candidate<S>,
    series: Vec<S>,
}

impl<S: Real> Arm<Dual<S, 1>> {
    /// The generated switch `(below, series)`: every term of the series, below `below`.
    fn generated((below, series): (S, &[S])) -> Self {
        Self {
            candidate: Candidate {
                terms: series.len(),
                switch_z: Dual::constant(below),
            },
            series: series.iter().map(|&x| Dual::constant(x)).collect(),
        }
    }
}

/// The scalar input `key` of `record` at `S`: at `f32` only if it is exactly binary32 (a lossless
/// cast), so a corpus `@f32` stratum reaches the subject unrounded. Every coefficient id's input
/// goes through this, swept or not.
pub(crate) fn exact_input<S: Real>(record: &Record, key: &str) -> Option<S> {
    let x = *record.input(key)?.first()?;
    match S::PRECISION {
        Precision::F64 => Some(S::lit(x)),
        Precision::F32 => exact_f32(x).map(|x| S::lit(f64::from(x))),
    }
}

/// The arguments of `id` in `record` at `S`: `z = fl(θ·θ)` (`n²` for `r`) formed at `S`, as the
/// sweep forms it, from [`exact_input`].
pub(crate) fn input<S: Real>(id: Swept, record: &Record) -> Option<Input<S>> {
    let get = |key: &str| exact_input::<S>(record, key);
    Some(match id {
        Swept::R => {
            let n = get("n")?;
            Input {
                z: n * n,
                w: get("w")?,
            }
        }
        _ => {
            let theta = get("theta")?;
            Input {
                z: theta * theta,
                w: S::one(),
            }
        }
    })
}

pub(crate) struct Seeded {
    name: String,
    version: String,
    defect: Option<Defect>,
    /// How the `b W²` term of `J_l` is formed. `seeded:correct` and every defect run
    /// `W2::Product`; the twin of 0037 (draft) runs `W2::Identity` and changes nothing else.
    form: se3::W2,
    /// By [`Swept::index`].
    arms: [Arm<D1>; 8],
    /// By [`Swept::index`], at `f32`, or why the subject has none.
    arms32: Result<[Arm<D1F32>; 8], String>,
}

impl Seeded {
    /// The correct kernels at the generated switches, in [`Swept::ALL`] order.
    pub(crate) fn generated() -> Self {
        let switches = [
            K_F64.parts(),
            A_F64.parts(),
            B_F64.parts(),
            C_F64.parts(),
            D_F64.parts(),
            E_F64.parts(),
            COS_HALF_F64.parts(),
            R_F64.parts(),
        ];
        let switches32 = [
            K_F32.parts(),
            A_F32.parts(),
            B_F32.parts(),
            C_F32.parts(),
            D_F32.parts(),
            E_F32.parts(),
            COS_HALF_F32.parts(),
            R_F32.parts(),
        ];
        Self {
            name: "seeded:correct".to_string(),
            version: "generated".to_string(),
            defect: None,
            form: se3::W2::Product,
            arms: switches.map(Arm::generated),
            arms32: Ok(switches32.map(Arm::generated)),
        }
    }

    /// The correct kernels with `J_l`'s `W²` by the algebraic identity: the one-variable twin of
    /// 0037 (draft) for an assembly form, as `host::Twin` is for a library. Planted, so no plain
    /// run and no bar reads it.
    pub(crate) fn w2_identity() -> Self {
        Self {
            name: TWIN_W2.to_string(),
            version: "generated@w2-identity".to_string(),
            form: se3::W2::Identity,
            arms32: Err(
                "the `W²` form twin is a binary64 subject; the `@f32` strata are the \
                         coefficient ids' own"
                    .to_string(),
            ),
            ..Self::generated()
        }
    }

    /// The correct kernels, all six at one `candidate` over the corpus's series: how a test drives
    /// the harness at a candidate the sweep also scores, the `D12` prior included.
    #[cfg(test)]
    pub(crate) fn uniform(series: &Series<D1>, candidate: Candidate<D1>) -> Self {
        Self {
            name: "seeded:correct".to_string(),
            version: format!("{}terms-z{}", candidate.terms, candidate.switch_z.v),
            defect: None,
            form: se3::W2::Product,
            arms: Swept::ALL.map(|id| Arm {
                candidate,
                series: series.swept(id).to_vec(),
            }),
            arms32: Err("a uniform candidate is a binary64 test subject".to_string()),
        }
    }

    /// [`Self::generated`] with its `f32` kernels all eight at one `candidate` over the corpus's
    /// `f32` series: [`Self::uniform`] at `f32`, how a test drives the harness at a candidate the
    /// `f32` sweep also scores.
    #[cfg(test)]
    pub(crate) fn uniform_f32(series: &Series<D1F32>, candidate: Candidate<D1F32>) -> Self {
        let mut subject = Self::generated();
        subject.arms32 = Ok(Swept::ALL.map(|id| Arm {
            candidate,
            series: series.swept(id).to_vec(),
        }));
        subject
    }

    /// The generated kernels with `defect` planted.
    pub(crate) fn planted(defect: Defect) -> Self {
        let mut seeded = Self::generated();
        seeded.name = format!("seeded:{}", defect.name());
        seeded.defect = Some(defect);
        if defect == Defect::CTwoTermsEarly {
            let (terms, z) = C_PLANTED;
            seeded.arms[Coeff::C.index()].candidate = Candidate {
                terms,
                switch_z: D1::constant(z),
            };
            seeded.arms32 = Err(format!(
                "its `c` ({terms} terms below z = {z:e}) is a candidate of the binary64 sweep, \
                 which the f32 kernels do not plant"
            ));
        }
        seeded
    }

    /// The series length and switch this subject runs `c` with.
    pub(crate) fn candidate(&self, c: Coeff) -> Candidate<D1> {
        self.arms[c.index()].candidate
    }

    /// The series length, switch and terms (widened exactly) this subject runs `id` with at
    /// `precision`; `None` for a subject with no `f32` kernel.
    #[cfg(test)]
    pub(crate) fn switch(&self, id: Swept, at: Precision) -> Option<(usize, f64, Vec<f64>)> {
        fn parts<S: Real + Into<f64>>(a: &Arm<Dual<S, 1>>) -> (usize, f64, Vec<f64>) {
            let (cand, series) = (a.candidate, a.series.iter());
            (
                cand.terms,
                cand.switch_z.v.into(),
                series.map(|x| x.v.into()).collect(),
            )
        }
        match at {
            Precision::F64 => Some(parts(&self.arms[id.index()])),
            Precision::F32 => self.arms32.as_ref().ok().map(|a| parts(&a[id.index()])),
        }
    }

    /// The series this subject's `c` takes its terms from.
    #[cfg(test)]
    pub(crate) fn series(&self, c: Coeff) -> &[D1] {
        &self.arms[c.index()].series
    }

    /// The generated kernels of all six coefficients, as `se3` reads them: the correct subject's,
    /// a coefficient defect being planted by `kernel`, which no `sen3_*` id reaches (`supports`).
    fn kernels(&self) -> se3::Kernels<'_, D1> {
        kernels_of(&self.arms)
    }

    pub(crate) fn registered(self) -> Registered {
        Registered {
            version: self.version.clone(),
            version_f32: self.version.clone(),
            no_f32: self.arms32.as_ref().err().cloned(),
            planted: self.defect.is_some() || self.form != se3::W2::Product,
            subject: Box::new(self),
        }
    }
}

impl Subject for Seeded {
    fn name(&self) -> &str {
        &self.name
    }

    fn supports(&self, fn_id: &str) -> bool {
        let coefficient = Coeff::of_fn(fn_id).is_some();
        let sen3 = se3::parse(fn_id);
        match self.defect {
            None => {
                Swept::of_fn(fn_id).is_some()
                    || sen3.is_some()
                    || matches!(fn_id, "so3_exp" | "so3_log" | "so3_jr" | "so3_jl")
            }
            Some(Defect::LogAcos | Defect::LogNoFlip) => fn_id == "so3_log",
            Some(Defect::Se3ExpTranslationFirst) => sen3.is_some_and(|(op, _)| op == se3::Op::Exp),
            Some(Defect::QMinusHalf) => sen3.is_some_and(|(op, _)| op != se3::Op::Exp),
            Some(_) => coefficient,
        }
    }

    /// An empty answer is an error in the harness, never a score. An id `supports` refuses is
    /// answered with nothing, so a planted subject never answers as the correct kernel would; at
    /// `f32` only the coefficient ids are answered (`0016`: no other id has an `@f32` stratum), only
    /// from an input that is exactly a binary32, and only by a subject with an `f32` kernel
    /// (`Registered::no_f32`, which the harness reads first, for its reason).
    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        if !self.supports(fn_id) {
            return Output::new();
        }
        match (Swept::of_fn(fn_id), precision) {
            (Some(id), Precision::F64) => answer(self.defect, id, &self.arms[id.index()], record),
            (Some(id), Precision::F32) => match &self.arms32 {
                Ok(arms) => answer(self.defect, id, &arms[id.index()], record),
                Err(_) => Output::new(),
            },
            (None, Precision::F64) if fn_id == "so3_exp" => self.exp(record),
            (None, Precision::F64) if fn_id == "so3_log" => self.log(record),
            (None, Precision::F64) if matches!(fn_id, "so3_jr" | "so3_jl") => {
                self.so3_jac(fn_id, record)
            }
            (None, Precision::F64) => self.sen3(fn_id, record),
            (None, Precision::F32) => Output::new(),
        }
    }
}

/// The kernel of `id` at `x`, the correct one or `defect`'s if it plants this coefficient.
fn kernel<S: Real>(defect: Option<Defect>, id: Swept, x: Input<S>, arm: &Arm<S>) -> S {
    let Arm { candidate, series } = arm;
    match (defect, id) {
        (Some(Defect::BNoSeries), Swept::Coeff(Coeff::B)) => b_no_series(x.z),
        (Some(Defect::KSqrtUnsafe), Swept::Coeff(Coeff::K)) => {
            k_sqrt_unsafe(x.z, *candidate, series)
        }
        _ => evaluate(id, x, *candidate, series),
    }
}

/// The adapter of both precisions: the arguments at the scalar `S` (`input`), the value and
/// `d/dz` of one `Dual<S, 1>` evaluation, widened exactly to binary64. Nothing when the record
/// holds no usable input.
fn answer<S: Real + Into<f64>>(
    defect: Option<Defect>,
    id: Swept,
    arm: &Arm<Dual<S, 1>>,
    record: &Record,
) -> Output {
    let Some(x) = input::<S>(id, record) else {
        return Output::new();
    };
    let r = kernel(defect, id, x.seed(), arm);
    Output::from([
        ("value".to_string(), vec![r.v.into()]),
        ("d_branch".to_string(), r.d.map(Into::into).to_vec()),
    ])
}

/// `so3_exp` at the scalar `S`: `q` of one evaluation on the generated `k` of `arm`.
fn exp_at<S: Real + Into<f64> + From<f64>>(record: &Record, arm: &Arm<Dual<S, 1>>) -> Output {
    let Some(&[x, y, z]) = record.input("phi").and_then(|p| p.first_chunk::<3>()) else {
        return Output::new();
    };
    let Arm { candidate, series } = arm;
    let phi = [x, y, z].map(|c| Dual::constant(S::from(c)));
    let q: [f64; 4] = so3::exp(phi, *candidate, series).map(|c| c.v.into());
    Output::from([("q".to_string(), q.to_vec())])
}

/// `so3_log` at the scalar `S`: `phi` of `so3::log`, the quaternion `atan2` (D5).
fn log_at<S: Real + Into<f64> + From<f64>>(record: &Record) -> Output {
    let Some(&q) = record.input("q").and_then(|q| q.first_chunk::<4>()) else {
        return Output::new();
    };
    let q: [Dual<S, 1>; 4] = q.map(|c| Dual::constant(S::from(c)));
    let phi: [f64; 3] = so3::log(q).map(|c| c.v.into());
    Output::from([("phi".to_string(), phi.to_vec())])
}

/// `so3_j{r,l}` at the scalar `S`: the dense column-major `J` of the value of one `Dual<S, 1>`
/// evaluation, which is `se3`'s rotation block alone (`se3::so3_jacobian`).
fn so3_jac_at<S: Real + Into<f64> + From<f64>>(
    fn_id: &str,
    record: &Record,
    kernels: &se3::Kernels<'_, Dual<S, 1>>,
    form: se3::W2,
) -> Output {
    let Some(&[x, y, z]) = record.input("phi").and_then(|p| p.first_chunk::<3>()) else {
        return Output::new();
    };
    // The side is read from the id by name, both spelled: `fn_id == "so3_jr"` would answer every
    // other id with `J_l` under its own name, and an id this subject does not compute must answer
    // with nothing (`Subject::eval`).
    let right = match fn_id {
        "so3_jr" => true,
        "so3_jl" => false,
        _ => return Output::new(),
    };
    let phi = [x, y, z].map(|c| Dual::constant(S::from(c)));
    let j = se3::so3_jacobian(phi, right, kernels, form);
    let j: Vec<f64> = j.iter().map(|c| c.v.into()).collect();
    Output::from([("J".to_string(), j)])
}

/// `sen3_{exp,jr,jl}_n<N>` at the scalar `S`: `q` and `x`, or the dense `J`, of the value of one
/// `Dual<S, 1>` evaluation. `order` and `half` carry the two planted SE(3) defects; the correct
/// kernel and the twin pass `NUMERICS.md`'s values.
fn sen3_at<S: Real + Into<f64> + From<f64>>(
    fn_id: &str,
    record: &Record,
    kernels: &se3::Kernels<'_, Dual<S, 1>>,
    order: se3::Order,
    half: f64,
    form: se3::W2,
) -> Output {
    let (Some((op, n)), Some(tau)) = (se3::parse(fn_id), record.input("tau")) else {
        return Output::new();
    };
    if tau.len() != 3 + 3 * n {
        return Output::new();
    }
    let tau: Vec<Dual<S, 1>> = tau.iter().map(|&t| Dual::constant(S::from(t))).collect();
    let value = |v: &[Dual<S, 1>]| v.iter().map(|c| c.v.into()).collect::<Vec<f64>>();
    match op {
        se3::Op::Exp => {
            let Some((q, x)) = se3::exp(&tau, kernels, order, form) else {
                return Output::new();
            };
            Output::from([("q".to_string(), value(&q)), ("x".to_string(), value(&x))])
        }
        se3::Op::Jr | se3::Op::Jl => {
            let half = Dual::constant(S::from(half));
            let j = se3::jacobian(&tau, op == se3::Op::Jr, kernels, half, form);
            j.map_or_else(Output::new, |j| {
                Output::from([("J".to_string(), value(&j))])
            })
        }
    }
}

/// The six coefficient kernels of `arms`, as `se3` reads them: [`Swept`] orders the coefficients
/// first, so the first six are `Coeff::ALL` by [`Coeff::index`].
fn kernels_of<S: Real>(arms: &[Arm<Dual<S, 1>>; 8]) -> se3::Kernels<'_, Dual<S, 1>> {
    se3::Kernels {
        arms: std::array::from_fn(|i| {
            let Arm { candidate, series } = &arms[i];
            (*candidate, series.as_slice())
        }),
    }
}

impl Seeded {
    fn exp(&self, record: &Record) -> Output {
        exp_at(record, &self.arms[Coeff::K.index()])
    }

    /// `so3_j{r,l}`, with this subject's `W²` form. No planted defect answers these ids
    /// (`supports`): `Q`'s `½` and the tangent order are not in the rotation block.
    fn so3_jac(&self, fn_id: &str, record: &Record) -> Output {
        so3_jac_at(fn_id, record, &self.kernels(), self.form)
    }

    /// `sen3_{exp,jr,jl}_n<N>`, with this subject's two planted SE(3) defects as the arguments of
    /// [`sen3_at`]: the tangent order of `Exp`, and the sign of `Q`'s `½`.
    fn sen3(&self, fn_id: &str, record: &Record) -> Output {
        let order = match self.defect {
            Some(Defect::Se3ExpTranslationFirst) => se3::Order::TranslationFirst,
            _ => se3::Order::RotationFirst,
        };
        let half = if self.defect == Some(Defect::QMinusHalf) {
            -0.5
        } else {
            0.5
        };
        sen3_at(fn_id, record, &self.kernels(), order, half, self.form)
    }

    fn log(&self, record: &Record) -> Output {
        let Some(&q) = record.input("q").and_then(|q| q.first_chunk::<4>()) else {
            return Output::new();
        };
        let phi = match self.defect {
            Some(Defect::LogAcos) => so3::log_acos(q),
            Some(Defect::LogNoFlip) => so3::log_no_flip(q.map(D1::constant)).map(|c| c.v),
            _ => return log_at::<f64>(record),
        };
        Output::from([("phi".to_string(), phi.to_vec())])
    }
}

/// The subjects `just conformance` knows: the correct kernels, every defect, and the host-`std`
/// twin ([`host`]). Everything but the correct kernels is planted: a run that names no subject
/// skips them.
pub(crate) fn registry() -> Vec<Registered> {
    let mut all = vec![Seeded::generated().registered()];
    let defects = Defect::COEFFICIENT
        .into_iter()
        .chain(Defect::LOG)
        .chain(Defect::SE3);
    all.extend(defects.map(|d| Seeded::planted(d).registered()));
    all.push(host::Twin::registered());
    all.push(Seeded::w2_identity().registered());
    all
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::corpus_dir;
    use crate::conformance::testkit::record;

    fn subject(defect: Option<Defect>) -> Seeded {
        match defect {
            None => Seeded::generated(),
            Some(d) => Seeded::planted(d),
        }
    }

    /// `d/dz` of the 16-term series, `Σ j s_j z^(j-1)`, in plain `f64`: independent of the kernels'
    /// exact arms and of `Dual`.
    fn series_derivative(s: &[f64], z: f64) -> f64 {
        let step = |acc: f64, (j, &s): (usize, &f64)| j as f64 * s + z * acc;
        s.iter().enumerate().skip(1).rev().fold(0.0, step)
    }

    #[test]
    fn a_subject_answers_the_value_and_the_derivative_in_the_branch_variable() -> Result<(), String>
    {
        let (s, theta) = (subject(None), 0.5);
        let rec = record(&[("theta", &[theta])], &[])?;
        let out = s.eval("coeff_a", &rec, Precision::F64);
        assert_eq!(out.keys().collect::<Vec<_>>(), ["d_branch", "value"]);
        // a = (1 - cos(theta)) / theta^2 and d/dz a = (theta sin(theta) - 2 (1 - cos(theta))) / (2 theta^4).
        let (sin, cos) = theta.sin_cos();
        let z = theta * theta;
        assert!((out["value"][0] - (1.0 - cos) / z).abs() < 1e-15);
        assert!(
            (out["d_branch"][0] - (theta * sin - 2.0 * (1.0 - cos)) / (2.0 * z * z)).abs() < 1e-13
        );
        Ok(())
    }

    #[test]
    fn the_derivative_is_in_z_not_in_theta_on_both_sides_of_the_switch() -> Result<(), String> {
        // At theta = 0.5, 2 theta = 1 and d/dz = d/dtheta: any other theta tells them apart.
        let (s, series) = (subject(None), Series::<f64>::load(&corpus_dir()?)?);
        for c in Coeff::ALL {
            let fn_id = format!("coeff_{}", c.name());
            for theta in [0.02, 0.06, 0.3, 0.7, 2.0, 3.0] {
                let rec = record(&[("theta", &[theta])], &[])?;
                let out = s.eval(&fn_id, &rec, Precision::F64);
                let (z, s16) = (theta * theta, series.of(c));
                let want = series_derivative(s16, z);
                let got = out["d_branch"][0];
                assert!(
                    (got - want).abs() <= 1e-7 * want.abs(),
                    "{fn_id} at {theta}: {got} vs {want}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn a_subject_supports_the_eight_coefficients_and_answers_nothing_it_cannot(
    ) -> Result<(), String> {
        let s = subject(None);
        for id in Swept::ALL {
            assert!(s.supports(&format!("coeff_{}", id.name())), "{id:?}");
        }
        for id in ["coeff_", "coeff_f", "so3_act", "coeff_series"] {
            assert!(!s.supports(id), "{id}");
        }
        // `cos θ/2` and `r` belong to the correct kernel; a planted defect runs the six.
        let planted = subject(Some(Defect::KSqrtUnsafe));
        assert!(planted.supports("coeff_e") && !planted.supports("coeff_r"));
        assert!(!planted.supports("coeff_cos_half"));
        // `r` reads `n` and `w`, the others `theta`.
        let rec = record(&[("theta", &[0.5])], &[])?;
        assert!(s.eval("coeff_r", &rec, Precision::F64).is_empty());
        let nw = record(&[("n", &[0.5]), ("w", &[1.0])], &[])?;
        assert!(s.eval("coeff_k", &nw, Precision::F64).is_empty());
        assert_eq!(s.eval("coeff_r", &nw, Precision::F64).len(), 2);
        for precision in [Precision::F64, Precision::F32] {
            let none = record(&[], &[])?;
            assert!(s.eval("coeff_k", &none, precision).is_empty());
            assert_eq!(s.eval("coeff_k", &rec, precision).len(), 2);
        }
        // An `f32` subject receives a binary32 or nothing: 0.1 is no binary32, and is not rounded.
        let tenth = record(&[("theta", &[0.1])], &[])?;
        assert!(s.eval("coeff_k", &tenth, Precision::F32).is_empty());
        assert_eq!(s.eval("coeff_k", &tenth, Precision::F64).len(), 2);
        Ok(())
    }

    #[test]
    fn at_f32_the_adapter_answers_the_binary32_kernel_widened_exactly() -> Result<(), String> {
        let s = subject(None);
        let dual = |x: f64| Dual::<f32, 1>::constant(x as f32);
        for id in Swept::ALL {
            let (_, below, terms) = s.switch(id, Precision::F32).ok_or("no f32 kernel")?;
            let cand = Candidate {
                terms: terms.len(),
                switch_z: dual(below),
            };
            let series: Vec<_> = terms.into_iter().map(dual).collect();
            for theta in [0.0f32, 1e-30, 1e-3, 0.05, 0.5, 2.0, 3.0] {
                // `n = θ` and `w = 1` for `r`, `θ` for the rest.
                let t = f64::from(theta);
                let rec = record(&[("theta", &[t]), ("n", &[t]), ("w", &[1.0])], &[])?;
                let out = s.eval(&format!("coeff_{}", id.name()), &rec, Precision::F32);
                let x = Input {
                    z: Dual::variable(theta * theta, 0),
                    w: dual(1.0),
                };
                let want = evaluate(id, x, cand, &series);
                let got = (out["value"][0].to_bits(), out["d_branch"][0].to_bits());
                let want = (f64::from(want.v).to_bits(), f64::from(want.d[0]).to_bits());
                assert_eq!(got, want, "{id:?} at {theta}");
            }
        }
        Ok(())
    }

    #[test]
    fn the_registry_holds_the_correct_kernel_and_every_defect_planted() {
        let all = registry();
        let names: Vec<(&str, bool)> = all.iter().map(|r| (r.subject.name(), r.planted)).collect();
        let want = [
            ("seeded:correct", false),
            ("seeded:b-no-series", true),
            ("seeded:k-sqrt-unsafe", true),
            ("seeded:c-two-terms-1e-8", true),
            ("seeded:log-acos", true),
            ("seeded:log-no-flip", true),
            ("seeded:se3-exp-translation-first", true),
            ("seeded:q-minus-half", true),
            ("seeded:host-std", true),
            ("seeded:w2-identity", true),
        ];
        assert_eq!(names, want);
        // Every subject runs the generated kernels; a twin says in its version what it changed,
        // since that is the only thing it changes (`0037`, draft; `0038`, draft).
        let kernels: Vec<&str> = all.iter().map(|r| r.version.as_str()).collect();
        assert_eq!(kernels[..8], ["generated"; 8]);
        assert_eq!(
            kernels[8..],
            ["generated@host-std", "generated@w2-identity"]
        );
        assert!(all
            .iter()
            .all(|r| r.version_at(Precision::F32) == r.version));
        // Three subjects have no `f32` kernel: the planted `c`, a candidate of the binary64 sweep,
        // and the two twins, whose subjects are binary64 readings.
        let none = all.iter().filter(|r| r.no_f32.is_some());
        let none: Vec<&str> = none.map(|r| r.subject.name()).collect();
        assert_eq!(
            none,
            [
                "seeded:c-two-terms-1e-8",
                "seeded:host-std",
                "seeded:w2-identity"
            ]
        );
    }

    #[test]
    fn a_subject_the_f32_kernels_do_not_model_answers_nothing_at_f32() -> Result<(), String> {
        let rec = record(&[("theta", &[0.5])], &[])?;
        let c = Seeded::planted(Defect::CTwoTermsEarly);
        for id in ["coeff_c", "coeff_b"] {
            assert!(c.eval(id, &rec, Precision::F32).is_empty(), "{id}");
            assert_eq!(c.eval(id, &rec, Precision::F64).len(), 2, "{id}");
        }
        // A uniform candidate is binary64's too; the correct kernel has one.
        let series = Series::<D1>::load(&corpus_dir()?)?;
        let uniform = Seeded::uniform(&series, d12(&series)?);
        assert!(uniform.eval("coeff_b", &rec, Precision::F32).is_empty());
        assert_eq!(subject(None).eval("coeff_b", &rec, Precision::F32).len(), 2);
        // A defect it does model is planted there: `b` without its series is another answer than
        // the correct kernel's at a small binary32 `θ`, where the exact arm cancels.
        let small = record(&[("theta", &[f64::from(1e-3f32)])], &[])?;
        let value = |d| subject(d).eval("coeff_b", &small, Precision::F32)["value"][0];
        let (good, bad) = (value(None), value(Some(Defect::BNoSeries)));
        assert_ne!(good.to_bits(), bad.to_bits());
        Ok(())
    }

    #[test]
    fn a_defect_changes_one_coefficient_of_the_generated_kernel_and_no_other() {
        let generated = Seeded::generated();
        let planted = Seeded::planted(Defect::CTwoTermsEarly);
        for c in Coeff::ALL {
            let (a, b) = (generated.candidate(c), planted.candidate(c));
            let same = (a.terms, a.switch_z.v.to_bits()) == (b.terms, b.switch_z.v.to_bits());
            assert_eq!(same, c != Coeff::C, "{c:?}");
        }
        // Two terms below 1e-16: §10's 1e-8 as θ, and the grid's first point.
        let c = planted.candidate(Coeff::C);
        assert_eq!((c.terms, c.switch_z.v.to_bits()), (2, 1e-16f64.to_bits()));
        assert_eq!(C_PLANTED, (2, 1e-16));
        // Its series is the generated one, of which two terms are used: the same length the
        // correct kernel holds, which is the sweep's term cap and not a number to type here.
        let i = Coeff::C.index();
        assert_eq!(planted.arms[i].series.len(), generated.arms[i].series.len());
    }

    #[test]
    fn so3_ids_belong_to_the_correct_subject_and_the_log_defects_take_so3_log_only(
    ) -> Result<(), String> {
        // `so3_act` is the control: the corpus holds it and this subject has no action.
        let ids = [
            "so3_exp", "so3_log", "coeff_k", "so3_act", "so3_jr", "so3_jl",
        ];
        let supported = |s: &Seeded| ids.map(|id| s.supports(id));
        assert_eq!(
            supported(&subject(None)),
            [true, true, true, false, true, true]
        );
        for d in Defect::LOG {
            assert_eq!(
                supported(&subject(Some(d))),
                [false, true, false, false, false, false]
            );
        }
        for d in Defect::COEFFICIENT {
            assert_eq!(
                supported(&subject(Some(d))),
                [false, false, true, false, false, false]
            );
        }
        let s = subject(None);
        let exp = s.eval(
            "so3_exp",
            &record(&[("phi", &[0.0, 0.0, 2.0])], &[])?,
            Precision::F64,
        );
        let log = s.eval(
            "so3_log",
            &record(&[("q", &[1.0, 0.0, 0.0, 0.0])], &[])?,
            Precision::F64,
        );
        assert_eq!(exp.keys().collect::<Vec<_>>(), ["q"]);
        assert_eq!(log.keys().collect::<Vec<_>>(), ["phi"]);
        assert_eq!((exp["q"].len(), log["phi"].len()), (4, 3));
        // A record without its input, or at `f32`, is answered with nothing.
        assert!(s
            .eval("so3_exp", &record(&[], &[])?, Precision::F64)
            .is_empty());
        let unit = record(&[("q", &[1.0, 0.0, 0.0, 0.0])], &[])?;
        assert!(s.eval("so3_log", &unit, Precision::F32).is_empty());
        // `so3_j{r,l}`: the dense `3 x 3`, and `J_l(φ)` the transpose of `J_r(φ)` to the bit. The
        // sign of the tangent is exact and `W²` is symmetric in the same products, so the two
        // sides are one program read twice and the corpus's column-major order is pinned here.
        let turn = record(&[("phi", &[0.3, -1.1, 0.7])], &[])?;
        let (jr, jl) = (
            s.eval("so3_jr", &turn, Precision::F64),
            s.eval("so3_jl", &turn, Precision::F64),
        );
        assert_eq!((jr["J"].len(), jl["J"].len()), (9, 9));
        let transposed: Vec<f64> = (0..3)
            .flat_map(|c| (0..3).map(move |r| (r, c)))
            .map(|(r, c)| jr["J"][r * 3 + c])
            .collect();
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&jl["J"]), bits(&transposed));
        assert!(s
            .eval("so3_jr", &record(&[], &[])?, Precision::F64)
            .is_empty());
        assert!(s.eval("so3_jl", &turn, Precision::F32).is_empty());
        Ok(())
    }

    #[test]
    fn the_sen3_ids_belong_to_the_correct_subject_and_each_se3_defect_answers_only_its_own(
    ) -> Result<(), String> {
        let ids = [
            "sen3_exp_n1",
            "sen3_exp_n2",
            "sen3_jr_n3",
            "sen3_jl_n1",
            "sen3_exp_n4",
            "sen3_log_n1",
            "sen3_jr_inv_n1",
        ];
        let supported = |s: &Seeded| ids.map(|id| s.supports(id));
        let (t, f) = (true, false);
        assert_eq!(supported(&subject(None)), [t, t, t, t, f, f, f]);
        let first = subject(Some(Defect::Se3ExpTranslationFirst));
        assert_eq!(supported(&first), [t, t, f, f, f, f, f]);
        assert_eq!(
            supported(&subject(Some(Defect::QMinusHalf))),
            [f, f, t, t, f, f, f]
        );
        for d in Defect::COEFFICIENT.into_iter().chain(Defect::LOG) {
            assert_eq!(supported(&subject(Some(d))), [f; 7]);
        }
        let s = subject(None);
        let answer = |id: &str, n: usize, precision| -> Result<Output, String> {
            let tau = vec![0.1; 3 + 3 * n];
            Ok(s.eval(id, &record(&[("tau", &tau)], &[])?, precision))
        };
        let exp = answer("sen3_exp_n2", 2, Precision::F64)?;
        assert_eq!(exp.keys().collect::<Vec<_>>(), ["q", "x"]);
        assert_eq!((exp["q"].len(), exp["x"].len()), (4, 6));
        let jr = answer("sen3_jr_n3", 3, Precision::F64)?;
        assert_eq!(jr.keys().collect::<Vec<_>>(), ["J"]);
        assert_eq!(jr["J"].len(), 144);
        // A tangent of the wrong length, no tangent, and `f32` are answered with nothing.
        assert!(answer("sen3_exp_n1", 2, Precision::F64)?.is_empty());
        assert!(answer("sen3_jl_n1", 1, Precision::F32)?.is_empty());
        let none = record(&[], &[])?;
        assert!(s.eval("sen3_jl_n1", &none, Precision::F64).is_empty());
        Ok(())
    }

    #[test]
    fn a_planted_subject_answers_nothing_it_does_not_support() -> Result<(), String> {
        // Every input a `theta`, `phi`, `q` or `tau` id reads, so only `supports` can refuse.
        let rec = record(
            &[
                ("theta", &[0.5]),
                ("phi", &[0.1, 0.2, 0.3]),
                ("q", &[1.0, 0.0, 0.0, 0.0]),
                ("tau", &[0.1; 6]),
            ],
            &[],
        )?;
        let ids = [
            "coeff_k",
            "so3_exp",
            "so3_log",
            "so3_jr",
            "sen3_exp_n1",
            "sen3_jr_n1",
        ];
        let defects = Defect::COEFFICIENT.into_iter().chain(Defect::LOG);
        for d in defects.chain(Defect::SE3) {
            let s = subject(Some(d));
            let answered = ids.map(|id| !s.eval(id, &rec, Precision::F64).is_empty());
            assert_eq!(answered, ids.map(|id| s.supports(id)), "{d:?}");
        }
        Ok(())
    }
}
