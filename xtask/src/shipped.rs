//! The `helicoid` conformance subject (`docs/PHASE1.md` §5) over the eight `coeff_*` ids: the
//! shipped groups of `helicoid::coeffs` through the hidden `__sweep` feature, so the code scored is
//! the code shipped (`0004` item 4), at both precisions on `Dual<S, 1>` seeded at the branch variable.
//!
//! The arguments are the seeded subject's (`seeded::input`: `z = fl(θ·θ)`, `n²` and `w` for `r`,
//! formed at the precision, `@f32` inputs cast losslessly) and the answer its fields, `value` and
//! `d_branch`. Each id reads the group that holds it at the call sites of `docs/PHASE3.md` §3:
//! `k` and `cos θ/2` from `exp_coeffs`, `a` and `b` from `jr_coeffs`, `c` from `jr_inv_coeff`, `d`
//! and `e` from `q_coeffs`, `r` from `log_ratio`. `b` also lives in `q_coeffs`; the two agree to
//! the bit at both precisions (`b_is_the_same_in_both_groups_that_hold_it` and its `_at_f32` twin,
//! sweep tests), so one is scored. The group ids, the geodesics and `0056`'s routines follow, each
//! family an enum; `0056`'s are answered at both precisions.
//!
//! What these tests do not cover: the boundary (`s = switch`, `w = 0`) and the safe-argument
//! pattern. A scalar mask makes `S::branch` lazy, so the exact arm never sees an unsafe argument
//! and a finite score is no evidence for it; the `helicoid` crate's lane and boundary tests are.
//!
//! `subject_version` is the workspace version: there is no other version of the shipped kernel.

use helicoid::__sweep as k;
use helicoid::{Jac, LieGroup, Quat, SEn3, SEn3Jac, SEn3Tangent, Tangent, SO3};
use helicoid_linalg::{
    chol, chol_solve, eig3, solve_cubic, Dual, Mat3, Matrix, Precision, Real, StridedMut, Vec3,
    Vector,
};

use crate::conformance::corpus::{exact_f32, Record};
use crate::conformance::subject::{Output, Registered, Subject};
use crate::seeded::{input, Coeff, Host, Input, Swept};

/// The group of `id` at `x`, as the shipped kernel evaluates it.
pub(crate) fn shipped<S: Real>(id: Swept, x: Input<Dual<S, 1>>) -> Dual<S, 1> {
    match id {
        Swept::Coeff(Coeff::K) => k::exp_coeffs(x.z).0,
        Swept::CosHalf => k::exp_coeffs(x.z).1,
        Swept::Coeff(Coeff::A) => k::jr_coeffs(x.z).0,
        Swept::Coeff(Coeff::B) => k::jr_coeffs(x.z).1,
        Swept::Coeff(Coeff::C) => k::jr_inv_coeff(x.z),
        Swept::Coeff(Coeff::D) => k::q_coeffs(x.z).1,
        Swept::Coeff(Coeff::E) => k::q_coeffs(x.z).2,
        Swept::R => k::log_ratio(x.z, x.w),
    }
}

/// An `so3_*` id this subject answers.
///
/// An enum with an exhaustive `answer`, not a `&str` match with a `_` arm. The `&str` form had two
/// levels of catch-all: an id added to the supported list without its own arm fell through `eval`
/// into the Jacobian helper and through *that* into `jl_inv`, so the harness would write a scored,
/// committed conformance row labelled with the new id and holding `J_l⁻¹` — no error, no empty
/// output, nothing to notice. Here the same mistake is a non-exhaustive `match`.
#[derive(Clone, Copy)]
enum So3 {
    Exp,
    Log,
    Act,
    FromMatrix,
    Jr,
    Jl,
    JrInv,
    JlInv,
}

impl So3 {
    /// Every `so3_*` id the corpus holds (`PHASE1.md` §4.3), with its name.
    const ALL: [(&'static str, So3); 8] = [
        ("so3_exp", So3::Exp),
        ("so3_log", So3::Log),
        ("so3_act", So3::Act),
        ("so3_from_matrix", So3::FromMatrix),
        ("so3_jr", So3::Jr),
        ("so3_jl", So3::Jl),
        ("so3_jr_inv", So3::JrInv),
        ("so3_jl_inv", So3::JlInv),
    ];

    fn of_fn(fn_id: &str) -> Option<Self> {
        Self::ALL
            .iter()
            .find(|(name, _)| *name == fn_id)
            .map(|&(_, id)| id)
    }

    /// The shipped answer at `S`; nothing when the record holds no usable input.
    ///
    /// `f64` only, by the `From<f64>` bound: no vector id has an `@f32` stratum until a record
    /// extends `0016`, and that bound — which `f32` does not implement — makes the `f32`
    /// instantiation statically unreachable, as it does for the seeded subject's `so3_*` arms.
    fn answer<S: Real + Into<f64> + From<f64>>(self, record: &Record) -> Output {
        match self {
            So3::Exp => {
                let Some(phi) = vec3_of::<S>(record, "phi") else {
                    return Output::new();
                };
                // Through `read_dense`, so the dense order the corpus writes is the one §1 states
                // and not this file's reading of it (`API.md` R3).
                let tau = <SO3<S> as LieGroup<S>>::Tangent::read_dense(&phi.0);
                let q = SO3::<S>::exp(&tau).quat();
                quat_out(&q)
            }
            So3::Log => {
                let Some(q) = quat_of::<S>(record, "q") else {
                    return Output::new();
                };
                let mut phi = [S::zero(); 3];
                SO3::from_quat_unchecked(q).log().write_dense(&mut phi);
                Output::from([("phi".to_string(), phi.map(Into::into).to_vec())])
            }
            So3::Act => {
                let (Some(q), Some(p)) = (quat_of::<S>(record, "q"), vec3_of::<S>(record, "p"))
                else {
                    return Output::new();
                };
                // §3.3 as written, so a non-unit `q` gives the scaled rotation §1 defines rather
                // than its normalization; `SO3::act`'s *Domain* has the measured gap.
                let rp = SO3::from_quat_unchecked(q).act(p);
                Output::from([("Rp".to_string(), rp.0.map(Into::into).to_vec())])
            }
            So3::FromMatrix => {
                let Some(&m) = record.input("R").and_then(|r| r.first_chunk::<9>()) else {
                    return Output::new();
                };
                // Column-major, as the record's sibling `shape` says (`PHASE1.md` §4.3).
                let cols = core::array::from_fn(|c| {
                    Vector(core::array::from_fn(|r| S::from(m[c * 3 + r])))
                });
                quat_out(&SO3::<S>::from_matrix(&Matrix::from_cols(cols)).quat())
            }
            So3::Jr | So3::Jl | So3::JrInv | So3::JlInv => {
                let Some(phi) = vec3_of::<S>(record, "phi") else {
                    return Output::new();
                };
                let tau = <SO3<S> as LieGroup<S>>::Tangent::read_dense(&phi.0);
                let j = match self {
                    So3::Jr => SO3::<S>::jr(&tau),
                    So3::Jl => SO3::<S>::jl(&tau),
                    So3::JrInv => SO3::<S>::jr_inv(&tau),
                    _ => SO3::<S>::jl_inv(&tau),
                };
                jac_out(&j)
            }
        }
    }
}

/// The quaternion a record holds under `key`, built from its fields.
///
/// Not through `Quat::from_wxyz_unchecked`, whose `debug_assert!` would reject the `q:nonunit` and
/// `q:w0` strata that `so3_log` and `so3_act` carry on purpose, and would make the scale
/// invariance §3.2 states unreachable from the corpus.
fn quat_of<S: Real + From<f64>>(record: &Record, key: &str) -> Option<Quat<S>> {
    let &[w, x, y, z] = record.input(key)?.first_chunk::<4>()?;
    Some(Quat {
        w: S::from(w),
        x: S::from(x),
        y: S::from(y),
        z: S::from(z),
    })
}

/// A `3`-vector a record holds under `key`.
fn vec3_of<S: Real + From<f64>>(record: &Record, key: &str) -> Option<Vector<S, 3>> {
    let &[a, b, c] = record.input(key)?.first_chunk::<3>()?;
    Some(Vector([a, b, c].map(S::from)))
}

/// `q` as the corpus holds it, `w` first (`NUMERICS.md` §1).
fn quat_out<S: Real + Into<f64>>(q: &Quat<S>) -> Output {
    Output::from([(
        "q".to_string(),
        [q.w, q.x, q.y, q.z].map(Into::into).to_vec(),
    )])
}

/// `J` as the corpus holds it: the dense `3 x 3`, column-major (`PHASE1.md` §4.3), through
/// `Jac::write_dense` -- the shipped path a consumer takes, structural zeros included.
fn jac_out<S: Real + Into<f64>>(j: &Mat3<S>) -> Output {
    let mut buf = [S::zero(); 9];
    Jac::<S, <SO3<S> as LieGroup<S>>::Tangent>::write_dense(
        j,
        &mut StridedMut::col_major(&mut buf, 3, 3),
    );
    Output::from([("J".to_string(), buf.map(Into::into).to_vec())])
}

/// The subject's answer at `S`: the value and `d/dz`, widened exactly to binary64; nothing when the
/// record holds no usable input.
fn answer<S: Real + Into<f64>>(id: Swept, record: &Record) -> Output {
    let Some(x) = input::<S>(id, record) else {
        return Output::new();
    };
    let r = shipped(id, x.seed());
    Output::from([
        ("value".to_string(), vec![r.v.into()]),
        ("d_branch".to_string(), r.d.map(Into::into).to_vec()),
    ])
}

/// What a `sen3_*` id computes (`docs/PHASE3.md` §5).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Sen3 {
    Exp,
    Log,
    Ad,
    Jr,
    Jl,
    JrInv,
    JlInv,
}

impl Sen3 {
    /// `(op, N)` of `sen3_<op>_n<N>`; `None` for every other id. The same spelling the seeded
    /// subject parses (`seeded::se3::parse`), extended by the ids only the library answers.
    fn of_fn(fn_id: &str) -> Option<(Self, usize)> {
        let (op, n) = fn_id.strip_prefix("sen3_")?.rsplit_once("_n")?;
        let op = match op {
            "exp" => Self::Exp,
            "log" => Self::Log,
            "ad" => Self::Ad,
            "jr" => Self::Jr,
            "jl" => Self::Jl,
            "jr_inv" => Self::JrInv,
            "jl_inv" => Self::JlInv,
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
}

/// The group a record holds under the keys `(q, x)`, at `N`.
///
/// The keys are parameters because the geodesic ids hold two poses per record, `(q0, x0)` and
/// `(q1, x1)` (`PHASE1.md` §4.3), and one builder for all three is one place for the
/// `from_quat_unchecked` reading below to live.
fn sen3_of<S: Real + From<f64>, const N: usize>(
    record: &Record,
    qk: &str,
    xk: &str,
) -> Option<SEn3<S, N>> {
    let q = quat_of::<S>(record, qk)?;
    let x = record.input(xk).filter(|x| x.len() == 3 * N)?;
    let cols = core::array::from_fn(|i| Vector(core::array::from_fn(|r| S::from(x[3 * i + r]))));
    // Not through `Quat::from_wxyz_unchecked`'s assert, as `quat_of` says: the `q:*` strata carry
    // quaternions that are unit only to rounding on purpose.
    Some(SEn3::from_parts(SO3::from_quat_unchecked(q), cols))
}

/// The dense `(3 + 3N)`-square matrix of a `SEn3Jac` under `field`, column-major as the corpus
/// holds it (`PHASE1.md` §4.3), through `Jac::write_dense` -- the shipped path, structural zeros
/// included.
///
/// `field` is a parameter because the corpus names the adjoint `Ad` and the Jacobians `J`: building
/// one key and renaming it would answer `{"Ad": []}` on a mismatch, which is non-empty, so the
/// harness would score a committed row holding no numbers -- the silent-mislabel failure this
/// module's `So3` note is about.
fn sen3_jac_out<S: Real + Into<f64>, const N: usize>(field: &str, j: &SEn3Jac<S, N>) -> Output {
    let d = 3 + 3 * N;
    let mut buf = vec![S::zero(); d * d];
    Jac::<S, SEn3Tangent<S, N>>::write_dense(j, &mut StridedMut::col_major(&mut buf, d, d));
    Output::from([(
        field.to_string(),
        buf.into_iter().map(Into::into).collect::<Vec<f64>>(),
    )])
}

/// The shipped answer of `op` at `N`, or nothing when the record holds no usable input.
fn sen3_answer<S: Real + Into<f64> + From<f64>, const N: usize>(
    op: Sen3,
    record: &Record,
) -> Output {
    let flat = |v: Vec<Vec3<S>>| {
        v.into_iter()
            .flat_map(|c| c.0)
            .map(Into::into)
            .collect::<Vec<f64>>()
    };
    match op {
        Sen3::Exp => {
            let Some(tau) = record.input("tau").filter(|t| t.len() == 3 + 3 * N) else {
                return Output::new();
            };
            // Through `read_dense`, so the dense order is §1's and not this file's reading of it.
            let tau: Vec<S> = tau.iter().map(|&t| S::from(t)).collect();
            let x = SEn3::<S, N>::exp(&SEn3Tangent::read_dense(&tau));
            let (r, cols) = x.parts();
            let mut out = quat_out(&r.quat());
            out.insert("x".to_string(), flat(cols.to_vec()));
            out
        }
        Sen3::Log => {
            let Some(x) = sen3_of::<S, N>(record, "q", "x") else {
                return Output::new();
            };
            let mut tau = vec![S::zero(); 3 + 3 * N];
            x.log().write_dense(&mut tau);
            Output::from([(
                "tau".to_string(),
                tau.into_iter().map(Into::into).collect::<Vec<f64>>(),
            )])
        }
        Sen3::Ad => {
            let Some(x) = sen3_of::<S, N>(record, "q", "x") else {
                return Output::new();
            };
            sen3_jac_out("Ad", &x.adjoint())
        }
        Sen3::Jr | Sen3::Jl | Sen3::JrInv | Sen3::JlInv => {
            let Some(tau) = record.input("tau").filter(|t| t.len() == 3 + 3 * N) else {
                return Output::new();
            };
            let tau: Vec<S> = tau.iter().map(|&t| S::from(t)).collect();
            let tau = SEn3Tangent::<S, N>::read_dense(&tau);
            // Every arm spelled, no `_`: a variant added to the outer pattern without its own arm
            // here would otherwise be answered with `J_l⁻¹` under its own name, which is the
            // mislabel the `So3` note above describes. A non-exhaustive match is a build error.
            let j = match op {
                Sen3::Jr => SEn3::<S, N>::jr(&tau),
                Sen3::Jl => SEn3::<S, N>::jl(&tau),
                Sen3::JrInv => SEn3::<S, N>::jr_inv(&tau),
                Sen3::JlInv => SEn3::<S, N>::jl_inv(&tau),
                Sen3::Exp | Sen3::Log | Sen3::Ad => return Output::new(),
            };
            sen3_jac_out("J", &j)
        }
    }
}

/// [`sen3_answer`] at the `N` the id names: `N` is a const parameter, so the three widths the
/// corpus holds are three instantiations and the dispatch is this match.
fn sen3_at<S: Real + Into<f64> + From<f64>>(op: Sen3, n: usize, record: &Record) -> Output {
    match n {
        1 => sen3_answer::<S, 1>(op, record),
        2 => sen3_answer::<S, 2>(op, record),
        3 => sen3_answer::<S, 3>(op, record),
        _ => Output::new(),
    }
}

/// A geodesic id this subject answers (`PHASE4.md` §4).
///
/// An enum with an exhaustive `answer`, for the reason the [`So3`] note gives: a `_` arm once
/// scored a committed row under the wrong id's name.
#[derive(Clone, Copy)]
enum Geodesic {
    So3,
    Se3,
}

impl Geodesic {
    const ALL: [(&'static str, Geodesic); 2] = [
        ("so3_geodesic", Geodesic::So3),
        ("se3_geodesic", Geodesic::Se3),
    ];

    fn of_fn(fn_id: &str) -> Option<Self> {
        Self::ALL
            .iter()
            .find(|(name, _)| *name == fn_id)
            .map(|&(_, id)| id)
    }

    /// The shipped answer at `S`; nothing when the record holds no usable input.
    ///
    /// The poses reach the group as the corpus holds them, through `from_quat_unchecked`: the
    /// reference is the geodesic between the two rotations the records *denote*, so it normalizes
    /// both, and a program carrying the stored quaternion reads about half an `u` of `‖q‖² − 1`
    /// for it. That difference is the shipped behaviour and is what the row measures.
    fn answer<S: Real + Into<f64> + From<f64>>(self, record: &Record) -> Output {
        let Some(&[t]) = record.input("t").and_then(|v| v.first_chunk::<1>()) else {
            return Output::new();
        };
        let t = S::from(t);
        match self {
            Geodesic::So3 => {
                let (Some(q0), Some(q1)) = (quat_of::<S>(record, "q0"), quat_of::<S>(record, "q1"))
                else {
                    return Output::new();
                };
                let (x0, x1) = (SO3::from_quat_unchecked(q0), SO3::from_quat_unchecked(q1));
                quat_out(&SO3::geodesic(&x0, &x1, t).quat())
            }
            Geodesic::Se3 => {
                let (Some(x0), Some(x1)) = (
                    sen3_of::<S, 1>(record, "q0", "x0"),
                    sen3_of::<S, 1>(record, "q1", "x1"),
                ) else {
                    return Output::new();
                };
                let (r, cols) = SEn3::<S, 1>::geodesic(&x0, &x1, t).parts();
                let mut out = quat_out(&r.quat());
                out.insert("x".to_string(), cols[0].0.map(Into::into).to_vec());
                out
            }
        }
    }
}

/// An id of `0056` this subject answers: the routines D7 did not reach, at both precisions.
///
/// An enum with an exhaustive `answer`, for the reason the [`So3`] note gives.
#[derive(Clone, Copy)]
enum Linalg {
    SolveCubic,
    Eig3,
    Chol3,
    Chol6,
    CholSolve3,
    CholSolve6,
    Renormalize,
    Sqrt,
    Cbrt,
    Acos,
    SinCos,
    Atan2,
    Div,
}

impl Linalg {
    const ALL: [(&'static str, Linalg); 13] = [
        ("solve_cubic", Linalg::SolveCubic),
        ("eig3", Linalg::Eig3),
        ("chol_n3", Linalg::Chol3),
        ("chol_n6", Linalg::Chol6),
        ("chol_solve_n3", Linalg::CholSolve3),
        ("chol_solve_n6", Linalg::CholSolve6),
        ("quat_renormalize", Linalg::Renormalize),
        ("real_sqrt", Linalg::Sqrt),
        ("real_cbrt", Linalg::Cbrt),
        ("real_acos", Linalg::Acos),
        ("real_sin_cos", Linalg::SinCos),
        ("real_atan2", Linalg::Atan2),
        ("real_div", Linalg::Div),
    ];

    fn of_fn(fn_id: &str) -> Option<Self> {
        Self::ALL
            .iter()
            .find(|(name, _)| *name == fn_id)
            .map(|&(_, id)| id)
    }

    /// The shipped answer at `S`; nothing when the record holds no usable input (an `@f32` record
    /// whose input is not exactly a binary32 included).
    fn answer<S: Real<Mask = bool> + Into<f64>>(self, record: &Record) -> Output {
        let one = |key: &str| scalars::<S>(record, key).and_then(|v| v.first().copied());
        let flat = |v: &[S]| v.iter().map(|&x| x.into()).collect::<Vec<f64>>();
        let bit = |b: bool| if b { 1.0 } else { 0.0 };
        match self {
            Linalg::SolveCubic => {
                let (Some(a), Some(b), Some(c), Some(d)) = (one("a"), one("b"), one("c"), one("d"))
                else {
                    return Output::new();
                };
                let (r, ok) = solve_cubic(a, b, c, d);
                Output::from([
                    ("roots".to_string(), flat(&r.0)),
                    ("valid".to_string(), ok.map(bit).to_vec()),
                ])
            }
            Linalg::Eig3 => {
                let Some(a) = matrix::<S, 3>(record, "A") else {
                    return Output::new();
                };
                let (lambda, v) = eig3(&a);
                Output::from([
                    ("lambda".to_string(), flat(&lambda.0)),
                    ("V".to_string(), cols_out(&v)),
                ])
            }
            Linalg::Chol3 => chol_answer::<S, 3>(record),
            Linalg::Chol6 => chol_answer::<S, 6>(record),
            Linalg::CholSolve3 => chol_solve_answer::<S, 3>(record),
            Linalg::CholSolve6 => chol_solve_answer::<S, 6>(record),
            Linalg::Renormalize => {
                let Some(&[w, x, y, z]) = scalars::<S>(record, "q").as_deref() else {
                    return Output::new();
                };
                // Not through `Quat::from_wxyz_unchecked`: its assert is the very drift the strata hold.
                let mut q = Quat { w, x, y, z };
                q.renormalize();
                quat_out(&q)
            }
            Linalg::Sqrt | Linalg::Cbrt | Linalg::Acos => {
                let Some(x) = one("x") else {
                    return Output::new();
                };
                let x = Dual::<S, 1>::variable(x, 0);
                let r = match self {
                    Linalg::Sqrt => x.sqrt(),
                    Linalg::Cbrt => x.cbrt(),
                    _ => x.acos(),
                };
                Output::from([
                    ("value".to_string(), vec![r.v.into()]),
                    ("d".to_string(), vec![r.d[0].into()]),
                ])
            }
            Linalg::SinCos => {
                let Some(x) = one("x") else {
                    return Output::new();
                };
                let (sin, cos) = Dual::<S, 1>::variable(x, 0).sin_cos();
                Output::from([
                    ("sin".to_string(), vec![sin.v.into()]),
                    ("cos".to_string(), vec![cos.v.into()]),
                    ("d_sin".to_string(), vec![sin.d[0].into()]),
                    ("d_cos".to_string(), vec![cos.d[0].into()]),
                ])
            }
            Linalg::Atan2 | Linalg::Div => {
                let (a, b, names) = match self {
                    Linalg::Atan2 => (one("y"), one("x"), ["d_y", "d_x"]),
                    _ => (one("n"), one("d"), ["d_n", "d_d"]),
                };
                let (Some(a), Some(b)) = (a, b) else {
                    return Output::new();
                };
                let (a, b) = (Dual::<S, 2>::variable(a, 0), Dual::<S, 2>::variable(b, 1));
                let r = match self {
                    Linalg::Atan2 => a.atan2(b),
                    _ => a / b,
                };
                Output::from([
                    ("value".to_string(), vec![r.v.into()]),
                    (names[0].to_string(), vec![r.d[0].into()]),
                    (names[1].to_string(), vec![r.d[1].into()]),
                ])
            }
        }
    }
}

/// The values a record holds under `key` at `S`: binary64 as stored, or at `f32` each exactly a
/// binary32 (`0016`), else nothing.
fn scalars<S: Real>(record: &Record, key: &str) -> Option<Vec<S>> {
    record
        .input(key)?
        .iter()
        .map(|&x| match S::PRECISION {
            Precision::F64 => Some(S::lit(x)),
            Precision::F32 => exact_f32(x).map(|x| S::lit(f64::from(x))),
        })
        .collect()
}

/// The `N x N` matrix a record holds under `key`, column-major (`PHASE1.md` §4.3).
fn matrix<S: Real, const N: usize>(record: &Record, key: &str) -> Option<Matrix<S, N, N>> {
    let m = scalars::<S>(record, key).filter(|m| m.len() == N * N)?;
    Some(Matrix::from_cols(core::array::from_fn(|c| {
        Vector(core::array::from_fn(|r| m[c * N + r]))
    })))
}

/// A square matrix column-major, as the corpus holds one.
fn cols_out<S: Real + Into<f64>, const N: usize>(m: &Matrix<S, N, N>) -> Vec<f64> {
    (0..N).flat_map(|c| m.col(c).0).map(Into::into).collect()
}

fn chol_answer<S: Real<Mask = bool> + Into<f64>, const N: usize>(record: &Record) -> Output {
    let Some(a) = matrix::<S, N>(record, "A") else {
        return Output::new();
    };
    let (l, ok) = chol(&a);
    Output::from([
        ("L".to_string(), cols_out(&l)),
        ("valid".to_string(), vec![if ok { 1.0 } else { 0.0 }]),
    ])
}

/// `chol` then `chol_solve`, the composition a consumer runs (`NUMERICS.md` §14).
fn chol_solve_answer<S: Real<Mask = bool> + Into<f64>, const N: usize>(record: &Record) -> Output {
    let (Some(a), Some(b)) = (
        matrix::<S, N>(record, "A"),
        scalars::<S>(record, "b").filter(|b| b.len() == N),
    ) else {
        return Output::new();
    };
    let x = chol_solve(&chol(&a).0, Vector(core::array::from_fn(|i| b[i])));
    Output::from([("x".to_string(), x.0.map(Into::into).to_vec())])
}

pub(crate) struct Helicoid;

impl Subject for Helicoid {
    fn name(&self) -> &str {
        "helicoid"
    }

    fn supports(&self, fn_id: &str) -> bool {
        Swept::of_fn(fn_id).is_some()
            || So3::of_fn(fn_id).is_some()
            || Sen3::of_fn(fn_id).is_some()
            || Geodesic::of_fn(fn_id).is_some()
            || Linalg::of_fn(fn_id).is_some()
    }

    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        if let Some(id) = Linalg::of_fn(fn_id) {
            return match precision {
                Precision::F64 => id.answer::<f64>(record),
                Precision::F32 => id.answer::<f32>(record),
            };
        }
        if let Some(id) = So3::of_fn(fn_id) {
            // A plain `f32` run skips every `so3_*` id; asking for one by name at `f32` is the
            // harness's error, not this subject's (`So3::answer` says why).
            return match precision {
                Precision::F64 => id.answer::<f64>(record),
                Precision::F32 => Output::new(),
            };
        }
        if let Some(id) = Geodesic::of_fn(fn_id) {
            // `f64` as for the other vector ids: the geodesic strata have no `@f32` twin.
            return match precision {
                Precision::F64 => id.answer::<f64>(record),
                Precision::F32 => Output::new(),
            };
        }
        if let Some((op, n)) = Sen3::of_fn(fn_id) {
            // `f32` as for the `so3_*` ids: no vector id has an `@f32` stratum until a record
            // extends `0016`.
            return match precision {
                Precision::F64 => sen3_at::<f64>(op, n, record),
                Precision::F32 => Output::new(),
            };
        }
        let Some(id) = Swept::of_fn(fn_id) else {
            return Output::new();
        };
        match precision {
            Precision::F64 => answer::<f64>(id, record),
            Precision::F32 => answer::<f32>(id, record),
        }
    }
}

/// The library's own host-`std` twin: this subject's program at [`Host`], whose transcendentals
/// are Rust `std`'s — the host's libm, so glibc on Linux — where a default build routes them
/// through the `libm` crate (D16).
///
/// `judge_with` reads a twin to attribute a domination failure to D16, and the twin it read could
/// only ever be `seeded:host-std`: the *stand-in's* program, which answers `exp`, `jr` and `jl`
/// and nothing else of `PHASE1.md` §10's. That left every `sen3_log`, `sen3_jr_inv` and
/// `sen3_jl_inv` failure unattributed, and attributed the rest through a program that is not the
/// one being judged. This twin is the candidate's own, for every id it answers, which is what
/// `judge_with`'s own words ask for — "the candidate's own program with one variable changed".
///
/// It is planted, so a plain `just conformance` skips it and no bar reads its rows: the envelope
/// comparing `helicoid` with a differently-rounded copy of itself would measure nothing.
pub(crate) struct HostStd;

impl Subject for HostStd {
    fn name(&self) -> &str {
        "helicoid:host-std"
    }

    /// Every id the candidate answers, so D16's price is measured wherever the candidate is scored.
    fn supports(&self, fn_id: &str) -> bool {
        Helicoid.supports(fn_id)
    }

    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        // `f64` alone, as `seeded:host-std` is: `std`'s `f32` transcendentals are not the host's
        // `f32` libm on every target (LLVM may widen), so an `f32` twin would change two things.
        if precision != Precision::F64 {
            return Output::new();
        }
        if let Some(id) = Linalg::of_fn(fn_id) {
            return id.answer::<Host>(record);
        }
        if let Some(id) = So3::of_fn(fn_id) {
            return id.answer::<Host>(record);
        }
        if let Some(id) = Geodesic::of_fn(fn_id) {
            return id.answer::<Host>(record);
        }
        if let Some((op, n)) = Sen3::of_fn(fn_id) {
            return sen3_at::<Host>(op, n, record);
        }
        let Some(id) = Swept::of_fn(fn_id) else {
            return Output::new();
        };
        answer::<Host>(id, record)
    }
}

pub(crate) fn host_std() -> Registered {
    let version = env!("CARGO_PKG_VERSION");
    Registered {
        version: format!("{version}@host-std"),
        version_f32: format!("{version}@host-std"),
        no_f32: Some(
            "the twin measures the host's binary64 transcendentals; `0016`'s `@f32` strata are \
             the coefficient ids' own"
                .to_string(),
        ),
        planted: true,
        subject: Box::new(HostStd),
    }
}

pub(crate) fn registered() -> Registered {
    let version = env!("CARGO_PKG_VERSION").to_string();
    Registered {
        version: version.clone(),
        version_f32: version,
        no_f32: None,
        subject: Box::new(Helicoid),
        planted: false,
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;
    use crate::conformance::corpus::{self, Record};
    use crate::conformance::corpus_dir;
    use crate::conformance::metric::{rule, Score};
    use crate::conformance::subject::registry;
    use crate::conformance::testkit::record;
    use crate::seeded::{d12, Seeded, Series, D1};

    type F32 = Dual<f32, 1>;

    /// Every record of `coeff_<id>` at `precision`: the strata without the suffix at binary64, the
    /// `@f32` ones at binary32.
    fn records(id: Swept, precision: Precision) -> Result<Vec<Record>, String> {
        let dir = corpus_dir()?;
        let fn_id = format!("coeff_{}", id.name());
        let entry = corpus::manifest(&dir)?
            .into_iter()
            .find(|e| e.fn_id == fn_id);
        let all = corpus::read(&dir, &entry.ok_or("no corpus file")?)?;
        let f32 = precision == Precision::F32;
        Ok(all
            .into_iter()
            .filter(|r| r.is_f32_stratum() == f32)
            .collect())
    }

    /// The largest score per stratum of `subject` on `id`, `None` for a non-finite record.
    fn maxima(
        subject: &dyn Subject,
        id: Swept,
        precision: Precision,
    ) -> Result<BTreeMap<String, Option<f64>>, String> {
        let fn_id = format!("coeff_{}", id.name());
        let rule = rule(&fn_id).ok_or("no rule")?;
        let mut out: BTreeMap<String, Option<f64>> = BTreeMap::new();
        for rec in records(id, precision)? {
            let answer = subject.eval(&fn_id, &rec, precision);
            let score = match rule.score(&rec, &answer, precision)? {
                Score::Finite(u) => Some(u),
                _ => None,
            };
            let slot = out.entry(rec.stratum.clone()).or_insert(Some(0.0));
            *slot = slot.zip(score).map(|(a, b)| a.max(b));
        }
        Ok(out)
    }

    /// Equal to the bit, but for the payload of a NaN, which Rust leaves unspecified.
    fn same(a: &Output, b: &Output) -> bool {
        let bits = |o: &Output| -> Vec<(String, Vec<u64>)> {
            let nan = f64::NAN.to_bits();
            let map = |x: &f64| if x.is_nan() { nan } else { x.to_bits() };
            o.iter()
                .map(|(k, v)| (k.clone(), v.iter().map(map).collect()))
                .collect()
        };
        bits(a) == bits(b)
    }

    #[test]
    fn the_registry_holds_helicoid_as_a_plain_subject_with_both_precisions() {
        let all = registry();
        let h = &all[0];
        assert_eq!((h.subject.name(), h.planted), ("helicoid", false));
        assert_eq!(h.version, env!("CARGO_PKG_VERSION"));
        assert_eq!(h.version_at(Precision::F32), h.version_at(Precision::F64));
        assert!(h.no_f32.is_none());
        assert_eq!(
            all.iter()
                .filter(|r| r.subject.name() == "helicoid")
                .count(),
            1
        );
    }

    /// The library's twin answers the candidate's ids, at `f64` only, and its answers **differ** —
    /// so the swap reaches the group code and is a change, not a relabelling.
    ///
    /// The seeded twin's `the_scalar_differs_from_the_libm_crate_in_the_transcendentals_only`
    /// proves that of the scalar; this proves it of the path the groups take to it.
    #[test]
    fn the_library_twin_answers_every_candidate_id_and_differs_somewhere() -> Result<(), String> {
        let twin = HostStd;
        for (id, _) in So3::ALL {
            assert!(twin.supports(id), "{id}");
        }
        for id in Swept::ALL.map(|i| format!("coeff_{}", i.name())) {
            assert!(twin.supports(&id), "{id}");
        }
        for op in ["exp", "log", "ad", "jr", "jl", "jr_inv", "jl_inv"] {
            let id = format!("sen3_{op}_n2");
            assert!(twin.supports(&id), "{id}");
        }
        assert!(!twin.supports("so3_jl_jr"));
        // Sampled, not one point: two libms agreeing to the bit at a given argument is ordinary,
        // so a single quaternion would make this test a property of the host's libc. The claim is
        // that *some* `θ` differs, which is what `host::the_scalar_differs_from_the_libm_crate_
        // in_the_transcendentals_only` asserts of the scalar; this asserts it of the group path.
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let mut next = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let bits = |o: &Output| o["tau"].iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        let mut differed = 0usize;
        for _ in 0..2000 {
            // A unit quaternion at a random `θ` in `(0, π)` about a random axis, so `Log` is on
            // its exact arm and reaches `atan2`.
            let theta = core::f64::consts::PI * next();
            let v = [next() - 0.5, next() - 0.5, next() - 0.5];
            let n = v
                .iter()
                .map(|c| c * c)
                .sum::<f64>()
                .sqrt()
                .max(f64::MIN_POSITIVE);
            let (sin, cos) = (theta * 0.5).sin_cos();
            let q = [cos, sin * v[0] / n, sin * v[1] / n, sin * v[2] / n];
            let rec = record(&[("q", &q), ("x", &[0.5, -2.0, 0.25])], &[])?;
            let (mine, theirs) = (
                Helicoid.eval("sen3_log_n1", &rec, Precision::F64),
                twin.eval("sen3_log_n1", &rec, Precision::F64),
            );
            assert_eq!(mine.keys().collect::<Vec<_>>(), ["tau"]);
            assert_eq!(theirs.keys().collect::<Vec<_>>(), ["tau"]);
            differed += usize::from(bits(&mine) != bits(&theirs));
        }
        assert!(differed > 0, "the swap changed nothing over 2000 draws");
        let rec = record(
            &[("q", &[1.0, 0.0, 0.0, 0.0]), ("x", &[0.5, -2.0, 0.25])],
            &[],
        )?;
        // `f64` only, for the reason the subject's doc gives.
        assert!(twin.eval("sen3_log_n1", &rec, Precision::F32).is_empty());
        assert!(twin.eval("so3_exp", &rec, Precision::F32).is_empty());
        Ok(())
    }

    #[test]
    fn helicoid_supports_the_coefficient_and_so3_ids_and_answers_nothing_it_cannot(
    ) -> Result<(), String> {
        for id in Swept::ALL {
            assert!(Helicoid.supports(&format!("coeff_{}", id.name())), "{id:?}");
        }
        // Every `so3_*` id the corpus holds, and every `sen3_*` one (`PHASE3.md` §5).
        for (id, _) in So3::ALL {
            assert!(Helicoid.supports(id), "{id}");
        }
        for op in ["exp", "log", "ad", "jr", "jl", "jr_inv", "jl_inv"] {
            for n in 1..=3 {
                let id = format!("sen3_{op}_n{n}");
                assert!(Helicoid.supports(&id), "{id}");
            }
        }
        for id in [
            "coeff_",
            "coeff_f",
            "coeff_series",
            "so3_jl_jr",
            // A width the corpus does not hold, and an op that is not one.
            "sen3_exp_n4",
            "sen3_exp_n0",
            "sen3_adjoint_n1",
            "sen3_exp",
        ] {
            assert!(!Helicoid.supports(id), "{id}");
            let rec = record(&[("theta", &[0.5]), ("n", &[0.5]), ("w", &[1.0])], &[])?;
            assert!(Helicoid.eval(id, &rec, Precision::F64).is_empty(), "{id}");
        }
        // An `so3_*` id reads `phi`, `q`, `p` or `R`; a record without them is answered with
        // nothing, and every one of them is answered with nothing at `f32`, where no `@f32`
        // stratum exists for a vector id.
        let empty = record(&[("theta", &[0.5])], &[])?;
        for (id, _) in So3::ALL {
            assert!(Helicoid.eval(id, &empty, Precision::F64).is_empty(), "{id}");
            assert!(Helicoid.eval(id, &empty, Precision::F32).is_empty(), "{id}");
        }
        // Every `sen3_*` id: the fields the corpus names, the dense width, nothing from a record
        // without its input or at `f32`, and -- for `Ad` at the identity rotation, where the
        // adjoint is `I + ε[x]_×` -- the column-major order and the block the corpus puts it in.
        for n in 1..=3usize {
            let d = 3 + 3 * n;
            let tau: Vec<f64> = (0..d).map(|i| 0.1 * (i as f64 + 1.0)).collect();
            let x: Vec<f64> = (0..3 * n).map(|i| 0.25 * (i as f64 + 1.0)).collect();
            let rec = record(
                &[("tau", &tau), ("q", &[1.0, 0.0, 0.0, 0.0]), ("x", &x)],
                &[],
            )?;
            for (op, field, len) in [
                ("exp", "q", 4),
                ("log", "tau", d),
                ("ad", "Ad", d * d),
                ("jr", "J", d * d),
                ("jl", "J", d * d),
                ("jr_inv", "J", d * d),
                ("jl_inv", "J", d * d),
            ] {
                let id = format!("sen3_{op}_n{n}");
                let out = Helicoid.eval(&id, &rec, Precision::F64);
                let keys: Vec<&str> = out.keys().map(String::as_str).collect();
                let want: Vec<&str> = match op {
                    "exp" => vec!["q", "x"],
                    _ => vec![field],
                };
                assert_eq!(keys, want, "{id}");
                assert_eq!(out[field].len(), len, "{id}");
                assert!(Helicoid.eval(&id, &rec, Precision::F32).is_empty(), "{id}");
                assert!(
                    Helicoid.eval(&id, &empty, Precision::F64).is_empty(),
                    "{id}"
                );
            }
            // `[x₁]_×` sits at block row 1, block column 0: entries `(4, 0)` and `(5, 0)` of the
            // dense matrix are `x_z` and `−x_y`, at flat indices `4` and `5` column-major.
            let ad = Helicoid.eval(&format!("sen3_ad_n{n}"), &rec, Precision::F64);
            assert_eq!(ad["Ad"][4].to_bits(), x[2].to_bits(), "n{n}");
            assert_eq!(ad["Ad"][5].to_bits(), (-x[1]).to_bits(), "n{n}");
            assert_eq!(ad["Ad"][0].to_bits(), 1.0_f64.to_bits(), "n{n}");
        }
        // `r` reads `n` and `w`, the others `theta`; a record without them is answered with nothing.
        let theta = record(&[("theta", &[0.5])], &[])?;
        let nw = record(&[("n", &[0.5]), ("w", &[1.0])], &[])?;
        for precision in [Precision::F64, Precision::F32] {
            assert!(Helicoid.eval("coeff_r", &theta, precision).is_empty());
            assert!(Helicoid.eval("coeff_k", &nw, precision).is_empty());
            let out = Helicoid.eval("coeff_r", &nw, precision);
            assert_eq!(out.keys().collect::<Vec<_>>(), ["d_branch", "value"]);
            let out = Helicoid.eval("coeff_e", &theta, precision);
            assert_eq!(out.keys().collect::<Vec<_>>(), ["d_branch", "value"]);
        }
        // An `f32` subject receives a binary32 or nothing: 0.1 is no binary32, and is not rounded.
        let tenth = record(&[("theta", &[0.1])], &[])?;
        assert!(Helicoid.eval("coeff_k", &tenth, Precision::F32).is_empty());
        assert_eq!(Helicoid.eval("coeff_k", &tenth, Precision::F64).len(), 2);
        Ok(())
    }

    #[test]
    fn helicoid_is_the_seeded_correct_kernel_bit_for_bit_at_both_precisions() -> Result<(), String>
    {
        // The two kernels share every candidate (`the_shipped_arms_are_the_seeded_arms_bit_for_bit`
        // and the two committed sweeps), so they agree on every record, in value and derivative.
        let correct = Seeded::generated();
        let mut n = 0;
        for precision in [Precision::F64, Precision::F32] {
            for id in Swept::ALL {
                let fn_id = format!("coeff_{}", id.name());
                for rec in records(id, precision)? {
                    let (h, c) = (
                        Helicoid.eval(&fn_id, &rec, precision),
                        correct.eval(&fn_id, &rec, precision),
                    );
                    assert!(
                        !h.is_empty() && same(&h, &c),
                        "{fn_id} {precision:?} {}",
                        rec.stratum
                    );
                    n += 1;
                }
            }
        }
        // Seven ids of 3420 records and `coeff_r` of 3426, over both precisions.
        assert_eq!(n, 7 * 3420 + 3426);
        Ok(())
    }

    #[test]
    fn no_record_of_any_stratum_scores_non_finite_at_either_precision() -> Result<(), String> {
        for precision in [Precision::F64, Precision::F32] {
            for id in Swept::ALL {
                for (stratum, max) in maxima(&Helicoid, id, precision)? {
                    let max = max.ok_or(format!("{id:?} {precision:?} {stratum}: non-finite"))?;
                    assert!(max.is_finite() && max >= 0.0, "{id:?} {stratum}: {max}");
                }
            }
        }
        Ok(())
    }

    /// The strata on which the prior (four terms below `z = 0.01`) has a smaller maximum than
    /// `helicoid`, as `(id, precision, stratum)`. It is `tf_tree`'s D12 for `a`, `b`, `c`; for `k`,
    /// `d`, `e` and `cos θ/2` D12 applied to a coefficient it was not defined for, and for `r`
    /// `(4 terms, s < 0.01)`, `θ = 0.2` and not D12's (0014 (draft) questions 10 and 29).
    fn prior_wins() -> Result<Vec<(&'static str, Precision, String)>, String> {
        let dir = corpus_dir()?;
        let (s64, s32) = (Series::<D1>::load(&dir)?, Series::<F32>::load(&dir)?);
        let prior64 = Seeded::uniform(&s64, d12(&s64)?);
        let prior32 = Seeded::uniform_f32(&s32, d12(&s32)?);
        let mut wins = Vec::new();
        for precision in [Precision::F64, Precision::F32] {
            for id in Swept::ALL {
                let prior = match precision {
                    Precision::F64 => &prior64,
                    Precision::F32 => &prior32,
                };
                let (h, p) = (
                    maxima(&Helicoid, id, precision)?,
                    maxima(prior, id, precision)?,
                );
                for (stratum, h) in h {
                    if p[&stratum].zip(h).is_some_and(|(p, h)| p < h) {
                        wins.push((id.name(), precision, stratum));
                    }
                }
            }
        }
        Ok(wins)
    }

    /// **Two strata, both at `θ = π − 0.1`, and no other.** `0039` lifted the sweep's grid to span
    /// the domain and the term cap to the corpus's series length, and the nine strata where the
    /// `tf_tree` D12 prior used to beat `helicoid` — all of them `cos θ/2`, where the old switch
    /// sat at `θ < 7.5e-8` because every candidate under `theta:1e0`'s objective tied — are gone:
    /// `cos θ/2`'s switch is now `z = 5.23` at `f64`, so the series arm serves the small-angle
    /// strata the prior used to win.
    ///
    /// What is left is near `π`, where the prior's four-term series below `θ = 0.1` is not what
    /// decides: both candidates are on their exact arms there and the margin is the assembly's.
    /// A stratum's maximum is what the bars compare (D8), the sweep's objective one maximum over
    /// all of them: open, 0014 (draft) question 30.
    #[test]
    fn the_prior_beats_helicoid_near_pi_on_two_strata_and_on_no_other() -> Result<(), String> {
        let want = [
            ("b", Precision::F64, "theta:pi-1e-1"),
            ("d", Precision::F32, "theta:pi-1e-1@f32"),
        ];
        let mut got = prior_wins()?;
        got.sort_by(|a, b| (a.0, a.2.clone()).cmp(&(b.0, b.2.clone())));
        let mut want = want.map(|(c, p, s)| (c, p, s.to_string())).to_vec();
        want.sort_by(|a, b| (a.0, a.2.clone()).cmp(&(b.0, b.2.clone())));
        assert_eq!(got, want);
        Ok(())
    }
}
