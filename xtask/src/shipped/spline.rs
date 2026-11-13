//! The spline ids of `0116`, answered by `helicoid-spline` (`0114` step 4).
//!
//! A record's window length is its order: one monomorphization per `K` in `2..=6` and derivative
//! order `D`, so each id runs the code a consumer of that order runs. Every Jacobian is the right
//! one, through `JacSegment::kinematics_jacobians` at the id's own `D`. A derivative output is
//! scaled to unit size as its id states, `ω / h⁻¹` and `ω̇ / (h⁻¹ h⁻¹)`. The record readers here
//! are shared with `seeded:spline`, so the two subjects read the same window.

use helicoid::{Jac, LieGroup, Quat, Right, Rn, SEn3, Tangent, SE3, SO3};
use helicoid_linalg::{Real, StridedMut, Vector};
use helicoid_spline::{derivative_gram, JacSegment, Order, Segment, ValidOrder, Weights};

use super::{lookup, scalars};
use crate::conformance::corpus::Record;
use crate::conformance::subject::Output;

/// The group a spline id evaluates on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum On {
    So3,
    Rn,
    Se3,
}

/// A spline id: a kinematic output (or its control Jacobians) at derivative order `order`, or the
/// Gram matrix.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SplineId {
    Eval { on: On, jac: bool, order: usize },
    Gram,
}

const fn eval(on: On, jac: bool, order: usize) -> SplineId {
    SplineId::Eval { on, jac, order }
}

impl SplineId {
    const ALL: [(&'static str, SplineId); 16] = [
        ("so3_spline_value", eval(On::So3, false, 0)),
        ("so3_spline_vel", eval(On::So3, false, 1)),
        ("so3_spline_acc", eval(On::So3, false, 2)),
        ("so3_spline_value_jac", eval(On::So3, true, 0)),
        ("so3_spline_vel_jac", eval(On::So3, true, 1)),
        ("so3_spline_acc_jac", eval(On::So3, true, 2)),
        ("rn_spline_value", eval(On::Rn, false, 0)),
        ("rn_spline_vel", eval(On::Rn, false, 1)),
        ("rn_spline_acc", eval(On::Rn, false, 2)),
        ("rn_spline_value_jac", eval(On::Rn, true, 0)),
        ("rn_spline_vel_jac", eval(On::Rn, true, 1)),
        ("rn_spline_acc_jac", eval(On::Rn, true, 2)),
        ("se3_spline_value", eval(On::Se3, false, 0)),
        ("se3_spline_vel", eval(On::Se3, false, 1)),
        ("se3_spline_acc", eval(On::Se3, false, 2)),
        ("rn_spline_gram", SplineId::Gram),
    ];

    pub(crate) fn of_fn(fn_id: &str) -> Option<Self> {
        lookup(&Self::ALL, fn_id)
    }

    /// The answer at `S`, nothing for a record it cannot read.
    pub(crate) fn answer<S: Real + Into<f64>>(self, record: &Record) -> Output {
        let out = match self {
            SplineId::Gram => gram::<S>(record),
            SplineId::Eval { on, jac, order } => {
                segment::<S>(record).and_then(|(u, inv_h)| match on {
                    On::So3 => by_len(&rotations::<S>(record)?, u, inv_h, jac, order),
                    On::Rn => by_len(&points::<S>(record)?, u, inv_h, jac, order),
                    On::Se3 => by_len(&poses::<S>(record)?, u, inv_h, jac, order),
                })
            }
        };
        out.map_or_else(Output::new, |out| {
            out.into_iter().map(|(k, v)| (k.to_string(), v)).collect()
        })
    }
}

/// The fields a group's outputs are written under, and its Jacobian's entries.
trait Answered<S: Real + Into<f64>>: LieGroup<S> {
    /// `ω`'s and `ω̇`'s fields.
    const TANGENT: [&'static str; 2];
    /// The Jacobian's field.
    const JAC: &'static str;
    fn value_out(&self) -> Vec<(&'static str, Vec<f64>)>;
    /// The entries of one control's Jacobian, as the id's schema holds them.
    fn jac_entries(j: &Self::Jac) -> Vec<S>;
}

/// `j`, column-major.
fn dense_jac<S: Real, G: LieGroup<S>>(j: &G::Jac) -> Vec<S> {
    let n = G::DOF;
    let mut buf = vec![S::zero(); n * n];
    Jac::<S, G::Tangent>::write_dense(j, &mut StridedMut::col_major(&mut buf, n, n));
    buf
}

impl<S: Real + Into<f64>> Answered<S> for SO3<S> {
    const TANGENT: [&'static str; 2] = ["omega", "alpha"];
    const JAC: &'static str = "jac";
    fn value_out(&self) -> Vec<(&'static str, Vec<f64>)> {
        vec![("q", quat_of(self))]
    }
    fn jac_entries(j: &Self::Jac) -> Vec<S> {
        dense_jac::<S, Self>(j)
    }
}

impl<S: Real + Into<f64>> Answered<S> for Rn<S, 3> {
    const TANGENT: [&'static str; 2] = ["v", "a"];
    const JAC: &'static str = "k";
    fn value_out(&self) -> Vec<(&'static str, Vec<f64>)> {
        vec![("p", self.0 .0.map(Into::into).to_vec())]
    }
    /// `k` of `k I`.
    fn jac_entries(j: &Self::Jac) -> Vec<S> {
        vec![dense_jac::<S, Self>(j)[0]]
    }
}

impl<S: Real + Into<f64>> Answered<S> for SE3<S> {
    const TANGENT: [&'static str; 2] = ["twist", "twist"];
    const JAC: &'static str = "jac";
    fn value_out(&self) -> Vec<(&'static str, Vec<f64>)> {
        let (r, t) = self.parts();
        vec![("q", quat_of(&r)), ("x", t[0].0.map(Into::into).to_vec())]
    }
    fn jac_entries(j: &Self::Jac) -> Vec<S> {
        dense_jac::<S, Self>(j)
    }
}

fn by_len<S: Real + Into<f64>, G: Answered<S>>(
    x: &[G],
    u: S,
    inv_h: S,
    jac: bool,
    order: usize,
) -> Option<Vec<(&'static str, Vec<f64>)>> {
    match x.len() {
        2 => by_order::<S, G, 2>(x.try_into().ok()?, u, inv_h, jac, order),
        3 => by_order::<S, G, 3>(x.try_into().ok()?, u, inv_h, jac, order),
        4 => by_order::<S, G, 4>(x.try_into().ok()?, u, inv_h, jac, order),
        5 => by_order::<S, G, 5>(x.try_into().ok()?, u, inv_h, jac, order),
        6 => by_order::<S, G, 6>(x.try_into().ok()?, u, inv_h, jac, order),
        _ => None,
    }
}

fn by_order<S: Real + Into<f64>, G: Answered<S>, const K: usize>(
    x: &[G; K],
    u: S,
    inv_h: S,
    jac: bool,
    order: usize,
) -> Option<Vec<(&'static str, Vec<f64>)>>
where
    Order<K>: ValidOrder,
{
    Some(match order {
        0 => at::<S, G, K, 0>(x, u, inv_h, jac),
        1 => at::<S, G, K, 1>(x, u, inv_h, jac),
        2 => at::<S, G, K, 2>(x, u, inv_h, jac),
        _ => return None,
    })
}

/// The id's output at derivative order `D`: the order's kinematic output, or its row of
/// Jacobians, one control after another.
fn at<S: Real + Into<f64>, G: Answered<S>, const K: usize, const D: usize>(
    x: &[G; K],
    u: S,
    inv_h: S,
    jac: bool,
) -> Vec<(&'static str, Vec<f64>)>
where
    Order<K>: ValidOrder,
{
    let w = Weights::<S, K, D>::at(u, inv_h);
    if jac {
        let (_, j) = JacSegment::new(x).kinematics_jacobians::<Right, D>(&w);
        let flat = [&j.value, &j.vel, &j.acc][D]
            .iter()
            .flat_map(G::jac_entries)
            .collect();
        return vec![(G::JAC, scaled(flat, D, inv_h))];
    }
    let k = Segment::new(x).kinematics(&w);
    match D {
        0 => k.value.value_out(),
        1 => vec![(G::TANGENT[0], scaled(dense(&k.vel, G::DOF), 1, inv_h))],
        _ => vec![(G::TANGENT[1], scaled(dense(&k.acc, G::DOF), 2, inv_h))],
    }
}

/// `G^(r)` at the record's `k`, `r` and `inv_h`, column-major.
fn gram<S: Real + Into<f64>>(record: &Record) -> Option<Vec<(&'static str, Vec<f64>)>> {
    let k = scalars::<S>(record, "k")?.first()?.value_f64();
    let r = scalars::<S>(record, "r")?.first()?.value_f64();
    let inv_h = *scalars::<S>(record, "inv_h")?.first()?;
    fn cols<S: Real + Into<f64>, const K: usize>(m: helicoid_linalg::Matrix<S, K, K>) -> Vec<f64> {
        (0..K)
            .flat_map(|t| (0..K).map(move |s| m.get(s, t).into()))
            .collect()
    }
    let g = match (k as usize, r as usize) {
        (2, 1) => cols(derivative_gram::<S, 2, 1>(inv_h)),
        (2, 2) => cols(derivative_gram::<S, 2, 2>(inv_h)),
        (3, 1) => cols(derivative_gram::<S, 3, 1>(inv_h)),
        (3, 2) => cols(derivative_gram::<S, 3, 2>(inv_h)),
        (4, 1) => cols(derivative_gram::<S, 4, 1>(inv_h)),
        (4, 2) => cols(derivative_gram::<S, 4, 2>(inv_h)),
        (5, 1) => cols(derivative_gram::<S, 5, 1>(inv_h)),
        (5, 2) => cols(derivative_gram::<S, 5, 2>(inv_h)),
        (6, 1) => cols(derivative_gram::<S, 6, 1>(inv_h)),
        (6, 2) => cols(derivative_gram::<S, 6, 2>(inv_h)),
        _ => return None,
    };
    Some(vec![("g", g)])
}

/// `u` and `inv_h` of a record.
pub(crate) fn segment<S: Real>(record: &Record) -> Option<(S, S)> {
    let u = *scalars::<S>(record, "u")?.first()?;
    let inv_h = *scalars::<S>(record, "inv_h")?.first()?;
    Some((u, inv_h))
}

/// The window's rotations: `q` holds one quaternion per control, `w` first.
pub(crate) fn rotations<S: Real>(record: &Record) -> Option<Vec<SO3<S>>> {
    let q = scalars::<S>(record, "q")?;
    let rotation = |c: &[S]| {
        SO3::from_quat_unchecked(Quat {
            w: c[0],
            x: c[1],
            y: c[2],
            z: c[3],
        })
    };
    (q.len() % 4 == 0).then(|| q.chunks(4).map(rotation).collect())
}

/// The window's poses: `rotations` with `x`'s translations.
pub(crate) fn poses<S: Real>(record: &Record) -> Option<Vec<SE3<S>>> {
    let r = rotations::<S>(record)?;
    let x = scalars::<S>(record, "x").filter(|x| x.len() == 3 * r.len())?;
    let pose = |(r, t): (&SO3<S>, &[S])| SEn3::from_parts(*r, [Vector([t[0], t[1], t[2]])]);
    Some(r.iter().zip(x.chunks(3)).map(pose).collect())
}

/// The window's points.
pub(crate) fn points<S: Real>(record: &Record) -> Option<Vec<Rn<S, 3>>> {
    let p = scalars::<S>(record, "p")?;
    (p.len() % 3 == 0).then(|| {
        p.chunks(3)
            .map(|c| Rn(Vector([c[0], c[1], c[2]])))
            .collect()
    })
}

pub(crate) fn dense<S: Real, T: Tangent<S>>(t: &T, n: usize) -> Vec<S> {
    let mut out = vec![S::zero(); n];
    t.write_dense(&mut out);
    out
}

/// `values`, each divided by `inv_h` once per derivative order, then widened.
pub(crate) fn scaled<S: Real + Into<f64>>(values: Vec<S>, order: usize, inv_h: S) -> Vec<f64> {
    let by = match order {
        0 => S::one(),
        1 => inv_h,
        _ => inv_h * inv_h,
    };
    values
        .into_iter()
        .map(|v| if order == 0 { v } else { v / by })
        .map(Into::into)
        .collect()
}

pub(crate) fn quat_of<S: Real + Into<f64>>(r: &SO3<S>) -> Vec<f64> {
    let q = r.quat();
    [q.w, q.x, q.y, q.z].map(Into::into).to_vec()
}

#[cfg(test)]
mod tests {
    use helicoid_linalg::Precision;

    use super::SplineId;
    use crate::conformance::corpus;
    use crate::conformance::corpus_dir;
    use crate::conformance::subject::{Output, Subject};
    use crate::seeded::spline::{Spline, IDS};
    use crate::shipped::Helicoid;

    /// Every corpus record of every spline id, at both precisions: the crate answers with the
    /// seeded evaluator's fields, finite, and within a few ulps of the output's size of it. The
    /// corpus run (`just conformance`) is the bar; this is its cheap regression.
    #[test]
    fn the_crate_answers_every_spline_record_as_section_21_does() -> Result<(), String> {
        let dir = corpus_dir()?;
        let manifest = corpus::manifest(&dir)?;
        let mut n = 0;
        for id in IDS {
            assert!(
                SplineId::of_fn(id).is_some() && Helicoid.supports(id),
                "{id}"
            );
            let entry = manifest.iter().find(|e| e.fn_id == id).ok_or(id)?;
            for record in corpus::read(&dir, entry)? {
                n += 1;
                let (precision, eps) = if record.is_f32_stratum() {
                    (Precision::F32, f64::from(f32::EPSILON))
                } else {
                    (Precision::F64, f64::EPSILON)
                };
                let (a, b): (Output, Output) = (
                    Helicoid.eval(id, &record, precision),
                    Spline(None).eval(id, &record, precision),
                );
                let fields = |o: &Output| {
                    o.iter()
                        .map(|(k, v)| (k.clone(), v.len()))
                        .collect::<Vec<_>>()
                };
                assert!(
                    !a.is_empty() && fields(&a) == fields(&b),
                    "{id} {}",
                    record.stratum
                );
                for (k, va) in &a {
                    let vb = &b[k];
                    let size = vb.iter().fold(1.0f64, |m, x| m.max(x.abs()));
                    for (x, y) in va.iter().zip(vb) {
                        assert!(x.is_finite(), "{id} {} {k}", record.stratum);
                        assert!(
                            (x - y).abs() <= 64.0 * eps * size,
                            "{id} {} {k}: {x:e} against {y:e}",
                            record.stratum
                        );
                    }
                }
            }
        }
        // The committed tier: 855 binary64 records and their 855 `@f32` twins.
        assert_eq!(n, 1710);
        Ok(())
    }
}
