//! `helicoid:composed-inv-act`: `0103`'s two ids answered the way a consumer composed them before
//! the record, `X.inverse() * p` and `act_jacobians(X⁻¹, p) · inverse_jacobian(X)`: the pair's
//! `NUMERICS.md` §14 twin, `reference::{so3, se3}_inv_act_jacobian`, kept as a subject so
//! `cargo xtask compare` reruns the record's *Context* over the corpus. The composition loses
//! `u‖t‖/‖x‖` on `tx:*` (the point and `Jr`) and `u‖t‖/‖p‖` on `tp:*` (`Jl`), the shipped pair
//! neither.
//!
//! Planted, as `helicoid:homogeneous-r` is: a twin is not a candidate, and no bar reads its rows.

use helicoid::{reference, Left, LieGroup, Right, SE3};
use helicoid_linalg::{Matrix, Point, Precision, Real, StridedMut, Vector};

use super::{inv_act_inputs, inv_act_out, Helicoid, InvActId};
use crate::conformance::corpus::Record;
use crate::conformance::subject::{Output, Registered, Subject};

/// The `3 x C` matrix a twin writes, read back.
fn written<S: Real, const C: usize>(write: impl FnOnce(&mut StridedMut<'_, S>)) -> Matrix<S, 3, C> {
    let mut buf = [[S::zero(); 3]; C];
    write(&mut StridedMut::col_major(buf.as_flattened_mut(), 3, C));
    Matrix::from_cols(buf.map(Vector))
}

fn answer<S: Real + Into<f64>>(id: InvActId, record: &Record) -> Output {
    let Some((r, t, p)) = inv_act_inputs::<S>(id, record) else {
        return Output::new();
    };
    match id {
        InvActId::So3 => {
            let ri = r.inverse();
            let jr = written::<S, 3>(|o| reference::so3_inv_act_jacobian::<S, Right>(&r, p, o));
            let jl = written::<S, 3>(|o| reference::so3_inv_act_jacobian::<S, Left>(&r, p, o));
            inv_act_out(id, ri.act(p), &jr, &jl, &ri.act_jacobians::<Right>(p).1)
        }
        InvActId::Se3 => {
            let (x, pt) = (SE3::from_rt(r, t), Point(p.0));
            let xi = x.inverse();
            let jr = written::<S, 6>(|o| reference::se3_inv_act_jacobian::<S, Right>(&x, pt, o));
            let jl = written::<S, 6>(|o| reference::se3_inv_act_jacobian::<S, Left>(&x, pt, o));
            let dp = xi.act_jacobians::<Right>(pt).1;
            inv_act_out(id, Vector((xi * pt).0), &jr, &jl, &dp)
        }
    }
}

pub(crate) struct ComposedInvAct;

impl Subject for ComposedInvAct {
    fn name(&self) -> &str {
        "helicoid:composed-inv-act"
    }

    fn supports(&self, fn_id: &str) -> bool {
        InvActId::of_fn(fn_id).is_some() && Helicoid.supports(fn_id)
    }

    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        let Some(id) = InvActId::of_fn(fn_id) else {
            return Output::new();
        };
        match precision {
            Precision::F64 => answer::<f64>(id, record),
            Precision::F32 => answer::<f32>(id, record),
        }
    }
}

pub(crate) fn registered() -> Registered {
    let version = format!("{}@composed-inv-act", env!("CARGO_PKG_VERSION"));
    Registered {
        version: version.clone(),
        version_f32: version,
        no_f32: None,
        subject: Box::new(ComposedInvAct),
        planted: true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::testkit::record;

    /// At the identity rotation and integer inputs every product of the composition is exact,
    /// so it is the shipped answer on every field, translation terms included: the right block
    /// cancels `[p]_× − [t]_×` to `[p − t]_×` and the left `[x]_× + [t]_×` to `[p]_×`, exactly.
    /// Both precisions, both ids.
    #[test]
    fn the_composition_is_the_shipped_answer_where_nothing_rounds() -> Result<(), String> {
        let rec = record(
            &[
                ("q", &[1.0, 0.0, 0.0, 0.0]),
                ("x", &[1.0, -2.0, 3.0]),
                ("p", &[1.5, -2.0, 0.25]),
            ],
            &[],
        )?;
        for precision in [Precision::F64, Precision::F32] {
            for (id, _) in InvActId::ALL {
                let (got, want) = (
                    ComposedInvAct.eval(id, &rec, precision),
                    Helicoid.eval(id, &rec, precision),
                );
                assert_eq!(
                    got.keys().collect::<Vec<_>>(),
                    want.keys().collect::<Vec<_>>()
                );
                for (k, v) in &want {
                    assert_eq!(&got[k], v, "{id} {k}");
                }
            }
        }
        assert!(!ComposedInvAct.supports("so3_act"));
        Ok(())
    }
}
